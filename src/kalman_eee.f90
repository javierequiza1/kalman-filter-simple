program kalman_complete_system
    implicit none
    
    ! --- PARÁMETROS GLOBALES ---
    integer, parameter :: dp = kind(0.d0)
    integer, parameter :: n_steps = 2000
    integer, parameter :: MAX_LAG = 500
    integer, parameter :: BUF_SIZE = MAX_LAG + 10
    
    ! --- VARIABLES DEL SISTEMA ---
    real(dp) :: dt, t_now, pi_val
    real(dp) :: F(2,2), Q(2,2), H(2), R_noise
    real(dp) :: m(2), P(2,2)           ! Estado Actual
    real(dp) :: m_pred(2), P_pred(2,2) ! Predicciones (Necesarias para Smoother)
    real(dp) :: z_meas
    
    ! --- BUFFERS CIRCULARES (HISTORIAL) ---
    real(dp) :: hist_m(2, BUF_SIZE), hist_P(2,2, BUF_SIZE)
    real(dp) :: hist_mp(2, BUF_SIZE), hist_Pp(2,2, BUF_SIZE)
    real(dp) :: hist_z(BUF_SIZE), hist_t(BUF_SIZE)
    
    ! --- VARIABLES DE CONTROL ---
    ! IMPORTANTE: He añadido 'k' aquí para el bucle de vaciado
    integer :: user_lag, i, k, p_head, p_target, steps_filled
    real(dp) :: output_state(2)
    logical :: use_smoother

    ! ==========================================================================
    ! 1. CONFIGURACIÓN
    ! ==========================================================================
    pi_val = 4.0_dp * atan(1.0_dp)
    dt = 2*pi_val / real(n_steps, dp)
    
    ! -- PREGUNTAR AL USUARIO --
    print *, "--- SISTEMA KALMAN COMPLETO ---"
    print *, "Introduce LAG (0 = Online Puro, >0 = Smoother/Retrasado):"
    read *, user_lag
    
    if (user_lag > MAX_LAG) user_lag = MAX_LAG
    if (user_lag < 0) user_lag = 0
    
    use_smoother = (user_lag > 0)
    print *, "Modo configurado: Lag =", user_lag

    ! -- INICIALIZACIÓN DE MATRICES --
    ! Modelo: Velocidad Constante
    F(1,1) = 1.0_dp; F(1,2) = dt
    F(2,1) = 0.0_dp; F(2,2) = 1.0_dp
    H = [1.0_dp, 0.0_dp]
    
    ! -- SINTONIZACIÓN (TUNING) AGRESIVA PARA CURVAS --
    Q = 0.0_dp
    Q(1,1) = 0.0001_dp    ! Incertidumbre en posición
    Q(2,2) = 0.01d0       ! ALTA Incertidumbre en velocidad
    R_noise = 0.05_dp     ! Confianza media en el sensor

    ! -- ESTADO INICIAL --
    m = 0.0_dp
    P = 0.0_dp; P(1,1)=1.0_dp; P(2,2)=1.0_dp
    t_now = 0.0_dp
    
    ! Punteros del buffer
    p_head = 0
    steps_filled = 0
    
    open(unit=10, file='kalman_output.txt', status='replace')

    ! ==========================================================================
    ! 2. BUCLE PRINCIPAL
    ! ==========================================================================
    do i = 1, n_steps
        
        ! A. SIMULACIÓN Y MEDICIÓN
        t_now = t_now + dt
        ! Señal real (Coseno) + Ruido
        z_meas = log(t_now) + sqrt(R_noise) * gaussian_noise()
        
        ! B. LLAMADA AL FILTRO ONLINE (Subrutina aislada)
        call run_kalman_step(m, P, m_pred, P_pred, z_meas, F, Q, H, R_noise)
        
        ! C. GUARDAR EN HISTORIAL (BUFFER CIRCULAR)
        p_head = mod(p_head, BUF_SIZE) + 1
        steps_filled = min(steps_filled + 1, BUF_SIZE)
        
        hist_m(:, p_head)   = m
        hist_P(:,:, p_head) = P
        hist_mp(:, p_head)  = m_pred
        hist_Pp(:,:, p_head)= P_pred
        hist_z(p_head)      = z_meas
        hist_t(p_head)      = t_now
        
        ! D. PROCESAR SALIDA (SMOOTHER O ONLINE)
        
        if (use_smoother) then
            ! Solo podemos suavizar si tenemos suficientes datos almacenados
            if (steps_filled > user_lag) then
                ! Calculamos el índice del pasado que queremos corregir
                p_target = mod(p_head - user_lag - 1 + BUF_SIZE, BUF_SIZE) + 1
                
                ! Llamada al Smoother (Calcula el valor mejorado para p_target)
                output_state = calculate_rts_smooth(p_head, user_lag, BUF_SIZE, &
                                                    hist_m, hist_P, hist_mp, hist_Pp, F)
                                                    
                ! Escribimos: Tiempo Antiguo, Medición Antigua, Estimación Suavizada
                write(10, *) hist_t(p_target), hist_z(p_target), output_state(1)
            end if
        else
            ! MODO ONLINE PURO: Escribimos lo que acaba de salir del filtro
            write(10, *) t_now, z_meas, m(1)
        end if
        
    end do
    
    ! ==========================================================================
    ! 3. VACIADO DEL BUFFER FINAL (FLUSH) -- [NUEVO BLOQUE]
    ! ==========================================================================
    ! Al terminar el bucle, quedan 'user_lag' datos en el buffer sin escribir.
    ! Los procesamos reduciendo el lag progresivamente hasta 0.
    
    if (use_smoother .and. steps_filled >= user_lag) then
        print *, "Simulacion terminada. Vaciando cola de suavizado..."
        
        ! Recorremos desde (Lag - 1) hasta 0 (el ultimo dato)
        do k = user_lag - 1, 0, -1
            
            ! Calculamos el índice objetivo (avanzando hacia el final)
            p_target = mod(p_head - k - 1 + BUF_SIZE, BUF_SIZE) + 1
            
            ! Llamamos al Smoother con un lag variable 'k'
            ! El punto de referencia 'future' sigue siendo p_head (el final de la simulación)
            output_state = calculate_rts_smooth(p_head, k, BUF_SIZE, &
                                                hist_m, hist_P, hist_mp, hist_Pp, F)
            
            write(10, *) hist_t(p_target), hist_z(p_target), output_state(1)
        end do
    end if

    close(10)
    print *, "Proceso terminado. Datos guardados en 'kalman_output.txt'"

contains

    ! ----------------------------------------------------------------------
    ! SUBRUTINA: UN PASO DEL FILTRO DE KALMAN (ONLINE)
    ! ----------------------------------------------------------------------
    subroutine run_kalman_step(x, P_cov, x_pred, P_prior, z, F_mat, Q_mat, H_vec, R_scalar)
        real(dp), intent(inout) :: x(2), P_cov(2,2)
        real(dp), intent(out)   :: x_pred(2), P_prior(2,2)
        real(dp), intent(in)    :: z, F_mat(2,2), Q_mat(2,2), H_vec(2), R_scalar
        
        real(dp) :: y_innov, S_cov, K_gain(2)
        real(dp) :: PHt(2)
        
        ! 1. Predicción (Time Update)
        x_pred = matmul(F_mat, x)
        P_prior = matmul(matmul(F_mat, P_cov), transpose(F_mat)) + Q_mat
        
        ! 2. Actualización (Measurement Update)
        y_innov = z - (x_pred(1)*H_vec(1) + x_pred(2)*H_vec(2))
        
        ! Auxiliar P * H^T
        PHt(1) = P_prior(1,1)*H_vec(1) + P_prior(1,2)*H_vec(2)
        PHt(2) = P_prior(2,1)*H_vec(1) + P_prior(2,2)*H_vec(2)
        
        S_cov = (H_vec(1)*PHt(1) + H_vec(2)*PHt(2)) + R_scalar
        
        K_gain(1) = PHt(1) / S_cov
        K_gain(2) = PHt(2) / S_cov
        
        ! Estado actualizado
        x = x_pred + (K_gain * y_innov)
        
        ! Covarianza actualizada: P = (I - KH)P_prior
        P_cov(1,1) = P_prior(1,1) - K_gain(1)*S_cov*K_gain(1)
        P_cov(1,2) = P_prior(1,2) - K_gain(1)*S_cov*K_gain(2)
        P_cov(2,1) = P_prior(2,1) - K_gain(2)*S_cov*K_gain(1)
        P_cov(2,2) = P_prior(2,2) - K_gain(2)*S_cov*K_gain(2)
        
    end subroutine run_kalman_step

    ! ----------------------------------------------------------------------
    ! FUNCIÓN: RTS SMOOTHER RETROSPECTIVO
    ! Devuelve el estado suavizado en t-lag dado el buffer actual
    ! ----------------------------------------------------------------------
    function calculate_rts_smooth(curr_idx, lag, dim, hm, hP, hmp, hPp, F_sys) result(x_smooth)
        integer, intent(in) :: curr_idx, lag, dim
        real(dp), intent(in) :: hm(2,dim), hP(2,2,dim), hmp(2,dim), hPp(2,2,dim), F_sys(2,2)
        real(dp) :: x_smooth(2)
        
        integer :: j, idx_now, idx_future
        real(dp) :: x_curr(2), P_curr(2,2)
        real(dp) :: P_fut_pred(2,2), x_fut_smooth(2)
        real(dp) :: C_mat(2,2), det_P, P_inv(2,2)
        real(dp) :: diff_x(2)
        
        ! Empezamos desde el "futuro" (la cabeza del buffer)
        x_fut_smooth = hm(:, curr_idx)
        
        ! Retrocedemos paso a paso hasta llegar al punto de lag
        do j = 0, lag-1
            ! Índices circulares
            idx_future = mod(curr_idx - j - 1 + dim, dim) + 1
            idx_now    = mod(curr_idx - j - 2 + dim, dim) + 1
            
            ! Recuperamos estado filtrado en t y predicción en t+1
            x_curr = hm(:, idx_now)
            P_curr = hP(:,:, idx_now)
            
            ! CORRECCIÓN IMPORTANTE: Necesitamos la P_pred del FUTURO (idx_future)
            ! En tu código original ponía 'idx_now', lo cual es incorrecto para RTS
            P_fut_pred = hPp(:,:, idx_future) 
            
            ! Inversa de P_pred (2x2 simplificada)
            det_P = P_fut_pred(1,1)*P_fut_pred(2,2) - P_fut_pred(1,2)*P_fut_pred(2,1)
            if (abs(det_P) < 1.0d-20) det_P = 1.0d-20
            
            P_inv(1,1) =  P_fut_pred(2,2)/det_P
            P_inv(1,2) = -P_fut_pred(1,2)/det_P
            P_inv(2,1) = -P_fut_pred(2,1)/det_P
            P_inv(2,2) =  P_fut_pred(1,1)/det_P
            
            ! Matriz de Suavizado C = P_filt * F^T * P_pred^-1
            C_mat = matmul(matmul(P_curr, transpose(F_sys)), P_inv)
            
            ! Estado suavizado: x_s = x_filt + C * (x_future_s - x_future_pred)
            diff_x = x_fut_smooth - hmp(:, idx_future)
            x_curr = x_curr + matmul(C_mat, diff_x)
            
            ! Actualizamos para la siguiente iteración
            x_fut_smooth = x_curr
        end do
        
        x_smooth = x_fut_smooth
    end function calculate_rts_smooth

    ! ----------------------------------------------------------------------
    ! GENERADOR DE RUIDO
    ! ----------------------------------------------------------------------
    function gaussian_noise() result(v)
        real(dp) :: u1, u2, v
        call random_number(u1)
        call random_number(u2)
        if (u1 < 1.0d-30) u1 = 1.0d-30
        v = sqrt(-2.0_dp*log(u1)) * cos(6.283185307_dp*u2)
    end function

end program