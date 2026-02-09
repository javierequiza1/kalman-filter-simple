program kalman_optimizer_phase_space
    implicit none 

    ! --- PARÁMETROS ---
    integer, parameter :: dp = kind(0.d0) ! Doble precisión (Real*8)
    integer, parameter :: n = 2000
    real(dp), parameter :: pi = 3.141592653589793_dp
    
    ! --- ARRAYS DE DATOS ---
    real(dp) :: t(n), z(n), true_pos(n), true_vel(n)
    
    ! --- KALMAN VARS ---
    real(dp) :: m(2), P(2,2), Q(2,2), F(2,2), R_noise
    real(dp) :: dt
    
    ! --- OPTIMIZACIÓN (GRID SEARCH) ---
    real(dp) :: q1_exp, q2_exp        ! Exponentes (Log10)
    real(dp) :: q1_val, q2_val        ! Valores Reales
    real(dp) :: mse, rmse             
    real(dp) :: err_pos, err_vel      ! Errores individuales
    
    ! --- BUCLES ---
    integer :: i1, i2      ! Contadores optimización
    integer :: steps = 50  ! Resolución del mapa de calor
    integer :: k           ! Contador simulación
    
    ! --- GENERACIÓN DE RUIDO ---
    real(dp) :: var_val 
    
    ! ====================================================================
    ! 1. CONFIGURACIÓN Y GENERACIÓN DE LA "VERDAD"
    ! ====================================================================
    var_val = 0.3_dp      
    R_noise = 0.09_dp     
    
    ! Generamos el eje temporal y la verdad física
    t(1) = 0.0_dp
    true_pos(1) = sin(t(1))
    true_vel(1) = cos(t(1)) ! La derivada exacta
    z(1) = signal_func(t(1), var_val)
    
    do k = 2, n
        t(k) = t(k-1) + (4.0_dp * pi / real(n, dp))
        
        ! VERDAD DE FÍSICA (Estado completo)
        true_pos(k) = sin(t(k))
        true_vel(k) = cos(t(k)) 
        
        ! MEDICIÓN (Solo posición + ruido)
        z(k) = signal_func(t(k), var_val)
    end do
    
    dt = t(2) - t(1)
    
    ! Matriz F (Modelo Velocidad Constante)
    F(1,1) = 1.0_dp; F(1,2) = dt
    F(2,1) = 0.0_dp; F(2,2) = 1.0_dp

    print *, "--- OPTIMIZADOR DE ESPACIO DE FASES (EEE) ---"
    print *, "Minimizando error conjunto: Posicion + Velocidad"
    
    open(20, file='error_surface.txt', status='replace')
    write(20, *) "Log_Q_Pos Log_Q_Vel RMSE_Combined" 

    ! ====================================================================
    ! 2. GRID SEARCH (BÚSQUEDA DE LA Q PERFECTA)
    ! ====================================================================
    ! Exploramos exponentes de -9.0 a 1.0
    
    do i1 = 0, steps
        q1_exp = -9.0_dp + (real(i1, dp) * 0.2_dp)
        
        do i2 = 0, steps
            q2_exp = -9.0_dp + (real(i2, dp) * 0.2_dp)
            
            ! A. Configurar Matriz Q candidata
            Q = 0.0_dp
            Q(1,1) = 10.0_dp**q1_exp
            Q(2,2) = 10.0_dp**q2_exp
            
            ! B. Reiniciar Filtro
            m = 0.0_dp      
            P = 0.0_dp; P(1,1) = 1.0_dp; P(2,2) = 1.0_dp
            mse = 0.0_dp    
            
            ! C. Simulación del Filtro
            do k = 1, n
                call kalman_step_cinematico(z(k), F, m, P, Q, R_noise)
                
                ! D. CÁLCULO DE ERROR (A partir del paso 50 para estabilizar)
                if (k > 50) then
                    err_pos = m(1) - true_pos(k)
                    err_vel = m(2) - true_vel(k) ! <--- ESTO ES LO NUEVO  !!el tema de esto es que conocemos las señales true generalmente esto no es posible habria que ver otros approaches para evaluar el error de velocidad, pero para este ejercicio es válido asumir que la verdad es conocida y así evaluar el error de ambos estados
                    
                    ! Sumamos cuadrados de ambos estados
                    mse = mse + (err_pos**2) + (err_vel**2)
                end if
            end do
            
            ! E. RMSE Combinado
            rmse = sqrt(mse / real(n-50, dp))
            
            ! Guardar resultado
            write(20, *) q1_exp, q2_exp, rmse
            
        end do
        
        ! Feedback visual en terminal para no desesperar
        if (mod(i1, 10) == 0) print *, "Procesando Q_Pos 10^", q1_exp
        write(20, *) "" ! Salto de línea para Gnuplot pm3d
    end do

    close(20)
    print *, "Calculo finalizado. Ejecuta 'heatmap.gp' para ver el optimo."

contains

    ! --- Función auxiliar: Genera señal ruidosa ---
    function signal_func(val, v) result(f)
        real(dp), intent(in) :: val, v
        real(dp) :: f
        f = sin(val) + v * gaussian_noise()
    end function signal_func

    ! --- Función auxiliar: Ruido Gaussiano ---
    function gaussian_noise() result(z_noise)
        real(dp) :: u1, u2, z_noise
        call random_number(u1)
        call random_number(u2)
        if (u1 < 1.0d-30) u1 = 1.0d-30
        z_noise = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * pi * u2)
    end function gaussian_noise

    ! --- SUBRUTINA: UN PASO DE KALMAN ---
    subroutine kalman_step_cinematico(y, F, m, P, Q, R)
        real(dp), intent(in)    :: y, R
        real(dp), intent(in)    :: F(2,2), Q(2,2)
        real(dp), intent(inout) :: m(2), P(2,2)
        
        real(dp) :: m_pred(2), P_pred(2,2)
        real(dp) :: H(2), S, y_res, K(2), PHt(2)
        
        ! 1. Predicción
        m_pred = matmul(F, m)
        P_pred = matmul(matmul(F, P), transpose(F)) + Q
        
        ! 2. Actualización
        H = [1.0_dp, 0.0_dp] ! Solo medimos posición
        
        y_res = y - m_pred(1)
        
        ! K = P*H' / S
        PHt(1) = P_pred(1,1)*H(1) + P_pred(1,2)*H(2)
        PHt(2) = P_pred(2,1)*H(1) + P_pred(2,2)*H(2)
        
        S = (H(1)*PHt(1) + H(2)*PHt(2)) + R
        
        K(1) = PHt(1) / S
        K(2) = PHt(2) / S
        
        m = m_pred + (K * y_res)
        
        ! P = P - K*S*K'
        P(1,1) = P_pred(1,1) - K(1)*S*K(1)
        P(1,2) = P_pred(1,2) - K(1)*S*K(2)
        P(2,1) = P_pred(2,1) - K(2)*S*K(1)
        P(2,2) = P_pred(2,2) - K(2)*S*K(2)
        
    end subroutine kalman_step_cinematico

end program kalman_optimizer_phase_space