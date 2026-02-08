program kalman_optimizer
    implicit none 

    ! Parámetros de la simulación
    integer, parameter :: n = 2000
    real(8), parameter :: pi = 3.141592653589793d0
    
    ! Arrays de datos
    real(8) :: t(n), z(n), true_signal(n)
    
    ! Variables del Filtro de Kalman
    real(8) :: m(2), P(2,2), Q(2,2), F(2,2), R_noise
    real(8) :: dt
    
    ! Variables para el Grid Search
    real(8) :: q1_exp, q2_exp       ! El valor real del exponente
    real(8) :: q1_val, q2_val       ! El valor real de Q (10^exp)
    real(8) :: mse, rmse            
    real(8) :: error_inst           
    
    ! --- NUEVAS VARIABLES PARA LOS BUCLES ---
    integer :: i1, i2              ! Contadores enteros para los bucles externos
    ! Rango: de -9.0 a 1.0 son 10 unidades. Paso 0.2. Total 50 pasos.
    integer, parameter :: steps = 50 
    
    ! Configuración del ruido
    real(8) :: var_val 
    integer :: i, k
    
    ! --- 1. CONFIGURACIÓN INICIAL ---
    var_val = 0.3d0          
    R_noise = 0.09d0         
    
    ! Generación de datos
    t(1) = 0.0d0
    true_signal(1) = sin(t(1))
    z(1) = signal_func(t(1), var_val)
    
    do i = 2, n
        t(i) = t(i-1) + (4.0d0 * pi / real(n, 8))
        true_signal(i) = sin(t(i))
        z(i) = signal_func(t(i), var_val)
    end do
    
    dt = t(2) - t(1)
    
    ! Matriz F Cinemática
    F(1,1) = 1.0d0; F(1,2) = dt
    F(2,1) = 0.0d0; F(2,2) = 1.0d0

    print *, "Iniciando Grid Search masivo (Corregido)..."
    open(20, file='error_surface.txt', status='replace')
    write(20, *) "Log_Q_Pos Log_Q_Vel RMSE" 

    ! --- 2. BUCLE DE OPTIMIZACIÓN CON ENTEROS ---
    ! Vamos de 0 a 50 pasos.
    ! Fórmula: valor = inicio + (paso * contador)
    ! Inicio = -9.0, Paso = 0.2
    
    do i1 = 0, steps
        ! Calculamos el exponente real aquí dentro
        q1_exp = -9.0d0 + (real(i1, 8) * 0.2d0)
        
        do i2 = 0, steps
            q2_exp = -9.0d0 + (real(i2, 8) * 0.2d0)
            
            ! A. Configurar Q
            Q = 0.0d0
            q1_val = 10.0d0**q1_exp
            q2_val = 10.0d0**q2_exp
            
            Q(1,1) = q1_val
            Q(2,2) = q2_val
            
            ! B. Reiniciar Filtro
            m = 0.0d0      
            P = 0.0d0
            P(1,1) = 1.0d0 
            P(2,2) = 1.0d0
            mse = 0.0d0    
            
            ! C. Simulación
            do k = 1, n
                call kalman_step_cinematico(z(k), F, m, P, Q, R_noise)
                
                if (k > 50) then
                    error_inst = m(1) - true_signal(k)
                    mse = mse + (error_inst**2)
                end if
            end do
            
            ! D. Calcular RMSE
            rmse = sqrt(mse / real(n-50, 8))
            
            ! E. Guardar
            write(20, *) q1_exp, q2_exp, rmse
            
        end do
        
        print *, "Procesado Q_Pos 10^", q1_exp
        write(20, *) "" 
    end do

    close(20)
    print *, "Calculo finalizado. Datos en error_surface.txt"

contains

    function signal_func(val, v) result(f)
        real(8), intent(in) :: val, v
        real(8) :: f
        f = sin(val) + v * gaussian_noise()
    end function signal_func

    function gaussian_noise() result(z_noise)
        real(8) :: u1, u2, z_noise
        call random_number(u1)
        call random_number(u2)
        z_noise = sqrt(-2.0d0 * log(u1)) * cos(2.0d0 * pi * u2)
    end function gaussian_noise

    subroutine kalman_step_cinematico(y, F, m, P, Q, R)
        real(8), intent(in)    :: y, R
        real(8), intent(in)    :: F(2,2), Q(2,2)
        real(8), intent(inout) :: m(2), P(2,2)
        
        real(8) :: m_pred(2), P_pred(2,2)
        real(8) :: H(2), S, y_res
        real(8) :: K(2), PHt(2)
        
        m_pred = matmul(F, m)
        P_pred = matmul(matmul(F, P), transpose(F)) + Q
        
        H(1) = 1.0d0
        H(2) = 0.0d0
        
        y_res = y - m_pred(1)
        
        PHt(1) = P_pred(1,1)*H(1) + P_pred(1,2)*H(2)
        PHt(2) = P_pred(2,1)*H(1) + P_pred(2,2)*H(2)
        
        S = (H(1)*PHt(1) + H(2)*PHt(2)) + R
        
        K(1) = PHt(1) / S
        K(2) = PHt(2) / S
        
        m = m_pred + (K * y_res)
        
        P(1,1) = P_pred(1,1) - K(1)*PHt(1)
        P(1,2) = P_pred(1,2) - K(1)*PHt(2)
        P(2,1) = P_pred(2,1) - K(2)*PHt(1)
        P(2,2) = P_pred(2,2) - K(2)*PHt(2)
        
    end subroutine kalman_step_cinematico

end program kalman_optimizer