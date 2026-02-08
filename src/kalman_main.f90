program kalman_continuous
    implicit none 

    integer, parameter :: n = 2000
    real(8), parameter :: pi = 3.141592653589793d0
    
    real(8) :: t(n), z(n), true_signal(n), error_val
    real(8) :: m(2), P(2,2), Q(2,2), F(2,2)
    
    ! Mantenemos R como vector por si quieres experimentar luego,
    ! pero ahora tendrá el mismo valor siempre.
    real(8) :: R_vec(n) 
    
    real(8) :: dt, var_val
    real(8) :: m_pred(2), P_pred(2,2)
    real(8) :: H(2), S, y_res, K(2), PHt(2)
    integer :: i

    ! 1. CONFIGURACIÓN (Tus valores actuales)
    Q(1,1) = 1.0d-8
    Q(1,2) = 0.0d0; Q(2,1) = 0.0d0
    Q(2,2) = 1.0d-1
    
    ! 2. GENERACIÓN DE DATOS
    var_val = 0.3d0
    t(1) = 0.0d0
    
    do i = 2, n
        t(i) = t(i-1) + (4.0d0 * pi / real(n, 8))
    end do
    
    do i = 1, n
        true_signal(i) = sin(t(i))
        z(i) = signal_func(t(i), var_val)
    end do
    dt = t(2) - t(1)

    ! 3. CONFIGURACIÓN DEL SENSOR (MODIFICADO)
    ! Asignamos el ruido estándar a TODO el vector.
    ! Hemos borrado el bloque 'WHERE' que creaba el apagón.
    R_vec = 0.09d0 

    ! Inicialización del filtro
    F(1,1) = 1.0d0; F(1,2) = dt
    F(2,1) = 0.0d0; F(2,2) = 1.0d0
    m = 0.0d0
    P = 0.0d0; P(1,1)=1.0d0; P(2,2)=1.0d0
    H(1) = 1.0d0; H(2) = 0.0d0

    open(40, file='kalman_continuous.txt', status='replace')

    ! 4. BUCLE PRINCIPAL
    do i = 1, n
        
        ! --- PREDICCIÓN ---
        m_pred = matmul(F, m)
        P_pred = matmul(matmul(F, P), transpose(F)) + Q

        ! --- ACTUALIZACIÓN ---
        y_res = z(i) - m_pred(1)
        
        PHt(1) = P_pred(1,1)*H(1) + P_pred(1,2)*H(2)
        PHt(2) = P_pred(2,1)*H(1) + P_pred(2,2)*H(2)
        
        ! R_vec(i) siempre es 0.09, así que siempre corrige
        S = (H(1)*PHt(1) + H(2)*PHt(2)) + R_vec(i)
        
        K(1) = PHt(1) / S
        K(2) = PHt(2) / S
        
        m = m_pred + (K * y_res)
        
        P(1,1) = P_pred(1,1) - K(1)*PHt(1)
        P(1,2) = P_pred(1,2) - K(1)*PHt(2)
        P(2,1) = P_pred(2,1) - K(2)*PHt(1)
        P(2,2) = P_pred(2,2) - K(2)*PHt(2)

        ! Calculamos y guardamos el error
        error_val = m(1) - true_signal(i)
        
        write(40, *) t(i), z(i), m(1), true_signal(i), error_val
    end do

    close(40)
    print *, "Calculo continuo completado. Archivo: kalman_continuous.txt"

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

end program kalman_continuous