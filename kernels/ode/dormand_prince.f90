! fortran/kernels/ode/dormand_prince.f90
!
! Dormand-Prince method (DOPRI5) - adaptive step-size Runge-Kutta solver.
! 5th-order method with 4th-order embedded error estimation.
!
! This is the default ODE solver in MATLAB (ode45) and SciPy (solve_ivp).
! Excellent for non-stiff problems with moderate accuracy requirements.
!
! For the ODE system: dy/dt = f(t, y)
!
! Usage:
!   use dormand_prince
!   type(dopri_config_t) :: config
!   real(wp) :: y(n), t_span(2)
!   real(wp), allocatable :: t_out(:), y_out(:,:)
!   
!   config = dopri_config_t(rtol=1e-6_wp, atol=1e-8_wp)
!   t_span = [0.0_wp, 10.0_wp]
!   call dopri_solve(ode_func, t_span, y, config, t_out, y_out, status)
!
module dormand_prince
    use kinds, only: wp, i64
    use constants, only: TOL_DEFAULT
    implicit none
    private
    
    ! Public interface
    public :: dopri_config_t, dopri_status_t
    public :: dopri_solve, dopri_step
    public :: DOPRI_SUCCESS, DOPRI_ERR_MAX_STEPS, DOPRI_ERR_INVALID_ARG
    public :: DOPRI_ERR_STEP_TOO_SMALL, DOPRI_ERR_FUNC_EVAL
    
    ! Status codes
    integer, parameter :: DOPRI_SUCCESS = 0
    integer, parameter :: DOPRI_ERR_MAX_STEPS = 1
    integer, parameter :: DOPRI_ERR_INVALID_ARG = 2
    integer, parameter :: DOPRI_ERR_STEP_TOO_SMALL = 3
    integer, parameter :: DOPRI_ERR_FUNC_EVAL = 4
    
    ! Dormand-Prince coefficients (Butcher tableau)
    real(wp), parameter :: A2 = 1.0_wp / 5.0_wp
    real(wp), parameter :: A3 = 3.0_wp / 10.0_wp
    real(wp), parameter :: A4 = 4.0_wp / 5.0_wp
    real(wp), parameter :: A5 = 8.0_wp / 9.0_wp
    
    real(wp), parameter :: B21 = 1.0_wp / 5.0_wp
    real(wp), parameter :: B31 = 3.0_wp / 40.0_wp
    real(wp), parameter :: B32 = 9.0_wp / 40.0_wp
    real(wp), parameter :: B41 = 44.0_wp / 45.0_wp
    real(wp), parameter :: B42 = -56.0_wp / 15.0_wp
    real(wp), parameter :: B43 = 32.0_wp / 9.0_wp
    real(wp), parameter :: B51 = 19372.0_wp / 6561.0_wp
    real(wp), parameter :: B52 = -25360.0_wp / 2187.0_wp
    real(wp), parameter :: B53 = 64448.0_wp / 6561.0_wp
    real(wp), parameter :: B54 = -212.0_wp / 729.0_wp
    real(wp), parameter :: B61 = 9017.0_wp / 3168.0_wp
    real(wp), parameter :: B62 = -355.0_wp / 33.0_wp
    real(wp), parameter :: B63 = 46732.0_wp / 5247.0_wp
    real(wp), parameter :: B64 = 49.0_wp / 176.0_wp
    real(wp), parameter :: B65 = -5103.0_wp / 18656.0_wp
    real(wp), parameter :: B71 = 35.0_wp / 384.0_wp
    real(wp), parameter :: B72 = 0.0_wp
    real(wp), parameter :: B73 = 500.0_wp / 1113.0_wp
    real(wp), parameter :: B74 = 125.0_wp / 192.0_wp
    real(wp), parameter :: B75 = -2187.0_wp / 6784.0_wp
    real(wp), parameter :: B76 = 11.0_wp / 84.0_wp
    
    ! 5th-order solution coefficients
    real(wp), parameter :: C1 = 35.0_wp / 384.0_wp
    real(wp), parameter :: C2 = 0.0_wp
    real(wp), parameter :: C3 = 500.0_wp / 1113.0_wp
    real(wp), parameter :: C4 = 125.0_wp / 192.0_wp
    real(wp), parameter :: C5 = -2187.0_wp / 6784.0_wp
    real(wp), parameter :: C6 = 11.0_wp / 84.0_wp
    real(wp), parameter :: C7 = 0.0_wp
    
    ! 4th-order solution coefficients (for error estimation)
    real(wp), parameter :: D1 = 5179.0_wp / 57600.0_wp
    real(wp), parameter :: D2 = 0.0_wp
    real(wp), parameter :: D3 = 7571.0_wp / 16695.0_wp
    real(wp), parameter :: D4 = 393.0_wp / 640.0_wp
    real(wp), parameter :: D5 = -92097.0_wp / 339200.0_wp
    real(wp), parameter :: D6 = 187.0_wp / 2100.0_wp
    real(wp), parameter :: D7 = 1.0_wp / 40.0_wp
    
    ! Error coefficients (difference between 5th and 4th order)
    real(wp), parameter :: E1 = C1 - D1
    real(wp), parameter :: E2 = C2 - D2
    real(wp), parameter :: E3 = C3 - D3
    real(wp), parameter :: E4 = C4 - D4
    real(wp), parameter :: E5 = C5 - D5
    real(wp), parameter :: E6 = C6 - D6
    real(wp), parameter :: E7 = C7 - D7
    
    ! Configuration type
    type :: dopri_config_t
        real(wp) :: rtol = 1.0e-6_wp       ! Relative tolerance
        real(wp) :: atol = 1.0e-8_wp       ! Absolute tolerance
        real(wp) :: dt_init = 0.0_wp       ! Initial step (0 = auto)
        real(wp) :: dt_max = huge(1.0_wp)  ! Maximum step size
        real(wp) :: dt_min = 0.0_wp        ! Minimum step size (0 = auto)
        real(wp) :: safety = 0.9_wp        ! Safety factor for step control
        integer :: max_steps = 100000      ! Maximum integration steps
        logical :: dense_output = .false.  ! Store all accepted steps
    end type dopri_config_t
    
    ! Status information
    type :: dopri_status_t
        integer :: code = DOPRI_SUCCESS
        integer :: steps_taken = 0
        integer :: steps_accepted = 0
        integer :: steps_rejected = 0
        integer :: func_evals = 0
        real(wp) :: final_time = 0.0_wp
        real(wp) :: final_step_size = 0.0_wp
    end type dopri_status_t
    
    ! Abstract interface for ODE function
    abstract interface
        subroutine ode_func_interface(t, y, dydt)
            import :: wp
            real(wp), intent(in) :: t
            real(wp), intent(in) :: y(:)
            real(wp), intent(out) :: dydt(:)
        end subroutine ode_func_interface
    end interface
    
contains

    !> Solve ODE system using adaptive Dormand-Prince method
    subroutine dopri_solve(func, t_span, y0, config, t_out, y_out, status)
        procedure(ode_func_interface) :: func
        real(wp), intent(in) :: t_span(2)
        real(wp), intent(in) :: y0(:)
        type(dopri_config_t), intent(in) :: config
        real(wp), allocatable, intent(out) :: t_out(:)
        real(wp), allocatable, intent(out) :: y_out(:, :)
        type(dopri_status_t), intent(out) :: status
        
        real(wp) :: t, tf, dt, dt_min
        real(wp) :: y(size(y0)), y_new(size(y0))
        real(wp) :: error_norm, scale_factor
        real(wp), allocatable :: t_temp(:), y_temp(:, :)
        integer :: n, step, output_idx, max_outputs
        logical :: step_accepted
        
        ! Initialize
        n = size(y0)
        t = t_span(1)
        tf = t_span(2)
        y = y0
        status = dopri_status_t()
        
        ! Validate inputs
        if (tf <= t .or. n <= 0 .or. config%rtol <= 0.0_wp .or. config%atol < 0.0_wp) then
            status%code = DOPRI_ERR_INVALID_ARG
            return
        end if
        
        ! Set minimum step size
        dt_min = config%dt_min
        if (dt_min <= 0.0_wp) then
            dt_min = 1.0e-12_wp * abs(tf - t)
        end if
        
        ! Estimate initial step size if not provided
        if (config%dt_init > 0.0_wp) then
            dt = config%dt_init
        else
            dt = estimate_initial_step(func, t, y, tf, config%rtol, config%atol)
        end if
        dt = min(dt, config%dt_max)
        
        ! Allocate temporary output storage (will resize later)
        max_outputs = min(10000, config%max_steps)
        allocate(t_temp(max_outputs))
        allocate(y_temp(n, max_outputs))
        
        ! Store initial condition
        output_idx = 1
        t_temp(output_idx) = t
        y_temp(:, output_idx) = y
        
        ! Integration loop
        do step = 1, config%max_steps
            ! Adjust step if we'd overshoot
            if (t + dt > tf) then
                dt = tf - t
            end if
            
            ! Take one DOPRI step and compute error
            call dopri_step(func, t, y, dt, y_new, error_norm)
            
            status%steps_taken = step
            status%func_evals = status%func_evals + 7
            
            ! Compute scale factor for step adjustment
            scale_factor = config%safety * (1.0_wp / max(error_norm, 1.0e-10_wp))**0.2_wp
            scale_factor = min(5.0_wp, max(0.2_wp, scale_factor))
            
            ! Check if step is accepted
            step_accepted = (error_norm <= 1.0_wp)
            
            if (step_accepted) then
                ! Accept step
                t = t + dt
                y = y_new
                status%steps_accepted = status%steps_accepted + 1
                
                ! Store output if dense output is requested
                if (config%dense_output) then
                    output_idx = output_idx + 1
                    if (output_idx > max_outputs) then
                        ! Resize arrays
                        call resize_output_arrays(t_temp, y_temp, max_outputs * 2)
                        max_outputs = max_outputs * 2
                    end if
                    t_temp(output_idx) = t
                    y_temp(:, output_idx) = y
                end if
                
                ! Check if we've reached final time
                if (abs(t - tf) < epsilon(1.0_wp)) then
                    status%code = DOPRI_SUCCESS
                    status%final_time = t
                    status%final_step_size = dt
                    exit
                end if
            else
                ! Reject step
                status%steps_rejected = status%steps_rejected + 1
            end if
            
            ! Adjust step size for next iteration
            dt = dt * scale_factor
            dt = min(dt, config%dt_max)
            
            ! Check if step size became too small
            if (dt < dt_min) then
                status%code = DOPRI_ERR_STEP_TOO_SMALL
                status%final_time = t
                status%final_step_size = dt
                exit
            end if
        end do
        
        ! Check if max steps exceeded
        if (status%code == DOPRI_SUCCESS .and. abs(t - tf) >= epsilon(1.0_wp)) then
            status%code = DOPRI_ERR_MAX_STEPS
            status%final_time = t
        end if
        
        ! Store final point if not already stored
        if (.not. config%dense_output .or. abs(t_temp(output_idx) - t) > epsilon(1.0_wp)) then
            output_idx = output_idx + 1
            if (output_idx > max_outputs) then
                call resize_output_arrays(t_temp, y_temp, output_idx)
            end if
            t_temp(output_idx) = t
            y_temp(:, output_idx) = y
        end if
        
        ! Copy to final output arrays
        allocate(t_out(output_idx))
        allocate(y_out(n, output_idx))
        t_out = t_temp(1:output_idx)
        y_out = y_temp(:, 1:output_idx)
        
        deallocate(t_temp, y_temp)
    end subroutine dopri_solve
    
    !> Perform single DOPRI5 step with error estimation
    pure subroutine dopri_step(func, t, y, dt, y_new, error_norm)
        procedure(ode_func_interface) :: func
        real(wp), intent(in) :: t
        real(wp), intent(in) :: y(:)
        real(wp), intent(in) :: dt
        real(wp), intent(out) :: y_new(:)
        real(wp), intent(out) :: error_norm
        
        real(wp) :: k1(size(y)), k2(size(y)), k3(size(y)), k4(size(y))
        real(wp) :: k5(size(y)), k6(size(y)), k7(size(y))
        real(wp) :: y_temp(size(y)), error(size(y)), scale(size(y))
        
        ! k1 = f(t, y)
        call func(t, y, k1)
        
        ! k2 = f(t + A2*dt, y + dt*(B21*k1))
        y_temp = y + dt * B21 * k1
        call func(t + A2 * dt, y_temp, k2)
        
        ! k3 = f(t + A3*dt, y + dt*(B31*k1 + B32*k2))
        y_temp = y + dt * (B31 * k1 + B32 * k2)
        call func(t + A3 * dt, y_temp, k3)
        
        ! k4 = f(t + A4*dt, y + dt*(B41*k1 + B42*k2 + B43*k3))
        y_temp = y + dt * (B41 * k1 + B42 * k2 + B43 * k3)
        call func(t + A4 * dt, y_temp, k4)
        
        ! k5 = f(t + A5*dt, y + dt*(B51*k1 + B52*k2 + B53*k3 + B54*k4))
        y_temp = y + dt * (B51 * k1 + B52 * k2 + B53 * k3 + B54 * k4)
        call func(t + A5 * dt, y_temp, k5)
        
        ! k6 = f(t + dt, y + dt*(B61*k1 + B62*k2 + B63*k3 + B64*k4 + B65*k5))
        y_temp = y + dt * (B61 * k1 + B62 * k2 + B63 * k3 + B64 * k4 + B65 * k5)
        call func(t + dt, y_temp, k6)
        
        ! 5th-order solution
        y_new = y + dt * (C1 * k1 + C3 * k3 + C4 * k4 + C5 * k5 + C6 * k6)
        
        ! k7 = f(t + dt, y_new) for FSAL property
        call func(t + dt, y_new, k7)
        
        ! Error estimate (difference between 5th and 4th order)
        error = dt * (E1 * k1 + E3 * k3 + E4 * k4 + E5 * k5 + E6 * k6 + E7 * k7)
        
        ! Compute error norm with mixed relative/absolute tolerance
        scale = max(abs(y), abs(y_new))
        error_norm = sqrt(sum((error / (1.0e-8_wp + 1.0e-6_wp * scale))**2) / size(y))
    end subroutine dopri_step
    
    !> Estimate good initial step size
    function estimate_initial_step(func, t0, y0, tf, rtol, atol) result(dt)
        procedure(ode_func_interface) :: func
        real(wp), intent(in) :: t0, y0(:), tf, rtol, atol
        real(wp) :: dt
        
        real(wp) :: f0(size(y0)), scale, d0, d1, d2
        real(wp) :: y1(size(y0)), f1(size(y0))
        
        ! Compute f0 = f(t0, y0)
        call func(t0, y0, f0)
        
        ! Compute scale
        scale = atol + rtol * maxval(abs(y0))
        
        ! d0 = ||y0||
        d0 = sqrt(sum((y0 / scale)**2) / size(y0))
        
        ! d1 = ||f0||
        d1 = sqrt(sum((f0 / scale)**2) / size(y0))
        
        ! Initial estimate
        if (d0 < 1.0e-5_wp .or. d1 < 1.0e-5_wp) then
            dt = 1.0e-6_wp * abs(tf - t0)
        else
            dt = 0.01_wp * d0 / d1
        end if
        
        ! Refine estimate with Euler step
        y1 = y0 + dt * f0
        call func(t0 + dt, y1, f1)
        
        d2 = sqrt(sum(((f1 - f0) / scale)**2) / size(y0)) / dt
        
        if (max(d1, d2) <= 1.0e-15_wp) then
            dt = max(1.0e-6_wp * abs(tf - t0), dt * 1.0e-3_wp)
        else
            dt = (0.01_wp / max(d1, d2))**(1.0_wp / 5.0_wp)
        end if
        
        ! Ensure reasonable bounds
        dt = min(100.0_wp * dt, abs(tf - t0))
    end function estimate_initial_step
    
    !> Resize output arrays (helper for dynamic allocation)
    subroutine resize_output_arrays(t_arr, y_arr, new_size)
        real(wp), allocatable, intent(inout) :: t_arr(:)
        real(wp), allocatable, intent(inout) :: y_arr(:, :)
        integer, intent(in) :: new_size
        
        real(wp), allocatable :: t_new(:), y_new(:, :)
        integer :: n, old_size
        
        old_size = size(t_arr)
        n = size(y_arr, 1)
        
        allocate(t_new(new_size))
        allocate(y_new(n, new_size))
        
        t_new(1:old_size) = t_arr
        y_new(:, 1:old_size) = y_arr
        
        call move_alloc(t_new, t_arr)
        call move_alloc(y_new, y_arr)
    end subroutine resize_output_arrays

end module dormand_prince