! fortran/kernels/ode/backward_euler.f90
!
! Backward (Implicit) Euler method for solving stiff ODEs.
! First-order implicit method with excellent stability properties.
!
! For the ODE system: dy/dt = f(t, y)
! Solves: y_{n+1} = y_n + dt * f(t_{n+1}, y_{n+1})
!
! Uses Newton iteration to solve the nonlinear system at each step.
! Optionally uses provided Jacobian or finite differences.
!
! Usage:
!   use backward_euler
!   type(beuler_config_t) :: config
!   real(wp) :: y(n), t_span(2)
!   real(wp), allocatable :: t_out(:), y_out(:,:)
!   
!   config = beuler_config_t(dt=0.01_wp, newton_tol=1e-6_wp)
!   t_span = [0.0_wp, 10.0_wp]
!   call beuler_solve(ode_func, t_span, y, config, t_out, y_out, status)
!
module backward_euler
    use kinds, only: wp, i64
    use constants, only: TOL_DEFAULT
    implicit none
    private
    
    ! Public interface
    public :: beuler_config_t, beuler_status_t
    public :: beuler_solve, beuler_step
    public :: BEULER_SUCCESS, BEULER_ERR_MAX_STEPS, BEULER_ERR_INVALID_ARG
    public :: BEULER_ERR_NEWTON_FAILED, BEULER_ERR_SINGULAR_JACOBIAN
    
    ! Status codes
    integer, parameter :: BEULER_SUCCESS = 0
    integer, parameter :: BEULER_ERR_MAX_STEPS = 1
    integer, parameter :: BEULER_ERR_INVALID_ARG = 2
    integer, parameter :: BEULER_ERR_NEWTON_FAILED = 3
    integer, parameter :: BEULER_ERR_SINGULAR_JACOBIAN = 4
    
    ! Configuration type
    type :: beuler_config_t
        real(wp) :: dt = 0.01_wp              ! Time step size
        real(wp) :: newton_tol = 1.0e-6_wp    ! Newton iteration tolerance
        real(wp) :: fd_epsilon = 1.0e-8_wp    ! Finite difference epsilon
        integer :: max_steps = 100000         ! Maximum integration steps
        integer :: max_newton_iter = 20       ! Max Newton iterations per step
        integer :: output_every = 1           ! Store output every N steps
        logical :: dense_output = .false.     ! Store all steps
        logical :: use_fd_jacobian = .true.   ! Use finite differences for Jacobian
    end type beuler_config_t
    
    ! Status information
    type :: beuler_status_t
        integer :: code = BEULER_SUCCESS
        integer :: steps_taken = 0
        integer :: newton_iters_total = 0
        integer :: func_evals = 0
        integer :: jacobian_evals = 0
        real(wp) :: final_time = 0.0_wp
    end type beuler_status_t
    
    ! Abstract interfaces
    abstract interface
        subroutine ode_func_interface(t, y, dydt)
            import :: wp
            real(wp), intent(in) :: t
            real(wp), intent(in) :: y(:)
            real(wp), intent(out) :: dydt(:)
        end subroutine ode_func_interface
        
        subroutine jacobian_interface(t, y, jac)
            import :: wp
            real(wp), intent(in) :: t
            real(wp), intent(in) :: y(:)
            real(wp), intent(out) :: jac(:, :)
        end subroutine jacobian_interface
    end interface
    
contains

    !> Solve ODE system using Backward Euler method
    subroutine beuler_solve(func, t_span, y0, config, t_out, y_out, status, jac)
        procedure(ode_func_interface) :: func
        real(wp), intent(in) :: t_span(2)
        real(wp), intent(in) :: y0(:)
        type(beuler_config_t), intent(in) :: config
        real(wp), allocatable, intent(out) :: t_out(:)
        real(wp), allocatable, intent(out) :: y_out(:, :)
        type(beuler_status_t), intent(out) :: status
        procedure(jacobian_interface), optional :: jac
        
        real(wp) :: t, tf, dt
        real(wp) :: y(size(y0))
        integer :: n, step, max_outputs, output_idx, newton_iters
        
        ! Initialize
        n = size(y0)
        t = t_span(1)
        tf = t_span(2)
        dt = config%dt
        y = y0
        status = beuler_status_t()
        
        ! Validate inputs
        if (dt <= 0.0_wp .or. tf <= t .or. n <= 0) then
            status%code = BEULER_ERR_INVALID_ARG
            return
        end if
        
        ! Estimate output size
        if (config%dense_output) then
            max_outputs = ceiling((tf - t) / dt) + 1
        else
            max_outputs = ceiling((tf - t) / (dt * config%output_every)) + 1
        end if
        
        ! Allocate output arrays
        allocate(t_out(max_outputs))
        allocate(y_out(n, max_outputs))
        
        ! Store initial condition
        output_idx = 1
        t_out(output_idx) = t
        y_out(:, output_idx) = y
        
        ! Integration loop
        do step = 1, config%max_steps
            ! Adjust final step if needed
            if (t + dt > tf) then
                dt = tf - t
            end if
            
            ! Take Backward Euler step
            if (present(jac)) then
                call beuler_step_analytical(func, jac, t, y, dt, config, y, &
                                           newton_iters, status%code)
            else
                call beuler_step_fd(func, t, y, dt, config, y, newton_iters, status%code)
            end if
            
            ! Check for Newton convergence failure
            if (status%code /= BEULER_SUCCESS) then
                status%final_time = t
                return
            end if
            
            t = t + dt
            status%steps_taken = step
            status%newton_iters_total = status%newton_iters_total + newton_iters
            
            ! Store output if needed
            if (config%dense_output .or. mod(step, config%output_every) == 0 &
                .or. abs(t - tf) < epsilon(1.0_wp)) then
                output_idx = output_idx + 1
                if (output_idx <= max_outputs) then
                    t_out(output_idx) = t
                    y_out(:, output_idx) = y
                end if
            end if
            
            ! Check if we've reached the final time
            if (abs(t - tf) < epsilon(1.0_wp)) then
                status%code = BEULER_SUCCESS
                status%final_time = t
                exit
            end if
        end do
        
        ! Trim output arrays
        if (output_idx < max_outputs) then
            t_out = t_out(1:output_idx)
            y_out = y_out(:, 1:output_idx)
        end if
        
        ! Check if max steps exceeded
        if (status%code /= BEULER_SUCCESS) then
            status%code = BEULER_ERR_MAX_STEPS
            status%final_time = t
        end if
    end subroutine beuler_solve
    
    !> Backward Euler step using analytical Jacobian
    subroutine beuler_step_analytical(func, jac, t, y, dt, config, y_new, &
                                     newton_iters, error_code)
        procedure(ode_func_interface) :: func
        procedure(jacobian_interface) :: jac
        real(wp), intent(in) :: t
        real(wp), intent(in) :: y(:)
        real(wp), intent(in) :: dt
        type(beuler_config_t), intent(in) :: config
        real(wp), intent(out) :: y_new(:)
        integer, intent(out) :: newton_iters
        integer, intent(out) :: error_code
        
        integer :: n, iter, info
        real(wp) :: f_new(size(y)), residual(size(y))
        real(wp) :: jacobian(size(y), size(y)), system_matrix(size(y), size(y))
        real(wp) :: delta(size(y))
        real(wp) :: res_norm
        integer, allocatable :: ipiv(:)
        
        n = size(y)
        allocate(ipiv(n))
        
        ! Initial guess: y_new = y (explicit Euler as predictor)
        call func(t + dt, y, f_new)
        y_new = y + dt * f_new
        
        error_code = BEULER_SUCCESS
        
        ! Newton iteration: solve F(y_new) = y_new - y - dt*f(t+dt, y_new) = 0
        do iter = 1, config%max_newton_iter
            ! Evaluate function at current guess
            call func(t + dt, y_new, f_new)
            
            ! Compute residual: F(y_new) = y_new - y - dt*f_new
            residual = y_new - y - dt * f_new
            res_norm = sqrt(sum(residual**2))
            
            ! Check convergence
            if (res_norm < config%newton_tol) then
                newton_iters = iter
                return
            end if
            
            ! Compute Jacobian: J = I - dt*df/dy
            call jac(t + dt, y_new, jacobian)
            system_matrix = -dt * jacobian
            do concurrent (iter = 1:n)
                system_matrix(iter, iter) = 1.0_wp + system_matrix(iter, iter)
            end do
            
            ! Solve: (I - dt*J) * delta = -residual
            delta = -residual
            call dgesv(n, 1, system_matrix, n, ipiv, delta, n, info)
            
            if (info /= 0) then
                error_code = BEULER_ERR_SINGULAR_JACOBIAN
                newton_iters = iter
                return
            end if
            
            ! Update: y_new = y_new + delta
            y_new = y_new + delta
        end do
        
        ! Newton iteration failed to converge
        error_code = BEULER_ERR_NEWTON_FAILED
        newton_iters = config%max_newton_iter
        
        deallocate(ipiv)
    end subroutine beuler_step_analytical
    
    !> Backward Euler step using finite difference Jacobian
    subroutine beuler_step_fd(func, t, y, dt, config, y_new, newton_iters, error_code)
        procedure(ode_func_interface) :: func
        real(wp), intent(in) :: t
        real(wp), intent(in) :: y(:)
        real(wp), intent(in) :: dt
        type(beuler_config_t), intent(in) :: config
        real(wp), intent(out) :: y_new(:)
        integer, intent(out) :: newton_iters
        integer, intent(out) :: error_code
        
        integer :: n, iter, i, info
        real(wp) :: f_new(size(y)), residual(size(y))
        real(wp) :: jacobian(size(y), size(y)), system_matrix(size(y), size(y))
        real(wp) :: delta(size(y)), y_perturbed(size(y))
        real(wp) :: f_perturbed(size(y)), eps
        real(wp) :: res_norm
        integer, allocatable :: ipiv(:)
        
        n = size(y)
        allocate(ipiv(n))
        eps = config%fd_epsilon
        
        ! Initial guess using explicit Euler
        call func(t + dt, y, f_new)
        y_new = y + dt * f_new
        
        error_code = BEULER_SUCCESS
        
        ! Newton iteration
        do iter = 1, config%max_newton_iter
            ! Evaluate function
            call func(t + dt, y_new, f_new)
            
            ! Compute residual
            residual = y_new - y - dt * f_new
            res_norm = sqrt(sum(residual**2))
            
            ! Check convergence
            if (res_norm < config%newton_tol) then
                newton_iters = iter
                deallocate(ipiv)
                return
            end if
            
            ! Compute Jacobian using finite differences
            do i = 1, n
                y_perturbed = y_new
                y_perturbed(i) = y_perturbed(i) + eps
                call func(t + dt, y_perturbed, f_perturbed)
                jacobian(:, i) = (f_perturbed - f_new) / eps
            end do
            
            ! System matrix: I - dt*J
            system_matrix = -dt * jacobian
            do concurrent (i = 1:n)
                system_matrix(i, i) = 1.0_wp + system_matrix(i, i)
            end do
            
            ! Solve linear system
            delta = -residual
            call dgesv(n, 1, system_matrix, n, ipiv, delta, n, info)
            
            if (info /= 0) then
                error_code = BEULER_ERR_SINGULAR_JACOBIAN
                newton_iters = iter
                deallocate(ipiv)
                return
            end if
            
            ! Update solution
            y_new = y_new + delta
        end do
        
        ! Newton failed to converge
        error_code = BEULER_ERR_NEWTON_FAILED
        newton_iters = config%max_newton_iter
        
        deallocate(ipiv)
    end subroutine beuler_step_fd
    
    !> Public interface for single step (uses finite differences)
    subroutine beuler_step(func, t, y, dt, config, y_new)
        procedure(ode_func_interface) :: func
        real(wp), intent(in) :: t
        real(wp), intent(in) :: y(:)
        real(wp), intent(in) :: dt
        type(beuler_config_t), intent(in) :: config
        real(wp), intent(out) :: y_new(:)
        
        integer :: newton_iters, error_code
        
        call beuler_step_fd(func, t, y, dt, config, y_new, newton_iters, error_code)
    end subroutine beuler_step

end module backward_euler