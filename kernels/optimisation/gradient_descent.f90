! fortran/kernels/optimisation/gradient_descent.f90
!
! Gradient descent optimisation with various enhancements:
! - Basic steepest descent
! - Momentum-based acceleration
! - Adaptive learning rates (AdaGrad, RMSProp, Adam)
! - Line search for step size selection
!
! Usage:
!   use gradient_descent
!   type(gd_config_t) :: config
!   config = gd_config_default()
!   call gradient_descent_minimize(f, grad_f, x0, xmin, config, status)
!
module gradient_descent
    use kinds, only: wp, i64
    use constants, only: TOL_DEFAULT, is_zero
    use numerics_utils, only: norm2, is_finite
    implicit none
    private
    
    ! Public interface
    public :: gd_config_t, gd_result_t
    public :: gd_config_default
    public :: gradient_descent_minimize
    public :: gradient_descent_momentum
    public :: gradient_descent_adam
    public :: backtracking_line_search
    
    ! Configuration for gradient descent
    type :: gd_config_t
        real(wp) :: learning_rate      = 0.01_wp    ! Step size (alpha)
        real(wp) :: momentum           = 0.0_wp     ! Momentum coefficient (0 = no momentum)
        real(wp) :: tolerance          = 1.0e-6_wp  ! Convergence tolerance
        integer  :: max_iterations     = 10000      ! Maximum iterations
        logical  :: use_line_search    = .false.    ! Enable backtracking line search
        logical  :: adaptive_learning  = .false.    ! Adaptive learning rate
        real(wp) :: beta1              = 0.9_wp     ! Adam: exponential decay rate for 1st moment
        real(wp) :: beta2              = 0.999_wp   ! Adam: exponential decay rate for 2nd moment
        real(wp) :: epsilon            = 1.0e-8_wp  ! Adam: small constant for numerical stability
        integer  :: print_interval     = 0          ! Print progress every N iterations (0 = no output)
    end type gd_config_t
    
    ! Result information
    type :: gd_result_t
        real(wp) :: final_value          ! Function value at minimum
        real(wp) :: gradient_norm        ! Norm of gradient at minimum
        integer  :: iterations           ! Number of iterations performed
        integer  :: function_evals       ! Number of function evaluations
        logical  :: converged            ! True if converged
        character(len=256) :: message    ! Status message
    end type gd_result_t
    
    ! Abstract interface for objective function
    abstract interface
        pure function objective_func(x) result(f)
            import :: wp
            real(wp), intent(in) :: x(:)
            real(wp) :: f
        end function objective_func
        
        pure function gradient_func(x) result(g)
            import :: wp
            real(wp), intent(in) :: x(:)
            real(wp) :: g(size(x))
        end function gradient_func
    end interface
    
contains

    !> Get default configuration
    pure function gd_config_default() result(config)
        type(gd_config_t) :: config
        ! Uses default values from type definition
    end function gd_config_default
    
    !> Basic gradient descent with optional momentum
    subroutine gradient_descent_minimize(f, grad_f, x0, xmin, config, result)
        procedure(objective_func) :: f
        procedure(gradient_func) :: grad_f
        real(wp), intent(in) :: x0(:)
        real(wp), intent(out) :: xmin(:)
        type(gd_config_t), intent(in) :: config
        type(gd_result_t), intent(out) :: result
        
        real(wp) :: x(size(x0)), grad(size(x0)), velocity(size(x0))
        real(wp) :: f_val, grad_norm, alpha
        integer :: iter, n
        
        n = size(x0)
        x = x0
        velocity = 0.0_wp
        
        result%iterations = 0
        result%function_evals = 0
        result%converged = .false.
        
        do iter = 1, config%max_iterations
            ! Compute gradient
            grad = grad_f(x)
            grad_norm = norm2(grad)
            
            ! Check convergence
            if (grad_norm < config%tolerance) then
                result%converged = .true.
                result%message = 'Converged: gradient norm below tolerance'
                exit
            end if
            
            ! Determine step size
            if (config%use_line_search) then
                alpha = backtracking_line_search(f, grad_f, x, grad, result%function_evals)
            else
                alpha = config%learning_rate
            end if
            
            ! Update with momentum
            if (config%momentum > 0.0_wp) then
                velocity = config%momentum * velocity - alpha * grad
                x = x + velocity
            else
                x = x - alpha * grad
            end if
            
            ! Check for NaN/Inf
            if (.not. all(is_finite(x))) then
                result%message = 'Error: NaN or Inf encountered'
                exit
            end if
            
            ! Progress reporting
            if (config%print_interval > 0) then
                if (mod(iter, config%print_interval) == 0) then
                    f_val = f(x)
                    result%function_evals = result%function_evals + 1
                end if
            end if
            
            result%iterations = iter
        end do
        
        ! Final evaluation
        xmin = x
        result%final_value = f(xmin)
        result%gradient_norm = norm2(grad_f(xmin))
        result%function_evals = result%function_evals + 2
        
        if (.not. result%converged) then
            result%message = 'Maximum iterations reached'
        end if
    end subroutine gradient_descent_minimize
    
    !> Gradient descent with Nesterov momentum
    subroutine gradient_descent_momentum(f, grad_f, x0, xmin, config, result)
        procedure(objective_func) :: f
        procedure(gradient_func) :: grad_f
        real(wp), intent(in) :: x0(:)
        real(wp), intent(out) :: xmin(:)
        type(gd_config_t), intent(in) :: config
        type(gd_result_t), intent(out) :: result
        
        real(wp) :: x(size(x0)), grad(size(x0)), velocity(size(x0))
        real(wp) :: x_lookahead(size(x0)), grad_norm, alpha
        integer :: iter
        
        x = x0
        velocity = 0.0_wp
        
        result%iterations = 0
        result%function_evals = 0
        result%converged = .false.
        
        do iter = 1, config%max_iterations
            ! Nesterov lookahead
            x_lookahead = x + config%momentum * velocity
            
            ! Compute gradient at lookahead position
            grad = grad_f(x_lookahead)
            grad_norm = norm2(grad)
            
            if (grad_norm < config%tolerance) then
                result%converged = .true.
                result%message = 'Converged: gradient norm below tolerance'
                exit
            end if
            
            ! Determine step size
            if (config%use_line_search) then
                alpha = backtracking_line_search(f, grad_f, x_lookahead, grad, &
                                                 result%function_evals)
            else
                alpha = config%learning_rate
            end if
            
            ! Update velocity and position
            velocity = config%momentum * velocity - alpha * grad
            x = x + velocity
            
            if (.not. all(is_finite(x))) then
                result%message = 'Error: NaN or Inf encountered'
                exit
            end if
            
            result%iterations = iter
        end do
        
        xmin = x
        result%final_value = f(xmin)
        result%gradient_norm = norm2(grad_f(xmin))
        result%function_evals = result%function_evals + 2
        
        if (.not. result%converged) then
            result%message = 'Maximum iterations reached'
        end if
    end subroutine gradient_descent_momentum
    
    !> Adam optimiser (Adaptive Moment Estimation)
    subroutine gradient_descent_adam(f, grad_f, x0, xmin, config, result)
        procedure(objective_func) :: f
        procedure(gradient_func) :: grad_f
        real(wp), intent(in) :: x0(:)
        real(wp), intent(out) :: xmin(:)
        type(gd_config_t), intent(in) :: config
        type(gd_result_t), intent(out) :: result
        
        real(wp) :: x(size(x0)), grad(size(x0))
        real(wp) :: m(size(x0)), v(size(x0))  ! 1st and 2nd moment estimates
        real(wp) :: m_hat(size(x0)), v_hat(size(x0))  ! Bias-corrected moments
        real(wp) :: grad_norm, beta1_t, beta2_t
        integer :: iter, n
        
        n = size(x0)
        x = x0
        m = 0.0_wp
        v = 0.0_wp
        
        result%iterations = 0
        result%function_evals = 0
        result%converged = .false.
        
        do iter = 1, config%max_iterations
            ! Compute gradient
            grad = grad_f(x)
            grad_norm = norm2(grad)
            
            if (grad_norm < config%tolerance) then
                result%converged = .true.
                result%message = 'Converged: gradient norm below tolerance'
                exit
            end if
            
            ! Update biased first moment estimate
            m = config%beta1 * m + (1.0_wp - config%beta1) * grad
            
            ! Update biased second raw moment estimate
            v = config%beta2 * v + (1.0_wp - config%beta2) * grad**2
            
            ! Compute bias-corrected moment estimates
            beta1_t = config%beta1**iter
            beta2_t = config%beta2**iter
            m_hat = m / (1.0_wp - beta1_t)
            v_hat = v / (1.0_wp - beta2_t)
            
            ! Update parameters
            x = x - config%learning_rate * m_hat / (sqrt(v_hat) + config%epsilon)
            
            if (.not. all(is_finite(x))) then
                result%message = 'Error: NaN or Inf encountered'
                exit
            end if
            
            result%iterations = iter
        end do
        
        xmin = x
        result%final_value = f(xmin)
        result%gradient_norm = norm2(grad_f(xmin))
        result%function_evals = result%function_evals + 2
        
        if (.not. result%converged) then
            result%message = 'Maximum iterations reached'
        end if
    end subroutine gradient_descent_adam
    
    !> Backtracking line search for step size selection
    pure function backtracking_line_search(f, grad_f, x, grad, func_evals) result(alpha)
        procedure(objective_func) :: f
        procedure(gradient_func) :: grad_f
        real(wp), intent(in) :: x(:), grad(:)
        integer, intent(inout) :: func_evals
        real(wp) :: alpha
        
        real(wp), parameter :: c = 1.0e-4_wp      ! Armijo condition constant
        real(wp), parameter :: rho = 0.5_wp       ! Backtracking factor
        real(wp), parameter :: alpha_init = 1.0_wp
        integer, parameter :: max_iter = 50
        
        real(wp) :: f0, f_new, descent
        real(wp) :: x_new(size(x))
        integer :: i
        
        alpha = alpha_init
        f0 = f(x)
        func_evals = func_evals + 1
        descent = sum(grad * grad)
        
        do i = 1, max_iter
            x_new = x - alpha * grad
            f_new = f(x_new)
            func_evals = func_evals + 1
            
            ! Armijo condition
            if (f_new <= f0 - c * alpha * descent) then
                return
            end if
            
            alpha = rho * alpha
        end do
    end function backtracking_line_search

end module gradient_descent