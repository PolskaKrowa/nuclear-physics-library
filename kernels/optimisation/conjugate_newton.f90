! fortran/kernels/optimisation/conjugate_gradient.f90
!
! Nonlinear conjugate gradient methods for unconstrained optimisation.
! Implements several update formulas:
! - Fletcher-Reeves
! - Polak-Ribière
! - Hestenes-Stiefel
! - Dai-Yuan
!
! Also includes linear CG for solving Ax = b with symmetric positive definite A.
!
! Usage:
!   use conjugate_gradient
!   type(cg_config_t) :: config
!   config = cg_config_default()
!   call conjugate_gradient_minimize(f, grad_f, x0, xmin, config, result)
!
module conjugate_gradient
    use kinds, only: wp
    use constants, only: TOL_DEFAULT
    use numerics_utils, only: norm2, is_finite, safe_divide
    implicit none
    private
    
    ! Public interface
    public :: cg_config_t, cg_result_t, cg_method_t
    public :: cg_config_default
    public :: conjugate_gradient_minimize
    public :: conjugate_gradient_linear
    
    ! CG method variants
    type :: cg_method_t
        integer :: id
        character(len=32) :: name
    end type cg_method_t
    
    ! Method identifiers
    integer, parameter, public :: CG_FLETCHER_REEVES = 1
    integer, parameter, public :: CG_POLAK_RIBIERE = 2
    integer, parameter, public :: CG_HESTENES_STIEFEL = 3
    integer, parameter, public :: CG_DAI_YUAN = 4
    
    ! Configuration
    type :: cg_config_t
        integer  :: method             = CG_POLAK_RIBIERE  ! CG variant
        real(wp) :: tolerance          = 1.0e-6_wp          ! Convergence tolerance
        integer  :: max_iterations     = 10000              ! Maximum iterations
        integer  :: restart_interval   = 50                 ! Restart every N iterations
        logical  :: use_line_search    = .true.             ! Enable line search
        real(wp) :: line_search_tol    = 1.0e-4_wp          ! Line search tolerance
        integer  :: print_interval     = 0                  ! Progress printing
    end type cg_config_t
    
    ! Result information
    type :: cg_result_t
        real(wp) :: final_value
        real(wp) :: gradient_norm
        integer  :: iterations
        integer  :: function_evals
        integer  :: restarts
        logical  :: converged
        character(len=256) :: message
    end type cg_result_t
    
    ! Abstract interfaces
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
    pure function cg_config_default() result(config)
        type(cg_config_t) :: config
        ! Uses default values from type definition
    end function cg_config_default
    
    !> Nonlinear conjugate gradient minimisation
    subroutine conjugate_gradient_minimize(f, grad_f, x0, xmin, config, result)
        procedure(objective_func) :: f
        procedure(gradient_func) :: grad_f
        real(wp), intent(in) :: x0(:)
        real(wp), intent(out) :: xmin(:)
        type(cg_config_t), intent(in) :: config
        type(cg_result_t), intent(out) :: result
        
        real(wp) :: x(size(x0)), x_new(size(x0))
        real(wp) :: grad(size(x0)), grad_new(size(x0))
        real(wp) :: direction(size(x0))
        real(wp) :: grad_norm, grad_norm_new, beta, alpha
        integer :: iter, n
        
        n = size(x0)
        x = x0
        grad = grad_f(x)
        direction = -grad
        grad_norm = norm2(grad)
        
        result%iterations = 0
        result%function_evals = 1
        result%restarts = 0
        result%converged = .false.
        
        do iter = 1, config%max_iterations
            ! Check convergence
            if (grad_norm < config%tolerance) then
                result%converged = .true.
                result%message = 'Converged: gradient norm below tolerance'
                exit
            end if
            
            ! Line search along direction
            if (config%use_line_search) then
                alpha = strong_wolfe_line_search(f, grad_f, x, direction, grad, &
                                                 config%line_search_tol, result%function_evals)
            else
                alpha = 1.0_wp
            end if
            
            ! Update position
            x_new = x + alpha * direction
            
            ! Compute new gradient
            grad_new = grad_f(x_new)
            grad_norm_new = norm2(grad_new)
            
            ! Check for numerical issues
            if (.not. all(is_finite(x_new))) then
                result%message = 'Error: NaN or Inf encountered'
                exit
            end if
            
            ! Compute beta using selected method
            beta = compute_beta(grad, grad_new, direction, config%method)
            
            ! Restart if beta is negative or at restart interval
            if (beta < 0.0_wp .or. mod(iter, config%restart_interval) == 0) then
                direction = -grad_new
                result%restarts = result%restarts + 1
            else
                ! Update search direction
                direction = -grad_new + beta * direction
            end if
            
            ! Prepare for next iteration
            x = x_new
            grad = grad_new
            grad_norm = grad_norm_new
            
            result%iterations = iter
        end do
        
        xmin = x
        result%final_value = f(xmin)
        result%gradient_norm = grad_norm
        result%function_evals = result%function_evals + 1
        
        if (.not. result%converged) then
            result%message = 'Maximum iterations reached'
        end if
    end subroutine conjugate_gradient_minimize
    
    !> Compute beta parameter based on selected method
    pure function compute_beta(grad_old, grad_new, direction, method) result(beta)
        real(wp), intent(in) :: grad_old(:), grad_new(:), direction(:)
        integer, intent(in) :: method
        real(wp) :: beta
        
        real(wp) :: grad_old_norm2, grad_new_norm2, y_dot_y, y_dot_d, y_dot_gnew
        real(wp) :: y(size(grad_old))
        
        y = grad_new - grad_old
        
        select case(method)
        case(CG_FLETCHER_REEVES)
            ! beta = ||g_new||^2 / ||g_old||^2
            grad_old_norm2 = sum(grad_old**2)
            grad_new_norm2 = sum(grad_new**2)
            beta = safe_divide(grad_new_norm2, grad_old_norm2)
            
        case(CG_POLAK_RIBIERE)
            ! beta = g_new^T * (g_new - g_old) / ||g_old||^2
            grad_old_norm2 = sum(grad_old**2)
            y_dot_gnew = sum(y * grad_new)
            beta = safe_divide(y_dot_gnew, grad_old_norm2)
            
        case(CG_HESTENES_STIEFEL)
            ! beta = g_new^T * (g_new - g_old) / d^T * (g_new - g_old)
            y_dot_gnew = sum(y * grad_new)
            y_dot_d = sum(y * direction)
            beta = safe_divide(y_dot_gnew, y_dot_d)
            
        case(CG_DAI_YUAN)
            ! beta = ||g_new||^2 / (d^T * (g_new - g_old))
            grad_new_norm2 = sum(grad_new**2)
            y_dot_d = sum(y * direction)
            beta = safe_divide(grad_new_norm2, y_dot_d)
            
        case default
            ! Default to Polak-Ribière
            grad_old_norm2 = sum(grad_old**2)
            y_dot_gnew = sum(y * grad_new)
            beta = safe_divide(y_dot_gnew, grad_old_norm2)
        end select
        
        ! Ensure beta is non-negative (Polak-Ribière+ variant)
        beta = max(beta, 0.0_wp)
    end function compute_beta
    
    !> Strong Wolfe line search
    pure function strong_wolfe_line_search(f, grad_f, x, direction, grad, &
                                          tol, func_evals) result(alpha)
        procedure(objective_func) :: f
        procedure(gradient_func) :: grad_f
        real(wp), intent(in) :: x(:), direction(:), grad(:), tol
        integer, intent(inout) :: func_evals
        real(wp) :: alpha
        
        real(wp), parameter :: c1 = 1.0e-4_wp
        real(wp), parameter :: c2 = 0.9_wp
        real(wp), parameter :: alpha_max = 2.0_wp
        integer, parameter :: max_iter = 20
        
        real(wp) :: alpha_low, alpha_high, f0, f_new, slope0, slope_new
        real(wp) :: x_new(size(x)), grad_new(size(x))
        integer :: i
        
        alpha = 1.0_wp
        alpha_low = 0.0_wp
        alpha_high = alpha_max
        
        f0 = f(x)
        slope0 = sum(grad * direction)
        func_evals = func_evals + 1
        
        do i = 1, max_iter
            x_new = x + alpha * direction
            f_new = f(x_new)
            func_evals = func_evals + 1
            
            ! Check Armijo condition
            if (f_new > f0 + c1 * alpha * slope0) then
                alpha_high = alpha
                alpha = 0.5_wp * (alpha_low + alpha_high)
                cycle
            end if
            
            ! Check curvature condition
            grad_new = grad_f(x_new)
            slope_new = sum(grad_new * direction)
            
            if (abs(slope_new) <= -c2 * slope0) then
                return  ! Strong Wolfe conditions satisfied
            end if
            
            if (slope_new >= 0.0_wp) then
                alpha_high = alpha
            else
                alpha_low = alpha
            end if
            
            alpha = 0.5_wp * (alpha_low + alpha_high)
        end do
    end function strong_wolfe_line_search
    
    !> Linear conjugate gradient for solving Ax = b
    !! A must be symmetric positive definite
    subroutine conjugate_gradient_linear(A, b, x0, x, config, result)
        real(wp), intent(in) :: A(:,:), b(:)
        real(wp), intent(in) :: x0(:)
        real(wp), intent(out) :: x(:)
        type(cg_config_t), intent(in) :: config
        type(cg_result_t), intent(out) :: result
        
        real(wp) :: r(size(b)), r_new(size(b))
        real(wp) :: p(size(b)), Ap(size(b))
        real(wp) :: alpha, beta, r_dot_r, r_dot_r_new
        integer :: iter, n
        
        n = size(b)
        x = x0
        
        ! Initial residual: r = b - Ax
        r = b - matmul(A, x)
        p = r
        r_dot_r = sum(r**2)
        
        result%iterations = 0
        result%converged = .false.
        
        do iter = 1, config%max_iterations
            ! Check convergence
            if (sqrt(r_dot_r) < config%tolerance) then
                result%converged = .true.
                result%message = 'Converged: residual norm below tolerance'
                exit
            end if
            
            ! Compute Ap
            Ap = matmul(A, p)
            
            ! Step size
            alpha = r_dot_r / sum(p * Ap)
            
            ! Update solution and residual
            x = x + alpha * p
            r_new = r - alpha * Ap
            r_dot_r_new = sum(r_new**2)
            
            ! Beta coefficient
            beta = r_dot_r_new / r_dot_r
            
            ! Update search direction
            p = r_new + beta * p
            
            ! Prepare for next iteration
            r = r_new
            r_dot_r = r_dot_r_new
            
            result%iterations = iter
        end do
        
        result%gradient_norm = sqrt(r_dot_r)
        
        if (.not. result%converged) then
            result%message = 'Maximum iterations reached'
        end if
    end subroutine conjugate_gradient_linear

end module conjugate_gradient