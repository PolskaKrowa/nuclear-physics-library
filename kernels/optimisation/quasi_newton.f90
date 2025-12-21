! fortran/kernels/optimisation/quasi_newton.f90
!
! Quasi-Newton methods for unconstrained optimisation.
! Implements BFGS and L-BFGS algorithms that approximate the Hessian
! using gradient information.
!
! BFGS: Maintains full inverse Hessian approximation
! L-BFGS: Limited-memory variant for large-scale problems
!
! Usage:
!   use quasi_newton
!   type(qn_config_t) :: config
!   config = qn_config_default()
!   call bfgs_minimize(f, grad_f, x0, xmin, config, result)
!
module quasi_newton
    use kinds, only: wp
    use constants, only: TOL_DEFAULT
    use numerics_utils, only: norm2, is_finite, safe_divide
    implicit none
    private
    
    ! Public interface
    public :: qn_config_t, qn_result_t
    public :: qn_config_default
    public :: bfgs_minimize
    public :: lbfgs_minimize
    
    ! Configuration
    type :: qn_config_t
        real(wp) :: tolerance        = 1.0e-6_wp   ! Convergence tolerance
        integer  :: max_iterations   = 1000        ! Maximum iterations
        integer  :: history_size     = 10          ! L-BFGS: number of correction pairs
        real(wp) :: line_search_tol  = 1.0e-4_wp   ! Line search tolerance
        logical  :: use_line_search  = .true.      ! Enable line search
        integer  :: print_interval   = 0           ! Progress printing
    end type qn_config_t
    
    ! Result information
    type :: qn_result_t
        real(wp) :: final_value
        real(wp) :: gradient_norm
        integer  :: iterations
        integer  :: function_evals
        logical  :: converged
        character(len=256) :: message
    end type qn_result_t
    
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
    pure function qn_config_default() result(config)
        type(qn_config_t) :: config
        ! Uses default values from type definition
    end function qn_config_default
    
    !> BFGS quasi-Newton method
    subroutine bfgs_minimize(f, grad_f, x0, xmin, config, result)
        procedure(objective_func) :: f
        procedure(gradient_func) :: grad_f
        real(wp), intent(in) :: x0(:)
        real(wp), intent(out) :: xmin(:)
        type(qn_config_t), intent(in) :: config
        type(qn_result_t), intent(out) :: result
        
        real(wp) :: x(size(x0)), x_new(size(x0))
        real(wp) :: grad(size(x0)), grad_new(size(x0))
        real(wp) :: H(size(x0), size(x0))  ! Inverse Hessian approximation
        real(wp) :: direction(size(x0))
        real(wp) :: s(size(x0)), y(size(x0))
        real(wp) :: grad_norm, alpha, rho, sy
        integer :: iter, n
        
        n = size(x0)
        x = x0
        grad = grad_f(x)
        grad_norm = norm2(grad)
        
        ! Initialise inverse Hessian as identity
        H = 0.0_wp
        do iter = 1, n
            H(iter, iter) = 1.0_wp
        end do
        
        result%iterations = 0
        result%function_evals = 1
        result%converged = .false.
        
        do iter = 1, config%max_iterations
            ! Check convergence
            if (grad_norm < config%tolerance) then
                result%converged = .true.
                result%message = 'Converged: gradient norm below tolerance'
                exit
            end if
            
            ! Compute search direction: d = -H * grad
            direction = -matmul(H, grad)
            
            ! Line search
            if (config%use_line_search) then
                alpha = wolfe_line_search(f, grad_f, x, direction, grad, &
                                         config%line_search_tol, result%function_evals)
            else
                alpha = 1.0_wp
            end if
            
            ! Update position
            x_new = x + alpha * direction
            grad_new = grad_f(x_new)
            
            ! Check for numerical issues
            if (.not. all(is_finite(x_new))) then
                result%message = 'Error: NaN or Inf encountered'
                exit
            end if
            
            ! Compute s and y for BFGS update
            s = x_new - x
            y = grad_new - grad
            sy = sum(s * y)
            
            ! BFGS update of inverse Hessian approximation
            if (sy > 1.0e-10_wp) then
                rho = 1.0_wp / sy
                call bfgs_update_hessian(H, s, y, rho)
            end if
            
            ! Prepare for next iteration
            x = x_new
            grad = grad_new
            grad_norm = norm2(grad)
            
            result%iterations = iter
        end do
        
        xmin = x
        result%final_value = f(xmin)
        result%gradient_norm = grad_norm
        result%function_evals = result%function_evals + 1
        
        if (.not. result%converged) then
            result%message = 'Maximum iterations reached'
        end if
    end subroutine bfgs_minimize
    
    !> Update inverse Hessian using BFGS formula
    pure subroutine bfgs_update_hessian(H, s, y, rho)
        real(wp), intent(inout) :: H(:,:)
        real(wp), intent(in) :: s(:), y(:), rho
        
        real(wp) :: Hy(size(s)), sT_H(size(s))
        real(wp) :: I_minus_rho_sy(size(s), size(s))
        real(wp) :: I(size(s), size(s))
        integer :: i, n
        
        n = size(s)
        
        ! Identity matrix
        I = 0.0_wp
        do i = 1, n
            I(i, i) = 1.0_wp
        end do
        
        ! Compute Hy
        Hy = matmul(H, y)
        
        ! BFGS update: H_new = (I - rho*s*y^T) * H * (I - rho*y*s^T) + rho*s*s^T
        ! Simplified: H_new = H + rho*rho*(y^T*H*y)*s*s^T - rho*(H*y*s^T + s*y^T*H)
        
        ! More numerically stable rank-two update:
        H = H - outer_product(Hy, s) * rho - outer_product(s, Hy) * rho &
            + (1.0_wp + rho * sum(y * Hy)) * outer_product(s, s) * rho
    end subroutine bfgs_update_hessian
    
    !> Limited-memory BFGS method
    subroutine lbfgs_minimize(f, grad_f, x0, xmin, config, result)
        procedure(objective_func) :: f
        procedure(gradient_func) :: grad_f
        real(wp), intent(in) :: x0(:)
        real(wp), intent(out) :: xmin(:)
        type(qn_config_t), intent(in) :: config
        type(qn_result_t), intent(out) :: result
        
        real(wp) :: x(size(x0)), x_new(size(x0))
        real(wp) :: grad(size(x0)), grad_new(size(x0))
        real(wp) :: direction(size(x0)), q(size(x0))
        real(wp), allocatable :: s_hist(:,:), y_hist(:,:), rho_hist(:)
        real(wp), allocatable :: alpha_hist(:)
        real(wp) :: grad_norm, alpha, gamma
        integer :: iter, n, m, k, i, bound
        
        n = size(x0)
        m = config%history_size
        
        ! Allocate history arrays
        allocate(s_hist(n, m), y_hist(n, m), rho_hist(m), alpha_hist(m))
        
        x = x0
        grad = grad_f(x)
        grad_norm = norm2(grad)
        
        result%iterations = 0
        result%function_evals = 1
        result%converged = .false.
        
        k = 0  ! Number of stored corrections
        
        do iter = 1, config%max_iterations
            ! Check convergence
            if (grad_norm < config%tolerance) then
                result%converged = .true.
                result%message = 'Converged: gradient norm below tolerance'
                exit
            end if
            
            ! Compute search direction using L-BFGS two-loop recursion
            q = grad
            bound = min(k, m)
            
            ! First loop
            do i = bound, 1, -1
                alpha_hist(i) = rho_hist(i) * sum(s_hist(:, i) * q)
                q = q - alpha_hist(i) * y_hist(:, i)
            end do
            
            ! Scaling
            if (k > 0) then
                gamma = sum(s_hist(:, 1) * y_hist(:, 1)) / sum(y_hist(:, 1)**2)
            else
                gamma = 1.0_wp
            end if
            
            direction = -gamma * q
            
            ! Second loop
            do i = 1, bound
                alpha = rho_hist(i) * sum(y_hist(:, i) * direction)
                direction = direction + s_hist(:, i) * (alpha_hist(i) - alpha)
            end do
            
            ! Line search
            if (config%use_line_search) then
                alpha = wolfe_line_search(f, grad_f, x, direction, grad, &
                                         config%line_search_tol, result%function_evals)
            else
                alpha = 1.0_wp
            end if
            
            ! Update position
            x_new = x + alpha * direction
            grad_new = grad_f(x_new)
            
            if (.not. all(is_finite(x_new))) then
                result%message = 'Error: NaN or Inf encountered'
                exit
            end if
            
            ! Update history
            call update_lbfgs_history(s_hist, y_hist, rho_hist, &
                                     x_new - x, grad_new - grad, k, m)
            
            ! Prepare for next iteration
            x = x_new
            grad = grad_new
            grad_norm = norm2(grad)
            
            result%iterations = iter
        end do
        
        xmin = x
        result%final_value = f(xmin)
        result%gradient_norm = grad_norm
        result%function_evals = result%function_evals + 1
        
        if (.not. result%converged) then
            result%message = 'Maximum iterations reached'
        end if
        
        deallocate(s_hist, y_hist, rho_hist, alpha_hist)
    end subroutine lbfgs_minimize
    
    !> Update L-BFGS history arrays
    pure subroutine update_lbfgs_history(s_hist, y_hist, rho_hist, s, y, k, m)
        real(wp), intent(inout) :: s_hist(:,:), y_hist(:,:), rho_hist(:)
        real(wp), intent(in) :: s(:), y(:)
        integer, intent(inout) :: k
        integer, intent(in) :: m
        
        real(wp) :: sy
        integer :: i
        
        sy = sum(s * y)
        
        if (sy > 1.0e-10_wp) then
            k = k + 1
            
            if (k > m) then
                ! Shift history to make room
                do i = 1, m - 1
                    s_hist(:, i) = s_hist(:, i + 1)
                    y_hist(:, i) = y_hist(:, i + 1)
                    rho_hist(i) = rho_hist(i + 1)
                end do
                s_hist(:, m) = s
                y_hist(:, m) = y
                rho_hist(m) = 1.0_wp / sy
                k = m
            else
                s_hist(:, k) = s
                y_hist(:, k) = y
                rho_hist(k) = 1.0_wp / sy
            end if
        end if
    end subroutine update_lbfgs_history
    
    !> Wolfe line search
    function wolfe_line_search(f, grad_f, x, direction, grad, &
                                   tol, func_evals) result(alpha)
        procedure(objective_func) :: f
        procedure(gradient_func) :: grad_f
        real(wp), intent(in) :: x(:), direction(:), grad(:), tol
        integer, intent(inout) :: func_evals
        real(wp) :: alpha
        
        real(wp), parameter :: c1 = 1.0e-4_wp
        real(wp), parameter :: c2 = 0.9_wp
        integer, parameter :: max_iter = 20
        
        real(wp) :: alpha_low, alpha_high, f0, f_new, slope0, slope_new
        real(wp) :: x_new(size(x)), grad_new(size(x))
        integer :: i
        
        alpha = 1.0_wp
        alpha_low = 0.0_wp
        alpha_high = 2.0_wp
        
        f0 = f(x)
        slope0 = sum(grad * direction)
        func_evals = func_evals + 1
        
        do i = 1, max_iter
            x_new = x + alpha * direction
            f_new = f(x_new)
            func_evals = func_evals + 1
            
            ! Armijo condition
            if (f_new > f0 + c1 * alpha * slope0) then
                alpha_high = alpha
                alpha = 0.5_wp * (alpha_low + alpha_high)
                cycle
            end if
            
            ! Curvature condition
            grad_new = grad_f(x_new)
            slope_new = sum(grad_new * direction)
            
            if (slope_new >= c2 * slope0) then
                return
            end if
            
            if (slope_new > 0.0_wp) then
                alpha_high = alpha
            else
                alpha_low = alpha
            end if
            
            if (alpha_high - alpha_low < tol) exit
            
            alpha = 0.5_wp * (alpha_low + alpha_high)
        end do
    end function wolfe_line_search
    
    !> Compute outer product of two vectors
    pure function outer_product(a, b) result(C)
        real(wp), intent(in) :: a(:), b(:)
        real(wp) :: C(size(a), size(b))
        integer :: i, j
        
        do j = 1, size(b)
            do i = 1, size(a)
                C(i, j) = a(i) * b(j)
            end do
        end do
    end function outer_product

end module quasi_newton