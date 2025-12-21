! fortran/kernels/optimisation/constrained.f90
!
! Constrained optimisation methods.
! Implements:
! - Projected gradient descent (box constraints)
! - Penalty methods (equality and inequality constraints)
! - Augmented Lagrangian method
! - Active-set methods for quadratic programming
!
! Usage:
!   use constrained
!   type(constrained_config_t) :: config
!   config = constrained_config_default()
!   call projected_gradient_minimize(f, grad_f, x0, lower, upper, xmin, config, result)
!
module constrained
    use kinds, only: wp
    use constants, only: TOL_DEFAULT
    use numerics_utils, only: norm2, is_finite, clip
    implicit none
    private
    
    ! Public interface
    public :: constrained_config_t, constrained_result_t
    public :: constrained_config_default
    public :: projected_gradient_minimize
    public :: penalty_method_minimize
    public :: augmented_lagrangian_minimize
    public :: project_onto_box
    
    ! Configuration
    type :: constrained_config_t
        real(wp) :: learning_rate     = 0.01_wp     ! Step size
        real(wp) :: tolerance         = 1.0e-6_wp   ! Convergence tolerance
        integer  :: max_iterations    = 10000       ! Maximum iterations
        real(wp) :: penalty_param     = 1.0_wp      ! Initial penalty parameter
        real(wp) :: penalty_increase  = 10.0_wp     ! Penalty increase factor
        integer  :: penalty_updates   = 10          ! Number of penalty updates
        logical  :: use_line_search   = .true.      ! Enable line search
        integer  :: print_interval    = 0           ! Progress printing
    end type constrained_config_t
    
    ! Result information
    type :: constrained_result_t
        real(wp) :: final_value
        real(wp) :: gradient_norm
        real(wp) :: constraint_violation
        integer  :: iterations
        integer  :: function_evals
        logical  :: converged
        character(len=256) :: message
    end type constrained_result_t
    
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
            real(wp) :: g(:)
        end function gradient_func
        
        function constraint_func(x) result(c)
            import :: wp
            real(wp), intent(in) :: x(:)
            real(wp), allocatable :: c(:)
        end function constraint_func
        
        function constraint_jacobian_func(x) result(J)
            import :: wp
            real(wp), intent(in) :: x(:)
            real(wp), allocatable :: J(:,:)
        end function constraint_jacobian_func
    end interface
    
contains

    !> Get default configuration
    pure function constrained_config_default() result(config)
        type(constrained_config_t) :: config
        ! Uses default values from type definition
    end function constrained_config_default
    
    !> Projected gradient descent for box constraints: lower <= x <= upper
    subroutine projected_gradient_minimize(f, grad_f, x0, lower, upper, xmin, &
                                          config, result)
        procedure(objective_func) :: f
        procedure(gradient_func) :: grad_f
        real(wp), intent(in) :: x0(:), lower(:), upper(:)
        real(wp), intent(out) :: xmin(:)
        type(constrained_config_t), intent(in) :: config
        type(constrained_result_t), intent(out) :: result
        
        real(wp) :: x(size(x0)), x_new(size(x0)), grad(size(x0))
        real(wp) :: grad_norm, alpha, projected_grad_norm
        integer :: iter, n
        
        n = size(x0)
        x = project_onto_box(x0, lower, upper)
        
        result%iterations = 0
        result%function_evals = 0
        result%converged = .false.
        result%constraint_violation = 0.0_wp
        
        do iter = 1, config%max_iterations
            ! Compute gradient
            grad = grad_f(x)
            grad_norm = norm2(grad)
            
            ! Project gradient (set to zero for active constraints)
            grad = project_gradient(x, grad, lower, upper)
            projected_grad_norm = norm2(grad)
            
            ! Check convergence
            if (projected_grad_norm < config%tolerance) then
                result%converged = .true.
                result%message = 'Converged: projected gradient norm below tolerance'
                exit
            end if
            
            ! Determine step size
            if (config%use_line_search) then
                alpha = projected_line_search(f, x, grad, lower, upper, &
                                             result%function_evals)
            else
                alpha = config%learning_rate
            end if
            
            ! Update with projection
            x_new = x - alpha * grad
            x_new = project_onto_box(x_new, lower, upper)
            
            if (.not. all(is_finite(x_new))) then
                result%message = 'Error: NaN or Inf encountered'
                exit
            end if
            
            x = x_new
            result%iterations = iter
        end do
        
        xmin = x
        result%final_value = f(xmin)
        result%gradient_norm = norm2(grad_f(xmin))
        result%function_evals = result%function_evals + 2
        
        if (.not. result%converged) then
            result%message = 'Maximum iterations reached'
        end if
    end subroutine projected_gradient_minimize
    
    !> Penalty method for general constraints
    subroutine penalty_method_minimize(f, grad_f, constraints, x0, xmin, &
                                      config, result)
        procedure(objective_func) :: f
        procedure(gradient_func) :: grad_f
        procedure(constraint_func) :: constraints
        real(wp), intent(in) :: x0(:)
        real(wp), intent(out) :: xmin(:)
        type(constrained_config_t), intent(in) :: config
        type(constrained_result_t), intent(out) :: result
        
        real(wp) :: x(size(x0)), grad(size(x0)), penalty
        real(wp), allocatable :: c(:)
        real(wp) :: grad_norm
        integer :: iter, k, n, nc
        
        n = size(x0)
        x = x0
        penalty = config%penalty_param
        
        result%iterations = 0
        result%function_evals = 0
        result%converged = .false.
        
        ! Outer loop: increase penalty
        do k = 1, config%penalty_updates
            ! Inner loop: minimise penalised objective
            do iter = 1, config%max_iterations / config%penalty_updates
                ! Compute penalised gradient
                if (.not. allocated(c)) allocate(c(size(constraints(x))))
                c = constraints(x)
                grad = grad_f(x) + penalty * penalty_gradient(x, c, constraints)
                grad_norm = norm2(grad)
                
                ! Check convergence
                result%constraint_violation = norm2(c)
                if (grad_norm < config%tolerance .and. &
                    result%constraint_violation < config%tolerance) then
                    result%converged = .true.
                    result%message = 'Converged: gradient and constraints satisfied'
                    exit
                end if
                
                ! Gradient step
                x = x - config%learning_rate * grad
                
                if (.not. all(is_finite(x))) then
                    result%message = 'Error: NaN or Inf encountered'
                    exit
                end if
                
                result%iterations = result%iterations + 1
            end do
            
            if (result%converged) exit
            
            ! Increase penalty
            penalty = penalty * config%penalty_increase
        end do
        
        xmin = x
        result%final_value = f(xmin)
        result%gradient_norm = norm2(grad_f(xmin))
        result%constraint_violation = norm2(constraints(xmin))
        
        if (.not. result%converged) then
            result%message = 'Maximum iterations reached'
        end if
    end subroutine penalty_method_minimize
    
    !> Augmented Lagrangian method
    subroutine augmented_lagrangian_minimize(f, grad_f, constraints, jac_constraints, &
                                            x0, xmin, config, result)
        procedure(objective_func) :: f
        procedure(gradient_func) :: grad_f
        procedure(constraint_func) :: constraints
        procedure(constraint_jacobian_func) :: jac_constraints
        real(wp), intent(in) :: x0(:)
        real(wp), intent(out) :: xmin(:)
        type(constrained_config_t), intent(in) :: config
        type(constrained_result_t), intent(out) :: result
        
        real(wp) :: x(size(x0)), grad(size(x0))
        real(wp), allocatable :: c(:), lambda(:), J(:,:)
        real(wp) :: grad_norm, penalty
        integer :: iter, k, n, nc, m
        
        n = size(x0)
        x = x0
        
        ! Initialise with one constraint evaluation to get size
        allocate(c, source=constraints(x))
        nc = size(c)
        allocate(lambda(nc), J(nc, n))
        lambda = 0.0_wp
        penalty = config%penalty_param
        
        result%iterations = 0
        result%function_evals = 0
        result%converged = .false.
        
        ! Outer loop: update multipliers
        do k = 1, config%penalty_updates
            ! Inner loop: minimise augmented Lagrangian
            do iter = 1, config%max_iterations / config%penalty_updates
                c = constraints(x)
                J = jac_constraints(x)
                
                ! Augmented Lagrangian gradient
                grad = grad_f(x) + matmul(transpose(J), lambda + penalty * c)
                grad_norm = norm2(grad)
                
                ! Check convergence
                result%constraint_violation = norm2(c)
                if (grad_norm < config%tolerance .and. &
                    result%constraint_violation < config%tolerance) then
                    result%converged = .true.
                    result%message = 'Converged: KKT conditions satisfied'
                    exit
                end if
                
                ! Gradient step
                x = x - config%learning_rate * grad
                
                if (.not. all(is_finite(x))) then
                    result%message = 'Error: NaN or Inf encountered'
                    exit
                end if
                
                result%iterations = result%iterations + 1
            end do
            
            if (result%converged) exit
            
            ! Update Lagrange multipliers
            c = constraints(x)
            lambda = lambda + penalty * c
            
            ! Increase penalty
            penalty = penalty * config%penalty_increase
        end do
        
        xmin = x
        result%final_value = f(xmin)
        result%gradient_norm = norm2(grad_f(xmin))
        result%constraint_violation = norm2(constraints(xmin))
        
        if (.not. result%converged) then
            result%message = 'Maximum iterations reached'
        end if
        
        deallocate(c, lambda, J)
    end subroutine augmented_lagrangian_minimize
    
    !> Project vector onto box constraints
    pure function project_onto_box(x, lower, upper) result(x_proj)
        real(wp), intent(in) :: x(:), lower(:), upper(:)
        real(wp) :: x_proj(size(x))
        integer :: i
        
        do i = 1, size(x)
            x_proj(i) = clip(x(i), lower(i), upper(i))
        end do
    end function project_onto_box
    
    !> Project gradient (zero out components for active constraints)
    pure function project_gradient(x, grad, lower, upper) result(grad_proj)
        real(wp), intent(in) :: x(:), grad(:), lower(:), upper(:)
        real(wp) :: grad_proj(size(grad))
        integer :: i
        real(wp), parameter :: active_tol = 1.0e-10_wp
        
        grad_proj = grad
        
        do i = 1, size(x)
            ! If at lower bound and gradient points downward, project to zero
            if (abs(x(i) - lower(i)) < active_tol .and. grad(i) < 0.0_wp) then
                grad_proj(i) = 0.0_wp
            end if
            
            ! If at upper bound and gradient points upward, project to zero
            if (abs(x(i) - upper(i)) < active_tol .and. grad(i) > 0.0_wp) then
                grad_proj(i) = 0.0_wp
            end if
        end do
    end function project_gradient
    
    !> Line search with projection
    function projected_line_search(f, x, grad, lower, upper, func_evals) result(alpha)
        procedure(objective_func) :: f
        real(wp), intent(in) :: x(:), grad(:), lower(:), upper(:)
        integer, intent(inout) :: func_evals
        real(wp) :: alpha
        
        real(wp), parameter :: c1 = 1.0e-4_wp
        real(wp), parameter :: rho = 0.5_wp
        integer, parameter :: max_iter = 20
        
        real(wp) :: f0, f_new, x_new(size(x))
        integer :: i
        
        alpha = 1.0_wp
        f0 = f(x)
        func_evals = func_evals + 1
        
        do i = 1, max_iter
            x_new = project_onto_box(x - alpha * grad, lower, upper)
            f_new = f(x_new)
            func_evals = func_evals + 1
            
            ! Sufficient decrease condition
            if (f_new <= f0 - c1 * alpha * sum(grad * grad)) then
                return
            end if
            
            alpha = rho * alpha
        end do
    end function projected_line_search
    
    !> Compute gradient of penalty term
    function penalty_gradient(x, c, constraints) result(grad_penalty)
        real(wp), intent(in) :: x(:), c(:)
        procedure(constraint_func) :: constraints
        real(wp) :: grad_penalty(size(x))
        
        real(wp) :: x_plus(size(x))
        real(wp), allocatable :: c_plus(:)
        real(wp), parameter :: h = 1.0e-8_wp
        integer :: i, n
        
        n = size(x)
        grad_penalty = 0.0_wp
        
        ! Finite difference approximation of Jacobian^T * c
        do i = 1, n
            x_plus = x
            x_plus(i) = x_plus(i) + h
            if (.not. allocated(c_plus)) allocate(c_plus, source=constraints(x_plus))
            c_plus = constraints(x_plus)
            grad_penalty(i) = sum((c_plus - c) / h * c)
        end do
    end function penalty_gradient

end module constrained