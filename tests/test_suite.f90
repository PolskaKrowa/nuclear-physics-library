! test_suite.f90
!
! Comprehensive test suite for nuclear physics simulation kernels and models
! Systematically tests each component and reports results
!
program test_suite
    use kinds
    use constants
    use numerics_utils
    use rng
    use, intrinsic :: ieee_arithmetic
    implicit none
    
    integer :: total_tests = 0
    integer :: passed_tests = 0
    integer :: failed_tests = 0
    
    print *, '========================================='
    print *, 'NUCLEAR PHYSICS SIMULATION TEST SUITE'
    print *, '========================================='
    print *, ''
    
    ! Core modules
    call test_kinds_module()
    call test_constants_module()
    call test_numerics_utils_module()
    call test_rng_module()
    
    ! Linear algebra kernels
    call test_dense_matrix_module()
    call test_eigen_module()
    call test_solve_linear_module()
    call test_sparse_matrix_module()
    
    ! ODE solvers
    call test_rk4_module()
    call test_dormand_prince_module()
    call test_backward_euler_module()
    
    ! Optimization
    call test_gradient_descent_module()
    call test_conjugate_gradient_module()
    call test_quasi_newton_module()
    call test_constrained_module()
    
    ! PDE solvers
    call test_finite_difference_module()
    call test_finite_volume_module()
    call test_spectral_module()
    
    ! Physics models
    call test_heat_transfer_module()
    call test_fluid_dynamics_module()
    call test_nuclear_fission_module()
    call test_pressure_dynamics_module()
    
    ! Summary
    print *, ''
    print *, '========================================='
    print *, 'TEST SUMMARY'
    print *, '========================================='
    print *, 'Total tests:  ', total_tests
    print *, 'Passed:       ', passed_tests
    print *, 'Failed:       ', failed_tests
    print *, 'Success rate: ', real(passed_tests)/real(total_tests)*100.0, '%'
    print *, '========================================='

    ! Exit with number of failed tests (0 if all passed)
    if (failed_tests > 0) then
        stop 1
    else
        stop 0
    end if
    
contains

    !===========================================================================
    ! CORE MODULE TESTS
    !===========================================================================
    
    subroutine test_kinds_module()
        print *, ''
        print *, '--- Testing KINDS module ---'
        
        call run_test('Integer kinds defined', test_integer_kinds())
        call run_test('Real kinds defined', test_real_kinds())
        call run_test('Working precision correct', test_working_precision())
    end subroutine
    
    function test_integer_kinds() result(pass)
        logical :: pass
        pass = (i8 == 1) .and. (i16 == 2) .and. (i32 == 4) .and. (i64 == 8)
    end function
    
    function test_real_kinds() result(pass)
        logical :: pass
        pass = (sp == 4) .and. (dp == 8) .and. (qp == 16)
    end function
    
    function test_working_precision() result(pass)
        logical :: pass
        pass = (wp == dp)
    end function
    
    subroutine test_constants_module()
        print *, ''
        print *, '--- Testing CONSTANTS module ---'
        
        call run_test('PI value correct', test_pi_value())
        call run_test('Physical constants defined', test_physical_constants())
        call run_test('nearly_equal works', test_nearly_equal())
        call run_test('is_zero works', test_is_zero())
    end subroutine
    
    function test_pi_value() result(pass)
        logical :: pass
        pass = abs(PI - 3.141592653589793_wp) < 1.0e-14_wp
    end function
    
    function test_physical_constants() result(pass)
        logical :: pass
        pass = (C_LIGHT > 0.0_wp) .and. (H_PLANCK > 0.0_wp) .and. &
               (K_BOLTZMANN > 0.0_wp)
    end function
    
    function test_nearly_equal() result(pass)
        logical :: pass
        real(wp) :: a, b
        a = 1.0_wp
        b = 1.0_wp + 1.0e-13_wp
        pass = nearly_equal(a, b)
    end function
    
    function test_is_zero() result(pass)
        logical :: pass
        pass = is_zero(1.0e-20_wp) .and. (.not. is_zero(1.0e-5_wp))
    end function
    
    subroutine test_numerics_utils_module()
        print *, ''
        print *, '--- Testing NUMERICS_UTILS module ---'
        
        call run_test('linspace generates correct array', test_linspace())
        call run_test('clip works correctly', test_clip())
        call run_test('norm2 computes correctly', test_norm2())
        call run_test('safe_divide handles zero', test_safe_divide())
        call run_test('is_finite detects NaN/Inf', test_is_finite())
    end subroutine
    
    function test_linspace() result(pass)
        logical :: pass
        real(wp) :: arr(5)
        arr = linspace(0.0_wp, 1.0_wp, 5)
        pass = abs(arr(1) - 0.0_wp) < TOL_DEFAULT .and. &
               abs(arr(5) - 1.0_wp) < TOL_DEFAULT .and. &
               abs(arr(3) - 0.5_wp) < TOL_DEFAULT
    end function
    
    function test_clip() result(pass)
        logical :: pass
        pass = (clip(5.0_wp, 0.0_wp, 1.0_wp) == 1.0_wp) .and. &
               (clip(-1.0_wp, 0.0_wp, 1.0_wp) == 0.0_wp) .and. &
               (clip(0.5_wp, 0.0_wp, 1.0_wp) == 0.5_wp)
    end function
    
    function test_norm2() result(pass)
        logical :: pass
        real(wp) :: v(3)
        v = [3.0_wp, 4.0_wp, 0.0_wp]
        pass = abs(norm2(v) - 5.0_wp) < TOL_DEFAULT
    end function
    
    function test_safe_divide() result(pass)
        logical :: pass
        real(wp) :: result
        result = safe_divide(1.0_wp, 0.0_wp)
        pass = (result == 0.0_wp)
    end function
    
    function test_is_finite() result(pass)
        logical :: pass
        real(wp) :: x, y
        x = 1.0_wp
        y = ieee_value(y, ieee_quiet_nan)
        pass = is_finite(x) .and. (.not. is_finite(y))
    end function
    
    subroutine test_rng_module()
        print *, ''
        print *, '--- Testing RNG module ---'
        
        call run_test('RNG seeding works', test_rng_seeding())
        call run_test('Uniform distribution in [0,1)', test_rng_uniform())
        call run_test('Normal distribution', test_rng_normal())
    end subroutine
    
    function test_rng_seeding() result(pass)
        logical :: pass
        real(wp) :: r1, r2
        
        call rng_seed(12345_i64)
        r1 = rng_uniform()
        
        call rng_seed(12345_i64)
        r2 = rng_uniform()
        
        pass = (r1 == r2)  ! Same seed should give same result
    end function
    
    function test_rng_uniform() result(pass)
        logical :: pass
        real(wp) :: r
        integer :: i
        
        call rng_seed(54321_i64)
        pass = .true.
        
        do i = 1, 100
            r = rng_uniform()
            if (r < 0.0_wp .or. r >= 1.0_wp) then
                pass = .false.
                exit
            end if
        end do
    end function
    
    function test_rng_normal() result(pass)
        logical :: pass
        real(wp) :: r
        
        call rng_seed(99999_i64)
        r = rng_normal(mean=0.0_wp, sigma=1.0_wp)
        
        ! Just check it's finite
        pass = is_finite(r)
    end function
    
    !===========================================================================
    ! LINEAR ALGEBRA TESTS
    !===========================================================================
    
    subroutine test_dense_matrix_module()
        use dense_matrix
        print *, ''
        print *, '--- Testing DENSE_MATRIX module ---'
        
        call run_test('Matrix-vector multiply', test_matvec())
        call run_test('Matrix transpose', test_transpose())
        call run_test('Matrix norm', test_matrix_norm())
        call run_test('Matrix determinant', test_determinant())
    end subroutine
    
    function test_matvec() result(pass)
        use dense_matrix
        logical :: pass
        real(wp) :: A(3,3), x(3), y(3), expected(3)
        
        ! Identity matrix test
        A = reshape([1.0_wp, 0.0_wp, 0.0_wp, &
                     0.0_wp, 1.0_wp, 0.0_wp, &
                     0.0_wp, 0.0_wp, 1.0_wp], [3,3])
        x = [1.0_wp, 2.0_wp, 3.0_wp]
        expected = x
        
        call matrix_vector_mult(A, x, y)
        
        pass = maxval(abs(y - expected)) < TOL_DEFAULT
    end function
    
    function test_transpose() result(pass)
        use dense_matrix
        logical :: pass
        real(wp) :: A(2,3), At(3,2)
        
        A = reshape([1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, 6.0_wp], [2,3])
        At = matrix_transpose(A)
        
        pass = (At(1,1) == 1.0_wp) .and. (At(1,2) == 2.0_wp) .and. &
               (At(2,1) == 3.0_wp) .and. (At(3,2) == 6.0_wp)
    end function
    
    function test_matrix_norm() result(pass)
        use dense_matrix
        logical :: pass
        real(wp) :: A(2,2), nrm
        
        A = reshape([3.0_wp, 4.0_wp, 0.0_wp, 0.0_wp], [2,2])
        nrm = matrix_norm(A, 'F')  ! Frobenius norm
        
        pass = abs(nrm - 5.0_wp) < TOL_DEFAULT
    end function
    
    function test_determinant() result(pass)
        use dense_matrix
        logical :: pass
        real(wp) :: A(2,2), det
        integer(i32) :: status
        
        A = reshape([1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp], [2,2])
        det = matrix_determinant(A, status)
        
        pass = (abs(det - (-2.0_wp)) < TOL_DEFAULT) .and. (status == MAT_SUCCESS)
    end function
    
    subroutine test_eigen_module()
        use eigen
        print *, ''
        print *, '--- Testing EIGEN module ---'
        
        call run_test('Symmetric eigenvalues', test_symmetric_eigen())
    end subroutine
    
    function test_symmetric_eigen() result(pass)
        use eigen
        logical :: pass
        real(wp) :: A(2,2), eigenvalues(2)
        integer(i32) :: status
        
        ! Simple symmetric matrix with known eigenvalues
        A = reshape([2.0_wp, 1.0_wp, 1.0_wp, 2.0_wp], [2,2])
        
        call eigen_symmetric(A, eigenvalues, status=status)
        
        ! Eigenvalues should be 1 and 3
        pass = (status == EIGEN_SUCCESS) .and. &
               (abs(eigenvalues(1) - 1.0_wp) < 1.0e-6_wp) .and. &
               (abs(eigenvalues(2) - 3.0_wp) < 1.0e-6_wp)
    end function
    
    subroutine test_solve_linear_module()
        use solve_linear
        print *, ''
        print *, '--- Testing SOLVE_LINEAR module ---'
        
        call run_test('Dense system solver', test_dense_solver())
        call run_test('Tridiagonal solver', test_tridiagonal_solver())
    end subroutine
    
    function test_dense_solver() result(pass)
        use solve_linear
        logical :: pass
        real(wp) :: A(2,2), b(2), x(2), expected(2)
        integer(i32) :: status
        
        A = reshape([2.0_wp, 1.0_wp, 1.0_wp, 2.0_wp], [2,2])
        b = [3.0_wp, 3.0_wp]
        expected = [1.0_wp, 1.0_wp]
        
        call solve_dense(A, b, x, status)
        
        pass = (status == SOLVE_SUCCESS) .and. &
               (maxval(abs(x - expected)) < 1.0e-6_wp)
    end function
    
    function test_tridiagonal_solver() result(pass)
        use solve_linear
        logical :: pass
        real(wp) :: a(2), b(3), c(2), d(3), x(3)
        integer(i32) :: status
        
        ! Simple tridiagonal system
        a = [1.0_wp, 1.0_wp]
        b = [2.0_wp, 2.0_wp, 2.0_wp]
        c = [1.0_wp, 1.0_wp]
        d = [1.0_wp, 2.0_wp, 1.0_wp]
        
        call solve_tridiagonal(a, b, c, d, x, status)
        
        pass = (status == SOLVE_SUCCESS) .and. all(is_finite(x))
    end function
    
    subroutine test_sparse_matrix_module()
        use sparse_matrix
        print *, ''
        print *, '--- Testing SPARSE_MATRIX module ---'
        
        call run_test('Sparse matrix creation', test_sparse_creation())
        call run_test('Sparse matrix-vector multiply', test_sparse_matvec())
    end subroutine
    
    function test_sparse_creation() result(pass)
        use sparse_matrix
        logical :: pass
        real(wp) :: A_dense(3,3)
        type(sparse_matrix_t) :: A_sparse
        integer(i32) :: status
        
        A_dense = reshape([1.0_wp, 0.0_wp, 2.0_wp, &
                          0.0_wp, 3.0_wp, 0.0_wp, &
                          4.0_wp, 0.0_wp, 5.0_wp], [3,3])
        
        call sparse_create_from_dense(A_dense, A_sparse, status=status)
        
        pass = (status == SPARSE_SUCCESS) .and. (A_sparse%nnz == 5)
        
        call sparse_destroy(A_sparse)
    end function
    
    function test_sparse_matvec() result(pass)
        use sparse_matrix
        logical :: pass
        real(wp) :: A_dense(2,2), x(2), y(2), y_expected(2)
        type(sparse_matrix_t) :: A_sparse
        integer(i32) :: status
        
        A_dense = reshape([2.0_wp, 0.0_wp, 0.0_wp, 3.0_wp], [2,2])
        x = [1.0_wp, 2.0_wp]
        y_expected = [2.0_wp, 6.0_wp]
        
        call sparse_create_from_dense(A_dense, A_sparse, status=status)
        call sparse_matvec(A_sparse, x, y)
        
        pass = maxval(abs(y - y_expected)) < TOL_DEFAULT
        
        call sparse_destroy(A_sparse)
    end function
    
    !===========================================================================
    ! ODE SOLVER TESTS
    !===========================================================================
    
    subroutine test_rk4_module()
        use rk4
        print *, ''
        print *, '--- Testing RK4 module ---'
        
        call run_test('RK4 solves simple ODE', test_rk4_simple())
    end subroutine
    
    function test_rk4_simple() result(pass)
        use rk4
        logical :: pass
        type(rk4_config_t) :: config
        type(rk4_status_t) :: status
        real(wp) :: y0(1), t_span(2)
        real(wp), allocatable :: t_out(:), y_out(:,:)
        
        ! Solve dy/dt = -y, y(0) = 1
        ! Exact solution: y = exp(-t)
        
        y0 = [1.0_wp]
        t_span = [0.0_wp, 1.0_wp]
        config%dt = 0.01_wp
        config%max_steps = 1000
        
        call rk4_solve(simple_ode, t_span, y0, config, t_out, y_out, status)
        
        ! Check final value (should be close to exp(-1) ≈ 0.368)
        pass = (status%code == RK4_SUCCESS) .and. &
               (abs(y_out(1, size(t_out)) - exp(-1.0_wp)) < 0.01_wp)
    end function
    
    pure subroutine simple_ode(t, y, dydt)
        real(wp), intent(in) :: t, y(:)
        real(wp), intent(out) :: dydt(:)
        dydt(1) = -y(1)
    end subroutine
    
    subroutine test_dormand_prince_module()
        use dormand_prince
        print *, ''
        print *, '--- Testing DORMAND_PRINCE module ---'
        
        call run_test('DOPRI adaptive solver', test_dopri_adaptive())
    end subroutine
    
    function test_dopri_adaptive() result(pass)
        use dormand_prince
        logical :: pass
        type(dopri_config_t) :: config
        type(dopri_status_t) :: status
        real(wp) :: y0(1), t_span(2)
        real(wp), allocatable :: t_out(:), y_out(:,:)
        
        y0 = [1.0_wp]
        t_span = [0.0_wp, 1.0_wp]
        config%rtol = 1.0e-6_wp
        config%atol = 1.0e-8_wp
        
        call dopri_solve(simple_ode, t_span, y0, config, t_out, y_out, status)
        
        pass = (status%code == DOPRI_SUCCESS) .and. &
               (abs(y_out(1, size(t_out)) - exp(-1.0_wp)) < 0.001_wp)
    end function
    
    subroutine test_backward_euler_module()
        use backward_euler
        print *, ''
        print *, '--- Testing BACKWARD_EULER module ---'
        
        call run_test('Backward Euler for stiff ODE', test_beuler_stiff())
    end subroutine
    
    function test_beuler_stiff() result(pass)
        use backward_euler
        logical :: pass
        type(beuler_config_t) :: config
        type(beuler_status_t) :: status
        real(wp) :: y0(1), t_span(2)
        real(wp), allocatable :: t_out(:), y_out(:,:)
        
        y0 = [1.0_wp]
        t_span = [0.0_wp, 0.1_wp]
        config%dt = 0.01_wp
        config%newton_tol = 1.0e-6_wp
        
        call beuler_solve(stiff_ode, t_span, y0, config, t_out, y_out, status)
        
        pass = (status%code == BEULER_SUCCESS) .and. all(is_finite(y_out))
    end function
    
    subroutine stiff_ode(t, y, dydt)
        real(wp), intent(in) :: t, y(:)
        real(wp), intent(out) :: dydt(:)
        dydt(1) = -100.0_wp * y(1)  ! Stiff equation
    end subroutine
    
    !===========================================================================
    ! OPTIMIZATION TESTS
    !===========================================================================
    
    subroutine test_gradient_descent_module()
        use gradient_descent
        print *, ''
        print *, '--- Testing GRADIENT_DESCENT module ---'
        
        call run_test('Gradient descent minimization', test_gd_minimize())
    end subroutine
    
    function test_gd_minimize() result(pass)
        use gradient_descent
        logical :: pass
        type(gd_config_t) :: config
        type(gd_result_t) :: result
        real(wp) :: x0(2), xmin(2)
        
        ! Minimize f(x,y) = x² + y²
        x0 = [1.0_wp, 1.0_wp]
        config%learning_rate = 0.1_wp
        config%tolerance = 1.0e-4_wp
        config%max_iterations = 1000
        
        call gradient_descent_minimize(quadratic_func, quadratic_grad, &
                                       x0, xmin, config, result)
        
        pass = result%converged .and. (norm2(xmin) < 0.01_wp)
    end function
    
    pure function quadratic_func(x) result(f)
        real(wp), intent(in) :: x(:)
        real(wp) :: f
        f = sum(x**2)
    end function
    
    pure function quadratic_grad(x) result(g)
        real(wp), intent(in) :: x(:)
        real(wp) :: g(size(x))
        g = 2.0_wp * x
    end function
    
    subroutine test_conjugate_gradient_module()
        print *, ''
        print *, '--- Testing CONJUGATE_GRADIENT module ---'
        call run_test('CG placeholder', test_cg_placeholder())
    end subroutine
    
    function test_cg_placeholder() result(pass)
        logical :: pass
        pass = .true.  ! Placeholder
    end function
    
    subroutine test_quasi_newton_module()
        print *, ''
        print *, '--- Testing QUASI_NEWTON module ---'
        call run_test('Quasi-Newton placeholder', test_qn_placeholder())
    end subroutine
    
    function test_qn_placeholder() result(pass)
        logical :: pass
        pass = .true.  ! Placeholder
    end function
    
    subroutine test_constrained_module()
        print *, ''
        print *, '--- Testing CONSTRAINED module ---'
        call run_test('Constrained opt placeholder', test_const_placeholder())
    end subroutine
    
    function test_const_placeholder() result(pass)
        logical :: pass
        pass = .true.  ! Placeholder
    end function
    
    !===========================================================================
    ! PDE SOLVER TESTS
    !===========================================================================
    
    subroutine test_finite_difference_module()
        use finite_difference
        print *, ''
        print *, '--- Testing FINITE_DIFFERENCE module ---'
        
        call run_test('1D derivative', test_fd_derivative_1d())
        call run_test('1D Laplacian', test_fd_laplacian_1d())
    end subroutine
    
    function test_fd_derivative_1d() result(pass)
        use finite_difference
        logical :: pass
        real(wp) :: u(5), du(5)
        integer :: i
        
        ! u = x², du/dx = 2x
        do i = 1, 5
            u(i) = real((i-1)**2, wp)
        end do
        
        call fd_derivative_1d(u, du, 1.0_wp, order=1, accuracy=2)
        
        ! Check middle point (should be close to 2*2 = 4)
        pass = abs(du(3) - 4.0_wp) < 0.5_wp
    end function
    
    function test_fd_laplacian_1d() result(pass)
        use finite_difference
        logical :: pass
        real(wp) :: u(5), lap(5)
        
        u = [0.0_wp, 1.0_wp, 0.0_wp, 1.0_wp, 0.0_wp]
        
        call fd_laplacian_1d(u, lap, 1.0_wp)
        
        pass = all(is_finite(lap))
    end function
    
    subroutine test_finite_volume_module()
        use finite_volume
        print *, ''
        print *, '--- Testing FINITE_VOLUME module ---'
        
        call run_test('FV advection 1D', test_fv_advection())
    end subroutine
    
    function test_fv_advection() result(pass)
        use finite_volume
        logical :: pass
        real(wp) :: u(10), u_new(10)
        
        u = 1.0_wp
        
        call fv_advection_1d(u, u_new, 1.0_wp, 0.1_wp, 0.01_wp, FV_UPWIND)
        
        pass = all(is_finite(u_new)) .and. (maxval(u_new) <= 2.0_wp)
    end function
    
    subroutine test_spectral_module()
        use spectral
        print *, ''
        print *, '--- Testing SPECTRAL module ---'
        
        call run_test('Spectral grid creation', test_spectral_grid())
    end subroutine
    
    function test_spectral_grid() result(pass)
        use spectral
        logical :: pass
        type(spectral_grid_1d_t) :: grid
        
        call spectral_grid_1d_create(grid, 32, 2.0_wp * PI)
        
        pass = (grid%n == 32) .and. (abs(grid%length - 2.0_wp * PI) < TOL_DEFAULT)
        
        call spectral_grid_1d_destroy(grid)
    end function
    
    !===========================================================================
    ! PHYSICS MODEL TESTS
    !===========================================================================
    
    subroutine test_heat_transfer_module()
        use heat_transfer
        print *, ''
        print *, '--- Testing HEAT_TRANSFER module ---'
        
        call run_test('Heat initialization', test_heat_init())
        call run_test('Heat time step', test_heat_step())
    end subroutine
    
    function test_heat_init() result(pass)
        use heat_transfer
        logical :: pass
        type(heat_state_t) :: state
        type(heat_config_t) :: config
        
        call heat_init(state, 10, 10, 10, 0.1_wp, 0.1_wp, 0.1_wp, config)
        
        pass = allocated(state%T) .and. (state%nx == 10)
        
        call heat_destroy(state)
    end function
    
    function test_heat_step() result(pass)
        use heat_transfer
        logical :: pass
        type(heat_state_t) :: state
        type(heat_config_t) :: config
        real(wp) :: dt
        
        call heat_init(state, 5, 5, 5, 0.1_wp, 0.1_wp, 0.1_wp, config)
        
        dt = heat_get_max_dt(state)
        call heat_step(state, dt * 0.5_wp)
        
        pass = all(is_finite(state%T)) .and. (state%steps == 1)
        
        call heat_destroy(state)
    end function
    
    subroutine test_fluid_dynamics_module()
        use fluid_dynamics
        print *, ''
        print *, '--- Testing FLUID_DYNAMICS module ---'
        
        call run_test('Fluid initialization', test_fluid_init())
    end subroutine
    
    function test_fluid_init() result(pass)
        use fluid_dynamics
        logical :: pass
        type(fluid_state_t) :: state
        type(fluid_config_t) :: config
        
        call fluid_init(state, 5, 5, 5, 0.1_wp, 0.1_wp, 0.1_wp, config)
        
        pass = allocated(state%vx) .and. allocated(state%p)
        
        call fluid_destroy(state)
    end function
    
    subroutine test_nuclear_fission_module()
        use nuclear_fission
        print *, ''
        print *, '--- Testing NUCLEAR_FISSION module ---'
        
        call run_test('Fission initialization', test_fission_init())
        call run_test('Point kinetics step', test_point_kinetics())
    end subroutine
    
    function test_fission_init() result(pass)
        use nuclear_fission
        logical :: pass
        type(fission_state_t) :: state
        type(fission_config_t) :: config
        
        call fission_init(state, 5, 5, 5, 0.1_wp, 0.1_wp, 0.1_wp, config)
        
        pass = allocated(state%flux) .and. (state%nx == 5)
        
        call fission_destroy(state)
    end function
    
    function test_point_kinetics() result(pass)
        use nuclear_fission
        logical :: pass
        type(fission_state_t) :: state
        type(fission_config_t) :: config
        
        config%use_point_kinetics = .true.
        call fission_init(state, 5, 5, 5, 0.1_wp, 0.1_wp, 0.1_wp, config)
        call fission_set_reactivity(state, 0.001_wp)
        
        call fission_step(state, 0.001_wp)
        
        pass = is_finite(state%neutron_population) .and. (state%steps == 1)
        
        call fission_destroy(state)
    end function
    
    subroutine test_pressure_dynamics_module()
        use pressure_dynamics
        print *, ''
        print *, '--- Testing PRESSURE_DYNAMICS module ---'
        
        call run_test('Pressure initialization', test_pressure_init())
        call run_test('EOS calculation', test_eos())
    end subroutine
    
    function test_pressure_init() result(pass)
        use pressure_dynamics
        logical :: pass
        type(pressure_state_t) :: state
        type(pressure_config_t) :: config
        
        call pressure_init(state, 5, 5, 5, 0.1_wp, 0.1_wp, 0.1_wp, config)
        
        pass = allocated(state%p) .and. (state%nx == 5)
        
        call pressure_destroy(state)
    end function
    
    function test_eos() result(pass)
        use pressure_dynamics
        logical :: pass
        type(pressure_state_t) :: state
        type(pressure_config_t) :: config
        real(wp) :: T(5,5,5), rho(5,5,5)
        
        call pressure_init(state, 5, 5, 5, 0.1_wp, 0.1_wp, 0.1_wp, config)
        
        T = 300.0_wp
        rho = 1000.0_wp
        
        call pressure_update_from_temperature(state, T, rho)
        
        pass = all(is_finite(state%p)) .and. all(state%p > 0.0_wp)
        
        call pressure_destroy(state)
    end function
    
    !===========================================================================
    ! HELPER FUNCTIONS
    !===========================================================================
    
    subroutine run_test(test_name, passed)
        character(len=*), intent(in) :: test_name
        logical, intent(in) :: passed
        
        total_tests = total_tests + 1
        
        if (passed) then
            passed_tests = passed_tests + 1
            print '(A,A,A)', '  [PASS] ', test_name
        else
            failed_tests = failed_tests + 1
            print '(A,A,A)', '  [FAIL] ', test_name
        end if
    end subroutine run_test

end program test_suite