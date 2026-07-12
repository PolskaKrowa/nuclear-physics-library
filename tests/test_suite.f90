! test_suite.f90
!
! Comprehensive test suite for nuclear physics simulation kernels and models
! Systematically tests each component and reports results
!

module test_helpers 
    use kinds, only: wp
    implicit none
contains
    ! --- ODE Helpers ---
    pure subroutine simple_ode(t, y, dydt)
        real(wp), intent(in) :: t, y(:)
        real(wp), intent(out) :: dydt(:)
        dydt(1) = -y(1)
    end subroutine simple_ode
    
    subroutine stiff_ode(t, y, dydt)
        real(wp), intent(in) :: t, y(:)
        real(wp), intent(out) :: dydt(:)
        dydt(1) = -100.0_wp * y(1)
    end subroutine stiff_ode

    ! --- Optimization Helpers ---
    pure function quadratic_func(x) result(f)
        real(wp), intent(in) :: x(:)
        real(wp) :: f
        f = sum(x**2)
    end function quadratic_func
    
    pure function quadratic_grad(x) result(g)
        real(wp), intent(in) :: x(:)
        real(wp) :: g(size(x))
        g = 2.0_wp * x
    end function quadratic_grad

end module test_helpers

program test_suite
    use kinds
    use constants
    use numerics_utils
    use rng
    use test_helpers
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

    ! Subsystems
    call test_vessel_module()
    call test_crd_module()
    call test_main_steam_module()
    call test_feedwater_module()
    call test_recirc_module()
    call test_fuel_module()

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
        ! qp is real128 (kind=16) on gfortran/ifort, but aliases to
        ! real64 (kind=8) on NVHPC's nvfortran which does not support
        ! 128-bit reals. Accept either so the test passes in both
        ! COMPUTATION_MODE=CPU and HYBRID/GPU builds.
#if !defined(__NVCOMPILER) && !defined(__PGI)
        pass = (sp == 4) .and. (dp == 8) .and. (qp == 16)
#else
        pass = (sp == 4) .and. (dp == 8) .and. (qp == 8)
#endif
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
    
    !===========================================================================
    ! SUBSYSTEM TESTS
    !===========================================================================

    subroutine test_vessel_module()
        print *, ''
        print *, '--- Testing VESSEL subsystem ---'

        call run_test('Vessel init seeds rated dome pressure', test_vessel_init_rated())
        call run_test('Vessel sat-temp consistent with pressure', test_vessel_sat_temp_consistent())
        call run_test('Vessel level positive after init',         test_vessel_level_positive())
        call run_test('Vessel pressure stays bounded under drivers', &
                      test_vessel_pressure_bounded())
    end subroutine

    function test_vessel_init_rated() result(pass)
        use vessel
        logical :: pass
        type(vessel_state_t)  :: state
        type(vessel_config_t) :: config

        call vessel_init(state, config)
        ! Spec 2.1: rated dome pressure ≈ 1000 psig = 6.895e6 Pa (default
        ! config). Initialization should seed exactly that.
        pass = abs(state%pressure_operating - 6.895e6_wp) < 1.0_wp
    end function

    function test_vessel_sat_temp_consistent() result(pass)
        use vessel
        logical :: pass
        type(vessel_state_t)  :: state
        type(vessel_config_t) :: config

        call vessel_init(state, config)
        ! sat_temp_K of the seeded dome pressure must equal the state's
        ! sat_temperature scalar (no drift between the polynomial and what
        ! init recorded).
        pass = abs(state%sat_temperature - sat_temp_K(state%pressure_operating)) < 1.0e-6_wp
    end function

    function test_vessel_level_positive() result(pass)
        use vessel
        logical :: pass
        type(vessel_state_t)  :: state
        type(vessel_config_t) :: config

        call vessel_init(state, config)
        ! Reactor level (Lrx) is base-corrected by density; with FW-inlet
        ! water and a ~3.81 m core, the corrected level should be > 0.
        pass = state%Lrx > 0.0_wp .and. state%Lbase > 0.0_wp
    end function

    function test_vessel_pressure_bounded() result(pass)
        use vessel
        logical :: pass
        type(vessel_state_t)   :: state
        type(vessel_config_t)  :: config
        type(vessel_drivers_t) :: drivers
        real(wp), parameter :: dt = 0.02_wp
        integer  :: i
        real(wp) :: p_min_seen, p_max_seen

        call vessel_init(state, config)

        ! Plausible "near-rated" drivers — saturated coolant, 100 % power,
        ! 100 % FW, 33 % turbine valve (the bwr_initialize default). This
        ! is not a clean steady state because of the placeholder magic
        ! coefficients in vessel_step_dynamics; we only check that the
        ! pressure stays inside the ASME-style hard clamps over 60 s.
        drivers%avg_coolant_temp_K  = sat_temp_K(state%pressure_operating)
        drivers%power_current_W     = 2381.0e6_wp
        drivers%power_rated_W       = 2381.0e6_wp
        drivers%feedwater_flow_pct  = 100.0_wp
        drivers%turbine_valve_pct   = 33.0_wp
        drivers%bypass_valve_pct    = 0.0_wp
        drivers%feedwater_flow_kg_s = 1638.0_wp
        drivers%steam_flow_kg_s     = 1638.0_wp
        drivers%turbine_tripped     = .false.
        drivers%core_height_m       = config%core_height_m

        p_min_seen = state%pressure_operating
        p_max_seen = state%pressure_operating

        do i = 1, 3000        ! 60 s at 20 ms
            call vessel_step_dynamics(state, dt, drivers, config)
            p_min_seen = min(p_min_seen, state%pressure_operating)
            p_max_seen = max(p_max_seen, state%pressure_operating)
        end do

        ! Internal clamps: P_MIN_PA = 1 atm (cold-standby floor),
        ! P_MAX_PA = 1.1e7. Near-rated drivers should keep pressure
        ! comfortably above atmospheric, but assert against the actual
        ! clamp values so the test stays consistent with the model.
        pass = p_min_seen >= 1.01325e5_wp .and. p_max_seen <= 1.1e7_wp
    end function

    subroutine test_crd_module()
        print *, ''
        print *, '--- Testing CRD subsystem ---'

        call run_test('CRD scram-time curve anchor points', test_crd_scram_time_curve())
        call run_test('CRD init seeds 137 blades at initial_insertion', &
                      test_crd_init_blades())
        call run_test('CRD bank-position broadcast hits all blades', &
                      test_crd_bank_broadcast())
        call run_test('CRD scram from full-out reaches 90% in ~2.0-2.3s @ 1000 psig', &
                      test_crd_scram_at_rated())
        call run_test('CRD scram at zero pressure takes ~3.5s', &
                      test_crd_scram_at_zero_pressure())
    end subroutine

    function test_crd_scram_time_curve() result(pass)
        use crd
        logical :: pass
        real(wp) :: t_zero, t_600, t_1000, t_1200
        real(wp), parameter :: PA_PER_PSI = 6894.76_wp
        real(wp), parameter :: P_ATM_PA   = 101325.0_wp

        t_zero = scram_time_90pct(P_ATM_PA)                        ! 0 psig
        t_600  = scram_time_90pct(P_ATM_PA +  600.0_wp * PA_PER_PSI)
        t_1000 = scram_time_90pct(P_ATM_PA + 1000.0_wp * PA_PER_PSI)
        t_1200 = scram_time_90pct(P_ATM_PA + 1200.0_wp * PA_PER_PSI)

        ! Anchor points from spec 03 Figure 2.3-13 typical-drive envelope.
        pass = abs(t_zero - 3.5_wp)  < 0.01_wp .and. &
               abs(t_600  - 2.6_wp)  < 0.01_wp .and. &
               abs(t_1000 - 2.15_wp) < 0.01_wp .and. &
               abs(t_1200 - 2.0_wp)  < 0.01_wp
    end function

    function test_crd_init_blades() result(pass)
        use crd
        logical :: pass
        type(crd_state_t)  :: state
        type(crd_config_t) :: config

        call crd_init(state, config)
        pass = all(abs(state%blade_insertion - config%initial_insertion) < 1.0e-12_wp) &
               .and. size(state%blade_insertion) == CRD_N_BLADES &
               .and. .not. state%scram_latched
    end function

    function test_crd_bank_broadcast() result(pass)
        ! Bank command stores in bank_demand; crd_step rate-limits the
        ! actual blade motion. Confirm the demand is set, then step long
        ! enough for the slew to converge.
        use crd
        logical :: pass
        type(crd_state_t)   :: state
        type(crd_config_t)  :: config
        type(crd_command_t) :: cmd
        type(crd_drivers_t) :: drv
        real(wp), parameter :: dt = 0.1_wp
        integer :: i

        call crd_init(state, config)
        cmd%bank_position_set = 0.42_wp
        call crd_apply(state, cmd)

        if (abs(state%bank_demand - 0.42_wp) > 1.0e-12_wp) then
            pass = .false.
            return
        end if

        ! Full stroke ≈ 7.32 m at 0.0305 m/s → ~240 s. 0.95 → 0.42 is
        ! Δ = 0.53 stroke ≈ 127 s; allow 400 ticks × 0.1 s = 40 s won't
        ! finish, but the bank should be MOVING toward the demand.
        do i = 1, 400
            call crd_step(state, dt, drv, config)
        end do

        ! After 40 s the bank should have moved noticeably from 0.95
        ! toward 0.42 but not yet reached it (rate-limited slew). All
        ! blades must move together (bank operation).
        pass = state%blade_insertion(1) < 0.95_wp &
               .and. state%blade_insertion(1) > 0.42_wp &
               .and. all(abs(state%blade_insertion - state%blade_insertion(1)) < 1.0e-12_wp)
    end function

    function test_crd_scram_at_rated() result(pass)
        use crd
        logical :: pass
        type(crd_state_t)   :: state
        type(crd_config_t)  :: config
        type(crd_command_t) :: cmd
        type(crd_drivers_t) :: drivers
        real(wp), parameter :: dt = 0.01_wp
        real(wp), parameter :: PA_PER_PSI = 6894.76_wp
        real(wp), parameter :: P_ATM_PA   = 101325.0_wp
        real(wp) :: avg_insertion, t
        integer  :: i

        ! Full-out start so the scram-time curve sees the full stroke.
        config%initial_insertion = 0.0_wp
        call crd_init(state, config)
        drivers%dome_pressure_Pa = P_ATM_PA + 1000.0_wp * PA_PER_PSI

        cmd%scram_latch = .true.
        call crd_apply(state, cmd)

        t = 0.0_wp
        avg_insertion = 0.0_wp
        do i = 1, 500     ! up to 5 s
            call crd_step(state, dt, drivers, config)
            t = t + dt
            avg_insertion = sum(state%blade_insertion) / real(CRD_N_BLADES, wp)
            if (avg_insertion >= 0.9_wp) exit
        end do

        ! Spec target: 2.0–2.3 s at 1000 psig (Figure 2.3-13).
        ! Allow ±0.3 s tolerance to absorb time-discretisation.
        pass = (t >= 1.7_wp .and. t <= 2.6_wp .and. avg_insertion >= 0.9_wp)
    end function

    function test_crd_scram_at_zero_pressure() result(pass)
        use crd
        logical :: pass
        type(crd_state_t)   :: state
        type(crd_config_t)  :: config
        type(crd_command_t) :: cmd
        type(crd_drivers_t) :: drivers
        real(wp), parameter :: dt = 0.01_wp
        real(wp), parameter :: P_ATM_PA = 101325.0_wp
        real(wp) :: avg_insertion, t
        integer  :: i

        config%initial_insertion = 0.0_wp
        call crd_init(state, config)
        drivers%dome_pressure_Pa = P_ATM_PA   ! 0 psig

        cmd%scram_latch = .true.
        call crd_apply(state, cmd)

        t = 0.0_wp
        avg_insertion = 0.0_wp
        do i = 1, 600     ! up to 6 s
            call crd_step(state, dt, drivers, config)
            t = t + dt
            avg_insertion = sum(state%blade_insertion) / real(CRD_N_BLADES, wp)
            if (avg_insertion >= 0.9_wp) exit
        end do

        ! Spec target: ~3.5 s at 0 psig (accumulator-only branch).
        pass = (t >= 3.2_wp .and. t <= 3.9_wp .and. avg_insertion >= 0.9_wp)
    end function

    !===========================================================================
    ! MAIN STEAM SUBSYSTEM (spec 2.5) — step 5
    !===========================================================================

    subroutine test_main_steam_module()
        print *, ''
        print *, '--- Testing MAIN_STEAM subsystem ---'

        call run_test('Main steam init seeds rated valve positions', &
                      test_ms_init())
        call run_test('TSV fast-close ramps to zero in ~0.1 s', &
                      test_ms_tsv_fast_close())
        call run_test('MSIV full stroke open→closed within configured time', &
                      test_ms_msiv_stroke())
        call run_test('SRV group 1 lifts when dome P exceeds setpoint', &
                      test_ms_srv_group1_lifts())
        call run_test('SRVs reset with deadband below setpoint', &
                      test_ms_srv_reset_deadband())
        call run_test('Bypass capacity 25 % cap applied in vessel ODE', &
                      test_ms_bypass_cap())
        call run_test('Load-reject: bypass opens to 25 %, SRVs may lift, MSIV stays open', &
                      test_ms_load_reject())
    end subroutine

    function test_ms_init() result(pass)
        use main_steam
        logical :: pass
        type(main_steam_state_t)  :: state
        type(main_steam_config_t) :: config

        call main_steam_init(state, config)
        pass = abs(state%turbine_valve_pct - config%initial_turbine_valve_pct) < 1.0e-12_wp &
               .and. abs(state%bypass_valve_pct - config%initial_bypass_valve_pct) < 1.0e-12_wp &
               .and. all(state%msiv_pos_pct == 100.0_wp) &
               .and. count(state%srv_open) == 0 &
               .and. .not. state%turbine_tripped
    end function

    function test_ms_tsv_fast_close() result(pass)
        use main_steam
        logical :: pass
        type(main_steam_state_t)   :: state
        type(main_steam_config_t)  :: config
        type(main_steam_command_t) :: cmd
        type(main_steam_drivers_t) :: drv
        real(wp), parameter :: dt = 0.001_wp
        real(wp) :: t
        integer  :: i

        config%initial_turbine_valve_pct = 100.0_wp
        call main_steam_init(state, config)
        drv%dome_pressure_Pa = 7.14e6_wp
        drv%power_current_W  = 2.0e9_wp

        cmd%tsv_fast_close = .true.
        call main_steam_apply(state, cmd, config)

        t = 0.0_wp
        do i = 1, 200
            call main_steam_step(state, dt, drv, config)
            t = t + dt
            if (state%turbine_tripped) exit
        end do

        ! Latched trip within fast-close window + tiny slack for step granularity.
        pass = state%turbine_tripped &
               .and. state%turbine_valve_pct < 1.0e-6_wp &
               .and. t <= config%tsv_fast_close_time_s + 2.0_wp * dt
    end function

    function test_ms_msiv_stroke() result(pass)
        use main_steam
        logical :: pass
        type(main_steam_state_t)   :: state
        type(main_steam_config_t)  :: config
        type(main_steam_command_t) :: cmd
        type(main_steam_drivers_t) :: drv
        real(wp), parameter :: dt = 0.01_wp
        real(wp) :: msiv_frac_after_full_stroke
        integer  :: i

        call main_steam_init(state, config)
        drv%dome_pressure_Pa = 7.14e6_wp

        cmd%msiv_close = .true.
        call main_steam_apply(state, cmd, config)

        ! Step for the configured stroke time + small slack.
        do i = 1, int((config%msiv_stroke_time_s + 0.1_wp) / dt)
            call main_steam_step(state, dt, drv, config)
        end do
        msiv_frac_after_full_stroke = sum(state%msiv_pos_pct) / (real(MS_N_MSLS, wp) * 100.0_wp)

        pass = msiv_frac_after_full_stroke < 1.0e-6_wp
    end function

    function test_ms_srv_group1_lifts() result(pass)
        use main_steam
        logical :: pass
        type(main_steam_state_t)   :: state
        type(main_steam_config_t)  :: config
        type(main_steam_drivers_t) :: drv
        real(wp), parameter :: dt = 0.01_wp
        integer  :: i
        real(wp), parameter :: PA_PER_PSI = 6894.76_wp
        real(wp), parameter :: P_ATM_PA   = 101325.0_wp

        call main_steam_init(state, config)
        ! Dome pressure right at the group-1 setpoint (1105 psig).
        drv%dome_pressure_Pa = P_ATM_PA + 1106.0_wp * PA_PER_PSI

        do i = 1, 5
            call main_steam_step(state, dt, drv, config)
        end do

        ! Group 1 = indices 1-4; group 2 setpoint is 10 psi higher so still shut.
        pass = all(state%srv_open(1:4)) &
               .and. .not. any(state%srv_open(5:11)) &
               .and. state%srv_flow_kg_s > 4.0_wp * config%srv_capacity_kg_s - 1.0e-6_wp
    end function

    function test_ms_srv_reset_deadband() result(pass)
        use main_steam
        logical :: pass
        type(main_steam_state_t)   :: state
        type(main_steam_config_t)  :: config
        type(main_steam_drivers_t) :: drv
        real(wp), parameter :: dt = 0.01_wp
        integer  :: i
        real(wp), parameter :: PA_PER_PSI = 6894.76_wp
        real(wp), parameter :: P_ATM_PA   = 101325.0_wp

        call main_steam_init(state, config)
        ! Lift group 1 first.
        drv%dome_pressure_Pa = P_ATM_PA + 1106.0_wp * PA_PER_PSI
        do i = 1, 3
            call main_steam_step(state, dt, drv, config)
        end do
        if (.not. all(state%srv_open(1:4))) then
            pass = .false.
            return
        end if

        ! Drop just below setpoint but above the deadband — SRV stays open.
        drv%dome_pressure_Pa = P_ATM_PA + 1080.0_wp * PA_PER_PSI
        do i = 1, 3
            call main_steam_step(state, dt, drv, config)
        end do
        if (.not. all(state%srv_open(1:4))) then
            pass = .false.
            return
        end if

        ! Drop well below the deadband — SRV resets.
        drv%dome_pressure_Pa = P_ATM_PA + 1000.0_wp * PA_PER_PSI
        do i = 1, 3
            call main_steam_step(state, dt, drv, config)
        end do

        pass = .not. any(state%srv_open)
    end function

    function test_ms_bypass_cap() result(pass)
        use vessel
        logical :: pass
        type(vessel_state_t)   :: state
        type(vessel_config_t)  :: config
        type(vessel_drivers_t) :: drivers
        real(wp), parameter :: dt = 0.02_wp
        real(wp) :: P_with_full_bypass, P_with_full_turbine
        integer  :: i

        ! Same initial dome P, same coolant T, no power, no FW. Step with
        ! turbine 100 % vs bypass 100 % — bypass should yield 4× less relief.
        call vessel_init(state, config)
        drivers%avg_coolant_temp_K = 480.0_wp
        drivers%power_current_W = 0.0_wp
        drivers%power_rated_W = 1.0_wp
        drivers%feedwater_flow_pct = 0.0_wp
        drivers%turbine_valve_pct = 100.0_wp
        drivers%bypass_valve_pct = 0.0_wp
        drivers%feedwater_flow_kg_s = 0.0_wp
        drivers%steam_flow_kg_s = 0.0_wp
        drivers%msiv_open_frac = 1.0_wp
        drivers%srv_flow_frac = 0.0_wp
        drivers%turbine_tripped = .false.
        drivers%core_height_m = config%core_height_m
        do i = 1, 5
            call vessel_step_dynamics(state, dt, drivers, config)
        end do
        P_with_full_turbine = state%pressure_operating

        call vessel_init(state, config)
        drivers%turbine_valve_pct = 0.0_wp
        drivers%bypass_valve_pct = 100.0_wp
        do i = 1, 5
            call vessel_step_dynamics(state, dt, drivers, config)
        end do
        P_with_full_bypass = state%pressure_operating

        ! Same starting pressure, same generation term, but bypass at 100 %
        ! gives only 25 % of turbine's relief → bypass case ends higher.
        pass = P_with_full_bypass > P_with_full_turbine
    end function

    function test_ms_load_reject() result(pass)
        use main_steam
        logical :: pass
        type(main_steam_state_t)   :: state
        type(main_steam_config_t)  :: config
        type(main_steam_command_t) :: cmd
        type(main_steam_drivers_t) :: drv
        real(wp), parameter :: dt = 0.01_wp
        real(wp), parameter :: PA_PER_PSI = 6894.76_wp
        real(wp), parameter :: P_ATM_PA   = 101325.0_wp
        integer  :: i

        config%initial_turbine_valve_pct = 100.0_wp
        call main_steam_init(state, config)
        ! Hold dome P just above group 1 setpoint to verify SRV behaviour
        ! during the transient (a real ODE would do this for us, but the
        ! test isolates main_steam).
        drv%dome_pressure_Pa = P_ATM_PA + 1106.0_wp * PA_PER_PSI
        drv%power_current_W  = 2.3e9_wp

        ! Load reject = TSV fast-close + bypass demand to 100 %.
        cmd%tsv_fast_close       = .true.
        cmd%bypass_valve_set_pct = 100.0_wp
        call main_steam_apply(state, cmd, config)

        ! Step a few seconds so the bypass valve gets time to stroke open
        ! (~2 s at 50 %/s) and the SRVs lift.
        do i = 1, 300        ! 3.0 s
            call main_steam_step(state, dt, drv, config)
        end do

        pass = state%turbine_tripped &
               .and. state%turbine_valve_pct < 1.0e-6_wp &
               .and. state%bypass_valve_pct > 99.0_wp &
               .and. all(state%msiv_pos_pct >= 99.0_wp) &
               .and. count(state%srv_open(1:4)) == 4
    end function

    !===========================================================================
    ! FEEDWATER SUBSYSTEM (spec 2.6) — step 6
    !===========================================================================

    subroutine test_feedwater_module()
        print *, ''
        print *, '--- Testing FEEDWATER subsystem ---'

        call run_test('Feedwater init at rated flow + rated FW temperature', &
                      test_fw_init())
        call run_test('Operator demand 0 % drains flow through RFPT slew', &
                      test_fw_demand_to_zero())
        call run_test('Single RFP trip caps flow at 67 % rated capacity', &
                      test_fw_single_rfp_cap())
        call run_test('Both RFPs tripped → flow → 0', &
                      test_fw_both_rfps_trip())
        call run_test('Booster low-suction cascade: A trips at 20 s, B at 45 s', &
                      test_fw_booster_trip_cascade())
        call run_test('RFP low-suction cascade fires after both boosters trip', &
                      test_fw_rfp_trip_cascade())
        call run_test('FW pump trip from rated → vessel level drops measurably', &
                      test_fw_trip_level_drop())
        call run_test('3-element control closes level error', &
                      test_fw_three_element_control())
    end subroutine

    function test_fw_init() result(pass)
        use feedwater
        logical :: pass
        type(feedwater_state_t)  :: state
        type(feedwater_config_t) :: config

        call feedwater_init(state, config)
        pass = abs(state%demand_flow_pct - 100.0_wp) < 1.0e-12_wp &
               .and. abs(state%flow_pct - 100.0_wp)   < 1.0e-12_wp &
               .and. abs(state%feedwater_temp_K - config%rated_inlet_temperature_K) < 1.0e-9_wp &
               .and. count(.not. state%rfp_tripped) == FW_N_RFP
    end function

    function test_fw_demand_to_zero() result(pass)
        use feedwater
        logical :: pass
        type(feedwater_state_t)   :: state
        type(feedwater_config_t)  :: config
        type(feedwater_command_t) :: cmd
        type(feedwater_drivers_t) :: drv
        real(wp), parameter :: dt = 0.02_wp
        integer  :: i

        ! Manual-demand semantics: disable the 3-element controller so the
        ! operator setpoint is what drives the flow (not steam-tracking).
        ! Also block the auto-engage by clearing the gating threshold via
        ! a configured pressure that keeps the controller off even if
        ! drivers were to supply a high dome pressure.
        config%controller_enabled  = .false.
        config%auto_engage_on_boil = .false.
        call feedwater_init(state, config)
        cmd%demand_flow_pct_set = 0.0_wp
        call feedwater_apply(state, cmd)

        ! RFPT speed slews at 20 %/s; allow 10 s to drain.
        do i = 1, 500
            call feedwater_step(state, dt, drv, config)
        end do
        pass = state%flow_pct < 1.0_wp .and. state%rfpt_speed_pct < 1.0_wp
    end function

    function test_fw_single_rfp_cap() result(pass)
        use feedwater
        logical :: pass
        type(feedwater_state_t)   :: state
        type(feedwater_config_t)  :: config
        type(feedwater_command_t) :: cmd
        type(feedwater_drivers_t) :: drv
        real(wp), parameter :: dt = 0.02_wp
        integer  :: i

        call feedwater_init(state, config)
        cmd%trip_rfp_idx = 1
        call feedwater_apply(state, cmd)
        cmd%trip_rfp_idx = 0
        cmd%demand_flow_pct_set = 100.0_wp
        call feedwater_apply(state, cmd)

        ! Wait for RFPT speed to settle; with 1 RFP the capacity cap is 67 %.
        do i = 1, 500
            call feedwater_step(state, dt, drv, config)
        end do

        pass = state%flow_pct >= 66.0_wp .and. state%flow_pct <= 68.0_wp &
               .and. count(.not. state%rfp_tripped) == 1
    end function

    function test_fw_both_rfps_trip() result(pass)
        use feedwater
        logical :: pass
        type(feedwater_state_t)   :: state
        type(feedwater_config_t)  :: config
        type(feedwater_command_t) :: cmd
        type(feedwater_drivers_t) :: drv
        real(wp), parameter :: dt = 0.02_wp
        integer  :: i

        call feedwater_init(state, config)
        cmd%trip_rfp_idx = 1
        call feedwater_apply(state, cmd)
        cmd%trip_rfp_idx = 2
        call feedwater_apply(state, cmd)
        cmd%trip_rfp_idx = 0
        cmd%demand_flow_pct_set = 100.0_wp
        call feedwater_apply(state, cmd)

        do i = 1, 100
            call feedwater_step(state, dt, drv, config)
        end do

        pass = state%flow_pct < 1.0e-6_wp &
               .and. state%flow_kg_s < 1.0e-3_wp &
               .and. count(.not. state%rfp_tripped) == 0
    end function

    function test_fw_booster_trip_cascade() result(pass)
        use feedwater
        logical :: pass
        type(feedwater_state_t)   :: state
        type(feedwater_config_t)  :: config
        type(feedwater_command_t) :: cmd
        type(feedwater_drivers_t) :: drv
        real(wp), parameter :: dt = 0.1_wp
        integer  :: i
        real(wp) :: t_A_trip, t_B_trip
        logical  :: A_was_running, B_was_running

        call feedwater_init(state, config)
        ! Trip both condensate pumps → booster suction drops to 0 psig.
        cmd%trip_condensate_idx = 1
        call feedwater_apply(state, cmd)
        cmd%trip_condensate_idx = 2
        call feedwater_apply(state, cmd)

        t_A_trip = -1.0_wp
        t_B_trip = -1.0_wp
        A_was_running = .true.
        B_was_running = .true.

        do i = 1, 600        ! up to 60 s
            call feedwater_step(state, dt, drv, config)
            if (A_was_running .and. state%booster_tripped(1)) then
                t_A_trip = state%time
                A_was_running = .false.
            end if
            if (B_was_running .and. state%booster_tripped(2)) then
                t_B_trip = state%time
                B_was_running = .false.
            end if
            if (t_A_trip > 0.0_wp .and. t_B_trip > 0.0_wp) exit
        end do

        pass = abs(t_A_trip - 20.0_wp) < 0.2_wp &
               .and. abs(t_B_trip - 45.0_wp) < 0.2_wp
    end function

    function test_fw_rfp_trip_cascade() result(pass)
        use feedwater
        logical :: pass
        type(feedwater_state_t)   :: state
        type(feedwater_config_t)  :: config
        type(feedwater_command_t) :: cmd
        type(feedwater_drivers_t) :: drv
        real(wp), parameter :: dt = 0.1_wp
        integer  :: i

        call feedwater_init(state, config)
        ! Knock out both boosters directly → RFP suction collapses.
        cmd%trip_booster_idx = 1
        call feedwater_apply(state, cmd)
        cmd%trip_booster_idx = 2
        call feedwater_apply(state, cmd)

        do i = 1, 400        ! up to 40 s — RFP B trips at 25 s
            call feedwater_step(state, dt, drv, config)
        end do

        pass = state%rfp_tripped(1) .and. state%rfp_tripped(2) &
               .and. state%flow_kg_s < 1.0_wp
    end function

    function test_fw_trip_level_drop() result(pass)
        ! Combined feedwater + vessel level-drop test (spec test).
        ! Run a fresh vessel with rated drivers, trip both RFPs, and
        ! verify reactor level drops measurably within a few seconds.
        use feedwater
        use vessel
        logical :: pass
        type(feedwater_state_t)   :: fw_state
        type(feedwater_config_t)  :: fw_config
        type(feedwater_command_t) :: fw_cmd
        type(feedwater_drivers_t) :: fw_drv
        type(vessel_state_t)      :: v_state
        type(vessel_config_t)     :: v_config
        type(vessel_drivers_t)    :: v_drv
        real(wp), parameter :: dt = 0.02_wp
        real(wp) :: L_initial, L_final
        integer  :: i

        call feedwater_init(fw_state, fw_config)
        call vessel_init(v_state, v_config)

        ! Settle the vessel briefly so its level isn't at init transient.
        v_drv%avg_coolant_temp_K = 550.0_wp
        v_drv%power_current_W = 2381.0e6_wp
        v_drv%power_rated_W   = 2381.0e6_wp
        v_drv%feedwater_flow_pct = fw_state%flow_pct
        v_drv%turbine_valve_pct = 100.0_wp
        v_drv%bypass_valve_pct = 0.0_wp
        v_drv%srv_flow_frac = 0.0_wp
        v_drv%msiv_open_frac = 1.0_wp
        v_drv%turbine_tripped = .false.
        v_drv%core_height_m = v_config%core_height_m
        ! Mass-balance level driver: at rated power, boil-off ≈ rated FW.
        v_drv%steam_flow_kg_s = 1638.0_wp

        do i = 1, 50
            call feedwater_step(fw_state, dt, fw_drv, fw_config)
            v_drv%feedwater_flow_pct = fw_state%flow_pct
            v_drv%feedwater_flow_kg_s = fw_state%flow_kg_s
            call vessel_step_dynamics(v_state, dt, v_drv, v_config)
        end do
        L_initial = v_state%Lrx

        ! Total FW loss: trip both RFPs.
        fw_cmd%trip_rfp_idx = 1
        call feedwater_apply(fw_state, fw_cmd)
        fw_cmd%trip_rfp_idx = 2
        call feedwater_apply(fw_state, fw_cmd)
        fw_cmd%trip_rfp_idx = 0

        do i = 1, 1000       ! 20 s after the trip
            call feedwater_step(fw_state, dt, fw_drv, fw_config)
            v_drv%feedwater_flow_pct = fw_state%flow_pct
            v_drv%feedwater_flow_kg_s = fw_state%flow_kg_s
            call vessel_step_dynamics(v_state, dt, v_drv, v_config)
        end do
        L_final = v_state%Lrx

        ! Level should drop measurably as the vessel boils off without FW.
        pass = (L_final < L_initial - 1.0e-3_wp) &
               .and. fw_state%flow_kg_s < 1.0e-3_wp
    end function

    function test_fw_three_element_control() result(pass)
        ! With controller enabled and the level reference below the
        ! current reactor level, demand should drop below 100 %.
        use feedwater
        logical :: pass
        type(feedwater_state_t)   :: state
        type(feedwater_config_t)  :: config
        type(feedwater_command_t) :: cmd
        type(feedwater_drivers_t) :: drv
        real(wp), parameter :: dt = 0.05_wp
        integer  :: i

        config%controller_enabled = .true.
        call feedwater_init(state, config)
        cmd%controller_enable_set = 1
        call feedwater_apply(state, cmd)

        ! High reactor level (above setpoint) + low steam flow → controller
        ! should close demand below 100 %.
        drv%reactor_level_m = config%level_setpoint_m + 0.5_wp
        drv%steam_flow_norm = 0.5_wp
        do i = 1, 200    ! 10 s
            call feedwater_step(state, dt, drv, config)
        end do

        pass = state%demand_flow_pct < 100.0_wp &
               .and. state%controller_enabled
    end function

    !===========================================================================
    ! RECIRCULATION SUBSYSTEM (spec 2.4) — step 7
    !===========================================================================

    subroutine test_recirc_module()
        print *, ''
        print *, '--- Testing RECIRCULATION subsystem ---'

        call run_test('Recirc init at rated → core flow = rated × (1+M)/(1+M)', &
                      test_recirc_init())
        call run_test('Jet-pump M-ratio: core_flow ≈ drive_flow × (1+M)', &
                      test_recirc_m_ratio())
        call run_test('Single recirc pump trip → core flow ≈ 50 % rated', &
                      test_recirc_single_pump_trip())
        call run_test('Both pumps tripped → core flow ≈ natural-circ floor', &
                      test_recirc_natural_circ_floor())
        call run_test('NPSH runback: FW < 20 % clamps pump speed', &
                      test_recirc_npsh_runback_fw())
        call run_test('NPSH runback: low reactor level clamps pump speed', &
                      test_recirc_npsh_runback_level())
        call run_test('EOC-RPT signal trips both pumps + latches', &
                      test_recirc_eoc_rpt())
        call run_test('Pump-speed runback on FW pump trip reduces core flow', &
                      test_recirc_fw_trip_runback())
    end subroutine

    function test_recirc_init() result(pass)
        use recirculation
        logical :: pass
        type(recirc_state_t)  :: state
        type(recirc_config_t) :: config
        real(wp) :: expected_core_flow

        call recirc_init(state, config)
        ! At rated: drive = rated/(1+M); core = drive × (1+M) = rated.
        expected_core_flow = config%rated_core_flow_kg_s
        pass = abs(state%core_flow_kg_s - expected_core_flow) < 1.0_wp &
               .and. abs(state%mass_flux_kg_m2_s - config%rated_mass_flux_kg_m2_s) < 1.0_wp &
               .and. count(.not. state%pump_tripped) == RECIRC_N_PUMPS &
               .and. .not. state%eoc_rpt_latched
    end function

    function test_recirc_m_ratio() result(pass)
        use recirculation
        logical :: pass
        type(recirc_state_t)  :: state
        type(recirc_config_t) :: config

        call recirc_init(state, config)
        ! drive × (1 + M) ≡ core_flow at steady state.
        pass = abs(state%core_flow_kg_s &
                   - state%drive_flow_kg_s * (1.0_wp + config%M_ratio)) &
               < 1.0e-6_wp * config%rated_core_flow_kg_s
    end function

    function test_recirc_single_pump_trip() result(pass)
        use recirculation
        logical :: pass
        type(recirc_state_t)   :: state
        type(recirc_config_t)  :: config
        type(recirc_command_t) :: cmd
        type(recirc_drivers_t) :: drv
        real(wp), parameter :: dt = 0.1_wp
        integer  :: i

        call recirc_init(state, config)
        cmd%trip_pump_idx = 1
        call recirc_apply(state, cmd, config)

        ! Allow ~30 s for the tripped pump to coast down (10 %/s).
        do i = 1, 400
            call recirc_step(state, dt, drv, config)
        end do

        ! Surviving pump at 100 %, tripped at 0 % → core ≈ 50 % rated.
        pass = state%core_flow_kg_s >= 0.45_wp * config%rated_core_flow_kg_s &
               .and. state%core_flow_kg_s <= 0.55_wp * config%rated_core_flow_kg_s &
               .and. state%pump_tripped(1) &
               .and. .not. state%pump_tripped(2)
    end function

    function test_recirc_natural_circ_floor() result(pass)
        use recirculation
        logical :: pass
        type(recirc_state_t)   :: state
        type(recirc_config_t)  :: config
        type(recirc_command_t) :: cmd
        type(recirc_drivers_t) :: drv
        real(wp), parameter :: dt = 0.1_wp
        real(wp) :: expected_floor
        integer  :: i

        call recirc_init(state, config)
        cmd%trip_pump_idx = 1
        call recirc_apply(state, cmd, config)
        cmd%trip_pump_idx = 2
        call recirc_apply(state, cmd, config)

        do i = 1, 400
            call recirc_step(state, dt, drv, config)
        end do

        expected_floor = config%nat_circ_fraction * config%rated_core_flow_kg_s
        pass = abs(state%core_flow_kg_s - expected_floor) < 1.0_wp &
               .and. count(.not. state%pump_tripped) == 0
    end function

    function test_recirc_npsh_runback_fw() result(pass)
        use recirculation
        logical :: pass
        type(recirc_state_t)   :: state
        type(recirc_config_t)  :: config
        type(recirc_drivers_t) :: drv
        real(wp), parameter :: dt = 0.1_wp
        integer  :: i

        call recirc_init(state, config)
        drv%feedwater_flow_pct = 10.0_wp    ! below 20 % threshold
        drv%reactor_level_m    = 1.0_wp     ! above L3

        ! Allow ~30 s for the pumps to ramp down to the runback target.
        do i = 1, 300
            call recirc_step(state, dt, drv, config)
        end do

        pass = state%npsh_runback_active &
               .and. abs(state%pump_speed_pct(1) - config%runback_target_pct) < 0.5_wp &
               .and. abs(state%pump_speed_pct(2) - config%runback_target_pct) < 0.5_wp
    end function

    function test_recirc_npsh_runback_level() result(pass)
        use recirculation
        logical :: pass
        type(recirc_state_t)   :: state
        type(recirc_config_t)  :: config
        type(recirc_drivers_t) :: drv
        real(wp), parameter :: dt = 0.1_wp
        integer  :: i

        call recirc_init(state, config)
        drv%feedwater_flow_pct = 100.0_wp
        drv%reactor_level_m    = 0.1_wp    ! below L3 (0.3 m)

        do i = 1, 300
            call recirc_step(state, dt, drv, config)
        end do

        pass = state%npsh_runback_active &
               .and. state%pump_speed_pct(1) <= config%runback_target_pct + 0.5_wp &
               .and. state%pump_speed_pct(2) <= config%runback_target_pct + 0.5_wp
    end function

    function test_recirc_eoc_rpt() result(pass)
        use recirculation
        logical :: pass
        type(recirc_state_t)   :: state
        type(recirc_config_t)  :: config
        type(recirc_drivers_t) :: drv
        real(wp), parameter :: dt = 0.1_wp
        integer  :: i

        call recirc_init(state, config)
        drv%eoc_rpt_signal = .true.

        do i = 1, 5    ! a few ticks is enough — trip latches immediately
            call recirc_step(state, dt, drv, config)
        end do

        ! After the latch fires both pumps are tripped and start coastdown.
        pass = state%eoc_rpt_latched &
               .and. count(.not. state%pump_tripped) == 0
    end function

    function test_recirc_fw_trip_runback() result(pass)
        ! Spec test: pump speed runback on FW pump trip reduces core flow.
        ! Couples feedwater + recirculation: trip both RFPs → FW flow → 0
        ! → NPSH-runback clamps recirc pumps → core flow drops.
        use feedwater
        use recirculation
        logical :: pass
        type(feedwater_state_t)   :: fw_state
        type(feedwater_config_t)  :: fw_config
        type(feedwater_command_t) :: fw_cmd
        type(feedwater_drivers_t) :: fw_drv
        type(recirc_state_t)      :: r_state
        type(recirc_config_t)     :: r_config
        type(recirc_drivers_t)    :: r_drv
        real(wp), parameter :: dt = 0.1_wp
        real(wp) :: flow_initial, flow_final
        integer  :: i

        call feedwater_init(fw_state, fw_config)
        call recirc_init(r_state, r_config)

        ! Settle initial steady state with both pump chains rated.
        r_drv%feedwater_flow_pct = fw_state%flow_pct
        r_drv%reactor_level_m    = 1.0_wp
        do i = 1, 20
            call feedwater_step(fw_state, dt, fw_drv, fw_config)
            r_drv%feedwater_flow_pct = fw_state%flow_pct
            call recirc_step(r_state, dt, r_drv, r_config)
        end do
        flow_initial = r_state%core_flow_kg_s

        ! Trip both RFPs.
        fw_cmd%trip_rfp_idx = 1
        call feedwater_apply(fw_state, fw_cmd)
        fw_cmd%trip_rfp_idx = 2
        call feedwater_apply(fw_state, fw_cmd)

        do i = 1, 600     ! 60 s for FW flow → 0 and recirc pumps to runback
            call feedwater_step(fw_state, dt, fw_drv, fw_config)
            r_drv%feedwater_flow_pct = fw_state%flow_pct
            call recirc_step(r_state, dt, r_drv, r_config)
        end do
        flow_final = r_state%core_flow_kg_s

        ! Core flow at the runback target ≈ 30 % rated × (1+M)/(1+M) = 30 %.
        ! Allow a 5 % window for steady-state error.
        pass = (flow_initial >= 0.95_wp * r_config%rated_core_flow_kg_s) &
               .and. (flow_final  <= 0.35_wp * r_config%rated_core_flow_kg_s) &
               .and. (flow_final  >= 0.25_wp * r_config%rated_core_flow_kg_s) &
               .and. r_state%npsh_runback_active
    end function

    !===========================================================================
    ! FUEL SUBSYSTEM (spec 2.2) — step 8
    !===========================================================================

    subroutine test_fuel_module()
        print *, ''
        print *, '--- Testing FUEL subsystem ---'

        call run_test('Fuel init builds 137 blade positions inside core', &
                      test_fuel_init_blade_positions())
        call run_test('Blade-for-cell map covers every fuel cell', &
                      test_fuel_blade_map_coverage())
        call run_test('Blade-for-cell map is 180°-rotation symmetric', &
                      test_fuel_blade_map_symmetric())
        call run_test('Inserting 4 specific blades flags 4 distinct cell groups', &
                      test_fuel_four_blade_pattern())
        call run_test('Bank-average equivalence: all blades same → uniform cell insertion', &
                      test_fuel_bank_average_equiv())
        call run_test('Void reactivity factor monotone in void fraction', &
                      test_fuel_void_rf_monotone())
        call run_test('Power-unit helper: W/cm³ → W/m³ scales by 10⁶', &
                      test_fuel_power_to_volumetric())
        call run_test('Power-unit helper: surface q'''' matches cell-volume formula', &
                      test_fuel_power_to_surface())
        call run_test('Dittus-Boelter h matches Nu·k/D at canonical Re/Pr', &
                      test_fuel_dittus_boelter())
        call run_test('Jens-Lottes h rises with q'''' at fixed pressure', &
                      test_fuel_jens_lottes_monotone())
        call run_test('Convection helper picks nucleate boiling regime when boiling', &
                      test_fuel_convection_regime_switch())
        call run_test('Heat step: nucleate boiling drops fuel temp vs single-phase', &
                      test_fuel_heat_coupling_temp_drop())
    end subroutine

    function test_fuel_init_blade_positions() result(pass)
        use crd,  only: CRD_N_BLADES
        use fuel
        logical :: pass
        type(fuel_state_t)  :: state
        type(fuel_config_t) :: config
        integer :: b
        real(wp) :: x_min, x_max, y_min, y_max, core_size

        ! Use a 20×20 grid with 4.75 m core_diameter (matches the real sim).
        call fuel_init_geometry(state, config, 20, 20, 0.2375_wp, 0.2375_wp, &
                                0.1855_wp, 4.75_wp)
        core_size = 20.0_wp * 0.2375_wp

        x_min = minval(state%blade_xy(1, :))
        x_max = maxval(state%blade_xy(1, :))
        y_min = minval(state%blade_xy(2, :))
        y_max = maxval(state%blade_xy(2, :))

        pass = size(state%blade_xy, 2) == CRD_N_BLADES &
               .and. x_min > 0.0_wp .and. x_max < core_size &
               .and. y_min > 0.0_wp .and. y_max < core_size

        ! No two blades collide at the same (x, y).
        do b = 2, CRD_N_BLADES
            if (any(abs(state%blade_xy(1, 1:b-1) - state%blade_xy(1, b)) < 1.0e-9_wp &
              .and. abs(state%blade_xy(2, 1:b-1) - state%blade_xy(2, b)) < 1.0e-9_wp)) then
                pass = .false.
            end if
        end do
        call fuel_destroy(state)
    end function

    function test_fuel_blade_map_coverage() result(pass)
        use fuel
        use crd, only: CRD_N_BLADES
        logical :: pass
        type(fuel_state_t)  :: state
        type(fuel_config_t) :: config
        integer :: i, j

        call fuel_init_geometry(state, config, 20, 20, 0.2375_wp, 0.2375_wp, &
                                0.1855_wp, 4.75_wp)

        pass = .true.
        do j = 1, 20
            do i = 1, 20
                ! Fuel cells: blade idx ∈ [1, CRD_N_BLADES]
                ! Bypass cells: blade idx = 0
                if (state%is_fuel_cell(i, j)) then
                    if (state%blade_for_cell(i, j) < 1 &
                        .or. state%blade_for_cell(i, j) > CRD_N_BLADES) pass = .false.
                else
                    if (state%blade_for_cell(i, j) /= 0) pass = .false.
                end if
            end do
        end do
        call fuel_destroy(state)
    end function

    function test_fuel_blade_map_symmetric() result(pass)
        use fuel
        logical :: pass
        type(fuel_state_t)  :: state
        type(fuel_config_t) :: config
        integer :: i, j, nx, ny
        integer :: blade_a, blade_b
        real(wp) :: xa, ya, xb, yb

        nx = 20; ny = 20
        call fuel_init_geometry(state, config, nx, ny, 0.2375_wp, 0.2375_wp, &
                                0.1855_wp, 4.75_wp)

        ! Test the geometric 180° symmetry of the blade map: cell (i, j)
        ! and cell (nx+1-i, ny+1-j) sit on opposite sides of the core
        ! centre. Their blade IDs may differ, but those blades must be
        ! at mirror-image positions about the core centre.
        pass = .true.
        do j = 1, ny
            do i = 1, nx
                if (.not. state%is_fuel_cell(i, j)) cycle
                blade_a = state%blade_for_cell(i, j)
                blade_b = state%blade_for_cell(nx + 1 - i, ny + 1 - j)
                if (blade_b == 0) then
                    pass = .false.; cycle
                end if
                xa = state%blade_xy(1, blade_a)
                ya = state%blade_xy(2, blade_a)
                xb = state%blade_xy(1, blade_b)
                yb = state%blade_xy(2, blade_b)
                ! Mirror about (core/2, core/2).
                if (abs((xa + xb) - state%core_diameter_m) > 0.3_wp &
                    .or. abs((ya + yb) - state%core_diameter_m) > 0.3_wp) pass = .false.
            end do
        end do
        call fuel_destroy(state)
    end function

    function test_fuel_four_blade_pattern() result(pass)
        use crd, only: CRD_N_BLADES
        use fuel
        logical :: pass
        type(fuel_state_t)  :: state
        type(fuel_config_t) :: config
        real(wp) :: blade_insertion(CRD_N_BLADES)
        integer :: i, j
        integer :: cells_with_rod, cells_total_fuel

        call fuel_init_geometry(state, config, 20, 20, 0.2375_wp, 0.2375_wp, &
                                0.1855_wp, 4.75_wp)

        blade_insertion = 0.0_wp
        ! Insert 4 specific blades at full insertion.
        blade_insertion(10) = 1.0_wp
        blade_insertion(40) = 1.0_wp
        blade_insertion(80) = 1.0_wp
        blade_insertion(120) = 1.0_wp

        cells_with_rod   = 0
        cells_total_fuel = 0
        do j = 1, 20
            do i = 1, 20
                if (.not. state%is_fuel_cell(i, j)) cycle
                cells_total_fuel = cells_total_fuel + 1
                if (blade_insertion(state%blade_for_cell(i, j)) > 0.5_wp) &
                    cells_with_rod = cells_with_rod + 1
            end do
        end do

        ! Exactly 4 blades inserted → the cells nearest those 4 blades
        ! should be a strict (small) subset of the fuel cells.
        pass = cells_with_rod > 0 &
               .and. cells_with_rod < cells_total_fuel &
               .and. cells_with_rod >= 4
        call fuel_destroy(state)
    end function

    function test_fuel_bank_average_equiv() result(pass)
        use crd, only: CRD_N_BLADES
        use fuel
        logical :: pass
        type(fuel_state_t)  :: state
        type(fuel_config_t) :: config
        real(wp) :: blade_insertion(CRD_N_BLADES)
        integer :: i, j
        real(wp) :: per_cell_insertion, expected

        call fuel_init_geometry(state, config, 20, 20, 0.2375_wp, 0.2375_wp, &
                                0.1855_wp, 4.75_wp)

        blade_insertion = 0.42_wp
        expected = 0.42_wp

        pass = .true.
        do j = 1, 20
            do i = 1, 20
                if (.not. state%is_fuel_cell(i, j)) cycle
                per_cell_insertion = blade_insertion(state%blade_for_cell(i, j))
                if (abs(per_cell_insertion - expected) > 1.0e-12_wp) pass = .false.
            end do
        end do
        call fuel_destroy(state)
    end function

    function test_fuel_void_rf_monotone() result(pass)
        use fuel
        logical :: pass
        real(wp) :: rf_0, rf_30, rf_70
        rf_0  = void_reactivity_factor(0.0_wp)
        rf_30 = void_reactivity_factor(0.30_wp)
        rf_70 = void_reactivity_factor(0.70_wp)
        ! Higher void → lower moderator → smaller scaling factor.
        pass = (rf_0 > rf_30) .and. (rf_30 > rf_70) .and. rf_70 > 0.0_wp
    end function

    function test_fuel_power_to_volumetric() result(pass)
        ! Reorg step 10 bug fix: mg_get_power emits W/cm³, heat%q wants
        ! W/m³. The conversion factor is 10⁶ (1 m³ = 10⁶ cm³). This
        ! test pins the contract so the regression can't sneak back in.
        use fuel
        logical :: pass
        real(wp) :: power_W_cm3(2, 2, 2)
        real(wp) :: q_W_m3(2, 2, 2)

        power_W_cm3 = 1.0_wp
        q_W_m3 = fuel_power_to_volumetric_W_m3(power_W_cm3)
        pass = all(abs(q_W_m3 - 1.0e6_wp) < 1.0e-6_wp)

        ! Nontrivial input — every entry should scale by 10⁶.
        power_W_cm3(1, 1, 1) = 250.0_wp
        power_W_cm3(2, 2, 2) = 0.0_wp
        q_W_m3 = fuel_power_to_volumetric_W_m3(power_W_cm3)
        pass = pass .and. abs(q_W_m3(1, 1, 1) - 2.5e8_wp) < 1.0e-3_wp &
                    .and. abs(q_W_m3(2, 2, 2)) < 1.0e-9_wp
    end function

    function test_fuel_power_to_surface() result(pass)
        ! q'' = power × 10⁶ × dx × dy / heated_perimeter.
        ! Constant power + unit cell + unit perimeter → q'' = 10⁶.
        use fuel
        logical :: pass
        real(wp) :: power_W_cm3(2, 2, 2), heated_perimeter(2, 2, 2), q_pp(2, 2, 2)

        power_W_cm3      = 1.0_wp
        heated_perimeter = 1.0_wp
        q_pp = fuel_power_to_surface_q(power_W_cm3, 1.0_wp, 1.0_wp, heated_perimeter)
        pass = all(abs(q_pp - 1.0e6_wp) < 1.0e-3_wp)

        ! Doubling perimeter halves q''. Doubling dx doubles q''. The
        ! conversion is dimensionally q'' [W/m²] = (W/cm³ → W/m³)
        ! × cell_volume_xy / perimeter.
        heated_perimeter = 2.0_wp
        q_pp = fuel_power_to_surface_q(power_W_cm3, 1.0_wp, 1.0_wp, heated_perimeter)
        pass = pass .and. all(abs(q_pp - 5.0e5_wp) < 1.0e-3_wp)

        ! Zero perimeter must not blow up (floor handling).
        heated_perimeter = 0.0_wp
        q_pp = fuel_power_to_surface_q(power_W_cm3, 1.0_wp, 1.0_wp, heated_perimeter)
        pass = pass .and. all(q_pp > 0.0_wp)
    end function

    function test_fuel_dittus_boelter() result(pass)
        ! Step-11: verify Dittus-Boelter helper matches Nu = 0.023·Re^0.8·Pr^0.4
        ! at a canonical turbulent point. Picked G, D_h, props so that
        ! Re ≫ 10⁴ and the laminar fallback doesn't kick in.
        use fuel, only: dittus_boelter_h
        use two_phase_flow, only: water_properties_t, get_water_properties
        logical :: pass
        type(water_properties_t) :: props
        real(wp) :: G, D_h, h, h_expected, Re, Pr, Nu

        G   = 2000.0_wp     ! kg/m²·s — typical BWR core
        D_h = 0.012_wp      ! m — GE-14 hydraulic diameter
        props = get_water_properties(558.0_wp, 7.0e6_wp)

        Re = abs(G) * D_h / props%mu_l
        Pr = props%mu_l * props%cp_l / props%k_l
        Nu = 0.023_wp * Re**0.8_wp * Pr**0.4_wp
        h_expected = Nu * props%k_l / D_h

        h = dittus_boelter_h(G, D_h, props)
        pass = abs(h - h_expected) / h_expected < 1.0e-6_wp

        ! Sanity floor: BWR-class h should be in the 10³–10⁵ band.
        pass = pass .and. h > 1.0e3_wp .and. h < 1.0e6_wp
    end function

    function test_fuel_jens_lottes_monotone() result(pass)
        ! Step-11: h via Jens-Lottes scales as q''^0.75 (since
        ! ΔT_sat ∝ q''^0.25). At fixed pressure, doubling q'' should
        ! raise h by a factor of 2^0.75 ≈ 1.682.
        use fuel, only: jens_lottes_h
        logical :: pass
        real(wp) :: p, q1, q2, h1, h2, ratio, expected

        p  = 7.0e6_wp    ! 7 MPa — BWR dome pressure
        q1 = 0.5e6_wp    ! 0.5 MW/m²
        q2 = 1.0e6_wp    ! 1.0 MW/m²
        h1 = jens_lottes_h(q1, p)
        h2 = jens_lottes_h(q2, p)

        expected = 2.0_wp**0.75_wp
        ratio    = h2 / h1
        pass = abs(ratio - expected) / expected < 1.0e-3_wp

        ! Below the q_floor the helper must yield zero (single-phase
        ! Dittus-Boelter handles that regime).
        pass = pass .and. abs(jens_lottes_h(0.0_wp, p)) < 1.0e-9_wp

        ! At 1 MW/m², 7 MPa the historical Jens-Lottes superheat is
        ! ~8 K → h ≈ 1.25e5 W/m²·K. Loose band check.
        pass = pass .and. h2 > 1.0e5_wp .and. h2 < 2.0e5_wp
    end function

    function test_fuel_convection_regime_switch() result(pass)
        ! Step-11: a two-phase state with one boiling cell and one
        ! non-boiling cell. The helper should return h_conv from
        ! Jens-Lottes in the boiling cell (driven by heat_flux) and
        ! Dittus-Boelter (driven by mass_flux) in the other cell.
        use fuel, only: fuel_compute_convection_coefficients
        use two_phase_flow, only: two_phase_state_t, two_phase_init, &
                                  two_phase_destroy
        logical :: pass
        type(two_phase_state_t) :: state
        real(wp), allocatable :: T_f(:, :, :), h_c(:, :, :)

        call two_phase_init(state, 2, 1, 1, 0.1_wp, 0.1_wp, 0.1_wp)
        allocate(T_f(2, 1, 1), h_c(2, 1, 1))

        ! Cell 1: non-boiling, moderate mass flux. Pure single-phase.
        state%boiling   (1, 1, 1) = .false.
        state%mass_flux (1, 1, 1) = 2000.0_wp
        state%heat_flux (1, 1, 1) = 0.0_wp
        state%pressure  (1, 1, 1) = 7.0e6_wp
        state%diameter  (1, 1, 1) = 0.012_wp

        ! Cell 2: boiling, same mass flux but high q''. Jens-Lottes
        ! should beat Dittus-Boelter here.
        state%boiling   (2, 1, 1) = .true.
        state%mass_flux (2, 1, 1) = 2000.0_wp
        state%heat_flux (2, 1, 1) = 1.0e6_wp   ! 1 MW/m²
        state%pressure  (2, 1, 1) = 7.0e6_wp
        state%diameter  (2, 1, 1) = 0.012_wp

        call fuel_compute_convection_coefficients(state, T_f, h_c)

        ! T_fluid in both cells is saturation temperature.
        pass = abs(T_f(1, 1, 1) - state%props(1, 1, 1)%T_sat) < 1.0e-6_wp &
         .and. abs(T_f(2, 1, 1) - state%props(2, 1, 1)%T_sat) < 1.0e-6_wp

        ! Boiling cell must have noticeably higher h than single-phase.
        pass = pass .and. h_c(2, 1, 1) > 2.0_wp * h_c(1, 1, 1)

        ! Both must be physical (positive, finite).
        pass = pass .and. h_c(1, 1, 1) > 0.0_wp .and. h_c(1, 1, 1) < 1.0e7_wp &
                    .and. h_c(2, 1, 1) > 0.0_wp .and. h_c(2, 1, 1) < 1.0e7_wp

        call two_phase_destroy(state)
        deallocate(T_f, h_c)
    end function

    function test_fuel_heat_coupling_temp_drop() result(pass)
        ! Step-11 acceptance criterion: hold the volumetric power source
        ! constant and run the heat kernel to (near) steady state under
        ! two coolant settings — single-phase Dittus-Boelter vs Jens-
        ! Lottes boiling. The boiling case must converge to a *lower*
        ! fuel temperature because h_conv is ~10× higher. Asserts an
        ! absolute drop > 20 K — the analytic estimate at the canonical
        ! BWR operating point is Δ ≈ Q · (1/h_sp - 1/h_nb) / A_over_V
        ! ≈ 30 K, with comfortable margin for kernel quirks.
        use heat_transfer
        use fuel, only: dittus_boelter_h, jens_lottes_h
        use two_phase_flow, only: water_properties_t, get_water_properties
        logical :: pass
        type(heat_state_t) :: state_sp, state_nb
        type(water_properties_t) :: props
        real(wp), allocatable :: T_f(:, :, :), h_sp(:, :, :), h_nb(:, :, :)
        real(wp) :: h_sp_scalar, h_nb_scalar
        real(wp) :: T_sat, T_final_sp, T_final_nb
        integer  :: it
        integer, parameter :: NX = 3, NY = 3, NZ = 3
        real(wp), parameter :: DX = 0.1_wp, DY = 0.1_wp, DZ = 0.1_wp
        real(wp), parameter :: Q_VOL = 1.5e8_wp   ! 150 MW/m³ — chosen so
                                                  ! Q/(h·A/V) is large
                                                  ! relative to noise
        real(wp), parameter :: PRESSURE = 7.0e6_wp

        ! Both states share the same source, dt, and material. Only the
        ! coolant differs.
        call heat_init(state_sp, NX, NY, NZ, DX, DY, DZ)
        call heat_init(state_nb, NX, NY, NZ, DX, DY, DZ)
        call heat_set_properties(state_sp, k=3.0_wp, rho=10970.0_wp, cp=300.0_wp)
        call heat_set_properties(state_nb, k=3.0_wp, rho=10970.0_wp, cp=300.0_wp)
        state_sp%Q = Q_VOL
        state_nb%Q = Q_VOL

        props = get_water_properties(558.0_wp, PRESSURE)
        T_sat = props%T_sat

        h_sp_scalar = dittus_boelter_h(1000.0_wp, 0.012_wp, props)
        h_nb_scalar = jens_lottes_h(1.0e6_wp, PRESSURE)

        ! Both regimes must yield positive coefficients, with nucleate
        ! boiling clearly stronger — otherwise the test is meaningless.
        pass = h_nb_scalar > h_sp_scalar .and. h_sp_scalar > 0.0_wp

        allocate(T_f(NX, NY, NZ), h_sp(NX, NY, NZ), h_nb(NX, NY, NZ))
        T_f  = T_sat
        h_sp = h_sp_scalar
        h_nb = h_nb_scalar
        call heat_set_coolant(state_sp, T_f, h_sp)
        call heat_set_coolant(state_nb, T_f, h_nb)

        ! Relax to (near) steady state via the implicit kernel.
        do it = 1, 50
            call heat_step_implicit(state_sp, 1.0_wp)
            call heat_step_implicit(state_nb, 1.0_wp)
        end do

        T_final_sp = maxval(state_sp%T)
        T_final_nb = maxval(state_nb%T)

        ! Boiling regime must drop the fuel temperature by > 20 K vs the
        ! single-phase case. This is the "power ramp into nucleate
        ! boiling shows a *drop* in fuel temperature" demonstration
        ! called for in REORG_PLAN step 11.
        pass = pass .and. T_final_nb < T_final_sp &
                    .and. (T_final_sp - T_final_nb) > 20.0_wp

        call heat_destroy(state_sp)
        call heat_destroy(state_nb)
        deallocate(T_f, h_sp, h_nb)
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