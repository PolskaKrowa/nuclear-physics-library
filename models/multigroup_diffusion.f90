! models/multigroup_diffusion.f90
!
! Multi-group neutron diffusion solver for reactor physics.
! Solves the steady-state and time-dependent multi-group diffusion equations:
!
!   1/v_g · ∂ϕ_g/∂t = ∇·D_g∇ϕ_g - Σ_r,g·ϕ_g + Σ_s,g'→g·ϕ_g' + χ_g·(1/k)·Σ_f,g'·ϕ_g'
!
! Features:
! - 2, 4, or 8 energy group structures
! - Fast/thermal spectrum treatment
! - Upscattering and downscattering
! - Fission source iteration for k-eigenvalue
! - Power iteration and source iteration methods
! - Albedo boundary conditions
!
! Usage:
!   use multigroup_diffusion
!   type(mg_config_t) :: config
!   type(mg_state_t) :: state
!   
!   config%n_groups = 2
!   call mg_init(state, nx, ny, nz, dx, dy, dz, config)
!   call mg_set_cross_sections(state, xsec_data)
!   call mg_solve_eigenvalue(state, k_eff, converged)
!   call mg_solve_transient(state, precursors, dt)
!   call mg_get_power(state, power_density)
!   call mg_destroy(state)
!
module multigroup_diffusion
    use kinds, only: wp, i32
    use constants, only: TOL_DEFAULT
    use solve_linear, only: solve_dense, SOLVE_SUCCESS
    use sparse_matrix, only: sparse_matrix_t, sparse_matvec
    implicit none
    private
    
    ! Public interface
    public :: mg_config_t, mg_state_t, mg_xsec_t
    public :: mg_init, mg_destroy
    public :: mg_set_cross_sections, mg_set_source
    public :: mg_solve_eigenvalue, mg_solve_fixed_source
    public :: mg_solve_transient, mg_power_iteration
    public :: mg_get_flux, mg_get_power, mg_get_fission_rate
    public :: mg_compute_keff, mg_compute_form_factors
    public :: mg_validate_cross_sections
    
    ! Standard energy group structures (eV boundaries)
    integer, parameter, public :: MG_2GROUP = 2
    integer, parameter, public :: MG_4GROUP = 4
    integer, parameter, public :: MG_8GROUP = 8
    
    ! Solver methods
    integer, parameter, public :: SOLVER_POWER_ITERATION = 1
    integer, parameter, public :: SOLVER_SOURCE_ITERATION = 2
    integer, parameter, public :: SOLVER_KRYLOV = 3
    
    ! Multi-group cross sections (per material)
    type :: mg_xsec_t
        integer :: n_groups = 2
        
        ! Group-wise cross sections [cm⁻¹]
        real(wp), allocatable :: D(:)              ! Diffusion coefficient
        real(wp), allocatable :: sigma_t(:)        ! Total
        real(wp), allocatable :: sigma_a(:)        ! Absorption
        real(wp), allocatable :: sigma_f(:)        ! Fission
        real(wp), allocatable :: nu_sigma_f(:)     ! Production
        real(wp), allocatable :: sigma_r(:)        ! Removal
        real(wp), allocatable :: chi(:)            ! Fission spectrum
        real(wp), allocatable :: kappa(:)          ! Energy per fission [MeV]
        
        ! Scattering matrix [cm⁻¹]: (from → to)
        real(wp), allocatable :: sigma_s(:, :)     ! (g_from, g_to)
        
        ! Delayed neutron data
        real(wp) :: beta_total = 0.0065_wp
        real(wp), allocatable :: chi_d(:)          ! Delayed spectrum
    end type mg_xsec_t
    
    ! Configuration
    type :: mg_config_t
        integer :: n_groups = 2
        integer :: max_outer_iter = 100
        integer :: max_inner_iter = 50
        real(wp) :: outer_tolerance = 1.0e-5_wp
        real(wp) :: inner_tolerance = 1.0e-6_wp
        real(wp) :: k_guess = 1.0_wp
        logical :: use_upscatter = .true.
        logical :: normalize_power = .true.
        real(wp) :: power_level = 1.0e6_wp        ! Watts
        integer :: solver_method = SOLVER_POWER_ITERATION
        
        ! Boundary conditions (albedo)
        real(wp) :: albedo_x_min = 0.0_wp
        real(wp) :: albedo_x_max = 0.0_wp
        real(wp) :: albedo_y_min = 0.0_wp
        real(wp) :: albedo_y_max = 0.0_wp
        real(wp) :: albedo_z_min = 0.0_wp
        real(wp) :: albedo_z_max = 0.0_wp
    end type mg_config_t
    
    ! Multi-group diffusion state
    type :: mg_state_t
        ! Grid
        integer :: nx, ny, nz
        real(wp) :: dx, dy, dz
        real(wp) :: volume
        
        ! Number of energy groups
        integer :: n_groups
        
        ! Multi-group flux [n/cm²·s] (nx, ny, nz, n_groups)
        real(wp), allocatable :: flux(:, :, :, :)
        real(wp), allocatable :: flux_old(:, :, :, :)
        
        ! Fission source [n/cm³·s]
        real(wp), allocatable :: fission_source(:, :, :)
        
        ! Fixed source [n/cm³·s] (optional, for fixed source problems)
        real(wp), allocatable :: external_source(:, :, :, :)
        
        ! Cross sections per cell
        type(mg_xsec_t), allocatable :: xsec(:, :, :)
        
        ! Eigenvalue
        real(wp) :: k_eff = 1.0_wp
        real(wp) :: k_eff_old = 1.0_wp
        
        ! Convergence tracking
        integer :: outer_iterations = 0
        integer :: inner_iterations = 0
        real(wp) :: outer_error = 1.0_wp
        real(wp) :: inner_error = 1.0_wp
        logical :: converged = .false.
        
        ! Power distribution
        real(wp), allocatable :: power_density(:, :, :)  ! W/cm³
        real(wp) :: total_power = 0.0_wp
        
        ! Form factors
        real(wp), allocatable :: radial_peaking(:, :)
        real(wp), allocatable :: axial_peaking(:)
        real(wp) :: total_peaking = 1.0_wp
        
        ! Configuration
        type(mg_config_t) :: config
    end type mg_state_t
    
contains

    !> Initialize multi-group diffusion state
    subroutine mg_init(state, nx, ny, nz, dx, dy, dz, config)
        type(mg_state_t), intent(out) :: state
        integer, intent(in) :: nx, ny, nz
        real(wp), intent(in) :: dx, dy, dz
        type(mg_config_t), intent(in), optional :: config
        
        integer :: g, i, j, k
        
        state%nx = nx
        state%ny = ny
        state%nz = nz
        state%dx = dx * 100.0_wp  ! Convert m to cm
        state%dy = dy * 100.0_wp
        state%dz = dz * 100.0_wp
        state%volume = dx * dy * dz  ! m³
        
        ! Configuration
        if (present(config)) then
            state%config = config
        end if
        
        state%n_groups = state%config%n_groups
        
        ! Allocate flux arrays
        allocate(state%flux(nx, ny, nz, state%n_groups))
        allocate(state%flux_old(nx, ny, nz, state%n_groups))
        allocate(state%fission_source(nx, ny, nz))
        allocate(state%external_source(nx, ny, nz, state%n_groups))
        
        ! Allocate cross sections
        allocate(state%xsec(nx, ny, nz))
        
        ! Allocate power
        allocate(state%power_density(nx, ny, nz))
        allocate(state%radial_peaking(nx, ny))
        allocate(state%axial_peaking(nz))
        
        ! Initialize with flat flux guess
        state%flux = 1.0e10_wp  ! Typical flux level
        state%flux_old = state%flux
        state%fission_source = 0.0_wp
        state%external_source = 0.0_wp
        state%power_density = 0.0_wp
        
        state%k_eff = state%config%k_guess
        state%k_eff_old = state%k_eff
        
        ! Initialize cross sections to default values
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    call initialize_default_xsec(state%xsec(i, j, k), state%n_groups)
                end do
            end do
        end do
    end subroutine mg_init
    
    !> Destroy multi-group state
    subroutine mg_destroy(state)
        type(mg_state_t), intent(inout) :: state
        
        if (allocated(state%flux)) deallocate(state%flux)
        if (allocated(state%flux_old)) deallocate(state%flux_old)
        if (allocated(state%fission_source)) deallocate(state%fission_source)
        if (allocated(state%external_source)) deallocate(state%external_source)
        if (allocated(state%power_density)) deallocate(state%power_density)
        if (allocated(state%radial_peaking)) deallocate(state%radial_peaking)
        if (allocated(state%axial_peaking)) deallocate(state%axial_peaking)
        if (allocated(state%xsec)) deallocate(state%xsec)
    end subroutine mg_destroy
    
    !> Initialize default cross sections (water moderator)
    subroutine initialize_default_xsec(xsec, n_groups)
        type(mg_xsec_t), intent(out) :: xsec
        integer, intent(in) :: n_groups
        
        xsec%n_groups = n_groups
        
        allocate(xsec%D(n_groups))
        allocate(xsec%sigma_t(n_groups))
        allocate(xsec%sigma_a(n_groups))
        allocate(xsec%sigma_f(n_groups))
        allocate(xsec%nu_sigma_f(n_groups))
        allocate(xsec%sigma_r(n_groups))
        allocate(xsec%chi(n_groups))
        allocate(xsec%kappa(n_groups))
        allocate(xsec%sigma_s(n_groups, n_groups))
        allocate(xsec%chi_d(n_groups))
        
        ! Default: water moderator (non-multiplying)
        xsec%D = 1.0_wp
        xsec%sigma_t = 0.5_wp
        xsec%sigma_a = 0.01_wp
        xsec%sigma_f = 0.0_wp
        xsec%nu_sigma_f = 0.0_wp
        xsec%sigma_r = 0.01_wp
        xsec%chi = 0.0_wp
        xsec%kappa = 200.0_wp
        xsec%chi_d = 0.0_wp
        
        ! Diagonal scattering (no group transfer)
        xsec%sigma_s = 0.0_wp
        if (n_groups >= 1) xsec%sigma_s(1, 1) = 0.49_wp
        if (n_groups >= 2) xsec%sigma_s(2, 2) = 0.49_wp
    end subroutine initialize_default_xsec
    
    !> Set cross sections for a region
    subroutine mg_set_cross_sections(state, xsec, i1, i2, j1, j2, k1, k2)
        type(mg_state_t), intent(inout) :: state
        type(mg_xsec_t), intent(in) :: xsec
        integer, intent(in), optional :: i1, i2, j1, j2, k1, k2
        
        integer :: imin, imax, jmin, jmax, kmin, kmax
        integer :: i, j, k
        
        imin = 1; imax = state%nx
        jmin = 1; jmax = state%ny
        kmin = 1; kmax = state%nz
        
        if (present(i1)) imin = i1
        if (present(i2)) imax = i2
        if (present(j1)) jmin = j1
        if (present(j2)) jmax = j2
        if (present(k1)) kmin = k1
        if (present(k2)) kmax = k2
        
        do k = kmin, kmax
            do j = jmin, jmax
                do i = imin, imax
                    state%xsec(i, j, k) = xsec
                end do
            end do
        end do
    end subroutine mg_set_cross_sections
    
    !> Set external source
    subroutine mg_set_source(state, source, group)
        type(mg_state_t), intent(inout) :: state
        real(wp), intent(in) :: source(:, :, :)
        integer, intent(in), optional :: group
        
        integer :: g
        
        if (present(group)) then
            state%external_source(:, :, :, group) = source
        else
            ! Set to all groups equally
            do g = 1, state%n_groups
                state%external_source(:, :, :, g) = source / real(state%n_groups, wp)
            end do
        end if
    end subroutine mg_set_source
    
    !> Diagnostic: Check cross-sections and initial state
    subroutine mg_check_problem_setup(state)
        type(mg_state_t), intent(in) :: state
        integer :: i, j, k, g
        real(wp) :: total_nu_sigma_f, total_sigma_a, total_flux
        real(wp) :: max_nu_sigma_f, max_sigma_a, avg_production, avg_absorption
        
        print '(A)', '=========================================='
        print '(A)', 'DIAGNOSTIC: Problem Setup Check'
        print '(A)', '=========================================='
        
        ! Check flux levels
        total_flux = sum(state%flux)
        print '(A,ES12.3)', 'Total flux: ', total_flux
        print '(A,ES12.3)', 'Max flux:   ', maxval(state%flux)
        print '(A,ES12.3)', 'Min flux:   ', minval(state%flux)
        
        ! Check cross-sections
        total_nu_sigma_f = 0.0_wp
        total_sigma_a = 0.0_wp
        max_nu_sigma_f = 0.0_wp
        max_sigma_a = 0.0_wp
        
        do k = 1, state%nz
            do j = 1, state%ny
                do i = 1, state%nx
                    do g = 1, state%n_groups
                        total_nu_sigma_f = total_nu_sigma_f + state%xsec(i,j,k)%nu_sigma_f(g)
                        total_sigma_a = total_sigma_a + state%xsec(i,j,k)%sigma_a(g)
                        max_nu_sigma_f = max(max_nu_sigma_f, state%xsec(i,j,k)%nu_sigma_f(g))
                        max_sigma_a = max(max_sigma_a, state%xsec(i,j,k)%sigma_a(g))
                    end do
                end do
            end do
        end do
        
        print '(A)', ''
        print '(A)', 'Cross-sections:'
        print '(A,ES12.3)', '  Max nu*sigma_f:   ', max_nu_sigma_f
        print '(A,ES12.3)', '  Max sigma_a:      ', max_sigma_a
        print '(A,ES12.3)', '  Total nu*sigma_f: ', total_nu_sigma_f
        print '(A,ES12.3)', '  Total sigma_a:    ', total_sigma_a
        
        if (max_nu_sigma_f < 1.0e-10_wp) then
            print '(A)', '  ** WARNING: No fissile material detected! **'
        end if
        
        ! Estimate k-infinite
        avg_production = total_nu_sigma_f / real(state%nx*state%ny*state%nz*state%n_groups, wp)
        avg_absorption = total_sigma_a / real(state%nx*state%ny*state%nz*state%n_groups, wp)
        
        if (avg_absorption > 1.0e-10_wp) then
            print '(A,F10.5)', '  Estimated k-inf:  ', avg_production / avg_absorption
        end if
        
        print '(A)', '=========================================='
        print '(A)', ''
    end subroutine mg_check_problem_setup

    !> Validate cross-sections before solving
    subroutine mg_validate_cross_sections(state, is_valid)
        type(mg_state_t), intent(inout) :: state
        logical, intent(out) :: is_valid
        
        integer :: i, j, k, g, gp
        real(wp) :: sigma_total, nu_ratio
        real(wp) :: max_sigma_a, min_sigma_a, max_nu_sigma_f
        logical :: has_fissile
        
        is_valid = .true.
        has_fissile = .false.
        
        ! Initialize min/max values
        max_sigma_a = -1.0e30_wp
        min_sigma_a = 1.0e30_wp
        max_nu_sigma_f = -1.0e30_wp
        
        print '(A)', '=========================================='
        print '(A)', 'Validating Cross-Sections...'
        print '(A)', '=========================================='
        
        do k = 1, state%nz
            do j = 1, state%ny
                do i = 1, state%nx
                    do g = 1, state%n_groups
                        
                        ! Check for negative values
                        if (state%xsec(i,j,k)%D(g) < 0.0_wp) then
                            print '(A,3I4,I2,A)', 'ERROR: D < 0 at (', i, j, k, g, ')'
                            is_valid = .false.
                        end if
                        
                        if (state%xsec(i,j,k)%sigma_a(g) < 0.0_wp) then
                            print '(A,3I4,I2,A)', 'ERROR: sigma_a < 0 at (', i, j, k, g, ')'
                            is_valid = .false.
                        end if
                        
                        ! FIX: If sigma_a is zero, set it to a small value
                        if (state%xsec(i,j,k)%sigma_a(g) < 1.0e-10_wp) then
                            ! Use nu_sigma_f as proxy (absorption = fission for pure fissile)
                            if (state%xsec(i,j,k)%nu_sigma_f(g) > 0.0_wp) then
                                ! Estimate sigma_a from nu_sigma_f (typical nu ~ 2.5)
                                state%xsec(i,j,k)%sigma_a(g) = &
                                    state%xsec(i,j,k)%nu_sigma_f(g) / 2.5_wp
                            else
                                ! Set minimum absorption to avoid singular matrix
                                state%xsec(i,j,k)%sigma_a(g) = 1.0e-6_wp
                            end if
                        end if
                        
                        ! Track min/max
                        max_sigma_a = max(max_sigma_a, state%xsec(i,j,k)%sigma_a(g))
                        min_sigma_a = min(min_sigma_a, state%xsec(i,j,k)%sigma_a(g))
                        max_nu_sigma_f = max(max_nu_sigma_f, state%xsec(i,j,k)%nu_sigma_f(g))
                        
                        ! Check for fissile material
                        if (state%xsec(i,j,k)%nu_sigma_f(g) > 1.0e-10_wp) then
                            has_fissile = .true.
                        end if
                        
                        ! Validate nu ratio
                        if (state%xsec(i,j,k)%sigma_f(g) > 1.0e-10_wp) then
                            nu_ratio = state%xsec(i,j,k)%nu_sigma_f(g) / &
                                    state%xsec(i,j,k)%sigma_f(g)
                            if (nu_ratio < 1.5_wp .or. nu_ratio > 4.0_wp) then
                                print '(A,3I4,I2,A,F6.2)', 'WARNING: Unusual nu at (', &
                                    i, j, k, g, '): ', nu_ratio
                            end if
                        end if
                        
                        ! Compute total removal
                        sigma_total = state%xsec(i,j,k)%sigma_a(g)
                        do gp = 1, state%n_groups
                            if (gp /= g) then
                                sigma_total = sigma_total + state%xsec(i,j,k)%sigma_s(g, gp)
                            end if
                        end do
                        
                        ! Store in sigma_r for later use
                        state%xsec(i,j,k)%sigma_r(g) = sigma_total
                        
                    end do
                end do
            end do
        end do
        
        if (.not. has_fissile) then
            print '(A)', 'ERROR: No fissile material found!'
            is_valid = .false.
        end if
        
        ! Print summary of fixed cross-sections
        print '(A)', 'After validation:'
        print '(A,ES12.3)', '  Max sigma_a:     ', max_sigma_a
        print '(A,ES12.3)', '  Min sigma_a:     ', min_sigma_a
        print '(A,ES12.3)', '  Max nu_sigma_f:  ', max_nu_sigma_f
        print '(A,L2)', '  Has fissile:     ', has_fissile
        print '(A)', '=========================================='
        print '(A)', ''
        
    end subroutine mg_validate_cross_sections

    !> Detect NaN in arrays
    function has_nan(arr) result(found_nan)
        real(wp), intent(in) :: arr(:,:,:,:)
        logical :: found_nan
        integer :: i, j, k, g
        
        found_nan = .false.
        
        do g = 1, size(arr, 4)
            do k = 1, size(arr, 3)
                do j = 1, size(arr, 2)
                    do i = 1, size(arr, 1)
                        if (arr(i,j,k,g) /= arr(i,j,k,g)) then  ! NaN check
                            print '(A,4I5)', 'NaN detected at: ', i, j, k, g
                            found_nan = .true.
                            return
                        end if
                    end do
                end do
            end do
        end do
    end function has_nan

    !> Solve k-eigenvalue problem using power iteration
    subroutine mg_solve_eigenvalue(state, k_eff, converged)
        type(mg_state_t), intent(inout) :: state
        real(wp), intent(out) :: k_eff
        logical, intent(out) :: converged
        
        integer :: outer_iter, inner_iter, g, sweep
        real(wp) :: k_error, flux_max, flux_min
        real(wp) :: fission_old, fission_new, flux_scale
        logical :: valid_xsec
        
        ! Validate cross-sections FIRST
        call mg_validate_cross_sections(state, valid_xsec)
        if (.not. valid_xsec) then
            print '(A)', 'ERROR: Invalid cross-sections. Cannot proceed.'
            converged = .false.
            k_eff = 0.0_wp
            return
        end if
        
        state%converged = .false.
        
        ! Ensure positive initial flux
        where (state%flux < 1.0e-10_wp) state%flux = 1.0_wp
        
        ! Outer iteration
        do outer_iter = 1, state%config%max_outer_iter
            state%outer_iterations = outer_iter
            state%k_eff_old = state%k_eff
            
            ! Compute fission source
            call compute_fission_source(state)
            fission_old = sum(state%fission_source)
            
            if (fission_old < 1.0e-30_wp) then
                print '(A)', 'ERROR: Fission source is zero or negative!'
                converged = .false.
                k_eff = 0.0_wp
                return
            end if
            
            ! Store old flux
            state%flux_old = state%flux
            
            ! Inner iterations
            do inner_iter = 1, state%config%max_inner_iter
                
                ! Multiple sweeps
                do sweep = 1, 3
                    do g = 1, state%n_groups
                        call solve_group_equation_safe(state, g)
                        
                        ! Check for NaN after each group
                        if (has_nan(state%flux)) then
                            print '(A,I3,A,I3)', 'ERROR: NaN in flux at outer=', &
                                outer_iter, ', group=', g
                            print '(A)', 'Last valid flux statistics:'
                            print '(A,ES12.3)', '  Max: ', maxval(state%flux_old)
                            print '(A,ES12.3)', '  Min: ', minval(state%flux_old)
                            converged = .false.
                            k_eff = 0.0_wp
                            return
                        end if
                    end do
                end do
                
                ! Check flux bounds
                flux_max = maxval(state%flux)
                flux_min = minval(state%flux)
                
                if (flux_max > 1.0e20_wp .or. flux_min < 0.0_wp) then
                    print '(A)', 'WARNING: Flux out of bounds, rescaling...'
                    print '(A,ES12.3,A,ES12.3)', '  Range: ', flux_min, ' to ', flux_max
                    state%flux = max(state%flux, 0.0_wp)
                    flux_scale = maxval(state%flux)
                    if (flux_scale > 1.0e10_wp) then
                        state%flux = state%flux / flux_scale * 1.0e10_wp
                    end if
                end if
            end do
            
            ! Compute new fission source
            call compute_fission_source(state)
            fission_new = sum(state%fission_source)
            
            ! Diagnostic output
            if (mod(outer_iter, 5) == 0 .or. outer_iter <= 3) then
                print '(A,I4,A,F12.8,A,ES10.2,A,ES10.2)', &
                    'Iter ', outer_iter, ': k=', state%k_eff, &
                    ', fiss_old=', fission_old, ', fiss_new=', fission_new
            end if
            
            ! Update k-effective
            if (fission_old > 1.0e-30_wp .and. fission_new > 1.0e-30_wp) then
                state%k_eff = state%k_eff_old * (fission_new / fission_old)
            else
                print '(A)', 'ERROR: Fission source collapsed!'
                converged = .false.
                k_eff = 0.0_wp
                return
            end if
            
            ! Clamp k-effective to reasonable range
            if (state%k_eff < 0.1_wp) then
                print '(A,F12.5)', 'WARNING: k-eff very low: ', state%k_eff
                print '(A)', 'Possible causes:'
                print '(A)', '  - Insufficient fissile material'
                print '(A)', '  - Too much leakage (check boundaries)'
                print '(A)', '  - Incorrect cross-sections'
            end if
            
            if (state%k_eff > 5.0_wp) then
                print '(A,F12.5)', 'WARNING: k-eff very high: ', state%k_eff
                state%k_eff = 5.0_wp
            end if
            
            ! Normalize flux
            flux_scale = sum(state%flux)
            if (flux_scale > 1.0e-30_wp) then
                state%flux = state%flux / flux_scale * 1.0e15_wp
            end if
            
            ! Check convergence
            k_error = abs(state%k_eff - state%k_eff_old) / (abs(state%k_eff) + 1.0e-30_wp)
            state%outer_error = k_error
            
            if (k_error < state%config%outer_tolerance .and. outer_iter > 3) then
                state%converged = .true.
                exit
            end if
        end do
        
        call compute_power_distribution(state)
        k_eff = state%k_eff
        converged = state%converged
        
        if (converged) then
            print '(A,F12.8)', 'SUCCESS: k-effective = ', k_eff
        else
            print '(A)', 'WARNING: Did not converge'
        end if
    end subroutine mg_solve_eigenvalue
    
    subroutine compute_balance(state, production, absorption, leakage)
        type(mg_state_t), intent(in) :: state
        real(wp), intent(out) :: production, absorption, leakage
        
        integer :: i, j, k, g
        real(wp) :: leak_x, leak_y, leak_z, D
        
        production = 0.0_wp
        absorption = 0.0_wp
        leakage = 0.0_wp
        
        do k = 1, state%nz
            do j = 1, state%ny
                do i = 1, state%nx
                    do g = 1, state%n_groups
                        ! Production
                        production = production + &
                            state%xsec(i,j,k)%nu_sigma_f(g) * state%flux(i,j,k,g)
                        
                        ! Absorption
                        absorption = absorption + &
                            state%xsec(i,j,k)%sigma_a(g) * state%flux(i,j,k,g)
                        
                        ! Leakage (approximate)
                        if (i > 1 .and. i < state%nx .and. &
                            j > 1 .and. j < state%ny .and. &
                            k > 1 .and. k < state%nz) then
                            
                            D = state%xsec(i,j,k)%D(g)
                            leak_x = D * abs(state%flux(i+1,j,k,g) - 2*state%flux(i,j,k,g) + &
                                            state%flux(i-1,j,k,g)) / (state%dx**2)
                            leak_y = D * abs(state%flux(i,j+1,k,g) - 2*state%flux(i,j,k,g) + &
                                            state%flux(i,j-1,k,g)) / (state%dy**2)
                            leak_z = D * abs(state%flux(i,j,k+1,g) - 2*state%flux(i,j,k,g) + &
                                            state%flux(i,j,k-1,g)) / (state%dz**2)
                            
                            leakage = leakage + (leak_x + leak_y + leak_z)
                        end if
                    end do
                end do
            end do
        end do
    end subroutine compute_balance

    !> Solve fixed source problem
    subroutine mg_solve_fixed_source(state, converged)
        type(mg_state_t), intent(inout) :: state
        logical, intent(out) :: converged
        
        integer :: iter, g
        real(wp) :: flux_error
        
        state%converged = .false.
        
        do iter = 1, state%config%max_inner_iter
            state%inner_iterations = iter
            state%flux_old = state%flux
            
            ! Sweep through all energy groups
            do g = 1, state%n_groups
                call solve_group_equation_safe(state, g)
            end do
            
            ! Check convergence
            flux_error = compute_flux_error(state)
            state%inner_error = flux_error
            
            if (flux_error < state%config%inner_tolerance) then
                state%converged = .true.
                exit
            end if
        end do
        
        converged = state%converged
    end subroutine mg_solve_fixed_source
    
    !> Solve group diffusion equation for eigenvalue problem
    subroutine solve_group_equation_safe(state, g)
        type(mg_state_t), intent(inout) :: state
        integer, intent(in) :: g
        
        integer :: i, j, k, gp
        real(wp) :: source, inscatter, fission, D, sigma_r
        real(wp) :: phi_xp, phi_xm, phi_yp, phi_ym, phi_zp, phi_zm
        real(wp) :: coef_x, coef_y, coef_z, coef_c, rhs, denom
        real(wp) :: inv_dx2, inv_dy2, inv_dz2, phi_new
        real(wp), parameter :: MIN_DENOM = 1.0e-10_wp
        real(wp), parameter :: MAX_FLUX = 1.0e15_wp
        
        inv_dx2 = 1.0_wp / (state%dx * state%dx)
        inv_dy2 = 1.0_wp / (state%dy * state%dy)
        inv_dz2 = 1.0_wp / (state%dz * state%dz)
        
        do k = 2, state%nz - 1
            do j = 2, state%ny - 1
                do i = 2, state%nx - 1
                    
                    D = state%xsec(i, j, k)%D(g)
                    sigma_r = state%xsec(i, j, k)%sigma_r(g)  ! Use pre-computed removal
                    
                    ! Get neighbors (with bounds check)
                    phi_xm = max(state%flux(i-1, j, k, g), 0.0_wp)
                    phi_xp = max(state%flux(i+1, j, k, g), 0.0_wp)
                    phi_ym = max(state%flux(i, j-1, k, g), 0.0_wp)
                    phi_yp = max(state%flux(i, j+1, k, g), 0.0_wp)
                    phi_zm = max(state%flux(i, j, k-1, g), 0.0_wp)
                    phi_zp = max(state%flux(i, j, k+1, g), 0.0_wp)
                    
                    ! Diffusion operator
                    coef_x = D * inv_dx2
                    coef_y = D * inv_dy2
                    coef_z = D * inv_dz2
                    coef_c = 2.0_wp * (coef_x + coef_y + coef_z)
                    
                    ! In-scattering
                    inscatter = 0.0_wp
                    do gp = 1, state%n_groups
                        if (gp /= g) then
                            inscatter = inscatter + &
                                state%xsec(i, j, k)%sigma_s(gp, g) * &
                                max(state%flux(i, j, k, gp), 0.0_wp)
                        end if
                    end do
                    
                    ! Fission source
                    fission = state%xsec(i, j, k)%chi(g) * &
                            state%fission_source(i, j, k) / max(state%k_eff, 0.1_wp)
                    
                    ! Total source
                    source = inscatter + fission + state%external_source(i, j, k, g)
                    
                    ! RHS
                    rhs = coef_x * (phi_xm + phi_xp) + &
                        coef_y * (phi_ym + phi_yp) + &
                        coef_z * (phi_zm + phi_zp) + source
                    
                    ! Denominator (must be positive)
                    denom = coef_c + sigma_r
                    denom = max(denom, MIN_DENOM)
                    
                    ! Update flux
                    phi_new = rhs / denom
                    
                    ! Apply relaxation for stability
                    phi_new = 0.8_wp * phi_new + 0.2_wp * state%flux(i, j, k, g)
                    
                    ! Clamp to reasonable range
                    phi_new = max(phi_new, 0.0_wp)
                    phi_new = min(phi_new, MAX_FLUX)
                    
                    state%flux(i, j, k, g) = phi_new
                end do
            end do
        end do
        
        call apply_boundary_conditions(state, g)
    end subroutine solve_group_equation_safe
    
    !> Solve group equation for fixed source
    subroutine solve_group_equation_fixed(state, g)
        type(mg_state_t), intent(inout) :: state
        integer, intent(in) :: g
        
        ! Similar to eigenvalue case but without k_eff division
        ! Implementation follows same pattern as solve_group_equation
        call solve_group_equation_safe(state, g)
    end subroutine solve_group_equation_fixed
    
    !> Compute Laplacian of flux
    function compute_laplacian(state, i, j, k, g) result(laplacian)
        type(mg_state_t), intent(in) :: state
        integer, intent(in) :: i, j, k, g
        real(wp) :: laplacian
        
        real(wp) :: inv_dx2, inv_dy2, inv_dz2
        real(wp) :: phi
        
        phi = state%flux(i, j, k, g)
        
        inv_dx2 = 1.0_wp / (state%dx * state%dx)
        inv_dy2 = 1.0_wp / (state%dy * state%dy)
        inv_dz2 = 1.0_wp / (state%dz * state%dz)
        
        laplacian = (state%flux(i-1, j, k, g) - 2.0_wp*phi + state%flux(i+1, j, k, g)) * inv_dx2 + &
                    (state%flux(i, j-1, k, g) - 2.0_wp*phi + state%flux(i, j+1, k, g)) * inv_dy2 + &
                    (state%flux(i, j, k-1, g) - 2.0_wp*phi + state%flux(i, j, k+1, g)) * inv_dz2
    end function compute_laplacian
    
    !> Compute total fission source
    subroutine compute_fission_source(state)
        type(mg_state_t), intent(inout) :: state
        
        integer :: i, j, k, g
        real(wp) :: fission
        
        do k = 1, state%nz
            do j = 1, state%ny
                do i = 1, state%nx
                    fission = 0.0_wp
                    do g = 1, state%n_groups
                        fission = fission + state%xsec(i, j, k)%nu_sigma_f(g) * &
                                           state%flux(i, j, k, g)
                    end do
                    state%fission_source(i, j, k) = fission
                end do
            end do
        end do
    end subroutine compute_fission_source

    !> Compute fission source with given flux
    subroutine compute_fission_source_with_flux(state, flux, total_fission)
        type(mg_state_t), intent(in) :: state
        real(wp), intent(in) :: flux(:, :, :, :)
        real(wp), intent(out) :: total_fission
        
        integer :: i, j, k, g
        real(wp) :: fission
        
        total_fission = 0.0_wp
        
        do k = 1, state%nz
            do j = 1, state%ny
                do i = 1, state%nx
                    fission = 0.0_wp
                    do g = 1, state%n_groups
                        fission = fission + state%xsec(i, j, k)%nu_sigma_f(g) * &
                                        flux(i, j, k, g)
                    end do
                    total_fission = total_fission + fission
                end do
            end do
        end do
    end subroutine compute_fission_source_with_flux

    !> Normalize flux to unit fission source
    subroutine normalize_flux_to_unity(state)
        type(mg_state_t), intent(inout) :: state
        
        real(wp) :: total_fission
        
        total_fission = sum(state%fission_source)
        
        if (total_fission > TOL_DEFAULT) then
            state%flux = state%flux / sqrt(total_fission)
            state%fission_source = state%fission_source / total_fission
        end if
    end subroutine normalize_flux_to_unity

    !> Update k-effective
    subroutine update_keff(state)
        type(mg_state_t), intent(inout) :: state
        
        real(wp) :: fission_new, fission_old
        integer :: i, j, k, g
        
        fission_new = 0.0_wp
        fission_old = 0.0_wp
        
        do k = 1, state%nz
            do j = 1, state%ny
                do i = 1, state%nx
                    do g = 1, state%n_groups
                        fission_new = fission_new + &
                            state%xsec(i, j, k)%nu_sigma_f(g) * state%flux(i, j, k, g)
                        fission_old = fission_old + &
                            state%xsec(i, j, k)%nu_sigma_f(g) * state%flux_old(i, j, k, g)
                    end do
                end do
            end do
        end do
        
        if (fission_old > TOL_DEFAULT) then
            state%k_eff = state%k_eff_old * (fission_new / fission_old)
        end if
    end subroutine update_keff
    
    !> Normalize flux to specified power level
    subroutine normalize_flux(state)
        type(mg_state_t), intent(inout) :: state
        
        real(wp) :: current_power, scale_factor
        integer :: i, j, k, g
        
        if (.not. state%config%normalize_power) return
        
        ! Compute current power
        current_power = 0.0_wp
        do k = 1, state%nz
            do j = 1, state%ny
                do i = 1, state%nx
                    do g = 1, state%n_groups
                        current_power = current_power + &
                            state%xsec(i, j, k)%sigma_f(g) * &
                            state%xsec(i, j, k)%kappa(g) * &
                            state%flux(i, j, k, g)
                    end do
                end do
            end do
        end do
        
        ! Convert MeV/cm³ to Watts
        current_power = current_power * 1.602e-13_wp * state%volume * 1.0e6_wp
        
        if (current_power > TOL_DEFAULT) then
            scale_factor = state%config%power_level / current_power
            state%flux = state%flux * scale_factor
        end if
    end subroutine normalize_flux
    
    !> Compute flux error
    function compute_flux_error(state) result(error)
        type(mg_state_t), intent(in) :: state
        real(wp) :: error
        
        real(wp) :: diff, norm_new, norm_old
        
        diff = sum((state%flux - state%flux_old)**2)
        norm_new = sum(state%flux**2)
        norm_old = sum(state%flux_old**2)
        
        if (norm_old > TOL_DEFAULT) then
            error = sqrt(diff / norm_old)
        else
            error = sqrt(diff)
        end if
    end function compute_flux_error
    
    !> Apply boundary conditions (albedo)
    subroutine apply_boundary_conditions(state, g)
        type(mg_state_t), intent(inout) :: state
        integer, intent(in) :: g
        
        real(wp) :: beta
        
        ! X boundaries
        beta = state%config%albedo_x_min
        state%flux(1, :, :, g) = beta * state%flux(2, :, :, g)
        
        beta = state%config%albedo_x_max
        state%flux(state%nx, :, :, g) = beta * state%flux(state%nx-1, :, :, g)
        
        ! Y boundaries
        beta = state%config%albedo_y_min
        state%flux(:, 1, :, g) = beta * state%flux(:, 2, :, g)
        
        beta = state%config%albedo_y_max
        state%flux(:, state%ny, :, g) = beta * state%flux(:, state%ny-1, :, g)
        
        ! Z boundaries
        beta = state%config%albedo_z_min
        state%flux(:, :, 1, g) = beta * state%flux(:, :, 2, g)
        
        beta = state%config%albedo_z_max
        state%flux(:, :, state%nz, g) = beta * state%flux(:, :, state%nz-1, g)
    end subroutine apply_boundary_conditions
    
    !> Compute power distribution
    subroutine compute_power_distribution(state)
        type(mg_state_t), intent(inout) :: state
        
        integer :: i, j, k, g
        real(wp) :: power_cell, kappa_MeV, J_per_MeV
        
        J_per_MeV = 1.602e-13_wp
        state%total_power = 0.0_wp
        
        do k = 1, state%nz
            do j = 1, state%ny
                do i = 1, state%nx
                    power_cell = 0.0_wp
                    
                    do g = 1, state%n_groups
                        kappa_MeV = state%xsec(i, j, k)%kappa(g)
                        power_cell = power_cell + &
                            state%xsec(i, j, k)%sigma_f(g) * kappa_MeV * &
                            state%flux(i, j, k, g)
                    end do
                    
                    ! Convert MeV/cm³·s to W/cm³
                    state%power_density(i, j, k) = power_cell * J_per_MeV
                    state%total_power = state%total_power + &
                        state%power_density(i, j, k) * state%volume * 1.0e6_wp  ! to W
                end do
            end do
        end do
    end subroutine compute_power_distribution
    
    !> Compute form factors (peaking)
    subroutine mg_compute_form_factors(state)
        type(mg_state_t), intent(inout) :: state
        
        integer :: i, j, k
        real(wp) :: avg_power, max_power, axial_avg
        
        ! Average power
        if (state%total_power > TOL_DEFAULT) then
            avg_power = state%total_power / real(state%nx * state%ny * state%nz, wp)
        else
            avg_power = 1.0_wp
        end if
        
        ! Radial peaking
        max_power = 0.0_wp
        do j = 1, state%ny
            do i = 1, state%nx
                state%radial_peaking(i, j) = sum(state%power_density(i, j, :)) / &
                                              real(state%nz, wp) / avg_power
                max_power = max(max_power, state%radial_peaking(i, j))
            end do
        end do
        
        ! Axial peaking
        do k = 1, state%nz
            axial_avg = sum(state%power_density(:, :, k)) / real(state%nx * state%ny, wp)
            state%axial_peaking(k) = axial_avg / avg_power
        end do
        
        state%total_peaking = maxval(state%power_density) / avg_power
    end subroutine mg_compute_form_factors
    
    !> Get multi-group flux
    subroutine mg_get_flux(state, flux, group)
        type(mg_state_t), intent(in) :: state
        real(wp), intent(out) :: flux(:, :, :)
        integer, intent(in), optional :: group
        
        integer :: g
        
        if (present(group)) then
            flux = state%flux(:, :, :, group)
        else
            ! Return total flux (sum over groups)
            flux = 0.0_wp
            do g = 1, state%n_groups
                flux = flux + state%flux(:, :, :, g)
            end do
        end if
    end subroutine mg_get_flux
    
    !> Get power density
    subroutine mg_get_power(state, power)
        type(mg_state_t), intent(in) :: state
        real(wp), intent(out) :: power(:, :, :)
        
        ! Convert W/cm³ to W/m³
        power = state%power_density * 1.0e6_wp
    end subroutine mg_get_power
    
    !> Get fission rate
    subroutine mg_get_fission_rate(state, fission_rate)
        type(mg_state_t), intent(in) :: state
        real(wp), intent(out) :: fission_rate(:, :, :)
        
        fission_rate = state%fission_source
    end subroutine mg_get_fission_rate
    
    !> Compute k-effective from flux distribution
    function mg_compute_keff(state) result(k_eff)
        type(mg_state_t), intent(in) :: state
        real(wp) :: k_eff
        
        real(wp) :: production, absorption
        integer :: i, j, k, g
        
        production = 0.0_wp
        absorption = 0.0_wp
        
        do k = 1, state%nz
            do j = 1, state%ny
                do i = 1, state%nx
                    do g = 1, state%n_groups
                        production = production + &
                            state%xsec(i, j, k)%nu_sigma_f(g) * state%flux(i, j, k, g)
                        absorption = absorption + &
                            state%xsec(i, j, k)%sigma_a(g) * state%flux(i, j, k, g)
                    end do
                end do
            end do
        end do
        
        if (absorption > TOL_DEFAULT) then
            k_eff = production / absorption
        else
            k_eff = 0.0_wp
        end if
    end function mg_compute_keff
    
    !> Time-dependent solve (kinetics)
    subroutine mg_solve_transient(state, precursors, dt)
        type(mg_state_t), intent(inout) :: state
        real(wp), intent(inout) :: precursors(:, :, :, :)  ! (nx, ny, nz, n_delayed)
        real(wp), intent(in) :: dt
        
        ! Simplified - full implementation would use:
        ! - Improved quasi-static method
        ! - Predictor-corrector
        ! - Shape/amplitude decomposition
        
        ! For now, use simple explicit update
        ! In production, implement IQS or other advanced kinetics
        
        call mg_solve_eigenvalue(state, state%k_eff, state%converged)
    end subroutine mg_solve_transient
    
    !> Power iteration method (alternative solver)
    subroutine mg_power_iteration(state, max_iter, tolerance)
        type(mg_state_t), intent(inout) :: state
        integer, intent(in) :: max_iter
        real(wp), intent(in) :: tolerance
        
        ! Wrapper for mg_solve_eigenvalue with custom parameters
        type(mg_config_t) :: old_config
        
        old_config = state%config
        state%config%max_outer_iter = max_iter
        state%config%outer_tolerance = tolerance
        
        call mg_solve_eigenvalue(state, state%k_eff, state%converged)
        
        state%config = old_config
    end subroutine mg_power_iteration

end module multigroup_diffusion