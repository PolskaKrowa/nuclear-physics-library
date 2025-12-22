! physics/nuclear_fission.f90
!
! Nuclear fission energy release model for reactor simulation.
! Implements:
! - Point kinetics equations
! - Multi-group neutron diffusion
! - Fission power density calculation
! - Decay heat generation
! - Delayed neutron precursors
!
! Usage:
!   use nuclear_fission
!   type(fission_config_t) :: config
!   type(fission_state_t) :: state
!   
!   call fission_init(state, nx, ny, nz, dx, dy, dz, config)
!   call fission_set_cross_sections(state, sigma_f, sigma_a, nu, kappa)
!   call fission_step(state, dt)
!   call fission_get_power_density(state, power)
!   call fission_destroy(state)
!
module nuclear_fission
    use kinds, only: wp, i32
    use constants, only: TOL_DEFAULT, E_CHARGE
    use finite_difference, only: fd_laplacian_3d
    use rk4, only: rk4_step, rk4_config_t
    implicit none
    private
    
    ! Public interface
    public :: fission_config_t, fission_state_t, nuclear_data_t
    public :: fission_init, fission_destroy
    public :: fission_set_cross_sections, fission_set_reactivity
    public :: fission_step, fission_step_point_kinetics
    public :: fission_get_power_density, fission_get_total_power
    public :: fission_compute_decay_heat
    
    ! Delayed neutron group parameters (6-group model)
    integer, parameter :: N_DELAYED_GROUPS = 6
    
    real(wp), parameter :: BETA_DELAYED(N_DELAYED_GROUPS) = &
        [0.000215_wp, 0.001424_wp, 0.001274_wp, 0.002568_wp, 0.000748_wp, 0.000273_wp]
    
    real(wp), parameter :: LAMBDA_DELAYED(N_DELAYED_GROUPS) = &
        [0.0124_wp, 0.0305_wp, 0.111_wp, 0.301_wp, 1.14_wp, 3.01_wp]  ! 1/s
    
    ! Nuclear data per material
    type :: nuclear_data_t
        real(wp) :: sigma_f         ! Fission cross section [cm⁻¹]
        real(wp) :: sigma_a         ! Absorption cross section [cm⁻¹]
        real(wp) :: sigma_s         ! Scattering cross section [cm⁻¹]
        real(wp) :: nu              ! Neutrons per fission
        real(wp) :: kappa           ! Energy per fission [MeV]
        real(wp) :: D               ! Diffusion coefficient [cm]
        real(wp) :: chi             ! Fission spectrum
        logical :: is_fuel          ! Is this a fuel region?
    end type nuclear_data_t
    
    ! Configuration
    type :: fission_config_t
        logical :: use_point_kinetics = .false.
        logical :: include_decay_heat = .true.
        logical :: include_delayed_neutrons = .true.
        integer :: num_energy_groups = 1         ! 1 for one-group
        real(wp) :: prompt_neutron_lifetime = 1.0e-4_wp  ! seconds
        real(wp) :: beta_total = 0.0065_wp       ! Total delayed neutron fraction
        real(wp) :: power_normalisation = 1.0e6_wp  ! Watts (1 MW)
    end type fission_config_t
    
    ! Fission state
    type :: fission_state_t
        ! Grid dimensions
        integer :: nx, ny, nz
        real(wp) :: dx, dy, dz
        real(wp) :: volume
        
        ! Neutron flux [n/cm²·s]
        real(wp), allocatable :: flux(:, :, :)
        real(wp), allocatable :: flux_old(:, :, :)
        
        ! Delayed neutron precursor concentrations [1/cm³]
        real(wp), allocatable :: precursors(:, :, :, :)  ! (nx, ny, nz, groups)
        
        ! Fission power density [W/m³]
        real(wp), allocatable :: power_density(:, :, :)
        
        ! Decay heat power [W/m³]
        real(wp), allocatable :: decay_heat(:, :, :)
        
        ! Nuclear data (per cell)
        type(nuclear_data_t), allocatable :: nuclear_data(:, :, :)
        
        ! Point kinetics variables
        real(wp) :: total_power = 0.0_wp
        real(wp) :: reactivity = 0.0_wp
        real(wp) :: neutron_population = 1.0_wp
        
        ! Configuration
        type(fission_config_t) :: config
        
        ! Time tracking
        real(wp) :: time = 0.0_wp
        integer :: steps = 0
        
        ! Cumulative fission energy [J]
        real(wp) :: cumulative_energy = 0.0_wp
    end type fission_state_t
    
contains

    !> Initialise fission state
    subroutine fission_init(state, nx, ny, nz, dx, dy, dz, config)
        type(fission_state_t), intent(out) :: state
        integer, intent(in) :: nx, ny, nz
        real(wp), intent(in) :: dx, dy, dz
        type(fission_config_t), intent(in), optional :: config
        
        state%nx = nx
        state%ny = ny
        state%nz = nz
        state%dx = dx * 0.01_wp  ! Convert m to cm
        state%dy = dy * 0.01_wp
        state%dz = dz * 0.01_wp
        state%volume = dx * dy * dz  ! m³
        
        ! Allocate fields
        allocate(state%flux(nx, ny, nz))
        allocate(state%flux_old(nx, ny, nz))
        allocate(state%precursors(nx, ny, nz, N_DELAYED_GROUPS))
        allocate(state%power_density(nx, ny, nz))
        allocate(state%decay_heat(nx, ny, nz))
        allocate(state%nuclear_data(nx, ny, nz))
        
        ! Initialise fields
        state%flux = 1.0e10_wp           ! Initial flux guess
        state%flux_old = 1.0e10_wp
        state%precursors = 0.0_wp
        state%power_density = 0.0_wp
        state%decay_heat = 0.0_wp
        
        ! Default nuclear data (non-fuel)
        state%nuclear_data(:, :, :)%sigma_f = 0.0_wp
        state%nuclear_data(:, :, :)%sigma_a = 0.01_wp
        state%nuclear_data(:, :, :)%sigma_s = 0.1_wp
        state%nuclear_data(:, :, :)%nu = 0.0_wp
        state%nuclear_data(:, :, :)%kappa = 0.0_wp
        state%nuclear_data(:, :, :)%D = 1.0_wp
        state%nuclear_data(:, :, :)%chi = 1.0_wp
        state%nuclear_data(:, :, :)%is_fuel = .false.
        
        ! Configuration
        if (present(config)) then
            state%config = config
        end if
    end subroutine fission_init
    
    !> Destroy fission state
    subroutine fission_destroy(state)
        type(fission_state_t), intent(inout) :: state
        
        if (allocated(state%flux)) deallocate(state%flux)
        if (allocated(state%flux_old)) deallocate(state%flux_old)
        if (allocated(state%precursors)) deallocate(state%precursors)
        if (allocated(state%power_density)) deallocate(state%power_density)
        if (allocated(state%decay_heat)) deallocate(state%decay_heat)
        if (allocated(state%nuclear_data)) deallocate(state%nuclear_data)
    end subroutine fission_destroy
    
    !> Set nuclear cross sections for a region
    subroutine fission_set_cross_sections(state, sigma_f, sigma_a, nu, kappa, &
                                          i1, i2, j1, j2, k1, k2)
        type(fission_state_t), intent(inout) :: state
        real(wp), intent(in) :: sigma_f, sigma_a, nu, kappa
        integer, intent(in), optional :: i1, i2, j1, j2, k1, k2
        
        integer :: imin, imax, jmin, jmax, kmin, kmax
        real(wp) :: D
        
        imin = 1; imax = state%nx
        jmin = 1; jmax = state%ny
        kmin = 1; kmax = state%nz
        
        if (present(i1)) imin = i1
        if (present(i2)) imax = i2
        if (present(j1)) jmin = j1
        if (present(j2)) jmax = j2
        if (present(k1)) kmin = k1
        if (present(k2)) kmax = k2
        
        ! Diffusion coefficient: D = 1/(3·Σ_tr)
        D = 1.0_wp / (3.0_wp * (sigma_a + state%nuclear_data(imin, jmin, kmin)%sigma_s))
        
        state%nuclear_data(imin:imax, jmin:jmax, kmin:kmax)%sigma_f = sigma_f
        state%nuclear_data(imin:imax, jmin:jmax, kmin:kmax)%sigma_a = sigma_a
        state%nuclear_data(imin:imax, jmin:jmax, kmin:kmax)%nu = nu
        state%nuclear_data(imin:imax, jmin:jmax, kmin:kmax)%kappa = kappa
        state%nuclear_data(imin:imax, jmin:jmax, kmin:kmax)%D = D
        state%nuclear_data(imin:imax, jmin:jmax, kmin:kmax)%is_fuel = (sigma_f > 0.0_wp)
    end subroutine fission_set_cross_sections
    
    !> Set reactivity for point kinetics
    subroutine fission_set_reactivity(state, rho)
        type(fission_state_t), intent(inout) :: state
        real(wp), intent(in) :: rho
        
        state%reactivity = rho
    end subroutine fission_set_reactivity
    
    !> Main time step
    subroutine fission_step(state, dt)
        type(fission_state_t), intent(inout) :: state
        real(wp), intent(in) :: dt
        
        if (state%config%use_point_kinetics) then
            call fission_step_point_kinetics(state, dt)
        else
            call fission_step_diffusion(state, dt)
        end if
        
        ! Update power density
        call fission_compute_power_density(state)
        
        ! Update decay heat
        if (state%config%include_decay_heat) then
            call fission_compute_decay_heat(state, dt)
        end if
        
        state%time = state%time + dt
        state%steps = state%steps + 1
    end subroutine fission_step
    
    !> Point kinetics equations (0D approximation)
    subroutine fission_step_point_kinetics(state, dt)
        type(fission_state_t), intent(inout) :: state
        real(wp), intent(in) :: dt
        
        real(wp) :: n, c(N_DELAYED_GROUPS), dndt, dcdt(N_DELAYED_GROUPS)
        real(wp) :: beta_total, lambda_eff
        integer :: i
        
        n = state%neutron_population
        beta_total = sum(BETA_DELAYED)
        
        ! Get precursor concentrations (spatially averaged)
        do i = 1, N_DELAYED_GROUPS
            c(i) = sum(state%precursors(:, :, :, i)) / &
                  real(state%nx * state%ny * state%nz, wp)
        end do
        
        ! Point kinetics equations:
        ! dn/dt = [(ρ-β)/Λ]·n + Σᵢ λᵢ·cᵢ
        ! dcᵢ/dt = (βᵢ/Λ)·n - λᵢ·cᵢ
        
        ! Neutron population
        dndt = ((state%reactivity - beta_total) / state%config%prompt_neutron_lifetime) * n
        do i = 1, N_DELAYED_GROUPS
            dndt = dndt + LAMBDA_DELAYED(i) * c(i)
        end do
        
        ! Precursor concentrations
        do i = 1, N_DELAYED_GROUPS
            dcdt(i) = (BETA_DELAYED(i) / state%config%prompt_neutron_lifetime) * n - &
                     LAMBDA_DELAYED(i) * c(i)
        end do
        
        ! Explicit Euler update
        state%neutron_population = n + dt * dndt
        
        do i = 1, N_DELAYED_GROUPS
            c(i) = c(i) + dt * dcdt(i)
            state%precursors(:, :, :, i) = c(i)
        end do
        
        ! Update flux based on power
        state%flux = state%neutron_population * 1.0e10_wp
        
        ! Total power
        state%total_power = state%neutron_population * state%config%power_normalisation
    end subroutine fission_step_point_kinetics
    
    !> Spatial neutron diffusion (simplified one-group)
    subroutine fission_step_diffusion(state, dt)
        type(fission_state_t), intent(inout) :: state
        real(wp), intent(in) :: dt
        
        real(wp), allocatable :: laplacian(:, :, :)
        real(wp), allocatable :: source(:, :, :)
        integer :: i, j, k, g
        real(wp) :: D, sigma_a, nu_sigma_f, production, absorption, leakage
        real(wp) :: prompt_source, delayed_source
        
        allocate(laplacian(state%nx, state%ny, state%nz))
        allocate(source(state%nx, state%ny, state%nz))
        
        state%flux_old = state%flux
        
        ! Compute Laplacian with variable diffusion coefficient
        call compute_weighted_laplacian(state, laplacian)
        
        ! Build source term
        source = 0.0_wp
        do k = 1, state%nz
            do j = 1, state%ny
                do i = 1, state%nx
                    ! Prompt fission source
                    nu_sigma_f = state%nuclear_data(i, j, k)%nu * &
                                state%nuclear_data(i, j, k)%sigma_f
                    prompt_source = (1.0_wp - state%config%beta_total) * &
                                   nu_sigma_f * state%flux(i, j, k)
                    
                    ! Delayed neutron source
                    delayed_source = 0.0_wp
                    if (state%config%include_delayed_neutrons) then
                        do g = 1, N_DELAYED_GROUPS
                            delayed_source = delayed_source + &
                                LAMBDA_DELAYED(g) * state%precursors(i, j, k, g)
                        end do
                    end if
                    
                    source(i, j, k) = prompt_source + delayed_source
                end do
            end do
        end do
        
        ! Time-dependent diffusion equation:
        ! (1/v)·∂ϕ/∂t = D·∇²ϕ - Σₐ·ϕ + source
        ! Simplified: assume v·dt << 1, so explicit update
        
        do k = 1, state%nz
            do j = 1, state%ny
                do i = 1, state%nx
                    D = state%nuclear_data(i, j, k)%D
                    sigma_a = state%nuclear_data(i, j, k)%sigma_a
                    
                    leakage = D * laplacian(i, j, k)
                    absorption = -sigma_a * state%flux(i, j, k)
                    production = source(i, j, k)
                    
                    state%flux(i, j, k) = state%flux(i, j, k) + dt * &
                        (leakage + absorption + production) / state%config%prompt_neutron_lifetime
                    
                    ! Ensure non-negative
                    state%flux(i, j, k) = max(state%flux(i, j, k), 0.0_wp)
                end do
            end do
        end do
        
        ! Update precursors
        if (state%config%include_delayed_neutrons) then
            call update_precursors(state, dt)
        end if
        
        deallocate(laplacian, source)
    end subroutine fission_step_diffusion
    
    !> Update delayed neutron precursors
    subroutine update_precursors(state, dt)
        type(fission_state_t), intent(inout) :: state
        real(wp), intent(in) :: dt
        
        integer :: i, j, k, g
        real(wp) :: production, decay, nu_sigma_f
        
        do k = 1, state%nz
            do j = 1, state%ny
                do i = 1, state%nx
                    if (.not. state%nuclear_data(i, j, k)%is_fuel) cycle
                    
                    nu_sigma_f = state%nuclear_data(i, j, k)%nu * &
                                state%nuclear_data(i, j, k)%sigma_f
                    
                    do g = 1, N_DELAYED_GROUPS
                        production = (BETA_DELAYED(g) / state%config%beta_total) * &
                                    nu_sigma_f * state%flux(i, j, k)
                        decay = LAMBDA_DELAYED(g) * state%precursors(i, j, k, g)
                        
                        state%precursors(i, j, k, g) = state%precursors(i, j, k, g) + &
                            dt * (production - decay)
                        
                        state%precursors(i, j, k, g) = max(state%precursors(i, j, k, g), 0.0_wp)
                    end do
                end do
            end do
        end do
    end subroutine update_precursors
    
    !> Compute power density from neutron flux
    subroutine fission_compute_power_density(state)
        type(fission_state_t), intent(inout) :: state
        
        integer :: i, j, k
        real(wp) :: kappa_joules, energy_per_fission
        
        ! Convert MeV to Joules: 1 MeV = 1.602e-13 J
        
        do k = 1, state%nz
            do j = 1, state%ny
                do i = 1, state%nx
                    if (.not. state%nuclear_data(i, j, k)%is_fuel) then
                        state%power_density(i, j, k) = 0.0_wp
                        cycle
                    end if
                    
                    ! Energy per fission in Joules
                    kappa_joules = state%nuclear_data(i, j, k)%kappa * 1.602e-13_wp
                    
                    ! Power density [W/m³] = Σf·ϕ·κ
                    ! Convert from cm⁻¹·(n/cm²·s)·J to W/m³
                    state%power_density(i, j, k) = &
                        state%nuclear_data(i, j, k)%sigma_f * &
                        state%flux(i, j, k) * kappa_joules * 1.0e6_wp
                end do
            end do
        end do
        
        ! Update cumulative energy
        state%cumulative_energy = state%cumulative_energy + &
            sum(state%power_density) * state%volume * &
            (state%time / real(state%steps, wp))
    end subroutine fission_compute_power_density
    
    !> Compute decay heat using ANS standard
    subroutine fission_compute_decay_heat(state, dt)
        type(fission_state_t), intent(inout) :: state
        real(wp), intent(in) :: dt
        
        integer :: i, j, k
        real(wp) :: decay_fraction
        
        ! Simplified decay heat model: fraction of fission power
        ! ANS-5.1 standard: approximately 6-7% immediately after shutdown
        
        do k = 1, state%nz
            do j = 1, state%ny
                do i = 1, state%nx
                    ! Simple model: exponential decay from fission power
                    if (state%flux(i, j, k) > 0.0_wp) then
                        ! During operation: ~7% of fission power
                        decay_fraction = 0.07_wp
                    else
                        ! After shutdown: decays with time
                        decay_fraction = 0.07_wp * exp(-0.1_wp * state%time)
                    end if
                    
                    state%decay_heat(i, j, k) = decay_fraction * state%power_density(i, j, k)
                end do
            end do
        end do
    end subroutine fission_compute_decay_heat
    
    !> Get total fission power
    function fission_get_total_power(state) result(power)
        type(fission_state_t), intent(in) :: state
        real(wp) :: power
        
        ! Sum power density over all cells [W]
        power = sum(state%power_density) * state%volume
    end function fission_get_total_power
    
    !> Get power density array
    subroutine fission_get_power_density(state, power)
        type(fission_state_t), intent(in) :: state
        real(wp), intent(out) :: power(:, :, :)
        
        power = state%power_density + state%decay_heat
    end subroutine fission_get_power_density
    
    !> Compute weighted Laplacian for variable diffusion coefficient
    subroutine compute_weighted_laplacian(state, laplacian)
        type(fission_state_t), intent(in) :: state
        real(wp), intent(out) :: laplacian(:, :, :)
        
        integer :: i, j, k
        real(wp) :: D_xm, D_xp, D_ym, D_yp, D_zm, D_zp
        real(wp) :: inv_dx2, inv_dy2, inv_dz2
        
        inv_dx2 = 1.0_wp / (state%dx * state%dx)
        inv_dy2 = 1.0_wp / (state%dy * state%dy)
        inv_dz2 = 1.0_wp / (state%dz * state%dz)
        
        do k = 2, state%nz - 1
            do j = 2, state%ny - 1
                do i = 2, state%nx - 1
                    ! Harmonic mean of diffusion coefficients at interfaces
                    D_xm = 2.0_wp * state%nuclear_data(i - 1, j, k)%D * &
                           state%nuclear_data(i, j, k)%D / &
                           (state%nuclear_data(i - 1, j, k)%D + state%nuclear_data(i, j, k)%D)
                    D_xp = 2.0_wp * state%nuclear_data(i, j, k)%D * &
                           state%nuclear_data(i + 1, j, k)%D / &
                           (state%nuclear_data(i, j, k)%D + state%nuclear_data(i + 1, j, k)%D)
                    
                    D_ym = 2.0_wp * state%nuclear_data(i, j - 1, k)%D * &
                           state%nuclear_data(i, j, k)%D / &
                           (state%nuclear_data(i, j - 1, k)%D + state%nuclear_data(i, j, k)%D)
                    D_yp = 2.0_wp * state%nuclear_data(i, j, k)%D * &
                           state%nuclear_data(i, j + 1, k)%D / &
                           (state%nuclear_data(i, j, k)%D + state%nuclear_data(i, j + 1, k)%D)
                    
                    D_zm = 2.0_wp * state%nuclear_data(i, j, k - 1)%D * &
                           state%nuclear_data(i, j, k)%D / &
                           (state%nuclear_data(i, j, k - 1)%D + state%nuclear_data(i, j, k)%D)
                    D_zp = 2.0_wp * state%nuclear_data(i, j, k)%D * &
                           state%nuclear_data(i, j, k + 1)%D / &
                           (state%nuclear_data(i, j, k)%D + state%nuclear_data(i, j, k + 1)%D)
                    
                    laplacian(i, j, k) = &
                        (D_xp * state%flux(i + 1, j, k) - &
                         (D_xp + D_xm) * state%flux(i, j, k) + &
                         D_xm * state%flux(i - 1, j, k)) * inv_dx2 + &
                        (D_yp * state%flux(i, j + 1, k) - &
                         (D_yp + D_ym) * state%flux(i, j, k) + &
                         D_ym * state%flux(i, j - 1, k)) * inv_dy2 + &
                        (D_zp * state%flux(i, j, k + 1) - &
                         (D_zp + D_zm) * state%flux(i, j, k) + &
                         D_zm * state%flux(i, j, k - 1)) * inv_dz2
                end do
            end do
        end do
        
        ! Boundaries (zero flux or reflective)
        laplacian(1, :, :) = 0.0_wp
        laplacian(state%nx, :, :) = 0.0_wp
        laplacian(:, 1, :) = 0.0_wp
        laplacian(:, state%ny, :) = 0.0_wp
        laplacian(:, :, 1) = 0.0_wp
        laplacian(:, :, state%nz) = 0.0_wp
    end subroutine compute_weighted_laplacian

end module nuclear_fission