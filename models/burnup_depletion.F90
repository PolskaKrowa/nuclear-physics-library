! models/burnup_depletion.f90
!
! Fuel burnup and isotopic depletion tracking.
! Implements:
! - Xenon-135 and Samarium-149 dynamics
! - Major actinide depletion (U-235, U-238, Pu-239, Pu-241)
! - Fission product poisoning
! - Burnup accumulation
! - Decay chains
!
! Xenon-135 chain:
!   Te-135 → I-135 → Xe-135 → Cs-135 (stable)
!           γ_I↗    ↙ γ_Xe
!          Fission
!
! Samarium-149 chain:
!   Nd-149 → Pm-149 → Sm-149 (stable)
!          ↗ γ_Pm
!      Fission
!
! Usage:
!   use burnup_depletion
!   type(burnup_state_t) :: burnup
!   
!   call burnup_init(burnup, nx, ny, nz, dx, dy, dz)
!   call burnup_set_initial_composition(burnup, enrichment)
!   call burnup_step(burnup, flux, power, dt)
!   call burnup_get_xenon(burnup, Xe_concentration)
!   call burnup_get_samarium(burnup, Sm_concentration)
!   call burnup_destroy(burnup)
!
module burnup_depletion
    use kinds, only: wp, i32
    use constants, only: TOL_DEFAULT, N_AVOGADRO
    use compute_mode, only: use_gpu, gpu_min_workload
#ifdef COMPUTATION_MODE_GPU
    use burnup_depletion_gpu, only: burnup_gpu_state_t, burnup_gpu_init, &
        burnup_gpu_cleanup, burnup_gpu_step, burnup_gpu_copy_from_device
#endif
#ifdef COMPUTATION_MODE_HYBRID
    use burnup_depletion_gpu, only: burnup_gpu_state_t, burnup_gpu_init, &
        burnup_gpu_cleanup, burnup_gpu_step, burnup_gpu_copy_from_device
#endif
    implicit none
    private
    
    ! Public interface
    public :: burnup_state_t, burnup_config_t, isotope_data_t
    public :: burnup_init, burnup_destroy
    public :: burnup_set_initial_composition
    public :: burnup_step, burnup_step_predictor_corrector
    public :: burnup_get_xenon, burnup_get_samarium, burnup_get_burnup
    public :: burnup_compute_reactivity_effect
    
    ! Decay constants [s⁻¹]
    real(wp), parameter :: LAMBDA_I135 = log(2.0_wp) / (6.57_wp * 3600.0_wp)   ! 6.57 h
    real(wp), parameter :: LAMBDA_XE135 = log(2.0_wp) / (9.14_wp * 3600.0_wp)  ! 9.14 h
    real(wp), parameter :: LAMBDA_PM149 = log(2.0_wp) / (53.1_wp * 3600.0_wp)  ! 53.1 h
    
    ! Fission yields (thermal fission of U-235)
    real(wp), parameter :: GAMMA_I135 = 0.061_wp    ! Direct I-135 yield
    real(wp), parameter :: GAMMA_XE135 = 0.003_wp   ! Direct Xe-135 yield
    real(wp), parameter :: GAMMA_PM149 = 0.0113_wp  ! Pm-149 yield
    real(wp), parameter :: GAMMA_SM149 = 0.0_wp     ! Direct Sm-149 (negligible)
    
    ! Microscopic cross sections [barns]
    real(wp), parameter :: SIGMA_A_XE135 = 2.65e6_wp
    real(wp), parameter :: SIGMA_A_SM149 = 4.08e4_wp
    real(wp), parameter :: SIGMA_A_U235 = 680.9_wp
    real(wp), parameter :: SIGMA_F_U235 = 585.1_wp
    real(wp), parameter :: SIGMA_A_U238 = 2.68_wp
    real(wp), parameter :: SIGMA_A_PU239 = 1017.7_wp
    real(wp), parameter :: SIGMA_F_PU239 = 748.1_wp
    real(wp), parameter :: SIGMA_A_PU241 = 1377.0_wp
    real(wp), parameter :: SIGMA_F_PU241 = 1012.0_wp
    
    ! Isotope data
    type :: isotope_data_t
        character(len=16) :: name
        real(wp) :: N                ! Number density [atoms/barn-cm]
        real(wp) :: sigma_a          ! Absorption xs [barns]
        real(wp) :: sigma_f          ! Fission xs [barns]
        real(wp) :: lambda           ! Decay constant [s⁻¹]
        real(wp) :: yield_fission    ! Fission product yield
    end type isotope_data_t
    
    ! Configuration
    type :: burnup_config_t
        logical :: track_xenon = .true.
        logical :: track_samarium = .true.
        logical :: track_actinides = .true.
        logical :: use_predictor_corrector = .true.
        real(wp) :: thermal_flux_fraction = 0.8_wp  ! For thermal XS
        integer :: depletion_substeps = 1
    end type burnup_config_t
    
    ! Burnup and depletion state
    type :: burnup_state_t
        ! Grid
        integer :: nx, ny, nz
        real(wp) :: dx, dy, dz
        real(wp) :: volume
        
        ! Burnup [MWd/kgU]
        real(wp), allocatable :: burnup(:, :, :)
        
        ! Poison concentrations [atoms/barn-cm]
        real(wp), allocatable :: I135(:, :, :)
        real(wp), allocatable :: Xe135(:, :, :)
        real(wp), allocatable :: Pm149(:, :, :)
        real(wp), allocatable :: Sm149(:, :, :)
        
        ! Major actinides [atoms/barn-cm]
        real(wp), allocatable :: U235(:, :, :)
        real(wp), allocatable :: U238(:, :, :)
        real(wp), allocatable :: Pu239(:, :, :)
        real(wp), allocatable :: Pu241(:, :, :)
        
        ! Fission product lumped poison
        real(wp), allocatable :: FP_poison(:, :, :)
        
        ! Integrated quantities
        real(wp), allocatable :: cumulative_fissions(:, :, :)
        real(wp), allocatable :: time_at_power(:, :, :)
        
        ! Configuration
        type(burnup_config_t) :: config
        
        ! Time tracking
        real(wp) :: time = 0.0_wp
        integer :: steps = 0
    end type burnup_state_t
    
contains

    !> Initialize burnup state
    subroutine burnup_init(state, nx, ny, nz, dx, dy, dz, config)
        type(burnup_state_t), intent(out) :: state
        integer, intent(in) :: nx, ny, nz
        real(wp), intent(in) :: dx, dy, dz
        type(burnup_config_t), intent(in), optional :: config
        
        state%nx = nx
        state%ny = ny
        state%nz = nz
        state%dx = dx
        state%dy = dy
        state%dz = dz
        state%volume = dx * dy * dz
        
        ! Allocate arrays
        allocate(state%burnup(nx, ny, nz))
        allocate(state%I135(nx, ny, nz))
        allocate(state%Xe135(nx, ny, nz))
        allocate(state%Pm149(nx, ny, nz))
        allocate(state%Sm149(nx, ny, nz))
        allocate(state%U235(nx, ny, nz))
        allocate(state%U238(nx, ny, nz))
        allocate(state%Pu239(nx, ny, nz))
        allocate(state%Pu241(nx, ny, nz))
        allocate(state%FP_poison(nx, ny, nz))
        allocate(state%cumulative_fissions(nx, ny, nz))
        allocate(state%time_at_power(nx, ny, nz))
        
        ! Initialize to zero
        state%burnup = 0.0_wp
        state%I135 = 0.0_wp
        state%Xe135 = 0.0_wp
        state%Pm149 = 0.0_wp
        state%Sm149 = 0.0_wp
        state%U235 = 0.0_wp
        state%U238 = 0.0_wp
        state%Pu239 = 0.0_wp
        state%Pu241 = 0.0_wp
        state%FP_poison = 0.0_wp
        state%cumulative_fissions = 0.0_wp
        state%time_at_power = 0.0_wp
        
        ! Configuration
        if (present(config)) then
            state%config = config
        end if
    end subroutine burnup_init
    
    !> Destroy burnup state
    subroutine burnup_destroy(state)
        type(burnup_state_t), intent(inout) :: state
        
        if (allocated(state%burnup)) deallocate(state%burnup)
        if (allocated(state%I135)) deallocate(state%I135)
        if (allocated(state%Xe135)) deallocate(state%Xe135)
        if (allocated(state%Pm149)) deallocate(state%Pm149)
        if (allocated(state%Sm149)) deallocate(state%Sm149)
        if (allocated(state%U235)) deallocate(state%U235)
        if (allocated(state%U238)) deallocate(state%U238)
        if (allocated(state%Pu239)) deallocate(state%Pu239)
        if (allocated(state%Pu241)) deallocate(state%Pu241)
        if (allocated(state%FP_poison)) deallocate(state%FP_poison)
        if (allocated(state%cumulative_fissions)) deallocate(state%cumulative_fissions)
        if (allocated(state%time_at_power)) deallocate(state%time_at_power)
    end subroutine burnup_destroy
    
    !> Set initial fuel composition
    subroutine burnup_set_initial_composition(state, enrichment, fuel_density)
        type(burnup_state_t), intent(inout) :: state
        real(wp), intent(in), optional :: enrichment
        real(wp), intent(in), optional :: fuel_density  ! g/cm³
        
        real(wp) :: enrich, rho_fuel, M_UO2, f_U235, f_U238
        real(wp) :: N_tot, N_U
        
        ! Default: 3.5% enrichment, UO2 at 10.97 g/cm³
        enrich = 0.035_wp
        if (present(enrichment)) enrich = enrichment
        
        rho_fuel = 10.97_wp
        if (present(fuel_density)) rho_fuel = fuel_density
        
        ! Molecular weight of UO2: ~270 g/mol
        M_UO2 = 270.0_wp
        
        ! Number density of UO2 molecules [molecules/cm³]
        N_tot = rho_fuel * N_AVOGADRO / M_UO2
        
        ! One uranium atom per UO2 molecule
        N_U = N_tot
        
        ! Convert to atoms/barn-cm (1 barn = 1e-24 cm²)
        N_U = N_U * 1.0e-24_wp
        
        ! Initial composition
        f_U235 = enrich
        f_U238 = 1.0_wp - enrich
        
        state%U235 = f_U235 * N_U
        state%U238 = f_U238 * N_U
        state%Pu239 = 0.0_wp
        state%Pu241 = 0.0_wp
        
        ! Fresh fuel has no poisons
        state%I135 = 0.0_wp
        state%Xe135 = 0.0_wp
        state%Pm149 = 0.0_wp
        state%Sm149 = 0.0_wp
        state%FP_poison = 0.0_wp
    end subroutine burnup_set_initial_composition
    
    !> Main depletion time step
    subroutine burnup_step(state, flux, power, dt)
        type(burnup_state_t), intent(inout) :: state
        real(wp), intent(in) :: flux(:, :, :)       ! Thermal flux [n/cm²·s]
        real(wp), intent(in) :: power(:, :, :)      ! Power density [W/cm³]
        real(wp), intent(in) :: dt                   ! Time step [s]
        
        integer :: i, j, k
        real(wp) :: phi_th, P
        
#if defined(COMPUTATION_MODE_GPU) || defined(COMPUTATION_MODE_HYBRID)
        ! GPU dispatch: offload per-cell depletion to GPU when large enough.
        if (use_gpu()) then
            if (state%nx * state%ny * state%nz >= gpu_min_workload()) then
                call burnup_step_gpu(state, flux, power, dt)
                return
            end if
        end if
#endif
        
        ! Note: deplete_cell handles its own substepping internally
        ! (it divides dt by n_substeps). The previous code here ALSO
        ! divided dt by n_substeps and then passed the already-divided
        ! value to deplete_cell, which divided AGAIN — so each cell
        ! only advanced by dt/n_substeps² instead of dt. We now pass
        ! the full dt to deplete_cell and let it do the substepping.
        !
        ! OpenMP: each cell's isotope concentrations are independent,
        ! so deplete_cell calls on different (i,j,k) are race-free.
        !$OMP PARALLEL DO PRIVATE(i, j, k, phi_th, P) SCHEDULE(STATIC)
        do k = 1, state%nz
            do j = 1, state%ny
                do i = 1, state%nx
                    phi_th = flux(i, j, k) * state%config%thermal_flux_fraction
                    P = power(i, j, k)
                    
                    ! Deplete each cell (full dt; deplete_cell substeps internally)
                    call deplete_cell(state, i, j, k, phi_th, P, dt, &
                                     state%config%depletion_substeps)
                end do
            end do
        end do
        !$OMP END PARALLEL DO
        
        state%time = state%time + dt
        state%steps = state%steps + 1
    end subroutine burnup_step
    
    !> Predictor-corrector depletion step.
    !! Saves the initial state BEFORE the predictor, then averages the
    !! initial and predicted concentrations. (Previously the copy was
    !! taken AFTER burnup_step had already advanced the state, so the
    !! "corrector" averaged the state with itself — a no-op.)
    !!
    !! Performance note: the previous implementation deep-copied the entire
    !! burnup_state_t TWICE (state_init + state_pred). Each copy duplicates
    !! all 12+ allocatable isotope arrays — O(Nx·Ny·Nz·8) bytes each. We
    !! now keep only the single state_init copy and average in-place into
    !! state (which already holds the predicted values after burnup_step),
    !! eliminating one full deep-copy.
    subroutine burnup_step_predictor_corrector(state, flux, power, dt)
        type(burnup_state_t), intent(inout) :: state
        real(wp), intent(in) :: flux(:, :, :)
        real(wp), intent(in) :: power(:, :, :)
        real(wp), intent(in) :: dt
        
        type(burnup_state_t) :: state_init
        integer :: i, j, k
        
        ! Save initial state BEFORE the predictor step.
        state_init = state
        
        ! Predictor: full step with initial conditions. After this call,
        ! 'state' holds the predicted concentrations.
        call burnup_step(state, flux, power, dt)
        
        ! Corrector: average initial and predicted concentrations in-place.
        ! OpenMP: each cell's averaging is independent.
        !$OMP PARALLEL DO PRIVATE(i, j, k) SCHEDULE(STATIC)
        do k = 1, state%nz
            do j = 1, state%ny
                do i = 1, state%nx
                    state%Xe135(i, j, k) = 0.5_wp * (state_init%Xe135(i, j, k) + &
                                                     state%Xe135(i, j, k))
                    state%I135(i, j, k)  = 0.5_wp * (state_init%I135(i, j, k)  + &
                                                     state%I135(i, j, k))
                    state%Sm149(i, j, k) = 0.5_wp * (state_init%Sm149(i, j, k) + &
                                                     state%Sm149(i, j, k))
                    state%Pm149(i, j, k) = 0.5_wp * (state_init%Pm149(i, j, k) + &
                                                     state%Pm149(i, j, k))
                    state%U235(i, j, k)  = 0.5_wp * (state_init%U235(i, j, k)  + &
                                                     state%U235(i, j, k))
                    state%U238(i, j, k)  = 0.5_wp * (state_init%U238(i, j, k)  + &
                                                     state%U238(i, j, k))
                    state%Pu239(i, j, k) = 0.5_wp * (state_init%Pu239(i, j, k) + &
                                                     state%Pu239(i, j, k))
                    state%Pu241(i, j, k) = 0.5_wp * (state_init%Pu241(i, j, k) + &
                                                     state%Pu241(i, j, k))
                end do
            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine burnup_step_predictor_corrector
    
    !> Deplete single cell
    subroutine deplete_cell(state, i, j, k, phi_thermal, power, dt, n_substeps)
        type(burnup_state_t), intent(inout) :: state
        integer, intent(in) :: i, j, k
        real(wp), intent(in) :: phi_thermal, power, dt
        integer, intent(in) :: n_substeps
        
        integer :: step
        real(wp) :: dt_sub
        
        dt_sub = dt / real(n_substeps, wp)
        
        do step = 1, n_substeps
            ! Update poisons
            if (state%config%track_xenon) then
                call update_xenon(state, i, j, k, phi_thermal, dt_sub)
            end if
            
            if (state%config%track_samarium) then
                call update_samarium(state, i, j, k, phi_thermal, dt_sub)
            end if
            
            ! Update actinides
            if (state%config%track_actinides) then
                call update_actinides(state, i, j, k, phi_thermal, dt_sub)
            end if
            
            ! Update burnup
            call update_burnup(state, i, j, k, power, dt_sub)
        end do
    end subroutine deplete_cell
    
    !> Update xenon-135 concentration
    subroutine update_xenon(state, i, j, k, phi, dt)
        type(burnup_state_t), intent(inout) :: state
        integer, intent(in) :: i, j, k
        real(wp), intent(in) :: phi, dt
        
        real(wp) :: I135, Xe135, Sigma_f, rate_fission
        real(wp) :: dI_dt, dXe_dt
        
        I135 = state%I135(i, j, k)
        Xe135 = state%Xe135(i, j, k)
        
        ! Fission rate [fissions/cm³·s]
        Sigma_f = SIGMA_F_U235 * state%U235(i, j, k) + &
                  SIGMA_F_PU239 * state%Pu239(i, j, k) + &
                  SIGMA_F_PU241 * state%Pu241(i, j, k)
        rate_fission = Sigma_f * phi * 1.0e-24_wp  ! Convert barns to cm²
        
        ! Iodine-135 balance:
        ! dI/dt = γ_I·Σ_f·φ - λ_I·I
        dI_dt = GAMMA_I135 * rate_fission - LAMBDA_I135 * I135
        
        ! Xenon-135 balance:
        ! dXe/dt = γ_Xe·Σ_f·φ + λ_I·I - λ_Xe·Xe - σ_a,Xe·φ·Xe
        dXe_dt = GAMMA_XE135 * rate_fission + LAMBDA_I135 * I135 - &
                 LAMBDA_XE135 * Xe135 - &
                 SIGMA_A_XE135 * 1.0e-24_wp * phi * Xe135
        
        ! Explicit Euler update
        state%I135(i, j, k) = I135 + dt * dI_dt
        state%Xe135(i, j, k) = Xe135 + dt * dXe_dt
        
        ! Ensure non-negative
        state%I135(i, j, k) = max(state%I135(i, j, k), 0.0_wp)
        state%Xe135(i, j, k) = max(state%Xe135(i, j, k), 0.0_wp)
    end subroutine update_xenon
    
    !> Update samarium-149 concentration
    subroutine update_samarium(state, i, j, k, phi, dt)
        type(burnup_state_t), intent(inout) :: state
        integer, intent(in) :: i, j, k
        real(wp), intent(in) :: phi, dt
        
        real(wp) :: Pm149, Sm149, Sigma_f, rate_fission
        real(wp) :: dPm_dt, dSm_dt
        
        Pm149 = state%Pm149(i, j, k)
        Sm149 = state%Sm149(i, j, k)
        
        ! Fission rate
        Sigma_f = SIGMA_F_U235 * state%U235(i, j, k) + &
                  SIGMA_F_PU239 * state%Pu239(i, j, k)
        rate_fission = Sigma_f * phi * 1.0e-24_wp
        
        ! Promethium-149 balance:
        ! dPm/dt = γ_Pm·Σ_f·φ - λ_Pm·Pm
        dPm_dt = GAMMA_PM149 * rate_fission - LAMBDA_PM149 * Pm149
        
        ! Samarium-149 balance:
        ! dSm/dt = λ_Pm·Pm - σ_a,Sm·φ·Sm
        ! (Sm-149 is stable, but burns out)
        dSm_dt = LAMBDA_PM149 * Pm149 - SIGMA_A_SM149 * 1.0e-24_wp * phi * Sm149
        
        ! Update
        state%Pm149(i, j, k) = Pm149 + dt * dPm_dt
        state%Sm149(i, j, k) = Sm149 + dt * dSm_dt
        
        ! Ensure non-negative
        state%Pm149(i, j, k) = max(state%Pm149(i, j, k), 0.0_wp)
        state%Sm149(i, j, k) = max(state%Sm149(i, j, k), 0.0_wp)
    end subroutine update_samarium
    
    !> Update actinide concentrations
    subroutine update_actinides(state, i, j, k, phi, dt)
        type(burnup_state_t), intent(inout) :: state
        integer, intent(in) :: i, j, k
        real(wp), intent(in) :: phi, dt
        
        real(wp) :: U235, U238, Pu239, Pu241
        real(wp) :: dU235_dt, dU238_dt, dPu239_dt, dPu241_dt
        real(wp) :: phi_barn
        
        U235 = state%U235(i, j, k)
        U238 = state%U238(i, j, k)
        Pu239 = state%Pu239(i, j, k)
        Pu241 = state%Pu241(i, j, k)
        
        phi_barn = phi * 1.0e-24_wp  ! Convert flux: n/cm²·s → n/barn·s
        
        ! U-235 depletion:
        ! dU235/dt = -σ_a,U235·φ·U235
        dU235_dt = -SIGMA_A_U235 * phi_barn * U235
        
        ! U-238 depletion and Pu-239 production:
        ! dU238/dt = -σ_a,U238·φ·U238
        dU238_dt = -SIGMA_A_U238 * phi_barn * U238
        
        ! Pu-239 production and depletion:
        ! dPu239/dt = 0.9·σ_a,U238·φ·U238 - σ_a,Pu239·φ·Pu239
        ! (0.9 accounts for U-239 → Np-239 → Pu-239 chain efficiency)
        dPu239_dt = 0.9_wp * SIGMA_A_U238 * phi_barn * U238 - &
                   SIGMA_A_PU239 * phi_barn * Pu239
        
        ! Pu-241 (simplified):
        ! dPu241/dt = 0.5·σ_a,Pu239·φ·Pu239 - σ_a,Pu241·φ·Pu241
        dPu241_dt = 0.5_wp * (SIGMA_A_PU239 - SIGMA_F_PU239) * phi_barn * Pu239 - &
                   SIGMA_A_PU241 * phi_barn * Pu241
        
        ! Update
        state%U235(i, j, k) = U235 + dt * dU235_dt
        state%U238(i, j, k) = U238 + dt * dU238_dt
        state%Pu239(i, j, k) = Pu239 + dt * dPu239_dt
        state%Pu241(i, j, k) = Pu241 + dt * dPu241_dt
        
        ! Ensure non-negative
        state%U235(i, j, k) = max(state%U235(i, j, k), 0.0_wp)
        state%U238(i, j, k) = max(state%U238(i, j, k), 0.0_wp)
        state%Pu239(i, j, k) = max(state%Pu239(i, j, k), 0.0_wp)
        state%Pu241(i, j, k) = max(state%Pu241(i, j, k), 0.0_wp)
        
        ! Track fissions
        state%cumulative_fissions(i, j, k) = state%cumulative_fissions(i, j, k) + &
            dt * phi_barn * (SIGMA_F_U235 * U235 + SIGMA_F_PU239 * Pu239 + SIGMA_F_PU241 * Pu241)
    end subroutine update_actinides
    
    !> Update burnup
    subroutine update_burnup(state, i, j, k, power, dt)
        type(burnup_state_t), intent(inout) :: state
        integer, intent(in) :: i, j, k
        real(wp), intent(in) :: power, dt  ! power in W/cm³
        
        real(wp) :: energy_MWd, mass_U_kg
        
        if (power < TOL_DEFAULT) return
        
        ! Energy released [MWd] = Power [W/cm³] × Volume [cm³] × Time [s] / (1e6 × 86400)
        energy_MWd = power * 1.0e6_wp * dt / (1.0e6_wp * 86400.0_wp)  ! per cm³
        
        ! Mass of uranium [kg/cm³]
        ! N_U [atoms/barn-cm] × A_U [g/mol] / N_A [1/mol] × 1e-24 [barn→cm²] × 1e-3 [g→kg]
        mass_U_kg = (state%U235(i, j, k) + state%U238(i, j, k)) * 238.0_wp / N_AVOGADRO * 1.0e-27_wp
        
        if (mass_U_kg > TOL_DEFAULT) then
            ! Burnup [MWd/kgU]
            state%burnup(i, j, k) = state%burnup(i, j, k) + energy_MWd / mass_U_kg
        end if
        
        ! Track time at power
        if (power > 1.0_wp) then
            state%time_at_power(i, j, k) = state%time_at_power(i, j, k) + dt
        end if
    end subroutine update_burnup
    
    !> Get xenon-135 concentration
    subroutine burnup_get_xenon(state, Xe_conc)
        type(burnup_state_t), intent(in) :: state
        real(wp), intent(out) :: Xe_conc(:, :, :)
        
        Xe_conc = state%Xe135
    end subroutine burnup_get_xenon
    
    !> Get samarium-149 concentration
    subroutine burnup_get_samarium(state, Sm_conc)
        type(burnup_state_t), intent(in) :: state
        real(wp), intent(out) :: Sm_conc(:, :, :)
        
        Sm_conc = state%Sm149
    end subroutine burnup_get_samarium
    
    !> Get burnup
    subroutine burnup_get_burnup(state, burnup)
        type(burnup_state_t), intent(in) :: state
        real(wp), intent(out) :: burnup(:, :, :)
        
        burnup = state%burnup
    end subroutine burnup_get_burnup
    
    !> Compute reactivity effect of poisons
    function burnup_compute_reactivity_effect(state, phi_thermal) result(rho)
        type(burnup_state_t), intent(in) :: state
        real(wp), intent(in) :: phi_thermal  ! Average thermal flux
        real(wp) :: rho  ! Reactivity [pcm]
        
        real(wp) :: Xe_avg, Sm_avg, worth_Xe, worth_Sm
        
        ! Average concentrations
        Xe_avg = sum(state%Xe135) / real(state%nx * state%ny * state%nz, wp)
        Sm_avg = sum(state%Sm149) / real(state%nx * state%ny * state%nz, wp)
        
        ! Reactivity worth (simplified)
        ! ρ = -Σ_a,poison / Σ_a,fuel
        worth_Xe = -Xe_avg * SIGMA_A_XE135 / (SIGMA_A_U235 * 0.01_wp) * 1.0e5_wp  ! pcm
        worth_Sm = -Sm_avg * SIGMA_A_SM149 / (SIGMA_A_U235 * 0.01_wp) * 1.0e5_wp
        
        rho = worth_Xe + worth_Sm
    end function burnup_compute_reactivity_effect

#if defined(COMPUTATION_MODE_GPU) || defined(COMPUTATION_MODE_HYBRID)
    ! -----------------------------------------------------------------
    ! GPU dispatch wrapper for the burnup depletion step.
    ! Creates a temporary GPU state, runs the per-cell depletion kernel,
    ! and copies the updated isotope concentrations back.
    ! -----------------------------------------------------------------
    subroutine burnup_step_gpu(state, flux, power, dt)
        type(burnup_state_t), intent(inout) :: state
        real(wp), intent(in) :: flux(state%nx, state%ny, state%nz)
        real(wp), intent(in) :: power(state%nx, state%ny, state%nz)
        real(wp), intent(in) :: dt
        type(burnup_gpu_state_t) :: g
        real(wp), allocatable :: thermal_flux(:,:,:)

        ! Scale flux by thermal fraction (matching the CPU path)
        allocate(thermal_flux(state%nx, state%ny, state%nz))
        thermal_flux = flux * state%config%thermal_flux_fraction

        ! Initialise and run GPU kernel
        call burnup_gpu_init(g, state%nx, state%ny, state%nz)
        call burnup_gpu_step(g, thermal_flux, power, dt, state%config%depletion_substeps)

        ! Copy updated concentrations back
        call burnup_gpu_copy_from_device(g, state%Xe135, state%I135, state%Sm149, &
            state%Pm149, state%U235, state%U238, state%Pu239, state%Pu241, &
            state%burnup)

        call burnup_gpu_cleanup(g)
        deallocate(thermal_flux)

        ! Update state metadata
        state%time = state%time + dt
        state%steps = state%steps + 1
    end subroutine burnup_step_gpu
#endif

end module burnup_depletion