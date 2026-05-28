! bwr_core_simulator.f90
!
! Non-OpenGL port of bwr_core_simulator_opengl.f90.
! All physics modules are identical to the original; the OpenGL/GLUT
! rendering layer has been replaced by a text-based interactive loop.
!
! Controls (type at the "sim>" prompt):
!   <Enter>      - advance STEPS_PER_PROMPT time steps and print summary
!   p            - pause / resume
!   1-5          - display focus (Power/Temp/Void/Burnup/T_sat ratio)
!   +/-          - increase / decrease power setpoint  5 %
!   r/R          - insert +20 / -20 pcm reactivity
!   c/C          - increase / decrease coolant flow  5 %
!   i/I          - insert / withdraw control rods  0.5 %
!   t/T          - open / close turbine valve  2 %
!   b/B          - open / close bypass valve   2 %
!   s            - print detailed steady-state summary
!   h/?          - print this help
!   q / quit     - quit

program bwr_core_simulator

    use kinds, only: wp

    ! Physics modules (unchanged from original)
    use multigroup_diffusion
    use cross_sections
    use burnup_depletion
    use heat_transfer
    use two_phase_flow
    use rng, only: rng_seed

    implicit none

    ! ── Simulation state ────────────────────────────────────────
    type :: simulation_t
        ! Grid
        integer  :: nx, ny, nz
        real(wp) :: dx, dy, dz
        real(wp) :: core_height, core_diameter
        real(wp) :: rod_bank_position  ! 0.0 (All Out) to 1.0 (All In) – BWR bottom entry

        ! Physics modules
        type(mg_state_t)         :: neutronics
        type(heat_state_t)       :: heat
        type(two_phase_state_t)  :: thermalhydraulics
        type(burnup_state_t)     :: burnup
        type(xsec_library_t)     :: xslib

        ! Operating conditions
        real(wp) :: power_rated, power_current
        real(wp) :: pressure_operating    ! dynamic [Pa]
        real(wp) :: mass_flux_core
        real(wp) :: inlet_temperature

        ! Core parameters
        real(wp) :: k_eff, reactivity_pcm
        real(wp) :: max_fuel_temp, max_clad_temp
        real(wp) :: avg_void_fraction, max_void_fraction
        real(wp) :: min_chfr, avg_burnup

        ! Feedback coefficients
        real(wp) :: alpha_doppler, alpha_void

        ! Time
        real(wp) :: time
        integer  :: n_steps

        ! Initialisation flag
        logical :: initialized

        ! Saturation and level
        real(wp) :: sat_temperature    ! Saturation temperature at Prx [K]
        real(wp) :: Lbase              ! Effective water inventory [m equivalent]
        real(wp) :: Lrx                ! Water level above core floor [m]

        ! 3-point smoothing buffers for coolant temperature
        real(wp) :: avg_coolant_temp
        real(wp) :: coolant_T_prev1
        real(wp) :: coolant_T_prev2

        ! 3-point smoothing for pressure
        real(wp) :: P_prev1
        real(wp) :: P_prev2

        ! Flash steam excess above saturation
        real(wp) :: T_excess

        ! Reactor period
        real(wp) :: neutrons_prev      ! Previous power level [W] for period calc
        real(wp) :: reactor_period     ! Reactor period [s]; +∞ = stable

        ! APRM instrumentation
        real(wp) :: aprm               ! Average Power Range Monitor [% rated]

        ! Void reactivity factor
        real(wp) :: rf_void

        ! Turbine
        real(wp) :: turbine_valve      ! Steam admission valve [0-100 %]
        real(wp) :: bypass_valve       ! Bypass valve [0-100 %]
        real(wp) :: steam_flow         ! Normalised total steam flow
        real(wp) :: turbine_speed      ! Rotational speed [rpm]
        real(wp) :: turbine_power      ! Net electrical output [W]
        logical  :: turbine_tripped    ! Turbine trip / isolation flag

        ! Feedwater temperature
        real(wp) :: feedwater_temp     ! [K]
    end type simulation_t

    ! ── Global state ────────────────────────────────────────────
    type(simulation_t), target, save :: g_sim
    real(wp), save :: g_dt            = 0.01_wp
    logical,  save :: g_converged     = .false.
    logical,  save :: g_paused        = .false.   ! replaces g_viz%paused
    integer,  save :: g_display_mode  = 1         ! 1=power,2=temp,3=void,4=burnup,5=T/T_sat

    ! Steps to advance before prompting the user
    integer, parameter :: STEPS_PER_PROMPT = 100

    ! ── Local variables ──────────────────────────────────────────
    integer            :: step, ios
    character(len=80)  :: cmd
    logical            :: running

    ! ── Banner ──────────────────────────────────────────────────
    print *, "=================================================="
    print *, "  NPL Example BWR Simulation Program"
    print *, "  (text-mode; no OpenGL)"
    print *, "=================================================="
    print *, ""

    call rng_seed(123456789_8)

    print *, "Setting up simulation..."
    call setup_bwr_simulation(g_sim)
    g_sim%initialized = .true.

    print *, "Computing initial steady-state..."
    call solve_steady_state(g_sim)
    call print_steady_state_summary(g_sim)

    call print_help()

    ! ── Main simulation loop ─────────────────────────────────────
    running = .true.

    do while (running)

        if (.not. g_paused) then
            ! Advance STEPS_PER_PROMPT physics steps
            do step = 1, STEPS_PER_PROMPT
                call coupled_time_step(g_sim, g_dt)
                call apply_automatic_controls(g_sim, g_dt)
            end do
            call print_transient_summary(g_sim, g_sim%n_steps, g_sim%time)
        end if

        ! Command prompt
        write(*, '(A)', advance='no') "sim> "
        read(*, '(A)', iostat=ios) cmd
        if (ios /= 0) exit   ! EOF (e.g. piped input exhausted)

        cmd = adjustl(cmd)
        call handle_command(cmd, running)

    end do

    print *, ""
    print *, "Shutting down..."
    call cleanup_simulation(g_sim)

contains

    ! ================================================================
    !  COMMAND HANDLER  (replaces keyboard_callback)
    ! ================================================================

    subroutine handle_command(cmd, running)
        character(len=*), intent(in)    :: cmd
        logical,          intent(inout) :: running

        character :: ch

        if (len_trim(cmd) == 0) return   ! bare Enter – just advance

        ch = cmd(1:1)

        select case (ch)

        case ('p', 'P')
            g_paused = .not. g_paused
            print *, merge("PAUSED  ", "RESUMED ", g_paused)

        case ('1'); g_display_mode = 1; print *, "Focus: Power distribution"
        case ('2'); g_display_mode = 2; print *, "Focus: Temperature"
        case ('3'); g_display_mode = 3; print *, "Focus: Void fraction"
        case ('4'); g_display_mode = 4; print *, "Focus: Burnup"
        case ('5'); g_display_mode = 5; print *, "Focus: T / T_sat ratio"

        case ('+', '=')
            g_sim%power_rated = g_sim%power_rated * 1.05_wp
            print '(A,F8.1,A)', " Power setpoint: ", g_sim%power_rated/1.0e6_wp, " MW"

        case ('-')
            g_sim%power_rated = g_sim%power_rated * 0.95_wp
            print '(A,F8.1,A)', " Power setpoint: ", g_sim%power_rated/1.0e6_wp, " MW"

        case ('r')
            call apply_reactivity_insertion(g_sim, 20.0_wp)
            print *, " +20 pcm reactivity inserted"

        case ('R')
            call apply_reactivity_insertion(g_sim, -20.0_wp)
            print *, " -20 pcm reactivity inserted"

        case ('c')
            g_sim%mass_flux_core = g_sim%mass_flux_core * 1.05_wp
            print '(A,F7.1,A)', " Flow: ", g_sim%mass_flux_core, " kg/m2·s"

        case ('C')
            g_sim%mass_flux_core = g_sim%mass_flux_core * 0.95_wp
            print '(A,F7.1,A)', " Flow: ", g_sim%mass_flux_core, " kg/m2·s"

        case ('i')
            g_sim%rod_bank_position = min(1.0_wp, g_sim%rod_bank_position + 0.005_wp)
            print '(A,F6.1,A)', " Rods: ", g_sim%rod_bank_position * 100.0_wp, " % inserted"

        case ('I')
            g_sim%rod_bank_position = max(0.0_wp, g_sim%rod_bank_position - 0.005_wp)
            print '(A,F6.1,A)', " Rods: ", g_sim%rod_bank_position * 100.0_wp, " % inserted"

        case ('t')
            g_sim%turbine_valve = min(100.0_wp, g_sim%turbine_valve + 2.0_wp)
            print '(A,F5.1,A)', " Turbine valve: ", g_sim%turbine_valve, " %"

        case ('T')
            g_sim%turbine_valve = max(0.0_wp, g_sim%turbine_valve - 2.0_wp)
            print '(A,F5.1,A)', " Turbine valve: ", g_sim%turbine_valve, " %"

        case ('b')
            g_sim%bypass_valve = min(100.0_wp, g_sim%bypass_valve + 2.0_wp)
            print '(A,F5.1,A)', " Bypass valve:  ", g_sim%bypass_valve, " %"

        case ('B')
            g_sim%bypass_valve = max(0.0_wp, g_sim%bypass_valve - 2.0_wp)
            print '(A,F5.1,A)', " Bypass valve:  ", g_sim%bypass_valve, " %"

        case ('s', 'S')
            call print_steady_state_summary(g_sim)

        case ('h', '?')
            call print_help()

        case ('q', 'Q')
            if (trim(adjustl(cmd)) == 'q'  .or. &
                trim(adjustl(cmd)) == 'Q'  .or. &
                trim(adjustl(cmd)) == 'quit') then
                running = .false.
            end if

        case default
            ! Silently ignore unknown single characters; check for "quit"
            if (trim(adjustl(cmd)) == 'quit' .or. &
                trim(adjustl(cmd)) == 'exit') then
                running = .false.
            end if

        end select

    end subroutine handle_command

    ! ================================================================
    !  HELP TEXT
    ! ================================================================

    subroutine print_help()
        print *, ""
        print *, "=== COMMANDS (type at 'sim>' prompt) ==="
        print *, "  <Enter>  - advance", STEPS_PER_PROMPT, "steps"
        print *, "  p        - pause / resume"
        print *, "  1-5      - display focus (Power/Temp/Void/Burnup/T_sat)"
        print *, "  +/-      - increase / decrease power setpoint 5 %"
        print *, "  r/R      - insert +20 / -20 pcm reactivity"
        print *, "  c/C      - increase / decrease coolant flow 5 %"
        print *, "  i/I      - insert / withdraw control rods 0.5 %"
        print *, "  t/T      - open / close turbine valve 2 %"
        print *, "  b/B      - open / close bypass valve 2 %"
        print *, "  s        - steady-state summary"
        print *, "  h / ?    - this help"
        print *, "  q / quit - quit"
        print *, "========================================"
        print *, ""
    end subroutine print_help

    ! ================================================================
    !  PURE PHYSICS FUNCTIONS  (unchanged from original)
    ! ================================================================

    !> Saturation temperature as a function of pressure.
    !  Polynomial fit to steam tables.  Returns temperature in Kelvin.
    pure function sat_temp_K(P_Pa) result(T_K)
        real(wp), intent(in) :: P_Pa
        real(wp) :: T_K, p
        p   = P_Pa / 1.0e6_wp   ! convert to MPa
        T_K = (-4.0964e-3_wp * p**6  &
              +  0.141738_wp  * p**5  &
              -  1.943057_wp  * p**4  &
              + 13.50428_wp   * p**3  &
              - 51.27379_wp   * p**2  &
              + 118.6854_wp   * p     &
              + 100.1542_wp) + 273.15_wp   ! °C → K
    end function sat_temp_K

    !> Liquid water density polynomial
    pure function coolant_density_gl(T_K) result(rho)
        real(wp), intent(in) :: T_K
        real(wp) :: rho, T_C
        T_C = T_K - 273.15_wp
        rho = -4.467711e-6_wp   * T_C**3  &
              - 5.60288485e-4_wp * T_C**2  &
              - 0.429148844451_wp* T_C     &
              + 1010.035413387815_wp
        rho = max(0.1_wp, rho)
    end function coolant_density_gl

    !> Reactor water level above core floor
    pure function calc_reactor_level(Lbase, T_K, H_core) result(Lrx)
        real(wp), intent(in) :: Lbase, T_K, H_core
        real(wp) :: Lrx, rho
        rho = coolant_density_gl(T_K)
        Lrx = (Lbase / (rho / 1000.0_wp)) - H_core
    end function calc_reactor_level

    !> Void reactivity factor from effective moderator density.
    pure function void_reactivity_factor(avg_void_frac) result(rf)
        real(wp), intent(in) :: avg_void_frac
        real(wp) :: rf, dens2
        dens2 = (1.0_wp - min(1.0_wp, max(0.0_wp, avg_void_frac))) * 1000.0_wp
        rf = -5.84e-10_wp   * dens2**3  &
             + 1.35422e-7_wp * dens2**2  &
             + 1.358252042e-3_wp * dens2  &
             + 0.090469568418_wp
    end function void_reactivity_factor

    ! ================================================================
    !  GEOMETRY SETUP  (unchanged)
    ! ================================================================

    subroutine setup_bwr_geometry_realistic(sim)
        type(simulation_t), intent(inout) :: sim

        integer  :: i, j, k
        real(wp) :: x, y, r, x_center, y_center, pitch, pin_radius
        integer  :: pin_i, pin_j
        real(wp) :: dx_from_pin, dy_from_pin, r_from_pin
        logical  :: in_fuel, found_water, found_fuel

        type(xsec_material_t) :: xsec_fuel, xsec_water

        call xslib_get_material(sim%xslib, "UO2_35", xsec_fuel, found_fuel)
        if (.not. found_fuel) then
            print *, "Fatal: material UO2_35 not found."
            call xslib_list_materials(sim%xslib); stop 1
        end if
        call xslib_get_material(sim%xslib, "H2O", xsec_water, found_water)
        if (.not. found_water) then
            print *, "Fatal: material H2O not found."
            call xslib_list_materials(sim%xslib); stop 1
        end if

        pitch      = 0.0163_wp  ! 1.63 cm pin pitch (BWR typical)
        pin_radius = 0.0052_wp  ! 5.2 mm fuel pellet radius

        print *, "BWR lattice parameters:"
        print '(A,F7.4,A)', "  Pin pitch:         ", pitch*100, " cm"
        print '(A,F7.4,A)', "  Pin radius:        ", pin_radius*100, " cm"
        print '(A,F7.2)',   "  Fuel-to-mod ratio: ", &
            (3.14159_wp*pin_radius**2) / (pitch**2 - 3.14159_wp*pin_radius**2)

        do k = 1, sim%nz
            do j = 1, sim%ny
                do i = 1, sim%nx
                    x = (real(i, wp) - 0.5_wp) * sim%dx
                    y = (real(j, wp) - 0.5_wp) * sim%dy

                    pin_i = nint(x / pitch);  pin_j = nint(y / pitch)
                    x_center = real(pin_i, wp) * pitch
                    y_center = real(pin_j, wp) * pitch

                    dx_from_pin = x - x_center
                    dy_from_pin = y - y_center
                    r_from_pin  = sqrt(dx_from_pin**2 + dy_from_pin**2)

                    in_fuel = (r_from_pin < pin_radius)

                    ! Cylindrical core boundary
                    r = sqrt((x - sim%core_diameter/2)**2 + (y - sim%core_diameter/2)**2)
                    if (r > sim%core_diameter / 2.2_wp) in_fuel = .false.

                    if (in_fuel) then
                        call mg_set_cross_sections(sim%neutronics, xsec_fuel%xsec_base, &
                            i, i, j, j, k, k)
                    else
                        call mg_set_cross_sections(sim%neutronics, xsec_water%xsec_base, &
                            i, i, j, j, k, k)
                    end if
                end do
            end do
        end do
        print *, "Pin-by-pin geometry set up successfully"
    end subroutine setup_bwr_geometry_realistic

    ! ================================================================
    !  SIMULATION SETUP  (unchanged)
    ! ================================================================

    subroutine setup_bwr_simulation(sim)
        type(simulation_t), intent(out) :: sim

        type(mg_config_t)        :: neutron_config
        type(heat_config_t)      :: heat_config
        type(two_phase_config_t) :: th_config
        type(burnup_config_t)    :: burnup_config
        type(xsec_material_t)    :: xsec_fuel, xsec_water

        ! ── Core geometry ─────────────────────────────────────────
        sim%nx = 20;  sim%ny = 20;  sim%nz = 20
        sim%core_height   = 3.71_wp   ! m
        sim%core_diameter = 4.75_wp   ! m
        sim%dx = sim%core_diameter / real(sim%nx, wp)
        sim%dy = sim%core_diameter / real(sim%ny, wp)
        sim%dz = sim%core_height   / real(sim%nz, wp)

        ! ── Operating conditions ──────────────────────────────────
        sim%power_rated         = 2381.0e6_wp  ! W thermal
        sim%pressure_operating  = 7.14e6_wp    ! Pa (dynamic)
        sim%mass_flux_core      = 1500.0_wp    ! kg/m²·s
        sim%inlet_temperature   = 551.0_wp     ! K

        ! ── Feedback coefficients ─────────────────────────────────
        sim%alpha_doppler = -3.5_wp   ! pcm/K
        sim%alpha_void    = -80.0_wp  ! pcm/(% void)

        sim%rod_bank_position = 0.95_wp  ! start 95 % inserted

        ! ── Physics initialisations ──────────────────
        sim%sat_temperature  = sat_temp_K(sim%pressure_operating)
        ! Normal BWR setpoint: ~0.5 m above top of active fuel.
        sim%Lbase            = (0.5_wp + sim%core_height) &
                             * (coolant_density_gl(sim%inlet_temperature) / 1000.0_wp)
        sim%Lrx              = calc_reactor_level(sim%Lbase, sim%inlet_temperature, sim%core_height)
        sim%avg_coolant_temp = sim%inlet_temperature
        sim%coolant_T_prev1  = sim%inlet_temperature
        sim%coolant_T_prev2  = sim%inlet_temperature
        sim%P_prev1          = sim%pressure_operating
        sim%P_prev2          = sim%pressure_operating
        sim%T_excess         = 0.0_wp

        sim%neutrons_prev  = 0.0_wp
        sim%reactor_period = 1.0e6_wp
        sim%aprm           = 0.0_wp

        sim%rf_void = void_reactivity_factor(0.0_wp)

        ! Turbine
        sim%turbine_valve   = 33.0_wp
        sim%bypass_valve    = 0.0_wp
        sim%steam_flow      = 0.0_wp
        sim%turbine_speed   = 3600.0_wp
        sim%turbine_power   = 0.0_wp
        sim%turbine_tripped = .false.
        sim%feedwater_temp  = sim%inlet_temperature

        print *, "Initialising physics modules..."

        ! ── Cross-section library ─────────────────────────────────
        call xslib_init(sim%xslib, n_groups=2)

        ! Fuel cross-sections
        xsec_fuel%name       = "UO2_35"
        xsec_fuel%n_groups   = 2
        xsec_fuel%is_fuel    = .true.
        xsec_fuel%enrichment = 0.035_wp
        xsec_fuel%T_ref      = 900.0_wp
        xsec_fuel%rho_ref    = 10.97_wp

        call xslib_create_two_group_fuel(xsec_fuel%xsec_base, enrichment=0.035_wp)

        allocate(xsec_fuel%alpha_D(2),   xsec_fuel%alpha_mod(2))
        allocate(xsec_fuel%alpha_rho(2), xsec_fuel%alpha_void(2))
        xsec_fuel%alpha_D    = [-2.0e-5_wp, -3.0e-5_wp]
        xsec_fuel%alpha_mod  = [ 0.0_wp,    0.0_wp    ]
        xsec_fuel%alpha_rho  = [ 0.0_wp,    0.0_wp    ]
        xsec_fuel%alpha_void = [ 0.0_wp,    0.0_wp    ]
        call xslib_add_material(sim%xslib, xsec_fuel)

        ! Water cross-sections
        xsec_water%name    = "H2O"
        xsec_water%n_groups= 2
        xsec_water%is_fuel = .false.
        xsec_water%T_ref   = 560.0_wp
        xsec_water%rho_ref = 0.74_wp

        call xslib_create_two_group_moderator(xsec_water%xsec_base)

        allocate(xsec_water%alpha_D(2),   xsec_water%alpha_mod(2))
        allocate(xsec_water%alpha_rho(2), xsec_water%alpha_void(2))
        xsec_water%alpha_D    = [ 0.0_wp,    0.0_wp    ]
        xsec_water%alpha_mod  = [ 0.0_wp,    1.0e-4_wp ]
        xsec_water%alpha_rho  = [-10.0_wp,  -50.0_wp   ]
        xsec_water%alpha_void = [-10.0_wp, -100.0_wp   ]
        call xslib_add_material(sim%xslib, xsec_water)

        ! ── Neutronics ────────────────────────────────────────────
        neutron_config%n_groups         = 2
        neutron_config%max_outer_iter   = 100
        neutron_config%outer_tolerance  = 1.0e-5_wp
        neutron_config%power_level      = sim%power_rated
        neutron_config%normalize_power  = .true.
        call mg_init(sim%neutronics, sim%nx, sim%ny, sim%nz, &
             sim%dx*100.0_wp, sim%dy*100.0_wp, sim%dz*100.0_wp, neutron_config)
        call setup_bwr_geometry_realistic(sim)

        ! ── Heat transfer ─────────────────────────────────────────
        heat_config%include_convection = .true.
        call heat_init(sim%heat, sim%nx, sim%ny, sim%nz, &
                       sim%dx, sim%dy, sim%dz, heat_config)
        call heat_set_properties(sim%heat, k=3.0_wp, rho=10970.0_wp, cp=300.0_wp, &
            i1=1, i2=int(0.9*sim%nx), j1=1, j2=int(0.9*sim%ny))
        call heat_set_properties(sim%heat, k=0.6_wp, rho=738.0_wp, cp=5200.0_wp, &
            i1=int(0.9*sim%nx)+1, i2=sim%nx, &
            j1=int(0.9*sim%ny)+1, j2=sim%ny)

        ! ── Two-phase thermal-hydraulics ──────────────────────────
        th_config%void_correlation          = VOID_CHEXAL_LELLOUCHE_ID
        th_config%chf_correlation           = CHF_GROENEVELD_ID
        th_config%include_subcooled_boiling = .true.
        call two_phase_init(sim%thermalhydraulics, sim%nx, sim%ny, sim%nz, &
                            sim%dx, sim%dy, sim%dz, th_config)
        call two_phase_set_geometry(sim%thermalhydraulics, diameter=0.012_wp)

        ! ── Burnup / depletion ────────────────────────────────────
        burnup_config%track_xenon     = .true.
        burnup_config%track_samarium  = .true.
        burnup_config%track_actinides = .true.
        call burnup_init(sim%burnup, sim%nx, sim%ny, sim%nz, &
                         sim%dx, sim%dy, sim%dz, burnup_config)
        call burnup_set_initial_composition(sim%burnup, enrichment=0.035_wp)

        ! ── Time counters ─────────────────────────────────────────
        sim%time      = 0.0_wp
        sim%n_steps   = 0
        sim%initialized = .false.

        print *, "  Grid:          ", sim%nx, "x", sim%ny, "x", sim%nz
        print *, "  Core diameter: ", sim%core_diameter, " m"
        print *, "  Core height:   ", sim%core_height,   " m"
        print *, "  Rated power:   ", sim%power_rated/1.0e6_wp, " MW"
        print '(A,F6.1,A)', "  Sat. temp at ", sim%pressure_operating/1.0e6_wp, &
            " MPa: " // trim(adjustl(fmt_real(sim%sat_temperature - 273.15_wp))) // " degC"
        print *, ""

    end subroutine setup_bwr_simulation

    ! ================================================================
    !  PRESSURE DYNAMICS  (unchanged)
    ! ================================================================

    subroutine update_pressure_dynamics(sim, dt)
        type(simulation_t), intent(inout) :: sim
        real(wp), intent(in) :: dt

        real(wp) :: T_sat, T_mean_smooth, rel_T, dP_gen, steam_val, dP_loss
        real(wp), parameter :: dP_COEF    = 5000.0_wp
        real(wp), parameter :: P_LOSS_COEF= 1.125e6_wp
        real(wp), parameter :: P_NOMINAL  = 7.14e6_wp

        sim%coolant_T_prev2 = sim%coolant_T_prev1
        sim%coolant_T_prev1 = sim%avg_coolant_temp
        T_mean_smooth = (sim%avg_coolant_temp + sim%coolant_T_prev1 + sim%coolant_T_prev2) / 3.0_wp

        T_sat = sat_temp_K(sim%pressure_operating)
        sim%sat_temperature = T_sat

        sim%T_excess = max(0.0_wp, T_mean_smooth - T_sat)
        rel_T = (T_mean_smooth - 373.15_wp) * (1.0_wp - sim%pressure_operating / 10.0e6_wp)
        dP_gen = dP_COEF * (rel_T + 30.0_wp * sim%T_excess) * dt

        if (.not. sim%turbine_tripped) then
            steam_val = (sim%turbine_valve + sim%bypass_valve) / 100.0_wp &
                        * (sim%pressure_operating / P_NOMINAL)
        else
            steam_val = 0.0_wp
        end if
        dP_loss = -P_LOSS_COEF * steam_val * dt

        sim%P_prev2 = sim%P_prev1
        sim%P_prev1 = sim%pressure_operating
        sim%pressure_operating = max(3.0e6_wp, min(11.0e6_wp, &
            sim%pressure_operating + dP_gen + dP_loss))

        sim%Lbase = sim%Lbase &
                  - (sim%power_current / max(1.0_wp, 2381.0e6_wp)) * 0.002_wp * dt
        sim%Lbase = sim%Lbase + (sim%mass_flux_core / 1500.0_wp) * 0.0015_wp * dt
        sim%Lbase = max(0.5_wp, min(sim%core_height + 3.0_wp, sim%Lbase))

        sim%Lrx = calc_reactor_level(sim%Lbase, sim%avg_coolant_temp, sim%core_height)

        sim%feedwater_temp = sim%inlet_temperature &
                           + (sim%avg_coolant_temp - sim%inlet_temperature) * 0.1_wp

    end subroutine update_pressure_dynamics

    ! ================================================================
    !  INSTRUMENTATION  (unchanged)
    ! ================================================================

    subroutine update_instrumentation(sim)
        type(simulation_t), intent(inout) :: sim
        real(wp) :: ratio

        sim%aprm = (sim%power_current / max(1.0_wp, sim%power_rated)) * 100.0_wp
        sim%aprm = max(0.0_wp, min(200.0_wp, sim%aprm))

        if (sim%neutrons_prev > 1.0e3_wp .and. sim%power_current > 1.0e3_wp) then
            ratio = sim%power_current / sim%neutrons_prev
            if (ratio > 0.0_wp .and. abs(ratio - 1.0_wp) > 1.0e-9_wp) then
                sim%reactor_period = 1.0_wp / log(ratio)
                sim%reactor_period = max(-999.0_wp, min(999.0_wp, sim%reactor_period))
            else
                sim%reactor_period = 9999.0_wp
            end if
        end if
        sim%neutrons_prev = sim%power_current

        sim%rf_void = void_reactivity_factor(sim%avg_void_fraction / 100.0_wp)

    end subroutine update_instrumentation

    ! ================================================================
    !  TURBINE MODEL  (unchanged)
    ! ================================================================

    subroutine update_turbine(sim, dt)
        type(simulation_t), intent(inout) :: sim
        real(wp), intent(in) :: dt

        real(wp) :: eta, ePres, steam_to_turbine
        real(wp), parameter :: ETA_CARNOT = 0.34_wp
        real(wp), parameter :: RATED_RPM  = 3600.0_wp
        real(wp), parameter :: P_NOMINAL  = 7.14e6_wp

        if (sim%turbine_tripped) then
            sim%turbine_speed = max(0.0_wp, sim%turbine_speed - 30.0_wp * dt)
            sim%turbine_power = 0.0_wp
            sim%steam_flow    = 0.0_wp
            return
        end if

        ePres = min(1.0_wp, 6.1e-8_wp * sim%pressure_operating + 0.567_wp)

        steam_to_turbine = sim%turbine_valve / 100.0_wp &
                         * 2000.0_wp * sim%pressure_operating / 10.0e6_wp
        sim%steam_flow = steam_to_turbine

        if (sim%turbine_speed > 3400.0_wp) then
            sim%turbine_speed = RATED_RPM
            eta = ETA_CARNOT * ePres * (sim%turbine_valve / 100.0_wp)
            sim%turbine_power = eta * sim%power_current
        else
            sim%turbine_speed = sim%turbine_speed + &
                (RATED_RPM * steam_to_turbine / 2000.0_wp - sim%turbine_speed) / 30.0_wp * dt
            sim%turbine_speed = max(0.0_wp, sim%turbine_speed)
            eta = ETA_CARNOT * ePres * (sim%turbine_speed / RATED_RPM) * (sim%turbine_valve / 100.0_wp)
            sim%turbine_power = eta * sim%power_current
        end if

        if (sim%turbine_tripped) sim%turbine_power = 0.0_wp

    end subroutine update_turbine

    ! ================================================================
    !  CROSS-SECTION FEEDBACK  (unchanged)
    ! ================================================================

    subroutine update_cross_sections_feedback(sim, T, rho)
        type(simulation_t), intent(inout) :: sim
        real(wp), intent(in) :: T(:, :, :)
        real(wp), intent(in) :: rho(:, :, :)

        integer  :: i, j, k
        type(mg_xsec_t) :: xsec
        real(wp) :: T_fuel, rho_mod, rf, rf_ref, rho_mod_corrected
        real(wp) :: node_bottom, node_top, rod_tip, inserted_fraction, H_core

        rf_ref = void_reactivity_factor(0.30_wp)
        rf     = sim%rf_void

        do k = 1, sim%nz
            do j = 1, sim%ny
                do i = 1, sim%nx
                    T_fuel  = T(i, j, k)
                    rho_mod = rho(i, j, k)

                    rho_mod_corrected = rho_mod * (rf / max(1.0e-9_wp, rf_ref))

                    call xslib_get_xsec(sim%xslib, "UO2_35", &
                        T_fuel, rho_mod_corrected, sim%burnup%burnup(i, j, k), xsec, &
                        Xe_conc = sim%burnup%Xe135(i, j, k), &
                        Sm_conc = sim%burnup%Sm149(i, j, k))

                    H_core      = real(sim%nz, wp) * sim%dz
                    node_bottom = real(k-1, wp) * sim%dz
                    node_top    = real(k,   wp) * sim%dz
                    rod_tip     = sim%rod_bank_position * H_core

                    inserted_fraction = 0.0_wp
                    if (rod_tip >= node_top) then
                        inserted_fraction = 1.0_wp
                    else if (rod_tip > node_bottom) then
                        inserted_fraction = (rod_tip - node_bottom) / sim%dz
                    end if

                    if (inserted_fraction > 0.0_wp) &
                        call xslib_apply_control_rod(xsec, inserted_fraction)

                    call mg_set_cross_sections(sim%neutronics, xsec, i, i, j, j, k, k)
                end do
            end do
        end do
    end subroutine update_cross_sections_feedback

    ! ================================================================
    !  COUPLED TIME STEP  (unchanged)
    ! ================================================================

    subroutine coupled_time_step(sim, dt)
        type(simulation_t), intent(inout) :: sim
        real(wp), intent(in) :: dt

        real(wp) :: power(sim%nx, sim%ny, sim%nz)
        real(wp) :: flux(sim%nx, sim%ny, sim%nz)
        real(wp) :: temperature(sim%nx, sim%ny, sim%nz)
        real(wp) :: void_fraction(sim%nx, sim%ny, sim%nz)
        real(wp) :: density(sim%nx, sim%ny, sim%nz)
        real(wp), parameter :: rho_liquid = 738.0_wp
        real(wp), parameter :: rho_vapor  = 0.038_wp
        logical  :: converged
        real(wp) :: v_coolant, Re, h_conv
        real(wp), parameter :: D_h     = 0.01_wp
        real(wp), parameter :: mu      = 0.0001_wp
        real(wp), parameter :: k_fluid = 0.6_wp
        real(wp), parameter :: Pr      = 0.9_wp

        v_coolant = sim%mass_flux_core / rho_liquid
        if (allocated(sim%heat%vz)) sim%heat%vz = v_coolant

        Re     = rho_liquid * v_coolant * D_h / mu
        h_conv = 0.023_wp * Re**0.8_wp * Pr**0.4_wp * k_fluid / D_h

        call mg_solve_eigenvalue(sim%neutronics, sim%k_eff, converged)
        call mg_get_power(sim%neutronics, power)
        call mg_get_flux(sim%neutronics, flux)
        sim%power_current = sim%neutronics%total_power

        call burnup_step(sim%burnup, flux, power, dt)
        sim%avg_burnup = sum(sim%burnup%burnup) / size(sim%burnup%burnup)

        sim%heat%q = power
        call heat_step(sim%heat, dt)
        temperature       = sim%heat%T
        sim%max_fuel_temp = maxval(temperature)

        sim%avg_coolant_temp = sum(temperature) / real(size(temperature), wp)

        call two_phase_step(sim%thermalhydraulics, &
            temperature, &
            sim%pressure_operating + 0.0_wp * temperature, &
            sim%mass_flux_core     + 0.0_wp * temperature, &
            power / (sim%dx * sim%dy), dt)

        call two_phase_get_void_fraction(sim%thermalhydraulics, void_fraction)

        density = (1.0_wp - void_fraction) * rho_liquid + void_fraction * rho_vapor
        density = max(density, rho_vapor)
        density = min(density, rho_liquid)

        sim%avg_void_fraction = sum(void_fraction) / size(void_fraction) * 100.0_wp
        sim%max_void_fraction = maxval(void_fraction) * 100.0_wp

        if (allocated(sim%thermalhydraulics%chf_ratio)) then
            sim%min_chfr = minval(sim%thermalhydraulics%chf_ratio)
        else
            sim%min_chfr = 999.0_wp
        end if

        call update_cross_sections_feedback(sim, temperature, density)
        call update_instrumentation(sim)
        call update_pressure_dynamics(sim, dt)
        call update_turbine(sim, dt)

        sim%reactivity_pcm = (sim%k_eff - 1.0_wp) / sim%k_eff * 1.0e5_wp

        sim%time    = sim%time + dt
        sim%n_steps = sim%n_steps + 1

    end subroutine coupled_time_step

    ! ================================================================
    !  STEADY-STATE SOLVE  (unchanged)
    ! ================================================================

    subroutine solve_steady_state(sim)
        type(simulation_t), intent(inout) :: sim

        integer  :: iter
        real(wp) :: power(sim%nx, sim%ny, sim%nz)
        real(wp) :: temperature(sim%nx, sim%ny, sim%nz)
        real(wp) :: void_fraction(sim%nx, sim%ny, sim%nz)
        real(wp) :: density(sim%nx, sim%ny, sim%nz)
        real(wp) :: error
        real(wp), parameter :: rho_liquid = 0.738_wp
        real(wp), parameter :: rho_vapor  = 0.038_wp
        logical  :: converged

        do iter = 1, 50
            call mg_solve_eigenvalue(sim%neutronics, sim%k_eff, converged)
            if (.not. converged) &
                print *, "Warning: neutronics did not converge at iteration", iter

            call mg_get_power(sim%neutronics, power)
            sim%power_current = sim%neutronics%total_power

            sim%heat%q = power
            call heat_step_implicit(sim%heat, 1.0_wp)
            temperature       = sim%heat%T
            sim%max_fuel_temp = maxval(temperature)

            sim%avg_coolant_temp = sum(temperature) / real(size(temperature), wp)

            call two_phase_step(sim%thermalhydraulics, &
                temperature, &
                sim%pressure_operating + 0.0_wp * temperature, &
                sim%mass_flux_core     + 0.0_wp * temperature, &
                power / (sim%dx * sim%dy), 1.0_wp)

            call two_phase_get_void_fraction(sim%thermalhydraulics, void_fraction)

            density = (1.0_wp - void_fraction) * rho_liquid + void_fraction * rho_vapor
            density = max(density, rho_vapor)
            density = min(density, rho_liquid)

            sim%avg_void_fraction = sum(void_fraction) / size(void_fraction) * 100.0_wp
            sim%max_void_fraction = maxval(void_fraction) * 100.0_wp

            call update_cross_sections_feedback(sim, temperature, density)

            error = abs(sim%k_eff - 1.0_wp) + &
                    abs(sim%power_current - sim%power_rated) / sim%power_rated
            if (error < 1.0e-4_wp .and. iter > 3) exit
        end do

        sim%reactivity_pcm = (sim%k_eff - 1.0_wp) / sim%k_eff * 1.0e5_wp

        if (allocated(sim%thermalhydraulics%chf_ratio)) then
            sim%min_chfr = minval(sim%thermalhydraulics%chf_ratio)
        else
            sim%min_chfr = 999.0_wp
        end if

        call update_instrumentation(sim)
        call update_pressure_dynamics(sim, 0.0_wp)
        call update_turbine(sim, 0.0_wp)

    end subroutine solve_steady_state

    ! ================================================================
    !  REACTIVITY INSERTION  (unchanged)
    ! ================================================================

    subroutine apply_reactivity_insertion(sim, rho_pcm)
        type(simulation_t), intent(inout) :: sim
        real(wp), intent(in) :: rho_pcm

        real(wp) :: factor
        integer  :: i, j, k, g

        factor = 1.0_wp + rho_pcm / 1.0e5_wp
        do k = 1, sim%nz
            do j = 1, sim%ny
                do i = 1, sim%nx
                    do g = 1, sim%neutronics%n_groups
                        sim%neutronics%xsec(i, j, k)%nu_sigma_f(g) = &
                            sim%neutronics%xsec(i, j, k)%nu_sigma_f(g) * factor
                    end do
                end do
            end do
        end do
    end subroutine apply_reactivity_insertion

    ! ================================================================
    !  AUTOMATIC CONTROLS  (unchanged; uses g_paused instead of g_viz%paused)
    ! ================================================================

    subroutine apply_automatic_controls(sim, dt)
        type(simulation_t), intent(inout) :: sim
        real(wp), intent(in) :: dt

        real(wp) :: power_fraction, rod_adjustment
        logical, save :: warning_printed = .false.

        real(wp), parameter :: SCRAM_TEMP        = 1673.15_wp
        real(wp), parameter :: HIGH_TEMP         = 1473.15_wp
        real(wp), parameter :: SCRAM_POWER       = 1.30_wp
        real(wp), parameter :: HIGH_POWER        = 1.15_wp
        real(wp), parameter :: SCRAM_PERIOD      = 5.0_wp
        real(wp), parameter :: SCRAM_LEVEL       = -1.5_wp
        real(wp), parameter :: SCRAM_PRESSURE_HI = 9.5e6_wp
        real(wp), parameter :: SCRAM_PRESSURE_LO = 4.0e6_wp

        power_fraction = sim%power_current / max(sim%power_rated, 1.0e6_wp)

        if (sim%max_fuel_temp > SCRAM_TEMP) then
            call scram(sim, "FUEL TEMPERATURE LIMIT", &
                "Max fuel: " // trim(deg_str(sim%max_fuel_temp)) // "degC  Limit: 1400degC")
            return
        end if

        if (power_fraction > SCRAM_POWER) then
            call scram(sim, "POWER EXCEEDED LIMIT", &
                "Power: " // trim(pct_str(power_fraction*100)) // "%  Limit: 130%")
            return
        end if

        if (sim%reactor_period > 0.0_wp .and. sim%reactor_period < SCRAM_PERIOD &
            .and. sim%power_current > 1.0e6_wp) then
            call scram(sim, "REACTOR PERIOD TOO SHORT", &
                "Period: " // trim(per_str(sim%reactor_period)) // " s  Limit: " // &
                trim(per_str(SCRAM_PERIOD)) // " s")
            return
        end if

        if (sim%Lrx < SCRAM_LEVEL) then
            call scram(sim, "REACTOR WATER LEVEL LOW", &
                "Level: " // trim(lev_str(sim%Lrx)) // " m  Limit: " // &
                trim(lev_str(SCRAM_LEVEL)) // " m")
            return
        end if

        if (sim%pressure_operating > SCRAM_PRESSURE_HI) then
            call scram(sim, "HIGH REACTOR PRESSURE", &
                "Pressure: " // trim(mpa_str(sim%pressure_operating/1.0e6_wp)) // &
                " MPa  Limit: 9.5 MPa")
            return
        end if

        if (sim%pressure_operating < SCRAM_PRESSURE_LO .and. sim%power_current > 0.1e6_wp) then
            call scram(sim, "LOW REACTOR PRESSURE", &
                "Pressure: " // trim(mpa_str(sim%pressure_operating/1.0e6_wp)) // &
                " MPa  Limit: 4.0 MPa")
            return
        end if

        ! ── Automatic power runback ───────────────────────────────
        if (power_fraction > HIGH_POWER) then
            rod_adjustment = (power_fraction - HIGH_POWER) * 0.02_wp * dt / 0.05_wp
            sim%rod_bank_position = min(1.0_wp, sim%rod_bank_position + rod_adjustment)
            if (.not. warning_printed) then
                print *, "  AUTO CONTROL: Power high, inserting rods"
                warning_printed = .true.
            end if
        else if (power_fraction < HIGH_POWER - 0.02_wp) then
            warning_printed = .false.
        end if

        if (sim%max_fuel_temp > HIGH_TEMP) then
            rod_adjustment = (sim%max_fuel_temp - HIGH_TEMP) / 500.0_wp * 0.05_wp * dt / 0.05_wp
            sim%rod_bank_position = min(1.0_wp, sim%rod_bank_position + rod_adjustment)
        end if

        if (sim%reactivity_pcm > 500.0_wp) then
            rod_adjustment = (sim%reactivity_pcm - 500.0_wp) / 10000.0_wp * dt / 0.05_wp
            sim%rod_bank_position = min(1.0_wp, sim%rod_bank_position + rod_adjustment)
        end if

        if (.not. sim%turbine_tripped) then
            call auto_pressure_control(sim)
        end if

    end subroutine apply_automatic_controls

    !> Auto pressure control (PD controller).
    subroutine auto_pressure_control(sim)
        type(simulation_t), intent(inout) :: sim
        real(wp) :: vel, dist, chg
        real(wp), parameter :: P_TARGET = 7.14e6_wp

        vel  = (sim%pressure_operating - (sim%P_prev1 + sim%P_prev2) / 2.0_wp) * 7.0_wp
        dist = sim%pressure_operating - P_TARGET
        chg  = max(-1.0_wp, min(1.0_wp, (dist + vel) / 1.0e5_wp))

        if (chg > 0.0_wp) then
            if (sim%turbine_valve < 100.0_wp) then
                sim%turbine_valve = min(100.0_wp, sim%turbine_valve + chg)
            else
                sim%bypass_valve  = min(100.0_wp, sim%bypass_valve + chg)
            end if
        else
            if (sim%bypass_valve > 0.0_wp) then
                sim%bypass_valve = max(0.0_wp, sim%bypass_valve + chg)
            else
                sim%turbine_valve = max(0.0_wp, sim%turbine_valve + chg)
            end if
        end if
    end subroutine auto_pressure_control

    ! ================================================================
    !  SAFETY LIMITS CHECK  (unchanged)
    ! ================================================================

    function check_safety_limits(sim) result(approaching_limits)
        type(simulation_t), intent(in) :: sim
        logical :: approaching_limits

        approaching_limits = .false.

        if (sim%max_fuel_temp > 1373.0_wp) then
            approaching_limits = .true.
            print *, "  High fuel temperature: ", sim%max_fuel_temp, " K"
        end if
        if (sim%min_chfr < 1.5_wp) then
            approaching_limits = .true.
            print *, "  Low CHFR: ", sim%min_chfr
        end if
        if (sim%power_current > 1.2_wp * sim%power_rated) then
            approaching_limits = .true.
            print *, "  High power: ", sim%power_current / 1.0e6_wp, " MW"
        end if
        if (sim%Lrx < 0.0_wp) then
            approaching_limits = .true.
            print '(A,F6.2,A)', "  Low reactor level: ", sim%Lrx, " m"
        end if
        if (sim%pressure_operating > 9.0e6_wp) then
            approaching_limits = .true.
            print '(A,F6.3,A)', "  High pressure: ", sim%pressure_operating/1.0e6_wp, " MPa"
        end if
    end function check_safety_limits

    ! ================================================================
    !  REPORTING  (unchanged)
    ! ================================================================

    subroutine print_steady_state_summary(sim)
        type(simulation_t), intent(in) :: sim

        print *, ""
        print *, "========== Steady-State Results =========="
        print *, "k_eff:                ", sim%k_eff
        print *, "Reactivity:           ", sim%reactivity_pcm, " pcm"
        print *, "Core power:           ", sim%power_current / 1.0e6_wp, " MW"
        print *, "Max fuel temp:        ", sim%max_fuel_temp - 273.15_wp, " degC"
        print *, "Avg void fraction:    ", sim%avg_void_fraction, " %"
        print *, "Max void fraction:    ", sim%max_void_fraction, " %"
        print *, "Min CHFR:             ", sim%min_chfr
        print *, "Sat. temperature:     ", sim%sat_temperature - 273.15_wp, " degC"
        print *, "Reactor level (Lrx):  ", sim%Lrx, " m above core floor"
        print *, "APRM:                 ", sim%aprm, " %"
        print *, "Void reactivity (rf): ", sim%rf_void
        print *, "Turbine valve:        ", sim%turbine_valve, " %"
        print *, "Turbine power:        ", sim%turbine_power / 1.0e6_wp, " MWe"
        print *, "Pressure:             ", sim%pressure_operating / 1.0e6_wp, " MPa"
        print *, "=========================================="
        print *, ""
    end subroutine print_steady_state_summary

    subroutine print_transient_summary(sim, step, t)
        type(simulation_t), intent(in) :: sim
        integer,  intent(in) :: step
        real(wp), intent(in) :: t

        write(*, '(A,I6,A,F8.2,A)') "Step ", step, ", t = ", t, " s"
        write(*, '(A,F10.5,A,F8.1,A)') &
            "  k_eff = ", sim%k_eff, ", rho = ", sim%reactivity_pcm, " pcm"
        write(*, '(A,F8.1,A,A,F5.1,A)') &
            "  Power = ", sim%power_current/1.0e6_wp, " MW", &
            "  APRM = ", sim%aprm, " %"
        write(*, '(A,F7.1,A)') &
            "  Max T_fuel = ", sim%max_fuel_temp - 273.15_wp, " degC"
        write(*, '(A,F6.2,A,A,F5.2)') &
            "  Avg void = ", sim%avg_void_fraction, " %", "  Min CHFR = ", sim%min_chfr
        write(*, '(A,F6.2,A,A,F6.3,A)') &
            "  Level = ", sim%Lrx, " m", &
            "  Pressure = ", sim%pressure_operating/1.0e6_wp, " MPa"
        write(*, '(A,F8.1,A,A,F8.2,A)') &
            "  Turbine = ", sim%turbine_power/1.0e6_wp, " MWe", &
            "  Period = ", sim%reactor_period, " s"
        write(*, '(A,F6.2)') &
            "  rf_void = ", sim%rf_void
        print *, ""
    end subroutine print_transient_summary

    ! ================================================================
    !  CLEANUP  (unchanged)
    ! ================================================================

    subroutine cleanup_simulation(sim)
        type(simulation_t), intent(inout) :: sim

        if (.not. sim%initialized) return

        call mg_destroy(sim%neutronics)
        call heat_destroy(sim%heat)
        call two_phase_destroy(sim%thermalhydraulics)
        call burnup_destroy(sim%burnup)
        call xslib_destroy(sim%xslib)

        sim%initialized = .false.
    end subroutine cleanup_simulation

    ! ================================================================
    !  HELPER PROCEDURES  (unchanged)
    ! ================================================================

    function fmt_real(v) result(s)
        real(wp), intent(in) :: v
        character(len=8) :: s
        write(s,'(F8.1)') v
    end function fmt_real

    !> Emergency SCRAM: full rod insertion + turbine trip + pause.
    !  Uses module-level g_paused instead of g_viz%paused.
    subroutine scram(sim, reason, detail)
        type(simulation_t), intent(inout) :: sim
        character(len=*),   intent(in)    :: reason, detail
        print *, ""
        print *, "╔════════════════════════════════════════════╗"
        print *, "║   *** EMERGENCY REACTOR SCRAM ***          ║"
        print *, "║   CAUSE: " // reason
        print *, "║   " // detail
        print *, "║   ACTION: Full rod insertion               ║"
        print *, "╚════════════════════════════════════════════╝"
        print *, ""
        sim%rod_bank_position = 1.0_wp
        sim%turbine_tripped   = .true.
        g_paused              = .true.   ! pause for operator review
    end subroutine scram

    function deg_str(v) result(s)
        real(wp), intent(in) :: v; character(len=8) :: s
        write(s,'(F8.1)') v - 273.15_wp
    end function deg_str

    function pct_str(v) result(s)
        real(wp), intent(in) :: v; character(len=8) :: s
        write(s,'(F8.1)') v
    end function pct_str

    function per_str(v) result(s)
        real(wp), intent(in) :: v; character(len=8) :: s
        write(s,'(F8.2)') v
    end function per_str

    function lev_str(v) result(s)
        real(wp), intent(in) :: v; character(len=8) :: s
        write(s,'(F8.2)') v
    end function lev_str

    function mpa_str(v) result(s)
        real(wp), intent(in) :: v; character(len=8) :: s
        write(s,'(F8.3)') v
    end function mpa_str

end program bwr_core_simulator