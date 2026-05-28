! Residual Heat Removal subsystem.
!
! Models the four-loop RHR system in lumped-parameter form. Each loop
! has a centrifugal pump (rated 7700 gpm @ 220 ft TDH for BWR/4) and
! draws from one of:
!   * recirculation suction (shutdown cooling mode, SDC)
!   * suppression pool (pool cooling, containment spray, LPCI test)
! and discharges through a heat exchanger (RHR HX, UA ≈ 8.5e6 W/K
! per HX at design service-water temp) to:
!   * feedwater line / reactor (SDC, LPCI)
!   * suppression pool spray header (pool cooling)
!   * drywell + wetwell spray nozzles (containment spray)
!   * recirc loop return (LPCI normal alignment)
!
! Modes (one global mode applies to all running loops; in real plant
! each loop can be aligned independently, but for the sim we fold them
! into a single operator-selectable mode):
!   * RHR_MODE_STANDBY: pumps idle, valves lined up for LPCI auto-init.
!   * RHR_MODE_SHUTDOWN_COOLING: SDC. Reactor coolant inventory cooled
!     directly via the RHR HX. Permissive: dome P < 135 psig (~1 MPa)
!     so HX shell side is not over-pressurised.
!   * RHR_MODE_SUPP_POOL_COOLING: pool water through HX back to pool.
!   * RHR_MODE_CONTAINMENT_SPRAY: pool water through HX into spray
!     headers (drywell + wetwell). Available only with drywell-pressure
!     permissive (post-LOCA).
!   * RHR_MODE_LPCI: low-pressure coolant injection. Auto-initiated by
!     ECCS logic on Level 1 + drywell high P, or manual. Injects from
!     suppression pool into the recirc loop discharge.
!
! Heat removal rate per running loop:
!     Q_loop = UA * (T_inlet − T_service)
! where T_service is the service-water (RHRSW) supply temp (default
! 305 K = 90 °F service water). T_inlet depends on mode:
!     SDC                        → reactor coolant T (avg_coolant_temp_K)
!     SUPP_POOL / CONT_SPRAY     → suppression pool T (driver)
!     LPCI                       → suppression pool T (initially)
!     STANDBY                    → 0 (no flow)
!
! The Q_loop sink reduces vessel power or pool temperature as
! appropriate; the orchestrator wires Q_total to the right state field.
module rhr
    use kinds, only: wp
    implicit none
    private

    public :: rhr_state_t, rhr_config_t
    public :: rhr_command_t, rhr_observation_t
    public :: rhr_drivers_t
    public :: rhr_init, rhr_destroy
    public :: rhr_step, rhr_observe
    public :: rhr_apply
    public :: RHR_N_LOOPS
    public :: RHR_MODE_STANDBY, RHR_MODE_SHUTDOWN_COOLING
    public :: RHR_MODE_SUPP_POOL_COOLING, RHR_MODE_CONTAINMENT_SPRAY, RHR_MODE_LPCI

    integer, parameter :: RHR_N_LOOPS = 4

    integer, parameter :: RHR_MODE_STANDBY              = 0
    integer, parameter :: RHR_MODE_SHUTDOWN_COOLING     = 1
    integer, parameter :: RHR_MODE_SUPP_POOL_COOLING    = 2
    integer, parameter :: RHR_MODE_CONTAINMENT_SPRAY    = 3
    integer, parameter :: RHR_MODE_LPCI                 = 4

    real(wp), parameter :: PA_PER_PSI = 6894.76_wp
    real(wp), parameter :: P_ATM_PA   = 101325.0_wp

    type :: rhr_config_t
        real(wp) :: rated_flow_kg_s_per_loop = 486.0_wp  ! 7700 gpm × 1 kg/L × 3.785 L/gal / 60
        real(wp) :: UA_W_per_K_per_loop      = 8.5e6_wp  ! design RHR HX duty
        real(wp) :: service_water_T_K        = 305.0_wp  ! 90 °F RHRSW intake
        ! SDC permissive: reactor pressure must drop below this before
        ! SDC can be aligned. Real plant: ~135 psig isolation valve
        ! protects HX tubes. Treat as a hard mode-rejection threshold.
        real(wp) :: sdc_max_pressure_Pa      = 135.0_wp * PA_PER_PSI + P_ATM_PA
        ! Containment spray permissive: drywell pressure above this
        ! arms the spray-header valves. Operator can override below.
        real(wp) :: cs_drywell_arm_Pa        = 1.7_wp * PA_PER_PSI + P_ATM_PA
        ! LPCI auto-init thresholds: Level 1 reactor level + drywell hi-P
        ! (spec ch. 10.3 owns the trip logic; here we just expose a
        ! manual flag because phase-3 ECCS hasn't landed yet).
        logical  :: lpci_auto_init           = .false.
    end type rhr_config_t

    type :: rhr_drivers_t
        real(wp) :: dome_pressure_Pa     = 0.0_wp
        real(wp) :: reactor_coolant_T_K  = 293.15_wp
        real(wp) :: suppression_pool_T_K = 305.0_wp
        real(wp) :: drywell_pressure_Pa  = 101325.0_wp
        real(wp) :: reactor_level_m      = 1.0_wp
    end type rhr_drivers_t

    type :: rhr_command_t
        ! Sentinel −1 for "no change". 0..4 selects a mode; trip_loop_idx
        ! 1..4 trips a specific loop; reset_all_loops recovers from trips.
        integer :: mode_set            = -1
        integer :: trip_loop_idx       = 0
        logical :: reset_all_loops     = .false.
        ! Per-loop pump enable. -1 = no change. 0 = stop. 1 = start.
        integer :: pump_set(RHR_N_LOOPS) = -1
    end type rhr_command_t

    type :: rhr_observation_t
        real(wp) :: time              = 0.0_wp
        integer  :: n_steps           = 0
        integer  :: mode              = RHR_MODE_STANDBY
        integer  :: n_pumps_running   = 0
        real(wp) :: total_flow_kg_s   = 0.0_wp
        real(wp) :: total_heat_W      = 0.0_wp   ! removed from coolant or pool
        real(wp) :: hx_inlet_T_K      = 0.0_wp
        real(wp) :: hx_outlet_T_K     = 0.0_wp
        logical  :: sdc_permissive    = .false.  ! true when SDC alignment allowed
        logical  :: cs_armed          = .false.
    end type rhr_observation_t

    type :: rhr_state_t
        real(wp) :: time     = 0.0_wp
        integer  :: n_steps  = 0
        integer  :: mode     = RHR_MODE_STANDBY
        logical  :: pump_running(RHR_N_LOOPS) = .false.
        logical  :: pump_tripped(RHR_N_LOOPS) = .false.
        ! Per-tick cached observables (also returned by rhr_observe).
        real(wp) :: total_heat_W      = 0.0_wp
        real(wp) :: total_flow_kg_s   = 0.0_wp
        real(wp) :: hx_inlet_T_K      = 0.0_wp
        real(wp) :: hx_outlet_T_K     = 0.0_wp
    end type rhr_state_t

contains

    subroutine rhr_init(state, config)
        type(rhr_state_t),  intent(out) :: state
        type(rhr_config_t), intent(in)  :: config
        ! config is intent(in) only to satisfy the contract; no-op for now.
        if (config%UA_W_per_K_per_loop < 0.0_wp) return
        state%time           = 0.0_wp
        state%n_steps        = 0
        state%mode           = RHR_MODE_STANDBY
        state%pump_running   = .false.
        state%pump_tripped   = .false.
        state%total_heat_W   = 0.0_wp
        state%total_flow_kg_s = 0.0_wp
        state%hx_inlet_T_K   = config%service_water_T_K
        state%hx_outlet_T_K  = config%service_water_T_K
    end subroutine rhr_init

    subroutine rhr_step(state, dt, drivers, config)
        type(rhr_state_t),   intent(inout) :: state
        real(wp),            intent(in)    :: dt
        type(rhr_drivers_t), intent(in)    :: drivers
        type(rhr_config_t),  intent(in)    :: config

        integer  :: i, n_run
        real(wp) :: T_in, Q_per_loop
        real(wp), parameter :: CP_WATER = 4186.0_wp   ! J/kg·K

        ! ── SDC permissive: kick mode back to STANDBY if dome P > limit.
        if (state%mode == RHR_MODE_SHUTDOWN_COOLING .and. &
            drivers%dome_pressure_Pa > config%sdc_max_pressure_Pa) then
            state%mode = RHR_MODE_STANDBY
        end if

        ! ── Per-loop pump state (no auto-start logic in standby for now;
        !    operator must hit pump_set + mode_set). Tripped pumps stay
        !    tripped until reset.
        do i = 1, RHR_N_LOOPS
            if (state%pump_tripped(i)) state%pump_running(i) = .false.
        end do

        n_run = count(state%pump_running)
        state%total_flow_kg_s = real(n_run, wp) * config%rated_flow_kg_s_per_loop

        ! ── Per-mode heat-exchanger inlet temperature. ────────────────────
        select case (state%mode)
        case (RHR_MODE_SHUTDOWN_COOLING)
            T_in = drivers%reactor_coolant_T_K
        case (RHR_MODE_SUPP_POOL_COOLING, RHR_MODE_CONTAINMENT_SPRAY, RHR_MODE_LPCI)
            T_in = drivers%suppression_pool_T_K
        case default
            T_in = config%service_water_T_K
        end select
        state%hx_inlet_T_K = T_in

        ! ── Heat removal per running loop. ────────────────────────────────
        ! Q = UA · (T_in − T_service). Outlet T comes from energy balance:
        ! Q = ṁ · cp · (T_in − T_out)
        if (n_run > 0 .and. state%mode /= RHR_MODE_STANDBY) then
            Q_per_loop = config%UA_W_per_K_per_loop * (T_in - config%service_water_T_K)
            Q_per_loop = max(0.0_wp, Q_per_loop)
            state%total_heat_W = real(n_run, wp) * Q_per_loop
            state%hx_outlet_T_K = T_in - state%total_heat_W &
                                       / max(1.0_wp, state%total_flow_kg_s * CP_WATER)
            state%hx_outlet_T_K = max(config%service_water_T_K, state%hx_outlet_T_K)
        else
            state%total_heat_W  = 0.0_wp
            state%hx_outlet_T_K = T_in
        end if

        state%time    = state%time + dt
        state%n_steps = state%n_steps + 1
    end subroutine rhr_step

    subroutine rhr_apply(state, command, config, drivers)
        type(rhr_state_t),   intent(inout) :: state
        type(rhr_command_t), intent(in)    :: command
        type(rhr_config_t),  intent(in)    :: config
        type(rhr_drivers_t), intent(in)    :: drivers

        integer :: i

        if (command%mode_set >= RHR_MODE_STANDBY .and. command%mode_set <= RHR_MODE_LPCI) then
            ! Permissive checks reject the mode set silently. Frontend can
            ! poll observation to confirm the mode actually changed.
            select case (command%mode_set)
            case (RHR_MODE_SHUTDOWN_COOLING)
                if (drivers%dome_pressure_Pa <= config%sdc_max_pressure_Pa) &
                    state%mode = command%mode_set
            case (RHR_MODE_CONTAINMENT_SPRAY)
                if (drivers%drywell_pressure_Pa >= config%cs_drywell_arm_Pa) &
                    state%mode = command%mode_set
            case default
                state%mode = command%mode_set
            end select
        end if

        if (command%trip_loop_idx >= 1 .and. command%trip_loop_idx <= RHR_N_LOOPS) &
            state%pump_tripped(command%trip_loop_idx) = .true.

        if (command%reset_all_loops) then
            state%pump_tripped = .false.
            state%pump_running = .false.
        end if

        do i = 1, RHR_N_LOOPS
            if (command%pump_set(i) == 0) then
                state%pump_running(i) = .false.
            else if (command%pump_set(i) == 1 .and. .not. state%pump_tripped(i)) then
                state%pump_running(i) = .true.
            end if
        end do
    end subroutine rhr_apply

    function rhr_observe(state, config, drivers) result(obs)
        type(rhr_state_t),   intent(in) :: state
        type(rhr_config_t),  intent(in) :: config
        type(rhr_drivers_t), intent(in) :: drivers
        type(rhr_observation_t) :: obs

        obs%time            = state%time
        obs%n_steps         = state%n_steps
        obs%mode            = state%mode
        obs%n_pumps_running = count(state%pump_running)
        obs%total_flow_kg_s = state%total_flow_kg_s
        obs%total_heat_W    = state%total_heat_W
        obs%hx_inlet_T_K    = state%hx_inlet_T_K
        obs%hx_outlet_T_K   = state%hx_outlet_T_K
        obs%sdc_permissive  = drivers%dome_pressure_Pa <= config%sdc_max_pressure_Pa
        obs%cs_armed        = drivers%drywell_pressure_Pa >= config%cs_drywell_arm_Pa
    end function rhr_observe

    subroutine rhr_destroy(state)
        type(rhr_state_t), intent(inout) :: state
        state%time           = 0.0_wp
        state%n_steps        = 0
        state%mode           = RHR_MODE_STANDBY
        state%pump_running   = .false.
        state%pump_tripped   = .false.
        state%total_heat_W   = 0.0_wp
        state%total_flow_kg_s = 0.0_wp
        state%hx_inlet_T_K   = 0.0_wp
        state%hx_outlet_T_K  = 0.0_wp
    end subroutine rhr_destroy

end module rhr