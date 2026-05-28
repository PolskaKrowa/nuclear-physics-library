! Vessel subsystem.
!
! The dome-pressure + level dynamics that used to live inside
! `update_pressure_dynamics` now run as `vessel_step_dynamics`, operating
! on a typed `vessel_state_t` that owns its ODE history (P_prev1/2,
! coolant_T_prev1/2, T_excess, sat_temperature, Lbase, Lrx,
! pressure_operating).
!
! NOTE on magic coefficients: `dP_COEF`, `P_LOSS_COEF`, and `P_NOMINAL` are
! kept verbatim from the pre-reorg code for behavioural continuity. They
! don't actually balance at the rated operating point (the ODE drifts under
! rated drivers); replacing the model with a proper steam-mass balance is
! a phase-2 task — see AUDIT §3 B5.
module vessel
    use kinds, only: wp
    implicit none
    private

    public :: vessel_state_t, vessel_config_t
    public :: vessel_command_t, vessel_observation_t
    public :: vessel_drivers_t
    public :: vessel_init, vessel_destroy
    public :: vessel_step, vessel_observe
    public :: vessel_apply
    public :: vessel_step_dynamics
    public :: sat_temp_K, coolant_density_gl, calc_reactor_level

    ! Magic-coefficient ODE constants. See module-header note.
    real(wp), parameter :: dP_COEF     = 5000.0_wp
    real(wp), parameter :: P_LOSS_COEF = 1.125e6_wp
    ! Reference dome pressure used by the steam-loss term. Stays at the
    ! original 7.14 MPa for behavioural continuity even though
    ! `rated_pressure_Pa` defaults to 1000 psig (6.895 MPa) per spec 2.1.
    real(wp), parameter :: P_NOMINAL   = 7.14e6_wp

    ! Hard clamps that used to live inline in update_pressure_dynamics.
    ! P_MIN_PA dropped from 3.0e6 to 1 atm so the vessel model can
    ! represent cold standby (no boiling -> dome at atmospheric). The
    ! near-rated dynamics never approach this floor; lowering it does
    ! not affect at-power behaviour.
    real(wp), parameter :: P_MIN_PA   = 1.01325e5_wp
    real(wp), parameter :: P_MAX_PA   = 11.0e6_wp
    real(wp), parameter :: L_BASE_MIN = 0.5_wp
    real(wp), parameter :: L_BASE_HEADROOM = 3.0_wp

    type :: vessel_config_t
        ! Rated values from spec 2.1 (BWR/4 GE-14 class, Browns Ferry /
        ! Peach Bottom). Used by `vessel_init` to seed the ODE state.
        real(wp) :: rated_pressure_Pa     = 6.895e6_wp     ! ~1000 psig (spec 2.1, p.16)
        real(wp) :: rated_core_flow_kg_s  = 9700.0_wp      ! 77.0e6 lb/hr (spec 2.1, p.15)
        ! FW inlet temperature: spec 2.1 says ~420 °F at vessel entry,
        ! spec 2.6 says ~300 °F at the heater string outlet. Plan calls
        ! out this discrepancy; vessel takes the spec 2.1 value as the
        ! authoritative vessel-inlet temperature.
        real(wp) :: rated_fw_inlet_temp_K = 488.7_wp       ! 420 °F (spec 2.1, p.16)
        real(wp) :: core_height_m         = 3.81_wp        ! 12.5 ft active fuel
        ! Effective downcomer area used by the level mass balance: only
        ! ~half the 218" ID cross-section (≈24 m²) contributes to indicated
        ! level since the shroud+separators displace the centerline volume.
        real(wp) :: downcomer_area_m2     = 12.0_wp
    end type vessel_config_t

    type :: vessel_command_t
        ! Vessel has no operator commands in phase 1; the dome is a passive
        ! pressure boundary. SRV opens / MSIV closures are commanded into
        ! `main_steam` (step 5).
        integer :: reserved = 0
    end type vessel_command_t

    ! Magic relief coefficient: each unit of SRV "flow fraction"
    ! (srv_flow_kg_s / rated_steam_flow_kg_s) takes the same chunk out of
    ! the dome-pressure ODE as a fully-open turbine line at rated dome P.
    ! Same units as P_LOSS_COEF so the two terms add.
    real(wp), parameter :: SRV_LOSS_COEF = 1.125e6_wp

    ! Bypass capacity cap: hardware-limited to 25 % of rated steam flow
    ! at 100 % valve open (spec 2.5 p.11). The visible bypass_valve_pct
    ! is the actual stroke; the relief term applies this factor.
    real(wp), parameter :: BYPASS_CAPACITY_FRAC = 0.25_wp

    type :: vessel_drivers_t
        ! Upstream state read each tick. Bundled into one struct so
        ! `vessel_step_dynamics` has a stable signature as more subsystems
        ! come online.
        real(wp) :: avg_coolant_temp_K = 0.0_wp
        real(wp) :: power_current_W    = 0.0_wp
        real(wp) :: power_rated_W      = 1.0_wp
        real(wp) :: feedwater_flow_pct = 0.0_wp
        real(wp) :: turbine_valve_pct  = 0.0_wp
        real(wp) :: bypass_valve_pct   = 0.0_wp
        ! Mass-flow drivers for the level ODE. Both come from peer subsystem
        ! observations: ṁ_FW from feedwater (rated 1638 kg/s at 100 %),
        ! ṁ_steam from main_steam (rated 1638 kg/s at 100 % TCV+BPV at
        ! rated pressure). Their difference is the rate of change of vessel
        ! water inventory; vessel_step_dynamics converts to a level rate.
        real(wp) :: feedwater_flow_kg_s = 0.0_wp
        real(wp) :: steam_flow_kg_s     = 0.0_wp
        ! SRV relief as a fraction of rated steam flow (already scaled by
        ! main_steam_observe). Set by the orchestrator each tick.
        real(wp) :: srv_flow_frac      = 0.0_wp
        ! MSIV gating: fraction of MSL paths still open (0=all isolated).
        ! The turbine/bypass loss term is gated by this since steam can't
        ! reach the equalising header through a closed MSIV.
        real(wp) :: msiv_open_frac     = 1.0_wp
        logical  :: turbine_tripped    = .false.
        real(wp) :: core_height_m      = 0.0_wp
    end type vessel_drivers_t

    type :: vessel_observation_t
        real(wp) :: time             = 0.0_wp
        integer  :: n_steps          = 0
        real(wp) :: pressure_Pa      = 0.0_wp
        real(wp) :: sat_temperature  = 0.0_wp
        real(wp) :: reactor_level_m  = 0.0_wp
        real(wp) :: T_excess         = 0.0_wp
    end type vessel_observation_t

    type :: vessel_state_t
        real(wp) :: time         = 0.0_wp
        integer  :: n_steps      = 0

        ! Dome-pressure ODE state
        real(wp) :: pressure_operating = 0.0_wp
        real(wp) :: sat_temperature    = 0.0_wp
        real(wp) :: T_excess           = 0.0_wp
        real(wp) :: P_prev1            = 0.0_wp
        real(wp) :: P_prev2            = 0.0_wp

        ! Coolant temperature smoothing buffer
        real(wp) :: coolant_T_prev1 = 0.0_wp
        real(wp) :: coolant_T_prev2 = 0.0_wp

        ! Level ODE state
        real(wp) :: Lbase = 0.0_wp
        real(wp) :: Lrx   = 0.0_wp
    end type vessel_state_t

contains

    ! ── Pure helpers (moved verbatim from bwr_c_interface.f90) ───────────

    pure function sat_temp_K(P_Pa) result(T_K)
        real(wp), intent(in) :: P_Pa
        real(wp) :: T_K, p
        p   = P_Pa / 1.0e6_wp
        T_K = (-4.0964e-3_wp * p**6  &
              +  0.141738_wp  * p**5  &
              -  1.943057_wp  * p**4  &
              + 13.50428_wp   * p**3  &
              - 51.27379_wp   * p**2  &
              + 118.6854_wp   * p     &
              + 100.1542_wp) + 273.15_wp
    end function sat_temp_K

    pure function coolant_density_gl(T_K) result(rho)
        real(wp), intent(in) :: T_K
        real(wp) :: rho, T_C
        T_C = T_K - 273.15_wp
        rho = -4.467711e-6_wp    * T_C**3  &
              - 5.60288485e-4_wp  * T_C**2  &
              - 0.429148844451_wp * T_C     &
              + 1010.035413387815_wp
        rho = max(0.1_wp, rho)
    end function coolant_density_gl

    pure function calc_reactor_level(Lbase, T_K, H_core) result(Lrx)
        real(wp), intent(in) :: Lbase, T_K, H_core
        real(wp) :: Lrx, rho
        rho = coolant_density_gl(T_K)
        Lrx = (Lbase / (rho / 1000.0_wp)) - H_core
    end function calc_reactor_level

    ! ── 5-call contract ──────────────────────────────────────────────────

    subroutine vessel_init(state, config)
        type(vessel_state_t),  intent(out) :: state
        type(vessel_config_t), intent(in)  :: config

        state%time    = 0.0_wp
        state%n_steps = 0

        ! Seed the ODE history at the rated design point.
        state%pressure_operating = config%rated_pressure_Pa
        state%P_prev1            = config%rated_pressure_Pa
        state%P_prev2            = config%rated_pressure_Pa
        state%sat_temperature    = sat_temp_K(config%rated_pressure_Pa)
        state%T_excess           = 0.0_wp

        state%coolant_T_prev1 = config%rated_fw_inlet_temp_K
        state%coolant_T_prev2 = config%rated_fw_inlet_temp_K

        ! Lbase is a mass-equivalent level scalar; the original init formula
        ! (`0.5 + core_height` corrected by density) is retained verbatim.
        state%Lbase = (0.5_wp + config%core_height_m) &
                    * (coolant_density_gl(config%rated_fw_inlet_temp_K) / 1000.0_wp)
        state%Lrx   = calc_reactor_level(state%Lbase, &
                                         config%rated_fw_inlet_temp_K, &
                                         config%core_height_m)
    end subroutine vessel_init

    subroutine vessel_step(state, dt)
        ! Bare contract entry point: bumps tick counters. Real physics is
        ! `vessel_step_dynamics` which needs upstream drivers; the
        ! orchestrator (`bwr_c_interface`) calls that directly.
        type(vessel_state_t), intent(inout) :: state
        real(wp),             intent(in)    :: dt
        state%time    = state%time + dt
        state%n_steps = state%n_steps + 1
    end subroutine vessel_step

    subroutine vessel_step_dynamics(state, dt, drivers, config)
        type(vessel_state_t),   intent(inout) :: state
        real(wp),               intent(in)    :: dt
        type(vessel_drivers_t), intent(in)    :: drivers
        type(vessel_config_t),  intent(in)    :: config

        real(wp) :: T_sat, T_mean_smooth, rel_T, dP_gen, steam_val, dP_loss
        real(wp) :: net_kg_s, dLbase_dt, area_m2
        real(wp), parameter :: RHO_REF = 1000.0_wp  ! reference density used by calc_reactor_level

        ! Coolant temperature smoothing (3-tap rolling average).
        state%coolant_T_prev2 = state%coolant_T_prev1
        state%coolant_T_prev1 = drivers%avg_coolant_temp_K
        T_mean_smooth = (drivers%avg_coolant_temp_K &
                        + state%coolant_T_prev1 + state%coolant_T_prev2) / 3.0_wp

        T_sat = sat_temp_K(state%pressure_operating)
        state%sat_temperature = T_sat
        state%T_excess = max(0.0_wp, T_mean_smooth - T_sat)

        ! Dome-pressure generation term (magic coefficients).
        rel_T  = (T_mean_smooth - 373.15_wp) &
               * (1.0_wp - state%pressure_operating / 10.0e6_wp)
        dP_gen = dP_COEF * (rel_T + 30.0_wp * state%T_excess) * dt

        ! Dome-pressure loss to turbine + bypass valves. Bypass valve %
        ! is the raw stroke position; the 25 % hardware capacity cap
        ! (spec 2.5 p.11) is applied here. Both flow paths are gated by
        ! the MSIV open fraction since they sit downstream of the MSIVs.
        if (.not. drivers%turbine_tripped) then
            steam_val = (drivers%turbine_valve_pct &
                         + BYPASS_CAPACITY_FRAC * drivers%bypass_valve_pct) &
                        / 100.0_wp &
                        * drivers%msiv_open_frac &
                        * (state%pressure_operating / P_NOMINAL)
        else
            steam_val = BYPASS_CAPACITY_FRAC * (drivers%bypass_valve_pct / 100.0_wp) &
                        * drivers%msiv_open_frac &
                        * (state%pressure_operating / P_NOMINAL)
        end if
        dP_loss = -P_LOSS_COEF * steam_val * dt &
                  - SRV_LOSS_COEF * drivers%srv_flow_frac * dt

        state%P_prev2 = state%P_prev1
        state%P_prev1 = state%pressure_operating
        state%pressure_operating = max(P_MIN_PA, min(P_MAX_PA, &
            state%pressure_operating + dP_gen + dP_loss))

        ! Level dynamics: derived mass balance over vessel inventory.
        ! dLbase/dt = (ṁ_FW − ṁ_steam) / (A_eff · ρ_ref)
        ! where Lbase is the mass-equivalent height tracked by the level
        ! observer (calc_reactor_level then converts to actual level using
        ! coolant_density_gl(T)). At rated balance (1638 kg/s each side)
        ! the rate is zero. Decoupling from the 0.002 m/s placeholder lets
        ! the operator overdrive / underdrive feedwater and see the right
        ! level slope (e.g., FW pump trip drops level at ≈ -0.14 m/s while
        ! rated boil-off continues, matching plant Tech Specs).
        area_m2 = max(1.0_wp, config%downcomer_area_m2)
        net_kg_s = drivers%feedwater_flow_kg_s - drivers%steam_flow_kg_s
        dLbase_dt = net_kg_s / (area_m2 * RHO_REF)
        state%Lbase = state%Lbase + dLbase_dt * dt
        state%Lbase = max(L_BASE_MIN, min(drivers%core_height_m + L_BASE_HEADROOM, state%Lbase))

        state%Lrx = calc_reactor_level(state%Lbase, drivers%avg_coolant_temp_K, &
                                       drivers%core_height_m)

        state%time    = state%time + dt
        state%n_steps = state%n_steps + 1
    end subroutine vessel_step_dynamics

    subroutine vessel_apply(state, command)
        type(vessel_state_t),   intent(inout) :: state
        type(vessel_command_t), intent(in)    :: command
        ! No vessel-level operator commands in phase 1. Touch both args to
        ! satisfy -Wextra without introducing fake side-effects.
        if (command%reserved /= 0) state%n_steps = state%n_steps
    end subroutine vessel_apply

    function vessel_observe(state) result(obs)
        type(vessel_state_t), intent(in) :: state
        type(vessel_observation_t) :: obs
        obs%time             = state%time
        obs%n_steps          = state%n_steps
        obs%pressure_Pa      = state%pressure_operating
        obs%sat_temperature  = state%sat_temperature
        obs%reactor_level_m  = state%Lrx
        obs%T_excess         = state%T_excess
    end function vessel_observe

    subroutine vessel_destroy(state)
        type(vessel_state_t), intent(inout) :: state
        state%time              = 0.0_wp
        state%n_steps           = 0
        state%pressure_operating = 0.0_wp
        state%sat_temperature   = 0.0_wp
        state%T_excess          = 0.0_wp
        state%P_prev1           = 0.0_wp
        state%P_prev2           = 0.0_wp
        state%coolant_T_prev1   = 0.0_wp
        state%coolant_T_prev2   = 0.0_wp
        state%Lbase             = 0.0_wp
        state%Lrx               = 0.0_wp
    end subroutine vessel_destroy

end module vessel