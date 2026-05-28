! Feedwater subsystem.
!
! State ownership transfers from bwr_state_t / vessel:
!   bwr_state_t%feedwater_flow_pct → feedwater%demand_flow_pct
!   vessel%feedwater_temp          → feedwater%feedwater_temp_K
!
! Convention: `flow_pct` is fraction of rated steam flow (100 % at rated
! power, > 100 % = above-rated FW makeup). Rated FW mass flow tracks
! rated steam flow at full power so vessel inventory stays balanced.
!
! Pump topology (spec 2.6 p.5-14):
!   Hotwell → 2 condensate pumps (3-stage centrifugal, 12,153 gpm @ 200
!   psig each) → SJAE/SPE/demineralisers → 2 booster pumps (12,153 gpm
!   @ 600 psig each, low-suction trip @ 35 psig with 20/45 s delays) →
!   LP heaters → 2 RFPs (14,000 gpm @ 1130 psig each, 67 % capacity, low-
!   suction trip @ 250 psig with 8/25 s delays) → HP heaters → 4 vessel
!   spargers.
!
! Pump-chain suction propagation: discharge pressure ≈ rated × running/2.
! When upstream stage discharge drops below the next stage's trip
! setpoint, the trip timer starts; on expiry the pump latches tripped.
!
! 3-element control (spec 2.6 p.211, controller owned by chapter 3.3):
!   demand_pct = 100 + K_L · (L_ref - L_current)
!                    + K_S · (steam_flow_pct - fw_flow_pct)
! The controller is OFF by default (manual operator demand); the C ABI
! `bwr_set_feedwater_flow_pct` writes the operator demand directly.
!
! Heater-train outlet temperature: simple linear model — at zero flow
! the train is at cold-inlet temperature (~46 °C); at rated flow it
! reaches `rated_fw_inlet_temp_K` (488.7 K = 420 °F per spec 2.1).
! The fw_inlet_T discrepancy between spec 2.1 (420 °F) and spec 2.6
! (300 °F) is resolved by `feedwater_config_t%rated_inlet_temperature_K`:
! single source of truth, default 488.7 K (BWR/4 Browns Ferry rating).
module feedwater
    use kinds, only: wp
    implicit none
    private

    public :: feedwater_state_t, feedwater_config_t
    public :: feedwater_command_t, feedwater_observation_t
    public :: feedwater_drivers_t
    public :: feedwater_init, feedwater_destroy
    public :: feedwater_step, feedwater_observe
    public :: feedwater_apply
    public :: FW_N_CONDENSATE, FW_N_BOOSTER, FW_N_RFP

    integer, parameter :: FW_N_CONDENSATE = 2
    integer, parameter :: FW_N_BOOSTER    = 2
    integer, parameter :: FW_N_RFP        = 2

    real(wp), parameter :: PA_PER_PSI = 6894.76_wp
    real(wp), parameter :: P_ATM_PA   = 101325.0_wp

    type :: feedwater_config_t
        ! Rated FW mass flow ≈ rated steam flow at full power. The
        ! 10.5 × 10⁶ lb/hr value from spec 2.1 = 1323 kg/s; the spec 2.5
        ! ~13 × 10⁶ lb/hr quote gives 1638 kg/s. Pick the 2.5 value to
        ! match main_steam (single source-of-truth).
        real(wp) :: rated_flow_kg_s          = 1638.0_wp
        ! Heater-train hot outlet target. Two specs quote different
        ! numbers (2.1: 420 °F = 488.7 K; 2.6: 300 °F = 422.0 K). 488.7 K
        ! reflects plant-rated operation; switch via this single config
        ! to flip to the 2.6 training-deck value.
        real(wp) :: rated_inlet_temperature_K = 488.7_wp
        ! Condenser hotwell water temperature (cold-suction temperature).
        real(wp) :: cold_inlet_temperature_K  = 319.0_wp   ! ~46 °C / 115 °F
        ! Reference vessel level for 3-element control. Spec 3.1 normal
        ! water level ≈ +37 in above top of active fuel = 0.94 m.
        real(wp) :: level_setpoint_m         = 1.0_wp
        ! 3-element control gains (outer-loop feedforward form). The
        ! spec's (ṁ_steam − ṁ_FW) mismatch is the inner-loop RFPT-speed
        ! trim; here flow tracks RFPT speed directly so the explicit FW
        ! feedback is folded into the slewed pump speed.
        ! When chapter 3.3 (FWC) lands these gains get pinned to FSAR.
        real(wp) :: K_level             = 30.0_wp ! pct / m of level error
        real(wp) :: K_steam_feedforward = 1.0_wp  ! demand / steam-norm
        logical  :: controller_enabled  = .true.
        ! Auto-engage threshold. Once dome pressure climbs above this
        ! (i.e., the reactor is boiling and producing real steam flow),
        ! the controller flips on if it wasn't already, so the operator
        ! doesn't have to remember to enable it manually during startup.
        ! Cold-preset boot disables explicitly and zeros the demand; the
        ! auto-engage only fires after pressure rises past the threshold.
        logical  :: auto_engage_on_boil = .true.
        real(wp) :: auto_engage_pressure_Pa = 2.0e6_wp  ! ~290 psig
        ! Manual operator demand persists across the auto-engage flip:
        ! when controller is OFF, demand_flow_pct is taken as a manual
        ! setpoint; when ON, the controller overwrites it each tick.
        ! RFPT speed dynamics — turbine-driven pump speed lag.
        real(wp) :: rfpt_speed_slew_pct_s = 20.0_wp
        ! Low-suction trip setpoints + delays (spec 2.6 p.8 / p.11).
        real(wp) :: booster_trip_setpoint_psig = 35.0_wp
        real(wp) :: rfp_trip_setpoint_psig     = 250.0_wp
        real(wp) :: booster_trip_delay_s(FW_N_BOOSTER) = [20.0_wp, 45.0_wp]
        real(wp) :: rfp_trip_delay_s(FW_N_RFP)         = [ 8.0_wp, 25.0_wp]
        ! Rated discharge pressures (used for suction propagation).
        real(wp) :: condensate_discharge_psig = 200.0_wp
        real(wp) :: booster_discharge_psig    = 600.0_wp
        real(wp) :: rfp_discharge_psig        = 1130.0_wp
    end type feedwater_config_t

    type :: feedwater_drivers_t
        ! Reactor vessel water level [m above TAF] — from vessel obs.
        real(wp) :: reactor_level_m   = 1.0_wp
        ! Steam flow [normalised to rated, 0..1+] — from main_steam obs.
        real(wp) :: steam_flow_norm   = 1.0_wp
        ! Dome pressure [Pa] — gates the auto-engage flip below.
        real(wp) :: dome_pressure_Pa  = 0.0_wp
    end type feedwater_drivers_t

    type :: feedwater_command_t
        ! Sentinel −1.0 means "no change" for operator demand. Positive
        ! values are interpreted as a flow % command; the controller
        ! flag, if true, ignores this and runs 3-element control.
        real(wp) :: demand_flow_pct_set = -1.0_wp
        ! Trip controls. trip_*_idx ∈ [1, n], 0 means "no action".
        integer  :: trip_condensate_idx = 0
        integer  :: trip_booster_idx    = 0
        integer  :: trip_rfp_idx        = 0
        ! Reset all tripped pumps + restart all chains.
        logical  :: reset_all_pumps     = .false.
        ! Toggle 3-element control on/off (-1 = no change, 0 = off, 1 = on).
        integer  :: controller_enable_set = -1
    end type feedwater_command_t

    type :: feedwater_observation_t
        real(wp) :: time                = 0.0_wp
        integer  :: n_steps             = 0
        real(wp) :: flow_pct            = 0.0_wp
        real(wp) :: flow_kg_s           = 0.0_wp
        real(wp) :: feedwater_temp_K    = 0.0_wp
        real(wp) :: rfpt_speed_pct      = 0.0_wp
        real(wp) :: demand_pct          = 0.0_wp
        integer  :: condensate_running  = 0
        integer  :: booster_running     = 0
        integer  :: rfp_running         = 0
        real(wp) :: booster_suction_psig = 0.0_wp
        real(wp) :: rfp_suction_psig     = 0.0_wp
    end type feedwater_observation_t

    type :: feedwater_state_t
        real(wp) :: time     = 0.0_wp
        integer  :: n_steps  = 0

        ! Operator demand + actual flow.
        real(wp) :: demand_flow_pct = 0.0_wp
        real(wp) :: flow_pct        = 0.0_wp     ! after pump-chain + slew
        real(wp) :: flow_kg_s       = 0.0_wp
        real(wp) :: feedwater_temp_K = 0.0_wp

        ! RFPT speed (single shared scalar — both RFPTs follow same demand).
        real(wp) :: rfpt_speed_pct  = 0.0_wp

        ! Pump-chain state.
        logical  :: condensate_tripped(FW_N_CONDENSATE) = .false.
        logical  :: booster_tripped(FW_N_BOOSTER)       = .false.
        logical  :: rfp_tripped(FW_N_RFP)               = .false.

        ! Suction-trip latch timers (seconds of continuous low-suction).
        real(wp) :: booster_low_suction_t(FW_N_BOOSTER) = 0.0_wp
        real(wp) :: rfp_low_suction_t(FW_N_RFP)         = 0.0_wp

        logical  :: controller_enabled = .false.
    end type feedwater_state_t

contains

    pure function n_running(tripped) result(n)
        logical, intent(in) :: tripped(:)
        integer :: n
        n = count(.not. tripped)
    end function n_running

    pure function clamp01(x, lo, hi) result(y)
        real(wp), intent(in) :: x, lo, hi
        real(wp) :: y
        y = max(lo, min(hi, x))
    end function clamp01

    subroutine feedwater_init(state, config)
        type(feedwater_state_t),  intent(out) :: state
        type(feedwater_config_t), intent(in)  :: config

        state%time            = 0.0_wp
        state%n_steps         = 0

        state%demand_flow_pct = 100.0_wp
        state%flow_pct        = 100.0_wp
        state%flow_kg_s       = 1.0_wp * config%rated_flow_kg_s
        state%feedwater_temp_K = config%rated_inlet_temperature_K
        state%rfpt_speed_pct  = 100.0_wp

        state%condensate_tripped = .false.
        state%booster_tripped    = .false.
        state%rfp_tripped        = .false.
        state%booster_low_suction_t = 0.0_wp
        state%rfp_low_suction_t     = 0.0_wp
        state%controller_enabled = config%controller_enabled
    end subroutine feedwater_init

    subroutine feedwater_step(state, dt, drivers, config)
        type(feedwater_state_t),   intent(inout) :: state
        real(wp),                  intent(in)    :: dt
        type(feedwater_drivers_t), intent(in)    :: drivers
        type(feedwater_config_t),  intent(in)    :: config

        integer  :: i, n_cond, n_boost, n_rfp
        real(wp) :: cond_disch_psig, boost_disch_psig, rfp_disch_psig
        real(wp) :: boost_suction_psig, rfp_suction_psig
        real(wp) :: demand, target_rfpt, max_step
        real(wp) :: pump_capacity_factor
        real(wp) :: t_demanded

        ! ── Pump-chain suction pressure propagation. ─────────────────────
        ! Discharge pressure scales linearly with the running pump count
        ! (each pump is rated for half of system capacity at design speed).
        n_cond  = n_running(state%condensate_tripped)
        n_boost = n_running(state%booster_tripped)
        n_rfp   = n_running(state%rfp_tripped)

        cond_disch_psig  = config%condensate_discharge_psig &
                         * real(n_cond, wp) / real(FW_N_CONDENSATE, wp)
        boost_disch_psig = config%booster_discharge_psig &
                         * real(n_boost, wp) / real(FW_N_BOOSTER, wp)

        ! Suction = upstream discharge minus header/heater ΔP at full flow.
        ! Approximate with a flat ~50 psi heater drop in each direction.
        boost_suction_psig = max(0.0_wp, cond_disch_psig - 50.0_wp)
        rfp_suction_psig   = max(0.0_wp, boost_disch_psig - 50.0_wp)

        ! ── Booster low-suction trip cascade. ────────────────────────────
        do i = 1, FW_N_BOOSTER
            if (state%booster_tripped(i)) cycle
            if (boost_suction_psig < config%booster_trip_setpoint_psig) then
                state%booster_low_suction_t(i) = state%booster_low_suction_t(i) + dt
                if (state%booster_low_suction_t(i) >= config%booster_trip_delay_s(i)) then
                    state%booster_tripped(i) = .true.
                end if
            else
                state%booster_low_suction_t(i) = 0.0_wp
            end if
        end do

        ! Re-evaluate downstream pressures after booster trips this tick.
        n_boost = n_running(state%booster_tripped)
        boost_disch_psig = config%booster_discharge_psig &
                         * real(n_boost, wp) / real(FW_N_BOOSTER, wp)
        rfp_suction_psig = max(0.0_wp, boost_disch_psig - 50.0_wp)

        ! ── RFP low-suction trip cascade. ────────────────────────────────
        do i = 1, FW_N_RFP
            if (state%rfp_tripped(i)) cycle
            if (rfp_suction_psig < config%rfp_trip_setpoint_psig) then
                state%rfp_low_suction_t(i) = state%rfp_low_suction_t(i) + dt
                if (state%rfp_low_suction_t(i) >= config%rfp_trip_delay_s(i)) then
                    state%rfp_tripped(i) = .true.
                end if
            else
                state%rfp_low_suction_t(i) = 0.0_wp
            end if
        end do
        n_rfp = n_running(state%rfp_tripped)

        rfp_disch_psig = config%rfp_discharge_psig &
                       * real(n_rfp, wp) / real(FW_N_RFP, wp)

        ! ── Auto-engage flip. ────────────────────────────────────────────
        ! Once dome pressure rises past the boiling-active threshold,
        ! flip the controller on. This converts the cold-preset "manual
        ! zero demand" into automatic ṁ_FW = ṁ_steam tracking once the
        ! plant is actually making steam. The flip is one-way per boot
        ! (a SCRAM that drops pressure does not disengage); to take back
        ! manual control the operator calls bwr_set_fw_controller_enabled(0).
        if (config%auto_engage_on_boil .and. .not. state%controller_enabled &
            .and. drivers%dome_pressure_Pa > config%auto_engage_pressure_Pa) then
            state%controller_enabled = .true.
        end if

        ! ── 3-element control (optional). ────────────────────────────────
        ! Outer loop: feedforward steam flow (matches boil-off rate at
        ! steady state) plus proportional level trim. The third element
        ! (FW flow) closes the inner loop by reading RFPT speed; the
        ! min(rfpt × cap, demand) clip below ties the two together
        ! without introducing the demand↔flow oscillation that the naive
        ! `K_S · (steam − FW)` form gives at large steam-flow excursions.
        if (state%controller_enabled) then
            demand = config%K_steam_feedforward * drivers%steam_flow_norm * 100.0_wp &
                   + config%K_level * (config%level_setpoint_m - drivers%reactor_level_m)
            demand = clamp01(demand, 0.0_wp, 200.0_wp)
            state%demand_flow_pct = demand
        end if

        ! ── RFPT speed slew toward the controller / operator demand. ─────
        target_rfpt = clamp01(state%demand_flow_pct, 0.0_wp, 100.0_wp)
        max_step = config%rfpt_speed_slew_pct_s * dt
        if (target_rfpt > state%rfpt_speed_pct) then
            state%rfpt_speed_pct = min(target_rfpt, state%rfpt_speed_pct + max_step)
        else
            state%rfpt_speed_pct = max(target_rfpt, state%rfpt_speed_pct - max_step)
        end if

        ! ── FW mass-flow: pump-running fraction × demanded RFPT speed. ───
        ! Each RFP is rated 67 % of system capacity (spec 2.6 p.11).
        ! Two pumps + 100 % RFPT speed = 100 % rated FW. One pump caps at
        ! 67 % regardless of demand; zero pumps = zero flow.
        if (n_rfp == 0) then
            pump_capacity_factor = 0.0_wp
        else if (n_rfp == 1) then
            pump_capacity_factor = 0.67_wp
        else
            pump_capacity_factor = 1.0_wp
        end if

        state%flow_pct = (state%rfpt_speed_pct / 100.0_wp) &
                       * pump_capacity_factor * 100.0_wp
        ! Operator demand above 100 % is only achievable if pump capacity
        ! allows; clip flow_pct to demand to preserve operator intent.
        state%flow_pct = min(state%flow_pct, state%demand_flow_pct)
        state%flow_kg_s = (state%flow_pct / 100.0_wp) * config%rated_flow_kg_s

        ! ── Heater-train outlet temperature. ─────────────────────────────
        ! Approximation: rated FW flow + extraction-steam heaters bring
        ! the train to rated_inlet_temperature_K. At very low flow, the
        ! heaters cool toward cold-inlet (no extraction steam). Linear
        ! ramp keyed off flow_pct; saturates at 30 % flow per spec p.15
        ! (RFP-startup handoff implies heaters effective at ≥ 30 % power).
        t_demanded = config%cold_inlet_temperature_K &
                   + (config%rated_inlet_temperature_K - config%cold_inlet_temperature_K) &
                     * clamp01(state%flow_pct / 30.0_wp, 0.0_wp, 1.0_wp)
        ! First-order lag (~5 s time constant) to soften step changes.
        state%feedwater_temp_K = state%feedwater_temp_K &
                              + (t_demanded - state%feedwater_temp_K) &
                                * clamp01(dt / 5.0_wp, 0.0_wp, 1.0_wp)

        state%time    = state%time + dt
        state%n_steps = state%n_steps + 1
    end subroutine feedwater_step

    subroutine feedwater_apply(state, command)
        type(feedwater_state_t),   intent(inout) :: state
        type(feedwater_command_t), intent(in)    :: command

        if (command%demand_flow_pct_set >= 0.0_wp) &
            state%demand_flow_pct = clamp01(command%demand_flow_pct_set, 0.0_wp, 200.0_wp)

        if (command%trip_condensate_idx >= 1 .and. &
            command%trip_condensate_idx <= FW_N_CONDENSATE) &
            state%condensate_tripped(command%trip_condensate_idx) = .true.

        if (command%trip_booster_idx >= 1 .and. &
            command%trip_booster_idx <= FW_N_BOOSTER) &
            state%booster_tripped(command%trip_booster_idx) = .true.

        if (command%trip_rfp_idx >= 1 .and. command%trip_rfp_idx <= FW_N_RFP) &
            state%rfp_tripped(command%trip_rfp_idx) = .true.

        if (command%reset_all_pumps) then
            state%condensate_tripped    = .false.
            state%booster_tripped       = .false.
            state%rfp_tripped           = .false.
            state%booster_low_suction_t = 0.0_wp
            state%rfp_low_suction_t     = 0.0_wp
        end if

        if (command%controller_enable_set == 0) then
            state%controller_enabled = .false.
        else if (command%controller_enable_set == 1) then
            state%controller_enabled = .true.
        end if
    end subroutine feedwater_apply

    function feedwater_observe(state, config) result(obs)
        type(feedwater_state_t),  intent(in) :: state
        type(feedwater_config_t), intent(in) :: config
        type(feedwater_observation_t) :: obs

        integer  :: n_cond, n_boost, n_rfp
        real(wp) :: cond_disch_psig, boost_disch_psig

        n_cond  = n_running(state%condensate_tripped)
        n_boost = n_running(state%booster_tripped)
        n_rfp   = n_running(state%rfp_tripped)

        cond_disch_psig  = config%condensate_discharge_psig &
                         * real(n_cond, wp) / real(FW_N_CONDENSATE, wp)
        boost_disch_psig = config%booster_discharge_psig &
                         * real(n_boost, wp) / real(FW_N_BOOSTER, wp)

        obs%time              = state%time
        obs%n_steps           = state%n_steps
        obs%flow_pct          = state%flow_pct
        obs%flow_kg_s         = state%flow_kg_s
        obs%feedwater_temp_K  = state%feedwater_temp_K
        obs%rfpt_speed_pct    = state%rfpt_speed_pct
        obs%demand_pct        = state%demand_flow_pct
        obs%condensate_running = n_cond
        obs%booster_running    = n_boost
        obs%rfp_running        = n_rfp
        obs%booster_suction_psig = max(0.0_wp, cond_disch_psig - 50.0_wp)
        obs%rfp_suction_psig     = max(0.0_wp, boost_disch_psig - 50.0_wp)
    end function feedwater_observe

    subroutine feedwater_destroy(state)
        type(feedwater_state_t), intent(inout) :: state
        state%time     = 0.0_wp
        state%n_steps  = 0
        state%demand_flow_pct = 0.0_wp
        state%flow_pct        = 0.0_wp
        state%flow_kg_s       = 0.0_wp
        state%feedwater_temp_K = 0.0_wp
        state%rfpt_speed_pct  = 0.0_wp
        state%condensate_tripped    = .false.
        state%booster_tripped       = .false.
        state%rfp_tripped           = .false.
        state%booster_low_suction_t = 0.0_wp
        state%rfp_low_suction_t     = 0.0_wp
        state%controller_enabled = .false.
    end subroutine feedwater_destroy

end module feedwater