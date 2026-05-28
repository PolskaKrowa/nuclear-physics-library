! Main steam subsystem.
!
! State ownership transfers from bwr_state_t:
!   turbine_valve   → main_steam%turbine_valve_pct  (actual TSV position)
!   bypass_valve    → main_steam%bypass_valve_pct
!   turbine_speed   → main_steam%turbine_speed_rpm
!   turbine_power   → main_steam%turbine_power_W
!   turbine_tripped → main_steam%turbine_tripped
!
! Conventions:
!   *_pct fields are valve fractions in [0, 100]; 0 = closed, 100 = open.
!   msiv_pos_pct(4) tracks the four MSL inboard/outboard pairs combined
!   (we collapse the 2-MSIV-per-line pair into one stroke state since
!   both swing together on an NSSSS isolation signal).
!   srv_open(11) is the per-valve open/closed flag. SRV lift is the
!   pilot-operated safety mode (dome pressure ≥ setpoint).
!
! SRV layout (3 staggered safety setpoint groups + ADS subset):
!   Indices 1-4 : group 1, setpoint 7.72e6 Pa (~1105 psig). 4 valves.
!   Indices 5-8 : group 2, setpoint 7.79e6 Pa (~1115 psig). 4 valves.
!   Indices 9-11: group 3, setpoint 7.93e6 Pa (~1135 psig). 3 valves.
!   ADS-designated: indices 1-7 (the 7 lowest-setpoint valves).
! Per-valve capacity 106 kg/s (~840 klbm/hr at lift) — typical BWR/4.
! Numbers are flagged `[BWR/4 typical]` in spec 05 since the deck is silent.
!
! TSV fast closure (p.12 of spec 05): full close in 0.1 s on turbine trip.
! In normal valve motion the operator setpoint slews at ~100 %/s.
!
! MSIV stroke (spec 05 "[BWR/4 typical 3-5 seconds full closure]"):
! 4.0 s default; configurable.
!
! Bypass capacity (p.11): 25 % of rated steam flow even at 100 % valve
! position. The valve % output from this subsystem is the raw stroke
! position; vessel.f90 applies the 0.25 capacity factor when computing
! its pressure-loss term.
!
! Turbine / generator routing (rework):
!   Steam leaves the equalising header through two parallel paths.
!     - TCV path: drives the turbine. Only this flow imparts torque on
!       the shaft, so only TCV flow accelerates the rotor.
!     - BPV path: dumps to the condenser. Drains the dome but does NOT
!       drive the turbine (no acceleration from bypass).
!   The shaft is a free spinner below sync (no electrical load) and a
!   speed-locked grid-tied machine once the operator closes the generator
!   breaker. Real BWR/4 tandem-compound 4-flow turbines run at 1800 rpm
!   into a 4-pole 60 Hz generator. Manual sync window: ±15 rpm around
!   1800. Below sync, all electrical output is zero; once synced, RPM is
!   pinned and turbine_power_W reflects mechanical power × generator η.
!   Load reject (open breaker) drops back to free-spinner with whatever
!   steam is still admitted, which can overspeed unless the operator
!   tracks TCV down (or governor / fast-close intervenes).
module main_steam
    use kinds, only: wp
    implicit none
    private

    public :: main_steam_state_t, main_steam_config_t
    public :: main_steam_command_t, main_steam_observation_t
    public :: main_steam_drivers_t
    public :: main_steam_init, main_steam_destroy
    public :: main_steam_step, main_steam_observe
    public :: main_steam_apply
    public :: MS_N_SRVS, MS_N_MSLS, MS_N_ADS

    integer, parameter :: MS_N_SRVS = 11
    integer, parameter :: MS_N_MSLS = 4
    integer, parameter :: MS_N_ADS  = 7

    real(wp), parameter :: P_ATM_PA  = 101325.0_wp
    real(wp), parameter :: PA_PER_PSI = 6894.76_wp

    ! Default SRV group setpoints (gauge pressure → absolute Pa).
    !   Group 1: 1105 psig → 7.72e6 Pa
    !   Group 2: 1115 psig → 7.79e6 Pa
    !   Group 3: 1135 psig → 7.93e6 Pa
    real(wp), parameter :: SRV_SP_G1 = 1105.0_wp * PA_PER_PSI + P_ATM_PA
    real(wp), parameter :: SRV_SP_G2 = 1115.0_wp * PA_PER_PSI + P_ATM_PA
    real(wp), parameter :: SRV_SP_G3 = 1135.0_wp * PA_PER_PSI + P_ATM_PA

    type :: main_steam_config_t
        real(wp) :: rated_steam_flow_kg_s    = 1638.0_wp   ! ~13e6 lbm/hr at 100 % (spec 2.5 open Q)
        real(wp) :: rated_pressure_Pa        = 7.14e6_wp   ! dome P_nominal (matches vessel)
        real(wp) :: msiv_stroke_time_s       = 4.0_wp      ! 3-5 s typical (spec 2.5 p.20 open Q)
        real(wp) :: tsv_normal_slew_pct_s    = 100.0_wp    ! demand-following rate
        real(wp) :: tsv_fast_close_time_s    = 0.1_wp      ! p.12 "fast closure"
        real(wp) :: bypass_slew_pct_s        = 50.0_wp     ! BPV stroke rate (no explicit deck value)
        real(wp) :: srv_capacity_kg_s        = 106.0_wp    ! ~840 klbm/hr per valve (spec open Q)
        real(wp) :: srv_reset_deadband_psi   = 50.0_wp     ! lift→reset deadband (typical SRV hysteresis)
        real(wp) :: srv_setpoint_Pa(MS_N_SRVS) = &
            [SRV_SP_G1, SRV_SP_G1, SRV_SP_G1, SRV_SP_G1, &      ! group 1 (× 4)
             SRV_SP_G2, SRV_SP_G2, SRV_SP_G2, SRV_SP_G2, &      ! group 2 (× 4)
             SRV_SP_G3, SRV_SP_G3, SRV_SP_G3]                   ! group 3 (× 3)
        logical  :: srv_is_ads(MS_N_SRVS) = &
            [.true., .true., .true., .true., &                  ! indices 1-7 are ADS
             .true., .true., .true., &
             .false., .false., .false., .false.]
        real(wp) :: initial_turbine_valve_pct = 33.0_wp
        real(wp) :: initial_bypass_valve_pct  = 0.0_wp

        ! Turbine / generator electromechanical model.
        real(wp) :: rated_speed_rpm        = 1800.0_wp  ! tandem-compound 4-flow @ 60 Hz, 4-pole gen
        real(wp) :: sync_window_rpm        = 15.0_wp    ! ±RPM tolerance for breaker closure
        real(wp) :: free_speed_overshoot   = 0.10_wp    ! +10 % above rated at full TCV with no load
        real(wp) :: speed_lag_time_s       = 8.0_wp     ! first-order RPM response time constant
        real(wp) :: speed_decel_rpm_s      = 60.0_wp    ! coast-down rate after trip / no steam
        real(wp) :: eta_turbine            = 0.34_wp    ! mechanical conversion (≈ Carnot fraction)
        real(wp) :: eta_generator          = 0.985_wp   ! synchronous-machine + mech losses
        real(wp) :: bearing_loss_W         = 5.0e5_wp   ! fixed bearing + windage drag at speed
        ! Roll speed at which the turbine first reaches sync. Sets the
        ! "no-load steady-state" curve via target_rpm = rated × (tcv/roll)^(1/3),
        ! clamped by free_speed_overshoot. Operators roll at ~1 % TCV typical.
        real(wp) :: roll_tcv_norm          = 0.05_wp
    end type main_steam_config_t

    type :: main_steam_drivers_t
        real(wp) :: dome_pressure_Pa = 7.14e6_wp
        real(wp) :: power_current_W  = 0.0_wp     ! drives turbine eta
    end type main_steam_drivers_t

    type :: main_steam_command_t
        ! Sentinel −1 means "no change" for the analogue setpoints, so a
        ! default-initialised command is a no-op.
        real(wp) :: turbine_valve_set_pct = -1.0_wp
        real(wp) :: bypass_valve_set_pct  = -1.0_wp
        logical  :: tsv_fast_close        = .false.   ! latches turbine trip
        logical  :: tsv_trip_reset        = .false.
        logical  :: msiv_close            = .false.   ! latches isolation signal
        logical  :: msiv_open             = .false.   ! clear isolation
        ! Generator breaker controls (manual sync). sync_close demands the
        ! generator breaker close — only honoured if the rotor is inside
        ! the ±sync_window_rpm band around rated speed. sync_open trips the
        ! breaker (load reject) and drops the turbine back to free-spinning.
        logical  :: sync_close            = .false.
        logical  :: sync_open             = .false.
    end type main_steam_command_t

    type :: main_steam_observation_t
        real(wp) :: time                = 0.0_wp
        integer  :: n_steps             = 0
        real(wp) :: turbine_valve_pct   = 0.0_wp
        real(wp) :: bypass_valve_pct    = 0.0_wp
        real(wp) :: turbine_speed_rpm   = 0.0_wp
        real(wp) :: turbine_power_W     = 0.0_wp     ! electrical out (0 unless synced)
        real(wp) :: turbine_mech_W      = 0.0_wp     ! shaft mechanical power
        real(wp) :: steam_flow_norm     = 0.0_wp     ! normalised to rated (kg/s / rated_kg/s)
        real(wp) :: tcv_flow_norm       = 0.0_wp     ! steam to the turbine
        real(wp) :: bpv_flow_norm       = 0.0_wp     ! steam to the bypass (condenser)
        real(wp) :: srv_flow_kg_s       = 0.0_wp
        real(wp) :: srv_flow_norm       = 0.0_wp     ! normalised to rated
        logical  :: turbine_tripped     = .false.
        logical  :: generator_synced    = .false.
        logical  :: sync_ready          = .false.    ! inside sync window, can close breaker
        integer  :: n_srvs_open         = 0
        real(wp) :: msiv_open_frac      = 1.0_wp     ! sum(msiv_pos_pct)/400
    end type main_steam_observation_t

    type :: main_steam_state_t
        real(wp) :: time     = 0.0_wp
        integer  :: n_steps  = 0

        ! Operator demands (held; positions slew toward these).
        real(wp) :: turbine_valve_demand_pct = 0.0_wp
        real(wp) :: bypass_valve_demand_pct  = 0.0_wp

        ! Actual physical positions.
        real(wp) :: turbine_valve_pct        = 0.0_wp
        real(wp) :: bypass_valve_pct         = 0.0_wp
        real(wp) :: msiv_pos_pct(MS_N_MSLS)  = 100.0_wp
        logical  :: msiv_close_latched       = .false.

        ! Turbine dynamics.
        real(wp) :: turbine_speed_rpm   = 0.0_wp
        real(wp) :: turbine_power_W     = 0.0_wp   ! electrical
        real(wp) :: turbine_mech_W      = 0.0_wp   ! mechanical shaft power
        real(wp) :: tcv_flow_norm       = 0.0_wp   ! steam through control valve
        real(wp) :: bpv_flow_norm       = 0.0_wp   ! steam through bypass
        real(wp) :: steam_flow_norm     = 0.0_wp   ! total = tcv + bpv
        logical  :: turbine_tripped     = .false.
        logical  :: generator_synced    = .false.
        logical  :: tsv_fast_close_latched = .false.
        real(wp) :: tsv_fast_close_t       = 0.0_wp
        real(wp) :: tsv_fast_close_initial = 0.0_wp

        ! SRVs.
        logical  :: srv_open(MS_N_SRVS) = .false.
        real(wp) :: srv_flow_kg_s       = 0.0_wp
    end type main_steam_state_t

contains

    pure function clamp01(x, lo, hi) result(y)
        real(wp), intent(in) :: x, lo, hi
        real(wp) :: y
        y = max(lo, min(hi, x))
    end function clamp01

    subroutine main_steam_init(state, config)
        type(main_steam_state_t),  intent(out) :: state
        type(main_steam_config_t), intent(in)  :: config

        state%time     = 0.0_wp
        state%n_steps  = 0

        state%turbine_valve_demand_pct = clamp01(config%initial_turbine_valve_pct, 0.0_wp, 100.0_wp)
        state%bypass_valve_demand_pct  = clamp01(config%initial_bypass_valve_pct,  0.0_wp, 100.0_wp)
        state%turbine_valve_pct        = state%turbine_valve_demand_pct
        state%bypass_valve_pct         = state%bypass_valve_demand_pct

        state%msiv_pos_pct        = 100.0_wp
        state%msiv_close_latched  = .false.

        state%turbine_speed_rpm    = 0.0_wp
        state%turbine_power_W      = 0.0_wp
        state%turbine_mech_W       = 0.0_wp
        state%tcv_flow_norm        = 0.0_wp
        state%bpv_flow_norm        = 0.0_wp
        state%steam_flow_norm      = 0.0_wp
        state%turbine_tripped      = .false.
        state%generator_synced     = .false.
        state%tsv_fast_close_latched = .false.
        state%tsv_fast_close_t       = 0.0_wp
        state%tsv_fast_close_initial = 0.0_wp

        state%srv_open       = .false.
        state%srv_flow_kg_s  = 0.0_wp
    end subroutine main_steam_init

    subroutine main_steam_step(state, dt, drivers, config)
        type(main_steam_state_t),    intent(inout) :: state
        real(wp),                    intent(in)    :: dt
        type(main_steam_drivers_t),  intent(in)    :: drivers
        type(main_steam_config_t),   intent(in)    :: config

        integer  :: i, n_open
        real(wp) :: stroke_pct, max_step, target, ePres, msiv_frac
        real(wp) :: p_dome, p_reset_deadband_Pa, p_ratio
        real(wp) :: target_rpm, free_speed_max, accel_lag
        real(wp) :: net_torque_W, alpha

        p_dome = drivers%dome_pressure_Pa
        ! Pressure-ratio for the choked-flow valve model: steam only flows
        ! when there is a real ΔP between dome and condenser (~atmospheric
        ! sink). Using absolute dome P would give a phantom ~1.4 % rated
        ! flow at cold standby (1 atm dome). The differential form drives
        ! flow to zero at 1 atm and to 1.0 at rated dome P.
        p_ratio = max(0.0_wp, (p_dome - P_ATM_PA) &
                            / max(1.0_wp, config%rated_pressure_Pa - P_ATM_PA))

        ! ── TSV: fast-close override beats the normal demand slew. ────────
        if (state%tsv_fast_close_latched) then
            state%tsv_fast_close_t = state%tsv_fast_close_t + dt
            if (state%tsv_fast_close_t >= config%tsv_fast_close_time_s) then
                state%turbine_valve_pct = 0.0_wp
                state%turbine_tripped   = .true.
            else
                state%turbine_valve_pct = state%tsv_fast_close_initial &
                    * (1.0_wp - state%tsv_fast_close_t / config%tsv_fast_close_time_s)
            end if
        else
            max_step = config%tsv_normal_slew_pct_s * dt
            target   = state%turbine_valve_demand_pct
            if (target > state%turbine_valve_pct) then
                state%turbine_valve_pct = min(target, state%turbine_valve_pct + max_step)
            else
                state%turbine_valve_pct = max(target, state%turbine_valve_pct - max_step)
            end if
        end if

        ! ── Bypass valve: simple slew to demand. ──────────────────────────
        max_step = config%bypass_slew_pct_s * dt
        target   = state%bypass_valve_demand_pct
        if (target > state%bypass_valve_pct) then
            state%bypass_valve_pct = min(target, state%bypass_valve_pct + max_step)
        else
            state%bypass_valve_pct = max(target, state%bypass_valve_pct - max_step)
        end if

        ! ── MSIVs: stroke from current position toward latched setpoint. ──
        if (config%msiv_stroke_time_s > 0.0_wp) then
            stroke_pct = 100.0_wp / config%msiv_stroke_time_s * dt
        else
            stroke_pct = 100.0_wp
        end if
        do i = 1, MS_N_MSLS
            if (state%msiv_close_latched) then
                state%msiv_pos_pct(i) = max(0.0_wp, state%msiv_pos_pct(i) - stroke_pct)
            else
                state%msiv_pos_pct(i) = min(100.0_wp, state%msiv_pos_pct(i) + stroke_pct)
            end if
        end do
        msiv_frac = sum(state%msiv_pos_pct) / (real(MS_N_MSLS, wp) * 100.0_wp)

        ! ── SRV lift logic with deadband. ─────────────────────────────────
        p_reset_deadband_Pa = config%srv_reset_deadband_psi * PA_PER_PSI
        do i = 1, MS_N_SRVS
            if (.not. state%srv_open(i)) then
                if (p_dome >= config%srv_setpoint_Pa(i)) state%srv_open(i) = .true.
            else
                if (p_dome <= config%srv_setpoint_Pa(i) - p_reset_deadband_Pa) &
                    state%srv_open(i) = .false.
            end if
        end do
        n_open = count(state%srv_open)
        state%srv_flow_kg_s = real(n_open, wp) * config%srv_capacity_kg_s

        ! ── Steam-flow routing. ───────────────────────────────────────────
        ! Both flows are downstream of the MSIVs, so closing them stops
        ! delivery to both turbine and bypass. The bypass is hardware-
        ! limited to 25 % of rated flow (BYPASS_CAPACITY_FRAC applied in
        ! vessel) but here we expose its raw stroke as a fraction of rated.
        if (state%turbine_tripped) then
            state%tcv_flow_norm = 0.0_wp
        else
            state%tcv_flow_norm = (state%turbine_valve_pct / 100.0_wp) &
                                * msiv_frac * p_ratio
        end if
        state%bpv_flow_norm = (state%bypass_valve_pct / 100.0_wp) &
                            * 0.25_wp * msiv_frac * p_ratio
        state%steam_flow_norm = state%tcv_flow_norm + state%bpv_flow_norm

        ! ── Turbine mechanical power input. ──────────────────────────────
        ! Only TCV flow imparts torque on the shaft; bypass dumps direct
        ! to the condenser and contributes no shaft work. Mechanical power
        ! tracks the thermal power available to the turbine, scaled by
        ! eta_turbine and the pressure-ratio factor that proxies enthalpy
        ! drop across the stages.
        ePres = min(1.0_wp, 6.1e-8_wp * p_dome + 0.567_wp)
        state%turbine_mech_W = config%eta_turbine * ePres &
                             * state%tcv_flow_norm * drivers%power_current_W

        ! ── Speed dynamics. ──────────────────────────────────────────────
        ! Three regimes:
        !   (a) tripped → coast down at fixed deceleration
        !   (b) synced  → RPM pinned to rated; mechanical power becomes
        !                 electrical output (minus generator + bearing loss)
        !   (c) free-spinning → RPM relaxes toward target = f(TCV) with
        !                 first-order lag. No electrical output.
        if (state%turbine_tripped) then
            state%turbine_speed_rpm = max(0.0_wp, &
                state%turbine_speed_rpm - config%speed_decel_rpm_s * dt)
            state%turbine_power_W   = 0.0_wp
            state%generator_synced  = .false.

        else if (state%generator_synced) then
            ! Sync locks shaft to grid frequency. Real load follows steam
            ! supply; any TCV drop translates to dropped electrical out.
            state%turbine_speed_rpm = config%rated_speed_rpm
            net_torque_W = state%turbine_mech_W - config%bearing_loss_W
            state%turbine_power_W = max(0.0_wp, config%eta_generator * net_torque_W)
            ! Loss-of-steam reverse-power protection: if mech_W < bearing
            ! loss for any duration, generator would motor — auto-trip.
            if (state%turbine_mech_W < 0.5_wp * config%bearing_loss_W) then
                state%generator_synced = .false.
            end if

        else
            ! Free-spinning: no-load steady-state RPM grows with TCV flow.
            ! Cube root models the cubic friction law on a free turbine;
            ! at TCV = roll_tcv_norm the no-load speed reaches rated.
            ! Clamped above by free_speed_overshoot to limit the runaway.
            free_speed_max = config%rated_speed_rpm * (1.0_wp + config%free_speed_overshoot)
            if (state%tcv_flow_norm > 1.0e-4_wp) then
                target_rpm = config%rated_speed_rpm &
                           * (state%tcv_flow_norm / max(1.0e-4_wp, config%roll_tcv_norm))**(1.0_wp/3.0_wp)
                target_rpm = min(free_speed_max, target_rpm)
            else
                target_rpm = 0.0_wp
            end if
            accel_lag = max(1.0e-3_wp, config%speed_lag_time_s)
            alpha = 1.0_wp - exp(-dt / accel_lag)
            state%turbine_speed_rpm = state%turbine_speed_rpm &
                                   + alpha * (target_rpm - state%turbine_speed_rpm)
            state%turbine_speed_rpm = max(0.0_wp, state%turbine_speed_rpm)
            state%turbine_power_W = 0.0_wp
        end if

        state%time    = state%time + dt
        state%n_steps = state%n_steps + 1
    end subroutine main_steam_step

    subroutine main_steam_apply(state, command, config)
        type(main_steam_state_t),   intent(inout) :: state
        type(main_steam_command_t), intent(in)    :: command
        type(main_steam_config_t),  intent(in)    :: config

        real(wp) :: rpm_err

        if (command%turbine_valve_set_pct >= 0.0_wp) &
            state%turbine_valve_demand_pct = clamp01(command%turbine_valve_set_pct, 0.0_wp, 100.0_wp)

        if (command%bypass_valve_set_pct >= 0.0_wp) &
            state%bypass_valve_demand_pct = clamp01(command%bypass_valve_set_pct, 0.0_wp, 100.0_wp)

        if (command%tsv_fast_close .and. .not. state%tsv_fast_close_latched) then
            state%tsv_fast_close_latched = .true.
            state%tsv_fast_close_t       = 0.0_wp
            state%tsv_fast_close_initial = state%turbine_valve_pct
            state%generator_synced       = .false.   ! breaker opens on trip
        else if (command%tsv_trip_reset) then
            state%tsv_fast_close_latched = .false.
            state%tsv_fast_close_t       = 0.0_wp
            state%turbine_tripped        = .false.
        end if

        if (command%msiv_close) then
            state%msiv_close_latched = .true.
        else if (command%msiv_open) then
            state%msiv_close_latched = .false.
        end if

        ! Manual generator sync. Only honoured inside the ±sync_window
        ! around rated RPM and only if not tripped. Outside the window
        ! the command silently fails (frontend should observe sync_ready
        ! before sending close).
        if (command%sync_close .and. .not. state%turbine_tripped) then
            rpm_err = abs(state%turbine_speed_rpm - config%rated_speed_rpm)
            if (rpm_err <= config%sync_window_rpm) then
                state%generator_synced  = .true.
                state%turbine_speed_rpm = config%rated_speed_rpm  ! snap to grid
            end if
        end if
        if (command%sync_open) state%generator_synced = .false.
    end subroutine main_steam_apply

    function main_steam_observe(state, config) result(obs)
        type(main_steam_state_t),  intent(in) :: state
        type(main_steam_config_t), intent(in) :: config
        type(main_steam_observation_t) :: obs

        real(wp) :: rpm_err

        obs%time              = state%time
        obs%n_steps           = state%n_steps
        obs%turbine_valve_pct = state%turbine_valve_pct
        obs%bypass_valve_pct  = state%bypass_valve_pct
        obs%turbine_speed_rpm = state%turbine_speed_rpm
        obs%turbine_power_W   = state%turbine_power_W
        obs%turbine_mech_W    = state%turbine_mech_W
        obs%steam_flow_norm   = state%steam_flow_norm
        obs%tcv_flow_norm     = state%tcv_flow_norm
        obs%bpv_flow_norm     = state%bpv_flow_norm
        obs%srv_flow_kg_s     = state%srv_flow_kg_s
        obs%srv_flow_norm     = state%srv_flow_kg_s / max(1.0_wp, config%rated_steam_flow_kg_s)
        obs%turbine_tripped   = state%turbine_tripped
        obs%generator_synced  = state%generator_synced
        rpm_err = abs(state%turbine_speed_rpm - config%rated_speed_rpm)
        obs%sync_ready        = (.not. state%turbine_tripped) &
                              .and. (.not. state%generator_synced) &
                              .and. (rpm_err <= config%sync_window_rpm)
        obs%n_srvs_open       = count(state%srv_open)
        obs%msiv_open_frac    = sum(state%msiv_pos_pct) / (real(MS_N_MSLS, wp) * 100.0_wp)
    end function main_steam_observe

    subroutine main_steam_destroy(state)
        type(main_steam_state_t), intent(inout) :: state
        state%time     = 0.0_wp
        state%n_steps  = 0
        state%turbine_valve_demand_pct = 0.0_wp
        state%bypass_valve_demand_pct  = 0.0_wp
        state%turbine_valve_pct        = 0.0_wp
        state%bypass_valve_pct         = 0.0_wp
        state%msiv_pos_pct             = 0.0_wp
        state%msiv_close_latched       = .false.
        state%turbine_speed_rpm        = 0.0_wp
        state%turbine_power_W          = 0.0_wp
        state%turbine_mech_W           = 0.0_wp
        state%tcv_flow_norm            = 0.0_wp
        state%bpv_flow_norm            = 0.0_wp
        state%steam_flow_norm          = 0.0_wp
        state%turbine_tripped          = .false.
        state%generator_synced         = .false.
        state%tsv_fast_close_latched   = .false.
        state%tsv_fast_close_t         = 0.0_wp
        state%tsv_fast_close_initial   = 0.0_wp
        state%srv_open                 = .false.
        state%srv_flow_kg_s            = 0.0_wp
    end subroutine main_steam_destroy

end module main_steam