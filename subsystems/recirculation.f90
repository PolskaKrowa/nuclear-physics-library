! Recirculation subsystem.
!
! State ownership transfers from bwr_state_t:
!   mass_flux_core → recirc%mass_flux_kg_m2_s
!
! Jet-pump model (spec 2.4 p.5):
!   At rated, drive flow ≈ ⅓ × core flow, driven flow ≈ ⅔ × core flow.
!   M = ṁ_driven / ṁ_drive ≈ 2.0 ⇒ core flow = drive × (1 + M) = 3 × drive.
!   This module aggregates the 20 jet pumps into a single linear model
!   (each loop's 10 jet pumps share their loop's drive flow uniformly);
!   detailed nozzle/throat geometry stays a phase-2 open question.
!
! Pump-speed dynamics:
!   * Two motor-driven variable-speed recirc pumps (spec range 30 - 102 %).
!   * Demand slews toward operator setpoint at `speed_slew_pct_s`.
!   * Trip → speed coasts to zero at `coastdown_pct_s`.
!   * NPSH protection: FW flow < 20 % rated OR reactor water below L3
!     clamps effective speed to `runback_target_pct` (spec p.3).
!   * EOC-RPT latch: external signal trips both pumps immediately, stays
!     latched until operator-reset.
!
! Mass-flux conversion: at 100 % core flow the kernel sees
! `rated_mass_flux_kg_m2_s`; below rated the mass flux scales linearly
! with `core_flow_kg_s / rated_core_flow_kg_s`. Detailed core flow area
! / fuel-bundle hydraulics live in `fuel` (step 8).
module recirculation
    use kinds, only: wp
    implicit none
    private

    public :: recirc_state_t, recirc_config_t
    public :: recirc_command_t, recirc_observation_t
    public :: recirc_drivers_t
    public :: recirc_init, recirc_destroy
    public :: recirc_step, recirc_observe
    public :: recirc_apply
    public :: RECIRC_N_PUMPS, RECIRC_N_JET_PUMPS

    integer, parameter :: RECIRC_N_PUMPS     = 2
    integer, parameter :: RECIRC_N_JET_PUMPS = 20

    type :: recirc_config_t
        ! Rated core mass-flow (spec 2.1, 77.0 × 10⁶ lb/hr = 9700 kg/s).
        real(wp) :: rated_core_flow_kg_s     = 9700.0_wp
        ! Rated core mass flux for the legacy two_phase kernel knob.
        real(wp) :: rated_mass_flux_kg_m2_s  = 1500.0_wp
        ! Jet-pump M-ratio at rated (spec 2.4 p.5): driven/drive ≈ 2.0.
        real(wp) :: M_ratio                  = 2.0_wp
        ! Natural-circulation core-flow fraction with both pumps tripped.
        ! BWR/4 sees ~25 % rated core flow on density-difference head
        ! alone (chapter 2.4: "can't reach rated power on natural
        ! circulation"); used as the floor for jet-pump-less operation.
        real(wp) :: nat_circ_fraction        = 0.25_wp
        ! Operator speed limits (spec p.3).
        real(wp) :: min_pump_speed_pct       = 30.0_wp
        real(wp) :: max_pump_speed_pct       = 102.0_wp
        ! Variable-speed-drive slew rate (not in deck; abstract).
        real(wp) :: speed_slew_pct_s         = 5.0_wp
        ! Trip coastdown rate (not in deck; abstract — real pump inertia
        ! gives a multi-second coastdown).
        real(wp) :: coastdown_pct_s          = 10.0_wp
        ! NPSH-protection clamp (spec p.3).
        real(wp) :: runback_target_pct       = 30.0_wp
        real(wp) :: fw_npsh_threshold_pct    = 20.0_wp
        real(wp) :: level_npsh_threshold_m   = 0.3_wp  ! L3 setpoint
        real(wp) :: initial_pump_speed_pct   = 100.0_wp
    end type recirc_config_t

    type :: recirc_drivers_t
        ! Inputs from peer subsystems — set by orchestrator each tick.
        real(wp) :: feedwater_flow_pct = 100.0_wp
        real(wp) :: reactor_level_m    = 1.0_wp
        ! External EOC-RPT signal (main_steam turbine trip wired here by
        ! the orchestrator). Edge-detected against the latch in state.
        logical  :: eoc_rpt_signal     = .false.
    end type recirc_drivers_t

    type :: recirc_command_t
        ! −1 sentinel = no change. Operator demand applies to both pumps
        ! uniformly when `pump_speed_set_both` is set; for asymmetric
        ! operation use `pump_speed_set_pct(2)`.
        real(wp) :: pump_speed_set_both    = -1.0_wp
        real(wp) :: pump_speed_set_pct(RECIRC_N_PUMPS) = -1.0_wp
        integer  :: trip_pump_idx          = 0    ! 1..2; 0 = no action
        logical  :: reset_all_pumps        = .false.
        logical  :: eoc_rpt_reset          = .false.
    end type recirc_command_t

    type :: recirc_observation_t
        real(wp) :: time                 = 0.0_wp
        integer  :: n_steps              = 0
        real(wp) :: pump_speed_pct(RECIRC_N_PUMPS) = 0.0_wp
        real(wp) :: drive_flow_kg_s      = 0.0_wp
        real(wp) :: core_flow_kg_s       = 0.0_wp
        real(wp) :: core_flow_pct        = 0.0_wp
        real(wp) :: mass_flux_kg_m2_s    = 0.0_wp
        integer  :: pumps_running        = 0
        logical  :: eoc_rpt_latched      = .false.
        logical  :: npsh_runback_active  = .false.
    end type recirc_observation_t

    type :: recirc_state_t
        real(wp) :: time     = 0.0_wp
        integer  :: n_steps  = 0

        ! Pump-speed dynamics.
        real(wp) :: pump_speed_pct(RECIRC_N_PUMPS)        = 0.0_wp
        real(wp) :: pump_speed_demand_pct(RECIRC_N_PUMPS) = 0.0_wp
        logical  :: pump_tripped(RECIRC_N_PUMPS)          = .false.

        ! EOC-RPT (end-of-cycle recirc pump trip) latch.
        logical  :: eoc_rpt_latched = .false.

        ! NPSH runback flag — set when speed is being clamped.
        logical  :: npsh_runback_active = .false.

        ! Derived flow scalars (recomputed every step).
        real(wp) :: drive_flow_kg_s   = 0.0_wp
        real(wp) :: core_flow_kg_s    = 0.0_wp
        real(wp) :: mass_flux_kg_m2_s = 0.0_wp
    end type recirc_state_t

contains

    pure function clamp01(x, lo, hi) result(y)
        real(wp), intent(in) :: x, lo, hi
        real(wp) :: y
        y = max(lo, min(hi, x))
    end function clamp01

    subroutine recirc_init(state, config)
        type(recirc_state_t),  intent(out) :: state
        type(recirc_config_t), intent(in)  :: config

        real(wp) :: speed0
        integer  :: i

        speed0 = clamp01(config%initial_pump_speed_pct, 0.0_wp, config%max_pump_speed_pct)
        state%time    = 0.0_wp
        state%n_steps = 0
        state%pump_speed_pct        = speed0
        state%pump_speed_demand_pct = speed0
        state%pump_tripped          = .false.
        state%eoc_rpt_latched       = .false.
        state%npsh_runback_active   = .false.

        ! Seed the flow scalars so the orchestrator's first read sees the
        ! steady-state values rather than zero.
        state%drive_flow_kg_s = 0.0_wp
        do i = 1, RECIRC_N_PUMPS
            state%drive_flow_kg_s = state%drive_flow_kg_s &
                + (state%pump_speed_pct(i) / 100.0_wp) &
                * (config%rated_core_flow_kg_s / (1.0_wp + config%M_ratio)) &
                / real(RECIRC_N_PUMPS, wp)
        end do
        state%core_flow_kg_s = max(state%drive_flow_kg_s * (1.0_wp + config%M_ratio), &
                                   config%nat_circ_fraction * config%rated_core_flow_kg_s)
        state%mass_flux_kg_m2_s = state%core_flow_kg_s / config%rated_core_flow_kg_s &
                                * config%rated_mass_flux_kg_m2_s
    end subroutine recirc_init

    subroutine recirc_step(state, dt, drivers, config)
        type(recirc_state_t),   intent(inout) :: state
        real(wp),               intent(in)    :: dt
        type(recirc_drivers_t), intent(in)    :: drivers
        type(recirc_config_t),  intent(in)    :: config

        integer  :: i
        real(wp) :: target_speed, max_step, effective_demand
        real(wp) :: drive_per_pump_rated, drive_sum
        logical  :: npsh_protection

        ! ── EOC-RPT latch (edge-triggered on the external signal). ──────
        if (drivers%eoc_rpt_signal .and. .not. state%eoc_rpt_latched) then
            state%eoc_rpt_latched = .true.
            state%pump_tripped    = .true.
        end if

        ! ── NPSH protection (spec p.3): FW < 20 % OR L < L3 → runback. ──
        npsh_protection = (drivers%feedwater_flow_pct < config%fw_npsh_threshold_pct) &
                     .or. (drivers%reactor_level_m    < config%level_npsh_threshold_m)
        state%npsh_runback_active = npsh_protection

        ! ── Per-pump speed dynamics. ────────────────────────────────────
        do i = 1, RECIRC_N_PUMPS
            if (state%pump_tripped(i)) then
                ! Coastdown: speed drifts to zero at coastdown rate.
                max_step = config%coastdown_pct_s * dt
                state%pump_speed_pct(i) = max(0.0_wp, state%pump_speed_pct(i) - max_step)
                cycle
            end if

            ! Effective demand is the operator setpoint clamped by NPSH.
            effective_demand = state%pump_speed_demand_pct(i)
            if (npsh_protection) &
                effective_demand = min(effective_demand, config%runback_target_pct)
            ! Clamp to operating envelope.
            effective_demand = clamp01(effective_demand, 0.0_wp, config%max_pump_speed_pct)
            ! Below min: a running pump can't operate stably below the
            ! turndown limit; collapse to zero.
            if (effective_demand > 0.0_wp .and. effective_demand < config%min_pump_speed_pct) &
                effective_demand = config%min_pump_speed_pct

            target_speed = effective_demand
            max_step     = config%speed_slew_pct_s * dt
            if (target_speed > state%pump_speed_pct(i)) then
                state%pump_speed_pct(i) = min(target_speed, &
                                              state%pump_speed_pct(i) + max_step)
            else
                state%pump_speed_pct(i) = max(target_speed, &
                                              state%pump_speed_pct(i) - max_step)
            end if
        end do

        ! ── Aggregated jet-pump flow model. ─────────────────────────────
        ! Each loop's pump drives 10 jet pumps with a linear flow-vs-speed
        ! contract. Total drive flow at rated = rated_core_flow / (1 + M).
        drive_per_pump_rated = config%rated_core_flow_kg_s / (1.0_wp + config%M_ratio) &
                             / real(RECIRC_N_PUMPS, wp)
        drive_sum = 0.0_wp
        do i = 1, RECIRC_N_PUMPS
            drive_sum = drive_sum + (state%pump_speed_pct(i) / 100.0_wp) * drive_per_pump_rated
        end do
        state%drive_flow_kg_s = drive_sum

        ! Core flow = drive × (1 + M); floored by natural-circ contribution.
        state%core_flow_kg_s = max(state%drive_flow_kg_s * (1.0_wp + config%M_ratio), &
                                   config%nat_circ_fraction * config%rated_core_flow_kg_s)

        state%mass_flux_kg_m2_s = state%core_flow_kg_s / config%rated_core_flow_kg_s &
                                * config%rated_mass_flux_kg_m2_s

        state%time    = state%time + dt
        state%n_steps = state%n_steps + 1
    end subroutine recirc_step

    subroutine recirc_apply(state, command, config)
        type(recirc_state_t),   intent(inout) :: state
        type(recirc_command_t), intent(in)    :: command
        type(recirc_config_t),  intent(in)    :: config

        integer :: i
        real(wp) :: s

        if (command%pump_speed_set_both >= 0.0_wp) then
            s = clamp01(command%pump_speed_set_both, 0.0_wp, config%max_pump_speed_pct)
            state%pump_speed_demand_pct = s
        end if

        do i = 1, RECIRC_N_PUMPS
            if (command%pump_speed_set_pct(i) >= 0.0_wp) then
                state%pump_speed_demand_pct(i) = clamp01(command%pump_speed_set_pct(i), &
                                                          0.0_wp, config%max_pump_speed_pct)
            end if
        end do

        if (command%trip_pump_idx >= 1 .and. command%trip_pump_idx <= RECIRC_N_PUMPS) &
            state%pump_tripped(command%trip_pump_idx) = .true.

        if (command%reset_all_pumps) then
            state%pump_tripped = .false.
        end if

        if (command%eoc_rpt_reset) then
            state%eoc_rpt_latched = .false.
        end if
    end subroutine recirc_apply

    function recirc_observe(state, config) result(obs)
        type(recirc_state_t),  intent(in) :: state
        type(recirc_config_t), intent(in) :: config
        type(recirc_observation_t) :: obs
        obs%time                = state%time
        obs%n_steps             = state%n_steps
        obs%pump_speed_pct      = state%pump_speed_pct
        obs%drive_flow_kg_s     = state%drive_flow_kg_s
        obs%core_flow_kg_s      = state%core_flow_kg_s
        obs%core_flow_pct       = state%core_flow_kg_s / config%rated_core_flow_kg_s * 100.0_wp
        obs%mass_flux_kg_m2_s   = state%mass_flux_kg_m2_s
        obs%pumps_running       = count(.not. state%pump_tripped)
        obs%eoc_rpt_latched     = state%eoc_rpt_latched
        obs%npsh_runback_active = state%npsh_runback_active
    end function recirc_observe

    subroutine recirc_destroy(state)
        type(recirc_state_t), intent(inout) :: state
        state%time     = 0.0_wp
        state%n_steps  = 0
        state%pump_speed_pct        = 0.0_wp
        state%pump_speed_demand_pct = 0.0_wp
        state%pump_tripped          = .false.
        state%eoc_rpt_latched       = .false.
        state%npsh_runback_active   = .false.
        state%drive_flow_kg_s   = 0.0_wp
        state%core_flow_kg_s    = 0.0_wp
        state%mass_flux_kg_m2_s = 0.0_wp
    end subroutine recirc_destroy

end module recirculation