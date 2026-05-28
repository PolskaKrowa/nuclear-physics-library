! Control rod drive subsystem.
!
! State ownership: this module now owns the rod position. The old
! scalar `bwr_state_t%rod_bank_position` (which fanned out to every cell
! in `update_cross_sections_feedback`) is gone; `bwr_c_interface` reads
! `sum(crd%blade_insertion) / CRD_N_BLADES` as the bank average for the
! single-position 1-D rod model still used by the cross-section feedback
! path. Per-cell control-rod XS using per-blade insertions lands in
! reorg step 8 (fuel subsystem).
!
! Convention: blade_insertion(i) ∈ [0, 1].
!   0.0 → blade fully withdrawn (tip below the active fuel)
!   1.0 → blade fully inserted  (tip at the top of the active fuel)
! Notch index follows the same monotonic relationship: notch 0 = fully
! withdrawn, notch 48 = fully inserted (6 in / 0.1524 m per notch,
! 48 × 6 in = 288 in = 24 ft of stroke per spec 2.3 p.26).
module crd
    use kinds, only: wp
    implicit none
    private

    public :: crd_state_t, crd_config_t
    public :: crd_command_t, crd_observation_t
    public :: crd_drivers_t
    public :: crd_init, crd_destroy
    public :: crd_step, crd_observe
    public :: crd_apply
    public :: scram_time_90pct
    public :: CRD_N_BLADES, CRD_N_NOTCHES, CRD_NOTCH_HEIGHT_M

    integer,  parameter :: CRD_N_BLADES        = 137
    integer,  parameter :: CRD_N_NOTCHES       = 48
    real(wp), parameter :: CRD_NOTCH_HEIGHT_M  = 0.1524_wp        ! 6 in per notch (spec 2.3 p.26)

    ! Pressure conversion constants used by the scram-time curve.
    real(wp), parameter :: PA_PER_PSI   = 6894.76_wp
    real(wp), parameter :: P_ATM_PA     = 101325.0_wp

    type :: crd_config_t
        ! Continuous-withdraw speed limit. Spec 2.3 doesn't tabulate this
        ! directly; the BWR/4 RMCS allows continuous motion at ~3 in/s
        ! (0.0762 m/s) during gang withdrawal, but operationally crews
        ! pull notch-by-notch with ~5 s dwell between notches → effective
        ! ~1.2 in/s (~0.03 m/s) during approach to critical. We use the
        ! slower number here: a 0.95 → 0.30 startup pull (Δ = 0.65 ×
        ! 24 ft stroke = 15.6 ft) at 1.2 in/s takes ~156 s ≈ 2.6 min,
        ! comfortably under the "startup < 1 hour" budget while leaving
        ! enough time for Doppler + xenon feedback to absorb each step
        ! and prevent the fuel-temp excursion the user flagged.
        real(wp) :: continuous_speed_m_s = 0.0305_wp              ! 1.2 in/s effective
        ! Active fuel height; passed through to keep `crd_config_t`
        ! self-contained even though geometry truly lives in `fuel`.
        real(wp) :: core_height_m        = 3.81_wp
        ! Initial bank insertion at reactor cold-iron. Matches the
        ! pre-reorg `bwr_initialize` value of 0.95 (mostly inserted —
        ! reactor startup configuration).
        real(wp) :: initial_insertion    = 0.95_wp
        ! Rod Block Monitor / Rod Worth Minimizer period setpoint.
        ! BWR/4 RBM blocks further withdrawal when the IRM/APRM period
        ! drops below ~30 s during the approach-to-critical (rule of thumb:
        ! 30 s gives the operator enough margin to detect overshoot
        ! before SCRAM-on-short-period at 5 s). Block applies to
        ! withdrawal only (insertion is always permitted). Period <= 0
        ! (the post-init default) is treated as "not yet measured" and
        ! does not block, otherwise the bank would never start moving
        ! on the first tick.
        real(wp) :: rod_block_period_s   = 30.0_wp
        ! Approach-to-critical block: freeze withdrawal once k_eff
        ! crosses this threshold. The period-based block above is one
        ! tick behind the solver (period is computed *after* the rod
        ! step decides whether to move), so during a fast power
        ! excursion period reads short only after the cliff. k_eff is
        ! computed by the same solver each tick — using it as the
        ! gating signal catches the approach before the excursion
        ! happens. Once at-power with stable feedback, the period-based
        ! block takes over.
        real(wp) :: rod_block_keff       = 0.99_wp
    end type crd_config_t

    type :: crd_drivers_t
        ! Reactor dome pressure feeds the scram-time curve via the
        ! "typical drive" envelope in Figure 2.3-13.
        real(wp) :: dome_pressure_Pa = 6.895e6_wp
        ! Current reactor period [s], from `update_instrumentation`. The
        ! Rod Block Monitor freezes withdrawal when this drops below
        ! `rod_block_period_s`. Negative period (decreasing flux) and
        ! period ≥ 1000 s (essentially stable / source-dominated
        ! subcritical) are treated as safe.
        real(wp) :: reactor_period_s = 9999.0_wp
        ! Current k_eff from the multigroup-diffusion solver. Used by
        ! the approach-to-critical block (faster trigger than period
        ! since k is available the same tick the rods would move).
        real(wp) :: k_eff = 0.0_wp
    end type crd_drivers_t

    type :: crd_command_t
        ! `scram_latch` latches the SCRAM (idempotent — subsequent latches
        ! while already scrammed are no-ops). `scram_reset` clears it.
        ! `bank_position_set >= 0` broadcasts that insertion fraction to
        ! all 137 blades (preserves the legacy C ABI).
        logical  :: scram_latch        = .false.
        logical  :: scram_reset        = .false.
        real(wp) :: bank_position_set  = -1.0_wp
    end type crd_command_t

    type :: crd_observation_t
        real(wp) :: time              = 0.0_wp
        integer  :: n_steps           = 0
        real(wp) :: bank_position_avg = 0.0_wp
        real(wp) :: bank_demand       = 0.0_wp
        logical  :: scrammed          = .false.
        logical  :: rod_block_active  = .false.
        real(wp) :: scram_t           = 0.0_wp
    end type crd_observation_t

    type :: crd_state_t
        real(wp) :: time                          = 0.0_wp
        integer  :: n_steps                       = 0

        ! Actual blade positions (what XS feedback reads).
        real(wp) :: blade_insertion(CRD_N_BLADES) = 0.0_wp
        integer  :: blade_notch(CRD_N_BLADES)     = 0
        ! Operator demand. `crd_apply` writes here; `crd_step` ramps
        ! blade_insertion toward this at continuous_speed_m_s. SCRAM
        ! drives the actual position directly and overrides the demand.
        real(wp) :: bank_demand                   = 0.95_wp

        ! Rod Block Monitor latch — true when withdrawal is currently
        ! held due to short reactor period. Exposed via observation so
        ! the frontend can light a "ROD BLOCK" lamp.
        logical  :: rod_block_active              = .false.

        ! Scram dynamics state.
        logical  :: scram_latched                 = .false.
        real(wp) :: scram_t                       = 0.0_wp
        real(wp) :: scram_initial(CRD_N_BLADES)   = 0.0_wp
    end type crd_state_t

contains

    pure function scram_time_90pct(p_dome_Pa) result(t90)
        ! Piecewise-linear fit to the "typical drive" envelope of
        ! Figure 2.3-13 (spec 03 chapter 2.3 p.23). Returns the time to
        ! 90 % insertion in seconds as a function of reactor dome
        ! pressure. Anchor points (from the figure):
        !
        !     pressure (psig)   t90 (s)
        !     ───────────────   ───────
        !     0                 3.5
        !     600               2.6
        !     1000              2.15
        !     1200              2.0
        !
        ! Above 1200 psig we hold the floor at 2.0 s — the figure plateaus.
        real(wp), intent(in) :: p_dome_Pa
        real(wp) :: t90, p_psig

        p_psig = max(0.0_wp, (p_dome_Pa - P_ATM_PA) / PA_PER_PSI)

        if (p_psig < 600.0_wp) then
            t90 = 3.5_wp + (2.6_wp - 3.5_wp) * (p_psig / 600.0_wp)
        else if (p_psig < 1000.0_wp) then
            t90 = 2.6_wp + (2.15_wp - 2.6_wp) * ((p_psig - 600.0_wp) / 400.0_wp)
        else if (p_psig < 1200.0_wp) then
            t90 = 2.15_wp + (2.0_wp - 2.15_wp) * ((p_psig - 1000.0_wp) / 200.0_wp)
        else
            t90 = 2.0_wp
        end if
    end function scram_time_90pct

    subroutine crd_init(state, config)
        type(crd_state_t),  intent(out) :: state
        type(crd_config_t), intent(in)  :: config

        state%time            = 0.0_wp
        state%n_steps         = 0
        state%blade_insertion = max(0.0_wp, min(1.0_wp, config%initial_insertion))
        state%blade_notch     = nint(state%blade_insertion(1) * real(CRD_N_NOTCHES, wp))
        state%bank_demand     = state%blade_insertion(1)
        state%scram_latched   = .false.
        state%scram_t         = 0.0_wp
        state%scram_initial   = state%blade_insertion
    end subroutine crd_init

    subroutine crd_step(state, dt, drivers, config)
        ! Advance the rod motion. Two regimes:
        !   (a) SCRAM latched → blades drive to fully-inserted following the
        !       Figure 2.3-13 envelope (pressure-dependent t90 → t_full).
        !   (b) Normal motion → ramp blade_insertion toward bank_demand at
        !       the continuous_speed_m_s rate. Operator commands change the
        !       demand instantly; the hardware stroke is what's rate-limited.
        type(crd_state_t),   intent(inout) :: state
        real(wp),            intent(in)    :: dt
        type(crd_drivers_t), intent(in)    :: drivers
        type(crd_config_t),  intent(in)    :: config

        real(wp) :: t90, t_full, frac, max_step, target, current_avg

        if (state%scram_latched) then
            state%scram_t = state%scram_t + dt

            ! Scram-time curve gives time to 90 % insertion. Extrapolate
            ! linearly to 100 %: full insertion at t90 / 0.9.
            t90    = scram_time_90pct(drivers%dome_pressure_Pa)
            t_full = max(1.0e-3_wp, t90 / 0.9_wp)

            frac = min(1.0_wp, state%scram_t / t_full)
            state%blade_insertion = state%scram_initial &
                                  + (1.0_wp - state%scram_initial) * frac
            state%blade_notch     = nint(state%blade_insertion * real(CRD_N_NOTCHES, wp))
        else
            ! Normal motion: rate-limited slew toward operator demand.
            ! Convert continuous_speed_m_s to fraction-per-second using the
            ! 24 ft total stroke (CRD_N_NOTCHES × CRD_NOTCH_HEIGHT_M ≈ 7.32 m).
            max_step = config%continuous_speed_m_s * dt &
                     / max(1.0e-3_wp, &
                           real(CRD_N_NOTCHES, wp) * CRD_NOTCH_HEIGHT_M)
            target = max(0.0_wp, min(1.0_wp, state%bank_demand))
            current_avg = state%blade_insertion(1)

            ! Approach-to-critical throttle: scale the withdraw rate by
            ! (1 - k_eff). Real BWR operators slow down as they approach
            ! critical for exactly this reason — a fixed rate means the
            ! reactivity insertion per second grows hyperbolically as the
            ! rod-worth curve steepens near the critical position. The
            ! gain of 5 means full speed below k=0.80, ~25 % speed at
            ! k=0.95, ~10 % at k=0.98, zero at k=1.0. Applies to withdrawal
            ! only — insertion stays at full speed for safe shutdown.
            if (target < current_avg .and. drivers%k_eff > 0.80_wp) then
                max_step = max_step &
                         * max(0.0_wp, min(1.0_wp, 5.0_wp * (1.0_wp - drivers%k_eff)))
            end if

            ! Rod Block Monitor: freeze withdrawal (target < current — rods
            ! are being pulled out) on either trigger:
            !   (1) k_eff already approaching critical from below — catches
            !       the criticality cliff before the prompt jump happens,
            !       since k is computed this same tick the rods would move
            !   (2) period below the configured safe minimum — backup
            !       trigger once at-power (post-criticality, k ~ 1, period
            !       becomes the operationally meaningful indicator)
            ! Insertion is always allowed; the block only halts continued
            ! withdrawal until either trigger releases. Negative and very
            ! large periods (no flux growth) bypass the period trigger —
            ! they mean the reactor is decaying or source-dominated.
            state%rod_block_active = (target < current_avg) .and. ( &
                (drivers%k_eff >= config%rod_block_keff) .or. &
                ((drivers%reactor_period_s > 0.0_wp) &
                    .and. (drivers%reactor_period_s < config%rod_block_period_s)))

            ! All 137 blades move together as a bank under bank-position
            ! commands. Per-blade addressing arrives with RMCS in phase 3.
            if (state%rod_block_active) then
                ! Hold position — operator must wait for period to recover.
                continue
            else if (target > current_avg) then
                current_avg = min(target, current_avg + max_step)
            else
                current_avg = max(target, current_avg - max_step)
            end if
            state%blade_insertion = current_avg
            state%blade_notch     = nint(current_avg * real(CRD_N_NOTCHES, wp))
        end if

        state%time    = state%time + dt
        state%n_steps = state%n_steps + 1
    end subroutine crd_step

    subroutine crd_apply(state, command)
        type(crd_state_t),   intent(inout) :: state
        type(crd_command_t), intent(in)    :: command

        if (command%scram_latch .and. .not. state%scram_latched) then
            state%scram_latched = .true.
            state%scram_t       = 0.0_wp
            state%scram_initial = state%blade_insertion
        else if (command%scram_reset .and. state%scram_latched) then
            state%scram_latched = .false.
            state%scram_t       = 0.0_wp
        end if

        ! Operator commands change the *demand* only — `crd_step` slews the
        ! actual blade position toward this at the hardware speed. Snapping
        ! blade_insertion directly (the pre-rate-limit behaviour) gave the
        ! operator an unrealistic step change in reactivity and was the
        ! root cause of the fuel-temp excursion when rods crossed the
        ! criticality cliff around 0.80 insertion.
        if (command%bank_position_set >= 0.0_wp) then
            state%bank_demand = max(0.0_wp, min(1.0_wp, command%bank_position_set))
        end if
    end subroutine crd_apply

    function crd_observe(state) result(obs)
        type(crd_state_t), intent(in) :: state
        type(crd_observation_t) :: obs
        obs%time              = state%time
        obs%n_steps           = state%n_steps
        obs%bank_position_avg = sum(state%blade_insertion) / real(CRD_N_BLADES, wp)
        obs%bank_demand       = state%bank_demand
        obs%scrammed          = state%scram_latched
        obs%rod_block_active  = state%rod_block_active
        obs%scram_t           = state%scram_t
    end function crd_observe

    subroutine crd_destroy(state)
        type(crd_state_t), intent(inout) :: state
        state%time            = 0.0_wp
        state%n_steps         = 0
        state%blade_insertion = 0.0_wp
        state%blade_notch     = 0
        state%scram_latched   = .false.
        state%scram_t         = 0.0_wp
        state%scram_initial   = 0.0_wp
    end subroutine crd_destroy

end module crd