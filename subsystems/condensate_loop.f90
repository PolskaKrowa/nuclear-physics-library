! Condensate loop subsystem.
!
! The feedwater subsystem (subsystems/feedwater.f90) models the RFP
! turbine speed + low-suction trip cascade + heater train outlet T.
! This module owns the water inventory that sits *between* the hotwell
! and the RFP suction: hotwell level, CST level, makeup-flow controls,
! and the dissolved-oxygen state of the deaerating condenser. It feeds
! suction conditions back into feedwater (RFPT trip logic) and exposes
! the operator-visible water-loop indicators (hotwell level, CST level,
! make-up flow rate, dissolved-O₂ status).
!
! Hotwell (spec 2.6 p.5):
!   * 2 condenser shells; 24" equalising line ties hotwells.
!   * 2-min retention from normal to minimum level for N-16 decay.
!   * LCV-03 makeup from CST, LCV-02 supplemental makeup, LCV-04 reject
!     (high level) to CST via the demin outlet.
!   * Normal operating hotwell level: ~6 ft above suction strainer
!     (mid-band of 2 ft … 10 ft normal/operating range).
!
! CST (spec 2.6 p.13):
!   * 50 ft dia × 35 ft tall, 550,000 gal effective.
!   * 150,000 gal HPCI/RCIC reserve (guaranteed by tank-tap elevation).
!   * Vented to atmosphere; cold-water source for makeup + ECCS suction
!     (HPCI + RCIC).
!
! Deaeration in the BWR/4 is performed by the deaerating condenser
! itself (no separate deaerator vessel). What the operator monitors is
! dissolved-O₂ in the condensate (spec 8/9 ppb at normal vacuum + spray).
! We track that as a single scalar: it rises rapidly when the condenser
! breaks vacuum or air in-leaks, and falls slowly when condensate flow
! resumes (SJAEs pull non-condensibles out).
!
! Condensate Transfer Pumps (spec 2.6 p.13):
!   * Rated 1000 gpm @ 225 ft TDH each (1+1 standby).
!   * Auto-start at 330 gpm system makeup demand.
!   * Used to refill the hotwell from CST.
!
! Mass balance:
!   dHotwell/dt = ṁ_RFP_return_kg_s_negative + ṁ_makeup_from_CST - ṁ_RFP_suction
!   where ṁ_RFP_suction = feedwater.flow_kg_s in steady state.
!   In practice over a closed loop (LP-turbine exhaust → condenser →
!   condensate → demin → booster → RFP → vessel → MSL → turbine), the
!   hotwell level only drifts due to net inventory shift between the
!   reactor and the condenser. We model that drift as:
!     Δhotwell ∝ -(ṁ_FW - ṁ_steam)  (water moving up into the vessel
!       drops the hotwell; condensing steam returning from the
!       turbine refills it).
!
! For phase-1 we keep the model lumped — operator sees one hotwell
! level, one CST level. Future phases can split into per-shell or
! per-condenser-bank.
module condensate_loop
    use kinds, only: wp
    implicit none
    private

    public :: cond_state_t, cond_config_t
    public :: cond_command_t, cond_observation_t
    public :: cond_drivers_t
    public :: cond_init, cond_destroy
    public :: cond_step, cond_observe
    public :: cond_apply
    public :: COND_N_TRANSFER_PUMPS

    integer, parameter :: COND_N_TRANSFER_PUMPS = 2

    type :: cond_config_t
        ! Hotwell geometry (per-shell, but lumped into one effective vol).
        ! 2-min retention at rated FW flow ≈ 1638 kg/s → ~200 m³ inventory.
        real(wp) :: hotwell_volume_m3        = 200.0_wp
        real(wp) :: hotwell_cross_area_m2    = 50.0_wp     ! cumulative across shells
        real(wp) :: hotwell_normal_level_m   = 2.0_wp      ! mid-band normal
        real(wp) :: hotwell_low_level_m      = 0.5_wp      ! transfer-pump auto-start
        real(wp) :: hotwell_hi_level_m       = 3.0_wp      ! reject-to-CST trigger

        ! CST geometry (spec 2.6 p.13): 550,000 gal effective ≈ 2080 m³.
        real(wp) :: cst_volume_m3            = 2080.0_wp
        real(wp) :: cst_cross_area_m2        = 182.0_wp    ! 50 ft dia → π·(7.62)²
        real(wp) :: cst_initial_level_m      = 9.0_wp      ! near-full normal startup
        real(wp) :: cst_hpci_reserve_level_m = 2.45_wp     ! 150 kgal floor

        ! Transfer pumps (spec 2.6 p.13): 1000 gpm @ 225 ft TDH each.
        real(wp) :: transfer_pump_kg_s       = 63.0_wp     ! 1000 gpm
        real(wp) :: transfer_auto_start_kg_s = 20.8_wp     ! 330 gpm

        ! Deaeration (dissolved-O₂ in ppb):
        ! Normal: 8 ppb; degraded vacuum: rises; sustained good vacuum
        ! pulls toward normal. We expose this as a slow scalar driven by
        ! a (condenser_vacuum_ok ? -decay : +ingress) term.
        real(wp) :: o2_normal_ppb            = 8.0_wp
        real(wp) :: o2_ingress_rate_ppb_s    = 0.5_wp      ! when condenser vac lost
        real(wp) :: o2_decay_rate_ppb_s      = 0.05_wp     ! SJAE removal rate
        real(wp) :: o2_high_alarm_ppb        = 50.0_wp     ! operator alarm
    end type cond_config_t

    type :: cond_drivers_t
        ! Net mass flow from feedwater into the reactor (kg/s). Equals
        ! ṁ_FW (positive into vessel) − ṁ_steam_returning_from_turbine
        ! (positive into hotwell). At steady balance both are equal and
        ! the hotwell level is constant.
        real(wp) :: feedwater_flow_kg_s = 0.0_wp
        real(wp) :: steam_flow_kg_s     = 0.0_wp
        ! Condenser vacuum status. True = SJAEs maintaining vacuum →
        ! O₂ decays toward normal. False = vacuum broken / lost →
        ! O₂ rises. The MSS subsystem owns the turbine bypass dump
        ! that holds vacuum; condenser pressure is a future addition.
        logical  :: condenser_vacuum_ok = .true.
    end type cond_drivers_t

    type :: cond_command_t
        ! Manual override knobs. -1 = no change.
        integer :: transfer_pump_set(COND_N_TRANSFER_PUMPS) = -1
        integer :: trip_transfer_pump_idx                   = 0
        logical :: reset_transfer_pumps                     = .false.
        ! Makeup demand override (manual reject/fill). -1 = auto.
        real(wp) :: makeup_demand_kg_s_set                  = -1.0_wp
    end type cond_command_t

    type :: cond_observation_t
        real(wp) :: time                = 0.0_wp
        integer  :: n_steps             = 0
        real(wp) :: hotwell_level_m     = 0.0_wp
        real(wp) :: cst_level_m         = 0.0_wp
        real(wp) :: makeup_flow_kg_s    = 0.0_wp
        real(wp) :: reject_flow_kg_s    = 0.0_wp
        real(wp) :: o2_ppb              = 0.0_wp
        integer  :: transfer_running    = 0
        logical  :: hotwell_low_alarm   = .false.
        logical  :: cst_hpci_reserve_breached = .false.
        logical  :: o2_high_alarm       = .false.
    end type cond_observation_t

    type :: cond_state_t
        real(wp) :: time     = 0.0_wp
        integer  :: n_steps  = 0

        ! Levels.
        real(wp) :: hotwell_level_m = 0.0_wp
        real(wp) :: cst_level_m     = 0.0_wp

        ! Flows.
        real(wp) :: makeup_flow_kg_s = 0.0_wp
        real(wp) :: reject_flow_kg_s = 0.0_wp

        ! Dissolved-O₂ tracker.
        real(wp) :: o2_ppb = 0.0_wp

        ! Transfer pumps (2 redundant, 1 normally running).
        logical :: transfer_running(COND_N_TRANSFER_PUMPS) = .false.
        logical :: transfer_tripped(COND_N_TRANSFER_PUMPS) = .false.
    end type cond_state_t

contains

    subroutine cond_init(state, config)
        type(cond_state_t),  intent(out) :: state
        type(cond_config_t), intent(in)  :: config

        state%time              = 0.0_wp
        state%n_steps           = 0
        state%hotwell_level_m   = config%hotwell_normal_level_m
        state%cst_level_m       = config%cst_initial_level_m
        state%makeup_flow_kg_s  = 0.0_wp
        state%reject_flow_kg_s  = 0.0_wp
        state%o2_ppb            = config%o2_normal_ppb
        state%transfer_running  = .false.
        state%transfer_tripped  = .false.
        ! One transfer pump auto-runs continuously per spec p.13.
        state%transfer_running(1) = .true.
    end subroutine cond_init

    subroutine cond_step(state, dt, drivers, config)
        type(cond_state_t),   intent(inout) :: state
        real(wp),             intent(in)    :: dt
        type(cond_drivers_t), intent(in)    :: drivers
        type(cond_config_t),  intent(in)    :: config

        real(wp) :: net_to_vessel_kg_s, makeup_demand, reject_demand
        real(wp) :: dhotwell, dcst
        real(wp), parameter :: RHO_REF = 1000.0_wp

        ! Net mass leaving the hotwell toward the vessel = ṁ_FW. The
        ! steam returning from the turbine condenses back into the
        ! hotwell, contributing +ṁ_steam. Net hotwell drain rate:
        net_to_vessel_kg_s = drivers%feedwater_flow_kg_s - drivers%steam_flow_kg_s

        ! ── Makeup logic. ─────────────────────────────────────────────────
        ! If hotwell drops below the low-level threshold, demand makeup
        ! from CST. If above high-level, reject to CST. Lazy hysteresis.
        if (state%hotwell_level_m < config%hotwell_low_level_m) then
            makeup_demand = config%transfer_pump_kg_s
        else if (state%hotwell_level_m < config%hotwell_normal_level_m) then
            ! Partial makeup proportional to deviation below normal.
            makeup_demand = config%transfer_auto_start_kg_s &
                          * (config%hotwell_normal_level_m - state%hotwell_level_m) &
                          / max(0.01_wp, &
                                config%hotwell_normal_level_m - config%hotwell_low_level_m)
        else
            makeup_demand = 0.0_wp
        end if

        if (state%hotwell_level_m > config%hotwell_hi_level_m) then
            reject_demand = config%transfer_auto_start_kg_s
        else
            reject_demand = 0.0_wp
        end if

        ! Honour transfer-pump status (running + not tripped).
        if (count(state%transfer_running .and. .not. state%transfer_tripped) == 0) then
            makeup_demand = 0.0_wp
        end if
        ! CST reserve floor: can't draw below HPCI/RCIC reserve.
        if (state%cst_level_m <= config%cst_hpci_reserve_level_m) makeup_demand = 0.0_wp

        state%makeup_flow_kg_s = makeup_demand
        state%reject_flow_kg_s = reject_demand

        ! ── Hotwell mass balance. ─────────────────────────────────────────
        dhotwell = (state%makeup_flow_kg_s - state%reject_flow_kg_s - net_to_vessel_kg_s) &
                 / (max(1.0_wp, config%hotwell_cross_area_m2) * RHO_REF)
        state%hotwell_level_m = max(0.0_wp, state%hotwell_level_m + dhotwell * dt)

        ! ── CST mass balance. ─────────────────────────────────────────────
        dcst = (state%reject_flow_kg_s - state%makeup_flow_kg_s) &
             / (max(1.0_wp, config%cst_cross_area_m2) * RHO_REF)
        state%cst_level_m = max(0.0_wp, state%cst_level_m + dcst * dt)

        ! ── Dissolved-O₂ tracker. ─────────────────────────────────────────
        if (drivers%condenser_vacuum_ok) then
            state%o2_ppb = state%o2_ppb &
                         - config%o2_decay_rate_ppb_s * dt &
                         * (state%o2_ppb - config%o2_normal_ppb)
        else
            state%o2_ppb = state%o2_ppb + config%o2_ingress_rate_ppb_s * dt
        end if
        state%o2_ppb = max(0.0_wp, state%o2_ppb)

        state%time    = state%time + dt
        state%n_steps = state%n_steps + 1
    end subroutine cond_step

    subroutine cond_apply(state, command)
        type(cond_state_t),   intent(inout) :: state
        type(cond_command_t), intent(in)    :: command

        integer :: i

        do i = 1, COND_N_TRANSFER_PUMPS
            if (command%transfer_pump_set(i) == 0) then
                state%transfer_running(i) = .false.
            else if (command%transfer_pump_set(i) == 1 &
                     .and. .not. state%transfer_tripped(i)) then
                state%transfer_running(i) = .true.
            end if
        end do

        if (command%trip_transfer_pump_idx >= 1 .and. &
            command%trip_transfer_pump_idx <= COND_N_TRANSFER_PUMPS) &
            state%transfer_tripped(command%trip_transfer_pump_idx) = .true.

        if (command%reset_transfer_pumps) then
            state%transfer_tripped = .false.
            state%transfer_running = .false.
            state%transfer_running(1) = .true.
        end if

        if (command%makeup_demand_kg_s_set >= 0.0_wp) &
            state%makeup_flow_kg_s = command%makeup_demand_kg_s_set
    end subroutine cond_apply

    function cond_observe(state, config) result(obs)
        type(cond_state_t),  intent(in) :: state
        type(cond_config_t), intent(in) :: config
        type(cond_observation_t) :: obs

        obs%time             = state%time
        obs%n_steps          = state%n_steps
        obs%hotwell_level_m  = state%hotwell_level_m
        obs%cst_level_m      = state%cst_level_m
        obs%makeup_flow_kg_s = state%makeup_flow_kg_s
        obs%reject_flow_kg_s = state%reject_flow_kg_s
        obs%o2_ppb           = state%o2_ppb
        obs%transfer_running = count(state%transfer_running .and. .not. state%transfer_tripped)
        obs%hotwell_low_alarm = state%hotwell_level_m < config%hotwell_low_level_m
        obs%cst_hpci_reserve_breached = state%cst_level_m < config%cst_hpci_reserve_level_m
        obs%o2_high_alarm = state%o2_ppb > config%o2_high_alarm_ppb
    end function cond_observe

    subroutine cond_destroy(state)
        type(cond_state_t), intent(inout) :: state
        state%time              = 0.0_wp
        state%n_steps           = 0
        state%hotwell_level_m   = 0.0_wp
        state%cst_level_m       = 0.0_wp
        state%makeup_flow_kg_s  = 0.0_wp
        state%reject_flow_kg_s  = 0.0_wp
        state%o2_ppb            = 0.0_wp
        state%transfer_running  = .false.
        state%transfer_tripped  = .false.
    end subroutine cond_destroy

end module condensate_loop