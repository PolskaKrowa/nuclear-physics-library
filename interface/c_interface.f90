! interface/c_interface.f90
!
! C interface for the full BWR simulation.
! Wraps coupled multigroup-diffusion / heat-transfer / two-phase / burnup
! physics derived from bwr_core_simulator_opengl.f90 for use from C/C++.
!
! Usage from C:
!   #include "nuclear_physics_bwr.h"
!   BwrHandle bwr = bwr_reactor_create(20, 20, 20, 0.2375, 0.2375, 0.1855,
!                                       2381e6, 551.0, 7.14e6, 1500.0,
!                                       0.035, 3.71, 4.75);
!   bwr_solve_steady_state(bwr);
!   while (running) {
!       bwr_step(bwr, 0.02);   // 50 Hz
!       double keff = bwr_get_keff(bwr);
!       ...
!   }
!   bwr_reactor_destroy(bwr);
!
! Field arrays are laid out in Fortran column-major order (i varies fastest).
! Element at cell (i,j,k) (1-based) is out[i-1 + nx*(j-1 + ny*(k-1))].
!
module bwr_c_interface
    use iso_c_binding
    use kinds,             only: wp
    use multigroup_diffusion
    use cross_sections
    use burnup_depletion
    use heat_transfer
    use two_phase_flow
    use rng,               only: rng_seed
    use vessel,            only: vessel_state_t, vessel_config_t, vessel_drivers_t, &
                                 vessel_init, vessel_step_dynamics, &
                                 sat_temp_K, coolant_density_gl, calc_reactor_level
    use crd,               only: crd_state_t, crd_config_t, crd_command_t, &
                                 crd_drivers_t, crd_init, crd_step, crd_apply, &
                                 crd_observe, CRD_N_BLADES
    use main_steam,        only: main_steam_state_t, main_steam_config_t, &
                                 main_steam_command_t, main_steam_drivers_t, &
                                 main_steam_init, main_steam_step, main_steam_apply, &
                                 MS_N_MSLS
    use feedwater,         only: feedwater_state_t, feedwater_config_t, &
                                 feedwater_command_t, feedwater_drivers_t, &
                                 feedwater_init, feedwater_step, feedwater_apply, &
                                 FW_N_RFP
    use recirculation,     only: recirc_state_t, recirc_config_t, &
                                 recirc_command_t, recirc_drivers_t, &
                                 recirc_init, recirc_step, recirc_apply, &
                                 RECIRC_N_PUMPS
    use fuel,              only: fuel_state_t, fuel_config_t, &
                                 fuel_init_geometry, fuel_destroy, &
                                 fuel_setup_geometry, fuel_apply_to_xs, &
                                 fuel_power_to_volumetric_W_m3, &
                                 fuel_power_to_surface_q, &
                                 fuel_compute_convection_coefficients, &
                                 void_reactivity_factor
    use rhr,               only: rhr_state_t, rhr_config_t, &
                                 rhr_command_t, rhr_drivers_t, &
                                 rhr_init, rhr_step, rhr_apply, rhr_destroy, &
                                 RHR_N_LOOPS, &
                                 RHR_MODE_STANDBY, RHR_MODE_SHUTDOWN_COOLING, &
                                 RHR_MODE_SUPP_POOL_COOLING, RHR_MODE_CONTAINMENT_SPRAY, &
                                 RHR_MODE_LPCI
    use condensate_loop,   only: cond_state_t, cond_config_t, &
                                 cond_command_t, cond_drivers_t, &
                                 cond_init, cond_step, cond_apply, cond_destroy, &
                                 COND_N_TRANSFER_PUMPS
    implicit none
    private

    ! ── BWR simulation state ───────────────────────────────────────────────────
    type :: bwr_state_t
        integer  :: nx, ny, nz
        real(wp) :: dx, dy, dz
        real(wp) :: core_height, core_diameter

        type(mg_state_t)         :: neutronics
        type(heat_state_t)       :: heat
        type(two_phase_state_t)  :: thermalhydraulics
        type(burnup_state_t)     :: burnup
        type(xsec_library_t)     :: xslib

        real(wp) :: power_rated, power_current
        real(wp) :: inlet_temperature

        ! Vessel subsystem (spec 2.1). Owns dome pressure, level, sat-temp,
        ! ODE history, FW-inlet temp filter. See subsystems/vessel.f90.
        type(vessel_state_t)  :: vessel
        type(vessel_config_t) :: vessel_config

        real(wp) :: k_eff, reactivity_pcm
        real(wp) :: max_fuel_temp
        real(wp) :: avg_void_fraction, max_void_fraction
        real(wp) :: min_chfr, avg_burnup

        real(wp) :: time
        integer  :: n_steps

        real(wp) :: avg_coolant_temp

        real(wp) :: neutrons_prev, reactor_period
        real(wp) :: aprm

        ! Fuel subsystem (spec 2.2). Owns the lattice geometry, the
        ! per-cell XS-feedback path, the 137-blade fan-out, and the
        ! operator reactivity perturbation. See subsystems/fuel.f90.
        type(fuel_state_t)  :: fuel
        type(fuel_config_t) :: fuel_config

        ! Control rod drive subsystem (spec 2.3). Owns the 137 blade
        ! insertions and the scram-time dynamics. See subsystems/crd.f90.
        type(crd_state_t)  :: crd
        type(crd_config_t) :: crd_config

        ! Main steam subsystem (spec 2.5). Owns MSIV/SRV/bypass + turbine
        ! throttle. See subsystems/main_steam.f90.
        type(main_steam_state_t)  :: ms
        type(main_steam_config_t) :: ms_config

        ! Feedwater subsystem (spec 2.6). Owns the three pump chains,
        ! staggered suction-trip cascade, RFPT speed dynamics, and the
        ! heater-train outlet temperature. See subsystems/feedwater.f90.
        type(feedwater_state_t)   :: fw
        type(feedwater_config_t)  :: fw_config

        ! Recirculation subsystem (spec 2.4). Owns the two-loop pump
        ! topology, the aggregated 20 jet-pump M-ratio model, NPSH
        ! runback, and the EOC-RPT latch. See subsystems/recirculation.f90.
        type(recirc_state_t)      :: recirc
        type(recirc_config_t)     :: recirc_config

        ! RHR subsystem (spec 10.4). Four loops, multi-mode (SDC, pool
        ! cooling, containment spray, LPCI). See subsystems/rhr.f90.
        type(rhr_state_t)         :: rhr
        type(rhr_config_t)        :: rhr_config

        ! Condensate loop subsystem (spec 2.6 water-side). Owns hotwell
        ! level, CST level, makeup/reject, transfer pumps, dissolved-O₂.
        ! See subsystems/condensate_loop.f90.
        type(cond_state_t)        :: cond
        type(cond_config_t)       :: cond_config

        ! Approximate suppression-pool temperature [K]. Phase-1 placeholder
        ! until containment subsystem lands. Driven by SRV discharge +
        ! RHR pool-cooling heat removal.
        real(wp) :: supp_pool_T_K
        real(wp) :: drywell_pressure_Pa

        ! Sustained operator reactivity insertion [pcm]. Reapplied to nu*sigma_f
        ! every step by update_cross_sections_feedback. bwr_apply_reactivity adds
        ! to this; bwr_reset_reactivity zeros it.
        real(wp) :: reactivity_perturbation_pcm

        ! Delayed-neutron precursor concentrations, 6 groups (U-235 standard).
        ! Advanced by mg_solve_transient each tick. Seeded to equilibrium
        ! C_d = beta_d * F / lambda_d at the end of do_solve_steady_state.
        real(wp), allocatable :: precursors(:, :, :, :)
    end type bwr_state_t

    ! ── Safety trip thresholds ─────────────────────────────────────────────────
    ! Single source of truth for the conditions that would SCRAM a real BWR.
    ! Values are RPS scram setpoints from BWR/4 Tech Specs (Browns Ferry /
    ! Peach Bottom class), not ASME design margins. bwr_get_trip_flags()
    ! returns a bitmask of which of these are exceeded; the GDExtension uses
    ! that to latch the SCRAM and to drive pre-trip alarms.
    real(wp), parameter :: TRIP_FUEL_TEMP_MAX_K = 1673.15_wp  ! 1400 °C peak fuel
    real(wp), parameter :: TRIP_POWER_OVERSHOOT = 1.30_wp     ! × rated thermal
    real(wp), parameter :: TRIP_PERIOD_MIN_S    = 5.0_wp      ! prompt-supercritical risk
    ! L3 reactor scram setpoint: roughly +12 in above top of active fuel.
    real(wp), parameter :: TRIP_LEVEL_MIN_M     = 0.3_wp
    ! High-pressure reactor scram (~1080 psig); below ASME 110 % limit.
    real(wp), parameter :: TRIP_PRESSURE_MAX_PA = 7.45e6_wp
    ! Low-MSL-pressure MSIV isolation in RUN mode (~798 psig).
    real(wp), parameter :: TRIP_PRESSURE_MIN_PA = 5.50e6_wp
    real(wp), parameter :: TRIP_CHFR_MIN        = 1.0_wp

    integer(c_int), parameter :: BWR_TRIP_FUEL_TEMP_HIGH = 1    ! 1 << 0
    integer(c_int), parameter :: BWR_TRIP_POWER_HIGH     = 2    ! 1 << 1
    integer(c_int), parameter :: BWR_TRIP_SHORT_PERIOD   = 4    ! 1 << 2
    integer(c_int), parameter :: BWR_TRIP_LEVEL_LOW      = 8    ! 1 << 3
    integer(c_int), parameter :: BWR_TRIP_PRESSURE_HIGH  = 16   ! 1 << 4
    integer(c_int), parameter :: BWR_TRIP_PRESSURE_LOW   = 32   ! 1 << 5
    integer(c_int), parameter :: BWR_TRIP_CHFR_LOW       = 64   ! 1 << 6

contains

    ! ══════════════════════════════════════════════════════════════════════════
    !  INTERNAL PHYSICS SUBROUTINES
    ! ══════════════════════════════════════════════════════════════════════════
    !
    ! `sat_temp_K`, `coolant_density_gl`, `calc_reactor_level` now live in
    ! `subsystems/vessel.f90` and are re-exported here through `use vessel`.
    ! `void_reactivity_factor`, the per-cell XS-feedback path, the
    ! lattice/core-mask geometry, and the operator reactivity perturbation
    ! moved to `subsystems/fuel.f90` in reorg step 8.

    pure function avg_coolant_temp_K(T_field, is_fuel_cell_xy) result(T_avg)
        ! Bulk coolant temperature — averages only the non-fuel cells of
        ! the heat grid. Previously we averaged the full grid, which folded
        ! the fuel temperature into the "coolant" reading. That created a
        ! positive-feedback runaway at cold standby: any tiny heat in the
        ! fuel raised the apparent T_bulk, the convection helper then drove
        ! the fuel toward that higher T_bulk, ad infinitum. Restricting the
        ! average to coolant cells breaks the loop. Falls back to the field
        ! average if no coolant cells are flagged (defensive — should never
        ! happen for a valid lattice).
        real(wp), intent(in) :: T_field(:, :, :)
        logical,  intent(in) :: is_fuel_cell_xy(:, :)
        real(wp) :: T_avg

        integer  :: i, j, k, n_cool
        real(wp) :: t_sum
        t_sum  = 0.0_wp
        n_cool = 0
        do k = 1, size(T_field, 3)
            do j = 1, size(T_field, 2)
                do i = 1, size(T_field, 1)
                    if (.not. is_fuel_cell_xy(i, j)) then
                        t_sum  = t_sum + T_field(i, j, k)
                        n_cool = n_cool + 1
                    end if
                end do
            end do
        end do
        if (n_cool > 0) then
            T_avg = t_sum / real(n_cool, wp)
        else
            T_avg = sum(T_field) / real(size(T_field), wp)
        end if
    end function avg_coolant_temp_K

    subroutine update_cross_sections_feedback(sim, T, rho)
        ! Thin wrapper that hands off to the fuel subsystem. Kept as a
        ! named subroutine because both do_solve_steady_state and do_step
        ! call it twice with the same argument shape.
        type(bwr_state_t), intent(inout) :: sim
        real(wp), intent(in) :: T(:, :, :)
        real(wp), intent(in) :: rho(:, :, :)

        call fuel_apply_to_xs(sim%fuel, sim%fuel_config, sim%xslib, sim%neutronics, &
                              sim%burnup%burnup, sim%burnup%Xe135, sim%burnup%Sm149, &
                              sim%crd%blade_insertion, T, rho, &
                              sim%avg_void_fraction / 100.0_wp, &
                              sim%reactivity_perturbation_pcm)
    end subroutine update_cross_sections_feedback

    subroutine update_instrumentation(sim)
        type(bwr_state_t), intent(inout) :: sim
        real(wp) :: ratio

        sim%aprm = (sim%power_current / max(1.0_wp, sim%power_rated)) * 100.0_wp
        sim%aprm = max(0.0_wp, min(200.0_wp, sim%aprm))

        ! Period detector — only meaningful in the IRM/APRM power range.
        ! Real BWR period meters drive off the linear IRM signal, which has
        ! a calibrated floor (≥ ~mWt-equivalent count rate). Below that
        ! floor the reactor is in the SRM regime (counts, not period).
        ! Without the floor here, microscopic flux variations from
        ! numerical noise at cold standby produce wildly oscillating
        ! negative periods (-30 s, -200 s, ...) that look like a runaway
        ! to the operator but reflect nothing physical — flux changing
        ! by parts-per-million between ticks at sub-watt power.
        ! Floor at 1 kW: source-multiplied cold standby (≪ 1 W) shows
        ! the default "stable" period; once the reactor is actually
        ! climbing through the kW range, the indication kicks in.
        ! Additionally, require |log(ratio)| > 1e-5 so a near-stable
        ! reactor at any power doesn't produce a giant noisy period
        ! magnitude from the 1/log singularity.
        block
            real(wp) :: log_ratio
            real(wp), parameter :: PERIOD_POWER_FLOOR_W = 1.0e3_wp
            real(wp), parameter :: PERIOD_LOG_FLOOR     = 1.0e-5_wp
            if (sim%neutrons_prev > PERIOD_POWER_FLOOR_W &
                .and. sim%power_current > PERIOD_POWER_FLOOR_W) then
                ratio = sim%power_current / sim%neutrons_prev
                if (ratio > 0.0_wp) then
                    log_ratio = log(ratio)
                    if (abs(log_ratio) > PERIOD_LOG_FLOOR) then
                        sim%reactor_period = 1.0_wp / log_ratio
                        sim%reactor_period = max(-999.0_wp, min(999.0_wp, sim%reactor_period))
                    else
                        sim%reactor_period = 9999.0_wp
                    end if
                end if
            else
                ! Below the indicator's calibrated range — report the
                ! default rather than holding the last computed value
                ! from when power was higher (which would be misleading).
                sim%reactor_period = 9999.0_wp
            end if
        end block
        sim%neutrons_prev = sim%power_current
    end subroutine update_instrumentation

    subroutine update_pressure_dynamics(sim, dt)
        ! Thin orchestration wrapper. Packs the upstream drivers the vessel
        ! ODE needs, invokes vessel_step_dynamics on the subsystem state.
        ! Valve positions + SRV blowdown are taken from main_steam's
        ! previous-tick state (it runs after vessel each tick); at 50 Hz
        ! the one-tick lag is well below SRV/MSIV stroke timescales.
        type(bwr_state_t), intent(inout) :: sim
        real(wp), intent(in) :: dt

        type(vessel_drivers_t) :: drivers

        drivers%avg_coolant_temp_K  = sim%avg_coolant_temp
        drivers%power_current_W     = sim%power_current
        drivers%power_rated_W       = sim%power_rated
        drivers%feedwater_flow_pct  = sim%fw%flow_pct
        drivers%turbine_valve_pct   = sim%ms%turbine_valve_pct
        drivers%bypass_valve_pct    = sim%ms%bypass_valve_pct
        drivers%feedwater_flow_kg_s = sim%fw%flow_kg_s
        drivers%steam_flow_kg_s     = sim%ms%steam_flow_norm &
                                    * sim%ms_config%rated_steam_flow_kg_s &
                                    + sim%ms%srv_flow_kg_s
        drivers%srv_flow_frac       = sim%ms%srv_flow_kg_s &
                                    / max(1.0_wp, sim%ms_config%rated_steam_flow_kg_s)
        drivers%msiv_open_frac      = sum(sim%ms%msiv_pos_pct) &
                                    / (real(MS_N_MSLS, wp) * 100.0_wp)
        drivers%turbine_tripped     = sim%ms%turbine_tripped
        drivers%core_height_m       = sim%core_height

        call vessel_step_dynamics(sim%vessel, dt, drivers, sim%vessel_config)
    end subroutine update_pressure_dynamics

    subroutine update_feedwater(sim, dt)
        ! Wrap feedwater_step into the orchestrator's tick. Runs after
        ! main_steam each tick — reads the freshly-stepped vessel level
        ! and the main-steam steam-flow norm for the 3-element controller.
        type(bwr_state_t), intent(inout) :: sim
        real(wp), intent(in) :: dt

        type(feedwater_drivers_t) :: fw_drivers
        fw_drivers%reactor_level_m = sim%vessel%Lrx
        fw_drivers%steam_flow_norm = sim%ms%steam_flow_norm
        fw_drivers%dome_pressure_Pa = sim%vessel%pressure_operating
        call feedwater_step(sim%fw, dt, fw_drivers, sim%fw_config)
    end subroutine update_feedwater

    subroutine update_main_steam(sim, dt)
        ! Wrap main_steam_step into the orchestrator's tick. Runs after
        ! vessel each tick — uses the freshly-stepped dome pressure to
        ! drive the SRV lift logic and the turbine eta.
        type(bwr_state_t), intent(inout) :: sim
        real(wp), intent(in) :: dt

        type(main_steam_drivers_t) :: ms_drivers
        ms_drivers%dome_pressure_Pa = sim%vessel%pressure_operating
        ms_drivers%power_current_W  = sim%power_current
        call main_steam_step(sim%ms, dt, ms_drivers, sim%ms_config)
    end subroutine update_main_steam

    subroutine update_rhr(sim, dt)
        ! Tick the RHR subsystem. Pool temperature climbs from SRV
        ! discharge and falls under pool-cooling RHR; SDC pulls heat
        ! from the reactor coolant directly.
        type(bwr_state_t), intent(inout) :: sim
        real(wp), intent(in) :: dt
        type(rhr_drivers_t) :: rhr_drivers
        real(wp), parameter :: CP_WATER  = 4186.0_wp
        ! Suppression pool BWR/4: ~3000 m³ water at 86 °F design. We use
        ! a fixed thermal capacity for the heat-up / cool-down balance.
        real(wp), parameter :: M_POOL_KG = 3.0e6_wp
        real(wp) :: q_srv_W, q_rhr_pool_W, q_rhr_core_W

        rhr_drivers%dome_pressure_Pa     = sim%vessel%pressure_operating
        rhr_drivers%reactor_coolant_T_K  = sim%avg_coolant_temp
        rhr_drivers%suppression_pool_T_K = sim%supp_pool_T_K
        rhr_drivers%drywell_pressure_Pa  = sim%drywell_pressure_Pa
        rhr_drivers%reactor_level_m      = sim%vessel%Lrx
        call rhr_step(sim%rhr, dt, rhr_drivers, sim%rhr_config)

        ! Heat-into-pool from SRV discharge (estimate: ṁ × Δh_fg ≈
        ! ṁ × 1.8e6 J/kg from sat steam to subcooled pool water).
        q_srv_W = sim%ms%srv_flow_kg_s * 1.8e6_wp

        ! Heat removed by RHR depends on mode:
        !   SDC          → out of reactor coolant
        !   pool cooling → out of suppression pool
        !   containment spray, LPCI → out of pool (treat as pool sink)
        q_rhr_pool_W = 0.0_wp
        q_rhr_core_W = 0.0_wp
        select case (sim%rhr%mode)
        case (RHR_MODE_SHUTDOWN_COOLING)
            q_rhr_core_W = sim%rhr%total_heat_W
        case (RHR_MODE_SUPP_POOL_COOLING, RHR_MODE_CONTAINMENT_SPRAY, RHR_MODE_LPCI)
            q_rhr_pool_W = sim%rhr%total_heat_W
        end select

        ! Pool temperature ODE.
        sim%supp_pool_T_K = sim%supp_pool_T_K &
            + (q_srv_W - q_rhr_pool_W) / (M_POOL_KG * CP_WATER) * dt
        sim%supp_pool_T_K = max(280.0_wp, min(420.0_wp, sim%supp_pool_T_K))

        ! SDC removes power from the reactor; reflect this by subtracting
        ! from the power-balance proxy. At cold standby this term is
        ! near-zero (no T difference). The vessel pressure ODE already
        ! reacts to coolant T, so reducing avg_coolant_temp here closes
        ! the loop. We use a small fraction so the SDC effect on
        ! avg_coolant_temp is metered and doesn't fight the heat solver
        ! every tick.
        if (q_rhr_core_W > 0.0_wp) then
            ! Approximate vessel water inventory ≈ 200 m³ → 1.5e5 kg
            sim%avg_coolant_temp = sim%avg_coolant_temp &
                - q_rhr_core_W / (1.5e5_wp * CP_WATER) * dt
        end if
    end subroutine update_rhr

    subroutine update_condensate_loop(sim, dt)
        ! Tick the condensate loop: hotwell level, CST level, makeup,
        ! reject, dissolved-O₂.
        type(bwr_state_t), intent(inout) :: sim
        real(wp), intent(in) :: dt
        type(cond_drivers_t) :: cond_drivers
        cond_drivers%feedwater_flow_kg_s = sim%fw%flow_kg_s
        cond_drivers%steam_flow_kg_s     = sim%ms%steam_flow_norm &
                                         * sim%ms_config%rated_steam_flow_kg_s
        ! Phase-1: assume condenser vacuum is held as long as the turbine
        ! is not tripped and MSIVs are open. Future condenser subsystem
        ! will own this properly.
        cond_drivers%condenser_vacuum_ok = .not. sim%ms%turbine_tripped &
                                           .and. sim%ms%msiv_pos_pct(1) > 50.0_wp
        call cond_step(sim%cond, dt, cond_drivers, sim%cond_config)
    end subroutine update_condensate_loop

    ! ── Physics initialisation ────────────────────────────────────────────────

    subroutine init_physics(sim, enrichment)
        type(bwr_state_t), intent(inout) :: sim
        real(wp), intent(in) :: enrichment

        type(mg_config_t)        :: neutron_config
        type(heat_config_t)      :: heat_config
        type(two_phase_config_t) :: th_config
        type(burnup_config_t)    :: burnup_config
        type(xsec_material_t)    :: xsec_fuel, xsec_water

        ! ── Cross-section library ──────────────────────────────────────────────
        call xslib_init(sim%xslib, n_groups=2)

        xsec_fuel%name       = "UO2_35"
        xsec_fuel%n_groups   = 2
        xsec_fuel%is_fuel    = .true.
        xsec_fuel%enrichment = enrichment
        xsec_fuel%T_ref      = 900.0_wp
        xsec_fuel%rho_ref    = 10.97_wp
        call xslib_create_two_group_fuel(xsec_fuel%xsec_base, enrichment=enrichment)
        allocate(xsec_fuel%alpha_D(2),   xsec_fuel%alpha_mod(2))
        allocate(xsec_fuel%alpha_rho(2), xsec_fuel%alpha_void(2))
        xsec_fuel%alpha_D    = [-2.0e-5_wp, -3.0e-5_wp]
        xsec_fuel%alpha_mod  = [ 0.0_wp,    0.0_wp    ]
        xsec_fuel%alpha_rho  = [ 0.0_wp,    0.0_wp    ]
        xsec_fuel%alpha_void = [ 0.0_wp,    0.0_wp    ]
        call xslib_add_material(sim%xslib, xsec_fuel)

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

        ! ── Neutronics ─────────────────────────────────────────────────────────
        neutron_config%n_groups        = 2
        neutron_config%max_outer_iter  = 100
        neutron_config%outer_tolerance = 1.0e-5_wp
        neutron_config%power_level     = sim%power_rated
        neutron_config%normalize_power = .true.
        call mg_init(sim%neutronics, sim%nx, sim%ny, sim%nz, &
             sim%dx, sim%dy, sim%dz, neutron_config)
        call fuel_setup_geometry(sim%fuel, sim%fuel_config, sim%xslib, &
                                 sim%neutronics, sim%nz)

        ! Delayed-neutron precursors (6-group U-235). Start at zero; the
        ! steady-state seed sets these to equilibrium once a flux shape exists.
        allocate(sim%precursors(sim%nx, sim%ny, sim%nz, 6))
        sim%precursors = 0.0_wp

        ! Startup neutron source. BWR/4 plants carry a Cf-252 / Sb-Be source
        ! plus spontaneous-fission and (γ,n) contributions from the core that
        ! together produce ~10⁷–10⁸ n/s. The source ensures subcritical flux
        ! stays finite (source × 1/(1-k) ≈ 25× at k=0.96) so SRM/IRM/period
        ! indications respond during approach to critical rather than the
        ! solver decaying to numerical zero. Magnitude chosen so subcritical
        ! power at k=0.96 lands a few μW (≫ 1 mW period-detector floor when
        ! near critical, ≪ 1 W so the at-power steady state is unperturbed).
        block
            real(wp) :: startup_source(sim%nx, sim%ny, sim%nz)
            integer  :: i, j, k
            startup_source = 0.0_wp
            ! Apply to fuel cells only via the 2-D fuel mask. Magnitude
            ! tuned down from initial 10 n/cm³·s — at that level the
            ! source-multiplied flux at cold standby accumulated enough
            ! over ~70 s to seed a localised flux hotspot that ran away
            ! to ~330 K in a single cell, even with zero net rod motion.
            ! 0.1 n/cm³·s ≈ 1e6 n/s over the active core volume — still
            ! enough for SRM counts to register during the approach to
            ! critical, but not enough to bootstrap a flux excursion at
            ! deep subcritical (k≈0.03). Thermal-group only (group 2).
            do k = 1, sim%nz
                do j = 1, sim%ny
                    do i = 1, sim%nx
                        if (sim%fuel%is_fuel_cell(i, j)) startup_source(i, j, k) = 0.1_wp
                    end do
                end do
            end do
            call mg_set_source(sim%neutronics, startup_source, group=2)
        end block

        ! ── Heat transfer ──────────────────────────────────────────────────────
        heat_config%include_convection = .true.
        call heat_init(sim%heat, sim%nx, sim%ny, sim%nz, &
                       sim%dx, sim%dy, sim%dz, heat_config)
        call heat_set_properties(sim%heat, k=3.0_wp, rho=10970.0_wp, cp=300.0_wp, &
            i1=1, i2=int(0.9_wp*sim%nx), j1=1, j2=int(0.9_wp*sim%ny))
        call heat_set_properties(sim%heat, k=0.6_wp, rho=738.0_wp, cp=5200.0_wp, &
            i1=int(0.9_wp*sim%nx)+1, i2=sim%nx, &
            j1=int(0.9_wp*sim%ny)+1, j2=sim%ny)
        ! Warm start: initialise fuel at inlet temperature so the first SS
        ! iteration sees the correct thermal-hydraulic baseline rather than
        ! room temperature, which would give a spuriously cold avg_coolant_temp
        ! and therefore a negative reactor level on the first vessel tick.
        sim%heat%T     = sim%inlet_temperature
        sim%heat%T_old = sim%inlet_temperature

        ! ── Two-phase thermalhydraulics ────────────────────────────────────────
        th_config%void_correlation          = VOID_CHEXAL_LELLOUCHE_ID
        th_config%chf_correlation           = CHF_GROENEVELD_ID
        th_config%include_subcooled_boiling = .true.
        call two_phase_init(sim%thermalhydraulics, sim%nx, sim%ny, sim%nz, &
                            sim%dx, sim%dy, sim%dz, th_config)
        call two_phase_set_geometry(sim%thermalhydraulics, diameter=0.012_wp)

        ! ── Burnup / depletion ─────────────────────────────────────────────────
        burnup_config%track_xenon     = .true.
        burnup_config%track_samarium  = .true.
        burnup_config%track_actinides = .true.
        call burnup_init(sim%burnup, sim%nx, sim%ny, sim%nz, &
                         sim%dx, sim%dy, sim%dz, burnup_config)
        call burnup_set_initial_composition(sim%burnup, enrichment=enrichment)

        ! Prime neutronics XS with the initial rod/temperature state so that the
        ! first SS eigenvalue solve sees rod absorption rather than bare fuel XS
        ! (k_inf≈1.75). Without this, Chebyshev acceleration drives keff to the
        ! upper clamp before any rod effects are applied.
        block
            real(wp) :: T_init(sim%nx, sim%ny, sim%nz)
            real(wp) :: rho_init(sim%nx, sim%ny, sim%nz)
            T_init   = sim%inlet_temperature
            rho_init = 0.738_wp  ! nominal liquid density, g/cm³
            call update_cross_sections_feedback(sim, T_init, rho_init)
        end block
    end subroutine init_physics

    ! ── Steady-state solve ────────────────────────────────────────────────────

    subroutine do_solve_steady_state(sim)
        type(bwr_state_t), intent(inout) :: sim

        integer  :: iter
        real(wp) :: power(sim%nx, sim%ny, sim%nz)
        real(wp) :: temperature(sim%nx, sim%ny, sim%nz)
        real(wp) :: void_fraction(sim%nx, sim%ny, sim%nz)
        real(wp) :: density(sim%nx, sim%ny, sim%nz)
        real(wp) :: error
        real(wp), parameter :: rho_liquid = 0.738_wp
        real(wp), parameter :: rho_vapor  = 0.038_wp
        ! Cold-startup guard: dome pressure well below operating means the
        ! reactor is in cold standby (no boiling, no power). Calling the SS
        ! solver here is meaningless — there is no power-rated equilibrium
        ! to converge to while rods are inserted, so the loop would burn
        ! its iteration budget for nothing. 1 MPa is a clean divider:
        ! cold preset boots at 1 atm; hot operating dome sits at ~7 MPa.
        real(wp), parameter :: COLD_SS_THRESHOLD_PA = 1.0e6_wp
        logical  :: converged

        if (sim%vessel%pressure_operating < COLD_SS_THRESHOLD_PA) then
            print '(A,F8.3,A)', "  SS skipped: vessel at ", &
                sim%vessel%pressure_operating / 1.0e5_wp, &
                " bar (cold standby). Use bwr_step to drive the startup transient."
            ! Leave state untouched: heat%T is cold, flux/power are zero,
            ! precursors are zero, k_eff is whatever init_physics left it.
            ! The transient solver in do_step will pick up from here.
            return
        end if

        ! Diagnostic: show center-cell XS before first eigenvalue solve to verify
        ! rod-insertion priming. Should show sigma_a(2)>>0.2 if 95% rods applied.
        block
            integer :: ci, cj, ck
            ci = sim%nx / 2;  cj = sim%ny / 2;  ck = sim%nz / 2
            print '(A,F8.4,A,F8.4,A,F8.4,A,F8.4)', &
                "  SS pre-solve  sigma_a2=", sim%neutronics%xsec(ci,cj,ck)%sigma_a(2), &
                "  nu_sf2=",  sim%neutronics%xsec(ci,cj,ck)%nu_sigma_f(2), &
                "  k_eff0=",  sim%neutronics%k_eff, &
                "  blade1=",  sim%crd%blade_insertion(1)
        end block

        ! Disable flux normalization during SS — it forces flux to rated power even
        ! for deeply subcritical startups, corrupting the eigenvalue convergence.
        sim%neutronics%config%normalize_power = .false.

        do iter = 1, 50
            call mg_solve_eigenvalue(sim%neutronics, sim%k_eff, converged)
            if (iter <= 3) print '(A,I3,A,F7.4,A,L1)', &
                "  SS iter ", iter, "  k_eff=", sim%k_eff, "  conv=", converged
            if (.not. converged .and. iter > 3) &
                print '(A,I3)', "  Warning: neutronics did not converge, iter ", iter

            call mg_get_power(sim%neutronics, power)
            sim%power_current = sim%neutronics%total_power
            if (sim%k_eff < 0.98_wp) then
                sim%power_current = 0.0_wp
                power = 0.0_wp
            end if
            print '(A,I3,A,F7.4,A,F8.1,A,L1)', &
                "  SS iter ", iter, &
                "  k_eff=", sim%k_eff, &
                "  Tmax_K=", sim%max_fuel_temp, &
                "  gated=", sim%k_eff < 0.98_wp

            ! mg_get_power → W/cm³; heat%q needs W/m³. Reorg step 10
            ! pinned this conversion into `fuel_power_to_volumetric_W_m3`
            ! so the unit contract is checked in unit tests.
            sim%heat%q = fuel_power_to_volumetric_W_m3(power)

            ! Step-11: feed per-cell T_fluid + h_conv from the previous
            ! iteration's two-phase state into the heat kernel. First
            ! iteration reads from the defaults set by two_phase_init
            ! (sat at 7 MPa, mass flux 1000 kg/m²·s, no boiling → pure
            ! Dittus-Boelter). Later iterations reflect the actual
            ! regime mix as the steady-state loop converges.
            block
                real(wp) :: T_f(sim%nx, sim%ny, sim%nz)
                real(wp) :: h_c(sim%nx, sim%ny, sim%nz)
                call fuel_compute_convection_coefficients(sim%thermalhydraulics, T_f, h_c, &
                                                          T_bulk=sim%avg_coolant_temp)
                call heat_set_coolant(sim%heat, T_f, h_c)
            end block

            call heat_step_implicit(sim%heat, 1.0_wp)
            temperature       = sim%heat%T
            sim%max_fuel_temp = maxval(temperature)

            sim%avg_coolant_temp = avg_coolant_temp_K(temperature, sim%fuel%is_fuel_cell)

            call two_phase_step(sim%thermalhydraulics, &
                temperature, &
                sim%vessel%pressure_operating + 0.0_wp * temperature, &
                sim%recirc%mass_flux_kg_m2_s + 0.0_wp * temperature, &
                fuel_power_to_surface_q(power, sim%dx, sim%dy, &
                                        sim%thermalhydraulics%heated_perimeter), &
                1.0_wp)

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

        ! Subcritical startups: the eigenvalue loop may have driven flux to
        ! the MAX_FLUX=1e15 clamp because normalize_flux skips power
        ! scaling below k=0.95. Renormalize so the transient solver
        ! starts from a sane amplitude. For near-critical SS exits the
        ! eigenvalue solver already scaled flux to rated power; leave it.
        if (sim%k_eff < 0.95_wp) then
            call normalize_flux_to_unity(sim%neutronics)
        end if

        sim%neutronics%config%normalize_power = .true.

        ! Seed delayed-neutron precursors at equilibrium
        ! C_d = beta_d * F / lambda_d, with F = Sum_g nu_sigma_f,g * phi_g
        ! (the delayed-neutron production rate; mirrors the per-cell
        ! fission_total used by mg_solve_transient). Without this seed
        ! the first transient steps would see a delayed source of zero
        ! and the prompt-only flux would slump on a 1/v timescale.
        block
            real(wp), parameter :: lambda(6) = [0.0124_wp, 0.0305_wp, 0.111_wp, &
                                                0.301_wp,  1.14_wp,   3.01_wp]
            real(wp), parameter :: beta(6)   = [0.000215_wp, 0.001424_wp, 0.001274_wp, &
                                                0.002568_wp, 0.000748_wp, 0.000273_wp]
            integer  :: i, j, k, d, g
            real(wp) :: F_cell
            do k = 1, sim%nz
                do j = 1, sim%ny
                    do i = 1, sim%nx
                        F_cell = 0.0_wp
                        do g = 1, sim%neutronics%n_groups
                            F_cell = F_cell + &
                                sim%neutronics%xsec(i,j,k)%nu_sigma_f(g) * &
                                sim%neutronics%flux(i,j,k,g)
                        end do
                        do d = 1, 6
                            sim%precursors(i,j,k,d) = beta(d) * F_cell / lambda(d)
                        end do
                    end do
                end do
            end do
        end block

        print '(A,F8.1,A,F8.1)', &
            "  SS done: max_fuel_temp_K=", sim%max_fuel_temp, &
            "  heat_T_max_K=", maxval(sim%heat%T)

        sim%reactivity_pcm = (sim%k_eff - 1.0_wp) / sim%k_eff * 1.0e5_wp

        if (allocated(sim%thermalhydraulics%chf_ratio)) then
            sim%min_chfr = minval(sim%thermalhydraulics%chf_ratio)
        else
            sim%min_chfr = 999.0_wp
        end if
        if (sim%power_current <= 0.0_wp) sim%min_chfr = 999.0_wp

        call update_instrumentation(sim)
        call update_pressure_dynamics(sim, 0.0_wp)
        call update_main_steam(sim, 0.0_wp)
        call update_feedwater(sim, 0.0_wp)
        call update_rhr(sim, 0.0_wp)
        call update_condensate_loop(sim, 0.0_wp)
    end subroutine do_solve_steady_state

    ! ── Coupled time step ─────────────────────────────────────────────────────

    subroutine do_step(sim, dt)
        type(bwr_state_t), intent(inout) :: sim
        real(wp), intent(in) :: dt

        real(wp) :: power(sim%nx, sim%ny, sim%nz)
        real(wp) :: flux(sim%nx, sim%ny, sim%nz)
        real(wp) :: temperature(sim%nx, sim%ny, sim%nz)
        real(wp) :: void_fraction(sim%nx, sim%ny, sim%nz)
        real(wp) :: density(sim%nx, sim%ny, sim%nz)
        ! rho_liquid_vel: kg/m³ — used for coolant velocity (mass_flux / density)
        ! rho_liquid_xs / rho_vapor_xs: g/cm³ — expected by fuel_apply_to_xs
        !   (must match the units used in do_solve_steady_state)
        real(wp), parameter :: rho_liquid_vel = 738.0_wp
        real(wp), parameter :: rho_liquid_xs  = 0.738_wp
        real(wp), parameter :: rho_vapor_xs   = 0.038_wp
        logical  :: converged
        real(wp) :: v_coolant

        type(crd_drivers_t)    :: crd_drivers
        type(recirc_drivers_t) :: recirc_drivers

        ! Tick order per REORG_PLAN §5: crd first — rod motion / scram
        ! dynamics update the blade insertions before cross-sections are
        ! rebuilt in update_cross_sections_feedback. Uses last-step's dome
        ! pressure (vessel hasn't run this tick yet); on a 20 ms tick at
        ! ~1000 psig this lag is well below the scram-curve sensitivity.
        crd_drivers%dome_pressure_Pa = sim%vessel%pressure_operating
        crd_drivers%reactor_period_s = sim%reactor_period
        ! k_eff from last tick's solver. CRD uses it for the approach-to-
        ! critical throttle (slew rate ∝ 1-k) and the hard-stop above
        ! rod_block_keff. One-tick lag is acceptable because the throttle
        ! already slows the motion by 10-25× near critical.
        crd_drivers%k_eff = sim%k_eff
        call crd_step(sim%crd, dt, crd_drivers, sim%crd_config)

        ! Recirculation runs before two_phase so the kernel sees the
        ! freshly-stepped core mass flux. NPSH protection reads last-tick
        ! FW flow + reactor level, and EOC-RPT consumes last-tick
        ! turbine_tripped. The cross-tick lag is harmless at 50 Hz.
        recirc_drivers%feedwater_flow_pct = sim%fw%flow_pct
        recirc_drivers%reactor_level_m    = sim%vessel%Lrx
        recirc_drivers%eoc_rpt_signal     = sim%ms%turbine_tripped
        call recirc_step(sim%recirc, dt, recirc_drivers, sim%recirc_config)

        ! Disable the heat solver's central-difference advection term.
        ! That term is `-v · ∇T` and is unconditionally unstable when the
        ! gradient is computed with central differences in an explicit
        ! Euler step (which is what heat_step does). At cold standby any
        ! tiny temperature non-uniformity from the source field amplified
        ! over ~70 s into a localised cell shooting up to 320 K+ — the
        ! "fuel goes infinite at cold standby" symptom.
        ! The per-cell convection sink q_conv = h · (T − T_fluid) already
        ! encodes the actual heat removal physics correctly and is stable;
        ! the gridded advection is redundant for the lumped sim. v_coolant
        ! is still computed so future use cases (e.g. velocity-dependent
        ! h_conv via Dittus-Boelter) can pick it up, but not pushed to
        ! the unstable solver path.
        v_coolant = sim%recirc%mass_flux_kg_m2_s / rho_liquid_vel
        if (allocated(sim%heat%vz)) sim%heat%vz = 0.0_wp
        if (allocated(sim%heat%vx)) sim%heat%vx = 0.0_wp
        if (allocated(sim%heat%vy)) sim%heat%vy = 0.0_wp

        ! Transient multigroup diffusion with 6-group delayed-neutron precursors.
        ! Replaces the per-step steady-state eigenvalue solve: instead of
        ! re-converging M*phi = (1/k)*F*phi from scratch each tick (10-50
        ! outer x 5-15 inner sweeps), we advance the flux + precursors one
        ! step with Crank-Nicolson on a warm-start flux, which converges in
        ! 1-3 iterations for the small XS perturbations between ticks.
        ! `converged` and `sim%k_eff` are updated diagnostically inside the
        ! solver (instantaneous k from production/absorption ratio).
        call mg_solve_transient(sim%neutronics, sim%precursors, dt)
        converged = .true.
        sim%k_eff = sim%neutronics%k_eff
        call mg_get_power(sim%neutronics, power)
        call mg_get_flux(sim%neutronics, flux)
        sim%power_current = sim%neutronics%total_power
        ! NOTE: the old "if k < 0.98 → zero power" gate is gone. It was a
        ! kludge from the eigenvalue-solver era when the solver produced
        ! spurious subcritical flux. mg_solve_transient + precursor coupling
        ! handles the subcritical regime correctly — source-multiplied flux
        ! at k=0.96 is ~25× the precursor decay floor, and the operator
        ! needs to *see* APRM / period evolving during approach to critical.
        ! Zeroing power here was the reason APRM stuck at 0 % with rods at
        ! 45 % withdrawn and the period readout pinned at the 1e6 s default.
        if (sim%n_steps < 5 .or. sim%max_fuel_temp > 1000.0_wp) then
            print '(A,I6,A,F7.4,A,F10.1,A,F10.1)', &
                "  step ", sim%n_steps, &
                "  k=", sim%k_eff, &
                "  Tmax_K=", sim%max_fuel_temp, &
                "  Q_max=", maxval(sim%heat%q)
        end if

        call burnup_step(sim%burnup, flux, power, dt)
        sim%avg_burnup = sum(sim%burnup%burnup) / size(sim%burnup%burnup)

        ! mg_get_power → W/cm³; heat%q needs W/m³. See note in
        ! do_solve_steady_state — `fuel_power_to_volumetric_W_m3` pins
        ! the contract.
        sim%heat%q = fuel_power_to_volumetric_W_m3(power)

        ! Step-11: per-cell T_fluid / h_conv from last-tick two-phase
        ! state. See do_solve_steady_state for the rationale; for the
        ! transient path the cross-tick lag is harmless at 50 Hz.
        block
            real(wp) :: T_f(sim%nx, sim%ny, sim%nz)
            real(wp) :: h_c(sim%nx, sim%ny, sim%nz)
            ! T_bulk = last-tick avg_coolant_temp. One-tick lag at 50 Hz
            ! is negligible. The convection helper clamps to T_sat so a
            ! fuel-dominated avg can never push T_fluid above saturation.
            call fuel_compute_convection_coefficients(sim%thermalhydraulics, T_f, h_c, &
                                                      T_bulk=sim%avg_coolant_temp)
            call heat_set_coolant(sim%heat, T_f, h_c)
            if (mod(sim%n_steps, 50) == 0 .or. sim%max_fuel_temp > 293.5_wp) then
                print '(A,I6,A,F8.2,A,F8.2,A,ES10.2,A,ES10.2,A,ES10.2)', &
                    "  DBG step ", sim%n_steps, &
                    "  Tfuel_max=", sim%max_fuel_temp, &
                    "  T_bulk=", sim%avg_coolant_temp, &
                    "  P_W=", sim%power_current, &
                    "  Q_max=", maxval(sim%heat%q), &
                    "  k=", sim%k_eff
            end if
        end block

        call heat_step(sim%heat, dt)
        temperature       = sim%heat%T
        sim%max_fuel_temp = maxval(temperature)

        sim%avg_coolant_temp = avg_coolant_temp_K(temperature, sim%fuel%is_fuel_cell)

        call two_phase_step(sim%thermalhydraulics, &
            temperature, &
            sim%vessel%pressure_operating + 0.0_wp * temperature, &
            sim%recirc%mass_flux_kg_m2_s + 0.0_wp * temperature, &
            fuel_power_to_surface_q(power, sim%dx, sim%dy, &
                                    sim%thermalhydraulics%heated_perimeter), &
            dt)

        call two_phase_get_void_fraction(sim%thermalhydraulics, void_fraction)

        density = (1.0_wp - void_fraction) * rho_liquid_xs + void_fraction * rho_vapor_xs
        density = max(density, rho_vapor_xs)
        density = min(density, rho_liquid_xs)

        sim%avg_void_fraction = sum(void_fraction) / size(void_fraction) * 100.0_wp
        sim%max_void_fraction = maxval(void_fraction) * 100.0_wp

        if (allocated(sim%thermalhydraulics%chf_ratio)) then
            sim%min_chfr = minval(sim%thermalhydraulics%chf_ratio)
        else
            sim%min_chfr = 999.0_wp
        end if
        ! CHF is undefined at zero heat flux — prevent false trips during subcritical operation.
        if (sim%power_current <= 0.0_wp) sim%min_chfr = 999.0_wp

        call update_cross_sections_feedback(sim, temperature, density)
        call update_instrumentation(sim)
        call update_pressure_dynamics(sim, dt)
        call update_main_steam(sim, dt)
        call update_feedwater(sim, dt)
        call update_rhr(sim, dt)
        call update_condensate_loop(sim, dt)

        sim%reactivity_pcm = (sim%k_eff - 1.0_wp) / sim%k_eff * 1.0e5_wp
        sim%time    = sim%time + dt
        sim%n_steps = sim%n_steps + 1
    end subroutine do_step

    ! ══════════════════════════════════════════════════════════════════════════
    !  C-BOUND FUNCTIONS
    ! ══════════════════════════════════════════════════════════════════════════

    ! ── Lifecycle ─────────────────────────────────────────────────────────────

    function bwr_reactor_create(nx, ny, nz, dx, dy, dz, &
                                 power_rated_W, inlet_temp_K, pressure_Pa, &
                                 mass_flux, enrichment, &
                                 core_height_m, core_diameter_m) &
                                 result(handle) bind(C, name="bwr_reactor_create")
        integer(c_int),  value :: nx, ny, nz
        real(c_double),  value :: dx, dy, dz
        real(c_double),  value :: power_rated_W, inlet_temp_K, pressure_Pa
        real(c_double),  value :: mass_flux, enrichment
        real(c_double),  value :: core_height_m, core_diameter_m
        type(c_ptr) :: handle

        type(bwr_state_t), pointer :: sim

        ! ── Cold-startup preset constants ──────────────────────────────────────
        ! The reactor boots in cold standby. Rated values from the C ABI
        ! (pressure_Pa, inlet_temp_K, mass_flux, power_rated_W) are kept
        ! as setpoints and references; only the *initial state* is cold.
        real(wp), parameter :: COLD_T_K     = 293.15_wp     ! ~20 °C ambient
        real(wp), parameter :: ATM_PA       = 1.01325e5_wp  ! sea-level atm
        real(wp), parameter :: COLD_RHO_GCC = 0.998_wp      ! water @20 °C, g/cm³
        real(wp), parameter :: FW_COLD_T_K  = 319.0_wp      ! hotwell suction (~46 °C)

        call rng_seed(123456789_8)
        allocate(sim)

        sim%nx = nx;  sim%ny = ny;  sim%nz = nz
        sim%dx = dx;  sim%dy = dy;  sim%dz = dz
        sim%core_height   = core_height_m
        sim%core_diameter = core_diameter_m

        sim%power_rated        = power_rated_W
        ! sim%inlet_temperature is only used inside init_physics to seed
        ! heat%T and the XS prime, *not* as a rated reference (rated FW
        ! inlet temp lives in vessel_config%rated_fw_inlet_temp_K and
        ! fw_config%rated_inlet_temperature_K). Override to cold standby
        ! so the heat field and XS feedback boot consistent with the rest
        ! of the cold preset.
        sim%inlet_temperature  = COLD_T_K

        ! Seed the fuel subsystem geometry: builds the 137-blade →
        ! per-cell index map and caches dx/dy/dz so the XS-feedback
        ! path doesn't need to redrive geometry calls each tick.
        call fuel_init_geometry(sim%fuel, sim%fuel_config, sim%nx, sim%ny, &
                                sim%dx, sim%dy, sim%dz, sim%core_diameter)

        ! Seed the CRD subsystem. Bank starts at 95 % inserted (reactor
        ! startup configuration), matching the pre-reorg rod_bank_position
        ! default. Active-fuel height is taken from the caller-supplied
        ! core_height_m so the notch tracking stays in sync.
        sim%crd_config%core_height_m     = core_height_m
        sim%crd_config%initial_insertion = 0.95_wp
        call crd_init(sim%crd, sim%crd_config)

        ! Seed the vessel subsystem. Rated values track the C ABI so
        ! downstream setpoint logic (P_NOMINAL ratios, sat_temp scaling)
        ! still references the rated operating point. The vessel state
        ! itself is then stomped to cold standby below.
        sim%vessel_config%rated_pressure_Pa     = pressure_Pa
        sim%vessel_config%rated_fw_inlet_temp_K = inlet_temp_K
        sim%vessel_config%core_height_m         = core_height_m
        call vessel_init(sim%vessel, sim%vessel_config)

        ! Cold vessel state: atmospheric dome, ambient coolant. The ODE
        ! history (P_prev, coolant_T_prev) is seeded coherent with the
        ! initial state so the first vessel_step_dynamics doesn't see a
        ! phantom transient. Lbase is recomputed with cold-water density
        ! so Lrx lands on the right level given the heavier coolant.
        sim%vessel%pressure_operating = ATM_PA
        sim%vessel%P_prev1            = ATM_PA
        sim%vessel%P_prev2            = ATM_PA
        sim%vessel%sat_temperature    = sat_temp_K(ATM_PA)
        sim%vessel%T_excess           = 0.0_wp
        sim%vessel%coolant_T_prev1    = COLD_T_K
        sim%vessel%coolant_T_prev2    = COLD_T_K
        sim%vessel%Lbase = (0.5_wp + core_height_m) &
                        * (coolant_density_gl(COLD_T_K) / 1000.0_wp)
        sim%vessel%Lrx   = calc_reactor_level(sim%vessel%Lbase, COLD_T_K, core_height_m)

        sim%avg_coolant_temp = COLD_T_K

        sim%neutrons_prev  = 0.0_wp
        sim%reactor_period = 1.0e6_wp
        sim%aprm           = 0.0_wp

        sim%reactivity_perturbation_pcm = 0.0_wp

        ! Main steam in cold preset: turbine isolated (0 % valve), bypass
        ! fully open (100 %). Until the dome boils there is no steam to
        ! pass; once boiling starts, the bypass routes all of it to the
        ! condenser so the turbine stays warm and idle.
        sim%ms_config%rated_pressure_Pa         = pressure_Pa
        sim%ms_config%initial_turbine_valve_pct = 0.0_wp
        sim%ms_config%initial_bypass_valve_pct  = 100.0_wp
        call main_steam_init(sim%ms, sim%ms_config)

        ! Feedwater in cold preset: zero flow. An unboiling reactor does
        ! not consume coolant; running rated FW would fill the vessel
        ! solid in seconds. The 3-element controller stays disabled (its
        ! level-error term would otherwise demand flow against an already
        ! full vessel).
        sim%fw_config%rated_inlet_temperature_K = inlet_temp_K
        call feedwater_init(sim%fw, sim%fw_config)
        sim%fw%demand_flow_pct  = 0.0_wp
        sim%fw%flow_pct         = 0.0_wp
        sim%fw%flow_kg_s        = 0.0_wp
        sim%fw%rfpt_speed_pct   = 0.0_wp
        sim%fw%feedwater_temp_K = FW_COLD_T_K
        sim%fw%controller_enabled = .false.

        ! Recirc pumps idle in cold preset. The pumps' purpose is to
        ! redistribute steam voids in the downcomer/jet-pump loop; with
        ! no boiling there are no voids to manage. Natural circulation
        ! provides the floor flow via recirc%nat_circ_fraction inside
        ! recirc_step.
        sim%recirc_config%rated_mass_flux_kg_m2_s = real(mass_flux, wp)
        sim%recirc_config%initial_pump_speed_pct  = 0.0_wp
        call recirc_init(sim%recirc, sim%recirc_config)

        ! RHR in cold preset: all four pumps stopped, mode STANDBY.
        ! Operator brings RHR online as needed (SDC for cooldown, pool
        ! cooling if pool heats up, LPCI on demand). Suppression pool at
        ! 305 K (≈ 90 °F design intake) — heats up only via SRV discharge.
        call rhr_init(sim%rhr, sim%rhr_config)
        sim%supp_pool_T_K       = sim%rhr_config%service_water_T_K
        sim%drywell_pressure_Pa = ATM_PA

        ! Condensate loop: hotwell at normal level, CST near full, one
        ! transfer pump running per spec.
        call cond_init(sim%cond, sim%cond_config)

        sim%power_current     = 0.0_wp
        sim%k_eff             = 1.0_wp
        sim%reactivity_pcm    = 0.0_wp
        sim%max_fuel_temp     = COLD_T_K
        sim%avg_void_fraction = 0.0_wp
        sim%max_void_fraction = 0.0_wp
        sim%min_chfr          = 999.0_wp
        sim%avg_burnup        = 0.0_wp
        sim%time              = 0.0_wp
        sim%n_steps           = 0

        call init_physics(sim, real(enrichment, wp))

        ! Re-prime XS feedback with cold-water density. init_physics
        ! seeds the XS at sim%inlet_temperature (now cold) but with a
        ! hardcoded rho=0.738 g/cm³ that reflects hot-channel density.
        ! Cold liquid is denser (~0.998 g/cm³), which strengthens
        ! moderation and shifts the cold-shutdown k slightly.
        block
            real(wp) :: T_cold(sim%nx, sim%ny, sim%nz)
            real(wp) :: rho_cold(sim%nx, sim%ny, sim%nz)
            T_cold   = COLD_T_K
            rho_cold = COLD_RHO_GCC
            call update_cross_sections_feedback(sim, T_cold, rho_cold)
        end block

        ! Two-phase model defaults to 7 MPa / 558 K (sat at rated). On
        ! the first do_step tick, fuel_compute_convection_coefficients
        ! runs *before* two_phase_step, so the heat solver would read
        ! the stale defaults and convect the fuel to 558 K in one step,
        ! kicking off a vessel-pressure runaway. Stomp the cached
        ! pressure / temperature / props field with cold values so the
        ! first heat_step sees a sub-saturation T_fluid.
        block
            type(water_properties_t) :: w_cold
            w_cold = get_water_properties(COLD_T_K, ATM_PA)
            sim%thermalhydraulics%pressure    = ATM_PA
            sim%thermalhydraulics%temperature = COLD_T_K
            sim%thermalhydraulics%props       = w_cold
            sim%thermalhydraulics%void_fraction = 0.0_wp
            sim%thermalhydraulics%quality       = 0.0_wp
            sim%thermalhydraulics%boiling       = .false.
        end block

        handle = c_loc(sim)
    end function bwr_reactor_create

    subroutine bwr_reactor_destroy(handle) bind(C, name="bwr_reactor_destroy")
        type(c_ptr), value :: handle
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        call mg_destroy(sim%neutronics)
        call heat_destroy(sim%heat)
        call two_phase_destroy(sim%thermalhydraulics)
        call burnup_destroy(sim%burnup)
        call xslib_destroy(sim%xslib)
        call fuel_destroy(sim%fuel)
        call rhr_destroy(sim%rhr)
        call cond_destroy(sim%cond)
        if (allocated(sim%precursors)) deallocate(sim%precursors)
        deallocate(sim)
    end subroutine bwr_reactor_destroy

    ! ── Debug / diagnostics ──────────────────────────────────────────────────

    ! Returns thermal (group-2) absorption XS at the core center and the
    ! mg_state k_eff before the first SS solve — used to verify XS priming
    ! and Chebyshev-disable fixes are taking effect.
    subroutine bwr_debug_xs(handle, sigma_a2, mg_keff) &
            bind(C, name="bwr_debug_xs")
        type(c_ptr),         value    :: handle
        real(c_double), intent(out)   :: sigma_a2, mg_keff
        type(bwr_state_t), pointer    :: sim
        integer :: ci, cj, ck
        call c_f_pointer(handle, sim)
        ci = sim%nx / 2;  cj = sim%ny / 2;  ck = sim%nz / 2
        sigma_a2 = real(sim%neutronics%xsec(ci, cj, ck)%sigma_a(2), c_double)
        mg_keff  = real(sim%neutronics%k_eff, c_double)
    end subroutine bwr_debug_xs

    ! ── Simulation control ────────────────────────────────────────────────────

    subroutine bwr_solve_steady_state(handle) bind(C, name="bwr_solve_steady_state")
        type(c_ptr), value :: handle
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        call do_solve_steady_state(sim)
    end subroutine bwr_solve_steady_state

    subroutine bwr_step(handle, dt) bind(C, name="bwr_step")
        type(c_ptr),    value :: handle
        real(c_double), value :: dt
        type(bwr_state_t), pointer :: sim
        real(wp) :: dt_clamped
        real(wp), parameter :: DT_MAX = 0.1_wp     ! 5× recommended 0.02 s; beyond this the implicit solvers risk divergence.
        real(wp), parameter :: DT_MIN = 1.0e-6_wp  ! 1 µs floor; smaller is numerical noise.

        call c_f_pointer(handle, sim)
        if (.not. (dt > 0.0_wp)) return            ! skip on NaN / <= 0
        dt_clamped = min(real(dt, wp), DT_MAX)
        dt_clamped = max(dt_clamped, DT_MIN)
        call do_step(sim, dt_clamped)
    end subroutine bwr_step

    ! ── Operator controls ─────────────────────────────────────────────────────

    subroutine bwr_set_control_rod_position(handle, insertion) &
            bind(C, name="bwr_set_control_rod_position")
        ! Broadcasts the requested insertion fraction to all 137 blades.
        ! This is the legacy single-knob rod-bank handle preserved for C-ABI
        ! compatibility. For a real scram with the Figure 2.3-13 timing,
        ! the caller should use bwr_scram() instead.
        type(c_ptr),    value :: handle
        real(c_double), value :: insertion
        type(bwr_state_t), pointer :: sim
        type(crd_command_t) :: cmd
        call c_f_pointer(handle, sim)
        cmd%bank_position_set = max(0.0_wp, min(1.0_wp, real(insertion, wp)))
        call crd_apply(sim%crd, cmd)
    end subroutine bwr_set_control_rod_position

    subroutine bwr_scram(handle) bind(C, name="bwr_scram")
        ! Trigger reactor protection system scram. The 137 control blades
        ! begin inserting following the "typical drive" envelope from
        ! spec 2.3 Figure 2.3-13 (time vs dome pressure). Idempotent —
        ! a second call while already scrammed is a no-op.
        type(c_ptr), value :: handle
        type(bwr_state_t), pointer :: sim
        type(crd_command_t) :: cmd
        call c_f_pointer(handle, sim)
        cmd%scram_latch = .true.
        call crd_apply(sim%crd, cmd)
    end subroutine bwr_scram

    subroutine bwr_scram_reset(handle) bind(C, name="bwr_scram_reset")
        ! Clear the latched scram condition. Does not retract the blades;
        ! operator must withdraw them via bwr_set_control_rod_position.
        type(c_ptr), value :: handle
        type(bwr_state_t), pointer :: sim
        type(crd_command_t) :: cmd
        call c_f_pointer(handle, sim)
        cmd%scram_reset = .true.
        call crd_apply(sim%crd, cmd)
    end subroutine bwr_scram_reset

    ! Legacy core mass-flux dial [kg/m²·s]. Now routes to the recirc pump
    ! speed setpoint: kg_m2_s is divided by rated_mass_flux to derive a
    ! 0..max_pump_speed_pct demand applied to both pumps. NPSH protection
    ! and the M-ratio model then determine the actual core mass flux —
    ! so the value the kernel sees may not match the operator's number
    ! during transients. For direct pump-speed control prefer the new
    ! `bwr_set_recirc_pump_speed_pct` API.
    subroutine bwr_set_mass_flux(handle, kg_m2_s) bind(C, name="bwr_set_mass_flux")
        type(c_ptr),    value :: handle
        real(c_double), value :: kg_m2_s
        type(bwr_state_t), pointer :: sim
        type(recirc_command_t) :: cmd
        real(wp) :: speed_pct
        call c_f_pointer(handle, sim)
        speed_pct = real(kg_m2_s, wp) / max(1.0_wp, sim%recirc_config%rated_mass_flux_kg_m2_s) &
                  * 100.0_wp
        cmd%pump_speed_set_both = max(0.0_wp, &
                                      min(sim%recirc_config%max_pump_speed_pct, speed_pct))
        call recirc_apply(sim%recirc, cmd, sim%recirc_config)
    end subroutine bwr_set_mass_flux

    subroutine bwr_set_recirc_pump_speed_pct(handle, pct) &
            bind(C, name="bwr_set_recirc_pump_speed_pct")
        ! Set both recirculation pump speed demands (in % of rated speed).
        ! Spec 2.4 p.3 operating range 30 – 102 %. Values outside are
        ! clamped; below the turndown threshold the pump-speed dynamics
        ! collapse to zero.
        type(c_ptr),    value :: handle
        real(c_double), value :: pct
        type(bwr_state_t), pointer :: sim
        type(recirc_command_t) :: cmd
        call c_f_pointer(handle, sim)
        cmd%pump_speed_set_both = real(pct, wp)
        call recirc_apply(sim%recirc, cmd, sim%recirc_config)
    end subroutine bwr_set_recirc_pump_speed_pct

    subroutine bwr_trip_recirc_pump(handle, idx) bind(C, name="bwr_trip_recirc_pump")
        ! Trip a specific recirc pump (idx ∈ {1, 2}). The pump coasts to
        ! zero speed; the unaffected loop keeps running. Both pumps off
        ! drops core flow to the natural-circulation floor (~25 % rated).
        type(c_ptr),    value :: handle
        integer(c_int), value :: idx
        type(bwr_state_t), pointer :: sim
        type(recirc_command_t) :: cmd
        call c_f_pointer(handle, sim)
        cmd%trip_pump_idx = int(idx)
        call recirc_apply(sim%recirc, cmd, sim%recirc_config)
    end subroutine bwr_trip_recirc_pump

    subroutine bwr_reset_recirc_pumps(handle) bind(C, name="bwr_reset_recirc_pumps")
        type(c_ptr), value :: handle
        type(bwr_state_t), pointer :: sim
        type(recirc_command_t) :: cmd
        call c_f_pointer(handle, sim)
        cmd%reset_all_pumps = .true.
        cmd%eoc_rpt_reset   = .true.
        call recirc_apply(sim%recirc, cmd, sim%recirc_config)
    end subroutine bwr_reset_recirc_pumps

    ! Feedwater flow [% of rated steam flow]. The reactor water level holds at
    ! rated power when this is 100 %; below drains the level, above floods it.
    ! Clamped to [0, 200].
    subroutine bwr_set_feedwater_flow_pct(handle, pct) &
            bind(C, name="bwr_set_feedwater_flow_pct")
        type(c_ptr),    value :: handle
        real(c_double), value :: pct
        type(bwr_state_t), pointer :: sim
        type(feedwater_command_t) :: cmd
        call c_f_pointer(handle, sim)
        cmd%demand_flow_pct_set = max(0.0_wp, min(200.0_wp, real(pct, wp)))
        call feedwater_apply(sim%fw, cmd)
    end subroutine bwr_set_feedwater_flow_pct

    function bwr_get_feedwater_flow_pct(handle) result(v) &
            bind(C, name="bwr_get_feedwater_flow_pct")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%fw%flow_pct, c_double)
    end function bwr_get_feedwater_flow_pct

    function bwr_get_feedwater_temp_K(handle) result(v) &
            bind(C, name="bwr_get_feedwater_temp_K")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%fw%feedwater_temp_K, c_double)
    end function bwr_get_feedwater_temp_K

    function bwr_get_rfp_running_count(handle) result(v) &
            bind(C, name="bwr_get_rfp_running_count")
        type(c_ptr), value :: handle
        integer(c_int) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = int(count(.not. sim%fw%rfp_tripped), c_int)
    end function bwr_get_rfp_running_count

    subroutine bwr_trip_feedwater_pump(handle, idx) bind(C, name="bwr_trip_feedwater_pump")
        ! Trip a specific RFP (1..2). Idempotent — re-tripping an already
        ! tripped pump is a no-op. Out-of-range indices are ignored.
        type(c_ptr),    value :: handle
        integer(c_int), value :: idx
        type(bwr_state_t), pointer :: sim
        type(feedwater_command_t) :: cmd
        call c_f_pointer(handle, sim)
        cmd%trip_rfp_idx = int(idx)
        call feedwater_apply(sim%fw, cmd)
    end subroutine bwr_trip_feedwater_pump

    subroutine bwr_reset_feedwater_pumps(handle) bind(C, name="bwr_reset_feedwater_pumps")
        type(c_ptr), value :: handle
        type(bwr_state_t), pointer :: sim
        type(feedwater_command_t) :: cmd
        call c_f_pointer(handle, sim)
        cmd%reset_all_pumps = .true.
        call feedwater_apply(sim%fw, cmd)
    end subroutine bwr_reset_feedwater_pumps

    subroutine bwr_set_turbine_valve(handle, pct) bind(C, name="bwr_set_turbine_valve")
        type(c_ptr),    value :: handle
        real(c_double), value :: pct
        type(bwr_state_t), pointer :: sim
        type(main_steam_command_t) :: cmd
        call c_f_pointer(handle, sim)
        cmd%turbine_valve_set_pct = max(0.0_wp, min(100.0_wp, real(pct, wp)))
        call main_steam_apply(sim%ms, cmd, sim%ms_config)
    end subroutine bwr_set_turbine_valve

    subroutine bwr_set_bypass_valve(handle, pct) bind(C, name="bwr_set_bypass_valve")
        type(c_ptr),    value :: handle
        real(c_double), value :: pct
        type(bwr_state_t), pointer :: sim
        type(main_steam_command_t) :: cmd
        call c_f_pointer(handle, sim)
        cmd%bypass_valve_set_pct = max(0.0_wp, min(100.0_wp, real(pct, wp)))
        call main_steam_apply(sim%ms, cmd, sim%ms_config)
    end subroutine bwr_set_bypass_valve

    subroutine bwr_trip_turbine(handle) bind(C, name="bwr_trip_turbine")
        ! Latch a TSV "fast closure" turbine trip. The four turbine stop
        ! valves ramp shut over the configured tsv_fast_close_time_s
        ! (default 0.1 s; spec 2.5 p.12) and the turbine_tripped flag
        ! latches once they reach zero. Idempotent.
        type(c_ptr), value :: handle
        type(bwr_state_t), pointer :: sim
        type(main_steam_command_t) :: cmd
        call c_f_pointer(handle, sim)
        cmd%tsv_fast_close = .true.
        call main_steam_apply(sim%ms, cmd, sim%ms_config)
    end subroutine bwr_trip_turbine

    subroutine bwr_reset_turbine_trip(handle) bind(C, name="bwr_reset_turbine_trip")
        type(c_ptr), value :: handle
        type(bwr_state_t), pointer :: sim
        type(main_steam_command_t) :: cmd
        call c_f_pointer(handle, sim)
        cmd%tsv_trip_reset = .true.
        call main_steam_apply(sim%ms, cmd, sim%ms_config)
    end subroutine bwr_reset_turbine_trip

    subroutine bwr_close_msivs(handle) bind(C, name="bwr_close_msivs")
        ! Latch the NSSSS isolation signal. All four MSIVs stroke shut
        ! over the configured msiv_stroke_time_s (default 4 s; spec 2.5
        ! "[BWR/4 typical 3-5 seconds full closure]"). Idempotent.
        type(c_ptr), value :: handle
        type(bwr_state_t), pointer :: sim
        type(main_steam_command_t) :: cmd
        call c_f_pointer(handle, sim)
        cmd%msiv_close = .true.
        call main_steam_apply(sim%ms, cmd, sim%ms_config)
    end subroutine bwr_close_msivs

    subroutine bwr_open_msivs(handle) bind(C, name="bwr_open_msivs")
        type(c_ptr), value :: handle
        type(bwr_state_t), pointer :: sim
        type(main_steam_command_t) :: cmd
        call c_f_pointer(handle, sim)
        cmd%msiv_open = .true.
        call main_steam_apply(sim%ms, cmd, sim%ms_config)
    end subroutine bwr_open_msivs

    !> Manual generator-breaker close. Honoured only if turbine speed is
    !> within ±sync_window_rpm of rated (1800 RPM at 60 Hz). Outside the
    !> window the call silently fails; poll bwr_get_sync_ready first.
    subroutine bwr_sync_generator(handle) bind(C, name="bwr_sync_generator")
        type(c_ptr), value :: handle
        type(bwr_state_t), pointer :: sim
        type(main_steam_command_t) :: cmd
        call c_f_pointer(handle, sim)
        cmd%sync_close = .true.
        call main_steam_apply(sim%ms, cmd, sim%ms_config)
    end subroutine bwr_sync_generator

    !> Open the generator breaker (load reject / desync). Drops the
    !> turbine back to free-spinning; without operator action on the TCV
    !> the rotor will overspeed.
    subroutine bwr_unsync_generator(handle) bind(C, name="bwr_unsync_generator")
        type(c_ptr), value :: handle
        type(bwr_state_t), pointer :: sim
        type(main_steam_command_t) :: cmd
        call c_f_pointer(handle, sim)
        cmd%sync_open = .true.
        call main_steam_apply(sim%ms, cmd, sim%ms_config)
    end subroutine bwr_unsync_generator

    function bwr_get_generator_synced(handle) result(v) &
            bind(C, name="bwr_get_generator_synced")
        type(c_ptr), value :: handle
        integer(c_int) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = merge(1, 0, sim%ms%generator_synced)
    end function bwr_get_generator_synced

    !> True (1) when RPM is inside the sync window and breaker can close.
    !> Frontend should latch a "SYNC READY" lamp on this.
    function bwr_get_sync_ready(handle) result(v) bind(C, name="bwr_get_sync_ready")
        type(c_ptr), value :: handle
        integer(c_int) :: v
        type(bwr_state_t), pointer :: sim
        real(wp) :: rpm_err
        call c_f_pointer(handle, sim)
        rpm_err = abs(sim%ms%turbine_speed_rpm - sim%ms_config%rated_speed_rpm)
        v = merge(1, 0, (.not. sim%ms%turbine_tripped) &
                   .and. (.not. sim%ms%generator_synced) &
                   .and. (rpm_err <= sim%ms_config%sync_window_rpm))
    end function bwr_get_sync_ready

    function bwr_get_turbine_mech_W(handle) result(v) bind(C, name="bwr_get_turbine_mech_W")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%ms%turbine_mech_W, c_double)
    end function bwr_get_turbine_mech_W

    function bwr_get_tcv_flow_norm(handle) result(v) bind(C, name="bwr_get_tcv_flow_norm")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%ms%tcv_flow_norm, c_double)
    end function bwr_get_tcv_flow_norm

    function bwr_get_bpv_flow_norm(handle) result(v) bind(C, name="bwr_get_bpv_flow_norm")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%ms%bpv_flow_norm, c_double)
    end function bwr_get_bpv_flow_norm

    function bwr_get_grid_frequency_Hz(handle) result(v) bind(C, name="bwr_get_grid_frequency_Hz")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        ! Sync ratio × 60 Hz. Free-spinning is open-loop speed / rated.
        v = real(sim%ms%turbine_speed_rpm / max(1.0_wp, sim%ms_config%rated_speed_rpm) &
              * 60.0_wp, c_double)
    end function bwr_get_grid_frequency_Hz

    !> Toggle the 3-element feedwater controller. 0 = manual demand,
    !> 1 = auto-track (ṁ_FW follows ṁ_steam + level trim).
    subroutine bwr_set_fw_controller_enabled(handle, on) &
            bind(C, name="bwr_set_fw_controller_enabled")
        type(c_ptr),    value :: handle
        integer(c_int), value :: on
        type(bwr_state_t), pointer :: sim
        type(feedwater_command_t) :: cmd
        call c_f_pointer(handle, sim)
        if (on == 0) then
            cmd%controller_enable_set = 0
        else
            cmd%controller_enable_set = 1
        end if
        call feedwater_apply(sim%fw, cmd)
    end subroutine bwr_set_fw_controller_enabled

    function bwr_get_fw_controller_enabled(handle) result(v) &
            bind(C, name="bwr_get_fw_controller_enabled")
        type(c_ptr), value :: handle
        integer(c_int) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = merge(1, 0, sim%fw%controller_enabled)
    end function bwr_get_fw_controller_enabled

    !> Select RHR mode. 0=STANDBY, 1=SHUTDOWN_COOLING, 2=POOL_COOLING,
    !> 3=CONTAINMENT_SPRAY, 4=LPCI. Permissive checks (SDC requires dome
    !> P < 135 psig, etc.) may silently reject the change; poll the
    !> getter to confirm.
    subroutine bwr_set_rhr_mode(handle, mode) bind(C, name="bwr_set_rhr_mode")
        type(c_ptr),    value :: handle
        integer(c_int), value :: mode
        type(bwr_state_t), pointer :: sim
        type(rhr_command_t) :: cmd
        type(rhr_drivers_t) :: drv
        call c_f_pointer(handle, sim)
        cmd%mode_set = int(mode)
        drv%dome_pressure_Pa     = sim%vessel%pressure_operating
        drv%reactor_coolant_T_K  = sim%avg_coolant_temp
        drv%suppression_pool_T_K = sim%supp_pool_T_K
        drv%drywell_pressure_Pa  = sim%drywell_pressure_Pa
        drv%reactor_level_m      = sim%vessel%Lrx
        call rhr_apply(sim%rhr, cmd, sim%rhr_config, drv)
    end subroutine bwr_set_rhr_mode

    function bwr_get_rhr_mode(handle) result(v) bind(C, name="bwr_get_rhr_mode")
        type(c_ptr), value :: handle
        integer(c_int) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = int(sim%rhr%mode, c_int)
    end function bwr_get_rhr_mode

    !> Start (1) or stop (0) a specific RHR loop pump (idx 1..4).
    subroutine bwr_set_rhr_pump(handle, idx, on) bind(C, name="bwr_set_rhr_pump")
        type(c_ptr),    value :: handle
        integer(c_int), value :: idx, on
        type(bwr_state_t), pointer :: sim
        type(rhr_command_t) :: cmd
        type(rhr_drivers_t) :: drv
        integer :: i
        call c_f_pointer(handle, sim)
        if (idx < 1 .or. idx > RHR_N_LOOPS) return
        do i = 1, RHR_N_LOOPS
            cmd%pump_set(i) = -1
        end do
        if (on == 0) then
            cmd%pump_set(int(idx)) = 0
        else
            cmd%pump_set(int(idx)) = 1
        end if
        drv%dome_pressure_Pa     = sim%vessel%pressure_operating
        drv%reactor_coolant_T_K  = sim%avg_coolant_temp
        drv%suppression_pool_T_K = sim%supp_pool_T_K
        drv%drywell_pressure_Pa  = sim%drywell_pressure_Pa
        drv%reactor_level_m      = sim%vessel%Lrx
        call rhr_apply(sim%rhr, cmd, sim%rhr_config, drv)
    end subroutine bwr_set_rhr_pump

    function bwr_get_rhr_pumps_running(handle) result(v) &
            bind(C, name="bwr_get_rhr_pumps_running")
        type(c_ptr), value :: handle
        integer(c_int) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = int(count(sim%rhr%pump_running .and. .not. sim%rhr%pump_tripped), c_int)
    end function bwr_get_rhr_pumps_running

    function bwr_get_rhr_total_heat_W(handle) result(v) &
            bind(C, name="bwr_get_rhr_total_heat_W")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%rhr%total_heat_W, c_double)
    end function bwr_get_rhr_total_heat_W

    function bwr_get_rhr_hx_outlet_T_K(handle) result(v) &
            bind(C, name="bwr_get_rhr_hx_outlet_T_K")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%rhr%hx_outlet_T_K, c_double)
    end function bwr_get_rhr_hx_outlet_T_K

    function bwr_get_supp_pool_T_K(handle) result(v) bind(C, name="bwr_get_supp_pool_T_K")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%supp_pool_T_K, c_double)
    end function bwr_get_supp_pool_T_K

    function bwr_get_hotwell_level_m(handle) result(v) bind(C, name="bwr_get_hotwell_level_m")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%cond%hotwell_level_m, c_double)
    end function bwr_get_hotwell_level_m

    function bwr_get_cst_level_m(handle) result(v) bind(C, name="bwr_get_cst_level_m")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%cond%cst_level_m, c_double)
    end function bwr_get_cst_level_m

    function bwr_get_cond_makeup_kg_s(handle) result(v) &
            bind(C, name="bwr_get_cond_makeup_kg_s")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%cond%makeup_flow_kg_s, c_double)
    end function bwr_get_cond_makeup_kg_s

    function bwr_get_cond_o2_ppb(handle) result(v) bind(C, name="bwr_get_cond_o2_ppb")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%cond%o2_ppb, c_double)
    end function bwr_get_cond_o2_ppb

    ! Add to the sustained operator reactivity perturbation [pcm]. The
    ! perturbation is reapplied every step (see apply_reactivity_perturbation),
    ! so the effect persists until cleared with bwr_reset_reactivity.
    subroutine bwr_apply_reactivity(handle, rho_pcm) bind(C, name="bwr_apply_reactivity")
        type(c_ptr),    value :: handle
        real(c_double), value :: rho_pcm
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        sim%reactivity_perturbation_pcm = sim%reactivity_perturbation_pcm &
                                        + real(rho_pcm, wp)
    end subroutine bwr_apply_reactivity

    subroutine bwr_reset_reactivity(handle) bind(C, name="bwr_reset_reactivity")
        type(c_ptr), value :: handle
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        sim%reactivity_perturbation_pcm = 0.0_wp
    end subroutine bwr_reset_reactivity

    function bwr_get_applied_reactivity_pcm(handle) result(v) &
            bind(C, name="bwr_get_applied_reactivity_pcm")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%reactivity_perturbation_pcm, c_double)
    end function bwr_get_applied_reactivity_pcm

    ! ── Scalar instrument readings ────────────────────────────────────────────

    function bwr_get_keff(handle) result(v) bind(C, name="bwr_get_keff")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%k_eff, c_double)
    end function bwr_get_keff

    function bwr_get_reactivity_pcm(handle) result(v) bind(C, name="bwr_get_reactivity_pcm")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%reactivity_pcm, c_double)
    end function bwr_get_reactivity_pcm

    function bwr_get_total_power_W(handle) result(v) bind(C, name="bwr_get_total_power_W")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%power_current, c_double)
    end function bwr_get_total_power_W

    function bwr_get_max_fuel_temp_K(handle) result(v) bind(C, name="bwr_get_max_fuel_temp_K")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%max_fuel_temp, c_double)
    end function bwr_get_max_fuel_temp_K

    function bwr_get_avg_void_fraction(handle) result(v) bind(C, name="bwr_get_avg_void_fraction")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%avg_void_fraction, c_double)
    end function bwr_get_avg_void_fraction

    function bwr_get_min_chfr(handle) result(v) bind(C, name="bwr_get_min_chfr")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%min_chfr, c_double)
    end function bwr_get_min_chfr

    function bwr_get_pressure_Pa(handle) result(v) bind(C, name="bwr_get_pressure_Pa")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%vessel%pressure_operating, c_double)
    end function bwr_get_pressure_Pa

    function bwr_get_reactor_level_m(handle) result(v) bind(C, name="bwr_get_reactor_level_m")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%vessel%Lrx, c_double)
    end function bwr_get_reactor_level_m

    function bwr_get_sat_temperature_K(handle) result(v) bind(C, name="bwr_get_sat_temperature_K")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%vessel%sat_temperature, c_double)
    end function bwr_get_sat_temperature_K

    function bwr_get_avg_coolant_temp_K(handle) result(v) &
            bind(C, name="bwr_get_avg_coolant_temp_K")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%avg_coolant_temp, c_double)
    end function bwr_get_avg_coolant_temp_K

    !> Dump a verbose two-line state summary to stdout. Useful when the
    !> frontend (Godot, etc.) can pipe Fortran stdout. The first line covers
    !> reactor thermal/neutronic state; the second covers level balance
    !> (steam evaporation vs feedwater inflow) and valve positions.
    subroutine bwr_print_state_summary(handle) bind(C, name="bwr_print_state_summary")
        type(c_ptr), value :: handle
        type(bwr_state_t), pointer :: sim
        real(wp) :: steam_kg_s, fw_kg_s
        call c_f_pointer(handle, sim)
        steam_kg_s = sim%ms%steam_flow_norm * sim%ms_config%rated_steam_flow_kg_s
        fw_kg_s    = sim%fw%flow_kg_s
        print '(A,F8.2,A,F8.5,A,ES10.3,A,F8.2,A,F8.2,A,F8.2,A,ES10.2,A,F6.2,A,F7.2,A,F6.3)', &
            'BWR t=', sim%time, &
            ' k=',     sim%k_eff, &
            ' P_MW=',  sim%power_current / 1.0e6_wp, &
            ' Tfuel=', sim%max_fuel_temp, &
            ' Tcool=', sim%avg_coolant_temp, &
            ' Tsat=',  sim%vessel%sat_temperature, &
            ' CHFR=',  sim%min_chfr, &
            ' void%=', sim%avg_void_fraction, &
            ' p_bar=', sim%vessel%pressure_operating / 1.0e5_wp, &
            ' Lrx=',   sim%vessel%Lrx
        print '(A,F7.1,A,F7.1,A,SP,F7.1,SS,A,F5.1,A,F5.1,A,F5.1,A,F5.1,A,F6.2)', &
            '    steam_kg_s=', steam_kg_s, &
            ' fw_kg_s=',       fw_kg_s, &
            ' net_kg_s=',      fw_kg_s - steam_kg_s, &
            ' TV%=', sim%ms%turbine_valve_pct, &
            ' BV%=', sim%ms%bypass_valve_pct, &
            ' FW%=', sim%fw%flow_pct, &
            ' rod%=', sim%crd%blade_insertion(1) * 100.0_wp, &
            ' pump%=', sum(sim%recirc%pump_speed_pct) / real(size(sim%recirc%pump_speed_pct), wp)
        print '(A,F7.1,A,F6.2,A,F7.2,A,F7.2,A,L2,A,I2,A,F6.1,A,F5.2,A,F6.2)', &
            '    rpm=',  sim%ms%turbine_speed_rpm, &
            ' Hz=',      sim%ms%turbine_speed_rpm / max(1.0_wp, sim%ms_config%rated_speed_rpm) * 60.0_wp, &
            ' Pmech_MW=', sim%ms%turbine_mech_W / 1.0e6_wp, &
            ' Pgen_MW=',  sim%ms%turbine_power_W / 1.0e6_wp, &
            ' synced=', sim%ms%generator_synced, &
            ' rhr_mode=', sim%rhr%mode, &
            ' poolT_K=', sim%supp_pool_T_K, &
            ' hotwell=', sim%cond%hotwell_level_m, &
            ' cst=',     sim%cond%cst_level_m
        flush(6)
    end subroutine bwr_print_state_summary

    function bwr_get_aprm_pct(handle) result(v) bind(C, name="bwr_get_aprm_pct")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%aprm, c_double)
    end function bwr_get_aprm_pct

    function bwr_get_reactor_period_s(handle) result(v) bind(C, name="bwr_get_reactor_period_s")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%reactor_period, c_double)
    end function bwr_get_reactor_period_s

    function bwr_get_turbine_power_W(handle) result(v) bind(C, name="bwr_get_turbine_power_W")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%ms%turbine_power_W, c_double)
    end function bwr_get_turbine_power_W

    function bwr_get_turbine_speed_rpm(handle) result(v) bind(C, name="bwr_get_turbine_speed_rpm")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%ms%turbine_speed_rpm, c_double)
    end function bwr_get_turbine_speed_rpm

    function bwr_get_turbine_valve_pct(handle) result(v) &
            bind(C, name="bwr_get_turbine_valve_pct")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%ms%turbine_valve_pct, c_double)
    end function bwr_get_turbine_valve_pct

    function bwr_get_bypass_valve_pct(handle) result(v) &
            bind(C, name="bwr_get_bypass_valve_pct")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%ms%bypass_valve_pct, c_double)
    end function bwr_get_bypass_valve_pct

    function bwr_get_msiv_open_frac(handle) result(v) &
            bind(C, name="bwr_get_msiv_open_frac")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sum(sim%ms%msiv_pos_pct) / (real(MS_N_MSLS, wp) * 100.0_wp), c_double)
    end function bwr_get_msiv_open_frac

    function bwr_get_srv_count_open(handle) result(v) &
            bind(C, name="bwr_get_srv_count_open")
        type(c_ptr), value :: handle
        integer(c_int) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = int(count(sim%ms%srv_open), c_int)
    end function bwr_get_srv_count_open

    function bwr_get_srv_flow_kg_s(handle) result(v) &
            bind(C, name="bwr_get_srv_flow_kg_s")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%ms%srv_flow_kg_s, c_double)
    end function bwr_get_srv_flow_kg_s

    !> Steam flow leaving the dome, normalised to rated [0..1+]. This is
    !> the *evaporation rate* the operator sees — turbine_valve × pressure
    !> ratio × MSIV open fraction. Pair with feedwater flow to monitor the
    !> inflow/outflow balance that drives reactor level.
    function bwr_get_steam_flow_norm(handle) result(v) &
            bind(C, name="bwr_get_steam_flow_norm")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%ms%steam_flow_norm, c_double)
    end function bwr_get_steam_flow_norm

    function bwr_get_steam_flow_kg_s(handle) result(v) &
            bind(C, name="bwr_get_steam_flow_kg_s")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%ms%steam_flow_norm * sim%ms_config%rated_steam_flow_kg_s, c_double)
    end function bwr_get_steam_flow_kg_s

    !> Feedwater mass-flow [kg/s] — the inflow that balances boil-off.
    !> Pair with bwr_get_steam_flow_kg_s to track level balance.
    function bwr_get_feedwater_flow_kg_s(handle) result(v) &
            bind(C, name="bwr_get_feedwater_flow_kg_s")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%fw%flow_kg_s, c_double)
    end function bwr_get_feedwater_flow_kg_s

    !> Simulated time [s] since reactor_create.
    function bwr_get_time_s(handle) result(v) bind(C, name="bwr_get_time_s")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%time, c_double)
    end function bwr_get_time_s

    !> Number of do_step calls since reactor_create.
    function bwr_get_step_count(handle) result(v) bind(C, name="bwr_get_step_count")
        type(c_ptr), value :: handle
        integer(c_int) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = int(sim%n_steps, c_int)
    end function bwr_get_step_count

    !> Rated thermal power [W] passed to reactor_create. Frontends should
    !> use this as the denominator for relative-power displays so the bar
    !> graph scales independently of the operating point.
    function bwr_get_power_rated_W(handle) result(v) bind(C, name="bwr_get_power_rated_W")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%power_rated, c_double)
    end function bwr_get_power_rated_W

    !> Turbine trip latch (1 = tripped, 0 = running).
    function bwr_get_turbine_tripped(handle) result(v) bind(C, name="bwr_get_turbine_tripped")
        type(c_ptr), value :: handle
        integer(c_int) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = merge(1, 0, sim%ms%turbine_tripped)
    end function bwr_get_turbine_tripped

    !> SCRAM latch (1 = scrammed, rods being driven in / fully in; 0 = normal).
    function bwr_get_scram_latched(handle) result(v) bind(C, name="bwr_get_scram_latched")
        type(c_ptr), value :: handle
        integer(c_int) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = merge(1, 0, sim%crd%scram_latched)
    end function bwr_get_scram_latched

    !> Rod Block Monitor latch (1 = withdrawal held due to short period,
    !> 0 = motion permitted). Insertion is never blocked.
    function bwr_get_rod_block_active(handle) result(v) &
            bind(C, name="bwr_get_rod_block_active")
        type(c_ptr), value :: handle
        integer(c_int) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = merge(1, 0, sim%crd%rod_block_active)
    end function bwr_get_rod_block_active

    !> Rod demand fraction — the operator's requested position, slewed by
    !> crd_step toward the actual blade insertion. Useful as a leading
    !> indicator: when demand differs from position, rods are still moving
    !> (or held by the RBM).
    function bwr_get_rod_demand(handle) result(v) bind(C, name="bwr_get_rod_demand")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%crd%bank_demand, c_double)
    end function bwr_get_rod_demand

    !> MSIV close latch (1 = isolation commanded / in progress, 0 = open path).
    function bwr_get_msiv_close_latched(handle) result(v) &
            bind(C, name="bwr_get_msiv_close_latched")
        type(c_ptr), value :: handle
        integer(c_int) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = merge(1, 0, sim%ms%msiv_close_latched)
    end function bwr_get_msiv_close_latched

    function bwr_get_mass_flux(handle) result(v) bind(C, name="bwr_get_mass_flux")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%recirc%mass_flux_kg_m2_s, c_double)
    end function bwr_get_mass_flux

    function bwr_get_core_flow_kg_s(handle) result(v) bind(C, name="bwr_get_core_flow_kg_s")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%recirc%core_flow_kg_s, c_double)
    end function bwr_get_core_flow_kg_s

    function bwr_get_core_flow_pct(handle) result(v) bind(C, name="bwr_get_core_flow_pct")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%recirc%core_flow_kg_s &
               / max(1.0_wp, sim%recirc_config%rated_core_flow_kg_s) * 100.0_wp, c_double)
    end function bwr_get_core_flow_pct

    function bwr_get_recirc_pump_speed_pct(handle, idx) result(v) &
            bind(C, name="bwr_get_recirc_pump_speed_pct")
        ! idx is 1-based (1 or 2). Out-of-range returns NaN-equivalent.
        type(c_ptr),    value :: handle
        integer(c_int), value :: idx
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        if (idx >= 1 .and. idx <= RECIRC_N_PUMPS) then
            v = real(sim%recirc%pump_speed_pct(int(idx)), c_double)
        else
            v = -1.0_c_double
        end if
    end function bwr_get_recirc_pump_speed_pct

    function bwr_get_control_rod_position(handle) result(v) &
            bind(C, name="bwr_get_control_rod_position")
        ! Returns the average insertion fraction across all 137 blades.
        ! For uniform-bank operation (the common case) this is the same
        ! as any individual blade's position.
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sum(sim%crd%blade_insertion) / real(CRD_N_BLADES, wp), c_double)
    end function bwr_get_control_rod_position

    function bwr_get_avg_burnup(handle) result(v) bind(C, name="bwr_get_avg_burnup")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%avg_burnup, c_double)
    end function bwr_get_avg_burnup

    ! ── Safety trips ──────────────────────────────────────────────────────────
    ! Returns a bitmask of currently-tripped thresholds; see BWR_TRIP_* in
    ! nuclear_physics_bwr.h. Non-zero means at least one SCRAM condition is
    ! met. Safe to poll between steps for pre-trip warning UI.
    function bwr_get_trip_flags(handle) result(flags) bind(C, name="bwr_get_trip_flags")
        type(c_ptr), value :: handle
        integer(c_int) :: flags
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        flags = 0
        if (sim%max_fuel_temp > TRIP_FUEL_TEMP_MAX_K) &
            flags = ior(flags, BWR_TRIP_FUEL_TEMP_HIGH)
        if (sim%power_current > sim%power_rated * TRIP_POWER_OVERSHOOT) &
            flags = ior(flags, BWR_TRIP_POWER_HIGH)
        if (sim%reactor_period > 0.0_wp .and. sim%reactor_period < TRIP_PERIOD_MIN_S) &
            flags = ior(flags, BWR_TRIP_SHORT_PERIOD)
        if (sim%vessel%Lrx < TRIP_LEVEL_MIN_M) &
            flags = ior(flags, BWR_TRIP_LEVEL_LOW)
        if (sim%vessel%pressure_operating > TRIP_PRESSURE_MAX_PA) &
            flags = ior(flags, BWR_TRIP_PRESSURE_HIGH)
        ! Low MSL pressure trip is the RUN-mode MSIV isolation: it detects a
        ! steam-line depressurization while the reactor is making steam for
        ! the turbine. In SHUTDOWN/STARTUP modes the mode switch bypasses
        ! it. We approximate that bypass by requiring meaningful power
        ! output (>1 % rated) — at cold standby the dome sits at 1 atm
        ! and there is no steam line to protect.
        if (sim%vessel%pressure_operating < TRIP_PRESSURE_MIN_PA .and. &
            sim%power_current > 0.01_wp * sim%power_rated) &
            flags = ior(flags, BWR_TRIP_PRESSURE_LOW)
        if (sim%min_chfr < TRIP_CHFR_MIN) then
            ! Log spurious zero-power CHFR events to aid root-cause diagnosis.
            if (sim%power_current < 1.0_wp) then
                block
                    integer :: ios
                    open(unit=87, file="/tmp/reaktor_chfr_debug.txt", &
                         position="append", action="write", iostat=ios)
                    if (ios == 0) then
                        write(87, '(2(A,ES12.4),A,F8.3,A,I0)') &
                            "CHFR<1 at zero power: min_chfr=", sim%min_chfr, &
                            "  power=", sim%power_current, &
                            "  Lrx=", sim%vessel%Lrx, &
                            "  steps=", sim%n_steps
                        close(87)
                    end if
                end block
            else
                flags = ior(flags, BWR_TRIP_CHFR_LOW)
            end if
        end if
    end function bwr_get_trip_flags

    ! ── 3-D field getters ─────────────────────────────────────────────────────
    ! Output arrays are laid out in Fortran column-major order (i varies fastest).
    ! Element at 1-based cell (i,j,k) is out[i-1 + nx*(j-1 + ny*(k-1))].

    subroutine bwr_get_power_density(handle, out) bind(C, name="bwr_get_power_density")
        type(c_ptr), value :: handle
        real(c_double), intent(out) :: out(*)
        type(bwr_state_t), pointer :: sim
        integer :: i, j, k, idx
        call c_f_pointer(handle, sim)
        idx = 1
        do k = 1, sim%nz
            do j = 1, sim%ny
                do i = 1, sim%nx
                    ! neutronics%power_density is W/cm³; convert to W/m³
                    out(idx) = real(sim%neutronics%power_density(i,j,k) * 1.0e6_wp, c_double)
                    idx = idx + 1
                end do
            end do
        end do
    end subroutine bwr_get_power_density

    subroutine bwr_get_temperature(handle, out) bind(C, name="bwr_get_temperature")
        type(c_ptr), value :: handle
        real(c_double), intent(out) :: out(*)
        type(bwr_state_t), pointer :: sim
        integer :: i, j, k, idx
        call c_f_pointer(handle, sim)
        idx = 1
        do k = 1, sim%nz
            do j = 1, sim%ny
                do i = 1, sim%nx
                    out(idx) = real(sim%heat%T(i,j,k), c_double)
                    idx = idx + 1
                end do
            end do
        end do
    end subroutine bwr_get_temperature

    subroutine bwr_get_void_fraction(handle, out) bind(C, name="bwr_get_void_fraction")
        type(c_ptr), value :: handle
        real(c_double), intent(out) :: out(*)
        type(bwr_state_t), pointer :: sim
        integer :: i, j, k, idx
        call c_f_pointer(handle, sim)
        idx = 1
        do k = 1, sim%nz
            do j = 1, sim%ny
                do i = 1, sim%nx
                    out(idx) = real(sim%thermalhydraulics%void_fraction(i,j,k), c_double)
                    idx = idx + 1
                end do
            end do
        end do
    end subroutine bwr_get_void_fraction

    subroutine bwr_get_burnup_field(handle, out) bind(C, name="bwr_get_burnup_field")
        type(c_ptr), value :: handle
        real(c_double), intent(out) :: out(*)
        type(bwr_state_t), pointer :: sim
        integer :: i, j, k, idx
        call c_f_pointer(handle, sim)
        idx = 1
        do k = 1, sim%nz
            do j = 1, sim%ny
                do i = 1, sim%nx
                    out(idx) = real(sim%burnup%burnup(i,j,k), c_double)
                    idx = idx + 1
                end do
            end do
        end do
    end subroutine bwr_get_burnup_field

    subroutine bwr_get_xenon(handle, out) bind(C, name="bwr_get_xenon")
        type(c_ptr), value :: handle
        real(c_double), intent(out) :: out(*)
        type(bwr_state_t), pointer :: sim
        integer :: i, j, k, idx
        call c_f_pointer(handle, sim)
        idx = 1
        do k = 1, sim%nz
            do j = 1, sim%ny
                do i = 1, sim%nx
                    out(idx) = real(sim%burnup%Xe135(i,j,k), c_double)
                    idx = idx + 1
                end do
            end do
        end do
    end subroutine bwr_get_xenon

    ! ── Grid dimensions ───────────────────────────────────────────────────────

    function bwr_get_nx(handle) result(v) bind(C, name="bwr_get_nx")
        type(c_ptr), value :: handle
        integer(c_int) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = int(sim%nx, c_int)
    end function bwr_get_nx

    function bwr_get_ny(handle) result(v) bind(C, name="bwr_get_ny")
        type(c_ptr), value :: handle
        integer(c_int) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = int(sim%ny, c_int)
    end function bwr_get_ny

    function bwr_get_nz(handle) result(v) bind(C, name="bwr_get_nz")
        type(c_ptr), value :: handle
        integer(c_int) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = int(sim%nz, c_int)
    end function bwr_get_nz

end module bwr_c_interface