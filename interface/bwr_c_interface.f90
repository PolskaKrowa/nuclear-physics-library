! interface/bwr_c_interface.f90
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
        real(wp) :: pressure_operating
        real(wp) :: mass_flux_core
        real(wp) :: inlet_temperature

        real(wp) :: k_eff, reactivity_pcm
        real(wp) :: max_fuel_temp
        real(wp) :: avg_void_fraction, max_void_fraction
        real(wp) :: min_chfr, avg_burnup

        real(wp) :: alpha_doppler, alpha_void

        real(wp) :: time
        integer  :: n_steps

        real(wp) :: sat_temperature
        real(wp) :: Lbase, Lrx
        real(wp) :: avg_coolant_temp
        real(wp) :: coolant_T_prev1, coolant_T_prev2
        real(wp) :: P_prev1, P_prev2
        real(wp) :: T_excess

        real(wp) :: neutrons_prev, reactor_period
        real(wp) :: aprm
        real(wp) :: rf_void

        real(wp) :: rod_bank_position

        real(wp) :: turbine_valve, bypass_valve
        real(wp) :: steam_flow, turbine_speed, turbine_power
        logical  :: turbine_tripped
        real(wp) :: feedwater_temp
    end type bwr_state_t

contains

    ! ══════════════════════════════════════════════════════════════════════════
    !  PURE PHYSICS HELPERS  (identical to bwr_core_simulator_opengl.f90)
    ! ══════════════════════════════════════════════════════════════════════════

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

    pure function void_reactivity_factor(avg_void_frac) result(rf)
        real(wp), intent(in) :: avg_void_frac
        real(wp) :: rf, dens2
        dens2 = (1.0_wp - min(1.0_wp, max(0.0_wp, avg_void_frac))) * 1000.0_wp
        rf = -5.84e-10_wp    * dens2**3  &
             + 1.35422e-7_wp  * dens2**2  &
             + 1.358252042e-3_wp * dens2  &
             + 0.090469568418_wp
    end function void_reactivity_factor

    ! ══════════════════════════════════════════════════════════════════════════
    !  INTERNAL PHYSICS SUBROUTINES
    ! ══════════════════════════════════════════════════════════════════════════

    subroutine update_cross_sections_feedback(sim, T, rho)
        type(bwr_state_t), intent(inout) :: sim
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

    subroutine update_instrumentation(sim)
        type(bwr_state_t), intent(inout) :: sim
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

    subroutine update_pressure_dynamics(sim, dt)
        type(bwr_state_t), intent(inout) :: sim
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

    subroutine update_turbine(sim, dt)
        type(bwr_state_t), intent(inout) :: sim
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
            eta = ETA_CARNOT * ePres * (sim%turbine_speed / RATED_RPM) &
                                     * (sim%turbine_valve / 100.0_wp)
            sim%turbine_power = eta * sim%power_current
        end if
    end subroutine update_turbine

    ! ── Geometry: pin-by-pin BWR lattice ──────────────────────────────────────

    subroutine setup_geometry(sim)
        type(bwr_state_t), intent(inout) :: sim

        integer  :: i, j, k
        real(wp) :: x, y, r, x_center, y_center
        integer  :: pin_i, pin_j
        real(wp) :: dx_from_pin, dy_from_pin, r_from_pin
        logical  :: in_fuel, found_water, found_fuel
        type(xsec_material_t) :: xsec_fuel, xsec_water

        real(wp), parameter :: pitch      = 0.0163_wp   ! 1.63 cm pin pitch
        real(wp), parameter :: pin_radius = 0.0052_wp   ! 5.2 mm pellet radius

        call xslib_get_material(sim%xslib, "UO2_35", xsec_fuel, found_fuel)
        if (.not. found_fuel) then
            print *, "bwr_c_interface: UO2_35 material not found"; stop 1
        end if
        call xslib_get_material(sim%xslib, "H2O", xsec_water, found_water)
        if (.not. found_water) then
            print *, "bwr_c_interface: H2O material not found"; stop 1
        end if

        do k = 1, sim%nz
            do j = 1, sim%ny
                do i = 1, sim%nx
                    x = (real(i, wp) - 0.5_wp) * sim%dx
                    y = (real(j, wp) - 0.5_wp) * sim%dy

                    pin_i = nint(x / pitch)
                    pin_j = nint(y / pitch)
                    x_center = real(pin_i, wp) * pitch
                    y_center = real(pin_j, wp) * pitch
                    dx_from_pin = x - x_center
                    dy_from_pin = y - y_center
                    r_from_pin  = sqrt(dx_from_pin**2 + dy_from_pin**2)

                    in_fuel = (r_from_pin < pin_radius)

                    r = sqrt((x - sim%core_diameter/2.0_wp)**2 + &
                             (y - sim%core_diameter/2.0_wp)**2)
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
    end subroutine setup_geometry

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
             sim%dx*100.0_wp, sim%dy*100.0_wp, sim%dz*100.0_wp, neutron_config)
        call setup_geometry(sim)

        ! ── Heat transfer ──────────────────────────────────────────────────────
        heat_config%include_convection = .true.
        call heat_init(sim%heat, sim%nx, sim%ny, sim%nz, &
                       sim%dx, sim%dy, sim%dz, heat_config)
        call heat_set_properties(sim%heat, k=3.0_wp, rho=10970.0_wp, cp=300.0_wp, &
            i1=1, i2=int(0.9_wp*sim%nx), j1=1, j2=int(0.9_wp*sim%ny))
        call heat_set_properties(sim%heat, k=0.6_wp, rho=738.0_wp, cp=5200.0_wp, &
            i1=int(0.9_wp*sim%nx)+1, i2=sim%nx, &
            j1=int(0.9_wp*sim%ny)+1, j2=sim%ny)

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
        logical  :: converged

        do iter = 1, 50
            call mg_solve_eigenvalue(sim%neutronics, sim%k_eff, converged)
            if (.not. converged) &
                print '(A,I3)', "  Warning: neutronics did not converge, iter ", iter

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
        real(wp), parameter :: rho_liquid = 738.0_wp
        real(wp), parameter :: rho_vapor  = 0.038_wp
        logical  :: converged
        real(wp) :: v_coolant
        real(wp), parameter :: D_h     = 0.01_wp
        real(wp), parameter :: mu      = 0.0001_wp
        real(wp), parameter :: k_fluid = 0.6_wp
        real(wp), parameter :: Pr      = 0.9_wp

        v_coolant = sim%mass_flux_core / rho_liquid
        if (allocated(sim%heat%vz)) sim%heat%vz = v_coolant

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

        call rng_seed(123456789_8)
        allocate(sim)

        sim%nx = nx;  sim%ny = ny;  sim%nz = nz
        sim%dx = dx;  sim%dy = dy;  sim%dz = dz
        sim%core_height   = core_height_m
        sim%core_diameter = core_diameter_m

        sim%power_rated        = power_rated_W
        sim%pressure_operating = pressure_Pa
        sim%mass_flux_core     = mass_flux
        sim%inlet_temperature  = inlet_temp_K

        sim%alpha_doppler     = -3.5_wp
        sim%alpha_void        = -80.0_wp
        sim%rod_bank_position = 0.95_wp

        sim%sat_temperature  = sat_temp_K(pressure_Pa)
        sim%Lbase            = core_height_m + 1.5_wp
        sim%Lrx              = 1.5_wp
        sim%avg_coolant_temp = inlet_temp_K
        sim%coolant_T_prev1  = inlet_temp_K
        sim%coolant_T_prev2  = inlet_temp_K
        sim%P_prev1          = pressure_Pa
        sim%P_prev2          = pressure_Pa
        sim%T_excess         = 0.0_wp

        sim%neutrons_prev  = 0.0_wp
        sim%reactor_period = 1.0e6_wp
        sim%aprm           = 0.0_wp
        sim%rf_void        = void_reactivity_factor(0.0_wp)

        sim%turbine_valve   = 33.0_wp
        sim%bypass_valve    = 0.0_wp
        sim%steam_flow      = 0.0_wp
        sim%turbine_speed   = 3600.0_wp
        sim%turbine_power   = 0.0_wp
        sim%turbine_tripped = .false.
        sim%feedwater_temp  = inlet_temp_K

        sim%power_current     = 0.0_wp
        sim%k_eff             = 1.0_wp
        sim%reactivity_pcm    = 0.0_wp
        sim%max_fuel_temp     = inlet_temp_K
        sim%avg_void_fraction = 0.0_wp
        sim%max_void_fraction = 0.0_wp
        sim%min_chfr          = 999.0_wp
        sim%avg_burnup        = 0.0_wp
        sim%time              = 0.0_wp
        sim%n_steps           = 0

        call init_physics(sim, real(enrichment, wp))

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
        deallocate(sim)
    end subroutine bwr_reactor_destroy

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
        call c_f_pointer(handle, sim)
        call do_step(sim, real(dt, wp))
    end subroutine bwr_step

    ! ── Operator controls ─────────────────────────────────────────────────────

    subroutine bwr_set_control_rod_position(handle, insertion) &
            bind(C, name="bwr_set_control_rod_position")
        type(c_ptr),    value :: handle
        real(c_double), value :: insertion
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        sim%rod_bank_position = max(0.0_wp, min(1.0_wp, real(insertion, wp)))
    end subroutine bwr_set_control_rod_position

    subroutine bwr_set_mass_flux(handle, kg_m2_s) bind(C, name="bwr_set_mass_flux")
        type(c_ptr),    value :: handle
        real(c_double), value :: kg_m2_s
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        sim%mass_flux_core = max(100.0_wp, real(kg_m2_s, wp))
    end subroutine bwr_set_mass_flux

    subroutine bwr_set_turbine_valve(handle, pct) bind(C, name="bwr_set_turbine_valve")
        type(c_ptr),    value :: handle
        real(c_double), value :: pct
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        sim%turbine_valve = max(0.0_wp, min(100.0_wp, real(pct, wp)))
    end subroutine bwr_set_turbine_valve

    subroutine bwr_set_bypass_valve(handle, pct) bind(C, name="bwr_set_bypass_valve")
        type(c_ptr),    value :: handle
        real(c_double), value :: pct
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        sim%bypass_valve = max(0.0_wp, min(100.0_wp, real(pct, wp)))
    end subroutine bwr_set_bypass_valve

    subroutine bwr_apply_reactivity(handle, rho_pcm) bind(C, name="bwr_apply_reactivity")
        type(c_ptr),    value :: handle
        real(c_double), value :: rho_pcm
        type(bwr_state_t), pointer :: sim
        integer :: i, j, k, g
        real(wp) :: factor
        call c_f_pointer(handle, sim)
        factor = 1.0_wp + real(rho_pcm, wp) / 1.0e5_wp
        do k = 1, sim%nz
            do j = 1, sim%ny
                do i = 1, sim%nx
                    do g = 1, sim%neutronics%n_groups
                        sim%neutronics%xsec(i,j,k)%nu_sigma_f(g) = &
                            sim%neutronics%xsec(i,j,k)%nu_sigma_f(g) * factor
                    end do
                end do
            end do
        end do
    end subroutine bwr_apply_reactivity

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
        v = real(sim%pressure_operating, c_double)
    end function bwr_get_pressure_Pa

    function bwr_get_reactor_level_m(handle) result(v) bind(C, name="bwr_get_reactor_level_m")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%Lrx, c_double)
    end function bwr_get_reactor_level_m

    function bwr_get_sat_temperature_K(handle) result(v) bind(C, name="bwr_get_sat_temperature_K")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%sat_temperature, c_double)
    end function bwr_get_sat_temperature_K

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
        v = real(sim%turbine_power, c_double)
    end function bwr_get_turbine_power_W

    function bwr_get_turbine_speed_rpm(handle) result(v) bind(C, name="bwr_get_turbine_speed_rpm")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%turbine_speed, c_double)
    end function bwr_get_turbine_speed_rpm

    function bwr_get_control_rod_position(handle) result(v) &
            bind(C, name="bwr_get_control_rod_position")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%rod_bank_position, c_double)
    end function bwr_get_control_rod_position

    function bwr_get_avg_burnup(handle) result(v) bind(C, name="bwr_get_avg_burnup")
        type(c_ptr), value :: handle
        real(c_double) :: v
        type(bwr_state_t), pointer :: sim
        call c_f_pointer(handle, sim)
        v = real(sim%avg_burnup, c_double)
    end function bwr_get_avg_burnup

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
