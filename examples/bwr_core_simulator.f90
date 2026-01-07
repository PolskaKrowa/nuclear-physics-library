! examples/bwr_core_simulator.f90
!
! Comprehensive Boiling Water Reactor (BWR) Core Simulator
! Demonstrates full coupling of the nuclear physics library modules:
! - Multi-group neutron diffusion
! - Two-phase thermal-hydraulics
! - Heat transfer (fuel to coolant)
! - Cross-section feedback
! - Burnup and depletion
! - Xenon/Samarium poisoning
!
! Simulates a typical BWR/4 reactor:
! - 764 fuel assemblies (8x8 rod array per assembly)
! - Thermal power: 2,381 MWth
! - Operating pressure: 7.14 MPa
! - Core flow: ~13,500 kg/s
! - Active height: 3.71 m
!
program bwr_core_simulator
    use kinds, only: wp
    
    ! Physics modules
    use multigroup_diffusion
    use cross_sections
    use burnup_depletion
    use heat_transfer
    use two_phase_flow
    use fluid_dynamics
    use pressure_dynamics
    
    ! Utilities
    use rng, only: rng_seed
    
    implicit none
    
    ! Simulation state
type :: simulation_t
        ! Grid
        integer :: nx, ny, nz
        real(wp) :: dx, dy, dz
        real(wp) :: core_height, core_diameter
        
        ! Physics modules
        type(mg_state_t) :: neutronics
        type(heat_state_t) :: heat
        type(two_phase_state_t) :: thermalhydraulics
        type(fluid_state_t) :: coolant
        type(pressure_state_t) :: pressure
        type(burnup_state_t) :: burnup
        type(xsec_library_t) :: xslib
        
        ! Operating conditions
        real(wp) :: power_rated        ! Rated thermal power [W]
        real(wp) :: power_current      ! Current power [W]
        real(wp) :: pressure_operating ! Operating pressure [Pa]
        real(wp) :: mass_flux_core     ! Core mass flux [kg/m²·s]
        real(wp) :: inlet_temperature  ! Inlet coolant temperature [K]
        
        ! Core parameters
        real(wp) :: k_eff              ! Effective multiplication factor
        real(wp) :: reactivity_pcm     ! Reactivity [pcm]
        real(wp) :: max_fuel_temp      ! Maximum fuel temperature [K]
        real(wp) :: max_clad_temp      ! Maximum cladding temperature [K]
        real(wp) :: avg_void_fraction  ! Average void fraction
        real(wp) :: max_void_fraction  ! Maximum void fraction
        real(wp) :: min_chfr           ! Minimum CHFR
        real(wp) :: avg_burnup         ! Average burnup [MWd/kgU]
        
        ! Feedback coefficients
        real(wp) :: alpha_doppler      ! Doppler coefficient [pcm/K]
        real(wp) :: alpha_void         ! Void coefficient [pcm/%void]
        
        ! Time
        real(wp) :: time
        integer :: n_steps
    end type simulation_t
    type(simulation_t) :: sim
    
    ! Time stepping
    real(wp) :: t_end, dt, t
    integer :: n_steps, step
    
    ! Output
    character(len=256) :: output_dir

    logical :: converged
    converged = .false.
    
    print *, "=============================================="
    print *, "  BWR Core Simulator"
    print *, "  Nuclear Physics Library Demonstration"
    print *, "=============================================="
    print *, ""
    
    ! Initialize random number generator
    call rng_seed(123456789_8)
    
    ! Setup simulation
    call setup_bwr_simulation(sim)
    
    ! Initial steady-state solve
    print *, "Computing initial steady-state..."
    call solve_steady_state(sim)
    call print_steady_state_summary(sim)
    
    ! Time-dependent simulation
    print *, ""
    print *, "Starting transient simulation..."
    print *, ""
    
    t = 0.0_wp
    t_end = 100.0_wp  ! 100 seconds
    n_steps = 1000
    dt = t_end / real(n_steps, wp)
    
    output_dir = "bwr_output"
    call execute_command_line("mkdir -p " // trim(output_dir))
    
    ! Write initial state
    call write_output(sim, 0, output_dir)
    
    ! Time integration loop
    do step = 1, n_steps
        t = t + dt
        
        ! Apply transient (example: control rod withdrawal)
        if (t >= 10.0_wp .and. t < 15.0_wp) then
            call apply_reactivity_insertion(sim, 50.0_wp)  ! +50 pcm
        end if
        
        ! Coupled physics time step
        call coupled_time_step(sim, dt)
        
        ! Monitor and output
        if (mod(step, 50) == 0) then
            call print_transient_summary(sim, step, t)
            call write_output(sim, step, output_dir)
        end if
        
        ! Safety checks
        if (check_safety_limits(sim)) then
            print *, "WARNING: Safety limits approached at t = ", t
        end if
    end do
    
    print *, ""
    print *, "Simulation complete!"
    print *, "Output written to: ", trim(output_dir)
    
    ! Cleanup
    call cleanup_simulation(sim)
    
contains

    !> Complete BWR simulation state
        
    !> Setup BWR simulation
    subroutine setup_bwr_simulation(sim)
        type(simulation_t), intent(out) :: sim
        
        type(mg_config_t) :: neutron_config
        type(heat_config_t) :: heat_config
        type(two_phase_config_t) :: th_config
        type(burnup_config_t) :: burnup_config
        
        integer :: i, j, k
        type(xsec_material_t) :: xsec_fuel, xsec_water
        real(wp) :: r, r_fuel, r_water
        
        ! Core geometry (simplified: cylindrical core)
        sim%nx = 20  ! Radial zones
        sim%ny = 20
        sim%nz = 20  ! Axial zones
        
        sim%core_height = 3.71_wp      ! meters
        sim%core_diameter = 4.75_wp    ! meters (approximate)
        
        sim%dx = sim%core_diameter / real(sim%nx, wp)
        sim%dy = sim%core_diameter / real(sim%ny, wp)
        sim%dz = sim%core_height / real(sim%nz, wp)
        
        ! Operating conditions (typical BWR/4)
        sim%power_rated = 2381.0e6_wp          ! 2,381 MWth
        sim%pressure_operating = 7.14e6_wp     ! 7.14 MPa
        sim%mass_flux_core = 1500.0_wp         ! kg/m²·s
        sim%inlet_temperature = 551.0_wp       ! ~278°C
        
        ! Feedback coefficients (typical BWR values)
        sim%alpha_doppler = -2.5_wp     ! pcm/K (negative)
        sim%alpha_void = -50.0_wp       ! pcm/%void (strongly negative)
        
        ! Initialize physics modules
        print *, "Initializing physics modules..."
        
        ! 1. Cross-section library
        call xslib_init(sim%xslib, n_groups=2)
        
        ! Create fuel cross sections (3.5% enriched UO2)
        call xslib_create_two_group_fuel(xsec_fuel, enrichment=0.035_wp)
        call xslib_add_material(sim%xslib, xsec_fuel)
        
        ! Create water/steam cross sections
        call xslib_create_two_group_moderator(xsec_water)
        call xslib_add_material(sim%xslib, xsec_water)
        
        ! 2. Neutronics (2-group diffusion)
        neutron_config%n_groups = 2
        neutron_config%max_outer_iter = 100
        neutron_config%outer_tolerance = 1.0e-5_wp
        neutron_config%power_level = sim%power_rated
        neutron_config%normalize_power = .true.
        
        call mg_init(sim%neutronics, sim%nx, sim%ny, sim%nz, &
                    sim%dx, sim%dy, sim%dz, neutron_config)
        
        ! Set fuel region (cylindrical core)
        r_fuel = 0.45_wp * sim%core_diameter  ! 90% of diameter
        r_water = 0.5_wp * sim%core_diameter
        
        do k = 1, sim%nz
            do j = 1, sim%ny
                do i = 1, sim%nx
                    r = sqrt((i - sim%nx/2.0_wp)**2 + (j - sim%ny/2.0_wp)**2) * sim%dx
                    
                    if (r < r_fuel) then
                        call mg_set_cross_sections(sim%neutronics, xsec_fuel%xsec_base, i, i, j, j, k, k)
                    else
                        call mg_set_cross_sections(sim%neutronics, xsec_water%xsec_base, i, i, j, j, k, k)
                    end if
                end do
            end do
        end do
        
        ! 3. Heat transfer
        heat_config%include_convection = .true.
        call heat_init(sim%heat, sim%nx, sim%ny, sim%nz, &
                      sim%dx, sim%dy, sim%dz, heat_config)
        
        ! Fuel properties (UO2)
        call heat_set_properties(sim%heat, &
            k=3.0_wp, rho=10970.0_wp, cp=300.0_wp, &  ! Fuel
            i1=1, i2=int(0.9*sim%nx), j1=1, j2=int(0.9*sim%ny))
        
        ! Water properties
        call heat_set_properties(sim%heat, &
            k=0.6_wp, rho=738.0_wp, cp=5200.0_wp)  ! Water/steam
        
        ! 4. Two-phase flow
        th_config%void_correlation = VOID_CHEXAL_LELLOUCHE_ID
        th_config%chf_correlation = CHF_GROENEVELD_ID
        th_config%include_subcooled_boiling = .true.
        
        call two_phase_init(sim%thermalhydraulics, sim%nx, sim%ny, sim%nz, &
                           sim%dx, sim%dy, sim%dz, th_config)
        
        ! Set channel geometry (typical BWR fuel assembly)
        call two_phase_set_geometry(sim%thermalhydraulics, diameter=0.012_wp)
        
        ! 5. Burnup tracking
        burnup_config%track_xenon = .true.
        burnup_config%track_samarium = .true.
        burnup_config%track_actinides = .true.
        
        call burnup_init(sim%burnup, sim%nx, sim%ny, sim%nz, &
                        sim%dx, sim%dy, sim%dz, burnup_config)
        
        call burnup_set_initial_composition(sim%burnup, enrichment=0.035_wp)
        
        ! Initialize time
        sim%time = 0.0_wp
        sim%n_steps = 0
        
        print *, "  Grid: ", sim%nx, "x", sim%ny, "x", sim%nz
        print *, "  Core diameter: ", sim%core_diameter, " m"
        print *, "  Core height: ", sim%core_height, " m"
        print *, "  Rated power: ", sim%power_rated/1.0e6_wp, " MW"
        print *, ""
    end subroutine setup_bwr_simulation
    
    !> Solve for initial steady state
    subroutine solve_steady_state(sim)
        type(simulation_t), intent(inout) :: sim
        
        integer :: iter
        real(wp) :: power(sim%nx, sim%ny, sim%nz)
        real(wp) :: temperature(sim%nx, sim%ny, sim%nz)
        real(wp) :: density(sim%nx, sim%ny, sim%nz)
        real(wp) :: error
        
        ! Iterative coupling for steady state
        do iter = 1, 50
            ! 1. Solve neutronics with current T, ρ
            call mg_solve_eigenvalue(sim%neutronics, sim%k_eff, converged)
            
            if (.not. converged) then
                print *, "Warning: Neutronics did not converge"
            end if
            
            ! 2. Get power distribution
            call mg_get_power(sim%neutronics, power)
            sim%power_current = sim%neutronics%total_power
            
            ! 3. Solve heat transfer
            sim%heat%q = power
            call heat_step(sim%heat, 1.0_wp)  ! Large pseudo-timestep
            
            temperature = sim%heat%T
            sim%max_fuel_temp = maxval(temperature)
            
            ! 4. Solve two-phase flow
            density = 738.0_wp  ! Simplified
            call two_phase_step(sim%thermalhydraulics, &
                temperature, sim%pressure_operating * sim%heat%T / sim%heat%T, &
                sim%mass_flux_core * sim%heat%T / sim%heat%T, &
                power / (sim%dx * sim%dy), 1.0_wp)
            
            ! Get void fraction
            call two_phase_get_void_fraction(sim%thermalhydraulics, density)
            sim%avg_void_fraction = sum(density) / size(density) * 100.0_wp
            sim%max_void_fraction = maxval(density) * 100.0_wp
            
            ! 5. Update cross sections with feedback
            call update_cross_sections_feedback(sim, temperature, density)
            
            ! Check convergence
            error = abs(sim%k_eff - 1.0_wp) + &
                   abs(sim%power_current - sim%power_rated) / sim%power_rated
            
            if (mod(iter, 10) == 0) then
                print *, "  Iteration", iter, ": k_eff =", sim%k_eff, &
                        ", error =", error
            end if
            
            if (error < 1.0e-4_wp) exit
        end do
        
        ! Final calculations
        sim%reactivity_pcm = (sim%k_eff - 1.0_wp) / sim%k_eff * 1.0e5_wp
        sim%min_chfr = minval(sim%thermalhydraulics%chf_ratio)
        
    end subroutine solve_steady_state
    
    !> Coupled physics time step
    subroutine coupled_time_step(sim, dt)
        type(simulation_t), intent(inout) :: sim
        real(wp), intent(in) :: dt
        
        real(wp) :: power(sim%nx, sim%ny, sim%nz)
        real(wp) :: flux(sim%nx, sim%ny, sim%nz)
        real(wp) :: temperature(sim%nx, sim%ny, sim%nz)
        real(wp) :: void_fraction(sim%nx, sim%ny, sim%nz)
        
        ! 1. Point kinetics for fast neutronics response
        ! (Could use full space-time kinetics for accuracy)
        call mg_solve_eigenvalue(sim%neutronics, sim%k_eff, converged)
        
        ! 2. Power distribution
        call mg_get_power(sim%neutronics, power)
        call mg_get_flux(sim%neutronics, flux)
        sim%power_current = sim%neutronics%total_power
        
        ! 3. Update burnup
        call burnup_step(sim%burnup, flux, power, dt)
        sim%avg_burnup = sum(sim%burnup%burnup) / size(sim%burnup%burnup)
        
        ! 4. Heat transfer
        sim%heat%q = power
        call heat_step(sim%heat, dt)
        temperature = sim%heat%T
        sim%max_fuel_temp = maxval(temperature)
        
        ! 5. Two-phase thermal-hydraulics
        call two_phase_step(sim%thermalhydraulics, &
            temperature, &
            sim%pressure_operating + 0.0_wp * temperature, &
            sim%mass_flux_core + 0.0_wp * temperature, &
            power / (sim%dx * sim%dy), dt)
        
        call two_phase_get_void_fraction(sim%thermalhydraulics, void_fraction)
        sim%avg_void_fraction = sum(void_fraction) / size(void_fraction) * 100.0_wp
        sim%max_void_fraction = maxval(void_fraction) * 100.0_wp
        sim%min_chfr = minval(sim%thermalhydraulics%chf_ratio)
        
        ! 6. Update cross sections with feedback
        call update_cross_sections_feedback(sim, temperature, &
            738.0_wp * (1.0_wp - void_fraction))
        
        ! Update reactivity
        sim%reactivity_pcm = (sim%k_eff - 1.0_wp) / sim%k_eff * 1.0e5_wp
        
        sim%time = sim%time + dt
        sim%n_steps = sim%n_steps + 1
        
    end subroutine coupled_time_step
    
    !> Update cross sections with T and ρ feedback
    subroutine update_cross_sections_feedback(sim, T, rho)
        type(simulation_t), intent(inout) :: sim
        real(wp), intent(in) :: T(:, :, :)
        real(wp), intent(in) :: rho(:, :, :)
        
        integer :: i, j, k
        type(mg_xsec_t) :: xsec
        real(wp) :: T_fuel, rho_mod
        
        do k = 1, sim%nz
            do j = 1, sim%ny
                do i = 1, sim%nx
                    T_fuel = T(i, j, k)
                    rho_mod = rho(i, j, k)
                    
                    ! Get updated cross sections with feedback
                    call xslib_get_xsec(sim%xslib, "UO2_35", &
                        T_fuel, rho_mod, sim%burnup%burnup(i, j, k), xsec, &
                        Xe_conc=sim%burnup%Xe135(i, j, k), &
                        Sm_conc=sim%burnup%Sm149(i, j, k))
                    
                    ! Update neutronics
                    call mg_set_cross_sections(sim%neutronics, xsec, i, i, j, j, k, k)
                end do
            end do
        end do
    end subroutine update_cross_sections_feedback
    
    !> Apply reactivity insertion
    subroutine apply_reactivity_insertion(sim, rho_pcm)
        type(simulation_t), intent(inout) :: sim
        real(wp), intent(in) :: rho_pcm
        
        ! Simplified: multiply nu by (1 + rho)
        real(wp) :: factor
        integer :: i, j, k, g
        
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
    
    !> Check safety limits
    function check_safety_limits(sim) result(approaching_limits)
        type(simulation_t), intent(in) :: sim
        logical :: approaching_limits
        
        approaching_limits = .false.
        
        ! Fuel temperature limit: 1200°C (1473 K)
        if (sim%max_fuel_temp > 1373.0_wp) then
            approaching_limits = .true.
            print *, "  High fuel temperature: ", sim%max_fuel_temp, " K"
        end if
        
        ! CHFR limit: > 1.3
        if (sim%min_chfr < 1.5_wp) then
            approaching_limits = .true.
            print *, "  Low CHFR: ", sim%min_chfr
        end if
        
        ! Power limit: 120% rated
        if (sim%power_current > 1.2_wp * sim%power_rated) then
            approaching_limits = .true.
            print *, "  High power: ", sim%power_current / 1.0e6_wp, " MW"
        end if
    end function check_safety_limits
    
    !> Print steady-state summary
    subroutine print_steady_state_summary(sim)
        type(simulation_t), intent(in) :: sim
        
        print *, ""
        print *, "========== Steady-State Results =========="
        print *, "k_eff:                ", sim%k_eff
        print *, "Reactivity:           ", sim%reactivity_pcm, " pcm"
        print *, "Core power:           ", sim%power_current / 1.0e6_wp, " MW"
        print *, "Max fuel temp:        ", sim%max_fuel_temp - 273.15_wp, " °C"
        print *, "Avg void fraction:    ", sim%avg_void_fraction, " %"
        print *, "Max void fraction:    ", sim%max_void_fraction, " %"
        print *, "Min CHFR:             ", sim%min_chfr
        print *, "=========================================="
        print *, ""
    end subroutine print_steady_state_summary
    
    !> Print transient summary
    subroutine print_transient_summary(sim, step, t)
        type(simulation_t), intent(in) :: sim
        integer, intent(in) :: step
        real(wp), intent(in) :: t
        
        write(*, '(A,I6,A,F8.2,A)') "Step ", step, ", t = ", t, " s"
        write(*, '(A,F10.5,A,F8.1,A)') "  k_eff = ", sim%k_eff, &
            ", ρ = ", sim%reactivity_pcm, " pcm"
        write(*, '(A,F8.1,A)') "  Power = ", sim%power_current / 1.0e6_wp, " MW"
        write(*, '(A,F7.1,A)') "  Max T_fuel = ", sim%max_fuel_temp - 273.15_wp, " °C"
        write(*, '(A,F6.2,A)') "  Avg void = ", sim%avg_void_fraction, " %"
        write(*, '(A,F6.2)') "  Min CHFR = ", sim%min_chfr
        print *, ""
    end subroutine print_transient_summary
    
    !> Write output files
    subroutine write_output(sim, step, output_dir)
        type(simulation_t), intent(in) :: sim
        integer, intent(in) :: step
        character(len=*), intent(in) :: output_dir
        
        character(len=256) :: filename
        integer :: unit, i, j, k
        
        ! Write summary file
        write(filename, '(A,A,I6.6,A)') trim(output_dir), "/summary_", step, ".txt"
        open(newunit=unit, file=trim(filename), status='replace')
        write(unit, *) "Time: ", sim%time
        write(unit, *) "k_eff: ", sim%k_eff
        write(unit, *) "Power: ", sim%power_current
        write(unit, *) "Max_Fuel_Temp: ", sim%max_fuel_temp
        write(unit, *) "Avg_Void: ", sim%avg_void_fraction
        write(unit, *) "Min_CHFR: ", sim%min_chfr
        write(unit, *) "Avg_Burnup: ", sim%avg_burnup
        close(unit)
        
        ! Write 2D power distribution (midplane)
        write(filename, '(A,A,I6.6,A)') trim(output_dir), "/power_", step, ".dat"
        open(newunit=unit, file=trim(filename), status='replace')
        k = sim%nz / 2
        do j = 1, sim%ny
            do i = 1, sim%nx
                write(unit, *) i, j, sim%neutronics%power_density(i, j, k)
            end do
        end do
        close(unit)
    end subroutine write_output
    
    !> Cleanup
    subroutine cleanup_simulation(sim)
        type(simulation_t), intent(inout) :: sim
        
        call mg_destroy(sim%neutronics)
        call heat_destroy(sim%heat)
        call two_phase_destroy(sim%thermalhydraulics)
        call burnup_destroy(sim%burnup)
        call xslib_destroy(sim%xslib)
    end subroutine cleanup_simulation

end program bwr_core_simulator