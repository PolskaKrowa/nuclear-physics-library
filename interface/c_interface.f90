! interface/c_interface.f90
!
! C interface for the nuclear physics library
! Provides ISO_C_BINDING compatible wrappers for C/C++ programs
!
! Usage from C:
!   #include "nuclear_physics.h"
!   void* reactor = reactor_create(nx, ny, nz, dx, dy, dz);
!   reactor_set_fuel_region(reactor, i1, i2, j1, j2, k1, k2, power);
!   reactor_step(reactor, dt);
!   reactor_get_temperature(reactor, T_array);
!   reactor_destroy(reactor);
!
module c_interface
    use iso_c_binding
    use kinds, only: wp
    use heat_transfer
    use fluid_dynamics
    use nuclear_fission
    use pressure_dynamics
    implicit none
    private
    
    ! Reactor state combining all physics
    type :: reactor_state_t
        type(heat_state_t) :: heat
        type(fluid_state_t) :: fluid
        type(fission_state_t) :: fission
        type(pressure_state_t) :: pressure
        logical :: initialized = .false.
        integer :: nx, ny, nz
    end type reactor_state_t
    
contains

    !> Create reactor simulation
    function reactor_create(nx, ny, nz, dx, dy, dz) result(handle) bind(C, name="reactor_create")
        integer(c_int), value :: nx, ny, nz
        real(c_double), value :: dx, dy, dz
        type(c_ptr) :: handle
        
        type(reactor_state_t), pointer :: reactor
        type(heat_config_t) :: heat_cfg
        type(fluid_config_t) :: fluid_cfg
        type(fission_config_t) :: fission_cfg
        type(pressure_config_t) :: pressure_cfg
        
        allocate(reactor)
        
        reactor%nx = nx
        reactor%ny = ny
        reactor%nz = nz
        
        ! Configure for stability
        fission_cfg%use_point_kinetics = .true.
        fission_cfg%power_normalisation = 1.0e6_wp  ! 1 MW nominal
        
        ! Initialize all physics modules
        call heat_init(reactor%heat, nx, ny, nz, dx, dy, dz, heat_cfg)
        call fluid_init(reactor%fluid, nx, ny, nz, dx, dy, dz, fluid_cfg)
        call fission_init(reactor%fission, nx, ny, nz, dx, dy, dz, fission_cfg)
        call pressure_init(reactor%pressure, nx, ny, nz, dx, dy, dz, pressure_cfg)
        
        ! Start at critical (zero reactivity)
        call fission_set_reactivity(reactor%fission, 0.0_wp)
        
        reactor%initialized = .true.
        handle = c_loc(reactor)
    end function reactor_create
    
    !> Destroy reactor simulation
    subroutine reactor_destroy(handle) bind(C, name="reactor_destroy")
        type(c_ptr), value :: handle
        type(reactor_state_t), pointer :: reactor
        
        call c_f_pointer(handle, reactor)
        if (.not. reactor%initialized) return
        
        call heat_destroy(reactor%heat)
        call fluid_destroy(reactor%fluid)
        call fission_destroy(reactor%fission)
        call pressure_destroy(reactor%pressure)
        
        deallocate(reactor)
    end subroutine reactor_destroy
    
    !> Set fuel region with nuclear properties
    subroutine reactor_set_fuel_region(handle, i1, i2, j1, j2, k1, k2, &
                                       power_density, enrichment) &
        bind(C, name="reactor_set_fuel_region")
        type(c_ptr), value :: handle
        integer(c_int), value :: i1, i2, j1, j2, k1, k2
        real(c_double), value :: power_density, enrichment
        type(reactor_state_t), pointer :: reactor
        
        real(wp) :: sigma_f, sigma_a, nu, kappa
        real(wp) :: k_fuel, rho_fuel, cp_fuel
        
        call c_f_pointer(handle, reactor)
        
        ! Typical UO2 fuel properties
        k_fuel = 3.5_wp      ! W/m·K
        rho_fuel = 10970.0_wp  ! kg/m³
        cp_fuel = 300.0_wp   ! J/kg·K
        
        ! Nuclear cross sections (simplified, enrichment-dependent)
        sigma_f = 0.008_wp * enrichment  ! Fission cross section
        sigma_a = 0.010_wp * enrichment  ! Absorption
        nu = 2.43_wp                     ! Neutrons per fission
        kappa = 200.0_wp                 ! MeV per fission
        
        ! Set material properties
        call heat_set_properties(reactor%heat, k_fuel, rho_fuel, cp_fuel, &
                                i1, i2, j1, j2, k1, k2)
        call heat_set_source(reactor%heat, power_density, &
                            i1, i2, j1, j2, k1, k2)
        
        call fission_set_cross_sections(reactor%fission, sigma_f, sigma_a, &
                                        nu, kappa, i1, i2, j1, j2, k1, k2)
        
        call fluid_set_properties(reactor%fluid, rho_fuel, 1.0e-3_wp, 0.0_wp, &
                                 i1, i2, j1, j2, k1, k2)
    end subroutine reactor_set_fuel_region
    
    !> Set coolant (water) region
    subroutine reactor_set_coolant_region(handle, i1, i2, j1, j2, k1, k2, &
                                         inlet_temp, mass_flow) &
        bind(C, name="reactor_set_coolant_region")
        type(c_ptr), value :: handle
        integer(c_int), value :: i1, i2, j1, j2, k1, k2
        real(c_double), value :: inlet_temp, mass_flow
        type(reactor_state_t), pointer :: reactor
        
        real(wp) :: k_water, rho_water, cp_water, mu_water
        
        call c_f_pointer(handle, reactor)
        
        ! Water properties at typical PWR/BWR conditions
        k_water = 0.6_wp       ! W/m·K
        rho_water = 700.0_wp   ! kg/m³ (subcooled)
        cp_water = 4180.0_wp   ! J/kg·K
        mu_water = 1.0e-4_wp   ! Pa·s
        
        call heat_set_properties(reactor%heat, k_water, rho_water, cp_water, &
                                i1, i2, j1, j2, k1, k2)
        
        call fluid_set_properties(reactor%fluid, rho_water, mu_water, 0.0_wp, &
                                 i1, i2, j1, j2, k1, k2)
        
        ! Set inlet temperature
        reactor%heat%T(i1:i2, j1:j2, k1:k2) = inlet_temp
    end subroutine reactor_set_coolant_region
    
    !> Set control rod position (affects reactivity)
    subroutine reactor_set_control_rods(handle, insertion_fraction) &
        bind(C, name="reactor_set_control_rods")
        type(c_ptr), value :: handle
        real(c_double), value :: insertion_fraction
        type(reactor_state_t), pointer :: reactor
        
        real(wp) :: reactivity, rho_worth
        
        call c_f_pointer(handle, reactor)
        
        ! Control rod worth: -$10 at full insertion, 0 at full withdrawal
        ! Reactivity swing: 0 to -10 dollars
        rho_worth = -10.0_wp * 0.0065_wp  ! -$10 in absolute units
        
        ! Map insertion fraction to reactivity
        ! 0.0 = withdrawn (slightly positive for power)
        ! 1.0 = inserted (very negative)
        reactivity = 0.001_wp + rho_worth * insertion_fraction
        
        ! Clamp to safe bounds
        reactivity = max(-0.065_wp, min(0.002_wp, reactivity))
        
        call fission_set_reactivity(reactor%fission, reactivity)
    end subroutine reactor_set_control_rods
    
    !> Perform one time step
    subroutine reactor_step(handle, dt) bind(C, name="reactor_step")
        type(c_ptr), value :: handle
        real(c_double), value :: dt
        type(reactor_state_t), pointer :: reactor
        
        real(wp), allocatable :: power_density(:,:,:)
        
        call c_f_pointer(handle, reactor)
        
        allocate(power_density(reactor%nx, reactor%ny, reactor%nz))
        
        ! Update fission power
        call fission_step(reactor%fission, dt)
        call fission_get_power_density(reactor%fission, power_density)
        
        ! Update heat transfer with fission power source
        reactor%heat%Q = power_density
        call heat_step(reactor%heat, dt)
        
        ! Update fluid dynamics with temperature field
        call fluid_set_temperature(reactor%fluid, reactor%heat%T)
        call fluid_step(reactor%fluid, dt)
        
        ! Update pressure based on temperature and density
        call pressure_update_from_temperature(reactor%pressure, &
                                             reactor%heat%T, &
                                             reactor%fluid%properties(:,:,:)%density)
        
        deallocate(power_density)
    end subroutine reactor_step
    
    !> Get temperature field
    subroutine reactor_get_temperature(handle, T_out) bind(C, name="reactor_get_temperature")
        type(c_ptr), value :: handle
        real(c_double), intent(out) :: T_out(*)
        type(reactor_state_t), pointer :: reactor
        
        integer :: i, j, k, idx
        
        call c_f_pointer(handle, reactor)
        
        idx = 1
        do k = 1, reactor%nz
            do j = 1, reactor%ny
                do i = 1, reactor%nx
                    T_out(idx) = reactor%heat%T(i, j, k)
                    idx = idx + 1
                end do
            end do
        end do
    end subroutine reactor_get_temperature
    
    !> Get pressure field
    subroutine reactor_get_pressure(handle, p_out) bind(C, name="reactor_get_pressure")
        type(c_ptr), value :: handle
        real(c_double), intent(out) :: p_out(*)
        type(reactor_state_t), pointer :: reactor
        
        integer :: i, j, k, idx
        
        call c_f_pointer(handle, reactor)
        
        idx = 1
        do k = 1, reactor%nz
            do j = 1, reactor%ny
                do i = 1, reactor%nx
                    p_out(idx) = reactor%pressure%p(i, j, k)
                    idx = idx + 1
                end do
            end do
        end do
    end subroutine reactor_get_pressure
    
    !> Get power density field
    subroutine reactor_get_power(handle, power_out) bind(C, name="reactor_get_power")
        type(c_ptr), value :: handle
        real(c_double), intent(out) :: power_out(*)
        type(reactor_state_t), pointer :: reactor
        
        integer :: i, j, k, idx
        
        call c_f_pointer(handle, reactor)
        
        idx = 1
        do k = 1, reactor%nz
            do j = 1, reactor%ny
                do i = 1, reactor%nx
                    power_out(idx) = reactor%fission%power_density(i, j, k)
                    idx = idx + 1
                end do
            end do
        end do
    end subroutine reactor_get_power
    
    !> Get total reactor power
    function reactor_get_total_power(handle) result(power) bind(C, name="reactor_get_total_power")
        type(c_ptr), value :: handle
        real(c_double) :: power
        type(reactor_state_t), pointer :: reactor
        
        call c_f_pointer(handle, reactor)
        power = fission_get_total_power(reactor%fission)
    end function reactor_get_total_power
    
    !> Get velocity field components
    subroutine reactor_get_velocity(handle, vx_out, vy_out, vz_out) &
        bind(C, name="reactor_get_velocity")
        type(c_ptr), value :: handle
        real(c_double), intent(out) :: vx_out(*), vy_out(*), vz_out(*)
        type(reactor_state_t), pointer :: reactor
        
        integer :: i, j, k, idx
        
        call c_f_pointer(handle, reactor)
        
        idx = 1
        do k = 1, reactor%nz
            do j = 1, reactor%ny
                do i = 1, reactor%nx
                    ! Average staggered velocities to cell centers
                    vx_out(idx) = 0.5_wp * (reactor%fluid%vx(i, j, k) + &
                                           reactor%fluid%vx(i+1, j, k))
                    vy_out(idx) = 0.5_wp * (reactor%fluid%vy(i, j, k) + &
                                           reactor%fluid%vy(i, j+1, k))
                    vz_out(idx) = 0.5_wp * (reactor%fluid%vz(i, j, k) + &
                                           reactor%fluid%vz(i, j, k+1))
                    idx = idx + 1
                end do
            end do
        end do
    end subroutine reactor_get_velocity
    
    !> Get maximum stable time step
    function reactor_get_max_dt(handle) result(dt_max) bind(C, name="reactor_get_max_dt")
        type(c_ptr), value :: handle
        real(c_double) :: dt_max
        type(reactor_state_t), pointer :: reactor
        
        real(wp) :: dt_heat, dt_fluid, dt_fission
        
        call c_f_pointer(handle, reactor)
        
        dt_heat = heat_get_max_dt(reactor%heat)
        dt_fluid = fluid_get_max_dt(reactor%fluid)
        dt_fission = 0.01_wp  ! Typical neutronic time scale
        
        dt_max = min(dt_heat, dt_fluid, dt_fission)
    end function reactor_get_max_dt

end module c_interface