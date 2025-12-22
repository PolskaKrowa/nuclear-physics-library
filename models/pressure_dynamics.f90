! physics/pressure_dynamics.f90
!
! Pressure dynamics model for nuclear reactor simulation.
! Computes pressure evolution using equation of state and conservation laws.
!
! Supports:
! - Ideal gas equation of state
! - Liquid/vapour phase transitions
! - Pressure wave propagation
! - Compressibility effects
!
! Usage:
!   use pressure_dynamics
!   type(pressure_config_t) :: config
!   type(pressure_state_t) :: state
!   
!   call pressure_init(state, nx, ny, nz, dx, dy, dz, config)
!   call pressure_update_from_temperature(state, temperature, density)
!   call pressure_step(state, dt)
!   call pressure_destroy(state)
!
module pressure_dynamics
    use kinds, only: wp, i32
    use constants, only: R_GAS, TOL_DEFAULT
    implicit none
    private
    
    ! Public interface
    public :: pressure_config_t, pressure_state_t, eos_properties_t
    public :: pressure_init, pressure_destroy
    public :: pressure_update_from_temperature
    public :: pressure_update_from_energy
    public :: pressure_step
    public :: pressure_get_sound_speed
    
    ! Equation of state types
    integer, parameter, public :: EOS_IDEAL_GAS = 1
    integer, parameter, public :: EOS_INCOMPRESSIBLE = 2
    integer, parameter, public :: EOS_STIFFENED_GAS = 3
    integer, parameter, public :: EOS_WATER_STEAM = 4  ! Simplified water/steam
    
    ! Equation of state properties
    type :: eos_properties_t
        integer :: eos_type = EOS_IDEAL_GAS
        real(wp) :: molar_mass = 0.018_wp      ! kg/mol (water)
        real(wp) :: gamma = 1.4_wp             ! Heat capacity ratio (air)
        real(wp) :: p_ref = 101325.0_wp        ! Reference pressure
        real(wp) :: rho_ref = 1000.0_wp        ! Reference density
        real(wp) :: bulk_modulus = 2.2e9_wp    ! K [Pa] for liquids
        real(wp) :: p_infinity = 0.0_wp        ! Stiffened gas parameter
    end type eos_properties_t
    
    ! Configuration
    type :: pressure_config_t
        logical :: include_compressibility = .false.
        logical :: include_phase_change = .false.
        real(wp) :: saturation_temperature = 373.15_wp  ! K (water at 1 atm)
        real(wp) :: latent_heat = 2.26e6_wp             ! J/kg (water)
    end type pressure_config_t
    
    ! Pressure dynamics state
    type :: pressure_state_t
        ! Grid dimensions
        integer :: nx, ny, nz
        real(wp) :: dx, dy, dz
        
        ! Pressure field [Pa]
        real(wp), allocatable :: p(:, :, :)
        real(wp), allocatable :: p_old(:, :, :)
        
        ! Density field [kg/m³]
        real(wp), allocatable :: rho(:, :, :)
        
        ! Temperature field [K]
        real(wp), allocatable :: T(:, :, :)
        
        ! Internal energy [J/kg]
        real(wp), allocatable :: e(:, :, :)
        
        ! Equation of state properties
        type(eos_properties_t), allocatable :: eos(:, :, :)
        
        ! Velocity for acoustic wave propagation
        real(wp), allocatable :: vx(:, :, :)
        real(wp), allocatable :: vy(:, :, :)
        real(wp), allocatable :: vz(:, :, :)
        
        ! Configuration
        type(pressure_config_t) :: config
        
        ! Time tracking
        real(wp) :: time = 0.0_wp
        integer :: steps = 0
    end type pressure_state_t
    
contains

    !> Initialise pressure dynamics state
    subroutine pressure_init(state, nx, ny, nz, dx, dy, dz, config)
        type(pressure_state_t), intent(out) :: state
        integer, intent(in) :: nx, ny, nz
        real(wp), intent(in) :: dx, dy, dz
        type(pressure_config_t), intent(in), optional :: config
        
        state%nx = nx
        state%ny = ny
        state%nz = nz
        state%dx = dx
        state%dy = dy
        state%dz = dz
        
        ! Allocate fields
        allocate(state%p(nx, ny, nz))
        allocate(state%p_old(nx, ny, nz))
        allocate(state%rho(nx, ny, nz))
        allocate(state%T(nx, ny, nz))
        allocate(state%e(nx, ny, nz))
        allocate(state%eos(nx, ny, nz))
        
        if (config%include_compressibility) then
            allocate(state%vx(nx, ny, nz))
            allocate(state%vy(nx, ny, nz))
            allocate(state%vz(nx, ny, nz))
            state%vx = 0.0_wp
            state%vy = 0.0_wp
            state%vz = 0.0_wp
        end if
        
        ! Initialise fields
        state%p = 101325.0_wp       ! Atmospheric pressure
        state%p_old = 101325.0_wp
        state%rho = 1000.0_wp        ! Water density
        state%T = 300.0_wp           ! Room temperature
        state%e = 0.0_wp
        
        ! Default EOS (liquid water)
        state%eos(:, :, :)%eos_type = EOS_INCOMPRESSIBLE
        state%eos(:, :, :)%rho_ref = 1000.0_wp
        state%eos(:, :, :)%bulk_modulus = 2.2e9_wp
        
        ! Configuration
        if (present(config)) then
            state%config = config
        end if
    end subroutine pressure_init
    
    !> Destroy pressure dynamics state
    subroutine pressure_destroy(state)
        type(pressure_state_t), intent(inout) :: state
        
        if (allocated(state%p)) deallocate(state%p)
        if (allocated(state%p_old)) deallocate(state%p_old)
        if (allocated(state%rho)) deallocate(state%rho)
        if (allocated(state%T)) deallocate(state%T)
        if (allocated(state%e)) deallocate(state%e)
        if (allocated(state%eos)) deallocate(state%eos)
        if (allocated(state%vx)) deallocate(state%vx)
        if (allocated(state%vy)) deallocate(state%vy)
        if (allocated(state%vz)) deallocate(state%vz)
    end subroutine pressure_destroy
    
    !> Update pressure from temperature and density using EOS
    subroutine pressure_update_from_temperature(state, T, rho)
        type(pressure_state_t), intent(inout) :: state
        real(wp), intent(in), optional :: T(:, :, :)
        real(wp), intent(in), optional :: rho(:, :, :)
        
        integer :: i, j, k
        
        if (present(T)) state%T = T
        if (present(rho)) state%rho = rho
        
        ! Apply equation of state
        do k = 1, state%nz
            do j = 1, state%ny
                do i = 1, state%nx
                    state%p(i, j, k) = compute_pressure_eos( &
                        state%T(i, j, k), state%rho(i, j, k), state%eos(i, j, k))
                end do
            end do
        end do
    end subroutine pressure_update_from_temperature
    
    !> Update pressure from internal energy
    subroutine pressure_update_from_energy(state, e, rho)
        type(pressure_state_t), intent(inout) :: state
        real(wp), intent(in) :: e(:, :, :)
        real(wp), intent(in), optional :: rho(:, :, :)
        
        integer :: i, j, k
        real(wp) :: cv, T_local
        
        state%e = e
        if (present(rho)) state%rho = rho
        
        do k = 1, state%nz
            do j = 1, state%ny
                do i = 1, state%nx
                    ! Convert internal energy to temperature
                    cv = compute_cv(state%eos(i, j, k))
                    T_local = state%e(i, j, k) / cv
                    state%T(i, j, k) = T_local
                    
                    ! Apply EOS
                    state%p(i, j, k) = compute_pressure_eos( &
                        T_local, state%rho(i, j, k), state%eos(i, j, k))
                end do
            end do
        end do
    end subroutine pressure_update_from_energy
    
    !> Time step for compressible flow (acoustic wave propagation)
    subroutine pressure_step(state, dt)
        type(pressure_state_t), intent(inout) :: state
        real(wp), intent(in) :: dt
        
        if (.not. state%config%include_compressibility) then
            ! Incompressible: pressure determined by incompressibility constraint
            ! (handled by fluid solver)
            state%time = state%time + dt
            state%steps = state%steps + 1
            return
        end if
        
        ! Compressible flow: solve acoustic equations
        call solve_acoustic_equations(state, dt)
        
        state%time = state%time + dt
        state%steps = state%steps + 1
    end subroutine pressure_step
    
    !> Solve acoustic wave equations
    subroutine solve_acoustic_equations(state, dt)
        type(pressure_state_t), intent(inout) :: state
        real(wp), intent(in) :: dt
        
        real(wp), allocatable :: dp_dx(:, :, :), dp_dy(:, :, :), dp_dz(:, :, :)
        real(wp), allocatable :: div_v(:, :, :)
        integer :: i, j, k
        real(wp) :: c2, inv_dx, inv_dy, inv_dz
        
        allocate(dp_dx(state%nx, state%ny, state%nz))
        allocate(dp_dy(state%nx, state%ny, state%nz))
        allocate(dp_dz(state%nx, state%ny, state%nz))
        allocate(div_v(state%nx, state%ny, state%nz))
        
        state%p_old = state%p
        
        inv_dx = 1.0_wp / state%dx
        inv_dy = 1.0_wp / state%dy
        inv_dz = 1.0_wp / state%dz
        
        ! Compute pressure gradients
        do k = 2, state%nz - 1
            do j = 2, state%ny - 1
                do i = 2, state%nx - 1
                    dp_dx(i, j, k) = (state%p(i + 1, j, k) - state%p(i - 1, j, k)) * &
                        0.5_wp * inv_dx
                    dp_dy(i, j, k) = (state%p(i, j + 1, k) - state%p(i, j - 1, k)) * &
                        0.5_wp * inv_dy
                    dp_dz(i, j, k) = (state%p(i, j, k + 1) - state%p(i, j, k - 1)) * &
                        0.5_wp * inv_dz
                end do
            end do
        end do
        
        ! Momentum equation: ∂v/∂t = -1/ρ · ∇p
        do k = 2, state%nz - 1
            do j = 2, state%ny - 1
                do i = 2, state%nx - 1
                    state%vx(i, j, k) = state%vx(i, j, k) - &
                        dt * dp_dx(i, j, k) / state%rho(i, j, k)
                    state%vy(i, j, k) = state%vy(i, j, k) - &
                        dt * dp_dy(i, j, k) / state%rho(i, j, k)
                    state%vz(i, j, k) = state%vz(i, j, k) - &
                        dt * dp_dz(i, j, k) / state%rho(i, j, k)
                end do
            end do
        end do
        
        ! Compute velocity divergence
        do k = 2, state%nz - 1
            do j = 2, state%ny - 1
                do i = 2, state%nx - 1
                    div_v(i, j, k) = &
                        (state%vx(i + 1, j, k) - state%vx(i - 1, j, k)) * 0.5_wp * inv_dx + &
                        (state%vy(i, j + 1, k) - state%vy(i, j - 1, k)) * 0.5_wp * inv_dy + &
                        (state%vz(i, j, k + 1) - state%vz(i, j, k - 1)) * 0.5_wp * inv_dz
                end do
            end do
        end do
        
        ! Continuity equation: ∂p/∂t = -ρ·c² · ∇·v
        do k = 2, state%nz - 1
            do j = 2, state%ny - 1
                do i = 2, state%nx - 1
                    c2 = pressure_get_sound_speed(state, i, j, k)**2
                    state%p(i, j, k) = state%p(i, j, k) - &
                        dt * state%rho(i, j, k) * c2 * div_v(i, j, k)
                end do
            end do
        end do
        
        deallocate(dp_dx, dp_dy, dp_dz, div_v)
    end subroutine solve_acoustic_equations
    
    !> Compute pressure from temperature and density using EOS
    function compute_pressure_eos(T, rho, eos) result(p)
        real(wp), intent(in) :: T, rho
        type(eos_properties_t), intent(in) :: eos
        real(wp) :: p
        
        real(wp) :: R_specific
        
        select case(eos%eos_type)
        case(EOS_IDEAL_GAS)
            ! p = ρ·R·T (R = R_gas/M)
            R_specific = R_GAS / eos%molar_mass
            p = rho * R_specific * T
            
        case(EOS_INCOMPRESSIBLE)
            ! Tait equation: p = p_ref + K·((ρ/ρ_ref)^γ - 1)
            p = eos%p_ref + eos%bulk_modulus * &
                ((rho / eos%rho_ref)**eos%gamma - 1.0_wp)
            
        case(EOS_STIFFENED_GAS)
            ! p = (γ-1)·ρ·e - γ·p_∞
            ! Simplified version assuming constant e
            R_specific = R_GAS / eos%molar_mass
            p = (eos%gamma - 1.0_wp) * rho * R_specific * T - &
                eos%gamma * eos%p_infinity
            
        case(EOS_WATER_STEAM)
            ! Simplified water/steam model
            if (T < 373.15_wp) then
                ! Liquid phase (incompressible)
                p = eos%p_ref + eos%bulk_modulus * (rho / eos%rho_ref - 1.0_wp)
            else
                ! Vapour phase (ideal gas)
                R_specific = R_GAS / eos%molar_mass
                p = rho * R_specific * T
            end if
            
        case default
            ! Default to ideal gas
            R_specific = R_GAS / eos%molar_mass
            p = rho * R_specific * T
        end select
    end function compute_pressure_eos
    
    !> Get local sound speed
    function pressure_get_sound_speed(state, i, j, k) result(c)
        type(pressure_state_t), intent(in) :: state
        integer, intent(in) :: i, j, k
        real(wp) :: c
        
        real(wp) :: gamma, K, rho, R_specific
        
        select case(state%eos(i, j, k)%eos_type)
        case(EOS_IDEAL_GAS)
            ! c = sqrt(γ·R·T)
            gamma = state%eos(i, j, k)%gamma
            R_specific = R_GAS / state%eos(i, j, k)%molar_mass
            c = sqrt(gamma * R_specific * state%T(i, j, k))
            
        case(EOS_INCOMPRESSIBLE)
            ! c = sqrt(K/ρ)
            K = state%eos(i, j, k)%bulk_modulus
            rho = state%rho(i, j, k)
            c = sqrt(K / rho)
            
        case(EOS_STIFFENED_GAS)
            gamma = state%eos(i, j, k)%gamma
            R_specific = R_GAS / state%eos(i, j, k)%molar_mass
            c = sqrt(gamma * R_specific * state%T(i, j, k))
            
        case(EOS_WATER_STEAM)
            if (state%T(i, j, k) < 373.15_wp) then
                ! Liquid
                K = state%eos(i, j, k)%bulk_modulus
                rho = state%rho(i, j, k)
                c = sqrt(K / rho)
            else
                ! Vapour
                gamma = state%eos(i, j, k)%gamma
                R_specific = R_GAS / state%eos(i, j, k)%molar_mass
                c = sqrt(gamma * R_specific * state%T(i, j, k))
            end if
            
        case default
            c = 1500.0_wp  ! Default (water at room temp)
        end select
        
        ! Ensure positive, finite value
        c = max(c, 100.0_wp)
        c = min(c, 10000.0_wp)
    end function pressure_get_sound_speed
    
    !> Compute specific heat at constant volume
    function compute_cv(eos) result(cv)
        type(eos_properties_t), intent(in) :: eos
        real(wp) :: cv
        
        real(wp) :: R_specific
        
        select case(eos%eos_type)
        case(EOS_IDEAL_GAS)
            ! cv = R/(γ-1)
            R_specific = R_GAS / eos%molar_mass
            cv = R_specific / (eos%gamma - 1.0_wp)
            
        case(EOS_WATER_STEAM)
            ! Water: approximately 4180 J/kg·K
            cv = 4180.0_wp
            
        case default
            cv = 4180.0_wp  ! Default to water
        end select
    end function compute_cv

end module pressure_dynamics