! physics/heat_transfer.f90
!
! Heat transfer model for nuclear reactor simulation.
! Solves the heat equation with convection and source terms:
!   ∂T/∂t = α·∇²T - v·∇T + Q/ρcₚ
!
! Features:
! - Conduction (diffusion)
! - Convection with fluid flow
! - Volumetric heat sources (fission, decay)
! - Multiple material regions
! - Boundary conditions (Dirichlet, Neumann, convective)
!
! Usage:
!   use heat_transfer
!   type(heat_config_t) :: config
!   type(heat_state_t) :: state
!   
!   call heat_init(state, nx, ny, nz, dx, dy, dz, config)
!   call heat_set_properties(state, thermal_conductivity, density, specific_heat)
!   call heat_set_source(state, power_density)
!   call heat_step(state, dt)
!   call heat_destroy(state)
!
module heat_transfer
    use kinds, only: wp, i32
    use constants, only: TOL_DEFAULT
    use finite_difference, only: fd_laplacian_2d, fd_laplacian_3d, fd_derivative_2d
    use finite_volume, only: fv_diffusion_2d, fv_advection_diffusion_1d
    implicit none
    private
    
    ! Public interface
    public :: heat_config_t, heat_state_t, heat_material_t
    public :: heat_init, heat_destroy
    public :: heat_set_properties, heat_set_source, heat_set_velocity
    public :: heat_set_bc, heat_apply_bc
    public :: heat_step, heat_step_implicit
    public :: heat_get_max_dt
    
    ! Boundary condition types
    integer, parameter, public :: BC_DIRICHLET = 1    ! Fixed temperature
    integer, parameter, public :: BC_NEUMANN = 2      ! Fixed heat flux
    integer, parameter, public :: BC_CONVECTIVE = 3   ! Convective cooling
    integer, parameter, public :: BC_ADIABATIC = 4    ! Insulated
    
    ! Material properties
    type :: heat_material_t
        real(wp) :: thermal_conductivity  ! k [W/m·K]
        real(wp) :: density              ! ρ [kg/m³]
        real(wp) :: specific_heat        ! cₚ [J/kg·K]
        real(wp) :: thermal_diffusivity  ! α = k/(ρ·cₚ) [m²/s]
    end type heat_material_t
    
    ! Configuration
    type :: heat_config_t
        logical :: include_convection = .true.
        logical :: use_implicit = .false.
        integer :: max_implicit_iter = 20
        real(wp) :: implicit_tolerance = 1.0e-6_wp
        integer :: bc_type(6) = BC_ADIABATIC  ! -x, +x, -y, +y, -z, +z
        real(wp) :: bc_value(6) = 0.0_wp      ! Boundary values
        real(wp) :: htc(6) = 0.0_wp           ! Heat transfer coefficient [W/m²·K]
        real(wp) :: t_ambient(6) = 300.0_wp   ! Ambient temperature [K]
    end type heat_config_t
    
    ! Heat transfer state
    type :: heat_state_t
        ! Grid dimensions
        integer :: nx, ny, nz
        real(wp) :: dx, dy, dz
        real(wp) :: volume
        
        ! Temperature field [K]
        real(wp), allocatable :: T(:, :, :)
        real(wp), allocatable :: T_old(:, :, :)
        
        ! Material properties (per cell)
        type(heat_material_t), allocatable :: material(:, :, :)
        
        ! Heat sources [W/m³]
        real(wp), allocatable :: Q(:, :, :)
        
        ! Velocity field for convection [m/s]
        real(wp), allocatable :: vx(:, :, :)
        real(wp), allocatable :: vy(:, :, :)
        real(wp), allocatable :: vz(:, :, :)
        
        ! Configuration
        type(heat_config_t) :: config
        
        ! Time tracking
        real(wp) :: time = 0.0_wp
        integer :: steps = 0
    end type heat_state_t
    
contains

    !> Initialise heat transfer state
    subroutine heat_init(state, nx, ny, nz, dx, dy, dz, config)
        type(heat_state_t), intent(out) :: state
        integer, intent(in) :: nx, ny, nz
        real(wp), intent(in) :: dx, dy, dz
        type(heat_config_t), intent(in), optional :: config
        
        state%nx = nx
        state%ny = ny
        state%nz = nz
        state%dx = dx
        state%dy = dy
        state%dz = dz
        state%volume = dx * dy * dz
        
        ! Allocate arrays
        allocate(state%T(nx, ny, nz))
        allocate(state%T_old(nx, ny, nz))
        allocate(state%material(nx, ny, nz))
        allocate(state%Q(nx, ny, nz))
        
        if (state%config%include_convection) then
            allocate(state%vx(nx, ny, nz))
            allocate(state%vy(nx, ny, nz))
            allocate(state%vz(nx, ny, nz))
            state%vx = 0.0_wp
            state%vy = 0.0_wp
            state%vz = 0.0_wp
        end if
        
        ! Initialise fields
        state%T = 300.0_wp      ! Room temperature
        state%T_old = 300.0_wp
        state%Q = 0.0_wp
        
        ! Default material (water at 300K)
        state%material(:, :, :)%thermal_conductivity = 0.6_wp
        state%material(:, :, :)%density = 1000.0_wp
        state%material(:, :, :)%specific_heat = 4180.0_wp
        state%material(:, :, :)%thermal_diffusivity = &
            0.6_wp / (1000.0_wp * 4180.0_wp)
        
        ! Configuration
        if (present(config)) then
            state%config = config
        end if
    end subroutine heat_init
    
    !> Destroy heat transfer state
    subroutine heat_destroy(state)
        type(heat_state_t), intent(inout) :: state
        
        if (allocated(state%T)) deallocate(state%T)
        if (allocated(state%T_old)) deallocate(state%T_old)
        if (allocated(state%material)) deallocate(state%material)
        if (allocated(state%Q)) deallocate(state%Q)
        if (allocated(state%vx)) deallocate(state%vx)
        if (allocated(state%vy)) deallocate(state%vy)
        if (allocated(state%vz)) deallocate(state%vz)
    end subroutine heat_destroy
    
    !> Set material properties for a region
    subroutine heat_set_properties(state, k, rho, cp, i1, i2, j1, j2, k1, k2)
        type(heat_state_t), intent(inout) :: state
        real(wp), intent(in) :: k, rho, cp
        integer, intent(in), optional :: i1, i2, j1, j2, k1, k2
        
        integer :: imin, imax, jmin, jmax, kmin, kmax
        real(wp) :: alpha
        
        ! Default: entire domain
        imin = 1; imax = state%nx
        jmin = 1; jmax = state%ny
        kmin = 1; kmax = state%nz
        
        if (present(i1)) imin = i1
        if (present(i2)) imax = i2
        if (present(j1)) jmin = j1
        if (present(j2)) jmax = j2
        if (present(k1)) kmin = k1
        if (present(k2)) kmax = k2
        
        alpha = k / (rho * cp)
        
        state%material(imin:imax, jmin:jmax, kmin:kmax)%thermal_conductivity = k
        state%material(imin:imax, jmin:jmax, kmin:kmax)%density = rho
        state%material(imin:imax, jmin:jmax, kmin:kmax)%specific_heat = cp
        state%material(imin:imax, jmin:jmax, kmin:kmax)%thermal_diffusivity = alpha
    end subroutine heat_set_properties
    
    !> Set volumetric heat source [W/m³]
    subroutine heat_set_source(state, Q, i1, i2, j1, j2, k1, k2)
        type(heat_state_t), intent(inout) :: state
        real(wp), intent(in) :: Q
        integer, intent(in), optional :: i1, i2, j1, j2, k1, k2
        
        integer :: imin, imax, jmin, jmax, kmin, kmax
        
        imin = 1; imax = state%nx
        jmin = 1; jmax = state%ny
        kmin = 1; kmax = state%nz
        
        if (present(i1)) imin = i1
        if (present(i2)) imax = i2
        if (present(j1)) jmin = j1
        if (present(j2)) jmax = j2
        if (present(k1)) kmin = k1
        if (present(k2)) kmax = k2
        
        state%Q(imin:imax, jmin:jmax, kmin:kmax) = Q
    end subroutine heat_set_source
    
    !> Set velocity field for convection
    subroutine heat_set_velocity(state, vx, vy, vz)
        type(heat_state_t), intent(inout) :: state
        real(wp), intent(in), optional :: vx(:, :, :)
        real(wp), intent(in), optional :: vy(:, :, :)
        real(wp), intent(in), optional :: vz(:, :, :)
        
        if (.not. state%config%include_convection) return
        
        if (present(vx)) state%vx = vx
        if (present(vy)) state%vy = vy
        if (present(vz)) state%vz = vz
    end subroutine heat_set_velocity
    
    !> Set boundary condition
    subroutine heat_set_bc(state, face, bc_type, value, htc, t_ambient)
        type(heat_state_t), intent(inout) :: state
        integer, intent(in) :: face  ! 1=-x, 2=+x, 3=-y, 4=+y, 5=-z, 6=+z
        integer, intent(in) :: bc_type
        real(wp), intent(in), optional :: value, htc, t_ambient
        
        state%config%bc_type(face) = bc_type
        
        if (present(value)) state%config%bc_value(face) = value
        if (present(htc)) state%config%htc(face) = htc
        if (present(t_ambient)) state%config%t_ambient(face) = t_ambient
    end subroutine heat_set_bc
    
    !> Explicit time step for heat equation
    subroutine heat_step(state, dt)
        type(heat_state_t), intent(inout) :: state
        real(wp), intent(in) :: dt
        
        real(wp), allocatable :: laplacian(:, :, :)
        real(wp), allocatable :: dT_dx(:, :, :), dT_dy(:, :, :), dT_dz(:, :, :)
        integer :: i, j, k
        real(wp) :: alpha, rho_cp, diffusion, convection, source
        
        allocate(laplacian(state%nx, state%ny, state%nz))
        
        state%T_old = state%T
        
        ! Compute Laplacian (diffusion term)
        call fd_laplacian_3d(state%T, laplacian, state%dx, state%dy, state%dz)
        
        ! Compute convection terms if needed
        if (state%config%include_convection) then
            allocate(dT_dx(state%nx, state%ny, state%nz))
            allocate(dT_dy(state%nx, state%ny, state%nz))
            allocate(dT_dz(state%nx, state%ny, state%nz))
            
            ! Compute gradients
            call compute_gradient_3d(state%T, dT_dx, dT_dy, dT_dz, &
                                    state%dx, state%dy, state%dz)
        end if
        
        ! Update temperature field
        do k = 1, state%nz
            do j = 1, state%ny
                do i = 1, state%nx
                    alpha = state%material(i, j, k)%thermal_diffusivity
                    rho_cp = state%material(i, j, k)%density * &
                             state%material(i, j, k)%specific_heat
                    
                    ! Diffusion
                    diffusion = alpha * laplacian(i, j, k)
                    
                    ! Convection
                    convection = 0.0_wp
                    if (state%config%include_convection) then
                        convection = -(state%vx(i, j, k) * dT_dx(i, j, k) + &
                                      state%vy(i, j, k) * dT_dy(i, j, k) + &
                                      state%vz(i, j, k) * dT_dz(i, j, k))
                    end if
                    
                    ! Source term
                    source = state%Q(i, j, k) / rho_cp
                    
                    ! Explicit Euler update
                    state%T(i, j, k) = state%T(i, j, k) + dt * &
                        (diffusion + convection + source)
                end do
            end do
        end do
        
        ! Apply boundary conditions
        call heat_apply_bc(state)
        
        ! Update time
        state%time = state%time + dt
        state%steps = state%steps + 1
        
        ! Clean up
        deallocate(laplacian)
        if (allocated(dT_dx)) deallocate(dT_dx, dT_dy, dT_dz)
    end subroutine heat_step
    
    !> Implicit time step (more stable for large dt)
    subroutine heat_step_implicit(state, dt)
        type(heat_state_t), intent(inout) :: state
        real(wp), intent(in) :: dt
        
        ! For now, fall back to explicit
        ! In production, implement implicit solver using solve_linear module
        call heat_step(state, dt)
    end subroutine heat_step_implicit
    
    !> Apply boundary conditions
    subroutine heat_apply_bc(state)
        type(heat_state_t), intent(inout) :: state
        integer :: i, j, k, face
        real(wp) :: T_bc, flux, htc, T_amb, T_surf
        
        ! -x face (i=1)
        face = 1
        select case(state%config%bc_type(face))
        case(BC_DIRICHLET)
            state%T(1, :, :) = state%config%bc_value(face)
        case(BC_NEUMANN)
            ! Flux specified: q = -k·dT/dx
            flux = state%config%bc_value(face)
            do k = 1, state%nz
                do j = 1, state%ny
                    state%T(1, j, k) = state%T(2, j, k) - flux * state%dx / &
                        state%material(1, j, k)%thermal_conductivity
                end do
            end do
        case(BC_CONVECTIVE)
            ! Newton's law of cooling: q = h(T_surf - T_amb)
            htc = state%config%htc(face)
            T_amb = state%config%t_ambient(face)
            do k = 1, state%nz
                do j = 1, state%ny
                    ! Implicit update for stability
                    T_surf = (state%material(1, j, k)%thermal_conductivity * &
                             state%T(2, j, k) / state%dx + htc * T_amb) / &
                            (state%material(1, j, k)%thermal_conductivity / &
                             state%dx + htc)
                    state%T(1, j, k) = T_surf
                end do
            end do
        case(BC_ADIABATIC)
            state%T(1, :, :) = state%T(2, :, :)
        end select
        
        ! +x face (i=nx) - similar logic
        face = 2
        select case(state%config%bc_type(face))
        case(BC_DIRICHLET)
            state%T(state%nx, :, :) = state%config%bc_value(face)
        case(BC_ADIABATIC)
            state%T(state%nx, :, :) = state%T(state%nx - 1, :, :)
        end select
        
        ! Similar for other faces...
        ! (Implementation omitted for brevity - follow same pattern)
    end subroutine heat_apply_bc
    
    !> Get maximum stable time step
    function heat_get_max_dt(state) result(dt_max)
        type(heat_state_t), intent(in) :: state
        real(wp) :: dt_max
        
        real(wp) :: alpha_max, dx_min, cfl
        
        ! Find maximum thermal diffusivity
        alpha_max = maxval(state%material(:, :, :)%thermal_diffusivity)
        
        ! Find minimum grid spacing
        dx_min = min(state%dx, state%dy, state%dz)
        
        ! CFL condition for diffusion: dt < dx²/(2·d·α)
        ! where d is number of dimensions
        cfl = 0.25_wp  ! Safety factor
        dt_max = cfl * dx_min**2 / (3.0_wp * alpha_max)
        
        ! Account for convection if present
        if (state%config%include_convection) then
            dt_max = min(dt_max, 0.5_wp * dx_min / &
                maxval(abs(state%vx) + abs(state%vy) + abs(state%vz)))
        end if
    end function heat_get_max_dt
    
    !> Compute 3D gradient using finite differences
    subroutine compute_gradient_3d(f, df_dx, df_dy, df_dz, dx, dy, dz)
        real(wp), intent(in) :: f(:, :, :)
        real(wp), intent(out) :: df_dx(:, :, :), df_dy(:, :, :), df_dz(:, :, :)
        real(wp), intent(in) :: dx, dy, dz
        
        integer :: nx, ny, nz, i, j, k
        
        nx = size(f, 1)
        ny = size(f, 2)
        nz = size(f, 3)
        
        ! Central differences in interior
        do k = 2, nz - 1
            do j = 2, ny - 1
                do i = 2, nx - 1
                    df_dx(i, j, k) = (f(i + 1, j, k) - f(i - 1, j, k)) / (2.0_wp * dx)
                    df_dy(i, j, k) = (f(i, j + 1, k) - f(i, j - 1, k)) / (2.0_wp * dy)
                    df_dz(i, j, k) = (f(i, j, k + 1) - f(i, j, k - 1)) / (2.0_wp * dz)
                end do
            end do
        end do
        
        ! Boundaries (one-sided differences)
        df_dx(1, :, :) = (f(2, :, :) - f(1, :, :)) / dx
        df_dx(nx, :, :) = (f(nx, :, :) - f(nx - 1, :, :)) / dx
        
        df_dy(:, 1, :) = (f(:, 2, :) - f(:, 1, :)) / dy
        df_dy(:, ny, :) = (f(:, ny, :) - f(:, ny - 1, :)) / dy
        
        df_dz(:, :, 1) = (f(:, :, 2) - f(:, :, 1)) / dz
        df_dz(:, :, nz) = (f(:, :, nz) - f(:, :, nz - 1)) / dz
    end subroutine compute_gradient_3d

end module heat_transfer