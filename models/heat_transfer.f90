! models/heat_transfer.f90
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
    use solve_linear, only: solve_tridiagonal, SOLVE_SUCCESS
    implicit none
    private
    
    ! Public interface
    public :: heat_config_t, heat_state_t, heat_material_t
    public :: heat_init, heat_destroy
    public :: heat_set_properties, heat_set_source, heat_set_velocity
    public :: heat_set_bc, heat_apply_bc, heat_set_coolant
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

        ! Per-cell coolant coupling (reorg step 11). When allocated,
        ! heat_step / heat_step_implicit use these instead of the legacy
        ! T_FLUID = 558 K, H_CONV = 30000 W/m²·K constants. They are
        ! populated by the orchestrator from two_phase_state_t via
        ! fuel_compute_convection_coefficients so that the convective
        ! sink reflects local flow regime (single-phase forced convection
        ! vs nucleate boiling vs DNB).
        real(wp), allocatable :: T_fluid(:, :, :)   ! [K]
        real(wp), allocatable :: h_conv(:, :, :)    ! [W/m²·K]

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
        if (allocated(state%T_fluid)) deallocate(state%T_fluid)
        if (allocated(state%h_conv)) deallocate(state%h_conv)
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
    
    !> Set per-cell coolant temperature and convective coefficient.
    !!
    !! Lazy-allocates state%T_fluid and state%h_conv on first call. Both
    !! arrays must match the kernel grid shape. When these are allocated,
    !! heat_step / heat_step_implicit use per-cell values; otherwise they
    !! fall back to the legacy T_FLUID = 558 K, H_CONV = 30000 W/m²·K
    !! constants. See fuel_compute_convection_coefficients (subsystems/
    !! fuel.f90) for the regime-aware producer.
    subroutine heat_set_coolant(state, T_fluid, h_conv)
        type(heat_state_t), intent(inout) :: state
        real(wp), intent(in) :: T_fluid(:, :, :)
        real(wp), intent(in) :: h_conv(:, :, :)

        if (.not. allocated(state%T_fluid)) &
            allocate(state%T_fluid(state%nx, state%ny, state%nz))
        if (.not. allocated(state%h_conv)) &
            allocate(state%h_conv(state%nx, state%ny, state%nz))

        state%T_fluid = T_fluid
        state%h_conv  = h_conv
    end subroutine heat_set_coolant

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
        real(wp) :: h_conv_ij, T_fluid_ij, q_conv  ! per-cell convective cooling
        real(wp) :: A_over_V                       ! Surface area / volume [1/m]

        real(wp), parameter :: D_EQUIV     = 0.01_wp      ! [m]
        real(wp), parameter :: H_CONV_DEF  = 30000.0_wp   ! [W/m²·K]
        real(wp), parameter :: T_FLUID_DEF = 558.0_wp     ! BWR sat @ 7 MPa [K]
        logical :: use_per_cell

        allocate(laplacian(state%nx, state%ny, state%nz))

        state%T_old = state%T

        ! Compute Laplacian (diffusion term)
        call fd_laplacian_3d(state%T, laplacian, state%dx, state%dy, state%dz)

        ! Compute convection terms if needed
        if (state%config%include_convection) then
            allocate(dT_dx(state%nx, state%ny, state%nz))
            allocate(dT_dy(state%nx, state%ny, state%nz))
            allocate(dT_dz(state%nx, state%ny, state%nz))

            call compute_gradient_3d(state%T, dT_dx, dT_dy, dT_dz, &
                                    state%dx, state%dy, state%dz)
        end if

        ! Geometry parameters for convective cooling. A_over_V uses the
        ! legacy 10 mm hydraulic diameter; the heated_perimeter that
        ! drives the two-phase kernel is per-cell but the BWR fuel-bundle
        ! A/V is comparable, so this is acceptable for the current step.
        A_over_V = 4.0_wp / D_EQUIV

        ! Step-11: pick per-cell coolant if the orchestrator wired it in,
        ! else fall back to the legacy constants.
        use_per_cell = allocated(state%T_fluid) .and. allocated(state%h_conv)

        ! Update temperature field
        do k = 1, state%nz
            do j = 1, state%ny
                do i = 1, state%nx
                    alpha = state%material(i, j, k)%thermal_diffusivity
                    rho_cp = state%material(i, j, k)%density * &
                            state%material(i, j, k)%specific_heat

                    ! Diffusion
                    diffusion = alpha * laplacian(i, j, k)

                    ! Convection (fluid flow)
                    convection = 0.0_wp
                    if (state%config%include_convection) then
                        convection = -(state%vx(i, j, k) * dT_dx(i, j, k) + &
                                    state%vy(i, j, k) * dT_dy(i, j, k) + &
                                    state%vz(i, j, k) * dT_dz(i, j, k))
                    end if

                    if (use_per_cell) then
                        h_conv_ij  = state%h_conv(i, j, k)
                        T_fluid_ij = state%T_fluid(i, j, k)
                    else
                        h_conv_ij  = H_CONV_DEF
                        T_fluid_ij = T_FLUID_DEF
                    end if

                    ! Convective heat removal to coolant [W/m³]
                    q_conv = h_conv_ij * A_over_V * (state%T(i, j, k) - T_fluid_ij)

                    ! Source term (includes heat generation minus cooling)
                    source = state%Q(i, j, k) / rho_cp - q_conv / rho_cp
                    
                    ! Explicit Euler update
                    state%T(i, j, k) = state%T(i, j, k) + dt * &
                        (diffusion + convection + source)
                    
                    ! Prevent unphysical temperatures. Floor at 273.15 K so
                    ! cold standby (293 K coolant) is representable; the
                    ! previous 300 K floor silently warmed the cold preset
                    ! by ~7 K every tick.
                    state%T(i, j, k) = max(state%T(i, j, k), 273.15_wp)
                    state%T(i, j, k) = min(state%T(i, j, k), 3000.0_wp)
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
    
    !> Implicit time step (unconditionally stable for any dt).
    !!
    !! Backward-Euler diffusion + analytic implicit reaction, advection-free.
    !! Diffusion is split into three 1-D ADI sweeps (Lie splitting), each
    !! solved with a tridiagonal Thomas pass; the convective sink to T_fluid
    !! is handled point-wise in closed form before the sweeps. Variable
    !! material properties are accommodated via the harmonic mean of the
    !! cell-centred diffusivities at each face.
    !!
    !! Intended for the steady-state relaxation loop and any future call
    !! site that needs to take a large dt without CFL constraints. The
    !! solver assumes velocities are zero (true during steady-state solve);
    !! for transient runs with non-zero coolant flow, use heat_step.
    subroutine heat_step_implicit(state, dt)
        type(heat_state_t), intent(inout) :: state
        real(wp), intent(in) :: dt

        real(wp), parameter :: T_FLUID_DEF = 558.0_wp   ! BWR saturation @ 7 MPa [K]
        real(wp), parameter :: H_CONV_DEF  = 30000.0_wp ! [W/m²·K]
        real(wp), parameter :: D_EQUIV     = 0.01_wp    ! [m]
        real(wp), parameter :: A_over_V    = 4.0_wp / D_EQUIV
        real(wp), parameter :: ALPHA_FLOOR = 1.0e-20_wp

        real(wp), allocatable :: alpha(:,:,:)
        real(wp), allocatable :: a_diag(:), b_diag(:), c_diag(:), rhs(:), sol(:)
        real(wp) :: rho_cp, beta, alpha_face, coef
        real(wp) :: h_conv_ij, T_fluid_ij
        logical  :: use_per_cell
        integer  :: i, j, k, n_max
        integer(i32) :: status

        if (.not. allocated(state%T)) return

        n_max = max(state%nx, state%ny, state%nz)
        allocate(alpha(state%nx, state%ny, state%nz))
        allocate(a_diag(n_max - 1), b_diag(n_max), c_diag(n_max - 1))
        allocate(rhs(n_max), sol(n_max))

        alpha = state%material(:, :, :)%thermal_diffusivity
        state%T_old = state%T

        ! Step-11: pick per-cell coolant if the orchestrator wired it in.
        use_per_cell = allocated(state%T_fluid) .and. allocated(state%h_conv)

        ! ── Reaction + source: analytic implicit step on the point-wise ODE ─────
        !   dT/dt = Q/(ρcp) - β·(T - T_fluid),  β = h_conv·A/V/(ρcp)
        !   T* = (T + dt·(Q/(ρcp) + β·T_fluid)) / (1 + dt·β)
        do k = 1, state%nz
            do j = 1, state%ny
                do i = 1, state%nx
                    rho_cp = state%material(i,j,k)%density * &
                             state%material(i,j,k)%specific_heat
                    if (use_per_cell) then
                        h_conv_ij  = state%h_conv(i, j, k)
                        T_fluid_ij = state%T_fluid(i, j, k)
                    else
                        h_conv_ij  = H_CONV_DEF
                        T_fluid_ij = T_FLUID_DEF
                    end if
                    beta = h_conv_ij * A_over_V / max(rho_cp, ALPHA_FLOOR)
                    state%T(i,j,k) = (state%T(i,j,k) + &
                        dt * (state%Q(i,j,k) / max(rho_cp, ALPHA_FLOOR) + beta * T_fluid_ij)) &
                        / (1.0_wp + dt * beta)
                end do
            end do
        end do

        ! ── ADI sweep 1: x-direction implicit diffusion ────────────────────────
        do k = 1, state%nz
            do j = 1, state%ny
                do i = 1, state%nx - 1
                    alpha_face = 2.0_wp * alpha(i,j,k) * alpha(i+1,j,k) / &
                                 max(alpha(i,j,k) + alpha(i+1,j,k), ALPHA_FLOOR)
                    coef = -dt * alpha_face / (state%dx * state%dx)
                    a_diag(i) = coef
                    c_diag(i) = coef
                end do
                b_diag(1) = 1.0_wp - c_diag(1)
                do i = 2, state%nx - 1
                    b_diag(i) = 1.0_wp - a_diag(i-1) - c_diag(i)
                end do
                b_diag(state%nx) = 1.0_wp - a_diag(state%nx - 1)
                rhs(1:state%nx) = state%T(1:state%nx, j, k)
                call solve_tridiagonal(a_diag(1:state%nx-1), b_diag(1:state%nx), &
                                       c_diag(1:state%nx-1), rhs(1:state%nx), &
                                       sol(1:state%nx), status)
                if (status == SOLVE_SUCCESS) state%T(:, j, k) = sol(1:state%nx)
            end do
        end do

        ! ── ADI sweep 2: y-direction implicit diffusion ────────────────────────
        do k = 1, state%nz
            do i = 1, state%nx
                do j = 1, state%ny - 1
                    alpha_face = 2.0_wp * alpha(i,j,k) * alpha(i,j+1,k) / &
                                 max(alpha(i,j,k) + alpha(i,j+1,k), ALPHA_FLOOR)
                    coef = -dt * alpha_face / (state%dy * state%dy)
                    a_diag(j) = coef
                    c_diag(j) = coef
                end do
                b_diag(1) = 1.0_wp - c_diag(1)
                do j = 2, state%ny - 1
                    b_diag(j) = 1.0_wp - a_diag(j-1) - c_diag(j)
                end do
                b_diag(state%ny) = 1.0_wp - a_diag(state%ny - 1)
                rhs(1:state%ny) = state%T(i, 1:state%ny, k)
                call solve_tridiagonal(a_diag(1:state%ny-1), b_diag(1:state%ny), &
                                       c_diag(1:state%ny-1), rhs(1:state%ny), &
                                       sol(1:state%ny), status)
                if (status == SOLVE_SUCCESS) state%T(i, :, k) = sol(1:state%ny)
            end do
        end do

        ! ── ADI sweep 3: z-direction implicit diffusion ────────────────────────
        do j = 1, state%ny
            do i = 1, state%nx
                do k = 1, state%nz - 1
                    alpha_face = 2.0_wp * alpha(i,j,k) * alpha(i,j,k+1) / &
                                 max(alpha(i,j,k) + alpha(i,j,k+1), ALPHA_FLOOR)
                    coef = -dt * alpha_face / (state%dz * state%dz)
                    a_diag(k) = coef
                    c_diag(k) = coef
                end do
                b_diag(1) = 1.0_wp - c_diag(1)
                do k = 2, state%nz - 1
                    b_diag(k) = 1.0_wp - a_diag(k-1) - c_diag(k)
                end do
                b_diag(state%nz) = 1.0_wp - a_diag(state%nz - 1)
                rhs(1:state%nz) = state%T(i, j, 1:state%nz)
                call solve_tridiagonal(a_diag(1:state%nz-1), b_diag(1:state%nz), &
                                       c_diag(1:state%nz-1), rhs(1:state%nz), &
                                       sol(1:state%nz), status)
                if (status == SOLVE_SUCCESS) state%T(i, j, :) = sol(1:state%nz)
            end do
        end do

        call heat_apply_bc(state)
        state%T = max(state%T, 273.15_wp)
        state%T = min(state%T, 3000.0_wp)

        state%time  = state%time + dt
        state%steps = state%steps + 1

        deallocate(alpha, a_diag, b_diag, c_diag, rhs, sol)
    end subroutine heat_step_implicit
    
    !> Apply boundary conditions on all 6 faces.
    !! The previous implementation only handled -x and +x; the y and z
    !! faces were silently skipped (the comment said "omitted for brevity"),
    !! which corrupted PDE dynamics for any non-adiabatic y/z BC.
    subroutine heat_apply_bc(state)
        type(heat_state_t), intent(inout) :: state
        integer :: i, j, k, face, nx, ny, nz
        real(wp) :: T_bc, flux, htc, T_amb, T_surf, k_therm
        
        nx = state%nx;  ny = state%ny;  nz = state%nz
        
        ! ---- -x face (i=1) ----
        face = 1
        select case(state%config%bc_type(face))
        case(BC_DIRICHLET)
            state%T(1, :, :) = state%config%bc_value(face)
        case(BC_NEUMANN)
            flux = state%config%bc_value(face)
            do k = 1, nz; do j = 1, ny
                state%T(1, j, k) = state%T(2, j, k) - flux * state%dx / &
                    state%material(1, j, k)%thermal_conductivity
            end do; end do
        case(BC_CONVECTIVE)
            htc = state%config%htc(face);  T_amb = state%config%t_ambient(face)
            do k = 1, nz; do j = 1, ny
                k_therm = state%material(1, j, k)%thermal_conductivity
                state%T(1, j, k) = (k_therm * state%T(2, j, k) / state%dx + htc * T_amb) / &
                                   (k_therm / state%dx + htc)
            end do; end do
        case(BC_ADIABATIC)
            state%T(1, :, :) = state%T(2, :, :)
        end select
        
        ! ---- +x face (i=nx) ----
        face = 2
        select case(state%config%bc_type(face))
        case(BC_DIRICHLET)
            state%T(nx, :, :) = state%config%bc_value(face)
        case(BC_NEUMANN)
            flux = state%config%bc_value(face)
            do k = 1, nz; do j = 1, ny
                state%T(nx, j, k) = state%T(nx-1, j, k) + flux * state%dx / &
                    state%material(nx, j, k)%thermal_conductivity
            end do; end do
        case(BC_CONVECTIVE)
            htc = state%config%htc(face);  T_amb = state%config%t_ambient(face)
            do k = 1, nz; do j = 1, ny
                k_therm = state%material(nx, j, k)%thermal_conductivity
                state%T(nx, j, k) = (k_therm * state%T(nx-1, j, k) / state%dx + htc * T_amb) / &
                                    (k_therm / state%dx + htc)
            end do; end do
        case(BC_ADIABATIC)
            state%T(nx, :, :) = state%T(nx-1, :, :)
        end select
        
        ! ---- -y face (j=1) ----
        face = 3
        select case(state%config%bc_type(face))
        case(BC_DIRICHLET)
            state%T(:, 1, :) = state%config%bc_value(face)
        case(BC_NEUMANN)
            flux = state%config%bc_value(face)
            do k = 1, nz; do i = 1, nx
                state%T(i, 1, k) = state%T(i, 2, k) - flux * state%dy / &
                    state%material(i, 1, k)%thermal_conductivity
            end do; end do
        case(BC_CONVECTIVE)
            htc = state%config%htc(face);  T_amb = state%config%t_ambient(face)
            do k = 1, nz; do i = 1, nx
                k_therm = state%material(i, 1, k)%thermal_conductivity
                state%T(i, 1, k) = (k_therm * state%T(i, 2, k) / state%dy + htc * T_amb) / &
                                   (k_therm / state%dy + htc)
            end do; end do
        case(BC_ADIABATIC)
            state%T(:, 1, :) = state%T(:, 2, :)
        end select
        
        ! ---- +y face (j=ny) ----
        face = 4
        select case(state%config%bc_type(face))
        case(BC_DIRICHLET)
            state%T(:, ny, :) = state%config%bc_value(face)
        case(BC_NEUMANN)
            flux = state%config%bc_value(face)
            do k = 1, nz; do i = 1, nx
                state%T(i, ny, k) = state%T(i, ny-1, k) + flux * state%dy / &
                    state%material(i, ny, k)%thermal_conductivity
            end do; end do
        case(BC_CONVECTIVE)
            htc = state%config%htc(face);  T_amb = state%config%t_ambient(face)
            do k = 1, nz; do i = 1, nx
                k_therm = state%material(i, ny, k)%thermal_conductivity
                state%T(i, ny, k) = (k_therm * state%T(i, ny-1, k) / state%dy + htc * T_amb) / &
                                    (k_therm / state%dy + htc)
            end do; end do
        case(BC_ADIABATIC)
            state%T(:, ny, :) = state%T(:, ny-1, :)
        end select
        
        ! ---- -z face (k=1) ----
        face = 5
        select case(state%config%bc_type(face))
        case(BC_DIRICHLET)
            state%T(:, :, 1) = state%config%bc_value(face)
        case(BC_NEUMANN)
            flux = state%config%bc_value(face)
            do j = 1, ny; do i = 1, nx
                state%T(i, j, 1) = state%T(i, j, 2) - flux * state%dz / &
                    state%material(i, j, 1)%thermal_conductivity
            end do; end do
        case(BC_CONVECTIVE)
            htc = state%config%htc(face);  T_amb = state%config%t_ambient(face)
            do j = 1, ny; do i = 1, nx
                k_therm = state%material(i, j, 1)%thermal_conductivity
                state%T(i, j, 1) = (k_therm * state%T(i, j, 2) / state%dz + htc * T_amb) / &
                                   (k_therm / state%dz + htc)
            end do; end do
        case(BC_ADIABATIC)
            state%T(:, :, 1) = state%T(:, :, 2)
        end select
        
        ! ---- +z face (k=nz) ----
        face = 6
        select case(state%config%bc_type(face))
        case(BC_DIRICHLET)
            state%T(:, :, nz) = state%config%bc_value(face)
        case(BC_NEUMANN)
            flux = state%config%bc_value(face)
            do j = 1, ny; do i = 1, nx
                state%T(i, j, nz) = state%T(i, j, nz-1) + flux * state%dz / &
                    state%material(i, j, nz)%thermal_conductivity
            end do; end do
        case(BC_CONVECTIVE)
            htc = state%config%htc(face);  T_amb = state%config%t_ambient(face)
            do j = 1, ny; do i = 1, nx
                k_therm = state%material(i, j, nz)%thermal_conductivity
                state%T(i, j, nz) = (k_therm * state%T(i, j, nz-1) / state%dz + htc * T_amb) / &
                                    (k_therm / state%dz + htc)
            end do; end do
        case(BC_ADIABATIC)
            state%T(:, :, nz) = state%T(:, :, nz-1)
        end select
    end subroutine heat_apply_bc
    
    !> Get maximum stable time step.
    !! Guards against division by zero when all velocities are zero
    !! (the previous code unconditionally divided by
    !! maxval(abs(vx)+abs(vy)+abs(vz)) which is 0 in pure-conduction mode).
    function heat_get_max_dt(state) result(dt_max)
        type(heat_state_t), intent(in) :: state
        real(wp) :: dt_max
        
        real(wp) :: alpha_max, dx_min, cfl, v_max
        
        ! Find maximum thermal diffusivity (guard against zero)
        alpha_max = maxval(state%material(:, :, :)%thermal_diffusivity)
        if (alpha_max < TOL_DEFAULT) alpha_max = TOL_DEFAULT
        
        ! Find minimum grid spacing
        dx_min = min(state%dx, state%dy, state%dz)
        
        ! CFL condition for diffusion: dt < dx²/(2·d·α)
        ! where d is number of dimensions
        cfl = 0.25_wp  ! Safety factor
        dt_max = cfl * dx_min**2 / (3.0_wp * alpha_max)
        
        ! Account for convection if present — guard against zero velocity
        ! (e.g. pure-conduction standby mode where vx=vy=vz=0).
        if (state%config%include_convection) then
            v_max = maxval(abs(state%vx) + abs(state%vy) + abs(state%vz))
            if (v_max > TOL_DEFAULT) then
                dt_max = min(dt_max, 0.5_wp * dx_min / v_max)
            end if
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