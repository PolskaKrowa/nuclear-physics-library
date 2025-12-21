! fortran/kernels/pde/finite_volume.f90
!
! Finite volume methods for conservation laws and advection-diffusion equations.
! Provides flux computation, reconstruction, and time integration schemes.
!
! Usage:
!   use finite_volume
!   call fv_advection_1d(u, flux, velocity, dx, dt, scheme=FV_UPWIND)
!   call fv_compute_fluxes_2d(u, fx, fy, velocity_x, velocity_y, scheme=FV_LAX_WENDROFF)
!
module finite_volume
    use kinds, only: wp, dp
    use constants, only: TOL_DEFAULT
    implicit none
    private
    
    ! Flux scheme types
    integer, parameter, public :: FV_UPWIND = 1
    integer, parameter, public :: FV_LAX_FRIEDRICHS = 2
    integer, parameter, public :: FV_LAX_WENDROFF = 3
    integer, parameter, public :: FV_GODUNOV = 4
    integer, parameter, public :: FV_MUSCL = 5
    
    ! Limiter types for high-order schemes
    integer, parameter, public :: LIMITER_NONE = 0
    integer, parameter, public :: LIMITER_MINMOD = 1
    integer, parameter, public :: LIMITER_VAN_LEER = 2
    integer, parameter, public :: LIMITER_SUPERBEE = 3
    
    ! Public interfaces
    public :: fv_advection_1d
    public :: fv_advection_2d
    public :: fv_diffusion_1d
    public :: fv_diffusion_2d
    public :: fv_advection_diffusion_1d
    public :: fv_compute_fluxes_1d
    public :: fv_compute_fluxes_2d
    public :: fv_reconstruct
    public :: fv_apply_limiter
    
contains

    !> 1D advection using finite volume method
    !! Solves ∂u/∂t + a·∂u/∂x = 0
    subroutine fv_advection_1d(u, u_new, velocity, dx, dt, scheme)
        real(wp), intent(in) :: u(:)
        real(wp), intent(out) :: u_new(:)
        real(wp), intent(in) :: velocity
        real(wp), intent(in) :: dx, dt
        integer, intent(in), optional :: scheme
        
        integer :: n, i, flux_scheme
        real(wp), allocatable :: flux(:)
        real(wp) :: cfl
        
        n = size(u)
        flux_scheme = FV_UPWIND
        if (present(scheme)) flux_scheme = scheme
        
        ! CFL condition check
        cfl = abs(velocity) * dt / dx
        if (cfl > 1.0_wp) then
            ! Warning: CFL condition violated
            u_new = u
            return
        end if
        
        allocate(flux(n + 1))
        
        ! Compute fluxes at cell interfaces
        call fv_compute_fluxes_1d(u, flux, velocity, flux_scheme)
        
        ! Update cell averages (conservative form)
        do i = 1, n
            u_new(i) = u(i) - (dt / dx) * (flux(i + 1) - flux(i))
        end do
        
        deallocate(flux)
    end subroutine fv_advection_1d
    
    !> 2D advection using finite volume method
    subroutine fv_advection_2d(u, u_new, vx, vy, dx, dy, dt, scheme)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: u_new(:, :)
        real(wp), intent(in) :: vx, vy
        real(wp), intent(in) :: dx, dy, dt
        integer, intent(in), optional :: scheme
        
        integer :: nx, ny, i, j, flux_scheme
        real(wp), allocatable :: flux_x(:, :), flux_y(:, :)
        real(wp) :: cfl_x, cfl_y
        
        nx = size(u, 1)
        ny = size(u, 2)
        flux_scheme = FV_UPWIND
        if (present(scheme)) flux_scheme = scheme
        
        ! CFL check
        cfl_x = abs(vx) * dt / dx
        cfl_y = abs(vy) * dt / dy
        if (cfl_x + cfl_y > 1.0_wp) then
            u_new = u
            return
        end if
        
        allocate(flux_x(nx + 1, ny), flux_y(nx, ny + 1))
        
        ! Compute fluxes
        call fv_compute_fluxes_2d(u, flux_x, flux_y, vx, vy, flux_scheme)
        
        ! Update
        do j = 1, ny
            do i = 1, nx
                u_new(i, j) = u(i, j) - (dt / dx) * (flux_x(i + 1, j) - flux_x(i, j)) &
                                       - (dt / dy) * (flux_y(i, j + 1) - flux_y(i, j))
            end do
        end do
        
        deallocate(flux_x, flux_y)
    end subroutine fv_advection_2d
    
    !> 1D diffusion using finite volume method
    !! Solves ∂u/∂t = D·∂²u/∂x²
    subroutine fv_diffusion_1d(u, u_new, diffusivity, dx, dt)
        real(wp), intent(in) :: u(:)
        real(wp), intent(out) :: u_new(:)
        real(wp), intent(in) :: diffusivity
        real(wp), intent(in) :: dx, dt
        
        integer :: n, i
        real(wp) :: alpha, inv_dx
        
        n = size(u)
        alpha = diffusivity * dt / (dx * dx)
        
        ! Stability check (alpha <= 0.5 for explicit)
        if (alpha > 0.5_wp) then
            u_new = u
            return
        end if
        
        inv_dx = 1.0_wp / dx
        
        ! Interior points
        do i = 2, n - 1
            u_new(i) = u(i) + alpha * (u(i + 1) - 2.0_wp * u(i) + u(i - 1))
        end do
        
        ! Periodic boundaries
        u_new(1) = u(1) + alpha * (u(2) - 2.0_wp * u(1) + u(n))
        u_new(n) = u(n) + alpha * (u(1) - 2.0_wp * u(n) + u(n - 1))
    end subroutine fv_diffusion_1d
    
    !> 2D diffusion using finite volume method
    subroutine fv_diffusion_2d(u, u_new, diffusivity, dx, dy, dt)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: u_new(:, :)
        real(wp), intent(in) :: diffusivity
        real(wp), intent(in) :: dx, dy, dt
        
        integer :: nx, ny, i, j
        real(wp) :: alpha_x, alpha_y
        
        nx = size(u, 1)
        ny = size(u, 2)
        alpha_x = diffusivity * dt / (dx * dx)
        alpha_y = diffusivity * dt / (dy * dy)
        
        ! Stability check
        if (alpha_x + alpha_y > 0.5_wp) then
            u_new = u
            return
        end if
        
        ! Interior points
        do j = 2, ny - 1
            do i = 2, nx - 1
                u_new(i, j) = u(i, j) + alpha_x * (u(i + 1, j) - 2.0_wp * u(i, j) + u(i - 1, j)) &
                                       + alpha_y * (u(i, j + 1) - 2.0_wp * u(i, j) + u(i, j - 1))
            end do
        end do
        
        ! Simplified periodic boundaries
        u_new(1, :) = u(1, :)
        u_new(nx, :) = u(nx, :)
        u_new(:, 1) = u(:, 1)
        u_new(:, ny) = u(:, ny)
    end subroutine fv_diffusion_2d
    
    !> Combined advection-diffusion in 1D
    subroutine fv_advection_diffusion_1d(u, u_new, velocity, diffusivity, dx, dt, scheme)
        real(wp), intent(in) :: u(:)
        real(wp), intent(out) :: u_new(:)
        real(wp), intent(in) :: velocity, diffusivity
        real(wp), intent(in) :: dx, dt
        integer, intent(in), optional :: scheme
        
        real(wp), allocatable :: u_adv(:), u_diff(:)
        integer :: n
        
        n = size(u)
        allocate(u_adv(n), u_diff(n))
        
        ! Operator splitting: advection then diffusion
        call fv_advection_1d(u, u_adv, velocity, dx, dt, scheme)
        call fv_diffusion_1d(u_adv, u_new, diffusivity, dx, dt)
        
        deallocate(u_adv, u_diff)
    end subroutine fv_advection_diffusion_1d
    
    !> Compute numerical fluxes in 1D
    subroutine fv_compute_fluxes_1d(u, flux, velocity, scheme)
        real(wp), intent(in) :: u(:)
        real(wp), intent(out) :: flux(:)
        real(wp), intent(in) :: velocity
        integer, intent(in) :: scheme
        
        integer :: n, i
        real(wp) :: u_left, u_right
        
        n = size(u)
        
        select case(scheme)
        case(FV_UPWIND)
            ! First-order upwind
            do i = 1, n + 1
                if (velocity >= 0.0_wp) then
                    if (i == 1) then
                        flux(i) = velocity * u(n)  ! Periodic
                    else
                        flux(i) = velocity * u(i - 1)
                    end if
                else
                    if (i == n + 1) then
                        flux(i) = velocity * u(1)  ! Periodic
                    else
                        flux(i) = velocity * u(i)
                    end if
                end if
            end do
            
        case(FV_LAX_FRIEDRICHS)
            ! Lax-Friedrichs flux
            do i = 2, n
                u_left = u(i - 1)
                u_right = u(i)
                flux(i) = 0.5_wp * velocity * (u_left + u_right) - &
                          0.5_wp * abs(velocity) * (u_right - u_left)
            end do
            ! Periodic boundaries
            flux(1) = 0.5_wp * velocity * (u(n) + u(1)) - &
                      0.5_wp * abs(velocity) * (u(1) - u(n))
            flux(n + 1) = flux(1)
            
        case(FV_LAX_WENDROFF)
            ! Lax-Wendroff (second-order)
            do i = 2, n
                flux(i) = velocity * (u(i - 1) + u(i)) / 2.0_wp - &
                          velocity**2 * (u(i) - u(i - 1)) / 2.0_wp
            end do
            flux(1) = velocity * (u(n) + u(1)) / 2.0_wp - &
                      velocity**2 * (u(1) - u(n)) / 2.0_wp
            flux(n + 1) = flux(1)
            
        case default
            ! Default to upwind
            flux = 0.0_wp
        end select
    end subroutine fv_compute_fluxes_1d
    
    !> Compute numerical fluxes in 2D
    subroutine fv_compute_fluxes_2d(u, flux_x, flux_y, vx, vy, scheme)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: flux_x(:, :), flux_y(:, :)
        real(wp), intent(in) :: vx, vy
        integer, intent(in) :: scheme
        
        integer :: nx, ny, i, j
        real(wp) :: u_left, u_right, u_down, u_up
        
        nx = size(u, 1)
        ny = size(u, 2)
        
        ! X-direction fluxes
        select case(scheme)
        case(FV_UPWIND)
            do j = 1, ny
                do i = 1, nx + 1
                    if (vx >= 0.0_wp) then
                        if (i == 1) then
                            flux_x(i, j) = vx * u(nx, j)
                        else
                            flux_x(i, j) = vx * u(i - 1, j)
                        end if
                    else
                        if (i == nx + 1) then
                            flux_x(i, j) = vx * u(1, j)
                        else
                            flux_x(i, j) = vx * u(i, j)
                        end if
                    end if
                end do
            end do
            
        case(FV_LAX_FRIEDRICHS)
            do j = 1, ny
                do i = 2, nx
                    u_left = u(i - 1, j)
                    u_right = u(i, j)
                    flux_x(i, j) = 0.5_wp * vx * (u_left + u_right) - &
                                   0.5_wp * abs(vx) * (u_right - u_left)
                end do
                ! Periodic
                flux_x(1, j) = 0.5_wp * vx * (u(nx, j) + u(1, j)) - &
                               0.5_wp * abs(vx) * (u(1, j) - u(nx, j))
                flux_x(nx + 1, j) = flux_x(1, j)
            end do
            
        case default
            flux_x = 0.0_wp
        end select
        
        ! Y-direction fluxes
        select case(scheme)
        case(FV_UPWIND)
            do i = 1, nx
                do j = 1, ny + 1
                    if (vy >= 0.0_wp) then
                        if (j == 1) then
                            flux_y(i, j) = vy * u(i, ny)
                        else
                            flux_y(i, j) = vy * u(i, j - 1)
                        end if
                    else
                        if (j == ny + 1) then
                            flux_y(i, j) = vy * u(i, 1)
                        else
                            flux_y(i, j) = vy * u(i, j)
                        end if
                    end if
                end do
            end do
            
        case(FV_LAX_FRIEDRICHS)
            do i = 1, nx
                do j = 2, ny
                    u_down = u(i, j - 1)
                    u_up = u(i, j)
                    flux_y(i, j) = 0.5_wp * vy * (u_down + u_up) - &
                                   0.5_wp * abs(vy) * (u_up - u_down)
                end do
                ! Periodic
                flux_y(i, 1) = 0.5_wp * vy * (u(i, ny) + u(i, 1)) - &
                               0.5_wp * abs(vy) * (u(i, 1) - u(i, ny))
                flux_y(i, ny + 1) = flux_y(i, 1)
            end do
            
        case default
            flux_y = 0.0_wp
        end select
    end subroutine fv_compute_fluxes_2d
    
    !> MUSCL reconstruction for higher-order accuracy
    subroutine fv_reconstruct(u, u_left, u_right, limiter_type)
        real(wp), intent(in) :: u(:)
        real(wp), intent(out) :: u_left(:), u_right(:)
        integer, intent(in), optional :: limiter_type
        
        integer :: n, i, lim_type
        real(wp) :: slope, r, phi
        
        n = size(u)
        lim_type = LIMITER_MINMOD
        if (present(limiter_type)) lim_type = limiter_type
        
        do i = 2, n - 1
            ! Compute slope ratio
            if (abs(u(i) - u(i - 1)) > TOL_DEFAULT) then
                r = (u(i + 1) - u(i)) / (u(i) - u(i - 1))
            else
                r = 0.0_wp
            end if
            
            ! Apply limiter
            phi = fv_apply_limiter(r, lim_type)
            
            ! Reconstruct values at interfaces
            slope = phi * (u(i) - u(i - 1)) / 2.0_wp
            u_left(i) = u(i) - slope
            u_right(i) = u(i) + slope
        end do
        
        ! Boundary values
        u_left(1) = u(1)
        u_right(1) = u(1)
        u_left(n) = u(n)
        u_right(n) = u(n)
    end subroutine fv_reconstruct
    
    !> Apply flux limiter function
    pure function fv_apply_limiter(r, limiter_type) result(phi)
        real(wp), intent(in) :: r
        integer, intent(in) :: limiter_type
        real(wp) :: phi
        
        select case(limiter_type)
        case(LIMITER_NONE)
            phi = 1.0_wp
            
        case(LIMITER_MINMOD)
            phi = max(0.0_wp, min(1.0_wp, r))
            
        case(LIMITER_VAN_LEER)
            phi = (r + abs(r)) / (1.0_wp + abs(r))
            
        case(LIMITER_SUPERBEE)
            phi = max(0.0_wp, min(1.0_wp, 2.0_wp * r), min(2.0_wp, r))
            
        case default
            phi = 1.0_wp
        end select
    end function fv_apply_limiter

end module finite_volume