! fortran/kernels/pde/finite_difference.f90
!
! Finite difference methods for spatial derivatives and PDE operators.
! Provides various accuracy orders and boundary condition handling.
!
! Usage:
!   use finite_difference
!   call fd_derivative_1d(u, du_dx, dx, order=2, accuracy=2)
!   call fd_laplacian_2d(u, laplacian, dx, dy)
!
module finite_difference
    use kinds, only: wp, dp
    use constants, only: TOL_DEFAULT
    implicit none
    private
    
    ! Public interfaces
    public :: fd_derivative_1d
    public :: fd_derivative_2d
    public :: fd_laplacian_1d
    public :: fd_laplacian_2d
    public :: fd_laplacian_3d
    public :: fd_gradient_2d
    public :: fd_divergence_2d
    public :: fd_curl_2d
    
    ! Boundary condition types
    integer, parameter, public :: BC_PERIODIC = 1
    integer, parameter, public :: BC_DIRICHLET = 2
    integer, parameter, public :: BC_NEUMANN = 3
    integer, parameter, public :: BC_EXTRAPOLATE = 4
    
    ! Stencil coefficients for various accuracy orders
    real(wp), parameter :: C2_CENTERED(3) = [1.0_wp, -2.0_wp, 1.0_wp]
    real(wp), parameter :: C4_CENTERED(5) = [-1.0_wp, 16.0_wp, -30.0_wp, 16.0_wp, -1.0_wp] / 12.0_wp
    
contains

    !> Compute 1D derivative using finite differences
    !! order: 1 (first derivative), 2 (second derivative)
    !! accuracy: 2 (second order), 4 (fourth order)
    subroutine fd_derivative_1d(u, du, dx, order, accuracy, bc_type)
        real(wp), intent(in) :: u(:)
        real(wp), intent(out) :: du(:)
        real(wp), intent(in) :: dx
        integer, intent(in) :: order
        integer, intent(in), optional :: accuracy
        integer, intent(in), optional :: bc_type
        
        integer :: n, i, acc, bc
        real(wp) :: inv_dx, inv_dx2
        
        n = size(u)
        acc = 2
        bc = BC_PERIODIC
        if (present(accuracy)) acc = accuracy
        if (present(bc_type)) bc = bc_type
        
        inv_dx = 1.0_wp / dx
        inv_dx2 = 1.0_wp / (dx * dx)
        
        if (order == 1) then
            ! First derivative
            if (acc == 2) then
                call fd_first_derivative_o2(u, du, inv_dx, bc)
            else if (acc == 4) then
                call fd_first_derivative_o4(u, du, inv_dx, bc)
            else
                du = 0.0_wp  ! Unsupported accuracy
            end if
        else if (order == 2) then
            ! Second derivative
            if (acc == 2) then
                call fd_second_derivative_o2(u, du, inv_dx2, bc)
            else if (acc == 4) then
                call fd_second_derivative_o4(u, du, inv_dx2, bc)
            else
                du = 0.0_wp
            end if
        else
            du = 0.0_wp  ! Unsupported order
        end if
    end subroutine fd_derivative_1d
    
    !> First derivative, second-order accurate
    subroutine fd_first_derivative_o2(u, du, inv_dx, bc)
        real(wp), intent(in) :: u(:)
        real(wp), intent(out) :: du(:)
        real(wp), intent(in) :: inv_dx
        integer, intent(in) :: bc
        integer :: n, i
        
        n = size(u)
        
        ! Interior points (centered difference)
        do i = 2, n - 1
            du(i) = 0.5_wp * (u(i + 1) - u(i - 1)) * inv_dx
        end do
        
        ! Boundary points
        call apply_bc_first_deriv(u, du, 1, inv_dx, bc)
        call apply_bc_first_deriv(u, du, n, inv_dx, bc)
    end subroutine fd_first_derivative_o2
    
    !> First derivative, fourth-order accurate
    subroutine fd_first_derivative_o4(u, du, inv_dx, bc)
        real(wp), intent(in) :: u(:)
        real(wp), intent(out) :: du(:)
        real(wp), intent(in) :: inv_dx
        integer, intent(in) :: bc
        integer :: n, i
        
        n = size(u)
        
        ! Interior points
        do i = 3, n - 2
            du(i) = (u(i - 2) - 8.0_wp * u(i - 1) + 8.0_wp * u(i + 1) - u(i + 2)) / 12.0_wp * inv_dx
        end do
        
        ! Near-boundary points (fall back to second order)
        if (n >= 4) then
            du(2) = 0.5_wp * (u(3) - u(1)) * inv_dx
            du(n - 1) = 0.5_wp * (u(n) - u(n - 2)) * inv_dx
        end if
        
        ! Boundary points
        call apply_bc_first_deriv(u, du, 1, inv_dx, bc)
        call apply_bc_first_deriv(u, du, n, inv_dx, bc)
    end subroutine fd_first_derivative_o4
    
    !> Second derivative, second-order accurate
    subroutine fd_second_derivative_o2(u, du, inv_dx2, bc)
        real(wp), intent(in) :: u(:)
        real(wp), intent(out) :: du(:)
        real(wp), intent(in) :: inv_dx2
        integer, intent(in) :: bc
        integer :: n, i
        
        n = size(u)
        
        ! Interior points
        do i = 2, n - 1
            du(i) = (u(i - 1) - 2.0_wp * u(i) + u(i + 1)) * inv_dx2
        end do
        
        ! Boundary points
        call apply_bc_second_deriv(u, du, 1, inv_dx2, bc)
        call apply_bc_second_deriv(u, du, n, inv_dx2, bc)
    end subroutine fd_second_derivative_o2
    
    !> Second derivative, fourth-order accurate
    subroutine fd_second_derivative_o4(u, du, inv_dx2, bc)
        real(wp), intent(in) :: u(:)
        real(wp), intent(out) :: du(:)
        real(wp), intent(in) :: inv_dx2
        integer, intent(in) :: bc
        integer :: n, i
        
        n = size(u)
        
        ! Interior points
        do i = 3, n - 2
            du(i) = (-u(i - 2) + 16.0_wp * u(i - 1) - 30.0_wp * u(i) + &
                     16.0_wp * u(i + 1) - u(i + 2)) / 12.0_wp * inv_dx2
        end do
        
        ! Near-boundary (fall back to second order)
        if (n >= 4) then
            du(2) = (u(1) - 2.0_wp * u(2) + u(3)) * inv_dx2
            du(n - 1) = (u(n - 2) - 2.0_wp * u(n - 1) + u(n)) * inv_dx2
        end if
        
        ! Boundary points
        call apply_bc_second_deriv(u, du, 1, inv_dx2, bc)
        call apply_bc_second_deriv(u, du, n, inv_dx2, bc)
    end subroutine fd_second_derivative_o4
    
    !> Apply boundary conditions for first derivative
    subroutine apply_bc_first_deriv(u, du, i, inv_dx, bc)
        real(wp), intent(in) :: u(:)
        real(wp), intent(inout) :: du(:)
        integer, intent(in) :: i, bc
        real(wp), intent(in) :: inv_dx
        integer :: n
        
        n = size(u)
        
        select case(bc)
        case(BC_PERIODIC)
            if (i == 1) then
                du(1) = 0.5_wp * (u(2) - u(n)) * inv_dx
            else
                du(n) = 0.5_wp * (u(1) - u(n - 1)) * inv_dx
            end if
            
        case(BC_DIRICHLET)
            if (i == 1) then
                du(1) = (u(2) - u(1)) * inv_dx  ! Forward difference
            else
                du(n) = (u(n) - u(n - 1)) * inv_dx  ! Backward difference
            end if
            
        case(BC_NEUMANN)
            du(i) = 0.0_wp  ! Derivative specified at boundary
            
        case(BC_EXTRAPOLATE)
            if (i == 1) then
                du(1) = 2.0_wp * du(2) - du(3)
            else
                du(n) = 2.0_wp * du(n - 1) - du(n - 2)
            end if
        end select
    end subroutine apply_bc_first_deriv
    
    !> Apply boundary conditions for second derivative
    subroutine apply_bc_second_deriv(u, du, i, inv_dx2, bc)
        real(wp), intent(in) :: u(:)
        real(wp), intent(inout) :: du(:)
        integer, intent(in) :: i, bc
        real(wp), intent(in) :: inv_dx2
        integer :: n
        
        n = size(u)
        
        select case(bc)
        case(BC_PERIODIC)
            if (i == 1) then
                du(1) = (u(n) - 2.0_wp * u(1) + u(2)) * inv_dx2
            else
                du(n) = (u(n - 1) - 2.0_wp * u(n) + u(1)) * inv_dx2
            end if
            
        case(BC_DIRICHLET)
            if (i == 1) then
                du(1) = (2.0_wp * u(1) - 5.0_wp * u(2) + 4.0_wp * u(3) - u(4)) * inv_dx2
            else
                du(n) = (2.0_wp * u(n) - 5.0_wp * u(n - 1) + 4.0_wp * u(n - 2) - u(n - 3)) * inv_dx2
            end if
            
        case(BC_NEUMANN)
            ! Use one-sided differences
            if (i == 1) then
                du(1) = (u(1) - 2.0_wp * u(2) + u(3)) * inv_dx2
            else
                du(n) = (u(n) - 2.0_wp * u(n - 1) + u(n - 2)) * inv_dx2
            end if
            
        case(BC_EXTRAPOLATE)
            if (i == 1) then
                du(1) = 2.0_wp * du(2) - du(3)
            else
                du(n) = 2.0_wp * du(n - 1) - du(n - 2)
            end if
        end select
    end subroutine apply_bc_second_deriv
    
    !> Compute 2D derivative (partial with respect to x or y)
    subroutine fd_derivative_2d(u, du, dx, dy, direction, order, accuracy)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: du(:, :)
        real(wp), intent(in) :: dx, dy
        integer, intent(in) :: direction  ! 1 = x, 2 = y
        integer, intent(in) :: order
        integer, intent(in), optional :: accuracy
        
        integer :: nx, ny, i, j, acc
        real(wp), allocatable :: u_slice(:), du_slice(:)
        real(wp) :: h
        
        nx = size(u, 1)
        ny = size(u, 2)
        acc = 2
        if (present(accuracy)) acc = accuracy
        
        if (direction == 1) then
            ! Derivative with respect to x
            h = dx
            allocate(u_slice(nx), du_slice(nx))
            do j = 1, ny
                u_slice = u(:, j)
                call fd_derivative_1d(u_slice, du_slice, h, order, acc, BC_PERIODIC)
                du(:, j) = du_slice
            end do
        else
            ! Derivative with respect to y
            h = dy
            allocate(u_slice(ny), du_slice(ny))
            do i = 1, nx
                u_slice = u(i, :)
                call fd_derivative_1d(u_slice, du_slice, h, order, acc, BC_PERIODIC)
                du(i, :) = du_slice
            end do
        end if
        
        deallocate(u_slice, du_slice)
    end subroutine fd_derivative_2d
    
    !> Compute 1D Laplacian
    pure subroutine fd_laplacian_1d(u, laplacian, dx, bc_type)
        real(wp), intent(in) :: u(:)
        real(wp), intent(out) :: laplacian(:)
        real(wp), intent(in) :: dx
        integer, intent(in), optional :: bc_type
        
        integer :: n, i, bc
        real(wp) :: inv_dx2
        
        n = size(u)
        bc = BC_PERIODIC
        if (present(bc_type)) bc = bc_type
        
        inv_dx2 = 1.0_wp / (dx * dx)
        
        ! Interior points
        do i = 2, n - 1
            laplacian(i) = (u(i - 1) - 2.0_wp * u(i) + u(i + 1)) * inv_dx2
        end do
        
        ! Boundaries
        if (bc == BC_PERIODIC) then
            laplacian(1) = (u(n) - 2.0_wp * u(1) + u(2)) * inv_dx2
            laplacian(n) = (u(n - 1) - 2.0_wp * u(n) + u(1)) * inv_dx2
        else
            laplacian(1) = 0.0_wp
            laplacian(n) = 0.0_wp
        end if
    end subroutine fd_laplacian_1d
    
    !> Compute 2D Laplacian: ∇²u = ∂²u/∂x² + ∂²u/∂y²
    pure subroutine fd_laplacian_2d(u, laplacian, dx, dy)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: laplacian(:, :)
        real(wp), intent(in) :: dx, dy
        
        integer :: nx, ny, i, j
        real(wp) :: inv_dx2, inv_dy2
        
        nx = size(u, 1)
        ny = size(u, 2)
        inv_dx2 = 1.0_wp / (dx * dx)
        inv_dy2 = 1.0_wp / (dy * dy)
        
        ! Interior points (standard 5-point stencil)
        do j = 2, ny - 1
            do i = 2, nx - 1
                laplacian(i, j) = (u(i - 1, j) - 2.0_wp * u(i, j) + u(i + 1, j)) * inv_dx2 + &
                                  (u(i, j - 1) - 2.0_wp * u(i, j) + u(i, j + 1)) * inv_dy2
            end do
        end do
        
        ! Periodic boundaries
        ! Bottom and top edges
        j = 1
        do i = 2, nx - 1
            laplacian(i, j) = (u(i - 1, j) - 2.0_wp * u(i, j) + u(i + 1, j)) * inv_dx2 + &
                              (u(i, ny) - 2.0_wp * u(i, j) + u(i, j + 1)) * inv_dy2
        end do
        j = ny
        do i = 2, nx - 1
            laplacian(i, j) = (u(i - 1, j) - 2.0_wp * u(i, j) + u(i + 1, j)) * inv_dx2 + &
                              (u(i, j - 1) - 2.0_wp * u(i, j) + u(i, 1)) * inv_dy2
        end do
        
        ! Left and right edges
        i = 1
        do j = 2, ny - 1
            laplacian(i, j) = (u(nx, j) - 2.0_wp * u(i, j) + u(i + 1, j)) * inv_dx2 + &
                              (u(i, j - 1) - 2.0_wp * u(i, j) + u(i, j + 1)) * inv_dy2
        end do
        i = nx
        do j = 2, ny - 1
            laplacian(i, j) = (u(i - 1, j) - 2.0_wp * u(i, j) + u(1, j)) * inv_dx2 + &
                              (u(i, j - 1) - 2.0_wp * u(i, j) + u(i, j + 1)) * inv_dy2
        end do
        
        ! Corners
        laplacian(1, 1) = (u(nx, 1) - 2.0_wp * u(1, 1) + u(2, 1)) * inv_dx2 + &
                          (u(1, ny) - 2.0_wp * u(1, 1) + u(1, 2)) * inv_dy2
        laplacian(nx, 1) = (u(nx - 1, 1) - 2.0_wp * u(nx, 1) + u(1, 1)) * inv_dx2 + &
                           (u(nx, ny) - 2.0_wp * u(nx, 1) + u(nx, 2)) * inv_dy2
        laplacian(1, ny) = (u(nx, ny) - 2.0_wp * u(1, ny) + u(2, ny)) * inv_dx2 + &
                           (u(1, ny - 1) - 2.0_wp * u(1, ny) + u(1, 1)) * inv_dy2
        laplacian(nx, ny) = (u(nx - 1, ny) - 2.0_wp * u(nx, ny) + u(1, ny)) * inv_dx2 + &
                            (u(nx, ny - 1) - 2.0_wp * u(nx, ny) + u(nx, 1)) * inv_dy2
    end subroutine fd_laplacian_2d
    
    !> Compute 3D Laplacian
    pure subroutine fd_laplacian_3d(u, laplacian, dx, dy, dz)
        real(wp), intent(in) :: u(:, :, :)
        real(wp), intent(out) :: laplacian(:, :, :)
        real(wp), intent(in) :: dx, dy, dz
        
        integer :: nx, ny, nz, i, j, k
        real(wp) :: inv_dx2, inv_dy2, inv_dz2
        
        nx = size(u, 1)
        ny = size(u, 2)
        nz = size(u, 3)
        inv_dx2 = 1.0_wp / (dx * dx)
        inv_dy2 = 1.0_wp / (dy * dy)
        inv_dz2 = 1.0_wp / (dz * dz)
        
        ! Interior points (7-point stencil)
        do k = 2, nz - 1
            do j = 2, ny - 1
                do i = 2, nx - 1
                    laplacian(i, j, k) = &
                        (u(i - 1, j, k) - 2.0_wp * u(i, j, k) + u(i + 1, j, k)) * inv_dx2 + &
                        (u(i, j - 1, k) - 2.0_wp * u(i, j, k) + u(i, j + 1, k)) * inv_dy2 + &
                        (u(i, j, k - 1) - 2.0_wp * u(i, j, k) + u(i, j, k + 1)) * inv_dz2
                end do
            end do
        end do
        
        ! Boundaries set to zero (simplified)
        laplacian(1, :, :) = 0.0_wp
        laplacian(nx, :, :) = 0.0_wp
        laplacian(:, 1, :) = 0.0_wp
        laplacian(:, ny, :) = 0.0_wp
        laplacian(:, :, 1) = 0.0_wp
        laplacian(:, :, nz) = 0.0_wp
    end subroutine fd_laplacian_3d
    
    !> Compute 2D gradient: (∂u/∂x, ∂u/∂y)
    subroutine fd_gradient_2d(u, grad_x, grad_y, dx, dy)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: grad_x(:, :), grad_y(:, :)
        real(wp), intent(in) :: dx, dy
        
        call fd_derivative_2d(u, grad_x, dx, dy, direction=1, order=1, accuracy=2)
        call fd_derivative_2d(u, grad_y, dx, dy, direction=2, order=1, accuracy=2)
    end subroutine fd_gradient_2d
    
    !> Compute 2D divergence: ∇·F = ∂Fx/∂x + ∂Fy/∂y
    subroutine fd_divergence_2d(fx, fy, div, dx, dy)
        real(wp), intent(in) :: fx(:, :), fy(:, :)
        real(wp), intent(out) :: div(:, :)
        real(wp), intent(in) :: dx, dy
        
        real(wp), allocatable :: dfx_dx(:, :), dfy_dy(:, :)
        integer :: nx, ny
        
        nx = size(fx, 1)
        ny = size(fx, 2)
        allocate(dfx_dx(nx, ny), dfy_dy(nx, ny))
        
        call fd_derivative_2d(fx, dfx_dx, dx, dy, direction=1, order=1, accuracy=2)
        call fd_derivative_2d(fy, dfy_dy, dx, dy, direction=2, order=1, accuracy=2)
        
        div = dfx_dx + dfy_dy
        
        deallocate(dfx_dx, dfy_dy)
    end subroutine fd_divergence_2d
    
    !> Compute 2D curl: ∇×F = ∂Fy/∂x - ∂Fx/∂y
    subroutine fd_curl_2d(fx, fy, curl, dx, dy)
        real(wp), intent(in) :: fx(:, :), fy(:, :)
        real(wp), intent(out) :: curl(:, :)
        real(wp), intent(in) :: dx, dy
        
        real(wp), allocatable :: dfy_dx(:, :), dfx_dy(:, :)
        integer :: nx, ny
        
        nx = size(fx, 1)
        ny = size(fx, 2)
        allocate(dfy_dx(nx, ny), dfx_dy(nx, ny))
        
        call fd_derivative_2d(fy, dfy_dx, dx, dy, direction=1, order=1, accuracy=2)
        call fd_derivative_2d(fx, dfx_dy, dx, dy, direction=2, order=1, accuracy=2)
        
        curl = dfy_dx - dfx_dy
        
        deallocate(dfy_dx, dfx_dy)
    end subroutine fd_curl_2d

end module finite_difference