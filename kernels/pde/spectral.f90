! fortran/kernels/pde/spectral.f90
!
! Spectral methods for solving PDEs.
! Provides FFT-based spatial derivatives and spectral solvers.
!
! Usage:
!   use spectral
!   call spectral_derivative_1d(u, du_dx, domain_length, derivative_order)
!   call spectral_solve_poisson_2d(rhs, solution, lx, ly)
!
module spectral
    use kinds, only: wp, dp
    use constants, only: PI, TWO_PI
    implicit none
    private
    
    ! Public interfaces
    public :: spectral_derivative_1d
    public :: spectral_derivative_2d
    public :: spectral_solve_poisson_1d
    public :: spectral_solve_poisson_2d
    public :: spectral_filter
    
    ! Spectral grid type
    type, public :: spectral_grid_1d_t
        integer :: n                    ! Number of grid points
        real(wp) :: length              ! Domain length
        real(wp) :: dx                  ! Grid spacing
        real(wp), allocatable :: x(:)   ! Physical grid points
        real(wp), allocatable :: k(:)   ! Wavenumbers
    end type spectral_grid_1d_t
    
    type, public :: spectral_grid_2d_t
        integer :: nx, ny               ! Number of grid points
        real(wp) :: lx, ly              ! Domain lengths
        real(wp) :: dx, dy              ! Grid spacings
        real(wp), allocatable :: x(:), y(:)    ! Physical coordinates
        real(wp), allocatable :: kx(:), ky(:)  ! Wavenumbers
    end type spectral_grid_2d_t
    
    ! Grid constructors
    public :: spectral_grid_1d_create, spectral_grid_2d_create
    public :: spectral_grid_1d_destroy, spectral_grid_2d_destroy
    
contains

    !> Create 1D spectral grid
    subroutine spectral_grid_1d_create(grid, n, length)
        type(spectral_grid_1d_t), intent(out) :: grid
        integer, intent(in) :: n
        real(wp), intent(in) :: length
        integer :: i
        
        grid%n = n
        grid%length = length
        grid%dx = length / real(n, wp)
        
        allocate(grid%x(n))
        allocate(grid%k(n))
        
        ! Physical grid (periodic)
        do i = 1, n
            grid%x(i) = real(i - 1, wp) * grid%dx
        end do
        
        ! Wavenumbers (FFT convention)
        do i = 1, n/2 + 1
            grid%k(i) = TWO_PI * real(i - 1, wp) / length
        end do
        do i = n/2 + 2, n
            grid%k(i) = TWO_PI * real(i - 1 - n, wp) / length
        end do
    end subroutine spectral_grid_1d_create
    
    !> Create 2D spectral grid
    subroutine spectral_grid_2d_create(grid, nx, ny, lx, ly)
        type(spectral_grid_2d_t), intent(out) :: grid
        integer, intent(in) :: nx, ny
        real(wp), intent(in) :: lx, ly
        integer :: i
        
        grid%nx = nx
        grid%ny = ny
        grid%lx = lx
        grid%ly = ly
        grid%dx = lx / real(nx, wp)
        grid%dy = ly / real(ny, wp)
        
        allocate(grid%x(nx), grid%y(ny))
        allocate(grid%kx(nx), grid%ky(ny))
        
        ! Physical grids
        do i = 1, nx
            grid%x(i) = real(i - 1, wp) * grid%dx
        end do
        do i = 1, ny
            grid%y(i) = real(i - 1, wp) * grid%dy
        end do
        
        ! Wavenumbers
        do i = 1, nx/2 + 1
            grid%kx(i) = TWO_PI * real(i - 1, wp) / lx
        end do
        do i = nx/2 + 2, nx
            grid%kx(i) = TWO_PI * real(i - 1 - nx, wp) / lx
        end do
        
        do i = 1, ny/2 + 1
            grid%ky(i) = TWO_PI * real(i - 1, wp) / ly
        end do
        do i = ny/2 + 2, ny
            grid%ky(i) = TWO_PI * real(i - 1 - ny, wp) / ly
        end do
    end subroutine spectral_grid_2d_create
    
    !> Destroy 1D grid
    subroutine spectral_grid_1d_destroy(grid)
        type(spectral_grid_1d_t), intent(inout) :: grid
        if (allocated(grid%x)) deallocate(grid%x)
        if (allocated(grid%k)) deallocate(grid%k)
    end subroutine spectral_grid_1d_destroy
    
    !> Destroy 2D grid
    subroutine spectral_grid_2d_destroy(grid)
        type(spectral_grid_2d_t), intent(inout) :: grid
        if (allocated(grid%x)) deallocate(grid%x)
        if (allocated(grid%y)) deallocate(grid%y)
        if (allocated(grid%kx)) deallocate(grid%kx)
        if (allocated(grid%ky)) deallocate(grid%ky)
    end subroutine spectral_grid_2d_destroy
    
    !> Compute spectral derivative in 1D
    !! Uses FFT to compute d^n u / dx^n
    subroutine spectral_derivative_1d(u, du, length, order)
        real(wp), intent(in) :: u(:)
        real(wp), intent(out) :: du(:)
        real(wp), intent(in) :: length
        integer, intent(in) :: order
        
        integer :: n, i
        complex(wp), allocatable :: u_hat(:), du_hat(:)
        real(wp) :: k
        
        n = size(u)
        allocate(u_hat(n), du_hat(n))
        
        ! Forward FFT (simplified - in practice use FFTW)
        call simple_fft_forward(u, u_hat)
        
        ! Apply spectral derivative
        do i = 1, n
            if (i <= n/2 + 1) then
                k = TWO_PI * real(i - 1, wp) / length
            else
                k = TWO_PI * real(i - 1 - n, wp) / length
            end if
            du_hat(i) = (cmplx(0.0_wp, 1.0_wp, wp) * k)**order * u_hat(i)
        end do
        
        ! Inverse FFT
        call simple_fft_backward(du_hat, du)
        
        deallocate(u_hat, du_hat)
    end subroutine spectral_derivative_1d
    
    !> Compute spectral derivative in 2D
    subroutine spectral_derivative_2d(u, du_dx, du_dy, lx, ly)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: du_dx(:, :), du_dy(:, :)
        real(wp), intent(in) :: lx, ly
        
        integer :: nx, ny, i, j
        complex(wp), allocatable :: u_hat(:, :)
        real(wp) :: kx, ky
        
        nx = size(u, 1)
        ny = size(u, 2)
        allocate(u_hat(nx, ny))
        
        ! 2D FFT
        call simple_fft_2d_forward(u, u_hat)
        
        ! x-derivative
        do j = 1, ny
            do i = 1, nx
                if (i <= nx/2 + 1) then
                    kx = TWO_PI * real(i - 1, wp) / lx
                else
                    kx = TWO_PI * real(i - 1 - nx, wp) / lx
                end if
                u_hat(i, j) = cmplx(0.0_wp, 1.0_wp, wp) * kx * u_hat(i, j)
            end do
        end do
        call simple_fft_2d_backward(u_hat, du_dx)
        
        ! Recompute u_hat for y-derivative
        call simple_fft_2d_forward(u, u_hat)
        
        ! y-derivative
        do j = 1, ny
            if (j <= ny/2 + 1) then
                ky = TWO_PI * real(j - 1, wp) / ly
            else
                ky = TWO_PI * real(j - 1 - ny, wp) / ly
            end if
            do i = 1, nx
                u_hat(i, j) = cmplx(0.0_wp, 1.0_wp, wp) * ky * u_hat(i, j)
            end do
        end do
        call simple_fft_2d_backward(u_hat, du_dy)
        
        deallocate(u_hat)
    end subroutine spectral_derivative_2d
    
    !> Solve Poisson equation in 1D: d²u/dx² = f
    subroutine spectral_solve_poisson_1d(f, u, length)
        real(wp), intent(in) :: f(:)
        real(wp), intent(out) :: u(:)
        real(wp), intent(in) :: length
        
        integer :: n, i
        complex(wp), allocatable :: f_hat(:), u_hat(:)
        real(wp) :: k, k2
        
        n = size(f)
        allocate(f_hat(n), u_hat(n))
        
        ! Forward FFT of RHS
        call simple_fft_forward(f, f_hat)
        
        ! Solve in Fourier space: -k² û = f̂
        do i = 1, n
            if (i <= n/2 + 1) then
                k = TWO_PI * real(i - 1, wp) / length
            else
                k = TWO_PI * real(i - 1 - n, wp) / length
            end if
            
            k2 = k * k
            if (abs(k2) > tiny(1.0_wp)) then
                u_hat(i) = -f_hat(i) / k2
            else
                u_hat(i) = cmplx(0.0_wp, 0.0_wp, wp)  ! Zero mean
            end if
        end do
        
        ! Inverse FFT
        call simple_fft_backward(u_hat, u)
        
        deallocate(f_hat, u_hat)
    end subroutine spectral_solve_poisson_1d
    
    !> Solve Poisson equation in 2D: ∇²u = f
    subroutine spectral_solve_poisson_2d(f, u, lx, ly)
        real(wp), intent(in) :: f(:, :)
        real(wp), intent(out) :: u(:, :)
        real(wp), intent(in) :: lx, ly
        
        integer :: nx, ny, i, j
        complex(wp), allocatable :: f_hat(:, :), u_hat(:, :)
        real(wp) :: kx, ky, k2
        
        nx = size(f, 1)
        ny = size(f, 2)
        allocate(f_hat(nx, ny), u_hat(nx, ny))
        
        ! 2D FFT of RHS
        call simple_fft_2d_forward(f, f_hat)
        
        ! Solve in Fourier space
        do j = 1, ny
            if (j <= ny/2 + 1) then
                ky = TWO_PI * real(j - 1, wp) / ly
            else
                ky = TWO_PI * real(j - 1 - ny, wp) / ly
            end if
            
            do i = 1, nx
                if (i <= nx/2 + 1) then
                    kx = TWO_PI * real(i - 1, wp) / lx
                else
                    kx = TWO_PI * real(i - 1 - nx, wp) / lx
                end if
                
                k2 = kx*kx + ky*ky
                if (abs(k2) > tiny(1.0_wp)) then
                    u_hat(i, j) = -f_hat(i, j) / k2
                else
                    u_hat(i, j) = cmplx(0.0_wp, 0.0_wp, wp)
                end if
            end do
        end do
        
        ! Inverse 2D FFT
        call simple_fft_2d_backward(u_hat, u)
        
        deallocate(f_hat, u_hat)
    end subroutine spectral_solve_poisson_2d
    
    !> Apply spectral filter to remove high frequencies
    subroutine spectral_filter(u, cutoff_fraction)
        real(wp), intent(inout) :: u(:)
        real(wp), intent(in) :: cutoff_fraction
        
        integer :: n, i, cutoff_mode
        complex(wp), allocatable :: u_hat(:)
        
        n = size(u)
        cutoff_mode = nint(cutoff_fraction * real(n / 2, wp))
        allocate(u_hat(n))
        
        call simple_fft_forward(u, u_hat)
        
        ! Zero out high frequencies
        do i = cutoff_mode + 1, n - cutoff_mode
            u_hat(i) = cmplx(0.0_wp, 0.0_wp, wp)
        end do
        
        call simple_fft_backward(u_hat, u)
        
        deallocate(u_hat)
    end subroutine spectral_filter
    
    !> Simple DFT forward transform (placeholder for FFTW)
    !! In production, replace with FFTW calls
    subroutine simple_fft_forward(x, x_hat)
        real(wp), intent(in) :: x(:)
        complex(wp), intent(out) :: x_hat(:)
        integer :: n, k, j
        real(wp) :: theta
        
        n = size(x)
        x_hat = cmplx(0.0_wp, 0.0_wp, wp)
        
        do k = 1, n
            do j = 1, n
                theta = -TWO_PI * real((k - 1) * (j - 1), wp) / real(n, wp)
                x_hat(k) = x_hat(k) + x(j) * cmplx(cos(theta), sin(theta), wp)
            end do
        end do
    end subroutine simple_fft_forward
    
    !> Simple DFT backward transform (placeholder for FFTW)
    subroutine simple_fft_backward(x_hat, x)
        complex(wp), intent(in) :: x_hat(:)
        real(wp), intent(out) :: x(:)
        integer :: n, k, j
        real(wp) :: theta
        complex(wp) :: sum_val
        
        n = size(x)
        
        do j = 1, n
            sum_val = cmplx(0.0_wp, 0.0_wp, wp)
            do k = 1, n
                theta = TWO_PI * real((k - 1) * (j - 1), wp) / real(n, wp)
                sum_val = sum_val + x_hat(k) * cmplx(cos(theta), sin(theta), wp)
            end do
            x(j) = real(sum_val, wp) / real(n, wp)
        end do
    end subroutine simple_fft_backward
    
    !> Simple 2D DFT forward (placeholder for FFTW)
    subroutine simple_fft_2d_forward(x, x_hat)
        real(wp), intent(in) :: x(:, :)
        complex(wp), intent(out) :: x_hat(:, :)
        integer :: nx, ny, i
        complex(wp), allocatable :: temp(:, :)
        
        nx = size(x, 1)
        ny = size(x, 2)
        allocate(temp(nx, ny))
        
        ! Transform rows
        do i = 1, ny
            call simple_fft_forward(x(:, i), temp(:, i))
        end do
        
        ! Transform columns
        do i = 1, nx
            call simple_fft_forward(real(temp(i, :), wp), x_hat(i, :))
        end do
        
        deallocate(temp)
    end subroutine simple_fft_2d_forward
    
    !> Simple 2D DFT backward (placeholder for FFTW)
    subroutine simple_fft_2d_backward(x_hat, x)
        complex(wp), intent(in) :: x_hat(:, :)
        real(wp), intent(out) :: x(:, :)
        integer :: nx, ny, i
        complex(wp), allocatable :: temp(:, :)
        real(wp), allocatable :: temp_real(:)
        
        nx = size(x, 1)
        ny = size(x, 2)
        allocate(temp(nx, ny), temp_real(ny))
        
        ! Transform columns
        do i = 1, nx
            call simple_fft_backward(x_hat(i, :), temp_real)
            temp(i, :) = cmplx(temp_real, 0.0_wp, wp)
        end do
        
        ! Transform rows
        do i = 1, ny
            call simple_fft_backward(temp(:, i), x(:, i))
        end do
        
        deallocate(temp, temp_real)
    end subroutine simple_fft_2d_backward

end module spectral