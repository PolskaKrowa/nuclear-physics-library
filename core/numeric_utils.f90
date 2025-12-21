! fortran/core/numerics_utils.f90
!
! Common numerical utilities and helper functions.
! Provides array operations, interpolation, and numerical checks.
!
! Usage:
!   use numerics_utils
!   y = linspace(0.0_wp, 10.0_wp, 100)
!   z = clip(x, 0.0_wp, 1.0_wp)
!
module numerics_utils
    use kinds, only: wp, dp
    use constants, only: TOL_DEFAULT, is_zero, nearly_equal
    implicit none
    private
    
    ! Public utilities
    public :: linspace, logspace, arange
    public :: clip, sign_safe
    public :: interp1d, trapz
    public :: norm2, normalize
    public :: meshgrid_2d
    public :: safe_divide, safe_sqrt, safe_log
    public :: is_finite, is_nan
    
contains

    !> Generate linearly spaced array from start to end
    pure function linspace(start, end, n) result(array)
        real(wp), intent(in) :: start, end
        integer, intent(in) :: n
        real(wp) :: array(n)
        integer :: i
        real(wp) :: step
        
        if (n == 1) then
            array(1) = start
        else
            step = (end - start) / real(n - 1, wp)
            do i = 1, n
                array(i) = start + real(i - 1, wp) * step
            end do
            array(n) = end  ! Ensure exact end value
        end if
    end function linspace
    
    !> Generate logarithmically spaced array from 10^start to 10^end
    pure function logspace(start, end, n) result(array)
        real(wp), intent(in) :: start, end
        integer, intent(in) :: n
        real(wp) :: array(n)
        integer :: i
        real(wp) :: step
        
        if (n == 1) then
            array(1) = 10.0_wp**start
        else
            step = (end - start) / real(n - 1, wp)
            do i = 1, n
                array(i) = 10.0_wp**(start + real(i - 1, wp) * step)
            end do
        end if
    end function logspace
    
    !> Generate array from start to end with given step
    pure function arange(start, end, step) result(array)
        real(wp), intent(in) :: start, end, step
        real(wp), allocatable :: array(:)
        integer :: n, i
        
        n = max(1, ceiling((end - start) / step))
        allocate(array(n))
        
        do i = 1, n
            array(i) = start + real(i - 1, wp) * step
        end do
    end function arange
    
    !> Clip value to range [vmin, vmax]
    elemental function clip(x, vmin, vmax) result(y)
        real(wp), intent(in) :: x, vmin, vmax
        real(wp) :: y
        
        y = max(vmin, min(vmax, x))
    end function clip
    
    !> Safe sign function (returns 0 for near-zero values)
    elemental function sign_safe(x, tol) result(s)
        real(wp), intent(in) :: x
        real(wp), intent(in), optional :: tol
        real(wp) :: s
        real(wp) :: tolerance
        
        tolerance = TOL_DEFAULT
        if (present(tol)) tolerance = tol
        
        if (abs(x) < tolerance) then
            s = 0.0_wp
        else
            s = sign(1.0_wp, x)
        end if
    end function sign_safe
    
    !> Linear interpolation at point x
    pure function interp1d(xp, fp, x) result(f)
        real(wp), intent(in) :: xp(:), fp(:)
        real(wp), intent(in) :: x
        real(wp) :: f
        integer :: n, i
        real(wp) :: t
        
        n = size(xp)
        
        ! Handle boundary cases
        if (x <= xp(1)) then
            f = fp(1)
            return
        else if (x >= xp(n)) then
            f = fp(n)
            return
        end if
        
        ! Find interpolation interval
        do i = 1, n - 1
            if (x >= xp(i) .and. x <= xp(i + 1)) then
                t = (x - xp(i)) / (xp(i + 1) - xp(i))
                f = (1.0_wp - t) * fp(i) + t * fp(i + 1)
                return
            end if
        end do
        
        ! Should never reach here
        f = fp(n)
    end function interp1d
    
    !> Trapezoidal integration
    pure function trapz(y, x) result(integral)
        real(wp), intent(in) :: y(:)
        real(wp), intent(in), optional :: x(:)
        real(wp) :: integral
        integer :: n, i
        real(wp) :: dx
        
        n = size(y)
        integral = 0.0_wp
        
        if (n < 2) return
        
        if (present(x)) then
            ! Non-uniform spacing
            do i = 1, n - 1
                dx = x(i + 1) - x(i)
                integral = integral + 0.5_wp * dx * (y(i) + y(i + 1))
            end do
        else
            ! Uniform spacing (dx = 1)
            do i = 1, n - 1
                integral = integral + 0.5_wp * (y(i) + y(i + 1))
            end do
        end if
    end function trapz
    
    !> Compute L2 norm of vector
    pure function norm2(v) result(n)
        real(wp), intent(in) :: v(:)
        real(wp) :: n
        
        n = sqrt(sum(v**2))
    end function norm2
    
    !> Normalize vector to unit length
    pure function normalize(v) result(u)
        real(wp), intent(in) :: v(:)
        real(wp) :: u(size(v))
        real(wp) :: n
        
        n = norm2(v)
        if (n > tiny(1.0_wp)) then
            u = v / n
        else
            u = 0.0_wp
        end if
    end function normalize
    
    !> Create 2D meshgrid from x and y vectors
    pure subroutine meshgrid_2d(x, y, xx, yy)
        real(wp), intent(in) :: x(:), y(:)
        real(wp), intent(out) :: xx(:, :), yy(:, :)
        integer :: i, j, nx, ny
        
        nx = size(x)
        ny = size(y)
        
        do j = 1, ny
            do i = 1, nx
                xx(i, j) = x(i)
                yy(i, j) = y(j)
            end do
        end do
    end subroutine meshgrid_2d
    
    !> Safe division (returns 0 if denominator near zero)
    elemental function safe_divide(numerator, denominator, tol) result(ratio)
        real(wp), intent(in) :: numerator, denominator
        real(wp), intent(in), optional :: tol
        real(wp) :: ratio
        real(wp) :: tolerance
        
        tolerance = TOL_DEFAULT
        if (present(tol)) tolerance = tol
        
        if (abs(denominator) > tolerance) then
            ratio = numerator / denominator
        else
            ratio = 0.0_wp
        end if
    end function safe_divide
    
    !> Safe square root (returns 0 for negative values)
    elemental function safe_sqrt(x) result(y)
        real(wp), intent(in) :: x
        real(wp) :: y
        
        if (x >= 0.0_wp) then
            y = sqrt(x)
        else
            y = 0.0_wp
        end if
    end function safe_sqrt
    
    !> Safe logarithm (clamps to minimum value)
    elemental function safe_log(x, minval) result(y)
        real(wp), intent(in) :: x
        real(wp), intent(in), optional :: minval
        real(wp) :: y
        real(wp) :: xmin
        
        xmin = tiny(1.0_wp)
        if (present(minval)) xmin = minval
        
        y = log(max(x, xmin))
    end function safe_log
    
    !> Check if value is finite (not NaN or Inf)
    elemental function is_finite(x) result(finite)
        real(wp), intent(in) :: x
        logical :: finite

        ! NaN is the only value for which (x /= x) is true.
        ! Infinity will have magnitude greater than huge(1.0_wp).
        finite = .not. (x /= x) .and. (abs(x) <= huge(1.0_wp))
    end function is_finite

    !> Check if value is NaN
    elemental function is_nan(x) result(nan)
        real(wp), intent(in) :: x
        logical :: nan

        ! NaN compares unequal to itself
        nan = (x /= x)
    end function is_nan

end module numerics_utils