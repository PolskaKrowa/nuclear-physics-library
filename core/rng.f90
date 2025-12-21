! fortran/core/rng.f90
!
! Random number generation with deterministic seeding.
! Provides uniform, normal, and other distributions.
!
! CRITICAL: This module maintains internal state for reproducibility.
! Always seed explicitly with rng_seed() for deterministic results.
!
! Usage:
!   use rng
!   call rng_seed(12345_i64)  ! Explicit seed
!   x = rng_uniform()
!   y = rng_normal(mean=0.0_wp, sigma=1.0_wp)
!
module rng
    use kinds, only: wp, i64
    use constants, only: PI, SQRT2, TWO_PI
    implicit none
    private
    
    ! Public interface
    public :: rng_seed, rng_state_save, rng_state_restore
    public :: rng_uniform, rng_uniform_array
    public :: rng_normal, rng_normal_array
    public :: rng_exponential, rng_integer
    
    ! RNG state (for reproducibility)
    type, public :: rng_state_t
        integer :: seed_size
        integer, allocatable :: seed(:)
    end type rng_state_t
    
    ! Internal state flag
    logical, save :: is_seeded = .false.
    
contains

    !> Seed the random number generator with explicit value
    !! This MUST be called before any RNG usage for deterministic results
    subroutine rng_seed(seed_value)
        integer(i64), intent(in) :: seed_value
        integer :: seed_size, i
        integer, allocatable :: seed_array(:)
        
        ! Get size of seed array
        call random_seed(size=seed_size)
        allocate(seed_array(seed_size))
        
        ! Fill seed array from input value
        ! Use different values for each element to avoid patterns
        do i = 1, seed_size
            seed_array(i) = int(ieor(seed_value, int(i * 2654435761_i64, i64)))
        end do
        
        ! Set the seed
        call random_seed(put=seed_array)
        deallocate(seed_array)
        
        is_seeded = .true.
    end subroutine rng_seed
    
    !> Save current RNG state for checkpointing
    subroutine rng_state_save(state)
        type(rng_state_t), intent(out) :: state
        
        call random_seed(size=state%seed_size)
        allocate(state%seed(state%seed_size))
        call random_seed(get=state%seed)
    end subroutine rng_state_save
    
    !> Restore RNG state from checkpoint
    subroutine rng_state_restore(state)
        type(rng_state_t), intent(in) :: state
        
        call random_seed(put=state%seed)
        is_seeded = .true.
    end subroutine rng_state_restore
    
    !> Generate uniform random number in [0, 1)
    function rng_uniform() result(r)
        real(wp) :: r
        
        if (.not. is_seeded) then
            ! Emergency fallback - seed from system time
            call rng_seed(123456789_i64)
        end if
        
        call random_number(r)
    end function rng_uniform
    
    !> Generate array of uniform random numbers
    subroutine rng_uniform_array(array)
        real(wp), intent(out) :: array(:)
        
        if (.not. is_seeded) then
            call rng_seed(123456789_i64)
        end if
        
        call random_number(array)
    end subroutine rng_uniform_array
    
    !> Generate normal (Gaussian) random number using Box-Muller transform
    function rng_normal(mean, sigma) result(r)
        real(wp), intent(in), optional :: mean, sigma
        real(wp) :: r
        real(wp) :: u1, u2, z
        real(wp) :: mu, sig
        
        ! Default parameters
        mu = 0.0_wp
        sig = 1.0_wp
        if (present(mean)) mu = mean
        if (present(sigma)) sig = sigma
        
        ! Box-Muller transform
        u1 = rng_uniform()
        u2 = rng_uniform()
        
        ! Avoid log(0)
        u1 = max(u1, tiny(1.0_wp))
        
        z = sqrt(-2.0_wp * log(u1)) * cos(TWO_PI * u2)
        r = mu + sig * z
    end function rng_normal
    
    !> Generate array of normal random numbers
    subroutine rng_normal_array(array, mean, sigma)
        real(wp), intent(out) :: array(:)
        real(wp), intent(in), optional :: mean, sigma
        integer :: i
        
        do i = 1, size(array)
            array(i) = rng_normal(mean, sigma)
        end do
    end subroutine rng_normal_array
    
    !> Generate exponential random number with given rate parameter
    function rng_exponential(lambda) result(r)
        real(wp), intent(in) :: lambda
        real(wp) :: r, u
        
        u = rng_uniform()
        u = max(u, tiny(1.0_wp))  ! Avoid log(0)
        r = -log(u) / lambda
    end function rng_exponential
    
    !> Generate random integer in range [imin, imax]
    function rng_integer(imin, imax) result(r)
        integer, intent(in) :: imin, imax
        integer :: r
        real(wp) :: u
        
        u = rng_uniform()
        r = imin + floor(u * real(imax - imin + 1, wp))
        r = min(r, imax)  ! Handle edge case where u is very close to 1
    end function rng_integer

end module rng