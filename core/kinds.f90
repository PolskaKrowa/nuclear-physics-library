! fortran/core/kinds.f90
!
! Defines standard kind parameters for all numerical types used in the project.
! This ensures portability and numerical consistency across different platforms.
!
! Usage:
!   use kinds
!   real(dp) :: x
!   integer(i64) :: count
!
module kinds
    use iso_fortran_env, only: int8, int16, int32, int64, real32, real64, real128
    implicit none
    private
    
    ! Integer kinds
    public :: i8, i16, i32, i64
    integer, parameter :: i8  = int8   ! 1 byte integer
    integer, parameter :: i16 = int16  ! 2 byte integer
    integer, parameter :: i32 = int32  ! 4 byte integer
    integer, parameter :: i64 = int64  ! 8 byte integer
    
    ! Real kinds
    public :: sp, dp, qp
    integer, parameter :: sp = real32   ! Single precision (32-bit)
    integer, parameter :: dp = real64   ! Double precision (64-bit)
    integer, parameter :: qp = real128  ! Quadruple precision (128-bit)
    
    ! Complex kinds
    public :: spc, dpc, qpc
    integer, parameter :: spc = sp      ! Single precision complex
    integer, parameter :: dpc = dp      ! Double precision complex
    integer, parameter :: qpc = qp      ! Quadruple precision complex
    
    ! Default working precision (double precision)
    public :: wp
    integer, parameter :: wp = dp
    
    ! String length parameters
    public :: SHORT_STR, MEDIUM_STR, LONG_STR
    integer, parameter :: SHORT_STR  = 64
    integer, parameter :: MEDIUM_STR = 256
    integer, parameter :: LONG_STR   = 1024
    
    ! Numerical properties for each kind
    public :: get_epsilon, get_huge, get_tiny, get_digits
    
contains

    !> Returns machine epsilon for given kind
    pure function get_epsilon(kind_param) result(eps)
        integer, intent(in) :: kind_param
        real(dp) :: eps
        
        select case(kind_param)
        case(sp)
            eps = real(epsilon(1.0_sp), dp)
        case(dp)
            eps = epsilon(1.0_dp)
        case(qp)
            eps = real(epsilon(1.0_qp), dp)
        case default
            eps = epsilon(1.0_dp)
        end select
    end function get_epsilon
    
    !> Returns largest representable number for given kind
    pure function get_huge(kind_param) result(big)
        integer, intent(in) :: kind_param
        real(dp) :: big
        
        select case(kind_param)
        case(sp)
            big = real(huge(1.0_sp), dp)
        case(dp)
            big = huge(1.0_dp)
        case(qp)
            big = huge(1.0_dp)
        case default
            big = huge(1.0_dp)
        end select
    end function get_huge
    
    !> Returns smallest positive normalized number for given kind
    pure function get_tiny(kind_param) result(small)
        integer, intent(in) :: kind_param
        real(dp) :: small
        
        select case(kind_param)
        case(sp)
            small = real(tiny(1.0_sp), dp)
        case(dp)
            small = tiny(1.0_dp)
        case(qp)
            small = real(tiny(1.0_qp), dp)
        case default
            small = tiny(1.0_dp)
        end select
    end function get_tiny
    
    !> Returns number of significant decimal digits for given kind
    pure function get_digits(kind_param) result(dig)
        integer, intent(in) :: kind_param
        integer :: dig
        
        select case(kind_param)
        case(sp)
            dig = digits(1.0_sp)
        case(dp)
            dig = digits(1.0_dp)
        case(qp)
            dig = digits(1.0_qp)
        case default
            dig = digits(1.0_dp)
        end select
    end function get_digits

end module kinds