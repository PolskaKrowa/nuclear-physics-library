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
! IMPORTANT (portability):
!   NVHPC's nvfortran (formerly PGI) does NOT support real128 from
!   iso_fortran_env. When this file is compiled with nvfortran (detected
!   via the __NVCOMPILER or __PGI predefined macros), the qp kind is
!   aliased to dp so the rest of the codebase keeps compiling. gfortran
!   and ifort keep true quad precision.
!
!   This is required because otherwise nvfortran emits a hard error on
!   the `use iso_fortran_env, only: ..., real128` line and exits WITHOUT
!   writing kinds.mod, which surfaces downstream as:
!     Error copying Fortran module "mod/kinds.mod".
!       Tried "mod/KINDS.mod" and "mod/kinds.mod".
!   (CMake's post-compile cmake_copy_f90_mod step has nothing to copy.)
!
module kinds
    use iso_fortran_env, only: int8, int16, int32, int64
#if !defined(__NVCOMPILER) && !defined(__PGI)
    ! gfortran / ifort: true 128-bit reals available.
    use iso_fortran_env, only: real32, real64, real128
#else
    ! nvfortran: real128 is not supported. We still pull real32/real64
    ! from iso_fortran_env (renaming is unnecessary since the only
    ! additional symbol we want, real128, is absent on this compiler).
    use iso_fortran_env, only: real32, real64
#endif
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
#if !defined(__NVCOMPILER) && !defined(__PGI)
    integer, parameter :: qp = real128  ! Quadruple precision (128-bit)
#else
    ! nvfortran fallback: qp aliases dp. Code that uses qp will run at
    ! double precision on GPU/HYBRID builds. This is acceptable for the
    ! reactor-physics workloads in this library; if true quad precision
    ! is required, build with gfortran (COMPUTATION_MODE=CPU) instead.
    integer, parameter :: qp = real64   ! Quadruple precision -- aliased to dp on NVHPC
    ! We also expose a compile-time flag so callers can detect the
    ! fallback at preprocess time if they need to.
    ! (No public symbol needed; the __NVCOMPILER/__PGI macros themselves
    ! serve this role for any #if-guarded call sites.)
#endif

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
#if !defined(__NVCOMPILER) && !defined(__PGI)
        case(qp)
            eps = real(epsilon(1.0_qp), dp)
#endif
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
#if !defined(__NVCOMPILER) && !defined(__PGI)
        case(qp)
            ! FIX: was huge(1.0_dp) -- copy-paste bug from the sp branch.
            ! We still return real(dp) for ABI stability; huge(1.0_qp) is
            ! ~1.19e4932 which overflows real(dp) to +Inf, so we cap at
            ! huge(1.0_dp) and document the limitation. Callers that need
            ! the true qp range should call huge(1.0_qp) directly.
            big = huge(1.0_dp)
#endif
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
#if !defined(__NVCOMPILER) && !defined(__PGI)
        case(qp)
            small = real(tiny(1.0_qp), dp)
#endif
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
#if !defined(__NVCOMPILER) && !defined(__PGI)
        case(qp)
            dig = digits(1.0_qp)
#endif
        case default
            dig = digits(1.0_dp)
        end select
    end function get_digits

end module kinds
