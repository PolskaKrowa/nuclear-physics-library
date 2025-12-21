! fortran/core/constants.f90
!
! Provides mathematical and physical constants at working precision.
! All constants are defined with sufficient accuracy for double precision.
!
! Usage:
!   use constants
!   real(wp) :: circumference = 2.0_wp * PI * radius
!
module constants
    use kinds, only: wp, dp
    implicit none
    private
    
    ! Mathematical constants
    public :: PI, TWO_PI, HALF_PI, QUARTER_PI
    public :: E, SQRT2, SQRT3, SQRT5
    public :: PHI, LN2, LN10
    public :: DEG_TO_RAD, RAD_TO_DEG
    
    ! Pi and related
    real(wp), parameter :: PI         = 3.141592653589793238462643383279502884197_wp
    real(wp), parameter :: TWO_PI     = 6.283185307179586476925286766559005768394_wp
    real(wp), parameter :: HALF_PI    = 1.570796326794896619231321691639751442099_wp
    real(wp), parameter :: QUARTER_PI = 0.785398163397448309615660845819875721049_wp
    
    ! Euler's number
    real(wp), parameter :: E = 2.718281828459045235360287471352662497757_wp
    
    ! Common square roots
    real(wp), parameter :: SQRT2 = 1.414213562373095048801688724209698078570_wp
    real(wp), parameter :: SQRT3 = 1.732050807568877293527446341505872366943_wp
    real(wp), parameter :: SQRT5 = 2.236067977499789696409173668731276235441_wp
    
    ! Golden ratio
    real(wp), parameter :: PHI = 1.618033988749894848204586834365638117720_wp
    
    ! Natural logarithms
    real(wp), parameter :: LN2  = 0.693147180559945309417232121458176568076_wp
    real(wp), parameter :: LN10 = 2.302585092994045684017991454684364207601_wp
    
    ! Angle conversions
    real(wp), parameter :: DEG_TO_RAD = PI / 180.0_wp
    real(wp), parameter :: RAD_TO_DEG = 180.0_wp / PI
    
    ! Physical constants (SI units)
    public :: C_LIGHT, G_GRAV, H_PLANCK, HBAR, K_BOLTZMANN
    public :: E_CHARGE, M_ELECTRON, M_PROTON, M_NEUTRON
    public :: N_AVOGADRO, R_GAS, EPS0, MU0
    
    ! Speed of light (m/s)
    real(wp), parameter :: C_LIGHT = 299792458.0_wp
    
    ! Gravitational constant (m³/kg/s²)
    real(wp), parameter :: G_GRAV = 6.67430e-11_wp
    
    ! Planck constant (J·s)
    real(wp), parameter :: H_PLANCK = 6.62607015e-34_wp
    real(wp), parameter :: HBAR = H_PLANCK / TWO_PI
    
    ! Boltzmann constant (J/K)
    real(wp), parameter :: K_BOLTZMANN = 1.380649e-23_wp
    
    ! Elementary charge (C)
    real(wp), parameter :: E_CHARGE = 1.602176634e-19_wp
    
    ! Particle masses (kg)
    real(wp), parameter :: M_ELECTRON = 9.1093837015e-31_wp
    real(wp), parameter :: M_PROTON   = 1.67262192369e-27_wp
    real(wp), parameter :: M_NEUTRON  = 1.67492749804e-27_wp
    
    ! Avogadro's number (1/mol)
    real(wp), parameter :: N_AVOGADRO = 6.02214076e23_wp
    
    ! Gas constant (J/mol/K)
    real(wp), parameter :: R_GAS = 8.314462618_wp
    
    ! Vacuum permittivity (F/m)
    real(wp), parameter :: EPS0 = 8.8541878128e-12_wp
    
    ! Vacuum permeability (H/m)
    real(wp), parameter :: MU0 = 1.25663706212e-6_wp
    
    ! Numerical tolerance defaults
    public :: TOL_SP, TOL_DP, TOL_QP, TOL_DEFAULT
    
    real(wp), parameter :: TOL_SP      = 1.0e-6_wp
    real(wp), parameter :: TOL_DP      = 1.0e-12_wp
    real(wp), parameter :: TOL_QP      = 1.0e-24_wp
    real(wp), parameter :: TOL_DEFAULT = TOL_DP
    
    ! Utility functions
    public :: nearly_equal, is_zero
    
contains

    !> Test if two floating-point numbers are nearly equal
    pure function nearly_equal(a, b, tol) result(equal)
        real(wp), intent(in) :: a, b
        real(wp), intent(in), optional :: tol
        logical :: equal
        real(wp) :: tolerance, rel_diff
        
        tolerance = TOL_DEFAULT
        if (present(tol)) tolerance = tol
        
        if (a == b) then
            equal = .true.
        else
            rel_diff = abs(a - b) / max(abs(a), abs(b), tiny(1.0_wp))
            equal = rel_diff < tolerance
        end if
    end function nearly_equal
    
    !> Test if a number is effectively zero
    pure function is_zero(x, tol) result(zero)
        real(wp), intent(in) :: x
        real(wp), intent(in), optional :: tol
        logical :: zero
        real(wp) :: tolerance
        
        tolerance = TOL_DEFAULT
        if (present(tol)) tolerance = tol
        
        zero = abs(x) < tolerance
    end function is_zero

end module constants