! models/cross_sections.f90
!
! Cross-section library with temperature and density feedback.
! Implements:
! - Doppler broadening (fuel temperature effects)
! - Moderator density effects
! - Xenon/Samarium poisoning
! - Soluble boron worth
! - Control rod worth curves
! - Burnup-dependent cross sections
!
! Usage:
!   use cross_sections
!   type(xsec_library_t) :: xslib
!   type(mg_xsec_t) :: xsec
!   
!   call xslib_init(xslib, n_groups=2)
!   call xslib_load_material(xslib, "UO2_35", enrichment=0.035_wp)
!   call xslib_get_xsec(xslib, "UO2_35", T_fuel, rho_mod, burnup, xsec)
!   call xslib_apply_doppler(xsec, T_fuel, T_ref)
!   call xslib_apply_moderator_density(xsec, rho_mod, rho_ref)
!   call xslib_destroy(xslib)
!
module cross_sections
    use kinds, only: wp, i32
    use constants, only: TOL_DEFAULT, N_AVOGADRO
    use multigroup_diffusion, only: mg_xsec_t
    implicit none
    private
    
    ! Public interface
    public :: xsec_library_t, xsec_material_t
    public :: xslib_init, xslib_destroy
    public :: xslib_add_material, xslib_get_xsec
    public :: xslib_apply_doppler, xslib_apply_moderator_density
    public :: xslib_apply_xenon, xslib_apply_boron
    public :: xslib_create_two_group_fuel, xslib_create_two_group_moderator
    public :: compute_doppler_coefficient, compute_moderator_coefficient
    public :: compute_xenon_worth, compute_samarium_worth
    
    ! Maximum materials in library
    integer, parameter :: MAX_MATERIALS = 100
    
    ! Material data at reference conditions
    type :: xsec_material_t
        character(len=64) :: name
        integer :: n_groups
        logical :: is_fuel
        
        ! Reference conditions
        real(wp) :: T_ref = 900.0_wp          ! K (fuel)
        real(wp) :: rho_ref = 0.7_wp          ! g/cm³ (moderator)
        real(wp) :: burnup_ref = 0.0_wp       ! MWd/kgU
        
        ! Base cross sections (at reference)
        type(mg_xsec_t) :: xsec_base
        
        ! Temperature coefficients
        real(wp), allocatable :: alpha_D(:)       ! Doppler coefficient [pcm/K]
        real(wp), allocatable :: alpha_mod(:)     ! Moderator temp coeff [pcm/K]
        
        ! Density coefficients
        real(wp), allocatable :: alpha_rho(:)     ! Density coefficient [pcm/(g/cm³)]
        real(wp), allocatable :: alpha_void(:)    ! Void coefficient [pcm/%void]
        
        ! Burnup data (polynomial fit)
        real(wp), allocatable :: burnup_coeff(:, :)  ! (n_groups, order)
        
        ! Isotopic data (for detailed tracking)
        real(wp) :: N_U235 = 0.0_wp           ! atoms/barn-cm
        real(wp) :: N_U238 = 0.0_wp
        real(wp) :: N_Pu239 = 0.0_wp
        real(wp) :: N_Pu241 = 0.0_wp
        real(wp) :: enrichment = 0.035_wp
    end type xsec_material_t
    
    ! Cross-section library
    type :: xsec_library_t
        integer :: n_materials = 0
        integer :: n_groups = 2
        type(xsec_material_t) :: materials(MAX_MATERIALS)
        
        ! Global constants
        real(wp) :: boron_worth = -8.0_wp     ! pcm/ppm
        real(wp) :: xenon_yield = 0.061_wp    ! Fission yield
        real(wp) :: sm_yield = 0.011_wp       ! Fission yield
    end type xsec_library_t
    
contains

    !> Initialize cross-section library
    subroutine xslib_init(library, n_groups)
        type(xsec_library_t), intent(out) :: library
        integer, intent(in), optional :: n_groups
        
        library%n_materials = 0
        
        if (present(n_groups)) then
            library%n_groups = n_groups
        end if
        
        ! Add default materials
        call add_default_materials(library)
    end subroutine xslib_init
    
    !> Destroy library
    subroutine xslib_destroy(library)
        type(xsec_library_t), intent(inout) :: library
        
        library%n_materials = 0
    end subroutine xslib_destroy
    
    !> Add default materials (UO2, water, etc.)
    subroutine add_default_materials(library)
        type(xsec_library_t), intent(inout) :: library
        
        type(xsec_material_t) :: mat
        
        ! UO2 fuel at 3.5% enrichment (2-group)
        if (library%n_groups == 2) then
            call xslib_create_two_group_fuel(mat, enrichment=0.035_wp)
            mat%name = "UO2_35"
            call xslib_add_material(library, mat)
            
            ! Water moderator
            call xslib_create_two_group_moderator(mat)
            mat%name = "H2O"
            call xslib_add_material(library, mat)
        end if
    end subroutine add_default_materials
    
    !> Add material to library
    subroutine xslib_add_material(library, material)
        type(xsec_library_t), intent(inout) :: library
        type(xsec_material_t), intent(in) :: material
        
        if (library%n_materials < MAX_MATERIALS) then
            library%n_materials = library%n_materials + 1
            library%materials(library%n_materials) = material
        else
            print *, "Warning: Cross-section library full"
        end if
    end subroutine xslib_add_material
    
    !> Get cross sections with feedback applied
    subroutine xslib_get_xsec(library, material_name, T_fuel, rho_mod, burnup, &
                              xsec, Xe_conc, Sm_conc, boron_ppm)
        type(xsec_library_t), intent(in) :: library
        character(len=*), intent(in) :: material_name
        real(wp), intent(in) :: T_fuel
        real(wp), intent(in) :: rho_mod
        real(wp), intent(in) :: burnup
        type(mg_xsec_t), intent(out) :: xsec
        real(wp), intent(in), optional :: Xe_conc, Sm_conc, boron_ppm
        
        type(xsec_material_t) :: mat
        integer :: i
        logical :: found
        
        ! Find material
        found = .false.
        do i = 1, library%n_materials
            if (trim(library%materials(i)%name) == trim(material_name)) then
                mat = library%materials(i)
                found = .true.
                exit
            end if
        end do
        
        if (.not. found) then
            print *, "Error: Material not found: ", trim(material_name)
            return
        end if
        
        ! Start with base cross sections
        xsec = mat%xsec_base
        
        ! Apply temperature feedback
        if (mat%is_fuel) then
            call xslib_apply_doppler(xsec, T_fuel, mat%T_ref, mat%alpha_D)
        else
            call xslib_apply_moderator_temperature(xsec, T_fuel, mat%T_ref, mat%alpha_mod)
        end if
        
        ! Apply density feedback
        call xslib_apply_moderator_density(xsec, rho_mod, mat%rho_ref, mat%alpha_rho)
        
        ! Apply burnup
        if (burnup > TOL_DEFAULT) then
            call apply_burnup_effects(xsec, burnup, mat)
        end if
        
        ! Apply xenon
        if (present(Xe_conc)) then
            call xslib_apply_xenon(xsec, Xe_conc)
        end if
        
        ! Apply samarium
        if (present(Sm_conc)) then
            call xslib_apply_samarium(xsec, Sm_conc)
        end if
        
        ! Apply boron
        if (present(boron_ppm)) then
            call xslib_apply_boron(xsec, boron_ppm, library%boron_worth)
        end if
    end subroutine xslib_get_xsec
    
    !> Apply Doppler broadening (fuel temperature effect)
    subroutine xslib_apply_doppler(xsec, T_fuel, T_ref, alpha_D)
        type(mg_xsec_t), intent(inout) :: xsec
        real(wp), intent(in) :: T_fuel, T_ref
        real(wp), intent(in), optional :: alpha_D(:)
        
        integer :: g
        real(wp) :: delta_T, factor, alpha_doppler
        
        delta_T = T_fuel - T_ref
        
        ! Doppler broadening primarily affects resonance absorption
        ! σ_a(T) ≈ σ_a(T_ref) * sqrt(T_ref / T)
        
        do g = 1, xsec%n_groups
            if (g == 1) then
                ! Fast group - small Doppler effect
                alpha_doppler = -2.0e-5_wp  ! -2 pcm/K typical
            else
                ! Thermal group - larger Doppler effect
                alpha_doppler = -3.0e-5_wp  ! -3 pcm/K typical
            end if
            
            if (present(alpha_D)) then
                if (g <= size(alpha_D)) alpha_doppler = alpha_D(g)
            end if
            
            ! Temperature-dependent resonance absorption
            factor = sqrt(T_ref / max(T_fuel, 300.0_wp))
            xsec%sigma_a(g) = xsec%sigma_a(g) * factor * (1.0_wp + alpha_doppler * delta_T / 1.0e5_wp)
            
            ! Resonance fission also affected
            xsec%sigma_f(g) = xsec%sigma_f(g) * (1.0_wp + 0.5_wp * alpha_doppler * delta_T / 1.0e5_wp)
            xsec%nu_sigma_f(g) = xsec%sigma_f(g) * 2.43_wp
        end do
        
        ! Update removal
        do g = 1, xsec%n_groups
            xsec%sigma_r(g) = xsec%sigma_a(g) + sum(xsec%sigma_s(g, :)) - xsec%sigma_s(g, g)
        end do
    end subroutine xslib_apply_doppler
    
    !> Apply moderator temperature effect
    subroutine xslib_apply_moderator_temperature(xsec, T_mod, T_ref, alpha_mod)
        type(mg_xsec_t), intent(inout) :: xsec
        real(wp), intent(in) :: T_mod, T_ref
        real(wp), intent(in), optional :: alpha_mod(:)
        
        real(wp) :: delta_T, factor
        integer :: g
        
        delta_T = T_mod - T_ref
        
        ! Moderator temperature affects thermal scattering
        do g = 1, xsec%n_groups
            if (g == xsec%n_groups) then
                ! Thermal group - hardening of spectrum
                factor = 1.0_wp + 1.0e-4_wp * delta_T
                xsec%sigma_s(g, g) = xsec%sigma_s(g, g) * factor
            end if
        end do
    end subroutine xslib_apply_moderator_temperature
    
    !> Apply moderator density effect
    subroutine xslib_apply_moderator_density(xsec, rho_mod, rho_ref, alpha_rho)
        type(mg_xsec_t), intent(inout) :: xsec
        real(wp), intent(in) :: rho_mod, rho_ref
        real(wp), intent(in), optional :: alpha_rho(:)
        
        real(wp) :: ratio
        
        if (rho_ref < TOL_DEFAULT) return
        
        ! Density ratio
        ratio = rho_mod / rho_ref
        
        ! Scale moderator cross sections linearly with density
        ! (For more accurate treatment, use void correlation)
        xsec%sigma_a = xsec%sigma_a * ratio
        xsec%sigma_s = xsec%sigma_s * ratio
        xsec%sigma_r = xsec%sigma_r * ratio
        
        ! Diffusion coefficient inversely proportional to density
        xsec%D = xsec%D / ratio
    end subroutine xslib_apply_moderator_density
    
    !> Apply xenon-135 poisoning
    subroutine xslib_apply_xenon(xsec, Xe_conc)
        type(mg_xsec_t), intent(inout) :: xsec
        real(wp), intent(in) :: Xe_conc  ! atoms/barn-cm
        
        real(wp) :: sigma_a_Xe_thermal
        
        ! Xe-135 thermal absorption cross section: ~2.65e6 barns
        sigma_a_Xe_thermal = 2.65e6_wp * 1.0e-24_wp  ! Convert to cm²
        
        ! Add to thermal group absorption
        if (xsec%n_groups >= 2) then
            xsec%sigma_a(xsec%n_groups) = xsec%sigma_a(xsec%n_groups) + &
                                          sigma_a_Xe_thermal * Xe_conc
            xsec%sigma_r(xsec%n_groups) = xsec%sigma_r(xsec%n_groups) + &
                                          sigma_a_Xe_thermal * Xe_conc
        end if
    end subroutine xslib_apply_xenon
    
    !> Apply samarium-149 poisoning
    subroutine xslib_apply_samarium(xsec, Sm_conc)
        type(mg_xsec_t), intent(inout) :: xsec
        real(wp), intent(in) :: Sm_conc  ! atoms/barn-cm
        
        real(wp) :: sigma_a_Sm_thermal
        
        ! Sm-149 thermal absorption cross section: ~4.08e4 barns
        sigma_a_Sm_thermal = 4.08e4_wp * 1.0e-24_wp  ! Convert to cm²
        
        ! Add to thermal group absorption
        if (xsec%n_groups >= 2) then
            xsec%sigma_a(xsec%n_groups) = xsec%sigma_a(xsec%n_groups) + &
                                          sigma_a_Sm_thermal * Sm_conc
            xsec%sigma_r(xsec%n_groups) = xsec%sigma_r(xsec%n_groups) + &
                                          sigma_a_Sm_thermal * Sm_conc
        end if
    end subroutine xslib_apply_samarium
    
    !> Apply soluble boron
    subroutine xslib_apply_boron(xsec, boron_ppm, boron_worth)
        type(mg_xsec_t), intent(inout) :: xsec
        real(wp), intent(in) :: boron_ppm
        real(wp), intent(in), optional :: boron_worth
        
        real(wp) :: sigma_a_B10, N_B10, worth
        
        worth = -8.0_wp  ! Default: -8 pcm/ppm
        if (present(boron_worth)) worth = boron_worth
        
        ! B-10 thermal absorption: ~3840 barns
        sigma_a_B10 = 3840.0_wp * 1.0e-24_wp
        
        ! Number density from ppm (assuming water density 0.7 g/cm³)
        ! ppm = mg/kg, so for 1 ppm: ~6e-6 g/cm³ of boron
        N_B10 = boron_ppm * 6.0e-6_wp * 0.2_wp * N_AVOGADRO / 10.8_wp  ! 20% B-10 abundance
        
        ! Add to thermal absorption
        if (xsec%n_groups >= 2) then
            xsec%sigma_a(xsec%n_groups) = xsec%sigma_a(xsec%n_groups) + sigma_a_B10 * N_B10
            xsec%sigma_r(xsec%n_groups) = xsec%sigma_r(xsec%n_groups) + sigma_a_B10 * N_B10
        end if
    end subroutine xslib_apply_boron
    
    !> Apply burnup effects
    subroutine apply_burnup_effects(xsec, burnup, mat)
        type(mg_xsec_t), intent(inout) :: xsec
        real(wp), intent(in) :: burnup  ! MWd/kgU
        type(xsec_material_t), intent(in) :: mat
        
        integer :: g
        real(wp) :: factor, delta_burnup
        
        delta_burnup = burnup - mat%burnup_ref
        
        ! Simplified burnup model
        ! As fuel burns:
        ! - U-235 depletes → lower nu, sigma_f
        ! - Pu-239 builds up → partially compensates
        ! - Fission products accumulate → higher sigma_a
        
        do g = 1, xsec%n_groups
            ! Fission cross section decreases
            factor = exp(-delta_burnup / 50000.0_wp)  ! ~50 GWd/MTU half-life
            xsec%sigma_f(g) = xsec%sigma_f(g) * factor
            xsec%nu_sigma_f(g) = xsec%nu_sigma_f(g) * factor
            
            ! Absorption increases (fission products)
            factor = 1.0_wp + delta_burnup / 100000.0_wp
            xsec%sigma_a(g) = xsec%sigma_a(g) * factor
            xsec%sigma_r(g) = xsec%sigma_r(g) * factor
        end do
    end subroutine apply_burnup_effects
    
    !> Create 2-group fuel cross sections
    subroutine xslib_create_two_group_fuel(mat, enrichment)
        type(xsec_material_t), intent(out) :: mat
        real(wp), intent(in), optional :: enrichment
        
        real(wp) :: enrich
        
        enrich = 0.035_wp
        if (present(enrichment)) enrich = enrichment
        
        mat%name = "UO2_fuel"
        mat%n_groups = 2
        mat%is_fuel = .true.
        mat%enrichment = enrich
        mat%T_ref = 900.0_wp
        
        ! Allocate base cross sections
        mat%xsec_base%n_groups = 2
        allocate(mat%xsec_base%D(2))
        allocate(mat%xsec_base%sigma_t(2))
        allocate(mat%xsec_base%sigma_a(2))
        allocate(mat%xsec_base%sigma_f(2))
        allocate(mat%xsec_base%nu_sigma_f(2))
        allocate(mat%xsec_base%sigma_r(2))
        allocate(mat%xsec_base%chi(2))
        allocate(mat%xsec_base%kappa(2))
        allocate(mat%xsec_base%sigma_s(2, 2))
        allocate(mat%xsec_base%chi_d(2))
        
        ! Typical 2-group constants for LEU fuel (3.5% enrichment)
        ! Group 1: Fast (10 eV - 10 MeV)
        ! Group 2: Thermal (0 - 10 eV)
        
        ! Diffusion coefficient [cm]
        mat%xsec_base%D(1) = 1.5_wp
        mat%xsec_base%D(2) = 0.4_wp
        
        ! Absorption [cm⁻¹]
        mat%xsec_base%sigma_a(1) = 0.010_wp * (enrich / 0.035_wp)
        mat%xsec_base%sigma_a(2) = 0.085_wp * (enrich / 0.035_wp)
        
        ! Fission [cm⁻¹]
        mat%xsec_base%sigma_f(1) = 0.004_wp * (enrich / 0.035_wp)
        mat%xsec_base%sigma_f(2) = 0.065_wp * (enrich / 0.035_wp)
        
        ! Production
        mat%xsec_base%nu_sigma_f(1) = 0.010_wp * (enrich / 0.035_wp)
        mat%xsec_base%nu_sigma_f(2) = 0.158_wp * (enrich / 0.035_wp)
        
        ! Scattering matrix [cm⁻¹]
        mat%xsec_base%sigma_s(1, 1) = 0.30_wp   ! Fast → Fast
        mat%xsec_base%sigma_s(1, 2) = 0.020_wp  ! Fast → Thermal (slowing down)
        mat%xsec_base%sigma_s(2, 1) = 0.0_wp    ! Thermal → Fast (no upscatter for fuel)
        mat%xsec_base%sigma_s(2, 2) = 0.90_wp   ! Thermal → Thermal
        
        ! Removal
        mat%xsec_base%sigma_r(1) = mat%xsec_base%sigma_a(1) + mat%xsec_base%sigma_s(1, 2)
        mat%xsec_base%sigma_r(2) = mat%xsec_base%sigma_a(2)
        
        ! Fission spectrum
        mat%xsec_base%chi(1) = 1.0_wp
        mat%xsec_base%chi(2) = 0.0_wp
        
        ! Energy per fission [MeV]
        mat%xsec_base%kappa(1) = 200.0_wp
        mat%xsec_base%kappa(2) = 200.0_wp
        
        ! Delayed neutron spectrum (similar to prompt)
        mat%xsec_base%chi_d(1) = 1.0_wp
        mat%xsec_base%chi_d(2) = 0.0_wp
        
        ! Allocate feedback coefficients
        allocate(mat%alpha_D(2))
        mat%alpha_D(1) = -1.5_wp  ! pcm/K (fast Doppler)
        mat%alpha_D(2) = -3.0_wp  ! pcm/K (thermal Doppler)
    end subroutine xslib_create_two_group_fuel
    
    !> Create 2-group moderator cross sections
    subroutine xslib_create_two_group_moderator(mat)
        type(xsec_material_t), intent(out) :: mat
        
        mat%name = "H2O_moderator"
        mat%n_groups = 2
        mat%is_fuel = .false.
        mat%rho_ref = 0.7_wp  ! g/cm³ at operating conditions
        
        ! Allocate
        mat%xsec_base%n_groups = 2
        allocate(mat%xsec_base%D(2))
        allocate(mat%xsec_base%sigma_t(2))
        allocate(mat%xsec_base%sigma_a(2))
        allocate(mat%xsec_base%sigma_f(2))
        allocate(mat%xsec_base%nu_sigma_f(2))
        allocate(mat%xsec_base%sigma_r(2))
        allocate(mat%xsec_base%chi(2))
        allocate(mat%xsec_base%kappa(2))
        allocate(mat%xsec_base%sigma_s(2, 2))
        allocate(mat%xsec_base%chi_d(2))
        
        ! Water (non-multiplying)
        mat%xsec_base%D(1) = 1.2_wp
        mat%xsec_base%D(2) = 0.16_wp
        
        mat%xsec_base%sigma_a(1) = 0.0008_wp
        mat%xsec_base%sigma_a(2) = 0.0020_wp
        
        mat%xsec_base%sigma_f = 0.0_wp
        mat%xsec_base%nu_sigma_f = 0.0_wp
        
        ! Strong scattering (H-1)
        mat%xsec_base%sigma_s(1, 1) = 0.50_wp
        mat%xsec_base%sigma_s(1, 2) = 0.10_wp   ! Moderation
        mat%xsec_base%sigma_s(2, 1) = 0.001_wp  ! Small upscatter
        mat%xsec_base%sigma_s(2, 2) = 1.50_wp
        
        mat%xsec_base%sigma_r(1) = mat%xsec_base%sigma_a(1) + mat%xsec_base%sigma_s(1, 2)
        mat%xsec_base%sigma_r(2) = mat%xsec_base%sigma_a(2) + mat%xsec_base%sigma_s(2, 1)
        
        mat%xsec_base%chi = 0.0_wp
        mat%xsec_base%kappa = 0.0_wp
        mat%xsec_base%chi_d = 0.0_wp
        
        ! Moderator feedback
        allocate(mat%alpha_mod(2))
        allocate(mat%alpha_rho(2))
        
        mat%alpha_mod(1) = 0.0_wp
        mat%alpha_mod(2) = 1.0e-4_wp  ! Small positive (spectral hardening)
        
        mat%alpha_rho(1) = 50.0_wp    ! pcm/(g/cm³) - positive
        mat%alpha_rho(2) = 80.0_wp
    end subroutine xslib_create_two_group_moderator
    
    !> Compute Doppler coefficient
    function compute_doppler_coefficient(xsec_low, xsec_high, delta_T) result(alpha_D)
        type(mg_xsec_t), intent(in) :: xsec_low, xsec_high
        real(wp), intent(in) :: delta_T
        real(wp) :: alpha_D
        
        real(wp) :: delta_sigma_a
        
        ! Simplified: use thermal group absorption change
        delta_sigma_a = xsec_high%sigma_a(xsec_high%n_groups) - &
                       xsec_low%sigma_a(xsec_low%n_groups)
        
        alpha_D = (delta_sigma_a / xsec_low%sigma_a(xsec_low%n_groups)) / delta_T * 1.0e5_wp  ! pcm/K
    end function compute_doppler_coefficient
    
    !> Compute moderator temperature coefficient
    function compute_moderator_coefficient(k_low, k_high, delta_T) result(alpha_mod)
        real(wp), intent(in) :: k_low, k_high, delta_T
        real(wp) :: alpha_mod
        
        alpha_mod = (k_high - k_low) / k_low / delta_T * 1.0e5_wp  ! pcm/K
    end function compute_moderator_coefficient
    
    !> Compute xenon worth
    function compute_xenon_worth(k_clean, k_xenon) result(rho_xe)
        real(wp), intent(in) :: k_clean, k_xenon
        real(wp) :: rho_xe
        
        ! Reactivity in pcm
        rho_xe = (k_xenon - k_clean) / (k_clean * k_xenon) * 1.0e5_wp
    end function compute_xenon_worth
    
    !> Compute samarium worth
    function compute_samarium_worth(k_clean, k_sm) result(rho_sm)
        real(wp), intent(in) :: k_clean, k_sm
        real(wp) :: rho_sm
        
        rho_sm = (k_sm - k_clean) / (k_clean * k_sm) * 1.0e5_wp
    end function compute_samarium_worth

end module cross_sections