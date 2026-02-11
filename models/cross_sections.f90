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
    use constants, only: TOL_DEFAULT, N_AVOGADRO, PI
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
    public :: compute_xenon_worth, compute_samarium_worth, xslib_get_material
    public :: xslib_list_materials, xslib_apply_control_rod
    
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
            ! Initialize the material structure
            mat%name = "UO2_35"
            mat%n_groups = 2
            mat%is_fuel = .true.
            mat%enrichment = 0.035_wp
            mat%T_ref = 900.0_wp
            mat%rho_ref = 10.97_wp  ! UO2 density g/cm³
            
            ! Create the base cross-sections - pass xsec_base, not mat
            call xslib_create_two_group_fuel(mat%xsec_base, enrichment=0.035_wp)
            
            ! Allocate feedback coefficients
            allocate(mat%alpha_D(2))
            allocate(mat%alpha_mod(2))
            allocate(mat%alpha_rho(2))
            allocate(mat%alpha_void(2))
            
            ! Set typical feedback coefficients for UO2
            mat%alpha_D = [-2.0e-5_wp, -3.0e-5_wp]  ! Doppler (pcm/K)
            mat%alpha_mod = [0.0_wp, 0.0_wp]        ! Not applicable for fuel
            mat%alpha_rho = [0.0_wp, 0.0_wp]        ! Not applicable for fuel
            mat%alpha_void = [0.0_wp, 0.0_wp]       ! Not applicable for fuel
            
            call xslib_add_material(library, mat)
            
            ! Water moderator
            mat%name = "H2O"
            mat%n_groups = 2
            mat%is_fuel = .false.
            mat%T_ref = 560.0_wp    ! ~287°C
            mat%rho_ref = 0.74_wp   ! Water density at 7 MPa, 287°C
            
            ! Create the base cross-sections - pass xsec_base, not mat
            call xslib_create_two_group_moderator(mat%xsec_base)
            
            ! Allocate feedback coefficients
            if (.not. allocated(mat%alpha_D)) allocate(mat%alpha_D(2))
            if (.not. allocated(mat%alpha_mod)) allocate(mat%alpha_mod(2))
            if (.not. allocated(mat%alpha_rho)) allocate(mat%alpha_rho(2))
            if (.not. allocated(mat%alpha_void)) allocate(mat%alpha_void(2))
            
            ! Set typical feedback coefficients for water
            mat%alpha_D = [0.0_wp, 0.0_wp]          ! Not applicable for moderator
            mat%alpha_mod = [0.0_wp, 1.0e-4_wp]     ! Temperature effect (pcm/K)
            mat%alpha_rho = [-10.0_wp, -50.0_wp]    ! Density effect (pcm/(g/cm³))
            mat%alpha_void = [-10.0_wp, -100.0_wp]  ! Void effect (pcm/% void)
            
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
        real(wp), intent(in) :: Xe_conc  ! Input is atoms/barn-cm

        real(wp) :: sigma_a_Xe_thermal

        ! Xe-135 thermal absorption cross section: ~2.65e6 barns
        ! CORRECTION: Keep in barns to match atoms/barn-cm input
        sigma_a_Xe_thermal = 2.65e6_wp 

        ! Add to thermal group absorption
        if (xsec%n_groups >= 2) then
            xsec%sigma_a(xsec%n_groups) = xsec%sigma_a(xsec%n_groups) + &
                                          sigma_a_Xe_thermal * Xe_conc
            
            ! Also update removal cross-section (absorption is part of removal)
            xsec%sigma_r(xsec%n_groups) = xsec%sigma_r(xsec%n_groups) + &
                                          sigma_a_Xe_thermal * Xe_conc
        end if
    end subroutine xslib_apply_xenon

    !> Apply samarium-149 poisoning
    subroutine xslib_apply_samarium(xsec, Sm_conc)
        type(mg_xsec_t), intent(inout) :: xsec
        real(wp), intent(in) :: Sm_conc  ! Input is atoms/barn-cm

        real(wp) :: sigma_a_Sm_thermal

        ! Sm-149 thermal absorption cross section: ~4.08e4 barns
        ! CORRECTION: Keep in barns
        sigma_a_Sm_thermal = 4.08e4_wp 

        ! Add to thermal group absorption
        if (xsec%n_groups >= 2) then
            xsec%sigma_a(xsec%n_groups) = xsec%sigma_a(xsec%n_groups) + &
                                          sigma_a_Sm_thermal * Sm_conc
            xsec%sigma_r(xsec%n_groups) = xsec%sigma_r(xsec%n_groups) + &
                                          sigma_a_Sm_thermal * Sm_conc
        end if
    end subroutine xslib_apply_samarium
    
    !> Apply Control Rods
    !> Adds strong absorption to represent a B4C or Hafnium blade
    subroutine xslib_apply_control_rod(xsec, rod_fraction)
        type(mg_xsec_t), intent(inout) :: xsec
        real(wp), intent(in) :: rod_fraction ! 0.0 (empty) to 1.0 (full blade)

        real(wp) :: delta_sigma_a_fast, delta_sigma_a_thermal
        
        ! Macroscopic worth constants [cm^-1]
        ! These are "black" absorbers relative to fuel
        delta_sigma_a_fast    = 0.015_wp   ! Small effect in fast group
        delta_sigma_a_thermal = 0.150_wp   ! Massive effect in thermal group

        ! Linearly weight by how much rod is in this specific node
        if (rod_fraction > 0.0_wp) then
            ! Update Fast Group
            xsec%sigma_a(1) = xsec%sigma_a(1) + (delta_sigma_a_fast * rod_fraction)
            xsec%sigma_r(1) = xsec%sigma_r(1) + (delta_sigma_a_fast * rod_fraction)
            
            ! Update Thermal Group
            xsec%sigma_a(2) = xsec%sigma_a(2) + (delta_sigma_a_thermal * rod_fraction)
            xsec%sigma_r(2) = xsec%sigma_r(2) + (delta_sigma_a_thermal * rod_fraction)
        end if
    end subroutine xslib_apply_control_rod

    !> Apply soluble boron
    subroutine xslib_apply_boron(xsec, boron_ppm, boron_worth)
        type(mg_xsec_t), intent(inout) :: xsec
        real(wp), intent(in) :: boron_ppm
        real(wp), intent(in), optional :: boron_worth

        real(wp) :: sigma_a_B10, N_B10

        ! B-10 thermal absorption: ~3840 barns
        ! CORRECTION: Keep in barns
        sigma_a_B10 = 3840.0_wp 

        ! Calculate Number Density
        ! ppm = parts per million by mass (g B / 1e6 g Water)
        ! Density of water approx 0.74 g/cm^3 at operating conditions
        ! N = (ppm * 1e-6) * rho_water * Abundance * N_Avogadro / AtomicWeight
        
        ! Note: We multiply by 1.0e-24 to convert from atoms/cm^3 to atoms/barn-cm
        N_B10 = (boron_ppm * 1.0e-6_wp) * &        ! mass fraction
                0.74_wp * &                        ! water density (g/cm^3)
                0.20_wp * &                        ! B-10 abundance (20%)
                (N_AVOGADRO / 10.81_wp) * &        ! moles/g
                1.0e-24_wp                         ! CONVERSION: cm^-3 -> (barn-cm)^-1

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

    subroutine xslib_get_material(library, material_name, material, found)
        type(xsec_library_t), intent(in)  :: library
        character(len=*),     intent(in)  :: material_name
        type(xsec_material_t), intent(out) :: material
        logical,               intent(out) :: found

        integer :: i
        found = .false.

        do i = 1, library%n_materials
            if (trim(adjustl(library%materials(i)%name)) == trim(adjustl(material_name))) then
                material = library%materials(i)
                found = .true.
                return
            end if
        end do

        ! not found: print diagnostics (helpful during debugging) but do not leave material uninitialised
        print *, 'xslib_get_material: material "', trim(material_name), '" not found in library.'
        print *, 'Available materials:'
        do i = 1, library%n_materials
            print *, '  - ', trim(library%materials(i)%name)
        end do
    end subroutine xslib_get_material

    subroutine xslib_list_materials(library)
        type(xsec_library_t), intent(in) :: library
        integer :: i
        print *, 'xslib: material list (', library%n_materials, '):'
        do i = 1, library%n_materials
            print *, '  ', trim(library%materials(i)%name)
        end do
    end subroutine xslib_list_materials

    subroutine xslib_create_two_group_fuel(xsec, enrichment)
        type(mg_xsec_t), intent(out) :: xsec
        real(wp), intent(in) :: enrichment
        
        real(wp) :: f_U235, f_U238
        real(wp) :: V_fuel, V_mod, f_fuel  ! Volume fractions
        real(wp) :: pitch, r_fuel

        real(wp) :: sigma_a_fuel, sigma_f_fuel, N_U
        real(wp) :: sigma_a_water

        xsec%n_groups = 2
        
        allocate(xsec%D(2))
        allocate(xsec%sigma_t(2))
        allocate(xsec%sigma_a(2))
        allocate(xsec%sigma_f(2))
        allocate(xsec%nu_sigma_f(2))
        allocate(xsec%sigma_r(2))
        allocate(xsec%chi(2))
        allocate(xsec%kappa(2))
        allocate(xsec%sigma_s(2, 2))
        allocate(xsec%chi_d(2))
        
        ! Enrichment fractions
        f_U235 = enrichment
        f_U238 = 1.0_wp - enrichment
        
        ! Compute fuel volume fraction for homogenization
        pitch = 0.0163_wp    ! 1.63 cm BWR pin pitch
        r_fuel = 0.0052_wp   ! 5.2 mm fuel radius
        
        ! Fuel volume fraction in lattice
        ! (Requires PI to be imported or defined)
        V_fuel = PI * r_fuel**2
        V_mod = pitch**2 - V_fuel
        f_fuel = V_fuel / (V_fuel + V_mod)
        
        print '(A,F6.4)', "  Fuel volume fraction: ", f_fuel
        
        ! ========================================================================
        ! GROUP 1: FAST (E > 0.625 eV)
        ! ========================================================================
        
        ! Homogenized diffusion coefficient
        xsec%D(1) = f_fuel * 1.50_wp + (1.0_wp - f_fuel) * 1.30_wp
        
        ! Total cross-section
        xsec%sigma_t(1) = 0.28_wp
        
        ! CRITICAL: Homogenized absorption
        xsec%sigma_a(1) = f_fuel * 0.001_wp + (1.0_wp - f_fuel) * 0.0001_wp
        
        ! Fast fission
        xsec%sigma_f(1) = f_fuel * f_U235 * 0.0053_wp
        xsec%nu_sigma_f(1) = 2.50_wp * xsec%sigma_f(1)
        
        ! All fast neutrons born in fast group
        xsec%chi(1) = 1.0_wp
        xsec%kappa(1) = 200.0_wp
        
        ! ========================================================================
        ! GROUP 2: THERMAL (E < 0.625 eV)  
        ! ========================================================================
        
        ! Homogenized diffusion coefficient
        xsec%D(2) = f_fuel * 0.40_wp + (1.0_wp - f_fuel) * 0.16_wp
        
        xsec%sigma_t(2) = 1.50_wp
        
        ! CRITICAL FIX: Homogenized macroscopic cross-sections
        ! These are averaged over fuel + moderator volumes
        
        ! Fuel contribution (using proper number density)
        N_U = 0.02354_wp  ! atoms/barn-cm
        
        ! Macroscopic XS for pure fuel [cm^-1]
        sigma_a_fuel = N_U * (f_U235 * 681.0_wp + f_U238 * 2.7_wp)
        sigma_f_fuel = N_U * f_U235 * 584.0_wp
        
        ! Water absorption: σ_a,water ≈ 0.0196 cm^-1
        sigma_a_water = 0.0196_wp
        
        ! HOMOGENIZED cross-sections (volume-weighted)
        xsec%sigma_a(2) = f_fuel * sigma_a_fuel + (1.0_wp - f_fuel) * sigma_a_water
        xsec%sigma_f(2) = f_fuel * sigma_f_fuel
        
        ! Production
        xsec%nu_sigma_f(2) = 2.42_wp * xsec%sigma_f(2)
        
        xsec%chi(2) = 0.0_wp
        xsec%kappa(2) = 200.0_wp
        
        ! ========================================================================
        ! SCATTERING MATRIX [cm^-1]
        ! ========================================================================
        
        ! Fast-to-fast (homogenized)
        xsec%sigma_s(1, 1) = f_fuel * 0.260_wp + (1.0_wp - f_fuel) * 1.05_wp
        
        ! CRITICAL: Fast-to-thermal slowing down (mostly in water!)
        xsec%sigma_s(1, 2) = (1.0_wp - f_fuel) * 0.048_wp + f_fuel * 0.012_wp
        
        ! No upscatter
        xsec%sigma_s(2, 1) = 0.0_wp
        
        ! Thermal-to-thermal (mostly water)
        xsec%sigma_s(2, 2) = f_fuel * 1.35_wp + (1.0_wp - f_fuel) * 3.48_wp
        
        ! Compute removal cross-section
        xsec%sigma_r(1) = xsec%sigma_a(1) + xsec%sigma_s(1, 2)
        xsec%sigma_r(2) = xsec%sigma_a(2)
        
        ! Delayed neutron data
        xsec%beta_total = 0.0065_wp
        xsec%chi_d(1) = 1.0_wp
        xsec%chi_d(2) = 0.0_wp
        
        ! Verification printout
        print *, "Created HOMOGENIZED UO2 fuel cross-sections:"
        print '(A,F6.4)', "  Enrichment:        ", enrichment * 100.0_wp, "%"
        print '(A,ES11.4)', "  σ_a(thermal):      ", xsec%sigma_a(2), " cm^-1"
        print '(A,ES11.4)', "  σ_f(thermal):      ", xsec%sigma_f(2), " cm^-1"
        print '(A,ES11.4)', "  ν*σ_f(thermal):    ", xsec%nu_sigma_f(2), " cm^-1"
        print '(A,F6.3)', "  k-infinity est:    ", xsec%nu_sigma_f(2) / xsec%sigma_a(2)
        
    end subroutine xslib_create_two_group_fuel

    subroutine xslib_create_two_group_moderator(xsec)
        type(mg_xsec_t), intent(out) :: xsec
        
        xsec%n_groups = 2
        
        allocate(xsec%D(2))
        allocate(xsec%sigma_t(2))
        allocate(xsec%sigma_a(2))
        allocate(xsec%sigma_f(2))
        allocate(xsec%nu_sigma_f(2))
        allocate(xsec%sigma_r(2))
        allocate(xsec%chi(2))
        allocate(xsec%kappa(2))
        allocate(xsec%sigma_s(2, 2))
        allocate(xsec%chi_d(2))
        
        ! ========================================================================
        ! GROUP 1: FAST (E > 0.625 eV)
        ! ========================================================================
        
        xsec%D(1) = 1.30_wp
        xsec%sigma_t(1) = 1.10_wp
        
        ! Water absorption in fast group (very small)
        xsec%sigma_a(1) = 0.0001_wp  ! Mostly potential scattering
        
        ! No fission in water
        xsec%sigma_f(1) = 0.0_wp
        xsec%nu_sigma_f(1) = 0.0_wp
        xsec%chi(1) = 0.0_wp
        xsec%kappa(1) = 0.0_wp
        
        ! ========================================================================
        ! GROUP 2: THERMAL (E < 0.625 eV)
        ! ========================================================================
        
        xsec%D(2) = 0.16_wp
        xsec%sigma_t(2) = 3.50_wp
        
        ! Water thermal absorption (THIS IS CRITICAL!)
        ! H-1: σ_a ~ 0.332 barns at 0.0253 eV
        ! At liquid water density (0.74 g/cm³ at BWR conditions):
        !   N_H = 0.74 * 6.022e23 / 18 * 2 = 4.94e22 atoms/cm³
        !   σ_a = N * σ = 4.94e22 * 0.332e-24 = 0.0164 cm⁻¹
        xsec%sigma_a(2) = 0.0196_wp  ! Includes temperature correction
        
        ! No fission in water
        xsec%sigma_f(2) = 0.0_wp
        xsec%nu_sigma_f(2) = 0.0_wp
        xsec%chi(2) = 0.0_wp
        xsec%kappa(2) = 0.0_wp
        
        ! ========================================================================
        ! SCATTERING MATRIX
        ! ========================================================================
        
        ! Fast-to-fast (elastic scattering)
        xsec%sigma_s(1, 1) = 1.05_wp
        
        ! Fast-to-thermal (slowing down - VERY IMPORTANT for moderation!)
        xsec%sigma_s(1, 2) = 0.048_wp
        
        ! Thermal-to-fast (negligible upscatter)
        xsec%sigma_s(2, 1) = 0.0_wp
        
        ! Thermal-to-thermal (elastic + inelastic)
        xsec%sigma_s(2, 2) = 3.48_wp
        
        ! Compute removal
        xsec%sigma_r(1) = xsec%sigma_a(1) + xsec%sigma_s(1, 2)
        xsec%sigma_r(2) = xsec%sigma_a(2)
        
        ! No delayed neutrons from water
        xsec%beta_total = 0.0_wp
        xsec%chi_d = 0.0_wp
        
        print *, "Created water moderator cross-sections:"
        print '(A,ES11.4)', "  σ_a(fast):         ", xsec%sigma_a(1), " cm⁻¹"
        print '(A,ES11.4)', "  σ_a(thermal):      ", xsec%sigma_a(2), " cm⁻¹"
        print '(A,ES11.4)', "  σ_s(1→2):          ", xsec%sigma_s(1, 2), " cm⁻¹"
        
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