! Fuel subsystem.
!
! Replaces the legacy 1-D bank-average rod fan-out:
!   bank_avg = mean(blade_insertion)        ← all cells see the same insertion
! with per-blade spatial fan-out:
!   blade_for_cell(i, j) ∈ [1, CRD_N_BLADES]
!   per-cell insertion = blade_insertion(blade_for_cell(i, j))
!
! Blade layout: 137 control rod drives placed on a 12 × 12 = 144 grid,
! the seven outermost corner positions dropped to hit the BWR/4 count
! (4 outer corners + 3 from the second ring; deterministic). Each fuel
! cell maps to its nearest blade by Euclidean distance. Cells outside
! the core radial mask (r > D / `core_mask_factor`) have blade_for_cell
! = 0 (water bypass; no rod XS contribution).
!
! NOTE on `void_reactivity_factor`: moved here from `bwr_c_interface`
! because it's a reactivity coefficient (not a vessel concern). The
! `alpha_doppler` / `alpha_void` scalars from the legacy god-state were
! never actually consumed by the XS feedback path and are dropped.
!
! NOTE on operator reactivity perturbation: `apply_reactivity_perturbation`
! also lives here — it directly mutates the live `nu_sigma_f` array each
! tick, which is squarely fuel-physics scope.
module fuel
    use kinds, only: wp
    use multigroup_diffusion, only: mg_state_t, mg_xsec_t, mg_set_cross_sections
    use cross_sections, only: xsec_library_t, xsec_material_t, &
                              xslib_get_material, xslib_get_xsec, &
                              xslib_apply_control_rod
    use crd, only: CRD_N_BLADES
    use two_phase_flow, only: two_phase_state_t, water_properties_t
    implicit none
    private

    public :: fuel_state_t, fuel_config_t
    public :: fuel_command_t, fuel_observation_t
    public :: fuel_init, fuel_init_geometry, fuel_destroy
    public :: fuel_step, fuel_observe, fuel_apply
    public :: fuel_setup_geometry, fuel_apply_to_xs
    public :: fuel_apply_reactivity_perturbation
    public :: void_reactivity_factor
    public :: fuel_power_to_volumetric_W_m3
    public :: fuel_power_to_surface_q
    public :: fuel_compute_convection_coefficients
    public :: dittus_boelter_h, jens_lottes_h
    public :: FUEL_BLADE_GRID_N

    ! 12 × 12 = 144 blade-position grid, less 7 outermost corners to
    ! land at 137. The trim mask is deterministic — see `init_blade_xy`.
    integer, parameter :: FUEL_BLADE_GRID_N = 12

    type :: fuel_config_t
        ! GE-14 10×10 BWR/4 fuel lattice (spec 02).
        real(wp) :: rod_pitch_m       = 0.01295_wp ! 12.95 mm rod pitch
        real(wp) :: pin_radius_m      = 0.00437_wp ! 4.37 mm pellet radius
        ! Cells with r > core_diameter / core_mask_factor are bypass
        ! water (no fuel, no control rod). Magic value 2.2 retained from
        ! the legacy setup_geometry for behavioural continuity.
        real(wp) :: core_mask_factor  = 2.2_wp
        ! Reactivity-coefficient placeholders. These were on the legacy
        ! bwr_state_t but never consumed by the XS-feedback path; kept
        ! here for forward use when chapter 7 (RFC) lands.
        real(wp) :: alpha_doppler     = -3.5_wp
        real(wp) :: alpha_void        = -80.0_wp
        ! Reference void fraction the rf_void scaling normalises against.
        real(wp) :: void_ref_fraction = 0.30_wp
    end type fuel_config_t

    type :: fuel_command_t
        ! No fuel-level operator commands in phase 1. Rod motion arrives
        ! through crd, fuel/coolant T arrives through heat + two_phase.
        integer :: reserved = 0
    end type fuel_command_t

    type :: fuel_observation_t
        real(wp) :: time            = 0.0_wp
        integer  :: n_steps         = 0
        real(wp) :: avg_fuel_temp_K = 0.0_wp
        real(wp) :: max_fuel_temp_K = 0.0_wp
    end type fuel_observation_t

    type :: fuel_state_t
        real(wp) :: time     = 0.0_wp
        integer  :: n_steps  = 0

        ! Cached lattice/grid geometry.
        integer  :: nx = 0
        integer  :: ny = 0
        real(wp) :: dx = 0.0_wp
        real(wp) :: dy = 0.0_wp
        real(wp) :: dz = 0.0_wp
        real(wp) :: core_diameter_m = 0.0_wp

        ! Per-cell blade index. 0 means water-bypass (no rod XS).
        integer, allocatable :: blade_for_cell(:, :)
        ! 2-D blade positions in metres. Shape (2, CRD_N_BLADES); row 1 = x, 2 = y.
        real(wp), allocatable :: blade_xy(:, :)
        ! Fuel/water mask (true = fuel cell), for both fuel_setup_geometry
        ! and downstream observability.
        logical,  allocatable :: is_fuel_cell(:, :)

        ! Optional caches updated each tick.
        real(wp) :: avg_fuel_temp_K = 0.0_wp
        real(wp) :: max_fuel_temp_K = 0.0_wp
    end type fuel_state_t

contains

    ! ── Power unit-conversion helpers (reorg step 10) ────────────────────
    ! mg_get_power returns power density in W/cm³. Downstream kernels
    ! want different units — and the legacy bwr_c_interface used to feed
    ! the raw W/cm³ array straight into `heat%q` (which the heat-transfer
    ! docstring declares as W/m³), suppressing fuel heating by a factor
    ! of 10⁶. These helpers make the conversion explicit at the call
    ! sites and let unit tests pin the contract.

    pure function fuel_power_to_volumetric_W_m3(power_W_cm3) result(q_W_m3)
        ! W/cm³ → W/m³. Factor 10⁶ (1 m³ = 10⁶ cm³). The heat kernel's
        ! `heat%q` field expects this unit per the heat_transfer.f90
        ! docstring.
        real(wp), intent(in) :: power_W_cm3(:, :, :)
        real(wp) :: q_W_m3(size(power_W_cm3, 1), size(power_W_cm3, 2), &
                           size(power_W_cm3, 3))
        q_W_m3 = power_W_cm3 * 1.0e6_wp
    end function fuel_power_to_volumetric_W_m3

    pure function fuel_power_to_surface_q(power_W_cm3, dx, dy, heated_perimeter) &
                                                                  result(q_pp_W_m2)
        ! W/cm³ → W/m² at the heated surface for the two-phase model:
        !   q'' = (power_density × cell_volume) / heated_surface_area
        !       = power_density [W/cm³] × 10⁶ × dx [m] × dy [m]
        !         / heated_perimeter [m]
        ! Guard against zero perimeter (cells with no heated surface).
        real(wp), intent(in) :: power_W_cm3(:, :, :)
        real(wp), intent(in) :: dx, dy
        real(wp), intent(in) :: heated_perimeter(:, :, :)
        real(wp) :: q_pp_W_m2(size(power_W_cm3, 1), size(power_W_cm3, 2), &
                              size(power_W_cm3, 3))
        real(wp), parameter :: P_FLOOR = 1.0e-9_wp
        q_pp_W_m2 = power_W_cm3 * 1.0e6_wp * dx * dy &
                  / max(heated_perimeter, P_FLOOR)
    end function fuel_power_to_surface_q

    ! ── Heat ↔ two-phase coupling helpers (reorg step 11) ────────────────
    ! Computes per-cell convective coefficient and bulk fluid temperature
    ! from the two-phase state so the heat kernel can stop hard-coding
    ! T_FLUID = 558 K and H_CONV = 30000 W/m²·K. Two-regime model:
    !   • single-phase forced convection  → Dittus-Boelter
    !   • saturated nucleate boiling      → Jens-Lottes
    ! The DNB transition (post-CHF film boiling) is deliberately not
    ! modelled at this step — TODO when CHF coupling lands.

    pure function dittus_boelter_h(G, D_h, props) result(h)
        ! Single-phase forced-convection h.
        !   Nu = 0.023 · Re^0.8 · Pr^0.4    (heating)
        !   h  = Nu · k_l / D_h
        ! For 0 < Re < ~10⁴ falls back to Nu = 4.36 (fully-developed laminar
        ! pipe flow with constant heat flux) so the formula stays positive
        ! and bounded at low mass flux. Re ≥ ~10⁴ recovers the standard
        ! Dittus-Boelter relation.
        real(wp), intent(in) :: G        ! mass flux [kg/m²·s]
        real(wp), intent(in) :: D_h      ! hydraulic diameter [m]
        type(water_properties_t), intent(in) :: props
        real(wp) :: h

        real(wp) :: Re, Pr, Nu
        real(wp), parameter :: MU_FLOOR = 1.0e-9_wp
        real(wp), parameter :: K_FLOOR  = 1.0e-9_wp
        real(wp), parameter :: D_FLOOR  = 1.0e-9_wp
        real(wp), parameter :: NU_LAMINAR = 4.36_wp
        real(wp), parameter :: RE_LAM_TURB = 1.0e4_wp

        Re = abs(G) * D_h / max(props%mu_l, MU_FLOOR)
        Pr = props%mu_l * props%cp_l / max(props%k_l, K_FLOOR)

        if (Re < RE_LAM_TURB) then
            Nu = NU_LAMINAR
        else
            Nu = 0.023_wp * Re**0.8_wp * Pr**0.4_wp
        end if

        h = Nu * props%k_l / max(D_h, D_FLOOR)
    end function dittus_boelter_h

    pure function jens_lottes_h(q_pp, p) result(h)
        ! Saturated nucleate-boiling h via Jens-Lottes:
        !   ΔT_sat = 25 · (q''/1e6)^0.25 · exp(-p / 6.2e6)   [K]
        !   h      = q'' / ΔT_sat
        ! q'' in W/m², p in Pa. The exp factor encodes the pressure
        ! dependence; at p = 6.89 MPa (1000 psia, typical BWR dome) it
        ! gives the historical Jens-Lottes calibration. At low q'' the
        ! formula goes singular — clamp q'' to a small positive value
        ! and let Dittus-Boelter handle the no-boiling regime.
        real(wp), intent(in) :: q_pp     ! [W/m²]
        real(wp), intent(in) :: p        ! [Pa]
        real(wp) :: h

        real(wp) :: dT_sat
        real(wp), parameter :: Q_FLOOR  = 1.0_wp   ! W/m²; below this h≈0
        real(wp), parameter :: DT_FLOOR = 1.0e-3_wp

        if (q_pp < Q_FLOOR) then
            h = 0.0_wp
            return
        end if

        dT_sat = 25.0_wp * (q_pp / 1.0e6_wp)**0.25_wp &
                * exp(-p / 6.2e6_wp)
        h = q_pp / max(dT_sat, DT_FLOOR)
    end function jens_lottes_h

    subroutine fuel_compute_convection_coefficients(thermalhydraulics, T_fluid, h_conv, &
                                                    T_bulk)
        ! Per-cell coolant temperature and convective coefficient. The
        ! boiling flag from two_phase (set when quality ≥ 0) picks the
        ! regime; on each side of the regime boundary we max with a
        ! Dittus-Boelter floor so that the coefficient never drops below
        ! the single-phase value even in low-q'' nucleate-boiling cells.
        !
        ! T_bulk (optional) is the bulk coolant temperature used in the
        ! single-phase (subcooled) regime. Without it the helper falls
        ! back to T_sat — the legacy behaviour, which is correct for hot
        ! at-power operation where the bulk fluid sits just below sat,
        ! but is *wrong* for cold standby where the actual coolant is
        ! deeply subcooled. Boiling cells always use T_sat regardless.
        type(two_phase_state_t), intent(in) :: thermalhydraulics
        real(wp), intent(out) :: T_fluid(thermalhydraulics%nx, &
                                         thermalhydraulics%ny, &
                                         thermalhydraulics%nz)
        real(wp), intent(out) :: h_conv (thermalhydraulics%nx, &
                                         thermalhydraulics%ny, &
                                         thermalhydraulics%nz)
        real(wp), intent(in), optional :: T_bulk

        integer :: i, j, k
        real(wp) :: h_sp, h_nb
        real(wp) :: T_sat_cell

        do k = 1, thermalhydraulics%nz
            do j = 1, thermalhydraulics%ny
                do i = 1, thermalhydraulics%nx
                    T_sat_cell = thermalhydraulics%props(i, j, k)%T_sat

                    if (thermalhydraulics%boiling(i, j, k)) then
                        ! Saturated boiling regime: mixture sits at T_sat.
                        T_fluid(i, j, k) = T_sat_cell
                    else if (present(T_bulk)) then
                        ! Subcooled single-phase: bulk fluid is below sat.
                        ! Clamp to T_sat as a safety against callers that
                        ! pass an averaged temp accidentally above sat.
                        T_fluid(i, j, k) = min(T_bulk, T_sat_cell)
                    else
                        ! Legacy fallback (no bulk temp supplied).
                        T_fluid(i, j, k) = T_sat_cell
                    end if

                    h_sp = dittus_boelter_h( &
                        thermalhydraulics%mass_flux(i, j, k), &
                        thermalhydraulics%diameter(i, j, k),  &
                        thermalhydraulics%props(i, j, k))

                    if (thermalhydraulics%boiling(i, j, k)) then
                        h_nb = jens_lottes_h( &
                            thermalhydraulics%heat_flux(i, j, k), &
                            thermalhydraulics%pressure(i, j, k))
                        h_conv(i, j, k) = max(h_sp, h_nb)
                    else
                        h_conv(i, j, k) = h_sp
                    end if
                end do
            end do
        end do
    end subroutine fuel_compute_convection_coefficients

    pure function void_reactivity_factor(avg_void_frac) result(rf)
        ! Legacy polynomial fit from the pre-reorg `bwr_c_interface`.
        ! Inputs: void fraction in [0, 1]. Output is a moderator-density
        ! scaling factor used in fuel_apply_to_xs.
        real(wp), intent(in) :: avg_void_frac
        real(wp) :: rf, dens2
        dens2 = (1.0_wp - min(1.0_wp, max(0.0_wp, avg_void_frac))) * 1000.0_wp
        rf = -5.84e-10_wp        * dens2**3  &
             + 1.35422e-7_wp     * dens2**2  &
             + 1.358252042e-3_wp * dens2     &
             + 0.090469568418_wp
    end function void_reactivity_factor

    subroutine init_blade_xy(state, config)
        ! Lay out the 137 control rod blades on a 12 × 12 grid centred on
        ! the core. Drop the 7 outermost corner positions (4 corners +
        ! 3 second-ring corner-adjacent) so the count matches CRD_N_BLADES.
        ! Positions are in metres relative to the (0, 0) cell origin —
        ! matches the cell-centre convention used by fuel_setup_geometry.
        type(fuel_state_t),  intent(inout) :: state
        type(fuel_config_t), intent(in)    :: config

        integer  :: ig, jg, b, dropped
        real(wp) :: blade_pitch, x_origin, y_origin
        logical  :: skip
        integer  :: drop_positions(7, 2)

        ! Drop these (1-based grid (ig, jg)) — chosen for deterministic
        ! symmetry: 4 outer corners + 3 adjacent ring positions.
        drop_positions = reshape([ 1,  1,  &
                                  12,  1,  &
                                   1, 12,  &
                                  12, 12,  &
                                   2,  1,  &
                                  11,  1,  &
                                   1,  2], [7, 2], order=[2, 1])

        if (allocated(state%blade_xy)) deallocate(state%blade_xy)
        allocate(state%blade_xy(2, CRD_N_BLADES))

        blade_pitch = real(state%nx, wp) * state%dx / real(FUEL_BLADE_GRID_N, wp)
        x_origin    = 0.5_wp * blade_pitch
        y_origin    = 0.5_wp * blade_pitch

        b = 0
        dropped = 0
        do jg = 1, FUEL_BLADE_GRID_N
            do ig = 1, FUEL_BLADE_GRID_N
                ! Skip dropped positions until we have 144 − 7 = 137 left.
                skip = any(drop_positions(:, 1) == ig .and. drop_positions(:, 2) == jg)
                if (skip) then
                    dropped = dropped + 1
                    cycle
                end if
                b = b + 1
                if (b > CRD_N_BLADES) exit
                state%blade_xy(1, b) = x_origin + real(ig - 1, wp) * blade_pitch
                state%blade_xy(2, b) = y_origin + real(jg - 1, wp) * blade_pitch
            end do
            if (b >= CRD_N_BLADES) exit
        end do

        ! Silence -Wunused-parameter when running this on a config-only
        ! init (config presently has no blade-layout knobs).
        if (config%core_mask_factor < 0.0_wp) state%blade_xy(:, 1) = 0.0_wp
    end subroutine init_blade_xy

    subroutine init_blade_for_cell(state, config)
        ! For each cell (i, j) inside the core mask, find the nearest
        ! blade by Euclidean distance. Bypass-water cells get index 0.
        type(fuel_state_t),  intent(inout) :: state
        type(fuel_config_t), intent(in)    :: config

        integer  :: i, j, b, best_b
        real(wp) :: x, y, dx_, dy_, r_core, d2, best_d2, r_from_centre

        if (allocated(state%blade_for_cell)) deallocate(state%blade_for_cell)
        if (allocated(state%is_fuel_cell))  deallocate(state%is_fuel_cell)
        allocate(state%blade_for_cell(state%nx, state%ny))
        allocate(state%is_fuel_cell(state%nx, state%ny))

        r_core = state%core_diameter_m / config%core_mask_factor

        do j = 1, state%ny
            do i = 1, state%nx
                x = (real(i, wp) - 0.5_wp) * state%dx
                y = (real(j, wp) - 0.5_wp) * state%dy
                r_from_centre = sqrt((x - state%core_diameter_m / 2.0_wp)**2 &
                                   + (y - state%core_diameter_m / 2.0_wp)**2)
                state%is_fuel_cell(i, j) = (r_from_centre <= r_core)

                if (.not. state%is_fuel_cell(i, j)) then
                    state%blade_for_cell(i, j) = 0
                    cycle
                end if

                best_b  = 1
                best_d2 = huge(1.0_wp)
                do b = 1, CRD_N_BLADES
                    dx_ = x - state%blade_xy(1, b)
                    dy_ = y - state%blade_xy(2, b)
                    d2  = dx_ * dx_ + dy_ * dy_
                    if (d2 < best_d2) then
                        best_d2 = d2
                        best_b  = b
                    end if
                end do
                state%blade_for_cell(i, j) = best_b
            end do
        end do
    end subroutine init_blade_for_cell

    subroutine fuel_init_geometry(state, config, nx, ny, dx, dy, dz, core_diameter_m)
        ! Public init variant that takes the grid dimensions directly
        ! (used by tests and the orchestrator). The 5-call contract's
        ! `fuel_init` calls into here with sensible defaults.
        type(fuel_state_t),  intent(out) :: state
        type(fuel_config_t), intent(in)  :: config
        integer,             intent(in)  :: nx, ny
        real(wp),            intent(in)  :: dx, dy, dz, core_diameter_m

        state%time     = 0.0_wp
        state%n_steps  = 0
        state%nx       = nx
        state%ny       = ny
        state%dx       = dx
        state%dy       = dy
        state%dz       = dz
        state%core_diameter_m = core_diameter_m
        state%avg_fuel_temp_K = 0.0_wp
        state%max_fuel_temp_K = 0.0_wp

        call init_blade_xy(state, config)
        call init_blade_for_cell(state, config)
    end subroutine fuel_init_geometry

    subroutine fuel_init(state, config)
        ! Bare 5-call contract entry. Sets up a 1×1×1 default geometry
        ! (suitable for unit tests that never call into the kernels);
        ! the orchestrator should use fuel_init_geometry instead.
        type(fuel_state_t),  intent(out) :: state
        type(fuel_config_t), intent(in)  :: config
        call fuel_init_geometry(state, config, 1, 1, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp)
    end subroutine fuel_init

    subroutine fuel_setup_geometry(state, config, xslib, neutronics, nz)
        ! Materialise the initial fuel/water cross-section pattern in
        ! the multigroup-diffusion state. Cells inside the core radial
        ! mask get UO2_35; outside get H2O. Replaces the legacy
        ! `setup_geometry` in bwr_c_interface.
        type(fuel_state_t),    intent(in)    :: state
        type(fuel_config_t),   intent(in)    :: config
        type(xsec_library_t),  intent(in)    :: xslib
        type(mg_state_t),      intent(inout) :: neutronics
        integer,               intent(in)    :: nz

        integer :: i, j, k
        logical :: found_fuel, found_water
        type(xsec_material_t) :: xsec_fuel, xsec_water

        call xslib_get_material(xslib, "UO2_35", xsec_fuel,  found_fuel)
        if (.not. found_fuel) then
            print *, "fuel_setup_geometry: UO2_35 material not found"; stop 1
        end if
        call xslib_get_material(xslib, "H2O",    xsec_water, found_water)
        if (.not. found_water) then
            print *, "fuel_setup_geometry: H2O material not found"; stop 1
        end if

        do k = 1, nz
            do j = 1, state%ny
                do i = 1, state%nx
                    if (state%is_fuel_cell(i, j)) then
                        call mg_set_cross_sections(neutronics, xsec_fuel%xsec_base, &
                            i, i, j, j, k, k)
                    else
                        call mg_set_cross_sections(neutronics, xsec_water%xsec_base, &
                            i, i, j, j, k, k)
                    end if
                end do
            end do
        end do

        ! Touch config so -Wunused-dummy-argument stays quiet when this
        ! routine's parameter list grows but config isn't yet consumed.
        if (config%alpha_doppler == 0.0_wp .and. config%alpha_void == 0.0_wp) return
    end subroutine fuel_setup_geometry

    subroutine fuel_apply_to_xs(state, config, xslib, neutronics, &
                                burnup_field, xenon_field, samarium_field, &
                                blade_insertion, T, rho, avg_void_fraction, &
                                reactivity_perturbation_pcm)
        ! Per-cell XS feedback (Doppler, moderator, burnup, Xe, Sm) +
        ! per-blade control-rod fan-out. Replaces the legacy
        ! `update_cross_sections_feedback` in bwr_c_interface.
        type(fuel_state_t),    intent(inout) :: state
        type(fuel_config_t),   intent(in)    :: config
        type(xsec_library_t),  intent(in)    :: xslib
        type(mg_state_t),      intent(inout) :: neutronics
        real(wp),              intent(in)    :: burnup_field(:, :, :)
        real(wp),              intent(in)    :: xenon_field(:, :, :)
        real(wp),              intent(in)    :: samarium_field(:, :, :)
        real(wp),              intent(in)    :: blade_insertion(:)
        real(wp),              intent(in)    :: T(:, :, :)
        real(wp),              intent(in)    :: rho(:, :, :)
        real(wp),              intent(in)    :: avg_void_fraction   ! [0, 1]
        real(wp),              intent(in)    :: reactivity_perturbation_pcm

        integer  :: i, j, k, nz
        type(mg_xsec_t) :: xsec
        real(wp) :: T_fuel, rho_mod, rf, rf_ref, rho_mod_corrected
        real(wp) :: node_bottom, node_top, rod_tip, inserted_fraction, H_core
        real(wp) :: per_cell_insertion

        nz = size(T, 3)

        rf_ref = void_reactivity_factor(config%void_ref_fraction)
        rf     = void_reactivity_factor(avg_void_fraction)

        do k = 1, nz
            do j = 1, state%ny
                do i = 1, state%nx
                    T_fuel  = T(i, j, k)
                    rho_mod = rho(i, j, k)
                    rho_mod_corrected = rho_mod * (rf / max(1.0e-9_wp, rf_ref))

                    if (.not. state%is_fuel_cell(i, j)) then
                        ! Out-of-core cell: apply moderator XS, not fuel.
                        ! Matches fuel_setup_geometry's initial assignment.
                        call xslib_get_xsec(xslib, "H2O", T_fuel, rho_mod, 0.0_wp, xsec)
                        call mg_set_cross_sections(neutronics, xsec, i, i, j, j, k, k)
                        cycle
                    end if

                    call xslib_get_xsec(xslib, "UO2_35", &
                        T_fuel, rho_mod_corrected, burnup_field(i, j, k), xsec, &
                        Xe_conc = xenon_field(i, j, k), &
                        Sm_conc = samarium_field(i, j, k))

                    ! Per-blade fan-out: each cell asks its nearest blade
                    ! for the rod insertion.
                    per_cell_insertion = blade_insertion(state%blade_for_cell(i, j))

                    ! Axial overlap of the inserted-rod tip with this
                    ! node. Bottom-entry rods: rod_tip rises with
                    ! `per_cell_insertion ∈ [0, 1]`.
                    H_core      = real(nz, wp) * state%dz
                    node_bottom = real(k - 1, wp) * state%dz
                    node_top    = real(k,     wp) * state%dz
                    rod_tip     = per_cell_insertion * H_core

                    inserted_fraction = 0.0_wp
                    if (rod_tip >= node_top) then
                        inserted_fraction = 1.0_wp
                    else if (rod_tip > node_bottom) then
                        inserted_fraction = (rod_tip - node_bottom) / state%dz
                    end if

                    if (inserted_fraction > 0.0_wp) &
                        call xslib_apply_control_rod(xsec, inserted_fraction)

                    call mg_set_cross_sections(neutronics, xsec, i, i, j, j, k, k)
                end do
            end do
        end do

        ! Operator-driven sustained reactivity perturbation. Mutates
        ! nu_sigma_f in place each tick; persists until cleared.
        call fuel_apply_reactivity_perturbation(neutronics, reactivity_perturbation_pcm)

        ! Cache instantaneous fuel-temp scalars for observability.
        state%avg_fuel_temp_K = sum(T) / real(size(T), wp)
        state%max_fuel_temp_K = maxval(T)
    end subroutine fuel_apply_to_xs

    subroutine fuel_apply_reactivity_perturbation(neutronics, perturbation_pcm)
        ! Scale ν·Σ_f by (1 + perturbation/1e5) — applied AFTER the
        ! per-cell XS rebuild so the perturbation persists tick over tick.
        type(mg_state_t), intent(inout) :: neutronics
        real(wp),         intent(in)    :: perturbation_pcm

        integer  :: i, j, k, g
        real(wp) :: factor

        if (abs(perturbation_pcm) < 1.0e-12_wp) return
        factor = 1.0_wp + perturbation_pcm / 1.0e5_wp
        do k = 1, neutronics%nz
            do j = 1, neutronics%ny
                do i = 1, neutronics%nx
                    do g = 1, neutronics%n_groups
                        neutronics%xsec(i,j,k)%nu_sigma_f(g) = &
                            neutronics%xsec(i,j,k)%nu_sigma_f(g) * factor
                    end do
                end do
            end do
        end do
    end subroutine fuel_apply_reactivity_perturbation

    subroutine fuel_step(state, dt)
        type(fuel_state_t), intent(inout) :: state
        real(wp),           intent(in)    :: dt
        state%time    = state%time + dt
        state%n_steps = state%n_steps + 1
    end subroutine fuel_step

    subroutine fuel_apply(state, command)
        type(fuel_state_t),   intent(inout) :: state
        type(fuel_command_t), intent(in)    :: command
        ! No fuel-level commands in phase 1. Touch both args so -Wextra
        ! doesn't flag the placeholder.
        if (command%reserved /= 0) state%n_steps = state%n_steps
    end subroutine fuel_apply

    function fuel_observe(state) result(obs)
        type(fuel_state_t), intent(in) :: state
        type(fuel_observation_t) :: obs
        obs%time            = state%time
        obs%n_steps         = state%n_steps
        obs%avg_fuel_temp_K = state%avg_fuel_temp_K
        obs%max_fuel_temp_K = state%max_fuel_temp_K
    end function fuel_observe

    subroutine fuel_destroy(state)
        type(fuel_state_t), intent(inout) :: state
        state%time     = 0.0_wp
        state%n_steps  = 0
        state%avg_fuel_temp_K = 0.0_wp
        state%max_fuel_temp_K = 0.0_wp
        if (allocated(state%blade_for_cell)) deallocate(state%blade_for_cell)
        if (allocated(state%blade_xy))       deallocate(state%blade_xy)
        if (allocated(state%is_fuel_cell))   deallocate(state%is_fuel_cell)
    end subroutine fuel_destroy

end module fuel