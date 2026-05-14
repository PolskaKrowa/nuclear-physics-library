! bwr_core_simulator_opengl.f90

program bwr_core_simulator_opengl
    use kinds, only: wp
    use iso_c_binding

    ! Physics modules
    use multigroup_diffusion
    use cross_sections
    use burnup_depletion
    use heat_transfer
    use two_phase_flow
    use fluid_dynamics
    use pressure_dynamics
    use rng, only: rng_seed

    implicit none

    ! ── OpenGL / GLUT interfaces ────────────────────────────────
    interface
        subroutine glutInit(argcp, argv) bind(C, name='glutInit')
            import :: c_ptr, c_int
            type(c_ptr), value :: argcp, argv
        end subroutine
        subroutine glutInitDisplayMode(mode) bind(C, name='glutInitDisplayMode')
            import :: c_int
            integer(c_int), value :: mode
        end subroutine
        subroutine glutInitWindowSize(w, h) bind(C, name='glutInitWindowSize')
            import :: c_int
            integer(c_int), value :: w, h
        end subroutine
        function glutCreateWindow(title) bind(C, name='glutCreateWindow')
            import :: c_int, c_char
            character(kind=c_char) :: title(*)
            integer(c_int) :: glutCreateWindow
        end function
        subroutine glutDisplayFunc(func) bind(C, name='glutDisplayFunc')
            import :: c_funptr
            type(c_funptr), value :: func
        end subroutine
        subroutine glutReshapeFunc(func) bind(C, name='glutReshapeFunc')
            import :: c_funptr
            type(c_funptr), value :: func
        end subroutine
        subroutine glutKeyboardFunc(func) bind(C, name='glutKeyboardFunc')
            import :: c_funptr
            type(c_funptr), value :: func
        end subroutine
        subroutine glutIdleFunc(func) bind(C, name='glutIdleFunc')
            import :: c_funptr
            type(c_funptr), value :: func
        end subroutine
        subroutine glutMainLoop() bind(C, name='glutMainLoop')
        end subroutine
        subroutine glutPostRedisplay() bind(C, name='glutPostRedisplay')
        end subroutine
        subroutine glutSwapBuffers() bind(C, name='glutSwapBuffers')
        end subroutine
        subroutine glClear(mask) bind(C, name='glClear')
            import :: c_int
            integer(c_int), value :: mask
        end subroutine
        subroutine glClearColor(r, g, b, a) bind(C, name='glClearColor')
            import :: c_float
            real(c_float), value :: r, g, b, a
        end subroutine
        subroutine glEnable(cap) bind(C, name='glEnable')
            import :: c_int
            integer(c_int), value :: cap
        end subroutine
        subroutine glLoadIdentity() bind(C, name='glLoadIdentity')
        end subroutine
        subroutine glMatrixMode(mode) bind(C, name='glMatrixMode')
            import :: c_int
            integer(c_int), value :: mode
        end subroutine
        subroutine glViewport(x, y, w, h) bind(C, name='glViewport')
            import :: c_int
            integer(c_int), value :: x, y, w, h
        end subroutine
        subroutine gluPerspective(fovy, aspect, zNear, zFar) bind(C, name='gluPerspective')
            import :: c_double
            real(c_double), value :: fovy, aspect, zNear, zFar
        end subroutine
        subroutine glTranslatef(x, y, z) bind(C, name='glTranslatef')
            import :: c_float
            real(c_float), value :: x, y, z
        end subroutine
        subroutine glRotatef(angle, x, y, z) bind(C, name='glRotatef')
            import :: c_float
            real(c_float), value :: angle, x, y, z
        end subroutine
        subroutine glBegin(mode) bind(C, name='glBegin')
            import :: c_int
            integer(c_int), value :: mode
        end subroutine
        subroutine glEnd() bind(C, name='glEnd')
        end subroutine
        subroutine glVertex3f(x, y, z) bind(C, name='glVertex3f')
            import :: c_float
            real(c_float), value :: x, y, z
        end subroutine
        subroutine glColor3f(r, g, b) bind(C, name='glColor3f')
            import :: c_float
            real(c_float), value :: r, g, b
        end subroutine
        subroutine glRasterPos2f(x, y) bind(C, name='glRasterPos2f')
            import :: c_float
            real(c_float), value :: x, y
        end subroutine
        subroutine glutBitmapCharacter(font, character) bind(C, name='glutBitmapCharacter')
            import :: c_ptr, c_int
            type(c_ptr), value :: font
            integer(c_int), value :: character
        end subroutine
        subroutine glBlendFunc(sfactor, dfactor) bind(C, name='glBlendFunc')
            import :: c_int
            integer(c_int), value :: sfactor, dfactor
        end subroutine
        subroutine glColor4f(r, g, b, a) bind(C, name='glColor4f')
            import :: c_float
            real(c_float), value :: r, g, b, a
        end subroutine
    end interface

    ! ── OpenGL constants ────────────────────────────────────────
    integer(c_int), parameter :: GLUT_RGB               = int(z'0000', c_int)
    integer(c_int), parameter :: GLUT_DOUBLE            = int(z'0002', c_int)
    integer(c_int), parameter :: GLUT_DEPTH             = int(z'0010', c_int)
    integer(c_int), parameter :: GL_COLOR_BUFFER_BIT    = int(z'00004000', c_int)
    integer(c_int), parameter :: GL_DEPTH_BUFFER_BIT    = int(z'00000100', c_int)
    integer(c_int), parameter :: GL_DEPTH_TEST          = int(z'0B71', c_int)
    integer(c_int), parameter :: GL_PROJECTION          = int(z'1701', c_int)
    integer(c_int), parameter :: GL_MODELVIEW           = int(z'1700', c_int)
    integer(c_int), parameter :: GL_QUADS               = int(z'0007', c_int)
    integer(c_int), parameter :: GL_LINES               = int(z'0001', c_int)
    integer(c_int), parameter :: GL_BLEND               = int(z'0BE2', c_int)
    integer(c_int), parameter :: GL_SRC_ALPHA           = int(z'0302', c_int)
    integer(c_int), parameter :: GL_ONE_MINUS_SRC_ALPHA = int(z'0303', c_int)

    ! ── Simulation state ────────────────────────────────────────
    type :: simulation_t
        ! Grid
        integer  :: nx, ny, nz
        real(wp) :: dx, dy, dz
        real(wp) :: core_height, core_diameter
        real(wp) :: rod_bank_position  ! 0.0 (All Out) to 1.0 (All In) – BWR bottom entry

        ! Physics modules
        type(mg_state_t)         :: neutronics
        type(heat_state_t)       :: heat
        type(two_phase_state_t)  :: thermalhydraulics
        type(burnup_state_t)     :: burnup
        type(xsec_library_t)     :: xslib

        ! Operating conditions
        real(wp) :: power_rated, power_current
        real(wp) :: pressure_operating    ! Now dynamic [Pa]
        real(wp) :: mass_flux_core
        real(wp) :: inlet_temperature

        ! Core parameters
        real(wp) :: k_eff, reactivity_pcm
        real(wp) :: max_fuel_temp, max_clad_temp
        real(wp) :: avg_void_fraction, max_void_fraction
        real(wp) :: min_chfr, avg_burnup

        ! Feedback coefficients
        real(wp) :: alpha_doppler, alpha_void

        ! Time
        real(wp) :: time
        integer  :: n_steps

        ! Initialisation flag
        logical :: initialized

        ! Saturation and level
        real(wp) :: sat_temperature    ! Saturation temperature at Prx [K]
        real(wp) :: Lbase              ! Effective water inventory [m equivalent]
        real(wp) :: Lrx                ! Water level above core floor [m]

        ! 3-point smoothing buffers for coolant temperature
        real(wp) :: avg_coolant_temp   ! Spatial average of heat%T [K]
        real(wp) :: coolant_T_prev1    ! Previous avg coolant temp
        real(wp) :: coolant_T_prev2

        ! 3-point smoothing for pressure
        real(wp) :: P_prev1
        real(wp) :: P_prev2

        ! Flash steam excess above saturation
        real(wp) :: T_excess

        ! Reactor period
        real(wp) :: neutrons_prev      ! Previous power level [W] for period calc
        real(wp) :: reactor_period     ! Reactor period [s]; +∞ = stable

        ! APRM instrumentation
        real(wp) :: aprm               ! Average Power Range Monitor [% rated]

        ! Void reactivity factor
        real(wp) :: rf_void

        ! Turbine
        real(wp) :: turbine_valve      ! Steam admission valve [0-100 %]
        real(wp) :: bypass_valve       ! Bypass valve [0-100 %]
        real(wp) :: steam_flow         ! Normalised total steam flow
        real(wp) :: turbine_speed      ! Rotational speed [rpm]
        real(wp) :: turbine_power      ! Net electrical output [W]
        logical  :: turbine_tripped    ! Turbine trip / isolation flag

        ! Feedwater temperature
        real(wp) :: feedwater_temp     ! [K]
    end type simulation_t

    ! ── Visualisation state ─────────────────────────────────────
    type :: viz_state_t
        integer  :: display_mode  ! 1=power,2=temp,3=void,4=burnup,5=T/T_sat ratio
        real(wp) :: rotation_x, rotation_y
        real(wp) :: zoom
        logical  :: paused
        logical  :: show_help
        integer  :: update_counter
    end type viz_state_t

    ! ── Global state (must be SAVE for OpenGL callbacks) ────────
    type(simulation_t), target, save :: g_sim
    type(viz_state_t),          save :: g_viz
    real(wp), save :: g_dt    = 0.01_wp
    logical,  save :: g_converged = .false.
    integer(c_int), save :: g_window

    ! ── Main program ─────────────────────────────────────────────
    print *, "=================================================="
    print *, "  NPL Example BWR Simulation Program"
    print *, "=================================================="
    print *, ""

    call rng_seed(123456789_8)

    print *, "Setting up simulation..."
    call setup_bwr_simulation(g_sim)
    g_sim%initialized = .true.

    print *, "Computing initial steady-state..."
    call solve_steady_state(g_sim)
    call print_steady_state_summary(g_sim)

    ! Initialise visualisation state
    g_viz%display_mode  = 1
    g_viz%rotation_x    = 30.0_wp
    g_viz%rotation_y    = 45.0_wp
    g_viz%zoom          = -15.0_wp
    g_viz%paused        = .false.
    g_viz%show_help     = .true.
    g_viz%update_counter= 0

    print *, ""
    print *, "=== CONTROLS ==="
    print *, "  SPACE    - Pause / Resume"
    print *, "  1-5      - Display (Power/Temp/Void/Burnup/T_sat ratio)"
    print *, "  +/-      - Increase / Decrease power setpoint"
    print *, "  r/R      - Insert +/- 20 pcm reactivity"
    print *, "  c/C      - Increase / Decrease coolant flow"
    print *, "  i/o      - Insert / Withdraw control rods"
    print *, "  t/T      - Open / Close turbine valve"
    print *, "  b/B      - Open / Close bypass valve"
    print *, "  w/s      - Rotate up / down"
    print *, "  a/d      - Rotate left / right"
    print *, "  z/x      - Zoom in / out"
    print *, "  h        - Toggle help"
    print *, "  q / ESC  - Quit"
    print *, "================"
    print *, ""
    print *, "Starting OpenGL visualisation..."

    call init_opengl()
    call glutMainLoop()

    call cleanup_simulation(g_sim)

contains

    ! ================================================================
    !  PURE PHYSICS FUNCTIONS
    ! ================================================================

    !> Saturation temperature as a function of pressure.
    !  Polynomial fit to steam tables.
    !  Returns temperature in Kelvin.
    pure function sat_temp_K(P_Pa) result(T_K)
        real(wp), intent(in) :: P_Pa
        real(wp) :: T_K, p
        p   = P_Pa / 1.0e6_wp   ! convert to MPa
        T_K = (-4.0964e-3_wp * p**6  &
              +  0.141738_wp  * p**5  &
              -  1.943057_wp  * p**4  &
              + 13.50428_wp   * p**3  &
              - 51.27379_wp   * p**2  &
              + 118.6854_wp   * p     &
              + 100.1542_wp) + 273.15_wp   ! °C → K
    end function sat_temp_K

    !> Liquid water density polynomial
    pure function coolant_density_gl(T_K) result(rho)
        real(wp), intent(in) :: T_K
        real(wp) :: rho, T_C
        T_C = T_K - 273.15_wp
        rho = -4.467711e-6_wp   * T_C**3  &
              - 5.60288485e-4_wp * T_C**2  &
              - 0.429148844451_wp* T_C     &
              + 1010.035413387815_wp
        rho = max(0.1_wp, rho)  ! guard against out-of-range temperatures
    end function coolant_density_gl

    !> Reactor water level above core floor
    pure function calc_reactor_level(Lbase, T_K, H_core) result(Lrx)
        real(wp), intent(in) :: Lbase, T_K, H_core
        real(wp) :: Lrx, rho
        rho = coolant_density_gl(T_K)
        Lrx = (Lbase / (rho / 1000.0_wp)) - H_core
    end function calc_reactor_level

    !> Void reactivity factor from effective moderator density.
    pure function void_reactivity_factor(avg_void_frac) result(rf)
        real(wp), intent(in) :: avg_void_frac
        real(wp) :: rf, dens2
        ! Effective moderator density: 0 (full void) → 1000 kg/m³ (full liquid)
        dens2 = (1.0_wp - min(1.0_wp, max(0.0_wp, avg_void_frac))) * 1000.0_wp
        rf = -5.84e-10_wp   * dens2**3  &
             + 1.35422e-7_wp * dens2**2  &
             + 1.358252042e-3_wp * dens2  &
             + 0.090469568418_wp
    end function void_reactivity_factor

    ! ================================================================
    !  OPENGL INITIALISATION
    ! ================================================================

    subroutine init_opengl()
        character(kind=c_char, len=25) :: title
        integer(c_int), target :: argc = 0
        type(c_ptr),    target :: argv = c_null_ptr

        title = "BWR Core Simulator"//c_null_char

        call glutInit(c_loc(argc), c_loc(argv))
        call glutInitDisplayMode(ior(ior(GLUT_RGB, GLUT_DOUBLE), GLUT_DEPTH))
        call glutInitWindowSize(1200_c_int, 900_c_int)
        g_window = glutCreateWindow(title)

        call glutDisplayFunc(c_funloc(display_callback))
        call glutReshapeFunc(c_funloc(reshape_callback))
        call glutKeyboardFunc(c_funloc(keyboard_callback))
        call glutIdleFunc(c_funloc(idle_callback))

        call glEnable(GL_BLEND)
        call glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
        call glClearColor(0.1_c_float, 0.1_c_float, 0.15_c_float, 1.0_c_float)
        call glEnable(GL_DEPTH_TEST)
    end subroutine init_opengl

    ! ================================================================
    !  DISPLAY CALLBACK
    ! ================================================================

    subroutine display_callback() bind(C)
        integer  :: i, j, k, skip
        real(wp) :: x, y, z, value, r, g, b
        real(wp) :: val_min, val_max

        if (.not. g_sim%initialized) return

        call glClear(ior(GL_COLOR_BUFFER_BIT, GL_DEPTH_BUFFER_BIT))
        call glMatrixMode(GL_MODELVIEW)
        call glLoadIdentity()

        call glTranslatef(0.0_c_float, 0.0_c_float, real(g_viz%zoom, c_float))
        call glRotatef(real(g_viz%rotation_x, c_float), 1.0_c_float, 0.0_c_float, 0.0_c_float)
        call glRotatef(real(g_viz%rotation_y, c_float), 0.0_c_float, 1.0_c_float, 0.0_c_float)

        call draw_axes()

        ! Determine colour map range for active display mode
        val_min = 0.0_wp; val_max = 1.0_wp

        select case (g_viz%display_mode)
        case (1)
            if (allocated(g_sim%neutronics%power_density)) then
                val_min = minval(g_sim%neutronics%power_density)
                val_max = maxval(g_sim%neutronics%power_density)
            end if
        case (2)
            if (allocated(g_sim%heat%T)) then
                val_min = minval(g_sim%heat%T) - 273.15_wp
                val_max = maxval(g_sim%heat%T) - 273.15_wp
            end if
        case (3)
            val_min = 0.0_wp;  val_max = 0.8_wp
        case (4)
            val_min = 0.0_wp
            if (allocated(g_sim%burnup%burnup)) &
                val_max = max(maxval(g_sim%burnup%burnup), 1.0_wp)
        case (5)
            ! T / T_sat ratio: 1.0 = at saturation, >1 = superheated steam
            val_min = 0.5_wp;  val_max = 1.1_wp
        end select

        skip = max(1, g_sim%nx / 15)

        do k = 1, g_sim%nz, skip
            do j = 1, g_sim%ny, skip
                do i = 1, g_sim%nx, skip
                    value = 0.0_wp

                    select case (g_viz%display_mode)
                    case (1)
                        if (allocated(g_sim%neutronics%power_density)) &
                            value = g_sim%neutronics%power_density(i, j, k)
                    case (2)
                        if (allocated(g_sim%heat%T)) &
                            value = g_sim%heat%T(i, j, k) - 273.15_wp
                    case (3)
                        if (allocated(g_sim%thermalhydraulics%void_fraction)) &
                            value = g_sim%thermalhydraulics%void_fraction(i, j, k)
                    case (4)
                        if (allocated(g_sim%burnup%burnup)) &
                            value = g_sim%burnup%burnup(i, j, k)
                    case (5)
                        ! T/T_sat ratio per cell
                        if (allocated(g_sim%heat%T)) &
                            value = g_sim%heat%T(i, j, k) / max(1.0_wp, g_sim%sat_temperature)
                    end select

                    call value_to_colour(value, val_min, val_max, r, g, b)

                    x = (real(i, wp) - real(g_sim%nx, wp) / 2.0_wp) * g_sim%dx
                    y = (real(j, wp) - real(g_sim%ny, wp) / 2.0_wp) * g_sim%dy
                    z = (real(k, wp) - real(g_sim%nz, wp) / 2.0_wp) * g_sim%dz

                    call draw_cube(x, y, z, g_sim%dx * 0.4_wp * skip, r, g, b, 0.2_wp)
                end do
            end do
        end do

        ! Draw water level indicator
        call draw_water_level()

        call glutSwapBuffers()
    end subroutine display_callback

    ! ================================================================
    !  RESHAPE CALLBACK
    ! ================================================================

    subroutine reshape_callback(width, height) bind(C)
        integer(c_int), value :: width, height
        real(c_double) :: aspect

        call glViewport(0_c_int, 0_c_int, width, height)
        call glMatrixMode(GL_PROJECTION)
        call glLoadIdentity()
        aspect = real(width, c_double) / real(max(height, 1), c_double)
        call gluPerspective(45.0_c_double, aspect, 0.1_c_double, 100.0_c_double)
        call glMatrixMode(GL_MODELVIEW)
    end subroutine reshape_callback

    ! ================================================================
    !  KEYBOARD CALLBACK
    ! ================================================================

    subroutine keyboard_callback(key, x, y) bind(C)
        integer(c_int), value :: key, x, y

        select case (key)
        case (32)  ! SPACE – pause / resume
            g_viz%paused = .not. g_viz%paused
            print *, merge("PAUSED  ", "RESUMED ", g_viz%paused)

        ! Display modes
        case (49); g_viz%display_mode = 1; print *, "Display: Power Distribution"
        case (50); g_viz%display_mode = 2; print *, "Display: Temperature"
        case (51); g_viz%display_mode = 3; print *, "Display: Void Fraction"
        case (52); g_viz%display_mode = 4; print *, "Display: Burnup"
        case (53); g_viz%display_mode = 5; print *, "Display: T / T_sat Ratio"

        ! Power setpoint
        case (43, 61)  ! '+' or '='
            g_sim%power_rated = g_sim%power_rated * 1.05_wp
            print '(A,F8.1,A)', " Power setpoint: ", g_sim%power_rated/1.0e6_wp, " MW"
        case (45)      ! '-'
            g_sim%power_rated = g_sim%power_rated * 0.95_wp
            print '(A,F8.1,A)', " Power setpoint: ", g_sim%power_rated/1.0e6_wp, " MW"

        ! Reactivity insertions
        case (114)  ! 'r'
            call apply_reactivity_insertion(g_sim, 20.0_wp)
            print *, " +20 pcm reactivity inserted"
        case (82)   ! 'R'
            call apply_reactivity_insertion(g_sim, -20.0_wp)
            print *, " -20 pcm reactivity inserted"

        ! Coolant flow
        case (99)   ! 'c'
            g_sim%mass_flux_core = g_sim%mass_flux_core * 1.05_wp
            print '(A,F7.1,A)', " Flow: ", g_sim%mass_flux_core, " kg/m²·s"
        case (67)   ! 'C'
            g_sim%mass_flux_core = g_sim%mass_flux_core * 0.95_wp
            print '(A,F7.1,A)', " Flow: ", g_sim%mass_flux_core, " kg/m²·s"

        ! Control rods (BWR bottom-entry, 'i' inserts from bottom)
        case (ichar('i'), ichar('I'))
            g_sim%rod_bank_position = min(1.0_wp, g_sim%rod_bank_position + 0.005_wp)
            print '(A,F6.1,A)', " Rods inserting: ", g_sim%rod_bank_position * 100.0_wp, " %"
        case (ichar('o'), ichar('O'))
            g_sim%rod_bank_position = max(0.0_wp, g_sim%rod_bank_position - 0.005_wp)
            print '(A,F6.1,A)', " Rods withdrawing: ", g_sim%rod_bank_position * 100.0_wp, " %"

        ! ── Turbine valve ────────────────
        case (116)  ! 't' – open turbine valve
            g_sim%turbine_valve = min(100.0_wp, g_sim%turbine_valve + 2.0_wp)
            print '(A,F5.1,A)', " Turbine valve: ", g_sim%turbine_valve, " %"
        case (84)   ! 'T' – close turbine valve
            g_sim%turbine_valve = max(0.0_wp, g_sim%turbine_valve - 2.0_wp)
            print '(A,F5.1,A)', " Turbine valve: ", g_sim%turbine_valve, " %"
        case (98)   ! 'b' – open bypass valve
            g_sim%bypass_valve = min(100.0_wp, g_sim%bypass_valve + 2.0_wp)
            print '(A,F5.1,A)', " Bypass valve:  ", g_sim%bypass_valve, " %"
        case (66)   ! 'B' – close bypass valve
            g_sim%bypass_valve = max(0.0_wp, g_sim%bypass_valve - 2.0_wp)
            print '(A,F5.1,A)', " Bypass valve:  ", g_sim%bypass_valve, " %"

        ! Visualisation controls
        case (119); g_viz%rotation_x = g_viz%rotation_x + 5.0_wp  ! w
        case (115); g_viz%rotation_x = g_viz%rotation_x - 5.0_wp  ! s
        case (97);  g_viz%rotation_y = g_viz%rotation_y - 5.0_wp  ! a
        case (100); g_viz%rotation_y = g_viz%rotation_y + 5.0_wp  ! d
        case (122); g_viz%zoom = g_viz%zoom + 1.0_wp               ! z
        case (120); g_viz%zoom = g_viz%zoom - 1.0_wp               ! x
        case (104); g_viz%show_help = .not. g_viz%show_help         ! h

        ! Quit
        case (27, 113)  ! ESC or 'q'
            print *, "Shutting down..."
            call cleanup_simulation(g_sim)
            stop
        end select

        call glutPostRedisplay()
    end subroutine keyboard_callback

    ! ================================================================
    !  IDLE CALLBACK  (physics update)
    ! ================================================================

    subroutine idle_callback() bind(C)
        if (.not. g_viz%paused .and. g_sim%initialized) then
            call coupled_time_step(g_sim, g_dt)
            call apply_automatic_controls(g_sim, g_dt)

            g_viz%update_counter = g_viz%update_counter + 1

            if (mod(g_viz%update_counter, 100) == 0) then
                call print_transient_summary(g_sim, g_sim%n_steps, g_sim%time)
            end if
        end if

        call glutPostRedisplay()
    end subroutine idle_callback

    ! ================================================================
    !  RENDERING HELPERS
    ! ================================================================

    subroutine draw_cube(x, y, z, size, r, g, b, alpha)
        real(wp), intent(in) :: x, y, z, size, r, g, b, alpha
        real(c_float) :: xf, yf, zf, sf

        xf = real(x, c_float); yf = real(y, c_float)
        zf = real(z, c_float); sf = real(size, c_float)

        call glColor4f(real(r, c_float), real(g, c_float), real(b, c_float), real(alpha, c_float))
        call glBegin(GL_QUADS)
        ! Front
        call glVertex3f(xf-sf, yf-sf, zf+sf); call glVertex3f(xf+sf, yf-sf, zf+sf)
        call glVertex3f(xf+sf, yf+sf, zf+sf); call glVertex3f(xf-sf, yf+sf, zf+sf)
        ! Back
        call glVertex3f(xf-sf, yf-sf, zf-sf); call glVertex3f(xf-sf, yf+sf, zf-sf)
        call glVertex3f(xf+sf, yf+sf, zf-sf); call glVertex3f(xf+sf, yf-sf, zf-sf)
        ! Top
        call glVertex3f(xf-sf, yf+sf, zf-sf); call glVertex3f(xf-sf, yf+sf, zf+sf)
        call glVertex3f(xf+sf, yf+sf, zf+sf); call glVertex3f(xf+sf, yf+sf, zf-sf)
        ! Bottom
        call glVertex3f(xf-sf, yf-sf, zf-sf); call glVertex3f(xf+sf, yf-sf, zf-sf)
        call glVertex3f(xf+sf, yf-sf, zf+sf); call glVertex3f(xf-sf, yf-sf, zf+sf)
        ! Right
        call glVertex3f(xf+sf, yf-sf, zf-sf); call glVertex3f(xf+sf, yf+sf, zf-sf)
        call glVertex3f(xf+sf, yf+sf, zf+sf); call glVertex3f(xf+sf, yf-sf, zf+sf)
        ! Left
        call glVertex3f(xf-sf, yf-sf, zf-sf); call glVertex3f(xf-sf, yf-sf, zf+sf)
        call glVertex3f(xf-sf, yf+sf, zf+sf); call glVertex3f(xf-sf, yf+sf, zf-sf)
        call glEnd()
    end subroutine draw_cube

    subroutine draw_axes()
        real(c_float), parameter :: len = 3.0_c_float
        call glBegin(GL_LINES)
        call glColor3f(1.0_c_float, 0.0_c_float, 0.0_c_float)  ! X – red
        call glVertex3f(0.0_c_float, 0.0_c_float, 0.0_c_float)
        call glVertex3f(len, 0.0_c_float, 0.0_c_float)
        call glColor3f(0.0_c_float, 1.0_c_float, 0.0_c_float)  ! Y – green
        call glVertex3f(0.0_c_float, 0.0_c_float, 0.0_c_float)
        call glVertex3f(0.0_c_float, len, 0.0_c_float)
        call glColor3f(0.0_c_float, 0.0_c_float, 1.0_c_float)  ! Z – blue
        call glVertex3f(0.0_c_float, 0.0_c_float, 0.0_c_float)
        call glVertex3f(0.0_c_float, 0.0_c_float, len)
        call glEnd()
    end subroutine draw_axes

    !> Draw a translucent horizontal quad at the reactor water level.
    subroutine draw_water_level()
        real(c_float) :: R, half_d, level_z
        real(wp) :: frac

        if (.not. g_sim%initialized) return

        ! Lrx is the level above the core floor; map to scene z-coordinate
        ! Core extends from -core_height/2 to +core_height/2 in scene space.
        frac   = g_sim%Lrx / max(0.01_wp, g_sim%core_height + 2.0_wp)
        level_z = real(-g_sim%core_height / 2.0_wp + g_sim%Lrx, c_float)
        half_d  = real(g_sim%core_diameter * 0.6_wp, c_float)

        ! Semi-transparent cyan surface for water level
        call glColor4f(0.0_c_float, 0.7_c_float, 1.0_c_float, 0.25_c_float)
        call glBegin(GL_QUADS)
        call glVertex3f(-half_d, -half_d, level_z)
        call glVertex3f( half_d, -half_d, level_z)
        call glVertex3f( half_d,  half_d, level_z)
        call glVertex3f(-half_d,  half_d, level_z)
        call glEnd()
    end subroutine draw_water_level

    !> Convert a scalar value to a heat-map colour (blue→cyan→green→yellow→red).
    subroutine value_to_colour(value, vmin, vmax, r, g, b)
        real(wp), intent(in)  :: value, vmin, vmax
        real(wp), intent(out) :: r, g, b
        real(wp) :: t

        if (abs(vmax - vmin) > 1.0e-10_wp) then
            t = (value - vmin) / (vmax - vmin)
        else
            t = 0.5_wp
        end if
        t = max(0.0_wp, min(1.0_wp, t))

        if      (t < 0.25_wp) then
            r = 0.0_wp; g = 4.0_wp*t;           b = 1.0_wp
        else if (t < 0.50_wp) then
            r = 0.0_wp; g = 1.0_wp;              b = 1.0_wp - 4.0_wp*(t-0.25_wp)
        else if (t < 0.75_wp) then
            r = 4.0_wp*(t-0.5_wp); g = 1.0_wp;  b = 0.0_wp
        else
            r = 1.0_wp; g = 1.0_wp - 4.0_wp*(t-0.75_wp); b = 0.0_wp
        end if
    end subroutine value_to_colour

    ! ================================================================
    !  GEOMETRY SETUP
    ! ================================================================

    subroutine setup_bwr_geometry_realistic(sim)
        type(simulation_t), intent(inout) :: sim

        integer  :: i, j, k
        real(wp) :: x, y, r, x_center, y_center, pitch, pin_radius
        integer  :: pin_i, pin_j
        real(wp) :: dx_from_pin, dy_from_pin, r_from_pin
        logical  :: in_fuel, found_water, found_fuel

        type(xsec_material_t) :: xsec_fuel, xsec_water

        call xslib_get_material(sim%xslib, "UO2_35", xsec_fuel, found_fuel)
        if (.not. found_fuel) then
            print *, "Fatal: material UO2_35 not found."
            call xslib_list_materials(sim%xslib); stop 1
        end if
        call xslib_get_material(sim%xslib, "H2O", xsec_water, found_water)
        if (.not. found_water) then
            print *, "Fatal: material H2O not found."
            call xslib_list_materials(sim%xslib); stop 1
        end if

        pitch      = 0.0163_wp  ! 1.63 cm pin pitch (BWR typical)
        pin_radius = 0.0052_wp  ! 5.2 mm fuel pellet radius

        print *, "BWR lattice parameters:"
        print '(A,F7.4,A)', "  Pin pitch:         ", pitch*100, " cm"
        print '(A,F7.4,A)', "  Pin radius:        ", pin_radius*100, " cm"
        print '(A,F7.2)',   "  Fuel-to-mod ratio: ", &
            (3.14159_wp*pin_radius**2) / (pitch**2 - 3.14159_wp*pin_radius**2)

        do k = 1, sim%nz
            do j = 1, sim%ny
                do i = 1, sim%nx
                    x = (real(i, wp) - 0.5_wp) * sim%dx
                    y = (real(j, wp) - 0.5_wp) * sim%dy

                    pin_i = nint(x / pitch);  pin_j = nint(y / pitch)
                    x_center = real(pin_i, wp) * pitch
                    y_center = real(pin_j, wp) * pitch

                    dx_from_pin = x - x_center
                    dy_from_pin = y - y_center
                    r_from_pin  = sqrt(dx_from_pin**2 + dy_from_pin**2)

                    in_fuel = (r_from_pin < pin_radius)

                    ! Cylindrical core boundary
                    r = sqrt((x - sim%core_diameter/2)**2 + (y - sim%core_diameter/2)**2)
                    if (r > sim%core_diameter / 2.2_wp) in_fuel = .false.

                    if (in_fuel) then
                        call mg_set_cross_sections(sim%neutronics, xsec_fuel%xsec_base, &
                            i, i, j, j, k, k)
                    else
                        call mg_set_cross_sections(sim%neutronics, xsec_water%xsec_base, &
                            i, i, j, j, k, k)
                    end if
                end do
            end do
        end do
        print *, "Pin-by-pin geometry set up successfully"
    end subroutine setup_bwr_geometry_realistic

    ! ================================================================
    !  SIMULATION SETUP
    ! ================================================================

    subroutine setup_bwr_simulation(sim)
        type(simulation_t), intent(out) :: sim

        type(mg_config_t)        :: neutron_config
        type(heat_config_t)      :: heat_config
        type(two_phase_config_t) :: th_config
        type(burnup_config_t)    :: burnup_config
        type(xsec_material_t)    :: xsec_fuel, xsec_water

        ! ── Core geometry ─────────────────────────────────────────
        sim%nx = 20;  sim%ny = 20;  sim%nz = 20
        sim%core_height   = 3.71_wp   ! m
        sim%core_diameter = 4.75_wp   ! m
        sim%dx = sim%core_diameter / real(sim%nx, wp)
        sim%dy = sim%core_diameter / real(sim%ny, wp)
        sim%dz = sim%core_height   / real(sim%nz, wp)

        ! ── Operating conditions ──────────────────────────────────
        sim%power_rated         = 2381.0e6_wp  ! W thermal
        sim%pressure_operating  = 7.14e6_wp    ! Pa (now dynamic)
        sim%mass_flux_core      = 1500.0_wp    ! kg/m²·s
        sim%inlet_temperature   = 551.0_wp     ! K

        ! ── Feedback coefficients ─────────────────────────────────
        sim%alpha_doppler = -3.5_wp   ! pcm/K
        sim%alpha_void    = -80.0_wp  ! pcm/(% void)

        sim%rod_bank_position = 0.95_wp  ! start 95 % inserted

        ! ── Physics initialisations ──────────────────
        ! Saturation temperature at nominal pressure (sat_temp_K at 7.14 MPa ≈ 286°C = 559 K)
        sim%sat_temperature  = sat_temp_K(sim%pressure_operating)
        ! Water level: Lbase is effective water inventory; initialise to core_height + 1.5 m
        sim%Lbase            = sim%core_height + 1.5_wp   ! m
        sim%Lrx              = 1.5_wp                     ! m above core top
        sim%avg_coolant_temp = sim%inlet_temperature
        sim%coolant_T_prev1  = sim%inlet_temperature
        sim%coolant_T_prev2  = sim%inlet_temperature
        sim%P_prev1          = sim%pressure_operating
        sim%P_prev2          = sim%pressure_operating
        sim%T_excess         = 0.0_wp

        ! Period / APRM
        sim%neutrons_prev  = 0.0_wp
        sim%reactor_period = 1.0e6_wp  ! effectively infinite at start
        sim%aprm           = 0.0_wp

        ! Void reactivity factor at zero void (full moderator density)
        sim%rf_void = void_reactivity_factor(0.0_wp)  ! ≈ 1.0

        ! Turbine
        sim%turbine_valve   = 33.0_wp
        sim%bypass_valve    = 0.0_wp
        sim%steam_flow      = 0.0_wp
        sim%turbine_speed   = 3600.0_wp ! rpm – synchronised to grid
        sim%turbine_power   = 0.0_wp
        sim%turbine_tripped = .false.
        sim%feedwater_temp  = sim%inlet_temperature

        print *, "Initialising physics modules..."

        ! ── Cross-section library ─────────────────────────────────
        call xslib_init(sim%xslib, n_groups=2)

        ! Fuel cross-sections
        xsec_fuel%name       = "UO2_35"
        xsec_fuel%n_groups   = 2
        xsec_fuel%is_fuel    = .true.
        xsec_fuel%enrichment = 0.035_wp
        xsec_fuel%T_ref      = 900.0_wp
        xsec_fuel%rho_ref    = 10.97_wp

        call xslib_create_two_group_fuel(xsec_fuel%xsec_base, enrichment=0.035_wp)

        allocate(xsec_fuel%alpha_D(2),   xsec_fuel%alpha_mod(2))
        allocate(xsec_fuel%alpha_rho(2), xsec_fuel%alpha_void(2))
        xsec_fuel%alpha_D    = [-2.0e-5_wp, -3.0e-5_wp]
        xsec_fuel%alpha_mod  = [ 0.0_wp,    0.0_wp    ]
        xsec_fuel%alpha_rho  = [ 0.0_wp,    0.0_wp    ]
        xsec_fuel%alpha_void = [ 0.0_wp,    0.0_wp    ]
        call xslib_add_material(sim%xslib, xsec_fuel)

        ! Water cross-sections
        xsec_water%name    = "H2O"
        xsec_water%n_groups= 2
        xsec_water%is_fuel = .false.
        xsec_water%T_ref   = 560.0_wp
        xsec_water%rho_ref = 0.74_wp

        call xslib_create_two_group_moderator(xsec_water%xsec_base)

        allocate(xsec_water%alpha_D(2),   xsec_water%alpha_mod(2))
        allocate(xsec_water%alpha_rho(2), xsec_water%alpha_void(2))
        xsec_water%alpha_D    = [ 0.0_wp,    0.0_wp    ]
        xsec_water%alpha_mod  = [ 0.0_wp,    1.0e-4_wp ]
        xsec_water%alpha_rho  = [-10.0_wp,  -50.0_wp   ]
        xsec_water%alpha_void = [-10.0_wp, -100.0_wp   ]
        call xslib_add_material(sim%xslib, xsec_water)

        ! ── Neutronics ────────────────────────────────────────────
        neutron_config%n_groups         = 2
        neutron_config%max_outer_iter   = 100
        neutron_config%outer_tolerance  = 1.0e-5_wp
        neutron_config%power_level      = sim%power_rated
        neutron_config%normalize_power  = .true.
        call mg_init(sim%neutronics, sim%nx, sim%ny, sim%nz, &
             sim%dx*100.0_wp, sim%dy*100.0_wp, sim%dz*100.0_wp, neutron_config)
        call setup_bwr_geometry_realistic(sim)

        ! ── Heat transfer ─────────────────────────────────────────
        heat_config%include_convection = .true.
        call heat_init(sim%heat, sim%nx, sim%ny, sim%nz, &
                       sim%dx, sim%dy, sim%dz, heat_config)
        call heat_set_properties(sim%heat, k=3.0_wp, rho=10970.0_wp, cp=300.0_wp, &
            i1=1, i2=int(0.9*sim%nx), j1=1, j2=int(0.9*sim%ny))
        call heat_set_properties(sim%heat, k=0.6_wp, rho=738.0_wp, cp=5200.0_wp, &
            i1=int(0.9*sim%nx)+1, i2=sim%nx, &
            j1=int(0.9*sim%ny)+1, j2=sim%ny)

        ! ── Two-phase thermal-hydraulics ──────────────────────────
        th_config%void_correlation       = VOID_CHEXAL_LELLOUCHE_ID
        th_config%chf_correlation        = CHF_GROENEVELD_ID
        th_config%include_subcooled_boiling = .true.
        call two_phase_init(sim%thermalhydraulics, sim%nx, sim%ny, sim%nz, &
                            sim%dx, sim%dy, sim%dz, th_config)
        call two_phase_set_geometry(sim%thermalhydraulics, diameter=0.012_wp)

        ! ── Burnup / depletion ────────────────────────────────────
        burnup_config%track_xenon     = .true.
        burnup_config%track_samarium  = .true.
        burnup_config%track_actinides = .true.
        call burnup_init(sim%burnup, sim%nx, sim%ny, sim%nz, &
                         sim%dx, sim%dy, sim%dz, burnup_config)
        call burnup_set_initial_composition(sim%burnup, enrichment=0.035_wp)

        ! ── Time counters ─────────────────────────────────────────
        sim%time      = 0.0_wp
        sim%n_steps   = 0
        sim%initialized = .false.

        print *, "  Grid:          ", sim%nx, "×", sim%ny, "×", sim%nz
        print *, "  Core diameter: ", sim%core_diameter, " m"
        print *, "  Core height:   ", sim%core_height,   " m"
        print *, "  Rated power:   ", sim%power_rated/1.0e6_wp, " MW"
        print '(A,F6.1,A)', "  Sat. temp at ", sim%pressure_operating/1.0e6_wp, &
            " MPa: " // trim(adjustl(fmt_real(sim%sat_temperature - 273.15_wp))) // " °C"
        print *, ""

    end subroutine setup_bwr_simulation

    ! ================================================================
    !  PRESSURE DYNAMICS
    ! ================================================================
    subroutine update_pressure_dynamics(sim, dt)
        type(simulation_t), intent(inout) :: sim
        real(wp), intent(in) :: dt

        real(wp) :: T_sat, T_mean_smooth, rel_T, dP_gen, steam_val, dP_loss
        real(wp), parameter :: dP_COEF    = 5000.0_wp   ! Pa K⁻¹ s⁻¹
        real(wp), parameter :: P_LOSS_COEF= 1.125e6_wp  ! Pa s⁻¹ at 100 % valve, nominal P
        real(wp), parameter :: P_NOMINAL  = 7.14e6_wp   ! Pa

        ! 3-point smoothing of coolant temperature
        sim%coolant_T_prev2 = sim%coolant_T_prev1
        sim%coolant_T_prev1 = sim%avg_coolant_temp
        T_mean_smooth = (sim%avg_coolant_temp + sim%coolant_T_prev1 + sim%coolant_T_prev2) / 3.0_wp

        ! Saturation temperature at current pressure (sat_temp polynomial)
        T_sat = sat_temp_K(sim%pressure_operating)
        sim%sat_temperature = T_sat

        ! Flash steam excess: temperature above saturation
        sim%T_excess = max(0.0_wp, T_mean_smooth - T_sat)

        ! Relative temperature term
        !   (T_mean - 100°C) * (1 - Prx/10 MPa)
        rel_T = (T_mean_smooth - 373.15_wp) * (1.0_wp - sim%pressure_operating / 10.0e6_wp)

        ! Pressure rise from reactor heat
        dP_gen = dP_COEF * (rel_T + 30.0_wp * sim%T_excess) * dt

        ! Steam consumers: turbine + bypass reduce pressure
        if (.not. sim%turbine_tripped) then
            steam_val = (sim%turbine_valve + sim%bypass_valve) / 100.0_wp &
                        * (sim%pressure_operating / P_NOMINAL)
        else
            steam_val = 0.0_wp   ! isolation on turbine trip
        end if
        dP_loss = -P_LOSS_COEF * steam_val * dt

        ! Pressure history smoothing
        sim%P_prev2 = sim%P_prev1
        sim%P_prev1 = sim%pressure_operating
        sim%pressure_operating = max(3.0e6_wp, min(11.0e6_wp, &
            sim%pressure_operating + dP_gen + dP_loss))

        ! Update reactor level
        ! Lbase adjusts slightly as steam is produced (level drop) or condensate returns
        ! Net change: production proportional to power, return proportional to feed flow
        sim%Lbase = sim%Lbase &
                  - (sim%power_current / max(1.0_wp, 2381.0e6_wp)) * 0.002_wp * dt  ! steam out
        sim%Lbase = sim%Lbase + (sim%mass_flux_core / 1500.0_wp) * 0.0015_wp * dt   ! feed in
        sim%Lbase = max(0.5_wp, min(sim%core_height + 3.0_wp, sim%Lbase))

        sim%Lrx = calc_reactor_level(sim%Lbase, sim%avg_coolant_temp, sim%core_height)

        ! Feedwater temperature
        sim%feedwater_temp = sim%inlet_temperature &
                           + (sim%avg_coolant_temp - sim%inlet_temperature) * 0.1_wp

    end subroutine update_pressure_dynamics

    ! ================================================================
    !  INSTRUMENTATION: APRM, period, void reactivity factor
    ! ================================================================

    subroutine update_instrumentation(sim)
        type(simulation_t), intent(inout) :: sim

        real(wp) :: ratio

        ! APRM: percentage of rated power
        sim%aprm = (sim%power_current / max(1.0_wp, sim%power_rated)) * 100.0_wp
        sim%aprm = max(0.0_wp, min(200.0_wp, sim%aprm))

        ! Reactor period from successive power levels
        if (sim%neutrons_prev > 1.0e3_wp .and. sim%power_current > 1.0e3_wp) then
            ratio = sim%power_current / sim%neutrons_prev
            if (ratio > 0.0_wp .and. abs(ratio - 1.0_wp) > 1.0e-9_wp) then
                sim%reactor_period = 1.0_wp / log(ratio)
                ! Clamp to display range
                sim%reactor_period = max(-999.0_wp, min(999.0_wp, sim%reactor_period))
            else
                sim%reactor_period = 9999.0_wp  ! effectively stable
            end if
        end if
        sim%neutrons_prev = sim%power_current

        ! Void reactivity factor from volume-averaged void fraction
        sim%rf_void = void_reactivity_factor(sim%avg_void_fraction / 100.0_wp)

    end subroutine update_instrumentation

    ! ================================================================
    !  TURBINE MODEL
    ! ================================================================
    !
    ! The Fortran turbine is synchronised to the grid at 3600 rpm.
    ! Electrical output is computed from thermal power × thermodynamic
    ! efficiency × valve fraction (Rankine cycle approximation).
    ! On a trip the turbine isolates and steam flow ceases.
    !
    subroutine update_turbine(sim, dt)
        type(simulation_t), intent(inout) :: sim
        real(wp), intent(in) :: dt

        real(wp) :: eta, ePres, steam_to_turbine
        real(wp), parameter :: ETA_CARNOT  = 0.34_wp   ! ~BWR thermal efficiency
        real(wp), parameter :: RATED_RPM   = 3600.0_wp
        real(wp), parameter :: P_NOMINAL   = 7.14e6_wp

        if (sim%turbine_tripped) then
            ! Coastdown after trip
            sim%turbine_speed = max(0.0_wp, sim%turbine_speed - 30.0_wp * dt)
            sim%turbine_power = 0.0_wp
            sim%steam_flow    = 0.0_wp
            return
        end if

        ! Pressure efficiency factor
        ePres = min(1.0_wp, 6.1e-8_wp * sim%pressure_operating + 0.567_wp)

        ! Steam flow to turbine (stmT2 = tSP * 2000 * Prx/10e6 / 100)
        steam_to_turbine = sim%turbine_valve / 100.0_wp &
                         * 2000.0_wp * sim%pressure_operating / 10.0e6_wp
        sim%steam_flow = steam_to_turbine

        ! Turbine speed – assume grid synchronisation when above 3400 rpm
        ! (sw_Sync path, tSpd locked to 3600)
        if (sim%turbine_speed > 3400.0_wp) then
            sim%turbine_speed = RATED_RPM
            eta = ETA_CARNOT * ePres * (sim%turbine_valve / 100.0_wp)
            sim%turbine_power = eta * sim%power_current
        else
            ! Off-speed (startup / transient rundown)
            sim%turbine_speed = sim%turbine_speed + &
                (RATED_RPM * steam_to_turbine / 2000.0_wp - sim%turbine_speed) / 30.0_wp * dt
            sim%turbine_speed = max(0.0_wp, sim%turbine_speed)
            eta = ETA_CARNOT * ePres * (sim%turbine_speed / RATED_RPM) * (sim%turbine_valve / 100.0_wp)
            sim%turbine_power = eta * sim%power_current
        end if

        ! Safety: turbine trip on high condenser back-pressure (simplified)
        ! S.tTrip fires on genT>110, oilQ<30, or high vibration
        if (sim%turbine_tripped) sim%turbine_power = 0.0_wp

    end subroutine update_turbine

    ! ================================================================
    !  CROSS-SECTION FEEDBACK
    ! ================================================================
    !
    ! The rf_void factor applies a bulk moderator-density correction on
    ! top of the local cross-section feedback already handled by
    ! xslib_get_xsec.  This captures the effect of reactor water level
    ! changes (e.g. bulk voiding) that are not captured by local void
    ! fraction alone.
    !
    subroutine update_cross_sections_feedback(sim, T, rho)
        type(simulation_t), intent(inout) :: sim
        real(wp), intent(in) :: T(:, :, :)
        real(wp), intent(in) :: rho(:, :, :)

        integer  :: i, j, k
        type(mg_xsec_t) :: xsec
        real(wp) :: T_fuel, rho_mod, rf, rf_ref, rho_mod_corrected
        real(wp) :: node_bottom, node_top, rod_tip, inserted_fraction, H_core

        ! Reference rf at nominal void (30 % void) used for normalisation
        ! so the correction is unity at design conditions.
        rf_ref = void_reactivity_factor(0.30_wp)

        ! Current bulk void reactivity factor
        rf = sim%rf_void   ! updated in update_instrumentation

        do k = 1, sim%nz
            do j = 1, sim%ny
                do i = 1, sim%nx
                    T_fuel  = T(i, j, k)
                    rho_mod = rho(i, j, k)

                    ! scale effective moderator density
                    ! by the ratio rf/rf_ref.  rf < 1 at high bulk void →
                    ! reduced moderation density seen by the cross-section library.
                    rho_mod_corrected = rho_mod * (rf / max(1.0e-9_wp, rf_ref))

                    call xslib_get_xsec(sim%xslib, "UO2_35", &
                        T_fuel, rho_mod_corrected, sim%burnup%burnup(i, j, k), xsec, &
                        Xe_conc = sim%burnup%Xe135(i, j, k), &
                        Sm_conc = sim%burnup%Sm149(i, j, k))

                    ! Control rod insertion (BWR bottom-entry)
                    H_core      = real(sim%nz, wp) * sim%dz
                    node_bottom = real(k-1, wp) * sim%dz
                    node_top    = real(k,   wp) * sim%dz
                    rod_tip     = sim%rod_bank_position * H_core

                    inserted_fraction = 0.0_wp
                    if (rod_tip >= node_top) then
                        inserted_fraction = 1.0_wp
                    else if (rod_tip > node_bottom) then
                        inserted_fraction = (rod_tip - node_bottom) / sim%dz
                    end if

                    if (inserted_fraction > 0.0_wp) &
                        call xslib_apply_control_rod(xsec, inserted_fraction)

                    call mg_set_cross_sections(sim%neutronics, xsec, i, i, j, j, k, k)
                end do
            end do
        end do
    end subroutine update_cross_sections_feedback

    ! ================================================================
    !  COUPLED TIME STEP
    ! ================================================================

    subroutine coupled_time_step(sim, dt)
        type(simulation_t), intent(inout) :: sim
        real(wp), intent(in) :: dt

        real(wp) :: power(sim%nx, sim%ny, sim%nz)
        real(wp) :: flux(sim%nx, sim%ny, sim%nz)
        real(wp) :: temperature(sim%nx, sim%ny, sim%nz)
        real(wp) :: void_fraction(sim%nx, sim%ny, sim%nz)
        real(wp) :: density(sim%nx, sim%ny, sim%nz)
        real(wp), parameter :: rho_liquid = 738.0_wp
        real(wp), parameter :: rho_vapor  = 0.038_wp
        logical  :: converged
        real(wp) :: v_coolant, Re, h_conv
        real(wp), parameter :: D_h     = 0.01_wp
        real(wp), parameter :: mu      = 0.0001_wp
        real(wp), parameter :: k_fluid = 0.6_wp
        real(wp), parameter :: Pr      = 0.9_wp

        ! ── Convective flow setup ─────────────────────────────────
        v_coolant = sim%mass_flux_core / rho_liquid
        if (allocated(sim%heat%vz)) sim%heat%vz = v_coolant

        Re     = rho_liquid * v_coolant * D_h / mu
        h_conv = 0.023_wp * Re**0.8_wp * Pr**0.4_wp * k_fluid / D_h

        ! ── Neutronics ────────────────────────────────────────────
        call mg_solve_eigenvalue(sim%neutronics, sim%k_eff, converged)
        call mg_get_power(sim%neutronics, power)
        call mg_get_flux(sim%neutronics, flux)
        sim%power_current = sim%neutronics%total_power

        ! ── Burnup / depletion ────────────────────────────────────
        call burnup_step(sim%burnup, flux, power, dt)
        sim%avg_burnup = sum(sim%burnup%burnup) / size(sim%burnup%burnup)

        ! ── Heat transfer ─────────────────────────────────────────
        sim%heat%q = power
        call heat_step(sim%heat, dt)
        temperature       = sim%heat%T
        sim%max_fuel_temp = maxval(temperature)

        ! Update average coolant temperature
        sim%avg_coolant_temp = sum(temperature) / real(size(temperature), wp)

        ! ── Two-phase thermal-hydraulics ──────────────────────────
        call two_phase_step(sim%thermalhydraulics, &
            temperature, &
            sim%pressure_operating + 0.0_wp * temperature, &
            sim%mass_flux_core     + 0.0_wp * temperature, &
            power / (sim%dx * sim%dy), dt)

        call two_phase_get_void_fraction(sim%thermalhydraulics, void_fraction)

        density = (1.0_wp - void_fraction) * rho_liquid + void_fraction * rho_vapor
        density = max(density, rho_vapor)
        density = min(density, rho_liquid)

        sim%avg_void_fraction = sum(void_fraction) / size(void_fraction) * 100.0_wp
        sim%max_void_fraction = maxval(void_fraction) * 100.0_wp

        if (allocated(sim%thermalhydraulics%chf_ratio)) then
            sim%min_chfr = minval(sim%thermalhydraulics%chf_ratio)
        else
            sim%min_chfr = 999.0_wp
        end if

        ! ── Cross-section feedback ──
        call update_cross_sections_feedback(sim, temperature, density)

        ! APRM, period, void reactivity factor
        call update_instrumentation(sim)

        ! Pressure dynamics + reactor level
        call update_pressure_dynamics(sim, dt)

        ! Turbine speed and electrical output
        call update_turbine(sim, dt)

        ! ── Reactivity ───────────────────────────────────────────
        sim%reactivity_pcm = (sim%k_eff - 1.0_wp) / sim%k_eff * 1.0e5_wp

        sim%time    = sim%time + dt
        sim%n_steps = sim%n_steps + 1

    end subroutine coupled_time_step

    ! ================================================================
    !  STEADY-STATE SOLVE
    ! ================================================================

    subroutine solve_steady_state(sim)
        type(simulation_t), intent(inout) :: sim

        integer  :: iter
        real(wp) :: power(sim%nx, sim%ny, sim%nz)
        real(wp) :: temperature(sim%nx, sim%ny, sim%nz)
        real(wp) :: void_fraction(sim%nx, sim%ny, sim%nz)
        real(wp) :: density(sim%nx, sim%ny, sim%nz)
        real(wp) :: error
        real(wp), parameter :: rho_liquid = 0.738_wp
        real(wp), parameter :: rho_vapor  = 0.038_wp
        logical  :: converged

        do iter = 1, 50
            call mg_solve_eigenvalue(sim%neutronics, sim%k_eff, converged)
            if (.not. converged) &
                print *, "Warning: neutronics did not converge at iteration", iter

            call mg_get_power(sim%neutronics, power)
            sim%power_current = sim%neutronics%total_power

            sim%heat%q = power
            call heat_step_implicit(sim%heat, 1.0_wp)
            temperature       = sim%heat%T
            sim%max_fuel_temp = maxval(temperature)

            ! Update average coolant temp
            sim%avg_coolant_temp = sum(temperature) / real(size(temperature), wp)

            call two_phase_step(sim%thermalhydraulics, &
                temperature, &
                sim%pressure_operating + 0.0_wp * temperature, &
                sim%mass_flux_core     + 0.0_wp * temperature, &
                power / (sim%dx * sim%dy), 1.0_wp)

            call two_phase_get_void_fraction(sim%thermalhydraulics, void_fraction)

            density = (1.0_wp - void_fraction) * rho_liquid + void_fraction * rho_vapor
            density = max(density, rho_vapor)
            density = min(density, rho_liquid)

            sim%avg_void_fraction = sum(void_fraction) / size(void_fraction) * 100.0_wp
            sim%max_void_fraction = maxval(void_fraction) * 100.0_wp

            call update_cross_sections_feedback(sim, temperature, density)

            error = abs(sim%k_eff - 1.0_wp) + &
                    abs(sim%power_current - sim%power_rated) / sim%power_rated
            if (error < 1.0e-4_wp .and. iter > 3) exit
        end do

        sim%reactivity_pcm = (sim%k_eff - 1.0_wp) / sim%k_eff * 1.0e5_wp

        if (allocated(sim%thermalhydraulics%chf_ratio)) then
            sim%min_chfr = minval(sim%thermalhydraulics%chf_ratio)
        else
            sim%min_chfr = 999.0_wp
        end if

        ! Initialise instrumentation after first steady state
        call update_instrumentation(sim)
        call update_pressure_dynamics(sim, 0.0_wp)   ! compute sat_temp, Lrx
        call update_turbine(sim, 0.0_wp)             ! initial turbine power

    end subroutine solve_steady_state

    ! ================================================================
    !  REACTIVITY INSERTION
    ! ================================================================

    subroutine apply_reactivity_insertion(sim, rho_pcm)
        type(simulation_t), intent(inout) :: sim
        real(wp), intent(in) :: rho_pcm

        real(wp) :: factor
        integer  :: i, j, k, g

        factor = 1.0_wp + rho_pcm / 1.0e5_wp
        do k = 1, sim%nz
            do j = 1, sim%ny
                do i = 1, sim%nx
                    do g = 1, sim%neutronics%n_groups
                        sim%neutronics%xsec(i, j, k)%nu_sigma_f(g) = &
                            sim%neutronics%xsec(i, j, k)%nu_sigma_f(g) * factor
                    end do
                end do
            end do
        end do
    end subroutine apply_reactivity_insertion

    ! ================================================================
    !  AUTOMATIC CONTROLS
    ! ================================================================

    subroutine apply_automatic_controls(sim, dt)
        type(simulation_t), intent(inout) :: sim
        real(wp), intent(in) :: dt

        real(wp) :: power_fraction, rod_adjustment
        logical, save :: warning_printed = .false.

        real(wp), parameter :: SCRAM_TEMP        = 1673.15_wp   ! 1400 °C
        real(wp), parameter :: HIGH_TEMP         = 1473.15_wp   ! 1200 °C
        real(wp), parameter :: SCRAM_POWER       = 1.30_wp      ! 130 % rated
        real(wp), parameter :: HIGH_POWER        = 1.15_wp      ! 115 % rated
        real(wp), parameter :: SCRAM_PERIOD      = 5.0_wp       ! s (prompt supercrit. guard)
        real(wp), parameter :: SCRAM_LEVEL       = -1.5_wp      ! m below core floor
        real(wp), parameter :: SCRAM_PRESSURE_HI = 9.5e6_wp     ! Pa (high P)
        real(wp), parameter :: SCRAM_PRESSURE_LO = 4.0e6_wp     ! Pa (low P)

        power_fraction = sim%power_current / max(sim%power_rated, 1.0e6_wp)

        ! ── Emergency SCRAM conditions ────────────────────────────

        if (sim%max_fuel_temp > SCRAM_TEMP) then
            call scram(sim, "FUEL TEMPERATURE LIMIT", &
                "Max fuel: " // trim(deg_str(sim%max_fuel_temp)) // "°C  Limit: 1400°C")
            return
        end if

        if (power_fraction > SCRAM_POWER) then
            call scram(sim, "POWER EXCEEDED LIMIT", &
                "Power: " // trim(pct_str(power_fraction*100)) // "%  Limit: 130%")
            return
        end if

        ! period scram (prompt critical guard)
        if (sim%reactor_period > 0.0_wp .and. sim%reactor_period < SCRAM_PERIOD &
            .and. sim%power_current > 1.0e6_wp) then
            call scram(sim, "REACTOR PERIOD TOO SHORT", &
                "Period: " // trim(per_str(sim%reactor_period)) // " s  Limit: " // &
                trim(per_str(SCRAM_PERIOD)) // " s")
            return
        end if

        ! Reactor level scram (low water)
        if (sim%Lrx < SCRAM_LEVEL) then
            call scram(sim, "REACTOR WATER LEVEL LOW", &
                "Level: " // trim(lev_str(sim%Lrx)) // " m  Limit: " // &
                trim(lev_str(SCRAM_LEVEL)) // " m")
            return
        end if

        ! High/low pressure scram (pressure model)
        if (sim%pressure_operating > SCRAM_PRESSURE_HI) then
            call scram(sim, "HIGH REACTOR PRESSURE", &
                "Pressure: " // trim(mpa_str(sim%pressure_operating/1.0e6_wp)) // &
                " MPa  Limit: 9.5 MPa")
            return
        end if
        if (sim%pressure_operating < SCRAM_PRESSURE_LO .and. sim%power_current > 0.1e6_wp) then
            call scram(sim, "LOW REACTOR PRESSURE", &
                "Pressure: " // trim(mpa_str(sim%pressure_operating/1.0e6_wp)) // &
                " MPa  Limit: 4.0 MPa")
            return
        end if

        ! ── Automatic power runback ───────────────────────────────
        if (power_fraction > HIGH_POWER) then
            rod_adjustment = (power_fraction - HIGH_POWER) * 0.02_wp * dt / 0.05_wp
            sim%rod_bank_position = min(1.0_wp, sim%rod_bank_position + rod_adjustment)
            if (.not. warning_printed) then
                print *, "⚠ AUTO CONTROL: Power high, inserting rods"
                warning_printed = .true.
            end if
        else if (power_fraction < HIGH_POWER - 0.02_wp) then
            warning_printed = .false.
        end if

        if (sim%max_fuel_temp > HIGH_TEMP) then
            rod_adjustment = (sim%max_fuel_temp - HIGH_TEMP) / 500.0_wp * 0.05_wp * dt / 0.05_wp
            sim%rod_bank_position = min(1.0_wp, sim%rod_bank_position + rod_adjustment)
        end if

        if (sim%reactivity_pcm > 500.0_wp) then
            rod_adjustment = (sim%reactivity_pcm - 500.0_wp) / 10000.0_wp * dt / 0.05_wp
            sim%rod_bank_position = min(1.0_wp, sim%rod_bank_position + rod_adjustment)
        end if

        ! ── Auto pressure control (autoPres PD controller) ──
        ! Adjust turbine valve to maintain operating pressure near 7.14 MPa.
        if (.not. sim%turbine_tripped) then
            call auto_pressure_control(sim)
        end if

    end subroutine apply_automatic_controls

    !> Auto pressure control (autoPres PD controller).
    !  Adjusts turbine valve to maintain 7.14 MPa.
    subroutine auto_pressure_control(sim)
        type(simulation_t), intent(inout) :: sim
        real(wp) :: vel, dist, chg
        real(wp), parameter :: P_TARGET = 7.14e6_wp

        ! PD controller: vel=derivative, dist=proportional
        vel  = (sim%pressure_operating - (sim%P_prev1 + sim%P_prev2) / 2.0_wp) * 7.0_wp
        dist = sim%pressure_operating - P_TARGET
        chg  = max(-1.0_wp, min(1.0_wp, (dist + vel) / 1.0e5_wp))

        ! Positive change (overpressure) → open turbine/bypass valve
        if (chg > 0.0_wp) then
            if (sim%turbine_valve < 100.0_wp) then
                sim%turbine_valve = min(100.0_wp, sim%turbine_valve + chg)
            else
                sim%bypass_valve  = min(100.0_wp, sim%bypass_valve + chg)
            end if
        else
            ! Negative change (underpressure) → close bypass first, then turbine
            if (sim%bypass_valve > 0.0_wp) then
                sim%bypass_valve = max(0.0_wp, sim%bypass_valve + chg)
            else
                sim%turbine_valve = max(0.0_wp, sim%turbine_valve + chg)
            end if
        end if
    end subroutine auto_pressure_control

    ! ================================================================
    !  SAFETY LIMITS CHECK
    ! ================================================================

    function check_safety_limits(sim) result(approaching_limits)
        type(simulation_t), intent(in) :: sim
        logical :: approaching_limits

        approaching_limits = .false.

        if (sim%max_fuel_temp > 1373.0_wp) then
            approaching_limits = .true.
            print *, "  High fuel temperature: ", sim%max_fuel_temp, " K"
        end if
        if (sim%min_chfr < 1.5_wp) then
            approaching_limits = .true.
            print *, "  Low CHFR: ", sim%min_chfr
        end if
        if (sim%power_current > 1.2_wp * sim%power_rated) then
            approaching_limits = .true.
            print *, "  High power: ", sim%power_current / 1.0e6_wp, " MW"
        end if
        if (sim%Lrx < 0.0_wp) then
            approaching_limits = .true.
            print '(A,F6.2,A)', "  Low reactor level: ", sim%Lrx, " m"
        end if
        if (sim%pressure_operating > 9.0e6_wp) then
            approaching_limits = .true.
            print '(A,F6.3,A)', "  High pressure: ", sim%pressure_operating/1.0e6_wp, " MPa"
        end if
    end function check_safety_limits

    ! ================================================================
    !  REPORTING
    ! ================================================================

    subroutine print_steady_state_summary(sim)
        type(simulation_t), intent(in) :: sim

        print *, ""
        print *, "========== Steady-State Results =========="
        print *, "k_eff:                ", sim%k_eff
        print *, "Reactivity:           ", sim%reactivity_pcm, " pcm"
        print *, "Core power:           ", sim%power_current / 1.0e6_wp, " MW"
        print *, "Max fuel temp:        ", sim%max_fuel_temp - 273.15_wp, " °C"
        print *, "Avg void fraction:    ", sim%avg_void_fraction, " %"
        print *, "Max void fraction:    ", sim%max_void_fraction, " %"
        print *, "Min CHFR:             ", sim%min_chfr
        print *, "Sat. temperature:     ", sim%sat_temperature - 273.15_wp, " °C"
        print *, "Reactor level (Lrx):  ", sim%Lrx, " m above core floor"
        print *, "APRM:                 ", sim%aprm, " %"
        print *, "Void reactivity (rf): ", sim%rf_void
        print *, "Turbine valve:        ", sim%turbine_valve, " %"
        print *, "Turbine power:        ", sim%turbine_power / 1.0e6_wp, " MWe"
        print *, "Pressure:             ", sim%pressure_operating / 1.0e6_wp, " MPa"
        print *, "=========================================="
        print *, ""
    end subroutine print_steady_state_summary

    subroutine print_transient_summary(sim, step, t)
        type(simulation_t), intent(in) :: sim
        integer,  intent(in) :: step
        real(wp), intent(in) :: t

        write(*, '(A,I6,A,F8.2,A)') "Step ", step, ", t = ", t, " s"
        write(*, '(A,F10.5,A,F8.1,A)') &
            "  k_eff = ", sim%k_eff, ", ρ = ", sim%reactivity_pcm, " pcm"
        write(*, '(A,F8.1,A,A,F5.1,A)') &
            "  Power = ", sim%power_current/1.0e6_wp, " MW", &
            "  APRM = ", sim%aprm, " %"
        write(*, '(A,F7.1,A)') &
            "  Max T_fuel = ", sim%max_fuel_temp - 273.15_wp, " °C"
        write(*, '(A,F6.2,A,A,F5.2)') &
            "  Avg void = ", sim%avg_void_fraction, " %", "  Min CHFR = ", sim%min_chfr
        write(*, '(A,F6.2,A,A,F6.3,A)') &
            "  Level = ", sim%Lrx, " m", &
            "  Pressure = ", sim%pressure_operating/1.0e6_wp, " MPa"
        write(*, '(A,F8.1,A,A,F8.2,A)') &
            "  Turbine = ", sim%turbine_power/1.0e6_wp, " MWe", &
            "  Period = ", sim%reactor_period, " s"
        write(*, '(A,F6.2)') &
            "  rf_void = ", sim%rf_void
        print *, ""
    end subroutine print_transient_summary

    ! ================================================================
    !  CLEANUP
    ! ================================================================

    subroutine cleanup_simulation(sim)
        type(simulation_t), intent(inout) :: sim

        if (.not. sim%initialized) return

        call mg_destroy(sim%neutronics)
        call heat_destroy(sim%heat)
        call two_phase_destroy(sim%thermalhydraulics)
        call burnup_destroy(sim%burnup)
        call xslib_destroy(sim%xslib)

        sim%initialized = .false.
    end subroutine cleanup_simulation

    ! ================================================================
    !  HELPER PROCEDURES
    ! ================================================================

    !> Tiny inline helper to format a real for print statements
    function fmt_real(v) result(s)
        real(wp), intent(in) :: v
        character(len=8) :: s
        write(s,'(F8.1)') v
    end function fmt_real

    !> Emergency SCRAM: full rod insertion + turbine trip + pause
    subroutine scram(sim, reason, detail)
        type(simulation_t), intent(inout) :: sim
        character(len=*),   intent(in)    :: reason, detail
        print *, ""
        print *, "╔════════════════════════════════════════════╗"
        print *, "║   *** EMERGENCY REACTOR SCRAM ***          ║"
        print *, "║   CAUSE: " // reason
        print *, "║   " // detail
        print *, "║   ACTION: Full rod insertion               ║"
        print *, "╚════════════════════════════════════════════╝"
        print *, ""
        sim%rod_bank_position = 1.0_wp
        sim%turbine_tripped   = .true.
        g_viz%paused          = .true.
    end subroutine scram

    ! Tiny string formatters for scram messages
    function deg_str(v) result(s)
        real(wp), intent(in) :: v
        character(len=8) :: s
        write(s,'(F8.1)') v - 273.15_wp
    end function deg_str

    function pct_str(v) result(s)
        real(wp), intent(in) :: v
        character(len=8) :: s
        write(s,'(F8.1)') v
    end function pct_str

    function per_str(v) result(s)
        real(wp), intent(in) :: v
        character(len=8) :: s
        write(s,'(F8.2)') v
    end function per_str

    function lev_str(v) result(s)
        real(wp), intent(in) :: v
        character(len=8) :: s
        write(s,'(F8.2)') v
    end function lev_str

    function mpa_str(v) result(s)
        real(wp), intent(in) :: v
        character(len=8) :: s
        write(s,'(F8.3)') v
    end function mpa_str

end program bwr_core_simulator_opengl
