! examples/bwr_core_simulator_opengl.f90
!
! Enhanced BWR Core Simulator with OpenGL Visualisation
! Integrates full physics simulation with interactive 3D display
!
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
    
    ! OpenGL interfaces (simplified and safer)
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
    
    ! OpenGL constants
    integer(c_int), parameter :: GLUT_RGB = int(z'0000', c_int)
    integer(c_int), parameter :: GLUT_DOUBLE = int(z'0002', c_int)
    integer(c_int), parameter :: GLUT_DEPTH = int(z'0010', c_int)
    integer(c_int), parameter :: GL_COLOR_BUFFER_BIT = int(z'00004000', c_int)
    integer(c_int), parameter :: GL_DEPTH_BUFFER_BIT = int(z'00000100', c_int)
    integer(c_int), parameter :: GL_DEPTH_TEST = int(z'0B71', c_int)
    integer(c_int), parameter :: GL_PROJECTION = int(z'1701', c_int)
    integer(c_int), parameter :: GL_MODELVIEW = int(z'1700', c_int)
    integer(c_int), parameter :: GL_QUADS = int(z'0007', c_int)
    integer(c_int), parameter :: GL_LINES = int(z'0001', c_int)
    integer(c_int), parameter :: GL_BLEND = int(z'0BE2', c_int)
    integer(c_int), parameter :: GL_SRC_ALPHA = int(z'0302', c_int)
    integer(c_int), parameter :: GL_ONE_MINUS_SRC_ALPHA = int(z'0303', c_int)
    
    ! Simulation state (same as original)
    type :: simulation_t
        ! Grid
        integer :: nx, ny, nz
        real(wp) :: dx, dy, dz
        real(wp) :: core_height, core_diameter
        real(wp) :: rod_bank_position  ! 0.0 (All Out) to 1.0 (All In)
        
        ! Physics modules
        type(mg_state_t) :: neutronics
        type(heat_state_t) :: heat
        type(two_phase_state_t) :: thermalhydraulics
        type(burnup_state_t) :: burnup
        type(xsec_library_t) :: xslib
        
        ! Operating conditions
        real(wp) :: power_rated, power_current
        real(wp) :: pressure_operating
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
        integer :: n_steps
        
        ! Initialization flag
        logical :: initialized
    end type simulation_t
    
    ! Visualisation state
    type :: viz_state_t
        integer :: display_mode  ! 1=power, 2=temperature, 3=void, 4=burnup
        real(wp) :: rotation_x, rotation_y
        real(wp) :: zoom
        logical :: paused
        logical :: show_help
        integer :: update_counter
    end type viz_state_t
    
    ! Global state (must be SAVE for OpenGL callbacks)
    type(simulation_t), target, save :: g_sim
    type(viz_state_t), save :: g_viz
    real(wp), save :: g_dt = 0.01_wp
    logical, save :: g_converged = .false.
    integer(c_int), save :: g_window
    
    print *, "=============================================="
    print *, "  BWR Core Simulator with OpenGL"
    print *, "  Interactive 3D Visualisation"
    print *, "=============================================="
    print *, ""
    
    ! Initialize RNG
    call rng_seed(123456789_8)
    
    ! Setup simulation - MUST be done before OpenGL
    print *, "Setting up simulation..."
    call setup_bwr_simulation(g_sim)
    g_sim%initialized = .true.
    
    ! Initial steady-state solve
    print *, "Computing initial steady-state..."
    call solve_steady_state(g_sim)
    call print_steady_state_summary(g_sim)
    
    ! Initialize visualisation state
    g_viz%display_mode = 1  ! Power distribution
    g_viz%rotation_x = 30.0_wp
    g_viz%rotation_y = 45.0_wp
    g_viz%zoom = -15.0_wp
    g_viz%paused = .false.
    g_viz%show_help = .true.
    g_viz%update_counter = 0
    
    print *, ""
    print *, "=== CONTROLS ==="
    print *, "  SPACE  - Pause/Resume"
    print *, "  1-4    - Change display (Power/Temp/Void/Burnup)"
    print *, "  +/-    - Increase/Decrease power"
    print *, "  r/R    - Insert positive/negative reactivity"
    print *, "  c/C    - Increase/Decrease coolant flow"
    print *, "  i/o    - Inseart/Withdraw control rods"
    print *, "  w/s    - Rotate up/down"
    print *, "  a/d    - Rotate left/right"
    print *, "  z/x    - Zoom in/out"
    print *, "  h      - Toggle help"
    print *, "  q/ESC  - Quit"
    print *, "================"
    print *, ""
    print *, "Starting OpenGL visualisation..."
    
    ! Initialize OpenGL
    call init_opengl()
    
    ! Start main loop (never returns)
    call glutMainLoop()
    
    ! Cleanup (only reached on error)
    call cleanup_simulation(g_sim)
    
contains

    !> Initialize OpenGL window and callbacks
    subroutine init_opengl()
        character(kind=c_char, len=25) :: title
        type(c_ptr) :: null_ptr

        ! 1. Define variables for command line arguments
        integer(c_int), target :: argc = 0       ! Must be 'target' to use c_loc
        type(c_ptr), target    :: argv = c_null_ptr
        
        title = "BWR Core Simulator"//c_null_char
        
        ! Initialize GLUT with null pointers (safer)
        null_ptr = c_null_ptr
        call glutInit(c_loc(argc), c_loc(argv))
        
        call glutInitDisplayMode(ior(ior(GLUT_RGB, GLUT_DOUBLE), GLUT_DEPTH))
        call glutInitWindowSize(1200_c_int, 900_c_int)
        g_window = glutCreateWindow(title)
        
        ! Set callbacks
        call glutDisplayFunc(c_funloc(display_callback))
        call glutReshapeFunc(c_funloc(reshape_callback))
        call glutKeyboardFunc(c_funloc(keyboard_callback))
        call glutIdleFunc(c_funloc(idle_callback))

        call glEnable(GL_BLEND)
        call glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
            
        ! OpenGL settings
        call glClearColor(0.1_c_float, 0.1_c_float, 0.15_c_float, 1.0_c_float)
        call glEnable(GL_DEPTH_TEST)
    end subroutine init_opengl
    
    !> Display callback - renders the scene
    subroutine display_callback() bind(C)
        integer :: i, j, k, skip
        real(wp) :: x, y, z, value, r, g, b
        real(wp) :: val_min, val_max
        
        ! Safety check
        if (.not. g_sim%initialized) return
        
        call glClear(ior(GL_COLOR_BUFFER_BIT, GL_DEPTH_BUFFER_BIT))
        call glMatrixMode(GL_MODELVIEW)
        call glLoadIdentity()
        
        ! Set up view
        call glTranslatef(0.0_c_float, 0.0_c_float, real(g_viz%zoom, c_float))
        call glRotatef(real(g_viz%rotation_x, c_float), 1.0_c_float, 0.0_c_float, 0.0_c_float)
        call glRotatef(real(g_viz%rotation_y, c_float), 0.0_c_float, 1.0_c_float, 0.0_c_float)
        
        ! Draw coordinate axes
        call draw_axes()
        
        ! Determine value range for colour mapping
        val_min = 0.0_wp
        val_max = 1.0_wp
        
        if (g_viz%display_mode == 1 .and. allocated(g_sim%neutronics%power_density)) then
            val_min = minval(g_sim%neutronics%power_density)
            val_max = maxval(g_sim%neutronics%power_density)
        else if (g_viz%display_mode == 2 .and. allocated(g_sim%heat%T)) then
            val_min = minval(g_sim%heat%T) - 273.15_wp
            val_max = maxval(g_sim%heat%T) - 273.15_wp
        else if (g_viz%display_mode == 3) then
            val_min = 0.0_wp
            val_max = 0.8_wp  ! Max void fraction
        else if (g_viz%display_mode == 4) then
            val_min = 0.0_wp
            if (allocated(g_sim%burnup%burnup)) then
                val_max = max(maxval(g_sim%burnup%burnup), 1.0_wp)
            end if
        end if
        
        ! Draw reactor core (skip some cells for performance)
        skip = max(1, g_sim%nx / 15)  ! Reduce resolution for performance
        
        do k = 1, g_sim%nz, skip
            do j = 1, g_sim%ny, skip
                do i = 1, g_sim%nx, skip
                    ! Get value for this cell with safety checks
                    value = 0.0_wp
                    
                    select case(g_viz%display_mode)
                    case(1)  ! Power
                        if (allocated(g_sim%neutronics%power_density)) then
                            value = g_sim%neutronics%power_density(i, j, k)
                        end if
                    case(2)  ! Temperature
                        if (allocated(g_sim%heat%T)) then
                            value = g_sim%heat%T(i, j, k) - 273.15_wp
                        end if
                    case(3)  ! Void
                        if (allocated(g_sim%thermalhydraulics%void_fraction)) then
                            value = g_sim%thermalhydraulics%void_fraction(i, j, k)
                        end if
                    case(4)  ! Burnup
                        if (allocated(g_sim%burnup%burnup)) then
                            value = g_sim%burnup%burnup(i, j, k)
                        end if
                    end select
                    
                    ! Convert to colour
                    call value_to_colour(value, val_min, val_max, r, g, b)
                    
                    ! Calculate position (centred)
                    x = (real(i, wp) - real(g_sim%nx, wp) / 2.0_wp) * g_sim%dx
                    y = (real(j, wp) - real(g_sim%ny, wp) / 2.0_wp) * g_sim%dy
                    z = (real(k, wp) - real(g_sim%nz, wp) / 2.0_wp) * g_sim%dz
                    
                    ! Draw voxel
                    call draw_cube(x, y, z, g_sim%dx * 0.4_wp * skip, r, g, b, 0.2_wp)
                end do
            end do
        end do
        
        ! Draw text overlay
        call draw_text_overlay()
        
        call glutSwapBuffers()
    end subroutine display_callback
    
    !> Reshape callback
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
    
    !> Keyboard callback
    subroutine keyboard_callback(key, x, y) bind(C)
        integer(c_int), value :: key, x, y
        
        select case(key)
        case(32)  ! SPACE
            g_viz%paused = .not. g_viz%paused
            print *, merge("PAUSED  ", "RESUMED ", g_viz%paused)
            
        case(49)  ! '1'
            g_viz%display_mode = 1
            print *, "Display: Power Distribution"
            
        case(50)  ! '2'
            g_viz%display_mode = 2
            print *, "Display: Temperature"
            
        case(51)  ! '3'
            g_viz%display_mode = 3
            print *, "Display: Void Fraction"
            
        case(52)  ! '4'
            g_viz%display_mode = 4
            print *, "Display: Burnup"
            
        case(43, 61)  ! '+' or '='
            g_sim%power_rated = g_sim%power_rated * 1.05_wp
            print '(A,F8.1,A)', " Power setpoint: ", g_sim%power_rated/1e6_wp, " MW"
            
        case(45)  ! '-'
            g_sim%power_rated = g_sim%power_rated * 0.95_wp
            print '(A,F8.1,A)', " Power setpoint: ", g_sim%power_rated/1e6_wp, " MW"
            
        case(114)  ! 'r'
            call apply_reactivity_insertion(g_sim, 20.0_wp)
            print *, " +20 pcm reactivity inserted"
            
        case(82)  ! 'R'
            call apply_reactivity_insertion(g_sim, -20.0_wp)
            print *, " -20 pcm reactivity inserted"
            
        case(99)  ! 'c'
            g_sim%mass_flux_core = g_sim%mass_flux_core * 1.05_wp
            print '(A,F7.1,A)', " Flow: ", g_sim%mass_flux_core, " kg/m²·s"
            
        case(67)  ! 'C'
            g_sim%mass_flux_core = g_sim%mass_flux_core * 0.95_wp
            print '(A,F7.1,A)', " Flow: ", g_sim%mass_flux_core, " kg/m²·s"
            
        case(119)  ! 'w'
            g_viz%rotation_x = g_viz%rotation_x + 5.0_wp
            
        case(115)  ! 's'
            g_viz%rotation_x = g_viz%rotation_x - 5.0_wp
            
        case(97)  ! 'a'
            g_viz%rotation_y = g_viz%rotation_y - 5.0_wp
            
        case(100)  ! 'd'
            g_viz%rotation_y = g_viz%rotation_y + 5.0_wp
            
        case(122)  ! 'z'
            g_viz%zoom = g_viz%zoom + 1.0_wp
            
        case(120)  ! 'x'
            g_viz%zoom = g_viz%zoom - 1.0_wp
            
        case(104)  ! 'h'
            g_viz%show_help = .not. g_viz%show_help
        
        ! [I]nsert Rods (Move Up)
        case (ichar('i'), ichar('I'))
            g_sim%rod_bank_position = min(1.0_wp, g_sim%rod_bank_position + 0.005_wp)
            print *, "Rods Inserting: ", g_sim%rod_bank_position * 100.0, "%"

        ! [O]ut / Withdraw Rods (Move Down)
        case (ichar('o'), ichar('O'))
            g_sim%rod_bank_position = max(0.0_wp, g_sim%rod_bank_position - 0.005_wp)
            print *, "Rods Withdrawing: ", g_sim%rod_bank_position * 100.0, "%"
            
        case(27, 113)  ! ESC or 'q'
            print *, "Shutting down..."
            call cleanup_simulation(g_sim)
            stop
        end select
        
        call glutPostRedisplay()
    end subroutine keyboard_callback
    
    !> Idle callback - physics update
    subroutine idle_callback() bind(C)
        if (.not. g_viz%paused .and. g_sim%initialized) then
            ! Run physics time step
            call coupled_time_step(g_sim, g_dt)
            call apply_automatic_controls(g_sim, g_dt)
            
            g_viz%update_counter = g_viz%update_counter + 1
            
            ! Periodic console output
            if (mod(g_viz%update_counter, 100) == 0) then
                call print_transient_summary(g_sim, g_sim%n_steps, g_sim%time)
            end if
            
            ! Safety check
            ! disabled because funny
            !if (check_safety_limits(g_sim)) then
            !    g_viz%paused = .true.
            !    print *, "*** SIMULATION PAUSED - Safety limits approached ***"
            !end if
        end if
        
        call glutPostRedisplay()
    end subroutine idle_callback
    
    !> Draw a coloured cube
    subroutine draw_cube(x, y, z, size, r, g, b, alpha)
        real(wp), intent(in) :: x, y, z, size, r, g, b, alpha
        real(c_float) :: xf, yf, zf, sf
        
        xf = real(x, c_float)
        yf = real(y, c_float)
        zf = real(z, c_float)
        sf = real(size, c_float)
        
        call glColor4f(real(r, c_float), real(g, c_float), real(b, c_float), real(alpha, c_float))
        call glBegin(GL_QUADS)
        
        ! Front
        call glVertex3f(xf - sf, yf - sf, zf + sf)
        call glVertex3f(xf + sf, yf - sf, zf + sf)
        call glVertex3f(xf + sf, yf + sf, zf + sf)
        call glVertex3f(xf - sf, yf + sf, zf + sf)
        
        ! Back
        call glVertex3f(xf - sf, yf - sf, zf - sf)
        call glVertex3f(xf - sf, yf + sf, zf - sf)
        call glVertex3f(xf + sf, yf + sf, zf - sf)
        call glVertex3f(xf + sf, yf - sf, zf - sf)
        
        ! Top
        call glVertex3f(xf - sf, yf + sf, zf - sf)
        call glVertex3f(xf - sf, yf + sf, zf + sf)
        call glVertex3f(xf + sf, yf + sf, zf + sf)
        call glVertex3f(xf + sf, yf + sf, zf - sf)
        
        ! Bottom
        call glVertex3f(xf - sf, yf - sf, zf - sf)
        call glVertex3f(xf + sf, yf - sf, zf - sf)
        call glVertex3f(xf + sf, yf - sf, zf + sf)
        call glVertex3f(xf - sf, yf - sf, zf + sf)
        
        ! Right
        call glVertex3f(xf + sf, yf - sf, zf - sf)
        call glVertex3f(xf + sf, yf + sf, zf - sf)
        call glVertex3f(xf + sf, yf + sf, zf + sf)
        call glVertex3f(xf + sf, yf - sf, zf + sf)
        
        ! Left
        call glVertex3f(xf - sf, yf - sf, zf - sf)
        call glVertex3f(xf - sf, yf - sf, zf + sf)
        call glVertex3f(xf - sf, yf + sf, zf + sf)
        call glVertex3f(xf - sf, yf + sf, zf - sf)
        
        call glEnd()
    end subroutine draw_cube
    
    !> Draw coordinate axes
    subroutine draw_axes()
        real(c_float), parameter :: len = 3.0_c_float
        
        call glBegin(GL_LINES)
        
        ! X axis (red)
        call glColor3f(1.0_c_float, 0.0_c_float, 0.0_c_float)
        call glVertex3f(0.0_c_float, 0.0_c_float, 0.0_c_float)
        call glVertex3f(len, 0.0_c_float, 0.0_c_float)
        
        ! Y axis (green)
        call glColor3f(0.0_c_float, 1.0_c_float, 0.0_c_float)
        call glVertex3f(0.0_c_float, 0.0_c_float, 0.0_c_float)
        call glVertex3f(0.0_c_float, len, 0.0_c_float)
        
        ! Z axis (blue)
        call glColor3f(0.0_c_float, 0.0_c_float, 1.0_c_float)
        call glVertex3f(0.0_c_float, 0.0_c_float, 0.0_c_float)
        call glVertex3f(0.0_c_float, 0.0_c_float, len)
        
        call glEnd()
    end subroutine draw_axes
    
    !> Convert value to heat-map colour
    subroutine value_to_colour(value, vmin, vmax, r, g, b)
        real(wp), intent(in) :: value, vmin, vmax
        real(wp), intent(out) :: r, g, b
        real(wp) :: t
        
        ! Normalise to [0, 1]
        if (abs(vmax - vmin) > 1.0e-10_wp) then
            t = (value - vmin) / (vmax - vmin)
        else
            t = 0.5_wp
        end if
        t = max(0.0_wp, min(1.0_wp, t))
        
        ! Heat map: blue -> cyan -> green -> yellow -> red
        if (t < 0.25_wp) then
            r = 0.0_wp
            g = 4.0_wp * t
            b = 1.0_wp
        else if (t < 0.5_wp) then
            r = 0.0_wp
            g = 1.0_wp
            b = 1.0_wp - 4.0_wp * (t - 0.25_wp)
        else if (t < 0.75_wp) then
            r = 4.0_wp * (t - 0.5_wp)
            g = 1.0_wp
            b = 0.0_wp
        else
            r = 1.0_wp
            g = 1.0_wp - 4.0_wp * (t - 0.75_wp)
            b = 0.0_wp
        end if
    end subroutine value_to_colour
    
    !> Draw text overlay with reactor status
    subroutine draw_text_overlay()
        ! Simplified text rendering - prints status to console instead
        ! Full implementation would use glutBitmapString
    end subroutine draw_text_overlay
    
    subroutine setup_bwr_geometry_realistic(sim)
        type(simulation_t), intent(inout) :: sim
        
        integer :: i, j, k
        real(wp) :: x, y, r
        real(wp) :: x_center, y_center
        real(wp) :: pitch, pin_radius
        integer :: pin_i, pin_j
        real(wp) :: dx_from_pin, dy_from_pin, r_from_pin
        logical :: in_fuel
        logical :: found_water, found_fuel
        
        type(xsec_material_t) :: xsec_fuel, xsec_water
        
        ! Get materials from library
        call xslib_get_material(sim%xslib, "UO2_35", xsec_fuel, found_fuel)
        if (.not. found_fuel) then
            print *, 'Fatal: material not found - aborting.'
            call xslib_list_materials(sim%xslib)  ! optional
            stop 1
        end if
        call xslib_get_material(sim%xslib, "H2O", xsec_water, found_water)
        if (.not. found_water) then
            print *, 'Fatal: material not found - aborting.'
            call xslib_list_materials(sim%xslib)  ! optional
            stop 1
        end if
        
        ! BWR fuel assembly parameters
        pitch = 0.0163_wp  ! 1.63 cm pin pitch (typical for BWR)
        pin_radius = 0.0052_wp  ! 5.2 mm fuel pellet radius
        
        ! Fuel-to-moderator ratio check
        print *, "BWR lattice parameters:"
        print '(A,F7.4,A)', "  Pin pitch:           ", pitch * 100, " cm"
        print '(A,F7.4,A)', "  Pin radius:          ", pin_radius * 100, " cm"
        print '(A,F7.2,A)', "  Fuel-to-mod ratio:   ", &
            (3.14159_wp * pin_radius**2) / (pitch**2 - 3.14159_wp * pin_radius**2)
        
        ! Loop through all cells and determine if fuel or water
        do k = 1, sim%nz
            do j = 1, sim%ny
                do i = 1, sim%nx
                    
                    ! Position of cell center
                    x = (real(i, wp) - 0.5_wp) * sim%dx
                    y = (real(j, wp) - 0.5_wp) * sim%dy
                    
                    ! Find nearest pin center
                    pin_i = nint(x / pitch)
                    pin_j = nint(y / pitch)
                    
                    x_center = real(pin_i, wp) * pitch
                    y_center = real(pin_j, wp) * pitch
                    
                    ! Distance from pin center
                    dx_from_pin = x - x_center
                    dy_from_pin = y - y_center
                    r_from_pin = sqrt(dx_from_pin**2 + dy_from_pin**2)
                    
                    ! Determine if inside fuel pin
                    in_fuel = (r_from_pin < pin_radius)
                    
                    ! Also check if inside core boundary (cylindrical)
                    r = sqrt((x - sim%core_diameter/2)**2 + (y - sim%core_diameter/2)**2)
                    if (r > sim%core_diameter / 2.2_wp) then
                        in_fuel = .false.  ! Outside core = reflector (water)
                    end if
                    
                    ! Set cross-sections
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
    
    !> Setup BWR simulation
    subroutine setup_bwr_simulation(sim)
        type(simulation_t), intent(out) :: sim
        
        type(mg_config_t) :: neutron_config
        type(heat_config_t) :: heat_config
        type(two_phase_config_t) :: th_config
        type(burnup_config_t) :: burnup_config
        
        integer :: i, j, k
        type(xsec_material_t) :: xsec_fuel, xsec_water
        real(wp) :: r, r_fuel, r_water
        
        ! Core geometry
        sim%nx = 20
        sim%ny = 20
        sim%nz = 20
        
        sim%core_height = 3.71_wp
        sim%core_diameter = 4.75_wp
        
        sim%dx = sim%core_diameter / real(sim%nx, wp)
        sim%dy = sim%core_diameter / real(sim%ny, wp)
        sim%dz = sim%core_height / real(sim%nz, wp)
        
        ! Operating conditions
        sim%power_rated = 2381.0e6_wp
        sim%pressure_operating = 7.14e6_wp
        sim%mass_flux_core = 1500.0_wp
        sim%inlet_temperature = 551.0_wp
        
        ! Feedback coefficients
        sim%alpha_doppler = -3.5_wp
        sim%alpha_void = -80.0_wp

        sim%rod_bank_position = 1.0_wp ! Start 100% inserted
        
        print *, "Initializing physics modules..."
        
        ! Cross-section library
        call xslib_init(sim%xslib, n_groups=2)
        
        ! Fuel cross sections
        xsec_fuel%name = "UO2_35"
        xsec_fuel%n_groups = 2
        xsec_fuel%is_fuel = .true.
        xsec_fuel%enrichment = 0.035_wp
        xsec_fuel%T_ref = 900.0_wp
        xsec_fuel%rho_ref = 10.97_wp
        
        call xslib_create_two_group_fuel(xsec_fuel%xsec_base, enrichment=0.035_wp)
        
        allocate(xsec_fuel%alpha_D(2))
        allocate(xsec_fuel%alpha_mod(2))
        allocate(xsec_fuel%alpha_rho(2))
        allocate(xsec_fuel%alpha_void(2))
        
        xsec_fuel%alpha_D = [-2.0e-5_wp, -3.0e-5_wp]
        xsec_fuel%alpha_mod = [0.0_wp, 0.0_wp]
        xsec_fuel%alpha_rho = [0.0_wp, 0.0_wp]
        xsec_fuel%alpha_void = [0.0_wp, 0.0_wp]
        
        call xslib_add_material(sim%xslib, xsec_fuel)
        
        ! Water cross sections
        xsec_water%name = "H2O"
        xsec_water%n_groups = 2
        xsec_water%is_fuel = .false.
        xsec_water%T_ref = 560.0_wp
        xsec_water%rho_ref = 0.74_wp
        
        call xslib_create_two_group_moderator(xsec_water%xsec_base)
        
        allocate(xsec_water%alpha_D(2))
        allocate(xsec_water%alpha_mod(2))
        allocate(xsec_water%alpha_rho(2))
        allocate(xsec_water%alpha_void(2))
        
        xsec_water%alpha_D = [0.0_wp, 0.0_wp]
        xsec_water%alpha_mod = [0.0_wp, 1.0e-4_wp]
        xsec_water%alpha_rho = [-10.0_wp, -50.0_wp]
        xsec_water%alpha_void = [-10.0_wp, -100.0_wp]
        
        call xslib_add_material(sim%xslib, xsec_water)
        
        ! Neutronics
        neutron_config%n_groups = 2
        neutron_config%max_outer_iter = 100
        neutron_config%outer_tolerance = 1.0e-5_wp
        neutron_config%power_level = sim%power_rated
        neutron_config%normalize_power = .true.
        
        call mg_init(sim%neutronics, sim%nx, sim%ny, sim%nz, &
                    sim%dx, sim%dy, sim%dz, neutron_config)
        
        ! Setup geometry
        call setup_bwr_geometry_realistic(sim)
        
        ! Heat transfer
        heat_config%include_convection = .true.
        call heat_init(sim%heat, sim%nx, sim%ny, sim%nz, &
                      sim%dx, sim%dy, sim%dz, heat_config)
        
        call heat_set_properties(sim%heat, &
            k=3.0_wp, rho=10970.0_wp, cp=300.0_wp, &
            i1=1, i2=int(0.9*sim%nx), j1=1, j2=int(0.9*sim%ny))
        
        call heat_set_properties(sim%heat, &
            k=0.6_wp, rho=738.0_wp, cp=5200.0_wp)
        
        ! Two-phase flow
        th_config%void_correlation = VOID_CHEXAL_LELLOUCHE_ID
        th_config%chf_correlation = CHF_GROENEVELD_ID
        th_config%include_subcooled_boiling = .true.
        
        call two_phase_init(sim%thermalhydraulics, sim%nx, sim%ny, sim%nz, &
                           sim%dx, sim%dy, sim%dz, th_config)
        
        call two_phase_set_geometry(sim%thermalhydraulics, diameter=0.012_wp)
        
        ! Burnup
        burnup_config%track_xenon = .true.
        burnup_config%track_samarium = .true.
        burnup_config%track_actinides = .true.
        
        call burnup_init(sim%burnup, sim%nx, sim%ny, sim%nz, &
                        sim%dx, sim%dy, sim%dz, burnup_config)
        
        call burnup_set_initial_composition(sim%burnup, enrichment=0.035_wp)
        
        ! Initialize time
        sim%time = 0.0_wp
        sim%n_steps = 0
        sim%initialized = .false.
        
        print *, "  Grid: ", sim%nx, "x", sim%ny, "x", sim%nz
        print *, "  Core diameter: ", sim%core_diameter, " m"
        print *, "  Core height: ", sim%core_height, " m"
        print *, "  Rated power: ", sim%power_rated/1.0e6_wp, " MW"
        print *, ""
    end subroutine setup_bwr_simulation
    
    !> Automatic reactor protection system
    !> Prevents power excursions and overtemperature conditions
    !> Add this subroutine to the 'contains' section
    subroutine apply_automatic_controls(sim, dt)
        type(simulation_t), intent(inout) :: sim
        real(wp), intent(in) :: dt
        
        real(wp) :: power_fraction, power_error, rod_adjustment
        real(wp) :: temp_celsius, max_temp_celsius
        real(wp) :: reactivity_pcm
        logical, save :: warning_printed = .false.
        
        ! Temperature limits
        real(wp), parameter :: SCRAM_TEMP = 1673.15_wp     ! 1400°C - automatic scram
        real(wp), parameter :: HIGH_TEMP = 1473.15_wp      ! 1200°C - runback starts
        real(wp), parameter :: NORMAL_TEMP = 1273.15_wp    ! 1000°C - normal max
        
        ! Power limits
        real(wp), parameter :: SCRAM_POWER = 1.30_wp       ! 130% rated power
        real(wp), parameter :: HIGH_POWER = 1.15_wp        ! 115% rated power
        
        real(wp) :: power_rate, max_rate

        ! Get current state
        power_fraction = sim%power_current / max(sim%power_rated, 1.0e6_wp)
        max_temp_celsius = maxval(sim%heat%T) - 273.15_wp
        reactivity_pcm = sim%reactivity_pcm
        
        ! ================================================================
        ! EMERGENCY SCRAM CONDITIONS
        ! ================================================================
        
        ! Temperature scram
        if (sim%max_fuel_temp > SCRAM_TEMP) then
            print *, ""
            print *, "╔════════════════════════════════════════════╗"
            print *, "║   *** EMERGENCY REACTOR SCRAM ***          ║"
            print *, "║   CAUSE: FUEL TEMPERATURE LIMIT            ║"
            print '(A,F7.1,A)', " ║   Temperature: ", max_temp_celsius, " °C                  ║"
            print *, "║   Limit: 1400 °C                           ║"
            print *, "║   ACTION: Full rod insertion               ║"
            print *, "╚════════════════════════════════════════════╝"
            print *, ""
            sim%rod_bank_position = 1.0_wp
            g_viz%paused = .true.
            return
        end if
        
        ! Power scram
        if (power_fraction > SCRAM_POWER) then
            print *, ""
            print *, "╔════════════════════════════════════════════╗"
            print *, "║   *** EMERGENCY REACTOR SCRAM ***          ║"
            print *, "║   CAUSE: POWER EXCEEDED LIMIT              ║"
            print '(A,F6.1,A)', " ║   Power: ", power_fraction * 100.0_wp, " %                        ║"
            print *, "║   Limit: 130%                              ║"
            print *, "║   ACTION: Full rod insertion               ║"
            print *, "╚════════════════════════════════════════════╝"
            print *, ""
            sim%rod_bank_position = 1.0_wp
            g_viz%paused = .true.
            return
        end if
        
        ! ================================================================
        ! AUTOMATIC POWER RUNBACK
        ! ================================================================
        
        ! High power condition - gradually insert rods
        if (power_fraction > HIGH_POWER) then
            rod_adjustment = (power_fraction - HIGH_POWER) * 0.02_wp * dt / 0.05_wp
            sim%rod_bank_position = min(1.0_wp, sim%rod_bank_position + rod_adjustment)
            
            if (.not. warning_printed) then
                print *, ""
                print *, "⚠ AUTO CONTROL: Power high, inserting rods"
                warning_printed = .true.
            end if
        else if (power_fraction < HIGH_POWER - 0.02_wp) then
            warning_printed = .false.
        end if
        
        ! High temperature condition - more aggressive rod insertion
        if (sim%max_fuel_temp > HIGH_TEMP) then
            rod_adjustment = (sim%max_fuel_temp - HIGH_TEMP) / 500.0_wp * 0.05_wp * dt / 0.05_wp
            sim%rod_bank_position = min(1.0_wp, sim%rod_bank_position + rod_adjustment)
            
            if (mod(int(sim%time * 10), 50) == 0) then  ! Print every 5 seconds
                print '(A,F7.1,A)', " ⚠ AUTO: High temp (", max_temp_celsius, &
                    " °C), inserting rods"
            end if
        end if
        
        ! ================================================================
        ! REACTIVITY LIMITING
        ! ================================================================
        
        ! If reactivity is too positive, insert rods
        if (reactivity_pcm > 500.0_wp) then
            rod_adjustment = (reactivity_pcm - 500.0_wp) / 10000.0_wp * dt / 0.05_wp
            sim%rod_bank_position = min(1.0_wp, sim%rod_bank_position + rod_adjustment)
            
            if (mod(int(sim%time * 10), 50) == 0) then
                print '(A,F7.1,A)', " ⚠ AUTO: High reactivity (", reactivity_pcm, &
                    " pcm), limiting"
            end if
        end if
        
        ! ================================================================
        ! POWER RATE LIMITING
        ! ================================================================
        
        ! Limit rate of power increase (prevent prompt critical excursions)
        if (sim%n_steps > 1) then
            
            ! Maximum allowed power increase rate: 10% per second
            max_rate = 0.10_wp * sim%power_rated * dt
            power_rate = (sim%power_current - sim%neutronics%total_power) / dt
            
            if (power_rate > max_rate) then
                rod_adjustment = 0.001_wp * dt / 0.05_wp
                sim%rod_bank_position = min(1.0_wp, sim%rod_bank_position + rod_adjustment)
            end if
        end if
        
    end subroutine apply_automatic_controls

    subroutine coupled_time_step(sim, dt)
        type(simulation_t), intent(inout) :: sim
        real(wp), intent(in) :: dt
        
        real(wp) :: power(sim%nx, sim%ny, sim%nz)
        real(wp) :: flux(sim%nx, sim%ny, sim%nz)
        real(wp) :: temperature(sim%nx, sim%ny, sim%nz)
        real(wp) :: void_fraction(sim%nx, sim%ny, sim%nz)
        real(wp) :: density(sim%nx, sim%ny, sim%nz)
        real(wp), parameter :: rho_liquid = 0.738_wp
        real(wp), parameter :: rho_vapor = 0.038_wp
        logical :: converged
        
        real(wp) :: v_coolant
        integer :: i, j, k
        real(wp) :: Re, h_conv
        real(wp), parameter :: D_h = 0.01_wp       ! hydraulic diameter
        real(wp), parameter :: mu = 0.0001_wp      ! viscosity
        real(wp), parameter :: k_fluid = 0.6_wp    ! thermal conductivity
        real(wp), parameter :: Pr = 0.9_wp         ! Prandtl number for water
        
        ! Set coolant velocity field for heat transfer
        ! Approximate axial flow velocity from mass flux
        ! v_z = G / rho  where G is mass flux [kg/m^2.s]
        
        v_coolant = sim%mass_flux_core / rho_liquid
        
        ! Set velocity in heat transfer module
        if (allocated(sim%heat%vz)) then
            sim%heat%vz = v_coolant
        end if
        
        ! Also update convective boundary conditions
        ! Higher flow = higher heat transfer coefficient
        
        ! Reynolds number
        Re = rho_liquid * v_coolant * D_h / mu
        
        ! Dittus-Boelter correlation for heat transfer coefficient
        h_conv = 0.023_wp * Re**0.8_wp * Pr**0.4_wp * k_fluid / D_h
        
        ! Update boundary conditions in heat transfer
        do k = 1, sim%nz
            do j = 1, sim%ny
                do i = 1, sim%nx
                    ! Apply convective cooling in heat source term
                    
                end do
            end do
        end do
        
        ! ========================================================================
        ! Continue with original coupling
        ! ========================================================================
        
        ! Neutronics
        call mg_solve_eigenvalue(sim%neutronics, sim%k_eff, converged)
        
        ! Power distribution
        call mg_get_power(sim%neutronics, power)
        call mg_get_flux(sim%neutronics, flux)
        sim%power_current = sim%neutronics%total_power
        
        ! Burnup
        call burnup_step(sim%burnup, flux, power, dt)
        sim%avg_burnup = sum(sim%burnup%burnup) / size(sim%burnup%burnup)
        
        ! Heat transfer
        sim%heat%q = power
        call heat_step(sim%heat, dt)
        temperature = sim%heat%T
        sim%max_fuel_temp = maxval(temperature)
        
        ! Two-phase thermal-hydraulics
        call two_phase_step(sim%thermalhydraulics, &
            temperature, &
            sim%pressure_operating + 0.0_wp * temperature, &
            sim%mass_flux_core + 0.0_wp * temperature, &
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
        
        ! Update cross sections with feedback
        call update_cross_sections_feedback(sim, temperature, density)
        
        sim%reactivity_pcm = (sim%k_eff - 1.0_wp) / sim%k_eff * 1.0e5_wp
        
        sim%time = sim%time + dt
        sim%n_steps = sim%n_steps + 1
        
    end subroutine coupled_time_step

    !> Solve for initial steady state
    subroutine solve_steady_state(sim)
        type(simulation_t), intent(inout) :: sim
        
        integer :: iter
        real(wp) :: power(sim%nx, sim%ny, sim%nz)
        real(wp) :: temperature(sim%nx, sim%ny, sim%nz)
        real(wp) :: void_fraction(sim%nx, sim%ny, sim%nz)
        real(wp) :: density(sim%nx, sim%ny, sim%nz)
        real(wp) :: error
        real(wp), parameter :: rho_liquid = 0.738_wp
        real(wp), parameter :: rho_vapor = 0.038_wp
        logical :: converged
        
        do iter = 1, 50
            call mg_solve_eigenvalue(sim%neutronics, sim%k_eff, converged)
            
            if (.not. converged) then
                print *, "Warning: Neutronics did not converge at iteration", iter
            end if
            
            call mg_get_power(sim%neutronics, power)
            sim%power_current = sim%neutronics%total_power
            
            sim%heat%q = power
            call heat_step(sim%heat, 1.0_wp)
            
            temperature = sim%heat%T
            sim%max_fuel_temp = maxval(temperature)
            
            call two_phase_step(sim%thermalhydraulics, &
                temperature, &
                sim%pressure_operating + 0.0_wp * temperature, &
                sim%mass_flux_core + 0.0_wp * temperature, &
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
        
    end subroutine solve_steady_state
    
    !> Update cross sections with feedback
    subroutine update_cross_sections_feedback(sim, T, rho)
        type(simulation_t), intent(inout) :: sim
        real(wp), intent(in) :: T(:, :, :)
        real(wp), intent(in) :: rho(:, :, :)
        
        integer :: i, j, k
        type(mg_xsec_t) :: xsec
        real(wp) :: T_fuel, rho_mod
        
        do k = 1, sim%nz
            do j = 1, sim%ny
                do i = 1, sim%nx
                    T_fuel = T(i, j, k)
                    rho_mod = rho(i, j, k)
                    
                    ! Get cross-sections with feedback
                    call xslib_get_xsec(sim%xslib, "UO2_35", &
                        T_fuel, rho_mod, sim%burnup%burnup(i, j, k), xsec, &
                        Xe_conc=sim%burnup%Xe135(i, j, k), &
                        Sm_conc=sim%burnup%Sm149(i, j, k))
                    
                    ! ================================================================
                    ! CRITICAL FIX: Apply control rod BEFORE setting cross-sections!
                    ! ================================================================
                    block
                        real(wp) :: node_bottom, node_top, rod_tip, inserted_fraction
                        real(wp) :: H_core

                        H_core = real(sim%nz, wp) * sim%dz 
                        
                        ! Calculate physical height of this node
                        node_bottom = real(k-1, wp) * sim%dz
                        node_top    = real(k, wp) * sim%dz
                        
                        ! BWR rods come from the BOTTOM
                        ! Position 0.0 = fully withdrawn (bottom)
                        ! Position 1.0 = fully inserted (top)
                        rod_tip = sim%rod_bank_position * H_core
                        
                        inserted_fraction = 0.0_wp
                        
                        if (rod_tip >= node_top) then
                            ! Rod is completely covering this node
                            inserted_fraction = 1.0_wp
                        else if (rod_tip > node_bottom) then
                            ! Rod tip is partially inside this node
                            inserted_fraction = (rod_tip - node_bottom) / sim%dz
                        end if

                        ! Apply control rod to LOCAL xsec variable
                        ! BEFORE calling mg_set_cross_sections!
                        if (inserted_fraction > 0.0_wp) then
                            call xslib_apply_control_rod(xsec, inserted_fraction)
                        end if
                    end block

                    ! Now set the modified cross-sections
                    call mg_set_cross_sections(sim%neutronics, xsec, i, i, j, j, k, k)
                end do
            end do
        end do
    end subroutine update_cross_sections_feedback
    
    !> Apply reactivity insertion
    subroutine apply_reactivity_insertion(sim, rho_pcm)
        type(simulation_t), intent(inout) :: sim
        real(wp), intent(in) :: rho_pcm
        
        real(wp) :: factor
        integer :: i, j, k, g
        
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
    
    !> Check safety limits
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
    end function check_safety_limits
    
    !> Print steady-state summary
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
        print *, "=========================================="
        print *, ""
    end subroutine print_steady_state_summary
    
    !> Print transient summary
    subroutine print_transient_summary(sim, step, t)
        type(simulation_t), intent(in) :: sim
        integer, intent(in) :: step
        real(wp), intent(in) :: t
        
        write(*, '(A,I6,A,F8.2,A)') "Step ", step, ", t = ", t, " s"
        write(*, '(A,F10.5,A,F8.1,A)') "  k_eff = ", sim%k_eff, &
            ", ρ = ", sim%reactivity_pcm, " pcm"
        write(*, '(A,F8.1,A)') "  Power = ", sim%power_current / 1.0e6_wp, " MW"
        write(*, '(A,F7.1,A)') "  Max T_fuel = ", sim%max_fuel_temp - 273.15_wp, " °C"
        write(*, '(A,F6.2,A)') "  Avg void = ", sim%avg_void_fraction, " %"
        write(*, '(A,F6.2)') "  Min CHFR = ", sim%min_chfr
        print *, ""
    end subroutine print_transient_summary
    
    !> Cleanup
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

end program bwr_core_simulator_opengl