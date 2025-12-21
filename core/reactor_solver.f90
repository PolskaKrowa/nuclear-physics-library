! fortran/c_interface/reactor_solver.f90
!
! ISO_C_BINDING wrapper for ODE solvers to be called from C
!
module reactor_solver_c
    use, intrinsic :: iso_c_binding
    use kinds, only: wp, dp
    use rk4, only: rk4_config_t, rk4_status_t, rk4_solve
    use dormand_prince, only: dopri_config_t, dopri_status_t, dopri_solve
    implicit none
    
    ! C-compatible structures
    type, bind(C) :: c_rk4_config
        real(c_double) :: dt
        integer(c_int) :: max_steps
        integer(c_int) :: output_every
        integer(c_int) :: dense_output
    end type c_rk4_config
    
    type, bind(C) :: c_rk4_status
        integer(c_int) :: code
        integer(c_int) :: steps_taken
        integer(c_int) :: func_evals
        real(c_double) :: final_time
    end type c_rk4_status
    
    type, bind(C) :: c_dopri_config
        real(c_double) :: rtol
        real(c_double) :: atol
        real(c_double) :: dt_init
        real(c_double) :: dt_max
        real(c_double) :: dt_min
        real(c_double) :: safety
        integer(c_int) :: max_steps
        integer(c_int) :: dense_output
    end type c_dopri_config
    
    type, bind(C) :: c_dopri_status
        integer(c_int) :: code
        integer(c_int) :: steps_taken
        integer(c_int) :: steps_accepted
        integer(c_int) :: steps_rejected
        integer(c_int) :: func_evals
        real(c_double) :: final_time
        real(c_double) :: final_step_size
    end type c_dopri_status
    
    ! Store function pointer from C
    abstract interface
        subroutine c_ode_function(t, y, dydt, n) bind(C)
            use, intrinsic :: iso_c_binding
            real(c_double), intent(in) :: t
            integer(c_int), value :: n
            real(c_double), intent(in) :: y(n)
            real(c_double), intent(out) :: dydt(n)
        end subroutine c_ode_function
    end interface
    
    ! Module variable to store C function pointer
    procedure(c_ode_function), pointer :: c_func_ptr => null()
    
contains

    !> Wrapper to convert C function to Fortran pure subroutine
    pure subroutine fortran_ode_wrapper(t, y, dydt)
        real(wp), intent(in) :: t
        real(wp), intent(in) :: y(:)
        real(wp), intent(out) :: dydt(:)
        
        integer(c_int) :: n
        real(c_double) :: t_c, y_c(size(y)), dydt_c(size(y))
        
        n = int(size(y), c_int)
        t_c = real(t, c_double)
        y_c = real(y, c_double)
        
        ! This won't work with pure - need to redesign
        ! For now, we'll use a module variable approach
        dydt_c = 0.0_c_double
        
        dydt = real(dydt_c, wp)
    end subroutine fortran_ode_wrapper
    
    !> RK4 solver callable from C
    subroutine c_rk4_solve(func, t_span, y0, n, config, &
                          t_out, y_out, n_out, status) bind(C, name='rk4_solve_c')
        type(c_funptr), value :: func
        real(c_double), intent(in) :: t_span(2)
        integer(c_int), value :: n
        real(c_double), intent(in) :: y0(n)
        type(c_rk4_config), intent(in) :: config
        type(c_ptr), intent(out) :: t_out
        type(c_ptr), intent(out) :: y_out
        integer(c_int), intent(out) :: n_out
        type(c_rk4_status), intent(out) :: status
        
        ! Fortran types
        type(rk4_config_t) :: f_config
        type(rk4_status_t) :: f_status
        real(wp), allocatable :: f_t_out(:), f_y_out(:,:)
        real(wp) :: f_t_span(2), f_y0(n)
        
        ! Convert C types to Fortran types
        f_config%dt = real(config%dt, wp)
        f_config%max_steps = config%max_steps
        f_config%output_every = config%output_every
        f_config%dense_output = (config%dense_output /= 0)
        
        f_t_span = real(t_span, wp)
        f_y0 = real(y0, wp)
        
        ! Call Fortran solver (simplified - needs proper function pointer handling)
        ! call rk4_solve(fortran_ode_wrapper, f_t_span, f_y0, f_config, &
        !                f_t_out, f_y_out, f_status)
        
        ! For now, allocate dummy output
        allocate(f_t_out(100), f_y_out(n, 100))
        n_out = 100
        
        ! Convert back to C types
        status%code = f_status%code
        status%steps_taken = f_status%steps_taken
        status%func_evals = f_status%func_evals
        status%final_time = real(f_status%final_time, c_double)
        
        ! Allocate C-compatible arrays and copy data
        ! (Caller must free these)
        t_out = c_loc(f_t_out(1))
        y_out = c_loc(f_y_out(1,1))
    end subroutine c_rk4_solve
    
    !> Dormand-Prince solver callable from C
    subroutine c_dopri_solve(func, t_span, y0, n, config, &
                            t_out, y_out, n_out, status) bind(C, name='dopri_solve_c')
        type(c_funptr), value :: func
        real(c_double), intent(in) :: t_span(2)
        integer(c_int), value :: n
        real(c_double), intent(in) :: y0(n)
        type(c_dopri_config), intent(in) :: config
        type(c_ptr), intent(out) :: t_out
        type(c_ptr), intent(out) :: y_out
        integer(c_int), intent(out) :: n_out
        type(c_dopri_status), intent(out) :: status
        
        type(dopri_config_t) :: f_config
        type(dopri_status_t) :: f_status
        real(wp), allocatable :: f_t_out(:), f_y_out(:,:)
        real(wp) :: f_t_span(2), f_y0(n)
        
        ! Convert C config to Fortran config
        f_config%rtol = real(config%rtol, wp)
        f_config%atol = real(config%atol, wp)
        f_config%dt_init = real(config%dt_init, wp)
        f_config%dt_max = real(config%dt_max, wp)
        f_config%dt_min = real(config%dt_min, wp)
        f_config%safety = real(config%safety, wp)
        f_config%max_steps = config%max_steps
        f_config%dense_output = (config%dense_output /= 0)
        
        f_t_span = real(t_span, wp)
        f_y0 = real(y0, wp)
        
        ! Call Fortran solver (simplified)
        ! call dopri_solve(fortran_ode_wrapper, f_t_span, f_y0, f_config, &
        !                  f_t_out, f_y_out, f_status)
        
        ! Dummy allocation
        allocate(f_t_out(100), f_y_out(n, 100))
        n_out = 100
        
        ! Convert status back
        status%code = f_status%code
        status%steps_taken = f_status%steps_taken
        status%steps_accepted = f_status%steps_accepted
        status%steps_rejected = f_status%steps_rejected
        status%func_evals = f_status%func_evals
        status%final_time = real(f_status%final_time, c_double)
        status%final_step_size = real(f_status%final_step_size, c_double)
        
        t_out = c_loc(f_t_out(1))
        y_out = c_loc(f_y_out(1,1))
    end subroutine c_dopri_solve
    
    !> Simple single-step RK4 for easier C integration
    subroutine c_rk4_step(func, t, y, dt, y_new, n) bind(C, name='rk4_step_c')
        type(c_funptr), value :: func
        real(c_double), value :: t, dt
        integer(c_int), value :: n
        real(c_double), intent(in) :: y(n)
        real(c_double), intent(out) :: y_new(n)
        
        procedure(c_ode_function), pointer :: f_ptr
        real(c_double) :: k1(n), k2(n), k3(n), k4(n)
        real(c_double) :: temp(n), t_half, t_end
        
        ! Convert C function pointer to Fortran procedure pointer
        call c_f_procpointer(func, f_ptr)
        
        ! RK4 algorithm
        call f_ptr(t, y, k1, n)
        
        temp = y + 0.5_c_double * dt * k1
        t_half = t + 0.5_c_double * dt
        call f_ptr(t_half, temp, k2, n)
        
        temp = y + 0.5_c_double * dt * k2
        call f_ptr(t_half, temp, k3, n)
        
        temp = y + dt * k3
        t_end = t + dt
        call f_ptr(t_end, temp, k4, n)
        
        y_new = y + (dt / 6.0_c_double) * (k1 + 2.0_c_double*k2 + &
                                            2.0_c_double*k3 + k4)
    end subroutine c_rk4_step

end module reactor_solver_c