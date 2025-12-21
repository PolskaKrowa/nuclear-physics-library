! fortran/kernels/ode/rk4.f90
!
! Fourth-order Runge-Kutta (RK4) method for solving ODEs.
! Classic fixed-step explicit integrator with fourth-order accuracy.
!
! For the ODE system: dy/dt = f(t, y)
!
! Usage:
!   use rk4
!   type(rk4_config_t) :: config
!   real(wp) :: y(n), t_span(2)
!   real(wp), allocatable :: t_out(:), y_out(:,:)
!   
!   config = rk4_config_t(dt=0.01_wp, max_steps=10000)
!   t_span = [0.0_wp, 10.0_wp]
!   call rk4_solve(ode_func, t_span, y, config, t_out, y_out, status)
!
module rk4
    use kinds, only: wp, i64
    use constants, only: TOL_DEFAULT
    implicit none
    private
    
    ! Public interface
    public :: rk4_config_t, rk4_status_t
    public :: rk4_solve, rk4_step
    public :: RK4_SUCCESS, RK4_ERR_MAX_STEPS, RK4_ERR_INVALID_ARG
    public :: RK4_ERR_FUNC_EVAL
    
    ! Status codes
    integer, parameter :: RK4_SUCCESS = 0
    integer, parameter :: RK4_ERR_MAX_STEPS = 1
    integer, parameter :: RK4_ERR_INVALID_ARG = 2
    integer, parameter :: RK4_ERR_FUNC_EVAL = 3
    
    ! Configuration type
    type :: rk4_config_t
        real(wp) :: dt = 0.01_wp           ! Time step size
        integer :: max_steps = 100000      ! Maximum integration steps
        integer :: output_every = 1        ! Store output every N steps
        logical :: dense_output = .false.  ! Store all steps vs. endpoints only
    end type rk4_config_t
    
    ! Status information
    type :: rk4_status_t
        integer :: code = RK4_SUCCESS      ! Status code
        integer :: steps_taken = 0         ! Number of steps executed
        integer :: func_evals = 0          ! Function evaluations
        real(wp) :: final_time = 0.0_wp    ! Final integration time reached
    end type rk4_status_t
    
    ! Abstract interface for ODE function
    abstract interface
        pure subroutine ode_func_interface(t, y, dydt)
            import :: wp
            real(wp), intent(in) :: t
            real(wp), intent(in) :: y(:)
            real(wp), intent(out) :: dydt(:)
        end subroutine ode_func_interface
    end interface
    
contains

    !> Solve ODE system using RK4 method
    !!
    !! @param func     ODE function f(t, y) = dy/dt
    !! @param t_span   Integration interval [t0, tf]
    !! @param y0       Initial conditions
    !! @param config   Solver configuration
    !! @param t_out    Output time points (allocated by routine)
    !! @param y_out    Output solution (allocated by routine)
    !! @param status   Status information
    subroutine rk4_solve(func, t_span, y0, config, t_out, y_out, status)
        procedure(ode_func_interface) :: func
        real(wp), intent(in) :: t_span(2)
        real(wp), intent(in) :: y0(:)
        type(rk4_config_t), intent(in) :: config
        real(wp), allocatable, intent(out) :: t_out(:)
        real(wp), allocatable, intent(out) :: y_out(:, :)
        type(rk4_status_t), intent(out) :: status
        
        real(wp) :: t, tf, dt
        real(wp) :: y(size(y0))
        integer :: n, step, max_outputs, output_idx
        
        ! Initialize
        n = size(y0)
        t = t_span(1)
        tf = t_span(2)
        dt = config%dt
        y = y0
        status = rk4_status_t()
        
        ! Validate inputs
        if (dt <= 0.0_wp .or. tf <= t .or. n <= 0) then
            status%code = RK4_ERR_INVALID_ARG
            return
        end if
        
        ! Estimate output size
        if (config%dense_output) then
            max_outputs = ceiling((tf - t) / dt) + 1
        else
            max_outputs = ceiling((tf - t) / (dt * config%output_every)) + 1
        end if
        
        ! Allocate output arrays
        allocate(t_out(max_outputs))
        allocate(y_out(n, max_outputs))
        
        ! Store initial condition
        output_idx = 1
        t_out(output_idx) = t
        y_out(:, output_idx) = y
        
        ! Integration loop
        do step = 1, config%max_steps
            ! Adjust final step if needed
            if (t + dt > tf) then
                dt = tf - t
            end if
            
            ! Take RK4 step
            call rk4_step(func, t, y, dt, y)
            t = t + dt
            
            status%steps_taken = step
            status%func_evals = status%func_evals + 4
            
            ! Store output if needed
            if (config%dense_output .or. mod(step, config%output_every) == 0 &
                .or. abs(t - tf) < epsilon(1.0_wp)) then
                output_idx = output_idx + 1
                if (output_idx <= max_outputs) then
                    t_out(output_idx) = t
                    y_out(:, output_idx) = y
                end if
            end if
            
            ! Check if we've reached the final time
            if (abs(t - tf) < epsilon(1.0_wp)) then
                status%code = RK4_SUCCESS
                status%final_time = t
                exit
            end if
        end do
        
        ! Trim output arrays to actual size
        if (output_idx < max_outputs) then
            t_out = t_out(1:output_idx)
            y_out = y_out(:, 1:output_idx)
        end if
        
        ! Check if max steps exceeded
        if (status%code /= RK4_SUCCESS) then
            status%code = RK4_ERR_MAX_STEPS
            status%final_time = t
        end if
    end subroutine rk4_solve
    
    !> Perform single RK4 integration step
    !!
    !! @param func   ODE function f(t, y) = dy/dt
    !! @param t      Current time
    !! @param y      Current state (input)
    !! @param dt     Time step size
    !! @param y_new  New state (output)
    pure subroutine rk4_step(func, t, y, dt, y_new)
        procedure(ode_func_interface) :: func
        real(wp), intent(in) :: t
        real(wp), intent(in) :: y(:)
        real(wp), intent(in) :: dt
        real(wp), intent(out) :: y_new(:)
        
        real(wp) :: k1(size(y)), k2(size(y)), k3(size(y)), k4(size(y))
        real(wp) :: y_temp(size(y))
        
        ! k1 = f(t, y)
        call func(t, y, k1)
        
        ! k2 = f(t + dt/2, y + dt*k1/2)
        y_temp = y + 0.5_wp * dt * k1
        call func(t + 0.5_wp * dt, y_temp, k2)
        
        ! k3 = f(t + dt/2, y + dt*k2/2)
        y_temp = y + 0.5_wp * dt * k2
        call func(t + 0.5_wp * dt, y_temp, k3)
        
        ! k4 = f(t + dt, y + dt*k3)
        y_temp = y + dt * k3
        call func(t + dt, y_temp, k4)
        
        ! y_new = y + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
        y_new = y + (dt / 6.0_wp) * (k1 + 2.0_wp * k2 + 2.0_wp * k3 + k4)
    end subroutine rk4_step

end module rk4