! core/compute_mode.f90
!
! Computation mode dispatch infrastructure.
!
! This module provides the compile-time and runtime framework for selecting
! between three execution modes:
!
!   CPU    — All computation on the CPU (gfortran + OpenMP + OpenBLAS).
!   GPU    — All GPU-eligible computation on the GPU (nvfortran + CUDA Fortran).
!            CPU is used only for orchestration and inherently serial work.
!   HYBRID — Work is split between GPU and CPU based on problem size.
!            Large parallel workloads go to the GPU; small/serial tasks
!            stay on the CPU. CUDA streams can overlap the two.
!
! The mode is selected at CMake configure time via -DCOMPUTATION_MODE=CPU|GPU|HYBRID
! and baked in as preprocessor macros. Callers query the mode via the
! get_compute_mode() function or the USE_GPU / USE_HYBRID preprocessor macros.
!
! Usage in application code:
!
!   use compute_mode
!   if (use_gpu()) then
!       call gpu_kernel(...)
!   else
!       call cpu_kernel(...)
!   end if
!
! For hybrid mode, the gpu_partition() function returns the fraction of
! work (0.0 to 1.0) that should be offloaded to the GPU, allowing the
! caller to split arrays and dispatch to both CPU and GPU concurrently.

module compute_mode
    use kinds, only: wp, i32
    implicit none
    private

    public :: COMPUTE_MODE_CPU, COMPUTE_MODE_GPU, COMPUTE_MODE_HYBRID
    public :: get_compute_mode, get_compute_mode_name
    public :: use_gpu, use_hybrid
    public :: gpu_partition, set_gpu_partition
    public :: gpu_min_workload, set_gpu_min_workload
    public :: get_num_threads, set_num_threads

    !> Computation mode identifiers.
    integer(i32), parameter :: COMPUTE_MODE_CPU    = 1
    integer(i32), parameter :: COMPUTE_MODE_GPU    = 2
    integer(i32), parameter :: COMPUTE_MODE_HYBRID = 3

    ! Runtime-adjustable parameters (mainly for hybrid mode).
    ! Default partition: 70% of work to GPU, 30% to CPU.
    real(wp), save :: gpu_work_fraction = 0.7_wp
    ! Minimum workload (number of cells / matrix size) below which GPU
    ! offloading is skipped even in GPU/HYBRID mode — the H2D transfer
    ! overhead dominates for small problems.
    integer(i32), save :: min_gpu_workload = 4096
    ! CPU thread count for hybrid mode (0 = use OMP default).
    integer(i32), save :: hybrid_num_threads = 0

contains

    !> Return the compile-time-selected computation mode.
    pure function get_compute_mode() result(mode)
        integer(i32) :: mode
#ifdef COMPUTATION_MODE_GPU
        mode = COMPUTE_MODE_GPU
#elif defined(COMPUTATION_MODE_HYBRID)
        mode = COMPUTE_MODE_HYBRID
#else
        mode = COMPUTE_MODE_CPU
#endif
    end function get_compute_mode

    !> Return a human-readable name for the current computation mode.
    pure function get_compute_mode_name() result(name)
        character(len=8) :: name
        integer(i32) :: mode
        mode = get_compute_mode()
        select case (mode)
        case (COMPUTE_MODE_GPU)
            name = 'GPU'
        case (COMPUTE_MODE_HYBRID)
            name = 'HYBRID'
        case default
            name = 'CPU'
        end select
    end function get_compute_mode_name

    !> Return .true. if GPU acceleration is compiled in (GPU or HYBRID mode).
    pure function use_gpu() result(flag)
        logical :: flag
#ifdef COMPUTATION_MODE_GPU
        flag = .true.
#elif defined(COMPUTATION_MODE_HYBRID)
        flag = .true.
#else
        flag = .false.
#endif
    end function use_gpu

    !> Return .true. if hybrid GPU/CPU mode is active.
    pure function use_hybrid() result(flag)
        logical :: flag
#ifdef COMPUTATION_MODE_HYBRID
        flag = .true.
#else
        flag = .false.
#endif
    end function use_hybrid

    !> Return the fraction of work to offload to the GPU in hybrid mode.
    !! In pure GPU mode this returns 1.0; in CPU mode it returns 0.0.
    pure function gpu_partition() result(frac)
        real(wp) :: frac
        integer(i32) :: mode
        mode = get_compute_mode()
        select case (mode)
        case (COMPUTE_MODE_GPU)
            frac = 1.0_wp
        case (COMPUTE_MODE_HYBRID)
            frac = gpu_work_fraction
        case default
            frac = 0.0_wp
        end select
    end function gpu_partition

    !> Override the default GPU work fraction (hybrid mode only).
    !! Value must be in [0.0, 1.0]; values outside are clamped.
    subroutine set_gpu_partition(frac)
        real(wp), intent(in) :: frac
        gpu_work_fraction = max(0.0_wp, min(1.0_wp, frac))
    end subroutine set_gpu_partition

    !> Return the minimum workload size below which GPU offloading is
    !! skipped. Used by callers to avoid H2D transfer overhead on small
    !! problems. Typical threshold: a few thousand cells.
    pure function gpu_min_workload() result(n)
        integer(i32) :: n
        n = min_gpu_workload
    end function gpu_min_workload

    !> Override the default minimum GPU workload threshold.
    subroutine set_gpu_min_workload(n)
        integer(i32), intent(in) :: n
        min_gpu_workload = max(0, n)
    end subroutine set_gpu_min_workload

    !> Return the number of CPU threads to use in hybrid mode.
    !! Returns 0 to indicate "use OpenMP default" (OMP_NUM_THREADS env var).
    pure function get_num_threads() result(n)
        integer(i32) :: n
        n = hybrid_num_threads
    end function get_num_threads

    !> Set the number of CPU threads for hybrid mode (0 = OMP default).
    subroutine set_num_threads(n)
        integer(i32), intent(in) :: n
        hybrid_num_threads = max(0, n)
    end subroutine set_num_threads

end module compute_mode
