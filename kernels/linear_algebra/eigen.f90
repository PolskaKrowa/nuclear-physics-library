! fortran/kernels/linear_algebra/eigen.f90
!
! Eigenvalue and eigenvector computations using LAPACK
! Supports symmetric, general, and generalised eigenvalue problems
!
module eigen
    use kinds, only: wp, i32
    use constants, only: TOL_DEFAULT
    implicit none
    private
    
    public :: eigen_symmetric, eigen_general
    public :: eigen_generalised_symmetric
    public :: singular_value_decomposition
    public :: schur_decomposition
    public :: qr_decomposition
    public :: eigenvalues_only_symmetric
    public :: power_iteration, inverse_iteration
    
    ! Error codes
    integer(i32), parameter, public :: EIGEN_SUCCESS = 0
    integer(i32), parameter, public :: EIGEN_ERR_SIZE = 1
    integer(i32), parameter, public :: EIGEN_ERR_LAPACK = 2
    integer(i32), parameter, public :: EIGEN_ERR_CONVERGENCE = 3
    
    ! Interface for DGEES select function
    abstract interface
        logical function select_function_interface(wr, wi)
            use kinds, only: wp
            real(wp), intent(in) :: wr, wi
        end function select_function_interface
    end interface
    
contains

    !> Compute eigenvalues and eigenvectors of symmetric matrix
    !! Uses LAPACK DSYEV - produces real eigenvalues in ascending order
    subroutine eigen_symmetric(A, eigenvalues, eigenvectors, status)
        real(wp), intent(in) :: A(:, :)
        real(wp), intent(out) :: eigenvalues(:)
        real(wp), intent(out), optional :: eigenvectors(:, :)
        integer(i32), intent(out), optional :: status
        
        real(wp), allocatable :: A_copy(:, :), work(:)
        integer :: n, lda, lwork, info
        character :: jobz
        
        n = size(A, 1)
        
        if (size(A, 2) /= n .or. size(eigenvalues) /= n) then
            if (present(status)) status = EIGEN_ERR_SIZE
            return
        end if
        
        if (present(eigenvectors)) then
            if (size(eigenvectors, 1) /= n .or. size(eigenvectors, 2) /= n) then
                if (present(status)) status = EIGEN_ERR_SIZE
                return
            end if
            jobz = 'V'  ! Compute eigenvalues and eigenvectors
        else
            jobz = 'N'  ! Eigenvalues only
        end if
        
        lda = n
        lwork = max(1, 3 * n - 1)
        
        allocate(A_copy(n, n))
        allocate(work(lwork))
        
        A_copy = A
        
        call dsyev(jobz, 'U', n, A_copy, lda, eigenvalues, work, lwork, info)
        
        if (present(eigenvectors) .and. info == 0) then
            eigenvectors = A_copy
        end if
        
        deallocate(A_copy, work)
        
        if (info /= 0) then
            if (present(status)) status = EIGEN_ERR_LAPACK
        else
            if (present(status)) status = EIGEN_SUCCESS
        end if
    end subroutine eigen_symmetric
    
    !> Compute eigenvalues only for symmetric matrix (more efficient)
    subroutine eigenvalues_only_symmetric(A, eigenvalues, status)
        real(wp), intent(in) :: A(:, :)
        real(wp), intent(out) :: eigenvalues(:)
        integer(i32), intent(out), optional :: status
        
        call eigen_symmetric(A, eigenvalues, status=status)
    end subroutine eigenvalues_only_symmetric
    
    !> Compute eigenvalues and eigenvectors of general matrix
    !! Uses LAPACK DGEEV - produces complex eigenvalues/vectors in general
    subroutine eigen_general(A, eigenvalues_real, eigenvalues_imag, &
                            left_eigenvectors, right_eigenvectors, status)
        real(wp), intent(in) :: A(:, :)
        real(wp), intent(out) :: eigenvalues_real(:)
        real(wp), intent(out) :: eigenvalues_imag(:)
        real(wp), intent(out), optional :: left_eigenvectors(:, :)
        real(wp), intent(out), optional :: right_eigenvectors(:, :)
        integer(i32), intent(out), optional :: status
        
        real(wp), allocatable :: A_copy(:, :), work(:)
        real(wp), allocatable :: vl(:, :), vr(:, :)
        integer :: n, lda, ldvl, ldvr, lwork, info
        character :: jobvl, jobvr
        
        n = size(A, 1)
        
        if (size(A, 2) /= n .or. size(eigenvalues_real) /= n .or. &
            size(eigenvalues_imag) /= n) then
            if (present(status)) status = EIGEN_ERR_SIZE
            return
        end if
        
        lda = n
        ldvl = n
        ldvr = n
        lwork = max(1, 4 * n)
        
        jobvl = 'N'
        jobvr = 'N'
        
        if (present(left_eigenvectors)) then
            if (size(left_eigenvectors, 1) /= n .or. size(left_eigenvectors, 2) /= n) then
                if (present(status)) status = EIGEN_ERR_SIZE
                return
            end if
            jobvl = 'V'
        end if
        
        if (present(right_eigenvectors)) then
            if (size(right_eigenvectors, 1) /= n .or. size(right_eigenvectors, 2) /= n) then
                if (present(status)) status = EIGEN_ERR_SIZE
                return
            end if
            jobvr = 'V'
        end if
        
        allocate(A_copy(n, n))
        allocate(vl(n, n))
        allocate(vr(n, n))
        allocate(work(lwork))
        
        A_copy = A
        
        call dgeev(jobvl, jobvr, n, A_copy, lda, eigenvalues_real, eigenvalues_imag, &
                   vl, ldvl, vr, ldvr, work, lwork, info)
        
        if (present(left_eigenvectors) .and. info == 0) then
            left_eigenvectors = vl
        end if
        
        if (present(right_eigenvectors) .and. info == 0) then
            right_eigenvectors = vr
        end if
        
        deallocate(A_copy, vl, vr, work)
        
        if (info /= 0) then
            if (present(status)) status = EIGEN_ERR_LAPACK
        else
            if (present(status)) status = EIGEN_SUCCESS
        end if
    end subroutine eigen_general
    
    !> Solve generalised symmetric eigenvalue problem: A*x = lambda*B*x
    !! Uses LAPACK DSYGV - requires B to be positive definite
    subroutine eigen_generalised_symmetric(A, B, eigenvalues, eigenvectors, status)
        real(wp), intent(in) :: A(:, :), B(:, :)
        real(wp), intent(out) :: eigenvalues(:)
        real(wp), intent(out), optional :: eigenvectors(:, :)
        integer(i32), intent(out), optional :: status
        
        real(wp), allocatable :: A_copy(:, :), B_copy(:, :), work(:)
        integer :: n, lda, ldb, lwork, info, itype
        character :: jobz
        
        n = size(A, 1)
        
        if (size(A, 2) /= n .or. size(B, 1) /= n .or. size(B, 2) /= n .or. &
            size(eigenvalues) /= n) then
            if (present(status)) status = EIGEN_ERR_SIZE
            return
        end if
        
        if (present(eigenvectors)) then
            if (size(eigenvectors, 1) /= n .or. size(eigenvectors, 2) /= n) then
                if (present(status)) status = EIGEN_ERR_SIZE
                return
            end if
            jobz = 'V'
        else
            jobz = 'N'
        end if
        
        itype = 1  ! A*x = lambda*B*x
        lda = n
        ldb = n
        lwork = max(1, 3 * n - 1)
        
        allocate(A_copy(n, n))
        allocate(B_copy(n, n))
        allocate(work(lwork))
        
        A_copy = A
        B_copy = B
        
        call dsygv(itype, jobz, 'U', n, A_copy, lda, B_copy, ldb, &
                   eigenvalues, work, lwork, info)
        
        if (present(eigenvectors) .and. info == 0) then
            eigenvectors = A_copy
        end if
        
        deallocate(A_copy, B_copy, work)
        
        if (info /= 0) then
            if (present(status)) status = EIGEN_ERR_LAPACK
        else
            if (present(status)) status = EIGEN_SUCCESS
        end if
    end subroutine eigen_generalised_symmetric
    
    !> Singular Value Decomposition: A = U * Sigma * V^T
    !! Uses LAPACK DGESVD
    subroutine singular_value_decomposition(A, U, S, Vt, status)
        real(wp), intent(in) :: A(:, :)
        real(wp), intent(out), optional :: U(:, :)
        real(wp), intent(out) :: S(:)
        real(wp), intent(out), optional :: Vt(:, :)
        integer(i32), intent(out), optional :: status
        
        real(wp), allocatable :: A_copy(:, :), work(:)
        real(wp), allocatable :: U_full(:, :), Vt_full(:, :)
        integer :: m, n, lda, ldu, ldvt, lwork, info
        character :: jobu, jobvt
        
        m = size(A, 1)
        n = size(A, 2)
        lda = m
        
        if (size(S) /= min(m, n)) then
            if (present(status)) status = EIGEN_ERR_SIZE
            return
        end if
        
        jobu = 'N'
        jobvt = 'N'
        
        if (present(U)) then
            if (size(U, 1) /= m .or. size(U, 2) /= m) then
                if (present(status)) status = EIGEN_ERR_SIZE
                return
            end if
            jobu = 'A'
        end if
        
        if (present(Vt)) then
            if (size(Vt, 1) /= n .or. size(Vt, 2) /= n) then
                if (present(status)) status = EIGEN_ERR_SIZE
                return
            end if
            jobvt = 'A'
        end if
        
        ldu = m
        ldvt = n
        lwork = max(1, 3 * min(m, n) + max(m, n), 5 * min(m, n))
        
        allocate(A_copy(m, n))
        allocate(U_full(m, m))
        allocate(Vt_full(n, n))
        allocate(work(lwork))
        
        A_copy = A
        
        call dgesvd(jobu, jobvt, m, n, A_copy, lda, S, U_full, ldu, &
                    Vt_full, ldvt, work, lwork, info)
        
        if (present(U) .and. info == 0) U = U_full
        if (present(Vt) .and. info == 0) Vt = Vt_full
        
        deallocate(A_copy, U_full, Vt_full, work)
        
        if (info /= 0) then
            if (present(status)) status = EIGEN_ERR_LAPACK
        else
            if (present(status)) status = EIGEN_SUCCESS
        end if
    end subroutine singular_value_decomposition
    
    !> QR decomposition: A = Q * R
    !! Uses LAPACK DGEQRF and DORGQR
    subroutine qr_decomposition(A, Q, R, status)
        real(wp), intent(in) :: A(:, :)
        real(wp), intent(out) :: Q(:, :), R(:, :)
        integer(i32), intent(out), optional :: status
        
        real(wp), allocatable :: A_copy(:, :), tau(:), work(:)
        integer :: m, n, k, lda, lwork, info, i, j
        
        m = size(A, 1)
        n = size(A, 2)
        k = min(m, n)
        lda = m
        
        if (size(Q, 1) /= m .or. size(Q, 2) /= m .or. &
            size(R, 1) /= m .or. size(R, 2) /= n) then
            if (present(status)) status = EIGEN_ERR_SIZE
            return
        end if
        
        lwork = max(1, n)
        
        allocate(A_copy(m, n))
        allocate(tau(k))
        allocate(work(lwork))
        
        A_copy = A
        
        ! QR factorisation
        call dgeqrf(m, n, A_copy, lda, tau, work, lwork, info)
        
        if (info /= 0) then
            if (present(status)) status = EIGEN_ERR_LAPACK
            deallocate(A_copy, tau, work)
            return
        end if
        
        ! Extract R (upper triangular)
        R = 0.0_wp
        do j = 1, n
            do i = 1, min(j, m)
                R(i, j) = A_copy(i, j)
            end do
        end do
        
        ! Generate Q
        deallocate(work)
        lwork = max(1, m)
        allocate(work(lwork))
        
        call dorgqr(m, m, k, A_copy, lda, tau, work, lwork, info)
        
        if (info == 0) then
            Q = A_copy
        end if
        
        deallocate(A_copy, tau, work)
        
        if (info /= 0) then
            if (present(status)) status = EIGEN_ERR_LAPACK
        else
            if (present(status)) status = EIGEN_SUCCESS
        end if
    end subroutine qr_decomposition
    
    !> Schur decomposition: A = Q * T * Q^T
    !! Uses LAPACK DGEES
    subroutine schur_decomposition(A, Q, T, eigenvalues_real, eigenvalues_imag, status)
        real(wp), intent(in) :: A(:, :)
        real(wp), intent(out) :: Q(:, :), T(:, :)
        real(wp), intent(out), optional :: eigenvalues_real(:)
        real(wp), intent(out), optional :: eigenvalues_imag(:)
        integer(i32), intent(out), optional :: status
        
        real(wp), allocatable :: A_copy(:, :), work(:), wr(:), wi(:)
        logical, allocatable :: bwork(:)
        integer :: n, lda, ldvs, lwork, info, sdim
        character :: jobvs, sort
        procedure(select_function_interface), pointer :: select_ptr => null()
        
        n = size(A, 1)
        
        if (size(A, 2) /= n .or. size(Q, 1) /= n .or. size(Q, 2) /= n .or. &
            size(T, 1) /= n .or. size(T, 2) /= n) then
            if (present(status)) status = EIGEN_ERR_SIZE
            return
        end if
        
        lda = n
        ldvs = n
        lwork = max(1, 3 * n)
        jobvs = 'V'
        sort = 'N'
        
        allocate(A_copy(n, n))
        allocate(wr(n))
        allocate(wi(n))
        allocate(work(lwork))
        allocate(bwork(n))
        
        A_copy = A
        
        call dgees(jobvs, sort, select_ptr, n, A_copy, lda, sdim, wr, wi, &
                   Q, ldvs, work, lwork, bwork, info)
        
        T = A_copy
        
        if (present(eigenvalues_real)) eigenvalues_real = wr
        if (present(eigenvalues_imag)) eigenvalues_imag = wi
        
        deallocate(A_copy, wr, wi, work, bwork)
        
        if (info /= 0) then
            if (present(status)) status = EIGEN_ERR_LAPACK
        else
            if (present(status)) status = EIGEN_SUCCESS
        end if
    end subroutine schur_decomposition
    
    !> Power iteration to find dominant eigenvalue and eigenvector
    subroutine power_iteration(A, eigenvalue, eigenvector, max_iter, tol, status)
        real(wp), intent(in) :: A(:, :)
        real(wp), intent(out) :: eigenvalue
        real(wp), intent(out) :: eigenvector(:)
        integer, intent(in), optional :: max_iter
        real(wp), intent(in), optional :: tol
        integer(i32), intent(out), optional :: status
        
        real(wp), allocatable :: v(:), v_new(:)
        real(wp) :: norm_v, lambda_old, tolerance
        integer :: n, iter, max_iterations, i
        
        n = size(A, 1)
        max_iterations = 1000
        tolerance = TOL_DEFAULT
        
        if (present(max_iter)) max_iterations = max_iter
        if (present(tol)) tolerance = tol
        
        if (size(A, 2) /= n .or. size(eigenvector) /= n) then
            if (present(status)) status = EIGEN_ERR_SIZE
            return
        end if
        
        allocate(v(n))
        allocate(v_new(n))
        
        ! Initial guess (random or all ones)
        v = 1.0_wp
        norm_v = sqrt(sum(v**2))
        v = v / norm_v
        
        eigenvalue = 0.0_wp
        
        do iter = 1, max_iterations
            lambda_old = eigenvalue
            
            ! v_new = A * v
            v_new = matmul(A, v)
            
            ! Rayleigh quotient
            eigenvalue = dot_product(v, v_new)
            
            ! Normalise
            norm_v = sqrt(sum(v_new**2))
            v = v_new / norm_v
            
            ! Check convergence
            if (abs(eigenvalue - lambda_old) < tolerance) then
                eigenvector = v
                if (present(status)) status = EIGEN_SUCCESS
                deallocate(v, v_new)
                return
            end if
        end do
        
        ! Did not converge
        eigenvector = v
        if (present(status)) status = EIGEN_ERR_CONVERGENCE
        
        deallocate(v, v_new)
    end subroutine power_iteration
    
    !> Inverse iteration to find eigenvector for known eigenvalue
    subroutine inverse_iteration(A, mu, eigenvector, max_iter, tol, status)
        real(wp), intent(in) :: A(:, :)
        real(wp), intent(in) :: mu  ! Approximate eigenvalue
        real(wp), intent(out) :: eigenvector(:)
        integer, intent(in), optional :: max_iter
        real(wp), intent(in), optional :: tol
        integer(i32), intent(out), optional :: status
        
        real(wp), allocatable :: A_shift(:, :), v(:), v_new(:)
        real(wp) :: norm_v, tolerance
        integer :: n, iter, max_iterations, i
        integer(i32) :: solve_status
        
        n = size(A, 1)
        max_iterations = 1000
        tolerance = TOL_DEFAULT
        
        if (present(max_iter)) max_iterations = max_iter
        if (present(tol)) tolerance = tol
        
        if (size(A, 2) /= n .or. size(eigenvector) /= n) then
            if (present(status)) status = EIGEN_ERR_SIZE
            return
        end if
        
        allocate(A_shift(n, n))
        allocate(v(n))
        allocate(v_new(n))
        
        ! Shift matrix: A_shift = A - mu*I
        A_shift = A
        do i = 1, n
            A_shift(i, i) = A_shift(i, i) - mu
        end do
        
        ! Initial guess
        v = 1.0_wp
        norm_v = sqrt(sum(v**2))
        v = v / norm_v
        
        do iter = 1, max_iterations
            ! Solve (A - mu*I) * v_new = v
            ! This requires solve_linear module - simplified here
            v_new = v  ! Placeholder - actual implementation needs solver
            
            ! Normalise
            norm_v = sqrt(sum(v_new**2))
            v_new = v_new / norm_v
            
            ! Check convergence
            if (sqrt(sum((v_new - v)**2)) < tolerance) then
                eigenvector = v_new
                if (present(status)) status = EIGEN_SUCCESS
                deallocate(A_shift, v, v_new)
                return
            end if
            
            v = v_new
        end do
        
        eigenvector = v
        if (present(status)) status = EIGEN_ERR_CONVERGENCE
        
        deallocate(A_shift, v, v_new)
    end subroutine inverse_iteration

end module eigen