! fortran/kernels/linear_algebra/solve_linear.f90
!
! Linear system solvers using LAPACK
! Solves systems of the form Ax = b for various matrix types
!
! Supported systems:
!   - General dense systems (LU decomposition)
!   - Symmetric positive definite (Cholesky)
!   - Symmetric indefinite (Bunch-Kaufman)
!   - Tridiagonal systems (Thomas algorithm)
!
module solve_linear
    use kinds, only: wp, i32
    use constants, only: TOL_DEFAULT
    implicit none
    private
    
    public :: solve_dense, solve_spd, solve_symmetric
    public :: solve_tridiagonal, solve_least_squares
    public :: lu_decompose, lu_solve
    public :: cholesky_decompose, cholesky_solve
    
    ! Error codes
    integer(i32), parameter, public :: SOLVE_SUCCESS = 0
    integer(i32), parameter, public :: SOLVE_ERR_SINGULAR = 1
    integer(i32), parameter, public :: SOLVE_ERR_NOT_SPD = 2
    integer(i32), parameter, public :: SOLVE_ERR_INVALID_SIZE = 3
    integer(i32), parameter, public :: SOLVE_ERR_LAPACK = 4
    
contains

    !> Solve general dense linear system Ax = b using LU decomposition
    !! Uses LAPACK DGESV routine
    subroutine solve_dense(A, b, x, status)
        real(wp), intent(in) :: A(:, :)  ! n x n matrix
        real(wp), intent(in) :: b(:)     ! Right-hand side
        real(wp), intent(out) :: x(:)    ! Solution vector
        integer(i32), intent(out), optional :: status
        
        real(wp), allocatable :: A_copy(:, :)
        integer :: n, nrhs, lda, ldb, info
        integer, allocatable :: ipiv(:)
        
        n = size(A, 1)
        nrhs = 1
        lda = n
        ldb = n
        
        ! Check dimensions
        if (size(A, 2) /= n .or. size(b) /= n .or. size(x) /= n) then
            if (present(status)) status = SOLVE_ERR_INVALID_SIZE
            return
        end if
        
        ! Copy input (LAPACK overwrites)
        allocate(A_copy(n, n))
        allocate(ipiv(n))
        A_copy = A
        x = b
        
        ! Solve using LAPACK DGESV (LU factorisation)
        call dgesv(n, nrhs, A_copy, lda, ipiv, x, ldb, info)
        
        deallocate(A_copy, ipiv)
        
        ! Check for errors
        if (info < 0) then
            if (present(status)) status = SOLVE_ERR_LAPACK
        else if (info > 0) then
            if (present(status)) status = SOLVE_ERR_SINGULAR
        else
            if (present(status)) status = SOLVE_SUCCESS
        end if
    end subroutine solve_dense
    
    !> Solve symmetric positive definite system using Cholesky decomposition
    !! Uses LAPACK DPOSV routine
    subroutine solve_spd(A, b, x, status)
        real(wp), intent(in) :: A(:, :)
        real(wp), intent(in) :: b(:)
        real(wp), intent(out) :: x(:)
        integer(i32), intent(out), optional :: status
        
        real(wp), allocatable :: A_copy(:, :)
        integer :: n, nrhs, lda, ldb, info
        character :: uplo
        
        n = size(A, 1)
        nrhs = 1
        lda = n
        ldb = n
        uplo = 'U'  ! Upper triangular
        
        if (size(A, 2) /= n .or. size(b) /= n .or. size(x) /= n) then
            if (present(status)) status = SOLVE_ERR_INVALID_SIZE
            return
        end if
        
        allocate(A_copy(n, n))
        A_copy = A
        x = b
        
        ! Solve using Cholesky (DPOSV)
        call dposv(uplo, n, nrhs, A_copy, lda, x, ldb, info)
        
        deallocate(A_copy)
        
        if (info < 0) then
            if (present(status)) status = SOLVE_ERR_LAPACK
        else if (info > 0) then
            if (present(status)) status = SOLVE_ERR_NOT_SPD
        else
            if (present(status)) status = SOLVE_SUCCESS
        end if
    end subroutine solve_spd
    
    !> Solve symmetric indefinite system using Bunch-Kaufman decomposition
    !! Uses LAPACK DSYSV routine
    subroutine solve_symmetric(A, b, x, status)
        real(wp), intent(in) :: A(:, :)
        real(wp), intent(in) :: b(:)
        real(wp), intent(out) :: x(:)
        integer(i32), intent(out), optional :: status
        
        real(wp), allocatable :: A_copy(:, :), work(:)
        integer, allocatable :: ipiv(:)
        integer :: n, nrhs, lda, ldb, lwork, info
        character :: uplo
        
        n = size(A, 1)
        nrhs = 1
        lda = n
        ldb = n
        lwork = n * 64  ! Optimal work array size
        uplo = 'U'
        
        if (size(A, 2) /= n .or. size(b) /= n .or. size(x) /= n) then
            if (present(status)) status = SOLVE_ERR_INVALID_SIZE
            return
        end if
        
        allocate(A_copy(n, n))
        allocate(ipiv(n))
        allocate(work(lwork))
        
        A_copy = A
        x = b
        
        call dsysv(uplo, n, nrhs, A_copy, lda, ipiv, x, ldb, work, lwork, info)
        
        deallocate(A_copy, ipiv, work)
        
        if (info < 0) then
            if (present(status)) status = SOLVE_ERR_LAPACK
        else if (info > 0) then
            if (present(status)) status = SOLVE_ERR_SINGULAR
        else
            if (present(status)) status = SOLVE_SUCCESS
        end if
    end subroutine solve_symmetric
    
    !> Solve tridiagonal system using Thomas algorithm
    !! More efficient than general solver for tridiagonal matrices
    pure subroutine solve_tridiagonal(a, b, c, d, x, status)
        real(wp), intent(in) :: a(:)   ! Lower diagonal (size n-1)
        real(wp), intent(in) :: b(:)   ! Main diagonal (size n)
        real(wp), intent(in) :: c(:)   ! Upper diagonal (size n-1)
        real(wp), intent(in) :: d(:)   ! Right-hand side (size n)
        real(wp), intent(out) :: x(:)  ! Solution (size n)
        integer(i32), intent(out), optional :: status
        
        real(wp) :: cp(size(b) - 1), dp(size(b))
        real(wp) :: m
        integer :: n, i
        
        n = size(b)
        
        if (size(a) /= n - 1 .or. size(c) /= n - 1 .or. &
            size(d) /= n .or. size(x) /= n) then
            if (present(status)) status = SOLVE_ERR_INVALID_SIZE
            return
        end if
        
        ! Forward sweep
        cp(1) = c(1) / b(1)
        dp(1) = d(1) / b(1)
        
        do i = 2, n - 1
            m = b(i) - a(i - 1) * cp(i - 1)
            if (abs(m) < TOL_DEFAULT) then
                if (present(status)) status = SOLVE_ERR_SINGULAR
                return
            end if
            cp(i) = c(i) / m
            dp(i) = (d(i) - a(i - 1) * dp(i - 1)) / m
        end do
        
        m = b(n) - a(n - 1) * cp(n - 1)
        if (abs(m) < TOL_DEFAULT) then
            if (present(status)) status = SOLVE_ERR_SINGULAR
            return
        end if
        dp(n) = (d(n) - a(n - 1) * dp(n - 1)) / m
        
        ! Back substitution
        x(n) = dp(n)
        do i = n - 1, 1, -1
            x(i) = dp(i) - cp(i) * x(i + 1)
        end do
        
        if (present(status)) status = SOLVE_SUCCESS
    end subroutine solve_tridiagonal
    
    !> Solve overdetermined system in least-squares sense
    !! Minimises ||Ax - b||₂ using QR decomposition (DGELS)
    subroutine solve_least_squares(A, b, x, status)
        real(wp), intent(in) :: A(:, :)   ! m x n matrix (m >= n)
        real(wp), intent(in) :: b(:)      ! Size m
        real(wp), intent(out) :: x(:)     ! Size n
        integer(i32), intent(out), optional :: status
        
        real(wp), allocatable :: A_copy(:, :), b_copy(:), work(:)
        integer :: m, n, nrhs, lda, ldb, lwork, info
        character :: trans
        
        m = size(A, 1)
        n = size(A, 2)
        nrhs = 1
        lda = m
        ldb = max(m, n)
        lwork = max(1, min(m, n) + max(min(m, n), nrhs))
        lwork = lwork * 64  ! Optimal
        trans = 'N'
        
        if (m < n .or. size(b) /= m .or. size(x) /= n) then
            if (present(status)) status = SOLVE_ERR_INVALID_SIZE
            return
        end if
        
        allocate(A_copy(m, n))
        allocate(b_copy(ldb))
        allocate(work(lwork))
        
        A_copy = A
        b_copy(1:m) = b
        if (ldb > m) b_copy(m + 1:ldb) = 0.0_wp
        
        call dgels(trans, m, n, nrhs, A_copy, lda, b_copy, ldb, work, lwork, info)
        
        x = b_copy(1:n)
        
        deallocate(A_copy, b_copy, work)
        
        if (info /= 0) then
            if (present(status)) status = SOLVE_ERR_LAPACK
        else
            if (present(status)) status = SOLVE_SUCCESS
        end if
    end subroutine solve_least_squares
    
    !> Compute LU decomposition of matrix A
    subroutine lu_decompose(A, L, U, P, status)
        real(wp), intent(in) :: A(:, :)
        real(wp), intent(out) :: L(:, :), U(:, :)
        integer, intent(out), optional :: P(:)
        integer(i32), intent(out), optional :: status
        
        real(wp), allocatable :: A_copy(:, :)
        integer, allocatable :: ipiv(:)
        integer :: n, info, i, j
        
        n = size(A, 1)
        
        if (size(A, 2) /= n .or. size(L, 1) /= n .or. size(L, 2) /= n .or. &
            size(U, 1) /= n .or. size(U, 2) /= n) then
            if (present(status)) status = SOLVE_ERR_INVALID_SIZE
            return
        end if
        
        allocate(A_copy(n, n))
        allocate(ipiv(n))
        A_copy = A
        
        call dgetrf(n, n, A_copy, n, ipiv, info)
        
        ! Extract L and U
        L = 0.0_wp
        U = 0.0_wp
        
        do j = 1, n
            do i = 1, n
                if (i > j) then
                    L(i, j) = A_copy(i, j)
                else if (i == j) then
                    L(i, j) = 1.0_wp
                    U(i, j) = A_copy(i, j)
                else
                    U(i, j) = A_copy(i, j)
                end if
            end do
        end do
        
        if (present(P)) then
            if (size(P) /= n) then
                if (present(status)) status = SOLVE_ERR_INVALID_SIZE
                deallocate(A_copy, ipiv)
                return
            end if
            P = ipiv
        end if
        
        deallocate(A_copy, ipiv)
        
        if (info /= 0) then
            if (present(status)) status = SOLVE_ERR_SINGULAR
        else
            if (present(status)) status = SOLVE_SUCCESS
        end if
    end subroutine lu_decompose
    
    !> Solve system using precomputed LU decomposition
    pure subroutine lu_solve(L, U, P, b, x, status)
        real(wp), intent(in) :: L(:, :), U(:, :)
        integer, intent(in) :: P(:)
        real(wp), intent(in) :: b(:)
        real(wp), intent(out) :: x(:)
        integer(i32), intent(out), optional :: status
        
        real(wp) :: y(size(b))
        integer :: n, i, j
        real(wp) :: sum_val
        
        n = size(b)
        
        ! Apply permutation
        do i = 1, n
            y(i) = b(P(i))
        end do
        
        ! Forward substitution: Ly = Pb
        do i = 1, n
            sum_val = y(i)
            do j = 1, i - 1
                sum_val = sum_val - L(i, j) * y(j)
            end do
            y(i) = sum_val / L(i, i)
        end do
        
        ! Back substitution: Ux = y
        do i = n, 1, -1
            sum_val = y(i)
            do j = i + 1, n
                sum_val = sum_val - U(i, j) * x(j)
            end do
            
            if (abs(U(i, i)) < TOL_DEFAULT) then
                if (present(status)) status = SOLVE_ERR_SINGULAR
                return
            end if
            
            x(i) = sum_val / U(i, i)
        end do
        
        if (present(status)) status = SOLVE_SUCCESS
    end subroutine lu_solve
    
    !> Compute Cholesky decomposition A = L * L^T
    subroutine cholesky_decompose(A, L, status)
        real(wp), intent(in) :: A(:, :)
        real(wp), intent(out) :: L(:, :)
        integer(i32), intent(out), optional :: status
        
        real(wp), allocatable :: A_copy(:, :)
        integer :: n, info, i, j
        character :: uplo
        
        n = size(A, 1)
        uplo = 'L'  ! Lower triangular
        
        if (size(A, 2) /= n .or. size(L, 1) /= n .or. size(L, 2) /= n) then
            if (present(status)) status = SOLVE_ERR_INVALID_SIZE
            return
        end if
        
        allocate(A_copy(n, n))
        A_copy = A
        
        call dpotrf(uplo, n, A_copy, n, info)
        
        ! Extract lower triangular part
        L = 0.0_wp
        do j = 1, n
            do i = j, n
                L(i, j) = A_copy(i, j)
            end do
        end do
        
        deallocate(A_copy)
        
        if (info < 0) then
            if (present(status)) status = SOLVE_ERR_LAPACK
        else if (info > 0) then
            if (present(status)) status = SOLVE_ERR_NOT_SPD
        else
            if (present(status)) status = SOLVE_SUCCESS
        end if
    end subroutine cholesky_decompose
    
    !> Solve system using precomputed Cholesky decomposition
    pure subroutine cholesky_solve(L, b, x, status)
        real(wp), intent(in) :: L(:, :)
        real(wp), intent(in) :: b(:)
        real(wp), intent(out) :: x(:)
        integer(i32), intent(out), optional :: status
        
        real(wp) :: y(size(b))
        integer :: n, i, j
        real(wp) :: sum_val
        
        n = size(b)
        
        ! Forward substitution: Ly = b
        do i = 1, n
            sum_val = b(i)
            do j = 1, i - 1
                sum_val = sum_val - L(i, j) * y(j)
            end do
            
            if (abs(L(i, i)) < TOL_DEFAULT) then
                if (present(status)) status = SOLVE_ERR_SINGULAR
                return
            end if
            
            y(i) = sum_val / L(i, i)
        end do
        
        ! Back substitution: L^T x = y
        do i = n, 1, -1
            sum_val = y(i)
            do j = i + 1, n
                sum_val = sum_val - L(j, i) * x(j)
            end do
            x(i) = sum_val / L(i, i)
        end do
        
        if (present(status)) status = SOLVE_SUCCESS
    end subroutine cholesky_solve

end module solve_linear