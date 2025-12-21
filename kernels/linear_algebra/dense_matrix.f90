! fortran/kernels/linear_algebra/dense_matrix.f90
!
! Dense matrix operations using BLAS and LAPACK
! Provides matrix-vector, matrix-matrix operations, and matrix properties
!
module dense_matrix
    use kinds, only: wp, i32
    use constants, only: TOL_DEFAULT
    implicit none
    private
    
    public :: matrix_vector_mult, matrix_matrix_mult
    public :: matrix_transpose, matrix_add, matrix_scale
    public :: matrix_norm, matrix_trace, matrix_determinant
    public :: matrix_rank, matrix_condition_number
    public :: matrix_inverse, matrix_power
    public :: is_symmetric, is_positive_definite
    public :: outer_product, kronecker_product
    
    ! Error codes
    integer(i32), parameter, public :: MAT_SUCCESS = 0
    integer(i32), parameter, public :: MAT_ERR_SIZE = 1
    integer(i32), parameter, public :: MAT_ERR_SINGULAR = 2
    integer(i32), parameter, public :: MAT_ERR_LAPACK = 3
    
    ! BLAS/LAPACK interfaces
    interface
        subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
            use kinds, only: wp
            character, intent(in) :: trans
            integer, intent(in) :: m, n, lda, incx, incy
            real(wp), intent(in) :: alpha, beta
            real(wp), intent(in) :: a(:, :), x(:)
            real(wp), intent(inout) :: y(:)
        end subroutine dgemv
        
        subroutine dgemm(transA, transB, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
            use kinds, only: wp
            character, intent(in) :: transA, transB
            integer, intent(in) :: m, n, k, lda, ldb, ldc
            real(wp), intent(in) :: alpha, beta
            real(wp), intent(in) :: a(:, :), b(:, :)
            real(wp), intent(inout) :: c(:, :)
        end subroutine dgemm
        
        real(wp) function dlange(norm, m, n, a, lda, work)
            use kinds, only: wp
            character, intent(in) :: norm
            integer, intent(in) :: m, n, lda
            real(wp), intent(in) :: a(:, :)
            real(wp), intent(inout) :: work(:)
        end function dlange
        
        subroutine dgetrf(m, n, a, lda, ipiv, info)
            use kinds, only: wp
            integer, intent(in) :: m, n, lda
            real(wp), intent(inout) :: a(:, :)
            integer, intent(out) :: ipiv(:), info
        end subroutine dgetrf
        
        subroutine dgetri(n, a, lda, ipiv, work, lwork, info)
            use kinds, only: wp
            integer, intent(in) :: n, lda, lwork
            real(wp), intent(inout) :: a(:, :), work(:)
            integer, intent(in) :: ipiv(:)
            integer, intent(out) :: info
        end subroutine dgetri
        
        subroutine dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
            use kinds, only: wp
            character, intent(in) :: jobu, jobvt
            integer, intent(in) :: m, n, lda, ldu, ldvt, lwork
            real(wp), intent(inout) :: a(:, :)
            real(wp), intent(out) :: s(:), u(:, :), vt(:, :), work(:)
            integer, intent(out) :: info
        end subroutine dgesvd
        
        subroutine dpotrf(uplo, n, a, lda, info)
            use kinds, only: wp
            character, intent(in) :: uplo
            integer, intent(in) :: n, lda
            real(wp), intent(inout) :: a(:, :)
            integer, intent(out) :: info
        end subroutine dpotrf
    end interface
    
contains

    !> Matrix-vector multiplication: y = alpha * A * x + beta * y
    !! Uses BLAS DGEMV
    subroutine matrix_vector_mult(A, x, y, alpha, beta, trans)
        real(wp), intent(in) :: A(:, :)
        real(wp), intent(in) :: x(:)
        real(wp), intent(inout) :: y(:)
        real(wp), intent(in), optional :: alpha, beta
        logical, intent(in), optional :: trans
        
        real(wp) :: alpha_val, beta_val
        character :: trans_char
        integer :: m, n, lda, incx, incy
        
        m = size(A, 1)
        n = size(A, 2)
        lda = m
        incx = 1
        incy = 1
        
        alpha_val = 1.0_wp
        beta_val = 0.0_wp
        trans_char = 'N'
        
        if (present(alpha)) alpha_val = alpha
        if (present(beta)) beta_val = beta
        if (present(trans)) then
            if (trans) trans_char = 'T'
        end if
        
        call dgemv(trans_char, m, n, alpha_val, A, lda, x, incx, beta_val, y, incy)
    end subroutine matrix_vector_mult
    
    !> Matrix-matrix multiplication: C = alpha * A * B + beta * C
    !! Uses BLAS DGEMM
    subroutine matrix_matrix_mult(A, B, C, alpha, beta, transA, transB)
        real(wp), intent(in) :: A(:, :), B(:, :)
        real(wp), intent(inout) :: C(:, :)
        real(wp), intent(in), optional :: alpha, beta
        logical, intent(in), optional :: transA, transB
        
        real(wp) :: alpha_val, beta_val
        character :: transA_char, transB_char
        integer :: m, n, k, lda, ldb, ldc
        
        alpha_val = 1.0_wp
        beta_val = 0.0_wp
        transA_char = 'N'
        transB_char = 'N'
        
        if (present(alpha)) alpha_val = alpha
        if (present(beta)) beta_val = beta
        if (present(transA)) then
            if (transA) transA_char = 'T'
        end if
        if (present(transB)) then
            if (transB) transB_char = 'T'
        end if
        
        ! Determine dimensions
        if (transA_char == 'N') then
            m = size(A, 1)
            k = size(A, 2)
        else
            m = size(A, 2)
            k = size(A, 1)
        end if
        
        if (transB_char == 'N') then
            n = size(B, 2)
        else
            n = size(B, 1)
        end if
        
        lda = size(A, 1)
        ldb = size(B, 1)
        ldc = size(C, 1)
        
        call dgemm(transA_char, transB_char, m, n, k, alpha_val, &
                   A, lda, B, ldb, beta_val, C, ldc)
    end subroutine matrix_matrix_mult
    
    !> Transpose a matrix
    pure function matrix_transpose(A) result(At)
        real(wp), intent(in) :: A(:, :)
        real(wp) :: At(size(A, 2), size(A, 1))
        integer :: i, j
        
        do j = 1, size(A, 2)
            do i = 1, size(A, 1)
                At(j, i) = A(i, j)
            end do
        end do
    end function matrix_transpose
    
    !> Matrix addition: C = A + B
    pure function matrix_add(A, B) result(C)
        real(wp), intent(in) :: A(:, :), B(:, :)
        real(wp) :: C(size(A, 1), size(A, 2))
        
        C = A + B
    end function matrix_add
    
    !> Matrix scaling: B = alpha * A
    pure function matrix_scale(alpha, A) result(B)
        real(wp), intent(in) :: alpha
        real(wp), intent(in) :: A(:, :)
        real(wp) :: B(size(A, 1), size(A, 2))
        
        B = alpha * A
    end function matrix_scale
    
    !> Compute matrix norm
    !! norm_type: '1' = 1-norm, 'I' = infinity norm, 'F' = Frobenius norm
    function matrix_norm(A, norm_type) result(nrm)
        real(wp), intent(in) :: A(:, :)
        character(len=1), intent(in), optional :: norm_type
        real(wp) :: nrm
        character(len=1) :: norm_char
        integer :: m, n, lda
        real(wp), allocatable :: work(:)
        
        m = size(A, 1)
        n = size(A, 2)
        lda = m
        
        norm_char = 'F'  ! Default: Frobenius norm
        if (present(norm_type)) norm_char = norm_type
        
        if (norm_char == 'I') then
            allocate(work(m))
        else
            allocate(work(1))
        end if
        
        nrm = dlange(norm_char, m, n, A, lda, work)
        
        deallocate(work)
    end function matrix_norm
    
    !> Compute matrix trace (sum of diagonal elements)
    pure function matrix_trace(A) result(tr)
        real(wp), intent(in) :: A(:, :)
        real(wp) :: tr
        integer :: i, n
        
        n = min(size(A, 1), size(A, 2))
        tr = 0.0_wp
        
        do i = 1, n
            tr = tr + A(i, i)
        end do
    end function matrix_trace
    
    !> Compute matrix determinant using LU decomposition
    function matrix_determinant(A, status) result(det)
        real(wp), intent(in) :: A(:, :)
        integer(i32), intent(out), optional :: status
        real(wp) :: det
        
        real(wp), allocatable :: A_copy(:, :)
        integer, allocatable :: ipiv(:)
        integer :: n, info, i, parity
        
        n = size(A, 1)
        
        if (size(A, 2) /= n) then
            det = 0.0_wp
            if (present(status)) status = MAT_ERR_SIZE
            return
        end if
        
        allocate(A_copy(n, n))
        allocate(ipiv(n))
        
        A_copy = A
        
        ! LU factorisation
        call dgetrf(n, n, A_copy, n, ipiv, info)
        
        if (info /= 0) then
            det = 0.0_wp
            if (present(status)) status = MAT_ERR_SINGULAR
            deallocate(A_copy, ipiv)
            return
        end if
        
        ! Determinant is product of diagonal of U, times sign of permutation
        det = 1.0_wp
        parity = 0
        
        do i = 1, n
            det = det * A_copy(i, i)
            if (ipiv(i) /= i) parity = parity + 1
        end do
        
        ! Account for permutation sign
        if (mod(parity, 2) /= 0) det = -det
        
        deallocate(A_copy, ipiv)
        
        if (present(status)) status = MAT_SUCCESS
    end function matrix_determinant
    
    !> Compute matrix rank using SVD
    function matrix_rank(A, tol) result(rnk)
        real(wp), intent(in) :: A(:, :)
        real(wp), intent(in), optional :: tol
        integer :: rnk
        
        real(wp), allocatable :: A_copy(:, :), s(:), u(:, :), vt(:, :), work(:)
        integer :: m, n, lda, ldu, ldvt, lwork, info, i
        real(wp) :: tolerance
        
        m = size(A, 1)
        n = size(A, 2)
        lda = m
        ldu = m
        ldvt = n
        lwork = max(3 * min(m, n) + max(m, n), 5 * min(m, n))
        
        tolerance = TOL_DEFAULT
        if (present(tol)) tolerance = tol
        
        allocate(A_copy(m, n))
        allocate(s(min(m, n)))
        allocate(u(m, m))
        allocate(vt(n, n))
        allocate(work(lwork))
        
        A_copy = A
        
        ! Compute SVD
        call dgesvd('A', 'A', m, n, A_copy, lda, s, u, ldu, vt, ldvt, work, lwork, info)
        
        ! Count singular values above tolerance
        rnk = 0
        do i = 1, min(m, n)
            if (s(i) > tolerance * s(1)) then
                rnk = rnk + 1
            end if
        end do
        
        deallocate(A_copy, s, u, vt, work)
    end function matrix_rank
    
    !> Compute matrix condition number (ratio of largest to smallest singular value)
    function matrix_condition_number(A, norm_type) result(cond)
        real(wp), intent(in) :: A(:, :)
        character(len=1), intent(in), optional :: norm_type
        real(wp) :: cond
        
        real(wp), allocatable :: A_copy(:, :), s(:), work(:)
        integer :: m, n, lda, lwork, info
        character(len=1) :: norm_char
        
        m = size(A, 1)
        n = size(A, 2)
        lda = m
        lwork = max(1, 3 * min(m, n) + max(m, n), 5 * min(m, n))
        
        norm_char = '2'  ! 2-norm (spectral norm)
        if (present(norm_type)) norm_char = norm_type
        
        allocate(A_copy(m, n))
        allocate(s(min(m, n)))
        allocate(work(lwork))
        
        A_copy = A
        
        ! Compute singular values
        call dgesvd('N', 'N', m, n, A_copy, lda, s, A_copy, 1, A_copy, 1, work, lwork, info)
        
        if (info == 0 .and. s(min(m, n)) > TOL_DEFAULT) then
            cond = s(1) / s(min(m, n))
        else
            cond = huge(1.0_wp)
        end if
        
        deallocate(A_copy, s, work)
    end function matrix_condition_number
    
    !> Compute matrix inverse using LU decomposition
    subroutine matrix_inverse(A, Ainv, status)
        real(wp), intent(in) :: A(:, :)
        real(wp), intent(out) :: Ainv(:, :)
        integer(i32), intent(out), optional :: status
        
        integer, allocatable :: ipiv(:)
        real(wp), allocatable :: work(:)
        integer :: n, info, lwork
        
        n = size(A, 1)
        
        if (size(A, 2) /= n .or. size(Ainv, 1) /= n .or. size(Ainv, 2) /= n) then
            if (present(status)) status = MAT_ERR_SIZE
            return
        end if
        
        lwork = n * 64
        allocate(ipiv(n))
        allocate(work(lwork))
        
        Ainv = A
        
        ! LU factorisation
        call dgetrf(n, n, Ainv, n, ipiv, info)
        
        if (info /= 0) then
            if (present(status)) status = MAT_ERR_SINGULAR
            deallocate(ipiv, work)
            return
        end if
        
        ! Invert
        call dgetri(n, Ainv, n, ipiv, work, lwork, info)
        
        deallocate(ipiv, work)
        
        if (info /= 0) then
            if (present(status)) status = MAT_ERR_SINGULAR
        else
            if (present(status)) status = MAT_SUCCESS
        end if
    end subroutine matrix_inverse
    
    !> Compute matrix power A^k using repeated squaring
    function matrix_power(A, k, status) result(Ak)
        real(wp), intent(in) :: A(:, :)
        integer, intent(in) :: k
        integer(i32), intent(out), optional :: status
        real(wp) :: Ak(size(A, 1), size(A, 2))
        
        real(wp) :: temp(size(A, 1), size(A, 2))
        real(wp) :: base(size(A, 1), size(A, 2))
        integer :: n, power, i
        
        n = size(A, 1)
        
        if (size(A, 2) /= n .or. k < 0) then
            Ak = 0.0_wp
            if (present(status)) status = MAT_ERR_SIZE
            return
        end if
        
        if (k == 0) then
            ! Identity matrix
            Ak = 0.0_wp
            do i = 1, n
                Ak(i, i) = 1.0_wp
            end do
            if (present(status)) status = MAT_SUCCESS
            return
        end if
        
        base = A
        Ak = 0.0_wp
        do i = 1, n
            Ak(i, i) = 1.0_wp
        end do
        
        power = k
        
        ! Binary exponentiation
        do while (power > 0)
            if (mod(power, 2) == 1) then
                temp = Ak
                call matrix_matrix_mult(temp, base, Ak)
            end if
            
            if (power > 1) then
                temp = base
                call matrix_matrix_mult(temp, temp, base)
            end if
            
            power = power / 2
        end do
        
        if (present(status)) status = MAT_SUCCESS
    end function matrix_power
    
    !> Check if matrix is symmetric
    pure function is_symmetric(A, tol) result(sym)
        real(wp), intent(in) :: A(:, :)
        real(wp), intent(in), optional :: tol
        logical :: sym
        integer :: n, i, j
        real(wp) :: tolerance
        
        n = size(A, 1)
        sym = .false.
        
        if (size(A, 2) /= n) return
        
        tolerance = TOL_DEFAULT
        if (present(tol)) tolerance = tol
        
        do j = 1, n
            do i = j + 1, n
                if (abs(A(i, j) - A(j, i)) > tolerance) return
            end do
        end do
        
        sym = .true.
    end function is_symmetric
    
    !> Check if matrix is positive definite using Cholesky
    function is_positive_definite(A) result(pd)
        real(wp), intent(in) :: A(:, :)
        logical :: pd
        
        real(wp), allocatable :: A_copy(:, :)
        integer :: n, info
        
        n = size(A, 1)
        pd = .false.
        
        if (size(A, 2) /= n) return
        
        allocate(A_copy(n, n))
        A_copy = A
        
        ! Try Cholesky factorisation
        call dpotrf('L', n, A_copy, n, info)
        
        pd = (info == 0)
        
        deallocate(A_copy)
    end function is_positive_definite
    
    !> Compute outer product: C = a * b^T
    pure function outer_product(a, b) result(C)
        real(wp), intent(in) :: a(:), b(:)
        real(wp) :: C(size(a), size(b))
        integer :: i, j
        
        do j = 1, size(b)
            do i = 1, size(a)
                C(i, j) = a(i) * b(j)
            end do
        end do
    end function outer_product
    
    !> Compute Kronecker product: C = A ⊗ B
    pure function kronecker_product(A, B) result(C)
        real(wp), intent(in) :: A(:, :), B(:, :)
        real(wp) :: C(size(A, 1) * size(B, 1), size(A, 2) * size(B, 2))
        integer :: i, j, k, l, ma, na, mb, nb
        integer :: row_offset, col_offset
        
        ma = size(A, 1)
        na = size(A, 2)
        mb = size(B, 1)
        nb = size(B, 2)
        
        do j = 1, na
            do i = 1, ma
                row_offset = (i - 1) * mb
                col_offset = (j - 1) * nb
                
                do l = 1, nb
                    do k = 1, mb
                        C(row_offset + k, col_offset + l) = A(i, j) * B(k, l)
                    end do
                end do
            end do
        end do
    end function kronecker_product

end module dense_matrix