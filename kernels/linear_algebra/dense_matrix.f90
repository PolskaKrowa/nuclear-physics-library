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
    
contains

    !> Matrix-vector multiplication: y = alpha * A * x + beta * y
    subroutine matrix_vector_mult(A, x, y, alpha, beta, trans)
        real(wp), intent(in) :: A(:, :)
        real(wp), intent(in) :: x(:)
        real(wp), intent(inout) :: y(:)
        real(wp), intent(in), optional :: alpha, beta
        logical, intent(in), optional :: trans
        
        real(wp) :: alpha_val, beta_val
        character :: trans_char
        integer :: m, n
        real(wp), allocatable :: A_work(:, :), x_work(:), y_work(:)
        
        external :: dgemv
        
        m = size(A, 1)
        n = size(A, 2)
        
        alpha_val = 1.0_wp
        beta_val = 0.0_wp
        trans_char = 'N'
        
        if (present(alpha)) alpha_val = alpha
        if (present(beta)) beta_val = beta
        if (present(trans)) then
            if (trans) trans_char = 'T'
        end if
        
        ! Create contiguous working copies
        allocate(A_work(m, n))
        allocate(x_work(size(x)))
        allocate(y_work(size(y)))
        
        A_work = A
        x_work = x
        y_work = y
        
        ! Call BLAS dgemv
        call dgemv(trans_char, m, n, alpha_val, A_work, m, x_work, 1, beta_val, y_work, 1)
        
        ! Copy result back
        y = y_work
        
        deallocate(A_work, x_work, y_work)
    end subroutine matrix_vector_mult
    
    !> Matrix-matrix multiplication: C = alpha * A * B + beta * C
    subroutine matrix_matrix_mult(A, B, C, alpha, beta, transA, transB)
        real(wp), intent(in) :: A(:, :), B(:, :)
        real(wp), intent(inout) :: C(:, :)
        real(wp), intent(in), optional :: alpha, beta
        logical, intent(in), optional :: transA, transB
        
        real(wp) :: alpha_val, beta_val
        character :: transA_char, transB_char
        integer :: m, n, k, lda, ldb, ldc
        real(wp), allocatable :: A_work(:, :), B_work(:, :), C_work(:, :)
        
        external :: dgemm
        
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
        
        ! Create contiguous working copies
        allocate(A_work(size(A, 1), size(A, 2)))
        allocate(B_work(size(B, 1), size(B, 2)))
        allocate(C_work(size(C, 1), size(C, 2)))
        
        A_work = A
        B_work = B
        C_work = C
        
        call dgemm(transA_char, transB_char, m, n, k, alpha_val, &
                   A_work, lda, B_work, ldb, beta_val, C_work, ldc)
        
        ! Copy result back
        C = C_work
        
        deallocate(A_work, B_work, C_work)
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
    function matrix_norm(A, norm_type) result(nrm)
        real(wp), intent(in) :: A(:, :)
        character(len=1), intent(in), optional :: norm_type
        real(wp) :: nrm
        character(len=1) :: norm_char
        integer :: m, n
        real(wp), allocatable :: work(:), A_work(:, :)
        
        double precision, external :: dlange
        
        m = size(A, 1)
        n = size(A, 2)
        
        norm_char = 'F'  ! Default: Frobenius norm
        if (present(norm_type)) norm_char = norm_type
        
        allocate(A_work(m, n))
        A_work = A
        
        if (norm_char == 'I' .or. norm_char == 'i') then
            allocate(work(m))
        else
            allocate(work(1))
        end if
        
        nrm = dlange(norm_char, m, n, A_work, m, work)
        
        deallocate(work, A_work)
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
        
        external :: dgetrf
        
        n = size(A, 1)
        
        if (size(A, 2) /= n) then
            det = 0.0_wp
            if (present(status)) status = MAT_ERR_SIZE
            return
        end if
        
        allocate(A_copy(n, n))
        allocate(ipiv(n))
        
        A_copy = A
        
        ! LU factorization with LAPACK
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
        integer :: m, n, lwork, info, i
        real(wp) :: tolerance
        
        external :: dgesvd
        
        m = size(A, 1)
        n = size(A, 2)
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
        call dgesvd('A', 'A', m, n, A_copy, m, s, u, m, vt, n, work, lwork, info)
        
        ! Count singular values above tolerance
        rnk = 0
        if (info == 0) then
            do i = 1, min(m, n)
                if (s(i) > tolerance * s(1)) then
                    rnk = rnk + 1
                end if
            end do
        end if
        
        deallocate(A_copy, s, u, vt, work)
    end function matrix_rank
    
    !> Compute matrix condition number
    function matrix_condition_number(A, norm_type) result(cond)
        real(wp), intent(in) :: A(:, :)
        character(len=1), intent(in), optional :: norm_type
        real(wp) :: cond
        
        real(wp), allocatable :: A_copy(:, :), s(:), work(:), u(:, :), vt(:, :)
        integer :: m, n, lwork, info
        
        external :: dgesvd
        
        m = size(A, 1)
        n = size(A, 2)
        lwork = max(1, 3 * min(m, n) + max(m, n), 5 * min(m, n))
        
        allocate(A_copy(m, n))
        allocate(s(min(m, n)))
        allocate(u(1, 1))  ! Not computed
        allocate(vt(1, 1)) ! Not computed
        allocate(work(lwork))
        
        A_copy = A
        
        ! Compute singular values only
        call dgesvd('N', 'N', m, n, A_copy, m, s, u, 1, vt, 1, work, lwork, info)
        
        if (info == 0 .and. s(min(m, n)) > TOL_DEFAULT) then
            cond = s(1) / s(min(m, n))
        else
            cond = huge(1.0_wp)
        end if
        
        deallocate(A_copy, s, u, vt, work)
    end function matrix_condition_number
    
    !> Compute matrix inverse using LU decomposition
    subroutine matrix_inverse(A, Ainv, status)
        real(wp), intent(in) :: A(:, :)
        real(wp), intent(out) :: Ainv(:, :)
        integer(i32), intent(out), optional :: status
        
        integer, allocatable :: ipiv(:)
        real(wp), allocatable :: work(:)
        integer :: n, info, lwork
        
        external :: dgetrf, dgetri
        
        n = size(A, 1)
        
        if (size(A, 2) /= n .or. size(Ainv, 1) /= n .or. size(Ainv, 2) /= n) then
            if (present(status)) status = MAT_ERR_SIZE
            return
        end if
        
        lwork = n * 64
        allocate(ipiv(n))
        allocate(work(lwork))
        
        Ainv = A
        
        ! LU factorization
        call dgetrf(n, n, Ainv, n, ipiv, info)
        
        if (info /= 0) then
            if (present(status)) status = MAT_ERR_SINGULAR
            deallocate(ipiv, work)
            return
        end if
        
        ! Invert using LU
        call dgetri(n, Ainv, n, ipiv, work, lwork, info)
        
        deallocate(ipiv, work)
        
        if (info /= 0) then
            if (present(status)) status = MAT_ERR_SINGULAR
        else
            if (present(status)) status = MAT_SUCCESS
        end if
    end subroutine matrix_inverse
    
    !> Compute matrix power A^k
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
    
    !> Check if matrix is positive definite
    function is_positive_definite(A) result(pd)
        real(wp), intent(in) :: A(:, :)
        logical :: pd
        
        real(wp), allocatable :: A_copy(:, :)
        integer :: n, info
        
        external :: dpotrf
        
        n = size(A, 1)
        pd = .false.
        
        if (size(A, 2) /= n) return
        
        allocate(A_copy(n, n))
        A_copy = A
        
        ! Try Cholesky factorization
        call dpotrf('L', n, A_copy, n, info)
        
        pd = (info == 0)
        
        deallocate(A_copy)
    end function is_positive_definite
    
    !> Compute outer product
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
    
    !> Compute Kronecker product
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