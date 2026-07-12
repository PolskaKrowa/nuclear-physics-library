! fortran/kernels/linear_algebra/sparse_matrix.f90
!
! Sparse matrix operations using CSR (Compressed Sparse Row) format
! Provides memory-efficient storage and operations for sparse matrices
!
module sparse_matrix
    use kinds, only: wp, i32, i64
    use constants, only: TOL_DEFAULT
    implicit none
    private
    
    ! Sparse matrix in CSR format
    type, public :: sparse_matrix_t
        integer :: nrows = 0        ! Number of rows
        integer :: ncols = 0        ! Number of columns
        integer(i64) :: nnz = 0     ! Number of non-zero elements
        real(wp), allocatable :: values(:)      ! Non-zero values
        integer, allocatable :: col_indices(:)  ! Column indices
        integer(i64), allocatable :: row_ptr(:) ! Row pointers
    end type sparse_matrix_t
    
    public :: sparse_create_from_dense, sparse_create_empty
    public :: sparse_destroy, sparse_copy
    public :: sparse_set_value, sparse_get_value
    public :: sparse_matvec, sparse_transpose
    public :: sparse_to_dense, sparse_nnz
    public :: sparse_diagonal, sparse_add_scaled
    public :: sparse_matrix_mult
    
    integer(i32), parameter, public :: SPARSE_SUCCESS = 0
    integer(i32), parameter, public :: SPARSE_ERR_SIZE = 1
    integer(i32), parameter, public :: SPARSE_ERR_ALLOC = 2
    integer(i32), parameter, public :: SPARSE_ERR_INDEX = 3
    
contains

    !> Create sparse matrix from dense matrix
    subroutine sparse_create_from_dense(A_dense, A_sparse, tol, status)
        real(wp), intent(in) :: A_dense(:, :)
        type(sparse_matrix_t), intent(out) :: A_sparse
        real(wp), intent(in), optional :: tol
        integer(i32), intent(out), optional :: status
        
        real(wp) :: tolerance
        integer :: i, j, count
        integer(i64) :: idx
        
        tolerance = TOL_DEFAULT
        if (present(tol)) tolerance = tol
        
        A_sparse%nrows = size(A_dense, 1)
        A_sparse%ncols = size(A_dense, 2)
        
        ! Count non-zeros
        A_sparse%nnz = 0
        do j = 1, A_sparse%ncols
            do i = 1, A_sparse%nrows
                if (abs(A_dense(i, j)) > tolerance) then
                    A_sparse%nnz = A_sparse%nnz + 1
                end if
            end do
        end do
        
        ! Allocate arrays
        allocate(A_sparse%values(A_sparse%nnz))
        allocate(A_sparse%col_indices(A_sparse%nnz))
        allocate(A_sparse%row_ptr(A_sparse%nrows + 1))
        
        ! Fill CSR structure
        idx = 1
        A_sparse%row_ptr(1) = 1
        
        do i = 1, A_sparse%nrows
            count = 0
            do j = 1, A_sparse%ncols
                if (abs(A_dense(i, j)) > tolerance) then
                    A_sparse%values(idx) = A_dense(i, j)
                    A_sparse%col_indices(idx) = j
                    idx = idx + 1
                    count = count + 1
                end if
            end do
            A_sparse%row_ptr(i + 1) = A_sparse%row_ptr(i) + count
        end do
        
        if (present(status)) status = SPARSE_SUCCESS
    end subroutine sparse_create_from_dense
    
    !> Create empty sparse matrix with pre-allocated space
    subroutine sparse_create_empty(A_sparse, nrows, ncols, nnz_estimate, status)
        type(sparse_matrix_t), intent(out) :: A_sparse
        integer, intent(in) :: nrows, ncols
        integer(i64), intent(in) :: nnz_estimate
        integer(i32), intent(out), optional :: status
        
        A_sparse%nrows = nrows
        A_sparse%ncols = ncols
        A_sparse%nnz = 0
        
        allocate(A_sparse%values(nnz_estimate))
        allocate(A_sparse%col_indices(nnz_estimate))
        allocate(A_sparse%row_ptr(nrows + 1))
        
        A_sparse%row_ptr = 1
        
        if (present(status)) status = SPARSE_SUCCESS
    end subroutine sparse_create_empty
    
    !> Deallocate sparse matrix
    subroutine sparse_destroy(A_sparse)
        type(sparse_matrix_t), intent(inout) :: A_sparse
        
        if (allocated(A_sparse%values)) deallocate(A_sparse%values)
        if (allocated(A_sparse%col_indices)) deallocate(A_sparse%col_indices)
        if (allocated(A_sparse%row_ptr)) deallocate(A_sparse%row_ptr)
        
        A_sparse%nrows = 0
        A_sparse%ncols = 0
        A_sparse%nnz = 0
    end subroutine sparse_destroy
    
    !> Copy sparse matrix
    subroutine sparse_copy(A_src, A_dest, status)
        type(sparse_matrix_t), intent(in) :: A_src
        type(sparse_matrix_t), intent(out) :: A_dest
        integer(i32), intent(out), optional :: status
        
        A_dest%nrows = A_src%nrows
        A_dest%ncols = A_src%ncols
        A_dest%nnz = A_src%nnz
        
        allocate(A_dest%values(A_src%nnz))
        allocate(A_dest%col_indices(A_src%nnz))
        allocate(A_dest%row_ptr(A_src%nrows + 1))
        
        A_dest%values = A_src%values
        A_dest%col_indices = A_src%col_indices
        A_dest%row_ptr = A_src%row_ptr
        
        if (present(status)) status = SPARSE_SUCCESS
    end subroutine sparse_copy
    
    !> Set value in sparse matrix (inefficient for many updates)
    subroutine sparse_set_value(A_sparse, row, col, value, status)
        type(sparse_matrix_t), intent(inout) :: A_sparse
        integer, intent(in) :: row, col
        real(wp), intent(in) :: value
        integer(i32), intent(out), optional :: status
        
        integer(i64) :: start_idx, end_idx, i
        logical :: found
        
        if (row < 1 .or. row > A_sparse%nrows .or. &
            col < 1 .or. col > A_sparse%ncols) then
            if (present(status)) status = SPARSE_ERR_INDEX
            return
        end if
        
        start_idx = A_sparse%row_ptr(row)
        end_idx = A_sparse%row_ptr(row + 1) - 1
        
        found = .false.
        do i = start_idx, end_idx
            if (A_sparse%col_indices(i) == col) then
                A_sparse%values(i) = value
                found = .true.
                exit
            end if
        end do
        
        if (.not. found .and. present(status)) then
            status = SPARSE_ERR_INDEX  ! Element not in sparsity pattern
            return
        end if
        
        if (present(status)) status = SPARSE_SUCCESS
    end subroutine sparse_set_value
    
    !> Get value from sparse matrix
    pure function sparse_get_value(A_sparse, row, col) result(value)
        type(sparse_matrix_t), intent(in) :: A_sparse
        integer, intent(in) :: row, col
        real(wp) :: value
        
        integer(i64) :: start_idx, end_idx, i
        
        value = 0.0_wp
        
        if (row < 1 .or. row > A_sparse%nrows .or. &
            col < 1 .or. col > A_sparse%ncols) return
        
        start_idx = A_sparse%row_ptr(row)
        end_idx = A_sparse%row_ptr(row + 1) - 1
        
        do i = start_idx, end_idx
            if (A_sparse%col_indices(i) == col) then
                value = A_sparse%values(i)
                return
            end if
        end do
    end function sparse_get_value
    
    !> Sparse matrix-vector multiplication: y = A * x
    pure subroutine sparse_matvec(A_sparse, x, y)
        type(sparse_matrix_t), intent(in) :: A_sparse
        real(wp), intent(in) :: x(:)
        real(wp), intent(out) :: y(:)
        
        integer :: i
        integer(i64) :: j, start_idx, end_idx
        real(wp) :: sum_val
        
        y = 0.0_wp
        
        do i = 1, A_sparse%nrows
            sum_val = 0.0_wp
            start_idx = A_sparse%row_ptr(i)
            end_idx = A_sparse%row_ptr(i + 1) - 1
            
            do j = start_idx, end_idx
                sum_val = sum_val + A_sparse%values(j) * x(A_sparse%col_indices(j))
            end do
            
            y(i) = sum_val
        end do
    end subroutine sparse_matvec
    
    !> Transpose sparse matrix
    subroutine sparse_transpose(A_sparse, At_sparse, status)
        type(sparse_matrix_t), intent(in) :: A_sparse
        type(sparse_matrix_t), intent(out) :: At_sparse
        integer(i32), intent(out), optional :: status
        
        integer, allocatable :: col_counts(:)
        integer(i64) :: i, j, idx, new_idx
        integer :: row, col
        
        At_sparse%nrows = A_sparse%ncols
        At_sparse%ncols = A_sparse%nrows
        At_sparse%nnz = A_sparse%nnz
        
        allocate(At_sparse%values(At_sparse%nnz))
        allocate(At_sparse%col_indices(At_sparse%nnz))
        allocate(At_sparse%row_ptr(At_sparse%nrows + 1))
        allocate(col_counts(At_sparse%nrows))
        
        ! Count elements per row in transposed matrix
        col_counts = 0
        do i = 1, A_sparse%nnz
            col = A_sparse%col_indices(i)
            col_counts(col) = col_counts(col) + 1
        end do
        
        ! Build row pointers
        At_sparse%row_ptr(1) = 1
        do i = 1, At_sparse%nrows
            At_sparse%row_ptr(i + 1) = At_sparse%row_ptr(i) + col_counts(i)
        end do
        
        ! Fill transposed matrix
        col_counts = 0
        do row = 1, A_sparse%nrows
            do idx = A_sparse%row_ptr(row), A_sparse%row_ptr(row + 1) - 1
                col = A_sparse%col_indices(idx)
                new_idx = At_sparse%row_ptr(col) + col_counts(col)
                
                At_sparse%values(new_idx) = A_sparse%values(idx)
                At_sparse%col_indices(new_idx) = row
                
                col_counts(col) = col_counts(col) + 1
            end do
        end do
        
        deallocate(col_counts)
        
        if (present(status)) status = SPARSE_SUCCESS
    end subroutine sparse_transpose
    
    !> Convert sparse to dense matrix
    pure subroutine sparse_to_dense(A_sparse, A_dense)
        type(sparse_matrix_t), intent(in) :: A_sparse
        real(wp), intent(out) :: A_dense(:, :)
        
        integer :: i
        integer(i64) :: j, start_idx, end_idx
        
        A_dense = 0.0_wp
        
        do i = 1, A_sparse%nrows
            start_idx = A_sparse%row_ptr(i)
            end_idx = A_sparse%row_ptr(i + 1) - 1
            
            do j = start_idx, end_idx
                A_dense(i, A_sparse%col_indices(j)) = A_sparse%values(j)
            end do
        end do
    end subroutine sparse_to_dense
    
    !> Get number of non-zero elements
    pure function sparse_nnz(A_sparse) result(nnz)
        type(sparse_matrix_t), intent(in) :: A_sparse
        integer(i64) :: nnz
        
        nnz = A_sparse%nnz
    end function sparse_nnz
    
    !> Extract diagonal of sparse matrix
    pure subroutine sparse_diagonal(A_sparse, diag)
        type(sparse_matrix_t), intent(in) :: A_sparse
        real(wp), intent(out) :: diag(:)
        
        integer :: i
        integer(i64) :: j, start_idx, end_idx
        
        diag = 0.0_wp
        
        do i = 1, min(A_sparse%nrows, A_sparse%ncols)
            start_idx = A_sparse%row_ptr(i)
            end_idx = A_sparse%row_ptr(i + 1) - 1
            
            do j = start_idx, end_idx
                if (A_sparse%col_indices(j) == i) then
                    diag(i) = A_sparse%values(j)
                    exit
                end if
            end do
        end do
    end subroutine sparse_diagonal
    
    !> Scaled addition: C = alpha * A + beta * B (both sparse)
    !! Note: Only adds elements in A's sparsity pattern
    subroutine sparse_add_scaled(A, B, C, alpha, beta, status)
        type(sparse_matrix_t), intent(in) :: A, B
        type(sparse_matrix_t), intent(out) :: C
        real(wp), intent(in) :: alpha, beta
        integer(i32), intent(out), optional :: status
        
        integer :: i
        integer(i64) :: j, start_idx, end_idx
        
        if (A%nrows /= B%nrows .or. A%ncols /= B%ncols) then
            if (present(status)) status = SPARSE_ERR_SIZE
            return
        end if
        
        ! Copy structure from A
        call sparse_copy(A, C)
        
        ! Scale A's values by alpha
        C%values = alpha * A%values
        
        ! Add beta * B where elements overlap
        do i = 1, A%nrows
            start_idx = A%row_ptr(i)
            end_idx = A%row_ptr(i + 1) - 1
            
            do j = start_idx, end_idx
                C%values(j) = C%values(j) + beta * sparse_get_value(B, i, A%col_indices(j))
            end do
        end do
        
        if (present(status)) status = SPARSE_SUCCESS
    end subroutine sparse_add_scaled
    
    !> Sparse matrix-matrix multiplication: C = A * B
    !! Result may have different sparsity pattern.
    !!
    !! Performance note: implements Gustavson's algorithm with a Sparse
    !! Accumulator (SPA). The previous code allocated a fully dense
    !! A%nrows × B%ncols intermediate — for two 10⁵×10⁵ matrices with
    !! ~10 nnz/row that allocated ~80 GB. We now accumulate each row of C
    !! in a dense workspace of size B%ncols (typically a few KB), then
    !! extract only the non-zeros into C's CSR structure. Memory usage is
    !! O(nnz(C) + B%ncols) instead of O(A%nrows × B%ncols).
    subroutine sparse_matrix_mult(A, B, C, status)
        type(sparse_matrix_t), intent(in) :: A, B
        type(sparse_matrix_t), intent(out) :: C
        integer(i32), intent(out), optional :: status
        
        ! SPA workspace
        real(wp), allocatable :: acc_values(:)      ! dense accumulator, size B%ncols
        integer, allocatable :: acc_marker(:)       ! marks which cols were touched, size B%ncols
        integer, allocatable :: acc_cols(:)         ! list of touched columns per row
        
        ! CSR build buffers. C_row_ptr and row_offset MUST be integer(i64)
        ! to match the sparse_matrix_t%row_ptr field kind (i64) — otherwise
        ! move_alloc fails with a kind mismatch. C_col_indices is integer
        ! (default kind) to match sparse_matrix_t%col_indices.
        integer(i64), allocatable :: C_row_ptr(:)
        real(wp), allocatable :: C_values(:)
        integer, allocatable :: C_col_indices(:)
        integer, allocatable :: row_nnz(:)          ! nnz per row of C
        integer(i64), allocatable :: row_offset(:)  ! prefix sum of row_nnz
        
        ! Loop indices. k/a_start/a_end/b_start/b_end come from A%row_ptr
        ! and B%row_ptr which are integer(i64), so they must be i64 too.
        ! kk_inner is i64 (indexes B%values/col_indices via row_ptr);
        ! kk_extract is plain integer (iterates 1..pos for extraction).
        integer :: i, j, col_b, pos, kk_extract
        integer(i64) :: k, a_start, a_end, b_start, b_end, nnz_C, kk_inner
        integer, parameter :: UNTOUCHED = -1
        
        if (A%ncols /= B%nrows) then
            if (present(status)) status = SPARSE_ERR_SIZE
            return
        end if
        
        ! --- Symbolic pass: count nnz per row of C ---
        allocate(acc_values(B%ncols))
        allocate(acc_marker(B%ncols))
        allocate(acc_cols(B%ncols))
        allocate(row_nnz(A%nrows))
        row_nnz = 0
        acc_marker = UNTOUCHED
        
        do i = 1, A%nrows
            a_start = A%row_ptr(i)
            a_end = A%row_ptr(i + 1) - 1
            row_nnz(i) = 0
            
            do k = a_start, a_end
                j = A%col_indices(k)
                b_start = B%row_ptr(j)
                b_end = B%row_ptr(j + 1) - 1
                
                do kk_inner = b_start, b_end
                    col_b = B%col_indices(kk_inner)
                    if (acc_marker(col_b) /= i) then
                        acc_marker(col_b) = i
                        row_nnz(i) = row_nnz(i) + 1
                    end if
                end do
            end do
        end do
        
        ! --- Allocate CSR arrays via prefix sum ---
        allocate(C_row_ptr(A%nrows + 1))
        allocate(row_offset(A%nrows))
        C_row_ptr(1) = 1
        do i = 1, A%nrows
            C_row_ptr(i + 1) = C_row_ptr(i) + row_nnz(i)
            row_offset(i) = C_row_ptr(i)
        end do
        nnz_C = C_row_ptr(A%nrows + 1) - 1
        
        allocate(C_values(nnz_C))
        allocate(C_col_indices(nnz_C))
        C_values = 0.0_wp
        
        ! --- Numerical pass: accumulate and extract ---
        acc_marker = UNTOUCHED
        acc_values = 0.0_wp
        
        do i = 1, A%nrows
            a_start = A%row_ptr(i)
            a_end = A%row_ptr(i + 1) - 1
            pos = 0  ! number of non-zeros in this row so far
            
            do k = a_start, a_end
                j = A%col_indices(k)
                b_start = B%row_ptr(j)
                b_end = B%row_ptr(j + 1) - 1
                
                do kk_inner = b_start, b_end
                    col_b = B%col_indices(kk_inner)
                    if (acc_marker(col_b) /= i) then
                        acc_marker(col_b) = i
                        pos = pos + 1
                        acc_cols(pos) = col_b
                        acc_values(col_b) = A%values(k) * B%values(kk_inner)
                    else
                        acc_values(col_b) = acc_values(col_b) + A%values(k) * B%values(kk_inner)
                    end if
                end do
            end do
            
            ! Extract non-zeros into C's CSR. Sort column indices for
            ! canonical CSR ordering (required by many downstream solvers).
            call sort_int_real_pair(acc_cols, acc_values, pos)
            do kk_extract = 1, pos
                C_col_indices(row_offset(i) + kk_extract - 1) = acc_cols(kk_extract)
                C_values(row_offset(i) + kk_extract - 1) = acc_values(acc_cols(kk_extract))
            end do
        end do
        
        ! --- Build result ---
        C%nrows = A%nrows
        C%ncols = B%ncols
        C%nnz = nnz_C
        call move_alloc(C_row_ptr, C%row_ptr)
        call move_alloc(C_col_indices, C%col_indices)
        call move_alloc(C_values, C%values)
        
        deallocate(acc_values, acc_marker, acc_cols, row_nnz, row_offset)
        
        if (present(status)) status = SPARSE_SUCCESS
    contains
        ! Simple insertion sort for small column-index lists (typical row
        ! has < 50 non-zeros). Keeps C in canonical sorted-column order.
        subroutine sort_int_real_pair(idx, vals, n)
            integer, intent(inout) :: idx(:)
            real(wp), intent(inout) :: vals(:)
            integer, intent(in) :: n
            integer :: a, b, tmp_idx
            do a = 2, n
                tmp_idx = idx(a)
                b = a - 1
                do while (b >= 1 .and. idx(b) > tmp_idx)
                    idx(b + 1) = idx(b)
                    b = b - 1
                end do
                idx(b + 1) = tmp_idx
            end do
            ! vals is indexed by the column number, so it doesn't need
            ! to be rearranged — only the index list is sorted.
        end subroutine sort_int_real_pair
    end subroutine sparse_matrix_mult

end module sparse_matrix