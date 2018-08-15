module array_utilities
!
!!
!!    Array utilities module
!!    Written by Sarai D. Folkstad and Eirik F. Kjønstad, May 2018
!!
!!    This module contains routines that perform various operations on arrays
!!    and that do not belong elsewhere (such as in the reordering module,
!!    where all such routines are gathered for convenience).
!!
!
   use index
   use kinds
   use memory_manager_class
   use dpstrf_eT
!
   implicit none
!
contains
!
!
   integer(i15) function get_max_index(x, dim)
!!
!!    Get max index
!!    Written by Eirik F. Kjønstad, 2018
!!
      implicit none
!
      integer(i15), intent(in) :: dim
      real(dp), dimension(dim,1), intent(in) :: x
!
      integer(i15) :: I
      real(dp)     :: maxval
!
      get_max_index = 1
      maxval = x(1,1)
      do I = 2, dim
!
         if (x(I,1) .gt. maxval) then
!
            get_max_index = I
            maxval = x(I,1)
!
         endif
!
      enddo
!
   end function get_max_index
!
!
   real(dp) function dot_product(x, y, n)
!!
!!    Calculate dot product
!!    Written by Eirik F. Kjønstad, June 2018
!!
!!    Returns the dot product of x and y, two vectors of length n
!!
      implicit none
!
      integer(i15), intent(in) :: n
!
      real(dp), dimension(:,:), intent(in) :: x
      real(dp), dimension(:,:), intent(in) :: y
!
      real(dp) :: ddot
!
      dot_product = ddot(n, x, 1, y, 1)
!
   end function dot_product
!
!
   logical function is_significant(vec, dim, threshold, screening)
!!
!!    Is vector significant ?
!!    Written by Eirik F. Kjønstad and Sarai D. Folkstad, June 2018
!!
!!    Returns true if all elements are below threshold
!!
      implicit none
!
      integer(i15), intent(in) :: dim
!
      real(dp), dimension(dim,1), intent(in)  :: vec
      real(dp), dimension(dim,1), intent(in), optional  :: screening
!
      real(dp), intent(in)  :: threshold
!
      integer(i15) :: i = 0
!
      is_significant = .false.
!
      if (present(screening)) then
         do i = 1, dim
!
            if (abs(vec(i, 1)*screening(i, 1)) .gt. threshold) then
!
               is_significant = .true.
               return
!
            endif
!
         enddo
      else
!
         do i = 1, dim
!
            if (abs(vec(i, 1)) .gt. threshold) then
!
               is_significant = .true.
               return
!
            endif
!
         enddo
      endif
!
   end function is_significant
!
!
   integer(i15) function n_significant(vec, dim, threshold)
!!
!!    Number of significant in vector
!!    Written by Eirik F. Kjønstad and Sarai D. Folkstad, June 2018
!!
!!    Returns the number of elements in vector larger than threshold
!!
      implicit none
!
      integer(i15), intent(in) :: dim
!
      real(dp), dimension(dim,1), intent(in)  :: vec
!
      real(dp), intent(in)  :: threshold
!
      integer(i15) :: i = 0
!
      n_significant = 0
!
      do i = 1, dim
!
         if (abs(vec(i, 1)) .gt. threshold) then
!
            n_significant = n_significant + 1
!
         endif
!
      enddo
!
   end function n_significant
!
!
   subroutine reduce_vector(vec, vec_reduced, block_firsts, block_significant, n_blocks, dim, dim_reduced)
!!
!!    Reduce vector
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Cuts the significant blocks out of a vector and places them in a reduced size
!!    vector
!!
      implicit none
!
      integer(i15) :: dim, dim_reduced, n_blocks
!
      logical, dimension(n_blocks, 1) :: block_significant
      integer(i15), dimension(n_blocks + 1, 1) :: block_firsts
!
      real(dp), dimension(dim, 1) :: vec
      real(dp), dimension(dim_reduced, 1) :: vec_reduced
!
      integer(i15) :: block, current_pos, first, last, size
!
      current_pos = 1
!
      do block = 1, n_blocks
!
         if (block_significant(block, 1)) then
!
            first = block_firsts(block, 1)
            last  = block_firsts(block + 1, 1) - 1
            size  = last - first + 1
!
            vec_reduced(current_pos : current_pos + size - 1, 1) = vec(first : last, 1)
            current_pos = current_pos + size
!
         endif
!
      enddo
!
   end subroutine reduce_vector
!
!
   subroutine reduce_vector_int(vec, vec_reduced, block_firsts, block_significant, n_blocks, dim, dim_reduced)
!!
!!    Reduce vector
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Cuts the significant blocks out of a vector and places them in a reduced size
!!    vector
!!
      implicit none
!
      integer(i15) :: dim, dim_reduced, n_blocks
!
      logical, dimension(n_blocks, 1) :: block_significant
      integer(i15), dimension(n_blocks + 1, 1) :: block_firsts
!
      integer(i15), dimension(dim, 1) :: vec
      integer(i15), dimension(dim_reduced, 1) :: vec_reduced
!
      integer(i15) :: block, current_pos, first, last, size
!
      current_pos = 1
!
      do block = 1, n_blocks
!
         if (block_significant(block, 1)) then
!
            first = block_firsts(block, 1)
            last  = block_firsts(block + 1, 1) - 1
            size  = last - first + 1
!
            vec_reduced(current_pos : current_pos + size - 1, 1) = vec(first : last, 1)
            current_pos = current_pos + size
!
         endif
!
      enddo
!
   end subroutine reduce_vector_int
!
!
   subroutine reduce_array(array, array_reduced, block_firsts, block_significant, n_blocks, dim, dim_reduced, columns)
!!
!!    Reduce array
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Cuts the significant row blocks out of a array and places them in a reduced size
!!    array
!!
      implicit none
!
      integer(i15) :: dim, dim_reduced, n_blocks, columns
!
      logical, dimension(n_blocks, 1) :: block_significant
      integer(i15), dimension(n_blocks + 1, 1) :: block_firsts
!
      real(dp), dimension(dim, columns) :: array
      real(dp), dimension(dim_reduced, columns) :: array_reduced
!
      integer(i15) :: block, current_pos, first, last, size, I
!
!$omp parallel do schedule(static) private(I, current_pos, block, first, last, size)
      do I = 1, columns
!
         current_pos = 1
!
         do block = 1, n_blocks
!
            if (block_significant(block, 1)) then
!
               first = block_firsts(block, 1)
               last  = block_firsts(block + 1, 1) - 1
               size  = last - first + 1
!
               array_reduced(current_pos : current_pos + size - 1, I) = array(first : last, I)
!
               current_pos = current_pos + size
!
            endif
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine reduce_array
!
!
   subroutine reduce_array_column(array, array_reduced, block_firsts, block_significant, n_blocks, dim, dim_reduced, rows)
!!
!!    Reduce array column
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Cuts the significant column blocks out of a array and places them in a reduced size array
!!
      implicit none
!
      integer(i15) :: dim, dim_reduced, n_blocks, rows
!
      logical, dimension(n_blocks, 1) :: block_significant
      integer(i15), dimension(n_blocks + 1, 1) :: block_firsts
!
      real(dp), dimension(rows, dim) :: array
      real(dp), dimension(rows, dim_reduced) :: array_reduced
!
      integer(i15) :: block, current_pos, first, last, size
!
      current_pos = 1
!
      do block = 1, n_blocks
!
         if (block_significant(block, 1)) then
!
            first = block_firsts(block, 1)
            last  = block_firsts(block + 1, 1) - 1
            size  = last - first + 1
!
            array_reduced(:, current_pos : current_pos + size - 1) = array(:, first : last)
         !   array_reduced(current_pos : current_pos + size - 1, 1 : columns) = array(first : last, 1 : columns)
!
            current_pos = current_pos + size
!
         endif
!
      enddo
!
   end subroutine reduce_array_column
!
!
   subroutine reduce_array_int(array, array_reduced, block_firsts, block_significant, n_blocks, dim, dim_reduced, columns)
!!
!!    Reduce array
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Cuts the significant row blocks out of a array and places them in a reduced size
!!    array
!!
      implicit none
!
      integer(i15) :: dim, dim_reduced, n_blocks, columns
!
      logical, dimension(n_blocks, 1) :: block_significant
      integer(i15), dimension(n_blocks + 1, 1) :: block_firsts
!
      integer(i15), dimension(dim, columns) :: array
      integer(i15), dimension(dim_reduced, columns) :: array_reduced
!
      integer(i15) :: block, current_pos, first, last, size
!
      current_pos = 1
!
      do block = 1, n_blocks
!
         if (block_significant(block, 1)) then
!
            first = block_firsts(block, 1)
            last  = block_firsts(block + 1, 1) - 1
            size  = last - first + 1
!
            array_reduced(current_pos : (current_pos + size - 1), :) = array(first : last, :)
!
            current_pos = current_pos + size
!
         endif
!
      enddo
!
   end subroutine reduce_array_int
!
!
   subroutine full_cholesky_decomposition(matrix, cholesky_vectors, dim, n_vectors,&
                                        threshold, used_diag)
!!
!!    Cholesky decomposition,
!!    Written by Sarai Dery Folkestad, June 2017
!!
!!
      implicit none
!
      integer(i15), intent(in) :: dim
      integer(i15), intent(out) :: n_vectors
!
      real(dp), intent(in) :: threshold
!
      real(dp), dimension(dim, dim), intent(inout) :: matrix
      real(dp), dimension(dim, dim), intent(out) :: cholesky_vectors
!
      integer(i15), dimension(dim, 1), optional, intent(out) :: used_diag
!
      integer(i15) :: i, j, k, index_max
      real(dp) :: max_diagonal
!
      real(dp), parameter :: tolerance = 1.0d-10
!
      if (present(used_diag)) used_diag = 0
!
!     Looping over the number of cholesky vectors
!
      do i = 1, dim
         n_vectors = i
!
!        Find the maximum diagonal
!
         index_max = 0
         max_diagonal = 0.0d0
!
         do j = 1, dim
!
            if (abs(matrix(j, j)) .gt. abs(max_diagonal)) then
!
               max_diagonal = matrix(j,j)
               index_max    = j
!
            endif
!
         enddo
!
!        Check against threshold and whether diagonal is negative
!
         if (max_diagonal .lt. 0.0d0) then
            if (abs(max_diagonal) .gt. tolerance) then
!
               write(output%unit,*)'Error: Found negative diagonal in cholesky decomposition.'
               stop
!
            endif
         endif
!
         if (abs(max_diagonal) .lt. threshold) then
!
            n_vectors = n_vectors - 1
            return
         else
            if (present(used_diag)) used_diag(n_vectors, 1) = index_max
         endif
!
!        Cholesky vectors
!
         do j = 1, dim
!
            cholesky_vectors(j,i) = matrix(j, index_max)/sqrt(max_diagonal)
!
         enddo
!
!        Subtract from matrix
!
         do j = 1, dim
            do k = 1, dim
!
               matrix(k,j) = matrix(k,j) - cholesky_vectors(k,i)*cholesky_vectors(j,i)
!
            enddo
         enddo
!
         do j = 1, dim
            matrix(j,index_max) = 0.0D0
            matrix(index_max,j) = 0.0D0
         enddo
!
      enddo

   end subroutine full_cholesky_decomposition
!
!
   subroutine full_cholesky_decomposition_effective(matrix, cholesky_vectors, dim, n_vectors,&
                                        threshold, used_diag)
!!
!!    Cholesky decomposition,
!!    Written by Sarai Dery Folkestad, June 2017.
!!
!!
      implicit none
!
      integer(i15), intent(in) :: dim
      integer(i15), intent(out) :: n_vectors
!
      real(dp), intent(in) :: threshold
!
      real(dp), dimension(dim, dim), intent(inout) :: matrix
      real(dp), dimension(dim, dim), intent(out) :: cholesky_vectors
!
      integer(i15), dimension(dim, 1), intent(out) :: used_diag
!
      integer(i15) :: i, j, k, index_max
      real(dp) :: max_diagonal, min_diagonal
!
      real(dp), dimension(:,:), allocatable :: diagonal, temp_cholesky_vector
!
      real(dp), parameter :: tolerance = 1.0d-10
!
      used_diag = 0
!
!     Looping over the number of cholesky vectors
!
      call mem%alloc(diagonal, dim, 1)
!
      do i = 1, dim
!
         diagonal(i, 1) = matrix(i, i)
!
      enddo
!
      do i = 1, dim
!
         n_vectors = i
!
!        Find the maximum diagonal
!
         index_max = 0
         max_diagonal = 0.0d0
!
         do j = 1, dim
!
            if (abs(diagonal(j, 1)) .gt. abs(max_diagonal)) then
!
               max_diagonal = diagonal(j, 1)
               index_max    = j
!
            endif
!
         enddo
!
!        Check against threshold and whether diagonal is negative
!
         if (max_diagonal .lt. 0.0d0) then
            if (abs(max_diagonal) .gt. tolerance) then
!
               write(output%unit,*)'Error: Found negative diagonal in cholesky decomposition.'
               stop
!
            endif
         endif
!
         if (abs(max_diagonal) .lt. threshold) then
!
            n_vectors = n_vectors - 1
!
            min_diagonal = 1.0D10
!
            do j = 1, dim
!
               if (diagonal(j, 1) .lt. min_diagonal) min_diagonal = diagonal(j, 1)

!
            enddo
!
            write(output%unit, '(t3, a46, e12.4)') 'The smallest diagonal after decomposition is: ', min_diagonal
            call mem%dealloc(diagonal, dim, 1)
!
            return
!
         else
!
            used_diag(n_vectors, 1) = index_max
!
         endif
!
!        Cholesky vectors
!
         cholesky_vectors(:, n_vectors) = matrix(:, index_max)
!
         if (n_vectors .gt. 1) then
!
            call mem%alloc(temp_cholesky_vector, 1, n_vectors - 1)
            temp_cholesky_vector(1, :) = cholesky_vectors(index_max, 1 : n_vectors - 1)
!
            call dgemm('N', 'T',                         &
                        dim,                             &
                        1,                               &
                        n_vectors - 1,                   &
                        -one,                            &
                        cholesky_vectors,                &
                        dim,                             &
                        temp_cholesky_vector,            &
                        1,                               &
                        one,                             &
                        cholesky_vectors(1, n_vectors),  &
                        dim)
!
            call mem%dealloc(temp_cholesky_vector, 1, n_vectors - 1)
!
         endif
!
         do j = 1, n_vectors - 1
!
            cholesky_vectors(used_diag(j,1), n_vectors) = zero
!
         enddo
!
         call dscal(dim, one/sqrt(max_diagonal), cholesky_vectors(1, n_vectors), 1)
!
         do j = 1, dim
!
            diagonal(j, 1) = diagonal(j, 1) - cholesky_vectors(j, n_vectors)**2
!
         enddo
!
         diagonal(index_max, 1) = zero
!
         do j = 1, dim
!
            matrix(j,index_max) = 0.0D0
            matrix(index_max,j) = 0.0D0
!
         enddo
!
      enddo
!
      min_diagonal = 1.0D10
!
      do j = 1, dim
!
            if (diagonal(j, 1) .lt. min_diagonal) min_diagonal = diagonal(j, 1)

!
      enddo
!
      write(output%unit, '(t3, a46, e12.4)') 'The smallest diagonal after decomposition is: ', min_diagonal
!
      call mem%dealloc(diagonal, dim, 1)
!
   end subroutine full_cholesky_decomposition_effective
!
!
   subroutine cholesky_decomposition_limited_diagonal(matrix, cholesky_vectors, dim, &
                                                     n_vectors, threshold, n_included_diagonals, &
                                                     included_diagonals, n_vectors_requested)
!!
!!    Cholesky decomposition reduced diagonal,
!!    Written by Sarai Dery Folkestad, June 2017.
!!
!!    On exit matrix_xy = matrix_xy - sum_J L_xJ*LyJ
!!
      implicit none
!
      integer(i15), intent(in) :: dim, n_included_diagonals
      integer(i15), intent(out) :: n_vectors
!
      real(dp), intent(in) :: threshold
!
      integer(i15), intent(in) :: n_vectors_requested
!
      real(dp), dimension(dim, dim), intent(inout) :: matrix
      real(dp), dimension(dim, n_included_diagonals), intent(out) :: cholesky_vectors
!
      integer(i15), dimension(n_included_diagonals, 1), intent(in) :: included_diagonals
!
      integer(i15), dimension(:, :), allocatable :: used_diag
!
      integer(i15) :: i, j, k, index_max
      real(dp) :: max_diagonal, min_diagonal
!
      real(dp), dimension(:,:), allocatable :: diagonal, temp_cholesky_vector
!
      real(dp), parameter :: tolerance = 1.0d-10
!
      call mem%alloc(diagonal, dim, 1)
!
      do i = 1, dim
!
         diagonal(i, 1) = matrix(i, i)
!
      enddo
!
      call mem%alloc_int(used_diag, dim, 1)
      used_diag = 0
!
      do i = 1, n_vectors_requested
!
         n_vectors = i
!
!        Find the maximum diagonal
!
         index_max = 0
         max_diagonal = 0.0d0
!
         do j = 1, n_included_diagonals
!
            if (abs(diagonal(included_diagonals(j, 1), 1)) .gt. abs(max_diagonal)) then
!
               max_diagonal = diagonal(included_diagonals(j, 1), 1)
               index_max    = included_diagonals(j, 1)
!
            endif
!
         enddo
!
!        Check against threshold and whether diagonal is negative
!
         if (max_diagonal .lt. 0.0d0) then
            if (abs(max_diagonal) .gt. tolerance) then
!
               write(output%unit,*)'Error: Found negative diagonal in cholesky decomposition.'
               stop
!
            endif
         endif
!
         if (abs(max_diagonal) .lt. threshold) then
!
            n_vectors = n_vectors - 1
!
            call mem%dealloc_int(used_diag, dim, 1)
!
!           On exit, cholesky vectors subtracted from matrix
!
            call dgemm('N', 'T',          &
                        dim,              &
                        dim,              &
                        n_vectors,        &
                        one,              &
                        cholesky_vectors, &
                        dim,              &
                        cholesky_vectors, &
                        dim,              &
                        -one,             &
                        matrix,           &
                        dim)
!
            return
!
         else
!
            used_diag(n_vectors, 1) = index_max
!
         endif
!
!        Cholesky vectors
!
         cholesky_vectors(:, n_vectors) = matrix(:, index_max)
!
         if (n_vectors .gt. 1) then
!
            call mem%alloc(temp_cholesky_vector, 1, n_vectors - 1)
            temp_cholesky_vector(1, :) = cholesky_vectors(index_max, 1 : n_vectors - 1)
!
            call dgemm('N', 'T',                         &
                        dim,                             &
                        1,                               &
                        n_vectors - 1,                   &
                        -one,                            &
                        cholesky_vectors,                &
                        dim,                             &
                        temp_cholesky_vector,            &
                        1,                               &
                        one,                             &
                        cholesky_vectors(1, n_vectors),  &
                        dim)
!
            call mem%dealloc(temp_cholesky_vector, 1, n_vectors - 1)
!
         endif
!
         do j = 1, n_vectors - 1
!
            cholesky_vectors(used_diag(j,1), n_vectors) = zero
!
         enddo
!
         call dscal(dim, one/sqrt(max_diagonal), cholesky_vectors(1, n_vectors), 1)
!
         do j = 1, dim
!
            diagonal(j, 1) = diagonal(j, 1) - cholesky_vectors(j, n_vectors)**2
!
         enddo
!
         diagonal(index_max, 1) = zero
!
      enddo
!
      call mem%dealloc_int(used_diag, dim, 1)
!
!     On exit, cholesky vectors subtracted from matrix
!
      call dgemm('N', 'T',          &
                  dim,              &
                  dim,              &
                  dim,              &
                  -one,             &
                  cholesky_vectors, &
                  dim,              &
                  cholesky_vectors, &
                  dim,              &
                  one,              &
                  matrix,           &
                  dim)
!
   end subroutine cholesky_decomposition_limited_diagonal
!
!
 subroutine full_cholesky_decomposition_system(matrix, cholesky_vectors, dim, n_vectors, &
                                                   threshold, used_diag)
!!
!!    Cholesky decomposition,
!!    Written by Sarai Dery Folkestad, June 2017.
!!
!!
      implicit none
!
      integer(i15), intent(in) :: dim
      integer(i15), intent(out) :: n_vectors
!
      real(dp), intent(in) :: threshold
!
      real(dp), dimension(dim, dim), intent(in)  :: matrix
      real(dp), dimension(dim, dim), intent(out) :: cholesky_vectors
!
      integer(kind=4), dimension(dim) :: used_diag
!
      integer(kind=4) :: n_vectors_4
!
      real(dp), dimension(:), allocatable :: work  ! work array for LAPACK
!
      integer(kind=4) :: info
      integer(i15) :: I, J
!
      cholesky_vectors = matrix
!
      allocate(work(2*dim))
!
!     DPSTRF computes the Cholesky factorization with complete pivoting
!     of a real symmetric positive semidefinite matrix.
!
      call dpstrf_e('L',      &
            dim,              &
            cholesky_vectors, &
            dim,              &
            used_diag,        &
            n_vectors_4,      &
            threshold,        &
            work,             &
            info)
!
      deallocate(work)
!
      n_vectors = int(n_vectors_4, kind=i15)
!
      do I = 1, dim ! Zero upper unreferenced triangle
         do J = 1, I - 1
!
            cholesky_vectors(J, I) = zero
!
         enddo
      enddo
!
      if (info .lt. 0) then
         write(*,*) info
         stop 'Cholesky decomposition failed! Something wrong in call to dpstrf'
      end if
!
   end subroutine full_cholesky_decomposition_system
!
!
   subroutine inv(Ainv, A, n)
!!
!!    Invert matrix A
!!
      implicit none
!
      integer(i15), intent(in) :: n
!
      real(dp), dimension(n,n), intent(in) :: A
      real(dp), dimension(n,n), intent(out) :: Ainv

      real(dp), dimension(n) :: work  ! work array for LAPACK
      integer(i15), dimension(n) :: ipiv   ! pivot indices
      integer(kind=4) :: info
!
!     Store A in Ainv to prevent it from being overwritten by LAPACK
!
      Ainv = A
!
!     DGETRF computes an LU factorization of a general M-by-N matrix A
!     using partial pivoting with row interchanges.
!
      call DGETRF(n, n, Ainv, n, ipiv, info)
!
      if (info /= 0) then
         stop 'Matrix is numerically singular!'
      end if
!
!     DGETRI computes the inverse of a matrix using the LU factorization
!     computed by DGETRF.
!
      call DGETRI(n, Ainv, n, ipiv, work, n, info)
!
      if (info /= 0) then
         stop 'Matrix inversion failed!'
      end if
!
   end subroutine inv
!
!
   subroutine inv_lower_tri(Ainv, A, n)
!!
!!    Invert lower triagonal matrix A
!!
      implicit none
!
      integer(i15), intent(in) :: n
!
      real(dp), dimension(n, n), intent(in) :: A
      real(dp), dimension(n, n), intent(out) :: Ainv
!
      integer(kind=4) :: info
!
!     Store A in Ainv to prevent it from being overwritten by LAPACK
!
      Ainv = A
!
!     DTRTRI computes the inverse of a real upper or lower triangular
!     matrix A.
!
      call DTRTRI('l','n', n, Ainv, n, info)
!
      if (info /= 0) then
         write(output%unit, *) 'Error: matrix inversion failed!', info
         stop
      end if
!
   end subroutine inv_lower_tri
!
!
   subroutine symmetrize(M, n)
!!
!!    Symmetrize matrix
!!    Written by Eirik F. Kjønstad, 2018
!!
!!     M <- 1/2 * (M + M^T)
!!
      implicit none
!
      integer(i15), intent(in) :: n
!
      integer(i15) :: i, j
!
      real(dp), dimension(n, n) :: M
      real(dp), dimension(:, :), allocatable :: MT
!
      call mem%alloc(MT, n, n)
!
      MT = transpose(M)
!
      do j = 1, n
         do i = 1, n
!
            M(i, j) = half*(M(i, j) + MT(i, j))
!
         enddo
      enddo
!
      call mem%dealloc(MT, n, n)
!
   end subroutine symmetrize
!
!
   subroutine anti_symmetrize(M, n)
!!
!!    Antisymmetrize matrix
!!    Written by Eirik F. Kjønstad, 2018
!!
!!     M <- 1/2 * (M - M^T)
!!
      implicit none
!
      integer(i15), intent(in) :: n
!
      integer(i15) :: i, j
!
      real(dp), dimension(n, n) :: M
      real(dp), dimension(:, :), allocatable :: MT
!
      call mem%alloc(MT, n, n)
!
      MT = transpose(M)
!
      do j = 1, n
         do i = 1, n
!
            M(i, j) = half*(M(i, j) - MT(i, j))
!
         enddo
      enddo
!
      call mem%dealloc(MT, n, n)
!
   end subroutine anti_symmetrize
!
!
   real(dp) function get_abs_max(X, n)
!!
!!    Get absolute maximum value of vector
!!    Written by Eirik F. Kjønstad, 2018
!!
!!    Returns the largest element of |X|, i.e. the largest 
!!    element in terms of absolute value
!!
      implicit none 
!
      integer(i15), intent(in) :: n
!
      real(dp), dimension(n, 1), intent(in) :: X
!
      integer(i15) :: i
!
      get_abs_max = zero
!
      do i = 1, n
!
         if (abs(X(i, 1)) .gt. get_abs_max) get_abs_max = abs(X(i, 1))
!
      enddo
!
   end function get_abs_max

!
!
   subroutine sandwich(X, A, B, n, left)
!!
!!    Sandwich
!!    Written by Eirik F. Kjønstad, 2018
!!
!!    Overwrites the n-by-n matrix X with a sandwich-ed X:
!!
!!       X <- A^T X B
!!
      implicit none
!
      integer(i15), intent(in) :: n
!
      real(dp), dimension(n, n), intent(in) :: A
      real(dp), dimension(n, n), intent(in) :: B
!
      real(dp), dimension(n, n) :: X
!
      logical, optional, intent(in) :: left 
!
      real(dp), dimension(:, :), allocatable :: tmp
!
      call mem%alloc(tmp, n, n)
!
      if (present(left)) then 
!
         if (left) then ! Transpose the left factor, X <- A^T X B
!
            call dgemm('N', 'N', &
                        n,       &
                        n,       &
                        n,       &
                        one,     &
                        X,       &
                        n,       &
                        B,       &
                        n,       &
                        zero,    &
                        tmp,     & ! tmp = X B
                        n)
!
            call dgemm('T', 'N', &
                        n,       &
                        n,       &
                        n,       &
                        one,     &
                        A,       &
                        n,       &
                        tmp,     &
                        n,       &
                        zero,    &
                        X,       & ! X = A^T tmp = A^T X B
                        n)
!
         else ! Transpose the right factor, X <- A X B^T
!
            call dgemm('N', 'T', &
                        n,       &
                        n,       &
                        n,       &
                        one,     &
                        X,       &
                        n,       &
                        B,       &
                        n,       &
                        zero,    &
                        tmp,     & ! tmp = X B^T
                        n)
!
            call dgemm('N', 'N', &
                        n,       &
                        n,       &
                        n,       &
                        one,     &
                        A,       &
                        n,       &
                        tmp,     &
                        n,       &
                        zero,    &
                        X,       & ! X = A tmp = A X B^T
                        n)
!
            endif
!
         else ! Transpose left factor (standard)
!  
            call dgemm('N', 'N', &
                        n,       &
                        n,       &
                        n,       &
                        one,     &
                        X,       &
                        n,       &
                        B,       &
                        n,       &
                        zero,    &
                        tmp,     & ! tmp = X B
                        n)
!
            call dgemm('T', 'N', &
                        n,       &
                        n,       &
                        n,       &
                        one,     &
                        A,       &
                        n,       &
                        tmp,     &
                        n,       &
                        zero,    &
                        X,       & ! X = A^T tmp = A^T X B
                        n)
!
         endif
!

!
      call mem%dealloc(tmp, n, n)
!
   end subroutine sandwich
!
!
   subroutine symmetric_sandwich(Xr, X, A, m, n)
!!
!!    Symmetric sandwich
!!    Written by Eirik F. Kjønstad, 2018
!!
!!    Overwrites the n-by-n matrix X with a sandwich-ed X:
!!
!!       X <- A^T X A
!!
      implicit none
!
      integer(i15), intent(in) :: m, n
!
      real(dp), dimension(m, n), intent(in) :: A
!
      real(dp), dimension(m, m) :: X
      real(dp), dimension(n, n) :: Xr
!
      real(dp), dimension(:, :), allocatable :: tmp
!
      call mem%alloc(tmp, m, n)
!
      call dgemm('N', 'N', &
                  m,       &
                  n,       &
                  m,       &
                  one,     &
                  X,       &
                  m,       &
                  A,       &
                  m,       &
                  zero,    &
                  tmp,     & ! tmp = X A
                  m)
!
      call dgemm('T', 'N', &
                  n,       &
                  n,       &
                  m,       &
                  one,     &
                  A,       &
                  m,       &
                  tmp,     &
                  m,       &
                  zero,    &
                  Xr,      & ! X = A^T tmp = A^T X B
                  n)
!
      call mem%dealloc(tmp, m, n)
!
   end subroutine symmetric_sandwich
!
!
   subroutine symmetric_sandwich_right(Xr, X, A, m, n)
!!
!!    Symmetric sandwich using right transposition
!!    Written by Eirik F. Kjønstad, 2018
!!
!!    Overwrites the n-by-n matrix X with a sandwich-ed X:
!!
!!       X <- A X A^T
!!
      implicit none
!
      integer(i15), intent(in) :: m, n
!
      real(dp), dimension(m, n), intent(in) :: A
!
      real(dp), dimension(n, n) :: X
      real(dp), dimension(m, m) :: Xr
!
      real(dp), dimension(:, :), allocatable :: tmp
!
      call mem%alloc(tmp, n, m)
!
      call dgemm('N', 'T', &
                  n,       &
                  m,       &
                  n,       &
                  one,     &
                  X,       &
                  n,       &
                  A,       &
                  m,       &
                  zero,    &
                  tmp,     & ! tmp = X A^T
                  n)
!
      call dgemm('N', 'N', &
                  m,       &
                  m,       &
                  n,       &
                  one,     &
                  A,       &
                  m,       &
                  tmp,     &
                  n,       &
                  zero,    &
                  Xr,      & ! X = A tmp = A X A^T 
                  m)
!
      call mem%dealloc(tmp, n, m)
!
   end subroutine symmetric_sandwich_right
!
!
end module array_utilities
