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
   logical function is_significant(vec, dim, threshold)
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
!
      real(dp), intent(in)  :: threshold
!
      integer(i15) :: i = 0
!
      is_significant = .false.
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
            array_reduced(current_pos : current_pos + size - 1, 1 : columns) = array(first : last, 1 : columns)
            current_pos = current_pos + size
!
         endif
!
      enddo
!
   end subroutine reduce_array
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
   module subroutine full_cholesky_decomposition(matrix, cholesky_vectors, dim, n_vectors,&
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
            if (matrix(j,j) .gt. max_diagonal) then
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
!        Subtract from density
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
   function inv(A, dim) result(Ainv)
!!
!!    Invert matrix A
!!
      implicit none
!
      integer(i15), intent(in) :: dim
!
      real(dp), dimension(dim,dim), intent(in) :: A
      real(dp), dimension(dim,dim) :: Ainv

      real(dp), dimension(dim) :: work  ! work array for LAPACK
      integer(kind=4), dimension(dim) :: ipiv   ! pivot indices
      integer(kind=4) :: n, info
!
!     Store A in Ainv to prevent it from being overwritten by LAPACK
!
      Ainv = A
      n = int(dim, kind=4)
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
   end function inv
!
!
end module array_utilities
