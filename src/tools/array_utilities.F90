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
   subroutine reduce_vector(vec, reduced_vec, block_firsts, block_significant, n_blocks, dim, dim_reduced)
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
      integer(i15), dimension(n_blocks, 1) :: block_significant
!
      real(dp), dimension(dim, 1) :: vec
      real(dp), dimension(dim_reduced, 1) :: vec_reduced
!
      current_pos = 1
!
      do block = 1, n_blocks
!
         if (block_significant(block, 1)) then
!
            first = block_firsts(block)
            last  = block_firsts(block + 1) - 1
            size  = last - first + 1
!
            reduce_vec(current_pos : current_pos + size - 1, 1) = vec(first : last, 1)
            current_pos = current_pos + size
!
         endif
!
      enddo
!
   end subroutine reduce_vector
!
!
   subroutine reduce_array(array, array_vec, block_firsts, block_significant, n_blocks, dim, dim_reduced, columns)
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
      integer(i15), dimension(n_blocks, 1) :: block_significant
!
      real(dp), dimension(dim, :) :: array
      real(dp), dimension(dim_reduced, :) :: array_reduced
!
      current_pos = 1
!
      do block = 1, n_blocks
!
         if (block_significant(block, 1)) then
!
            first = block_firsts(block)
            last  = block_firsts(block + 1) - 1
            size  = last - first + 1
!
            reduce_vec(current_pos : current_pos + size - 1, 1 : columns) = vec(first : last, 1 : columns)
            current_pos = current_pos + size
!
         endif
!
      enddo
!
   end subroutine reduce_array
!
!
end module array_utilities
