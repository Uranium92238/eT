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
end module array_utilities
