module interval_class
!!
!!    Interval class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
   use kinds
!
   implicit none
!
   type interval
!
      integer(i15) :: first = -1
      integer(i15) :: last  = -1
      integer(i15) :: size  = -1
!
   contains
!
   end type interval
!
end module interval_class
