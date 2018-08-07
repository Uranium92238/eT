module active_atoms_class
!
!!
!!    Active atoms class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   use kinds
!
!
   implicit none
!
   type :: active_atoms
!
      integer(i15) :: n_acitve_atoms
!
      real(dp), dimension(:,:), allocatable :: atoms
!
   contains
!
!
   end type active_atoms
!
!
contains
!
!     Sanity checks for active atoms
!
end module active_atoms_class
