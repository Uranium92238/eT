module active_atoms_info_class
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
   type :: active_atoms_info
!
      integer(i15) :: n_active_atoms
!
      integer(i15), dimension(:,:), allocatable :: atoms
!
   contains
!
!
   end type active_atoms_info
!
!
contains
!
!     Sanity checks for active atoms
!
end module active_atoms_info_class
