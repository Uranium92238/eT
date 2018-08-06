module active_atoms_class
!
!!
!!    Active atoms class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
!
   implicit none
!
   type :: active_atoms
!
      integer(i15) :: n_acitve_atoms
!
      real(15), dimension(:,:), allocatable :: active_atoms
!
   contains
!
!
   end type active_atoms
!
!
contains
!
!
end module active_atoms_class
