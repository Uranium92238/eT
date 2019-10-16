module active_atoms_class
!
!!
!!    Active atom space class module
!!    Written by Sarai D. Folkestad, Apr 2019
!!
!
   implicit none
!
   type :: active_atoms
!
      character(len=100) :: level
!
      integer :: first_atom, last_atom
!
   end type active_atoms
!
contains
!
!
end module active_atoms_class