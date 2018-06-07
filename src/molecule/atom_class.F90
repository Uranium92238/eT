module atom_class
!
!!
!!    Atom class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   use kinds
!
   implicit none
!
   type :: atom
!
      character(len=2) :: type
!
      character(len=40) :: basis
!
      real(dp), dimension(3, 1) :: xyz
!
   contains
!
   end type atom
!
end module atom_class
