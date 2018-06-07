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
      character(len=100) :: basis
!
      real(dp) :: x
      real(dp) :: y
      real(dp) :: z
!
   contains
!
   end type atom
!
end module atom_class
