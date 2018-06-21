module shell_class
!
!!
!!    Shell class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   use kinds
!
   implicit none
!
   type :: shell
!
      integer(i15) :: size  ! The number of basis functions in the shell
      integer(i15) :: first ! The first AO index of the shell
      integer(i15) :: l     ! The angular momentum of the shell, e.g. l = 1 denotes p-orbitals
!
   contains
!
!     No routines yet
!
   end type shell
!
contains
!
!  No routines yet
!
end module
