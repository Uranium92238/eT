module file_class
!
!!
!!		File class module
!!		Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
	use kinds
!
	implicit none
!
! 	Definition of the file class
!
   type :: file
!
!     Filename
!
      character(len=40) :: name = 'no_name'
!
!     Unit identifier
!
      integer(i15) :: unit = -1
!
!     File size (in bytes)
!
      integer(i15) :: size = -1
!
!     Logical for whether the file is currently opened or not
!
      logical :: opened = .false.
!
   contains
!
! 		No routines yet
!
   end type file
!
! 	Declaration of the main output file used throughout eT
!
   type(file) :: output ! Main output file
!
end module file_class
