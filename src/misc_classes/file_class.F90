module file_class
!
!!
!!                               File class module                                 
!!             Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2018         
!!
!! 
!
!  :::::::::::::::::::::::::::::::::::
!  -::- Modules used by the class -::-
!  :::::::::::::::::::::::::::::::::::
!
   use types
   use input_output
!
!  ::::::::::::::::::::::::::::::::::::::
!  -::- Definition of the file class -::-
!  ::::::::::::::::::::::::::::::::::::::
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
   end type file                                                                             
!
!
contains
!
!
end module file_class
