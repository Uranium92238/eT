module calc_settings_class 
!
!!
!!                         Calculation settings class module                                 
!!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017         
!!
!
!  :::::::::::::::::::::::::::::::::::
!  -::- Modules used by the class -::-
!  :::::::::::::::::::::::::::::::::::
!
   use types
!
!  :::::::::::::::::::::::::::::::::::::::::::::::
!  -::- Definition of the calc_settings class -::-
!  :::::::::::::::::::::::::::::::::::::::::::::::
!
   type :: calc_settings 
!
      real(dp) :: disk_space = 100D0      ! 100 gigabytes; can be modified by user
      real(dp) :: memory     = 700000000  ! In words (= ca 5 gb) !! This should be hooked up to workspace!!
!
      character(len=40) :: print_level = 'developer'
!
   end type calc_settings                                                                            
!
contains
!
!  No procedures yet.
!
end module calc_settings_class