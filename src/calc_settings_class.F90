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
      real(dp) :: energy_threshold   = 1.0D-07 ! Threshold for energies 
      real(dp) :: equation_threshold = 1.0D-07 ! Threshold for equation residuals
!
      integer(i15) :: ground_state_max_iterations = 50! Ground state maximum number of iterations 
!
      logical :: restart = .false.
!
   end type calc_settings                                                                            
!
contains
!
!  No procedures yet.
!
end module calc_settings_class