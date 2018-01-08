module ground_state_specs_class
!
!!
!!                       Ground State Specifictions class module                                 
!!            Written by Sarai D. Folkestad and Eirik F. Kjønstad, Nov 2017         
!!
!!    The class contains specifications for the ground state calculation.
!!    
!
!  :::::::::::::::::::::::::::::::::::
!  -::- Modules used by the class -::-
!  :::::::::::::::::::::::::::::::::::
!
   use types
!
!  :::::::::::::::::::::::::::::::::::::::::::::::::
!  -::- Definition of the ground_state_specs class -::-
!  :::::::::::::::::::::::::::::::::::::::::::::::::
!
   type :: ground_state_specs
!
!     Settings for the ground state calculation
!
      real(dp) :: energy_threshold   = 1.0D-6 ! Threshold for energies 
      real(dp) :: residual_threshold = 1.0D-6 ! Threshold for equation residuals
!
      integer(i15) :: max_iterations = 50 ! Ground state maximum number of iterations 
!
      logical :: restart = .false.
!
   end type ground_state_specs                                                                          
!
contains
!
end module ground_state_specs_class