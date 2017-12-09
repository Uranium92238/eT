module response_calc_specs_class
!
!!
!!                      Response calculation specifictions class module                                 
!!                Written by Sarai D. Folkestad and Eirik F. Kjønstad, Nov 2017         
!!
!!    The class contains specifications for the response calculation, and necesary variables 
!!    provided by user.
!! 
!
!  :::::::::::::::::::::::::::::::::::
!  -::- Modules used by the class -::-
!  :::::::::::::::::::::::::::::::::::
!
   use types
!
!  :::::::::::::::::::::::::::::::::::::::::::::::::::::
!  -::- Definition of the response_calc_specs class -::-
!  :::::::::::::::::::::::::::::::::::::::::::::::::::::
!
   type :: response_calc_specs
!
!     Settings for the excited state calculation
!
      real(dp) :: energy_threshold   = 1.0D-6 ! Threshold for energies 
      real(dp) :: residual_threshold = 1.0D-6 ! Threshold for equation residuals
!
      integer(i15) :: max_iterations = 50 ! excited state maximum number of iterations 
!
      logical :: restart = .false.
!
   end type response_calc_specs                                                                          
!
contains
!
!   Nothing yet
!
end module response_calc_specs_class
