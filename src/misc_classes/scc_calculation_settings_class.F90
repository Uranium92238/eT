module scc_calculation_settings_class
!
!! 
!! 						SCC calculation settings class 
!! 					Written by Eirik F. Kjønstad, Dec 2017
!!
!
!  :::::::::::::::::::::::::::::::::::
!  -::- Modules used by the class -::-
!  :::::::::::::::::::::::::::::::::::
!
   use types
!
!  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  -::- Definition of the scc_calculation_settings class -::-
!  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
   type :: scc_calculation_settings 
!
!     Whether to do a restart of an SCCSD calculations 
!
      logical :: restart = .false.
!
!     Threshold for convergence of zero overlap 
!
      real(dp) :: overlap_threshold = 1.0D-6
!
   contains
!
!     Nothing yet
!
   end type scc_calculation_settings                                                                             
!
contains
!
!  Nothing yet
!
end module scc_calculation_settings_class
