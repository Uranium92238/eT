module mlcc_calculation_settings_class
!
!! 
!
!  :::::::::::::::::::::::::::::::::::
!  -::- Modules used by the class -::-
!  :::::::::::::::::::::::::::::::::::
!
   use types
!
!  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  -::- Definition of the mlcc_calculation_settings class -::-
!  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
   type :: mlcc_calculation_settings 
!
!  Levels of the hierarchy
!
   logical :: CCS  = .false.
   logical :: CC2  = .false.
   logical :: CCSD = .false.
   logical :: CC3  = .false.
!
   contains
!
   end type mlcc_calculation_settings                                                                             
!
contains
!
!
end module mlcc_calculation_settings_class
