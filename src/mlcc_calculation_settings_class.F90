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
!  Localization method used
!
   logical :: cholesky = .false.
   logical :: cnto     = .false.
!
   real(dp) :: delta_o = 1.0D-06
   real(dp) :: delta_v = 1.0D-06
!
!  Cholesky specific variables
!
   contains
!
   end type mlcc_calculation_settings                                                                             
!
contains
!
!
end module mlcc_calculation_settings_class