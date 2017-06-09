module mlcc2_class
!
!!
!!                 Multi-level CC2 (MLCC2) class module                                
!!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017         
!!                                                                           
!
!  :::::::::::::::::::::::::::::::::::
!  -::- Modules used by the class -::-
!  :::::::::::::::::::::::::::::::::::
!
!  General tools
!
   use types
   use utils
   use workspace
   use input_output
   use input_reader
   use mlcc_calculation_settings_class
!
!  The ancestor class module (CCS)
!
   use ccs_class
!
   implicit none 
!
!  :::::::::::::::::::::::::::::::::::::::
!  -::- Definition of the MLCC2 class -::-
!  ::::::::::::::::::::::::::::::::::::::: 
!
   type, extends(ccs) :: mlcc2
!
      type(mlcc_calculation_settings) :: mlcc_settings
!
   contains
!
!     Initialization and driver routines
!
      procedure :: init => init_mlcc2
!
   end type mlcc2
!
!
contains
!
!
   subroutine init_mlcc2(wf)
!!
!!
      implicit none
!
      class(mlcc2) :: wf
!
!     Set model name 
!
      wf%name = 'MLCC2'
!
!     MLCC sanity check
!
      if(wf%mlcc_settings%CC3 .or. wf%mlcc_settings%CCSD) then
         write(unit_output,*)'WARNING: CC3 and CCSD active spaces not available for MLCC2'
         stop
      endif
!
!     Set implemented methods
!
!
!     Read Hartree-Fock info
!
     call wf%read_hf_info
!
!     Orbital partitioning
!
   end subroutine init_mlcc2
!
!
end module mlcc2_class