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
      integer(i15) :: n_active_spaces
!
   contains
!
!     Initialization and driver routines
!
      procedure :: init => init_mlcc2
!
      procedure :: orbital_partitioning         => orbital_partitioning_mlcc2
      procedure :: cholesky_decomposition       => cholesky_decomposition_mlcc2
      procedure ::cholesky_orbital_localization => cholesky_orbital_localization_mlcc2
!
   end type mlcc2
!
   interface
!
!
      module subroutine orbital_partitioning_mlcc2(wf)
!!
!!
         implicit none
!
         class(mlcc2) :: wf
!
      end subroutine orbital_partitioning_mlcc2
!
!
      module subroutine cholesky_orbital_localization_mlcc2(wf)
!!
!!
         implicit none
!
         class(mlcc2) :: wf
!
      end subroutine cholesky_orbital_localization_mlcc2
!
!
   module subroutine cholesky_decomposition_mlcc2(wf, density, cholesky_vectors, n_active_aos, active_ao_index_list, n_vectors)
!!
!!
      implicit none
!
      class(mlcc2)                                   :: wf
      integer(i15)                                   :: n_active_aos
      integer(i15)                                   :: n_vectors
      integer(i15), dimension( n_active_aos,1)       :: active_ao_index_list
      real(dp), dimension(wf%n_ao,wf%n_ao)           :: density
      real(dp), dimension(wf%n_ao, wf%n_ao)          :: cholesky_vectors
   !
      end subroutine cholesky_decomposition_mlcc2
!
!
      module function get_number_of_active_spaces(unit_cholesky_decomp)
!!
!!
      implicit none
!
      integer(i15) :: get_number_of_active_spaces
      integer(i15) :: unit_cholesky_decomp
!
   end function get_number_of_active_spaces
!
   module function get_number_of_active_atoms(unit_cholesky_decomp, active_space, ml_level)
!!
!!
      implicit none
!
      integer(i15)      :: get_number_of_active_atoms
      integer(i15)      :: unit_cholesky_decomp
      integer(i15)      :: active_space
      character(len=5) :: ml_level
!
   end function get_number_of_active_atoms
!
   module subroutine get_active_atoms(unit_cholesky_decomp, active_atoms, n_active_atoms, active_space, ml_level)
!!
!!
      implicit none
!
      integer(i15)      :: unit_cholesky_decomp
      integer(i15)      :: active_space, n_active_atoms
      character(len=5)  :: ml_level
!
      integer(i15), dimension(n_active_atoms,1) :: active_atoms
!
      end subroutine get_active_atoms
!
!
      module subroutine construct_active_ao_index_list(active_ao_index_list, n_active_aos, active_atoms, &
                                                   n_active_atoms, ao_center_info, n_ao)
!!
!!
      implicit none
!
      integer(i15)                               :: n_active_aos
      integer(i15)                               :: n_active_atoms
      integer(i15)                               :: n_ao
      integer(i15), dimension(n_active_aos, 1)   :: active_ao_index_list
      integer(i15), dimension(n_active_atoms, 1) :: active_atoms
      integer(i15), dimension(n_ao, 2)           :: ao_center_info
!
      end subroutine construct_active_ao_index_list
!
!
   end interface
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
      if (wf%mlcc_settings%CCS) then
         call wf%orbital_partitioning
      else
         write(unit_output,*)'Full CC2 requested, orbital partitioning skipped'
      endif
!
   end subroutine init_mlcc2
!
!
end module mlcc2_class