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
      integer(i15)                              :: n_CCS_o = 0
      integer(i15)                              :: n_CCS_v = 0
!          
      integer(i15), dimension(:,:), allocatable :: n_CC2_o
      integer(i15), dimension(:,:), allocatable :: n_CC2_v
!
   contains
!
!     Initialization and driver routines
!
      procedure :: init                    => init_mlcc2
      procedure :: initialize_orbital_info => initialize_orbital_info_mlcc2
      procedure :: destruct_orbital_info   => destruct_orbital_info_mlcc2
!
      procedure :: orbital_partitioning          => orbital_partitioning_mlcc2
      procedure :: cholesky_decomposition        => cholesky_decomposition_mlcc2
      procedure :: cholesky_localization         => cholesky_localization_mlcc2
      procedure :: cholesky_orbitals             => cholesky_orbitals_mlcc2
      procedure :: cholesky_orbital_drv          => cholesky_orbital_drv_mlcc2
!
!     Omega
!
      procedure :: omega_mlcc2_a1_mlcc2 => omega_mlcc2_a1_mlcc2
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
   module subroutine cholesky_localization_mlcc2(wf, orbitals, orbital_energies,&
                                              n_nuclei, ao_center_info, n_ao_on_center, unit_cholesky_decomp)
!!
!!
      implicit none
!
      class(mlcc2) :: wf
      real(dp), dimension(wf%n_ao, wf%n_mo) :: orbitals
      real(dp), dimension(wf%n_mo, 1)       :: orbital_energies
      integer(i15)                          :: n_nuclei
      integer(i15), dimension(wf%n_ao,2)    :: ao_center_info
      integer(i15), dimension(n_nuclei, 1)  :: n_ao_on_center
      integer(i15)                          :: unit_cholesky_decomp
!
      end subroutine cholesky_localization_mlcc2
!
!
   module subroutine cholesky_decomposition_mlcc2(wf, density, cholesky_vectors,&
                                                     n_vectors, selection, n_active_aos, active_ao_index_list)
!!
!!
      implicit none
!
      class(mlcc2)                                       :: wf
      integer(i15)                                       :: n_active_aos
      integer(i15)                                       :: n_vectors
      real(dp), dimension(wf%n_ao,wf%n_ao)               :: density
      real(dp), dimension(wf%n_ao, wf%n_ao)              :: cholesky_vectors
      logical                                            :: selection  
      integer(i15), dimension( n_active_aos,1), optional :: active_ao_index_list
   !
      end subroutine cholesky_decomposition_mlcc2
!
!
      module subroutine cholesky_orbitals_mlcc2(wf, cholesky_vectors, n_vectors, orbitals, orbital_energies, ao_fock)
!!
!!
!!
      implicit none
!
      class(mlcc2)                              :: wf
      real(dp), dimension(wf%n_ao, wf%n_ao)     :: cholesky_vectors
      real(dp), dimension(wf%n_ao, n_vectors)   :: orbitals
      real(dp), dimension(n_vectors, 1)         :: orbital_energies
      integer(i15)                              :: n_vectors
      real(dp), dimension(wf%n_ao,wf%n_ao)      :: ao_fock
!
      end subroutine cholesky_orbitals_mlcc2
!
   module subroutine cholesky_orbital_drv_mlcc2(wf, orbitals, orbital_energies, offset, ao_fock, density, n_vectors,&
                              selection, n_active_aos, active_ao_index_list)
!!
!!
      implicit none
!
      class(mlcc2)                                       :: wf
      real(dp), dimension(wf%n_ao, wf%n_mo)              :: orbitals
      real(dp), dimension(wf%n_mo, 1)                    :: orbital_energies
      real(dp), dimension(wf%n_ao, wf%n_ao)              :: ao_fock
      integer(i15)                                       :: n_active_aos, offset
      integer(i15)                                       :: n_vectors
      real(dp), dimension(wf%n_ao,wf%n_ao)               :: density
      logical                                            :: selection
      integer(i15), dimension( n_active_aos,1), optional :: active_ao_index_list
!
   end subroutine cholesky_orbital_drv_mlcc2
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
      module subroutine read_atom_info(n_nuclei, n_ao)
!!
!!
!!
         implicit none
!
         integer(i15) :: n_nuclei,n_ao 
!
   end subroutine read_atom_info
!
   module subroutine read_center_info(n_nuclei, n_ao, n_ao_on_center, ao_center_info)
!!
!!
!!
      implicit none
!
      integer(i15) :: n_nuclei
      integer(i15) :: n_ao
      integer, dimension(n_nuclei, 1)  :: n_ao_on_center
      integer, dimension(n_ao, 2)      :: ao_center_info
!
   end subroutine read_center_info
!
!
    module subroutine omega_mlcc2_a1_mlcc2(wf, active_space)
! 
!     Omega A1
!     Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!   
!     Calculates the A1 term of omega, 
!   
!     A1: sum_ckd g_adkc * u_ki^cd,
!  
!     and adds it to the projection vector (omega1) of
!     the wavefunction object wf
! 
!     u_ki^cd = 2*s_ki^cd - s_ik^cd 
! 
      implicit none
!
      class(mlcc2)   :: wf
      integer(i15) :: active_space
!
   end subroutine omega_mlcc2_a1_mlcc2
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
!
         write(unit_output,*)'Full CC2 requested, orbital partitioning skipped'
         wf%n_active_spaces = 1
         call wf%initialize_orbital_info
         wf%n_CC2_o = wf%n_o
         wf%n_CC2_v = wf%n_v
!
      endif
!
!     Initialize amplitudes and associated attributes
!
      call wf%initialize_amplitudes
!
!     Set the number of parameters in the wavefunction
!     (that are solved for in the ground and excited state solvers) 
!
      wf%n_parameters = wf%n_t1am
!
!     Read Cholesky AO integrals and transform to MO basis
!
      call wf%read_transform_cholesky
!
!     Initialize fock matrix
!
      call wf%initialize_fock_matrix
!
   end subroutine init_mlcc2
!
   subroutine initialize_orbital_info_mlcc2(wf)
!!
!!
      implicit none
!
      class(mlcc2) :: wf
!
      if (.not. allocated(wf%n_CC2_o)) call allocator_int(wf%n_CC2_o, wf%n_active_spaces, 1)
      if (.not. allocated(wf%n_CC2_v)) call allocator_int(wf%n_CC2_v, wf%n_active_spaces, 1)
      wf%n_CC2_o = 0
      wf%n_CC2_v = 0
!
   end subroutine initialize_orbital_info_mlcc2
!
!
   subroutine destruct_orbital_info_mlcc2(wf)
!!
!!
      implicit none
!
      class(mlcc2) :: wf
!
      if (allocated(wf%n_CC2_o)) call deallocator_int(wf%n_CC2_o, wf%n_active_spaces, 1)
      if (allocated(wf%n_CC2_v)) call deallocator_int(wf%n_CC2_v, wf%n_active_spaces, 1)
!
   end subroutine destruct_orbital_info_mlcc2
!
end module mlcc2_class