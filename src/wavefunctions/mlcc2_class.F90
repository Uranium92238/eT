module mlcc2_class
!
!!
!!
!!
!!         Multilevel Coupled cluster singles and perturbative doubles (MLCC2) class module                                 
!!                Written by Sarai D. Folkestad and Eirik F. Kjønstad, Jun 2017         
!! 
!!
!!    This module contains the definition of the multilevel Coupled cluster singles and perturbative doubles (MLCC2)
!!    wavefunction class. It is structured into four sections:
!!
!!       1. Modules used by the class: 
!!
!!             Basic utilities and the ancestor class
!!
!!       2. Definition of the class: 
!!
!!             Non-inherited variables, followed by non-inherited or overridden procedures
!!
!!       3. Interfaces to submodules:
!!
!!             The procedures in the class are grouped according to functionality, with
!!             detailed definitions given in the following class submodules:
!!                
!!                - Input Reader
!!                - Orbital partitioning
!!                - Ground state
!!                - Omega
!!                - Excited state 
!!                - Jacobian (right transformation)
!!                 
!!
!!             The interfaces shows incoming variables and their type, but contains 
!!             no information of the procedure itself. The procedure is shown in full 
!!             in the corresponding submodule. 
!!
!!       4. Class module routines (i.e., non-submodule procedures). These include
!!          the initialization and driver routines of the class, along with procedures that
!!          are not (yet, at least) easily gathered in a submodule.
!!                                                                    
!
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
   use mlcc_orbitals_class
   use mlcc_calculation_settings_class
!
!  The ancestor class module (CCS)
!
   use ccs_class
!
   implicit none 
!
!
!  ::::::::::::::::::::::::::::::::::::::::::
!  -::- 2. Definition of the MLCC2 class -::-
!  ::::::::::::::::::::::::::::::::::::::::::
!
!
   type, extends(ccs) :: mlcc2
!
      type(mlcc_calculation_settings)  :: mlcc_settings
!
!     ML variables
!
      integer(i15) :: n_CCS_o = 0
      integer(i15) :: n_CCS_v = 0
!
      integer(i15) :: first_CCS_o = 0
      integer(i15) :: first_CCS_v = 0
!          
      integer(i15) :: n_CC2_o = 0
      integer(i15) :: n_CC2_v = 0
!
      integer(i15) :: first_CC2_o = 0
      integer(i15) :: first_CC2_v = 0
!
      type(mlcc_orbitals) :: CC2_orbitals
!
!     Excited state variables
!
      integer(i15)                           :: n_x2am = 0 
      real(dp), dimension(:,:), allocatable  :: x2am 
!
   contains
!
!
!     -::- Initialization routine -::-
!     --------------------------------
!
      procedure :: init => init_mlcc2
!
!
!     -::- Other class routine pointers not located in submodules -::-
!     ----------------------------------------------------------------
!
      procedure :: get_CC2_active_indices       => get_CC2_active_indices_mlcc2
      procedure :: get_CC2_n_active             => get_CC2_n_active_mlcc2
!
      procedure :: calc_energy                  => calc_energy_mlcc2
!
!     Routines to allocate amplitudes
!
      procedure :: initialize_amplitudes              => initialize_amplitudes_mlcc2
      procedure :: initialize_cc2_double_amplitudes   => initialize_cc2_double_amplitudes_mlcc2
!
!     Routines to deallocate amplitudes
!
      procedure :: destruct_amplitudes             => destruct_amplitudes_mlcc2
      procedure :: destruct_cc2_double_amplitudes  => destruct_cc2_double_amplitudes_mlcc2
!
!     Routine to save the amplitudes to disk 
!
      procedure :: save_amplitudes              => save_amplitudes_mlcc2
!
!!    Routines to read the amplitudes from disk (and allocate if necessary)
!
      procedure :: read_amplitudes              => read_amplitudes_mlcc2
      procedure :: read_cc2_double_amplitudes   => read_cc2_double_amplitudes_mlcc2
!
!     -::- Input reader submodule routine pointers -::-
!     -------------------------------------------------
!
      procedure :: mlcc_reader         => mlcc_reader_mlcc2
      procedure :: read_orbital_info   => read_orbital_info_mlcc2
!
!
!     -::- Orbital partitioning submodule routine pointers -::-
!     ---------------------------------------------------------
!
      procedure :: orbital_partitioning         => orbital_partitioning_mlcc2
      procedure :: cholesky_decomposition       => cholesky_decomposition_mlcc2
      procedure :: cholesky_localization_drv    => cholesky_localization_drv_mlcc2
      procedure :: cholesky_orbitals            => cholesky_orbitals_mlcc2
      procedure :: cholesky_orbital_constructor => cholesky_orbital_constructor_mlcc2
!
      procedure :: cnto_orbital_drv             => cnto_orbital_drv_mlcc2
      procedure :: cnto_lower_level_method      => cnto_lower_level_method_mlcc2
      procedure :: cnto_orbitals                => cnto_orbitals_mlcc2
      procedure :: cnto_init_ccs                => cnto_init_ccs_mlcc2
!
      procedure :: print_orbital_info           => print_orbital_info_mlcc2
!
!
!     -::- Ground state submodule routine pointers -::-
!     -------------------------------------------------
!
      procedure :: ground_state_preparations => ground_state_preparations_mlcc2
      procedure :: ground_state_cleanup      => ground_state_cleanup_mlcc2
!
!
!     -::- Omega submodule routine pointers -::-
!     ------------------------------------------
!
      procedure :: omega_mlcc2_a1                     => omega_mlcc2_a1_mlcc2
      procedure :: omega_mlcc2_b1                     => omega_mlcc2_b1_mlcc2
      procedure :: construct_omega                    => construct_omega_mlcc2
      procedure :: get_s2am                           => get_s2am_mlcc2
!
!
!     -::- Excited state submodule routine pointers -::-
!     --------------------------------------------------
!
      procedure :: excited_state_preparations        => excited_state_preparations_mlcc2
!      procedure :: excited_state_cleanup             => excited_state_cleanup_mlcc2
!
      procedure :: calculate_orbital_differences      => calculate_orbital_differences_mlcc2
      procedure :: transform_trial_vectors            => transform_trial_vectors_mlcc2
!
      procedure :: cvs_residual_projection            => cvs_residual_projection_mlcc2 
!
      procedure :: analyze_double_excitation_vector   => analyze_double_excitation_vector_mlcc2
      procedure :: summary_excited_state_info         => summary_excited_state_info_mlcc2
      procedure :: print_excitation_vector            => print_excitation_vector_mlcc2
!
!
!     -::- Jacobian submodule routine pointers -::-
!     ---------------------------------------------
!
      procedure :: jacobian_mlcc2_transformation   => jacobian_mlcc2_transformation_mlcc2
      procedure :: jacobian_mlcc2_a1               => jacobian_mlcc2_a1_mlcc2
      procedure :: jacobian_mlcc2_b1               => jacobian_mlcc2_b1_mlcc2
      procedure :: jacobian_mlcc2_a2               => jacobian_mlcc2_a2_mlcc2
      procedure :: jacobian_mlcc2_b2               => jacobian_mlcc2_b2_mlcc2
!
      procedure :: cvs_rho_aibj_projection         => cvs_rho_aibj_projection_mlcc2
!
   end type mlcc2
!
!
!  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  -::- 3. Interfaces to the submodules of the CCS class -::- 
!  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
   interface
!
!
!     -::- Input reader submodule routine interface -::-
!     --------------------------------------------------
!
      module subroutine mlcc_reader_mlcc2(wf, unit_input)
!!
!!
!!
         implicit none
!
         integer(i15)      :: unit_input
!
         class(mlcc2)      :: wf
!
      end subroutine mlcc_reader_mlcc2
!
!
      module subroutine read_orbital_info_mlcc2(wf, unit_input)
!!
!!
         implicit none
!
         integer(i15)      :: unit_input
!
         class(mlcc2)      :: wf
!
      end subroutine read_orbital_info_mlcc2
!
!
   end interface
!
!
   interface
!
!     -::- Orbital partitioning submodule interface -::-
!     ::::::::::::::::::::::::::::::::::::::::::::::::::
!
      module subroutine orbital_partitioning_mlcc2(wf)
!!
!!       Orbital partitioning,
!!       Written by Sarai D. Folkestad, June 2017
!!
         implicit none
!
         class(mlcc2) :: wf
!
      end subroutine orbital_partitioning_mlcc2
!
!
      module subroutine cholesky_localization_drv_mlcc2(wf)
!!
!!       Cholesky orbital localization. driver,
!!       Written by Sarai D. Folkestad, June 2017
!!
         implicit none
!
         class(mlcc2) :: wf
!
      end subroutine cholesky_localization_drv_mlcc2
!
!
!
      module subroutine cholesky_decomposition_mlcc2(wf, density, cholesky_vectors,&
                                                     n_vectors, selection, n_active_aos, active_ao_index_list)
!!
!!       Cholesky decomposition, 
!!       Written by Sarai dery Folkestad, June 2017.
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
!!       Cholesky orbitals,
!!       Written by Sarai Dery Folkestad, June 2017
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
      module subroutine cholesky_orbital_constructor_mlcc2(wf, orbitals, orbital_energies, offset, ao_fock, density, n_vectors,&
                              selection, n_active_aos, active_ao_index_list)
!!
!!       Cholesky orbital constructor,
!!       Written by Sarai Dery Folkestad, June 2017
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
      end subroutine cholesky_orbital_constructor_mlcc2
!
!
      module subroutine construct_active_ao_index_list(active_ao_index_list, n_active_aos, active_atoms, &
                                                   n_active_atoms, ao_center_info, n_ao)
!!
!!       Construct active ao index list,
!!       Written by Sarai Dery Folkestad, June 2017.
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
      module subroutine cnto_orbital_drv_mlcc2(wf)
!!
!!       CNTO orbital driver,
!!       Written by Sarai D. Folkestad, June 2017.
!!
         implicit none 
!
         class(mlcc2) :: wf
!
      end subroutine cnto_orbital_drv_mlcc2
!
!
      module subroutine cnto_lower_level_method_mlcc2(wf)
!!
!!       CNTO lower level calculation (MLCC2),
!!       Written by Sarai D. Folkestad, June 2017.
!!
         implicit none 
!
         class(mlcc2) :: wf
!
      end subroutine cnto_lower_level_method_mlcc2
!
!
      module subroutine cnto_orbitals_mlcc2(wf)
!!
!!       CNTO Oritals (MLCC2),
!!       Written by Sarai D. Folkestad Aug. 2017
!!   
         implicit none
!
         class(mlcc2) :: wf
!
      end subroutine cnto_orbitals_mlcc2
!
!
      module subroutine print_orbital_info_mlcc2(wf)
!!
!!       Print CNTO info, 
!!       Written by Sarai D. Folkestad, Aug. 2017
!!
         implicit none 
!
         class(mlcc2) :: wf
!
      end subroutine print_orbital_info_mlcc2
!
!
      module subroutine cnto_init_ccs_mlcc2(wf, ccs_wf)
!!
!!       CNTO Oritals (MLCC2),
!!       Written by Sarai D. Folkestad Aug. 2017
!!   
         implicit none
!
         class(mlcc2)   :: wf
         class(ccs)      :: ccs_wf
!
      end subroutine cnto_init_ccs_mlcc2
!
!
   end interface
!
!
   interface
! 
!     -::- Omega submodule interface -::-
!     :::::::::::::::::::::::::::::::::::
!
      module subroutine omega_mlcc2_a1_mlcc2(wf)
!! 
!!       Omega A1
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!  
         implicit none
!
         class(mlcc2)   :: wf
!
      end subroutine omega_mlcc2_a1_mlcc2
!
!
      module subroutine omega_mlcc2_b1_mlcc2(wf)
!!
!!       Omega B1
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
         implicit none
!
         class(mlcc2)   :: wf
      
      end subroutine omega_mlcc2_b1_mlcc2
!
!
      module subroutine construct_omega_mlcc2(wf)
!! 
!!       Construct Omega (CC2)
!!       Written by Eirik F. Kjønstad and Sarai Folkestad, Apr 2017
!! 
         implicit none 
!
         class(mlcc2) :: wf
!
      end subroutine construct_omega_mlcc2
!
      module subroutine get_s2am_mlcc2(wf, s_ai_bj)
!!
!!       Get S_2 amplitudes, 
!!       Written by Sarai D. Folkestad, July 2017 
!!
         implicit none
!
         class(mlcc2) :: wf
! 
         real(dp), dimension((wf%n_CC2_v)*(wf%n_CC2_o), (wf%n_CC2_v)*(wf%n_CC2_o)) :: s_ai_bj
!
      end subroutine get_s2am_mlcc2
!
!
   end interface 
!
!
   interface
!
!
!     -::- Ground state submodule routine  interface -::-
!     ---------------------------------------------------
!
      module subroutine ground_state_preparations_mlcc2(wf)
!!
!!       Ground State Preparations (mlcc2)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
         class(mlcc2) :: wf 
!
      end subroutine ground_state_preparations_mlcc2
!
!
      module subroutine ground_state_cleanup_mlcc2(wf)
!!
!!       Ground State Cleanup (mlcc2)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
         class(mlcc2) :: wf 
!
      end subroutine ground_state_cleanup_mlcc2
!  
!
   end interface
!
!
   interface 
!
!     -::- Excited state submodule interface -::-
!     :::::::::::::::::::::::::::::::::::::::::::
!
!
      module subroutine excited_state_preparations_mlcc2(wf)
!!
!!       Excited State Preparations (MLCC2)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
         class(mlcc2) :: wf 
!
!        Do nothing for mlcc2
!
      end subroutine excited_state_preparations_mlcc2
!  
!
      module subroutine excited_state_cleanup_mlcc2(wf)
!!
!!       Excited State Cleanup (MLCC2)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
         implicit none
!
         class(mlcc2) :: wf 
!
      end subroutine excited_state_cleanup_mlcc2
!  
!
      module subroutine calculate_orbital_differences_mlcc2(wf, orbital_diff)
!!
!!       Calculate Orbital Differences (MLCC2)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad May 2017
!!
         implicit none
!
         class(mlcc2) :: wf
!
         real(dp), dimension(wf%n_parameters, 1) :: orbital_diff
!
      end subroutine calculate_orbital_differences_mlcc2
!
!
      module subroutine transform_trial_vectors_mlcc2(wf, first_trial, last_trial)
!!
!!       Transformation Trial Vectors (MLCC2)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
         implicit none
!
         class(mlcc2) :: wf
!
         integer(i15), intent(in) :: first_trial, last_trial ! Which trial_vectors we are to transform
!
      end subroutine transform_trial_vectors_mlcc2
!
!
      module subroutine cvs_residual_projection_mlcc2(wf, residual)
!!
!!       Residual projection (MLCC2), 
!!       Written by Sarai D. Folkestad Aug. 2017    
!!
         implicit none
!
         class(mlcc2) :: wf
         real(dp), dimension(wf%n_parameters, 1) :: residual
!
      end subroutine cvs_residual_projection_mlcc2
!
!
     module subroutine print_excitation_vector_mlcc2(wf, vec, unit_id)
!!
!!       Print excitation vector,
!!       Written by Eirik F. Kjønstad and Sarai D. Folekstad,Jun 2017
!!
         implicit none
!  
         class(mlcc2) :: wf
!
         real(dp), dimension(wf%n_parameters, 1) :: vec
!
         integer(i15) :: unit_id     
!
      end subroutine print_excitation_vector_mlcc2
!
!
      module subroutine analyze_double_excitation_vector_mlcc2(wf, vec, n, sorted_short_vec, index_list)
!!
!!       Analze double excitation vector,
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jun 2017
!!
         implicit none
!  
         class(mlcc2) :: wf
!
         real(dp), dimension(wf%n_x2am, 1) :: vec    
!
         integer(i15) :: a = 0, i = 0, ai = 0, b = 0, j = 0, bj = 0, aibj = 0, k = 0
!
         integer(i15) :: n    ! Number of elements wanted
!  
         real(dp), dimension(n, 1)    :: sorted_short_vec
!  
         integer(i15), dimension(n, 4) ::index_list
!
      end subroutine analyze_double_excitation_vector_mlcc2
!
!
      module subroutine summary_excited_state_info_mlcc2(wf, energies)
!!
!!       Summary of excited state information,
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad
!!
         implicit none
!  
         class(mlcc2) :: wf
!
         real(dp), dimension(wf%excited_state_specifications%n_singlet_states,1) :: energies
!
      end subroutine summary_excited_state_info_mlcc2
!
!
   end interface
!
!
   interface
!
!     -::- Jacobian transformation submodule -::-
!     :::::::::::::::::::::::::::::::::::::::::::
!
      module subroutine jacobian_mlcc2_transformation_mlcc2(wf, c_a_i, c_aibj)
!!
!!       Jacobian transformation (MLCC2)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, June 2017
!!
         implicit none
!
         class(mlcc2) :: wf 
!
!        Incoming vector c 
!
         real(dp), dimension(wf%n_v, wf%n_o) :: c_a_i  ! c_ai 
         real(dp), dimension(wf%n_x2am, 1)   :: c_aibj ! c_aibj  
!  
      end subroutine jacobian_mlcc2_transformation_mlcc2
!
!
      module subroutine jacobian_mlcc2_a1_mlcc2(wf, rho_a_i, c_a_i)
!!
!!       Jacobian tem A1
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, June 2017
!!
         implicit none
!
         class(mlcc2) :: wf 
!
!        Incoming vectors c and rho 
!
         real(dp), dimension(wf%n_v, wf%n_o) :: c_a_i 
         real(dp), dimension(wf%n_v, wf%n_o) :: rho_a_i 
!
      end subroutine jacobian_mlcc2_a1_mlcc2
!
!
      module subroutine jacobian_mlcc2_b1_mlcc2(wf, rho_a_i, c_ai_bj)
!!
!!       Jacobian tem B1
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, June 2017
!!
         implicit none
!  
         class(mlcc2) :: wf
!
         real(dp), dimension(:,:)            :: c_ai_bj 
         real(dp), dimension(wf%n_v, wf%n_o) :: rho_a_i
!
      end subroutine jacobian_mlcc2_b1_mlcc2
!
!
      module subroutine jacobian_mlcc2_a2_mlcc2(wf, rho_ai_bj, c_a_i)
!!
!!       Jacobian tem A2
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, June 2017
!!
         implicit none
!     
         class(mlcc2) :: wf
!  
         real(dp), dimension(:,:)            :: rho_ai_bj 
         real(dp), dimension(wf%n_v, wf%n_o) :: c_a_i 
!
      end subroutine jacobian_mlcc2_a2_mlcc2
!
!
      module subroutine jacobian_mlcc2_b2_mlcc2(wf, rho_ai_bj, c_ai_bj)
!!
!!       Jacobian tem B2
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, June 2017
!!
         implicit none
!  
         class(mlcc2) :: wf
!
         real(dp), dimension(:,:)   :: c_ai_bj 
         real(dp), dimension(:,:)   :: rho_ai_bj
!
!
      end subroutine jacobian_mlcc2_b2_mlcc2
!
!
      module subroutine cvs_rho_aibj_projection_mlcc2(wf, vec_aibj)
!!
!!       Rho projection for CVS (MLCC2),
!!       Written by Sarai D. Folkestad, Aug. 2017
!!
         implicit none
!
         class(mlcc2) :: wf
         real(dp), dimension(:, :) :: vec_aibj
!
      end subroutine cvs_rho_aibj_projection_mlcc2
!
!
   end interface
! 
!
contains
!
!  -::- MLCC2 initialization routine -::-
!  ::::::::::::::::::::::::::::::::::::::
!
   subroutine init_mlcc2(wf)
!!
!!
      implicit none
!
      class(mlcc2) :: wf
!
      integer(i15) :: i, j
!
      integer(i15) :: unit_input = -1
!
!     Set model name 
!
      wf%name = 'MLCC2'
!
!     Open input file eT.inp
!
      call generate_unit_identifier(unit_input)
      open(unit=unit_input, file='eT.inp', status='old', form='formatted')
      rewind(unit_input)
!
!     Read general specifications (memory and diskspace for calculation)
!
      call wf%general_specs_reader(unit_input)
!
!     Read MLCC info 
!
      call wf%mlcc_reader(unit_input)
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
      wf%implemented%ground_state = .true.
      wf%implemented%excited_state = .true.
      wf%implemented%core_excited_state = .true.
!
!     Read calculation tasks from input file eT.inp
!     
      call wf%calculation_reader(unit_input)
!
!     Read orbital info 
!
      call wf%read_orbital_info(unit_input)
!
!     Close input file
!
      close(unit_input)
!
!     Read Hartree-Fock info
!
      call wf%read_hf_info
!
!     Orbital partitioning - only if we have CCS region
!
      if (wf%mlcc_settings%CCS) then
!
         call wf%orbital_partitioning
!
      else
!
!        Do full space CC2 calculation
!
         write(unit_output,'(/t6,a50/)')'Full CC2 requested, orbital partitioning skipped.'
!
         wf%n_CC2_o = wf%n_o
         wf%n_CC2_v = wf%n_v
!
         wf%first_CC2_o = 1
         wf%first_CC2_v = 1         
!
      endif
!
!     Initialize amplitudes and associated attributes
!
      wf%n_t1am = (wf%n_o)*(wf%n_v)
!
      call wf%initialize_single_amplitudes
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
      call wf%destruct_single_amplitudes
!
   end subroutine init_mlcc2
!
!
!  :::::::::::::::::::::::::::::::::::::::::
!  -::- Class subroutines and functions -::- 
!  :::::::::::::::::::::::::::::::::::::::::
!
!  -::- Energy routine -::-
!  ::::::::::::::::::::::::
!
   subroutine calc_energy_mlcc2(wf)
!!
!!    Calculate Energy (MLCC2)
!!
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!    Calculates the MLCC2 energy, 
!!
!!    E_CC2 = E_HF + sum_aibj L_iajb*(s_ij^ab + t_i^a*t_j^b),
!!   
!!    with s_ij^ab = - g_aibj/(e_a + e_b - e_i - e_j) where 
!!    g_aibj are T1-transformed integrals.
!!    Batching over a.
!! 
!!
   implicit none
!
      class(mlcc2) :: wf
!
      logical :: debug = .false.
!
!     Integrals
!
      real(dp), dimension(:,:), allocatable :: L_IA_J
      real(dp), dimension(:,:), allocatable :: L_ai_J
      real(dp), dimension(:,:), allocatable :: g_IA_JB ! = g_aibj
!
!     s2 amplitudes
!
      real(dp), dimension(:,:), allocatable :: s_ia_jb ! = g_aibj/(e_a + e_b - e_i - e_j)
!
!     Batching variables
!  
      integer(i15) :: a_batch, a_first, a_last, a_length
      integer(i15) :: required, available, n_batch, batch_dimension, max_batch_length, offset
!
!     Indices
!
      integer(i15) :: a = 0, b = 0
      integer(i15) :: i = 0, j = 0
!
      integer(i15) :: ai = 0, bj = 0
      integer(i15) :: ia = 0, ib = 0, jb = 0, ja = 0
      integer(i15) :: IA_full = 0, IB_full = 0, JB_full = 0, JA_full = 0
!
      integer(i15) :: aibj = 0
!
!     ML variables
!
      integer(i15) :: n_active_o = 0, n_active_v = 0
      integer(i15) :: first_active_o ! first active occupied index 
      integer(i15) :: first_active_v ! first active virtual index
      integer(i15) :: last_active_o ! first active occupied index 
      integer(i15) :: last_active_v ! first active virtual index
!
!     :: t1 contribution ::
!
!     sum_aibj t_ai*t_bj*L_ia_jb
!
!
!     Allocate the Cholesky vector L_ia_J = L_ia^J and set to zero 
!
      call wf%mem%alloc(L_IA_J, (wf%n_o)*(wf%n_v), wf%n_J)
      L_IA_J = zero
!
!     Get the Cholesky vector L_ia_J 
!
      call wf%get_cholesky_ia(L_IA_J)
!
!     Allocate g_ia_jb = g_iajb and set it to zero
!
      call wf%mem%alloc(g_IA_JB, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      g_IA_JB = zero
!
!     Calculate the integrals g_ia_jb from the Cholesky vector L_ia_J 
!
      call dgemm('N','T',                                   &
                  (wf%n_o)*(wf%n_v),                        &
                  (wf%n_o)*(wf%n_v),                        &
                  wf%n_J,                                   &
                  one,                                      &
                  L_IA_J,                                   &
                  (wf%n_o)*(wf%n_v),                        &
                  L_IA_J,                                   &
                  (wf%n_o)*(wf%n_v),                        &
                  zero,                                     &
                  g_IA_JB,                                  &
                  (wf%n_o)*(wf%n_v))
!
!     Deallocate the Cholesky vector L_ia_J 
!
      call wf%mem%dealloc(L_IA_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!     Set the initial value of the energy 
!
      wf%energy = wf%scf_energy
!
!
!     Add the correlation energy E = E + sum_aibj (t_ij^ab + t_i^a t_j^b) L_iajb
!
      do I = 1, wf%n_o
         do A = 1, wf%n_v
!
            IA = index_two(I, A, wf%n_o)
!
            do J = 1, wf%n_o
!
               JA = index_two(J, A, wf%n_o)
!
              do B = 1, wf%n_v
! 
                  JB = index_two(J, B, wf%n_o)
                  IB = index_two(I, B, wf%n_o)
!
!                 Add the correlation energy 
!
                  wf%energy = wf%energy &
                            + (two*g_IA_JB(IA,JB) - g_IA_JB(JA, IB))*(wf%t1am(A,I))*(wf%t1am(B,J))
                           
!
               enddo
            enddo
         enddo
      enddo
!
!     :: s2 contribution ::
!
!     sum_aibj s_ai_bj*L_ia_jb
!
!     Set ML variables
!
      call wf%get_CC2_active_indices(first_active_o, first_active_v)
!
      call wf%get_CC2_n_active(n_active_o, n_active_v)
!
      last_active_o = first_active_o + n_active_o - 1
      last_active_v = first_active_v + n_active_v - 1 
!
!     Prepare for batching over index a
!  
      required = (2*(wf%n_v)*(wf%n_o)*(wf%n_J) &                            ! Needed for g_aibj  
                  + 2*((wf%n_v)**2)*(wf%n_J) &                              ! and 's2' amplitudes  
                  + 2*(wf%n_v)**2*(wf%n_o)**2)                              !
!         
      required = 4*required ! In words
      available = get_available()
!
      batch_dimension  = n_active_v ! Batch over the virtual index a
      max_batch_length = 0          ! Initilization of unset variables 
      n_batch          = 0
!
      call num_batch(required, available, max_batch_length, n_batch, batch_dimension)       
!
!     Loop over the number of a batches 
!
      do a_batch = 1, n_batch
!
!        For each batch, get the limits for the a index 
!
         call batch_limits(a_first, a_last, a_batch, max_batch_length, batch_dimension)
!
!        a is active index, and thus a_first and a_last must be displaced
!
         a_first  = a_first + (first_active_v - 1)
         a_last   = a_last  + (first_active_v - 1)
!
         if (a_last .gt. last_active_v) a_last = last_active_v
!
         a_length = a_last - a_first + 1 

!        :: Calculate cc2 doubles amplitudes ::
!
         call wf%mem%alloc(L_ai_J, n_active_o*n_active_v, wf%n_J)
         L_ai_J = zero
!    
         call wf%get_cholesky_ai(L_ai_J, first_active_v, last_active_v, first_active_o, last_active_o)
!
         call wf%mem%alloc(L_ia_J, (n_active_o)*(n_active_v), wf%n_J)
!
!        reorder and constrain L_bi_J
!
         do a = 1, n_active_v
            do i = 1, n_active_o
!
               ia = index_two(i, a, n_active_o)
               ai = index_two(a, i, n_active_v)
!
               do J = 1, wf%n_J
!
                  L_ia_J(ia, J) = L_ai_J(ai, J) 
!
               enddo
            enddo
         enddo
!
         call wf%mem%dealloc(L_ai_J, n_active_o*n_active_v, wf%n_J)
!
         call wf%mem%alloc(s_ia_jb, (n_active_o)*a_length, (n_active_o)*n_active_v)
!
         offset = index_two(1, a_first, n_active_o)
!
         call dgemm('N', 'T',                    &
                     (n_active_o)*a_length,      &
                     (n_active_o)*(n_active_v),  &
                     (wf%n_J),                   &
                     one,                        &
                     L_ia_J(offset,1),           &
                     (n_active_o)*(n_active_v),  &
                     L_ia_J,                     &
                     (n_active_o)*(n_active_v),  &
                     zero,                       &
                     s_ia_jb,                    &
                     (n_active_o)*a_length)
!
!
         call wf%mem%dealloc(L_ia_J, (n_active_o)*(n_active_v), wf%n_J)
!
!        Add the rest of the correlation energy E = E + sum_aibj (s_ij^ab ) L_iajb
!
         do a = 1, a_length
            do i = 1, n_active_o
!
               ia = index_two(i, a, n_active_o)
               IA_full = index_two(i + first_active_o - 1, a + a_first - 1, wf%n_o)
!
               do b = 1, n_active_v
!
                     IB_full = index_two(i + first_active_o - 1, b + first_active_v - 1, wf%n_o)
!
                  do j = 1, n_active_o
!
                     jb = index_two(j, b, n_active_o)
                     JB_full = index_two(j + first_active_o - 1, b + first_active_v - 1, wf%n_o)
                     JA_full = index_two(j + first_active_o - 1, a + a_first - 1, wf%n_o)
!
                     wf%energy = wf%energy + (two*g_ia_jb(IA_full, JB_full) - g_ia_jb(JA_full, IB_full))*((s_ia_jb(ia,jb))&
                                             /(wf%fock_diagonal(i + first_active_o - 1 ,1)&
                                              +wf%fock_diagonal(j + first_active_o - 1 ,1) &
                                             - wf%fock_diagonal(wf%n_o + a + a_first - 1 ,1)&
                                             - wf%fock_diagonal(wf%n_o + b + first_active_v - 1 ,1)))
!
                  enddo
               enddo
            enddo
         enddo
!  
         call wf%mem%dealloc(s_ia_jb, a_length*n_active_o, n_active_o*n_active_v)
!
      enddo ! End of batching
!
      call wf%mem%dealloc(g_ia_jb, wf%n_o*wf%n_v, wf%n_o*wf%n_v)
!
   end subroutine calc_energy_mlcc2
!
!  -::- ML helper routines -::-
!  ::::::::::::::::::::::::::::
!
   subroutine get_CC2_active_indices_mlcc2(wf, first_o, first_v)
!!
!!    Get CC2 active indices,
!!    Written by Sarai D. Folkestad, June 2017
!!
!!    Returns the first active occupied and virtual indices 
!!    of the active space.
!!
      implicit none
!
      class(mlcc2) :: wf
      integer(i15) :: first_o
      integer(i15) :: first_v
!
      first_o = wf%first_CC2_o
      first_v = wf%first_CC2_v
!
   end subroutine get_CC2_active_indices_mlcc2
!
!
   subroutine get_CC2_n_active_mlcc2(wf, n_active_o, n_active_v)
!!
!!    Get CC2 active indices,
!!    Written by Sarai D. Folkestad, June 2017
!!
!!    Returns the first active occupied and virtual indices 
!!    of the active space.
!!
      implicit none
!
      class(mlcc2) :: wf
      integer(i15) :: n_active_o
      integer(i15) :: n_active_v
!
      n_active_o = wf%n_CC2_o
      n_active_v = wf%n_CC2_v
!
   end subroutine get_CC2_n_active_mlcc2
!
!  -::- Amplitude routines -::-
!  ::::::::::::::::::::::::::::
!
   subroutine save_amplitudes_mlcc2(wf)
!!
!!    Save Amplitudes (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjøsntad, May 2017
!!
!!    Store the amplitudes to disk (T1AM)
!!
      implicit none 
!
      class(mlcc2) :: wf
!
      integer(i15) :: unit_x1am = -1 
      integer(i15) :: unit_x2am = -1 
!
      real(dp), dimension(:,:), allocatable :: s_ai_bj
      real(dp), dimension(:,:), allocatable :: s2am
!
      integer(i15) :: a = 0, i = 0, ai = 0, ia = 0, b = 0, j = 0, bj = 0, jb = 0
      integer(i15) :: aibj = 0
!
!     Open amplitude files
!
      call generate_unit_identifier(unit_x1am)
      open(unit_x1am, file='t1am', status='unknown', form='unformatted')
      rewind(unit_x1am)
!
      call generate_unit_identifier(unit_x2am)
      open(unit_x2am, file='x2am', status='unknown', form='unformatted')
      rewind(unit_x2am)
!
!     Write t1 amplitudes
!
      write(unit_x1am) wf%t1am
!
!     Construct s2 amplitudes
!
      call wf%mem%alloc(s_ai_bj, (wf%n_CC2_v)*(wf%n_CC2_o), (wf%n_CC2_v)*(wf%n_CC2_o))
      call wf%get_s2am(s_ai_bj)
!  
!     Reorder and pack in
!
      call wf%mem%alloc(s2am, (wf%n_CC2_v)*(wf%n_CC2_o)*((wf%n_CC2_v)*(wf%n_CC2_o)+1)/2, 1)
!
      do i = 1, wf%n_CC2_o
         do a = 1, wf%n_CC2_v
            ai = index_two(a, i, wf%n_CC2_v)
            do j = 1, wf%n_CC2_o
               do b = 1, wf%n_CC2_v
!
                  bj = index_two(b, j, wf%n_CC2_v)
!
                  aibj = index_packed(ai, bj)
!
                  s2am(aibj, 1) = s_ai_bj(ai, bj)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(s_ai_bj, (wf%n_CC2_v)*(wf%n_CC2_o), (wf%n_CC2_v)*(wf%n_CC2_o))
!
!     Write s2 amplitudes 
!
      write(unit_x2am) s2am
!
      call wf%mem%dealloc(s2am, (wf%n_CC2_v)*(wf%n_CC2_o)*((wf%n_CC2_v)*(wf%n_CC2_o)+1)/2, 1)
!
!     Close amplitude file
!
      close(unit_x1am)
      close(unit_x2am)
!
   end subroutine save_amplitudes_mlcc2
!
!
   subroutine read_amplitudes_mlcc2(wf)
!!
!!    Read Amplitudes (MLCC2)
!!    Written by Sarai D. Folkestad and Eirik F. Kjøsntad, May 2017
!!
!!    Reads the amplitudes from disk (T1AM)
!!
      implicit none 
!
      class(mlcc2) :: wf
!
      call wf%read_single_amplitudes
      call wf%read_cc2_double_amplitudes
!
   end subroutine read_amplitudes_mlcc2
!
   subroutine read_cc2_double_amplitudes_mlcc2(wf)
!!
!!    Read Amplitudes (MLCC2)
!!    Written by Sarai D. Folkestad and Eirik F. Kjøsntad, May 2017
!!
!!    Reads the amplitudes from disk (T1AM, S2AM)
!!
      implicit none 
!
      class(mlcc2) :: wf
!
      integer(i15) :: unit_x2am = -1 
!
      logical :: file_exists = .false.
!
!     Check to see whether file exists
!
      inquire(file='x2am',exist=file_exists)
!
      if (file_exists) then 
!
!        Open amplitude files if they exist
!
         call generate_unit_identifier(unit_x2am)
!
         open(unit_x2am, file='x2am', status='unknown', form='unformatted')
!
         rewind(unit_x2am)
!
!        Read from file & close
!
         wf%n_x2am = ((wf%n_CC2_v)*(wf%n_CC2_o))&
                   *((wf%n_CC2_v )*(wf%n_CC2_o)+1)/2 
!
         if (.not. allocated(wf%x2am)) call wf%mem%alloc(wf%x2am, wf%n_x2am, 1) 
         read(unit_x2am) wf%x2am
!
         close(unit_x2am)
!
      else
!
         write(unit_output,'(t3,a)') 'Error: amplitude files do not exist.'
         stop
!
      endif
!
   end subroutine read_cc2_double_amplitudes_mlcc2
!
   subroutine destruct_cc2_double_amplitudes_mlcc2(wf)
!!
!!
!!
      implicit none
!
      class(mlcc2) :: wf

!
      integer(i15) :: n_active_v, n_active_o 
!
      call wf%get_CC2_n_active(n_active_o, n_active_v)
       wf%n_x2am = ((n_active_v)*(n_active_o))&
                   *((n_active_v)*(n_active_o)+1)/2  
!
      if (allocated(wf%x2am)) call wf%mem%dealloc(wf%x2am, wf%n_x2am, 1)
!
   end subroutine destruct_cc2_double_amplitudes_mlcc2
!
!
   subroutine initialize_amplitudes_mlcc2(wf)
!!
!!    Initialize Amplitudes (MLCC2)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!    Allocates the amplitudes, sets them to zero.
!!
      implicit none 
!
      class(mlcc2) :: wf
!
      write(unit_output,*) 'Error: do not use initialize_amplitudes for ML'
      stop
!
   end subroutine initialize_amplitudes_mlcc2
!
!
   subroutine initialize_cc2_double_amplitudes_mlcc2(wf)
!!
!!    Initialize double amplitudes (MLCC2)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!    Allocates the doubles amplitudes, sets them to zero.
!!
      implicit none 
!
      class(mlcc2) :: wf

!
      integer(i15) :: n_active_v, n_active_o 
!
!     Allocate the doubles amplitudes and set to zero
!
      if (.not. allocated(wf%x2am)) then
!
         call wf%get_CC2_n_active(n_active_o, n_active_v)
         wf%n_x2am = ((n_active_v)*(n_active_o))&
                   *((n_active_v)*(n_active_o)+1)/2 
!
         call wf%mem%alloc(wf%x2am, wf%n_x2am, 1)
         wf%x2am = zero
!
      else
!
         write(unit_output,'(t3,a)') 'Warning: attempted to allocate and zero already allocated x2am'
!
      endif
!
   end subroutine initialize_cc2_double_amplitudes_mlcc2
!
!
   subroutine destruct_amplitudes_mlcc2(wf)
!!
!!    Destruct Amplitudes (MLCC2)
!!    Written by Sarai D. Folkestad and Eirik F. Kjøsntad, May 2017
!!
!!    Deallocates the (doubles) amplitudes.
!!
      implicit none
!
      class(mlcc2) :: wf
!
      write(unit_output,*) 'Error: do not use destruct_amplitudes for ML'
      stop
!
   end subroutine destruct_amplitudes_mlcc2
!
!
end module mlcc2_class
