module sccsd_class
!
!!
!!                Similarity constrained coupled cluster singles and doubles (SCCSD) class module                                 
!!                            Written by Eirik F. Kjønstad, June 2017         
!!                                                                           
!!
!!    This module contains the definition of the similarity constrained coupled cluster singles
!!    and doubles (SCCSD) wavefunction class. It is structured into four sections:
!!
!!       1. Modules used by the class: 
!!
!!             Basic utilities and the ancestor class (CCSD)
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
!!                - Excited state 
!!                - Omega 
!!                - Jacobian (right transformation)
!!                - Jacobian transpose (left transformation)
!!                - Input reader 
!!                - Metric (overlap between constrained states)
!!
!!             The interfaces shows incoming variables and their type, but contains 
!!             no information of the procedure itself. The procedure is shown in full 
!!             in the corresponding submodule. 
!!
!!       4. Class module routines (i.e., non-submodule procedures). These include
!!          the initialization and driver routines of the class, along with procedures that
!!          are not (yet, at least) easily gathered in a submodule.
!!         
!! 
!  ::::::::::::::::::::::::::::::::::::::
!  -::- 1. Modules used by the class -::-
!  ::::::::::::::::::::::::::::::::::::::
!
!  General tools
!
   use types
   use utils
   use workspace
   use input_output
   use scc_calculation_settings_class
!
!  Ancestor class module (CCSD)
!
   use ccsd_class
!
   implicit none 
!
!  ::::::::::::::::::::::::::::::::::::::::::
!  -::- 2. Definition of the SCCSD class -::-
!  ::::::::::::::::::::::::::::::::::::::::::
!
   type, extends(ccsd) :: sccsd
!
!     The triple amplitude that enforces similarity (= t_IJK^ABC)
!
      real(dp) :: triples = zero 
!
!     The associated triple excitation operator (= tau_IJK^ABC)
!
      integer(i15) :: I = 0
      integer(i15) :: J = 0
      integer(i15) :: K = 0
!
      integer(i15) :: A = 0
      integer(i15) :: B = 0
      integer(i15) :: C = 0
!
!     The excited states that are constrained to be similar
!
      integer(i15) :: state_A = 0
      integer(i15) :: state_B = 0
!
!     The generalized overlap between the two constrained states
!
      real(dp) :: overlap = zero
!
!     Calculation settings object for SCC calculations
!
      type(scc_calculation_settings) :: scc_settings
!
   contains
!
!
!     -::- Initialization routines -::-
!     ---------------------------------
!
      procedure :: init => init_sccsd
!
!
!     -::- Excited state submodule routine pointers -::-
!     --------------------------------------------------
!
      procedure :: excited_state_driver             => excited_state_driver_sccsd 
      procedure :: excited_state_intersection_cycle => excited_state_intersection_cycle_sccsd
      procedure :: eigenvector_controller           => eigenvector_controller_sccsd
!
      procedure :: ground_state_intersection_cycle     => ground_state_intersection_cycle_sccsd 
      procedure :: ground_state_eigenvector_controller => ground_state_eigenvector_controller_sccsd
!
!
!     -::- Omega submodule routine pointers -::-
!     ------------------------------------------
!
      procedure :: construct_omega => construct_omega_sccsd
!
      procedure :: omega_sccsd_a1  => omega_sccsd_a1_sccsd
!
!
!     -::- Jacobian submodule routine pointers -::-
!     ---------------------------------------------
!
      procedure :: jacobian_ccsd_transformation => jacobian_ccsd_transformation_sccsd 
!
      procedure :: jacobian_sccsd_a2 => jacobian_sccsd_a2_sccsd
      procedure :: jacobian_sccsd_b2 => jacobian_sccsd_b2_sccsd 
      procedure :: jacobian_sccsd_c2 => jacobian_sccsd_c2_sccsd
!
!
!     -::- Jacobian transpose submodule routine pointers -::-
!     -------------------------------------------------------
!
      procedure :: jacobian_transpose_ccsd_transformation => jacobian_transpose_ccsd_transformation_sccsd
!
      procedure :: jacobian_transpose_sccsd_a1 => jacobian_transpose_sccsd_a1_sccsd
      procedure :: jacobian_transpose_sccsd_b1 => jacobian_transpose_sccsd_b1_sccsd 
      procedure :: jacobian_transpose_sccsd_c1 => jacobian_transpose_sccsd_c1_sccsd
!
!
!     -::- Metric submodule routine pointers -::-
!     -------------------------------------------
!
      procedure :: metric_transformation => metric_transformation_sccsd
!
!     Helper routines for the metric transformation 
!
      procedure :: construct_q      => construct_q_sccsd       ! q_mu
      procedure :: Q_transformation => Q_transformation_sccsd  ! Q_mu,nu
      procedure :: S_transformation => S_transformation_sccsd  ! S_mu,nu  
!
!     Routine to calculate the overlap between the similarity constrained states
!
      procedure :: calc_overlap              => calc_overlap_sccsd
      procedure :: calc_ground_state_overlap => calc_ground_state_overlap_sccsd
!
!
!     -::- Input reader submodule routine pointers -::-
!     -------------------------------------------------
!
      procedure :: scc_reader => scc_reader_sccsd
!
!
!     -::- Other class routine pointers not located in submodules -::-
!     ----------------------------------------------------------------
!
      procedure :: read_triples => read_triples_sccsd
!
      procedure :: sccsd_diis => sccsd_diis_sccsd
!
   end type sccsd 
!
!
!  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  -::- 3. Interfaces to the submodules of the SCCSD class -::- 
!  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
   interface 
!
!
!     -::- Omega submodule interface -::-
!     -----------------------------------
!
      module subroutine construct_omega_sccsd(wf)
!!
!!       Construct Omega (SCCSD)
!!       Written by Eirik F. Kjønstad, June 2017
!!
         implicit none 
!
         class(sccsd) :: wf
!
      end subroutine construct_omega_sccsd
!
!
      module subroutine omega_sccsd_a1_sccsd(wf)
!!
!!       Omega SCCSD A1 
!!       Written by Eirik F. Kjønstad, June 2017
!!
         implicit none 
!
         class(sccsd) :: wf 
!
      end subroutine omega_sccsd_a1_sccsd
!
!
   end interface
!
!
   interface 
!
!
!     -::- Jacobian submodule interface -::-
!     --------------------------------------
!
      module subroutine jacobian_ccsd_transformation_sccsd(wf, c_a_i, c_aibj)
!!
!!       Jacobian transformation (SCCSD)
!!       Written by Eirik F. Kjønstad, May 2017
!!
         implicit none 
!
         class(sccsd) :: wf
!
         real(dp), dimension(wf%n_v, wf%n_o) :: c_a_i  
         real(dp), dimension(wf%n_t2am, 1)   :: c_aibj 
!
      end subroutine jacobian_ccsd_transformation_sccsd
!
!
      module subroutine jacobian_sccsd_a2_sccsd(wf, rho, X)
!!
!!       Jacobian SCCSD A2 
!!       Written by Eirik F. Kjønstad, June 2017
!!
         implicit none 
!
         class(sccsd) :: wf 
!
         real(dp), dimension(wf%n_v, wf%n_o) :: X                         
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho 
!
      end subroutine jacobian_sccsd_a2_sccsd
!
!
      module subroutine jacobian_sccsd_b2_sccsd(wf, rho, Y)
!!
!!       Jacobian SCCSD B2
!!       Written by Eirik F. Kjønstad, June 2017
!!
         implicit none 
!
         class(sccsd) :: wf 
!
         real(dp), dimension(wf%n_v, wf%n_o) :: Y 
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho 
!
      end subroutine jacobian_sccsd_b2_sccsd
!
!
      module subroutine jacobian_sccsd_c2_sccsd(wf, rho, Z)
!!
!!       Jacobian SCCSD C2
!!       Written by Eirik F. Kjønstad, June 2017
!!
         implicit none 
!
         class(sccsd) :: wf 
!
         real(dp), dimension((wf%n_v)**2, (wf%n_o)*(wf%n_v)) :: Z  
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho 
!
      end subroutine jacobian_sccsd_c2_sccsd
!
!
   end interface 
!
!
   interface 
!
!
!     -::- Jacobian transpose submodule interface -::-
!     ------------------------------------------------
!
      module subroutine jacobian_transpose_ccsd_transformation_sccsd(wf, b_a_i, b_aibj)
!!
!!       Jacobian transpose CCSD transformation (SCCSD)
!!       Written by Eirik F. Kjønstad, June 2017
!!
         implicit none 
!
         class(sccsd) :: wf 
!
         real(dp), dimension(wf%n_v, wf%n_o) :: b_a_i 
         real(dp), dimension(wf%n_t2am, 1)   :: b_aibj 
!
      end subroutine jacobian_transpose_ccsd_transformation_sccsd
!
!
      module subroutine jacobian_transpose_sccsd_a1_sccsd(wf, sigma, b2am, L)
!!
!!       Jacobian transpose SCCSD A1
!!       Written by Eirik F. Kjønstad, June 2017
!!
         implicit none 
!
         class(sccsd) :: wf 
!
         real(dp), dimension(wf%n_v, wf%n_o) :: sigma                      
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: L    
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: b2am 
!
      end subroutine jacobian_transpose_sccsd_a1_sccsd
!
!
      module subroutine jacobian_transpose_sccsd_b1_sccsd(wf, sigma, b2am, g)
!!
!!       Jacobian transpose SCCSD B1
!!       Written by Eirik F. Kjønstad, June 2017
!!
         implicit none 
!
         class(sccsd) :: wf 
!
         real(dp), dimension(wf%n_v, wf%n_o) :: sigma                       
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: g    
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: b2am 
!
      end subroutine jacobian_transpose_sccsd_b1_sccsd 
!
!
      module subroutine jacobian_transpose_sccsd_c1_sccsd(wf, sigma, b2am, g)
!!
!!       Jacobian transpose SCCSD C1
!!       Written by Eirik F. Kjønstad, June 2017
!!
         implicit none 
!
         class(sccsd) :: wf 
!
         real(dp), dimension(wf%n_v, wf%n_o) :: sigma                      
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: g    
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: b2am 
!
      end subroutine jacobian_transpose_sccsd_c1_sccsd  
!
!  
   end interface
!
!
   interface 
!
!
!     -::- Metric submodule interface -::-
!     ------------------------------------
!
      module subroutine calc_overlap_sccsd(wf)
!!
!!       Calculate Overlap (SCCSD)
!!       Written by Eirik F. Kjønstad, June 2017
!!
         implicit none 
!
         class(sccsd) :: wf 
!
      end subroutine calc_overlap_sccsd
!
!
      module subroutine calc_ground_state_overlap_sccsd(wf)
!!
!!       Calculate ground state overlap (SCCSD)
!!       Written by Eirik F. Kjønstad, Feb 2018
!!
         implicit none 
!
         class(sccsd) :: wf 
!
      end subroutine calc_ground_state_overlap_sccsd
!
!
      module subroutine metric_transformation_sccsd(wf, r_a_i, r_aibj)
!!
!!       Metric Transformation (SCCSD)
!!       Written by Eirik F. Kjønstad, June 2017
!!
         implicit none
!
         class(sccsd) :: wf 
!
         real(dp), dimension(wf%n_v, wf%n_o) :: r_a_i  ! r_ai 
         real(dp), dimension(wf%n_t2am, 1)   :: r_aibj ! r_aibj 
!
      end subroutine metric_transformation_sccsd
!
!
      module subroutine construct_q_sccsd(wf, q_a_i, q_aibj)
!!
!!       Construct q (SCCSD)
!!       Written by Eirik F. Kjønstad, June 2017
!!
         implicit none 
!
         class(sccsd) :: wf 
!
         real(dp), dimension(wf%n_v, wf%n_o) :: q_a_i  ! q_ai 
         real(dp), dimension(wf%n_t2am, 1)   :: q_aibj ! q_aibj 
!
      end subroutine construct_q_sccsd
!
!
      module subroutine Q_transformation_sccsd(wf, r_a_i, r_aibj, transpose)
!!
!!       Q transformation (SCCSD)
!!       Written by Eirik F. Kjønstad, June 2017
!!
         implicit none 
!
         class(sccsd) :: wf 
!
         real(dp), dimension(wf%n_v, wf%n_o) :: r_a_i  ! r_ai 
         real(dp), dimension(wf%n_t2am, 1)   :: r_aibj ! r_aibj
!
         logical :: transpose ! If true, transform by Q^T instead of Q  
!
      end subroutine Q_transformation_sccsd    
!
!
      module subroutine S_transformation_sccsd(wf, r_a_i, r_aibj)
!!
!!       S Transformation (SCCSD)
!!       Written by Eirik F. Kjønstad, June 2017
!!
         class(sccsd) :: wf
!
         real(dp), dimension(wf%n_v, wf%n_o) :: r_a_i 
         real(dp), dimension(wf%n_t2am, 1)   :: r_aibj
!
      end subroutine S_transformation_sccsd
!
!
   end interface 
!
!
   interface 
!
!
!     -::- Excited state submodule interface -::-
!     -------------------------------------------
!
      module subroutine excited_state_driver_sccsd(wf)
!!
!!       Excited state driver (SCCSD)
!!       Written by Eirik F. Kjønstad, June 2017
!!
         implicit none 
!
         class(sccsd) :: wf 
!
      end subroutine excited_state_driver_sccsd
!
!
      module subroutine excited_state_intersection_cycle_sccsd(wf, iteration)
!!
!!       Excited state intersection cycle (SCCSD)
!!       Written by Eirik F. Kjønstad, June 2017
!!
         implicit none 
!
         class(sccsd) :: wf 
!
         integer(i15) :: iteration 
!
      end subroutine excited_state_intersection_cycle_sccsd
!
!
      module subroutine ground_state_intersection_cycle_sccsd(wf, iteration)
!!
!!       Ground state intersection cycle (SCCSD)
!!       Written by Eirik F. Kjønstad, June 2017
!!
         implicit none 
!
         class(sccsd) :: wf 
!
         integer(i15) :: iteration 
!
      end subroutine ground_state_intersection_cycle_sccsd
!
!
      module subroutine eigenvector_controller_sccsd(wf, iteration)
!!
!!       Eigenvector controller (SCCSD)
!!       Written by Eirik F. Kjønstad, Dec 2017
!!
         implicit none 
!
         class(sccsd) :: wf 
!
         integer(i15), intent(in) :: iteration 
!
      end subroutine eigenvector_controller_sccsd
!
!
      module subroutine ground_state_eigenvector_controller_sccsd(wf, iteration)
!!
!!       Ground state eigenvector controller (SCCSD)
!!       Written by Eirik F. Kjønstad, Dec 2017
!!
         implicit none 
!
         class(sccsd) :: wf 
!
         integer(i15), intent(in) :: iteration 
!
      end subroutine ground_state_eigenvector_controller_sccsd
!
!
      module subroutine sccsd_diis_sccsd(wf, dt, t_dt, iteration)
!!
!!       SCCSD DIIS routine
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
         implicit none 
!
         class(sccsd) :: wf 
!
         integer(i15), intent(in) :: iteration
!
         real(dp) :: dt 
         real(dp) :: t_dt 
!
      end subroutine sccsd_diis_sccsd
!
!
!
   end interface
!
!
   interface
!
!
!     -::- Input reader submodule interface -::-
!     ------------------------------------------
!
      module subroutine scc_reader_sccsd(wf, unit_input)
!!
!!       SCC reader (SCCSD)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Dec 2017
!!
         implicit none
!
         integer(i15) :: unit_input
!
         class(sccsd) :: wf
!
      end subroutine scc_reader_sccsd
!
!
   end interface
!
!
contains
!
!
!  ::::::::::::::::::::::::::::::::::::::::::::
!  -::- 4. Class subroutines and functions -::- 
!  ::::::::::::::::::::::::::::::::::::::::::::
!
!
   subroutine init_sccsd(wf)
!!
!!    Initialize SCCSD object
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
      implicit none 
!
      class(sccsd) :: wf
!
      integer(i15) :: unit_input = -1
!
!     Set model name 
!
      wf%name = 'SCCSD'
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
!     Set implemented methods
!
      wf%implemented%ground_state       = .true.
      wf%implemented%excited_state      = .true.
      wf%implemented%core_excited_state = .false.
      wf%implemented%ionized_state      = .false.
      wf%implemented%multipliers        = .false.
!
!     Read calculation tasks from input file eT.inp
!     
      call wf%calculation_reader(unit_input)
!
!     Initialize the triples amplitude 
!
      wf%triples = zero
!
!     Standard output setting is minimalist
!
   !   wf%settings%print_level = 'minimal'
!
!     Read SCC specific information 
!
      call wf%scc_reader(unit_input)
!
!     Close input file
!
      close(unit_input)
!
!     The thresholds for energies and residuals should be, at least, as strict as the 
!     overlap threshold 
!
      if (wf%excited_state_specifications%energy_threshold .gt. wf%scc_settings%overlap_threshold .or. &
          wf%excited_state_specifications%residual_threshold .gt. wf%scc_settings%overlap_threshold) then 
!
         write(unit_output,'(/t3,a/)') 'Warning: energy or residual threshold too low; changed to equal overlap threshold.'
         flush(unit_output)
!
         wf%excited_state_specifications%energy_threshold   = wf%scc_settings%overlap_threshold
         wf%excited_state_specifications%residual_threshold = wf%scc_settings%overlap_threshold
!
      endif
!
!     Set maximum number of iterations to equal that of the excited state (temporary solution)
!
      wf%ground_state_specifications%max_iterations = wf%excited_state_specifications%max_iterations
!
!     Set ground state & response thresholds to equal the
!     excited state thresholds set by user (these should be consistent,
!     otherwise accuracy is limited by the lowest threshold)
!
      wf%ground_state_specifications%energy_threshold   = wf%excited_state_specifications%energy_threshold
      wf%ground_state_specifications%residual_threshold = wf%excited_state_specifications%residual_threshold
!
      wf%response_specifications%energy_threshold       = wf%excited_state_specifications%energy_threshold
      wf%response_specifications%residual_threshold     = wf%excited_state_specifications%residual_threshold
!
!     Print SCC-specific settings to output file 
!
      write(unit_output,'(/t3,a)') 'General settings for SCC calculation:'
!
      write(unit_output,'(/t6,a20,e9.2)') 'Overlap threshold:  ', wf%scc_settings%overlap_threshold
      write(unit_output,'(t6,a20,e9.2)')  'Energy thresholds:  ', wf%excited_state_specifications%energy_threshold
      write(unit_output,'(t6,a20,e9.2)')  'Residual thresholds:', wf%excited_state_specifications%residual_threshold
!
      write(unit_output, '(/t6,a,/t6,a,i1,a,i1,a)') &
                                     'Constraining the Jacobian to be nondefective in', &
                                     'the subspace spanned by the states ', wf%state_A, ' and ', wf%state_B, '.'
!
      write(unit_output,'(/t6,a)') 'Triple excitation:'
!
      write(unit_output,'(/t9,a,3i3)') 'Occupied indices:', wf%I, wf%J, wf%K
      write(unit_output,'(t9,a,3i3/)') 'Virtual indices: ', wf%A, wf%B, wf%C
!
!     Read Hartree-Fock info from SIRIUS
!
      call wf%read_hf_info
!
!     Read Cholesky AO integrals and transform to MO basis
!
      call wf%read_transform_cholesky 
!
!     Initialize (singles and doubles) amplitudes
!
      wf%n_t1am = (wf%n_o)*(wf%n_v) 
      wf%n_t2am = (wf%n_t1am)*(wf%n_t1am + 1)/2 
!
!     Set the number of parameters solved for in the ground state 
!     and excited state equations
!
      wf%n_parameters = wf%n_t1am + wf%n_t2am
!
!     Initialize the Fock matrix (allocate and construct given the initial amplitudes)
!
      call wf%initialize_single_amplitudes ! t1am = zero
      call wf%initialize_fock_matrix
      call wf%destruct_single_amplitudes
!
      if (wf%excited_state_specifications%restart) then
!
!           Set restarts to true 
!
            wf%ground_state_specifications%restart = .true.
            wf%excited_state_specifications%restart = .true.
            wf%response_specifications%restart = .true.
!
      endif
!
   end subroutine init_sccsd
!
!
   subroutine read_triples_sccsd(wf)
!!
!!    Read triples (SCCSD)
!!    Written by Eirik F. Kjønstad, June 2017
!!
!!    Reads the triple amplitude from file, if the file exists. 
!!    Note: the assumption is that if it exists it should and will be read.
!!
      implicit none 
!
      class(sccsd) :: wf 
!
      logical :: file_exists = .false.
!
      integer(i15) :: unit_triples = -1
!
!     Check to see whether file exists  
!
      inquire(file='triples',exist=file_exists)
!
      if (file_exists) then 
!
         call generate_unit_identifier(unit_triples)
         open(unit_triples, file='triples', status='unknown', form='unformatted')
         rewind(unit_triples)
         read(unit_triples) wf%triples 
!
         write(unit_output,'(/t3,a,/t3,a)') &
            'Found a triple amplitude on file.', 'Restarting with stored triple amplitude.'
!
         close(unit_triples)
!
      endif
!
   end subroutine read_triples_sccsd
!
!
end module sccsd_class
