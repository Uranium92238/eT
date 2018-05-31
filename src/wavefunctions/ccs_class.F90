module ccs_class
!
!!
!!
!!
!!                      Coupled cluster singles (CCS) class module
!!                Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!
!!
!!    This module contains the definition of the coupled cluster singles (CCS)
!!    wavefunction class. It is structured into four sections:
!!
!!       1. Modules used by the class:
!!
!!             Basic utilities and the ancestor class (HF)
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
!!                - Ground state
!!                - Excited state
!!                - Response
!!                - Input Reader
!!                - Cholesky
!!                - Integrals
!!                - Fock
!!                - Jacobian (right transformation)
!!                - Jacobian Transpose (left transformation)
!!                - Ionized State
!!                - CVS
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
!!
!
!
!  ::::::::::::::::::::::::::::::::::::::
!  -::- 1. Modules used by the class -::-
!  ::::::::::::::::::::::::::::::::::::::
!
!
!  General tools
!
   use types
   use utils
   use workspace
   use input_output
!
!  The ancestor class module (HF)
!
   use hf_class
!
   implicit none
!
!
!  ::::::::::::::::::::::::::::::::::::::::
!  -::- 2. Definition of the CCS class -::-
!  ::::::::::::::::::::::::::::::::::::::::
!
!
   type, extends(hf) :: ccs
!
!     Cluster amplitudes
!
      integer(i15) :: n_t1am       = 0 ! Number of single excitation amplitudes
      integer(i15) :: n_parameters = 0 ! Number of parameters in the wavefunction
!
      real(dp), dimension(:,:), allocatable :: t1am ! Single excitation amplitudes vector, t_i^a
!
!     The omega, or projection, vector < mu | exp(-T) H exp(T) | R >
!
      real(dp), dimension(:,:), allocatable :: omega1 ! Singles projection vector
!
!     The T1-transformed Fock matrix
!
      real(dp), dimension(:,:), allocatable :: fock_ij
      real(dp), dimension(:,:), allocatable :: fock_ia
      real(dp), dimension(:,:), allocatable :: fock_ai
      real(dp), dimension(:,:), allocatable :: fock_ab
!
!     Excitation energies
!
      real(dp), dimension(:,:), allocatable :: excitation_energies
!
   contains
!
!
!     -::- Initialization and driver routines -::-
!     --------------------------------------------
!
      procedure :: init => init_ccs ! Initialization of class
      procedure :: drv  => drv_ccs  ! Driver of class
!
!
!     -::- Ground state submodule routine pointers -::-
!     -------------------------------------------------
!
!     Driver and solver
!
      procedure :: ground_state_driver       => ground_state_driver_ccs
      procedure :: ground_state_solver       => ground_state_solver_ccs
!
!     Preparations and cleanup routines (before and after solver)
!
      procedure :: ground_state_preparations => ground_state_preparations_ccs
      procedure :: ground_state_cleanup      => ground_state_cleanup_ccs
!
!     Ground state restart routine
!
      procedure :: ground_state_restart      => ground_state_restart_ccs
!
!     DIIS component of solver, with helper routines
!
      procedure, non_overridable :: diis     => diis_ccs
!
      procedure :: new_amplitudes            => new_amplitudes_ccs
      procedure :: calc_ampeqs               => calc_ampeqs_ccs
      procedure :: calc_ampeqs_norm          => calc_ampeqs_norm_ccs
      procedure :: calc_quasi_Newton_singles => calc_quasi_Newton_singles_ccs
!
!
!     -::- Excited state submodule routine pointers -::-
!     --------------------------------------------------
!
!     Driver and solver
!
      procedure                  :: excited_state_driver      => excited_state_driver_ccs
      procedure, non_overridable :: excited_state_solver      => excited_state_solver_ccs
      procedure, non_overridable :: excited_state_solver_diis => excited_state_solver_diis_ccs
!
!     Preparations and cleanup routines (before and after solver)
!
      procedure :: excited_state_preparations => excited_state_preparations_ccs
      procedure :: excited_state_cleanup      => excited_state_cleanup_ccs
!
!     Excited state restart routine
!
      procedure :: excited_state_restart      => excited_state_restart_ccs
!
!     Helper routines for excited state solver
!
      procedure :: transform_trial_vectors          => transform_trial_vectors_ccs
      procedure :: calculate_orbital_differences    => calculate_orbital_differences_ccs
!
      procedure :: precondition_residual            => precondition_residual_ccs
      procedure :: precondition_residual_valence    => precondition_residual_valence_ccs
!
      procedure :: print_excited_state_info         => print_excited_state_info_ccs
      procedure :: print_excitation_vector          => print_excitation_vector_ccs
!
      procedure :: analyze_single_excitation_vector => analyze_single_excitation_vector_ccs
      procedure :: summary_excited_state_info       => summary_excited_state_info_ccs
!
!     Valence excited states specific routines
!
      procedure, non_overridable :: initialize_trial_vectors            => initialize_trial_vectors_ccs
      procedure, non_overridable :: initialize_trial_vectors_valence    => initialize_trial_vectors_valence_ccs
      procedure, non_overridable :: find_start_trial_indices            => find_start_trial_indices_ccs
      procedure, non_overridable :: trial_vectors_from_stored_solutions => trial_vectors_from_stored_solutions_ccs
!
!     Core excited states specific routines
!
      procedure                  :: precondition_residual_core          => precondition_residual_core_ccs
      procedure, non_overridable :: find_start_trial_indices_core       => find_start_trial_indices_core_ccs
      procedure, non_overridable :: find_core_mo                        => find_core_mo_ccs
      procedure, non_overridable :: initialize_trial_vectors_core       => initialize_trial_vectors_core_ccs
!
      procedure, non_overridable :: solve_reduced_eigenvalue_equation   => solve_reduced_eigenvalue_equation_ccs
      procedure, non_overridable :: construct_next_trial_vectors        => construct_next_trial_vectors_ccs
!
!
!     -::- Response submodule routine pointers -::-
!     ---------------------------------------------
!
!     Driver and solver
!
      procedure :: response_driver => response_driver_ccs
      procedure :: response_solver => response_solver_ccs
!
!     Helper routines for solver
!
      procedure :: response_preparations                 => response_preparations_ccs
      procedure :: response_cleanup                      => response_cleanup_ccs
      procedure :: initialize_response                   => initialize_response_ccs
      procedure :: solve_reduced_response_equation       => solve_reduced_response_equation_ccs
      procedure :: construct_reduced_matrix              => construct_reduced_matrix_ccs
      procedure :: construct_reduced_gradient            => construct_reduced_gradient_ccs
      procedure :: construct_next_response_trial_vectors => construct_next_response_trial_vectors_ccs
      procedure :: construct_gradient_vector             => construct_gradient_vector_ccs
!
!
!     -::- Input reader submodule routine pointers -::-
!     -------------------------------------------------
!
      procedure :: general_specs_reader      => general_specs_reader_ccs
      procedure :: calculation_reader        => calculation_reader_ccs
      procedure :: read_ground_state_specs   => read_ground_state_specs_ccs
      procedure :: read_excited_state_specs  => read_excited_state_specs_ccs
      procedure :: read_property_specs       => read_property_specs_ccs
!
!
!     -::- Cholesky submodule routine pointers -::-
!     ---------------------------------------------
!
      procedure, non_overridable :: get_cholesky_ij => get_cholesky_ij_ccs
      procedure, non_overridable :: get_cholesky_ia => get_cholesky_ia_ccs
      procedure, non_overridable :: get_cholesky_ai => get_cholesky_ai_ccs
      procedure, non_overridable :: get_cholesky_ab => get_cholesky_ab_ccs
!
!
!     -::- Integral submodule routine pointers -::-
!     ---------------------------------------------
!
!     Get integral routines
!
      procedure :: get_oo_oo => get_oo_oo_ccs
      procedure :: get_oo_ov => get_oo_ov_ccs
      procedure :: get_ov_oo => get_ov_oo_ccs
      procedure :: get_oo_vo => get_oo_vo_ccs
      procedure :: get_vo_oo => get_vo_oo_ccs
      procedure :: get_oo_vv => get_oo_vv_ccs
      procedure :: get_vv_oo => get_vv_oo_ccs
      procedure :: get_ov_ov => get_ov_ov_ccs
!
      procedure :: get_vo_vo => get_vo_vo_ccs
      procedure :: get_ov_vo => get_ov_vo_ccs
      procedure :: get_vo_ov => get_vo_ov_ccs
      procedure :: get_ov_vv => get_ov_vv_ccs
      procedure :: get_vv_ov => get_vv_ov_ccs
      procedure :: get_vo_vv => get_vo_vv_ccs
      procedure :: get_vv_vo => get_vv_vo_ccs
      procedure :: get_vv_vv => get_vv_vv_ccs
!
      procedure :: get_oo_oo_electronic_repulsion => get_oo_oo_electronic_repulsion_ccs
      procedure :: get_oo_ov_electronic_repulsion => get_oo_ov_electronic_repulsion_ccs
      procedure :: get_ov_oo_electronic_repulsion => get_ov_oo_electronic_repulsion_ccs
      procedure :: get_oo_vo_electronic_repulsion => get_oo_vo_electronic_repulsion_ccs
      procedure :: get_vo_oo_electronic_repulsion => get_vo_oo_electronic_repulsion_ccs
      procedure :: get_oo_vv_electronic_repulsion => get_oo_vv_electronic_repulsion_ccs
      procedure :: get_vv_oo_electronic_repulsion => get_vv_oo_electronic_repulsion_ccs
      procedure :: get_ov_ov_electronic_repulsion => get_ov_ov_electronic_repulsion_ccs
!
      procedure :: get_vo_vo_electronic_repulsion => get_vo_vo_electronic_repulsion_ccs
      procedure :: get_ov_vo_electronic_repulsion => get_ov_vo_electronic_repulsion_ccs
      procedure :: get_vo_ov_electronic_repulsion => get_vo_ov_electronic_repulsion_ccs
      procedure :: get_ov_vv_electronic_repulsion => get_ov_vv_electronic_repulsion_ccs
      procedure :: get_vv_ov_electronic_repulsion => get_vv_ov_electronic_repulsion_ccs
      procedure :: get_vo_vv_electronic_repulsion => get_vo_vv_electronic_repulsion_ccs
      procedure :: get_vv_vo_electronic_repulsion => get_vv_vo_electronic_repulsion_ccs
      procedure :: get_vv_vv_electronic_repulsion => get_vv_vv_electronic_repulsion_ccs
!
!     Routines to store, read, and T1-transform electronic repulsion integrals
!
      procedure :: store_vv_vv_electronic_repulsion    => store_vv_vv_electronic_repulsion_ccs
!
      procedure :: read_vv_vv_electronic_repulsion     => read_vv_vv_electronic_repulsion_ccs
      procedure :: t1_transform_vv_vv                  => t1_transform_vv_vv_ccs
!
      procedure :: store_t1_vv_vv_electronic_repulsion => store_t1_vv_vv_electronic_repulsion_ccs
      procedure :: store_t1_vo_ov_electronic_repulsion => store_t1_vo_ov_electronic_repulsion_ccs
      procedure :: store_t1_vv_vo_electronic_repulsion => store_t1_vv_vo_electronic_repulsion_ccs
      procedure :: store_t1_vv_ov_electronic_repulsion => store_t1_vv_ov_electronic_repulsion_ccs
!
      procedure :: read_t1_vv_vo_electronic_repulsion  => read_t1_vv_vo_electronic_repulsion_ccs
      procedure :: read_t1_vv_vv_electronic_repulsion  => read_t1_vv_vv_electronic_repulsion_ccs
      procedure :: read_t1_vo_ov_electronic_repulsion  => read_t1_vo_ov_electronic_repulsion_ccs
      procedure :: read_t1_vv_ov_electronic_repulsion  => read_t1_vv_ov_electronic_repulsion_ccs
!
      procedure :: get_vvvv_required_mem               => get_vvvv_required_mem_ccs
      procedure :: get_vvvo_required_mem               => get_vvvo_required_mem_ccs
      procedure :: get_vvov_required_mem               => get_vvov_required_mem_ccs
      procedure :: get_vvoo_required_mem               => get_vvoo_required_mem_ccs
!
!
!     -::- Fock submodule routine pointers -::-
!     -----------------------------------------
!
      procedure, non_overridable :: initialize_fock_matrix => initialize_fock_matrix_ccs
      procedure, non_overridable :: construct_fock         => construct_fock_ccs
      procedure, non_overridable :: one_electron_t1        => one_electron_t1_ccs
!
!
!     -::- Jacobian submodule routine pointers -::-
!     ---------------------------------------------
!
      procedure :: jacobian_ccs_transformation      => jacobian_ccs_transformation_ccs
!
      procedure, non_overridable :: jacobian_ccs_a1 => jacobian_ccs_a1_ccs
      procedure, non_overridable :: jacobian_ccs_b1 => jacobian_ccs_b1_ccs
!
      procedure :: jacobi_test => jacobi_test_ccs ! A debug routine
!
!
!     -::- Jacobian transpose submodule routine pointers -::-
!     -------------------------------------------------------
!
      procedure :: jacobian_transpose_ccs_transformation      => jacobian_transpose_ccs_transformation_ccs
!
      procedure, non_overridable :: jacobian_transpose_ccs_a1 => jacobian_transpose_ccs_a1_ccs
      procedure, non_overridable :: jacobian_transpose_ccs_b1 => jacobian_transpose_ccs_b1_ccs
!
!
!     -::- Ionized state submodule routine pointers -::-
!     --------------------------------------------------
!
      procedure :: ionized_state_driver                           => ionized_state_driver_ccs
      procedure :: initialize_trial_vectors_core_ionization       => initialize_trial_vectors_core_ionization_ccs
      procedure :: initialize_trial_vectors_valence_ionization    => initialize_trial_vectors_valence_ionization_ccs
      procedure :: precondition_residual_valence_ionization       => precondition_residual_valence_ionization_ccs
      procedure :: ionization_residual_projection                 => ionization_residual_projection_ccs
      procedure :: ionization_rho_a_i_projection                  => ionization_rho_a_i_projection_ccs
      procedure :: precondition_residual_core_ionization          => precondition_residual_core_ionization_ccs
!
!
!     -::- CVS submodule routine pointers -::-
!     ----------------------------------------
!
      procedure :: cvs_rho_a_i_projection             => cvs_rho_a_i_projection_ccs
      procedure :: cvs_residual_projection            => cvs_residual_projection_ccs
!
!
!     -::- Other class routine pointers not located in submodules -::-
!     ----------------------------------------------------------------
!
!     Routine to construct omega, or projection, vector
!
      procedure :: construct_omega => construct_omega_ccs
      procedure :: omega_ccs_a1    => omega_ccs_a1_ccs
!
!     Ground state energy calculation routine
!
      procedure :: calc_energy => calc_energy_ccs
!
!     Routine to construct eta, or right projection, vector
!
      procedure :: construct_eta => construct_eta_ccs
!
!     Routines to allocate amplitudes
!
      procedure :: initialize_amplitudes        => initialize_single_amplitudes_ccs ! Note: points to the same
      procedure :: initialize_single_amplitudes => initialize_single_amplitudes_ccs ! routine in CCS, for obvious reasons
!
!     Routines to deallocate amplitudes
!
      procedure :: destruct_amplitudes        => destruct_single_amplitudes_ccs ! Note: points to the same
      procedure :: destruct_single_amplitudes => destruct_single_amplitudes_ccs ! routine in CCS, for obvious reasons
!
!     Routine to save the amplitudes to disk
!
      procedure :: save_amplitudes        => save_amplitudes_ccs
!
!     Routines to read the amplitudes from disk (and allocate if necessary)
!
      procedure :: read_amplitudes        => read_single_amplitudes_ccs ! Note: points to the same
      procedure :: read_single_amplitudes => read_single_amplitudes_ccs ! routine in CCS, for obvious reasons
!
!     Routines to allocate and deallocate omega
!
      procedure :: initialize_omega => initialize_omega_ccs ! Allocate and zero projection vector
      procedure :: destruct_omega   => destruct_omega_ccs   ! Deallocate projection vector
!
!
   end type ccs
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
!     -::- Ground state submodule interface -::-
!     ------------------------------------------
!
      module subroutine ground_state_driver_ccs(wf)
!!
!!       Ground State Driver (CCS)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
         implicit none
!
         class(ccs) :: wf
!
      end subroutine ground_state_driver_ccs
!
!
      module subroutine ground_state_preparations_ccs(wf)
!!
!!       Ground State Preparations (CCS)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
         class(ccs) :: wf
!
      end subroutine ground_state_preparations_ccs
!
!
      module subroutine ground_state_solver_ccs(wf)
!!
!!       Ground State Solver
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
         implicit none
!
         class(ccs) :: wf
!
      end subroutine ground_state_solver_ccs
!
!
      module subroutine ground_state_cleanup_ccs(wf)
!!
!!       Ground State Cleanup (CCS)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
         class(ccs) :: wf
!
      end subroutine ground_state_cleanup_ccs
!
!
      module subroutine ground_state_restart_ccs(wf)
!!
!!       Ground state restart (CCS)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Dec 2017
!!
         class(ccs) :: wf
!
      end subroutine ground_state_restart_ccs
!
!
      module subroutine calc_ampeqs_ccs(wf)
!!
!!       Calculate Amplitude Equations (CCS)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
         implicit none
!
         class(ccs) :: wf
!
      end subroutine calc_ampeqs_ccs
!
!
      module subroutine calc_ampeqs_norm_ccs(wf, ampeqs_norm)
!!
!!       Calculate Amplitude Equations Norm (CCS)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
         implicit none
!
         class(ccs) :: wf
         real(dp)   :: ampeqs_norm
!
      end subroutine calc_ampeqs_norm_ccs
!
!
      module subroutine new_amplitudes_ccs(wf, diis_ground_state)
!!
!!       New Amplitudes (CCS)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         class(diis) :: diis_ground_state
!
      end subroutine new_amplitudes_ccs
!
!
      module subroutine calc_quasi_Newton_singles_ccs(wf,dt)
!!
!!       Calculate quasi-Newton Singles (CCS)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         real(dp), dimension(wf%n_parameters, 1) :: dt
!
      end subroutine calc_quasi_Newton_singles_ccs
!
!
      module subroutine diis_ccs(wf,dt,t_dt)
!!
!!       DIIS (CCS)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
         implicit none
!
         class(ccs), intent(in)   :: wf
!
         real(dp), dimension(wf%n_parameters, 1) :: dt
         real(dp), dimension(wf%n_parameters, 1) :: t_dt
!
      end subroutine diis_ccs
!
!
   end interface
!
!
   interface
!
!     -::- Excited state submodule interface -::-
!     -------------------------------------------
!
      module subroutine initialize_trial_vectors_ccs(wf)
!!
!!       Initialize Trial Vectors (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
         implicit none
!
         class(ccs) :: wf
!
      end subroutine initialize_trial_vectors_ccs
!
!
      module subroutine find_start_trial_indices_ccs(wf, index_list)
!!
!!       Find Start Trial Indices (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
         implicit none
!
         class(ccs) :: wf
         integer(i15), dimension(wf%excited_state_specifications%n_singlet_states,1), intent(inout) :: index_list
!
      end subroutine find_start_trial_indices_ccs
!
!
      module subroutine transform_trial_vectors_ccs(wf, first_trial, last_trial)
!!
!!       Transform trial vectors (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
         implicit none
!
         class(ccs) :: wf
         integer(i15), intent(in) :: first_trial, last_trial
!
      end subroutine transform_trial_vectors_ccs
!
!
      module subroutine excited_state_driver_ccs(wf)
!!
!!       Excited state driver (CCS)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
         implicit none
!
         class(ccs) :: wf
!
      end subroutine excited_state_driver_ccs
!
!
      module subroutine excited_state_preparations_ccs(wf)
!!
!!       Excited State Preparations (CCS)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
         class(ccs) :: wf
!
      end subroutine excited_state_preparations_ccs
!
!
      module subroutine excited_state_cleanup_ccs(wf)
!!
!!       Excited State Cleanup (CCS)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
         class(ccs) :: wf
!
      end subroutine excited_state_cleanup_ccs
!
!
      module subroutine excited_state_solver_ccs(wf)
!!
!!       Excited state solver
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
         implicit none
!
         class(ccs) :: wf
!
      end subroutine excited_state_solver_ccs
!
!
      module subroutine excited_state_solver_diis_ccs(wf)
!!
!!       Excited state diis solver
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2018
!!
         implicit none
!
         class(ccs) :: wf
!
      end subroutine excited_state_solver_diis_ccs
!
!
      module subroutine excited_state_restart_ccs(wf)
!!
!!       Excited state restart (CCS)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Dec 2017
!!
         class(ccs) :: wf
!
      end subroutine excited_state_restart_ccs
!
!
      module subroutine solve_reduced_eigenvalue_equation_ccs(wf, eigenvalues_Re, eigenvalues_Im, &
                                                               solution_vectors_reduced, reduced_dim, n_new_trials)
!!
!!       Solve Reduced Eigenvalue Equation
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         integer(i15) :: reduced_dim, n_new_trials
!
         real(dp), dimension(wf%excited_state_specifications%n_singlet_states,1) :: eigenvalues_Re
         real(dp), dimension(wf%excited_state_specifications%n_singlet_states,1) :: eigenvalues_Im
!
         real(dp), dimension(reduced_dim, wf%excited_state_specifications%n_singlet_states) :: solution_vectors_reduced
!
      end subroutine solve_reduced_eigenvalue_equation_ccs
!
!
     module subroutine construct_next_trial_vectors_ccs(wf, eigenvalues_Re, eigenvalues_Im, &
                                                   solution_vectors_reduced, &
                                                   reduced_dim, n_new_trials)
!!
!!       Construct Next Trial Vectors
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         real(dp), dimension(wf%excited_state_specifications%n_singlet_states,1) :: eigenvalues_Re
         real(dp), dimension(wf%excited_state_specifications%n_singlet_states,1) :: eigenvalues_Im
!
         integer(i15) :: reduced_dim
         integer(i15) :: n_new_trials
         real(dp), dimension(reduced_dim, wf%excited_state_specifications%n_singlet_states) :: solution_vectors_reduced
!
      end subroutine construct_next_trial_vectors_ccs
!
!
      module subroutine initialize_trial_vectors_valence_ccs(wf)
!!
!!       Initialize Trial Vectors Valence
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2017
!!
         implicit none
!
         class(ccs) :: wf
!
      end subroutine initialize_trial_vectors_valence_ccs
!
!
!
      module subroutine initialize_trial_vectors_core_ccs(wf)
!!
!!       Initialize trial vectors, for CVS calculation
!!       Written by Sarai D. Folkestad, Aug 2017
!!
         implicit none
!
         class(ccs) :: wf
!
      end subroutine initialize_trial_vectors_core_ccs
!
!
      module subroutine trial_vectors_from_stored_solutions_ccs(wf)
!!
!!    Trial Vectors from Old Solutions
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, June 2017
!!
      implicit none
!
      class(ccs) :: wf
!
      end subroutine trial_vectors_from_stored_solutions_ccs
!
!
      module subroutine find_start_trial_indices_core_ccs(wf, index_list)
!!
!!       Find Start Trial Indices Core
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2017
!!
         implicit none
!
         class(ccs) :: wf
         integer(i15), dimension(wf%excited_state_specifications%n_singlet_states,1), intent(inout) :: index_list
!
      end subroutine find_start_trial_indices_core_ccs
!
!
      module subroutine find_core_mo_ccs(wf)
!!
!!       Find Core MO
!!       Written by Sarai D. Folkestad, Aug 2017
!!
         implicit none
!
         class(ccs) :: wf
!
      end subroutine find_core_mo_ccs
!
!
      module subroutine calculate_orbital_differences_ccs(wf,orbital_diff)
!!
!!       Calculate Orbital Differences
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
         implicit none
!
         class(ccs) :: wf
         real(dp), dimension(wf%n_parameters, 1) :: orbital_diff
!
      end subroutine calculate_orbital_differences_ccs
!
!
      module subroutine initialize_excited_states_ccs(wf)
!!
!!       Initialize Excited States
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, June 2017
!!
         implicit none
!
         class(ccs) :: wf
!
      end subroutine initialize_excited_states_ccs
!
!
      module subroutine precondition_residual_ccs(wf, residual)
!!
!!       Precondition Residual
!!       Written by Sarai D. Folkestad, Aug 2017
!!
         implicit none
!
         class(ccs) :: wf
         real(dp), dimension(wf%n_parameters ,1) :: residual
!
      end subroutine precondition_residual_ccs
!
!
      module subroutine precondition_residual_valence_ccs(wf, residual)
!!
!!       Precondition Residual Valence
!!       Written by Sarai D. Folkestad, Aug 2017
!!
         implicit none
!
         class(ccs) :: wf
         real(dp), dimension(wf%n_parameters ,1) :: residual
!
      end subroutine precondition_residual_valence_ccs
!
!
      module subroutine precondition_residual_core_ccs(wf, residual)
!!
!!       Precondition Residual Core
!!       Written by Sarai D. Folkestad, Aug 2017
!!
         implicit none
!
         class(ccs) :: wf
         real(dp), dimension(wf%n_parameters ,1) :: residual
!
!
      end subroutine precondition_residual_core_ccs
!
!
      module subroutine print_excited_state_info_ccs(wf)
!!
!!       Print Excited State Information
!!       Written by Sarai D. Folkestad, Aug 2017
!!
         implicit none
!
         class(ccs) :: wf
!
      end subroutine print_excited_state_info_ccs
!
!
      module subroutine print_excitation_vector_ccs(wf, vec, unit_id)
!!
!!       Print Excitation Vector
!!       Written by Sarai D. Folkestad, Aug 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         real(dp), dimension(wf%n_parameters, 1) :: vec
!
         integer(i15) :: unit_id
!
      end subroutine print_excitation_vector_ccs
!
!
      module subroutine analyze_single_excitation_vector_ccs(wf, vec, n, sorted_short_vec, index_list)
!!
!!       Analyze Single Excitation Vector (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Nov 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         real(dp), dimension(wf%n_o*wf%n_v, 1) :: vec
!
         integer(i15) :: a = 0, i = 0, ai = 0
!
         integer(i15) :: n    ! Number of elements wanted
!
         real(dp), dimension(n, 1)    :: sorted_short_vec
!
         integer(i15), dimension(n, 2) ::index_list
!
      end subroutine analyze_single_excitation_vector_ccs
!
!
      module subroutine summary_excited_state_info_ccs(wf, energies)
!!
!!       Summary Excited State Info (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Nov 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         real(dp), dimension(wf%excited_state_specifications%n_singlet_states,1) :: energies
!
      end subroutine summary_excited_state_info_ccs
!
!
   end interface
!
!
   interface
!
!     -::- Response submodule interface -::-
!     --------------------------------------
!
      module subroutine response_driver_ccs(wf)
!!
!!       Response Driver (CCS)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
         implicit none
!
         class(ccs) :: wf
!
      end subroutine response_driver_ccs
!
!
      module subroutine response_solver_ccs(wf)
!!
!!       Response Solver (CCS)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
         implicit none
!
         class(ccs) :: wf
!
      end subroutine response_solver_ccs
!
!
      module subroutine response_preparations_ccs(wf)
!!
!!       Response calculation preparations (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Dec 2017
!!
         implicit none
!
         class(ccs) :: wf
!
      end subroutine response_preparations_ccs
!
!
      module subroutine response_cleanup_ccs(wf)
!!
!!       Response cleanup (CCS)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Dec 2017
!!
         implicit none
!
         class(ccs) :: wf
!
      end subroutine response_cleanup_ccs
!
!
      module subroutine initialize_response_ccs(wf)
!!
!!       Initialize Response (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, June 2017
!!
         implicit none
!
         class(ccs) :: wf
!
      end subroutine initialize_response_ccs
!
!
      module subroutine solve_reduced_response_equation_ccs(wf, solution_vector_reduced, reduced_dim, n_new_trials)
!!
!!       Solve Reduced Response Equation
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, June 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         integer(i15) :: reduced_dim, n_new_trials
!
         real(dp), dimension(reduced_dim, 1) :: solution_vector_reduced
!
      end subroutine solve_reduced_response_equation_ccs
!
!
      module subroutine construct_reduced_matrix_ccs(wf, A_reduced, reduced_dim, n_new_trials)
!!
!!       Construct Reduced Matrix (CCS)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         integer(i15) :: reduced_dim, n_new_trials
!
         real(dp), dimension(reduced_dim, reduced_dim) :: A_reduced
!
      end subroutine construct_reduced_matrix_ccs
!
!
      module subroutine construct_reduced_gradient_ccs(wf, F_reduced, reduced_dim, n_new_trials)
!!
!!       Construct Reduced Gradient (CCS)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         integer(i15) :: reduced_dim, n_new_trials
!
         real(dp), dimension(reduced_dim, 1) :: F_reduced
!
      end subroutine construct_reduced_gradient_ccs
!
!
      module subroutine construct_gradient_vector_ccs(wf)
!!
!!       Construct Gradient Vector (CCS)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
         implicit none
!
         class(ccs) :: wf
!
      end subroutine construct_gradient_vector_ccs
!
!
      module subroutine construct_next_response_trial_vectors_ccs(wf, solution_vector_reduced, reduced_dim, n_new_trials)
!!
!!       Construct Next Response Trial Vectors (CCS)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         integer(i15) :: reduced_dim, n_new_trials
!
         real(dp), dimension(reduced_dim, 1) :: solution_vector_reduced ! X_reduced
!
      end subroutine construct_next_response_trial_vectors_ccs
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
      module subroutine general_specs_reader_ccs(wf, unit_input)
!!
!!       General Specifications Reader
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Nov. 2017
!!
         implicit none
!
         class(ccs)   :: wf
         integer(i15) :: unit_input
!
         integer :: memory = 0, disk_space = 0
!
      end subroutine general_specs_reader_ccs
!
!
      module subroutine calculation_reader_ccs(wf, unit_input)
!!
!!       Calculation Reader
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Nov. 2017
!!
         implicit none
!
         integer(i15) :: unit_input
!
         class(ccs) :: wf
!
      end subroutine calculation_reader_ccs
!
!
      module subroutine read_ground_state_specs_ccs(wf, unit_input)
!!
!!       Read Ground State Specifications
!!       Written by Eirik F. Kjønstad and Sarai Dery Folkestad, Nov. 2017
!!
         implicit none
!
         integer :: unit_input
!
         class(ccs) :: wf
!
      end subroutine read_ground_state_specs_ccs
!
!
      module subroutine read_excited_state_specs_ccs(wf, unit_input)
!!
!!       Read Excited State Specifications
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Nov. 2017
!!
         implicit none
!
         integer(i15)      :: unit_input
!
         class(ccs)        :: wf
!
      end subroutine read_excited_state_specs_ccs
!
!
      module subroutine read_property_specs_ccs(wf, unit_input)
!!
!!       Read Property Specifications
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Nov. 2017
!!
         implicit none
!
         integer(i15)      :: unit_input
!
         class(ccs)        :: wf
!
      end subroutine read_property_specs_ccs
!
!
   end interface
!
!
   interface
!
!
!     -::- Cholesky submodule interface -::-
!     --------------------------------------
!
      module subroutine get_cholesky_ij_ccs(wf, L_ij_J, i_first, i_last, j_first, j_last)
!!
!!       Get Cholesky IJ
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
         implicit none
!
         class(ccs)               :: wf
         integer(i15), optional   :: i_first, j_first ! First index (can differ from 1 when batching or for mlcc)
         integer(i15), optional   :: i_last, j_last   ! Last index (can differ from n_o when batching or for mlcc)
         real(dp), dimension(:,:) :: L_ij_J
!
      end subroutine get_cholesky_ij_ccs
!
!
      module subroutine get_cholesky_ia_ccs(wf, L_ia_J, i_first, i_last, a_first, a_last)
!!
!!       Get Cholesky IA
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
         implicit none
!
         class(ccs)               :: wf
         integer(i15), optional   :: i_first, a_first   ! First index (can differ from 1 when batching or for mlcc)
         integer(i15), optional   :: i_last, a_last     ! Last index (can differ from n_v/n_o when batching or for mlcc)
         real(dp), dimension(:,:) :: L_ia_J
!
      end subroutine get_cholesky_ia_ccs
!
!
      module subroutine get_cholesky_ai_ccs(wf, L_ai_J, a_first, a_last, i_first, i_last)
!!
!!       Get Cholesky AI
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
         implicit none
!
         class(ccs)               :: wf
         integer(i15), optional   :: i_first, a_first     ! First index (can differ from 1 when batching or for mlcc)
         integer(i15), optional   :: i_last, a_last       ! Last index (can differ from n_o/n_v when batching or for mlcc)
         real(dp), dimension(:,:) :: L_ai_J
!
      end subroutine get_cholesky_ai_ccs
!
!
      module subroutine get_cholesky_ab_ccs(wf, L_ab_J, a_first, a_last, b_first, b_last)
!!
!!       Get Cholesky AB
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
         implicit none
!
         class(ccs)               :: wf
         integer(i15), intent(in) :: a_first, b_first   ! First index (can differ from 1 when batching or for mlcc)
         integer(i15), intent(in) :: a_last, b_last     ! Last index (can differ from n_v when batching or for mlcc)
         real(dp), dimension(((b_last - b_first + 1)*(a_last - a_first + 1)), wf%n_J) :: L_ab_J ! L_ab^J
!
      end subroutine get_cholesky_ab_ccs
!
!
   end interface
!
!

   interface
!
!     -::- Integral submodule interface -::-
!     --------------------------------------
!
      module subroutine store_vv_vv_electronic_repulsion_ccs(wf)
!!
!!       Store vvvv Electronic Repulsion
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
         implicit none
!
         class(ccs) :: wf
!
      end subroutine store_vv_vv_electronic_repulsion_ccs
!
!
      module subroutine store_t1_vv_vv_electronic_repulsion_ccs(wf)
!!
!!       Store t1 vvvv Electronic Repulsion Integrals
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
         implicit none
!
         class(ccs) :: wf
!
      end subroutine store_t1_vv_vv_electronic_repulsion_ccs
!
!
      module subroutine store_t1_vo_ov_electronic_repulsion_ccs(wf)
!!
!!       Store t1 voov Electronic Repulsion Integrals
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
         implicit none
!
         class(ccs) :: wf
!
      end subroutine store_t1_vo_ov_electronic_repulsion_ccs
!
!
      module subroutine store_t1_vv_ov_electronic_repulsion_ccs(wf)
!!
!!       Store t1 voov Electronic Repulsion Integrals
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
         implicit none
!
         class(ccs) :: wf
!
      end subroutine store_t1_vv_ov_electronic_repulsion_ccs
!
!
      module subroutine store_t1_vv_vo_electronic_repulsion_ccs(wf)
!!
!!       Store t1 vvvo Electronic Repulsion Integrals
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
         implicit none
!
         class(ccs) :: wf
!
      end subroutine store_t1_vv_vo_electronic_repulsion_ccs
!
!
      module subroutine read_vv_vv_electronic_repulsion_ccs(wf, x_vv_vv,    &
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!!
!!       Read vvvv Electronic Repulsion
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         real(dp), dimension(:,:) :: x_vv_vv
!
         integer(i15) :: index1_first, index1_last
         integer(i15) :: index2_first, index2_last
         integer(i15) :: index3_first, index3_last
         integer(i15) :: index4_first, index4_last
!
      end subroutine read_vv_vv_electronic_repulsion_ccs
!
!
      module subroutine read_t1_vv_vv_electronic_repulsion_ccs(wf, x_vv_vv,    &
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!!
!!       Read t1 vvvv Electronic Repulsion
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         real(dp), dimension(:,:) :: x_vv_vv
!
         integer(i15) :: index1_first, index1_last
         integer(i15) :: index2_first, index2_last
         integer(i15) :: index3_first, index3_last
         integer(i15) :: index4_first, index4_last
!
      end subroutine read_t1_vv_vv_electronic_repulsion_ccs
!
!
      module subroutine read_t1_vo_ov_electronic_repulsion_ccs(wf, x_vo_ov,    &
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!!
!!       Read t1 voov Electronic Repulsion
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         real(dp), dimension(:,:) :: x_vo_ov
!
         integer(i15) :: index1_first, index1_last
         integer(i15) :: index2_first, index2_last
         integer(i15) :: index3_first, index3_last
         integer(i15) :: index4_first, index4_last
!
      end subroutine read_t1_vo_ov_electronic_repulsion_ccs
!
!
      module subroutine read_t1_vv_vo_electronic_repulsion_ccs(wf, x_vv_vo,    &
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!!
!!       Read t1 vvvo Electronic Repulsion Integrals
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
         implicit none
!
         class(ccs) :: wf
!
!        Integral
!
         real(dp), dimension(:,:) :: x_vv_vo
!
         integer(i15) :: index1_first, index1_last
         integer(i15) :: index2_first, index2_last
         integer(i15) :: index3_first, index3_last
         integer(i15) :: index4_first, index4_last
!
      end subroutine read_t1_vv_vo_electronic_repulsion_ccs
!
!
      module subroutine read_t1_vv_ov_electronic_repulsion_ccs(wf, x_vv_ov,    &
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!!
!!       Read t1 vvov Electronic Repulsion Integrals
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
         implicit none
!
         class(ccs) :: wf
!
!        Integral
!
         real(dp), dimension(:,:) :: x_vv_ov
!
         integer(i15) :: index1_first, index1_last
         integer(i15) :: index2_first, index2_last
         integer(i15) :: index3_first, index3_last
         integer(i15) :: index4_first, index4_last
!
      end subroutine read_t1_vv_ov_electronic_repulsion_ccs
!
!
      module subroutine get_oo_oo_ccs(wf, integral_type, x_oo_oo,    &
                                          index1_first, index1_last, &
                                          index2_first, index2_last, &
                                          index3_first, index3_last, &
                                          index4_first, index4_last)
!!
!!       Get x_oo,oo integral (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         character(len=40) :: integral_type
!
         real(dp), dimension(:, :) :: x_oo_oo
!
         integer(i15), optional :: index1_first, index1_last
         integer(i15), optional :: index2_first, index2_last
         integer(i15), optional :: index3_first, index3_last
         integer(i15), optional :: index4_first, index4_last
!
      end subroutine get_oo_oo_ccs
!
!
      module subroutine get_oo_ov_ccs(wf, integral_type, x_oo_ov,    &
                                          index1_first, index1_last, &
                                          index2_first, index2_last, &
                                          index3_first, index3_last, &
                                          index4_first, index4_last)
!!
!!       Get x_oo,ov integral (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         character(len=40) :: integral_type
!
         real(dp), dimension(:,:) :: x_oo_ov
!
         integer(i15), optional :: index1_first, index1_last
         integer(i15), optional :: index2_first, index2_last
         integer(i15), optional :: index3_first, index3_last
         integer(i15), optional :: index4_first, index4_last
!
      end subroutine get_oo_ov_ccs
!
!
      module subroutine get_ov_oo_ccs(wf, integral_type, x_ov_oo,    &
                                          index1_first, index1_last, &
                                          index2_first, index2_last, &
                                          index3_first, index3_last, &
                                          index4_first, index4_last)
!!
!!       Get x_ov,oo integral (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         character(len=40) :: integral_type
!
         real(dp), dimension(:,:) :: x_ov_oo
!
         integer(i15), optional :: index1_first, index1_last
         integer(i15), optional :: index2_first, index2_last
         integer(i15), optional :: index3_first, index3_last
         integer(i15), optional :: index4_first, index4_last
!
      end subroutine get_ov_oo_ccs
!
!
      module subroutine get_oo_vo_ccs(wf, integral_type, x_oo_vo,    &
                                          index1_first, index1_last, &
                                          index2_first, index2_last, &
                                          index3_first, index3_last, &
                                          index4_first, index4_last)
!!
!!       Get x_oo,vo integral (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         character(len=40) :: integral_type
!
         real(dp), dimension(:,:) :: x_oo_vo
!
         integer(i15), optional :: index1_first, index1_last
         integer(i15), optional :: index2_first, index2_last
         integer(i15), optional :: index3_first, index3_last
         integer(i15), optional :: index4_first, index4_last
!
      end subroutine get_oo_vo_ccs
!
!
      module subroutine get_vo_oo_ccs(wf, integral_type, x_vo_oo,    &
                                          index1_first, index1_last, &
                                          index2_first, index2_last, &
                                          index3_first, index3_last, &
                                          index4_first, index4_last)
!!
!!       Get x_vo,oo integral (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         character(len=40) :: integral_type
!
         real(dp), dimension(:,:) :: x_vo_oo
!
         integer(i15), optional :: index1_first, index1_last
         integer(i15), optional :: index2_first, index2_last
         integer(i15), optional :: index3_first, index3_last
         integer(i15), optional :: index4_first, index4_last
!
      end subroutine get_vo_oo_ccs
!
!
      module subroutine get_oo_vv_ccs(wf, integral_type, x_oo_vv,    &
                                          index1_first, index1_last, &
                                          index2_first, index2_last, &
                                          index3_first, index3_last, &
                                          index4_first, index4_last)
!!
!!       Get x_oo,vv integral (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         character(len=40) :: integral_type
!
         real(dp), dimension(:,:) :: x_oo_vv
!
         integer(i15), optional :: index1_first, index1_last
         integer(i15), optional :: index2_first, index2_last
         integer(i15), optional :: index3_first, index3_last
         integer(i15), optional :: index4_first, index4_last
!
      end subroutine get_oo_vv_ccs
!
!
      module subroutine get_vv_oo_ccs(wf, integral_type, x_vv_oo,    &
                                          index1_first, index1_last, &
                                          index2_first, index2_last, &
                                          index3_first, index3_last, &
                                          index4_first, index4_last)
!!
!!       Get x_vv,oo integral (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         character(len=40) :: integral_type
!
         real(dp), dimension(:,:) :: x_vv_oo
!
         integer(i15), optional :: index1_first, index1_last
         integer(i15), optional :: index2_first, index2_last
         integer(i15), optional :: index3_first, index3_last
         integer(i15), optional :: index4_first, index4_last
!
      end subroutine get_vv_oo_ccs
!
!
      module subroutine get_ov_ov_ccs(wf, integral_type, x_ov_ov,    &
                                          index1_first, index1_last, &
                                          index2_first, index2_last, &
                                          index3_first, index3_last, &
                                          index4_first, index4_last)
!!
!!       Get x_ov,ov integral (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         character(len=40) :: integral_type
!
         real(dp), dimension(:,:) :: x_ov_ov
!
         integer(i15), optional :: index1_first, index1_last
         integer(i15), optional :: index2_first, index2_last
         integer(i15), optional :: index3_first, index3_last
         integer(i15), optional :: index4_first, index4_last
!
      end subroutine get_ov_ov_ccs
!
!
      module subroutine get_vo_vo_ccs(wf, integral_type, x_vo_vo,    &
                                          index1_first, index1_last, &
                                          index2_first, index2_last, &
                                          index3_first, index3_last, &
                                          index4_first, index4_last)
!!
!!       Get x_vo,vo integral (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         character(len=40) :: integral_type
!
         real(dp), dimension(:,:) :: x_vo_vo
!
         integer(i15), optional :: index1_first, index1_last
         integer(i15), optional :: index2_first, index2_last
         integer(i15), optional :: index3_first, index3_last
         integer(i15), optional :: index4_first, index4_last
!
      end subroutine get_vo_vo_ccs
!
!
      module subroutine get_ov_vo_ccs(wf, integral_type, x_ov_vo,    &
                                          index1_first, index1_last, &
                                          index2_first, index2_last, &
                                          index3_first, index3_last, &
                                          index4_first, index4_last)
!!
!!       Get x_ov,vo integral (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         character(len=40) :: integral_type
!
         real(dp), dimension(:,:) :: x_ov_vo
!
         integer(i15), optional :: index1_first, index1_last
         integer(i15), optional :: index2_first, index2_last
         integer(i15), optional :: index3_first, index3_last
         integer(i15), optional :: index4_first, index4_last
!
      end subroutine get_ov_vo_ccs
!
!
      module subroutine get_vo_ov_ccs(wf, integral_type, x_vo_ov,    &
                                          index1_first, index1_last, &
                                          index2_first, index2_last, &
                                          index3_first, index3_last, &
                                          index4_first, index4_last)
!!
!!       Get x_vo,ov integral (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         character(len=40) :: integral_type
!
         real(dp), dimension(:,:) :: x_vo_ov
!
         integer(i15), optional :: index1_first, index1_last
         integer(i15), optional :: index2_first, index2_last
         integer(i15), optional :: index3_first, index3_last
         integer(i15), optional :: index4_first, index4_last
!
      end subroutine get_vo_ov_ccs
!
!
      module subroutine get_ov_vv_ccs(wf, integral_type, x_ov_vv,    &
                                          index1_first, index1_last, &
                                          index2_first, index2_last, &
                                          index3_first, index3_last, &
                                          index4_first, index4_last)
!!
!!       Get x_ov,vv integral (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         character(len=40) :: integral_type
!
         real(dp), dimension(:,:) :: x_ov_vv
!
         integer(i15), optional :: index1_first, index1_last
         integer(i15), optional :: index2_first, index2_last
         integer(i15), optional :: index3_first, index3_last
         integer(i15), optional :: index4_first, index4_last
!
      end subroutine get_ov_vv_ccs
!
!
      module subroutine get_vv_ov_ccs(wf, integral_type, x_vv_ov,    &
                                          index1_first, index1_last, &
                                          index2_first, index2_last, &
                                          index3_first, index3_last, &
                                          index4_first, index4_last)
!!
!!       Get x_vv,ov integral (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         character(len=40) :: integral_type
!
         real(dp), dimension(:,:) :: x_vv_ov
!
         integer(i15), optional :: index1_first, index1_last
         integer(i15), optional :: index2_first, index2_last
         integer(i15), optional :: index3_first, index3_last
         integer(i15), optional :: index4_first, index4_last
!
      end subroutine get_vv_ov_ccs
!
!
      module subroutine get_vo_vv_ccs(wf, integral_type, x_vo_vv,    &
                                          index1_first, index1_last, &
                                          index2_first, index2_last, &
                                          index3_first, index3_last, &
                                          index4_first, index4_last)
!!
!!       Get x_vo,vv integral (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         character(len=40) :: integral_type
!
         real(dp), dimension(:,:) :: x_vo_vv
!
         integer(i15), optional :: index1_first, index1_last
         integer(i15), optional :: index2_first, index2_last
         integer(i15), optional :: index3_first, index3_last
         integer(i15), optional :: index4_first, index4_last
!
      end subroutine get_vo_vv_ccs
!
!
      module subroutine get_vv_vo_ccs(wf, integral_type, x_vv_vo,    &
                                          index1_first, index1_last, &
                                          index2_first, index2_last, &
                                          index3_first, index3_last, &
                                          index4_first, index4_last)
!!
!!       Get x_vv,vo integral (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         character(len=40) :: integral_type
!
         real(dp), dimension(:,:) :: x_vv_vo
!
         integer(i15), optional :: index1_first, index1_last
         integer(i15), optional :: index2_first, index2_last
         integer(i15), optional :: index3_first, index3_last
         integer(i15), optional :: index4_first, index4_last
!
      end subroutine get_vv_vo_ccs
!
!
      module subroutine get_vv_vv_ccs(wf, integral_type, x_vv_vv,    &
                                          index1_first, index1_last, &
                                          index2_first, index2_last, &
                                          index3_first, index3_last, &
                                          index4_first, index4_last)
!!
!!       Get x_vv,vv integral (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         character(len=40) :: integral_type
!
         real(dp), dimension(:,:) :: x_vv_vv
!
         integer(i15), optional :: index1_first, index1_last
         integer(i15), optional :: index2_first, index2_last
         integer(i15), optional :: index3_first, index3_last
         integer(i15), optional :: index4_first, index4_last
!
      end subroutine get_vv_vv_ccs
!
!
      module subroutine get_oo_oo_electronic_repulsion_ccs(wf, x_oo_oo,             &
                                                         index1_first, index1_last, &
                                                         index2_first, index2_last, &
                                                         index3_first, index3_last, &
                                                         index4_first, index4_last)
!!
!!       Get oooo Electronic Repulsion (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         real(dp), dimension(:, :) :: x_oo_oo
!
         integer(i15) :: index1_first, index1_last
         integer(i15) :: index2_first, index2_last
         integer(i15) :: index3_first, index3_last
         integer(i15) :: index4_first, index4_last
!
         real(dp), dimension(:,:), allocatable :: L_ij_J, L_kl_J
!
         integer(i15) :: length_1 = 0, length_2 = 0, length_3 = 0, length_4 = 0
!
      end subroutine get_oo_oo_electronic_repulsion_ccs
!
!
      module subroutine get_oo_ov_electronic_repulsion_ccs(wf, x_oo_ov,             &
                                                         index1_first, index1_last, &
                                                         index2_first, index2_last, &
                                                         index3_first, index3_last, &
                                                         index4_first, index4_last)
!!
!!       Get ooov Electronic Repulsion (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         real(dp), dimension(:, :) :: x_oo_ov
!
         integer(i15) :: index1_first, index1_last
         integer(i15) :: index2_first, index2_last
         integer(i15) :: index3_first, index3_last
         integer(i15) :: index4_first, index4_last
!
         real(dp), dimension(:,:), allocatable :: L_ij_J, L_ka_J
!
         integer(i15) :: length_1 = 0, length_2 = 0, length_3 = 0, length_4 = 0
!
      end subroutine get_oo_ov_electronic_repulsion_ccs
!
!
      module subroutine get_ov_oo_electronic_repulsion_ccs(wf, x_ov_oo,             &
                                                         index1_first, index1_last, &
                                                         index2_first, index2_last, &
                                                         index3_first, index3_last, &
                                                         index4_first, index4_last)
!!
!!       Get ovoo Electronic Repulsion (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         real(dp), dimension(:, :) :: x_ov_oo
!
         integer(i15) :: index1_first, index1_last
         integer(i15) :: index2_first, index2_last
         integer(i15) :: index3_first, index3_last
         integer(i15) :: index4_first, index4_last
!
         real(dp), dimension(:,:), allocatable :: L_ia_J, L_jk_J
!
         integer(i15) :: length_1 = 0, length_2 = 0, length_3 = 0, length_4 = 0
!
      end subroutine get_ov_oo_electronic_repulsion_ccs
!
!
      module subroutine get_oo_vo_electronic_repulsion_ccs(wf, x_oo_vo,             &
                                                         index1_first, index1_last, &
                                                         index2_first, index2_last, &
                                                         index3_first, index3_last, &
                                                         index4_first, index4_last)
!!
!!       Get oovo Electronic Repulsion (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         real(dp), dimension(:, :) :: x_oo_vo
!
         integer(i15) :: index1_first, index1_last
         integer(i15) :: index2_first, index2_last
         integer(i15) :: index3_first, index3_last
         integer(i15) :: index4_first, index4_last
!
         real(dp), dimension(:,:), allocatable :: L_ij_J, L_ak_J
!
         integer(i15) :: length_1 = 0, length_2 = 0, length_3 = 0, length_4 = 0
!
      end subroutine get_oo_vo_electronic_repulsion_ccs
!
!
      module subroutine get_oo_vv_electronic_repulsion_ccs(wf, x_oo_vv,             &
                                                         index1_first, index1_last, &
                                                         index2_first, index2_last, &
                                                         index3_first, index3_last, &
                                                         index4_first, index4_last)
!!
!!       Get oovv Electronic Repulsion (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         real(dp), dimension(:, :) :: x_oo_vv
!
         integer(i15) :: index1_first, index1_last
         integer(i15) :: index2_first, index2_last
         integer(i15) :: index3_first, index3_last
         integer(i15) :: index4_first, index4_last
!
         real(dp), dimension(:,:), allocatable :: L_ij_J, L_ab_J
!
         integer(i15) :: length_1 = 0, length_2 = 0, length_3 = 0, length_4 = 0
!
      end subroutine get_oo_vv_electronic_repulsion_ccs
!
!
      module subroutine get_vv_oo_electronic_repulsion_ccs(wf, x_vv_oo,             &
                                                         index1_first, index1_last, &
                                                         index2_first, index2_last, &
                                                         index3_first, index3_last, &
                                                         index4_first, index4_last)
!!
!!       Get vvoo Electronic Repulsion (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         real(dp), dimension(:, :) :: x_vv_oo
!
         integer(i15) :: index1_first, index1_last
         integer(i15) :: index2_first, index2_last
         integer(i15) :: index3_first, index3_last
         integer(i15) :: index4_first, index4_last
!
         real(dp), dimension(:,:), allocatable :: L_ab_J, L_ij_J
!
         integer(i15) :: length_1 = 0, length_2 = 0, length_3 = 0, length_4 = 0
!
      end subroutine get_vv_oo_electronic_repulsion_ccs
!
!
      module subroutine get_ov_ov_electronic_repulsion_ccs(wf, x_ov_ov,             &
                                                         index1_first, index1_last, &
                                                         index2_first, index2_last, &
                                                         index3_first, index3_last, &
                                                         index4_first, index4_last)
!!
!!       Get ovov Electronic Repulsion (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         real(dp), dimension(:, :) :: x_ov_ov
!
         integer(i15) :: index1_first, index1_last
         integer(i15) :: index2_first, index2_last
         integer(i15) :: index3_first, index3_last
         integer(i15) :: index4_first, index4_last
!
         real(dp), dimension(:,:), allocatable :: L_ia_J, L_jb_J
!
         integer(i15) :: length_1 = 0, length_2 = 0, length_3 = 0, length_4 = 0
!
      end subroutine get_ov_ov_electronic_repulsion_ccs
!
!
      module subroutine get_vo_oo_electronic_repulsion_ccs(wf, x_vo_oo,             &
                                                         index1_first, index1_last, &
                                                         index2_first, index2_last, &
                                                         index3_first, index3_last, &
                                                         index4_first, index4_last)
!!
!!       Get vooo Electronic Repulsion (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         real(dp), dimension(:, :) :: x_vo_oo
!
         integer(i15) :: index1_first, index1_last
         integer(i15) :: index2_first, index2_last
         integer(i15) :: index3_first, index3_last
         integer(i15) :: index4_first, index4_last
!
         real(dp), dimension(:,:), allocatable :: L_ai_J, L_jk_J
!
         integer(i15) :: length_1 = 0, length_2 = 0, length_3 = 0, length_4 = 0
!
      end subroutine get_vo_oo_electronic_repulsion_ccs
!
!
      module subroutine get_vo_vo_electronic_repulsion_ccs(wf, x_vo_vo, &
                                          index1_first, index1_last,    &
                                          index2_first, index2_last,    &
                                          index3_first, index3_last,    &
                                          index4_first, index4_last)
!!
!!       Get vovo Electronic Repulsion (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         real(dp), dimension(:,:) :: x_vo_vo
!
         integer(i15) :: index1_first, index1_last
         integer(i15) :: index2_first, index2_last
         integer(i15) :: index3_first, index3_last
         integer(i15) :: index4_first, index4_last
!
      end subroutine get_vo_vo_electronic_repulsion_ccs
!
!
      module subroutine get_ov_vo_electronic_repulsion_ccs(wf, x_ov_vo, &
                                          index1_first, index1_last,    &
                                          index2_first, index2_last,    &
                                          index3_first, index3_last,    &
                                          index4_first, index4_last)
!!
!!       Get ovvo Electronic Repulsion (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         real(dp), dimension(:,:) :: x_ov_vo
!
         integer(i15) :: index1_first, index1_last
         integer(i15) :: index2_first, index2_last
         integer(i15) :: index3_first, index3_last
         integer(i15) :: index4_first, index4_last
!
      end subroutine get_ov_vo_electronic_repulsion_ccs
!
!
      module subroutine get_vo_ov_electronic_repulsion_ccs(wf, x_vo_ov, &
                                          index1_first, index1_last,    &
                                          index2_first, index2_last,    &
                                          index3_first, index3_last,    &
                                          index4_first, index4_last)
!!
!!       Get voov Electronic Repulsion (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         real(dp), dimension(:,:) :: x_vo_ov
!
         integer(i15) :: index1_first, index1_last
         integer(i15) :: index2_first, index2_last
         integer(i15) :: index3_first, index3_last
         integer(i15) :: index4_first, index4_last
!
      end subroutine get_vo_ov_electronic_repulsion_ccs
!
!
      module subroutine get_ov_vv_electronic_repulsion_ccs(wf, x_ov_vv, &
                                          index1_first, index1_last,    &
                                          index2_first, index2_last,    &
                                          index3_first, index3_last,    &
                                          index4_first, index4_last)
!!
!!       Get ovvv Electronic Repulsion (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         real(dp), dimension(:,:) :: x_ov_vv
!
         integer(i15) :: index1_first, index1_last
         integer(i15) :: index2_first, index2_last
         integer(i15) :: index3_first, index3_last
         integer(i15) :: index4_first, index4_last
!
      end subroutine get_ov_vv_electronic_repulsion_ccs
!
!
      module subroutine get_vv_ov_electronic_repulsion_ccs(wf, x_vv_ov, &
                                          index1_first, index1_last,    &
                                          index2_first, index2_last,    &
                                          index3_first, index3_last,    &
                                          index4_first, index4_last)
!!
!!       Get vvov Electronic Repulsion (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         real(dp), dimension(:,:) :: x_vv_ov
!
         integer(i15) :: index1_first, index1_last
         integer(i15) :: index2_first, index2_last
         integer(i15) :: index3_first, index3_last
         integer(i15) :: index4_first, index4_last
!
      end subroutine get_vv_ov_electronic_repulsion_ccs
!
!
      module subroutine get_vo_vv_electronic_repulsion_ccs(wf, x_vo_vv, &
                                          index1_first, index1_last,    &
                                          index2_first, index2_last,    &
                                          index3_first, index3_last,    &
                                          index4_first, index4_last)
!!
!!       Get vovv Electronic Repulsion (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         real(dp), dimension(:,:) :: x_vo_vv
!
         integer(i15) :: index1_first, index1_last
         integer(i15) :: index2_first, index2_last
         integer(i15) :: index3_first, index3_last
         integer(i15) :: index4_first, index4_last
!
      end subroutine get_vo_vv_electronic_repulsion_ccs
!
!
      module subroutine get_vv_vo_electronic_repulsion_ccs(wf, x_vv_vo, &
                                          index1_first, index1_last,    &
                                          index2_first, index2_last,    &
                                          index3_first, index3_last,    &
                                          index4_first, index4_last)
!!
!!       Get vvvo Electronic Repulsion (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         real(dp), dimension(:,:) :: x_vv_vo
!
         integer(i15) :: index1_first, index1_last
         integer(i15) :: index2_first, index2_last
         integer(i15) :: index3_first, index3_last
         integer(i15) :: index4_first, index4_last
!
      end subroutine get_vv_vo_electronic_repulsion_ccs
!
!
      module subroutine get_vv_vv_electronic_repulsion_ccs(wf, x_vv_vv, &
                                          index1_first, index1_last,    &
                                          index2_first, index2_last,    &
                                          index3_first, index3_last,    &
                                          index4_first, index4_last)
!!
!!       Get vvvv Electronic Repulsion (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         real(dp), dimension(:,:) :: x_vv_vv
!
         integer(i15) :: index1_first, index1_last
         integer(i15) :: index2_first, index2_last
         integer(i15) :: index3_first, index3_last
         integer(i15) :: index4_first, index4_last
!
      end subroutine get_vv_vv_electronic_repulsion_ccs
!
!
      module subroutine t1_transform_vv_vv_ccs(wf, g_vv_vv,                &
                                             index1_first, index1_last, &
                                             index2_first, index2_last, &
                                             index3_first, index3_last, &
                                             index4_first, index4_last)
!!
!!       T1 Transformation of g_vv_vv Integrals (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Oct 2017.
!!
         implicit none
!
         class(ccs) :: wf
!
         real(dp), dimension(:, :) :: g_vv_vv
!
         integer(i15) :: index1_first, index1_last
         integer(i15) :: index2_first, index2_last
         integer(i15) :: index3_first, index3_last
         integer(i15) :: index4_first, index4_last
!
         end subroutine t1_transform_vv_vv_ccs
!
!
      module function get_vvvv_required_mem_ccs(wf, dim_1, dim_2, dim_3, dim_4)
!!
!!       Get vvvv required memory (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Dec 2017
!!
         implicit none
!
         class(ccs), intent(in)              :: wf
!
         integer(i15), intent(in), optional  :: dim_1, dim_2, dim_3, dim_4
!
         integer(i15) :: get_vvvv_required_mem_ccs
!
      end function get_vvvv_required_mem_ccs
!
!
      module function get_vvvo_required_mem_ccs(wf, dim_1, dim_2, dim_3, dim_4)
!
!        Get vvvo required memory (CCS)
!        Written by Eirik F. Kjønstad and Sarai D. Folkestad, Dec 2017
!
         implicit none
!
         class(ccs), intent(in)              :: wf
!
         integer(i15), intent(in), optional  :: dim_1, dim_2, dim_3, dim_4
!
         integer(i15) :: get_vvvo_required_mem_ccs
!
      end function get_vvvo_required_mem_ccs
!
!
      module function get_vvov_required_mem_ccs(wf, dim_1, dim_2, dim_3, dim_4)
!
!        Get vvov required memory (CCS)
!        Written by Eirik F. Kjønstad and Sarai D. Folkestad, Dec 2017
!
         implicit none

         class(ccs), intent(in)              :: wf
!
         integer(i15), intent(in), optional  :: dim_1, dim_2, dim_3, dim_4
!
         integer(i15) :: get_vvov_required_mem_ccs
!
      end function get_vvov_required_mem_ccs
!
!
      module function get_vvoo_required_mem_ccs(wf, dim_1, dim_2, dim_3, dim_4)
!
!        Get vvoo required memory (CCS)
!        Written by Eirik F. Kjønstad and Sarai D. Folkestad, Dec 2017
!
         implicit none
!
         class(ccs), intent(in)              :: wf
!
         integer(i15), intent(in), optional  :: dim_1, dim_2, dim_3, dim_4
!
         integer(i15) :: get_vvoo_required_mem_ccs
!
      end function get_vvoo_required_mem_ccs
!
   end interface
!
!
   interface
!
!
!    -::- Fock submodule interface -::-
!    ----------------------------------
!
      module subroutine initialize_fock_matrix_ccs(wf)
!!
!!       Initialize Fock Matrix
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
         implicit none
!
         class(ccs) :: wf
!
      end subroutine initialize_fock_matrix_ccs
!
!
      module subroutine construct_fock_ccs(wf)
!!
!!       Construct Fock
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
         implicit none
!
         class(ccs) :: wf
!
      end subroutine construct_fock_ccs
!
!
      module subroutine one_electron_t1_ccs(wf, h1 ,h1_T1)
!!
!!       One-electron T1
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         real(dp), dimension(wf%n_mo, wf%n_mo) :: h1
         real(dp), dimension(wf%n_mo, wf%n_mo) :: h1_T1
!
      end subroutine one_electron_t1_ccs
!
!
   end interface
!
!
   interface
!
!
!    -::- Ionized state submodule interface -::-
!    -------------------------------------------
!
      module subroutine ionized_state_driver_ccs(wf)
!!
!!       Ionized State Driver (CCS)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
         implicit none
!
         class(ccs) :: wf
!
      end subroutine ionized_state_driver_ccs
!
!
      module subroutine initialize_trial_vectors_valence_ionization_ccs(wf)
!!
!!       Initialize Trial Vectors for Valence Ionization
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, June 2017
!!
         implicit none
!
         class(ccs) :: wf
!
      end subroutine initialize_trial_vectors_valence_ionization_ccs
!
!
      module subroutine initialize_trial_vectors_core_ionization_ccs(wf)
!!
!!       Initialize Trial Vectors for Core Ionization
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, June 2017
!!
         implicit none
!
         class(ccs) :: wf
!
      end subroutine initialize_trial_vectors_core_ionization_ccs
!
!
      module subroutine precondition_residual_valence_ionization_ccs(wf, residual)
!!
!!       Precondition Residual Valence
!!       Written by Sarai D. Folkestad, Aug 2017
!!
         implicit none
!
         class(ccs) :: wf
         real(dp), dimension(wf%n_parameters ,1) :: residual
!
      end subroutine precondition_residual_valence_ionization_ccs
!
!
      module subroutine ionization_residual_projection_ccs(wf, residual)
!!
!!       Core Ionization Residual Projection (CCS)
!!       Written by Sarai D. Folkestad, Aug 2017
!!
         implicit none
!
         class(ccs) :: wf
         real(dp), dimension(wf%n_parameters, 1) :: residual
!
      end subroutine ionization_residual_projection_ccs
!
!
     module subroutine ionization_rho_a_i_projection_ccs(wf, rho_a_i)
!!
!!       Ionization rho_a_i Projection (CCS)
!!       Written by Sarai D. Folkestad, Aug 2017
!!
         implicit none
!
         class(ccs) :: wf
         real(dp), dimension(wf%n_v, wf%n_o) :: rho_a_i
!
      end subroutine ionization_rho_a_i_projection_ccs
!
!
      module subroutine precondition_residual_core_ionization_ccs(wf, residual)
!!
!!       Precondition Residual Valence
!!       Written by Sarai D. Folkestad, Aug 2017
!!
         implicit none
!
         class(ccs) :: wf
         real(dp), dimension(wf%n_parameters ,1) :: residual
!
      end subroutine precondition_residual_core_ionization_ccs
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
      module subroutine jacobian_ccs_transformation_ccs(wf, c_a_i)
!!
!!       Jacobian CCS transformation
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
         implicit none
!
         class(ccs) :: wf
         real(dp), dimension(wf%n_v, wf%n_o)   :: c_a_i
!
      end subroutine jacobian_ccs_transformation_ccs
!
!
      module subroutine jacobian_ccs_a1_ccs(wf,rho,c1)
!!
!!       Jacobian CCS A1
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
         implicit none
!
         class(ccs) :: wf
         real(dp), dimension(wf%n_o,wf%n_v) :: c1
         real(dp), dimension(wf%n_o,wf%n_v) :: rho
!
      end subroutine jacobian_ccs_a1_ccs
!
!
      module subroutine jacobian_ccs_b1_ccs(wf,rho,c1)
!!
!!       Jacobian CCS B1
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
         implicit none
!
         class(ccs) :: wf
         real(dp), dimension(wf%n_o,wf%n_v) :: c1
         real(dp), dimension(wf%n_o,wf%n_v) :: rho
!
      end subroutine jacobian_ccs_b1_ccs
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
      module subroutine jacobian_transpose_ccs_transformation_ccs(wf, b_a_i)
!!
!!       Jacobian transpose transformation (CCS)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         real(dp), dimension(wf%n_v, wf%n_o) :: b_a_i
!
      end subroutine jacobian_transpose_ccs_transformation_ccs
!
!
      module subroutine jacobian_transpose_ccs_a1_ccs(wf, sigma_a_i, b_a_i)
!!
!!       Jacobian transpose A1 (CCS)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         real(dp), dimension(wf%n_v, wf%n_o) :: sigma_a_i
         real(dp), dimension(wf%n_v, wf%n_o) :: b_a_i
!
      end subroutine jacobian_transpose_ccs_a1_ccs
!
!
      module subroutine jacobian_transpose_ccs_b1_ccs(wf, sigma_a_i, b_a_i)
!!
!!       Jacobian transpose B1 (CCS)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
         implicit none
!
         class(ccs) :: wf
!
         real(dp), dimension(wf%n_v, wf%n_o) :: sigma_a_i
         real(dp), dimension(wf%n_v, wf%n_o) :: b_a_i
!
      end subroutine jacobian_transpose_ccs_b1_ccs
!
!
   end interface
!
!
   interface
!
!     -::- CVS submodule interface -::-
!     ---------------------------------
!
      module subroutine cvs_rho_a_i_projection_ccs(wf, vec_a_i)
!!
!!       CVS rho_a_i Projection
!!       Written by Sarai D. Folkestad, Aug 2017
!!
         implicit none
!
         class(ccs) :: wf
         real(dp), dimension(wf%n_o,wf%n_v) :: vec_a_i
!
      end subroutine cvs_rho_a_i_projection_ccs
!
!
      module subroutine cvs_residual_projection_ccs(wf, residual)
!!
!!       CVS Residual Projection
!!       Written by Sarai D. Folkestad, Aug. 2017
!!
         implicit none
!
         class(ccs) :: wf
         real(dp), dimension(wf%n_parameters, 1) :: residual
!
      end subroutine cvs_residual_projection_ccs
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
   subroutine init_ccs(wf)
!!
!!    Initialize CCS object
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!    Performs the following tasks
!!
!!    - Sets HF orbital and energy information by reading from file
!!    - Transforms AO Cholesky vectors to MO basis and saves to file
!!    - Allocates the singles amplitudes and sets them to zero, and sets associated properties
!!    - Allocates the omega vector and sets it to zero
!!    - Initializes the Fock matrix and sets it to zero
!!
      implicit none
!
      class(ccs) :: wf
!
      integer(i15) :: unit_input = -1
!
!     Set model name
!
      wf%name = 'CCS'
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
!     Set implemented generic methods
!
      wf%implemented%ground_state         = .true.
      wf%implemented%excited_state        = .true.
      wf%implemented%ionized_state        = .true.
      wf%implemented%core_excited_state   = .true.
      wf%implemented%core_ionized_state   = .true.
      wf%implemented%multipliers          = .false.
!
!     Read calculation tasks from input file eT.inp
!
      call wf%calculation_reader(unit_input)
!
!     Close input file
!
      close(unit_input)
!
!     Figure out the size of the calculation folder, and update
!     the available disk space by subtracting it
!
      call wf%disk%subtract_folder_size
!
!     Read Hartree-Fock info from SIRIUS
!
      call wf%read_hf_info
!
!     Read Cholesky AO integrals and transform to MO basis
!
      call wf%read_transform_cholesky
!
!     Set amplitude attributes
!
      wf%n_t1am = (wf%n_o)*(wf%n_v)
      wf%n_parameters = wf%n_t1am ! The number of variables solved for by the solvers
!
!     Allocate Fock matrix and set to zero
!
      call wf%initialize_single_amplitudes ! t1am allocated, then set to zero
      call wf%initialize_fock_matrix
!
   end subroutine init_ccs
!
!
   subroutine drv_ccs(wf)
!!
!!    CCS Driver
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!    The driver for CCS is written so as to be inherited unaltered.
!!    It finds which calculations are requested by the user, and controls
!!    that the calculation can be done. If the method is implemented, it
!!    calls the driver for that particular calculation (e.g., ground state
!!    energy).
!!
      implicit none
!
      class(ccs) :: wf
!
      if (wf%tasks%ground_state) then
!
!        Ground state calculation requested
!
         if (wf%implemented%ground_state) then
!
            call wf%ground_state_driver
!
         else
!
            write(unit_output,'(t3,a,a)') &
               'Error: ground state solver not implemented for ',trim(wf%name)
            stop
!
         endif
      endif
!
      if (wf%tasks%excited_state) then
!
         if (wf%implemented%excited_state) then
!
            call wf%excited_state_driver
!
         else
!
            write(unit_output,'(t3,a,a)') &
               'Error: excited state solver not implemented for ',trim(wf%name)
            flush(unit_output)
            stop
!
         endif
!
      endif
!
      if (wf%tasks%core_excited_state) then
!
!        Excited state calculation requested
!
         if (wf%implemented%core_excited_state) then
!
            call wf%excited_state_driver
!
         else
!
            write(unit_output,'(t3,a,a)') &
               'Error: core excited state solver not implemented for ',trim(wf%name)
            flush(unit_output)
            stop
!
         endif
!
      endif
!
      if (wf%tasks%ionized_state) then
!
         if (wf%implemented%ionized_state) then
!
            call wf%ionized_state_driver
!
         else
!
            write(unit_output,'(t3,a,a)') &
               'Error: ionized state solver not implemented for ',trim(wf%name)
            flush(unit_output)
            stop
!
         endif
!
      endif
!
      if (wf%tasks%core_ionized_state) then
!
         if (wf%implemented%core_ionized_state) then
!
            call wf%ionized_state_driver
!
         else
!
            write(unit_output,'(t3,a,a)') &
               'Error: core ionized state solver not implemented for ',trim(wf%name)
            flush(unit_output)
            stop
!
         endif
!
      endif
!
      if (wf%tasks%multipliers) then
!
!        Multipliers calculation requested
!
         if (wf%implemented%multipliers) then
!
            call wf%response_driver
!
         else
!
            write(unit_output,'(t3,a,a)')'Error: Multipliers not implemented for ',trim(wf%name)
            stop
!
         endif
!
      endif
!
   end subroutine drv_ccs
!
!
   subroutine initialize_single_amplitudes_ccs(wf)
!!
!!    Initialize single amplitudes (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!    Allocates the singles amplitudes, sets them to zero.
!!
      implicit none
!
      class(ccs) :: wf
!
!     Allocate the singles amplitudes and set to zero
!     (which is also the value that solves the projected Scrödinger eq.)
!
      if (.not. allocated(wf%t1am)) then
!
         call wf%mem%alloc(wf%t1am, wf%n_v, wf%n_o)
         wf%t1am = zero
!
      else
!
         write(unit_output,'(t3,a)') 'Warning: attempted to allocate and zero already allocated t1am'
!
      endif
!
   end subroutine initialize_single_amplitudes_ccs
!
!
   subroutine initialize_omega_ccs(wf)
!!
!!    Initialize Omega (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!    Allocates and sets the projection vector to zero (which is
!!    also its correct value, by Brillouin)
!!
      implicit none
!
      class(ccs) :: wf
!
      if (.not. allocated(wf%omega1)) then
!
         call wf%mem%alloc(wf%omega1, wf%n_v, wf%n_o)
         wf%omega1 = zero
!
      else
!
         write(unit_output,'(t3,a)') 'Warning: attempted to allocate and zero already allocated omega1'
!
      endif
!
   end subroutine initialize_omega_ccs
!
!
   subroutine calc_energy_ccs(wf)
!!
!!    Calculate Energy (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
      implicit none
!
      class(ccs) :: wf
!
      wf%energy = wf%scf_energy
!
   end subroutine calc_energy_ccs
!
!
   subroutine construct_omega_ccs(wf)
!!
!!    Construct Omega (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
      implicit none
!
      class(ccs) :: wf
!
      wf%omega1 = zero ! Brillouin
!
   end subroutine construct_omega_ccs
!
!
   subroutine omega_ccs_a1_ccs(wf)
!!
!!    Omega D1
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, March 2017
!!
!!    Omega_ai^D1 = F_ai_T1
!!
      implicit none
!
      class(ccs) :: wf
!
!     Add F_a_i to omega
!
      call daxpy((wf%n_o)*(wf%n_v), one, wf%fock_ai, 1, wf%omega1, 1)
!
   end subroutine omega_ccs_a1_ccs
!
!
   subroutine construct_eta_ccs(wf, eta)
!!
!!    Construct Eta (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
!!    Note: the routine assumes that the Fock matrix
!!    has been constructed.
!!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_parameters, 1) :: eta ! eta_ai
!
      integer(i15) :: i = 0, a = 0, ai = 0
!
      eta = zero
!
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            ai = index_two(a, i, wf%n_v)
            eta(ai, 1) = two*(wf%fock_ia(i, a)) ! eta_ai = 2 F_ia
!
         enddo
      enddo
!
   end subroutine construct_eta_ccs
!
!
   subroutine save_amplitudes_ccs(wf)
!!
!!    Save amplitudes (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjøsntad, May 2017
!!
!!    Store the amplitudes to disk (T1AM)
!!
      implicit none
!
      class(ccs) :: wf
!
      integer(i15) :: unit_t1am = -1
!
!     Open amplitude file
!
      call generate_unit_identifier(unit_t1am)
      open(unit_t1am, file='t1am', status='unknown', form='unformatted')
      rewind(unit_t1am)
!
      write(unit_t1am) wf%t1am
!
!     Close amplitude file
!
      close(unit_t1am)
!
   end subroutine save_amplitudes_ccs
!
!
   subroutine read_single_amplitudes_ccs(wf)
!!
!!    Read single amplitudes (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjøsntad, May 2017
!!
!!    Reads the single amplitudes from disk. If the amplitudes are not allocated,
!!    the routine automatically allocates them.
!!
      implicit none
!
      class(ccs) :: wf
!
      integer(i15) :: unit_t1am = -1
!
      logical :: file_exists = .false.
!
!     Check to see whether file exists
!
      inquire(file='t1am',exist=file_exists)
!
      if (file_exists) then
!
!        Open amplitude files if they exist
!
         call generate_unit_identifier(unit_t1am)
!
         open(unit_t1am, file='t1am', status='unknown', form='unformatted')
!
         rewind(unit_t1am)
!
!        Allocate amplitudes if they aren't allocated
!
         if (.not. allocated(wf%t1am)) then
!
            call wf%mem%alloc(wf%t1am, wf%n_v, wf%n_o)
            wf%t1am = zero
!
         endif
!
!        Read from file
!
         read(unit_t1am) wf%t1am
!
!        Close file
!
         close(unit_t1am)
!
      else
!
         write(unit_output,'(t3,a)') 'Error: amplitude file t1am does not exist.'
         stop
!
      endif
!
   end subroutine read_single_amplitudes_ccs
!
!
   subroutine destruct_amplitudes_ccs(wf)
!!
!!    Destruct amplitudes (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!    Deallocates the amplitudes.
!!
      implicit none
!
      class(ccs) :: wf
!
      if (allocated(wf%t1am)) call wf%mem%dealloc(wf%t1am, wf%n_v, wf%n_o)
!
   end subroutine destruct_amplitudes_ccs
!
!
   subroutine destruct_single_amplitudes_ccs(wf)
!!
!!    Destruct single amplitudes (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!    Deallocates the single amplitudes.
!!
      implicit none
!
      class(ccs) :: wf
!
      if (allocated(wf%t1am)) call wf%mem%dealloc(wf%t1am, wf%n_v, wf%n_o)
!
   end subroutine destruct_single_amplitudes_ccs
!
!
   subroutine destruct_omega_ccs(wf)
!!
!!    Destruct omega (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!    Deallocates the projection vector.
!!
      implicit none
!
      class(ccs) :: wf
!
      if (allocated(wf%omega1)) call wf%mem%dealloc(wf%omega1, wf%n_v, wf%n_o)
!
   end subroutine destruct_omega_ccs
!
!
   subroutine jacobi_test_ccs(wf)
!
      implicit none
!
      class(ccs) :: wf
!
   end subroutine jacobi_test_ccs
!
!
   subroutine read_atom_info(n_nuclei, n_ao)
!!
!!    Read atom info,
!!    Written by Sarai Dery Folkestad, June 2017.
!!
!!    Reads atom info from DALTON generated file:
!!    Reads:
!!       - Number of nuclei
!!       - Number of AO's
!!
      implicit none
!
      integer(i15) :: n_nuclei,n_ao
!
      integer(i15) :: unit_center = 0
      integer(i15) :: ioerror     = 0
!
      call generate_unit_identifier(unit_center)
      open(unit=unit_center, file='center_info', status='unknown', form='unformatted', iostat=ioerror)
      if (ioerror .ne. 0) write(unit_output,*)'WARNING: Error while opening center_info'
      rewind(unit_center)
!
!     Read number of nuclei and aos
!
      read(unit_center) n_nuclei, n_ao
!
      close(unit_center)
!
   end subroutine read_atom_info
!
   subroutine read_center_info(n_nuclei, n_ao, n_ao_on_center, ao_center_info)
!!
!!    Read center info,
!!    Written by Sarai Dery Folkestad, June 2017.
!!
!!    Reads atom info from DALTON generated file:
!!    Reads:
!!       - Information of which ao's belong to which nuclei
!!
      implicit none
!
      integer(i15) :: n_nuclei
      integer(i15) :: n_ao
      integer, dimension(n_nuclei, 1)  :: n_ao_on_center
      integer, dimension(n_ao, 2)      :: ao_center_info
!
      integer(i15) :: offset = 0, nucleus = 0
      integer(i15) :: i = 0
      integer(i15) :: unit_center = 0, ioerror = 0
!
!     Read number of aos on each center
!
      call generate_unit_identifier(unit_center)
      open(unit=unit_center, file='center_info', status='unknown', form='unformatted', iostat=ioerror)
      if (ioerror .ne. 0) write(unit_output,*)'WARNING: Error while opening center_info'
      rewind(unit_center)
!
!     Empty read
!
      read(unit_center)
!
      read(unit_center, iostat=ioerror) n_ao_on_center
      if (ioerror .ne. 0) write(unit_output,*)'WARNING: Error while reading center_info'
!
      offset = 1
!
      do nucleus = 1, n_nuclei
!
         do i = 1, n_ao_on_center(nucleus, 1)
            ao_center_info(offset + i - 1,1)=nucleus
         enddo
!
!        Read which aos are on which centers
!
         read(unit_center, iostat=ioerror) (ao_center_info(offset + i - 1, 2), i = 1, n_ao_on_center(nucleus, 1))
         if (ioerror .ne. 0) write(unit_output,*)'WARNING: Error while reading center_info'
!
         offset = offset + n_ao_on_center(nucleus,1)
!
      enddo
!
      close(unit_center)
!
   end subroutine read_center_info
!
!
end module ccs_class
