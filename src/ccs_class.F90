module ccs_class
!
!!
!!              Coupled cluster singles (CCS) class module                                 
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
!
!  The ancestor class module (HF)
!
   use hf_class
!
   implicit none 
!
!  :::::::::::::::::::::::::::::::::::::
!  -::- Definition of the CCS class -::-
!  ::::::::::::::::::::::::::::::::::::: 
!
   type, extends(hf) :: ccs
!
!     Amplitude attributes
!
      integer(i15) :: n_t1am = 0                    ! Number of singles amplitudes
      real(dp), dimension(:,:), allocatable :: t1am ! Singles amplitude vector
!
      integer(i15) :: n_parameters = 0 ! Number of parameters in the wavefunction
!
!     Projection vector < mu | exp(-T) H exp(T) | R > (the omega vector)
! 
      real(dp), dimension(:,:), allocatable :: omega1 ! Singles vector 
!
!     The T1-transformed Fock matrix (in vir-occ block form)
!
      real(dp), dimension(:,:), allocatable :: fock_ij ! occ-occ block
      real(dp), dimension(:,:), allocatable :: fock_ia ! occ-vir block
      real(dp), dimension(:,:), allocatable :: fock_ai ! vir-occ block
      real(dp), dimension(:,:), allocatable :: fock_ab ! vir-vir block 
!
!     Variables that keep track of which response task is being performed 
!
      character(len=40) :: response_task 
      character(len=40) :: excited_state_task
      character(len=40) :: current_task = 'ground_state'
!
!     The excitation energies (omega_1, omega_2, ...)
!
      real(dp), dimension(:,:), allocatable :: excited_state_energies
!
   contains 
!
!     Initialization and driver routines
!
      procedure :: init => init_ccs
      procedure :: drv  => drv_ccs
!
!     Initialization routine for the amplitudes and omega 
!      
      procedure :: initialize_amplitudes => initialize_amplitudes_ccs
      procedure :: initialize_omega      => initialize_omega_ccs
!
!     Initialization routine for the Fock matrix, and a T1 Fock matrix constructor
!
      procedure, non_overridable :: construct_fock         => construct_fock_ccs
      procedure, non_overridable :: initialize_fock_matrix => initialize_fock_matrix_ccs
!
      procedure, non_overridable :: one_electron_t1 => one_electron_t1_ccs ! T1-transf. of h_pq
!
!     Routine to calculate the energy
!
      procedure :: calc_energy => calc_energy_ccs
!
!     get Cholesky routines to calculate the occ/vir-occ/vir blocks of the 
!     T1-transformed MO Cholesky vectors
!
      procedure, non_overridable :: get_cholesky_ij => get_cholesky_ij_ccs ! occ-occ
      procedure, non_overridable :: get_cholesky_ia => get_cholesky_ia_ccs ! occ-vir
      procedure, non_overridable :: get_cholesky_ai => get_cholesky_ai_ccs ! vir-occ
      procedure, non_overridable :: get_cholesky_ab => get_cholesky_ab_ccs ! vir-vir
!
!     Routine to construct projection vector (omega)
!
      procedure :: construct_omega => construct_omega_ccs
      procedure :: omega_ccs_a1    => omega_ccs_a1_ccs
!
!     Routine to construct right projection vector (eta)
!
      procedure :: construct_eta => construct_eta_ccs 
!
!     Ground state driver routine (and helpers)
!
!     Note: while this solver is uneccessary for CCS, where the solution is trivial, 
!     it is inherited mostly unaltered by descendants (CCSD, CC2, etc.).
!
      procedure :: ground_state_driver => ground_state_driver_ccs
!
!     Solver preparations and cleanup routines plus solver routine and its helpers
!
      procedure :: ground_state_preparations => ground_state_preparations_ccs
      procedure :: ground_state_solver       => ground_state_solver_ccs
      procedure :: ground_state_cleanup      => ground_state_cleanup_ccs
!
      procedure :: initialize_ground_state   => initialize_ground_state_ccs
      procedure :: destruct_ground_state     => destruct_ground_state_ccs
      procedure :: new_amplitudes            => new_amplitudes_ccs
      procedure :: calc_ampeqs               => calc_ampeqs_ccs
      procedure :: calc_ampeqs_norm          => calc_ampeqs_norm_ccs
      procedure :: calc_quasi_Newton_singles => calc_quasi_Newton_singles_ccs
!
      procedure, non_overridable :: diis     => diis_ccs
!
!     Routine to save and read the amplitudes (to/from disk)
!
      procedure :: save_amplitudes        => save_amplitudes_ccs
!
      procedure :: read_amplitudes        => read_amplitudes_ccs
      procedure :: read_single_amplitudes => read_single_amplitudes_ccs
!
!     Routines to destroy amplitudes and omega 
!
      procedure :: destruct_amplitudes => destruct_amplitudes_ccs
      procedure :: destruct_omega      => destruct_omega_ccs
!
!     Coupled cluster Jacobian transformation routine
!
      procedure :: jacobian_ccs_transformation => jacobian_ccs_transformation_ccs
!
      procedure, non_overridable :: jacobian_ccs_a1 => jacobian_ccs_a1_ccs 
      procedure, non_overridable :: jacobian_ccs_b1 => jacobian_ccs_b1_ccs
!
      procedure :: jacobi_test => jacobi_test_ccs ! A debug routine for A transformation 
!
!     Core-valence separation approximation routines
!
      procedure :: cvs_rho_a_i_projection             => cvs_rho_a_i_projection_ccs
      procedure :: cvs_residual_projection            => cvs_residual_projection_ccs
!
!     Coupled cluster Jacobian transpose transformation routine
!
      procedure :: jacobian_transpose_ccs_transformation => jacobian_transpose_ccs_transformation_ccs
!
      procedure, non_overridable :: jacobian_transpose_ccs_a1 => jacobian_transpose_ccs_a1_ccs
      procedure, non_overridable :: jacobian_transpose_ccs_b1 => jacobian_transpose_ccs_b1_ccs
!
!     Excited state driver & solver 
!
      procedure                  :: excited_state_preparations => excited_state_preparations_ccs
      procedure                  :: excited_state_driver       => excited_state_driver_ccs 
      procedure                  :: excited_state_cleanup      => excited_state_cleanup_ccs
      procedure, non_overridable :: excited_state_solver       => excited_state_solver_ccs
!
!     Helper routines for excited state solver 
!
      procedure :: initialize_excited_states     => initialize_excited_states_ccs
      procedure :: transform_trial_vectors       => transform_trial_vectors_ccs
      procedure :: calculate_orbital_differences => calculate_orbital_differences_ccs ! Must be overwritten for CCSD 
!
      procedure :: precondition_residual           => precondition_residual_ccs
      procedure :: precondition_residual_valence   => precondition_residual_valence_ccs
!
      procedure :: print_excited_state_info  => print_excited_state_info_ccs      
      procedure :: print_excitation_vector   => print_excitation_vector_ccs
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
!     Response driver & solver 
!
      procedure :: response_driver => response_driver_ccs
      procedure :: response_solver => response_solver_ccs
!
!     Helper routines for response solver 
!
      procedure :: initialize_response                   => initialize_response_ccs
      procedure :: solve_reduced_response_equation       => solve_reduced_response_equation_ccs
      procedure :: construct_reduced_matrix              => construct_reduced_matrix_ccs
      procedure :: construct_reduced_gradient            => construct_reduced_gradient_ccs
      procedure :: construct_next_response_trial_vectors => construct_next_response_trial_vectors_ccs
      procedure :: construct_gradient_vector             => construct_gradient_vector_ccs
!
!     Ionized state by super diffuse orbital
!
      procedure :: ionized_state_driver                           => ionized_state_driver_ccs
      procedure :: initialize_trial_vectors_core_ionization       => initialize_trial_vectors_core_ionization_ccs
      procedure :: initialize_trial_vectors_valence_ionization    => initialize_trial_vectors_valence_ionization_ccs
      procedure :: precondition_residual_valence_ionization       => precondition_residual_valence_ionization_ccs
      procedure :: ionization_residual_projection                 => ionization_residual_projection_ccs
      procedure :: ionization_rho_a_i_projection                  => ionization_rho_a_i_projection_ccs
      procedure :: precondition_residual_core_ionization          => precondition_residual_core_ionization_ccs
!
!     Integral routines (o: occupied index, v: virtual index)
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
!     Routine to store, read, and T1-transform read electronic repulsion integrals (g_abcd)
!
      procedure, non_overridable :: store_vv_vv_electronic_repulsion    => store_vv_vv_electronic_repulsion_ccs
!
      procedure, non_overridable :: read_vv_vv_electronic_repulsion     => read_vv_vv_electronic_repulsion_ccs
      procedure, non_overridable :: t1_transform_vv_vv                  => t1_transform_vv_vv_ccs
!
      procedure, non_overridable :: store_t1_vv_vv_electronic_repulsion => store_t1_vv_vv_electronic_repulsion_ccs
      procedure, non_overridable :: store_t1_vo_ov_electronic_repulsion => store_t1_vo_ov_electronic_repulsion_ccs
      procedure, non_overridable :: store_t1_vv_vo_electronic_repulsion => store_t1_vv_vo_electronic_repulsion_ccs
      procedure, non_overridable :: store_t1_vv_ov_electronic_repulsion => store_t1_vv_ov_electronic_repulsion_ccs
!
      procedure, non_overridable :: read_t1_vv_vo_electronic_repulsion  => read_t1_vv_vo_electronic_repulsion_ccs
      procedure, non_overridable :: read_t1_vv_vv_electronic_repulsion  => read_t1_vv_vv_electronic_repulsion_ccs
      procedure, non_overridable :: read_t1_vo_ov_electronic_repulsion  => read_t1_vo_ov_electronic_repulsion_ccs
      procedure, non_overridable :: read_t1_vv_ov_electronic_repulsion  => read_t1_vv_ov_electronic_repulsion_ccs
!
   end type ccs
!
!  ::::::::::::::::::::::::::::::::::::::::::::::::::::
!  -::- Interface to the submodule routines of CCS -::- 
!  ::::::::::::::::::::::::::::::::::::::::::::::::::::
!
   interface 
!
!
!     -::- Cholesky submodule interface -::-
!     :::::::::::::::::::::::::::::::::::::: 
!
      module subroutine get_cholesky_ij_ccs(wf, L_ij_J, i_first, i_last, j_first, j_last)
!!
!!       Get Cholesky IJ
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!       Reads and T1-transforms the IA Cholesky vectors:
!!     
!!          L_ij_J_T1 = L_ij_J + sum_a t_aj * L_ia_J
!!
!!       Memory required in routine:
!!
!!          2*n_J*(i_length)*n_v     -> for reading L_ia_J contribution and reordering
!!          i_length = i_last - i_first + 1
!!
!!       Optional arguments: i_first, i_last, j_first, j_last can be used in order to restrict indices
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
!!       Reads and T1-transforms IA Cholesky vectors
!!
!!          L_ia_J_T1 = L_ia_J (only reading necessary)
!!
!!       Memory required in routine:
!!
!!          No additional memory
!!
!!       Optional arguments: i_first, i_last, a_first, a_last can be used in order to restrict indices
!!
         implicit none 
!
         class(ccs)               :: wf
         integer(i15), optional   :: i_first, a_first   ! First index (can differ from 1 when batching or for mlcc)
         integer(i15), optional   :: i_last, a_last    ! Last index (can differ from n_v/n_o when batching or for mlcc)
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
!!       Read and T1-transform Cholesky AI vectors:
!!       
!!          L_ai_J_T1 = L_ia_J - sum_j  t_aj*L_ji_J 
!!                             + sum_b  t_bi*L_ab_J
!!                             - sum_bj t_aj*t_bi*L_jb_J
!!
!!       Allocations in routine:
!!
!!         (1) n_J*(i_length)*(a_length) + 2*n_J*(a_length)*batch_length  ->  for L_ab_J contribution (batches of b)
!!         (2) n_J*(i_length)*n_v + 2*n_J*n_o*(i_length)                  ->  for L_ij_J contribution
!!         (3) 2*n_J*n_o*n_v                                              ->  for L_jb_J contribution
!!
!!         i_length = i_last - i_first + 1          
!!         a_length = a_last - a_first + 1          
!!
!!         (1) determines memory requirement. 
!!
!!       Optional arguments: i_first, i_last, a_first, a_last can be used in order to restrict indices
!!
         implicit none 
!
         class(ccs)               :: wf
         integer(i15), optional   :: i_first, a_first     ! First index (can differ from 1 when batching or for mlcc)
         integer(i15), optional   :: i_last, a_last      ! Last index (can differ from n_o/n_v when batching or for mlcc)
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
!!       Reads and T1-transforms the IA Cholesky vectors:
!!
!!          L_ab_J_T1 = L_ab_J - sum_i t_ai*L_ib_J
!!
!!
!!       Required memory: 
!!
!!          n_J*b_length*a_length       ->   For reordering of L_ab_J / L_ba_J
!!          2*b_length*n_o*n_J          ->   For L_ib_J contribution
!!
!!         a_length = a_last - a_first + 1          
!!         b_length = b_last - b_first + 1 
!!

         implicit none
!
         class(ccs)               :: wf
         integer(i15), intent(in) :: a_first, b_first   ! First index (can differ from 1 when batching or for mlcc)
         integer(i15), intent(in) :: a_last, b_last    ! Last index (can differ from n_v when batching or for mlcc)
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
!    -::- Fock submodule interface -::-
!    ::::::::::::::::::::::::::::::::::
!
!
      module subroutine initialize_fock_matrix_ccs(wf)
!!  
!!       Initialize Fock matrix
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!       Allocates and sets Fock matrix blocks (ij, ia, ai, ab) to zero
!!       before calling the Fock matrix constructor.
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
!!       Constructs the T1-transformed Fock matrix blocks (occ/vir-occ/vir),
!!       and saves the result in the class variables fock_pq.  
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
!!       T1-transforms the one-electron MO integrals h_pq
!!
!!           h_p_q_T1 = sum_st x_p_s * y_q_t * h_s_t,
!!
!!       where
!!
!!           x = I - t1,
!!           y = I - t1^T.
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
!     -::- Ground state submodule interface -::-
!     ::::::::::::::::::::::::::::::::::::::::::
!
      module subroutine ground_state_driver_ccs(wf)
!!
!!       Ground state driver (CCS)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
!!       Directs the solution of the ground state problem for CCS. The
!!       routine is written so as to be inherited unaltered in the CC hierarchy. 
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
!!       A routine for preparation tasks (if any). Can be overwritten
!!       in descendants if other preparations prove necessary.    
!!
         class(ccs) :: wf 
!
      end subroutine ground_state_preparations_ccs
!
!
      module subroutine ground_state_cleanup_ccs(wf)
!!
!!       Ground State Cleanup (CCS)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
!!       A routine for cleanup tasks (if any). Can be overwritten
!!       in descendants if other cleanups prove necessary.    
!!
         class(ccs) :: wf 
!
      end subroutine ground_state_cleanup_ccs
!
!
      module subroutine ground_state_solver_ccs(wf)
!!
!!       Ground State Solver 
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!       Directs the solution of the ground state amplitude equations
!!       using a DIIS algorithm. The problem the routine solves is 
!!
!!          X_mu(t) = 0, where t = { t_mu }_mu 
!!
!!       For standard coupled cluster theories, the vector X is the
!!       projection vector (omega).
!!
         implicit none 
!
         class(ccs) :: wf 
!
      end subroutine ground_state_solver_ccs
!
!
      module subroutine calc_ampeqs_ccs(wf)
!!
!!       Calculate Amplitude Equations (CCS)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!       Constructs the amplitude equations vector (the projection vector 
!!       in CCS) for the amplitudes of the current iteration of the ground state
!!       solver. It also calculates the norm of the amplitude equations, which 
!!       is zero when the equations are exactly solved.
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
      module subroutine new_amplitudes_ccs(wf)
!!
!!       New Amplitudes (CCS)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!       Directs the calculation of the quasi-Newton estimate Δ t_i, 
!!       and t_i + Δ t_i, and calls the DIIS routine to save & get 
!!       the amplitudes for the next iteration.
!!
         implicit none 
!
         class(ccs) :: wf 
!
      end subroutine new_amplitudes_ccs
!
!
      module subroutine calc_quasi_Newton_singles_ccs(wf,dt)
!!
!!       Calculate quasi-Newton estimate (CCS)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!       Calculates the quasi-Newton estimate Δ t_i (singles part)
!!       and places the contribution in the dt vector (of length n_parameters,
!!       with singles first, then doubles, etc. if inherited)
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
!!       DIIS routine (CCS)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad
!!
!!       The next amplitudes are 
!!
!!          t_n+1 = sum_k w_k (t_k + dt_k), 
!! 
!!       where the weights w_k in front of the quasi-Newton estimate dt_k
!!       are determined so as to minimize 
!!
!!          f(w_k) = sum_k w_k dt_k, 
!!
!!       with the constraint that g(w_k) = sum_k w_k - 1 = 0.
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
      module subroutine initialize_ground_state_ccs(wf)
!!
!!       Initialize Ground State (CCS)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!       Initializes the amplitudes and the projection vector. This routine 
!!       can be inherited unaltered by standard CC methods.
!!
         implicit none 
!
         class(ccs) :: wf
!
      end subroutine initialize_ground_state_ccs
!
!
      module subroutine destruct_ground_state_ccs(wf)
!!
!!       Destruct Ground State (CCS)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!       Deallocates the amplitudes and the projection vector. This routine 
!!       can be inherited unaltered by standard CC methods.
!!
         implicit none
!
         class(ccs) :: wf
!
      end subroutine destruct_ground_state_ccs
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
      module subroutine initialize_trial_vectors_ccs(wf)
!!
!!       Initialize trial vectors (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!       Initializes start trial vectors for the calculation of 
!!       singlet excited states and writes them to file 'trial_vecs'.
!!
!!       n start vectors are constructed by finding the n lowest orbital differences,      
!!       where n = n_singlet_states. Vector i has a 1.0D0 at the element corresponding to the i'th lowest
!!       orbital difference and 0.0d0 everywhere else
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
         integer(i15), dimension(wf%tasks%n_singlet_states,1), intent(inout) :: index_list
! 
      end subroutine find_start_trial_indices_ccs
!
!
      module subroutine transform_trial_vectors_ccs(wf, first_trial, last_trial)
!!
!!       Transform trial vectors (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!       Each trial vector in first_trial to last_trial is read from file and
!!       transformed before the transformed vector is written to file.
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
!!       Directs the solution of the excited state problem for CCS. The
!!       routine is inherited is to be inherited unaltered in the CC hierarchy. 
!!
!!       Note: it is only necessary to alter this routine if the excited states are 
!!       solved for by a different algorithm (such as in similarity constrained CC, 
!!       where the excited states and ground state are determined simultaneously).
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
!!       A routine for preparation tasks (if any). Can be overwritten
!!       in descendants if other preparations prove necessary.    
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
!!       A routine for cleanup tasks (if any). Can be overwritten
!!       in descendants if other cleanups prove necessary.    
!!
         class(ccs) :: wf 
!
      end subroutine excited_state_cleanup_ccs
!
!
      module subroutine excited_state_solver_ccs(wf)
!!
!!       Excited State Solver
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!       Directs the solution of the excited states using a Davidson algorithm.
!!       The routine aims to find the right eigenvectors of the Jacobian matrix
!!
!!          AX = eX,  
!!
!!       and the eigenvalues which corresponds to the excitation energies.
!!
!!       The problem is solved in reduced space. To find n roots, n start trial vectors {c_i}_i=1, 
!!       n are generated according to the lowest orbital differences. Then a reduced space Jacobian 
!!       is constructed,
!! 
!!          A_red_ij = c_i^T * A c_j,
!!
!!       and the eigenvalues e and eigenvectors x of this matrix are found.
!!       The full space vectors {X_j}_j=1,n are then given by
!!
!!          X_j = sum_i x_j_i*c_i, 
!! 
!!       and the j'th residual vector is given by
!!
!!          R_j = (A*X_j - e*X_j)/|X_j|.
!!
!!       If the norm of this residual is sufficiently small (and the excitation energies 
!!       are converged within a given threshold), convergence is reached. If not, new trial 
!!       vectors will be generated by orthogonalizing the residual vector against the previous 
!!       trial vectors and then normalizing them, thereby expanding the dimension 
!!       of the reduced space for the next iteration.
!!
!!       The linear system (equivalently, the residual) is preconditioned with a diagonal 
!!       matrix with elements equal to the inverse orbital differences.
!!   
         implicit none
!  
         class(ccs) :: wf
!
      end subroutine excited_state_solver_ccs
!
!
      module subroutine solve_reduced_eigenvalue_equation_ccs(wf, eigenvalues_Re, eigenvalues_Im, &
                                                               solution_vectors_reduced, reduced_dim, n_new_trials)
!!
!!       Solve Reduced Eigenvalue Equation
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!       Constructs the reduced A matrix, solves its eigenvalue equation,
!!       and returns its first n eigenvalues and eigenvectors (reduced space
!!       solution vectors).
!!
         implicit none
!
         class(ccs) :: wf
!
         integer(i15) :: reduced_dim, n_new_trials
!
         real(dp), dimension(wf%tasks%n_singlet_states,1) :: eigenvalues_Re
         real(dp), dimension(wf%tasks%n_singlet_states,1) :: eigenvalues_Im
! 
         real(dp), dimension(reduced_dim, wf%tasks%n_singlet_states) :: solution_vectors_reduced
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
!!       Constructs the next eigenvectors by constructing the residual vectors
!!    
!!          R_j = (A*X_j - e*X_j)/|X_j|,
!!
!!       orthogonalizing them against the other trial vectors.
!!
!!       Residual vectors are preconditioned before orthogonalization.
!!       This is done by dividing by the orbital differences.
!!    
!!       If norm of orthogonal vector is very small 
!!       (i.e. high degree of linear dependence on previous trial vectors)
!!       it is scrapped. If norm sufficiently large, vector is normalized and
!!       stored in trial_vec file, to be used in the next iteration.
!!
!!       Routine also constructs full space solution vectors and stores them
!!       in file solution_vectors 
!! 
         implicit none
!
         class(ccs) :: wf
!
         real(dp), dimension(wf%tasks%n_singlet_states,1) :: eigenvalues_Re
         real(dp), dimension(wf%tasks%n_singlet_states,1) :: eigenvalues_Im
!
         integer(i15) :: reduced_dim
         integer(i15) :: n_new_trials
         real(dp), dimension(reduced_dim, wf%tasks%n_singlet_states) :: solution_vectors_reduced
!
      end subroutine construct_next_trial_vectors_ccs
!
!
      module subroutine initialize_trial_vectors_valence_ccs(wf)
!!
!!       Initialize trial vectors valence
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad
!!
!!       Initializes start trial vectors for the calculation of 
!!       singlet excited states and writes them to file 'trial_vecs'.
!!
!!       n start vectors are constructed by finding the n lowest orbital differences,      
!!       where n = n_singlet_states. Vector i has a 1.0D0 at the element corresponding to the i'th lowest
!!       orbital difference and 0.0d0 everywhere else
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
!!       Written by Sarai D. Folkestad, Aug. 2017
!!
!!       Finds correct core MO, and selects the n_singlet_state lowest 
!!       orbital differences where one of the occupied indices corresponds to the 
!!       core MO
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
!!    Trial vectors from old solutions,
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, June 2017
!!
!!    Restart: Use old solutions as trial vectors
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
!!       Find indices for lowest orbital differences core excitation
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad
!!
         implicit none
!
         class(ccs) :: wf
         integer(i15), dimension(wf%tasks%n_singlet_states,1), intent(inout) :: index_list
!
      end subroutine find_start_trial_indices_core_ccs
!
!
      module subroutine find_core_mo_ccs(wf)
!!
!!       Find which mo are core mos
!!       Written by Sarai D. Folkestad, Aug. 2017
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
!!       Calculate and return orbital differences
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad May 2017
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
!!       Initialize excited states
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, June 2017
!!
         implicit none 
!    
         class(ccs) :: wf
!
      end subroutine initialize_excited_states_ccs
!
!
      module subroutine cvs_residual_projection_ccs(wf, residual)
!!
!!       Residual projection for CVS
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
      module subroutine precondition_residual_ccs(wf, residual)
!!
!!       Precondition residual
!!       Written by Sarai D. Folkestad, Aug. 2017
!!
!!       Calls precondition_residual_valence for normal excited state calculation
!!       Calls precondition_residual_core for cvs calculation
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
!!       Precondition residual valence
!!       Written by Sarai D. Folkestad, Aug. 2017
!!
!!       Divide elements of residual by orbital difference       
!!
         implicit none
!
         class(ccs) :: wf
         real(dp), dimension(wf%n_parameters ,1) :: residual
!   
!
      end subroutine precondition_residual_valence_ccs
!
!
!
      module subroutine precondition_residual_core_ccs(wf, residual)
!!
!!       Precondition residual core
!!       Written by Sarai D. Folkestad, Aug. 2017
!!
!!       Project out elements not corresponding to the core excitation
!!       Divide elements of residual by orbital difference
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
!!
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
!!
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
   end interface
!
!
   interface
!
!    -::- Ionized state submodule interface -::-
!    :::::::::::::::::::::::::::::::::::::::::::
!

      module subroutine ionized_state_driver_ccs(wf)
!!
!!       Ionized state driver (CCS)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
!!       Directs the solution of the eionized state problem for CCS. The
!!       routine is inherited is to be inherited unaltered in the CC hierarchy. 
!!
!!       Note: it is only necessary to alter this routine if the ionized states are 
!!       solved for by a different algorithm (such as in similarity constrained CC, 
!!       where the excited states and ground state are determined simultaneously).
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
!!
!!       Initialize trial vectors for valence ionized state
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad
!!
!!       Initializes start trial vectors for the calculation of 
!!       singlet excited states and writes them to file 'trial_vecs'.
!!
!!       n start vectors are constructed by finding the n lowest orbital differences,      
!!       where n = n_singlet_states. Vector i has a 1.0D0 at the element corresponding to the i'th lowest
!!       orbital difference and 0.0d0 everywhere else
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
!!       Initialize trial vectors for core ionized state
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad
!!
!!       Initializes start trial vectors for the calculation of 
!!       singlet excited states and writes them to file 'trial_vecs'.
!!
!!       n start vectors are constructed by finding the n lowest orbital differences,      
!!       where n = n_singlet_states. Vector i has a 1.0D0 at the element corresponding to the i'th lowest
!!       orbital difference and 0.0d0 everywhere else
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
!!       Precondition residual valence
!!       Written by Sarai D. Folkestad, Aug. 2017
!!
!!       Projects out contamination from regular excitations
!!       Divide elements of residual by orbital difference          
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
!!       Residual projection for core ionization (CCS)
!!       Written by Sarai D. Folkestad, Aug. 2017
!!
!!       Projects out contaminations from regular excitations
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
!!       Residual projection for CVS
!!       Written by Sarai D. Folkestad, Aug. 2017
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
!     
!        Precondition residual valence
!        Written by Sarai D. Folkestad, Aug. 2017
!     
!        Divide elements of residual by orbital difference       
!     
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
!    -::- Jacobian submodule interface -::-
!    ::::::::::::::::::::::::::::::::::::::
!
      module subroutine core_jacobian_ccs_transformation_ccs(wf, c_a_i)
!!
!!       Jacobian transformation for CVS calculation
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
         implicit none
!
         class(ccs) :: wf 
         real(dp), dimension(wf%n_v, wf%n_o)   :: c_a_i    
!
      end subroutine core_jacobian_ccs_transformation_ccs
!
      module subroutine jacobian_ccs_a1_ccs(wf,rho,c1)
!!
!!       A1 contribution to right transform of Jacobian
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad
!!
!!       Calculates the A1 term of the right transform of the
!!       Jacobian,
!!
!!       A1: sum_b F_ab*c_bi + sum_j F_ji*c_aj
!!
!!       and adds it to the rho vector.
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
!!       B1 contribution to right transform of Jacobian
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad
!!
!!       Calculates the B1 term of the right transform of the
!!       Jacobian,
!!
!!       B1: sum_bj L_aijb*c_bj
!!
!!       and adds it to the rho vector.
!!
         implicit none
!
         class(ccs) :: wf
         real(dp), dimension(wf%n_o,wf%n_v) :: c1
         real(dp), dimension(wf%n_o,wf%n_v) :: rho                
!
      end subroutine jacobian_ccs_b1_ccs
!
!
      module subroutine cvs_rho_a_i_projection_ccs(wf, vec_a_i)
!!
!!       Rho projection for CVS
!!       Written by Sarai D. Folkestad, Aug. 2017
!!
         implicit none
!
         class(ccs) :: wf
         real(dp), dimension(wf%n_o,wf%n_v) :: vec_a_i
!
      end subroutine cvs_rho_a_i_projection_ccs
!
!
   end interface 
!
!
   interface 
!
!     -::- Jacobian submodule interface -::-
!     ::::::::::::::::::::::::::::::::::::::
!
      module subroutine jacobian_ccs_transformation_ccs(wf, c_a_i)
!!
!!       Jacobian CCS transformation
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!       Directs the transformation by the CCSD Jacobi matrix,
!!
!!          A_mu,nu = < mu | exp(-T) [H, tau_nu] exp(T) | nu >. 
!!
!!       In particular,
!!
!!          rho_mu = (A c)_mu = sum_ck A_mu,ck c_ck.
!! 
!!       On exit, c is overwritten by rho. 
!!
         implicit none
!
         class(ccs) :: wf 
         real(dp), dimension(wf%n_v, wf%n_o)   :: c_a_i       
!
      end subroutine jacobian_ccs_transformation_ccs
!
!
   end interface 
!
!
   interface 
!
!     -::- Jacobian transpose submodule interface -::-
!     ::::::::::::::::::::::::::::::::::::::::::::::::
!
      module subroutine jacobian_transpose_ccs_transformation_ccs(wf, b_a_i)
!!
!!       Jacobian transpose transformation (CCS)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
!!       Calculates the transpose Jacobian transformation, i.e., the transformation 
!!       by the transpose of the Jacobian matrix
!!
!!          A_mu,nu = < mu | exp(-T) [H, tau_nu] exp(T) | R >.
!!
!!       The transformation is performed as sigma^T = b^T A, where b is the vector
!!       sent to the routine. On exit, the vector b is equal to sigma (the transformed
!!       vector).
!
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
!!       Calculates the A1 term,
!!
!!          sum_c b_ci F_ca - sum_k b_ak F_ik,
!!
!!       and adds it to the sigma-vector (b^T -> sigma^T = b^T A).
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
!!       Calculates the B1 term,
!!
!!          sum_ck L_ckia b_ck
!!
!!       and adds it to the sigma-vector (b^T -> sigma^T = b^T A).
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
!     -::- Response submodule interface -::-
!     ::::::::::::::::::::::::::::::::::::::
!
      module subroutine response_driver_ccs(wf)
!!
!!       Response driver (CCS)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
!!       Directs the solution of molculear properties for CCS. The
!!       routine is inherited is to be inherited unaltered in the CC hierarchy. 
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
!!       Solves the response equation
!!
!!          A X = F,
!!
!!       where F is the gradient vector and A is either the Jacobian (A) or 
!!       the transposed Jacobian (A^T) matrix. 
!!
!!       The equation is solved by constructing a set of trial vectors, c_i, and 
!!       solving the projected/reduced equations associated with the reduced A and
!!       reduced F vectors:
!!
!!          A_ij = c_i^T A c_j,    F_i = c_i^T F 
!!
!!       The space of trial vectors is expanded according to a Davidson algorithm.
!!
         implicit none 
!
         class(ccs) :: wf
!
      end subroutine response_solver_ccs
!
!
      module subroutine initialize_response_ccs(wf)
!!
!!       Initialize response (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, June 2017
!!
!!       Performs two tasks:
!!
!!       1. Initializes start trial vector for response solver. We use the 
!!       the diagonal approximation D of A (and A^T) to form the trial vector:
!!
!!          A X = F => X ~ D^-1 F 
!!
!!       The diagonal approximation of A consists of the orbital energy differences.
!!
!!       2. Constructs the vector F and saves it to file.
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
!!       Constructs the reduced A matrix and the reduced F vector,
!!       and solves the (reduced space) linear equation A X = F for X. 
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
!!       Constructs A_ij = c_i^T A c_j by reading from file and constructing the missing elements.
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
!!       Constructs F_i = c_i^T F by reading from file and constructing the missing elements.
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
!!       Constructs the gradient vector, given the current value of "response_task",
!!       and stores the vector on disk for use by the solver. 
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
!!       Constructs the next trial vector by constructing the residual vector
!!    
!!          R = (A*X - F)/|X|,
!!
!!       and orthogonalizing it against the previous trial vectors.
!!
!!       Residual vectors are preconditioned before orthogonalization.
!!       This is done by dividing by the orbital differences.
!!    
!!       If the norm of orthogonal vector is very small 
!!       (i.e. high degree of linear dependence on previous trial vectors)
!!       it is scrapped. If norm sufficiently large, the vector is normalized and
!!       stored in trial_vec file, to be used in the next iteration.
!!
!!       The routine also constructs full space solution vectors and stores them
!!       in file solution_vectors 
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
!     -::- Integral submodule interface -::-
!     ::::::::::::::::::::::::::::::::::::::
!
      module subroutine store_vv_vv_electronic_repulsion_ccs(wf)
!
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
!!       Tests whether it is possible to store t1-transformed vir-vir-vir-vir integrals and,
!!       if possible, writes the integrals to disk 
!!
         implicit none 
!
         class(ccs) :: wf 
!
      end subroutine store_t1_vv_vv_electronic_repulsion_ccs
!
      module subroutine store_t1_vo_ov_electronic_repulsion_ccs(wf)
!!
!!       Store t1 voov Electronic Repulsion Integrals 
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!! 
!!       Tests whether it is possible to store t1-transformed vir-occ-occ-vir integrals and,
!!       if possible, writes the integrals to disk 
!! 
         implicit none 
!  
         class(ccs) :: wf 
!
      end subroutine store_t1_vo_ov_electronic_repulsion_ccs
!
      module subroutine store_t1_vv_ov_electronic_repulsion_ccs(wf)
!!
!!       Store t1 voov Electronic Repulsion Integrals 
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!! 
!!       Tests whether it is possible to store t1-transformed vir-vir-occ-vir integrals and,
!!       if possible, writes the integrals to disk 
!! 
         implicit none 
!  
         class(ccs) :: wf 
!
      end subroutine store_t1_vv_ov_electronic_repulsion_ccs
!
      module subroutine store_t1_vv_vo_electronic_repulsion_ccs(wf)
!!
!!       Store t1 vvvo Electronic Repulsion Integrals 
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
!!       Tests whether it is possible to store t1-transformed vir-vir-vir-occ integrals and,
!!       if possible, writes the integrals to disk 
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
!
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
!
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
      module subroutine read_t1_vo_ov_electronic_repulsion_ccs(wf, x_vo_ov,    & 
                                       index1_first, index1_last, &
                                       index2_first, index2_last, &
                                       index3_first, index3_last, &
                                       index4_first, index4_last)
!
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
!!       Reads the T1-transformed vir-vir-occ-vir integrals from file.
!!
!!       The integrals are stored on file as (a, bci) = (a, :), where a 
!!       is the record number and : denotes all the bci elements.
!!
!!       The recommended use is therefore to batch over the a index,
!!       as this will involve the minimum amount of wasteful read statements
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
!!       Reads the T1-transformed vir-vir-occ-vir integrals from file.
!!
!!       The integrals are stored on file as (a, bci) = (a, :), where a 
!!       is the record number and : denotes all the bci elements.
!!
!!       The recommended use is therefore to batch over the a index,
!!       as this will involve the minimum amount of wasteful read statements
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
!
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
!  
!  
!  
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
!  
!  
!  
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
!  
!  
!  
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
!  
!  
!  
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
!  
!  
!  
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
!  
!  
!  
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
!  
!  
!  
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
!  
!  
!  
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
!  
!  
!  
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
!  
!  
!  
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
!  
!  
!  
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
!  
!  
!  
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
!  
!  
!  
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
!  
!  
!  
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
!  
!  
!  
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
!!
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
!!
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
!!
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
!!
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
!!
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
!!
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
!!
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
!!
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
!!  
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
!  
!  
!  
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
!  
!  
!  
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
!  
!  
!  
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
!  
!  
!  
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
!  
!  
!  
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
!  
!  
!  
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
!  
!  
!  
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
   module subroutine t1_transform_vv_vv_ccs(wf, g_vv_vv,                & 
                                             index1_first, index1_last, &
                                             index2_first, index2_last, &
                                             index3_first, index3_last, &
                                             index4_first, index4_last)
!!
!!       T1 transformation of g_vv_vv integrals (CCS)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, Oct 2017. 
!!
!!       g_ab_cd_T1 = g_ab_cd + sum_(J) sum_(i) t_a_i * L_ib_J * L_cd_J
!!                            + sum_(J) sum_(k) t_c_k * L_kd_J * L_ab_J
!!                            + sum_(J) sum_(ki) t_a_i * t_c_k * L_kd_J * L_ib_J
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
   end interface
!
!
contains
!
!  ::::::::::::::::::::::::::::::::::::::::::::
!  -::- Initialization and driver routines -::- 
!  ::::::::::::::::::::::::::::::::::::::::::::
!
   subroutine init_ccs(wf)
!!
!!    Initialize CCS object
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!    Performs the following tasks
!!
!!    1. Sets HF orbital and energy information by reading from file
!!    2. Transforms AO Cholesky vectors to MO basis and saves to file 
!!    3. Allocates the singles amplitudes and sets them to zero, and sets associated properties 
!!    4. Allocates the omega vector and sets it to zero
!!    5. Initializes the Fock matrix and sets it to zero 
!!
      implicit none 
!
      class(ccs) :: wf
!
!     Set model name 
!
      wf%name = 'CCS'
!
!     Set implemented generic methods
!
      wf%implemented%ground_state = .true.
      wf%implemented%excited_state = .true.
      wf%implemented%ionized_state = .true.
      wf%implemented%core_excited_state = .true.
      wf%implemented%core_ionized_state = .true.
      wf%implemented%properties = .false.
!
!     Read Hartree-Fock info from SIRIUS
!
      call wf%read_hf_info
!
!     Read Cholesky AO integrals and transform to MO basis
!
      call wf%read_transform_cholesky
!
!     Test for the possibility of storing vir-vir-vir-vir
!     electronic repulsion integrals (g_abcd), storing the
!     integrals if possible
!
     ! call wf%store_vvvv_electronic_repulsion_integrals Should this occur in ccs? I don't think so
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
!     Initialize the projection vector 
!
      call wf%initialize_omega
!
!     Allocate Fock matrix and set to zero
!
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
            wf%current_task = 'ground_state'
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
            wf%excited_state_task = 'right_valence'
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
            wf%excited_state_task = 'right_core'
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
            wf%excited_state_task = 'right_valence'
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
            wf%excited_state_task = 'right_core'
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
      if (wf%tasks%properties) then
!
!        Properties calculation requested
!
         if (wf%implemented%properties) then 
!
            wf%current_task = 'response'
            wf%response_task = 'multipliers'
            call wf%response_driver
!
         else
!
            write(unit_output,'(t3,a,a)') &
               'Error: properties not implemented for ',trim(wf%name)
            stop
!
         endif
!
      endif
!
   end subroutine drv_ccs
!
!  :::::::::::::::::::::::::::::::::::::::::
!  -::- Class subroutines and functions -::- 
!  :::::::::::::::::::::::::::::::::::::::::
!
   subroutine initialize_amplitudes_ccs(wf)
!!
!!    Initialize Amplitudes (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!    Allocates the singles amplitudes, sets them to zero, and calculates
!!    the number of singles amplitudes.
!!
      implicit none 
!
      class(ccs) :: wf
!
!     Calculate the number of singles amplitudes
!
      wf%n_t1am = (wf%n_o)*(wf%n_v) 
!
!     Allocate the singles amplitudes and set to zero
!     (which is also the value that solves the projected Scrödinger eq.)
!
      if (.not. allocated(wf%t1am)) call allocator(wf%t1am, wf%n_v, wf%n_o)
      wf%t1am = zero
!
   end subroutine initialize_amplitudes_ccs
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
      if (.not. allocated(wf%omega1)) call allocator(wf%omega1, wf%n_v, wf%n_o)
      wf%omega1 = zero
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
!!    Save Amplitudes (CCS)
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
   subroutine read_amplitudes_ccs(wf)
!!
!!    Read Amplitudes (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjøsntad, May 2017
!!
!!    Reads the amplitudes from disk (T1AM)
!!
      implicit none 
!
      class(ccs) :: wf
!
      call wf%read_single_amplitudes
!
   end subroutine read_amplitudes_ccs
!
   subroutine read_single_amplitudes_ccs(wf)
!!
!!    Read Amplitudes (CCSD)
!!    Written by Sarai D. Folkestad and Eirik F. Kjøsntad, May 2017
!!
!!    Reads the amplitudes from disk (T1AM, T2AM)
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
!        Read from file & close
!
         wf%t1am = zero
!
         read(unit_t1am) wf%t1am 
!  
         close(unit_t1am)
!
      else
!
         write(unit_output,'(t3,a)') 'Error: amplitude files do not exist.'
         stop
!
      endif
!
   end subroutine read_single_amplitudes_ccs
!
!
   subroutine destruct_amplitudes_ccs(wf)
!!
!!    Destruct Amplitudes (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjøsntad, May 2017
!!
!!    Deallocates the amplitudes.
!!
      implicit none
!
      class(ccs) :: wf
!
      if (allocated(wf%t1am)) then
         call deallocator(wf%t1am, wf%n_v, wf%n_o)
      endif
!
   end subroutine destruct_amplitudes_ccs
!
!
   subroutine destruct_omega_ccs(wf)
!!
!!    Destruct Omega (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjøsntad, May 2017
!!
!!    Deallocates the projection vector.
!!
      implicit none
!
      class(ccs) :: wf
!
      if (allocated(wf%omega1)) then
         call deallocator(wf%omega1, wf%n_v, wf%n_o)
      endif
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
