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
      character(len=40)                     :: excited_state_task     = 'right_eigenvectors' ! The task being performed at the moment 
      real(dp), dimension(:,:), allocatable :: excited_state_energies
!
   contains 
!
!     Initialization and driver routines
!
      procedure :: init => init_ccs
      procedure :: drv  => drv_ccs
!
!     Initialization routine for the (singles) amplitudes
!      
      procedure :: initialize_amplitudes => initialize_amplitudes_ccs
!
!     Routine to initialize omega (allocate and set to zero)
!
      procedure :: initialize_omega => initialize_omega_ccs
!
!     Initialization routine for the Fock matrix, and a Fock matrix constructor
!     (for the given T1 amplitudes)
!
      procedure                  :: initialize_fock_matrix => initialize_fock_matrix_ccs
      procedure, non_overridable :: construct_fock         => construct_fock_ccs
!
      procedure, non_overridable :: one_electron_t1        => one_electron_t1_ccs ! T1-transf. of h_pq
!
!     Routine to calculate the energy (trivial: it is the SCF energy)
!
      procedure :: calc_energy => calc_energy_ccs
!
!     get Cholesky routines to calculate the occ/vir-occ/vir  blocks of the 
!     T1-transformed Cholesky vectors
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
!     Ground state solver routines (and helpers)
!
!     Note: while this solver is strictly uneccessary for CCS, where the solution
!           is trivial, it is inherited mostly unaltered by descendants (CCSD, CC2, 
!           etc.), where it serves an actual role.
!
      procedure :: ground_state_solver       => ground_state_solver_ccs
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
      procedure :: save_amplitudes => save_amplitudes_ccs
      procedure :: read_amplitudes => read_amplitudes_ccs
      procedure :: read_single_amplitudes => read_single_amplitudes_ccs
!
!     Routines to destroy amplitudes and omega 
!
      procedure :: destruct_amplitudes   => destruct_amplitudes_ccs
      procedure :: destruct_omega        => destruct_omega_ccs
!
!     Jacobian transformation routine 
!
      procedure :: jacobian_ccs_transformation => jacobian_ccs_transformation_ccs
!
!     Helper routines
!
!     Non-overridable, they will be used for contributions
!     to linear of higher order coupled cluster methods
!
      procedure, non_overridable :: jacobian_ccs_a1 => jacobian_ccs_a1_ccs 
      procedure, non_overridable :: jacobian_ccs_b1 => jacobian_ccs_b1_ccs
!
      procedure :: jacobi_test => jacobi_test_ccs
!
!     Jacobian transpose transformation routine (b^T -> b^T A, i.e., b -> A^T b)
!
      procedure :: jacobian_transpose_ccs_transformation => jacobian_transpose_ccs_transformation_ccs
!
!     Helper routines 
!
!     Non-overridable, they will be used for contributions
!     to linear of higher order coupled cluster methods
!
      procedure, non_overridable :: jacobian_transpose_ccs_a1 => jacobian_transpose_ccs_a1_ccs
      procedure, non_overridable :: jacobian_transpose_ccs_b1 => jacobian_transpose_ccs_b1_ccs

!     Excited state driver & solver 
!
      procedure                  :: excited_state_driver => excited_state_driver_ccs 
      procedure, non_overridable :: excited_state_solver => excited_state_solver_ccs
!
!     Helper routines 
!
      procedure :: transform_trial_vectors       => transform_trial_vectors_ccs
      procedure :: calculate_orbital_differences => calculate_orbital_differences_ccs ! Must be overwritten for CCSD 
!
      procedure, non_overridable :: find_start_trial_indices => find_start_trial_indices_ccs
      procedure, non_overridable :: initialize_trial_vectors => initialize_trial_vectors_ccs
      procedure, non_overridable :: trial_vectors_from_stored_solutions => trial_vectors_from_stored_solutions_ccs
!
      procedure, non_overridable :: solve_reduced_eigenvalue_equation => solve_reduced_eigenvalue_equation_ccs
      procedure, non_overridable :: construct_next_trial_vectors      => construct_next_trial_vectors_ccs
!
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
   module subroutine get_cholesky_ij_ccs(wf, L_ij_J, i_first, i_last, j_first, j_last)
!!
!!    Get Cholesky IJ
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!    Reads and T1-transforms the IA Cholesky vectors:
!!     
!!       L_ij_J_T1 = L_ij_J + sum_a t_aj * L_ia_J
!!
!!    Memory required in routine:
!!
!!       2*n_J*(i_length)*n_v     -> for reading L_ia_J contribution and reordering
!!       i_length = i_last - i_first + 1
!!
!!    Optional arguments: i_first, i_last, j_first, j_last can be used in order to restrict indices
!!       
      implicit none 
!
      class(ccs)               :: wf
      integer(i15), optional   :: i_first, j_first     ! First index (can differ from 1 when batching or for mlcc)
      integer(i15), optional   :: i_last, j_last      ! Last index (can differ from n_o when batching or for mlcc)
      real(dp), dimension(:,:) :: L_ij_J
!
      end subroutine get_cholesky_ij_ccs
!
!
   module subroutine get_cholesky_ia_ccs(wf, L_ia_J, i_first, i_last, a_first, a_last)
!!
!!    Get Cholesky IA
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!    Reads and T1-transforms IA Cholesky vectors
!!
!!       L_ia_J_T1 = L_ia_J (only reading necessary)
!!
!!    Memory required in routine:
!!
!!       No additional memory
!!
!!    Optional arguments: i_first, i_last, a_first, a_last can be used in order to restrict indices
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
      module subroutine get_cholesky_ai_ccs(wf, L_ai_J, i_first, i_last, a_first, a_last)
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
      module subroutine initialize_trial_vectors_ccs(wf)
!!
!!       Initialize trial vectors
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad
!!
!!       Initializes start trial vectors for the calculation of 
!!       singlet excited states and writes them to file 'trial_vecs'.
!!
!!       n start vectors are constructed by finding the n lowest orbital differences,      
!!       where n = n_singlet_states. Vector i has a 1.0D0 at the element corresponding to the i'th lowest
!!       orbital difference and 0.0d0 everywhere else
         implicit none
!
         class(ccs) :: wf
!
! 
      end subroutine initialize_trial_vectors_ccs
!
!
      module subroutine trial_vectors_from_stored_solutions_ccs(wf)
!!
!!
!!
      implicit none
!
      class(ccs) :: wf
!
      end subroutine trial_vectors_from_stored_solutions_ccs
!
!
      module subroutine find_start_trial_indices_ccs(wf, index_list)
!!
!!       Find indices for lowest orbital differences
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad
!!
!!
         implicit none
!
         class(ccs) :: wf
         integer(i15), dimension(wf%tasks%n_singlet_states,1), intent(inout) :: index_list
! 
      end subroutine find_start_trial_indices_ccs
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
      module subroutine transform_trial_vectors_ccs(wf, first_trial, last_trial)
!!
!!       Construct Right Transform of Jacobian
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad
!!
!!
         implicit none
!
         class(ccs) :: wf 
         integer(i15), intent(in) :: first_trial, last_trial       
!
      end subroutine transform_trial_vectors_ccs
!
!
      module subroutine jacobian_ccs_transformation_ccs(wf, c_a_i)
!!
!!       Jacobian transformation
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad
!!
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
         real(dp), dimension(reduced_dim, wf%tasks%n_singlet_states) :: solution_vectors_reduced
!
         integer(i15) :: reduced_dim
         integer(i15) :: n_new_trials
!
      end subroutine construct_next_trial_vectors_ccs
!
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
!     Set implemented methods
!
      wf%implemented%ground_state = .true.
      wf%implemented%excited_state = .true.
!
!     Read Hartree-Fock info from SIRIUS
!
      call wf%read_hf_info
!
!     Read Cholesky AO integrals and transform to MO basis
!
      call wf%read_transform_cholesky
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
            call wf%ground_state_solver
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
!        Excited state calculation requested
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
      if (wf%tasks%properties) then
!
!        Properties calculation requested
!
         if (wf%implemented%properties) then 
!
          !  call wf%properties
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
   subroutine destruct_amplitudes_ccs(wf)
!!
!!    Destruct Amplitudes (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjøsntad, May 2017
!!
!!    Deallocates the singles amplitudes.
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
!!    Deallocates the singles projection vector.
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
end module ccs_class
