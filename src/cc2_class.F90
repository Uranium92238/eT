module cc2_class
!
!!
!!            Coupled cluster perturbative doubles (CC2) class module                                 
!!         Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017         
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
!  The ancestor class module (CCS)
!
   use ccs_class
!
!
   implicit none 
!
!  :::::::::::::::::::::::::::::::::::::
!  -::- Definition of the CC2 class -::-
!  ::::::::::::::::::::::::::::::::::::: 
!
   type, extends(ccs) :: cc2
!
      real(dp), dimension(:,:), allocatable :: s2am
!
      integer(dp) :: n_s2am
!
   contains 
!
!     Initialization and driver routines
!
      procedure :: init => init_cc2
!
!     Routines to construct the projection vector (omega)
!
      procedure :: construct_omega => construct_omega_cc2
!
!     Helper routines for construct_omega
!
      procedure :: omega_cc2_a1 => omega_cc2_a1_cc2 
      procedure :: omega_cc2_b1 => omega_cc2_b1_cc2 
      procedure :: get_s2am => get_s2am_cc2  
!
      procedure ::read_cc2_double_amplitudes => read_cc2_double_amplitudes_cc2  
      procedure ::read_amplitudes            => read_amplitudes_cc2  
      procedure ::save_amplitudes            => save_amplitudes_cc2  
      procedure ::destruct_s2am              => destruct_s2am_cc2 
!
!     Ground state solver helper routines
!
      procedure :: calc_energy => calc_energy_cc2
!
!     Jacobian
!
      procedure :: jacobian_cc2_transformation      => jacobian_cc2_transformation_cc2
      procedure :: cvs_jacobian_cc2_transformation  => cvs_jacobian_cc2_transformation_cc2
!
      procedure :: cvs_rho_ai_bj_projection           => cvs_rho_ai_bj_projection_cc2
!
      procedure :: jacobian_cc2_a1                  => jacobian_cc2_a1_cc2
      procedure :: jacobian_cc2_b1                  => jacobian_cc2_b1_cc2
      procedure :: jacobian_cc2_a2                  => jacobian_cc2_a2_cc2
      procedure :: jacobian_cc2_b2                  => jacobian_cc2_b2_cc2
!
!     Excited states
!
      procedure :: initialize_excited_states     => initialize_excited_states_cc2
      procedure :: calculate_orbital_differences => calculate_orbital_differences_cc2
      procedure :: transform_trial_vectors       => transform_trial_vectors_cc2
      procedure :: cvs_residual_projection       => cvs_residual_projection_cc2 
!
   end type cc2
!
!  ::::::::::::::::::::::::::::::::::::::::::::::::::::
!  -::- Interface to the submodule routines of CCS -::- 
!  ::::::::::::::::::::::::::::::::::::::::::::::::::::
!
   interface
!
!     -::- Omega interface -::-
!     :::::::::::::::::::::::::
!
      module subroutine construct_omega_cc2(wf)
!!
!!        Construct Omega 
!!        Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!        Directs the construction of the projection vector < mu | exp(-T) H exp(T) | R >
!!        for the current amplitudes of the object wf 
!!
         implicit none
!
         class(cc2) :: wf 
!
      end subroutine construct_omega_cc2
!
!
      module subroutine omega_cc2_a1_cc2(wf)
!!
!!        Omega A1
!!        Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!  
!!        Calculates the A1 term of omega, 
!!  
!!        A1: sum_ckd g_adkc * u_ki^cd,
!! 
!!        and adds it to the projection vector (omega1) of
!!        the wavefunction object wf
!!
!!        u_ki^cd = 2*t_ki^cd - t_ik^cd 
!!
         implicit none
!
         class(cc2) :: wf
!
      end subroutine omega_cc2_a1_cc2
!
!
      module subroutine omega_cc2_b1_cc2(wf)
!!
!!        Omega B1
!!        Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!        Calculates the B1 term of omega, 
!!
!!        B1: - sum_ckl u_kl^ac * g_kilc,
!! 
!!        and adds it to the projection vector (omeg1) of
!!        the wavefunction object wf
!!
!!        u_kl^ac = 2*t_kl^ac - t_lk^ac 
!!
         implicit none
!
         class(cc2) :: wf 
!
      end subroutine omega_cc2_b1_cc2
!
!
      module subroutine get_s2am_cc2(wf, s_ia_jb, b_first, b_length)
!!
!!       Get S_2 amplitudes, 
!!       Written by Sarai D. Folkestad, July 2017 
!!
!!       Construct
!!
!!          s_ai_bj = - 1/ε_ij^ab * g_aibj,
!!
!!       while batching over b.
!!
         implicit none
!
         class(cc2) :: wf
! 
         integer(i15) :: b_first, b_length
         real(dp), dimension((wf%n_v)*(wf%n_o), b_length*(wf%n_o)) :: s_ia_jb
!
      end subroutine get_s2am_cc2
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
      module subroutine jacobian_cc2_transformation_cc2(wf, c_a_i, c_aibj)
!!
!!       Jacobian transformation (CC2)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, June 2017
!!
!!       Directs the transformation by the CCSD Jacobi matrix,
!!
!!          A_mu,nu = < mu | exp(-T) [H, tau_nu] exp(T) | nu >,
!!
!!       where the basis employed for the brackets is biorthonormal. 
!!       The transformation is rho = A c, i.e., 
!!
!!          rho_mu = (A c)_mu = sum_ck A_mu,ck c_ck 
!!                  + 1/2 sum_ckdl A_mu,ckdl c_ckdl (1 + delta_ck,dl).
!!
!!       On exit, c is overwritten by rho. That is, c_a_i = rho_a_i,
!!       and c_aibj = rho_aibj. 
!!
         implicit none
!
         class(cc2) :: wf 
!
!        Incoming vector c 
!
         real(dp), dimension(wf%n_v, wf%n_o) :: c_a_i  ! c_ai 
         real(dp), dimension(wf%n_s2am, 1)   :: c_aibj ! c_aibj  
!  
      end subroutine jacobian_cc2_transformation_cc2
!
!
      module subroutine cvs_jacobian_cc2_transformation_cc2(wf, c_a_i, c_aibj)
!!
!!    Jacobian transformation (CC2)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, June 2017
!!
!!    Directs the transformation by the CCSD Jacobi matrix,
!!
!!       A_mu,nu = < mu | exp(-T) [H, tau_nu] exp(T) | nu >,
!!
!!    where the basis employed for the brackets is biorthonormal. 
!!    The transformation is rho = A c, i.e., 
!!
!!       rho_mu = (A c)_mu = sum_ck A_mu,ck c_ck 
!!                  + 1/2 sum_ckdl A_mu,ckdl c_ckdl (1 + delta_ck,dl).
!!
!!    On exit, c is overwritten by rho. That is, c_a_i = rho_a_i,
!!    and c_aibj = rho_aibj. 
!!
      implicit none
!
      class(cc2) :: wf 
!
!     Incoming vector c 
!
      real(dp), dimension(wf%n_v, wf%n_o) :: c_a_i  ! c_ai 
      real(dp), dimension(wf%n_s2am, 1)   :: c_aibj ! c_aibj     
!
   end subroutine cvs_jacobian_cc2_transformation_cc2
!
!
      module subroutine jacobian_cc2_a1_cc2(wf, rho_a_i, c_a_i)
!!
!!       Jacobian tem A1
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, June 2017
!!
!!       Calculates the A1 contribution to the jacobi transformation,
!!
!!          A1: 2*sum_BJck u_ik^ac*g_kc,JB*c_BJ - sum_Bjck u_kj^ca*g_kc,jB*c_BI
!!            - sum_BJck u_ik^ac*g_kB,Jc*c_BJ - sum_bJck u_ki^cb*g_kc,Jb*c_AJ, 
!!
!!       with, 
!!
!!       u_ik^ac = 2*s_ik^ac - 2*s_ik^ca,
!!
!!       which is constructed while batching over c
!!
         implicit none
!
         class(cc2) :: wf 
!
!        Incoming vectors c and rho 
!
         real(dp), dimension(wf%n_v, wf%n_o) :: c_a_i 
         real(dp), dimension(wf%n_v, wf%n_o) :: rho_a_i 
!
      end subroutine jacobian_cc2_a1_cc2
!
!
      module subroutine jacobian_cc2_b1_cc2(wf, rho_a_i, c_ai_bj)
!!
!!       Jacobian tem B1
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, June 2017
!!
!!       Calculates the B1 contribution to the jacobi transformation,
!!
!!       B1:   sum_ck F_kc*(2c_ai,ck - c_ak,ci) 
!!           - sum_ckj L_jIkc * c_aj,ck + sum_cbk L_Abkc * c_bi,ck
!!
!!
!!       L_Abkc is constructed while batching over A.
!!
         implicit none
!  
         class(cc2) :: wf
!
         real(dp), dimension(:,:)            :: c_ai_bj 
         real(dp), dimension(wf%n_v, wf%n_o) :: rho_a_i
!
      end subroutine jacobian_cc2_b1_cc2
!
!
      module subroutine jacobian_cc2_a2_cc2(wf, rho_ai_bj, c_a_i)
!!
!!       Jacobian tem A2
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, June 2017
!! 
!!       Calculates the A2 contribution to the jacobi transformation,
!! 
!!          A2:   sum_C g_ai,bC * c_Cj - sum_K g_ai,Kj * C_bK.
!! 
!!       g_ai,bC is constructed in batches of C.
!! 
         implicit none
!     
         class(cc2) :: wf
!  
         real(dp), dimension(:,:)            :: rho_ai_bj 
         real(dp), dimension(wf%n_v, wf%n_o) :: c_a_i 
!
      end subroutine jacobian_cc2_a2_cc2
!
!
      module subroutine jacobian_cc2_b2_cc2(wf, rho_ai_bj, c_ai_bj)
!!
!!       Jacobian tem B2
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, June 2017
!!
!!       Calculates the B2 contribution to the jacobi transformation,
!!
!!          B2:   ε_ij^ab*c_ai,bj.
!!
!!
         implicit none
!  
         class(cc2) :: wf
!
         real(dp), dimension(:,:)   :: c_ai_bj 
         real(dp), dimension(:,:)   :: rho_ai_bj
!
!
      end subroutine jacobian_cc2_b2_cc2
!
!
      module subroutine cvs_rho_ai_bj_projection_cc2(wf, vec_ai_bj)
!!
!!       Rho projection for CVS (CC2),
!!       Written by Sarai D. Folkestad, Aug. 2017
!!
         implicit none
!
         class(cc2) :: wf
         real(dp), dimension(:, :) :: vec_ai_bj
!
      end subroutine cvs_rho_ai_bj_projection_cc2
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
      module subroutine initialize_excited_states_cc2(wf)
!!
!!       Initialize excited states
!!       Written by Sarai D. Folkestad, June 2017
!!
!!       Calculates and sets n_s2am, and updates n_parameters
!!       for excited state calculation
!! 
         implicit none 
!    
         class(cc2) :: wf
!
      end subroutine initialize_excited_states_cc2
!  
!

      module subroutine calculate_orbital_differences_cc2(wf, orbital_diff)
!!
!!       Calculate Orbital Differences (CCSD)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad May 2017
!!
!!       Calculates orbital differences
!!
!!          1) ε_I^A = ε_A - ε_I
!!          2) ε_ij^ab = ε_a + ε_b - ε_i - ε_j (for active spaces only)
!!
!!       and puts them in orbital_diff, which is a vector of length n_parameters.        
!!
         implicit none
!
         class(cc2) :: wf
!
         real(dp), dimension(wf%n_parameters, 1) :: orbital_diff
!
      end subroutine calculate_orbital_differences_cc2
!
!
      module subroutine transform_trial_vectors_cc2(wf, first_trial, last_trial)
!!
!!       Transformation Trial Vectors (CC2)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!! 
!!       Each trial vector in first_trial to last_trial is read from file and
!!       transformed before the transformed vector is written to file.
!! 
!!       Singles and doubles part of the transformed vectors are written to 
!!       the same record in file transformed_vec, record length is n_parameters long.
!!
         implicit none
!
         class(cc2) :: wf
!
         integer(i15), intent(in) :: first_trial, last_trial ! Which trial_vectors we are to transform
!
      end subroutine transform_trial_vectors_cc2
!
!
      module subroutine cvs_residual_projection_cc2(wf, residual)
!!
!!       Residual projection (CC2), 
!!       Written by Sarai D. Folkestad Aug. 2017    
!!
         implicit none
!
         class(cc2) :: wf
         real(dp), dimension(wf%n_parameters, 1) :: residual
!
      end subroutine cvs_residual_projection_cc2
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
   subroutine init_cc2(wf)
!!
!!     Initialize CC2 object
!!     Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!     Performs the following tasks
!!
!!     1. Sets HF orbital and energy information by reading from file (read_hf_info)
!!     2. Transforms AO Cholesky vectors to MO basis and saves to file (read_transform_cholesky)
!!     3. Allocates the Fock matrix and sets it to zero (note: the matrix is constructed in 
!!        the descendant classes) 
!!     4. Allocates the singles amplitudes and sets them to zero, and sets associated properties 
!!     5. Allocate Omega vector
!!
      implicit none
!
      class(cc2)  :: wf
!
      integer(i15) :: unit_input = -1
!
!     Set model name
!
      wf%name = 'CC2'
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
!     Set implemented
!
      wf%implemented%ground_state         = .true.
      wf%implemented%excited_state        = .true.
      wf%implemented%core_excited_state   = .true.
!
!     Read calculation tasks from input file eT.inp
!     
      call wf%calculation_reader(unit_input)
!
!     Close input file
!
      close(unit_input)
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
!     Allocate Fock matrix and set to zero
!
      call wf%initialize_fock_matrix
!
      wf%n_parameters = wf%n_t1am
!
!     Initialize omega vector
!
      call wf%initialize_omega
!
   end subroutine init_cc2
!
!
!  :::::::::::::::::::::::::::::::::::::::::
!  -::- Class subroutines and functions -::- 
!  :::::::::::::::::::::::::::::::::::::::::
!
   subroutine calc_energy_cc2(wf)
!!
!!    Calculate Energy (CC2)
!!
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!    Calculates the CC2 energy, 
!!
!!    E_CC2 = E_HF + sum_aibj L_iajb*(t_ij^ab + t_i^a*t_j^b),
!!   
!!    with t_ij^ab = - g_aibj/(e_a + e_b - e_i - e_j) where 
!!    g_aibj are T1-transformed integrals.
!!    Batching over a.
!!
   implicit none
!
   class(cc2) :: wf
!
!     Integrals
!
      real(dp), dimension(:,:), allocatable :: L_bj_J
      real(dp), dimension(:,:), allocatable :: L_ia_J
      real(dp), dimension(:,:), allocatable :: g_ia_bj ! = g_aibj
      real(dp), dimension(:,:), allocatable :: g_ia_jb 
!
!     t2 amplitudes
!
      real(dp), dimension(:,:), allocatable :: t_ia_bj ! = g_aibj/(e_a + e_b - e_i - e_j)
!
!     Batching variables
!  
      integer(i15) :: a_batch, a_first, a_last, a_length
      integer(i15) :: required, available, n_batch, batch_dimension, max_batch_length
!
!     Indices
!
      integer(i15) :: a = 0, b = 0
      integer(i15) :: i = 0, j = 0
!
      integer(i15) :: ai = 0, bj = 0
      integer(i15) :: ia = 0, ib = 0, jb = 0, ja = 0
!
      integer(i15) :: aibj = 0
!
!     Prepare for batching over index a
!  
      required = (2*(wf%n_v)*(wf%n_o)*(wf%n_J) &                           !    
               + 2*(wf%n_v)*(wf%n_o)*(wf%n_J) &                            ! Needed for g_aibj  
               + 2*((wf%n_v)**2)*(wf%n_J) + ((wf%n_o)**2)*(wf%n_J) &       ! and 't2' amplitudes  
               + 2*(wf%n_v)**2*(wf%n_o)**2)                                !
!      
      required = 4*required ! In words
      available = get_available()
!
      batch_dimension  = wf%n_v ! Batch over the virtual index a
      max_batch_length = 0      ! Initilization of unset variables 
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
         a_length = a_last - a_first + 1 
!
!        :: Calculate cc2 doubles amplitudes ::
!
!        Allocate L_bj_J and L_ia_J (= reordering of L_bj_J constrained to the batch)
!
         call allocator(L_bj_J, (wf%n_v)*(wf%n_o), wf%n_J)
         call allocator(L_ia_J, a_length*(wf%n_o), wf%n_J)
         L_bj_J = zero
         L_ia_J = zero
!
         call wf%get_cholesky_ai(L_bj_J)
!
!        Create L_ia_J
!
         do a = 1, a_length
            do i = 1, wf%n_o
               do J = 1, wf%n_J
!
!                 Calculate compound indices
!
                  ia = index_two(i, a, wf%n_o)
                  ai = index_two(a + a_first - 1, i, wf%n_v)
!
                  L_ia_J(ia, J) = L_bj_J(ai, J)
!
               enddo
            enddo
         enddo
!
!        Allocate g_ia_bj
!
         call allocator(g_ia_bj, a_length*(wf%n_o), (wf%n_o)*(wf%n_v))
!
!        Construct integral g_ia_bj (= g_aibj for the batch)
!
         call dgemm('N','T',            &
                     a_length*(wf%n_o), &
                     (wf%n_o)*(wf%n_v), &
                     wf%n_J,            &
                     one,               &
                     L_ia_J,            &
                     a_length*(wf%n_o), &
                     L_bj_J,            &
                     (wf%n_o)*(wf%n_v), &
                     zero,              &
                     g_ia_bj,           &
                     a_length*(wf%n_o))
!
!        L_bj_J and L_ia_J
!
         call deallocator(L_bj_J, (wf%n_v)*(wf%n_o), wf%n_J)
         call deallocator(L_ia_J, a_length*(wf%n_o), wf%n_J)
!
!        :: Construct the needed integrals for the enegry ::
!
!        Allocate t_ia_bj
!
         call allocator(t_ia_bj, a_length*(wf%n_o), (wf%n_o)*(wf%n_v))
!
!        Create t2 amplitudes
!
         do a = 1, a_length
            do i = 1, wf%n_o
!
               ia = index_two(i, a, wf%n_o)
!
               do j = 1, wf%n_o
                  do b = 1, wf%n_v
!
                     bj = index_two(b, j, wf%n_v)
!
                     t_ia_bj(ia, bj) = - g_ia_bj(ia, bj)/(wf%fock_diagonal(a + wf%n_o, 1) &
                                          + wf%fock_diagonal(b + wf%n_o, 1) &
                                          - wf%fock_diagonal(i, 1) - wf%fock_diagonal(j, 1))
!
                  enddo
               enddo
            enddo
         enddo
!
!        Deallocate g_ia_bj
!
         call deallocator(g_ia_bj, a_length*(wf%n_o), (wf%n_o)*(wf%n_v))
!
!        Allocate the Cholesky vector L_ia_J = L_ia^J and set to zero 
!
         call allocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
         L_ia_J = zero
!
!         Get the Cholesky vector L_ia_J 
!
         call wf%get_cholesky_ia(L_ia_J)
!
!        Allocate g_ia_jb = g_iajb and set it to zero
!
         call allocator(g_ia_jb, (wf%n_o)*a_length, (wf%n_o)*(wf%n_v))
         g_ia_jb = zero
!
!        Calculate the integrals g_ia_jb from the Cholesky vector L_ia_J 
!
         call dgemm('N','T',                                   &
                     (wf%n_o)*a_length,                        &
                     (wf%n_o)*(wf%n_v),                        &
                     wf%n_J,                                   &
                     one,                                      &
                     L_ia_J(index_two(1, a_first, wf%n_o), 1), &
                     (wf%n_o)*(wf%n_v),                        &
                     L_ia_J,                                   &
                     (wf%n_o)*(wf%n_v),                        &
                     zero,                                     &
                     g_ia_jb,                                  &
                     (wf%n_o)*a_length)
!
!     Deallocate the Cholesky vector L_ia_J 
!
      call deallocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!     Set the initial value of the energy 
!
      wf%energy = wf%scf_energy
!
!     Add the correlation energy E = E + sum_aibj (t_ij^ab + t_i^a t_j^b) L_iajb
!
         do i = 1, wf%n_o
            do a = 1, a_length
!
               ai = index_two(a, i, a_length)
               ia = index_two(i, a, wf%n_o)
!
               do j = 1, wf%n_o
!
                  ja = index_two(j, a, wf%n_o)
!
                 do b = 1, wf%n_v
! 
                     bj = index_two(b, j, wf%n_v)
                     jb = index_two(j, b, wf%n_o)
                     ib = index_two(i, b, wf%n_o)
!
                     aibj = index_packed(ai, bj)
!
!                    Add the correlation energy 
!
                     wf%energy = wf%energy + &
                     (two*g_ia_jb(ia,jb) - g_ia_jb(ja, ib))*(t_ia_bj(ia, bj) + (wf%t1am(a,i))*(wf%t1am(b,j)))
                              
!
                  enddo
               enddo
            enddo
         enddo
!
!        Deallocate g_ia_jb
!
         call deallocator(g_ia_jb, (wf%n_o)*a_length, (wf%n_o)*(wf%n_v))
!
!        Deallocate t_ia_bj
!  
         call deallocator(t_ia_bj, a_length*(wf%n_o), (wf%n_o)*(wf%n_v))
!
      enddo ! End of batching
!
   end subroutine calc_energy_cc2
!
!
   subroutine save_amplitudes_cc2(wf)
!!
!!    Save Amplitudes (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjøsntad, May 2017
!!
!!    Store the amplitudes to disk (T1AM)
!!
      implicit none 
!
      class(cc2) :: wf
!
      integer(i15) :: unit_t1am = -1 
      integer(i15) :: unit_s2am = -1 
!
      real(dp), dimension(:,:), allocatable :: s_ia_jb
      real(dp), dimension(:,:), allocatable :: s2am
!
      integer(i15) :: a = 0, i = 0, ai = 0, ia = 0, b = 0, j = 0, bj = 0, jb = 0
      integer(i15) :: aibj = 0
!
!     Open amplitude files
!
      call generate_unit_identifier(unit_t1am)
      open(unit_t1am, file='t1am', status='unknown', form='unformatted')
      rewind(unit_t1am)
!
      call generate_unit_identifier(unit_s2am)
      open(unit_s2am, file='s2am', status='unknown', form='unformatted')
      rewind(unit_s2am)
!
!     Write t1 amplitudes
!
      write(unit_t1am) wf%t1am
!
!     Construct s2 amplitudes
!
      call allocator(s_ia_jb, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      call wf%get_s2am(s_ia_jb, 1, wf%n_v)
!  
!     Reorder and pack in
!
      call allocator(s2am, (wf%n_v)*(wf%n_o)*((wf%n_v)*(wf%n_o)+1)/2, 1)
!
      do i = 1, wf%n_o
         do a = 1, wf%n_v
            ai = index_two(a, i, wf%n_v)
            ia = index_two(i, a, wf%n_o)
            do j = 1, wf%n_o
               do b = 1, wf%n_v
!
                  bj = index_two(b, j, wf%n_v)
                  jb = index_two(j, b, wf%n_o)
!
                  aibj = index_packed(ai, bj)
!
                  s2am(aibj, 1) = s_ia_jb(ia, jb)
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(s_ia_jb, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     Write s2 amplitudes 
!
      write(unit_s2am) s2am
!
      call deallocator(s2am, (wf%n_v)*(wf%n_o)*((wf%n_v)*(wf%n_o)+1)/2, 1)
!
!     Close amplitude file
!
      close(unit_t1am)
      close(unit_s2am)
!
   end subroutine save_amplitudes_cc2
!
!
   subroutine destruct_s2am_cc2(wf)
!!
!!
   implicit none
!
   class(cc2) :: wf
!
   if (allocated(wf%s2am)) call deallocator(wf%s2am, wf%n_s2am, 1)
!
   end subroutine destruct_s2am_cc2
!
   subroutine read_amplitudes_cc2(wf)
!!
!!    Read Amplitudes (CC2)
!!    Written by Sarai D. Folkestad and Eirik F. Kjøsntad, May 2017
!!
!!    Reads the amplitudes from disk (T1AM)
!!
      implicit none 
!
      class(cc2) :: wf
!
      call wf%read_single_amplitudes
      call wf%read_cc2_double_amplitudes
!
   end subroutine read_amplitudes_cc2
!
!
   subroutine read_cc2_double_amplitudes_cc2(wf)
!!
!!    Read Amplitudes (CC2)
!!    Written by Sarai D. Folkestad and Eirik F. Kjøsntad, May 2017
!!
!!    Reads the amplitudes from disk (T1AM, S2AM)
!!
      implicit none 
!
      class(cc2) :: wf
!
      integer(i15) :: unit_s2am = -1 
!
      logical :: file_exists = .false.
!
!     Check to see whether file exists
!
      inquire(file='s2am',exist=file_exists)
!
      if (file_exists) then 
!
!        Open amplitude files if they exist
!
         call generate_unit_identifier(unit_s2am)
!
         open(unit_s2am, file='s2am', status='unknown', form='unformatted')
!
         rewind(unit_s2am)
!
!        Read from file & close
!
         wf%n_s2am = + ((wf%n_v)*(wf%n_o))&
                   *((wf%n_v )*(wf%n_o)+1)/2 
!
         if (.not. allocated(wf%s2am)) call allocator(wf%s2am, wf%n_s2am, 1) 
         read(unit_s2am) wf%s2am
!
         close(unit_s2am)
!
      else
!
         write(unit_output,'(t3,a)') 'Error: amplitude files do not exist.'
         stop
!
      endif
!
   end subroutine read_cc2_double_amplitudes_cc2
!
end module cc2_class