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
      integer(i15) :: n_s2am = 0 ! Number of CC2 double amplitudes, only needed for excited state calculation.
!
      integer(i15) :: n_total_active_o = 0 
      integer(i15) :: n_total_active_v = 0 
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
!     ML helper routines
!
      procedure :: get_CC2_active_indices => get_CC2_active_indices_mlcc2
!
!     Omega
!
      procedure :: omega_mlcc2_a1   => omega_mlcc2_a1_mlcc2
      procedure :: omega_mlcc2_b1   => omega_mlcc2_b1_mlcc2
      procedure :: construct_omega  => construct_omega_mlcc2
!
      procedure :: calc_energy => calc_energy_mlcc2
!
!     Jacobian
!
      procedure :: initialize_excited_states     => initialize_excited_states_mlcc2
      procedure :: calculate_orbital_differences => calculate_orbital_differences_mlcc2
      procedure :: transform_trial_vectors       => transform_trial_vectors_mlcc2
      procedure :: jacobian_mlcc2_transformation => jacobian_mlcc2_transformation_mlcc2
      procedure :: jacobian_mlcc2_a1             => jacobian_mlcc2_a1_mlcc2
      procedure :: jacobian_mlcc2_b1             => jacobian_mlcc2_b1_mlcc2
      procedure :: jacobian_mlcc2_a2             => jacobian_mlcc2_a2_mlcc2
      procedure :: jacobian_mlcc2_b2             => jacobian_mlcc2_b2_mlcc2
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
   module subroutine cholesky_localization_mlcc2(wf)
!!
!!
      implicit none
!
      class(mlcc2) :: wf
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
    module subroutine omega_mlcc2_b1_mlcc2(wf, active_space)
! 
!     Omega B1
!     Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!   
!     Calculates the B1 term of omega, 
!   
!     B1: - sum_bjk u_jk^ab*g_jikb
!
      implicit none
!
      class(mlcc2)   :: wf
      integer(i15)   :: active_space 
      
   end subroutine omega_mlcc2_b1_mlcc2
!
!
   module subroutine construct_omega_mlcc2(wf)
!  
!        Construct Omega (CC2)
!        Written by Eirik F. Kjønstad and Sarai Folkestad, Apr 2017
!  
!        Constructs t2-amplitudes on the fly, according to the CC2
!        expression for the doubles amplitudes,
!  
!        t_ij^ab = - g_ai_bj / (e_a + e_b - e_i - e_j),
!  
!        where g_ai_bj are T1-transformed two-electron integrals 
!        and e_x is the orbital enegy of orbital x.
!         
!        The routine also sets up timing variables.    
!  
         implicit none 
!
         class(mlcc2) :: wf
!
   end subroutine construct_omega_mlcc2
!
!
   module subroutine initialize_excited_states_mlcc2(wf)
!!
      implicit none 
!    
      class(mlcc2) :: wf
!
   end subroutine initialize_excited_states_mlcc2
!  
!
   module subroutine calculate_orbital_differences_mlcc2(wf, orbital_diff)
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
         class(mlcc2) :: wf
!
         real(dp), dimension(wf%n_parameters, 1) :: orbital_diff
!
   end subroutine calculate_orbital_differences_mlcc2
!
!
   module subroutine transform_trial_vectors_mlcc2(wf, first_trial, last_trial)
!!
!!    Transformation Trial Vectors (CCSD)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!    Each trial vector in first_trial to last_trial is read from file and
!!    transformed before the transformed vector is written to file.
!!
!!    Singles and doubles part of the transformed vectors are written to 
!!    the same record in file transformed_vec, record length is n_parameters long.
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
   module subroutine jacobian_mlcc2_transformation_mlcc2(wf, c_a_i, c_aibj)
!!
!!    Jacobian transformation (MLCC2)
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
      class(mlcc2) :: wf 
!
!     Incoming vector c 
!
      real(dp), dimension(wf%n_v, wf%n_o) :: c_a_i  ! c_ai 
      real(dp), dimension(wf%n_s2am, 1)   :: c_aibj ! c_aibj    
   end subroutine jacobian_mlcc2_transformation_mlcc2
!
!
   module subroutine jacobian_mlcc2_a1_mlcc2(wf, rho_a_i, c_a_i)
!!
!!
!!    A1: 2*sum_BJck u_ik^ac*g_kc,JB*c_BJ - sum_Bjck u_kj^ca*g_kc,jB*c_BI
!!        - sum_BJck u_ik^ac*g_kB,Jc*c_BJ - sum_bJck u_ki^cb*g_kc,Jb*c_AJ 
!!
      implicit none
!
      class(mlcc2) :: wf 
!
!     Incoming vectors c and rho 
!
      real(dp), dimension(wf%n_v, wf%n_o) :: c_a_i 
      real(dp), dimension(wf%n_v, wf%n_o) :: rho_a_i 
!
   end subroutine jacobian_mlcc2_a1_mlcc2
!
!
   module subroutine jacobian_mlcc2_b1_mlcc2(wf, rho_a_i, c_ai_bj, n_active_o, &
                                              n_active_v, first_active_o, first_active_v)
!!
!!
!!    B1:   sum_ck F_kc*(2c_ai,ck - c_ak,ci) 
!!        - sum_ckj L_jIkc * c_aj,ck + sum_cbk L_Abkc * c_bi,ck
!!
!!
      implicit none
!  
      class(mlcc2) :: wf
!
      real(dp), dimension((n_active_v)*(n_active_o), (n_active_v)*(n_active_o))   :: c_ai_bj 
      real(dp), dimension(wf%n_v, wf%n_o)                                                                             :: rho_a_i
      integer(i15) :: n_active_o, n_active_v
      integer(i15) :: first_active_o, first_active_v
!
   end subroutine jacobian_mlcc2_b1_mlcc2
!
!
   module subroutine jacobian_mlcc2_a2_mlcc2(wf, rho_ai_bj, c_a_i,&
                                  n_active_o, n_active_v, first_active_o, first_active_v)
!!
!!
!!    A2:   sum_C g_aibC*c_Cj - sum_K g_aiKb*C_bK
!!
!!
      implicit none
!  
      class(mlcc2) :: wf
!
      real(dp), dimension((n_active_v)*(n_active_o), (n_active_v)*(n_active_o))   :: rho_ai_bj 
      real(dp), dimension(wf%n_v, wf%n_o)                                                                             :: c_a_i 
      integer(i15) :: n_active_o, n_active_v
      integer(i15) :: first_active_o, first_active_v
!
   end subroutine jacobian_mlcc2_a2_mlcc2
!
!
   module subroutine jacobian_mlcc2_b2_mlcc2(wf, rho_ai_bj, c_ai_bj,&
                                  n_active_o, n_active_v, first_active_o, first_active_v)
!!
!!
!!    B2:   ε_ij^ab*c_ai,bj
!!
!!
      implicit none
!  
      class(mlcc2) :: wf
!
      real(dp), dimension((n_active_v)*(n_active_o), (n_active_v)*(n_active_o))   :: c_ai_bj 
      real(dp), dimension((n_active_v)*(n_active_o), (n_active_v)*(n_active_o))   :: rho_ai_bj
      integer(i15) :: n_active_o, n_active_v
      integer(i15) :: first_active_o, first_active_v 
!
!
   end subroutine jacobian_mlcc2_b2_mlcc2
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
      integer(i15) :: i, j
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
      wf%implemented%ground_state = .true.
      wf%implemented%excited_state = .true.
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
      call wf%initialize_omega
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
!
!  :::::::::::::::::::::::::::::::::::::::::
!  -::- Class subroutines and functions -::- 
!  :::::::::::::::::::::::::::::::::::::::::
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
      integer(i15) :: active_space
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
      call allocator(L_IA_J, (wf%n_o)*(wf%n_v), wf%n_J)
      L_IA_J = zero
!
!     Get the Cholesky vector L_ia_J 
!
      call wf%get_cholesky_ia(L_IA_J)
!
!     Allocate g_ia_jb = g_iajb and set it to zero
!
      call allocator(g_IA_JB, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
      call deallocator(L_IA_J, (wf%n_o)*(wf%n_v), wf%n_J)
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
!     Loop over active spaces
!
      do active_space = 1, wf%n_active_spaces
!
!        Set ML variables
!
         call wf%get_CC2_active_indices(first_active_v, first_active_o, active_space)
!
         n_active_o = wf%n_CC2_o(active_space,1) 
         n_active_v = wf%n_CC2_v(active_space,1)
!
         last_active_o = first_active_o + n_active_o - 1
         last_active_v = first_active_v + n_active_v - 1 
!
!        Prepare for batching over index a
!  
         required = (2*(wf%n_v)*(wf%n_o)*(wf%n_J) &                            ! Needed for g_aibj  
                  + 2*((wf%n_v)**2)*(wf%n_J) &               ! and 's2' amplitudes  
                  + 2*(wf%n_v)**2*(wf%n_o)**2)                                !
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
!        Loop over the number of a batches 
!
         do a_batch = 1, n_batch
!
!           For each batch, get the limits for the a index 
!
            call batch_limits(a_first, a_last, a_batch, max_batch_length, batch_dimension)
!
!           a is active index, and thus a_first and a_last must be displaced
!
            a_first  = a_first + (first_active_v - 1)
            a_last   = a_last  + (first_active_v - 1)
!
            if (a_last .gt. last_active_v) a_last = last_active_v
!
            a_length = a_last - a_first + 1 

!           :: Calculate cc2 doubles amplitudes ::
!
            call allocator(L_ai_J, n_active_o*n_active_v, wf%n_J)
            L_ai_J = zero
!
            call wf%get_cholesky_ai(L_ai_J, first_active_v, last_active_v, first_active_o, last_active_o)
!
            call allocator(L_ia_J, (n_active_o)*(n_active_v), wf%n_J)
!
!           reorder and constrain L_bi_J
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
            call deallocator(L_ai_J, n_active_o*n_active_v, wf%n_J)
!
            call allocator(s_ia_jb, (n_active_o)*a_length, (n_active_o)*n_active_v)
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
            call deallocator(L_ia_J, (n_active_o)*(n_active_v), wf%n_J)
!
!           Add the rest of the correlation energy E = E + sum_aibj (s_ij^ab ) L_iajb
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
            call deallocator(s_ia_jb, a_length*n_active_o, n_active_o*n_active_v)
!
         enddo ! End of batching
      enddo ! End loop over active spaces
!
      call deallocator(g_ia_jb, wf%n_o*wf%n_v, wf%n_o*wf%n_v)
!
   end subroutine calc_energy_mlcc2
!
!
   subroutine initialize_orbital_info_mlcc2(wf)
!!
!!    Initialize orbital information
!!    Written by Sarai D. Folkestad, June 2017
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
!!    Destruct orbital info
!!    Written by Sarai D. Folkestad, June 2017
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
!
   subroutine get_CC2_active_indices_mlcc2(wf, first_o, first_v, active_space)
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
   integer(i15) :: active_space
!
   integer(i15) :: i
!
   first_o = 1
   first_v = 1
!
   do i = 1, active_space - 1
      first_o = first_o + wf%n_CC2_o(i,1)
      first_v = first_v + wf%n_CC2_v(i,1)
   enddo
!
   end subroutine get_CC2_active_indices_mlcc2
!
!
end module mlcc2_class