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
!     Excited state variables
!
      integer(i15)                           :: n_x2am = 0 
      real(dp), dimension(:,:), allocatable  :: x2am 
!
   contains
!
!     Initialization and driver routines
!
      procedure :: init => init_mlcc2
!
!     Orbital partitioning
!
      procedure :: orbital_partitioning          => orbital_partitioning_mlcc2
      procedure :: cholesky_decomposition        => cholesky_decomposition_mlcc2
      procedure :: cholesky_localization_drv     => cholesky_localization_drv_mlcc2
      procedure :: cholesky_orbitals             => cholesky_orbitals_mlcc2
      procedure :: cholesky_orbital_constructor  => cholesky_orbital_constructor_mlcc2
!
      procedure :: cnto_orbital_drv => cnto_orbital_drv_mlcc2
!
!     ML helper routines
!
      procedure :: get_CC2_active_indices => get_CC2_active_indices_mlcc2
      procedure :: get_CC2_n_active       => get_CC2_n_active_mlcc2
!
!     Amplitude routines
!
      procedure :: save_amplitudes              => save_amplitudes_mlcc2
      procedure :: read_amplitudes              => read_amplitudes_mlcc2
      procedure :: read_cc2_double_amplitudes   => read_cc2_double_amplitudes_mlcc2
      procedure :: destruct_x2am                => destruct_x2am_mlcc2
!
!     Omega
!
      procedure :: omega_mlcc2_a1   => omega_mlcc2_a1_mlcc2
      procedure :: omega_mlcc2_b1   => omega_mlcc2_b1_mlcc2
      procedure :: construct_omega  => construct_omega_mlcc2
      procedure :: get_s2am         => get_s2am_mlcc2
!
!     Ground state energy calculation
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
!     -::- Orbital partitioning submodule interface -::-
!     ::::::::::::::::::::::::::::::::::::::::::::::::::
!
      module subroutine orbital_partitioning_mlcc2(wf)
!!
!!       Orbital partitioning,
!!       Written by Sarai D. Folkestad, June 2017
!!
!!       Directs the partitioning for mlcc calculations.
!!
!!       So far only Cholesky decomposition is available. 
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
!!       Driver for Cholesky density decomposition.  
!!
!!       - Collects atom and ao-basis information.
!!       - Constructs occupied and vacant densities.
!!       - Constructs AO Fock matrix.  (This is currently an N^5 operation, should be optimized/removed)
!!       - By looping over active spaces, the occupied and virtual densities are Cholesky decomposed
!!         and the cholesky vectors are used to generate new localized MO's.
!!       - New orbitals are tested for orthonormality (Not implemented yet, only need overlap matrix from DALTON)    
!!
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
!!       Cholesky decomposes the density (occupied/virtual).
!!       Pivoting elements are chosen according to the active ao-list if it is pressent.
!!       If not, maximum diagonal elements are chosen as pivoting elements.
!!
!!       The Cholesky vectors are subtracted from the incoming density matrix, and the returned density can be furteh decomposed
!!       for inactive region.
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
!!       Makes the new MO's fromthe Cholesky vectors   
!!       - Transforms the AO fock matrix by the Cholesky vectors
!!       - Diagonalize the MO fock to get orbital energies and new orbitals.
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
!!       Constructs new localized orbitals (occupied/virtual) by 
!!       - Decomposing the density (occupied/virual)
!!       - Transforming Fock matrix with Cholesky vectors, and diagonalizing it.
!!         New orbitals are eigenvectors, orbital energies are eigenvectors.
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
      module function get_number_of_active_atoms(unit_cholesky_decomp, ml_level)
!!
!!       Get number of active atoms
!!       Written by Sarai D. Folkestad June 2017
!! 
!!       Reads cholesky.inp, and returns the number of atoms treated at ml_level of the CC hierarchy
!!       in active space of question.
!! 
         implicit none
!  
         integer(i15)      :: get_number_of_active_atoms
         integer(i15)      :: unit_cholesky_decomp
         character(len=5)  :: ml_level
!
      end function get_number_of_active_atoms
!
      module subroutine get_active_atoms(unit_cholesky_decomp, active_atoms, n_active_atoms, ml_level)
!!
!!       Get active atoms
!!       Written by Sarai D. Folkestad June 2017
!!
!!       Reads cholesky.inp, and returns the indices of the active atoms treated at ml_level of the CC hierarchy
!!       in active space of question.
!!
         implicit none
!
         integer(i15)      :: unit_cholesky_decomp
         integer(i15)      :: n_active_atoms
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
!!       Construct active ao index list,
!!       Written by Sarai Dery Folkestad, June 2017.
!!
!!       Constructs list of active ao's for cholesky decomposition.
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
!!       Read atom info,
!!       Written by Sarai Dery Folkestad, June 2017.
!!
!!       Reads atom info from DALTON generated file:
!!       Reads:
!!       - Number of nuclei
!!       - Number of AO's
!!
         implicit none
!
         integer(i15) :: n_nuclei,n_ao 
!
      end subroutine read_atom_info
!
!
      module subroutine read_center_info(n_nuclei, n_ao, n_ao_on_center, ao_center_info)
!!
!!       Read center info,
!!       Written by Sarai Dery Folkestad, June 2017.
!!
!!       Reads atom info from DALTON generated file:
!!       Reads:
!!          - Information of which ao's belong to which nuclei
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
      module subroutine cnto_orbital_drv_mlcc2(wf)
!!
!!       CNTO orbital driver,
!!       Written by Sarai D. Folkestad, June 2017.
!!
!!       A CCS calculation ground state and excited states is performed.
!!       The M and N matrices are then constructed, 
!! 
!!          M_ij = sum_a R1_ai*R1_aj + sum_a R2_ai*R2_aj + ...
!!          N_ab = sum_i R1_ai*R1_bi + sum_a R2_ai*R2_bi + ...
!!   
!!       where Ri_ai is the i'th single excitation vector obtained from the CCS calculation. 
!!       The transformation matrices for the occupied and virtual part
!!       are constructed by diagonalizing M and N. The number of active occupied
!!       and virtual orbitals are determined from δ_o and δ_v
!!
!!          1 - sum_i λ^o_i < δ_o
!!          1 - sum_i λ^v_i < δ_v
!!
!!       Where the orbitals of highest eigenvalues λ^o/λ^v are selected first.
!!
!!       Fock matrix is block diagonalized in active and inactive blocks in order to obtain 
!!       the orbitals and orbital energies used in the CC2 calculation.
!!

         implicit none 
!
         class(mlcc2) :: wf
!
      end subroutine cnto_orbital_drv_mlcc2
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
! 
!        Omega A1
!        Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!   
!        Calculates the A1 term of omega, 
!   
!        A1: sum_ckd g_adkc * u_ki^cd,
!  
!        and adds it to the projection vector (omega1) of
!        the wavefunction object wf
! 
!        u_ki^cd = 2*s_ki^cd - s_ik^cd 
! 
         implicit none
!
         class(mlcc2)   :: wf
!
      end subroutine omega_mlcc2_a1_mlcc2
!
!
      module subroutine omega_mlcc2_b1_mlcc2(wf)
! 
!        Omega B1
!        Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!   
!        Calculates the B1 term of omega, 
!   
!        B1: - sum_bjk u_jk^ab*g_jikb
!
         implicit none
!
         class(mlcc2)   :: wf
      
      end subroutine omega_mlcc2_b1_mlcc2
!
!
      module subroutine construct_omega_mlcc2(wf)
!  
!           Construct Omega (CC2)
!           Written by Eirik F. Kjønstad and Sarai Folkestad, Apr 2017
!  
!           Constructs t2-amplitudes on the fly, according to the CC2
!           expression for the doubles amplitudes,
!  
!           t_ij^ab = - g_ai_bj / (e_a + e_b - e_i - e_j),
!  
!           where g_ai_bj are T1-transformed two-electron integrals 
!           and e_x is the orbital enegy of orbital x.
!            
!           The routine also sets up timing variables.    
!  
            implicit none 
!
            class(mlcc2) :: wf
!
      end subroutine construct_omega_mlcc2
!
      module subroutine get_s2am_mlcc2(wf, s_ia_jb, b_first, b_length)
!!
!!       Batching over b
!!
         implicit none
!
         class(mlcc2) :: wf
! 
         integer(i15) :: b_first, b_length
         real(dp), dimension((wf%n_CC2_v)*(wf%n_CC2_o), b_length*(wf%n_CC2_o)) :: s_ia_jb
!
      end subroutine get_s2am_mlcc2
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
!!       Transformation Trial Vectors (CCSD)
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
         class(mlcc2) :: wf
!
         integer(i15), intent(in) :: first_trial, last_trial ! Which trial_vectors we are to transform
!
      end subroutine transform_trial_vectors_mlcc2
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
!!       Calculates the A2 contribution to the jacobi transformation,
!! 
!!          A2:   sum_C g_ai,bC * c_Cj - sum_K g_ai,Kj * C_bK.
!! 
!!       g_ai,bC is constructed in batches of C.
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
!!       Calculates the B2 contribution to the jacobi transformation,
!!
!!          B2:   ε_ij^ab*c_ai,bj.
!!
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
!        Do full space CC2 calculation
!
         write(unit_output,'(/t3,a50)')'Full CC2 requested, orbital partitioning skipped'
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
         call allocator(L_ai_J, n_active_o*n_active_v, wf%n_J)
         L_ai_J = zero
!    
         call wf%get_cholesky_ai(L_ai_J, first_active_v, last_active_v, first_active_o, last_active_o)
!
         call allocator(L_ia_J, (n_active_o)*(n_active_v), wf%n_J)
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
         call deallocator(s_ia_jb, a_length*n_active_o, n_active_o*n_active_v)
!
      enddo ! End of batching
!
      call deallocator(g_ia_jb, wf%n_o*wf%n_v, wf%n_o*wf%n_v)
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
      real(dp), dimension(:,:), allocatable :: s_ia_jb
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
      call allocator(s_ia_jb, (wf%n_CC2_v)*(wf%n_CC2_o), (wf%n_CC2_v)*(wf%n_CC2_o))
      call wf%get_s2am(s_ia_jb, wf%first_CC2_v, wf%first_CC2_v + wf%n_CC2_v - 1)
!  
!     Reorder and pack in
!
      call allocator(s2am, (wf%n_CC2_v)*(wf%n_CC2_o)*((wf%n_CC2_v)*(wf%n_CC2_o)+1)/2, 1)
!
      do i = 1, wf%n_CC2_o
         do a = 1, wf%n_CC2_v
            ai = index_two(a, i, wf%n_CC2_v)
            ia = index_two(i, a, wf%n_CC2_o)
            do j = 1, wf%n_CC2_o
               do b = 1, wf%n_CC2_v
!
                  bj = index_two(b, j, wf%n_CC2_v)
                  jb = index_two(j, b, wf%n_CC2_o)
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
      call deallocator(s_ia_jb, (wf%n_CC2_v)*(wf%n_CC2_o), (wf%n_CC2_v)*(wf%n_CC2_o))
!
!     Write s2 amplitudes 
!
      write(unit_x2am) s2am
!
      call deallocator(s2am, (wf%n_CC2_v)*(wf%n_CC2_o)*((wf%n_CC2_v)*(wf%n_CC2_o)+1)/2, 1)
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
         wf%n_x2am = + ((wf%n_CC2_v)*(wf%n_CC2_o))&
                   *((wf%n_CC2_v )*(wf%n_CC2_o)+1)/2 
!
         if (.not. allocated(wf%x2am)) call allocator(wf%x2am, wf%n_x2am, 1) 
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
   subroutine destruct_x2am_mlcc2(wf)
!!
!!
!!
   implicit none
!
   class(mlcc2) :: wf
!
   if (allocated(wf%x2am)) call deallocator(wf%x2am, wf%n_x2am, 1)
!
   end subroutine destruct_x2am_mlcc2
!
end module mlcc2_class