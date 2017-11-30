module mlccsd_class
!
!!
!!                 Multi-level CCSD (MLCCSD) class module                                
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
!  The ancestor class module (MLCC2)
!
   use mlcc2_class
!
   implicit none 
!
!  :::::::::::::::::::::::::::::::::::::::
!  -::- Definition of the MLCCSD class -::-
!  ::::::::::::::::::::::::::::::::::::::: 
!
   type, extends(mlcc2) :: mlccsd
!
!     ML-variables
!
      integer(i15) :: n_CCSD_o = 0
      integer(i15) :: n_CCSD_v = 0
!
      integer(i15) :: first_CCSD_o = 0
      integer(i15) :: first_CCSD_v = 0
!
      type(mlcc_orbitals) :: CCSD_orbitals
!
      real(dp), dimension(:,:), allocatable :: mo_coef_cc2_ccs ! MO coefficient matrix for mlccsd basis
      real(dp), dimension(:,:), allocatable :: T_o ! Occupied MO transformation matrix for mlccsd basis
      real(dp), dimension(:,:), allocatable :: T_v ! Virtual MO transformation matrix for mlccsd basis
!
!     Fock matrix variables
!
      real(dp), dimension(:,:), allocatable :: fock_diagonal_cc2_ccs  ! diagonal vector for mlccsd basis
!
!     Amplitude variables
!
      integer(i15) :: n_t2am = 0                    ! Number CCSD of doubles amplitudes
!
      real(dp), dimension(:,:), allocatable :: t2am ! Doubles amplitude vector
!
!     Schrödinger equation projection vector (the omega vector)
! 
!        < mu | exp(-T) H exp(T) | R >
!
      real(dp), dimension(:,:), allocatable :: omega2 ! Doubles vector
!
   contains
!
      procedure :: mlcc_reader       => mlcc_reader_mlccsd
      procedure :: read_orbital_info => read_orbital_info_mlccsd
!
      procedure :: init => init_mlccsd
!
!     Initialization routines
!
      procedure :: initialize_amplitudes   => initialize_amplitudes_mlccsd
      procedure :: initialize_omega        => initialize_omega_mlccsd
!
!     Energy calculation routine
!
      procedure :: calc_energy             => calc_energy_mlccsd
!
!     Orbital partitioning
!
      procedure :: orbital_partitioning               => orbital_partitioning_mlccsd
      procedure :: cholesky_localization_drv          => cholesky_localization_drv_mlccsd
      procedure :: cholesky_localization_CCSD_CC2_CCS => cholesky_localization_CCSD_CC2_CCS_mlccsd
      procedure :: cholesky_localization_CCSD_CCS     => cholesky_localization_CCSD_CCS_mlccsd
      procedure :: cholesky_localization_CCSD_CC2     => cholesky_localization_CCSD_CC2_mlccsd
      procedure :: construct_perturbative_doubles     => construct_perturbative_doubles_mlccsd
!
      procedure :: construct_MO_transformation_matrix => construct_MO_transformation_matrix_mlccsd
!
      procedure :: cnto_orbital_drv                   => cnto_orbital_drv_mlccsd
      procedure :: ccsd_cnto_lower_level_method       => ccsd_cnto_lower_level_method_mlccsd
      procedure :: ccsd_cnto_orbitals                 => ccsd_cnto_orbitals_mlccsd
      procedure :: print_orbital_info                 => print_orbital_info_mlccsd
!
!     Cholesky vector routines for CC2 space
!
      procedure :: read_transform_cholesky_for_CC2_amplitude   => read_transform_cholesky_for_CC2_amplitude_mlccsd
      procedure :: read_cholesky_ai_for_cc2_amplitudes         => read_cholesky_ai_for_cc2_amplitudes_mlccsd
      procedure :: read_cholesky_ia_for_cc2_amplitudes         => read_cholesky_ia_for_cc2_amplitudes_mlccsd
      procedure :: read_cholesky_ij_for_cc2_amplitudes         => read_cholesky_ij_for_cc2_amplitudes_mlccsd
      procedure :: read_cholesky_ab_for_cc2_amplitudes         => read_cholesky_ab_for_cc2_amplitudes_mlccsd
      procedure :: get_cholesky_ai_for_cc2_amplitudes          => get_cholesky_ai_for_cc2_amplitudes_mlccsd
!
!     ML helper routines
!
      procedure :: get_CC2_active_indices    => get_CC2_active_indices_mlccsd
      procedure :: get_CC2_n_active          => get_CC2_n_active_mlccsd
      procedure :: get_CCSD_active_indices   => get_CCSD_active_indices_mlccsd
      procedure :: get_CCSD_n_active         => get_CCSD_n_active_mlccsd
!
!     Omega routines
!
      procedure :: construct_omega  => construct_omega_mlccsd
!
      procedure :: omega_mlccsd_a1   => omega_mlccsd_a1_mlccsd
      procedure :: omega_mlccsd_b1   => omega_mlccsd_b1_mlccsd
!
      procedure :: get_mlccsd_x2am  => get_mlccsd_x2am_mlccsd
!
      procedure :: omega_mlccsd_a2  => omega_mlccsd_a2_mlccsd
      procedure :: omega_mlccsd_b2  => omega_mlccsd_b2_mlccsd
      procedure :: omega_mlccsd_c2  => omega_mlccsd_c2_mlccsd
      procedure :: omega_mlccsd_d2  => omega_mlccsd_d2_mlccsd
      procedure :: omega_mlccsd_e2  => omega_mlccsd_e2_mlccsd
!
!     Ground state solver routines
!
      procedure :: calc_ampeqs_norm          => calc_ampeqs_norm_mlccsd
      procedure :: new_amplitudes            => new_amplitudes_mlccsd
      procedure :: calc_quasi_Newton_doubles => calc_quasi_Newton_doubles_mlccsd
      procedure :: initialize_ground_state   => initialize_ground_state_mlccsd
!
!     Helper routines for onstruction, destruction, saving and reading
!
      procedure :: construct_orbital_coef_CC2_CCS     => construct_orbital_coef_CC2_CCS_mlccsd
      procedure :: construct_orbital_energy_CC2_CCS   => construct_orbital_energy_CC2_CCS_mlccsd
!
      procedure :: destruct_orbital_coef_CC2_CCS     => destruct_orbital_coef_CC2_CCS_mlccsd
      procedure :: destruct_orbital_energy_CC2_CCS   => destruct_orbital_energy_CC2_CCS_mlccsd
      procedure :: destruct_amplitudes               => destruct_amplitudes_mlccsd
      procedure :: destruct_double_amplitudes        => destruct_double_amplitudes_mlccsd
!
      procedure :: save_amplitudes                 => save_amplitudes_mlccsd
      procedure :: read_amplitudes                 => read_amplitudes_mlccsd
      procedure :: read_mlccsd_double_amplitudes   => read_mlccsd_double_amplitudes_mlccsd
      procedure :: read_double_amplitudes          => read_double_amplitudes_mlccsd
!
!     Jacobian 
!
      procedure :: jacobian_mlccsd_transformation => jacobian_mlccsd_transformation_mlccsd
!
      procedure :: jacobian_mlccsd_b2 => jacobian_mlccsd_b2_mlccsd
      procedure :: jacobian_mlccsd_c2 => jacobian_mlccsd_c2_mlccsd
      procedure :: jacobian_mlccsd_d2 => jacobian_mlccsd_d2_mlccsd
      procedure :: jacobian_mlccsd_e2 => jacobian_mlccsd_e2_mlccsd
      procedure :: jacobian_mlccsd_f2 => jacobian_mlccsd_f2_mlccsd
      procedure :: jacobian_mlccsd_g2 => jacobian_mlccsd_g2_mlccsd
      procedure :: jacobian_mlccsd_h2 => jacobian_mlccsd_h2_mlccsd
      procedure :: jacobian_mlccsd_i2 => jacobian_mlccsd_i2_mlccsd
      procedure :: jacobian_mlccsd_j2 => jacobian_mlccsd_j2_mlccsd
      procedure :: jacobian_mlccsd_k2 => jacobian_mlccsd_k2_mlccsd
!
!     Excited states
!
      procedure :: transform_trial_vectors   => transform_trial_vectors_mlccsd
      procedure :: initialize_excited_states => initialize_excited_states_mlccsd
      procedure :: print_excitation_vector   => print_excitation_vector_mlccsd
!
   end type mlccsd
!
!
   interface
!
!    -::- Input reader submodule interface -::-
!    :::::::::::::::::::::::::::::::::::::::::: 
!
      module subroutine mlcc_reader_mlccsd(wf, unit_input)
!!
!!
!!
         implicit none
!
         integer(i15) :: unit_input
!
         class(mlccsd)  :: wf
!
      end subroutine mlcc_reader_mlccsd
!
!
      module subroutine read_orbital_info_mlccsd(wf, unit_input)
!!
!!
          implicit none
!
          integer(i15)      :: unit_input
!
          class(mlccsd)      :: wf
!
      end subroutine read_orbital_info_mlccsd
!
   end interface
!
!
   interface
!
!    -::- Orbital partitioning submodule interface -::-
!    :::::::::::::::::::::::::::::::::::::::::::::::::: 
!
!
      module subroutine orbital_partitioning_mlccsd(wf)
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
         class(mlccsd) :: wf
!  
      end subroutine orbital_partitioning_mlccsd
!
!
      module subroutine cholesky_localization_drv_mlccsd(wf)
!!
!!       Cholesky orbital localization driver,
!!       Written by Sarai D. Folkestad, July 2017.
!!
!!       Driver for Cholesky density decomposition 
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
         class(mlccsd) :: wf
 !
      end subroutine cholesky_localization_drv_mlccsd
!
!
      module subroutine cholesky_localization_CCSD_CC2_CCS_mlccsd(wf, ao_center_info, n_ao_on_center,&
                                                       ao_fock, n_nuclei, unit_cholesky_decomp)
!!
!!       Cholesky orbital localization CCS/CC2/CCSD,
!!       Written by Sarai D. Folkestad, July 2017
!!
!!       Cholesky partitiining routine for CCS/CC2/CCSD calculation
!!
         implicit none
!
!        Input arguments
!
         class(mlccsd) :: wf
!
         integer(i15), dimension(wf%n_ao, 2)  :: ao_center_info
         integer(i15), dimension(n_nuclei, 1) :: n_ao_on_center
!
         real(dp), dimension(:,:) :: ao_fock
!
         integer(i15) :: n_nuclei
         integer(i15) :: unit_cholesky_decomp
!
      end subroutine cholesky_localization_CCSD_CC2_CCS_mlccsd
!
!
      module subroutine cholesky_localization_CCSD_CCS_mlccsd(wf, ao_center_info, n_ao_on_center,&
                                                       ao_fock, n_nuclei, unit_cholesky_decomp)
!!
!!       Cholesky orbital localization CCS/CCSD,
!!       Written by Sarai D. Folkestad, July 2017
!!
!!       Cholesky partitiining routine for CCS/CCSD calculation
!!
         implicit none
!
!        Input arguments
!
         class(mlccsd) :: wf
!
         integer(i15), dimension(wf%n_ao, 2)  :: ao_center_info
         integer(i15), dimension(n_nuclei, 1) :: n_ao_on_center
!
         real(dp), dimension(:,:) :: ao_fock
!
         integer(i15) :: n_nuclei
         integer(i15) :: unit_cholesky_decomp
!
      end subroutine cholesky_localization_CCSD_CCS_mlccsd
!
!
      module subroutine cholesky_localization_CCSD_CC2_mlccsd(wf, ao_center_info, n_ao_on_center,&
                                                       ao_fock, n_nuclei, unit_cholesky_decomp)
!!
!!       Cholesky orbital localization CC2/CCSD
!!       Written by Sarai D. Folkestad, July 2017
!!
!!       Cholesky partitiining routine for CC2/CCSD calculation
!!
         implicit none
!
!        Input arguments
!
         class(mlccsd) :: wf
!
         integer(i15), dimension(wf%n_ao, 2)  :: ao_center_info
         integer(i15), dimension(n_nuclei, 1) :: n_ao_on_center
!
         real(dp), dimension(:,:) :: ao_fock
!
         integer(i15) :: n_nuclei
         integer(i15) :: unit_cholesky_decomp
!
      end subroutine cholesky_localization_CCSD_CC2_mlccsd
!
!
      module subroutine construct_MO_transformation_matrix_mlccsd(wf)
!!
!!       Construct MO transformation matrix,
!!       Written by Sarai D. Fokestad, July 2017
!!
!!       Constructs transformation matrix,
!!
!!          T = (C_CCSD)^T * S * C_CC2
!!
!!       between CC2 basis and CCSD basis. 
!!       Needed for transforming s_ij_ab from CC2 to CCSD basis.
!!
         implicit none
!
         class(mlccsd) :: wf
!
      end subroutine construct_MO_transformation_matrix_mlccsd
!
!
      module subroutine cnto_orbital_drv_mlccsd(wf)
!!
!!       CNTO orbital driver,
!!       Written by Sarai D. Folkestad, July 2017.
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
         class(mlccsd) :: wf
!
      end subroutine cnto_orbital_drv_mlccsd
!
!
      module subroutine ccsd_cnto_lower_level_method_mlccsd(wf, cc2_n_parameters, cc2_n_x2am, n_cc2_o, n_cc2_v)
!!
!!    CNTO constructor (MLCCSD),
!!    Written by Sarai D. Folkestad, Aug. 2017
!!
!!    Constructs CNTOs and partitions orbital space
!!
      implicit none 
!
      class(mlccsd) :: wf
!
      integer(i15) :: cc2_n_x2am, cc2_n_parameters
      integer(i15) :: n_CC2_o, n_CC2_v
!
      end subroutine ccsd_cnto_lower_level_method_mlccsd
!
!
      module subroutine ccsd_cnto_orbitals_mlccsd(wf, cc2_n_parameters, cc2_n_x2am, n_cc2_o, n_cc2_v)
!!
!!       CNTO Oritals (MLCCSD),
!!       Written by Sarai D. Folkestad Aug. 2017
!!
!!       Constructs the CNTO orbitals based on exitation vectors from lower level method
!!    
         implicit none
!
         class(mlccsd) :: wf
!
         integer(i15) :: cc2_n_x2am, cc2_n_parameters
         integer(i15) :: n_CC2_o, n_CC2_v
!
      end subroutine ccsd_cnto_orbitals_mlccsd
!
!
      module subroutine print_orbital_info_mlccsd(wf)
!!
!!       Print CNTO info, 
!!       Written by Sarai D. Folkestad, Aug. 2017
!!
!!       Prints information on CNTO partitioning
!!
         implicit none 
!
         class(mlccsd) :: wf
!
      end subroutine print_orbital_info_mlccsd
!
!
   end interface
!
!
   interface
!
!     -::- Cholesky submodule interface -::-
!     ::::::::::::::::::::::::::::::::::::::
!
      module subroutine read_transform_cholesky_for_CC2_amplitude_mlccsd(wf)
!
         implicit none
!
         class(mlccsd) :: wf
!
      end subroutine read_transform_cholesky_for_CC2_amplitude_mlccsd
!
!
      module subroutine read_cholesky_ai_for_cc2_amplitudes_mlccsd(wf,L_ai_J, a_first, a_last, i_first, i_last)
!!
!!       Read Cholesky IA 
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!       Reads the MO Cholesky IA (occ-vir) vectors from file and
!!       places them in the incoming L_ia_J matrix
!!
!!
         implicit none
!
         class(mlccsd)            :: wf
         real(dp), dimension(:,:) :: L_ai_J
         integer(i15)             :: i_first, a_first     ! First index (can differ from 1 when batching or for mlcc)
         integer(i15)             :: i_last, a_last      ! Last index (can differ from n_o/n_v when batching or for mlcc)
!
!
      end subroutine read_cholesky_ai_for_cc2_amplitudes_mlccsd
!
!
      module subroutine get_cholesky_ai_for_cc2_amplitudes_mlccsd(wf, L_ai_J, a_first, a_last, i_first, i_last)
!!
!!       Get Cholesky AI
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!       Read and T1-transform Cholesky AI vectors:
!!    
!!        L_ai_J_T1 = L_ia_J - sum_j  t_aj*L_ji_J 
!!                            + sum_b  t_bi*L_ab_J
!!                            - sum_bj t_aj*t_bi*L_jb_J
!!
!!       Allocations in routine:
!!
!!       (1) n_J*(i_length)*(a_length) + 2*n_J*(a_length)*batch_length  ->  for L_ab_J contribution (batches of b)
!!       (2) n_J*(i_length)*n_v + 2*n_J*n_o*(i_length)                  ->  for L_ij_J contribution
!!       (3) 2*n_J*n_o*n_v                                              ->  for L_jb_J contribution
!!
!!       i_length = i_last - i_first + 1          
!!       a_length = a_last - a_first + 1          
!!
!!       (1) determines memory requirement. 
!!
         implicit none 
!
         class(mlccsd)            :: wf
         real(dp), dimension(:,:) :: L_ai_J
         integer(i15)             :: i_first, a_first     ! First index (can differ from 1 when batching or for mlcc)
         integer(i15)             :: i_last, a_last      ! Last index (can differ from n_o/n_v when batching or for mlcc)
!
      end subroutine get_cholesky_ai_for_cc2_amplitudes_mlccsd
!
!
    module subroutine read_cholesky_ab_for_cc2_amplitudes_mlccsd(wf, L_ab_J, a_first, a_last, b_first, b_last)
!!
!!    Read Cholesky AB 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!    Reads the MO Cholesky AB (vir-vir) vectors from file and
!!    places them in the incoming L_ab_J matrix, with batching 
!!    if necessary
!!
!!    Optional arguments: b_first, b_last, a_first, a_last can be used in order to restrict indices
!!
      implicit none
!
      class(mlccsd)                :: wf
      integer(i15), intent(in) :: a_first, b_first   ! First index (can differ from 1 when batching  or for mlcc)
      integer(i15), intent(in) :: a_last, b_last    ! Last index  (can differ from n_v when batching or for mlcc)
      real(dp), dimension(((a_last - a_first + 1)*(b_last - b_first + 1)), wf%n_J) :: L_ab_J ! L_ab^J
!
      end subroutine read_cholesky_ab_for_cc2_amplitudes_mlccsd
!
!!   
   module subroutine read_cholesky_ia_for_cc2_amplitudes_mlccsd(wf,L_ia_J, i_first, i_last, a_first, a_last)
!!
!!    Read Cholesky IA 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!    Reads the MO Cholesky IA (occ-vir) vectors from file and
!!    places them in the incoming L_ia_J matrix
!!
!!
!!    Optional arguments: i_first, i_last, a_first, a_last can be used in order to restrict indices
!!
      implicit none
!
      class(mlccsd)                 :: wf    
      integer(i15)                  :: i_first, a_first     ! First index (can differ from 1 when batching or for mlcc) 
      integer(i15)                  :: i_last, a_last      ! Last index (can differ from n_o when batching or for mlcc) 
      real(dp), dimension(:,:)      :: L_ia_J ! L_ia^J
!
   end subroutine read_cholesky_ia_for_cc2_amplitudes_mlccsd
!
!
   module subroutine read_cholesky_ij_for_cc2_amplitudes_mlccsd(wf,L_ij_J , i_first, i_last, j_first, j_last)
!!
!!    Read Cholesky IJ 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!    Reads the MO Cholesky IJ (occ-occ) vectors from file and 
!!    places them in the incoming L_ij_J matrix
!!
!!    Optional arguments: i_first, i_last, j_first, j_last can be used in order to restrict indices
!!
      implicit none
!
      class(mlccsd)            :: wf
      integer(i15)             :: i_first, j_first     ! First index (can differ from 1 when batching or for mlcc) 
      integer(i15)             :: i_last, j_last      ! Last index (can differ from n_o when batching or for mlcc)   
      real(dp), dimension(:,:) :: L_ij_J ! L_ij^J
!
   end subroutine read_cholesky_ij_for_cc2_amplitudes_mlccsd
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
      module subroutine initialize_omega_mlccsd(wf)
!
!        Initialize Omega (MLCCSD)
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!
!        Allocates the projection vector (omega1, omega2) and sets it
!        to zero.
!
         implicit none 
!
         class(mlccsd) :: wf
!
      end subroutine initialize_omega_mlccsd
!
!
      module subroutine construct_omega_mlccsd(wf)
!! 
!!       Construct Omega (MLCCSD)
!!       Written by Eirik F. Kjønstad and Sarai Folkestad, Apr 2017
!!
!!       Constructs the MlCC2 omega.
!! 
!!       s2-amplitudes are constructed on the fly, according to the CC2
!!       expression for the doubles amplitudes. 
!!
!!       Calculated by looping over active spaces, 
!!       Adding the omega contribution from each active space in turn.
!! 
         implicit none 
!
         class(mlccsd) :: wf
!
      end subroutine construct_omega_mlccsd
!
!
!
   end interface
!
!
   interface
!
!     -::- Omega CC2 submodule interface -::-
!     :::::::::::::::::::::::::::::::::::::::
!
      module subroutine get_mlccsd_x2am_mlccsd(wf, x_ia_jb)
!!
!!
         implicit none
!
         class(mlccsd) :: wf
!
         real(dp), dimension(:,:) :: x_ia_jb
!
     end subroutine get_mlccsd_x2am_mlccsd
!
!
    module subroutine omega_mlccsd_a1_mlccsd(wf, x_ib_jc)
!! 
!!    Omega A1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!   
!!    Calculates the A1 term of omega for the active space, 
!!   
!!    A1: sum_bcj g_Abjc * u_Ij^bc,
!!
!!    where upper case letters indicate CCS space, and 
!!    lower case letters are the combined CC2/CCSD spaces.
!!  
!!    A1 is added to the projection vector (omega1) of
!!    the wavefunction object wf.
!! 
!!    u_ij^bc = 2*x_ij^bc - x_ij^cb 
!!
!!    Batching over A.
!!
!! 
      implicit none
!
      class(mlccsd)            :: wf
      real(dp), dimension(:,:) :: x_ib_jc
!
      end subroutine omega_mlccsd_a1_mlccsd
!
!
      module subroutine omega_mlccsd_b1_mlccsd(wf, x_ja_kb)
!! 
!!       Omega B1
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!   
!!       Calculates the B1 term of omega, 
!!   
!!       B1: - sum_bjk u_jk^ab*g_kbjI + sum_bj u_ij^ab F_jb,
!!
!!       with u_ij^ab = 2*x_ij^ab - x_ij^ba. 
!!
         implicit none
!
         class(mlccsd)            :: wf
         real(dp), dimension(:,:) :: x_ja_kb  
!
      end subroutine omega_mlccsd_b1_mlccsd
!
!
      module subroutine omega_mlccsd_a2_mlccsd(wf, x_IC_JD)
!
!        Omega A2 term: Omega A2 = sum_(cd)g_aC_bD * x_Ci_Dj
!
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, 10 Mar 2017
!
!        Structure: Batching over both a and b for A2.2.
!                   x^+_Ci_Dj = x_Ci_Dj + x_Di_Cj
!                   x^-_Ci_Dj = x_Ci_Dj - x_Di_Cj
!                   g^+_aC_bD = g_aC_bD + g_bC_aD 
!                   g^-_aC_bD = g_aC_bD - g_bC_aD 
! 
!                   omega_A2_ai_bj = 1/4*(g^+_aC_bD*x^+_Ci_Dj + g^-_aC_bD*x^-_Ci_Dj)
!                   omega_A2_aj_bi = 1/4*(g^+_aC_bD*x^+_Ci_Dj - g^-_aC_bD*x^-_Ci_Dj)
!
         implicit none
!
         class(mlccsd)  :: wf
         real(dp), dimension(:,:) :: x_IC_JD
!
      end subroutine omega_mlccsd_a2_mlccsd
!
!
      module subroutine omega_mlccsd_b2_mlccsd(wf, x_kc_ld)
!!
!!       Omega B2
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, 11 Mar 2017
!! 
!!       Omega B2 = sum_(kl) x_ak_bl*(g_kilj + sum_(cd) x_ci_dj * g_kc_ld)
!!
!!       Structure: g_kilj is constructed first and reordered as g_kl_ij. 
!!       Then the contraction over cd is performed, and the results added to g_kl_ij.
!!       x_ka_lb is then reordered as x_ab_kl and the contraction over kl is performed.
!!
         implicit none
!
         class(mlccsd) :: wf 
         real(dp), dimension(:,:) :: x_kc_ld 
!
      end subroutine omega_mlccsd_b2_mlccsd
!
!
      module subroutine omega_mlccsd_c2_mlccsd(wf, x_lc_kd)
!!
!!       Omega C2 
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2017
!!       
!!       Z_ai_bj = - sum_(ck) x^bc_kj*( g_ki,ac - 1/2 * sum_(dl) x^ad_li*g_kd,cl )
!!
!!       Omega_ai_bj = P_ij^ab (1/2 Z_ai_bj + Z_aj_bi) = 1/2 Z_ai_bj + 1/2 Z_bj_ai + Z_aj_bi+ Z_bi_aj
!!
!!       
         implicit none
!
         class(mlccsd)            :: wf
         real(dp), dimension(:,:) :: x_lc_kd
!
      end subroutine omega_mlccsd_c2_mlccsd
!
!
      module subroutine omega_mlccsd_d2_mlccsd(wf, x_KC_LD)
!!
!!       Omega D2 
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!       Calculates the D2 term,
!!
!!        D2: sum_ck u_jk^bc g_aikc 
!!          - 1/2 * sum_ck u_jk^bc g_acki 
!!          + 1/4 * sum_ck u_jk^bc sum_dl L_ldkc u_il^ad,
!!
!!       where 
!!
!!          u_jk^bc = 2 * t_jk^bc - t_kj^bc,
!!          L_ldkc  = 2 * g_ldkc  - g_lckd.
!!
!!       The first, second, and third terms are referred to as D2.1, D2.2, and D2.3, 
!!       and comes out ordered as (ai,bj). All terms are added to the omega vector of the 
!!       wavefunction object wf.
!
!        The routine adds the terms in the following order: D2.3, D2.1, D2.2
!
         implicit none 
!
         class(mlccsd) :: wf 
         real(dp), dimension(:,:) :: x_KC_LD
!
      end subroutine omega_mlccsd_d2_mlccsd
!
!
      module subroutine omega_mlccsd_e2_mlccsd(wf, x_kc_ld)
!!
!!     Omega E2
!!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!     Calculates the E2 term,
!!
!!      E2: sum_c t_ij^ac (F_bc - sum_dkl g_ldkc u_kl^bd) 
!!        - sum_k t_ik^ab (F_kj + sum_cdl g_ldkc u_lj^dc),
!!
!!     where
!!
!!        u_kl^bc = 2 * t_kl^bc - t_lk^bc.
!!
!!     The first term is referred to as the E2.1 term, and comes out ordered as (b,jai).
!!     The second term is referred to as the E2.2 term, and comes out ordered as (aib,j).
!!
!!     Both are permuted added to the projection vector element omega2(ai,bj) of
!!     the wavefunction object wf.
!!
       implicit none 
!
      class(mlccsd) :: wf 
      real(dp), dimension(:,:) :: x_kc_ld
!
   end subroutine omega_mlccsd_e2_mlccsd
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
      module subroutine calc_ampeqs_norm_mlccsd(wf, ampeqs_norm)
!
!        Calculate Amplitude Equations Norm (MLCCSD)
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!
         implicit none 
!
         class(mlccsd) :: wf 
!
         real(dp) :: ampeqs_norm 
!
      end subroutine calc_ampeqs_norm_mlccsd
!
!
   module subroutine new_amplitudes_mlccsd(wf)
!
!     New Amplitudes (MLCCSD)
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!
!     Directs the calculation of the quasi-Newton estimate Δ t_i, 
!     and t_i + Δ t_i, and calls the DIIS routine to save & get 
!     the amplitudes for the next iteration.
!
      implicit none 
!
      class(mlccsd) :: wf 
!
   end subroutine new_amplitudes_mlccsd
!
!
   module subroutine calc_quasi_Newton_doubles_mlccsd(wf,dt)
!
!     Calculate quasi-Newtoni doubles estimate (CCSD)
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!
!     Calculates the quasi-Newton estimate Δ t_i (doubbles part)
!     and places the contribution in the dt vector (of length n_parameters,
!     with singles first, then doubles, etc. if inherited)
!
      implicit none 
!
      class(mlccsd) :: wf 
!
      real(dp), dimension(wf%n_parameters, 1) :: dt
!
!
   end subroutine calc_quasi_Newton_doubles_mlccsd
!
!
   module subroutine initialize_ground_state_mlccsd(wf)
!!
!!    Initialize Ground State (CCSD)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!    Initializes the amplitudes and the projection vector for the ground
!!    state solver.
!!
      implicit none 
!
      class(mlccsd) :: wf
!
   end subroutine initialize_ground_state_mlccsd
!
!
   end interface
!
!
   interface
!
!     -::- Jacobi submodule interface -::-
!     ::::::::::::::::::::::::::::::::::::
!
   module subroutine jacobian_mlccsd_transformation_mlccsd(wf, c_a_i, c_aibj)
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
      class(mlccsd) :: wf 
!
!      Incoming vector c 
!
      real(dp), dimension(wf%n_v, wf%n_o) :: c_a_i  ! c_ai 
      real(dp), dimension(wf%n_x2am, 1)   :: c_aibj ! c_aibj  
!
      end subroutine jacobian_mlccsd_transformation_mlccsd
!
!
      module subroutine jacobian_mlccsd_b2_mlccsd(wf, rho_ai_bj, c_a_i)
!!
!!    Jacobian CCSD B2 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!    rho_ai_bj^B2 = - sum_kc (F_kc t_ij^ac c_bk + F_kc t_ik^ab c_cj)
!!
!!    The term is added as rho_ai_bj(ai,bj) = rho_ai_bj(ai,bj) + rho_ai_bj^B2,
!!    where c_a_i(a,i) = c_ai above.
!!
      implicit none 
!
      class(mlccsd) :: wf 
!
      real(dp), dimension(wf%n_v, wf%n_o)                       :: c_a_i
      real(dp), dimension(:, :) :: rho_ai_bj
!
!
      end subroutine jacobian_mlccsd_b2_mlccsd
!
!
      module subroutine jacobian_mlccsd_c2_mlccsd(wf, rho_ai_bj, c_a_i)
!!
!!    Jacobian CCSD C2 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2017
!!
!!    rho_ai_bj^C2 = sum_kcl (g_ljkc * t_ki^ac * c_bl) + (g_ljkc * t_li^bc * c_ak) 
!!                         + (g_ljkc * t_lk^ba * c_ci) 
!!                         - (L_ljkc * t_ik^ac * c_bl)- (L_ljkc * t_il^ab * c_ck)
!!                
!!
      implicit none 
!
      class(mlccsd) :: wf 
!
      real(dp), dimension(wf%n_v, wf%n_o) :: c_a_i
      real(dp), dimension(:,:)            :: rho_ai_bj
!
      end subroutine jacobian_mlccsd_c2_mlccsd
!
!
      module subroutine jacobian_mlccsd_d2_mlccsd(wf, rho_ai_bj, c_a_i)
!!
!!    Jacobian CCSD D2 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!    rho_ai_bj^D2 = - sum_kcd g_kcbd (x_ij^cd c_ak + x_kj^ad c_ci + x_ik^ca c_dj)
!!                       + sum_kcd L_kcbd (x_ik^ac c_dj + x_ij^ad c_ck)
!!
!!    Note: the code is structured so that we batch over the index b,
!!          where the integrals are made as g_kc_db = g_kcbd and held
!!          in some ordering or other throughout a given batch (i.e.,
!!          all five terms are constructed gradually in the batches).
!!
!!    The term is added as rho_ai_bj(ai,bj) = rho_ai_bj(ai,bj) + rho_ai_bj^D2,
!!    where c_a_i(a,i) = c_ai above.
!!
      implicit none 
!
      class(mlccsd) :: wf 
!
      real(dp), dimension(wf%n_v, wf%n_o)                       :: c_a_i
      real(dp), dimension(:,:) :: rho_ai_bj
!
      end subroutine jacobian_mlccsd_d2_mlccsd
!
!
      module subroutine jacobian_mlccsd_e2_mlccsd(wf, rho_ai_bj, c_ai_ck)
!!
!!    Jacobian MLCCSD E2 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!    rho_ai_bj^E2 = 2 sum_dlck x_bj,DL * L_KC,LD * c_ai,CK 
!!
      implicit none 
!
      class(mlccsd) :: wf 
!
      real(dp), dimension(:,:) :: rho_ai_bj
      real(dp), dimension(:,:) :: c_ai_ck
!
      end subroutine jacobian_mlccsd_e2_mlccsd
!
!
      module subroutine jacobian_mlccsd_f2_mlccsd(wf, rho_ai_bj, c_ai_bj)
!!
!!       Jacobian MLCCSD F2 
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!       rho_ai_bj^F2 =   - sum_(CKDL) x_ai,CK * L_KC,LD * c_bL,Dj 
!!                        - sum_(CKDL) x_ai,Dj * L_KC,LD * c_bL,CK
!!                        - sum_(CKDL) x_ai_bL * L_KC,LD * c_CK,Dj
!!
!!
         implicit none 
!
         class(mlccsd) :: wf 
!
         real(dp), dimension(:,:) :: rho_ai_bj
         real(dp), dimension(:,:) :: c_ai_bj
!
      end subroutine jacobian_mlccsd_f2_mlccsd
!
!
     module subroutine jacobian_mlccsd_g2_mlccsd(wf, rho_ai_bj, c_ai_bj)
!!
!!       Jacobian MLCCSD G2 
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!       rho_ai_bj^G2 =  - sum_ckdl x_bL,Dj * L_KC,LD * c_ai,CK 
!!                       - sum_ckdl x_CK_bL * L_KC,LD * c_ai,Dj 
!!                       - sum_ckld x_CK,Dj * L_KC,LD * c_ai,bL 
!!
!!
         implicit none 
!
         class(mlccsd) :: wf 
!
         real(dp), dimension(:,:) :: rho_ai_bj
         real(dp), dimension(:,:) :: c_ai_bj
!
      end subroutine jacobian_mlccsd_g2_mlccsd
!
!
     module subroutine jacobian_mlccsd_h2_mlccsd(wf, rho_ai_bj, c_ai_bj)
!!
!!       Jacobian MLCCSD H2 
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!       rho_ai_bj^H2 =  sum_CKDL x_Ci,aK * g_KC,LD * c_bL,Dj 
!!                     + sum_CKDL x_Cj,aL * g_KC,LD * c_bK,Di
!!                
         implicit none 
!
         class(mlccsd) :: wf 
!
         real(dp), dimension(:,:) :: rho_ai_bj
         real(dp), dimension(:,:) :: c_ai_bj
!
   end subroutine jacobian_mlccsd_h2_mlccsd
!
!
      module subroutine jacobian_mlccsd_i2_mlccsd(wf, rho_ai_bj, c_ai_bj)
!!
!!       Jacobian MLCCSD I2 
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!       rho_ai_bj^I2 =  sum_C F_bC * c_ai,Cj - sum_K F_jK * c_ai,bK
!!                     + sum_ck L_bj,KC * c_ai,CK 
!!                     - sum_ck ( g_KC,bj * c_aK,Ci + g_Ki,bC * c_aK,Cj ) 
!!                
!!       Batch over c to construct  g_ki_bC
!! 
         implicit none 
!
         class(mlccsd) :: wf 
!
         real(dp), dimension(:,:) :: rho_ai_bj
         real(dp), dimension(:,:) :: c_ai_bj
!
   end subroutine jacobian_mlccsd_i2_mlccsd
!
!
      module subroutine jacobian_mlccsd_j2_mlccsd(wf, rho_ab_ij, c_ab_ij)
!!
!!       Jacobian CCSD J2 
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!       rho_ab_ij^J2 =    sum_ckld t_ci,dj * g_kc,ld * c_ak,bl 
!!                       + sum_ckdl t_ak,bl * g_kc,ld * c_ci,dj
!!             
         implicit none 
!
         class(mlccsd) :: wf 
!
         real(dp), dimension(:,:) :: rho_ab_ij
         real(dp), dimension(:,:) :: c_ab_ij
!
      end subroutine jacobian_mlccsd_j2_mlccsd
!
!
      module subroutine jacobian_mlccsd_k2_mlccsd(wf, rho_ab_ij, c_ab_ij)
!!
!!       Jacobian MLCCSD K2 
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!       rho_ab_ij^K2 =    sum_kl g_Ki,Lj * c_aK,bL 
!!                       + sum_cd g_aC,bD * c_Ci,Dj
!! 
!!       For the last term we batch over a and b and 
!!       add each batch to rho_ai_bj 
!!               
         implicit none 
!
         class(mlccsd) :: wf 
!
         real(dp), dimension(:,:) :: rho_ab_ij
         real(dp), dimension(:,:) :: c_ab_ij
!
      end subroutine jacobian_mlccsd_k2_mlccsd
!
!
!
   end interface
!
!
   interface
!
!
!     -::- Excited state submodule interface -::-
!     :::::::::::::::::::::::::::::::::::::::::::
!
      module subroutine transform_trial_vectors_mlccsd(wf, first_trial, last_trial)
!!
!!    Transformation of Trial Vectors (MLCC2)
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
      class(mlccsd) :: wf
!
      integer(i15), intent(in) :: first_trial, last_trial ! Which trial_vectors we are to transform
!
      end subroutine transform_trial_vectors_mlccsd
!
!
      module subroutine initialize_excited_states_mlccsd(wf)
!!
!!       Initialize excited states
!!       Written by Sarai D. Folkestad, Aug 2017
!!
!!       Calculates and sets n_s2am, and updates n_parameters
!!       for excited state calculation
!!
         implicit none 
!       
         class(mlccsd) :: wf
      end subroutine initialize_excited_states_mlccsd
!
!
      module subroutine print_excitation_vector_mlccsd(wf, vec, unit_id)
!!
!!
!!
         implicit none
!  
         class(mlccsd) :: wf
!
         real(dp), dimension(wf%n_parameters, 1) :: vec
!
         integer(i15) :: unit_id     
!
      END subroutine print_excitation_vector_mlccsd
!
!
   end interface
contains
!
!
  subroutine init_mlccsd(wf)
!!
!!
      implicit none
!
      class(mlccsd) :: wf
!
      integer(i15) :: i, j, active_space
!
      integer(i15) :: unit_input = -1
!
!
      write(unit_output,'(//t3,a/)')    ':: Initialization' 
!
!
!     Set model name 
!
      wf%name = 'MLCCSD'
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
      if(wf%mlcc_settings%CC3) then
         write(unit_output,*)'WARNING: CC3 active spaces not available for MLCCSD'
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
!     Orbital partitioning - only if we have CCS/CC2 region
!
      if (wf%mlcc_settings%CC2) then
!
            call wf%construct_orbital_coef_CC2_CCS
            call wf%construct_orbital_energy_CC2_CCS
!
      endif
!
      if (wf%mlcc_settings%CCS .or. wf%mlcc_settings%CC2) then
!
         call wf%orbital_partitioning
!
      else
!
         write(unit_output,'(/t3,a)')'Full CCSD requested, orbital partitioning skipped.'
         flush(unit_output)
!
!        Do full space CCSD calculation
!
         wf%n_CCSD_o = wf%n_o
         wf%n_CCSD_v = wf%n_v
!
         wf%first_CCSD_o = 1
         wf%first_CCSD_v = 1
!
         wf%n_CC2_o = 0
         wf%n_CC2_v = 0
!
      endif
!
!     Initialize amplitudes and associated attributes
!
!     Set number of amplitudes for CCSD active space
!
      wf%n_t2am = ((wf%n_CCSD_v)*(wf%n_CCSD_o))*((wf%n_CCSD_v )*(wf%n_CCSD_o)+1)/2

!
      call wf%initialize_amplitudes
      call wf%initialize_omega
!
!     Set the number of parameters in the wavefunction
!     (that are solved for in the ground and excited state solvers) 
!
      wf%n_parameters = wf%n_t1am + wf%n_t2am
!
!     Read Cholesky AO integrals and transform to MO basis
!
      if (wf%mlcc_settings%CC2) then
!
         call wf%read_transform_cholesky_for_CC2_amplitude
         call wf%construct_MO_transformation_matrix
         call wf%destruct_orbital_coef_CC2_CCS
!
      endif
!
      call wf%read_transform_cholesky
!
!     Initialize fock matrix
!
      call wf%initialize_fock_matrix
!
   end subroutine init_mlccsd
!
!
   subroutine initialize_amplitudes_mlccsd(wf)
!!
!!     Initialize Amplitudes (MLCCSD)
!!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!     Allocates the amplitudes, sets them to zero, and calculates
!!     the number of amplitudes.
!!
      implicit none 
!
      class(mlccsd) :: wf
!
!     Allocate the doubles amplitudes and set to zero
!
      wf%n_t1am = (wf%n_o)*(wf%n_v)
      if (.not. allocated(wf%t1am)) call allocator(wf%t1am, wf%n_v, wf%n_o)
      wf%t1am = zero
!
      wf%n_t2am = ((wf%n_CCSD_v)*(wf%n_CCSD_o))*((wf%n_CCSD_v )*(wf%n_CCSD_o)+1)/2
!
      if (.not. allocated(wf%t2am)) call allocator(wf%t2am, wf%n_t2am, 1)
      wf%t2am = zero
!
   end subroutine initialize_amplitudes_mlccsd
!
   subroutine construct_perturbative_doubles_mlccsd(wf)
!!
!!    Construct Perturbative Doubles (CCSD)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!    Sets the doubles amplitudes (t2am) to its MP2 estimate. This is
!!    the initial guess used in the solver for the ground state amplitude 
!!    equations.
!!
      implicit none 
!
      class(mlccsd) :: wf
!
      real(dp), dimension(:,:), allocatable :: L_ai_J
      real(dp), dimension(:,:), allocatable :: g_ai_bj
!
      integer(i15) :: i = 0, j = 0, a = 0, b = 0
      integer(i15) :: ai = 0, bj = 0, ia = 0, jb = 0, aibj = 0 
!
!     Active space variables
!
      integer(i15) :: n_active_o = 0, n_active_v = 0
      integer(i15) :: first_active_o ! first active occupied index 
      integer(i15) :: first_active_v ! first active virtual index
      integer(i15) :: last_active_o ! first active occupied index 
      integer(i15) :: last_active_v ! first active virtual index
!
!     Calculate first/last indeces
! 
      call wf%get_CCSD_active_indices(first_active_o, first_active_v)
      call wf%get_CCSD_n_active(n_active_o, n_active_v)
!
      last_active_o = first_active_o + n_active_o - 1
      last_active_v = first_active_v + n_active_v - 1 
!
!     Allocate L_ia_J and g_ia_jb
!
      call allocator(L_ai_J, (n_active_o)*(n_active_v), wf%n_J)
      call allocator(g_ai_bj, (n_active_o)*(n_active_v), (n_active_o)*(n_active_v))
!
      L_ai_J = zero
      g_ai_bj = zero
!
!     Get the Cholesky IA vector 
!
      call wf%get_cholesky_ai(L_ai_J, first_active_v, last_active_v, first_active_o, last_active_o)
!
!     Calculate g_ia_jb = g_iajb
!
      call dgemm('N','T',                 &
                  n_active_o*n_active_v,  & 
                  n_active_o*n_active_v,  &
                  wf%n_J,                 &
                  one,                    &
                  L_ai_J,                 &
                  n_active_o*n_active_v,  &
                  L_ai_J,                 &
                  n_active_o*n_active_v,  &
                  zero,                   &
                  g_ai_bj,                &
                  n_active_o*n_active_v)
!
!     Set the doubles amplitudes
!
      do i = 1, n_active_o
         do a = 1, n_active_v
!
            ai = index_two(a, i, n_active_v)
!
            do j = 1, n_active_o
               do b = 1, n_active_v
!     
                  bj = index_two(b, j, n_active_v)
!
!                 Set the doubles amplitudes
!
                  if (ai .le. bj) then ! To avoid setting the same element twice
!
                     aibj = index_packed(ai,bj)
!
                     wf%t2am(aibj, 1) = - g_ai_bj(ai,bj)/(wf%fock_diagonal(wf%n_o + a + first_active_v - 1, 1) + &
                                                            wf%fock_diagonal(wf%n_o + b + first_active_v - 1, 1) - &
                                                            wf%fock_diagonal(i + first_active_o - 1, 1) - &
                                                            wf%fock_diagonal(j + first_active_o - 1, 1))
!
                  endif
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocations
!
      call deallocator(L_ai_J, (n_active_o)*(n_active_v), wf%n_J)
      call deallocator(g_ai_bj, (n_active_o)*(n_active_v), (n_active_o)*(n_active_v)) 
!
   end subroutine construct_perturbative_doubles_mlccsd
!
   subroutine get_CCSD_active_indices_mlccsd(wf, first_o, first_v)
!!
!!    Get CC2 active indices,
!!    Written by Sarai D. Folkestad, June 2017
!!
!!    Returns the first active occupied and virtual indices 
!!    of the active space.
!!
   implicit none
!
   class(mlccsd) :: wf
   integer(i15) :: first_o
   integer(i15) :: first_v
! 
   first_o = wf%first_CCSD_o
   first_v = wf%first_CCSD_v
!
   end subroutine get_CCSD_active_indices_mlccsd
!
!
   subroutine get_CC2_active_indices_mlccsd(wf, first_o, first_v)
!!
!!    Get CC2 active indices,
!!    Written by Sarai D. Folkestad, June 2017
!!
!!    Returns the first active occupied and virtual indices 
!!    of the active space.
!!
      implicit none
!
      class(mlccsd) :: wf
      integer(i15) :: first_o
      integer(i15) :: first_v
!
      first_o = wf%first_CCSD_o
      first_v = wf%first_CCSD_v
!
   end subroutine get_CC2_active_indices_mlccsd
!
!
   subroutine get_CC2_n_active_mlccsd(wf, n_active_o, n_active_v)
!!
!!    Get CC2 active indices,
!!    Written by Sarai D. Folkestad, June 2017
!!
!!    Returns the first active occupied and virtual indices 
!!    of the active space.
!!
      implicit none
!
      class(mlccsd) :: wf
      integer(i15) :: n_active_o
      integer(i15) :: n_active_v
!
      n_active_o = wf%n_CC2_o + wf%n_CCSD_o
      n_active_v = wf%n_CC2_v + wf%n_CCSD_v
!
   end subroutine get_CC2_n_active_mlccsd
!
!
   subroutine get_CCSD_n_active_mlccsd(wf, n_active_o, n_active_v)
!!
!!    Get CC2 active indices,
!!    Written by Sarai D. Folkestad, June 2017
!!
!!    Returns the first active occupied and virtual indices 
!!    of the active space.
!!
      implicit none
!
      class(mlccsd) :: wf
      integer(i15) :: n_active_o
      integer(i15) :: n_active_v
!
      n_active_o = wf%n_CCSD_o
      n_active_v = wf%n_CCSD_v
!
   end subroutine get_CCSD_n_active_mlccsd
!
!
   subroutine calc_energy_mlccsd(wf)
!!
!!     Calculate Energy (MLCCSD)
!!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!     Calculates the MLCCSD energy for the wavefunction's current amplitudes.
!!
!!       E = E_scf + sum_AIBJ (t^A_I * t^B_J)*L_IAJB
!!                 + sum aibj (s^ab_ij + t^ab_ij)L_ia_jb
!!
!!     where lower case letters indicate indices restricted to the CC2+CCSD regions.
!!
!!     s^ab_ij is zero in the CCSD region, and t^ab_ij is zero in the CC2 region.
!!
!!
      implicit none 
!
      class(mlccsd) :: wf 
!
      real(dp), dimension(:,:), allocatable :: L_ia_J  ! L_ia^J
      real(dp), dimension(:,:), allocatable :: g_ia_jb ! g_iajb
      real(dp), dimension(:,:), allocatable :: x_ia_jb ! s_ij^ab
!
      integer(i15) :: a = 0, i = 0, b = 0, j = 0, ai = 0
      integer(i15) :: bj = 0, aibj = 0, ia = 0, jb = 0, ib = 0, ja = 0
      integer(i15) :: ia_full = 0, jb_full = 0
!
!     Active space variables
!
      integer(i15) :: n_active_o = 0, n_active_v = 0
      integer(i15) :: first_active_o ! first active occupied index 
      integer(i15) :: first_active_v ! first active virtual index
      integer(i15) :: last_active_o ! first active occupied index 
      integer(i15) :: last_active_v ! first active virtual index
!
!     Allocate the Cholesky vector L_ia_J = L_ia^J and set to zero 
!
      call allocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
      L_ia_J = zero
!
!     Get the Cholesky vector L_ia_J 
!
      call wf%get_cholesky_ia(L_ia_J)
!
!     Allocate g_ia_jb = g_iajb and set it to zero
!
      call allocator(g_ia_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      g_ia_jb = zero
!
!     Calculate the integrals g_ia_jb from the Cholesky vector L_ia_J 
!
      call dgemm('N','T',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  wf%n_J,            &
                  one,               &
                  L_ia_J,            &
                  (wf%n_o)*(wf%n_v), &
                  L_ia_J,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  g_ia_jb,           &
                  (wf%n_o)*(wf%n_v))
!
!     Deallocate the Cholesky vector L_ia_J 
!
      call deallocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!     Set the initial value of the energy 
!
      wf%energy = wf%scf_energy
!
!     t1 amplitude contribution
!
!
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            ia = index_two(i, a, wf%n_o)
!
            do j = 1, wf%n_o
!
               ja = index_two(j, a, wf%n_o)
!
               do b = 1, wf%n_v
!
                  jb = index_two(j, b, wf%n_o)
                  ib = index_two(i, b, wf%n_o)
!
!                 Add the correlation energy from the single amplitudes
!
                  wf%energy = wf%energy +                       & 
                                 (wf%t1am(a,i))*(wf%t1am(b,j))* &
                                 (two*g_ia_jb(ia,jb) - g_ia_jb(ib,ja))
!
               enddo
            enddo
         enddo
      enddo
!
!     t2 amplitude contribution to the energy
!
!     Calculate first/last indeces
! 
      call wf%get_CCSD_active_indices(first_active_o, first_active_v)
      call wf%get_CC2_n_active(n_active_o, n_active_v)
!
      call allocator(x_ia_jb, n_active_v*n_active_o, n_active_v*n_active_o)
      call wf%get_mlccsd_x2am(x_ia_jb)
!
      do i = 1, n_active_o
         do a = 1, n_active_v
!
            ia = index_two(i, a, n_active_o)
            ia_full = index_two(i + first_active_o - 1, a + first_active_v - 1, wf%n_o)
!
            do j = 1, n_active_o
!
               ja = index_two(j + first_active_o - 1, a + first_active_v - 1, wf%n_o)
!
               do b = 1, n_active_v
!
                  jb_full = index_two(j + first_active_o - 1, b + first_active_v - 1, wf%n_o)
                  jb = index_two(j , b , n_active_o)
                  ib = index_two(i + first_active_o - 1, b + first_active_v - 1, wf%n_o)
!

!
!                 Add the correlation energy from the double amplitudes
!
                  wf%energy = wf%energy +              & 
                                 (x_ia_jb(ia, jb))*    &
                                 (two*g_ia_jb(ia_full, jb_full) - g_ia_jb(ib,ja))
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate g_ia_jb
!
      call deallocator(g_ia_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call deallocator(x_ia_jb, n_active_v*n_active_o, n_active_v*n_active_o)
!
   end subroutine calc_energy_mlccsd
!
!
   subroutine construct_orbital_coef_CC2_CCS_mlccsd(wf)
!!
!!
      implicit none
!
      class(mlccsd) :: wf
!
      if (.not. allocated(wf%mo_coef_cc2_ccs)) call allocator(wf%mo_coef_cc2_ccs, wf%n_ao, wf%n_mo)
      wf%mo_coef_cc2_ccs = zero
!
   end subroutine construct_orbital_coef_CC2_CCS_mlccsd
!
!
   subroutine destruct_orbital_coef_CC2_CCS_mlccsd(wf)
!!
!!
      implicit none
!
      class(mlccsd) :: wf
!
      if (allocated(wf%mo_coef_cc2_ccs)) call deallocator(wf%mo_coef_cc2_ccs, wf%n_ao, wf%n_mo)
!
   end subroutine destruct_orbital_coef_CC2_CCS_mlccsd
!
!
   subroutine construct_orbital_energy_CC2_CCS_mlccsd(wf)
!!
!!
      implicit none
!
      class(mlccsd) :: wf
!
      if (.not. allocated(wf%fock_diagonal_cc2_ccs)) call allocator(wf%fock_diagonal_cc2_ccs, wf%n_mo, 1)
      wf%mo_coef_cc2_ccs = zero
!
   end subroutine construct_orbital_energy_CC2_CCS_mlccsd
!
   subroutine destruct_orbital_energy_CC2_CCS_mlccsd(wf)
!!
!!
      implicit none
!
      class(mlccsd) :: wf
!
      if (allocated(wf%fock_diagonal_cc2_ccs)) call deallocator(wf%fock_diagonal_cc2_ccs, wf%n_mo, 1)
!
   end subroutine destruct_orbital_energy_CC2_CCS_mlccsd
!
!
   subroutine destruct_double_amplitudes_mlccsd(wf)
!!
!!
      implicit none
!
      class(mlccsd) :: wf
!
      if (allocated(wf%t2am)) call deallocator(wf%t2am, wf%n_t2am, 1)
!
   end subroutine destruct_double_amplitudes_mlccsd
!
   subroutine destruct_amplitudes_mlccsd(wf)
!!
!!
      implicit none
!
      class(mlccsd) :: wf
!
      call wf%destruct_double_amplitudes
      call destruct_amplitudes_ccs(wf)
!
   end subroutine destruct_amplitudes_mlccsd
!
!
   subroutine save_amplitudes_mlccsd(wf)
!!
!!    Save Amplitudes (CCSD)
!!    Written by Sarai D. Folkestad and Eirik F. Kjøsntad, May 2017
!!
!!    Store the amplitudes to disk (T1AM, T2AM)
!!
      implicit none 
!
      class(mlccsd) :: wf
!
      integer(i15) :: unit_x1am = -1
      integer(i15) :: unit_x2am = -1
      integer(i15) :: unit_t2am = -1
      integer(i15) :: n_active_o, n_active_v
!
      integer(i15) ::  i = 0, a = 0, b = 0, j = 0
      integer(i15) ::  ia = 0, ai = 0, jb = 0, bj = 0
      integer(i15) ::  aibj = 0
!
      real(dp), dimension(:,:), allocatable :: x_ia_jb
!
!     Open amplitude files
!
      call generate_unit_identifier(unit_x1am)
      call generate_unit_identifier(unit_x2am)
      call generate_unit_identifier(unit_t2am)
!
      open(unit_x1am, file='t1am', status='unknown', form='unformatted')
      open(unit_t2am, file='t2am', status='unknown', form='unformatted')
      open(unit_x2am, file='x2am', status='unknown', form='unformatted')
!
      rewind(unit_x1am)
      rewind(unit_x2am)
      rewind(unit_t2am)
!
!     Write amplitudes to files
!
      write(unit_x1am) wf%t1am
      write(unit_t2am) wf%t2am
!
      call wf%get_CC2_n_active(n_active_o, n_active_v)
      wf%n_x2am = (n_active_v)*(n_active_o)*((n_active_v)*(n_active_o)+1)/2
!
      call allocator(x_ia_jb, (n_active_v)*(n_active_o), (n_active_v)*(n_active_o)) 
      call wf%get_mlccsd_x2am(x_ia_jb)
!
!     Reorder and pack in
!
      call allocator(wf%x2am, wf%n_x2am, 1)
!
      do i = 1, n_active_o
         do a = 1, n_active_v
            ai = index_two(a, i, n_active_v)
            ia = index_two(i, a, n_active_o)
            do j = 1, n_active_o
               do b = 1, n_active_v
!
                  bj = index_two(b, j, n_active_v)
                  jb = index_two(j, b, n_active_o)
!
                  aibj = index_packed(ai, bj)
!
                  wf%x2am(aibj, 1) = x_ia_jb(ia, jb)
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(x_ia_jb, (n_active_v)*(n_active_o), (n_active_v)*(n_active_o)) 
!
      write(unit_x2am) wf%x2am
      call deallocator(wf%x2am, wf%n_x2am, 1)
!
!     Close amplitude files
!
      close(unit_x1am)
      close(unit_x2am)
      close(unit_t2am)
!
   end subroutine save_amplitudes_mlccsd
!
   subroutine read_amplitudes_mlccsd(wf)
!!
!!    Read Amplitudes (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjøsntad, May 2017
!!
!!    Reads the amplitudes from disk (T1AM)
!!
      implicit none 
!
      class(mlccsd) :: wf
!
      call wf%read_single_amplitudes
      call wf%read_mlccsd_double_amplitudes
!
   end subroutine read_amplitudes_mlccsd
!
   subroutine read_mlccsd_double_amplitudes_mlccsd(wf)
!!
!!    Read Amplitudes (CCSD)
!!    Written by Sarai D. Folkestad and Eirik F. Kjøsntad, May 2017
!!
!!    Reads the amplitudes from disk (T1AM, T2AM)
!!
      implicit none 
!
      class(mlccsd) :: wf
!
      integer(i15) :: unit_x2am = -1 
      integer(i15) :: n_active_v, n_active_o 
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
         call wf%get_CC2_n_active(n_active_o, n_active_v)
         wf%n_x2am = (n_active_v)*(n_active_o)*((n_active_v)*(n_active_o)+1)/2
!
         if (.not. allocated(wf%x2am)) call allocator(wf%x2am, wf%n_x2am, 1)
!
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
   end subroutine read_mlccsd_double_amplitudes_mlccsd
!
!
   subroutine read_double_amplitudes_mlccsd(wf)
!!
!!    Read Amplitudes (CCSD)
!!    Written by Sarai D. Folkestad and Eirik F. Kjøsntad, May 2017
!!
!!    Reads the amplitudes from disk (T1AM, T2AM)
!!
      implicit none 
!
      class(mlccsd) :: wf
!
      integer(i15) :: unit_t2am = -1 
      integer(i15) :: n_active_v, n_active_o 
!
      logical :: file_exists = .false.
!
!     Check to see whether file exists
!
      inquire(file='t2am',exist=file_exists)
!
      if (file_exists) then 
!
!        Open amplitude files if they exist
!
         call generate_unit_identifier(unit_t2am)
!
         open(unit_t2am, file='t2am', status='unknown', form='unformatted')
!
         rewind(unit_t2am)
!
!        Read from file & close
!
         call wf%get_CCSD_n_active(n_active_o, n_active_v)
         wf%n_t2am = (n_active_v)*(n_active_o)*((n_active_v)*(n_active_o)+1)/2
!
         if (.not. allocated(wf%t2am)) call allocator(wf%t2am, wf%n_t2am, 1)
!
         read(unit_t2am) wf%t2am
!
         close(unit_t2am)
!
      else
!
         write(unit_output,'(t3,a)') 'Error: amplitude files do not exist.'
         stop
!
      endif
!
   end subroutine read_double_amplitudes_mlccsd
end module mlccsd_class