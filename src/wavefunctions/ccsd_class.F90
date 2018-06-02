module ccsd_class
!
!!
!!
!!                Coupled cluster singles and doubles (CCSD) class module
!!              Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!
!!    This module contains the definition of the coupled cluster singles
!!    and doubles (CCSD) wavefunction class. It is structured into four sections:
!!
!!       1. Modules used by the class:
!!
!!             Basic utilities and the ancestor class (CCS)
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
!!                - Omega
!!                - Excited state
!!                - Jacobian (right transformation)
!!                - Jacobian transpose (left transformation)
!!                - Ionized state
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
!
!
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
!
!  Ancestor class module (CCS)
!
   use ccs_class
!
   implicit none
!
!
!  :::::::::::::::::::::::::::::::::::::::::
!  -::- 2. Definition of the CCSD class -::-
!  :::::::::::::::::::::::::::::::::::::::::
!
!
   type, extends(ccs) :: ccsd
!
!     Cluster amplitudes
!
      integer(i15) :: n_t2am = 0 ! Number of doubles excitation amplitudes
!
      real(dp), dimension(:,:), allocatable :: t2am ! Doubles amplitude vector, t_ij^ab
!
!     The omega, or projection, vector < mu | exp(-T) H exp(T) | R >
!
      real(dp), dimension(:,:), allocatable :: omega2 ! Doubles projection vector
!
   contains
!
!
!     -::- Initialization routine -::-
!     ---------------------------------
!
      procedure :: init => init_ccsd
!
!
!     -::- Ground state submodule routine pointers -::-
!     -------------------------------------------------
!
      procedure :: calc_ampeqs_norm          => calc_ampeqs_norm_ccsd
      procedure :: new_amplitudes            => new_amplitudes_ccsd
      procedure :: calc_quasi_Newton_doubles => calc_quasi_Newton_doubles_ccsd
!
      procedure :: ground_state_preparations => ground_state_preparations_ccsd
      procedure :: summary_ground_state_info => summary_ground_state_info_ccsd
!
!
!     -::- Omega submodule routine pointers -::-
!     ------------------------------------------
!
      procedure :: destruct_omega   => destruct_omega_ccsd
      procedure :: initialize_omega => initialize_omega_ccsd
!
      procedure :: construct_omega  => construct_omega_ccsd
!
      procedure :: omega_ccsd_a1    => omega_ccsd_a1_ccsd
      procedure :: omega_ccsd_b1    => omega_ccsd_b1_ccsd
      procedure :: omega_ccsd_c1    => omega_ccsd_c1_ccsd
!
      procedure :: omega_ccsd_a2    => omega_ccsd_a2_ccsd
      procedure :: omega_ccsd_b2    => omega_ccsd_b2_ccsd
      procedure :: omega_ccsd_c2    => omega_ccsd_c2_ccsd
      procedure :: omega_ccsd_d2    => omega_ccsd_d2_ccsd
      procedure :: omega_ccsd_e2    => omega_ccsd_e2_ccsd
!
!
!     -::- Excited state submodule routine pointers -::-
!     --------------------------------------------------
!
      procedure :: calculate_orbital_differences    => calculate_orbital_differences_ccsd
      procedure :: transform_trial_vectors          => transform_trial_vectors_ccsd
      procedure :: print_excitation_vector          => print_excitation_vector_ccsd
      procedure :: summary_excited_state_info       => summary_excited_state_info_ccsd
      procedure :: excited_state_preparations       => excited_state_preparations_ccsd
!
!
!     -::- Jacobian submodule routine pointers -::-
!     ---------------------------------------------
!
      procedure :: jacobian_ccsd_transformation => jacobian_ccsd_transformation_ccsd
!
      procedure :: jacobian_ccsd_a1 => jacobian_ccsd_a1_ccsd
      procedure :: jacobian_ccsd_b1 => jacobian_ccsd_b1_ccsd
      procedure :: jacobian_ccsd_c1 => jacobian_ccsd_c1_ccsd
      procedure :: jacobian_ccsd_d1 => jacobian_ccsd_d1_ccsd
!
      procedure :: jacobian_ccsd_a2 => jacobian_ccsd_a2_ccsd
      procedure :: jacobian_ccsd_b2 => jacobian_ccsd_b2_ccsd
      procedure :: jacobian_ccsd_c2 => jacobian_ccsd_c2_ccsd
      procedure :: jacobian_ccsd_d2 => jacobian_ccsd_d2_ccsd
      procedure :: jacobian_ccsd_e2 => jacobian_ccsd_e2_ccsd
      procedure :: jacobian_ccsd_f2 => jacobian_ccsd_f2_ccsd
      procedure :: jacobian_ccsd_g2 => jacobian_ccsd_g2_ccsd
      procedure :: jacobian_ccsd_h2 => jacobian_ccsd_h2_ccsd
      procedure :: jacobian_ccsd_i2 => jacobian_ccsd_i2_ccsd
      procedure :: jacobian_ccsd_j2 => jacobian_ccsd_j2_ccsd
      procedure :: jacobian_ccsd_k2 => jacobian_ccsd_k2_ccsd
!
      procedure :: jacobi_test      => jacobi_test_ccsd      ! A debug routine
!
!
!     -::- Jacobian transpose submodule routine pointers -::-
!     -------------------------------------------------------
!
      procedure :: jacobian_transpose_ccsd_transformation => jacobian_transpose_ccsd_transformation_ccsd
!
      procedure :: jacobian_transpose_ccsd_a1 => jacobian_transpose_ccsd_a1_ccsd
      procedure :: jacobian_transpose_ccsd_b1 => jacobian_transpose_ccsd_b1_ccsd
      procedure :: jacobian_transpose_ccsd_c1 => jacobian_transpose_ccsd_c1_ccsd
      procedure :: jacobian_transpose_ccsd_d1 => jacobian_transpose_ccsd_d1_ccsd
      procedure :: jacobian_transpose_ccsd_e1 => jacobian_transpose_ccsd_e1_ccsd
      procedure :: jacobian_transpose_ccsd_f1 => jacobian_transpose_ccsd_f1_ccsd
      procedure :: jacobian_transpose_ccsd_g1 => jacobian_transpose_ccsd_g1_ccsd
!
      procedure :: jacobian_transpose_ccsd_a2 => jacobian_transpose_ccsd_a2_ccsd
      procedure :: jacobian_transpose_ccsd_b2 => jacobian_transpose_ccsd_b2_ccsd
      procedure :: jacobian_transpose_ccsd_c2 => jacobian_transpose_ccsd_c2_ccsd
      procedure :: jacobian_transpose_ccsd_d2 => jacobian_transpose_ccsd_d2_ccsd
      procedure :: jacobian_transpose_ccsd_e2 => jacobian_transpose_ccsd_e2_ccsd
      procedure :: jacobian_transpose_ccsd_f2 => jacobian_transpose_ccsd_f2_ccsd
      procedure :: jacobian_transpose_ccsd_g2 => jacobian_transpose_ccsd_g2_ccsd
      procedure :: jacobian_transpose_ccsd_h2 => jacobian_transpose_ccsd_h2_ccsd
      procedure :: jacobian_transpose_ccsd_i2 => jacobian_transpose_ccsd_i2_ccsd
!
!
!     -::- Ionized state submodule routine pointers -::-
!     --------------------------------------------------
!
      procedure :: ionization_residual_projection => ionization_residual_projection_ccsd
      procedure :: ionization_rho_aibj_projection => ionization_rho_aibj_projection_ccsd
!
!
!     -::- CVS submodule routine pointers -::-
!     ----------------------------------------
!
      procedure :: cvs_rho_aibj_projection => cvs_rho_aibj_projection_ccsd
      procedure :: cvs_residual_projection => cvs_residual_projection_ccsd
!
!
!     -::- Other class routine pointers not located in submodules -::-
!     ----------------------------------------------------------------
!
!     Routine to construct right projection vector (eta)
!
      procedure :: construct_eta => construct_eta_ccsd
!
!     Routine to allocate & zero the amplitudes
!
      procedure :: initialize_amplitudes        => initialize_amplitudes_ccsd
      procedure :: initialize_double_amplitudes => initialize_double_amplitudes_ccsd
!
!     Routines to deallocate amplitudes
!
      procedure :: destruct_amplitudes        => destruct_amplitudes_ccsd
      procedure :: destruct_double_amplitudes => destruct_double_amplitudes_ccsd
!
!     Routine to save and read the amplitudes
!
      procedure :: save_amplitudes        => save_amplitudes_ccsd
      procedure :: save_double_amplitudes => save_double_amplitudes_ccsd
!
      procedure :: read_amplitudes        => read_amplitudes_ccsd
      procedure :: read_double_amplitudes => read_double_amplitudes_ccsd
!
!     Set the double amplitudes to the MP2 guess
!
      procedure :: construct_perturbative_doubles => construct_perturbative_doubles_ccsd
!
!     Routine to calculate the energy
!
      procedure :: calc_energy => calc_energy_ccsd
!
   end type ccsd
!
!
!  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  -::- 3. Interfaces to the submodules of the CCSD class -::-
!  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
   interface
!
!
!     -::- Ground state submodule interface -::-
!     ------------------------------------------
!
      module subroutine calc_ampeqs_norm_ccsd(wf, ampeqs_norm)
!!
!!       Calculate amplitude equations norm (CCSD)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp) :: ampeqs_norm
!
      end subroutine calc_ampeqs_norm_ccsd
!
!
      module subroutine new_amplitudes_ccsd(wf, diis_ground_state)
!!
!!       New amplitudes (CCSD)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         class(diis) :: diis_ground_state
!
      end subroutine new_amplitudes_ccsd
!
!
      module subroutine calc_quasi_Newton_doubles_ccsd(wf,dt)
!!
!!       Calculate quasi-Newton estimate (CCSD)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension(wf%n_parameters, 1) :: dt
!
      end subroutine calc_quasi_Newton_doubles_ccsd
!
!
      module subroutine ground_state_preparations_ccsd(wf)
!!
!!       Ground state preparations (CCSD)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
         class(ccsd) :: wf
!
      end subroutine ground_state_preparations_ccsd
!
!
      module subroutine summary_ground_state_info_ccsd(wf, time)
!!
!!       Summary ground state info (CCSD)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2018
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp) :: time
!
      end subroutine summary_ground_state_info_ccsd
!
!
   end interface
!
!
   interface
!
!
!     -::- Omega submodule interface -::-
!     -----------------------------------
!
      module subroutine initialize_omega_ccsd(wf)
!!
!!       Initialize omega (CCSD)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
      end subroutine initialize_omega_ccsd
!
!
      module subroutine destruct_omega_ccsd(wf)
!!
!!       Destruct Omega (CCSD)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
      end subroutine destruct_omega_ccsd
!
!
      module subroutine construct_omega_ccsd(wf)
!!
!!       Construct omega (CCSD)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
      end subroutine construct_omega_ccsd
!
!
      module subroutine omega_ccsd_a1_ccsd(wf)
!!
!!       Omega CCSD A1 term
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
      end subroutine omega_ccsd_a1_ccsd
!
!
      module subroutine omega_ccsd_b1_ccsd(wf)
!!
!!       Omega CCSD B1
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
      end subroutine omega_ccsd_b1_ccsd
!
!
      module subroutine omega_ccsd_c1_ccsd(wf)
!!
!!       Omega CCSD C1
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
      end subroutine omega_ccsd_c1_ccsd
!
!
      module subroutine omega_ccsd_a2_ccsd(wf)
!!
!!       Omega CCSD A2
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
      end subroutine omega_ccsd_a2_ccsd
!
!
      module subroutine omega_ccsd_b2_ccsd(wf)
!!
!!       Omega CCSD B2
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
      end subroutine omega_ccsd_b2_ccsd
!
!
      module subroutine omega_ccsd_c2_ccsd(wf)
!!
!!       Omega CCSD C2
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
      end subroutine omega_ccsd_c2_ccsd
!
!
      module subroutine omega_ccsd_d2_ccsd(wf)
!!
!!       Omega CCSD D2
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
      end subroutine omega_ccsd_d2_ccsd
!
!
      module subroutine omega_ccsd_e2_ccsd(wf)
!!
!!       Omega E2
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
      end subroutine omega_ccsd_e2_ccsd
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
      module subroutine calculate_orbital_differences_ccsd(wf,orbital_diff)
!!
!!       Calculate orbital differences (CCSD)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
         implicit none
!
         class(ccsd) :: wf
         real(dp), dimension(wf%n_parameters, 1) :: orbital_diff
!
      end subroutine calculate_orbital_differences_ccsd
!
!
      module subroutine transform_trial_vectors_ccsd(wf, first_trial, last_trial)
!!
!!       Transformation trial vectors (CCSD)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         integer(i15), intent(in) :: first_trial, last_trial ! Which trial_vectors we are to transform
!
      end subroutine transform_trial_vectors_ccsd
!
!
      module subroutine print_excitation_vector_ccsd(wf, vec, unit_id)
!!
!!       Print excitation vector
!!       Written by Sarai D. Folkestad, Oct 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension(wf%n_parameters, 1) :: vec
!
         integer(i15) :: unit_id
!
      end subroutine print_excitation_vector_ccsd
!
!
      module subroutine cvs_residual_projection_ccsd(wf, residual)
!!
!!       CVS residual projection (CCSD)
!!       Written by Sarai D. Folkestad, Aug 2017
!!
         implicit none
!
         class(ccsd) :: wf
         real(dp), dimension(wf%n_parameters, 1) :: residual
!
      end subroutine cvs_residual_projection_ccsd
!
!
      module subroutine excited_state_preparations_ccsd(wf)
!!
!!       Excited State Preparations (CCSD)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
         class(ccsd) :: wf
!
      end subroutine excited_state_preparations_ccsd
!
!
      module subroutine summary_excited_state_info_ccsd(wf, energies)
!!
!!       Summary of excited state info (CCSD)
!!       Written by Sarai D. Folkestad, Oct 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension(wf%excited_state_specifications%n_singlet_states,1) :: energies
!
      end subroutine summary_excited_state_info_ccsd
!
!
   end interface
!
!
   interface
!
!
!     -::- Ionized state submodule interface -::-
!     -------------------------------------------
!
      module subroutine ionization_residual_projection_ccsd(wf, residual)
!!
!!       Residual projection for CVS
!!       Written by Sarai D. Folkestad, Aug 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension(wf%n_parameters, 1) :: residual
!
      end subroutine ionization_residual_projection_ccsd
!
!
      module subroutine ionization_rho_aibj_projection_ccsd(wf, rho_aibj)
!!
!!       Residual projection for CVS
!!       Written by Sarai D. Folkestad, Aug 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension(:,:) :: rho_aibj
!
      end subroutine ionization_rho_aibj_projection_ccsd
!
!
      module subroutine ionization_jacobian_ccsd_transformation_ccsd(wf, c_a_i, c_aibj)
!!
!!      Ionization Jacobian transformation (CCSD)
!!      Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
        implicit none
!
        class(ccsd) :: wf
!
        real(dp), dimension(wf%n_v, wf%n_o) :: c_a_i
        real(dp), dimension(wf%n_t2am, 1)   :: c_aibj
!
      end subroutine ionization_jacobian_ccsd_transformation_ccsd
!
!
      module subroutine core_ionization_jacobian_ccsd_transformation_ccsd(wf, c_a_i, c_aibj)
!!
!!       Core ionization Jacobian transformation (CCSD)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension(wf%n_v, wf%n_o) :: c_a_i
         real(dp), dimension(wf%n_t2am, 1)   :: c_aibj
!
      end subroutine core_ionization_jacobian_ccsd_transformation_ccsd
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
      module subroutine jacobian_ccsd_transformation_ccsd(wf, c_a_i, c_aibj)
!!
!!       Jacobian CCSD transformation
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension(wf%n_v, wf%n_o) :: c_a_i  ! c_ai
         real(dp), dimension(wf%n_t2am, 1)   :: c_aibj ! c_aibj
!
      end subroutine jacobian_ccsd_transformation_ccsd
!
!
      module subroutine cvs_jacobian_ccsd_transformation_ccsd(wf, c_a_i, c_aibj)
!!
!!       CVS Jacobian transformation (CCSD)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension(wf%n_v, wf%n_o) :: c_a_i  ! c_ai
         real(dp), dimension(wf%n_t2am, 1)   :: c_aibj ! c_aibj
!
      end subroutine cvs_jacobian_ccsd_transformation_ccsd
!
!
      module subroutine jacobian_ccsd_a1_ccsd(wf, rho_a_i, c_a_i)
!!
!!       Jacobian CCSD A1
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_a_i   ! c_ai
         real(dp), dimension(wf%n_v, wf%n_o) :: rho_a_i ! rho_ai
!
      end subroutine jacobian_ccsd_a1_ccsd
!
!
      module subroutine jacobian_ccsd_b1_ccsd(wf, rho_a_i, c_ai_bj)
!!
!!       Jacobian CCSD B1
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)), intent(in) :: c_ai_bj   ! c_aibj
         real(dp), dimension(wf%n_v, wf%n_o) :: rho_a_i ! rho_ai
!
      end subroutine jacobian_ccsd_b1_ccsd
!
!
      module subroutine jacobian_ccsd_c1_ccsd(wf, rho_a_i, c_ai_bj)
!!
!!       Jacobian CCSD C1
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)), intent(in) :: c_ai_bj
         real(dp), dimension(wf%n_v, wf%n_o) :: rho_a_i ! rho_ai
!
      end subroutine jacobian_ccsd_c1_ccsd
!
!
      module subroutine jacobian_ccsd_d1_ccsd(wf, rho_a_i, c_bi_cj)
!!
!!       Jacobian CCSD D1
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)), intent(in) :: c_bi_cj
         real(dp), dimension(wf%n_v, wf%n_o) :: rho_a_i ! rho_ai
!
      end subroutine jacobian_ccsd_d1_ccsd
!
!
      module subroutine jacobian_ccsd_a2_ccsd(wf, rho_ai_bj, c_a_i)
!!
!!       Jacobian CCSD A2
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
         real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_a_i
!
      end subroutine jacobian_ccsd_a2_ccsd
!
!
      module subroutine jacobian_ccsd_b2_ccsd(wf, rho_ai_bj, c_a_i)
!!
!!       Jacobian CCSD B2
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_a_i
!
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
!
      end subroutine jacobian_ccsd_b2_ccsd
!
!
      module subroutine jacobian_ccsd_c2_ccsd(wf, rho_ai_bj, c_a_i)
!!
!!       Jacobian CCSD C2
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_a_i
!
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
!
      end subroutine jacobian_ccsd_c2_ccsd
!
!
      module subroutine jacobian_ccsd_d2_ccsd(wf, rho_ai_bj, c_a_i)
!!
!!       Jacobian CCSD D2
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_a_i
!
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
!
      end subroutine jacobian_ccsd_d2_ccsd
!
!
      module subroutine jacobian_ccsd_e2_ccsd(wf, rho_ai_bj, c_ai_ck)
!!
!!       Jacobian CCSD E2
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
!
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)), intent(in) :: c_ai_ck
!
      end subroutine jacobian_ccsd_e2_ccsd
!
!
      module subroutine jacobian_ccsd_f2_ccsd(wf, rho_ai_bj, c_ai_bj)
!!
!!       Jacobian CCSD F2
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
!
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)), intent(in) :: c_ai_bj
!
      end subroutine jacobian_ccsd_f2_ccsd
!
!
      module subroutine jacobian_ccsd_g2_ccsd(wf, rho_ai_bj, c_ai_bj)
!!
!!       Jacobian CCSD G2
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)), intent(in) :: c_ai_bj
!
      end subroutine jacobian_ccsd_g2_ccsd
!
!
      module subroutine jacobian_ccsd_h2_ccsd(wf, rho_ai_bj, c_ai_bj)
!!
!!       Jacobian CCSD H2
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)), intent(in) :: c_ai_bj
!
      end subroutine jacobian_ccsd_h2_ccsd
!
!
      module subroutine jacobian_ccsd_i2_ccsd(wf, rho_ai_bj, c_ai_bj)
!!
!!       Jacobian CCSD I2
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)), intent(in) :: c_ai_bj
!
      end subroutine jacobian_ccsd_i2_ccsd
!
!
      module subroutine jacobian_ccsd_j2_ccsd(wf, rho_ab_ij, c_ab_ij)
!!
!!       Jacobian CCSD J2
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension((wf%n_v)**2, (wf%n_o)**2) :: rho_ab_ij
         real(dp), dimension((wf%n_v)**2, (wf%n_o)**2), intent(in) :: c_ab_ij
!
      end subroutine jacobian_ccsd_j2_ccsd
!
!
      module subroutine jacobian_ccsd_k2_ccsd(wf, rho_ab_ij, c_ab_ij)
!!
!!       Jacobian CCSD K2
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension((wf%n_v)**2, (wf%n_o)**2) :: rho_ab_ij
         real(dp), dimension((wf%n_v)**2, (wf%n_o)**2), intent(in) :: c_ab_ij
!
      end subroutine jacobian_ccsd_k2_ccsd
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
      module subroutine jacobian_transpose_ccsd_transformation_ccsd(wf, b_a_i, b_aibj)
!!
!!       Jacobian transpose transformation (CCSD)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension(wf%n_v, wf%n_o) :: b_a_i
         real(dp), dimension(wf%n_t2am, 1)   :: b_aibj
!
      end subroutine jacobian_transpose_ccsd_transformation_ccsd
!
!
      module subroutine jacobian_transpose_ccsd_a1_ccsd(wf, sigma_a_i, b_a_i)
!!
!!       Jacobian transpose CCSD A1
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension(wf%n_v, wf%n_o) :: b_a_i
         real(dp), dimension(wf%n_v, wf%n_o) :: sigma_a_i
!
      end subroutine jacobian_transpose_ccsd_a1_ccsd
!
!
      module subroutine jacobian_transpose_ccsd_b1_ccsd(wf, sigma_a_i, b_a_i)
!!
!!       Jacobian transpose CCSD B1
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension(wf%n_v, wf%n_o) :: b_a_i
         real(dp), dimension(wf%n_v, wf%n_o) :: sigma_a_i
!
      end subroutine jacobian_transpose_ccsd_b1_ccsd
!
!
      module subroutine jacobian_transpose_ccsd_c1_ccsd(wf, sigma_a_i, b_ai_bj)
!!
!!       Jacobian transpose CCSD C1
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension(wf%n_v, wf%n_o)                       :: sigma_a_i
         real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: b_ai_bj
!
      end subroutine jacobian_transpose_ccsd_c1_ccsd
!
!
      module subroutine jacobian_transpose_ccsd_d1_ccsd(wf, sigma_a_i, b_ai_bj)
!!
!!       Jacobian transpose CCSD D1
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension(wf%n_v, wf%n_o)                       :: sigma_a_i
         real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: b_ai_bj
!
      end subroutine jacobian_transpose_ccsd_d1_ccsd
!
!
      module subroutine jacobian_transpose_ccsd_e1_ccsd(wf, sigma_a_i, b_ai_bj)
!!
!!       Jacobian transpose CCSD E1
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension(wf%n_v, wf%n_o)                       :: sigma_a_i
         real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: b_ai_bj
!
      end subroutine jacobian_transpose_ccsd_e1_ccsd
!
!
      module subroutine jacobian_transpose_ccsd_f1_ccsd(wf, sigma_a_i, b_ai_bj)
!!
!!       Jacobian transpose CCSD F1
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension(wf%n_v, wf%n_o)                       :: sigma_a_i
         real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: b_ai_bj
!
      end subroutine jacobian_transpose_ccsd_f1_ccsd
!
!
      module subroutine jacobian_transpose_ccsd_g1_ccsd(wf, sigma_a_i, b_ai_bj)
!!
!!       Jacobian transpose CCSD G1
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension(wf%n_v, wf%n_o)                       :: sigma_a_i
         real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: b_ai_bj
!
      end subroutine jacobian_transpose_ccsd_g1_ccsd
!
!
      module subroutine jacobian_transpose_ccsd_a2_ccsd(wf, sigma_ai_bj, b_a_i)
!!
!!       Jacobian transpose CCSD A2
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension(wf%n_v, wf%n_o)                       :: b_a_i
         real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: sigma_ai_bj
!
      end subroutine jacobian_transpose_ccsd_a2_ccsd
!
!
      module subroutine jacobian_transpose_ccsd_b2_ccsd(wf, sigma_ai_bj, b_ai_bj)
!!
!!       Jacobian transpose CCSD B2
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: b_ai_bj
         real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: sigma_ai_bj
!
      end subroutine jacobian_transpose_ccsd_b2_ccsd
!
!
      module subroutine jacobian_transpose_ccsd_c2_ccsd(wf, sigma_ai_bj, b_ai_bj)
!!
!!       Jacobian transpose CCSD C2
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: b_ai_bj
         real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: sigma_ai_bj
!
      end subroutine jacobian_transpose_ccsd_c2_ccsd
!
!
      module subroutine jacobian_transpose_ccsd_d2_ccsd(wf, sigma_ai_bj, b_ai_bj)
!!
!!       Jacobian transpose CCSD D2
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: b_ai_bj
         real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: sigma_ai_bj
!
      end subroutine jacobian_transpose_ccsd_d2_ccsd
!
!
      module subroutine jacobian_transpose_ccsd_e2_ccsd(wf, sigma_ai_bj, b_ai_bj)
!!
!!       Jacobian transpose CCSD E2
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: b_ai_bj
         real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: sigma_ai_bj
!
      end subroutine jacobian_transpose_ccsd_e2_ccsd
!
!
      module subroutine jacobian_transpose_ccsd_f2_ccsd(wf, sigma_ai_bj, b_ai_bj)
!!
!!       Jacobian transpose CCSD F2
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: b_ai_bj
         real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: sigma_ai_bj
!
      end subroutine jacobian_transpose_ccsd_f2_ccsd
!
!
      module subroutine jacobian_transpose_ccsd_g2_ccsd(wf, sigma_ai_bj, b_ai_bj)
!!
!!       Jacobian transpose CCSD G2
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: b_ai_bj
         real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: sigma_ai_bj
!
      end subroutine jacobian_transpose_ccsd_g2_ccsd
!
!
      module subroutine jacobian_transpose_ccsd_h2_ccsd(wf, sigma_ab_ij, b_ab_ij)
!!
!!       Jacobian transpose CCSD H2
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension((wf%n_v)**2, (wf%n_o)**2) :: b_ab_ij
         real(dp), dimension((wf%n_v)**2, (wf%n_o)**2) :: sigma_ab_ij
!
      end subroutine jacobian_transpose_ccsd_h2_ccsd
!
!
      module subroutine jacobian_transpose_ccsd_i2_ccsd(wf, sigma_ab_ij, b_ab_ij)
!!
!!       Jacobian transpose CCSD I2
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension((wf%n_v)**2, (wf%n_o)**2) :: b_ab_ij
         real(dp), dimension((wf%n_v)**2, (wf%n_o)**2) :: sigma_ab_ij
!
      end subroutine jacobian_transpose_ccsd_i2_ccsd
!
!
      module subroutine cvs_rho_aibj_projection_ccsd(wf, vec_aibj)
!!
!!       CVS rho_aibj projection (CCSD)
!!       Written by Sarai D. Folkestad, Aug 2017
!!
         implicit none
!
         class(ccsd) :: wf
         real(dp), dimension(:, :) :: vec_aibj
!
      end subroutine cvs_rho_aibj_projection_ccsd
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
   subroutine init_ccsd(wf)
!!
!!    Initialize CCSD object
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!    Wavefunction initialization. The tasks performed are, in broad terms:
!!
!!    a) Read the input file to determine general specifications
!!       (model type, memory, disk space, etc.), the requested
!!       calculations (ground state, excited state, etc.) with associated
!!       settings (thresholds, number of iterations, etc.). For more information,
!!       see the input reader submodule for the reader-type routines.
!!
!!    b) Set wavefunction settings. This includes the name of the object,
!!       used mainly for printing, and also which methods are implemented,
!!       which is used to stop calculations that are not possible to execute
!!       given what is currently implemented. Moreover, a couple of central
!!       attributes are set (number of occupied and virtual orbitals, single
!!       and double amplitudes, and the number of variables, which is used
!!       by solvers in particular)
!!
!!    c) Hartree-Fock Dalton interface. This includes reading Hartree-Fock
!!       results (orbital coefficients and orbital attributes, see b)) and
!!       Cholesky decomposition files (for atomic Cholesky vectors).
!!
!!    One or more of a)-c) are different for different wavefunctions,
!!    so most new wavefunction has to have an overwritten initialization
!!    routine. This is in contrast to the driver subroutines, which are
!!    mostly kept in the CCS form for the entire hierarchy.
!!
      implicit none
!
      class(ccsd) :: wf
!
      type(file) :: input_file
!
!     Set model name
!
      wf%name = 'CCSD'
!
!     Open input file eT.inp
!
      input_file%name = 'eT.inp'
      call wf%disk%open_file(input_file, 'formatted', 'read', 'sequential')
!
!     Read general specifications (memory and diskspace for calculation)
!
      call wf%general_specs_reader(input_file%unit)
!
!     Set implemented methods
!
      wf%implemented%ground_state       = .true.
      wf%implemented%excited_state      = .true.
      wf%implemented%core_excited_state = .true.
      wf%implemented%ionized_state      = .true.
      wf%implemented%multipliers        = .true.
!
!     Read calculation tasks from input file eT.inp
!
      call wf%calculation_reader(input_file%unit)
!
!     Close input file
!
      call wf%disk%close_file(input_file)
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
!     Set single and double amplitude attributes & the total number of variables
!
      wf%n_t1am = (wf%n_o)*(wf%n_v)
      wf%n_t2am = (wf%n_t1am)*(wf%n_t1am + 1)/2
!
      wf%n_parameters = wf%n_t1am + wf%n_t2am
!
!     Initialize the Fock matrix (allocate and construct given the initial amplitudes)
!
      call wf%initialize_single_amplitudes
!
      call wf%initialize_fock_matrix
!
      call wf%destruct_single_amplitudes
!
   end subroutine init_ccsd
!
!
   subroutine initialize_amplitudes_ccsd(wf)
!!
!!    Initialize amplitudes (CCSD)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!    Allocates the amplitudes and sets them to zero.
!!
      implicit none
!
      class(ccsd) :: wf
!
      call wf%initialize_single_amplitudes
      call wf%initialize_double_amplitudes
!
   end subroutine initialize_amplitudes_ccsd
!
!
   subroutine initialize_double_amplitudes_ccsd(wf)
!!
!!    Initialize double amplitudes (CCSD)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!    Allocates the doubles amplitudes if they are not allocated,
!!    then sets them to zero.
!!
      implicit none
!
      class(ccsd) :: wf
!
!     Allocate the doubles amplitudes and set to zero
!
      if (.not. allocated(wf%t2am)) then
!
         call wf%mem%alloc(wf%t2am, wf%n_t2am, 1)
!
      else
!
         write(unit_output,'(/t3,a)') 'Warning: attempted to allocate and zero already allocated t2am.'
         write(unit_output,'(t3,a)')  'This is a bug. Be aware that the vector is zeroed anyway.'
!
      endif
!
      wf%t2am = zero
!
   end subroutine initialize_double_amplitudes_ccsd
!
!
   subroutine construct_perturbative_doubles_ccsd(wf)
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
      class(ccsd) :: wf
!
      real(dp), dimension(:,:), allocatable :: g_ia_jb
!
      integer(i15) :: i = 0, j = 0, a = 0, b = 0
      integer(i15) :: ai = 0, bj = 0, ia = 0, jb = 0, aibj = 0
!
      character(len=40) :: integral_type
!
!     Form g_ia_jb
!
      call wf%mem%alloc(g_ia_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      integral_type = 'electronic_repulsion'
      call wf%get_ov_ov(integral_type, g_ia_jb)
!
!     Set the doubles amplitudes
!
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            ai = index_two(a, i, wf%n_v)
            ia = index_two(i, a, wf%n_o)
!
            do j = 1, wf%n_o
               do b = 1, wf%n_v
!
                  jb = index_two(j, b, wf%n_o)
                  bj = index_two(b, j, wf%n_v)
!
!                 Set the doubles amplitudes
!
                  if (ai .le. bj) then ! To avoid setting the same element twice
!
                     aibj = index_packed(ai,bj)
!
                     wf%t2am(aibj, 1) = - g_ia_jb(ia,jb)/(wf%fock_diagonal(wf%n_o + a, 1) + &
                                                          wf%fock_diagonal(wf%n_o + b, 1) - &
                                                          wf%fock_diagonal(i, 1) -          &
                                                          wf%fock_diagonal(j, 1))
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
      call wf%mem%dealloc(g_ia_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
   end subroutine construct_perturbative_doubles_ccsd
!
!
   subroutine calc_energy_ccsd(wf)
!!
!!     Calculate energy (CCSD)
!!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!     Calculates the CCSD energy. This is only equal to the actual
!!     energy when the ground state equations are solved, of course.
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(:,:), allocatable :: L_ia_J  ! L_ia^J
      real(dp), dimension(:,:), allocatable :: g_ia_jb ! g_iajb
!
      integer(i15) :: a = 0, i = 0, b = 0, j = 0, ai = 0
      integer(i15) :: bj = 0, aibj = 0, ia = 0, jb = 0, ib = 0, ja = 0
!
      character(len=40) :: integral_type
!
!     Get g_ia_jb = g_iajb
!
      call wf%mem%alloc(g_ia_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      integral_type = 'electronic_repulsion'
      call wf%get_ov_ov(integral_type, g_ia_jb)
!
!     Set the initial value of the energy
!
      wf%energy = wf%scf_energy
!
!     Add the correlation energy E = E + sum_aibj (t_ij^ab + t_i^a t_j^b) L_iajb
!
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            ai = index_two(a, i, wf%n_v)
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
!                 Add the correlation energy
!
                  wf%energy = wf%energy +                                           &
                                 (wf%t2am(aibj,1) + (wf%t1am(a,i))*(wf%t1am(b,j)))* &
                                 (two*g_ia_jb(ia,jb) - g_ia_jb(ib,ja))
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate g_ia_jb
!
      call wf%mem%dealloc(g_ia_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
   end subroutine calc_energy_ccsd
!
!
   subroutine destruct_amplitudes_ccsd(wf)
!!
!!    Destruct amplitudes (CCSD)
!!    Written by Sarai D. Folkestad and Eirik F. Kjøsntad, May 2017
!!
!!    Deallocates the amplitudes. The reason for including a special
!!    destruct routine is that sometimes one may want to make sure the
!!    amplitudes are deallocated, without wanting to actually deallocate
!!    an existing vector (introducing errors in available memory estimation)
!!
      implicit none
!
      class(ccsd) :: wf
!
      call wf%destruct_single_amplitudes
      call wf%destruct_double_amplitudes
!
   end subroutine destruct_amplitudes_ccsd
!
!
   subroutine destruct_double_amplitudes_ccsd(wf)
!!
!!    Destruct double amplitudes (CCSD)
!!    Written by Sarai D. Folkestad and Eirik F. Kjøsntad, May 2017
!!
!!    Deallocates the doubles amplitudes if they are allocated.
!!
      implicit none
!
      class(ccsd) :: wf
!
      if (allocated(wf%t2am)) call wf%mem%dealloc(wf%t2am, wf%n_t2am, 1)
!
   end subroutine destruct_double_amplitudes_ccsd
!
!
   subroutine save_amplitudes_ccsd(wf)
!!
!!    Save amplitudes (CCSD)
!!    Written by Sarai D. Folkestad and Eirik F. Kjøsntad, May 2017
!!
!!    Store the amplitudes to disk (in the files t1am and t2am)
!!
      implicit none
!
      class(ccsd) :: wf
!
      call wf%save_single_amplitudes
      call wf%save_double_amplitudes
!
   end subroutine save_amplitudes_ccsd
!
!
   subroutine save_double_amplitudes_ccsd(wf)
!!
!!    Save double amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjøsntad, May 2017
!!
!!    Stores the double amplitudes to disk
!!
      implicit none
!
      class(ccsd) :: wf
!
      type(file) :: t2am_file
!
!     Open amplitude file
!
      t2am_file%name = 't2am'
      call wf%disk%open_file(t2am_file, 'unformatted', 'write', 'sequential')
!
!     Write amplitudes to file
!
      write(t2am_file%unit) wf%t2am
!
!     Close amplitude files
!
      call wf%disk%close_file(t2am_file)
!
   end subroutine save_double_amplitudes_ccsd
!
!
   subroutine read_amplitudes_ccsd(wf)
!!
!!    Read amplitudes (CCSD)
!!    Written by Sarai D. Folkestad and Eirik F. Kjøsntad, May 2017
!!
!!    Reads the single and double amplitudes (t1am and t2am) from disk.
!!    If the arrays are not allocated, the routine also allocates them.
!!
      implicit none
!
      class(ccsd) :: wf
!
!     Read the amplitudes from disk
!
      call wf%read_single_amplitudes
      call wf%read_double_amplitudes
!
   end subroutine read_amplitudes_ccsd
!
!
   subroutine read_double_amplitudes_ccsd(wf)
!!
!!    Read double amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017 / June 2018
!!
!!    Reads the double amplitudes from disk
!!
      implicit none
!
      class(ccsd) :: wf
!
      type(file) :: t2am_file
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
         t2am_file%name = 't2am'
         call wf%disk%open_file(t2am_file, 'unformatted', 'read', 'sequential')
!
!        Allocate doubles amplitudes if they aren't allocated
!
         if (.not. allocated(wf%t2am)) then
!
            call wf%mem%alloc(wf%t2am, wf%n_t2am, 1)
            wf%t2am = zero
!
         endif
!
!        Read from file & close
!
         read(t2am_file%unit) wf%t2am
!
         call wf%disk%close_file(t2am_file)
!
      else
!
         write(unit_output,'(t3,a)') 'Error: attempted to read double amplitudes from file,'
         write(unit_output,'(t3,a)') 'but the file does not appear to exist.'
         stop
!
      endif
!
   end subroutine read_double_amplitudes_ccsd
!
!
   subroutine jacobi_test_ccsd(wf)
!!
!!    Jacobian test (CCSD)
!!    Written by Eirik F. Kjønstad, June 2017
!!
!!    Calculates the Jacobian matrix by numerically differentiating
!!    the projection vector:
!!
!!       A_mu,nu = d Omega_mu / d t_nu.
!!
!!    Used to debug Jacobian transformation.
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(:,:), allocatable :: c_a_i
      real(dp), dimension(:,:), allocatable :: c_aibj
!
      integer(i15) :: a = 0, i = 0, b = 0, j = 0, ai = 0, bj = 0, aibj = 0
      integer(i15) :: c = 0, k = 0, ck = 0, ckdl = 0, l = 0, d = 0, dl = 0
!
      real(dp), dimension(:,:), allocatable :: r1am
      real(dp), dimension(:,:), allocatable :: r2am
!
      real(dp) :: displacement
!
!     Calculate the transformation of the t1 amplitudes
!
      call wf%read_double_amplitudes
!
      call wf%mem%alloc(r1am, wf%n_v, wf%n_o)
      call wf%mem%alloc(r2am, wf%n_t2am, 1)
!
      r1am = wf%t1am
      r2am = wf%t2am
!
!     Make sure fock matrix is up to date
!
      call wf%construct_fock
!
   !   r1am = zero
   !   r2am = zero
!
      call wf%destruct_amplitudes
!
      write(unit_output,*) 'T1AM'
!
      do i = 1, wf%n_o
         do a = 1, wf%n_v
            ai = index_two(a, i, wf%n_v)
            write(unit_output,*) ai,r1am(a,i)
         enddo
      enddo
!
      write(unit_output,*) 'T2AM'
!
      do j = 1, 15
         write(unit_output,*) j, r2am(j,1)
      enddo
!
   !   call wf%jacobian_ccsd_transformation(r1am,r2am)
      call wf%jacobian_transpose_ccsd_transformation(r1am,r2am)
!
      write(unit_output,*) 'TRF(SINGLES)'
!
      do i = 1, wf%n_o
         do a = 1, wf%n_v
            ai = index_two(a, i, wf%n_v)
            write(unit_output,*) ai,r1am(a,i)
         enddo
      enddo
!
      write(unit_output,*) 'TRF(DOUBLES)'
!
      do j = 1, 500
         write(unit_output,*) j, r2am(j,1)
      enddo
!
!     We wish to calculate A_mu,nu = d(omega)_mu / dt_nu,
!     in two different ways:
!
!        - By transforming c_tau = delta_tau,nu with A. Then (A c)_mu = sum_tau A_mu,tau c_tau = A_mu,nu
!        - By calculating A_mu,nu = (omega(t+Dt_nu)_mu - omega(t)_mu)/Dt_nu.
!
!     :: First approach: transform by A ::
!
      call wf%mem%alloc(c_a_i, wf%n_v, wf%n_o)
      call wf%mem%alloc(c_aibj, wf%n_t2am, 1)
!
      c_a_i  = zero
      c_aibj = zero
!
      write(unit_output,*) 'A_mu,nu, singles, by transformation'
!
      do j = 1, wf%n_o
         do b = 1, wf%n_v
!
!           Let c have a one for the excitation bj, zero otherwise
!
            c_a_i(b, j) = one
!
!           Transform c by A. The result should be A_mu,bj
!
            call wf%jacobian_ccsd_transformation(c_a_i,c_aibj)
!
!           Print these elements of A.
!
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  ai = index_two(a, i, wf%n_v)
                  bj = index_two(b, j, wf%n_v)
!
                  write(unit_output,*) 'ai, bj, A_ai,bj', ai, bj, c_a_i(a, i)
!
               enddo
            enddo
!
            c_a_i = zero
            c_aibj = zero
!
         enddo
      enddo
!
      write(unit_output,*) 'A_ckdl,bj, singles, by transformation ZZ'
!
      do j = 1, wf%n_o
         do b = 1, wf%n_v
!
!           Let c have a one for the excitation bj, zero otherwise
!
            c_a_i(b, j) = one
!
!           Transform c by A. The result should be A_mu,bj
!
            call wf%jacobian_ccsd_transformation(c_a_i,c_aibj)
!
!           Print these elements of A.
!
            do k = 1, wf%n_o
               do c = 1, wf%n_v
                  do l = 1, wf%n_o
                     do d = 1, wf%n_v
!
                        ck = index_two(c, k, wf%n_v)
                        dl = index_two(d, l, wf%n_v)
!
                        ckdl = index_packed(ck, dl)
!
                        bj = index_two(b, j, wf%n_v)
!
                        write(unit_output,*) 'ckdl, bj, A_ckdl,bj', ckdl, bj, c_aibj(ckdl,1)
!
                     enddo
                  enddo
!
               enddo
            enddo
!
            c_a_i = zero
            c_aibj = zero
!
         enddo
      enddo
!
      c_a_i  = zero
      c_aibj = zero
!
      write(unit_output,*) 'A_ck,aibj, doubles, by transformation XX'
!
      do j = 1, wf%n_o
         do b = 1, wf%n_v
!
            bj = index_two(b, j, wf%n_v)
!
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  ai = index_two(a, i, wf%n_v)
!
                  aibj = index_packed(ai, bj)
!
                  c_aibj(aibj,1) = one
!
!                 Transform c by A. The result should be A_mu,aibj
!
                  call wf%jacobian_ccsd_transformation(c_a_i,c_aibj)
!
!                 Print singles-doubles block
!
                  do k = 1, wf%n_o
                     do c = 1, wf%n_v
!
                           ck = index_two(c, k, wf%n_v)
!
                           write(unit_output,*) 'ck, aibj, A_ck, aibj', ck, aibj, c_a_i(c, k)
!
                     enddo
                  enddo
!
                  c_a_i = zero
                  c_aibj = zero
!
               enddo
            enddo
         enddo
      enddo
!
      write(unit_output,*) 'A_ckdl,aibj, doubles, by transformation YY'
!
      do j = 1, wf%n_o
         do b = 1, wf%n_v
!
            bj = index_two(b, j, wf%n_v)
!
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  ai = index_two(a, i, wf%n_v)
!
                  aibj = index_packed(ai, bj)
!
                  c_aibj(aibj,1) = one
!
!                 Transform c by A. The result should be A_mu,aibj
!
                  call wf%jacobian_ccsd_transformation(c_a_i,c_aibj)
!
!                 Print singles-doubles block
!
                  do k = 1, wf%n_o
                     do c = 1, wf%n_v
                        do l = 1, wf%n_o
                           do d = 1, wf%n_v
!
                              ck = index_two(c, k, wf%n_v)
                              dl = index_two(d, l, wf%n_v)
!
                              ckdl = index_packed(ck, dl)
!
                              write(unit_output,*) 'ckdl, aibj, A_ckdl, aibj', ckdl, aibj, c_aibj(ckdl,1)
!
                           enddo
                        enddo
                     enddo
                  enddo
!
                  c_a_i = zero
                  c_aibj = zero
!
               enddo
            enddo
         enddo
      enddo
!
!     :: Second approach: differentaition of omega ::
!
      write(unit_output,*) 'A_mu,nu, singles, by derivation'
!
      call wf%initialize_amplitudes
      call wf%read_double_amplitudes
!
      displacement = 1.0D-10
!
      call wf%initialize_omega
!
      do j = 1, wf%n_o
         do b = 1, wf%n_v
!
!           Construct omega and save in c
!
            wf%omega1 = zero
            wf%omega2 = zero
!
            call wf%construct_fock
            call wf%construct_omega
!
            c_a_i  = wf%omega1
            c_aibj = wf%omega2
!
!           Shift the j,b amplitude
!
            wf%t1am(b,j) = wf%t1am(b,j) + displacement
!
!           Construct omega
!
            wf%omega1 = zero
            wf%omega2 = zero
!
            call wf%construct_fock
            call wf%construct_omega
!
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
!                 Calculate A_ai,bj = (omega(t)_ai - omega(t+Dt_bj)_ai)/Dt_bj.
!
                  ai = index_two(a, i, wf%n_v)
                  bj = index_two(b, j, wf%n_v)
!
                  write(unit_output,*) 'ai, bj, A_ai,bj', ai, bj, (wf%omega1(a,i)-c_a_i(a,i))/displacement
!
               enddo
            enddo
!
            wf%t1am(b,j) = wf%t1am(b,j) - displacement
!
         enddo
      enddo
!
      write(unit_output,*) 'A_ck,aibj, doubles, by differentiation XX'
!
      do j = 1, wf%n_o
         do b = 1, wf%n_v
!
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  ai = index_two(a, i, wf%n_v)
                  bj = index_two(b, j, wf%n_v)
!
                  aibj = index_packed(ai, bj)
!
!                 Construct omega and save in c
!
                  wf%omega1 = zero
                  wf%omega2 = zero
!
                  call wf%construct_fock
                  call wf%construct_omega
!
                  c_a_i  = wf%omega1
                  c_aibj = wf%omega2
!
!                 Shift the aibj amplitude
!
                  wf%t2am(aibj,1) = wf%t2am(aibj,1) + displacement
!
!                 Construct omega
!
                  wf%omega1 = zero
                  wf%omega2 = zero
!
                  call wf%construct_fock
                  call wf%construct_omega
!
                  do k = 1, wf%n_o
                     do c = 1, wf%n_v
!
                        ck = index_two(c, k, wf%n_v)
!
                        write(unit_output,*) 'ck, aibj, A_ck, aibj', ck, aibj, (wf%omega1(c,k)-c_a_i(c,k))/displacement
!
                     enddo
                  enddo
!
                  wf%t2am(aibj,1) = wf%t2am(aibj,1) - displacement
!
               enddo
            enddo
!
         enddo
      enddo
!
!
      write(unit_output,*) 'A_ckdl,aibj, doubles, by differentiation YY'
!
      do j = 1, wf%n_o
         do b = 1, wf%n_v
!
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  ai = index_two(a, i, wf%n_v)
                  bj = index_two(b, j, wf%n_v)
!
                  aibj = index_packed(ai, bj)
!
!                 Construct omega and save in c
!
                  wf%omega1 = zero
                  wf%omega2 = zero
!
                  call wf%construct_fock
                  call wf%construct_omega
!
                  c_a_i  = wf%omega1
                  c_aibj = wf%omega2
!
!                 Shift the aibj amplitude
!
                  wf%t2am(aibj,1) = wf%t2am(aibj,1) + displacement
!
!                 Construct omega
!
                  wf%omega1 = zero
                  wf%omega2 = zero
!
                  call wf%construct_fock
                  call wf%construct_omega
!
                  do k = 1, wf%n_o
                     do c = 1, wf%n_v
                        do l = 1, wf%n_o
                           do d = 1, wf%n_v
!
                              ck = index_two(c, k, wf%n_v)
                              dl = index_two(d, l, wf%n_v)
!
                              ckdl = index_packed(ck, dl)
!
                              if (c .eq. d .and. k .eq. l) then
                                 write(unit_output,*) 'ckdl, aibj, A_ckdl, aibj', &
                                          ckdl, aibj, half*(wf%omega2(ckdl,1)-c_aibj(ckdl,1))/displacement
                              else
!
                                 write(unit_output,*) 'ckdl, aibj, A_ckdl, aibj', &
                                          ckdl, aibj, (wf%omega2(ckdl,1)-c_aibj(ckdl,1))/displacement
                              endif
!
                           enddo
                        enddo
                     enddo
                  enddo
!
                  wf%t2am(aibj,1) = wf%t2am(aibj,1) - displacement
!
               enddo
            enddo
!
         enddo
      enddo
!
     write(unit_output,*) 'A_ckdl,bj, singles, by derivation ZZ'
!
      do j = 1, wf%n_o
         do b = 1, wf%n_v
!
!           Construct omega and save in c
!
            wf%omega1 = zero
            wf%omega2 = zero
!
            call wf%construct_fock
            call wf%construct_omega
!
            c_a_i  = wf%omega1
            c_aibj = wf%omega2
!
!           Shift the j,b amplitude
!
            wf%t1am(b,j) = wf%t1am(b,j) + displacement
!
!           Construct omega
!
            wf%omega1 = zero
            wf%omega2 = zero
!
            call wf%construct_fock
            call wf%construct_omega
!
            do k = 1, wf%n_o
               do c = 1, wf%n_v
                  do l = 1, wf%n_o
                     do d = 1, wf%n_v
!
!                       Calculate A_ckdl,bj = (omega(t)_ckdl - omega(t+Dt_bj)_ckdl)/Dt_bj.
!
                        ck = index_two(c, k, wf%n_v)
                        dl = index_two(d, l, wf%n_v)
                        bj = index_two(b, j, wf%n_v)
!
                        ckdl = index_packed(ck, dl)
!
                        if (c .eq. d .and. k .eq. l) then
!
                           write(unit_output,*) 'ckdl, bj, A_ckdl,bj', ckdl, bj, &
                                    half*(wf%omega2(ckdl,1)-c_aibj(ckdl,1))/displacement
!
                        else
!
                           write(unit_output,*) 'ckdl, bj, A_ckdl,bj', ckdl, bj, &
                                    (wf%omega2(ckdl,1)-c_aibj(ckdl,1))/displacement
                        endif
!
                     enddo
                  enddo
               enddo
            enddo
!
            wf%t1am(b,j) = wf%t1am(b,j) - displacement
!
         enddo
      enddo
!
   end subroutine jacobi_test_ccsd
!
!
   subroutine construct_eta_ccsd(wf, eta)
!!
!!    Construct eta (CCSD)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
!!    Note: the routine assumes that eta is initialized and that the Fock matrix
!!    has been constructed.
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_parameters, 1) :: eta ! eta = ( eta_ai eta_aibj )
!
      real(dp), dimension(:,:), allocatable :: g_ia_jb
      real(dp), dimension(:,:), allocatable :: eta_ai_bj
!
      character(len=40) :: integral_type
!
      integer(i15) :: i = 0, a = 0, j = 0, b = 0, aibj = 0
      integer(i15) :: ib = 0, ja = 0, jb = 0, ia = 0, bj = 0, ai = 0
!
      call wf%construct_fock
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
!     Form g_ia_jb = g_iajb
!
      call wf%mem%alloc(g_ia_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      integral_type = 'electronic_repulsion'
      call wf%get_ov_ov(integral_type, g_ia_jb)
!
!     Form eta_ai_bj = 2* L_iajb = 2 * ( 2 * g_iajb - g_ibja)
!                                = 4 * g_ia_jb(ia,jb) - 2 * g_ia_jb(ib,ja)
!
      call wf%mem%alloc(eta_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      eta_ai_bj = zero
!
      do j = 1, wf%n_o
         do b = 1, wf%n_v
!
            bj = index_two(b, j, wf%n_v)
            jb = index_two(j, b, wf%n_o)
!
            do i = 1, wf%n_o
!
               ib = index_two(i, b, wf%n_o)
!
               do a = 1, wf%n_v
!
                  ai = index_two(a, i, wf%n_v)
                  ia = index_two(i, a, wf%n_o)
                  ja = index_two(j, a, wf%n_o)
!
                  eta_ai_bj(ai, bj) = four*g_ia_jb(ia, jb) - two*g_ia_jb(ib, ja) ! 2 * L_iajb
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(g_ia_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Pack vector into doubles eta
!
      do j = 1, wf%n_o
         do b = 1, wf%n_v
!
            bj = index_two(b, j, wf%n_v)
!
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  ai = index_two(a, i, wf%n_v)
!
                  aibj = index_packed(ai, bj)
!
                  eta(wf%n_t1am + aibj, 1) = eta_ai_bj(ai, bj)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%mem%dealloc(eta_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
   end subroutine construct_eta_ccsd
!
!
end module ccsd_class
!
