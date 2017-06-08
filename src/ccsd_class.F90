module ccsd_class
!
!!
!!           Coupled cluster singles and doubles (CCSD) class module                                 
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
!  Ancestor class module (CCS)
!
   use ccs_class
!
   implicit none 
!
!  ::::::::::::::::::::::::::::::::::::::
!  -::- Definition of the CCSD class -::-
!  ::::::::::::::::::::::::::::::::::::::
!
   type, extends(ccs) :: ccsd
!
!     Amplitude variables
!
      integer(i15) :: n_t2am = 0                    ! Number of doubles amplitudes
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
!     Initialization routine (driver is inherited)
!
      procedure :: init => init_ccsd
!
!     Initialization routine for the (singles, doubles) amplitudes,
!     and routine to set MP2 guess for the doubles amplitudes 
!
      procedure :: initialize_amplitudes          => initialize_amplitudes_ccsd
      procedure :: construct_perturbative_doubles => construct_perturbative_doubles_ccsd
!
!     Routine to calculate the energy from the current amplitudes
!
      procedure :: calc_energy => calc_energy_ccsd 
!
!     Routine to initialize omega (allocate and set to zero)
!
      procedure :: initialize_omega => initialize_omega_ccsd
!
!     Routines to construct the projection vector (omega)
!
      procedure :: construct_omega => construct_omega_ccsd
!
!     Helper routines for construct_omega
!
      procedure :: omega_a1 => omega_a1_ccsd 
      procedure :: omega_b1 => omega_b1_ccsd 
      procedure :: omega_c1 => omega_c1_ccsd
!
      procedure :: omega_a2 => omega_a2_ccsd 
      procedure :: omega_b2 => omega_b2_ccsd 
      procedure :: omega_c2 => omega_c2_ccsd 
      procedure :: omega_d2 => omega_d2_ccsd 
      procedure :: omega_e2 => omega_e2_ccsd   
!
!     Ground state solver routine (helpers only, see CCS for the rest)
!
      procedure :: initialize_ground_state   => initialize_ground_state_ccsd
      procedure :: calc_ampeqs_norm          => calc_ampeqs_norm_ccsd
      procedure :: new_amplitudes            => new_amplitudes_ccsd
      procedure :: calc_quasi_Newton_doubles => calc_quasi_Newton_doubles_ccsd
!
!     Routine to save and read the amplitudes (to/from disk)
!
      procedure :: save_amplitudes => save_amplitudes_ccsd
      procedure :: read_amplitudes => read_amplitudes_ccsd
      procedure :: read_double_amplitudes => read_double_amplitudes_ccsd
!
!     Jacobian transformation routine 
!
      procedure :: calculate_orbital_differences => calculate_orbital_differences_ccsd
      procedure :: jacobian_ccsd_transformation  => jacobian_ccsd_transformation_ccsd
      procedure :: transform_trial_vectors       => transform_trial_vectors_ccsd
!
!     Helper routines for Jacobian transformation 
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
      procedure :: jacobi_test => jacobi_test_ccsd
!
!     Routines to destroy amplitudes and omega 
!
      procedure :: destruct_amplitudes => destruct_amplitudes_ccsd
      procedure :: destruct_omega      => destruct_omega_ccsd
!
   end type ccsd
!
!  :::::::::::::::::::::::::::::::::::::::::::::::::::::
!  -::- Interface to the submodule routines of CCSD -::- 
!  :::::::::::::::::::::::::::::::::::::::::::::::::::::
!
   interface
!
!
      module subroutine initialize_omega_ccsd(wf)
!!
!!       Initialize Omega (CCSD)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!       Allocates the projection vector (omega1, omega2) and sets it
!!       to zero.
!!
         implicit none 
!
         class(ccsd) :: wf
!
      end subroutine initialize_omega_ccsd
!
!
      module subroutine construct_omega_ccsd(wf)
!!
!!       Construct Omega (CCSD)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!       Directs the construction of the projection vector < mu | exp(-T) H exp(T) | R >
!!       for the current amplitudes of the object wfn 
!!
         implicit none 
!
         class(ccsd) :: wf 
!
      end subroutine construct_omega_ccsd
!
!
      module subroutine omega_a1_ccsd(wf)
!!
!!       Omega A1 term
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!  
!!       Calculates the A1 term, 
!!  
!!       A1: sum_ckd g_adkc * u_ki^cd,
!!  
!!       and adds it to the singles projection vector (omeg1) of
!!       the wavefunction object wf.
!!
         implicit none 
!
         class(ccsd) :: wf
!
      end subroutine omega_a1_ccsd
!
!
      module subroutine omega_b1_ccsd(wf)
!!
!!       Omega B1
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!       Calculates the B1 term, 
!!
!!       B1: - sum_ckl u_kl^ac * g_kilc,
!! 
!!       and adds it to the singles projection vector (omeg1) of
!!       the wavefunction object wf
!!
         implicit none 
!
         class(ccsd) :: wf 
!
      end subroutine omega_b1_ccsd
!
!
      module subroutine omega_c1_ccsd(wf)
!!  
!!     Omega C1
!!     Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!  
!!     Calculates the C1 term of omega,
!!  
!!     C1: sum_ck F_kc*u_ai_ck,
!!  
!!     and adds it to the projection vector (omega1) of    
!!     the wavefunction object wf                           
!!  
!!     u_ai_ck = 2*t_ck_ai - t_ci_ak
!! 
         implicit none 
!
         class(ccsd) :: wf 
!
      end subroutine omega_c1_ccsd
!
!
      module subroutine omega_a2_ccsd(wf)
!!
!!     Omega A2 term: Omega A2 = g_ai_bj + sum_(cd)g_ac_bd * t_ci_dj = A2.1 + A.2.2
!!
!!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, 10 Mar 2017
!!
!!     Structure: Batching over both a and b for A2.2.
!!                t^+_ci_dj = t_ci_dj + t_di_cj
!!                t^-_ci_dj = t_ci_dj - t_di_cj
!!                g^+_ac_bd = g_ac_bd + g_bc_ad 
!!                g^-_ac_bd = g_ac_bd - g_bc_ad 
!! 
!!                omega_A2.2_ai_bj = 1/4*(g^+_ac_bd*t^+_ci_dj + g^-_ac_bd*t^-_ci_dj) = omega_A2.2_bj_ai
!!                omega_A2.2_aj_bi = 1/4*(g^+_ac_bd*t^+_ci_dj - g^-_ac_bd*t^-_ci_dj) = omega_A2.2_bi_aj
!!
         implicit none 
!
         class(ccsd) :: wf
!
      end subroutine omega_a2_ccsd
!
!
      module subroutine omega_b2_ccsd(wf)
!!
!!       Omega B2
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, 11 Mar 2017
!! 
!!       Omega B2 = sum_(kl) t_ak_bl*(g_kilj + sum_(cd) t_ci_dj * g_kc_ld)
!!
!!       Structure: g_kilj is constructed first and reordered as g_kl_ij. 
!!       Then the contraction over cd is performed, and the results added to g_kl_ij.
!!       t_ak_bl is then reordered as t_ab_kl and the contraction over kl is performed.
!!
         implicit none 
!
         class(ccsd) :: wf
!
      end subroutine omega_b2_ccsd
!
!
      module subroutine omega_c2_ccsd(wf)
!!
!!       Omega C2 
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2017
!!     
!!       Omega C2 = -1/2* sum_(ck)t_bk_cj*(g_ki_ac -1/2 sum_(dl)t_al_di * g_kd_lc)
!!                                  - sum_(ck) t_bk_ci (g_kj_ac-sum_(dl)t_al_dj*g_kd_lc)
!!
         implicit none 
!
         class(ccsd) :: wf
!
      end subroutine omega_c2_ccsd
!
!
      module subroutine omega_d2_ccsd(wf)
!!
!!       Omega D2 
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!       Calculates the D2 term,
!!
!!          sum_ck u_jk^bc g_aikc 
!!           - 1/2 * sum_ck u_jk^bc g_acki 
!!           + 1/4 * sum_ck u_jk^bc sum_dl L_ldkc u_il^ad,
!!
!!       where 
!!
!!          u_jk^bc = 2 * t_jk^bc - t_kj^bc,
!!          L_ldkc  = 2 * g_ldkc  - g_lckd.
!!
!!       The first, second, and third terms are referred to as D2.1, D2.2, and D2.3, 
!!       and comes out ordered as (ai,bj). All terms are added to the omega vector of the 
!!       wavefunction object wf.
!!
!!       The routine adds the terms in the following order: D2.3, D2.1, D2.2
!!
         implicit none 
!
         class(ccsd) :: wf
!
      end subroutine omega_d2_ccsd
!
!
      module subroutine omega_e2_ccsd(wf)
!!
!!       Omega E2
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!       Calculates the E2 term,
!!
!!           sum_c t_ij^ac (F_bc - sum_dkl g_ldkc u_kl^bd) 
!!           - sum_k t_ik^ab (F_kj + sum_cdl g_ldkc u_lj^dc),
!!
!!       where
!!
!!          u_kl^bc = 2 * t_kl^bc - t_lk^bc.
!!
!!       The first term is referred to as the E2.1 term, and comes out ordered as (b,jai).
!!       The second term is referred to as the E2.2 term, and comes out ordered as (aib,j).
!!
!!       Both are permuted added to the projection vector element omega2(ai,bj) of
!!       the wavefunction object wf.
!!
         implicit none 
!
         class(ccsd) :: wf
!
      end subroutine omega_e2_ccsd
!
!
      module subroutine calc_ampeqs_norm_ccsd(wf, ampeqs_norm)
!!
!!       Calculate Amplitude Equations Norm (CCSD)
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
      module subroutine new_amplitudes_ccsd(wf)
!!
!!       New Amplitudes (CCSD)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!       Directs the calculation of the quasi-Newton estimate Δ t_i, 
!!       and t_i + Δ t_i, and calls the DIIS routine to save & get 
!!       the amplitudes for the next iteration. 
!!
         implicit none 
!
         class(ccsd) :: wf 
!
      end subroutine new_amplitudes_ccsd
!
!
      module subroutine calc_quasi_Newton_doubles_ccsd(wf,dt)
!!
!!       Calculate quasi-Newton estimate (CCSD)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!       Calculates the quasi-Newton estimate Δ t_i (doubbles part)
!!       and places the contribution in the dt vector (of length n_parameters,
!!       with singles first, then doubles, etc. if inherited)
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
      module subroutine initialize_ground_state_ccsd(wf)
!!
!!       Initialize Ground State (CCSD)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!       Initializes the amplitudes and the projection vector for the ground
!!       state solver. 
!!
         implicit none 
!
         class(ccsd) :: wf
!
      end subroutine initialize_ground_state_ccsd
!
!
      module subroutine calculate_orbital_differences_ccsd(wf,orbital_diff)
!!
!!       Calculate Orbital Differences (CCSD)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad May 2017
!!
!!       Calculates orbital differences
!!
!!          1) ε_i^a = ε_a - ε_i
!!          2) ε_ij^ab = ε_a + ε_b - ε_i - ε_j
!!
!!       and puts them in orbital_diff, which is a vector of length n_parameters.        
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
         class(ccsd) :: wf
!
         integer(i15), intent(in) :: first_trial, last_trial ! Which trial_vectors we are to transform
!
      end subroutine transform_trial_vectors_ccsd
!
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
      module subroutine jacobian_ccsd_a1_ccsd(wf, rho_a_i, c_a_i)
!!
!!       Jacobian CCSD A1
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017 
!!
!!       rho_ai^A1 = sum_ckdl L_lckd (u_li^ca c_dk  - t_li^cd c_ak - t_lk^ad c_ci)
!!
         implicit none 
!
         class(ccsd) :: wf
!
         real(dp), dimension(wf%n_v, wf%n_o) :: c_a_i   ! c_ai 
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
!!       rho_ai^B1 = sum_bj F_jb (2*c_ai_bj  -  c_aj_bi) 
!!              = sum_bj F_jb v_ai_bj
!!
!!       The term is added as rho_a_i(a,i) = rho_a_i(a,i) + rho_ai^A1,
!!       where c_a_i(a,i) = c_ai above. 
!!
         implicit none 
!
         class(ccsd) :: wf
!
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: c_ai_bj   ! c_aibj 
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
!!       rho_ai^C1 = - sum_bjk L_jikb c_aj_bk 
!!                 = - sum_bjk (2*g_jikb - g_kijb) c_aj_bk 
!!
!!       The term is added as rho_a_i(a,i) = rho_a_i(a,i) + rho_ai^A1,
!!       where c_ai_bj(ai,bj) = c_aibj above. 
!!
         implicit none 
!
         class(ccsd) :: wf
!
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: c_ai_bj
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
!!       rho_ai^D1 =  sum_bcj L_abjc c_bicj
!!
!!       The term is added as rho_a_i(a,i) = rho_a_i(a,i) + rho_ai^A1,
!!       where c_bi_cj(bi,cj) = c_bicj above. 
!!
         implicit none 
!
         class(ccsd) :: wf
!
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: c_bi_cj
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
!!       rho_ai_bj^A2 = sum_c g_aibc c_cj - sum_k g_aikj c_bk 
!!
         implicit none 
!
         class(ccsd) :: wf 
!
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
         real(dp), dimension(wf%n_v, wf%n_o) :: c_a_i
!
      end subroutine jacobian_ccsd_a2_ccsd
!
!
      module subroutine jacobian_ccsd_b2_ccsd(wf, rho_ai_bj, c_a_i)
!!
!!       Jacobian CCSD B2 
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!       rho_ai_bj^B2 = sum_kc (F_kc t_ij^ac c_bk + F_kc t_ik^ab c_cj)
!!
         implicit none 
!
         class(ccsd) :: wf 
!
         real(dp), dimension(wf%n_v, wf%n_o) :: c_a_i
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
!!       rho_ai_bj^C2 = sum_kcl g_ljkc (t_ki^ac c_bl + t_li^bc c_ak + t_lk^ba c_ci)
!!                    - sum_kcl L_ljkc (t_il^ab c_ck + t_ik^ac c_bl)
!!
         implicit none 
!
         class(ccsd) :: wf 
!
         real(dp), dimension(wf%n_v, wf%n_o) :: c_a_i
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
!!       rho_ai_bj^D2 = sum_kcd g_kcbd (t_ij^cd c_ak + t_kj^ad c_ci + t_ik^ca c_dj)
!!                    - sum_kcd L_kcbd (t_ik^ac c_dj + t_ij^ad c_ck)
!!
!!       Note: the code is structured so that we batch over the index b,
!!             where the integrals are made as g_kc_db = g_kcbd and held
!!             in some ordering or other throughout a given batch (i.e.,
!!             all five terms are constructed gradually in the batches).
!!
         implicit none 
!
         class(ccsd) :: wf 
!
         real(dp), dimension(wf%n_v, wf%n_o) :: c_a_i
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
!!       rho_ai_bj^E2 = 2 sum_dlck t_bj,dl * L_kc,ld * c_ai,ck 
!!
         implicit none 
!
         class(ccsd) :: wf 
!
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
!
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: c_ai_ck
!
      end subroutine jacobian_ccsd_e2_ccsd
!
!
      module subroutine jacobian_ccsd_f2_ccsd(wf, rho_ai_bj, c_ai_bj)
!!
!!       Jacobian CCSD F2 
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!       rho_ai_bj^F2 =   - sum_ckld t_ai,ck * L_kc,ld * c_bl,dj 
!!                        - sum_ckdl t_ai,dj * L_kc,ld * c_bl,ck
!!                        - sum_ckdl t_ai_bl * L_kc,ld * c_ck,dj
!!
!!       L_kc,ld = 2*g_kc,ld - g_kd,lc = 2*g_kc_ld(kc,ld) - 2*g_kc_ld(kd,lc)
!!
         implicit none 
!
         class(ccsd) :: wf 
!
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
!
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: c_ai_bj
!
      end subroutine jacobian_ccsd_f2_ccsd
!
!
      module subroutine jacobian_ccsd_g2_ccsd(wf, rho_ai_bj, c_ai_bj)
!!
!!       Jacobian CCSD G2 
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!       rho_ai_bj^G2 =  P_(ij)^(ab) (- sum_ckld t_ck,dj * L_kc,ld * c_ai,bl 
!!                                    - sum_ckdl t_bl,dj * L_kc,ld * c_ai,ck
!!                                    - sum_ckdl t_ck_bl * L_kc,ld * c_ai,dj )
!!
         implicit none 
!
         class(ccsd) :: wf 
!
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: c_ai_bj
!
      end subroutine jacobian_ccsd_g2_ccsd
!
!
      module subroutine jacobian_ccsd_h2_ccsd(wf, rho_ai_bj, c_ai_bj)
!!
!!       Jacobian CCSD H2 
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!       rho_ai_bj^H2 =  P_(ij)^(ab) ( sum_ckld t_ci,ak * g_kc,ld * c_bl,dj 
!!                                   + sum_ckdl t_cj,al * g_kc,ld * c_bk,di)
!!                
!!
         implicit none 
!
         class(ccsd) :: wf 
!
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: c_ai_bj
!
      end subroutine jacobian_ccsd_h2_ccsd
!
!
      module subroutine jacobian_ccsd_i2_ccsd(wf, rho_ai_bj, c_ai_bj)
!!
!!       Jacobian CCSD I2 
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!       rho_ai_bj^I2 =  sum_c F_bc * c_ai,cj - sum_k F_jk * c_ai,bk
!!                     + sum_ck L_bj,kc * c_ai,ck 
!!                     - sum_ck ( g_kc,bj * c_ak,ci + g_ki,bc * c_ak,cj ) 
!!                
!!       Batch over c to construct  g_ki_bc
!! 
         implicit none 
!
         class(ccsd) :: wf 
!
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: c_ai_bj
!
      end subroutine jacobian_ccsd_i2_ccsd
!
!
      module subroutine jacobian_ccsd_j2_ccsd(wf, rho_ab_ij, c_ab_ij)
!!
!!       Jacobian CCSD J2 
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!       rho_ab_ij^J2 =    sum_ckld t_ci,dj * g_kc,ld * c_ak,bl 
!!                       + sum_ckdl t_ak,bl * g_kc,ld * c_ci,dj
!!                
!!
         implicit none 
!
         class(ccsd) :: wf 
!
         real(dp), dimension((wf%n_v)**2, (wf%n_o)**2) :: rho_ab_ij
         real(dp), dimension((wf%n_v)**2, (wf%n_o)**2) :: c_ab_ij
!
      end subroutine jacobian_ccsd_j2_ccsd
!
!
      module subroutine jacobian_ccsd_k2_ccsd(wf, rho_ab_ij, c_ab_ij)
!!
!!       Jacobian CCSD K2 
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!       rho_ab_ij^K2 =    sum_kl g_ki,lj * c_ak,bl 
!!                       + sum_cd g_ac,bd * c_ci,dj
!! 
!!       For the last term we batch over a and b and 
!!       add each batch to rho_ai_bj 
!!   
         implicit none 
!
         class(ccsd) :: wf 
!
         real(dp), dimension((wf%n_v)**2, (wf%n_o)**2) :: rho_ab_ij
         real(dp), dimension((wf%n_v)**2, (wf%n_o)**2) :: c_ab_ij
!
!
      end subroutine jacobian_ccsd_k2_ccsd
!
!
   end interface
!
!
contains
!
!
!  ::::::::::::::::::::::::::::::::::::::::::::
!  -::- Initialization and driver routines -::-
!  ::::::::::::::::::::::::::::::::::::::::::::
!
!
   subroutine init_ccsd(wf)
!!
!!    Initialize CCSD object
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!    Performs the following tasks:
!!
!!       1. Sets HF orbital and energy information by reading from file (read_hf_info)
!!       2. Transforms AO Cholesky vectors to MO basis and saves to file (read_transform_cholesky)
!!       3. Allocates the Fock matrix and sets it to zero
!!       4. Initializes the amplitudes (sets their initial values and associated variables)
!!
!!    Note: this routine does not calculate the energy, which is postponed until the wavefunction
!!    is passed to the ground-state solver.
!!
      implicit none 
!
      class(ccsd) :: wf
!
!     Set model name 
!
      wf%name = 'CCSD'
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
!     Initialize (singles and doubles) amplitudes
!
      wf%n_t1am = (wf%n_o)*(wf%n_v) 
      wf%n_t2am = (wf%n_t1am)*(wf%n_t1am + 1)/2 
!
      wf%n_parameters = wf%n_t1am + wf%n_t2am
!
!     Initialize the Fock matrix (allocate and construct given the initial amplitudes)
!
      if (.not. allocated(wf%t1am)) call allocator(wf%t1am, wf%n_v, wf%n_o)
      wf%t1am = zero
!
      call wf%initialize_fock_matrix
!
   end subroutine init_ccsd
!
!
!  :::::::::::::::::::::::::::::::::::::::::
!  -::- Class subroutines and functions -::- 
!  :::::::::::::::::::::::::::::::::::::::::
!
   subroutine initialize_amplitudes_ccsd(wf)
!!
!!     Initialize Amplitudes (CCSD)
!!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!     Allocates the amplitudes, sets them to zero, and calculates
!!     the number of amplitudes.
!!
      implicit none 
!
      class(ccsd) :: wf
!
!     Allocate the doubles amplitudes and set to zero
!
      if (.not. allocated(wf%t2am)) call allocator(wf%t2am, wf%n_t2am, 1)
      wf%t2am = zero
!
   end subroutine initialize_amplitudes_ccsd
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
      real(dp), dimension(:,:), allocatable :: L_ia_J
      real(dp), dimension(:,:), allocatable :: g_ia_jb
!
      integer(i15) :: i = 0, j = 0, a = 0, b = 0
      integer(i15) :: ai = 0, bj = 0, ia = 0, jb = 0, aibj = 0 
!
!     Allocate L_ia_J and g_ia_jb
!
      call allocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
      call allocator(g_ia_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      L_ia_J = zero
      g_ia_jb = zero
!
!     Get the Cholesky IA vector 
!
      call wf%get_cholesky_ia(L_ia_J)
!
!     Calculate g_ia_jb = g_iajb
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
                                                            wf%fock_diagonal(i, 1) - &
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
      call deallocator(L_ia_J, (wf%n_o)*(wf%n_v), (wf%n_J))
      call deallocator(g_ia_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) 
!
   end subroutine construct_perturbative_doubles_ccsd
!
!
   subroutine calc_energy_ccsd(wf)
!!
!!     Calculate Energy (CCSD)
!!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!     Calculates the CCSD energy for the wavefunction's current amplitudes.
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
      call deallocator(g_ia_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
   end subroutine calc_energy_ccsd
!
!
   subroutine destruct_amplitudes_ccsd(wf)
!
      implicit none
!
      class(ccsd) :: wf
!
      if (allocated(wf%t2am)) then
         call deallocator(wf%t2am, wf%n_t2am, 1)
      endif
!
   end subroutine destruct_amplitudes_ccsd
!
!
   subroutine destruct_omega_ccsd(wf)
!
      implicit none
!
      class(ccsd) :: wf
!
      if (allocated(wf%omega1)) then
         call deallocator(wf%omega1, wf%n_v, wf%n_o)
      endif
      if (allocated(wf%omega2)) then
         call deallocator(wf%omega2, wf%n_t2am, 1)
      endif
!
   end subroutine destruct_omega_ccsd
!
!
   subroutine save_amplitudes_ccsd(wf)
!!
!!    Save Amplitudes (CCSD)
!!    Written by Sarai D. Folkestad and Eirik F. Kjøsntad, May 2017
!!
!!    Store the amplitudes to disk (T1AM, T2AM)
!!
      implicit none 
!
      class(ccsd) :: wf
!
      integer(i15) :: unit_t1am = -1
      integer(i15) :: unit_t2am = -1
!
!     Open amplitude files
!
      call generate_unit_identifier(unit_t1am)
      call generate_unit_identifier(unit_t2am)
!
      open(unit_t1am, file='t1am', status='unknown', form='unformatted')
      open(unit_t2am, file='t2am', status='unknown', form='unformatted')
!
      rewind(unit_t1am)
      rewind(unit_t2am)
!
!     Write amplitudes to files
!
      write(unit_t1am) wf%t1am 
      write(unit_t2am) wf%t2am
!
!     Close amplitude files
!
      close(unit_t1am)
      close(unit_t2am)
!
   end subroutine save_amplitudes_ccsd
!
   subroutine read_amplitudes_ccsd(wf)
!!
!!    Read Amplitudes (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjøsntad, May 2017
!!
!!    Reads the amplitudes from disk (T1AM)
!!
      implicit none 
!
      class(ccsd) :: wf
!
      call wf%read_single_amplitudes
      call wf%read_double_amplitudes
!
   end subroutine read_amplitudes_ccsd
!
   subroutine read_double_amplitudes_ccsd(wf)
!!
!!    Read Amplitudes (CCSD)
!!    Written by Sarai D. Folkestad and Eirik F. Kjøsntad, May 2017
!!
!!    Reads the amplitudes from disk (T1AM, T2AM)
!!
      implicit none 
!
      class(ccsd) :: wf
!
      integer(i15) :: unit_t2am = -1 
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
         wf%t2am = zero
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
   end subroutine read_double_amplitudes_ccsd
!
!
!
!
   subroutine jacobi_test_ccsd(wf)
!
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
      call wf%initialize_amplitudes
      call wf%read_double_amplitudes
!
      call allocator(r1am, wf%n_v, wf%n_o)
      call allocator(r2am, wf%n_t2am, 1)
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
      call wf%jacobian_ccsd_transformation(r1am,r2am)
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
      do j = 1, 15
         write(unit_output,*) j, r2am(j,1)
      enddo
!
!     We wish to calculate A_mu,nu = d(omega)_mu / dt_nu, 
!     in two different ways:
!
!        1. By transforming c_tau = delta_tau,nu with A. Then (A c)_mu = sum_tau A_mu,tau c_tau = A_mu,nu
!        2. By calculating A_mu,nu = (omega(t+Dt_nu)_mu - omega(t)_mu)/Dt_nu.
!
!     :: First approach: transform by A :: 
!
      call allocator(c_a_i, wf%n_v, wf%n_o)
      call allocator(c_aibj, wf%n_t2am, 1)
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
end module ccsd_class
!
