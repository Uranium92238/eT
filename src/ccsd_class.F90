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
!     Initialization routine for the (singles, doubles) amplitudes
!
      procedure :: initialize_amplitudes => initialize_amplitudes_ccsd
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
      procedure :: omega_d1 => omega_d1_ccsd
!
      procedure :: omega_a2 => omega_a2_ccsd 
      procedure :: omega_b2 => omega_b2_ccsd 
      procedure :: omega_c2 => omega_c2_ccsd 
      procedure :: omega_d2 => omega_d2_ccsd 
      procedure :: omega_e2 => omega_e2_ccsd       
!
!     Ground state solver routine (helpers only, see CCS for the rest)
!
      procedure :: calc_ampeqs_norm          => calc_ampeqs_norm_ccsd
      procedure :: new_amplitudes            => new_amplitudes_ccsd
      procedure :: calc_quasi_Newton_doubles => calc_quasi_Newton_doubles_ccsd
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
      module subroutine omega_d1_ccsd(wf)
!!
!!       Omega D1 
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, March 2017
!!
!!       Omega_ai^D1 = F_ai_T1
!!
         implicit none 
!
         class(ccsd) :: wf
!
      end subroutine omega_d1_ccsd
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
      module subroutine calc_quasi_Newton_doubles_ccsd(wf,dt,n_variables)
!!
!!       Calculate quasi-Newton estimate (CCSD)
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!       Calculates the quasi-Newton estimate Δ t_i (doubbles part)
!!       and places the contribution in the dt vector (of length n_variables,
!!       with singles first, then doubles, etc. if inherited)
!!
         implicit none 
!
         class(ccsd) :: wf 
!
         integer(i15), intent(in) :: n_variables
         real(dp), dimension(n_variables, 1) :: dt
!
      end subroutine calc_quasi_Newton_doubles_ccsd
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
      call wf%initialize_amplitudes
!
!     Initialize the Fock matrix (allocate and construct given the initial amplitudes)
!
      call wf%initialize_fock_matrix
!
!     Initialize the projection vector (omega)
!
      call wf%initialize_omega
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
!!     Allocates the amplitudes, sets them to zero, calculates
!!     the number of amplitudes, and sets the doubles amplitudes
!!     to the perturbative MP2 estimate.
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
!     Calculate the number of singles and doubles amplitudes
!
      wf%n_t1am = (wf%n_o)*(wf%n_v) 
      wf%n_t2am = (wf%n_t1am)*(wf%n_t1am + 1)/2
!
!     Allocate the singles amplitudes and set to zero
!
      call allocator(wf%t1am, wf%n_v, wf%n_o)
      wf%t1am = zero
!
!     Allocate the doubles amplitudes and set to zero
!
      call allocator(wf%t2am, wf%n_t2am, 1)
      wf%t2am = zero
!
!
!     :: Initialize the doubles amplitudes to the MP2 estimate ::
!
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
   end subroutine initialize_amplitudes_ccsd
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
end module ccsd_class
