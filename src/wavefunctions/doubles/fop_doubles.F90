!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2019 the authors of eT
!
!  eT is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  eT is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program. If not, see <https://www.gnu.org/licenses/>.
!
!
submodule (doubles_class) fop_doubles
!
!!
!!    First order properties submodule
!!
!!    Routines for construction of the right-hand-side, η^X, and left-hand-side, ξ^X
!!    vectors and the left-hand-side (ρ^L) and right-hand-side (ρ^R) transition densities
!!    for transition moments.
!!
!!    Equation-of-motion (EOM):
!!
!!    (Following Koch, H., Kobayashi, R., Sanches de Merás, A., and Jørgensen, P.,
!!    J. Chem. Phys. 100, 4393 (1994))
!!
!!       η_μ^X,EOM =  < Λ | [X, τ_μ] | CC > + (< Λ | τ_μ X | CC >  - tbar_μ < Λ | X | CC > )
!!                 = η^{X,0} + η^{X,corr}
!!
!!    Where the last two terms are called the EOM-corrections and the first term also 
!!    appears in LR-CC.
!!
!!    The left-hand-side vector is the same in EOM-CC and LR-CC:   
!!
!!       ξ^X_μ = < μ | exp(-T) X exp(T)| R >
!!
!!    The transition density matrices are construct as follows:
!!
!!       ρ^L_pq = < k | E_pq | CC >
!!       ρ^R_pq = < Λ | E_pq | k >
!!
!!    where |k> and <k| are the eigenvectors of the Jacobian with amplitudes R_μ, L_μ
!!
!!       | k > = sum_μ (τ_μ | CC > R_{k,μ} - tbar_μ | CC > R_{k,μ})
!!       < k | = sum_μ L_{k,μ} < μ | e^-T
!!
!
   implicit none
!
!
contains
!
!
   module subroutine construct_left_transition_density_doubles(wf, state)
!!
!!    Construct left one-electron transition density for the state k
!!    Written by Alexander Paul, June 2019
!!
!!          ρ^L_pq = < k | E_pq | CC >
!!
!!    where <k| is the left eigenvector of the Jacobian
!!    with amplitudes L_μ
!!
!!          < k | = sum_μ L_{k,μ} < μ | e^-T
!!
      implicit none
!
      class(doubles) :: wf
!
      integer, intent(in) :: state
!
      real(dp), dimension(:), allocatable :: L_k
!
      real(dp), dimension(:,:), allocatable :: L_ai
      real(dp), dimension(:,:,:,:), allocatable :: L_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: t_aibj
!
      type(timings) :: L_TDM_timer
!
      call L_TDM_timer%turn_on()
!
      L_TDM_timer = timings('Left transition density')
!
      call mem%alloc(L_k, wf%n_es_amplitudes)
      call wf%read_excited_state(L_k, state, 'left')
!
      call zero_array(wf%left_transition_density, wf%n_mo**2)
!
!     Allocate the singles part of the excitation vector
!
      call mem%alloc(L_ai, wf%n_v, wf%n_o)
!
      call dcopy(wf%n_t1, L_k, 1, L_ai, 1)
!
      call wf%gs_one_el_density_ccs_vo(wf%left_transition_density, L_ai)
!
      call mem%alloc(t_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_aibj, (wf%n_v)*(wf%n_o))
!
      call wf%gs_one_el_density_doubles_ov(wf%left_transition_density, L_ai, t_aibj)
!
      call mem%dealloc(L_ai, wf%n_v, wf%n_o)
!
!     Allocate and unpack doubles part of the excitation vector
!
      call mem%alloc(L_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call squareup(L_k(wf%n_t1 + 1 : wf%n_es_amplitudes), L_aibj, wf%n_t1)
!
      call mem%dealloc(L_k, wf%n_es_amplitudes)
!
      call wf%gs_one_el_density_doubles_oo(wf%left_transition_density, L_aibj, t_aibj)
      call wf%gs_one_el_density_doubles_vv(wf%left_transition_density, L_aibj, t_aibj)
!
      call mem%dealloc(t_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(L_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call L_TDM_timer%turn_off()
!      
   end subroutine construct_left_transition_density_doubles
!
!
   module subroutine construct_right_transition_density_doubles(wf, state)
!!
!!    Construct right one-electron transition density for the state k
!!    Written by Alexander Paul, June 2019
!!
!!          ρ^R_pq = < Λ | E_pq | k >
!!
!!    where |k> is the right eigenvector of the Jacobian
!!    with amplitudes R_μ
!!
!!          | k > = sum_μ (τ_μ | CC > R_{k,μ} - tbar_μ | CC > R_{k,μ}) 
!!
      implicit none
!
      class(doubles) :: wf
!
      integer, intent(in) :: state
!
      real(dp), dimension(:), allocatable :: R_k
!
      real(dp), dimension(:,:), allocatable :: R_ai
      real(dp), dimension(:,:,:,:), allocatable :: R_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: tbar_aibj
!
      real(dp), dimension(:,:), allocatable :: rho_corr
!
      real(dp) :: scaling_factor
      real(dp) :: ddot
!
      integer :: a, i
!
      type(timings) :: R_TDM_timer
!
      R_TDM_timer = timings('Right transition density')
!
      call R_TDM_timer%turn_on()
!
      call zero_array(wf%right_transition_density, (wf%n_mo)**2)
!
      call mem%alloc(R_k, wf%n_es_amplitudes)
      call wf%read_excited_state(R_k, state, 'right')
!
      call mem%alloc(R_ai, wf%n_v, wf%n_o)
!
      call dcopy(wf%n_t1, R_k, 1, R_ai, 1)
!
      scaling_factor = -one*ddot(wf%n_v*wf%n_o, R_ai, 1, wf%t1bar, 1)
!
      call wf%right_transition_density_ccs_oo(wf%t1bar, R_ai)
      call wf%right_transition_density_ccs_ov(R_ai)
      call wf%right_transition_density_ccs_vv(wf%t1bar, R_ai)
!
      call mem%alloc(tbar_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2bar, tbar_aibj, (wf%n_v)*(wf%n_o))
!
      call wf%right_transition_density_doubles_ov(tbar_aibj, R_ai)
      call wf%right_transition_density_doubles_vo(tbar_aibj, R_ai)
!
      call mem%dealloc(R_ai, wf%n_v, wf%n_o)
!
!     Allocate and unpack doubles part of the excitation vector
!
      call mem%alloc(R_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call squareup(R_k(wf%n_t1 + 1 : wf%n_es_amplitudes), R_aibj, wf%n_t1)
!
      call mem%dealloc(R_k, wf%n_es_amplitudes)
!
!     Scale the doubles vector by 1 + δ_ai,bj, i.e.
!     redefine to c_ckdl = c_ckdl (1 + δ_ck,dl)
!
!$omp parallel do schedule(static) private(a,i)
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            R_aibj(a,i,a,i) = two*R_aibj(a,i,a,i)
!
         enddo
      enddo
!$omp end parallel do
!
      scaling_factor = scaling_factor &
                        - half * ddot((wf%n_v)**2*(wf%n_o)**2, R_aibj, 1, tbar_aibj, 1)
!
      call wf%gs_one_el_density_doubles_oo(wf%right_transition_density, tbar_aibj, R_aibj)
      call wf%gs_one_el_density_doubles_vv(wf%right_transition_density, tbar_aibj, R_aibj)
!
      call mem%dealloc(tbar_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%gs_one_el_density_doubles_ov(wf%right_transition_density, wf%t1bar, R_aibj)
!
      call mem%dealloc(R_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Right transition density, contribution from the ground state density
!     ρ^R_pq -= sum_μν R_{k,μ}tbar_μ tbar_ν < ν |e^-T E_pq e^T| HF >
!            -= sum_μ R_{k,μ} tbar_μ (D_GS - D_HF)
!            -= sum_μ R_{k,μ} tbar_μ ρ_corr
!
      call mem%alloc(rho_corr, wf%n_mo, wf%n_mo)
      call dcopy(wf%n_mo**2, wf%density, 1, rho_corr, 1)
!
!$omp parallel do private(i)
      do i = 1, wf%n_o
!
         rho_corr(i,i) = rho_corr(i,i) - two  
!
      enddo
!$omp end parallel do
!
      call daxpy(wf%n_mo**2, scaling_factor, rho_corr, 1, wf%right_transition_density, 1)
!
      call mem%dealloc(rho_corr, wf%n_mo, wf%n_mo)
!
      call R_TDM_timer%turn_off()
!
   end subroutine construct_right_transition_density_doubles
!
!
   module subroutine right_transition_density_doubles_ov_doubles(wf, tbar_aibj, R_ai)
!!
!!    Right transition density ov contribution (CCSD)
!!    from the singles part of the excitation vector
!!    Written by Alexander Paul, June 2019
!!
!!    ρ^R_kc += sum_abij R^a_i tbar^ab_ij (2t^bc_jk - t^bc_kj)
!!             -sum_abij tbar^ab_ij (R^b_k t^ac_ij + R^c_j t^ab_ik)
!!      
      implicit none
!
      class(doubles) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: tbar_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: R_ai
!
      real(dp), dimension(:,:,:,:), allocatable :: u_bjck, t_aick, t_aicj
      real(dp), dimension(:,:), allocatable :: X_bj, X_ck, X_jk, X_bc
!
      integer :: k, c
!
!     :: Term 1: sum_abij R^a_i tbar^ab_ij (2t^bc_jk - t^bc_kj) ::
!
      call mem%alloc(X_bj, wf%n_v, wf%n_o)
!
      call dgemv('T',            &
                  wf%n_o*wf%n_v, &
                  wf%n_o*wf%n_v, &
                  one,           &
                  tbar_aibj,     & ! tbar_bj_ai
                  wf%n_o*wf%n_v, &
                  R_ai,          & ! R_ai
                  1,             &
                  zero,          &
                  X_bj,          &
                  1)
!
      call mem%alloc(t_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_aick, wf%n_v*wf%n_o)
!
      call mem%alloc(u_bjck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dcopy((wf%n_v)**2*(wf%n_o)**2, t_aick, 1, u_bjck, 1)
      call dscal((wf%n_v)**2*(wf%n_o)**2, two, u_bjck, 1)
      call add_1432_to_1234(-one, t_aick, u_bjck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(t_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(X_ck, wf%n_v, wf%n_o)
!
      call dgemv('T',            &
                  wf%n_o*wf%n_v, &
                  wf%n_o*wf%n_v, &
                  one,           &
                  u_bjck,        & ! u_ck_bj
                  wf%n_o*wf%n_v, &
                  X_bj,          & ! X_bj
                  1,             &
                  zero,          &
                  X_ck,          &
                  1)
!
      call mem%dealloc(X_bj, wf%n_v, wf%n_o)
      call mem%dealloc(u_bjck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     :: Term 2: -sum_abij tbar^ab_ij (R^b_k t^ac_ij + R^c_j t^ab_ik) ::
!
      call mem%alloc(t_aicj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_aicj, wf%n_v*wf%n_o)
!
      call mem%alloc(X_bc, wf%n_v, wf%n_v)
!
!     X_bc = sum_aij tbar_aibj t_aicj = sum_abi tbar_bjai t_cjai
!
      call dgemm('N', 'T',          &
                  wf%n_v,           &
                  wf%n_v,           &
                  wf%n_v*wf%n_o**2, &
                  one,              &
                  tbar_aibj,        & ! tbar_b_jai
                  wf%n_v,           &
                  t_aicj,           & ! t_c_jai
                  wf%n_v,           &
                  zero,             &
                  X_bc,             &
                  wf%n_v)
!
!     X_jk = sum_aib tbar_aibj t_aibk
!
      call mem%alloc(X_jk, wf%n_o, wf%n_o)
!
      call dgemm('T', 'N',          &
                  wf%n_o,           &
                  wf%n_o,           &
                  wf%n_o*wf%n_v**2, &
                  one,              &
                  tbar_aibj,        & ! tbar_aib_j
                  wf%n_o*wf%n_v**2, &
                  t_aicj,           & ! t_aib_k
                  wf%n_o*wf%n_v**2, &
                  zero,             &
                  X_jk,             &
                  wf%n_o)
!
      call mem%dealloc(t_aicj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T','N',  &
                  wf%n_v,  &
                  wf%n_o,  &
                  wf%n_v,  &
                  -one,    &
                  X_bc,    & ! X_c_b
                  wf%n_v,  &
                  R_ai,    & ! R_b_k
                  wf%n_v,  &
                  one,     &
                  X_ck,    &
                  wf%n_v)
!
      call mem%dealloc(X_bc, wf%n_v, wf%n_v)
!
      call dgemm('N','N',  &
                  wf%n_v,  &
                  wf%n_o,  &
                  wf%n_o,  &
                  -one,    &
                  R_ai,    & ! R_c_j
                  wf%n_v,  &
                  X_jk,    & ! X_j_k
                  wf%n_o,  &
                  one,     &
                  X_ck,    &
                  wf%n_v)
!
      call mem%dealloc(X_jk, wf%n_o, wf%n_o)
!
!$omp parallel do private(c, k)
      do c = 1, wf%n_v
         do k = 1, wf%n_o
!
            wf%right_transition_density(k, wf%n_o + c) = X_ck(c, k) &
                        + wf%right_transition_density(k, wf%n_o + c)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(X_ck, wf%n_v, wf%n_o)
!
   end subroutine right_transition_density_doubles_ov_doubles
!
!
   module subroutine right_transition_density_doubles_vo_doubles(wf, tbar_aibj, R_ai)
!!
!!    Right transition density ov contribution (CCSD)
!!    Written by Alexander Paul, June 2019
!!
!!    ρ^R_bj += sum_ai R^a_i tbar^ab_ij
!!      
      implicit none
!
      class(doubles) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: tbar_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: R_ai
!
      real(dp), dimension(:,:), allocatable :: rho_vo
!
      integer :: i, a
!
      call mem%alloc(rho_vo, wf%n_v, wf%n_o)
!
      call dgemv('T',            &
                  wf%n_o*wf%n_v, &
                  wf%n_o*wf%n_v, &
                  one,           &
                  tbar_aibj,     & ! tbar_bj_ai
                  wf%n_o*wf%n_v, &
                  R_ai,          & ! R_ai
                  1,             &
                  zero,          &
                  rho_vo,        &
                  1)
!
!$omp parallel do private(a, i)
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            wf%right_transition_density(wf%n_o + a, i) = rho_vo(a, i) &
                        + wf%right_transition_density(wf%n_o + a, i)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(rho_vo, wf%n_v, wf%n_o)
!
   end subroutine right_transition_density_doubles_vo_doubles
!
   module subroutine construct_eom_etaX_doubles(wf, X, csiX, etaX)
!!
!!    Construct EOM etaX
!!    Written by Sarai D. Folkestad, May 2019
!!
!!    Constructs the EOM effective etaX vector, adding the EOM
!!    correction to etaX. 
!!
      implicit none
!
      class(doubles), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: csiX
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: etaX
!
      call wf%construct_etaX(X, etaX)
!
      call wf%etaX_eom_a(etaX, csiX)
!
      call wf%etaX_eom_doubles_a1(X, etaX(1:wf%n_t1))
!
   end subroutine construct_eom_etaX_doubles
!
!
   module subroutine construct_etaX_doubles(wf, X, etaX)
!!
!!    Construct η^X
!!    Written by Josefine H. Andersen, 2019
!!
!!    Adapted by Sarai D. Folekstad, Apr 2019
!!
!!    Constructs left-hand-side vector etaX:
!!
!!       η^X_μ = < Λ | [X, τ_μ] | CC >
!!
      implicit none
!
      class(doubles), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!      
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: etaX
!
      real(dp), dimension(:,:), allocatable :: etaX_ai
      real(dp), dimension(:,:,:,:), allocatable :: etaX_aibj
!
      integer :: a, i, ai
!
      call zero_array(etaX, wf%n_es_amplitudes)
!
!     etaX_ai:
!
      call mem%alloc(etaX_ai, wf%n_v, wf%n_o)
      call zero_array(etaX_ai, (wf%n_o*wf%n_v))
!
      call wf%etaX_ccs_a1(X, etaX_ai)
      call wf%etaX_ccs_b1(X, etaX_ai)
!
      call wf%etaX_doubles_a1(X, etaX_ai)
!
!$omp parallel do private (a, i, ai)
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            ai = wf%n_v*(i - 1) + a
!
            etaX(ai) = etaX_ai(a,i) 
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(etaX_ai, wf%n_v, wf%n_o)
!
!     etaX_aibj:
!
      call mem%alloc(etaX_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call zero_array(etaX_aibj, (wf%n_o*wf%n_v)**2)
!
      call wf%etaX_doubles_a2(X, etaX_aibj)
      call wf%etaX_doubles_b2(X, etaX_aibj)
!
      call symmetrize_and_add_to_packed(etaX(wf%n_t1 + 1 : wf%n_es_amplitudes), etaX_aibj, (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(etaX_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine construct_etaX_doubles
!
!
   module subroutine etaX_doubles_a1_doubles(wf, X, etaX_ai)
!!
!!    etaX CCSD A1
!!    Written by Josefine H. Andersen, Feb 2019
!!
!!    Adapted by Sarai D. Folkestad, Apr 2019
!!
!!       A1 = - sum_ckdl (tb_ckal X_id t_ckdl + tb_ckdi X_la t_ckdl)
!!
      implicit none
!
      class(doubles), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!      
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: etaX_ai
!        
      real(dp), dimension(:,:), allocatable :: X_id ! X_la
!
      real(dp), dimension(:,:,:,:), allocatable :: tb_ckal ! tb_ckdi
      real(dp), dimension(:,:,:,:), allocatable :: t_lckd
!
      real(dp), dimension(:,:), allocatable :: I_ad    ! intermediate, first term
      real(dp), dimension(:,:), allocatable :: I_li    ! intermediate, second term
!
      integer :: d, i
!
      call mem%alloc(X_id, wf%n_o, wf%n_v)
!
!$omp parallel do private(d, i)
      do d = 1, wf%n_v
         do i = 1,wf%n_o
!
            X_id(i, d) = X(i, d + wf%n_o)
!
         enddo
      enddo
!$omp end parallel do 
!
!     Squareup multipliers
!
      call mem%alloc(tb_ckal, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2bar, tb_ckal, (wf%n_v)*(wf%n_o))
!      
!     Read amplitudes and order as t_lck_d = t_kl^cd
!
      call mem%alloc(t_lckd, (wf%n_o), (wf%n_v), (wf%n_o), wf%n_v)
      call squareup_and_sort_1234_to_2341(wf%t2, t_lckd, (wf%n_v), (wf%n_o), (wf%n_v), (wf%n_o))
!      
!     :: First term: - sum_ckdl tb_ckal X_id t_ckdl
!
!     I_a_d = sum_ckl tb_a_lck t_lck_d = sum_ckl tb_ckal t_kl^cd
!  
      call mem%alloc(I_ad, wf%n_v, wf%n_v)
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  tb_ckal,              & 
                  wf%n_v,               &
                  t_lckd,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  zero,                 &
                  I_ad,                 &
                  wf%n_v)
!
!     Add   - sum_ckdl tb_ckal X_id t_kl^cd
!         = - sum_d I_a_d X_id 
!         = - sum_d I_a_d X_i_a^T(d,i)
!
      call dgemm('N','T',     &
                  wf%n_v,     &
                  wf%n_o,     &
                  wf%n_v,     &
                  -one,       &
                  I_ad,       &
                  wf%n_v,     &
                  X_id,       &
                  wf%n_o,     &
                  one,        &
                  etaX_ai,    &
                  wf%n_v)          
!
      call mem%dealloc(I_ad, wf%n_v, wf%n_v)
!
!     :: Second term: sum_ckdl tb_ckdi X_la t_ckdl
!
!     X_l_i = sum_ckd t_l_ckd tb_ckd_i  = sum_ckd tb_ckdi t_kl^cd
!
      call mem%alloc(I_li, wf%n_o, wf%n_o)
!
      call dgemm('N','N',               &
                  wf%n_o,               &
                  wf%n_o,               &
                  (wf%n_o)*(wf%n_v)**2, &
                  one,                  &
                  t_lckd,               & ! t_l_ckd
                  (wf%n_o),             &
                  tb_ckal,              & ! tbar_ckd_i
                  (wf%n_o)*(wf%n_v)**2, &
                  zero,                 &
                  I_li,                 &
                  wf%n_o)
!
      call mem%dealloc(t_lckd, (wf%n_o), (wf%n_v), (wf%n_o), wf%n_v)
      call mem%dealloc(tb_ckal, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Add - sum_ckdl b_ckdi X_la t_kl^cd = - sum_l X_la I_l_i = - sum_l X_i_a^T(a,l) I_l_i(l,i)
!
      call dgemm('T','N',     &
                  wf%n_v,     &
                  wf%n_o,     &
                  wf%n_o,     &
                  -one,       &
                  X_id,       &
                  wf%n_o,     &
                  I_li,       &
                  wf%n_o,     &
                  one,        &
                  etaX_ai,    &
                  wf%n_v)
!
      call mem%dealloc(I_li, wf%n_o, wf%n_o)
      call mem%dealloc(X_id, wf%n_o, wf%n_v)
!
   end subroutine etaX_doubles_a1_doubles
!
   module subroutine etaX_doubles_a2_doubles(wf, X, etaX_aibj)
!!
!!    etaX CCSD A2
!!    Written by Josefine H. Andersen, Feb 2019
!!
!!    Adapted by Sarai D. Folkestad, Apr 2019
!!   
!!    Constructs the A2 term of the etaX vector
!! 
!!       A2 = 2 X_jb tb_ai - X_ib tb_aj
!!
      implicit none
!
      class(doubles), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout) :: etaX_aibj
!
      integer :: i, a, j, b
!  
!$omp parallel do private(a, i, b, j)
       do j = 1, wf%n_o
         do b = 1, wf%n_v
!  
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  etaX_aibj(a, i, b, j) = etaX_aibj(a, i, b, j) - (X(i, b + wf%n_o))*(wf%t1bar(a, j)) &
                                       + two*(X(j, b + wf%n_o))*(wf%t1bar(a, i))
!
               enddo
            enddo
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine etaX_doubles_a2_doubles
!
!
   module subroutine etaX_doubles_b2_doubles(wf, X, etaX_aibj)
!!
!!    etaX CCSD B2
!!    Written by Josefine H. Andersen, Feb 2019
!!
!!    Adapted by Sarai D. Folkestad, Apr 2019
!!   
!!    Constructs the B2 term of the etaX vector
!! 
!!       B2 = sum_c tb_aicj X_cb - sum_k tb_aibk X_jk
!!
      implicit none
!
      class(doubles), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout) :: etaX_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: etaX_aijb
!
      real(dp), dimension(:,:,:,:), allocatable :: tb_aibj
      real(dp), dimension(:,:,:,:), allocatable :: tb_aijc
!
      real(dp), dimension(:,:), allocatable :: X_cb
      real(dp), dimension(:,:), allocatable :: X_jk
!
      integer :: b, c, k, j
!
!     Get and squareup multipliers
!
      call mem%alloc(tb_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2bar, tb_aibj, (wf%n_v)*(wf%n_o))
!
!     Reorder multipiers to tb_aijc
!
      call mem%alloc(tb_aijc, (wf%n_v), (wf%n_o), (wf%n_o), (wf%n_v))
      call sort_1234_to_1243(tb_aibj, tb_aijc, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     :: First term: sum_c tb_aicj X_cb 
!
      call mem%alloc(X_cb, wf%n_v, wf%n_v) 
!
!$omp parallel do private(b, c)
      do b = 1, wf%n_v
         do c = 1, wf%n_v
!
            X_cb(c,b) = X(c + wf%n_o, b + wf%n_o)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%alloc(etaX_aijb, (wf%n_v), (wf%n_o), (wf%n_o), (wf%n_v))
!
      call dgemm('N','N',               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_v,               &
                  wf%n_v,               &
                  one,                  &
                  tb_aijc,              &
                  (wf%n_v)*(wf%n_o)**2, &
                  X_cb,                 &
                  wf%n_v,               &
                  zero,                 &
                  etaX_aijb,            &
                  (wf%n_v)*(wf%n_o)**2)
!
      call mem%dealloc(tb_aijc, (wf%n_v), (wf%n_o), (wf%n_o), (wf%n_v))      
      call mem%dealloc(X_cb, wf%n_v, wf%n_v)
!
      call add_1243_to_1234(one, etaX_aijb, etaX_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(etaX_aijb, (wf%n_v), (wf%n_o), (wf%n_o), (wf%n_v))
!
!     :: Second term: -sum_k tb_aick X_jk
!
      call mem%alloc(X_jk, wf%n_o, wf%n_o)
!
!$omp parallel do private(j, k)
      do k = 1, wf%n_o
         do j = 1, wf%n_o
!
            X_jk(j,k) = X(j,k)
!
         enddo
      enddo
!$omp end parallel do
!
      call dgemm('N','T',                 &
                  (wf%n_o)*(wf%n_v)**2,   &
                  wf%n_o,                 &
                  wf%n_o,                 &
                  -one,                   &
                  tb_aibj,                &
                  (wf%n_o)*(wf%n_v)**2,   &
                  X_jk,                   &
                  wf%n_o,                 &
                  one,                    &
                  etaX_aibj,              &
                  (wf%n_o)*(wf%n_v)**2)
!          
      call mem%dealloc(tb_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X_jk, wf%n_o, wf%n_o)
!
   end subroutine etaX_doubles_b2_doubles
!
!
   module subroutine construct_csiX_doubles(wf, X, csiX)
!!
!!    Construct csiX
!!    Written by Josefine H. Andersen, 2019
!!
!!    Adapted by Sarai D. Folkestad
!!
!!    Constructs ξ^X_μ :
!!
!!       ξ^X_μ = < μ | exp(-T) X exp(T)| R >
!!
      implicit none
!
      class(doubles), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!      
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: csiX
!
      real(dp), dimension(:,:), allocatable :: csiX_ai
      real(dp), dimension(:,:,:,:), allocatable :: csiX_aibj
!
      integer :: a, i, ai
!
      call zero_array(csiX, wf%n_es_amplitudes)
!
!     csiX_ai
!
      call mem%alloc(csiX_ai, wf%n_v, wf%n_o)
      call zero_array(csiX_ai, (wf%n_o*wf%n_v))
!
      call wf%csiX_ccs_a1(X, csiX_ai)
      call wf%csiX_doubles_a1(X, csiX_ai)
!
!$omp parallel do private (a, i, ai)
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            ai = wf%n_v*(i - 1) + a
!        
            csiX(ai) = csiX_ai(a,i) 
!
         enddo
      enddo
!$omp end parallel do
!
!     csiX_aibj
!      
      call mem%alloc(csiX_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call zero_array(csiX_aibj, (wf%n_o*wf%n_v)**2)
!
      call wf%csiX_doubles_a2(X, csiX_aibj)
!
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            csiX_aibj(a, i, a, i) = half*csiX_aibj(a, i, a, i)
!
         enddo
      enddo
!
      call symmetrize_and_add_to_packed(csiX(wf%n_t1 + 1 : wf%n_es_amplitudes), csiX_aibj, (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(csiX_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine construct_csiX_doubles
!
!
   module subroutine csiX_doubles_a1_doubles(wf, X, csiX_ai)
!!
!!    csiX CCSD A1
!!    Written by Josefine H. Andersen, Feb 2019
!!
!!    Adapted by Sarai D. Folkestad, Apr 2010
!!
!!    Constructs the A1 term of csiX
!!
!!       A1 = sum_ck u_aick X_kc,
!!    
!!    where u_aick = 2t_ckai - t_ciak
!!
      implicit none
!
      class(doubles), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!      
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: csiX_ai
!      
      real(dp), dimension(:,:,:,:), allocatable   :: u_aick
      real(dp), dimension(:,:,:,:), allocatable   :: t_aick
!
      real(dp), dimension(:,:), allocatable   :: X_ck
!
      integer :: k, c
!
!     X_kc ordered as X_ck
!
      call mem%alloc(X_ck, wf%n_v, wf%n_o)
!
!$omp parallel do private(k, c)
      do k = 1, wf%n_o
         do c = 1, wf%n_v
!
            X_ck(c, k) = X(k, c + wf%n_o)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%alloc(t_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_aick,  wf%n_o*wf%n_v)
!
!     Form u_aick = 2 t_ai_ck - t_ak_ci
!
      call mem%alloc(u_aick, (wf%n_v), (wf%n_o), (wf%n_v), (wf%n_o))
      call zero_array(u_aick, (wf%n_o*wf%n_v)**2)
!
      call add_1432_to_1234(-one, t_aick, u_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o) 
      call daxpy((wf%n_o)**2 * (wf%n_v)**2, two, t_aick, 1, u_aick, 1)
!
      call mem%dealloc(t_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     sum_ck u_ai_ck X_kc
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  1,                 &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  u_aick,            &
                  (wf%n_o)*(wf%n_v), &
                  X_ck,              &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  csiX_ai,           &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(u_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X_ck, wf%n_v, wf%n_o)
!
   end subroutine csiX_doubles_a1_doubles
!
!
   module subroutine csiX_doubles_a2_doubles(wf, X, csiX_aibj)
!!
!!    CsiX CCSD A2
!!    Written by Josefine H. Andersen, Feb 2019
!!
!!    Adapted by Sarai D. Folekstad
!!
!!    Construct csiX A2 contribution:
!!
!!       A2 = sum_c t_aicj X_bc - sum_k t_aibk X_kj
!!
      implicit none
!
      class(doubles), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!      
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout) :: csiX_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: t_cjai
!
      real(dp), dimension(:,:), allocatable :: X_bc
      real(dp), dimension(:,:), allocatable :: X_kj
!
      integer :: b, c, j, k
!
      call mem%alloc(t_cjai, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_cjai, wf%n_v*wf%n_o)
!
!     :: First term: sum_c t_aicj X_bc
!
      call mem%alloc(X_bc, wf%n_v, wf%n_v)
!
!$omp parallel do private(b, c)
      do c = 1, wf%n_v
         do b = 1, wf%n_v
!
            X_bc(b, c) = X(b + wf%n_o, c + wf%n_o)
!
         enddo
      enddo
!$omp end parallel do
!
      call dgemm('N','N',            &
                 wf%n_v,             &
                 wf%n_v*(wf%n_o)**2, &
                 wf%n_v,             &
                 one,                &
                 X_bc,               &
                 wf%n_v,             &
                 t_cjai,             &
                 wf%n_v,             &
                 one,                &
                 csiX_aibj,          & ! csiX_bjai but it does not matter since we will symmetrize anyhow
                 wf%n_v)
!         
      call mem%dealloc(X_bc, wf%n_v, wf%n_v)
!      
!     :: Second term: -sum_k t_ai_bk X_kj
!
      call mem%alloc(X_kj, wf%n_o, wf%n_o)
!
!$omp parallel do private(b, c)
      do j = 1, wf%n_o
         do k = 1, wf%n_o
!
            X_kj(k, j) = X(k, j)
!
         enddo
      enddo
!$omp end parallel do
!      
      call dgemm('N','N',              &
                 (wf%n_o)*(wf%n_v)**2, &
                 wf%n_o,               &
                 wf%n_o,               &
                 -one,                 &
                 t_cjai,               &
                 wf%n_o*(wf%n_v)**2,   &
                 X_kj,                 &
                 wf%n_o,               &
                 one,                  &
                 csiX_aibj,            &
                 wf%n_o*(wf%n_v)**2)
!
      call mem%dealloc(X_kj, wf%n_o, wf%n_o)
      call mem%dealloc(t_cjai, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine csiX_doubles_a2_doubles
!
!
   module subroutine etaX_eom_doubles_a1_doubles(wf, X, etaX_ai)
!!
!!    etaX EOM CCSD A1
!!    Written by Josefine H. Andersen, Feb 2019
!!
!!    Adapted by Sarai D. Folkestad
!!
!!    Constructs the A1 correction term to η^X for EOM
!!
!!       A1 = sum_ck tb_aick X_ck + sum_ckdl tb_aick u_ckdl X_ld,
!!
!!    where u_ckdl = 2*t_ckdl - t_cldk
!!
      implicit none
!
      class(doubles), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: etaX_ai
!
      real(dp), dimension(:,:,:,:), allocatable :: tb_aibj
      real(dp), dimension(:,:,:,:), allocatable :: t_ckdl
      real(dp), dimension(:,:,:,:), allocatable :: u_ckdl
!
      real(dp), dimension(:,:,:,:), allocatable :: I_aidl
!
      real(dp), dimension(:,:), allocatable :: X_ck
      real(dp), dimension(:,:), allocatable :: X_dl
!
      integer :: c, k, d, l
!
!     :: First term: sum_ck tb_aick X_ck
!
      call mem%alloc(tb_aibj, wf%n_v, wf%n_o,  wf%n_v, wf%n_o)
      call squareup(wf%t2bar, tb_aibj, wf%n_v*wf%n_o)
!
      call mem%alloc(X_ck, wf%n_v, wf%n_o)
!
!$omp parallel do private(k, c)
      do k = 1, wf%n_o
         do c = 1, wf%n_v
!
            X_ck(c, k) = X(c + wf%n_o, k)
!
         enddo
      enddo
!$omp end parallel do
!      
      call dgemm('N','N',           &
                 (wf%n_v)*(wf%n_o), &
                 1,                 &
                 (wf%n_v)*(wf%n_o), &
                 one,               &
                 tb_aibj,           &
                 (wf%n_v)*(wf%n_o), &
                 X_ck,              &
                 (wf%n_v)*(wf%n_o), &
                 one,               &
                 etaX_ai,           &
                 (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(X_ck, wf%n_v, wf%n_o)
!
!     :: Second term: sum_ckdl tb_aick u_ckdl X_ld
!
!     X_ld ordered as X_dl
! 
      call mem%alloc(X_dl, wf%n_v, wf%n_o)
!
!$omp parallel do private(l, d)
      do l = 1, wf%n_o
         do d = 1, wf%n_v
!
            X_dl(d, l) = X(l, d + wf%n_o)
!
         enddo
      enddo
!$omp end parallel do
!
!     Form u_aick
!
      call mem%alloc(t_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_ckdl, wf%n_v*wf%n_o)
!
      call mem%alloc(u_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call zero_array(u_ckdl, (wf%n_o*wf%n_v)**2)
!
      call add_1432_to_1234(-one, t_ckdl, u_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call daxpy((wf%n_o)**2*(wf%n_v)**2, two, t_ckdl, 1, u_ckdl, 1)
!
      call mem%dealloc(t_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Form the intermediate I_ai_dl = sum_ck tb_ai_ck u_ck_dl
!
      call mem%alloc(I_aidl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!      
      call dgemm('N','N',           &
                 (wf%n_v)*(wf%n_o), &
                 (wf%n_v)*(wf%n_o), &
                 (wf%n_v)*(wf%n_o), &
                 one,               &
                 tb_aibj,           &
                 (wf%n_v)*(wf%n_o), &
                 u_ckdl,            &
                 (wf%n_v)*(wf%n_o), &
                 zero,              &
                 I_aidl,            &
                 (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(tb_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(u_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Form sum_dl I_ai_dl X_ld^T
!
      call dgemm('N','N',           &
                 (wf%n_v)*(wf%n_o), &
                 1,                 &
                 (wf%n_v)*(wf%n_o), &
                 one,               &
                 I_aidl,            &
                 (wf%n_v)*(wf%n_o), &
                 X_dl,              &
                 (wf%n_v)*(wf%n_o), &
                 one,               &
                 etaX_ai,           &
                 (wf%n_v)*(wf%n_o))
!
            call mem%dealloc(I_aidl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
            call mem%dealloc(X_dl, wf%n_v, wf%n_o)
!
   end subroutine etaX_eom_doubles_a1_doubles
!
!
   module subroutine etaX_eom_a_doubles(wf, etaX, csiX)
!!
!!    Get eom contribution
!!    Written by Josefine H. Andersen, Feb 2019
!!
!!    Add EOM contribution to etaX vector
!!
!!       EOM correction:  η^X,corr_μ += tbar_μ (ξ * tbar) 
!!
      implicit none
!
      class(doubles), intent(in) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: etaX
      real(dp), dimension(wf%n_es_amplitudes), intent(in)    :: csiX
!
      real(dp) :: X_cc
      real(dp) :: ddot
!
      real(dp), dimension(:), allocatable :: multipliers
!
      call mem%alloc(multipliers, wf%n_es_amplitudes)
!
      call dcopy(wf%n_t1, wf%t1bar, 1, multipliers(1:wf%n_t1), 1)
      call dcopy(wf%n_t2, wf%t2bar, 1, multipliers(wf%n_t1 + 1: wf%n_es_amplitudes), 1)
!
      X_cc = ddot(wf%n_es_amplitudes, multipliers, 1, csiX, 1)
!
      call daxpy(wf%n_es_amplitudes, -X_cc, multipliers, 1, etaX, 1)
!
      call mem%dealloc(multipliers, wf%n_es_amplitudes)
!
   end subroutine etaX_eom_a_doubles
!
!
!
!
end submodule fop_doubles
