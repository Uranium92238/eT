!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2022 the authors of eT
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
submodule (ccsd_class) jacobian_ccsd
!
!!
!!    Jacobian submodule
!!
!!    Routines for the linear transform of trial
!!    vectors by the Jacobian matrix
!!
!!    ρ_i = A * c_i,
!!
!!    where
!!
!!    A_μ,ν = < μ |exp(-T) [H, τ_ν] exp(T) | R >.
!!
!!
!!
!!    (Transfered to the current eT program from the first version
!!    of eT by Andreas Skeidsvoll and Sarai D. Folkestad, 2018.)
!!
!
   implicit none
!
!
contains
!
!
   module subroutine prepare_for_jacobian_ccsd(wf)
!!
!!    Prepare for jacobian
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2019
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
!
      type(timings), allocatable :: prep_timer
!
      prep_timer = timings("Prepare for Jacobian CCSD", pl='normal')
      call prep_timer%turn_on()
!
      call wf%doubles%prepare_for_jacobian()
!
      call wf%save_jacobian_c2_intermediates()
      call wf%save_jacobian_d2_intermediate()
      call wf%save_jacobian_e2_intermediate()
      call wf%save_jacobian_g2_intermediates()
      call wf%save_jacobian_h2_intermediate()
      call wf%save_jacobian_j2_intermediate()
!
      call prep_timer%turn_off()
!
   end subroutine prepare_for_jacobian_ccsd
!
!
   module subroutine jacobian_transformation_ccsd(wf, c, rho)
!!
!!    Jacobian transformation
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!!    Directs the transformation by the CCSD Jacobi matrix,
!!
!!       A_μ,ν = < μ | exp(-T) [H, τ_ν] exp(T) | R >,
!!
!!    where the basis employed for the brackets is biorthonormal.
!!    The transformation is rho = A c, i.e.,
!!
!!       rho_mu = (A c)_mu = sum_ck A_mu,ck c_ck
!!                  + 1/2 sum_ckdl A_mu,ckdl c_ckdl (1 + delta_ck,dl).
!!
!!    On exit, c is overwritten by rho. That is, c_ai = rho_ai,
!!    and c_aibj = rho_aibj.
!!
      use array_initialization, only: zero_array
      use array_utilities, only: scale_diagonal
      use reordering, only: symmetric_sum, sort_1234_to_1324, squareup
      use reordering, only: packin_and_add_from_1324_order
!
      implicit none
!
      class(ccsd), intent(inout) :: wf
!
      real(dp), dimension(wf%n_t1 + wf%n_t2), intent(in)  :: c
      real(dp), dimension(wf%n_t1 + wf%n_t2), intent(out) :: rho
!
      real(dp), dimension(:,:,:,:), allocatable :: c_aibj, c_abij
!
      real(dp), dimension(:,:,:,:), allocatable :: rho_aibj, rho_abij
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian CCSD transformation', pl='normal')
      call timer%turn_on()
!
      call zero_array(rho, wf%n_t1 + wf%n_t2)
!
!     Doubles contributions
!
      call wf%doubles%jacobian_transformation(c, rho)
!
!     CCSD contributions
!
      call mem%alloc(c_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(c(wf%n_t1+1:), c_aibj, wf%n_t1)
      call scale_diagonal(two, c_aibj, wf%n_t1)
!
      call mem%alloc(rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o, set_zero=.true.)
!
!     Contributions from singles vector c
!
      call wf%jacobian_ccsd_b2(rho_aibj, c(1 : wf%n_t1))
      call wf%jacobian_ccsd_c2(rho_aibj, c(1 : wf%n_t1))
      call wf%jacobian_ccsd_d2(rho_aibj, c(1 : wf%n_t1))
!
!     Contributions from doubles vector c
!
      call wf%jacobian_ccsd_e2(rho_aibj, c_aibj)
      call wf%jacobian_ccsd_f2(rho_aibj, c_aibj)
      call wf%jacobian_ccsd_g2(rho_aibj, c_aibj)
      call wf%jacobian_ccsd_h2(rho_aibj, c_aibj)
      call wf%jacobian_ccsd_i2(rho_aibj, c_aibj)
!
!     Last three terms are already symmetric (J2, K2, and L2). Perform the symmetrization
!     rho_aibj = P_ij^ab rho_aibj now, for convenience
!
      call symmetric_sum(rho_aibj, wf%n_t1)
!
!     In preparation for last two terms, reorder
!     rho_aibj to rho_abij, and c_aibj to c_abij
!
      call mem%alloc(rho_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(c_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call sort_1234_to_1324(c_aibj, c_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_1324(rho_aibj, rho_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(c_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%jacobian_ccsd_j2(rho_abij, c_abij)
      call wf%jacobian_ccsd_k2(rho_abij, c_abij)
!
      call mem%dealloc(c_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
!     Scale diagonal by 1/(1 + delta_ai,bj), packin, and add non-A2
!     CCSD doubles contributions
!
      call scale_diagonal(half, rho_abij, wf%n_v, wf%n_o)
!
      call packin_and_add_from_1324_order(rho(wf%n_t1+1:), rho_abij, wf%n_v, wf%n_o)
!
      call mem%dealloc(rho_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
!     Add CCSD-Omega-A2 contribution
!
      call wf%omega_ccsd_a2(rho(wf%n_t1+1:), c(wf%n_t1+1:), right=.true., diagonal_factor=two)
!
      call timer%turn_off()
!
   end subroutine jacobian_transformation_ccsd
!
!
   module subroutine approximate_Jacobian_transform_ccsd(wf, r_or_l, X, R, w)
!!
!!    Approximate Jacobian transform
!!    Written by Eirik F. Kjønstad, Mar 2021
!!
!!    Wrapper for a lower-level Jacobian transformation that is the best approximation
!!    with a lower computational scaling.
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
!
      character(len=*), intent(in) :: r_or_l
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in)  :: X
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: R
!
      real(dp), intent(in), optional :: w
!
      if (trim(r_or_l) == 'left') &
         call output%error_msg('Approximate CCSD Jacobian transpose not &
                                &yet supported.')
!
      call wf%doubles%construct_Jacobian_transform(r_or_l, X, R, w)
!
   end subroutine approximate_Jacobian_transform_ccsd
!
!
   module subroutine prepare_for_approximate_Jacobians_ccsd(wf, r_or_l)
!!
!!    Prepare for approximate Jacobians
!!    Written by Eirik F. Kjønstad, Mar 2021
!!
!!    Wrapper for preparations to a lower-level Jacobian transformation that is
!!    the best approximation with a lower computational scaling.
!!
!!    r_or_l: 'left', 'right', or 'both'
!!            (prepares for A^T, A, or both A^T and A)
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
!
      character(len=*), intent(in) :: r_or_l
!
      if (trim(r_or_l) == 'left') &
         call output%error_msg('Approximate CCSD Jacobian transpose not &
                                &yet supported.')
!
      call wf%doubles%prepare_for_Jacobians(r_or_l)
!
   end subroutine prepare_for_approximate_Jacobians_ccsd
!
!
   module subroutine jacobian_ccsd_b2(wf, rho_aibj, c_ai)
!!
!!    Jacobian CCSD B2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    rho_aibj^B2 = - sum_kc (F_kc t_ij^ac c_bk + F_kc t_ik^ab c_cj)
!!
      use reordering, only: squareup
!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)       :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o)   :: rho_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: t_aicj   ! t_ij^ac
      real(dp), dimension(:,:,:,:), allocatable :: X_kjai   ! An intermediate
!
      real(dp), dimension(:,:), allocatable :: X_kj         ! An intermediate
!
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian CCSD B2 transformation', pl='verbose')
      call timer%turn_on()
!
!     :: Term 1. - sum_kc F_kc t_ij^ac c_bk ::
!
!     Order the amplitudes as t_ca_ij = t_ij^ac
!
      call mem%alloc(t_aicj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call squareup(wf%t2, t_aicj, (wf%n_o)*(wf%n_v))
!
!     Form the intermediate X_k_jai = sum_c F_k_c t_c_jai
!
      call mem%alloc(X_kjai, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',              &
                  wf%n_o,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_v,               &
                  one,                  &
                  wf%fock_ia,           & ! F_k,c
                  wf%n_o,               &
                  t_aicj,               & ! t_cjai = t_c,jai
                  wf%n_v,               &
                  zero,                 &
                  X_kjai,               &
                  wf%n_o)
!
!     Form rho_b_jai = sum_k c_ai(b,k) X_kjai(k,jai)
!
      call dgemm('N', 'N',              &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_o,               &
                  -one,                 &
                  c_ai,                 & ! c_b,k
                  wf%n_v,               &
                  X_kjai,               & ! X_k_jai
                  wf%n_o,               &
                  one,                  &
                  rho_aibj,             & ! rho_b,jai -> will be (bj,ai)-symmetrized
                  wf%n_v)
!
      call mem%dealloc(X_kjai, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     :: Term 2. - sum_kc F_kc t_ik^ab c_cj ::
!
!     Form X_kj = sum_c F_kc c_cj = sum_c fock_ia(k,c) c_ai(c,j)
!
      call mem%alloc(X_kj, wf%n_o, wf%n_o)
!
      call dgemm('N','N',     &
                  wf%n_o,     &
                  wf%n_o,     &
                  wf%n_v,     &
                  one,        &
                  wf%fock_ia, & ! F_k,c
                  wf%n_o,     &
                  c_ai,       & ! c_c,j
                  wf%n_v,     &
                  zero,       &
                  X_kj,       &
                  wf%n_o)
!
!     Form rho_aib_j = - sum_k t_aib_k X_k_j
!
      call dgemm('N','N',               &
                  (wf%n_o)*(wf%n_v)**2, &
                  wf%n_o,               &
                  wf%n_o,               &
                  -one,                 &
                  t_aicj,               & ! t_aib,k
                  (wf%n_o)*(wf%n_v)**2, &
                  X_kj,                 &
                  wf%n_o,               &
                  one,                  &
                  rho_aibj,             & ! rho_aib,j
                  (wf%n_o)*(wf%n_v)**2)
!
      call mem%dealloc(X_kj, wf%n_o, wf%n_o)
      call mem%dealloc(t_aicj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_ccsd_b2
!
!
   module subroutine jacobian_ccsd_c2(wf, rho_aibj, c_ai)
!!
!!    Jacobian CCSD C2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    rho_aibj^C2 = sum_kcl g_ljkc (t_ki^ac c_bl + t_li^bc c_ak + t_lk^ba c_ci)
!!                 - sum_kcl L_ljkc (t_il^ab c_ck + t_ik^ac c_bl)
!!
      use reordering, only: add_1432_to_1234, sort_1234_to_3142
      use reordering, only: sort_1234_to_1342, add_3124_to_1234
      use reordering, only: add_1243_to_1234, add_4213_to_1234
      use reordering, only: sort_1234_to_1423, squareup
!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)           :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: rho_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ljkc ! g_ljkc
      real(dp), dimension(:,:,:,:), allocatable :: L_ljck ! L_ljkc
!
      real(dp), dimension(:,:,:,:), allocatable :: X_ljai  ! An intermediate, term 1
      real(dp), dimension(:,:,:,:), allocatable :: X_kjbi  ! An intermediate, term 2
      real(dp), dimension(:,:,:,:), allocatable :: X_ljki  ! An intermediate, term 3
      real(dp), dimension(:,:,:,:), allocatable :: X_klij  ! X_kj_li reordered
      real(dp), dimension(:,:), allocatable     :: X_lj    ! An intermediate, term 4
      real(dp), dimension(:,:,:,:), allocatable :: Y_ljai  ! An intermediate, term 5
!
      real(dp), dimension(:,:,:,:), allocatable :: t_akci ! t_ki^ac
      real(dp), dimension(:,:,:,:), allocatable :: t_bakl ! t_lk^ba
      real(dp), dimension(:,:,:,:), allocatable :: t_aibl ! t_il^ab
!
      real(dp), dimension(:,:,:,:), allocatable :: rho_ajbi ! rho_aibj, term 2
      real(dp), dimension(:,:,:,:), allocatable :: rho_baij ! rho_aibj, term 3
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian CCSD C2 transformation', pl='verbose')
      call timer%turn_on()
!
!     :: Term 1. sum_kcl g_ljkc t_ki^ac c_bl ::
!
!     Get intermediate X_lj_ai = sum_ck g_lj_kc t_kc_ai
!
      call mem%alloc(X_ljai, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call wf%jacobian_c2_intermediate_oovo_1%open_('read', 'rewind')
      call wf%jacobian_c2_intermediate_oovo_1%read_(X_ljai, (wf%n_o)**3*(wf%n_v))
      call wf%jacobian_c2_intermediate_oovo_1%close_('keep')
!
!     Calculate rho_b_jai = sum_l c_bl X_ljai
!
      call dgemm('N', 'N',              &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_o,               &
                  one,                  &
                  c_ai,                 & ! c_b_l
                  wf%n_v,               &
                  X_ljai,               & ! X_l_jai
                  wf%n_o,               &
                  one,                  &
                  rho_aibj,             & ! rho_b_jai but we will symmetrize later
                  wf%n_v)
!
      call mem%dealloc(X_ljai, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     :: Term 2. sum_kcl g_ljkc t_li^bc c_ak ::
!
!     Get the intermediate X_kjbi = sum_lc g_kjlc t_lcbi
!
      call mem%alloc(X_kjbi, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call wf%jacobian_c2_intermediate_oovo_2%open_('read', 'rewind')
      call wf%jacobian_c2_intermediate_oovo_2%read_(X_kjbi, (wf%n_o)**3*(wf%n_v))
      call wf%jacobian_c2_intermediate_oovo_2%close_('keep')
!
!     Calculate rho_a_jbi = sum_k c_ak X_kjbi
!
      call mem%alloc(rho_ajbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',              &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_o,               &
                  one,                  &
                  c_ai,                 & ! c_a_k
                  wf%n_v,               &
                  X_kjbi,               & ! X_k_jbi
                  wf%n_o,               &
                  zero,                 &
                  rho_ajbi,             & ! rho_a_jbi
                  wf%n_v)
!
      call mem%dealloc(X_kjbi, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     Add rho_ajbi to rho_aibj
!
      call add_1432_to_1234(one, rho_ajbi, rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(rho_ajbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     :: Term 3. sum_kcl g_ljkc t_lk^ba c_ci ::
!
!     Form the intermediate X_ljk_i = sum_c g_ljkc c_ci
!
      call mem%alloc(g_ljkc, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call wf%eri_t1%get('ooov', g_ljkc)
!
      call mem%alloc(X_ljki, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call dgemm('N','N',      &
                  (wf%n_o)**3, &
                  wf%n_o,      &
                  wf%n_v,      &
                  one,         &
                  g_ljkc,      & ! g_ljk_c
                  (wf%n_o)**3, &
                  c_ai,        & ! c_c_i
                  wf%n_v,      &
                  zero,        &
                  X_ljki,      & ! X_ljk_i
                  (wf%n_o)**3)
!
!     Reorder to X_kl_ij = X_kj_li
!
      call mem%alloc(X_klij, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call sort_1234_to_3142(X_ljki, X_klij, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call mem%dealloc(X_ljki, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
!     t_bakl = t_lk^ba = t_blak
!
      call mem%alloc(t_akci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_akci, (wf%n_o)*(wf%n_v))
!
      call mem%alloc(t_bakl, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_1234_to_1342(t_akci, t_bakl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(t_akci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Calculate rho_ba_ij = sum_kcl g_ljkc t_lk^ba c_ci
!                         = sum kl ( sum_c g_ljkc c_ci ) t_lk^ba
!                         = sum_kl t_ba_kl X_kl_ij
!
      call mem%alloc(rho_baij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N','N',      &
                  (wf%n_v)**2, &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  one,         &
                  t_bakl,      & ! t_ba_kl
                  (wf%n_v)**2, &
                  X_klij,      & ! X_kl_ij
                  (wf%n_o)**2, &
                  zero,        &
                  rho_baij,    & ! rho_ba_ij
                  (wf%n_v)**2)
!
      call mem%dealloc(X_klij, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
!     Add rho_baij into rho_aibj
!
      call add_3124_to_1234(one, rho_baij, rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(rho_baij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
!     :: Term 4. - sum_kcl L_ljkc t_il^ab c_ck ::
!
!     Form L_lj_ck(lj,ck) = L_ljkc = 2 * g_ljkc - g_lckj
!                  1234   = 2 * g_ljkc - g_kjlc = 2* g_kj_lc(kj,lc) - g_kj_lc(lj,kc)
!                                                            4213             1243
!
      call mem%alloc(L_ljck, wf%n_o, wf%n_o, wf%n_v, wf%n_o, set_zero=.true.)
!
      call add_1243_to_1234(two, g_ljkc, L_ljck, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call add_4213_to_1234(-one, g_ljkc, L_ljck, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_ljkc, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
!     Calculate the intermediate X_lj = sum_ck L_ljck c_ck
!
      call mem%alloc(X_lj, wf%n_o, wf%n_o)
!
      call dgemv('N',                &
                  (wf%n_o)**2,       &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  L_ljck,            & ! L_lj_ck
                  (wf%n_o)**2,       &
                  c_ai,              & ! c_ck
                  1,                 &
                  zero,              &
                  X_lj,              & ! X_lj
                  1)
!
      call mem%dealloc(L_ljck, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     Order the amplitudes as t_ai_bl = t_il^ab
!
      call mem%alloc(t_aibl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     t_lk^ba = t_ba_kl(ba,kl) => t_il^ab = t_bakl(ab,li) = t_aibl(ai,bl)
!                                                  1234            1423
!
      call sort_1234_to_1423(t_bakl, t_aibl, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call mem%dealloc(t_bakl, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
!     Form rho_aibj =+ - sum_l t_il^ab X_lj = - sum_l t_aib_l X_lj
!
      call dgemm('N','N',               &
                  (wf%n_o)*(wf%n_v)**2, &
                  wf%n_o,               &
                  wf%n_o,               &
                  -one,                 &
                  t_aibl,               & ! t_aib_l
                  (wf%n_o)*(wf%n_v)**2, &
                  X_lj,                 & ! X_l_j
                  wf%n_o,               &
                  one,                  &
                  rho_aibj,             & ! rho_aib_j
                  (wf%n_o)*(wf%n_v)**2)
!
      call mem%dealloc(t_aibl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X_lj, wf%n_o, wf%n_o)
!
!     :: Term 5. - sum_kcl L_ljkc t_ik^ac c_bl ::
!
!     t_il^ab = t_ai_bl(ai,bl) => t_ai_bl(ck,ai) = t_ki^ca = t_ik^ac
!
!     Get the intermediate Y_lj_ai = sum_kc L_ljkc t_ik^ac = sum_kc L_lj_ck t_ck_ai
!
      call mem%alloc(Y_ljai, wf%n_o, wf%n_o, wf%n_v,wf%n_o)
      call wf%jacobian_c2_intermediate_oovo_3%open_('read', 'rewind')
      call wf%jacobian_c2_intermediate_oovo_3%read_(Y_ljai, (wf%n_o)**3*(wf%n_v))
      call wf%jacobian_c2_intermediate_oovo_3%close_('keep')
!
!     Calculate rho_b_jai =+ - sum_l c_bl Y_lj_ai
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_o,               &
                  -one,                 &
                  c_ai,                 & ! c_b_l
                  wf%n_v,               &
                  Y_ljai,               & ! Y_l_jai
                  wf%n_o,               &
                  one,                  &
                  rho_aibj,             & ! rho_b_jai we will symmetrize later
                  wf%n_v)
!
      call mem%dealloc(Y_ljai, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_ccsd_c2
!
!
   module subroutine jacobian_ccsd_d2(wf, rho_aibj, c_ai)
!!
!!    Jacobian CCSD D2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018 & 2021
!!
!!    rho_aibj^D2 = - sum_kcd g_kcbd (t_ij^cd c_ak + t_kj^ad c_ci + t_ik^ca c_dj)
!!                       + sum_kcd L_kcbd (t_ik^ac c_dj + t_ij^ad c_ck)
!!
!!    The first term is constructed via a pre-calculated intermediate X:
!!
!!       - sum_kcd g_kcbd t_ij^cd c_ak = - c_ak X_kijb
!!
!!    The remaining terms are:
!!
!!       -g_kcbd c_ci t_kj^ad = - (L_J,kc c_ci) L_J,bd t_akdj = - W_J,ki L_J,bd t_akdj
!!       -g_kcbd c_dj t_ik^ca = - L_J,kc (L_J,bd c_dj) t_akci = - L_J,kc V_J,bj t_akci
!!      2 g_kcbd c_dj t_ik^ac = 2 L_J,kc (L_J,bd c_dj) t_aick = 2 L_J,kc V_J,bj t_aick
!!       -g_kdbc c_dj t_ik^ac = - (L_J,kd c_dj) L_J,bc t_aick = - W_J,kj L_J,bc t_aick
!!      2 g_kcbd c_ck t_ij^ad = 2 (L_J,kc c_ck) L_J,bd t_aidj = 2 Y_J    L_J,bd t_aidj
!!       -g_kdbc c_ck t_ij^ad = - L_J,kd (L_J,bc c_ck) t_aidj = - L_J,kd V_J,bk t_aidj
!!
!!    We perform the final contractions as:
!!
!!     1:     -  W_J,ki L_J,bd t_akdj   = -WL_kibd t_akdj (o3v3)
!!     2:     - (L_J,kc t_akci) V_J,bj  = - M_J,ai V_J,bj (N5)
!!     3:      2 L_J,kc V_J,bj t_aick   = 2 (L_J,kc t_aick) V_J,bj = 2 Z_J,ai V_J,bj (N5)
!!     4:      - W_J,kj L_J,bd t_aidk   = - WL_kjbd t_aidk (o3v3)
!!     5 and 6:  2 (Y_J L_J,bd) t_aidj - (L_J,kd V_J,bk) t_aidj
!!                                      = (2 N_bd - Z_bd) t_aidj (N5)
!!                                      = X_bd t_aidj
!!
      use batching_index_class, only: batching_index
      use array_initialization, only: copy_and_scale
      use reordering, only: sort_12_to_21, sort_123_to_213, sort_1234_to_4132
      use reordering, only: sort_1234_to_1432, sort_1234_to_3241
      use reordering, only: sort_1234_to_1243, sort_1234_to_1423
      use reordering, only: add_1243_to_1234, add_1432_to_1234, squareup
!
      implicit none
!
      class(ccsd), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout) :: rho_aibj
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
!
      type(timings), allocatable :: timer
!
      real(dp), dimension(:,:,:,:), allocatable ::  X_kijb, rho_aijb
!
      real(dp), dimension(:), allocatable :: Y_J
!
      real(dp), dimension(:,:), allocatable :: c_kc, N_vv, X_vv, Z_vv
!
      real(dp), dimension(:,:,:), allocatable :: M_J_vo, V_J_vo, V_vJo
!
      real(dp), dimension(:,:,:), pointer :: L_J_ov, L_J_vv
!
      real(dp), dimension(:,:,:), allocatable :: W_J_oo, Z_J_vo
!
      real(dp), dimension(:,:,:,:), allocatable :: rho_ajbi, t_aijd, t_aikc, t_ajdk
      real(dp), dimension(:,:,:,:), allocatable :: t_vovo, WL_oovv, WL_bidk, WL_dkbj
!
      integer :: req0, req1, batch
      type(batching_index), allocatable :: batch_v
!
      timer = timings('Jacobian CCSD D2 transformation', pl='verbose')
      call timer%turn_on()
!
!     Term 1: - sum_kcd g_kcbd t_ij^cd c_ak = - c_ak X_kijb
!
      call mem%alloc(X_kijb, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call wf%jacobian_d2_intermediate%open_('read','rewind')
      call wf%jacobian_d2_intermediate%read_(X_kijb, (wf%n_o)**3*(wf%n_v))
      call wf%jacobian_d2_intermediate%close_('keep')
!
      call mem%alloc(rho_aijb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N', 'N',                &
                  wf%n_v,                 &
                  (wf%n_v)*(wf%n_o)**2,   &
                  wf%n_o,                 &
                  -one,                   &
                  c_ai,                   & ! c_a,k
                  wf%n_v,                 &
                  X_kijb,                 & ! X_k,ijb
                  wf%n_o,                 &
                  zero,                   &
                  rho_aijb,               & ! rho_a,ijb
                  wf%n_v)
!
      call add_1243_to_1234(one, rho_aijb, rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(rho_aijb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call mem%dealloc(X_kijb, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
!     Construct intermediates used in terms 2 to 6
!
!     W_J,ki = L_J,kc c_ci
!
      call wf%L_t1%load_block(L_J_ov, 1, wf%n_o, wf%n_o + 1, wf%n_mo)
!
      call mem%alloc(W_J_oo, wf%L_t1%n_J, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N',             &
                 wf%L_t1%n_J*wf%n_o,   &
                 wf%n_o,               &
                 wf%n_v,               &
                 one,                  &
                 L_J_ov,               &
                 wf%L_t1%n_J*wf%n_o,   &
                 c_ai,                 &
                 wf%n_v,               &
                 zero,                 &
                 W_J_oo,               &
                 wf%L_t1%n_J*wf%n_o)
!
!     Y_J = L_J,kc c_ck
!
      call mem%alloc(Y_J, wf%L_t1%n_J)
!
      call mem%alloc(c_kc, wf%n_o, wf%n_v)
      call sort_12_to_21(c_ai, c_kc, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',          &
                  wf%L_t1%n_J,      &
                  1,                &
                  wf%n_o * wf%n_v,  &
                  one,              &
                  L_J_ov,           &
                  wf%L_t1%n_J,      &
                  c_kc,             &
                  wf%n_o * wf%n_v,  &
                  zero,             &
                  Y_J,              &
                  wf%L_t1%n_J)
!
      call mem%dealloc(c_kc, wf%n_o, wf%n_v)
!
!     Make intermediates that require L_J_vv
!
      call mem%alloc(V_J_vo, wf%L_t1%n_J, wf%n_v, wf%n_o, set_zero=.true.)
!
      call mem%alloc(WL_oovv, wf%n_o, wf%n_o, wf%n_v, wf%n_v, set_zero=.true.)
!
      call mem%alloc(N_vv, wf%n_v, wf%n_v, set_zero=.true.)
!
      req0 = 0
      req1 = wf%L_t1%load_memory_estimate(wf%n_o + 1, wf%n_mo, wf%n_o + 1, wf%n_o + 1)
!
      batch_v = batching_index(wf%n_v)
!
      call mem%batch_setup(batch_v, req0, req1, 'jacobian_ccsd_d2')
!
      do batch = 1, batch_v%num_batches
!
         call batch_v%determine_limits(batch)
!
         call wf%L_t1%load_block(L_J_vv,                      &
                                 wf%n_o + 1,                  &
                                 wf%n_mo,                     &
                                 wf%n_o + batch_v%first,      &
                                 wf%n_o + batch_v%get_last())
!
!        V_J,bj = L_J,bd c_dj
!
         call dgemm('N', 'N',                      &
                     wf%L_t1%n_J * wf%n_v,         &
                     wf%n_o,                       &
                     batch_v%length,               &
                     one,                          &
                     L_J_vv(:,:,1:batch_v%length), &! L_Jb,d
                     wf%L_t1%n_J * wf%n_v,         &
                     c_ai(batch_v%first, 1),       &! c_d,j
                     wf%n_v,                       &
                     one,                          &
                     V_J_vo,                       &
                     wf%L_t1%n_J * wf%n_v)
!
!        WL_kibd = W_J,ki L_J,bd
!
         call dgemm('T', 'N',                                           &
                     wf%n_o**2,                                         &
                     wf%n_v * batch_v%length,                           &
                     wf%L_t1%n_J,                                       &
                     one,                                               &
                     W_J_oo,                                            &
                     wf%L_t1%n_J,                                       &
                     L_J_vv(:,:,1:batch_v%length),                      &
                     wf%L_t1%n_J,                                       &
                     one,                                               &
                     WL_oovv(:,:,:,batch_v%first:batch_v%get_last()),   &
                     wf%n_o**2)
!
!        N_bd = Y_J L_J,bd
!
         call dgemm('T','N',                          &
                     wf%n_v * batch_v%length,         &
                     1,                               &
                     wf%L_t1%n_J,                     &
                     one,                             &
                     L_J_vv(:,:,1:batch_v%length),    &
                     wf%L_t1%n_J,                     &
                     Y_J,                             &
                     wf%L_t1%n_J,                     &
                     one,                             &
                     N_vv(1, batch_v%first),          &
                     wf%n_v * batch_v%length)
!
         call wf%L_t1%offload_block(wf%n_o + 1,                  &
                                    wf%n_mo,                     &
                                    wf%n_o + batch_v%first,      &
                                    wf%n_o + batch_v%get_last())
!
      enddo
!
      call mem%batch_finalize()
!
!     Z_bd = L_J,kd V_J,bk
!
!     Sort V_Jbk as V_bJk
!
      call mem%alloc(V_vJo, wf%n_v, wf%L_t1%n_J, wf%n_o)
!
      call sort_123_to_213(V_J_vo, V_vJo, wf%L_t1%n_J, wf%n_v, wf%n_o)
!
      call mem%alloc(Z_vv, wf%n_v, wf%n_v)
!
      call dgemm('N', 'N',             &
                  wf%n_v,              &
                  wf%n_v,              &
                  wf%L_t1%n_J * wf%n_o,&
                  one,                 &
                  V_vJo,               &! V_b,Jk
                  wf%n_v,              &
                  L_J_ov,              &! L_Jk,d
                  wf%L_t1%n_J * wf%n_o,&
                  zero,                &
                  Z_vv,                &
                  wf%n_v)
!
      call mem%dealloc(V_vJo, wf%n_v, wf%L_t1%n_J, wf%n_o)
!
!     X_bd = two * N_bd - Z_bd
!
      call mem%alloc(X_vv, wf%n_v, wf%n_v)
!
      call copy_and_scale(two, N_vv, X_vv, wf%n_v**2)
      call daxpy(wf%n_v**2, -one, Z_vv, 1, X_vv, 1)
!
      call mem%dealloc(N_vv, wf%n_v, wf%n_v)
      call mem%dealloc(Z_vv, wf%n_v, wf%n_v)
!
!     Z_J,ai = L_J,kc t_aick
!
      call mem%alloc(t_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call squareup(wf%t2, t_vovo, wf%n_t1)
!
!     Sort t_aick to t_aikc
!
      call mem%alloc(t_aikc, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call sort_1234_to_1243(t_vovo, t_aikc, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(Z_J_vo, wf%L_t1%n_J, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T',          &
                  wf%L_t1%n_J,      &
                  wf%n_v * wf%n_o,  &
                  wf%n_v * wf%n_o,  &
                  one,              &
                  L_J_ov,           &! L_J,kc
                  wf%L_t1%n_J,      &
                  t_aikc,           &! t_ai,kc
                  wf%n_v * wf%n_o,  &
                  zero,             &
                  Z_J_vo,           &
                  wf%L_t1%n_J)
!
!     M_J,ai = L_J,kc t_akci
!
!     Sort t_akci to t_aikc
!
      call sort_1234_to_1423(t_vovo, t_aikc, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(M_J_vo, wf%L_t1%n_J, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T',          &
                  wf%L_t1%n_J,      &
                  wf%n_o * wf%n_v,  &
                  wf%n_o * wf%n_v,  &
                  one,              &
                  L_J_ov,           &! L_J,kc
                  wf%L_t1%n_J,      &
                  t_aikc,           &! t_ai,kc
                  wf%n_o * wf%n_v,  &
                  zero,             &
                  M_J_vo,           &
                  wf%L_t1%n_J)
!

      call wf%L_t1%offload_block(1, wf%n_o, wf%n_o + 1, wf%n_mo)
!
      call mem%dealloc(t_aikc, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
!     Term 1: -WL_kibd t_akdj
!
      call mem%alloc(t_ajdk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call sort_1234_to_1432(t_vovo, t_ajdk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(WL_bidk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call sort_1234_to_3241(WL_oovv, WL_bidk, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
      call mem%alloc(rho_ajbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T',          &
                  wf%n_v * wf%n_o,  &
                  wf%n_v * wf%n_o,  &
                  wf%n_v * wf%n_o,  &
                  -one,             &
                  t_ajdk,           &
                  wf%n_v * wf%n_o,  &
                  WL_bidk,          &
                  wf%n_v * wf%n_o,  &
                  zero,             &
                  rho_ajbi,         &
                  wf%n_v * wf%n_o)
!
      call mem%dealloc(WL_bidk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call add_1432_to_1234(one, rho_ajbi, rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(t_ajdk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_ajbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Term 4. - WL_kjbd t_aidk
!
      call mem%alloc(WL_dkbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call sort_1234_to_4132(WL_oovv, WL_dkbj, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
      call dgemm('N', 'N',          &
                  wf%n_v * wf%n_o,  &
                  wf%n_v * wf%n_o,  &
                  wf%n_v * wf%n_o,  &
                  -one,             &
                  t_vovo,           & ! t_ai,dk
                  wf%n_v * wf%n_o,  &
                  WL_dkbj,          &
                  wf%n_v * wf%n_o,  &
                  one,              &
                  rho_aibj,         &
                  wf%n_v * wf%n_o)
!
      call mem%dealloc(WL_dkbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Term 2. - M_J,ai V_J,bj
!
      call dgemm('T', 'N',          &
                  wf%n_v * wf%n_o,  &
                  wf%n_v * wf%n_o,  &
                  wf%L_t1%n_J,      &
                  -one,             &
                  M_J_vo,           &
                  wf%L_t1%n_J,      &
                  V_J_vo,           &
                  wf%L_t1%n_J,      &
                  one,              &
                  rho_aibj,         &
                  wf%n_v * wf%n_o)
!
      call mem%dealloc(M_J_vo, wf%L_t1%n_J, wf%n_v, wf%n_o)
!
!     Term 3. 2 Z_J,ai V_J,bj
!
      call dgemm('T', 'N',          &
                  wf%n_v * wf%n_o,  &
                  wf%n_v * wf%n_o,  &
                  wf%L_t1%n_J,      &
                  two,              &
                  Z_J_vo,           &
                  wf%L_t1%n_J,      &
                  V_J_vo,           &
                  wf%L_t1%n_J,      &
                  one,              &
                  rho_aibj,         &
                  wf%n_v * wf%n_o)
!
!     Terms 5 and 6. X_bd t_aidj
!
      call mem%alloc(t_aijd, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call sort_1234_to_1243(t_vovo, t_aijd, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(rho_aijb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N', 'T',             &
                  wf%n_o**2 * wf%n_v,  &
                  wf%n_v,              &
                  wf%n_v,              &
                  one,                 &
                  t_aijd,              &
                  wf%n_o**2 * wf%n_v,  &
                  X_vv,                &
                  wf%n_v,              &
                  zero,                &
                  rho_aijb,            &
                  wf%n_o**2 * wf%n_v)
!
      call add_1243_to_1234(one, rho_aijb, rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(rho_aijb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call mem%dealloc(t_aijd, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call mem%dealloc(t_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(WL_oovv, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
      call mem%dealloc(V_J_vo, wf%L_t1%n_J, wf%n_v, wf%n_o)
      call mem%dealloc(W_J_oo, wf%L_t1%n_J, wf%n_o, wf%n_o)
!
      call mem%dealloc(Y_J, wf%L_t1%n_J)
      call mem%dealloc(X_vv, wf%n_v, wf%n_v)
      call mem%dealloc(Z_J_vo, wf%L_t1%n_J, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_ccsd_d2
!
!
   module subroutine jacobian_ccsd_e2(wf, rho_aibj, c_aick)
!!
!!    Jacobian CCSD E2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    rho_aibj^E2 = 2 sum_dlck t_bj,dl * L_kc,ld * c_ai,ck
!!                  - sum_ckld t_bj,dl * L_kc,ld * c_ak,ci
!!                = sum_ck Y_bjck (2c_aick - c_akci)
!!                = sum_ck Y_bjck v_aick
!!
!!    with
!!
!!       Y_bjck = t_bj,dl * L_kc,ld
!!       v_aick = (2c_aick - c_akci)
!!
!!    which is constructed in prepare_for_jacobian
!!
      use array_initialization, only: copy_and_scale
      use reordering, only: add_1432_to_1234
!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: rho_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_aick
!
      real(dp), dimension(:,:,:,:), allocatable :: Y_bjck
      real(dp), dimension(:,:,:,:), allocatable :: v_aick
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian CCSD E2 transformation', pl='verbose')
      call timer%turn_on()
!
      call mem%alloc(Y_bjck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%jacobian_e2_intermediate%open_('read', 'rewind')
      call wf%jacobian_e2_intermediate%read_(Y_bjck, (wf%n_v**2)*(wf%n_o**2))
      call wf%jacobian_e2_intermediate%close_('keep')
!
      call mem%alloc(v_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call copy_and_scale(two, c_aick, v_aick, (wf%n_v**2)*(wf%n_o**2))
      call add_1432_to_1234(-one, c_aick, v_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T',             &
                  (wf%n_v)*(wf%n_o),   &
                  (wf%n_v)*(wf%n_o),   &
                  (wf%n_v)*(wf%n_o),   &
                  one,                 &
                  v_aick,              &
                  (wf%n_v)*(wf%n_o),   &
                  Y_bjck,              &
                  (wf%n_v)*(wf%n_o),   &
                  one,                 &
                  rho_aibj,            &
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(v_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(Y_bjck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_ccsd_e2
!
!
   module subroutine jacobian_ccsd_f2(wf, rho_aibj, c_aibj)
!!
!!    Jacobian CCSD F2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!       rho_aibj^F2 =   - sum_ckdl t_ai,dj * L_kc,ld * c_bl,ck
!!                       - sum_ckdl t_ai_bl * L_kc,ld * c_ck,dj
!!
!!    L_kc,ld = 2*g_kc,ld - g_kd,lc = 2*g_kcld(kc,ld) - 2*g_kcld(kd,lc)
!!
      use reordering, only: add_4123_to_1234, add_4321_to_1234, squareup
!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: rho_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: g_kcld
      real(dp), dimension(:,:,:,:), allocatable :: L_dlck
!
      real(dp), dimension(:,:,:,:), allocatable :: t_djai
 !
      real(dp), dimension(:,:), allocatable :: Y_d_b
      real(dp), dimension(:,:), allocatable :: Z_jl
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian CCSD F2 transformation', pl='verbose')
      call timer%turn_on()
!
!     :: Term 1: - sum_ckdl t_aidj * L_kcld * c_blck
!
      call mem%alloc(g_kcld, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%eri_t1%get('ovov', g_kcld)
!
!     Construct L_ckdl reordered as L_dlck
!
      call mem%alloc(L_dlck, wf%n_v, wf%n_o, wf%n_v, wf%n_o, set_zero=.true.)
!
!     L_dlck(d,l,c,k) =- g_kcld(k,d,l,c) (4123)
!
      call add_4123_to_1234(-one, g_kcld, L_dlck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     L_dlck(d,l,c,k) =+ 2*g_kcld(k,c,l,d) (4321)
!
      call add_4321_to_1234(two, g_kcld, L_dlck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_kcld, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     Y_d_b = sum_clk L_d_lck * c_b_lck
!
      call mem%alloc(Y_d_b, wf%n_v, wf%n_v)
!
      call dgemm('N', 'T',              &
                  wf%n_v,               &
                  wf%n_v,               &
                  ((wf%n_o)**2)*wf%n_v, &
                  one,                  &
                  L_dlck,               & ! L_d_lck
                  wf%n_v,               &
                  c_aibj,               & ! c_b_lck
                  wf%n_v,               &
                  zero,                 &
                  Y_d_b,                & ! Y_d_b
                  wf%n_v)
!
!     rho_aij_b =- sum_d Y_d_b t_djai
!
      call mem%alloc(t_djai, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call squareup(wf%t2, t_djai, wf%n_t1)
!
      call dgemm('T', 'N',                &
                  wf%n_v,                 &
                  (wf%n_v)*(wf%n_o)**2,   &
                  wf%n_v,                 &
                  -one,                   &
                  Y_d_b,                  & ! Y_d,b
                  wf%n_v,                 &
                  t_djai,                 & ! t_d,jai
                  wf%n_v,                 &
                  one,                    &
                  rho_aibj,               & ! rho_b,jai -> will be (bj,ai)-symmtrized
                  wf%n_v)
!
      call mem%dealloc(Y_d_b, wf%n_v, wf%n_v)
!
!     :: Term 3: - sum_ckdl t_aibl * L_kcld * c_ckdj ::
!
!     Note: Using symmetry L_dlck = L_ckdl (L_ldkc = L_kcld)
!
!     Z_jl = sum_ckd L_kdlc c_ckd_j
!
      call mem%alloc(Z_jl, wf%n_o, wf%n_o)
!
      call dgemm('T', 'N',              &
                  wf%n_o,               &
                  wf%n_o,               &
                  ((wf%n_v)**2)*wf%n_o, &
                  one,                  &
                  c_aibj,               & ! c_ckd_j
                  ((wf%n_v)**2)*wf%n_o, &
                  L_dlck,               &
                  ((wf%n_v)**2)*wf%n_o, &
                  zero,                 &
                  Z_jl,                 &
                  wf%n_o)
!
      call mem%dealloc(L_dlck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     rho_aibj = sum_l t_aibl Z_l_j
!
      call dgemm('N','T',                 &
                  ((wf%n_v)**2)*(wf%n_o), &
                  wf%n_o,                 &
                  wf%n_o,                 &
                  -one,                   &
                  t_djai,                 & ! t_aib,l
                  ((wf%n_v)**2)*(wf%n_o), &
                  Z_jl,                   &
                  wf%n_o,                 &
                  one,                    &
                  rho_aibj,               &
                  ((wf%n_v)**2)*(wf%n_o))
!
      call mem%dealloc(Z_jl, wf%n_o, wf%n_o)
      call mem%dealloc(t_djai, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_ccsd_f2
!
!
   module subroutine jacobian_ccsd_g2(wf, rho_aibj, c_aibj)
!!
!!    Jacobian CCSD G2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    rho_aibj^G2 =  - sum_ckdl t_bl,dj * L_kc,ld * c_ai,ck
!!                       - sum_ckdl t_ck_bl * L_kc,ld * c_ai,dj
!!                       - sum_ckld t_ck,dj * L_kc,ld * c_ai,bl
!!                = - sum_ck Y_bjck c_aick
!!                  - sum_d Y_bd c_aidj
!!                  - sum_l Y_jl c_aibl
!!
!!    The intermediates are constructed once in prepare_for_jacobian
!!    in the routine save_jacobian_g2_intermediates
!!
!!    L_kc,ld = 2*g_kc,ld - g_kd,lc = 2*g_kcld(kc,ld) - 2*g_kcld(kd,lc)
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o)               :: rho_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)   :: c_aibj
!
      real(dp), dimension(:,:), allocatable :: Y_bd
      real(dp), dimension(:,:), allocatable :: Y_jl
      real(dp), dimension(:,:,:,:), allocatable :: Y_bjck
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian CCSD G2 transformation', pl='verbose')
      call timer%turn_on()
!
!     Term 1: - Y_bjck c_aick
!
      call mem%alloc(Y_bjck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%jacobian_g2_intermediate_vovo%open_('read', 'rewind')
      call wf%jacobian_g2_intermediate_vovo%read_(Y_bjck, wf%n_t1**2)
      call wf%jacobian_g2_intermediate_vovo%close_('keep')
!
      call dgemm('N', 'T',             &
                  (wf%n_v)*(wf%n_o),   &
                  (wf%n_v)*(wf%n_o),   &
                  (wf%n_v)*(wf%n_o),   &
                  -one,                &
                  c_aibj,              & ! c_aick
                  (wf%n_v)*(wf%n_o),   &
                  Y_bjck,              &
                  (wf%n_v)*(wf%n_o),   &
                  one,                 &
                  rho_aibj,            &
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(Y_bjck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Term 2: - sum_d Y_bd c_aidj
!
      call mem%alloc(Y_bd, wf%n_v, wf%n_v)
!
      call wf%jacobian_a1_intermediate_vv%open_('read', 'rewind')
      call wf%jacobian_a1_intermediate_vv%read_(Y_bd, wf%n_v**2)
      call wf%jacobian_a1_intermediate_vv%close_('keep')
!
      call dgemm('N', 'N',             &
                  wf%n_v,              &
                  wf%n_v*(wf%n_o**2),  &
                  wf%n_v,              &
                  -one,                &
                  Y_bd,                &
                  wf%n_v,              &
                  c_aibj,              & ! c_djai
                  wf%n_v,              &
                  one,                 &
                  rho_aibj,            & ! rho_bjai but we will symmetrize anyway
                  wf%n_v)
!
      call mem%dealloc(Y_bd, wf%n_v, wf%n_v)
!
!     Term 3: - Y_jl c_aibl
!
      call mem%alloc(Y_jl, wf%n_o, wf%n_o)
!
      call wf%jacobian_a1_intermediate_oo%open_('read', 'rewind')
      call wf%jacobian_a1_intermediate_oo%read_(Y_jl, wf%n_o**2)
      call wf%jacobian_a1_intermediate_oo%close_('keep')
!
      call dgemm('N', 'T',             &
                  (wf%n_v**2)*wf%n_o,  &
                  wf%n_o,              &
                  wf%n_o,              &
                  -one,                &
                  c_aibj,              & ! c_aibl
                  (wf%n_v**2)*wf%n_o,  &
                  Y_jl,                &
                  wf%n_o,              &
                  one,                 &
                  rho_aibj,            &
                  (wf%n_v**2)*wf%n_o)
!
      call mem%dealloc(Y_jl, wf%n_o, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_ccsd_g2
!
!
   module subroutine jacobian_ccsd_h2(wf, rho_aibj, c_aibj)
!!
!!    Jacobian CCSD H2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!       rho_aibj^H2 =  sum_ckdl t_ci,ak * g_kc,ld * c_bl,dj
!!                     + sum_ckdl t_cj,al * g_kc,ld * c_bk,di
!!                   = sum_dl Y_aild c_bldj
!!                     + sum_dk Y_ajkd c_bkdi
!!
      use reordering, only: sort_1234_to_2314, sort_1234_to_3241
      use reordering, only: squareup, add_1432_to_1234
!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o)               :: rho_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)   :: c_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: c_ldbj
!
      real(dp), dimension(:,:,:,:), allocatable :: Y_ajkd
!
      real(dp), dimension(:,:,:,:), allocatable :: rho_ajbi
!
      real(dp), dimension(:,:,:,:), allocatable :: t_vovo, t_aikc
!
      real(dp), dimension(:,:,:), allocatable :: L_J_ov, V_J_vo, W_J_vo
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian CCSD H2 transformation', pl='verbose')
      call timer%turn_on()
!
!     Term 1: Y_aild c_bldj
!
!     t_ciak g_kcld c_bldj
!     = (L_kc^J t_ciak) (L_ld^J c_bldj)
!     = W_ai^J V_bj^J
!
      call mem%alloc(L_J_ov, wf%eri_t1%n_J, wf%n_o, wf%n_v)
      call wf%L_t1%get(L_J_ov, &
                       1, wf%n_o, &
                       wf%n_o + 1, wf%n_mo)
!
!     t_ciak -> t_aikc
!
      call mem%alloc(t_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_vovo, wf%n_t1)
!
      call mem%alloc(t_aikc, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_3241(t_vovo, t_aikc, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(t_vovo, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(W_J_vo, wf%eri_t1%n_J, wf%n_v, wf%n_o)
!
      call dgemm('N', 'T',          &
                  wf%eri_t1%n_J,    &
                  wf%n_v * wf%n_o,  &
                  wf%n_v * wf%n_o,  &
                  one,              &
                  L_J_ov,           & ! L_J,kc
                  wf%eri_t1%n_J,    &
                  t_aikc,           &
                  wf%n_v * wf%n_o,  &
                  zero,             &
                  W_J_vo,           &
                  wf%eri_t1%n_J)
!
      call mem%dealloc(t_aikc, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     c_bldj -> c_ldbj
!
      call mem%alloc(c_ldbj, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call sort_1234_to_2314(c_aibj, c_ldbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(V_J_vo, wf%eri_t1%n_J, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',          &
                  wf%eri_t1%n_J,    &
                  wf%n_v * wf%n_o,  &
                  wf%n_v * wf%n_o,  &
                  one,              &
                  L_J_ov,           &  ! L_J,ld
                  wf%eri_t1%n_J,    &
                  c_ldbj,           &
                  wf%n_v * wf%n_o,  &
                  zero,             &
                  V_J_vo,           &
                  wf%eri_t1%n_J)
!
      call mem%dealloc(L_J_ov, wf%eri_t1%n_J, wf%n_o, wf%n_v)
!
      call dgemm('T', 'N',          &
                  wf%n_v * wf%n_o,  &
                  wf%n_v * wf%n_o,  &
                  wf%eri_t1%n_J,    &
                  one,              &
                  W_J_vo,           &
                  wf%eri_t1%n_J,    &
                  V_J_vo,           &
                  wf%eri_t1%n_J,    &
                  one,              &
                  rho_aibj,         &
                  wf%n_v * wf%n_o)
!
      call mem%dealloc(V_J_vo, wf%eri_t1%n_J, wf%n_v, wf%n_o)
      call mem%dealloc(W_J_vo, wf%eri_t1%n_J, wf%n_v, wf%n_o)
!
!     Term 2: Y_ajkd c_bkdi
!
!     Pretend c_ldbj is c_kdbi
!
      call mem%alloc(Y_ajkd, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call wf%jacobian_h2_intermediate%open_('read', 'rewind')
      call wf%jacobian_h2_intermediate%read_(Y_ajkd, wf%n_t1**2)
      call wf%jacobian_h2_intermediate%close_('keep')
!
      call mem%alloc(rho_ajbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',       &
                  wf%n_v*wf%n_o, &
                  wf%n_v*wf%n_o, &
                  wf%n_v*wf%n_o, &
                  one,           &
                  Y_ajkd,        &
                  wf%n_v*wf%n_o, &
                  c_ldbj,        & ! c_kdbi
                  wf%n_v*wf%n_o, &
                  zero,          &
                  rho_ajbi,      &
                  wf%n_v*wf%n_o)
!
      call mem%dealloc(Y_ajkd, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call mem%dealloc(c_ldbj, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      call add_1432_to_1234(one, rho_ajbi, rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(rho_ajbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_ccsd_h2
!
!
   module subroutine jacobian_ccsd_i2(wf, rho_aibj, c_aibj)
!!
!!    Jacobian CCSD I2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    rho_aibj^I2 = sum_ck L_bj,kc * c_ai,ck
!!                - sum_ck ( g_kc,bj * c_ak,ci + g_ki,bc * c_ak,cj )
!!
      use reordering, only: sort_1234_to_2314, sort_1234_to_1432
      use reordering, only: add_3421_to_1234, add_3124_to_1234, add_3214_to_1234
!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o)               :: rho_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)   :: c_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: rho_biaj
!
      real(dp), dimension(:,:,:,:), allocatable :: L_ckbj
!
      real(dp), dimension(:,:,:,:), allocatable :: c_kcai
      real(dp), dimension(:,:,:,:), allocatable :: g_bikc
      real(dp), dimension(:,:,:,:), allocatable :: g_bjkc
      real(dp), dimension(:,:,:,:), allocatable :: g_bckj
!
      real(dp), dimension(:,:,:), allocatable :: L_J_ov, L_J_vo, W_J_vo
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian CCSD I2 transformation', pl='verbose')
      call timer%turn_on()
!
!     rho_bjai =- g_bjkc c_akci = g_bj,kc c_kc,ai
!
!     L_bj^J L_kc^J c_kc,ai = L_bj^J W_ai^J
!
      call mem%alloc(c_kcai, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call sort_1234_to_2314(c_aibj, c_kcai, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(L_J_ov, wf%eri_t1%n_J, wf%n_o, wf%n_v)
!
      call wf%L_t1%get(L_J_ov,       &
                       1, wf%n_o,    &
                       wf%n_o + 1, wf%n_mo)
!
      call mem%alloc(W_J_vo, wf%eri_t1%n_J, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',             &
                  wf%eri_t1%n_J,       &
                  wf%n_o * wf%n_v,     &
                  wf%n_o * wf%n_v,     &
                  one,                 &
                  L_J_ov,              &
                  wf%eri_t1%n_J,       &
                  c_kcai,              &
                  wf%n_o * wf%n_v,     &
                  zero,                &
                  W_J_vo,              &
                  wf%eri_t1%n_J)
!

      call mem%dealloc(L_J_ov, wf%eri_t1%n_J, wf%n_o, wf%n_v)
!
      call mem%alloc(L_J_vo, wf%eri_t1%n_J, wf%n_v, wf%n_o)
!
      call wf%L_t1%get(L_J_vo,                &
                       wf%n_o + 1, wf%n_mo,   &
                       1, wf%n_o)
!
      call dgemm('T', 'N',          &
                  wf%n_o * wf%n_v,  &
                  wf%n_o * wf%n_v,  &
                  wf%eri_t1%n_J,    &
                  -one,             &
                  W_J_vo,           &
                  wf%eri_t1%n_J,    &
                  L_J_vo,           &
                  wf%eri_t1%n_J,    &
                  one,              &
                  rho_aibj,         &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(W_J_vo, wf%eri_t1%n_J, wf%n_v, wf%n_o)
      call mem%dealloc(L_J_vo, wf%eri_t1%n_J, wf%n_v, wf%n_o)
!
!     rho_aibj =+ c_aick L_bjkc = c_aick L_ck,bj
!
      call mem%alloc(g_bjkc, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%eri_t1%get('voov', g_bjkc)
!
!     Construct L_bjkc ordered as L_ckbj = 2 g_bjkc - g_bckj
!
      call mem%alloc(L_ckbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o, set_zero=.true.)
!
      call add_3421_to_1234(two, g_bjkc, L_ckbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_bjkc, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(g_bckj, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call wf%eri_t1%get('vvoo', g_bckj)
!
      call add_3124_to_1234(-one, g_bckj, L_ckbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     rho_aibj =- g_bcki c_akcj =- g_bi,kc c_kc,aj
!
      call mem%alloc(g_bikc, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_1234_to_1432(g_bckj, g_bikc, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call mem%dealloc(g_bckj, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call mem%alloc(rho_biaj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',             &
                  (wf%n_o)*(wf%n_v),   &
                  (wf%n_o)*(wf%n_v),   &
                  (wf%n_o)*(wf%n_v),   &
                  -one,                &
                  g_bikc,              &
                  (wf%n_o)*(wf%n_v),   &
                  c_kcai,              & ! c_kc,aj
                  (wf%n_o)*(wf%n_v),   &
                  zero,                &
                  rho_biaj,            &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(g_bikc, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call mem%dealloc(c_kcai, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      call add_3214_to_1234(one, rho_biaj, rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(rho_biaj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     rho_aibj =+ c_aick L_ck,bj
!
      call dgemm('N', 'N',             &
                  (wf%n_o)*(wf%n_v),   &
                  (wf%n_o)*(wf%n_v),   &
                  (wf%n_o)*(wf%n_v),   &
                  one,                 &
                  c_aibj,              & ! c_ai,ck
                  (wf%n_o)*(wf%n_v),   &
                  L_ckbj,              &
                  (wf%n_o)*(wf%n_v),   &
                  one,                 &
                  rho_aibj,            &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(L_ckbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_ccsd_i2
!
!
   module subroutine jacobian_ccsd_j2(wf, rho_abij, c_abij)
!!
!!    Jacobian CCSD J2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!       rho_abij^J2 =    sum_ckld t_ci,dj * g_kc,ld * c_ak,bl
!!                       + sum_ckdl t_ak,bl * g_kc,ld * c_ci,dj
!!                   =    sum_ckld Y_klij * c_ak,bl
!!                       + sum_ckdl t_ak,bl * g_kc,ld * c_ci,dj
      use reordering, only: sort_1234_to_1324, squareup_and_sort_1234_to_1324
!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o) :: rho_abij
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: c_abij
!
      real(dp), dimension(:,:,:,:), allocatable :: g_kcld
      real(dp), dimension(:,:,:,:), allocatable :: g_klcd
!
      real(dp), dimension(:,:,:,:), allocatable :: t_abij
!
      real(dp), dimension(:,:,:,:), allocatable :: Y_klij
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian CCSD J2 transformation', pl='verbose')
      call timer%turn_on()
!
!     Y_kl_ij = g_kl_cd * t_cd_ij
!
      call mem%alloc(Y_klij, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call wf%jacobian_j2_intermediate_oooo%open_('read', 'rewind')
!
      call wf%jacobian_j2_intermediate_oooo%read_(Y_klij, wf%n_o**4)
!
      call wf%jacobian_j2_intermediate_oooo%close_('keep')
!
!     rho_abij += c_ab_kl * X_kl_ij
!
      call dgemm('N', 'N',     &
                  (wf%n_v)**2, &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  one,         &
                  c_abij,      & ! c_ab_kl
                  (wf%n_v)**2, &
                  Y_klij,      & ! X_kl_ij
                  (wf%n_o)**2, &
                  one,         &
                  rho_abij,    & ! rho_ab_ij
                  (wf%n_v)**2)
!
!     Y_kl_ij = g_kl_cd * c_cd_ij
!
!     Construct g_kcld ordered as g_klcd
!
      call mem%alloc(g_kcld, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%eri_t1%get('ovov', g_kcld)
!
      call mem%alloc(g_klcd, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
      call sort_1234_to_1324(g_kcld, g_klcd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%dealloc(g_kcld, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call dgemm('N', 'N',     &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  (wf%n_v)**2, &
                  one,         &
                  g_klcd,      & ! g_kl_cd
                  (wf%n_o)**2, &
                  c_abij,      & ! c_cd_ij
                  (wf%n_v)**2, &
                  zero,        &
                  Y_klij,      & ! X_kl_ij
                  (wf%n_o)**2)
!
      call mem%dealloc(g_klcd, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
!     rho_abij += t_abkl * X_klij
!
      call mem%alloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call squareup_and_sort_1234_to_1324(wf%t2, t_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',     &
                  (wf%n_v)**2, &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  one,         &
                  t_abij,      & ! t_ab_kl
                  (wf%n_v)**2, &
                  Y_klij,      & ! X_kl_ij
                  (wf%n_o)**2, &
                  one,         &
                  rho_abij,    & ! rho_ab_ij
                  (wf%n_v)**2)
!
      call mem%dealloc(Y_klij, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
      call mem%dealloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_ccsd_j2
!
!
   module subroutine jacobian_ccsd_k2(wf, rho_abij, c_abij)
!!
!!    Jacobian CCSD K2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    rho_abij^K2 =    sum_kl g_kilj * c_akbl
!!
      use reordering, only: sort_1234_to_1324
!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o) :: rho_abij
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: c_abij
!
      real(dp), dimension(:,:,:,:), allocatable :: g_kilj
      real(dp), dimension(:,:,:,:), allocatable :: g_klij
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian CCSD K2 transformation', pl='verbose')
      call timer%turn_on()
!
!     Construct g_kilj ordered as g_klij
!
      call mem%alloc(g_kilj, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call wf%eri_t1%get('oooo', g_kilj)
!
      call mem%alloc(g_klij, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call sort_1234_to_1324(g_kilj, g_klij, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call mem%dealloc(g_kilj, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
!     rho_abij += sum_kl g_kilj * c_akbl = sum_kl c_abij(a,b,k,l) g_klij
!
      call dgemm('N', 'N',     &
                  (wf%n_v)**2, &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  one,         &
                  c_abij,      & ! c_ab_kl
                  (wf%n_v)**2, &
                  g_klij,      & ! g_kl_ij
                  (wf%n_o)**2, &
                  one,         &
                  rho_abij,    & ! rho_ab_ij
                  (wf%n_v)**2)
!
      call mem%dealloc(g_klij, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_ccsd_k2
!
!
   module subroutine save_jacobian_c2_intermediates(wf)
!!
!!    Save jacobian c2 intermediates
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2019
!!
!!    Constructs the intermediates for Jacobian C2:
!!
!!       X_ljai = sum_ck g_ljkc t_ki^ac
!!       X_kjbi = sum_lc g_ljkc t_li^bc
!!       Y_ljai = sum_kc L_ljkc t_ik^ac
!!
!!    used in the c2-term. This is done only once in prepare_for_jacobian
!!    and the intermediate is stored in the file jacobian_c2_intermediate
!!    which is a wf variable.
!!
      use reordering, only: sort_1234_to_2314, sort_1234_to_3214, squareup
      use reordering, only: add_3214_to_1234, sort_1234_to_4312
!
      implicit none
!
      class(ccsd) :: wf
!
      type(timings), allocatable :: timer
!
      real(dp), dimension(:,:,:,:), allocatable :: X_ljai
      real(dp), dimension(:,:,:,:), allocatable :: X_kjbi
      real(dp), dimension(:,:,:,:), allocatable :: Y_ljai
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ljkc, t_akci, t_kcai, g_kjlc, L_ljkc
!
      timer = timings('Jacobian CCSD C2 intermediates construction', pl='verbose')
      call timer%turn_on()
!
!     Intermediate X_ljai
!
!     Form g_ljkc
!
      call mem%alloc(g_ljkc, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call wf%eri_t1%get('ooov', g_ljkc)
!
!     Square up (t_ak_ci = t_ki^ac)
!
      call mem%alloc(t_akci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_akci, (wf%n_o)*(wf%n_v))
!
!     Order as t_kc_ai = t_ki^ac
!
      call mem%alloc(t_kcai, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
!     Reorder t_ak_ci to t_kc_ai
!
      call sort_1234_to_2314(t_akci, t_kcai, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Form X_lj_ai = sum_ck g_lj_kc t_kc_ai
!
      call mem%alloc(X_ljai, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',            &
                  (wf%n_o)**2,       &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  g_ljkc,            & ! g_lj_kc
                  (wf%n_o)**2,       &
                  t_kcai,            & ! t_kc_ai
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_ljai,            & ! X_lj_ai
                  (wf%n_o)**2)
!
      wf%jacobian_c2_intermediate_oovo_1 = sequential_file('jacobian_c2_intermediate_oovo_1_ccsd')
      call wf%jacobian_c2_intermediate_oovo_1%open_('write', 'rewind')
      call wf%jacobian_c2_intermediate_oovo_1%write_(X_ljai, (wf%n_o)**3*(wf%n_v))
      call wf%jacobian_c2_intermediate_oovo_1%close_('keep')
!
      call mem%dealloc(X_ljai, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     Intermediate X_kjbi
!
!     Reorder to g_kj_lc = g_lj_kc = g_ljkc
!                  3214      1234
!
      call mem%alloc(g_kjlc, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call sort_1234_to_3214(g_ljkc, g_kjlc, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
!     Form the intermediate X_kjbi = sum_lc g_kjlc t_lcbi
!
      call mem%alloc(X_kjbi, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',           &
                  (wf%n_o)**2,       &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  g_kjlc,            & ! g_kj_lc
                  (wf%n_o)**2,       &
                  t_kcai,            & ! t_lc_bi
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_kjbi,            & ! X_kj_bi
                  (wf%n_o)**2)
!
      call mem%dealloc(t_kcai, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call mem%dealloc(g_kjlc, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      wf%jacobian_c2_intermediate_oovo_2 = sequential_file('jacobian_c2_intermediate_oovo_2_ccsd')
      call wf%jacobian_c2_intermediate_oovo_2%open_('write', 'rewind')
      call wf%jacobian_c2_intermediate_oovo_2%write_(X_kjbi, (wf%n_o)**3*(wf%n_v))
      call wf%jacobian_c2_intermediate_oovo_2%close_('keep')
!
      call mem%dealloc(X_kjbi, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     Intermediate Y_ljai = sum_kc L_ljkc t_ik^ac = sum_kc L_lj,kc t_kc,ai
!
!     L_lj,kc = 2 g_ljkc - g_kjlc
!
      call mem%alloc(L_ljkc, wf%N_o, wf%n_o, wf%n_o, wf%n_v)
!
      L_ljkc = two*g_ljkc
      call add_3214_to_1234(-one, g_ljkc, L_ljkc, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
!     t_kc,ai = t_ik^ac = t_aick
!
      call mem%alloc(t_kcai, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call sort_1234_to_4312(t_akci, t_kcai, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(Y_ljai, wf%n_o, wf%n_o, wf%n_v,wf%n_o)
!
      call dgemm('N','N',            &
                  (wf%n_o)**2,       &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  L_ljkc,            &
                  (wf%n_o)**2,       &
                  t_kcai,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  Y_ljai,            &
                  (wf%n_o)**2)
!
      call mem%dealloc(L_ljkc, wf%N_o, wf%n_o, wf%n_o, wf%n_v)
      call mem%dealloc(t_kcai, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      wf%jacobian_c2_intermediate_oovo_3 = sequential_file('jacobian_c2_intermediate_oovo_3_ccsd')
      call wf%jacobian_c2_intermediate_oovo_3%open_('write', 'rewind')
      call wf%jacobian_c2_intermediate_oovo_3%write_(Y_ljai, (wf%n_o)**3*(wf%n_v))
      call wf%jacobian_c2_intermediate_oovo_3%close_('keep')
!
      call mem%dealloc(Y_ljai, wf%n_o, wf%n_o, wf%n_v,wf%n_o)
      call mem%dealloc(g_ljkc, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call mem%dealloc(t_akci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine save_jacobian_c2_intermediates
!
!
   module subroutine save_jacobian_d2_intermediate(wf)
!!
!!    Save jacobian d2 intermediate
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2019
!!
!!    Constructs the intermediate
!!
!!       X_kbij = sum_dl g_kcbd t_ij^cd
!!
!!    used in the d2-term. This is done only once in prepare_for_jacobian
!!    and the intermediate is stored in the file jacobian_d2_intermediate
!!    which is a wf variable.
!!
      use batching_index_class, only: batching_index
      use reordering, only: squareup_and_sort_1234_to_2413
      use reordering, only: sort_1234_to_4231, sort_1234_to_3124
!
      implicit none
!
      class(ccsd) :: wf
!
      type(timings), allocatable :: timer
!
      integer :: req1, req0
      integer :: current_b_batch
!
      type(batching_index) :: batch_b
!
      real(dp), dimension(:,:,:,:), allocatable :: g_bdkc, g_cdkb, X_ij_kb, X_k_ijb, X_kijb_full, t_ijcd
!
      timer = timings('Jacobian CCSD D2 intermediate construction', pl='verbose')
      call timer%turn_on()
!
      call mem%alloc(X_kijb_full, wf%n_o, wf%n_o, wf%n_o, wf%n_v, set_zero=.true.)
!
!     Order amplitudes as t_ij_cd = t_ij^cd = t_ci_dj
!
      call mem%alloc(t_ijcd, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call squareup_and_sort_1234_to_2413(wf%t2, t_ijcd, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Initialize batching variable
!
      req0 = wf%n_v*wf%n_o*wf%eri_t1%n_J
      req1 = wf%n_v*wf%eri_t1%n_J + 2*(wf%n_o)*(wf%n_v)**2 + 2*(wf%n_o)**3
!
      batch_b = batching_index(wf%n_v)
      call mem%batch_setup(batch_b, req0, req1, 'save_jacobian_d2_intermediate')
!
!     Start looping over b-batches
!
      do current_b_batch = 1, batch_b%num_batches
!
!        Get batching limits for current b-batch
!
         call batch_b%determine_limits(current_b_batch)
!
!        Form g_kc_db = g_kcbd
!
         call mem%alloc(g_bdkc, batch_b%length, wf%n_v, wf%n_o, wf%n_v)
!
         call wf%eri_t1%get('vvov', g_bdkc, first_p=batch_b%first, last_p=batch_b%get_last())
!
!        Reorder g_bd_kc to g_cd_kb (= g_kcbd), i.e. 1234 to 4231
!
         call mem%alloc(g_cdkb, wf%n_v, wf%n_v, wf%n_o, batch_b%length)
!
         call sort_1234_to_4231(g_bdkc, g_cdkb, batch_b%length, wf%n_v, wf%n_o, wf%n_v)
!
         call mem%dealloc(g_bdkc, batch_b%length, wf%n_v, wf%n_o, wf%n_v)

!
!        Form intermediate X_ij_kb = sum_cd g_kcdb t_ij^cd
!                                  = sum_cd t_ij_cd g_cd_kb
!
         call mem%alloc(X_ij_kb, wf%n_o, wf%n_o, wf%n_o, batch_b%length)
!
         call dgemm('N', 'N',                   &
                     (wf%n_o)**2,               &
                     (wf%n_o)*(batch_b%length), &
                     (wf%n_v)**2,               &
                     one,                       &
                     t_ijcd,                    & ! t_ij_cd
                     (wf%n_o)**2,               &
                     g_cdkb,                    & ! g_cd_kb
                     (wf%n_v)**2,               &
                     zero,                      &
                     X_ij_kb,                   & ! X_ij_kb
                     (wf%n_o)**2)
!
         call mem%dealloc(g_cdkb, wf%n_v, wf%n_v, wf%n_o, batch_b%length)
!
!        Reorder to X_k_ijb = X_ij_kb & add to the full space X
!
         call mem%alloc(X_k_ijb, wf%n_o, wf%n_o, wf%n_o, batch_b%length)
         call sort_1234_to_3124(X_ij_kb, X_k_ijb, wf%n_o, wf%n_o, wf%n_o, (batch_b%length))
         call mem%dealloc(X_ij_kb, wf%n_o, wf%n_o, wf%n_o, batch_b%length)
!
         call daxpy((wf%n_o)**3*(batch_b%length), one, X_k_ijb, 1, X_kijb_full(1,1,1,batch_b%first), 1)
         call mem%dealloc(X_k_ijb, wf%n_o, wf%n_o, wf%n_o, batch_b%length)
!
      enddo
!
      call mem%batch_finalize()
!
!     Store intermediate to file
!
      wf%jacobian_d2_intermediate = sequential_file('jacobian_d2_intermediate_ccsd')
      call wf%jacobian_d2_intermediate%open_('write', 'rewind')
!
      call wf%jacobian_d2_intermediate%write_(X_kijb_full, (wf%n_o)**3*(wf%n_v))
!
      call mem%dealloc(X_kijb_full, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call mem%dealloc(t_ijcd, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
      call wf%jacobian_d2_intermediate%close_('keep')
!
      call timer%turn_off()
!
   end subroutine save_jacobian_d2_intermediate
!
!
   module subroutine save_jacobian_e2_intermediate(wf)
!!
!!    Save jacobian e2 intermediate
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2019
!!
!!    Constructs the intermediate
!!
!!       Y_bjck = sum_dl t_bjdl L_ldkc
!!
!!    used in the e2-term. This is done only once in prepare_for_jacobian
!!    and the intermediate is stored in the file jacobian_e2_intermediate
!!    which is a wf variable
!!
      use reordering, only: add_2143_to_1234, add_2341_to_1234, squareup
!
      implicit none
!
      class(ccsd) :: wf
!
      type(timings), allocatable :: timer
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ldkc
      real(dp), dimension(:,:,:,:), allocatable :: L_dlck
      real(dp), dimension(:,:,:,:), allocatable :: t_bjdl
      real(dp), dimension(:,:,:,:), allocatable :: Y_bjck
!
      timer = timings('Jacobian CCSD E2 intermediate construction', pl='verbose')
      call timer%turn_on()
!
      call mem%alloc(g_ldkc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call wf%eri_t1%get('ovov', g_ldkc)
!
      call mem%alloc(L_dlck, wf%n_v, wf%n_o, wf%n_v, wf%n_o, set_zero=.true.)
      call add_2143_to_1234(two, g_ldkc, L_dlck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call add_2341_to_1234(-one, g_ldkc, L_dlck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_ldkc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%alloc(t_bjdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call squareup(wf%t2, t_bjdl, wf%n_t1)
!
      call mem%alloc(Y_bjck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',            &
                (wf%n_v)*(wf%n_o),   &
                (wf%n_v)*(wf%n_o),   &
                (wf%n_v)*(wf%n_o),   &
                one,                 &
                t_bjdl,              &
                (wf%n_v)*(wf%n_o),   &
                L_dlck,              &
                (wf%n_v)*(wf%n_o),   &
                zero,                &
                Y_bjck,              &
                (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(L_dlck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(t_bjdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      wf%jacobian_e2_intermediate = sequential_file('jacobian_e2_intermediate_ccsd')
      call wf%jacobian_e2_intermediate%open_('write', 'rewind')
!
      call wf%jacobian_e2_intermediate%write_(Y_bjck, wf%n_t1**2)
!
      call mem%dealloc(Y_bjck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%jacobian_e2_intermediate%close_('keep')
!
      call timer%turn_off()
!
   end subroutine save_jacobian_e2_intermediate
!
!
   module subroutine save_jacobian_g2_intermediates(wf)
!!
!!    Save jacobian g2 intermediates
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2019
!!
!!    Constructs the intermediate
!!
!!       Y_bjck = t_bldj L_kcdl
!!
!!    The intermediates are stored in the file:
!!
!!       jacobian_g2_intermediate_vovo
!!
!!    G2 also needs intermediates formed for A1 term:
!!
!!       Y_jl   = t_ckdj L_kcld
!!       Y_bd   = t_blck L_kcld
!!
!!    Which are constructed in save_jacobian_a1_intermediates
!!    and stored on files
!!
!!       jacobian_a1_intermediate_oo
!!       jacobian_a1_intermediate_vv
!!
!!    which are wavefunction variables
!!
      use reordering, only: add_2143_to_1234, add_2341_to_1234
      use reordering, only: squareup, sort_1234_to_1432
!
      implicit none
!
      class(ccsd) :: wf
!
      type(timings), allocatable :: timer
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ldkc
      real(dp), dimension(:,:,:,:), allocatable :: L_dlck
      real(dp), dimension(:,:,:,:), allocatable :: Y_bjck
      real(dp), dimension(:,:,:,:), allocatable :: t_blck
      real(dp), dimension(:,:,:,:), allocatable :: t_bjdl
!
      timer = timings('Jacobian CCSD G2 intermediate construction', pl='verbose')
      call timer%turn_on()
!
      call mem%alloc(g_ldkc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call wf%eri_t1%get('ovov', g_ldkc)
!
      call mem%alloc(L_dlck, wf%n_v, wf%n_o, wf%n_v, wf%n_o, set_zero=.true.)
!
      call add_2143_to_1234(two, g_ldkc, L_dlck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call add_2341_to_1234(-one, g_ldkc, L_dlck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_ldkc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%alloc(t_blck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_blck, wf%n_t1)
!
!     Y_bjck = t_bldj L_kcdl
!
!     Reorder t_bldj to t_bjdl
!
      call mem%alloc(t_bjdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_1432(t_blck, t_bjdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(t_blck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(Y_bjck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',       &
                  wf%n_v*wf%n_o, &
                  wf%n_v*wf%n_o, &
                  wf%n_v*wf%n_o, &
                  one,           &
                  t_bjdl,        &
                  wf%n_v*wf%n_o, &
                  L_dlck,        &
                  wf%n_v*wf%n_o, &
                  zero,          &
                  Y_bjck,        &
                  wf%n_v*wf%n_o)
!
      call mem%dealloc(t_bjdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(L_dlck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      wf%jacobian_g2_intermediate_vovo = sequential_file('jacobian_g2_intermediate_vovo_ccsd')
      call wf%jacobian_g2_intermediate_vovo%open_('write', 'rewind')
!
      call wf%jacobian_g2_intermediate_vovo%write_(Y_bjck, wf%n_t1**2)
!
      call mem%dealloc(Y_bjck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%jacobian_g2_intermediate_vovo%close_('keep')
!
      call timer%turn_off()
!
   end subroutine save_jacobian_g2_intermediates
!
!
   module subroutine save_jacobian_h2_intermediate(wf)
!!
!!    Save jacobian h2 intermediates
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2019
!!
!!    Constructs the intermediates
!!
!!       Y_ajkd = t_cjal g_kcld
!!
!!    The intermediate is stored in the file
!!
!!       jacobian_h2_intermediate_voov
!!
      use reordering, only: squareup_and_sort_1234_to_1423
      use reordering, only: sort_1234_to_1432
!
      implicit none
!
      class(ccsd) :: wf
!
      type(timings), allocatable :: timer
!
      real(dp), dimension(:,:,:,:), allocatable :: g_kcld
      real(dp), dimension(:,:,:,:), allocatable :: g_lckd
      real(dp), dimension(:,:,:,:), allocatable :: t_ajcl
      real(dp), dimension(:,:,:,:), allocatable :: Y_ajkd

      timer = timings('Jacobian CCSD H2 intermediate construction', pl='verbose')
      call timer%turn_on()
!
!     Y_ajkd = t_alcj g_kcld
!
      call mem%alloc(g_kcld, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call wf%eri_t1%get('ovov', g_kcld)
!
      call mem%alloc(t_ajcl, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call squareup_and_sort_1234_to_1423(wf%t2, t_ajcl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(g_lckd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call sort_1234_to_1432(g_kcld, g_lckd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%dealloc(g_kcld, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%alloc(Y_ajkd, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N', 'N',       &
                  wf%n_v*wf%n_o, &
                  wf%n_v*wf%n_o, &
                  wf%n_v*wf%n_o, &
                  one,           &
                  t_ajcl,        &
                  wf%n_v*wf%n_o, &
                  g_lckd,        &
                  wf%n_v*wf%n_o, &
                  zero,          &
                  Y_ajkd,        &
                  wf%n_v*wf%n_o)
!
      call mem%dealloc(g_lckd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%dealloc(t_ajcl, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      wf%jacobian_h2_intermediate = sequential_file('jacobian_h2_intermediate_ccsd')
      call wf%jacobian_h2_intermediate%open_('write', 'rewind')
!
      call wf%jacobian_h2_intermediate%write_(Y_ajkd, wf%n_t1**2)
!
      call mem%dealloc(Y_ajkd, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call wf%jacobian_h2_intermediate%close_('keep')
!
      call timer%turn_off()
!
   end subroutine save_jacobian_h2_intermediate
!
!
   module subroutine save_jacobian_j2_intermediate(wf)
!!
!!    Save jacobian j2 intermediate
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2019
!!
!!    Constructs the intermediate
!!
!!       Y_klij = t_cidj g_kcld
!!
!!    The intermediates are stored in the files:
!!
!!       jacobian_j2_intermediate_oooo
!!
!!    which are wavefunction variables
!!
      use reordering, only: sort_1234_to_1324
      use reordering, only: squareup_and_sort_1234_to_1324
!
      implicit none
!
      class(ccsd) :: wf
!
      type(timings), allocatable :: timer
!
      real(dp), dimension(:,:,:,:), allocatable :: g_kcld
      real(dp), dimension(:,:,:,:), allocatable :: g_klcd
      real(dp), dimension(:,:,:,:), allocatable :: t_cdij
      real(dp), dimension(:,:,:,:), allocatable :: Y_klij

      timer = timings('Jacobian CCSD J2 intermediate construction', pl='verbose')
      call timer%turn_on()
!
      call mem%alloc(g_kcld, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call wf%eri_t1%get('ovov', g_kcld)
!
      call mem%alloc(g_klcd, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
      call sort_1234_to_1324(g_kcld, g_klcd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%dealloc(g_kcld, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%alloc(t_cdij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call squareup_and_sort_1234_to_1324(wf%t2, t_cdij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(Y_klij, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N',    &
                  wf%n_o**2,  &
                  wf%n_o**2,  &
                  wf%n_v**2,  &
                  one,        &
                  g_klcd,     &
                  wf%n_o**2,  &
                  t_cdij,     &
                  wf%n_v**2,  &
                  zero,       &
                  Y_klij,     &
                  wf%n_o**2)
!
      call mem%dealloc(t_cdij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(g_klcd, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
      wf%jacobian_j2_intermediate_oooo = sequential_file('jacobian_j2_intermediate_oooo_ccsd')
      call wf%jacobian_j2_intermediate_oooo%open_('write', 'rewind')
!
      call wf%jacobian_j2_intermediate_oooo%write_(Y_klij, wf%n_o**4)
!
      call mem%dealloc(Y_klij, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call wf%jacobian_j2_intermediate_oooo%close_('keep')
!
      call timer%turn_off()
!
   end subroutine save_jacobian_j2_intermediate
!
!
end submodule jacobian_ccsd
