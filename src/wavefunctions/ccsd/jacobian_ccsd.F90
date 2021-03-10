!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2021 the authors of eT
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
!!    A_μ,ν = < μ | exp(-T) [H, τ_ν] exp(T) | R >.
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
      prep_timer = timings("Time preparing for CCSD Jacobian", pl='normal')
      call prep_timer%turn_on()
!
      call wf%save_jacobian_a1_intermediates()
      call wf%save_jacobian_c2_intermediates()
      call wf%save_jacobian_d2_intermediate()
      call wf%save_jacobian_e2_intermediate()
      call wf%save_jacobian_g2_intermediates()
      call wf%save_jacobian_h2_intermediates()
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
!     Allocate and zero the transformed vector
!
      call zero_array(rho, wf%n_t1 + wf%n_t2)
!
!     :: CCS contributions to the singles c vector ::
!
      call wf%ccs%jacobian_transformation(c(1 : wf%n_t1), rho(1 : wf%n_t1))
!
!     :: CCSD contributions to the transformed singles vector ::
!
      call wf%jacobian_doubles_a1(rho(1 : wf%n_t1), c(1 : wf%n_t1))
!
      call mem%alloc(c_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call squareup(c(wf%n_t1+1:), c_aibj, wf%n_t1)
!
!     Scale the doubles vector by 1 + delta_ai,bj, i.e.
!     redefine to c_ckdl = c_ckdl (1 + delta_ck,dl)
!
      call scale_diagonal(two, c_aibj, wf%n_t1)
      call wf%jacobian_doubles_b1(rho(1 : wf%n_t1), c_aibj)
      call wf%jacobian_doubles_c1(rho(1 : wf%n_t1), c_aibj)
      call wf%jacobian_doubles_d1(rho(1 : wf%n_t1), c_aibj)
!
!     :: CCSD contributions to the transformed doubles vector ::
!
!     Allocate unpacked transformed vector
!
      call mem%alloc(rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call zero_array(rho_aibj, wf%n_t1**2)
!
!     Contributions from singles vector c
!
      call wf%jacobian_doubles_a2(rho_aibj, c(1 : wf%n_t1))
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
      call symmetric_sum(rho_aibj, (wf%n_v)*(wf%n_o))
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
!     Scale diagonal by 1/(1 + delta_ai,bj) and packin
!
      call scale_diagonal(half, rho_abij, wf%n_v, wf%n_o)
!
!     Pack in the rho vector
!
      call packin(rho(wf%n_t1+1:), rho_abij, wf%n_v, wf%n_o)
!
      call mem%dealloc(rho_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call wf%omega_ccsd_a2(rho(wf%n_t1+1:), c(wf%n_t1+1:), right=.true., diagonal_factor=two)
!
      call timer%turn_off()
!
   end subroutine jacobian_transformation_ccsd
!
!
   module subroutine jacobian_ccsd_b2_ccsd(wf, rho_aibj, c_ai)
!!
!!    Jacobian CCSD B2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    rho_aibj^B2 = - sum_kc (F_kc t_ij^ac c_bk + F_kc t_ik^ab c_cj)
!!
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
   end subroutine jacobian_ccsd_b2_ccsd
!
!
   module subroutine jacobian_ccsd_c2_ccsd(wf, rho_aibj, c_ai)
!!
!!    Jacobian CCSD C2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    rho_aibj^C2 = sum_kcl g_ljkc (t_ki^ac c_bl + t_li^bc c_ak + t_lk^ba c_ci)
!!                 - sum_kcl L_ljkc (t_il^ab c_ck + t_ik^ac c_bl)
!!
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
      real(dp), dimension(:,:,:,:), allocatable :: X_ljai   ! An intermediate, term 1
      real(dp), dimension(:,:,:,:), allocatable :: X_kjbi   ! An intermediate, term 2
      real(dp), dimension(:,:,:,:), allocatable :: X_ljki   ! An intermediate, term 3
      real(dp), dimension(:,:,:,:), allocatable :: X_klij   ! X_kj_li reordered
      real(dp), dimension(:,:), allocatable     :: X_lj     ! An intermediate, term 4
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
      call wf%eri%get_eri_t1('ooov', g_ljkc)
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
      call mem%alloc(L_ljck, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call zero_array(L_ljck, (wf%n_o**3)*wf%n_v)
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
   end subroutine jacobian_ccsd_c2_ccsd
!
!
  module subroutine jacobian_ccsd_d2_ccsd(wf, rho_aibj, c_ai)
!!
!!    Jacobian CCSD D2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    rho_aibj^D2 = - sum_kcd g_kcbd (t_ij^cd c_ak + t_kj^ad c_ci + t_ik^ca c_dj)
!!                       + sum_kcd L_kcbd (t_ik^ac c_dj + t_ij^ad c_ck)
!!
!!    Note: the code is structured so that we batch over the index b,
!!          where the integrals are made as g_kc_db = g_kcbd and held
!!          in some ordering or other throughout a given batch (i.e.,
!!          all five terms are constructed gradually in the batches).
!!
!!    The first term is constructed from the D2 intermediate: 
!!
!!       - sum_kcd g_kcbd t_ij^cd c_ak = - sum_d c_ak X_kijb  
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)     :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: rho_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: g_bdkc ! g_kcbd
      real(dp), dimension(:,:,:,:), allocatable :: g_ckbd ! g_kcbd reordered
      real(dp), dimension(:,:,:,:), allocatable :: L_ckbd ! L_kcbd = 2 g_kcbd - g-kdbc
!
      real(dp), dimension(:,:,:,:), allocatable :: t_dkaj ! t_kj^ad
      real(dp), dimension(:,:,:,:), allocatable :: t_aick ! t_ik^ca
      real(dp), dimension(:,:,:,:), allocatable :: t_aijd ! t_ij^ad
!
      real(dp), dimension(:,:,:,:), allocatable :: X_k_ijb ! Intermediate, term 1
      real(dp), dimension(:,:,:,:), allocatable :: X_bdki  ! Intermediate, term 2
      real(dp), dimension(:,:,:,:), allocatable :: X_ib_dk ! The above intermediate, reordered
      real(dp), dimension(:,:,:,:), allocatable :: X_ckb_j ! An intermediate, term 3
      real(dp), dimension(:,:,:,:), allocatable :: Y_ckb_j ! An intermediate, term 4
!      
      real(dp), dimension(:,:), allocatable :: X_bd    ! An intermediate, term 5
!
      real(dp), dimension(:,:,:,:), allocatable :: rho_aijb ! rho_aibj, batching over b
      real(dp), dimension(:,:,:,:), allocatable :: rho_ibaj ! rho_aibj, batching over b
      real(dp), dimension(:,:,:,:), allocatable :: rho_aib_j ! rho_aibj, batching over b
      real(dp), dimension(:,:,:,:), allocatable :: rho_baij ! rho_aibj, batching over b
!
      integer :: req1, req0
      integer :: current_b_batch
!
      type(batching_index) :: batch_b

      integer :: b, i, j, a
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian CCSD D2 transformation', pl='verbose')
      call timer%turn_on()
!
!     Use intermediate to form term 1
!
      call mem%alloc(X_k_ijb, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call wf%jacobian_d2_intermediate%open_('read','rewind')
      call wf%jacobian_d2_intermediate%read_(X_k_ijb, (wf%n_o)**3*(wf%n_v))
      call wf%jacobian_d2_intermediate%close_('keep')
!
      call mem%alloc(rho_aijb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N', 'N',                &
                  wf%n_v,                 &
                  (wf%n_v)*(wf%n_o)**2,   &
                  wf%n_o,                 &
                  -one,                   &
                  c_ai,                   & ! c_a_k
                  wf%n_v,                 &
                  X_k_ijb,                & ! X_k_ijb
                  wf%n_o,                 &
                  zero,                   &
                  rho_aijb,               & ! rho_aijb
                  wf%n_v)   
!
      call add_1243_to_1234(one, rho_aijb, rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o) 
      call mem%dealloc(rho_aijb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call mem%dealloc(X_k_ijb, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
!     Initialize batching variable
!
      req0 = (wf%n_o**2)*(wf%n_v**2) + wf%n_v*wf%n_o*wf%eri%n_J
      req1 = wf%n_v**2*wf%n_o                         &
            + wf%n_v*wf%eri%n_J                       &
            + max(2*(wf%n_o**3),                      &
                  (wf%n_o**3) + (wf%n_o**2)*(wf%n_v), &
                  2*(wf%n_o)*(wf%n_v**2),             &
                  2*(wf%n_o**2)*(wf%n_v))

!
      batch_b = batching_index(wf%n_v)
      call mem%batch_setup(batch_b, req0, req1)
!
!     Start looping over b-batches
!
      do current_b_batch = 1, batch_b%num_batches
!
!        Get batching limits for current b-batch
!
         call batch_b%determine_limits(current_b_batch)
!
!        Form g_kc_db = g_kcbd & reorder 
!
         call mem%alloc(g_bdkc, batch_b%length, wf%n_v, wf%n_o, wf%n_v)
!
         call wf%eri%get_eri_t1('vvov', g_bdkc, first_p=batch_b%first, last_p=batch_b%last)
!
!        :: Term 2. - sum_kcd g_kcbd t_kj^ad c_ci ::
!
!        X_bdki = g_bdkc c_ci 
!
         call mem%alloc(X_bdki, batch_b%length, wf%n_v, wf%n_o, wf%n_o)
!
         call dgemm('N','N',                             &
                     (wf%n_v)*(wf%n_o)*(batch_b%length), &
                     wf%n_o,                             &
                     wf%n_v,                             &
                     one,                                &
                     g_bdkc,                             & ! g_bdk,c
                     (wf%n_v)*(wf%n_o)*(batch_b%length), &
                     c_ai,                               & ! c_c,i 
                     wf%n_v,                             &
                     zero,                               &
                     X_bdki,                             & ! X_bdk,i
                     (wf%n_v)*(wf%n_o)*(batch_b%length))         
!
!        sum_kcd g_kcbd t_kj^ad c_ci = sum_kd (sum_c c_ci g_kcbd) t_kj^ad
!                                    = sum_kd X_bdki t_jk^da
!                                    = sum_kd X_ib_dk t_dk_aj
!
!        Reorder X_bdki to X_ibdk 
!
          call mem%alloc(X_ib_dk, wf%n_o, batch_b%length, wf%n_v, wf%n_o)
          call sort_1234_to_4123(X_bdki, X_ib_dk, batch_b%length, wf%n_v, wf%n_o, wf%n_o)
!
          call mem%dealloc(X_bdki, wf%n_o, wf%n_v, wf%n_o, batch_b%length)
!
!        Reorder t_djak to t_dkaj
!
         call mem%alloc(t_dkaj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
         call squareup_and_sort_1234_to_1432(wf%t2, t_dkaj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!        Calculate rho_ib_aj = - sum_dk X_ib_dk t_dk_aj
!
         call mem%alloc(rho_ibaj, wf%n_o, batch_b%length, wf%n_v, wf%n_o)
!
         call dgemm('N','N',                    &
                     (wf%n_o)*(batch_b%length), &
                     (wf%n_o)*(wf%n_v),         &
                     (wf%n_o)*(wf%n_v),         &
                     -one,                      &
                     X_ib_dk,                   & ! X_ib_dk
                     (wf%n_o)*(batch_b%length), &
                     t_dkaj,                    & ! t_dk_aj
                     (wf%n_o)*(wf%n_v),         &
                     zero,                      &
                     rho_ibaj,                  & ! rho_ib_aj
                     (wf%n_o)*(batch_b%length))
!
         call mem%dealloc(X_ib_dk, wf%n_o, batch_b%length, wf%n_v, wf%n_o)
         call mem%dealloc(t_dkaj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!        Add rho_ibaj (batch over b) ro rho_aibj (full space)
!
!$omp parallel do private(b, j, i, a)
         do a = 1, wf%n_v
            do j = 1, wf%n_o
               do b = 1, batch_b%length
                  do i = 1, wf%n_o
!
                     rho_aibj(a,i,batch_b%first + b - 1,j) = rho_aibj(a,i,batch_b%first + b - 1,j) &
                                                            + rho_ibaj(i,b,a,j)
!
                  enddo
               enddo
            enddo
         enddo
!$omp end parallel do
!
         call mem%dealloc(rho_ibaj, wf%n_o, batch_b%length, wf%n_v, wf%n_o)
!
!        :: Term 3. - sum_kcd g_kcbd t_ik^ca c_dj ::
!
!        sum_d g_bdkc c_dj = sum_d g_ckb_d c_d_j
!
!        Reorder g_bdkc to g_ckbd 
!
         call mem%alloc(g_ckbd, wf%n_v, wf%n_o, batch_b%length, wf%n_v)
!
         call sort_1234_to_4312(g_bdkc, g_ckbd, batch_b%length, wf%n_v, wf%n_o, wf%n_v)
!
         call mem%dealloc(g_bdkc, batch_b%length, wf%n_v, wf%n_o, wf%n_v)
!
!        Form the intermediate X_ckb_j = sum_d g_ckb_d c_d_j
!
         call mem%alloc(X_ckb_j, wf%n_v, wf%n_o, batch_b%length, wf%n_o)
!
         call dgemm('N','N',                             &
                     (wf%n_v)*(wf%n_o)*(batch_b%length), &
                     wf%n_o,                             &
                     wf%n_v,                             &
                     one,                                &
                     g_ckbd,                             & ! g_ckb_d
                     (wf%n_v)*(wf%n_o)*(batch_b%length), &
                     c_ai,                               & ! c_d_j
                     wf%n_v,                             &
                     zero,                               &
                     X_ckb_j,                            & ! X_ckb_j
                     (wf%n_v)*(wf%n_o)*(batch_b%length))
!
!        Order amplitudes as t_ai_ck = t_ik^ca = t_ki^ac (t_ak,ci)
!
         call mem%alloc(t_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
         call squareup_and_sort_1234_to_1432(wf%t2, t_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!        Form rho_aib_j = -sum_kcd g_kcbd t_ik^ca c_dj = sum_ck t_ai_ck X_ckb_j
!
         call mem%alloc(rho_aib_j, wf%n_v, wf%n_o, batch_b%length, wf%n_o)
!
         call dgemm('N','N',                    &
                     (wf%n_v)*(wf%n_o),         &
                     (wf%n_o)*(batch_b%length), &
                     (wf%n_v)*(wf%n_o),         &
                     -one,                      &
                     t_aick,                    & ! t_ai_ck
                     (wf%n_v)*(wf%n_o),         &
                     X_ckb_j,                   & ! X_ck_bj
                     (wf%n_v)*(wf%n_o),         &
                     zero,                      &
                     rho_aib_j,                 & ! rho_aibj
                     (wf%n_v)*(wf%n_o))
!
         call mem%dealloc(X_ckb_j, wf%n_v, wf%n_o, batch_b%length, wf%n_o)
         call mem%dealloc(t_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!        Add rho_aibj to rho_aibj
!
!$omp parallel do private(b, j, i, a)
         do j = 1, wf%n_o
            do b = 1, batch_b%length
               do i = 1, wf%n_o
                  do a = 1, wf%n_v
!
                     rho_aibj(a,i,batch_b%first + b - 1,j) = rho_aibj(a,i,batch_b%first + b - 1,j) &
                                                            + rho_aib_j(a,i,b,j)
!
                 enddo
               enddo
            enddo
         enddo
!$omp end parallel do
!
!        Deallocations for term 3 (keep g_ckb_d = g_kcbd)
!
         call mem%dealloc(rho_aib_j, wf%n_v, wf%n_o, batch_b%length, wf%n_o)
!
!        :: Term 4.  sum_kcd L_kcbd t_ik^ac c_dj ::
!
!        sum_d L_kcbd c_dj = sum_d L_ckb_d c_dj
!
!        Form L_ckb_d = L_kcbd = 2 * g_kcbd - g_kdbc = 2 * g_ckb_d(ckb, d) - g_ckb_d(dkb, c)
!
         call mem%alloc(L_ckbd, wf%n_v, wf%n_o, batch_b%length, wf%n_v)
         call zero_array(L_ckbd, (wf%n_v**2)*wf%n_o*batch_b%length)
!
!        Note: exchange g_ck_bd(dk,bc) -> 4231
!
         call add_4231_to_1234(-one, g_ckbd, L_ckbd, wf%n_v, wf%n_o, batch_b%length, wf%n_v)
         call daxpy((wf%n_v)**2*(wf%n_o)*(batch_b%length), two, g_ckbd, 1, L_ckbd, 1)
!
         call mem%dealloc(g_ckbd, wf%n_v, wf%n_o, batch_b%length, wf%n_v)
!
!        Form the intermediate Y_ckb_j = sum_d L_kcbd c_dj = sum_d L_ckb_d c_dj
!
         call mem%alloc(Y_ckb_j, wf%n_v, wf%n_o, batch_b%length, wf%n_o)
!
         call dgemm('N','N',                             &
                     (wf%n_v)*(wf%n_o)*(batch_b%length), &
                     wf%n_o,                             &
                     wf%n_v,                             &
                     one,                                &
                     L_ckbd,                             & ! L_ckb_d
                     (wf%n_v)*(wf%n_o)*(batch_b%length), &
                     c_ai,                               & ! c_d_j
                     wf%n_v,                             &
                     zero,                               &
                     Y_ckb_j,                            & ! Y_ckb_j
                     (wf%n_v)*(wf%n_o)*(batch_b%length))
!
!        Order amplitudes as t_ai_ck = t_ik^ac = t2(aick)
!
         call mem%alloc(t_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
         call squareup(wf%t2, t_aick, (wf%n_o)*(wf%n_v))
!
!        Form rho_aib_j =  sum_ck t_ai_ck Y_ckb_j
!
         call mem%alloc(rho_aib_j, wf%n_v, wf%n_o, batch_b%length, wf%n_o)
!
         call dgemm('N','N',                    &
                     (wf%n_o)*(wf%n_v),         &
                     (wf%n_o)*(batch_b%length), &
                     (wf%n_o)*(wf%n_v),         &
                     one,                       &
                     t_aick,                    & ! t_ai_ck
                     (wf%n_o)*(wf%n_v),         &
                     Y_ckb_j,                   & ! Y_ck_bj
                     (wf%n_o)*(wf%n_v),         &
                     zero,                      &
                     rho_aib_j,                 & ! rho_ai_bj
                     (wf%n_o)*(wf%n_v))
!
         call mem%dealloc(Y_ckb_j, wf%n_v, wf%n_o, batch_b%length, wf%n_o)
         call mem%dealloc(t_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!        Add rho_aib_j to rho_aibj
!
!$omp parallel do private(b, j, i, a)
         do j = 1, wf%n_o
            do b = 1, batch_b%length
               do i = 1, wf%n_o
                  do a = 1, wf%n_v
!
                     rho_aibj(a,i,batch_b%first + b - 1,j) = rho_aibj(a,i,batch_b%first + b - 1,j) &
                                                            + rho_aib_j(a,i,b,j)
!
                  enddo
               enddo
            enddo
         enddo
!$omp end parallel do
!
!        Deallocations for term 4 (keep L_ckb_d = L_kcbd)
!
         call mem%dealloc(rho_aib_j, wf%n_v, wf%n_o, batch_b%length, wf%n_o)
!
!        :: Term 5.  sum_kcd L_kcbd t_ij^ad c_ck ::
!
!        Form the intermediate X_1,bd = sum_ck c_ck L_kcbd = sum_ck c_1,ck L_ckb_d
!
!        Note: c_ai is interpreted as c_1,ai in the matrix multiplication
!
         call mem%alloc(X_bd, 1, (batch_b%length)*(wf%n_v))
!
         call dgemv('T',                        &
                     (wf%n_o)*(wf%n_v),         &
                     (batch_b%length)*(wf%n_v), &
                     one,                       &
                     L_ckbd,                    & ! L_ck_bd
                     (wf%n_o)*(wf%n_v),         &
                     c_ai,                      & ! c_ck
                     1,                         &
                     zero,                      &
                     X_bd,                      & ! X_bd
                     1)
!
!        Order amplitudes as t_ai_jd = t_ij^ad
!
         call mem%alloc(t_aijd, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
         call squareup_and_sort_1234_to_1243(wf%t2, t_aijd, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!        Form rho_b_aij =  sum_kcd L_kcbd t_ij^ad c_ck
!                       =  sum_d X_bd t_d_aij
!
!        Note: X_bd is interpreted as X_b_d in the matrix multiplication
!
         call mem%alloc(rho_baij, batch_b%length, wf%n_v, wf%n_o, wf%n_o)
!
         call dgemm('N','T',               &
                     batch_b%length,       &
                     (wf%n_v)*(wf%n_o)**2, &
                     wf%n_v,               &
                     one,                  &
                     X_bd,                 & ! X_b_d
                     batch_b%length,       &
                     t_aijd,               & ! t_aij_d
                     (wf%n_v)*(wf%n_o)**2, &
                     zero,                 &
                     rho_baij,             & ! rho_b_aij
                     batch_b%length)
!
         call mem%dealloc(X_bd, 1, (batch_b%length)*(wf%n_v))
         call mem%dealloc(t_aijd, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
!        Add rho_baij to rho_aibj
!
!$omp parallel do private(b, j, i, a)
         do j = 1, wf%n_o
            do i = 1, wf%n_o
               do a = 1, wf%n_v
                  do b = 1, batch_b%length
!
                     rho_aibj(a,i,batch_b%first + b - 1,j) = rho_aibj(a,i,batch_b%first + b - 1,j) &
                                                            + rho_baij(b,a,i,j)
!
                  enddo
               enddo
            enddo
         enddo
!$omp end parallel do
!
!        Final deallocations in batching loop
!
         call mem%dealloc(rho_baij, batch_b%length, wf%n_v, wf%n_o, wf%n_o)
         call mem%dealloc(L_ckbd, wf%n_v, wf%n_o,batch_b%length, wf%n_v)
!
      enddo ! End of batches over b
!
      call timer%turn_off()
!
   end subroutine jacobian_ccsd_d2_ccsd
!
!
    module subroutine jacobian_ccsd_e2_ccsd(wf, rho_aibj, c_aick)
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
   end subroutine jacobian_ccsd_e2_ccsd
!
!
   module subroutine jacobian_ccsd_f2_ccsd(wf, rho_aibj, c_aibj)
!!
!!    Jacobian CCSD F2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!       rho_aibj^F2 =   - sum_ckdl t_ai,dj * L_kc,ld * c_bl,ck
!!                       - sum_ckdl t_ai_bl * L_kc,ld * c_ck,dj
!!
!!    L_kc,ld = 2*g_kc,ld - g_kd,lc = 2*g_kcld(kc,ld) - 2*g_kcld(kd,lc)
!!
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
      call wf%eri%get_eri_t1('ovov', g_kcld)
!
!     Construct L_ckdl reordered as L_dlck
!
      call mem%alloc(L_dlck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call zero_array(L_dlck, wf%n_t1**2)
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
   end subroutine jacobian_ccsd_f2_ccsd
!
!
   module subroutine jacobian_ccsd_g2_ccsd(wf, rho_aibj, c_aibj)
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
!!    in the routine save_jacobian_g2_intermediates_ccsd
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
   end subroutine jacobian_ccsd_g2_ccsd
!
!
   module subroutine jacobian_ccsd_h2_ccsd(wf, rho_aibj, c_aibj)
!!
!!    Jacobian CCSD H2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!       rho_aibj^H2 =  sum_ckdl t_ci,ak * g_kc,ld * c_bl,dj
!!                     + sum_ckdl t_cj,al * g_kc,ld * c_bk,di
!!                   = sum_dl Y_aild c_bldj
!!                     + sum_dk Y_ajkd c_bkdi 
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o)               :: rho_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)   :: c_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: c_ldbj
!
      real(dp), dimension(:,:,:,:), allocatable :: Y_aild
      real(dp), dimension(:,:,:,:), allocatable :: Y_ajkd
!
      real(dp), dimension(:,:,:,:), allocatable :: rho_ajbi
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian CCSD H2 transformation', pl='verbose')
      call timer%turn_on()
!
!     Term 1: Y_aild c_bldj
!
      call mem%alloc(Y_aild, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call wf%jacobian_h2_intermediate_voov_1%open_('read', 'rewind')
      call wf%jacobian_h2_intermediate_voov_1%read_(Y_aild, wf%n_t1**2)
      call wf%jacobian_h2_intermediate_voov_1%close_('keep')
!
      call mem%alloc(c_ldbj, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call sort_1234_to_2314(c_aibj, c_ldbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N', &
                  wf%n_v*wf%n_o, &
                  wf%n_v*wf%n_o, &
                  wf%n_v*wf%n_o, &
                  one,           &
                  Y_aild,        &
                  wf%n_v*wf%n_o, &
                  c_ldbj,        &
                  wf%n_v*wf%n_o, &
                  one,           &
                  rho_aibj,      &
                  wf%n_v*wf%n_o)
!
      call mem%dealloc(Y_aild, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
!     Term 2: Y_ajkd c_bkdi
!
!     Pretend c_ldbj is c_kdbi
!
      call mem%alloc(Y_ajkd, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call wf%jacobian_h2_intermediate_voov_2%open_('read', 'rewind')
      call wf%jacobian_h2_intermediate_voov_2%read_(Y_ajkd, wf%n_t1**2)
      call wf%jacobian_h2_intermediate_voov_2%close_('keep')
!
      call mem%alloc(rho_ajbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N', &
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
   end subroutine jacobian_ccsd_h2_ccsd
!
!
   module subroutine jacobian_ccsd_i2_ccsd(wf, rho_aibj, c_aibj)
!!
!!    Jacobian CCSD I2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    rho_aibj^I2 =  sum_c F_bc * c_ai,cj - sum_k F_jk * c_ai,bk
!!                   + sum_ck L_bj,kc * c_ai,ck
!!                   - sum_ck ( g_kc,bj * c_ak,ci + g_ki,bc * c_ak,cj )
!!
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
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian CCSD I2 transformation', pl='verbose')
      call timer%turn_on()
!
!     :: sum_c F_bc c_ai,cj ::
!
!     rho_bjai =+ F_bc c_cjai
!
      call dgemm('N','N',                 &
                  wf%n_v,                 &
                  (wf%n_v)*(wf%n_o)**2,   &
                  wf%n_v,                 &
                  one,                    &
                  wf%fock_ab,             & ! F_b,c
                  wf%n_v,                 &
                  c_aibj,                 & ! c_c,jai
                  wf%n_v,                 &
                  one,                    &
                  rho_aibj,               & ! rho_b,jai -> will be (ai,bj)-symmetrized 
                  wf%n_v)
!
!     :: - sum_k F_jk * c_aibk  ::
!
!     rho_aibj += - sum_k F_jk * c_aibk = - sum_k c_aibk F_ij(k,j)^T
!
      call dgemm('N', 'N',                &
                  (wf%n_o)*((wf%n_v)**2), &
                  wf%n_o,                 &
                  wf%n_o,                 &
                  -one,                   &
                  c_aibj,                 & ! c_aib_k
                  (wf%n_o)*((wf%n_v)**2), &
                  wf%fock_ij,             & ! F_k_j
                  wf%n_o,                 &
                  one,                    &
                  rho_aibj,               & ! rho_aib_j
                  (wf%n_o)*((wf%n_v)**2))
!
!     rho_aibj =+ c_aick L_bjkc = c_aick L_ck,bj 
!
      call mem%alloc(g_bjkc, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%eri%get_eri_t1('voov', g_bjkc)
!
!     rho_bjai =- g_bjkc c_akci = g_bj,kc c_kc,ai
!
      call mem%alloc(c_kcai, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call sort_1234_to_2314(c_aibj, c_kcai, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',             &
                  (wf%n_o)*(wf%n_v),   &
                  (wf%n_o)*(wf%n_v),   &
                  (wf%n_o)*(wf%n_v),   &
                  -one,                &
                  g_bjkc,              &
                  (wf%n_o)*(wf%n_v),   &
                  c_kcai,              &
                  (wf%n_o)*(wf%n_v),   &
                  one,                 &
                  rho_aibj,            & ! rho_bj,ai -> will be (ai,bj)-symmetrized
                  (wf%n_o)*(wf%n_v))
!
!     Construct L_bjkc ordered as L_ckbj = 2 g_bjkc - g_bckj 
!
      call mem%alloc(L_ckbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call zero_array(L_ckbj, (wf%n_t1)**2)
!
      call add_3421_to_1234(two, g_bjkc, L_ckbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_bjkc, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(g_bckj, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call wf%eri%get_eri_t1('vvoo', g_bckj)
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
   end subroutine jacobian_ccsd_i2_ccsd
!
!
   module subroutine jacobian_ccsd_j2_ccsd(wf, rho_abij, c_abij)
!!
!!    Jacobian CCSD J2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!       rho_abij^J2 =    sum_ckld t_ci,dj * g_kc,ld * c_ak,bl
!!                       + sum_ckdl t_ak,bl * g_kc,ld * c_ci,dj
!!                   =    sum_ckld Y_klij * c_ak,bl
!!                       + sum_ckdl t_ak,bl * g_kc,ld * c_ci,dj
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
      call wf%eri%get_eri_t1('ovov', g_kcld)
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
   end subroutine jacobian_ccsd_j2_ccsd
!
!
   module subroutine jacobian_ccsd_k2_ccsd(wf, rho_abij, c_abij)
!!
!!    Jacobian CCSD K2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    rho_abij^K2 =    sum_kl g_kilj * c_akbl
!!
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
      call wf%eri%get_eri_t1('oooo', g_kilj)
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
   end subroutine jacobian_ccsd_k2_ccsd
!
!
   module subroutine save_jacobian_c2_intermediates_ccsd(wf)
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
      call wf%eri%get_eri_t1('ooov', g_ljkc)
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
   end subroutine save_jacobian_c2_intermediates_ccsd
!
!
   module subroutine save_jacobian_d2_intermediate_ccsd(wf)
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
      call mem%alloc(X_kijb_full, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call zero_array(X_kijb_full, (wf%n_o**3)*wf%n_v)
!
!     Order amplitudes as t_ij_cd = t_ij^cd = t_ci_dj
!
      call mem%alloc(t_ijcd, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call squareup_and_sort_1234_to_2413(wf%t2, t_ijcd, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Initialize batching variable
!
      req0 = wf%n_v*wf%n_o*wf%eri%n_J
      req1 = wf%n_v*wf%eri%n_J + 2*(wf%n_o)*(wf%n_v)**2 + 2*(wf%n_o)**3
!
      batch_b = batching_index(wf%n_v)
      call mem%batch_setup(batch_b, req0, req1)
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
         call wf%eri%get_eri_t1('vvov', g_bdkc, first_p=batch_b%first, last_p=batch_b%last)
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
   end subroutine save_jacobian_d2_intermediate_ccsd
!
!
   module subroutine save_jacobian_e2_intermediate_ccsd(wf)
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
      call wf%eri%get_eri_t1('ovov', g_ldkc)
!
      call mem%alloc(L_dlck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call zero_array(L_dlck, (wf%n_v**2)*(wf%n_o**2))
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
   end subroutine save_jacobian_e2_intermediate_ccsd
!
!
   module subroutine save_jacobian_g2_intermediates_ccsd(wf)
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
      call wf%eri%get_eri_t1('ovov', g_ldkc)
!
      call mem%alloc(L_dlck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call zero_array(L_dlck, (wf%n_v**2)*(wf%n_o**2))
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
   end subroutine save_jacobian_g2_intermediates_ccsd
!
!
   module subroutine save_jacobian_h2_intermediates_ccsd(wf)
!!
!!    Save jacobian h2 intermediates
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2019
!!
!!    Constructs the intermediates 
!!
!!       Y_aild = t_ciak g_kcld
!!       Y_ajkd = t_cjal g_kcld
!!
!!    The intermediates are stored in the files:
!!
!!       jacobian_h2_intermediate_voov_1
!!       jacobian_h2_intermediate_voov_2
!!
!!    which are wavefunction variables
!!
      implicit none
!
      class(ccsd) :: wf
!
      type(timings), allocatable :: timer
!
      real(dp), dimension(:,:,:,:), allocatable :: g_kcld
      real(dp), dimension(:,:,:,:), allocatable :: g_lckd
      real(dp), dimension(:,:,:,:), allocatable :: t_aikc
      real(dp), dimension(:,:,:,:), allocatable :: Y_aild
      real(dp), dimension(:,:,:,:), allocatable :: Y_ajkd
     
      timer = timings('Jacobian CCSD H2 intermediate construction', pl='verbose')
      call timer%turn_on()
!
      call mem%alloc(t_aikc, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call squareup_and_sort_1234_to_1423(wf%t2, t_aikc, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(g_kcld, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call wf%eri%get_eri_t1('ovov', g_kcld)
!
!     Y_aild = t_ciak g_kcld
!
      call mem%alloc(Y_aild, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N', 'N', &
                  wf%n_v*wf%n_o, &
                  wf%n_v*wf%n_o, &
                  wf%n_v*wf%n_o, &
                  one,           &
                  t_aikc,        &
                  wf%n_v*wf%n_o, &
                  g_kcld,        &
                  wf%n_v*wf%n_o, &
                  zero,          &
                  Y_aild,        &
                  wf%n_v*wf%n_o)
!
      wf%jacobian_h2_intermediate_voov_1 = sequential_file('jacobian_h2_intermediate_voov_1_ccsd')
      call wf%jacobian_h2_intermediate_voov_1%open_('write', 'rewind')
!
      call wf%jacobian_h2_intermediate_voov_1%write_(Y_aild, wf%n_t1**2)
!
      call mem%dealloc(Y_aild, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call wf%jacobian_h2_intermediate_voov_1%close_('keep')
!
!     Y_ajkd = t_alcj g_kcld
!
!     Note: pretend that t_aikc is t_ajcl
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
                  t_aikc,        & ! t_ajlc
                  wf%n_v*wf%n_o, &
                  g_lckd,        & 
                  wf%n_v*wf%n_o, &
                  zero,          &
                  Y_ajkd,        &
                  wf%n_v*wf%n_o)       
!
      call mem%dealloc(g_lckd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%dealloc(t_aikc, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      wf%jacobian_h2_intermediate_voov_2 = sequential_file('jacobian_h2_intermediate_voov_2_ccsd')
      call wf%jacobian_h2_intermediate_voov_2%open_('write', 'rewind')
!
      call wf%jacobian_h2_intermediate_voov_2%write_(Y_ajkd, wf%n_t1**2)
!
      call mem%dealloc(Y_ajkd, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call wf%jacobian_h2_intermediate_voov_2%close_('keep')
!
      call timer%turn_off()
!
   end subroutine save_jacobian_h2_intermediates_ccsd
!
!
   module subroutine save_jacobian_j2_intermediate_ccsd(wf)
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
      call wf%eri%get_eri_t1('ovov', g_kcld)
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
   end subroutine save_jacobian_j2_intermediate_ccsd
!
!
end submodule jacobian_ccsd
