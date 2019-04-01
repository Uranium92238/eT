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
submodule (ccsd_class) jacobian_transpose_ccsd
!
!!
!!    Jacobian transpose submodule (CCSD)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Routines for the linear transform of trial
!!    vectors by the transpose of the Jacobian matrix
!!
!!    σ_i = A^T * b_i,
!!
!!    where
!!
!!    A_μ,ν = < μ | exp(-T) [H, τ_ν] exp(T) | R >.
!!
!!    Transfered to the current eT program from the first version
!!    of eT by Andreas Skeidsvoll and Sarai D. Folkestad, 2018.
!!
!
   implicit none
!
!
contains
!
!
   module subroutine jacobian_transpose_transform_trial_vector_ccsd(wf, c_i)
!!
!!    Jacobian transform trial vector
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep 2018
!!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes) :: c_i
!
      call wf%jacobian_transpose_ccsd_transformation(c_i)
!
   end subroutine jacobian_transpose_transform_trial_vector_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_transformation_ccsd(wf, b)
!!
!!    Jacobian transpose transformation (CCSD)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Calculates the transpose Jacobian transformation, i.e., the transformation
!!    by the transpose of the Jacobian matrix
!!
!!       A_mu,nu = < mu | exp(-T) [H, tau_nu] exp(T) | R >.
!!
!!    The transformation is performed as sigma^T = b^T A, where b is the vector
!!    sent to the routine. On exit, the vector b is equal to sigma (the transformed
!!    vector).
!!
      implicit none
!
      class(ccsd) :: wf
!
!     Incoming vector b
!
      real(dp), dimension(wf%n_es_amplitudes) :: b
!
      real(dp), dimension(:,:), allocatable :: b_ai
!
!     Local unpacked and reordered vectors
!
      real(dp), dimension(:,:,:,:), allocatable :: b_aibj ! Unpacked b_aibj
      real(dp), dimension(:,:,:,:), allocatable :: b_abij ! b_aibj, reordered
!
      real(dp), dimension(:,:,:,:), allocatable :: sigma_abij     ! sigma_aibj, reordered
!
      real(dp), dimension(:,:), allocatable :: sigma_ai
      real(dp), dimension(:,:,:,:), allocatable :: sigma_aibj
!
!     Indices
!
      type(timings) :: jacobian_transpose_timer
!
      call jacobian_transpose_timer%init('jacobian transpose')
      call jacobian_transpose_timer%start()
!
      call mem%alloc(sigma_ai, wf%n_v, wf%n_o)
      call mem%alloc(b_ai, wf%n_v, wf%n_o)
!
      call dcopy(wf%n_t1, b(1:wf%n_t1), 1, b_ai, 1)
!
!     Calculate and add the CCS contributions to the
!     singles transformed vector

      call wf%jacobian_transpose_ccs_a1(sigma_ai, b_ai)
      call wf%jacobian_transpose_ccs_b1(sigma_ai, b_ai)
!
!     Calculate and add the CCSD contributions to the
!     singles transformed vector
!
      call wf%jacobian_transpose_ccsd_a1(sigma_ai, b_ai)
      call wf%jacobian_transpose_ccsd_b1(sigma_ai, b_ai)
!
      call mem%alloc(b_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(b(wf%n_t1+1:wf%n_t1+wf%n_t2), b_aibj, wf%n_v*wf%n_o)
!
      call wf%jacobian_transpose_ccsd_c1(sigma_ai, b_aibj)
      call wf%jacobian_transpose_ccsd_d1(sigma_ai, b_aibj)
      call wf%jacobian_transpose_ccsd_e1(sigma_ai, b_aibj)
      call wf%jacobian_transpose_ccsd_f1(sigma_ai, b_aibj)
      call wf%jacobian_transpose_ccsd_g1(sigma_ai, b_aibj)
!
      call dcopy(wf%n_t1, sigma_ai, 1, b(1:wf%n_t1), 1)
      call mem%dealloc(sigma_ai, wf%n_v, wf%n_o)
!
!     Add the CCSD contributions to the doubles vector arising from
!     the incoming singles vector
!
      call mem%alloc(sigma_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%jacobian_transpose_ccsd_a2(sigma_aibj, b_ai)
!
      call mem%dealloc(b_ai, wf%n_v, wf%n_o)
!
      call wf%jacobian_transpose_ccsd_b2(sigma_aibj, b_aibj)
      call wf%jacobian_transpose_ccsd_c2(sigma_aibj, b_aibj)
      call wf%jacobian_transpose_ccsd_d2(sigma_aibj, b_aibj)
      call wf%jacobian_transpose_ccsd_e2(sigma_aibj, b_aibj)
      call wf%jacobian_transpose_ccsd_f2(sigma_aibj, b_aibj)
      call wf%jacobian_transpose_ccsd_g2(sigma_aibj, b_aibj)
!
!     Last two terms are already symmetric (h2 and i2). Perform the symmetrization
!     sigma_aibj = P_ij^ab sigma_aibj now, for convenience
!
      call symmetric_sum(sigma_aibj, wf%n_v*wf%n_o)
!
!     In preparation for last two terms, reorder b_aibj to b_abij
!
      call mem%alloc(b_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_1234_to_1324(b_aibj, b_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(b_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(sigma_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
!     Add the last two terms
!
      call wf%jacobian_transpose_ccsd_h2(sigma_abij, b_abij)
      call wf%jacobian_transpose_ccsd_i2(sigma_abij, b_abij)
!
      call mem%dealloc(b_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call add_1324_to_1234(one, sigma_abij, sigma_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(sigma_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
!     Overwrite the incoming doubles b vector
!
      call packin(b(wf%n_t1+1:wf%n_t1+wf%n_t2), sigma_aibj, wf%n_v*wf%n_o)
!
      call mem%dealloc(sigma_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call jacobian_transpose_timer%freeze()
      call jacobian_transpose_timer%switch_off()
!
   end subroutine jacobian_transpose_ccsd_transformation_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_a1_ccsd(wf, sigma_ai, b_ai)
!!
!!    Jacobian transpose CCSD A1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Calculates the A1 term,
!!
!!       sum_ckdl b_ck L_iald u_kl^cd,
!!
!!    abd adds it to the transformed vector sigma_ai.
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)     :: b_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)  :: sigma_ai
!
      real(dp), dimension(:,:,:,:), allocatable :: u_ldck
!
      real(dp), dimension(:,:,:,:), allocatable :: t_ldck
!
      real(dp), dimension(:), allocatable :: X_ld ! An intermediate, see below
!
      real(dp), dimension(:,:,:,:), allocatable :: g_iald ! g_iald
      real(dp), dimension(:,:,:,:), allocatable :: L_aild
!
!     Form u_ldck = u_kl^cd = 2 * t_kl^cd - t_lk^cd = 2 * t_ldck(l,d,c,k) - t_ldck(k,d,c,l)
!
      call mem%alloc(t_ldck, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      call squareup_and_sort_1234_to_4312(wf%t2, t_ldck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(u_ldck, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      u_ldck = zero
!
      call daxpy((wf%n_v**2)*(wf%n_o**2), two, t_ldck, 1, u_ldck, 1)
      call add_4231_to_1234(-one, t_ldck, u_ldck, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      call mem%dealloc(t_ldck, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
!     Form the intermediate X_ld = sum_ck u_ld_ck b_ck
!
      call mem%alloc(X_ld, (wf%n_v)*(wf%n_o))
!
      call dgemm('N','N',            &
                  (wf%n_v)*(wf%n_o), &
                  1,                 &
                  (wf%n_v)*(wf%n_o), &
                  one,               &
                  u_ldck,            & ! u_ld_ck
                  (wf%n_v)*(wf%n_o), &
                  b_ai,              & ! "b_ck"
                  (wf%n_v)*(wf%n_o), &
                  zero,              &
                  X_ld,              &
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(u_ldck, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
!     Form L_iald = 2 * g_iald - g_idla
!                 = 2 * g_iald(i,a,l,d) - g_iald(i,d,l,a)
!
      call mem%alloc(g_iald, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_ovov(g_iald)
!
      call mem%alloc(L_aild, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      L_aild = zero
!
      call add_2134_to_1234(two, g_iald, L_aild, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call add_2431_to_1234(-one, g_iald, L_aild, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call mem%dealloc(g_iald, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     Calculate and add sum_ckdl b_ck L_iald u_kl^cd
!                       = sum_ld L_ai_ld X_ld
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  1,                 &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  L_aild,            &
                  (wf%n_o)*(wf%n_v), &
                  X_ld,              &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  sigma_ai,          & ! "sigma_ai"
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(L_aild, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call mem%dealloc(X_ld, (wf%n_o)*(wf%n_v))
!
   end subroutine jacobian_transpose_ccsd_a1_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_b1_ccsd(wf, sigma_ai, b_ai)
!!
!!    Jacobian transpose CCSD B1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Calculates the B1 term,
!!
!!       - sum_ckdl (b_al L_kcid t_kl^cd + b_ci L_ldka t_kl^cd),
!!
!!    abd adds it to the transformed vector sigma_ai.
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: b_ai
!
      real(dp), dimension(:,:,:,:), allocatable :: g_kcid
      real(dp), dimension(:,:,:,:), allocatable :: L_kcdi
!
      real(dp), dimension(:,:,:,:), allocatable :: t_lkcd ! t_kl^cd
!
      real(dp), dimension(:,:,:,:), allocatable :: L_aldk ! L_ldka
      real(dp), dimension(:,:,:,:), allocatable :: t_ldkc ! t_kl^cd
!
      real(dp), dimension(:,:), allocatable :: X_li ! An intermediate, term 1
      real(dp), dimension(:,:), allocatable :: X_ac ! An intermediate, term 2
!
!
!     :: Term 1. - sum_ckdl b_al L_kcid t_kl^cd ::
!
!     Form L_kcdi = L_kcid = 2 * g_kcid - g_kdic
!                          = 2 * g_kcid(k,c,i,d) - g_kcid(k,d,i,c)
!
      call mem%alloc(g_kcid, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_ovov(g_kcid)
!
      call mem%alloc(L_kcdi, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      L_kcdi = zero
!
      call add_1243_to_1234(two, g_kcid, L_kcdi, wf%n_o, wf%n_v, wf%n_v, wf% n_o)
      call add_1342_to_1234(-one, g_kcid, L_kcdi, wf%n_o, wf%n_v, wf%n_v, wf% n_o)
!
      call mem%dealloc(g_kcid, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     Form t_lkcd = t_kl^cd
!
      call mem%alloc(t_lkcd, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
      call squareup_and_sort_1234_to_4213(wf%t2, t_lkcd, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Calculate the intermediate X_li = sum_kcd t_l_kcd L_kcd_i
!
      call mem%alloc(X_li, wf%n_o, wf%n_o)
!
      call dgemm('N','N',               &
                  wf%n_o,               &
                  wf%n_o,               &
                  (wf%n_o)*(wf%n_v)**2, &
                  one,                  &
                  t_lkcd,               & ! t_l_kcd
                  wf%n_o,               &
                  L_kcdi,               & ! L_kcd_i
                  (wf%n_o)*(wf%n_v)**2, &
                  zero,                 &
                  X_li,                 &
                  wf%n_o)
!
!     Add - sum_ckdl b_al L_kcid t_kl^cd = sum_l b_al X_li
!
      call dgemm('N','N',   &
                  wf%n_v,   &
                  wf%n_o,   &
                  wf%n_o,   &
                  -one,     &
                  b_ai,     & ! b_al
                  wf%n_v,   &
                  X_li,     &
                  wf%n_o,   &
                  one,      &
                  sigma_ai, &
                  wf%n_v)
!
      call mem%dealloc(X_li, wf%n_o, wf%n_o)
!
!     :: Term 2. - sum_ckdl b_ci L_ldka t_kl^cd ::
!
!     Form L_aldk = L_ldka = L_kcdi(l,d,a,k)
!
      call mem%alloc(L_aldk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call sort_1234_to_3124(L_kcdi, L_aldk, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      call mem%dealloc(L_kcdi, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
!     Form t_ldkc = t_kl^cd = t_lkcd(l,k,c,d)
!
      call mem%alloc(t_ldkc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call sort_1234_to_1423(t_lkcd, t_ldkc, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
      call mem%dealloc(t_lkcd, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
!     Calculate the intermediate X_ac = sum_ldk L_a_ldk t_ldk_c
!
      call mem%alloc(X_ac, wf%n_v, wf%n_v)
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  L_aldk,               & ! L_a_ldk
                  wf%n_v,               &
                  t_ldkc,               & ! t_ldk_c
                  (wf%n_v)*(wf%n_o)**2, &
                  zero,                 &
                  X_ac,                 &
                  wf%n_v)
!
      call mem%dealloc(L_aldk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(t_ldkc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     Add - sum_ckdl b_ci L_ldka t_kl^cd = - sum_c X_ac b_ci
!
      call dgemm('N','N',   &
                  wf%n_v,   &
                  wf%n_o,   &
                  wf%n_v,   &
                  -one,     &
                  X_ac,     &
                  wf%n_v,   &
                  b_ai,     & ! b_ci
                  wf%n_v,   &
                  one,      &
                  sigma_ai, &
                  wf%n_v)
!
      call mem%dealloc(X_ac, wf%n_v, wf%n_v)
!
   end subroutine jacobian_transpose_ccsd_b1_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_c1_ccsd(wf, sigma_ai, b_aibj)
!!
!!    Jacobian transpose CCSD C1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Calculates the C1 term,
!!
!!       sum_cdl b_cidl g_dlca - sum_kdl b_akdl g_dlik,
!!
!!    and adds it to the transformed vector sigma_ai.
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o)                 :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: g_dlca
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ikdl ! g_dlik = g_ikdl
!
!     Batching variables
!
      integer :: rec0, rec1
      integer :: current_a_batch
!
      type(batching_index) :: batch_a
!
!     :: Term 1. sum_cdl b_cidl g_dlca ::
!
!     Prepare batching over index a
!
      rec0 = wf%n_o*wf%n_v*wf%integrals%n_J
      rec1 = wf%n_v*wf%integrals%n_J + wf%n_v**2*wf%n_o
!
      call batch_a%init(wf%n_v)
      call mem%batch_setup(batch_a, rec0, rec1)
!
!     Loop over the a-batches
!
      do current_a_batch = 1, batch_a%num_batches
!
!        For each batch, get the limits for the a index
!
         call batch_a%determine_limits(current_a_batch)
!
!        Form g_dlca
!
         call mem%alloc(g_dlca, wf%n_v, wf%n_o, wf%n_v, batch_a%length)
!
         call wf%get_vovv(g_dlca,         &
                           1,             &
                           wf%n_v,        &
                           1,             &
                           wf%n_o,        &
                           1,             &
                           wf%n_v,        &
                           batch_a%first, &
                           batch_a%last)
!
!        Add sum_dlc g_dlc_a^T b_dlc_i
!
         call dgemm('T','N',                    &
                     batch_a%length,            &
                     wf%n_o,                    &
                     (wf%n_o)*(wf%n_v)**2,      &
                     one,                       &
                     g_dlca,                    & ! g_dlc_a
                     (wf%n_o)*(wf%n_v)**2,      &
                     b_aibj,                    & ! b_dlc_i
                     (wf%n_o)*(wf%n_v)**2,      &
                     one,                       &
                     sigma_ai(batch_a%first,1), &
                     wf%n_v)
!
         call mem%dealloc(g_dlca, wf%n_v, wf%n_o, wf%n_v, batch_a%length)
!
      enddo ! End of batches over a
!
!     :: Term 2. - sum_kdl b_akdl g_dlik = - sum_kdl b_akdl g_ikdl ::
!
      call mem%alloc(g_ikdl, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%get_oovo(g_ikdl)
!
      call dgemm('N','T',               &
                  wf%n_v,               &
                  wf%n_o,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  -one,                 &
                  b_aibj,               & ! b_a_kdl
                  wf%n_v,               &
                  g_ikdl,               & ! g_i_kdl
                  (wf%n_o),             &
                  one,                  &
                  sigma_ai,             &
                  wf%n_v)
!
      call mem%dealloc(g_ikdl, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine jacobian_transpose_ccsd_c1_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_d1_ccsd(wf, sigma_ai, b_aibj)
!!
!!    Jacobian transpose CCSD D1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Calculates the D1 term,
!!
!!       - sum_ckdl (b_ckal F_id t_kl^cd + b_ckdi F_la t_kl^cd),
!!
!!    and adds it to the transformed vector sigma_ai.
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o)                 :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: t_lckd ! t_kl^cd
!
      real(dp), dimension(:,:), allocatable :: X_ad ! An intermediate, term 1
      real(dp), dimension(:,:), allocatable :: X_li ! An intermediate, term 2
!
!     :: Term 1. - sum_ckdl b_ckal F_id t_kl^cd ::
!
!     Read amplitudes and order as t_lckd = t_kl^cd
!
      call mem%alloc(t_lckd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call squareup_and_sort_1234_to_2341(wf%t2, t_lckd, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Form the intermediate X_ad = sum_ckl b_a_lck t_lck_d = sum_ckl b_ckal t_kl^cd
!
      call mem%alloc(X_ad, wf%n_v, wf%n_v)
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  b_aibj,               &  ! b_a_lck = b_alck = b_ckal
                  wf%n_v,               &
                  t_lckd,               & ! t_lck_d
                  (wf%n_v)*(wf%n_o)**2, &
                  zero,                 &
                  X_ad,                 &
                  wf%n_v)
!
!     Add - sum_ckdl b_ckal F_id t_kl^cd
!           = - sum_d X_ad F_id
!           = - sum_d X_ad F_i_a^T(d,i)
!
      call dgemm('N','T',     &
                  wf%n_v,     &
                  wf%n_o,     &
                  wf%n_v,     &
                  -one,       &
                  X_ad,       &
                  wf%n_v,     &
                  wf%fock_ia, & ! F_i_a
                  wf%n_o,     &
                  one,        &
                  sigma_ai,   &
                  wf%n_v)
!
      call mem%dealloc(X_ad, wf%n_v, wf%n_v)
!
!     :: Term 2. - sum_ckdl b_ckdi F_la t_kl^cd
!
!     Form the intermediate X_li = sum_ckd t_l_ckd b_ckd_i  = sum_ckd b_ckdi t_kl^cd
!
!     Note: we interpret b_aibj as b_aib_j, such that b_aib_j(ckd, i) = b_ckdi
!           we interpret t_lckd as t_l_ckd, such that t_l_ckd(l,ckd) = t_kl^cd
!
      call mem%alloc(X_li, wf%n_o, wf%n_o)
!
      call dgemm('N','N',               &
                  wf%n_o,               &
                  wf%n_o,               &
                  (wf%n_o)*(wf%n_v)**2, &
                  one,                  &
                  t_lckd,               & ! t_l_ckd
                  (wf%n_o),             &
                  b_aibj,               & ! b_ckd_i
                  (wf%n_o)*(wf%n_v)**2, &
                  zero,                 &
                  X_li,                 &
                  wf%n_o)
!
      call mem%dealloc(t_lckd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     Add - sum_ckdl b_ckdi F_la t_kl^cd = - sum_l F_la X_li = - sum_l F_i_a^T(a,l) X_li(l,i)
!
      call dgemm('T','N',     &
                  wf%n_v,     &
                  wf%n_o,     &
                  wf%n_o,     &
                  -one,       &
                  wf%fock_ia, &
                  wf%n_o,     &
                  X_li,       &
                  wf%n_o,     &
                  one,        &
                  sigma_ai,   &
                  wf%n_v)
!
      call mem%dealloc(X_li, wf%n_o, wf%n_o)
!
   end subroutine jacobian_transpose_ccsd_d1_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_e1_ccsd(wf, sigma_ai, b_aibj)
!!
!!    Jacobian transpose CCSD E1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Calculates the E1 term,
!!
!!       sum_ckdle (b_ckdi L_dale t_kl^ce + b_ckdl L_deia t_kl^ce)
!!      -sum_ckdlm (b_ckal L_ilmd t_km^cd + b_ckdl L_mlia t_km^cd)
!!
!!    and adds it to the transformed vector sigma_ai.
!!
!!    The routine adds the third and forth terms first.
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o)                 :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: t_dmck ! t_km^cd
      real(dp), dimension(:,:,:,:), allocatable :: t_elck ! t_lk^ec
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ilmd
      real(dp), dimension(:,:,:,:), allocatable :: g_mlia
!
      real(dp), dimension(:,:,:,:), allocatable :: L_ildm
      real(dp), dimension(:,:,:,:), allocatable :: L_aiml ! L_mlia
!
      real(dp), dimension(:,:,:,:), allocatable :: X_ilck ! An intermediate, term 3
!
      real(dp), dimension(:,:), allocatable :: X_ml ! An intermediate, term 4
!
      real(dp), dimension(:,:,:,:), allocatable :: g_dale
      real(dp), dimension(:,:,:,:), allocatable :: L_aeld ! L_dale
!
      real(dp), dimension(:,:,:,:), allocatable :: g_deia
      real(dp), dimension(:,:,:,:), allocatable :: L_aied ! L_deia
!
      real(dp), dimension(:,:,:,:), allocatable :: X_eldi ! An intermediate, term 1
!
      real(dp), dimension(:,:), allocatable :: X_de ! An intermediate, term 2
      real(dp), dimension(:,:), allocatable :: X_ed ! Reordered intermediate, term 2
!
!     Batching variables
!
      integer :: rec0, rec1!, offset_eld
!
      integer :: current_d_batch
!
      type(batching_index) :: batch_d
!
!     :: Term 3. - sum_ckdlm b_ckal L_ilmd t_km^cd ::
!
      call mem%alloc(t_dmck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_dmck, (wf%n_o)*(wf%n_v))
!
      call mem%alloc(g_ilmd, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call wf%get_ooov(g_ilmd)
!
!     Form L_ildm = L_ilmd = 2 * g_ilmd - g_mlid
!                          = 2 * g_ilmd(i,l,m,d) - g_ilmd(m,l,i,d)
!
      call mem%alloc(L_ildm, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      L_ildm = zero
!
      call add_1243_to_1234(two, g_ilmd, L_ildm, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call add_4213_to_1234(-one, g_ilmd, L_ildm, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_ilmd, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
!     Form the intermediate X_ilck = sum_md L_ilmd t_mk^dc
!                                  = sum_md L_il_dm t_dm_ck
!
      call mem%alloc(X_ilck, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',            &
                  (wf%n_o)**2,       &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  L_ildm,            & ! L_il_dm
                  (wf%n_o)**2,       &
                  t_dmck,            & ! t_dm_ck
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_ilck,            & ! X_il_ck
                  (wf%n_o)**2)
!
      call mem%dealloc(L_ildm, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     Add - sum_ckdlm b_ckal L_ilmd t_km^cd
!         = - sum_ckl b_a_lck X_i_lck^T
!
      call dgemm('N','T',               &
                  wf%n_v,               &
                  wf%n_o,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  -one,                 &
                  b_aibj,               & ! b_a_lck (= b_al_ck = b_aibj)
                  wf%n_v,               &
                  X_ilck,               & ! X_i_lck
                  wf%n_o,               &
                  one,                  &
                  sigma_ai,             &
                  wf%n_v)
!
      call mem%dealloc(X_ilck, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     :: Term 4. - sum_ckdlm b_ckdl L_mlia t_km^cd ::
!
!     Form the intermediate X_ml = sum_ckd t_km^cd b_ckdl
!                                 = sum_ckd t_ckd_m^T b_ckd_l
!
      call mem%alloc(X_ml, wf%n_o, wf%n_o)
!
      call dgemm('T','N',               &
                  wf%n_o,               &
                  wf%n_o,               &
                  (wf%n_o)*(wf%n_v)**2, &
                  one,                  &
                  t_dmck,               & ! t_ckd_m
                  (wf%n_o)*(wf%n_v)**2, &
                  b_aibj,               & ! b_aib_j
                  (wf%n_o)*(wf%n_v)**2, &
                  zero,                 &
                  X_ml,                 &
                  wf%n_o)
!
      call mem%dealloc(t_dmck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Form g_mlia = g_mlia
!
      call mem%alloc(g_mlia, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call wf%get_ooov(g_mlia)
!
!     Form L_aiml = L_mlia = 2 * g_mlia - g_mail
!                          = 2 * g_mlia - g_ilma
!
      call mem%alloc(L_aiml, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      L_aiml = zero
!
      call add_3421_to_1234(two, g_mlia, L_aiml, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call add_2431_to_1234(-one, g_mlia, L_aiml, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call mem%dealloc(g_mlia, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
!     Add - sum_ckdlm b_ckdl L_mlia t_km^cd = - sum_lm L_ai_ml X_ml
!
      call dgemm('N','N',            &
                  (wf%n_v)*(wf%n_o), &
                  1,                 &
                  (wf%n_o)**2,       &
                  -one,              &
                  L_aiml,            & ! L_ai_ml
                  (wf%n_o)*(wf%n_v), &
                  X_ml,              & ! X_ml
                  (wf%n_o)**2,       &
                  one,               &
                  sigma_ai,          & ! sigma_ai
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(L_aiml, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%dealloc(X_ml, wf%n_o, wf%n_o)
!
!     :: Term 1. sum_ckdle b_ckdi L_dale t_kl^ce ::
!
!     Read amplitudes from disk
!
      call mem%alloc(t_elck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_elck, (wf%n_o)*(wf%n_v))
!
!     Form the intermediate X_eldi = sum_ck t_lk^ec b_ckdi
!                                  = sum_ck t_el_ck b_ck_di
!
      call mem%alloc(X_eldi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  t_elck,            & ! t_el_ck
                  (wf%n_o)*(wf%n_v), &
                  b_aibj,            & ! b_ck_di
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_eldi,            & ! X_el_di
                  (wf%n_o)*(wf%n_v))
!
!     sum_dle L_dale X_eldi
!
!     Prepare batching over index a
!
      rec0 = wf%n_v*wf%n_o*wf%integrals%n_J
!
      rec1 = wf%n_v*wf%integrals%n_J + wf%n_v**2*wf%n_o
!
      call batch_d%init(wf%n_v)
      call mem%batch_setup(batch_d, rec0, rec1)
!
      do current_d_batch = 1, batch_d%num_batches
!
         call batch_d%determine_limits(current_d_batch)
!
         call mem%alloc(g_dale, batch_d%length, wf%n_v, wf%n_o, wf%n_v)
!
            call wf%get_vvov(g_dale,      &
                           batch_d%first, &
                           batch_d%last,  &
                           1,             &
                           wf%n_v,        &
                           1,             &
                           wf%n_o,        &
                           1,             &
                           wf%n_v)
!
!        Form  L_aeld = L_dale = 2 * g_dale - g_dela
!                              = 2 * g_dale(d,a,l,e) - g_dale(d,e,l,a)
!
            call mem%alloc(L_aeld, wf%n_v, wf%n_v, wf%n_o, batch_d%length)
            L_aeld = zero
!
            call add_4132_to_1234(two, g_dale, L_aeld, wf%n_v, wf%n_v, wf%n_o, batch_d%length)
!
            call add_4231_to_1234(-one, g_dale, L_aeld, wf%n_v,  wf%n_v, wf%n_o, batch_d%length)
!
            call mem%dealloc(g_dale, batch_d%length, wf%n_v, wf%n_o, wf%n_v)
!
!           Add sum_ckdle b_ckdi L_dale t_kl^ce
!               = sum_eld L_a_eld X_eld_i
!
            call dgemm('N','N',                             &
                        wf%n_v,                             &
                        wf%n_o,                             &
                        (wf%n_o)*(wf%n_v)*batch_d%length,   &
                        one,                                &
                        L_aeld,                             & ! L_a_eld
                        wf%n_v,                             &
                        X_eldi(1, 1, batch_d%first, 1),     & ! X_eld_i
                        (wf%n_o)*(wf%n_v**2),               &
                        one,                                &
                        sigma_ai,                           &
                        wf%n_v)
!
            call mem%dealloc(L_aeld, wf%n_v, wf%n_v, wf%n_o, batch_d%length)
!
      enddo ! End of batches over a
!
      call mem%dealloc(X_eldi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     :: Term 2. sum_ckdle b_ckdl L_deia t_kl^ce ::
!
!     Form the intermediate X_de = sum_ckl b_ckdl t_kl^ce = sum_ckl b_d_lck t_e_lck^T
!
      call mem%alloc(X_de, wf%n_v, wf%n_v)
!
      call dgemm('N','T',               &
                  wf%n_v,               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  b_aibj,               & ! b_d_lck = b_dlck = b_ckdl
                  wf%n_v,               &
                  t_elck,               & ! t_e_lck
                  wf%n_v,               &
                  zero,                 &
                  X_de,                 &
                  wf%n_v)
!
      call mem%dealloc(t_elck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(X_ed, wf%n_v, wf%n_v)
!
      call sort_12_to_21(X_de, X_ed, wf%n_v, wf%n_v)
!
      call mem%dealloc(X_de, wf%n_v, wf%n_v)
!
!     sum_ckdle b_ckdl L_deia t_kl^ce = sum_de L_deia X_ed
!
!     Prepare batching over index d
!
      rec0 = wf%n_v*wf%n_o*wf%integrals%n_J
      rec1 = (wf%n_v**2)*(wf%n_o) + wf%n_v*wf%integrals%n_J
!
      call batch_d%init(wf%n_v)
      call mem%batch_setup(batch_d, rec0, rec1)
!
      do current_d_batch = 1, batch_d%num_batches
!
!        For each batch, get the limits for the d index
!
         call batch_d%determine_limits(current_d_batch)
!
!        Form g_deia
!
         call mem%alloc(g_deia, batch_d%length, wf%n_v, wf%n_o, wf%n_v)
!
         call wf%get_vvov(g_deia,         &
                           batch_d%first, &
                           batch_d%last,  &
                           1,             &
                           wf%n_v,        &
                           1,             &
                           wf%n_o,        &
                           1,             &
                           wf%n_v)
!
!        Form L_aied = L_deia = 2 * g_deia - g_daie
!                             = 2 * g_deia(d,e,i,a) - g_deia(d,a,i,e)
!
         call mem%alloc(L_aied, wf%n_v, wf%n_o, wf%n_v, batch_d%length)
         L_aied = zero
!
         call add_4321_to_1234(two, g_deia, L_aied, wf%n_v, wf%n_o, wf%n_v, batch_d%length)
         call add_4123_to_1234(-one, g_deia, L_aied, wf%n_v, wf%n_o, wf%n_v, batch_d%length)
!
         call mem%dealloc(g_deia, batch_d%length, wf%n_v, wf%n_o, wf%n_v)
!
         call dgemm('N','N',                    &
                     (wf%n_v)*(wf%n_o),         &
                     1,                         &
                     (wf%n_v)*(batch_d%length), &
                     one,                       &
                     L_aied,                    & ! L_ai_ed
                     (wf%n_v)*(wf%n_o),         &
                     X_ed(1,batch_d%first),     &
                     (wf%n_v)*(batch_d%length), &
                     one,                       &
                     sigma_ai,                  &
                     (wf%n_v)*(wf%n_o))
!
         call mem%dealloc(L_aied, wf%n_v, wf%n_o, wf%n_v, batch_d%length)
!
      enddo ! End of batches over d
!
      call mem%dealloc(X_ed, wf%n_v, wf%n_v)
!
   end subroutine jacobian_transpose_ccsd_e1_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_f1_ccsd(wf, sigma_ai, b_aibj)
!!
!!    Jacobian transpose CCSD F1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Calculates the F1 term,
!!
!!       sum_ckdlm (b_akdl t_lm^cd g_ikmc + b_ckal t_ml^cd g_mkid + b_ckdi t_ml^cd g_mkla)
!!
!!    and adds it to the transformed vector sigma_ai.
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o)                 :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: t_mcdl
      real(dp), dimension(:,:,:,:), allocatable :: t_cldm
      real(dp), dimension(:,:,:,:), allocatable :: t_cdml
!
      real(dp), dimension(:,:,:,:), allocatable :: b_akcl
      real(dp), dimension(:,:,:,:), allocatable :: b_kicd
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ikmc
      real(dp), dimension(:,:,:,:), allocatable :: g_kdmi
      real(dp), dimension(:,:,:,:), allocatable :: g_amkl
!
      real(dp), dimension(:,:,:,:), allocatable :: X_ikdl, X_kdli
!
      real(dp), dimension(:,:,:,:), allocatable :: X_akdm
!
      real(dp), dimension(:,:,:,:), allocatable :: X_kiml
      real(dp), dimension(:,:,:,:), allocatable :: X_mkli
!
!     :: Term 1. sum_ckdlm b_akdl t_lm^cd g_ikmc ::
!
!     X_ikdl = sum_mc t_lm^cd g_ikmc = sum_mc g_ik_mc t_mc_dl
!
!     Order amplitudes as t_mcdl = t_lm^cd
!
      call mem%alloc(t_mcdl, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      call squareup_and_sort_1234_to_4132(wf%t2, t_mcdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Form the integral g_ik_mc
!
      call mem%alloc(g_ikmc, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call wf%get_ooov(g_ikmc)
!
!     Form the intermediate X_ikdl = sum_mc t_lm^cd g_ikmc = sum_mc g_ik_mc t_mc_dl
!
      call mem%alloc(X_ikdl, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',            &
                  (wf%n_o)**2,       &
                  (wf%n_v)*(wf%n_o), &
                  (wf%n_v)*(wf%n_o), &
                  one,               &
                  g_ikmc,            & ! g_ik_mc
                  (wf%n_o)**2,       &
                  t_mcdl,            & ! t_mc_dl
                  (wf%n_v)*(wf%n_o), &
                  zero,              &
                  X_ikdl,            & ! X_ik_dl
                  (wf%n_o)**2)
!
!     Add sum_ckdlm b_akdl t_lm^cd g_ikmc
!         = sum_kdl b_a_kdl X_i_kdl^T
!
      call mem%alloc(X_kdli, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call sort_1234_to_2341(X_ikdl, X_kdli, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X_ikdl, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  wf%n_o,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  b_aibj,               & ! b_a_kdl
                  wf%n_v,               &
                  X_kdli,               & ! X_kdl_i
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  sigma_ai,             &
                  wf%n_v)
!
      call mem%dealloc(X_kdli, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
!     :: Term 2. sum_ckdlm b_ckal t_ml^cd g_mkid ::
!
!     X_akdm = sum_cl b_ckal t_ml^cd
!            = sum_cl b_ak_cl t_cl_dm
!
!     We have t_mcdl(m,c,d,l) = t_lm^cd
!     Reorder t_cldm(c,l,d,m) = t_mcdl(l,c,d,m) = t_ml^cd
!
      call mem%alloc(t_cldm, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call sort_1234_to_2134(t_mcdl, t_cldm, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      call mem%dealloc(t_mcdl, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
!     Reorder to b_akcl = b_ckal
!
      call mem%alloc(b_akcl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call sort_1234_to_3214(b_aibj, b_akcl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Form the intermediate X_akdm = sum_cl b_ak_cl t_cl_dm
!
      call mem%alloc(X_akdm, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',            &
                  (wf%n_v)*(wf%n_o), &
                  (wf%n_v)*(wf%n_o), &
                  (wf%n_v)*(wf%n_o), &
                  one,               &
                  b_akcl,            & ! b_ak_cl
                  (wf%n_v)*(wf%n_o), &
                  t_cldm,            & ! t_cl_dm
                  (wf%n_v)*(wf%n_o), &
                  zero,              &
                  X_akdm,            &
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(b_akcl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     sum_ckdlm b_ckal t_ml^cd g_mkid = sum_kdm X_akdm g_mkid
!
!     We have g_ikmc(i,k,m,c)
!     Reorder to g_kdmi(k,d,m,i) = g_mkid = g_ikmc(m,k,i,d)
!
      call mem%alloc(g_kdmi, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call sort_1234_to_2413(g_ikmc, g_kdmi, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call mem%dealloc(g_ikmc, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
!     Add sum_ckdlm b_ckal t_ml^cd g_mkid = sum_kdm X_akdm g_mkid
!                                         = sum_kdm X_akdm g_kdm_i
!
!     Note: we interpret X_akdm as X_a_kdm
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  wf%n_o,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  X_akdm,               & ! X_a_kdm
                  wf%n_v,               &
                  g_kdmi,               & ! g_kdm_i
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  sigma_ai,             &
                  wf%n_v)
!
      call mem%dealloc(X_akdm, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     :: Term 3. sum_ckdlm b_ckdi t_ml^cd g_mkla ::
!
!     X_kiml = sum_cd b_ckdi t_ml^cd
!
!     We have t_cldm(c,l,d,m) = t_ml^cd
!     Reorder to t_cdml(c,d,m,l) = t_cldm(c,l,d,m) = t_ml^cd
!
      call mem%alloc(t_cdml, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call sort_1234_to_1342(t_cldm, t_cdml, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(t_cldm, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Reorder to b_kicd = b_ckdi
!
      call mem%alloc(b_kicd, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
      call sort_1234_to_2413(b_aibj, b_kicd, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Form intermediate X_kiml = sum_cd b_ckdi t_ml^cd = sum_cd b_ki_cd t_cd_ml
!
      call mem%alloc(X_kiml, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call dgemm('N','N',            &
                  (wf%n_o)**2,       &
                  (wf%n_o)**2,       &
                  (wf%n_v)**2,       &
                  one,               &
                  b_kicd,            & ! b_ki_cd
                  (wf%n_o)**2,       &
                  t_cdml,            & ! t_cd_ml
                  (wf%n_v)**2,       &
                  zero,              &
                  X_kiml,            & ! X_ki_ml
                  (wf%n_o)**2)
!
      call mem%dealloc(t_cdml, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(b_kicd, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
!     sum_ckdlm b_ckdi t_ml^cd g_mkla = sum_klm g_mkla X_ki_ml
!
!     Reorder to X_mkli
!
      call mem%alloc(X_mkli, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call sort_1234_to_3142(X_kiml, X_mkli, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call mem%dealloc(X_kiml, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
!     We have g_kdm_i(kdm,i) = g_mkid
!     Reorder to g_amkl(a,m,k,l) = g_mkla = g_kdm_i(kam,l)
!
      call mem%alloc(g_amkl, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call sort_1234_to_2314(g_kdmi, g_amkl, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call mem%dealloc(g_kdmi, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
!     Add sum_ckdlm b_ckdi t_ml^cd g_mkla = sum_klm g_a_mkl X_mkl_i
!
      call dgemm('N','N',      &
                  wf%n_v,      &
                  wf%n_o,      &
                  (wf%n_o)**3, &
                  one,         &
                  g_amkl,      & ! g_a_mkl
                  wf%n_v,      &
                  X_mkli,      & ! X_mkl_i
                  (wf%n_o)**3, &
                  one,         &
                  sigma_ai,    &
                  wf%n_v)
!
      call mem%dealloc(g_amkl, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%dealloc(X_mkli, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
   end subroutine jacobian_transpose_ccsd_f1_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_g1_ccsd(wf, sigma_ai, b_aibj)
!!
!!    Jacobian transpose CCSD G1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Calculates the G1 term,
!!
!!       - sum_ckdle (b_akdl t_kl^ce g_icde + b_cidl t_kl^ce g_kade + b_cldi t_kl^ce g_keda)
!!
!!    and adds it to the transformed vector sigma_ai.
!!
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o)                 :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: b_dicl ! b_cidl
!
      real(dp), dimension(:,:,:,:), allocatable :: X_diek ! An intermediate, term 2
      real(dp), dimension(:,:,:,:), allocatable :: X_kdei ! Reordered intermediate, term 2
!
      real(dp), dimension(:,:,:,:), allocatable :: X_kedi ! Reordered intermediate, term 3
!
      real(dp), dimension(:,:,:,:), allocatable :: X_idkl ! An intermediate, term 1
      real(dp), dimension(:,:,:,:), allocatable :: X_kdli ! Reordered intermediate, term 1
!
      real(dp), dimension(:,:,:,:), allocatable :: t_clek ! t_kl^ce
      real(dp), dimension(:,:,:,:), allocatable :: t_cekl ! t_kl^ce
!
      real(dp), dimension(:,:,:,:), allocatable :: g_kade
      real(dp), dimension(:,:,:,:), allocatable :: g_akde
      real(dp), dimension(:,:,:,:), allocatable :: g_keda
      real(dp), dimension(:,:,:,:), allocatable :: g_icde
      real(dp), dimension(:,:,:,:), allocatable :: g_idce ! g_icde
!
!     Batching variables
!
      integer :: current_a_batch = 0
      integer :: current_d_batch = 0
      integer :: current_e_batch = 0
!
      type(batching_index) :: batch_a
      type(batching_index) :: batch_d
      type(batching_index) :: batch_e
!
      integer :: rec1, rec0
!
!     :: Term 2. - sum_ckdle b_cidl t_kl^ce g_kade ::
!
!     X_diek = sum_cl b_cidl t_kl^ce = sum_cl b_di_cl t_cl_ek
!
!     Reorder to b_dicl = b_cidl
!
      call mem%alloc(b_dicl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_3214(b_aibj, b_dicl, wf%n_v, wf%n_o, wf%n_v, wf%n_o) ! b_aibj = b_cidl
!
!     Order amplitudes as t_clek = t_kl^ce
!
      call mem%alloc(t_clek, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup_and_sort_1234_to_1432(wf%t2, t_clek, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Form the intermediate X_diek = sum_cl b_di_cl t_cl_ek
!
      call mem%alloc(X_diek, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  b_dicl,            & ! b_di_cl
                  (wf%n_o)*(wf%n_v), &
                  t_clek,            & ! t_cl_ek
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_diek,            &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(b_dicl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     - sum_ckdle b_cidl t_kl^ce g_kade = sum_kde g_kade X_diek
!
!     Reorder X_diek to X_kdei
!
      call mem%alloc(X_kdei, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call sort_1234_to_4132(X_diek, X_kdei, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(X_diek, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Prepare batching over index e
!
      rec0 = wf%n_v*wf%n_o*wf%integrals%n_J
      rec1 = 2*(wf%n_v**2)*(wf%n_o) + wf%n_v*wf%integrals%n_J
!
!     Initialize batching variable
!
      call batch_e%init(wf%n_v)
      call mem%batch_setup(batch_e, rec0, rec1)
!
!     Loop over the e-batches
!
      do current_e_batch = 1, batch_e%num_batches
!
!        For each batch, get the limits for the e index
!
         call batch_e%determine_limits(current_e_batch)
!
!        Form g_kade
!
         call mem%alloc(g_kade, wf%n_o, wf%n_v, wf%n_v, batch_e%length)
!
         call wf%get_ovvv(g_kade,         &
                           1,             &
                           wf%n_o,        &
                           1,             &
                           wf%n_v,        &
                           1,             &
                           wf%n_v,        &
                           batch_e%first, &
                           batch_e%last)
!
!        Reorder to g_akde = g_kade
!
         call mem%alloc(g_akde, wf%n_v, wf%n_o, wf%n_v, batch_e%length)
!
         call sort_1234_to_2134(g_kade, g_akde, wf%n_o, wf%n_v, wf%n_v, batch_e%length)
!
         call mem%dealloc(g_kade, wf%n_o, wf%n_v, wf%n_v, batch_e%length)
!
         call dgemm('N','N',                             &
                     wf%n_v,                             &
                     wf%n_o,                             &
                     (wf%n_o)*(wf%n_v)*(batch_e%length), &
                     -one,                               &
                     g_akde,                             &
                     wf%n_v,                             &
                     X_kdei(1, 1, batch_e%first,1),      &
                     (wf%n_o)*(wf%n_v)**2,               &
                     one,                                &
                     sigma_ai,                           &
                     wf%n_v)
!
         call mem%dealloc(g_akde, wf%n_v, wf%n_o, wf%n_v, batch_e%length)
!
      enddo ! End of batches over e
!
      call mem%dealloc(X_kdei, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
!     :: Term 3. - sum_ckdle b_cldi t_kl^ce g_keda ::
!
!     X_diek = sum_cl b_cldi t_kl^ce = sum_cl b_di_cl t_cl_ek
!
      call mem%alloc(X_diek, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T','N',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  b_aibj,            & ! b_cl_di
                  (wf%n_o)*(wf%n_v), &
                  t_clek,            & ! t_cl_ek
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_diek,            &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(t_clek, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     - sum_kde X_diek g_keda
!
!     Reorder X_diek to X_kedi
!
      call mem%alloc(X_kedi, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call sort_1234_to_4312(X_diek, X_kedi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(X_diek, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Prepare batching over a
!
      rec0 = wf%n_v*wf%n_o*wf%integrals%n_J
      rec1 = (wf%n_v**2)*(wf%n_o) + wf%n_v*wf%integrals%n_J
!
!     Initialize batching variable
!
      call batch_a%init(wf%n_v)
      call mem%batch_setup(batch_a, rec0, rec1)
!
!     Loop over the a-batches
!
      do current_a_batch = 1, batch_a%num_batches
!
!        For each batch, get the limits for the a index
!
         call batch_a%determine_limits(current_a_batch)
!
!        Form g_keda
!
         call mem%alloc(g_keda, wf%n_o, wf%n_v, wf%n_v, batch_a%length)
!
         call wf%get_ovvv(g_keda,         &
                           1,             &
                           wf%n_o,        &
                           1,             &
                           wf%n_v,        &
                           1,             &
                           wf%n_v,        &
                           batch_a%first, &
                           batch_a%last)
!
         call dgemm('T','N',                    &
                     batch_a%length,            &
                     wf%n_o,                    &
                     (wf%n_o)*(wf%n_v)**2,      &
                     -one,                      &
                     g_keda,                    & ! g_ked_a
                     (wf%n_o)*(wf%n_v)**2,      &
                     X_kedi,                    & ! X_ked_i
                     (wf%n_o)*(wf%n_v)**2,      &
                     one,                       &
                     sigma_ai(batch_a%first,1), &
                     wf%n_v)
!
         call mem%dealloc(g_keda, wf%n_o, wf%n_v, wf%n_v, batch_a%length)
!
      enddo ! End of batches over a
!
      call mem%dealloc(X_kedi, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
!     :: Term 1. - sum_ckdle b_akdl t_kl^ce g_icde ::
!
!     X_idkl = sum_ce t_kl^ce g_icde = sum_ce g_id_ce t_ce_kl
!
!     Reorder to t_cekl = t_kl^ce
!
      call mem%alloc(t_cekl, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call squareup_and_sort_1234_to_1324(wf%t2, t_cekl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(X_idkl, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      X_idkl = zero
!
!     Prepare for batching over d
!
      rec0 = wf%n_v*wf%n_o*wf%integrals%n_J
      rec1 = 2*(wf%n_v**2)*(wf%n_o) + wf%n_v*wf%integrals%n_J
!
      call batch_d%init(wf%n_v)
      call mem%batch_setup(batch_d, rec0, rec1)
!
!     Loop over the d-batches
!
      do current_d_batch = 1, batch_d%num_batches
!
         call batch_d%determine_limits(current_d_batch)
!
         call mem%alloc(g_icde, wf%n_o, wf%n_v, batch_d%length, wf%n_v)
!
         call wf%get_ovvv(g_icde,         &
                           1,             &
                           wf%n_o,        &
                           1,             &
                           wf%n_v,        &
                           batch_d%first, &
                           batch_d%last,  &
                           1,             &
                           wf%n_v)
!
!        Reorder to g_id_ce = g_ic_de
!
         call mem%alloc(g_idce, wf%n_o, batch_d%length, wf%n_v, wf%n_v)
!
         call sort_1234_to_1324(g_icde, g_idce, wf%n_o, wf%n_v, (batch_d%length), wf%n_v)
!
         call mem%dealloc(g_icde, wf%n_o, wf%n_v, batch_d%length, wf%n_v)
!
         call dgemm('N','N',                      &
                     (wf%n_o)*(batch_d%length),   &
                     (wf%n_o)**2,                 &
                     (wf%n_v)**2,                 &
                     one,                         &
                     g_idce,                      &
                     (wf%n_o)*(batch_d%length),   &
                     t_cekl,                      &
                     (wf%n_v)**2,                 &
                     one,                         &
                     X_idkl(1,batch_d%first,1,1), &
                     (wf%n_o)*(wf%n_v))
!
         call mem%dealloc(g_idce, wf%n_o, batch_d%length, wf%n_v, wf%n_v)
!
      enddo ! End of batches over d
!
      call mem%dealloc(t_cekl, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
!     - sum_ckdle b_akdl t_kl^ce g_icde = sum_kdl b_akdl X_id_kl
!
      call mem%alloc(X_kdli, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call sort_1234_to_3241(X_idkl, X_kdli, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call mem%dealloc(X_idkl, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
!     Add - sum_ckdle b_akdl t_kl^ce g_icde = - sum_dkl b_a_kdl X_kdl_i
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  wf%n_o,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  -one,                 &
                  b_aibj,               & ! b_aibj
                  wf%n_v,               &
                  X_kdli,               & ! X_kdl_i
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  sigma_ai,             &
                  wf%n_v)
!
      call mem%dealloc(X_kdli, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
   end subroutine jacobian_transpose_ccsd_g1_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_a2_ccsd(wf, sigma_aibj, b_ai)
!!
!!    Jacobian transpose CCSD A2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Calculates the A2 term,
!!
!!       2 F_jb b_ai - F_ib b_aj - sum_k L_ikjb b_ak + sum_c L_cajb b_ci
!!
!!    and adds it to the transformed vector sigma_aibj.
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: sigma_aibj
      real(dp), dimension(wf%n_v, wf%n_o)                 :: b_ai
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ikjb
      real(dp), dimension(:,:,:,:), allocatable :: L_kibj ! L_ikjb
!
      real(dp), dimension(:,:,:,:), allocatable :: g_cajb
      real(dp), dimension(:,:,:,:), allocatable :: g_cbja
!
      real(dp), dimension(:,:,:,:), allocatable :: sigma_iajb ! sigma_aibj contribution
      real(dp), dimension(:,:,:,:), allocatable :: sigma_ibja ! sigma_aibj contribution
!
      integer :: i, a, j, b
!
!     Batching variables
!
      integer :: rec0, rec1
!
      integer :: current_a_batch = 0
      integer :: current_b_batch = 0
!
      type(batching_index) :: batch_a
      type(batching_index) :: batch_b
!
!     :: Term 1 & 2. 2 F_jb b_ai - F_ib b_aj ::
!
!$omp parallel do private(i, a, j, b)
      do j = 1, wf%n_o
         do b = 1, wf%n_v
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  sigma_aibj(a,i,b,j) = -(wf%fock_ia(i, b))*b_ai(a, j) + two*(wf%fock_ia(j, b))*b_ai(a, i)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do

!     :: Term 3. - sum_k L_ikjb b_ak ::
!
!     Form g_ikjb
!
!
      call mem%alloc(g_ikjb, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call wf%get_ooov(g_ikjb)
!
!     Form L_kibj = L_ikjb = 2 * g_ikjb - g_jkib
!                          = 2 * g_ikjb(i,k,j,b) - g_ikjb(j,k,i,b)
!
      call mem%alloc(L_kibj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      L_kibj = zero
!
      call add_2143_to_1234(two, g_ikjb, L_kibj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call add_4123_to_1234(-one, g_ikjb, L_kibj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_ikjb, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
!     Add - sum_k L_ikjb b_ak = - sum_k b_a_k L_k_ibj
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_o,               &
                  -one,                 &
                  b_ai,                 & ! "b_a_k"
                  wf%n_v,               &
                  L_kibj,               & ! L_k_ibj
                  wf%n_o,               &
                  one,                  &
                  sigma_aibj,           & ! "sigma_aibj"
                  wf%n_v)
!
      call mem%dealloc(L_kibj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     :: Term 4. 2 sum_c g_cajb b_ci - sum_c g_cbja b_ci ::
!
!     2 sum_c g_cajb
!
!     Prepare for batching over a
!
      rec0 = wf%n_v*wf%n_o*wf%integrals%n_J
      rec1 = (wf%n_v**2)*(wf%n_o) + (wf%n_v)*(wf%n_o**2) + wf%n_v*wf%integrals%n_J
!
!     Initialize batching variable
!
      call batch_a%init(wf%n_v)
      call mem%batch_setup(batch_a, rec0, rec1)
!
!     Loop over the a-batches
!
      do current_a_batch = 1, batch_a%num_batches
!
!        For each batch, get the limits for the a index
!
         call batch_a%determine_limits(current_a_batch)
!
!        Form g_cajb
!
         call mem%alloc(g_cajb, wf%n_v, batch_a%length, wf%n_o, wf%n_v)
!
         call wf%get_vvov(g_cajb,                        &
                           1, wf%n_v,                    &
                           batch_a%first, batch_a%last,  &
                           1, wf%n_o,                    &
                           1, wf%n_v)
!
!        Add 2 sum_c g_cajb b_ci = 2 sum_c b_c_i^T(i,c) g_c_ajb
!
         call mem%alloc(sigma_iajb, wf%n_o, batch_a%length, wf%n_o, wf%n_v)
!
         call dgemm('T','N',                             &
                     wf%n_o,                             &
                     (wf%n_o)*(wf%n_v)*(batch_a%length), &
                     wf%n_v,                             &
                     two,                                &
                     b_ai,                               & ! b_c_i
                     wf%n_v,                             &
                     g_cajb,                             & ! g_c_ajb
                     wf%n_v,                             &
                     zero,                               &
                     sigma_iajb,                         &
                     wf%n_o)
!
         call mem%dealloc(g_cajb, wf%n_v, batch_a%length, wf%n_o, wf%n_v)
!
         call add_2143_to_1234(one, sigma_iajb, sigma_aibj, batch_a%length, wf%n_o, wf%n_v, wf%n_o)
!
         call mem%dealloc(sigma_iajb, wf%n_o, batch_a%length, wf%n_o, wf%n_v)
!
      enddo ! End of batches over a
!
!     - sum_c g_cbja b_ci
!
!     Prepare for batching over b
!
      rec0 = wf%n_v*wf%n_o*wf%integrals%n_J
      rec1 = (wf%n_v**2)*(wf%n_o) + (wf%n_v)*(wf%n_o**2) + wf%n_v*wf%integrals%n_J
!
!     Initialize batching variable
!
      call batch_b%init(wf%n_v)
      call mem%batch_setup(batch_b, rec0, rec1)
!
!     Loop over the number of b batches
!
      do current_b_batch = 1, batch_b%num_batches
!
!        For each batch, get the limits for the b index
!
         call batch_b%determine_limits(current_b_batch)
!
!        Form g_cbja
!
         call mem%alloc(g_cbja, wf%n_v, batch_b%length, wf%n_o, wf%n_v)
!
         call wf%get_vvov(g_cbja,                        &
                           1, wf%n_v,                    &
                           batch_b%first, batch_b%last,  &
                           1, wf%n_o,                    &
                           1, wf%n_v)
!
!         Form - sum_c g_cbja b_ci = - sum_c b_ci^T(i,c) g_c_bja
!
         call mem%alloc(sigma_ibja, wf%n_o, batch_b%length, wf%n_o, wf%n_v)
!
         call dgemm('T','N',                             &
                     wf%n_o,                             &
                     (wf%n_o)*(wf%n_v)*(batch_b%length), &
                     wf%n_v,                             &
                     -one,                               &
                     b_ai,                               & ! b_c_i
                     wf%n_v,                             &
                     g_cbja,                             & ! g_c_bja
                     wf%n_v,                             &
                     zero,                               &
                     sigma_ibja,                         &
                     wf%n_o)
!
         call mem%dealloc(g_cbja, wf%n_v, batch_b%length, wf%n_o, wf%n_v)
!
!        Add it to sigma_aibj
!
         call add_2341_to_1234(one, sigma_ibja, sigma_aibj, wf%n_v, wf%n_o, batch_b%length, wf%n_o)
!
         call mem%dealloc(sigma_ibja, wf%n_o, batch_b%length, wf%n_o, wf%n_v)
!
      enddo ! End of batches over b
!
   end subroutine jacobian_transpose_ccsd_a2_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_b2_ccsd(wf, sigma_aibj, b_aibj)
!!
!!    Jacobian transpose CCSD B2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Calculates the B2 term,
!!
!!       sum_c b_aicj F_cb - sum_k b_aibk F_jk + sum_ck b_aick L_ckjb
!!
!!    and adds it to the transformed vector sigma_aibj.
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: sigma_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: b_aijc ! b_aicj
!
      real(dp), dimension(:,:,:,:), allocatable :: sigma_aijb ! sigma_aibj contribution
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ckjb
      real(dp), dimension(:,:,:,:), allocatable :: g_ckbj ! g_ckjb & g_cbjk
!
      real(dp), dimension(:,:,:,:), allocatable :: g_cbjk_restricted ! g_cbjk, batch over b
!
      integer :: k, j, b, c
!
!     Batching variables
!
      integer :: rec0, rec1
      integer :: current_b_batch = 0
!
      type(batching_index) :: batch_b
!
!     :: Term 1. sum_c b_aicj F_cb ::
!
!     Reorder to b_aijc = b_aicj
!
      call mem%alloc(b_aijc, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      b_aijc = zero
!
      call sort_1234_to_1243(b_aibj, b_aijc, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Calculate and add sum_c b_aicj F_cb = sum_c b_aij_c F_c_b
!
      call mem%alloc(sigma_aijb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N','N',               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_v,               &
                  wf%n_v,               &
                  one,                  &
                  b_aijc,               & ! b_aij_c
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%fock_ab,           & ! F_c_b
                  wf%n_v,               &
                  zero,                 &
                  sigma_aijb,           &
                  (wf%n_v)*(wf%n_o)**2)
!
      call mem%dealloc(b_aijc, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call add_1243_to_1234(one, sigma_aijb, sigma_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(sigma_aijb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
!     :: Term 2. - sum_k b_aibk F_jk ::
!
!     - sum_k b_aibk F_jk = - sum_k b_aib_k F_jk^T(k,j)
!
      call dgemm('N','T',               &
                  (wf%n_o)*(wf%n_v)**2, &
                  wf%n_o,               &
                  wf%n_o,               &
                  -one,                 &
                  b_aibj,               & ! b_aib_k
                  (wf%n_o)*(wf%n_v)**2, &
                  wf%fock_ij,           & ! F_j_k
                  wf%n_o,               &
                  one,                  &
                  sigma_aibj,           & ! sigma_aib_j
                  (wf%n_o)*(wf%n_v)**2)
!
!     :: Term 3. sum_ck b_aick L_ckjb ::
!
!     Form g_ckjb
!
      call mem%alloc(g_ckjb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call wf%get_voov(g_ckjb)
!
!     Reorder to g_ckbj = g_ckjb
!
      call mem%alloc(g_ckbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call sort_1234_to_1243(g_ckjb, g_ckbj, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call mem%dealloc(g_ckjb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
!     Add 2 * sum_ck b_aick g_ckjb = 2 * sum_ck b_aick g_ckbj
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  two,               &
                  b_aibj,            & ! "b_ai_ck"
                  (wf%n_o)*(wf%n_v), &
                  g_ckbj,            & ! g_ck_bj
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  sigma_aibj,        &
                  (wf%n_o)*(wf%n_v))
!
!     - sum_ck b_aick g_cbjk
!
!     Prepare to batch over b to make g_cb_jk = g_cbjk successively
!
      g_ckbj = zero ! g_cbjk reordered
!
      rec0 = wf%n_o**2*wf%integrals%n_J
      rec1 = wf%n_v*wf%integrals%n_J  + (wf%n_o**2)*(wf%n_v)
!
!     Initialize batching variable
!
      call batch_b%init(wf%n_v)
      call mem%batch_setup(batch_b, rec0, rec1)
!
!     Loop over the number of b batches
!
      do current_b_batch = 1, batch_b%num_batches
!
!        For each batch, get the limits for the b index
!
         call batch_b%determine_limits(current_b_batch)
!
!        Form g_cb_jk = g_cbjk
!
         call mem%alloc(g_cbjk_restricted, wf%n_v, batch_b%length, wf%n_o, wf%n_o)
!
         call wf%get_vvoo(g_cbjk_restricted, &
                           1,                &
                           wf%n_v,           &
                           batch_b%first,    &
                           batch_b%last,     &
                           1,                &
                           wf%n_o,           &
                           1,                &
                           wf%n_o)
!
!        Place in reordered full space vector and deallocate restricted vector
!
!$omp parallel do schedule(static) private(k,j,b,c)
         do k = 1, wf%n_o
            do j = 1, wf%n_o
               do b = batch_b%first, batch_b%last
                  do c = 1, wf%n_v
!
                     g_ckbj(c,k,b,j) = g_cbjk_restricted(c,b-batch_b%first+1,j,k)
!
                  enddo
               enddo
            enddo
         enddo
!$omp end parallel do
!
         call mem%dealloc(g_cbjk_restricted, wf%n_v, batch_b%length, wf%n_o, wf%n_o)
!
      enddo ! End of batches over b
!
!     Add  - sum_ck b_aick g_cbjk = - sum_ck b_ai_ck g_ck_bj
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  -one,              &
                  b_aibj,            & ! "b_ai_ck"
                  (wf%n_o)*(wf%n_v), &
                  g_ckbj,            & ! g_ck_bj
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  sigma_aibj,        &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(g_ckbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine jacobian_transpose_ccsd_b2_ccsd
!

!
   module subroutine jacobian_transpose_ccsd_c2_ccsd(wf, sigma_aibj, b_aibj)
!!
!!    Jacobian transpose CCSD C2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Calculates the C2 term,
!!
!!       - sum_ck (b_ajck g_ibck + b_akcj g_ikcb)
!!
!!    and adds it to the transformed vector sigma_aibj.
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: sigma_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ibck
      real(dp), dimension(:,:,:,:), allocatable :: g_cbik
      real(dp), dimension(:,:,:,:), allocatable :: g_ckbi ! g_cbik
!
      real(dp), dimension(:,:,:,:), allocatable :: sigma_ajib ! sigma_aibj contribution
      real(dp), dimension(:,:,:,:), allocatable :: sigma_ajbi ! sigma_aibj contribution
!
      real(dp), dimension(:,:,:,:), allocatable :: b_ajck ! b_akcj
!
      integer :: k, i, b, c
!
!     Batching variables
!
      integer :: rec1, rec0
      integer :: current_b_batch = 0
!
      type(batching_index) :: batch_b
!
!     :: Term 1. - sum_ck b_ajck g_ibck ::
!
!     Form g_ibck
!
      call mem%alloc(g_ibck, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      call wf%get_ovvo(g_ibck)
!
!     Calculate and add - sum_ck b_ajck g_ibck = - sum_ck b_aj_ck g_ib_ck^T(c,k,i,b)
!
      call mem%alloc(sigma_ajib, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N','T',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  -one,              &
                  b_aibj,            & ! b_aj_ck
                  (wf%n_o)*(wf%n_v), &
                  g_ibck,            & ! g_ib_ck
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  sigma_ajib,        &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(g_ibck, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      call add_1423_to_1234(one, sigma_ajib, sigma_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(sigma_ajib, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
!     :: Term 2. - sum_ck b_akcj g_ikcb = - sum_ck b_akcj g_cbik ::
!
!     Make g_ckbi = g_cbik in batches over b
!
      call mem%alloc(g_ckbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      g_ckbi = zero
!
      rec0 = wf%n_o**2*wf%integrals%n_J
      rec1 = wf%n_v*wf%integrals%n_J  + (wf%n_o**2)*(wf%n_v)
!
!     Initialize batching variable
!
      call batch_b%init(wf%n_v)
      call mem%batch_setup(batch_b, rec0, rec1)
!
!     Loop over the b-batches
!
      do current_b_batch = 1, batch_b%num_batches
!
!        For each batch, get the limits for the b index
!
         call batch_b%determine_limits(current_b_batch)
!
!        Form g_cbik
!
         call mem%alloc(g_cbik, wf%n_v, batch_b%length, wf%n_o, wf%n_o)
!
         call wf%get_vvoo(g_cbik,         &
                           1,             &
                           wf%n_v,        &
                           batch_b%first, &
                           batch_b%last,  &
                           1,             &
                           wf%n_o,        &
                           1,             &
                           wf%n_o)
!
!        Place in reordered integral g_ckbi = g_cbik
!
!$omp parallel do schedule(static) private(k,i,b,c)
         do k = 1, wf%n_o
            do i = 1, wf%n_o
               do b = batch_b%first, batch_b%last
                  do c = 1, wf%n_v
!
                     g_ckbi(c,k,b,i) = g_cbik(c,b-batch_b%first+1,i,k)
!
                  enddo
               enddo
            enddo
         enddo
!$omp end parallel do
!
         call mem%dealloc(g_cbik, wf%n_v, batch_b%length, wf%n_o, wf%n_o)
!
      enddo ! End of batches over b
!
!     Reorder to b_ajck = b_akcj
!
      call mem%alloc(b_ajck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call sort_1234_to_1432(b_aibj, b_ajck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Form and add - sum_ck b_akcj g_cbik = - sum_ck b_ajck g_ck_bi
!
      call mem%alloc(sigma_ajbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  -one,              &
                  b_ajck,            & ! b_aj_ck
                  (wf%n_o)*(wf%n_v), &
                  g_ckbi,            & ! g_ck_bi
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  sigma_ajbi,        &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(b_ajck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(g_ckbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call add_1432_to_1234(one, sigma_ajbi, sigma_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(sigma_ajbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine jacobian_transpose_ccsd_c2_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_d2_ccsd(wf, sigma_aibj, b_aibj)
!!
!!    Jacobian transpose CCSD D2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Calculates the D2 term,
!!
!!       2 * sum_ckdl b_aick L_jbld t_kl^cd
!!
!!    and adds it to the transformed vector sigma_aibj.
!!
   implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: sigma_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: t_ckdl ! t_kl^cd
!
      real(dp), dimension(:,:,:,:), allocatable :: X_ckbj ! An intermediate
!
      real(dp), dimension(:,:,:,:), allocatable :: g_jbld ! g_jbld
      real(dp), dimension(:,:,:,:), allocatable :: L_dlbj ! L_jbld
!
!     Form g_jbld
!
      call mem%alloc(g_jbld, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_ovov(g_jbld)
!
!     Form L_dlbj = L_jbld = 2 * g_jbld - g_jdlb
!                          = 2 * g_jbld(j,b,l,d) - g_jbld(j,d,l,b)
!
      call mem%alloc(L_dlbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      L_dlbj = zero
!
      call add_4321_to_1234(two, g_jbld, L_dlbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call add_4123_to_1234(-one, g_jbld, L_dlbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_jbld, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     Form t_ckdl = t_kl^cd
!
      call mem%alloc(t_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      t_ckdl = zero
!
      call squareup(wf%t2, t_ckdl, (wf%n_o)*(wf%n_v))
!
!     Form the intermediate X_ckbj = sum_dl t_ck_dl L_dl_bj
!
      call mem%alloc(X_ckbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  t_ckdl,            & ! t_ck_dl
                  (wf%n_o)*(wf%n_v), &
                  L_dlbj,            & ! L_dl_bj
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_ckbj,            &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(t_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(L_dlbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Add 2 * sum_ckdl b_aick L_jbld t_kl^cd = 2 * sum_ck b_ai_ck X_ck_bj
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  two,               &
                  b_aibj,            & ! "b_ai_ck"
                  (wf%n_o)*(wf%n_v), &
                  X_ckbj,            & ! X_ck_bj
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  sigma_aibj,        &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(X_ckbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine jacobian_transpose_ccsd_d2_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_e2_ccsd(wf, sigma_aibj, b_aibj)
!!
!!    Jacobian transpose CCSD E2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Calculates the E2 term,
!!
!!       - sum_ckdl (b_aibl t_kl^cd L_kcjd + b_aicl t_kl^cd L_jbkd + b_aicj t_kl^cd L_ldkb)
!!
!!    and adds it to the transformed vector sigma_aibj.
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: sigma_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: g_kcjd
      real(dp), dimension(:,:,:,:), allocatable :: L_jckd ! L_kcjd
      real(dp), dimension(:,:,:,:), allocatable :: L_dkbj ! L_jbkd
      real(dp), dimension(:,:,:,:), allocatable :: L_ldkb
!
      real(dp), dimension(:,:,:,:), allocatable :: t_ckdl ! t_kl^cd
      real(dp), dimension(:,:,:,:), allocatable :: t_cldk ! t_kl^cd
!
      real(dp), dimension(:,:), allocatable :: X_jl   ! An intermediate, term 1
      real(dp), dimension(:,:,:,:), allocatable :: X_clbj ! An intermediate, term 2
      real(dp), dimension(:,:), allocatable :: X_cb   ! An intermediate, term 3
!
      real(dp), dimension(:,:,:,:), allocatable :: sigma_aijb ! sigma_aibj contribution
      real(dp), dimension(:,:,:,:), allocatable :: b_aijc     ! b_aicj
!
!     :: Term 1. - sum_ckdl b_aibl t_kl^cd L_kcjd ::
!
!     X_jl = sum_kcd L_kcjd t_kl^cd = sum_kcd L_j_ckd t_ckd_l
!
!     Form g_kcjd
!
      call mem%alloc(g_kcjd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_ovov(g_kcjd)
!
!     Form L_jckd = L_kcjd = 2 * g_kcjd - g_kdjc
!                          = 2 * g_kcjd(k,c,j,d) - g_kcjd(k,d,j,c)
!
      call mem%alloc(L_jckd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
         L_jckd = zero
!
      call add_3214_to_1234(two, g_kcjd, L_jckd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call add_3412_to_1234(-one, g_kcjd, L_jckd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%dealloc(g_kcjd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     Form t_ckdl = t_kl^cd
!
      call mem%alloc(t_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call squareup(wf%t2, t_ckdl, (wf%n_o)*(wf%n_v))
!
!     Form the intermediate X_jl = sum_kcd L_kcjd t_kl^cd
!                                 = sum_kcd L_j_ckd t_ckd_l
!
      call mem%alloc(X_jl, wf%n_o, wf%n_o)
!
      call dgemm('N','N',               &
                  wf%n_o,               &
                  wf%n_o,               &
                  (wf%n_o)*(wf%n_v)**2, &
                  one,                  &
                  L_jckd,               & ! L_j_ckd
                  wf%n_o,               &
                  t_ckdl,               & ! "t_ckd_l"
                  (wf%n_o)*(wf%n_v)**2, &
                  zero,                 &
                  X_jl,                 &
                  wf%n_o)
!
!     Add - sum_ckdl b_aibl t_kl^cd L_kcjd
!         = - sum_l b_aib_l X_jl^T(l,j)
!
      call dgemm('N','T',               &
                  (wf%n_o)*(wf%n_v)**2, &
                  wf%n_o,               &
                  wf%n_o,               &
                  -one,                 &
                  b_aibj,               & ! "b_aib_l"
                  (wf%n_o)*(wf%n_v)**2, &
                  X_jl,                 &
                  wf%n_o,               &
                  one,                  &
                  sigma_aibj,           & ! "sigma_aib_j"
                  (wf%n_o)*(wf%n_v)**2)
!
      call mem%dealloc(X_jl, wf%n_o, wf%n_o)
!
!     :: Term 2. -sum_ckdl b_aicl t_kl^cd L_jbkd ::
!
!     X_clbj = sum_kd t_kl^cd L_jbkd = sum_kd t_cl_dk L_dk_bj
!
!     We have L_jckd = L_kcjd   => L_jckd(k,b,j,d) = L_jbkd
!     We have t_ckd_l = t_kl^cd
!
!     Reorder to L_dkbj(d,k,b,j) = L_jbkd = L_jckd(k,b,j,d)
!
      call mem%alloc(L_dkbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call sort_1234_to_4123(L_jckd, L_dkbj, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%dealloc(L_jckd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     Reorder to t_cldk = t_kl^cd = t_ckdl
!
      call mem%alloc(t_cldk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call sort_1234_to_1432(t_ckdl, t_cldk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(t_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Form the intermediate X_clbj = sum_kd t_kl^cd L_jbkd
!                                  = sum_kd t_cl_dk L_dk_bj
!
      call mem%alloc(X_clbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  t_cldk,            & ! t_cl_dk
                  (wf%n_o)*(wf%n_v), &
                  L_dkbj,            & ! L_dk_bj
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_clbj,            &
                  (wf%n_o)*(wf%n_v))
!
!     Add - sum_ckdl b_aicl t_kl^cd L_jbkd
!           = - sum_cl b_ai_cl X_cl_bj
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  -one,              &
                  b_aibj,            & ! "b_ai_cl"
                  (wf%n_o)*(wf%n_v), &
                  X_clbj,            & ! X_cl_bj
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  sigma_aibj,        &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(X_clbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     :: Term 3. - sum_ckdl b_aicj t_kl^cd L_ldkb ::
!
!     - sum_c b_aijc X_cb,   X_cb = sum_kdl t_kl^cd L_ldkb
!                                    = sum_kdl t_c_ldk L_ldk_b
!
!     We have L_dkbj(d,k,b,j) = L_jbkd => L_dkbj(b,k,d,l) = L_ldkb
!
!     Reorder to L_ldkb = L_dkbj(b,k,d,l)
!
      call mem%alloc(L_ldkb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call sort_1234_to_4321(L_dkbj, L_ldkb, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(L_dkbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Form the intermediate X_cb = sum_kdl t_kl^cd L_ldkb
!                                 = sum_kdl t_c_ldk L_ldk_b
!
      call mem%alloc(X_cb, wf%n_v, wf%n_v)
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  t_cldk,               & ! t_c_ldk
                  wf%n_v,               &
                  L_ldkb,               & ! L_ldk_b
                  (wf%n_v)*(wf%n_o)**2, &
                  zero,                 &
                  X_cb,                 &
                  wf%n_v)
!
      call mem%dealloc(t_cldk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(L_ldkb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     Reorder to b_aijc = b_aicj
!
      call mem%alloc(b_aijc, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call sort_1234_to_1243(b_aibj, b_aijc, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Form and add - sum_ckdl b_aicj t_kl^cd L_ldkb
!                  = - sum_c b_aij_c X_cb
!
      call mem%alloc(sigma_aijb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N','N',               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_v,               &
                  wf%n_v,               &
                  -one,                 &
                  b_aijc,               & ! b_aij_c
                  (wf%n_v)*(wf%n_o)**2, &
                  X_cb,                 &
                  wf%n_v,               &
                  zero,                 &
                  sigma_aijb,           &
                  (wf%n_v)*(wf%n_o)**2)
!
      call mem%dealloc(X_cb, wf%n_v, wf%n_v)
      call mem%dealloc(b_aijc, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call add_1243_to_1234(one, sigma_aijb, sigma_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(sigma_aijb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
   end subroutine jacobian_transpose_ccsd_e2_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_f2_ccsd(wf, sigma_aibj, b_aibj)
!!
!!    Jacobian transpose CCSD F2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Calculates the F2 term,
!!
!!       - sum_ckdl (b_alck t_kl^cd L_jbid + b_ajck t_kl^cd L_ldib + b_djck t_kl^cd L_ialb)
!!
!!    and adds it to the transformed vector sigma_aibj.
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: sigma_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: t_lckd ! t_kl^cd
      real(dp), dimension(:,:,:,:), allocatable :: t_ckdl ! t_kl^cd
!
      real(dp), dimension(:,:,:,:), allocatable :: g_jbid
      real(dp), dimension(:,:,:,:), allocatable :: L_dibj ! L_jbid
      real(dp), dimension(:,:,:,:), allocatable :: L_dlbi ! L_ldib
!
      real(dp), dimension(:,:), allocatable :: X_ad   ! An intermediate, term 1
      real(dp), dimension(:,:,:,:), allocatable :: X_ckbi ! An intermediate, term 2
      real(dp), dimension(:,:), allocatable :: X_lj   ! An intermediate, term 3
!
      real(dp), dimension(:,:,:,:), allocatable :: sigma_ajbi ! sigma_aibj contribution
!
!     :: Term 1. - sum_ckdl b_alck t_kl^cd L_jbid ::
!
!     X_ad = b_a_lck t_lckd
!
!     Form t_lckd = t_kl^cd
!
      call mem%alloc(t_lckd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call squareup_and_sort_1234_to_4123(wf%t2, t_lckd, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Form the intermediate X_ad = sum_lck b_a_lck t_lck_d
!
      call mem%alloc(X_ad, wf%n_v, wf%n_v)
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  b_aibj,               & ! "b_a_lck"
                  wf%n_v,               &
                  t_lckd,               & ! t_lck_d
                  (wf%n_v)*(wf%n_o)**2, &
                  zero,                 &
                  X_ad,                 &
                  wf%n_v)
!
!     Form g_jbid
!
      call mem%alloc(g_jbid, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_ovov(g_jbid)
!
!     Form L_dibj = L_jbid = 2 * g_jbid - g_jdib
!                          = 2 * g_jbid(j,b,i,d) - g_jbid(j,d,i,b)
!
      call mem%alloc(L_dibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      L_dibj = zero
!
      call add_4321_to_1234(two, g_jbid, L_dibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call add_4123_to_1234(-one, g_jbid, L_dibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_jbid, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     Add - sum_ckdl b_alck t_kl^cd L_jbid
!         = - sum_d X_ad L_d_ibj
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_v,               &
                  -one,                 &
                  X_ad,                 &
                  wf%n_v,               &
                  L_dibj,               & ! L_d_ibj
                  wf%n_v,               &
                  one,                  &
                  sigma_aibj,           & ! "sigma_aibj"
                  wf%n_v)
!
      call mem%dealloc(X_ad, wf%n_v, wf%n_v)
!
!     :: Term 2. - sum_ckdl b_ajck t_kl^cd L_ldib ::
!
!     X_ckbi = sum_dl t_ck_dl L_dl_bi
!
!     Reorder to t_ckdl = t_lckd = t_kl^cd
!
      call mem%alloc(t_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call sort_1234_to_2341(t_lckd, t_ckdl, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%dealloc(t_lckd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     We have L_dibj = L_jbid => L_dibj(b,i,d,l) = L_ldib
!
!     Form L_dlbi = L_ldib = L_dibj(b,i,d,l)
!
      call mem%alloc(L_dlbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      L_dlbi = zero
!
      call sort_1234_to_3412(L_dibj, L_dlbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(L_dibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Form the intermediate X_ckbi = sum_dl t_ck_dl L_dl_bi
!
      call mem%alloc(X_ckbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',            &
                  (wf%n_v)*(wf%n_o), &
                  (wf%n_v)*(wf%n_o), &
                  (wf%n_v)*(wf%n_o), &
                  one,               &
                  t_ckdl,            & ! t_ck_dl
                  (wf%n_v)*(wf%n_o), &
                  L_dlbi,            & ! L_dl_bi
                  (wf%n_v)*(wf%n_o), &
                  zero,              &
                  X_ckbi,            &
                  (wf%n_v)*(wf%n_o))
!
!     Form and add - sum_ckdl b_ajck t_kl^cd L_ldib = - sum_ck b_aj_ck X_ck_bi
!
      call mem%alloc(sigma_ajbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',            &
                  (wf%n_v)*(wf%n_o), &
                  (wf%n_v)*(wf%n_o), &
                  (wf%n_v)*(wf%n_o), &
                  -one,              &
                  b_aibj,            & ! b_aj_ck
                  (wf%n_v)*(wf%n_o), &
                  X_ckbi,            & ! X_ck_bi
                  (wf%n_v)*(wf%n_o), &
                  zero,              &
                  sigma_ajbi,        &
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(X_ckbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call add_1432_to_1234(one, sigma_ajbi, sigma_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(sigma_ajbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     :: Term 3. - sum_ckdl b_djck t_kl^cd L_ialb = - sum_ckdl b_ckdj t_kl^cd L_ialb ::
!
!     X_lj = sum_ckd t_l_ckd b_ckd_j
!
!     Reorder to t_ckdl = t_kl^cd
!
      call mem%alloc(t_lckd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call sort_1234_to_4123(t_ckdl, t_lckd, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(t_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Form the intermediate X_lj = sum_ckd t_kl^cd b_ckdj
!                                 = sum_ckd t_l_ckd b_ckd_j
!
      call mem%alloc(X_lj, wf%n_o, wf%n_o)
!
      call dgemm('N','N',               &
                  wf%n_o,               &
                  wf%n_o,               &
                  (wf%n_o)*(wf%n_v)**2, &
                  one,                  &
                  t_lckd,               & ! t_l_ckd
                  wf%n_o,               &
                  b_aibj,               & ! b_ckd_j
                  (wf%n_o)*(wf%n_v)**2, &
                  zero,                 &
                  X_lj,                 &
                  wf%n_o)
!
      call mem%dealloc(t_lckd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     - sum_ckdl b_ckdj t_kl^cd L_ialb = - sum_l L_ialb X_lj
!
!     We have L_dlbi = (L_ldib) = L_aib_l
!     Add - sum_ckdl b_ckdj t_kl^cd L_ialb = - sum_l L_aib_l X_lj
!
      call dgemm('N','N',               &
                  (wf%n_o)*(wf%n_v)**2, &
                  wf%n_o,               &
                  wf%n_o,               &
                  -one,                 &
                  L_dlbi,               & ! L_dl_bi
                  (wf%n_o)*(wf%n_v)**2, &
                  X_lj,                 &
                  wf%n_o,               &
                  one,                  &
                  sigma_aibj,          & ! sigma_aib_j
                  (wf%n_o)*(wf%n_v)**2)
!
      call mem%dealloc(L_dlbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X_lj, wf%n_o, wf%n_o)
!
   end subroutine jacobian_transpose_ccsd_f2_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_g2_ccsd(wf, sigma_aibj, b_aibj)
!!
!!    Jacobian transpose CCSD G2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Calculates the G2 term,
!!
!!       sum_ckdl (b_alcj t_kl^cd g_kbid + b_ajcl t_kl^cd g_kdib)
!!
!!    and adds it to the transformed vector sigma_aibj.
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: sigma_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: t_cldk ! t_kl^cd
      real(dp), dimension(:,:,:,:), allocatable :: t_clkd ! t_kl^cd
!
      real(dp), dimension(:,:,:,:), allocatable :: g_kbid
      real(dp), dimension(:,:,:,:), allocatable :: g_dkbi ! g_kbid
      real(dp), dimension(:,:,:,:), allocatable :: g_kdib
!
      real(dp), dimension(:,:,:,:), allocatable :: b_ajcl     ! b_alcj
      real(dp), dimension(:,:,:,:), allocatable :: sigma_ajbi ! sigma_aibj contribution
      real(dp), dimension(:,:,:,:), allocatable :: sigma_ajib ! sigma_aibj contribution
!
      real(dp), dimension(:,:,:,:), allocatable :: X_clbi ! An intermediate, term 1
      real(dp), dimension(:,:,:,:), allocatable :: X_ajkd ! An intermediate, term 2
!
!     :: Term 1. sum_ckdl b_alcj t_kl^cd g_kbid ::
!
!     Form t_cldk = t_kl^cd
!
      call mem%alloc(t_cldk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call squareup_and_sort_1234_to_1432(wf%t2, t_cldk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Form g_kbid
!
      call mem%alloc(g_kbid, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_ovov(g_kbid)
!
!     Reorder to g_dkbi = g_kbid
!
      call mem%alloc(g_dkbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call sort_1234_to_4123(g_kbid, g_dkbi, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%dealloc(g_kbid, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     Form the intermediate X_clbi = sum_dk t_cl_dk g_dk_bi
!
      call mem%alloc(X_clbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',            &
                  (wf%n_v)*(wf%n_o), &
                  (wf%n_v)*(wf%n_o), &
                  (wf%n_v)*(wf%n_o), &
                  one,               &
                  t_cldk,            & ! t_cl_dk
                  (wf%n_v)*(wf%n_o), &
                  g_dkbi,            & ! g_dk_bi
                  (wf%n_v)*(wf%n_o), &
                  zero,              &
                  X_clbi,            &
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(t_cldk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(g_dkbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Reorder to b_ajcl = b_alcj
!
      call mem%alloc(b_ajcl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call sort_1234_to_1432(b_aibj, b_ajcl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Calculate and add sum_ckdl b_alcj t_kl^cd g_kbid
!                       = sum_cl b_aj_cl X_cl_bi
!
      call mem%alloc(sigma_ajbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',            &
                  (wf%n_v)*(wf%n_o), &
                  (wf%n_v)*(wf%n_o), &
                  (wf%n_v)*(wf%n_o), &
                  one,               &
                  b_ajcl,            & ! b_aj_cl
                  (wf%n_v)*(wf%n_o), &
                  X_clbi,            & ! X_cl_bi
                  (wf%n_v)*(wf%n_o), &
                  zero,              &
                  sigma_ajbi,        &
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(X_clbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(b_ajcl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call add_1432_to_1234(one, sigma_ajbi, sigma_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(sigma_ajbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     :: Term 2. sum_ckdl b_ajcl t_kl^cd g_kdib ::
!
!     Form t_clkd = t_kl^cd
!
      call mem%alloc(t_clkd, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call squareup_and_sort_1234_to_1423(wf%t2, t_clkd, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Form the intermediate X_ajkd = sum_cl b_aj_cl t_cl_kd
!
      call mem%alloc(X_ajkd, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N','N',            &
                  (wf%n_v)*(wf%n_o), &
                  (wf%n_v)*(wf%n_o), &
                  (wf%n_v)*(wf%n_o), &
                  one,               &
                  b_aibj,            & ! b_aj_cl
                  (wf%n_v)*(wf%n_o), &
                  t_clkd,            & ! t_cl_kd
                  (wf%n_v)*(wf%n_o), &
                  zero,              &
                  X_ajkd,            &
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(t_clkd, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
!     Form g_kdib
!
      call mem%alloc(g_kdib, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_ovov(g_kdib)
!
!     Form and add sum_ckdl b_ajcl t_kl^cd g_kdib
!                  = sum_kd X_aj_kd g_kd_ib
!
      call mem%alloc(sigma_ajib, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  X_ajkd,            & ! X_aj_kd
                  (wf%n_o)*(wf%n_v), &
                  g_kdib,            & ! g_kd_ib
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  sigma_ajib,        &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(g_kdib, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%dealloc(X_ajkd, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call add_1423_to_1234(one, sigma_ajib, sigma_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(sigma_ajib, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
   end subroutine jacobian_transpose_ccsd_g2_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_h2_ccsd(wf, sigma_abij, b_abij)
!!
!!    Jacobian transpose CCSD H2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Calculates the H2 term,
!!
!!       sum_kl b_akbl g_ikjl + sum_cd b_cidj g_cadb
!!
!!    and adds it to the transformed vector sigma_abij.
!!
!!    In this routine, the b and sigma vectors are ordered as
!!
!!       b_abij = b_aibj
!!       sigma_abij = sigma_abij
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o) :: sigma_abij
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o) :: b_abij
!
      real(dp), dimension(:,:,:,:), allocatable :: sigma_abij_batch
!
      real(dp), dimension(:,:,:,:), allocatable :: g_cadb
      real(dp), dimension(:,:,:,:), allocatable :: g_abcd ! g_cadb
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ikjl
      real(dp), dimension(:,:,:,:), allocatable :: g_klij ! g_ikjl
!
!     Batching variables
!
      integer :: rec2, rec0, rec1_a, rec1_b
!
      integer :: current_a_batch = 0
      integer :: current_b_batch = 0
!
      type(batching_index) :: batch_a
      type(batching_index) :: batch_b
!
!     :: Term 1. sum_kl b_akbl g_ikjl ::
!
!     Form g_ikjl
!
      call mem%alloc(g_ikjl, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call wf%get_oooo(g_ikjl)
!
!     Reorder to g_klij = g_ikjl
!
      call mem%alloc(g_klij, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call sort_1234_to_2413(g_ikjl, g_klij, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call mem%dealloc(g_ikjl, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
!     Add sum_kl b_akbl g_ikjl = sum_kl b_ab_kl g_kl_ij
!
      call dgemm('N','N',      &
                  (wf%n_v)**2, &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  one,         &
                  b_abij,      & ! b_ab_kl
                  (wf%n_v)**2, &
                  g_klij,      & ! g_kl_ij
                  (wf%n_o)**2, &
                  zero,        &
                  sigma_abij,  &
                  (wf%n_v)**2)
!
      call mem%dealloc(g_klij, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
!     :: Term 2. sum_cd b_cidj g_cadb ::
!
!     sum_cd b_cidj g_cadb = sum_cd g_cadb b_cd_ij
!                          = sum_cd g_ab_cd b_cd_ij
!
!     Prepare batching over a and b
!
      rec0 = 0
      rec1_a = wf%n_v*wf%integrals%n_J
      rec1_b = wf%n_v*wf%integrals%n_J
      rec2 = 2*wf%n_v**2 + wf%n_o**2
!
!     Initialize batching indices
!
      call batch_a%init(wf%n_v)
      call batch_b%init(wf%n_v)
!
      call mem%batch_setup(batch_a, batch_b, rec0, rec1_a, rec1_b, rec2)
!
!     Loop over a-batches
!
      do current_a_batch = 1, batch_a%num_batches
!
!        For each a batch, get the limits for the a index
!
         call batch_a%determine_limits(current_a_batch)
!
!        Loop over b-batches
!
         do current_b_batch = 1, batch_b%num_batches
!
!           For each b batch, get the limits for the b index
!
            call batch_b%determine_limits(current_b_batch)
!
!           Form g_cadb
!
            call mem%alloc(g_cadb, wf%n_v, batch_a%length, wf%n_v, batch_b%length)
!
            call wf%get_vvvv(g_cadb,         &
                              1,             &
                              wf%n_v,        &
                              batch_a%first, &
                              batch_a%last,  &
                              1,             &
                              wf%n_v,        &
                              batch_b%first, &
                              batch_b%last)
!
!           Reorder to g_abcd = g_cadb
!
            call mem%alloc(g_abcd, batch_a%length, batch_b%length, wf%n_v, wf%n_v)
!
            call sort_1234_to_2413(g_cadb, g_abcd, wf%n_v, batch_a%length, wf%n_v, batch_b%length)
!
            call mem%dealloc(g_cadb, wf%n_v, batch_a%length, wf%n_v, batch_b%length)
!
!           Calculate sigma_abij_batch = sum_cd g_abcd b_cd_ij
!           and add it to the full space sigma vector
!
            call mem%alloc(sigma_abij_batch, batch_a%length, batch_b%length, wf%n_o, wf%n_o)
!
            call dgemm('N','N',                            &
                        (batch_a%length)*(batch_b%length), &
                        (wf%n_o)**2,                       &
                        (wf%n_v)**2,                       &
                        one,                               &
                        g_abcd,                            & ! g_ab_cd
                        (batch_a%length)*(batch_b%length), &
                        b_abij,                            & ! "b_cd_ij"
                        (wf%n_v)**2,                       &
                        zero,                              &
                        sigma_abij_batch,                  &
                        (batch_a%length)*(batch_b%length))
!
            call mem%dealloc(g_abcd, batch_a%length, batch_b%length, wf%n_v, wf%n_v)
!
            call daxpy(batch_a%length*batch_b%length*wf%n_o**2, one, sigma_abij_batch, 1, sigma_abij, 1)
!
            call mem%dealloc(sigma_abij_batch, batch_a%length, batch_b%length, wf%n_o, wf%n_o)
!
         enddo ! End of batches over b
      enddo ! End of batches over a
!
   end subroutine jacobian_transpose_ccsd_h2_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_i2_ccsd(wf, sigma_abij, b_abij)
!!
!!    Jacobian transpose CCSD I2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Calculates the I2 term,
!!
!!       sum_ckdl b_cidj t_kl^cd g_kalb + sum_ckdl b_akbl t_kl^cd g_icjd
!!
!!    and adds it to the transformed vector sigma_abij.
!!
!!    In this routine, the b and sigma vectors are ordered as
!!
!!       b_abij = b_aibj
!!       sigma_abij = sigma_abij
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o) :: sigma_abij
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o) :: b_abij
!
      real(dp), dimension(:,:,:,:), allocatable :: t_klcd ! t_kl^cd
!
      real(dp), dimension(:,:,:,:), allocatable :: g_kalb ! g_kalb
      real(dp), dimension(:,:,:,:), allocatable :: g_abkl ! g_kalb
!
      real(dp), dimension(:,:,:,:), allocatable :: X_klij ! An intermediate, terms 1 & 2
!
!     :: Term 1. sum_ckdl b_cidj t_kl^cd g_kalb ::
!
!     sum_ckdl t_kl_cd b_cd_ij
!
!     Form t_kl_cd = t_kl^cd
!
      call mem%alloc(t_klcd, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
      call squareup_and_sort_1234_to_2413(wf%t2, t_klcd, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Form the intermediate X_klij = sum_cd t_kl_cd b_cd_ij
!
      call mem%alloc(X_klij, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call dgemm('N','N',      &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  (wf%n_v)**2, &
                  one,         &
                  t_klcd,      & ! t_kl_cd
                  (wf%n_o)**2, &
                  b_abij,      & ! "b_cd_ij"
                  (wf%n_v)**2, &
                  zero,        &
                  X_klij,      &
                  (wf%n_o)**2)
!
!     Form g_kalb
!
      call mem%alloc(g_kalb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_ovov(g_kalb)
!
!     Reorder to g_abkl = g_kalb
!
      call mem%alloc(g_abkl, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      g_abkl = zero
!
      call sort_1234_to_2413(g_kalb, g_abkl, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%dealloc(g_kalb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     Add sum_ckdl b_cidj t_kl^cd g_kalb
!         = sum_kl g_ab_kl X_kl_ij
!
      call dgemm('N','N',      &
                  (wf%n_v)**2, &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  one,         &
                  g_abkl,      &
                  (wf%n_v)**2, &
                  X_klij,      &
                  (wf%n_o)**2, &
                  one,         &
                  sigma_abij,  &
                  (wf%n_v)**2)
!
!     :: Term 2. sum_ckdl b_akbl t_kl^cd g_icjd ::
!
!     Repurpose X_kl_ij to make sum_cd t_kl^cd g_icjd
!                               = sum_cd t_kl_cd g_cd_ij
!                               = sum_cd t_kl_cd g_ab_kl(cd,ij)
!
      call dgemm('N','N',      &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  (wf%n_v)**2, &
                  one,         &
                  t_klcd,      & ! t_kl_cd
                  (wf%n_o)**2, &
                  g_abkl,      & ! "g_cd_ij"
                  (wf%n_v)**2, &
                  zero,        &
                  X_klij,      &
                  (wf%n_o)**2)
!
      call mem%dealloc(g_abkl, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(t_klcd, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
!     Add sum_ckdl b_akbl t_kl^cd g_icjd
!         = sum_kl b_ab_kl X_kl_ij
!
      call dgemm('N','N',      &
                  (wf%n_v)**2, &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  one,         &
                  b_abij,      & ! "b_ab_kl"
                  (wf%n_v)**2, &
                  X_klij,      & ! X_kl_ij
                  (wf%n_o)**2, &
                  one,         &
                  sigma_abij,  &
                  (wf%n_v)**2)
!
      call mem%dealloc(X_klij, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
   end subroutine jacobian_transpose_ccsd_i2_ccsd
!
!
end submodule jacobian_transpose_ccsd
