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
submodule (cc2_class) jacobian_transpose_cc2
!
!!
!!    Jacobian transpose submodule (CC2)
!!    Written by Sarai D. Folkestad and Alexander Paul, Feb 2019
!!
!!    Routines for the linear transform of trial
!!    vectors by the transpose of the Jacobian matrix
!!
!!    σ_i = A^T * c_i,
!!
!!    where
!!
!!    A_μ,ν = < μ | exp(-T) [H, τ_ν] exp(T) | R >.
!!
!!
!
   implicit none
!
!
contains
!
!
   module subroutine prepare_for_jacobian_transpose_cc2(wf)
!!
!!    Jacobian transpose submodule (CC2)
!!    Written by Sarai D. Folkestad and Alexander Paul, Feb 2019
!!
      implicit none
!
      class(cc2), intent(inout) :: wf
!
      call wf%initialize_u()
      call wf%construct_u()
!
   end subroutine prepare_for_jacobian_transpose_cc2
!
!
!
   module subroutine jacobian_transpose_transform_trial_vector_cc2(wf, c_i)
!!
!!    Jacobian transform trial vector
!!    Written by Sarai D. Folkestad and Alexander Paul, Feb 2019
!!
      class(cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes) :: c_i
!
      call wf%jacobian_transpose_cc2_transformation(c_i)
!
   end subroutine jacobian_transpose_transform_trial_vector_cc2
!
!
   module subroutine jacobian_transpose_cc2_transformation_cc2(wf, c)
!!
!!    Jacobian transpose transformation (CC2)
!!    Written by Sarai D. Folkestad and Alexander Paul, Feb 2019
!!
!!    Calculates the transpose Jacobian transformation, i.e., the transformation
!!    by the transpose of the Jacobian matrix
!!
!!       A_mu,nu = < mu | exp(-T) [H, tau_nu] exp(T) | R >.
!!
!!    The transformation is performed as sigma^T = c^T A, where c is the vector
!!    sent to the routine. On exit, the vector c is equal to sigma (the transformed
!!    vector).
!!
      implicit none
!
      class(cc2) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes) :: c
!
      real(dp), dimension(:,:), allocatable :: c_ai
      real(dp), dimension(:,:,:,:), allocatable :: c_aibj
!
      real(dp), dimension(:,:), allocatable :: sigma_ai
      real(dp), dimension(:,:,:,:), allocatable :: sigma_aibj
!
      integer :: i, j, a, b, ai, bj, aibj ! Index
!
!     Allocate and zero the transformed vecotr (singles part)
!
      call mem%alloc(sigma_ai, wf%n_v, wf%n_o)
      sigma_ai = zero
!
      call mem%alloc(c_ai, wf%n_v, wf%n_o)
!
!$omp parallel do schedule(static) private(a, i, ai)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            ai = wf%n_v*(i - 1) + a
!
            c_ai(a, i) = c(ai)
!
         enddo
      enddo
!$omp end parallel do
!
!     :: CCS contributions to the singles c vector ::
!
      call wf%jacobian_transpose_ccs_a1(sigma_ai, c_ai)
      call wf%jacobian_transpose_ccs_b1(sigma_ai, c_ai)
!
!     :: CC2 contributions to the transformed singles vector ::
!
      call wf%jacobian_transpose_cc2_a1(sigma_ai, c_ai)
!
!     Allocate the incoming unpacked doubles vector
!
      call mem%alloc(c_aibj, (wf%n_v), (wf%n_o), (wf%n_v), (wf%n_o))
!
!$omp parallel do schedule(static) private(a, i, b, j, ai, bj, aibj)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            ai = wf%n_v*(i - 1) + a
!
            do j = 1, wf%n_o
               do b = 1, wf%n_v
!
                  bj = wf%n_v*(j - 1) + b
!
                  if (ai .ge. bj) then
!
                     aibj = ai*(ai-3)/2 + ai + bj
!
                     c_aibj(a, i, b, j) = c(wf%n_o*wf%n_v + aibj)
                     c_aibj(b, j, a, i) = c(wf%n_o*wf%n_v + aibj)
!
                  endif
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call wf%jacobian_transpose_cc2_b1(sigma_ai, c_aibj)
!
!     Done with singles vector c; overwrite it with
!     transformed vector for exit
!
!$omp parallel do schedule(static) private(a, i, ai)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            ai = wf%n_v*(i - 1) + a
!
            c(ai) = sigma_ai(a, i)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(sigma_ai, wf%n_v, wf%n_o)
!
!     :: CC2 contributions to the transformed doubles vector ::
!
!     Allocate unpacked transformed vector
!
      call mem%alloc(sigma_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      sigma_aibj = zero
!
!     Contributions from singles vector c
!
      call wf%jacobian_transpose_cc2_a2(sigma_aibj, c_ai)
!
      call mem%dealloc(c_ai, wf%n_v, wf%n_o)
!
!     Contributions from doubles vector c
!
      call wf%jacobian_transpose_cc2_b2(sigma_aibj, c_aibj)
!
      call mem%dealloc(c_aibj, (wf%n_v), (wf%n_o), (wf%n_v), (wf%n_o))
!
!     Overwrite the incoming doubles c vector & pack in
!
!$omp parallel do schedule(static) private(a, i, b, j, ai, bj, aibj)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            ai = wf%n_v*(i - 1) + a
!
            do j = 1, wf%n_o
               do b = 1, wf%n_v
!
                  bj = wf%n_v*(j - 1) + b
!
                  if (ai .ge. bj) then
!
                     aibj = ai*(ai-3)/2 + ai + bj
!
                     c((wf%n_o)*(wf%n_v) + aibj) = sigma_aibj(a, i, b, j)
!
                  endif
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(sigma_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine jacobian_transpose_cc2_transformation_cc2
!
!
   module subroutine jacobian_transpose_cc2_a1_cc2(wf, sigma_ai, c_bj)
!!
!!    Jacobian transpose CC2 A1
!!    Written by Sarai D. Folkestad and Alexander Paul, Feb 2019
!!
!!    Calculates the A1 term,
!!
!!    u_kcbj = u^bc_jk =  2 t^bc_jk - t^bc_kj = -(2 g_bjck - g_bkcj)/eps^bc_jk
!!
!!    sigma_ai =+ sum_bjck u^bc_jk (c_bj L_iakc - c_ak g_jbic - c_ci g_jbka)
!!             =+ sum_ck X_kc L_iakc - sum_k c_ak Y_ik - sum_c Y_ca c_ci
!!
      implicit none
!
      class(cc2) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)     :: c_bj
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)  :: sigma_ai
!
      real(dp), dimension(:,:,:,:), allocatable :: g_iakc, g_icjb, g_jbka
      real(dp), dimension(:,:,:,:), allocatable :: L_aick
      real(dp), dimension(:,:,:,:), allocatable :: u_cjbk
!
      real(dp), dimension(:,:), allocatable :: X_ck, Y_ik, Y_ca
!
!     :: Term 1: sigma_ai =+ sum_bjck c_bj u^bc_jk L_iakc
!
      call mem%alloc(X_ck, wf%n_v, wf%n_o)
!
!     X_ck = sum_bj u_ckbj * c_bj
!
      call dgemm('N', 'N',             &
                  (wf%n_o)*(wf%n_v),   &
                  1,                   &
                  (wf%n_o)*(wf%n_v),   &
                  one,                 &
                  wf%u,                & ! u_ck_bj
                  (wf%n_o)*(wf%n_v),   &
                  c_bj,                & ! c_bj
                  (wf%n_o)*(wf%n_v),   &
                  zero,                &
                  X_ck,                &
                  (wf%n_o)*(wf%n_v))
!
!     L_iakc = 2 g_iakc - g_icka
!
      call mem%alloc(g_iakc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_ovov(g_iakc)
!
      call mem%alloc(L_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      L_aick = zero
      call add_2143_to_1234(two, g_iakc, L_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call add_2341_to_1234(-one, g_iakc, L_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_iakc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     sigma_ai =+ sum_kc L_iakc * X_ck
!
      call dgemm('N', 'N',             &
                  (wf%n_v)*(wf%n_o),   &
                  1,                   &
                  (wf%n_v)*(wf%n_o),   &
                  one,                 &
                  L_aick,              & ! L_ai_ck
                  (wf%n_v)*(wf%n_o),   &
                  X_ck,                & ! X_ck
                  (wf%n_v)*(wf%n_o),   &
                  one,                 &
                  sigma_ai,            &
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(L_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X_ck, wf%n_v, wf%n_o)
!
!     :: Term 2: sigma_ai =+ sum_bjck - c_ak u^bc_jk g_jbic
!
      call mem%alloc(g_icjb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_ovov(g_icjb)
!
      call mem%alloc(u_cjbk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_3214(wf%u, u_cjbk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(Y_ik, wf%n_o, wf%n_o)
!
!     Y_ik = sum_cjb g_jbic * u_cjbk = sum_cjb g_icjb * u_cjbk
!
      call dgemm('N', 'N',                &
                  wf%n_o,                 &
                  wf%n_o,                 &
                  (wf%n_o)*(wf%n_v)**2,   &
                  one,                    &
                  g_icjb,                 & ! g_i_cjb
                  wf%n_o,                 &
                  u_cjbk,                 & ! u_cjb_k
                  (wf%n_o)*(wf%n_v)**2,   &
                  zero,                   &
                  Y_ik,                   & ! Y_ik
                  wf%n_o)
!
      call mem%dealloc(g_icjb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     sigma_ai =+ - c_ak * Y_ik
!
      call dgemm('N', 'T',    & ! transpose of Y_ik
                  wf%n_v,     &
                  wf%n_o,     &
                  wf%n_o,     &
                  -one,       &
                  c_bj,       & ! c_a_k
                  wf%n_v,     &
                  Y_ik,       & ! Y_k_i
                  wf%n_o,     &
                  one,        &
                  sigma_ai,   &
                  wf%n_v)
!
      call mem%dealloc(Y_ik, wf%n_o, wf%n_o)
!
!     :: Term 3: sigma_ai =+ sum_bjck - c_ci u^bc_jk g_jbka
!
      call mem%alloc(g_jbka, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_ovov(g_jbka)
!
      call mem%alloc(Y_ca, wf%n_v, wf%n_v)
!
!     Y_ca = sum_jbk u_cjbk * g_jbka
!
      call dgemm('N', 'N',                &
                  wf%n_v,                 &
                  wf%n_v,                 &
                  (wf%n_v)*(wf%n_o)**2,   &
                  one,                    &
                  u_cjbk,                 & ! u_c_jbk
                  wf%n_v,                 &
                  g_jbka,                 & ! g_jbk_a
                  (wf%n_v)*(wf%n_o)**2,   &
                  zero,                   &
                  Y_ca,                   & ! Y_ca
                  wf%n_v)
!
      call mem%dealloc(g_jbka, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%dealloc(u_cjbk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     sigma_ai =+ - c_ci * Y_ca
!
      call dgemm('T', 'N',    & ! transpose of Y_ik
                  wf%n_v,     &
                  wf%n_o,     &
                  wf%n_v,     &
                  -one,       &
                  Y_ca,       & ! Y_a_c
                  wf%n_v,     &
                  c_bj,       & ! c_c_i
                  wf%n_v,     &
                  one,        &
                  sigma_ai,   &
                  wf%n_v)
!
      call mem%dealloc(Y_ca, wf%n_v, wf%n_v)
!
   end subroutine jacobian_transpose_cc2_a1_cc2
!
!
   module subroutine jacobian_transpose_cc2_b1_cc2(wf, sigma_ai, c_bjck)
!!
!!    Jacobian transpose CC2 B1
!!    Written by Sarai D. Folkestad and Alexander Paul, Feb 2019
!!
!!    Calculates the B1 term,
!!
!!    sigma_ai =+ sum_bjc c_bjci g_bjca - c_akbj g_bjik
!!
      implicit none
!
      class(cc2) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_bjck
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)  :: sigma_ai
!
      real(dp), dimension(:,:,:,:), allocatable :: g_bjca, g_ikbj
!
      type(batching_index) :: batch_a
!
      integer :: req0, req1, current_a_batch
!
!     :: Term 1: sigma_ai =+ sum_bjc c_bjci g_bjca = sum_bjc (g_bjca)^T c_bjci
!
      req0 = (wf%n_v)*(wf%n_o)*(wf%integrals%n_J)
      req1 = max((wf%n_v)*(wf%integrals%n_J) + (wf%n_o)*(wf%n_v)**2, 2*(wf%n_o)*(wf%n_v)**2)
!
      call batch_a%init(wf%n_v)
!
      call mem%batch_setup(batch_a, req0, req1)
!
      do current_a_batch = 1, batch_a%num_batches
!
         call batch_a%determine_limits(current_a_batch)
!
         call mem%alloc(g_bjca, wf%n_v, wf%n_o, wf%n_v, batch_a%length)
!
         call wf%get_vovv(g_bjca,                        &
                           1, wf%n_v,                    &
                           1, wf%n_o,                    &
                           1, wf%n_v,                    &
                           batch_a%first, batch_a%last)
!
!        sigma_ai =+ sum_bjc g_abjc * c_bjci
!
         call dgemm('T', 'N',                    & ! transposed g_bjca
                     wf%n_v,                     &
                     wf%n_o,                     &
                     (wf%n_o)*(wf%n_v)**2,       &
                     one,                        &
                     g_bjca,                     & ! g_a_bjc
                     (wf%n_o)*(wf%n_v)**2,       &
                     c_bjck,                     & ! c_bjc_i
                     (wf%n_o)*(wf%n_v)**2,       &
                     one,                        &
                     sigma_ai(batch_a%first, 1), &
                     wf%n_v)
!
         call mem%dealloc(g_bjca, wf%n_v, wf%n_o, wf%n_v, batch_a%length)
!
      enddo ! batch_a
!
!     :: Term 2: sigma_ai =+ sum_bjc c_akbj g_bjik = sum_bjc c_akbj (g_ikbj)^T
!
      call mem%alloc(g_ikbj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%get_oovo(g_ikbj)
!
!     sigma_ai =+ sum_bjk - c_akbj g_ikbj
!
      call dgemm('N', 'T',                & ! transposed g_ikbj
                  wf%n_v,                 &
                  wf%n_o,                 &
                  (wf%n_v)*(wf%n_o)**2,   &
                  -one,                   &
                  c_bjck,                 & ! c_a_kbj
                  wf%n_v,                 &
                  g_ikbj,                 & ! g_kbj_i
                  wf%n_o,                 &
                  one,                    &
                  sigma_ai,               &
                  wf%n_v)
!
      call mem%dealloc(g_ikbj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine jacobian_transpose_cc2_b1_cc2
!
!
   module subroutine jacobian_transpose_cc2_a2_cc2(wf, sigma_aibj, c_ai)
!!
!!    Jacobian transpose CC2 A2
!!    Written by Sarai D. Folkestad and Alexander Paul, Feb 2019
!!
!!    Calculates the A2 term,
!!
!!    sigma_aibj =+ P_ai,bj (2F_jb c_ai - F_ib c_aj - L_ikjb c_ak + L_cajb c_ci)
!!
      implicit none
!
      class(cc2) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                      :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout)   :: sigma_aibj
!
!     Local variables
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ikjb, g_cajb
      real(dp), dimension(:,:,:,:), allocatable :: L_kibj, L_cajb
      real(dp), dimension(:,:,:,:), allocatable :: sigma_iajb
!
      type(batching_index) :: batch_c
!
      integer :: req0, req1, current_c_batch
!
      integer :: a, i, b, j
!
!     Term 1: (2F_jb c_ai - F_ib c_aj)
!
      sigma_aibj = zero
!
!$omp parallel do private(a, i, b, j)
      do b = 1, wf%n_v
         do j = 1, wf%n_o
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  sigma_aibj(a, i, b, j) = sigma_aibj(a, i, b, j) &
                                          +  two*wf%fock_ia(j, b)*c_ai(a, i)&
                                          -  wf%fock_ia(i, b)*c_ai(a, j)

!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
!     Term 3: - L_ikjb c_ak
!
!     L_ikjb = 2 g_ikjb - g_jkib (ordered as g_kibj)
!
      call mem%alloc(g_ikjb, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call wf%get_ooov(g_ikjb)
!
      call mem%alloc(L_kibj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      L_kibj = zero
      call add_2143_to_1234(two, g_ikjb, L_kibj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call add_4123_to_1234(-one, g_ikjb, L_kibj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_ikjb, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N', 'N',             &
                  wf%n_v,              &
                  (wf%n_o**2)*wf%n_v,  &
                  wf%n_o,              &
                  -one,                &
                  c_ai,                & ! c_a_k
                  wf%n_v,              &
                  L_kibj,              & ! L_k_ibj
                  wf%n_o,              &
                  one,                 &
                  sigma_aibj,          & ! sigma_a_ibj
                  wf%n_v)
!
      call mem%dealloc(L_kibj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     Term 4: L_cajb c_ci
!
      call mem%alloc(sigma_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      sigma_iajb = zero
!
      req0 = (wf%n_v)*(wf%n_o)*(wf%integrals%n_J)
      req1 = (wf%n_v)*(wf%integrals%n_J) + (wf%n_o)*(wf%n_v)**2
!
      call batch_c%init(wf%n_v)
!
      call mem%batch_setup(batch_c, req0, req1)
!
      do current_c_batch = 1, batch_c%num_batches
!
         call batch_c%determine_limits(current_c_batch)
!
!        L_cajb = 2 g_cajb - g_cbja
!
         call mem%alloc(g_cajb, batch_c%length, wf%n_v, wf%n_o, wf%n_v)
!
         call wf%get_vvov(g_cajb, &
                           batch_c%first, batch_c%last,   &
                           1, wf%n_v,                     &
                           1, wf%n_o,                     &
                           1, wf%n_v)
!
         call mem%alloc(L_cajb, batch_c%length, wf%n_v, wf%n_o, wf%n_v)
         call dcopy((batch_c%length)*(wf%n_v**2)*(wf%n_o), g_cajb, 1, L_cajb, 1)
         call dscal((batch_c%length)*(wf%n_v**2)*(wf%n_o), two, L_cajb, 1)
         call add_1432_to_1234(-one, g_cajb, L_cajb, batch_c%length, wf%n_v, wf%n_o, wf%n_v)
!
         call mem%dealloc(g_cajb, batch_c%length, wf%n_v, wf%n_o, wf%n_v)
!
         call dgemm('T', 'N',                &
                     wf%n_o,                 &
                     (wf%n_v**2)*(wf%n_o),   &
                     batch_c%length,         &
                     one,                    &
                     c_ai(batch_c%first, 1), & ! c_c_i
                     wf%n_v,                 &
                     L_cajb,                 & ! L_c_ajb
                     batch_c%length,         &
                     one,                    &
                     sigma_iajb,             &
                     wf%n_o)
!
         call mem%dealloc(L_cajb, batch_c%length, wf%n_v, wf%n_o, wf%n_v)
!
      enddo ! batch_c
!
      call add_2143_to_1234(one, sigma_iajb, sigma_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(sigma_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     Symmetrize
!
      call symmetric_sum(sigma_aibj, (wf%n_o)*(wf%n_v))
!
   end subroutine jacobian_transpose_cc2_a2_cc2
!
!
   module subroutine jacobian_transpose_cc2_b2_cc2(wf, sigma_aibj, c_aibj)
!!
!!    Jacobian transpose CC2 B2
!!    Written by Sarai D. Folkestad and Alexander Paul, Feb 2019
!!
!!    Calculates the A2 term,
!!
!!    sigma_aibj =+ ε_aibj c_aibj
!!
      implicit none
!
      class(cc2) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)      :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout)   :: sigma_aibj
!
      integer :: a, i, b, j
!
!$omp parallel do private(a, i, b, j)
      do b = 1, wf%n_v
         do j = 1, wf%n_o
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  sigma_aibj(a, i, b, j) = sigma_aibj(a, i, b, j) + c_aibj(a,i,b,j) &
                                          * (wf%orbital_energies(a + wf%n_o) &
                                          + wf%orbital_energies(b + wf%n_o) &
                                          - wf%orbital_energies(i) &
                                          - wf%orbital_energies(j))
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine jacobian_transpose_cc2_b2_cc2
!
!
end submodule jacobian_transpose_cc2
