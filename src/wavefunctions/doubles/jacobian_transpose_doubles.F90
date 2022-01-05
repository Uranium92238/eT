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
submodule (doubles_class) jacobian_transpose_doubles
!
!!
!!    Jacobian transpose submodule
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
   module subroutine save_jacobian_transpose_a1_intermediates_doubles(wf, u_bjck)
!!
!!    Save jacobian transpose A1 intermediates
!!    Written by by E. F. Kjønstad, S. D. Folkestad and Alexander C. Paul
!!
!!    Calculates the intermediates,
!!
!!       Y_ik = sum_cjb g_icjb * u_bjck
!!       Y_ca = sum_jbk u_bjck * g_jbka
!!
!!    and saves them into
!!
!!       jacobian_transpose_a1_intermdiate_oo
!!       jacobian_transpose_a1_intermdiate_vv
!!
!!    u_bjck = u^bc_jk =  2 t^bc_jk - t^bc_kj
!!           = -(2 g_bjck - g_bkcj)/eps^bc_jk
!!
!!    Adapted by Tor S. Haugland, Oct 2019
!!
!!    Isolated the intermediates from the
!!    jacobian_transpose_doubles_a1_doubles and wrote them to file.
!!
      use reordering, only: sort_1234_to_3214
!
      implicit none
!
      class(doubles), intent(inout) :: wf
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: u_bjck
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ovov
      real(dp), dimension(:,:,:,:), allocatable :: u_cjbk
!
      real(dp), dimension(:,:),     allocatable :: Y_ik
      real(dp), dimension(:,:),     allocatable :: Y_ca
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian transpose A1 intermediate', pl='verbose')
      call timer%turn_on()
!
!     Both intermediates need g_ovov and u_bjck -> u_cjbk
!
      call mem%alloc(g_ovov, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call wf%eri%get_eri_t1('ovov', g_ovov)
!
      call mem%alloc(u_cjbk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_3214(u_bjck, u_cjbk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!        Y_ik = sum_cjb g_i_cjb * u_cjb_k
!
      call mem%alloc(Y_ik, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N',                &
                  wf%n_o,                 &
                  wf%n_o,                 &
                  (wf%n_o)*(wf%n_v)**2,   &
                  one,                    &
                  g_ovov,                 & ! g_i_cjb
                  wf%n_o,                 &
                  u_cjbk,                 & ! u_cjb_k
                  (wf%n_o)*(wf%n_v)**2,   &
                  zero,                   &
                  Y_ik,                   & ! Y_ik
                  wf%n_o)
!
!     Save Y_ik
!
      wf%jacobian_transpose_a1_intermediate_oo = &
                                          stream_file('jacobian_transpose_intermediate_a1_oo')
      call wf%jacobian_transpose_a1_intermediate_oo%open_('write', 'rewind')
!
      call wf%jacobian_transpose_a1_intermediate_oo%write_(Y_ik, wf%n_o**2)
!
      call wf%jacobian_transpose_a1_intermediate_oo%close_('keep')
!
      call mem%dealloc(Y_ik, wf%n_o, wf%n_o)
!
!        Y_ca = sum_jbk u_c_jbk g_jbk_a
!
      call mem%alloc(Y_ca, wf%n_v, wf%n_v)
!
      call dgemm('N', 'N',                &
                  wf%n_v,                 &
                  wf%n_v,                 &
                  (wf%n_v)*(wf%n_o)**2,   &
                  one,                    &
                  u_cjbk,                 & ! u_c_jbk
                  wf%n_v,                 &
                  g_ovov,                 & ! g_jbk_a
                  (wf%n_v)*(wf%n_o)**2,   &
                  zero,                   &
                  Y_ca,                   & ! Y_ca
                  wf%n_v)
!
!     Save Y_ca
!
      wf%jacobian_transpose_a1_intermediate_vv = &
                                          stream_file('jacobian_transpose_intermediate_a1_vv')
      call wf%jacobian_transpose_a1_intermediate_vv%open_('write', 'rewind')
!
      call wf%jacobian_transpose_a1_intermediate_vv%write_(Y_ca, wf%n_v**2)
!
      call wf%jacobian_transpose_a1_intermediate_vv%close_('keep')
!
      call mem%dealloc(Y_ca, wf%n_v, wf%n_v)
!
!     Cleanup
!
      call mem%dealloc(g_ovov, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%dealloc(u_cjbk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine save_jacobian_transpose_a1_intermediates_doubles
!
!
   module subroutine jacobian_transpose_doubles_a1_doubles(wf, sigma_ai, c_bj, u)
!!
!!    Jacobian transpose doubles A1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!    and Alexander C. Paul, Feb 2019
!!
!!    Calculates the A1 term,
!!
!!    u_ckbj = u^bc_jk =  2 t^bc_jk - t^bc_kj = -(2 g_bjck - g_bkcj)/eps^bc_jk
!!
!!    sigma_ai += sum_bjck u^bc_jk (c_bj L_iakc - c_ak g_jbic - c_ci g_jbka)
!!             += sum_ck X_kc L_iakc - sum_k c_ak Y_ik - sum_c Y_ca c_ci
!!
!!    Modified by Tor S. Haugland, Oct 2019
!!
!!    Use saved intermediates to construct Y_ik and Y_ca.
!!    Create intermediate X_kc using transpose to save time re-ordering g_iakc
!!
      use array_utilities, only: copy_and_scale
      use reordering, only: sort_12_to_21, add_1432_to_1234, sort_12_to_21
!
      implicit none
!
      class(doubles) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)     :: c_bj
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)  :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: u
!
      real(dp), dimension(:,:,:,:), allocatable :: g_iakc
      real(dp), dimension(:,:,:,:), allocatable :: L_iakc
!
      real(dp), dimension(:,:), allocatable :: X_ck, X_kc, Y_ik, Y_ca, sigma_ia, sigma_ai_temp
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian transpose doubles A1', pl='verbose')
      call timer%turn_on()
!
!     :: Term 1: sigma_ai += sum_bjck c_bj u^bc_jk L_iakc
!
!     Intermediate X_kc = (X_ck)^T = ( sum_bj u_ckbj * c_bj )^T
!
      call mem%alloc(X_ck, wf%n_v, wf%n_o)
!
      call dgemv('N',                  &
                  (wf%n_o)*(wf%n_v),   &
                  (wf%n_o)*(wf%n_v),   &
                  one,                 &
                  u,                   & ! u_ck_bj
                  (wf%n_o)*(wf%n_v),   &
                  c_bj,                & ! c_bj
                  1,                   &
                  zero,                &
                  X_ck,                &
                  1)
!
      call mem%alloc(X_kc, wf%n_o, wf%n_v)
!
      call sort_12_to_21(X_ck, X_kc, wf%n_v, wf%n_o)
!
      call mem%dealloc(X_ck, wf%n_v, wf%n_o)
!
!     L_iakc = 2 g_iakc - g_icka
!
      call mem%alloc(g_iakc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%alloc(L_iakc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call wf%eri%get_eri_t1('ovov', g_iakc)
!
      call copy_and_scale(two, g_iakc, L_iakc, wf%n_t1**2)
      call add_1432_to_1234(-one, g_iakc, L_iakc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%dealloc(g_iakc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     sigma_ai += sum_kc L_iakc * X_kc
!
      call mem%alloc(sigma_ia, wf%n_o, wf%n_v)
!
      call dgemm('N', 'N',             &
                  (wf%n_v)*(wf%n_o),   &
                  1,                   &
                  (wf%n_v)*(wf%n_o),   &
                  one,                 &
                  L_iakc,              & ! L_ia_kc
                  (wf%n_v)*(wf%n_o),   &
                  X_kc,                & ! X_kc
                  (wf%n_v)*(wf%n_o),   &
                  zero,                &
                  sigma_ia,            &
                  (wf%n_v)*(wf%n_o))
!
      call mem%alloc(sigma_ai_temp, wf%n_v, wf%n_o)
      call sort_12_to_21(sigma_ia, sigma_ai_temp, wf%n_o, wf%n_v)
!
      call daxpy(wf%n_v * wf%n_o, one, sigma_ai_temp, 1, sigma_ai, 1)
!
      call mem%dealloc(sigma_ai_temp, wf%n_v, wf%n_o)
!
!     Cleanup
!
      call mem%dealloc(X_kc, wf%n_o, wf%n_v)
      call mem%dealloc(L_iakc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%dealloc(sigma_ia, wf%n_o, wf%n_v)
!
!     :: Term 2: sigma_ai -= c_ak Y_ik
!
!     Read Y_ik from file
!
      call mem%alloc(Y_ik, wf%n_o, wf%n_o)
!
      call wf%jacobian_transpose_a1_intermediate_oo%open_('read', 'rewind')
      call wf%jacobian_transpose_a1_intermediate_oo%read_(Y_ik, wf%n_o**2)
      call wf%jacobian_transpose_a1_intermediate_oo%close_()
!
!     sigma_ai -= c_ak * Y_ik
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
!     :: Term 3: sigma_ai += - c_ci Y_ca
!
      call mem%alloc(Y_ca, wf%n_v, wf%n_v)
!
!     Read Y_ca from file
!
      call wf%jacobian_transpose_a1_intermediate_vv%open_('read', 'rewind')
      call wf%jacobian_transpose_a1_intermediate_vv%read_(Y_ca, wf%n_v**2)
      call wf%jacobian_transpose_a1_intermediate_vv%close_()
!
!     sigma_ai -= c_ci * Y_ca
!
      call dgemm('T', 'N',    & ! transpose of Y_ca
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
      call timer%turn_off()
!
   end subroutine jacobian_transpose_doubles_a1_doubles
!
!
  module subroutine jacobian_transpose_doubles_b1_doubles(wf, sigma_ai, c_bjck)
!!
!!    Jacobian transpose doubles B1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!    and Alexander C. Paul, Feb 2019
!!
!!    Calculates the B1 term,
!!
!!    sigma_ai =+ sum_bjc c_bjci g_bjca - c_akbj g_bjik
!!
      use batching_index_class, only: batching_index
!
      implicit none
!
      class(doubles) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_bjck
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)  :: sigma_ai
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ikbj
!
      real(dp), dimension(:,:,:), allocatable :: W_J_vo, L_J_vo, L_J_vv
!
      integer :: req0, req1, batch
!
      type(batching_index), allocatable :: batch_a
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian transpose doubles B1', pl='verbose')
      call timer%turn_on()
!
!     :: Term 2: sigma_ai =+ c_bjci g_bjca = (L_J_bj c_bjci) L_J_ca = L_J_ca W_J_ci
!
      call mem%alloc(L_J_vo, wf%eri%n_J, wf%n_v, wf%n_o)
      call wf%eri%get_cholesky_t1(L_J_vo, wf%n_o + 1, wf%n_mo, 1, wf%n_o)
!
      call mem%alloc(W_J_vo, wf%eri%n_J, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',          &
                  wf%eri%n_J,       &
                  wf%n_v * wf%n_o,  &
                  wf%n_v * wf%n_o,  &
                  one,              &
                  L_J_vo,           & ! L_J,bj
                  wf%eri%n_J,       &
                  c_bjck,           & ! c_bj,ci
                  wf%n_v * wf%n_o,  &
                  zero,             &
                  W_J_vo,           & ! W_J,ci
                  wf%eri%n_J)
!
      call mem%dealloc(L_J_vo, wf%eri%n_J, wf%n_v, wf%n_o)
!
      req0 = 0
      req1 = wf%n_v*wf%eri%n_J
!
      batch_a = batching_index(wf%n_v)
      call mem%batch_setup(batch_a, req0, req1, 'jacobian_transpose_doubles_b1')
!
      call mem%alloc(L_J_vv, wf%eri%n_J, wf%n_v, batch_a%max_length)
!
      do batch = 1, batch_a%num_batches
!
         call batch_a%determine_limits(batch)
!
         call wf%eri%get_cholesky_t1(L_J_vv, wf%n_o + 1, wf%n_mo, &
                                     wf%n_o + batch_a%first,     &
                                     wf%n_o + batch_a%get_last())
!
         call dgemm('T', 'N',                      &
                     batch_a%get_length(),         &
                     wf%n_o,                       &
                     wf%n_v * wf%eri%n_J,          &
                     one,                          &
                     L_J_vv,                       & ! L_Jc,a
                     wf%n_v * wf%eri%n_J,          &
                     W_J_vo,                       & ! W_Jc,i
                     wf%n_v * wf%eri%n_J,          &
                     one,                          &
                     sigma_ai(batch_a%first, 1),   &
                     wf%n_v)
!
      enddo
!
      call mem%dealloc(L_J_vv, wf%eri%n_J, wf%n_v, batch_a%max_length)
!
      call mem%batch_finalize()
!
      call mem%dealloc(W_J_vo, wf%eri%n_J, wf%n_v, wf%n_o)
!
!     :: Term 2: sigma_ai =+ sum_bjc c_akbj g_bjik = sum_bjc c_akbj (g_ikbj)^T
!
      call mem%alloc(g_ikbj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%eri%get_eri_t1('oovo', g_ikbj)
!
!     sigma_ai =- sum_bjk c_akbj g_ikbj
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
      call timer%turn_off()
!
   end subroutine jacobian_transpose_doubles_b1_doubles
!
!
  module subroutine jacobian_transpose_doubles_a2_doubles(wf, sigma_aibj, c_ai)
!!
!!    Jacobian transpose CC2 A2
!!    Written by Sarai D. Folkestad and Alexander C. Paul, Feb 2019
!!
!!    Calculates the A2 term,
!!
!!    sigma_aibj =+ (2F_jb c_ai - F_ib c_aj - L_ikjb c_ak + L_cajb c_ci)
!!
!!    Modified by Tor S. Haugland, Nov 2019
!!
!!    Now uses BLAS dger for outer-product instead of for-loops.
!!
      use batching_index_class, only: batching_index
      use array_utilities, only: zero_array, copy_and_scale
      use reordering, only: sort_12_to_21, add_1432_to_1234
      use reordering, only: add_2143_to_1234, add_4123_to_1234
!
      implicit none
!
      class(doubles) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                      :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout)   :: sigma_aibj
!
!     Local variables
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ikjb, g_cajb
      real(dp), dimension(:,:,:,:), allocatable :: L_kibj, L_cajb
      real(dp), dimension(:,:,:,:), allocatable :: sigma_ajbi
      real(dp), dimension(:,:,:,:), allocatable :: sigma_aibj_temp
      real(dp), dimension(:,:),     allocatable :: F_ai
!
      type(batching_index) :: batch_c
!
      integer :: req0, req1, current_c_batch
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian transpose doubles A2', pl='verbose')
      call timer%turn_on()
!
!     Term 1: (2F_jb c_ai - F_ib c_aj)
!
!     Sort F_ia to F_ai
!
      call mem%alloc(F_ai, wf%n_v, wf%n_o)
      call sort_12_to_21(wf%fock_ia, F_ai, wf%n_o, wf%n_v)
!
      call mem%alloc(sigma_aibj_temp, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call zero_array(sigma_aibj_temp, wf%n_o**2 * wf%n_v**2)
!
      call dger(wf%n_v * wf%n_o,    &
                wf%n_v * wf%n_o,    &
                one,                &
                c_ai,               & ! c_ai
                1,                  &
                F_ai,               & ! F_jb
                1,                  &
                sigma_aibj_temp,    & ! sigma_aibj
                wf%n_v * wf%n_o)
!
      call daxpy(wf%n_v**2 * wf%n_o**2, two, sigma_aibj_temp, 1, sigma_aibj, 1)
      call add_1432_to_1234(-one, sigma_aibj_temp, sigma_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(F_ai, wf%n_v, wf%n_o)
      call mem%dealloc(sigma_aibj_temp, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Term 3: - L_ikjb c_ak
!
!     L_ikjb = 2 g_ikjb - g_jkib (ordered as g_kibj)
!
      call mem%alloc(g_ikjb, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call wf%eri%get_eri_t1('ooov', g_ikjb)
!
      call mem%alloc(L_kibj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call zero_array(L_kibj, (wf%n_o**3)*wf%n_v)
!
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
      call mem%alloc(sigma_ajbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call zero_array(sigma_ajbi, wf%n_t1**2)
!
      req0 = (wf%n_v)*(wf%n_o)*(wf%eri%n_J)
      req1 = max((wf%n_v)*(wf%eri%n_J) + (wf%n_o)*(wf%n_v)**2, 2*(wf%n_o)*(wf%n_v)**2)
!
      batch_c = batching_index(wf%n_v)
!
      call mem%batch_setup(batch_c, req0, req1, 'jacobian_transpose_doubles_a2')
!
      do current_c_batch = 1, batch_c%num_batches
!
         call batch_c%determine_limits(current_c_batch)
!
!        L_cajb = 2 g_cajb - g_cbja
!
         call mem%alloc(g_cajb, batch_c%length, wf%n_v, wf%n_o, wf%n_v)
!
         call wf%eri%get_eri_t1('vvov', g_cajb, first_p=batch_c%first, last_p=batch_c%get_last())
!
         call mem%alloc(L_cajb, batch_c%length, wf%n_v, wf%n_o, wf%n_v)
!
         call copy_and_scale(two, g_cajb, L_cajb, (batch_c%length)*(wf%n_v**2)*(wf%n_o))
         call add_1432_to_1234(-one, g_cajb, L_cajb, batch_c%length, wf%n_v, wf%n_o, wf%n_v)
!
         call mem%dealloc(g_cajb, batch_c%length, wf%n_v, wf%n_o, wf%n_v)
!
         call dgemm('T', 'N',                &
                     (wf%n_v**2)*wf%n_o,     &
                     wf%n_o,                 &
                     batch_c%length,         &
                     one,                    &
                     L_cajb,                 &
                     batch_c%length,         &
                     c_ai(batch_c%first,1),  & ! c_ci
                     wf%n_v,                 &
                     one,                    &
                     sigma_ajbi,             &
                     (wf%n_v**2)*(wf%n_o))
!
         call mem%dealloc(L_cajb, batch_c%length, wf%n_v, wf%n_o, wf%n_v)
!
      enddo ! batch_c
!
      call mem%batch_finalize()
!
      call add_1432_to_1234(one, sigma_ajbi, sigma_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(sigma_ajbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_transpose_doubles_a2_doubles
!
!
end submodule jacobian_transpose_doubles
