!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2020 the authors of eT
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
submodule (doubles_class) jacobian_transpose_doubles_complex
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
   module subroutine save_jacobian_transpose_a1_intermediates_doubles_complex(wf, u_bjck)
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
!!    jacobian_transpose_doubles_a1_doubles_complex and wrote them to file. 
!!
      class(doubles), intent(inout) :: wf
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: u_bjck
!
      complex(dp), dimension(:,:,:,:), allocatable :: g_ovov
      complex(dp), dimension(:,:,:,:), allocatable :: u_cjbk
!
      complex(dp), dimension(:,:),     allocatable :: Y_ik
      complex(dp), dimension(:,:),     allocatable :: Y_ca
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian transpose A1 intermediate', pl='verbose')
      call timer%turn_on()
!
!     Both intermediates need g_ovov and u_bjck -> u_cjbk
!
      call mem%alloc(g_ovov, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call wf%get_ovov_complex(g_ovov)
!
      call mem%alloc(u_cjbk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_3214(u_bjck, u_cjbk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!        Y_ik = sum_cjb g_i_cjb * u_cjb_k
!
      call mem%alloc(Y_ik, wf%n_o, wf%n_o)
!
      call zgemm('N', 'N',                &
                  wf%n_o,                 &
                  wf%n_o,                 &
                  (wf%n_o)*(wf%n_v)**2,   &
                  one_complex,                    &
                  g_ovov,                 & ! g_i_cjb
                  wf%n_o,                 &
                  u_cjbk,                 & ! u_cjb_k
                  (wf%n_o)*(wf%n_v)**2,   &
                  zero_complex,                   &
                  Y_ik,                   & ! Y_ik
                  wf%n_o)
!
!     Save Y_ik
!
      wf%jacobian_transpose_a1_intermediate_oo = &
                                          sequential_file('jacobian_transpose_intermediate_a1_oo')
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
      call zgemm('N', 'N',                &
                  wf%n_v,                 &
                  wf%n_v,                 &
                  (wf%n_v)*(wf%n_o)**2,   &
                  one_complex,                    &
                  u_cjbk,                 & ! u_c_jbk
                  wf%n_v,                 &
                  g_ovov,                 & ! g_jbk_a
                  (wf%n_v)*(wf%n_o)**2,   &
                  zero_complex,                   &
                  Y_ca,                   & ! Y_ca
                  wf%n_v)
!
!     Save Y_ca
!
      wf%jacobian_transpose_a1_intermediate_vv = &
                                          sequential_file('jacobian_transpose_intermediate_a1_vv')
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
   end subroutine save_jacobian_transpose_a1_intermediates_doubles_complex
!
!
   module subroutine jacobian_transpose_doubles_a1_doubles_complex(wf, sigma_ai, c_bj, u)
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
      implicit none
!
      class(doubles) :: wf
!
      complex(dp), dimension(wf%n_v, wf%n_o), intent(in)     :: c_bj
      complex(dp), dimension(wf%n_v, wf%n_o), intent(inout)  :: sigma_ai
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: u
!
      complex(dp), dimension(:,:,:,:), allocatable :: g_iakc
      complex(dp), dimension(:,:,:,:), allocatable :: L_iakc
!
      complex(dp), dimension(:,:), allocatable :: X_ck, X_kc, Y_ik, Y_ca, sigma_ia, sigma_ai_temp
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
      call zgemv('N',                  &
                  (wf%n_o)*(wf%n_v),   &
                  (wf%n_o)*(wf%n_v),   &
                  one_complex,                 &
                  u,                   & ! u_ck_bj
                  (wf%n_o)*(wf%n_v),   &
                  c_bj,                & ! c_bj
                  1,                   &
                  zero_complex,                &
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
      call wf%get_ovov_complex(g_iakc)
!
      call copy_and_scale_complex(two_complex, g_iakc, L_iakc, wf%n_t1**2)
      call add_1432_to_1234(-one_complex, g_iakc, L_iakc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%dealloc(g_iakc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     sigma_ai += sum_kc L_iakc * X_kc
!
      call mem%alloc(sigma_ia, wf%n_o, wf%n_v)
!
      call zgemm('N', 'N',             &
                  (wf%n_v)*(wf%n_o),   &
                  1,                   &
                  (wf%n_v)*(wf%n_o),   &
                  one_complex,                 &
                  L_iakc,              & ! L_ia_kc
                  (wf%n_v)*(wf%n_o),   &
                  X_kc,                & ! X_kc
                  (wf%n_v)*(wf%n_o),   &
                  zero_complex,                &
                  sigma_ia,            &
                  (wf%n_v)*(wf%n_o))
!
      call mem%alloc(sigma_ai_temp, wf%n_v, wf%n_o)
      call sort_12_to_21(sigma_ia, sigma_ai_temp, wf%n_o, wf%n_v)
!
      call zaxpy(wf%n_v * wf%n_o, one_complex, sigma_ai_temp, 1, sigma_ai, 1)
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
      call zgemm('N', 'T',    & ! transpose of Y_ik
                  wf%n_v,     &
                  wf%n_o,     &
                  wf%n_o,     &
                  -one_complex,       &
                  c_bj,       & ! c_a_k
                  wf%n_v,     &
                  Y_ik,       & ! Y_k_i
                  wf%n_o,     &
                  one_complex,        &
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
      call zgemm('T', 'N',    & ! transpose of Y_ca
                  wf%n_v,     &
                  wf%n_o,     &
                  wf%n_v,     &
                  -one_complex,       &
                  Y_ca,       & ! Y_a_c
                  wf%n_v,     &
                  c_bj,       & ! c_c_i
                  wf%n_v,     &
                  one_complex,        &
                  sigma_ai,   &
                  wf%n_v)
!
      call mem%dealloc(Y_ca, wf%n_v, wf%n_v)
!
      call timer%turn_off()
!
   end subroutine jacobian_transpose_doubles_a1_doubles_complex
!
!
  module subroutine jacobian_transpose_doubles_b1_doubles_complex(wf, sigma_ai, c_bjck)
!!
!!    Jacobian transpose doubles B1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad 
!!    and Alexander C. Paul, Feb 2019
!!
!!    Calculates the B1 term,
!!
!!    sigma_ai =+ sum_bjc c_bjci g_bjca - c_akbj g_bjik
!!
      implicit none
!
      class(doubles) :: wf
!
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_bjck
      complex(dp), dimension(wf%n_v, wf%n_o), intent(inout)  :: sigma_ai
!
      complex(dp), dimension(:,:,:,:), allocatable :: g_bjca, g_ikbj
!
      type(batching_index) :: batch_a
!
      integer :: req0, req1, current_a_batch
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian transpose doubles B1', pl='verbose')
      call timer%turn_on()
!
!     :: Term 1: sigma_ai =+ sum_bjc c_bjci g_bjca = sum_bjc (g_bjca)^T c_bjci
!
      req0 = (wf%n_v)*(wf%n_o)*(wf%integrals%n_J)
      req1 = max((wf%n_v)*(wf%integrals%n_J) + (wf%n_o)*(wf%n_v)**2, 2*(wf%n_o)*(wf%n_v)**2)
!
      batch_a = batching_index(wf%n_v)
!
      call mem%batch_setup(batch_a, req0, req1)
!
      do current_a_batch = 1, batch_a%num_batches
!
         call batch_a%determine_limits(current_a_batch)
!
         call mem%alloc(g_bjca, wf%n_v, wf%n_o, wf%n_v, batch_a%length)
!
         call wf%get_vovv_complex(g_bjca,                        &
                           1, wf%n_v,                    &
                           1, wf%n_o,                    &
                           1, wf%n_v,                    &
                           batch_a%first, batch_a%last)
!
!        sigma_ai =+ sum_bjc g_abjc * c_bjci
!
         call zgemm('T', 'N',                    & ! transposed g_bjca
                     batch_a%length,             &
                     wf%n_o,                     &
                     (wf%n_o)*(wf%n_v)**2,       &
                     one_complex,                        &
                     g_bjca,                     & ! g_a_bjc
                     (wf%n_o)*(wf%n_v)**2,       &
                     c_bjck,                     & ! c_bjc_i
                     (wf%n_o)*(wf%n_v)**2,       &
                     one_complex,                        &
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
      call wf%get_oovo_complex(g_ikbj)
!
!     sigma_ai =- sum_bjk c_akbj g_ikbj
!
      call zgemm('N', 'T',                & ! transposed g_ikbj
                  wf%n_v,                 &
                  wf%n_o,                 &
                  (wf%n_v)*(wf%n_o)**2,   &
                  -one_complex,                   &
                  c_bjck,                 & ! c_a_kbj
                  wf%n_v,                 &
                  g_ikbj,                 & ! g_kbj_i
                  wf%n_o,                 &
                  one_complex,                    &
                  sigma_ai,               &
                  wf%n_v)
!
      call mem%dealloc(g_ikbj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_transpose_doubles_b1_doubles_complex
!
!
  module subroutine jacobian_transpose_doubles_a2_doubles_complex(wf, sigma_aibj, c_ai)
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
!!    Now uses BLAS zgeru for outer-product instead of for-loops.
!!
      implicit none
!
      class(doubles) :: wf
!
      complex(dp), dimension(wf%n_v, wf%n_o), intent(in)                      :: c_ai
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout)   :: sigma_aibj
!
!     Local variables
!
      complex(dp), dimension(:,:,:,:), allocatable :: g_ikjb, g_cajb
      complex(dp), dimension(:,:,:,:), allocatable :: L_kibj, L_cajb
      complex(dp), dimension(:,:,:,:), allocatable :: sigma_ajbi
      complex(dp), dimension(:,:,:,:), allocatable :: sigma_aibj_temp
      complex(dp), dimension(:,:),     allocatable :: F_ai
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
      call sort_12_to_21(wf%fock_ia_complex, F_ai, wf%n_o, wf%n_v)
!
      call mem%alloc(sigma_aibj_temp, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call zero_array_complex(sigma_aibj_temp, wf%n_o**2 * wf%n_v**2)
!
      call zgeru(wf%n_v * wf%n_o,    &
                wf%n_v * wf%n_o,    &
                one_complex,                &
                c_ai,               & ! c_ai
                1,                  &
                F_ai,               & ! F_jb
                1,                  &
                sigma_aibj_temp,    & ! sigma_aibj
                wf%n_v * wf%n_o)
!
      call zaxpy(wf%n_v**2 * wf%n_o**2, two_complex, sigma_aibj_temp, 1, sigma_aibj, 1)
      call add_1432_to_1234(-one_complex, sigma_aibj_temp, sigma_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(F_ai, wf%n_v, wf%n_o)
      call mem%dealloc(sigma_aibj_temp, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Term 3: - L_ikjb c_ak
!
!     L_ikjb = 2 g_ikjb - g_jkib (ordered as g_kibj)
!
      call mem%alloc(g_ikjb, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call wf%get_ooov_complex(g_ikjb)
!
      call mem%alloc(L_kibj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call zero_array_complex(L_kibj, (wf%n_o**3)*wf%n_v)
!
      call add_2143_to_1234(two_complex, g_ikjb, L_kibj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call add_4123_to_1234(-one_complex, g_ikjb, L_kibj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_ikjb, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call zgemm('N', 'N',             &
                  wf%n_v,              &
                  (wf%n_o**2)*wf%n_v,  &
                  wf%n_o,              &
                  -one_complex,                &
                  c_ai,                & ! c_a_k
                  wf%n_v,              &
                  L_kibj,              & ! L_k_ibj
                  wf%n_o,              &
                  one_complex,                 &
                  sigma_aibj,          & ! sigma_a_ibj
                  wf%n_v)
!
      call mem%dealloc(L_kibj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     Term 4: L_cajb c_ci
!
      call mem%alloc(sigma_ajbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call zero_array_complex(sigma_ajbi, wf%n_t1**2)
!
      req0 = (wf%n_v)*(wf%n_o)*(wf%integrals%n_J)
      req1 = max((wf%n_v)*(wf%integrals%n_J) + (wf%n_o)*(wf%n_v)**2, 2*(wf%n_o)*(wf%n_v)**2)
!
      batch_c = batching_index(wf%n_v)
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
         call wf%get_vvov_complex(g_cajb, &
                           batch_c%first, batch_c%last,   &
                           1, wf%n_v,                     &
                           1, wf%n_o,                     &
                           1, wf%n_v)
!
         call mem%alloc(L_cajb, batch_c%length, wf%n_v, wf%n_o, wf%n_v)
!
         call copy_and_scale_complex(two_complex, g_cajb, L_cajb, (batch_c%length)*(wf%n_v**2)*(wf%n_o))
         call add_1432_to_1234(-one_complex, g_cajb, L_cajb, batch_c%length, wf%n_v, wf%n_o, wf%n_v)
!
         call mem%dealloc(g_cajb, batch_c%length, wf%n_v, wf%n_o, wf%n_v)
!
         call zgemm('T', 'N',                &
                     (wf%n_v**2)*wf%n_o,     &
                     wf%n_o,                 &
                     batch_c%length,         &
                     one_complex,                    &
                     L_cajb,                 &
                     batch_c%length,         &
                     c_ai(batch_c%first,1),  & ! c_ci
                     wf%n_v,                 &
                     one_complex,                    &
                     sigma_ajbi,             &
                     (wf%n_v**2)*(wf%n_o))
!
         call mem%dealloc(L_cajb, batch_c%length, wf%n_v, wf%n_o, wf%n_v)
!
      enddo ! batch_c
!
      call add_1432_to_1234(one_complex, sigma_ajbi, sigma_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(sigma_ajbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_transpose_doubles_a2_doubles_complex
!
!
end submodule jacobian_transpose_doubles_complex
