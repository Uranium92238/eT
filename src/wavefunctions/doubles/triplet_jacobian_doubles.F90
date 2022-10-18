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
submodule (doubles_class) triplet_jacobian_doubles
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
!
   implicit none
!
!
contains
!
!
   module subroutine triplet_jacobian_s_s_b_doubles(wf, rho_ai, c_ai, t_aibj)
!!
!!    Triplet jacobian singles-singles B
!!    Written by Sarai D. Folkestad, Feb 2022
!!
!!    Computes
!!
!!      -rho_ai +=  L_kcld t_al,ck c_di - L_kcld t_ck,di c_al + g_lckd t_ak,ci c_dl
!!
!
      use array_initialization, only: copy_and_scale
      use reordering, only: sort_1234_to_1432, add_3214_to_1234
!
      implicit none
!
      class(doubles), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)  :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: t_aibj
!
      real(dp), dimension(:,:,:), pointer :: L_Jvo_mo
!
      real(dp), dimension(:,:,:,:), allocatable :: g_vovo_MO
      real(dp), dimension(:,:,:,:), allocatable :: L_ckdl, g_ckdl, t_aick
      real(dp), dimension(:,:), allocatable :: X_ck, X_ad, X_li
!
      call mem%alloc(g_vovo_MO, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%L_mo%load_block(L_Jvo_mo, wf%n_o + 1, wf%n_mo, 1, wf%n_o)
!
      call dgemm('T', 'N', &
                  wf%n_v*wf%n_o, &
                  wf%n_v*wf%n_o, &
                  wf%L_mo%n_J,   &
                  one,           &
                  L_Jvo_mo,      &
                  wf%L_mo%n_J,   &
                  L_Jvo_mo,      &
                  wf%L_mo%n_J,   &
                  zero,          &
                  g_vovo_MO,     &
                  wf%n_v*wf%n_o)
!
      call wf%L_mo%offload_block(wf%n_o + 1, wf%n_mo, 1, wf%n_o)
!
!     Compute g_lckd t_ak,ci c_dl
!
!     reorder g_cldk to g_ckdl
      call mem%alloc(g_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_1432(g_vovo_MO, g_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(X_ck, wf%n_v, wf%n_o)
      call dgemm('N', 'N',       &
                  wf%n_v*wf%n_o, &
                  1,             &
                  wf%n_v*wf%n_o, &
                  one,           &
                  g_ckdl,        &
                  wf%n_v*wf%n_o, &
                  c_ai,          & ! c_dl
                  wf%n_v*wf%n_o, &
                  zero,          &
                  X_ck,          &
                  wf%n_v*wf%n_o)
!
      call mem%dealloc(g_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     reorder t_akci to t_aick
      call mem%alloc(t_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_1432(t_aibj, t_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N', &
                  wf%n_v*wf%n_o, &
                  1,             &
                  wf%n_v*wf%n_o, &
                  one,           &
                  t_aick,        &
                  wf%n_v*wf%n_o, &
                  X_ck,          &
                  wf%n_v*wf%n_o, &
                  one,           &
                  rho_ai,        &
                  wf%n_v*wf%n_o)
!
      call mem%dealloc(t_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X_ck, wf%n_v, wf%n_o)
!
!     Compute - L_kcld t_al,ck c_di - L_kcld t_ck,di c_al
!
      call mem%alloc(L_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call copy_and_scale(two, g_vovo_MO, L_ckdl, wf%n_o**2*wf%n_v**2)
      call add_3214_to_1234(-one, g_vovo_MO, L_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_vovo_MO, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(X_ad, wf%n_v, wf%n_v)
!
      call dgemm('N', 'T',          &
                  wf%n_v,           &
                  wf%n_v,           &
                  wf%n_o**2*wf%n_v, &
                  one,              &
                  t_aibj,           & ! t_alck
                  wf%n_v,           &
                  L_ckdl,           & ! L_dlck
                  wf%n_v,           &
                  zero,             &
                  X_ad,             &
                  wf%n_v)
!
      call mem%alloc(X_li, wf%n_o, wf%n_o)
      call dgemm('T', 'N',          &
                  wf%n_o,           &
                  wf%n_o,           &
                  wf%n_v**2*wf%n_o, &
                  one,              &
                  L_ckdl,           &
                  wf%n_v**2*wf%n_o, &
                  t_aibj,           & ! t_ckdi
                  wf%n_v**2*wf%n_o, &
                  zero,             &
                  X_li,             &
                  wf%n_o)
!
      call mem%dealloc(L_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N', &
                 wf%n_v,  &
                 wf%n_o,  &
                 wf%n_v,  &
                 -one,    &
                 X_ad,    &
                 wf%n_v,  &
                 c_ai,    & ! c_di
                 wf%n_v,  &
                 one,     &
                 rho_ai,  &
                 wf%n_v)
!
!
      call dgemm('N', 'N', &
                  wf%n_v,  &
                  wf%n_o,  &
                  wf%n_o,  &
                  -one,    &
                  c_ai,    & ! c_al
                  wf%n_v,  &
                  X_li,    &
                  wf%n_o,  &
                  one,     &
                  rho_ai,  &
                  wf%n_v)
!
      call mem%dealloc(X_ad, wf%n_v, wf%n_v)
      call mem%dealloc(X_li, wf%n_o, wf%n_o)
!
   end subroutine triplet_jacobian_s_s_b_doubles
!
!
   module subroutine triplet_jacobian_s_d_a_doubles(wf, rho_ai, c_aibj)
!!
!!    Triplet Jacobian singles-doubles
!!    Written by Sarai D. Folkestad, Feb 2022
!!
!!    Computes
!!
!!       rho_ai += 2 F_kc c_aick
!!
!
      use reordering, only: sort_12_to_21
!
      implicit none
!
      class(doubles), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)  :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
!
      real(dp), dimension(:,:), allocatable :: F_ck
!
      call mem%alloc(F_ck, wf%n_v, wf%n_o)
      call sort_12_to_21(wf%fock_ia, F_ck, wf%n_o, wf%n_v)
!
      call dgemm('N', 'N',       &
                  wf%n_v*wf%n_o, &
                  1,             &
                  wf%n_v*wf%n_o, &
                  two,           &
                  c_aibj,        & ! c_aick
                  wf%n_v*wf%n_o, &
                  F_ck,          &
                  wf%n_v*wf%n_o, &
                  one,           &
                  rho_ai,        &
                  wf%n_v*wf%n_o)
!
      call mem%dealloc(F_ck, wf%n_v, wf%n_o)
!
   end subroutine triplet_jacobian_s_d_a_doubles
!
!
   module subroutine triplet_jacobian_s_d_b_doubles(wf, rho_ai, c_dick)
!!
!!
!!    Triplet Jacobian singles-doubles C
!!    Written by Sarai D. Folkestad, Feb 2022
!!
!!    Computes
!!
!!       rho_ai += 2 g_adkc c_dick
!!
!!    Note: we make no assuption on symmetries of C
!!
!
      use reordering, only: sort_1234_to_1432
      use batching_index_class, only: batching_index
!
      implicit none
!
      class(doubles), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)  :: c_dick
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
!
      real(dp), dimension(:,:,:,:), allocatable :: g_adkc, c_dkci
!
      integer :: req_0, req_1
      integer :: i, c, k, d
!
      integer :: current_d_batch
!
      type(batching_index) :: batch_d
!
      integer, dimension(2) :: memory
!
      memory = wf%eri_t1%get_memory_estimate('vvov', wf%n_v, 1, wf%n_o, wf%n_v)
!
      req_0 = memory(2)
      req_1 = wf%n_v*wf%n_o**2 + wf%n_v**2*wf%n_o + memory(1)
!
      batch_d = batching_index(wf%n_v)
!
      call mem%batch_setup(batch_d, req_0, req_1, tag='triplet_jacobian_s_d_b_doubles')
!
      do current_d_batch = 1, batch_d%num_batches
!
         call batch_d%determine_limits(current_d_batch)
!
         call mem%alloc(g_adkc, wf%n_v, batch_d%length, wf%n_o, wf%n_v)
!
         call wf%eri_t1%get('vvov', g_adkc,                    &
                            1, wf%n_v,                         &
                            batch_d%first, batch_d%get_last(), &
                            1, wf%n_o,                         &
                            1, wf%n_v)
!
         call mem%alloc(c_dkci, batch_d%length, wf%n_o, wf%n_v, wf%n_o)
!
!$omp parallel do private(i, c, k, d)
         do i = 1, wf%n_o
            do c = 1, wf%n_v
               do k = 1, wf%n_o
                  do d = 1, batch_d%length
!
                     c_dkci(d, k, c, i) = c_dick(d + batch_d%first - 1, i, c, k)
!
                  enddo
               enddo
            enddo
         enddo
!$omp end parallel do
!
         call dgemm('N', 'N',                      &
                     wf%n_v,                       &
                     wf%n_o,                       &
                     wf%n_v*wf%n_o*batch_d%length, &
                     two,                          &
                     g_adkc,                       &
                     wf%n_v,                       &
                     c_dkci,                       &
                     wf%n_v*wf%n_o*batch_d%length, &
                     one,                          &
                     rho_ai,                       &
                     wf%n_v)
!
         call mem%dealloc(c_dkci, batch_d%length, wf%n_o, wf%n_v, wf%n_o)
         call mem%dealloc(g_adkc, wf%n_v, batch_d%length, wf%n_o, wf%n_v)
!
      enddo
!
      call mem%batch_finalize()
!
   end subroutine triplet_jacobian_s_d_b_doubles
!
!
   module subroutine triplet_jacobian_s_d_c_doubles(wf, rho_ai, c_alck)
!!
!!
!!    Triplet Jacobian singles-doubles C
!!    Written by Sarai D. Folkestad, Feb 2022
!!
!!    Computes
!!
!!       rho_ai += - 2 g_kcli c_alck
!!
!!    Note: we make no assuption on symmetries of C
!!
!
      use reordering, only: sort_1234_to_3214
!
      implicit none
!
      class(doubles), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)    :: c_alck
      real(dp), dimension(wf%n_v, wf%n_o),                 intent(inout) :: rho_ai
!
      real(dp), dimension(:,:,:,:), allocatable :: g_kcli, g_lcki
!
      call mem%alloc(g_kcli, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call wf%eri_t1%get('ovoo', g_kcli)
!
      call mem%alloc(g_lcki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call sort_1234_to_3214(g_kcli, g_lcki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(g_kcli, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N',          &
                  wf%n_v,           &
                  wf%n_o,           &
                  wf%n_o**2*wf%n_v, &
                  -two,             &
                  c_alck,           &
                  wf%n_v,           &
                  g_lcki,           &
                  wf%n_o**2*wf%n_v, &
                  one,              &
                  rho_ai,           &
                  wf%n_v)
!
      call mem%dealloc(g_lcki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
   end subroutine triplet_jacobian_s_d_c_doubles
!
!
   module subroutine triplet_jacobian_d_s_a_doubles(wf, rho_aibj, c_ai)
!!
!!    Triplet Jacobian doubles-singles A
!!    Written by Sarai D. Folkestad, Feb 2022
!!
!!    Computes
!!
!!       rho_aibj += 1/2 g_aibc c_cj - 1/2 g_bjki c_ak
!!
      implicit none
!
      class(doubles), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o),                  intent(in)    :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o),  intent(inout) :: rho_aibj
!
      real(dp), dimension(:,:,:), pointer :: L_Jvo
!
      call wf%L_t1%load_block(L_Jvo, wf%n_o + 1, wf%n_mo, 1, wf%n_o)
!
      call wf%triplet_jacobian_d_s_a_1(rho_aibj, c_ai, L_Jvo)
      call wf%triplet_jacobian_d_s_a_2(rho_aibj, c_ai, L_Jvo)
!
      call wf%L_t1%offload_block(wf%n_o + 1, wf%n_mo, 1, wf%n_o)
!
   end subroutine triplet_jacobian_d_s_a_doubles
!
!
   module subroutine triplet_jacobian_d_s_a_1_doubles(wf, rho_aibj, c_ci, L_Jbj)
!!
!!    Triplet Jacobian doubles-singles A1
!!    Written by Sarai D. Folkestad, Feb 2022
!!
!!    Computes
!!
!!       rho_aibj += 1/2 g_bjac c_ci
!
!
      use batching_index_class, only: batching_index
      use array_initialization, only: zero_array
!
      implicit none
!
      class(doubles), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o),                  intent(in)     :: c_ci
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o),  intent(inout)  :: rho_aibj
      real(dp), dimension(wf%L_t1%n_J, wf%n_v, wf%n_o),     intent(in)     :: L_Jbj
!
      real(dp), dimension(:,:,:), pointer :: L_Jac
      real(dp), dimension(:,:,:), allocatable :: L_Jai
!
      integer :: req_0, req_1
!
      integer :: current_c_batch
!
      type(batching_index) :: batch_c
!
      call mem%alloc(L_Jai, wf%L_t1%n_J, wf%n_v, wf%n_o)
      call zero_array(L_Jai, wf%L_t1%n_J*wf%n_v*wf%n_o)
!
      req_0 = 0
      req_1 = wf%L_t1%load_memory_estimate(wf%n_o + 1, wf%n_mo, wf%n_o + 1, wf%n_o + 1)
!
      batch_c = batching_index(wf%n_v)
!
      call mem%batch_setup(batch_c, req_0, req_1, tag='triplet_jacobian_d_s_a_1_doubles')
!
      do current_c_batch = 1, batch_c%num_batches
!
         call batch_c%determine_limits(current_c_batch)

         call wf%L_t1%load_block(L_Jac,                  &
                                 wf%n_o + 1, wf%n_mo,    &
                                 wf%n_o + batch_c%first, &
                                 wf%n_o + batch_c%get_last())
!
         call dgemm('N', 'N',                      &
                     wf%L_t1%n_J*wf%n_v,           &
                     wf%n_o,                       &
                     batch_c%length,               &
                     one,                          &
                     L_Jac,                        &
                     wf%L_t1%n_J*wf%n_v,           &
                     c_ci(batch_c%first, 1),       &
                     wf%n_v,                       &
                     one,                          &
                     L_Jai,                        &
                     wf%L_t1%n_J*wf%n_v)
!
         call wf%L_t1%offload_block(wf%n_o + 1, wf%n_mo,    &
                                 wf%n_o + batch_c%first,    &
                                 wf%n_o + batch_c%get_last())
      enddo
!
      call mem%batch_finalize()
!
      call dgemm('T', 'N',            &
                  wf%n_o*wf%n_v,      &
                  wf%n_o*wf%n_v,      &
                  wf%L_t1%n_J,        &
                  half,               &
                  L_Jai,              &
                  wf%L_t1%n_J,        &
                  L_Jbj,              &
                  wf%L_t1%n_J,        &
                  one,                &
                  rho_aibj,           &
                  wf%n_o*wf%n_v)
!
      call mem%dealloc(L_Jai, wf%L_t1%n_J, wf%n_v, wf%n_o)
!
   end subroutine triplet_jacobian_d_s_a_1_doubles
!
!
   module subroutine triplet_jacobian_d_s_a_2_doubles(wf, rho_aibj, c_ak, L_Jbj)
!!
!!    Triplet Jacobian doubles-singles A2
!!    Written by Sarai D. Folkestad, Feb 2022
!!
!!    Computes
!!
!!       rho_aibj -=-1/2 g_bjki c_ak
!!
!
      use reordering, only: sort_123_to_132
!
      implicit none
!
      class(doubles), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ak
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout) :: rho_aibj
      real(dp), dimension(wf%L_t1%n_J, wf%n_v, wf%n_o), intent(in) :: L_Jbj
!
      real(dp), dimension(:,:,:), pointer :: L_Jki
      real(dp), dimension(:,:,:), allocatable :: L_Jai, L_Jik, L_Jia
!
      call wf%L_t1%load_block(L_Jki, 1, wf%n_o, 1, wf%n_o)
!
      call mem%alloc(L_Jik, wf%L_t1%n_J, wf%n_o, wf%n_o)
      call sort_123_to_132(L_Jki, L_Jik, wf%L_t1%n_J, wf%n_o, wf%n_o)
!
      call wf%L_t1%offload_block(1, wf%n_o, 1, wf%n_o)
!
      call mem%alloc(L_Jia, wf%L_t1%n_J, wf%n_o, wf%n_v)
      call dgemm('N', 'T',            &
                  wf%L_t1%n_J*wf%n_o, &
                  wf%n_v,             &
                  wf%n_o,             &
                  one,                &
                  L_Jik,              &
                  wf%L_t1%n_J*wf%n_o, &
                  c_ak,               &
                  wf%n_v,             &
                  zero,               &
                  L_Jia,              &
                  wf%L_t1%n_J*wf%n_o)
!
      call mem%dealloc(L_Jik, wf%L_t1%n_J, wf%n_o, wf%n_o)
!
      call mem%alloc(L_Jai, wf%L_t1%n_J, wf%n_v, wf%n_o)
      call sort_123_to_132(L_Jia, L_Jai, wf%L_t1%n_J, wf%n_o, wf%n_v)
      call mem%dealloc(L_Jia, wf%L_t1%n_J, wf%n_o, wf%n_v)
!
      call dgemm('T', 'N',            &
                  wf%n_o*wf%n_v,      &
                  wf%n_o*wf%n_v,      &
                  wf%L_t1%n_J,        &
                  -half,              &
                  L_Jai,              &
                  wf%L_t1%n_J,        &
                  L_Jbj,              &
                  wf%L_t1%n_J,        &
                  one,                &
                  rho_aibj,           &
                  wf%n_o*wf%n_v)
!
      call mem%dealloc(L_Jai, wf%L_t1%n_J, wf%n_v, wf%n_o)
!
   end subroutine triplet_jacobian_d_s_a_2_doubles
!
!
   module subroutine packout_triplet_d(wf, x, x_aibj_p, x_aibj_m)
!!
!!    Packout triplet doubles
!!    Written by Sarai D. Folkestad, Feb 2022
!!
!!    Doubles part of the triplet vector X^T = (Y, Z)^T contains
!!
!!       Y_aibj for a > b and i > j
!!       Z_aibj for ai > bj
!!
!!
      use array_initialization, only: zero_array
      use reordering, only: squareup_anti
      implicit none
!
      class(doubles), intent(in) :: wf
!
      real(dp), dimension(wf%n_triplet_amplitudes), intent(in) :: x
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(out) :: x_aibj_p, x_aibj_m
!
      integer :: i, j, a, b, ab, ij, abij, n
!
      call zero_array(x_aibj_p, wf%n_v**2*wf%n_o**2)
!
!$omp parallel do private(a, i, b, j, ab, ij, abij)
      do a = 2, wf%n_v
         do b = 1, a - 1
!
            ab = (max(a-1,b)*(max(a-1,b)-3)/2) + a-1 + b
!
            do i = 2, wf%n_o
               do j = 1, i - 1
!
                  ij = (max(i-1,j)*(max(i-1,j)-3)/2) + i-1 + j
                  abij = (wf%n_v*(wf%n_v-1)/2)*(ij - 1) + ab
!
                  x_aibj_p(a, i, b, j) = x(wf%n_t1 + abij)
                  x_aibj_p(b, j, a, i) = x(wf%n_t1 + abij)
                  x_aibj_p(b, i, a, j) = -x(wf%n_t1 + abij)
                  x_aibj_p(a, j, b, i) = -x(wf%n_t1 + abij)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      n = wf%n_t1 + 1 +(wf%n_v*(wf%n_v-1)/2)*(wf%n_o*(wf%n_o-1)/2)
      call squareup_anti(x(n:wf%n_triplet_amplitudes), x_aibj_m, wf%n_o*wf%n_v)
!
   end subroutine packout_triplet_d
!
!
   module subroutine packin_triplet_d(wf, x, x_aibj_p, x_aibj_m)
!!
!!    Packin triplet doubles
!!    Written by Sarai D. Folkestad, Feb 2022
!!
!!    Doubles part of the triplet vector X^T = (Y, Z)^T contains
!!
!!       Y_aibj for a > b and i > j
!!       Z_aibj for ai > bj
!!
!!
      use reordering, only: packin_anti
      implicit none
!
      class(doubles), intent(in) :: wf
!
      real(dp), dimension(wf%n_triplet_amplitudes), intent(inout) :: x
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: x_aibj_p, x_aibj_m
!
      integer :: i, j, a, b, ab, ij, abij, n
!
!$omp parallel do private(a, i, b, j, ab, ij, abij)
      do a = 2, wf%n_v
         do b = 1, a - 1
!
            ab = (max(a-1,b)*(max(a-1,b)-3)/2) + a-1 + b
!
            do i = 2, wf%n_o
               do j = 1, i - 1
!
                  ij = (max(i-1,j)*(max(i-1,j)-3)/2) + i-1 + j
                  abij = wf%n_v*(wf%n_v-1)/2*(ij - 1) + ab
!
                  x(wf%n_t1 + abij) = x_aibj_p(a, i, b, j)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      n = wf%n_t1 + (wf%n_v*(wf%n_v-1)/2)*(wf%n_o*(wf%n_o-1)/2) + 1

      call packin_anti(x(n:), x_aibj_m, wf%n_o*wf%n_v)
!
   end subroutine packin_triplet_d
!
!
end submodule triplet_jacobian_doubles
