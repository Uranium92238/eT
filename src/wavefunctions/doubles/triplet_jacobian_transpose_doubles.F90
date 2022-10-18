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
submodule (doubles_class) triplet_jacobian_transpose_doubles
!
!!
!!    Triplet Jacobian transpose submodule
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
   module subroutine triplet_jacobian_transpose_s_s_b_doubles(wf, sigma_ai, c_ai, t_aibj)
!!
!!
!!    Triplet Jacobian transpose singles-singles B
!!    Written by Sarai D. Folkestad, Feb 2022
!!
!!    Computes
!!
!!      sigma_ai += - L_kcla t_dl,ck c_di - L_kcid t_ck,dl c_al + g_icka t_dk,cl c_dl
!!
!
      use array_initialization, only: copy_and_scale
      use reordering, only: sort_1234_to_1432, add_1432_to_1234
!
      implicit none
!
      class(doubles), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)  :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: t_aibj
!
      real(dp), dimension(:,:,:), pointer :: L_Jvo_mo
!
      real(dp), dimension(:,:), allocatable :: X_ck, X_ad, X_li
!
      real(dp), dimension(:,:,:,:), allocatable :: g_vovo_MO, t_ckdl, g_aick, L_alck
!
!     Reorder t_cldk to t_ckdl
      call mem%alloc(t_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_1432(t_aibj, t_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(X_ck, wf%n_v, wf%n_o)
      call dgemm('N', 'N',       &
                  wf%n_v*wf%n_o, &
                  1,             &
                  wf%n_v*wf%n_o, &
                  one,           &
                  t_ckdl,        &
                  wf%n_v*wf%n_o, &
                  c_ai,          & ! c_dl
                  wf%n_v*wf%n_o, &
                  zero,          &
                  X_ck,          &
                  wf%n_v*wf%n_o)
      call mem%dealloc(t_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%L_mo%load_block(L_Jvo_mo, wf%n_o + 1, wf%n_mo, 1, wf%n_o)
!
      call mem%alloc(g_vovo_MO, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
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
      call wf%L_mo%offload_block(wf%n_o + 1, wf%n_mo, 1, wf%n_o)
!
!     Reorder g_akci to g_aick
      call mem%alloc(g_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_1432(g_vovo_MO, g_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N', &
                  wf%n_v*wf%n_o, &
                  1,             &
                  wf%n_v*wf%n_o, &
                  one,           &
                  g_aick,        &
                  wf%n_v*wf%n_o, &
                  X_ck,          &
                  wf%n_v*wf%n_o, &
                  one,           &
                  sigma_ai,        &
                  wf%n_o*wf%n_v)
      call mem%dealloc(g_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X_ck, wf%n_v, wf%n_o)
!
      call mem%alloc(L_alck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call copy_and_scale(two, g_vovo_MO, L_alck, wf%n_v**2*wf%n_o**2)
      call add_1432_to_1234(-one, g_vovo_MO, L_alck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_vovo_MO, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(X_ad, wf%n_v, wf%n_v)
      call dgemm('N', 'T',          &
                  wf%n_v,           &
                  wf%n_v,           &
                  wf%n_v*wf%n_o**2, &
                  one,              &
                  L_alck,           &
                  wf%n_v,           &
                  t_aibj,           & ! t_dlck
                  wf%n_v,           &
                  zero,             &
                  X_ad,             &
                  wf%n_v)
!
      call mem%alloc(X_li, wf%n_o, wf%n_o)
      call dgemm('T', 'N',          &
                  wf%n_o,           &
                  wf%n_o,           &
                  wf%n_o*wf%n_v**2, &
                  one,              &
                  t_aibj,           & ! t_ckdl
                  wf%n_o*wf%n_v**2, &
                  L_alck,           & ! L_ckdi
                  wf%n_o*wf%n_v**2, &
                  zero,             &
                  X_li,             &
                  wf%n_o)
!
      call mem%dealloc(L_alck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
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
                  sigma_ai,  &
                  wf%n_v)
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
                 sigma_ai,  &
                 wf%n_v)
!
      call mem%dealloc(X_ad, wf%n_v, wf%n_v)
      call mem%dealloc(X_li, wf%n_o, wf%n_o)
!
   end subroutine triplet_jacobian_transpose_s_s_b_doubles
!
!
   module subroutine triplet_jacobian_transpose_d_s_a_doubles(wf, sigma_ai, c_dlei)
!!
!!    Triplet Jacobian transpose doubles-singles A
!!    Written by Sarai D. Folkestad, Feb 2022
!!
!!    Computes
!!
!!       sigma_ai += c_eidl g_dlea - c_amdl g_imdl
!!
!
      use reordering, only: sort_123_to_132
      use batching_index_class, only: batching_index
!
      implicit none
!
      class(doubles), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)  :: c_dlei
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
!
      real(dp), dimension(:,:,:), pointer :: L_Jvo, L_Jea, L_Jim
      real(dp), dimension(:,:,:), allocatable :: X_Jvo, X_Jma, L_Jmi
!
      integer :: req_0, req_1
!
      integer :: current_a_batch
!
      type(batching_index) :: batch_a
!
      call wf%L_t1%load_block(L_Jvo, wf%n_o + 1, wf%n_mo, 1, wf%n_o)
!
      call mem%alloc(X_Jvo, wf%L_t1%n_J, wf%n_v, wf%n_o)
      call dgemm('N', 'T',       &
                  wf%L_t1%n_J,   &
                  wf%n_v*wf%n_o, &
                  wf%n_v*wf%n_o, &
                  one,           &
                  L_Jvo,         &
                  wf%L_t1%n_J,   &
                  c_dlei,        & ! c_ei_dl
                  wf%n_v*wf%n_o, &
                  zero,          &
                  X_Jvo,         &
                  wf%L_t1%n_J)
!
      call wf%L_t1%offload_block(wf%n_o + 1, wf%n_mo, 1, wf%n_o)
!
!     Compute c_eidl g_dlea
!
!
      req_0 = 0
      req_1 = wf%L_t1%load_memory_estimate(wf%n_o + 1,   &
                                          wf%n_mo,       &
                                          wf%n_o + 1,    &
                                          wf%n_o + 1)
!
      batch_a = batching_index(wf%n_v)
!
      call mem%batch_setup(batch_a, req_0, req_1, tag='triplet_jacobian_transpose_d_s_a_doubles')
!
      do current_a_batch = 1, batch_a%num_batches
!
         call batch_a%determine_limits(current_a_batch)
!
         call wf%L_t1%load_block(L_Jea, wf%n_o + 1, &
                                        wf%n_mo,    &
                                        wf%n_o + batch_a%first, &
                                        wf%n_o + batch_a%get_last())
         call dgemm('T', 'N',                  &
                     batch_a%length,           &
                     wf%n_o,                   &
                     wf%L_t1%n_J*wf%n_v,       &
                     one,                      &
                     L_Jea,                    &
                     wf%L_t1%n_J*wf%n_v,       &
                     X_Jvo,                    & ! X_Jei
                     wf%L_t1%n_J*wf%n_v,       &
                     one,                      &
                     sigma_ai(batch_a%first, 1), &
                     wf%n_v)
!
         call wf%L_t1%offload_block(wf%n_o + 1, &
                                    wf%n_mo,    &
                                    wf%n_o + batch_a%first, &
                                    wf%n_o + batch_a%get_last())
!
      enddo
!
      call mem%batch_finalize()
!
!     Compute - c_amdl g_imdl
!
      call mem%alloc(X_Jma, wf%L_t1%n_J, wf%n_o, wf%n_v)
      call sort_123_to_132(X_Jvo, X_Jma, wf%L_t1%n_J, wf%n_v, wf%n_o)
!
      call mem%dealloc(X_Jvo, wf%L_t1%n_J, wf%n_v, wf%n_o)
!
      call wf%L_t1%load_block(L_Jim, 1, wf%n_o, 1, wf%n_o)
!
      call mem%alloc(L_Jmi, wf%L_t1%n_J, wf%n_o, wf%n_o)
      call sort_123_to_132(L_Jim, L_Jmi, wf%L_t1%n_J, wf%n_o, wf%n_o)
      call wf%L_t1%offload_block(1, wf%n_o, 1, wf%n_o)

      call dgemm('T', 'N',             &
                  wf%n_v,              &
                  wf%n_o,              &
                  wf%L_t1%n_J*wf%n_o,  &
                  -one,                &
                  X_Jma,               &
                  wf%L_t1%n_J*wf%n_o,  &
                  L_Jmi,               &
                  wf%L_t1%n_J*wf%n_o,  &
                  one,                 &
                  sigma_ai,              &
                  wf%n_v)
!
      call mem%dealloc(X_Jma, wf%L_t1%n_J, wf%n_o, wf%n_v)
      call mem%dealloc(L_Jmi, wf%L_t1%n_J, wf%n_o, wf%n_o)
!
   end subroutine triplet_jacobian_transpose_d_s_a_doubles
!
!
   module subroutine triplet_jacobian_transpose_s_d_a_doubles(wf, sigma_aibj, c_ai)
!!
!!    Triplet Jacobian transpose singles-doubles A
!!    Written by Sarai D. Folkestad, Feb 2022
!!
!!    Computes
!!
!!       sigma_aibj += F_jb c_ai
!!
      implicit none
!
      class(doubles), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)  :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout) :: sigma_aibj
!
      integer :: a, i, b, j
!
!$omp parallel do private (a, i, b, j)
      do j = 1, wf%n_o
         do b = 1, wf%n_v
            do i = 1, wf%n_o
               do a = 1,wf%n_v
!
                  sigma_aibj(a, i, b, j) = sigma_aibj(a, i, b, j) + wf%fock_ia(j,b)*c_ai(a,i)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine triplet_jacobian_transpose_s_d_a_doubles
!
!
   module subroutine triplet_jacobian_transpose_s_d_b_doubles(wf, sigma_aibj, c_ai)
!!
!!    Triplet Jacobian transpose singles-doubles B
!!    Written by Sarai D. Folkestad, Feb 2022
!!
!!    Computes
!!
!!       sigma_aibj += -  g_jbil c_al
!!
!
      use reordering, only: add_4321_to_1234
!
      implicit none
!
      class(doubles), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)  :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout) :: sigma_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: g_jbil, g_jbia
!
      call mem%alloc(g_jbil, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call wf%eri_t1%get('ovoo', g_jbil)
!
      call mem%alloc(g_jbia, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call dgemm('N', 'T',          &
                  wf%n_o**2*wf%n_v, &
                  wf%n_v,           &
                  wf%n_o,           &
                  one,              &
                  g_jbil,           &
                  wf%n_o**2*wf%n_v, &
                  c_ai,             & ! c_al
                  wf%n_v,           &
                  zero,             &
                  g_jbia,           &
                  wf%n_o**2*wf%n_v)
!
      call mem%dealloc(g_jbil, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call add_4321_to_1234(-one, g_jbia, sigma_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_jbia, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
   end subroutine triplet_jacobian_transpose_s_d_b_doubles
!
!
   module subroutine triplet_jacobian_transpose_s_d_c_doubles(wf, sigma_aibj, c_ai)
!!
!!    Triplet Jacobian transpose singles-doubles C
!!    Written by Sarai D. Folkestad, Feb 2022
!!
!!    Computes
!!
!!       sigma_aibj +=  g_dajb c_di
!!
!
      use array_initialization, only: zero_array
      use reordering, only: add_4321_to_1234, sort_123_to_132
      use batching_index_class, only: batching_index
!
      implicit none
!
      class(doubles), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)  :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout) :: sigma_aibj
!
      real(dp), dimension(:,:,:), pointer :: L_Jda, L_Jbj
      real(dp), dimension(:,:,:), allocatable :: L_Jad, X_Jai
!
      integer :: req_0, req_1
!
      integer :: current_d_batch
!
      type(batching_index) :: batch_d
!
      call mem%alloc(X_Jai, wf%L_t1%n_J, wf%n_v, wf%n_o)
      call zero_array(X_Jai, wf%L_t1%n_J*wf%n_v*wf%n_o)
!
      req_0 = 0
      req_1 = wf%n_v*wf%L_t1%n_J &
            + wf%L_t1%load_memory_estimate(wf%n_o + 1, wf%n_o + 1, wf%n_o + 1, wf%n_mo)
!
      batch_d = batching_index(wf%n_v)
!
      call mem%batch_setup(batch_d, req_0, req_1, tag='triplet_jacobian_transpose_s_d_c_doubles')
!
      do current_d_batch = 1, batch_d%num_batches
!
         call batch_d%determine_limits(current_d_batch)
!
         call wf%L_t1%load_block(L_Jda, wf%n_o + batch_d%first, &
                                        wf%n_o + batch_d%get_last(), &
                                        wf%n_o + 1, wf%n_mo)
!
         call mem%alloc(L_Jad, wf%L_t1%n_J, wf%n_v, batch_d%length)
         call sort_123_to_132(L_Jda, L_Jad, wf%L_t1%n_J, batch_d%length, wf%n_v)
!
         call wf%L_t1%offload_block(wf%n_o + batch_d%first, &
                                    wf%n_o + batch_d%get_last(), &
                                    wf%n_o + 1, wf%n_mo)
!
         call dgemm('N', 'N',                &
                    wf%L_t1%n_J*wf%n_v,      &
                    wf%n_o,                  &
                    batch_d%length,          &
                    one,                     &
                    L_Jad,                   &
                    wf%L_t1%n_J*wf%n_v,      &
                    c_ai(batch_d%first, 1),  & ! c_di
                    wf%n_v,                  &
                    one,                     &
                    X_Jai,                   &
                    wf%L_t1%n_J*wf%n_v)
!
         call mem%dealloc(L_Jad, wf%L_t1%n_J, wf%n_v, batch_d%length)
!
      enddo
!
      call mem%batch_finalize()
!
!     using L_T1_jb = L_MO_jb = L_MO_bj
      call wf%L_mo%load_block(L_Jbj, wf%n_o + 1, wf%n_mo, 1, wf%n_o)
!
      call dgemm('T', 'N',       &
                  wf%n_v*wf%n_o, &
                  wf%n_v*wf%n_o, &
                  wf%L_t1%n_J,   &
                  one,           &
                  X_Jai,         &
                  wf%L_t1%n_J,   &
                  L_Jbj,         &
                  wf%L_t1%n_J,   &
                  one,           &
                  sigma_aibj,      &
                  wf%n_v*wf%n_o)
!
      call mem%dealloc(X_Jai, wf%L_t1%n_J, wf%n_v, wf%n_o)
      call wf%L_t1%offload_block(wf%n_o + 1, wf%n_mo, 1, wf%n_o)
!
   end subroutine triplet_jacobian_transpose_s_d_c_doubles
!
!
end submodule triplet_jacobian_transpose_doubles
