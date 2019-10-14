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
submodule (abstract_doubles_class) jacobian_transpose_abstract_doubles
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
   module subroutine jacobian_transpose_doubles_a1_abstract_doubles(wf, sigma_ai, c_bj, u)
!!
!!    Jacobian transpose doubles A1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad 
!!    and Alexander Paul, Feb 2019
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
      class(abstract_doubles) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)     :: c_bj
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)  :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: u
!
      real(dp), dimension(:,:,:,:), allocatable :: g_iakc
      real(dp), dimension(:,:,:,:), allocatable :: L_aick
      real(dp), dimension(:,:,:,:), allocatable :: u_cjbk
!
      real(dp), dimension(:,:), allocatable :: X_ck, Y_ik, Y_ca
!
      type(timings) :: timer
!
      timer = timings('jacobian transpose doubles a1')
      call timer%turn_on()
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
                  u,                   & ! u_ck_bj
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
      call zero_array(L_aick, wf%n_t1**2)
!
      call add_2143_to_1234(two, g_iakc, L_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call add_2341_to_1234(-one, g_iakc, L_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
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
!      Note: Pretend g_iakc is g_icjb
!
      call mem%alloc(u_cjbk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_3214(u, u_cjbk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
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
                  g_iakc,                 & ! g_i_cjb
                  wf%n_o,                 &
                  u_cjbk,                 & ! u_cjb_k
                  (wf%n_o)*(wf%n_v)**2,   &
                  zero,                   &
                  Y_ik,                   & ! Y_ik
                  wf%n_o)
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
!     Note: we pretend that g_iakc is g_jbka 
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
                  g_iakc,                 & ! g_jbk_a
                  (wf%n_v)*(wf%n_o)**2,   &
                  zero,                   &
                  Y_ca,                   & ! Y_ca
                  wf%n_v)
!
      call mem%dealloc(g_iakc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
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
      call timer%turn_off()
!
   end subroutine jacobian_transpose_doubles_a1_abstract_doubles
!
!
  module subroutine jacobian_transpose_doubles_b1_abstract_doubles(wf, sigma_ai, c_bjck)
!!
!!    Jacobian transpose doubles B1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad 
!!    and Alexander Paul, Feb 2019
!!
!!    Calculates the B1 term,
!!
!!    sigma_ai =+ sum_bjc c_bjci g_bjca - c_akbj g_bjik
!!
      implicit none
!
      class(abstract_doubles) :: wf
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
      type(timings) :: timer
!
      timer = timings('jacobian transpose doubles b1')
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
         call wf%get_vovv(g_bjca,                        &
                           1, wf%n_v,                    &
                           1, wf%n_o,                    &
                           1, wf%n_v,                    &
                           batch_a%first, batch_a%last)
!
!        sigma_ai =+ sum_bjc g_abjc * c_bjci
!
         call dgemm('T', 'N',                    & ! transposed g_bjca
                     batch_a%length,             &
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
   end subroutine jacobian_transpose_doubles_b1_abstract_doubles
!
!
  module subroutine jacobian_transpose_doubles_a2_abstract_doubles(wf, sigma_aibj, c_ai)
!!
!!    Jacobian transpose CC2 A2
!!    Written by Sarai D. Folkestad and Alexander Paul, Feb 2019
!!
!!    Calculates the A2 term,
!!
!!    sigma_aibj =+ (2F_jb c_ai - F_ib c_aj - L_ikjb c_ak + L_cajb c_ci)
!!
      implicit none
!
      class(abstract_doubles) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                      :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout)   :: sigma_aibj
!
!     Local variables
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ikjb, g_cajb
      real(dp), dimension(:,:,:,:), allocatable :: L_kibj, L_cajb
      real(dp), dimension(:,:,:,:), allocatable :: sigma_ajbi
!
      type(batching_index) :: batch_c
!
      integer :: req0, req1, current_c_batch
!
      integer :: a, i, b, j
!
      type(timings) :: timer
!
      timer = timings('jacobian transpose doubles a2')
      call timer%turn_on()
!
!     Term 1: (2F_jb c_ai - F_ib c_aj)
!
      call zero_array(sigma_aibj, (wf%n_o*wf%n_v)**2)
!
!$omp parallel do private(a, i, b, j)
      do j = 1, wf%n_o
         do b = 1, wf%n_v
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
      req0 = (wf%n_v)*(wf%n_o)*(wf%integrals%n_J)
      req1 = (wf%n_v)*(wf%integrals%n_J) + (wf%n_o)*(wf%n_v)**2
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
         call wf%get_vvov(g_cajb, &
                           batch_c%first, batch_c%last,   &
                           1, wf%n_v,                     &
                           1, wf%n_o,                     &
                           1, wf%n_v)
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
      call add_1432_to_1234(one, sigma_ajbi, sigma_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(sigma_ajbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_transpose_doubles_a2_abstract_doubles
!
end submodule jacobian_transpose_abstract_doubles
