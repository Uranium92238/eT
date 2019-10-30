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
submodule (mlcc2_class) jacobian_transpose_mlcc2
!
!!
!!    Jacobian transpose submodule (MLCC2)
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
!!    Note that all routines are adapted from jacobian_transpose_cc2.F90
!!    written by Sarai D. Folkestad and Alexander Paul, Feb 2019
!!
!
   implicit none
!
!
contains
!
!
   module subroutine prepare_for_jacobian_transpose_mlcc2(wf)
!!
!!    Jacobian transpose submodule (MLCC2)
!!    Written by Sarai D. Folkestad and Alexander Paul, Feb 2019
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      call wf%initialize_u()
      call wf%construct_u()
!
   end subroutine prepare_for_jacobian_transpose_mlcc2
!
!
   module subroutine jacobian_transpose_transformation_mlcc2(wf, b)
!!
!!    Jacobian transpose transformation (MLCC2)
!!    Adapted by Sarai D. Folkestad, Jul 2019
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
      class(mlcc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: b
!
      real(dp), dimension(:,:), allocatable :: b_ai
      real(dp), dimension(:,:,:,:), allocatable :: b_aibj
!
      real(dp), dimension(:,:), allocatable :: sigma_ai
      real(dp), dimension(:,:,:,:), allocatable :: sigma_aibj
!
!     Allocate and zero the transformed vecotr (singles part)
!
      call mem%alloc(sigma_ai, wf%n_v, wf%n_o)
      call zero_array(sigma_ai, wf%n_t1)
!
      call mem%alloc(b_ai, wf%n_v, wf%n_o)
!
      call dcopy(wf%n_t1, b, 1, b_ai, 1)
!
!     CCS contributions to the singles c vector
!
      call wf%jacobian_transpose_ccs_a1(sigma_ai, b_ai)
      call wf%jacobian_transpose_ccs_b1(sigma_ai, b_ai)
!
!     CC2 contributions to the transformed singles vector
!
      call wf%jacobian_transpose_cc2_a1(sigma_ai, b_ai, wf%n_cc2_o, wf%n_cc2_v, &
                                 wf%first_cc2_o, wf%first_cc2_v)
!
!     Allocate the incoming unpacked doubles vector
!
      call mem%alloc(b_aibj, (wf%n_cc2_v), (wf%n_cc2_o), (wf%n_cc2_v), (wf%n_cc2_o))
!
      call squareup(b(wf%n_t1 + 1 : wf%n_es_amplitudes), &
                    b_aibj, wf%n_cc2_v*wf%n_cc2_o)
!
      call wf%jacobian_transpose_cc2_b1(sigma_ai, b_aibj, wf%n_cc2_o, wf%n_cc2_v, &
                                       wf%first_cc2_o, wf%first_cc2_v, wf%last_cc2_o, wf%last_cc2_v)
!
!     Done with singles vector c; overwrite it with
!     transformed vector for exit
!
      call dcopy(wf%n_t1, sigma_ai, 1, b, 1)
!
      call mem%dealloc(sigma_ai, wf%n_v, wf%n_o)
!
!     CC2 contributions to the transformed doubles vector
!
!     Allocate unpacked transformed vector
!
      call mem%alloc(sigma_aibj, wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o)
      sigma_aibj = zero
!
!     Contributions from singles vector c
!
     call wf%jacobian_transpose_cc2_a2(sigma_aibj, b_ai, wf%n_cc2_o, wf%n_cc2_v, &
                                 wf%first_cc2_o, wf%first_cc2_v, wf%last_cc2_o, wf%last_cc2_v)

!     Symmetrize
!
      call symmetric_sum(sigma_aibj, (wf%n_cc2_o)*(wf%n_cc2_v))
!
      call mem%dealloc(b_ai, wf%n_v, wf%n_o)
!
!     Contributions from doubles vector c
!
      call wf%jacobian_transpose_cc2_b2(sigma_aibj, b_aibj)
!
      call mem%dealloc(b_aibj, (wf%n_cc2_v), (wf%n_cc2_o), (wf%n_cc2_v), (wf%n_cc2_o))
!
!     Overwrite the incoming doubles c vector & pack in
!
      call packin(b(wf%n_t1 + 1 : wf%n_es_amplitudes), &
                  sigma_aibj, wf%n_cc2_v*wf%n_cc2_o)
!
      call mem%dealloc(sigma_aibj, wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o)
!
   end subroutine jacobian_transpose_transformation_mlcc2
!
!
   module subroutine jacobian_transpose_cc2_a1_mlcc2(wf, sigma_ai, c_ai, n_cc2_o, n_cc2_v, &
                                                      first_o, first_v)
!!
!!    Jacobian transpose MLCC2 A1
!!    Written by Sarai D. Folkestad, Apr 2019
!!
!!    Adapted from jacobian_transpose_cc2.F90
!!    written by Sarai D. Folkestad and Alexander Paul, Feb 2019
!!
!!    Calculates the A1 term,
!!
!!       A1: sum_bcjk u^bc_jk (c_bj L_iakc - c_ak g_jbic - c_ci g_jbka)
!!
!!    and adds it to sigma_ai.
!!
!!    Index restrictions:
!!
!!       b, c, j, k : CC2 orbitals
!!
!!       a, i : unrestricted
!!
      implicit none
!
      class(mlcc2) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)     :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)  :: sigma_ai
!
      integer, intent(in) :: n_cc2_o, n_cc2_v, first_o, first_v
!
      real(dp), dimension(:,:,:,:), allocatable :: g_icjb, g_jbka, g_iajb
      real(dp), dimension(:,:,:,:), allocatable :: L_aick
      real(dp), dimension(:,:,:,:), allocatable :: u_cjbk
!
      real(dp), dimension(:,:), allocatable :: X_ck, Y_ik, Y_ca
      real(dp), dimension(:,:), allocatable :: c_bj_active
!
      integer :: b, j, A, I, c, k
!
!     Term 1: sum_bjck c_bj u_bjck L_iakc
!
!     X_ck = sum_bj u_ckbj * c_bj
!
      call mem%alloc(X_ck, n_cc2_v, n_cc2_o)
      call mem%alloc(c_bj_active, n_cc2_v, n_cc2_o)
!
!$omp parallel do private (b, j) collapse(2)
      do j = 1, n_cc2_o
         do b = 1, n_cc2_v
!
            c_bj_active(b, j) = c_ai(b + first_v - 1, j + first_o - 1)
!
         enddo
      enddo
!$omp end parallel do
!
      call dgemm('N', 'N',             &
                  (n_cc2_o)*(n_cc2_v), &
                  1,                   &
                  (n_cc2_o)*(n_cc2_v), &
                  one,                 &
                  wf%u,                & ! u_ck_bj
                  (n_cc2_o)*(n_cc2_v), &
                  c_bj_active,         & ! c_bj
                  (n_cc2_o)*(n_cc2_v), &
                  zero,                &
                  X_ck,                &
                  (n_cc2_o)*(n_cc2_v))
!

      call mem%dealloc(c_bj_active, n_cc2_v, n_cc2_o)
!
      call mem%alloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call wf%get_ovov(g_iajb)
!
!     L_iakc = 2 g_iakc - g_icka
!
      call mem%alloc(L_aick, wf%n_v, wf%n_o, n_cc2_v, n_cc2_o)
!
!$omp parallel do private(k, c, i, a) collapse(2)
       do k = 1, n_cc2_o
         do c = 1, n_cc2_v
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  L_aick(a,i, c, k) = two*g_iajb(i,a, k + first_o - 1, c + first_v - 1) &
                                       - g_iajb(i, c + first_v - 1, k + first_o - 1, a)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
!     sum_kc L_iakc X_ck
!
      call dgemm('N', 'N',             &
                  (wf%n_v)*(wf%n_o),   &
                  1,                   &
                  (n_cc2_v)*(n_cc2_o), &
                  one,                 &
                  L_aick,              & ! L_ai_ck
                  (wf%n_v)*(wf%n_o),   &
                  X_ck,                & ! X_ck
                  (n_cc2_v)*(n_cc2_o), &
                  one,                 &
                  sigma_ai,            &
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(L_aick, wf%n_v, wf%n_o, n_cc2_v, n_cc2_o)
      call mem%dealloc(X_ck, n_cc2_v, n_cc2_o)
!
!     Term 2: - sum_bjck c_ak u_bjck g_jbic
!
      call mem%alloc(g_icjb, wf%n_o, n_cc2_v, n_cc2_o, n_cc2_v)
!
!$omp parallel do private(i,c,j,b) collapse(2)
      do b = 1, n_cc2_v
         do j = 1, n_cc2_o
            do c = 1, n_cc2_v
               do i = 1, wf%n_o              
!
                  g_icjb(i, c, j, b) = g_iajb(i, c + first_v - 1, j + first_o - 1, b + first_v - 1)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%alloc(u_cjbk, n_cc2_v, n_cc2_o, n_cc2_v, n_cc2_o)
      call sort_1234_to_3214(wf%u, u_cjbk, n_cc2_v, n_cc2_o, n_cc2_v, n_cc2_o)
!
      call mem%alloc(Y_ik, wf%n_o, n_cc2_o)
!
!     Y_ik = sum_cjb g_jbic u_bjck
!
      call dgemm('N', 'N',                &
                  wf%n_o,                 &
                  n_cc2_o,                &
                  (n_cc2_o)*(n_cc2_v**2), &
                  one,                    &
                  g_icjb,                 & 
                  wf%n_o,                 &
                  u_cjbk,                 & 
                  (n_cc2_o)*(n_cc2_v**2), &
                  zero,                   &
                  Y_ik,                   &
                  wf%n_o)
!
      call mem%dealloc(g_icjb, wf%n_o, n_cc2_v, n_cc2_o, n_cc2_v)
!
!     sigma_ai -= c_ak Y_ik
!
      call dgemm('N', 'T',          & 
                  wf%n_v,           &
                  wf%n_o,           &
                  n_cc2_o,          &
                  -one,             &
                  c_ai(1, first_o), & ! c_a_k
                  wf%n_v,           &
                  Y_ik,             & ! Y_i_k
                  wf%n_o,           &
                  one,              &
                  sigma_ai,         &
                  wf%n_v)
!
      call mem%dealloc(Y_ik, wf%n_o, n_cc2_o)
!
!     Term 3: - sum_bjck c_ci u^bc_jk g_jbka
!
      call mem%alloc(g_jbka, n_cc2_o, n_cc2_v, n_cc2_o, wf%n_v)
!
!$omp parallel do private(j,b,k,a) collapse(2)
      do a = 1, wf%n_v
         do k = 1, n_cc2_o
            do b = 1, n_cc2_v
               do j = 1, n_cc2_o
!
                 g_jbka(j, b, k, a) = g_iajb(j + first_o - 1, b + first_v - 1, k + first_o - 1, a)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%alloc(Y_ca, n_cc2_v, wf%n_v)
!
!     Y_ca = sum_jbk u_cjbk * g_jbka
!
      call dgemm('N', 'N',                &
                  n_cc2_v,                &
                  wf%n_v,                 &
                  (n_cc2_v)*(n_cc2_o**2), &
                  one,                    &
                  u_cjbk,                 & ! u_c_jbk
                  n_cc2_v,                &
                  g_jbka,                 & ! g_jbk_a
                  (n_cc2_v)*(n_cc2_o**2), &
                  zero,                   &
                  Y_ca,                   & ! Y_ca
                  n_cc2_v)
!
      call mem%dealloc(g_jbka, n_cc2_o, n_cc2_v, n_cc2_o, wf%n_v)
      call mem%dealloc(u_cjbk, n_cc2_v, n_cc2_o, n_cc2_v, n_cc2_o)
!
!     sigma_ai -= c_ci * Y_ca
!
      call dgemm('T', 'N',          & 
                  wf%n_v,           &
                  wf%n_o,           &
                  n_cc2_v,          &
                  -one,             &
                  Y_ca,             & 
                  n_cc2_v,          &
                  c_ai(first_v, 1), & ! c_c_i
                  wf%n_v,           &
                  one,              &
                  sigma_ai,         &
                  wf%n_v)
!
      call mem%dealloc(Y_ca, n_cc2_v, wf%n_v)

   end subroutine jacobian_transpose_cc2_a1_mlcc2
!
!
   module subroutine jacobian_transpose_cc2_b1_mlcc2(wf, sigma_ai, c_aibj, n_cc2_o, n_cc2_v, &
                                                      first_o, first_v, last_o, last_v)
!!
!!    Jacobian transpose MLCC2 B1
!!    Written by Sarai D. Folkestad, Apr 2019
!!
!!    Adapted from jacobian_transpose_cc2.F90
!!    written by Sarai D. Folkestad and Alexander Paul, Feb 2019
!!
!!    Calculates the B1 term,
!!
!!       B1: sum_bjc c_bjci g_bjca - sum_bjk c_akbj g_bjik
!!
!!    and adds it to sigma_ai.
!!
!!    Index restrictions:
!!
!!       Term 1: 
!!
!!          b, j, c, i : CC2 orbitals
!!
!!          a : unretricted
!!
!!       Term 2: 
!!
!!          a, k, b, j : CC2 orbitals
!!
!!          i : unretricted
!!
      implicit none
!
      class(mlcc2) :: wf
!
      integer, intent(in) :: n_cc2_o, n_cc2_v, first_o, first_v, last_o, last_v
!
      real(dp), dimension(n_cc2_v, n_cc2_o, n_cc2_v, n_cc2_o), intent(in) :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)  :: sigma_ai
!
      real(dp), dimension(:,:,:,:), allocatable :: g_bjca, g_ikbj
!
      type(batching_index) :: batch_a
!
      integer :: req0, req1, current_a_batch
!
!     Term 1: sum_bjc c_bjci g_bjca
!
!        b, j, c, i : CC2 orbitals
!
!        a : unretricted
!
      req0 = (n_cc2_v)*(n_cc2_o)*(wf%integrals%n_j)
      req1 = (n_cc2_v)*(wf%integrals%n_j) + (n_cc2_o)*(n_cc2_v**2)
!
      batch_a = batching_index(wf%n_v)
!
      call mem%batch_setup(batch_a, req0, req1)
!
      do current_a_batch = 1, batch_a%num_batches
!
         call batch_a%determine_limits(current_a_batch)
!
         call mem%alloc(g_bjca, n_cc2_v, n_cc2_o, n_cc2_v, batch_a%length)
!
         call wf%get_vovv(g_bjca,             &
                           first_v, last_v,   &
                           first_o, last_o,   &
                           first_v, last_v,   &
                           batch_a%first, batch_a%last)
!
!        sigma_ai =+ sum_bjc g_abjc c_bjci
!
         call dgemm('T', 'N',                                  & 
                     batch_a%length,                           &
                     n_cc2_o,                                  &
                     (n_cc2_o)*(n_cc2_v**2),                   &
                     one,                                      &
                     g_bjca,                                   & ! g_bjc_a
                     (n_cc2_o)*(n_cc2_v**2),                   &
                     c_aibj,                                   & ! c_bjc_i
                     (n_cc2_o)*(n_cc2_v**2),                   &
                     one,                                      &
                     sigma_ai(batch_a%first, first_o),  &
                     wf%n_v)
!
         call mem%dealloc(g_bjca, n_cc2_v, n_cc2_o, n_cc2_v, batch_a%length)
!
      enddo ! batch_a
!
!     Term 2: sum_bjk c_akbj g_bjik
!
!        a, k, b, j : CC2 orbitals
!
!        i : unretricted
!
      call mem%alloc(g_ikbj, wf%n_o, n_cc2_o, n_cc2_v, n_cc2_o)
!
      call wf%get_oovo(g_ikbj,             &
                        1, wf%n_o,         &
                        first_o, last_o,   &
                        first_v, last_v,   &
                        first_o, last_o)
!
!     - sum_bjk c_akbj g_ikbj
!
      call dgemm('N', 'T',               & 
                  n_cc2_v,               &
                  wf%n_o,                &
                  (n_cc2_v)*(n_cc2_o)**2,&
                  -one,                  &
                  c_aibj,                & ! c_a_kbj
                  n_cc2_v,               &
                  g_ikbj,                & ! g_i_kbj
                  wf%n_o,                &
                  one,                   &
                  sigma_ai(first_v, 1),  &
                  wf%n_v)
!
      call mem%dealloc(g_ikbj, wf%n_o, n_cc2_o, n_cc2_v, n_cc2_o)
!
   end subroutine jacobian_transpose_cc2_b1_mlcc2
!
!
   module subroutine jacobian_transpose_cc2_a2_mlcc2(wf, sigma_aibj, c_ai, n_cc2_o, n_cc2_v, &
                                                      first_o, first_v, last_o, last_v)
!!
!!    Jacobian transpose MLCC2 A2
!!    Written by Sarai D. Folkestad, Apr 2019
!!
!!    Adapted from jacobian_transpose_cc2.F90
!!    written by Sarai D. Folkestad and Alexander Paul, Feb 2019
!!
!!    Calculates the A2 term,
!!
!!       A2:  2F_jb c_ai - F_ib c_aj - L_ikjb c_ak + L_cajb c_ci
!!
!!    and adds it to sigma_aibj.
!!
!!    Term 4 is calculated in batches of index c.
!!
!!    Index restrictions:
!!
!!       a, i, b, j : CC2 orbitals
!!
!!       Term 3:
!!
!!          k : unrestricted 
!!
!!       Term 4:
!!
!!          c : unrestricted 
!!
      implicit none
!
      class(mlcc2) :: wf
      integer, intent(in) :: n_cc2_o, n_cc2_v, first_o, first_v, last_o, last_v
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                          :: c_ai
      real(dp), dimension(n_cc2_v, n_cc2_o, n_cc2_v, n_cc2_o), intent(inout)   :: sigma_aibj
!
!     Local variables
!
      real(dp), dimension(:,:,:,:), allocatable :: g_iKjb, g_Cajb
      real(dp), dimension(:,:,:,:), allocatable :: L_Kibj, L_Cajb
      real(dp), dimension(:,:,:,:), allocatable :: sigma_iajb
!
      type(batching_index) :: batch_c
!
      integer :: req0, req1, current_c_batch
!
      integer :: a, i, b, j
!
!     Terms 1 and 2: (2F_jb c_ai - F_ib c_aj)
!
!        a, i, b, j : CC2 orbitals
!
!$omp parallel do private(a, i, b, j)
      do b = 1, n_cc2_v
         do j = 1, n_cc2_o
            do i = 1, n_cc2_o
               do a = 1, n_cc2_v
!
                  sigma_aibj(a, i, b, j) = sigma_aibj(a, i, b, j) &
                                          +  two*wf%fock_ia(j + first_o - 1, b + first_v - 1)&
                                                *c_ai(a + first_v - 1, i + first_o - 1)&
                                          -  wf%fock_ia(i + first_o - 1, b + first_v - 1)&
                                                *c_ai(a + first_v - 1, j + first_o - 1)

!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
!     Term 3: - L_ikjb c_ak
!
!        a, i, b, j : CC2 orbitals
!
!        k : unrestricted
!
!     L_ikjb = 2 g_ikjb - g_jkib (ordered as g_kibj)
!
      call mem%alloc(g_ikjb, n_cc2_o, wf%n_o, n_cc2_o, n_cc2_v)
!
      call wf%get_ooov(g_ikjb,            &
                        first_o, last_o,  &
                        1, wf%n_o,        &
                        first_o, last_o,  &
                        first_v, last_v)
!
      call mem%alloc(L_kibj, wf%n_o, n_cc2_o, n_cc2_v, n_cc2_o)
      call zero_array(L_kibj, wf%n_o*(n_cc2_o**2)*n_cc2_v)
!
      call add_2143_to_1234(two, g_ikjb, L_kibj, wf%n_o, n_cc2_o, n_cc2_v, n_cc2_o)
      call add_4123_to_1234(-one, g_ikjb, L_kibj, wf%n_o, n_cc2_o, n_cc2_v, n_cc2_o)
!
      call mem%dealloc(g_ikjb, n_cc2_o, wf%n_o, n_cc2_o, n_cc2_v)
!
      call dgemm('N', 'N',                &
                  n_cc2_v,                &
                  (n_cc2_o**2)*n_cc2_v,   &
                  wf%n_o,                 &
                  -one,                   &
                  c_ai(first_v, 1),       & ! c_a_k
                  wf%n_v,                 &
                  L_kibj,                 & ! L_k_ibj
                  wf%n_o,                 &
                  one,                    &
                  sigma_aibj,             & ! sigma_a_ibj
                  n_cc2_v)
!
      call mem%dealloc(L_kibj, wf%n_o, n_cc2_o, n_cc2_v, n_cc2_o)
!
!     Term 4: L_cajb c_ci
!
!        a, i, b, j : CC2 orbitals
!
!        c : unrestricted
!
      call mem%alloc(sigma_iajb, n_cc2_o, n_cc2_v, n_cc2_o, n_cc2_v)
      call zero_array(sigma_iajb, (n_cc2_o**2)*(n_cc2_v**2))
!
      req0 = (n_cc2_v)*(n_cc2_o)*(wf%integrals%n_j)
      req1 = max((n_cc2_v)*(wf%integrals%n_j) + (n_cc2_o)*(n_cc2_v)**2, 2*(n_cc2_o)*(n_cc2_v)**2)
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
         call mem%alloc(g_cajb, batch_c%length, n_cc2_v, n_cc2_o, n_cc2_v)
!
         call wf%get_vvov(g_cajb, &
                           batch_c%first, batch_c%last,   &
                           first_v, last_v, &
                           first_o, last_o, &
                           first_v, last_v)
!
         call mem%alloc(L_cajb, batch_c%length, n_cc2_v, n_cc2_o, n_cc2_v)
         call dcopy((batch_c%length)*(n_cc2_v**2)*(n_cc2_o), g_cajb, 1, L_cajb, 1)
         call dscal((batch_c%length)*(n_cc2_v**2)*(n_cc2_o), two, L_cajb, 1)
         call add_1432_to_1234(-one, g_cajb, L_cajb, batch_c%length, n_cc2_v, n_cc2_o, n_cc2_v)
!
         call mem%dealloc(g_cajb, batch_c%length, n_cc2_v, n_cc2_o, n_cc2_v)
!
         call dgemm('T', 'N',                         &
                     n_cc2_o,                         &
                     (n_cc2_v**2)*(n_cc2_o),          &
                     batch_c%length,                  &
                     one,                             &
                     c_ai(batch_c%first, first_o),    & ! c_c_i
                     wf%n_v,                          &
                     L_cajb,                          & ! L_c_ajb
                     batch_c%length,                  &
                     one,                             &
                     sigma_iajb,                      &
                     n_cc2_o)
!
         call mem%dealloc(L_cajb, batch_c%length, n_cc2_v, n_cc2_o, n_cc2_v)
!
      enddo ! batch_c
!
      call add_2143_to_1234(one, sigma_iajb, sigma_aibj, n_cc2_v, n_cc2_o, n_cc2_v, n_cc2_o)
!
      call mem%dealloc(sigma_iajb, n_cc2_o, n_cc2_v, n_cc2_o, n_cc2_v)
!
   end subroutine jacobian_transpose_cc2_a2_mlcc2
!
!
   module subroutine jacobian_transpose_cc2_b2_mlcc2(wf, sigma_aibj, c_aibj)
!!
!!    Jacobian transpose CC2 B2
!!    Written by Sarai D. Folkestad, Apr 2019
!!
!!    Adapted from jacobian_transpose_cc2.F90
!!    written by Sarai D. Folkestad and Alexander Paul, Feb 2019
!!
!!    Calculates the B2 term,
!!
!!       B2: ε_aibj c_aibj
!!
!!    and adds it to sigma_aibj.
!!
!!    Index restrictions:
!!
!!       a, i, b, j : CC2 orbitals 
!!
      implicit none
!
      class(mlcc2) :: wf
!
      real(dp), dimension(wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o), intent(in)      :: c_aibj
      real(dp), dimension(wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o), intent(inout)   :: sigma_aibj
!
      integer :: a, i, b, j
!
!$omp parallel do private(a, i, b, j) collapse(2)
      do j = 1, wf%n_cc2_o
         do b = 1, wf%n_cc2_v
            do i = 1, wf%n_cc2_o
               do a = 1, wf%n_cc2_v
!
                  sigma_aibj(a, i, b, j) = sigma_aibj(a, i, b, j) + c_aibj(a,i,b,j) &
                                          *(wf%orbital_energies(a + wf%first_cc2_v - 1 + wf%n_o) &
                                          + wf%orbital_energies(b + wf%first_cc2_v - 1 + wf%n_o) &
                                          - wf%orbital_energies(i + wf%first_cc2_o - 1) &
                                          - wf%orbital_energies(j + wf%first_cc2_o - 1))
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine jacobian_transpose_cc2_b2_mlcc2
!
!
end submodule jacobian_transpose_mlcc2
