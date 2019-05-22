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
submodule (cc2_class) jacobian_cc2
!
!!
!!    Jacobian submodule (CC2)
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 2018
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
! 
   implicit none
!
!
contains
!
!
   module subroutine prepare_for_jacobian_cc2(wf)
!!
!!    Prepare for Jacobian
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2019
!!
      implicit none 
!
      class(cc2), intent(inout) :: wf 
!
      call wf%initialize_u()
      call wf%construct_u()
!
   end subroutine prepare_for_jacobian_cc2
!
!
   module subroutine jacobian_transform_trial_vector_cc2(wf, c_i)
!!
!!    Jacobian transform trial vector 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep 2018
!!
      class(cc2), intent(in) :: wf 
!
      real(dp), dimension(wf%n_es_amplitudes) :: c_i
!
      call wf%jacobian_cc2_transformation(c_i)
!
   end subroutine jacobian_transform_trial_vector_cc2
!
!
   module subroutine jacobian_cc2_transformation_cc2(wf, c)
!!
!!    Jacobian transformation (cc2)
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, Jan 2019
!!
!!    Directs the transformation by the cc2 Jacobi matrix,
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
      class(cc2) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes)   :: c
!
      real(dp), dimension(:,:), allocatable     :: c_ai
      real(dp), dimension(:,:,:,:), allocatable :: c_aibj
!
      real(dp), dimension(:,:), allocatable     :: rho_ai   
      real(dp), dimension(:,:,:,:), allocatable :: rho_aibj
!
!
      integer :: i, j, a, b, ai, bj, aibj ! Index
!
      type(timings) :: timer
!
      timer = new_timer('Jacobian transformation CC2')
      call timer%turn_on()
!
!     Allocate and zero the transformed vector (singles part)
!
      call mem%alloc(rho_ai, wf%n_v, wf%n_o)
      call zero_array(rho_ai, wf%n_t1)
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
      call wf%jacobian_ccs_a1(rho_ai, c_ai)
      call wf%jacobian_ccs_b1(rho_ai, c_ai)
!
!     :: CC2 contributions to the transformed singles vector ::
!
      call wf%jacobian_cc2_a1(rho_ai, c_ai)
!
!     Allocate the incoming unpacked doubles vector
!
      call mem%alloc(c_aibj, (wf%n_v), (wf%n_o), (wf%n_v), (wf%n_o))
!
!$omp parallel do schedule(static) private(a, i, b, j, ai, bj, aibj) 
      do j = 1, wf%n_o
         do b = 1, wf%n_v
!
            bj = wf%n_v*(j - 1) + b
!  
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  ai = wf%n_v*(i - 1) + a
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
!     Scale the doubles vector by 1 + delta_ai,bj, i.e.
!     redefine to c_ckdl = c_ckdl (1 + delta_ck,dl)
!
!$omp parallel do schedule(static) private(a, i) 
      do i =1, wf%n_o
         do a = 1, wf%n_v
!
            c_aibj(a, i, a, i) = two*c_aibj(a, i, a, i)
!
         enddo
      enddo
!$omp end parallel do
!
      call wf%jacobian_cc2_b1(rho_ai, c_aibj)
!
!     Done with singles vector c; overwrite it with
!     transformed vector for exit
!
!$omp parallel do schedule(static) private(a, i, ai) 
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            ai = wf%n_v*(i - 1) + a
!
            c(ai) = rho_ai(a, i)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(rho_ai, wf%n_v, wf%n_o)
!
!     :: CC2 contributions to the transformed doubles vector ::
!
!     Allocate unpacked transformed vector
!
      call mem%alloc(rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call zero_array(rho_aibj, wf%n_t1**2)
!
!     Contributions from singles vector c
!
      call wf%jacobian_cc2_a2(rho_aibj, c_ai)
!
      call mem%dealloc(c_ai, wf%n_v, wf%n_o)
!
!     Contributions from doubles vector c      
!
      call wf%jacobian_cc2_b2(rho_aibj, c_aibj)
!
      call mem%dealloc(c_aibj, (wf%n_v), (wf%n_o), (wf%n_v), (wf%n_o))
!
!     Overwrite the incoming doubles c vector & pack in
!
!$omp parallel do schedule(static) private(a, i, b, j, ai, bj, aibj) 
      do j = 1, wf%n_o
         do b = 1, wf%n_v
!
            bj = wf%n_v*(j - 1) + b
!  
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  ai = wf%n_v*(i - 1) + a
!
                  if (ai .ge. bj) then
!
                     aibj = ai*(ai-3)/2 + ai + bj
!
                     c((wf%n_o)*(wf%n_v) + aibj) = rho_aibj(a, i, b, j)
!
                  endif
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_cc2_transformation_cc2
!
!
   module subroutine jacobian_cc2_a1_cc2(wf, rho_ai, c_ai)
!!
!!    Jacobian CC2 A1
!!    Written by Sarai D. Folkestad Eirik F. Kjønstad Jan 2019
!!
!!    rho_ai^A1 = sum_bjck (L_kcjb u_aick c_bj - g_kcjb (u_ckbi c_aj + u_ckaj c_bi))
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(out)   :: rho_ai
!
      real(dp), dimension(:,:,:,:), allocatable :: g_jbkc, u_bkci, X_kcji, L_ckbj
      real(dp), dimension(:,:), allocatable     :: X_ji, X_ck
!
      type(timings) :: timer
!
      timer = new_timer('jacobian cc2 a1')
      call timer%turn_on()
!
      call mem%alloc(g_jbkc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_ovov(g_jbkc)
!
!     Reorder u_ckbi as u_bkci
!
      call mem%alloc(u_bkci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_3214(wf%u, u_bkci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(X_ji, wf%n_o, wf%n_o)
!
!     X_ji = sum_kcb g_kcjb u_ckbi 
!
      call dgemm('N', 'N',                &
                  wf%n_o,                 &
                  wf%n_o,                 &
                  (wf%n_v**2)*(wf%n_o),   &
                  one,                    &
                  g_jbkc,                 & ! g_j_bkc
                  wf%n_o,                 &
                  u_bkci,                 & ! u_bkc_i 
                  (wf%n_v**2)*(wf%n_o),   &
                  zero,                   &
                  X_ji,                   &
                  wf%n_o)
!
      call dgemm('N', 'N', &
                  wf%n_v,  &
                  wf%n_o,  &
                  wf%n_o,  &
                  -one,    &
                  c_ai,    & ! c_a_j
                  wf%n_v,  &
                  X_ji,    & ! X_j_i
                  wf%n_o,  &
                  one,     &
                  rho_ai,  &
                  wf%n_v)
!
      call mem%dealloc(X_ji, wf%n_o, wf%n_o)
!
!     NOTE! We will now pretend that u_bkci = u(c, k, b, i) is  u_akcj = u(c, k, a, j)
!     NOTE! We will now pretend that g_jbkc = g_kcjb
!
      call mem%alloc(X_kcji, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N',                &
                  (wf%n_o**2)*(wf%n_v),   &
                  wf%n_o,                 &
                  wf%n_v,                 &
                  one,                    &
                  g_jbkc,                 & ! g_kcj_b
                  (wf%n_o**2)*(wf%n_v),   &
                  c_ai,                   & ! c_b_i
                  wf%n_v,                 &
                  zero,                   &
                  X_kcji,                 & ! X_kcj_i
                  (wf%n_o**2)*(wf%n_v))       
!
      call dgemm('N', 'N',                &
                  wf%n_v,                 &
                  wf%n_o,                 &
                  (wf%n_o**2)*(wf%n_v),   &
                  -one,                   &
                  u_bkci,                 & ! u_a_kcj
                  wf%n_v,                 &
                  X_kcji,                 & ! X_kcj_i
                  (wf%n_o**2)*(wf%n_v),   &
                  one,                    &
                  rho_ai,                 &
                  wf%n_v)
!
      call mem%dealloc(u_bkci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X_kcji, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
!     L_kcjb ordered as L_ckbj
!
      call mem%alloc(L_ckbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call zero_array(L_ckbj, wf%n_t1**2)
!
      call add_2143_to_1234(two, g_jbkc, L_ckbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o) ! Still pretending that g_jbkc is g_kcjb although this is not necessary
      call add_4123_to_1234(-one, g_jbkc, L_ckbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_jbkc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%alloc(X_ck, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',       &
                  wf%n_o*wf%n_v, &
                  1,             &
                  wf%n_o*wf%n_v, &
                  one,           &
                  L_ckbj,        & ! L_ck_bj
                  wf%n_o*wf%n_v, &
                  c_ai,          & ! c_bj
                  wf%n_o*wf%n_v, &
                  zero,          &
                  X_ck,          &
                  wf%n_o*wf%n_v)
!
      call mem%dealloc(L_ckbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',       &
                  wf%n_o*wf%n_v, &
                  1,             &
                  wf%n_o*wf%n_v, &
                  one,           &
                  wf%u,          & ! u_ai_ck
                  wf%n_o*wf%n_v, &
                  X_ck,          &
                  wf%n_o*wf%n_v, &
                  one,           &
                  rho_ai,        &
                  wf%n_o*wf%n_v)
!
      call mem%dealloc(X_ck, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_cc2_a1_cc2
!
!
   module subroutine jacobian_cc2_b1_cc2(wf, rho_ai, c_aibj)
!!
!!    Jacobian CC2 B1
!!    Written by Sarai D. Folkestad Eirik F. Kjønstad Jan 2019
!!
!!    rho_ai^B1 = 2 sum_bj F_jb c_aibj - F_jb c_ajbi  - sum_kjb L_kijb c_akbj + sum_bkc L_abkc c_bick
!!
      implicit none
!
      class(cc2) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)   :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(out)                  :: rho_ai   
!
      real(dp), dimension(:,:,:,:), allocatable :: X_aijb 
      real(dp), dimension(:,:,:,:), allocatable :: g_jbki 
      real(dp), dimension(:,:,:,:), allocatable :: g_abkc
      real(dp), dimension(:,:,:,:), allocatable :: L_abkc
      real(dp), dimension(:,:,:,:), allocatable :: L_kbji
      real(dp), dimension(:,:,:,:), allocatable :: c_bkci 
!
      type(batching_index) :: batch_a 
      integer         :: req0, req1, current_a_batch
!
      type(timings) :: timer
!
      timer = new_timer('jacobian cc2 b1')
      call timer%turn_on()
!
!     Make X_aijb = 2 c_aibj - c_ajbi 
!
      call mem%alloc(X_aijb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call zero_array(X_aijb, wf%n_t1**2)
!
      call add_1243_to_1234(two, c_aibj, X_aijb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call add_1342_to_1234(-one, c_aibj, X_aijb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
!     2 sum_bj F_jb c_aibj - F_jb c_ajbi = sum_bj X_aijb F_jb 
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  1,                 &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  X_aijb,            & ! X_ai,jb 
                  (wf%n_o)*(wf%n_v), &
                  wf%fock_ia,        & ! F_jb 
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  rho_ai,            &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(X_aijb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
!     - sum_kjb c_akbj L_jbki
!
!     Make L_jbki ordered as L_kbji (= 2 g_jbki - g_kbji)
!
      call mem%alloc(g_jbki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call wf%get_ovoo(g_jbki)
!
      call mem%alloc(L_kbji, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call copy_and_scale(-one, g_jbki, L_kbji, (wf%n_o**3)*wf%n_v)
!      
      call add_3214_to_1234(two, g_jbki, L_kbji, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call mem%dealloc(g_jbki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  wf%n_o,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  -one,                 &
                  c_aibj,               & ! c_a,kbj
                  wf%n_v,               &
                  L_kbji,               & ! L_kbj,i 
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  rho_ai,               & ! rho_a,i
                  wf%n_v)
!
      call mem%dealloc(L_kbji, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
!     sum_bkc L_abkc c_ckbi, batch over a 
!
!     Order c_ckbi as c_bkci
!
      call mem%alloc(c_bkci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_3214(c_aibj, c_bkci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      req0 = (wf%n_o)*(wf%n_v)*(wf%integrals%n_J)
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
         call mem%alloc(g_abkc, batch_a%length, wf%n_v, wf%n_o, wf%n_v)
!
         call wf%get_vvov(g_abkc,                       &
                           batch_a%first, batch_a%last, &
                           1, wf%n_v,                   &   
                           1, wf%n_o,                   &
                           1, wf%n_v)
!
         call mem%alloc(L_abkc, batch_a%length, wf%n_v, wf%n_o, wf%n_v)
!
         call copy_and_scale(two, g_abkc, L_abkc, (wf%n_v**2)*(batch_a%length)*wf%n_o)
!
         call add_1432_to_1234(-one, g_abkc, L_abkc, batch_a%length, wf%n_v, wf%n_o, wf%n_v)
!
         call mem%dealloc(g_abkc, batch_a%length, wf%n_v, wf%n_o, wf%n_v)
!
         call dgemm('N','N',                   &
                     batch_a%length,           &
                     wf%n_o,                   &
                     (wf%n_o)*(wf%n_v)**2,     &
                     one,                      &
                     L_abkc,                   & ! L_a,bkc 
                     batch_a%length,           &
                     c_bkci,                   & ! c_bkc,i
                     (wf%n_o)*(wf%n_v)**2,     &
                     one,                      &
                     rho_ai(batch_a%first, 1), &
                     wf%n_v)
!
         call mem%dealloc(L_abkc, batch_a%length, wf%n_v, wf%n_o, wf%n_v)
!
      enddo ! batch_a
!
      call mem%dealloc(c_bkci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_cc2_b1_cc2
!
!
   module subroutine jacobian_cc2_a2_cc2(wf, rho_aibj, c_ai)
!!
!!    Jacobian CC2 A2
!!    Written by Sarai D. Folkestad Eirik F. Kjønstad Jan 2019
!!
!!    rho_aibj^A2 = (1/Δ_aibj)P_aibj (sum_c g_aibc c_cj - sum_k g_aikj c_bk), 
!!
      implicit none
!
      class(cc2) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                   :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(out)  :: rho_aibj   
!
      real(dp), dimension(:,:,:,:), allocatable :: rho_bjai
      real(dp), dimension(:,:,:,:), allocatable :: g_kjai 
      real(dp), dimension(:,:,:,:), allocatable :: g_aibc
!
      integer :: a, i
!
      type(batching_index) :: batch_c 
      integer         :: req0, req1, current_c_batch 
!
      type(timings) :: timer
!
      timer = new_timer('jacobian cc2 a2')
      call timer%turn_on()
!
!     (1/Δ_aibj)P_aibj (- sum_k c_bk g_kjai)
!
      call mem%alloc(g_kjai, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%get_oovo(g_kjai)
!
      call mem%alloc(rho_bjai, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_o,               &
                  -one,                 &
                  c_ai,                 & ! c_b,k
                  wf%n_v,               &
                  g_kjai,               & ! g_k,jai
                  wf%n_o,               &
                  zero,                 &
                  rho_bjai,             & ! rho_b,jai 
                  wf%n_v)
!
      call mem%dealloc(g_kjai, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     (1/Δ_aibj)P_aibj (sum_c g_aibc c_cj)
!
      req0 = (wf%integrals%n_J)*(wf%n_o)*(wf%n_v)
      req1 = (wf%n_o)*(wf%n_v)**2 + (wf%integrals%n_J)*(wf%n_v)
!
      call batch_c%init(wf%n_v)
!
      call mem%batch_setup(batch_c, req0, req1)
!
      do current_c_batch = 1, batch_c%num_batches
!
         call batch_c%determine_limits(current_c_batch)
!
         call mem%alloc(g_aibc, wf%n_v, wf%n_o, wf%n_v, batch_c%length)
!
         call wf%get_vovv(g_aibc,                        &
                           1, wf%n_v,                    &
                           1, wf%n_o,                    &
                           1, wf%n_v,                    &
                           batch_c%first, batch_c%last)
!
         call dgemm('N','N',                 &
                     (wf%n_o)*(wf%n_v)**2,   &
                     wf%n_o,                 &
                     batch_c%length,         &
                     one,                    &
                     g_aibc,                 & ! g_aib,c 
                     (wf%n_o)*(wf%n_v)**2,   &
                     c_ai(batch_c%first, 1), & ! c_c,j
                     wf%n_v,                 &
                     one,                    &
                     rho_bjai,               & ! rho_aib,j But we will symmetrize and may place in rho_bjai
                     (wf%n_o)*(wf%n_v)**2)
!
         call mem%dealloc(g_aibc, wf%n_v, wf%n_o, wf%n_v, batch_c%length)
!
      enddo  ! batch over c 
!
!$omp parallel do private(a, i)
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            rho_bjai(a, i, a, i) = half*rho_bjai(a, i, a, i)
!
         enddo
      enddo
!$omp end parallel do 
!
      call daxpy((wf%n_o)**2*(wf%n_v)**2, one, rho_bjai, 1, rho_aibj, 1)
      call add_21_to_12(one, rho_bjai, rho_aibj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))     
!
      call mem%dealloc(rho_bjai, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_cc2_a2_cc2
!
!
   module subroutine jacobian_cc2_b2_cc2(wf, rho_aibj, c_aibj)
!!
!!    Jacobian CC2 B2
!!    Written by Sarai D. Folkestad Eirik F. Kjønstad Jan 2019
!!
!!    rho_aibj^B2 = ε_aibj c_aibj/(1/Δ_aibj) 
!!
      implicit none
!
      class(cc2) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout)   :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(out)     :: rho_aibj   
!
      integer :: i, j, a, b
!
      type(timings) :: timer
!
      timer = new_timer('jacobian cc2 b2')
      call timer%turn_on()
!
!     c_aibj/(1/Δ_aibj) 
!
!$omp parallel do private(a, i)
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            c_aibj(a, i, a, i) = half*c_aibj(a, i, a, i) 
!
         enddo
      enddo
!$omp end parallel do
!
!
!$omp parallel do private(a, i, b, j)
      do j = 1, wf%n_o
         do b = 1, wf%n_v
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  rho_aibj(a,i,b,j) = rho_aibj(a,i,b,j) + c_aibj(a,i,b,j)*&
                                          (- wf%orbital_energies(i) &
                                           - wf%orbital_energies(j) &
                                           + wf%orbital_energies(wf%n_o + a) &
                                           + wf%orbital_energies(wf%n_o + b) )
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call timer%turn_off()
!
   end subroutine jacobian_cc2_b2_cc2
!
!
end submodule jacobian_cc2
