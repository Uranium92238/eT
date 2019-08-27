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
submodule (lowmem_cc2_class) jacobian
!
!!
!!    Jacobian submodule (CC2)
!!    Written by Eirik F. Kjønstad and Sarai Dery Folkestad
!!    Linda Goletto, and Alexander Paul, Dec 2018
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
   module subroutine effective_jacobian_transformation_lowmem_cc2(wf, omega, c)
!!
!!    Effective jacobian transformation
!!    Written by Eirik F. Kjønstad and Sarai Dery Folkestad
!!    Linda Goletto, and Alexander Paul, Dec 2018
!!
      implicit none
!
      class(lowmem_cc2) :: wf
!
      real(dp), intent(in) :: omega
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: c
!
      real(dp), dimension(:,:), allocatable :: c_a_i
!
      real(dp), dimension(:,:), allocatable :: rho_a_i
!
      real(dp), dimension(:), allocatable :: eps_o
      real(dp), dimension(:), allocatable :: eps_v
!
!     Allocate and zero the transformed vector (singles part)
!
      call mem%alloc(rho_a_i, wf%n_v, wf%n_o)
      call zero_array(rho_a_i, wf%n_o*wf%n_v)
!
      call mem%alloc(c_a_i, wf%n_v, wf%n_o)
      call dcopy(wf%n_t1, c, 1, c_a_i, 1)
!
!     :: CCS contributions to the singles c vector ::
!
      call wf%jacobian_ccs_a1(rho_a_i, c_a_i)
!
!     :: CC2 contributions to the transformed singles vector ::
!
      call mem%alloc(eps_o, wf%n_o)
      call mem%alloc(eps_v, wf%n_v)
!
      eps_o = wf%orbital_energies(1:wf%n_o)
      eps_v = wf%orbital_energies(wf%n_o + 1 : wf%n_mo)
!
      call wf%jacobian_cc2_a1(rho_a_i, c_a_i)
      call wf%jacobian_cc2_b1(rho_a_i, c_a_i, eps_o, eps_v)
!
      call wf%effective_jacobian_cc2_a1(omega, rho_a_i, c_a_i, eps_o, eps_v)
      call wf%effective_jacobian_cc2_b1(omega, rho_a_i, c_a_i, eps_o, eps_v)
      call wf%effective_jacobian_cc2_c1(omega, rho_a_i, c_a_i, eps_o, eps_v)
      call wf%effective_jacobian_cc2_d1(omega, rho_a_i, c_a_i, eps_o, eps_v)
      call wf%effective_jacobian_cc2_e1(omega, rho_a_i, c_a_i, eps_o, eps_v)
      call wf%effective_jacobian_cc2_f1(omega, rho_a_i, c_a_i, eps_o, eps_v)
!
      call mem%dealloc(eps_o, wf%n_o)
      call mem%dealloc(eps_v, wf%n_v)
!
      call dcopy(wf%n_es_amplitudes, rho_a_i, 1, c, 1)
!
      call mem%dealloc(c_a_i, wf%n_v, wf%n_o)
      call mem%dealloc(rho_a_i, wf%n_v, wf%n_o)
!
   end subroutine effective_jacobian_transformation_lowmem_cc2
!
!
   module subroutine jacobian_cc2_a1_lowmem_cc2(wf, rho_ai, c_bj)
!!
!!    Jacobian CC2 A1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,
!!    Linda Goletto, and Alexander Paul, Dec 2018
!!
!!    rho_ai =+ sum_bj (2 g_aijb - g_abji) * c_bj
!!
!!    Separate calculation of both terms due to batching
!!
      implicit none
!
      class(lowmem_cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_bj
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
!
      real(dp), dimension(:,:), allocatable :: c_jb
!
      real(dp), dimension(:,:,:,:), allocatable :: g_aijb, g_abji
!
      integer :: req0, req1_i, req1_j, req1_b, req2
!
      integer :: current_i_batch, current_j_batch, current_b_batch
!
      integer :: rho_offset, j, b
!
      type(batching_index) :: batch_i, batch_j, batch_b
!
!     :: Term 1: rho_ai = sum_bj 2 g_aijb * c_bj ::
!
      req0 = 0
!
      req1_i = (wf%integrals%n_J)*(wf%n_v)
      req1_j = (wf%integrals%n_J)*(wf%n_v)
!
      req2 = (wf%n_v)**2
!
      call batch_i%init(wf%n_o)
      call batch_j%init(wf%n_o)
!
      call mem%batch_setup(batch_i, batch_j, req0, req1_i, req1_j, req2)
!
      do current_i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(current_i_batch)
!
         do current_j_batch = 1, batch_j%num_batches
!
            call batch_j%determine_limits(current_j_batch)
!
            call mem%alloc(g_aijb, wf%n_v, (batch_i%length), (batch_j%length), wf%n_v)
!
            call wf%get_voov(g_aijb,                        &
                              1, wf%n_v,                    &
                              batch_i%first, batch_i%last,  &
                              batch_j%first, batch_j%last,  &
                              1, wf%n_v)
!
            call mem%alloc(c_jb, (batch_j%length), wf%n_v)
!
!$omp parallel do private(j, b)
            do b = 1, wf%n_v
               do j = 1, batch_j%length
!
                  c_jb(j, b) = c_bj(b, j + batch_j%first - 1)
!
               enddo
            enddo
!$omp end parallel do
!
!           rho_a_i = rho_a_i + sum_bj 2 g_aijb * c_bj
!
            rho_offset = wf%n_v*(batch_i%first - 1) + 1
!
            call dgemm('N', 'N',                   &
                        (wf%n_v)*(batch_i%length), &
                        1,                         &
                        (batch_j%length)*(wf%n_v), &
                        two,                       &
                        g_aijb,                    & ! g_ai_jb
                        (wf%n_v)*(batch_i%length), &
                        c_jb,                      & ! c_jb
                        (batch_j%length)*(wf%n_v), &
                        one,                       &
                        rho_ai(rho_offset, 1),     &
                        (wf%n_v)*(wf%n_o))
!
            call mem%dealloc(c_jb, (batch_j%length), wf%n_v)
            call mem%dealloc(g_aijb, wf%n_v,(batch_i%length), (batch_j%length), wf%n_v)
!
         enddo ! batch_j
      enddo !batch_i
!
!     :: Term 2 rho_ai = - g_abji * c_bj::
!
      req1_i = max((wf%integrals%n_J)*(wf%n_v), (wf%integrals%n_J)*(wf%n_o))
      req1_b = max((wf%integrals%n_J)*(wf%n_v), (wf%integrals%n_J)*(wf%n_o))
!
      req2 = 2*(wf%n_o)*(wf%n_v)
!
      call batch_i%init(wf%n_o)
      call batch_b%init(wf%n_v)
!
      call mem%batch_setup(batch_i, batch_b, req0, req1_i, req1_b, req2)
!
      do current_i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(current_i_batch)
!
         do current_b_batch = 1, batch_b%num_batches
!
            call batch_b%determine_limits(current_b_batch)
!
            call mem%alloc(g_abji, wf%n_v, (batch_b%length), wf%n_o, (batch_i%length))
!
            call wf%get_vvoo(g_abji,                        &
                              1, wf%n_v,                    &
                              batch_b%first, batch_b%last,  &
                              1, wf%n_o,                    &
                              batch_i%first, batch_i%last)
!
!           Sort g_abji(a,b,j,i) as g_abji(a,i,j,b)
!
            call mem%alloc(g_aijb, wf%n_v, (batch_i%length), wf%n_o, (batch_b%length))
            call sort_1234_to_1432(g_abji, g_aijb, wf%n_v, (batch_b%length), wf%n_o, (batch_i%length))
!
            call mem%dealloc(g_abji, wf%n_v, (batch_b%length), wf%n_o, (batch_i%length))
!
            call mem%alloc(c_jb, wf%n_o, (batch_b%length))
!
!$omp parallel do private(j, b)
            do j = 1, wf%n_o
               do b = 1, batch_b%length
!
                  c_jb(j, b) = c_bj(b + batch_b%first - 1, j)
!
               enddo
            enddo
!$omp end parallel do
!
!           rho_a_i = rho_a_i - sum_bj g_aijb * c_jb
!
            rho_offset = wf%n_v*(batch_i%first - 1) + 1
!
            call dgemm('N', 'N',                   &
                        (wf%n_v)*(batch_i%length), &
                        1,                         &
                        (wf%n_o)*(batch_b%length), &
                        -one,                      &
                        g_aijb,                    & ! g_ai_jb
                        (wf%n_v)*(batch_i%length), &
                        c_jb,                      & ! c_jb
                        (wf%n_o)*(batch_b%length), &
                        one,                       & ! rho_ai
                        rho_ai(rho_offset, 1),     &
                        (wf%n_v)*(wf%n_o))
!
            call mem%dealloc(g_aijb, wf%n_v, (batch_i%length), wf%n_o,(batch_b%length))
            call mem%dealloc(c_jb, wf%n_o, (batch_b%length))
!
         enddo ! batch_b
      enddo ! batch_i
!
end subroutine jacobian_cc2_a1_lowmem_cc2
!
!
   module subroutine jacobian_cc2_b1_lowmem_cc2(wf, rho_ai, c_bj, eps_o, eps_v)
!!
!!    Jacobian CC2 B1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,
!!    Linda Goletto, and Alexander Paul, Dec 2018
!!
!!    rho_ai^B1 = L_kcjb c_bj (2 t^ac_ik - t^ac_ki)
!!                - L_kcjb t^cb_ki c_aj - L_kcjb t^ca_kj c_bi
!!
      implicit none
!
      class(lowmem_cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_bj
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
      real(dp), dimension(:,:), allocatable :: X_kc, X_ji, X_ab, X_kc_batch
      real(dp), dimension(:,:), allocatable :: rho_ai_batch
!
      real(dp), dimension(:,:,:,:), allocatable :: g_kcjb, g_ckbi, g_ckaj, g_aick
      real(dp), dimension(:,:,:,:), allocatable :: L_kcbj, L_kcjb, L_jckb
      real(dp), dimension(:,:,:,:), allocatable :: u_aikc, t_akcj, g_jckb
!
      integer :: i, j, k, a, b, c, bj_offset, kc_offset
!
      integer :: req0, req1_j, req1_k, req2, req1_a, req1_c
      integer :: req1_i, req2_ji, req2_ki, req2_kj, req3
!
      integer :: current_j_batch, current_k_batch, current_i_batch
      integer :: current_a_batch, current_c_batch
!
      type(batching_index) :: batch_j, batch_k, batch_a, batch_c, batch_i
!
!     :: Term 1: L_kcjb * c_bj * (2 t^ac_ik - t^ac_ki)  ::
!                L_kcjb * c_bj * u_aick
!
!     X_kc = sum_jb L_kcbj * c_bj
!     In batches of k and j
!
      call mem%alloc(X_kc, wf%n_o, wf%n_v)
      X_kc=zero
!
      req0 = 0
!
      req1_j = (wf%integrals%n_J)*(wf%n_v)
      req1_k = (wf%integrals%n_J)*(wf%n_v)
!
      req2 = 2*(wf%n_v)**2
!
      call batch_j%init(wf%n_o)
      call batch_k%init(wf%n_o)
!
      call mem%batch_setup(batch_j, batch_k, req0, req1_j, req1_k, req2)
!
      do current_j_batch = 1, batch_j%num_batches
!
         call batch_j%determine_limits(current_j_batch)
!
         do current_k_batch = 1, batch_k%num_batches
!
            call batch_k%determine_limits(current_k_batch)
!
!           L_kcjb = 2 g_kcjb - g_kbjc, ordered as L_kcbj
!
            call mem%alloc(g_kcjb, batch_k%length, wf%n_v, batch_j%length, wf%n_v)
!
            call wf%get_ovov(g_kcjb,                        &
                              batch_k%first, batch_k%last,  &
                              1, wf%n_v,                    &
                              batch_j%first, batch_j%last,  &
                              1, wf%n_v)
!
            call mem%alloc(L_kcbj, batch_k%length, wf%n_v, wf%n_v, batch_j%length)
!
            call zero_array(L_kcbj, batch_k%length*batch_j%length*wf%n_v**2)
!
            call add_1243_to_1234(two, g_kcjb, L_kcbj, batch_k%length, wf%n_v, wf%n_v, batch_j%length)
            call add_1342_to_1234(-one, g_kcjb, L_kcbj, batch_k%length, wf%n_v, wf%n_v, batch_j%length)
!
            call mem%dealloc(g_kcjb, batch_k%length, wf%n_v, batch_j%length, wf%n_v)
!
            call mem%alloc(X_kc_batch, batch_k%length, wf%n_v)
!
            bj_offset = wf%n_v*(batch_j%first - 1) + 1
!
!           X_kc = L_kcjb * c_bj
!
            call dgemm('N', 'N',                   &
                        (batch_k%length)*(wf%n_v), &
                        1,                         &
                        (wf%n_v)*(batch_j%length), &
                        one,                       &
                        L_kcbj,                    & ! L_kc_bj
                        (batch_k%length)*(wf%n_v), &
                        c_bj(bj_offset, 1),        & ! c_bj
                        (wf%n_v)*(wf%n_o),         &
                        zero,                      &
                        X_kc_batch,                & ! X_kc
                        batch_k%length*(wf%n_v))
!
!$omp parallel do private(k, c)
            do c = 1, wf%n_v
               do k = 1, batch_k%length
!
                  X_kc(k + batch_k%first - 1, c) = X_kc(k + batch_k%first - 1, c) &
                                                 + X_kc_batch(k, c)
!
               enddo
            enddo
!$omp end parallel do
!
            call mem%dealloc(X_kc_batch, batch_k%length, wf%n_v)
!
            call mem%dealloc(L_kcbj, batch_k%length, wf%n_v, wf%n_v, batch_j%length)
!
         enddo ! batch_k
      enddo ! batch_j
!
      req0 = 0
!
      req1_a = (wf%integrals%n_J)*(wf%n_o)
      req1_c = (wf%integrals%n_J)*(wf%n_o)
!
      req2 = 2*(wf%n_o)**2
!
      call batch_a%init(wf%n_v)
      call batch_c%init(wf%n_v)
!
      call mem%batch_setup(batch_a, batch_c, req0, req1_a, req1_c, req2)
!
      do current_a_batch = 1, batch_a%num_batches
!
         call batch_a%determine_limits(current_a_batch)
!
         do current_c_batch = 1, batch_c%num_batches
!
            call batch_c%determine_limits(current_c_batch)
!
!           u_aick = 2 t^ac_ik - t^ac_ki = - (2 g_aick - g_akci)/ε^{ac}_{ik}
!
            call mem%alloc(g_aick, batch_a%length, wf%n_o, batch_c%length, wf%n_o)
!
            call wf%get_vovo(g_aick,                     &
                           batch_a%first, batch_a%last,  &
                           1, wf%n_o,                    &
                           batch_c%first, batch_c%last,  &
                           1, wf%n_o)
!
            call mem%alloc(u_aikc, batch_a%length, wf%n_o, wf%n_o, batch_c%length)
!
!$omp parallel do private(c,k,i,a)
            do k = 1, wf%n_o
               do c = 1, batch_c%length
                  do i = 1, wf%n_o
                     do a = 1, batch_a%length

!
                        u_aikc(a,i,k,c) = - (two*g_aick(a,i,c,k)- g_aick(a,k,c,i))  &
                                          /(eps_v(a + batch_a%first - 1)            &
                                          + eps_v(c + batch_c%first - 1)            &
                                          - eps_o(i) - eps_o(k))
!
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
            call mem%dealloc(g_aick, batch_a%length, wf%n_o, batch_c%length, wf%n_o)
!
!           rho_ai = rho_ai + sum_ck u_aikc X_kc
!
            call mem%alloc(rho_ai_batch, batch_a%length, wf%n_o)
!
            kc_offset = wf%n_o*(batch_c%first - 1) + 1
!
            call dgemm('N', 'N',                   &
                        (batch_a%length)*(wf%n_o), &
                        1,                         &
                        (wf%n_o)*(batch_c%length), &
                        one,                       &
                        u_aikc,                    & ! u_ai_kc
                        (batch_a%length)*(wf%n_o), &
                        X_kc(kc_offset, 1),        & ! X_kc
                        (wf%n_o)*(wf%n_v),         &
                        zero,                      &
                        rho_ai_batch,              & ! rho_ai
                        (batch_a%length)*(wf%n_o))
!
            call mem%dealloc(u_aikc, batch_a%length, wf%n_o, wf%n_o, batch_c%length)
!
!$omp parallel do private(a, i)
            do i = 1, wf%n_o
               do a = 1, batch_a%length
!
                  rho_ai(a + batch_a%first - 1, i) = rho_ai(a + batch_a%first - 1, i) &
                                                 + rho_ai_batch(a, i)
!
               enddo
            enddo
!$omp end parallel do
!
            call mem%dealloc(rho_ai_batch, batch_a%length, wf%n_o)
!
         enddo ! batch_c
      enddo ! batch_a
!
      call mem%dealloc(X_kc, wf%n_o, wf%n_v)
!
!     :: Term 2: - L_kcjb t^cb_ki c_aj ::
!
      call mem%alloc(X_ji, wf%n_o, wf%n_o)
      call zero_array(X_ji, wf%n_o**2)
!
      req0 = 0
!
      req1_k = (wf%integrals%n_J)*(wf%n_v)
      req1_j = (wf%integrals%n_J)*(wf%n_v)
      req1_i = (wf%integrals%n_J)*(wf%n_v)
!
      req2_kj = 2*(wf%n_v**2)
      req2_ki = (wf%n_v**2)
      req2_ji = 0
!
      req3 = 0
!
      call batch_k%init(wf%n_o)
      call batch_j%init(wf%n_o)
      call batch_i%init(wf%n_o)
!
      call mem%batch_setup(batch_k, batch_j, batch_i, req0, req1_k, req1_j, req1_i, &
            req2_kj, req2_ki, req2_ji, req3)
!
      do current_k_batch = 1, batch_k%num_batches
!
         call batch_k%determine_limits(current_k_batch)
!
         do current_j_batch = 1, batch_j%num_batches
!
            call batch_j%determine_limits(current_j_batch)
!
               do current_i_batch = 1, batch_j%num_batches
!
                  call batch_i%determine_limits(current_i_batch)
!
!                 L_kcjb = 2 g_kcjb - g_jckb  (ordered as L_jckb)
!
                  call mem%alloc(g_jckb, batch_j%length, wf%n_v, batch_k%length, wf%n_v)
!
                  call wf%get_ovov(g_jckb,                     &
                                 batch_j%first, batch_j%last,  &
                                 1, wf%n_v,                    &
                                 batch_k%first, batch_k%last,  &
                                 1, wf%n_v)
!
                  call mem%alloc(L_jckb, batch_j%length, wf%n_v, batch_k%length, wf%n_v)
!
                  call dcopy(batch_j%length*(wf%n_v**2)*batch_k%length, g_jckb, 1, L_jckb, 1)
                  call dscal(batch_j%length*(wf%n_v**2)*batch_k%length, -one, L_jckb, 1)
!
                  call add_1432_to_1234(two, g_jckb, L_jckb, batch_j%length, wf%n_v, batch_k%length, wf%n_v)
                  call mem%dealloc(g_jckb, batch_j%length, wf%n_v, batch_k%length, wf%n_v)
!
!                 t_ckbi = - g_ckbi/ε^{cb}_{ik}
!
                  call mem%alloc(g_ckbi, wf%n_v, batch_k%length, wf%n_v, batch_i%length)
!
                  call wf%get_vovo(g_ckbi,                        &
                                    1, wf%n_v,                    &
                                    batch_k%first, batch_k%last,  &
                                    1, wf%n_v,                    &
                                    batch_i%first, batch_i%last)
!
!$omp parallel do private(i,b,k,c)
                  do c = 1, wf%n_v
                     do i = 1, batch_i%length
                        do k = 1, batch_k%length
                           do b = 1, wf%n_v
!
                              g_ckbi(c,k,b,i) = - g_ckbi(c,k,b,i) &
                                                /(eps_v(c) + eps_v(b) &
                                                - eps_o(i + batch_i%first - 1) &
                                                - eps_o(k + batch_k%first - 1))
                           enddo
                        enddo
                     enddo
                  enddo
!$omp end parallel do
!

            call dgemm('N', 'N',                                     &
                        batch_j%length,                              &
                        batch_i%length,                              &
                        (wf%n_v**2)*(batch_k%length),                &
                        one,                                         &
                        L_jckb,                                      & ! L_j_ckb
                        batch_j%length,                              &
                        g_ckbi,                                      & ! g_ckb_i
                        (wf%n_v**2)*(batch_k%length),                &
                        one,                                         &
                        X_ji(batch_j%first, batch_i%first),          & ! X_ji
                        (wf%n_o))
!
            call mem%dealloc(L_jckb, batch_j%length, wf%n_v, batch_k%length, wf%n_v)
            call mem%dealloc(g_ckbi, wf%n_v, batch_k%length, wf%n_v, batch_i%length)
!
            enddo ! batch_i
         enddo ! batch_j
      enddo ! batch_k
!
!     rho_ai = rho_ai - c_aj X_ji
!
      call dgemm('N', 'N',  &
                  (wf%n_v), &
                  (wf%n_o), &
                  (wf%n_o), &
                  -one,     &
                  c_bj,     & ! c_a_j
                  (wf%n_v), &
                  X_ji,     & ! X_j_i
                  (wf%n_o), &
                  one,      &
                  rho_ai,   & ! rho_a_i
                  (wf%n_v))
!
      call mem%dealloc(X_ji, wf%n_o, wf%n_o)
!
!     :: Term 3: - L_kcjb t^ca_kj c_bi ::
!
!     X_ab = t_akcj L_kcjb
!
      call mem%alloc(X_ab, (wf%n_v), (wf%n_v))
      call zero_array(X_ab, wf%n_v**2)
!
      req0 = 0
!
      req1_k = (wf%integrals%n_J)*(wf%n_v)
      req1_j = (wf%integrals%n_J)*(wf%n_v)
!
      req2 = 3*(wf%n_v)**2
!
      call batch_k%init(wf%n_o)
      call batch_j%init(wf%n_o)
!
      call mem%batch_setup(batch_k, batch_j, req0, req1_k, req1_j, req2)
!
      do current_k_batch = 1, batch_k%num_batches
!
         call batch_k%determine_limits(current_k_batch)
!
         do current_j_batch = 1, batch_j%num_batches
!
            call batch_j%determine_limits(current_j_batch)
!
!           L_kcjb = 2 g_kcjb - g_kbjc
!
            call mem%alloc(g_kcjb, batch_k%length, wf%n_v, batch_j%length, wf%n_v)
!
            call wf%get_ovov(g_kcjb,                        &
                              batch_k%first, batch_k%last,  &
                              1, wf%n_v,                    &
                              batch_j%first, batch_j%last,  &
                              1, wf%n_v)
!
            call mem%alloc(L_kcjb, batch_k%length, wf%n_v, batch_j%length, wf%n_v)
!
            call dcopy(batch_k%length*(wf%n_v**2)*batch_j%length, g_kcjb, 1, L_kcjb, 1)
            call dscal(batch_k%length*(wf%n_v**2)*batch_j%length, two, L_kcjb, 1)
!
            call add_1432_to_1234(-one, g_kcjb, L_kcjb, &
                                 batch_k%length, wf%n_v, batch_j%length, wf%n_v)
!
            call mem%dealloc(g_kcjb, batch_k%length, wf%n_v, batch_j%length, wf%n_v)
!
!           t_akcj = - g_ckaj/ε^{ca}_{jk}
!
            call mem%alloc(g_ckaj, wf%n_v, batch_k%length, wf%n_v, batch_j%length)
!
            call wf%get_vovo(g_ckaj,                        &
                              1, wf%n_v,                    &
                              batch_k%first, batch_k%last,  &
                              1, wf%n_v,                    &
                              batch_j%first, batch_j%last)
!
            call mem%alloc(t_akcj, wf%n_v, batch_k%length, wf%n_v, batch_j%length)
!
!$omp parallel do private(j,c,k,a)
            do c = 1, wf%n_v
               do j = 1, batch_j%length
                  do k = 1, batch_k%length
                     do a = 1, wf%n_v
!
                        t_akcj(a,k,c,j) = - g_ckaj(c,k,a,j) &
                                  /(eps_v(a) + eps_v(c) &
                                    - eps_o(j + batch_j%first - 1) &
                                    - eps_o(k + batch_k%first - 1))
!
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
            call mem%dealloc(g_ckaj, wf%n_v, batch_k%length, wf%n_v, batch_j%length)
!
            call dgemm('N', 'N',                                     &
                        (wf%n_v),                                    &
                        (wf%n_v),                                    &
                        (wf%n_v)*(batch_k%length)*(batch_j%length),  &
                        one,                                         &
                        t_akcj,                                      & ! t_a_kcj
                        (wf%n_v),                                    &
                        L_kcjb,                                      & ! L_kcj_b
                        (wf%n_v)*(batch_k%length)*(batch_j%length),  &
                        one,                                         &
                        X_ab,                                        & ! X_a_b
                        (wf%n_v))
!
            call mem%dealloc(t_akcj, wf%n_v, batch_k%length, wf%n_v, batch_j%length)
            call mem%dealloc(L_kcjb, batch_k%length, wf%n_v, batch_j%length, wf%n_v)
!
         enddo ! batch_k
      enddo ! batch_j
!
!     rho_ai = rho_ai - X_ab c_bi
!
      call dgemm('N', 'N',  &
                  (wf%n_v), &
                  (wf%n_o), &
                  (wf%n_v), &
                  -one,     &
                  X_ab,     & ! X_a_b
                  (wf%n_v), &
                  c_bj,     & ! c_b_i
                  (wf%n_v), &
                  one,      &
                  rho_ai ,  & ! rho_a_i
                  (wf%n_v))
!
      call mem%dealloc(X_ab, (wf%n_v), (wf%n_v))
!
   end subroutine jacobian_cc2_b1_lowmem_cc2
!
!
   module subroutine effective_jacobian_cc2_a1_lowmem_cc2(wf, omega, rho_ai, c_ai, eps_o, eps_v)
!!
!!    Jacobian CC2 A1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,
!!    Linda Goletto, and Alexander Paul, Dec 2018
!!
!!    rho_ai =+ F_kc * (-eps_ai,ck + w)^-1 * (2 g_aicd c_dk + 2 g_ckad c_di - g_akcd c_di - g_ciad c_dk)
!!           =+ F_kc * (-eps_ai,ck + w)^-1 * (2 X_aick - X_akci + 2 X_ckai - X_ciak)
!!           =+ F_kc * (Y_aick + Y_ckai)
!!
!!    The term is calculated in batches over the a and c indices.
!!
      implicit none
!
      class(lowmem_cc2), intent(in) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
      real(dp), dimension(:,:,:,:), allocatable :: X_aick
      real(dp), dimension(:,:,:,:), allocatable :: Y_aick
!
      real(dp), dimension(:,:,:,:), allocatable :: g_aicd
!
      real(dp), dimension(:,:), allocatable :: F_ck
      real(dp), dimension(:,:), allocatable :: reduced_rho_ai
!
      integer :: a, i, c, k
!
      type(batching_index) :: batch_a, batch_c
!
      integer :: current_a_batch, current_c_batch
      integer :: req0, req1_a, req1_c, req2
!
      req0   = 0
      req1_a = wf%integrals%n_J*wf%n_o
      req1_c = wf%integrals%n_J*wf%n_v
      req2   = wf%n_v*wf%n_o + wf%n_o**2
!
      call batch_a%init(wf%n_v)
      call batch_c%init(wf%n_v)
!
      call mem%batch_setup(batch_a, batch_c, req0, req1_a, req1_c, req2)
!
      do current_c_batch = 1, batch_c%num_batches
!
         call batch_c%determine_limits(current_c_batch)
!
         do current_a_batch = 1, batch_a%num_batches
!
            call batch_a%determine_limits(current_a_batch)
!
!           X_aick = sum_d g_aicd c_dk
!
            call mem%alloc(g_aicd, batch_a%length, wf%n_o, batch_c%length, wf%n_v)
!
            call wf%get_vovv(g_aicd,                       &
                              batch_a%first, batch_a%last, &
                              1, wf%n_o,                   &
                              batch_c%first, batch_c%last, &
                              1, wf%n_v)
!
            call mem%alloc(X_aick, batch_a%length, wf%n_o, batch_c%length, wf%n_o)
!
            call dgemm('N', 'N',                                  &
                        wf%n_o*(batch_a%length)*(batch_c%length), &
                        wf%n_o,                                   &
                        wf%n_v,                                   &
                        one,                                      &
                        g_aicd,                                   & ! g_aic_d
                        wf%n_o*(batch_a%length)*(batch_c%length), &
                        c_ai,                                     & ! c_d_k
                        wf%n_v,                                   &
                        zero,                                     &
                        X_aick,                                   & ! X_aic_k
                        wf%n_o*(batch_a%length)*(batch_c%length))
!
            call mem%dealloc(g_aicd, batch_a%length, wf%n_o, batch_c%length, wf%n_v)
!
!           Y_aick = (-eps_ai,ck + w)^-1 * (2 X_aick - X_akci)
!
            call mem%alloc(Y_aick, batch_a%length, wf%n_o, batch_c%length, wf%n_o)
!
!$omp parallel do private(k,c,i,a)
            do k = 1, wf%n_o
               do c = 1, batch_c%length
                  do i = 1, wf%n_o
                     do a = 1, batch_a%length
!
                        Y_aick(a, i, c, k) = (two*X_aick(a, i, c, k) - X_aick(a, k, c, i))/&
                              (- eps_v(a + batch_a%first - 1) - eps_v(c + batch_c%first - 1) &
                                 + eps_o(i) + eps_o(k) + omega)
!
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
            call mem%dealloc(X_aick, batch_a%length, wf%n_o, batch_c%length, wf%n_o)
!
!           Reorder occ-vir Fock matrix, then compute & add the term,
!           2 F_ck * (Y_aick + Y_ckai) to the transformed vector rho_ai
!
            call mem%alloc(F_ck, batch_c%length, wf%n_o)
!
!$omp parallel do private(c,k)
            do k = 1, wf%n_o
               do c = 1, batch_c%length
!
                  F_ck(c, k) = wf%fock_ia(k, c + batch_c%first - 1)
!
               enddo
            enddo
!$omp end parallel do
!
            call mem%alloc(reduced_rho_ai, batch_a%length, wf%n_o)
!
            call dgemm('N', 'N',                    &
                        wf%n_o*(batch_a%length),    &
                        1,                          &
                        wf%n_o*(batch_c%length),    &
                        one,                        &
                        Y_aick,                     & ! Y_ai_ck
                        wf%n_o*(batch_a%length),    &
                        F_ck,                       & ! F_ck
                        wf%n_o*(batch_c%length),    &
                        zero,                       &
                        reduced_rho_ai,             & ! reduced_rho_ai
                        wf%n_o*(batch_a%length))
!
!$omp parallel do private(i,a)
            do i = 1, wf%n_o
               do a = 1, batch_a%length
!
                  rho_ai(a + batch_a%first - 1, i) = rho_ai(a + batch_a%first - 1, i) &
                                                   + reduced_rho_ai(a, i)
!
               enddo
            enddo
!$omp end parallel do
!
            call mem%dealloc(reduced_rho_ai, batch_a%length, wf%n_o)
            call mem%dealloc(F_ck, batch_c%length, wf%n_o)
!
!           Now we pretend that ck is ai and vice versa s.t. Y_aick = Y_ckai
!
            call mem%alloc(F_ck, batch_a%length, wf%n_o)
!
!$omp parallel do private(c,k)
            do k = 1, wf%n_o
               do c = 1, batch_a%length
!
                  F_ck(c, k) = wf%fock_ia(k, c + batch_a%first - 1)
!
               enddo
            enddo
!$omp end parallel do
!
            call mem%alloc(reduced_rho_ai, batch_c%length, wf%n_o)
!
            call dgemm('T', 'N',                   &
                        wf%n_o*(batch_c%length),   &
                        1,                         &
                        wf%n_o*(batch_a%length),   &
                        one,                       &
                        Y_aick,                    & ! Pretend it is Y_ck_ai^T = Y_ai_ck
                        wf%n_o*(batch_a%length),   &
                        F_ck,                      & ! F_ck
                        wf%n_o*(batch_a%length),   &
                        zero,                      &
                        reduced_rho_ai,            &
                        wf%n_o*(batch_c%length))
!
            call mem%dealloc(Y_aick, batch_a%length, wf%n_o, batch_c%length, wf%n_o)
            call mem%dealloc(F_ck, batch_a%length, wf%n_o)
!
!
!$omp parallel do private(i,a)
            do i = 1, wf%n_o
               do a = 1, batch_c%length
!
                  rho_ai(a + batch_c%first - 1, i) = rho_ai(a + batch_c%first - 1, i) &
                                                   + reduced_rho_ai(a, i)
!
               enddo
            enddo
!$omp end parallel do
!
         call mem%dealloc(reduced_rho_ai, batch_c%length, wf%n_o)
!
         enddo ! batch_a
      enddo ! batch_c
!
   end subroutine effective_jacobian_cc2_a1_lowmem_cc2
!
!
   module subroutine effective_jacobian_cc2_b1_lowmem_cc2(wf, omega, rho_ai, c_ai, eps_o, eps_v)
!!
!!    Effective jacobian B1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!    Linda Goletto, and Alexander Paul, Jan 2019
!!
!!    Effective B1 = - 2 sum_{kcl} F_kc (1/(ε_{aick} + ω)) * (g_ailk c_cl + g_ckli c_al)
!!                     + sum_{kcl} F_kc (1/(ε_{akci} + ω)) * (g_akli c_cl + g_cilk c_al)
!!                 =   2 sum_{kcl} F_kc (- 2*X_ckai - 2*X_aick + X_ciak + X_akci)
!!
!!
      implicit none
!
      class(lowmem_cc2), intent(in) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
!
      real(dp), dimension(wf%n_o), intent(in)  :: eps_o
      real(dp), dimension(wf%n_v), intent(in)  :: eps_v
!
      real(dp), dimension(:,:,:,:), allocatable :: g_lkai, X_ckai
!
      integer :: a, c, i, k
!
      integer :: req0, req1_i, req1_k, req2
!
      integer :: current_i_batch, current_k_batch
!
      type(batching_index) :: batch_i, batch_k
!
      req0 = 0
!
      req1_i = (wf%integrals%n_J)*(wf%n_v)
      req1_k = (wf%integrals%n_J)*(wf%n_o)
!
      req2 = 2*(wf%n_o)*(wf%n_v)
!
      call batch_i%init(wf%n_o)
      call batch_k%init(wf%n_o)
!
      call mem%batch_setup(batch_i, batch_k, req0, req1_i, req1_k, req2)
!
      do current_i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(current_i_batch)
!
         do current_k_batch = 1, batch_k%num_batches
!
            call batch_k%determine_limits(current_k_batch)
!
!           X_c_kai = sum_l c_cl g_lkai
!
            call mem%alloc(g_lkai, wf%n_o, batch_k%length, wf%n_v, batch_i%length)
!
            call wf%get_oovo(g_lkai,                        &
                              1, wf%n_o,                    &
                              batch_k%first, batch_k%last,  &
                              1, wf%n_v,                    &
                              batch_i%first, batch_i%last)
!
            call mem%alloc(X_ckai, wf%n_v, batch_k%length, wf%n_v, batch_i%length)
!
            call dgemm('N', 'N',                                  &
                        wf%n_v,                                   &
                        (batch_i%length)*(batch_k%length)*wf%n_v, &
                        wf%n_o,                                   &
                        one,                                      &
                        c_ai,                                     & ! c_c_l
                        wf%n_v,                                   &
                        g_lkai,                                   & ! g_l_kai
                        wf%n_o,                                   &
                        zero,                                     &
                        X_ckai,                                   & ! X_c_kai
                        wf%n_v)
!
            call mem%dealloc(g_lkai, wf%n_o, batch_k%length, wf%n_v, batch_i%length)
!
!
!$omp parallel do private(a,i,k,c)
            do a = 1, wf%n_v
               do i = 1, batch_i%length
                  do k = 1, batch_k%length
                     do c = 1, wf%n_v
!
                        X_ckai(c, k, a, i) = X_ckai(c, k, a, i)/(eps_o(k + batch_k%first - 1) &
                                                               + eps_o(i + batch_i%first - 1) &
                                                               - eps_v(a) - eps_v(c) + omega)
!
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
!           Add all four terms while considering batched indices i and k
!
!$omp parallel do private(i,a,k,c)
            do a = 1, wf%n_v
               do i = 1, batch_i%length
                  do k = 1, batch_k%length
                     do c = 1, wf%n_v
!
                        rho_ai(a, i + batch_i%first - 1) = rho_ai(a, i + batch_i%first - 1)        &
                                 - two*X_ckai(c, k, a, i)*(wf%fock_ia(k + batch_k%first - 1, c))   &
                                 + X_ckai(a, k, c, i)*(wf%fock_ia(k + batch_k%first - 1, c))
!
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
!$omp parallel do private(a,k,i,c)
            do a = 1, wf%n_v
               do k = 1, batch_k%length
                  do i = 1, batch_i%length
                     do c = 1, wf%n_v
!
                        rho_ai(a, k + batch_k%first - 1) = rho_ai(a, k + batch_k%first - 1)        & !(k <-> i)
                                 - two*X_ckai(a, k, c, i)*(wf%fock_ia(i + batch_i%first - 1, c))   &
                                 + X_ckai(c, k, a, i)*(wf%fock_ia(i + batch_i%first - 1, c))
!
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
            call mem%dealloc(X_ckai, wf%n_v, batch_k%length, wf%n_v, batch_i%length)
!
         enddo ! batch_k
      enddo ! batch_i
!
   end subroutine effective_jacobian_cc2_b1_lowmem_cc2
!
!
   module subroutine effective_jacobian_cc2_c1_lowmem_cc2(wf, omega, rho_ai, c_cj, eps_o, eps_v)
!!
!!    Jacobian CC2 C1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,
!!    Linda Goletto, and Alexander Paul, Dec 2018
!!
!!    Implicit calculation of the doubles vector
!!    rho_ai^C1 =+ sum_ckbj - L_kijb  (g_akbc * c_cj + g_bjac * c_ck) (omega - ε_akbj)^-1
!!              =+ sum_kjb - L_kijb  (X_akbj + X_bjak)
!!              =+ sum_kjb - L_kijb Y_a_kjb
!!
      implicit none
!
      class(lowmem_cc2), intent(in) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_cj
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
      real(dp), dimension(:,:,:,:), allocatable :: g_akbc, g_bjac, g_kijb
      real(dp), dimension(:,:,:,:), allocatable :: L_kjbi
!
      real(dp), dimension(:,:,:,:), allocatable :: X_akbj, X_bjak, Y_akjb
!
      integer :: j, k, b, a
!
      type(batching_index) :: batch_i, batch_a, batch_b
!
      integer :: current_i_batch, current_a_batch, current_b_batch
      integer :: req0, req1_i, req1_a, req1_b, req2_ia, req2_ib, req2_ab, req3
!
      req0 = 0
!
      req1_i = (wf%integrals%n_J)*(wf%n_o)
      req1_a = max((wf%integrals%n_J)*(wf%n_o), (wf%integrals%n_J)*(wf%n_v))
      req1_b = max((wf%integrals%n_J)*(wf%n_o), (wf%integrals%n_J)*(wf%n_v))
!
      req2_ia = 0
      req2_ib = 2*(wf%n_o)**2
      req2_ab = max(3*(wf%n_o)**2, 2*(wf%n_o)**2 + (wf%n_o)*(wf%n_v))
!
      req3 = 0
!
      call batch_i%init(wf%n_o)
      call batch_a%init(wf%n_v)
      call batch_b%init(wf%n_v)
!
      call mem%batch_setup(batch_i, batch_a, batch_b, req0, req1_i, req1_a, req1_b, &
                           req2_ia, req2_ib, req2_ab, req3)
!
      do current_i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(current_i_batch)
!
         do current_a_batch = 1, batch_a%num_batches
!
               call batch_a%determine_limits(current_a_batch)
!
            do current_b_batch = 1, batch_b%num_batches
!
               call batch_b%determine_limits(current_b_batch)
!
!              X_akbj = sum_c g_akbc * c_cj
!
               call mem%alloc(g_akbc, batch_a%length, wf%n_o, batch_b%length, wf%n_v)
!
               call wf%get_vovv(g_akbc,                        &
                                 batch_a%first, batch_a%last,  &
                                 1, wf%n_o,                    &
                                 batch_b%first, batch_b%last,  &
                                 1, wf%n_v)
!
               call mem%alloc(X_akbj, batch_a%length, wf%n_o, batch_b%length, wf%n_o)
!
               call dgemm('N', 'N',                                     &
                           (batch_a%length)*(wf%n_o)*(batch_b%length),  &
                           (wf%n_o),                                    &
                           (wf%n_v),                                    &
                           one,                                         &
                           g_akbc,                                      & ! g_akb_c
                           (batch_a%length)*(wf%n_o)*(batch_b%length),  &
                           c_cj,                                        & ! c_c_j
                           (wf%n_v),                                    &
                           zero,                                        &
                           X_akbj,                                      & ! X_akb_j
                           (batch_a%length)*(wf%n_o)*(batch_b%length))
!
               call mem%dealloc(g_akbc, batch_a%length, wf%n_o, batch_b%length, wf%n_v)
!
!              X_bjak = sum_c g_bjac * c_ck
!
               call mem%alloc(g_bjac, batch_b%length, wf%n_o, batch_a%length, wf%n_v)
!
               call wf%get_vovv(g_bjac,                        &
                                 batch_b%first, batch_b%last,  &
                                 1, wf%n_o,                    &
                                 batch_a%first, batch_a%last,  &
                                 1, wf%n_v)
!
               call mem%alloc(X_bjak, batch_b%length, wf%n_o, batch_a%length, wf%n_o)
!
               call dgemm('N', 'N',                                     &
                           (batch_b%length)*(wf%n_o)*(batch_a%length),  &
                           (wf%n_o),                                    &
                           (wf%n_v),                                    &
                           one,                                         &
                           g_bjac,                                      & ! g_bja_c
                           (batch_b%length)*(wf%n_o)*(batch_a%length),  &
                           c_cj,                                        & ! c_c_k
                           (wf%n_v),                                    &
                           zero,                                        &
                           X_bjak,                                      & ! X_bja_k
                           (batch_b%length)*(wf%n_o)*(batch_a%length))
!
               call mem%dealloc(g_bjac, batch_b%length, wf%n_o, batch_a%length, wf%n_v)
!
               call mem%alloc(Y_akjb, batch_a%length, wf%n_o, wf%n_o, batch_b%length)
!
!              Y_akjb = (X_akbj + X_bjak)  /(ε_akbj + omega)
!
!$omp parallel do private(j,b,k,a)
               do j = 1, wf%n_o
                  do b = 1, batch_b%length
                     do k = 1, wf%n_o
                        do a = 1, batch_a%length
!
                           Y_akjb(a,k,j,b) = (X_akbj(a,k,b,j) + X_bjak(b,j,a,k)) &
                                             /(- eps_v(a + batch_a%first - 1)    &
                                             - eps_v(b + batch_b%first - 1)      &
                                             + eps_o(j) + eps_o(k) + omega)
!
                        enddo
                     enddo
                  enddo
               enddo
!$omp end parallel do
!
               call mem%dealloc(X_akbj, batch_a%length, wf%n_o, batch_b%length, wf%n_o)
               call mem%dealloc(X_bjak, batch_b%length, wf%n_o, batch_a%length, wf%n_o)
!
!              L_kijb = 2 g_kijb - g_jikb, ordered as L_kjbi
!
               call mem%alloc(g_kijb, wf%n_o, batch_i%length, wf%n_o, batch_b%length)
!
               call wf%get_ooov(g_kijb,                        &
                                 1, wf%n_o,                    &
                                 batch_i%first, batch_i%last,  &
                                 1, wf%n_o,                    &
                                 batch_b%first, batch_b%last)
!
               call mem%alloc(L_kjbi, wf%n_o, wf%n_o, batch_b%length, batch_i%length)
!
               call zero_array(L_kjbi, wf%n_v*wf%n_o**3)
!
               call add_1423_to_1234(two, g_kijb, L_kjbi, wf%n_o, wf%n_o, batch_b%length, batch_i%length)
               call add_2413_to_1234(-one, g_kijb, L_kjbi, wf%n_o, wf%n_o, batch_b%length, batch_i%length)
!
               call mem%dealloc(g_kijb, wf%n_o, batch_i%length, wf%n_o, batch_b%length)
!
!              rho_ai = rho_ai - sum_jbk Y_akjb * L_kjbi
!
               call dgemm('N', 'N',                               &
                           (batch_a%length),                      &
                           (batch_i%length),                      &
                           (batch_b%length)*(wf%n_o)**2,          &
                           -one,                                  &
                           Y_akjb,                                & ! Y_a_kjb
                           (batch_a%length),                      &
                           L_kjbi,                                & ! L_kjb_i
                           (batch_b%length)*(wf%n_o)**2,          &
                           one,                                  &
                           rho_ai(batch_a%first, batch_i%first),  & ! rho_ai
                           (wf%n_v))
!
               call mem%dealloc(Y_akjb, batch_a%length, wf%n_o, wf%n_o, batch_b%length)
               call mem%dealloc(L_kjbi, wf%n_o, wf%n_o, batch_b%length, batch_i%length)
!
            enddo ! batch_b
         enddo ! batch_a
      enddo ! batch_i
!
!
   end subroutine effective_jacobian_cc2_c1_lowmem_cc2
!
!
   module subroutine effective_jacobian_cc2_d1_lowmem_cc2(wf, omega, rho_ai, c_bl, eps_o, eps_v)
!!
!!    Jacobian CC2 D1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,
!!    Linda Goletto, and Alexander Paul, Dec 2018
!!
!!    Implicit calculation of the doubles vector
!!    rho_ai^D1 =+ sum_ckbj - L_kijb  (- g_aklj * c_bl - g_bjlk * c_al) (omega - ε_akbj)^-1 
!!              =+ sum_kjb L_kijb  (X_bjak + X_akbj)
!!              =+ sum_kjb L_kijb  Y_ajbk
!!
      implicit none
!
      class(lowmem_cc2), intent(in) :: wf
!
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_bl
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ljak, g_jbki, g_jikb, g_lkbj
      real(dp), dimension(:,:,:,:), allocatable :: L_jbki
!
      real(dp), dimension(:,:,:,:), allocatable :: X_akbj, X_bjak, Y_ajbk
!
      integer :: j, k, a, b
!
      integer :: req0, req1_j, req1_k, req2
!
      integer :: current_j_batch, current_k_batch
!
      type(batching_index) :: batch_j, batch_k
!
!     X_bjak = sum_l g_aklj * c_bl
!
      req0 = 0
      req1_j = max((wf%integrals%n_J)*(wf%n_v), (wf%integrals%n_J)*(wf%n_o))
      req1_k = max((wf%integrals%n_J)*(wf%n_v), (wf%integrals%n_J)*(wf%n_o))
      req2 = max((wf%n_v)**2 + 2*(wf%n_v)*(wf%n_o), 3*(wf%n_v)**2)
!
      call batch_j%init(wf%n_o)
      call batch_k%init(wf%n_o)
!
      call mem%batch_setup(batch_j, batch_k, req0, req1_j, req1_k, req2)
!
      do current_j_batch = 1, batch_j%num_batches
!
         call batch_j%determine_limits(current_j_batch)
!
         do current_k_batch = 1, batch_k%num_batches
!
            call batch_k%determine_limits(current_k_batch)
!
            call mem%alloc(g_ljak, wf%n_o, (batch_j%length), wf%n_v, (batch_k%length))
!
            call wf%get_oovo(g_ljak,                        &
                              1, wf%n_o,                    &
                              batch_j%first, batch_j%last,  &
                              1, wf%n_v,                    &
                              batch_k%first, batch_k%last)
!
            call mem%alloc(X_bjak, wf%n_v, (batch_j%length), wf%n_v, (batch_k%length))
!
            call dgemm('N', 'N',                                     &
                        (wf%n_v),                                    &
                        (batch_j%length)*(wf%n_v)*(batch_k%length),  &
                        (wf%n_o),                                    &
                        one,                                         &
                        c_bl,                                        & ! c_b_l
                        (wf%n_v),                                    &
                        g_ljak,                                      & ! g_l_jak
                        (wf%n_o),                                    &
                        zero,                                        &
                        X_bjak,                                      & ! X_b_jak
                        (wf%n_v))
!
            call mem%dealloc(g_ljak, wf%n_o, (batch_j%length), wf%n_v, (batch_k%length))
!
!           X_akbj = sum_c g_bjlk * c_al
!
            call mem%alloc(g_lkbj, wf%n_o, batch_k%length, wf%n_v, batch_j%length)
!
            call wf%get_oovo(g_lkbj,                        &
                              1, wf%n_o,                    &
                              batch_k%first, batch_k%last,  &
                              1, wf%n_v,                    &
                              batch_j%first, batch_j%last)
!
            call mem%alloc(X_akbj, wf%n_v, batch_k%length, wf%n_v, batch_j%length)
!
            call dgemm('N', 'N',                                     &
                        (wf%n_v),                                    &
                        (batch_k%length)*(batch_j%length)*(wf%n_v),  &
                        (wf%n_o),                                    &
                        one,                                         &
                        c_bl,                                        & ! c_al
                        (wf%n_v),                                    &
                        g_lkbj,                                      & ! g_l_kbj
                        (wf%n_o),                                    &
                        zero,                                        &
                        X_akbj,                                      & ! X_a_kbj
                        (wf%n_v))
!
            call mem%dealloc(g_lkbj, wf%n_o, batch_k%length, wf%n_v, batch_j%length)
!
!           Y_ajbk = (X_bjak + X_akbj)/(ε_bjak + ω)
!
            call mem%alloc(Y_ajbk, wf%n_v, (batch_j%length), wf%n_v, (batch_k%length))
!
!$omp parallel do private(k,b,j,a)
            do b = 1, wf%n_v
               do k = 1, batch_k%length
                  do j = 1, batch_j%length
                     do a = 1, wf%n_v
!
                        Y_ajbk(a,j,b,k) = (X_bjak(b,j,a,k) + X_akbj(a,k,b,j))  &
                                       /(omega - eps_v(a) - eps_v(b)          &
                                       + eps_o(j + batch_j%first - 1)         &
                                       + eps_o(k + batch_k%first - 1))
!
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
            call mem%dealloc(X_bjak, wf%n_v, batch_j%length, wf%n_v, batch_k%length)
            call mem%dealloc(X_akbj, wf%n_v, batch_k%length, wf%n_v, batch_j%length)
!
!           L_jbki = 2 g_kijb - g_kbji
!
            call mem%alloc(g_jbki,batch_j%length, wf%n_v, batch_k%length, wf%n_o)
!
            call wf%get_ovoo(g_jbki,                        &
                              batch_j%first, batch_j%last,  &
                              1, wf%n_v,                    &
                              batch_k%first, batch_k%last,  &
                              1, wf%n_o)
!
            call mem%alloc(L_jbki, batch_j%length, wf%n_v, batch_k%length, wf%n_o)
!
            call dcopy(batch_j%length*wf%n_v*batch_k%length*wf%n_o, g_jbki, 1, L_jbki, 1)
            call dscal(batch_j%length*wf%n_v*batch_k%length*wf%n_o, two, L_jbki, 1)
!
            call mem%dealloc(g_jbki, batch_j%length, wf%n_v, batch_k%length, wf%n_o)
!
            call mem%alloc(g_jikb, batch_j%length, wf%n_o, batch_k%length, wf%n_v)
!
            call wf%get_ooov(g_jikb,                        &
                              batch_j%first, batch_j%last,  &
                              1, wf%n_o,                    &
                              batch_k%first, batch_k%last,  &
                              1, wf%n_v)
!
            call add_1432_to_1234(-one, g_jikb, L_jbki, batch_j%length, wf%n_v, batch_k%length, wf%n_o)
!
            call mem%dealloc(g_jikb, batch_j%length, wf%n_o, batch_k%length, wf%n_v)
!
!           rho_ai = rho_ai + sum_jbk Y_ajbk * L_jbki
!
            call dgemm('N', 'N',                                     &
                        (wf%n_v),                                    &
                        (wf%n_o),                                    &
                        (wf%n_v)*(batch_j%length)*(batch_k%length),  &
                        one,                                         &
                        Y_ajbk,                                      & ! X_a_jbk
                        (wf%n_v),                                    &
                        L_jbki,                                      & ! L_jbk_i
                        (wf%n_v)*(batch_j%length)*(batch_k%length),  &
                        one,                                         &
                        rho_ai,                                      & ! rho_ai
                        (wf%n_v))
!
            call mem%dealloc(Y_ajbk, wf%n_v, batch_j%length, wf%n_v, batch_k%length)
            call mem%dealloc(L_jbki, batch_j%length, wf%n_v, batch_k%length, wf%n_o)
!
         enddo ! batch k
      enddo ! batch j
!
   end subroutine effective_jacobian_cc2_d1_lowmem_cc2
!
!
   module subroutine effective_jacobian_cc2_e1_lowmem_cc2(wf, omega, rho_ai, c_dk, eps_o, eps_v)
!!
!!    Jacobian CC2 E1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,
!!    Linda Goletto, and Alexander Paul, Jan 2019
!!
!!    Implicit calculation of the doubles vector
!!    rho_ai^E1 =+ sum_bckd L_abkc  (g_bicd * c_dk + g_ckbd * c_di) (omega - ε_bick)^-1
!!              =+ sum_bck L_abkc  (X_bick + X_ckbi)
!!
!!
      implicit none
!
      class(lowmem_cc2), intent(in) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_dk
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
      real(dp), dimension(:,:,:,:), allocatable :: g_bicd, g_ckbd, g_abkc, g_kbac
      real(dp), dimension(:,:,:,:), allocatable :: L_abkc
!
      real(dp), dimension(:,:,:,:), allocatable :: X_bick, X_ckbi, Y_bkci
!
      integer :: b, c, i, k
!
      integer :: req0, req1_b, req1_c, req2
!
      integer :: current_b_batch, current_c_batch
!
      type(batching_index) :: batch_b, batch_c
!
!
      req0 = 0
      req1_b = max((wf%integrals%n_J)*(wf%n_v),(wf%integrals%n_J)*(wf%n_o))
      req1_c = max((wf%integrals%n_J)*(wf%n_v),(wf%integrals%n_J)*(wf%n_o))
      req2 = max((wf%n_o)**2 + 2*(wf%n_v)*(wf%n_o), 3*(wf%n_o)**2)
!
      call batch_b%init(wf%n_v)
      call batch_c%init(wf%n_v)
!
      call mem%batch_setup(batch_b, batch_c, req0, req1_b, req1_c, req2)
!
      do current_b_batch = 1, batch_b%num_batches
!
         call batch_b%determine_limits(current_b_batch)
!
         do current_c_batch = 1, batch_c%num_batches
!
            call batch_c%determine_limits(current_c_batch)
!
!           X_bick = sum_d g_bicd * c_dk
!
            call mem%alloc(g_bicd, batch_b%length, wf%n_o, batch_c%length, wf%n_v)
!
            call wf%get_vovv(g_bicd,                        &
                              batch_b%first, batch_b%last,  &
                              1, wf%n_o,                    &
                              batch_c%first, batch_c%last,  &
                              1, wf%n_v)
!
            call mem%alloc(X_bick, batch_b%length, wf%n_o, batch_c%length, wf%n_o)
!
            call dgemm('N', 'N',                                     &
                        (batch_b%length)*(wf%n_o)*(batch_c%length),  &
                        (wf%n_o),                                    &
                        (wf%n_v),                                    &
                        one,                                         &
                        g_bicd,                                      & ! g_bic_d
                        (batch_b%length)*(wf%n_o)*(batch_c%length),  &
                        c_dk,                                        & ! c_d_k
                        (wf%n_v),                                    &
                        zero,                                        &
                        X_bick,                                      & ! X_bic_k
                        (batch_b%length)*(wf%n_o)*(batch_c%length))
!
            call mem%dealloc(g_bicd, batch_b%length, wf%n_o, batch_c%length, wf%n_v)
!
!           X_ckbi = sum_d g_ckbd * c_di
!
            call mem%alloc(g_ckbd, batch_c%length, wf%n_o, batch_b%length, wf%n_v)
!
            call wf%get_vovv(g_ckbd,                        &
                              batch_c%first, batch_c%last,  &
                              1, wf%n_o,                    &
                              batch_b%first, batch_b%last,  &
                              1, wf%n_v)
!
            call mem%alloc(X_ckbi, batch_c%length, wf%n_o, batch_b%length, wf%n_o)
!
            call dgemm('N', 'N',                                     &
                        (batch_c%length)*(wf%n_o)*(batch_b%length),  &
                        (wf%n_o),                                    &
                        (wf%n_v),                                    &
                        one,                                         &
                        g_ckbd,                                      & ! g_ckb_d
                        (batch_c%length)*(wf%n_o)*(batch_b%length),  &
                        c_dk,                                        & ! c_d_i
                        (wf%n_v),                                    &
                        zero,                                        &
                        X_ckbi,                                      & ! X_ckb_i
                        (batch_c%length)*(wf%n_o)*(batch_b%length))
!
            call mem%dealloc(g_ckbd, batch_c%length, wf%n_o, batch_b%length, wf%n_v)
!
!          Y_bkci = (X_bick + X_ckbi)/(ε_bick + ω)
!
            call mem%alloc(Y_bkci, batch_b%length, wf%n_o, batch_c%length, wf%n_o)
!
!$omp parallel do private(i,c,k,b)
            do i = 1, wf%n_o
               do c = 1, batch_c%length
                  do k = 1, wf%n_o
                     do b = 1, batch_b%length
!
                        Y_bkci(b,k,c,i) = (X_bick(b,i,c,k) + X_ckbi(c,k,b,i)) &
                                       /(omega + eps_o(i) + eps_o(k)          &
                                       - eps_v(b + batch_b%first - 1)         &
                                       - eps_v(c + batch_c%first - 1))
!
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
            call mem%dealloc(X_bick, batch_b%length, wf%n_o, batch_c%length, wf%n_o)
            call mem%dealloc(X_ckbi, batch_c%length, wf%n_o, batch_b%length, wf%n_o)
!
!           L_abkc = 2 g_abkc - g_kbac
!
            call mem%alloc(g_abkc, wf%n_v, batch_b%length, wf%n_o, batch_c%length)
!
            call wf%get_vvov(g_abkc,                        &
                              1, wf%n_v,                    &
                              batch_b%first, batch_b%last,  &
                              1, wf%n_o,                    &
                              batch_c%first, batch_c%last)
!
            call mem%alloc(L_abkc, wf%n_v, batch_b%length, wf%n_o, batch_c%length)
!
            call dcopy(wf%n_v*(batch_b%length)*(wf%n_o)*(batch_c%length), g_abkc, 1, L_abkc, 1)
!
            call mem%dealloc(g_abkc, wf%n_v, batch_b%length, wf%n_o, batch_c%length)
!
            call dscal(wf%n_v*(batch_b%length)*(wf%n_o)*(batch_c%length), two, L_abkc, 1)
!
            call mem%alloc(g_kbac, wf%n_o, batch_b%length, wf%n_v, batch_c%length)
!
            call wf%get_ovvv(g_kbac,                        &
                              1, wf%n_o,                    &
                              batch_b%first, batch_b%last,  &
                              1, wf%n_v,                    &
                              batch_c%first, batch_c%last)
!
            call add_3214_to_1234(-one, g_kbac, L_abkc, &
                                 wf%n_v, batch_b%length, wf%n_o, batch_c%length)
!
            call mem%dealloc(g_kbac, wf%n_o, batch_b%length, wf%n_v, batch_c%length)
!
!           rho_ai = rho_ai + sum_bkc L_abkc * Y_bkci
!
            call dgemm('N', 'N',                                     &
                        (wf%n_v),                                    &
                        (wf%n_o),                                    &
                        (wf%n_o)*(batch_b%length)*(batch_c%length),  &
                        one,                                         &
                        L_abkc,                                      & ! L_a_bkc
                        (wf%n_v),                                    &
                        Y_bkci,                                      & ! Y_bkc_i
                        (wf%n_o)*(batch_b%length)*(batch_c%length),  &
                        one,                                         &
                        rho_ai,                                      & ! rho_ai
                        (wf%n_v))
!
            call mem%dealloc(Y_bkci, batch_b%length, wf%n_o, batch_c%length, wf%n_o)
            call mem%dealloc(L_abkc, wf%n_v, batch_b%length, wf%n_o, batch_c%length)
!
         enddo ! batch_c
      enddo ! batch_b
!
   end subroutine effective_jacobian_cc2_e1_lowmem_cc2
!
!
   module subroutine effective_jacobian_cc2_f1_lowmem_cc2(wf, omega, rho_ai, c_ai, eps_o, eps_v)
!!
!!    Jacobian CC2 effective F1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,
!!    Linda Goletto, and Alexander Paul, Jan 2019
!!
!!    Effective F1 = - L_abkc (1/(ε_{bick} + ω) * (g_lkbi c_cl + g_lick c_bl))
!!
      implicit none
!
      class(lowmem_cc2), intent(in) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
      real(dp), dimension(:,:,:,:), allocatable :: X_ckbi, X_bick, Y_bcki
      real(dp), dimension(:,:,:,:), allocatable :: L_abck, g_abkc, g_lkbi, g_lick
!
      integer :: b, c, i, k
!
      type(batching_index) :: batch_i, batch_k, batch_a
!
      integer :: current_i_batch, current_k_batch, current_a_batch
      integer :: req0, req1_i, req1_k, req1_a, req2_ik, req2_ia, req2_ka, req3
!
      req0 = 0
!
      req1_i = max((wf%integrals%n_J)*(wf%n_o), (wf%integrals%n_J)*(wf%n_v))
      req1_k = max((wf%integrals%n_J)*(wf%n_o), (wf%integrals%n_J)*(wf%n_v))
      req1_a = (wf%integrals%n_J)*(wf%n_v)
!
      req2_ik = max(3*(wf%n_v)**2, 2*(wf%n_v)**2 + (wf%n_o)*(wf%n_v))
      req2_ia = 0
      req2_ka = 2*(wf%n_v)**2
!
      req3 = 0
!
      call batch_i%init(wf%n_o)
      call batch_k%init(wf%n_o)
      call batch_a%init(wf%n_v)
!
      call mem%batch_setup(batch_i, batch_k, batch_a, req0, req1_i, req1_k, req1_a, &
                           req2_ik, req2_ia, req2_ka, req3)
!
      do current_i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(current_i_batch)
!
         do current_k_batch = 1, batch_k%num_batches
!
            call batch_k%determine_limits(current_k_batch)
!
            do current_a_batch = 1, batch_a%num_batches
!
               call batch_a%determine_limits(current_a_batch)
!
!              X_ckbi = sum_c g_lkbi * c_cl
!
               call mem%alloc(g_lkbi, wf%n_o, batch_k%length, wf%n_v, batch_i%length)
!
               call wf%get_oovo(g_lkbi,                        &
                                 1, wf%n_o,                    &
                                 batch_k%first, batch_k%last,  &
                                 1, wf%n_v,                    &
                                 batch_i%first, batch_i%last)
!
               call mem%alloc(X_ckbi, wf%n_v, batch_k%length, wf%n_v, batch_i%length)
!
               call dgemm('N', 'N',                                     &
                           wf%n_v,                                      &
                           (batch_k%length)*(wf%n_v)*(batch_i%length),  &
                           wf%n_o,                                      &
                           one,                                         &
                           c_ai,                                        & ! c_c_l
                           wf%n_v,                                      &
                           g_lkbi,                                      & ! g_l_kbi
                           wf%n_o,                                      &
                           zero,                                        &
                           X_ckbi,                                      & ! X_c_kbi
                           wf%n_v)
!
               call mem%dealloc(g_lkbi, wf%n_o, batch_k%length, wf%n_v, batch_i%length)
!
!              X_bick = sum_l c_bl * g_lick
!
               call mem%alloc(g_lick, wf%n_o, batch_i%length, wf%n_v, batch_k%length)
!
               call wf%get_oovo(g_lick,                        &
                                 1, wf%n_o,                    &
                                 batch_i%first, batch_i%last,  &
                                 1, wf%n_v,                    &
                                 batch_k%first, batch_k%last)
!
               call mem%alloc(X_bick, wf%n_v, batch_i%length, wf%n_v, batch_k%length)
!
               call dgemm('N', 'N',                                     &
                           wf%n_v,                                      &
                           (batch_i%length)*(wf%n_v)*(batch_k%length),  &
                           wf%n_o,                                      &
                           one,                                         &
                           c_ai,                                        & ! c_b_l
                           wf%n_v,                                      &
                           g_lick,                                      & ! g_l_ick
                           wf%n_o,                                      &
                           zero,                                        &
                           X_bick,                                      & ! X_b_ick
                           wf%n_v)
!
               call mem%dealloc(g_lick, wf%n_o, batch_i%length, wf%n_v, batch_k%length)
!
!              Reorder and scale
!
               call mem%alloc(Y_bcki, wf%n_v, wf%n_v, batch_k%length, batch_i%length)
!
!              Y_bcki = (X_ckbi + X_bick)/(ε_ckbi + ω)
!
!$omp parallel do private(c,i,k,b)
               do c = 1, wf%n_v
                  do i = 1, batch_i%length
                     do k = 1, batch_k%length
                        do b = 1, wf%n_v
!
                           Y_bcki(b,c,k,i) = (X_ckbi(c,k,b,i) + X_bick(b,i,c,k)) &
                                             /(omega - eps_v(b) - eps_v(c)       &
                                             + eps_o(i + batch_i%first - 1)      &
                                             + eps_o(k + batch_k%first - 1))
!
                        enddo
                     enddo
                  enddo
               enddo
!$omp end parallel do
!
               call mem%dealloc(X_ckbi, wf%n_v, batch_k%length, wf%n_v, batch_i%length)
               call mem%dealloc(X_bick, wf%n_v, batch_i%length, wf%n_v, batch_k%length)
!
!              L_abkc = 2 g_abkc - g_ackb ordered as L_abck
!
               call mem%alloc(g_abkc, batch_a%length, wf%n_v, batch_k%length, wf%n_v)
!
               call wf%get_vvov(g_abkc,                        &
                                 batch_a%first, batch_a%last,  &
                                 1, wf%n_v,                    &
                                 batch_k%first, batch_k%last,  &
                                 1, wf%n_v)
!
               call mem%alloc(L_abck, batch_a%length, wf%n_v, wf%n_v, batch_k%length)
!
               call zero_array(L_abck, wf%n_o*wf%n_v**3)
!
               call add_1243_to_1234(two, g_abkc, L_abck, batch_a%length, wf%n_v, wf%n_v, batch_k%length)
               call add_1342_to_1234(-one, g_abkc, L_abck, batch_a%length, wf%n_v, wf%n_v, batch_k%length)
!
               call mem%dealloc(g_abkc, batch_a%length, wf%n_v, batch_k%length, wf%n_v)
!
!              rho_ai = rho_ai - sum_bkc L_kjbi * Y_akjb
!
               call dgemm('N', 'N',                               &
                           (batch_a%length),                      &
                           (batch_i%length),                      &
                           (batch_k%length)*(wf%n_v)**2,          &
                           -one,                                  &
                           L_abck,                                & ! L_a_bck
                           (batch_a%length),                      &
                           Y_bcki,                                & ! Y_bck_i
                           (batch_k%length)*(wf%n_v)**2,          &
                           one,                                   &
                           rho_ai(batch_a%first, batch_i%first),  & ! rho_ai
                           (wf%n_v))
!
               call mem%dealloc(Y_bcki, wf%n_v, wf%n_v, batch_k%length, batch_i%length)
               call mem%dealloc(L_abck, batch_a%length, wf%n_v, wf%n_v, batch_k%length)
!
            enddo ! batch_a
         enddo ! batch_k
      enddo ! batch_i
!
   end subroutine effective_jacobian_cc2_f1_lowmem_cc2
!
!
end submodule jacobian
