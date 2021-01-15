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
submodule (ccs_class) jacobian_ccs
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
!!    A_μ,ν = < μ | exp(-T) [H, τ_ν] exp(T) | R >.
!!
!
   implicit none
!
!
contains
!
!
   module subroutine prepare_for_jacobian_ccs(wf)
!!
!!    Prepare for jacobian
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2019
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
!     For now, do nothing.
!
      call output%printf('v', '- No preparations for the ' // trim(wf%name_) // &
                         ' excited state equation.', fs='(/t3,a)')
!
   end subroutine prepare_for_jacobian_ccs
!
!
   module subroutine jacobian_transformation_ccs(wf, c, rho)
!!
!!    Jacobian transformation
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!    Directs the transformation by the CCS Jacobi matrix,
!!
!!       A_μ,ν = < μ | exp(-T) [H, τ_ν] exp(T) | R >.
!!
!!    In particular,
!!
!!       rho_mu = (A c)_mu = sum_ck A_mu,ck c_ck.
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(wf%n_t1), intent(in) :: c
      real(dp), dimension(wf%n_t1), intent(out) :: rho
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian CCS transformation', pl='normal')
      call timer%turn_on()
!
      call zero_array(rho, wf%n_t1)
!
      call wf%jacobian_ccs_a1(rho, c)
      call wf%jacobian_ccs_b1(rho, c)
!
      call timer%turn_off()
!
   end subroutine jacobian_transformation_ccs
!
!
   module subroutine jacobian_ccs_a1_ccs(wf, rho1, c1)
!!
!!    Jacobian CCS A1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!    Calculates the A1 term,
!!
!!       sum_b F_ab c_bi - sum_j F_ji c_aj,
!!
!!    and adds it to the rho vector.
!!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c1
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho1
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian CCS A1', pl='verbose')
      call timer%turn_on()
!
!     sum_b F_a_b c_b_i 
!
      call dgemm('N', 'N',     &
                  wf%n_v,      &
                  wf%n_o,      &
                  wf%n_v,      &
                  one,         &
                  wf%fock_ab,  &
                  wf%n_v,      &
                  c1,          &
                  wf%n_v,      &
                  one,         &
                  rho1,        &
                  wf%n_v)
!
!     - sum_j c_a_j F_j_i
!
      call dgemm('N','N',      &
                  wf%n_v,      &
                  wf%n_o,      &
                  wf%n_o,      &
                  -one,        &
                  c1,          &
                  wf%n_v,      &
                  wf%fock_ij,  &
                  wf%n_o,      &
                  one,         &
                  rho1,        &
                  wf%n_v)
!
      
      call timer%turn_off()
!
   end subroutine jacobian_ccs_a1_ccs
!
!
   module subroutine jacobian_ccs_b1_ccs(wf, rho_ai, c_bj)
!!
!!    Jacobian CCS B1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,
!!    Linda Goletto, and Alexander C. Paul, Dec 2018
!!
!!    Calculates the A1 term
!!
!!       A1: sum_bj (2 g_aijb - g_abji) * c_bj
!!
!!    and adds it to rho_ai.
!!
!!    Separate calculation of both terms due to batching.
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_bj
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
!
      real(dp), dimension(:,:), allocatable :: rho_ai_batch
      real(dp), dimension(:), allocatable :: c_jb
!
      real(dp), dimension(:,:,:), allocatable :: L_J_jb, L_J_ai, L_J_ab, L_J_ji, X_Jaj, X_aJj
      real(dp), dimension(:), allocatable :: X_J
!
      integer :: req0, req1_i, req1_j, req1_a, req2
!
      integer :: current_i_batch, current_a_batch, current_j_batch
!
      integer :: a, i, j, b, jb
!
      type(batching_index) :: batch_i, batch_a, batch_j
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian CCS B1', pl='verbose')
      call timer%turn_on()
!
!     :: Term 1: rho_ai = sum_bj 2 g_aijb * c_bj ::
!
      req0 = 0
!
      req1_j = wf%n_v*wf%eri%n_J + wf%n_v
!
      batch_j = batching_index(wf%n_o)
!
      call mem%batch_setup(batch_j, req0, req1_j)
!
      call mem%alloc(X_J, wf%eri%n_J)
      call zero_array(X_J, wf%eri%n_J)
!
      call mem%alloc(L_J_jb, wf%eri%n_J, batch_j%max_length, wf%n_v)
      call mem%alloc(c_jb, batch_j%max_length*wf%n_v)
!
      do current_j_batch = 1, batch_j%num_batches
!
         call batch_j%determine_limits(current_j_batch)
         call wf%eri%get_cholesky_t1(L_J_jb, batch_j%first, batch_j%last,&
                                                   wf%n_o + 1, wf%n_mo)
!
!$omp parallel do private (b, j)
         do b = 1, wf%n_v
            do j = 1, batch_j%length
!
               jb = batch_j%length*(b - 1) + j
               c_jb(jb) = c_bj(b, j + batch_j%first - 1)
!
            enddo
         enddo
!$omp end parallel do
!
         call dgemv('N',                     &
                     wf%eri%n_J,             &
                     wf%n_v*batch_j%length,  &
                     one,                    &
                     L_J_jb,                 &
                     wf%eri%n_J,             &
                     c_jb,                   &
                     1,                      &
                     one,                    &
                     X_J,                    &
                     1)
!
      enddo !batch_j
!
      call mem%dealloc(L_J_jb, wf%eri%n_J, batch_j%max_length, wf%n_v)
      call mem%dealloc(c_jb, wf%n_v*batch_j%max_length)
!
      req0 = 0
!
      req1_i = wf%n_v*wf%eri%n_J

      batch_i = batching_index(wf%n_o)
!
      call mem%batch_setup(batch_i, req0, req1_i)
!
      call mem%alloc(L_J_ai, wf%eri%n_J, wf%n_v, batch_i%max_length)
!
      do current_i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(current_i_batch)
!
         call wf%eri%get_cholesky_t1(L_J_ai, wf%n_o + 1, wf%n_mo ,&
                                                   batch_i%first, batch_i%last)
!
         call dgemm('T', 'N',                   &
                     wf%n_v*batch_i%length,     &
                     1,                         &
                     wf%eri%n_J,                &
                     two,                       &
                     L_J_ai,                    &
                     wf%eri%n_J,                &
                     X_J,                       &
                     wf%eri%n_J,                &
                     one,                       &
                     rho_ai(1, batch_i%first),  &
                     wf%n_v*wf%n_o)
!
      enddo !batch_i
!
      call mem%dealloc(L_J_ai, wf%eri%n_J, wf%n_v, batch_i%max_length)
      call mem%dealloc(X_J, wf%eri%n_J)
!
!     :: Term 2 rho_ai = - g_abji * c_bj::
!
      req0 = 0
!
      req1_i = wf%eri%n_J*wf%n_o
      req1_a = max(wf%eri%n_J*wf%n_v +  wf%eri%n_J*wf%n_o, 2*wf%eri%n_J*wf%n_o)
!
      req2 = 1
!
      batch_i = batching_index(wf%n_o)
      batch_a = batching_index(wf%n_v)
!
      call mem%batch_setup(batch_i, batch_a, req0, req1_i, req1_a, req2)
!
      do current_i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(current_i_batch)
!
         do current_a_batch = 1, batch_a%num_batches
!
            call batch_a%determine_limits(current_a_batch)
!
            call mem%alloc(L_J_ab, wf%eri%n_J, batch_a%length, wf%n_v)
!
            call wf%eri%get_cholesky_t1(L_J_ab, wf%n_o + batch_a%first, &
                                                      wf%n_o + batch_a%last,  &
                                                      wf%n_o + 1, &
                                                      wf%n_mo)
!
            call mem%alloc(X_Jaj, wf%eri%n_J, batch_a%length, wf%n_o)
!
            call dgemm('N', 'N',                   &
                        wf%eri%n_J*batch_a%length, &
                        wf%n_o,                    &
                        wf%n_v,                    &
                        one,                       &
                        L_J_ab,                    & ! L_Ja_b
                        wf%eri%n_J*batch_a%length, &
                        c_bj,                      &
                        wf%n_v,                    &
                        zero,                      &
                        X_Jaj,                     &
                        wf%eri%n_J*batch_a%length)
!
            call mem%dealloc(L_J_ab, wf%eri%n_J, batch_a%length, wf%n_v)
!
            call mem%alloc(X_aJj, batch_a%length, wf%eri%n_J, wf%n_o)
            call sort_123_to_213(X_Jaj, X_aJj, wf%eri%n_J, batch_a%length, wf%n_o)
            call mem%dealloc(X_Jaj, wf%eri%n_J, batch_a%length, wf%n_o)
!
            call mem%alloc(L_J_ji, wf%eri%n_J, wf%n_o, batch_i%length)
!
            call wf%eri%get_cholesky_t1(L_J_ji, 1, wf%n_o, &
                                                      batch_i%first, batch_i%last)
!
            call mem%alloc(rho_ai_batch, batch_a%length, batch_i%length)
!
            call dgemm('N', 'N',             &
                        batch_a%length,      &
                        batch_i%length,      &
                        wf%eri%n_J*wf%n_o,   &
                        -one,                &
                        X_aJj,               & ! X_a_Jj
                        batch_a%length,      &
                        L_J_ji,              & ! L_Jj_i
                        wf%eri%n_J*wf%n_o,   &
                        zero,                &
                        rho_ai_batch,        &
                        batch_a%length)
!
            call mem%dealloc(L_J_ji, wf%eri%n_J, wf%n_o, batch_i%length)
            call mem%dealloc(X_aJj, batch_a%length, wf%eri%n_J, wf%n_o)
!
!$omp parallel do private (a, i)
            do i = 1, batch_i%length
               do a = 1, batch_a%length
!
                  rho_ai(a + batch_a%first - 1, i + batch_i%first - 1) = &
                     rho_ai(a + batch_a%first - 1, i + batch_i%first - 1) + rho_ai_batch(a, i)
!
               enddo
            enddo
!$omp end parallel do
!
            call mem%dealloc(rho_ai_batch, batch_a%length, batch_i%length)
!
         enddo ! batch_b
      enddo ! batch_i
!
      call timer%turn_off()
!
   end subroutine jacobian_ccs_b1_ccs
!
!
end submodule jacobian_ccs
