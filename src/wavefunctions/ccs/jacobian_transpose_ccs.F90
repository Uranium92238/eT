!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2021 the authors of eT
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
submodule (ccs_class) jacobian_transpose_ccs
!
!!
!!    Jacobian transpose submodule
!!
!!    Routines for the linear transform of trial
!!    vectors by the transpose of the Jacobian matrix
!!
!!    σ_i = A^T * b_i,
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
   module subroutine prepare_for_jacobian_transpose_ccs(wf)
!!
!!    Prepare for jacobian transpose
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
   end subroutine prepare_for_jacobian_transpose_ccs
!
!
   module subroutine jacobian_transpose_transformation_ccs(wf, b, sigma)
!!
!!    Jacobian transpose transformation
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
!!    Calculates the transpose Jacobian transformation, i.e., the transformation
!!    by the transpose of the Jacobian matrix
!!
!!       A_mu,nu = < mu | exp(-T) [H, τ_nu] exp(T) | R >.
!!
!!    In particular,
!!
!!       sigma_mu = (b^T A)_mu = sum_ck b_ck A_ck,mu.
!!
      use array_utilities, only: zero_array
!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(wf%n_t1), intent(in)  :: b
      real(dp), dimension(wf%n_t1), intent(out) :: sigma
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian transpose CCS', pl='normal')
      call timer%turn_on()
!
      call zero_array(sigma, wf%n_t1)
!
      call wf%jacobian_transpose_ccs_a1(sigma, b)
      call wf%jacobian_transpose_ccs_b1(sigma, b)
!
      call timer%turn_off()
!
   end subroutine jacobian_transpose_transformation_ccs
!
!
   module subroutine jacobian_transpose_ccs_a1_ccs(wf, sigma_ai, b_ai)
!!
!!    Jacobian transpose A1 (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
!!    Calculates the A1 term,
!!
!!       sum_c b_ci F_ca - sum_k b_ak F_ik,
!!
!!    and adds it to the sigma-vector (b^T -> sigma^T = b^T A).
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: b_ai
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian transpose CCS A1', pl='verbose')
      call timer%turn_on()
!
!     Add sum_c F_ca b_ci = sum_c F_ac^T b_ci
!
      call dgemm('T','N',     &
                  wf%n_v,     &
                  wf%n_o,     &
                  wf%n_v,     &
                  one,        &
                  wf%fock_ab, &
                  wf%n_v,     &
                  b_ai,       &
                  wf%n_v,     &
                  one,        &
                  sigma_ai,   &
                  wf%n_v)
!
!     Add - sum_k b_ak F_ik = - sum_k b_ak F_ki^T
!
      call dgemm('N','T',     &
                  wf%n_v,     &
                  wf%n_o,     &
                  wf%n_o,     &
                  -one,       &
                  b_ai,       &
                  wf%n_v,     &
                  wf%fock_ij, &
                  wf%n_o,     &
                  one,        &
                  sigma_ai,   &
                  wf%n_v)
!
      call timer%turn_off()
!
   end subroutine jacobian_transpose_ccs_a1_ccs
!
!
   module subroutine jacobian_transpose_ccs_b1_ccs(wf, sigma_ai, b_bj)
!!
!!    Jacobian transpose B1 (CCS)
!!    Written by Sarai D. Folkestad, Jun 2019
!!
!!    Calculates the (CCS) B1 term of the Jacobian transpose
!!    transfromation.
!!
!!       B1 = sum_bj L_bjia b_bj
!!          = sum_bj (2 g_bjia b_bj - g_baij b_bj)
!!          = sum_bjJ (2 L^J_bj L^J_ia b_bj - L^J_baL^J_ij b_bj)
!!
      use batching_index_class, only: batching_index
      use array_utilities, only: zero_array
      use reordering, only: sort_123_to_132
!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: b_bj
!
      integer :: req, req_i, req_a, req_ai
      integer :: current_a_batch, current_j_batch, current_i_batch
!
      type(batching_index) batch_a, batch_j, batch_i
!
      real(dp), dimension(:,:,:), allocatable :: L_Jbj, L_Jij
      real(dp), dimension(:,:,:), allocatable :: L_Jba, Y_Jib, Y_Jbi, L_Jia
      real(dp), dimension(:,:), allocatable :: sigma_ia
      real(dp), dimension(:), allocatable :: X_J
!
      type(timings), allocatable :: timer
!
      integer :: i, a
!
      timer = timings('Jacobian transpose CCS B1', pl='verbose')
      call timer%turn_on()
!
      call mem%alloc(X_J, wf%eri%n_J)
      call zero_array(X_J, wf%eri%n_J)
!
      batch_j = batching_index(wf%n_o)
!
      req = wf%eri%n_J*wf%n_v
      call mem%batch_setup(batch_j, 0, req, tag='jacobian_transpose_ccs_b1_ccs 1')
!
      call mem%alloc(L_Jbj, wf%eri%n_J, wf%n_v, batch_j%max_length)
!
      do current_j_batch = 1, batch_j%num_batches
!
         call batch_j%determine_limits(current_j_batch)
!
         call wf%eri%get_cholesky_t1(L_Jbj, wf%n_o + 1, wf%n_mo, batch_j%first, batch_j%get_last())
!
         call dgemm('N', 'N',                &
                     wf%eri%n_J,             &
                     1,                      &
                     batch_j%length*wf%n_v,  &
                     one,                    &
                     L_Jbj,                  &
                     wf%eri%n_J,             &
                     b_bj(1,batch_j%first),  &
                     batch_j%length*wf%n_v,  &
                     one,                    &
                     X_J,                    &
                     wf%eri%n_J)
!
      enddo
!
      call mem%dealloc(L_Jbj, wf%eri%n_J, wf%n_v, batch_j%max_length)
!
      call mem%batch_finalize()
!
      batch_i = batching_index(wf%n_o)
!
      req = 2*wf%eri%n_J*wf%n_v + wf%n_v
      call mem%batch_setup(batch_i, 0, req, tag='jacobian_transpose_ccs_b1_ccs 2')
!
      do current_i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(current_i_batch)

         call mem%alloc(L_Jia, wf%eri%n_J, batch_i%length, wf%n_v)
         call wf%eri%get_cholesky_t1(L_Jia, batch_i%first, batch_i%get_last(), wf%n_o + 1, wf%n_mo)
!
         call mem%alloc(sigma_ia, batch_i%length, wf%n_v)
         call dgemm('T', 'N',                   &
                     wf%n_v*batch_i%length,     &
                     1,                         &
                     wf%eri%n_J,                &
                     two,                       &
                     L_Jia,                     &
                     wf%eri%n_J,                &
                     X_J,                       &
                     wf%eri%n_J,                &
                     zero,                      &
                     sigma_ia,                  &
                     wf%n_v*batch_i%length)
!
         call mem%dealloc(L_Jia, wf%eri%n_J, batch_i%length, wf%n_v)
!
!$omp parallel do private (i, a)
         do i = batch_i%first, batch_i%get_last()
            do a = 1, wf%n_v
               sigma_ai(a, i) = sigma_ai(a, i) + sigma_ia(i - batch_i%first + 1, a)
            enddo
         enddo
!$omp end parallel do
!
         call mem%dealloc(sigma_ia, batch_i%length, wf%n_v)
!
      enddo
!
      call mem%dealloc(X_J, wf%eri%n_J)
!
      call mem%batch_finalize()
!
      batch_i = batching_index(wf%n_o)
      batch_a = batching_index(wf%n_v)
!
      req = 0
      req_a = 2*wf%eri%n_J*wf%n_v
      req_i = 2*wf%eri%n_J*wf%n_o
      req_ai = wf%eri%n_J
!
      call mem%batch_setup(batch_a, batch_i, req, req_a, req_i, req_ai, &
                           tag='jacobian_transpose_ccs_b1_ccs 3')
!
      call mem%alloc(L_Jba, wf%eri%n_J, wf%n_v, batch_a%max_length)
!
      do current_i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(current_i_batch)
!
         call mem%alloc(L_Jij, wf%eri%n_J, batch_i%length, wf%n_o)
         call mem%alloc(Y_Jbi, wf%eri%n_J, wf%n_v, batch_i%length)
         call mem%alloc(Y_Jib, wf%eri%n_J, batch_i%length, wf%n_v)
!
         call wf%eri%get_cholesky_t1(L_Jij, batch_i%first, batch_i%get_last(), 1, wf%n_o)
!
         call dgemm('N', 'T',                   &
                     wf%eri%n_J*batch_i%length, &
                     wf%n_v,                    &
                     wf%n_o,                    &
                     one,                       &
                     L_Jij,                     &
                     wf%eri%n_J*batch_i%length, &
                     b_bj,                      &
                     wf%n_v,                    &
                     zero,                      &
                     Y_Jib,                     &
                     wf%eri%n_J*batch_i%length)
!
         call mem%dealloc(L_Jij, wf%eri%n_J, batch_i%length, wf%n_o)
!
         call sort_123_to_132(Y_Jib, Y_Jbi, wf%eri%n_J, batch_i%length, wf%n_v)
         call mem%dealloc(Y_Jib, wf%eri%n_J, batch_i%length, wf%n_v)
!
         do current_a_batch = 1, batch_a%num_batches
!
            call batch_a%determine_limits(current_a_batch)
!
            call wf%eri%get_cholesky_t1(L_Jba, wf%n_o + 1, wf%n_mo, &
                                       wf%n_o + batch_a%first, wf%n_o + batch_a%get_last())
!
            call dgemm('T', 'N',                              &
                       batch_a%length,                        &
                       batch_i%length,                        &
                       wf%n_v*wf%eri%n_J,                     &
                       -one,                                  &
                       L_Jba,                                 &
                       wf%n_v*wf%eri%n_J,                     &
                       Y_Jbi,                                 &
                       wf%n_v*wf%eri%n_J,                     &
                       one,                                   &
                       sigma_ai(batch_a%first, batch_i%first),&
                       wf%n_v)
!
         enddo
!
         call mem%dealloc(Y_Jbi, wf%eri%n_J, wf%n_v, batch_i%length)
!
      enddo
!
      call mem%dealloc(L_Jba, wf%eri%n_J, wf%n_v, batch_a%max_length)

      call mem%batch_finalize()
!
      call timer%turn_off()
!
   end subroutine jacobian_transpose_ccs_b1_ccs
!
!
end submodule jacobian_transpose_ccs
