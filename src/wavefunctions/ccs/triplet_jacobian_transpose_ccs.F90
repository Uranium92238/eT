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
submodule (ccs_class) triplet_jacobian_transpose_ccs
!
!!
!!    Triplet Jacobian transpose submodule
!!
!!    Routines for the linear transform of trial
!!    vectors by the Jacobian matrix
!!
!!    sigma_i = A^T * c_i,
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
   module subroutine prepare_for_triplet_jacobian_transpose_ccs(wf)
!!
!!    Prepare for triplet Jacobian Transpose
!!    Written by Sarai D. Folkestad, Feb 2022
!!
!
      use warning_suppressor, only: do_nothing
!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      call do_nothing(wf)
!
   end subroutine prepare_for_triplet_jacobian_transpose_ccs
!
!
   module subroutine triplet_jacobian_transpose_transformation_ccs(wf, b, sigma)
!!
!!    Triplet Jacobian transformation
!!    Written by Sarai D. Folkestad, Feb 2022
!!
!!    Directs the transformation by the CCS Jacobi matrix,
!!
!!       A_μ,ν = < μ |exp(-T) [H, τ_ν] exp(T) | R >,
!!
!!    Where τ_ν is a triplet excitation operator,
!!    and < μ | = 1/2 < R | τ_μ^dagger
!!
!!    In particular,
!!
!!       sigma_mu = (A c)_mu = sum_ck A_mu,ck c_ck.
!!
      use array_initialization, only: zero_array
!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(wf%n_t1), intent(in) :: b
      real(dp), dimension(wf%n_t1), intent(out) :: sigma
!
      type(timings), allocatable :: timer
!
      timer = timings('Triplet Jacobian CCS transformation', pl='normal')
      call timer%turn_on()
!
      call zero_array(sigma, wf%n_t1)

      call wf%jacobian_transpose_ccs_a1(sigma, b)
      call wf%triplet_jacobian_transpose_s_s_a(sigma, b)
!
      call timer%turn_off()
!
   end subroutine triplet_jacobian_transpose_transformation_ccs
!
!
   module subroutine triplet_jacobian_transpose_s_s_a_ccs(wf, sigma_ai, c_bj)
!!
!!    Triplet Jacobian transpose singles-singles A
!!    Written by Sarai D. Folkestad, Feb 2022
!!
!!    Calculates the term,
!!
!!       sum_bj -g_ijba c_bj
!!
!!    and adds it to the rho vector.
!!
!
      use reordering, only: sort_123_to_132
      use array_initialization, only: zero_array
!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_bj
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
!
      real(dp), dimension(:,:,:), pointer :: L_Jba, L_Jij
!
      real(dp), dimension(:,:,:), allocatable :: L_Jbi, L_Jib
!
      type(timings), allocatable :: timer
!
      integer :: req0, req1_i, req1_a, req2
!
      integer :: current_a_batch, current_i_batch
!
      type(batching_index) :: batch_a, batch_i
!
      timer = timings('Triplet Jacobian CCS A', pl='verbose')
      call timer%turn_on()
!
      req0 = 0
!
      req1_i = wf%L_t1%n_J*wf%n_v &
               + wf%L_t1%load_memory_estimate(1, 1, 1, wf%n_o)
      req1_a = max(wf%L_t1%n_J*wf%n_v, wf%L_t1%n_J*wf%n_o) &
               + wf%L_t1%load_memory_estimate(wf%n_o + 1, wf%n_mo, wf%n_o + 1, wf%n_o + 1)
!
      req2 = 0
!
      batch_i = batching_index(wf%n_o)
      batch_a = batching_index(wf%n_v)
!
      call mem%batch_setup(batch_i, batch_a, req0, req1_i, req1_a, req2, &
                           tag='triplet_jacobian_transpose_ccs')
!
      do current_i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(current_i_batch)
         call wf%L_t1%load_block(L_Jij, batch_i%first, batch_i%get_last(), 1, wf%n_o)
!
         call mem%alloc(L_Jib, wf%L_t1%n_J, batch_i%length, wf%n_v)
         call dgemm('N', 'T',                   &
                    wf%L_t1%n_J*batch_i%length, &
                    wf%n_v,                     &
                    wf%n_o,                     &
                    one,                        &
                    L_Jij,                      &
                    wf%L_t1%n_J*batch_i%length, &
                    c_bj,                       &
                    wf%n_v,                     &
                    zero,                       &
                    L_Jib,                      &
                    wf%L_t1%n_J*batch_i%length)
!
         call wf%L_t1%offload_block(batch_i%first, batch_i%get_last(), 1, wf%n_o)
!
         call mem%alloc(L_Jbi, wf%L_t1%n_J, wf%n_v, batch_i%length)
         call sort_123_to_132(L_Jib, L_Jbi, wf%L_t1%n_J, batch_i%length, wf%n_v)
         call mem%dealloc(L_Jib, wf%L_t1%n_J, batch_i%length, wf%n_v)
!
         do current_a_batch = 1, batch_a%num_batches
!
            call batch_a%determine_limits(current_a_batch)
!
            call wf%L_t1%load_block(L_Jba, 1 + wf%n_o, wf%n_mo, &
                                    wf%n_o + batch_a%first,     &
                                    wf%n_o + batch_a%get_last())
!
            call dgemm('T', 'N',                               &
                        batch_a%length,                        &
                        batch_i%length,                        &
                        wf%n_v*wf%L_t1%n_J,                    &
                        -one,                                  &
                        L_Jba,                                 &
                        wf%n_v*wf%L_t1%n_J,                    &
                        L_Jbi,                                 &
                        wf%n_v*wf%L_t1%n_J,                    &
                        one,                                   &
                        sigma_ai(batch_a%first, batch_i%first),&
                        wf%n_v)
!
            call wf%L_t1%offload_block(1 + wf%n_o, wf%n_mo, &
                                    wf%n_o + batch_a%first,     &
                                    wf%n_o + batch_a%get_last())
!
         enddo

         call mem%dealloc(L_Jbi, wf%L_t1%n_J, wf%n_v, batch_i%length)
      enddo
      call mem%batch_finalize()
!
      call timer%turn_off()
!
   end subroutine triplet_jacobian_transpose_s_s_a_ccs
!
!
end submodule triplet_jacobian_transpose_ccs
