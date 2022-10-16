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
submodule (ccs_class) triplet_jacobian_ccs
!
!!
!!    Triplet Jacobian submodule
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
   module subroutine prepare_for_triplet_jacobian_ccs(wf)
!!
!!    Prepare for triplet Jacobian
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2019
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
   end subroutine prepare_for_triplet_jacobian_ccs
!
!
   module subroutine triplet_jacobian_transformation_ccs(wf, c, rho)
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
!!       rho_mu = (A c)_mu = sum_ck A_mu,ck c_ck.
!!
      use array_initialization, only: zero_array
!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(wf%n_t1), intent(in) :: c
      real(dp), dimension(wf%n_t1), intent(out) :: rho
!
      type(timings), allocatable :: timer
!
      timer = timings('Triplet Jacobian CCS transformation', pl='normal')
      call timer%turn_on()
!
      call zero_array(rho, wf%n_t1)

      call wf%jacobian_ccs_a1(rho, c)
      call wf%triplet_jacobian_s_s_a(rho, c)
!
      call timer%turn_off()
!
   end subroutine triplet_jacobian_transformation_ccs
!
!
!
   module subroutine triplet_jacobian_s_s_a_ccs(wf, rho_ai, c_bj)
!!
!!    Triplet Jacobian singles-singles A
!!    Written by Sarai D. Folkestad, Feb 2022
!!
!!    Calculates the term,
!!
!!       sum_bj -g_abji c_bj
!!
!!    and adds it to the rho vector.
!!
!
      use reordering, only: sort_123_to_132
      use array_initialization, only: zero_array
!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_bj
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
!
      real(dp), dimension(:,:,:), pointer :: L_Jab, L_Jji
!
      real(dp), dimension(:,:,:), allocatable :: L_Jaj, L_Jja
!
      integer :: req0, req1_j, req1_b, req2
!
      integer :: current_b_batch, current_j_batch
!
      type(batching_index) :: batch_b, batch_j
!
      type(timings), allocatable :: timer
!
      timer = timings('Triplet Jacobian CCS A', pl='verbose')
      call timer%turn_on()
!
      req0 = 0
!
      req1_j = 2*wf%L_t1%n_J*wf%n_v + wf%L_t1%load_memory_estimate(1,1,1,wf%n_o)
      req1_b = wf%L_t1%load_memory_estimate(1,wf%n_v,1,1)
!
      req2 = 0
!
      batch_j = batching_index(wf%n_o)
      batch_b = batching_index(wf%n_v)
!
      call mem%batch_setup(batch_j, batch_b, req0, req1_j, req1_b, req2, &
                           tag='triplet_jacobian_ccs')
!
      do current_j_batch = 1, batch_j%num_batches
!
         call batch_j%determine_limits(current_j_batch)
!
         call mem%alloc(L_Jaj, wf%L_t1%n_J, wf%n_v, batch_j%length)
         call zero_array(L_Jaj, wf%L_t1%n_J*wf%n_v*batch_j%length)
!
         do current_b_batch = 1, batch_b%num_batches
!
            call batch_b%determine_limits(current_b_batch)
!
            call wf%L_t1%load_block(L_Jab, wf%n_o + 1, wf%n_mo, &
                                           wf%n_o + batch_b%first, &
                                           wf%n_o + batch_b%get_last())
!
            call dgemm('N', 'N',              &
                        wf%L_t1%n_J*wf%n_v,   &
                        batch_j%length,       &
                        batch_b%length,       &
                        one,                  &
                        L_Jab,                &
                        wf%L_t1%n_J*wf%n_v,   &
                        c_bj(batch_b%first,batch_j%first),&
                        wf%n_v,               &
                        one,                  &
                        L_Jaj,                &
                        wf%L_t1%n_J*wf%n_v)
!
            call wf%L_t1%offload_block(wf%n_o + 1, wf%n_mo, &
                                       wf%n_o + batch_b%first, &
                                       wf%n_o + batch_b%get_last())
!
         enddo
!
         call mem%alloc(L_Jja, wf%L_t1%n_J, batch_j%length, wf%n_v)
         call sort_123_to_132(L_Jaj, L_Jja, wf%L_t1%n_J, wf%n_v, batch_j%length)
         call mem%dealloc(L_Jaj, wf%L_t1%n_J, wf%n_v, batch_j%length)
!
         call wf%L_t1%load_block(L_Jji, batch_j%first, batch_j%get_last(), 1, wf%n_o)
!
         call dgemm('T', 'N',                     &
                     wf%n_v,                      &
                     wf%n_o,                      &
                     wf%L_t1%n_J*batch_j%length,  &
                     -one,                        &
                     L_Jja,                       &
                     wf%L_t1%n_J*batch_j%length,  &
                     L_Jji,                       &
                     wf%L_t1%n_J*batch_j%length,  &
                     one,                         &
                     rho_ai,                      &
                     wf%n_v)
!
         call wf%L_t1%offload_block(batch_j%first, batch_j%get_last(), 1, wf%n_o)
         call mem%dealloc(L_Jja, wf%L_t1%n_J, batch_j%length, wf%n_v)
!
      enddo
!
      call mem%batch_finalize()
!
      call timer%turn_off()
!
   end subroutine triplet_jacobian_s_s_a_ccs
!
!
end submodule triplet_jacobian_ccs
