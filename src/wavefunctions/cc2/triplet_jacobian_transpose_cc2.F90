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
submodule (cc2_class) triplet_jacobian_transpose_cc2
!
!!
!!    Triplet Jacobian transpose submodule
!!
!!    Routines for the linear transform of trial
!!    vectors by the Jacobian matrix
!!
!!    sigma_i = A^T b_i
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
   module subroutine prepare_for_triplet_jacobian_transpose_cc2(wf)
!!
!!    Prepare for triplet Jacobian transpose
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2019
!!
      implicit none
!
      class(cc2), intent(inout) :: wf
!
      call wf%initialize_t2()
      call wf%construct_t2()
!
   end subroutine prepare_for_triplet_jacobian_transpose_cc2
!
!
   module subroutine triplet_jacobian_transpose_transformation_cc2(wf, b, sigma)
!!
!!    Triplet Jacobian transpose transformation
!!    Written by Sarai D. Folkestad, Feb 2022
!!
!!    Directs the transformation by the cc2 Jacobi matrix,
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
      use reordering, only: symmetric_sum, antisymmetric_sum, sort_1234_to_1432, squareup
!
      implicit none
!
      class(cc2), intent(inout) :: wf
!
      real(dp), dimension(wf%n_triplet_amplitudes), intent(in) :: b
      real(dp), dimension(wf%n_triplet_amplitudes), intent(out) :: sigma
!
      real(dp), dimension(:,:,:,:), allocatable :: t_aibj
      real(dp), dimension(:,:,:,:), allocatable :: b_aibj_p, b_aibj_m
      real(dp), dimension(:,:,:,:), allocatable :: b_aibj_pm
      real(dp), dimension(:,:,:,:), allocatable :: sigma_aibj_p, sigma_aibj_m, sigma_ajbi_p
!
      type(timings), allocatable :: timer
!
      timer = timings('Triplet Jacobian cc2 transformation', pl='normal')
      call timer%turn_on()
!
      call zero_array(sigma, wf%n_triplet_amplitudes)
!
      call wf%ccs%triplet_jacobian_transpose_transformation(b(1:wf%n_t1), sigma(1:wf%n_t1))
!
      call mem%alloc(t_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_aibj, wf%n_v*wf%n_o)
      call wf%triplet_jacobian_transpose_s_s_b(sigma(1:wf%n_t1), b(1:wf%n_t1), t_aibj)

      call mem%dealloc(t_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(b_aibj_p, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(b_aibj_m, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%packout_triplet_d(b, b_aibj_p, b_aibj_m)
!
      call mem%alloc(b_aibj_pm, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call dcopy(wf%n_o**2*wf%n_v**2, b_aibj_p, 1, b_aibj_pm, 1)
      call daxpy(wf%n_o**2*wf%n_v**2, one, b_aibj_m, 1, b_aibj_pm, 1)
!
      call wf%triplet_jacobian_transpose_d_s_a(sigma(1:wf%n_t1), b_aibj_pm)
!
      call mem%dealloc(b_aibj_pm, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(sigma_aibj_p, wf%n_v, wf%n_o, wf%n_v, wf%n_o, set_zero=.true.)
      call mem%alloc(sigma_aibj_m, wf%n_v, wf%n_o, wf%n_v, wf%n_o, set_zero=.true.)
!
      call wf%triplet_jacobian_transpose_s_d_a(sigma_aibj_p, b(1:wf%n_t1))
      call wf%triplet_jacobian_transpose_s_d_b(sigma_aibj_p, b(1:wf%n_t1))
      call wf%triplet_jacobian_transpose_s_d_c(sigma_aibj_p, b(1:wf%n_t1))
!
      call wf%triplet_jacobian_transpose_s_d_a(sigma_aibj_m, b(1:wf%n_t1))
      call wf%triplet_jacobian_transpose_s_d_b(sigma_aibj_m, b(1:wf%n_t1))
      call wf%triplet_jacobian_transpose_s_d_c(sigma_aibj_m, b(1:wf%n_t1))
!
      call symmetric_sum(sigma_aibj_p, wf%n_o*wf%n_v)
      call antisymmetric_sum(sigma_aibj_m, wf%n_o*wf%n_v)
!
      call mem%alloc(sigma_ajbi_p, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_1432(sigma_aibj_p, sigma_ajbi_p, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call daxpy(wf%n_v**2*wf%n_o**2, -one, sigma_ajbi_p, 1, sigma_aibj_p, 1)
      call mem%dealloc(sigma_ajbi_p, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%triplet_jacobian_d_d_a(sigma_aibj_p, b_aibj_p)
      call wf%triplet_jacobian_d_d_a(sigma_aibj_m, b_aibj_m)
!
      call mem%dealloc(b_aibj_p, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(b_aibj_m, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%packin_triplet_d(sigma, sigma_aibj_p, sigma_aibj_m)
!
      call mem%dealloc(sigma_aibj_p, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(sigma_aibj_m, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine triplet_jacobian_transpose_transformation_cc2
!
!
end submodule triplet_jacobian_transpose_cc2
