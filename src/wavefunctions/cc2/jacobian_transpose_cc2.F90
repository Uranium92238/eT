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
submodule (cc2_class) jacobian_transpose_cc2
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
!!    A_μ,ν = < μ |exp(-T) [H, τ_ν] exp(T) | R >.
!!
!!
!
   implicit none
!
!
contains
!
!
   module subroutine prepare_for_jacobian_transpose_cc2(wf)
!!
!!    Jacobian transpose submodule (CC2)
!!    Written by Sarai D. Folkestad and Alexander C. Paul, Feb 2019
!!
!!    Modified by Tor S. Haugland, Oct 2019
!!
!!    Removed construction of wf%u as it is constructed in the ground state.
!!
      implicit none
!
      class(cc2), intent(inout) :: wf
!
      type(timings), allocatable :: timer
!
      timer = timings('Prepare for Jacobian transpose CC2', pl='normal')
      call timer%turn_on()
!
      call wf%save_jacobian_transpose_a1_intermediates(wf%u_aibj)
!
      call timer%turn_off()
!
   end subroutine prepare_for_jacobian_transpose_cc2
!
!
   module subroutine jacobian_transpose_transformation_cc2(wf, b, sigma)
!!
!!    Jacobian transpose transformation
!!    Written by Sarai D. Folkestad and Alexander C. Paul, Feb 2019
!!
!!    Calculates the transpose Jacobian transformation, i.e., the transformation
!!    by the transpose of the Jacobian matrix
!!
!!       A_mu,nu = < mu | exp(-T) [H, tau_nu] exp(T) | R >.
!!
!!    The transformation is performed as sigma^T = c^T A, where c is the vector
!!    sent to the routine.
!!
      use array_utilities, only: zero_array
      use reordering, only: squareup, symmetric_sum, packin
!
      implicit none
!
      class(cc2), intent(inout) :: wf
!
      real(dp), dimension(wf%n_t1 + wf%n_t2), intent(in)  :: b
      real(dp), dimension(wf%n_t1 + wf%n_t2), intent(out) :: sigma
!
      real(dp), dimension(:,:,:,:), allocatable :: b_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: sigma_aibj
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian transpose CC2', pl='normal')
      call timer%turn_on()
!
      call zero_array(sigma, wf%n_t1 + wf%n_t2)
!
!     :: CCS contributions to the singles b vector ::
!
      call wf%ccs%jacobian_transpose_transformation(b(1 : wf%n_t1), sigma(1 : wf%n_t1))
!
!     :: CC2 contributions to the transformed singles vector ::
!
      call wf%jacobian_transpose_doubles_a1(sigma(1 : wf%n_t1), b(1 : wf%n_t1), wf%u_aibj)
!
      call mem%alloc(b_aibj, (wf%n_v), (wf%n_o), (wf%n_v), (wf%n_o))
!
      call squareup(b(wf%n_t1+1:), b_aibj, wf%n_t1)
!
      call wf%jacobian_transpose_doubles_b1(sigma(1 : wf%n_t1), b_aibj)
!
!     :: CC2 contributions to the transformed doubles vector ::
!
      call mem%alloc(sigma_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call zero_array(sigma_aibj, wf%n_t1**2)
!
!     Contributions from singles vector b
!
      call wf%jacobian_transpose_doubles_a2(sigma_aibj, b(1 : wf%n_t1))
      call symmetric_sum(sigma_aibj, wf%n_t1)
!
!     Contributions from doubles vector b
!
      call wf%jacobian_transpose_cc2_b2(sigma_aibj, b_aibj)
!
      call mem%dealloc(b_aibj, (wf%n_v), (wf%n_o), (wf%n_v), (wf%n_o))
!
!     Overwrite the incoming doubles b vector & pack in
!
      call packin(sigma(wf%n_t1+1:), sigma_aibj, wf%n_t1)
!
      call mem%dealloc(sigma_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_transpose_transformation_cc2
!
!
   module subroutine jacobian_transpose_cc2_b2_cc2(wf, sigma_aibj, c_aibj)
!!
!!    Jacobian transpose CC2 B2
!!    Written by Sarai D. Folkestad and Alexander C. Paul, Feb 2019
!!
!!    Calculates the A2 term,
!!
!!    sigma_aibj =+ ε_aibj c_aibj
!!
      implicit none
!
      class(cc2) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)      :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout)   :: sigma_aibj
!
      integer :: a, i, b, j
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian transpose CC2 B2', pl='verbose')
      call timer%turn_on()
!
!$omp parallel do private(a, i, b, j)
      do j = 1, wf%n_o
         do b = 1, wf%n_v
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  sigma_aibj(a, i, b, j) = sigma_aibj(a, i, b, j) + c_aibj(a,i,b,j) &
                                          * (wf%orbital_energies(a + wf%n_o) &
                                          + wf%orbital_energies(b + wf%n_o) &
                                          - wf%orbital_energies(i) &
                                          - wf%orbital_energies(j))
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call timer%turn_off()
!
   end subroutine jacobian_transpose_cc2_b2_cc2
!
!
end submodule jacobian_transpose_cc2
