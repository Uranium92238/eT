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
submodule (cc2_class) jacobian_transpose_cc2
!
!!
!!    Jacobian transpose submodule (CC2)
!!    Written by Sarai D. Folkestad and Alexander C. Paul, Feb 2019
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
!
   module subroutine prepare_for_jacobian_transpose_cc2(wf)
!!
!!    Jacobian transpose submodule (CC2)
!!    Written by Sarai D. Folkestad and Alexander C. Paul, Feb 2019
!!
      implicit none
!
      class(cc2), intent(inout) :: wf
!
      call wf%initialize_u()
      call wf%construct_u()
!
   end subroutine prepare_for_jacobian_transpose_cc2
!
!
   module subroutine jacobian_transpose_transformation_cc2(wf, b)
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
!!    sent to the routine. On exit, the vector c is equal to sigma (the transformed
!!    vector).
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: b
!
      real(dp), dimension(:,:), allocatable :: b_ai
      real(dp), dimension(:,:,:,:), allocatable :: b_aibj
!
      real(dp), dimension(:,:), allocatable :: sigma_ai
      real(dp), dimension(:,:,:,:), allocatable :: sigma_aibj
!
      type(timings) :: timer
!
      timer = timings('Jacobian transpose transformation CC2')
      call timer%turn_on()
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
!     :: CCS contributions to the singles b vector ::
!
      call wf%jacobian_transpose_ccs_a1(sigma_ai, b_ai)
      call wf%jacobian_transpose_ccs_b1(sigma_ai, b_ai)
!
!     :: CC2 contributions to the transformed singles vector ::
!
      call wf%jacobian_transpose_doubles_a1(sigma_ai, b_ai, wf%u)
!
!     Allocate the incoming unpacked doubles vector
!
      call mem%alloc(b_aibj, (wf%n_v), (wf%n_o), (wf%n_v), (wf%n_o))
!
      call squareup(b(wf%n_t1 + 1 : wf%n_es_amplitudes), b_aibj, wf%n_t1)
!
      call wf%jacobian_transpose_doubles_b1(sigma_ai, b_aibj)
!
!     Done with singles vector c; overwrite it with
!     transformed vector for exit
!
      call dcopy(wf%n_t1, sigma_ai, 1, b, 1)
!
      call mem%dealloc(sigma_ai, wf%n_v, wf%n_o)
!
!     :: CC2 contributions to the transformed doubles vector ::
!
!     Allocate unpacked transformed vector
!
      call mem%alloc(sigma_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call zero_array(sigma_aibj, wf%n_t1**2)
!
!     Contributions from singles vector b
!
      call wf%jacobian_transpose_doubles_a2(sigma_aibj, b_ai)
      call symmetric_sum(sigma_aibj, wf%n_t1)
!
      call mem%dealloc(b_ai, wf%n_v, wf%n_o)
!
!     Contributions from doubles vector b
!
      call wf%jacobian_transpose_cc2_b2(sigma_aibj, b_aibj)
!
      call mem%dealloc(b_aibj, (wf%n_v), (wf%n_o), (wf%n_v), (wf%n_o))
!
!     Overwrite the incoming doubles b vector & pack in
!
      call packin(b(wf%n_t1 + 1 : wf%n_es_amplitudes), sigma_aibj, wf%n_t1)
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
      type(timings) :: timer
!
      timer = timings('jacobian transpose b2 cc2')
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
