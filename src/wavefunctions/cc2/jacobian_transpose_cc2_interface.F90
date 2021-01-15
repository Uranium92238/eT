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
      implicit none
!
      class(cc2), intent(inout) :: wf
      real(dp), dimension(wf%n_t1 + wf%n_t2), intent(in)  :: b
      real(dp), dimension(wf%n_t1 + wf%n_t2), intent(out) :: sigma
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
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)      :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout)   :: sigma_aibj
!
   end subroutine jacobian_transpose_cc2_b2_cc2
