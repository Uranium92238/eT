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
   module subroutine prepare_for_jacobian_transpose_cc2(wf)
!!
!!    Jacobian transpose transform trial vector
!!    Written by Sarai D. Folkestad and Alexander Paul, Feb 2019
!!
      implicit none
!
      class(cc2), intent(inout) :: wf
!
   end subroutine prepare_for_jacobian_transpose_cc2
!
!
   module subroutine jacobian_transpose_transform_trial_vector_cc2(wf, c_i)
!!
!!    Jacobian transpose transform trial vector
!!    Written by Sarai D. Folkestad and Alexander Paul, Feb 2019
!!
      class(cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes) :: c_i
!
   end subroutine jacobian_transpose_transform_trial_vector_cc2
!
!
   module subroutine jacobian_transpose_cc2_transformation_cc2(wf, c)
!!
!!    Jacobian transpose transformation (CC2)
!!    Written by Sarai D. Folkestad and Alexander Paul, Feb 2019
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
      class(cc2) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes) :: c
!
!
   end subroutine jacobian_transpose_cc2_transformation_cc2
!
!
   module subroutine jacobian_transpose_cc2_a1_cc2(wf, sigma_ai, c_bj)
!!
!!    Jacobian transpose CC2 A1
!!    Written by Sarai D. Folkestad and Alexander Paul, Feb 2019
!!
      implicit none
!
      class(cc2) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_bj
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
!
   end subroutine jacobian_transpose_cc2_a1_cc2
!
!
   module subroutine jacobian_transpose_cc2_b1_cc2(wf, sigma_ai, c_bjck)
!!
!!    Jacobian transpose CC2 B1
!!    Written by Sarai D. Folkestad and Alexander Paul, Feb 2019
!!
      implicit none
!
      class(cc2) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)                :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)   :: c_bjck
!
   end subroutine jacobian_transpose_cc2_b1_cc2
!
!
   module subroutine jacobian_transpose_cc2_a2_cc2(wf, sigma_aibj, c_ai)
!!
!!    Jacobian transpose CC2 A2
!!    Written by Sarai D. Folkestad and Alexander Paul, Feb 2019
!!
      implicit none
!
      class(cc2) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                      :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout)   :: sigma_aibj
!
   end subroutine jacobian_transpose_cc2_a2_cc2
!
!
   module subroutine jacobian_transpose_cc2_b2_cc2(wf, sigma_aibj, c_aibj)
!!
!!    Jacobian transpose CC2 B2
!!    Written by Sarai D. Folkestad and Alexander Paul, Feb 2019
!!
      implicit none
!
      class(cc2) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout)   :: sigma_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)      :: c_aibj
!
   end subroutine jacobian_transpose_cc2_b2_cc2
!
!
