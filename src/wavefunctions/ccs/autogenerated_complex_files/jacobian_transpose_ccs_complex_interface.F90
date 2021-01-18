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
   module subroutine prepare_for_jacobian_transpose_ccs_complex(wf)
!!
!!    Prepare for jacobian transpose
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2019
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
   end subroutine prepare_for_jacobian_transpose_ccs_complex
!
!
   module subroutine jacobian_transpose_transformation_ccs_complex(wf, b, sigma)
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
      implicit none
!
      class(ccs), intent(inout) :: wf
      complex(dp), dimension(wf%n_t1), intent(in)  :: b
      complex(dp), dimension(wf%n_t1), intent(out) :: sigma
!
   end subroutine jacobian_transpose_transformation_ccs_complex
!
!
   module subroutine jacobian_transpose_ccs_a1_ccs_complex(wf, sigma_ai, b_ai)
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
      complex(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      complex(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: b_ai
!
   end subroutine jacobian_transpose_ccs_a1_ccs_complex
!
!
   module subroutine jacobian_transpose_ccs_b1_ccs_complex(wf, sigma_ai, b_ai)
!!
!!    Jacobian transpose B1 (CCS)
!!    Written by Sarai D. Folkestad, Jun 2019
!!
!!    Calculates the (CCS) B1 term of the Jacobian transpose 
!!    transfromation. 
!!
!!       B1 = sum_bj L_bjia b_bj
!!          = sum_bj (2 g_bjia b_bj - g_baij b_bj)
!!    
      implicit none
!
      class(ccs), intent(inout) :: wf
      complex(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      complex(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: b_ai
!
   end subroutine jacobian_transpose_ccs_b1_ccs_complex
