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
   module subroutine save_jacobian_transpose_a1_intermediates_doubles_complex(wf, u_bjck)
!!
!!    Save jacobian transpose A1 intermediates
!!    Written by by E. F. Kjønstad, S. D. Folkestad and Alexander C. Paul
!!
!!    Calculates the intermediates,
!!
!!       Y_ik = sum_cjb g_icjb * u_bjck
!!       Y_ca = sum_jbk u_bjck * g_jbka
!!
!!    and saves them into
!!
!!       jacobian_transpose_a1_intermdiate_oo
!!       jacobian_transpose_a1_intermdiate_vv
!!
!!    u_bjck = u^bc_jk =  2 t^bc_jk - t^bc_kj 
!!           = -(2 g_bjck - g_bkcj)/eps^bc_jk
!!    
!!    Adapted by Tor S. Haugland, Oct 2019
!!
!!    Isolated the intermediates from the
!!    jacobian_transpose_doubles_a1_doubles_complex and wrote them to file. 
!!
      class(doubles), intent(inout) :: wf
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: u_bjck
!
   end subroutine save_jacobian_transpose_a1_intermediates_doubles_complex
!
!
   module subroutine jacobian_transpose_doubles_a1_doubles_complex(wf, sigma_ai, c_bj, u)
!!
!!    Jacobian transpose doubles A1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad 
!!    and Alexander C. Paul, Feb 2019
!!
!!    Calculates the A1 term,
!!
!!    u_ckbj = u^bc_jk =  2 t^bc_jk - t^bc_kj = -(2 g_bjck - g_bkcj)/eps^bc_jk
!!
!!    sigma_ai += sum_bjck u^bc_jk (c_bj L_iakc - c_ak g_jbic - c_ci g_jbka)
!!             += sum_ck X_kc L_iakc - sum_k c_ak Y_ik - sum_c Y_ca c_ci
!!
!!    Modified by Tor S. Haugland, Oct 2019
!!
!!    Use saved intermediates to construct Y_ik and Y_ca.
!!    Create intermediate X_kc using transpose to save time re-ordering g_iakc
!!
      implicit none
!
      class(doubles) :: wf
      complex(dp), dimension(wf%n_v, wf%n_o), intent(in)     :: c_bj
      complex(dp), dimension(wf%n_v, wf%n_o), intent(inout)  :: sigma_ai
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: u
!
   end subroutine jacobian_transpose_doubles_a1_doubles_complex
!
!
  module subroutine jacobian_transpose_doubles_b1_doubles_complex(wf, sigma_ai, c_bjck)
!!
!!    Jacobian transpose doubles B1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad 
!!    and Alexander C. Paul, Feb 2019
!!
!!    Calculates the B1 term,
!!
!!    sigma_ai =+ sum_bjc c_bjci g_bjca - c_akbj g_bjik
!!
      implicit none
!
      class(doubles) :: wf
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_bjck
      complex(dp), dimension(wf%n_v, wf%n_o), intent(inout)  :: sigma_ai
!
   end subroutine jacobian_transpose_doubles_b1_doubles_complex
!
!
  module subroutine jacobian_transpose_doubles_a2_doubles_complex(wf, sigma_aibj, c_ai)
!!
!!    Jacobian transpose CC2 A2
!!    Written by Sarai D. Folkestad and Alexander C. Paul, Feb 2019
!!
!!    Calculates the A2 term,
!!
!!    sigma_aibj =+ (2F_jb c_ai - F_ib c_aj - L_ikjb c_ak + L_cajb c_ci)
!!
!!    Modified by Tor S. Haugland, Nov 2019
!!
!!    Now uses BLAS zgeru for outer-product instead of for-loops.
!!
      implicit none
!
      class(doubles) :: wf
      complex(dp), dimension(wf%n_v, wf%n_o), intent(in)                      :: c_ai
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout)   :: sigma_aibj
!
   end subroutine jacobian_transpose_doubles_a2_doubles_complex
