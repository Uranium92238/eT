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
   module subroutine prepare_for_jacobian_cc2(wf)
!!
!!    Prepare for Jacobian
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2019
!!
      implicit none 
!
      class(cc2), intent(inout) :: wf 
!
   end subroutine prepare_for_jacobian_cc2
!
!
   module subroutine jacobian_transform_trial_vector_cc2(wf, c_i)
!!
!!    Jacobian transform trial vector 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep 2018
!!
      class(cc2), intent(in) :: wf 
!
      real(dp), dimension(wf%n_es_amplitudes) :: c_i
!
   end subroutine jacobian_transform_trial_vector_cc2
!
!
   module subroutine jacobian_cc2_transformation_cc2(wf, c)
!!
!!    Jacobian transformation (cc2)
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, Jan 2019
!!
!!    Directs the transformation by the cc2 Jacobi matrix,
!!
!!       A_mu,nu = < mu | exp(-T) [H, tau_nu] exp(T) | R >,
!!
!!    where the basis employed for the brackets is biorthonormal.
!!    The transformation is rho = A c, i.e.,
!!
!!       rho_mu = (A c)_mu = sum_ck A_mu,ck c_ck
!!                  + 1/2 sum_ckdl A_mu,ckdl c_ckdl (1 + delta_ck,dl).
!!
!!    On exit, c is overwritten by rho. That is, c_a_i = rho_a_i,
!!    and c_aibj = rho_aibj.
!!
      implicit none
!
      class(cc2) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes) :: c
!
   end subroutine jacobian_cc2_transformation_cc2
!
!
   module subroutine jacobian_cc2_a1_cc2(wf, rho_ai, c_ai)
!!
!!    Jacobian CC2 A1
!!    Written by Sarai D. Folkestad Eirik F. Kjønstad Jan 2019
!!
!!    rho_ai^A1 = sum_bjck (L_kcjb u_aick c_bj - g_kcjb (u_ckbi c_aj + u_ckaj c_bi))
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(out)   :: rho_ai   
!
   end subroutine jacobian_cc2_a1_cc2
!
!
   module subroutine jacobian_cc2_b1_cc2(wf, rho_ai, c_aibj)
!!
!!    Jacobian CC2 B1
!!    Written by Sarai D. Folkestad Eirik F. Kjønstad Jan 2019
!!
!!    rho_ai^B1 = 2 sum_bj F_jb c_aibj - F_jb c_ajbi ) - sum_kjb L_kijb c_akbj
!!
      implicit none
!
      class(cc2) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)   :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(out)                  :: rho_ai   
!
   end subroutine jacobian_cc2_b1_cc2
!
!
   module subroutine jacobian_cc2_a2_cc2(wf, rho_aibj, c_ai)
!!
!!    Jacobian CC2 A2
!!    Written by Sarai D. Folkestad Eirik F. Kjønstad Jan 2019
!!
!!    rho_aibj^A2 = /(1/Δ_aibj)P_aibj sum_c g_abkc c_cj - sum_k g_aikj c_bk, 
!!
      implicit none
!
      class(cc2) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                   :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(out)  :: rho_aibj   
!
   end subroutine jacobian_cc2_a2_cc2
!
!
   module subroutine jacobian_cc2_b2_cc2(wf, rho_aibj, c_aibj)
!!
!!    Jacobian CC2 B2
!!    Written by Sarai D. Folkestad Eirik F. Kjønstad Jan 2019
!!
!!    rho_aibj^B2 = ε_aibj c_aibj/(1/Δ_aibj) 
!!
      implicit none
!
      class(cc2) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout)   :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(out)  :: rho_aibj   
!
   end subroutine jacobian_cc2_b2_cc2
!
