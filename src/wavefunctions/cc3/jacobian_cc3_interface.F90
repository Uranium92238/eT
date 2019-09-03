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
   module subroutine prepare_for_jacobian_cc3(wf)
!!
!!    Prepare for jacobian
!!    Written by Rolf Heilemann Myhre, April 2019
!!
      implicit none
!
      class(cc3), intent(inout) :: wf
!
   end subroutine prepare_for_jacobian_cc3
!
!
   module subroutine effective_jacobian_transformation_cc3(wf, omega, c)
!!
!!    Jacobian transformation (CC3)
!!    Alexander Paul and Rolf H. Myhre Feb 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: c
!
   end subroutine effective_jacobian_transformation_cc3
!
!
   module subroutine jacobian_cc3_t3_a2_cc3(wf, c_ai, rho_abij)
!!
!!    Computes the first contribution of the T3 amplitudes to rho_2
!!
!!    Reads in the intermediates X_abid and X_ajil prepared in 
!!    prepare_jacobian_transform contracts with c_ai and adds to rho_abij
!!
!!    rho_abil += sum_abi X_abid * C_dl
!!    rho_daji += sum_aik C_dl * X_ajil
!!
!!    where: X_abid = - sum_jck (2 t^abc_ijk - t^cba_ijk - t^acb_ijk) * g_kcjd
!!           X_ajil = - sum_bck (2 t^abc_ijk - t^cba_ijk - t^acb_ijk) * g_lbkc
!!
!!    Written by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(out) :: rho_abij
!
   end subroutine jacobian_cc3_t3_a2_cc3
!
!
   module subroutine jacobian_cc3_t3_b2_cc3(wf, c_ai, rho_abij)
!!
!!    Computes the second contribution of the T3 amplitudes to rho_2
!!
!!    sigma_abij +=  sum_ckdl C^d_l L_kcld (t^abc_ijk - t^bac_ijk)
!!               +=  sum_ck F_kc_c1 * (t^abc_ijk - t^bac_ijk)
!!
!!    Constructs t^abc_ijk for fixed ijk contracts with 
!!    the c1-transformed Fock matrix
!!    
!!    Written by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout) :: rho_abij
!
   end subroutine jacobian_cc3_t3_b2_cc3
!
!
   module subroutine jacobian_cc3_c3_a_cc3(wf, omega, c_ai, c_abij, rho_ai, rho_abij)
!!
!!    Construct C^abc_ijk in single batches of ijk and compute the contributions
!!    to the singles and doubles part of the outgoing vector
!!
!!    The triples amplitudes are expressed in terms of doubles amplitudes:
!!    C_3 = (omega - ε^abc_ijk)^-1 (< μ3 | [H,C_2] | HF > + < μ3 | [[H,C_1],T_2] | HF >)
!!
!!    rho1 += < μ1 | [H,C_3] | R >
!!    rho2 += < μ2 | [H,C_3] | R >
!!
!!    Based on omega_cc3_a_cc3 written by Rolf H. Myhre
!!    Modified by Alexander Paul and Rolf H. Myhre, Feb 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: c_abij
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout) :: rho_abij
!
   end subroutine jacobian_cc3_c3_a_cc3
!
!
   module subroutine jacobian_cc3_A_cc3(wf, omega, c_ai, c_abij, rho_ai, rho_abij)
!!
!!
!!    CC3 jacobian terms
!!
!!    The triples amplitudes are expressed in terms of doubles amplitudes:
!!    C_3 = (omega - ε^abc_ijk)^-1 (< mu3 | [H,C_2] | HF > + < mu3 | [[H,C_1],T_2] | HF >)
!!    T_3 = (omega - ε^abc_ijk)^-1 < mu3 | [H,T_2] | HF >
!!
!!    They are then used to compute the contributions 
!!    to the singles and doubles part of the transformed vector
!!
!!    rho1 = rho1(CCSD) + < mu1 | [H,C_3] | HF >
!!    rho2 = rho2(CCSD) + < mu2 | [H,C_3] | HF > + < mu2 | [[H,C_1],T_3] | HF >
!!
!!    Alex C. Paul and Rolf H. Myhre, Feb 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: c_abij
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout) :: rho_abij
!
   end subroutine jacobian_cc3_A_cc3
!
!
   module subroutine construct_c1_integrals_cc3(wf, c_ai)
!!
!!    Construct c1 transformed integrals needed in CC3 Jacobian
!!    Alexander Paul and Rolf H. Myhre Feb 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
!
   end subroutine construct_c1_integrals_cc3
!
!
   module subroutine construct_c1_fock_cc3(wf, c_ai, F_ia_c1)
!!
!!    Calculates C1 transformed elements of the Fock matrix required for the CC3 jacobian
!!    Rolf H. Myhre and Alexander Paul, Feb 2019
!!
!!    F_ia_c1 = sum_j L_iajj' = sum_j 2 g_iajj' - g_ij'ja
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(out) :: F_ia_c1
!
   end subroutine construct_c1_fock_cc3
!
!
   module subroutine jacobian_cc3_b2_fock_cc3(wf, i, j, k, t_abc, u_abc, rho_abij, F_kc)
!!
!!    Calculate the Fock contribution to rho2 for fixed i,j and k
!!
!!    rho_abji =+ P^ab_ij sum_kc (C^abc_ijk - C^cba_ijk) F_kc
!!
!!    Alexander Paul and Rolf H. Myhre, Feb 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)              :: t_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)             :: u_abc
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout)   :: rho_abij
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                      :: F_kc
!
   end subroutine jacobian_cc3_b2_fock_cc3
!
!
