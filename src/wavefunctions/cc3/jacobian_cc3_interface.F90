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
   module subroutine effective_jacobian_transformation_cc3(wf, omega, c, rho)
!!
!!    Effective Jacobian transformation (CC3)
!!    Written by Alexander C. Paul and Rolf H. Myhre, Feb 2019
!!
!!    Directs the transformation by the CC3 Jacobi matrix,
!!
!!       A_mu,nu = < mu| exp(-T) [H, tau_nu] exp(T) | R >,
!!
!!    where the basis employed for the brackets is biorthonormal.
!!    The transformation is rho = A c, i.e.,
!!
!!       rho_mu = (A c)_mu 
!!       = sum_ck A_mu,ck c_ck + 1/2 sum_ckdl A_mu,ckdl c_ckdl (1 + delta_ck,dl)
!!
!!    On exit, c is overwritten by rho. That is, c(ai) = rho_a_i,
!!    and c(aibj) = rho_aibj.
!!
      implicit none
!
      class(cc3) :: wf
      real(dp), intent(in) :: omega
      real(dp), dimension(wf%n_t1 + wf%n_t2), intent(in)  :: c
      real(dp), dimension(wf%n_t1 + wf%n_t2), intent(out) :: rho
!
   end subroutine effective_jacobian_transformation_cc3
!
!
   module subroutine jacobian_cc3_t3_a2_cc3(wf, c_ai, rho_abij)
!!
!!    Jacobian CC3 A2
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
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
      implicit none
!
      class(cc3) :: wf
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(out) :: rho_abij
!
   end subroutine jacobian_cc3_t3_a2_cc3
!
!
   module subroutine jacobian_cc3_t3_b2_cc3(wf, c_ai, rho_abij)
!!
!!    Jacobian CC3 B2
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!!    sigma_abij +=  sum_ckdl C^d_l L_kcld (t^abc_ijk - t^bac_ijk)
!!               +=  sum_ck F_kc_c1 * (t^abc_ijk - t^bac_ijk)
!!
!!    Constructs t^abc_ijk for fixed ijk contracts with
!!    the c1-transformed Fock matrix
!!
      implicit none
!
      class(cc3) :: wf
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout) :: rho_abij
!
   end subroutine jacobian_cc3_t3_b2_cc3
!
!
   module subroutine jacobian_cc3_c3_a_cc3(wf, omega, c_ai, c_abij, rho_ai, rho_abij)
!!
!!    Contribution of the C3/R3 terms
!!    Written by Alexander C. Paul and Rolf H. Myhre, Feb 2019
!!
!!    Construct C^abc_ijk in single batches of ijk and compute the contributions
!!    to the singles and doubles part of the outgoing vector
!!
!!    The triples amplitudes are expressed in terms of doubles amplitudes:
!!    C_3 = (omega - epsilon_mu3)^-1 (< mu3| [H,C_2] | HF > 
!!                                  + < mu3| [[H,C_1],T_2] |HF >)
!!
!!    c^abc = (omega - epsilon^abc_ijk)^-1 * P^abc_ijk 
!!             (sum_d c^ad_ij g_ckbd - sum_l c^ab_il g_cklj
!!            + sum_d t^ad_ij g'_bdck - sum_l t^ab_il g'_cklj
!!
!!    rho1 += < mu1| [H,C_3] |R >
!!    rho2 += < mu2| [H,C_3] |R >
!!
!!    Based on omega_cc3_a_cc3 written by Rolf H. Myhre
!!
      implicit none
!
      class(cc3) :: wf
      real(dp), intent(in) :: omega
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: c_abij
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout) :: rho_abij
!
   end subroutine jacobian_cc3_c3_a_cc3
!
!
   module subroutine construct_c1_integrals_cc3(wf, c_ai)
!!
!!    Construct c1 transformed integrals
!!    Alexander C. Paul and Rolf H. Myhre Feb 2019
!!
!!    g'_bdck = (b'd|ck) + (bd|c'k) + (bd|ck')   ordered as dbc,k
!!    g'_ljck = (lj'|ck) + (lj|ck') + (lj|c'k)   ordered as lc,jk
!!
!!    Based on omega_cc3_integrals_cc3 written by Rolf H. Myhre
!!
      implicit none
!
      class(cc3) :: wf
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
!
   end subroutine construct_c1_integrals_cc3
!
!
   module subroutine construct_c1_fock_cc3(wf, c_ai, F_ia_c1)
!!
!!    Construct C1-transformed fock matrix ov-block
!!    Written by Alexander C. Paul, Feb 2019
!!
!!    Calculates C1-transformed occupied-virtual elements of the Fock matrix
!!    required for the CC3 jacobian and returns it ordered as n_v, n_o
!!
!!    F_ia_c1 = sum_j L_iajj' = sum_j 2 g_iajj' - g_ij'ja
!!
      implicit none
!
      class(cc3) :: wf
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(out) :: F_ia_c1
!
   end subroutine construct_c1_fock_cc3
!
!
   module subroutine jacobian_cc3_b2_fock_cc3(wf, i, j, k, t_abc, u_abc, rho_abij, F_ov_ck)
!!
!!    Jacobian CC3 contribution c1-transformed fock matrix
!!    Written by Alexander C. Paul and Rolf H. Myhre, Feb 2019
!!
!!    rho_2 =+ P^{ab}_{ij} sum_kc (t^abc_ijk - t^cba_ijk) F_kc
!!
!!    The permutations of i,j,k are necessary 
!!    due to the index restrictions in the batching loops
!!
      implicit none
!
      class(cc3) :: wf
      integer, intent(in) :: i, j, k
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)              :: t_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)             :: u_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout)   :: rho_abij
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                      :: F_ov_ck
!
   end subroutine jacobian_cc3_b2_fock_cc3
