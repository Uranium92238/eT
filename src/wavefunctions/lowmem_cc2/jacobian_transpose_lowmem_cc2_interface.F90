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
   module subroutine effective_jacobian_transpose_transformation_lowmem_cc2(wf, omega, b, cvs)
!!
!!    Effective Jacobian transpose transformation
!!    Written by Sarai Dery Folkestad, Jun 2019
!!
      implicit none
!
      class(lowmem_cc2) :: wf
      real(dp), intent(in) :: omega
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: b
      logical, intent(in) :: cvs
!
   end subroutine effective_jacobian_transpose_transformation_lowmem_cc2
!
!
   module subroutine jacobian_transpose_cc2_a1_lowmem_cc2(wf, sigma_ai, b_ai, eps_o, eps_v)
!!
!!    Jacobian transpose A1 (CC2)
!!    Written by Sarai D. Folkestad, Jun 2019
!!
!!    Constructs the Jacobian transpose A1 term
!!
!!       A1 = sum_ckbj u_bjck L_iakc b_bj - t_bjck L_kcja b_bi
!!
!!    where
!!
!!       u_bjck = 2t_bjck - t_bkcj
!!
!!    and 
!!
!!       t_bjck = - g_bjck/ε_bjck
!!
!!    Batching over j and k, we will construct the intermediates
!!
!!       X_ck = sum_bj u_bjck b_bj = - 2 sum_bj g_bjck b_bj / ε_bjck
!!                                       + sum_bj g_cjbk b_bj / ε_bjck
!!
!!       Y_ba = sum_bj t_bjck L_kcja = - sum_bj g_bjck L_kcja / ε_bjck
!!
      implicit none
!
      class(lowmem_cc2), intent(in) :: wf
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: b_ai
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
   end subroutine jacobian_transpose_cc2_a1_lowmem_cc2
!
!
   module subroutine jacobian_transpose_cc2_b1_lowmem_cc2(wf, sigma_ai, b_ai, eps_o, eps_v)
!!
!!    Jacobian transpose B1 (CC2)
!!    Written by Sarai D. Folkestad, Jun 2019
!!
!!    Constructs the Jacobian transpose A1 term
!!
!!       B1 =  - sum_ckbj t_bjck L_kcib b_aj
!!
!!    where
!!
!!       t_bjck = - g_bjck/ε_bjck
!!
!!    Batching over b and c, we will construct the intermediate
!!
!!       X_ij = sum_bck t_ckbj L_ibkc
!!
      implicit none
!
      class(lowmem_cc2), intent(in) :: wf
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: b_ai
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
   end subroutine jacobian_transpose_cc2_b1_lowmem_cc2
!
!
   module subroutine effective_jacobian_transpose_cc2_a1_lowmem_cc2(wf, omega, sigma_ai, b_ai, eps_o, eps_v)
!!
!!    Effective jacobian transpose CC2 A1
!!    Written by Sarai D. Folkestad, Jul 2019
!!
!!    sigma_ai =+ g_bjca [P_bj,ci( 2 F_jb b_ci - F_ib b_cj)] 1/(ω - ε_cibj)
!!
!!    The term is calculated in batches over the b and c.
!!
      implicit none
!
      class(lowmem_cc2), intent(in) :: wf
      real(dp), intent(in) :: omega
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: b_ai
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
   end subroutine effective_jacobian_transpose_cc2_a1_lowmem_cc2
!
!
   module subroutine effective_jacobian_transpose_cc2_b1_lowmem_cc2(wf, omega, sigma_ai, b_ai, eps_o, eps_v)
!!
!!    Effective jacobian transpose CC2 B1
!!    Written by Sarai D. Folkestad, Jul 2019
!!
!!    sigma_ai =+ g_bjca [P_bj,ci L_dcjb b_di] 1/(ω - ε_cibj)
!!
!!    The term is calculated in batches over the b and c.
!!
      implicit none
!
      class(lowmem_cc2), intent(in) :: wf
      real(dp), intent(in) :: omega
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: b_ai
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
   end subroutine effective_jacobian_transpose_cc2_b1_lowmem_cc2
!
!
   module subroutine effective_jacobian_transpose_cc2_c1_lowmem_cc2(wf, omega, sigma_ai, b_ai, eps_o, eps_v)
!!
!!    Effective jacobian transpose CC2 C1
!!    Written by Sarai D. Folkestad, Jul 2019
!!
!!    sigma_ai =- g_bjca [P_bj,ci L_ikjb b_ck] 1/(ω - ε_cibj)
!!
!!    The term is calculated in batches over the k, b, c
!!
      implicit none
!
      class(lowmem_cc2), intent(in) :: wf
      real(dp), intent(in) :: omega
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: b_ai
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
   end subroutine effective_jacobian_transpose_cc2_c1_lowmem_cc2
!
!
   module subroutine effective_jacobian_transpose_cc2_d1_lowmem_cc2(wf, omega, sigma_ai, b_ai, eps_o, eps_v)
!!
!!    Effective jacobian transpose CC2 D1
!!    Written by Sarai D. Folkestad, Jul 2019
!!
!!    sigma_ai =- g_bjik [P_bj,ak( 2 F_jb b_ak - F_kb b_aj)] 1/(ω - ε_bjak)
!!
!!    The term is calculated in batches over the k and j.
!!
      implicit none
!
      class(lowmem_cc2), intent(in) :: wf
      real(dp), intent(in) :: omega
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: b_ai
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
   end subroutine effective_jacobian_transpose_cc2_d1_lowmem_cc2
!
!
   module subroutine effective_jacobian_transpose_cc2_e1_lowmem_cc2(wf, omega, sigma_ai, b_ai, eps_o, eps_v)
!!
!!    Effective jacobian transpose CC2 E1
!!    Written by Sarai D. Folkestad, Jul 2019
!!
!!    sigma_ai =+ g_bjik [P_bj,ak( L_jbkl b_al)] 1/(ω - ε_bjak)
!!
!!    The term is calculated in batches over the k and j.
!!
      implicit none
!
      class(lowmem_cc2), intent(in) :: wf
      real(dp), intent(in) :: omega
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: b_ai
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
   end subroutine effective_jacobian_transpose_cc2_e1_lowmem_cc2
!
!
   module subroutine effective_jacobian_transpose_cc2_f1_lowmem_cc2(wf, omega, sigma_ai, b_ai, eps_o, eps_v)
!!
!!    Effective jacobian transpose CC2 F1
!!    Written by Sarai D. Folkestad, Jul 2019
!!
!!    sigma_ai =- g_bjik [P_bj,ak( L_cajb b_ck)] 1/(ω - ε_bjak)
!!
!!    The term is calculated in batches over the k, j and c indices.
!!
      implicit none
!
      class(lowmem_cc2), intent(in) :: wf
      real(dp), intent(in) :: omega
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: b_ai
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
   end subroutine effective_jacobian_transpose_cc2_f1_lowmem_cc2
