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
   module subroutine effective_jacobian_transpose_transformation_cc3(wf, omega, c)
!!
!!    Jacobian transpose transformation (CC3)
!!    Alexander Paul and Rolf H. Myhre, March 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: c
!
   end subroutine effective_jacobian_transpose_transformation_cc3
!
!
   module subroutine jacobian_transpose_cc3_sigma1_t3_A1_cc3(wf, c_abij, sigma_ai)
!!
!!    Computes first contribution of the T3 amplitudes to sigma_1
!!
!!    Reads in the intermediates X_abid and Y_akil prepared in prepare_jacobian_transpose
!!    contracts with c_abij and adds to sigma_ai
!!
!!    sigma_dl =  sum_abi X_abid * C_abil + sum_aik C_daki * Y_akil
!!    
!!    Written by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: c_abij
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
!
   end subroutine jacobian_transpose_cc3_sigma1_t3_A1_cc3
!
!
   module subroutine jacobian_transpose_cc3_sigma1_t3_B1_cc3(wf, c_abij, sigma_ai)
!!
!!    Computes first contribution of the T3 amplitudes to sigma_1
!!
!!    Constructs t^abc_ijk for fixed ijk and contracts with c_abij
!!    The intermediate X_ck is then contracted with L_kcld
!!
!!    sigma_dl =  sum_abcijk C^ab_ij (t^abc_ijk - t^acb_ijk) L_kcld
!!    
!!    Written by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: c_abij
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
!
   end subroutine jacobian_transpose_cc3_sigma1_t3_B1_cc3
!
!
   module subroutine jacobian_transpose_cc3_X_ck_calc_cc3(wf, i, j, k, t_abc, u_abc, X_ck, c_abij)
!!
!!    Constructs the intermediate X_ck
!!
!!    X_ck =  sum_abcijk C^ab_ij (t^abc_ijk - t^acb_ijk)
!!
!!    All permutations for i,j,k have to be considered due to the restrictions in the i,j,k loops
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)           :: t_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: u_abc
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)                :: X_ck
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in)   :: c_abij
!
   end subroutine jacobian_transpose_cc3_X_ck_calc_cc3
!
!
   module subroutine jacobian_transpose_cc3_C3_terms_cc3(wf, omega, c_ai, c_abij, sigma_ai, sigma_abij)
!!
!!    Construct C^abc_ijk in single batches of ijk and compute the contributions
!!    to the singles and doubles part of the outgoing vector
!!
!!    The construction of C3 is split into contributions 
!!    from outer products and matrix multiplications
!!
!!    1 array for each Permutation of C_abc will be used 
!!    to reduce the amount of N^7-contractions and sorting
!!
!!    c_μ3 = (ω - ε^abc_ijk)^-1 (c_μ1 < μ1 | [H,tau_ν3] | R > + c_μ2 < μ2 | [H,tau_ν3] | R >
!!
!!    σ1 += c_μ3 < μ3 | [[H,T_2],tau_ν1] | R >
!!    σ2 += c_μ3 < μ3 | [H,tau_ ν2] | R >
!!
!!    Written by Alexander Paul and Rolf H. Myhre, April 2019
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
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout) :: sigma_abij
!
   end subroutine jacobian_transpose_cc3_C3_terms_cc3
!
!
   module subroutine jacobian_transpose__cc3_calc_c3_matmul_cc3(wf, i, j, k, c_abij, c_abc, c_bac, & 
                                                               c_cba, c_acb, c_cab, c_bca, u_abc,  &
                                                               g_dbic, g_dbjc, g_dbkc,             &
                                                               g_jlic, g_klic, g_kljc,             &
                                                               g_iljc, g_ilkc, g_jlkc)
!!
!!    Calculate the contributions from matrix multiplications 
!!    to the  C3 amplitudes for fixed indices i,j,k
!!
!!    C^abc_ijk 
!!    = (ω - ε^abc_ijk)^-1 P^abc_ijk (C_ai*L_jbkc - C_ak*L_jbic + Cabij*F_kc - C_abik*F_jc
!!    + sum_l (C_ablk g_iljc - C_abil L_jlkc) - sum_d (C_adjk g_ibdc - C_adij L_dbkc)
!!
!!    Contibutions:
!!    sum_l (C_ablk g_iljc + C_abil g_jckl - 2 C_abil g_jlkc) 
!!    - sum_d (C_adjk g_ibdc + C_adij g_dckb - 2 C_adij g_dbkc)
!!
!!    Written by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: c_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: c_bac
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: c_cba
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: c_acb
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: c_cab
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: c_bca
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: u_abc
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: c_abij
!
!     g_dbkc ordered bcd,k
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in) :: g_dbic
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in) :: g_dbjc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in) :: g_dbkc
!
!     g_jlkc ordered cl,jk
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_jlic
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_klic
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_kljc
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_iljc
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_ilkc
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_jlkc
!
   end subroutine jacobian_transpose__cc3_calc_c3_matmul_cc3
!
!