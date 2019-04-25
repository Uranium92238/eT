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
!!    The intermediate X_ai is then contracted with L_kcld
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
   module subroutine jacobian_transpose_cc3_X_ai_calc_cc3(wf, i, j, k, t_abc, u_abc, X_ai, c_bcjk)
!!
!!    Constructs the intermediate X_ai
!!
!!    X_ai =  sum_bcjk (t^abc_ijk - t^bac_ijk)*C^bc_jk
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
      real(dp), dimension(wf%n_v, wf%n_o), intent(out)                  :: X_ai
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in)   :: c_bcjk
!
   end subroutine jacobian_transpose_cc3_X_ai_calc_cc3
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
   module subroutine jacobian_transpose_cc3_calc_outer_cc3(wf, i, j, k, c_ai, c_abij,     & 
                                                            c_abc, c_bac, c_cba, c_acb,   &
                                                            c_bca, u_abc, F_kc,           &
                                                            L_jbic, L_kbic, L_jbkc)
!!
!!    Calculate the contributions from outer products 
!!    to the  C3 amplitudes for fixed indices i,j,k
!!
!!    C^abc_ijk 
!!    = (ω - ε^abc_ijk)^-1 P^abc_ijk (C_ai*L_jbkc - C_ak*L_jbic + Cabij*F_kc - C_abik*F_jc)
!!    + sum_l (C_ablk g_iljc - C_abil L_jlkc) - sum_d (C_adjk g_ibdc - C_adij L_dbkc)
!!
!!    Contibutions in this routine:
!!    P^abc_ijk (C_ai*L_jbkc - C_ak*L_jbic + Cabij*F_kc - C_abik*F_jc)
!!
!!    L_jlkc and L_dbkc split up to reduce the amount of N^7 contractions
!!    but 6 arrays for c_abc needed (for all permutations of abc)
!!
!!    Written by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: c_abij
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: c_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: c_bac
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: c_cba
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: c_acb
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: c_bca
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: u_abc
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: F_kc
!
      real(dp), dimension(wf%n_v, wf%n_v), intent(in) :: L_jbic
      real(dp), dimension(wf%n_v, wf%n_v), intent(in) :: L_kbic
      real(dp), dimension(wf%n_v, wf%n_v), intent(in) :: L_jbkc
!
   end subroutine jacobian_transpose_cc3_calc_outer_cc3
!
!
   module subroutine jacobian_transpose_cc3_calc_c3_matmul_cc3(wf, i, j, k, c_abij, c_abc, c_bac, & 
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
   end subroutine jacobian_transpose_cc3_calc_c3_matmul_cc3
!
!
   module subroutine jacobian_transpose_cc3_collect_c3_cc3(wf, omega, i, j, k, c_abc, c_bac, & 
                                                            c_cba, c_acb, c_cab, c_bca)
!!
!!    Adds up the contributions from all permutations of the indices abc to c_abc
!!    from the matrix multiplications and outer products
!!
!!    Divides by (ω - ε^abc_ijk)
!!
!!    Written by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), intent(in) :: omega
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
   end subroutine jacobian_transpose_cc3_collect_c3_cc3
!
!
   module subroutine jacobian_transpose_cc3_sigma2_cc3(wf, i, j, k, c_abc, u_abc, sigma_abij,   &
                                                      g_bdci, g_bdcj, g_bdck, g_ljci, g_lkci,   &
                                                      g_lkcj, g_licj, g_lick, g_ljck)
!!
!!    Calculates triples contribution to sigma2
!!
!!    sigma_adij =   sum_ckd c^abc_ijk g_bdck
!!    sigma_abil = - sum_cki c^bac_ijk g_lick
!!
!!    All permutations for i,j,k have to be considered due to the restrictions in the i,j,k loops   
!!
!!    Written by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)              :: c_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)             :: u_abc
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout)   :: sigma_abij
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)              :: g_bdci
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)              :: g_bdcj
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)              :: g_bdck
!
      real(dp), dimension(wf%n_o, wf%n_v), intent(in)                      :: g_ljci
      real(dp), dimension(wf%n_o, wf%n_v), intent(in)                      :: g_lkci
      real(dp), dimension(wf%n_o, wf%n_v), intent(in)                      :: g_lkcj
      real(dp), dimension(wf%n_o, wf%n_v), intent(in)                      :: g_licj
      real(dp), dimension(wf%n_o, wf%n_v), intent(in)                      :: g_lick
      real(dp), dimension(wf%n_o, wf%n_v), intent(in)                      :: g_ljck
!
   end subroutine jacobian_transpose_cc3_sigma2_cc3
!
!
   module subroutine construct_intermediates_c3_cc3(wf, i, j, k, c_abc, u_abc, t_abij, Y_cmjk,   &
                                                   X_bcei, X_bcej, X_bcek)
!!
!!    Constructs the intermediates X_bcei and Y_cmjk used to compute the c3 contributions to sigma_ai
!!
!!    X_bcei = sum_aij c^abc_ijk * t^ae_ij
!!    Y_cmjk = sum_abj c^bac_ijk * t^ba_mj
!!
!!    All permutations for i,j,k have to be considered due to the restrictions in the i,j,k loops
!!
!!    Written by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)              :: c_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)             :: u_abc
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in)      :: t_abij
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_o, wf%n_o), intent(inout)   :: Y_cmjk
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(inout)           :: X_bcei
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(inout)           :: X_bcej
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(inout)           :: X_bcek
!
   end subroutine construct_intermediates_c3_cc3
!
!
   module subroutine jacobian_transpose_cc3_sigma1_C3_A1_cc3(wf, sigma_ai, Y_cmjk)
!!
!!    Computes the contribution of the intermediate Y_cmjk to sigma_1
!!
!!    sigma1 += sum_mjk Y_cmjk * g_mjlk
!!    sigma1 += sum_cmj g_mjcd * Y_cmjk
!!    sigma1 += sum_cmk g_leck * Y_cmjk
!!    
!!    Written by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_o, wf%n_o), intent(in) :: Y_cmjk
!
   end subroutine jacobian_transpose_cc3_sigma1_C3_A1_cc3
!
!
   module subroutine jacobian_transpose_cc3_sigma1_C3_B1_cc3(wf, sigma_ai)
!!
!!    Computes the contribution of the intermediate X_bcek to sigma_1
!!
!!    sigma1 += sum_bec g_becd * X_bcek
!!    sigma1 += sum_cek X_bcek * g_leck
!!    sigma1 += sum_bek X_bcek * g_lkbe
!!    
!!    Written by Alexander Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
!
   end subroutine jacobian_transpose_cc3_sigma1_C3_B1_cc3
