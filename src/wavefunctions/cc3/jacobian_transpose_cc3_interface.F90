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
   module subroutine effective_jacobian_transpose_transformation_cc3(wf, omega, b, sigma, cvs)
!!
!!    Effective Jacobian transpose transformation (CC3)
!!    Alexander C. Paul and Rolf H. Myhre, March 2019
!!
!!    Directs the transformation by the transpose of the  CC3 Jacobi matrix,
!!
!!       A_mu,nu = < mu| exp(-T) [H, tau_nu] exp(T)|R >,
!!
!!    The transformation is performed as sigma^T = c^T A, where c is the vector
!!    sent to the routine. 
!!
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      class(cc3) :: wf
      real(dp), intent(in) :: omega
      real(dp), dimension(wf%n_t1 + wf%n_t2), intent(in)  :: b
      real(dp), dimension(wf%n_t1 + wf%n_t2), intent(out) :: sigma
      logical, intent(in) :: cvs
!
   end subroutine effective_jacobian_transpose_transformation_cc3
!
!
   module subroutine jacobian_transpose_cc3_t3_a1_cc3(wf, c_abij, sigma_ai)
!!
!!    Jacobian transpose T3 A1 term
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!!    Computes the first contribution of the T3 amplitudes to sigma_1
!!
!!    Reads in the intermediates X_abid and X_ajil prepared in prepare_jacobian_transpose
!!    contracts with c_abij and adds to sigma_ai
!!
!!    sigma_dl +=  sum_abi X_abid * C_abil + sum_aik C_daji * X_ajil
!!
!!    where: X_abid = - sum_jck (2t^abc_ijk - t^cba_ijk - t^acb_ijk) * g_kcjd
!!           X_ajil = - sum_bck (2t^abc_ijk - t^cba_ijk - t^acb_ijk) * g_lbkc
!!    
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      class(cc3) :: wf
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: c_abij
      real(dp), dimension(wf%n_v, wf%n_o), intent(out) :: sigma_ai
!
   end subroutine jacobian_transpose_cc3_t3_a1_cc3
!
!
   module subroutine jacobian_transpose_cc3_t3_b1_cc3(wf, c_abij, sigma_ai, cvs)
!!
!!    Jacobian transpose T3 B1 term
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!!    Constructs t^abc_ijk for fixed ijk and contracts with c_abij
!!    The intermediate X_ai is then contracted with L_iald
!!
!!    sigma_dl +=  sum_abcijk C^bc_jk (t^abc_ijk - t^bac_ijk) L_iald
!!             +=  sum_ai X_ai * L_iald
!!
      implicit none
!
      class(cc3) :: wf
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: c_abij
      logical, intent(in) :: cvs
!
   end subroutine jacobian_transpose_cc3_t3_b1_cc3
!
!
   module subroutine construct_x_ai_intermediate_cc3(wf, i, j, k, t_abc, u_abc, x_ai, c_bcjk)
!!
!!    Constructs the intermediate X_ai
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!!    X_ai = sum_abcijk (t^abc_ijk - t^bac_ijk) C^bc_jk 
!!
!!    All permutations for i,j,k have to be considered due to the restrictions in the i,j,k loops
!!
      implicit none
!
      class(cc3) :: wf
      integer, intent(in) :: i, j, k
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)           :: t_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: u_abc
      real(dp), dimension(wf%n_v, wf%n_o), intent(out)                  :: X_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in)   :: c_bcjk
!
   end subroutine construct_x_ai_intermediate_cc3
!
!
   module subroutine jacobian_transpose_cc3_c3_a_cc3(wf, omega,            &
                                                     c_ai, c_abij,         &
                                                     sigma_ai, sigma_abij, &
                                                     cvs)
!!
!!    Contributions of the c3/L3 terms
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
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
!!    c_mu3 = (omega - epsilon^abc_ijk)^-1 (c_mu1 < mu1| [H,tau_nu3] |R > 
!!                                        + c_mu2 < mu2| [H,tau_nu3] |R >)
!!
!!    sigma_1 += c_mu3 < mu3| [[H,T_2],tau_nu1] |R >
!!    sigma_2 += c_mu3 < mu3| [H,tau_ nu2] |R >
!!
      implicit none
!
      class(cc3) :: wf
      real(dp), intent(in) :: omega
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: c_abij
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout) :: sigma_abij
      logical, intent(in) :: cvs
!
   end subroutine jacobian_transpose_cc3_c3_a_cc3
!
!
   module subroutine jacobian_transpose_cc3_c3_calc_cc3(wf, i, j ,k, c_ai, c_abij, &
                                                        c_abc, u_abc,              &
                                                        v_abc, F_ov_ck,               &
                                                        L_ibjc, L_ibkc, L_jbkc,    &
                                                        g_dbic, g_dbjc, g_dbkc,    &
                                                        g_jlic, g_klic, g_kljc,    &
                                                        g_iljc, g_ilkc, g_jlkc)
!!
!!    Construct c3 amplitudes for fixed i,j,k
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!!    C^abc_ijk 
!!    = (omega - epsilon^abc_ijk)^-1 P^abc_ijk (C__ai*L_jbkc - C__ak*L_jbic 
!!                                  + C_abij*F_kc - C__abik*F_jc
!!                                  + sum_l (C__ablk g_iljc - C__abil L_jlkc) 
!!                                  - sum_d (C__adjk g_ibdc - C__adij L_dbkc))
!!
!!    Contributions from outer products:
!!    P^abc_ijk (C_ai*L_jbkc - C_ak*L_jbic + Cabij*F_kc - C_abik*F_jc)
!!
!!    Contributions from matrix multiplication:
!!      sum_l P^abc_ijk (C_ablk g_iljc + C_abil g_jckl - 2 C_abil g_jlkc) 
!!    - sum_d P^abc_ijk (C_adjk g_ibdc + C_adij g_dckb - 2 C_adij g_dbkc)
!!
      implicit none
!
      class(cc3) :: wf
      integer, intent(in) :: i, j, k
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: c_abij
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: c_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: u_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: v_abc
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: F_ov_ck
      real(dp), dimension(wf%n_v, wf%n_v), intent(in) :: L_ibjc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in) :: L_ibkc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in) :: L_jbkc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in) :: g_dbic
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in) :: g_dbjc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in) :: g_dbkc
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_jlic
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_klic
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_kljc
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_iljc
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_ilkc
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_jlkc
!
   end subroutine jacobian_transpose_cc3_c3_calc_cc3
!
!
   module subroutine jacobian_transpose_cc3_a_n7_cc3(wf, i, j, k, c_abc, u_abc, sigma_abij,  &
                                                     g_bdci, g_bdcj, g_bdck, g_ljci, g_lkci, &
                                                     g_lkcj, g_licj, g_lick, g_ljck)
!!
!!    Jacobian transpose A2 term
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!!    sigma_adij =   sum_ckd c^abc_ijk g_bdck
!!    sigma_abil = - sum_cki c^bac_ijk g_lick
!!
!!    All permutations for i,j,k have to be considered 
!!    due to the restrictions in the i,j,k loops
!!
      implicit none
!
      class(cc3) :: wf
      integer, intent(in) :: i, j, k
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)           :: c_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: u_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(out)  :: sigma_abij
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)           :: g_bdci
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)           :: g_bdcj
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)           :: g_bdck
      real(dp), dimension(wf%n_o, wf%n_v), intent(in)                   :: g_ljci
      real(dp), dimension(wf%n_o, wf%n_v), intent(in)                   :: g_lkci
      real(dp), dimension(wf%n_o, wf%n_v), intent(in)                   :: g_lkcj
      real(dp), dimension(wf%n_o, wf%n_v), intent(in)                   :: g_licj
      real(dp), dimension(wf%n_o, wf%n_v), intent(in)                   :: g_lick
      real(dp), dimension(wf%n_o, wf%n_v), intent(in)                   :: g_ljck
!
   end subroutine jacobian_transpose_cc3_a_n7_cc3
!
!
   module subroutine construct_y_intermediates_cc3(wf, i, j, k, c_abc, u_abc, t_abij, y_cmjk,   &
                                                   Y_bcei, Y_bcej, Y_bcek)
!!
!!    Construct Y intermediates
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!!    Y_bcei = sum_aij c^abc_ijk * t^ae_ij
!!    Y_cmjk = sum_abj c^bac_ijk * t^ba_mj
!!
!!    All permutations for i,j,k have to be considered due to the restrictions in the i,j,k loops
!!
      implicit none
!
      class(cc3) :: wf
      integer, intent(in) :: i, j, k
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)           :: c_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: u_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in)   :: t_abij
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_o, wf%n_o), intent(out)  :: Y_cmjk
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: Y_bcei
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: Y_bcej
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: Y_bcek
!
   end subroutine construct_y_intermediates_cc3
!
!
   module subroutine jacobian_transpose_cc3_c3_a1_y_o_cc3(wf, sigma_ai, Y_cmjk)
!!
!!    Jacobian tranpose contribution of Y_vooo
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!!    sigma_1 += sum_mjk Y_cmjk * g_mjlk
!!    sigma_1 += sum_cmj Y_cmjk * g_mjcd
!!    sigma_1 += sum_cmk Y_cmjk * g_mdck
!!   
      implicit none
!
      class(cc3) :: wf
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_o, wf%n_o), intent(in) :: Y_cmjk
!
   end subroutine jacobian_transpose_cc3_c3_a1_y_o_cc3
!
!
   module subroutine jacobian_transpose_cc3_c3_a1_y_v_cc3(wf, sigma_ai)
!!
!!    Jacobian transpose contribution Y_vvvo
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!!    sigma_dk += sum_bec Y_bcek * g_becd
!!    sigma_cl += sum_bek Y_bcek * g_lkbe
!!    sigma_bl += sum_cek Y_bcek * g_leck
!!
      implicit none
!
      class(cc3) :: wf
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
!
   end subroutine jacobian_transpose_cc3_c3_a1_y_v_cc3
