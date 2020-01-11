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
   module subroutine construct_left_transition_density_cc3(wf, state)
!!
!!    Construct left one-electron transition density
!!    Written by Alexander C. Paul, June 2019
!!
!!          D^L_pq = < k| E_pq |CC >
!!
!!    where < k| is the left eigenvector of the Jacobian
!!    with amplitudes L_mu
!!
!!          < k| = sum_mu L_{k,mu} < mu| e^-T
!!
      implicit none
!
      class(cc3) :: wf
      integer, intent(in) :: state
!
   end subroutine construct_left_transition_density_cc3
!
!
   module subroutine construct_right_transition_density_cc3(wf, state)
!!
!!    Construct right one-electron transition density
!!    Written by Alexander C. Paul, June 2019
!!
!!          rho^R_pq = < Lambda| E_pq |k >
!!
!!    where |k > is the right eigenvector of the Jacobian
!!    with amplitudes R_mu
!!
!!          |k > = sum_mu (tau_mu |CC > R_{k,mu} - tbar_mu |CC > R_{k,mu}) 
!!
      implicit none
!
      class(cc3) :: wf
      integer, intent(in) :: state
!
   end subroutine construct_right_transition_density_cc3
!
!
   module subroutine density_cc3_mu_nu_ov_cc3(wf, density_ov, R_ai)
!!
!!    One electron density (EOM) excited-determinant/excited-determinant ov-term 
!!    Written by Alexander C. Paul, August 2019
!!
!!    Computes terms of the form:
!!
!!          D_pq += sum_mu,nu X_mu < mu| e^(-T) E_pq e^T |nu > Y_nu
!!
!!    explicit term in this routine:
!!
!!       D^R_ld = -1/2 sum_{abcijk} tbar^abc_ijk(R^c_l t^abd_ijk + R^d_k t^abc_ijl)
!!                = sum_{abcijk}( -1/2 tbar^abc_ijk t^abc_ijl R^d_k 
!!                                -1/2 tbar^abc_ijk t^abd_ijk R^c_l)
!!                = sum_k D_lk R^d_k - sum_c D_cd R^c_l
!!
!!    Calculates the contribution of the oo- and vv-blocks
!!    of the GS density to the ov-block of the right TDM
!!
      implicit none
!
      class(cc3) :: wf
      real(dp), dimension(wf%n_o, wf%n_v), intent(inout) :: density_ov
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: R_ai
!
   end subroutine density_cc3_mu_nu_ov_cc3
!
!
   module subroutine density_cc3_mu_nu_oo_ov_vv_cc3(wf, density_oo, density_ov, &
                                                    density_vv, R_ai, R_aibj)
!!
!!    One electron density (EOM) excited-determinant/excited-determinant 
!!    oo-, ov- and vv-term 
!!    Written by Alexander C. Paul, August 2019
!!
!!    Computes terms of the form:
!!
!!          D_pq += sum_mu,nu X_mu < mu| e^(-T) E_pq e^T |nu > Y_nu
!!
!!    explicit terms in this routine:
!!
!!          D^R_lk -= sum{abcij}  tbar^abc_ijk t^bc_jl R^a_i
!!                 -= sum{ai}  Y_alki R^a_i
!!
!!          D^R_ld -= sum{abcijk} tbar^abc_ijk (t^bd_jk R^ac_il + t^bc_jl R^ad_ik)
!!                 -= sum{aci} Y_cadi R^ac_il + sum{aik} Y_alki R^ad_ik
!!
!!          D^R_cd += sum{abijk}  tbar^abc_ijk t^bd_jk R^a_i
!!                 += sum{ai}  Y_cadi R^a_i
!!
!!    Calculates the contribution of the Y_clik- and Y_bcdk-intermediates
!!    to the oo-, ov- and vv-block of the right TDM
!!
!!    The intermediates are constructed and written to file during
!!    the construction of the GS-density and while solving for the multipliers
!!
!!          Y_bcek = sum_aij tbar^abc_ijk * t^ae_ij
!!          Y_clik = sum_abj tbar^abc_ijk * t^ab_lj
!!
      implicit none
!
      class(cc3) :: wf
      real(dp), dimension(wf%n_o, wf%n_o), intent(inout) :: density_oo
      real(dp), dimension(wf%n_o, wf%n_v), intent(inout) :: density_ov
      real(dp), dimension(wf%n_v, wf%n_v), intent(inout) :: density_vv
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: R_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: R_aibj
!
   end subroutine density_cc3_mu_nu_oo_ov_vv_cc3
!
!
   module subroutine density_cc3_mu_nu_ijk_cc3(wf, density_ov, density_vo, &
                                               density_vv, omega,          &
                                               tbar_ai, tbar_abij,         &
                                               R_ai, R_abij, tbar_R_overlap)
!!
!!    One electron density excited-determinant/excited-determinant term 
!!    in batches of the occupied orbitals i,j,k
!!    Written by Alexander C. Paul, July 2019
!!
!!    Computes terms of the form:
!!
!!          D_pq += sum_mu,nu X_mu * < mu| e^(-T) E_pq e^T |nu > Y_nu
!!
!!    where X_mu and Y_nu are general amplitude (tbar or L)
!!
!!
!!    Construct R^abc_ijk and tbar^abc_ijk in batches of i,j,k and compute
!!    the contributions to the vo- and vv-block of the right TDM.
!!
!!    explicit terms in this routine
!!
!!    R_mu3 = (omega - eps_mu3)^-1 (< mu3| [H,R_2] |HF > 
!!                                + < mu3| [[H,R_1],T_2] |HF >)
!!    tbar_mu3 = (- eps_mu3)^-1 (tbar_mu1 < mu1| [H,tau_nu3] |R > 
!!                             + tbar_mu2 < mu2| [H,tau_nu3] |R >
!!
!!    vo-part:
!!          D^R_ck += 1/2 sum{abij} tbar^abc_ijk R^ab_ij
!!
!!    vv-part:
!!          D^R_cd += 1/2 sum_{abijk} tbar^abc_ijk R^abd_ijk
!!
!!    Also construct the intermediate Z_bcjk needed for the ov-block
!!          Z_bcjk = sum{ai} tbar^abc_ijk R^a_i
!!
      implicit none
!
      class(cc3) :: wf
      real(dp), intent(in) :: omega
      real(dp), intent(inout) :: tbar_R_overlap
      real(dp), dimension(wf%n_o, wf%n_v), intent(inout) :: density_ov
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: density_vo
      real(dp), dimension(wf%n_v, wf%n_v), intent(inout) :: density_vv
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: tbar_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: tbar_abij
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: R_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: R_abij
!
   end subroutine density_cc3_mu_nu_ijk_cc3
!
!
   module subroutine density_cc3_mu_nu_vo_cc3(wf, i, j, k, tbar_abc,   &
                                              v_abc, density_vo, R_abij)
!!
!!    One electron density (EOM) excited-determinant/excited-determinant vo-term 
!!    Written by Alexander C. Paul, August 2019
!!
!!    Computes terms of the form:
!!
!!          D_pq += sum_mu,nu X_mu < mu| e^(-T) E_pq e^T |nu > Y_nu
!!
!!    explicit term in this routine:
!!
!!          D^R_ck += sum_{abij} tbar^abc_ijk R^ab_ij
!!
!!    All permutations for i,j,k have to be considered 
!!    due to the restrictions in the i,j,k loops
!!
!!    based on construct_x_ai_intermediate_cc3
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      class(cc3) :: wf
      integer, intent(in) :: i, j, k
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)         :: tbar_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)        :: v_abc
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)              :: density_vo
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: R_abij
!
   end subroutine density_cc3_mu_nu_vo_cc3
!
!
   module subroutine construct_Z_intermediate_cc3(wf, i, j, k, tbar_abc, &
                                                  v_abc, Z_bcjk, R_ai)
!!
!!    Constructs Z-intermediate 
!!    Written by Alexander C. Paul, August 2019
!!
!!    based on jacobian_cc3_b2_fock_cc3
!!    Written by Alexander C. Paul and Rolf H. Myhre, Feb 2019
!!
!!    needed for the ov-block of the right TDM
!!
!!       Z_bcjk += sum_{ai} tbar^abc_ijk R^a_i
!!
!!    All permutations for i,j,k have to be considered 
!!    due to the restrictions in the i,j,k loops
!!      
      implicit none
!
      class(cc3) :: wf
      integer, intent(in) :: i, j, k
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)           :: tbar_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: v_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(out)  :: Z_bcjk
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)   :: R_ai
!
   end subroutine construct_Z_intermediate_cc3
!
!
   module subroutine density_cc3_mu3_nu2_ov_cc3(wf, density_ov, density_vo)
!!
!!    One electron density (EOM) triples-determinant/doubles-determinant ov-term 
!!    Written by Alexander C. Paul, August 2019
!!
!!    Computes terms of the form:
!!
!!          D_pq += sum_mu,nu X_mu < mu| e^(-T) E_pq e^T |nu > Y_nu
!!
!!    explicit term in this routine:
!!
!!          D^R_ld += sum{abcijk} 1/2 tbar^abc_ijk R^ab_ij(2t^cd_kl-t^cd_lk)
!!                 += sum{ck} D^R_ck (2t^cd_kl-t^cd_lk)
!!
!!    Calculates the contribution of the vo-block of the right TDM to
!!    the ov-block of the right TDM
!!
      implicit none
!
      class(cc3) :: wf
      real(dp), dimension(wf%n_o, wf%n_v), intent(inout) :: density_ov
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: density_vo
!
   end subroutine density_cc3_mu3_nu2_ov_cc3
!
!
   module subroutine density_cc3_mu_nu_abc_cc3(wf, density_oo, omega, &
                                               tbar_ia, tbar_ijab, R_ijab)
!!
!!    One electron density excited-determinant/excited-determinant term 
!!    in batches of the virtual orbitals a,b,c
!!    Written by Alexander C. Paul, August 2019
!!
!!    Computes terms of the form:
!!
!!          D_pq += sum_mu,nu X_mu * < mu| e^(-T) E_pq e^T |nu > Y_nu
!!
!!    where X_mu and Y_nu are general amplitude (tbar or L)
!!
!!    Construct R^abc_ijk and tbar^abc_ijk in batches of a,b,c
!!    and compute the contribution to the oo-block of the right TDM.
!!
!!    explicit terms in this routine
!!
!!    R_mu3 = (omega - eps_mu3)^-1 (< mu3| [H,R_2] |HF > 
!!                                + < mu3| [[H,R_1],T_2] |HF >)
!!
!!    tbar_mu3 = (- eps_mu3)^-1 (tbar_mu1 < mu1| [H,tau_nu3] |R > 
!!                             + tbar_mu2 < mu2| [H,tau_nu3] |R >
!!
!!    oo-part:
!!          D^R_kl -= 1/2 sum_{abcij} tbar^abc_ijl R^abc_ijk
!!
      implicit none
!
      class(cc3) :: wf
      real(dp), intent(in) :: omega
      real(dp), dimension(wf%n_o, wf%n_o), intent(inout) :: density_oo
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: tbar_ia
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: tbar_ijab
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: R_ijab
!
   end subroutine density_cc3_mu_nu_abc_cc3
!
!
   module subroutine density_cc3_mu_nu_ov_Z_term_cc3(wf, density_ov)
!!
!!    One electron density (EOM) excited-determinant/excited-determinant 
!!    ov-term from Z_intermediate
!!    Written by Alexander C. Paul, August 2019
!!
!!    Computes terms of the form:
!!
!!          D_pq += sum_mu,nu X_mu < mu| e^(-T) E_pq e^T |nu > Y_nu
!!
!!    explicit term in this routine:
!!
!!       D^R_ld += sum_{abcijk} tbar^abc_ijk R^a_i (t^bcd_jkl - t^bcd_jlk)
!!              += sum_{bcjk} Z_bcjk (t^bcd_jkl - t^bcd_jlk)
!!
!!    Construct t3-amplitudes in batches of i,j,k and calculate the
!!    contribution to the ov-block involving the Z-intermediate
!!
!!       t_mu3 = -< mu3| [U,T2] |HF > (eps_mu3)^-1
!!
!!    The Z-intermediate was constructed in density_cc3_mu_nu_ijk_cc3
!!
      implicit none
!
      class(cc3) :: wf
      real(dp), dimension(wf%n_o, wf%n_v), intent(inout) :: density_ov
!
   end subroutine density_cc3_mu_nu_ov_Z_term_cc3
