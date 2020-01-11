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
   module subroutine prepare_for_density_cc3(wf)
!!
!!    Prepare for ground state density
!!    Written by Alexander C. Paul and Rolf H. Myhre, March 2019
!!
!!    Sets up the integral files needed to construct t3 and tbar3
!!    batching over the virtual indices.
!!    g_vvvo ordered 2413
!!    g_oovo ordered 1243
!!    g_vvov ordered 1324
!!    g_ooov ordered 2134
!!    L_ovov ordered 1324
!!
      implicit none
!
      class(cc3), intent(inout) :: wf
!
   end subroutine prepare_for_density_cc3
!
!
   module subroutine construct_gs_density_cc3(wf)
!!
!!    Construct one-electron density
!!    Written by Alexander C. Paul
!!    based on construct_gs_density_ccsd by Sarai D. Folkestad
!!
!!    Constructs the one-electron density matrix in the T1 basis
!!
!!    D_pq = < Lambda| E_pq |CC >
!!
!!    Contributions to the density are split up as follows:
!!    D_pq = D_pq(ref-ref) + sum_mu tbar_mu D_pq(mu-ref)
!!
      implicit none
!
      class(cc3) :: wf
!
   end subroutine construct_gs_density_cc3
!
!
   module subroutine density_cc3_mu_ref_abc_cc3(wf, density_oo, omega, tbar_ia, &
                                                tbar_ijab, t_ijab, cvs)
!!
!!    One electron density excited-determinant/reference term 
!!    in batches of the virtual orbitals a,b,c
!!    Written by Alexander C. Paul, July 2019
!!
!!    Computes terms of the form:
!!
!!          D_pq += sum_mu X_mu * < mu| e^(-T) E_pq e^T |HF >
!!
!!    where X_mu is a general amplitude (tbar or L)
!!
!!    Construct t^abc_ijk and tbar^abc_ijk in batches of a,b,c.
!!    Calls routines calculating the actual density contributions
!!    to the oo-block
!!
!!    t_mu3 = -< mu3|{U,T2}|HF > (eps_mu3)^-1
!!    tbar_mu3 = (- eps_mu3)^-1 (tbar_mu1 < mu1| [H,tau_nu3] |R > 
!!                             + tbar_mu2 < mu2| [H,tau_nu3] |R >
!!
!!    oo-block:
!!          D_kl += -1/2 sum_ij,abc t^abc_ijk tbar^abc_ijl    
!!
      implicit none
!
      class(cc3) :: wf
      real(dp), dimension(wf%n_o, wf%n_o), intent(inout) :: density_oo
      real(dp), intent(in) :: omega
      logical, intent(in)  :: cvs
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: tbar_ia
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_v, wf%n_v), intent(in) :: tbar_ijab
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_v, wf%n_v), intent(in) :: t_ijab
!
   end subroutine density_cc3_mu_ref_abc_cc3
!
!
   module subroutine density_cc3_mu_ref_oo_cc3(wf, a, b, c, density_oo, t_ijk, &
                                                   u_ijk, tbar_ijk, v_ijk)
!!
!!    One electron density excited-determinant/reference oo-term 
!!    Written by Alexander C. Paul, August 2019
!!
!!    Calculates CC3 terms of the form:
!!
!!          D_pq += sum_mu X_mu * < mu| e^(-T) E_pq e^T |HF >
!!
!!    where X_mu is a general amplitude (tbar or L)
!!
!!    explicit Term:
!!          D_kl += -1/2 sum_ij,abc t^abc_ijk tbar^abc_ijl
!!
!!    All permutations for a,b,c have to be considered 
!!    due to the restrictions in the a,b,c loops   
!!
      implicit none
!
      class(cc3) :: wf
      integer, intent(in) :: a, b, c
      real(dp), dimension(wf%n_o, wf%n_o), intent(inout)  :: density_oo
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(in)  :: t_ijk
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(out) :: u_ijk
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(in)  :: tbar_ijk
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(out) :: v_ijk
!
   end subroutine density_cc3_mu_ref_oo_cc3
!
!
   module subroutine density_cc3_mu_ref_ijk_cc3(wf, density_ov, density_vv, omega, &
                                                tbar_ai, tbar_abij, t_abij,        &
                                                cvs, keep_Y)
!!
!!    One electron density excited-determinant/reference term 
!!    in batches of the occupied orbitals i,j,k
!!    Written by Alexander C. Paul, July 2019
!!
!!    Computes terms of the form:
!!
!!          D_pq += sum_mu X_mu * < mu| e^(-T) E_pq e^T |HF >
!!
!!    where X_mu is a general amplitude (tbar or L)
!!
!!    Construct t^abc_ijk and tbar^abc_ijk in batches of a,b,c.
!!    Calls routines calculating the actual density contributions
!!    to the ov- and vv-blocks.
!!
!!    t_mu3 = -< mu3| [U,T2] |HF > (eps_mu3)^-1
!!    tbar_mu3 = (- eps_mu3)^-1 (tbar_mu1 < mu1 | [H,tau_nu3] |R > 
!!                             + tbar_mu2 < mu2 | [H,tau_nu3] |R >
!!
!!    ov-block:
!!       rho^L_kc += sum_{abij} tbar^ab_ij (t^abc_ijk - t^bac_ijk)
!!       rho^L_ld -= sum_{abijk} tbar^abc_ijk t^ab_lj t^dc_ik
!!
!!    vv-block:
!!       rho^L_cd += 1/2 sum_{abijk} tbar^abc_ijk t^abd_ijk
!!
      implicit none
!
      class(cc3) :: wf
      real(dp), intent(in) :: omega
      logical, intent(in)  :: cvs
      logical, intent(in), optional :: keep_Y
      real(dp), dimension(wf%n_o, wf%n_v), intent(inout) :: density_ov
      real(dp), dimension(wf%n_v, wf%n_v), intent(inout) :: density_vv
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: tbar_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: tbar_abij
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: t_abij
!
   end subroutine density_cc3_mu_ref_ijk_cc3
!
!
   module subroutine density_cc3_mu_ref_vv_cc3(wf, i, j, k, density_vv, t_abc, &
                                               u_abc, tbar_abc, v_abc)
!!
!!    One electron density excited-determinant/reference vv-term 
!!    Written by Alexander C. Paul, August 2019
!!
!!    Calculates CC3 terms of the form:
!!
!!          D_pq += sum_mu X_mu * < mu| e^(-T) E_pq e^T |HF >
!!
!!    where X_mu is a general amplitude (tbar or L)
!!
!!    explicit Term:
!!          D_cd += 1/2 sum_ab,ijk tbar^abc_ijk t^abd_ijk
!!
!!    All permutations for i,j,k have to be considered 
!!    due to the restrictions in the i,j,k loops   
!!
      implicit none
!
      class(cc3) :: wf
      integer, intent(in) ::  i, j, k
      real(dp), dimension(wf%n_v, wf%n_v), intent(inout)  :: density_vv
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)  :: t_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: u_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)  :: tbar_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: v_abc
!
   end subroutine density_cc3_mu_ref_vv_cc3
!
!
   module subroutine construct_y_intermediate_vo3_cc3(wf, i, j, k, tbar_abc, &
                                                      u_abc, t_abij, Y_clik)
!!
!!    Construct Y_vooo intermediate  
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!!    Y_clik = sum_abj tbar^abc_ijk * t^ab_lj
!!    used to compute the tbar3 contributions to the ov part of the GS density
!!
!!    All permutations for i,j,k have to be considered 
!!    due to the restrictions in the i,j,k loops
!!
!!    equal to some of the terms in jacobian_transpose construct_y_intermediates
!!
      implicit none
!
      class(cc3) :: wf
      integer, intent(in) :: i, j, k
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)           :: tbar_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: u_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in)   :: t_abij
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_o, wf%n_o), intent(out)  :: Y_clik
!
   end subroutine construct_y_intermediate_vo3_cc3
