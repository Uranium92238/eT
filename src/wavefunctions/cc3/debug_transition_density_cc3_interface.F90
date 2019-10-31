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
   module subroutine left_TDM_debug_cc3(wf, state, l_TDM)
!!
!!    Construct CC3 left transition density matrix
!!    Written by Alexander C. Paul, Sep 2019
!!
!!    NB: only for debugging
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: state
!
      real(dp), dimension(wf%n_mo,wf%n_mo), intent(out) :: l_TDM
!
   end subroutine left_TDM_debug_cc3
!
!
   module subroutine right_TDM_debug_cc3(wf, state, r_TDM)
!!
!!    Construct CC3 right transition density matrix
!!    Written by Alexander C. Paul, Sep 2019
!!
!!    NB: only for debugging
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: state
!
      real(dp), dimension(wf%n_mo,wf%n_mo), intent(inout) :: r_TDM
!
   end subroutine right_TDM_debug_cc3
!
!
   module subroutine construct_full_R3_cc3(wf, R_abij, R_abcijk, omega)
!!
!!    Construct full R3
!!    Written by Alexander C. Paul, Sep 2019
!!
!!    constructs full v3o3 array containing the triples excitation vector 
!!    using unwrapped loops
!!
!!     NB: the amplitudes are scaled by the biorthonormal factor
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_o,wf%n_o), intent(in) :: R_abij
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_v,wf%n_o,wf%n_o,wf%n_o), intent(inout) :: R_abcijk
!
      real(dp), intent(in) :: omega
!
   end subroutine construct_full_R3_cc3
!
!
   module subroutine construct_full_t3_cc3(wf, t_abcijk)
!!
!!    Construct full t3
!!    Written by Alexander C. Paul, Sep 2019
!!
!!    constructs full v3o3 array containing 
!!    the triples ground state amplitudes
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_v,wf%n_o,wf%n_o,wf%n_o), intent(inout) :: t_abcijk
!
   end subroutine construct_full_t3_cc3
!
!
   module subroutine construct_full_tbar3_cc3(wf, tbar_ai, tbar_abij, tbar_abcijk, omega, cvs)
!!
!!    Construct full tbar3 debug
!!    Written by Alexander C. Paul, Sep 2019
!!
!!    constructs full v3o3 array containing the triples multipliers
!!    using unwrapped loops
!!
!!    NB: can also be used for left excitation vector
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: tbar_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: tbar_abij
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v, wf%n_o, wf%n_o, wf%n_o), intent(inout) :: tbar_abcijk
!
      real(dp), intent(in) :: omega
!
      logical, intent(in) :: cvs
!
   end subroutine construct_full_tbar3_cc3
!
!
   module subroutine debug_left_oo_cc3(wf, density_oo, L_abcijk, t_abcijk)
!!
!!    CC3 contribution to the oo-block of the left TDM
!!    Written by Alexander C. Paul, September 2019
!!
!!    rho_kl -= 1/2 sum_abcijk L^abc_ijk*t^abc_ijk
!! 
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_o), intent(inout) :: density_oo
!
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_v,wf%n_o,wf%n_o,wf%n_o), intent(in) :: L_abcijk
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_v,wf%n_o,wf%n_o,wf%n_o), intent(in) :: t_abcijk
!
   end subroutine debug_left_oo_cc3
!
!
   module subroutine debug_left_ov_N7_cc3(wf, density_ov, L_abcijk)
!!
!!    CC3 contribution to the ov-block of the left TDM (N7 scaling)
!!    Written by Alexander C. Paul, September 2019
!!
!!    rho_ld -= sum_abcijk L^abc_ijk*t^ad_ik*t^bc_jl
!!    rho_ld -= sum_aik Y_alik*t^ad_ik
!!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_v), intent(inout) :: density_ov
!
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_v,wf%n_o,wf%n_o,wf%n_o), intent(in) :: L_abcijk
!
   end subroutine debug_left_ov_N7_cc3
!
!
   module subroutine debug_left_ov_N6_cc3(wf, density_ov, L_abij, t_abcijk)
!!
!!    CC3 contribution to the ov-block of the left TDM N6 scaling
!!    Written by Alexander C. Paul, September 2019
!!
!!    rho_ld += sum_abcijk L^ab_ij*(t^abd_ijl-t^adb_ijl)
!!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_v), intent(inout) :: density_ov
!
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_o,wf%n_o), intent(in) :: L_abij
!
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_v,wf%n_o,wf%n_o,wf%n_o), intent(in) :: t_abcijk
!
   end subroutine debug_left_ov_N6_cc3
!
!
   module subroutine debug_left_vv_cc3(wf, density_vv, L_abcijk, t_abcijk)
!!
!!    CC3 contribution to the vv-block of the left TDM
!!    Written by Alexander C. Paul, September 2019
!!
!!    rho_cd += 1/2 sum_abcijk L^abc_ijk*t^abd_ijk
!!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_v), intent(inout) :: density_vv
!
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_v,wf%n_o,wf%n_o,wf%n_o), intent(in) :: L_abcijk
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_v,wf%n_o,wf%n_o,wf%n_o), intent(in) :: t_abcijk
!
   end subroutine debug_left_vv_cc3
!
!
   module subroutine debug_right_vo_cc3(wf, density_vo, R_abij, tbar_abcijk)
!!
!!    CC3 contribution to the vo-block of the right TDM
!!    Written by Alexander C. Paul, September 2019
!!
!!    rho_ck += sum_abij 1/2 tbar^abc_ijk*R^ab_ij
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: density_vo
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: R_abij
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v, wf%n_o, wf%n_o, wf%n_o), intent(in) :: tbar_abcijk
!
   end subroutine debug_right_vo_cc3
!
!
   module subroutine debug_right_ov_R3_cc3(wf, density_ov, R_abcijk)
!!
!!    CC3 contribution to the ov-block of the right TDM (R3 contribution)
!!    Written by Alexander C. Paul, September 2019
!!
!!    rho_ld += sum_abij tbar^ab_ij*(R^abd_ijl - R^abd_ilj)
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_v), intent(inout) :: density_ov
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v, wf%n_o, wf%n_o, wf%n_o), intent(in) :: R_abcijk
!
   end subroutine debug_right_ov_R3_cc3
!
!
   module subroutine debug_right_ov_t3_cc3(wf, density_ov, density_vo, R_ai, tbar_abcijk, t_abcijk)
!!
!!    CC3 contribution to the ov-block of the right TDM (t3 contribution)
!!    Written by Alexander C. Paul, September 2019
!!
!!    rho_ld += -1/2 tbar^abc_ijk*t^abc_ijl*R^d_k
!!
!!    rho_ld += sum_abcijk 1/2 tbar^abc_ijk R^ab_ij (2t^cd_kl - t^cd_lk)
!!    rho_ld += sum_ck rho_ck (2t^cd_kl - t^cd_lk)
!!
!!    rho_ld += sum_abcijk tbar^abc_ijk*R^a_i*(t^bcd_jkl - t^bcd_jlk)
!!    rho_ld += sum_bcjk Z^bc_jk*(t^bcd_jkl - t^bcd_jlk)
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_v), intent(inout) :: density_ov
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: density_vo
!
      real(dp), dimension(wf%n_v, wf%n_o),intent(in) :: R_ai
!
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_v,wf%n_o,wf%n_o,wf%n_o), intent(in) :: tbar_abcijk, t_abcijk
!
   end subroutine debug_right_ov_t3_cc3
!
!
   module subroutine debug_right_ov_Y_contribution_cc3(wf, density_ov, R_abij, tbar_abcijk)
!!
!!    CC3 contribution to the ov-block of the right TDM (Y intermediates)
!!    Written by Alexander C. Paul, September 2019
!!
!!    rho_ld -= sum_abcijk tbar^abc_ijk*R^ad_ik*t^bc_jl
!!    rho_ld -= sum_aik Y_alik*R^ad_ik
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_v), intent(inout) :: density_ov
!
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_o,wf%n_o), intent(in) :: R_abij
!
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_v,wf%n_o,wf%n_o,wf%n_o), intent(in) :: tbar_abcijk
!
   end subroutine debug_right_ov_Y_contribution_cc3
!
!
   module subroutine debug_right_oo_cc3(wf, density_oo, R_ai, tbar_abcijk, R_abcijk)
!!
!!    CC3 contribution to the oo-block of the right TDM
!!    Written by Alexander C. Paul, September 2019
!!
!!    rho_lk -= sum_abcij (tbar^abc_ijk*R^a_i*t^bc_jl 
!!                      + half tbar^abc_ijk*R^abc_ijl*delta_ai,bj,cl)
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_o), intent(inout) :: density_oo
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: R_ai
!
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_v,wf%n_o,wf%n_o,wf%n_o), intent(in) :: tbar_abcijk, R_abcijk
!
   end subroutine debug_right_oo_cc3
!
!
   module subroutine debug_right_vv_cc3(wf, density_vv, R_ai, tbar_abcijk, R_abcijk)
!!
!!    CC3 contribution to the vv-block of the right TDM
!!    Written by Alexander C. Paul, September 2019
!!
!!    rho_cd += sum_abijk tbar^abc_ijk*R^a_i*t^bd_jk 
!!    rho_cd += half sum_abcijk tbar^abc_ijk*R^abd_ijk*Δ_ai,bj,dk
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_o), intent(inout) :: density_vv
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: R_ai
!
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_v,wf%n_o,wf%n_o,wf%n_o), intent(in) :: tbar_abcijk, R_abcijk
!
   end subroutine debug_right_vv_cc3