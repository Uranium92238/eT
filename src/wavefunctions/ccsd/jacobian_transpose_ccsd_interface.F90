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
   module subroutine prepare_for_jacobian_transpose_ccsd(wf)
!!
!!    Prepare for jacobian transpose
!!    Written by Tor S. Haugland, Oct 2019
!!
!!    Based on prepare_for_jacobian_ccsd by E. F. Kjønstad and S. D. Folkestad
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
!
   end subroutine prepare_for_jacobian_transpose_ccsd
!
!
   module subroutine jacobian_transpose_transformation_ccsd(wf, b)
!!
!!    Jacobian transpose transformation
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: b
!
   end subroutine jacobian_transpose_transformation_ccsd
!
!
   module subroutine save_jacobian_transpose_d1_intermediates_ccsd(wf, t_aibj)
!!
!!    Save Jacobian transpose D1 intermediates
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, and 
!!    Tor S. Haugland, Oct 2019
!!
      implicit none
!
      class(ccsd) :: wf
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: t_aibj
!
   end subroutine save_jacobian_transpose_d1_intermediates_ccsd
!
!
   module subroutine save_jacobian_transpose_e1_intermediates_ccsd(wf, t_aibj, L_ilmd)
!!
!!    Save Jacobian transpose E1 intermediates
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, and 
!!    Tor S. Haugland, Oct 2019
!!
      implicit none
!
      class(ccsd) :: wf
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: t_aibj
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o, wf%n_v), intent(in) :: L_ilmd
!
   end subroutine save_jacobian_transpose_e1_intermediates_ccsd
!
!
   module subroutine save_jacobian_transpose_f1_intermediates_ccsd(wf, t_aibj, g_ikmc)
!!
!!    Save Jacobian transpose F1 intermediates
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, and
!!    Andreas Skeidsvoll, Oct 2019
!!
      implicit none
!
      class(ccsd) :: wf
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(in) :: t_aibj
      real(dp), dimension(wf%n_o,wf%n_o,wf%n_o,wf%n_v), intent(in) :: g_ikmc
!
   end subroutine save_jacobian_transpose_f1_intermediates_ccsd
!
!
   module subroutine save_jacobian_transpose_g1_intermediates_ccsd(wf, t_aibj)
!!
!!    Save Jacobian transpose D2 intermediates
!!    Written by Andreas Skeidsvoll, Tor S. Haugland,
!!    Sarai D. Folkestad and Eirik F. Kjønstad , Nov 2019
!!
      implicit none
!
      class(ccsd) :: wf
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(in) :: t_aibj
!
   end subroutine save_jacobian_transpose_g1_intermediates_ccsd
!
!
   module subroutine save_jacobian_transpose_d2_intermediates_ccsd(wf, u_ckdl, L_dlbj)
!!
!!    Save Jacobian transpose D2 intermediates
!!    Written by Andreas Skeidsvoll, Tor S. Haugland,
!!    Sarai D. Folkestad and Eirik F. Kjønstad , Nov 2019
!!
      implicit none
!
      class(ccsd) :: wf
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: u_ckdl
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: L_dlbj ! L_jbld
!
   end subroutine save_jacobian_transpose_d2_intermediates_ccsd
!
!
   module subroutine save_jacobian_transpose_f2_intermediates_ccsd(wf, t_ckdl, L_dlbi)
!!
!!    Save Jacobian transpose f2 intermediates
!!    Written by Eirik F. Kjønstad, Tor S. Haugland and 
!!    Sarai D. Folkestad, Nov 2019
!!
      implicit none
!
      class(ccsd) :: wf
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: t_ckdl
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: L_dlbi
!
   end subroutine save_jacobian_transpose_f2_intermediates_ccsd
!
!
   module subroutine save_jacobian_transpose_g2_intermediates_ccsd(wf, t_aibj, g_kdib)
!!
!!    Save Jacobian transpose g2 intermediates
!!    Written by Tor S. Haugland, Eirik F. Kjønstad and 
!!    S. D. Folkestad, Nov 2019
!!
      implicit none
!
      class(ccsd) :: wf
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(in) :: t_aibj
      real(dp), dimension(wf%n_o,wf%n_v,wf%n_o,wf%n_v), intent(in) :: g_kdib
!
   end subroutine save_jacobian_transpose_g2_intermediates_ccsd
!
!
   module subroutine save_jacobian_transpose_i2_intermediates_ccsd(wf, t_aibj, g_ovov)
!!
!!    Save Jacobian transpose i2 intermediates
!!    Written by Tor S. Haugland, Eirik F. Kjønstad and 
!!    Sarai D. Folkestad, Nov 2019
!!
      implicit none
!
      class(ccsd) :: wf
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(in) :: t_aibj
      real(dp), dimension(wf%n_o,wf%n_v,wf%n_o,wf%n_v), intent(in) :: g_ovov
!
   end subroutine save_jacobian_transpose_i2_intermediates_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_d1_ccsd(wf, sigma_ai, b_aibj)
!!
!!    Jacobian transpose CCSD D1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o)                 :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
   end subroutine jacobian_transpose_ccsd_d1_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_e1_ccsd(wf, sigma_ai, b_aibj)
!!
!!    Jacobian transpose CCSD E1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o)                 :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
   end subroutine jacobian_transpose_ccsd_e1_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_f1_ccsd(wf, sigma_ai, b_aibj)
!!
!!    Jacobian transpose CCSD F1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o)                 :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
   end subroutine jacobian_transpose_ccsd_f1_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_g1_ccsd(wf, sigma_ai, b_aibj)
!!
!!    Jacobian transpose CCSD G1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o)                 :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
   end subroutine jacobian_transpose_ccsd_g1_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_b2_ccsd(wf, sigma_aibj, b_aibj)
!!
!!    Jacobian transpose CCSD B2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: sigma_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
   end subroutine jacobian_transpose_ccsd_b2_ccsd
!

!
   module subroutine jacobian_transpose_ccsd_c2_ccsd(wf, sigma_aibj, b_aibj)
!!
!!    Jacobian transpose CCSD C2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: sigma_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
   end subroutine jacobian_transpose_ccsd_c2_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_d2_ccsd(wf, sigma_aibj, b_aibj)
!!
!!    Jacobian transpose CCSD D2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
   implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: sigma_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
   end subroutine jacobian_transpose_ccsd_d2_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_e2_ccsd(wf, sigma_aibj, b_aibj)
!!
!!    Jacobian transpose CCSD E2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: sigma_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
   end subroutine jacobian_transpose_ccsd_e2_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_f2_ccsd(wf, sigma_aibj, b_aibj)
!!
!!    Jacobian transpose CCSD F2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: sigma_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
   end subroutine jacobian_transpose_ccsd_f2_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_g2_ccsd(wf, sigma_aibj, b_aibj)
!!
!!    Jacobian transpose CCSD G2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: sigma_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: b_aibj
!
   end subroutine jacobian_transpose_ccsd_g2_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_h2_ccsd(wf, sigma_abij, b_abij)
!!
!!    Jacobian transpose CCSD H2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o) :: sigma_abij
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o) :: b_abij
!
   end subroutine jacobian_transpose_ccsd_h2_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_i2_ccsd(wf, sigma_abij, b_abij)
!!
!!    Jacobian transpose CCSD I2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o) :: sigma_abij
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o) :: b_abij
!
   end subroutine jacobian_transpose_ccsd_i2_ccsd
!
!
   module subroutine save_jacobian_transpose_e2_oo_intermediate_ccsd(wf, t_ckdl, L_ckdj)
!!
!!    Save Jacobian transpose e2 oo intermediate
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018    
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: t_ckdl
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: L_ckdj ! L_kcjd
!
   end subroutine save_jacobian_transpose_e2_oo_intermediate_ccsd
!
!
   module subroutine save_jacobian_transpose_e2_vv_intermediate_ccsd(wf, t_ckdl, L_bkdl)
!!
!!    Save Jacobian transpose e2 vv intermediate
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018    
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: t_ckdl
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: L_bkdl ! L_kbld
!
   end subroutine save_jacobian_transpose_e2_vv_intermediate_ccsd
