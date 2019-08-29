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
   module subroutine prepare_for_density_cc2(wf)
!!
!!    Prepare for the construction of density matrices
!!    Written by Sarai D. Folekstad, May 2019
!!
      implicit none
!
      class(cc2), intent(inout) :: wf
!
   end subroutine prepare_for_density_cc2
!
!
   module subroutine construct_gs_density_cc2(wf)
!!
!!    Construct one-electron density
!!    Written by Sarai Dery Folkestad, 2019
!!
      implicit none
!
      class(cc2) :: wf
!
   end subroutine construct_gs_density_cc2
!
!
   module subroutine gs_one_el_density_cc2_oo_cc2(wf, density, tbar_akbj, t_akbi)
!!
!!    One electron density oo
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!
      class(cc2) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: density
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: tbar_akbj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: t_akbi
!
   end subroutine gs_one_el_density_cc2_oo_cc2
!
!
   module subroutine gs_one_el_density_cc2_vv_cc2(wf, density, tbar_ajci, t_bjci)
!!
!!    One electron density vv
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!
      class(cc2) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: density
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: tbar_ajci
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: t_bjci
!
   end subroutine gs_one_el_density_cc2_vv_cc2
!
!
   module subroutine gs_one_el_density_cc2_ov_cc2(wf, density, tbar_ai, t_aibj)
!!
!!    One electron density ov
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!
      class(cc2) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: density
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: tbar_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: t_aibj
!
   end subroutine gs_one_el_density_cc2_ov_cc2
!
!
