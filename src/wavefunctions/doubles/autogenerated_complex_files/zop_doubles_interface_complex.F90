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
   module subroutine construct_gs_density_doubles_complex(wf)
!!
!!    Construct one_complex-electron density
!!    Written by Sarai Dery Folkestad, 2019
!!
      implicit none
!
      class(doubles) :: wf
!
      complex(dp), dimension(:,:,:,:), allocatable :: tbar_aibj
      complex(dp), dimension(:,:,:,:), allocatable :: t_aibj
!
   end subroutine construct_gs_density_doubles_complex
!
!
   module subroutine gs_one_el_density_doubles_oo_doubles_complex(wf, density, tbar_akbj, t_akbi)
!!
!!    One electron density oo
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!
      class(doubles) :: wf
!
      complex(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: density
!
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: tbar_akbj
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: t_akbi
!
   end subroutine gs_one_el_density_doubles_oo_doubles_complex
!
!
   module subroutine gs_one_el_density_doubles_vv_doubles_complex(wf, density, tbar_ajci, t_bjci)
!!
!!    One electron density vv
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!
      class(doubles) :: wf
!
      complex(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: density
!
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: tbar_ajci
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: t_bjci
!
   end subroutine gs_one_el_density_doubles_vv_doubles_complex
!
!
   module subroutine gs_one_el_density_doubles_ov_doubles_complex(wf, density, tbar_ai, t_aibj)
!!
!!    One electron density ov
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!
      class(doubles) :: wf
!
      complex(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: density
!
      complex(dp), dimension(wf%n_v, wf%n_o), intent(in) :: tbar_ai
      complex(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: t_aibj
!
   end subroutine gs_one_el_density_doubles_ov_doubles_complex
!
