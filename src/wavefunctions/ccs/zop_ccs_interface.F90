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
   module subroutine prepare_for_density_ccs(wf)
!!
!!    Prepare for the construction of density matrices
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2019
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
   end subroutine prepare_for_density_ccs
!
!
   module subroutine construct_gs_density_ccs(wf)
!!
!!    Construct one-electron density
!!    Written by Sarai Dery Folkestad, 2019
!!
      implicit none
!
      class(ccs) :: wf
!
   end subroutine construct_gs_density_ccs
!
!
   module subroutine gs_one_el_density_ccs_oo_ccs(wf, density)
!!
!!    One electron density oo
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: density
!
   end subroutine gs_one_el_density_ccs_oo_ccs
!
!
   module subroutine gs_one_el_density_ccs_vo_ccs(wf, density, tbar_ai)
!!
!!    One electron density vo
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: density
      real(dp), dimension(wf%n_v, wf%n_o) :: tbar_ai
!
   end subroutine gs_one_el_density_ccs_vo_ccs
!
!
   module function calculate_expectation_value_ccs(wf, A, density) result(expectation_value)
!!
!!    Calculate expectation value
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!  
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: A
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: density
!
      real(dp) :: expectation_value
!
   end function calculate_expectation_value_ccs
!
!
   module subroutine calculate_energy_ccs(wf)
!!
!!    Calculate energy
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
   end subroutine calculate_energy_ccs