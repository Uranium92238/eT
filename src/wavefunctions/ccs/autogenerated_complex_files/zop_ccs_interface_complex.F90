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
   module subroutine prepare_for_density_ccs_complex(wf)
!!
!!    Prepare for the construction of density matrices
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2019
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
   end subroutine prepare_for_density_ccs_complex
!
!
   module subroutine construct_gs_density_ccs_complex(wf)
!!
!!    Construct one_complex-electron density
!!    Written by Sarai Dery Folkestad, 2019
!!
      implicit none
!
      class(ccs) :: wf
!
   end subroutine construct_gs_density_ccs_complex
!
!
   module subroutine gs_one_el_density_ccs_oo_ccs_complex(wf, density)
!!
!!    One electron density oo
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!
      class(ccs) :: wf
!
      complex(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: density
!
   end subroutine gs_one_el_density_ccs_oo_ccs_complex
!
!
   module subroutine gs_one_el_density_ccs_vo_ccs_complex(wf, density, tbar_ai)
!!
!!    One electron density vo
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!
      class(ccs) :: wf
!
      complex(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: density
      complex(dp), dimension(wf%n_v, wf%n_o) :: tbar_ai
!
   end subroutine gs_one_el_density_ccs_vo_ccs_complex
!
!
   module function calculate_expectation_value_ccs_complex(wf, A, density) result(expectation_value)
!!
!!    Calculate expectation value
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none
!  
      class(ccs), intent(in) :: wf
!
      complex(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: A
!
      complex(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: density
!
      complex(dp) :: expectation_value
!
   end function calculate_expectation_value_ccs_complex
!
!
   module subroutine calculate_energy_ccs_complex(wf)
!!
!!    Calculate energy_complex
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
   end subroutine calculate_energy_ccs_complex
!
!
   module subroutine calculate_energy_omega_term_ccs_complex(wf)
!!
!!    Calculate energy_complex omega term (CCS)
!!    Written by Andreas Skeidsvoll, Jan 2019
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
   end subroutine calculate_energy_omega_term_ccs_complex
!
!
   module subroutine calculate_energy_length_dipole_term_ccs_complex(wf, electric_field)
!!
!!    Calculate energy_complex length dipole term (CCS)
!!    Written by Andreas Skeidsvoll, Jan 2019
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      complex(dp), dimension(3), intent(in) :: electric_field
!
   end subroutine calculate_energy_length_dipole_term_ccs_complex