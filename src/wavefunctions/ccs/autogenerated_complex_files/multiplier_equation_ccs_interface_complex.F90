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
   module subroutine prepare_for_multiplier_equation_ccs_complex(wf)
!!
!!    Prepare for the construction of the multipliers
!!    Written by Alexander C. Paul, July 2019
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
   end subroutine prepare_for_multiplier_equation_ccs_complex
!
!
   module subroutine construct_multiplier_equation_ccs_complex(wf, equation)
!!
!!    Construct multiplier equation
!!    Written by Eirik F. Kjønstad, Oct 2018
!!
!!    Constructs
!!
!!       t-bar^T A + eta,
!!
!!    and places the result in 'equation'.
!!
      implicit none
!
      class(ccs), intent(in) :: wf
      complex(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: equation
!
   end subroutine construct_multiplier_equation_ccs_complex
!
!
   module subroutine construct_eta_ccs_complex(wf, eta)
!!
!!    Construct eta
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
      implicit none
!
      class(ccs), intent(in) :: wf
      complex(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: eta
!
   end subroutine construct_eta_ccs_complex
