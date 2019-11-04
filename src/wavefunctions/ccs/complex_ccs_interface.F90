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
   module subroutine make_complex_ccs(wf)
!!
!!    Make complex (CCS)
!!    Written by Andreas Skeidsvoll, Sep 2019
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
   end subroutine make_complex_ccs
!
!
   module subroutine make_ccs_complex_ccs(wf)
!!
!!    Make CCS complex (CCS)
!!    Written by Andreas Skeidsvoll, Sep 2019
!!
!!    Allocates complex CCS variables, puts the real variables into the real part of the complex
!!    variables and zero into the imaginary part, and deallocates the real variables.
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
   end subroutine make_ccs_complex_ccs
!
!
   module subroutine construct_complex_time_derivative_ccs(wf, ddt_amplitudes_multipliers)
!!
!!    Construct time derivative of amplitudes and multipliers (CCS)
!!    Written by Alice Balbi and Andreas Skeidsvoll, Oct 2018
!!
!!    Returns the time derivative of the CCS amplitudes and multipliers. Does not have to be
!!    overwritten in descendants.
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
      complex(dp), intent(out), dimension(2*wf%n_gs_amplitudes) :: ddt_amplitudes_multipliers
!
   end subroutine construct_complex_time_derivative_ccs
!
!
   module subroutine construct_complex_time_derivative_amplitudes_ccs(wf, ddt_amplitudes)
!!
!!    Construct time derivative of amplitudes (CCS)
!!    Written by Alice Balbi and Andreas Skeidsvoll, Oct 2018
!!
!!    Returns the time derivative of the CCS amplitudes (complex). The time derivative
!!    of the amplitudes is, according to Koch and Jørgensen (1990), given by
!!
!!       ddt_amplitudes = -i*omega
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
      complex(dp), dimension(wf%n_gs_amplitudes), intent(out) :: ddt_amplitudes
!
   end subroutine construct_complex_time_derivative_amplitudes_ccs
!
!
   module subroutine construct_complex_time_derivative_multipliers_ccs(wf, ddt_multipliers)
!!
!!    Construct time derivative multipliers (CCS)
!!    Written by Alice Balbi and Andreas Skeidsvoll, Oct 2018
!!
!!    Returns the time derivative of the CCS multipliers (complex). The time derivative
!!    of the multipliers is, according to Koch and Jørgensen (1990), given by
!!
!!       ddt_multipliers = i*multiplier_equation
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
      complex(dp), dimension(wf%n_gs_amplitudes), intent(out) :: ddt_multipliers
!
   end subroutine construct_complex_time_derivative_multipliers_ccs