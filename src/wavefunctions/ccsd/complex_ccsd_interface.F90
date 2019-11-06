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
   module subroutine make_complex_ccsd(wf)
!!
!!    Make complex (CCSD)
!!    Written by Andreas Skeidsvoll, Sep 2019
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
!
   end subroutine make_complex_ccsd
!
!
   module subroutine construct_complex_time_derivative_amplitudes_ccsd(wf, ddt_amplitudes)
!!
!!    Construct complex time derivative amplitudes (CCSD)
!!    Written by Alice Balbi and Andreas Skeidsvoll, Oct 2018
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
      complex(dp), dimension(wf%n_gs_amplitudes), intent(out) :: ddt_amplitudes
!
   end subroutine construct_complex_time_derivative_amplitudes_ccsd
!
!
   module subroutine construct_complex_time_derivative_multipliers_ccsd(wf, ddt_multipliers)
!!
!!    Construct complex time derivative multipliers (CCSD)
!!    Written by Alice Balbi and Andreas Skeidsvoll, Oct 2018
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
      complex(dp), dimension(wf%n_gs_amplitudes), intent(out) :: ddt_multipliers
!
   end subroutine construct_complex_time_derivative_multipliers_ccsd