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
   module subroutine prepare_for_multiplier_equation_cc2(wf)
!!
!!    Prepare for jacobian transpose transformation
!!    Written by Tor S. Haugland, Oct 2019
!!
      implicit none
!
      class(cc2), intent(inout) :: wf
!
   end subroutine prepare_for_multiplier_equation_cc2
!
!
   module subroutine construct_multiplier_equation_cc2(wf, equation)
!!
!!    Construct multiplier equation 
!!    Written by Sarai D. Folkestad, Feb 2019
!!
      implicit none 
!
      class(cc2), intent(in) :: wf 
!
      real(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: equation 
!
   end subroutine construct_multiplier_equation_cc2