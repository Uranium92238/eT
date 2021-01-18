!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2021 the authors of eT
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
   module subroutine prepare_for_multiplier_equation_ccsd_complex(wf)
!!
!!    Prepare for multiplier equation
!!    Written by Tor S. Haugland, 2019
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
!
   end subroutine prepare_for_multiplier_equation_ccsd_complex
!
!
   module subroutine construct_eta_ccsd_complex(wf, eta)
!!
!!    Construct eta (CCSD)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
!!    Note: the routine assumes that eta is initialized and that the Fock matrix
!!    has been constructed.
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
      complex(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: eta
!
   end subroutine construct_eta_ccsd_complex
