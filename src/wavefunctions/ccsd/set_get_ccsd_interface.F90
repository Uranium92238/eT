!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2020 the authors of eT
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
   module subroutine set_amplitudes_ccsd(wf, amplitudes)
!!
!!    Set amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(ccsd) :: wf
      real(dp), dimension(wf%n_gs_amplitudes), intent(in) :: amplitudes
!
   end subroutine set_amplitudes_ccsd
!
!
   module subroutine get_amplitudes_ccsd(wf, amplitudes)
!!
!!    Get amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
      real(dp), dimension(wf%n_gs_amplitudes) :: amplitudes
!
   end subroutine get_amplitudes_ccsd
!
!
   module subroutine set_multipliers_ccsd(wf, multipliers)
!!
!!    Set multipliers
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Nov 2018
!!
      implicit none
!
      class(ccsd) :: wf
      real(dp), dimension(wf%n_gs_amplitudes), intent(in) :: multipliers
!
   end subroutine set_multipliers_ccsd
!
!
   module subroutine get_multipliers_ccsd(wf, multipliers)
!!
!!    Get multipliers
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Nov 2018
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
      real(dp), dimension(wf%n_gs_amplitudes) :: multipliers
!
   end subroutine get_multipliers_ccsd
