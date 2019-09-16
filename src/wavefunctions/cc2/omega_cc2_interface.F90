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
   module subroutine construct_omega_cc2(wf, omega)
!!
!!    Construct omega
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, Jan 2019
!!
!!    Directs the construction of the projection vector < mu | exp(-T) H exp(T) | R >
!!    for the current wavefunction amplitudes.
!!
      implicit none
!
      class(cc2), intent(inout) :: wf
!
      real(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: omega
!
   end subroutine construct_omega_cc2
!
!
   module subroutine construct_omega2_cc2(wf, omega2, t)
!!
!!    Construct Omega2
!!    Written by Sarai D. Folkestad, Sep 2019
!!
!!    Constructs the doubles part of omega for MLCC2
!!
      implicit none
!
      class(cc2) :: wf
!
      real(dp), dimension(wf%n_t2), intent(out) :: omega2
      real(dp), dimension(wf%n_t2), intent(in) :: t
!
   end subroutine construct_omega2_cc2
