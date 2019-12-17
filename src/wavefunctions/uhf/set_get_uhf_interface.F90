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
   module subroutine set_ao_fock_uhf(wf, F)
!!
!!    Set AO Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Sets the AO Fock.
!!
      implicit none
!
      class(uhf) :: wf
      real(dp), dimension(wf%n_ao*(wf%n_ao + 1)/2, wf%n_densities), intent(in) :: F ! Packed
!
   end subroutine set_ao_fock_uhf
!
!
   module subroutine get_ao_fock_uhf(wf, F)
!!
!!    Get AO Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Returns the AO Fock
!!
      implicit none
!
      class(uhf), intent(in) :: wf
      real(dp), dimension(wf%n_ao*(wf%n_ao+1)/2, wf%n_densities), intent(inout) :: F ! Packed
!
   end subroutine get_ao_fock_uhf
!
!
   module subroutine get_ao_density_sq_uhf(wf, D)
!!
!!    Get AO density squared
!!    Written by Eirik F. Kjønstad, Nov 2018
!!
!!    Returns the unpacked AO density matrix D
!!    (or density matrices in descendants, see overwriting routines)
!!
      implicit none
!
      class(uhf), intent(in) :: wf
      real(dp), dimension(wf%n_ao**2, wf%n_densities), intent(inout) :: D
!
   end subroutine get_ao_density_sq_uhf
