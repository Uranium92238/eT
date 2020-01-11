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
   module subroutine save_orbital_coefficients_uhf(wf)
!!
!!    Save orbital coefficients
!!    Written by Eirik F. Kjønstad, Oct 2018
!!
      implicit none
!
      class(uhf), intent(inout) :: wf
!
   end subroutine save_orbital_coefficients_uhf
!
!
   module subroutine read_orbital_coefficients_uhf(wf)
!!
!!    Read orbital coefficients
!!    Written by Eirik F. Kjønstad, Oct 2018
!!
      implicit none
!
      class(uhf), intent(inout) :: wf
!
   end subroutine read_orbital_coefficients_uhf
!
!
   module subroutine save_orbital_energies_uhf(wf)
!!
!!    Save orbital energies
!!    Written by Eirik F. Kjønstad, Mar 2019
!!
      implicit none
!
      class(uhf), intent(inout) :: wf
!
   end subroutine save_orbital_energies_uhf
!
!
   module subroutine read_orbital_energies_uhf(wf)
!!
!!    Save orbital energies
!!    Written by Eirik F. Kjønstad, Mar 2019
!!
      implicit none
!
      class(uhf), intent(inout) :: wf
!
   end subroutine read_orbital_energies_uhf
!
!
   module subroutine save_ao_density_uhf(wf)
!!
!!    Save AO density
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Save the AO density (or densities, if unrestricted) based
!!    on the current orbital coefficient matrix (or matrices).
!!
      implicit none
!
      class(uhf) :: wf
!
   end subroutine save_ao_density_uhf
