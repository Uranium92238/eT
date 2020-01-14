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
   module subroutine calculate_energy_lowmem_cc2(wf)
!!
!!     Calculate energy (lowmem CC2)
!!     Written by Sarai D. Folkestad, Eirik F. Kjønstad,
!!     Andreas Skeidsvoll, 2018
!!
!!     Calculates the lowmem CC2 energy. This is only equal to the actual
!!     energy when the ground state equations are solved, of course.
!!
      implicit none
!
      class(lowmem_cc2), intent(inout) :: wf
!
   end subroutine calculate_energy_lowmem_cc2
