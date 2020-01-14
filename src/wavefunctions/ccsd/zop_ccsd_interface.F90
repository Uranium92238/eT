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
   module subroutine calculate_energy_ccsd(wf)
!!
!!    Calculate energy (CCSD)
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad,
!!    Andreas Skeidsvoll, 2018
!!
!!    Calculates the CCSD energy. This is only equal to the actual
!!    energy when the ground state equations are solved, of course.
!!
!!       E = E_hf + sum_aibj (t_ij^ab + t_i^a t_j^b) L_iajb
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
!
   end subroutine calculate_energy_ccsd
