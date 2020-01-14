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
   module subroutine calculate_energy_mp2(wf)
!!
!!    Calculate energy
!!    Written by Andreas Skeidsvoll, 2018
!!
!!    Calculates the MP2 energy from HF energy, E_HF, vovo integrals, g_aibj, 
!!    and the orbital energies, eps. The total MP2 energy is calculated as
!!
!!       E = E_HF - sum_aibj g_aibj*L_aibj/(eps(a)+eps(b)-eps(i)-eps(j))
!!
!!    where
!!
!!       L_aibj = 2*g_aibj - g_ajbi.
!!
!!    On entry, it is assumed that the energy is equal to the HF energy 
!!    (i.e. the routine only adds the correction to the energy variable).
!!
      implicit none
!
      class(mp2), intent(inout) :: wf 
!
   end subroutine calculate_energy_mp2
