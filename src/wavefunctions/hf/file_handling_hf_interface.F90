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
   module subroutine save_orbital_coefficients_hf(wf)
!!
!!    Save orbital coefficients
!!    Written by Eirik F. Kjønstad, Oct 2018
!!
      implicit none
!
      class(hf), intent(inout) :: wf
!
   end subroutine save_orbital_coefficients_hf
!
!
   module subroutine read_orbital_coefficients_hf(wf)
!!
!!    Save orbital coefficients
!!    Written by Eirik F. Kjønstad, Oct 2018
!!
      implicit none
!
      class(hf), intent(inout) :: wf
!
   end subroutine read_orbital_coefficients_hf
!
!
   module subroutine save_orbital_energies_hf(wf)
!!
!!    Save orbital energies
!!    Written by Eirik F. Kjønstad, Oct 2018
!!
      implicit none
!
      class(hf), intent(inout) :: wf
!
   end subroutine save_orbital_energies_hf
!
!
   module subroutine read_orbital_energies_hf(wf)
!!
!!    Save orbital energies
!!    Written by Eirik F. Kjønstad, Oct 2018
!!
      implicit none
!
      class(hf), intent(inout) :: wf
!
   end subroutine read_orbital_energies_hf
!
!
   module subroutine write_scf_restart_hf(wf)
!!
!!    Write HF restart file
!!    Written by Linda Goletto, Oct 2019
!!
!!    Writes a file used for consistency checks when restarting
!!
      implicit none
!
      class(hf) :: wf
!
   end subroutine write_scf_restart_hf
!
!
   module subroutine write_orbital_information_hf(wf)
!!
!!    Write HF information file
!!    Written by Linda Goletto, Oct 2019
!!
!!    Writes orbital information
!!    Used for CC and MP2 and for restarting MLHF
!!
      implicit none
!
      class(hf) :: wf
!
   end subroutine write_orbital_information_hf
!
!
   module subroutine save_ao_density_hf(wf)
!!
!!    Save AO density
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Save the AO density based
!!    on the current orbital coefficient matrix (or matrices).
!!
      implicit none
!
      class(hf) :: wf
!
   end subroutine save_ao_density_hf
