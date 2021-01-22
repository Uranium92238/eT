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
submodule (hf_class) file_handling_hf
!
!!
!!    File handling submodule
!!
!!    Gathers routines that save wavefunction parameters to file,
!!    and reads them from file, plus other routines related to the 
!!    handling of the files that belong to the wavefunction.
!!
!
   implicit none
!
!
contains
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
      call wf%orbital_coefficients_file%open_('write', 'rewind')
!
      call wf%orbital_coefficients_file%write_(wf%orbital_coefficients, wf%ao%n*wf%n_mo)
!
      call wf%orbital_coefficients_file%close_
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
      call wf%orbital_coefficients_file%open_('read', 'rewind')
!
      call wf%orbital_coefficients_file%read_(wf%orbital_coefficients, wf%ao%n*wf%n_mo)
!
      call wf%orbital_coefficients_file%close_
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
      call wf%orbital_energies_file%open_('write', 'rewind')
!
      call wf%orbital_energies_file%write_(wf%orbital_energies, wf%n_mo)
!
      call wf%orbital_energies_file%close_
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
      call wf%orbital_energies_file%open_('read', 'rewind')
!
      call wf%orbital_energies_file%read_(wf%orbital_energies, wf%n_mo)
!
      call wf%orbital_energies_file%close_
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
      call wf%restart_file%open_('write', 'rewind')
!
      call wf%restart_file%write_(wf%ao%n)
      call wf%restart_file%write_(wf%n_densities)
      call wf%restart_file%write_(wf%ao%get_n_electrons())
!
      call wf%restart_file%close_
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
      type(sequential_file) :: CC_orbitals_file, CC_orbital_energies_file
!
      wf%orbital_information_file = sequential_file('orbital_information')
      call wf%orbital_information_file%open_('write', 'rewind')
!
      call wf%orbital_information_file%write_(wf%n_o)
      call wf%orbital_information_file%write_(wf%n_v)
      call wf%orbital_information_file%write_(wf%ao%n)
      call wf%orbital_information_file%write_(wf%n_mo)
      call wf%orbital_information_file%write_(wf%energy)
!
      call wf%orbital_information_file%close_
!
      CC_orbitals_file = sequential_file('cc_orbital_coefficients')
      call CC_orbitals_file%open_('write', 'rewind')
!
      call CC_orbitals_file%write_(wf%orbital_coefficients, wf%ao%n*wf%n_mo)
!
      call CC_orbitals_file%close_('keep')
!
      CC_orbital_energies_file = sequential_file('cc_orbital_energies')
      call CC_orbital_energies_file%open_('write', 'rewind')
!
      call CC_orbital_energies_file%write_(wf%orbital_energies, wf%n_mo)
!
      call CC_orbital_energies_file%close_('keep')
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
      type(sequential_file) :: ao_density_file
!
      ao_density_file = sequential_file('ao_density')
      call ao_density_file%open_('write', 'rewind')
!
      call ao_density_file%write_(wf%ao_density, wf%ao%n*wf%ao%n)
!
      call ao_density_file%close_
!
   end subroutine save_ao_density_hf
!
!
end submodule file_handling_hf