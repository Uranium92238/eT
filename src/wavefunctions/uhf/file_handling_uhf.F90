!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2022 the authors of eT
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
submodule (uhf_class) file_handling_uhf
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
   module subroutine save_orbitals_uhf(wf)
!!
!!    Save orbital coefficients
!!    Written by Eirik F. Kjønstad, Oct 2018
!!
      implicit none
!
      class(uhf), intent(inout) :: wf
!
      call wf%orbital_file%open_('write')
      call wf%orbital_file%write_(int(wf%ao%n,kind=i64))
      call wf%orbital_file%write_(int(wf%n_mo,kind=i64))
      call wf%orbital_file%write_(wf%orbital_energies_a, wf%n_mo)
      call wf%orbital_file%write_(wf%orbital_coefficients_a, wf%ao%n*wf%n_mo)
      call wf%orbital_file%write_(wf%orbital_energies_b, wf%n_mo)
      call wf%orbital_file%write_(wf%orbital_coefficients_b, wf%ao%n*wf%n_mo)
      call wf%orbital_file%close_('keep')
!
   end subroutine save_orbitals_uhf
!
!
   module subroutine read_orbitals_uhf(wf)
!!
!!    Read orbital coefficients
!!    Written by Eirik F. Kjønstad, Oct 2018
!!
      implicit none
!
      class(uhf), intent(inout) :: wf
      integer :: n_ao, n_mo
      integer(i64) :: n
!
      call wf%orbital_file%open_('read', 'rewind')
!
      call wf%orbital_file%read_(n)
      n_ao = int(n)
      call wf%orbital_file%read_(n)
      n_mo = int(n)
!
      if (n_ao .ne. wf%ao%n) call output%error_msg('Number of AOs does not match')
      if (n_mo .ne. wf%n_mo) call output%error_msg('Number of MOs does not match')
!
      call wf%orbital_file%read_(wf%orbital_energies_a, wf%n_mo)
      call wf%orbital_file%read_(wf%orbital_coefficients_a, wf%ao%n*wf%n_mo)
      call wf%orbital_file%read_(wf%orbital_energies_b, wf%n_mo)
      call wf%orbital_file%read_(wf%orbital_coefficients_b, wf%ao%n*wf%n_mo)
!
      call wf%orbital_file%close_('keep')
!
   end subroutine read_orbitals_uhf
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
      type(sequential_file) :: ao_density_file
      type(sequential_file) :: ao_density_file_a
      type(sequential_file) :: ao_density_file_b
!
      ao_density_file   = sequential_file('ao_density')
      ao_density_file_a = sequential_file('ao_density_a')
      ao_density_file_b = sequential_file('ao_density_b')
!
      call ao_density_file%open_('write', 'rewind')
      call ao_density_file%write_(wf%ao_density, wf%ao%n*wf%ao%n)
      call ao_density_file%close_
!
      call ao_density_file_a%open_('write', 'rewind')
      call ao_density_file_a%write_(wf%ao_density_a, wf%ao%n*wf%ao%n)
      call ao_density_file_a%close_
!
      call ao_density_file_b%open_('write', 'rewind')
      call ao_density_file_b%write_(wf%ao_density_b, wf%ao%n*wf%ao%n)
      call ao_density_file_b%close_
!
   end subroutine save_ao_density_uhf
!
!
end submodule file_handling_uhf
