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
   module subroutine save_orbitals_hf(wf)
!!
!!    Save orbitals
!!    Written by Eirik F. Kjønstad, Oct 2018
!!
      implicit none
!
      class(hf), intent(inout) :: wf
!
      call wf%orbital_file%open_('write')
      call wf%orbital_file%write_(int(wf%ao%n,kind=i64))
      call wf%orbital_file%write_(int(wf%n_mo,kind=i64))
      call wf%orbital_file%write_(wf%orbital_energies, wf%n_mo)
      call wf%orbital_file%write_(wf%orbital_coefficients, wf%ao%n*wf%n_mo)
      call wf%orbital_file%close_('keep')
!
   end subroutine save_orbitals_hf
!
!
   module subroutine read_orbitals_hf(wf)
!!
!!    Save orbitals
!!    Written by Eirik F. Kjønstad, Oct 2018
!!
      implicit none
!
      class(hf), intent(inout) :: wf
!
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
      call wf%orbital_file%read_(wf%orbital_energies, wf%n_mo)
      call wf%orbital_file%read_(wf%orbital_coefficients, wf%ao%n*wf%n_mo)
      call wf%orbital_file%close_('keep')
!
   end subroutine read_orbitals_hf
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
!
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
   module subroutine save_tdhf_vector_hf(wf, vector, energy, I)
!!
!!    Save TDHF vector
!!    Written by Sarai D. Folkestad, 2021
!!
      implicit none
!
      class(hf), intent(inout) :: wf
!
      real(dp), dimension(wf%n_o*wf%n_v), intent(in) :: vector
!
      real(dp), intent(in) :: energy
!
      integer, intent(in) :: I
!
      integer(i64) :: n
!
      n = int(wf%n_o*wf%n_v, kind=i64)
!
      call wf%tdhf_files(I)%open_('write', 'rewind')
!
      call wf%tdhf_files(I)%write_(energy)
      call wf%tdhf_files(I)%write_(n)
      call wf%tdhf_files(I)%write_(vector, wf%n_o*wf%n_v)
!
      call wf%tdhf_files(I)%close_
!
   end subroutine save_tdhf_vector_hf
!
!
   module subroutine read_tdhf_vector_hf(wf, vector, I)
!!
!!    Read TDHF vector
!!    Written By Sarai D. Folkestad, 2021
!!
      implicit none
!
      class(hf), intent(inout) :: wf
      real(dp), dimension(wf%n_o*wf%n_v), intent(out) :: vector
!
      integer, intent(in) :: I
      integer(i64) :: n
!
      call wf%tdhf_files(I)%open_('read', 'rewind')
      call wf%tdhf_files(I)%read_(n, dp + 1)
!
      if (int(n) /= wf%n_o*wf%n_v) &
            call output%error_msg('Wrong number of parameters in (a0): (i0)', &
                                  chars=[wf%tdhf_files(I)%get_name()], ints=[int(n)])
!
      call wf%tdhf_files(I)%read_(vector, wf%n_o*wf%n_v)
!
      call wf%tdhf_files(I)%close_
!
   end subroutine read_tdhf_vector_hf
!
!
   module subroutine read_tdhf_energy_hf(wf, I)
!!
!!    Read TDHF energy
!!    Written by Sarai D. Folkestad, 2021
!!
      implicit none
!
      class(hf), intent(inout) :: wf
!
      integer, intent(in) :: I
!
      call wf%tdhf_files(I)%open_('read', 'rewind')
!
      call wf%tdhf_files(I)%read_(wf%tdhf_excitation_energies(I))
!
      call wf%tdhf_files(I)%close_
!
   end subroutine read_tdhf_energy_hf
!
!
   module subroutine initialize_tdhf_files_hf(wf)
!!
!!    Initialize TDHF files
!!    Written by Sarai D. Folkestad, 2021
!!
!
      implicit none
!
      class(hf), intent(inout) :: wf
!
      character(len=21) :: file_name
      integer :: state
!
      allocate(wf%tdhf_files(wf%n_tdhf_vectors))
!
      do state = 1, wf%n_tdhf_vectors
!
         write(file_name,'(a,i3.3)') 'tdhf_eigenvectors_', state
         wf%tdhf_files(state) = stream_file(trim(file_name))
!
      end do
!
   end subroutine initialize_tdhf_files_hf
!
!
   module subroutine write_molden_file_hf(wf)
!!
!!    Write Molden file
!!    Written by Alexander C. Paul, May 2021
!!
!!    Writes file readable by "molden" (https://www3.cmbi.umcn.nl/molden/)
!!
      implicit none
!
      class(hf), intent(in) :: wf
      type(output_file), allocatable :: molden
!
!     Map aos from libint ordering to molden ordering
      integer, dimension(:), allocatable  :: ao_map
      real(dp), dimension(:), allocatable :: scaling_factor
      integer p, q, i
!
      real(dp) :: occupation
!
      call mem%alloc(ao_map, wf%ao%n)
      call mem%alloc(scaling_factor, wf%ao%n)
!
      ao_map = wf%ao%get_molden_ao_indices()
      scaling_factor = wf%ao%get_ao_normalization_factors()
!
      molden = output_file('eT.molden')
      call molden%open_
!
!     Print header
!
      call molden%printf('m', '[Molden Format]', fs='(t1,a)')
!
!     Print Geometry
!
      call molden%printf('m', '[ATOMS] AU', fs='(t1,a)')
      call wf%ao%print_molden_geometry(molden)
!
!     Print basis set information per atom
!     Atom_number (as in geometry)
!     shell_label number_of_primitives
!     exponent coefficient
!
      call molden%printf('m', '[GTO]', fs='(t1,a)')
      call wf%ao%print_basis_set_molden(molden)
!
!     Print MO coefficients
!
      call molden%printf('m', '[MO]', fs='(t1,a)')
!
      do p = 1, wf%n_mo
!
         if (p .le. wf%n_o) then
            occupation = two
         else
            occupation = zero
         end if
!
         call molden%printf('m', 'Sym=X', fs='(t1,a)')
         call molden%printf('m', 'Ene=(f17.10)', fs='(t1,a)', reals=[wf%orbital_energies(p)])
         call molden%printf('m', 'Spin=Alpha', fs='(t1,a)')
         call molden%printf('m', 'Occup=(f6.4)', reals=[occupation], fs='(t1,a)')
!
         do q = 1, wf%ao%n
!
            i = ao_map(q)
!
            call molden%printf('m', '(i6)   (f17.12)', fs='(t1,a)', ints=[q], &
                               reals=[scaling_factor(i)*wf%orbital_coefficients(i,p)])
!
         end do
      end do
!
      call mem%dealloc(ao_map, wf%ao%n)
      call mem%dealloc(scaling_factor, wf%ao%n)
      call molden%close_
!
   end subroutine write_molden_file_hf
!
!
end submodule file_handling_hf
