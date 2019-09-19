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
interface
!
   subroutine set_coulomb_precision_c(prec) bind(C, name='set_coulomb_precision')
!
      use iso_c_binding
      implicit none
!
      real(c_double) :: prec
!
   end subroutine set_coulomb_precision_c
!
   subroutine initialize_basis_c(basisset,filename,cartesian_gaussians) bind(C, name='initialize_basis')
!
      use iso_c_binding, only: C_CHAR, C_LONG, C_BOOL
      implicit none
!
      character(kind = c_char) :: basisset(*)
      character(kind = c_char) :: filename(*)
      logical(kind = c_bool), value :: cartesian_gaussians
!
   end subroutine initialize_basis_c
!
   subroutine reset_basis_c() bind(C, name='reset_basis')
!
      implicit none
!
   end subroutine reset_basis_c
!
   subroutine initialize_shell2atom_c() bind(C, name='initialize_shell2atom')
!
      implicit none
!
   end subroutine initialize_shell2atom_c
!
   subroutine initialize_atoms_c(name) bind(C, name='initialize_atoms')
!
      use iso_c_binding, only: C_CHAR
      implicit none
!
      character(kind = c_char) :: name(*)
!
   end subroutine initialize_atoms_c
!
   subroutine initialize_libint_c() bind(C, name='initialize_libint')
!
      use iso_c_binding
      implicit none
!
   end subroutine initialize_libint_c
!
   subroutine finalize_libint_c() bind(C, name='finalize_libint')
!
      use iso_c_binding
      implicit none
!
   end subroutine finalize_libint_c
!
   subroutine initialize_coulomb_c() bind(C, name='initialize_coulomb')
!
      use iso_c_binding
      implicit none
!
   end subroutine initialize_coulomb_c
!
   subroutine initialize_kinetic_c() bind(C, name='initialize_kinetic')
!
      use iso_c_binding
      implicit none
!
   end subroutine initialize_kinetic_c
!
   subroutine initialize_overlap_c() bind(C, name='initialize_overlap')
!
      use iso_c_binding
      implicit none
!
   end subroutine initialize_overlap_c
!
   subroutine initialize_dipole_c() bind(C, name='initialize_dipole')
!
      use iso_c_binding
      implicit none
!
   end subroutine initialize_dipole_c
!
   subroutine initialize_quadrupole_c() bind(C, name='initialize_quadrupole')
!
      use iso_c_binding
      implicit none
!
   end subroutine initialize_quadrupole_c
!
   subroutine initialize_nuclear_c() bind(C, name='initialize_nuclear')
!
      use iso_c_binding
      implicit none
!
   end subroutine initialize_nuclear_c
!
   subroutine initialize_potential_c(charges,coordinates,n_points) bind(C, name='initialize_potential')
!
      use iso_c_binding
      implicit none
      real(c_double)  :: charges(*)
      real(c_double)  :: coordinates(*)
      integer(c_int)  :: n_points
!
   end subroutine initialize_potential_c
!
end interface
