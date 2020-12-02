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
module libint_initialization
!
   use kinds
   use iso_c_binding
!
   include "libint_initialization_cdef.F90"
!
contains
!
   subroutine set_coulomb_precision(prec)
!
      implicit none
!
      real(dp) :: prec
!
      call set_coulomb_precision_c(prec)
!
   end subroutine set_coulomb_precision
!
   subroutine initialize_libint()
!
      implicit none
!
      call initialize_libint_c()
!
   end subroutine initialize_libint
!
   subroutine finalize_libint()
!
      implicit none
!
      call finalize_libint_c()
!
   end subroutine finalize_libint
!
   subroutine initialize_overlap()
!
      implicit none
!
      call initialize_overlap_c()
!
   end subroutine initialize_overlap
!
   subroutine initialize_dipole()
!
      implicit none
!
      call initialize_dipole_c()
!
   end subroutine initialize_dipole
!
   subroutine initialize_quadrupole()
!
      implicit none
!
      call initialize_quadrupole_c()
!
   end subroutine initialize_quadrupole
!
   subroutine initialize_coulomb()
!
      implicit none
!
      call initialize_coulomb_c()
!
   end subroutine initialize_coulomb
!
   subroutine initialize_kinetic()
!
      implicit none
!
      call initialize_kinetic_c()
!
   end subroutine initialize_kinetic
!
   subroutine initialize_nuclear()
!
      implicit none
!
      call initialize_nuclear_c()
!
   end subroutine initialize_nuclear
!
   subroutine export_geometry_and_basis_to_libint(atoms)
!!
!!    Export geometry and basis
!!    Written by Rolf Heilemann Myhre, Mar. 2020
!!
!!    Uses the input array of atoms to generate
!!    a C integer array with atomic numbers,
!!    a C double array with coordinates,
!!    a C char array with basis set names,
!!    a C pointer array to the basis set names, and
!!    a C integer array which is 1 if basis on atom is cartesian and 0 otherwise
!!    and pass them on to the C side
!!
!
      use atomic_class, only: atomic
      use parameters, only: angstrom_to_bohr
      use iso_c_binding, only: c_int, c_double, c_char, c_null_char
!
      implicit none
!
      integer(c_int), parameter :: max_length = 100
!
      type(atomic), dimension(:), intent(in) :: atoms
!
      integer(c_int), dimension(:), allocatable :: atomic_numbers_c
      real(c_double), dimension(:,:), allocatable :: atomic_coordinates_c
!
      integer :: i
!
      integer :: n_atoms
      integer(c_int) :: n_atoms_c
!
      character(len=max_length, kind=c_char), dimension(:), allocatable :: basis_sets_c
!
      integer(kind=c_int), dimension(:), allocatable :: cartesians_c
!
      n_atoms = size(atoms)
      n_atoms_c = int(n_atoms, c_int)
!
      allocate(atomic_numbers_c(n_atoms))
      allocate(atomic_coordinates_c(3, n_atoms))
      allocate(basis_sets_c(n_atoms))
      allocate(cartesians_c(n_atoms))
!
!
!     Construct C arrays with atom numbers and coordinates
      do i =1,n_atoms
!
         call atoms(i)%set_number()
         atomic_numbers_c(i) = int(atoms(i)%number_, c_int)
!
         atomic_coordinates_c(1,i) = real(atoms(i)%x*angstrom_to_bohr, c_double)
         atomic_coordinates_c(2,i) = real(atoms(i)%y*angstrom_to_bohr, c_double)
         atomic_coordinates_c(3,i) = real(atoms(i)%z*angstrom_to_bohr, c_double)
!
         basis_sets_c(i) = trim(atoms(i)%basis)//c_null_char
!
         if (atoms(i)%cartesian) then
            cartesians_c(i) = int(1, c_int)
         else
            cartesians_c(i) = int(0, c_int)
         endif
!
      enddo
!
!     Export the geometry and basis
!     We only pass on the pointers to basis_sets_c because chars are difficult in C
      call export_geometry_and_basis_to_libint_c(n_atoms_c, &
                                                 atomic_numbers_c, &
                                                 atomic_coordinates_c, &
                                                 basis_sets_c, &
                                                 max_length, &
                                                 cartesians_c)
!
!
   end subroutine export_geometry_and_basis_to_libint
!
   subroutine initialize_potential(charges,coordinates,n_points)
!
      use iso_c_binding, only: c_int
      implicit none
      real(dp) :: charges(*)
      real(dp) :: coordinates(*)
      integer(c_int) :: n_points
!
      call initialize_potential_c(charges,coordinates,n_points)
!
   end subroutine initialize_potential
!
end module libint_initialization
