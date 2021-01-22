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
module libint_initialization
!
   use kinds
   use iso_c_binding
!
   include "libint_initialization_cdef.F90"
!
contains
!
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
      use atomic_center_class, only: atomic_center
      use parameters, only: angstrom_to_bohr
      use iso_c_binding, only: c_int, c_double, c_char, c_null_char
!
      implicit none
!
      integer(c_int), parameter :: max_length = 100
!
      type(atomic_center), dimension(:), intent(in) :: atoms
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
         atomic_numbers_c(i) = int(atoms(i)%number_, c_int)
!
         atomic_coordinates_c(1,i) = real(atoms(i)%coordinates(1)*angstrom_to_bohr, c_double)
         atomic_coordinates_c(2,i) = real(atoms(i)%coordinates(2)*angstrom_to_bohr, c_double)
         atomic_coordinates_c(3,i) = real(atoms(i)%coordinates(3)*angstrom_to_bohr, c_double)
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
!
end module libint_initialization
