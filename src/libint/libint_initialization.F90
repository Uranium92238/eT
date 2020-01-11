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
   subroutine initialize_atoms(mol_name)
!
      use iso_c_binding, only: c_char, c_null_char
      implicit none
!
      character(len=*) :: mol_name
      character(len=40) :: mol_name_temp
      character(kind=c_char),dimension(40) :: cpp_temp_mol_name
!
      integer :: j
!
      mol_name_temp = trim(mol_name)//c_null_char
!
      do j=1,len_trim(mol_name_temp)
!
         cpp_temp_mol_name(j) = mol_name_temp(j:j)  
!
      enddo
!
      call initialize_atoms_c(cpp_temp_mol_name)
!
   end subroutine initialize_atoms
!
   subroutine initialize_basis(basis_set, file_name, cartesian_gaussians_int)
!
      use iso_c_binding, only: c_char, c_null_char, c_int
      implicit none
!
      character(len=*) :: basis_set
      character(len=*) :: file_name
!
      character(len=40) :: basis_set_temp
      character(len=40) :: file_name_temp
!
      character(kind=c_char),dimension(40) :: cpp_temp_basis
      character(kind=c_char),dimension(40) :: cpp_temp_file
!
      integer :: cartesian_gaussians_int
      integer(c_int) :: cartesian_gaussians_int_c
!
      integer :: j
!
      basis_set_temp = trim(basis_set)//c_null_char
!
      do j=1,len_trim(basis_set_temp)
!
         cpp_temp_basis(j) = basis_set_temp(j:j)
!
      enddo
!
      file_name_temp = trim(file_name)//c_null_char
!
      do j=1,len_trim(file_name_temp)
!
         cpp_temp_file(j) = file_name_temp(j:j)  
!
      enddo
!
      cartesian_gaussians_int_c = int(cartesian_gaussians_int, kind=c_int)
!     
      call initialize_basis_c(cpp_temp_basis, cpp_temp_file, cartesian_gaussians_int_c)
!
   end subroutine initialize_basis
!
   subroutine initialize_potential(charges,coordinates,n_points)
!
      use iso_c_binding
      implicit none
      real(dp) :: charges(*)
      real(dp) :: coordinates(*)
      integer(i6) :: n_points
!
      call initialize_potential_c(charges,coordinates,n_points)
!
   end subroutine initialize_potential
!
end module libint_initialization
