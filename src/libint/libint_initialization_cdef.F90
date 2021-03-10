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
interface
!
!
   subroutine export_geometry_and_basis_to_libint_c(nAtoms, &
                                                    atomicNumbers, &
                                                    atomicCoordinates, &
                                                    basisSets, &
                                                    maxLen, &
                                                    cartesians) &
               bind(C, name='export_geometry_and_basis_to_libint')
!
      use iso_c_binding, only: c_int, c_double, c_char
      implicit none
!
      integer(c_int), value                        :: nAtoms
      integer(c_int), dimension(*), intent(in)     :: atomicNumbers
      real(c_double), dimension(*), intent(in)     :: atomicCoordinates
      character(c_char), dimension(*), intent(in)  :: basisSets
      integer(c_int), value                        :: maxLen
      integer(c_int), dimension(*), intent(in)     :: cartesians
!
   end subroutine export_geometry_and_basis_to_libint_c
!
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
   subroutine initialize_eri_c(eri_precision) bind(C, name='initialize_eri')
!
      use iso_c_binding
      implicit none
      real(c_double), value :: eri_precision
!
   end subroutine initialize_eri_c
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
   subroutine initialize_coulomb_external_charges_c(charges,coordinates,n_points) &
            bind(C, name='initialize_coulomb_external_charges')
!
      use iso_c_binding
      implicit none
      real(c_double), dimension(*), intent(in) :: charges
      real(c_double), dimension(*), intent(in) :: coordinates
      integer(c_int), value  :: n_points
!
   end subroutine initialize_coulomb_external_charges_c
!
   subroutine initialize_coulomb_external_unit_charges_c(coordinates,n_points) &
            bind(C, name='initialize_coulomb_external_unit_charges')
!
      use iso_c_binding
      implicit none
      real(c_double), dimension(*), intent(in) :: coordinates
      integer(c_int), value :: n_points
!
   end subroutine initialize_coulomb_external_unit_charges_c
!
   subroutine set_eri_precision_c(prec) bind(C, name='set_eri_precision')
!
      use iso_c_binding
      implicit none
!
      real(c_double), value :: prec
!
   end subroutine set_eri_precision_c
!
!
end interface
