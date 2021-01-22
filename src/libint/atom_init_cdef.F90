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
   subroutine get_shell_numbers_c(atom, sn) bind(C, name='get_shell_numbers')
!
      use iso_c_binding
!
      implicit none
!
      integer(c_int) :: atom
      integer(c_int), dimension(1,1) :: sn
!
   end subroutine get_shell_numbers_c
!
   subroutine get_first_ao_in_shells_c(atom, faois) bind(C, name='get_first_ao_in_shells')
!
      use iso_c_binding
!
      implicit none
!
      integer(c_int) :: atom
      integer(c_int), dimension(1,1) :: faois
!
   end subroutine get_first_ao_in_shells_c
!
   subroutine get_n_aos_in_shell_c(atom, n_aos_in_shell) bind(C, name='get_n_aos_in_shell')
!
      use iso_c_binding
!
      implicit none
!
      integer(c_int) :: atom
      integer(c_int), dimension(1,1) :: n_aos_in_shell
!
   end subroutine get_n_aos_in_shell_c
!
   subroutine get_n_shells_on_atoms_c(nsoa) bind(C, name='get_n_shells_on_atoms')
!
      use iso_c_binding
!
      implicit none
!
      integer(c_int), dimension(1,1) :: nsoa
!
   end subroutine get_n_shells_on_atoms_c
!
   subroutine get_n_shells_on_atom_c(atom, n_shells) bind(C, name='get_n_shells_on_atom')
!
      use iso_c_binding
!
      implicit none
!
      integer(c_int) :: atom, n_shells
!
   end subroutine get_n_shells_on_atom_c
!
   subroutine initialize_atom_to_shell_list_c() bind(C, name='initialize_atom_to_shell_list')
!
      use iso_c_binding
!
      implicit none
!
   end subroutine initialize_atom_to_shell_list_c
!
   subroutine initialize_shell_to_first_ao_c() bind(C, name='initialize_shell_to_first_ao')
!
      use iso_c_binding
!
      implicit none
!
   end subroutine initialize_shell_to_first_ao_c
!
end interface
