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
   subroutine get_n_ri_shells_c(n_shells) bind(C, name='get_n_ri_shells')
!
      use iso_c_binding
!
      implicit none
!
      integer(c_int), dimension(1) :: n_shells
!
   end subroutine get_n_ri_shells_c
!
   subroutine get_n_ri_ao_c(n_ao) bind(C, name='get_n_ri_ao')
!
      use iso_c_binding
!
      implicit none
!
      integer(c_int), dimension(1) :: n_ao
!
   end subroutine get_n_ri_ao_c
!
   subroutine get_ri_shell_size_c(shell, n_ao) bind(C, name='get_ri_shell_size')
!
      use iso_c_binding
!
      implicit none
!
      integer(c_int), dimension(1) :: n_ao, shell
!
   end subroutine get_ri_shell_size_c
!
!
end interface
