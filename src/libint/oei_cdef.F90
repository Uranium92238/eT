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
interface
!
!
   subroutine get_oei_c(type_, x, s1, s2) bind(C, name='get_oei')
!
      use iso_c_binding
!
      implicit none
!
      character(c_char) :: type_
      real(c_double) :: x(*)
      integer(c_int) :: s1, s2
!
   end subroutine get_oei_c
!
!
   subroutine get_oei_1der_c(type_, x, s1, s2) bind(C, name='get_oei_1der')
!
      use iso_c_binding
!
      implicit none
!
      character(c_char) :: type_
      real(c_double) :: x(*)
      integer(c_int) :: s1, s2
!
   end subroutine get_oei_1der_c
!
!
   subroutine add_nuclear_h_1der_c(h_wxqk, s1, s2, n_ao) &
                           bind(C, name='add_nuclear_h_1der')
!
      use iso_c_binding
!
      implicit none
!
      real(c_double) :: h_wxqk(*)
      integer(c_int) :: s1, s2, n_ao
!
   end subroutine add_nuclear_h_1der_c
!
!
end interface
