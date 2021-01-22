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
   subroutine get_eri_c(g, s1, s2, s3, s4, epsilon, &
               skip, n1, n2, n3, n4) bind(C, name='get_eri')
!
      use iso_c_binding
!
      implicit none
!
      real(c_double) :: g(*)
      real(c_double) :: epsilon
      integer(c_int) :: s1, s2, s3, s4
      integer(c_int) :: skip, n1, n2, n3, n4
!
   end subroutine get_eri_c
!
!
   subroutine set_eri_precision_c(prec) bind(C, name='set_eri_precision')
!
      use iso_c_binding
      implicit none
!
      real(c_double) :: prec
!
   end subroutine set_eri_precision_c
!
!
   subroutine get_eri_1der_c(g_wxyzqk, s1, s2, s3, s4) bind(C, name='get_eri_1der')
!
      use iso_c_binding
!
      implicit none
!
      real(c_double) :: g_wxyzqk(*)
      integer(c_int) :: s1, s2, s3, s4
!
   end subroutine get_eri_1der_c
!
!
end interface
