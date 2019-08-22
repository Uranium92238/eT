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
!
   subroutine construct_ao_g_wxyz_c(g, s1, s2, s3, s4) bind(C, name='construct_ao_g_wxyz')
!
      use iso_c_binding
!
      implicit none
!
      real(c_double) :: g(*)
      integer(c_int) :: s1, s2, s3, s4
!
   end subroutine construct_ao_g_wxyz_c
!
!
   subroutine construct_ao_g_wxyz_1der_c(g_wxyzqk, s1, s2, s3, s4) bind(C, name='construct_ao_g_wxyz_1der')
!
      use iso_c_binding
!
      implicit none
!
      real(c_double) :: g_wxyzqk(*)
      integer(c_int) :: s1, s2, s3, s4
!
   end subroutine construct_ao_g_wxyz_1der_c
!
!
   subroutine construct_ao_g_wxyz_epsilon_c(g, s1, s2, s3, s4, epsilon, &
               thread, skip, n1, n2, n3, n4) bind(C, name='construct_ao_g_wxyz_epsilon')
!
      use iso_c_binding
!
      implicit none
!
      real(c_double) :: g(*)
      real(c_double) :: epsilon
      integer(c_int) :: s1, s2, s3, s4, thread
      integer(c_int) :: skip, n1, n2, n3, n4
!
   end subroutine construct_ao_g_wxyz_epsilon_c
!
!
end interface
