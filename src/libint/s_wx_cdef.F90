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
interface
!
!
   subroutine construct_ao_s_wx_c(s, s1, s2) bind(C, name='construct_ao_s_wx')
!
      use iso_c_binding
!
      implicit none
!
      real(c_double) :: s(*)
      integer(c_int) :: s1, s2
!
   end subroutine construct_ao_s_wx_c
!
!
   subroutine construct_ao_s_wx_1der_c(s_1x, s_1y, s_1z, s_2x, s_2y, s_2z, s1, s2) bind(C, name='construct_ao_s_wx_1der')
!
      use iso_c_binding
!
      implicit none
!
      real(c_double) :: s_1x(*), s_1y(*), s_1z(*), s_2x(*), s_2y(*), s_2z(*)
      integer(c_int) :: s1, s2
!
   end subroutine construct_ao_s_wx_1der_c
!
!
end interface
