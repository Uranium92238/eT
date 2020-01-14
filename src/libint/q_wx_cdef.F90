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
   subroutine construct_ao_q_wx_c(q_xx, q_xy, q_xz, q_yy, q_yz, q_zz, s1, s2) bind(C, name='construct_ao_q_wx')
!
      use iso_c_binding
!
      implicit none
!
      real(c_double) :: q_xx(*), q_xy(*), q_xz(*), q_yy(*), q_yz(*), q_zz(*)
      integer(c_int) :: s1, s2
!
   end subroutine construct_ao_q_wx_c
!
end interface
