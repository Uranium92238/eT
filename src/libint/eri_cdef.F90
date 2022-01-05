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
   subroutine get_eri_c(g, s1, s2, s3, s4, &
                        epsilon_, skip,    &
                        n1, n2, n3, n4) bind(C, name='get_eri')
!
      use iso_c_binding
!
      implicit none
!
      real(c_double), dimension(*), intent(out) :: g
      integer(c_int), value :: s1, s2, s3, s4
      real(c_double), value :: epsilon_
      integer(c_int), intent(out) :: skip
      integer(c_int), value :: n1, n2, n3, n4
!
   end subroutine get_eri_c
!
!
   subroutine get_eri_3c_c(g, J, s3, s4, &
                        epsilon_, skip,    &
                        nJ, n3, n4) bind(C, name='get_eri_3c')
!
      use iso_c_binding
!
      implicit none
!
      real(c_double), dimension(*), intent(out) :: g
      integer(c_int), value :: J, s3, s4
      real(c_double), value :: epsilon_
      integer(c_int), intent(out) :: skip
      integer(c_int), value :: nJ, n3, n4
!
   end subroutine get_eri_3c_c
!
!
   subroutine get_eri_2c_c(g, J, K, &
                        epsilon_, skip,    &
                        nJ, nK) bind(C, name='get_eri_2c')
!
      use iso_c_binding
!
      implicit none
!
      real(c_double), dimension(*), intent(out) :: g
      integer(c_int), value :: J, K
      real(c_double), value :: epsilon_
      integer(c_int), intent(out) :: skip
      integer(c_int), value :: nJ, nK
!
   end subroutine get_eri_2c_c
!
!
   subroutine get_eri_1der_c(g_wxyzqk, s1, s2, s3, s4) bind(C, name='get_eri_1der')
!
      use iso_c_binding
!
      implicit none
!
      real(c_double), dimension(*), intent(out) :: g_wxyzqk
      integer(c_int), value :: s1, s2, s3, s4
!
   end subroutine get_eri_1der_c
!
!
end interface
