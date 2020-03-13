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
module angular_momentum
!
!!
!!    angular momentum module
!!    Written by Alexander C. Paul, Dec 2019
!!
!!    Defines labels for angular momentum functions
!!
!!       In declarations:   "real(dp) :: foo", "integer(i64) :: foo_int", etc.
!!       In record-lengths: "recl=dp*200" (200 double precision numbers per record)
!!
!
   implicit none
!
   character(len=1), parameter :: s = 's'
!
   character(len=3), dimension(3), parameter :: p_cart = ['p x', 'p y', 'p z']
   character(len=4), dimension(3), parameter :: p = ['p  1', 'p  0', 'p -1']
!
   character(len=4), dimension(6), parameter :: d_cart = ['d xx','d xy','d xz',&
                                                          'd yy','d yz','d zz']
   character(len=4), dimension(5), parameter :: d = ['d  2','d  1','d  0','d -1','d -2']
!
   character(len=5), dimension(10), parameter :: f_cart = ['f xxx','f xxy','f xxz','f xyy','f xyz',&
                                                           'f xzz','f yyy','f yyz','f yzz','f zzz']
   character(len=4), dimension(7), parameter :: f = ['f  3','f  2','f  1','f  0','f -1','f -2','f -3']
!
   character(len=1), parameter :: g = 'g'
!
   character(len=1), parameter :: h = 'h'
!
end module angular_momentum
