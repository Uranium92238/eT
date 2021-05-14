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
module angular_momentum
!
!!
!!    angular momentum module
!!    Written by Alexander C. Paul, Dec 2019
!!
!!    Defines labels for angular momentum functions 
!!    as defined by the default ordering in libint2
!!
!!    Libint ordering described in Appendix B of 
!!    https://doi.org/10.1002/jcc.20815
!!
!
   implicit none
!
   character(len=1), parameter :: s = 's'
!
   character(len=1), dimension(3), parameter :: p_cart = ['x', 'y', 'z']
   character(len=2), dimension(3), parameter :: p = [' 1', ' 0', '-1']
!
   character(len=2), dimension(6), parameter :: d_cart = ['xx','xy','xz',&
                                                          'yy','yz','zz']
   character(len=2), dimension(5), parameter :: d = [' 2',' 1',' 0','-1','-2']
!
   character(len=3), dimension(10), parameter :: f_cart = ['xxx','xxy','xxz','xyy','xyz',&
                                                           'xzz','yyy','yyz','yzz','zzz']
   character(len=2), dimension(7), parameter :: f = [' 3',' 2',' 1',' 0','-1','-2','-3']
!
   character(len=4), dimension(15), parameter :: g_cart = ['xxxx','xxxy','xxxz','xxyy','xxyz',&
                                                           'xxzz','xyyy','xyyz','xyzz','xzzz',&
                                                           'yyyy','yyyz','yyzz','yzzz','zzzz']
   character(len=2), dimension(9), parameter  :: g = [' 4',' 3',' 2',' 1',' 0','-1','-2','-3','-4']
!
!
!  To write molden files the AOs have to be reordered
!  within a shell the n-th ao has to have the number (first + offset(n) - 1)
!
!                           instead of   xx, xy, xz, yy, yz, zz
!                           Molden wants xx, yy, zz, xy, xz, yz
   integer, dimension(6) :: d_offsets_cart = [1,4,6,2,3,5]
!
!                            instead of   xxx, xxy, xxz, xyy, xyz, xzz, yyy, yyz, yzz, zzz
!                            Molden wants xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz
   integer, dimension(10) :: f_offsets_cart = [1,7,10,4,2,3,6,9,8,5]
!
!                            instead of   xxxx, xxxy, xxxz, xxyy, xxyz, xxzz, xyyy, xyyz, xyzz, xzzz,
!                                         yyyy, yyyz, yyzz, yzzz, zzzz
!                            Molden wants xxxx, yyyy, zzzz, xxxy, xxxz, yyyx, yyyz, zzzx, zzzy, xxyy,
!                                         xxzz, yyzz, xxyz, yyxz, zzxy
   integer, dimension(15) :: g_offsets_cart = [1,11,15,2,3,7,12,10,14,4,6,13,5,8,9]
!
!                           instead of    2, 1, 0,-1,-2
!                           Molden wants  0, 1,-1, 2,-2
   integer, dimension(5) :: d_offsets = [3,2,4,1,5]
!
!                           instead of    3, 2, 1, 0,-1,-2,-3
!                           Molden wants  0, 1,-1, 2,-2, 3,-3
   integer, dimension(7) :: f_offsets = [4,3,5,2,6,1,7]
!
!                           instead of    4, 3, 2, 1, 0,-1,-2,-3,-4
!                           Molden wants  0, 1,-1, 2,-2, 3,-3, 4,-4
   integer, dimension(9) :: g_offsets = [5,4,6,3,7,2,8,1,9]
!
end module angular_momentum
