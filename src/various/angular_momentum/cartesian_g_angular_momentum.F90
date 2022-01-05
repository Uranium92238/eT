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
module cartesian_g_angular_momentum_class
!
!!
!!    Cartesian gangular momentum class
!!    Written by Alexander C. Paul, June 2021
!!
!!    Stores information about cartesian angular momentum
!!
!!    Libint ordering and normalization described in:
!!    https://github.com/evaleev/libint/wiki/using-modern-CPlusPlus-API
!!
!
   use abstract_angular_momentum_class
!
!
   type, extends(abstract_angular_momentum) :: cartesian_g_angular_momentum
   end type
!
   interface cartesian_g_angular_momentum
!
      procedure :: new_cartesian_g_angular_momentum
!
   end interface cartesian_g_angular_momentum
!
!
contains
!
!
   function new_cartesian_g_angular_momentum() result(this)
!!
!!    New cartesian g angular momentum
!!    Written by Alexander C. Paul, June 2021
!!
      implicit none
!
      type(cartesian_g_angular_momentum) :: this
!
      this%l = 4
!
      this%l_letter = 'g'
!
      allocate(this%offsets(15))
      allocate(this%components(15))
      allocate(this%normalization(15))
!
!     instead of   xxxx, xxxy, xxxz, xxyy, xxyz, xxzz, xyyy, xyyz,
!                  xyzz, xzzz, yyyy, yyyz, yyzz, yzzz, zzzz
!     Molden needs xxxx, yyyy, zzzz, xxxy, xxxz, yyyx, yyyz, zzzx,
!                  zzzy, xxyy, xxzz, yyzz, xxyz, yyxz, zzxy
      this%offsets = [1,11,15,2,3,7,12,10,14,4,6,13,5,8,9]
!
      this%components = ['xxxx','xxxy','xxxz','xxyy','xxyz','xxzz','xyyy','xyyz', &
                         'xyzz','xzzz','yyyy','yyyz','yyzz','yzzz','zzzz']
!
      this%normalization = [one,  inv_sqrt_7, inv_sqrt_7, &
                            inv_sqrt_7*inv_sqrt_5*sqrt_3, inv_sqrt_7*inv_sqrt_5, &
                            inv_sqrt_7*inv_sqrt_5*sqrt_3, inv_sqrt_7, &
                            inv_sqrt_7*inv_sqrt_5, inv_sqrt_7*inv_sqrt_5, &
                            inv_sqrt_7, one, inv_sqrt_7, &
                            inv_sqrt_7*inv_sqrt_5*sqrt_3, inv_sqrt_7, one]
!
   end function new_cartesian_g_angular_momentum
!
!
end module cartesian_g_angular_momentum_class
