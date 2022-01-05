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
module cartesian_d_angular_momentum_class
!
!!
!!    Cartesian d angular momentum class
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
   type, extends(abstract_angular_momentum) :: cartesian_d_angular_momentum
   end type
!
   interface cartesian_d_angular_momentum
!
      procedure :: new_cartesian_d_angular_momentum
!
   end interface cartesian_d_angular_momentum
!
!
contains
!
!
   function new_cartesian_d_angular_momentum() result(this)
!!
!!    New cartesian d angular momentum
!!    Written by Alexander C. Paul, June 2021
!!
      implicit none
!
      type(cartesian_d_angular_momentum) :: this
!
      this%l = 2
!
      this%l_letter = 'd'
!
      allocate(this%offsets(6))
      allocate(this%components(6))
      allocate(this%normalization(6))
!
!     instead of   xx, xy, xz, yy, yz, zz
!     Molden needs xx, yy, zz, xy, xz, yz
      this%offsets = [1,4,6,2,3,5]
!
      this%components = ['xx','xy','xz','yy','yz','zz']
!
      this%normalization = [one, inv_sqrt_3, inv_sqrt_3, one, inv_sqrt_3, one]
!
   end function new_cartesian_d_angular_momentum
!
!
end module cartesian_d_angular_momentum_class
