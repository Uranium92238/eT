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
module cartesian_f_angular_momentum_class
!
!!
!!    Cartesian f angular momentum class
!!    Written by Alexander C. Paul, June 2021
!!
!!    Stores information about cartesian angular momentum
!!
!!    Libint ordering and normalization described in:
!!    https://github.com/evaleev/libint/wiki/using-modern-CPlusPlus-API
!!
!
   use abstract_cartesian_angular_momentum_class
!
!
   type, extends(abstract_cartesian_angular_momentum) :: cartesian_f_angular_momentum
   end type
!
   interface cartesian_f_angular_momentum
!
      procedure :: new_cartesian_f_angular_momentum
!
   end interface cartesian_f_angular_momentum
!
!
contains
!
!
   function new_cartesian_f_angular_momentum() result(this)
!!
!!    New cartesian f angular momentum
!!    Written by Alexander C. Paul, June 2021
!!
      implicit none
!
      type(cartesian_f_angular_momentum) :: this
!
      this%l = 3
!
      this%l_letter = 'f'
!
      this%n_functions = 10
!
      allocate(this%offsets(this%n_functions))
      allocate(this%components(this%n_functions))
      allocate(this%normalization(this%n_functions))
!
!     instead of   xxx, xxy, xxz, xyy, xyz, xzz, yyy, yyz, yzz, zzz
!     Molden needs xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz
      this%offsets = [1,7,10,4,2,3,6,9,8,5]
!
      this%components = ['xxx','xxy','xxz','xyy','xyz','xzz','yyy','yyz','yzz','zzz']
!
      this%normalization = [one, inv_sqrt_5, inv_sqrt_5, inv_sqrt_5, &
                            inv_sqrt_5*inv_sqrt_3, inv_sqrt_5, &
                            one, inv_sqrt_5, inv_sqrt_5, one]
!
   end function new_cartesian_f_angular_momentum
!
!
end module cartesian_f_angular_momentum_class
