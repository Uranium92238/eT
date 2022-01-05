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
module s_angular_momentum_class
!
!!
!!    s angular momentum class
!!    Written by Alexander C. Paul, June 2021
!!
!!    Stores information about s angular momentum
!!
!!    Libint ordering described in:
!!    https://github.com/evaleev/libint/wiki/using-modern-CPlusPlus-API
!!
!
   use abstract_angular_momentum_class
!
   type, extends(abstract_angular_momentum) :: s_angular_momentum
   end type
!
   interface s_angular_momentum
!
      procedure :: new_s_angular_momentum
!
   end interface s_angular_momentum
!
!
contains
!
!
   function new_s_angular_momentum() result(this)
!!
!!    New s angular momentum
!!    Written by Alexander C. Paul, June 2021
!!
      implicit none
!
      type(s_angular_momentum) :: this
!
      this%l = 0
      this%l_letter = 's'
!
      allocate(this%offsets(1))
      allocate(this%components(1))
      allocate(this%normalization(1))
!
      this%offsets = 1
      this%components = ''
      this%normalization = one
!
   end function new_s_angular_momentum
!
!
end module s_angular_momentum_class
