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
module spherical_i_angular_momentum_class
!
!!
!!    Spherical i angular momentum class
!!    Written by Alexander C. Paul, June 2021
!!
!!    Stores information about spherical angular momentum
!!
!!    Libint ordering described in:
!!    https://github.com/evaleev/libint/wiki/using-modern-CPlusPlus-API
!!
!
   use abstract_angular_momentum_class
!
   type, extends(abstract_angular_momentum) :: spherical_i_angular_momentum
   end type
!
!
   interface spherical_i_angular_momentum
!
      procedure :: new_spherical_i_angular_momentum
!
   end interface spherical_i_angular_momentum
!
!
contains
!
!
   function new_spherical_i_angular_momentum() result(this)
!!
!!    New spherical i angular momentum
!!    Written by Alexander C. Paul, June 2021
!!
      implicit none
!
      type(spherical_i_angular_momentum) :: this
!
      integer :: i
!
      this%l = 6
!
      this%l_letter = 'i'
!
      allocate(this%offsets(13))
      allocate(this%components(13))
      allocate(this%normalization(13))
!
!     Molden cannot handle l > 4
      this%offsets = [(i, i=1, 13)]
!
      this%components = ['-6','-5','-4','-3','-2','-1',' 0',' 1',' 2',' 3',' 4', ' 5', ' 6']
!
      this%normalization = one
!
   end function new_spherical_i_angular_momentum
!
!
end module spherical_i_angular_momentum_class
