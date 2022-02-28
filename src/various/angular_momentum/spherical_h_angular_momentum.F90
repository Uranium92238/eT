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
module spherical_h_angular_momentum_class
!
!!
!!    Spherical h angular momentum class
!!    Written by Alexander C. Paul, June 2021
!!
!!    Stores information about spherical angular momentum
!!
!!    Libint ordering described in:
!!    https://github.com/evaleev/libint/wiki/using-modern-CPlusPlus-API
!!
!
   use abstract_spherical_angular_momentum_class
!
   type, extends(abstract_spherical_angular_momentum) :: spherical_h_angular_momentum
   end type
!
!
   interface spherical_h_angular_momentum
!
      procedure :: new_spherical_h_angular_momentum
!
   end interface spherical_h_angular_momentum
!
!
contains
!
!
   function new_spherical_h_angular_momentum() result(this)
!!
!!    New spherical h angular momentum
!!    Written by Alexander C. Paul, June 2021
!!
      implicit none
!
      type(spherical_h_angular_momentum) :: this
!
      integer :: i
!
      this%l = 5
!
      this%l_letter = 'h'
!
      this%n_functions = 11
!
      allocate(this%offsets(this%n_functions))
      allocate(this%components(this%n_functions))
      allocate(this%normalization(this%n_functions))
!
!     Molden cannot handle l > 4
      this%offsets = [(i, i=1, 11)]
!
      this%components = ['-5','-4','-3','-2','-1',' 0',' 1',' 2',' 3',' 4', ' 5']
!
      this%normalization = one
!
   end function new_spherical_h_angular_momentum
!
!
end module spherical_h_angular_momentum_class
