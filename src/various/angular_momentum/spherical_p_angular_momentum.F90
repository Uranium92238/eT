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
module spherical_p_angular_momentum_class
!
!!
!!    Spherical p angular momentum class
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
   type, extends(abstract_spherical_angular_momentum) :: spherical_p_angular_momentum
!
      contains
!
      procedure :: get_angular_part &
                => get_angular_part_spherical_p_angular_momentum
!
   end type
!
!
   interface spherical_p_angular_momentum
!
      procedure :: new_spherical_p_angular_momentum
!
   end interface spherical_p_angular_momentum
!
!
contains
!
!
   function new_spherical_p_angular_momentum() result(this)
!!
!!    New spherical p angular momentum
!!    Written by Alexander C. Paul, June 2021
!!
      implicit none
!
      type(spherical_p_angular_momentum) :: this
!
      this%l = 1
!
      this%l_letter = 'p'
!
      this%n_functions = 3
!
      allocate(this%offsets(this%n_functions))
      allocate(this%components(this%n_functions))
      allocate(this%normalization(this%n_functions))
!
!     instead of  -1, 0, 1
!     Molden needs 1, -1, 0 (which is identical to x, y, z)
      this%offsets = [3,1,2]
!
      this%components = ['-1',' 0',' 1']
      ! this%components = [' y',' z',' x']
!
      this%normalization = one
!
   end function new_spherical_p_angular_momentum
!
!
   function get_angular_part_spherical_p_angular_momentum(this, x, y ,z) &
                                                   result(angular_part)
!!
!!    Get angular part
!!    Written by Alexander C. Paul, June 2021
!!
!!    Return angular part of the AOs evaluated at x, y, z
!!    See MEST, p 211, Table 6.3
!!
      implicit none
!
      class(spherical_p_angular_momentum), intent(in) :: this
!
      real(dp), intent(in) :: x, y, z
!
      real(dp), dimension(this%n_functions) :: angular_part
!
!     ml = -1, angular part: y
      angular_part(1) = y
!
!     ml = 0, angular part: z
      angular_part(2) = z
!
!     ml = 1, angular part: x
      angular_part(3) = x
!
   end function get_angular_part_spherical_p_angular_momentum
!
!
end module spherical_p_angular_momentum_class
