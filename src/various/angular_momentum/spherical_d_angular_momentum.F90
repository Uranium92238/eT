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
module spherical_d_angular_momentum_class
!
!!
!!    Spherical d angular momentum class
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
   type, extends(abstract_spherical_angular_momentum) :: spherical_d_angular_momentum
!
      contains
!
      procedure :: get_angular_part &
                => get_angular_part_spherical_d_angular_momentum
!
   end type
!
!
   interface spherical_d_angular_momentum
!
      procedure :: new_spherical_d_angular_momentum
!
   end interface spherical_d_angular_momentum
!
!
contains
!
!
   function new_spherical_d_angular_momentum() result(this)
!!
!!    New spherical d angular momentum
!!    Written by Alexander C. Paul, June 2021
!!
      implicit none
!
      type(spherical_d_angular_momentum) :: this
!
      this%l = 2
!
      this%l_letter = 'd'
!
      this%n_functions = 5
!
      allocate(this%offsets(this%n_functions))
      allocate(this%components(this%n_functions))
      allocate(this%normalization(this%n_functions))
!
!     instead of  -2,-1, 0, 1, 2
!     Molden needs 0, 1,-1, 2,-2
      this%offsets = [3,4,2,5,1]
!
      this%components = ['-2','-1',' 0',' 1',' 2']
!
      this%normalization = one
!
   end function new_spherical_d_angular_momentum
!
!
   function get_angular_part_spherical_d_angular_momentum(this, x, y ,z) &
                                                   result(angular_part)
!!
!!    Get angular part
!!    Written by Alexander C. Paul, June 2021
!!
!!    Return angular part of the AOs evaluated at x, y, z
!!
      implicit none
!
      class(spherical_d_angular_momentum), intent(in) :: this
!
      real(dp), intent(in) :: x, y, z
!
      real(dp), dimension(this%n_functions) :: angular_part
!
      real(dp) :: r_squared
!
      r_squared = x**2 + y**2 + z**2
!
!     ml = -2, angular part: sqrt(3)xy
      angular_part(1)   = sqrt_3 * x*y
!
!     ml = -1, angular part: sqrt(3)yz
      angular_part(2) = sqrt_3 * y*z
!
!     ml = 0, angular part: 1.5z^2 - 0.5r^2
      angular_part(3) = half * (three*z**2 - r_squared)
!
!     ml = 1, angular part: sqrt(3)xz
      angular_part(4) = sqrt_3 * x*z
!
!     ml = 2, angular part: 1/2 sqrt(3)(x^2 - y^2)
      angular_part(5) = half*sqrt_3 * (x**2 - y**2)
!
end function get_angular_part_spherical_d_angular_momentum
!
!
end module spherical_d_angular_momentum_class
