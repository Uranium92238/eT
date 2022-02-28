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
module spherical_f_angular_momentum_class
!
!!
!!    Spherical f angular momentum class
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
   type, extends(abstract_spherical_angular_momentum) :: spherical_f_angular_momentum
!
      contains
!
      procedure :: get_angular_part &
                => get_angular_part_spherical_f_angular_momentum
!
   end type
!
!
   interface spherical_f_angular_momentum
!
      procedure :: new_spherical_f_angular_momentum
!
   end interface spherical_f_angular_momentum
!
!
contains
!
!
   function new_spherical_f_angular_momentum() result(this)
!!
!!    New spherical f angular momentum
!!    Written by Alexander C. Paul, June 2021
!!
      implicit none
!
      type(spherical_f_angular_momentum) :: this
!
      this%l = 3
!
      this%l_letter = 'f'
!
      this%n_functions = 7
!
      allocate(this%offsets(this%n_functions))
      allocate(this%components(this%n_functions))
      allocate(this%normalization(this%n_functions))
!
!     instead of  -3,-2,-1, 0, 1, 2, 3
!     Molden needs 0, 1,-1, 2,-2, 3,-3
      this%offsets = [4,5,3,6,2,7,1]
!
      this%components = ['-3','-2','-1',' 0',' 1',' 2',' 3']
!
      this%normalization = one
!
   end function new_spherical_f_angular_momentum
!
!
   function get_angular_part_spherical_f_angular_momentum(this, x, y ,z) &
                                                   result(angular_part)
!!
!!    Get angular part
!!    Written by Alexander C. Paul, June 2021
!!
!!    Return angular part of the AOs evaluated at x, y, z
!!
      implicit none
!
      class(spherical_f_angular_momentum), intent(in) :: this
!
      real(dp), intent(in) :: x, y, z
!
      real(dp), dimension(this%n_functions) :: angular_part
!
      real(dp) :: r_squared
!
      real(dp), parameter :: sqrt_3_half = sqrt(three*half)
      real(dp), parameter :: sqrt_15     = sqrt(15.0d0)
      real(dp), parameter :: sqrt_5_half = sqrt(five*half)
!
      r_squared = x**2 + y**2 + z**2
!
!     ml = -3, angular part: 1/2 sqrt(5/2) (3x^2 - y^2)y
      angular_part(1) = half*sqrt_5_half * (three*x**2 - y**2)*y
!
!     ml = -2, angular part: sqrt(15) xyz
      angular_part(2) = sqrt_15 * x*y*z
!
!     ml = -1, angular part: 1/2 sqrt(3/2) (5*z**2 - r**2)y
!     ml =  1, angular part: 1/2 sqrt(3/2) (5*z**2 - r**2)x
      angular_part(3) = half*sqrt_3_half * (five*z**2 - r_squared)
      angular_part(5) = angular_part(3) * x
      angular_part(3) = angular_part(3) * y
!
!     ml = 0, angular part: 1/2(5z**2 - 3r^2)z
      angular_part(4) = half * (five*z**2 - three*r_squared)*z
!
!     ml = 2, angular part: 1/2 sqrt(15)(x^2-y^2)z
      angular_part(6) = half*sqrt_15 * (x**2 - y**2)*z
!
!     ml = 3, angular part: 1/2 sqrt(5/2)(x^2-3y^2)x
      angular_part(7) = half*sqrt_5_half * (x**2 - three*y**2)*x
!
   end function get_angular_part_spherical_f_angular_momentum
!
!
end module spherical_f_angular_momentum_class
