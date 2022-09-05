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
module abstract_cartesian_angular_momentum_class
!
!!
!!    Abstract cartesian angular momentum class
!!    Written by Alexander C. Paul, Feb 2022
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
   type, abstract, extends(abstract_angular_momentum) :: abstract_cartesian_angular_momentum
!
   contains
!
      procedure :: get_angular_part &
                => get_angular_part_abstract_cartesian_angular_momentum
!
   end type
!
!
contains
!
!
   function get_angular_part_abstract_cartesian_angular_momentum(this, x, y ,z) &
                                                                 result(angular_part)
!!
!!    Get angular part
!!    Written by Alexander C. Paul, June 2021
!!
!!    Return angular part of the AOs evaluated at x, y, z
!!
      use array_initialization, only: constant_array
!
      implicit none
!
      class(abstract_cartesian_angular_momentum), intent(in) :: this
!
      real(dp), intent(in) :: x, y, z
!
      real(dp), dimension(this%n_functions) :: angular_part
!
      integer :: i, j
      character(len=6) :: component
!
      call constant_array(angular_part, this%n_functions, one)
!
      do i = 1, this%n_functions
!
         component = this%components(i)
!
         do j = 1, 6
!
            select case (component(j:j))
               case ('x')
                  angular_part(i) = angular_part(i) * x
               case ('y')
                  angular_part(i) = angular_part(i) * y
               case ('z')
                  angular_part(i) = angular_part(i) * z
            end select
!
         end do
      end do
!
   end function get_angular_part_abstract_cartesian_angular_momentum
!
!
end module abstract_cartesian_angular_momentum_class
