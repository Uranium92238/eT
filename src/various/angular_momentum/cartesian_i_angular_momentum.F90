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
module cartesian_i_angular_momentum_class
!
!!
!!    Cartesian i angular momentum class
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
   type, extends(abstract_angular_momentum) :: cartesian_i_angular_momentum
   end type
!
   interface cartesian_i_angular_momentum
!
      procedure :: new_cartesian_i_angular_momentum
!
   end interface cartesian_i_angular_momentum
!
!
contains
!
!
   function new_cartesian_i_angular_momentum() result(this)
!!
!!    New cartesian i angular momentum
!!    Written by Alexander C. Paul, June 2021
!!
      implicit none
!
      type(cartesian_i_angular_momentum) :: this
!
      integer :: i
!
      this%l = 6
!
      this%l_letter = 'i'
!
      allocate(this%offsets(28))
      allocate(this%components(28))
      allocate(this%normalization(28))
!
!     Molden cannot handle l > 4
      this%offsets = [(i, i=1, 28)]
!
      this%components = ['xxxxxx','xxxxxy','xxxxxz','xxxxyy','xxxxyz','xxxxzz','xxxyyy', &
                         'xxxyyz','xxxyzz','xxxzzz','xxyyyy','xxyyyz','xxyyzz','xxyzzz', &
                         'xxzzzz','xyyyyy','xyyyyz','xyyyzz','xyyzzz','xyzzzz','xzzzzz', &
                         'yyyyyy','yyyyyz','yyyyzz','yyyzzz','yyzzzz','yzzzzz','zzzzzz']
!
      this%normalization = [one, &
                            inv_sqrt_11, &
                            inv_sqrt_11, &
                            inv_sqrt_3*inv_sqrt_11, &
                            third*inv_sqrt_11, &
                            inv_sqrt_3*inv_sqrt_11, &
                            inv_sqrt_3*inv_sqrt_7*inv_sqrt_11*sqrt_5, &
                            inv_sqrt_3*inv_sqrt_7*inv_sqrt_11, &
                            inv_sqrt_3*inv_sqrt_7*inv_sqrt_11, &
                            inv_sqrt_3*inv_sqrt_7*inv_sqrt_11*sqrt_5, &
                            inv_sqrt_3*inv_sqrt_11, &
                            inv_sqrt_3*inv_sqrt_7*inv_sqrt_11, &
                            inv_sqrt_5*inv_sqrt_7*inv_sqrt_11, &
                            inv_sqrt_3*inv_sqrt_7*inv_sqrt_11, &
                            inv_sqrt_3*inv_sqrt_11, &
                            inv_sqrt_11, &
                            third*inv_sqrt_11, &
                            inv_sqrt_3*inv_sqrt_7*inv_sqrt_11, &
                            inv_sqrt_3*inv_sqrt_7*inv_sqrt_11, &
                            third*inv_sqrt_11, &
                            inv_sqrt_11, &
                            one, &
                            inv_sqrt_11, &
                            inv_sqrt_3*inv_sqrt_11, &
                            inv_sqrt_3*inv_sqrt_7*inv_sqrt_11*sqrt_5, &
                            inv_sqrt_3*inv_sqrt_11, &
                            inv_sqrt_11, &
                            one]
!
   end function new_cartesian_i_angular_momentum
!
!
end module cartesian_i_angular_momentum_class
