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
module abstract_angular_momentum_class
!
!!
!!    Abstract angular momenta class
!!    Written by Alexander C. Paul, June 2021
!!
!!    Stores and returns information about spherical angular momenta
!!
!!    Libint ordering described in:
!!    https://github.com/evaleev/libint/wiki/using-modern-CPlusPlus-API
!!
!
   use parameters
!
   type, abstract :: abstract_angular_momentum
!
      integer :: l
      character(len=1) :: l_letter
      character(len=6), dimension(:), allocatable :: components
!
      integer, dimension(:), allocatable :: offsets
!
      real(dp), dimension(:), allocatable :: normalization
!
   contains
!
      procedure, non_overridable :: get_label
      procedure, non_overridable :: get_molden_offset
      procedure, non_overridable :: get_normalization_factor
      procedure, non_overridable :: get_l_integer
      procedure, non_overridable :: get_l_string
!
   end type
!
!
contains
!
!
   pure function get_l_integer(this) result(l)
!!
!!    Get l integer
!!    Written by Alexander C. Paul, June 2021
!!
!!    Return angular momentum as integer
!!
      implicit none
!
      class(abstract_angular_momentum), intent(in) :: this
!
      integer :: l
!
      l = this%l
!
   end function get_l_integer
!
!
   pure function get_l_string(this) result(l)
!!
!!    Get l string
!!    Written by Alexander C. Paul, June 2021
!!
!!    Return angular momentum as letter
!!
      implicit none
!
      class(abstract_angular_momentum), intent(in) :: this
!
      character(len=1) :: l
!
      l = this%l_letter
!
   end function get_l_string
!
!
   pure function get_label(this, index_) result(label)
!!
!!    Get label
!!    Written by Alexander C. Paul, June 2021
!!
!!    Returns string containing the angular momenta and its spatial component
!!
      implicit none
!
      class(abstract_angular_momentum), intent(in) :: this
!
      integer, intent(in) :: index_
!
      character(len=:), allocatable :: label
!
      label = this%l_letter // ' ' // trim(this%components(index_))
!
   end function get_label
!
!
   pure function get_molden_offset(this, index_) result(offset)
!!
!!    Get Molden offset
!!    Written by Alexander C. Paul, June 2021
!!
!!    Return offset from the first AO index_ in the shell
!!    to obtain the AO order molden expects.
!!
      implicit none
!
      class(abstract_angular_momentum), intent(in) :: this
!
      integer, intent(in)  ::index_
!
      integer :: offset
!
      offset = this%offsets(index_)
!
   end function get_molden_offset
!
!
   pure function get_normalization_factor(this, index_) result(factor)
!!
!!    Get normalization factor
!!    Written by Alexander C. Paul, June 2021
!!
!!    Return factor to normalize the function determined
!!    by the angular momentum l and the index_
!!
      implicit none
!
      class(abstract_angular_momentum), intent(in) :: this
!
      integer, intent(in) :: index_
!
      real(dp) :: factor
!
      factor = this%normalization(index_)
!
   end function get_normalization_factor
!
!
end module abstract_angular_momentum_class
