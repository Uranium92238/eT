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
module named_range_class
!!
!!    Named range class
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2020
!!
!!    Extends the range class with a name
!!
   use range_class
!
   implicit none
!
   type, extends(range_) :: named_range
!
      character(len=200), private :: name_ = 'no_name'
!
   contains
!
      procedure, public :: set_range => set_range_named_range
      procedure, public :: get_name  => get_name_named_range
!
   end type named_range
!
   interface named_range
!
      procedure :: new_named_range
      procedure :: new_uninitialized_named_range
      procedure :: copy_named_range
!
   end interface named_range
!
contains
!
!
   pure function new_named_range(name_, first, length, step) result(this)
!!
!!    New named range
!!    Written by Eirik F. Kjønstad, 2020
!!
      implicit none
!
      character(len=*), intent(in) :: name_
!
      integer, intent(in) :: first, length
      integer, optional, intent(in) :: step
!
      type(named_range) :: this
!
      this%range_ = range_(first, length, step)
!
      this%name_ = trim(name_)
!
   end function new_named_range
!
!
   pure function new_uninitialized_named_range(name_) result(this)
!!
!!    New uninitialized range
!!    Written by Rolf H. Myhre, Jun 2021
!!
      implicit none
!
      character(len=*), intent(in) :: name_
!
      type(named_range) :: this
!
      this%name_ = trim(name_)
!
   end function new_uninitialized_named_range
!
!
   pure function copy_named_range(that) result(this)
!!
!!    Copy named range
!!    Written by Eirik F. Kjønstad, 2020
!!
      implicit none
!
      class(named_range), intent(in) :: that
!
      type(named_range) :: this
!
      this%range_ = range_(that)
!
      this%name_ = trim(that%name_)
!
   end function copy_named_range
!
!
   pure subroutine set_range_named_range(this, first, length, step)
!!
!!    Set range
!!    Written by Rolf H. Myhre, Apr 2021
!!
      implicit none
!
      class(named_range), intent(inout) :: this
      integer, intent(in)               :: first, length
      integer, optional, intent(in)     :: step
!
      this%range_ = range_(first, length, step)
!
   end subroutine set_range_named_range
!
!
   pure function get_name_named_range(this) result(name_)
!!
!!    Get name
!!    Written by Eirik F. Kjønstad, 2020
!!
      implicit none
!
      class(named_range), intent(in) :: this
!
      character(len=len_trim(this%name_)) :: name_
!
      name_ = trim(this%name_)
!
   end function get_name_named_range
!
!
end module named_range_class
