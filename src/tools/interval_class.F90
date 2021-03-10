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
module interval_class
!!
!!    Interval class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
   implicit none
!
   type interval
!
      integer :: first  = -1
      integer :: last   = -1
      integer :: length = 0
!
   contains
!
      procedure :: is_subset              => is_subset_interval
      procedure :: empty_intersection     => empty_intersection_interval 
      procedure :: element_is_member      => element_is_member_interval
      procedure :: set_interval           => set_interval_interval
!
   end type interval
!
   interface interval
!
      procedure :: new_interval 
      procedure :: new_interval_from_template 
!
   end interface interval 
!
contains 
!
!
   pure function new_interval(first, last) result(this)
!!
!!    New interval 
!!    Written by Eirik F. Kjønstad, 2019 
!!
      implicit none 
!
      integer, intent(in) :: first, last 
!
      type(interval) :: this 
!
      call this%set_interval(first, last)
!
   end function new_interval
!
!
   pure function new_interval_from_template(that) result(this) 
!!
!!    New interval from template 
!!    Written by Eirik F. Kjønstad, 2020
!!
      implicit none 
!
      class(interval), intent(in) :: that
!
      type(interval) :: this 
!
      call this%set_interval(that%first, that%last)
!
   end function new_interval_from_template
!
!
   pure function is_subset_interval(this, first, last) result(is_subset)
!!
!!    Is subset 
!!    Written by Eirik F. Kjønstad, Apr 2020 
!!
!!    Checks if [first, last] is a subset of the interval.
!!
      implicit none 
!
      class(interval), intent(in) :: this 
!
      integer, intent(in) :: first, last 
!
      logical :: is_subset
!
      is_subset = .false.
      if (first .ge. this%first .and. last .le. this%last) is_subset = .true.
!
   end function is_subset_interval
!
!
   pure function empty_intersection_interval(this, first, last) result(empty)
!!
!!    Empty intersection 
!!    Written by Eirik F. Kjønstad, Apr 2020 
!!
!!    Does [first, last] have no elements in common with interval?
!!
      implicit none 
!
      class(interval), intent(in) :: this 
!
      integer, intent(in) :: first, last 
!
      logical :: empty 
!
      empty = .true.
!
!     Is first or last in the interval? 
!
      if (first .ge. this%first .and. first .le. this%last .or. &
          last  .ge. this%first .and. last  .le. this%last) then 
!
         empty = .false.
         return 
!
      endif 
!
!     Is interval contained in [first, last]? 
!
      if (first .lt. this%first .and. last .gt. this%last) empty = .false.
!
   end function empty_intersection_interval
!
!
   pure function element_is_member_interval(this, n) result(member)
!!
!!    Element is member 
!!    Written by Eirik F. Kjønstad, Mar 2020
!!
!!    Returns true if n is in the interval.
!!
      implicit none 
!
      class(interval), intent(in) :: this 
!
      integer, intent(in) :: n  
!
      logical :: member
!
      member = .false.
      if (n .ge. this%first .and. n .le. this%last) member = .true.
!
   end function element_is_member_interval
!
!
   pure subroutine set_interval_interval(this, first, last)
!!
!!    Set interval
!!    Written by Eirik F. Kjønstad, Jan 2021
!!
!!    Sets the interval (first, last, length) from specified first and last.
!!
      implicit none 
!
      class(interval), intent(inout) :: this 
!
      integer, intent(in) :: first, last 
!
      this%first  = first 
      this%last   = last 
      this%length = last - first + 1 
!
   end subroutine set_interval_interval
!
!
end module interval_class
