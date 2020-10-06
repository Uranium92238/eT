!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2020 the authors of eT
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
   use kinds
!
   implicit none
!
   type interval
!
      integer :: first = -1
      integer :: last  = -1
      integer :: length = -1
!
   contains
!
      procedure :: is_subset              => is_subset_interval
      procedure :: empty_intersection     => empty_intersection_interval 
      procedure :: element_is_member      => element_is_member_interval
!
   end type interval
!
   interface interval
!
      procedure :: new_interval 
      procedure :: new_interval_copy
!
   end interface interval 
!
contains 
!
!
   function new_interval(first, last) result(intval)
!!
!!    New interval 
!!    Written by Eirik F. Kjønstad, 2019 
!!
      implicit none 
!
      integer, intent(in) :: first, last 
!
      type(interval) :: intval 
!
      intval%first = first 
      intval%last  = last 
      intval%length = last - first + 1
!
   end function new_interval
!
!
   function new_interval_copy(intval_copy) result(intval)
!!
!!    New interval 
!!    Written by Eirik F. Kjønstad, 2019 
!!
      implicit none 
!
      class(interval), intent(in) :: intval_copy 
!
      type(interval) :: intval 
!
      intval%first  = intval_copy%first 
      intval%last   = intval_copy%last  
      intval%length = intval_copy%length
!
   end function new_interval_copy
!
!
   pure function is_subset_interval(intval, first, last) result(is_subset)
!!
!!    Is subset 
!!    Written by Eirik F. Kjønstad, Apr 2020 
!!
!!    Checks if [first, last] is a subset of the interval.
!!
      implicit none 
!
      class(interval), intent(in) :: intval 
!
      integer, intent(in) :: first, last 
!
      logical :: is_subset
!
      is_subset = .false.
      if (first .ge. intval%first .and. last .le. intval%last) is_subset = .true.
!
   end function is_subset_interval
!
!
   pure function empty_intersection_interval(intval, first, last) result(empty)
!!
!!    Empty intersection 
!!    Written by Eirik F. Kjønstad, Apr 2020 
!!
!!    Does [first, last] have no elements in common with interval?
!!
      implicit none 
!
      class(interval), intent(in) :: intval 
!
      integer, intent(in) :: first, last 
!
      logical :: empty 
!
      empty = .true.
!
!     Is first or last in the interval? 
!
      if (first .ge. intval%first .and. first .le. intval%last .or. &
          last  .ge. intval%first .and. last .le. intval%last) then 
!
         empty = .false.
         return 
!
      endif 
!
!     Is interval contained in [first, last]? 
!
      if (first .lt. intval%first .and. last .gt. intval%last) empty = .false.
!
   end function empty_intersection_interval
!
!
   pure function element_is_member_interval(intval, n) result(member)
!!
!!    Element is member 
!!    Written by Eirik F. Kjønstad, Mar 2020
!!
!!    Returns true if n is in the interval.
!!
      implicit none 
!
      class(interval), intent(in) :: intval 
!
      integer, intent(in) :: n  
!
      logical :: member
!
      member = .false.
      if (n .ge. intval%first .and. n .le. intval%last) member = .true.
!
   end function element_is_member_interval
!
!
end module interval_class
