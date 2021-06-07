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
module range_class
!!
!!    Range class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!    Changed to range by Rolf H. Myhre, Jun 2020
!!
   implicit none
!
   type range_
!
      integer, public :: first  = 0
      integer, public :: length = 0
      integer, public :: step   = 0
!
   contains
!
      procedure, public    :: shift          => shift_range
      procedure, public    :: shift_first_to => shift_first_to_range
      procedure, public    :: shift_last_to  => shift_last_to_range
!
      procedure, public    :: get_first      => get_first_range
      procedure, public    :: get_last       => get_last_range
      procedure, public    :: get_length     => get_length_range
      procedure, public    :: get_step       => get_step_range
!
      procedure, public    :: get_min        => get_min_range
      procedure, public    :: get_max        => get_max_range
      procedure, public    :: get_abs_step   => get_abs_step_range
!
      procedure, public    :: overlaps       => overlaps_range
      procedure, public    :: get_overlap    => get_overlap_range
!
      procedure, private   :: contains_range
      procedure, private   :: contains_integer
      generic, public      :: contains_     => contains_range, &
                                            contains_integer
!
      procedure, private   :: is_equal_range
      generic, public      :: operator(.eq.) => is_equal_range
!
   end type range_
!
   interface range_
!
      procedure :: new_range
      procedure :: copy_range
!
   end interface range_
!
contains
!
!
   elemental function new_range(first, length, step) result(this)
!!
!!    New range
!!    Written by Rolf H. Myhre, Apr 2021
!!
      implicit none
!
      integer, intent(in) :: first, length
      integer, optional, intent(in) :: step
!
      type(range_) :: this
!
      integer :: step_
!
      step_ = 1
      if (present(step)) step_ = step
!
      if (length .gt. 0 .and. step_ .ne. 0) then
         this%first  = first
         this%length = length
         this%step   = step_
      else
         error stop "length < 1 or step == 0 in new_range"
      endif
!
   end function new_range
!
!
   elemental function copy_range(that) result(this)
!!
!!    Copy
!!    Written by Rolf H. Myhre, Apr 2021
!!
      implicit none
!
      class(range_), intent(in)   :: that
!
      type(range_) :: this
!
      this%first  = that%first
      this%length = that%length
      this%step   = that%step
!
   end function copy_range
!
!
   elemental function range_from_first_last(first, last, step) result(this)
!!
!!    Range from first and last
!!    Written by Rolf H. Myhre, Apr 2021
!!
      implicit none
!
      integer, intent(in) :: first, last
      integer, optional, intent(in) :: step
!
      type(range_) :: this
!
      integer :: step_
!
      step_ = 1
      if (present(step)) step_ = step
!
      if (((last - first) .ge. 0 .and. step_ .gt. 0) .or. &
          ((last - first) .lt. 0 .and. step_ .lt. 0)) then
         this%first  = first
         this%length = (last - first)/step_ + 1
         this%step   = step_
      else
         error stop "Inconsistent first, last, and step in range_from_first_last"
      endif
!
   end function range_from_first_last
!
!
   elemental function is_equal_range(this, that) result(is_equal)
!!
!!    Is equal
!!    Written by Rolf H. Myhre, Apr 2021
!!
      implicit none
!
      class(range_), intent(in) :: this
      class(range_), intent(in) :: that
!
      logical :: is_equal
!
      is_equal = ((this%first  .eq. that%first)  .and. &
                  (this%length .eq. that%length) .and. &
                  (this%step   .eq. that%step))
!
   end function is_equal_range
!
!
   elemental subroutine shift_range(this, n)
!!
!!    Shift
!!    Written by Rolf H. Myhre, Apr 2021
!!
      implicit none
!
      class(range_), intent(inout) :: this
      integer, intent(in)          :: n
!
      this%first = this%first + n
!
   end subroutine shift_range
!
!
   elemental subroutine shift_first_to_range(this, n)
!!
!!    Set first to
!!    Written by Rolf H. Myhre, Apr 2021
!!
      implicit none
!
      class(range_), intent(inout) :: this
      integer, intent(in)          :: n
!
      call this%shift(n - this%first)
!
   end subroutine shift_first_to_range
!
!
   elemental subroutine shift_last_to_range(this, n)
!!
!!    Set last to
!!    Written by Rolf H. Myhre, Apr 2021
!!
      implicit none
!
      class(range_), intent(inout) :: this
      integer, intent(in)          :: n
!
      call this%shift(n - this%get_last())
!
   end subroutine shift_last_to_range
!
!
   elemental function get_first_range(this) result(first)
!!
!!    Get first
!!    Written by Rolf H. Myhre, Apr 2021
!!
      implicit none
!
      class(range_), intent(in) :: this
!
      integer :: first
!
      first = this%first
!
   end function get_first_range
!
!
   elemental function get_length_range(this) result(length)
!!
!!    Get length
!!    Written by Rolf H. Myhre, Apr 2021
!!
      implicit none
!
      class(range_), intent(in) :: this
!
      integer :: length
!
      length = this%length
!
   end function get_length_range
!
!
   elemental function get_step_range(this) result(step)
!!
!!    Get step
!!    Written by Rolf H. Myhre, Apr 2021
!!
      implicit none
!
      class(range_), intent(in) :: this
!
      integer :: step
!
      step = this%step
!
   end function get_step_range
!
!
   elemental function get_last_range(this) result(last)
!!
!!    Get last
!!    Written by Rolf H. Myhre, Apr 2021
!!
      implicit none
!
      class(range_), intent(in) :: this
!
      integer :: last
!
      last = this%first + (this%length-1)*this%step
!
   end function get_last_range
!
!
   elemental function get_min_range(this) result(min_val)
!!
!!    Get min
!!    Written by Rolf H. Myhre, Apr 2021
!!
      implicit none
!
      class(range_), intent(in) :: this
!
      integer :: min_val
!
      if (this%step .gt. 0) then
         min_val = this%first
      else
         min_val = this%get_last()
      endif
!
   end function get_min_range
!
!
   elemental function get_max_range(this) result(max_val)
!!
!!    Get max
!!    Written by Rolf H. Myhre, Apr 2021
!!
      implicit none
!
      class(range_), intent(in) :: this
!
      integer :: max_val
!
      if (this%step .gt. 0) then
         max_val = this%get_last()
      else
         max_val = this%first
      endif
!
   end function get_max_range
!
!
   elemental function get_abs_step_range(this) result(abs_step)
!!
!!    Get absolute step
!!    Written by Rolf H. Myhre, Apr 2021
!!
      implicit none
!
      class(range_), intent(in) :: this
!
      integer :: abs_step
!
      abs_step = abs(this%step)
!
   end function get_abs_step_range
!
!
   elemental function overlaps_range(this, that) result(overlaps)
!!
!!    Overlaps
!!    Written by Rolf H. Myhre, Apr 2021
!!
!!    Does 'that' have elements in common with 'this'?
!!
!!    Determines whether the ranges
!!
!!    range1 = f1 + n1 * s1
!!    range2 = f2 + n2 * s2
!!
!!    overlap by checking if the integer equation
!!
!!    n1*s1 - n2*s2 = f2 - f1
!!
!!    has any solutions.
!!
      use math_utilities, only: gcd, least_positive_n1s1_n2s2
!
      implicit none
!
      class(range_), intent(in) :: this
      class(range_), intent(in) :: that
!
      logical :: overlaps
!
      integer :: g, l, offset
      integer :: n1, n2
!
!     First check that the first and last values overlap.
      if (this%get_min() .le. that%get_max() .and. &
          that%get_min() .le. this%get_max()) then
!
!        If offset is not dividable by greatest common divisor,
!        the ranges will not overlap.
         g = gcd(this%step, that%step)
         offset = that%get_min() - this%get_min()
         if (mod(offset, g) .eq. 0) then
!
!           If the overlap is greater than the least common multiple,
!           there will be an overlap.
            l = (this%get_abs_step()/g)*that%get_abs_step()
            if ((min(this%get_max(),  that%get_max()) -          &
                 max(this%get_min(), that%get_min()) + 1) .ge. l) then
               overlaps = .true.
               return
            endif
!
!           Find first positive n1 corresponding to a positive n2
            n1 = least_positive_n1s1_n2s2(this%get_abs_step(), that%get_abs_step(), offset)
            n2 = (n1*this%get_abs_step() - offset)/that%get_abs_step()
!
!           Check if both n1 and n2 are less than length
            if ((n1 .lt. this%length) .and. (n2 .lt. that%length)) then
               overlaps = .true.
               return
            endif
         endif
      endif
!
      overlaps = .false.
!
   end function overlaps_range
!
!
   elemental function contains_range(this, that) result(contains_)
!!
!!    Contains
!!    Written by Rolf H. Myhre, Apr 2021
!!
!!    Is that contained in this?
!!
!!    Check if the ranges overlap, if they do,
!!    check that that%min is greater than or equal to this%min
!!    and that%max .le. this%max
!!    and that%abs_step is a multiple of this%abs_step or that%abs_step is 1
!!
      implicit none
!
      class(range_), intent(in) :: this
      class(range_), intent(in) :: that
!
      logical :: contains_
!
      if (this%overlaps(that)) then
!
         if (that%get_min() .ge. this%get_min() .and. &
             that%get_max() .le. this%get_max()) then
!
            if ((that%get_abs_step() .ge. this%get_abs_step() .and. &
                 mod(that%get_abs_step(), this%get_abs_step()) .eq. 0) .or. &
                (that%length .eq. 1)) then
               contains_ = .true.
               return
            endif
!
         endif
!
      endif
!
      contains_ = .false.
!
   end function contains_range
!
!
   elemental function contains_integer(this, an_integer) result(contains_)
!!
!!    Contains integer
!!    Written by Rolf H. Myhre, Apr 2021
!!
!!    Is 'an_integer' contained in 'this'?
!!
      implicit none
!
      class(range_), intent(in) :: this
      integer, intent(in)       :: an_integer
!
      logical :: contains_
!
      if (an_integer .ge. this%get_min() .and. &
          an_integer .le. this%get_max()) then
         if (mod(an_integer - this%get_min(), this%get_abs_step()) .eq. 0) then
            contains_ = .true.
            return
         endif
      endif
!
      contains_ = .false.
!
   end function contains_integer
!
!
   elemental function get_overlap_range(this, that) result(overlap)
!!
!!    Get overlap
!!    Written by Rolf H. Myhre, Apr 2021
!!
!!    Returns a range with the overlap of the two ranges with positive steplength
!!
!!    first is the first common element
!!    step is the least common multiple
!!    length is (min(this%last, that%last) - first)/step + 1)
!!
!!    If the ranges don't overlap, length will be zero
!!
      use math_utilities, only: lcm, least_positive_n1s1_n2s2
!
      implicit none
!
      class(range_), intent(in) :: this
      class(range_), intent(in) :: that
!
      type(range_) :: overlap
!
      integer :: n1
      integer :: first, length, step
!
      if (this%overlaps(that)) then
!
         n1 = least_positive_n1s1_n2s2(this%get_abs_step(), that%get_abs_step(), &
                                       that%get_min() - this%get_min())
!
         first  = this%get_min() + n1*this%get_abs_step()
         step   = lcm(this%get_abs_step(), that%get_abs_step())
         length = (min(this%get_max(), that%get_max()) - first)/step + 1
!
         overlap = range_(first, length, step)
!
      endif
!
   end function get_overlap_range
!
!
end module range_class
