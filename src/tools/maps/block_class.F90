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
module block_class
!
!!
!!    Block class
!!    Written by Rolf H. Myhre, Apr 2021
!!
!!    Class with two ranges interpreted as the dimensions of a matrix.
!
   use range_class
!
   implicit none
!
   type :: block_
!
      type(range_) :: p_range
      type(range_) :: q_range
!
   contains
!
      procedure :: overlaps      => overlaps_block
!
      procedure :: contains_     => contains_block
!
      procedure :: get_overlap   => get_overlap_block
!
      procedure :: get_area      => get_area_block
      procedure :: get_extent    => get_extent_block
      procedure :: get_position  => get_position_block
!
      procedure :: get_p_range   => get_p_range_block
      procedure :: get_q_range   => get_q_range_block
!
      procedure :: is_diagonal   => is_diagonal_block
!
      procedure :: is_equal_block
      generic, public :: operator(.eq.) => is_equal_block
!
   end type block_
!
!
   interface block_
!
      procedure :: new_integers_block
      procedure :: new_ranges_block
      procedure :: copy_block
!
   end interface block_
!
!
contains
!
!
   elemental function new_integers_block(first_p, last_p, first_q, last_q) result(this)
!!
!!    New integers block
!!    Written by Rolf H. Myhre, Apr 2021
!!
      implicit none
!
      integer, intent(in) :: first_p, last_p, first_q, last_q
!
      type(block_) :: this
!
      this%p_range = range_(first_p, last_p)
      this%q_range = range_(first_q, last_q)
!
   end function new_integers_block
!
!
   elemental function new_ranges_block(p_range, q_range) result(this)
!!
!!    New ranges block
!!    Written by Rolf H. Myhre, Apr 2021
!!
      implicit none
!
      class(range_), intent(in) :: p_range, q_range
!
      type(block_) :: this
!
      this%p_range = p_range
      this%q_range = q_range
!
   end function new_ranges_block
!
!
   elemental function copy_block(that) result(this)
!!
!!    Copy block
!!    Written by Rolf H. Myhre, Apr 2021
!!
      implicit none
!
      class(block_), intent(in) :: that
!
      type(block_) :: this
!
      this%p_range = that%get_p_range()
      this%q_range = that%get_q_range()
!
   end function copy_block
!
!
   elemental function is_equal_block(this, that) result(is_equal)
!!
!!    Is equal
!!    Written by Rolf H. Myhre, Apr 2021
!!
      implicit none
!
      class(block_), intent(in) :: this
      class(block_), intent(in) :: that
!
      logical :: is_equal
!
      is_equal = ((this%p_range .eq. that%p_range) .and. (this%q_range .eq. that%q_range))
!
   end function is_equal_block
!
!
   elemental function overlaps_block(this, that) result(overlaps)
!!
!!    Overlaps
!!    Written by Rolf H. Myhre, Apr 2021
!!
!!    Does that overlap with this
!!
      implicit none
!
      class(block_), intent(in) :: this
      class(block_), intent(in) :: that
!
      logical :: overlaps
!
      overlaps = this%p_range%overlaps(that%p_range) .and. &
                 this%q_range%overlaps(that%q_range)
!
   end function overlaps_block
!
!
   elemental function get_overlap_block(this, that) result(overlap)
!!
!!    Get overlap
!!    Written by Rolf H. Myhre, Apr 2021
!!
      implicit none
!
      class(block_), intent(in) :: this
      class(block_), intent(in) :: that
!
      type(block_) :: overlap
!
      if (this%overlaps(that)) then
         overlap = block_(this%p_range%get_overlap(that%get_p_range()), &
                          this%q_range%get_overlap(that%get_q_range()))
      endif
!
   end function get_overlap_block
!
!
   elemental function contains_block(this, that) result(contains_)
!!
!!    Contains
!!    Written by Rolf H. Myhre, Apr 2021
!!
!!    Determine if that is contained in this
!!
      implicit none
!
      class(block_), intent(in) :: this
      class(block_), intent(in) :: that
!
      logical :: contains_
!
      contains_ = this%p_range%contains_(that%p_range) .and. &
                  this%q_range%contains_(that%q_range)
!
   end function contains_block
!
!
   elemental function get_area_block(this) result(area)
!!
!!    Get area
!!    Written by Rolf H. Myhre, Apr 2021
!!
      implicit none
!
      class(block_), intent(in) :: this
!
      integer :: area
!
      area = this%p_range%length*this%q_range%length
!
   end function get_area_block
!
!
   elemental function get_extent_block(this) result(extent)
!!
!!    Get extent
!!    Written by Rolf H. Myhre, Apr 2021
!!
      implicit none
!
      class(block_), intent(in) :: this
!
      integer :: extent
!
      extent = this%p_range%get_extent()*this%q_range%get_extent()
!
   end function get_extent_block
!
!
   elemental function get_position_block(this, p, q) result(position_)
!!
!!    Get position
!!    Written by Rolf H. Myhre, Apr 2021
!!
      implicit none
!
      class(block_), intent(in) :: this
!
      integer, intent(in) :: p, q
      integer :: position_
!
      position_ = this%p_range%get_position(p) + &
                 (this%q_range%get_position(q) - 1) * this%p_range%get_extent()
!
   end function get_position_block
!
!
   elemental function is_diagonal_block(this) result(diagonal)
!!
!!    Is diagonal
!!    Written by Rolf H. Myhre, Apr 2021
!!
      implicit none
!
      class(block_), intent(in) :: this
!
      logical :: diagonal
!
      diagonal = (this%p_range .eq. this%q_range)
!
   end function is_diagonal_block
!
!
   elemental function get_p_range_block(this) result(p_range)
!!
!!    Get p range
!!    Written by Rolf H. Myhre, Apr 2021
!!
      implicit none
!
      class(block_), intent(in) :: this
!
      type(range_) :: p_range
!
      p_range = this%p_range
!
   end function get_p_range_block
!
!
   elemental function get_q_range_block(this) result(q_range)
!!
!!    Get q range
!!    Written by Rolf H. Myhre, Apr 2021
!!
      implicit none
!
      class(block_), intent(in) :: this
!
      type(range_) :: q_range
!
      q_range = this%q_range
!
   end function get_q_range_block
!
!
end module block_class
