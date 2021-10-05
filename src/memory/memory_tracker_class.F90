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
module memory_tracker_class
!
!!
!!    Memory tracker
!!    Written by Eirik F. Kjønstad, June 2021
!!
!!    Keeps track of memory used and causes error stops when the memory exceeds the
!!    specified maximum. The tracker object must be notified about allocations and
!!    deallocations (via "update").
!!
!
   use kinds
   use global_out, only: output
!
   type :: memory_tracker
!
      integer(i64), private :: current
      integer(i64), private :: max_allowed
!
   contains
!
      procedure, public :: update
!
      final :: destructor
!
   end type memory_tracker
!
!
   interface memory_tracker
!
      procedure :: new_memory_tracker
!
   end interface memory_tracker
!
!
contains
!
!
   function new_memory_tracker(max_allowed) result(this)
!!
!!    New memory tracker
!!    Written by Eirik F. Kjønstad, June 2021
!!
!!    max_allowed: the maximum memory usage allowed (in bytes)
!!
      implicit none
!
      type(memory_tracker) :: this
!
      integer(i64), intent(in) :: max_allowed
!
      this%current = 0
      this%max_allowed = max_allowed
!
      call output%printf('debug', 'Memory tracker initialized - allowed memory: (i0) B', &
                         ints=[int(this%max_allowed)])
!
   end function new_memory_tracker
!
!
   subroutine update(this, bytes_allocated)
!!
!!    Update
!!    Written by Eirik F. Kjønstad, June 2021
!!
      implicit none
!
      class(memory_tracker), intent(inout) :: this
!
      integer(i64), intent(in) :: bytes_allocated
!
      this%current = this%current + bytes_allocated
!
      call output%printf('debug', 'Current memory used ((i0)/(i0))', &
                           ints=[int(this%current), int(this%max_allowed)])
!
      if (this%current .gt. this%max_allowed) then
!
         call output%printf('v','Exceeded expected memory in memory tracker. &
                           &(i0) B used out of (i0) B allowed.', &
                            ints=[int(this%current), int(this%max_allowed)])
!
      endif
!
   end subroutine update
!
!
   subroutine destructor(this)
!!
!!    Destructor
!!    Written by Eirik F. Kjønstad, June 2021
!!
      implicit none
!
      type(memory_tracker) :: this
!
      call output%printf('debug', 'Memory tracker finalized - allowed memory: (i0) B', &
                         ints=[int(this%max_allowed)])
!
   end subroutine destructor
!
!
end module memory_tracker_class
