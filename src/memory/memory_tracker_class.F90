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
      integer(i64), private :: max_allowed, max_used
!
      character(len=:), allocatable :: tag
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
   function new_memory_tracker(max_allowed, tag) result(this)
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
      character(len=*), intent(in) :: tag
!
      this%current = 0
      this%max_used = 0
      this%max_allowed = max_allowed
!
      this%tag = trim(tag)
!
      call output%printf('debug', 'Memory tracker initialized - ' // this%tag, ll=100)
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
      character(len=17) :: current_string, max_string
!
      this%current = this%current + bytes_allocated
!
      if (this%current > this%max_used) this%max_used = this%current
!
      if (this%current .gt. this%max_allowed) then
!
         write(current_string, '(i0)') this%current
         write(max_string, '(i0)') this%max_allowed
!
         call output%printf('v','Exceeded expected memory in memory tracker: (a0)', &
                             chars=[this%tag])
         call output%printf('v', '(a0) B used out of (a0) B allowed.', &
                            chars=[current_string, max_string])
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
      character(len=17) :: difference, allowed_string, used_string
!
      write(allowed_string, '(i0)') this%max_allowed
      write(used_string, '(i0)') this%current
      write(difference,'(i0)') this%max_allowed - this%max_used
!
      call output%printf('debug', 'Memory tracker finalized - ' // this%tag // ':')
      call output%printf('debug', 'allowed: (i0) B  used: (i0) B  difference: (i0) B', &
                         chars=[allowed_string, used_string, difference], ll=100)
!
   end subroutine destructor
!
!
end module memory_tracker_class
