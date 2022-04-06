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
module observable_class
!
!!
!! Observable class
!! Written by Sarai D. Folkestad, Oct 2021
!!
!! Part of the general implementation of the
!! Observer design pattern
!!
!! An 'observable' keeps a list of many 'observers' which it can 'notify'
!! when the observable changes. The observers 'update' when
!! they are notified.
!!
!
   use parameters
!
   use observer_list_class, only: observer_list
   use observer_class, only: observer
!
   implicit none
!
   type, abstract :: observable
!
      type(observer_list) :: observers
!
   contains
!
      procedure :: add_observer
      procedure :: remove_observer
      procedure :: notify_observers
!
   end type observable
!
contains
!
   subroutine add_observer(this, tag, o)
!!
!!    Add observer
!!    Written by Sarai D. Folkestad, Oct 2021
!!
      implicit none
!
      class(observable),   intent(inout) :: this
!
      class(observer),     intent(inout)   :: o
      character(len=*),    intent(in)  :: tag
!
      call this%observers%add(tag, o)
!
   end subroutine add_observer
!
!
   subroutine remove_observer(this, tag)
!!
!!    Remove observer
!!    Written by Sarai D. Folkestad, Oct 2021
!!
      implicit none
!
      class(observable), intent(inout) :: this
      character(len=*),  intent(in)  :: tag
!
      call this%observers%remove(tag)
!
   end subroutine remove_observer
!
!
   subroutine notify_observers(this)
!!
!!    Notify observers
!!    Written by Sarai D. Folkestad, Oct 2021
!!
      implicit none
!
      class(observable),   intent(inout) :: this
!
      call this%observers%update()
!
   end subroutine notify_observers
!
!
end module observable_class
