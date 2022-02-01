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
module observer_class
!
!!
!! Observer class
!! Written by Sarai D. Folkestad, Oct 2021
!!
!! Part of the general implementation of the
!! Observer design pattern
!!
!! A concrete observer must implement an 'update'
!! which is used by the 'observable'.
!!
!
   use parameters
!
   implicit none
!
   type, abstract :: observer
!
   contains
!
      procedure(update_observer), deferred :: update
!
   end type observer
!
   abstract interface
!
      subroutine update_observer(this)
!!
         import observer
!
         implicit none
!
         class(observer), intent(inout) :: this
!
      end subroutine
!
   end interface
!
end module observer_class
