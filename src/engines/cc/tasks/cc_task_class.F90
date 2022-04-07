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
module cc_task_class
!
!!
!! CC task class
!! Written by Eirik F. Kjønstad, 2021
!!
!
   use ccs_class, only: ccs
   use banner_printer_class, only: banner_printer
   use timings_class, only: timings
!
   implicit none
!
   type, abstract :: cc_task
!
      character(len=200) :: name_
      type(banner_printer) :: banner
      type(timings) :: timer
!
   contains
!
      procedure(execute), public, deferred :: execute
!
      procedure, public :: print_header
!
      procedure, public :: start_timer
      procedure, public :: end_timer
!
   end type cc_task
!
!
   abstract interface
!
      subroutine execute(this, wf)
!
         import :: ccs, cc_task
!
         implicit none
!
         class(cc_task), intent(inout) :: this
!
         class(ccs), intent(inout), target :: wf
!
      end subroutine execute
!
   end interface
!
!
contains
!
!
   subroutine start_timer(this)
!!
!!    Start timer
!!    Written by Eirik F. Kjønstad, 2022
!!
      implicit none
!
      class(cc_task) :: this
!
      this%timer = timings(this%name_, pl='m')
      call this%timer%turn_on()
!
   end subroutine start_timer
!
!
   subroutine end_timer(this)
!!
!!    End timer
!!    Written by Eirik F. Kjønstad, 2022
!!
      implicit none
!
      class(cc_task) :: this
!
      call this%timer%turn_off()
!
   end subroutine end_timer
!
!
   subroutine print_header(this)
!!
!!    Print header
!!    Written by Eirik F. Kjønstad, 2022
!!
      implicit none
!
      class(cc_task) :: this
!
      call this%banner%print_(this%name_)
!
   end subroutine print_header
!
!
end module cc_task_class
