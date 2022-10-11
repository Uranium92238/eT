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
module abstract_file_class
!
!!
!! Abstract file class module
!! Written by Sarai D. Folkestad, Eirik F. Kjønstad,
!! Rolf H. Myhre and Alexander C. Paul, 2018-2022
!!
!
   use kinds
!
   type, abstract :: abstract_file
!
!     These variables cannot be private because of problems in the constructors
      character(len=255) :: name_   = 'no_name'
      character(len=12)  :: access_ = 'unknown'
      character(len=12)  :: format_ = 'unknown'
!
!     Unit identifier will be assigned by new_unit
      integer :: unit_ = -1
!
      logical :: is_open = .false.
!
      character(len=12) :: status_ = 'unknown'
      character(len=12) :: action_ = 'unknown'
!
   contains
!
      procedure, public :: open_
      procedure, public :: close_
      procedure, public :: delete_
!
      procedure, public :: rewind_
      procedure, public :: exists => exists_
      procedure, public :: determine_status
!
      procedure, public :: get_name
!
      procedure, public :: check_io_status
      procedure, public :: error_if_closed
      procedure, public :: error_if_open
!
      procedure(io_error), deferred :: io_error
!
   end type abstract_file
!
!
   abstract interface
!
      subroutine io_error(this, io_status, io_message, task)
!
         import abstract_file
         implicit none
!
         class(abstract_file), intent(in) :: this
!
         integer, intent(in) :: io_status
         character(len=*), intent(in) :: io_message
         character(len=*), intent(in) :: task
!
      end subroutine
!
   end interface
!
!
contains
!
!
   subroutine open_(this, new_position)
!!
!!    Open
!!    Written by Rolf H. Myhre Nov. 2019
!!
!!    Opens the file
!!    new_position: optional character string specifying position
!!                  rewind, append or asis. Default: asis
!!
      implicit none
!
      class(abstract_file), intent(inout) :: this
!
      character(len=*), optional, intent(in) :: new_position
!
      integer :: io_status
      character(len=10) :: position_
      character(len = 100) :: io_message
!
      call this%error_if_open('open')
!
      position_ = 'asis'
      if(present(new_position)) position_ = new_position
!
!     We let open handle the error if position is wrong
!
      open(newunit=this%unit_, file=this%get_name(), access=this%access_, &
           form=this%format_, action=this%action_, status=this%status_, &
           position=position_, iostat=io_status, iomsg=io_message)
!
      call this%check_io_status(io_status, io_message, task='open')
!
      this%is_open = .true.
      this%status_ = 'old'
!
   end subroutine open_
!
!
   subroutine close_(this, destiny)
!!
!!    Close
!!    Written by Rolf H. Myhre Nov. 2019
!!
!!    Closes the file
!!    Destiny: Optional character string that decides whether to keep
!!             or delete. Deafault: keep
!!
      implicit none
!
      class(abstract_file), intent(inout) :: this
!
      character(len=*), intent(in), optional :: destiny
!
      integer            :: io_status
      character(len=100) :: io_message
      character(len=10)  :: local_destiny
!
      call this%error_if_closed('close')
!
      local_destiny = 'keep'
      if(present(destiny)) local_destiny = destiny
!
!     We let close handle the error if local_destiny is wrong
      close(this%unit_, status=local_destiny, iostat=io_status, iomsg=io_message)
!
      call this%check_io_status(io_status, io_message, task='close')
!
      this%is_open = .false.
      this%unit_ = -1
!
      if (local_destiny == 'keep') then
         this%status_ = 'old'
      else
         this%status_ = 'new'
      endif
!
   end subroutine close_
!
!
   subroutine delete_(this)
!!
!!    Delete
!!    Written by Rolf H. Myhre, Aug 2019
!!
      implicit none
!
      class(abstract_file) :: this
!
      if(.not. this%is_open) call this%open_()
      call this%close_("delete")
!
   end subroutine delete_
!
!
   function get_name(this) result(name_)
!!
!!    Get name
!!    Written by and Alexander C. Paul, 2022
!!
      implicit none
!
      class(abstract_file), intent(in) :: this
      character(len=:), allocatable :: name_
!
      name_ = trim(this%name_)
!
   end function get_name
!
!
   subroutine rewind_(this)
!!
!!    Rewind
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(abstract_file), intent(inout) :: this
!
      integer :: io_status
      character(len=100) :: io_msg
!
      call this%error_if_closed('rewind')
      rewind(this%unit_, iostat=io_status, iomsg=io_msg)
      call this%check_io_status(io_status, io_msg, task='rewind')
!
   end subroutine rewind_
!
!
   function exists_(this) result(exists)
!!
!!    Exists
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad and Rolf H. Myhre 2018
!!
      implicit none
!
      class(abstract_file), intent(in) :: this
!
      logical :: exists
!
      inquire(file=trim(this%name_), exist=exists)
!
   end function exists_
!
!
   subroutine determine_status(this)
!!
!!    Determine status
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none
!
      class(abstract_file), intent(inout) :: this
!
      if (this%exists()) then
         this%status_ = 'old'
      else
         this%status_ = 'new'
      end if
!
   end subroutine determine_status
!
!
   subroutine check_io_status(this, io_status, io_message, task, status_)
!!
!!    Check io status
!!    Written by Alexander C. Paul, Oct 2020
!!
!!    Checks IO status and prints error message according to task and io_message
!!
!!    task: string specifying which task is checked: open, close, read from, write to
!!
!!    status_: Optionally returns the io_status to be handled outside
!!             e.g. for checking if the end of file was reached.
!!
!!    io_status < 0: end of file or end of record
!!    io_status > 0: error occured that prevents further execution
!!
      implicit none
!
      class(abstract_file), intent(in) :: this
!
      integer, intent(in) :: io_status
!
      character(len=*), intent(in) :: io_message
      character(len=*), intent(in) :: task
!
      integer, intent(out), optional :: status_
!
!
      if(present(status_) .and. io_status <= 0) then
!
         status_ = io_status
!
      else if(io_status > 0) then
!
         call this%io_error(io_status, io_message, task)
!
      endif
!
   end subroutine check_io_status
!
!
   subroutine error_if_closed(this, task)
!!
!!    Error if closed
!!    Written by Alexander C. Paul, Oct 2020
!!
!!    task: task for which an open file is expected (e.g. close, read, write)
!!
      implicit none
!
      class(abstract_file), intent(in) :: this
      character(len=*), intent(in) :: task
!
      if(.not. this%is_open) then
!
         print '(5a)', 'Trying to ' // task // ' file ' // this%get_name() // &
                      'which is closed.'
         error stop
!
      end if
!
   end subroutine error_if_closed
!
!
   subroutine error_if_open(this, task)
!!
!!    Error if open
!!    Written by Alexander C. Paul, Oct 2020
!!
!!    task: task for which a closed file is expected (e.g. open)
!!
      implicit none
!
      class(abstract_file), intent(in) :: this
!
      character(len=*), intent(in) :: task
!
      if(this%is_open) then
!
         print '(5a)', 'Trying to ' // task // ' file ' // this%get_name() // &
                      'which is open.'
         error stop
!
      end if
!
   end subroutine error_if_open
!
!
end module abstract_file_class
