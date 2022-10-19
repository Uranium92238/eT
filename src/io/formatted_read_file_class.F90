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
module formatted_read_file_class
!
!!
!!    formatted read file class module
!!    Written by Rolf H. Myhre and Alexander C. Paul, 2019-2022
!!
!
   use kinds
   use abstract_file_class, only : abstract_file
   use global_out, only : output
!
   type, extends(abstract_file) :: formatted_read_file
!
   contains
!
      procedure, public :: skip => skip_formatted_read_file
      procedure, public :: io_error => io_error_formatted_read_file
!
      procedure, private :: read_r_formatted_read_file
      procedure, private :: read_r_1_formatted_read_file
      procedure, private :: read_r_2_formatted_read_file
      procedure, private :: read_character_formatted_read_file
!
      generic, public :: read_ => read_r_formatted_read_file,   &
                                  read_r_1_formatted_read_file, &
                                  read_r_2_formatted_read_file, &
                                  read_character_formatted_read_file
!
   end type formatted_read_file
!
   interface formatted_read_file
!
      procedure new_formatted_read_file
!
   end interface formatted_read_file
!
contains
!
!
   function new_formatted_read_file(name_) result(this)
!!
!!    new formatted read file
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      type(formatted_read_file) :: this
!
      character(len=*), intent(in) :: name_
!
      this%name_ = name_
!
      this%access_ = 'sequential'
      this%format_ = 'formatted'
      this%status_ = 'unknown'
      this%action_ = 'read'
!
      this%is_open = .false.
      this%unit_ = -1
!
   end function new_formatted_read_file
!
!
   subroutine skip_formatted_read_file(this, n_lines)
!!
!!    Skip
!!    Written by Rolf H. Myhre, May 2019
!!    Skip number of lines. Default: 1
!!
      implicit none
!
      class(formatted_read_file), intent(inout) :: this
      integer, optional, intent(in) :: n_lines
!
      integer :: i, n_lines_to_skip
!
      integer :: io_status
      character(len=100) :: io_msg
!
      n_lines_to_skip = 1
      if(present(n_lines)) n_lines_to_skip = n_lines
!
      do i = 1, n_lines_to_skip
!
         read(this%unit_, iostat=io_status, iomsg=io_msg)
         call this%check_io_status(io_status, io_msg, task='skip lines')
!
      enddo
!
   end subroutine skip_formatted_read_file
!
!
   subroutine read_r_formatted_read_file(this, scalar, io_stat)
!!
!!    read real(dp) scalar
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(formatted_read_file), intent(in)  :: this
!
      real(dp), intent(out)          :: scalar
      integer, intent(out), optional :: io_stat
!
      integer :: io_status
      character(len=100) :: io_msg
!
      read(this%unit_, *, iostat=io_status, iomsg=io_msg) scalar
      call this%check_io_status(io_status, io_msg, task='read', status_=io_stat)
!
   end subroutine read_r_formatted_read_file
!
!
   subroutine read_r_1_formatted_read_file(this, array, n, io_stat)
!!
!!    read real(dp) 1 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(formatted_read_file), intent(in)  :: this
!
      integer, intent(in)                 :: n
      real(dp), dimension(n), intent(out) :: array
      integer, intent(out), optional      :: io_stat
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      read(this%unit_, *, iostat=io_error, iomsg=io_msg) array
      call this%check_io_status(io_error, io_msg, task='read', status_=io_stat)
!
   end subroutine read_r_1_formatted_read_file
!
!
   subroutine read_r_2_formatted_read_file(this, array, n, io_stat)
!!
!!    read real(dp) 2 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(formatted_read_file), intent(in)     :: this
      integer, intent(in)                    :: n
      real(dp), dimension(:,:), intent(out)  :: array
      integer, intent(out), optional         :: io_stat
!
      call this%read_r_1_formatted_read_file(array, n, io_stat)
!
   end subroutine read_r_2_formatted_read_file
!
!
   subroutine read_character_formatted_read_file(this, string, format_string, io_stat)
!!
!!    read character
!!    Written by Rolf H. Myhre and Andreas Skeidsvol, 2019
!!
!!    string:  character string to read in to
!!    format_string: optional format string
!!    io_stat: optional integer to return io_error if it is less than or equal to 0
!!
      implicit none
!
      class(formatted_read_file), intent(in) :: this
!
      character(len=*), intent(out)          :: string
      character(len=*), intent(in), optional :: format_string
      integer, intent(out), optional         :: io_stat
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      if (present(format_string) .and. trim(format_string) .ne. '*') then
         read(this%unit_, trim(format_string), iostat=io_error, iomsg=io_msg) string
      else
         read(this%unit_, *, iostat=io_error, iomsg=io_msg) string
      endif
!
      call this%check_io_status(io_error, io_msg, task='read', status_=io_stat)
!
   end subroutine read_character_formatted_read_file
!
!
   subroutine io_error_formatted_read_file(this, io_status, io_message, task)
!!
!!    print io error
!!    Written by Alexander C. Paul, Sep 2022
!!
!!    task: string specifying which task is checked: open, close, read from, write to
!!
      implicit none
!
      class(formatted_read_file), intent(in) :: this
!
      integer, intent(in) :: io_status
!
      character(len=*), intent(in) :: io_message
      character(len=*), intent(in) :: task
!
      call output%error_msg('Failed to '// task // 'formatted read file (a0), status &
                            &is (i0) and error message is: ' // trim(io_message), &
                            chars = [this%get_name()], ints = [io_status])
!
   end subroutine io_error_formatted_read_file
!
!
end module formatted_read_file_class
