!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2019 the authors of eT
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
module abstract_stream_class
!
!!
!!    Abstract stream class module
!!    Written by Rolf H. Myhre, Feb. 2020
!!
!
   use kinds
   use global_out, only: output
!
   type, abstract :: abstract_stream
!
!     Filename
!
      character(len=255), private :: name_ = 'no_name'
      character(len=10), private  :: status_ = 'no_status'
!
!     Unit identifier
!
      integer, private :: unit_ = 1
!
!     Logical for whether the file is currently opened or not
!
      logical, private :: is_open = .false.
!
   contains
!
      procedure :: set_name         => set_name_abstract_stream
      procedure :: get_name         => get_name_abstract_stream
! 
      procedure :: set_status       => set_status_abstract_stream
! 
      procedure :: get_open         => get_open_abstract_stream
! 
      procedure :: open_            => open_abstract_stream
      procedure :: close_           => close_abstract_stream
      procedure :: delete_          => delete_abstract_stream
! 
      procedure :: exists           => exists_abstract_stream
      procedure :: get_file_size    => get_file_size_abstract_stream
      procedure :: copy             => copy_abstract_stream
!
!     Error handling
      procedure :: check_io_status  => check_io_status_abstract_stream
      procedure :: should_be_open   => should_be_open_abstract_stream
      procedure :: should_be_closed => should_be_closed_abstract_stream
!
      procedure :: read_0_real_dp_abstract_stream
      procedure :: read_0_real_sp_abstract_stream
      procedure :: read_0_complex_dp_abstract_stream
      procedure :: read_0_int_32_abstract_stream
      procedure :: read_0_int_64_abstract_stream
      procedure :: read_0_log_abstract_stream
!
      procedure :: read_1_real_dp_abstract_stream
      procedure :: read_1_real_sp_abstract_stream
      procedure :: read_1_complex_dp_abstract_stream
      procedure :: read_1_int_32_abstract_stream
      procedure :: read_1_int_64_abstract_stream
      procedure :: read_1_log_abstract_stream
!
      procedure :: write_0_real_dp_abstract_stream
      procedure :: write_0_real_sp_abstract_stream
      procedure :: write_0_complex_dp_abstract_stream
      procedure :: write_0_int_32_abstract_stream
      procedure :: write_0_int_64_abstract_stream
      procedure :: write_0_log_abstract_stream
!
      procedure :: write_1_real_dp_abstract_stream
      procedure :: write_1_real_sp_abstract_stream
      procedure :: write_1_complex_dp_abstract_stream
      procedure :: write_1_int_32_abstract_stream
      procedure :: write_1_int_64_abstract_stream
      procedure :: write_1_log_abstract_stream
!
   end type abstract_stream
!
!
contains
!
!
   subroutine set_name_abstract_stream(the_file, new_name)
!!
!!    Set name
!!    Written by Rolf H. Myhre Nov. 2019
!!
!!    Sets the private file name variable
!!
      implicit none
!
      class(abstract_stream), intent(inout) :: the_file
!
      character(len=*), intent(in) :: new_name
!
      if (len_trim(new_name) .le. 0) then
!
         call output%error_msg("Blank string in set_name")
!
      end if
!
      the_file%name_ = trim(new_name)
!
   end subroutine set_name_abstract_stream
!
!
   function get_name_abstract_stream(the_file) result(name_)
!!
!!    Get name
!!    Written by Rolf H. Myhre Feb. 2020
!!
      implicit none
!
      class(abstract_stream), intent(in) :: the_file
!
      character(len=255), allocatable :: name_
!
      name_ = trim(the_file%name_)
!
   end function get_name_abstract_stream
!
!
   subroutine set_status_abstract_stream(the_file, new_status)
!!
!!    Set status
!!    Written by Rolf H. Myhre Nov. 2019
!!
!!    Sets the private file status variable
!!
      implicit none
!
      class(abstract_stream), intent(inout) :: the_file
!
      character(len=*), intent(in) :: new_status
!
      if (trim(new_status) .eq. 'new' .or. &
          trim(new_status) .eq. 'old' .or. &
          trim(new_status) .eq. 'unknown') then
!
         the_file%status_ = trim(new_status)
!
      else
!
         call output%error_msg('(a0) is not an acceptable status', chars=[new_status])
!
      end if
!
   end subroutine set_status_abstract_stream
!
!
   function get_open_abstract_stream(the_file) result(is_open)
!!
!!    Get open status
!!    Written by Rolf H. Myhre Feb. 2020
!!
      implicit none
!
      class(abstract_stream) :: the_file
!
      logical :: is_open
!
      is_open = the_file%is_open
!
   end function get_open_abstract_stream
!
!
   subroutine open_abstract_stream(the_file, new_action, new_position)
!!
!!    Open
!!    Written by Rolf H. Myhre Nov. 2019
!!
!!    Opens the file
!!
!!    new_action:   optional character string specifying actions
!!                  read, write or readwrite. Default: readwrite
!!
!!    new_position: optional character string specifying position
!!                  rewind, append or asis. Default: asis
!!
      implicit none
!
      class(abstract_stream), intent(inout) :: the_file
!
      character(len=*), optional, intent(in) :: new_action
      character(len=*), optional, intent(in) :: new_position
!
      integer              :: io_status
      character(len = 100) :: io_message
!
      character(len=10) :: action_
      character(len=10) :: position_
!
      call the_file%should_be_closed('open')
!
!     Check if new_action or new_position is present
!     We let open handle the error if they are messed up
!
      action_ = 'readwrite'
      if(present(new_action)) then
         action_ = new_action
      end if
!
      position_ = 'asis'
      if(present(new_position)) then
         position_ = new_position
      end if
!
      open(newunit = the_file%unit_, file=trim(the_file%name_), access = 'stream', &
           form = 'unformatted', action = action_, status = the_file%status_, &
           position = position_, iostat = io_status, iomsg = io_message)
!
      call the_file%check_io_status(io_status, trim(io_message), 'open')
!
      the_file%is_open = .true.
      the_file%status_ = 'old'
!
   end subroutine open_abstract_stream
!
!
   subroutine close_abstract_stream(the_file, new_destiny)
!!
!!    Close
!!    Written by Rolf H. Myhre Nov. 2019
!!
!!    Closes the file
!!
!!    Destiny: Optional character string that decides whether to keep
!!             or delete. Deafault: keep
!!
      implicit none
!
      class(abstract_stream), intent(inout) :: the_file
!
      character(len=*), intent(in), optional :: new_destiny
!
      integer              :: io_status
      character(len = 100) :: io_message
      character(len = 10)  :: destiny
!
      call the_file%should_be_open('close')
!
!     Check if destiny is present
!     We let close handle the error if destiny is messed up
!
      if(.not. present(new_destiny)) then
         destiny = 'keep'
      else
         destiny = new_destiny
      end if
!
!     Things seems fine, close the file
!
      close(the_file%unit_, status = destiny, iostat = io_status, iomsg = io_message)
!
!     Was it fine?
!
      call the_file%check_io_status(io_status, trim(io_message), 'close')
!
      the_file%is_open = .false.
      the_file%unit_ = 1
!
      if(destiny .eq. 'keep') then
         the_file%status_ = 'old'
      else
         the_file%status_ = 'new'
      endif
!
   end subroutine close_abstract_stream
!
!
   function exists_abstract_stream(the_file) result(it_exists)
!!
!!    File exists
!!    Written by Rolf H. Myhre Nov. 2019
!!
!!    Checks if file exists
!!
      implicit none
!
      class(abstract_stream), intent(inout) :: the_file
!
      logical :: it_exists
!
      inquire(file=trim(the_file%name_), exist=it_exists)
!
      if(it_exists) then
!
         the_file%status_ = 'old'
!
      else
!
         the_file%status_ = 'new'
!
      end if
!
!
   end function exists_abstract_stream
!
!
   function get_file_size_abstract_stream(the_file) result(file_size)
!!
!!    Get file Size
!!    Written by Rolf H. Myhre, Mar. 2020
!!
!!    Returns the number of the last byte written
!!
!!    Returns -1 if the number can't be determined,
!!    typically because the file does not exist.
!!
      implicit none
!
      class(abstract_stream), intent(inout) :: the_file
!
      integer :: file_size
!
      inquire(file=trim(the_file%name_), size=file_size)
!
   end function get_file_size_abstract_stream
!
!
   subroutine read_0_real_dp_abstract_stream(the_file, scalar, position_, status_)
!!
!!    read rank 0 real dp
!!    Written by Rolf H. Myhre, Feb. 2020
!!
!!    Modified by Alexander C. Paul, Mar. 2020
!!    Check if position_ is present
!!
!!    scalar: real double precision scalar
!!
!!    position_:  optional integer, position to read from in file
!!                positions counted in bytes, starting at 1
!!                default: current file pointer position
!!
!!    status_: optional integer, returns iostat e.g. to check for end of file
!!             on the outside
!!
      implicit none
!
      class(abstract_stream), intent(in) :: the_file
!
      integer, intent(in), optional  :: position_
      integer, intent(out), optional :: status_
!
      real(dp), intent(out) :: scalar
!
      integer ::            io_status = 1
      character(len=100) :: io_message
!
      call the_file%should_be_open('read from')
!
      if(.not. the_file%is_open) then
!
         call output%error_msg('Tried to read from file (a0), which is not open', &
                               chars = [the_file%name_])
      end if
!
      if(present(position_)) then
!
         read(the_file%unit_, pos=position_, iostat=io_status, iomsg=io_message) scalar
!
      else
!
         read(the_file%unit_, iostat=io_status, iomsg=io_message) scalar
!
      end if
!
      call the_file%check_io_status(io_status, trim(io_message), 'read from', status_)
!
   end subroutine read_0_real_dp_abstract_stream
!
!
   subroutine read_0_real_sp_abstract_stream(the_file, scalar, position_, status_)
!!
!!    Read rank 0 real(sp)
!!    Written by Rolf H. Myhre, Feb. 2020
!!
!!    Modified by Alexander C. Paul, Mar. 2020
!!    Check if position_ is present
!!
!!    scalar:     real(sp)
!!
!!    position_:  optional integer, position to read from in file
!!                positions counted in bytes, starting at 1
!!                default: current file pointer position
!!
!!    status_: optional integer, returns iostat e.g. to check for end of file
!!             on the outside
!!
      implicit none
!
      class(abstract_stream), intent(in) :: the_file
!
      integer, intent(in), optional  :: position_
      integer, intent(out), optional :: status_
!
      real(sp), intent(out) :: scalar
!
      integer::             io_status = 1
      character(len=100) :: io_message
!
      call the_file%should_be_open('read from')
!
      if(present(position_)) then
!
         read(the_file%unit_, pos=position_, iostat=io_status, iomsg=io_message) scalar
!
      else
!
         read(the_file%unit_, iostat=io_status, iomsg=io_message) scalar
!
      end if
!
      call the_file%check_io_status(io_status, trim(io_message), 'read from', status_)
!
   end subroutine read_0_real_sp_abstract_stream
!
!
   subroutine read_0_complex_dp_abstract_stream(the_file, scalar, position_, status_)
!!
!!    read rank 0 complex dp
!!    Written by Rolf H. Myhre, Feb. 2020
!!
!!    Modified by Alexander C. Paul, Mar. 2020
!!    Check if position_ is present
!!
!!    scalar: complex double precision scalar
!!
!!    position_:  optional integer, position to read from in file
!!                positions counted in bytes, starting at 1
!!                default: current file pointer position
!!
!!    status_: optional integer, returns iostat e.g. to check for end of file
!!             on the outside
!!
      implicit none
!
      class(abstract_stream), intent(in) :: the_file
!
      integer, intent(in), optional :: position_
      integer, intent(out), optional :: status_
!
      complex(dp), intent(out) :: scalar
!
      integer::             io_status = 1
      character(len=100) :: io_message
!
      call the_file%should_be_open('read from')
!
      if(present(position_)) then
!
         read(the_file%unit_, pos=position_, iostat=io_status, iomsg=io_message) scalar
!
      else
!
         read(the_file%unit_, iostat=io_status, iomsg=io_message) scalar
!
      end if
!
      call the_file%check_io_status(io_status, trim(io_message), 'read from', status_)
!
   end subroutine read_0_complex_dp_abstract_stream
!
!
   subroutine read_0_int_32_abstract_stream(the_file, scalar, position_, status_)
!!
!!    read rank 0 int 32
!!    Written by Rolf H. Myhre, Feb. 2020
!!
!!    Modified by Alexander C. Paul, Mar. 2020
!!    Check if position_ is present
!!
!!    scalar: 32 bit integer
!!
!!    position_:  optional integer, position to read from in file
!!                positions counted in bytes, starting at 1
!!                default: current file pointer position
!!
!!    status_: optional integer, returns iostat e.g. to check for end of file
!!             on the outside
!!
      implicit none
!
      class(abstract_stream), intent(in) :: the_file
!
      integer, intent(in), optional  :: position_
      integer, intent(out), optional :: status_
!
      integer(i32), intent(out) :: scalar
!
      integer::             io_status = 1
      character(len=100) :: io_message
!
      call the_file%should_be_open('read from')
!
      if(present(position_)) then
!
         read(the_file%unit_, pos=position_, iostat=io_status, iomsg=io_message) scalar
!
      else
!
         read(the_file%unit_, iostat=io_status, iomsg=io_message) scalar
!
      end if
!
      call the_file%check_io_status(io_status, trim(io_message), 'read from', status_)
!
   end subroutine read_0_int_32_abstract_stream
!
!
   subroutine read_0_int_64_abstract_stream(the_file, scalar, position_, status_)
!!
!!    read rank 0 int 64
!!    Written by Rolf H. Myhre, Feb. 2020
!!
!!    Modified by Alexander C. Paul, Mar. 2020
!!    Check if position_ is present
!!
!!    scalar: 64 bit integer
!!
!!    position_: optional integer, position to read from in file
!!               positions counted in bytes, starting at 1
!!               default: current file pointer position
!!
!!    status_: optional integer, returns iostat e.g. to check for end of file
!!             on the outside
!!
      implicit none
!
      class(abstract_stream), intent(in) :: the_file
!
      integer, intent(in), optional  :: position_
      integer, intent(out), optional :: status_
!
      integer(i64), intent(out) :: scalar
!
      integer::             io_status = 1
      character(len=100) :: io_message
!
      call the_file%should_be_open('read from')
!
      if(present(position_)) then
!
         read(the_file%unit_, pos=position_, iostat=io_status, iomsg=io_message) scalar
!
      else
!
         read(the_file%unit_, iostat=io_status, iomsg=io_message) scalar
!
      end if
!
      call the_file%check_io_status(io_status, trim(io_message), 'read from', status_)
!
   end subroutine read_0_int_64_abstract_stream
!
!
   subroutine read_0_log_abstract_stream(the_file, scalar, position_, status_)
!!
!!    read rank 0 logical
!!    Written by Rolf H. Myhre and Alexander C. Paul, Mar. 2020
!!
!!    scalar: logical
!!
!!    position_: optional integer, position to read from in file
!!               positions counted in bytes, starting at 1
!!               default: current file pointer position
!!
!!    status_: optional integer, returns iostat e.g. to check for end of file
!!             on the outside
!!
      implicit none
!
      class(abstract_stream), intent(in) :: the_file
!
      integer, intent(in), optional  :: position_
      integer, intent(out), optional :: status_
!
      logical, intent(out) :: scalar
!
      integer::             io_status = 1
      character(len=100) :: io_message
!
      call the_file%should_be_open('read from')
!
      if(present(position_)) then
!
         read(the_file%unit_, pos=position_, iostat=io_status, iomsg=io_message) scalar
!
      else
!
         read(the_file%unit_, iostat=io_status, iomsg=io_message) scalar
!
      end if
!
      call the_file%check_io_status(io_status, trim(io_message), 'read from', status_)
!
   end subroutine read_0_log_abstract_stream
!
!
   subroutine read_1_real_dp_abstract_stream(the_file, array, n, position_, status_)
!!
!!    read rank 1 real dp
!!    Written by Rolf H. Myhre, Feb. 2020
!!
!!    Modified by Alexander C. Paul, Mar. 2020
!!    Check if position_ is present
!!
!!    array: real double precision array of length n
!!
!!    n: length of array
!!
!!    position_:  optional integer, position to read from in file
!!                positions counted in bytes, starting at 1
!!                default: current file pointer position
!!
!!    status_: optional integer, returns iostat e.g. to check for end of file
!!             on the outside
!!
      implicit none
!
      class(abstract_stream), intent(in) :: the_file
!
      integer, intent(in) :: n
      integer, intent(in), optional  :: position_
      integer, intent(out), optional :: status_
!
      real(dp), dimension(n), intent(out) :: array
!
      integer::             io_status = 1
      character(len=100) :: io_message
!
      call the_file%should_be_open('read from')
!
      if(present(position_)) then
!
         read(the_file%unit_, pos=position_, iostat=io_status, iomsg=io_message) array
!
      else
!
         read(the_file%unit_, iostat=io_status, iomsg=io_message) array
!
      end if
!
      call the_file%check_io_status(io_status, trim(io_message), 'read from', status_)
!
   end subroutine read_1_real_dp_abstract_stream
!
!
   subroutine read_1_real_sp_abstract_stream(the_file, array, n, position_, status_)
!!
!!    Read rank 1 real(sp)
!!    Written by Rolf H. Myhre, Feb. 2020
!!
!!    Modified by Alexander C. Paul, Mar. 2020
!!    Check if position_ is present
!!
!!    array: real(sp) array of length n
!!
!!    n: length of array
!!
!!    position_:  optional integer, position to read from in file
!!                positions counted in bytes, starting at 1
!!                default: current file pointer position
!!
!!    status_: optional integer, returns iostat e.g. to check for end of file
!!             on the outside
!!
      implicit none
!
      class(abstract_stream), intent(in) :: the_file
!
      integer, intent(in) :: n
      integer, intent(in), optional  :: position_
      integer, intent(out), optional :: status_
!
      real(sp), dimension(n), intent(out) :: array
!
      integer::             io_status = 1
      character(len=100) :: io_message
!
      call the_file%should_be_open('read from')
!
      if(present(position_)) then
!
         read(the_file%unit_, pos=position_, iostat=io_status, iomsg=io_message) array
!
      else
!
         read(the_file%unit_, iostat=io_status, iomsg=io_message) array
!
      end if
!
      call the_file%check_io_status(io_status, trim(io_message), 'read from', status_)
!
   end subroutine read_1_real_sp_abstract_stream
!
!
   subroutine read_1_complex_dp_abstract_stream(the_file, array, n, position_, status_)
!!
!!    read rank 1 complex dp
!!    Written by Rolf H. Myhre, Feb. 2020
!!
!!    Modified by Alexander C. Paul, Mar. 2020
!!    Check if position_ is present
!!
!!    array: complex double precision array of length n
!!
!!    n: length of array
!!
!!    position_:  optional integer, position to read from in file
!!                positions counted in bytes, starting at 1
!!                default: current file pointer position
!!
!!    status_: optional integer, returns iostat e.g. to check for end of file
!!             on the outside
!!
      implicit none
!
      class(abstract_stream), intent(in) :: the_file
!
      integer, intent(in) :: n
      integer, intent(in), optional  :: position_
      integer, intent(out), optional :: status_
!
      complex(dp), dimension(n), intent(out) :: array
!
      integer::             io_status = 1
      character(len=100) :: io_message
!
      call the_file%should_be_open('read from')
!
      if(present(position_)) then
!
         read(the_file%unit_, pos=position_, iostat=io_status, iomsg=io_message) array
!
      else
!
         read(the_file%unit_, iostat=io_status, iomsg=io_message) array
!
      end if
!
      call the_file%check_io_status(io_status, trim(io_message), 'read from', status_)
!
   end subroutine read_1_complex_dp_abstract_stream
!
!
   subroutine read_1_int_32_abstract_stream(the_file, array, n, position_, status_)
!!
!!    read rank 1 int 32
!!    Written by Rolf H. Myhre, Feb. 2020
!!
!!    Modified by Alexander C. Paul, Mar. 2020
!!    Check if position_ is present
!!
!!    array: 32 bit integer array of length n
!!
!!    n: length of array
!!
!!    position_:  optional integer, position to read from in file
!!                positions counted in bytes, starting at 1
!!                default: current file pointer position
!!
!!    status_: optional integer, returns iostat e.g. to check for end of file
!!             on the outside
!!
      implicit none
!
      class(abstract_stream), intent(in) :: the_file
!
      integer, intent(in) :: n
      integer, intent(in), optional  :: position_
      integer, intent(out), optional :: status_
!
      integer(i32), dimension(n), intent(out) :: array
!
      integer::             io_status = 1
      character(len=100) :: io_message
!
      call the_file%should_be_open('read from')
!
      if(present(position_)) then
!
         read(the_file%unit_, pos=position_, iostat=io_status, iomsg=io_message) array
!
      else
!
         read(the_file%unit_, iostat=io_status, iomsg=io_message) array
!
      end if
!
      call the_file%check_io_status(io_status, trim(io_message), 'read from', status_)
!
   end subroutine read_1_int_32_abstract_stream
!
!
   subroutine read_1_int_64_abstract_stream(the_file, array, n, position_, status_)
!!
!!    read rank 1 int 64
!!    Written by Rolf H. Myhre, Feb. 2020
!!
!!    Modified by Alexander C. Paul, Mar. 2020
!!    Check if position_ is present
!!
!!    array: 64 bit integer array of length n
!!
!!    n: length of array
!!
!!    position_:  optional integer, position to read from in file
!!                positions counted in bytes, starting at 1
!!                default: current file pointer position
!!
!!    status_: optional integer, returns iostat e.g. to check for end of file
!!             on the outside
!!
      implicit none
!
      class(abstract_stream), intent(in) :: the_file
!
      integer, intent(in) :: n
      integer, intent(in), optional  :: position_
      integer, intent(out), optional :: status_
!
      integer(i64), dimension(n), intent(out) :: array
!
      integer::             io_status = 1
      character(len=100) :: io_message
!
      call the_file%should_be_open('read from')
!
      if(present(position_)) then
!
         read(the_file%unit_, pos=position_, iostat=io_status, iomsg=io_message) array
!
      else
!
         read(the_file%unit_, iostat=io_status, iomsg=io_message) array
!
      end if
!
      call the_file%check_io_status(io_status, trim(io_message), 'read from', status_)
!
   end subroutine read_1_int_64_abstract_stream
!
!
   subroutine read_1_log_abstract_stream(the_file, array, n, position_, status_)
!!
!!    read rank 1 logical
!!    Written by Rolf H. Myhre and Alexander C. Paul, Mar. 2020
!!
!!    array: array of logicals with length n
!!
!!    n: length of array
!!
!!    position_:  optional integer, position to read from in file
!!                positions counted in bytes, starting at 1
!!                default: current file pointer position
!!
!!    status_: optional integer, returns iostat e.g. to check for end of file
!!             on the outside
!!
      implicit none
!
      class(abstract_stream), intent(in) :: the_file
!
      integer, intent(in) :: n
      integer, intent(in), optional  :: position_
      integer, intent(out), optional :: status_
!
      logical, dimension(n), intent(out) :: array
!
      integer::             io_status = 1
      character(len=100) :: io_message
!
      call the_file%should_be_open('read from')
!
      if(present(position_)) then
!
         read(the_file%unit_, pos=position_, iostat=io_status, iomsg=io_message) array
!
      else
!
         read(the_file%unit_, iostat=io_status, iomsg=io_message) array
!
      end if
!
      call the_file%check_io_status(io_status, trim(io_message), 'read from', status_)
!
   end subroutine read_1_log_abstract_stream
!
!
   subroutine write_0_real_dp_abstract_stream(the_file, scalar, position_)
!!
!!    write rank 0 real dp
!!    Written by Rolf H. Myhre, Feb. 2020
!!
!!    Modified by Alexander C. Paul, Mar. 2020
!!    Check if position_ is present
!!
!!    scalar: real double precision scalar
!!
!!    position_:  optional integer, position to read from in file
!!                positions counted in bytes, starting at 1
!!                default: current file pointer position
!!
      implicit none
!
      class(abstract_stream), intent(in) :: the_file
!
      integer, intent(in), optional :: position_
!
      real(dp), intent(in) :: scalar
!
      integer::             io_status = 1
      character(len=100) :: io_message
!
      call the_file%should_be_open('write to')
!
      if(present(position_)) then
!
         write(the_file%unit_, pos=position_, iostat=io_status, iomsg=io_message) scalar
!
      else
!
         write(the_file%unit_, iostat=io_status, iomsg=io_message) scalar
!
      end if
!
      call the_file%check_io_status(io_status, trim(io_message), 'write to')
!
   end subroutine write_0_real_dp_abstract_stream
!
!
   subroutine write_0_real_sp_abstract_stream(the_file, scalar, position_)
!!
!!    write rank 0 real(sp)
!!    Written by Rolf H. Myhre, Feb. 2020
!!
!!    Modified by Alexander C. Paul, Mar. 2020
!!    Check if position_ is present
!!
!!    scalar: real(sp) scalar
!!
!!    position_:  optional integer, position to read from in file
!!                positions counted in bytes, starting at 1
!!                default: current file pointer position
!!
      implicit none
!
      class(abstract_stream), intent(in) :: the_file
!
      integer, intent(in), optional :: position_
!
      real(sp), intent(in) :: scalar
!
      integer::             io_status = 1
      character(len=100) :: io_message
!
      call the_file%should_be_open('write to')
!
      if(present(position_)) then
!
         write(the_file%unit_, pos=position_, iostat=io_status, iomsg=io_message) scalar
!
      else
!
         write(the_file%unit_, iostat=io_status, iomsg=io_message) scalar
!
      end if
!
      call the_file%check_io_status(io_status, trim(io_message), 'write to')
!
   end subroutine write_0_real_sp_abstract_stream
!
!
   subroutine write_0_complex_dp_abstract_stream(the_file, scalar, position_)
!!
!!    write rank 0 complex dp
!!    Written by Rolf H. Myhre, Feb. 2020
!!
!!    Modified by Alexander C. Paul, Mar. 2020
!!    Check if position_ is present
!!
!!    scalar: complex double precision scalar
!!
!!    position_:  optional integer, position to read from in file
!!                positions counted in bytes, starting at 1
!!                default: current file pointer position
!!
      implicit none
!
      class(abstract_stream), intent(in) :: the_file
!
      integer, intent(in), optional :: position_
!
      complex(dp), intent(in) :: scalar
!
      integer::             io_status = 1
      character(len=100) :: io_message
!
      call the_file%should_be_open('write to')
!
      if(present(position_)) then
!
         write(the_file%unit_, pos=position_, iostat=io_status, iomsg=io_message) scalar
!
      else
!
         write(the_file%unit_, iostat=io_status, iomsg=io_message) scalar
!
      end if
!
      call the_file%check_io_status(io_status, trim(io_message), 'write to')
!
   end subroutine write_0_complex_dp_abstract_stream
!
!
   subroutine write_0_int_32_abstract_stream(the_file, scalar, position_)
!!
!!    write rank 0 int 32
!!    Written by Rolf H. Myhre, Feb. 2020
!!
!!    Modified by Alexander C. Paul, Mar. 2020
!!    Check if position_ is present
!!
!!    scalar: 32 bit integer
!!
!!    position_:  optional integer, position to read from in file
!!                positions counted in bytes, starting at 1
!!                default: current file pointer position
!!
      implicit none
!
      class(abstract_stream), intent(in) :: the_file
!
      integer, intent(in), optional :: position_
!
      integer(i32), intent(in) :: scalar
!
      integer::             io_status = 1
      character(len=100) :: io_message
!
      call the_file%should_be_open('write to')
!
      if(present(position_)) then
!
         write(the_file%unit_, pos=position_, iostat=io_status, iomsg=io_message) scalar
!
      else
!
         write(the_file%unit_, iostat=io_status, iomsg=io_message) scalar
!
      end if
!
      call the_file%check_io_status(io_status, trim(io_message), 'write to')
!
   end subroutine write_0_int_32_abstract_stream
!
!
   subroutine write_0_int_64_abstract_stream(the_file, scalar, position_)
!!
!!    write rank 0 int 64
!!    Written by Rolf H. Myhre, Feb. 2020
!!
!!    Modified by Alexander C. Paul, Mar. 2020
!!    Check if position_ is present
!!
!!    scalar: 64 bit integer
!!
!!    position_:  optional integer, position to read from in file
!!                positions counted in bytes, starting at 1
!!                default: current file pointer position
!!
      implicit none
!
      class(abstract_stream), intent(in) :: the_file
!
      integer, intent(in), optional :: position_
!
      integer(i64), intent(in) :: scalar
!
      integer::             io_status = 1
      character(len=100) :: io_message
!
      call the_file%should_be_open('write to')
!
      if(present(position_)) then
!
         write(the_file%unit_, pos=position_, iostat=io_status, iomsg=io_message) scalar
!
      else
!
         write(the_file%unit_, iostat=io_status, iomsg=io_message) scalar
!
      end if
!
      call the_file%check_io_status(io_status, trim(io_message), 'write to')
!
   end subroutine write_0_int_64_abstract_stream
!
!
   subroutine write_0_log_abstract_stream(the_file, scalar, position_)
!!
!!    write rank 0 logical
!!    Written by Rolf H. Myhre and Alexander C. Paul, Mar. 2020
!!
!!    scalar: logical
!!
!!    position_:  optional integer, position to read from in file
!!                positions counted in bytes, starting at 1
!!                default: current file pointer position
!!
      implicit none
!
      class(abstract_stream), intent(in) :: the_file
!
      integer, intent(in), optional :: position_
!
      logical, intent(in) :: scalar
!
      integer::             io_status = 1
      character(len=100) :: io_message
!
      call the_file%should_be_open('write to')
!
      if(present(position_)) then
!
         write(the_file%unit_, pos=position_, iostat=io_status, iomsg=io_message) scalar
!
      else
!
         write(the_file%unit_, iostat=io_status, iomsg=io_message) scalar
!
      end if
!
      call the_file%check_io_status(io_status, trim(io_message), 'write to')
!
   end subroutine write_0_log_abstract_stream
!
!
   subroutine write_1_real_dp_abstract_stream(the_file, array, n, position_)
!!
!!    write rank 1 real dp
!!    Written by Rolf H. Myhre, Mar. 2020
!!
!!    Modified by Alexander C. Paul, Mar. 2020
!!    Check if position_ is present
!!
!!    array: real double precision array of length n
!!
!!    n: length of array
!!
!!    position_:  optional integer, position to read from in file
!!                positions counted in bytes, starting at 1
!!                default: current file pointer position
!!
      implicit none
!
      class(abstract_stream), intent(in) :: the_file
!
      integer, intent(in) :: n
      integer, intent(in), optional :: position_
!
      real(dp), dimension(n), intent(in) :: array
!
      integer::             io_status = 1
      character(len=100) :: io_message
!
      call the_file%should_be_open('write to')
!
      if(present(position_)) then
!
         write(the_file%unit_, pos=position_, iostat=io_status, iomsg=io_message) array
!
      else
!
         write(the_file%unit_, iostat=io_status, iomsg=io_message) array
!
      end if
!
      call the_file%check_io_status(io_status, trim(io_message), 'write to')
!
   end subroutine write_1_real_dp_abstract_stream
!
!
   subroutine write_1_real_sp_abstract_stream(the_file, array, n, position_)
!!
!!    write rank 1 real(sp)
!!    Written by Rolf H. Myhre, Mar. 2020
!!
!!    Modified by Alexander C. Paul, Mar. 2020
!!    Check if position_ is present
!!
!!    array: real(sp) array of length n
!!
!!    n: length of array
!!
!!    position_:  optional integer, position to read from in file
!!                positions counted in bytes, starting at 1
!!                default: current file pointer position
!!
      implicit none
!
      class(abstract_stream), intent(in) :: the_file
!
      integer, intent(in) :: n
      integer, intent(in), optional :: position_
!
      real(sp), dimension(n), intent(in) :: array
!
      integer::             io_status = 1
      character(len=100) :: io_message
!
      call the_file%should_be_open('write to')
!
      if(present(position_)) then
!
         write(the_file%unit_, pos=position_, iostat=io_status, iomsg=io_message) array
!
      else
!
         write(the_file%unit_, iostat=io_status, iomsg=io_message) array
!
      end if
!
      call the_file%check_io_status(io_status, trim(io_message), 'write to')
!
   end subroutine write_1_real_sp_abstract_stream
!
!
   subroutine write_1_complex_dp_abstract_stream(the_file, array, n, position_)
!!
!!    write rank 1 complex dp
!!    Written by Rolf H. Myhre, Mar. 2020
!!
!!    Modified by Alexander C. Paul, Mar. 2020
!!    Check if position_ is present
!!
!!    array: complex double precision array of length n
!!
!!    n: length of array
!!
!!    position_:  optional integer, position to read from in file
!!                positions counted in bytes, starting at 1
!!                default: current file pointer position
!!
      implicit none
!
      class(abstract_stream), intent(in) :: the_file
!
      integer, intent(in) :: n
      integer, intent(in), optional :: position_
!
      complex(dp), dimension(n), intent(in) :: array
!
      integer::             io_status = 1
      character(len=100) :: io_message
!
      call the_file%should_be_open('write to')
!
      if(present(position_)) then
!
         write(the_file%unit_, pos=position_, iostat=io_status, iomsg=io_message) array
!
      else
!
         write(the_file%unit_, iostat=io_status, iomsg=io_message) array
!
      end if
!
      call the_file%check_io_status(io_status, trim(io_message), 'write to')
!
   end subroutine write_1_complex_dp_abstract_stream
!
!
   subroutine write_1_int_32_abstract_stream(the_file, array, n, position_)
!!
!!    write rank 1 int 32
!!    Written by Rolf H. Myhre, Mar. 2020
!!
!!    Modified by Alexander C. Paul, Mar. 2020
!!    Check if position_ is present
!!
!!    array: 32 bit integer array of length n
!!
!!    n: length of array
!!
!!    position_:  optional integer, position to read from in file
!!                positions counted in bytes, starting at 1
!!                default: current file pointer position
!!
      implicit none
!
      class(abstract_stream), intent(in) :: the_file
!
      integer, intent(in) :: n
      integer, intent(in), optional :: position_
!
      integer(i32), dimension(n), intent(in) :: array
!
      integer::             io_status = 1
      character(len=100) :: io_message
!
      call the_file%should_be_open('write to')
!
      if(present(position_)) then
!
         write(the_file%unit_, pos=position_, iostat=io_status, iomsg=io_message) array
!
      else
!
         write(the_file%unit_, iostat=io_status, iomsg=io_message) array
!
      end if
!
      call the_file%check_io_status(io_status, trim(io_message), 'write to')
!
   end subroutine write_1_int_32_abstract_stream
!
!
   subroutine write_1_int_64_abstract_stream(the_file, array, n, position_)
!!
!!    write rank 1 int 64
!!    Written by Rolf H. Myhre, Mar. 2020
!!
!!    Modified by Alexander C. Paul, Mar. 2020
!!    Check if position_ is present
!!
!!    array: 64 bit integer array of length n
!!
!!    n: length of array
!!
!!    position_:  optional integer, position to read from in file
!!                positions counted in bytes, starting at 1
!!                default: current file pointer position
!!
      implicit none
!
      class(abstract_stream), intent(in) :: the_file
!
      integer, intent(in) :: n
      integer, intent(in), optional :: position_
!
      integer(i64), dimension(n), intent(in) :: array
!
      integer::             io_status = 1
      character(len=100) :: io_message
!
      call the_file%should_be_open('write to')
!
      if(present(position_)) then
!
         write(the_file%unit_, pos=position_, iostat=io_status, iomsg=io_message) array
!
      else
!
         write(the_file%unit_, iostat=io_status, iomsg=io_message) array
!
      end if
!
      call the_file%check_io_status(io_status, trim(io_message), 'write to')
!
   end subroutine write_1_int_64_abstract_stream
!
!
   subroutine write_1_log_abstract_stream(the_file, array, n, position_)
!!
!!    write rank 1 logical
!!    Written by Rolf H. Myhre and Alexander C. Paul, Mar. 2020
!!
!!    array: 64 bit integer array of length n
!!
!!    n: length of array
!!
!!    position_:  optional integer, position to read from in file
!!                positions counted in bytes, starting at 1
!!                default: current file pointer position
!!
      implicit none
!
      class(abstract_stream), intent(in) :: the_file
!
      integer, intent(in) :: n
      integer, intent(in), optional :: position_
!
      logical, dimension(n), intent(in) :: array
!
      integer::             io_status = 1
      character(len=100) :: io_message
!
      call the_file%should_be_open('write to')
!
      if(present(position_)) then
!
         write(the_file%unit_, pos=position_, iostat=io_status, iomsg=io_message) array
!
      else
!
         write(the_file%unit_, iostat=io_status, iomsg=io_message) array
!
      end if
!
      call the_file%check_io_status(io_status, trim(io_message), 'write to')
!
   end subroutine write_1_log_abstract_stream
!
!
   subroutine copy_abstract_stream(the_file, filename)
!!
!!    Copy abstract file
!!    Written by Alexander C. Paul and Rolf H. Myhre, September 2019
!!
!!    This routine is slow and should not be used for files of nontrivial size.
!!
      implicit none
!
      class(abstract_stream) :: the_file
!
      character(*), intent(in) :: filename
!
      integer              :: io_status
      character(len = 100) :: io_message
!
      logical :: was_closed
!
      integer :: copy_unit
!
!     One byte character
      integer(i8) :: byte
!
!     Open a new file
      open(newunit=copy_unit, file=trim(filename), access='stream', &
           form='unformatted', action='write', status='new', &
           iostat=io_status, iomsg=io_message)
!
      call the_file%check_io_status(io_status, trim(io_message), 'open')
!
!     Check if the file is open
!
      if(the_file%is_open) then
!
!        Rewind to start
         read(the_file%unit_, pos=1)
         was_closed = .false.
!
      else
!
!        Open the file for reading from start
         call the_file%open_('read', 'rewind')
         was_closed = .true.
!
      endif
!
!     Read byte by byte and write it to the new file
!
      do
!
         read(the_file%unit_, end=200) byte !Read until end of file, then go to 200
         write(copy_unit) byte             !Write whatever you just read
!
      enddo
      200 continue !End of file reached, should be done
!
!     Close the files
      close(copy_unit, status='keep', iostat=io_status, iomsg=io_message)
!
      call the_file%check_io_status(io_status, trim(io_message), 'close')
!
      if(was_closed) then
         call the_file%close_()
      end if
!
!
   end subroutine copy_abstract_stream
!
!
   subroutine delete_abstract_stream(the_file)
!!
!!    Delete file
!!    Written by Rolf H. Myhre, Aug 2019
!!
      implicit none
!
      class(abstract_stream) :: the_file
!
      integer            :: io_status
      character(len=100) :: io_message
!
      if(the_file%is_open) then
!
         close(the_file%unit_, iostat=io_status, iomsg=io_message, status='delete')
!
         call the_file%check_io_status(io_status, trim(io_message), 'delete')
!
      else
!
         call the_file%open_('write')
!
         close(the_file%unit_, iostat=io_status, iomsg=io_message, status='delete')
!
         call the_file%check_io_status(io_status, trim(io_message), 'delete')
!
      endif
!
      the_file%is_open = .false.
      the_file%status_ = 'new'
!
   end subroutine delete_abstract_stream
!
!
   subroutine check_io_status_abstract_stream(the_file, io_status, io_message, task, status_)
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
      implicit none
!
      class(abstract_stream), intent(in) :: the_file
!
      integer, intent(in) :: io_status
!
      character(len=*), intent(in) :: io_message
      character(len=*), intent(in) :: task
!
      integer, intent(out), optional :: status_
!
      if(present(status_) .and. io_status .le. 0) then
!
         status_ = io_status
!
      else if(io_status .ne. 0) then
!
         call output%error_msg('Failed to '// task // ' file (a0), status is &
                               &(i0) and error message is ' // io_message,   &
                               chars = [the_file%name_], ints = [io_status])
      endif
!
   end subroutine check_io_status_abstract_stream
!
!
   subroutine should_be_open_abstract_stream(the_file, task)
!!
!!    Should be open
!!    Written by Alexander C. Paul, Oct 2020
!!
!!    Checks if file is open and returns error if not
!!
!!    task: task for which an open file is expected (close, read, write)
!!
      implicit none
!
      class(abstract_stream), intent(in) :: the_file
!
      character(len=*), intent(in) :: task
!
      if(.not. the_file%is_open) then
!
         call output%error_msg('Trying to ' // task // ' file (a0), which is closed.', &
                                chars=[the_file%name_])
      end if
!
   end subroutine should_be_open_abstract_stream
!
!
   subroutine should_be_closed_abstract_stream(the_file, task)
!!
!!    Should be closed
!!    Written by Alexander C. Paul, Oct 2020
!!
!!    Checks if file is open and returns error if not
!!
!!    task: task for which an open file is expected (close, read, write)
!!
      implicit none
!
      class(abstract_stream), intent(in) :: the_file
!
      character(len=*), intent(in) :: task
!
      if(the_file%is_open) then
!
         call output%error_msg('Trying to ' // task // ' file (a0), which is open.', &
                                chars=[the_file%name_])
      end if
!
   end subroutine should_be_closed_abstract_stream
!
!
end module abstract_stream_class
