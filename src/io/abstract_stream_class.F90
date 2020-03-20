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
      procedure :: set_name                => set_name_abstract_stream
      procedure :: get_name                => get_name_abstract_stream
!
      procedure :: set_status              => set_status_abstract_stream
!
      procedure :: get_open                => get_open_abstract_stream
!
      procedure :: open_                   => open_abstract_stream
      procedure :: close_                  => close_abstract_stream
      procedure :: delete_                 => delete_abstract_stream
!
      procedure :: exists                  => exists_abstract_stream
      procedure :: get_file_size           => get_file_size_abstract_stream
      procedure :: copy                    => copy_abstract_stream
!
      procedure :: read_real_dp_scalar     => read_real_dp_scalar_abstract_stream
      procedure :: read_complex_dp_scalar  => read_complex_dp_scalar_abstract_stream
      procedure :: read_int_32_dp_scalar   => read_int_32_scalar_abstract_stream
      procedure :: read_int_64_dp_scalar   => read_int_64_scalar_abstract_stream
!
      procedure :: read_real_dp            => read_real_dp_abstract_stream
      procedure :: read_complex_dp         => read_complex_dp_abstract_stream
      procedure :: read_int_32             => read_int_32_abstract_stream
      procedure :: read_int_64             => read_int_64_abstract_stream
!
      procedure :: write_real_dp_scalar    => write_real_dp_scalar_abstract_stream
      procedure :: write_complex_dp_scalar => write_complex_dp_scalar_abstract_stream
      procedure :: write_int_32_dp_scalar  => write_int_32_scalar_abstract_stream
      procedure :: write_int_64_dp_scalar  => write_int_64_scalar_abstract_stream
!
      procedure :: write_real_dp           => write_real_dp_abstract_stream
      procedure :: write_complex_dp        => write_complex_dp_abstract_stream
      procedure :: write_int_32            => write_int_32_abstract_stream
      procedure :: write_int_64            => write_int_64_abstract_stream
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
      if(the_file%is_open) then
!
         call output%error_msg('Trying to open (a0), but it is already open.', &
                                chars=[the_file%name_])
      end if
!
!
!     Check if new_action or new_position is present
!     We let open handle the error if they is messed up
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
!
      open(newunit = the_file%unit_, file=trim(the_file%name_), access = 'stream', &
           form = 'unformatted', action = action_, status = the_file%status_, &
           position = position_, iostat = io_status, iomsg = io_message)
!
!
      if(io_status .ne. 0) then
!
         call output%error_msg('Failed to open file (a0), status is (i0) and &
                               &error message is (a0)', chars = [the_file%name_, io_message], &
                               ints = [io_status])
      end if
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
!
!     Check if the file is open
!
      if(.not. the_file%is_open) then
         call output%error_msg('Trying to close (a0), but it is already closed.', &
                                chars=[the_file%name_])
      end if
!
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
!
!     Things seems fine, close the file
!
      close(the_file%unit_, status = destiny, iostat = io_status, iomsg = io_message)
!
!
!     Was it fine?
!
      if(io_status .ne. 0) then
!
         call output%error_msg('Failed to close file (a0), status is (i0) and &
                               &error message is (a0)', chars = [the_file%name_, io_message], &
                               ints = [io_status])
      end if
!
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
      inquire(file=trim(the_file%name_), size=file_size )
!
   end function get_file_size_abstract_stream
!
!
   subroutine read_real_dp_scalar_abstract_stream(the_file, scalar, position_)
!!
!!    Read real dp scalar
!!    Written by Rolf H. Myhre, Feb. 2020
!!
!!    scalar: real double precision scalar
!!
!!    pos_:  optional integer, position to read from in file
!!           positions counted in bytes, starting at 1
!!           default: current file pointer position
!!
      implicit none
!
      class(abstract_stream), intent(in) :: the_file
!
      integer, intent(in), optional :: position_
!
      real(dp), intent(out) :: scalar
!
      integer::             io_status = 1
      character(len=100) :: io_message
!
      if(.not. the_file%is_open) then
!
         call output%error_msg('Tried to read from file (a0), which is not open', &
                               chars = [the_file%name_])
      end if
!
      read(the_file%unit_, pos=position_, iostat=io_status, iomsg=io_message) scalar
!
      if(io_status .ne. 0) then
!
         call output%error_msg('Failed to read from file (a0), status is (i0) and &
                               &error message is (a0)', chars = [the_file%name_, io_message], &
                               ints = [io_status])
      end if
!
   end subroutine read_real_dp_scalar_abstract_stream
!
!
   subroutine read_complex_dp_scalar_abstract_stream(the_file, scalar, position_)
!!
!!    Read complex dp scalar
!!    Written by Rolf H. Myhre, Feb. 2020
!!
!!    scalar: complex double precision scalar
!!
!!    pos_:  optional integer, position to read from in file
!!           positions counted in bytes, starting at 1
!!           default: current file pointer position
!!
      implicit none
!
      class(abstract_stream), intent(in) :: the_file
!
      integer, intent(in), optional :: position_
!
      complex(dp), intent(out) :: scalar
!
      integer::             io_status = 1
      character(len=100) :: io_message
!
      if(.not. the_file%is_open) then
!
         call output%error_msg('Tried to read from file (a0), which is not open', &
                               chars = [the_file%name_])
      end if
!
      read(the_file%unit_, pos=position_, iostat=io_status, iomsg=io_message) scalar
!
      if(io_status .ne. 0) then
!
         call output%error_msg('Failed to read from file (a0), status is (i0) and &
                               &error message is (a0)', chars = [the_file%name_, io_message], &
                               ints = [io_status])
      end if
!
   end subroutine read_complex_dp_scalar_abstract_stream
!
!
   subroutine read_int_32_scalar_abstract_stream(the_file, scalar, position_)
!!
!!    Read int 32 scalar
!!    Written by Rolf H. Myhre, Feb. 2020
!!
!!    scalar: real 32 bit integer
!!
!!    pos_:  optional integer, position to read from in file
!!           positions counted in bytes, starting at 1
!!           default: current file pointer position
!!
      implicit none
!
      class(abstract_stream), intent(in) :: the_file
!
      integer, intent(in), optional :: position_
!
      integer(i32), intent(out) :: scalar
!
      integer::             io_status = 1
      character(len=100) :: io_message
!
      if(.not. the_file%is_open) then
!
         call output%error_msg('Tried to read from file (a0), which is not open', &
                               chars = [the_file%name_])
      end if
!
      read(the_file%unit_, pos=position_, iostat=io_status, iomsg=io_message) scalar
!
      if(io_status .ne. 0) then
!
         call output%error_msg('Failed to read from file (a0), status is (i0) and &
                               &error message is (a0)', chars = [the_file%name_, io_message], &
                               ints = [io_status])
      end if
!
   end subroutine read_int_32_scalar_abstract_stream
!
!
   subroutine read_int_64_scalar_abstract_stream(the_file, scalar, position_)
!!
!!    Read int 64 scalar
!!    Written by Rolf H. Myhre, Feb. 2020
!!
!!    scalar: real 64 bit integer
!!
!!    pos_:  optional integer, position to read from in file
!!           positions counted in bytes, starting at 1
!!           default: current file pointer position
!!
      implicit none
!
      class(abstract_stream), intent(in) :: the_file
!
      integer, intent(in), optional :: position_
!
      integer(i64), intent(out) :: scalar
!
      integer::             io_status = 1
      character(len=100) :: io_message
!
      if(.not. the_file%is_open) then
!
         call output%error_msg('Tried to read from file (a0), which is not open', &
                               chars = [the_file%name_])
      end if
!
      read(the_file%unit_, pos=position_, iostat=io_status, iomsg=io_message) scalar
!
      if(io_status .ne. 0) then
!
         call output%error_msg('Failed to read from file (a0), status is (i0) and &
                               &error message is (a0)', chars = [the_file%name_, io_message], &
                               ints = [io_status])
      end if
!
   end subroutine read_int_64_scalar_abstract_stream
!
!
   subroutine read_real_dp_abstract_stream(the_file, array, n, position_)
!!
!!    Read real dp
!!    Written by Rolf H. Myhre, Feb. 2020
!!
!!    array: real double precision array of length n
!!
!!    n: length of array
!!
!!    pos_:  optional integer, position to read from in file
!!           positions counted in bytes, starting at 1
!!           default: current file pointer position
!!
      implicit none
!
      class(abstract_stream), intent(in) :: the_file
!
      integer, intent(in) :: n
      integer, intent(in), optional :: position_
!
      real(dp), dimension(n), intent(out) :: array
!
      integer::             io_status = 1
      character(len=100) :: io_message
!
      if(.not. the_file%is_open) then
!
         call output%error_msg('Tried to read from file (a0), which is not open', &
                               chars = [the_file%name_])
      end if
!
      read(the_file%unit_, pos=position_, iostat=io_status, iomsg=io_message) array
!
      if(io_status .ne. 0) then
!
         call output%error_msg('Failed to read from file (a0), status is (i0) and &
                               &error message is (a0)', chars = [the_file%name_, io_message], &
                               ints = [io_status])
      end if
!
   end subroutine read_real_dp_abstract_stream
!
!
   subroutine read_complex_dp_abstract_stream(the_file, array, n, position_)
!!
!!    Read complex dp
!!    Written by Rolf H. Myhre, Feb. 2020
!!
!!    array: complex double precision array of length n
!!
!!    n: length of array
!!
!!    pos_:  optional integer, position to read from in file
!!           positions counted in bytes, starting at 1
!!           default: current file pointer position
!!
      implicit none
!
      class(abstract_stream), intent(in) :: the_file
!
      integer, intent(in) :: n
      integer, intent(in), optional :: position_
!
      complex(dp), dimension(n), intent(out) :: array
!
      integer::             io_status = 1
      character(len=100) :: io_message
!
      if(.not. the_file%is_open) then
!
         call output%error_msg('Tried to read from file (a0), which is not open', &
                               chars = [the_file%name_])
      end if
!
      read(the_file%unit_, pos=position_, iostat=io_status, iomsg=io_message) array
!
      if(io_status .ne. 0) then
!
         call output%error_msg('Failed to read from file (a0), status is (i0) and &
                               &error message is (a0)', chars = [the_file%name_, io_message], &
                               ints = [io_status])
      end if
!
   end subroutine read_complex_dp_abstract_stream
!
!
   subroutine read_int_32_abstract_stream(the_file, array, n, position_)
!!
!!    Read int 32
!!    Written by Rolf H. Myhre, Feb. 2020
!!
!!    array: 32 bit integer array of length n
!!
!!    n: length of array
!!
!!    pos_:  optional integer, position to read from in file
!!           positions counted in bytes, starting at 1
!!           default: current file pointer position
!!
      implicit none
!
      class(abstract_stream), intent(in) :: the_file
!
      integer, intent(in) :: n
      integer, intent(in), optional :: position_
!
      integer(i32), dimension(n), intent(out) :: array
!
      integer::             io_status = 1
      character(len=100) :: io_message
!
      if(.not. the_file%is_open) then
!
         call output%error_msg('Tried to read from file (a0), which is not open', &
                               chars = [the_file%name_])
      end if
!
      read(the_file%unit_, pos=position_, iostat=io_status, iomsg=io_message) array
!
      if(io_status .ne. 0) then
!
         call output%error_msg('Failed to read from file (a0), status is (i0) and &
                               &error message is (a0)', chars = [the_file%name_, io_message], &
                               ints = [io_status])
      end if
!
   end subroutine read_int_32_abstract_stream
!
!
   subroutine read_int_64_abstract_stream(the_file, array, n, position_)
!!
!!    Read int 64
!!    Written by Rolf H. Myhre, Feb. 2020
!!
!!    array: 64 bit integer array of length n
!!
!!    n: length of array
!!
!!    pos_:  optional integer, position to read from in file
!!           positions counted in bytes, starting at 1
!!           default: current file pointer position
!!
      implicit none
!
      class(abstract_stream), intent(in) :: the_file
!
      integer, intent(in) :: n
      integer, intent(in), optional :: position_
!
      integer(i64), dimension(n), intent(out) :: array
!
      integer::             io_status = 1
      character(len=100) :: io_message
!
      if(.not. the_file%is_open) then
!
         call output%error_msg('Tried to read from file (a0), which is not open', &
                               chars = [the_file%name_])
      end if
!
      read(the_file%unit_, pos=position_, iostat=io_status, iomsg=io_message) array
!
      if(io_status .ne. 0) then
!
         call output%error_msg('Failed to read from file (a0), status is (i0) and &
                               &error message is (a0)', chars = [the_file%name_, io_message], &
                               ints = [io_status])
      end if
!
   end subroutine read_int_64_abstract_stream
!
!
   subroutine write_real_dp_scalar_abstract_stream(the_file, scalar, position_)
!!
!!    write real dp scalar
!!    Written by Rolf H. Myhre, Feb. 2020
!!
!!    scalar: real double precision scalar
!!
!!    pos_:  optional integer, position to write from in file
!!           positions counted in bytes, starting at 1
!!           default: current file pointer position
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
      if(.not. the_file%is_open) then
!
         call output%error_msg('Tried to write from file (a0), which is not open', &
                               chars = [the_file%name_])
      end if
!
      write(the_file%unit_, pos=position_, iostat=io_status, iomsg=io_message) scalar
!
      if(io_status .ne. 0) then
!
         call output%error_msg('Failed to write from file (a0), status is (i0) and &
                               &error message is (a0)', chars = [the_file%name_, io_message], &
                               ints = [io_status])
      end if
!
   end subroutine write_real_dp_scalar_abstract_stream
!
!
   subroutine write_complex_dp_scalar_abstract_stream(the_file, scalar, position_)
!!
!!    write complex dp scalar
!!    Written by Rolf H. Myhre, Feb. 2020
!!
!!    scalar: complex double precision scalar
!!
!!    pos_:  optional integer, position to write from in file
!!           positions counted in bytes, starting at 1
!!           default: current file pointer position
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
      if(.not. the_file%is_open) then
!
         call output%error_msg('Tried to write from file (a0), which is not open', &
                               chars = [the_file%name_])
      end if
!
      write(the_file%unit_, pos=position_, iostat=io_status, iomsg=io_message) scalar
!
      if(io_status .ne. 0) then
!
         call output%error_msg('Failed to write from file (a0), status is (i0) and &
                               &error message is (a0)', chars = [the_file%name_, io_message], &
                               ints = [io_status])
      end if
!
   end subroutine write_complex_dp_scalar_abstract_stream
!
!
   subroutine write_int_32_scalar_abstract_stream(the_file, scalar, position_)
!!
!!    write int 32 scalar
!!    Written by Rolf H. Myhre, Feb. 2020
!!
!!    scalar: real 32 bit integer
!!
!!    pos_:  optional integer, position to write from in file
!!           positions counted in bytes, starting at 1
!!           default: current file pointer position
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
      if(.not. the_file%is_open) then
!
         call output%error_msg('Tried to write from file (a0), which is not open', &
                               chars = [the_file%name_])
      end if
!
      write(the_file%unit_, pos=position_, iostat=io_status, iomsg=io_message) scalar
!
      if(io_status .ne. 0) then
!
         call output%error_msg('Failed to write from file (a0), status is (i0) and &
                               &error message is (a0)', chars = [the_file%name_, io_message], &
                               ints = [io_status])
      end if
!
   end subroutine write_int_32_scalar_abstract_stream
!
!
   subroutine write_int_64_scalar_abstract_stream(the_file, scalar, position_)
!!
!!    write int 64 scalar
!!    Written by Rolf H. Myhre, Feb. 2020
!!
!!    scalar: real 64 bit integer
!!
!!    pos_:  optional integer, position to write from in file
!!           positions counted in bytes, starting at 1
!!           default: current file pointer position
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
      if(.not. the_file%is_open) then
!
         call output%error_msg('Tried to write from file (a0), which is not open', &
                               chars = [the_file%name_])
      end if
!
      write(the_file%unit_, pos=position_, iostat=io_status, iomsg=io_message) scalar
!
      if(io_status .ne. 0) then
!
         call output%error_msg('Failed to write from file (a0), status is (i0) and &
                               &error message is (a0)', chars = [the_file%name_, io_message], &
                               ints = [io_status])
      end if
!
   end subroutine write_int_64_scalar_abstract_stream
!
!
   subroutine write_real_dp_abstract_stream(the_file, array, n, position_)
!!
!!    Write real dp
!!    Written by Rolf H. Myhre, Mar. 2020
!!
!!    array: real double precision array of length n
!!
!!    n: length of array
!!
!!    pos_:  optional integer, position to write to in file
!!           positions counted in bytes, starting at 1
!!           default: current file pointer position
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
      if(.not. the_file%is_open) then
!
         call output%error_msg('Tried to write to file (a0), which is not open', &
                               chars = [the_file%name_])
      end if
!
      write(the_file%unit_, pos=position_, iostat=io_status, iomsg=io_message) array
!
      if(io_status .ne. 0) then
!
         call output%error_msg('Failed to write to file (a0), status is (i0) and &
                               &error message is (a0)', chars = [the_file%name_, io_message], &
                               ints = [io_status])
      end if
!
   end subroutine write_real_dp_abstract_stream
!
!
   subroutine write_complex_dp_abstract_stream(the_file, array, n, position_)
!!
!!    Write complex dp
!!    Written by Rolf H. Myhre, Mar. 2020
!!
!!    array: complex double precision array of length n
!!
!!    n: length of array
!!
!!    pos_:  optional integer, position to write to in file
!!           positions counted in bytes, starting at 1
!!           default: current file pointer position
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
      if(.not. the_file%is_open) then
!
         call output%error_msg('Tried to write to file (a0), which is not open', &
                               chars = [the_file%name_])
      end if
!
      write(the_file%unit_, pos=position_, iostat=io_status, iomsg=io_message) array
!
      if(io_status .ne. 0) then
!
         call output%error_msg('Failed to write to file (a0), status is (i0) and &
                               &error message is (a0)', chars = [the_file%name_, io_message], &
                               ints = [io_status])
      end if
!
   end subroutine write_complex_dp_abstract_stream
!
!
   subroutine write_int_32_abstract_stream(the_file, array, n, position_)
!!
!!    Write int 32
!!    Written by Rolf H. Myhre, Mar. 2020
!!
!!    array: 32 bit integer array of length n
!!
!!    n: length of array
!!
!!    pos_:  optional integer, position to write to in file
!!           positions counted in bytes, starting at 1
!!           default: current file pointer position
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
      if(.not. the_file%is_open) then
!
         call output%error_msg('Tried to write to file (a0), which is not open', &
                               chars = [the_file%name_])
      end if
!
      write(the_file%unit_, pos=position_, iostat=io_status, iomsg=io_message) array
!
      if(io_status .ne. 0) then
!
         call output%error_msg('Failed to write to file (a0), status is (i0) and &
                               &error message is (a0)', chars = [the_file%name_, io_message], &
                               ints = [io_status])
      end if
!
   end subroutine write_int_32_abstract_stream
!
!
   subroutine write_int_64_abstract_stream(the_file, array, n, position_)
!!
!!    Write int 64
!!    Written by Rolf H. Myhre, Mar. 2020
!!
!!    array: 64 bit integer array of length n
!!
!!    n: length of array
!!
!!    pos_:  optional integer, position to write to in file
!!           positions counted in bytes, starting at 1
!!           default: current file pointer position
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
      if(.not. the_file%is_open) then
!
         call output%error_msg('Tried to write to file (a0), which is not open', &
                               chars = [the_file%name_])
      end if
!
      write(the_file%unit_, pos=position_, iostat=io_status, iomsg=io_message) array
!
      if(io_status .ne. 0) then
!
         call output%error_msg('Failed to write to file (a0), status is (i0) and &
                               &error message is (a0)', chars = [the_file%name_, io_message], &
                               ints = [io_status])
      end if
!
   end subroutine write_int_64_abstract_stream
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
      if(io_status .ne. 0) then
!
         call output%error_msg('Failed to open file (a0), status is (i0) and &
                               &error message is (a0)', chars = [filename, io_message], &
                               ints = [io_status])
      end if
!
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
      if(io_status .ne. 0) then
!
         call output%error_msg('Failed to close file (a0), status is (i0) and &
                               &error message is (a0)', chars = [filename, io_message], &
                               ints = [io_status])
      end if
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
      integer              :: io_error
      character(len=100)   :: io_msg
!
      if(the_file%is_open) then
!
         close(the_file%unit_, iostat=io_error, iomsg=io_msg, status='delete')
!
         if (io_error .ne. 0) then
            call output%error_msg('could not delete eT file '//trim(the_file%name_)//&
                                 &'. Error message: '//trim(io_msg))
         endif
!
      else
!
!        Open the file with unformatted stream access
!        Stream doesn't care about format and records and neither do we
         open(newunit=the_file%unit_, file=the_file%name_, access='stream', &
              form='unformatted', action='write', status='old', &
              iostat=io_error, iomsg=io_msg)
!
         if (io_error .ne. 0) then
            call output%error_msg('could not delete eT file '//trim(the_file%name_)//&
                                 &'. Error message: '//trim(io_msg))
         endif
!
         close(the_file%unit_, iostat=io_error, iomsg=io_msg, status='delete')
!
         if (io_error .ne. 0) then
            call output%error_msg('could not delete eT file '//trim(the_file%name_)//&
                                 &'. Error message: '//trim(io_msg))
         endif
!
      endif
!
      the_file%is_open = .false.
!
   end subroutine delete_abstract_stream
!
!
end module abstract_stream_class
