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
module abstract_stream_class
!
!!
!!    Abstract stream class module
!!    Written by Rolf H. Myhre and Alexander C. Paul, 2020-2022
!!
!
   use kinds
   use global_out, only: output
   use abstract_file_class, only: abstract_file
!
   type, abstract, extends(abstract_file) :: abstract_stream
!
   contains
!
      procedure, public :: initialize &
                        => initialize_abstract_stream
!
      procedure, public :: get_file_size &
                        => get_file_size_abstract_stream
!
      procedure, public :: read_0_real_dp_abstract_stream
      procedure, public :: read_0_int_64_abstract_stream
!
      procedure, public :: read_1_real_dp_abstract_stream
      procedure, public :: read_1_complex_dp_abstract_stream
      procedure, public :: read_1_int_32_abstract_stream
      procedure, public :: read_1_int_64_abstract_stream
      procedure, public :: read_1_log_abstract_stream
!
      procedure, public :: write_0_real_dp_abstract_stream
      procedure, public :: write_0_int_32_abstract_stream
      procedure, public :: write_0_int_64_abstract_stream
!
      procedure, public :: write_1_real_dp_abstract_stream
      procedure, public :: write_1_real_sp_abstract_stream
      procedure, public :: write_1_complex_dp_abstract_stream
      procedure, public :: write_1_int_32_abstract_stream
      procedure, public :: write_1_int_64_abstract_stream
      procedure, public :: write_1_log_abstract_stream
!
      procedure, public :: io_error &
                        => io_error_abstract_stream
!
   end type abstract_stream
!
!
contains
!
!
   subroutine initialize_abstract_stream(this, name_, new_status)
!!
!!    Initialize
!!    Written by Alexander C. Paul, Sep. 2022
!!
!!    new_action: Optional, read, write or readwrite (default)
!!    new_status: Optional, old, new, unkown (default)
!!
      implicit none
!
      class(abstract_stream), intent(inout) :: this
!
      character(len=*), intent(in) :: name_
      character(len=*), intent(in), optional :: new_status
!
      this%name_ = trim(name_)
      this%access_ = 'stream'
      this%format_ = 'unformatted'
!
      this%action_ = 'readwrite'
      this%status_ = 'unknown'
      if (present(new_status)) this%status_ = new_status
!
   end subroutine initialize_abstract_stream
!
!
   function get_file_size_abstract_stream(this) result(file_size)
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
      class(abstract_stream), intent(inout) :: this
!
      integer :: file_size
!
      inquire(file=this%get_name(), size=file_size)
!
   end function get_file_size_abstract_stream
!
!
   subroutine read_0_real_dp_abstract_stream(this, scalar, position_, status_)
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
      class(abstract_stream), intent(in) :: this
!
      integer, intent(in), optional  :: position_
      integer, intent(out), optional :: status_
!
      real(dp), intent(out) :: scalar
!
      integer ::            io_status = 1
      character(len=100) :: io_message
!
      call this%error_if_closed(task='read')
!
      if(present(position_)) then
!
         read(this%unit_, pos=position_, iostat=io_status, iomsg=io_message) scalar
!
      else
!
         read(this%unit_, iostat=io_status, iomsg=io_message) scalar
!
      end if
!
      call this%check_io_status(io_status, io_message, task='read', status_=status_)
!
   end subroutine read_0_real_dp_abstract_stream
!
!
   subroutine read_0_int_64_abstract_stream(this, scalar, position_, status_)
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
      class(abstract_stream), intent(in) :: this
!
      integer, intent(in), optional  :: position_
      integer, intent(out), optional :: status_
!
      integer(i64), intent(out) :: scalar
!
      integer::             io_status = 1
      character(len=100) :: io_message
!
      call this%error_if_closed(task='read')
!
      if(present(position_)) then
!
         read(this%unit_, pos=position_, iostat=io_status, iomsg=io_message) scalar
!
      else
!
         read(this%unit_, iostat=io_status, iomsg=io_message) scalar
!
      end if
!
      call this%check_io_status(io_status, io_message, task='read', status_=status_)
!
   end subroutine read_0_int_64_abstract_stream
!
!
   subroutine read_1_real_dp_abstract_stream(this, array, n, position_, status_)
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
      class(abstract_stream), intent(in) :: this
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
      call this%error_if_closed(task='read')
!
      if(present(position_)) then
!
         read(this%unit_, pos=position_, iostat=io_status, iomsg=io_message) array
!
      else
!
         read(this%unit_, iostat=io_status, iomsg=io_message) array
!
      end if
!
      call this%check_io_status(io_status, io_message, task='read', status_=status_)
!
   end subroutine read_1_real_dp_abstract_stream
!
!
   subroutine read_1_complex_dp_abstract_stream(this, array, n, position_, status_)
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
      class(abstract_stream), intent(in) :: this
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
      call this%error_if_closed(task='read')
!
      if(present(position_)) then
!
         read(this%unit_, pos=position_, iostat=io_status, iomsg=io_message) array
!
      else
!
         read(this%unit_, iostat=io_status, iomsg=io_message) array
!
      end if
!
      call this%check_io_status(io_status, io_message, task='read', status_=status_)
!
   end subroutine read_1_complex_dp_abstract_stream
!
!
   subroutine read_1_int_32_abstract_stream(this, array, n, position_, status_)
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
      class(abstract_stream), intent(in) :: this
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
      call this%error_if_closed(task='read')
!
      if(present(position_)) then
!
         read(this%unit_, pos=position_, iostat=io_status, iomsg=io_message) array
!
      else
!
         read(this%unit_, iostat=io_status, iomsg=io_message) array
!
      end if
!
      call this%check_io_status(io_status, io_message, task='read', status_=status_)
!
   end subroutine read_1_int_32_abstract_stream
!
!
   subroutine read_1_int_64_abstract_stream(this, array, n, position_, status_)
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
      class(abstract_stream), intent(in) :: this
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
      call this%error_if_closed(task='read')
!
      if(present(position_)) then
!
         read(this%unit_, pos=position_, iostat=io_status, iomsg=io_message) array
!
      else
!
         read(this%unit_, iostat=io_status, iomsg=io_message) array
!
      end if
!
      call this%check_io_status(io_status, io_message, task='read', status_=status_)
!
   end subroutine read_1_int_64_abstract_stream
!
!
   subroutine read_1_log_abstract_stream(this, array, n, position_, status_)
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
      class(abstract_stream), intent(in) :: this
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
      call this%error_if_closed(task='read')
!
      if(present(position_)) then
!
         read(this%unit_, pos=position_, iostat=io_status, iomsg=io_message) array
!
      else
!
         read(this%unit_, iostat=io_status, iomsg=io_message) array
!
      end if
!
      call this%check_io_status(io_status, io_message, task='read', status_=status_)
!
   end subroutine read_1_log_abstract_stream
!
!
   subroutine write_0_real_dp_abstract_stream(this, scalar, position_)
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
      class(abstract_stream), intent(in) :: this
!
      integer, intent(in), optional :: position_
!
      real(dp), intent(in) :: scalar
!
      integer::             io_status = 1
      character(len=100) :: io_message
!
      call this%error_if_closed(task='write')
!
      if(present(position_)) then
!
         write(this%unit_, pos=position_, iostat=io_status, iomsg=io_message) scalar
!
      else
!
         write(this%unit_, iostat=io_status, iomsg=io_message) scalar
!
      end if
!
      call this%check_io_status(io_status, io_message, task='write')
!
   end subroutine write_0_real_dp_abstract_stream
!
!
   subroutine write_0_int_32_abstract_stream(this, scalar, position_)
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
      class(abstract_stream), intent(in) :: this
!
      integer, intent(in), optional :: position_
!
      integer(i32), intent(in) :: scalar
!
      integer::             io_status = 1
      character(len=100) :: io_message
!
      call this%error_if_closed(task='write')
!
      if(present(position_)) then
!
         write(this%unit_, pos=position_, iostat=io_status, iomsg=io_message) scalar
!
      else
!
         write(this%unit_, iostat=io_status, iomsg=io_message) scalar
!
      end if
!
      call this%check_io_status(io_status, io_message, task='write')
!
   end subroutine write_0_int_32_abstract_stream
!
!
   subroutine write_0_int_64_abstract_stream(this, scalar, position_)
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
      class(abstract_stream), intent(in) :: this
!
      integer, intent(in), optional :: position_
!
      integer(i64), intent(in) :: scalar
!
      integer::             io_status = 1
      character(len=100) :: io_message
!
      call this%error_if_closed(task='write')
!
      if(present(position_)) then
!
         write(this%unit_, pos=position_, iostat=io_status, iomsg=io_message) scalar
!
      else
!
         write(this%unit_, iostat=io_status, iomsg=io_message) scalar
!
      end if
!
      call this%check_io_status(io_status, io_message, task='write')
!
   end subroutine write_0_int_64_abstract_stream
!
!
   subroutine write_1_real_dp_abstract_stream(this, array, n, position_)
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
      class(abstract_stream), intent(in) :: this
!
      integer, intent(in) :: n
      integer, intent(in), optional :: position_
!
      real(dp), dimension(n), intent(in) :: array
!
      integer::             io_status = 1
      character(len=100) :: io_message
!
      call this%error_if_closed(task='write')
!
      if(present(position_)) then
!
         write(this%unit_, pos=position_, iostat=io_status, iomsg=io_message) array
!
      else
!
         write(this%unit_, iostat=io_status, iomsg=io_message) array
!
      end if
!
      call this%check_io_status(io_status, io_message, task='write')
!
   end subroutine write_1_real_dp_abstract_stream
!
!
   subroutine write_1_real_sp_abstract_stream(this, array, n, position_)
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
      class(abstract_stream), intent(in) :: this
!
      integer, intent(in) :: n
      integer, intent(in), optional :: position_
!
      real(sp), dimension(n), intent(in) :: array
!
      integer::             io_status = 1
      character(len=100) :: io_message
!
      call this%error_if_closed(task='write')
!
      if(present(position_)) then
!
         write(this%unit_, pos=position_, iostat=io_status, iomsg=io_message) array
!
      else
!
         write(this%unit_, iostat=io_status, iomsg=io_message) array
!
      end if
!
      call this%check_io_status(io_status, io_message, task='write')
!
   end subroutine write_1_real_sp_abstract_stream
!
!
   subroutine write_1_complex_dp_abstract_stream(this, array, n, position_)
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
      class(abstract_stream), intent(in) :: this
!
      integer, intent(in) :: n
      integer, intent(in), optional :: position_
!
      complex(dp), dimension(n), intent(in) :: array
!
      integer::             io_status = 1
      character(len=100) :: io_message
!
      call this%error_if_closed(task='write')
!
      if(present(position_)) then
!
         write(this%unit_, pos=position_, iostat=io_status, iomsg=io_message) array
!
      else
!
         write(this%unit_, iostat=io_status, iomsg=io_message) array
!
      end if
!
      call this%check_io_status(io_status, io_message, task='write')
!
   end subroutine write_1_complex_dp_abstract_stream
!
!
   subroutine write_1_int_32_abstract_stream(this, array, n, position_)
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
      class(abstract_stream), intent(in) :: this
!
      integer, intent(in) :: n
      integer, intent(in), optional :: position_
!
      integer(i32), dimension(n), intent(in) :: array
!
      integer::             io_status = 1
      character(len=100) :: io_message
!
      call this%error_if_closed(task='write')
!
      if(present(position_)) then
!
         write(this%unit_, pos=position_, iostat=io_status, iomsg=io_message) array
!
      else
!
         write(this%unit_, iostat=io_status, iomsg=io_message) array
!
      end if
!
      call this%check_io_status(io_status, io_message, task='write')
!
   end subroutine write_1_int_32_abstract_stream
!
!
   subroutine write_1_int_64_abstract_stream(this, array, n, position_)
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
      class(abstract_stream), intent(in) :: this
!
      integer, intent(in) :: n
      integer, intent(in), optional :: position_
!
      integer(i64), dimension(n), intent(in) :: array
!
      integer::             io_status = 1
      character(len=100) :: io_message
!
      call this%error_if_closed(task='write')
!
      if(present(position_)) then
!
         write(this%unit_, pos=position_, iostat=io_status, iomsg=io_message) array
!
      else
!
         write(this%unit_, iostat=io_status, iomsg=io_message) array
!
      end if
!
      call this%check_io_status(io_status, io_message, task='write')
!
   end subroutine write_1_int_64_abstract_stream
!
!
   subroutine write_1_log_abstract_stream(this, array, n, position_)
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
      class(abstract_stream), intent(in) :: this
!
      integer, intent(in) :: n
      integer, intent(in), optional :: position_
!
      logical, dimension(n), intent(in) :: array
!
      integer::             io_status = 1
      character(len=100) :: io_message
!
      call this%error_if_closed(task='write')
!
      if(present(position_)) then
!
         write(this%unit_, pos=position_, iostat=io_status, iomsg=io_message) array
!
      else
!
         write(this%unit_, iostat=io_status, iomsg=io_message) array
!
      end if
!
      call this%check_io_status(io_status, io_message, task='write')
!
   end subroutine write_1_log_abstract_stream
!
!
   subroutine io_error_abstract_stream(this, io_status, io_message, task)
!!
!!    io error
!!    Written by Alexander C. Paul, Sep 2022
!!
!!    task: string specifying which task is checked: open, close, read, write
!!
      implicit none
!
      class(abstract_stream), intent(in) :: this
!
      integer, intent(in) :: io_status
!
      character(len=*), intent(in) :: io_message
      character(len=*), intent(in) :: task
!
      call output%error_msg('Failed to '// task // ' stream file (a0), status is &
                            &(i0) and error message is: ' // trim(io_message), &
                            chars = [this%get_name()], ints = [io_status])
!
   end subroutine io_error_abstract_stream
!
!
end module abstract_stream_class
