!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2020 the authors of eT
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
module stream_file_class
!
!!
!!    Stream file class module
!!    Writen by Alexander C. Paul, March 2020
!!
!!    Contains constructor and destructor for stream files and
!!    wrapper routines to handle arrays of higher ranks than 1.
!!
!
   use kinds
   use abstract_stream_class, only : abstract_stream
   use global_out, only : output
!
   type, extends(abstract_stream) :: stream_file
!
   contains
!
!     Read routines
!
!     Real double precision read
!
      procedure :: read_2_real_dp_stream_file
      procedure :: read_3_real_dp_stream_file
      procedure :: read_4_real_dp_stream_file
!
!     Complex double precision read
!
      procedure :: read_2_complex_dp_stream_file
      procedure :: read_3_complex_dp_stream_file
      procedure :: read_4_complex_dp_stream_file
!
!     32-bit integer read
!
      procedure :: read_2_int_32_stream_file
      procedure :: read_3_int_32_stream_file
      procedure :: read_4_int_32_stream_file
!
!     64-bit integer read
!
      procedure :: read_2_int_64_stream_file
      procedure :: read_3_int_64_stream_file
      procedure :: read_4_int_64_stream_file
!
!     Read generic
!
      generic :: read_ => read_real_dp_scalar,           &
                          read_real_dp,                  &
                          read_2_real_dp_stream_file,    &
                          read_3_real_dp_stream_file,    &
                          read_4_real_dp_stream_file,    &
                          read_complex_dp_scalar,        &
                          read_complex_dp,               &
                          read_2_complex_dp_stream_file, &
                          read_3_complex_dp_stream_file, &
                          read_4_complex_dp_stream_file, &
                          read_int_32_scalar,            &
                          read_int_32,                   &
                          read_2_int_32_stream_file,     &
                          read_3_int_32_stream_file,     &
                          read_4_int_32_stream_file,     &
                          read_int_64_scalar,            &
                          read_int_64,                   &
                          read_2_int_64_stream_file,     &
                          read_3_int_64_stream_file,     &
                          read_4_int_64_stream_file,     &
                          read_log_scalar,               &
                          read_log
!
!     Write routines
!
!
!     Real double precision write
!
      procedure :: write_2_real_dp_stream_file
      procedure :: write_3_real_dp_stream_file
      procedure :: write_4_real_dp_stream_file
!
!     Complex double precision write
!
      procedure :: write_2_complex_dp_stream_file
      procedure :: write_3_complex_dp_stream_file
      procedure :: write_4_complex_dp_stream_file
!
!     32-bit integer read
!
      procedure :: write_2_int_32_stream_file
      procedure :: write_3_int_32_stream_file
      procedure :: write_4_int_32_stream_file
!
!     64-bit integer read
!
      procedure :: write_2_int_64_stream_file
      procedure :: write_3_int_64_stream_file
      procedure :: write_4_int_64_stream_file
!
!     Write generic
!
      generic :: write_ => write_real_dp_scalar,            &
                           write_real_dp,                   &
                           write_2_real_dp_stream_file,     &
                           write_3_real_dp_stream_file,     &
                           write_4_real_dp_stream_file,     &
                           write_complex_dp_scalar,         &
                           write_complex_dp,                &
                           write_2_complex_dp_stream_file,  &
                           write_3_complex_dp_stream_file,  &
                           write_4_complex_dp_stream_file,  &
                           write_int_32_scalar,             &
                           write_int_32,                    &
                           write_2_int_32_stream_file,      &
                           write_3_int_32_stream_file,      &
                           write_4_int_32_stream_file,      &
                           write_int_64_scalar,             &
                           write_int_64,                    &
                           write_2_int_64_stream_file,      &
                           write_3_int_64_stream_file,      &
                           write_4_int_64_stream_file,      &
                           write_log_scalar,                &
                           write_log
!
      final :: destructor
!
   end type stream_file
!
   interface stream_file
!
      procedure new_stream_file
!
   end interface stream_file
!
contains
!
!
   function new_stream_file(name_, status_) result(the_file)
!!
!!    Stream file constructer
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      type(stream_file) :: the_file
!
      character(len=*), intent(in)           :: name_
      character(len=*), intent(in), optional :: status_
!
      call the_file%set_name(name_)
!
      if (present(status_)) then
         call the_file%set_status(status_)
      else
         call the_file%set_status('unknown')
      endif
!
   end function new_stream_file
!
!
   subroutine destructor(the_file)
!!
!!    Destructor
!!    Written by Rolf H. Myhre, Feb. 2020
!!
      implicit none
!
      type(stream_file) :: the_file
!
      if (the_file%get_open()) then
         call output%error_msg('Destructor for file (a0) called &
                               &while the file is still open',  &
                               chars=[the_file%get_name()])
      endif
!
   end subroutine destructor
!
!!
!! Wrapper routines for read_real_dp_abstract_stream
!! that accept rank 2, 3, and 4 arrays
!!
!! see documentation below
!!
!
   subroutine read_2_real_dp_stream_file(the_file, array, n, position_)
!!
!!    Read 2 real dp
!!    Written by Rolf H. Myhre and Alexander C. Paul, Mar 2020
!!
!!    array: real double precision array of length n and rank 2
!!
!!    n: total length of array (Number of elements)
!!
!!    position_:  optional integer, position to read from in file
!!                positions counted in bytes, starting at 1
!!                default: current file pointer position
!!
!
      implicit none
!
      class(stream_file), intent(in) :: the_file
      integer, intent(in) :: n
      real(dp), dimension(:,:), intent(inout) :: array
      integer, intent(in), optional :: position_
!
      call the_file%read_real_dp(array, n, position_)
!
   end subroutine read_2_real_dp_stream_file
!
!
   subroutine read_3_real_dp_stream_file(the_file, array, n, position_)
!
      implicit none
!
      class(stream_file), intent(in) :: the_file
      integer, intent(in) :: n
      real(dp), dimension(:,:,:), intent(inout) :: array
      integer, intent(in), optional :: position_
!
      call the_file%read_real_dp(array, n, position_)
!
   end subroutine read_3_real_dp_stream_file
!
!
   subroutine read_4_real_dp_stream_file(the_file, array, n, position_)
!
      implicit none
!
      class(stream_file), intent(in) :: the_file
      integer, intent(in) :: n
      real(dp), dimension(:,:,:,:), intent(inout) :: array
      integer, intent(in), optional :: position_
!
      call the_file%read_real_dp(array, n, position_)
!
   end subroutine read_4_real_dp_stream_file
!
!  Complex arrays
!
   subroutine read_2_complex_dp_stream_file(the_file, array, n, position_)
!
      implicit none
!
      class(stream_file), intent(in) :: the_file
      integer, intent(in) :: n
      complex(dp), dimension(:,:), intent(inout) :: array
      integer, intent(in), optional :: position_
!
      call the_file%read_complex_dp(array, n, position_)
!
   end subroutine read_2_complex_dp_stream_file
!
!
   subroutine read_3_complex_dp_stream_file(the_file, array, n, position_)
!
      implicit none
!
      class(stream_file), intent(in) :: the_file
      integer, intent(in) :: n
      complex(dp), dimension(:,:,:), intent(inout) :: array
      integer, intent(in), optional :: position_
!
      call the_file%read_complex_dp(array, n, position_)
!
   end subroutine read_3_complex_dp_stream_file
!
!
   subroutine read_4_complex_dp_stream_file(the_file, array, n, position_)
!
      implicit none
!
      class(stream_file), intent(in) :: the_file
      integer, intent(in) :: n
      complex(dp), dimension(:,:,:,:), intent(inout) :: array
      integer, intent(in), optional :: position_
!
      call the_file%read_complex_dp(array, n, position_)
!
   end subroutine read_4_complex_dp_stream_file
!
!  64-bit integer arrays
!
   subroutine read_2_int_64_stream_file(the_file, array, n, position_)
!
      implicit none
!
      class(stream_file), intent(in) :: the_file
      integer, intent(in) :: n
      integer(i64), dimension(:,:), intent(inout) :: array
      integer, intent(in), optional :: position_
!
      call the_file%read_int_64(array, n, position_)
!
   end subroutine read_2_int_64_stream_file
!
!
   subroutine read_3_int_64_stream_file(the_file, array, n, position_)
!
      implicit none
!
      class(stream_file), intent(in) :: the_file
      integer, intent(in) :: n
      integer(i64), dimension(:,:,:), intent(inout) :: array
      integer, intent(in), optional :: position_
!
      call the_file%read_int_64(array, n, position_)
!
   end subroutine read_3_int_64_stream_file
!
!
   subroutine read_4_int_64_stream_file(the_file, array, n, position_)
!
      implicit none
!
      class(stream_file), intent(in) :: the_file
      integer, intent(in) :: n
      integer(i64), dimension(:,:,:,:), intent(inout) :: array
      integer, intent(in), optional :: position_
!
      call the_file%read_int_64(array, n, position_)
!
   end subroutine read_4_int_64_stream_file 
!
!  32-bit integer arrays
!
   subroutine read_2_int_32_stream_file(the_file, array, n, position_)
!
      implicit none
!
      class(stream_file), intent(in) :: the_file
      integer, intent(in) :: n
      integer(i32), dimension(:,:), intent(inout) :: array
      integer, intent(in), optional :: position_
!
      call the_file%read_int_32(array, n, position_)
!
   end subroutine read_2_int_32_stream_file
!
!
   subroutine read_3_int_32_stream_file(the_file, array, n, position_)
!
      implicit none
!
      class(stream_file), intent(in) :: the_file
      integer, intent(in) :: n
      integer(i32), dimension(:,:,:), intent(inout) :: array
      integer, intent(in), optional :: position_
!
      call the_file%read_int_32(array, n, position_)
!
   end subroutine read_3_int_32_stream_file
!
!
   subroutine read_4_int_32_stream_file(the_file, array, n, position_)
!
      implicit none
!
      class(stream_file), intent(in) :: the_file
      integer, intent(in) :: n
      integer(i32), dimension(:,:,:,:), intent(inout) :: array
      integer, intent(in), optional :: position_
!
      call the_file%read_int_32(array, n, position_)
!
   end subroutine read_4_int_32_stream_file
!
!!
!! Wrapper routines for write_real_dp_abstract_stream
!! that accept rank 2, 3, and 4 arrays
!!
!! see documentation below
!!
!
   subroutine write_2_real_dp_stream_file(the_file, array, n, position_)
!!
!!    Write 2 real dp
!!    Written by Rolf H. Myhre and Alexander C. Paul, Mar 2020
!!
!!    array: real double precision array of length n and rank 2
!!
!!    n: total length of array (Number of elements)
!!
!!    position_:  optional integer, position to read from in file
!!                positions counted in bytes, starting at 1
!!                default: current file pointer position
!!
!
      implicit none
!
      class(stream_file), intent(in) :: the_file
      integer, intent(in) :: n
      real(dp), dimension(:,:), intent(in) :: array
      integer, intent(in), optional :: position_
!
      call the_file%write_real_dp(array, n, position_)
!
   end subroutine write_2_real_dp_stream_file
!
!
   subroutine write_3_real_dp_stream_file(the_file, array, n, position_)
!
      implicit none
!
      class(stream_file), intent(in) :: the_file
      integer, intent(in) :: n
      real(dp), dimension(:,:,:), intent(in) :: array
      integer, intent(in), optional :: position_
!
      call the_file%write_real_dp(array, n, position_)
!
   end subroutine write_3_real_dp_stream_file
!
!
   subroutine write_4_real_dp_stream_file(the_file, array, n, position_)
!
      implicit none
!
      class(stream_file), intent(in) :: the_file
      integer, intent(in) :: n
      real(dp), dimension(:,:,:,:), intent(in) :: array
      integer, intent(in), optional :: position_
!
      call the_file%write_real_dp(array, n, position_)
!
   end subroutine write_4_real_dp_stream_file
!
!  Complex arrays
!
   subroutine write_2_complex_dp_stream_file(the_file, array, n, position_)
!
      implicit none
!
      class(stream_file), intent(in) :: the_file
      integer, intent(in) :: n
      complex(dp), dimension(:,:), intent(in) :: array
      integer, intent(in), optional :: position_
!
      call the_file%write_complex_dp(array, n, position_)
!
   end subroutine write_2_complex_dp_stream_file
!
!
   subroutine write_3_complex_dp_stream_file(the_file, array, n, position_)
!
      implicit none
!
      class(stream_file), intent(in) :: the_file
      integer, intent(in) :: n
      complex(dp), dimension(:,:,:), intent(in) :: array
      integer, intent(in), optional :: position_
!
      call the_file%write_complex_dp(array, n, position_)
!
   end subroutine write_3_complex_dp_stream_file
!
!
   subroutine write_4_complex_dp_stream_file(the_file, array, n, position_)
!
      implicit none
!
      class(stream_file), intent(in) :: the_file
      integer, intent(in) :: n
      complex(dp), dimension(:,:,:,:), intent(in) :: array
      integer, intent(in), optional :: position_
!
      call the_file%write_complex_dp(array, n, position_)
!
   end subroutine write_4_complex_dp_stream_file
!
!  64-bit integer arrays
!
   subroutine write_2_int_64_stream_file(the_file, array, n, position_)
!
      implicit none
!
      class(stream_file), intent(in) :: the_file
      integer, intent(in) :: n
      integer(i64), dimension(:,:), intent(in) :: array
      integer, intent(in), optional :: position_
!
      call the_file%write_int_64(array, n, position_)
!
   end subroutine write_2_int_64_stream_file
!
!
   subroutine write_3_int_64_stream_file(the_file, array, n, position_)
!
      implicit none
!
      class(stream_file), intent(in) :: the_file
      integer, intent(in) :: n
      integer(i64), dimension(:,:,:), intent(in) :: array
      integer, intent(in), optional :: position_
!
      call the_file%write_int_64(array, n, position_)
!
   end subroutine write_3_int_64_stream_file
!
!
   subroutine write_4_int_64_stream_file(the_file, array, n, position_)
!
      implicit none
!
      class(stream_file), intent(in) :: the_file
      integer, intent(in) :: n
      integer(i64), dimension(:,:,:,:), intent(in) :: array
      integer, intent(in), optional :: position_
!
      call the_file%write_int_64(array, n, position_)
!
   end subroutine write_4_int_64_stream_file 
!
!  32-bit integer arrays
!
   subroutine write_2_int_32_stream_file(the_file, array, n, position_)
!
      implicit none
!
      class(stream_file), intent(in) :: the_file
      integer, intent(in) :: n
      integer(i32), dimension(:,:), intent(in) :: array
      integer, intent(in), optional :: position_
!
      call the_file%write_int_32(array, n, position_)
!
   end subroutine write_2_int_32_stream_file
!
!
   subroutine write_3_int_32_stream_file(the_file, array, n, position_)
!
      implicit none
!
      class(stream_file), intent(in) :: the_file
      integer, intent(in) :: n
      integer(i32), dimension(:,:,:), intent(in) :: array
      integer, intent(in), optional :: position_
!
      call the_file%write_int_32(array, n, position_)
!
   end subroutine write_3_int_32_stream_file
!
!
   subroutine write_4_int_32_stream_file(the_file, array, n, position_)
!
      implicit none
!
      class(stream_file), intent(in) :: the_file
      integer, intent(in) :: n
      integer(i32), dimension(:,:,:,:), intent(in) :: array
      integer, intent(in), optional :: position_
!
      call the_file%write_int_32(array, n, position_)
!
   end subroutine write_4_int_32_stream_file
!
!
end module stream_file_class
