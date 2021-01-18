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
      procedure :: read_0_real_dp_stream_file
      procedure :: read_1_real_dp_stream_file
      procedure :: read_2_real_dp_stream_file
      procedure :: read_3_real_dp_stream_file
      procedure :: read_4_real_dp_stream_file
!
!     Real single precision read
!
      procedure :: read_0_real_sp_stream_file
      procedure :: read_1_real_sp_stream_file
!
!     Complex double precision read
!
      procedure :: read_0_complex_dp_stream_file
      procedure :: read_1_complex_dp_stream_file
      procedure :: read_2_complex_dp_stream_file
      procedure :: read_3_complex_dp_stream_file
      procedure :: read_4_complex_dp_stream_file
!
!     32-bit integer read
!
      procedure :: read_0_int_32_stream_file
      procedure :: read_1_int_32_stream_file
      procedure :: read_2_int_32_stream_file
      procedure :: read_3_int_32_stream_file
      procedure :: read_4_int_32_stream_file
!
!     64-bit integer read
!
      procedure :: read_0_int_64_stream_file
      procedure :: read_1_int_64_stream_file
      procedure :: read_2_int_64_stream_file
      procedure :: read_3_int_64_stream_file
      procedure :: read_4_int_64_stream_file
!
!     Logical read
!
      procedure :: read_0_log_stream_file
      procedure :: read_1_log_stream_file
!
!     Read generic
!
      generic :: read_ => read_0_real_dp_stream_file,    &
                          read_1_real_dp_stream_file,    &
                          read_2_real_dp_stream_file,    &
                          read_3_real_dp_stream_file,    &
                          read_4_real_dp_stream_file,    &
                          read_0_real_sp_stream_file,    &
                          read_1_real_sp_stream_file,    &
                          read_0_complex_dp_stream_file, &
                          read_1_complex_dp_stream_file, &
                          read_2_complex_dp_stream_file, &
                          read_3_complex_dp_stream_file, &
                          read_4_complex_dp_stream_file, &
                          read_0_int_32_stream_file,     &
                          read_1_int_32_stream_file,     &
                          read_2_int_32_stream_file,     &
                          read_3_int_32_stream_file,     &
                          read_4_int_32_stream_file,     &
                          read_0_int_64_stream_file,     &
                          read_1_int_64_stream_file,     &
                          read_2_int_64_stream_file,     &
                          read_3_int_64_stream_file,     &
                          read_4_int_64_stream_file,     &
                          read_0_log_stream_file,        &
                          read_1_log_stream_file
!
!     Write routines
!
!     Real double precision write
!
      procedure :: write_0_real_dp_stream_file
      procedure :: write_1_real_dp_stream_file
      procedure :: write_2_real_dp_stream_file
      procedure :: write_3_real_dp_stream_file
      procedure :: write_4_real_dp_stream_file
!
!     Real single precision write
!
      procedure :: write_0_real_sp_stream_file
      procedure :: write_1_real_sp_stream_file
!
!     Complex double precision write
!
      procedure :: write_0_complex_dp_stream_file
      procedure :: write_1_complex_dp_stream_file
      procedure :: write_2_complex_dp_stream_file
      procedure :: write_3_complex_dp_stream_file
      procedure :: write_4_complex_dp_stream_file
!
!     32-bit integer write
!
      procedure :: write_0_int_32_stream_file
      procedure :: write_1_int_32_stream_file
      procedure :: write_2_int_32_stream_file
      procedure :: write_3_int_32_stream_file
      procedure :: write_4_int_32_stream_file
!
!     64-bit integer write
!
      procedure :: write_0_int_64_stream_file
      procedure :: write_1_int_64_stream_file
      procedure :: write_2_int_64_stream_file
      procedure :: write_3_int_64_stream_file
      procedure :: write_4_int_64_stream_file
!
!     Logical write
!
      procedure :: write_0_log_stream_file
      procedure :: write_1_log_stream_file
!
!     write generic
!
      generic :: write_ => write_0_real_dp_stream_file,    &
                           write_1_real_dp_stream_file,    &
                           write_2_real_dp_stream_file,    &
                           write_3_real_dp_stream_file,    &
                           write_4_real_dp_stream_file,    &
                           write_0_real_sp_stream_file,    &
                           write_1_real_sp_stream_file,    &
                           write_0_complex_dp_stream_file, &
                           write_1_complex_dp_stream_file, &
                           write_2_complex_dp_stream_file, &
                           write_3_complex_dp_stream_file, &
                           write_4_complex_dp_stream_file, &
                           write_0_int_32_stream_file,     &
                           write_1_int_32_stream_file,     &
                           write_2_int_32_stream_file,     &
                           write_3_int_32_stream_file,     &
                           write_4_int_32_stream_file,     &
                           write_0_int_64_stream_file,     &
                           write_1_int_64_stream_file,     &
                           write_2_int_64_stream_file,     &
                           write_3_int_64_stream_file,     &
                           write_4_int_64_stream_file,     &
                           write_0_log_stream_file,        &
                           write_1_log_stream_file
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
!! Wrapper routines for read_x_real_dp_abstract_stream
!! that accept scalars and rank 1, 2, 3, and 4 arrays
!!
!! The wrapper for scalars and rank 1 are needed because ifort does not permit 
!! a generic procedure pointer to a routine of the parent class
!!
!! see documentation below
!!
!
   subroutine read_0_real_dp_stream_file(the_file, scalar, position_, status_)
!!
!!    Read rank 0 real dp
!!    Written by Rolf H. Myhre and Alexander C. Paul, Mar 2020
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
      class(stream_file), intent(in) :: the_file
      real(dp), intent(out)          :: scalar
      integer, intent(in), optional  :: position_
      integer, intent(out), optional :: status_
!
      call the_file%read_0_real_dp_abstract_stream(scalar, position_, status_)
!
   end subroutine read_0_real_dp_stream_file
!
   subroutine read_1_real_dp_stream_file(the_file, array, n, position_, status_)
!!
!!    Read rank 1 real dp
!!    Written by Rolf H. Myhre and Alexander C. Paul, Mar 2020
!!
!!    array: real double precision array of length n
!!
!!    n: total length of array (Number of elements)
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
      class(stream_file), intent(in)      :: the_file
      integer, intent(in)                 :: n
      real(dp), dimension(n), intent(out) :: array
      integer, intent(in), optional       :: position_
      integer, intent(out), optional      :: status_
!
      call the_file%read_1_real_dp_abstract_stream(array, n, position_, status_)
!
   end subroutine read_1_real_dp_stream_file
!
   subroutine read_2_real_dp_stream_file(the_file, array, n, position_, status_)
!!
!!    Read rank 2 real dp
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
!!    status_: optional integer, returns iostat e.g. to check for end of file
!!             on the outside
!!
!
      implicit none
!
      class(stream_file), intent(in)        :: the_file
      real(dp), dimension(:,:), intent(out) :: array
      integer, intent(in)                   :: n
      integer, intent(in), optional         :: position_
      integer, intent(out), optional        :: status_
!
      call the_file%read_1_real_dp_stream_file(array, n, position_, status_)
!
   end subroutine read_2_real_dp_stream_file
!
   subroutine read_3_real_dp_stream_file(the_file, array, n, position_, status_)
!
      implicit none
!
      class(stream_file), intent(in)          :: the_file
      real(dp), dimension(:,:,:), intent(out) :: array
      integer, intent(in)                     :: n
      integer, intent(in), optional           :: position_
      integer, intent(out), optional          :: status_
!
      call the_file%read_1_real_dp_stream_file(array, n, position_, status_)
!
   end subroutine read_3_real_dp_stream_file
!
   subroutine read_4_real_dp_stream_file(the_file, array, n, position_, status_)
!
      implicit none
!
      class(stream_file), intent(in)            :: the_file
      real(dp), dimension(:,:,:,:), intent(out) :: array
      integer, intent(in)                       :: n
      integer, intent(in), optional             :: position_
      integer, intent(out), optional            :: status_
!
      call the_file%read_1_real_dp_stream_file(array, n, position_, status_)
!
   end subroutine read_4_real_dp_stream_file
!
!  Real single precision
!
   subroutine read_0_real_sp_stream_file(the_file, scalar, position_, status_)
!
      implicit none
!
      class(stream_file), intent(in) :: the_file
      real(sp), intent(out)          :: scalar
      integer, intent(in), optional  :: position_
      integer, intent(out), optional :: status_
!
      call the_file%read_0_real_sp_abstract_stream(scalar, position_, status_)
!
   end subroutine read_0_real_sp_stream_file
!
   subroutine read_1_real_sp_stream_file(the_file, array, n, position_, status_)
!
      implicit none
!
      class(stream_file), intent(in)      :: the_file
      integer, intent(in)                 :: n
      real(sp), dimension(n), intent(out) :: array
      integer, intent(in), optional       :: position_
      integer, intent(out), optional      :: status_
!
      call the_file%read_1_real_sp_abstract_stream(array, n, position_, status_)
!
   end subroutine read_1_real_sp_stream_file
!
!  Complex
!
   subroutine read_0_complex_dp_stream_file(the_file, scalar, position_, status_)
!
      implicit none
!
      class(stream_file), intent(in) :: the_file
      complex(dp), intent(out)       :: scalar
      integer, intent(in), optional  :: position_
      integer, intent(out), optional :: status_
!
      call the_file%read_0_complex_dp_abstract_stream(scalar, position_, status_)
!
   end subroutine read_0_complex_dp_stream_file
!
   subroutine read_1_complex_dp_stream_file(the_file, array, n, position_, status_)
!
      implicit none
!
      class(stream_file), intent(in)         :: the_file
      integer, intent(in)                    :: n
      complex(dp), dimension(n), intent(out) :: array
      integer, intent(in), optional          :: position_
      integer, intent(out), optional         :: status_
!
      call the_file%read_1_complex_dp_abstract_stream(array, n, position_, status_)
!
   end subroutine read_1_complex_dp_stream_file
!
   subroutine read_2_complex_dp_stream_file(the_file, array, n, position_, status_)
!
      implicit none
!
      class(stream_file), intent(in)           :: the_file
      complex(dp), dimension(:,:), intent(out) :: array
      integer, intent(in)                      :: n
      integer, intent(in), optional            :: position_
      integer, intent(out), optional           :: status_
!
      call the_file%read_1_complex_dp_stream_file(array, n, position_, status_)
!
   end subroutine read_2_complex_dp_stream_file
!
   subroutine read_3_complex_dp_stream_file(the_file, array, n, position_, status_)
!
      implicit none
!
      class(stream_file), intent(in)             :: the_file
      complex(dp), dimension(:,:,:), intent(out) :: array
      integer, intent(in)                        :: n
      integer, intent(in), optional              :: position_
      integer, intent(out), optional             :: status_
!
      call the_file%read_1_complex_dp_stream_file(array, n, position_, status_)
!
   end subroutine read_3_complex_dp_stream_file
!
   subroutine read_4_complex_dp_stream_file(the_file, array, n, position_, status_)
!
      implicit none
!
      class(stream_file), intent(in)               :: the_file
      complex(dp), dimension(:,:,:,:), intent(out) :: array
      integer, intent(in)                          :: n
      integer, intent(in), optional                :: position_
      integer, intent(out), optional               :: status_
!
      call the_file%read_1_complex_dp_stream_file(array, n, position_, status_)
!
   end subroutine read_4_complex_dp_stream_file
!
!  64-bit integers
!
   subroutine read_0_int_64_stream_file(the_file, scalar, position_, status_)
!
      implicit none
!
      class(stream_file), intent(in) :: the_file
      integer(i64), intent(out)      :: scalar
      integer, intent(in), optional  :: position_
      integer, intent(out), optional :: status_
!
      call the_file%read_0_int_64_abstract_stream(scalar, position_, status_)
!
   end subroutine read_0_int_64_stream_file
!
   subroutine read_1_int_64_stream_file(the_file, array, n, position_, status_)
!
      implicit none
!
      class(stream_file), intent(in)          :: the_file
      integer, intent(in)                     :: n
      integer(i64), dimension(n), intent(out) :: array
      integer, intent(in), optional           :: position_
      integer, intent(out), optional          :: status_
!
      call the_file%read_1_int_64_abstract_stream(array, n, position_, status_)
!
   end subroutine read_1_int_64_stream_file
!
   subroutine read_2_int_64_stream_file(the_file, array, n, position_, status_)
!
      implicit none
!
      class(stream_file), intent(in)            :: the_file
      integer(i64), dimension(:,:), intent(out) :: array
      integer, intent(in)                       :: n
      integer, intent(in), optional             :: position_
      integer, intent(out), optional            :: status_
!
      call the_file%read_1_int_64_stream_file(array, n, position_, status_)
!
   end subroutine read_2_int_64_stream_file
!
   subroutine read_3_int_64_stream_file(the_file, array, n, position_, status_)
!
      implicit none
!
      class(stream_file), intent(in)              :: the_file
      integer(i64), dimension(:,:,:), intent(out) :: array
      integer, intent(in)                         :: n
      integer, intent(in), optional               :: position_
      integer, intent(out), optional              :: status_
!
      call the_file%read_1_int_64_stream_file(array, n, position_, status_)
!
   end subroutine read_3_int_64_stream_file
!
   subroutine read_4_int_64_stream_file(the_file, array, n, position_, status_)
!
      implicit none
!
      class(stream_file), intent(in)                :: the_file
      integer(i64), dimension(:,:,:,:), intent(out) :: array
      integer, intent(in)                           :: n
      integer, intent(in), optional                 :: position_
      integer, intent(out), optional                :: status_
!
      call the_file%read_1_int_64_stream_file(array, n, position_, status_)
!
   end subroutine read_4_int_64_stream_file 
!
!  32-bit integers
!
   subroutine read_0_int_32_stream_file(the_file, scalar, position_, status_)
!
      implicit none
!
      class(stream_file), intent(in) :: the_file
      integer(i32), intent(out)      :: scalar
      integer, intent(in), optional  :: position_
      integer, intent(out), optional :: status_
!
      call the_file%read_0_int_32_abstract_stream(scalar, position_, status_)
!
   end subroutine read_0_int_32_stream_file
!
   subroutine read_1_int_32_stream_file(the_file, array, n, position_, status_)
!
      implicit none
!
      class(stream_file), intent(in)          :: the_file
      integer, intent(in)                     :: n
      integer(i32), dimension(n), intent(out) :: array
      integer, intent(in), optional           :: position_
      integer, intent(out), optional          :: status_
!
      call the_file%read_1_int_32_abstract_stream(array, n, position_, status_)
!
   end subroutine read_1_int_32_stream_file
!
   subroutine read_2_int_32_stream_file(the_file, array, n, position_, status_)
!
      implicit none
!
      class(stream_file), intent(in)            :: the_file
      integer(i32), dimension(:,:), intent(out) :: array
      integer, intent(in)                       :: n
      integer, intent(in), optional             :: position_
      integer, intent(out), optional            :: status_
!
      call the_file%read_1_int_32_stream_file(array, n, position_, status_)
!
   end subroutine read_2_int_32_stream_file
!
   subroutine read_3_int_32_stream_file(the_file, array, n, position_, status_)
!
      implicit none
!
      class(stream_file), intent(in)              :: the_file
      integer(i32), dimension(:,:,:), intent(out) :: array
      integer, intent(in)                         :: n
      integer, intent(in), optional               :: position_
      integer, intent(out), optional              :: status_
!
      call the_file%read_1_int_32_stream_file(array, n, position_, status_)
!
   end subroutine read_3_int_32_stream_file
!
   subroutine read_4_int_32_stream_file(the_file, array, n, position_, status_)
!
      implicit none
!
      class(stream_file), intent(in)                :: the_file
      integer(i32), dimension(:,:,:,:), intent(out) :: array
      integer, intent(in)                           :: n
      integer, intent(in), optional                 :: position_
      integer, intent(out), optional                :: status_
!
      call the_file%read_1_int_32_stream_file(array, n, position_, status_)
!
   end subroutine read_4_int_32_stream_file
!
!  Logicals
!
   subroutine read_0_log_stream_file(the_file, scalar, position_, status_)
!
      implicit none
!
      class(stream_file), intent(in) :: the_file
      logical, intent(out)           :: scalar
      integer, intent(in), optional  :: position_
      integer, intent(out), optional :: status_
!
      call the_file%read_0_log_abstract_stream(scalar, position_, status_)
!
   end subroutine read_0_log_stream_file
!
   subroutine read_1_log_stream_file(the_file, array, n, position_, status_)
!
      implicit none
!
      class(stream_file), intent(in)     :: the_file
      integer, intent(in)                :: n
      logical, dimension(n), intent(out) :: array
      integer, intent(in), optional      :: position_
      integer, intent(out), optional     :: status_
!
      call the_file%read_1_log_abstract_stream(array, n, position_, status_)
!
   end subroutine read_1_log_stream_file
!
!!
!! Wrapper routines for write_real_dp_abstract_stream
!! that accept scalars and rank 1, 2, 3, and 4 arrays
!!
!! The wrapper for scalars and rank 1 are needed because ifort does not permit 
!! a generic procedure pointer to a routine of the parent class
!!
!! see documentation below
!!
!
   subroutine write_0_real_dp_stream_file(the_file, scalar, position_)
!!
!!    write rank 0 real dp
!!    Written by Rolf H. Myhre and Alexander C. Paul, Mar 2020
!!
!!    scalar: real double precision scalar
!!
!!    position_:  optional integer, position to write to in file
!!                positions counted in bytes, starting at 1
!!                default: current file pointer position
!!
      implicit none
!
      class(stream_file), intent(in) :: the_file
      real(dp), intent(in)           :: scalar
      integer, intent(in), optional  :: position_
!
      call the_file%write_0_real_dp_abstract_stream(scalar, position_)
!
   end subroutine write_0_real_dp_stream_file
!
   subroutine write_1_real_dp_stream_file(the_file, array, n, position_)
!!
!!    write rank 1 real dp
!!    Written by Rolf H. Myhre and Alexander C. Paul, Mar 2020
!!
!!    array: real double precision array of length n
!!
!!    n: total length of array (Number of elements)
!!
!!    position_:  optional integer, position to write to in file
!!                positions counted in bytes, starting at 1
!!                default: current file pointer position
!!
      implicit none
!
      class(stream_file), intent(in)     :: the_file
      integer, intent(in)                :: n
      real(dp), dimension(n), intent(in) :: array
      integer, intent(in), optional      :: position_
!
      call the_file%write_1_real_dp_abstract_stream(array, n, position_)
!
   end subroutine write_1_real_dp_stream_file
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
!!    position_:  optional integer, position to write to in file
!!                positions counted in bytes, starting at 1
!!                default: current file pointer position
!!
!
      implicit none
!
      class(stream_file), intent(in)       :: the_file
      real(dp), dimension(:,:), intent(in) :: array
      integer, intent(in)                  :: n
      integer, intent(in), optional        :: position_
!
      call the_file%write_1_real_dp_stream_file(array, n, position_)
!
   end subroutine write_2_real_dp_stream_file
!
   subroutine write_3_real_dp_stream_file(the_file, array, n, position_)
!
      implicit none
!
      class(stream_file), intent(in)         :: the_file
      real(dp), dimension(:,:,:), intent(in) :: array
      integer, intent(in)                    :: n
      integer, intent(in), optional          :: position_
!
      call the_file%write_1_real_dp_stream_file(array, n, position_)
!
   end subroutine write_3_real_dp_stream_file
!
   subroutine write_4_real_dp_stream_file(the_file, array, n, position_)
!
      implicit none
!
      class(stream_file), intent(in)           :: the_file
      real(dp), dimension(:,:,:,:), intent(in) :: array
      integer, intent(in)                      :: n
      integer, intent(in), optional            :: position_
!
      call the_file%write_1_real_dp_stream_file(array, n, position_)
!
   end subroutine write_4_real_dp_stream_file
!
!  Real single precision
!
   subroutine write_0_real_sp_stream_file(the_file, scalar, position_)
!
      implicit none
!
      class(stream_file), intent(in) :: the_file
      real(sp), intent(in)           :: scalar
      integer, intent(in), optional  :: position_
!
      call the_file%write_0_real_sp_abstract_stream(scalar, position_)
!
   end subroutine write_0_real_sp_stream_file
!
   subroutine write_1_real_sp_stream_file(the_file, array, n, position_)
!
      implicit none
!
      class(stream_file), intent(in)     :: the_file
      integer, intent(in)                :: n
      real(sp), dimension(n), intent(in) :: array
      integer, intent(in), optional      :: position_
!
      call the_file%write_1_real_sp_abstract_stream(array, n, position_)
!
   end subroutine write_1_real_sp_stream_file
!
!  Complex
!
   subroutine write_0_complex_dp_stream_file(the_file, scalar, position_)
!
      implicit none
!
      class(stream_file), intent(in) :: the_file
      complex(dp), intent(in)        :: scalar
      integer, intent(in), optional  :: position_
!
      call the_file%write_0_complex_dp_abstract_stream(scalar, position_)
!
   end subroutine write_0_complex_dp_stream_file
!
   subroutine write_1_complex_dp_stream_file(the_file, array, n, position_)
!
      implicit none
!
      class(stream_file), intent(in)        :: the_file
      integer, intent(in)                   :: n
      complex(dp), dimension(n), intent(in) :: array
      integer, intent(in), optional         :: position_
!
      call the_file%write_1_complex_dp_abstract_stream(array, n, position_)
!
   end subroutine write_1_complex_dp_stream_file
!
   subroutine write_2_complex_dp_stream_file(the_file, array, n, position_)
!
      implicit none
!
      class(stream_file), intent(in)          :: the_file
      complex(dp), dimension(:,:), intent(in) :: array
      integer, intent(in)                     :: n
      integer, intent(in), optional           :: position_
!
      call the_file%write_1_complex_dp_stream_file(array, n, position_)
!
   end subroutine write_2_complex_dp_stream_file
!
   subroutine write_3_complex_dp_stream_file(the_file, array, n, position_)
!
      implicit none
!
      class(stream_file), intent(in)            :: the_file
      complex(dp), dimension(:,:,:), intent(in) :: array
      integer, intent(in)                       :: n
      integer, intent(in), optional             :: position_
!
      call the_file%write_1_complex_dp_stream_file(array, n, position_)
!
   end subroutine write_3_complex_dp_stream_file
!
!
   subroutine write_4_complex_dp_stream_file(the_file, array, n, position_)
!
      implicit none
!
      class(stream_file), intent(in)              :: the_file
      complex(dp), dimension(:,:,:,:), intent(in) :: array
      integer, intent(in)                         :: n
      integer, intent(in), optional               :: position_
!
      call the_file%write_1_complex_dp_stream_file(array, n, position_)
!
   end subroutine write_4_complex_dp_stream_file
!
!  64-bit integer
!
   subroutine write_0_int_64_stream_file(the_file, scalar, position_)
!
      implicit none
!
      class(stream_file), intent(in) :: the_file
      integer(i64), intent(in)       :: scalar
      integer, intent(in), optional  :: position_
!
      call the_file%write_0_int_64_abstract_stream(scalar, position_)
!
   end subroutine write_0_int_64_stream_file
!
   subroutine write_1_int_64_stream_file(the_file, array, n, position_)
!
      implicit none
!
      class(stream_file), intent(in)         :: the_file
      integer, intent(in)                    :: n
      integer(i64), dimension(n), intent(in) :: array
      integer, intent(in), optional          :: position_
!
      call the_file%write_1_int_64_abstract_stream(array, n, position_)
!
   end subroutine write_1_int_64_stream_file
!
   subroutine write_2_int_64_stream_file(the_file, array, n, position_)
!
      implicit none
!
      class(stream_file), intent(in)           :: the_file
      integer(i64), dimension(:,:), intent(in) :: array
      integer, intent(in)                      :: n
      integer, intent(in), optional            :: position_
!
      call the_file%write_1_int_64_stream_file(array, n, position_)
!
   end subroutine write_2_int_64_stream_file
!
   subroutine write_3_int_64_stream_file(the_file, array, n, position_)
!
      implicit none
!
      class(stream_file), intent(in)             :: the_file
      integer(i64), dimension(:,:,:), intent(in) :: array
      integer, intent(in)                        :: n
      integer, intent(in), optional              :: position_
!
      call the_file%write_1_int_64_stream_file(array, n, position_)
!
   end subroutine write_3_int_64_stream_file
!
   subroutine write_4_int_64_stream_file(the_file, array, n, position_)
!
      implicit none
!
      class(stream_file), intent(in)               :: the_file
      integer(i64), dimension(:,:,:,:), intent(in) :: array
      integer, intent(in)                          :: n
      integer, intent(in), optional                :: position_
!
      call the_file%write_1_int_64_stream_file(array, n, position_)
!
   end subroutine write_4_int_64_stream_file 
!
!  32-bit integer
!
   subroutine write_0_int_32_stream_file(the_file, scalar, position_)
!
      implicit none
!
      class(stream_file), intent(in) :: the_file
      integer(i32), intent(in)       :: scalar
      integer, intent(in), optional  :: position_
!
      call the_file%write_0_int_32_abstract_stream(scalar, position_)
!
   end subroutine write_0_int_32_stream_file
!
   subroutine write_1_int_32_stream_file(the_file, array, n, position_)
!
      implicit none
!
      class(stream_file), intent(in)         :: the_file
      integer, intent(in)                    :: n
      integer(i32), dimension(n), intent(in) :: array
      integer, intent(in), optional          :: position_
!
      call the_file%write_1_int_32_abstract_stream(array, n, position_)
!
   end subroutine write_1_int_32_stream_file
!
   subroutine write_2_int_32_stream_file(the_file, array, n, position_)
!
      implicit none
!
      class(stream_file), intent(in)           :: the_file
      integer(i32), dimension(:,:), intent(in) :: array
      integer, intent(in)                      :: n
      integer, intent(in), optional            :: position_
!
      call the_file%write_1_int_32_stream_file(array, n, position_)
!
   end subroutine write_2_int_32_stream_file
!
   subroutine write_3_int_32_stream_file(the_file, array, n, position_)
!
      implicit none
!
      class(stream_file), intent(in)             :: the_file
      integer(i32), dimension(:,:,:), intent(in) :: array
      integer, intent(in)                        :: n
      integer, intent(in), optional              :: position_
!
      call the_file%write_1_int_32_stream_file(array, n, position_)
!
   end subroutine write_3_int_32_stream_file
!
   subroutine write_4_int_32_stream_file(the_file, array, n, position_)
!
      implicit none
!
      class(stream_file), intent(in)               :: the_file
      integer(i32), dimension(:,:,:,:), intent(in) :: array
      integer, intent(in)                          :: n
      integer, intent(in), optional                :: position_
!
      call the_file%write_1_int_32_stream_file(array, n, position_)
!
   end subroutine write_4_int_32_stream_file
!
!  Logicals
!
   subroutine write_0_log_stream_file(the_file, scalar, position_)
!
      implicit none
!
      class(stream_file), intent(in) :: the_file
      logical, intent(in)            :: scalar
      integer, intent(in), optional  :: position_
!
      call the_file%write_0_log_abstract_stream(scalar, position_)
!
   end subroutine write_0_log_stream_file
!
   subroutine write_1_log_stream_file(the_file, array, n, position_)
!
      implicit none
!
      class(stream_file), intent(in)    :: the_file
      integer, intent(in)               :: n
      logical, dimension(n), intent(in) :: array
      integer, intent(in), optional     :: position_
!
      call the_file%write_1_log_abstract_stream(array, n, position_)
!
   end subroutine write_1_log_stream_file
!
!
end module stream_file_class
