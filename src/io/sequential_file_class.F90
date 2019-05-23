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
module sequential_file_class
!
!!
!!    Sequential access file class module
!!    Written by Rolf H. Myhre, May 2019
!!
!!
!
   use kinds    
   use abstract_file_class      
   use output_file_class      
   use disk_manager_class
!
   type, extends(abstract_file) :: sequential_file
!
      integer, private  :: word_size      ! Size of a word, default is double precision
!
   contains
!
!     Open and close
!
      procedure, public :: open_file => open_file_sequential_file
      procedure, public :: close_file => close_file_sequential_file
      procedure, public :: rewind_file => rewind_file_sequential_file
!
!     Writer routines
!
      procedure, public :: writer_dp_sequential_file
      procedure, public :: writer_i6_sequential_file
      procedure, public :: writer_i15_sequential_file
      generic           :: writer => writer_dp_sequential_file, &
                                     writer_i6_sequential_file, &
                                     writer_i15_sequential_file
!
!     Reader routines
!
      procedure, public :: reader_dp_sequential_file
      procedure, public :: reader_i6_sequential_file
      procedure, public :: reader_i15_sequential_file
      generic           :: reader => reader_dp_sequential_file, &
                                     reader_i6_sequential_file, &
                                     reader_i15_sequential_file
!
   end type sequential_file
!
   interface sequential_file
!
      procedure new_sequential_file
!
   end interface sequential_file
!
contains
!
!
   module function new_sequential_file(file_name, w_size) result(the_file)
!!
!!    Sequential file constructer
!!    Writen by Rolf H. Myhre, May 2019
!!
!!    rec_dim is number of words in each record
!!    w_size (optional) is the size of each word, default is double precision
!!    record length is rec_dim*w_size
!!
      implicit none
!
      type(sequential_file) :: the_file
!
      character(len=*), intent(in) :: file_name
      integer, intent(in), optional :: w_size
!
      if (present(w_size)) then
         if (w_size .gt. 0) then
            the_file%word_size = w_size
         else
            call output%error_msg("Word size less than zero for file "//trim(file_name))
         endif
      else
         the_file%word_size = dp
      endif
!
      the_file%file_name = file_name
!
      the_file%file_access = 'sequential'
      the_file%file_format = 'unformatted'
!
   end function
!
!
   subroutine open_file_sequential_file(the_file, file_action, file_pos)
!!
!!    Open eT sequential file
!!    Written by Rolf Heilemann Myhre, May 2019
!!
      implicit none
!
      class(sequential_file)                 :: the_file
      character(len=*), optional, intent(in) :: file_action
      character(len=*), optional, intent(in) :: file_pos
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      character(len=20)    :: act, pos
!
      if(present(file_action)) then
         act = trim(file_action)
      else
         act = 'readwrite'
      endif 
!
      if(present(file_pos)) then
         pos = trim(file_pos)
      else
         pos = 'rewind'
      endif 
!
      open(newunit=the_file%unit, file=the_file%file_name, access=the_file%file_access, &
           action=trim(act), status='unknown', form=the_file%file_format, position=pos, &
           iostat=io_error, iomsg=io_msg)
!
      if (io_error .ne. 0) then 
         call output%error_msg('Error: could not open eT sequential file '//trim(the_file%file_name)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
      the_file%file_opened = .true.
!
      call the_file%set_open_file_size()
!
   end subroutine open_file_sequential_file
!
!
   subroutine close_file_sequential_file(the_file, file_status)
!!
!!    Open the sequential file
!!    Written by Rolf Heilemann Myhre, May 2019
!!
      implicit none
!
      class(sequential_file)                 :: the_file
      character(len=*), optional, intent(in) :: file_status
!
      integer  :: file_change
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      character(len=20)    :: stat
!
      if(present(file_status)) then
         stat = trim(file_status)
      else
         stat = 'keep'
      endif 
!
      close(the_file%unit, iostat=io_error, iomsg=io_msg, status=trim(stat))
!
      if (io_error .ne. 0) then 
         call output%error_msg('Error: could not close eT file '//trim(the_file%file_name)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
      the_file%file_opened = .false.
!
      file_change = the_file%get_file_change()
      call disk%update(file_change, the_file%file_name)
!
   end subroutine close_file_sequential_file
!
!
   subroutine rewind_file_sequential_file(the_file)
!!
!!    Rewind the sequential file
!!    Written by Rolf Heilemann Myhre, May 2019
!!
      implicit none
!
      class(sequential_file)                 :: the_file
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      rewind(the_file%unit, iostat=io_error, iomsg=io_msg)
!
      if (io_error .ne. 0) then 
         call output%error_msg('Error: could not rewind eT file '//trim(the_file%file_name)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine rewind_file_sequential_file
!
!
   module subroutine writer_dp_sequential_file(the_file, array, n)
!!
!!    Sequential file writer, real double precision
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(sequential_file) :: the_file
!
      integer, intent(in)                 :: n
      real(dp), dimension(n), intent(in)  :: array
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      write(the_file%unit, iostat=io_error, iomsg=io_msg) array
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to write to file: '//the_file%file_name//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine writer_dp_sequential_file
!
!
   module subroutine writer_i6_sequential_file(the_file, array, n)
!!
!!    Sequential file writer, integer 32
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(sequential_file) :: the_file
!
      integer, intent(in)                    :: n
      integer(i6), dimension(n), intent(in)  :: array
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      write(the_file%unit, iostat=io_error, iomsg=io_msg) array
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to write to file: '//the_file%file_name//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine writer_i6_sequential_file
!
!
   module subroutine writer_i15_sequential_file(the_file, array, n)
!!
!!    Sequential file writer, integer 64
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(sequential_file) :: the_file
!
      integer, intent(in)                    :: n
      integer(i15), dimension(n), intent(in) :: array
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      write(the_file%unit, iostat=io_error, iomsg=io_msg) array
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to write to file: '//trim(the_file%file_name)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine writer_i15_sequential_file
!
!
   module subroutine reader_dp_sequential_file(the_file, array, n)
!!
!!    Sequential file reader, real double precision
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(sequential_file) :: the_file
!
      integer, intent(in)                 :: n
      real(dp), dimension(n), intent(out) :: array
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      read(the_file%unit, iostat=io_error, iomsg=io_msg) array
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to read from file: '//trim(the_file%file_name)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine reader_dp_sequential_file
!
!
   module subroutine reader_i6_sequential_file(the_file, array, n)
!!
!!    Sequential file reader, integer 32
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(sequential_file) :: the_file
!
      integer, intent(in)                    :: n
      integer(i6), dimension(n), intent(out) :: array
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      read(the_file%unit, iostat=io_error, iomsg=io_msg) array
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to read from file: '//trim(the_file%file_name)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine reader_i6_sequential_file
!
!
   module subroutine reader_i15_sequential_file(the_file, array, n)
!!
!!    Sequential file reader, integer 64
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(sequential_file) :: the_file
!
      integer, intent(in)                     :: n
      integer(i15), dimension(n), intent(out) :: array
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      read(the_file%unit, iostat=io_error, iomsg=io_msg) array
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to read from file: '//trim(the_file%file_name)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine reader_i15_sequential_file
!
!
end module sequential_file_class
