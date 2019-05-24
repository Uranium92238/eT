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
module direct_file_class
!
!!
!!    Direct access file class module
!!    Written by Rolf H. Myhre, May 2019
!!
!!
!
   use kinds    
   use abstract_file_class, only : abstract_file
   use output_file_class, only : output
   use disk_manager_class, only : disk
!
   type, extends(abstract_file) :: direct_file
!
      integer, private  :: record_dim     ! Number of words per record
      integer, private  :: word_size      ! Size of a word, default is double precision
      integer, private  :: record_length  ! record_dim*word_size
!
   contains
!
!     Open and close
!
      procedure, public :: open_file => open_file_direct_file
      procedure, public :: close_file => close_file_direct_file
!
!     Writer routines
!
      procedure, public :: writer_dp_direct_file
      procedure, public :: writer_dp_dim_direct_file
      procedure, public :: writer_i_direct_file
      procedure, public :: writer_i_dim_direct_file
      generic           :: writer => writer_dp_direct_file, &
                                     writer_dp_dim_direct_file, &
                                     writer_i_direct_file, &
                                     writer_i_dim_direct_file
!
!     Reader routines
!
      procedure, public :: reader_dp_direct_file
      procedure, public :: reader_dp_dim_direct_file
      procedure, public :: reader_i_direct_file
      procedure, public :: reader_i_dim_direct_file
      generic           :: reader => reader_dp_direct_file, &
                                     reader_dp_dim_direct_file, &
                                     reader_i_direct_file, &
                                     reader_i_dim_direct_file
!
   end type direct_file
!
   interface direct_file
!
      procedure new_direct_file
!
   end interface direct_file
!
contains
!
!
   module function new_direct_file(file_name, rec_dim, w_size) result(the_file)
!!
!!    Direct file constructer
!!    Writen by Rolf H. Myhre, May 2019
!!
!!    rec_dim is number of words in each record
!!    w_size (optional) is the size of each word, default is double precision
!!    record length is rec_dim*w_size
!!
      implicit none
!
      type(direct_file) :: the_file
!
      character(len=*), intent(in) :: file_name
      integer, intent(in) :: rec_dim
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
      the_file%file_access = 'direct'
      the_file%file_format = 'unformatted'
!
      if (rec_dim .le. 0) then
         call output%error_msg("Record dimension less than zero for file "//file_name)
      endif
!
      the_file%record_dim = rec_dim
      the_file%record_length = rec_dim*the_file%word_size
!
   end function
!
!
   subroutine open_file_direct_file(the_file, file_action)
!!
!!    Open eT direct file
!!    Written by Rolf Heilemann Myhre, May 2019
!!
      implicit none
!
      class(direct_file)                     :: the_file
      character(len=*), optional, intent(in) :: file_action
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      character(len=20)    :: act
!
      if(present(file_action)) then
         act = trim(file_action)
      else
         act = 'readwrite'
      endif 
!
      open(newunit=the_file%unit, file=the_file%file_name, access=the_file%file_access, &
           action=trim(act), recl=the_file%record_length, status='unknown', form=the_file%file_format, &
           iostat=io_error, iomsg=io_msg)
!
      if (io_error .ne. 0) then 
         call output%error_msg('Error: could not open eT direct file '//trim(the_file%file_name)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
      the_file%file_opened = .true.
!
      call the_file%set_open_file_size()
!
   end subroutine open_file_direct_file
!
!
   subroutine close_file_direct_file(the_file, file_status)
!!
!!    Open the output file
!!    Written by Rolf Heilemann Myhre, May 2019
!!
      implicit none
!
      class(direct_file)                     :: the_file
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
   end subroutine close_file_direct_file
!
!
   module subroutine writer_dp_direct_file(the_file, scalar, record)
!!
!!    Direct file writer, real(dp) scalar
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(direct_file) :: the_file
!
      real(dp), intent(in) :: scalar
      integer, intent(in)  :: record
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      write(the_file%unit, rec=record, iostat=io_error, iomsg=io_msg) scalar
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to write to file: '//the_file%file_name//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine writer_dp_direct_file
!
!
   module subroutine writer_dp_dim_direct_file(the_file, array, record)
!!
!!    Direct file writer, real(dp) array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(direct_file) :: the_file
!
      real(dp), dimension(the_file%record_dim), intent(in) :: array
      integer, intent(in) :: record
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      write(the_file%unit, rec=record, iostat=io_error, iomsg=io_msg) array
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to write to file: '//the_file%file_name//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine writer_dp_dim_direct_file
!
!
   module subroutine writer_i_direct_file(the_file, scalar, record)
!!
!!    Direct file writer, integer scalar
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(direct_file) :: the_file
!
      integer, intent(in) :: scalar
      integer, intent(in) :: record
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      write(the_file%unit, rec=record, iostat=io_error, iomsg=io_msg) scalar
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to write to file: '//the_file%file_name//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine writer_i_direct_file
!
!
   module subroutine writer_i_dim_direct_file(the_file, array, record)
!!
!!    Direct file writer, integer array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(direct_file) :: the_file
!
      integer, dimension(the_file%record_dim), intent(in) :: array
      integer, intent(in) :: record
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      write(the_file%unit, rec=record, iostat=io_error, iomsg=io_msg) array
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to write to file: '//trim(the_file%file_name)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine writer_i_dim_direct_file
!
!
   module subroutine reader_dp_direct_file(the_file, scalar, record)
!!
!!    Direct file reader, real(dp) scalar
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(direct_file) :: the_file
!
      real(dp), intent(out)   :: scalar
      integer, intent(in)     :: record
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      read(the_file%unit, rec=record, iostat=io_error, iomsg=io_msg) scalar
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to read from file: '//trim(the_file%file_name)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine reader_dp_direct_file
!
!
   module subroutine reader_dp_dim_direct_file(the_file, array, record)
!!
!!    Direct file reader, real(dp) array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(direct_file) :: the_file
!
      real(dp), dimension(the_file%record_dim), intent(out) :: array
      integer, intent(in) :: record
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      read(the_file%unit, rec=record, iostat=io_error, iomsg=io_msg) array
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to read from file: '//trim(the_file%file_name)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine reader_dp_dim_direct_file
!
!
   module subroutine reader_i_direct_file(the_file, scalar, record)
!!
!!    Direct file reader, integer scalar
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(direct_file) :: the_file
!
      integer, intent(out) :: scalar
      integer, intent(in)  :: record
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      read(the_file%unit, rec=record, iostat=io_error, iomsg=io_msg) scalar
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to read from file: '//trim(the_file%file_name)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine reader_i_direct_file
!
!
   module subroutine reader_i_dim_direct_file(the_file, array, record)
!!
!!    Direct file reader, integer array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(direct_file) :: the_file
!
      integer, dimension(the_file%record_dim), intent(out) :: array
      integer, intent(in) :: record
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      read(the_file%unit, rec=record, iostat=io_error, iomsg=io_msg) array
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to read from file: '//trim(the_file%file_name)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine reader_i_dim_direct_file
!
!
end module direct_file_class
