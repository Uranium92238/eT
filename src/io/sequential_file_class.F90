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
   use abstract_file_class, only : abstract_file
   use output_file_class, only : output
   use disk_manager_class, only : disk
!
   type, extends(abstract_file) :: sequential_file
!
      integer, private  :: word_size      ! Size of a word, default is double precision
!
   contains
!
!     Open and close
!
      procedure, public :: open_   => open__sequential_file
      procedure, public :: close_  => close__sequential_file
      procedure, public :: rewind_ => rewind__sequential_file
      procedure, public :: skip    => skip_sequential_file
!
!     Write routines
!
      procedure, public :: write_dp_sequential_file
      procedure, public :: write_dp_1_sequential_file
      procedure, public :: write_dp_2_sequential_file
      procedure, public :: write_dp_3_sequential_file
      procedure, public :: write_dp_4_sequential_file
      procedure, public :: write_i_sequential_file
      procedure, public :: write_i_1_sequential_file
      procedure, public :: write_i_2_sequential_file
      procedure, public :: write_i_3_sequential_file
      procedure, public :: write_i_4_sequential_file
      generic           :: write_ => write_dp_sequential_file,   &
                                     write_dp_1_sequential_file, &
                                     write_dp_2_sequential_file, &
                                     write_dp_3_sequential_file, &
                                     write_dp_4_sequential_file, &
                                     write_i_sequential_file,    &
                                     write_i_1_sequential_file,  &
                                     write_i_2_sequential_file,  &
                                     write_i_3_sequential_file,  &
                                     write_i_4_sequential_file
!
!     Read routines
!
      procedure, public :: read_dp_sequential_file
      procedure, public :: read_dp_1_sequential_file
      procedure, public :: read_dp_2_sequential_file
      procedure, public :: read_dp_3_sequential_file
      procedure, public :: read_dp_4_sequential_file
      procedure, public :: read_i_sequential_file
      procedure, public :: read_i_1_sequential_file
      procedure, public :: read_i_2_sequential_file
      procedure, public :: read_i_3_sequential_file
      procedure, public :: read_i_4_sequential_file
      generic           :: read_ => read_dp_sequential_file, &
                                     read_dp_1_sequential_file, &
                                     read_dp_2_sequential_file, &
                                     read_dp_3_sequential_file, &
                                     read_dp_4_sequential_file, &
                                     read_i_sequential_file, &
                                     read_i_1_sequential_file, &
                                     read_i_2_sequential_file, &
                                     read_i_3_sequential_file, &
                                     read_i_4_sequential_file
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
   module function new_sequential_file(name_) result(the_file)
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
      character(len=*), intent(in) :: name_
!
      the_file%name_ = name_
!
      the_file%access_ = 'sequential'
      the_file%format_ = 'unformatted'
!
   end function
!
!
   subroutine open__sequential_file(the_file, file_action, file_pos)
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
      open(newunit=the_file%unit, file=the_file%name_, access=the_file%access_, &
           action=trim(act), status='unknown', form=the_file%format_, position=pos, &
           iostat=io_error, iomsg=io_msg)
!
      if (io_error .ne. 0) then 
         call output%error_msg('Error: could not open eT sequential file '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
      the_file%opened = .true.
!
      call the_file%set_open_size()
!
   end subroutine open__sequential_file
!
!
   subroutine close__sequential_file(the_file, file_status)
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
         call output%error_msg('Error: could not close eT file '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
      the_file%opened = .false.
!
      file_change = the_file%get_change()
      call disk%update(file_change, the_file%name_)
!
   end subroutine close__sequential_file
!
!
   subroutine rewind__sequential_file(the_file)
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
         call output%error_msg('Error: could not rewind eT file '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine rewind__sequential_file
!
!
   subroutine skip_sequential_file(the_file, jump)
!!
!!    Skip sequential file
!!    Written by Rolf Heilemann Myhre, May 2019
!!    Skip jump number of lines. Default: 1
!!
      implicit none
!
      class(sequential_file), intent(in)  :: the_file
      integer, optional,  intent(in)      :: jump
!
      integer :: i,skips
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      if(present(jump)) then
         skips = jump
      else
         skips = 1
      endif
!
      do i = 1,skips
         read(the_file%unit, iostat=io_error, iomsg=io_msg) 
!
         if (io_error .ne. 0) then 
            call output%error_msg('Error: could not skip eT sequential file '//trim(the_file%name_)//&
                                 &'. Error message: '//trim(io_msg))
         endif
!
      enddo
!
   end subroutine skip_sequential_file
!
!
   module subroutine write_dp_sequential_file(the_file, scalar)
!!
!!    Sequential file write, real(dp0 scalar
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(sequential_file), intent(in) :: the_file
!
      real(dp), intent(in) :: scalar
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      write(the_file%unit, iostat=io_error, iomsg=io_msg) scalar
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to write to file: '//the_file%name_//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine write_dp_sequential_file
!
!
   module subroutine write_dp_1_sequential_file(the_file, array, n)
!!
!!    Sequential file write, real(dp) array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(sequential_file), intent(in)  :: the_file
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
         call output%error_msg('Failed to write to file: '//the_file%name_//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine write_dp_1_sequential_file
!
!
   module subroutine write_dp_2_sequential_file(the_file, array, n)
!!
!!    Sequential file write, real(dp) array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(sequential_file), intent(in)     :: the_file
      integer, intent(in)                    :: n
      real(dp), dimension(:,:), intent(in)   :: array
!
      call the_file%write_dp_1_sequential_file(array, n)
!
   end subroutine write_dp_2_sequential_file
!
!
   module subroutine write_dp_3_sequential_file(the_file, array, n)
!!
!!    Sequential file write, real(dp) array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(sequential_file), intent(in)     :: the_file
      integer, intent(in)                    :: n
      real(dp), dimension(:,:,:), intent(in) :: array
!
      call the_file%write_dp_1_sequential_file(array, n)
!
   end subroutine write_dp_3_sequential_file
!
!
   module subroutine write_dp_4_sequential_file(the_file, array, n)
!!
!!    Sequential file write, real(dp) array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(sequential_file), intent(in)        :: the_file
      integer, intent(in)                       :: n
      real(dp), dimension(:,:,:,:), intent(in)  :: array
!
      call the_file%write_dp_1_sequential_file(array, n)
!
   end subroutine write_dp_4_sequential_file
!
!
   module subroutine write_i_sequential_file(the_file, scalar)
!!
!!    Sequential file write, integer scalar
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(sequential_file), intent(in) :: the_file
!
      integer, intent(in)  :: scalar
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      write(the_file%unit, iostat=io_error, iomsg=io_msg) scalar
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to write to file: '//the_file%name_//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine write_i_sequential_file
!
!
   module subroutine write_i_1_sequential_file(the_file, array, n)
!!
!!    Sequential file write, integer array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(sequential_file), intent(in)  :: the_file
!
      integer, intent(in)                 :: n
      integer, dimension(n), intent(in)   :: array
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      write(the_file%unit, iostat=io_error, iomsg=io_msg) array
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to write to file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine write_i_1_sequential_file
!
!
   module subroutine write_i_2_sequential_file(the_file, array, n)
!!
!!    Sequential file write, integer array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(sequential_file), intent(in)  :: the_file
      integer, intent(in)                 :: n
      integer, dimension(:,:), intent(in) :: array
!
      call the_file%write_i_1_sequential_file(array, n)
!
   end subroutine write_i_2_sequential_file
!
!
   module subroutine write_i_3_sequential_file(the_file, array, n)
!!
!!    Sequential file write, integer array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(sequential_file), intent(in)     :: the_file
      integer, intent(in)                    :: n
      integer, dimension(:,:,:), intent(in)  :: array
!
      call the_file%write_i_1_sequential_file(array, n)
!
   end subroutine write_i_3_sequential_file
!
!
   module subroutine write_i_4_sequential_file(the_file, array, n)
!!
!!    Sequential file write, integer array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(sequential_file), intent(in)        :: the_file
      integer, intent(in)                       :: n
      integer, dimension(:,:,:,:), intent(in)   :: array
!
      call the_file%write_i_1_sequential_file(array, n)
!
   end subroutine write_i_4_sequential_file
!
!
   module subroutine read_dp_sequential_file(the_file, scalar)
!!
!!    Sequential file read, real(dp) scalar
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(sequential_file), intent(in)  :: the_file
!
      real(dp), intent(out) :: scalar
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      read(the_file%unit, iostat=io_error, iomsg=io_msg) scalar
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to read from file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine read_dp_sequential_file
!
!
   module subroutine read_dp_1_sequential_file(the_file, array, n)
!!
!!    Sequential file read, real(dp) array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(sequential_file), intent(in)  :: the_file
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
         call output%error_msg('Failed to read from file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine read_dp_1_sequential_file
!
!
   module subroutine read_dp_2_sequential_file(the_file, array, n)
!!
!!    Sequential file read, real(dp) array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(sequential_file), intent(in)     :: the_file
      integer, intent(in)                    :: n
      real(dp), dimension(:,:), intent(out)  :: array
!
      call the_file%read_dp_1_sequential_file(array, n)
!
   end subroutine read_dp_2_sequential_file
!
!
   module subroutine read_dp_3_sequential_file(the_file, array, n)
!!
!!    Sequential file read, real(dp) array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(sequential_file), intent(in)        :: the_file
      integer, intent(in)                       :: n
      real(dp), dimension(:,:,:), intent(out)   :: array
!
      call the_file%read_dp_1_sequential_file(array, n)
!
   end subroutine read_dp_3_sequential_file
!
!
   module subroutine read_dp_4_sequential_file(the_file, array, n)
!!
!!    Sequential file read, real(dp) array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(sequential_file), intent(in)        :: the_file
      integer, intent(in)                       :: n
      real(dp), dimension(:,:,:,:), intent(out) :: array
!
      call the_file%read_dp_1_sequential_file(array, n)
!
   end subroutine read_dp_4_sequential_file
!
!
   module subroutine read_i_sequential_file(the_file, scalar)
!!
!!    Sequential file read, integer scalar
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(sequential_file), intent(in) :: the_file
!
      integer, intent(out) :: scalar
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      read(the_file%unit, iostat=io_error, iomsg=io_msg) scalar
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to read from file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine read_i_sequential_file
!
!
   module subroutine read_i_1_sequential_file(the_file, array, n)
!!
!!    Sequential file read, integer array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(sequential_file), intent(in) :: the_file
!
      integer, intent(in)                 :: n
      integer, dimension(n), intent(out)  :: array
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      read(the_file%unit, iostat=io_error, iomsg=io_msg) array
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to read from file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine read_i_1_sequential_file
!
!
   module subroutine read_i_2_sequential_file(the_file, array, n)
!!
!!    Sequential file read, integer array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(sequential_file), intent(in)     :: the_file
      integer, intent(in)                    :: n
      integer, dimension(:,:), intent(out)   :: array
!
      call the_file%read_i_1_sequential_file(array, n)
!
   end subroutine read_i_2_sequential_file
!
!
   module subroutine read_i_3_sequential_file(the_file, array, n)
!!
!!    Sequential file read, integer array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(sequential_file), intent(in)     :: the_file
      integer, intent(in)                    :: n
      integer, dimension(:,:,:), intent(out) :: array
!
      call the_file%read_i_1_sequential_file(array, n)
!
   end subroutine read_i_3_sequential_file
!
!
   module subroutine read_i_4_sequential_file(the_file, array, n)
!!
!!    Sequential file read, integer array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(sequential_file), intent(in)        :: the_file
      integer, intent(in)                       :: n
      integer, dimension(:,:,:,:), intent(out)  :: array
!
      call the_file%read_i_1_sequential_file(array, n)
!
   end subroutine read_i_4_sequential_file
!
!
end module sequential_file_class
