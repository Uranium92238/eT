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
   use global_out, only : output
   use disk_manager_class, only : disk
!
   type, extends(abstract_file) :: sequential_file
!
   contains
!
!     Open and close
!
      procedure, public :: open_   => open_sequential_file
      procedure, public :: close_  => close_sequential_file
      procedure, public :: delete_ => delete_sequential_file
      procedure, public :: rewind_ => rewind_sequential_file
      procedure, public :: skip    => skip_sequential_file
!
!     Write routines
!
      procedure, public :: write_blank_sequential_file
!
      procedure, public :: write_r_sequential_file
      procedure, public :: write_r_1_sequential_file
      procedure, public :: write_r_2_sequential_file
      procedure, public :: write_r_3_sequential_file
      procedure, public :: write_r_4_sequential_file
!
      procedure, public :: write_c_sequential_file
      procedure, public :: write_c_1_sequential_file
!
      procedure, public :: write_i_sequential_file
      procedure, public :: write_i_i_sequential_file
      procedure, public :: write_i_1_sequential_file
      procedure, public :: write_i_2_sequential_file
      procedure, public :: write_i_3_sequential_file
      procedure, public :: write_i_4_sequential_file
!
      procedure, public :: write_l_sequential_file
      procedure, public :: write_l_1_sequential_file
!
      procedure, public :: write_char_sequential_file
!
      generic           :: write_ => write_blank_sequential_file, &
                                     write_r_sequential_file,     &
                                     write_r_1_sequential_file,   &
                                     write_r_2_sequential_file,   &
                                     write_r_3_sequential_file,   &
                                     write_r_4_sequential_file,   &
                                     write_c_sequential_file,     &
                                     write_c_1_sequential_file,   &
                                     write_i_sequential_file,     &
                                     write_i_i_sequential_file,   &
                                     write_i_1_sequential_file,   &
                                     write_i_2_sequential_file,   &
                                     write_i_3_sequential_file,   &
                                     write_i_4_sequential_file,   &
                                     write_l_sequential_file,     &
                                     write_l_1_sequential_file,   &
                                     write_char_sequential_file
!
!     Read routines
!
      procedure, public :: read_blank_sequential_file
!
      procedure, public :: read_r_sequential_file
      procedure, public :: read_r_1_sequential_file
      procedure, public :: read_r_2_sequential_file
      procedure, public :: read_r_3_sequential_file
      procedure, public :: read_r_4_sequential_file
!
      procedure, public :: read_c_sequential_file
      procedure, public :: read_c_1_sequential_file
      procedure, public :: read_c_2_sequential_file
      procedure, public :: read_c_3_sequential_file
      procedure, public :: read_c_4_sequential_file
!
      procedure, public :: read_i_sequential_file
      procedure, public :: read_i_i_sequential_file
      procedure, public :: read_i_1_sequential_file
      procedure, public :: read_i_2_sequential_file
      procedure, public :: read_i_3_sequential_file
      procedure, public :: read_i_4_sequential_file
!
      procedure, public :: read_char_sequential_file
!
      generic           :: read_ => read_blank_sequential_file, &
                                    read_r_sequential_file,     &
                                    read_r_1_sequential_file,   &
                                    read_r_2_sequential_file,   &
                                    read_r_3_sequential_file,   &
                                    read_r_4_sequential_file,   &
                                    read_c_sequential_file,     &
                                    read_c_1_sequential_file,   &
                                    read_c_2_sequential_file,   &
                                    read_c_3_sequential_file,   &
                                    read_c_4_sequential_file,   &
                                    read_i_sequential_file,     &
                                    read_i_i_sequential_file,   &
                                    read_i_1_sequential_file,   &
                                    read_i_2_sequential_file,   &
                                    read_i_3_sequential_file,   &
                                    read_i_4_sequential_file,   &
                                    read_char_sequential_file
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
   function new_sequential_file(name_, format_) result(the_file)
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
      character(len=*), optional, intent(in) :: format_
!
      the_file%name_ = name_
!
      the_file%access_ = 'sequential'
!
      if(present(format_)) then
         if(format_ .eq. 'formatted' .or. format_ .eq. 'unformatted') then
            the_file%format_ = format_
         else
            call output%error_msg('Wrong format specifier in eT sequential file '//the_file%name_ // &
                                 & ' ,format specifier: '//format_)
         endif
      else
         the_file%format_ = 'unformatted'
      endif
!
      the_file%is_open = .false.
      the_file%unit = -1
!
   end function new_sequential_file
!
!
   subroutine open_sequential_file(the_file, file_action, position_)
!!
!!    Open eT sequential file
!!    Written by Rolf Heilemann Myhre, May 2019
!!
      implicit none
!
      class(sequential_file)                 :: the_file
      character(len=*), optional, intent(in) :: file_action
      character(len=*), optional, intent(in) :: position_
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      character(len=20)    :: pos
!
      if(present(file_action)) then
         the_file%action_ = trim(file_action)
      else
         the_file%action_ = 'readwrite'
      endif 
!
      if(present(position_)) then
         pos = trim(position_)
      else
         pos = 'rewind'
      endif 
!
      if (the_file%is_open) then
!
         call output%error_msg(trim(the_file%name_)//' is already open')
!
      endif
!
      open(newunit=the_file%unit, file=the_file%name_, access=the_file%access_, &
           action=the_file%action_, status='unknown', form=the_file%format_, &
           position=pos, iostat=io_error, iomsg=io_msg)
!
      if (io_error .ne. 0) then 
         call output%error_msg('Error: could not open eT sequential file '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
      the_file%is_open = .true.
!
      call the_file%set_open_size()
!
   end subroutine open_sequential_file
!
!
   subroutine close_sequential_file(the_file, file_status)
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
      if (.not. the_file%is_open) then
         call output%error_msg(trim(the_file%name_)//' already closed')
      end if
!
      close(the_file%unit, iostat=io_error, iomsg=io_msg, status=trim(stat))
!
      if (io_error .ne. 0) then 
         call output%error_msg('could not close eT file '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
      file_change = the_file%get_change()
      call disk%update(file_change, the_file%name_)
!
      the_file%is_open = .false.
      the_file%unit = -1
      the_file%action_ = 'unknown'
!
   end subroutine close_sequential_file
!
!
   subroutine rewind_sequential_file(the_file)
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
   end subroutine rewind_sequential_file
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
   subroutine delete_sequential_file(the_file)
!!
!!    Delete file
!!    Written by Rolf Heilemann Myhre, Aug 2019
!!
      implicit none
!
      class(sequential_file) :: the_file
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      if(the_file%is_open) then
!
         close(the_file%unit, iostat=io_error, iomsg=io_msg, status='delete')
!
         if (io_error .ne. 0) then 
            call output%error_msg('Error: could not delete eT file '//trim(the_file%name_)//&
                                 &'. Error message: '//trim(io_msg))
         endif
!
      else
!
         open(newunit=the_file%unit, file=the_file%name_, access=the_file%access_, &
              action='write', iostat=io_error, iomsg=io_msg)
!
         if (io_error .ne. 0) then 
            call output%error_msg('Error: could not open eT file '//trim(the_file%name_)//&
                                 &'. Error message: '//trim(io_msg))
         endif
!
         close(the_file%unit, iostat=io_error, iomsg=io_msg, status='delete')
!
         if (io_error .ne. 0) then 
            call output%error_msg('Error: could not delete eT file '//trim(the_file%name_)//&
                                 &'. Error message: '//trim(io_msg))
         endif
!
      endif
!
   end subroutine delete_sequential_file
!
!
   subroutine write_blank_sequential_file(the_file)
!!
!!    Sequential file write, blank line
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(sequential_file), intent(in) :: the_file
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      if (the_file%format_ .eq. 'unformatted') then
         write(the_file%unit, iostat=io_error, iomsg=io_msg) 
      else
         write(the_file%unit, *, iostat=io_error, iomsg=io_msg) 
      endif
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to write to file: '//the_file%name_//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine write_blank_sequential_file
!
!
   subroutine write_r_sequential_file(the_file, scalar)
!!
!!    Sequential file write, real(dp) scalar
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
      if (the_file%format_ .eq. 'unformatted') then
         write(the_file%unit, iostat=io_error, iomsg=io_msg) scalar
      else
         write(the_file%unit, *, iostat=io_error, iomsg=io_msg) scalar
      endif
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to write to file: '//the_file%name_//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine write_r_sequential_file
!
!
   subroutine write_r_1_sequential_file(the_file, array, n)
!!
!!    Sequential file write, real(dp) 1 dimensional array
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
      if (the_file%format_ .eq. 'unformatted') then
         write(the_file%unit, iostat=io_error, iomsg=io_msg) array
      else
         write(the_file%unit, *, iostat=io_error, iomsg=io_msg) array
      endif
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to write to file: '//the_file%name_//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine write_r_1_sequential_file
!
!
   subroutine write_r_2_sequential_file(the_file, array, n)
!!
!!    Sequential file write, real(dp) 2 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(sequential_file), intent(in)     :: the_file
      integer, intent(in)                    :: n
      real(dp), dimension(:,:), intent(in)   :: array
!
      call the_file%write_r_1_sequential_file(array, n)
!
   end subroutine write_r_2_sequential_file
!
!
   subroutine write_r_3_sequential_file(the_file, array, n)
!!
!!    Sequential file write, real(dp) 3 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(sequential_file), intent(in)     :: the_file
      integer, intent(in)                    :: n
      real(dp), dimension(:,:,:), intent(in) :: array
!
      call the_file%write_r_1_sequential_file(array, n)
!
   end subroutine write_r_3_sequential_file
!
!
   subroutine write_r_4_sequential_file(the_file, array, n)
!!
!!    Sequential file write, real(dp) 4 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(sequential_file), intent(in)        :: the_file
      integer, intent(in)                       :: n
      real(dp), dimension(:,:,:,:), intent(in)  :: array
!
      call the_file%write_r_1_sequential_file(array, n)
!
   end subroutine write_r_4_sequential_file
!
!
   subroutine write_c_sequential_file(the_file, scalar)
!!
!!    Sequential file write, complex(dp) scalar
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(sequential_file), intent(in) :: the_file
!
      complex(dp), intent(in) :: scalar
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      if (the_file%format_ .eq. 'unformatted') then
         write(the_file%unit, iostat=io_error, iomsg=io_msg) scalar
      else
         write(the_file%unit, *, iostat=io_error, iomsg=io_msg) scalar
      endif
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to write to file: '//the_file%name_//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine write_c_sequential_file
!
!
   subroutine write_c_1_sequential_file(the_file, array, n)
!!
!!    Sequential file write, complex(dp) 1 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(sequential_file), intent(in)  :: the_file
!
      integer, intent(in)                 :: n
      complex(dp), dimension(n), intent(in)  :: array
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      if (the_file%format_ .eq. 'unformatted') then
         write(the_file%unit, iostat=io_error, iomsg=io_msg) array
      else
         write(the_file%unit, *, iostat=io_error, iomsg=io_msg) array
      endif
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to write to file: '//the_file%name_//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine write_c_1_sequential_file
!
!
   subroutine write_c_2_sequential_file(the_file, array, n)
!!
!!    Sequential file write, complex(dp) 2 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(sequential_file), intent(in)      :: the_file
      integer, intent(in)                     :: n
      complex(dp), dimension(:,:), intent(in) :: array
!
      call the_file%write_c_1_sequential_file(array, n)
!
   end subroutine write_c_2_sequential_file
!
!
   subroutine write_c_3_sequential_file(the_file, array, n)
!!
!!    Sequential file write, complex(dp) 3 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(sequential_file), intent(in)        :: the_file
      integer, intent(in)                       :: n
      complex(dp), dimension(:,:,:), intent(in) :: array
!
      call the_file%write_c_1_sequential_file(array, n)
!
   end subroutine write_c_3_sequential_file
!
!
   subroutine write_c_4_sequential_file(the_file, array, n)
!!
!!    Sequential file write, complex(dp) 4 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(sequential_file), intent(in)          :: the_file
      integer, intent(in)                         :: n
      complex(dp), dimension(:,:,:,:), intent(in) :: array
!
      call the_file%write_c_1_sequential_file(array, n)
!
   end subroutine write_c_4_sequential_file
!
!
   subroutine write_i_sequential_file(the_file, scalar)
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
      if (the_file%format_ .eq. 'unformatted') then
         write(the_file%unit, iostat=io_error, iomsg=io_msg) scalar
      else
         write(the_file%unit, *, iostat=io_error, iomsg=io_msg) scalar
      endif
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to write to file: '//the_file%name_//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine write_i_sequential_file
!
!
   subroutine write_i_i_sequential_file(the_file, scalar1, scalar2)
!!
!!    Sequential file write, two integer scalars
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(sequential_file), intent(in) :: the_file
!
      integer, intent(in)  :: scalar1, scalar2
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      if (the_file%format_ .eq. 'unformatted') then
         write(the_file%unit, iostat=io_error, iomsg=io_msg) scalar1, scalar2
      else
         write(the_file%unit, *, iostat=io_error, iomsg=io_msg) scalar1, scalar2
      endif
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to write to file: '//the_file%name_//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine write_i_i_sequential_file
!
!
   subroutine write_i_1_sequential_file(the_file, array, n)
!!
!!    Sequential file write, integer 1 dimensional array
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
      if (the_file%format_ .eq. 'unformatted') then
         write(the_file%unit, iostat=io_error, iomsg=io_msg) array
      else
         write(the_file%unit, *, iostat=io_error, iomsg=io_msg) array
      endif
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to write to file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine write_i_1_sequential_file
!
!
   subroutine write_i_2_sequential_file(the_file, array, n)
!!
!!    Sequential file write, integer 2 dimensional array
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
   subroutine write_i_3_sequential_file(the_file, array, n)
!!
!!    Sequential file write, integer 3 dimensional array
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
   subroutine write_i_4_sequential_file(the_file, array, n)
!!
!!    Sequential file write, integer 4 dimensional array
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
   subroutine write_l_sequential_file(the_file, scalar)
!!
!!    Sequential file write, logical scalar
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(sequential_file), intent(in) :: the_file
!
      logical, intent(in)  :: scalar
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      if (the_file%format_ .eq. 'unformatted') then
         write(the_file%unit, iostat=io_error, iomsg=io_msg) scalar
      else
         write(the_file%unit, *, iostat=io_error, iomsg=io_msg) scalar
      endif
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to write to file: '//the_file%name_//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine write_l_sequential_file
!
!
   subroutine write_l_1_sequential_file(the_file, array, n)
!!
!!    Sequential file write, logical array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(sequential_file), intent(in)  :: the_file
!
      integer, intent(in)                 :: n
      logical, dimension(n), intent(in)   :: array
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      if (the_file%format_ .eq. 'unformatted') then
         write(the_file%unit, iostat=io_error, iomsg=io_msg) array
      else
         write(the_file%unit, *, iostat=io_error, iomsg=io_msg) array
      endif
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to write to file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine write_l_1_sequential_file
!
!
   subroutine write_char_sequential_file(the_file, string, format_string)
!!
!!    Sequential file write, character
!!    Written by Rolf H. Myhre, May 2019
!!    string:  character string to write from
!!    format_string: optional format string
!!
      implicit none
!
      class(sequential_file), intent(in) :: the_file
!
      character(len=*), intent(in)           :: string
      character(len=*), intent(in), optional :: format_string
!
      integer              :: io_error
      character(len=100)   :: io_msg
      character(len=100)   :: fstring
!
      if(present(format_string)) then
         if (the_file%format_ .eq. 'unformatted') then
            call output%error_msg(trim(the_file%name_)//' is unformatted')
         endif
         fstring = trim(format_string)
      else
         fstring = '*'
      endif
!
      if (the_file%format_ .eq. 'unformatted') then
         write(the_file%unit, iostat=io_error, iomsg=io_msg) string
      else
         write(the_file%unit, trim(fstring), iostat=io_error, iomsg=io_msg) string
      endif
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to write to file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine write_char_sequential_file
!
!
   subroutine read_blank_sequential_file(the_file)
!!
!!    Sequential file read, blank line
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(sequential_file), intent(in)  :: the_file
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      if (the_file%format_ .eq. 'unformatted') then
         read(the_file%unit, iostat=io_error, iomsg=io_msg)
      else
         read(the_file%unit, *, iostat=io_error, iomsg=io_msg)
      endif
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to read from file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine read_blank_sequential_file
!
!
   subroutine read_r_sequential_file(the_file, scalar)
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
      if (the_file%format_ .eq. 'unformatted') then
         read(the_file%unit, iostat=io_error, iomsg=io_msg) scalar
      else
         read(the_file%unit, *, iostat=io_error, iomsg=io_msg) scalar
      endif
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to read from file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine read_r_sequential_file
!
!
   subroutine read_r_1_sequential_file(the_file, array, n)
!!
!!    Sequential file read, real(dp) 1 dimensional array
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
      if (the_file%format_ .eq. 'unformatted') then
         read(the_file%unit, iostat=io_error, iomsg=io_msg) array
      else
         read(the_file%unit, *, iostat=io_error, iomsg=io_msg) array
      endif
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to read from file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine read_r_1_sequential_file
!
!
   subroutine read_r_2_sequential_file(the_file, array, n)
!!
!!    Sequential file read, real(dp) 2 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(sequential_file), intent(in)     :: the_file
      integer, intent(in)                    :: n
      real(dp), dimension(:,:), intent(out)  :: array
!
      call the_file%read_r_1_sequential_file(array, n)
!
   end subroutine read_r_2_sequential_file
!
!
   subroutine read_r_3_sequential_file(the_file, array, n)
!!
!!    Sequential file read, real(dp) 3 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(sequential_file), intent(in)        :: the_file
      integer, intent(in)                       :: n
      real(dp), dimension(:,:,:), intent(out)   :: array
!
      call the_file%read_r_1_sequential_file(array, n)
!
   end subroutine read_r_3_sequential_file
!
!
   subroutine read_r_4_sequential_file(the_file, array, n)
!!
!!    Sequential file read, real(dp) 4 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(sequential_file), intent(in)        :: the_file
      integer, intent(in)                       :: n
      real(dp), dimension(:,:,:,:), intent(out) :: array
!
      call the_file%read_r_1_sequential_file(array, n)
!
   end subroutine read_r_4_sequential_file
!
!
   subroutine read_c_sequential_file(the_file, scalar)
!!
!!    Sequential file read, complex(dp) scalar
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(sequential_file), intent(in)  :: the_file
!
      complex(dp), intent(out) :: scalar
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      if (the_file%format_ .eq. 'unformatted') then
         read(the_file%unit, iostat=io_error, iomsg=io_msg) scalar
      else
         read(the_file%unit, *, iostat=io_error, iomsg=io_msg) scalar
      endif
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to read from file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine read_c_sequential_file
!
!
   subroutine read_c_1_sequential_file(the_file, array, n)
!!
!!    Sequential file read, complex(dp) 1 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(sequential_file), intent(in)  :: the_file
!
      integer, intent(in)                    :: n
      complex(dp), dimension(n), intent(out) :: array
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      if (the_file%format_ .eq. 'unformatted') then
         read(the_file%unit, iostat=io_error, iomsg=io_msg) array
      else
         read(the_file%unit, *, iostat=io_error, iomsg=io_msg) array
      endif
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to read from file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine read_c_1_sequential_file
!
!
   subroutine read_c_2_sequential_file(the_file, array, n)
!!
!!    Sequential file read, complex(dp) 2 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(sequential_file), intent(in)       :: the_file
      integer, intent(in)                      :: n
      complex(dp), dimension(:,:), intent(out) :: array
!
      call the_file%read_c_1_sequential_file(array, n)
!
   end subroutine read_c_2_sequential_file
!
!
   subroutine read_c_3_sequential_file(the_file, array, n)
!!
!!    Sequential file read, complex(dp) 3 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(sequential_file), intent(in)         :: the_file
      integer, intent(in)                        :: n
      complex(dp), dimension(:,:,:), intent(out) :: array
!
      call the_file%read_c_1_sequential_file(array, n)
!
   end subroutine read_c_3_sequential_file
!
!
   subroutine read_c_4_sequential_file(the_file, array, n)
!!
!!    Sequential file read, complex(dp) 4 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(sequential_file), intent(in)           :: the_file
      integer, intent(in)                          :: n
      complex(dp), dimension(:,:,:,:), intent(out) :: array
!
      call the_file%read_c_1_sequential_file(array, n)
!
   end subroutine read_c_4_sequential_file
!
!
   subroutine read_i_sequential_file(the_file, scalar)
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
      if (the_file%format_ .eq. 'unformatted') then
         read(the_file%unit, iostat=io_error, iomsg=io_msg) scalar
      else
         read(the_file%unit, *, iostat=io_error, iomsg=io_msg) scalar
      endif
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to read from file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine read_i_sequential_file
!
!
   subroutine read_i_i_sequential_file(the_file, scalar1, scalar2)
!!
!!    Sequential file read, two integer scalars
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(sequential_file), intent(in) :: the_file
!
      integer, intent(out) :: scalar1, scalar2
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      if (the_file%format_ .eq. 'unformatted') then
         read(the_file%unit, iostat=io_error, iomsg=io_msg) scalar1, scalar2
      else
         read(the_file%unit, *, iostat=io_error, iomsg=io_msg) scalar1, scalar2
      endif
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to read from file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine read_i_i_sequential_file
!
!
   subroutine read_i_1_sequential_file(the_file, array, n)
!!
!!    Sequential file read, integer 1 dimensional array
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
      if (the_file%format_ .eq. 'unformatted') then
         read(the_file%unit, iostat=io_error, iomsg=io_msg) array
      else
         read(the_file%unit, *, iostat=io_error, iomsg=io_msg) array
      endif
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to read from file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine read_i_1_sequential_file
!
!
   subroutine read_i_2_sequential_file(the_file, array, n)
!!
!!    Sequential file read, integer 2 dimensional array
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
   subroutine read_i_3_sequential_file(the_file, array, n)
!!
!!    Sequential file read, integer 3 dimensional array
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
   subroutine read_i_4_sequential_file(the_file, array, n)
!!
!!    Sequential file read, integer 4 dimensional array
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
   subroutine read_l_sequential_file(the_file, scalar)
!!
!!    Sequential file read, logical scalar
!!    Written by Rolf H. Myhre, September 2019
!!
      implicit none
!
      class(sequential_file), intent(in) :: the_file
!
      logical, intent(out) :: scalar
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      if (the_file%format_ .eq. 'unformatted') then
         read(the_file%unit, iostat=io_error, iomsg=io_msg) scalar
      else
         read(the_file%unit, *, iostat=io_error, iomsg=io_msg) scalar
      endif
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to read from file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine read_l_sequential_file
!
!
   subroutine read_l_1_sequential_file(the_file, array, n)
!!
!!    Sequential file read, logical 1 dimensional array
!!    Written by Rolf H. Myhre, September 2019
!!
      implicit none
!
      class(sequential_file), intent(in) :: the_file
!
      integer, intent(in)                 :: n
      logical, dimension(n), intent(out)  :: array
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      if (the_file%format_ .eq. 'unformatted') then
         read(the_file%unit, iostat=io_error, iomsg=io_msg) array
      else
         read(the_file%unit, *, iostat=io_error, iomsg=io_msg) array
      endif
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to read from file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine read_l_1_sequential_file
!
!
   subroutine read_char_sequential_file(the_file, string, format_string)
!!
!!    Sequential file read, character
!!    Written by Rolf H. Myhre, May 2019
!!    string:  character string to read in to
!!    format_string: optional format string
!!
      implicit none
!
      class(sequential_file), intent(in) :: the_file
!
      character(len=*), intent(out)          :: string
      character(len=*), intent(in), optional :: format_string
!
      integer              :: io_error
      character(len=100)   :: io_msg
      character(len=100)   :: fstring
!
      if(present(format_string)) then
         if (the_file%format_ .eq. 'unformatted') then
            call output%error_msg(trim(the_file%name_)//' is unformatted')
         endif
         fstring = trim(format_string)
      else
         fstring = '*'
      endif
!
      if (the_file%format_ .eq. 'unformatted') then
         read(the_file%unit, iostat=io_error, iomsg=io_msg) string
      else
         read(the_file%unit, trim(fstring), iostat=io_error, iomsg=io_msg) string
      endif
!
      if(io_error .ne. 0) then
         call output%error_msg('Failed to read from file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine read_char_sequential_file
!
!
end module sequential_file_class
