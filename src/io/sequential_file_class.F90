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
!
   type, extends(abstract_file) :: sequential_file
!
   contains
!
!     Open, skip and rewind
!
      procedure, public :: open_   => open_sequential_file
      procedure, public :: close_  => close_sequential_file
      procedure, public :: delete_ => delete_sequential_file
      procedure, public :: copy    => copy_sequential_file
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
      procedure, public :: write_c_2_sequential_file
      procedure, public :: write_c_3_sequential_file
      procedure, public :: write_c_4_sequential_file
!
      procedure, public :: write_i_sequential_file
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
                                     write_c_2_sequential_file,   &
                                     write_c_3_sequential_file,   &
                                     write_c_4_sequential_file,   &
                                     write_i_sequential_file,     &
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
      procedure, public :: read_blank => read_blank_sequential_file
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
      procedure, public :: read_i_1_sequential_file
      procedure, public :: read_i_2_sequential_file
      procedure, public :: read_i_3_sequential_file
      procedure, public :: read_i_4_sequential_file
!
      procedure, public :: read_l_sequential_file
      procedure, public :: read_l_1_sequential_file
!
      procedure, public :: read_char_sequential_file
!
      generic           :: read_ => read_r_sequential_file,     &
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
                                    read_i_1_sequential_file,   &
                                    read_i_2_sequential_file,   &
                                    read_i_3_sequential_file,   &
                                    read_i_4_sequential_file,   &
                                    read_l_sequential_file,     &
                                    read_l_1_sequential_file,   &
                                    read_char_sequential_file
!
      final :: destructor
!
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
!!    new sequential file
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
      the_file%unit_ = -1
!
   end function new_sequential_file
!
!
   subroutine open_sequential_file(the_file, file_action, position_)
!!
!!    Open
!!    Written by Rolf H. Myhre, May 2019
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
      open(newunit=the_file%unit_, file=the_file%name_, access=the_file%access_, &
           action=the_file%action_, status='unknown', form=the_file%format_, &
           position=pos, iostat=io_error, iomsg=io_msg)
!
      if (io_error .ne. 0) then 
         call output%error_msg('could not open eT sequential file '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
      the_file%is_open = .true.
!
   end subroutine open_sequential_file
!
!
   subroutine close_sequential_file(the_file, new_destiny)
!!
!!    Close
!!    Written by Rolf H. Myhre Nov. 2019
!!
!!    Destiny: Optional character string that decides whether to keep
!!             or delete. Default: keep
!!
      implicit none
!
      class(sequential_file), intent(inout) :: the_file
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
                               &error message is (a0)', &
                               chars = [character(len=255)::the_file%name_, io_message], &
                               ints = [io_status])
      end if
!
!
      the_file%is_open = .false.
      the_file%unit_ = 1
      the_file%action_ = 'unknown'
!
   end subroutine close_sequential_file
!
!
   subroutine delete_sequential_file(the_file)
!!
!!    Delete
!!    Written by Rolf H. Myhre, Aug 2019
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
   end subroutine delete_sequential_file
!
!
   subroutine copy_sequential_file(the_file, filename)
!!
!!    Copy 
!!    Written by Alexander C. Paul and Rolf H. Myhre, September 2019
!!
!!    Very similar to abstract copy, but with access to output
!!
      implicit none
!
      class(sequential_file) :: the_file
!
      character(*), intent(in) :: filename
!
      integer              :: copy_unit
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
!     Character to hold a byte
      character :: byte
!
!     Check that file is closed
      if(the_file%is_open) then 
!
         call output%error_msg(the_file%name_//' is not closed in copy.')
!
      endif
!
!     Open the file with stream unformatted access
      open(newunit=the_file%unit_, file=the_file%name_, access='stream', &
           form='unformatted', action='read', status='old', & 
           iostat=io_error, iomsg=io_msg)
!
      if(io_error .ne. 0) then 
!
         call output%error_msg('Failed to open '//trim(the_file%name_)//' in copy.'//&
                              &' io_msg: '//trim(io_msg))
!
      endif
!
!     Open a new file
      open(newunit=copy_unit, file=trim(filename), access='stream', &
           form='unformatted', action='write', status='new', &
           iostat=io_error, iomsg=io_msg)
!
      if(io_error .ne. 0) then 
!
         call output%error_msg('Failed to open '//trim(filename)//' in copy.'//&
                              &' io_msg: '//trim(io_msg))
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
      close(copy_unit, status='keep', iostat=io_error, iomsg=io_msg)
!
      if(io_error .ne. 0) then 
!
         call output%error_msg('Failed to close '//trim(filename)//' in copy.'//&
                              &' io_msg: '//trim(io_msg))
!
      endif
!
      close(the_file%unit_, status='keep', iostat=io_error, iomsg=io_msg)
!
      if(io_error .ne. 0) then 
!
         call output%error_msg('Failed to close '//trim(the_file%name_)//' in copy.'//&
                              &' io_msg: '//trim(io_msg))
!
      endif
!
   end subroutine copy_sequential_file
!
!
   subroutine rewind_sequential_file(the_file)
!!
!!    Rewind
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(sequential_file)                 :: the_file
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
!     Check if the file is open
      if (.not. the_file%is_open) then
!
         call output%error_msg(trim('in rewind '//the_file%name_)//' is not open')
!
      endif
!
      rewind(the_file%unit_, iostat=io_error, iomsg=io_msg)
!
      if (io_error .ne. 0) then 
         call output%error_msg('Could not rewind eT file '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine rewind_sequential_file
!
!
   subroutine skip_sequential_file(the_file, jump)
!!
!!    Skip
!!    Written by Rolf H. Myhre, May 2019
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
         read(the_file%unit_, iostat=io_error, iomsg=io_msg) 
!
         if (io_error .ne. 0) then 
            call output%error_msg('Could not skip eT sequential file '//trim(the_file%name_)//&
                                 &'. Error message: '//trim(io_msg))
         endif
!
      enddo
!
   end subroutine skip_sequential_file
!
!
   subroutine write_blank_sequential_file(the_file)
!!
!!    write blank
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
         write(the_file%unit_, iostat=io_error, iomsg=io_msg) 
      else
         write(the_file%unit_, *, iostat=io_error, iomsg=io_msg) 
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
         write(the_file%unit_, iostat=io_error, iomsg=io_msg) scalar
      else
         write(the_file%unit_, *, iostat=io_error, iomsg=io_msg) scalar
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
         write(the_file%unit_, iostat=io_error, iomsg=io_msg) array
      else
         write(the_file%unit_, *, iostat=io_error, iomsg=io_msg) array
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
         write(the_file%unit_, iostat=io_error, iomsg=io_msg) scalar
      else
         write(the_file%unit_, *, iostat=io_error, iomsg=io_msg) scalar
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
         write(the_file%unit_, iostat=io_error, iomsg=io_msg) array
      else
         write(the_file%unit_, *, iostat=io_error, iomsg=io_msg) array
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
         write(the_file%unit_, iostat=io_error, iomsg=io_msg) scalar
      else
         write(the_file%unit_, *, iostat=io_error, iomsg=io_msg) scalar
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
         write(the_file%unit_, iostat=io_error, iomsg=io_msg) array
      else
         write(the_file%unit_, *, iostat=io_error, iomsg=io_msg) array
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
         write(the_file%unit_, iostat=io_error, iomsg=io_msg) scalar
      else
         write(the_file%unit_, *, iostat=io_error, iomsg=io_msg) scalar
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
         write(the_file%unit_, iostat=io_error, iomsg=io_msg) array
      else
         write(the_file%unit_, *, iostat=io_error, iomsg=io_msg) array
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
!
      if (the_file%format_ .eq. 'unformatted') then
!
         if (present(format_string)) call output%error_msg(trim(the_file%name_)//' is unformatted')
         write(the_file%unit_, iostat=io_error, iomsg=io_msg) string
!
      else
!
         if (present(format_string) .and. trim(format_string) .ne. '*') then
            write(the_file%unit_, trim(format_string), iostat=io_error, iomsg=io_msg) string
         else
            write(the_file%unit_, *, iostat=io_error, iomsg=io_msg) string
         endif
!
      endif
!
      if (io_error .ne. 0) then
         call output%error_msg('Failed to write to file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine write_char_sequential_file
!
!
   subroutine read_blank_sequential_file(the_file, io_stat)
!!
!!    Sequential file read, blank line
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(sequential_file), intent(in)  :: the_file
!
      integer, optional :: io_stat 
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      if (the_file%format_ .eq. 'unformatted') then
         read(the_file%unit_, iostat=io_error, iomsg=io_msg)
      else
         read(the_file%unit_, *, iostat=io_error, iomsg=io_msg)
      endif
!
      if (present(io_stat) .and. io_error .le. 0) then 
!
         io_stat = io_error
!  
      elseif (io_error .ne. 0) then
!
         call output%error_msg('Failed to read from file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
!
      endif
!
   end subroutine read_blank_sequential_file
!
!
   subroutine read_r_sequential_file(the_file, scalar, io_stat)
!!
!!    Sequential file read, real(dp) scalar
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(sequential_file), intent(in)  :: the_file
!
      real(dp), intent(out)          :: scalar
      integer, intent(out), optional :: io_stat
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      if (the_file%format_ .eq. 'unformatted') then
         read(the_file%unit_, iostat=io_error, iomsg=io_msg) scalar
      else
         read(the_file%unit_, *, iostat=io_error, iomsg=io_msg) scalar
      endif
!
      if(present(io_stat) .and. io_error .le. 0) then
!
         io_stat = io_error
!
      elseif(io_error .ne. 0) then
!
         call output%error_msg('Failed to read from file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine read_r_sequential_file
!
!
   subroutine read_r_1_sequential_file(the_file, array, n, io_stat)
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
      integer, intent(out), optional      :: io_stat
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      if (the_file%format_ .eq. 'unformatted') then
         read(the_file%unit_, iostat=io_error, iomsg=io_msg) array
      else
         read(the_file%unit_, *, iostat=io_error, iomsg=io_msg) array
      endif
!
      if(present(io_stat) .and. io_error .le. 0) then
!
         io_stat = io_error
!
      elseif(io_error .ne. 0) then
!
         call output%error_msg('Failed to read from file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine read_r_1_sequential_file
!
!
   subroutine read_r_2_sequential_file(the_file, array, n, io_stat)
!!
!!    Sequential file read, real(dp) 2 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(sequential_file), intent(in)     :: the_file
      integer, intent(in)                    :: n
      real(dp), dimension(:,:), intent(out)  :: array
      integer, intent(out), optional         :: io_stat
!
      call the_file%read_r_1_sequential_file(array, n, io_stat)
!
   end subroutine read_r_2_sequential_file
!
!
   subroutine read_r_3_sequential_file(the_file, array, n, io_stat)
!!
!!    Sequential file read, real(dp) 3 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(sequential_file), intent(in)        :: the_file
      integer, intent(in)                       :: n
      real(dp), dimension(:,:,:), intent(out)   :: array
      integer, intent(out), optional            :: io_stat
!
      call the_file%read_r_1_sequential_file(array, n, io_stat)
!
   end subroutine read_r_3_sequential_file
!
!
   subroutine read_r_4_sequential_file(the_file, array, n, io_stat)
!!
!!    Sequential file read, real(dp) 4 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(sequential_file), intent(in)        :: the_file
      integer, intent(in)                       :: n
      real(dp), dimension(:,:,:,:), intent(out) :: array
      integer, intent(out), optional            :: io_stat
!
      call the_file%read_r_1_sequential_file(array, n, io_stat)
!
   end subroutine read_r_4_sequential_file
!
!
   subroutine read_c_sequential_file(the_file, scalar, io_stat)
!!
!!    Sequential file read, complex(dp) scalar
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(sequential_file), intent(in)  :: the_file
!
      complex(dp), intent(out)       :: scalar
      integer, intent(out), optional :: io_stat
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      if (the_file%format_ .eq. 'unformatted') then
         read(the_file%unit_, iostat=io_error, iomsg=io_msg) scalar
      else
         read(the_file%unit_, *, iostat=io_error, iomsg=io_msg) scalar
      endif
!
      if(present(io_stat) .and. io_error .le. 0) then
!
         io_stat = io_error
!
      elseif(io_error .ne. 0) then
!
         call output%error_msg('Failed to read from file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine read_c_sequential_file
!
!
   subroutine read_c_1_sequential_file(the_file, array, n, io_stat)
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
      integer, intent(out), optional         :: io_stat
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      if (the_file%format_ .eq. 'unformatted') then
         read(the_file%unit_, iostat=io_error, iomsg=io_msg) array
      else
         read(the_file%unit_, *, iostat=io_error, iomsg=io_msg) array
      endif
!
      if(present(io_stat) .and. io_error .le. 0) then
!
         io_stat = io_error
!
      elseif(io_error .ne. 0) then
!
         call output%error_msg('Failed to read from file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine read_c_1_sequential_file
!
!
   subroutine read_c_2_sequential_file(the_file, array, n, io_stat)
!!
!!    Sequential file read, complex(dp) 2 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(sequential_file), intent(in)       :: the_file
      integer, intent(in)                      :: n
      complex(dp), dimension(:,:), intent(out) :: array
      integer, intent(out), optional           :: io_stat
!
      call the_file%read_c_1_sequential_file(array, n, io_stat)
!
   end subroutine read_c_2_sequential_file
!
!
   subroutine read_c_3_sequential_file(the_file, array, n, io_stat)
!!
!!    Sequential file read, complex(dp) 3 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(sequential_file), intent(in)         :: the_file
      integer, intent(in)                        :: n
      complex(dp), dimension(:,:,:), intent(out) :: array
      integer, intent(out), optional             :: io_stat
!
      call the_file%read_c_1_sequential_file(array, n, io_stat)
!
   end subroutine read_c_3_sequential_file
!
!
   subroutine read_c_4_sequential_file(the_file, array, n, io_stat)
!!
!!    Sequential file read, complex(dp) 4 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(sequential_file), intent(in)           :: the_file
      integer, intent(in)                          :: n
      complex(dp), dimension(:,:,:,:), intent(out) :: array
      integer, intent(out), optional               :: io_stat
!
      call the_file%read_c_1_sequential_file(array, n, io_stat)
!
   end subroutine read_c_4_sequential_file
!
!
   subroutine read_i_sequential_file(the_file, scalar, io_stat)
!!
!!    Sequential file read, integer scalar
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(sequential_file), intent(in) :: the_file
!
      integer, intent(out) :: scalar
      integer, intent(out), optional :: io_stat
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      if (the_file%format_ .eq. 'unformatted') then
         read(the_file%unit_, iostat=io_error, iomsg=io_msg) scalar
      else
         read(the_file%unit_, *, iostat=io_error, iomsg=io_msg) scalar
      endif
!
      if(present(io_stat) .and. io_error .le. 0) then
!
         io_stat = io_error
!
      elseif(io_error .ne. 0) then
!
         call output%error_msg('Failed to read from file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine read_i_sequential_file
!
!
   subroutine read_i_1_sequential_file(the_file, array, n, io_stat)
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
      integer, intent(out), optional      :: io_stat
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      if (the_file%format_ .eq. 'unformatted') then
         read(the_file%unit_, iostat=io_error, iomsg=io_msg) array
      else
         read(the_file%unit_, *, iostat=io_error, iomsg=io_msg) array
      endif
!
      if(present(io_stat) .and. io_error .le. 0) then
!
         io_stat = io_error
!
      elseif(io_error .ne. 0) then
!
         call output%error_msg('Failed to read from file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine read_i_1_sequential_file
!
!
   subroutine read_i_2_sequential_file(the_file, array, n, io_stat)
!!
!!    Sequential file read, integer 2 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(sequential_file), intent(in)     :: the_file
      integer, intent(in)                    :: n
      integer, dimension(:,:), intent(out)   :: array
      integer, intent(out), optional         :: io_stat
!
      call the_file%read_i_1_sequential_file(array, n, io_stat)
!
   end subroutine read_i_2_sequential_file
!
!
   subroutine read_i_3_sequential_file(the_file, array, n, io_stat)
!!
!!    Sequential file read, integer 3 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(sequential_file), intent(in)     :: the_file
      integer, intent(in)                    :: n
      integer, dimension(:,:,:), intent(out) :: array
      integer, intent(out), optional         :: io_stat
!
      call the_file%read_i_1_sequential_file(array, n, io_stat)
!
   end subroutine read_i_3_sequential_file
!
!
   subroutine read_i_4_sequential_file(the_file, array, n, io_stat)
!!
!!    Sequential file read, integer 4 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(sequential_file), intent(in)        :: the_file
      integer, intent(in)                       :: n
      integer, dimension(:,:,:,:), intent(out)  :: array
      integer, intent(out), optional            :: io_stat
!
      call the_file%read_i_1_sequential_file(array, n, io_stat)
!
   end subroutine read_i_4_sequential_file
!
!
   subroutine read_l_sequential_file(the_file, scalar, io_stat)
!!
!!    Sequential file read, logical scalar
!!    Written by Rolf H. Myhre, September 2019
!!
!!    scalar  : scalar logical to read into
!!    io_stat : optional integer set to io_error if it is less than or equal to 0
!!
      implicit none
!
      class(sequential_file), intent(in) :: the_file
!
      logical, intent(out)           :: scalar
      integer, intent(out), optional :: io_stat
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      if (the_file%format_ .eq. 'unformatted') then
         read(the_file%unit_, iostat=io_error, iomsg=io_msg) scalar
      else
         read(the_file%unit_, *, iostat=io_error, iomsg=io_msg) scalar
      endif
!
      if(present(io_stat) .and. io_error .le. 0) then
!
         io_stat = io_error
!
      elseif(io_error .ne. 0) then
!
         call output%error_msg('Failed to read from file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine read_l_sequential_file
!
!
   subroutine read_l_1_sequential_file(the_file, array, n, io_stat)
!!
!!    Sequential file read, logical 1 dimensional array
!!    Written by Rolf H. Myhre, September 2019
!!
!!    array   : array of logicals to read into
!!    n       : number of elements to read
!!    io_stat : optional integer set to io_error if it is less than or equal to 0
!!
      implicit none
!
      class(sequential_file), intent(in) :: the_file
!
      integer, intent(in)                 :: n
      logical, dimension(n), intent(out)  :: array
      integer, intent(out), optional      :: io_stat
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      if (the_file%format_ .eq. 'unformatted') then
         read(the_file%unit_, iostat=io_error, iomsg=io_msg) array
      else
         read(the_file%unit_, *, iostat=io_error, iomsg=io_msg) array
      endif
!
      if(present(io_stat) .and. io_error .le. 0) then
!
         io_stat = io_error
!
      elseif(io_error .ne. 0) then
!
         call output%error_msg('Failed to read from file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine read_l_1_sequential_file
!
!
   subroutine read_char_sequential_file(the_file, string, format_string, io_stat)
!!
!!    Sequential file read, character
!!    Written by Rolf H. Myhre, May 2019
!!    string:  character string to read in to
!!    format_string: optional format string
!!
!!    Modified by Andreas Skeidsvoll, Dec 2019
!!    Added optional argument
!!    io_stat: optional integer set to io_error if it is less than or equal to 0
!!
      implicit none
!
      class(sequential_file), intent(in) :: the_file
!
      character(len=*), intent(out)          :: string
      character(len=*), intent(in), optional :: format_string
      integer, intent(out), optional         :: io_stat
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      if (the_file%format_ .eq. 'unformatted') then
!
         if (present(format_string)) call output%error_msg(trim(the_file%name_)//' is unformatted')
         read(the_file%unit_, iostat=io_error, iomsg=io_msg) string
!
      else
!
         if (present(format_string) .and. trim(format_string) .ne. '*') then
            read(the_file%unit_, trim(format_string), iostat=io_error, iomsg=io_msg) string
         else
            read(the_file%unit_, *, iostat=io_error, iomsg=io_msg) string
         endif
!
      endif
!
      if (present(io_stat) .and. io_error .le. 0) then
!
         io_stat = io_error
!
      elseif (io_error .ne. 0) then
!
         call output%error_msg('Failed to read from file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
!
      endif
!
   end subroutine read_char_sequential_file
!
!
   subroutine destructor(the_file)
!!
!!    Destructor 
!!    Written by Rolf H. Myhre, Sep 2019 
!!
      implicit none 
!
      type(sequential_file) :: the_file
!
      if (the_file%is_open) then
         call output%error_msg('Destructor for file (a0) called &
                               &while the file is still open', chars=[the_file%name_])
      endif
!
   end subroutine destructor
!
!
end module sequential_file_class
