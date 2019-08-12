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
      procedure, public :: open_    => open_direct_file
      procedure, public :: close_   => close_direct_file
!
!     Copy a file
!
      procedure, public :: copy  => copy_direct_file
!
!     Write routines
!
      procedure, public :: write_dp_direct_file
      procedure, public :: write_dp_1_direct_file
      procedure, public :: write_dp_2_direct_file
      procedure, public :: write_dp_3_direct_file
      procedure, public :: write_dp_4_direct_file
      procedure, public :: write_i_direct_file
      procedure, public :: write_i_1_direct_file
      procedure, public :: write_i_2_direct_file
      procedure, public :: write_i_3_direct_file
      procedure, public :: write_i_4_direct_file
      generic           :: write_ => write_dp_direct_file,    &
                                     write_dp_1_direct_file,  &
                                     write_dp_2_direct_file,  &
                                     write_dp_3_direct_file,  &
                                     write_dp_4_direct_file,  &
                                     write_i_direct_file,     &
                                     write_i_1_direct_file,   &
                                     write_i_2_direct_file,   &
                                     write_i_3_direct_file,   &
                                     write_i_4_direct_file
!
!     Read routines
!
      procedure, public :: read_dp_direct_file
      procedure, public :: read_dp_1_direct_file
      procedure, public :: read_dp_2_direct_file
      procedure, public :: read_dp_3_direct_file
      procedure, public :: read_dp_4_direct_file
      procedure, public :: read_i_direct_file
      procedure, public :: read_i_1_direct_file
      procedure, public :: read_i_2_direct_file
      procedure, public :: read_i_3_direct_file
      procedure, public :: read_i_4_direct_file
      generic           :: read_ => read_dp_direct_file,    &
                                    read_dp_1_direct_file,  &
                                    read_dp_2_direct_file,  &
                                    read_dp_3_direct_file,  &
                                    read_dp_4_direct_file,  &
                                    read_i_direct_file,     &
                                    read_i_1_direct_file,   &
                                    read_i_2_direct_file,   &
                                    read_i_3_direct_file,   &
                                    read_i_4_direct_file
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
   function new_direct_file(name_, rec_dim, w_size) result(the_file)
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
      character(len=*), intent(in) :: name_
      integer, intent(in) :: rec_dim
      integer, intent(in), optional :: w_size
!
      if (present(w_size)) then
         if (w_size .gt. 0) then
            the_file%word_size = w_size
         else
            call output%error_msg("Word size less than zero for file "//trim(name_))
         endif
      else
         the_file%word_size = dp
      endif
!
      the_file%name_ = name_
!
      the_file%access_ = 'direct'
      the_file%format_ = 'unformatted'
!
      if (rec_dim .le. 0) then
         call output%error_msg("Record dimension less than zero for file "//name_)
      endif
!
      the_file%record_dim = rec_dim
      the_file%record_length = rec_dim*the_file%word_size
!
   end function new_direct_file
!
!
   subroutine open_direct_file(the_file, file_action)
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
      open(newunit=the_file%unit, file=the_file%name_, access=the_file%access_, &
           action=trim(act), recl=the_file%record_length, status='unknown', form=the_file%format_, &
           iostat=io_error, iomsg=io_msg)
!
      if (io_error .ne. 0) then 
         call output%error_msg('Error: could not open eT direct file '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
      the_file%opened = .true.
!
      call the_file%set_open_size()
!
   end subroutine open_direct_file
!
!
   subroutine copy_direct_file(the_file, filename)
!!
!!    Copy eT direct file
!!    Written by Alexander Paul and Rolf H. Myhre, August 2019
!!
      implicit none
!
      class(direct_file) :: the_file
!
      character(*), intent(in) :: filename
!
      type(direct_file) :: copy
!
      integer              :: io_error_1, io_error_2 
      character(len=100)   :: io_msg_1, io_msg_2
!
      character(len=the_file%record_length) :: line
!
      integer :: record, next_record
!
      logical :: file_was_open
!
      file_was_open = the_file%opened
      if (.not. file_was_open) call the_file%open_('read')
!
      copy = direct_file(filename, the_file%record_dim, the_file%word_size)
      call copy%open_('write')
!
      record = 1
!
!     Read record by record and write it to the new file
!
      do
!
         read(the_file%unit, rec=record, iostat=io_error_1, iomsg=io_msg_1) line
!
!        Direct access files don't recognize end of file
!        Therefore check number of the next record
!        if the current record is equal to the next then EOF -> exit
!
         inquire(the_file%unit, nextrec=next_record)
!
         if (next_record .eq. record) then
!
            exit 
!
         else if (io_error_1 .ne. 0) then
!
            call output%error_msg('Failed to read from file: '//trim(the_file%name_)//&
                                 &'. Error message: '//trim(io_msg_1))
!
         end if
!
         write(copy%unit, rec=record, iostat=io_error_2, iomsg=io_msg_2) line
!
         record = record + 1
      enddo
!
      if(io_error_2 .ne. 0) then
         call output%error_msg('Failed to write to file: '//trim(copy%name_)//&
                              &'. Error message: '//trim(io_msg_2))
      endif
!
      call copy%close_()
!
!     Close copied file again if it was closed in the beginning
!
      if (.not. file_was_open) call the_file%close_()
!
   end subroutine copy_direct_file
!
!
   subroutine close_direct_file(the_file, file_status)
!!
!!    Close the output file
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
         call output%error_msg('Error: could not close eT file '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
      the_file%opened = .false.
!
      file_change = the_file%get_change()
      call disk%update(file_change, the_file%name_)
!
   end subroutine close_direct_file
!
!
   subroutine write_dp_direct_file(the_file, scalar, record)
!!
!!    Direct file write_, real(dp) scalar
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(direct_file), intent(in) :: the_file
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
         call output%error_msg('Failed to write to file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine write_dp_direct_file
!
!
   subroutine write_dp_1_direct_file(the_file, array, record)
!!
!!    Direct file write_, real(dp) array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(direct_file), intent(in) :: the_file
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
         call output%error_msg('Failed to write to file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine write_dp_1_direct_file
!
!
   subroutine write_dp_2_direct_file(the_file, array, record)
!!
!!    Direct file write_, real(dp) array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(direct_file), intent(in)         :: the_file
      real(dp), dimension(:,:), intent(in)   :: array
      integer, intent(in)                    :: record
!
      call the_file%write_dp_1_direct_file(array, record)
!
   end subroutine write_dp_2_direct_file
!
!
   subroutine write_dp_3_direct_file(the_file, array, record)
!!
!!    Direct file write_, real(dp) array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(direct_file), intent(in)         :: the_file
      real(dp), dimension(:,:,:), intent(in) :: array
      integer, intent(in)                    :: record
!
      call the_file%write_dp_1_direct_file(array, record)
!
   end subroutine write_dp_3_direct_file
!
!
   subroutine write_dp_4_direct_file(the_file, array, record)
!!
!!    Direct file write_, real(dp) array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(direct_file), intent(in)            :: the_file
      real(dp), dimension(:,:,:,:), intent(in)  :: array
      integer, intent(in)                       :: record
!
      call the_file%write_dp_1_direct_file(array, record)
!
   end subroutine write_dp_4_direct_file
!
!
   subroutine write_i_direct_file(the_file, scalar, record)
!!
!!    Direct file write_, integer scalar
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(direct_file), intent(in) :: the_file
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
         call output%error_msg('Failed to write to file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine write_i_direct_file
!
!
   subroutine write_i_1_direct_file(the_file, array, record)
!!
!!    Direct file write_, integer array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(direct_file), intent(in) :: the_file
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
         call output%error_msg('Failed to write to file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine write_i_1_direct_file
!
!
   subroutine write_i_2_direct_file(the_file, array, record)
!!
!!    Direct file write_, integer array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(direct_file), intent(in)      :: the_file
      integer, dimension(:,:), intent(in) :: array
      integer, intent(in)                 :: record
!
      call the_file%write_i_1_direct_file(array, record)
!
   end subroutine write_i_2_direct_file
!
!
   subroutine write_i_3_direct_file(the_file, array, record)
!!
!!    Direct file write_, integer array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(direct_file), intent(in)         :: the_file
      integer, dimension(:,:,:), intent(in)  :: array
      integer, intent(in)                    :: record
!
      call the_file%write_i_1_direct_file(array, record)
!
   end subroutine write_i_3_direct_file
!
!
   subroutine write_i_4_direct_file(the_file, array, record)
!!
!!    Direct file write_, integer array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(direct_file), intent(in)            :: the_file
      integer, dimension(:,:,:,:), intent(in)   :: array
      integer, intent(in)                       :: record
!
      call the_file%write_i_1_direct_file(array, record)
!
   end subroutine write_i_4_direct_file
!
!
   subroutine read_dp_direct_file(the_file, scalar, record)
!!
!!    Direct file read_, real(dp) scalar
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(direct_file), intent(in) :: the_file
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
         call output%error_msg('Failed to read from file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine read_dp_direct_file
!
!
   subroutine read_dp_1_direct_file(the_file, array, record)
!!
!!    Direct file read_, real(dp) array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(direct_file), intent(in) :: the_file
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
         call output%error_msg('Failed to read from file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine read_dp_1_direct_file
!
!
   subroutine read_dp_2_direct_file(the_file, array, record)
!!
!!    Direct file read_, real(dp) array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(direct_file), intent(in)         :: the_file
      real(dp), dimension(:,:), intent(out)  :: array
      integer, intent(in)                    :: record
!
      call the_file%read_dp_1_direct_file(array,record)
!
   end subroutine read_dp_2_direct_file
!
!
   subroutine read_dp_3_direct_file(the_file, array, record)
!!
!!    Direct file read_, real(dp) array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(direct_file), intent(in)            :: the_file
      real(dp), dimension(:,:,:), intent(out)   :: array
      integer, intent(in)                       :: record
!
      call the_file%read_dp_1_direct_file(array,record)
!
   end subroutine read_dp_3_direct_file
!
!
   subroutine read_dp_4_direct_file(the_file, array, record)
!!
!!    Direct file read_, real(dp) array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(direct_file), intent(in)            :: the_file
      real(dp), dimension(:,:,:,:), intent(out) :: array
      integer, intent(in)                       :: record
!
      call the_file%read_dp_1_direct_file(array,record)
!
   end subroutine read_dp_4_direct_file
!
!
   subroutine read_i_direct_file(the_file, scalar, record)
!!
!!    Direct file read_, integer scalar
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(direct_file), intent(in) :: the_file
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
         call output%error_msg('Failed to read from file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine read_i_direct_file
!
!
   subroutine read_i_1_direct_file(the_file, array, record)
!!
!!    Direct file read_, integer array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(direct_file), intent(in) :: the_file
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
         call output%error_msg('Failed to read from file: '//trim(the_file%name_)//&
                              &'. Error message: '//trim(io_msg))
      endif
!
   end subroutine read_i_1_direct_file
!
!
   subroutine read_i_2_direct_file(the_file, array, record)
!!
!!    Direct file read_, integer array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(direct_file), intent(in)         :: the_file
      integer, dimension(:,:), intent(out)   :: array
      integer, intent(in)                    :: record
!
      call the_file%read_i_1_direct_file(array, record)
!
   end subroutine read_i_2_direct_file
!
!
   subroutine read_i_3_direct_file(the_file, array, record)
!!
!!    Direct file read_, integer array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(direct_file), intent(in)         :: the_file
      integer, dimension(:,:,:), intent(out) :: array
      integer, intent(in)                    :: record
!
      call the_file%read_i_1_direct_file(array, record)
!
   end subroutine read_i_3_direct_file
!
!
   subroutine read_i_4_direct_file(the_file, array, record)
!!
!!    Direct file read_, integer array
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
      class(direct_file), intent(in)            :: the_file
      integer, dimension(:,:,:,:), intent(out)  :: array
      integer, intent(in)                       :: record
!
      call the_file%read_i_1_direct_file(array, record)
!
   end subroutine read_i_4_direct_file
!
!
end module direct_file_class
