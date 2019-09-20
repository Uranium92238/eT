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
module abstract_other_file_class
!
!!
!!    Abstract other file class module
!!    Written by Rolf H. Myhre, September 2019
!!
!!    Abstract file class which direct and sequential file
!!    inherits from. Necessary in order to give these files 
!!    access to disk and output without circular use statements
!!
!
   use kinds          
!
   use global_out, only : output
!
   use abstract_file_class, only : abstract_file
!
   use disk_manager_class, only : disk
!
   type, abstract, extends(abstract_file) :: abstract_other_file 
!
   contains
!
      procedure :: close_  => close_abstract_other_file
      procedure :: delete_ => delete_abstract_other_file
      procedure :: copy_   => copy_abstract_other_file
!
   end type abstract_other_file
!
!
contains
!
!
   subroutine close_abstract_other_file(the_file, file_status)
!!
!!    Close file
!!    Written by Rolf Heilemann Myhre, May 2019
!!
      implicit none
!
      class(abstract_other_file)             :: the_file
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
         print *, 'Error: '//trim(the_file%name_)//' already closed'
      end if
!
      close(the_file%unit, iostat=io_error, iomsg=io_msg, status=trim(stat))
!
      if (io_error .ne. 0) then 
         print *, 'Error: could not close eT file '//trim(the_file%name_)
         print *, 'Error message: '//trim(io_msg)
         stop
      endif
!
      file_change = the_file%get_change()
      call disk%update(file_change, the_file%name_)
!
      the_file%is_open = .false.
      the_file%unit = -1
      the_file%action_ = 'unknown'
!
   end subroutine close_abstract_other_file
!
!
   subroutine delete_abstract_other_file(the_file)
!!
!!    Delete file
!!    Written by Rolf Heilemann Myhre, Aug 2019
!!
      implicit none
!
      class(abstract_other_file) :: the_file
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      if(the_file%is_open) then
!
         close(the_file%unit, iostat=io_error, iomsg=io_msg, status='delete')
!
         if (io_error .ne. 0) then 
            call output%error_msg('could not delete eT file '//trim(the_file%name_)//&
                                 &'. Error message: '//trim(io_msg))
         endif
!
      else
!
         open(newunit=the_file%unit, file=the_file%name_, access=the_file%access_, &
              action='write', iostat=io_error, iomsg=io_msg)
!
         if (io_error .ne. 0) then 
            call output%error_msg('could not open eT file '//trim(the_file%name_)//&
                                 &'. Error message: '//trim(io_msg))
         endif
!
         close(the_file%unit, iostat=io_error, iomsg=io_msg, status='delete')
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
   end subroutine delete_abstract_other_file
!
!
   subroutine copy_abstract_other_file(the_file, filename)
!!
!!    Copy abstract file
!!    Written by Alexander Paul and Rolf H. Myhre, September 2019
!!
!!    Very similar to abstract copy, but with access to output
!!
      implicit none
!
      class(abstract_other_file) :: the_file
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
      open(newunit=the_file%unit, file=the_file%name_, access='stream', &
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
         read(the_file%unit, iostat=io_error, iomsg=io_msg) !Read until end of file
!                                                                         
         if(io_error .gt. 0) then 
!
            call output%error_msg('Failed to read '//trim(the_file%name_)//' in copy.'//&
                                 &' io_msg: '//trim(io_msg))
!
         elseif(io_error .lt. 0) then !Reached end
!
            exit
!
         endif
!
         write(copy_unit, iostat=io_error, iomsg=io_msg) byte !Write whatever you just read
!
         if(io_error .ne. 0) then 
!
            call output%error_msg('Failed to write '//trim(filename)//' in copy.'//&
                                 &' io_msg: '//trim(io_msg))
!
         endif
!
      enddo
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
      close(the_file%unit, status='keep', iostat=io_error, iomsg=io_msg)
!
      if(io_error .ne. 0) then 
!
         call output%error_msg('Failed to close '//trim(the_file%name_)//' in copy.'//&
                              &' io_msg: '//trim(io_msg))
!
      endif
!
   end subroutine copy_abstract_other_file
!
!
end module abstract_other_file_class
