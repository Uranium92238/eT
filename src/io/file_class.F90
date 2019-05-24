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
module file_class
!
!!
!!    File class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2018
!!
!!
!
   use kinds    
   use output_file_class, only : output
   use abstract_file_class      
!
   type, extends(abstract_file) :: file
!
      integer :: record_length = 0
!
!
   contains
!
      procedure :: init                   => init_file
!
      procedure :: prepare_to_read_line   => prepare_to_read_line_file
!
   end type file
!
!
contains
!
!
   subroutine init_file(the_file, file_name, file_access, file_format, record_length)
!!
!!    Initialize file
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2018
!!
!!    Initializes a file object
!!
!!       - file_access ('direct' or 'sequential')
!!       - file_format ('formatted' or unformatted)
!!       - record_length, optional argument which must be provided if
!!         the_file is direct access
!!
!!    File is closed, to open it use disk_manager
!!
      implicit none
!
      class(file) :: the_file
!
      character(len=*) :: file_name
      character(len=*) :: file_access
      character(len=*) :: file_format
!
      integer, optional :: record_length
!
!     Sanity checks
!
      the_file%file_name = file_name
!
      if (.not. present(record_length)) then
!
        if (file_access == 'direct') then
!
            call output%error_msg('for direct access files a record length must be specified.')
!
         endif
!
      elseif (file_access .ne. 'direct' .and. file_access .ne. 'sequential') then
!
         call output%error_msg('illegal access type specified for file: ' // file_name // file_access)
!
      elseif (file_format .ne. 'unformatted' .and. file_format .ne. 'formatted') then
!
         call output%error_msg('illegal format specified for file: ' // file_name)
!
      endif
!
      the_file%file_opened = .false.
!
      the_file%file_access = file_access
      the_file%file_format = file_format
!
      if (present(record_length)) then
!
         the_file%record_length = record_length
!
      else
!
        the_file%record_length = 0
!
      endif
!
   end subroutine init_file
!
!
   subroutine prepare_to_read_line_file(the_file, line)
!!
!!    Prepare to read line
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2018
!!
!!    Places cursor for reading a sequential file at the position line
!!    by rewinding and doing (line - 1) emply reads
!!
!!    Requires file to be initialized and open. Files are opened by the disk_manager
!!
      implicit none
!
      class(file) :: the_file
!
      integer :: line
!
      integer :: i = 0
!
!     Sanity checks
!
      if (.not. the_file%file_opened) then
!
         call output%error_msg('attempted to read unopened file:' // trim(the_file%file_name))
!
      elseif (the_file%file_access == 'direct') then
!
         call output%warning_msg('no need to prepare to read line for a direct access file.')
         return
!
      endif
!
      rewind(the_file%unit)
!
      if (the_file%file_format == 'unformatted') then
!
         do i = 1, (line - 1)
!
            read(the_file%unit)
!
         enddo
!
      elseif (the_file%file_format == 'formatted') then
!
         do i = 1, (line - 1)
!
            read(the_file%unit, *)
!
         enddo
!
      endif
!
   end subroutine prepare_to_read_line_file
!
!  
end module file_class
