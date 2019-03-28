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
   use abstract_file_class      
   use output_file_class      
!
   type, extends(abstract_file) :: file
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
!  The 'global' eT files 
!
   type(file) :: input
   type(file) :: timing 
!
contains
!
!
   subroutine init_file(the_file, name, access, format, record_length)
!!
!!    Initialize file
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2018
!!
!!    Initializes a file object
!!
!!       - access ('direct' or 'sequential')
!!       - format ('formatted' or unformatted)
!!       - record_length, optional argument which must be provided if
!!         the_file is direct access
!!
!!    File is closed, to open it use disk_manager
!!
      implicit none
!
      class(file) :: the_file
!
      character(len=*) :: name
      character(len=*) :: access
      character(len=*) :: format
!
      integer, optional :: record_length
!
!     Sanity checks
!
      the_file%name = name
!
      if (.not. present(record_length)) then
!
        if (access == 'direct') then
!
            call output%error_msg('for direct access files a record length must be specified.')
!
         endif
!
      elseif (access .ne. 'direct' .and. access .ne. 'sequential') then
!
         call output%error_msg('illegal access type specified for file: ' // name // access)
!
      elseif (format .ne. 'unformatted' .and. format .ne. 'formatted') then
!
         call output%error_msg('illegal format specified for file: ' // name)
!
      endif
!
      the_file%opened = .false.
!
      the_file%access = access
      the_file%format = format
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
      if (.not. the_file%opened) then
!
         call output%error_msg('attempted to read unopened file:' // trim(the_file%name))
!
      elseif (the_file%access == 'direct') then
!
         call output%warning_msg('no need to prepare to read line for a direct access file.')
         return
!
      endif
!
      rewind(the_file%unit)
!
      if (the_file%format == 'unformatted') then
!
         do i = 1, (line - 1)
!
            read(the_file%unit)
!
         enddo
!
      elseif (the_file%format == 'formatted') then
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
