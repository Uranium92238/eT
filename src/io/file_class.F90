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
!
   type :: file
!
!     Filename
!
      character(len=255) :: name = 'no_name'
!
!     Unit identifier
!
      integer :: unit = -1
!
!     File size (in bytes)
!
      integer, private :: file_size = -1
!
!     Logical for whether the file is currently opened or not
!
      logical :: opened = .false.
!
      character(len=40) :: access = 'unknown'
      character(len=40) :: format = 'unknown'
!
      integer :: record_length = 0
!
   contains
!
      procedure :: init                   => init_file
!
      procedure :: prepare_to_read_line   => prepare_to_read_line_file
!
      procedure :: error_msg              => error_msg_file
      procedure :: warning_msg            => warning_msg_file
!
      procedure :: determine_file_size    => determine_file_size_file
!
      procedure :: get_file_size          => get_file_size_file
      procedure :: file_exists            => file_exists_file
!
   end type file
!
!  The 'global' eT files 
!
   type(file) :: output
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
   subroutine error_msg_file(out_file, error_specs, error_int)
!!
!!    Error message
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(file) :: out_file
!
      character(len=*) :: error_specs
!
      integer, optional :: error_int 
!
      character(len=40) :: error_int_char = ' '
!
      if (present(error_int)) then
!
         write(error_int_char, '(i12)') error_int   
!
         write(out_file%unit, '(a)') 'Error: ' // trim(error_specs) // ' ' // error_int_char
!
      else
!
         write(out_file%unit, '(a)') 'Error: ' // trim(error_specs)
!
      endif
!
      stop
!
   end subroutine error_msg_file
!
!
   subroutine warning_msg_file(out_file, warning_specs)
!!
!!    Warning message
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(file) :: out_file
!
      character(len=*) :: warning_specs
!
      write(out_file%unit, '(a)') 'Warning: ' // trim(warning_specs)
!
   end subroutine warning_msg_file
!
!
   subroutine determine_file_size_file(the_file)
!!
!!    Determine file size
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2018
!!    Moved to file by Rolf H. Myhre Nov. 2018
!!
!!    The disk manager handles files. This routine is called by it
!!    and should never be called by the user (because it can lead to
!!    errors in the disk space estimates).
!!
      implicit none
!
      class(file) :: the_file
!
!     Inquire about the file size
!
      inquire(file=the_file%name, size=the_file%file_size)
!
!     Check whether the file size could be calculated
!
      if (the_file%file_size .eq. -1) then
!
         call output%error_msg('could not calculate file size of the file ' // trim(the_file%name))
!
      endif
!
   end subroutine determine_file_size_file
!
!
   function get_file_size_file(the_file)
!!
!!    Return private variable file_size
!!    Written by Rolf H. Myhre, 2018
!!    
!
      implicit none
!  
      class(file), intent(in) :: the_file
!
      integer :: get_file_size_file
!
      get_file_size_file = the_file%file_size
!
   end function get_file_size_file
!
!  
   function file_exists_file(the_file)
!!
!!    File exists
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!    Moved to file by Rolf H. Myhre Nov. 2018
!!    
      implicit none
!  
      class(file), intent(in) :: the_file
!
      logical :: file_exists_file
!
      inquire(file=the_file%name, exist=file_exists_file)
!
   end function file_exists_file
!
!  
end module file_class
