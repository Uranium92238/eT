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
      character(len=40) :: name = 'no_name'
!
!     Unit identifier 
!
      integer(i15) :: unit = -1     
!
!     File size (in bytes)
! 
      integer(i15) :: size = -1  
!
!     Logical for whether the file is currently opened or not 
!
      logical :: opened = .false.
!
      character(len=40) :: access = 'unknown'
      character(len=40) :: format = 'unknown'
!
      integer(i15) :: record_length = 0
!
   contains 
!
      procedure :: init => init_file
!
      procedure :: prepare_to_read_line => prepare_to_read_line_file
!
   end type file                                                                             
!
      type(file) :: output
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
      integer(i15), optional :: record_length 
!
!     Sanity checks
!
      if (.not. present(record_length)) then
        if (access == 'direct') then
!
            write(output%unit,'(/t3,a)') 'Error: for direct access files a record length must be specified.'
            stop
         endif
!
      elseif (access .ne. 'direct' .and. access .ne. 'sequential') then
!
         write(output%unit,'(/t3,a)') 'Error: illegal access type specified for file: ', name
         stop
!
      elseif (format .ne. 'unformatted' .and. format .ne. 'formatted') then
!
         write(output%unit,'(/t3,a)') 'Error: illegal format specified for file: ', name
         stop
!
      endif
!
      the_file%opened = .false.
!
      the_file%name = name
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
      integer(i15) :: line
!
      integer(i15) :: i = 0
!
!     Sanity checks
!
      if (.not. the_file%opened) then
!
         write(output%unit,'(/t3,a)') 'Error: attempted to read unopened file.'
         stop
!
      elseif (the_file%access == 'direct') then
!
         write(output%unit,'(/t3,a)') 'Warning: no need to prepare to read line for a direct access file.'
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
