module disk_manager_class
!
!!
!!    Disk manager class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2018
!!
!
   use kinds
   use file_class
   use io_utilities
!
   type :: disk_manager
!
!     The total amount of disk space specified by user (standard: 30 GB)
!
      integer(i15), private :: total = 3000000000000
!
!     The amount of disk space currently available, based on the files
!     currently stored on file
!
      integer(i15), private :: available = 3000000000000
!
   contains
!
!     Initialization routine (used if user specifies a disk space different from standard)
!
      procedure :: prepare => prepare_disk_manager
!
!     Account for disk space already used (the size of the calculation folder)
!
      procedure :: subtract_folder_size => subtract_folder_size_disk_manager
!
!     Routine to open and close files
!
      procedure :: open_file_sequential   => open_file_sequential_disk_manager
      procedure :: open_file_direct       => open_file_direct_disk_manager
      procedure :: open_file              => open_file_disk_manager
!
      procedure :: close_file => close_file_disk_manager
!
!     Routine to determine file size
!
      procedure, private :: determine_file_size => determine_file_size_disk_manager
!
      procedure :: read_settings  => read_settings_disk_manager
      procedure :: print_settings => print_settings_disk_manager
!
      procedure :: delete => delete_disk_manager
!
   end type disk_manager
!
   type(disk_manager) :: disk
!
contains
!
!
   subroutine subtract_folder_size_disk_manager(disk)
!!
!!    Subtract folder size
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2018
!!
      implicit none
!
      class(disk_manager) :: disk
!
      type(file) :: scratch_size
!
!     Strings for reading the size of the directory from file
!
      character(len=40) :: scratch_size_entry
      character(len=40) :: scratch_size_entry_only_size
!
      integer(i15) :: counter
!
      integer(i15) :: size_of_directory ! in MB
!
!     For calculations with restart, it is necessary to account for the files
!     already present in the folder. We calculate this and subtract from the available.
!
      call system('du -m . >> scratch_size') ! Ask for number of megabytes
!
!     Open file to read the size
!
      call scratch_size%init('scratch_size', 'sequential', 'formatted')
      call disk%open_file(output, 'readwrite')
!
      read(scratch_size%unit, *) scratch_size_entry
!
      rewind(scratch_size%unit)
!
      counter = 1
      do while (scratch_size_entry(counter:counter) /= '.')
!
         counter = counter + 1
!
      enddo
!
      scratch_size_entry_only_size = scratch_size_entry(1:counter-1)
!
!     Now write the size to the file, rewind, and read into integer
!
      write(scratch_size%unit, *) scratch_size_entry_only_size
!
      rewind(scratch_size%unit)
!
      read(scratch_size%unit, *) size_of_directory
!
      write(output%unit,'(/t3,a36,i15)') 'Size of calculation directory (MB): ',size_of_directory
!
      size_of_directory = size_of_directory*1000 ! Megabytes to bytes
!
!     Set the initial available space
!
      disk%available = disk%total - size_of_directory
!
      call disk%close_file(scratch_size)
!
   end subroutine subtract_folder_size_disk_manager
!
!
   subroutine prepare_disk_manager(disk)
!!
!!    Init (disk manager)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2018
!!
!!    Initializes the disk manager object by setting the
!!    total and initial available disk. This is only called
!!    if the user specifies a total disk different from the standard.
!!
      implicit none
!
      class(disk_manager) :: disk
!
      integer(i15) :: total ! in GB
!
      if (requested_section('disk')) then
!
         call disk%read_settings()
!
      else
!
!        Set default
!
         disk%total = 30*1000000000
!
      endif
!
      disk%available = disk%total
!
      call disk%print_settings()
!
   end subroutine prepare_disk_manager
!
!
 subroutine open_file_disk_manager(disk, the_file, permissions, pos)
!!
!!    Open file
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2018
!!
!!    Wrapper for opening files.
!!
!!    The routine takes the following arguments:
!!
!!       - the_file (an object of type "file"). It is assumed that a file name
!!         has been set: i.e., the_file%name = 'filename'.
!!       - permissions ('read', 'write', 'readwrite')
!!       - pos ('rewind', 'append'). Optional argument specified for overwriting or appending
!!         sequential file. Default is writing to current position where ever that might be.
!!
      implicit none
!
      class(disk_manager) :: disk
!
      class(file) :: the_file ! the file
!
      character(len=*) :: permissions
      character(len=*), optional :: pos
!
      integer(i15) :: io_error = -1
!
!     Sanity checks
!
      if ( present(pos)) then
         if (the_file%access  == 'direct') then
!
            write(output%unit,'(/t3,a)') 'Warning: position specifier is disregarded for direct access file.'
            stop
!
         endif
      elseif (the_file%access .ne. 'direct' .and. the_file%access .ne. 'sequential' ) then
!
         write(output%unit,'(/t3,a)') 'Error: illegal access type for file: ', the_file%name, the_file%access
         stop
!
      endif
!
      if (the_file%access  == 'direct') then
!
         call disk%open_file_direct(the_file, permissions)
!
      elseif (the_file%access  == 'sequential') then
!
         call disk%open_file_sequential(the_file, permissions, pos)
!
      endif
!
   end subroutine open_file_disk_manager
!
!
   subroutine open_file_sequential_disk_manager(disk, the_file, permissions, pos)
!!
!!    Open sequential file
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2018
!!
!!    Opens a sequential access file object
!!
!!    The routine takes the following arguments:
!!
!!       - the_file (an object of type "file"). It is assumed that a file name
!!         has been set: i.e., the_file%name = 'filename'.
!!       - permissions ('read', 'write', 'readwrite')
!!       - position ('rewind', 'append'). Optional argument specified for overwriting or appending.
!!                                        Default is writing to current position whatever that might be.
!!
      implicit none
!
      class(disk_manager) :: disk
!
      class(file) :: the_file ! the file
!
      character(len=*) :: permissions
      character(len=*), optional :: pos
!
      integer(i15) :: io_error = -1
!
!     Sanity checks
!
      if (the_file%name == 'no_name') then
!
         write(output%unit,'(/t3,a)') 'Error: to open a file, you must set the name of the file.'
         stop
!
      elseif (the_file%format == 'unknown') then
!
         write(output%unit,'(/t3,a)') 'Error: to open a file, you must set the format of the file.'
         stop
!
      elseif (the_file%access == 'direct') then

!
         write(output%unit,'(/t3,a)') 'Error: tried to open sequential access file as a direct access file.'
         stop
!
      endif
!
!     Open file
!
      if (present(pos)) then
!
         io_error = -1
!
         open(newunit=the_file%unit, file=the_file%name, access='sequential', &
              action=permissions, status='unknown', form=the_file%format, position=pos, iostat=io_error)
!
      else
!
         io_error = -1
!
         open(newunit=the_file%unit, file=the_file%name, access='sequential', &
              action=permissions, status='unknown', form=the_file%format, iostat=io_error)
!
      endif
!
!     Check whether file open was successful
!
      if (io_error .ne. 0) then
!
         write(output%unit,'(/t3,a)') 'Error: could not open file: ', the_file%name
         stop
!
      endif
!
!     Tell the file that it is now open, and ask it to calculate its own size
!
      the_file%opened = .true.
!
      call disk%determine_file_size(the_file)
!
!     If the intent is 'write' or 'readwrite' and the disk is entirely filled
!     (according to the specified available disk space), the calculation will stop:
!
      if (disk%available .lt. 0 .and. (permissions == 'write' .or. permissions == 'readwrite')) then
!
         write(output%unit,'(t3,a/,t3,a)') 'Error: the specified disk space is used up and', &
                                           'a file was opened with permission to write.'
         stop
!
      endif
!
   end subroutine open_file_sequential_disk_manager
!
!
   subroutine open_file_direct_disk_manager(disk, the_file, permissions)
!!
!!    Open direct access file
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2018
!!
!!    Opens a direct access file object.
!!
!!    The routine takes the following arguments:
!!
!!       - the_file (an object of type "file"). It is assumed that a file name
!!         has been set: i.e., the_file%name = 'filename'.
!!       - permissions ('read', 'write', 'readwrite')
!!
      implicit none
!
      class(disk_manager) :: disk
!
      class(file) :: the_file ! the file
!
      character(len=*) :: permissions
!
      integer(i15) :: io_error = -1
!
!     Sanity checks
!
      if (the_file%name == 'no_name') then
!
         write(output%unit,'(/t3,a)') 'Error: to open a file, you must set the name of the file.'
         stop
!
      elseif (the_file%format == 'unknown') then
!
         write(output%unit,'(/t3,a)') 'Error: to open a file, you must set the format of the file.'
         stop
!
      elseif (the_file%access == 'sequential') then
!
         write(output%unit,'(/t3,a)') 'Error: tried to open direct access file as a sequential access file.'
         stop
!
      elseif (the_file%record_length == 0) then
!
         write(output%unit,'(/t3,a)') 'Error: tried to open direct access file without a set record length.'
         stop
!
      endif
!
!     Open file
!
      io_error = -1
!
      open(newunit=the_file%unit, file=the_file%name, access='direct', &
           action=permissions, status='unknown', form=the_file%format, recl=the_file%record_length, iostat=io_error)
!
!
!     Check whether file open was successful
!
      if (io_error .ne. 0) then
!
         write(output%unit,'(/t3,a)') 'Error: could not open file: ', the_file%name
         stop
!
      endif
!
!     Tell the file that it is now open, and ask it to calculate its own size
!
      the_file%opened = .true.
!
      call disk%determine_file_size(the_file)
!
!     If the intent is 'write' or 'readwrite' and the disk is entirely filled
!     (according to the specified available disk space), the calculation will stop:
!
      if (disk%available .lt. 0 .and. (permissions == 'write' .or. permissions == 'readwrite')) then
!
         write(output%unit,'(t3,a/,t3,a)') 'Error: the specified disk space is used up and', &
                                           'a file was opened with permission to write.'
         stop
!
      endif
!
   end subroutine open_file_direct_disk_manager
!
!
   subroutine close_file_disk_manager(disk, the_file, destiny)
!!
!!    Close file
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2018
!!
      implicit none
!
      class(disk_manager) :: disk
!
      class(file) :: the_file
!
      character(len=*), optional :: destiny ! i.e. 'status' after close, can be 'keep' or 'delete'
!
      integer(i15) :: file_size_when_opened
      integer(i15) :: file_size_when_closed
!
      integer(i15) :: bytes_written_to_disk ! Negative if storage is freed up, via 'delete'
!
!
!     Sanity check
!
      if (.not. the_file%opened) then
!
         write(output%unit,'(t3,a)') 'Error: tried to close a file that has not been opened.'
         stop
!
      endif
!
!     Get the file size (both when opened & closed)
!
      file_size_when_opened = the_file%size
!
      call disk%determine_file_size(the_file)
!
      file_size_when_closed = the_file%size
!
!     Close file
!
      if (present(destiny)) then ! Either keep or delete
!
         if (.not. (destiny == 'keep' .or. destiny == 'delete')) then
!
            write(output%unit,'(t3,a)') 'Error: could not recognize status when closing file.'
            write(output%unit,'(t3,a)') 'Must equal keep or delete.'
            stop
!
         else
!
            if (destiny == 'keep') then
!
                  close(the_file%unit, status=destiny)
                  bytes_written_to_disk = file_size_when_closed - file_size_when_opened
!
            else ! destiny == 'delete'
!
                  close(the_file%unit, status=destiny)
                  bytes_written_to_disk = -file_size_when_closed
!
!                 Sanity check
!
                  if (file_size_when_opened .ne. file_size_when_closed) then
!
                     write(output%unit,'(t3,a)') 'Warning: deleting a file that has been written to since'
                     write(output%unit,'(t3,a)') 'it was opened. To avoid an apparent accumulation of storage space,'
                     write(output%unit,'(t3,a)') 'the estimated freed up space is taken to be the initial file size.'
                     bytes_written_to_disk = -file_size_when_opened
!
                  endif
!
            endif
!
!
         endif
!
      else ! Keep file
!
         close(the_file%unit)
         bytes_written_to_disk = file_size_when_closed - file_size_when_opened
!
      endif
!
!     Let the user know how much has been written, or freed up, and to which file
!
!       if (bytes_written_to_disk .gt. 0) then
! !
!          write(output%unit,'(/t3,a,a,a,i14)') 'Number of bytes written to file ', &
!                            trim(the_file%name), ': ', bytes_written_to_disk
! !
!       elseif (bytes_written_to_disk .lt. 0) then
! !
!          write(output%unit,'(/t3,a,a,a,i14)') 'File ', &
!                         trim(the_file%name),  ' modified or deleted, with bytes freed up: ', bytes_written_to_disk
! !
!       endif
!
!     Update the available disk space
!
      disk%available = disk%available - bytes_written_to_disk
!
      if (disk%available .lt. 0) then
!
         write(output%unit,'(t3,a/,t3,a)') 'Warning: the specified disk space is now used up. If any', &
                                           'more has to be stored, the calculation will stop.'
!
      endif
!
!     Set the status for the file to 'closed'
!
      the_file%opened = .false.
!
   end subroutine close_file_disk_manager
!
!
   subroutine determine_file_size_disk_manager(disk, the_file)
!!
!!    Determine file size
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2018
!!
!!    The disk manager handles files. This routine is called by it
!!    and should never be called by the user (because it can lead to
!!    errors in the disk space estimates).
!!
      implicit none
!
      class(disk_manager) :: disk
!
      class(file) :: the_file
!
!     Inquire about the file size
!
      inquire(file=the_file%name, size=the_file%size)
!
!     Check whether the file size could be calculated
!
      if (the_file%size .eq. -1) then
!
         write(output%unit,*) 'Error: Could not calculate file size of the file ', trim(the_file%name)
         stop
!
      endif
!
   end subroutine determine_file_size_disk_manager
!
!
   subroutine read_settings_disk_manager(disk)
!!
!!    Read settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!  
      class(disk_manager) :: disk
!
      integer(i15) :: n_specs, i
!
      character(len=100) :: line
!
      disk%total = 8.0d0
!
      call move_to_section('disk', n_specs)
!
      do i = 1, n_specs
!
         read(input%unit, '(a100)') line
         line = remove_preceding_blanks(line)
!
         if (line(1:10) == 'available:' ) then
!
            read(line(11:100), *) disk%total
!
         endif
!
      enddo
!
      disk%total = disk%total*1000000000
!
   end subroutine read_settings_disk_manager
!
!
   subroutine print_settings_disk_manager(disk)
!!
!!    Print settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!  
      class(disk_manager) :: disk
!
      write(output%unit, '(t3, a38, i3, a)') 'Disk space available for calculation: ',&
                                                 disk%total/1000000000, ' GB'
!
   end subroutine print_settings_disk_manager
!
!
   subroutine delete_disk_manager(disk, the_file)
!!
!!    Delete file
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!    
!
      implicit none
!  
      class(disk_manager) :: disk
!
      type(file) :: the_file
!
      logical :: file_exists
!
      inquire(file=the_file%name, exist=file_exists)
!
      if (file_exists) then
!
         call disk%open_file(the_file, 'write')
         call disk%close_file(the_file, 'delete')
!
      endif
!
   end subroutine delete_disk_manager
!
!
   function file_exists_disk_manager(disk, the_file)
!!
!!    File exists
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!    
!
      implicit none
!  
      class(disk_manager), intent(in) :: disk
!
      type(file), intent(in) :: the_file
!
      logical :: file_exists_disk_manager
!
      inquire(file=the_file%name, exist=file_exists_disk_manager)
!
   end function file_exists_disk_manager
!
end module disk_manager_class
