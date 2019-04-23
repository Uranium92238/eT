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
module disk_manager_class
!
!!
!!    Disk manager class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2018
!!
!
   use kinds
   use file_class
   use input_file_class
   !use io_utilities
!
   type :: disk_manager
!
!     The total amount of disk space specified by user (standard: 30 GB)
!
      integer(i15), private :: total
!
!     The amount of disk space currently available, based on the files
!     currently stored on file
!
      integer(i15), private :: available
!
   contains
!
      procedure :: prepare                      => prepare_disk_manager
!
      procedure :: open_file_sequential         => open_file_sequential_disk_manager
      procedure :: open_file_direct             => open_file_direct_disk_manager
      procedure :: open_file                    => open_file_disk_manager
!
      procedure :: close_file                   => close_file_disk_manager
!
      procedure :: read_settings                => read_settings_disk_manager
      procedure :: print_settings               => print_settings_disk_manager
!
      procedure :: delete                       => delete_disk_manager
!
   end type disk_manager
!
   type(disk_manager) :: disk
!
contains
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
      disk%total = 30000000000_i15
      call disk%read_settings()
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
!     Sanity checks
!
      if ( present(pos)) then
         if (the_file%access  == 'direct') then
!
            call output%warning_msg('position specifier is disregarded for direct access file.')
!
         endif
      elseif (the_file%access .ne. 'direct' .and. the_file%access .ne. 'sequential' ) then
!
         call output%error_msg('illegal access type for file: ' // trim(the_file%name) //',' // trim(the_file%access))
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
      integer :: io_error = -1
!
!     Sanity checks
!
      if (the_file%name == 'no_name') then
!
         call output%error_msg('to open a file, you must set the name of the file.')
!
      elseif (the_file%format == 'unknown') then
!
         call output%error_msg('to open a file, you must set the format of the file.')
!
      elseif (the_file%access == 'direct') then

!
         call output%error_msg('tried to open sequential access file as a direct access file.')
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
         call output%error_msg('could not open file: ' //  trim(the_file%name))
!
      endif
!
!     Tell the file that it is now open, and ask it to calculate its own size
!
      the_file%opened = .true.
!
      call the_file%determine_file_size()
!
!     If the intent is 'write' or 'readwrite' and the disk is entirely filled
!     (according to the specified available disk space), the calculation will stop:
!
      if (disk%available .lt. 0 .and. (permissions == 'write' .or. permissions == 'readwrite')) then
!
         call output%error_msg('the specified disk space is used up and' // &
                                           'a file was opened with permission to write.')
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
      integer :: io_error = -1
!
!     Sanity checks
!
      if (the_file%name == 'no_name') then
!
         call output%error_msg('to open a file, you must set the name of the file.')
!
      elseif (the_file%format == 'unknown') then
!
         call output%error_msg('to open a file, you must set the format of the file.')
!
      elseif (the_file%access == 'sequential') then
!
         call output%error_msg('tried to open direct access file as a sequential access file.')
!
      elseif (the_file%record_length == 0) then
!
         call output%error_msg('tried to open direct access file without a set record length.')
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
         call output%error_msg('could not open file: ' // trim(the_file%name))
!
      endif
!
!     Tell the file that it is now open, and ask it to calculate its own size
!
      the_file%opened = .true.
!
      call the_file%determine_file_size()
!
!     If the intent is 'write' or 'readwrite' and the disk is entirely filled
!     (according to the specified available disk space), the calculation will stop:
!
      if (disk%available .lt. 0 .and. (permissions == 'write' .or. permissions == 'readwrite')) then
!
         call output%error_msg('Error: the specified disk space is used up and ' // &
                                           'a file was opened with permission to write.')
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
      integer :: file_size_when_opened
      integer :: file_size_when_closed
!
      integer :: bytes_written_to_disk ! Negative if storage is freed up, via 'delete'
!
      bytes_written_to_disk = 0
!
!     Sanity check
!
      if (.not. the_file%opened) then
!
         call output%error_msg('tried to close a file that has not been opened.')
!
      endif
!
!     Get the file size (both when opened & closed)
!
      file_size_when_opened = the_file%get_file_size()
!
      call the_file%determine_file_size()
!
      file_size_when_closed = the_file%get_file_size()
!
!     Close file
!
      if (present(destiny)) then ! Either keep or delete
!
         if (.not. (destiny == 'keep' .or. destiny == 'delete')) then
!
            call output%error_msg('could not recognize status when closing file.')
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
                  !  write(output%unit,'(t3,a)') 'Warning: deleting a file that has been written to since'
                  !  write(output%unit,'(t3,a)') 'it was opened. To avoid an apparent accumulation of storage space,'
                  !  write(output%unit,'(t3,a)') 'the estimated freed up space is taken to be the initial file size.'
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
!     Update the available disk space
!
      disk%available = disk%available - bytes_written_to_disk
!
      if (disk%available .lt. 0) then
!
         call output%warning_msg('the specified disk space is now used up. If any ' //&
                                           'more has to be stored, the calculation will stop.')
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
   subroutine read_settings_disk_manager(disk)
!!
!!    Read settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!  
      class(disk_manager) :: disk
!
      if (input%requested_keyword_in_section('available','disk')) then
!
         call input%get_keyword_in_section('available', 'disk', disk%total)
         disk%total = disk%total*1000000000
!
      endif 
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
      write(output%unit, '(t3, a38, i5, a)') 'Disk space available for calculation: ',&
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
      class(file) :: the_file
!
      logical :: file_exists
!
      file_exists = the_file%file_exists()
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
end module disk_manager_class
