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
   use global_in,  only : input
   use global_out, only : output
!
   type :: disk_manager
!
!     The total amount of disk space specified by user (standard: 30 GB)
!
      integer(i15), private :: total = 300000
!
!     The amount of disk space currently available, based on the files
!     currently stored on file
!
      integer(i15), private :: available = 300000
!
   contains
!
      procedure :: prepare                      => prepare_disk_manager
!
      procedure :: check_available              => check_available_disk_manager
      procedure :: update                       => update_disk_manager
!
      procedure :: read_settings                => read_settings_disk_manager
      procedure :: print_settings               => print_settings_disk_manager
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
      disk%total = 300000000000_i15
      call disk%read_settings()
!
      disk%available = disk%total
!
      call disk%print_settings()
!
   end subroutine prepare_disk_manager
!
!
   function check_available_disk_manager(disk, disk_used) result(room)
!!
!!    Check available
!!    Written by Rolf Heilemann Myhre, May 2019
!!    Check if room for input integer on disk 
!!
      implicit none
!
      class(disk_manager), intent(in)  :: disk
      integer, intent(in)              :: disk_used
      logical                          :: room
!
      if(disk%available - disk_used .lt. 0) then
         room = .false.
      else
         room = .true.
      endif
!
   end function check_available_disk_manager
!
!
   subroutine update_disk_manager(disk, disk_used, name_)
!!
!!    Update disk manager
!!    Written by Rolf Heilemann Myhre, May 2019
!!    Check if disk space available and update if there is
!!
      implicit none
!
      class(disk_manager), intent(inout)  :: disk
      integer, intent(in)  :: disk_used
      character(len=*), intent(in) :: name_
!
      logical  :: room
!
      room = disk%check_available(disk_used)
!
      if (room) then
         disk%available = disk%available - disk_used
      else
         call output%error_msg('File '//trim(name_)//' has used too much disk space')
      endif
!
   end subroutine update_disk_manager
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
end module disk_manager_class
