module disk_manager_class
!
!!
!!                            Disk manager class module                                 
!!             Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2018         
!!
!! 
!
!  :::::::::::::::::::::::::::::::::::
!  -::- Modules used by the class -::-
!  :::::::::::::::::::::::::::::::::::
!
   use types
   use input_output
!
!  ::::::::::::::::::::::::::::::::::::::::::::::
!  -::- Definition of the disk_manager class -::-
!  ::::::::::::::::::::::::::::::::::::::::::::::
!
   type :: disk_manager 
!
!     The total amount of disk space specified by user (standard: 30 GB)
!
      integer(i15) :: total = 30000000000
!
!     The amount of disk space currently available, based on the files 
!     currently stored on file
!
      integer(i15) :: available = 30000000000
!
   contains 
!
!     Initialization routine (used if user specifies a disk space different from standard)
!
      procedure :: init => init_disk_manager 
!
   end type disk_manager                                                                            
!
!
contains
!
!
   subroutine init_disk_manager(disk, total)
!!
!!    Init (disk manager)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2018
!!
!!    Initializes the memory manager object by setting the 
!!    total and initial available memory. This is only called 
!!    if the user specifies a total memory different from the standard.
!! 
      implicit none 
!
      class(disk_manager) :: disk 
!
      integer(i15) :: total ! In GBs
!     
!     Set the specified total memory in bytes 
!
      disk%total = total*1.0D9
!
!     Update the initially available memory 
!
      disk%available = disk%total 
!
   end subroutine init_disk_manager
!
!
end module disk_manager_class
