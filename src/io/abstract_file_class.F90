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
module abstract_file_class
!
!!
!!    Abstract file class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2018
!!
!!
!
   use kinds          
!
   type, abstract :: abstract_file
!
!     Filename
!
      character(len=255) :: name_ = 'no_name'
!
!     Unit identifier
!
      integer :: unit = -1
!
!     File size (in bytes)
!
      integer, private :: open_size = -1
      integer, private :: current_size = -1
!
!     Logical for whether the file is currently opened or not
!
      logical :: opened = .false.
!
      character(len=40) :: access_ = 'unknown'
      character(len=40) :: format_ = 'unknown'
!
   contains
!
      procedure :: set_current_size       => set_current_size_abstract_file
      procedure :: set_open_size          => set_open_size_abstract_file
      procedure :: get_size               => get_size_abstract_file
      procedure :: get_change             => get_change_abstract_file
      procedure :: exists                 => exists_abstract_file
!
   end type abstract_file
!
!
contains
!
!
   subroutine set_current_size_abstract_file(the_file)
!!
!!    Determine current file size
!!    Written by Rolf H. Myhre May 2019
!!
      implicit none
!
      class(abstract_file) :: the_file
!
      the_file%current_size = the_file%get_size()
!
   end subroutine set_current_size_abstract_file
!
!
   subroutine set_open_size_abstract_file(the_file)
!!
!!    Determine current file size
!!    Written by Rolf H. Myhre May 2019
!!
      implicit none
!
      class(abstract_file) :: the_file
!
      the_file%open_size = the_file%get_size()
!
   end subroutine set_open_size_abstract_file
!
!
   function get_size_abstract_file(the_file) result(file_size)
!!
!!    Return private variable file_size
!!    Written by Rolf H. Myhre, 2018
!!    
!
      implicit none
!  
      class(abstract_file), intent(in) :: the_file
!
      integer :: file_size
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
!     Inquire about the file size
!
      inquire(file=trim(the_file%name_), size=file_size, &
              iostat=io_error, iomsg=io_msg)
!
!     Check whether the file size could be calculated
!
      if (io_error .ne. 0) then
!
          print *, 'Error: Could not determine size of file ' // trim(the_file%name_)
          stop
!
      endif
!
!     Inquire returns -1 if file is deleted, change to 0      
!
      if (file_size .eq. -1) file_size = 0
!
   end function get_size_abstract_file
!
!  
   function get_change_abstract_file(the_file) result(change)
!!
!!    Return private variable file_size
!!    Written by Rolf H. Myhre, 2018
!!    
!
      implicit none
!  
      class(abstract_file), intent(in) :: the_file
!
      integer :: change
!
      call the_file%set_current_size()
      change = the_file%current_size - the_file%open_size
!
   end function get_change_abstract_file
!
!  
   function exists_abstract_file(the_file)
!!
!!    File exists
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!    Moved to file by Rolf H. Myhre Nov. 2018
!!    
      implicit none
!  
      class(abstract_file), intent(in) :: the_file
!
      logical :: exists_abstract_file
!
      inquire(file=trim(the_file%name_), exist=exists_abstract_file)
!
   end function exists_abstract_file
!
!  
end module abstract_file_class
