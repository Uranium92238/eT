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
      procedure :: determine_file_size    => determine_file_size_abstract_file
      procedure :: get_file_size          => get_file_size_abstract_file
      procedure :: file_exists            => file_exists_abstract_file
!
   end type abstract_file
!
contains
!
!
   subroutine determine_file_size_abstract_file(the_file)
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
      class(abstract_file) :: the_file
!
!     Inquire about the file size
!
      inquire(file=the_file%name, size=the_file%file_size)
!
!     Check whether the file size could be calculated
!
      if (the_file%file_size .eq. -1) then
!
         write(*,'(/a, a/)') 'Error: Could not determine size of file '// trim(the_file%name)
         stop
!
      endif
!
   end subroutine determine_file_size_abstract_file
!
!
   function get_file_size_abstract_file(the_file)
!!
!!    Return private variable file_size
!!    Written by Rolf H. Myhre, 2018
!!    
!
      implicit none
!  
      class(abstract_file), intent(in) :: the_file
!
      integer :: get_file_size_abstract_file
!
      get_file_size_abstract_file = the_file%file_size
!
   end function get_file_size_abstract_file
!
!  
   function file_exists_abstract_file(the_file)
!!
!!    File exists
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!    Moved to file by Rolf H. Myhre Nov. 2018
!!    
      implicit none
!  
      class(abstract_file), intent(in) :: the_file
!
      logical :: file_exists_abstract_file
!
      inquire(file=the_file%name, exist=file_exists_abstract_file)
!
   end function file_exists_abstract_file
!
!  
end module abstract_file_class
