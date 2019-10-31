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
!!    Rewritten by Rolf H. Myhre, November 2018 and September 2019
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
!     Logical for whether the file is currently opened or not
!
      logical :: is_open = .false.
!
      character(len=40) :: access_ = 'unknown'
      character(len=40) :: format_ = 'unknown'
      character(len=40) :: action_ = 'unknown'
!
   contains
!
      procedure :: exists                 => exists_abstract_file
      procedure :: copy                   => copy_abstract_file
!
   end type abstract_file
!
!
contains
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
   subroutine copy_abstract_file(the_file, filename)
!!
!!    Copy abstract file
!!    Written by Alexander C. Paul and Rolf H. Myhre, September 2019
!!
      implicit none
!
      class(abstract_file) :: the_file
!
      character(*), intent(in) :: filename
!
      integer :: copy_unit
!
!     Character to hold a byte
      character :: byte
!
!     Check that file is closed
      if(the_file%is_open) then 
!
         print *, 'Error in copy: '//the_file%name_//' is not closed'
         stop
!
      endif
!
!     Open the file with stream unformatted access
      open(newunit=the_file%unit, file=the_file%name_, access='stream', &
           form='unformatted', action='read')
!
!     Open a new file
      open(newunit=copy_unit, file=trim(filename), access='stream', &
           form='unformatted', action='write', status='new')
!
!
!     Read byte by byte and write it to the new file
!
      do
!
         read(the_file%unit, end=200) byte !Read until end of file, then go to 200
         write(copy_unit) byte             !Write whatever you just read
!
      enddo
      200 continue !End of file reached, should be done
!
!     Close the files
      close(copy_unit, status='keep')
      close(the_file%unit, status='keep')
!
!
   end subroutine copy_abstract_file
!
!
end module abstract_file_class
