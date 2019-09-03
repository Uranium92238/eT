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
module section_class
!
!!
!!    Section class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
!
   use kinds    
   use string_utilities  
   use global_out, only: output
!
   type :: section
!
      character(len=:), allocatable :: name_
!
      logical :: required
!
      character(len=21), allocatable :: keywords(:)
!
   contains
!
      procedure :: print_keywords 
!
   end type section  
!
contains 
!
!
   subroutine print_keywords(the_section)
!!
!!    Print keywords 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
      implicit none 
!
      class(section) :: the_section
!
      integer :: k 
!
      write(output%unit, '(/t3,a,a,a/)') 'The valid keywords in the section named "', trim(the_section%name_), '" are:'
!
      do k = 1, size(the_section%keywords)
!
         write(output%unit, '(t6,a)') the_section%keywords(k)
!
      enddo
!
   end subroutine print_keywords
!
!
end module section_class
