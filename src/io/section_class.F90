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
   use output_file_class
   use string_utilities  
!
   type :: section
!
      character(len=:), allocatable :: name_
!
      character(len=21), dimension(:), allocatable :: keywords
!
   contains
   end type section  
!
end module section_class
