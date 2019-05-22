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
module direct_file_class
!
!!
!!    Direct access file class module
!!    Written by Rolf H. Myhre, May 2019
!!
!!
!
   use kinds    
   use abstract_file_class      
   use output_file_class      
!
   type, extends(abstract_file) :: direct_file
!
      integer :: record_length
!
   contains
!
   end type direct_file
!
   interface direct_file
!
      procedure new_direct_file
!
   end interface direct_file
!
contains
!
!
   module function new_direct_file(file_name, record_length) result(the_file)
!!
!!    Direct file constructer
!!    Writen by Rolf H. Myhre May 2019
!!
      type(direct_file) :: the_file
!
      character(len=*), intent(in) :: file_name
      integer, intent(in) :: record_length
!
      the_file%file_name = file_name
!
      the_file%file_access = 'direct'
      the_file%file_format = 'unformatted'
!
      if (record_length .le. 0) then
         call output%error_msg("Record length less than zero for file"//file_name)
      endif
!
      the_file%record_length = record_length 
!
   end function
!
!
end module direct_file_class
