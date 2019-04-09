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
module output_file_class
!
!!
!!    File class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2018
!!
!!
!
   use kinds    
   use abstract_file_class      
!
   type, extends(abstract_file) :: output_file
!
!
   contains
!
      procedure :: init                   => init_output_file
!
      procedure :: error_msg              => error_msg_output_file
      procedure :: warning_msg            => warning_msg_output_file
!
   end type output_file
!
   type(output_file) :: output, timing
!
contains
!
!
   subroutine init_output_file(the_file, name)
!!
!!    Initialize output file
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2018
!!
!!    Initializes an output file file object
!!
!!    Output file is formatted and sequantial file.
!!    Routine sets these, and sets tha file name    
!!
      implicit none
!
      class(output_file) :: the_file
!
      character(len=*) :: name
!
      the_file%name = name
!
      the_file%access = 'sequential'
      the_file%format = 'formatted'
!
   end subroutine init_output_file
!
!
   subroutine error_msg_output_file(out_file, error_specs, error_int)
!!
!!    Error message
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(output_file) :: out_file
!
      character(len=*) :: error_specs
!
      integer, optional :: error_int 
!
      character(len=40) :: error_int_char = ' '
!
      if (present(error_int)) then
!
         write(error_int_char, '(i12)') error_int   
!
         write(out_file%unit, '(/t3,a)') 'Error: ' // trim(error_specs) // ' ' // error_int_char
!
      else
!
         write(out_file%unit, '(/t3,a)') 'Error: ' // trim(error_specs)
!
      endif
!
      stop
!
   end subroutine error_msg_output_file
!
!
   subroutine warning_msg_output_file(out_file, warning_specs)
!!
!!    Warning message
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(output_file) :: out_file
!
      character(len=*) :: warning_specs
!
      write(out_file%unit, '(/t3,a)') 'Warning: ' // trim(warning_specs)
!
   end subroutine warning_msg_output_file
!
!  
end module output_file_class
