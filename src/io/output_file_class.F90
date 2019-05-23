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
   use abstract_file_class, only : abstract_file 
!
   type, extends(abstract_file) :: output_file
!
!
   contains
!
      procedure :: open_file              => open_file_output_file
      procedure :: close_file             => close_file_output_file
      procedure :: flush_file             => flush_file_output_file
!
      procedure :: error_msg              => error_msg_output_file
      procedure :: warning_msg            => warning_msg_output_file
!
   end type output_file
!
   interface output_file
!
      procedure new_output_file
!
   end interface output_file
!
   type(output_file) :: output
!
contains
!
!
   module function new_output_file(file_name) result(the_file)
!!
!!    Output file constructer
!!    Writen by Rolf H. Myhre May 2019
!!
!!    Output file is formatted and sequantial file.
!!    Routine sets these, and sets the file name    
!!
      implicit none
!
      type(output_file) :: the_file
!
      character(len=*), intent(in) :: file_name
!
      the_file%file_name = file_name
!
      the_file%file_access = 'sequential'
      the_file%file_format = 'formatted'
!
   end function
!
!
   subroutine open_file_output_file(the_file)
!!
!!    Open the output file
!!    Written by Rolf Heilemann Myhre, May 2019
!!
      implicit none
!
      class(output_file) :: the_file
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      open(newunit=the_file%unit, file=the_file%file_name, access=the_file%file_access, &
           action='write', status='unknown', form=the_file%file_format, iostat=io_error, iomsg=io_msg)
!
      if (io_error /= 0) stop 'Error: could not open eT output file '//trim(the_file%file_name)//&
                             &'. Error message: '//trim(io_msg)
!
      the_file%file_opened = .true.
!
      call the_file%set_open_file_size()
!
   end subroutine open_file_output_file
!
!
   subroutine close_file_output_file(the_file)
!!
!!    Close the output file
!!    Written by Rolf Heilemann Myhre, May 2019
!!
      implicit none
!
      class(output_file) :: the_file
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      close(the_file%unit, iostat=io_error, iomsg=io_msg, status='keep')
!
      if (io_error.ne. 0) then
          stop 'Error: could not close eT output file '//trim(the_file%file_name)//&
              &'error message: '//trim(io_msg)
      endif
!
      the_file%file_opened = .false.
!
   end subroutine close_file_output_file
!
!
   subroutine flush_file_output_file(the_file)
!!
!!    Flush the output file
!!    Written by Rolf Heilemann Myhre, May 2019
!!
      implicit none
!
      class(output_file) :: the_file
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      flush(the_file%unit, iostat=io_error, iomsg=io_msg)
!
      if (io_error /= 0) stop 'Error: could not flush eT output file '//trim(the_file%file_name)//&
                             &'error message: '//trim(io_msg)
!
   end subroutine flush_file_output_file
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
      stop "Something went wrong, check the .out file"
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
