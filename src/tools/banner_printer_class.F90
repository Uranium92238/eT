!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2022 the authors of eT
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
module banner_printer_class
!
!!
!!    Banner printer class
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, and Alexander C. Paul, 2022
!!
!
   use global_out, only: output, timing
!
   implicit none
!
   type :: banner_printer 
!
      character(len=1), private :: print_level = 'm'
!
   contains
!
      procedure, public :: print_ 
!
      procedure, private :: print_to_file
!
   end type banner_printer
!
contains
!
!
   subroutine print_(this, header)
!!
!!    Print
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, and Alexander C. Paul, 2022
!!
      implicit none 
!
      class(banner_printer), intent(in) :: this 
!
      character(len=*), intent(in) :: header 
!
      call this%print_to_file(output, header)
      call this%print_to_file(timing, header)
!
   end subroutine print_ 
!
!
   subroutine print_to_file(this, file_, header)
!!
!!    Print to file
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, and Alexander C. Paul, 2022
!!
      use output_file_class, only: output_file
!
      implicit none 
!
      class(banner_printer), intent(in) :: this 
!
      class(output_file), intent(inout) :: file_ 
!
      character(len=*), intent(in) :: header 
!
      call file_%printf(this%print_level, trim(header), fs='(//t3,a)')
      call file_%print_separator(this%print_level, len(trim(header)), '=')
!
   end subroutine print_to_file
!
!
end module banner_printer_class
