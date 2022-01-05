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
module timings_file_class
!
!!
!!    Timings file class module
!!    Written by Sarai D. Folkestad, Sep 2021
!!
!!
!
   use kinds
   use output_file_class, only: output_file
!
   type, extends(output_file) :: timings_file
!
   contains
!
      procedure, public :: print_banner => print_banner_timings_file
      procedure, public :: print_time => print_time_timings_file
      procedure, private :: print_formated_task_name
      procedure, nopass, private :: get_format
!
   end type timings_file
!
   interface timings_file
!
      procedure new_timings_file
!
   end interface timings_file
!
contains
!
!
   function new_timings_file(name_) result(the_file)
!!
!!    New
!!    Writen by Rolf H. Myhre May 2019
!!
      implicit none
!
      type(timings_file) :: the_file
!
      character(len=*), intent(in) :: name_
!
      call the_file%initialize(name_)
!
   end function new_timings_file
!
!
   subroutine print_banner_timings_file(the_file)
!!
!!    Print banner
!!    Writtan by Sarai D. Folkestad, Sep 2021
!!
      implicit none
!
      class(timings_file) :: the_file
!
      if (.not. the_file%is_open) call the_file%error_msg('tried to print to closed timings file')
!
!
      call the_file%printf('m', 'Description                                 wall[s]         &
                                 &cpu[s]    cpu/wall', fs='(/t3,a)', ll=200)
      call the_file%print_separator('m', 85, '-', fs='(t3,a)')
!
   end subroutine print_banner_timings_file
!
!
   subroutine print_time_timings_file(the_file, wall, cpu, name_, pl)
!!
!!    Print times
!!    Written by Eirik F. Kjønstad, Dec 2018
!!
      implicit none
!
      class(timings_file), intent(inout) :: the_file
!
      real(dp), intent(in) :: wall, cpu
      character(len=*), intent(in) :: name_
      character(len=100), intent(in) :: pl
!
      call print_formated_task_name(the_file, name_, pl)
!
      if (wall > 1.0d-8) then
!
         call the_file%printf(pl, ' (f10.2)  (f14.2)  (f8.2)', fs='(t3,a)', &
                        reals=[wall, cpu, cpu/wall], ll=100)
!
      else
!
         call the_file%printf(pl, ' (f10.2)  (f14.2)        --', fs='(t3,a)', &
                           reals=[wall, cpu], ll=100)
!
      endif
!
   end subroutine print_time_timings_file
!
!
   subroutine print_formated_task_name(the_file, name_, pl)
!!
!!    Print fromated task name 
!!    Written by Sarai D. Folkestad, Sep 2021
!!
!!    Prints the timing description (name_)
!!    If needed, over several lines. Line breaks at spaces are 
!!    preferred.
!!
      use string_utilities, only: n_instances_of_character, last_instance_of_character
!
      implicit none
!
      class(timings_file), intent(inout) :: the_file
      character(len=*), intent(in) :: name_
      character(len=100), intent(in) :: pl
!
      integer :: name_length, offset, chunk_length, cursor
!
      character(len=100) :: name_chunk
!
      character(len=20) :: format_string_name
!
      integer, parameter :: max_len = 35
!
      name_length = len(trim(name_))
!
      if (name_length .gt. max_len) then
!
         offset = 0
         name_chunk = name_
!
         do while (len_trim(name_chunk) .gt. max_len)
!
            cursor = max_len
            if (n_instances_of_character(trim(name_chunk(1:cursor)), ' ') > 0) &
               cursor = last_instance_of_character(trim(name_chunk(1:cursor)), ' ')
!
            name_chunk = name_chunk(1:cursor)
!
            chunk_length = len_trim(name_chunk)
            format_string_name = the_file%get_format(chunk_length, max_len)
            call the_file%printf(pl, trim(format_string_name), chars=[trim(name_chunk)])
!
            offset = offset + cursor
            name_chunk = adjustl(name_(offset + 1:))
!
         enddo
!
         chunk_length = len_trim(name_chunk)
         format_string_name = the_file%get_format(chunk_length, max_len)
         call the_file%printf(pl, trim(format_string_name), adv=.false., chars=[trim(name_chunk)])
!
      else
!
         format_string_name = the_file%get_format(name_length, max_len)
         call the_file%printf(pl, trim(format_string_name), adv=.false., chars=[name_])
!
      endif

!
   end subroutine print_formated_task_name
!
!
   function get_format(length, max_length) result(format_)
!!
!!    Get format
!!    Written by Sarai D. Folkestad, Sep 2021
!!
      implicit none
!
      integer, intent(in) :: length, max_length
!
      character(len=20) :: format_
!
      if (max_length - length .gt. 0) then
         write(format_, '(a2,i0,a4,i0,a1)') &
              '(a', length, ') (x', max_length-length,')'
      else
         write(format_, '(a2,i0,a1)') '(a', max_length, ')'
      endif
!
   end function get_format
!
end module timings_file_class
