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
      procedure, private :: print_formatted_task_name
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
!!    Written by Rolf H. Myhre May 2019
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
!!    Written by Sarai D. Folkestad, Sep 2021
!!
      implicit none
!
      class(timings_file) :: the_file
!
      if (.not. the_file%is_open) &
         call the_file%error_msg('tried to print to closed timings file')
!
      call the_file%printf('m', 'Description(x41)wall[s](x10)cpu[s]  cpu/wall', &
                           fs='(/t3,a)', ll=100)
!
      call the_file%print_separator('m', 85, '-', fs='(t3,a)')
!
   end subroutine print_banner_timings_file
!
!
   subroutine print_time_timings_file(the_file, wall, cpu, name_, pl)
!!
!!    Print time
!!    Written by Eirik F. Kjønstad, Dec 2018
!!
      implicit none
!
      class(timings_file), intent(inout) :: the_file
!
      real(dp), intent(in) :: wall, cpu
      character(len=*), intent(in) :: name_
      character(len=*), intent(in) :: pl
!
      call print_formatted_task_name(the_file, name_, pl)
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
   subroutine print_formatted_task_name(the_file, name_, pl)
!!
!!    Print formatted task name
!!    Written by Sarai D. Folkestad, Sep 2021
!!
!!    Prints the timing description (name_)
!!    If needed, over several lines. Line breaks at spaces are preferred.
!!
      use string_utilities, only: n_instances_of_character, last_instance_of_character
!
      implicit none
!
      class(timings_file), intent(inout) :: the_file
      character(len=*), intent(in) :: name_
      character(len=*), intent(in) :: pl
!
      integer :: name_length, offset, cursor
!
      character(len=100) :: remaining_name
!
      integer, parameter :: max_length = 46
      character(len=5) :: format_string, local_format_string
!
      write(format_string, '(a,i0,a)') "(b", max_length, ")"
!
      name_length = len_trim(name_)
!
      if (name_length <= max_length) then
         call the_file%printf(pl, format_string, adv=.false., chars=[name_])
         return
      end if
!
      offset = 1
      remaining_name = name_
!
      do while (len_trim(remaining_name) .gt. max_length)
!
         cursor = max_length
         if (n_instances_of_character(trim(remaining_name(1:cursor)), ' ') > 0) &
            cursor = last_instance_of_character(trim(remaining_name(1:cursor)), ' ')
!
         write(local_format_string, '(a,i0,a)') "(b", cursor, ")"
         call the_file%printf(pl, trim(local_format_string), &
                              chars=[trim(remaining_name(1:cursor))])
!
         offset = offset + cursor
         remaining_name = adjustl(name_(offset:))
!
      enddo
!
      call the_file%printf(pl, format_string, adv=.false., chars=[trim(remaining_name)])
!
   end subroutine print_formatted_task_name
!
!
end module timings_file_class
