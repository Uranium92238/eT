!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2021 the authors of eT
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
module task_list_class
!
!!
!! Task list class
!! Written by Eirik F. Kjønstad and Sarai D. Folkestad, Dec 2019   
!!
!! Handles task prints for the engines.
!!
!! Should be allocated in the 'set_printables' routine
!! of the engine. There all tasks should be added to the 
!! task list through the use of 'add_'
!!
!! In the banner, all tasks are printed using 'print_all'.
!!
!! In each subtask of the engine, the current task is printed
!!
!
   use global_out, only : output
!
   implicit none
!
   type task_list 
!
      private 
!
      integer :: n_tasks
!
      character(len=150), dimension(:), allocatable :: descriptions
      character(len=150), dimension(:), allocatable :: labels
!
      integer :: max_tasks
!
   contains 
!
      procedure :: add        => add_task_list
      procedure :: print_     => print_task_list
      procedure :: print_all  => print_all_task_list
!
      final :: destructor_task_list
!
   end type task_list 
!
   interface task_list
!
      procedure :: new_task_list
!
   end interface task_list 
!
contains
!
!
   function new_task_list() result(tasks)
!!
!!    New task list 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Dec 2019
!!
      implicit none
!
      type(task_list) :: tasks
!
      tasks%n_tasks = 0
      tasks%max_tasks = 20
!
      allocate(tasks%descriptions(tasks%max_tasks))
      allocate(tasks%labels(tasks%max_tasks))
!
   end function new_task_list
!
!
   subroutine add_task_list(tasks, label, description)
!!
!!    Add task
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Dec 2019
!!
!!    Adds a new task to the task list
!!
      implicit none
!
      class(task_list), intent(inout) :: tasks
!
      character(len=*), intent(in) :: label
      character(len=*), intent(in) :: description
!
      tasks%n_tasks = tasks%n_tasks + 1
!
      if (tasks%n_tasks .gt. tasks%max_tasks) &
         call output%error_msg('exceeded maximal number of tasks')
!
      tasks%labels(tasks%n_tasks) = label
      tasks%descriptions(tasks%n_tasks) = description
!
   end subroutine add_task_list
!
!
   subroutine print_task_list(tasks, label, append_string, append_fs)
!!
!!    Print
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Dec 2019
!!
!!    Prints the description matching the 
!!    passed 'label'
!!
!!    'append_string' : optional argumet to add a new line of text under 'description'
!!
!!    'append_fs' : optional argument for format of appended string
!!
      use timings_class, only: timing
!
      implicit none
!
      class(task_list), intent(in) :: tasks
!
      character(len=*), intent(in) :: label
!
      character(len=*), intent(in), optional :: append_string
      character(len=*), intent(in), optional :: append_fs
!
      integer :: i
!
      do i = 1, tasks%n_tasks
!
         if (label == trim(tasks%labels(i))) then
!
            call output%printf('m', '(i0)) ' // trim(tasks%descriptions(i)), &
                               ints=[i], fs='(//t3,a)')
!
            call timing%printf('m', '(i0)) ' // trim(tasks%descriptions(i)), &
                               ints=[i], fs='(//t3,a)')
!
            if (present(append_string)) then
!
               if (present(append_fs)) then
!
                  call output%printf('m', append_string, fs=append_fs)
                  call timing%printf('m', append_string, fs=append_fs)
!
               else
!
                  call output%printf('m', append_string, fs='(//t3,a)')
                  call timing%printf('m', append_string, fs='(//t3,a)')
!
               endif
!
            endif
!
            return
!
         endif
!
      enddo
!
      call output%error_msg('did not recognize task label')
!
   end subroutine print_task_list
!
!
   subroutine print_all_task_list(tasks)
!!
!!    Print all 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad
!!
!!    Prints all the tasks to the output file
!!
      implicit none
!
      class(task_list), intent(in) :: tasks
!
      integer :: i
!
      do i = 1, tasks%n_tasks
!
         call output%printf('m', '(i0)) ' // trim(tasks%descriptions(i)), &
                            ints=[i], fs='(t6,a)')
!
      enddo
!
   end subroutine print_all_task_list
!
!
   subroutine destructor_task_list(tasks)
!!
!!    Destructor
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Dec 2019
!!
      implicit none
!
      type(task_list) :: tasks
!
      deallocate(tasks%labels)
      deallocate(tasks%descriptions)
!
   end subroutine destructor_task_list
!
!
end module task_list_class
