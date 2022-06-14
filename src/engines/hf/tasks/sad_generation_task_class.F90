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
module sad_generation_task_class
!
!!
!! SAD generation task class
!! Written by Alexander C. Paul, Sarai D. Folkestad, and Eirik F. Kjønstad, 2022
!!
!
   use kinds
   use hf_class,      only: hf
   use hf_task_class, only: hf_task
!
   implicit none
!
   type, extends(hf_task) :: sad_generation_task
!
   contains
!
      procedure, public :: execute &
                        => execute_sad_generation_task
!
      procedure, private, nopass :: get_gradient_threshold
!
   end type sad_generation_task
!
!
   interface sad_generation_task
!
      procedure :: new_sad_generation_task
!
   end interface sad_generation_task
!
!
contains
!
!
   function new_sad_generation_task() result(this)
!!
!!    New
!!    Written by Eirik F. Kjønstad, 2022
!!
      implicit none
!
      type(sad_generation_task) :: this
!
      this%name_ = 'Generating initial SAD density'
!
   end function new_sad_generation_task
!
!
   subroutine execute_sad_generation_task(this, wf)
!!
!!    Execute
!!    Written by Sarai D. Folkestad, Oct 2021
!!
      use sad_tool_class, only: sad_tool
      use atomic_center_class, only: atomic_center
!
      implicit none
!
      class(sad_generation_task), intent(inout) :: this
      class(hf), target, intent(inout) :: wf
!
      type(atomic_center), dimension(:), allocatable :: centers
!
      type(sad_tool) :: sad
!
      call this%print_header()
      call this%start_timer()
!
      sad = sad_tool(this%get_gradient_threshold())
!
      allocate(centers(wf%ao%get_n_centers()))
      call wf%ao%get_centers(centers, 1, wf%ao%get_n_centers())
!
      call sad%generate(wf%ao%get_n_centers(), centers)
!
      deallocate(centers)
!
!     Libint is overwritten by SAD. Re-initialize.
      call wf%ao%export_centers_to_libint()
!
!     Re-determine status of a file because SAD may have deleted it
!     (so the status must go from "old" -> "new")
      call wf%orbital_file%determine_status()
!
      call this%end_timer()
!
   end subroutine execute_sad_generation_task
!
!
   function get_gradient_threshold() result(threshold)
!!
!!    Get gradient threshold
!!    Written by Alexander C. Paul, 2022
!!
      use global_in, only: input
!
      implicit none
!
      real(dp) :: threshold
!
      threshold = 1.0d-6
      call input%get_keyword('gradient threshold', 'solver scf', threshold)
!
   end function get_gradient_threshold
!
!
end module sad_generation_task_class
   