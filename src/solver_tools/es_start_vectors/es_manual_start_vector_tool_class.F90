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
module es_manual_start_vector_tool_class
!
!!
!!    Excited state manual start vector tool class module
!!    Written by Eirik F. Kjønstad, 2021
!!
!
   use ccs_class, only: ccs 
!
   use array_utilities, only: get_n_lowest
   use es_start_vector_tool_class
   use memory_manager_class, only: mem 
   use global_out, only: output 
   use global_in, only: input
!
   implicit none
!
   type, extends(es_start_vector_tool) :: es_manual_start_vector_tool
!
   contains
!
   end type es_manual_start_vector_tool
!
!
   interface es_manual_start_vector_tool
!
      procedure :: new_es_manual_start_vector_tool
!
   end interface es_manual_start_vector_tool
!
!
contains 
!
!
   function new_es_manual_start_vector_tool(wf, side, restart) result(this)
!!
!!    New ES manual start vector tool 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Dec 2021
!!
      implicit none 
!
      class(ccs), intent(in), target :: wf
      character(len=*), intent(in) :: side
      logical, intent(in) :: restart
!
      type(es_manual_start_vector_tool) :: this
!
      integer :: length_state_guesses, state
!
      integer, dimension(:), allocatable :: occupied, virtual
!
      this%n_vectors = wf%n_singlet_states
      this%n_parameters = wf%n_es_amplitudes
      this%side = side
      this%wf => wf
      this%restart = restart
!
      length_state_guesses = input%get_n_state_guesses()
!
      if (length_state_guesses /= this%n_vectors) &
         call output%error_msg('when manually specifying starting guesses for excited states, &
                               &you need to specify a guess for all states requested.')
!
      call mem%alloc(occupied, this%n_vectors)
      call mem%alloc(virtual, this%n_vectors)
!
      call input%get_state_guesses(occupied, virtual, this%n_vectors)
!
      allocate(this%indices(this%n_vectors))
!
      do state = 1, this%n_vectors
!
         this%indices(state) = wf%n_v*(occupied(state) - 1) + virtual(state)
!
      end do 
!
      call mem%dealloc(occupied, this%n_vectors)
      call mem%dealloc(virtual, this%n_vectors)
!
   end function new_es_manual_start_vector_tool
!
!
end module es_manual_start_vector_tool_class
