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
module es_valence_start_vector_tool_class
!
!!
!!    Excited state valence start vector tool class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018-2019
!!
!
   use ccs_class, only: ccs 
!
   use array_utilities, only: get_n_lowest
   use es_start_vector_tool_class
   use memory_manager_class, only: mem 
   use global_out, only: output 
!
   implicit none
!
   type, extends(es_start_vector_tool) :: es_valence_start_vector_tool
!
   contains
!
   end type es_valence_start_vector_tool
!
!
   interface es_valence_start_vector_tool
!
      procedure :: new_es_valence_start_vector_tool
!
   end interface es_valence_start_vector_tool
!
!
contains 
!
!
   function new_es_valence_start_vector_tool(wf, side, restart) result(this)
!!
!!    New ES valence start vector tool 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018-2019
!!
      implicit none 
!
      class(ccs), intent(in), target :: wf
      character(len=*), intent(in) :: side
      logical, intent(in) :: restart
!
      type(es_valence_start_vector_tool) :: this
!
      real(dp), dimension(:), allocatable :: eps         ! Orbital differences 
      real(dp), dimension(:), allocatable :: lowest_eps  ! Lowest orbital differences 
!
      if (wf%bath_orbital) &
         call output%error_msg('Bath orbitals can not be used in valence excitation calculation')
!
      this%n_vectors    = wf%n_singlet_states
      this%n_parameters = wf%n_es_amplitudes
      this%side         = side
      this%wf           => wf
      this%restart      = restart
!
      call mem%alloc(eps, this%n_parameters)
      call wf%get_orbital_differences(eps, this%n_parameters)
!
      call mem%alloc(lowest_eps, this%n_vectors)
!
      allocate(this%indices(this%n_vectors))
      call get_n_lowest(this%n_vectors, this%n_parameters, eps, lowest_eps, this%indices)
!
      call mem%dealloc(lowest_eps, this%n_vectors)
      call mem%dealloc(eps, this%n_parameters)
!
   end function new_es_valence_start_vector_tool
!
!
end module es_valence_start_vector_tool_class
