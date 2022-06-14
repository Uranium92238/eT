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
module fci_engine_class
!
!!
!! FCI engine class
!! Written by Sarai D. Folkestad, Eirik F. Kjønstad, and Alexander C. Paul, 2018-2022
!!
!
   use abstract_fci_engine_class, only: abstract_fci_engine
   use fci_class, only: fci
!
   use fci_eigenproblem_task_class, only: fci_eigenproblem_task
   use fci_mean_value_task_class, only: fci_mean_value_task
!
   implicit none
!
   type, extends(abstract_fci_engine) :: fci_engine
!
      type(fci_eigenproblem_task), allocatable, private :: eigenproblem
      type(fci_mean_value_task), allocatable, private :: mean_value
!
   contains
!
      procedure, public :: ignite => ignite_fci_engine
!
   end type fci_engine
!
contains
!
!
   subroutine ignite_fci_engine(this, wf)
!!
!!    Ignite
!!    Written by Alexander C. Paul, May 2022
!!
      use global_in, only: input
!
      implicit none
!
      class(fci_engine), intent(inout) :: this
      class(fci), intent(inout) :: wf
!
      this%eigenproblem = fci_eigenproblem_task()
      call this%eigenproblem%execute(wf)
!
      if (input%is_section_present('fci mean value')) then
!
         this%mean_value = fci_mean_value_task()
         call this%mean_value%execute(wf)
!
      end if
!
   end subroutine ignite_fci_engine
!
!
   end module fci_engine_class
   