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
module cc_es_start_vector_factory_class
!
!!
!! CC es start vector factory class
!! Written by Alexander C. Paul, May 2022
!!
!
   use kinds
   use global_out, only: output
   use global_in,  only: input
   use ccs_class,  only: ccs
!
   implicit none
!
   type:: cc_es_start_vector_factory
!
      character(len=:), allocatable :: es_type, transformation
      logical :: restart
!
   contains
!
      procedure, public :: create => create_cc_es_start_vector_factory
!
   end type cc_es_start_vector_factory
!
!
   interface cc_es_start_vector_factory
!
      procedure :: new_cc_es_start_vector_factory
!
   end interface
!
!
contains
!
!
   function new_cc_es_start_vector_factory(es_type, transformation, restart) result(this)
!!
!!    New
!!    Written by Alexander C. Paul, May 2022
!!
      implicit none
!
      character(len=*), intent(in) :: es_type, transformation
      logical, intent(in) :: restart
!
      type(cc_es_start_vector_factory) :: this
!
      this%es_type = es_type
      this%transformation = transformation
      this%restart = restart
!
   end function new_cc_es_start_vector_factory
!
!
   subroutine create_cc_es_start_vector_factory(this, wf, tool)
!!
!!    Create
!!    Written by Sarai D. Folkestad and Alexander C. Paul, 2020-2022
!!
      use global_in, only: input
!
      use start_vector_tool_class, only: start_vector_tool
      use es_ip_start_vector_tool_class, only: es_ip_start_vector_tool
      use es_cvs_start_vector_tool_class, only: es_cvs_start_vector_tool
      use es_manual_start_vector_tool_class, only: es_manual_start_vector_tool
      use es_valence_start_vector_tool_class, only: es_valence_start_vector_tool
!
      implicit none
!
      class(cc_es_start_vector_factory), intent(in)  :: this
      class(ccs), intent(in) :: wf
!
      class(start_vector_tool), allocatable :: tool
!
      if (input%is_keyword_present('state guesses', 'solver cc es')) then
!
         tool = es_manual_start_vector_tool(wf, &
                                            this%transformation, &
                                            this%restart)
!
         return
!
      end if
!
!
      select case(this%es_type)
!
         case ('valence')
!
            tool = es_valence_start_vector_tool(wf, &
                                                this%transformation, &
                                                this%restart)
!
         case ('core')
!
            tool = es_cvs_start_vector_tool(wf, &
                                            this%transformation, &
                                            this%restart)
         case ('ionize')
!
            tool = es_ip_start_vector_tool(wf, &
                                           this%transformation, &
                                           this%restart)
!
         case ('remove core')
!
            tool = es_valence_start_vector_tool(wf, &
                                                this%transformation, &
                                                this%restart)
!
         case default
!
            call output%error_msg('could not recognize excited state type in cc_es_start_vector_factory')
!
      end select
!
   end subroutine create_cc_es_start_vector_factory
!
!
end module cc_es_start_vector_factory_class
