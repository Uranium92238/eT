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
module cc_es_projector_factory_class
!
!!
!! CC es projector factory class
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
   type:: cc_es_projector_factory
!
      character(len=:), allocatable :: es_type
!
   contains
!
      procedure, public :: create => create_cc_es_projector_factory
!
   end type cc_es_projector_factory
!
   interface cc_es_projector_factory
!
      procedure :: new_cc_es_projector_factory
!
   end interface
!
contains
!
!
   function new_cc_es_projector_factory(es_type) result(this)
!!
!!    New
!!    Written by Alexander C. Paul, May 2022
!!
      implicit none
!
      character(len=*), intent(in) :: es_type
      type(cc_es_projector_factory) :: this
!
      this%es_type = es_type
!
   end function new_cc_es_projector_factory
!
!
   subroutine create_cc_es_projector_factory(this, wf, projector)
!!
!!    Create
!!    Written by Eirik F. Kjønstad, Sep 2019
!!
      use abstract_projection_tool_class, only: abstract_projection_tool
      use null_projection_tool_class, only: null_projection_tool
      use es_cvs_projection_tool_class, only: es_cvs_projection_tool
      use es_ip_projection_tool_class, only: es_ip_projection_tool
      use es_rm_core_projection_tool_class, only: es_rm_core_projection_tool
!
      implicit none
!
      class(cc_es_projector_factory), intent(in) :: this
      class(ccs), intent(in) :: wf
!
      class(abstract_projection_tool), allocatable, intent(out) :: projector
!
      select case (this%es_type)
!
         case ('valence')
!
            projector = null_projection_tool(wf%n_es_amplitudes)
!
         case ('core')
!
            projector = es_cvs_projection_tool(wf)
!
         case ('ionize')
!
            projector = es_ip_projection_tool(wf)
!
         case ('remove core')
!
            projector = es_rm_core_projection_tool(wf)
!
         case default
!
            call output%error_msg('could not recognize excited state type in cc_es_projector_factory')
!
      end select
!
   end subroutine create_cc_es_projector_factory
!
!
end module cc_es_projector_factory_class
