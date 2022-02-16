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
module cc_response_engine_factory_class
!
!!
!! CC response engine factory class
!! Written by Alexander C. Paul, Feb 2022
!!
!
   use global_in,  only: input
   use global_out, only: output
!
   implicit none
!
!
   type :: cc_response_engine_factory
!
   contains
!
      procedure, public :: create
!
      procedure, private, nopass :: create_transition_moment_engine
!
end type cc_response_engine_factory
!
!
contains
!
!
   subroutine create(this, response_engine)
!!
!!    Create
!!    Written by Alexander C. Paul, Feb 2022
!!
      use cc_engine_class, only: cc_engine
      use cc_polarizability_engine_class, only: cc_polarizability_engine
      use eom_transition_moment_engine_class, only: eom_transition_moment_engine
!
      implicit none
!
      class(cc_response_engine_factory), intent(inout) :: this
!
      class(cc_engine), allocatable :: response_engine
!
      if (input%is_keyword_present('polarizabilities', 'cc response')) then
!
         allocate(cc_polarizability_engine :: response_engine)
!
      else if (input%is_keyword_present('transition moments', 'cc response')) then
!
         call this%create_transition_moment_engine(response_engine)
!
      else if (input%is_keyword_present('permanent moments', 'cc response')) then
!
         allocate(eom_transition_moment_engine :: response_engine)
!
      else
!
         call output%error_msg("Either (a0), (a0), or (a0) need to be &
                               &specified for response calculation.", &
                               chars=[character(len=18) ::  "polarizabilities", &
                                                            "transition moments", &
                                                            "permanent moments"])
!
      end if
!
   end subroutine create
!
!
   subroutine create_transition_moment_engine(response_engine)
!!
!!    Create transition moment engine
!!    Written by Alexander C. Paul, Feb 2022
!!
      use cc_engine_class, only: cc_engine
      use eom_transition_moment_engine_class, only: eom_transition_moment_engine
      use lr_transition_moment_engine_class, only: lr_transition_moment_engine
!
      implicit none
!
      class(cc_engine), allocatable :: response_engine
!
      if (input%is_keyword_present('lr', 'cc response')) then
!
         allocate(lr_transition_moment_engine :: response_engine)
!
      else if (input%is_keyword_present('eom', 'cc response')) then
!
         allocate(eom_transition_moment_engine :: response_engine)
!
      else
!
         call output%error_msg("Either eom or lr need to be specified for response.")
!
      end if
!
   end subroutine create_transition_moment_engine
!
!
end module cc_response_engine_factory_class
