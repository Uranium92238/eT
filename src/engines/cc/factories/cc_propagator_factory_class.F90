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
module cc_propagator_factory_class
!
!!
!! CC propagator factory class
!! Written by Eirik F. Kjønstad, 2021
!!
!
   use ccs_class,  only: ccs
   use global_in,  only: input
   use global_out, only: output
!
   implicit none
!
!
   type :: cc_propagator_factory
!
      character(len=200), private :: integrator
!
   contains
!
      procedure, public :: create
!
   end type cc_propagator_factory
!
!
   interface cc_propagator_factory
!
      procedure :: new_cc_propagator_factory
!
   end interface cc_propagator_factory
!
!
contains
!
!
   function new_cc_propagator_factory() result(this)
!!
!!    New CC propagator factory
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none
!
      type(cc_propagator_factory) :: this
!
      this%integrator = 'rk4'
      call input%get_keyword('integrator', 'solver cc propagation', this%integrator)
!
   end function new_cc_propagator_factory
!
!
   subroutine create(this, wf, propagator)
!!
!!    Create
!!    Written by Eirik F. Kjønstad, 2021
!!
      use cc_propagation_class,           only: cc_propagation
!
      use euler_cc_propagation_class,     only: euler_cc_propagation
      use rk4_cc_propagation_class,       only: rk4_cc_propagation
      use gl2_cc_propagation_class,       only: gl2_cc_propagation
      use gl4_cc_propagation_class,       only: gl4_cc_propagation
      use gl6_cc_propagation_class,       only: gl6_cc_propagation
!
      implicit none
!
      class(cc_propagator_factory), intent(inout) :: this
!
      class(ccs), intent(in) :: wf
!
      class(cc_propagation), allocatable :: propagator
!
      if (this%integrator == 'euler') then
!
         propagator = euler_cc_propagation(wf)
!
      elseif (this%integrator == 'rk4') then
!
         propagator = rk4_cc_propagation(wf)
!
      elseif (this%integrator == 'gl2') then
!
         propagator = gl2_cc_propagation(wf)
!
      elseif (this%integrator == 'gl4') then
!
         propagator = gl4_cc_propagation(wf)
!
      elseif (this%integrator == 'gl6') then
!
         propagator = gl6_cc_propagation(wf)
!
      else
!
         call output%error_msg('Integrator for CC propagation not recognized!')
!
      endif
!
   end subroutine create
!
!
end module cc_propagator_factory_class
