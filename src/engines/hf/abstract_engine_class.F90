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
module abstract_engine_class
!
!!
!! Abstract engine class module
!! Written by Tor S. Haugland, 2019
!!
!! Each descendant must make an ignite-routine with the corresponding wavefunction.
!!
!
   use kinds
   use parameters
!
   use global_out,           only: output
   use timings_class,        only: timings
   use string_utilities,     only: convert_to_uppercase
   use memory_manager_class, only: mem
   use task_list_class,      only: task_list
!
   use wavefunction_class, only: wavefunction
!
   implicit none
!
   type, abstract :: abstract_engine
!
      character(len=200) :: name_
      character(len=200) :: tag
      character(len=200) :: description
!
      logical :: dipole
      logical :: quadrupole
      logical :: plot_density
!
      type(timings) :: timer ! Timer for engine. Obs! must be turned on in constructor
!
      type(task_list), allocatable :: tasks  ! The printed tasks of the engine.
                                             ! Should be set in constructor
!
   contains
!
      procedure :: print_banner           => print_banner_abstract_engine
      procedure :: print_timings          => print_timings_abstract_engine
!
   end type abstract_engine
!
!
contains
!
!
   subroutine print_banner_abstract_engine(engine, wf)
!!
!!    Print banner
!!    Written by Sarai D. Folkestad, May 2019
!!
!!    Prints:
!!
!!       - Engine name
!!       - Authors and date
!!       - Wavefunction type
!!       - Engine tasks
!!
!!    Dependencies:
!!
!!       - The printables of the engine must be set for each decendant (constructor)
!!
      use timings_class,          only: timing
!
      implicit none
!
      class(abstract_engine), intent(in) :: engine
      class(wavefunction),    intent(in) :: wf
!
      if (.not. allocated(engine%tasks)) then
         call output%error_msg('Tasks of engine was not set. Do this in prepare.')
      end if
!
      call output%printf('m', ":: (a0)", fs='(//t3,a)', chars=[engine%name_])
      call output%print_separator('minimal', len(trim(engine%name_))+6, '=')
!
      call timing%printf('m', ":: (a0)", fs='(/t3,a)', chars=[engine%name_])
      call timing%print_separator('minimal', len(trim(engine%name_))+6, '=')
!
      call output%printf('m', "(a0)", fs='(/t3,a)', chars=[engine%description])
!
      call output%printf('m', 'This is a (a0) (a0) calculation.', &
                         chars=[character(len=500)::convert_to_uppercase(wf%name_), &
                         engine%tag], ffs='(/t3,a)')
!
      call output%printf('m', 'The following tasks will be performed:', fs='(t3,a/)')
!
      call engine%tasks%print_all()
!
   end subroutine print_banner_abstract_engine
!
!
   subroutine print_timings_abstract_engine(engine, wf)
!!
!!    Print timings
!!    Written by Sarai D. Folkestad, May 2019
!!
!!    Prints the timings of the engine.
!!
      implicit none
!
      class(abstract_engine), intent(inout) :: engine
      class(wavefunction),    intent(in)    :: wf
!
      call output%printf('n', '- Timings for the (a0) (a0) calculation', &
                         chars=[character(len=500)::convert_to_uppercase(wf%name_), &
                         engine%tag], fs='(/t3, a)')
!
      call engine%timer%turn_off()
!
      call output%printf('n', 'Total wall time (sec): (f20.5)', &
                         reals=[engine%timer%get_elapsed_time('wall')], fs='(/t6, a)')
      call output%printf('n', 'Total cpu time (sec):  (f20.5)', &
                         reals=[engine%timer%get_elapsed_time('cpu')], fs='(t6, a)')
!
   end subroutine print_timings_abstract_engine
!
!
end module abstract_engine_class
