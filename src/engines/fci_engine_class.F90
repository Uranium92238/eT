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
!!
!! Full CI engine class module
!! Written by Enrico Ronca, 2020
!!
   use kinds
   use global_in,       only: input
   use global_out,      only: output
   use timings_class,   only: timings
   use task_list_class, only: task_list
!
   use abstract_engine_class, only: abstract_engine
   use fci_class,       only: fci
!
   type, extends(abstract_engine) :: fci_engine
!
      character(len=200) :: fci_algorithm
!
      logical :: fci_restart
!
   contains
!
      procedure :: ignite           => ignite_fci_engine
      procedure :: run              => run_fci_engine
      procedure :: read_settings    => read_settings_fci_engine
      procedure :: set_printables   => set_printables_fci_engine
!
   end type fci_engine
!
!
   interface fci_engine
!
      procedure :: new_fci_engine
!
   end interface fci_engine
!
!
contains
!
!
   function new_fci_engine() result(engine)
!!
!!    New FCI engine
!!    Written by Enrico Ronca, 2020
!!
      implicit none
!
      type(fci_engine) :: engine
!
!     Set defaults and then read if nonstandard
!
      engine%fci_algorithm = 'davidson'
      engine%fci_restart   = .false.
!
      call engine%read_settings()
!
      call engine%set_printables()
!
      engine%timer = timings(trim(engine%name_))
      call engine%timer%turn_on()
!
   end function new_fci_engine
!
!
   subroutine ignite_fci_engine(engine, wf)
!!
!!    Ignite
!!    Written by Enrico Ronca, 2020
!!
      implicit none
!
      class(fci_engine) :: engine
!
      class(fci) :: wf
!
      call engine%print_banner(wf)
      call engine%run(wf)
      call engine%print_timings(wf)
!
   end subroutine ignite_fci_engine
!
!
   subroutine read_settings_fci_engine(engine)
!!
!!    Read fci settings
!!    Written by Enrico Ronca, 2020
!!
      implicit none
!
      class(fci_engine) :: engine
!
      engine%fci_restart = input%is_keyword_present('restart', 'solver fci')
!
   end subroutine read_settings_fci_engine
!
!
   subroutine run_fci_engine(engine, wf)
!!
!!    Run
!!    Written by Enrico Ronca, 2021
!!
      use eigen_davidson_solver_class, only: eigen_davidson_solver
      use fci_solver_factory_class,   only: fci_solver_factory
      use warning_suppressor, only: do_nothing
!
      implicit none
!
      class(fci_engine)  :: engine
      class(fci)        :: wf
!
      class(eigen_davidson_solver), allocatable :: solver
!
      type(fci_solver_factory) :: solver_factory
!
      call do_nothing(engine)
!
      call solver_factory%create(wf, solver)
!
      call solver%run()
!
      call wf%print_fci_summary()
!
   end subroutine run_fci_engine
!
!
   subroutine set_printables_fci_engine(engine)
!!
!!    Set printables
!!    Written by Enrico Ronca, 2020
!!
      use string_utilities, only: convert_to_uppercase
!
      implicit none
!
      class(fci_engine) :: engine
!
      engine%name_  = 'Full CI engine'
!
      engine%tag = 'fci'
!
!     Prepare the list of tasks
!
      engine%tasks = task_list()
!
      call engine%tasks%add(label='fci solver',                              &
                        description='Diagonalization of the Hamiltonian ('//  &
                        trim((engine%fci_algorithm))//' algorithm)')
!
      engine%description  = 'Calculates the FCI states and energies'
!
   end subroutine set_printables_fci_engine
!
!
end module fci_engine_class
