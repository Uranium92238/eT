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
!  along with engine program. If not, see <https://www.gnu.org/licenses/>.
!
!
module tdhf_engine_class
!
!!
!! Time dependent Hartree-Fock engine class module
!! Written by Sarai D. Folkestad, May 2021
!!
!
   use kinds
!
   use reference_engine_class, only: reference_engine
!
   use global_in,              only: input
   use global_out,             only: output
   use timings_class,          only: timings
   use memory_manager_class,   only: mem
   use task_list_class,        only: task_list
!
   use hf_class,               only: hf
!
   type, extends(reference_engine) :: tdhf_engine
!
   contains
!
      procedure :: run &
                => run_tdhf_engine
!
      procedure :: set_printables &
                => set_printables_tdhf_engine
!
      procedure, private :: do_tdhf
!
   end type tdhf_engine
!
   interface tdhf_engine
!
      procedure :: new_tdhf_engine
!
   end interface tdhf_engine
!
contains
!
   function new_tdhf_engine() result(engine)
!!
!!    New TDHF engine
!!    Written by Sarai D. Folkestad, May 2021
!!
      implicit none
!
      type(tdhf_engine) :: engine
!
      engine%name_       = 'Time dependent Hartree-Fock engine'
!
      engine%description = 'Drives the calculation of the Hartree-Fock &
                           &excitation energies and properties. '
!
      engine%tag         = 'properties'
!
      engine%ao_density_guess = 'sad'
      engine%algorithm        = 'scf-diis'
      engine%restart          = .false.
      engine%dipole           = .false.
      engine%quadrupole       = .false.
      engine%plot_orbitals    = .false.
      engine%plot_density     = .false.
      engine%write_mo_info    = .false.
      engine%molden_file      = .false.
      engine%skip_scf         = .false.
!
      call engine%read_settings()
!
   end function new_tdhf_engine
!
!
   subroutine run_tdhf_engine(engine, wf)
!!
!!    Run
!!    Written by Sarai D. Folkestad, May 2021
!!
      implicit none
!
      class(tdhf_engine), intent(in)    :: engine
      class(hf)         , intent(inout) :: wf
!
      call engine%reference_engine%run(wf)
!
      call engine%do_tdhf(wf)
!
   end subroutine run_tdhf_engine
!
!
   subroutine do_tdhf(engine, wf)
!!
!!    Do TDHF
!!    Written by Sarai D. Folkestad, May 2021
!!
      use eigen_davidson_solver_class, only: eigen_davidson_solver
      use tdhf_solver_factory_class,   only: tdhf_solver_factory
!
      implicit none

      class(tdhf_engine), intent(in) :: engine
      class(hf), intent(inout)       :: wf
!
      class(eigen_davidson_solver), allocatable :: solver
      type(tdhf_solver_factory) :: solver_factory
!
      call engine%tasks%print_('excitation energies')
!
      call solver_factory%create(wf, solver)
!
      call wf%initialize_tdhf_quantities(solver%get_n_solutions())
!
      call solver%run()
      call wf%tdhf_summary()
!
   end subroutine do_tdhf
!
!
   subroutine set_printables_tdhf_engine(engine)
!!
!!    Set Printables
!!    Written by Sarai D. Folkestad, May 2021
!!
!!    Should be overwritten by descendants.
!
      use string_utilities, only : convert_to_uppercase
!
      implicit none
!
      class(tdhf_engine), intent(inout) :: engine
!
      call engine%reference_engine%set_printables()
!
      call engine%tasks%add(label='excitation energies', &
            description='Calculate TDHF excitation energies')
!
   end subroutine set_printables_tdhf_engine
!
!
end module tdhf_engine_class
