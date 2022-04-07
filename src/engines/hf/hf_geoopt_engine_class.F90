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
module hf_geoopt_engine_class
!
!!
!! Hartree-Fock geometry optimization engine class module
!! Written by Eirik F. Kjønstad, 2019
!!
!
   use reference_engine_class, only: reference_engine
!
   use global_out,      only: output
   use global_in,       only: input
   use timings_class,   only: timings
   use hf_class,        only: hf
   use task_list_class, only: task_list
!
   type, extends(reference_engine) :: hf_geoopt_engine
!
   contains
!
      procedure :: run              => run_hf_geoopt_engine
!
      procedure :: read_settings    => read_settings_hf_geoopt_engine
      procedure :: set_printables   => set_printables_hf_geoopt_engine
!
   end type hf_geoopt_engine
!
!
   interface hf_geoopt_engine
!
      procedure :: new_hf_geoopt_engine
!
   end interface hf_geoopt_engine
!
!
contains
!
!
   function new_hf_geoopt_engine() result(engine)
!!
!!    New HF engine
!!    Written by Eirik F. Kjønstad, 2019
!!
      implicit none
!
      type(hf_geoopt_engine) :: engine
!
      engine%name_       = 'Hartree-Fock geometry optimization engine'
!
      engine%description = 'Calculates the optimized geometry for the Hartree-Fock wavefunction.'
      engine%tag         = 'geometry optimization'
!
      engine%algorithm        = 'bfgs'
      engine%restart          = .false.
      engine%plot_orbitals    = .false.
      engine%plot_density     = .false.
      engine%write_mo_info    = .false.
!
      call engine%read_settings()
!
      call engine%set_printables()
!
      engine%timer = timings(trim(engine%name_))
      call engine%timer%turn_on()
!
   end function new_hf_geoopt_engine
!
!
   subroutine run_hf_geoopt_engine(engine, wf)
!!
!!    Run
!!    Written by Eirik F. Kjønstad, 2019
!!
      use kinds
      use hf_energy_function_class, only: hf_energy_function
      use redundant_internal_coords_class, only: redundant_internal_coords
      use bfgs_solver_class, only: bfgs_solver
!
      implicit none
!
      class(hf_geoopt_engine), intent(in)    :: engine
      class(hf),               intent(inout) :: wf
!
      type(hf_energy_function), allocatable :: energy_function
      type(redundant_internal_coords), allocatable :: internals
!
      class(bfgs_solver), allocatable :: solver
!
      integer :: max_iterations
      real(dp) :: max_step, energy_threshold, gradient_threshold
!
      if (wf%ao%has_ghost_atoms()) &
         call output%warning_msg("Ghosts are experimental in geometry optimization.")
!
      if (wf%embedded) &
         call output%error_msg('geometry optimization with embedding is not supported')
!
      if (trim(engine%algorithm) == 'bfgs') then
!
         call engine%tasks%print_('optimize geometry')
!
         max_iterations = 100
         call input%get_keyword('max iterations', 'solver scf geoopt', max_iterations)
!
         max_step = 0.5d0
         call input%get_keyword('max step', 'solver scf geoopt', max_step)
!
         energy_threshold = 1.0d-6
         call input%get_keyword('energy threshold', 'solver scf geoopt', energy_threshold)
!
         gradient_threshold = 3.0d-4
         call input%get_keyword('gradient threshold', 'solver scf geoopt', gradient_threshold)
!
         energy_function = hf_energy_function(wf)
!
         internals = redundant_internal_coords(3*wf%n_atomic_centers)
!
         call wf%initialize_redundant_internal_coordinates(internals)
!
         solver = bfgs_solver(energy_function,     &
                              internals,           &
                              max_iterations,      &
                              max_step,            &
                              energy_threshold,    &
                              gradient_threshold)
!
         call solver%initialize()
         call solver%run()
         call solver%finalize()
!
         deallocate(internals)
!
      else
!
         call output%error_msg('did not recognize hf geoopt algorithm: '// engine%algorithm)
!
      endif
!
   end subroutine run_hf_geoopt_engine
!
!
   subroutine read_settings_hf_geoopt_engine(engine)
!!
!!    Read settings
!!    Written by Eirik F. Kjønstad, 2019
!!
      implicit none
!
      class(hf_geoopt_engine), intent(inout) :: engine
!
      call input%get_keyword('algorithm', 'solver scf geoopt', engine%algorithm)
!
      if (input%is_keyword_present('restart', 'solver scf geoopt')) engine%restart = .true.
      if (input%is_keyword_present('restart', 'do')) engine%restart = .true.
!
      engine%write_mo_info = input%is_keyword_present('print orbitals', 'solver scf')
!
   end subroutine read_settings_hf_geoopt_engine
!
!
   subroutine set_printables_hf_geoopt_engine(engine)
!!
!!    Set Printables
!!    Written by Sarai D. Folkestad, May 2019
!!
!!    Should be overwritten by descendants.
!!
!
      use string_utilities, only : convert_to_uppercase
!
      implicit none
!
      class(hf_geoopt_engine), intent(inout) :: engine
!
      engine%tasks = task_list()
!
      if (trim(engine%ao_density_guess) == 'sad' .and. .not. engine%restart) &
         call engine%tasks%add(label='sad', description='Generate initial SAD density')
!
      call engine%tasks%add(label='optimize geometry', &
                            description='Calculation of optimized geometry (' //&
                            trim((engine%algorithm)) // ' algorithm)')
!
   end subroutine set_printables_hf_geoopt_engine
!
!
end module hf_geoopt_engine_class
