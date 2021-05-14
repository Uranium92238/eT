!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2021 the authors of eT
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
!!
!!    Hartree-Fock geometry optimization engine class module
!!    Written by Eirik F. Kjønstad, 2019
!!
   use reference_engine_class, only: reference_engine
!
   use global_out,      only: output
   use global_in,       only: input
   use timings_class,   only: timings
   use hf_class,        only: hf
   use task_list_class, only: task_list
!
   use bfgs_geoopt_hf_class, only: bfgs_geoopt_hf
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
      engine%algorithm        = 'bfgs'
      engine%ao_density_guess = 'sad'
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
      implicit none
!
      class(hf_geoopt_engine) :: engine
      class(hf)               :: wf
!
      type(bfgs_geoopt_hf) :: bfgs_geoopt
!
      if (wf%ao%has_ghost_atoms()) &
         call output%warning_msg("Ghosts are experimental in geometry optimization.")
!
      if (wf%embedded) &
         call output%error_msg('geometry optimization with embedding is not supported')
!
      if (.not. engine%restart .and. (trim(engine%ao_density_guess) == 'sad')) then
!
!        Generate SAD if requested
!
         call engine%generate_sad_density(wf)
!
      endif
!
      if (trim(engine%algorithm) == 'bfgs') then
!
         call engine%tasks%print_('optimize geometry')
!
         bfgs_geoopt = bfgs_geoopt_hf(engine%restart)
         call bfgs_geoopt%run(wf)
         call bfgs_geoopt%cleanup()
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
      class(hf_geoopt_engine) :: engine
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
      class(hf_geoopt_engine) :: engine
!
      engine%name_       = 'Hartree-Fock geometry optimization engine'
!
      engine%description = 'Calculates the optimized geometry for the Hartree-Fock wavefunction.'
      engine%tag         = 'geometry optimization'
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
