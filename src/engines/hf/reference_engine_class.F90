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
module reference_engine_class
!
!!
!! Hartree-Fock engine class module
!! Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!
   use parameters
!
   use abstract_engine_class, only: abstract_engine
!
   use global_in,                   only: input
   use global_out,                  only: output
   use timings_class,               only: timings
   use memory_manager_class,        only: mem
   use task_list_class,             only: task_list
!
   use hf_class,                    only: hf
!
   use scf_solver_class,            only: scf_solver
!
   type, extends(abstract_engine) :: reference_engine
!
      character(len=200) :: ao_density_guess
      character(len=200) :: algorithm
!
      logical :: restart
      logical :: requested_mean_value
!
      logical :: plot_orbitals
      logical :: write_mo_info, molden_file
!
      logical :: skip_scf
!
   contains
!
      procedure :: ignite                              => ignite_reference_engine
!
      procedure :: run                                 => run_reference_engine
      procedure :: read_settings                       => read_settings_reference_engine
      procedure :: read_mean_value_settings &
                    => read_mean_value_settings_reference_engine
!
      procedure :: calculate_mean_values               => calculate_mean_values_reference_engine
!
      procedure :: set_printables                      => set_printables_reference_engine
!
      procedure :: generate_sad_density                => generate_sad_density_reference_engine
!
      procedure :: do_visualization                    => do_visualization_reference_engine
      procedure, nopass :: do_orbital_plotting         => do_orbital_plotting_reference_engine
!
      procedure :: do_ground_state                     => do_ground_state_reference_engine
!
   end type reference_engine
!
!
   interface reference_engine
!
      procedure :: new_reference_engine
!
   end interface reference_engine
!
!
contains
!
!
   function new_reference_engine() result(engine)
!!
!!    New reference engine
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      type(reference_engine) :: engine
!
      engine%name_            = 'Hartree-Fock engine'
!
      engine%description      = 'Drives the calculation of the Hartree-Fock state. '
      engine%tag              = 'ground state'
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
   end function new_reference_engine
!
!
   subroutine ignite_reference_engine(engine, wf)
!!
!!    Ignite
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Apr 2019
!!
      implicit none
!
      class(reference_engine), intent(inout) :: engine
      class(hf),               intent(inout) :: wf
!
!     Overwrite restart if the corresponding files don't exist
!
      if (engine%restart) engine%restart = wf%is_restart_possible()
      call engine%set_printables()
!
      engine%timer = timings(trim(engine%name_))
      call engine%timer%turn_on()
!
      call engine%print_banner(wf)
      call engine%run(wf)
      call engine%print_timings(wf)
!
   end subroutine ignite_reference_engine
!
!
   subroutine run_reference_engine(engine, wf)
!!
!!    Run
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(reference_engine), intent(in)    :: engine
      class(hf),               intent(inout) :: wf
!
      if ((.not. engine%restart) .and.  &
          (.not. engine%skip_scf) .and. &
          (trim(engine%ao_density_guess) == 'sad')) then
!
!        Generate SAD if requested
!
         call engine%generate_sad_density(wf)
!
      endif
!
!     Solve equations
!
      call engine%do_ground_state(wf)
!
      if (.not. engine%skip_scf) call wf%flip_final_orbitals()
      call wf%print_summary(engine%write_mo_info)
      if (engine%molden_file) call wf%write_molden_file()
!
!     Plot orbitals and/or density
!
      if (engine%plot_orbitals .or. engine%plot_density) call engine%do_visualization(wf)
!
!     Calculate properties
!
      if(engine%requested_mean_value) call engine%calculate_mean_values(wf)
!
   end subroutine run_reference_engine
!
!
   subroutine read_settings_reference_engine(engine)
!!
!!    Read settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(reference_engine), intent(inout) :: engine
!
      call input%get_keyword('algorithm', 'solver scf', engine%algorithm)
!
      engine%restart = input%is_keyword_present('restart', 'solver scf')
      engine%skip_scf = input%is_keyword_present('skip', 'solver scf')
!
      call input%get_keyword('ao density guess', 'solver scf', engine%ao_density_guess)
!
      engine%write_mo_info = input%is_keyword_present('print orbitals', 'solver scf')
!
      engine%molden_file = input%is_keyword_present('write molden', 'solver scf')
!
      engine%plot_orbitals = input%is_keyword_present('plot hf orbitals', 'visualization')
!
      engine%plot_density = input%is_keyword_present('plot hf density', 'visualization')
!
!     Global restart
      if (input%is_keyword_present('restart', 'do')) then
         engine%restart = .true.
      end if
!
      call engine%read_mean_value_settings()
!
   end subroutine read_settings_reference_engine
!
!
   subroutine set_printables_reference_engine(engine)
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
      class(reference_engine), intent(inout) :: engine
!
!     Prepare the list of tasks
!
      engine%tasks = task_list()
!
      if (trim(engine%ao_density_guess) == 'sad' .and. .not. engine%restart) &
         call engine%tasks%add(label='sad', description='Generate initial SAD density')
!
      call engine%tasks%add(label='gs solver',                                &
                            description='Calculation of reference state (' // &
                                 trim(convert_to_uppercase(engine%algorithm)) // ' algorithm)')
!
      if (engine%plot_orbitals .or. engine%plot_density) &
         call engine%tasks%add(label='plotting', description='Plot orbitals and/or density')
!
      if (engine%dipole .or. engine%quadrupole) &
         call engine%tasks%add(label='expectation value', &
            description='Calculate dipole and/or quadrupole moments')
!
   end subroutine set_printables_reference_engine
!
!
   subroutine generate_sad_density_reference_engine(engine, wf)
!!
!!    Generate SAD density
!!    Written by Sarai D. Folkestad, Oct 2021
!
      use sad_tool_class, only: sad_tool
      use atomic_center_class, only: atomic_center
!
      implicit none
!
      class(reference_engine)          :: engine
      class(hf)                        :: wf
!
      type(atomic_center), dimension(:), allocatable :: centers
!
      type(sad_tool) :: sad
!
      real(dp) :: gradient_threshold
!
      call engine%tasks%print_('sad')
!
      gradient_threshold = 1.0D-6
      call input%get_keyword('gradient threshold', 'solver scf', gradient_threshold)
!
      sad = sad_tool(gradient_threshold)
!
      allocate(centers(wf%ao%get_n_centers()))
      call wf%ao%get_centers(centers, 1, wf%ao%get_n_centers())
!
      call sad%generate(wf%ao%get_n_centers(), centers)
!
      deallocate(centers)
!
!     Libint is overwritten by SAD. Re-initialize.
      call wf%ao%export_centers_to_libint()
!
!     Re-determine status of a file because SAD may have deleted it
!     (so the status must go from "old" -> "new")
      call wf%orbital_file%determine_status()
!
   end subroutine generate_sad_density_reference_engine
!
!
   subroutine do_visualization_reference_engine(engine, wf)
!!
!!    Do visualization
!!    Written by Sarai D. Folkestad, Oct 2019
!!
!!    Writes orbitals and/or density to .plt files
!!    which may be opened in Chimera, if requested
!!    on input.
!!
!
      use visualization_class, only : visualization
!
      implicit none
!
      class(reference_engine) :: engine
      class(hf) :: wf
!
      type(visualization), allocatable :: plotter
!
      character(len=200) :: density_file_tag
!
      type(timings), allocatable :: density_plotting_timer
!
      call engine%tasks%print_('plotting')
!
      if (trim(wf%name_) .eq. 'uhf') call output%error_msg('no plotting for UHF')
!
!     Initialize the plotter
!
      plotter = visualization(wf%ao)
      call plotter%initialize(wf%ao)
!
      if (engine%plot_orbitals) then
!
         call engine%do_orbital_plotting(plotter, wf)
!
      endif
!
      if (engine%plot_density) then
!
         density_plotting_timer = timings('Density plotting time', pl='normal')
         call density_plotting_timer%turn_on()
!
         density_file_tag = 'AO_density'
         call plotter%plot_density(wf%ao, wf%ao_density, density_file_tag)
!
         call density_plotting_timer%turn_off()
!
      endif
!
      call plotter%cleanup()
!
   end subroutine do_visualization_reference_engine
!
!
   subroutine do_orbital_plotting_reference_engine(plotter, wf)
!!
!!    Do orbital plotting
!!    Written by Sarai D. Folkestad, Oct 2019
!!
!!    Reads orbitals to plot, and extracts the
!!    corresponding MO coefficients. Prepares
!!    file tags for the .plt files and writes
!!    orbital plot files using the visualization
!!    tool.
!!
      use visualization_class, only : visualization
      use memory_manager_class, only : mem
!
      implicit none
!
      class(hf) :: wf
!
      type(visualization) :: plotter
!
      integer :: i, n_orbitals_to_plot
!
      integer, dimension(:), allocatable :: orbitals_to_plot
!
      real(dp), dimension(:,:), allocatable :: orbital_coefficients
!
      character(len=200), dimension(:), allocatable :: orbital_file_tags
!
      type(timings), allocatable :: timer
!
      timer = timings('Plotting orbitals', pl='normal')
      call timer%turn_on()
!
!     Read orbital plotting settings
!
      n_orbitals_to_plot = input%get_n_elements_for_keyword('plot hf orbitals', 'visualization')
!
      call mem%alloc(orbitals_to_plot, n_orbitals_to_plot)
      call input%get_array_for_keyword('plot hf orbitals', 'visualization', &
            n_orbitals_to_plot, orbitals_to_plot)
!
!     Extract the orbitals we will plot from the full orbital coefficient matrix
!
      call mem%alloc(orbital_coefficients, wf%ao%n, n_orbitals_to_plot)
!
      do i = 1, n_orbitals_to_plot
!
         call dcopy(wf%ao%n, wf%orbital_coefficients(1, orbitals_to_plot(i)), &
               1, orbital_coefficients(1, i), 1)
!
      enddo
!
      allocate(orbital_file_tags(n_orbitals_to_plot))
!
!     Set file tags
!
      do i = 1, n_orbitals_to_plot
!
         write(orbital_file_tags(i), '(a,i4.4)') 'MO_', orbitals_to_plot(i)
!
      enddo
!
      call mem%dealloc(orbitals_to_plot, n_orbitals_to_plot)
!
!     Plot orbitals
!
      call plotter%plot_orbitals(wf%ao, orbital_coefficients, &
                                 n_orbitals_to_plot, orbital_file_tags)
!
      call mem%dealloc(orbital_coefficients, wf%ao%n, n_orbitals_to_plot)
      deallocate(orbital_file_tags)
!
      call timer%turn_off()
!
   end subroutine do_orbital_plotting_reference_engine
!
!
   subroutine read_mean_value_settings_reference_engine(engine)
!!
!!    Read mean value settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!    Modified by Tommaso Giovannini and Linda Goletto
!!
      implicit none
!
      class(reference_engine) :: engine
!
      engine%requested_mean_value = input%is_section_present('hf mean value')
!
      if (engine%requested_mean_value) then
!
         if (input%is_keyword_present('dipole','hf mean value')) &
             engine%dipole = .true.
!
         if (input%is_keyword_present('quadrupole','hf mean value')) &
             engine%quadrupole = .true.
!
      endif
!
   end subroutine read_mean_value_settings_reference_engine
!
!
   subroutine calculate_mean_values_reference_engine(engine, wf)
!!
!!    Calculate expectation values
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
      implicit none
!
      class(reference_engine), intent(in) :: engine
!
      class(hf), intent(in) :: wf
!
      type(timings), allocatable :: timer
!
      timer = timings('Time to calculte dipole and/or quadrupole', pl='normal')
      call timer%turn_on()
!
      call engine%tasks%print_('expectation value')
!
      if(engine%dipole) call wf%calculate_and_print_dipole()
!
      if (engine%quadrupole) call wf%calculate_and_print_quadrupole()
!
      call timer%turn_off()
!
   end subroutine calculate_mean_values_reference_engine
!
!
   subroutine do_ground_state_reference_engine(engine, wf)
!!
!!    Do ground state
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Constructs the solver specified on input.
!!    Solves the ground state.
!!
!
      use scf_solver_factory_class, only: scf_solver_factory
!
      implicit none
!
      class(reference_engine), intent(in)       :: engine
      class(hf), intent(inout)                  :: wf
!
      class(scf_solver),  allocatable           :: scf
!
      type(scf_solver_factory) :: factory
!
      call engine%tasks%print_('gs solver')
!
      call wf%prepare_for_scf(engine%restart, engine%skip_scf)
!
      factory = scf_solver_factory()
      call factory%create(wf, scf, engine%skip_scf)
!
      call scf%run(wf)
!
      call wf%finalize_gs()
!
   end subroutine do_ground_state_reference_engine
!
!
end module reference_engine_class
