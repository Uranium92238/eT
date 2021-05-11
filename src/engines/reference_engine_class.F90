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
module reference_engine_class
!!
!!    Hartree-Fock engine class module 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
   use kinds
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
      logical :: print_mo_info
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
      procedure, nopass :: calculate_quadrupole_moment => calculate_quadrupole_moment_reference_engine
      procedure, nopass :: calculate_dipole_moment     => calculate_dipole_moment_reference_engine
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
      procedure, private :: check_algorithm            => check_algorithm_reference_engine
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
      engine%ao_density_guess = 'sad'
      engine%algorithm        = 'scf-diis'
      engine%restart          = .false.
      engine%dipole           = .false.
      engine%quadrupole       = .false.
      engine%plot_orbitals    = .false.
      engine%plot_density     = .false.
      engine%print_mo_info    = .false.
      engine%skip_scf         = .false.
!  
      call engine%read_settings()
      call engine%check_algorithm()
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
      class(reference_engine) :: engine 
      class(hf)               :: wf 
!
!     Overwrite restart if the corresponding files don't exist
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
      class(reference_engine)    :: engine 
      class(hf)                  :: wf
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
      call wf%print_summary(engine%print_mo_info)
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
      class(reference_engine) :: engine 
!
      call input%get_keyword('algorithm', 'solver scf', engine%algorithm)
!
      engine%restart = input%is_keyword_present('restart', 'solver scf')
      engine%skip_scf = input%is_keyword_present('skip', 'solver scf')
!
      call input%get_keyword('ao density guess', 'solver scf', engine%ao_density_guess)
!
      if (input%is_keyword_present('print orbitals', 'solver scf')) then
         engine%print_mo_info = .true.
      end if
!
      if (input%is_keyword_present('plot hf orbitals', 'visualization')) then
         engine%plot_orbitals = .true.
      end if
!
      if (input%is_keyword_present('plot hf density', 'visualization')) then 
         engine%plot_density = .true.
      end if
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
      class(reference_engine) :: engine
!
      engine%name_       = 'Hartree-Fock engine'
!
      engine%description = 'Drives the calculation of the Hartree-Fock state. '
      engine%tag         = 'ground state'
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
!!    Written by Tor S. Haugland, Sep 2019
!!
!!    Generates SAD density for every unique (atom, basis) pair.
!!
!!       1: a molecular system with only one atom for every (atom, basis) is created.
!!       2: an UHF wavefunction is created for that system.
!!       3: the wavefunction is run through a SCF solver.
!!       4: the density files are moved to where the wavefunction expects them to be.
!!
!!    Modified by Tor S. Haugland, Nov 2019
!!
!!    Removed deletion of files generated in SCF for SAD to fix a re-occuring bug.
!!    Added "Found SAD". Only loop over unique atoms.
!!
      use sequential_file_class,  only: sequential_file
!
      use atomic_center_class,    only: atomic_center
!
      use uhf_class,              only: uhf
!
      use string_utilities,       only: index_of_unique_strings
!
      use timings_class,          only: timing
!
      implicit none
!
      class(reference_engine)          :: engine
!
      class(hf)                        :: wf
!
      type(uhf),         allocatable         :: sad_wf
      type(scf_solver), allocatable          :: sad_solver
!
      character(len=200)    :: ao_density_guess
      real(dp)              :: energy_threshold
      real(dp)              :: gradient_threshold
      integer               :: max_iterations
!
      character(len=200)    :: name_
      integer               :: multiplicity
!
      character(len=200)    :: alpha_fname
      character(len=200)    :: beta_fname
      type(sequential_file) :: alpha_density_file
      type(sequential_file) :: beta_density_file
!
      integer :: I
      character(len=50), dimension(wf%n_atomic_centers) :: atom_and_basis
      integer,           dimension(wf%n_atomic_centers) :: unique_atom_index
!
      type(timings), allocatable :: sad_generation_timer
!
      type(atomic_center), allocatable :: center 
!
      call engine%tasks%print_('sad')
!
      sad_generation_timer = timings('SAD generation time', pl='normal')
      call sad_generation_timer%turn_on()
!
!     SAD solver settings
!
      ao_density_guess   = 'core'
      max_iterations     = 100
!
      energy_threshold   = 1.0D-6
      call input%get_keyword('energy threshold', 'solver scf', energy_threshold)
      energy_threshold   = min(1.0D-6, energy_threshold)
!
      gradient_threshold = 1.0D-6
      call input%get_keyword('gradient threshold', 'solver scf', gradient_threshold)
      gradient_threshold = min(1.0D-6, gradient_threshold)
!
!     Find atomic index of unique atom/basis combinations
!
      allocate(center)
!
      do I = 1, wf%n_atomic_centers
!
         call wf%ao%get_center(I, center)
!
         atom_and_basis(I) = trim(center%symbol) // trim(center%basis)
!
      enddo
!
      call index_of_unique_strings(unique_atom_index, wf%n_atomic_centers, atom_and_basis)
!
!     For every unique atom, generate SAD density to file
!
      call timing%mute()
!
      do I = 1, wf%n_atomic_centers
!
         call wf%ao%get_center(I, center)
!
!        Check for unique atoms and ghosts
         if (all(unique_atom_index /= I) .or. center%is_ghost()) cycle
!
         name_ = "sad_" // trim(center%basis) // "_" // trim(center%symbol)
!
         alpha_fname = trim(name_) // '_alpha'
         beta_fname  = trim(name_) // '_beta'
!
!        Prepare molecule of the chosen atom
!
         call output%mute()
!
         multiplicity = center%get_ground_state_multiplicity()
!
!        Prepare SAD wavefunction
!
         sad_wf = uhf(fractional_uniform_valence=.true., &
                      multiplicity=multiplicity)
!
         call sad_wf%prepare([center],  embedding=.false., charge=0)
!
!        Prepare and run solver
!
         sad_solver = scf_solver(restart=.false.,                 &
                           ao_density_guess=ao_density_guess,     &
                           max_iterations=max_iterations,         &
                           gradient_threshold=gradient_threshold, &
                           acceleration_type='none',              &
                           skip = .false.)
!
         call sad_solver%run(sad_wf)
!
!        Cleanup and generate ao_density_a and ao_density_b
!
         call sad_wf%orbital_file%delete_() 
         call sad_wf%cleanup()
!
         deallocate(sad_wf)
!
         call output%unmute()
!
         call output%printf('v', 'Generated atomic density for ' //  &
                            adjustl(center%symbol) // ' using UHF/(a0)', &
                            chars=[center%basis], fs='(t6,a)')
!
!        Move densities to where "set_ao_density_sad" can use them,
!        but first delete SAD if it already exists.
!
         alpha_density_file = sequential_file(alpha_fname)
         if (alpha_density_file%exists()) call alpha_density_file%delete_()
!
         beta_density_file  = sequential_file(beta_fname)
         if (beta_density_file%exists())  call beta_density_file%delete_()
!
         alpha_density_file = sequential_file('ao_density_a')
         call alpha_density_file%copy(alpha_fname)
         call alpha_density_file%delete_()
!
         beta_density_file  = sequential_file('ao_density_b')
         call beta_density_file%copy(beta_fname)
         call beta_density_file%delete_()
!
      enddo
!
      deallocate(center)
!
!     Libint is overwritten by SAD. Re-initialize.
!
      call wf%ao%export_centers_to_libint()
!
      call timing%unmute()
      call sad_generation_timer%turn_off()
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
      timer = timings('Orbital plotting time', pl='normal')
      call timer%turn_on()
!
!     Read orbital plotting settings
!
      n_orbitals_to_plot = input%get_n_elements_for_keyword('plot hf orbitals', &
                                                                        'visualization')
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
         write(orbital_file_tags(i), '(i4.4)') orbitals_to_plot(i)
         orbital_file_tags(i) = 'MO_' // trim(orbital_file_tags(i))
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
      real(dp), dimension(3) :: mu_electronic
      real(dp), dimension(3) :: mu_nuclear
      real(dp), dimension(3) :: mu_total
!
      real(dp), dimension(6) :: q_electronic
      real(dp), dimension(6) :: q_nuclear
      real(dp), dimension(6) :: q_total
!      
      character(len=4), dimension(:), allocatable :: components
!
      type(timings), allocatable :: timer 
!
      timer = timings('Time to calculte dipole and/or quadrupole', pl='normal')
      call timer%turn_on()
!
      call engine%tasks%print_('expectation value')
!
      if(engine%dipole) then 
!
         call engine%calculate_dipole_moment(wf, mu_electronic, mu_nuclear, mu_total)
!
         allocate(components(3))
!
         components = (/'x   ',&
                        'y   ',&
                        'z   '/)
!
         call engine%print_operator('dipole moment', mu_electronic, mu_nuclear, mu_total, &
                                    components, 3)
!
         deallocate(components)
!
      endif
!
      if (engine%quadrupole) then
!
         call engine%calculate_quadrupole_moment(wf, q_electronic, q_nuclear, q_total)
!
         allocate(components(6))
!
         components = (/ 'xx  ',   &
                         'xy  ',   &
                         'xz  ',   &
                         'yy  ',   &
                         'yz  ',   &
                         'zz  '    /)
!
         call engine%print_operator('quadrupole moment (with trace)', q_electronic, q_nuclear, q_total, &
                                    components, 6)
!
         call engine%remove_trace(q_electronic)
         call engine%remove_trace(q_nuclear)
!
         q_total = q_electronic + q_nuclear
!
         call output%printf('m', 'The traceless quadrupole is calculated as:', fs='(/t6,a)')
         call output%printf('m', 'Q_ij = 1/2[3*q_ij - tr(q)*delta_ij]',fs='(/t9,a)')
         call output%printf('m', 'where q_ij is the non-traceless matrix',fs='(/t6,a)')
!
         call engine%print_operator('traceless quadrupole moment', q_electronic, q_nuclear, q_total, &
                                    components, 6)
!
         deallocate(components)
!
      endif
!
      call timer%turn_off()
!
   end subroutine calculate_mean_values_reference_engine
!
!
   subroutine calculate_dipole_moment_reference_engine(wf, electronic, nuclear, total)
!!
!!    Calculate dipole moment
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
!!    Modified by Linda Goletto, Anders Hutcheson 
!!    and Tommaso Giovannini, Oct 2019
!!
!!    Calculates tr(D mu) in the AO basis; if the wf is mlhf,
!!    it only calculates tr(D Q) for the active space
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(3), intent(out) :: electronic
      real(dp), dimension(3), intent(out) :: nuclear
      real(dp), dimension(3), intent(out) :: total
!
      integer :: k
!
      real(dp), dimension(:,:,:), allocatable :: mu_pqk
!
      if(wf%name_.eq.'mlhf') &
         call output%warning_msg('dipole moments are size-extensive and&
                                 & are not well defined in MLHF.')
!
!     Get the integrals mu_pqk for components k = 1, 2, 3
!
      call mem%alloc(mu_pqk, wf%ao%n, wf%ao%n, 3)
      call wf%ao%get_oei('dipole', mu_pqk)
!
!     Get electronic expectation value contribution
!
      do k = 1, 3
!
         electronic(k) = wf%calculate_expectation_value(mu_pqk(:,:,k), wf%ao_density)
!
      enddo
!
      call mem%dealloc(mu_pqk, wf%ao%n, wf%ao%n, 3)
!
!     Get nuclear expectation value contribution, then sum the two
!
      nuclear = wf%get_nuclear_dipole()
!
      total = electronic + nuclear
!
   end subroutine calculate_dipole_moment_reference_engine
!
!
   subroutine calculate_quadrupole_moment_reference_engine(wf, electronic, nuclear, total)
!!
!!    Calculate quadrupole moment
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
!!    Modified by Linda Goletto, Anders Hutcheson 
!!    and Tommaso Giovannini, Oct 2019
!!
!!    Calculates tr(D Q) in the AO basis; if the wf is mlhf,
!!    it only calculates tr(D Q) for the active space
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(6), intent(out) :: electronic
      real(dp), dimension(6), intent(out) :: nuclear
      real(dp), dimension(6), intent(out) :: total
!
      integer :: k
!
      real(dp), dimension(:,:,:), allocatable :: q_pqk
!
!     Get the integrals q_pqk for components k = 1, 2, ..., 6 in the T1-transformed basis
!
      call mem%alloc(q_pqk, wf%ao%n, wf%ao%n, 6)
!
      call wf%ao%get_oei('quadrupole', q_pqk)
!
!     Get electronic expectation value contribution
!
      do k = 1, 6
!
         electronic(k) = wf%calculate_expectation_value(q_pqk(:,:,k), wf%ao_density)
!
      enddo
!
      call mem%dealloc(q_pqk, wf%ao%n, wf%ao%n, 6)
!
!     Get nuclear expectation value contribution, then sum the two
!
      if(wf%name_.eq.'mlhf') &
         call output%warning_msg('quadrupole moments are size-extensive&
                                 &and are not well defined in MLHF.')
!
      nuclear = wf%get_nuclear_quadrupole()
!
      total = electronic + nuclear
!
   end subroutine calculate_quadrupole_moment_reference_engine
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
      implicit none
!
      class(reference_engine), intent(in)       :: engine 
      class(hf), intent(inout)                  :: wf 
!
      class(scf_solver),  allocatable           :: scf
      character(len=200)                        :: acceleration_type
!
      call engine%tasks%print_('gs solver')
!
      acceleration_type = 'none'
!
      if (trim(engine%algorithm) == 'scf-diis' .or. &
          trim(engine%algorithm) == 'mo-scf-diis') acceleration_type = 'diis'
!
      scf = scf_solver(engine%restart, acceleration_type, engine%skip_scf)     
      call scf%run(wf)
!
   end subroutine do_ground_state_reference_engine
!
!
   subroutine check_algorithm_reference_engine(engine)
!!
!!    Check algorithm 
!!    Written by Sarai D. Folkestad, 2020
!!
      implicit none

      class(reference_engine), intent(in) :: engine
!
      if (trim(engine%algorithm) .ne. 'scf-diis'    .and. &
          trim(engine%algorithm) .ne. 'scf'         .and. &
          trim(engine%algorithm) .ne. 'mo-scf-diis') then
!
         call output%error_msg('did not recognize SCF algorithm')
!
      endif
!
   end subroutine check_algorithm_reference_engine
!
!
end module reference_engine_class
