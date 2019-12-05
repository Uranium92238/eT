!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2019 the authors of eT
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
   use global_in,            only: input
   use global_out,           only: output
   use timings_class,        only: timings
   use memory_manager_class, only: mem
   use task_list_class,      only: task_list
!
   use hf_class,          only: hf
!
   use scf_hf_class,      only: scf_hf
   use scf_diis_hf_class, only: scf_diis_hf
   use mo_scf_diis_class, only: mo_scf_diis
!
!
   type, extends(abstract_engine) :: reference_engine
!
      character(len=200) :: ao_density_guess
      character(len=200) :: algorithm 
!
      logical :: restart
      logical :: requested_zop
!
      logical :: plot_orbitals
!
   contains 
!
      procedure :: ignite                              => ignite_reference_engine
!
      procedure :: run                                 => run_reference_engine
      procedure :: read_settings                       => read_settings_reference_engine
      procedure :: read_zop_settings                   => read_zop_settings_reference_engine
!
      procedure :: calculate_expectation_values        => calculate_expectation_values_reference_engine
      procedure, nopass :: calculate_quadrupole_moment => calculate_quadrupole_moment_reference_engine
      procedure, nopass :: calculate_dipole_moment     => calculate_dipole_moment_reference_engine
!
      procedure :: set_printables                      => set_printables_reference_engine
!
      procedure :: generate_sad_density               => generate_sad_density_reference_engine
!
      procedure :: do_visualization                   => do_visualization_reference_engine
      procedure, nopass :: do_orbital_plotting        => do_orbital_plotting_reference_engine
!
      procedure :: do_ground_state                    => do_ground_state_reference_engine
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
!
      call engine%read_settings()
!
      call engine%set_printables()
!
      engine%timer = timings(trim(engine%name_))
      call engine%timer%turn_on()
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
      class(hf)        :: wf 
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
      class(reference_engine)  :: engine 
      class(hf)         :: wf 
!
      if (.not. engine%restart .and. (trim(engine%ao_density_guess) == 'sad')) then
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
!     Plot orbitals and/or density
!
      if (engine%plot_orbitals .or. engine%plot_density) call engine%do_visualization(wf)
!
!     Calculate the zeroth order properties
!
      if(engine%requested_zop) call engine%calculate_expectation_values(wf)
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
      call input%get_keyword_in_section('algorithm', 'solver scf', engine%algorithm)
      if (input%requested_keyword_in_section('restart', 'solver scf')) engine%restart = .true.
!
      call input%get_keyword_in_section('ao density guess', 'solver scf', engine%ao_density_guess)
!
      if (input%requested_keyword_in_section('plot hf orbitals', 'visualization')) engine%plot_orbitals = .true.
      if (input%requested_keyword_in_section('plot hf density', 'visualization')) engine%plot_density = .true.
!
      call engine%read_zop_settings()
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
      use atomic_class,           only: atomic
      use molecular_system_class, only: molecular_system
!
      use uhf_class,              only: uhf
      use scf_hf_class,           only: scf_hf
!
      use string_utilities,       only: index_of_unique_strings
!
      implicit none
!
      class(reference_engine)          :: engine
!
      class(hf)                        :: wf
!
      type(atomic)                     :: atom
      type(molecular_system)           :: sad_system
      type(uhf),         allocatable   :: sad_wf
      type(scf_hf), allocatable        :: sad_solver
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
      character(len=50), dimension(wf%system%n_atoms) :: atom_and_basis
      integer,           dimension(wf%system%n_atoms) :: unique_atom_index
!
      call engine%tasks%print_('sad')
!
!     SAD solver settings
!
      ao_density_guess   = 'core'
      max_iterations     = 100
!
      energy_threshold   = 1.0D-6
      call input%get_keyword_in_section('energy threshold', 'solver scf', energy_threshold)
      energy_threshold   = min(1.0D-6, energy_threshold)
!
      gradient_threshold = 1.0D-6
      call input%get_keyword_in_section('gradient threshold', 'solver scf', gradient_threshold)
      gradient_threshold = min(1.0D-6, gradient_threshold)
!
!     Find atomic index of unique atom/basis combinations
!
      do I = 1, wf%system%n_atoms
!
         atom_and_basis(I) = trim(wf%system%atoms(I)%symbol) // trim(wf%system%atoms(I)%basis)
!
      enddo
!
      call index_of_unique_strings(unique_atom_index, wf%system%n_atoms, atom_and_basis)
!
!     For every unique atom, generate SAD density to file
!
      do I = 1, wf%system%n_atoms
!
!        Check unique
!
         if ( all(unique_atom_index /= I)) cycle
!
         atom = wf%system%atoms(I)
!
         name_       = "sad_" // trim(atom%basis) // "_" // trim(atom%symbol)
         alpha_fname = trim(name_) // '_alpha'
         beta_fname  = trim(name_) // '_beta'
!
!        Prepare molecule of the chosen atom
!
         call output%mute()
!
         multiplicity = atom%get_multiplicity()
         sad_system   = molecular_system(atoms=[atom],              &
                                         name_=name_,               &
                                         charge=0,                  &
                                         multiplicity=multiplicity, &
                                         mm_calculation=.false.,    &
                                         pcm_calculation=.false.     )
!
!        Prepare SAD wavefunction
!
         sad_wf = uhf(sad_system, fractional_uniform_valence=.true.)
!
!        Prepare and run solver
!
         sad_solver = scf_hf(wf=sad_wf,                       &
                           restart=.false.,                   &
                           ao_density_guess=ao_density_guess, &
                           energy_threshold=energy_threshold, &
                           max_iterations=max_iterations,     &
                           gradient_threshold=gradient_threshold)
!
         call sad_solver%run(sad_wf)
!
!        Cleanup and generate ao_density_a and ao_density_b
!
         call sad_wf%cleanup()
         call sad_system%cleanup()
!
         call output%unmute()
!
         call output%printf('Generated atomic density for ' // adjustl(atom%symbol) // &
                            ' using UHF/(a0)', pl='verbose', fs='(t6,a)', chars=[atom%basis])
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
!     Libint is overwritten by SAD. Re-initialize.
!
      call wf%system%initialize_libint_atoms_and_bases()
      call wf%system%initialize_libint_integral_engines()
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
      character(len=200)      :: density_file_tag
!
      call engine%tasks%print_('plotting')
!
      if (trim(wf%name_) .eq. 'uhf') call output%error_msg('no plotting for UHF')
!
!     Initialize the plotter
!
      plotter = visualization(wf%system, wf%n_ao)
!
      if (engine%plot_orbitals) then
!
         call engine%do_orbital_plotting(plotter, wf)
!
      endif
!
      if (engine%plot_density) then
!
         density_file_tag = 'AO_density'
         call plotter%plot_density(wf%system,wf%ao_density, density_file_tag)
!
      endif
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
!
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
!     Read orbital plotting settings
!
      n_orbitals_to_plot = input%get_n_elements_for_keyword_in_section('plot hf orbitals', &
                                                                        'visualization')
!
      call mem%alloc(orbitals_to_plot, n_orbitals_to_plot)
      call input%get_array_for_keyword_in_section('plot hf orbitals', 'visualization', &
            n_orbitals_to_plot, orbitals_to_plot)
!
!     Extract the orbitals we will plot from the full orbital coefficient matrix
!
      call mem%alloc(orbital_coefficients, wf%n_ao, n_orbitals_to_plot)
!
      do i = 1, n_orbitals_to_plot
!
         call dcopy(wf%n_ao, wf%orbital_coefficients(1, orbitals_to_plot(i)), &
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
      call plotter%plot_orbitals(wf%system, orbital_coefficients, &
                                       n_orbitals_to_plot, orbital_file_tags)
!
      call mem%dealloc(orbital_coefficients, wf%n_ao, n_orbitals_to_plot)
      deallocate(orbital_file_tags)
!
   end subroutine do_orbital_plotting_reference_engine
!
!
   subroutine read_zop_settings_reference_engine(engine)
!!
!!    Read ZOP settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!    Modified by Tommaso Giovannini and Linda Goletto
!!
      implicit none
!
      class(reference_engine) :: engine
!
      engine%requested_zop = input%requested_section('hf zop')
!
      if (engine%requested_zop) then 
!
         if (input%requested_keyword_in_section('dipole','hf zop')) &
             engine%dipole = .true.
!
         if (input%requested_keyword_in_section('quadrupole','hf zop')) &
             engine%quadrupole = .true.
!
      endif
!
   end subroutine read_zop_settings_reference_engine
!
!
   subroutine calculate_expectation_values_reference_engine(engine, wf)
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
         call output%printf('The traceless quadrupole is calculated as:', pl='minimal',fs='(/t6,a)')
         call output%printf('Q_ij = 1/2[3*q_ij - tr(q)*delta_ij]', pl='minimal',fs='(/t9,a)')
         call output%printf('where q_ij is the non-traceless matrix', pl='minimal',fs='(/t6,a)')
!
         call engine%print_operator('traceless quadrupole moment', q_electronic, q_nuclear, q_total, &
                                    components, 6)
!
         deallocate(components)
!
      endif
!
   end subroutine calculate_expectation_values_reference_engine
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
!     Get the integrals mu_pqk for components k = 1, 2, 3
!
      call mem%alloc(mu_pqk, wf%n_ao, wf%n_ao, 3)
      call wf%get_ao_mu_wx(mu_pqk(:,:,1), mu_pqk(:,:,2), mu_pqk(:,:,3))
!
!     Get electronic expectation value contribution
!
      do k = 1, 3
!
         electronic(k) = wf%calculate_expectation_value(mu_pqk(:,:,k), wf%ao_density)
!
      enddo
!
      call mem%dealloc(mu_pqk, wf%n_ao, wf%n_ao, 3)
!
!     Get nuclear expectation value contribution, then sum the two
!
      if(wf%name_.eq.'mlhf') then
!
         call output%printf('In MLHF the dipole moment is computed only in the active space', pl='minimal',fs='(/t6,a)')
!
         call wf%system%get_nuclear_dipole_active(nuclear)
!
      else
!
         call wf%system%get_nuclear_dipole(nuclear)
!
      endif
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
      call mem%alloc(q_pqk, wf%n_ao, wf%n_ao, 6)
      call wf%get_ao_q_wx(q_pqk(:,:,1), q_pqk(:,:,2), q_pqk(:,:,3), &
                          q_pqk(:,:,4), q_pqk(:,:,5), q_pqk(:,:,6))
!
!     Get electronic expectation value contribution
!
      do k = 1, 6
!
         electronic(k) = wf%calculate_expectation_value(q_pqk(:,:,k), wf%ao_density)
!
      enddo
!
      call mem%dealloc(q_pqk, wf%n_ao, wf%n_ao, 6)
!
!     Get nuclear expectation value contribution, then sum the two
!
      if(wf%name_.eq.'mlhf') then
!
         call output%printf('In MLHF the quadrupole moment is computed only in the active space', pl='minimal',fs='(/t6,a)')
!
         call wf%system%get_nuclear_quadrupole_active(nuclear)
!
      else
!
         call wf%system%get_nuclear_quadrupole(nuclear)
!
      endif
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
      class(reference_engine), intent(in)    :: engine 
      class(hf), intent(inout)               :: wf 
!
      type(scf_hf),      allocatable :: scf
      type(scf_diis_hf), allocatable :: scf_diis
      type(mo_scf_diis), allocatable :: mo_scf_diis_
!
      call engine%tasks%print_('gs solver')
!
      if (trim(engine%algorithm) .eq. 'scf-diis' .and. trim(wf%name_) .eq. 'mlhf') then
!
         call output%error_msg('MLHF can not run with scf-diis, try mo-scf-diis.')
!
      elseif (trim(engine%algorithm) .eq. 'mo-scf-diis' .and. trim(wf%name_) .eq. 'uhf') then
!
         call output%error_msg('UHF can not run with mo-scf-diis, try scf-diis.')
!
      elseif (trim(engine%algorithm) == 'scf-diis') then
!
         scf_diis = scf_diis_hf(wf, engine%restart)
         call scf_diis%run(wf)
!
      elseif (trim(engine%algorithm) == 'mo-scf-diis') then
!
         mo_scf_diis_ = mo_scf_diis(wf, engine%restart)
         call mo_scf_diis_%run(wf)
!
      elseif (trim(engine%algorithm) == 'scf') then 
!
         scf = scf_hf(wf, engine%restart)
         call scf%run(wf)
!
      else
!
         call output%error_msg('did not recognize hf algorithm: '// engine%algorithm)
!
      endif
!
   end subroutine do_ground_state_reference_engine
!
!
end module reference_engine_class
