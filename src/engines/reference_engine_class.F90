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
   use global_in,         only: input
   use global_out,        only: output
   use timings_class,     only: timings
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
      logical :: restart
!
   contains 
!
      procedure :: ignite                       => ignite_reference_engine
!
      procedure :: run                          => run_reference_engine
      procedure :: read_settings                => read_settings_reference_engine
!
      procedure :: set_printables               => set_printables_reference_engine
!
      procedure, nopass :: generate_sad_density => generate_sad_density_reference_engine
!
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
      type(scf_hf),      allocatable :: scf
      type(scf_diis_hf), allocatable :: scf_diis
      type(mo_scf_diis), allocatable :: mo_scf_diis_
!
!     Generate SAD if requested
!
      if (.not. engine%restart .and. (trim(engine%ao_density_guess) == 'sad')) then
!
         call engine%generate_sad_density(wf)
!
      endif
!
!     Choose solver
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
         call scf_diis%cleanup(wf)
!
      elseif (trim(engine%algorithm) == 'mo-scf-diis') then
!
         mo_scf_diis_ = mo_scf_diis(wf, engine%restart)
         call mo_scf_diis_%run(wf)
         call mo_scf_diis_%cleanup(wf)
!
      elseif (trim(engine%algorithm) == 'scf') then 
!
         scf = scf_hf(wf, engine%restart)
         call scf%run(wf)
         call scf%cleanup(wf)
!
      else
!
         call output%error_msg('did not recognize hf algorithm: '// engine%algorithm)
!
      endif
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
      implicit none
!
      class(reference_engine) :: engine
!
      engine%name_       = 'Reference state engine'
      engine%author      = 'E. F. Kjønstad, S. D. Folkestad, 2018'
!
      engine%description = 'Calculates the reference wavefunction | R >.'
      engine%tag         = 'ground state'
!
      engine%tasks       = [character(len=150) ::                                                         &
                           'Generate initial density (' // trim(engine%ao_density_guess) // ')',          &
                           'Calculation of reference state (' // trim(engine%algorithm) // ' algorithm)', &
                           'Calculation of the ground state energy']
!
   end subroutine set_printables_reference_engine
!
!
   subroutine generate_sad_density_reference_engine(wf)
!!    
!!    Generate SAD density
!!    Written by Tor S. Haugland, 2019
!!
!!    Generates SAD density for every unique (atom, basis) pair and sets the wf density.
!!    First, a molecular system with only one atom for every (atom, basis) is created.
!!    Secondly, a UHF wavefunction is created for that system.
!!    Thirdly, the wavefunction is run through a HF solver.
!!    Finally, the density files are moved to where the wavefunction expects them to be.
!!
!!    The renaming of files and the cleanup and re-prepare of the full system
!!    is present to avoid overwriting files.
!!
!!    Modified by Tor S. Haugland, 2019
!!
!!    Removed deletion of files generated in SCF for SAD to fix a re-occuring bug. Restart files
!!    are now checked from a parameter list. Added "Found SAD". Only loop over unique atoms.
!!
!
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
      class(hf)                        :: wf
!
      type(atomic)                     :: atom
      type(molecular_system)           :: sad_system
      type(uhf),         allocatable   :: sad_wf
      type(scf_hf), allocatable        :: sad_solver
!
      type(sequential_file) :: alpha_density_file
      type(sequential_file) :: beta_density_file
      type(sequential_file) :: restart_file
!
      character(len=200) :: alpha_fname
      character(len=200) :: beta_fname
!
      character(len=200) :: ao_density_guess
      character(len=200) :: name_
      integer            :: multiplicity
!
      real(dp) :: energy_threshold
      real(dp) :: gradient_threshold
      integer  :: max_iterations
!
      integer :: I
!
      character(len=50), dimension(wf%system%n_atoms) :: atom_and_basis
      integer,           dimension(wf%system%n_atoms) :: unique_atom_index
!
      character(len=20), dimension(2), parameter :: restart_files = ["scf_restart_file   ", &
                                                                     "orbital_information"  ]
!
      call output%printf('- Generating SAD guess', pl='minimal', fs='(/t3,a)')
!
!     SAD DIIS solver settings
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
!     Rename Hartree-Fock restart files so that they are not overwritten:
!
      do I = 1, size(restart_files)
!
         restart_file = sequential_file(trim(restart_files(I)))
!
         if (restart_file%exists()) then
!
            call output%printf('Found restart, copy to temp. position: temp_(a0)', &
                              pl='debug', chars=[trim(restart_files(I))], fs='(t6,a)' )
!
            call restart_file%copy("temp_" // trim(restart_files(I)))
!
            call restart_file%delete_()
!
         endif
!
      enddo
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
!        if SAD already exist, skip to next atom
!
         alpha_density_file = sequential_file(alpha_fname)
         beta_density_file  = sequential_file(beta_fname)
!
         if (alpha_density_file%exists() .and. beta_density_file%exists()) then
!
            call output%printf('Found SAD for '// adjustl(atom%symbol) &
               // ' in ' // trim(atom%basis), pl='verbose', fs='(t6,a)')
!
            cycle
!
         elseif (alpha_density_file%exists() .or. beta_density_file%exists()) then
!
            call output%warning_msg("Could only find 1 of 2 density files for " // trim(name_) // &
                                  ". Deleting it and creating both SAD densities.")
!
            if (alpha_density_file%exists()) call alpha_density_file%delete_()
            if (beta_density_file%exists())  call beta_density_file%delete_()
!
         endif
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
                                         mm_calculation=.false.     )
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
         call sad_solver%cleanup(sad_wf)
         call sad_wf%cleanup()
         call sad_system%cleanup()
!
         call output%unmute()
!
         call output%printf('Generated SAD for '// adjustl(atom%symbol) //' in '&
                  // trim(atom%basis), pl='verbose', fs='(t6,a)')
!
!        Move densities to where "set_ao_density_sad" can use them
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
!     Rename Hartree-Fock restart files back to their original names
!
      do I = 1, size(restart_files)
!
         restart_file = sequential_file('temp_' // trim(restart_files(I)))
!
         if (restart_file%exists()) then
!
            call output%printf('Found temp. restart, copy it back:     (a0)', &
                              pl='debug', chars=[restart_files(I)], fs='(t6,a)' )
!
            call restart_file%copy(trim(restart_files(I)))
!
            call restart_file%delete_()
!
         endif
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
end module reference_engine_class
