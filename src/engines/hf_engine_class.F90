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
module hf_engine_class
!!
!!    Hartree-Fock engine class module 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
   use parameters
!
   use global_out,             only: output
   use global_in,              only: input
   use sequential_file_class,  only: sequential_file
   use memory_manager_class,   only: mem
!
   use atomic_class,           only: atomic
   use molecular_system_class, only: molecular_system
   use hf_class,               only: hf
   use uhf_class,              only: uhf
   use scf_diis_hf_class,      only: scf_diis_hf
   use scf_hf_class,           only: scf_hf
!
!
   type hf_engine 
!
      character(len=200) :: ao_density_guess
      character(len=200) :: algorithm 
      logical :: restart
!
   contains 
!
      procedure :: ignite                       => ignite_hf_engine
!
      procedure, private :: run                 => run_hf_engine
!
      procedure, private :: read_settings       => read_settings_hf_engine
!
      procedure, nopass :: generate_sad_density => generate_sad_density_hf_engine
!
   end type hf_engine 
!
!
   interface hf_engine
!
      procedure :: new_hf_engine 
!
   end interface hf_engine
!
!
contains
!
!
   function new_hf_engine() result(engine)
!!
!!    New HF engine 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
      implicit none 
!
      type(hf_engine) :: engine
!
      engine%ao_density_guess = 'sad'
      engine%algorithm        = 'scf-diis'
      engine%restart          = .false.
!
      call engine%read_settings()
!
   end function new_hf_engine
!
!
   subroutine run_hf_engine(engine, wf)
!!
!!    Run 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none 
!
      class(hf_engine)  :: engine 
      class(hf)         :: wf 
!
      type(scf_hf), allocatable :: scf
      type(scf_diis_hf), allocatable :: scf_diis
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
      if (trim(engine%algorithm) == 'scf-diis') then
!
         scf_diis = scf_diis_hf(wf, engine%restart)
         call scf_diis%run(wf)
         call scf_diis%cleanup(wf)
!
      elseif (trim(engine%algorithm) == 'scf') then 
!
         if (engine%restart) call output%error_msg('SCF does not support restart.')
!
         scf = scf_hf(wf)
         call scf%run(wf)
         call scf%cleanup(wf)
!
      else
!
         call output%error_msg('did not recognize hf algorithm: '// engine%algorithm)
!
      endif
!
   end subroutine run_hf_engine
!
!
   subroutine read_settings_hf_engine(engine)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
      implicit none
!
      class(hf_engine) :: engine 
!
      call input%get_keyword_in_section('algorithm', 'solver hf', engine%algorithm)
      if (input%requested_keyword_in_section('restart', 'solver hf')) engine%restart = .true.
!
      call input%get_keyword_in_section('ao density guess', 'solver hf', engine%ao_density_guess)
!
   end subroutine read_settings_hf_engine
!
!
   subroutine ignite_hf_engine(engine, wf)
!!
!!    Ignite
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Apr 2019
!!
      implicit none 
!
      class(hf_engine) :: engine 
      class(hf)        :: wf 
!
      call engine%run(wf)
!
   end subroutine ignite_hf_engine
!
!
   subroutine generate_sad_density_hf_engine(wf)
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
      type(sequential_file) :: other_file
!
      character(len=200) :: alpha_fname
      character(len=200) :: beta_fname
!
      character(len=200) :: ao_density_guess
      character(len=200) :: name
      integer            :: multiplicity
!
      real(dp) :: energy_threshold
      real(dp) :: gradient_threshold
      integer  :: max_iterations
!
      integer :: i
!
!     SAD DIIS solver settings
!
      ao_density_guess   = 'core'
      max_iterations     = 100
!
      energy_threshold   = 1.0D-6
      call input%get_keyword_in_section('energy threshold', 'solver hf', energy_threshold)
      energy_threshold   = min(1.0D-6, energy_threshold)
!
      gradient_threshold = 1.0D-6
      call input%get_keyword_in_section('gradient threshold', 'solver hf', gradient_threshold)
      gradient_threshold = min(1.0D-6, gradient_threshold)
!
!     Mute output
!
      call output%printf('- SAD generation', pl='normal', fs='(/t3,a)')
!
      call output%mute()
!
!     Rename Hartree-Fock restart files so that they are not overwritten
!
      restart_file = sequential_file("hf_restart_file")
!
      if (restart_file%exists()) then
         call restart_file%copy("temp_restart_file")
         call restart_file%delete_()
      endif
!
!     For every unique atom, generate SAD density to file
!
      do i = 1, wf%system%n_atoms
!
         atom = wf%system%atoms(i)
!
         name        = "sad_" // trim(atom%basis) // "_" // trim(atom%symbol)
         alpha_fname = trim(name) // '_alpha'
         beta_fname  = trim(name) // '_beta'
!
!        if SAD already exist, skip to next atom
!
         alpha_density_file = sequential_file(alpha_fname)
         beta_density_file  = sequential_file(beta_fname)
!
         if (alpha_density_file%exists() .AND. beta_density_file%exists()) cycle
!
!        Prepare molecule of the chosen atom
!
         multiplicity = atom%get_multiplicity()
         sad_system   = molecular_system(atoms=[atom],              &
                                         name=name,                 &
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
         sad_solver = scf_hf(wf=sad_wf, &
                           ao_density_guess=ao_density_guess, &
                           energy_threshold=energy_threshold, &
                           max_iterations=max_iterations,     &
                           gradient_threshold=gradient_threshold)
!
         call sad_solver%run(sad_wf)
!
!        Print which density was generated
!
         call output%unmute()
!
         call output%printf('Generated SAD in (a0) for (a0)', chars=[character(len=100) :: atom%basis, atom%symbol], &
                                             pl='normal', ffs='(/t6,a)', fs='(t6,a)')
!
         call output%mute()
!
!        Cleanup and generate ao_density_a and ao_density_b
!
         call sad_solver%cleanup(sad_wf)
         call sad_wf%cleanup()
         call sad_system%cleanup()
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
!     Remove junk from SAD-generation
!
      other_file = sequential_file("ao_density")
      call other_file%delete_()
!
      other_file = sequential_file("orbital_coefficients")
      call other_file%delete_()
!
      other_file = sequential_file("orbital_energies")
      call other_file%delete_()
!
      restart_file = sequential_file('hf_restart_file')
      call restart_file%delete_()
!
!     Rename Hartree-Fock restart files back to their original names
!
      restart_file = sequential_file('temp_restart_file')
!
      if (restart_file%exists()) then
         call restart_file%copy('hf_restart_file')
         call restart_file%delete_()
      endif
!
!     Libint is overwritten by SAD. Re-initialize.
!
      call wf%system%initialize_libint_atoms_and_bases()
      call wf%system%initialize_libint_integral_engines()
!
!
      call output%unmute()
!
   end subroutine generate_sad_density_hf_engine
!
!
end module hf_engine_class
