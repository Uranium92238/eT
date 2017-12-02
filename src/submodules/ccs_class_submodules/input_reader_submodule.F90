submodule (ccs_class) input_reader
!
!!
!!    Input reader submodule (CCS)
!!    Written by Eirik F. Kjønstad and Sarai Dery Folkestad, May 2017
!!
!!    Contains the following family of procedures of the CCS class:
!!
!!
   use input_reader
!
contains 
!
!
   module subroutine general_specs_reader_ccs(wf, unit_input)
!!
!!    General Specifications Reader
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Nov. 2017
!!
!!    Reads general settings 
!!       - memory
!!       - disk space
!!    from the eT.inp file.
!!
!!    Both variables are given in gb in the input file and must be converted to words.
!!
      implicit none
!
      class(ccs)   :: wf
      integer(i15) :: unit_input
!
      integer :: memory = 0, disk_space = 0
!
      character(len=40) :: line
!
!     Start at the beginning of eT.inp
!
      rewind(unit_input)
!
      do ! General do loop - ends when it reaches 'exit'. Either when method is specified or if no method is given.
!
         read(unit_input,'(a40)') line 
!
!        Remove blanks preceding text
!
         line = remove_preceding_blanks(line)
!
!        Find general input section
!
         if(trim(line) == 'CC' .or. trim(line) == 'MLCC' .or. trim(line) == 'SCC') then
!
            read(unit_input,'(a40)') line 
            line = remove_preceding_blanks(line)
!
            if (trim(line) == '{') then
!
               do ! General do loop - ends when it reaches 'exit'. 
!
                  read(unit_input,'(a40)') line
                  line = remove_preceding_blanks(line)
!
                  if (trim(line) == 'disk_space:') then
!
                     read(unit_input,*) disk_space
!
!                    Converting from gb to words
!
                     wf%settings%disk_space = disk_space*134500000
!
                     cycle
!
                  elseif (trim(line) == 'memory:') then
!
                     read(unit_input,*) memory
!
!                    Converting from gb to words
!
                     wf%settings%memory = memory*134500000
!
                     cycle
!
                  elseif (trim(line) == 'print_level:') then
!
                     read(unit_input,*) wf%settings%print_level
!
                     cycle             
!
                  elseif (trim(line) == 'end of eT input') then ! Use defaults for memory and disk_space
!
                     exit
!
                  endif 
!
               enddo ! End general do loop
!
               exit
!
            else       
!
               write(unit_output,*)'Error: no general input section found.'
               stop
            endif
!
         elseif (trim(line) == 'end of eT input') then
!
            write(unit_output,*)'Error: no general input section found.'
            stop
!
         endif
!
      enddo
      write(unit_output,'(t3,a,i4,a)')  'Memory available for calculation:     ',memory,' gb'
      write(unit_output,'(t3,a,i4,a/)') 'Disk space available for calculation: ',disk_space,' gb'
      flush(unit_output)
!
   end subroutine general_specs_reader_ccs
!
!
 
   module subroutine calculation_reader_ccs(wf, unit_input)
                                 
!!
!!    Calculation reader,
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Nov. 2017
!!
!!    Reads eT.inp to determine which calculation tasks
!!
!!       - ground state
!!       - excited state
!!       - ionized state
!!       - property calculation
!!
!!     to do, and then calls specific readers for each of these tasks to set task specific settings.
!!
      implicit none
!
      integer(i15) :: unit_input
!
      class(ccs) :: wf
!
      character(len=40) :: line
!
!     Prints
!
      write(unit_output,'(t3,a/)')'Calculations requested:'
!
!     Always do ground state calculation
!
      wf%tasks%ground_state = .true.
!
!     Start at the begining of eT.inp
!
      rewind(unit_input)
!
      do ! General do loop - ends when it reaches 'exit'. 
!
         read(unit_input,'(a40)') line 
!
!        Remove blanks preceding text
!
         line = remove_preceding_blanks(line)
!
!        Find ground state/excited state/ionized state/property input section
!
         if (trim(line) == 'ground state') then ! Ground State
!
!           Prints
!
            write(unit_output,'(t6,a,a)')  trim(wf%name),' ground state'
!
            read(unit_input,'(a40)') line
            line = remove_preceding_blanks(line)
!
            if (trim(line) == '{') then ! Specifications given for ground state
!
               call wf%read_ground_state_specs(unit_input)
!
               cycle
!
            else ! Use defaults for ground state calculation
!
               cycle
!
            endif
!
         elseif (trim(line) == 'excited state') then ! Excited state
!
!           Prints
!
            write(unit_output,'(t6,a,a)')  trim(wf%name),  ' excited state'
!
             wf%tasks%excited_state = .true. ! This will be set to false again in the case of core_excited_state = .true.
!
            read(unit_input,'(a40)') line
            line = remove_preceding_blanks(line)
!
            if (trim(line) == '{') then ! Specifications and calculation details given for excited state
!
               call wf%read_excited_state_specs(unit_input)
!
               cycle
!
            else ! No excited state calculation without any information in eT.inp
!
               write(unit_output,*)'Error: excited state calculation section of eT.inp is missing.'
               stop
!
            endif
!
         elseif (trim(line) == 'ionized state') then ! Ionized state
!
!           Prints
!
            write(unit_output,'(t6,a,a)')  trim(wf%name),  ' ionized state'
!
             wf%tasks%ionized_state = .true. ! This will be set to false again in the case of core_ionized_state = .true.
!
            read(unit_input,'(a40)') line
            line = remove_preceding_blanks(line)
!
            if (trim(line) == '{') then ! Specifications and calculation details given for ionized state
!
               call wf%read_excited_state_specs(unit_input)
!
               cycle
!
            else ! No ionized state calculation without any information in eT.inp
!
               write(unit_output,*)'Error: ionized state calculation section of eT.inp missing.'
               stop
!
            endif
!
         elseif (trim(line) == 'property') then ! Property calculation
!
!           Prints
!
            write(unit_output,'(t6,a)')  '- Property'
!
            read(unit_input,'(a40)') line
            line = remove_preceding_blanks(line)
!
            if (trim(line) == '{') then ! Specifications and calculation details given for property calculation
!
               call wf%read_property_specs(unit_input)
!
               cycle
!
            else ! No property calculation without any information in eT.inp
!
               write(unit_output,*)'Error: property calculation section of eT.inp missing.'
               stop
!
            endif
         elseif (trim(line) == 'end of eT input') then
            backspace(unit_input)
            exit
         endif
!
      enddo
!
   end subroutine calculation_reader_ccs
!
!
   module subroutine read_ground_state_specs_ccs(wf, unit_input)
!!
!!    Read ground state specifications.
!!    Written by Eirik F. Kjønstad and Sarai Dery Folkestad, Nov. 2017
!!
!!    Reads settings
!!
!!       - restart
!!       - energy and residual thresholds
!!       - max number of iterations
!!
!!    for the ground state calculation
!!
      implicit none
!
      integer :: unit_input
!
      class(ccs) :: wf
!
      character(len=40) :: line
!
!     File is open, and positioned at the top of the ground state section
!     thus we will NOT rewind.
!
      do ! General do loop - ends when it reaches 'exit'. 
!
         read(unit_input,'(a40)') line 
!
!        Remove blanks preceding text
!
         line = remove_preceding_blanks(line)
!
         if (trim(line) == 'restart') then ! Restart
!
            wf%ground_state_specifications%restart = .true.
            cycle
!
         elseif (trim(line) == 'energy_threshold:') then ! Energy threshold
!
            read(unit_input, *) wf%ground_state_specifications%energy_threshold
            cycle
!
         elseif (trim(line) == 'residual_threshold:') then ! Residual threshold
!
            read(unit_input, *) wf%ground_state_specifications%residual_threshold
            cycle
!
         elseif (trim(line) == 'max_iterations:') then ! Max number of iterations
!
            read(unit_input, *) wf%ground_state_specifications%max_iterations
            cycle
!
         elseif (trim(line) == '}') then ! Done 
!
            return
!
         elseif (trim(line) == 'enf of eT input') then
!
            write(unit_output,*)'Error: could not find the end of the ground state calculation section of eT.inp .'
            stop
!
         endif
!
      enddo
!
   end subroutine read_ground_state_specs_ccs
!
   module subroutine read_excited_state_specs_ccs(wf, unit_input)
!!
!!    Read excited state specifications,
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Nov. 2017
!!
!!    Reads settings
!!
!!       - restart
!!       - energy and residual thresholds
!!       - max number of iterations
!!
!!    for the excited state calculation and specifications
!!      
!!       - n_singlet_states
!!       - n_triplet_states
!!       - start_vectors
!!    
!!    and for core excited states 
!!
!!       - n_equivalent_cores
!!       - cores
!!
!!     is also read. These MUST be present (in this order) in eT.inp if core excited/ionized state is requested.
!!
      implicit none
!
      integer(i15)      :: unit_input
!
      class(ccs)        :: wf
!
      character(len=40) :: line
!
!     File is open, and positioned at the top of the ground state section
!     thus we will NOT rewind.
!
      do ! General do loop - ends when it reaches 'exit'. 
!
         read(unit_input,'(a40)') line 
!
!        Remove blanks preceding text
!
         line = remove_preceding_blanks(line)
!
         if (trim(line) == 'restart') then ! Restart
!
            wf%excited_state_specifications%restart = .true.
            cycle
!
         elseif (trim(line) == 'energy_threshold:') then ! Energy threshold
!
            read(unit_input, *) wf%excited_state_specifications%energy_threshold
            cycle
!
         elseif (trim(line) == 'residual_threshold:') then ! Residual threshold
!
            read(unit_input, *) wf%excited_state_specifications%residual_threshold
            cycle
!
         elseif (trim(line) == 'max_iterations:') then ! Max number of iterations
!
            read(unit_input, *) wf%excited_state_specifications%max_iterations
            cycle
!
         elseif (trim(line) == 'n_singlet_states:') then ! Number of singlets
!
            read(unit_input, *) wf%excited_state_specifications%n_singlet_states
            cycle
!
         elseif (trim(line) == 'n_triplet_states:') then ! Number of triplets
!
            read(unit_input, *) wf%excited_state_specifications%n_triplet_states
            cycle
!
         elseif (trim(line) == 'start_vectors:') then ! Start vector provided in eT.inp
!
            wf%excited_state_specifications%user_specified_start_vector = .true.
!
            if (wf%excited_state_specifications%n_singlet_states == 0) then
!
               write(unit_output,*)'Error: Number of singlet excited states must be specified before start vectors are given'
               stop
!
            endif
!
            call allocator_int(wf%excited_state_specifications%start_vectors, wf%excited_state_specifications%n_singlet_states, 1)
!
            read(unit_input, *) wf%excited_state_specifications%start_vectors
            cycle
!
         elseif (trim(line) == 'core excited state' .or. trim(line) == 'core ionized state') then ! Requested core excitations 
!
!           Set calculation tasks
!
            if (wf%tasks%excited_state) then
!
               wf%tasks%excited_state = .false.
               wf%tasks%core_excited_state = .true.
!
            elseif (wf%tasks%ionized_state) then
!
               wf%tasks%ionized_state = .false.
               wf%tasks%core_ionized_state = .true.
!
            endif
!
!           This input section har strict structure:
!           E.g. a calculation of core excitations from the second atom given MOLECULE.INP
!           when this atom is the only atom of the this type with exactly the same environment (no equivalent cores).
!
!           n_equivalent_cores:
!           1
!           cores:
!           2
!
            read(unit_input,'(a40)') line 
            line = remove_preceding_blanks(line)
!
            if (trim(line) == '{') then
!
               do
                  read(unit_input,'(a40)') line 
                  line = remove_preceding_blanks(line)
!
                  if (trim(line) == 'n_equivalent_cores:') then ! Number of equivalent cores 
!
                     read(unit_input,*) wf%core_excited_state_specifications%n_equivalent_cores
                     cycle
!
                  elseif (trim(line) == 'cores:') then ! Cores (given by the order of which it appears in MOLECULE.INP)
!
                     if (wf%core_excited_state_specifications%n_equivalent_cores == 0) then
                        write(unit_output,*)'Number of equivalent cores should be selected before atoms are selected.'
                        stop
                     endif
!
                     call allocator_int(wf%core_excited_state_specifications%cores, &
                                       wf%core_excited_state_specifications%n_equivalent_cores, 1)
!
                     read(unit_input,*) wf%core_excited_state_specifications%cores
                     cycle
!
                  elseif (trim(line) == '}')  then
!
                     exit
!
                  elseif (trim(line) == 'enf of eT input') then
                     backspace(unit_input)
                     exit
!
                  endif
               enddo
!
!              Sanity check for core excited states
!
               if (wf%core_excited_state_specifications%n_equivalent_cores ==0 &
                  .or. .not. allocated(wf%core_excited_state_specifications%cores)) then
                  write(unit_output,*)'Error: atoms for core excitation was not found in eT.inp'
               endif
!
            else
!
               write(unit_output,*)'Error: core excited state input section is incomplete.'
               stop
!
            endif
            cycle
!
!
         elseif (trim(line) == '}') then ! Done
            exit
         elseif (trim(line) == 'enf of eT input') then ! Done
!
            backspace(unit_input)
            exit
!
         endif
!
      enddo
!
!     Sanity check:
!
      if (wf%excited_state_specifications%n_singlet_states == 0 .and. &
          wf%excited_state_specifications%n_triplet_states == 0) then
!
         write(unit_output,*)'Error: please specify either n_singlet_states or n_triplet_states in eT.inp'
         stop
!
      endif
!
   end subroutine read_excited_state_specs_ccs
!
!
module subroutine read_property_specs_ccs(wf, unit_input)
!!
!!    Read property specifications,
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Nov. 2017
!!
!!    Reads which properties to calculate. 
!!
      implicit none
!
      integer(i15)      :: unit_input
!
      class(ccs)        :: wf
!
      character(len=40) :: line
!
!     File is open, and positioned at the top of the ground state section
!     thus we will NOT rewind.
!
      do ! General do loop - ends when it reaches 'exit'. 
!
         read(unit_input,'(a40)') line 
!
!        Remove blanks preceding text
!
         line = remove_preceding_blanks(line)
!
         if (trim(line) == 'multipliers') then
!
            wf%tasks%multipliers = .true.
!
         elseif (trim(line) == 'end of eT input') then
!
            backspace(unit_input)
            exit
!
         endif
!
      enddo
!
!     Sanity check:
!
      if (wf%excited_state_specifications%n_singlet_states == 0 .and. &
          wf%excited_state_specifications%n_triplet_states == 0) then
!
         write(unit_output,*)'Property calculation requires specification of excited state calculation in eT.inp'
         stop
!
      endif
!
   end subroutine read_property_specs_ccs
!
!
end submodule input_reader