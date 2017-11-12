module input_reader
!
!!
!!    Input reader module                                 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017 
!!
!
   use types
   use workspace
   use input_output
   use calc_procedures_class
   use calc_settings_class
   use mlcc_calculation_settings_class
!
   implicit none
!
contains
!
!
   subroutine method_reader(unit_input, method)
!!
!!    Method Reader
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
      implicit none 
!
      integer(i15), intent(in) :: unit_input
!
      character(len=40) :: method
!
      character(len=40) :: line 
!
      do ! General do loop - ends when it reaches 'exit'
!
         read(unit_input,'(a40)') line 
!
         do while (line(1:1) == '!' .or. trim(line) == '') ! Comment or blank line: read the next line
            read(unit_input,'(a40)') line 
         enddo
!
         if (trim(line) == 'Method:') then
!
            read(unit_input,'(a40)') line
!
            do while (line(1:1) == '!' .or. trim(line) == '') ! Comment or blank line: read the next line
               read(unit_input,'(a40)') line 
            enddo
!
            if (line(1:1) == '.') then 
!
               method = trim(line(2:40))
               exit ! Escape from general do loop
!
            else
!
               write(unit_output,*) 'Input error: expected method, not the line ',trim(line),'.'
               stop ! Terminate program
!
            endif
!
         else
!
            write(unit_output,*) 'Input error: expected method section, not the line ',trim(line),'.'
            stop ! Terminate program
!
         endif
!
      enddo ! End general do loop
!
   end subroutine method_reader
!
!
   subroutine calculation_reader(unit_input, tasks)
!!
!!    Calculation Reader
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!    Reads the calculation from the input file and initializes the 
!!    tasks requested of the wavefunction.
!!
      implicit none 
!
      type(calc_procedures) :: tasks 
!
      integer(i15), intent(in) :: unit_input
!
      character(len=40) :: calculation 
!
      character(len=40) :: line
!
      logical :: calculation_info = .false.
!
      integer(i15) :: i
!
      rewind(unit_input)

      do
         read(unit_input,'(a40)') line   
!
         do while (line(1:1) == '!' .or. trim(line) == '') ! Comment or blank line: read the next line
            read(unit_input,'(a40)') line 
         enddo
!
         if (trim(line) == 'Calculation:') then
            calculation_info = .true.
!
            do ! Read calculations 
!
               read(unit_input,'(a40)') line
!
               do while (line(1:1) == '!' .or. trim(line) == '') ! Comment or blank line: read the next line
                  read(unit_input,'(a40)') line 
               enddo
!
               if (line(1:1) == '.') then 
!
                  calculation = trim(line(2:40))
!
!                 Test for which type, set the logical in tasks, and cycle!
!
                  if (calculation == 'ground_state') then
!
                     tasks%ground_state = .true.
                     cycle
!
                  elseif (calculation == 'excited_state') then
!
                     tasks%excited_state = .true.
                     read(unit_input,'(i3,i3)') tasks%n_singlet_states, tasks%n_triplet_states
!
                     read(unit_input,'(a40)') line

                     if (trim(line) == 'start_vector') then
!
                        tasks%user_specified_start_vector = .true.
                        call allocator_int(tasks%start_vectors, tasks%n_singlet_states, 1)
                        read(unit_input,*) tasks%start_vectors
!
                     else
!
                        backspace(unit_input)
!
                     endif
!
                     cycle 
!
                  elseif (calculation == 'core_excited_state') then
!
                     tasks%core_excited_state = .true.
                     read(unit_input,'(i3,i3)') tasks%n_singlet_states, tasks%n_triplet_states
                     read(unit_input,'(i3)') tasks%n_cores
                     call allocator_int(tasks%cores, tasks%n_cores, 1)
                     read(unit_input,*) (tasks%cores(i, 1), i = 1, tasks%n_cores)
                     cycle 
!
                  elseif (calculation == 'ionized_state') then
!
                     tasks%ionized_state = .true.
                     read(unit_input,'(i3,i3)') tasks%n_singlet_states, tasks%n_triplet_states
!
                     cycle
!
                  elseif (calculation == 'core_ionized_state') then
!
                     tasks%core_ionized_state = .true.
                     read(unit_input,'(i3,i3)') tasks%n_singlet_states, tasks%n_triplet_states
                     read(unit_input,'(i3)') tasks%n_cores
                     call allocator_int(tasks%cores, tasks%n_cores, 1)
                     read(unit_input,*) (tasks%cores(i, 1), i = 1, tasks%n_cores)
                     cycle 
!
                  elseif (calculation == 'properties') then
!
                     tasks%properties = .true. 
                     cycle 
!
                  else
!
                     write(unit_output,*) 'Input error: calculation ',trim(line(2:40)),' not recognized.'
                     stop
!
                  endif
!
               elseif (trim(line) == 'Settings:' .or. trim(line) == 'MLCC:') then
!
                  exit ! escape the do loop
!
               elseif ( trim(line) == '#end of eT input') then
!
                  backspace(unit_input)
                  exit
!
               endif
!
            enddo
!
         elseif ( trim(line) == '#end of eT input') then
            backspace(unit_input)
            exit
         endif
      enddo
!
      if (.not. calculation_info) then
         write(unit_output,*)'WARNING: Calculation information expected.'
         stop
      endif
!
   end subroutine calculation_reader
!
!
   subroutine settings_reader(unit_input, settings)
!!
!!    Settings Reader
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!    Reads the calculation settings from the input file and initializes 
!!    the settings requested of the wavefunction.
!!
      implicit none 
!
      type(calc_settings) :: settings 
!
      integer(i15), intent(in) :: unit_input
!
      character(len=40) :: setting 
!
      character(len=40) :: line 
!
      rewind(unit_input)
!
      do
         read(unit_input,'(a40)') line  
!
         do while (line(1:1) == '!' .or. trim(line) == '') ! Comment or blank line: read the next line
                  read(unit_input,'(a40)') line 
         enddo 
!
         if (trim(line) == 'Settings:') then
!
            do ! Read settings 
!
               read(unit_input,'(a40)') line
!
               do while (line(1:1) == '!' .or. trim(line) == '') ! Comment or blank line: read the next line
                  read(unit_input,'(a40)') line 
               enddo
!
               if (line(1:1) == '.') then 
!
                  setting = trim(line(2:40))
!
!                 Test for which type, set the logical in tasks, and cycle!
!
                  if (setting == 'energy_threshold') then
!
                     read(unit_input,*) settings%energy_threshold
                     cycle
!
                  elseif (setting == 'equation_threshold') then 
!
                     read(unit_input,*) settings%equation_threshold
                     cycle
!
                  elseif (setting == 'memory') then
!
                     read(unit_input,*) mem 
                     cycle
!
                  elseif (setting == 'restart') then
!
                     settings%restart = .true.
                     cycle
!
                  elseif (setting == 'ground_state_max_iterations') then 
!
                     read(unit_input,*) settings%ground_state_max_iterations
                     cycle
!
                  else
!
                     write(unit_output,*) 'Input error: setting ',trim(line(2:40)),' not recognized.'
                     stop
!
                  endif
               elseif (trim(line) == 'Calculation:' .or. trim(line) == 'MLCC:') then
!
                  exit
!
               elseif (trim(line) == '#end of eT input') then
!
                  backspace(unit_input)
                  exit ! escape do loop 
!
               endif
!
            enddo
!
         elseif (trim(line) == '#end of eT input') then
!
            backspace(unit_input)
            exit
!
         endif
!
      enddo
!
   end subroutine settings_reader
!
   subroutine mlcc_reader(unit_input,mlcc_settings)
!!
!!
   implicit none
!
      integer(i15) :: unit_input
      type(mlcc_calculation_settings) :: mlcc_settings
!
      character(len=40) :: setting 
!
      character(len=40) :: line 
!
      logical :: mlcc_info = .false.
      logical :: read_new_line = .true.
!
      rewind(unit_input)
!
      do
!
         read(unit_input,'(a40)') line
!
         do while (line(1:1) == '!' .or. trim(line) == '') ! Comment or blank line: read the next line
!
               read(unit_input,'(a40)') line
!
         enddo
!
         if (trim(line) == 'MLCC:') then
            mlcc_info = .true.
            do
!
            read(unit_input,'(a40)') line
!
!           Reading settings
!
            do while (line(1:1) == '!' .or. trim(line) == '') ! Comment or blank line: read the next line
!
               read(unit_input,'(a40)') line
!
            enddo

!
               if (line(1:1) == '.') then
!
                  setting = trim(line(2:40))
!
!                 Test for which type, set the logical in tasks, and cycle!
!
                  if (setting == 'cholesky') then
!
                     mlcc_settings%cholesky = .true.
                     cycle
!
                  elseif (setting == 'cnto') then 
!
                     mlcc_settings%cnto = .true.
!
!                    Get thresholds if in input
!
                     read(unit_input,'(a40)') line
                     if(trim(line)=='cnto_threshold_o') then
!
                        read(unit_input,*) mlcc_settings%delta_o
!
                        read(unit_input,'(a40)') line
!
                        if(trim(line)=='cnto_threshold_v') then
                          read(unit_input,*) mlcc_settings%delta_v
                        endif
! 
                     elseif(trim(line)=='cnto_threshold_v') then
!
                        read(unit_input,*) mlcc_settings%delta_v
!
                        read(unit_input,'(a40)') line
!
                        if(trim(line)=='cnto_threshold_o') then
                          read(unit_input,*) mlcc_settings%delta_o
                        endif
! 
                     elseif (trim(line) == '#end of eT input') then
!
                        backspace(unit_input)
                        exit
!
                     endif
!
!                    Done
!
                     cycle
!
                  elseif (setting == 'CCS') then
!
                     mlcc_settings%CCS = .true.
                     cycle
                  elseif (setting == 'CC2') then
!
                     mlcc_settings%CC2 = .true.
                     cycle
!
                  elseif (setting == 'CCSD') then
!
                     mlcc_settings%CCSD = .true.
                     cycle
!
                  elseif (setting == 'CC3') then
!
                     mlcc_settings%CC3 = .true.
                     cycle
!
                  elseif (trim(line) == 'Calculation:' .or. trim(line) == 'Settings:') then
!
                     exit
!
                  elseif (trim(line) == '#end of eT input') then
!
                     backspace(unit_input)
                     exit
!
                  endif
!
               elseif (trim(line) == '#end of eT input') then
!
                  backspace(unit_input)
                  exit
!
               endif
!
            enddo
!           
         elseif ( trim(line) == '#end of eT input') then
            backspace(unit_input)
            exit
         endif
      enddo
!
      if (.not. mlcc_info) then 
         write(unit_output,*) 'Input error: expected mlcc settings section.'
         stop ! Terminate program
      endif
!

!
   end subroutine
!
end module input_reader
