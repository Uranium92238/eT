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
      read(unit_input,'(a40)') line   
!
      do while (line(1:1) == '!' .or. trim(line) == '') ! Comment or blank line: read the next line
!
         read(unit_input,'(a40)') line 
!
      enddo
!
      if (trim(line) == 'Calculation:') then
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
!              Test for which type, set the logical in tasks, and cycle!
!
               if (calculation == 'ground_state') then
!
                  tasks%ground_state = .true.
                  cycle
!
               elseif (calculation == 'excited_state') then
!
                  tasks%excited_state = .true.
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
            elseif (trim(line) == 'Settings:') then 
!
               backspace(unit_input)
               exit ! escape the do loop
!
            else
! 
               write(unit_output,*) 'Input error: line ',trim(line),' not recognized.'
               stop
!
            endif
!
         enddo
!
      else
!
         write(unit_output,*) 'Expected calculation settings, not ',trim(line),'.'
!
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
      read(unit_input,'(a40)') line   
!
      do while (line(1:1) == '!' .or. trim(line) == '') ! Comment or blank line: read the next line
!
         read(unit_input,'(a40)') line 
!
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
!              Test for which type, set the logical in tasks, and cycle!
!
               if (setting == 'energy_threshold') then
!
                  read(unit_input,*) settings%energy_threshold
                  cycle
!
               elseif (setting == 'ampeqs_threshold') then 
!
                  read(unit_input,*) settings%ampeqs_threshold
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
               elseif (trim(line) == '#end of eT input') then
!
                  exit ! escape do loop 
!
               else
!
                  write(unit_output,*) 'Input error: setting ',trim(line(2:40)),' not recognized.'
                  stop
!
               endif
!
            elseif (trim(line) == '#end of eT input') then
!
               exit
!
            else
!
               write(unit_output,*) 'Input error: unregonized line ',trim(line),'.'
               stop
!
            endif
!
         enddo
!
      else
!
         write(unit_output,*) 'Input error: expected settings section, not the line ',trim(line),'.'
         stop ! Terminate program
!
      endif
!
   end subroutine settings_reader
!
!
end module input_reader