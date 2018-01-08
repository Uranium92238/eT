submodule (mlccsd_class) input_reader
!
!!
!!    Input reader submodule (mlccsd)
!!    Written by Eirik F. Kjønstad and Sarai Dery Folkestad, May 2017
!!
!!    Contains the following family of procedures of the mlccsd class:
!!
!!
   use input_reader
!
contains 
!
!
      module subroutine mlcc_reader_mlccsd(wf, unit_input)
!!
!!
!!
      implicit none
!
      integer(i15) :: unit_input
!
      class(mlccsd)  :: wf
!
      character(len=40) :: line
!
!     Start at the begining of eT.inp      
!   
!
      rewind(unit_input)
!
      do ! General loop 1
!
         read(unit_input,'(a40)') line 
!
!        Remove blanks preceding text
!
         line = remove_preceding_blanks(line)
!
         if (trim(line) == 'MLCCSD') then
! 
            read(unit_input,'(a40)') line 
            line = remove_preceding_blanks(line) 
!
               if (trim(line) == '{') then
!
                  do ! General loop 2
!
                     read(unit_input,'(a40)') line 
                     line = remove_preceding_blanks(line)
!
                     if (trim(line) == 'CCSD') then
!
                        wf%mlcc_settings%CCSD = .true.
                        cycle
!
                      elseif (trim(line) == 'CC2') then
!
                        wf%mlcc_settings%CC2 = .true.
                        cycle
!
                     elseif (trim(line) == 'CCS') then
!
                        wf%mlcc_settings%CCS = .true.
                        cycle
!
                     elseif (trim(line) == '}') then
!
                        exit ! Exit general loop 2
!
                     elseif (trim(line) == 'end of eT input') then
!
                        write(unit_output,*)'Error: could not find end of MLCCSD section in eT.inp.'
                        stop  
!
                     endif
!
                  enddo
!
                  exit ! Exit general loop 1
!
               else
                  write(unit_output,*)'Error: missing "{" for MLCCSD section.'
                  stop
               endif 
!
         elseif (trim(line) == 'end of eT input') then
!
            write(unit_output,*)'Error: no MLCCSD input section found.'
            stop
!
         endif
! 
      enddo
!
   end subroutine mlcc_reader_mlccsd
!
!
   module subroutine read_orbital_info_mlccsd(wf, unit_input)
!!
!!
      implicit none
!
      integer(i15)      :: unit_input
!
      class(mlccsd)      :: wf
!
      character(len=40) :: line
!
!     Start at the begining of eT.inp      
!
      rewind(unit_input)
!
      do ! General loop 1
!
         read(unit_input,'(a40)') line 
!
!        Remove blanks preceding text
!
         line = remove_preceding_blanks(line)
!
         if (trim(line) == 'CC2 orbitals') then
!
            call wf%CC2_orbitals%orbital_reader(unit_input)
            cycle
!
         elseif (trim(line) == 'CCSD orbitals') then
!
            call wf%CCSD_orbitals%orbital_reader(unit_input)
            cycle
!
         elseif (trim(line) == 'end of eT input') then
!
            backspace(unit_input)
            exit
!
         endif
      enddo
!
   end subroutine read_orbital_info_mlccsd
!
!
!
end submodule input_reader