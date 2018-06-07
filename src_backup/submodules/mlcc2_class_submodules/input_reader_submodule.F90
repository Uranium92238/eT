submodule (mlcc2_class) input_reader
!
!!
!!    Input reader submodule (MLCC2)
!!    Written by Eirik F. Kjønstad and Sarai Dery Folkestad, May 2017
!!
!!    Contains the following family of procedures of the MLCC2 class:
!!
!!
   use input_reader
!
contains 
!
!
      module subroutine mlcc_reader_mlcc2(wf, unit_input)
!!
!!
!!
      implicit none
!
      integer(i15)      :: unit_input
!
      class(mlcc2)      :: wf
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
         if (trim(line) == 'MLCC2') then
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
                     if (trim(line) == 'CC2') then
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
                        write(unit_output,*)'Error: could not find end of MLCC2 section in eT.inp.'
                        stop  
!
                     endif
!
                  enddo
!
                  exit ! Exit general loop 1
!
               else
                  write(unit_output,*)'Error: missing "{" for MLCC2 section.'
                  stop
               endif 
!
         elseif (trim(line) == 'end of eT input') then
!
            write(unit_output,*)'Error: no MLCC2 input section found.'
            stop
!
         endif
! 
      enddo
!
   end subroutine mlcc_reader_mlcc2
!
!
   module subroutine read_orbital_info_mlcc2(wf, unit_input)
!!
!!
      implicit none
!
      integer(i15)      :: unit_input
!
      class(mlcc2)      :: wf
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
            exit
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
   end subroutine read_orbital_info_mlcc2
!
!
end submodule input_reader
