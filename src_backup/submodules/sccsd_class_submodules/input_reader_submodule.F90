submodule (sccsd_class) input_reader
!
!!
!!    Input reader submodule (SCCSD)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Dec 2017
!!
   use input_reader
!
contains 
!
!
   module subroutine scc_reader_sccsd(wf, unit_input)
!!
!!    SCC reader (SCCSD)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Dec 2017
!!
!!    Reads the SCCSD specific input section of input file.
!!
      implicit none
!
      integer(i15)      :: unit_input
!
      class(sccsd)      :: wf
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
         if (trim(line) == 'SCCSD') then
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
                     if (trim(line) == 'restart:') then
!
                        wf%scc_settings%restart = .true.
                        cycle
!
                     elseif (trim(line) == 'states:') then
!
                        read(unit_input,*) wf%state_A, wf%state_B
                        cycle
!
                     elseif (trim(line) == 'overlap_threshold:') then
!
                        read(unit_input,*) wf%scc_settings%overlap_threshold
                        cycle
!
                     elseif (trim(line) == 'triple_amplitude:') then
!
                        read(unit_input,*) wf%triples
                        cycle
!
                     elseif (trim(line) == 'virtual:') then
!
                        read(unit_input,*) wf%A, wf%B, wf%C 
                        cycle
!
                     elseif (trim(line) == 'occupied:') then
!
                        read(unit_input,*) wf%I, wf%J, wf%K 
                        cycle
!
                     elseif (trim(line) == '}') then
!
                        exit ! Exit general loop 2
!
                     elseif (trim(line) == 'end of eT input') then
!
                        write(unit_output,*) 'Error: could not find end of SCCSD section in eT.inp.'
                        stop  
!
                     endif
!
                  enddo
!
                  exit ! Exit general loop 1
!
               else
!
                  write(unit_output,*) 'Error: missing "{" for SCCSD section.'
                  stop
!
               endif 
!
         elseif (trim(line) == 'end of eT input') then
!
            write(unit_output,*) 'Error: no SCCSD input section found.'
            stop
!
         endif
! 
      enddo
!
   end subroutine scc_reader_sccsd
!
!
end submodule input_reader
