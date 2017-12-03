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
!
   implicit none
!
contains
!
!
   subroutine method_reader(unit_input, method)
!!
!!    Method Reader
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Nov 2017
!!
      implicit none 
!
      integer(i15), intent(in) :: unit_input
!
      character(len=40) :: method
!
      character(len=40) :: line
!
!     Start at the begining of eT.inp
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
         if (trim(line) == 'CC') then
!
            read(unit_input,'(a40)') line
            line = remove_preceding_blanks(line)
!
            if (trim(line) == '{') then
!
               do ! General do loop - ends when it reaches 'exit'. Either when method is specified or if no method is given.
!
                  read(unit_input,'(a40)') line
                  line = remove_preceding_blanks(line)
!
                  if (trim(line) == 'method:') then ! Set method 
!
                     read(unit_input,'(a40)') line
                     line = remove_preceding_blanks(line)
!
                     method = trim(line)
!
                     exit
!
                  elseif (trim(line) == '}') then
                     write(unit_output,*)'Error: method was not specified in eT.inp.'
                     stop
                  endif 
!
               enddo ! End general do loop
!
            else 
!
               write(unit_output,*)'Error: method was not specified in eT.inp.'
               stop
!
            endif
!
!           
            exit
!
         elseif (trim(line) == 'MLCC') then
!
!
!           Determine what type of MLCC method
!
            read(unit_input,'(a40)') line
            line = remove_preceding_blanks(line)
!
            if (trim(line) == '{') then
!
               do ! General do loop - ends when it reaches 'exit'. Either when method is specified or if no method is given.
                  read(unit_input,'(a40)') line
                  line = remove_preceding_blanks(line)
!
                  if (trim(line) == 'method:') then ! Set method
!
                     read(unit_input,'(a40)') line
                     line = remove_preceding_blanks(line)
!
!                    Set method
!
                     method = trim(line)
!
                     exit
!
                  elseif (trim(line) == '}') then
!
                     write(unit_output,*)'Error: method was not specified in eT.inp.'
                     stop
!
                  endif 
!
               enddo ! End general do loop
!
            else 
!
               write(unit_output,*)'Error: method was not specified in eT.inp.'
               stop
!
            endif
            exit
!
         elseif (trim(line) == 'SCC') then
!
!           Go to SCC specific reader
!
            write(unit_output, '(t3,a)') 'Error: similarity constrained CC calculation requested, but it has not been implemented'
            stop
!
         elseif (trim(line) == 'end of eT input') then
!
            write(unit_output,*)'Error: method was not specified in eT.inp.'
            stop
!
         endif
!
      enddo ! End general do loop
!
   end subroutine method_reader
!
!
   function remove_preceding_blanks(line)
!!
!!    Remove Preceding Blanks. 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Nov 2017
!!    
!!    Removes white spaces before text from line
!!
      implicit none
!
      character(len=40) :: line
!
      character(len=40) :: remove_preceding_blanks
!
      integer(i15) :: i = 0, j = 0, length = 0
!
      do i = 1, 40
         if (line(i:i) == ' ') then
            continue
         else
            length = 40 - (i - 1)
            remove_preceding_blanks(1:length) = line(i:40)
            remove_preceding_blanks(length+1:40) = ' '
            return
         endif
      enddo
!
   end function remove_preceding_blanks
!
!
end module input_reader
