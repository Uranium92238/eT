module input_reader
!
!                          Input reader module                                 
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017 
!
   use types
   use input_output
   use calculation_procedures
!
contains
!
!
   subroutine method_reader(unit_input, method)
!
!     Method Reader
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!
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
         do while (line(1:1) == '!') ! Comment; read the next line
            read(unit_input,'(a40)') line 
         enddo
!
         if (trim(line)) == 'Method:') then
!
            read(unit_input,'(a40)') line
!
            do while (line(1:1) == '!') ! Comment; read the next line
               read(unit_input,'(a40)') line 
            enddo
!
            if (line(1:1) == '.') then 
               method = trim(line(2:40))
               exit ! Escape from general do loop
            else
               write(unit_output,*) 'Input error: could not read method.'
               stop ! Terminate program
            exit
!
         else
!
            write(unit_output,*) 'Input error: method section must come first in input file.'
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
!
!     Calculation Reader
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!
!     Reads the calculation from the input file and initializes the 
!     tasks requested of the wavefunction.
!
      implicit none 
!
      type(calculation_procedures) :: tasks 
!
!
      integer(i15), intent(in) :: unit_input
!
      character(len=40) :: method, calculation 
!
      character(len=40) :: line 
!
      do ! General do loop - ends when it reaches 'exit'
!
         read(unit_input,'(a40)') line 
!
         do while (line(1:1) == '!') ! Comment; read the next line
            read(unit_input,'(a40)') line 
         enddo
!
         if (trim(line)) == 'Calculation:') then
!
            read(unit_input,'(a40)') line
!
            do while (line(1:1) == '!') ! Comment; read the next line
               read(unit_input,'(a40)') line 
            enddo
!
            if (line(1:1) == '.') then 
               calculation = trim(line(2:40))
               ! Test for which type, set the logical in tasks, and cycle!
            else
               write(unit_output,*) 'Input error: could not read calculation.'
               stop ! Terminate program
            exit
!
         else
!
            write(unit_output,*) 'Input error: method section must come first in input file.'
            stop ! Terminate program
!
         endif
!
      enddo ! End general do loop
!
   end subroutine calculation_reader
!
!
end module input_reader