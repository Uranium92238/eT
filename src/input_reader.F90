module input_reader
!
!                          Input reader module                                 
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017 
!
   use types
   use input_output
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
      integer(i15) :: position = 1, line_number = 0
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
               write(unit_output,*) 'Input error: could not identify method.'
               stop ! Terminate program
            exit
!
         else
!
            write(unit_output,*) 'Input error: "Method" section must come first in input file.'
            stop ! Terminate program
!
         endif
!
      enddo ! End general do loop
!
   end subroutine method_reader
!
!
end module input_reader