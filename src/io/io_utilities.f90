module io_utilities
!
!!
!!    IO utilities module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkstad, June 2018
!!
!
   use kinds
!
   contains
!
      function remove_preceding_blanks(line)
!     
!        Remove preceding blanks
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, Nov 2017
!     
!        Removes white spaces before text from line
!     
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
!
               continue
!
            else
!
               length = 40 - (i - 1)
               remove_preceding_blanks(1:length) = line(i:40)
               remove_preceding_blanks(length+1:40) = ' '
               return
!
            endif
         enddo
!
      end function remove_preceding_blanks
!
!
      function check_if_blank(line)
!     
!        Check if blank
!        Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2018
!     
!        Checks if line is blank
!     
         implicit none
!
         character(len=40) :: line
!
         logical :: check_if_blank
!
         check_if_blank = .false.
!
         if (line == '') then
             check_if_blank = .true.
             write(*,*)'halla'
         endif
!
      end function check_if_blank
!
   end module
!