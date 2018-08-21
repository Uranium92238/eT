module io_utilities
!
!!
!!    IO utilities module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkstad, June 2018
!!
!
   use kinds
   use file_class
!
contains
!
   function remove_preceding_blanks(line)
!  
!     Remove preceding blanks
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Nov 2017
!  
!     Removes white spaces before text from line
!  
      implicit none
!
      character(len=100) :: line
!
      character(len=100) :: remove_preceding_blanks
!
      integer(i15) :: i = 0, j = 0, length = 0
!
      remove_preceding_blanks = ' '
!
      do i = 1, 100
         if (line(i:i) == ' ') then
!
            continue
!
         else
!
            length = 100 - (i - 1)
            remove_preceding_blanks(1:length) = line(i:100)
            remove_preceding_blanks(length+1:100) = ' '
            return
!
         endif
      enddo
!
   end function remove_preceding_blanks
!
!
   logical function requested_section(string)
!!
!!
!!
      implicit none
!
      character(len=*) :: string
!
      character(len=100) :: line
!
      rewind(input%unit)
!
      requested_section = .false.
!
      do 
!
         read(input%unit, '(a100)') line
         line = remove_preceding_blanks(line)
!
         if (trim(line) .eq. 'end geometry') then
            backspace(input%unit)
            return
         endif
!
         if (trim(line) == string) then
!
            requested_section = .true.
            return
! 
         endif
!
      enddo
!
   end function requested_section
!
!
end module
