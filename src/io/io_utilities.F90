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
!!  
!!     Remove preceding blanks
!!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Nov 2017
!!  
!!     Removes white spaces before text from line
!!  
      implicit none
!
      character(len=100) :: line
!
      character(len=100) :: remove_preceding_blanks
!
      integer(i15) :: i = 0, length = 0
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
         if (trim(line) == 'end ' // string) then
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
   subroutine move_to_section(string, n_records)
!!
!!    Move to section,
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    Moves cursor to section given by string, and counts the number of records in the section
!!
      implicit none
!
      character(len=*), intent(in) :: string
!
      integer(i15) :: n_records
!
      character(len=100) :: line
!
      integer(i15) :: n_start, count_start, count_rec_end, count_rec_start, i
!
      rewind(input%unit)
!
      count_rec_end = 0
      n_start = 0
!
      count_start = 0
      count_rec_start = 0
!
      do 
!
         read(input%unit, '(a100)') line
         line = remove_preceding_blanks(line)
!
         count_rec_end = count_rec_end + 1
!
         if (trim(line) .eq. 'geometry') then
!
            backspace(input%unit)
            call output%error_msg('could not move to requested section: '// string)
!
         endif
!
         if (trim(line) == 'end ' // string) then        
!
            rewind(input%unit)
!
            do i = 1, count_rec_end - 1
!
               read(input%unit, '(a100)') line
               line = remove_preceding_blanks(line) 
!  
               if (trim(line) == string) then
!
                  n_start =  n_start + 1
!
               endif
!
            enddo 
!
            rewind(input%unit)
!
            do i = 1, count_rec_end - 1
!
              read(input%unit, '(a100)') line
              line = remove_preceding_blanks(line)
!  
              count_rec_start = count_rec_start + 1
!  
               if (trim(line) == string) then
!
                  count_start = count_start + 1
!
                  if (count_start == n_start) then
!
                     n_records = count_rec_end - count_rec_start - 1
!
                     return
!
                  endif
!
               endif
!
            enddo    
! 
         endif
!
      enddo
!
   end subroutine move_to_section
!
!
end module
