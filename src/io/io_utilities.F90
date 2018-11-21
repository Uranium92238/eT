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
   subroutine long_string_print(string,format_string,colons,fformat_string,lformat_string,line_length)
!!
!!    Long string print
!!    Written by Rolf H. Myhre, Nov 2018
!!
!!    Prints a reasonbly formatted long string over several lines
!!    format_string: optional format string
!!    fformat_string: optional format string for first printed line, will be used for single line if present
!!    lformat_string: optional format string for last printed line
!!    line_length: optional argument for length printed lines, routine adds extra length to not split up words
!!
!!
      implicit none
!
      character(len=*), intent(in) :: string
      character(len=*), intent(in), optional :: format_string,fformat_string,lformat_string
      integer(i15), intent(in), optional :: line_length
      logical, intent(in), optional :: colons
!
      character(90)  :: temp
      integer(i15)   :: l, l_left, lines, l_length
      integer(i15)   :: i,j, padd, printed 
      character(20)  :: fs,fstring,ffstring,lfstring
      logical        :: col 
!
!     Default line length
      l_length = 70
      if(present(line_length)) then
         l_length = line_length
      endif
!
!     Figure out the formatting
      fstring = '(t3,a)'
      if(present(format_string)) then
         fstring = format_string
      endif
!
      ffstring = fstring
      if(present(fformat_string)) then
         ffstring = fformat_string
      endif
!
      lfstring = fstring
      if(present(lformat_string)) then
         lfstring = lformat_string
      endif
!
!     First line format
      fs = ffstring
!
!     Fancy colons for banner headings
      col = .false.
      if(present(colons)) then
         col = colons
      endif
      if(col) then
         l_length = l_length - 3
      endif
!
      l = len_trim(string)      
      l_left = l
      lines = l/l_length + 1
      printed = 1
!
      do i = 1,lines
!
         if(i .ne. lines) then
!
            if(i .ne. 1) then
               fs = fstring
            endif
!
!           Add some extra padding to not split words
            do j = 1,18
               padd = j
               if(string(printed+l_length+j:printed+l_length+j) == ' ') then
                  exit
               endif
            enddo
!
!           Copy string to be printed. Add hyphen if word must be split
            if(padd == 18) then
               temp(1:l_length+padd+1) = string(printed:printed+l_length+padd)
               temp(l_length+padd+1:l_length+padd+1) = '-'
               printed = printed + l_length + padd
            else
               temp(1:l_length+padd+1) = string(printed:printed+l_length+padd+1)
               printed = printed + l_length + padd + 1
            endif
!
!           Print
            if(col) then
               write(output%unit, fs)  ':: '//temp(1:l_length+padd+1)
            else
               write(output%unit, fs)  temp(1:l_length+padd+1)
            endif
!
         else
!           Print the remaining string
            fs = lfstring
            if(col) then
               write(output%unit, fs)  ':: '//string(printed:l)
            else
               write(output%unit, fs)  string(printed:l)
            endif
            printed = l
         endif
!
      enddo
!
      flush(output%unit)
!
   end subroutine long_string_print
!
!
end module
