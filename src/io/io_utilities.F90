!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2019 the authors of eT
!
!  eT is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  eT is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program. If not, see <https://www.gnu.org/licenses/>.
!
!
module io_utilities
!
!!
!!    IO utilities module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkstad, June 2018
!!
!
   use output_file_class
!
contains
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
      integer, intent(in), optional :: line_length
      logical, intent(in), optional :: colons
!
      character(90)  :: temp
      integer        :: l, l_left, lines, l_length
      integer        :: i,j, padd, printed 
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
