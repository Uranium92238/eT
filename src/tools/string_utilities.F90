!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2021 the authors of eT
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
module string_utilities
!
!!
!!    String utilities module
!!
!!    Routines that manipulate and analyze strings.
!!
!
   use kinds
   use global_out, only: output
!
contains
!
!
   function get_n_elements_in_string(string) result(n_elements)
!!
!!    Get n elements in string
!!    Written by Sarai D. Folkstad and Eirik F. Kjønstad, Mar 2019
!!
!!    Gets the number of elements in range or list,
!!    To be used for reading of input.
!!
!!    Ranges should always be given as [a,b].
!!
!!    Lists should always be given as {a, b, c, d},
!!    that is, in set notation.
!!
      implicit none
!
      character(len=200), intent(inout) :: string
!
      integer :: n_elements
!
!     Local variables
!
      integer :: first, last, n_characters
      integer :: i
!
      n_elements = 0
!
      string = adjustl(string)
!
      n_characters = len_trim(string)
!
      if (string(1:1) == '[') then ! range given
!
!        Sanity check - Is set closed?
!
         if (string(n_characters:n_characters) /= ']') call output%error_msg('found open range in input file.')
!
         do i = 2, n_characters - 1
!
            if (string(i:i) == ',') exit
!
         enddo
!
!        Read first element
!
         read(string(2:i-1), *) first
!
!        Read last element
!
         read(string(i+1:n_characters - 1), *) last
!
!        Calculate number of elements
!
         n_elements = last - first + 1
!
      elseif (string(1:1)=='{') then ! list given
!
!        Sanity check - Is set closed?
!
         if (string(n_characters:n_characters) /= '}') call output%error_msg('found open set in input file.')
!
         n_elements = 1 ! Assuming that the set contains at least one element (otherwize why give list?) 
!
!        Loop through and count commas
!
         do i = 2, n_characters - 1
!
            if (string(i:i) == ',') n_elements = n_elements + 1
!
         enddo
!
      else ! Did not find list or 
!
         n_elements = 0
!
      endif
!
   end function get_n_elements_in_string
!
!
   subroutine get_elements_in_string(string, n_elements, elements)
!!
!!    Get elements
!!    Written by Sarai D. Folkstad and Eirik F. Kjønstad, Mar 2019
!!
!!    Gets the elements from range or list.
!!    To be used for reading of input.
!!
!!    Ranges should always be given as [a,b].
!!
!!    Lists should always be given as {a, b, c, d},
!!    that is, in set notation.
!!
      implicit none
!
      character(len=200), intent(inout) :: string
!
      integer, intent(in) :: n_elements
!
      integer, dimension(n_elements), intent(out) :: elements
!
!     Local variables
!
      integer :: first, last, n_characters, n_elements_found
      integer :: i, j
!
      string = adjustl(string)
!
      n_characters = len_trim(string)
!
      if (string(1:1) == '[') then ! range given
!
!        Sanity check - Is set closed?
!
         if (string(n_characters:n_characters) /= ']') call output%error_msg('found open range in input file.')
!
         do i = 2, n_characters - 1
!
            if (string(i:i) == ',') exit
!
         enddo
!
!        Read first element
!
         read(string(2:i-1), *) first
!
!        Read last element
!
         read(string(i+1:n_characters - 1), *) last
!
!        Sanity check - Is the number of elements found equal to n_elements
!
         if (n_elements .ne. last - first + 1) call output%error_msg('Mismatch in number of elements to be read.')
!
         do i = 1, n_elements
!
            elements(i) = first + i - 1
!
         enddo
!
      elseif (string(1:1)=='{') then ! list given
!
!        Sanity check - Is set closed?
!
         if (string(n_characters:n_characters) /= '}') call output%error_msg('found open set in input file.')
!
!        Loop through and set the elements
!
         first            = 2
         n_elements_found = 0
!
         do j = 1, n_elements
!
            do i = first, n_characters - 1
!
               if (string(i:i) == ',') exit
!
            enddo
!
            read(string(first:i-1), *) elements(j)
!
            n_elements_found = n_elements_found + 1
!
            first = i + 1
!
            if (first == n_characters) exit
!
         enddo
!
         if (n_elements_found .ne. n_elements) call output%error_msg('Mismatch in number of elements to be read.')
!
      else ! Did not find list or 
!
         call output%error_msg('neither list nor range was found.')
!
      endif
!
   end subroutine get_elements_in_string
!
!
   subroutine get_reals_in_string(string, n_elements, elements)
!!
!!    Get reals
!!    Written by Sarai D. Folkstad and Eirik F. Kjønstad, Mar 2019
!!    Modified by Andreas Skeidsvoll, Sep 2019: Reads reals instead of integers
!!
!!    Gets the reals from list.
!!    To be used for reading of input.
!!
!!    Lists should always be given as {a, b, c, d},
!!    that is, in set notation.
!!
      implicit none
!
      character(len=200), intent(inout) :: string
!
      integer, intent(in) :: n_elements
!
      real(dp), dimension(n_elements), intent(out) :: elements
!
!     Local variables
!
      integer :: first, n_characters, n_elements_found
      integer :: i, j
!
      string = adjustl(string)
!
      n_characters = len_trim(string)
!
      if (string(1:1)=='{') then ! list given
!
!        Sanity check - Is set closed?
!
         if (string(n_characters:n_characters) /= '}') call output%error_msg('found open set in input file.')
!
!        Loop through and set the elements
!
         first            = 2
         n_elements_found = 0
!
         do j = 1, n_elements
!
            do i = first, n_characters - 1
!
               if (string(i:i) == ',') exit
!
            enddo
!
            read(string(first:i-1), *) elements(j)
!
            n_elements_found = n_elements_found + 1
!
            first = i + 1
!
            if (first == n_characters) exit
!
         enddo
!
         if (n_elements_found .ne. n_elements) call output%error_msg('Mismatch in number of elements to be read.')
!
      else ! Did not find list or 
!
         call output%error_msg('neither list nor range was found.')
!
      endif
!
   end subroutine get_reals_in_string
!
!
   function set_cursor_to_character(string,final_character) result(cursor)
!!
!!    Set cursor to string
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!    Modified by Tommaso Giovannini, May 2019
!!
!!    Sets the cursor to the final character in the string
!!    More general than previous white space
!!
!!    final_character :: optional 
!!    if present --> look for final_character
!!    otherwise  --> look for whitespace
!!
      implicit none
!
      character(len=200), intent(inout) :: string
!
      character(len=1), intent(in), optional :: final_character
!
      integer :: cursor
!!
!!    internal variables
!!
      character(len=1) :: look_character
!      
      string = adjustl(string)
!
      if (present(final_character)) then
!      
         look_character = final_character
!         
      else 
!      
         look_character = ' '
!         
      endif
!
      cursor = 1
!
      do while (cursor .lt. 200)
         if (string(cursor:cursor) .eq. look_character) then
            exit
         else
            cursor = cursor + 1
            cycle
         endif
      enddo
!
   end function set_cursor_to_character
!
!
   subroutine convert_to_lowercase(string)
!!
!!    Convert to lowercase 
!!    Written by Eirik F. Kjønstad, Mar 2019 
!!
!!    Adapted from routines posted on the Stack-exchange.
!!    
!!    Assumes ASCII table for representing characters as integers, 
!!    where the lowercase letter is +32 relative to the uppercase letters.
!!
!!    Note: uppercase (65-90) and lowercase (97-122).
!!
      implicit none 
!
      character(len=*), intent(inout) :: string 
!
      integer :: character, current_character
!
      do character = 1, len(string)
!
!        Represent character as integer 
!  
         current_character = ichar(string(character : character)) 
!
!        Convert if character is in the range of uppercase characters 
!
         if (current_character >= 65 .and. current_character <= 90) then ! Between A and Z
!
            current_character = current_character + 32 
!
         endif
!
!        Replace the character by the (possibly) lowercased letter 
!
         string(character : character) = char(current_character)
!
      enddo
!
   end subroutine convert_to_lowercase
!
!
   function convert_char_to_uppercase(char_) result(char_out)
!!
!!    Convert to uppercase
!!    Written by Eirik F. Kjønstad, Mar 2019 
!!
!!    Adapted from routines posted on the Stack-exchange.
!!    
!!    Assumes ASCII table for representing characters as integers, 
!!    where the lowercase letter is +32 relative to the uppercase letters.
!!
!!    Note: uppercase (65-90) and lowercase (97-122).
!!
      implicit none 
!
      character, intent(in) :: char_ 
      character :: char_out 
!
      integer :: char_int
!
!     Represent character as integer 
!  
      char_int = ichar(char_) 
!
!     Convert if character is in the range of uppercase characters 
!
      if (char_int >= 97 .and. char_int <= 122) then ! Between a and z
!
         char_int = char_int - 32 
!
      endif
!
!     Replace the character by the (possibly) uppercased letter 
!
     char_out = char(char_int)
!
   end function convert_char_to_uppercase
!
!
   function convert_to_uppercase(string) result(string_out)
!!
!!    Convert to lowercase 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
!!    Adapted from routines posted on the Stack-exchange.
!!    
!!    Assumes ASCII table for representing characters as integers, 
!!    where the lowercase letter is +32 relative to the uppercase letters.
!!
!!    Note: uppercase (65-90) and lowercase (97-122).
!!
      implicit none 
!
      character(len=*), intent(in) :: string
      character(len=500) :: string_out 
!
      integer :: character, current_character
!
      if (len(string) .gt. 500) call output%error_msg('string is too long to convert.')
!
      string_out = repeat(' ', 500)
!
      do character = 1, len(string)
!
!        Represent character as integer 
!  
         current_character = ichar(string(character : character)) 
!
!        Convert if character is in the range of uppercase characters 
!
         if (current_character >= 97 .and. current_character <= 122) then ! Between a and z
!
            current_character = current_character - 32 
!
         endif
!
!        Replace the character by the (possibly) lowercased letter 
!
         string_out(character : character) = char(current_character)
!
      enddo
!
   end function convert_to_uppercase
!
!
   subroutine remove_spaces_etc(string)
!!
!!    Remove spaces etc.
!!    written by Rolf H. Myhre, Mar. 2020
!!
!!    Replaces characters that are not suitable for file names from string
!!    Trailing blanks are not touched
!!
!!    ' ' -> '_'
!!    '*' -> 's'
!!    '+' -> 'p'
!!
      implicit none
!
      character(len=*), intent(inout) :: string
!
      integer :: i
!
      do i = 1, len_trim(string)
!
         if (string(i:i) .eq. ' ') then 
            string(i:i) = '_'
         else if (string(i:i) .eq. '*') then
            string(i:i) = 's'
         else if (string(i:i) .eq. '+') then
            string(i:i) = 'p'
         endif
!
      enddo
!
   end subroutine remove_spaces_etc
!
!
   subroutine index_of_unique_strings(indices, n_strings, strings)
!!
!!    Index of unique strings
!!    Written by Tor S. Haugland, Oct 2019
!!
!!    Find index of unique strings in an array of strings.
!!
!!    indices:    on input an array of size n_strings.
!!                   on output the index of unique elements in strings. 0 elsewhere.
!!                   n_unique = count( indices /= 0 )
!!    n_strings:  check first `n_strings` values of strings
!!    strings:    array of strings to check for unique. 
!!
!!    Example: strings = ['a', 'b', 'a', 'c'], n_strings = 4
!!             will return indices = [1, 2, 4, 0]
!!    
      implicit none 
!
      integer,                                intent(in)   :: n_strings
      integer, dimension(n_strings),          intent(out)  :: indices
      character(len=*), dimension(n_strings), intent(in)   :: strings
!
      character(len=len(strings(1))), dimension(n_strings) :: unique_strings
      integer :: i, n_unique
!
      n_unique = 0
!
      indices(:) = 0
!
      do i = 1, n_strings
!
!        Try to add string to unique_strings. if failed: check next
!
         if (.not. any(unique_strings(1:n_unique) == strings(i) )) then
!
            n_unique = n_unique + 1
!
            unique_strings(n_unique) = strings(i)
            indices(n_unique) = i
!
         endif
!
      enddo
!
   end subroutine index_of_unique_strings
!
!
end module string_utilities
