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
!
contains
!
!
   function set_cursor_to_substring(string, substring) result(cursor)
!!
!!    Set cursor to substring
!!    Written by Sarai D. Folkestad, Dec 2020
!!
!!    Sets cursor to the final character of the first occurrence
!!    of the provided substring.
!!
      implicit none
!
      character(len=*), intent(in) :: string
!
      character(len=*), intent(in) :: substring
!
      integer :: cursor, start_, substring_length
!
      substring_length = len(trim(substring))
!
      start_ = 1
      cursor = start_ + substring_length - 1
!
      do while (cursor .lt. len(string))
!
         if (string(start_ : cursor) == trim(substring)) then
!
            return
!
         else
!
            start_ = start_ + 1
            cursor = cursor + 1
!
         endif
      enddo
!
      print *, "Error in set_cursor_to_substring"
      error stop
!
   end function set_cursor_to_substring
!
!
   pure function n_instances_of_character(string, char) result(n)
!!
!!    Number of instances of character
!!    Written by Sarai D. Folkestad, Sep 2021
!!
      implicit none
!
      character(len=*), intent(in) :: string
!
      character(len=1), intent(in) :: char
!
      integer :: n, cursor
!
      cursor = 1
      n = 0
!
      do while (cursor .lt. len(string))
!
         if (string(cursor : cursor) == char) n = n + 1
!
         cursor = cursor + 1
!
      enddo
!
   end function n_instances_of_character
!
!
   function last_instance_of_character(string, char) result(cursor)
!!
!!    Last instance of character
!!    Written by Sarai D. Folkestad, Dec 2020
!!
!!    Sets cursor to the final character of the first occurrence
!!    of the provided substring.
!!
      implicit none
!
      character(len=*), intent(in) :: string
!
      character(len=1), intent(in) :: char
!
      integer :: n, cursor
!
      cursor = len(string)
      n = 0
!
      do while (cursor .gt. 0)
!
         if (string(cursor : cursor) == char) return
!
         cursor = cursor - 1
!
      enddo
!
      print *, "Error in last_instance_of_character"
      error stop
!
   end function last_instance_of_character
!
!
   function first_instance_of_character(string, char) result(cursor)
!!
!!    First instance of character
!!    Written by Sarai D. Folkestad, Dec 2020
!!
      implicit none
!
      character(len=*), intent(in) :: string
!
      character(len=1), intent(in) :: char
!
      integer :: n, cursor
!
      cursor = 1
      n = 0
!
      do while (cursor .gt. 0)
!
         if (string(cursor : cursor) == char) return
!
         cursor = cursor + 1
!
      enddo
!
      print *, "Error in first_instance_of_character"
      error stop
!
   end function first_instance_of_character
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
      if (len(string) .gt. 500) then
         print *, 'Error: string is too long to convert.'
         error stop
      endif
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
   function is_substring_in_string(string, substring) result(substring_found)
!!
!!    Substring in string
!!    Written by Sarai D. Folkestad, Aug 2021
!!
      implicit none
!
      character(len=*), intent(in) :: string
      character(len=*), intent(in) :: substring
!
      logical :: substring_found
!
      integer :: start_, substring_length, cursor
!
      substring_length = len(trim(substring))
!
      start_ = 1
      cursor = start_ + substring_length - 1
!
      do while (cursor .le. len(string))
!
         if (string(start_ : cursor) == trim(substring)) then
!
            substring_found = .true.
            return
!
         else
!
            start_ = start_ + 1
            cursor = cursor + 1
!
         endif
      enddo
!
      substring_found = .false.
!
   end function is_substring_in_string
!
!
   subroutine split_at_delimiter(string, n, array_of_strings, delimiter)
!!
!!    Split at delimiter
!!    Written by Sarai D. Folkestad, 2021
!!
      implicit none
!
      character(len=*), intent(in) :: string
!
      integer, intent(in) :: n
!
      character(len=200), dimension(n), intent(out) :: array_of_strings
!
      character, intent(in) :: delimiter
!
      integer :: string_length, cursor_2, cursor_1, counter
!
      character(len=1000) :: left_aligned_string
!
      string_length = len_trim(adjustl(string))
      if (string_length > 1000) stop "String too long for split_at_delimiter routine"
!
      left_aligned_string = adjustl(string)
!
      cursor_2 = 1
      cursor_1 = 1
!
      counter = 0
!
      do while (cursor_2 .le. string_length)
!
         if (trim(left_aligned_string(cursor_1:)) == '') return
!
         if (left_aligned_string(cursor_2:cursor_2) == delimiter) then
!
            if (cursor_1 == cursor_2) cycle
!
            counter =  counter + 1
            if (counter > n) stop "too many substrings for split_at_delimiter routine"
!
            array_of_strings(counter) = trim(adjustl(left_aligned_string(cursor_1:cursor_2 - 1)))
            cursor_1 = cursor_2 + 1
!
         endif
!
         cursor_2 = cursor_2 + 1
!
      enddo
!
   end subroutine split_at_delimiter
!
!
   subroutine remove_delimiter_from_string(string, delimiter)
!!
!!    Remove delimiter from string
!!    Written by Sarai D. Folkestad, Dec 2021
!!
      implicit none
!
      character(len=*), intent(inout) :: string
!
      character, intent(in) :: delimiter
!
      integer :: i
!
      do i = 1, len(string)
!
         if (string(i:i) == delimiter) string(i:) = string(i+1:)
!
      enddo
!
   end subroutine remove_delimiter_from_string
!
!
end module string_utilities
