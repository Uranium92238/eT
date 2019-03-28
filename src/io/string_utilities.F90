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
module string_utilities
!
!!
!!    String utilities module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkstad, Mar 2019
!!
!
   use kinds
   use file_class
!
contains
!
!
   function get_n_elements(string) result(n_elements)
!!
!!    Get n elements
!!    Written by Sarai D. Folkstad and Eirik F. Kjønstad, Mar 2019
!!
!!    Gets the number of elements in range or list,
!!    To be used for reading of input.
!!
!!    Ranges should always be given as [a,b].
!!
!!    Lists should always be given as {a, b, c, d}
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
         i = 2
!
         n_elements = 1 ! Assuming that the set contains at least one element (otherwize why give list?) 
!
!        Loop through and count commas
!
         do while (i < n_characters)
!
            if (string(i:i) == ',') n_elements = n_elements + 1
!
            i = i + 1
!
         enddo
!
      else ! Did not find list or 
!
         call output%error_msg('neither list nor range was found.')
!
      endif
!
   end function get_n_elements
!
!
   subroutine get_elements(string, n_elements, elements)
!!
!!    Get elements
!!    Written by Sarai D. Folkstad and Eirik F. Kjønstad, Mar 2019
!!
!!    Gets the elements from range or list.
!!    To be used for reading of input.
!!
!!    Ranges should always be given as [a,b].
!!
!!    Lists should always be given as {a, b, c, d}
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
         i = 2
         first = 2
!
         n_elements_found = 0
!
         do j = 1, n_elements
!
            do while (string(i:i) /= ',' .and. i < n_characters)
!
               i = i + 1
!
            enddo
!
            read(string(first:i-1), *) elements(j)
!
            n_elements_found = n_elements_found + 1
!
            i = i + 1
            first = i
!
            if (i == n_characters) exit
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
   end subroutine get_elements
!
!
end module string_utilities
