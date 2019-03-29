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
module input_file_class
!
!!
!!    Input file class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
!
   use kinds   
   use section_class 
   use output_file_class
   use string_utilities  
!
   type, extends(abstract_file) :: input_file 
!
      type(section) :: cholesky_section
!
   contains
!
      procedure :: init => init_input_file
!
      procedure, nopass :: string_is_comment    => string_is_comment_input_file
!
      procedure :: check_for_illegal_keywords   => check_for_illegal_keywords_input_file
      procedure :: check_section_for_illegal_keywords   => check_section_for_illegal_keywords_input_file
!
      procedure, nopass :: extract_keyword_from_string  => extract_keyword_from_string_input_file
!
      procedure :: requested_section            => requested_section_input_file
      procedure :: requested_keyword_in_section => requested_keyword_in_section_input_file
!
      generic :: get_keyword_in_section => get_integer_keyword_in_section_input_file,   &
                                            get_string_keyword_in_section_input_file,    &
                                            get_dp_keyword_in_section_input_file
!
      generic :: get_required_keyword_in_section => get_required_string_keyword_in_section_input_file, &
                                                     get_required_integer_keyword_in_section_input_file, &
                                                     get_required_dp_keyword_in_section_input_file
!
      procedure :: get_integer_keyword_in_section_input_file
      procedure :: get_string_keyword_in_section_input_file
      procedure :: get_dp_keyword_in_section_input_file
      procedure :: get_required_string_keyword_in_section_input_file
      procedure :: get_required_integer_keyword_in_section_input_file
      procedure :: get_required_dp_keyword_in_section_input_file
!
      procedure, private :: get_string_keyword_in_section_wo_safety => get_string_keyword_in_section_wo_safety_input_file
!
      procedure :: move_to_section  => move_to_section_input_file
!
      procedure :: get_n_elements_for_keyword_in_section => get_n_elements_for_keyword_in_section_input_file
      procedure :: get_array_for_keyword_in_section      => get_array_for_keyword_in_section_input_file
!
   end type input_file
!
!
   type(input_file) :: input
!
!
contains
!
!
   subroutine init_input_file(the_file, name)
!!
!!    Initialize input file
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2018
!!
!!    Initializes an input file file object
!!
!!    Output file is formatted and sequantial file.
!!    Routine sets these, and sets tha file name    
!!
      implicit none
!
      class(input_file) :: the_file
!
      character(len=*) :: name
!
      the_file%name = name
!
      the_file%access = 'sequential'
      the_file%format = 'formatted'
!
      input%cholesky_section%name_ = 'cholesky'
!
      input%cholesky_section%keywords = [ 'threshold           ',    &
                                          'span                ',    &
                                          'batches             ',    &
                                          'qualified           ',    &
                                          'one center          ',    &
                                          'no vectors          '     ]
!
   !   input%ccgs_section%name_ = 'cc ground state'
!
   !   input%ccgs_section%keywords = [ '',    &]
!
   end subroutine init_input_file
!
!
   subroutine check_for_illegal_keywords_input_file(the_file)
!!
!!    Check for illegal keywords 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
!!    Checks each section in turn, stopping with an error if it finds a keyword 
!!    that is not recognized.
!!
      implicit none 
!
      class(input_file) :: the_file
!
      call the_file%check_section_for_illegal_keywords(the_file%cholesky_section)
!
   end subroutine check_for_illegal_keywords_input_file
!
!
   subroutine check_section_for_illegal_keywords_input_file(the_file, the_section)
!!
!!    Checks a particular section for illegal keywords 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
      implicit none 
!
      class(input_file) :: the_file
!
      class(section) :: the_section
!
      character(len=200) :: string, keyword 
!
      logical :: recognized 
!
      integer :: n_records, record, k
!
      if (the_file%requested_section(the_section%name_)) then 
!
         call the_file%move_to_section(the_section%name_, n_records)
!
         do record = 1, n_records
!
            recognized = .false.
!
            read(the_file%unit, '(a200)') string 
!
            if (.not. input%string_is_comment(string)) then 
!
               call input%extract_keyword_from_string(string, keyword)
!
               do k = 1, size(the_section%keywords)
!
                  if (trim(the_section%keywords(k)) == trim(keyword)) recognized = .true. 
!
               enddo
!
               if (.not. recognized) then 
!
                  write(output%unit, '(/t3,a,a,a,a,a)') 'Could not recognize keyword "', trim(keyword), &
                                                         '" in section "', trim(the_section%name_), '".'
!
                  call the_section%print_keywords()
!
                  call output%error_msg('Something is wrong in the input file. See above.')
!
               endif 
!
            endif 
!
         enddo
!
      endif 
!
   end subroutine check_section_for_illegal_keywords_input_file
!
!
   subroutine get_integer_keyword_in_section_input_file(the_file, keyword, section, keyword_value)
!!
!!    Read integer keyword in section 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
!!    If specified, reads keyword as an integer into keyword value.
!!
!!    Note: if the keyword is not present in the section, keyword_value will not be set,
!!    and no error occurs. In typical usage, the standard value of the keyword is set before
!!    this routine is called. Changes from standard are made only when the keyword is specified
!!    - hence no errors when it does not find the keyword. 
!!  
      implicit none 
!
      class(input_file), intent(in) :: the_file
!
      character(len=*), intent(in) :: keyword 
      character(len=*), intent(in) :: section  
!
      integer, intent(out) :: keyword_value 
!
      character(len=200) :: keyword_value_string
!
      if (the_file%requested_section(section)) then
!
         if (the_file%requested_keyword_in_section(keyword, section)) then 
!
!           Get the keyword value in string format 
!
            call the_file%get_string_keyword_in_section_wo_safety(keyword, section, keyword_value_string)
!
!           Extract the integer from the string
!
            read(keyword_value_string, *) keyword_value
!
         endif
!
      endif
!
   end subroutine get_integer_keyword_in_section_input_file
!
!
   subroutine get_required_integer_keyword_in_section_input_file(the_file, keyword, section, keyword_value)
!!
!!    Read requested integer keyword in section 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
!!    If keyword and section not specified
!!    ends with error because required keyword 
!!    not provided.
!!
!!  
      implicit none 
!
      class(input_file), intent(in) :: the_file
!
      character(len=*), intent(in) :: keyword 
      character(len=*), intent(in) :: section  
!
      integer, intent(out) :: keyword_value 
!
      character(len=200) :: keyword_value_string
!
      if (.not. the_file%requested_section(section)) & 
         call output%error_msg('could not find the required section: '// trim(section) // '.')
!
      if (.not. the_file%requested_keyword_in_section(keyword, section)) & 
         call output%error_msg('could not find the required keyword '// trim(keyword) // ' in section ' // trim(section))
!
!     Get the keyword value in string format 
!
      call the_file%get_string_keyword_in_section_wo_safety(keyword, section, keyword_value_string)
!
!     Extract the integer from the string
!
      read(keyword_value_string, *) keyword_value
!
   end subroutine get_required_integer_keyword_in_section_input_file
!
!
   subroutine get_dp_keyword_in_section_input_file(the_file, keyword, section, keyword_value)
!!
!!    Read double precision keyword in section 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
!!    If specified, reads keyword as a real double precision into keyword value.
!!
!!    Note: if the keyword is not present in the section, keyword_value will not be set,
!!    and no error occurs. In typical usage, the standard value of the keyword is set before
!!    this routine is called. Changes from standard are made only when the keyword is specified
!!    - hence no errors when it does not find the keyword. 
!!    
      implicit none 
!
      class(input_file), intent(in) :: the_file
!
      character(len=*), intent(in) :: keyword 
      character(len=*), intent(in) :: section  
!
      real(dp), intent(out) :: keyword_value 
!
      character(len=200) :: keyword_value_string
!
      if (the_file%requested_section(section)) then
!
         if (the_file%requested_keyword_in_section(keyword, section)) then 
!
!           Get the keyword value in string format 
!
            call the_file%get_string_keyword_in_section_wo_safety(keyword, section, keyword_value_string)
!
!           Extract the integer from the string
!
            read(keyword_value_string, *) keyword_value
!
         endif 
!
      endif
!
   end subroutine get_dp_keyword_in_section_input_file
!
!
   subroutine get_required_dp_keyword_in_section_input_file(the_file, keyword, section, keyword_value)
!!
!!    Read requested double precision keyword in section 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
!!    If specified, reads keyword as a real double precision into keyword value.
!!
!!    Note: if the keyword is not present in the section, keyword_value will not be set,
!!    and no error occurs. In typical usage, the standard value of the keyword is set before
!!    this routine is called. Changes from standard are made only when the keyword is specified
!!    - hence no errors when it does not find the keyword. 
!!    
      implicit none 
!
      class(input_file), intent(in) :: the_file
!
      character(len=*), intent(in) :: keyword 
      character(len=*), intent(in) :: section  
!
      real(dp), intent(out) :: keyword_value 
!
      character(len=200) :: keyword_value_string
!
      if (.not. the_file%requested_section(section)) & 
         call output%error_msg('could not find the required section: '// trim(section) // '.')
!
      if (.not. the_file%requested_keyword_in_section(keyword, section)) & 
         call output%error_msg('could not find the required keyword '// trim(keyword) // ' in section ' // trim(section))
!
!     Get the keyword value in string format 
!
      call the_file%get_string_keyword_in_section_wo_safety(keyword, section, keyword_value_string)
!
!     Extract the integer from the string
!
      read(keyword_value_string, *) keyword_value
!
   end subroutine get_required_dp_keyword_in_section_input_file
!
!
   subroutine get_string_keyword_in_section_input_file(the_file, keyword, section, keyword_value)
!!
!!    Read string keyword in section 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
!!    If specified, reads keyword as a string into keyword value.
!!
!!    Note: if the keyword is not present in the section, keyword_value will not be set,
!!    and no error occurs. In typical usage, the standard value of the keyword is set before
!!    this routine is called. Changes from standard are made only when the keyword is specified
!!    - hence no errors when it does not find the keyword. 
!!    
      implicit none 
!
      class(input_file), intent(in) :: the_file
!
      character(len=*), intent(in) :: keyword 
      character(len=*), intent(in) :: section  
!
      character(len=200) :: keyword_value 
!
      if (the_file%requested_section(section)) then
!
         if (the_file%requested_keyword_in_section(keyword, section)) then 
!
!           Get the keyword value in string format 
!
            call the_file%get_string_keyword_in_section_wo_safety(keyword, section, keyword_value)
!
         endif
!
      endif 
!
   end subroutine get_string_keyword_in_section_input_file
!
!
   subroutine get_required_string_keyword_in_section_input_file(the_file, keyword, section, keyword_value)
!!
!!    Read string keyword in section 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
!!    If specified, reads keyword as a string into keyword value.
!!
!!    Note: if the keyword is not present in the section, keyword_value will not be set,
!!    and no error occurs. In typical usage, the standard value of the keyword is set before
!!    this routine is called. Changes from standard are made only when the keyword is specified
!!    - hence no errors when it does not find the keyword. 
!!    
      implicit none 
!
      class(input_file), intent(in) :: the_file
!
      character(len=*), intent(in) :: keyword 
      character(len=*), intent(in) :: section  
!
      character(len=200) :: keyword_value 
!
      if (.not. the_file%requested_section(section)) & 
         call output%error_msg('could not find the required section: '// trim(section) // '.')
!
      if (.not. the_file%requested_keyword_in_section(keyword, section)) & 
         call output%error_msg('could not find the required keyword '// trim(keyword) // ' in section ' // trim(section))
!
!     Get the keyword value in string format 
!
      call the_file%get_string_keyword_in_section_wo_safety(keyword, section, keyword_value)
!
   end subroutine get_required_string_keyword_in_section_input_file
!
!
   subroutine get_string_keyword_in_section_wo_safety_input_file(the_file, keyword, section, keyword_value)
!!
!!    Read string keyword in section without safety 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
!!    Reads keyword in section, placing the result in the string keyword_value. This 
!!    routine gives an error if the keyword is not located. 
!!
      implicit none 
!
      class(input_file), intent(in) :: the_file
!
      character(len=*), intent(in) :: keyword 
      character(len=*), intent(in) :: section  
!
      character(len=200) :: keyword_value 
!
      integer :: n_records, record, len_line_keyword
!
      character(len=200) :: line
!
!     Move to the requested section & get the number of records in that section 
!
      call the_file%move_to_section(section, n_records)
!
!     Loop through records within the section to locate & get the keyword value 
!
      do record = 1, n_records
!
         read(the_file%unit, '(a200)') line 
!
         line = adjustl(line) 
!
         if (line(1 : 1) == '!') cycle ! Comment
!
         len_line_keyword = 0
         do while (line(len_line_keyword + 1 : len_line_keyword + 1) /= ':' .and. len_line_keyword < 199)
!
            len_line_keyword = len_line_keyword + 1
!
         enddo
!
         if (len_line_keyword == 0) then 
!
            call output%error_msg('Failed to read keyword ' // keyword // ' in section ' // section)
!
         endif 
!
         if (trim(line(1 : len_line_keyword)) == keyword) then 
!
            keyword_value = adjustl(line(len_line_keyword + 2 : 200))
            return
!
         endif 
!
      enddo 
!
!     If you are here, you have not returned, so you have not found the keyword! 
!
      call output%error_msg('Failed to read keyword ' // keyword // ' in section ' // section)
!
   end subroutine get_string_keyword_in_section_wo_safety_input_file
!
!
   logical function string_is_comment_input_file(string)
!!
!!    Is the string a comment?
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
      implicit none 
!
      character(len=200), intent(in) :: string
!
      character(len=200) :: tmp_string 
!
      tmp_string = adjustl(string)
!
      if (string(1:1) == '!') then 
!
         string_is_comment_input_file = .true.
!
      else
!
         string_is_comment_input_file = .false.
!
      endif
!
   end function string_is_comment_input_file
!
!
   subroutine extract_keyword_from_string_input_file(string, keyword)
!!
!!    Extract keyword from string 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!    
!!    Note: assumes that the string is not a comment. This routine is therefore 
!!    called after a call to "is comment" logical routine. 
!!
      implicit none 
!
      character(len=200), intent(in)   :: string 
      character(len=200), intent(out)  :: keyword 
!
      integer :: k, colon_position
!
      keyword = adjustl(string)
!
!     If there is a ":", we have to remove the ":" as well as the value(s) that follow 
!
      colon_position = -1 
!
      do k = 1, 200
!
         if (keyword(k : k) == ':') then 
!
            colon_position = k 
!
         endif 
!
      enddo
!
      if (colon_position .ne. -1) then 
!
         do k = colon_position, 200
!
            keyword(k : k) = ' '
!
         enddo
!
      endif 
!
   end subroutine extract_keyword_from_string_input_file
!
!
   logical function requested_keyword_in_section_input_file(the_file, keyword, section)
!!
!!    Is string keyword in section?
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
!!    Returns true if the keyword is in the section. 
!!
!!    Note: stops inside "move to section" if the section is not present!
!!    => this routine should only be called if you know the section in present.
!!    In typical usage, it is called after we have made sure the section exists; 
!!    see the "section exists" function.
!!
      implicit none 
!
      class(input_file), intent(in) :: the_file
!
      character(len=*), intent(in) :: keyword 
      character(len=*), intent(in) :: section  
!
      integer :: n_records, record, len_line_keyword
!
      character(len=200) :: line
!
!     Move to the requested section & get the number of records in that section 
!
      call the_file%move_to_section(section, n_records)
!
!     Loop through records within the section & try to locate the keyword 
!
      do record = 1, n_records
!
         read(the_file%unit, '(a200)') line 
!
         line = adjustl(line) 
!
         if (line(1 : 1) == '!') cycle ! Comment
!
         len_line_keyword = 0
         do while (line(len_line_keyword + 1 : len_line_keyword + 1) /= ':' .and. len_line_keyword < 199)
!
            len_line_keyword = len_line_keyword + 1
!
         enddo
!
         if (len_line_keyword == 0) then 
!
            call output%error_msg('Failed to read keyword ' // keyword // ' in section ' // section)
!
         endif 
!
         if (trim(line(1 : len_line_keyword)) == keyword) then 
!
            requested_keyword_in_section_input_file = .true.
            return
!
         endif 
!
      enddo 
!
!     If you are here, you have not returned, so you have not found the keyword! 
!
      requested_keyword_in_section_input_file = .false.
!
   end function requested_keyword_in_section_input_file
!
!
   logical function requested_section_input_file(the_file, section)
!!
!!    Does section exist? 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
!!    Returns true if the section exists, false if it doesn't 
!!
      implicit none 
!
      class(input_file), intent(in) :: the_file 
!
      character(len=*), intent(in) :: section 
!
      character(len=200) :: line 
!
      requested_section_input_file = .false.
!
      rewind(the_file%unit)
!
      do 
!
         read(the_file%unit, '(a200)') line
         line = adjustl(line)
!
         if (trim(line) == section) then 
!
            requested_section_input_file = .true.
            return 
!
         endif 
!
         if (trim(line) == 'geometry') then 
!
            requested_section_input_file = .false.
            return 
!
         endif
!
      enddo 
!
   end function requested_section_input_file
!
!
   subroutine move_to_section_input_file(the_file, string, n_records)
!!
!!    Move to section
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 / Mar 2019
!!
!!    Moves cursor to section given by string, and 
!!    counts & returns the number of records in that section.
!!
      implicit none
!
      class(input_file), intent(in) :: the_file
!
      character(len=*), intent(in) :: string
!
      integer :: n_records
!
      character(len=200) :: line
!
      integer :: n_start, count_start, count_rec_end, count_rec_start, i
!
      rewind(the_file%unit)
!
      count_rec_end = 0
      n_start = 0
!
      count_start = 0
      count_rec_start = 0
!
      do 
!
         read(the_file%unit, '(a200)') line
         line = adjustl(line)
!
         count_rec_end = count_rec_end + 1
!
         if (trim(line) .eq. 'geometry') then
!
            backspace(the_file%unit)
            call output%error_msg('could not move to requested section: '// string)
!
         endif
!
         if (trim(line) == 'end ' // string) then        
!
            rewind(the_file%unit)
!
            do i = 1, count_rec_end - 1
!
               read(the_file%unit, '(a200)') line
               line = adjustl(line) 
!  
               if (trim(line) == string) then
!
                  n_start =  n_start + 1
!
               endif
!
            enddo 
!
            rewind(the_file%unit)
!
            do i = 1, count_rec_end - 1
!
              read(the_file%unit, '(a200)') line
              line = adjustl(line)
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
   end subroutine move_to_section_input_file
!
!
   function get_n_elements_for_keyword_in_section_input_file(the_file, keyword, section) result(n_elements)
!!
!!    Get n elements for keyword in section 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
!!    Gets the number of elements of input variable array
!!    for keyword which is specified on input by either an
!!    integer range or list (of length n_elements).
!!
!!    Ranges should always be given as [a,b].
!!
!!    Lists should always be given as {a, b, c, d},
!!    that is, in set notation.
!!
!!    Function is called in preparation for
!!    get_array_for_keyword_in_section
!!
      implicit none 
!
      class(input_file), intent(in) :: the_file
!
      character(len=*), intent(in) :: keyword 
      character(len=*), intent(in) :: section  
!
      integer :: n_elements
!
      character(len=200) :: keyword_value_string
!
!     Get the keyword value in string format 
!
      call the_file%get_keyword_in_section(keyword, section, keyword_value_string)
!
!     Use string utility functionality to get n_elements   
!  
      n_elements = get_n_elements_in_string(keyword_value_string)
!
   end function get_n_elements_for_keyword_in_section_input_file
!
!
   subroutine get_array_for_keyword_in_section_input_file(the_file, keyword, section, n_elements, array_)
!!
!!    Get array for keyword in section 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
!!    Gets input variable array (array_) for keyword
!!    which is specified on input by either an
!!    integer range or list (of length n_elements).
!!
!!    Ranges should always be given as [a,b].
!!
!!    Lists should always be given as {a, b, c, d},
!!    that is, in set notation.
!!
!!    Routine should be called after the 
!!    get_n_elements_for_keyword_in_section is called
!!    in order to determine n_elements so that array_ 
!!    can be allocated.
!!
!!
      implicit none 
!
      class(input_file), intent(in) :: the_file
!
      character(len=*), intent(in) :: keyword 
      character(len=*), intent(in) :: section  
!
      integer, intent(in) :: n_elements
!
      integer, dimension(n_elements) :: array_
!
      character(len=200) :: keyword_value_string
!
!     Get the keyword value in string format 
!
      call the_file%get_keyword_in_section(keyword, section, keyword_value_string)
!
!     Use string utility functionality to get the array
!  
      call get_elements_in_string(keyword_value_string, n_elements, array_)
!
   end subroutine get_array_for_keyword_in_section_input_file
!
!
end module input_file_class
