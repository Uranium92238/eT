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
!!
!
   use kinds    
   use output_file_class      
!
   type, extends(abstract_file) :: input_file 
!
!
   contains
!
      procedure :: init                      => init_input_file
!
      procedure :: section_exists     => section_exists_input_file
!
      procedure :: keyword_is_in_section => keyword_is_in_section_input_file
!
      generic :: read_keyword_in_section => read_integer_keyword_in_section_input_file, &
                                             read_string_keyword_in_section_input_file, &
                                             read_dp_keyword_in_section_input_file
!
      procedure :: read_integer_keyword_in_section_input_file
      procedure :: read_string_keyword_in_section_input_file
      procedure :: read_dp_keyword_in_section_input_file
!
      procedure :: move_to_section           => move_to_section_input_file
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
   end subroutine init_input_file
!
!
   subroutine read_integer_keyword_in_section_input_file(the_file, keyword, section, keyword_value)
!!
!!    Read integer keyword in section 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
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
      if (the_file%keyword_is_in_section(keyword, section)) then 
!
!        Get the keyword value in string format 
!
         call the_file%read_keyword_in_section(keyword, section, keyword_value_string, .false.)
!
!        Extract the integer from the string
!
         read(keyword_value_string, *) keyword_value
!
      endif 
!
   end subroutine read_integer_keyword_in_section_input_file
!
!
   subroutine read_dp_keyword_in_section_input_file(the_file, keyword, section, keyword_value)
!!
!!    Read double precision keyword in section 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
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
      if (the_file%keyword_is_in_section(keyword, section)) then 
!
!        Get the keyword value in string format 
!
         call the_file%read_keyword_in_section(keyword, section, keyword_value_string, .false.)
!
!        Extract the integer from the string
!
         read(keyword_value_string, *) keyword_value
!
      endif 
!
   end subroutine read_dp_keyword_in_section_input_file
!
!
   subroutine read_string_keyword_in_section_input_file(the_file, keyword, section, keyword_value, safeguard)
!!
!!    Read string keyword in section 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
      implicit none 
!
      class(input_file), intent(in) :: the_file
!
      character(len=*), intent(in) :: keyword 
      character(len=*), intent(in) :: section  
!
      logical, optional, intent(in) :: safeguard
!
      logical :: local_safeguard, look_for_keyword
!
      character(len=200) :: keyword_value 
!
      integer :: n_records, record, len_line_keyword
!
      character(len=200) :: line
!
      local_safeguard = .true. 
!
      if (present(safeguard)) then 
!
         if (.not. safeguard) local_safeguard = .false.
!
      endif 
!
      if (local_safeguard) then 
!
         look_for_keyword = the_file%keyword_is_in_section(keyword, section)
!
      else 
!
         look_for_keyword = .true.
!
      endif
!
      if (look_for_keyword) then 
!
!        Move to the requested section & get the number of records in that section 
!
         call the_file%move_to_section(section, n_records)
!
!        Loop through records within the section to locate & get the keyword value 
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
!        If you are here, you have not returned, so you have not found the keyword! 
!
         call output%error_msg('Failed to read keyword ' // keyword // ' in section ' // section)
!
      endif 
!
   end subroutine read_string_keyword_in_section_input_file
!
!
   logical function keyword_is_in_section_input_file(the_file, keyword, section)
!!
!!    Is string keyword in section?
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
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
            keyword_is_in_section_input_file = .true.
            return
!
         endif 
!
      enddo 
!
!     If you are here, you have not returned, so you have not found the keyword! 
!
      keyword_is_in_section_input_file = .false.
!
   end function keyword_is_in_section_input_file
!
!
   logical function section_exists_input_file(the_file, section)
!!
!!    Does section exist? 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
      implicit none 
!
      class(input_file), intent(in) :: the_file 
!
      character(len=*), intent(in) :: section 
!
      character(len=200) :: line 
!
      section_exists_input_file = .false.
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
            section_exists_input_file = .true.
            return 
!
         endif 
!
         if (trim(line) == 'end geometry') then 
!
            section_exists_input_file = .false.
            return 
!
         endif
!
      enddo 
!
   end function section_exists_input_file
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
end module input_file_class
