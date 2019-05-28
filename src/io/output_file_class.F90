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
module output_file_class
!
!!
!!    Output ile class module
!!    Written by Rolf H. Myhre, May 2019
!!
!!
!
   use kinds    
   use abstract_file_class, only : abstract_file 
!
   type, extends(abstract_file) :: output_file
!
!
   contains
!
      procedure :: open_                  => open_output_file
      procedure :: close_                 => close_output_file
      procedure :: flush_                 => flush_output_file
!
      procedure :: error_msg              => error_msg_output_file
      procedure :: warning_msg            => warning_msg_output_file
!
      procedure, public :: printf         => printf_output_file
      procedure, public :: printd         => printd_output_file
      procedure, public :: author         => author_output_file
!
      procedure :: long_string_print      => long_string_print_output_file
!
   end type output_file
!
   interface output_file
!
      procedure new_output_file
!
   end interface output_file
!
   type(output_file) :: output
!
contains
!
!
   module function new_output_file(name_) result(the_file)
!!
!!    Output file constructer
!!    Writen by Rolf H. Myhre May 2019
!!
!!    Output file is a formatted and sequantial file.
!!    Routine sets these, and sets the file name    
!!
      implicit none
!
      type(output_file) :: the_file
!
      character(len=*), intent(in) :: name_
!
      the_file%name_ = name_
!
      the_file%access_ = 'sequential'
      the_file%format_ = 'formatted'
!
   end function
!
!
   subroutine open_output_file(the_file)
!!
!!    Open the output file
!!    Written by Rolf Heilemann Myhre, May 2019
!!
      implicit none
!
      class(output_file) :: the_file
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      open(newunit=the_file%unit, file=the_file%name_, access=the_file%access_, &
           action='write', status='unknown', form=the_file%format_, iostat=io_error, iomsg=io_msg)
!
      if (io_error /= 0) stop 'Error: could not open eT output file '//trim(the_file%name_)//&
                             &'. Error message: '//trim(io_msg)
!
      the_file%opened = .true.
!
      call the_file%set_open_size()
!
   end subroutine open_output_file
!
!
   subroutine close_output_file(the_file)
!!
!!    Close the output file
!!    Written by Rolf Heilemann Myhre, May 2019
!!
      implicit none
!
      class(output_file) :: the_file
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      close(the_file%unit, iostat=io_error, iomsg=io_msg, status='keep')
!
      if (io_error.ne. 0) then
          stop 'Error: could not close eT output file '//trim(the_file%name_)//&
              &'error message: '//trim(io_msg)
      endif
!
      the_file%opened = .false.
!
   end subroutine close_output_file
!
!
   subroutine flush_output_file(the_file)
!!
!!    Flush the output file
!!    Written by Rolf Heilemann Myhre, May 2019
!!
      implicit none
!
      class(output_file) :: the_file
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      flush(the_file%unit, iostat=io_error, iomsg=io_msg)
!
      if (io_error /= 0) stop 'Error: could not flush eT output file '//trim(the_file%name_)//&
                             &'error message: '//trim(io_msg)
!
   end subroutine flush_output_file
!
!
   subroutine error_msg_output_file(out_file, error_specs, error_int)
!!
!!    Error message
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(output_file) :: out_file
!
      character(len=*) :: error_specs
!
      integer, optional :: error_int 
!
      character(len=40) :: error_int_char = ' '
!
      if (present(error_int)) then
!
         write(error_int_char, '(i12)') error_int   
!
         write(out_file%unit, '(/t3,a)') 'Error: ' // trim(error_specs) // ' ' // error_int_char
!
      else
!
         write(out_file%unit, '(/t3,a)') 'Error: ' // trim(error_specs)
!
      endif
!
      stop "Something went wrong, check the .out file"
!
   end subroutine error_msg_output_file
!
!
   subroutine warning_msg_output_file(out_file, warning_specs)
!!
!!    Warning message
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(output_file), intent(in) :: out_file
!
      character(len=*) :: warning_specs
!
      write(out_file%unit, '(/t3,a)') 'Warning: ' // trim(warning_specs)
!
   end subroutine warning_msg_output_file
!
!  
   subroutine printf_output_file(the_file, string, reals, ints, fs, ll)
!!
!!    Printf 
!!    Written by Rolf Heilemann Myhre, May 2019
!!
!!    Prints any number of reals and integers formatted Python style.
!!
!!    string:   String of character that should be printed, 
!!              including formatting of reals and integers 
!!    reals:    Array of real(dp) to print - in the order specified by string 
!!    ints:     Array of integers to print - in the order specified by string 
!!    fs:       Specifies the format of 'string', e.g. fs='(/t6,a)' gives 
!!              a new line, then indentation 6, then the value of 'string'
!!              with reals and integers as specified.
!!    ll:       Integer specifying number of characters per line of print.
!!
!!    Example: call output%printf('Energy (a.u.): (f18.12)', reals=[wf%energy], fs='(/t6,a)')
!!
      implicit none
!
      class(output_file), intent(in) :: the_file
!
      character(len=*), intent(in)                 :: string 
      real(dp), dimension(:), intent(in), optional :: reals 
      integer, dimension(:), intent(in), optional  :: ints
      integer, intent(in), optional                :: ll
      character(len=*), optional, intent(in)       :: fs
!
      character(len=1000)  :: pstring 
      character(len=20)    :: fstring 
!
      integer :: i, p_pos, int_check, i_err
      integer :: int_len, real_len, string_len
      integer :: int_count, real_count
      integer :: print_pos, printed
!
      integer           :: l_length
      character(len=20) :: f_string
!
!     Set print string and format string blank
      pstring = ' '
      fstring = ' '
!
      string_len = len_trim(string)
!
!     Check how many reals are present
      if(present(reals)) then
         real_len = size(reals)
      else
         real_len = 0
      endif
!
!     Check how many ints are present
      if(present(ints)) then
         int_len = size(ints)
      else
         int_len = 0
      endif
!
!     Line length to send to long_string_print
      if(present(ll)) then
         l_length = ll
      else
         l_length = 70
      endif
!
!     Format string to send to long_string_print
      if(present(fs)) then
         f_string = fs
      else
         f_string = '(t3,a)'
      endif
!
      real_count = 0
      int_count  = 0
!
      print_pos = 1
      printed = 1
      i = 0
!
      do while (i .lt. string_len)
!
         i = i + 1
!
!        Look for (
         if(string(i:i) .eq. "(") then
!
!           Is ( followed by f, F, e or E?
            if(string(i+1:i+1) .eq. "f" .or. string(i+1:i+1) .eq. "F" .or. &
               string(i+1:i+1) .eq. "e" .or. string(i+1:i+1) .eq. "E") then
!
!              Is it followed by a number, if so, assume a format string
               read(string(i+2:),'(i1)',iostat=i_err) int_check
               if (i_err .eq. 0) then
!
                  real_count = real_count + 1
                  if (real_count .gt. real_len) call the_file%error_msg('Not enough reals in printf')
!
!                 Print everything between previous print and (
                  write(pstring(print_pos:),'(a)') string(printed:i-1)
                  print_pos = print_pos + i - printed
                  p_pos = i
!
!                 Get length of format string
                  do while (string(i:i) .ne. ")" .and. i .le. string_len)
                     i = i + 1
                  enddo
!
                  printed = i+1
!
!                 Copy format string to fstring and write the next real to pstring
                  fstring = string(p_pos:i)
                  write(pstring(print_pos:),fstring) reals(real_count)
!
!                 Set next position to print
                  print_pos = len_trim(pstring) + 1
!
               endif
!
!           Is ( followed by i or I?
            elseif(string(i+1:i+1) .eq. "i" .or. string(i+1:i+1) .eq. "I") then
!
!              Is it followed by a number, if so, assume a format string
               read(string(i+2:),'(i1)',iostat=i_err) int_check
               if (i_err .eq. 0) then
!
                  int_count = int_count + 1
                  if (int_count .gt. int_len) call the_file%error_msg('Not enough ints in printf')
!
!                 Print everything between previous print and (
                  write(pstring(print_pos:),'(a)') string(printed:i-1)
                  print_pos = print_pos + i - printed
                  p_pos = i
!
!                 Get length of format string
                  do while (string(i:i) .ne. ")" .and. i .le. string_len)
                     i = i + 1
                  enddo
!
                  printed = i + 1
!
!                 Copy format string to fstring and write the next real to pstring
                  fstring = string(p_pos:i)
                  write(pstring(print_pos:),fstring) ints(int_count)
!
!                 Set next position to print
                  print_pos = len_trim(pstring) + 1
!
               endif
            endif
!  
         elseif (i .eq. string_len) then
!
!           Reached the end, print the rest
            write(pstring(print_pos:),'(a)') string(printed:i)
!
         endif
!
      enddo 
!
      call the_file%long_string_print(pstring, fs=f_string, ll=l_length)
!
   end subroutine printf_output_file
!
!  
   subroutine printd_output_file(the_file, string)
!!
!!    Print formatted
!!    Written by Rolf Heilemann Myhre, May 2019
!!
      implicit none
!
      class(output_file), intent(in) :: the_file
!
      character(len=*), intent(in) :: string
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      write(the_file%unit, '(t3,a)', iostat=io_error, iomsg=io_msg) string
!
      if (io_error /= 0) stop 'Error: could not print to eT output file '//trim(the_file%name_)//&
                             &'error message: '//trim(io_msg)
!
   end subroutine printd_output_file
!
!  
   subroutine author_output_file(the_file, author, contribution)
!!
!!    Print formatted
!!    Written by Rolf Heilemann Myhre, May 2019
!!
      implicit none
!
      class(output_file), intent(in) :: the_file
!
      character(len=*), intent(in) :: author
      character(len=*), intent(in) :: contribution
!
      character(len=100) :: d_author
      character(len=100) :: d_contribution
!
      d_author = author
      d_contribution = contribution
!
      write(the_file%unit,'(t4,a23,a54)') d_author, d_contribution
!
   end subroutine author_output_file
!
!  
   subroutine long_string_print_output_file(the_file, string, fs, colons, &
                                            ffs, lfs, ll)
!!
!!    Long string print
!!    Written by Rolf H. Myhre, Nov 2018
!!
!!    Prints a reasonbly formatted long string over several lines
!!    fs: optional format string
!!    ffs: optional format string for first printed line, 
!!                    will be used for single line if present
!!    lfs: optional format string for last printed line
!!    ll: optional argument for length printed lines, 
!!                 routine adds extra length to not split up words
!!
!!
      implicit none
!
      class(output_file), intent(in) :: the_file
!
      character(len=*), intent(in) :: string
      character(len=*), intent(in), optional :: fs, ffs, lfs
      integer, intent(in), optional :: ll
      logical, intent(in), optional :: colons
!
      character(90)  :: temp
      integer        :: l, l_left, lines, l_length
      integer        :: i,j, padd, printed 
      character(20)  :: f_s,fstring,ffstring,lfstring
      logical        :: col 
!
!     Default line length
      l_length = 70
      if(present(ll)) then
         l_length = ll
      endif
!
!     Figure out the formatting
      fstring = '(t3,a)'
      if(present(fs)) then
         fstring = fs
      endif
!
      ffstring = fstring
      if(present(ffs)) then
         ffstring = ffs
      endif
!
      lfstring = fstring
      if(present(lfs)) then
         lfstring = lfs
      endif
!
!     First line format
      f_s = ffstring
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
      lines = (l-1)/l_length + 1
      printed = 1
!
      do i = 1,lines
!
         if(i .ne. lines) then
!
            if(i .ne. 1) then
               f_s = fstring
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
               write(the_file%unit,f_s) ':: '//temp(1:l_length+padd+1)
            else
               write(the_file%unit,f_s) temp(1:l_length+padd+1)
            endif
!
         else
!           Print the remaining string
            f_s = lfstring
            if(col) then
               write(the_file%unit,f_s) ':: '//string(printed:l)
            else
               write(the_file%unit,f_s) string(printed:l)
            endif
            printed = l
         endif
!
      enddo
!
      call the_file%flush_()
!
   end subroutine long_string_print_output_file
!
!
end module output_file_class
