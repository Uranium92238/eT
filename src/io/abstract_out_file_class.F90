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
module abstract_out_file_class
!
!!
!!    abstract out file class module
!!    Written by Rolf H. Myhre, Sep 2019
!!
!!
!
   use kinds
   use abstract_file_class, only : abstract_file
!
   type, abstract, extends(abstract_file) :: abstract_out_file
!
   character(len=1), dimension(5) :: simple_formats = ['a', 'b', 'i', 'l', 'x']
   character(len=1), dimension(4) :: float_formats =  ['f', 'e', 'c', 'z']
!
   contains
!
      procedure :: open_                     => open_abstract_out_file
      procedure :: close_                    => close_abstract_out_file
      procedure :: flush_                    => flush_abstract_out_file
!
      procedure, public :: preformat_print         => preformat_print_abstract_out_file
!
      procedure, public :: format_print            => format_print_abstract_out_file
!
      procedure, public :: format_print_matrix     => format_print_matrix_abstract_out_file
!
      procedure, public :: format_print_separator  => format_print_separator_abstract_out_file
!
      procedure, public :: format_print_vector     => format_print_vector_abstract_out_file
!
      procedure, private :: long_string_print      => long_string_print_abstract_out_file
!
      procedure, private :: get_format_length      => get_format_length_abstract_out_file
      procedure, private :: get_tab_length         => get_tab_length_abstract_out_file
      procedure, private :: is_format              => is_format_abstract_out_file
      procedure, private :: get_n_numbers          => get_n_numbers_abstract_out_file
      procedure, nopass, private :: is_number      => is_number_abstract_out_file
      procedure, nopass, private :: get_format_end => get_format_end_abstract_out_file
      procedure, nopass, private :: get_n_escapes  => get_n_escapes_abstract_out_file
!
      procedure, private :: modify_real_format_abstract_out_file
      procedure, private :: modify_complex_format_abstract_out_file
      procedure, private :: modify_int_format_abstract_out_file
      procedure, private :: modify_log_format_abstract_out_file
      procedure, private :: modify_char_format_abstract_out_file
      procedure, private :: modify_blank_format_abstract_out_file
!
      generic :: modify_format => modify_real_format_abstract_out_file,    &
                                  modify_complex_format_abstract_out_file, &
                                  modify_int_format_abstract_out_file,     &
                                  modify_log_format_abstract_out_file,     &
                                  modify_char_format_abstract_out_file,    &
                                  modify_blank_format_abstract_out_file
!
   end type abstract_out_file
!
contains
!
!
   subroutine open_abstract_out_file(the_file, position_)
!!
!!    Open the abstract_out file
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(abstract_out_file) :: the_file
!
      integer              :: io_error
      character(len=100)   :: io_msg
      character(len=*), optional, intent(in) :: position_
      character(len=20)    :: pos
!
!
      if(present(position_)) then
         pos = trim(position_)
      else
         pos = 'rewind'
      endif
!
      if (the_file%is_open) then
!
         print *, trim(the_file%name_)//' is already open'
         stop
!
      endif
!
      open(newunit=the_file%unit_, file=the_file%name_, access=the_file%access_, &
           action='write', status='unknown', form=the_file%format_, position=pos, &
           iostat=io_error, iomsg=io_msg)
!
      if (io_error /= 0) then
!
         print *, 'Error: could not open eT abstract_out file '//trim(the_file%name_)//&
                             &'. Error message: '//trim(io_msg)
         stop
!
      endif
!
      the_file%is_open = .true.
!
   end subroutine open_abstract_out_file
!
!
   subroutine close_abstract_out_file(the_file)
!!
!!    Close the output file
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(abstract_out_file) :: the_file
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      if (.not. the_file%is_open) then
         print *, trim(the_file%name_)//' already closed'
         stop
      end if
!
      close(the_file%unit_, iostat=io_error, iomsg=io_msg, status='keep')
!
      if (io_error.ne. 0) then
!
         print *, 'Error: could not close eT output file '//trim(the_file%name_)//&
              &'error message: '//trim(io_msg)
         stop
!
      endif
!
      the_file%is_open = .false.
      the_file%unit_ = -1
!
   end subroutine close_abstract_out_file
!
!
   subroutine flush_abstract_out_file(the_file)
!!
!!    Flush the output file
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(abstract_out_file) :: the_file
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      flush(the_file%unit_, iostat=io_error, iomsg=io_msg)
!
      if (io_error /= 0) then
!
         print *, 'Error: could not flush eT output file '//trim(the_file%name_)//&
                             &'error message: '//trim(io_msg)
         stop
!
      endif
!
   end subroutine flush_abstract_out_file
!
!
   subroutine preformat_print_abstract_out_file(the_file, string, reals, complexs, ints, chars, logs, &
                                                fs, ffs, ll, padd, adv)
!!
!!    Preformat print
!!    Written by Rolf H. Myhre, May 2019
!!
!!    Modify format strings in string before passing them on.
!!    If a number is found in front of a format string, the string is multiplied by that number
!!    If a * is found, the number of format strings is set equal to the number of variables
!!    Example: 4(f9.3) -> (9.3)(9.3)(9.3)(9.3)
!!    Multipliers can be escaped using "\" for example, "\12(i0), ints=[1]" will be printed as 121
!!
!!    string:  String of characters that should be printed,
!!             including formatting of reals, integers, characters and logicals
!!
!!    reals:   Optional array of reals to print - in the order specified by string
!!             Default: None
!!    ints:    Optional array of integers to print - in the order specified by string
!!             Default: None
!!    chars:   Optional array of character strings to print - in the order specified
!!             by string. Note that all the strings must be of same length in Fortran
!!             Default: None
!!    logs:    Optional array of logicals to print - in the order specified by string
!!             Default: None
!!
!!    fs:      Optional character string specifies the format of the entire string,
!!             e.g. fs='(/t6,a)' gives a new line, then indentation 5, then the value
!!             of 'string' with reals and integers as specified. Default: '(t3,a)'
!!    ffs:     Optional character string specifies the format of the first printed line if
!!             different from fs. Default: same as fs
!!
!!    ll:      Optional integer specifying number of characters to print per line before
!!             looking for a white space to add a line break after. Default: 70
!!    padd:    Optional integer specifies how many characters beyond ll to search for
!!             a white space. Default: 18
!!
!!    adv:     Optional logical specifies whether advance is 'yes' or 'no' for the last line.
!!             Default: .true.
!!
!!    Note that the number of characters to print per line will typically be between ll and
!!    ll + padd minus the number of blank spaces specified by the t specifier in the format
!!    string, assuming there are enough characters to print.
!!
      implicit none
!
      class(abstract_out_file), intent(in) :: the_file
!
!     Data to print
      character(len=*), intent(in)                          :: string
      real(dp), dimension(:), intent(in), optional          :: reals
      complex(dp), dimension(:), intent(in), optional       :: complexs
      integer, dimension(:), intent(in), optional           :: ints
      character(len=*), dimension(:), intent(in), optional  :: chars
      logical, dimension(:), intent(in), optional           :: logs
!
!     Parameters to pass to long_string_print
      integer, optional, intent(in)                         :: ll
      integer, optional, intent(in)                         :: padd
      character(len=*), optional, intent(in)                :: fs
      character(len=*), optional, intent(in)                :: ffs
!
      logical, intent(in), optional :: adv
!
      character(len=1000) :: pstring
      character(len=40)   :: fstring
!
      integer :: format_start, format_end, n_escapes, n_numbers
      integer :: string_len, factor, f_length
      integer :: i, j
      integer :: print_position, printed, print_end
!
!     Set print string and format string blank
      pstring = ''
      fstring = ''
!
      string_len = len_trim(string)
!
      print_position = 1
      printed = 0
      i = 1
      n_escapes = 0
!
      do while (i .lt. string_len)
!
         i = i + 1
!
!        Look for a format string with a number in front
         if(the_file%is_format(string(i:))) then
            if(the_file%is_number(string(i-1:i-1)) .or. string(i-1:i-1) .eq. "*") then
!
               format_start = i
!
               f_length = the_file%get_format_end(string(format_start:))
               format_end = format_start + f_length - 1
               fstring = string(format_start:format_end)
!
               if(string(i-1:i-1) .eq. "*") then
                  n_numbers = 1
                  factor = 0
                  if(fstring(2:2) .eq. "f" .or. fstring(2:2) .eq. "e") then
                     if(present(reals)) factor = size(reals)
                  elseif(fstring(2:2) .eq. "c" .or. fstring(2:2) .eq. "z") then
                     if(present(complexs)) factor = size(complexs)
                  elseif(fstring(2:2) .eq. "i") then
                     if(present(ints)) factor = size(ints)
                  elseif(fstring(2:2) .eq. "l") then
                     if(present(logs)) factor = size(logs)
                  else
                     if(present(chars)) factor = size(chars)
                  endif
               else
                  n_numbers = the_file%get_n_numbers(string(:format_start - 1))
                  read(string(i-n_numbers:i-1),*) factor
               endif
!
!              Print everything between previous print and number
               if(format_start .gt. 2) then
!
!                 Figure out the number of escape characters "\"
                  n_escapes = the_file%get_n_escapes(string(:format_start - n_numbers - 1))
                  print_end = format_start - n_numbers - n_escapes - 1
!
                  write(pstring(print_position:),'(a)') string(printed+1 : print_end)
                  print_position = print_position + print_end - printed
!
               endif
!
!              Process format multiplication if even number of escape characters
               if(mod(n_escapes,2) .eq. 0) then
!
                  do j = 1, n_escapes
                     write(pstring(print_position:),'(a)') "\"
                     print_position = print_position + 1
                  enddo
!
                  do j = 1,factor
                     write(pstring(print_position:),'(a)') fstring
                     print_position = print_position + f_length
                  enddo
!
                  printed = format_start - 1
!
               else !Escape
!
                  do j = 1, n_escapes/2
                     write(pstring(print_position:),'(a)') "\"
                     print_position = print_position + 1
                  enddo
!
                  write(pstring(print_position:),'(a)') &
                        string(format_start-n_numbers:format_start-1)
                  print_position = print_position + n_numbers
                  write(pstring(print_position:),'(a)') fstring
                  print_position = print_position + f_length
!
                  printed = format_start - n_numbers - 1
               endif
               i = format_end
            endif
         elseif (i .eq. string_len) then
!
!           Reached the end, print the rest
            write(pstring(print_position:),'(a)') string(printed+1:i)
!
         endif
!
      enddo
!
      call the_file%format_print(pstring, reals, complexs, ints, chars, logs, &
                                 fs, ffs, ll, padd, adv)
!
   end subroutine preformat_print_abstract_out_file
!
!
   subroutine format_print_abstract_out_file(the_file, string, reals, complexs, ints, chars, logs, &
                                             fs, ffs, ll, padd, adv)
!!
!!    Format print
!!    Written by Rolf H. Myhre, May 2019
!!
!!    Prints any number of reals, integers, characters and logicals formatted Python style.
!!
!!    string:  String of characters that should be printed,
!!             including formatting of reals, integers, characters and logicals
!!
!!    reals:   Optional array of reals to print - in the order specified by string
!!             Default: None
!!    ints:    Optional array of integers to print - in the order specified by string
!!             Default: None
!!    chars:   Optional array of character strings to print - in the order specified
!!             by string. Note that all the strings must be of same length in Fortran
!!             Default: None
!!    logs:    Optional array of logicals to print - in the order specified by string
!!             Default: None
!!
!!    fs:      Optional character string specifies the format of the entire string,
!!             e.g. fs='(/t6,a)' gives a new line, then indentation 5, then the value
!!             of 'string' with reals and integers as specified. Default: '(t3,a)'
!!    ffs:     Optional character string specifies the format of the first printed line if
!!             different from fs. Default: same as fs
!!
!!    ll:      Optional integer specifying number of characters to print per line before
!!             looking for a white space to add a line break after. Default: 70
!!    padd:    Optional integer specifies how many characters beyond ll to search for
!!             a white space. Default: 18
!!
!!    adv:     Optional logical specifies whether advance is 'yes' or 'no' for the last line.
!!             Default: .true.
!!
!!    Note that the number of characters to print per line will typically be between ll and
!!    ll + padd minus the number of blank spaces specified by the t specifier in the format
!!    string, assuming there are enough characters to print.
!!
      implicit none
!
      class(abstract_out_file), intent(in) :: the_file
!
!     Data to print
      character(len=*), intent(in)                          :: string
      real(dp), dimension(:), intent(in), optional          :: reals
      complex(dp), dimension(:), intent(in), optional       :: complexs
      integer, dimension(:), intent(in), optional           :: ints
      character(len=*), dimension(:), intent(in), optional  :: chars
      logical, dimension(:), intent(in), optional           :: logs
!
!     Parameters to pass to long_string_print
      integer, optional, intent(in)                         :: ll
      integer, optional, intent(in)                         :: padd
      character(len=*), optional, intent(in)                :: fs
      character(len=*), optional, intent(in)                :: ffs
!
      logical, intent(in), optional :: adv
!
      character(len=1000) :: pstring
      character(len=40)   :: fstring
      character(len=1)    :: ftype
!
      integer :: format_start, format_end, format_length, n_escapes
      integer :: string_len
      integer :: i, j
      integer :: nints, nreals, nlogs, nchars, ncomplexs
      integer :: int_count, real_count, log_count, char_count, complex_count
      integer :: print_position, printed, print_end
      integer :: length
!
!     Set print string and format string blank
      pstring = ' '
      fstring = ' '
!
      string_len = len_trim(string)
      format_length = 0
!
      real_count    = 0
      int_count     = 0
      log_count     = 0
      char_count    = 0
      complex_count = 0
!
      nreals    = 0
      nints     = 0
      nlogs     = 0
      nchars    = 0
      ncomplexs = 0
!
!     Set number of reals, ints, logicals or chars if present
      if(present(reals)) nreals       = size(reals)
      if(present(complexs)) ncomplexs = size(complexs)
      if(present(ints))  nints        = size(ints)
      if(present(logs))  nlogs        = size(logs)
      if(present(chars)) nchars       = size(chars)
!
      print_position = 1
      printed = 0
      i = 0
      n_escapes = 0
!
      do while (i .lt. string_len)
!
         i = i + 1
!
!        Look for a format string
         if(the_file%is_format(string(i:))) then
!
            format_start = i
            format_end = format_start + the_file%get_format_end(string(i:)) - 1
            fstring = string(format_start:format_end)
            ftype = fstring(2:2)
!
!           Print everything between previous print and "\" or "("
            if(format_start .gt. 1) then
!
!              Figure out the number of escape characters "\"
               n_escapes = the_file%get_n_escapes(string(:format_start - 1))
               print_end = format_start - n_escapes - 1
!
               write(pstring(print_position:),'(a)') string(printed+1 : print_end)
               print_position = print_position + print_end - printed
!
            endif
!
            do j = 1, n_escapes/2
               write(pstring(print_position:),'(a)') "\"
               print_position = print_position + 1
            enddo
!
!           Process format if even number of escape characters
            if(mod(n_escapes,2) .eq. 0) then
!
!              Is ( followed by f or e?
               if(ftype .eq. "f" .or. ftype .eq. "e") then
!
                  real_count = real_count + 1
                  if (real_count .gt. nreals) then
                     print *, 'Not enough reals in printf'
                     stop
                  endif
!
!                 Modify fstring if f0 and get format length
                  call the_file%modify_format(fstring, reals(real_count), format_length)
!
!                 Copy format string to fstring and write the next real to pstring
                  write(pstring(print_position:),fstring) reals(real_count)
!
!              Is ( followed by c or z?
               elseif(ftype .eq. "c" .or. ftype .eq. "z") then
!
                  complex_count = complex_count + 1
                  if (complex_count .gt. ncomplexs) then
                     print *, 'Not enough complexs in printf'
                     stop
                  endif
!
!                 Modify fstring if f0 and get format length
                  call the_file%modify_format(fstring, complexs(complex_count), format_length)
!
!                 Copy format string to fstring and write the next real to pstring
                  write(pstring(print_position:),fstring) real(complexs(complex_count)), &
                                                          abs(aimag(complexs(complex_count)))
!
!              Is ( followed by i?
               elseif(ftype .eq. "i") then
!
                  int_count = int_count + 1
                  if (int_count .gt. nints) then
                     print *, 'Not enough ints in printf'
                     stop
                  endif
!
!                 Modify fstring if i0
                  call the_file%modify_format(fstring, ints(int_count), format_length)
!
                  write(pstring(print_position:),fstring) ints(int_count)
!
                  format_length = the_file%get_format_length(fstring)
!
!              Is ( followed by l?
               elseif(ftype .eq. "l") then
!
                  log_count = log_count + 1
                  if (log_count .gt. nlogs) then
                     print *, 'Not enough logicals in printf'
                     stop
                  endif
!
!                 Modify fstring to a and figure out length if l0
                  call the_file%modify_format(fstring, logs(log_count), format_length)
!
                  if(logs(log_count)) then
                     write(pstring(print_position:), fstring) 'True'
                  else
                     write(pstring(print_position:), fstring) 'False'
                  endif
!
!              Is ( followed by a or b?
               elseif((ftype .eq. "a") .or. (ftype .eq. "b")) then
!
                  char_count = char_count + 1
                  if (char_count .gt. nchars) then
                     print *,  'Not enough chars in printf'
                     stop
                  endif
!
                  call the_file%modify_format(fstring, chars(char_count), format_length)
!
                  if (format_length .gt. 0) then
                     write(pstring(print_position:), fstring) trim(chars(char_count))
                  endif
!
!              Is ( followed by x?
               elseif(ftype .eq. "x") then
!
                  call the_file%modify_format(fstring, format_length)
!
                  write(pstring(print_position:), fstring)
                  format_length = the_file%get_format_length(fstring)
!
               endif
!
               i = format_end
               printed = format_end
               print_position = print_position + format_length
!
            else !Escape
!
               printed = format_start - 1
!
            endif
!
         elseif (i .eq. string_len) then
!
!           Reached the end, print the rest
            write(pstring(print_position:),'(a)') string(printed+1:i)
!
            print_position = print_position + i - printed
!
         endif
!
      enddo
!
      length = print_position - 1
!
      if(nreals .ne. real_count) then
         print *, "Too many reals in format_print"
         stop
      endif
      if(nints .ne. int_count) then
         print *, "Too many ints in format_print"
         stop
      endif
      if(ncomplexs .ne. complex_count) then
         print *, "Too many complexs in format_print"
      endif
      if(nlogs .ne. log_count) then
         print *, "Too many logs in format_print"
      endif
      if(nchars .ne. char_count) then
         print *, "Too many logs in format_print"
      endif
!
      call the_file%long_string_print(pstring, length, fs, ffs, ll, padd, adv)
!
   end subroutine format_print_abstract_out_file
!
!
   subroutine long_string_print_abstract_out_file(the_file, string, length, fs, ffs, ll, padd, adv)
!!
!!    Long string print
!!    Written by Rolf H. Myhre, Nov 2018
!!
!!    Prints a character string, possibly dividing it up over several lines
!!    based on ll, padd and format and prints them to the_file
!!
!!    string:  Character string to print
!!
!!    length:  Length of string
!!
!!    fs:      Optional character string specifies the format of the entire string,
!!             e.g. fs='(/t6,a)' gives a new line, then indentation 5, then the value
!!             of 'string' with reals and integers as specified. Default: '(t3,a)'
!!    ffs:     Optional character string specifies the format of the first printed line if
!!             different from fs. Default: same as fs
!!
!!    ll:      Optional integer specifying number of characters to print per line before
!!             looking for a white space to add a line break after. Default: 70
!!    padd:    Optional integer specifies how many characters beyond ll to search for
!!             a white space. Default: 18
!!
!!    adv:     Optional logical specifies whether advance is 'yes' or 'no' for the last line.
!!             Default: .true.
!!
!!    Note: The number of characters to print per line will typically be between ll and
!!    ll + padd minus the number of blank spaces specified by the t specifier in the format
!!    string, assuming there are enough characters to print.
!!
      implicit none
!
      class(abstract_out_file), intent(in) :: the_file
!
      character(len=*), intent(in)           :: string
      integer, intent(in)                    :: length
      character(len=*), intent(in), optional :: fs, ffs
      integer, intent(in), optional          :: ll
      logical, intent(in), optional          :: adv
      integer, intent(in), optional          :: padd
!
      character(200) :: temp
      integer        :: length_with_tab
      integer        :: temp_length, printed
      integer        :: line_length, fline_length, l_l
      integer        :: padding, max_padd
      character(20)  :: f_s, fstring, ffstring
      logical        :: do_advance
!
      character(len=3) :: advancing
!
!     Just print a new line and return if empty
      if (length .eq. 0) then
         write(the_file%unit_, *)
         return
      endif
!
!     Default line length
      length_with_tab = 70
      if(present(ll)) then
         length_with_tab = ll
      endif
!
      if(length_with_tab .le. 0) then
         print *, "Too short line length: ", length_with_tab
         stop
      endif
!
!     Figure out the formatting for the bulk
      fstring = '(t3,a)'
      if(present(fs)) then
         fstring = fs
      endif
      line_length = length_with_tab - the_file%get_tab_length(fstring)
!
      if(line_length .le. 0) then
         print *, "Too many tabs or too short line length: ", line_length
         stop
      endif
!
!     Figure out the formatting of the first line
      ffstring = fstring
      if(present(ffs)) then
         ffstring = ffs
      endif
      fline_length = length_with_tab - the_file%get_tab_length(ffstring)
!
      if(fline_length .le. 0) then
         print *, "Too many tabs or too short first line length: ", fline_length
         stop
      endif
!
!     How much padding?
      max_padd = 18
      if(present(padd)) then
         max_padd = padd
      endif
!
      if(max_padd .lt. 0) then
         print *, "Padding is less than zero: ", max_padd
         stop
      endif
!
!     advancing?
      do_advance = .true.
      if(present(adv)) then
         do_advance = adv
      endif
!
!     First line format and length
      f_s = ffstring
      l_l = fline_length
      advancing = 'yes'
!
      printed = 0
!
      do while (printed .lt. length)
!
         temp_length = l_l
!
!        See if we've reached the end
         if ((printed + temp_length) .ge. length) then
            temp_length = length - printed
!
!        Search for a white space or the end of the string within max_padd
!        Note that the white space will be included in the print
         else
            do padding = 0, max_padd
!
               if((string(printed + temp_length : printed + temp_length) .eq. ' ') .or. &
                  (printed + temp_length .eq. length)) then
                  exit
               endif
!
               temp_length = temp_length + 1
!
            enddo
         endif
!
!        Check if we reached the end and set advancing parameter and disable hyphen
         if ((printed + temp_length) .eq. length) then
            padding = -1
            if (.not. do_advance) then
               advancing = 'no'
            end if
         end if
!
!        Copy string to be printed
         temp(1 : temp_length) = string(printed + 1 : printed + temp_length)
         printed = printed + temp_length
!
!        Add hyphen if splitting words and max_padd is not zero
         if((padding .eq. max_padd) .and. (max_padd .ne. 0)) then
            temp_length = temp_length + 1
            temp(temp_length : temp_length) = '-'
         endif
!
!        Write to file
         write(the_file%unit_, f_s, advance=advancing) temp(1 : temp_length)
!
         f_s = fstring
         l_l = line_length
!
      enddo
!
   end subroutine long_string_print_abstract_out_file
!
!
   pure function get_format_length_abstract_out_file(the_file, fstring) result(length)
!
!!    Get printed length of format
!!    Written by Rolf H. Myhre, Sep 2019
!!
!!    Figure out how many printed characters a format string corresponds to
!!
!!    Note! The function assumes the string is a format string and does not work for 0 formats
!!
      implicit none
!
      class(abstract_out_file), intent(in) :: the_file
!
      character(len=*), intent(in)  :: fstring
!
      character(len=1)  :: endchar
      integer :: i, length
!
      endchar = ")"
!
!     Look for "." if floating point format
      if(any(fstring(2:2) .eq. the_file%float_formats)) then
         endchar = "."
      endif
!
      i = 3
      do while (fstring(i+1:i+1) .ne. endchar)
         i = i + 1
      enddo
!
      read(fstring(3:i),*) length
!
   end function get_format_length_abstract_out_file
!
!
   pure function get_n_escapes_abstract_out_file(string) result(n_escapes)
!
!!    Get n escapes
!!    Written by Rolf H. Myhre, Nov 2020
!!
!!    Returns the number of consecutive escape characters at the end of string.
!!    "f\sf\\\" would return 3
!!
      implicit none
!
      character(len=*), intent(in)  :: string
!
      integer :: n_escapes, i
!
      n_escapes = 0
      if (len(string) .eq. 0) return
!
      do i = len(string), 1, -1
         if (string(i:i) .eq. "\") then
            n_escapes = n_escapes + 1
         else
            return
         endif
      enddo
!
   end function get_n_escapes_abstract_out_file
!
!
   pure function get_n_numbers_abstract_out_file(the_file, string) result(n_numbers)
!
!!    Get n numbers
!!    Written by Rolf H. Myhre, Nov 2020
!!
!!    Returns number at the end of string if it exists,
!!    for example, "f2sf456" would return 3
!!
      implicit none
!
      class(abstract_out_file), intent(in) :: the_file
!
      character(len=*), intent(in)  :: string
!
      integer :: n_numbers, i
!
      n_numbers = 0
      if (len(string) .le. 1) return
!
      do i = len(string), 1, -1
         if (the_file%is_number(string(i:i))) then
            n_numbers = n_numbers + 1
         else
            return
         endif
      enddo
!
   end function get_n_numbers_abstract_out_file
!
!
   pure function get_format_end_abstract_out_file(string) result(i)
!
!!    Get format end
!!    Written by Rolf H. Myhre, Nov 2020
!!
!!    Returns the position of first ")"
!!
!!    Note! Routine assumes there is a ")" in the string
!!
      implicit none
!
      character(len=*), intent(in)  :: string
!
      integer :: i
!
      i = 1
      do while (string(i:i) .ne. ")")
         i = i + 1
      enddo
!
   end function get_format_end_abstract_out_file
!
!
   pure function is_format_abstract_out_file(the_file, string) result(is_format)
!
!!    is format
!!    Written by Rolf H. Myhre, Oct 2020
!!
!!    Check if string starts with a format string,
!!
!!
      implicit none
!
      class(abstract_out_file), intent(in) :: the_file
!
      character(len=*), intent(in)   :: string
!
      logical   :: is_format
      integer   :: i, length, maxlength
!
      is_format = .false.
      maxlength = 7
!
      length = len_trim(string)
      if (length .lt. 4) return
!
!     Check that the first character is a "("
!     and that the second is one of the formatting characters
!     and that the third is a number
      if(string(1:1) .eq. "(") then
         if (any(string(2:2) .eq. the_file%float_formats) .or. &
             any(string(2:2) .eq. the_file%simple_formats)) then
            if (the_file%is_number(string(3:3))) then
!
               i = 4
!
!              Look for "." in float formats
               if(any(string(2:2) .eq. the_file%float_formats)) then
!
!                 Check that there are only numbers until "."
                  do while (string(i:i) .ne. "." .and. i .lt. length-1 .and. i .lt. maxlength)
!
                     if (.not. the_file%is_number(string(i:i))) exit
                     i = i + 1
!
                  enddo
!
!                 Check that we have a "."
                  if ((i .gt. length - 2) .or. &
                      (.not. string(i:i) .eq. ".")) then
                     return
                  else
                     maxlength = maxlength + i
                     i = i + 1
                  endif
               endif
!
!              Check that there are only numbers until )
               do while (string(i:i) .ne. ")" .and. i .lt. length .and. i .lt. maxlength)
!
                  if (.not. the_file%is_number(string(i:i))) exit
                  i = i + 1
!
               enddo
!
!              Check that we have a ")"
               if ((i .gt. length) .or. &
                   (.not. string(i:i) .eq. ")")) then
                  return
               else
                  is_format = .true.
               endif
!
            endif
         endif
      endif
!
   end function is_format_abstract_out_file
!
!
   function get_tab_length_abstract_out_file(the_file, fstring) result(length)
!
!!    Get tab length
!!    Written by Rolf H. Myhre, Sep 2019
!!
!!    Figure out how many blank spaces are added in front by a format string
!!    For example fstring='(/t6,a)' should return 5
!
      implicit none
!
      class(abstract_out_file), intent(in) :: the_file
!
      character(len=*), intent(in)  :: fstring
!
      integer :: i, j
      integer :: length, string_length, stat
!
      string_length = len_trim(fstring)
!
      length = 0
!
!     If character, integer or logical
      if (fstring(1:1) .eq. '(') then
!
         if (fstring(string_length:string_length) .ne. ')') then
!
            print *, fstring//' does not end with )'
            stop
!
         endif
!
!        Loop over length of string
         do i = 2,string_length-2
!
!           If it doesn't start with a tab we don't care
            if(fstring(i:i) .ne. '/' .and. fstring(i:i) .ne. 't') then
               exit
!
!           If we find a t
            elseif(fstring(i:i) .eq. 't') then
!
!              Check that we're not at the end, then something is very wrong
               if (i .lt. string_length) then
!
!                 Count the numbers after the t
                  j = i
                  do while (the_file%is_number(fstring(j+1:j+1)))
                     j = j + 1
                  enddo
!
!                 If we found any numbers, read them into the length int
                  if (j .gt. i) then
!
                     read(fstring(i+1:j),*,iostat=stat) length
                     if (stat .ne. 0) then
                        print *, 'Failed to read numbers from '//fstring
                        stop
                     endif
!
!                 No numbers
                  else
                     print *, 'No numbers after t in '//fstring
                     stop
                  endif
!
               else
                  print *, fstring//' is not a format string'
                  stop
               endif
!
               length = length - 1
!
               if (length .lt. 0) then
                  print *, 'Tab length less than 0 for '//fstring
                  stop
               endif
!
               exit
!
            endif
         enddo
      else
         print *, fstring//' is not a format string'
         stop
      endif
!
   end function get_tab_length_abstract_out_file
!
!
   pure function is_number_abstract_out_file(check_char) result(it_is)
!
!!    is number
!!    Written by Rolf H. Myhre, Nov 2019
!!
!!    Check if check_char is a number or not
!!
!
      implicit none
!
      character, intent(in)   :: check_char
      logical                 :: it_is
!
      character(len=1), dimension(10)  :: zero_to_nine
!
      zero_to_nine=['0','1','2','3','4','5','6','7','8','9']
!
      it_is = any(check_char .eq. zero_to_nine)
!
   end function is_number_abstract_out_file
!
!
   pure subroutine modify_real_format_abstract_out_file(the_file, fstring, variable, format_length)
!
!!    modify real format
!!    Written by Rolf H. Myhre, Nov 2019
!!
!!    Set correct length for f0 and e0 and set format_length
!!
      implicit none
!
      class(abstract_out_file), intent(in) :: the_file
      character(len=*), intent(inout)      :: fstring
      real(dp), intent(in)                 :: variable
      integer, intent(out)                 :: format_length
!
      character(len=40) :: tempstring
      integer :: j, f_length, stringlength
!
      if (fstring(3:4) .eq. "0.") then
!
         stringlength = len_trim(fstring)
!
         read(fstring(5:stringlength-1),*) f_length
!
         if(fstring(2:2) .eq. "f") then
            j = 1
            do while (abs(variable)/(10**j) .gt. 1.0)
               j = j + 1
            enddo
         else
            j = 5
         endif
!
         f_length = f_length + j + 1
         if (variable .lt. 0.0) f_length = f_length + 1
!
         tempstring = fstring(4:stringlength)
         write(fstring(3:), '(i0,a)') f_length, tempstring(1:stringlength-3)
!
      endif
!
      format_length = the_file%get_format_length(fstring)
!
   end subroutine modify_real_format_abstract_out_file
!
!
   pure subroutine modify_complex_format_abstract_out_file(the_file, fstring, variable, format_length)
!
!!    modify complex format
!!    Written by Rolf H. Myhre, Nov 2019
!!
!!    Set complex format to two real formats and set format_length
!!
      implicit none
!
      class(abstract_out_file), intent(in) :: the_file
      character(len=*), intent(inout)      :: fstring
      complex(dp), intent(in)              :: variable
      integer, intent(out)                 :: format_length
!
      character(len=40) :: tempstring
      character(len=5)  :: plusminus
      character(len=4)  :: istring
      integer :: j, k, real_length
      logical :: zero_format
!
      istring = '"i")'
!
      plusminus = '" + "'
      if(aimag(variable) .lt. 0.0) plusminus(3:3) = '-'
      zero_format = .false.
      if(fstring(3:4) .eq. "0.") zero_format = .true.
!
      if(fstring(2:2) .eq. "c") then
         fstring(2:2) = "f"
      else
         fstring(2:2) = "e"
      endif
!
      if(zero_format) tempstring = trim(fstring)
!
      call the_file%modify_format(fstring, real(variable), real_length)
      j = len_trim(fstring)
      format_length = real_length + 3
!
      if (zero_format) then
!
         call the_file%modify_format(tempstring, abs(aimag(variable)), real_length)
         k = len_trim(tempstring)
!
         write(fstring(j:), '(a,a,a)') plusminus, tempstring(2:k-1), istring
!
      else
!
         write(fstring(j:), '(a,a,a)') plusminus, fstring(2:j-1), istring
!
      endif
!
      format_length = format_length + real_length + 1
!
   end subroutine modify_complex_format_abstract_out_file
!
!
   pure subroutine modify_int_format_abstract_out_file(the_file, fstring, variable, format_length)
!
!!    modify float format
!!    Written by Rolf H. Myhre, Nov 2019
!!
!!    Set correct length for i0 and set format_length
!!
      implicit none
!
      class(abstract_out_file), intent(in) :: the_file
      character(len=*), intent(inout)      :: fstring
      integer, intent(in)                  :: variable
      integer, intent(out)                 :: format_length
!
      integer :: j
!
      if (fstring(1:4) .eq. "(i0)") then
!
         j = 1
         do while (abs(variable)/(10**j) .gt. 0)
            j = j + 1
         enddo
!
         if (variable .lt. 0) j = j + 1
!
         write(fstring, '(a2,i0,a1)') "(i", j, ")"
!
      endif
!
      format_length = the_file%get_format_length(fstring)
!
   end subroutine modify_int_format_abstract_out_file
!
!
   pure subroutine modify_log_format_abstract_out_file(the_file, fstring, variable, format_length)
!
!!    modify log format
!!    Written by Rolf H. Myhre, Nov 2019
!!
!!    Change l formats to a and set format_length
!!
      implicit none
!
      class(abstract_out_file), intent(in) :: the_file
      character(len=*), intent(inout)      :: fstring
      logical, intent(in)                  :: variable
      integer, intent(out)                 :: format_length
!
      if (fstring(1:4) .eq. "(l0)") then
         if(variable) then
            fstring = "(a4)"
         else
            fstring = "(a5)"
         endif
      else
         fstring(2:2) = "a"
      endif
!
      format_length = the_file%get_format_length(fstring)
!
   end subroutine modify_log_format_abstract_out_file
!
!
   pure subroutine modify_char_format_abstract_out_file(the_file, fstring, variable, format_length)
!
!!    modify char format
!!    Written by Rolf H. Myhre, Nov 2019
!!
!!    Change b formats to a and set format_length
!!
      implicit none
!
      class(abstract_out_file), intent(in) :: the_file
      character(len=*), intent(inout)      :: fstring
      character(len=*), intent(in)         :: variable
      integer, intent(out)                 :: format_length
!
      integer :: len_trimmed
!
      format_length = the_file%get_format_length(fstring)
      len_trimmed = len_trim(variable)
!
      if (fstring(3:3) .eq. "0" .or. len_trimmed .eq. format_length) then
!
         write(fstring, '(a2,i0,a1)') "(a", len_trimmed, ")"
         format_length = len_trimmed
!
      else
         if (fstring(2:2) .eq. "b") then
            write(fstring, '(a2,i0,a1,i0,a2)') "(a", len_trimmed, &
                                               ",", format_length - len_trimmed, "x)"
         endif
      endif
!
   end subroutine modify_char_format_abstract_out_file
!
!
   pure subroutine modify_blank_format_abstract_out_file(the_file, fstring, format_length)
!
!!    modify blank format
!!    Written by Rolf H. Myhre, Nov 2019
!!
!!    Change format to fortran format and set format_length
!!
      implicit none
!
      class(abstract_out_file), intent(in) :: the_file
      character(len=*), intent(inout)      :: fstring
      integer, intent(out)                 :: format_length
!
      format_length = the_file%get_format_length(fstring)
      if(format_length .eq. 0) then
         write(fstring, '(a3)') "(a)"
      else
         write(fstring, '(a1,i0,a2)') '(', format_length, 'x)'
      endif
!
   end subroutine modify_blank_format_abstract_out_file
!
!
   subroutine format_print_matrix_abstract_out_file(the_file, name_, matrix, &
                                                    dim_1, dim_2, fs, columns)
!!
!!    Format print matrix
!!    Written by Tommaso Giovannini, Mar. 2019
!!
!!    Modified by Rolf H. Myhre, Oct. 2019
!!
!!    Moved to output file and added format string and number of columns
!!
!!    name_: Name to be printed above the matrix
!!
!!    Matrix to be printed with dimension dim_1 x dim_2
!!
!!    fs:      Optional format string for numbers, default is (f13.8)
!!    columns: Optional integer specifying number of columns to print per line, default is 5
!!
      implicit none
!
      class(abstract_out_file), intent(in)            :: the_file
!
      character(len=*), intent(in)                    :: name_
!
      integer, intent(in)                             :: dim_1
      integer, intent(in)                             :: dim_2
!
      real(dp), dimension(dim_1, dim_2), intent(in)   :: matrix
!
      character(len=*), optional, intent(in)          :: fs
      integer, intent(in), optional                   :: columns
!
      integer :: i, j, k
      integer :: n_columns, n_prints, columns_printed
      integer :: form_length, line_length, name_length, print_length
      integer :: row_int_length, col_int_length, first_col_int_length
      character(len=20)  :: form_string, print_fs
      character(len=20)  :: row_index_format, col_index_format, first_col_index_format
      character(len=400) :: int_string, real_string
!
      integer, dimension(:), allocatable  :: ints_to_print
      real(dp), dimension(:), allocatable :: reals_to_print
!
!     Set format if provided
      if (present(fs)) then
         form_string = trim(fs)
      else
         form_string = '(f13.8)'
      endif
!
!
!     Set number of columns if required
      if (present(columns)) then
         n_columns = min(dim_2, columns)
      else
         n_columns = min(dim_2, 5)
      endif
!
!     Tiny arrays and mem depends on output
      allocate(ints_to_print(n_columns))
      allocate(reals_to_print(n_columns))
!
!     Number of prints to devide the columns into
      n_prints = (dim_2-1)/n_columns + 1
!
!     Get various lengths
      name_length = len_trim(name_)
      form_length = the_file%get_format_length(form_string)
!
      write(int_string, '(i0)') dim_1
      row_int_length = len_trim(int_string)
!
      write(int_string, '(i0)') dim_2
      col_int_length = len_trim(int_string)
!
!     Calculate the line lengths
      line_length = (row_int_length + 1) + n_columns*(form_length+1)
!
!     Calculate the first line length, name should be in the middle
      print_length = (line_length + name_length)/2
!
!     Find format of first line and print name
      write(print_fs, '(a,i0,a)') '(/a', print_length, ')'
      call the_file%format_print(name_, fs=print_fs)
!
!     Set up string for column indices, assuming form_length .ge. col_int_length
      first_col_int_length = (row_int_length + 1) + form_length/2 + col_int_length/2 + 1
      col_int_length = form_length
!
      write(first_col_index_format,'(a,i0,a)') '(i', first_col_int_length, ')'
      write(col_index_format,'(a,i0,a)') '(i', col_int_length, ')'
      int_string = trim(first_col_index_format) // repeat(' '//trim(col_index_format), n_columns - 1)
!
!     Set up string for matrix elements
      write(row_index_format,'(a,i0,a)') '(i', row_int_length+1, ')'
      real_string = trim(row_index_format) // repeat(' '//trim(form_string), n_columns)
!
!     Print separator
      call the_file%format_print_separator(line_length+1, '=')
!
!     Start iterating over number of prints
!
      columns_printed = 0
!
      do i = 1,n_prints
!
!        Some space between the prints
         if (i .ne. 1) then
            call the_file%format_print('')
         end if
!
!        Calculate number of columns if last print and not first print and set up format strings
         if ((i .eq.  n_prints) .and. (i .ne. 1) .and. (mod(dim_2, n_columns) .ne. 0)) then
            n_columns = mod(dim_2, n_columns)
            int_string = trim(first_col_index_format) // &
                         repeat(' '//trim(col_index_format), n_columns - 1)
            real_string = trim(row_index_format) // repeat(' '//trim(form_string), n_columns)
         endif
!
!        Set up column indices to print
         do j = 1, n_columns
            ints_to_print(j) = columns_printed + j
         enddo
!
         call the_file%format_print(int_string, ints=ints_to_print(1:n_columns), ll=line_length)
!
         do k = 1, dim_1
!
            do j = 1, n_columns
               reals_to_print(j) = matrix(k,j + columns_printed)
            enddo
!
            call the_file%format_print(real_string, ints=[k], &
                                    reals=reals_to_print(1:n_columns), ll=line_length)
!
         enddo
!
!
         columns_printed = columns_printed + n_columns
!
      enddo
!
!     Print separator
      call the_file%format_print_separator(line_length+1, '=')
!
!
      deallocate(ints_to_print)
      deallocate(reals_to_print)
!
   end subroutine format_print_matrix_abstract_out_file
!
!
   subroutine format_print_separator_abstract_out_file(the_file, n, symbol, fs)
!!
!!    Format print separator
!!
!!    Written by Rolf H. Myhre, Oct. 2019
!!
!!    n: Number of symbols to print
!!    symbol: optional, what symbol to print. Default: '-'
!!    fs: optional format string. Default: '(t3,a)'
!!
      implicit none
!
      class(abstract_out_file), intent(in) :: the_file
!
      integer, intent(in) :: n
!
      character, intent(in), optional :: symbol
!
      character(len=*), optional, intent(in) :: fs
!
      character :: sym
!
      character(len=n) :: separator_line
!
      integer :: line_length
!
      if (present(symbol)) then
         sym = symbol
      else
         sym = '-'
      endif
!
      if (present(fs)) then
         line_length = n+the_file%get_tab_length(fs)
      else
         line_length = n + 2
      endif
!
      separator_line = repeat(sym, n)
!
      call the_file%format_print(separator_line, ll=line_length, fs=fs, padd=0)
!
   end subroutine format_print_separator_abstract_out_file
!
!
   subroutine format_print_vector_abstract_out_file(the_file, name_, dim_, vector, fs, columns, transpose_)
!!
!!    Format print vector
!!    Written by Tor S. Haugland, Oct 2019
!!
!!    Prints a given vector using formats and columns.
!!    Based on the now deprecated print_vector by Eirik F. Kjønstad.
!!
!!    name_:      Name of vector
!!    dim_:       Number of elements to be printed
!!    vector:     Vector to print
!!
!!       OPTIONAL
!!    fs:         Format string for the numbers.         Default '(f13.8)'
!!    columns:    Number of columns to print per row.    Default: 4
!!    transpose_: Transpose the vector printing.         Default: .false.
!!
!!
!!    By default vector will be printed vertically in 4 columns,
!!       1 vector(1)   4 vector(4)   7 vector(7)  10 vector(10)
!!       2 vector(2)   5 vector(5)   8 vector(8)
!!       3 vector(3)   6 vector(6)   9 vector(9)
!!
      implicit none
!
      class(abstract_out_file),  intent(in)           :: the_file
      character(len=*),          intent(in)           :: name_
      integer,                   intent(in)           :: dim_
      real(dp), dimension(dim_), intent(in)           :: vector
!
      character(len=*),          intent(in), optional :: fs
      integer,                   intent(in), optional :: columns
      logical,                   intent(in), optional :: transpose_
!
      integer :: I, index_
      integer :: form_length, line_length, int_length
      integer :: n_columns, n_rows, row, col
!
      character(len=20) :: form_string, int_string, index_format
!
      logical :: adv, transpose_local
!
!     Set defaults
!
      if (present(fs)) then
         form_string = trim(fs)
      else
         form_string = '(f13.8)'
      endif
!
      if (present(columns)) then
         n_columns = columns
      else
         n_columns = 5
      endif
!
      if (present(transpose_)) then
         transpose_local = transpose_
      else
         transpose_local = .false.
      endif
!
!     Initialize formatting variables
!
      form_length = the_file%get_format_length(form_string)
!
      write(int_string, '(i0)') dim_
      int_length = len_trim(int_string)
!
      write(index_format,'(a,i0,a)') '(i', int_length, ')'
!
      line_length = (3 + form_length + int_length) * n_columns - 1
!
!     Print header
!
      call the_file%format_print(name_, fs='(/t3,a/)')
!
      call the_file%format_print_separator(line_length, symbol='-')
!
!     Print content
!
      n_rows = 1 + (dim_-1) / n_columns ! integer division
!
      do I = 1, n_rows * n_columns
!
         if (.not. transpose_local) then
!
!           Find row and column of element
!
            row = (I-1) / n_columns     ! which row
            col = mod(I - 1, n_columns) ! which column
!
            index_ = 1 + row + col * n_rows
            adv = ( mod(I,n_columns) == 0 .or. index_ + n_rows > dim_)
!
         else
!
            index_ = I
            adv = ( mod(I,n_columns) == 0 .or. index_ == dim_)
!
         endif
!
!        When size(dim_) < n_rows * n_columns, the index can be out of bounds
!
         if (index_ > dim_) cycle
!
         call the_file%format_print( trim(index_format) // ' ' // trim(form_string), adv=adv, &
                                    ints=[index_], reals=[vector(index_)])
!
      enddo
!
      call the_file%format_print_separator(line_length, symbol='-')
!
   end subroutine format_print_vector_abstract_out_file
!
!
end module abstract_out_file_class
