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
!
   contains
!
      procedure :: open_                     => open_abstract_out_file
      procedure :: close_                    => close_abstract_out_file
      procedure :: flush_                    => flush_abstract_out_file
!
      procedure, public :: format_print            => format_print_abstract_out_file
!     
      procedure, public :: format_print_matrix     => format_print_matrix_abstract_out_file
!
      procedure, public :: format_print_separator  => format_print_separator_abstract_out_file
!
      procedure :: long_string_print         => long_string_print_abstract_out_file
!
      procedure, nopass, private :: get_format_length
!
   end type abstract_out_file
!
contains
!
!
   subroutine open_abstract_out_file(the_file, position_)
!!
!!    Open the abstract_out file
!!    Written by Rolf Heilemann Myhre, May 2019
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
      open(newunit=the_file%unit, file=the_file%name_, access=the_file%access_, &
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
!!    Written by Rolf Heilemann Myhre, May 2019
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
      end if
!
      close(the_file%unit, iostat=io_error, iomsg=io_msg, status='keep')
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
      the_file%unit = -1
!
   end subroutine close_abstract_out_file
!
!
   subroutine flush_abstract_out_file(the_file)
!!
!!    Flush the output file
!!    Written by Rolf Heilemann Myhre, May 2019
!!
      implicit none
!
      class(abstract_out_file) :: the_file
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      flush(the_file%unit, iostat=io_error, iomsg=io_msg)
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
   subroutine format_print_abstract_out_file(the_file, string, reals, ints, chars, logs, fs, ffs, lfs, ll, adv)
!!
!!    Format print 
!!    Written by Rolf Heilemann Myhre, May 2019
!!
!!    Prints any number of reals and integers formatted Python style.
!!
!!    string:   String of character that should be printed, 
!!              including formatting of reals and integers 
!!    reals:    Array of real(dp) to print - in the order specified by string 
!!    ints:     Array of integers to print - in the order specified by string 
!!    chars:    Array of strings to print - in the order specified by string
!!              Note that all the strings must be of same length in Fortran
!!
!!    fs:       Specifies the format of the entire string, e.g. fs='(/t6,a)' gives 
!!              a new line, then indentation 5, then the value of 'string'
!!              with reals and integers as specified. Default: '(t3,a)'
!!    ffs:      Specifies the format of the first printed line if different from fs. 
!!              Default: same as fs
!!    lfs:      Specifies the format of the last printed line if different from fs. 
!!              Default: same as fs
!!    ll:       Integer specifying number of characters per line of print.
!!    adv:      Logical specifies whether advance is 'yes' or 'no', default = 'yes'
!!
!!    Example: call output%printf('Energy (a.u.): (f19.12)', reals=[wf%energy], fs='(/t6,a)')
!!
      implicit none
!
      class(abstract_out_file), intent(in) :: the_file
!
!     Data to print
      character(len=*), intent(in)                          :: string 
      real(dp), dimension(:), intent(in), optional          :: reals 
      integer, dimension(:), intent(in), optional           :: ints
      character(len=*), dimension(:), intent(in), optional  :: chars
      logical, dimension(:), intent(in), optional           :: logs
!
!     Parameters to pass to long_string_print
      integer, intent(in), optional                         :: ll
      character(len=*), optional, intent(in)                :: fs
      character(len=*), optional, intent(in)                :: ffs
      character(len=*), optional, intent(in)                :: lfs
!
      logical, intent(in), optional :: adv
!
      character(len=1000)  :: pstring 
      character(len=20)    :: fstring 
!
      integer :: i, p_pos, int_check, i_err, add_pos
      integer :: int_len, real_len, log_len, char_len, string_len
      integer :: int_count, real_count, log_count, char_count
      integer :: print_pos, printed
!
      integer           :: l_length
      character(len=20) :: f_string
      character(len=20) :: ff_string
      character(len=20) :: lf_string
      logical           :: adva
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
!     Check how many logicals are present
      if(present(logs)) then
         log_len = size(logs)
      else
         log_len = 0
      endif
!
!     Check how many chars are present
      if(present(chars)) then
         char_len = size(chars)
      else
         char_len = 0
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
!     First line format string to send to long_string_print
      if(present(ffs)) then
         ff_string = ffs
      else
         ff_string = f_string
      endif
!
!     Last line format string to send to long_string_print
      if(present(lfs)) then
         lf_string = lfs
      else
         lf_string = f_string
      endif
!
!     Advance write?
      if(present(adv)) then
         adva = adv
      else
         adva = .True.
      endif
!
      real_count = 0
      int_count  = 0
      log_count  = 0
      char_count  = 0
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
                  if (real_count .gt. real_len) print *, 'Not enough reals in printf'
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
                  if (int_count .gt. int_len) print *, 'Not enough ints in printf'
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
!  
!           Is ( followed by l or L?
            elseif(string(i+1:i+1) .eq. "l" .or. string(i+1:i+1) .eq. "L") then
!
!              Is it followed by a number, if so, assume a format string
               read(string(i+2:),'(i1)',iostat=i_err) int_check
               if (i_err .eq. 0) then
!
                  log_count = log_count + 1
                  if (log_count .gt. log_len) print *, 'Not enough ints in printf'
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
                  if (fstring(2:3) .eq. 'l0') then
                     fstring = '(a)'
                     add_pos = len_trim(chars(char_count))
                  else
                     fstring(2:2) = 'a'
                     add_pos = the_file%get_format_length(fstring)
                  endif
!
                  if(logs(log_count)) then
                     write(pstring(print_pos:), fstring) 'True'     
                  else
                     write(pstring(print_pos:), fstring) 'False'     
                  endif
!
!                 Set next position to print
                  print_pos = print_pos + add_pos
!
               endif
!  
!           Is ( followed by a or A?
            elseif(string(i+1:i+1) .eq. "a" .or. string(i+1:i+1) .eq. "A") then
!
!              Is it followed by a number, if so, assume a format string
               read(string(i+2:),'(i1)',iostat=i_err) int_check
               if (i_err .eq. 0) then
!
!
                  char_count = char_count + 1
                  if (char_count .gt. char_len) print *,  'Not enough chars in printf'
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
!                 Copy format string to fstring and check if a0
                  fstring = string(p_pos:i)
                  if (fstring(2:3) .eq. 'a0') then
                     fstring = '(a)'
                     add_pos = len_trim(chars(char_count))
                  else
                     add_pos = the_file%get_format_length(fstring)
                  endif
!
                  write(pstring(print_pos:), fstring) trim(chars(char_count))
!
!                 Set next position to print
                  print_pos = print_pos + add_pos
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
      call the_file%long_string_print(pstring, fs=f_string, ffs=ff_string, lfs=lf_string, & 
                                      ll=l_length, adv=adva)
!
   end subroutine format_print_abstract_out_file
!
!  
   subroutine long_string_print_abstract_out_file(the_file, string, fs, colons, &
                                            ffs, lfs, ll, adv)
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
      class(abstract_out_file), intent(in) :: the_file
!
      character(len=*), intent(in) :: string
      character(len=*), intent(in), optional :: fs, ffs, lfs
      integer, intent(in), optional :: ll
      logical, intent(in), optional :: colons
      logical, intent(in), optional :: adv
!
      character(90)  :: temp
      integer        :: l, l_left, lines, l_length
      integer        :: i,j, padd, printed 
      character(20)  :: f_s,fstring,ffstring,lfstring
      logical        :: col 
!
      character(len=3) :: adva
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
      adva = 'yes'
      if(present(adv)) then
         if(adv) then
            adva = 'yes'
         else
            adva = 'no'
         endif
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
               write(the_file%unit,f_s, advance=trim(adva)) ':: '//string(printed:l)
            else
               write(the_file%unit,f_s, advance=trim(adva)) string(printed:l)
            endif
            printed = l
         endif
!
      enddo
!
      call the_file%flush_()
!
   end subroutine long_string_print_abstract_out_file
!
!
   function get_format_length(fstring) result(length)
!
!!    Get printed length of format
!!    Written by Rolf H. Myhre, Sep 2019
!!
!!    Figure out how many printed characters a format string corresponds to
!
      implicit none
!
      character(len=*), intent(in)  :: fstring
!
      integer i, length, string_len, stat
!
      string_len = len_trim(fstring)
!
!     If character, integer or logical
      if(fstring(1:1) .eq. '(') then 
         if(fstring(2:2) .eq. 'a' .or. fstring(2:2) .eq. 'A' .or. &
            fstring(2:2) .eq. 'l' .or. fstring(2:2) .eq. 'L' .or. &
            fstring(2:2) .eq. 'i' .or. fstring(2:2) .eq. 'I') then
!
            i = 3
            do while (fstring(i+1:i+1) .ne. ")" .and. i .lt. string_len)
               i = i + 1
            enddo
!
            read(fstring(3:i),*,iostat=stat) length
!
            if(stat .ne. 0) then
               print *, fstring//' is not an acceptable format string'
            endif 
!
!        If float
         elseif(fstring(2:2) .eq. 'f' .or. fstring(2:2) .eq. 'F' .or. &
                fstring(2:2) .eq. 'e' .or. fstring(2:2) .eq. 'E') then
!
            i = 3
            do while (fstring(i+1:i+1) .ne. "." .and. i .lt. string_len)
               i = i + 1
            enddo
!
            read(fstring(3:i),*,iostat=stat) length
!
            if(stat .ne. 0) then
               print *, fstring//' is not an acceptable format string'
            endif 
!
         endif
      else
         print *, fstring//' is not an acceptable format string'
      endif
!
   end function get_format_length
!
!
   subroutine format_print_matrix_abstract_out_file(the_file, name_, matrix, dim_1, dim_2, fs, columns)
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
      call the_file%format_print_separator(line_length+1)
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
         if ((i .eq.  n_prints) .and. (i .ne. 1)) then
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
               reals_to_print(j) = matrix(k,j)
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
      call the_file%format_print_separator(line_length+1)
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
!!    symbol: optional, what symbol to print, default is '='
!!    fs: optional format string
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
      if (present(symbol)) then
         sym = symbol
      else
         sym = '='
      endif
!
      separator_line = repeat(sym, n)
!
      call the_file%format_print(separator_line, ll=n, fs=fs)
!
   end subroutine format_print_separator_abstract_out_file
!
!
end module abstract_out_file_class
