!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2022 the authors of eT
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
!! Output file class module
!! Written by Rolf H. Myhre and Alexander C. Paul, 2019-2022
!!
!
   use kinds
   use abstract_out_file_class
!
   type, extends(abstract_out_file) :: output_file
!
!
      character(len=10), private :: global_print_level
      character(len=10), private :: local_print_level
      logical, private :: is_mute
!
      integer :: warning_counter
!
   contains
!
      procedure, public :: error_msg               => error_msg_output_file
      procedure, public :: warning_msg             => warning_msg_output_file
      procedure, public :: check_for_warnings      => check_for_warnings_output_file
!
      procedure, public :: printf                  => printf_output_file
      procedure, public :: print_matrix            => print_matrix_output_file
      procedure, public :: print_separator         => print_separator_output_file
      procedure, public :: print_vector            => print_vector_output_file
      procedure, public :: print_table             => print_table_output_file
!
      procedure, public :: newline                 => newline_output_file
!
      procedure, public :: set_global_print_level  => set_global_print_level_output_file
      procedure, public :: set_local_print_level   => set_local_print_level_output_file
      procedure, public :: reset_local_print_level => reset_local_print_level_output_file
!
      procedure, public :: mute   => mute_output_file
      procedure, public :: unmute => unmute_output_file
!
      procedure, public :: initialize => initialize_output_file
!
      procedure, private :: should_print => should_print_output_file
!
!     Generic write routines without formatstring
!
      procedure, public :: write_ => write_r_1
!
   end type output_file
!
   interface output_file
!
      procedure new_output_file
!
   end interface output_file
!
contains
!
!
   function new_output_file(name_) result(this)
!!
!!    Output file constructer
!!    Written by Rolf H. Myhre May 2019
!!
!!    Output file is a formatted and sequantial file.
!!    Routine sets these, and sets the file name
!!
      implicit none
!
      type(output_file) :: this
!
      character(len=*), intent(in) :: name_
!
      call this%initialize(name_)
!
   end function new_output_file
!
!
   subroutine initialize_output_file(this, name_)
!!
!!    Initialize
!!    Written by Sarai D. Folkestad, Sep 2021
!!
      implicit none
!
      class(output_file), intent(inout) :: this
!
      character(len=*), intent(in) :: name_
!
      this%name_ = trim(name_)
!
      this%access_ = 'sequential'
      this%format_ = 'formatted'
      this%action_ = 'write'
      this%status_ = 'unknown'
!
      this%global_print_level = 'normal'
      this%local_print_level = 'normal'
      this%is_mute = .false.
!
      this%is_open = .false.
      this%unit_ = -1
!
      this%warning_counter = 0
!
   end subroutine initialize_output_file
!
!
   subroutine set_global_print_level_output_file(this, print_level)
!!
!!    Set the global print level
!!    Written by Rolf H. Myhre, May 2019
!!
!!    Also sets the local print level
!!
      implicit none
!
      class(output_file) :: this
!
      character(len=*), intent(in) :: print_level
!
      if (trim(print_level) .eq. 'normal') then
         this%global_print_level = 'normal'
         this%local_print_level  = 'normal'
!
      elseif (trim(print_level) .eq. 'minimal') then
         this%global_print_level = 'minimal'
         this%local_print_level  = 'minimal'
!
      elseif (trim(print_level) .eq. 'verbose') then
         this%global_print_level = 'verbose'
         this%local_print_level  = 'verbose'
!
      elseif (trim(print_level) .eq. 'debug') then
         this%global_print_level = 'debug'
         this%local_print_level  = 'debug'
!
      else
         error stop 'Print level not minimal, normal, verbose or debug'
      endif
!
   end subroutine set_global_print_level_output_file
!
!
   subroutine set_local_print_level_output_file(this, print_level)
!!
!!    Set the local print level
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(output_file) :: this
!
      character(len=*), intent(in) :: print_level
!
      if (trim(print_level) .eq. 'normal') then
         this%local_print_level='normal'
!
      elseif (trim(print_level) .eq. 'minimal') then
         this%local_print_level='minimal'
!
      elseif (trim(print_level) .eq. 'verbose') then
         this%local_print_level='verbose'
!
      elseif (trim(print_level) .eq. 'debug') then
         this%local_print_level='debug'
!
      else
         error stop 'Print level not normal, minimal, verbose or debug'
      endif
!
   end subroutine set_local_print_level_output_file
!
!
   subroutine reset_local_print_level_output_file(this)
!!
!!    Set the local print level to global
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(output_file) :: this
!
      this%local_print_level = this%global_print_level
!
   end subroutine reset_local_print_level_output_file
!
!
   subroutine mute_output_file(this)
!!
!!    Set the file to mute
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(output_file) :: this
!
      this%is_mute = .true.
!
   end subroutine mute_output_file
!
!
   subroutine unmute_output_file(this)
!!
!!    Set the file to unmute
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(output_file) :: this
!
      this%is_mute = .false.
!
   end subroutine unmute_output_file
!
!
   subroutine error_msg_output_file(this, error_specs, &
                                    reals, complexs, ints, chars, logs, fs, ffs, ll, padd)
!!
!!    Error message
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Modified by Alexander C. Paul, Nov 2019
!!    Uses format_print for the error message
!!
!!    Modified by Rolf H. Myhre, Nov 2019
!!    Updated to reflect changes in format_print, plus update of documentation
!!
!!    error_specs: String of character that should be printed as error message,
!!                 including formatting of reals and integers
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
!!    Note: The number of characters to print per line will typically be between ll and
!!    ll + padd minus the number of blank spaces specified by the t specifier in the format
!!    string, assuming there are enough characters to print.
!!
      implicit none
!
      class(output_file) :: this
!
      character(len=*) :: error_specs
!
      real(dp)        , dimension(:), intent(in), optional  :: reals
      complex(dp)     , dimension(:), intent(in), optional  :: complexs
      integer         , dimension(:), intent(in), optional  :: ints
      character(len=*), dimension(:), intent(in), optional  :: chars
      logical         , dimension(:), intent(in), optional  :: logs
!
      character(len=*), optional, intent(in)                :: fs
      character(len=*), optional, intent(in)                :: ffs
!
      integer, optional, intent(in)                         :: ll
      integer, optional, intent(in)                         :: padd
!
      character(len=20) :: ff_string
      character(len=20) :: f_string
!
!     Default format: New line with t3 - Error message aligned after the colon
!
!     Format for the first line
      if(present(ffs)) then
         ff_string = ffs
      else
         ff_string = '(/t3,a)'
      endif
!
!     Format for the core message
      if(present(fs)) then
         f_string = fs
      else
         f_string = '(t10,a)'
      endif
!
!     Option advance from format_print shall always be true for errors
!
      call this%format_print('Error: ' // trim(error_specs),     &
                                 reals, complexs, ints, chars, logs, &
                                 fs = f_string, ffs = ff_string,     &
                                 ll = ll, padd = padd,               &
                                 adv = .true.)
!
      call this%flush_()
!
!     Cannot be an error stop because of the unit tests
      stop "Something went wrong, check the .out file"
!
   end subroutine error_msg_output_file
!
!
   subroutine warning_msg_output_file(this, warning_specs, &
                                      reals, complexs, ints, chars, logs, fs, ffs, ll, padd)
!!
!!    Warning message
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Modified by Alexander C. Paul, Nov 2019
!!    Uses format_print for the warning message
!!
!!    Modified by Rolf H. Myhre, Nov 2019
!!    Updated to reflect changes in format_print, plus update of documentation
!!
!!    error_specs: String of character that should be printed as error message,
!!                 including formatting of reals and integers
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
!!    Note: The number of characters to print per line will typically be between ll and
!!    ll + padd minus the number of blank spaces specified by the t specifier in the format
!!    string, assuming there are enough characters to print.
!!
      implicit none
!
      class(output_file) :: this
!
      character(len=*) :: warning_specs
!
      real(dp)        , dimension(:), intent(in), optional  :: reals
      complex(dp)     , dimension(:), intent(in), optional  :: complexs
      integer         , dimension(:), intent(in), optional  :: ints
      character(len=*), dimension(:), intent(in), optional  :: chars
      logical         , dimension(:), intent(in), optional  :: logs
!
      character(len=*), optional, intent(in)                :: fs
      character(len=*), optional, intent(in)                :: ffs
!
      integer, optional, intent(in)                         :: ll
      integer, optional, intent(in)                         :: padd
!
      character(len=20) :: ff_string
      character(len=20) :: f_string
!
      this%warning_counter = this%warning_counter + 1
!
!     Default format: New line with t3 - Warning message aligned after the colon
!
!     Format for the first line
      if(present(ffs)) then
         ff_string = ffs
      else
         ff_string = '(/t3,a)'
      endif
!
!     Format for the core message
      if(present(fs)) then
         f_string = fs
      else
         f_string = '(t12,a)'
      endif
!
!     Option advance from format_print shall always be true for warnings
!
      call this%format_print('Warning: ' // trim(warning_specs), &
                                 reals, complexs, ints, chars, logs, &
                                 fs = f_string, ffs = ff_string,     &
                                 ll = ll, padd = padd,               &
                                 adv = .true.)
!
      call this%flush_()
!
   end subroutine warning_msg_output_file
!
!
   subroutine check_for_warnings_output_file(this)
!!
!!    Check for warnings
!!    Written by Alexander C. Paul, Nov 2019
!!
      implicit none
!
      class(output_file) :: this
!
      if (this%warning_counter == 1) then
!
         call this%printf('m', ':: There was 1 warning during the execution of eT. ::', &
                              fs='(/t3,a)')
!
      else if(this%warning_counter > 1) then
!
         call this%printf('m', ':: There were (i0) warnings during the execution of eT. ::', &
                              ints=[this%warning_counter], fs='(/t3,a)')
!
      else if(this%warning_counter < 0) then
!
         call this%error_msg('Negative number of warning messages.' // &
                             'Something went horribly wrong.')
!
      end if
!
   end subroutine check_for_warnings_output_file
!
!
   subroutine printf_output_file(this, pl, string, &
                                 reals, complexs, ints, chars, logs, &
                                 fs, ffs, ll, padd, adv)
!!
!!    printf
!!    Written by Rolf H. Myhre, May 2019
!!
!!    Printf output_file wrapper that checks for print level and silence before calling
!!    format_print which prints any number of reals, integers, characters and logicals
!!    formatted Python style.
!!
!!    pl:      print level
!!             character string that is compared to the print level of the file with four
!!             allowed levels:
!!             'minimal' or 'm' : Will always be printed. Only banners, final results
!!                                like total energies or excitation energies, and solver
!!                                settings or other essential information
!!             'normal' or 'n'  : Will normally be printed, for example convergence iterations
!!             'verbose' or 'v' : Will only be printed if verbose output is specified in input,
!!                                for example extra norms and MO coefficients
!!             'debug'          : Print information only useful for developers such as extra tests
!!                                and index dimensions
!!
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
!!    Note: The number of characters to print per line will typically be between ll and
!!    ll + padd minus the number of blank spaces specified by the t specifier in the format
!!    string, assuming there are enough characters to print.
!!
!!    Example:
!!    call output%printf('(a0) ground state energy (a.u): (f19.12)', &
!!                       pl='minimal', reals=[wf%energy], chars=[wf%name_], fs='(/t6,a)')
!!
      implicit none
!
      class(output_file), intent(in) :: this
!
!     Data to print
      character(len=*), intent(in)                          :: string
      character(len=*), intent(in)                          :: pl
      real(dp), dimension(:), intent(in), optional          :: reals
      complex(dp), dimension(:), intent(in), optional       :: complexs
      integer, dimension(:), intent(in), optional           :: ints
      character(len=*), dimension(:), intent(in), optional  :: chars
      logical, dimension(:), intent(in), optional           :: logs
!
!     Parameters to pass to long_string_print
      integer, intent(in), optional                         :: ll
      character(len=*), optional, intent(in)                :: fs
      character(len=*), optional, intent(in)                :: ffs
      integer, optional, intent(in)                         :: padd
!
      logical, intent(in), optional :: adv
!
      if (this%is_mute) return
!
      if (this%should_print(pl)) then
!
         call this%preformat_print(string, reals, complexs, ints, chars, logs, &
                                       fs, ffs, ll, padd, adv)
!
         call this%flush_()
!
      endif

!
   end subroutine printf_output_file
!
!
   function should_print_output_file(this, print_level) result(should_print)
!!
!!    Should print
!!
!!    Written by Rolf H. Myhre, Oct. 2019
!!
!!    Checks if print_level is valid and compares it to local_print_level
!!
      implicit none
!
      class(output_file), intent(in) :: this
!
      character(len=*), intent(in)   :: print_level
!
      logical :: should_print
!
!     Default value is false
      should_print = .false.
!
!     Always print if print_level is minimal
      if ((trim(print_level) .eq. 'minimal') .or. (trim(print_level) .eq. 'm' )) then
!
         should_print = .true.
!
!
!     Print if print_level is normal and local print level not minimal
      elseif ((trim(print_level) .eq. 'normal') .or. (trim(print_level) .eq. 'n' )) then
!
         if (this%local_print_level .ne. 'minimal') then
            should_print = .true.
         endif
!
!
!     Print if print_level is verbose and local print level is verbose or debug
      elseif ((trim(print_level) .eq. 'verbose') .or. (trim(print_level) .eq. 'v' )) then
!
         if ((this%local_print_level .eq. 'verbose') .or. &
             (this%local_print_level .eq. 'debug')) then
            should_print = .true.
         endif
!
!
!     Print if print_level is debug and local print level is debug
      elseif (trim(print_level) .eq. 'debug') then
!
         if (this%local_print_level .eq. 'debug') then
            should_print = .true.
         endif
!
!
      else
!
         print '(3a)', "Error: ", trim(print_level), " is not an acceptable print level"
         error stop
!
      endif
!
!
   end function should_print_output_file
!
!
   subroutine print_matrix_output_file(this, pl, name_, matrix, dim_1, dim_2, fs, columns)
!!
!!    Print matrix
!!
!!    Written by Rolf H. Myhre, Oct. 2019
!!
!!    Calls format_print_matrix if appropriate print level
!!
!!    name_: Name to be printed above the matrix
!!
!!    Matrix to be printed with dimension dim_1 x dim_2
!!
!!    fs:      Optional format string for numbers, default is (f19.12)
!!    columns: Optional integer specifying number of columns to print per line, default is 5
!!
      implicit none
!
      class(output_file), intent(in) :: this
!
      character(len=*), intent(in)                    :: pl
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
      if (this%is_mute) return
!
      if (this%should_print(pl)) then
!
         call this%format_print_matrix(name_, matrix, dim_1, dim_2, fs, columns)
!
         call this%flush_()
!
      endif

!
   end subroutine print_matrix_output_file
!
!
   subroutine print_separator_output_file(this, pl, n, symbol, fs)
!!
!!    Print separator
!!
!!    Written by Rolf H. Myhre, Oct. 2019
!!
!!    Calls format_print_separator if appropriate print level
!!
!!    n: Number of symbols to print
!!    symbol: optional, what symbol to print, default is '='
!!    fs: optional format string
!!
      implicit none
!
      class(output_file), intent(in)         :: this
!
      character(len=*), intent(in)           :: pl
!
      integer, intent(in)                    :: n
!
      character, intent(in), optional        :: symbol
!
      character(len=*), optional, intent(in) :: fs
!
      if (this%is_mute) return
!
      if (this%should_print(pl)) then
!
         call this%format_print_separator(n, symbol, fs)
!
         call this%flush_()
!
      endif

!
   end subroutine print_separator_output_file
!
!
   subroutine print_vector_output_file(this, pl, name_, dim_, vector, fs, columns, transpose_)
!!
!!    Print vector
!!    Written by Tor S. Haugland, Oct 2019
!!
!!    Calls format_print_vector if appropriate print level
!!
!!    Based on print_matrix_output_file by Rolf H. Myhre.
!!
!!    name_     : name to be printed above the vector
!!    dim       : size(vector)
!!    vector    : vector to print
!!
!!       OPTIONAL
!!    fs        : format string, default '(f12.8)'
!!    columns   : number of columns to print per row, default 4
!!    transpose : transpose the printed vector
!!
      implicit none
!
      class(output_file),        intent(in)           :: this
      character(len=*),          intent(in)           :: pl
      character(len=*),          intent(in)           :: name_
      integer,                   intent(in)           :: dim_
      real(dp), dimension(dim_), intent(in)           :: vector
!
      character(len=*),          intent(in), optional :: fs
      integer,                   intent(in), optional :: columns
      logical,                   intent(in), optional :: transpose_
!
      if (this%is_mute) return
!
      if (this%should_print(pl)) then
!
         call this%format_print_vector(name_, dim_, vector, fs, columns, transpose_)
!
         call this%flush_()
!
      endif

!
   end subroutine print_vector_output_file
!
!
   subroutine newline_output_file(this, pl)
!!
!!    Newline
!!    Written by Sarai D. Folkestad, 2021
!!
      implicit none

      class(output_file), intent(inout) :: this
      character(len=*),   intent(in)    :: pl
!
      call this%printf(pl, '')
!
   end subroutine newline_output_file
!
!
   subroutine print_table_output_file(this, top_left_corner, row_labels, column_labels, &
                                      data_, n_rows, n_columns)
!!
!!    Print table
!!    Written by Alexander C. Paul, June 2020
!!
!!    Prints table of data with n_rows and n_columns
!!
      implicit none
!
      class(output_file), intent(inout) :: this
!
      integer, intent(in) :: n_rows, n_columns
!
      real(dp), dimension(n_rows, n_columns), intent(in) :: data_
!
      character(len=5), intent(in) :: top_left_corner
      character(len=5), dimension(n_rows), intent(in) :: row_labels
      character(len=18), dimension(n_columns), intent(in) :: column_labels
!
      integer :: n_tables, n_printed, n_to_print, ll
      integer :: row, table, first, last
!
      character(len=100) :: line
!
!     Print in tables with 4 columns
      n_tables = (n_columns - 1)/4 + 1
!
      n_printed = 0
!
      do table = 1, n_tables
!
         n_to_print = min(n_columns - n_printed, 4)
!
         first = n_printed + 1
         last = first + n_to_print - 1
!
         ll = 6 + 18*n_to_print
!
!        Print header of the table (first_label  Label_1  Label_2 ...)
!
         write(line,'(a,1x,i0,a)') top_left_corner, n_to_print, '(a18)'
!
         call this%printf('m', trim(line), fs='(/t6,a)', ll=95, &
                            chars=[column_labels(first : last)])
!
         call this%print_separator('m', ll,'-', fs='(t6,a)')
!
!        Print rows for the subsets
!
         do row = 1, n_rows
!
            write(line,'(a5,1x,i0,a)') row_labels(row), n_to_print, '(f18.10)'
!
            call this%printf('m', trim(line), fs='(t6,a)', ll=95, &
                               reals=[data_(row, first : last)])
!
         end do
!
         call this%print_separator('m', ll,'-', fs='(t6,a)')
!
         n_printed = n_printed + n_to_print
!
      end do
!
   end subroutine print_table_output_file
!
!
   subroutine write_r_1(this, array, n)
!!
!!    Write, real(dp) 1 dimensional array
!!    Written by Rolf H. Myhre, May 2019
!!
!!    Uses default format, for user-defined format see printf
!!
      implicit none
!
      class(output_file), intent(inout) :: this
!
      integer, intent(in)                :: n
      real(dp), dimension(n), intent(in) :: array
!
      integer :: io_status
      character(len=100) :: io_message
!
      write(this%unit_, *, iostat=io_status, iomsg=io_message) array
      call this%check_io_status(io_status, io_message, task="write")
!
   end subroutine write_r_1
!
!
end module output_file_class
