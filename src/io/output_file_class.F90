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
!!    Output file class module
!!    Written by Rolf H. Myhre, May 2019
!!
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
   contains
!
      procedure :: open_                     => open_output_file
!
      procedure :: error_msg                 => error_msg_output_file
      procedure :: warning_msg               => warning_msg_output_file
!
      procedure, public :: printf            => printf_output_file
      procedure, public :: print_matrix      => print_matrix_output_file
      procedure, public :: print_separator   => print_separator_output_file
      procedure, public :: print_vector      => print_vector_output_file
!
      procedure, public :: set_global_print_level  => set_global_print_level_output_file 
      procedure, public :: set_local_print_level   => set_local_print_level_output_file  
      procedure, public :: reset_local_print_level => reset_local_print_level_output_file  
      procedure, public :: check_print_level       => check_print_level_output_file  
!
      procedure, public :: mute   => mute_output_file
      procedure, public :: unmute => unmute_output_file
!
      procedure, private :: should_print => should_print_output_file
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
   function new_output_file(name_) result(the_file)
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
      the_file%action_ = 'write'
!
      the_file%global_print_level='normal'
      the_file%local_print_level='normal'
      the_file%is_mute = .false.
!
      the_file%is_open = .false.
      the_file%unit = -1
!
   end function new_output_file
!
!
   subroutine set_global_print_level_output_file(the_file, print_level)
!!
!!    Set the global print level
!!    Written by Rolf H. Myhre, May 2019
!!
!!    Also sets the local print level
!!
      implicit none
!
      class(output_file) :: the_file
!
      character(len=*), intent(in) :: print_level
!
      if (trim(print_level) .eq. 'normal') then
         the_file%global_print_level = 'normal'
         the_file%local_print_level  = 'normal'
!
      elseif (trim(print_level) .eq. 'minimal') then
         the_file%global_print_level = 'minimal'
         the_file%local_print_level  = 'minimal'
!
      elseif (trim(print_level) .eq. 'verbose') then
         the_file%global_print_level = 'verbose'
         the_file%local_print_level  = 'verbose'
!
      elseif (trim(print_level) .eq. 'debug') then
         the_file%global_print_level = 'debug'
         the_file%local_print_level  = 'debug'
!
      else
         print *, 'Error: Print level not normal, minimal, verbose or debug'
         stop
      endif
!
   end subroutine set_global_print_level_output_file
!
!
   subroutine set_local_print_level_output_file(the_file, print_level)
!!
!!    Set the local print level
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(output_file) :: the_file
!
      character(len=*), intent(in) :: print_level
!
      if (trim(print_level) .eq. 'normal') then
         the_file%local_print_level='normal'
!
      elseif (trim(print_level) .eq. 'minimal') then
         the_file%local_print_level='minimal'
!
      elseif (trim(print_level) .eq. 'verbose') then
         the_file%local_print_level='verbose'
!
      elseif (trim(print_level) .eq. 'debug') then
         the_file%local_print_level='debug'
!
      else
         print *, 'Error: Print level not normal, minimal, verbose or debug'
         stop
      endif
!
   end subroutine set_local_print_level_output_file
!
!
   subroutine reset_local_print_level_output_file(the_file)
!!
!!    Set the local print level to global
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(output_file) :: the_file
!
      the_file%local_print_level = the_file%global_print_level
!
   end subroutine reset_local_print_level_output_file
!
!
   subroutine check_print_level_output_file(the_file)
!!
!!    Check if local and global print level is the same
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(output_file) :: the_file
!
      if (the_file%local_print_level .ne. the_file%global_print_level) then
         print *, 'Warning: global and local print levels are not the same'
      endif
!
   end subroutine check_print_level_output_file
!
!
   subroutine mute_output_file(the_file)
!!
!!    Set the file to mute
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(output_file) :: the_file
!
      the_file%is_mute = .true.
!
   end subroutine mute_output_file
!
!
   subroutine unmute_output_file(the_file)
!!
!!    Set the file to unmute
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(output_file) :: the_file
!
      the_file%is_mute = .false.
!
   end subroutine unmute_output_file
!
!
   subroutine open_output_file(the_file, position_)
!!
!!    Open the output file
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(output_file) :: the_file
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
         print *, 'Error: could not open eT output file '//trim(the_file%name_)//&
                             &'. Error message: '//trim(io_msg)
         stop 
!
      endif 
!
      the_file%is_open = .true.
!
   end subroutine open_output_file
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
   subroutine printf_output_file(the_file, string, pl, reals, ints, chars, logs, fs, ffs, lfs, ll, adv)
!!
!!    printf
!!    Written by Rolf H. Myhre, May 2019
!!
!!    Printf output_file wrapper that checks for print level and silence
!!    Prints any number of reals and integers formatted Python style.
!!
!!    pl:       print level
!!              compared to stored print level variable and four allowed levels
!!              'minimal' or 'm' : Will always be printed. Only banners, final results 
!!                                 like total energies or excitation energies, and solver 
!!                                 settings or other essential information 
!!              'normal' or 'n'  : Will normally be printed, for example convergence iterations
!!              'verbose' or 'v' : Will only be printed if verbose output is specified in input,
!!                                 for example extra norms and MO coefficients 
!!              'debug'          : Print information only useful for developers such as extra tests 
!!                                 and index dimensions
!!
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
!!    Example: 
!!    call output%printf('Energy (a.u.): (f19.12)', pl='minimal', reals=[wf%energy], fs='(/t6,a)')
!!
!!
      implicit none
!
      class(output_file), intent(in) :: the_file
!
!     Data to print
      character(len=*), intent(in)                          :: string 
      character(len=*), intent(in), optional                :: pl
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
      character(len=10) :: plvl
!
      if(present(pl)) then !Set local plvl to pl if present
         plvl = trim(pl)
      else !Else default
         plvl = 'normal'
      endif
!
      if (the_file%is_mute) then !File is muted, make a quiet return
!
         return 
!
      endif
!
!
      if (the_file%should_print(plvl)) then
!
         call the_file%format_print(string, reals, ints, chars, logs, fs, ffs, lfs, ll, adv)
!
!        Flush if not a verbose print
         if (trim(plvl) .ne. 'verbose' .and. trim(plvl) .ne. 'v') then 
!
            call the_file%flush_()
!
         endif 
!
      endif
      
!
   end subroutine printf_output_file
!
!
   function should_print_output_file(the_file, print_level) result(should_print)
!!
!!    Should print
!!
!!    Written by Rolf H. Myhre, Oct. 2019
!!
!!    Checks if print_level is valid and compares it to local_print_level
!!
      implicit none
!
      class(output_file), intent(in) :: the_file
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
         if (the_file%local_print_level .ne. 'minimal') then
            should_print = .true.
         endif
!
!
!     Print if print_level is verbose and local print level is verbose or debug
      elseif ((trim(print_level) .eq. 'verbose') .or. (trim(print_level) .eq. 'v' )) then
!        
         if ((the_file%local_print_level .eq. 'verbose') .or. &
             (the_file%local_print_level .eq. 'debug')) then
            should_print = .true.
         endif
!
!
!     Print if print_level is debug and local print level is debug
      elseif (trim(print_level) .eq. 'debug') then
!
         if (the_file%local_print_level .eq. 'debug') then
            should_print = .true.
         endif
!
!
      else
!
         print *, 'Error: '//trim(print_level)// 'is not an acceptable print level'
         error stop 
!
      endif
!
!
   end function should_print_output_file
!
!
   subroutine print_matrix_output_file(the_file, pl, name_, matrix, dim_1, dim_2, fs, columns)
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
      class(output_file), intent(in) :: the_file
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
!
      if (the_file%is_mute) then !File is muted, make a quiet return
!
         return 
!
      endif
!
!
      if (the_file%should_print(pl)) then
!
         call the_file%format_print_matrix(name_, matrix, dim_1, dim_2, fs, columns)
!
         call the_file%flush_()
!
      endif
      
!
   end subroutine print_matrix_output_file
!
!
   subroutine print_separator_output_file(the_file, pl, n, symbol, fs)
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
      class(output_file), intent(in)         :: the_file
!
      character(len=*), intent(in)           :: pl
!
      integer, intent(in)                    :: n
!
      character, intent(in), optional        :: symbol
!
      character(len=*), optional, intent(in) :: fs
!
      if (the_file%is_mute) then !File is muted, make a quiet return
!
         return 
!
      endif
!
!
      if (the_file%should_print(pl)) then
!
         call the_file%format_print_separator(n, symbol, fs)
!
         call the_file%flush_()
!
      endif
      
!
   end subroutine print_separator_output_file
!
!
   subroutine print_vector_output_file(the_file, pl, name_, dim_, vector, fs, columns, transpose_)
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
      class(output_file),        intent(in)           :: the_file
      character(len=*),          intent(in)           :: pl
      character(len=*),          intent(in)           :: name_
      integer,                   intent(in)           :: dim_
      real(dp), dimension(dim_), intent(in)           :: vector
!
      character(len=*),          intent(in), optional :: fs
      integer,                   intent(in), optional :: columns
      logical,                   intent(in), optional :: transpose_
!
      if (the_file%is_mute) then !File is muted, make a quiet return
!
         return 
!
      endif
!
!
      if (the_file%should_print(pl)) then
!
         call the_file%format_print_vector(name_, dim_, vector, fs, columns, transpose_)
!
         call the_file%flush_()
!
      endif
      
!
   end subroutine print_vector_output_file
!
!
end module output_file_class
