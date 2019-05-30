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
   use section_class, only : section
   use string_utilities 
   use abstract_file_class, only : abstract_file
!
   type, extends(abstract_file) :: input_file 
!
      type(section), allocatable :: sections(:)
!
      character(len=21), allocatable :: rf_wfs(:)
      character(len=21), allocatable :: cc_wfs(:)
!
   contains
!
      procedure :: open_                                                => open__input_file
      procedure :: close_                                               => close__input_file
!
      procedure :: check_for_errors                                     => check_for_errors_input_file
      procedure :: requested_section                                    => requested_section_input_file
      procedure :: requested_keyword_in_section                         => requested_keyword_in_section_input_file
      procedure :: get_n_elements_for_keyword_in_section                => get_n_elements_for_keyword_in_section_input_file
      procedure :: get_array_for_keyword_in_section                     => get_array_for_keyword_in_section_input_file
      procedure :: get_n_atoms                                          => get_n_atoms_input_file
      procedure :: get_geometry                                         => get_geometry_input_file
      procedure :: get_reference_wf                                     => get_reference_wf_input_file
      procedure :: get_cc_wf                                            => get_cc_wf_input_file
!
      procedure :: requested_reference_calculation                      => requested_reference_calculation_input_file
      procedure :: requested_cc_calculation                             => requested_cc_calculation_input_file
!
      procedure, private :: get_string_keyword_in_section_wo_safety     => get_string_keyword_in_section_wo_safety_input_file
      procedure, private :: move_to_section                             => move_to_section_input_file
      procedure, private :: check_section_for_illegal_keywords          => check_section_for_illegal_keywords_input_file
      procedure, private :: check_for_illegal_sections                  => check_for_illegal_sections_input_file
      procedure, private :: print_sections                              => print_sections_input_file
      procedure, private :: read_adjustl_lower                          => read_adjustl_lower_input_file
!
      procedure, nopass, private :: string_is_comment                   => string_is_comment_input_file
      procedure, nopass, private :: extract_keyword_from_string         => extract_keyword_from_string_input_file
      procedure, nopass, private :: extract_keyword_value_from_string   => extract_keyword_value_from_string_input_file
!
      generic :: get_keyword_in_section                                 => get_integer_keyword_in_section_input_file,   &
                                                                           get_integer8_keyword_in_section_input_file,  &
                                                                           get_string_keyword_in_section_input_file,    &
                                                                           get_dp_keyword_in_section_input_file
!
      generic :: get_required_keyword_in_section                        => get_required_string_keyword_in_section_input_file,    &
                                                                           get_required_integer_keyword_in_section_input_file,   &
                                                                           get_required_integer8_keyword_in_section_input_file,  &
                                                                           get_required_dp_keyword_in_section_input_file
!
      procedure :: get_integer_keyword_in_section_input_file
      procedure :: get_integer8_keyword_in_section_input_file
      procedure :: get_string_keyword_in_section_input_file
      procedure :: get_dp_keyword_in_section_input_file
      procedure :: get_required_string_keyword_in_section_input_file
      procedure :: get_required_integer_keyword_in_section_input_file
      procedure :: get_required_integer8_keyword_in_section_input_file
      procedure :: get_required_dp_keyword_in_section_input_file
!
   end type input_file
!
   interface input_file
!
      procedure new_input_file
!
   end interface input_file
!
!
contains
!
!
   module function new_input_file(name_) result(the_file)
!!
!!    Initialize input file
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2018
!!
!!    Initializes the input file for opening by the disk manager
!!    and sets the valid sections and their keywords. 
!!
      implicit none
!
      type(input_file) :: the_file
!
      character(len=*), intent(in) :: name_
!
      integer :: k
!
      type(section) :: calculations
      type(section) :: system 
      type(section) :: memory 
      type(section) :: disk 
      type(section) :: cc_zop 
      type(section) :: cc_fop 
      type(section) :: method 
      type(section) :: solver_cholesky
      type(section) :: solver_hf
      type(section) :: solver_cc_gs
      type(section) :: solver_cc_es
      type(section) :: solver_cc_multipliers 
      type(section) :: active_atoms
!
!     Set input file name, access and format 
!
      the_file%name_ = name_
!
      the_file%access_ = 'sequential'
      the_file%format_ = 'formatted'
!
!     Set method section 
!
      the_file%rf_wfs = (/ 'hf                   ',   &
                           'uhf                  '    /)
!
      the_file%cc_wfs = (/ 'ccs                  ',   &
                           'mp2                  ',   &
                           'cc2                  ',   &
                           'lowmem-cc2           ',   &
                           'ccsd                 ',   &
                           'cc3                  '    /)
!
      method%name_    = 'method'
      method%required = .true.
!
      allocate(method%keywords(size(the_file%rf_wfs) + size(the_file%cc_wfs)))
!
      do k = 1, size(the_file%rf_wfs)
!
         method%keywords(k) = the_file%rf_wfs(k)
!
      enddo 
!
      do k = 1, size(the_file%cc_wfs)
!
         method%keywords(size(the_file%rf_wfs) + k) = the_file%cc_wfs(k)
!
      enddo 
!
!     Set other sections
!
      calculations%name_    = 'do'
      calculations%required = .true.
      calculations%keywords = (/ 'ground state         ',   &
                                 'excited state        ',   &
                                 'zop                  ',   &
                                 'fop                  ',   &
                                 'cholesky eri         ',   &
                                 'multipliers          '   /)
!
      system%name_    = 'system'
      system%required = .true.
      system%keywords = (/ 'name                 ',   &
                           'charge               ',   &
                           'multiplicity         '    /)
!
      memory%name_    = 'memory'
      memory%required = .false.
      memory%keywords = (/ 'available            ' /)
!
      disk%name_    = 'disk'
      disk%required = .false.
      disk%keywords = (/ 'available            ' /)
!
      cc_zop%name_    = 'cc zop'
      cc_zop%required = .false.
      cc_zop%keywords = (/ 'dipole               ',   &
                           'quadrupole           ' /)
!
      cc_fop%name_    = 'cc fop'
      cc_fop%required = .false.
      cc_fop%keywords = (/ 'dipole length        ',   &
                           'lr                   ',   &
                           'eom                  ' /)
!
      solver_cholesky%name_    = 'solver cholesky'
      solver_cholesky%required = .false.
      solver_cholesky%keywords = (/ 'threshold           ',    &
                                    'span                ',    &
                                    'batches             ',    &
                                    'qualified           ',    &
                                    'one center          ',    &
                                    'no vectors          '     /)
!
      solver_hf%name_    = 'solver hf'
      solver_hf%required = .false.
      solver_hf%keywords = (/ 'algorithm            ',   &
                              'energy threshold     ',   &
                              'gradient threshold   ',   &
                              'max iterations       ',   &
                              'diis dimension       ',   &
                              'restart              ',   &
                              'ao density guess     '    /)
!
      solver_cc_gs%name_    = 'solver cc gs'
      solver_cc_gs%required = .false.
      solver_cc_gs%keywords = (/ 'algorithm            ',   &
                                 'energy threshold     ',   &
                                 'omega threshold      ',   &
                                 'max micro iterations ',   &
                                 'rel micro threshold  ',   &
                                 'max iterations       ',   &
                                 'diis dimension       ',   &
                                 'restart              '    /)
!
      solver_cc_es%name_    = 'solver cc es'
      solver_cc_es%required = .false.
      solver_cc_es%keywords = (/ 'algorithm            ',   &
                                 'core excitation      ',   &
                                 'energy threshold     ',   &
                                 'residual threshold   ',   &
                                 'max iterations       ',   &
                                 'restart              ',   &
                                 'left eigenvectors    ',   &
                                 'right eigenvectors   ',   &
                                 'singlet states       ',   &
                                 'diis dimension       ',   &
                                 'max reduced dimension' /)
!
      solver_cc_multipliers%name_    = 'solver cc multipliers'
      solver_cc_multipliers%required = .false.
      solver_cc_multipliers%keywords = (/ 'algorithm            ',   &
                                          'threshold            ',   &
                                          'restart              ',   &
                                          'max iterations       '    /)
!
      active_atoms%name_    = 'active atoms'
      active_atoms%required = .false.
      active_atoms%keywords = (/ 'selection type       ', &
                                 'central atom         ', &
                                 'hf                   ', &
                                 'active basis         ' /)
!
!     Gather all sections into the file's section array 
!
      the_file%sections = [calculations,           &
                           system,                 &
                           memory,                 &
                           disk,                   &
                           method,                 &
                           cc_zop,                 &
                           cc_fop,                 &
                           solver_cholesky,        &
                           solver_hf,              &
                           solver_cc_gs,           &
                           solver_cc_es,           &
                           solver_cc_multipliers,  &
                           active_atoms]
!
   end function new_input_file
!
!
   subroutine open__input_file(the_file)
!!
!!    Open the input file
!!    Written by Rolf Heilemann Myhre, May 2019
!!
      implicit none
!
      class(input_file) :: the_file
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      open(newunit=the_file%unit, file=the_file%name_, access=the_file%access_, &
           action='read', status='unknown', form=the_file%format_, iostat=io_error, iomsg=io_msg)
!
      if (io_error /= 0) stop 'Error: could not open eT input file '//trim(the_file%name_)//&
                             &'error message: '//trim(io_msg)
!
      the_file%opened = .true.
!
   end subroutine open__input_file
!
!
   subroutine close__input_file(the_file)
!!
!!    Close the input file
!!    Written by Rolf Heilemann Myhre, May 2019
!!
      implicit none
!
      class(input_file) :: the_file
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      close(the_file%unit, iostat=io_error, iomsg=io_msg, status='keep')
!
      if (io_error /= 0) stop 'Error: could not close eT input file '//trim(the_file%name_)//&
                             &'error message: '//trim(io_msg)
!
      the_file%opened = .false.
!
   end subroutine close__input_file
!
!
   subroutine check_for_errors_input_file(the_file)
!!
!!    Check for errors
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
!!    Looks for errors the user may have made when writing the input file.
!!
      implicit none 
!
      class(input_file) :: the_file
!
      integer :: k 
!
!     Look for illegal sections 
!
      call the_file%check_for_illegal_sections()
!
!     For each legal section present, look for illegal keywords within that section 
!
      do k = 1, size(the_file%sections)
!
         if (the_file%requested_section(the_file%sections(k)%name_)) then
!
            call the_file%check_section_for_illegal_keywords(the_file%sections(k))
!
         else
!
            if (the_file%sections(k)%required) then 
!
               write(output%unit, '(/t3,a,a,a)') 'All calculations require the section "', trim(the_file%sections(k)%name_), &
                                                   '". It appears to be missing.'
!
               call output%error_msg('Something is wrong in the input file. See above.')
! 
            endif 
!
         endif
!
      enddo
!
   end subroutine check_for_errors_input_file
!
!
   subroutine check_for_illegal_sections_input_file(the_file)
!!
!!    Check for illegal sections 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
!!    Finds the sections in the input file and controls whether they are among 
!!    the valid sections in the program.
!!
!!    Tests only the end of each section. If the beginning of a section is incorrectly 
!!    specified, an error will occur when it checks the section for illegal keywords. 
!!
      implicit none 
!
      class(input_file) :: the_file
!
      character(len=200) :: line
!
      integer :: n_elements, k
!
      logical :: recognized 
!  
      rewind(the_file%unit)
!
      line = the_file%read_adjustl_lower()
!
      do while (trim(line) /= 'end geometry') 
!
         if (line(1 : 3) == 'end') then
!
!           Located the end of a section -
!           attempt to move to the beginning of that section => fails if inconsistent beginning and end
!
            call the_file%move_to_section(trim(adjustl(line(4 : 200))), n_elements)
!
!           Check whether section name is valid 
!
            recognized = .false.
!
            do k = 1, size(the_file%sections)
!
               if (trim(the_file%sections(k)%name_) == trim(adjustl(line(4 : 200)))) recognized = .true. 
!
            enddo
!
            if (.not. recognized) then 
!
               write(output%unit, '(/t3,a,a,a)') 'Could not recognize section named "', trim(adjustl(line(4 : 200))), '".'
!
               call the_file%print_sections()
!
               call output%error_msg('Something is wrong in the input file. See above.')
!
            endif 
!
!           Move to the end of the section again 
!
            do k = 1, n_elements + 1
!
               read(the_file%unit, *) 
!
            enddo
!
         endif 
!
         line = the_file%read_adjustl_lower()
!
      enddo
!
   end subroutine check_for_illegal_sections_input_file
!
!
   subroutine print_sections_input_file(the_file)
!!
!!    Print sections 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
      implicit none 
!
      class(input_file), intent(in) :: the_file
!
      integer :: k
!
      write(output%unit, '(/t3,a/)') 'The valid input sections are:'
!
      do k = 1, size(the_file%sections)
!
         write(output%unit, '(t6,a)') the_file%sections(k)%name_
!
      enddo
!
   end subroutine print_sections_input_file
!
!
   subroutine check_section_for_illegal_keywords_input_file(the_file, the_section)
!!
!!    Checks a section for illegal keywords 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
!!    If there is an illegal keyword in the section, it stops with an error. 
!!    If a legal keyword is specified more than once, it stops with an error.
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
      integer, dimension(:), allocatable :: keywords_instances
!
      allocate(keywords_instances(size(the_section%keywords)))
      keywords_instances = 0
!
      call the_file%move_to_section(the_section%name_, n_records)
!
      do record = 1, n_records
!
         recognized = .false.
!
         read(the_file%unit, '(a200)') string 
!
         if (.not. the_file%string_is_comment(string)) then 
!
            call the_file%extract_keyword_from_string(string, keyword)
!
            do k = 1, size(the_section%keywords)
!
               if (trim(the_section%keywords(k)) == trim(keyword)) then  
!
                  keywords_instances(k) = keywords_instances(k) + 1
                  recognized = .true. 
!
               endif 
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
      do k = 1, size(the_section%keywords)
!
         if (keywords_instances(k) .gt. 1) then 
!
            write(output%unit, '(/t3,a,i0,a,a,a,a,a)') 'Found ', keywords_instances(k), ' instances of the keyword "', &
                              trim(the_section%keywords(k)), '" in the section "', the_section%name_, '".'
!
            call output%error_msg('Something is wrong in the input file. See above.')
!
         endif
!
      enddo
!
      deallocate(keywords_instances)
!
   end subroutine check_section_for_illegal_keywords_input_file
!
!
   logical function requested_reference_calculation_input_file(the_file)
!!
!!    Requested reference calculation 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019 
!!
      implicit none 
!
      class(input_file), intent(in) :: the_file
!
      integer :: n_rf_wfs, k
!
      n_rf_wfs = 0
      do k = 1, size(the_file%rf_wfs)
!
         if (the_file%requested_keyword_in_section(the_file%rf_wfs(k), 'method')) then 
!
            n_rf_wfs = n_rf_wfs + 1
!
         endif 
!
      enddo 
!
      if (n_rf_wfs == 1) then 
!
         requested_reference_calculation_input_file = .true.
!
      elseif (n_rf_wfs > 1) then
!
         requested_reference_calculation_input_file = .false.
         call output%error_msg('Requested more than one reference wavefunction.')
!
      else
!
         requested_reference_calculation_input_file = .false.
!
      endif  
!
   end function requested_reference_calculation_input_file
!
!
   character(len=21) function get_reference_wf_input_file(the_file)
!!
!!    Get reference wavefunction 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019 
!!
      implicit none 
!
      class(input_file), intent(in) :: the_file
!
      integer :: k
!
      logical :: recognized 
!
      recognized = .false.   
!
      do k = 1, size(the_file%rf_wfs)
!
         if (the_file%requested_keyword_in_section(the_file%rf_wfs(k),'method')) then 
!
            get_reference_wf_input_file = the_file%rf_wfs(k)
            recognized = .true. 
!
         endif 
!
      enddo
!
      if (.not. recognized) call output%error_msg('Tried to read reference wavefunction, but could not find any.')
!
   end function get_reference_wf_input_file
!
!
   logical function requested_cc_calculation_input_file(the_file)
!!
!!    Requested CC calculation 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019 
!!
      implicit none 
!
      class(input_file), intent(in) :: the_file
!
      integer :: n_cc_wfs, k
!
      n_cc_wfs = 0
      do k = 1, size(the_file%cc_wfs)
!
         if (the_file%requested_keyword_in_section(the_file%cc_wfs(k), 'method')) then 
!
            n_cc_wfs = n_cc_wfs + 1
!
         endif 
!
      enddo 
!
      if (n_cc_wfs == 1) then 
!
         requested_cc_calculation_input_file = .true.
!
      elseif (n_cc_wfs > 1) then
!
         requested_cc_calculation_input_file = .false.
         call output%error_msg('Requested more than one CC wavefunction.')
!
      else
!
         requested_cc_calculation_input_file = .false.
!
      endif  
!
   end function requested_cc_calculation_input_file
!
!
   character(len=21) function get_cc_wf_input_file(the_file)
!!
!!    Get CC wavefunction 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019 
!!
      implicit none 
!
      class(input_file), intent(in) :: the_file
!
      integer :: k
!
      logical :: recognized 
!
      recognized = .false.   
!
      do k = 1, size(the_file%cc_wfs)
!
         if (the_file%requested_keyword_in_section(the_file%cc_wfs(k),'method')) then 
!
            get_cc_wf_input_file = the_file%cc_wfs(k)
            recognized = .true. 
!
         endif 
!
      enddo
!
      if (.not. recognized) call output%error_msg('Tried to read CC wavefunction, but could not find any.')
!
   end function get_cc_wf_input_file
!
!
   subroutine get_integer_keyword_in_section_input_file(the_file, keyword, section, keyword_value)
!!
!!    Read integer keyword in section 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
!!    If specified, reads keyword as an integer into keyword value.
!!  
      implicit none 
!
      class(input_file), intent(in) :: the_file
!
      character(len=*), intent(in) :: keyword 
      character(len=*), intent(in) :: section  
!
      integer(i6), intent(out) :: keyword_value 
!
      character(len=200) :: keyword_value_string
!
      if (the_file%requested_keyword_in_section(keyword, section)) then 
!
!        Get the keyword value in string format 
!
         call the_file%get_string_keyword_in_section_wo_safety(keyword, section, keyword_value_string)
!
!        Extract the integer from the string
!
         read(keyword_value_string, *) keyword_value
!
      endif
!
   end subroutine get_integer_keyword_in_section_input_file
!
!
   subroutine get_integer8_keyword_in_section_input_file(the_file, keyword, section, keyword_value)
!!
!!    Read integer keyword in section 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
!!    If specified, reads keyword as an integer into keyword value.
!!  
      implicit none 
!
      class(input_file), intent(in) :: the_file
!
      character(len=*), intent(in) :: keyword 
      character(len=*), intent(in) :: section  
!
      integer(i15), intent(out) :: keyword_value 
!
      character(len=200) :: keyword_value_string
!
      if (the_file%requested_keyword_in_section(keyword, section)) then 
!
!        Get the keyword value in string format 
!
         call the_file%get_string_keyword_in_section_wo_safety(keyword, section, keyword_value_string)
!
!        Extract the integer from the string
!
         read(keyword_value_string, *) keyword_value
!
      endif
!
   end subroutine get_integer8_keyword_in_section_input_file
!
!
   subroutine get_required_integer_keyword_in_section_input_file(the_file, keyword, section, keyword_value)
!!
!!    Read required integer keyword in section 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
!!    Reads keyword as an integer into keyword value. If the keyword or the 
!!    section is not specified, an error occurs because the keyword is "required".
!!
!!  
      implicit none 
!
      class(input_file), intent(in) :: the_file
!
      character(len=*), intent(in) :: keyword 
      character(len=*), intent(in) :: section  
!
      integer(i6), intent(out) :: keyword_value 
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
   subroutine get_required_integer8_keyword_in_section_input_file(the_file, keyword, section, keyword_value)
!!
!!    Read required integer keyword in section 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
!!    Reads keyword as an integer into keyword value. If the keyword or the 
!!    section is not specified, an error occurs because the keyword is "required".
!!
!!  
      implicit none 
!
      class(input_file), intent(in) :: the_file
!
      character(len=*), intent(in) :: keyword 
      character(len=*), intent(in) :: section  
!
      integer(i15), intent(out) :: keyword_value 
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
   end subroutine get_required_integer8_keyword_in_section_input_file
!
!
   subroutine get_dp_keyword_in_section_input_file(the_file, keyword, section, keyword_value)
!!
!!    Read double precision keyword in section 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
!!    If specified, reads keyword as a real double precision into keyword value.
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
      if (the_file%requested_keyword_in_section(keyword, section)) then 
!
!        Get the keyword value in string format 
!
         call the_file%get_string_keyword_in_section_wo_safety(keyword, section, keyword_value_string)
!
!        Extract the integer from the string
!
         read(keyword_value_string, *) keyword_value
!
      endif 
!
   end subroutine get_dp_keyword_in_section_input_file
!
!
   subroutine get_required_dp_keyword_in_section_input_file(the_file, keyword, section, keyword_value)
!!
!!    Read required double precision keyword in section 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
!!    Reads keyword as an double precision into keyword value. If the keyword or the 
!!    section is not specified, an error occurs because the keyword is "required".
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
      implicit none 
!
      class(input_file), intent(in) :: the_file
!
      character(len=*), intent(in) :: keyword 
      character(len=*), intent(in) :: section  
!
      character(len=200) :: keyword_value 
!
      if (the_file%requested_keyword_in_section(keyword, section)) then 
!
!        Get the keyword value in string format 
!
         call the_file%get_string_keyword_in_section_wo_safety(keyword, section, keyword_value)
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
      integer :: n_records, record
!
      character(len=200) :: line, local_keyword
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
         if (the_file%string_is_comment(line)) then
!
            cycle
!
         else  
!
            call the_file%extract_keyword_from_string(line, local_keyword)
!
            if (trim(local_keyword) == keyword) then 
!
               call the_file%extract_keyword_value_from_string(line, keyword_value)
               return
!
            endif 
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
      if (tmp_string(1:1) == '!') then 
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
!!    Note: assumes that the string is not a comment.
!!
      implicit none 
!
      character(len=200), intent(in)   :: string 
      character(len=200), intent(out)  :: keyword 
!
      integer :: k, colon_position
!
      keyword = adjustl(string)
      call convert_to_lowercase(keyword)
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
            exit 
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
   subroutine extract_keyword_value_from_string_input_file(string, keyword_value)
!!
!!    Extract keyword value from string 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!    
!!    Note: assumes that the string is not a comment. This routine is therefore 
!!    called after a call to "is comment" logical routine. 
!!
      implicit none 
!
      character(len=200), intent(in)   :: string 
      character(len=200), intent(out)  :: keyword_value 
!
      integer :: k, colon_position
!
!     If there is a ":", we have to remove the ":" as well as the value(s) that follow 
!
      keyword_value = adjustl(string)
      call convert_to_lowercase(keyword_value)
!
      colon_position = -1 
!
      do k = 1, 200
!
         if (keyword_value(k : k) == ':') then 
!
            colon_position = k 
            exit 
!
         endif 
!
      enddo
!
      if (colon_position .ne. -1) then 
!
         do k = 1, colon_position
!
            keyword_value(k : k) = ' '
!
         enddo
!
      endif 
!
      keyword_value = adjustl(keyword_value)
      call convert_to_lowercase(keyword_value)
!
   end subroutine extract_keyword_value_from_string_input_file
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
!!
      implicit none 
!
      class(input_file), intent(in) :: the_file
!
      character(len=*), intent(in) :: keyword 
      character(len=*), intent(in) :: section  
!
      integer :: n_records, record
!
      character(len=200) :: line, local_keyword
!
!     Move to the requested section & get the number of records in that section 
!
      if (.not. the_file%requested_section(section)) then
!
         requested_keyword_in_section_input_file = .false.
         return
!
      endif
!
      call the_file%move_to_section(section, n_records)
!
!     Loop through records within the section & try to locate the keyword 
!
      do record = 1, n_records
!
         read(the_file%unit, '(a200)') line 
!
         if (the_file%string_is_comment(line)) then
!
            cycle
!
         else  
!
            call the_file%extract_keyword_from_string(line, local_keyword)
!
            if (trim(local_keyword) == keyword) then 
!
               requested_keyword_in_section_input_file = .true.
               return
!
            endif 
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
!!    Requested section?
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
!!    Returns true if the section exists, false if it doesn't.
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
        line = the_file%read_adjustl_lower()
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
      integer, intent(out) :: n_records
!
      character(len=200) :: line 
!
      integer :: start_record, end_record 
!
      integer :: n_beginnings, n_ends 
!
!     Find the number of instances of the section,
!     stopping with an error if there are any inconsistencies 
!
      rewind(the_file%unit)
!
      line = repeat(' ', 200)
!
      n_ends = 0
      n_beginnings = 0
!
      do while (trim(line) /= 'end geometry')
!
         line = the_file%read_adjustl_lower()
!
         if (trim(line) == 'end ' // string) n_ends = n_ends + 1
         if (trim(line) == string) n_beginnings = n_beginnings + 1      
!
      enddo    
!
      if (n_beginnings > 1) call output%error_msg('Tried to move to section "' // string // '" with more than one starting clause.')
      if (n_ends > 1) call output%error_msg('Tried to move to section "' // string // '" with more than one ending clause.')
      if (n_ends == 0 .and. n_beginnings == 0) call output%error_msg('Tried to move to non-existent section "' // string // '".')
      if (n_ends < 1) call output%error_msg('Tried to move to section "' // string // '" with no end.')
      if (n_beginnings < 1) call output%error_msg('Tried to move to section "' // string // '" with no beginning.')
!
!     Find the end of the section 
!
      rewind(the_file%unit)
!
      line = the_file%read_adjustl_lower()
!
      end_record = 1
!
      do while (trim(line) /= 'end geometry' .and. trim(line) /= 'end ' // string) 
!
         end_record = end_record + 1
! 
        line = the_file%read_adjustl_lower()     
!
      enddo   
!
!     Find the beginning of the section;
!     this also places the pointer in the correct position
!
      rewind(the_file%unit)
!
      line = the_file%read_adjustl_lower()
!
      start_record = 1
!
      do while (trim(line) /= 'end geometry' .and. trim(line) /= string) 
!
         start_record = start_record + 1
! 
        line = the_file%read_adjustl_lower()       
!
      enddo
!
!     Set the number of records inside the section 
!
      n_records = end_record - start_record - 1
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
   function get_n_atoms_input_file(the_file) result(n_atoms)
!!
!!    Get n atoms 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
!!    Reads the geometry section of the input and
!!    counts the number of atoms
!!
      implicit none
!
      class(input_file), intent(in) :: the_file
!
      integer :: n_atoms
!
!     Local variables
!
      integer :: n_records, record
!
      character(len=200) :: string 
!
      n_atoms = 0
!
      call the_file%move_to_section('geometry', n_records)
!
      do record = 1, n_records
!
         string = the_file%read_adjustl_lower()
!
         if(string(1:6) /= 'basis:') n_atoms = n_atoms + 1
!
      enddo
!
   end function  get_n_atoms_input_file
!
!
   subroutine get_geometry_input_file(the_file, n_atoms, symbols, positions, basis_sets)
!!
!!    Get geometry
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
!!    Reads the geometry from the output file and sets it in the
!!    list of atoms.
!!
!!    Note: This routine should be called after a call to the 
!!    function get_n_atoms_input_file which finds the number of 
!!    atoms in the system
!!
      implicit none
!
      class(input_file), intent(in) :: the_file
!
      integer, intent(in) :: n_atoms
!
      character(len=2), dimension(n_atoms), intent(out)   :: symbols
      character(len=100), dimension(n_atoms), intent(out) :: basis_sets
!
      real(dp), dimension(n_atoms, 3), intent(out) :: positions ! x, y, z
!
!     Local variables
!
      integer :: n_records, record, cursor, current_atom, i
!
      character(len=200) :: string, coordinate
      character(len=100) :: current_basis
!
      call the_file%move_to_section('geometry', n_records)
!
!     Set initial basis -> Error if not
!
      string = the_file%read_adjustl_lower()
!
      if(string(1:6) /= 'basis:') call output%error_msg('did not find basis in geometry section.')
      current_basis = trim(adjustl(string(7:200)))
!
!     Loop through the rest of the geometry section to get atoms
!
      current_atom = 0
!
      do record = 1, n_records - 1
!
         string = the_file%read_adjustl_lower()
!
         if(string(1:6) == 'basis:') then
!
            current_basis = trim(adjustl(string(7:200)))
!
         else
!
            current_atom = current_atom + 1
!
            basis_sets(current_atom) = current_basis
            symbols(current_atom)    = trim(string(1:2))
!
            string = string(3:200)
!
            cursor = set_cursor_to_whitespace(string)
!
            coordinate = string(1:cursor)
            read(coordinate, '(f21.16)') positions(current_atom, 1)
!
            string = string(cursor + 1:200)
!
            cursor = set_cursor_to_whitespace(string)
!
            coordinate = string(1:cursor)
            read(coordinate, '(f21.16)') positions(current_atom, 2)
!
            coordinate = string(cursor + 1:200)
            coordinate = adjustl(coordinate)
!
            read(coordinate, '(f21.16)') positions(current_atom, 3)
!
         endif
!
      enddo
!
!     First character of symbol should be upper case
!
      do i = 1, n_atoms
!
         symbols(i)(1:1) = convert_char_to_uppercase(symbols(i)(1:1))
!
      enddo
!
   end subroutine get_geometry_input_file
!
!
   function read_adjustl_lower_input_file(the_file) result(line)
!!
!!    Read, adjustl and convert to lower case
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
      implicit none
!
      class(input_file), intent(in) :: the_file
!
      character(len=200) :: line
!  
      read(the_file%unit, '(a200)') line
      line = adjustl(line)
      call convert_to_lowercase(line)
!
   end function read_adjustl_lower_input_file

!
!
end module input_file_class
