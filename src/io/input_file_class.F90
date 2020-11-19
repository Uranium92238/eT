!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2020 the authors of eT
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
   use abstract_file_class, only: abstract_file
   use global_out, only: output
!
   type, extends(abstract_file) :: input_file
!
      type(section), allocatable :: sections(:)
!
      character(len=30), allocatable :: rf_wfs(:)
      character(len=30), allocatable :: cc_wfs(:)
!
      integer :: n_keyword_lines ! Number of lines excluding the geometry section
      integer :: n_geometry_lines    ! Number of QM atoms
      integer :: n_mm_atom_lines ! Number of MM atoms
!
      character(len=200), dimension(:), allocatable :: input_  ! Array of the input
                                                               ! lines excluding geometry
!
      character(len=200), dimension(:), allocatable :: geometry    ! Array of the QM input geometry
      character(len=200), dimension(:), allocatable :: mm_geometry ! Array of the MM input geometry
!
   contains
!
      procedure :: open_                                                => open_input_file
      procedure :: close_                                               => close_input_file
!
      procedure :: check_for_errors                                     => check_for_errors_input_file
!
      procedure :: requested_section                                    => requested_section_input_file
      procedure :: requested_keyword_in_section                         => requested_keyword_in_section_input_file
      procedure :: get_n_elements_for_keyword_in_section                => get_n_elements_for_keyword_in_section_input_file
      procedure :: get_integer_array_for_keyword_in_section             => get_integer_array_for_keyword_in_section_input_file
      procedure :: get_real_array_for_keyword_in_section                => get_real_array_for_keyword_in_section_input_file
      procedure :: get_n_atoms                                          => get_n_atoms_input_file
      procedure :: get_mm_n_atoms                                       => get_mm_n_atoms_input_file
      procedure :: get_geometry                                         => get_geometry_input_file
      procedure :: get_mm_geometry                                      => get_mm_geometry_input_file
      procedure :: get_reference_wf                                     => get_reference_wf_input_file
      procedure :: get_cc_wf                                            => get_cc_wf_input_file
!
      procedure :: requested_reference_calculation                      => requested_reference_calculation_input_file
      procedure :: requested_cc_calculation                             => requested_cc_calculation_input_file
!
      procedure, private :: get_string_keyword_in_section_wo_safety     => get_string_keyword_in_section_wo_safety_input_file
      procedure, private :: get_section_limits                          => get_section_limits_input_file
      procedure, private :: check_section_for_illegal_keywords          => check_section_for_illegal_keywords_input_file
      procedure, private :: check_for_illegal_sections                  => check_for_illegal_sections_input_file
      procedure, private :: print_sections                              => print_sections_input_file
!
      procedure, nopass, private :: string_is_comment                   => string_is_comment_input_file
      procedure, nopass, private :: extract_keyword_from_string         => extract_keyword_from_string_input_file
      procedure, nopass, private :: extract_keyword_value_from_string   => extract_keyword_value_from_string_input_file
!
      generic :: get_keyword_in_section                                 => get_integer4_keyword_in_section_input_file,  &
                                                                           get_integer8_keyword_in_section_input_file,  &
                                                                           get_string_keyword_in_section_input_file,    &
                                                                           get_dp_keyword_in_section_input_file
!
      generic :: get_required_keyword_in_section                        => get_required_string_keyword_in_section_input_file,    &
                                                                           get_required_integer4_keyword_in_section_input_file,  &
                                                                           get_required_integer8_keyword_in_section_input_file,  &
                                                                           get_required_dp_keyword_in_section_input_file
!
      generic :: get_array_for_keyword_in_section                       => get_integer_array_for_keyword_in_section, &
                                                                           get_real_array_for_keyword_in_section
!
      procedure :: get_integer4_keyword_in_section_input_file
      procedure :: get_integer8_keyword_in_section_input_file
      procedure :: get_string_keyword_in_section_input_file
      procedure :: get_dp_keyword_in_section_input_file
      procedure :: get_required_string_keyword_in_section_input_file
      procedure :: get_required_integer4_keyword_in_section_input_file
      procedure :: get_required_integer8_keyword_in_section_input_file
      procedure :: get_required_dp_keyword_in_section_input_file
!
      procedure :: is_string_in_cs_list   => is_string_in_cs_list_input_file
!
      procedure :: read_keywords_and_geometry => read_keywords_and_geometry_input_file
!
      procedure :: cleanup_geometry    => cleanup_geometry_input_file
      procedure :: cleanup_keywords    => cleanup_keywords_input_file
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
   function new_input_file(name_) result(the_file)
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
      type(section) :: calculations
      type(section) :: system
      type(section) :: memory
      type(section) :: cc_mean_value
      type(section) :: cc_response
      type(section) :: hf_mean_value
      type(section) :: cc_td
      type(section) :: method
      type(section) :: solver_cholesky
      type(section) :: solver_scf
      type(section) :: solver_scf_geoopt
      type(section) :: solver_cc_gs
      type(section) :: solver_cc_es
      type(section) :: solver_cc_multipliers
      type(section) :: solver_cc_response
      type(section) :: solver_cc_propagation
      type(section) :: solver_fft_dipole_moment
      type(section) :: solver_fft_electric_field
      type(section) :: electric_field
      type(section) :: active_atoms
      type(section) :: cc
      type(section) :: mlcc
      type(section) :: mm
      type(section) :: mlhf
      type(section) :: visualization
      type(section) :: pcm
      type(section) :: global_print
      type(section) :: frozen_orbitals
      type(section) :: integrals
!
!     Set input file name, access and format
!
      the_file%name_ = name_
!
      the_file%access_ = 'sequential'
      the_file%format_ = 'formatted'
      the_file%action_ = 'read'
!
!     Set method section
!
      method%name_    = 'method'
      method%required = .false.
!
      the_file%rf_wfs = [character(len=30) :: &
                        'hf',  &
                        'uhf', &
                        'mlhf']
!
      the_file%cc_wfs = [character(len=30) :: &
                        'ccs',                &
                        'mp2',                &
                        'cc2',                &
                        'lowmem-cc2',         &
                        'ccsd',               &
                        'cc3',                &
                        'ccsd(t)',            &
                        'mlcc2',              &
                        'mlccsd']
!
      allocate(method%keywords(size(the_file%rf_wfs) &
                             + size(the_file%cc_wfs)))
!
      method%keywords = [the_file%rf_wfs, the_file%cc_wfs]
!
!     Set other sections in alphabetical order
!
      active_atoms%name_    = 'active atoms'
      active_atoms%required = .false.
      active_atoms%keywords = [character(len=30) :: &
                              'selection type',     &
                              'central atom',       &
                              'hf',                 &
                              'ccs',                &
                              'cc2',                &
                              'ccsd',               &
                              'cc3',                &
                              'ccsd(t)',            &
                              'inactive basis',     &
                              'hf basis',           &
                              'ccs basis',          &
                              'cc2 basis',          &
                              'ccsd basis',         &
                              'cc3 basis',          &
                              'ccsd(t) basis ']
!
!
      cc%name_    = 'cc'
      cc%required = .false.
      cc%keywords = [character(len=30) :: 'bath orbital']
!
!
      cc_mean_value%name_    = 'cc mean value'
      cc_mean_value%required = .false.
      cc_mean_value%keywords = [character(len=30) :: &
                               'dipole',             &
                               'quadrupole']
!
!
      cc_response%name_    = 'cc response'
      cc_response%required = .false.
      cc_response%keywords = [character(len=30) :: &
                             'transition moments', &
                             'frequencies',        &
                             'polarizabilities',   &
                             'dipole length',      &
                             'lr',                 &
                             'eom']
!
!
      cc_td%name_    = 'cc td'
      cc_td%required = .false.
      cc_td%keywords = [character(len=30) :: &
                       'propagation',        &
                       'fft dipole moment',  &
                       'fft electric field']
!
!
      calculations%name_    = 'do'
      calculations%required = .true.
      calculations%keywords = [character(len=30) ::   &
                              'ground state',         &
                              'ground state geoopt',  &
                              'excited state',        &
                              'response',             &
                              'mean value',           &
                              'time dependent state', &
                              'cholesky eri',         &
                              'restart']
!
!
      electric_field%name_    = 'electric field'
      electric_field%required = .false.
      electric_field%keywords = [character(len=30) ::        &
                                'envelope',                  &
                                'x polarization',            &
                                'y polarization',            &
                                'z polarization',            &
                                'central time',              &
                                'width',                     &
                                'central angular frequency', &
                                'peak strength',             &
                                'phase shift',               &
                                'repetition',                &
                                'separation']
!
!
      hf_mean_value%name_    = 'hf mean value'
      hf_mean_value%required = .false.
      hf_mean_value%keywords = [character(len=30) :: &
                               'dipole',             &
                               'quadrupole']
!
!
      frozen_orbitals%name_    = 'frozen orbitals'
      frozen_orbitals%required = .false.
      frozen_orbitals%keywords = [character(len=30) :: &
                                 'hf', &
                                 'core']
!
!
      integrals%name_    = 'integrals'
      integrals%required = .false.
      integrals%keywords = [character(len=30) :: &
                           'cholesky storage',   &
                           'mo eri in memory',   &
                           't1 eri in memory',   &
                           'eri storage']
!
!
      memory%name_    = 'memory'
      memory%required = .false.
      memory%keywords = [character(len=30) :: &
                        'available', &
                        'unit']
!
      mlcc%name_    = 'mlcc'
      mlcc%required = .false.
      mlcc%keywords = [character(len=30) ::     &
                      'levels',                 &
                      'cc2 orbitals',           &
                      'ccsd orbitals',          &
                      'cholesky threshold',     &
                      'cnto restart',           &
                      'orbital restart',        &
                      'cnto occupied cc2',      &
                      'cnto virtual cc2',       &
                      'cnto occupied ccsd',     &
                      'cnto virtual ccsd',      &
                      'cnto states',            &
                      'nto states',             &
                      'nto occupied cc2',       &
                      'nto occupied ccsd',      &
                      'canonical virtual cc2',  &
                      'canonical virtual ccsd', &
                      'print ccs calculation',  &
                      'print cc2 calculation']
!
!
      mm%name_    = 'molecular mechanics'
      mm%required = .false.
      mm%keywords = [character(len=30) :: &
                    'forcefield',         &
                    'algorithm ']
!
!
      mlhf%name_    = 'multilevel hf'
      mlhf%required = .false.
      mlhf%keywords = [character(len=30) ::           &
                      'initial hf optimization',    &
                      'initial hf threshold',       &
                      'print initial hf',           &
                      'cholesky threshold',         &
                      'project on minimal basis',   &
                      'cholesky virtuals']
!
!
      pcm%name_    = 'pcm'
      pcm%required = .false.
      pcm%keywords = [character(len=30) :: &
                     'solvent',            &
                     'input',              &
                     'tesserae area',      &
                     'solver type']
!
!
      global_print%name_    = 'print'
      global_print%required = .false.
      global_print%keywords = [character(len=30) ::  &
                              'output print level ', &
                              'timing print level ']
!
!
      solver_cholesky%name_    = 'solver cholesky'
      solver_cholesky%required = .false.
      solver_cholesky%keywords = [character(len=30) :: &
                                 'threshold',          &
                                 'span',               &
                                 'batches',            &
                                 'qualified',          &
                                 'one center',         &
                                 'no vectors',         &
                                 'mo screening']
!
!
      solver_scf%name_    = 'solver scf'
      solver_scf%required = .false.
      solver_scf%keywords = [character(len=30) ::        &
                            'algorithm',                 &
                            'energy threshold',          &
                            'gradient threshold',        &
                            'storage',                   &
                            'crop',                      &
                            'cumulative fock threshold', &
                            'max iterations',            &
                            'coulomb threshold',         &
                            'exchange threshold',        &
                            'integral precision',        &
                            'integral cutoff',           &
                            'diis dimension',            &
                            'restart',                   &
                            'ao density guess',          &
                            'print orbitals']
!
!
      solver_scf_geoopt%name_    = 'solver scf geoopt'
      solver_scf_geoopt%required = .false.
      solver_scf_geoopt%keywords = [character(len=30) :: &
                                   'algorithm',          &
                                   'max step',           &
                                   'energy threshold',   &
                                   'gradient threshold', &
                                   'max iterations',     &
                                   'restart']
!
!
      solver_cc_es%name_    = 'solver cc es'
      solver_cc_es%required = .false.
      solver_cc_es%keywords = [character(len=30) ::       &
                              'algorithm',                &
                              'core excitation',          &
                              'ionization',               &
                              'energy threshold ',        &
                              'crop',                     &
                              'residual threshold',       &
                              'max iterations',           &
                              'restart',                  &
                              'left eigenvectors',        &
                              'right eigenvectors',       &
                              'storage',                  &
                              'singlet states',           &
                              'diis dimension',           &
                              'max micro iterations',     &
                              'max reduced dimension',    &
                              'davidson preconvergence',  &
                              'preconvergence threshold', &
                              'max micro iterations',     &
                              'rel micro threshold',      &
                              'chain length',             &
                              'lanczos normalization']
!
!
      solver_cc_gs%name_    = 'solver cc gs'
      solver_cc_gs%required = .false.
      solver_cc_gs%keywords = [character(len=30) ::   &
                              'algorithm',            &
                              'energy threshold',     &
                              'omega threshold',      &
                              'crop',                 &
                              'max micro iterations', &
                              'rel micro threshold',  &
                              'storage',              &
                              'max iterations',       &
                              'diis dimension',       &
                              'restart' ]
!
!
      solver_cc_multipliers%name_    = 'solver cc multipliers'
      solver_cc_multipliers%required = .false.
      solver_cc_multipliers%keywords = [character(len=30) ::    &
                                       'algorithm',             &
                                       'threshold',             &
                                       'storage',               &
                                       'crop',                  &
                                       'diis dimension',        &
                                       'restart',               &
                                       'max reduced dimension', &
                                       'max iterations']
!
!
      solver_cc_propagation%name_    = 'solver cc propagation'
      solver_cc_propagation%required = .false.
      solver_cc_propagation%keywords = [character(len=30) ::    &
                                       'initial time',          &
                                       'final time',            &
                                       'time step',             &
                                       'steps between output',  &
                                       'implicit threshold',    &
                                       'energy output',         &
                                       'dipole moment output',  &
                                       'electric field output', &
                                       'amplitudes output',     &
                                       'multipliers output',    &
                                       'density matrix output', &
                                       'integrator']
!
!
      solver_cc_response%name_    = 'solver cc response'
      solver_cc_response%required = .false.
      solver_cc_response%keywords = [character(len=30) ::    &
                                    'threshold            ', &
                                    'storage              ', &
                                    'max iterations       ']
!
!
      solver_fft_dipole_moment%name_    = 'solver fft dipole moment'
      solver_fft_dipole_moment%required = .false.
      solver_fft_dipole_moment%keywords = [character(len=30) :: & 
                                          'initial time',       &
                                          'final time',         &
                                          'time step']
!
!
      solver_fft_electric_field%name_    = 'solver fft electric field'
      solver_fft_electric_field%required = .false.
      solver_fft_electric_field%keywords = [character(len=30) :: &
                                           'initial time',       &
                                           'final time',         &
                                           'time step']
!
!
      system%name_    = 'system'
      system%required = .true.
      system%keywords = [character(len=30) ::'name', &
                        'cartesian gaussians',       &
                        'pure gaussians',            &
                        'charge',                    &
                        'multiplicity']
!
!
      visualization%name_    = 'visualization'
      visualization%required = .false.
      visualization%keywords = [character(len=30) ::        &
                               'grid spacing',              &
                               'grid buffer',               &
                               'plot cc density',           &
                               'plot hf orbitals',          &
                               'plot hf density',           &
                               'plot hf active density',    &
                               'plot transition densities', &
                               'states to plot']
!
!     Gather all sections into the file's section array
!
      the_file%sections = [active_atoms,              &
                           calculations,              &
                           cc,                        &
                           cc_mean_value,             &
                           cc_response,               &
                           cc_td,                     &
                           electric_field,            &
                           frozen_orbitals,           &
                           global_print,              &
                           hf_mean_value,             &
                           integrals,                 &
                           mlcc,                      &
                           mlhf,                      &
                           mm,                        &
                           memory,                    &
                           method,                    &
                           pcm,                       &
                           solver_cholesky,           &
                           solver_cc_gs,              &
                           solver_cc_es,              &
                           solver_cc_multipliers,     &
                           solver_cc_response,        &
                           solver_cc_propagation,     &
                           solver_fft_dipole_moment,  &
                           solver_fft_electric_field, &
                           solver_scf,                &
                           solver_scf_geoopt,         &
                           system,                    &
                           visualization]
!
      the_file%is_open = .false.
      the_file%unit_ = -1
!
   end function new_input_file
!
!
   subroutine open_input_file(the_file)
!!
!!    Open the input file
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(input_file) :: the_file
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      if (the_file%is_open) then
!
         call output%error_msg(trim(the_file%name_)//' is already open')
!
      endif
!
      open(newunit=the_file%unit_, file=the_file%name_, access=the_file%access_, &
           action='read', status='unknown', form=the_file%format_, iostat=io_error, iomsg=io_msg)
!
      if (io_error .ne. 0) then
!
         call output%error_msg('could not open eT input file '//trim(the_file%name_)//&
                              &'error message: '//trim(io_msg))
      endif
!
      the_file%is_open = .true.
!
   end subroutine open_input_file
!
!
   subroutine close_input_file(the_file)
!!
!!    Close the input file
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(input_file) :: the_file
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      if (.not. the_file%is_open) then
         call output%error_msg(trim(the_file%name_)//' already closed')
      end if
!
      close(the_file%unit_, iostat=io_error, iomsg=io_msg, status='keep')
!
      if (io_error .ne. 0) then
!
         call output%error_msg('could not close eT input file '//trim(the_file%name_)//&
                              &'error message: '//trim(io_msg))
      endif
!
      the_file%is_open = .false.
      the_file%unit_ = -1
!
   end subroutine close_input_file
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
               call output%printf('m', 'All calculations require the section "' // &
                                  trim(the_file%sections(k)%name_) // '". It &
                                  &appears to be missing.', fs='(/t3,a)')
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
      integer :: k, start_, end_, i
!
      logical :: recognized
!
      do i = 1, the_file%n_keyword_lines
!
         if (the_file%input_(i)(1 : 3) == 'end') then
!
!           Located the end of a section -
!           attempt to move to the beginning of that section => fails if inconsistent beginning and end
!
            call the_file%get_section_limits(trim(adjustl(the_file%input_(i)(4 : 200))), start_, end_)
!
!           Check whether section name is valid
!
            recognized = .false.
!
            do k = 1, size(the_file%sections)
!
               if (trim(the_file%sections(k)%name_) == trim(adjustl(the_file%input_(i)(4 : 200)))) &
                  recognized = .true.
!
            enddo
!
            if (.not. recognized) then
!
               call output%printf('m', 'Could not recognize section named "' //  &
                                   trim(adjustl(the_file%input_(i)(4 : 200))) // '".', fs='(/t3,a)')
!
               call the_file%print_sections()
!
               call output%error_msg('Something is wrong in the input file. See above.')
!
            endif
!
         endif
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
      call output%printf('m', 'The valid input sections are:', fs='(/t3,a/)')
!
      do k = 1, size(the_file%sections)
!
         call output%printf('m', '(a0)', chars=[trim(the_file%sections(k)%name_)], fs='(t6,a)')
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
      character(len=200) :: keyword
!
      logical :: recognized
!
      integer :: record, k, start_, end_
!
      integer, dimension(:), allocatable :: keywords_instances
!
      allocate(keywords_instances(size(the_section%keywords)))
      keywords_instances = 0
!
      call the_file%get_section_limits(the_section%name_, start_, end_)
!
      do record = start_ + 1, end_ - 1
!
         recognized = .false.
!
         if (.not. the_file%string_is_comment(the_file%input_(record))) then
!
            call the_file%extract_keyword_from_string(the_file%input_(record), keyword)
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
               call output%printf('m', 'Could not recognize keyword "' //  &
                                  trim(keyword) // '" in section "' //  &
                                  trim(the_section%name_) // '".', fs='(/t3,a)')

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
         call output%printf('m', 'Found (i0) instances of the keyword "' //  &
                            trim(the_section%keywords(k)) // '" in the section "' &
                            // the_section%name_ // '".', &
                            ints=[keywords_instances(k)], fs='(/t3,a)')
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
   character(len=30) function get_reference_wf_input_file(the_file)
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
   character(len=30) function get_cc_wf_input_file(the_file)
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
!
!
   subroutine get_integer4_keyword_in_section_input_file(the_file, keyword, section, keyword_value)
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
      integer(i32), intent(inout) :: keyword_value
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
   end subroutine get_integer4_keyword_in_section_input_file
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
      integer(i64), intent(inout) :: keyword_value
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
   subroutine get_required_integer4_keyword_in_section_input_file(the_file, keyword, section, keyword_value)
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
      integer(i32), intent(inout) :: keyword_value
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
   end subroutine get_required_integer4_keyword_in_section_input_file
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
      integer(i64), intent(out) :: keyword_value
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
      real(dp), intent(inout) :: keyword_value
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
      integer :: record, start_, end_
!
      character(len=200) :: local_keyword
!
!     Move to the requested section & get the number of records in that section
!
      call the_file%get_section_limits(section, start_, end_)
!
!     Loop through records within the section to locate & get the keyword value
!
      do record = start_ + 1, end_ - 1
!
         if (the_file%string_is_comment(the_file%input_(record))) then
!
            cycle
!
         else
!
            call the_file%extract_keyword_from_string(the_file%input_(record), local_keyword)
!
            if (trim(local_keyword) == keyword) then
!
               call the_file%extract_keyword_value_from_string(the_file%input_(record), &
                  keyword_value)
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
      integer :: record, start_, end_
!
      character(len=200) :: local_keyword
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
      call the_file%get_section_limits(section, start_, end_)
!
!     Loop through records within the section & try to locate the keyword
!
      do record = start_ + 1, end_ - 1
!
         if (the_file%string_is_comment(the_file%input_(record))) then
!
            cycle
!
         else
!
            call the_file%extract_keyword_from_string(the_file%input_(record), local_keyword)
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
   pure function requested_section_input_file(the_file, section) result(requested)
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
      logical :: requested
!
      integer :: i
!
      requested = .false.
!
      do i = 1, the_file%n_keyword_lines
!
         if (trim(the_file%input_(i)) == section) then
!
            requested = .true.
            return
!
         endif
!
      enddo
!
   end function requested_section_input_file
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
      n_elements = 0
!
      if (.not. the_file%requested_keyword_in_section(keyword, section)) return
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
   subroutine get_integer_array_for_keyword_in_section_input_file(the_file, keyword, section, n_elements, array_)
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
      if (.not. the_file%requested_keyword_in_section(keyword, section)) return
!
!     Get the keyword value in string format
!
      call the_file%get_keyword_in_section(keyword, section, keyword_value_string)
!
!     Use string utility functionality to get the array
!
      call get_elements_in_string(keyword_value_string, n_elements, array_)
!
   end subroutine get_integer_array_for_keyword_in_section_input_file
!
!
   subroutine get_real_array_for_keyword_in_section_input_file(the_file, keyword, section, n_elements, array_)
!!
!!    Get real array for keyword in section
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!    Modified by Andreas Skeidsvoll, Sep 2019: Reads real array instead of integer array
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
      real(dp), dimension(n_elements) :: array_
!
      character(len=200) :: keyword_value_string
!
!     Get the keyword value in string format
!
      call the_file%get_keyword_in_section(keyword, section, keyword_value_string)
!
!     Use string utility functionality to get the real array
!
      call get_reals_in_string(keyword_value_string, n_elements, array_)
!
   end subroutine get_real_array_for_keyword_in_section_input_file
!
!
   pure function get_n_atoms_input_file(the_file) result(n_atoms)
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
      integer :: record
!
      n_atoms = 0
!
      do record = 1, the_file%n_geometry_lines
!
         if (the_file%geometry(record)(1:6) .ne. 'basis:' &
            .and. (the_file%geometry(record)(1:6) .ne. 'units:')) &
            n_atoms = n_atoms + 1
!
      enddo
!
   end function  get_n_atoms_input_file
!
!
   pure function get_mm_n_atoms_input_file(the_file) result(n_atoms)
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
      n_atoms = the_file%n_mm_atom_lines
!
   end function  get_mm_n_atoms_input_file
!
!
   subroutine get_geometry_input_file(the_file, n_atoms, symbols, &
                                       positions, basis_sets, units_angstrom)
!!
!!    Get geometry
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!    Modified by Åsmund H. Tveten, Oct 2019. Generalized to Bohr units.
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
      logical, intent(out) :: units_angstrom ! True if units are Ångström/unspecified, false if Bohr
!
!     Local variables
!
      integer :: cursor, current_atom, i
!
      integer :: start_
!
      character(len=200) :: string, coordinate
      character(len=100) :: current_basis
!
      start_ = 1 ! Specifies the line of the first and required basis
!
!     Are units specified?
!     Note that units can only be specified as the first line of the geometry
!
      units_angstrom = .true. ! Default units are Angstrom
!
      if (the_file%geometry(1)(1:6) == 'units:') then
!
         string = (trim(adjustl(the_file%geometry(1)(7:200))))
!
         if (string(1:4) == 'bohr') then
!
            units_angstrom = .false.
!
         elseif (string(1:8) /= 'angstrom') then
!
            call output%error_msg('units of atom coordinates must be either angstrom or bohr')
!
         endif
!
         start_ = 2
!
      endif
!
!     Error if next line is not a basis set line
!
      if(the_file%geometry(start_)(1:6) /= 'basis:') &
            call output%error_msg('did not find basis in geometry section.')
!
!     Loop through geometry
!
      current_atom = 0
!
      do i = start_, the_file%n_geometry_lines
!
         if (the_file%geometry(i)(1:6) == 'units:') &
         call output%error_msg('Units must be specified as the first line in the geometry section.')
!
         if(the_file%geometry(i)(1:6) == 'basis:') then
!
            current_basis = trim(adjustl(the_file%geometry(i)(7:200)))
!
         else
!
            current_atom = current_atom + 1
!
            basis_sets(current_atom) = current_basis
            symbols(current_atom)    = trim(the_file%geometry(i)(1:2))
!
            string = the_file%geometry(i)(3:200)
!
            cursor = set_cursor_to_character(string)
!
            coordinate = string(1:cursor)
            read(coordinate, '(f25.16)') positions(current_atom, 1)
!
            string = string(cursor + 1:200)
!
            cursor = set_cursor_to_character(string)
!
            coordinate = string(1:cursor)
            read(coordinate, '(f25.16)') positions(current_atom, 2)
!
            coordinate = string(cursor + 1:200)
            coordinate = adjustl(coordinate)
!
            read(coordinate, '(f25.16)') positions(current_atom, 3)
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
   subroutine get_mm_geometry_input_file(the_file, n_atoms, symbols, mols, positions, charge, chi, eta)
!!
!!    Get MM geometry
!!    Written by Tommaso Giovanini, May 2019
!!
!!    Reads the geometry of the MM portion from the input file and
!!    sets it in the list of atoms.
!!
!!    Note: In order to be run, you need to know the number of MM atoms
!!          in the system
!!
      implicit none
!
      class(input_file), intent(in) :: the_file
!
      integer, intent(in) :: n_atoms
!
      character(len=2), dimension(n_atoms), intent(out)   :: symbols
!
      real(dp), dimension(3, n_atoms), intent(out) :: positions ! x, y, z
      real(dp), dimension(n_atoms), intent(out), optional :: charge ! charges for non-polarizable QMMM
      real(dp), dimension(n_atoms), intent(out), optional :: chi ! electronegativity for FQ
      real(dp), dimension(n_atoms), intent(out), optional :: eta ! chemical hardness for FQ
!
      integer, dimension(n_atoms), intent(out) :: mols
!
!     Local variables
!
      integer :: record, cursor, current_atom, i
!
      character(len=200) :: string, coordinate, imolecule, charge_read, chi_read, eta_read
!
!     Loop through the rest of the geometry section to get atoms
!
      current_atom = 0
!
      do record = 1, the_file%n_mm_atom_lines
!
         current_atom = current_atom + 1
!
         symbols(current_atom)    = trim(the_file%mm_geometry(record)(1:2))
!
         string = the_file%mm_geometry(record)(3:200)
!
         cursor = set_cursor_to_character(string,'=')
!
         string = string(cursor+1:200)
!
         cursor = set_cursor_to_character(string,']')
!
         imolecule = string(1:cursor-1)
         read(imolecule,'(i4)') mols(current_atom)
!
         string = string(cursor+1:200)
!
         cursor = set_cursor_to_character(string)
!
         coordinate = string(1:cursor)
         read(coordinate, '(f21.16)') positions(1, current_atom)
!
         string = string(cursor + 1:200)
!
         cursor = set_cursor_to_character(string)
!
         coordinate = string(1:cursor)
         read(coordinate, '(f21.16)') positions(2, current_atom)
!
         string = string(cursor + 1:200)
!
         cursor = set_cursor_to_character(string)
!
         coordinate = string(1:cursor)
         read(coordinate, '(f21.16)') positions(3, current_atom)
!
         if (present (charge)) then ! non-polarizable QM/MM [q=value]
!
             string = string(cursor + 1:200)
!
             cursor = set_cursor_to_character(string,'q')
!
             string = string(cursor+2:200)
!
             cursor = set_cursor_to_character(string,']')
!
             charge_read = string(1:cursor-1)
             read(charge_read,'(f21.16)') charge(current_atom)
!
             if(abs(charge(current_atom)).lt.1.0d-8) then
!
                call output%printf('m', 'Electrostatic Embedding QM/MM ', fs='(/t6,a)')
                call output%warning_msg('You put zero charge on atom = (i0)', &
                                         ints=[current_atom], fs='(t6,a)')
!
             endif
!
         else if (present (chi).and.present(eta)) then ! polarizable QM/FQ [chi=value,eta=value]
!
             string = string(cursor + 1:200)
!
             cursor = set_cursor_to_character(string,'=')
!
             string = string(cursor+1:200)
!
             cursor = set_cursor_to_character(string,',')
!
             chi_read = string(1:cursor-1)
             read(chi_read,'(f21.16)') chi(current_atom)
!
             string = string(cursor + 1:200)
!
             cursor = set_cursor_to_character(string,'=')
!
             string = string(cursor+1:200)
!
             cursor = set_cursor_to_character(string,']')
!
             eta_read = string(1:cursor-1)
             read(eta_read,'(f21.16)') eta(current_atom)
!
             if(abs(eta(current_atom)).lt.1.0d-8) then
!
               call output%printf('m', 'Polarizable QM/FQ ', fs='(/t6,a)')
               call output%error_msg('You have put zero chemical hardness on atom: (i0)', &
                                      ints=[current_atom])
!
             endif
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
   end subroutine get_mm_geometry_input_file
!
!
   function is_string_in_cs_list_input_file(the_file, keyword, section, string) result(found)
!!
!!    Is string in comma separated list
!!    Written by Sarai D. Folkestad, 2019-2020
!!
!!    Determines if the 'string' is one of the keyword values in the comma separated
!!    list given with 'keyword' in the given 'section' of the input file.
!!
!!    Examples of such a keyword is 'levels' from the 'mlcc' section:
!!
!!       levels: ccs, cc2, ccsd
!!
!!    This function assumes that the keyword is requested. If it is not, a test should be performed
!!    before the routine is called
!!
      implicit none
!
      class(input_file), intent(in) :: the_file
!
      character(len=*), intent(in) :: section
      character(len=*), intent(in) :: keyword
      character(len=*), intent(in) :: string
!
      logical :: found
!
      character(len=200) :: line
      character(len=200) :: current_snippet
!
      integer :: i, current_cursor
!
      found = .false.
!
      call the_file%get_required_keyword_in_section(keyword, section, line)
!
      current_cursor = 1
      current_snippet = repeat(' ', 200)
!
      do i = 1, len_trim(line)
!
         if (line(i:i) .ne. ',') then
!
            current_snippet(current_cursor:current_cursor) = line(i:i)
!
            current_cursor = current_cursor + 1
!
            if (i == len_trim(line)) then ! the end of the line
!
               if (trim(adjustl(current_snippet)) == string) then
!
                  found = .true.
                  return
!
               endif
!
            endif
!
         else ! We have found a comma
!
            current_cursor = 1 ! Reset cursor
!
            if (trim(adjustl(current_snippet)) == string) then
!
               found = .true.
               return
!
            endif
!
            current_snippet = repeat(' ', 200)
!
         endif
!
      enddo
!
   end function is_string_in_cs_list_input_file
!
!
   subroutine read_keywords_and_geometry_input_file(the_file)
!!
!!    Read all
!!    Written by sarai D. Folkestad, Mar 2020
!!
!!       Includes the former routine
!!
!!       print_to_output_input_file, witten by Eirik F. Kjønstad, Jan 2020
!!
!!    1) Prints the input file - except for the geometry specification - to the output file.
!!
!!    2) Places input file in memory
!!
!!       Places the input file, excluding the geometry,
!!       in to the array input_
!!
!!       Places the geometry, excluding the MM geometry,
!!       in to the array geometry
!!
!!       Places the MM geometry, if present,
!!       in to the array mm_geometry
!!
!!    3) Checks for illegal keywords and sections in the input
!!
      implicit none
!
      class(input_file), intent(inout) :: the_file
!
      character(len=200) :: line
!
      integer :: i, io_error
!
      logical :: QM_atoms
!
      call output%printf('m', ':: Input file', fs='(//t3,a)')
      call output%print_separator('m', 16, '=', fs='(t3,a/)')
!
      call output%printf('m', 'Note: geometry section is excluded from this print', fs='(t6,a/)')
!
      rewind(the_file%unit_)
!
      the_file%n_keyword_lines = 0
!
      do
!
         read(the_file%unit_, '(a)', iostat=io_error) line
!
         if (trim(adjustl(line)) == 'geometry') then
!
            exit
!
         elseif (io_error .ne. 0) then
!
            call output%error_msg("The 'geometry' section appears to be missing in the input file.")
!
         endif
!
         call output%printf('m', line, fs='(t6,a)', ll=120)
         the_file%n_keyword_lines = the_file%n_keyword_lines + 1
!
      enddo
!
!     Continue reading the geometry to count QM and MM atoms
!
      QM_atoms = .true.
      the_file%n_geometry_lines = 0
      the_file%n_mm_atom_lines = 0
!
      do
!
         read(the_file%unit_, '(a)', iostat=io_error) line
!
         if (trim(adjustl(line)) == 'end geometry') then
!
            exit
!
         endif
!
         if (trim(adjustl(line)) == '--') then ! QM and MM atoms separator
!
            QM_atoms = .false.
            cycle
!
         endif
!
         if (QM_atoms) then
!
            the_file%n_geometry_lines = the_file%n_geometry_lines + 1
!
         else
!
            the_file%n_mm_atom_lines = the_file%n_mm_atom_lines + 1
!
         endif
!
      enddo
!
!     Place keyword lines in memory
!
      allocate(the_file%input_(the_file%n_keyword_lines))
!
      rewind(the_file%unit_)
!
      do i = 1, the_file%n_keyword_lines
!
         read(the_file%unit_, '(a)', iostat=io_error) line
!
         line = adjustl(line)
         call convert_to_lowercase(line)
!
         the_file%input_(i) = line
!
      enddo
!
!     Place geometries in memory
!
      if (the_file%n_geometry_lines < 1) call output%error_msg('No QM atoms in geometry!')
!
      allocate(the_file%geometry(the_file%n_geometry_lines))
!
      read(the_file%unit_, '(a)', iostat=io_error) line ! Reads 'geometry' line
!
!     QM geometry
!
      do i = 1, the_file%n_geometry_lines
!
         read(the_file%unit_, '(a)', iostat=io_error) line
!
         line = adjustl(line)
         call convert_to_lowercase(line)
!
         the_file%geometry(i) = line
!
      enddo
!
!     MM geometry if present
!
      if (the_file%n_mm_atom_lines > 0) then
!
         read(the_file%unit_, '(a)', iostat=io_error) line ! Reads '--' line
!
         allocate(the_file%mm_geometry(the_file%n_mm_atom_lines))
!
         do i = 1, the_file%n_mm_atom_lines
!
            read(the_file%unit_, '(a)', iostat=io_error) line
!
            line = adjustl(line)
            call convert_to_lowercase(line)
!
            the_file%mm_geometry(i) = line
!
         enddo
      endif
!
      call the_file%check_for_errors() ! Check for incorrect/missing keywords/sections
!
   end subroutine read_keywords_and_geometry_input_file
!
!
   subroutine get_section_limits_input_file(the_file, string, start_, end_)
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
      integer, intent(out) :: start_, end_
!
      integer :: n_beginnings, n_ends, i
!
!     Find the number of instances of the section,
!     stopping with an error if there are any inconsistencies
!
      n_ends = 0
      n_beginnings = 0
!
      do i = 1, the_file%n_keyword_lines !while (trim(line) /= 'end geometry')
!
         if (trim(the_file%input_(i)) == 'end ' // string) then
!
            n_ends = n_ends + 1
            end_ = i
!
         endif

         if (trim(the_file%input_(i)) == string) then
!
            n_beginnings = n_beginnings + 1
            start_ = i
!
         endif
!
      enddo
!
      if (n_beginnings > 1)                                 &
         call output%error_msg('Tried to move to section "' &
            // string // '" with more than one starting clause.')
!
      if (n_ends > 1)                                       &
         call output%error_msg('Tried to move to section "' &
            // string // '" with more than one ending clause.')
!
      if (n_ends == 0 .and. n_beginnings == 0)                             &
         call output%error_msg('Tried to move to non-existent section "'   &
            // string // '".')
!
      if (n_ends < 1) &
         call output%error_msg('Tried to move to section "' // string // '" with no end.')
!
      if (n_beginnings < 1) &
         call output%error_msg('Tried to move to section "' // string // '" with no beginning.')
!
   end subroutine get_section_limits_input_file
!
!
   subroutine cleanup_geometry_input_file(the_file)
!!
!!    Cleanup geometries
!!    Written by Sarai D. Folkestad, Mar 2020
!!
!!    Deallocates the geometry arrays
!!
      implicit none
!
      class(input_file), intent(inout) :: the_file
!
      deallocate(the_file%geometry)
!
      if (allocated(the_file%mm_geometry)) deallocate(the_file%mm_geometry)
!
   end subroutine cleanup_geometry_input_file
!
!
   subroutine cleanup_keywords_input_file(the_file)
!!
!!    Cleanup keywords
!!    Written by Sarai D. Folkestad, Mar 2020
!!
!!    Deallocates the keyword arrays
!!
      implicit none
!
      class(input_file), intent(inout) :: the_file
!
      deallocate(the_file%input_)
!
   end subroutine cleanup_keywords_input_file
!
!
end module input_file_class
