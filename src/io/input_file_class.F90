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
      type(section), allocatable, private :: sections(:)
!
      character(len=30), allocatable, private :: rf_wfs(:)
      character(len=30), allocatable, private :: cc_wfs(:)
!
      integer, private :: n_keyword_lines ! Number of lines excluding the geometry section
      integer, private :: n_qm_atom_lines ! Number of QM atoms
      integer, private :: n_mm_atom_lines ! Number of MM atoms
!
!     Array of the input lines excluding geometry
!
      character(len=200), dimension(:), allocatable, private :: input_
!
!     Array of the QM and MM input geometries
!
      character(len=200), dimension(:), allocatable, private :: geometry
      character(len=200), dimension(:), allocatable, private :: mm_geometry
!
   contains
!
      procedure :: open_                                    &
                => open_input_file
!
      procedure :: close_                                   &
                => close_input_file
!
      generic :: get_keyword                                &
              => get_integer4_keyword,                      &
                 get_integer8_keyword,                      &
                 get_string_keyword,                        &
                 get_real_dp_keyword
!
      generic :: get_required_keyword                       &
              => get_required_string_keyword,               &
                 get_required_integer4_keyword,             &
                 get_required_integer8_keyword,             &
                 get_required_real_dp_keyword
!
      generic :: get_array_for_keyword                      &
              => get_integer_array_for_keyword,             &
                 get_real_array_for_keyword
!
      procedure :: check_for_errors                         &
                => check_for_errors_input_file
!
      procedure :: is_section_present                       &
                => is_section_present_input_file
!
      procedure :: is_keyword_present                       &
                => is_keyword_present_input_file
!
      procedure :: get_n_elements_for_keyword               &
                => get_n_elements_for_keyword
!
      procedure :: get_n_atoms                              &
                => get_n_atoms_input_file
!
      procedure :: get_n_mm_atoms                           &
                => get_n_mm_atoms_input_file
!
      procedure :: get_n_mm_molecules                       &
                => get_n_mm_molecules_input_file
!
      procedure :: get_geometry                             &
                => get_geometry_input_file
!
      procedure :: get_mm_geometry_fq                       &
                => get_mm_geometry_fq_input_file
!
      procedure :: get_mm_geometry_non_polarizable          &
                => get_mm_geometry_non_polarizable_input_file
!
      procedure :: get_reference_wf                         &
                => get_reference_wf_input_file
!
      procedure :: get_cc_wf                                &
                => get_cc_wf_input_file
!
      procedure :: requested_reference_calculation          &
                => requested_reference_calculation_input_file
!
      procedure :: requested_cc_calculation                 &
                => requested_cc_calculation_input_file
!
      procedure :: is_string_in_cs_list                     &
                => is_string_in_cs_list_input_file
!
      procedure :: read_keywords_and_geometry               &
                => read_keywords_and_geometry_input_file
!
      procedure :: cleanup_geometry                         &
                => cleanup_geometry_input_file
!
      procedure :: cleanup_keywords                         &
                => cleanup_keywords_input_file
!
      procedure :: is_embedding_on                          &
                => is_embedding_on_input_file
!
      procedure :: place_records_in_memory                  &
                => place_records_in_memory_input_file
!
      procedure, private :: get_wf
      procedure, private :: requested_calculation
!
      procedure, private :: get_string_keyword_wo_safety
      procedure, private :: get_section_limits
      procedure, private :: check_section_for_illegal_keywords
      procedure, private :: check_for_illegal_sections
      procedure, private :: print_sections
!
      procedure, nopass, private :: string_is_comment
      procedure, nopass, private :: extract_keyword_from_string
      procedure, nopass, private :: extract_keyword_value_from_string
!
      procedure, private :: get_integer4_keyword
      procedure, private :: get_integer8_keyword
      procedure, private :: get_string_keyword
      procedure, private :: get_real_dp_keyword
      procedure, private :: get_required_string_keyword
      procedure, private :: get_required_integer4_keyword
      procedure, private :: get_required_integer8_keyword
      procedure, private :: get_required_real_dp_keyword
!
      procedure, private :: get_integer_array_for_keyword
      procedure, private :: get_real_array_for_keyword
!
      procedure, private :: section_is_allowed
      procedure, private :: keyword_is_allowed
      procedure, private :: check_keyword_and_section
!
      procedure, public :: get_n_state_guesses
      procedure, public :: get_state_guesses
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
   function new_input_file(name_) result(this)
!!
!!    New input file
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2018
!!
      implicit none
!
      type(input_file) :: this
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
      type(section) :: qed
      type(section) :: tdhf
!
!     Set input file name, access and format
!
      this%name_ = name_
!
      this%access_ = 'sequential'
      this%format_ = 'formatted'
      this%action_ = 'read'
!
!     Set method section
!
      method%name_    = 'method'
      method%required = .false.
!
      this%rf_wfs = [character(len=30) :: &
                        'hf',     &
                        'uhf',    &
                        'mlhf',   &
                        'qed-hf', &
                        'cuhf',   &
                        'rohf']
!
      this%cc_wfs = [character(len=30) :: &
                        'ccs',            &
                        'mp2',            &
                        'cc2',            &
                        'lowmem-cc2',     &
                        'ccsd',           &
                        'cc3',            &
                        'ccsd(t)',        &
                        'mlcc2',          &
                        'mlccsd']
!
      allocate(method%keywords(size(this%rf_wfs) &
                             + size(this%cc_wfs)))
!
      method%keywords = [this%rf_wfs, this%cc_wfs]
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
                             'dipole length',      &
                             'eom',                &
                             'frequencies',        &
                             'initial states',     &
                             'lr',                 &
                             'polarizabilities',   &
                             'permanent moments',  &
                             'transition moments']
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
                              'cholesky eri',         &
                              'ground state',         &
                              'ground state geoopt',  &
                              'excited state',        &
                              'response',             &
                              'mean value',           &
                              'time dependent state', &
                              'restart',              &
                              'time dependent hf']
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
                           'eri storage',        &
                           'mo eri in memory',   &
                           'ri',                 &
                           't1 eri in memory']
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
                      'nto restart',            &
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
      qed%name_    = 'qed'
      qed%required = .false.
      qed%keywords = [character(len=30) :: &
                     'coherent state',     &
                     'coupling bilinear',  &
                     'coupling self',      &
                     'coupling',           &
                     'frequency',          &
                     'hf coherent state',  &
                     'modes',              &
                     'polarization',       &
                     'quadrupole oei',     &
                     'wavevector']
!
!
      global_print%name_    = 'print'
      global_print%required = .false.
      global_print%keywords = [character(len=30) ::  &
                              'output print level ', &
                              'timing print level ', &
                              'full references',     &
                              'z-matrix']
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
                                 'mo screening']
!
!
      solver_scf%name_    = 'solver scf'
      solver_scf%required = .false.
      solver_scf%keywords = [character(len=30) ::           &
                            'algorithm',                    &
                            'ao density guess',             &
                            'coulomb threshold',            &
                            'crop',                         &
                            'cumulative fock threshold',    &
                            'diis dimension',               &
                            'energy threshold',             &
                            'exchange threshold',           &
                            'gradient threshold',           &
                            'integral cutoff',              &
                            'integral precision',           &
                            'max iterations',               &
                            'one-electron integral cutoff', &
                            'print orbitals',               &
                            'restart',                      &
                            'rohf coupling parameters',     &
                            'skip',                         &
                            'storage',                      &
                            'write molden']
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
                              'max reduced dimension',    &
                              'davidson preconvergence',  &
                              'preconvergence threshold', &
                              'max micro iterations',     &
                              'rel micro threshold',      &
                              'chain length',             &
                              'lanczos normalization',    &
                              'remove core',              &
                              'state guesses']
!
!
      solver_cc_gs%name_    = 'solver cc gs'
      solver_cc_gs%required = .false.
      solver_cc_gs%keywords = [character(len=30) ::       &
                              'algorithm',                &
                              'energy threshold',         &
                              'omega threshold',          &
                              'crop',                     &
                              'micro iteration storage',  &
                              'max micro iterations',     &
                              'multimodel newton',        &
                              'rel micro threshold',      &
                              'storage',                  &
                              'max iterations',           &
                              'diis dimension',           &
                              'restart' ]
!
!
      solver_cc_multipliers%name_    = 'solver cc multipliers'
      solver_cc_multipliers%required = .false.
      solver_cc_multipliers%keywords = [character(len=30) ::       &
                                       'algorithm',                &
                                       'threshold',                &
                                       'storage',                  &
                                       'crop',                     &
                                       'micro iteration storage',  &
                                       'max micro iterations',     &
                                       'multimodel newton',        &
                                       'rel micro threshold',      &
                                       'diis dimension',           &
                                       'restart',                  &
                                       'max reduced dimension',    &
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
      system%keywords = [character(len=30) ::   &
                        'name',                 &
                        'cartesian gaussians',  &
                        'pure gaussians',       &
                        'charge',               &
                        'multiplicity']
!
!
      tdhf%name_    = 'solver tdhf'
      tdhf%required = .false.
      tdhf%keywords = [character(len=30) ::    &
                       'energy threshold',     &
                       'max iterations',       &
                       'max reduced dimension',&
                       'residual threshold',   &
                       'restart',              &
                       'storage',              &
                       'states',               &
                       'tamm-dancoff']
!
!
      visualization%name_    = 'visualization'
      visualization%required = .false.
      visualization%keywords = [character(len=30) ::        &
                               'file format',               &
                               'grid spacing',              &
                               'grid buffer',               &
                               'grid min',                  &
                               'grid max',                  &
                               'plot cc density',           &
                               'plot es densities',         &
                               'plot hf orbitals',          &
                               'plot hf density',           &
                               'plot hf active density',    &
                               'plot cntos',                &
                               'plot ntos',                 &
                               'nto threshold',             &
                               'plot transition densities', &
                               'states to plot']
!
!     Gather all sections into the file's section array
!
      this%sections = [active_atoms,              &
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
                       memory,                    &
                       method,                    &
                       mlcc,                      &
                       mlhf,                      &
                       mm,                        &
                       pcm,                       &
                       qed,                       &
                       solver_cc_es,              &
                       solver_cc_gs,              &
                       solver_cc_multipliers,     &
                       solver_cc_propagation,     &
                       solver_cc_response,        &
                       solver_cholesky,           &
                       solver_fft_dipole_moment,  &
                       solver_fft_electric_field, &
                       solver_scf,                &
                       solver_scf_geoopt,         &
                       system,                    &
                       tdhf,                      &
                       visualization]
!
      this%is_open = .false.
      this%unit_ = -1
!
   end function new_input_file
!
!
   subroutine open_input_file(this)
!!
!!    Open input file
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(input_file) :: this
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      if (this%is_open) then
!
         call output%error_msg(trim(this%name_)//' is already open')
!
      endif
!
      open(newunit=this%unit_, file=this%name_, access=this%access_, &
           action='read', status='unknown', form=this%format_, iostat=io_error, iomsg=io_msg)
!
      if (io_error .ne. 0) then
!
         call output%error_msg('could not open eT input file '//trim(this%name_)//&
                              &'error message: '//trim(io_msg))
      endif
!
      this%is_open = .true.
!
   end subroutine open_input_file
!
!
   subroutine close_input_file(this)
!!
!!    Close input file
!!    Written by Rolf H. Myhre, May 2019
!!
      implicit none
!
      class(input_file) :: this
!
      integer              :: io_error
      character(len=100)   :: io_msg
!
      if (.not. this%is_open) then
         call output%error_msg(trim(this%name_)//' already closed')
      end if
!
      close(this%unit_, iostat=io_error, iomsg=io_msg, status='keep')
!
      if (io_error .ne. 0) then
!
         call output%error_msg('could not close eT input file '//trim(this%name_)//&
                              &'error message: '//trim(io_msg))
      endif
!
      this%is_open = .false.
      this%unit_ = -1
!
   end subroutine close_input_file
!
!
   subroutine check_for_errors_input_file(this)
!!
!!    Check for errors
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
!!    Looks for errors the user may have made when writing the input file.
!!
      implicit none
!
      class(input_file) :: this
!
      integer :: k
!
!     Look for illegal sections
!
      call this%check_for_illegal_sections()
!
!     For each legal section present, look for illegal keywords within that section
!
      do k = 1, size(this%sections)
!
         if (this%is_section_present(this%sections(k)%name_)) then
!
            call this%check_section_for_illegal_keywords(this%sections(k))
!
         else
!
            if (this%sections(k)%required) then
!
               call output%printf('m', 'All calculations require the section "' // &
                                  trim(this%sections(k)%name_) // '". It &
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
   subroutine check_for_illegal_sections(this)
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
      class(input_file) :: this
!
      integer :: start_, end_, i
!
      logical :: recognized
!
      character(len=200) :: section
!
      do i = 1, this%n_keyword_lines
!
         if (this%input_(i)(1 : 3) == 'end') then
!
!           Located the end of a section -
!           attempt to move to the beginning of that section => fails if inconsistent beginning and end
!
            section = trim(adjustl(this%input_(i)(4 : 200)))
!
            call this%get_section_limits(section, start_, end_)
!
!           Check whether section name is valid
!
            recognized = this%section_is_allowed(section)
!
            if (.not. recognized) then
!
               call output%printf('m', 'Could not recognize section named "(a0)".', &
                                       chars=[section], fs='(/t3,a)')
!
               call this%print_sections()
!
               call output%error_msg('Something is wrong in the input file. See above.')
!
            endif
!
         endif
!
      enddo
!
   end subroutine check_for_illegal_sections
!
!
   subroutine print_sections(this)
!!
!!    Print sections
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
      implicit none
!
      class(input_file), intent(in) :: this
!
      integer :: k
!
      call output%printf('m', 'The valid input sections are:', fs='(/t3,a/)')
!
      do k = 1, size(this%sections)
!
         call output%printf('m', '(a0)', chars=[trim(this%sections(k)%name_)], fs='(t6,a)')
!
      enddo
!
   end subroutine print_sections
!
!
   subroutine check_section_for_illegal_keywords(this, the_section)
!!
!!    Checks a section for illegal keywords
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
!!    If there is an illegal keyword in the section, it stops with an error.
!!    If a legal keyword is specified more than once, it stops with an error.
!!
      implicit none
!
      class(input_file) :: this
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
      call this%get_section_limits(the_section%name_, start_, end_)
!
      do record = start_ + 1, end_ - 1
!
         recognized = .false.
!
         if (.not. this%string_is_comment(this%input_(record))) then
!
            call this%extract_keyword_from_string(this%input_(record), keyword)
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
   end subroutine check_section_for_illegal_keywords
!
!
   function requested_reference_calculation_input_file(this) result(requested)
!!
!!    Requested reference calculation
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
      implicit none
!
      class(input_file), intent(in) :: this
!
      logical :: requested
!
      requested = this%requested_calculation(this%rf_wfs)
!
   end function requested_reference_calculation_input_file
!
!
   function get_reference_wf_input_file(this) result(ref_wf)
!!
!!    Get reference wavefunction
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
      implicit none
!
      class(input_file), intent(in) :: this
!
      character(len=30) :: ref_wf
!
      ref_wf = this%get_wf(this%rf_wfs)
!
   end function get_reference_wf_input_file
!
!
   function requested_calculation(this, wfs) result(requested)
!!
!!    Requested calculation
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
!!    Looks for wfs(k) in method and returns true if it finds one such wf
!!    in the input file. Gives an error if it finds more than one of wfs(k).
!!
      implicit none
!
      class(input_file), intent(in) :: this
!
      character(len=30), dimension(:), intent(in) :: wfs
!
      logical :: requested
!
      integer :: n_wfs, k
!
      n_wfs = 0
      do k = 1, size(wfs)
!
         if (this%is_keyword_present(wfs(k), 'method')) then
!
            n_wfs = n_wfs + 1
!
         endif
!
      enddo
!
      if (n_wfs == 1) then
!
         requested = .true.
!
      elseif (n_wfs > 1) then
!
         requested = .false.
         call output%error_msg('Requested more than one method &
                               & of a single kind (e.g., two CC methods).')
!
      else
!
         requested = .false.
!
      endif
!
   end function requested_calculation
!
!
   function requested_cc_calculation_input_file(this) result(requested)
!!
!!    Requested CC calculation
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
      implicit none
!
      class(input_file), intent(in) :: this
!
      logical :: requested
!
      requested = this%requested_calculation(this%cc_wfs)
!
   end function requested_cc_calculation_input_file
!
!
   function get_cc_wf_input_file(this) result(cc_wf)
!!
!!    Get CC wavefunction
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
      implicit none
!
      class(input_file), intent(in) :: this
!
      character(len=30) :: cc_wf
!
      cc_wf = this%get_wf(this%cc_wfs)
!
   end function get_cc_wf_input_file
!
!
   function get_wf(this, wfs) result(wf)
!!
!!    Get wavefunction
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
!!    Looks for wfs(k) in the 'method' section of the input.
!!
      implicit none
!
      class(input_file), intent(in) :: this
!
      character(len=30), dimension(:), intent(in) :: wfs
!
      character(len=30) :: wf
!
      integer :: k
!
      do k = 1, size(wfs)
!
         if (this%is_keyword_present(wfs(k),'method')) then
!
            wf = wfs(k)
            return
!
         endif
!
      enddo
!
      call output%error_msg('Failed to read wavefunction from input. &
         & Did you forget to specify all wavefunctions for the calculation?')
!
   end function get_wf
!
!
!
!
   subroutine get_integer4_keyword(this, keyword, section, keyword_value)
!!
!!    Get integer(4) keyword
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
!!    If specified, reads keyword as an integer into keyword value.
!!
      implicit none
!
      class(input_file), intent(in) :: this
!
      character(len=*), intent(in) :: keyword
      character(len=*), intent(in) :: section
!
      integer(i32), intent(inout) :: keyword_value
!
      character(len=200) :: keyword_value_string
!
      call this%check_keyword_and_section(keyword, section)
!
      if (this%is_keyword_present(keyword, section)) then
!
!        Get the keyword value in string format
!
         call this%get_string_keyword_wo_safety(keyword, section, keyword_value_string)
!
!        Extract the integer from the string
!
         read(keyword_value_string, *) keyword_value
!
      endif
!
   end subroutine get_integer4_keyword
!
!
   subroutine get_integer8_keyword(this, keyword, section, keyword_value)
!!
!!    Get integer(8) keyword
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
!!    If specified, reads keyword as an integer into keyword value.
!!
      implicit none
!
      class(input_file), intent(in) :: this
!
      character(len=*), intent(in) :: keyword
      character(len=*), intent(in) :: section
!
      integer(i64), intent(inout) :: keyword_value
!
      character(len=200) :: keyword_value_string
!
      call this%check_keyword_and_section(keyword, section)
!
      if (this%is_keyword_present(keyword, section)) then
!
!        Get the keyword value in string format
!
         call this%get_string_keyword_wo_safety(keyword, section, keyword_value_string)
!
!        Extract the integer from the string
!
         read(keyword_value_string, *) keyword_value
!
      endif
!
   end subroutine get_integer8_keyword
!
!
   subroutine get_required_integer4_keyword(this, keyword, section, keyword_value)
!!
!!    Get required integer(4) keyword
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
!!    Reads keyword as an integer into keyword value. If the keyword or the
!!    section is not specified, an error occurs because the keyword is "required".
!!
      implicit none
!
      class(input_file), intent(in) :: this
!
      character(len=*), intent(in) :: keyword
      character(len=*), intent(in) :: section
!
      integer(i32), intent(inout) :: keyword_value
!
      character(len=200) :: keyword_value_string
!
      if (.not. this%is_section_present(section)) &
         call output%error_msg('could not find the required section: '// trim(section) // '.')
!
      if (.not. this%is_keyword_present(keyword, section)) &
         call output%error_msg('could not find the required keyword '// trim(keyword) // &
                               ' in section ' // trim(section))
!
!     Get the keyword value in string format
!
      call this%get_string_keyword_wo_safety(keyword, section, keyword_value_string)
!
!     Extract the integer from the string
!
      read(keyword_value_string, *) keyword_value
!
   end subroutine get_required_integer4_keyword
!
!
   subroutine get_required_integer8_keyword(this, keyword, section, keyword_value)
!!
!!    Get required integer(8) keyword
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
!!    Reads keyword as an integer into keyword value. If the keyword or the
!!    section is not specified, an error occurs because the keyword is "required".
!!
      implicit none
!
      class(input_file), intent(in) :: this
!
      character(len=*), intent(in) :: keyword
      character(len=*), intent(in) :: section
!
      integer(i64), intent(out) :: keyword_value
!
      character(len=200) :: keyword_value_string
!
      if (.not. this%is_section_present(section)) &
         call output%error_msg('could not find the required section: '// trim(section) // '.')
!
      if (.not. this%is_keyword_present(keyword, section)) &
         call output%error_msg('could not find the required keyword '// trim(keyword) // &
                               ' in section ' // trim(section))
!
!     Get the keyword value in string format
!
      call this%get_string_keyword_wo_safety(keyword, section, keyword_value_string)
!
!     Extract the integer from the string
!
      read(keyword_value_string, *) keyword_value
!
   end subroutine get_required_integer8_keyword
!
!
   subroutine get_real_dp_keyword(this, keyword, section, keyword_value)
!!
!!    Get real(dp) keyword
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
!!    If specified, reads keyword as a real double precision into keyword value.
!!
      implicit none
!
      class(input_file), intent(in) :: this
!
      character(len=*), intent(in) :: keyword
      character(len=*), intent(in) :: section
!
      real(dp), intent(inout) :: keyword_value
!
      character(len=200) :: keyword_value_string
!
      call this%check_keyword_and_section(keyword, section)
!
      if (this%is_keyword_present(keyword, section)) then
!
!        Get the keyword value in string format
!
         call this%get_string_keyword_wo_safety(keyword, section, keyword_value_string)
!
!        Extract the integer from the string
!
         read(keyword_value_string, *) keyword_value
!
      endif
!
   end subroutine get_real_dp_keyword
!
!
   subroutine get_required_real_dp_keyword(this, keyword, section, keyword_value)
!!
!!    Get required real(dp) keyword
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
!!    Reads keyword as an double precision into keyword value. If the keyword or the
!!    section is not specified, an error occurs because the keyword is "required".
!!
      implicit none
!
      class(input_file), intent(in) :: this
!
      character(len=*), intent(in) :: keyword
      character(len=*), intent(in) :: section
!
      real(dp), intent(out) :: keyword_value
!
      character(len=200) :: keyword_value_string
!
      if (.not. this%is_section_present(section)) &
         call output%error_msg('could not find the required section: '// trim(section) // '.')
!
      if (.not. this%is_keyword_present(keyword, section)) &
         call output%error_msg('could not find the required keyword '// trim(keyword) // &
                               ' in section ' // trim(section))
!
!     Get the keyword value in string format
!
      call this%get_string_keyword_wo_safety(keyword, section, keyword_value_string)
!
!     Extract the integer from the string
!
      read(keyword_value_string, *) keyword_value
!
   end subroutine get_required_real_dp_keyword
!
!
   subroutine get_string_keyword(this, keyword, section, keyword_value)
!!
!!    Get string keyword
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
!!    If specified, reads keyword as a string into keyword value.
!!
      implicit none
!
      class(input_file), intent(in) :: this
!
      character(len=*), intent(in) :: keyword
      character(len=*), intent(in) :: section
!
      character(len=200) :: keyword_value
!
      call this%check_keyword_and_section(keyword, section)
!
      if (this%is_keyword_present(keyword, section)) then
!
!        Get the keyword value in string format
!
         call this%get_string_keyword_wo_safety(keyword, section, keyword_value)
!
      endif
!
   end subroutine get_string_keyword
!
!
   subroutine get_required_string_keyword(this, keyword, section, keyword_value)
!!
!!    Get required string keyword
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
!!    If specified, reads keyword as a string into keyword value.
!!
      implicit none
!
      class(input_file), intent(in) :: this
!
      character(len=*), intent(in) :: keyword
      character(len=*), intent(in) :: section
!
      character(len=200) :: keyword_value
!
      if (.not. this%is_section_present(section)) &
         call output%error_msg('could not find the required section: '// trim(section) // '.')
!
      if (.not. this%is_keyword_present(keyword, section)) &
         call output%error_msg('could not find the required keyword '// trim(keyword) // ' in section ' // trim(section))
!
!     Get the keyword value in string format
!
      call this%get_string_keyword_wo_safety(keyword, section, keyword_value)
!
   end subroutine get_required_string_keyword
!
!
   subroutine get_string_keyword_wo_safety(this, keyword, section, keyword_value)
!!
!!    Get string keyword without safety
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
!!    Reads keyword in section, placing the result in the string keyword_value. This
!!    routine gives an error if the keyword is not located ("no safety").
!!
      implicit none
!
      class(input_file), intent(in) :: this
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
      call this%get_section_limits(section, start_, end_)
!
!     Loop through records within the section to locate & get the keyword value
!
      do record = start_ + 1, end_ - 1
!
         if (this%string_is_comment(this%input_(record))) then
!
            cycle
!
         else
!
            call this%extract_keyword_from_string(this%input_(record), local_keyword)
!
            if (trim(local_keyword) == keyword) then
!
               call this%extract_keyword_value_from_string(this%input_(record), &
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
   end subroutine get_string_keyword_wo_safety
!
!
   pure function string_is_comment(string) result(is_comment)
!!
!!    String is comment?
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
      implicit none
!
      character(len=200), intent(in) :: string
!
      logical :: is_comment
!
      character(len=200) :: tmp_string
!
      tmp_string = adjustl(string)
!
      if (tmp_string(1:1) == '!') then
!
         is_comment = .true.
!
      else
!
         is_comment = .false.
!
      endif
!
   end function string_is_comment
!
!
   subroutine extract_keyword_from_string(string, keyword)
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
   end subroutine extract_keyword_from_string
!
!
   subroutine extract_keyword_value_from_string(string, keyword_value)
!!
!!    Extract keyword value from string
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
!!    Note: assumes that the string is not a comment.
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
   end subroutine extract_keyword_value_from_string
!
!
   function is_keyword_present_input_file(this, keyword, section) result(is_present)
!!
!!    Is keyword present?
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
!!    Returns true if the keyword is in the section.
!!
      implicit none
!
      class(input_file), intent(in) :: this
!
      character(len=*), intent(in) :: keyword
      character(len=*), intent(in) :: section
!
      logical :: is_present
!
      integer :: record, start_, end_
!
      character(len=200) :: local_keyword
!
      if (.not. this%keyword_is_allowed(keyword, section)) then
         call output%error_msg('Trying to read keyword ((a0)) that is not defined &
                               &in section: ' // section,  chars=[keyword])
      end if
!
!     Move to the requested section & get the number of records in that section
!
      if (.not. this%is_section_present(section)) then
!
         is_present = .false.
         return
!
      endif
!
      call this%get_section_limits(section, start_, end_)
!
!     Loop through records within the section & try to locate the keyword
!
      do record = start_ + 1, end_ - 1
!
         if (this%string_is_comment(this%input_(record))) then
!
            cycle
!
         else
!
            call this%extract_keyword_from_string(this%input_(record), local_keyword)
!
            if (trim(local_keyword) == keyword) then
!
               is_present = .true.
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
      is_present = .false.
!
   end function is_keyword_present_input_file
!
!
   pure function is_section_present_input_file(this, section) result(is_present)
!!
!!    Is section present?
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
!!    Returns true if the section exists, false if it doesn't.
!!
      implicit none
!
      class(input_file), intent(in) :: this
!
      character(len=*), intent(in) :: section
!
      logical :: is_present
!
      integer :: i
!
      is_present = .false.
!
      do i = 1, this%n_keyword_lines
!
         if (trim(this%input_(i)) == section) then
!
            is_present = .true.
            return
!
         endif
!
      enddo
!
   end function is_section_present_input_file
!
!
   function get_n_elements_for_keyword(this, keyword, section) result(n_elements)
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
!!    get_array_for_keyword
!!
      implicit none
!
      class(input_file), intent(in) :: this
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
      if (.not. this%is_keyword_present(keyword, section)) return
!
!     Get the keyword value in string format
!
      call this%get_keyword(keyword, section, keyword_value_string)
!
!     Use string utility functionality to get n_elements
!
      n_elements = get_n_elements_in_string(keyword_value_string)
!
   end function get_n_elements_for_keyword
!
!
   subroutine get_integer_array_for_keyword(this, keyword, section, n_elements, array_)
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
!!    get_n_elements_for_keyword is called
!!    in order to determine n_elements so that array_
!!    can be allocated.
!!
      implicit none
!
      class(input_file), intent(in) :: this
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
      if (.not. this%is_keyword_present(keyword, section)) return
!
!     Get the keyword value in string format
!
      call this%get_keyword(keyword, section, keyword_value_string)
!
!     Use string utility functionality to get the array
!
      call get_elements_in_string(keyword_value_string, n_elements, array_)
!
   end subroutine get_integer_array_for_keyword
!
!
   subroutine get_real_array_for_keyword(this, keyword, section, n_elements, array_)
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
!!    get_n_elements_for_keyword is called
!!    in order to determine n_elements so that array_
!!    can be allocated.
!!
      implicit none
!
      class(input_file), intent(in) :: this
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
      call this%get_keyword(keyword, section, keyword_value_string)
!
!     Use string utility functionality to get the real array
!
      call get_reals_in_string(keyword_value_string, n_elements, array_)
!
   end subroutine get_real_array_for_keyword
!
!
   pure function get_n_atoms_input_file(this) result(n_atoms)
!!
!!    Get n atoms
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
!!    Reads the geometry section of the input and
!!    counts the number of atoms
!!
      implicit none
!
      class(input_file), intent(in) :: this
!
      integer :: n_atoms
!
!     Local variables
!
      integer :: record
!
      n_atoms = 0
!
      do record = 1, this%n_qm_atom_lines
!
         if (this%geometry(record)(1:6) .ne. 'basis:' &
            .and. (this%geometry(record)(1:6) .ne. 'units:') &
            .and. (this%geometry(record)(1:5) .ne. 'ghost')) &
            n_atoms = n_atoms + 1
!
      enddo
!
   end function  get_n_atoms_input_file
!
!
   pure function get_n_mm_atoms_input_file(this) result(n_atoms)
!!
!!    Get n MM atoms
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
!!    Reads the geometry section of the input and
!!    counts the number of atoms
!!
      implicit none
!
      class(input_file), intent(in) :: this
!
      integer :: n_atoms
!
      n_atoms = this%n_mm_atom_lines
!
   end function  get_n_mm_atoms_input_file
!
!
   function get_n_mm_molecules_input_file(this) result(n_molecules)
!!
!!    Get n MM molecules
!!    Written by Sarai D. Folkestad, Sep 2020
!!
!!    Reads the MM geometry section of the input and
!!    counts the number of molecules
!!
      implicit none
!
      class(input_file), intent(in) :: this
!
      integer :: n_molecules
!
!     Local variables
!
      integer :: cursor, record, current_molecule, previous_molecule
!
      character(len=200) :: string, imolecule
!
!     Loop through the MM geometry section
!
      previous_molecule = 1
      n_molecules       = 1
!
      do record = 1, this%n_mm_atom_lines
!
         string = (trim(adjustl(this%mm_geometry(record))))
         cursor = first_instance_of_character(string,'=')
!
         string = string(cursor+1:200)
!
         cursor = first_instance_of_character(string,']')
!
         imolecule = string(1:cursor-1)
         read(imolecule,'(i4)') current_molecule
!
         if (record == 1 .and. current_molecule .ne. 1) &
            call output%error_msg('Molecule 1 is missing!')
!
         if (current_molecule .eq. previous_molecule + 1) then
!
            n_molecules = n_molecules + 1
!
         elseif (current_molecule .ne. previous_molecule .and. &
               current_molecule .ne. previous_molecule + 1) then
!
               call output%error_msg('Molecule (i0) missing even though (i0) &
                                    &should be present in QM/MM calculation.', &
                                    ints=[previous_molecule + 1, n_molecules])
!
         endif
!
         previous_molecule = current_molecule
!
      enddo
!
   end function get_n_mm_molecules_input_file
!
!
   subroutine get_geometry_input_file(this, n_atoms, symbols,  &
                                      positions, basis_sets,   &
                                      charge, units_angstrom,  &
                                      is_ghost)
!!
!!    Get geometry
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!    Modified by Åsmund H. Tveten, Oct 2019. Generalized to Bohr units.
!!    Modified by Tor S. Haugland, May 2021. Added ghost atoms.
!!    Modified by SDF, Jun 2021. Added atomic charges.
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
      class(input_file), intent(in) :: this
!
      integer, intent(in) :: n_atoms
!
      character(len=2), dimension(n_atoms), intent(out)   :: symbols
      character(len=100), dimension(n_atoms), intent(out) :: basis_sets
      logical, dimension(n_atoms), intent(out) :: is_ghost
!
      real(dp), dimension(3, n_atoms), intent(out) :: positions ! x, y, z
      integer, dimension(n_atoms), intent(out) :: charge
!
      logical, intent(out) :: units_angstrom ! True if units are Ångström/unspecified, false if Bohr
!
!     Local variables
!
      integer :: current_atom, i
!
      integer :: start_
!
      character(len=200) :: string
      character(len=100) :: current_basis
      logical :: is_ghost_atom = .false.
!
      integer :: cursor
!
      start_ = 1 ! Specifies the line of the first and required basis
!
!     Are units specified?
!     Note that units can only be specified as the first line of the geometry
!
      units_angstrom = .true. ! Default units are Angstrom
      charge = 0 ! Default charge is zero
!
      if (this%geometry(1)(1:6) == 'units:') then
!
         string = (trim(adjustl(this%geometry(1)(7:200))))
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
      if(this%geometry(start_)(1:6) /= 'basis:') &
            call output%error_msg('did not find basis in geometry section.')
!
!     Loop through geometry
!
      current_atom = 0
!
      do i = start_, this%n_qm_atom_lines
!
         string = trim(adjustl(this%geometry(i)))
!
         if (string(1:6) == 'units:') &
            call output%error_msg('Units must be specified as the first line in the geometry section.')
!
         if(string(1:6) == 'basis:') then
!
            current_basis = trim(adjustl(string(7:200)))
!
         elseif (string(1:5) == 'ghost') then
!
!           Every atom after 'ghost' keyword are ghosts
            is_ghost_atom = .true.
!
         else
!
            current_atom = current_atom + 1
!
            basis_sets(current_atom) = current_basis
            symbols(current_atom)    = string(1:2)
            is_ghost(current_atom)   = is_ghost_atom
!
            string = adjustl(string(3:200))
!
            if (is_substring_in_string(string, 'q')) then
!
                  cursor = set_cursor_to_substring(string, 'q')
                  read(string(1 : cursor - 1), *) positions(:, current_atom)
!
                  string = adjustl(string(cursor+1:200))
!
                  if (string(1:1) .ne. '=') &
                        call output%error_msg('in asignment of charge to atom in the geometry')
!
                  string = adjustl(string(2:200))
                  read(string, *) charge(current_atom)
            else
!
                  read(string(1:), *) positions(:, current_atom)
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
   end subroutine get_geometry_input_file
!
!
   subroutine get_mm_geometry_fq_input_file(this, n_atoms, n_molecules, &
                                          n_atoms_per_molecule, symbols,   &
                                          positions, chi, eta)
!!
!!    Get MM geometry FQ
!!    Written by Tommaso Giovanini, May 2019
!!    Modified for robustness by Sarai D. Folkestad 2020
!!
!!    Reads the geometry of the MM portion from the input file in the case
!!    of an fq calculation and sets it in the list of atoms.
!!
!!    Note: In order to be run, you need to know the number of MM atoms and molecules
!!          in the system
!!
      implicit none
!
      class(input_file),   intent(in) :: this
      integer,             intent(in) :: n_atoms
      integer,             intent(in) :: n_molecules
!
      character(len=2), dimension(n_atoms),     intent(out) :: symbols
      real(dp),         dimension(3, n_atoms),  intent(out) :: positions
      real(dp),         dimension(n_atoms),     intent(out) :: chi
      real(dp),         dimension(n_atoms),     intent(out) :: eta
      integer,          dimension(n_molecules), intent(out) :: n_atoms_per_molecule
!
!     Local variables
!
      integer :: record, cursor, current_atom, i, current_molecule
!
      character(len=200) :: string, coordinate
!
      current_atom = 0
!
      current_molecule     = 0
      n_atoms_per_molecule = 0
      cursor = 0
!
!     Loop through the MM atoms specified on input
!
      do record = 1, this%n_mm_atom_lines
!
         current_atom = current_atom + 1
         string = trim(adjustl(this%mm_geometry(record)))
!
!        Determine symbol
!
         symbols(current_atom) = string(1:2)
         string = adjustl(string(3:200))
!
!        Determine molecule index
!
         if (is_substring_in_string(string, 'mol')) then
            cursor = set_cursor_to_substring(string, 'mol')
         else
            call output%error_msg('could not find mol in MM geometry')
         endif
!
         string = adjustl(string(cursor+1:200))
!
         if (string(1:1) .ne. '=') call output%error_msg('in mol specification in MM geometry')
         string = adjustl(string(2:200))
!
         cursor = first_instance_of_character(string, ']')
!
         read(string(1:cursor-1),*) current_molecule
         string = adjustl(string(cursor+1:200))
!
         n_atoms_per_molecule(current_molecule) = n_atoms_per_molecule(current_molecule) + 1
!
         cursor = first_instance_of_character(string, '[')
!
         coordinate = string(1:cursor-1)
         read(coordinate, *) positions(:, current_atom)
!
         string = adjustl(string(cursor+1:200))
!
!        Determine chi
!
         if (is_substring_in_string(string, 'chi')) then
            cursor = set_cursor_to_substring(string, 'chi')
         else
            call output%error_msg('could not find chi in MM geometry')
         endif
!
         string = adjustl(string(cursor+1:200))
!
         if (string(1:1) .ne. '=') call output%error_msg('in chi specification in MM geometry')
         string = adjustl(string(2:200))
!
         cursor = first_instance_of_character(string,',')
         read(string(1:cursor-1), * ) chi(current_atom)
         string = adjustl(string(cursor+1:200))
!
!        Determine eta
!
         if (is_substring_in_string(string, 'eta')) then
            cursor = set_cursor_to_substring(string, 'eta')
         else
            call output%error_msg('could not find eta in MM geometry')
         endif
!
         string = trim(adjustl(string(cursor+1:200)))
!
         if (string(1:1) .ne. '=') call output%error_msg('in chi specification in MM geometry')
         string = adjustl(string(2:200))
!
         cursor = first_instance_of_character(string,']')
         read(string(1:cursor-1), * ) eta(current_atom)
!
         if(abs(eta(current_atom)).lt.1.0d-8) then
!
            call output%error_msg('You have put zero chemical hardness on atom: (i0)', &
                                  ints=[current_atom])
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
   end subroutine get_mm_geometry_fq_input_file
!
!
   subroutine get_mm_geometry_non_polarizable_input_file(this, n_atoms, symbols,   &
                                                         positions, charge)
!!
!!    Get MM geometry non-polarizable
!!    Written by Tommaso Giovanini, May 2019
!!    Modified for robustness by Sarai D. Folkestad 2020
!!
!!    Reads the geometry of the MM portion from the input file and
!!    sets it in the list of atoms.
!!
!!    Note: In order to be run, you need to know the number of MM atoms
!!          there are in the system
!!
      implicit none
!
      class(input_file), intent(in) :: this
!
      integer, intent(in) :: n_atoms
!
      character(len=2), dimension(n_atoms),     intent(out) :: symbols
      real(dp),         dimension(3, n_atoms),  intent(out) :: positions
      real(dp),         dimension(n_atoms),     intent(out) :: charge
!
!     Local variables
!
      integer :: record, cursor, current_atom, i, current_molecule
!
      character(len=200) :: string, coordinate
!
!     Loop through the MM atoms specified on input
!
      current_atom      = 0
      current_molecule  = 0
      cursor            = 0
!
      do record = 1, this%n_mm_atom_lines
!
         current_atom = current_atom + 1
         string = trim(adjustl(this%mm_geometry(record)))
!
!        Determine symbol
!
         symbols(current_atom) = string(1:2)
         string = adjustl(string(3:200))
!
!        Determine molecule index
!
         if (is_substring_in_string(string, 'mol')) then
            cursor = set_cursor_to_substring(string, 'mol')
         else
            call output%error_msg('could not find mol in MM geometry')
         endif
!
         string = adjustl(string(cursor+1:200))
!
         if (string(1:1) .ne. '=') call output%error_msg('in mol specification in MM geometry')
         string = adjustl(string(2:200))

         cursor = first_instance_of_character(string, ']')
         read(string(1:cursor-1),*) current_molecule
         string = adjustl(string(cursor+1:200))
!
!        Determine position
!
         cursor = first_instance_of_character(string, '[')
!
         coordinate = string(1:cursor-1)
         read(coordinate, *) positions(:, current_atom)
!
         string = adjustl(string(cursor+1:200))
!
!        Determine charge
!
         cursor = first_instance_of_character(string, 'q')
         string = adjustl(string(cursor+1:200))
!
         if (string(1:1) .ne. '=') call output%error_msg('in charge specification in MM geometry')
         string = adjustl(string(2:200))
!
         cursor = first_instance_of_character(string,']')
         read(string(1:cursor-1), * ) charge(current_atom)
!
         if(abs(charge(current_atom)).lt.1.0d-8) then
!
            call output%warning_msg('You put zero charge on atom = (i0)', &
                                     ints=[current_atom], fs='(t6,a)')
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
   end subroutine get_mm_geometry_non_polarizable_input_file
!
!
   function is_string_in_cs_list_input_file(this, keyword, section, string) result(in_list)
!!
!!    Is string in comma separated list?
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
      class(input_file), intent(in) :: this
!
      character(len=*), intent(in) :: section
      character(len=*), intent(in) :: keyword
      character(len=*), intent(in) :: string
!
      logical :: in_list
!
      character(len=200) :: line
      character(len=200) :: current_snippet
!
      integer :: i, current_cursor
!
      in_list = .false.
!
      call this%get_required_keyword(keyword, section, line)
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
                  in_list = .true.
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
               in_list = .true.
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
   subroutine read_keywords_and_geometry_input_file(this)
!!
!!    Read keywords and geometry
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
      class(input_file), intent(inout) :: this
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
      rewind(this%unit_)
!
      this%n_keyword_lines = 0
!
      do
!
         read(this%unit_, '(a)', iostat=io_error) line
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
         this%n_keyword_lines = this%n_keyword_lines + 1
!
      enddo
!
!     Continue reading the geometry to count QM and MM atoms
!
      QM_atoms = .true.
      this%n_qm_atom_lines = 0
      this%n_mm_atom_lines = 0
!
      do
!
         read(this%unit_, '(a)', iostat=io_error) line
!
         if (trim(adjustl(line)) == 'end geometry') then
!
            exit
!
         elseif (io_error .ne. 0) then
!
            call output%error_msg("The 'geometry' section appears to have no end in the input file.")
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
            this%n_qm_atom_lines = this%n_qm_atom_lines + 1
!
         else
!
            this%n_mm_atom_lines = this%n_mm_atom_lines + 1
!
         endif
!
      enddo
!
!     Place keyword lines in memory
!
      allocate(this%input_(this%n_keyword_lines))
!
      rewind(this%unit_)
!
      do i = 1, this%n_keyword_lines
!
         read(this%unit_, '(a)', iostat=io_error) line
!
         line = adjustl(line)
         call convert_to_lowercase(line)
!
         this%input_(i) = line
!
      enddo
!
!     Place geometries in memory
!
      if (this%n_qm_atom_lines < 1) call output%error_msg('No QM atoms in geometry!')
!
      allocate(this%geometry(this%n_qm_atom_lines))
!
      read(this%unit_, '(a)', iostat=io_error) line ! Reads 'geometry' line
!
!     QM geometry
!
      do i = 1, this%n_qm_atom_lines
!
         read(this%unit_, '(a)', iostat=io_error) line
!
         line = adjustl(line)
         call convert_to_lowercase(line)
!
         this%geometry(i) = line
!
      enddo
!
!     MM geometry if present
!
      if (this%n_mm_atom_lines > 0) then
!
         read(this%unit_, '(a)', iostat=io_error) line ! Reads '--' line
!
         allocate(this%mm_geometry(this%n_mm_atom_lines))
!
         do i = 1, this%n_mm_atom_lines
!
            read(this%unit_, '(a)', iostat=io_error) line
!
            line = adjustl(line)
            call convert_to_lowercase(line)
!
            this%mm_geometry(i) = line
!
         enddo
      endif
!
      call this%check_for_errors() ! Check for incorrect/missing keywords/sections
!
   end subroutine read_keywords_and_geometry_input_file
!
!
   subroutine get_section_limits(this, section, start_, end_)
!!
!!    Get section limits
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 / Mar 2019
!!
!!    Finds the start and the end index of the section given by the string 'section'
!!
      implicit none
!
      class(input_file), intent(in) :: this
!
      character(len=*), intent(in) :: section
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
      do i = 1, this%n_keyword_lines
!
         if (trim(this%input_(i)) == 'end ' // section) then
!
            n_ends = n_ends + 1
            end_ = i
!
         endif

         if (trim(this%input_(i)) == section) then
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
            // section // '" with more than one starting clause.')
!
      if (n_ends > 1)                                       &
         call output%error_msg('Tried to move to section "' &
            // section // '" with more than one ending clause.')
!
      if (n_ends == 0 .and. n_beginnings == 0)                             &
         call output%error_msg('Tried to move to non-existent section "'   &
            // section // '".')
!
      if (n_ends < 1) &
         call output%error_msg('Tried to move to section "' // section // '" with no end.')
!
      if (n_beginnings < 1) &
         call output%error_msg('Tried to move to section "' // section // '" with no beginning.')
!
   end subroutine get_section_limits
!
!
   subroutine cleanup_geometry_input_file(this)
!!
!!    Cleanup geometry
!!    Written by Sarai D. Folkestad, Mar 2020
!!
!!    Deallocates the geometry arrays
!!
      implicit none
!
      class(input_file), intent(inout) :: this
!
      deallocate(this%geometry)
!
      if (allocated(this%mm_geometry)) deallocate(this%mm_geometry)
!
   end subroutine cleanup_geometry_input_file
!
!
   subroutine cleanup_keywords_input_file(this)
!!
!!    Cleanup keywords
!!    Written by Sarai D. Folkestad, Mar 2020
!!
!!    Deallocates the keyword arrays
!!
      implicit none
!
      class(input_file), intent(inout) :: this
!
      deallocate(this%input_)
!
   end subroutine cleanup_keywords_input_file
!
!
   pure function is_embedding_on_input_file(this) result(embedding)
!!
!!    Is embedding on?
!!    Written by Sarai D. Folkestad, Sep 2020
!!
!!    Returns true if wavefunction is embedded.
!!
      implicit none
!
      class(input_file), intent(in) :: this
!
      logical :: embedding
!
      embedding = this%is_section_present('molecular mechanics') .or. &
                  this%is_section_present('pcm')
!
   end function is_embedding_on_input_file
!
!
   subroutine place_records_in_memory_input_file(this, section, records_in_memory)
!!
!!    Place records in memory
!!    Written by Sarai D. Folkestad, 2021
!!
!!    Checks if storage for 'section' is in memory or on disk
!!
!!    records_in_mem is intent(inout) because it should be set
!!    to a default value before a call to this routine
!!
      implicit none
!
      class(input_file), intent(in) :: this
      character(len=*),  intent(in) :: section
!
      logical, intent(inout) :: records_in_memory
!
      character(len=200) :: storage
!
      if (this%is_keyword_present('storage', trim(section))) then
!
         call this%get_keyword('storage', trim(section), storage)
!
!        Determine whether to store records in memory or on file
!
         if (trim(storage) == 'memory') then
!
            records_in_memory = .true.
!
         elseif (trim(storage) == 'disk') then
!
            records_in_memory = .false.
!
         else
!
            call output%error_msg('Could not recognize keyword storage in solver: ' // &
                                    trim(storage))
!
         endif
!
      endif
!
   end subroutine place_records_in_memory_input_file
!
!
   function keyword_is_allowed(this, keyword, section) result(allowed)
!!
!!    Keyword is allowed
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2019-2021
!!
      implicit none

      class(input_file), intent(in) :: this
!
      character(len=*), intent(in) :: keyword
      character(len=*), intent(in) :: section
!
      logical :: allowed
!
      integer :: i
!
      allowed = .false.
!
      do i = 1, size(this%sections)
!
         if (this%sections(i)%name_ == trim(section)) then
!
            allowed = any(this%sections(i)%keywords == trim(keyword))
            return
!
         endif
      enddo
!
   end function keyword_is_allowed
!
!
   function section_is_allowed(this, section) result(allowed)
!!
!!    Section is allowed
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2019-2021
!!
      implicit none

      class(input_file), intent(in) :: this
!
      character(len=*), intent(in) :: section
!
      logical :: allowed
!
      integer :: k
!
      allowed = .false.
!
      do k = 1, size(this%sections)
!
         if (this%sections(k)%name_ == trim(section)) allowed = .true.
!
      enddo
!
   end function section_is_allowed
!
!
   subroutine check_keyword_and_section(this, keyword, section)
!!
!!    Check keyword and section
!!    Written by Sarai D. Folkestad, Jul 2021
!!
      implicit none

      class(input_file), intent(in) :: this
!
      character(len=*), intent(in) :: keyword
      character(len=*), intent(in) :: section
!
      if (.not. this%section_is_allowed(section)) &
         call output%error_msg('requested keyword from non-existing section (a0)', &
                               chars=[trim(section)])
!
      if (.not. this%keyword_is_allowed(keyword, section)) &
         call output%error_msg('requested illegal keyword (a0)', chars=[trim(keyword)])
!
   end subroutine check_keyword_and_section
!
!
   function get_n_sets_and_ranges(string) result(n_elements)
!!
!!    Get n sets and ranges in string
!!    Written by Sarai D. Folkstad and Eirik F. Kjønstad, Mar 2019
!!
!!    Gets the number of sets/ranges in a string,
!!    To be used for reading of input.
!!
!!    Ranges should always be given as [a,b].
!!
!!    Lists/ranges should always be given as {a, b, c, d},
!!    that is, in set notation.
!!
      implicit none
!
      character(len=200), intent(inout) :: string
!
      integer :: n_elements
!
      integer :: i
      integer :: n_set_start, n_set_end
      integer :: n_range_start, n_range_end
!
      n_set_start = 0
      n_set_end   = 0
!
      n_range_start = 0
      n_range_end   = 0
!
      string = adjustl(string)
!
      do i = 1, len_trim(string)
!
         if (string(i:i) == '{') n_set_start = n_set_start + 1
         if (string(i:i) == '}') n_set_end = n_set_end + 1
!
         if (string(i:i) == '[') n_range_start = n_range_start + 1
         if (string(i:i) == ']') n_range_end = n_range_end + 1
!
      enddo
!
      if (n_set_start /= n_set_end) call output%error_msg('found open set in input file')
      if (n_range_start /= n_range_end) call output%error_msg('found open range in input file')
!
      n_elements = n_set_start + n_range_start
!
   end function get_n_sets_and_ranges
!
!
   function get_n_elements_in_string(string) result(n_elements)
!!
!!    Get n elements in string
!!    Written by Sarai D. Folkstad and Eirik F. Kjønstad, Mar 2019
!!
!!    Gets the number of elements in range or list,
!!    To be used for reading of input.
!!
!!    Ranges should always be given as [a,b].
!!
!!    Lists should always be given as {a, b, c, d},
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
         n_elements = 1 ! Assuming that the set contains at least one element (otherwize why give list?)
!
!        Loop through and count commas
!
         do i = 2, n_characters - 1
!
            if (string(i:i) == ',') n_elements = n_elements + 1
!
         enddo
!
      else ! Did not find list or
!
         n_elements = 0
!
      endif
!
   end function get_n_elements_in_string
!
!
   subroutine get_elements_in_string(string, n_elements, elements)
!!
!!    Get elements
!!    Written by Sarai D. Folkstad and Eirik F. Kjønstad, Mar 2019
!!
!!    Gets the elements from range or list.
!!    To be used for reading of input.
!!
!!    Ranges should always be given as [a,b].
!!
!!    Lists should always be given as {a, b, c, d},
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
         first            = 2
         n_elements_found = 0
!
         do j = 1, n_elements
!
            do i = first, n_characters - 1
!
               if (string(i:i) == ',') exit
!
            enddo
!
            read(string(first:i-1), *) elements(j)
!
            n_elements_found = n_elements_found + 1
!
            first = i + 1
!
            if (first == n_characters) exit
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
   end subroutine get_elements_in_string
!
!
   subroutine get_reals_in_string(string, n_elements, elements)
!!
!!    Get reals
!!    Written by Sarai D. Folkstad and Eirik F. Kjønstad, Mar 2019
!!    Modified by Andreas Skeidsvoll, Sep 2019: Reads reals instead of integers
!!
!!    Gets the reals from list.
!!    To be used for reading of input.
!!
!!    Lists should always be given as {a, b, c, d},
!!    that is, in set notation.
!!
      implicit none
!
      character(len=200), intent(inout) :: string
!
      integer, intent(in) :: n_elements
!
      real(dp), dimension(n_elements), intent(out) :: elements
!
!     Local variables
!
      integer :: first, n_characters, n_elements_found
      integer :: i, j
!
      string = adjustl(string)
!
      n_characters = len_trim(string)
!
      if (string(1:1)=='{') then ! list given
!
!        Sanity check - Is set closed?
!
         if (string(n_characters:n_characters) /= '}') call output%error_msg('found open set in input file.')
!
!        Loop through and set the elements
!
         first            = 2
         n_elements_found = 0
!
         do j = 1, n_elements
!
            do i = first, n_characters - 1
!
               if (string(i:i) == ',') exit
!
            enddo
!
            read(string(first:i-1), *) elements(j)
!
            n_elements_found = n_elements_found + 1
!
            first = i + 1
!
            if (first == n_characters) exit
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
   end subroutine get_reals_in_string
!
!
   function get_n_state_guesses(this) result(n_guesses)
!!
!!
      implicit none
!
      class(input_file) :: this
!
      integer :: n_guesses
!
      character(len=200) :: state_guesses
!
      call this%get_keyword('state guesses', 'solver cc es', state_guesses)
!
      n_guesses = get_n_sets_and_ranges(state_guesses)
!
   end function get_n_state_guesses
!
!
   subroutine get_state_guesses(this, occupied, virtual, n_states)
!!
!!    Get state guesses
!!    Written by Sarai D. Folkestad, Dec 2021
!!
!!    state guesses: {i=1, a=2}, {i=2, a=2}, {i=1, a=1}
!!
!!    Provided as a (comma separated) list of sets.
!!
!!    A set comprises of two elements which are
!!    the occupied and virtual indices of the start vector
!!
!!    'i =' provides occupied index
!!    'a =' provides virtual index
!!
!
      implicit none
!
      class(input_file) :: this
!
      integer, intent(in) :: n_states
      integer, dimension(n_states), intent(out) :: virtual, occupied
!
      character(len=200) :: state_guesses
      character(len=200), dimension(:), allocatable :: split_string
!
      integer :: n_elements, cursor_a, cursor_i, i, n_occupied, n_virtual
!
      call this%get_keyword('state guesses', 'solver cc es', state_guesses)
      state_guesses = adjustl(state_guesses)
!
      n_occupied = n_instances_of_character(state_guesses, 'i')
      n_virtual = n_instances_of_character(state_guesses, 'a')
!
      if (n_virtual /= n_occupied) &
         call output%error_msg('in keyword (state guesses)')
!
      if (n_virtual /= n_states) &
         call output%error_msg('in keyword (state guesses)')
!
!     {i=1, a=2}, {i=2, a=2}, {i=1, a=1} -> i=1, a=2} i=2, a=2} i=1, a=1}
      call remove_delimiter_from_string(state_guesses, '{')
!
      n_elements = n_instances_of_character(state_guesses, '}')
!
      if (n_elements /= n_states) &
         call output%error_msg('in keyword (state guesses)')
!
      allocate(split_string(n_elements))
!
!     i=1, a=2} i=2, a=2} i=1, a=1} -> string array ("i=1, a=2", "i=2, a=2", "i=1, a=1")
      call split_at_delimiter(state_guesses, n_elements, split_string, '}')
!
      do i = 1, n_elements
!
         call remove_delimiter_from_string(split_string(i), '=')
         call remove_delimiter_from_string(split_string(i), ',')
!
         cursor_i = set_cursor_to_substring(split_string(i), 'i')
         cursor_a = set_cursor_to_substring(split_string(i), 'a')
!
         if (cursor_i < cursor_a) then
!
            read(split_string(i)(cursor_i + 1 : cursor_a -1), *) occupied(i)
            read(split_string(i)(cursor_a + 1:), *) virtual(i)
!
         else
!
            read(split_string(i)(cursor_a + 1 : cursor_i -1), *) virtual(i)
            read(split_string(i)(cursor_i + 1:), *) occupied(i)
!
         endif
!
      enddo
!
      deallocate(split_string)
!
   end subroutine get_state_guesses
!
!
!
end module input_file_class
