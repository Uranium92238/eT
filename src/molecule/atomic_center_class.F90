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
module atomic_center_class
!
!!
!!    Atomic center class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018 and 2020
!!
!
   use parameters
   use global_out, only : output
   use sequential_file_class, only : sequential_file
   use shell_class, only : shell
   use memory_manager_class, only : mem
!
   implicit none
!
!  NIST Atomic Spectra Database (ver. 5.8),
!  https://physics.nist.gov/asd, 2021, January 27
!
!  Multiplicity for lead from:
!  Dembcyński et al. Phys. Rev. A, 49, 745-754, 1994
!
   integer, dimension(105) :: atomic_multiplicities = (/ &
            2,                                                            1, & ! H,    He
            2,1,                                                2,3,4,3,2,1, & ! Li .. Ne
            2,1,                                                2,3,4,3,2,1, & ! Na .. Ar
            2,1,2,                            3,4,7,6,5,4,3,2,1,2,3,4,3,2,1, & ! K  .. Kr
            2,1,2,                            3,6,7,6,5,4,1,2,1,2,3,4,3,2,1, & ! Rb .. Xe
            2,1,2,2,1,4,5,6,7,8,9,6,5,4,3,2,1,3,4,5,6,5,4,3,2,1,2,3,4,3,2,1, & ! Cs .. Rn
            2,1,2,2,3,4,5,6,7,8,9,6,5,4,3,2,1,3,4 /)                           ! Fr .. Db
!
   type :: atomic_center
!
      character(len=2) :: symbol ! He, Fe, P, ...
!
      integer :: n_ao
!
      integer :: n_shells
      type(shell), dimension(:), allocatable :: shells
!
      character(len=100) :: basis
!
      real(dp), dimension(3) :: coordinates ! x, y, z
!
      integer :: number_        ! Atomic number Z
      integer :: charge
      integer :: nuclear_charge ! = number_ (if (ghost) then = 0)
!
      integer :: input_number  ! Atom # in the input file
      integer :: libint_number ! Center # in Libint
!
      logical :: cartesian ! If not, it is spherical
!
   contains
!
      procedure :: initialize_shells &
                => initialize_shells_atomic_center
!
      procedure :: symbol_to_number &
                => symbol_to_number_atomic_center
!
      procedure :: rename_core_valence_dunning_sets &
                => rename_core_valence_dunning_sets_atomic_center
!
      procedure :: get_basis_set_name &
                => get_basis_set_name_atomic_center
!
      procedure :: set_cartesian &
                => set_cartesian_atomic_center
!
      procedure :: read_atomic_uhf_density &
                => read_atomic_uhf_density_atomic_center
!
      procedure :: get_ground_state_multiplicity &
                => get_ground_state_multiplicity_atomic_center
!
      procedure :: get_ao_range &
                => get_ao_range_atomic_center
!
      procedure :: is_ghost &
                => is_ghost_atomic_center
!
      procedure :: get_molden_order &
                => get_molden_order_atomic_center
!
      procedure :: get_ao_normalization_factors &
                => get_ao_normalization_factors_atomic_center
!
      procedure :: get_identifier_string &
                => get_identifier_string_atomic_center
!
      procedure :: cleanup &
                => cleanup_atomic_center
!
      procedure, private :: read_basis_info & ! Reads primitives, exponents, coefficients
                         => read_basis_info_atomic_center
!
      procedure, private :: read_shell_basis_info &
                         => read_shell_basis_info_atomic_center
!
   end type atomic_center
!
   include "../libint/atom_init_cdef.F90"
!
   interface atomic_center
!
      procedure :: new_atomic_center
!
   end interface atomic_center
!
!
contains
!
!
   function new_atomic_center(libint_number,    &
                              input_number,     &
                              symbol,           &
                              coordinates,      &
                              basis,            &
                              basis_type_,      &
                              is_ghost) result(center)
!!
!!    New atomic center
!!    Written by Eirik F. Kjønstad, 2020
!!
      implicit none
!
      type(atomic_center) :: center
!
      integer, intent(in) :: libint_number, input_number
!
      character(len=2), intent(in) :: symbol
!
      real(dp), dimension(3), intent(in) :: coordinates
!
      character(len=100), intent(in) :: basis
!
      character(len=*), intent(in) :: basis_type_
!
      logical, intent(in) :: is_ghost
!
      center%libint_number = libint_number
      center%symbol        = symbol
      center%basis         = basis
      center%coordinates   = coordinates
      center%input_number  = input_number
!
      call center%symbol_to_number()
!
      center%charge = 0
      center%nuclear_charge = center%number_
!
      if (is_ghost) center%nuclear_charge = 0
!
      if (center%number_ .eq. -1) &
         call output%error_msg('illegal atomic symbol, check the eT.inp file ')
!
      call center%rename_core_valence_dunning_sets() ! Otherwise Libint looks for 'augmentation'
                                                     ! files that don't exist for these basis sets
!
      call center%set_cartesian(basis_type_)
!
   end function new_atomic_center
!
!
   subroutine rename_core_valence_dunning_sets_atomic_center(center)
!!
!!    Rename core valence Dunning sets
!!    Written by Eirik F. Kjønstad, Nov 2019
!!
!!    Renames the basis sets (if any) of the core-valence type:
!!
!!       aug-cc-pCVXZ -> _aug-cc-pCVXZ
!!
!!    This is necessary because the Libint files must be named
!!    with an "_" prefix to avoid it looking for an "augmentation" file.
!!
      implicit none
!
      class(atomic_center) :: center
!
      integer :: k
!
      integer, parameter :: n_renamings = 4
!
      character(len=20), dimension(n_renamings), parameter :: original_names = ['aug-cc-pcvdz', &
                                                                                'aug-cc-pcvtz', &
                                                                                'aug-cc-pcvqz', &
                                                                                'aug-cc-pcv5z']
!
      character(len=20), dimension(n_renamings), parameter :: new_names = ['_aug-cc-pcvdz', &
                                                                           '_aug-cc-pcvtz', &
                                                                           '_aug-cc-pcvqz', &
                                                                           '_aug-cc-pcv5z']
!
      do k = 1, n_renamings
!
         if (trim(center%basis) == trim(original_names(k))) then
!
            center%basis = trim(new_names(k))
!
         endif
!
      enddo
!
   end subroutine rename_core_valence_dunning_sets_atomic_center
!
!
   pure function get_basis_set_name_atomic_center(center) result(name_)
!!
!!    Get basis set name
!!    Written by Eirik F. Kjønstad, 2020
!!
      implicit none
!
      class(atomic_center), intent(in) :: center
!
      character(len=100) :: name_
!
      integer :: I
!
      I = 1
      if (trim(center%basis(1:1)) .eq. '_') I = 2
!
      name_ = center%basis(I:)
!
   end function get_basis_set_name_atomic_center
!
!
   pure subroutine symbol_to_number_atomic_center(center)
!!
!!    Symbol to atomic number
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!!    Uses the periodic table to determine atomic number from atomic symbol
!!
      use periodic_table
!
      implicit none
!
      class(atomic_center), intent(inout) :: center
!
      integer :: i
!
      center%number_ = -1
!
      do i = 1, size_periodic_table
!
         if (atomic_symbol(i) == center%symbol) then
!
            center%number_ = i
            return
!
         endif
      enddo
!
   end subroutine symbol_to_number_atomic_center
!
!
   subroutine initialize_shells_atomic_center(center)
!!
!!    Initialize shells
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!!    Allocates the shell array, determines the range information of each shell
!!    (first, length), as well as basis information for each shell (primitives,
!!    exponents, and coefficients).
!!
      use iso_c_binding, only: c_int
!
      implicit none
!
      class(atomic_center) :: center
!
      integer(c_int), dimension(:), allocatable :: n_aos_in_shell
      integer(c_int), dimension(:), allocatable :: first_ao_in_shell
      integer(c_int), dimension(:), allocatable :: shell_numbers
!
      integer :: j
!
      integer(c_int) :: center_number_c, n_shells_c
!
!     Set number of shells and allocate shell array
!
      center_number_c = int(center%libint_number, kind=c_int)
!
      call get_n_shells_on_atom_c(center_number_c, n_shells_c)
!
      center%n_shells = int(n_shells_c)
!
      allocate(center%shells(center%n_shells))
!
!     Get shell ranges and determine the total number of AOs on the center
!
      allocate(n_aos_in_shell(center%n_shells))
      allocate(shell_numbers(center%n_shells))
      allocate(first_ao_in_shell(center%n_shells))
!
      call get_n_aos_in_shell_c(center_number_c, n_aos_in_shell)
      call get_shell_numbers_c(center_number_c, shell_numbers)
      call get_first_ao_in_shells_c(center_number_c, first_ao_in_shell)
!
      center%n_ao = 0
!
      do j = 1, center%n_shells
!
         center%shells(j) = shell(first=int(first_ao_in_shell(j)), &
                                  length=int(n_aos_in_shell(j)),   &
                                  number_=int(shell_numbers(j)),   &
                                  cartesian=center%cartesian)
!
         center%n_ao = center%n_ao + center%shells(j)%length
!
      enddo
!
      deallocate(n_aos_in_shell)
      deallocate(shell_numbers)
      deallocate(first_ao_in_shell)
!
!     Read additional shell information from basis set file(s)
!     (number of primitives, exponents, coefficients)
!
      call center%read_basis_info()
!
   end subroutine initialize_shells_atomic_center
!
!
   subroutine read_atomic_uhf_density_atomic_center(center, atomic_D)
!!
!!    Read atomic density
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    Read the atomic density matrices from file and adds them together. By assumption, these
!!    densities are the result of atomic UHF calculations with the restriction that valence
!!    electrons are smeared out to ensure that the density is spherically symmetric (as required
!!    for a rotationally invariant SAD guess).
!!
      implicit none
!
      class(atomic_center), intent(in) :: center
!
      real(dp), dimension(center%n_ao, center%n_ao), intent(out) :: atomic_D
!
      real(dp), dimension(:,:), allocatable :: temporary
      character(len=100)                    :: alpha_fname, beta_fname, name_
!
      type(sequential_file), allocatable :: alpha_D_file
      type(sequential_file), allocatable :: beta_D_file
!
      atomic_D = zero
!
      name_ = "sad_" // trim(center%get_identifier_string())
!
      alpha_fname = trim(name_) // '_alpha'
      beta_fname  = trim(name_) // '_beta'
!
      alpha_D_file = sequential_file(trim(alpha_fname))
      beta_D_file  = sequential_file(trim(beta_fname))
!
      call alpha_D_file%open_('read', 'rewind')
      call beta_D_file%open_( 'read', 'rewind')
!
      call mem%alloc(temporary, center%n_ao, center%n_ao)
!
      call alpha_D_file%read_(temporary, center%n_ao**2)
      atomic_D = atomic_D + temporary
!
      call beta_D_file%read_(temporary, center%n_ao**2)
      atomic_D = atomic_D + temporary
!
      call mem%dealloc(temporary, center%n_ao, center%n_ao)
!
      call alpha_D_file%close_()
      call beta_D_file%close_()
!
   end subroutine read_atomic_uhf_density_atomic_center
!
!
   function get_ground_state_multiplicity_atomic_center(center) result(multiplicity)
!!
!!    Get multiplicity
!!    Written by Tor S. Haugland, 2019
!!
!!    Returns the multiplicity from the atomic number.
!!    Based on the table 5.1 in Griffith's "Introduction to Quantum Mechanics", 2. ed
!!
      implicit none
!
      class(atomic_center), intent(in) :: center
!
      integer :: multiplicity
!
      if (center%number_ > size(atomic_multiplicities)) &
         call output%error_msg("Atom multiplicity not supported for Z > (i0)", &
                                ints=[size(atomic_multiplicities)])
!
      multiplicity = atomic_multiplicities(center%number_)
!
   end function get_ground_state_multiplicity_atomic_center
!
!
  subroutine read_basis_info_atomic_center(center)
!!
!!    Read basis info
!!    Written by Sarai D. Folkestad, Dec 2018
!!
!!    Opens and reads the basis set file and calls routine read_shell_basis_info to set
!!    the information for each shell (n_primitives, exponents, and coefficients).
!!
!!    Moved and adapted to atomic center class, Eirik F. Kjønstad, 2020.
!!
      use string_utilities, only: convert_to_lowercase
!
      implicit none
!
      class(atomic_center), intent(inout) :: center
!
      character(len=200) :: basis
      character(len=200) :: libint_path
      character(len=200) :: filename
!
      type(sequential_file), allocatable :: basis_file
!
      integer :: sh
!
!     Initialize shell index to zero
!
      sh = 0
!
!     Do the shells in the non-augmented part first
!
      if (center%basis(1:3) .ne. 'aug') then
!
         write(basis,'(a)') center%basis
!
      else
!
         write(basis,'(a)') center%basis(5:100)
!
      endif
!
      call convert_to_lowercase(basis)
      call get_environment_variable("LIBINT_DATA_PATH", libint_path)
!
      filename = trim(libint_path) // '/' // trim(basis) // '.g94'
!
      basis_file = sequential_file(trim(filename), 'formatted')
!
      call basis_file%open_('read', 'rewind')
!
      call center%read_shell_basis_info(basis_file, sh)
!
      call basis_file%close_()
!
!     Then the shells in the augmented part
!
      if (center%basis(1:3) .eq. 'aug') then
!
         basis = 'augmentation-' // trim(basis)
!
         filename = trim(libint_path) // '/' // trim(basis) // '.g94'
!
         basis_file = sequential_file(trim(filename), 'formatted')
!
         call basis_file%open_('read', 'rewind')
!
         call center%read_shell_basis_info(basis_file, sh)
!
         call basis_file%close_()
!
      endif
!
   end subroutine read_basis_info_atomic_center
!
!
   subroutine read_shell_basis_info_atomic_center(center, basis_file, sh)
!!
!!    Read shell basis info
!!    Written by Sarai D. Folkestad, Dec 2018
!!
!!    Sets the number of primitives, and the coefficient and exponents of the primitives, in
!!    each of the shells of the atom.
!!
!!       basis_file: File associated with the basis set. Assumed to be opened
!!                   and rewinded when the routine is called.
!!
!!       sh:         Equal to zero on entry if there is only one file associated with the
!!                   basis set. For each read shell, this integer is incremented by one.
!!                   In the case of multiple basis set files, shell is first incremented
!!                   from zero when reading the first file and then further incremented
!!                   when reading subsequent (e.g., non-augmented and then augmented).
!!
!!    Modified by Marco Scavino, 2019
!!    Extended to read the SP primitives from Gaussian format
!!
!!    Moved and adapted to atomic center class, Eirik F. Kjønstad, 2020
!!
      implicit none
!
      class(atomic_center), intent(inout) :: center
!
      type(sequential_file), intent(in) :: basis_file
!
      integer, intent(inout) :: sh
!
      character(len=200) :: line
      logical            :: elm_found ! Element found?
!
!
      integer            :: n_primitive, primitive
      character(len=2)   :: ang_mom
      real(dp)           :: coefficient, coefficient_2, exponent_
!
!     Find coordinates of atom/element in the basis file
!     Look for:  ****
!     Symbol
!
      elm_found = .false.
!
      call basis_file%read_(line,'(a200)')
!
      do while (.not. elm_found)
!
         if (trim(line) == '****') then
!
            call basis_file%read_(line,'(a200)')
!
            if ((line(1:2)) == center%symbol) then
!
               elm_found = .true.
!
            endif
!
         endif
!
         if (.not. elm_found) then
!
            call basis_file%read_(line,'(a200)')
!
         endif
!
      enddo
!
!     Read angular momentum and number of primitives
!
      call basis_file%read_(line,'(a200)')
!
      do while (trim(line) .ne. '****')
!
!        Increment shell index and read angular momentum symbol and the number of primitive
!        basis functions in the shell
!
         sh = sh + 1
!
         read(line, *) ang_mom, n_primitive
!
!        Sanity check -> does shell exceed number of shells on atom
!
         if (sh .gt. center%n_shells) &
            call output%error_msg('Mismatch in number of shells in read_shell_basis_info')
!
!        Sanity check -> does shell have the correct angular momentum
!
         if (angular_momentum_from_symbol(ang_mom) .ne. center%shells(sh)%get_angular_momentum()) &
            call output%error_msg('Mismatch in angular momentum in read_shell_basis_info')
!
!        Set number of primitives and initialize exponents and coefficient array
!
         call center%shells(sh)%set_n_primitives(n_primitive)
!
         call center%shells(sh)%initialize_exponents()
         call center%shells(sh)%initialize_coefficients()
!
!        Loop over primitives and set coefficient and exponent
!
         if (ang_mom == "SP") then
!
!           In case of "SP" shell, split S and P coefficients
!
            call center%shells(sh+1)%set_n_primitives(n_primitive)
!
            call center%shells(sh+1)%initialize_exponents()
            call center%shells(sh+1)%initialize_coefficients()
!
            do primitive = 1, n_primitive
!
!              Read coefficients and exponent
!
               call basis_file%read_(line,'(a200)')
               read(line, *) exponent_, coefficient, coefficient_2
!
               call center%shells(sh)%set_exponent_i(primitive, exponent_)
               call center%shells(sh+1)%set_exponent_i(primitive, exponent_)
!
               call center%shells(sh)%set_coefficient_i(primitive, coefficient)
               call center%shells(sh+1)%set_coefficient_i(primitive, coefficient_2)
!
            enddo
!
!           Since we have now read two shells, we need to increment shell again
!
            sh = sh + 1
!
         else
!
            do primitive = 1, n_primitive
!
!              Read coefficient and exponent
!
               call basis_file%read_(line,'(a200)')
               read(line, *) exponent_, coefficient
!
               call center%shells(sh)%set_exponent_i(primitive, exponent_)
               call center%shells(sh)%set_coefficient_i(primitive, coefficient)

!
            enddo
!
         end if
!
         call basis_file%read_(line,'(a200)')
!
      enddo
!
   end subroutine read_shell_basis_info_atomic_center
!
!
   function angular_momentum_from_symbol(symbol) result(l)
!!
!!    Angular momentum from symbol
!!    Written by Sarai D. Folkestad, Dec 2018
!!
!!    From the symbol ('S', 'P', 'D', ...), the function returns the associated
!!    angular momentum (0, 1, 2, ...).
!!
      implicit none
!
      character(len=1), intent(in) :: symbol
!
      integer :: l
!
      l = 0
!
      if (symbol == 'S') then
!
         l = 0
!
      elseif (symbol == 'P') then
!
         l = 1
!
      elseif (symbol == 'D') then
!
         l = 2
!
      elseif (symbol == 'F') then
!
         l = 3
!
      elseif (symbol == 'G') then
!
         l = 4
!
      elseif (symbol == 'H') then
!
         l = 5
!
      elseif (symbol == 'I') then
!
         l = 6
!
      else
!
         call output%error_msg('no support for orbitals with l > 6')
!
      endif
!
   end function angular_momentum_from_symbol
!
!
   subroutine set_cartesian_atomic_center(center, basis_type)
!!
!!    Set cartesian
!!    Written by Eirik F. Kjønstad, 2020
!!
!!    Includes previous routine 'default_cartesian_basis' written by Rolf H. Myhre, 2020,
!!    for setting default basis.
!!
!!    Determines whether the basis is Cartesian or spherical.
!!
!!       basis_type: string specifying whether basis should be 'cartesian', 'spherical',
!!                   or 'default' (see list in if statement below)
!!
      implicit none
!
      class(atomic_center), intent(inout) :: center
!
      character(len=*), intent(in) :: basis_type
!
      if (trim(basis_type) == 'cartesian') then
!
         center%cartesian = .true.
!
      elseif (trim(basis_type) == 'spherical') then
!
         center%cartesian = .false.
!
      elseif (trim(basis_type) == 'default') then
!
         center%cartesian = .false.
!
         if (trim(center%basis) .eq. '3-21g' .or.         &
             trim(center%basis) .eq. '6-31g' .or.         &
             trim(center%basis) .eq. '6-31+g' .or.        &
             trim(center%basis) .eq. '6-31++g' .or.       &
             trim(center%basis) .eq. '6-31g*' .or.        &
             trim(center%basis) .eq. '6-31g**' .or.       &
             trim(center%basis) .eq. '6-31+g*' .or.       &
             trim(center%basis) .eq. '6-31+g**' .or.      &
             trim(center%basis) .eq. '6-31++g**' .or.     &
             trim(center%basis) .eq. '6-31g(d,p)' .or.    &
             trim(center%basis) .eq. '6-31g(2df,p)' .or.  &
             trim(center%basis) .eq. '6-31g(3df,3pd)') then
!
            center%cartesian = .true.
!
         endif
!
      endif
!
   end subroutine set_cartesian_atomic_center
!
!
   function get_ao_range_atomic_center(center) result(aos)
!!
!!    Get AO range
!!    Written by Eirik F. Kjønstad, 2020
!!
!!    Returns an range containing the first and last AOs on the center,
!!    as well as the length (in an range object 'aos')
!!
      use range_class
!
      implicit none
!
      class(atomic_center), intent(in) :: center
!
      type(range_) :: aos
!
      integer :: first, length
!
      first  = center%shells(1)%first
      length = center%shells(center%n_shells)%get_last() - first + 1
!
      aos = range_(first, length)
!
   end function get_ao_range_atomic_center
!
!
   pure function is_ghost_atomic_center(center) result(is_ghost)
!!
!!    Is ghost?
!!    Written by Tor S. Hauglan, May 2021
!!
!!    Ghosts are atoms with zero charge. Their only function is
!!    enlargening the basis.
!!
      implicit none
!
      class(atomic_center), intent(in) :: center
      logical :: is_ghost
!
      is_ghost = (center%nuclear_charge == 0)
!
   end function is_ghost_atomic_center
!
!
   function get_molden_order_atomic_center(center) result(map)
!!
!!    Get Molden order
!!    Written by Alexander C. Paul, May 2021
!!
!!    Return index list mapping AOs to the order molden expects
!!
      implicit none
!
      class(atomic_center), intent(in) :: center
!
      integer, dimension(center%n_ao) :: map
!
      integer :: s, offset
!
      offset = 1
!
      do s = 1, center%n_shells
!
         map(offset : offset+center%shells(s)%length-1) = &
            center%shells(s)%get_molden_order()
!
         offset = offset + center%shells(s)%length
!
      end do
!
   end function get_molden_order_atomic_center
!
!
   function get_ao_normalization_factors_atomic_center(center) result(factors)
!!
!!    Get AO normalization factors
!!    Written by Alexander C. Paul, May 2021
!!
      implicit none
!
      class(atomic_center), intent(in) :: center
!
      real(dp), dimension(center%n_ao) :: factors
!
      integer :: s, offset
!
      offset = 1
!
      do s = 1, center%n_shells
!
         factors(offset : offset+center%shells(s)%length-1) = &
            center%shells(s)%get_normalization_factor()
!
         offset = offset + center%shells(s)%length
!
      end do
!
   end function get_ao_normalization_factors_atomic_center
!
!
   pure function get_identifier_string_atomic_center(center) result(identifier)
!!
!!    Get identifier string
!!    Written by Sarai D. Folkestad, Jun 2021
!!
!!    Returns a string with symbol and basis set
!!
      implicit none
!
      class(atomic_center), intent(in) :: center
      character(len=50) :: identifier
!
      identifier = trim(center%symbol) // "_" // trim(center%basis)
!
   end function get_identifier_string_atomic_center
!
!
   subroutine cleanup_atomic_center(center)
!!
!!    Cleanup
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
      implicit none
!
      class(atomic_center) :: center
!
      integer :: sh
!
      do sh = 1, center%n_shells
!
         call center%shells(sh)%cleanup()
!
      enddo
!
   end subroutine cleanup_atomic_center
!
!
end module atomic_center_class
