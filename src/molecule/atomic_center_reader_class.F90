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
module atomic_center_reader_class
!
!!
!!    Atomic center reader
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2020-2021
!!
!!    Class responsible for reading the atomic centers from the input file.
!!
!!    The atomic center reader processes the information in the input file and initializes
!!    a set of atomic center objects ('centers') and an array of center subsets that specifies the
!!    active atom spaces ('nonempty sets').
!!
!!    The object processes information via the 'read_centers' subroutine.
!!
!!    The user can then get the centers and sets via get-routines.
!!
!
   use range_class,           only: range_
   use named_range_class,     only: named_range
   use global_in,             only: input
   use global_out,            only: output
   use kinds,                 only: dp
   use memory_manager_class,  only: mem
   use atomic_center_class,   only: atomic_center
!
   implicit none
!
   type :: atomic_center_reader
!
!     Atomic center objects and center subsets
      integer, private :: n_centers
      type(atomic_center), dimension(:), allocatable, private :: centers
!
      integer, private :: n_sets, n_nonempty_sets
      type(named_range), dimension(:), allocatable, private :: sets, nonempty_sets
!
!     ordered_wfs = wf-order used to order the center subsets
      character(len=200), dimension(6), private :: ordered_wfs
!
!     type = standard, spherical, or Gaussian
      character(len=200), private :: basis_type_
!
   contains
!
      procedure, public :: get_centers &
                        => get_centers_atomic_center_reader
!
      procedure, public :: get_n_centers &
                        => get_n_centers_atomic_center_reader
!
      procedure, public :: get_center_subsets &
                        => get_center_subsets_atomic_center_reader
!
      procedure, public :: get_n_center_subsets &
                        => get_n_center_subsets_atomic_center_reader
!
      procedure, public :: read_centers &
                        => read_centers_atomic_center_reader
!
      procedure, private :: create_centers                ! Sets up the array of atomic centers
      procedure, private :: order_center_indices_by_basis
      procedure, private :: create_nonempty_subsets       ! Removes empty atomic subsets
!
!     The routines below are used to set up atomic center subsets. The information extracted is the
!     name of each set ('unclassified', 'hf', ...) and the range information (first, last, length)
!
      procedure, private :: determine_n_subsets
      procedure, private :: determine_subset_names
      procedure, private :: determine_subset_bases
!
!     determine_subset_indices_and_ranges determines the atomic center indices
!     for each of the sets (subset-indices) and stores how many centers there
!     are in each set (subset-length). The indices are needed to reorder the centers
!     so that each set is contiguous (see create_centers).
!
!     How this is done depends on how the active spaces are specified: no active space
!     (unclassified), specified as {1,3,5,...} or [1-4, 6-8] (list/range), or as a central atom
!     with a surrounding radius (central_atom)
!
      procedure, private :: determine_subset_indices_and_ranges
!
      procedure, private :: determine_subset_indices_and_ranges_list_range
      procedure, private :: determine_subset_indices_and_ranges_central_atom
      procedure, private :: determine_subset_indices_and_ranges_unclassified
!
   end type atomic_center_reader
!
!
   interface atomic_center_reader
!
      procedure :: new_atomic_center_reader
!
   end interface atomic_center_reader
!
!
contains
!
!
   function new_atomic_center_reader(basis_type_) result(this)
!!
!!    New atomic center reader
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none
!
      type(atomic_center_reader) :: this
!
      character(len=200), intent(in) :: basis_type_ ! standard, spherical, or Gaussian
!
      this%basis_type_ = basis_type_
!
      this%ordered_wfs = [character(len=200) :: 'cc3',      &
                                                'ccsd(t)',  &
                                                'ccsd',     &
                                                'cc2',      &
                                                'ccs',      &
                                                'hf']
!
   end function new_atomic_center_reader
!
!
   subroutine get_centers_atomic_center_reader(this, centers)
!!
!!    Get centers
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none
!
      class(atomic_center_reader), intent(in) :: this
!
      type(atomic_center), dimension(this%n_centers), intent(out) :: centers
!
      centers = this%centers
!
   end subroutine get_centers_atomic_center_reader
!
!
   pure function get_n_centers_atomic_center_reader(this) result(n_centers)
!!
!!    Get number of centers
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none
!
      class(atomic_center_reader), intent(in) :: this
!
      integer :: n_centers
!
      n_centers = this%n_centers
!
   end function get_n_centers_atomic_center_reader
!
!
   pure subroutine get_center_subsets_atomic_center_reader(this, center_subsets)
!!
!!    Get center subsets
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none
!
      class(atomic_center_reader), intent(in) :: this
!
      type(named_range), dimension(this%n_nonempty_sets), intent(out) :: center_subsets
!
      center_subsets = this%nonempty_sets
!
   end subroutine get_center_subsets_atomic_center_reader
!
!
   pure function get_n_center_subsets_atomic_center_reader(this) result(n_center_subsets)
!!
!!    Get number of center subsets
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none
!
      class(atomic_center_reader), intent(in) :: this
!
      integer :: n_center_subsets
!
      n_center_subsets = this%n_nonempty_sets
!
   end function get_n_center_subsets_atomic_center_reader
!
!
   subroutine read_centers_atomic_center_reader(this)
!!
!!    Read centers
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2020
!!
!!    Reads the basis centers from input and initializes 'centers' objects as well as
!!    'center subsets' (which specifies the active atomic-center subspaces).
!!
!!    Based on routines with contributions from Sarai D. Folkestad, Eirik F. Kjønstad,
!!    Rolf H. Myhre, and Åsmund H. Tveten.
!!
      use parameters, only: bohr_to_angstrom
!
      implicit none
!
      class(atomic_center_reader), intent(inout) :: this
!
      character(len=2), dimension(:), allocatable   :: center_symbols
      character(len=100), dimension(:), allocatable :: center_bases
      real(dp), dimension(:,:), allocatable         :: center_positions
      integer, dimension(:), allocatable            :: center_indices
      logical, dimension(:), allocatable            :: center_is_ghosts
      integer, dimension(:), allocatable            :: center_charges
!
      logical :: units_angstrom
!
!     Read the number of centers, as well as associated symbols, positions and bases
!
      this%n_centers = input%get_n_atoms()
!
      call mem%alloc(center_positions, 3, this%n_centers)
      call mem%alloc(center_charges, this%n_centers)
!
      allocate(center_symbols(this%n_centers))
      allocate(center_bases(this%n_centers))
      allocate(center_is_ghosts(this%n_centers))
!
      call input%get_geometry(this%n_centers,   &
                              center_symbols,   &
                              center_positions, &
                              center_bases,     &
                              center_charges,   &
                              units_angstrom,   &
                              center_is_ghosts)
!
      if (.not. units_angstrom) then
!
         call dscal(3*this%n_centers,  &
                    bohr_to_angstrom,  &
                    center_positions,  &
                    1)
!
      endif
!
!     Calculate number of subsets and set their names
!
      call this%determine_n_subsets()
!
      allocate(this%sets(this%n_sets))
!
      call this%determine_subset_names()
!
!     Store center-indices and range (first, last, length) for each subset
!
      call mem%alloc(center_indices, this%n_centers)
!
      call this%determine_subset_indices_and_ranges(center_indices, center_positions)
!
      call this%determine_subset_bases(center_indices, center_bases)
!
!     Create the centers in the correct order and store offsets for the subsets
!
      allocate(this%centers(this%n_centers))
!
      call this%create_centers(center_indices, center_bases, &
                               center_symbols, center_positions, &
                               center_charges, center_is_ghosts)
!
      call mem%dealloc(center_indices, this%n_centers)
      call mem%dealloc(center_positions, 3, this%n_centers)
      call mem%dealloc(center_charges, this%n_centers)
!
      deallocate(center_symbols)
      deallocate(center_bases)
!
!     Set up array of non-empty center sets
!
      call this%create_nonempty_subsets()
!
   end subroutine read_centers_atomic_center_reader
!
!
   subroutine create_centers(this, center_indices, center_bases, &
                                   center_symbols, center_positions, &
                                   center_charges, center_is_ghosts)
!!
!!    Create centers
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2020-2021
!!
!!    Creates atomic centers in the order specified by ordered_wfs.
!!
!!    Based on routines with contributions from Sarai D. Folkestad, Eirik F. Kjønstad,
!!    Rolf H. Myhre, and Åsmund H. Tveten.
!!
      implicit none
!
      class(atomic_center_reader), intent(inout) :: this
!
      integer, dimension(this%n_centers), intent(in)            :: center_indices
      character(len=100), dimension(this%n_centers), intent(in) :: center_bases
      character(len=2), dimension(this%n_centers), intent(in)   :: center_symbols
      real(dp), dimension(3, this%n_centers), intent(in)        :: center_positions
      logical, dimension(this%n_centers), intent(in)            :: center_is_ghosts
      integer, dimension(this%n_centers), intent(in)            :: center_charges
!
      integer :: J, J_c
!
      integer, dimension(:), allocatable :: ordered_center_indices
!
!     Order center indices according to basis set within each subset
!
      call mem%alloc(ordered_center_indices, this%n_centers)
!
      call this%order_center_indices_by_basis(ordered_center_indices, center_indices, center_bases)
!
!     Initialize centers according to this reordering
!
      do J = 1, this%n_centers
!
         J_c = ordered_center_indices(J)
!
         this%centers(J) = atomic_center(J,                         &
                                         J_c,                       &
                                         center_symbols(J_c),       &
                                         center_positions(:, J_c),  &
                                         center_bases(J_c),         &
                                         this%basis_type_,          &
                                         center_charges(J_c),       &
                                         center_is_ghosts(J_c))
!
      enddo
!
      call mem%dealloc(ordered_center_indices, this%n_centers)
!
   end subroutine create_centers
!
!
   subroutine order_center_indices_by_basis(this, ordered_center_indices, &
                                                  center_indices, center_bases)
!!
!!    Order center indices by basis
!!    Written by Eirik F. Kjønstad, Jan 2021
!!
!!    Order the center indices by basis set within each subset.
!!
!!    On entry, center_indices are ordered according to subset but not according to
!!    basis within each subset.
!!
      implicit none
!
      class(atomic_center_reader), intent(in) :: this
!
      integer, dimension(this%n_centers), intent(out)           :: ordered_center_indices
      integer, dimension(this%n_centers), intent(in)            :: center_indices
      character(len=100), dimension(this%n_centers), intent(in) :: center_bases
!
      logical, dimension(:), allocatable :: center_reordered ! true for atomic centers that
                                                             ! have been placed in the ordered
                                                             ! center indices array
!
      integer :: center, set, J, J_c, K, K_c
!
      call mem%alloc(center_reordered, this%n_centers)
      center_reordered = .false.
!
      center = 0
!
      do set = 1, this%n_sets
!
         if (this%sets(set)%length .gt. 0) then
!
            do J = this%sets(set)%first, this%sets(set)%get_last()
!
               J_c = center_indices(J)
!
!              For each J, we find the K within the set that has the same basis set
!              and has not already been found and added to ordered_center_indices
!
               do K = J, this%sets(set)%get_last()
!
                  K_c = center_indices(K)
!
                  if (center_reordered(K_c)) cycle                    ! Already reordered? cycle
                  if (center_bases(K_c) .ne. center_bases(J_c)) cycle ! Not the current basis? cycle
!
                  center = center + 1
                  ordered_center_indices(center) = K_c
!
                  center_reordered(K_c) = .true.
!
               enddo
!
            enddo
!
         endif
!
      enddo
!
      call mem%dealloc(center_reordered, this%n_centers)
!
   end subroutine order_center_indices_by_basis
!
!
   subroutine create_nonempty_subsets(this)
!!
!!    Create nonempty subsets
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2020-2021
!!
!!    Allocates and sets an array where empty sets have been removed from
!!    the full list of subsets - in case there are sets that contain zero
!!    atoms.
!!
      implicit none
!
      class(atomic_center_reader), intent(inout) :: this
!
      integer :: set, I
!
      this%n_nonempty_sets = 0
!
      do set = 1, this%n_sets
!
         if (this%sets(set)%length .gt. 0) then
!
            this%n_nonempty_sets = this%n_nonempty_sets + 1
!
         endif
!
      enddo
!
      allocate(this%nonempty_sets(this%n_nonempty_sets))
!
      I = 0
!
      do set = 1, this%n_sets
!
         if (this%sets(set)%length .gt. 0) then
!
            I = I + 1
            this%nonempty_sets(I) = this%sets(set)
!
         endif
!
      enddo
!
   end subroutine create_nonempty_subsets
!
!
   subroutine determine_n_subsets(this)
!!
!!    Determine number of subsets
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2020-2021
!!
!!    Based on routines with contributions from Sarai D. Folkestad, Eirik F. Kjønstad,
!!    Rolf H. Myhre, and Åsmund H. Tveten.
!!
      implicit none
!
      class(atomic_center_reader), intent(inout) :: this
!
      integer :: k
!
      this%n_sets = 1
!
      do k = 1, size(this%ordered_wfs)
!
         if (input%is_keyword_present(trim(this%ordered_wfs(k)), 'active atoms')) then
!
            this%n_sets = this%n_sets + 1
!
         endif
!
      enddo
!
   end subroutine determine_n_subsets
!
!
   subroutine determine_subset_names(this)
!!
!!    Determine subset names
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2020-2021
!!
!!    Center-set names are set to the method of the active space ('ccsd', 'hf', etc.) or
!!    to 'unclassified' if no active spaces are specified in input.
!!
!!    First and length is set to zero here and determined later.
!!
!!    Based on routines with contributions from Sarai D. Folkestad, Eirik F. Kjønstad,
!!    Rolf H. Myhre, and Åsmund H. Tveten.
!!
      implicit none
!
      class(atomic_center_reader), intent(inout) :: this
!
      integer :: set, k
!
      set = 0
!
      do k = 1, size(this%ordered_wfs)
!
         if (input%is_keyword_present(trim(this%ordered_wfs(k)), 'active atoms')) then
!
            set = set + 1
            this%sets(set) = named_range(trim(this%ordered_wfs(k)))
!
         endif
!
      enddo
!
      this%sets(this%n_sets) = named_range('unclassified')
!
   end subroutine determine_subset_names
!
!
   subroutine determine_subset_indices_and_ranges(this, center_indices, center_positions)
!!
!!    Determine subset indices and ranges
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2020
!!
!!    Stores the atomic center indices, for each of the sets in the
!!    center indices arrays, as well as the number of atoms in each
!!    set (i.e., the subset-lengths).
!!
!!    Based on routines with contributions from Sarai D. Folkestad, Eirik F. Kjønstad,
!!    Rolf H. Myhre, and Åsmund H. Tveten.
!!
      implicit none
!
      class(atomic_center_reader), intent(inout) :: this
!
      integer, dimension(this%n_centers), intent(out)    :: center_indices
      real(dp), dimension(3, this%n_centers), intent(in) :: center_positions
!
      character(len=200) :: selection_type
!
      selection_type = 'none'
!
      call input%get_keyword('selection type',   &
                                        'active atoms',     &
                                        selection_type)
!
      if (trim(selection_type) == 'none') then
!
!        Only unclassified centers; one subset with all atoms
!
         if (this%n_sets .ne. 1) &
               call output%error_msg("You need to specify selection type for active atoms.")
!
         call this%determine_subset_indices_and_ranges_unclassified(center_indices)
!
      elseif (trim(selection_type) == 'list' .or. trim(selection_type) == 'range') then
!
!        Active atoms given by list {1,4,5} or range [1-4, 6-10]
!
         call this%determine_subset_indices_and_ranges_list_range(center_indices)
!
      elseif (trim(selection_type) == 'central atom') then
!
!        Central atoms with radii defining the different methods
!
         call this%determine_subset_indices_and_ranges_central_atom(center_indices, center_positions)
!
      endif
!
   end subroutine determine_subset_indices_and_ranges
!
!
   subroutine determine_subset_indices_and_ranges_list_range(this, center_indices)
!!
!!    Determine subset indices and ranges (list or range)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2020-2021
!!
!!    Determines the subset-indices and ranges when active atoms
!!    are given by a list {1,4,5} or a range [1-4, 6-10].
!!
!!    Based on routines with contributions from Sarai D. Folkestad, Eirik F. Kjønstad,
!!    Rolf H. Myhre, and Åsmund H. Tveten.
!!
      implicit none
!
      class(atomic_center_reader), intent(inout) :: this
!
      integer, dimension(this%n_centers), intent(out) :: center_indices
!
      integer :: i, first, length
!
!     Get center indices for atoms associated with methods
!
      first = 1
!
      do i = 1, this%n_sets - 1
!
         length = input%get_n_elements_for_keyword(this%sets(i)%get_name(), &
                                                   'active atoms')
!
         call this%sets(i)%set_range(first, length)
!
         call input%get_array_for_keyword(this%sets(i)%get_name(), &
                                          'active atoms',          &
                                          length,                  &
                                          center_indices(first:))
!
         first = first + length
!
      enddo
!
!     Store the inactive center indices & associated range
!
      call this%determine_subset_indices_and_ranges_unclassified(center_indices)
!
   end subroutine determine_subset_indices_and_ranges_list_range
!
!
   subroutine determine_subset_indices_and_ranges_central_atom(this, center_indices, center_positions)
!!
!!    Determine subset indices and ranges (central atom)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2020-2021
!!
!!    Determines the set-center-indices and range when active atoms
!!    are specified by a central atom with a radius.
!!
!!    Based on routines with contributions from Sarai D. Folkestad, Eirik F. Kjønstad,
!!    Rolf H. Myhre, and Åsmund H. Tveten.
!!
      use array_utilities, only: get_l2_norm
      use parameters, only: zero
!
      implicit none
!
      class(atomic_center_reader), intent(inout) :: this
!
      integer, dimension(this%n_centers), intent(out) :: center_indices
!
      real(dp), dimension(3, this%n_centers), intent(in) :: center_positions
!
      logical, dimension(:), allocatable :: unclassified
!
      integer :: central_atom, center, i, k, first, length
!
      real(dp), dimension(:), allocatable :: radii
!
      real(dp), dimension(3) :: r_vec
!
      real(dp) :: r
!
!     Get the central atom and the radii for different methods
!
      call input%get_keyword('central atom', &
                                        'active atoms', &
                                        central_atom)
!
      call mem%alloc(radii, this%n_sets)
!
      radii = zero
!
      do i = 1, this%n_sets - 1
!
         call input%get_keyword(this%sets(i)%get_name(), &
                                'active atoms',          &
                                radii(i))
!
      enddo
!
!     Set the center indices for the different methods & the number of indices
!
      call mem%alloc(unclassified, this%n_centers)
      unclassified = .true.
!
      first = 1
      center = 0
!
      do i = 1, this%n_sets - 1
!
         do k = 1, this%n_centers
!
            r_vec = center_positions(:, k) - center_positions(:, central_atom)
!
            r = get_l2_norm(r_vec, 3)
!
!           Is center (k) within radii(i) of central atom
!           & has the center not been designated to a previous method?
!
            if (r .lt. radii(i) .and. unclassified(k)) then
!
               center = center + 1
               center_indices(center) = k
!
               unclassified(k) = .false.
!
            endif
!
         enddo
!
         length = center - first + 1
         call this%sets(i)%set_range(first, length)
!
         first = first + center
!
      enddo
!
      call mem%dealloc(radii, this%n_sets)
      call mem%dealloc(unclassified, this%n_centers)
!
!     Store the inactive center indices & associated range
!
      call this%determine_subset_indices_and_ranges_unclassified(center_indices)
!
   end subroutine determine_subset_indices_and_ranges_central_atom
!
!
   subroutine determine_subset_indices_and_ranges_unclassified(this, center_indices)
!!
!!    Determine subset indices and ranges (unclassified)
!!    Written by Eirik F. Kjønstad, 2021
!!
!!    Sets the subset indices and range for unclassified centers
!!
      implicit none
!
      class(atomic_center_reader), intent(inout) :: this
!
      integer, dimension(this%n_centers), intent(inout) :: center_indices
!
      integer :: center, k, j, set, first
!
      logical :: k_is_unclassified
!
      set = this%n_sets ! by assumption, 'unclassified' is the last subset
!
      first = 1
      if (set .gt. 1) first = this%sets(set-1)%get_last() + 1
!
      center = 0
!
      do k = 1, this%n_centers
!
!        Figure out if k is unclassified and store center if this is the case
!
         k_is_unclassified = .true.
!
         do j = 1, first - 1
!
            if (k .eq. center_indices(j)) then
!
               k_is_unclassified = .false.
               exit
!
            endif
!
         enddo
!
         if (k_is_unclassified) then
!
            center = center + 1
            center_indices(first + center - 1) = k
!
         endif
!
      enddo
!
      if (center .gt. 0) call this%sets(set)%set_range(first, center)
!
   end subroutine determine_subset_indices_and_ranges_unclassified
!
!
   subroutine determine_subset_bases(this, center_indices, center_bases)
!!
!!    Determine subset bases
!!    Written by Sarai D. Folkestad, Jul 2021
!!
!!    Changes the basis set of an atomic subsystem, if a different basis set
!!    is specified for the susbet on input
!!
      implicit none
!
      class(atomic_center_reader), intent(inout) :: this
!
      integer, dimension(this%n_centers), intent(in) :: center_indices
      character(len=100), dimension(this%n_centers), intent(inout) :: center_bases
!
      integer :: set, center
      character(len=200) :: set_name
      character(len=200) :: set_basis
!
      do set = 1, this%n_sets
!
         set_name = this%sets(set)%get_name()
!
         if (trim(set_name) == 'unclassified') set_name = 'inactive'
!
         if (input%is_keyword_present(trim(set_name) // ' basis', 'active atoms')) then
!
            call input%get_keyword(trim(set_name) // ' basis', 'active atoms', set_basis)
!
            do center = this%sets(set)%first, this%sets(set)%get_last()
!
               center_bases(center_indices(center)) = trim(set_basis)
!
            enddo
!
         endif
!
      enddo
!
   end subroutine determine_subset_bases
!
!
end module atomic_center_reader_class
