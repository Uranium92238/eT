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
module eri_cd_class
!
!!
!!    Cholesky decomposition (CD) of electronic
!!    repulsion integrals (ERI) solver class
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!    Files updated by Rolf H. Myhre, September 2019
!!
!!    Handles the Cholesky decomposition of the ERIs
!!
!!       g_wxyz = (wx|yz) = sum_J L_wx_J L_yz_J
!!
!!    by first determining the Cholesky basis, i.e.,
!!    the decomposition pivots and then constructing the
!!    vectors.
!!
!!    After the basis is determined, the overlap
!!
!!       S_JK = (J|K),
!!
!!    is constructed for the pivot elements {J}. S is Cholesky
!!    decomposed,
!!
!!       S = Q^T Q,
!!
!!    and the factors Q are inverted.
!!
!!    Finally the Cholesky vectors are constructed
!!    as
!!
!!       L_wx_J = sum_K (wx|K)[Q^-T]_KJ.
!!
!!    For a more detailed description see
!!
!!       Folkestad, S. D., Kjønstad, E. F., and Koch, H.,
!!       JCP, 150(19), 194112 (2019).
!!
!!    There are options for one-center Cholesky decomposition
!!    method specific Cholesky decomposition or partitioned
!!    Cholesky decomposition (PCD).
!!
!
   use parameters
!
   use global_out, only: output
   use memory_manager_class, only: mem
   use range_class, only: range_
   use sequential_file_class, only: sequential_file
   use cholesky_array_list_class, only: cholesky_array_list
   use ao_tool_class, only: ao_tool
   use timings_class, only: timings
!
   implicit none
!
!
   type :: eri_cd
!
      character(len=100), private :: name_ = 'Cholesky decomposition of electronic &
                                            &repulsion integrals solver'

      character(len=500), private :: description1 = 'Performs a Cholesky decomposition &
                                             &of the two-electron electronic repulsion integrals &
                                             &in the atomic orbital basis,'

      character(len=500), private :: description2 = '(ab|cd) = sum_J L_ab^J L_cd^J.'

      character(len=500), private :: description3 = 'Once the Cholesky basis has been determined, &
                                            &the vectors L^J are constructed and stored to disk. &
                                            &These may either be used directly, &
                                            &or be transformed to the MO basis &
                                            &for use in post-HF calculations. &
                                            &For more information, &
                                            &see S. D. Folkestad, E. F. Kjønstad &
                                            &and H. Koch, JCP, 150(19), (2019)'
!
      real(dp), private :: threshold
      real(dp), private :: span
!
      integer, private :: max_qual
      integer, private :: iteration
!
      logical, private :: one_center
!
      type(sequential_file), private :: diagonal_info_cauchy_schwarz
      type(sequential_file), private :: Q
      type(sequential_file), private :: Q_inverse
      type(sequential_file), private :: diagonal_info_target
      type(sequential_file), private :: cholesky_basis_file
!
      integer, private :: n_cholesky
      integer, private :: n_shp_in_basis
      integer, private :: n_s, n_shp, n_ao, n_aop
!
      integer, private :: n_batches
!
      type(timings), private :: timer
!
   contains
!
      procedure, public :: run &
                        => run_eri_cd
!
      procedure :: construct_cholesky_mo_vectors_eri_cd
!
      generic, public :: construct_cholesky_mo_vectors => construct_cholesky_mo_vectors_eri_cd
!
      procedure, public :: diagonal_test &
                        => diagonal_test_eri_cd
!
      procedure, public :: get_n_cholesky &
                        => get_n_cholesky_eri_cd
!
      procedure, public :: cleanup &
                        => cleanup_eri_cd
!
      procedure, private :: construct_significant_diagonal
!
      procedure, private :: construct_shp_to_shells
      procedure, private :: determine_sig_shps_and_max_diagonal
      procedure, private :: determine_shps_to_construct_diagonal
      procedure, private :: calculate_n_sig_aops_and_shps
      procedure, private :: get_screening_vector
      procedure, private :: construct_sig_shp_to_shells_and_ao_offsets
      procedure, private :: construct_sig_diagonal_and_screening_vector
      procedure, private :: write_diagonal_info_target
      procedure, private :: write_diagonal_info_cauchy_schwarz
!
      procedure, private :: determine_cholesky_basis
!
      procedure, private :: determine_cholesky_basis_standard
      procedure, private :: determine_cholesky_basis_PCD
!
      procedure, private :: construct_diagonal_batches
      procedure, private :: construct_diagonal_from_batch_bases
!
      procedure, private :: invert_Q
      procedure, private :: construct_S
!
      procedure, private :: read_settings
      procedure, private :: print_banner
      procedure, private :: print_settings
!
   end type eri_cd
!
!
   interface eri_cd
!
      procedure :: new_eri_cd
!
   end interface eri_cd
!
!
contains
!
!
   function new_eri_cd(ao) result(this)
!!
!!    New ERI CD
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      use citation_printer_class, only: eT_citations
!
      implicit none
!
      type(eri_cd)  :: this
      type(ao_tool) :: ao
!
      this%timer = timings('Cholesky decomposition of ERIs')
      call this%timer%turn_on()
!
!     Set defaults
!
      this%threshold           = 1.0d-4
      this%span                = 1.0d-2
!
      this%max_qual            = 1000
      this%iteration           = 0
      this%n_batches           = 1
!
      this%one_center          = .false.
!
      this%n_cholesky          = 0
!
      call this%read_settings()
!
      this%n_aop   = ao%n*(ao%n+1)/2 ! Number of ao pairs packed
      this%n_ao    = ao%n
      this%n_s     = ao%n_sh
      this%n_shp   = this%n_s*(this%n_s + 1)/2 ! Number of shell pairs packed
!
!     Add citation for this implementation
!
      call eT_citations%add("Cholesky decomposition of ERIs")
!
!     Initialize files
!
      this%diagonal_info_cauchy_schwarz = sequential_file('cauchy_schwarz_diagonal_eri')
!
      this%Q = sequential_file('Q_eri')
      this%Q_inverse = sequential_file('Q_inverse_eri')
!
      this%diagonal_info_target = sequential_file('target_diagonal_eri')
      this%cholesky_basis_file = sequential_file('basis_shell_info')
!
      call this%print_banner()
      call this%print_settings()
!
!     Additional prints
!
      call output%printf('m', '- Cholesky decomposition ao details:', fs='(/t3, a)')
!
      call output%printf('m', 'Total number of AOs:         (i13)', &
                         ints=[this%n_ao], fs='(/t6,a)')
      call output%printf('m', 'Total number of shell pairs: (i13)', &
                         ints=[this%n_shp], fs='(t6, a)')
      call output%printf('m', 'Total number of AO pairs:    (i13)', &
                         ints=[this%n_aop], fs='(t6, a)')
!
   end function new_eri_cd
!
!
   subroutine run_eri_cd(this, ao, screening_vector)
!!
!!    Run
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(eri_cd) :: this
!
      type(ao_tool) :: ao
!
      real(dp), dimension(this%n_ao), optional :: screening_vector
!
      type(timings) :: basis_timer, inversion_timer
!
      basis_timer = timings("Cholesky: time to determine auxiliary basis")
      call basis_timer%turn_on()
!
      call this%construct_significant_diagonal(ao, screening_vector)
      call this%determine_cholesky_basis(ao)
!
      call basis_timer%turn_off()
!
      inversion_timer = timings("Cholesky: time to construct (J|K), decompose, and invert")
      call inversion_timer%turn_on()
!
      call this%construct_S(ao) ! S = (J | K), where J and K are pivots in the basis
      call this%invert_Q()      ! Q is the Cholesky factor: S = QQ^T
!
      call inversion_timer%turn_off()
!
   end subroutine run_eri_cd
!
!
   subroutine determine_cholesky_basis(this, ao)
!!
!!    Determine Cholesky basis
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Determine Cholesky pivots in batches (partitioned Cholesky decomposition, PCD) or
!!    in a single batch (standard Cholesky decomposition).
!!
      implicit none
!
      class(eri_cd), intent(inout) :: this
!
      class(ao_tool), intent(in) :: ao
!
      if (this%n_batches == 1) then
!
         call this%determine_cholesky_basis_standard(ao,                        &
                                                     this%diagonal_info_target, &
                                                     this%cholesky_basis_file)
!
      else
!
         call this%determine_cholesky_basis_PCD(ao)
!
      endif
!
   end subroutine determine_cholesky_basis
!
!
   pure function get_n_cholesky_eri_cd(this) result(n_J)
!!
!!    Get n Cholesky
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none
!
      class(eri_cd), intent(in) :: this
!
      integer :: n_J
!
      n_J = this%n_cholesky
!
   end function get_n_cholesky_eri_cd
!
!
   subroutine cleanup_eri_cd(this)
!!
!!    Cleanup
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(eri_cd), intent(inout) :: this
!
      call this%timer%turn_off()
!
      call output%printf('n', '- Finished decomposing the ERIs.', fs='(/t3, a)')
!
      call output%printf('n', 'Total wall time (sec): (f20.5)', &
                         reals=[this%timer%get_elapsed_time('wall')], fs='(/t6,a)')
      call output%printf('n', 'Total cpu time (sec):  (f20.5)', &
                         reals=[this%timer%get_elapsed_time('cpu')], fs='(t6, a)')
!
   end subroutine cleanup_eri_cd
!
!
   subroutine construct_significant_diagonal(this, ao, screening_vector)
!!
!!    Construct significant diagonal
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Constructs two screened diagonals
!!
!!       1. Screened diagonal for decomposition D_wx <= T   or   D_wx * V_w * V_x <= T
!!          The latter if optional screening vector is present
!!
!!       2. Screened diagonal for construction of cholesky vectros  sqrt( D_wx * D_max ) <= T
!!
!!    Writes all information to files target_diagonal and construct_diagonal
!!
      implicit none
!
      class(eri_cd)  :: this
      class(ao_tool) :: ao
!
      real(dp), dimension(this%n_ao*this%n_ao), target, optional :: screening_vector
!
      integer :: n_sig_aop, n_sig_shp, n_construct_aop, n_construct_shp
!
      integer, dimension(:,:), allocatable :: shp_to_shells, sig_shp_to_shells
!
      integer, dimension(:), allocatable  :: ao_offsets
!
      real(dp), dimension(:,:), allocatable :: screening_vector_local
!
      real(dp), dimension(:), allocatable :: screening_vector_reduced
      real(dp), dimension(:), allocatable :: D_xy
!
      logical, dimension(:), allocatable :: sig_shp, construct_shp
!
      real(dp) :: max_diagonal
!
      call output%printf('v', '- Preparing diagonal for decomposition:', fs='(/t3,a)')
!
!
      call mem%alloc(screening_vector_local, this%n_ao, this%n_ao)
      call mem%alloc(sig_shp, this%n_shp)
      call mem%alloc(shp_to_shells, this%n_shp, 2)
!
      call this%get_screening_vector(screening_vector_local, screening_vector)
      call this%construct_shp_to_shells(shp_to_shells)
!
      call this%determine_sig_shps_and_max_diagonal(ao,                     &
                                                    sig_shp,                & ! = D(AB) significant?
                                                    max_diagonal,           & ! = max(D(AB))
                                                    shp_to_shells,          &
                                                    screening_vector_local)
!
      call mem%alloc(construct_shp, this%n_shp)
!
      call this%determine_shps_to_construct_diagonal(ao,            &
                                                     construct_shp, &
                                                     shp_to_shells, &
                                                     max_diagonal)
!
      call this%calculate_n_sig_aops_and_shps(ao,                &
                                              n_sig_aop,         &
                                              n_sig_shp,         &
                                              sig_shp,           &
                                              shp_to_shells)
!
      call output%printf('n', 'Significant shell pairs: (i17)', ints=[n_sig_shp], fs='(/t6,a)')
      call output%printf('n', 'Significant AO pairs:    (i17)', ints=[n_sig_aop], fs='(t6, a)')
!
      call this%calculate_n_sig_aops_and_shps(ao,                &
                                              n_construct_aop,   &
                                              n_construct_shp,   &
                                              construct_shp,     &
                                              shp_to_shells)
!
      call output%printf('n', 'Construct shell pairs: (i19)', ints=[n_construct_shp], fs='(/t6,a)')
      call output%printf('n', 'Construct AO pairs:    (i19)', ints=[n_construct_aop], fs='(t6, a)')
!
!     Construct index lists needed for parallelized construction of the diagonal
!
      call mem%alloc(ao_offsets, n_sig_shp)
      call mem%alloc(sig_shp_to_shells, n_sig_shp, 2)
!
      call this%construct_sig_shp_to_shells_and_ao_offsets(ao,                    &
                                                           sig_shp_to_shells,     &
                                                           n_sig_shp,             &
                                                           ao_offsets,            &
                                                           shp_to_shells,         &
                                                           sig_shp)
!
      call mem%dealloc(shp_to_shells, this%n_shp, 2)
!
!     Construct significant diagonal and screening vector
!
      call mem%alloc(D_xy, n_sig_aop)
      call mem%alloc(screening_vector_reduced, n_sig_aop)
!
      call this%construct_sig_diagonal_and_screening_vector(ao,                        &
                                                            D_xy,                      &
                                                            screening_vector_reduced,  &
                                                            screening_vector_local,    &
                                                            sig_shp_to_shells,         &
                                                            n_sig_shp,                 &
                                                            n_sig_aop,                 &
                                                            ao_offsets)

!
      call mem%dealloc(screening_vector_local, this%n_ao, this%n_ao)
      call mem%dealloc(sig_shp_to_shells, n_sig_shp, 2)
      call mem%dealloc(ao_offsets, n_sig_shp)
!
      call this%write_diagonal_info_target(n_sig_shp, &
                                           n_sig_aop, &
                                           sig_shp,   &
                                           D_xy,      &
                                           screening_vector_reduced)
!
      call this%write_diagonal_info_cauchy_schwarz(n_construct_shp, &
                                                   n_construct_aop, &
                                                   construct_shp)
!
      call mem%dealloc(sig_shp, this%n_shp)
      call mem%dealloc(D_xy, n_sig_aop)
      call mem%dealloc(screening_vector_reduced, n_sig_aop)
!
      call mem%dealloc(construct_shp, this%n_shp)
!
   end subroutine construct_significant_diagonal
!
!
   subroutine write_diagonal_info_target(this, n_sig_shp, n_sig_aop, sig_shp, &
                                         D_xy, screening_vector_reduced)
!!
!!    Write diagonal info target
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Write info file for target diagonal containing
!!
!!       1. number of significant shell pairs, number of significant ao pairs
!!       2. sig_shp - vector of logicals to describe which shell pairs are significant
!!       3. D_xy = ( xy | xy ), the significant diagonal.
!!       4. Screening vector
!!
      implicit none
!
      class(eri_cd), intent(inout) :: this
!
      integer, intent(in) :: n_sig_shp, n_sig_aop
!
      logical, dimension(n_sig_shp), intent(in) :: sig_shp
!
      real(dp), dimension(n_sig_aop), intent(in) :: D_xy
!
      real(dp), dimension(n_sig_aop), intent(in) :: screening_vector_reduced
!
      call this%diagonal_info_target%open_('write', 'rewind')
!
      call this%diagonal_info_target%write_(n_sig_shp)
      call this%diagonal_info_target%write_(n_sig_aop)
      call this%diagonal_info_target%write_(sig_shp, this%n_shp)
      call this%diagonal_info_target%write_(D_xy, n_sig_aop)
      call this%diagonal_info_target%write_(screening_vector_reduced, n_sig_aop)
!
      call this%diagonal_info_target%close_()
!
   end subroutine write_diagonal_info_target
!
!
   subroutine write_diagonal_info_cauchy_schwarz(this, n_construct_shp, &
                                                     n_construct_aop, construct_shp)
!!
!!    Write diagonal info Cauchy-Schwarz
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Write info file for construct diagonal containing
!!
!!       1. number of shell pairs to construct, number of ao pairs to construct
!!       2. construct_shp - vector of logicals to describe which shell pairs are to be constructed
!!
      implicit none
!
      class(eri_cd), intent(inout) :: this
!
      integer, intent(in) :: n_construct_shp, n_construct_aop
!
      logical, dimension(n_construct_shp), intent(in) :: construct_shp
!
      call this%diagonal_info_cauchy_schwarz%open_('write', 'rewind')
!
      call this%diagonal_info_cauchy_schwarz%write_(n_construct_shp)
      call this%diagonal_info_cauchy_schwarz%write_(n_construct_aop)
      call this%diagonal_info_cauchy_schwarz%write_(construct_shp, this%n_shp)
!
      call this%diagonal_info_cauchy_schwarz%close_()
!
   end subroutine write_diagonal_info_cauchy_schwarz
!
!
   subroutine construct_sig_diagonal_and_screening_vector(this, ao, &
                                          D_xy, screening_vector_reduced, screening_vector_local, &
                                          sig_shp_to_shells, n_sig_shp, n_sig_aop, ao_offsets)
!!
!!    Construct significant diagonal and screening vector
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(eri_cd), intent(in) :: this
!
      class(ao_tool), intent(in) :: ao
!
      integer, intent(in) :: n_sig_shp, n_sig_aop
!
      real(dp), dimension(n_sig_aop), intent(out) :: D_xy
      real(dp), dimension(n_sig_aop), intent(out) :: screening_vector_reduced
!
      real(dp), dimension(this%n_ao, this%n_ao), intent(in) :: screening_vector_local
      integer, dimension(n_sig_shp, 2), intent(in) :: sig_shp_to_shells
      integer, dimension(n_sig_shp), intent(in) :: ao_offsets
!
      real(dp), dimension(:,:,:,:), pointer :: g_ABAB_p
!
      real(dp), dimension(ao%max_sh_size**4), target :: g_ABAB
!
      integer :: x, y, xy, xy_packed, A, B, I
!
      type(range_) :: A_range, B_range
!
!     Note: allocated with length n_significant_shp + 1, last element is used for n_significant_aop
!     This is convenient because significant_shp_to_first_significant_aop will be used to calculate lengths.
!
!$omp parallel do &
!$omp private(I, A, B, A_range, B_range, x, y, xy, xy_packed, g_ABAB, g_ABAB_p) &
!$omp shared(D_xy, screening_vector_reduced, ao_offsets) &
!$omp schedule(guided)
      do I = 1, n_sig_shp
!
         A = sig_shp_to_shells(I, 1)
         B = sig_shp_to_shells(I, 2)
!
         A_range = ao%shells(A)
         B_range = ao%shells(B)
!
         call ao%get_eri(g_ABAB, A, B, A, B)
!
         g_ABAB_p(1 : A_range%length, 1 : B_range%length, &
                  1 : A_range%length, 1 : B_range%length) &
                  => g_ABAB(1 : (A_range%length)**2*(B_range%length)**2)
!
         if (A .eq. B) then
!
            do x = 1, A_range%length
               do y = x, B_range%length
!
                  xy_packed = (max(x,y)*(max(x,y)-3)/2) + x + y
!
                  D_xy(xy_packed + ao_offsets(I)) = g_ABAB_p(x, y, x, y)
                  screening_vector_reduced(xy_packed + ao_offsets(I)) = &
                                          screening_vector_local(x + A_range%first - 1, &
                                                                 y + B_range%first - 1)
               enddo
            enddo
!
         else ! A ≠ B
!
            do x = 1, (A_range%length)
               do y = 1, (B_range%length)
!
                  xy = A_range%length*(y - 1) + x
                  D_xy(xy + ao_offsets(I)) = g_ABAB_p(x, y, x, y)
                  screening_vector_reduced(xy + ao_offsets(I)) = &
                                      screening_vector_local(x + A_range%first - 1, &
                                                             y + B_range%first - 1)
!
               enddo
            enddo
!
         endif
!
      enddo
!$omp end parallel do
!
   end subroutine construct_sig_diagonal_and_screening_vector
!
!
   subroutine construct_sig_shp_to_shells_and_ao_offsets(this, ao, &
                        sig_shp_to_shells, n_sig_shp, ao_offsets, shp_to_shells, sig_shp)
!!
!!    Construct sig shp to shells and AO offsets
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(eri_cd), intent(in) :: this
!
      class(ao_tool), intent(in) :: ao
!
      integer, intent(in) :: n_sig_shp
!
      integer, dimension(n_sig_shp, 2), intent(out) :: sig_shp_to_shells
!
      integer, dimension(n_sig_shp), intent(out) :: ao_offsets
!
      integer, dimension(this%n_shp, 2), intent(in) :: shp_to_shells
      logical, dimension(this%n_shp), intent(in) :: sig_shp
!
      integer :: I, A, B, current_sig_shp
!
      type(range_) :: A_range, B_range
!
      ao_offsets = 0
      current_sig_shp = 0
!
      do I = 1, this%n_shp
!
         if (sig_shp(I)) then
!
            current_sig_shp = current_sig_shp + 1
!
            A = shp_to_shells(I, 1)
            B = shp_to_shells(I, 2)
!
            A_range = ao%shells(A)
            B_range = ao%shells(B)
!
            sig_shp_to_shells(current_sig_shp, 1) = A
            sig_shp_to_shells(current_sig_shp, 2) = B
!
            if (current_sig_shp .lt. n_sig_shp) then
!
               ao_offsets(current_sig_shp + 1) = ao_offsets(current_sig_shp) + &
                           get_size_shp(A_range, B_range)
!
            endif
!
         endif
!
      enddo
!
   end subroutine construct_sig_shp_to_shells_and_ao_offsets
!
!
   subroutine get_screening_vector(this, screening_target, screening_source)
!!
!!    Get screening vector
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Sets target screening vector equal to source (if 'source' is present).
!!
!!    Otherwise, the target screening vector is set to 1.
!!
      implicit none
!
      class(eri_cd), intent(in) :: this
!
      real(dp), dimension(this%n_ao, this%n_ao), intent(out) :: screening_target
      real(dp), dimension(this%n_ao, this%n_ao), intent(in), optional :: screening_source
!
      integer :: i, j
!
      if (present(screening_source)) then
!
         call dcopy(this%n_ao**2, screening_source, 1, screening_target, 1)
!
      else
!
!$omp parallel do private (i, j)
         do i = 1, this%n_ao
            do j = 1, this%n_ao
!
               screening_target(i, j) = one
!
            enddo
         enddo
!$omp end parallel do
!
      endif
!
   end subroutine get_screening_vector
!
!
   subroutine calculate_n_sig_aops_and_shps(this, ao, n_aop, n_shp, sig_shp, shp_to_shells)
!!
!!    Calculate n sig AOPs and SHPs
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Compute number of AOPs and SHPs (I) that correspond to sig_shp(I) = .true.
!!
      use array_utilities, only: is_significant
!
      implicit none
!
      class(eri_cd), intent(in) :: this
!
      class(ao_tool), intent(in) :: ao
!
      integer, intent(out) :: n_aop, n_shp
!
      logical, dimension(this%n_shp), intent(in) :: sig_shp
!
      integer, dimension(this%n_shp, 2), intent(in) :: shp_to_shells
!
      integer :: I, A, B
!
      type(range_) :: A_range, B_range
!
      n_aop = 0 ! Number of sig AO pairs
      n_shp = 0 ! Number of sig shell pairs
!
      do I = 1, this%n_shp
!
         if (sig_shp(I)) then
!
            A = shp_to_shells(I, 1)
            B = shp_to_shells(I, 2)
!
            A_range = ao%shells(A)
            B_range = ao%shells(B)
!
            n_aop = n_aop + get_size_shp(A_range, B_range)
            n_shp = n_shp + 1
!
         endif
!
      enddo
!
   end subroutine calculate_n_sig_aops_and_shps
!
!
   subroutine determine_sig_shps_and_max_diagonal(this, ao, &
                        sig_shp, max_diagonal, shp_to_shells, &
                        screening_vector)
!!
!!    Determine sig SHPs and max diagonal
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Determines the shell pairs that are significant (sig_shp) and
!!    calculates the maximum value on the entire diagonal (max_diagonal).
!!
      use array_utilities, only: is_significant
!
      implicit none
!
      class(eri_cd), intent(in) :: this
!
      class(ao_tool), intent(in) :: ao
!
      logical, dimension(this%n_shp), intent(out) :: sig_shp
!
      real(dp), intent(out) :: max_diagonal
!
      integer, dimension(this%n_shp, 2), intent(in) :: shp_to_shells
!
      real(dp), dimension(this%n_ao, this%n_ao), intent(in) :: screening_vector
!
      real(dp), dimension(:,:,:,:), pointer :: g_ABAB_p
!
      real(dp), dimension(ao%max_sh_size**4), target :: g_ABAB
      real(dp), dimension(ao%max_sh_size**2) :: D_AB, D_AB_screen
!
      integer :: x, y, A, B, I, K
!
      type(range_) :: A_range, B_range
!
      real(dp), dimension(:), allocatable :: max_in_shp_diagonal
!
      sig_shp = .false.
!
      call mem%alloc(max_in_shp_diagonal, this%n_shp)
!
!$omp parallel do &
!$omp private(I, K, A, B, A_range, B_range, x, y, g_ABAB, g_ABAB_p, D_AB, D_AB_screen) &
!$omp shared(sig_shp,  max_in_shp_diagonal) &
!$omp schedule(guided)
      do I = 1, this%n_shp
!
         A = shp_to_shells(I, 1)
         B = shp_to_shells(I, 2)
!
         if (this%one_center .and. &
             ao%shell_to_center(B) .ne. ao%shell_to_center(A)) cycle
!
         A_range = ao%shells(A)
         B_range = ao%shells(B)
!
!        Construct diagonal D_AB for the given shell pair
!
         call ao%get_eri(g_ABAB, A, B, A, B)
!
         g_ABAB_p(1 : A_range%length, 1 : B_range%length, &
                  1 : A_range%length, 1 : B_range%length) &
                  => g_ABAB(1 : (A_range%length)**2*(B_range%length)**2)
!
         K = 0
         do x = 1, (A_range%length)
            do y = 1, (B_range%length)
!
               K = K + 1
               D_AB_screen(K) = g_ABAB_p(x, y, x, y)&
                           *screening_vector(x + A_range%first - 1, &
                                             y + B_range%first - 1)
!
               D_AB(K) = g_ABAB_p(x, y, x, y)
!
            enddo
         enddo
!
         sig_shp(I) = is_significant(D_AB_screen, A_range%length*B_range%length, this%threshold)
!
         max_in_shp_diagonal(I) = maxval(D_AB)
!
      enddo
!$omp end parallel do
!
      max_diagonal = maxval(max_in_shp_diagonal)
!
      call mem%dealloc(max_in_shp_diagonal, this%n_shp)
!
   end subroutine determine_sig_shps_and_max_diagonal
!
!
   subroutine determine_shps_to_construct_diagonal(this, ao, &
                                 construct_shp, shp_to_shells, max_diagonal)
!!
!!    Determine SHPs to construct diagonal
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Determines for which SHPs to construct the diagonal of the ERI matrix
!!    (construct_shp).
!!
      use array_utilities, only: is_significant
!
      implicit none
!
      class(eri_cd), intent(in) :: this
!
      class(ao_tool), intent(in) :: ao
!
      logical, dimension(this%n_shp), intent(out) :: construct_shp
!
      integer, dimension(this%n_shp, 2), intent(in) :: shp_to_shells
!
      real(dp), intent(in) :: max_diagonal
!
      real(dp), dimension(:,:,:,:), pointer :: g_ABAB_p
!
      real(dp), dimension(ao%max_sh_size**4), target :: g_ABAB
!
      integer :: I, K, A, B, x, y
!
      type(range_) :: A_range, B_range
!
      real(dp), dimension(ao%max_sh_size**2) :: construct_test
!
      construct_shp = .false.
!
!$omp parallel do &
!$omp private(I, K, A, B, A_range, B_range, x, y, g_ABAB, g_ABAB_p, construct_test) &
!$omp shared(construct_shp) &
!$omp schedule(guided)
      do I = 1, this%n_shp
!
         A = shp_to_shells(I, 1)
         B = shp_to_shells(I, 2)
!
         A_range = ao%shells(A)
         B_range = ao%shells(B)
!
         call ao%get_eri(g_ABAB, A, B, A, B)
!
         g_ABAB_p(1 : A_range%length, 1 : B_range%length, &
                  1 : A_range%length, 1 : B_range%length) &
                  => g_ABAB(1 : (A_range%length)**2*(B_range%length)**2)
!
         K = 0
         do x = 1, A_range%length
            do y = 1, B_range%length
!
               K = K + 1
               construct_test(K) = sqrt(g_ABAB_p(x, y, x, y)*max_diagonal)
!
            enddo
         enddo
!
         construct_shp(I) = is_significant(construct_test, &
                           (A_range%length)*(B_range%length), min(this%threshold,1.0d-8))
!
      enddo
!$omp end parallel do
!
   end subroutine determine_shps_to_construct_diagonal
!
!
   subroutine construct_shp_to_shells(this, shp_to_shells)
!!
!!    Construct shp to shells
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Creates map from a shell-pair AB to the individual shells [A, B].
!!
      implicit none
!
      class(eri_cd), intent(in) :: this
!
      integer, dimension(this%n_shp, 2), intent(out) :: shp_to_shells
!
      integer :: shp, A, B
!
      shp = 0
!
      do B = 1, this%n_s
         do A = B, this%n_s
!
            shp = shp + 1
!
            shp_to_shells(shp, 1) = A
            shp_to_shells(shp, 2) = B
!
         enddo
      enddo
!
   end subroutine construct_shp_to_shells
!
!
   subroutine construct_diagonal_batches(this, ao)
!!
!!    Construct diagonal batches
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Divides the significant diagonal into batches and prepares for
!!    partitioned decomposition
!!
      implicit none
!
      class(eri_cd) :: this
!
      type(ao_tool) :: ao
!
      integer :: n_sig_aop, n_sig_shp, n_sig_shp_batch, shp
!
      real(dp), dimension(:), allocatable :: D_xy, D_batch
!
      real(dp), dimension(:), allocatable :: screening_vector_batch, screening_vector
!
      logical, dimension(:), allocatable :: sig_shp, sig_shp_batch
!
      type(range_) :: A_range, B_range
!
      type(sequential_file) :: batch_file
!
      integer :: A, B, batch, batch_first, batch_last, batch_size, current_batch_size
      integer :: xy_first, xy_last
!
      character(len=100) :: temp_name
!
!     Read diagonal info file containing (name given as argument)
!
!        1. number of significant shell pairs, number of significant ao pairs
!        2. sig_shp - vector of logicals to describe which shell pairs are significant
!        3. D_xy = ( xy | xy ), the significant diagonal.
!
      call this%diagonal_info_target%open_('read', 'rewind')
!
      call this%diagonal_info_target%read_(n_sig_shp)
      call this%diagonal_info_target%read_(n_sig_aop)
!
      call mem%alloc(sig_shp, this%n_shp)
      call mem%alloc(D_xy, n_sig_aop)
      call mem%alloc(screening_vector, n_sig_aop)
!
      call this%diagonal_info_target%read_(sig_shp, this%n_shp)
      call this%diagonal_info_target%read_(D_xy, n_sig_aop)
      call this%diagonal_info_target%read_(screening_vector, n_sig_aop)
!
      call this%diagonal_info_target%close_()
!
!     Calculate size of batches
!
      batch_size = n_sig_aop/this%n_batches
!
      call mem%alloc(sig_shp_batch, (this%n_shp))
!
      batch_first = 1
      batch_last = batch_size
!
      do batch = 1, this%n_batches
!
!        Determine sig_shp_batch
!
         sig_shp_batch = .false.
!
         shp = 0        ! Shell pair number
         xy_first = 1
         xy_last = 0
         n_sig_shp_batch = 0
!
         do B = 1, this%n_s
            do A = B, this%n_s
!
               shp = shp + 1
!
               if (sig_shp(shp)) then
!
                  A_range = ao%shells(A)
                  B_range = ao%shells(B)
!
                  xy_last = xy_last + get_size_shp(A_range, B_range)
!
                  if ((xy_last .ge. batch_first) .and. (xy_first .le. batch_last)) then
!
                     sig_shp_batch(shp) = .true.
                     n_sig_shp_batch = n_sig_shp_batch + 1
!
                     if (xy_last .gt. batch_last) then
!
                        batch_last = xy_last
!
                     endif
!
                  endif
!
                  xy_first = xy_first + get_size_shp(A_range, B_range)
!
               endif
!
            enddo
         enddo
!
         current_batch_size = batch_last - batch_first + 1
!
         call mem%alloc(D_batch, current_batch_size)
         call mem%alloc(screening_vector_batch, current_batch_size)
!
         D_batch(:) = D_xy(batch_first : batch_last)
!
         screening_vector_batch(:) = screening_vector(batch_first : batch_last)
!
!        Write info file for batch diagonal containing
!
!        1. number of significant shell pairs, number of significant ao pairs
!        2. sig_shp - vector of logicals to describe which shell pairs are significant
!        3. D_xy = ( xy | xy ), the significant diagonal.
!        4. Screening vector
!
         write(temp_name, '(a14, i4.4)')'diagonal_info_', batch
         batch_file = sequential_file(trim(temp_name))
!
         call batch_file%open_('write', 'rewind')
!
         call output%printf('n', 'Significant AO and shell pairs in batch (i0):', &
                            ints=[batch], fs='(/t6,a)')
         call output%printf('n', 'Significant shell pairs: (i14)', &
                            ints=[n_sig_shp_batch], fs='(t9,a)')
         call output%printf('n', 'Significant AO pairs:    (i14)', &
                            ints=[current_batch_size], fs='(t9,a)')
!
         call batch_file%write_(n_sig_shp_batch)
         call batch_file%write_(current_batch_size)
         call batch_file%write_(sig_shp_batch, this%n_shp)
         call batch_file%write_(D_batch, current_batch_size)
         call batch_file%write_(screening_vector_batch, current_batch_size)
!
         call batch_file%close_()
!
         call mem%dealloc(D_batch, current_batch_size)
         call mem%dealloc(screening_vector_batch, current_batch_size)
!
         batch_first = batch_last + 1
         batch_last  = batch_size*(batch + 1)
!
         if ((batch + 1) == this%n_batches) batch_last = n_sig_aop
!
         if (batch_last .lt. batch_first) call output%error_msg('Batch size is too small.')
!
      enddo
!
      call mem%dealloc(sig_shp_batch, (this%n_shp))
      call mem%dealloc(sig_shp, (this%n_shp))
      call mem%dealloc(D_xy, n_sig_aop)
      call mem%dealloc(screening_vector, n_sig_aop)
!
   end subroutine construct_diagonal_batches
!
!
   subroutine construct_diagonal_from_batch_bases(this, ao, &
                                                  n_cholesky_batches, n_shp_in_basis_batches)
!!
!!    Construct diagonal from batch bases
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Nov 2018
!!
!!    Constructs the final diagonal from the bases obtained from diagonal batches.
!!    Called as preparation for final decomposition step in PCD.
!!
      use array_utilities, only: quicksort_with_index_ascending_int
!
      implicit none
!
      class(eri_cd) :: this
!
      type(ao_tool) :: ao
!
      integer, dimension(this%n_batches), intent(in) :: n_cholesky_batches
      integer, dimension(this%n_batches), intent(in) :: n_shp_in_basis_batches
!
      integer :: n_cholesky_total, n_shp_in_basis_total, J, I, n_sig_aop, n_sig_shp
      integer :: A_shell, B_shell, shp, alpha_in_A, beta_in_B, alpha_beta_in_AB, aop, batch
      integer :: n_basis_aop_in_AB_total, n_basis_aop_in_AB_offset, current_offset, current_offset_old
      integer :: count_sig, n_cholesky_offset, n_sig_aop_old, n_sig_shp_old, n_shp_in_basis_offset
!
      type(sequential_file) :: batch_file
!
      integer, dimension(:), allocatable :: alpha, beta, alpha_beta, sorted_alpha, sorted_beta, sorted_alpha_beta
      integer, dimension(:), allocatable :: index_alpha_beta, alpha_beta_offset, alpha_beta_offset_old
      integer, dimension(:), allocatable :: A, B, AB, sorted_A, sorted_B, sorted_AB
      integer, dimension(:), allocatable :: sorted_n_basis_aop_in_AB, n_basis_aop_in_AB, index_AB
!
      integer, dimension(:,:), allocatable :: basis_shell_info, cholesky_basis
!
      logical, dimension(:), allocatable :: sig_shp, sig_shp_old
!
      type(range_) :: A_range, B_range
!
      real(dp), dimension(:), allocatable :: D, D_old, screening_vector, screening_vector_old
!
      character(len=100) :: temp_name
!
      n_cholesky_total = 0
      n_shp_in_basis_total = 0
!
      do batch = 1, this%n_batches
!
         n_cholesky_total    = n_cholesky_total + n_cholesky_batches(batch)
         n_shp_in_basis_total = n_shp_in_basis_total + n_shp_in_basis_batches(batch)
!
      enddo
!
!     Read and paste together basis information
!     from the different batches
!
      call mem%alloc(alpha, n_cholesky_total)
      call mem%alloc(beta, n_cholesky_total)
      call mem%alloc(alpha_beta, n_cholesky_total)
!
      call mem%alloc(A, n_shp_in_basis_total)
      call mem%alloc(B, n_shp_in_basis_total)
      call mem%alloc(AB, n_shp_in_basis_total)
      call mem%alloc(n_basis_aop_in_AB, n_shp_in_basis_total)
!
      n_shp_in_basis_offset = 0
      n_cholesky_offset = 0
!
      do batch = 1, this%n_batches
!
!        Basis_shell_data file order:
!
!           1. number of shps in basis
!           2. basis_shell_info
!           3. cholesky_basis

         write(temp_name, '(a11, i4.4)') 'basis_info_', batch
         batch_file = sequential_file(trim(temp_name))
!
         call batch_file%open_('read', 'rewind')
!
         call mem%alloc(basis_shell_info, n_shp_in_basis_batches(batch), 4)
         call mem%alloc(cholesky_basis, n_cholesky_batches(batch), 3)
!
         call batch_file%read_blank()
         call batch_file%read_(basis_shell_info, n_shp_in_basis_batches(batch)*4)
         call batch_file%read_(cholesky_basis, n_cholesky_batches(batch)*3)
!
         call batch_file%close_()
!
         do J = 1, n_cholesky_batches(batch)
!
            alpha(n_cholesky_offset + J)      = cholesky_basis(J, 1)
            beta(n_cholesky_offset + J)       = cholesky_basis(J, 2)
            alpha_beta(n_cholesky_offset + J) = cholesky_basis(J, 3)
!
         enddo
!
         do shp = 1, n_shp_in_basis_batches(batch)
!
            A(n_shp_in_basis_offset + shp)                 = basis_shell_info(shp, 1)
            B(n_shp_in_basis_offset + shp)                 = basis_shell_info(shp, 2)
            AB(n_shp_in_basis_offset + shp)                = basis_shell_info(shp, 3)
            n_basis_aop_in_AB(n_shp_in_basis_offset + shp) = basis_shell_info(shp, 4)
!
         enddo
!
         n_shp_in_basis_offset = n_shp_in_basis_offset + n_shp_in_basis_batches(batch)
         n_cholesky_offset    = n_cholesky_offset + n_cholesky_batches(batch)
!
         call mem%dealloc(basis_shell_info, n_shp_in_basis_batches(batch), 4)
         call mem%dealloc(cholesky_basis, n_cholesky_batches(batch), 3)
!
      enddo
!
!     Sort the arrays according to an alphabeta and an AB ordering
!     from smallest to largest
!
      call mem%alloc(index_AB, n_shp_in_basis_total)
      call quicksort_with_index_ascending_int(AB, index_AB, n_shp_in_basis_total)
!
      call mem%alloc(index_alpha_beta, n_cholesky_total)
      call quicksort_with_index_ascending_int(alpha_beta, index_alpha_beta, n_cholesky_total)
!
      call mem%alloc(sorted_alpha, n_cholesky_total)
      call mem%alloc(sorted_beta, n_cholesky_total)
      call mem%alloc(sorted_alpha_beta, n_cholesky_total)
!
      sorted_alpha_beta = alpha_beta
      call mem%dealloc(alpha_beta, n_cholesky_total)
!
      do J = 1, n_cholesky_total
!
         sorted_alpha(J) = alpha(index_alpha_beta(J))
         sorted_beta(J)  = beta(index_alpha_beta(J))
!
      enddo
!
      call mem%dealloc(alpha, n_cholesky_total)
      call mem%dealloc(beta, n_cholesky_total)
!
      call mem%alloc(sorted_A, n_shp_in_basis_total)
      call mem%alloc(sorted_B, n_shp_in_basis_total)
      call mem%alloc(sorted_AB, n_shp_in_basis_total)
      call mem%alloc(sorted_n_basis_aop_in_AB, n_shp_in_basis_total)
!
      sorted_AB = AB
!
      call mem%dealloc(AB, n_shp_in_basis_total)
!
      do J = 1, n_shp_in_basis_total
!
         sorted_A(J) = A(index_AB(J))
         sorted_B(J) = B(index_AB(J))
         sorted_n_basis_aop_in_AB(J) = n_basis_aop_in_AB(index_AB(J))
!
      enddo
!
      call mem%dealloc(A, n_shp_in_basis_total)
      call mem%dealloc(B, n_shp_in_basis_total)
      call mem%dealloc(n_basis_aop_in_AB, n_shp_in_basis_total)
!
      call mem%dealloc(index_alpha_beta, n_cholesky_total)
      call mem%dealloc(index_AB, n_shp_in_basis_total)
!
!     Construct significant shell pair logical array,
!     and count the number of significant AO and shell pairs
!
      call mem%alloc(sig_shp, this%n_shp)
!
      sig_shp = .false.
!
      n_sig_shp = 0
      n_sig_aop = 0
!
      I = 1
      shp = 0
!
      do B_shell = 1, this%n_s
         do A_shell = B_shell, this%n_s
!
            shp = shp + 1
!
            if (shp == sorted_AB(I)) then
!
               A_range = ao%shells(A_shell)
               B_range = ao%shells(B_shell)
!
               sig_shp(shp) = .true.
               n_sig_shp = n_sig_shp + 1
               n_sig_aop = n_sig_aop + get_size_shp(A_range, B_range)
!
               I = I + 1
!
            endif
!
         enddo
      enddo
!
!     Read old diagonal from file, along with old sig_shp logical array, related info.,
!     and screening vector (old refers here to the initially screened diagonal)
!
      call mem%alloc(D, n_sig_aop)
      D = zero
      call mem%alloc(screening_vector, n_sig_aop)
!
      call this%diagonal_info_target%open_('read', 'rewind')
!
      call this%diagonal_info_target%read_(n_sig_shp_old)
      call this%diagonal_info_target%read_(n_sig_aop_old)
!
      call mem%alloc(sig_shp_old, this%n_shp)
      call mem%alloc(D_old, n_sig_aop_old)
      call mem%alloc(screening_vector_old, n_sig_aop_old)
!
      call this%diagonal_info_target%read_(sig_shp_old, this%n_shp)
      call this%diagonal_info_target%read_(D_old, n_sig_aop_old)
      call this%diagonal_info_target%read_(screening_vector_old, n_sig_aop_old)
!
      call this%diagonal_info_target%close_()
!
!     Copy the correct elements of the initial D into the new D using cholesky basis array.
!     We precalculate alpha beta offsets both old and new, then copy afterwards.
!
      call mem%alloc(alpha_beta_offset, n_sig_shp)
      call mem%alloc(alpha_beta_offset_old, n_sig_shp)
!
      count_sig = 0
      current_offset = 0
      current_offset_old = 0
!
      shp = 0
!
      do B_shell = 1, this%n_s
         do A_shell = B_shell, this%n_s
!
            shp = shp + 1
!
            if (sig_shp_old(shp)) then
!
               A_range = ao%shells(A_shell)
               B_range = ao%shells(B_shell)
!
               if (sig_shp(shp)) then
!
                  count_sig = count_sig + 1
!
                  alpha_beta_offset_old(count_sig) = current_offset_old
                  alpha_beta_offset(count_sig)     = current_offset
!
                  current_offset = current_offset + get_size_shp(A_range, B_range)
!
               endif
!
               current_offset_old = current_offset_old + get_size_shp(A_range, B_range)
!
            endif
!
         enddo
!
      enddo
!
      I = 0
      count_sig = 0
      n_basis_aop_in_AB_offset = 0
!
      do while (I .lt. n_shp_in_basis_total)
!
         I = I + 1
         count_sig = count_sig + 1
!
         n_basis_aop_in_AB_total = sorted_n_basis_aop_in_AB(I)
!
         A_range = ao%shells(sorted_A(I))
         B_range = ao%shells(sorted_B(I))
!
         do aop = 1, n_basis_aop_in_AB_total
!
            alpha_in_A = sorted_alpha(aop + n_basis_aop_in_AB_offset) - A_range%first + 1
            beta_in_B  = sorted_beta(aop + n_basis_aop_in_AB_offset) - B_range%first + 1
!
            if (sorted_A(I) == sorted_B(I)) then
!
               alpha_beta_in_AB = max(alpha_in_A, beta_in_B)*(max(alpha_in_A, beta_in_B) - 3)/2 + alpha_in_A + beta_in_B
!
            else
!
               alpha_beta_in_AB = A_range%length*(beta_in_B - 1) + alpha_in_A
!
            endif
!
            D(alpha_beta_in_AB + alpha_beta_offset(count_sig)) = &
                                    D_old(alpha_beta_in_AB + alpha_beta_offset_old(count_sig))
!
            screening_vector(alpha_beta_in_AB + alpha_beta_offset(count_sig)) = &
                                    screening_vector_old(alpha_beta_in_AB + alpha_beta_offset_old(count_sig))
!
         enddo
!
         n_basis_aop_in_AB_offset = n_basis_aop_in_AB_offset + n_basis_aop_in_AB_total
!
      enddo
!
      call mem%dealloc(alpha_beta_offset, n_sig_shp)
      call mem%dealloc(alpha_beta_offset_old, n_sig_shp)
!
      call mem%dealloc(D_old, n_sig_aop_old)
      call mem%dealloc(screening_vector_old, n_sig_aop_old)
!
      call mem%dealloc(sorted_alpha, n_cholesky_total)
      call mem%dealloc(sorted_beta, n_cholesky_total)
      call mem%dealloc(sorted_alpha_beta, n_cholesky_total)
!
      call mem%dealloc(sorted_A, n_shp_in_basis_total)
      call mem%dealloc(sorted_B, n_shp_in_basis_total)
      call mem%dealloc(sorted_AB, n_shp_in_basis_total)
      call mem%dealloc(sorted_n_basis_aop_in_AB, n_shp_in_basis_total)
!
      call mem%dealloc(sig_shp_old, this%n_shp)
!
!     Write info file for target diagonal containing
!
!        1. number of significant shell pairs, number of significant ao pairs
!        2. sig_shp - vector of logicals to describe which shell pairs are significant
!        3. D_xy = ( xy | xy ), the significant diagonal.
!        4. Screening vector
!
      call output%printf('n', 'Significant AO and shell pairs in final decomposition:', &
                         fs='(/t6, a)')
!
      call output%printf('n', 'Significant shell pairs: (i17)', ints=[n_sig_shp], fs='(t9,a)')
      call output%printf('n', 'Significant AO pairs:    (i17)', ints=[n_sig_aop], fs='(t9, a)')
!
!
      call this%diagonal_info_target%open_('write', 'rewind')
!
      call this%diagonal_info_target%write_(n_sig_shp)
      call this%diagonal_info_target%write_(n_sig_aop)
      call this%diagonal_info_target%write_(sig_shp, this%n_shp)
      call this%diagonal_info_target%write_(D, n_sig_aop)
      call this%diagonal_info_target%write_(screening_vector, n_sig_aop)
!
      call this%diagonal_info_target%close_()
!
      call mem%dealloc(sig_shp, this%n_shp)
!
      call mem%dealloc(D, n_sig_aop)
      call mem%dealloc(screening_vector, n_sig_aop)
!
   end subroutine construct_diagonal_from_batch_bases
!
!
   subroutine determine_cholesky_basis_PCD(this, ao)
!!
!!    Determine cholesky basis for partitioned Cholesky decomposition
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Determines the elements of the auxiliary basis by PCD
!!
      implicit none
!
      class(eri_cd), intent(inout) :: this
!
      class(ao_tool), intent(in) :: ao
!
      integer :: batch
!
      integer, dimension(:), allocatable :: n_cholesky_batches, n_shp_in_basis_batches
!
      type(sequential_file) :: batch_file_diag
      type(sequential_file) :: batch_file_basis
!
      character(len=100) :: temp_name
!
      call mem%alloc(n_cholesky_batches, this%n_batches)
      call mem%alloc(n_shp_in_basis_batches, this%n_batches)
!
      n_cholesky_batches = 0
      n_shp_in_basis_batches = 0
!
      call this%construct_diagonal_batches(ao)
!
      call output%printf('n', '- Decomposing batched diagonal', fs='(/t3,a)')
!
      do batch = 1, this%n_batches
!
         call output%printf('n', 'Batch (i0):', ints=[batch], fs='(/t3,a)')
!
         write(temp_name, '(a14, i4.4)')'diagonal_info_', batch
         batch_file_diag = sequential_file(trim(temp_name))

         write(temp_name, '(a11, i4.4)')'basis_info_', batch
         batch_file_basis = sequential_file(trim(temp_name))
!
         call this%determine_cholesky_basis_standard(ao, batch_file_diag, batch_file_basis)
!
         n_cholesky_batches(batch)     = this%n_cholesky
         n_shp_in_basis_batches(batch)  = this%n_shp_in_basis
!
      enddo
!
      call output%printf('n', '- Final decomposition step:', fs='(/t3,a)')
!
      call this%construct_diagonal_from_batch_bases(ao, n_cholesky_batches, n_shp_in_basis_batches)
      call this%determine_cholesky_basis_standard(ao, this%diagonal_info_target, this%cholesky_basis_file)
!
      call mem%dealloc(n_cholesky_batches, this%n_batches)
      call mem%dealloc(n_shp_in_basis_batches, this%n_batches)
!
   end subroutine determine_cholesky_basis_PCD
!
!
   subroutine determine_cholesky_basis_standard(this, ao, diagonal_info, basis_info)
!!
!!    Determine cholesky basis (standard)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Determines the elements of the Cholesky auxiliary basis.
!!
      use array_utilities, only: quicksort_with_index_descending
      use array_utilities, only: get_n_highest, is_significant
      use array_utilities, only: reduce_array_int, reduce_vector
!
      implicit none
!
      class(eri_cd), intent(inout) :: this
!
      class(ao_tool), intent(in) :: ao
!
      type(sequential_file), intent(inout) :: diagonal_info
      type(sequential_file), intent(inout) :: basis_info
!
!     Local variables
!
!     Integers
!
      integer :: n_sig_shp, n_sig_aop
      integer :: n_new_sig_shp, n_new_sig_aop, current_new_sig_shp
      integer :: n_qual_shp, n_qual_aop, n_previous_qual_aop, n_qual_aop_in_shp
      integer :: shp, current_sig_shp
      integer :: first_sig_aop, last_sig_aop, aop
      integer :: A, B, AB, AB_shp
      integer :: C, D, CD_shp
      integer :: I, J
      integer :: w, x, y, z
      integer :: xy, xy_packed, xy_max, wx, wx_packed
      integer :: sig_neg
      integer :: first, last
      integer :: first_x, first_y
      integer :: current_qual, qual
      integer :: n_cholesky_in_node, n_new_cholesky
      integer :: n_shp_in_basis
      integer :: sig_shp_counter
      integer :: shp_in_basis
!
!     Integer allocatable arrays
!
!     Map significant shell pair to first ao pair
      integer, dimension(:), allocatable :: sig_shp_to_first_sig_aop
!     Map significant shell pair to first ao pair
      integer, dimension(:), allocatable :: new_sig_shp_to_first_sig_aop
!     Index array for sorting shell pairs according to their maximum values
      integer, dimension(:), allocatable :: sorted_max_sig_shp
!     Index array for sorting the qualified ao pairs in shell pair
      integer, dimension(:), allocatable :: sorted_qual_aop_in_shp_indices
!     Offsets for omp-loop, number of qualified ao pairs in preceding shell pair
      integer, dimension(:), allocatable :: n_qual_aop_in_prev_shps
!     Index list containing order in which qualified diagonals are sig in decomposition
      integer, dimension(:), allocatable :: qual_max
!     Map significant shell pair indices to significant shell pair indices of last iteration
      integer, dimension(:), allocatable :: sig_shp_to_previous_sig_shp
!
!     Maps significant shell pair to shells
      integer, dimension(:,:), allocatable :: sig_shp_to_shells
!     Maps significant shell pair to shells
      integer, dimension(:,:), allocatable :: new_sig_shp_to_shells
!     Maps significant ao pair to aos
      integer, dimension(:,:), allocatable :: sig_aop_to_aos
!     Maps significant ao pair to aos
      integer, dimension(:,:), allocatable :: new_sig_aop_to_aos
!     List of qualified shell pairs
      integer, dimension(:,:), allocatable :: qual_shp
!     List of qualified shell pairs, copy used to reduce size
      integer, dimension(:,:), allocatable :: qual_shp_copy
!     List of qualified ao pairs
      integer, dimension(:,:), allocatable :: qual_aop
!     List of qualified ao pairs, copy used to reduce size
      integer, dimension(:,:), allocatable :: qual_aop_copy
!     ao and ao pair indices of the elements of the cholesky basis
      integer, dimension(:,:), allocatable :: cholesky_basis
!     ao and ao pair indices of the elements of the cholesky basis, written to file at the end
      integer, dimension(:,:), allocatable :: cholesky_basis_new
!     Info on shells containing elements of the basis
      integer, dimension(:,:), allocatable :: basis_shell_info_full
!     Info on shells containing elements of the basis, written to file at end of routine
      integer, dimension(:,:), allocatable :: basis_shell_info
!
!     Logicals
!
      logical :: done, write_warning, construct_more_choleskys, found
!
!     Logical allocatable arrays
!
      logical, dimension(:), allocatable :: sig_shp
      logical, dimension(:), allocatable :: new_sig_shp
!
!     Reals
!
      real(dp) :: D_max_full, D_max
!
!     Reals for timings
!
      real(dp) :: s_select_basis_time, e_select_basis_time
      real(dp) :: full_reduce_time, s_reduce_time, e_reduce_time
      real(dp) :: full_construct_time, s_construct_time, e_construct_time
!
!     Real allocatable arrays
!
!     Array for eri diagonal elements
      real(dp), dimension(:), allocatable :: D_xy
!     Array for eri diagonal elements, used for reduction
      real(dp), dimension(:), allocatable :: D_xy_new
!     Array for accumulating approximate diagonal
      real(dp), dimension(:), allocatable :: approx_diagonal_accumulative
!     Maximum in each significant shell pair
      real(dp), dimension(:), allocatable :: max_in_sig_shp
!     Sorted qualified ao pair in shell pair
      real(dp), dimension(:), allocatable :: sorted_qual_aop_in_shp
!     Screening vector for diagonal
      real(dp), dimension(:), allocatable :: screening_vector
!     Screening vector for diagonal, used for reduction
      real(dp), dimension(:), allocatable :: screening_vector_new
!
!     Array for eri
      real(dp), dimension(:,:), allocatable :: g_wxyz
!     Array used for dgemm, reordered copy of cholesky vectors of current batch of qualified
      real(dp), dimension(:,:), allocatable :: cholesky_tmp
!
!     Array for eri for shell pairs AB and CD
      real(dp), dimension(ao%max_sh_size**4), target :: g_ABCD
!     Pointer for eri for shell pairs AB and CD
      real(dp), dimension(:,:,:,:), pointer :: g_ABCD_p
!
!     Real pointers
!
      real(dp), dimension(:,:), pointer     :: cholesky
      real(dp), dimension(:,:), pointer     :: cholesky_new
!
!     Intervals
!
      type(range_) :: A_range, B_range, C_range, D_range
!
!     Cholesky linked list
!
      type(cholesky_array_list) :: cholesky_array
!
      call cpu_time(s_select_basis_time)

!     Read diagonal info file containing (name given as argument)
!
!        1. number of significant shell pairs, number of significant ao pairs
!        2. sig_shp - vector of logicals to describe which shell pairs are significant
!        3. D_xy = ( xy | xy ), the significant diagonal.
!
      call diagonal_info%open_('read','rewind')
!
      call diagonal_info%read_(n_sig_shp)
      call diagonal_info%read_(n_sig_aop)
!
      call mem%alloc(sig_shp, this%n_shp)
      call mem%alloc(D_xy, n_sig_aop)
      call mem%alloc(screening_vector, n_sig_aop)
!
      call diagonal_info%read_(sig_shp, this%n_shp)
      call diagonal_info%read_(D_xy, n_sig_aop)
      call diagonal_info%read_(screening_vector, n_sig_aop)
!
      call diagonal_info%close_()
!
!     Construct info arrays
!
      call mem%alloc(sig_shp_to_first_sig_aop, n_sig_shp + 1) ! Maps significant shell pair to first ao pair
      sig_shp_to_first_sig_aop = 0
!
!     Note: allocated with length n_significant_shp + 1, last element is used for n_sig_aop + 1
!     This is convenient because sig_shp_to_first_sig_aop will be used to calculate lengths.
!
      sig_shp_to_first_sig_aop(n_sig_shp + 1) = n_sig_aop + 1
!
      call mem%alloc(sig_shp_to_shells, n_sig_shp, 2) ! [A, B]
      sig_shp_to_shells = 0
!
      call mem%alloc(sig_aop_to_aos, n_sig_aop, 2) ! [alpha, beta]
      sig_aop_to_aos = 0
!
      shp              = 1
      current_sig_shp  = 1
      first_sig_aop    = 1
!
      do B = 1, this%n_s
!
         do A = B, this%n_s
!
            if (sig_shp(shp)) then
!
               sig_shp_to_first_sig_aop(current_sig_shp) = first_sig_aop
!
               A_range = ao%shells(A)
               B_range = ao%shells(B)
!
               sig_shp_to_shells(current_sig_shp, 1) = A
               sig_shp_to_shells(current_sig_shp, 2) = B
!
               if (A .eq. B) then
!
                  do x = 1, A_range%length
                     do y = 1, B_range%length
!
                        xy = A_range%length*(y - 1) + x
                        xy_packed = (max(x,y)*(max(x,y)-3)/2) + x + y
!
                        sig_aop_to_aos(xy_packed + first_sig_aop - 1, 1) = &
                                                      A_range%first + x - 1
!
                        sig_aop_to_aos(xy_packed + first_sig_aop - 1, 2) = &
                                                      B_range%first + y - 1
!
                     enddo
                  enddo
!
               else ! A ≠ B
!
                  do x = 1, (A_range%length)
                     do y = 1, (B_range%length)
!
                        xy = A_range%length*(y - 1) + x
!
                        sig_aop_to_aos(xy + first_sig_aop - 1, 1) = A_range%first + x - 1
                        sig_aop_to_aos(xy + first_sig_aop - 1, 2) = B_range%first + y - 1
!
                     enddo
                  enddo
!
               endif
!
               first_sig_aop = first_sig_aop + get_size_shp(A_range, B_range)
!
               current_sig_shp = current_sig_shp + 1
!
            endif ! End of if (significant)
!
            shp = shp + 1
!
         enddo
      enddo
!
!     Determining the basis
!
      call output%printf('n', 'Iter.  #Sign. ao pairs / shell pairs   Max &
                         &diagonal    #Qualified    #Cholesky    Cholesky array size', &
                         ll=103, fs='(/t3, a)')
      call output%print_separator('n', 103,'-', fs='(t3,a)')
!
      this%iteration   = 0
      this%n_cholesky  = 0
!
      sig_neg = 0 ! Counts number of significant neative diagonals (absolute value > 10^-10)
!
      full_reduce_time     = zero
      full_construct_time  = zero
!
      done = .false.
!
      cholesky_array = cholesky_array_list()
!
      do while (.not. done)
!
         write_warning = .true. ! Logical used to ensure warning of
!                                 significant negative diagonal only appears once per iteration
!
         this%iteration = this%iteration + 1
!
!        Shell maximums and shell maximums indices vectors
!
         call mem%alloc(max_in_sig_shp, n_sig_shp)
!
         max_in_sig_shp = zero
!
         do shp = 1, n_sig_shp
!
!           Get first and last indices of shell pair
!
            first = sig_shp_to_first_sig_aop(shp)
            last  = sig_shp_to_first_sig_aop(shp + 1) - 1
!
!           Determine the largest elements
!
            do I = first, last
!
               if (D_xy(I) .gt. max_in_sig_shp(shp)) then
!
                  max_in_sig_shp(shp) = D_xy(I)
!
               endif
!
            enddo
!
         enddo
!
!        Sort from largest to smallest and determine an index array of sorting
!
         call mem%alloc(sorted_max_sig_shp, n_sig_shp)
         sorted_max_sig_shp = 0
!
         call quicksort_with_index_descending(max_in_sig_shp, sorted_max_sig_shp, n_sig_shp)
!
         D_max_full  = max_in_sig_shp(1)
         n_qual_aop  = 0
         n_qual_shp   = 0
!
         call mem%dealloc(max_in_sig_shp, n_sig_shp)
!
         call mem%alloc(qual_aop, this%max_qual, 3)
         call mem%alloc(qual_shp, this%n_shp, 3)
         qual_shp = 0
         qual_aop = 0
!
!        Determine qualified shell pairs
!
         do shp = 1, n_sig_shp
!
            current_sig_shp = sorted_max_sig_shp(shp)
!
            first_sig_aop = sig_shp_to_first_sig_aop(current_sig_shp)
            last_sig_aop  = sig_shp_to_first_sig_aop(current_sig_shp + 1) - 1
!
            n_qual_aop_in_shp = 0
!
            do aop = first_sig_aop, last_sig_aop
!
               if ((D_xy(aop) .ge. this%span*D_max_full) .and. (n_qual_aop .lt. this%max_qual)) then
!
                  n_qual_aop_in_shp  = n_qual_aop_in_shp + 1
                  n_qual_aop        = n_qual_aop + 1
!
               endif
!
            enddo
!
            if (n_qual_aop_in_shp .ne. 0) then
!
               n_qual_shp = n_qual_shp + 1
!
               call mem%alloc(sorted_qual_aop_in_shp_indices, n_qual_aop_in_shp)
               call mem%alloc(sorted_qual_aop_in_shp, n_qual_aop_in_shp)
!
               call get_n_highest(n_qual_aop_in_shp, last_sig_aop - first_sig_aop + 1, &
                                 D_xy(first_sig_aop:last_sig_aop), sorted_qual_aop_in_shp, &
                                 sorted_qual_aop_in_shp_indices)
!
               n_previous_qual_aop = (n_qual_aop - n_qual_aop_in_shp)
!
               do aop = 1, n_qual_aop_in_shp
!
                  qual_aop(aop + n_previous_qual_aop, 1) = sig_aop_to_aos(sorted_qual_aop_in_shp_indices(aop) &
                                                               + first_sig_aop - 1, 1)
!
                  qual_aop(aop + n_previous_qual_aop, 2) = sig_aop_to_aos(sorted_qual_aop_in_shp_indices(aop) &
                                                               + first_sig_aop - 1, 2)
!
                  qual_aop(aop + n_previous_qual_aop, 3) = sorted_qual_aop_in_shp_indices(aop) &
                                                               + first_sig_aop - 1
!
               enddo
!
               first_x = sig_aop_to_aos(first_sig_aop, 1) ! alpha
               first_y = sig_aop_to_aos(first_sig_aop, 2) ! beta
!
               qual_shp(n_qual_shp, 1) = ao%ao_to_shell(first_x)
               qual_shp(n_qual_shp, 2) = ao%ao_to_shell(first_y)
               qual_shp(n_qual_shp, 3) = n_qual_aop_in_shp
!
               call mem%dealloc(sorted_qual_aop_in_shp_indices, n_qual_aop_in_shp)
               call mem%dealloc(sorted_qual_aop_in_shp, n_qual_aop_in_shp)
!
            endif
!
            if (n_qual_aop == this%max_qual) then
!
               exit
!
            endif
!
         enddo
!
         call mem%dealloc(sorted_max_sig_shp, n_sig_shp)
!
!        Cut out the qualified parts of the aop and shp lists
!
         call mem%alloc(qual_aop_copy, n_qual_aop, 3)
         call mem%alloc(qual_shp_copy, n_qual_shp, 3)
!
         qual_aop_copy(:, :) = qual_aop(1 : n_qual_aop, :)
         qual_shp_copy(:, :)  = qual_shp(1 : n_qual_shp, :)
!
         call mem%dealloc(qual_aop, this%max_qual, 3)
         call mem%dealloc(qual_shp, this%n_shp, 3)
!
         call mem%alloc(qual_aop, n_qual_aop, 3)
         call mem%alloc(qual_shp, n_qual_shp, 3)
!
         qual_aop    = qual_aop_copy
         qual_shp     = qual_shp_copy
!
         call mem%dealloc(qual_aop_copy, n_qual_aop, 3)
         call mem%dealloc(qual_shp_copy, n_qual_shp, 3)
!
!        Prepare to construct g_wxyz in parallelized loop
!
         call mem%alloc(n_qual_aop_in_prev_shps, n_qual_shp)
         n_qual_aop_in_prev_shps = 0
!
         do CD_shp = 1, n_qual_shp - 1
!
             n_qual_aop_in_prev_shps(CD_shp + 1) = n_qual_aop_in_prev_shps(CD_shp) + qual_shp(CD_shp, 3)
!
         enddo
!
!        Construct g_wxyz
!
         call mem%alloc(g_wxyz, n_sig_aop, n_qual_aop)
!
!$omp parallel do &
!$omp private(AB_shp, CD_shp, A, B, A_range, B_range, C, D, C_range, D_range, &
!$omp  aop, w, x, y, z, wx, wx_packed, g_ABCD, g_ABCD_p, n_qual_aop_in_shp) &
!$omp shared(g_wxyz, n_qual_aop_in_prev_shps, qual_aop) &
!$omp schedule(guided)
         do CD_shp = 1, n_qual_shp
!
            C                = qual_shp(CD_shp, 1)
            D                = qual_shp(CD_shp, 2)
            n_qual_aop_in_shp = qual_shp(CD_shp, 3)
!
            C_range = ao%shells(C)
            D_range = ao%shells(D)
!
!           Calculate the ({wx} | J) integrals,
!           where {wx} is the screened list of integrals
!
            do AB_shp = 1, n_sig_shp
!
               A = sig_shp_to_shells(AB_shp, 1)
               B = sig_shp_to_shells(AB_shp, 2)
!
               A_range = ao%shells(A)
               B_range = ao%shells(B)
!
               call ao%get_eri(g_ABCD, A, B, C, D)
!
               g_ABCD_p(1 : A_range%length, 1 : B_range%length, &
                        1 : C_range%length, 1 : D_range%length) &
                        => g_ABCD(1 : (A_range%length)*(B_range%length)*(C_range%length)*(D_range%length))
!
               do aop = 1, n_qual_aop_in_shp
!
                  y = qual_aop(aop + n_qual_aop_in_prev_shps(CD_shp), 1)
                  z = qual_aop(aop + n_qual_aop_in_prev_shps(CD_shp), 2)
!
                  if (A == B) then
!
                     do w = 1, A_range%length
                        do x = w, B_range%length
!
                           wx_packed = (max(w,x)*(max(w,x)-3)/2) + w + x
!
                           g_wxyz(sig_shp_to_first_sig_aop(AB_shp) + wx_packed - 1, aop + n_qual_aop_in_prev_shps(CD_shp)) &
                                   = g_ABCD_p(w, x, y - C_range%first + 1, z - D_range%first + 1)
!
                        enddo
                     enddo
!
                  else
!
                     do w = 1, A_range%length
                        do x = 1, B_range%length
!
                           wx = A_range%length*(x-1) + w
!
                           g_wxyz(sig_shp_to_first_sig_aop(AB_shp) + wx - 1, aop + n_qual_aop_in_prev_shps(CD_shp)) &
                                   = g_ABCD_p(w, x, y - C_range%first + 1, z - D_range%first + 1)
!
                        enddo
                     enddo
!
                  endif
!
               enddo
!
         enddo
!
      enddo ! cd_shp
!$omp end parallel do
!
         call mem%dealloc(n_qual_aop_in_prev_shps, n_qual_shp)
!
!        Subtract old cholesky vectors
!
         call cpu_time(s_construct_time)
!
         if (this%n_cholesky .ne. 0) then
!
            do I = 1, cholesky_array%n_nodes
!
               call cholesky_array%get_array(cholesky, I) ! Let cholesky point to node # I
!
               n_cholesky_in_node = cholesky_array%get_n_columns_element(I)
!
               call mem%alloc(cholesky_tmp, n_cholesky_in_node, n_qual_aop)
!
               do J = 1, n_qual_aop
!
                  cholesky_tmp(:, J) = cholesky(qual_aop(J, 3), :)
!
               enddo
!
               call dgemm('N', 'N',            &
                           n_sig_aop,          &
                           n_qual_aop,         &
                           n_cholesky_in_node, &
                           -one,               &
                           cholesky,           &
                           n_sig_aop,          &
                           cholesky_tmp,       &
                           n_cholesky_in_node, &
                           one,                &
                           g_wxyz,             &
                           n_sig_aop)
!
               call mem%dealloc(cholesky_tmp, n_cholesky_in_node, n_qual_aop)
!
            enddo
!
            call mem%alloc(cholesky_basis, this%n_cholesky + n_qual_aop, 3)
            cholesky_basis(1 : this%n_cholesky, :) = cholesky_basis_new(:, :)
            call mem%dealloc(cholesky_basis_new, this%n_cholesky, 3)
!
         else
!
            call mem%alloc(cholesky_basis, n_qual_aop, 3)
            cholesky_basis = 0
!
         endif
!
         call mem%alloc(approx_diagonal_accumulative, n_sig_aop)
         approx_diagonal_accumulative = zero
!
         call cholesky_array%push_back(n_sig_aop, n_qual_aop) ! Make new element at the tail
         call cholesky_array%get_array(cholesky_new, cholesky_array%n_nodes) ! Point cholesky_new to this last element
         cholesky_new = zero
!
         current_qual = 0
!
         construct_more_choleskys = .true.
!
         call mem%alloc(qual_max, n_qual_aop)
!
         do while ((current_qual .lt. n_qual_aop) .and. construct_more_choleskys)
!
            current_qual = current_qual + 1
!
            qual_max(current_qual) = 1
            xy_max = qual_aop(1, 3)
            D_max = D_xy(xy_max) - approx_diagonal_accumulative(xy_max)
!
            do qual = 1, n_qual_aop
!
               xy = qual_aop(qual, 3)
!
               if (D_xy(xy) - approx_diagonal_accumulative(xy) .gt. D_max) then
!
                  qual_max(current_qual) = qual
                  xy_max = xy
                  D_max = D_xy(xy) - approx_diagonal_accumulative(xy)
!
               endif
!
            enddo
!
            if ((D_max*screening_vector(xy_max) .gt. this%threshold) .and. &
               (D_max .gt. this%threshold)  .and. &
               (D_max .ge. this%span*D_max_full)) then
!
               cholesky_basis(this%n_cholesky + current_qual, 1) = qual_aop(qual_max(current_qual), 1)
               cholesky_basis(this%n_cholesky + current_qual, 2) = qual_aop(qual_max(current_qual), 2)
!
               A = ao%ao_to_shell(qual_aop(qual_max(current_qual), 1))
               B = ao%ao_to_shell(qual_aop(qual_max(current_qual), 2))
!
               cholesky_basis(this%n_cholesky + current_qual, 3) = get_shp_from_shells(A, B, this%n_s)
!
               cholesky_new(: , current_qual) = g_wxyz(:, qual_max(current_qual))
!
               if (current_qual .gt. 1) then
!
                  call mem%alloc(cholesky_tmp, 1, current_qual - 1)
!
                  cholesky_tmp(1, :) = cholesky_new(qual_aop(qual_max(current_qual), 3), 1 : current_qual - 1)
!
                  call dgemm('N', 'T',                               &
                           n_sig_aop,                                &
                           1,                                        &
                           current_qual - 1,                         &
                           -one,                                     &
                           cholesky_new,                             &
                           n_sig_aop,                                &
                           cholesky_tmp,                             &
                           1,                                        &
                           one,                                      &
                           cholesky_new(1:n_sig_aop, current_qual),  &
                           n_sig_aop)
!
                  call mem%dealloc(cholesky_tmp, 1, current_qual - 1)
!
               endif
!
               call dscal(n_sig_aop, one/sqrt(D_max), cholesky_new(1, current_qual), 1)
!
!$omp parallel do private(xy)
               do xy = 1, n_sig_aop
!
                  if (D_xy(xy) == zero) then
!
                     cholesky_new(xy, current_qual) = zero
!
                  endif
!
               enddo
!$omp end parallel do
!
!$omp parallel do private(xy)
               do xy = 1, n_sig_aop
!
                  approx_diagonal_accumulative(xy) = approx_diagonal_accumulative(xy) &
                                                   + cholesky_new(xy, current_qual)**2
!
               enddo
!$omp end parallel do
!
               D_xy(qual_aop(qual_max(current_qual), 3)) = zero
               approx_diagonal_accumulative(qual_aop(qual_max(current_qual), 3))  = zero
!
               do xy = 1, n_sig_aop
!
                  if (D_xy(xy) - approx_diagonal_accumulative(xy) .lt. zero) then
!
                     if (abs(D_xy(xy) - approx_diagonal_accumulative(xy)) .gt. 1.0d-10) then
                        if (write_warning) then
                           call output%warning_msg('Found significant negative diagonal!')
                           write_warning = .false.
                        endif
                        sig_neg = sig_neg + 1
                     endif
!
                     D_xy(xy) = zero
                     approx_diagonal_accumulative(xy) = zero
!
                  elseif ((D_xy(xy) - approx_diagonal_accumulative(xy))*screening_vector(xy) &
                          .lt. this%threshold .or. &
                         (D_xy(xy) - approx_diagonal_accumulative(xy)) .lt. this%threshold) then
!
                     D_xy(xy) = zero
                     approx_diagonal_accumulative(xy) = zero
!
                  endif
!
               enddo
!
            else
!
               construct_more_choleskys = .false.
               current_qual = current_qual - 1
!
            endif
!
         enddo ! End of decomposition
!
!$omp parallel do private(xy)
         do xy = 1, n_sig_aop
!
            D_xy(xy) = D_xy(xy) - approx_diagonal_accumulative(xy)
!
         enddo
!$omp end parallel do
!
         call mem%dealloc(approx_diagonal_accumulative, n_sig_aop)
!
         call mem%dealloc(qual_max, n_qual_aop)
!
         call cpu_time(e_construct_time)
         full_construct_time = full_construct_time + e_construct_time - s_construct_time
!
         call mem%dealloc(g_wxyz, n_sig_aop, n_qual_aop)
!
         n_new_cholesky = current_qual
!
!        Find new significant diagonals
!
         n_new_sig_shp = 0
         call mem%alloc(new_sig_shp, n_sig_shp)
!
         new_sig_shp = .false.
!
         sig_shp_counter = 0
!
         do shp = 1, this%n_shp
!
            if (sig_shp(shp)) then
!
               sig_shp_counter = sig_shp_counter + 1
!
               first = sig_shp_to_first_sig_aop(sig_shp_counter)
               last  = sig_shp_to_first_sig_aop(sig_shp_counter + 1) - 1
!
               new_sig_shp(sig_shp_counter) = (is_significant(D_xy(first:last), &
                                                last - first + 1, this%threshold, &
                                                screening_vector(first:last) ) .and. &
                                             is_significant(D_xy(first:last), &
                                                last - first + 1, this%threshold ))
!
               sig_shp(shp) = new_sig_shp(sig_shp_counter)
!
               if (new_sig_shp(sig_shp_counter)) then
!
                  n_new_sig_shp = n_new_sig_shp + 1
!
               endif
!
            endif
!
         enddo
!
         call cpu_time(s_reduce_time)
!
         if (n_new_sig_shp .gt. 0) then
!
!           Update index lists: shps -> aops, aops -> aos, and shps -> full shps
!
            call mem%alloc(new_sig_shp_to_first_sig_aop, n_new_sig_shp + 1)
            new_sig_shp_to_first_sig_aop = 0
!
            call mem%alloc(sig_shp_to_previous_sig_shp, n_sig_shp + 1) ! 1 2 3 4 ... n_sig_shp, n_sig_shp + 1
            sig_shp_to_previous_sig_shp(n_sig_shp + 1) = n_sig_shp + 1
!
            current_new_sig_shp    = 1
            n_new_sig_aop = 0
            first_sig_aop = 1
!
            do shp = 1, n_sig_shp
!
               sig_shp_to_previous_sig_shp(shp) = shp
!
               if (new_sig_shp(shp)) then
!
                  A = sig_shp_to_shells(shp, 1)
                  B = sig_shp_to_shells(shp, 2)
!
                  A_range = ao%shells(A)
                  B_range = ao%shells(B)
!
                  new_sig_shp_to_first_sig_aop(current_new_sig_shp) = first_sig_aop
!
                  first_sig_aop = first_sig_aop + get_size_shp(A_range, B_range)
                  n_new_sig_aop = first_sig_aop - 1
!
                  current_new_sig_shp    = current_new_sig_shp + 1

               endif
!
            enddo
!
            new_sig_shp_to_first_sig_aop(current_new_sig_shp) = n_new_sig_aop + 1
!
            call mem%alloc(new_sig_aop_to_aos, n_new_sig_aop, 2)
!
            call reduce_array_int(sig_aop_to_aos,        &
                               new_sig_aop_to_aos,       &
                               sig_shp_to_first_sig_aop, &
                               new_sig_shp,              &
                               n_sig_shp,                &
                               n_sig_aop,                &
                               n_new_sig_aop,            &
                               2)
!
            call mem%alloc(new_sig_shp_to_shells, n_new_sig_shp, 2)
            new_sig_shp_to_shells = 0
!
            call reduce_array_int(sig_shp_to_shells,           &
                                  new_sig_shp_to_shells,       &
                                  sig_shp_to_previous_sig_shp, &
                                  new_sig_shp,                 &
                                  n_sig_shp,                   &
                                  n_sig_shp,                   &
                                  n_new_sig_shp,               &
                                  2)
!
            call mem%dealloc(sig_shp_to_previous_sig_shp, n_sig_shp + 1)
            call mem%dealloc(sig_shp_to_shells, n_sig_shp, 2)
            call mem%alloc(sig_shp_to_shells, n_new_sig_shp, 2)
!
            sig_shp_to_shells = new_sig_shp_to_shells
            call mem%dealloc(new_sig_shp_to_shells, n_new_sig_shp, 2)
!
            call mem%alloc(D_xy_new, n_new_sig_aop)
!
           call reduce_vector(D_xy,                      &
                             D_xy_new,                   &
                             sig_shp_to_first_sig_aop,   &
                             new_sig_shp,                &
                             n_sig_shp,                  &
                             n_sig_aop,                  &
                             n_new_sig_aop)
!
            call mem%dealloc(D_xy, n_sig_aop)
            call mem%alloc(D_xy, n_new_sig_aop)
!
            call dcopy(n_new_sig_aop, D_xy_new, 1, D_xy, 1)
!
            call mem%dealloc(D_xy_new, n_new_sig_aop)
!
            call mem%alloc(screening_vector_new, n_new_sig_aop)
!
           call reduce_vector(screening_vector,          &
                             screening_vector_new,       &
                             sig_shp_to_first_sig_aop,   &
                             new_sig_shp,                &
                             n_sig_shp,                  &
                             n_sig_aop,                  &
                             n_new_sig_aop)
!
            call mem%dealloc(screening_vector, n_sig_aop)
            call mem%alloc(screening_vector, n_new_sig_aop)
!
            call dcopy(n_new_sig_aop, screening_vector_new, 1, screening_vector, 1)
!
            call mem%dealloc(screening_vector_new, n_new_sig_aop)
!
!           Remove the unused columns of cholesky new
!
            call cholesky_array%keep_columns(cholesky_array%n_nodes, 1, n_new_cholesky)
!
            call cholesky_array%reduce(sig_shp_to_first_sig_aop,  &
                                       new_sig_shp,               &
                                       n_sig_shp,                 &
                                       n_new_sig_aop)
!
            cholesky_new => null()
!
            call mem%alloc(cholesky_basis_new, this%n_cholesky + n_new_cholesky, 3)
            cholesky_basis_new(:, :) = cholesky_basis(1 : this%n_cholesky + n_new_cholesky, :)
            call mem%dealloc(cholesky_basis, this%n_cholesky + n_qual_aop, 3)
!
!           Deallocate old lists & reallocate + copy over new lists
!
            call mem%dealloc(new_sig_shp, (n_sig_shp))
!
            call mem%dealloc(sig_shp_to_first_sig_aop, n_sig_shp + 1)
            call mem%alloc(sig_shp_to_first_sig_aop, n_new_sig_shp + 1)
            sig_shp_to_first_sig_aop = new_sig_shp_to_first_sig_aop
            call mem%dealloc(new_sig_shp_to_first_sig_aop, n_new_sig_shp + 1)
!
            call mem%dealloc(sig_aop_to_aos, n_sig_aop, 2)
            call mem%alloc(sig_aop_to_aos, n_new_sig_aop, 2)
            sig_aop_to_aos = new_sig_aop_to_aos
            call mem%dealloc(new_sig_aop_to_aos, n_new_sig_aop, 2)
!
            n_sig_shp = n_new_sig_shp
            n_sig_aop = n_new_sig_aop
!
            this%n_cholesky = this%n_cholesky + n_new_cholesky
!
            call mem%dealloc(qual_aop, n_qual_aop, 3)
            call mem%dealloc(qual_shp, n_qual_shp, 3)
!
            call output%printf('n', '(i4)        (i10) /(i8)      (e12.5)    &
                               &(i8)        (i7)        (i10)', &
                               ints=[this%iteration, n_sig_aop, n_sig_shp, &
                               n_qual_aop, this%n_cholesky, &
                               this%n_cholesky*n_sig_aop], &
                               reals=[D_max_full], ll=103, fs='(t3,a)')
!
         else
!
            call cholesky_array%finalize()
!
            cholesky_new   => null()
            cholesky       => null()
!
            call mem%dealloc(D_xy, n_sig_aop)
!
            call mem%alloc(cholesky_basis_new, this%n_cholesky + n_new_cholesky, 3)
            cholesky_basis_new(:, :) = cholesky_basis(1 : this%n_cholesky + n_new_cholesky, :)
            call mem%dealloc(cholesky_basis, this%n_cholesky + n_qual_aop, 3)
            call mem%dealloc(screening_vector, n_sig_aop)
!
            call mem%dealloc(sig_shp_to_first_sig_aop, n_sig_shp + 1)
            call mem%dealloc(sig_shp_to_shells, n_sig_shp, 2)
            call mem%dealloc(sig_aop_to_aos, n_sig_aop, 2)
            call mem%dealloc(new_sig_shp, (n_sig_shp))
!
            call mem%dealloc(qual_aop, n_qual_aop, 3)
            call mem%dealloc(qual_shp, n_qual_shp, 3)
!
            this%n_cholesky = this%n_cholesky + n_new_cholesky
!
            done = .true.
!
            call output%printf('n', '(i4)        (i10) /(i8)      (e12.5)    &
                               &(i8)        (i7)        (i10)', &
                               ints=[this%iteration, 0, 0, n_qual_aop, &
                               this%n_cholesky, 0], reals=[D_max_full], &
                               ll=103, fs='(t3,a)')
!
!
         endif
!
         call cpu_time(e_reduce_time)
         full_reduce_time = full_reduce_time + e_reduce_time - s_reduce_time
!
      enddo ! while not done
!
      call output%print_separator('n', 103,'-', fs='(t3,a)')
!
!     Timings
!
      call cpu_time(e_select_basis_time)
!
      if (sig_neg .gt. 0) &
         call output%printf('n', 'Number of significant negative diagonals: (i0)', &
                             ints=[sig_neg], fs='(/t6,a)')
!
!     Prepare info on basis
!
!     Construct a list of all shell pairs (and shells) that contain elements of the basis
!     and how many elements of the basis they contain
!
      call mem%alloc(basis_shell_info_full, this%n_shp, 4) ! A, B, AB, n_basis_aops_in_shp
      basis_shell_info_full = 0
!
      n_shp_in_basis = 0
!
      do i = 1, this%n_cholesky
!
         A = ao%ao_to_shell(cholesky_basis_new(i, 1))
         B = ao%ao_to_shell(cholesky_basis_new(i, 2))
!
         AB = get_shp_from_shells(A, B, this%n_s)
!
         found = .false.
!
         do shp_in_basis = 1, n_shp_in_basis
!
            if (AB == basis_shell_info_full(shp_in_basis, 3)) then
               found = .true.
               basis_shell_info_full(shp_in_basis, 4) = basis_shell_info_full(shp_in_basis, 4) + 1
               exit
            endif
!
         enddo
!
         if(.not. found) then
!
            n_shp_in_basis = n_shp_in_basis + 1
!
            basis_shell_info_full(n_shp_in_basis, 1) = A
            basis_shell_info_full(n_shp_in_basis, 2) = B
            basis_shell_info_full(n_shp_in_basis, 3) = AB
            basis_shell_info_full(n_shp_in_basis, 4) = 1
!
         endif
!
      enddo
!
      call mem%alloc(basis_shell_info, n_shp_in_basis, 4)
      basis_shell_info(:, :) = basis_shell_info_full(1:n_shp_in_basis, :)
      call mem%dealloc(basis_shell_info_full, this%n_shp, 4)
!
!     Write basis_shell_data file containing
!
!        1. number shell pairs in basis
!        2. basis_shell_info
!        3. cholesky_basis
!
      call basis_info%open_('write','rewind')
!
      call basis_info%write_(n_shp_in_basis)
!
      call basis_info%write_(basis_shell_info, 4*n_shp_in_basis)
      call basis_info%write_(cholesky_basis_new, 3*this%n_cholesky)
!
      this%n_shp_in_basis = n_shp_in_basis
!
      call basis_info%close_()
!
      call mem%dealloc(basis_shell_info, n_shp_in_basis, 4)
      call mem%dealloc(cholesky_basis_new, this%n_cholesky, 3)
!
      call mem%dealloc(sig_shp, this%n_shp)
!
   end subroutine determine_cholesky_basis_standard
!
!
   subroutine construct_S(this, ao)
!!
!!    Construct S
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Constructs the overlap matrix (J|K) of the auxiliary basis
!!    and Cholesky decomposes it.
!!
      implicit none
!
      class(eri_cd) :: this
!
      type(ao_tool) :: ao
!
!     Local variables
!
!     Integers
!
      integer :: n_shp_in_basis, shp_in_basis
      integer :: n_vectors
      integer :: current_aop_in_shp
      integer :: A, B, C, D, AB, AB_shp, CD_shp
      integer :: I, J, K, L
      integer :: w, x, y, z, wx, yz
!
!     Integer allocatable arrays
!
!     Info on shells containing elements of the basis
      integer, dimension(:,:), allocatable :: basis_shell_info
      integer, dimension(:,:), allocatable :: basis_shell_info_full
!
!     ao and ao pair indices of the elements of the cholesky basis
      integer, dimension(:,:), allocatable :: cholesky_basis
      integer, dimension(:,:), allocatable :: cholesky_basis_updated
!
!     basis ao pairs in shell pair CD
      integer, dimension(3*ao%max_sh_size**2), target :: basis_aops_in_CD_shp
!
!     basis ao pairs in shell pair AB
      integer, dimension(3*ao%max_sh_size**2), target :: basis_aops_in_AB_shp
!
!     basis ao pairs in shell pair CD
      integer, dimension(:,:), pointer :: basis_aops_in_CD_shp_p
!
!     basis ao pairs in shell pair AB
      integer, dimension(:,:), pointer :: basis_aops_in_AB_shp_p
!
      integer, dimension(:,:), allocatable :: aops_in_basis
!
      integer, dimension(:), allocatable :: keep_vectors
!
!     Intervals
!
      type(range_) :: A_range, B_range, C_range, D_range
!
!     Reals
!
      real(dp) :: s_decomp_time, e_decomp_time, s_build_basis_time, e_build_basis_time
!
!     Real allocatable arrays
!
      real(dp), dimension(ao%max_sh_size**4), target :: g_AB_CD
      real(dp), dimension(:,:), pointer :: g_AB_CD_p
!
      real(dp), dimension(:,:), allocatable :: integrals_auxiliary
!
      real(dp), dimension(:), allocatable :: work  ! work array for LAPACK
!
      integer :: info, max_n_basis_aops_in_shp
!
!     Logicals
!
      logical :: found
!
      call cpu_time(s_build_basis_time)
!
!     Read basis_shell_data
!
      call this%cholesky_basis_file%open_('read','rewind')
!
      call this%cholesky_basis_file%read_(n_shp_in_basis)
!
      call mem%alloc(basis_shell_info, n_shp_in_basis, 4)
      call mem%alloc(cholesky_basis, this%n_cholesky, 3)
!
      call this%cholesky_basis_file%read_(basis_shell_info, 4*n_shp_in_basis)
      call this%cholesky_basis_file%read_(cholesky_basis, 3*this%n_cholesky)
!
      call this%cholesky_basis_file%close_('delete')
!
!     Construct integrals (J | J')
!
      call mem%alloc(integrals_auxiliary, this%n_cholesky, this%n_cholesky)
!
!$omp parallel do &
!$omp private(AB_shp, CD_shp, A, B, A_range, B_range, C, D, C_range, D_range, &
!$omp w, x, y, z, wx, yz, g_AB_CD, g_AB_CD_p, I, J, K, L, &
!$omp current_aop_in_shp, basis_aops_in_CD_shp, basis_aops_in_AB_shp, basis_aops_in_CD_shp_p, basis_aops_in_AB_shp_p) &
!$omp shared(integrals_auxiliary, cholesky_basis, basis_shell_info) &
!$omp schedule(guided)
      do CD_shp = 1, n_shp_in_basis
!
         C = basis_shell_info(CD_shp, 1)
         D = basis_shell_info(CD_shp, 2)
!
         C_range = ao%shells(C)
         D_range = ao%shells(D)
!
         basis_aops_in_CD_shp_p(1 : basis_shell_info(CD_shp, 4), 1 : 3) => &
               basis_aops_in_CD_shp(1 : basis_shell_info(CD_shp, 4)*3)
!
!        Determine which elements in the shell pair CD are elements of the basis
!
         current_aop_in_shp = 0
!
         do I = 1, this%n_cholesky
            if (cholesky_basis(I,3) == basis_shell_info(CD_shp, 3)) then
!
               current_aop_in_shp = current_aop_in_shp + 1
!
               basis_aops_in_CD_shp_p(current_aop_in_shp, 1) = cholesky_basis(I,1) - C_range%first + 1
               basis_aops_in_CD_shp_p(current_aop_in_shp, 2) = cholesky_basis(I,2) - D_range%first + 1
               basis_aops_in_CD_shp_p(current_aop_in_shp, 3) = I
!
            endif
         enddo
!
         do AB_shp = 1, n_shp_in_basis
!
            A = basis_shell_info(AB_shp, 1)
            B = basis_shell_info(AB_shp, 2)
!
            A_range = ao%shells(A)
            B_range = ao%shells(B)
!
            basis_aops_in_AB_shp_p(1 : basis_shell_info(AB_shp, 4), 1 : 3) => &
                  basis_aops_in_AB_shp(1 : basis_shell_info(AB_shp, 4)*3)
!
!           Determine which elements in the shell pair AB are elements of the basis
!
            current_aop_in_shp = 0
!
            do I = 1, this%n_cholesky
               if (cholesky_basis(I,3) == basis_shell_info(AB_shp, 3)) then
!
                  current_aop_in_shp = current_aop_in_shp + 1
!
                  basis_aops_in_AB_shp_p(current_aop_in_shp, 1) = cholesky_basis(I,1) - A_range%first + 1
                  basis_aops_in_AB_shp_p(current_aop_in_shp, 2) = cholesky_basis(I,2) - B_range%first + 1
                  basis_aops_in_AB_shp_p(current_aop_in_shp, 3) = I
!
               endif
            enddo
!
!           Construct integrals
!
            call ao%get_eri(g_AB_CD, A, B, C, D)
!
            g_AB_CD_p(1 : A_range%length*B_range%length, &
                        1 : C_range%length*D_range%length) &
                        => g_AB_CD(1 : A_range%length*B_range%length*C_range%length*D_range%length)
!
!           Only keep those that correspond to elements of the basis
!
            do I = 1, basis_shell_info(AB_shp, 4)
               do J = 1, basis_shell_info(CD_shp, 4)
!
                  y = basis_aops_in_CD_shp_p(J, 1)
                  z = basis_aops_in_CD_shp_p(J, 2)
                  yz = C_range%length*(z-1)+y
!
                  w = basis_aops_in_AB_shp_p(I, 1)
                  x = basis_aops_in_AB_shp_p(I, 2)
                  wx = A_range%length*(x-1)+w
!
                  K = basis_aops_in_AB_shp_p(I, 3)
                  L = basis_aops_in_CD_shp_p(J, 3)
!
                  integrals_auxiliary(K, L) = g_AB_CD_p(wx, yz)
                  integrals_auxiliary(L, K) = g_AB_CD_p(wx, yz)
!
               enddo
            enddo
!
         enddo ! AB
!
      enddo ! CD
!$omp end parallel do
!
      call mem%dealloc(basis_shell_info, n_shp_in_basis, 4)
!
      n_vectors = 0
      call mem%alloc(keep_vectors, (this%n_cholesky))
!
      call cpu_time(s_decomp_time)
!
      call mem%alloc(work, (2*this%n_cholesky))
!
!     DPSTRF computes the Cholesky factorization with complete pivoting
!     of a real symmetric positive semidefinite matrix.
!
      call dpstrf('L',                 &
            this%n_cholesky,         &
            integrals_auxiliary,       &
            this%n_cholesky,         &
            keep_vectors,              &
            n_vectors,                 &
            this%threshold*1.0d-1,   &
            work,                      &
            info)
!
      call mem%dealloc(work, (2*this%n_cholesky))
!
      call cpu_time(e_decomp_time)
!
      call mem%alloc(cholesky_basis_updated, n_vectors, 3)
!
!$omp parallel do private(I)
      do I = 1, n_vectors
!
         cholesky_basis_updated(I, :) = cholesky_basis(keep_vectors(I), :)
!
      enddo
!$omp end parallel do
!
      call mem%dealloc(cholesky_basis, this%n_cholesky, 3)
      call mem%dealloc(keep_vectors, (this%n_cholesky))
!
!     Write Q file containing
!
!        1. Cholesky factor Q, S = QQ^T
!
      call this%Q%open_('write', 'rewind')
!
      do i = 1,n_vectors
         call this%Q%write_(integrals_auxiliary(1:n_vectors, i), n_vectors)
      enddo
!
      call this%Q%close_()
!
      call mem%dealloc(integrals_auxiliary, this%n_cholesky, this%n_cholesky)
!
      this%n_cholesky = n_vectors
!
      call output%printf('n', '- Summary of Cholesky decomposition of &
                         &electronic repulsion integrals: ', fs='(/t3,a)')
      call output%printf('n', 'Final number of Cholesky vectors: (i0)', &
                         ints=[this%n_cholesky], fs='(/t6,a)')
!
!     Update the basis_shell_info array which contains information of which shell pairs (and shells)
!     contain elements of the basis and how many elements of the basis they contain.
!
      call mem%alloc(basis_shell_info_full, this%n_shp, 4) ! A, B, AB, n_basis_aops_in_shp
      basis_shell_info_full = 0
!
      n_shp_in_basis = 0
!
      do i = 1, this%n_cholesky
!
         A = ao%ao_to_shell(cholesky_basis_updated(i, 1))
         B = ao%ao_to_shell(cholesky_basis_updated(i, 2))
!
         AB = get_shp_from_shells(A, B, this%n_s)
!
         found = .false.
!
         do shp_in_basis = 1, n_shp_in_basis
!
              if (AB == basis_shell_info_full(shp_in_basis, 3)) then ! This shell pair is already
                                                                    ! in the list, so must only
                                                                    ! increment n_aops_in_shp.
!
                 found = .true.
                 basis_shell_info_full(shp_in_basis, 4) = basis_shell_info_full(shp_in_basis, 4) + 1
!
                 exit ! Loop over shp_in_basis
!
              endif
!
         enddo
!
         if(.not. found) then ! First aop in this shell pair was found,
                              ! set [A, B, AB, n_aops_in_shp = 1]
!
            n_shp_in_basis = n_shp_in_basis + 1
!
            basis_shell_info_full(n_shp_in_basis, 1) = A
            basis_shell_info_full(n_shp_in_basis, 2) = B
            basis_shell_info_full(n_shp_in_basis, 3) = AB
            basis_shell_info_full(n_shp_in_basis, 4) = 1
!
         endif
!
      enddo
!
!     Copy into correctly shaped array
!
      call mem%alloc(basis_shell_info, n_shp_in_basis, 4)
!
!$omp parallel do private (I, J) collapse(2)
      do J = 1, 4
         do I = 1, n_shp_in_basis

            basis_shell_info(I, J) = basis_shell_info_full(I, J)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(basis_shell_info_full, this%n_shp, 4)
!
!     Prepare array for RI-like expressions, either to construct Cholesky vectors
!     or for RI-like implementations using the Cholesky basis
!
      max_n_basis_aops_in_shp = maxval(basis_shell_info(1:n_shp_in_basis,4))
!
      call mem%alloc(aops_in_basis, n_shp_in_basis, &
                              3*max_n_basis_aops_in_shp)  ! [x, y] where x and y
                                                          ! are the indices of the
                                                          ! AOs in the shell they
                                                          ! belong to.
!
     aops_in_basis = 0
!
    do AB_shp = 1, n_shp_in_basis
!
      A  = basis_shell_info(AB_shp, 1)
      B  = basis_shell_info(AB_shp, 2)
      AB = basis_shell_info(AB_shp, 3)
!
      A_range = ao%shells(A)
      B_range = ao%shells(B)
!
!     Determine which elements in the shell pair AB are elements of the basis
!
      current_aop_in_shp = 0
!
      do I = 1, this%n_cholesky
         if (cholesky_basis_updated(I,3) == AB) then
!
            current_aop_in_shp = current_aop_in_shp + 1
!
            aops_in_basis(AB_shp, current_aop_in_shp) &
               = cholesky_basis_updated(I,1) - A_range%first + 1
!
            aops_in_basis(AB_shp, max_n_basis_aops_in_shp + current_aop_in_shp) &
               = cholesky_basis_updated(I,2) - B_range%first + 1
!
            aops_in_basis(AB_shp, 2*max_n_basis_aops_in_shp + current_aop_in_shp) &
               = I
!
         endif
      enddo
!
      if (current_aop_in_shp .ne. basis_shell_info(AB_shp, 4)) &
         call output%error_msg('something went wrong in construct_S.')
!
    enddo
!
!     Write basis_shell_data file containing
!
!        1. number shell pairs in basis
!        2. basis_shell_info
!        3. cholesky_basis
!
      call this%cholesky_basis_file%open_('write')
!
      call this%cholesky_basis_file%write_(n_shp_in_basis)
      call this%cholesky_basis_file%write_(basis_shell_info, 4*n_shp_in_basis)
      call this%cholesky_basis_file%write_(cholesky_basis_updated, 3*n_vectors)
      call this%cholesky_basis_file%write_(max_n_basis_aops_in_shp)
      call this%cholesky_basis_file%write_(aops_in_basis, n_shp_in_basis*3*max_n_basis_aops_in_shp)
!
      call this%cholesky_basis_file%close_()
!
      call mem%dealloc(basis_shell_info, n_shp_in_basis, 4)
      call mem%dealloc(cholesky_basis_updated, n_vectors, 3)
      call mem%dealloc(aops_in_basis, n_shp_in_basis, 3*max_n_basis_aops_in_shp)
!
      call cpu_time(e_build_basis_time)
!
   end subroutine construct_S
!
!
   subroutine invert_Q(this)
!!
!!    Invert Q
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, July 2018
!!
!!    Invert Cholesky factors Q of S
!!
      implicit none
!
      class(eri_cd) :: this
!
      real(dp):: s_invert_time, e_invert_time
!
      real(dp), dimension(:,:), allocatable :: cholesky_inverse
!
      integer :: I, info
!
      call cpu_time(s_invert_time)
!
!     Read Cholesky vectors of auxiliary basis overlap
!
      call this%Q%open_('read', 'rewind')
!
      call mem%alloc(cholesky_inverse, this%n_cholesky, this%n_cholesky)
!
      do i = 1,this%n_cholesky
         call this%Q%read_(cholesky_inverse(:,i), this%n_cholesky)
      enddo
!
      call this%Q%close_()
!
!     Invert cholesky vectors
!
      call dtrtri('l','n', this%n_cholesky, cholesky_inverse, this%n_cholesky, info)
!
      if (info /= 0) then
         call output%error_msg('Matrix inversion failed.' // &
                               ' "Dtrtri" finished with info: (i0)', ints=[info])
      end if
!
!     Write inverse Cholesky vectors of auxiliary basis overlap
!
!        1. n_cholesky: Number of elements of the basis
!        2. cholesky_inverse (n_cholesky, n_cholesky)
!
      call this%Q_inverse%open_('write', 'rewind')
!
!     Write out columns, but only the parts that are not rubbish left over from dpstrf
!
      do i = 1, this%n_cholesky
!
         call this%Q_inverse%write_(cholesky_inverse(1 + (i - 1) : this%n_cholesky, i), &
                                      this%n_cholesky + 1 - i)
!
      enddo
!
      call this%Q_inverse%close_()
!
      call mem%dealloc(cholesky_inverse, this%n_cholesky, this%n_cholesky)
!
!     Timings
!
      call cpu_time(e_invert_time)
!
!
   end subroutine invert_Q
!
!
   pure function get_size_shp(A_range, B_range) result(size_shp)
!!
!!    Get size shell pair
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Returns size of diagonal for a given shell pair.
!!
      implicit none
!
      type(range_), intent(in) :: A_range
      type(range_), intent(in) :: B_range
!
      integer :: size_shp
!
      if (A_range%first == B_range%first) then
!
         size_shp = A_range%length*(A_range%length + 1)/2
!
      else
!
         size_shp = (A_range%length)*(B_range%length)
!
      endif
!
   end function get_size_shp
!
!
   pure function get_shp_from_shells(s1, s2, n_s) result(shp)
!!
!!    Get shell pair from shells,
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      integer, intent(in) :: s1, s2, n_s
      integer :: shp
      integer :: A, B
!
      shp = 1
!
      if (s1 .ge. s2) then
!
         do B = 1, n_s
            do A = B, n_s
!
               if (s1 == A .and. s2 == B)  return
!
               shp = shp + 1
!
            enddo
         enddo
!
      else
!
         do B = 1, n_s
            do A = B, n_s
!
               if (s2 == A .and. s1 == B)  return
!
               shp = shp + 1
!
            enddo
         enddo
!
      endif
!
   end function get_shp_from_shells
!
!
   subroutine read_settings(this)
!!
!!    Read settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!     Read input if it is present:
!!
!!        this cholesky
!!           threshold: 1.0d-5
!!           span: 1.0d-2
!!           qualified: 1000
!!           one center
!!        end this cholesky
!!
      use global_in, only: input
!
      implicit none
!
      class(eri_cd) :: this
!
      call input%get_keyword('threshold', 'solver cholesky', this%threshold)
      call input%get_keyword('span', 'solver cholesky', this%span)
      call input%get_keyword('batches', 'solver cholesky', this%n_batches)
      call input%get_keyword('qualified', 'solver cholesky', this%max_qual)
!
      this%one_center = input%is_keyword_present('one center', 'solver cholesky')
!
   end subroutine read_settings
!
!
   subroutine print_banner(this)
!!
!!    Print banner
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(eri_cd) :: this
!
      call output%printf('m', ' - ' // trim(this%name_), fs='(/t3,a)')
      call output%print_separator('m', len(trim(this%name_)) + 6, '-')
!
      call output%printf('m', this%description1, ffs='(/t3,a)')
      call output%printf('m', this%description2, ffs='(/t3,a)')
      call output%printf('m', this%description3, ffs='(/t3,a)')
!
   end subroutine print_banner
!
!
   subroutine print_settings(this)
!!
!!    Print settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(eri_cd) :: this
!
      call output%printf('m', '- Cholesky decomposition settings:', fs='(/t3, a)')
!
      call output%printf('m', 'Target threshold is: (e10.2)', &
                         reals=[this%threshold], fs='(/t6, a)')
      call output%printf('m', 'Span factor:         (e10.2)', &
                         reals=[this%span], fs='(t6, a)')
      call output%printf('m', 'Max qual:            (i10)', &
                         ints=[this%max_qual], fs='(t6, a)')
!
      if (this%one_center) then
!
         call output%printf('m', 'Doing one-center Cholesky decomposition (1C-CD).', &
                            fs='(/t6, a)')
!
      endif
!
   end subroutine print_settings
!
!
   subroutine construct_cholesky_mo_vectors_eri_cd(this, ao, n_ao, n_mo, &
                                                   orbital_coefficients, cd_tool)
!!
!!    Construct Cholesky MO vectors
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Feb 2020
!!
!!    Constructs the Cholesky MO vectors
!!
!!       L_Jpq = C_wpCxq (wx|K)[Q^{-T}]_KJ,
!!
!!    where C_wp are the orbital coefficients,
!!    and Q is the Cholesky factor of the overlap of the elements
!!    of the Cholesky basis, S.
!!
!!    Passes the Cholesky vectors to the integrals (eri_tool)
!!
      use array_utilities, only: zero_array, transpose_
      use abstract_eri_cholesky_class, only: abstract_eri_cholesky
      use batching_index_class, only: batching_index
      use reordering, only: sort_123_to_213, sort_123_to_312
!
      implicit none
!
      class(eri_cd) :: this
!
      type(ao_tool) :: ao
!
      integer, intent(in) :: n_mo, n_ao
!
      real(dp), dimension(n_ao, n_mo), intent(in) :: orbital_coefficients
!
      class(abstract_eri_cholesky), intent(inout) :: cd_tool
!
      integer :: A, B, C, D, K_shp, AB_shp
      integer :: w, x, y, z
      integer :: q, p
      integer :: I, J
      integer :: shp
      integer :: n_construct_shp, n_construct_aop
      integer :: n_shp_in_basis
!
      integer, dimension(:,:), allocatable :: basis_shell_info ! Info on shells containing elements of the basis
      integer, dimension(:,:), allocatable :: cholesky_basis   ! Info on cholesky basis
      integer, dimension(:,:), allocatable :: aops_in_basis
      integer, dimension(:,:), allocatable   :: AB_info
!
      real(dp), dimension(ao%max_sh_size**4), target :: g_ABCD
      real(dp), dimension(:,:,:,:), pointer :: g_ABCD_p
!
      real(dp), dimension(:,:,:), allocatable :: g_wxK, g_pxK, g_xpK, g_qpK, g_qpK_red
      real(dp), dimension(:,:,:), allocatable :: L_qpJ, L_Jqp
      real(dp), dimension(:,:), allocatable :: aux_chol_inverse, aux_chol_inverse_transpose
!
!     Logical allocatable arrays
!
      logical, dimension(:), allocatable :: construct_shp
!
!     Intervals
!
      type(range_) :: A_range, B_range, C_range, D_range
!
      integer :: max_n_basis_aops_in_shp
!
!     Batching and memory handling
!
      integer :: req0, req1
      integer :: current_p_batch
!
      type(batching_index) :: batch_p
!
      type(timings) :: timer
!
      timer = timings('Construct and write Cholesky vectors', pl='n')
      call timer%turn_on()
!
!     Read diagonal info
!
      call mem%alloc(construct_shp, this%n_shp)
!
      call this%diagonal_info_cauchy_schwarz%open_('read', 'rewind')
!
      call this%diagonal_info_cauchy_schwarz%read_(n_construct_shp)
      call this%diagonal_info_cauchy_schwarz%read_(n_construct_aop)
      call this%diagonal_info_cauchy_schwarz%read_(construct_shp, this%n_shp)
!
      call this%diagonal_info_cauchy_schwarz%close_()
!
!     Read inverse of cholesky vectors of auxiliary overlap, Q^-1
!
      call this%Q_inverse%open_('read', 'rewind')
!
      call mem%alloc(aux_chol_inverse, this%n_cholesky, this%n_cholesky)
      call zero_array(aux_chol_inverse, this%n_cholesky**2)
!
      do I = 1, this%n_cholesky
!
         call this%Q_inverse%read_(aux_chol_inverse(1 + (i - 1) : this%n_cholesky, i), &
                                     this%n_cholesky + 1 - i)
!
      enddo
!
      call this%Q_inverse%close_()
!
!     Transpose Q^-1
!
      call mem%alloc(aux_chol_inverse_transpose, this%n_cholesky, this%n_cholesky)
      call transpose_(aux_chol_inverse, aux_chol_inverse_transpose, this%n_cholesky)
      call mem%dealloc(aux_chol_inverse, this%n_cholesky, this%n_cholesky)
!
!     Read Cholesky basis information
!
      call this%cholesky_basis_file%open_('read','rewind')
      call this%cholesky_basis_file%read_(n_shp_in_basis)
!
      call mem%alloc(basis_shell_info, n_shp_in_basis, 4)
      call mem%alloc(cholesky_basis, this%n_cholesky, 3)
!
      call this%cholesky_basis_file%read_(basis_shell_info, 4*n_shp_in_basis)
      call this%cholesky_basis_file%read_(cholesky_basis, 3*this%n_cholesky)
      call this%cholesky_basis_file%read_(max_n_basis_aops_in_shp)
!
      call mem%alloc(aops_in_basis, n_shp_in_basis, 3*max_n_basis_aops_in_shp)
!
      call this%cholesky_basis_file%read_(aops_in_basis, n_shp_in_basis*3*max_n_basis_aops_in_shp)
!
      call this%cholesky_basis_file%close_()
!
!     Prepare for OMP
!
      call mem%alloc(AB_info, n_construct_shp, 2) ! [A, B]
!
      shp = 0
      AB_shp = 0
!
      do B = 1, this%n_s
         do A = B, this%n_s
!
            AB_shp = AB_shp + 1
!
            if (construct_shp(AB_shp)) then
!
               shp = shp + 1
!
               AB_info(shp, 1) = A
               AB_info(shp, 2) = B
!
            endif
!
         enddo
      enddo
!
!
!     Prepare for batching
!
      req0 = n_ao**2*max_n_basis_aops_in_shp
      req1 = max(2*n_ao*max_n_basis_aops_in_shp     &
                  + n_mo*this%n_cholesky,        &
                 2*n_mo*this%n_cholesky)
!
      batch_p = batching_index(n_mo)
!
      call mem%batch_setup(batch_p, req0, req1, 'construct_cholesky_mo_vectors')
!
!     Loop over the number of a batches
!
      do current_p_batch = 1, batch_p%num_batches
!
         call batch_p%determine_limits(current_p_batch)
!
         call mem%alloc(g_qpK, n_mo, batch_p%length, this%n_cholesky)
!
         do K_shp = 1, n_shp_in_basis
!
            call mem%alloc(g_wxK, n_ao, n_ao, basis_shell_info(K_shp,4))
            call zero_array(g_wxK, n_ao**2*basis_shell_info(K_shp,4))
!
            C = basis_shell_info(K_shp, 1)
            D = basis_shell_info(K_shp, 2)
!
            C_range = ao%shells(C)
            D_range = ao%shells(D)
!
!$omp parallel do private(AB_shp, A, B, A_range, B_range, g_ABCD, g_ABCD_p, w, x, J, y, z)
            do AB_shp = 1, n_construct_shp
!
               A = AB_info(AB_shp, 1)
               B = AB_info(AB_shp, 2)
!
               B_range = ao%shells(B)
               A_range = ao%shells(A)
!
               call ao%get_eri(g_ABCD, A, B, C, D)
!
               g_ABCD_p(1 : A_range%length, 1 : B_range%length, &
                        1 : C_range%length, 1 : D_range%length) &
                        => g_ABCD(1 : A_range%length &
                                     *B_range%length &
                                     *C_range%length &
                                     *D_range%length)
!
               if (A == B) then
!
                  do w = A_range%first, A_range%get_last()
                     do x = B_range%first, B_range%get_last()
                        do J = 1, basis_shell_info(K_shp, 4)
!
                           y = aops_in_basis(K_shp, J)
                           z = aops_in_basis(K_shp, J + max_n_basis_aops_in_shp)

                           g_wxK(w, x, J) = g_ABCD_p(w - A_range%first + 1, &
                                                     x - B_range%first + 1, y, z)
!
                        enddo
                     enddo
                  enddo
!
               else
!
                  do w = A_range%first, A_range%get_last()
                     do x = B_range%first, B_range%get_last()
                        do J = 1, basis_shell_info(K_shp, 4)
!
                           y = aops_in_basis(K_shp, J)
                           z = aops_in_basis(K_shp, J + max_n_basis_aops_in_shp)

                           g_wxK(w, x, J) = g_ABCD_p(w - A_range%first + 1, &
                                                     x - B_range%first + 1, y, z)
                           g_wxK(x, w, J) = g_ABCD_p(w - A_range%first + 1, &
                                                     x - B_range%first + 1, y, z)
!
                        enddo
                     enddo
                  enddo
!
               endif

            enddo ! end loop over AB_shp
!$omp end parallel do
!
!           Transform wx -> pq
!
            call mem%alloc(g_pxK, batch_p%length, n_ao, basis_shell_info(K_shp,4))
!
            call dgemm('T', 'N',                               &
                        batch_p%length,                        &
                        n_ao*basis_shell_info(K_shp,4),        &
                        n_ao,                                  &
                        one,                                   &
                        orbital_coefficients(1,batch_p%first), & ! C_w_p
                        n_ao,                                  &
                        g_wxK,                                 & ! g_w_xK
                        n_ao,                                  &
                        zero,                                  &
                        g_pxK,                                 &
                        batch_p%length)
!
            call mem%dealloc(g_wxK, n_ao, n_ao, basis_shell_info(K_shp,4))
!
            call mem%alloc(g_xpK, n_ao, batch_p%length, basis_shell_info(K_shp,4))
!
            call sort_123_to_213(g_pxK, g_xpK, batch_p%length, n_ao, basis_shell_info(K_shp,4))
!
            call mem%dealloc(g_pxK, batch_p%length, n_ao, basis_shell_info(K_shp,4))
!
            call mem%alloc(g_qpK_red, n_mo, batch_p%length, basis_shell_info(K_shp,4))
!
            call dgemm('T', 'N',                                  &
                        n_mo,                                     &
                        batch_p%length*basis_shell_info(K_shp,4), &
                        n_ao,                                     &
                        one,                                      &
                        orbital_coefficients,                     & ! C_x_q
                        n_ao,                                     &
                        g_xpK,                                    &
                        n_ao,                                     &
                        zero,                                     &
                        g_qpK_red,                                &
                        n_mo)
!
            call mem%dealloc(g_xpK, n_ao, batch_p%length, basis_shell_info(K_shp,4))
!
!$omp parallel do private(p, q, J)
            do J = 1, basis_shell_info(K_shp,4)
               do p = 1, batch_p%length
                  do q = 1, n_mo

!
                     g_qpK(q,p, aops_in_basis(K_shp, J + max_n_basis_aops_in_shp*2)) = g_qpK_red(q,p,J)
!
                  enddo
               enddo
            enddo
!$omp end parallel do
!
            call mem%dealloc(g_qpK_red, n_mo, batch_p%length, basis_shell_info(K_shp,4))
!
         enddo ! end loop over K_shp
!
         call mem%alloc(L_qpJ, n_mo, batch_p%length, this%n_cholesky)
!
         call dgemm('N', 'N',                      &
                     n_mo*batch_p%length,          &
                     this%n_cholesky,            &
                     this%n_cholesky,            &
                     one,                          &
                     g_qpK,                        &
                     n_mo*batch_p%length,          &
                     aux_chol_inverse_transpose,   &
                     this%n_cholesky,            &
                     zero,                         &
                     L_qpJ,                        &
                     n_mo*batch_p%length)
!
         call mem%dealloc(g_qpK, n_mo, batch_p%length, this%n_cholesky)
!
         call mem%alloc(L_Jqp, this%n_cholesky, n_mo, batch_p%length)
         call sort_123_to_312(L_qpJ, L_Jqp, n_mo, batch_p%length, this%n_cholesky)
         call mem%dealloc(L_qpJ, n_mo, batch_p%length, this%n_cholesky)
!
         call cd_tool%set(L_Jqp,           &
                          1,               &
                          n_mo,            &
                          batch_p%first,   &
                          batch_p%get_last())
!
         call mem%dealloc(L_Jqp, this%n_cholesky, n_mo, batch_p%length)
!
      enddo
!
      call mem%batch_finalize()
!
      call mem%dealloc(AB_info, n_construct_shp, 2) ! [A, B]
!
      call mem%dealloc(cholesky_basis, this%n_cholesky, 3)
      call mem%dealloc(basis_shell_info, n_shp_in_basis, 4)
      call mem%dealloc(aops_in_basis, n_shp_in_basis, 3*max_n_basis_aops_in_shp)
!
      call mem%dealloc(aux_chol_inverse_transpose, this%n_cholesky, this%n_cholesky)
      call mem%dealloc(construct_shp, this%n_shp)
!
      call timer%turn_off()
!
   end subroutine construct_cholesky_mo_vectors_eri_cd
!
!
   subroutine diagonal_test_eri_cd(this, ao)
!!
!!    Cholesky vectors diagonal test
!!    Written by Sarai D. Folkestad, Feb 2020
!!
!!    Tests the decomposition by
!!
!!       1. finding the largest element of (D_sig - D_approx)
!!       2. finding the smallest (largest negative) element of (D_sig - D_approx)
!!
!
      use array_utilities, only: zero_array, transpose_, get_abs_max
      use batching_index_class, only: batching_index
!
      implicit none
!
      class(eri_cd) :: this
!
      type(ao_tool), intent(in) :: ao
!
      real(dp), dimension(ao%max_sh_size**4), target :: g_ABAB, g_ABCD
      real(dp), dimension(:,:,:,:), pointer :: g_ABAB_p, g_ABCD_p
!
      real(dp), dimension(:), allocatable :: D_xy
      real(dp), dimension(:,:), allocatable :: aux_chol_inverse, aux_chol_inverse_transpose
      real(dp), dimension(:,:), allocatable :: g_wxK, L_wxJ
!
      real(dp) :: max_diff, min_diff
!
      integer :: aop, n_construct_shp, n_construct_aop, current_construct_shp
      integer :: A, B, C, D
      integer :: AB_shp
      integer :: K, I, J
      integer :: x, y, xy, xy_packed, shp
      integer :: K_shp, w, wx, z
      integer :: n_shp_in_basis
!
      integer, dimension(:), allocatable :: ao_offsets
      integer, dimension(:,:), allocatable :: construct_shp_to_shells
!
      logical, dimension(:), allocatable :: construct_shp
!
      integer, dimension(:,:), allocatable :: basis_shell_info ! Info on shells containing elements of the basis
      integer, dimension(:,:), allocatable :: cholesky_basis   ! Info on cholesky basis
      integer, dimension(:,:), allocatable :: aops_in_basis
!
      type(range_) :: A_range, B_range, C_range, D_range
!
      integer :: max_n_basis_aops_in_shp
!
!     Batching and memory handling
!
      integer :: req0, req1
      integer :: current_J_batch
!
      type(batching_index) :: batch_J
!
      type(timings), allocatable :: timer
!
      timer = timings('Diagonal test', pl='n')
      call timer%turn_on()
!
!     Read diagonal info
!
      call mem%alloc(construct_shp, this%n_shp)
!
      call this%diagonal_info_cauchy_schwarz%open_('read', 'rewind')
!
      call this%diagonal_info_cauchy_schwarz%read_(n_construct_shp)
      call this%diagonal_info_cauchy_schwarz%read_(n_construct_aop)
      call this%diagonal_info_cauchy_schwarz%read_(construct_shp, this%n_shp)
!
      call this%diagonal_info_cauchy_schwarz%close_()
!
!     Read inverse of cholesky vectors of auxiliary overlap, Q^-1
!
      call this%Q_inverse%open_('read', 'rewind')
!
      call mem%alloc(aux_chol_inverse, this%n_cholesky, this%n_cholesky)
      call zero_array(aux_chol_inverse, this%n_cholesky**2)
!
      do I = 1, this%n_cholesky
!
         call this%Q_inverse%read_(aux_chol_inverse(1 + (i - 1) : this%n_cholesky, i), &
                                     this%n_cholesky + 1 - i)
!
      enddo
!
      call this%Q_inverse%close_()
!
!     Transpose Q^-1
!
      call mem%alloc(aux_chol_inverse_transpose, this%n_cholesky, this%n_cholesky)
      call transpose_(aux_chol_inverse, aux_chol_inverse_transpose, this%n_cholesky)
      call mem%dealloc(aux_chol_inverse, this%n_cholesky, this%n_cholesky)
!
!     Read Cholesky basis information
!
      call this%cholesky_basis_file%open_('read','rewind')
      call this%cholesky_basis_file%read_(n_shp_in_basis)
!
      call mem%alloc(basis_shell_info, n_shp_in_basis, 4)
      call mem%alloc(cholesky_basis, this%n_cholesky, 3)
!
      call this%cholesky_basis_file%read_(basis_shell_info, 4*n_shp_in_basis)
      call this%cholesky_basis_file%read_(cholesky_basis, 3*this%n_cholesky)
      call this%cholesky_basis_file%read_(max_n_basis_aops_in_shp)
!
      call mem%alloc(aops_in_basis, n_shp_in_basis, 3*max_n_basis_aops_in_shp)
!
      call this%cholesky_basis_file%read_(aops_in_basis, n_shp_in_basis*3*max_n_basis_aops_in_shp)
!
      call this%cholesky_basis_file%close_()
!
!     Prepare for construction of diagonal
!
      shp = 0        ! Shell pair number
!
      call mem%alloc(ao_offsets, n_construct_shp)
      ao_offsets = 0
!
      current_construct_shp = 0
!
      call mem%alloc(construct_shp_to_shells, n_construct_shp, 2)
!
      do B = 1, this%n_s
         do A = B, this%n_s
!
            shp = shp + 1
!
            if (construct_shp(shp)) then
!
               current_construct_shp = current_construct_shp + 1
!
               A_range = ao%shells(A)
               B_range = ao%shells(B)
!
               construct_shp_to_shells(current_construct_shp, 1) = A
               construct_shp_to_shells(current_construct_shp, 2) = B
!
               if (current_construct_shp .lt. n_construct_shp) then
!
                  ao_offsets(current_construct_shp + 1) = ao_offsets(current_construct_shp) + &
                           get_size_shp(A_range, B_range)
!
               endif
!
            endif
!
         enddo
      enddo
!
      call mem%dealloc(construct_shp, this%n_shp)
!
!     Construct significant diagonal
!
      call mem%alloc(D_xy, n_construct_aop)
!
!$omp parallel do &
!$omp private(I, A, B, A_range, B_range, x, y, xy, xy_packed, g_ABAB, g_ABAB_p) &
!$omp shared(D_xy, ao_offsets) &
!$omp schedule(guided)
      do AB_shp = 1, n_construct_shp
!
         A = construct_shp_to_shells(AB_shp, 1)
         B = construct_shp_to_shells(AB_shp, 2)
!
         A_range = ao%shells(A)
         B_range = ao%shells(B)
!
         call ao%get_eri(g_ABAB, A, B, A, B)
!
         g_ABAB_p(1 : A_range%length, &
                  1 : B_range%length, &
                  1 : A_range%length, &
                  1 : B_range%length) &
                        => g_ABAB(1 : (A_range%length)**2*(B_range%length)**2)
!
         if (A .eq. B) then
!
            do x = 1, A_range%length
               do y = x, B_range%length
!
                  xy_packed = (max(x,y)*(max(x,y)-3)/2) + x + y
!
                  D_xy(xy_packed + ao_offsets(AB_shp)) = g_ABAB_p(x, y, x, y)
!
               enddo
            enddo
!
         else ! A ≠ B
!
            do x = 1, (A_range%length)
               do y = 1, (B_range%length)
!
                  xy = A_range%length*(y - 1) + x
                  D_xy(xy + ao_offsets(AB_shp)) = g_ABAB_p(x, y, x, y)
!
               enddo
            enddo
!
         endif
!
      enddo
!$omp end parallel do
!
!     Prepare for batching
!
      req0 = n_construct_aop*max_n_basis_aops_in_shp
      req1 = n_construct_aop
!
      batch_J = batching_index(this%n_cholesky)
!
      call mem%batch_setup(batch_J, req0, req1, 'diagonal_test')
!
!     Loop over the number of a batches
!
      do current_J_batch = 1, batch_J%num_batches
!
         call batch_J%determine_limits(current_J_batch)
!
         call mem%alloc(L_wxJ, n_construct_aop, batch_J%length)
         call zero_array(L_wxJ, n_construct_aop*batch_J%length)
!
         do K_shp = 1, n_shp_in_basis
!
            C = basis_shell_info(K_shp, 1)
            D = basis_shell_info(K_shp, 2)
!
            C_range = ao%shells(C)
            D_range = ao%shells(D)
!
            call mem%alloc(g_wxK, n_construct_aop, basis_shell_info(K_shp,4))
!
!$omp parallel do private(AB_shp, A, B, A_range, B_range, g_ABCD, g_ABCD_p, w, x, K, wx, y, z)
            do AB_shp = 1, n_construct_shp
!
               A = construct_shp_to_shells(AB_shp, 1)
               B = construct_shp_to_shells(AB_shp, 2)
!
               A_range = ao%shells(A)
               B_range = ao%shells(B)
!
               call ao%get_eri(g_ABCD, A, B, C, D)
!
               g_ABCD_p(1 : A_range%length, 1 : B_range%length, &
                        1 : C_range%length, 1 : D_range%length) &
                        => g_ABCD(1 : A_range%length &
                                     *B_range%length &
                                     *C_range%length &
                                     *D_range%length)
!
               if (A == B) then
!
                  do w = A_range%first, A_range%get_last()
                     do x = B_range%first, B_range%get_last()
                        do K = 1, basis_shell_info(K_shp, 4)
!
                           y = aops_in_basis(K_shp, K)
                           z = aops_in_basis(K_shp, K + max_n_basis_aops_in_shp)
!
                           wx = (max(x - B_range%first + 1, w - A_range%first + 1)*&
                                (max(x - B_range%first + 1, w - A_range%first + 1)-3)/2) &
                              + x - B_range%first + 1 + w - A_range%first + 1
!
                           g_wxK(wx + ao_offsets(AB_shp), K) = &
                                          g_ABCD_p(w - A_range%first + 1, &
                                                   x - B_range%first + 1, y, z)

!
                        enddo
                     enddo
                  enddo
!
               else
!
                  do w = A_range%first, A_range%get_last()
                     do x = B_range%first, B_range%get_last()
                        do K = 1, basis_shell_info(K_shp, 4)
!
                           y = aops_in_basis(K_shp, K)
                           z = aops_in_basis(K_shp, K + max_n_basis_aops_in_shp)
!
                           wx = A_range%length*(x - B_range%first) + w - A_range%first + 1
!
                           g_wxK(wx + ao_offsets(AB_shp), K) = &
                                          g_ABCD_p(w - A_range%first + 1, &
                                                   x - B_range%first + 1, y, z)
!
                        enddo
                     enddo
                  enddo
!
               endif
            enddo ! end loop over AB_shp
!$omp end parallel do
!
!$omp parallel do private(J, wx, K)
            do J = 1, batch_J%length
               do wx = 1, n_construct_aop
                  do K = 1, basis_shell_info(K_shp, 4)

                     L_wxJ(wx, J) = L_wxJ(wx, J) + g_wxK(wx,K) &
                                       *aux_chol_inverse_transpose(aops_in_basis(K_shp, &
                                          K + max_n_basis_aops_in_shp*2), J + batch_J%first - 1)
                  enddo
               enddo
            enddo
!$omp end parallel do
!
            call mem%dealloc(g_wxK, n_construct_aop, basis_shell_info(K_shp,4))
!
         enddo ! K_shp
!
!$omp parallel do private(wx, J)
         do wx = 1, n_construct_aop
            do J = 1, batch_J%length
!
                  D_xy(wx) = D_xy(wx) - L_wxJ(wx, J)**2
!
            enddo
         enddo
!$omp end parallel do
!
         call mem%dealloc(L_wxJ, n_construct_aop, batch_J%length)
!
      enddo ! J batch
!
      call mem%batch_finalize()
!
      call mem%dealloc(ao_offsets, n_construct_shp)
!
      call mem%dealloc(construct_shp_to_shells, n_construct_shp, 2)
!
      call mem%dealloc(cholesky_basis, this%n_cholesky, 3)
      call mem%dealloc(basis_shell_info, n_shp_in_basis, 4)
      call mem%dealloc(aops_in_basis, n_shp_in_basis, 3*max_n_basis_aops_in_shp)
!
      call mem%dealloc(aux_chol_inverse_transpose, this%n_cholesky, this%n_cholesky)
!
!     Calculate maximal difference and minimal difference
!
      max_diff = get_abs_max(D_xy, n_construct_aop)
!
      min_diff = 1.0d5
!
      do aop = 1, n_construct_aop
         if (D_xy(aop) .lt. min_diff) min_diff = D_xy(aop)
      enddo
!
      call mem%dealloc(D_xy, n_construct_aop)
!
      call output%printf('m', '- Testing the Cholesky decomposition &
                         &decomposition electronic repulsion integrals:', &
                         ll=100, fs='(/t2,a)')
!
      call output%printf('m', 'Maximal difference between approximate and &
                         &actual diagonal: (e23.4)', reals=[max_diff], ll=100, fs='(/t6,a)')
      call output%printf('m', 'Minimal element of difference between &
                         &approximate and actual diagonal: (e12.4)', &
                         reals=[min_diff], ll=100, fs='( t6,a)')
!
      call timer%turn_off()
!
   end subroutine diagonal_test_eri_cd
!
!
end module eri_cd_class
