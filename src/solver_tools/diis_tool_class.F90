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
module diis_tool_class
!
!!
!!    DIIS tool class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2018
!!
!!    Objects of this class can be used to solve a non-linear set of equations e(x) = 0,
!!    where e and x are vectors of length n_equations and n_parameters.
!!
!!    Typical usage:
!!
!!       diis = diis_tool(...)
!!
!!       call this%initialize()
!!
!!       do while (.not. converged)
!!
!!          -> Calculate error vector e given the current parameters x
!!
!!          -> Set x to some estimated next value (x = x + dx), e.g. using Quasi-Newton
!!
!!          call this%update(e, x)
!!
!!          -> Now, x contains the extrapolated DIIS guess (*)
!!
!!       enddo
!!
!!       call this%finalize()
!!
!!    In DIIS, the squared norm of the averaged error vector
!!
!!       < e(x) > = sum_k w_k e(x)_k
!!
!!    is minimized with respect to the condition sum_k w_k = 1. The parameters
!!    are updated as
!!
!!       x = sum_k=1^n w_k x_k     (*)
!!
!!    where "x_k" is the current previous estimates for x (e.g., using quasi-Newton).
!!
!!    See Pulay, P. Convergence acceleration of iterative sequences.
!!    The case of SCF iteration. Chem. Phys. Lett. 1980, 73, 393−398.
!!
!
   use parameters
!
   use global_out,            only: output
   use array_initialization,  only: zero_array
   use sequential_file_class, only: sequential_file
   use record_storer_class,   only: record_storer
   use memory_manager_class,  only: mem
   use batching_index_class,  only: batching_index
   use timings_class,         only: timings
!
   type :: diis_tool
!
      character(len=40), private :: name_ ! determines the prefix of all DIIS files
!
      integer, private :: iteration       ! Increments by +1 when 'update' is called
!
      integer, private :: dimension_      ! Standard is 8, i.e. to keep at most 8 vectors
!
      class(record_storer), allocatable, private :: errors, parameters
!
      integer, private :: n_parameters    ! The length of the parameter (x) vectors
      integer, private :: n_equations     ! The length of the error (e) vectors
!
      logical, private :: crop            ! Conjugate residual with optimal trial vectors
                                          ! J. Chem. Theory Comput. 2015, 11, 4, 1518-1524
!
      integer, private :: dim_G
!
      real(dp), dimension(:,:), allocatable, private :: G     ! DIIS matrix, G_ij = e_i^T e_j
!
      real(dp), dimension(:), allocatable, private :: weights ! DIIS weights, G w = H
!
!     Mappings that keep track storage position for errors and parameters
!
!        - index:    DIIS matrix index
!        - record:   storage index
!
      integer, dimension(:), allocatable, private :: index_to_record
      integer, dimension(:), allocatable, private :: record_to_index
!
   contains
!
      procedure, public :: initialize &
                        => initialize_diis_tool
!
      procedure, public :: update &
                        => update_diis_tool
!
      procedure, public :: finalize &
                        => finalize_diis_tool
!
      procedure, private :: print_settings
!
      procedure, private :: calculate_diis_matrix_dimension
!
      procedure, private :: save_errors_and_parameters
      procedure, private :: modify_mappings_to_remove_oldest_record
!
      procedure, private :: update_diis_matrix
!
      procedure, private :: extract_old_G_elements
      procedure, private :: construct_new_G_elements
!
      procedure, private :: solve_diis_equation
!
      procedure, private :: construct_padded_G
      procedure, private :: construct_padded_H
!
      procedure, private :: extrapolate_parameters
!
      procedure, private :: construct_extrapolated_vector
      procedure, private :: update_crop_error_and_parameter
!
   end type diis_tool
!
!
   interface diis_tool
!
      procedure :: new_diis_tool
!
   end interface diis_tool
!
!
contains
!
!
   function new_diis_tool(name_, n_parameters, n_equations, &
                          dimension_, crop, records_in_memory) result(this)
!!
!!    New DIIS tool
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2018
!!
!!    name_:               Name of DIIS tool (used for temporary files)
!!
!!    n_parameters:        Number of parameters (dimensionality of the x vectors)
!!
!!    n_equations:         Number of equations (dimensionality of the e vectors)
!!
!!    dimension_:          Number of records to use in the DIIS extrapolation.
!!                         Default is to use the previous 8 records.
!!
!!    crop:                Use conjugate residual with optimal trial vectors
!!                         (see J. Chem. Theory Comput. 2015, 11, 4, 1518-1524). Optional.
!!                         Default is false.
!!
!!    records_in_memory:   Place records in memory. Default is false.
!!
      implicit none
!
      type(diis_tool) :: this
!
      character(len=*), intent(in) :: name_
!
      integer, intent(in) :: n_parameters
      integer, intent(in) :: n_equations
!
      integer, intent(in), optional :: dimension_
      logical, intent(in), optional :: crop
      logical, intent(in), optional :: records_in_memory
!
      logical :: records_in_memory_
!
      this%name_           = trim(name_)
      this%n_parameters    = n_parameters
      this%n_equations     = n_equations
 !
      this%iteration       = 1
      this%dimension_      = 8
!
      this%crop            = .false.
!
      if (present(dimension_)) then
!
         this%dimension_ = dimension_
!
      endif
!
      if (present(crop)) then
!
         this%crop = crop
!
      endif
!
      call this%print_settings()
!
      records_in_memory_ = .false.
      if (present(records_in_memory)) records_in_memory_ = records_in_memory
!
      this%errors = record_storer(trim(this%name_) // '_errors',  &
                                     this%n_equations,            &
                                     this%dimension_,             &
                                     records_in_memory_)
!
      this%parameters = record_storer(trim(this%name_) // '_parameters',   &
                                     this%n_parameters,                    &
                                     this%dimension_,                      &
                                     records_in_memory_)
!
   end function new_diis_tool
!
!
   subroutine print_settings(this)
!!
!!    Print settings
!!    Written by Sarai D. Folkestad, 2020
!!
      implicit none
!
      class(diis_tool), intent(in) :: this
!
      call output%printf('n', '- DIIS tool settings:', fs='(/t3,a)')
!
      call output%printf('n', 'DIIS dimension: (i3)', &
                         ints=[this%dimension_], fs='(/t6,a)')
!
      if (this%crop) then
!
         call output%printf('n', 'Enabled CROP in the DIIS algorithm.', fs='(/t6,a)')
!
      endif
!
      call output%newline('n')
!
   end subroutine print_settings
!
!
   subroutine initialize_diis_tool(this)
!!
!!    Initialize
!!    Written by Eirik F. Kjønstad, Nov 2019
!!
      implicit none
!
      class(diis_tool), intent(inout) :: this
!
      integer :: I
!
      call this%errors%initialize()
      call this%parameters%initialize()
!
      call mem%alloc(this%G, this%dimension_, this%dimension_)
!
      call mem%alloc(this%index_to_record, this%dimension_)
      call mem%alloc(this%record_to_index, this%dimension_)
!
      do I = 1, this%dimension_
!
         this%index_to_record(I) = I
         this%record_to_index(I) = I
!
      enddo
!
      call mem%alloc(this%weights, this%dimension_)
!
   end subroutine initialize_diis_tool
!
!
   subroutine finalize_diis_tool(this)
!!
!!    Finalize
!!    Written by Eirik F. Kjønstad, Nov 2019
!!
      implicit none
!
      class(diis_tool), intent(inout) :: this
!
      call this%errors%finalize()
      call this%parameters%finalize()
!
      call mem%dealloc(this%G, this%dimension_, this%dimension_)
!
      call mem%dealloc(this%index_to_record, this%dimension_)
      call mem%dealloc(this%record_to_index, this%dimension_)
!
      call mem%dealloc(this%weights, this%dimension_)
!
   end subroutine finalize_diis_tool
!
!
   subroutine update_diis_tool(this, e, x)
!!
!!    Update
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2021
!!
!!    e: current error
!!    x: current parameter
!!
!!    Determines the extrapolated DIIS vector <x> and sets x = <x>.
!!
!!    The <x> is a linear combination of the current and previous x's
!!    that minimizes an averaged error <e>.
!!
!!    To do this, the routine constructs the DIIS matrix G and solves the DIIS equation:
!!
!!       G w = H,   G_ij = e_i^T e_j,    H_i = 0
!!
!!    The G and H is padded with -1 and 0 to ensure normalization of the DIIS weights w.
!!
      implicit none
!
      class(diis_tool), intent(inout) :: this
!
      real(dp), dimension(this%n_equations), intent(in)     :: e
      real(dp), dimension(this%n_parameters), intent(inout) :: x
!
      type(timings), allocatable :: timer
!
      timer = timings('DIIS: update space and calculate extrapolated guess', 'v')
      call timer%turn_on()
!
      call this%calculate_diis_matrix_dimension()
!
      call this%save_errors_and_parameters(e, x)
!
      call this%update_diis_matrix(e)
!
      call this%extrapolate_parameters(x)
!
      if (this%crop) call this%update_crop_error_and_parameter(x)
!
      this%iteration = this%iteration + 1
!
      call timer%turn_off()
!
   end subroutine update_diis_tool
!
!
   pure subroutine calculate_diis_matrix_dimension(this)
!!
!!    Calculate DIIS matrix dimension
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
      implicit none
!
      class(diis_tool), intent(inout) :: this
!
      if (this%iteration .gt. this%dimension_) then
!
         this%dim_G = this%dimension_
!
      else
!
         this%dim_G = this%iteration - &
                 this%dimension_*((this%iteration-1)/this%dimension_)
!
      endif
!
   end subroutine calculate_diis_matrix_dimension
!
!
   subroutine save_errors_and_parameters(this, e, x)
!!
!!    Save errors and parameters
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
!!    Stores 'e' and 'x' as the current errors and parameters
!!
      implicit none
!
      class(diis_tool) :: this
!
      real(dp), dimension(this%n_equations), intent(in)  :: e
      real(dp), dimension(this%n_parameters), intent(in) :: x
!
      integer :: record
!
      if (this%iteration .gt. this%dimension_) then
!
         call this%modify_mappings_to_remove_oldest_record()
!
      endif
!
      record = this%index_to_record(this%dim_G)
!
      call this%errors%copy_record_in(e, record)
      call this%parameters%copy_record_in(x, record)
!
   end subroutine save_errors_and_parameters
!
!
   subroutine modify_mappings_to_remove_oldest_record(this)
!!
!!    Modify mappings to remove oldest record
!!    Written by Eirik F. Kjønstad, 2019
!!
!!    Reorders records as follows: ABCD -> BCDA
!!
!!    The oldest record (A) then appears last and will be overwritten next.
!!
      implicit none
!
      class(diis_tool), intent(inout) :: this
!
      integer :: I, first_record
!
      first_record = this%index_to_record(1)
!
      do I = 1, this%dimension_ - 1
!
         this%index_to_record(I) =  this%index_to_record(I + 1)
!
      enddo
!
      this%index_to_record(this%dimension_) = first_record
!
      do I = 1, this%dimension_
!
         this%record_to_index(this%index_to_record(I)) = I
!
      enddo
!
   end subroutine modify_mappings_to_remove_oldest_record
!
!
   subroutine update_diis_matrix(this, e)
!!
!!    Update DIIS matrix
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
!!    e: current error vector
!!
!!    Updates the DIIS matrix G by calculating the new elements (associated with 'e')
!!    and removing any old obsolete elements of G, where
!!
!!       G(i,j) = e_i^T e_j.
!!
      implicit none
!
      class(diis_tool), intent(inout) :: this
!
      real(dp), dimension(this%n_equations), intent(in) :: e
!
      real(dp), dimension(:,:), allocatable :: G
!
      call mem%alloc(G, this%dim_G, this%dim_G)
!
      if (this%dim_G .gt. 1) call this%extract_old_G_elements(G)
!
      call this%construct_new_G_elements(G,           &
                                         e,           &
                                         this%dim_G,  &
                                         this%dim_G)
!
      this%G(1 : this%dim_G, 1 : this%dim_G) = G(1 : this%dim_G, 1 : this%dim_G)
!
      call mem%dealloc(G, this%dim_G, this%dim_G)
!
   end subroutine update_diis_matrix
!
!
   subroutine extract_old_G_elements(this, G)
!!
!!    Extract old G elements
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none
!
      class(diis_tool), intent(in) :: this
!
      real(dp), dimension(this%dim_G, this%dim_G), intent(out) :: G
!
      if (this%iteration .gt. this%dimension_) then
!
!        Copy the parts that are associated with non-deleted records
!
         G(1 : this%dim_G - 1, &
           1 : this%dim_G - 1) = this%G(2 : this%dim_G, &
                                        2 : this%dim_G)

!
      else
!
!        Copy the whole matrix
!
         G(1 : this%dim_G - 1, &
           1 : this%dim_G - 1) = this%G(1 : this%dim_G - 1, &
                                        1 : this%dim_G - 1)
!
      endif
!
   end subroutine extract_old_G_elements
!
!
   subroutine construct_new_G_elements(this, G, e, element, dim_G)
!!
!!    Construct new G elements
!!    Written by Eirik F. Kjønstad, Mar 2020
!!
!!    Computes new elements of the DIIS matrix G:
!!
!!       G(i,element) = e_i^T e = G(element,i)          i = 1,2, ..., element
!!
!!    e: current error
!!
      implicit none
!
      class(diis_tool), intent(inout) :: this
!
      real(dp), dimension(this%n_equations), intent(in) :: e
!
      integer, intent(in) :: dim_G, element
!
      real(dp), dimension(dim_G, dim_G), intent(inout) :: G
!
      real(dp), dimension(:,:), pointer, contiguous :: prev_e
!
      integer :: i, current_i_batch, req_0, req_1, i_index
!
      type(batching_index), allocatable :: batch_i
!
      real(dp) :: ddot
!
      batch_i = batching_index(element)
!
      req_0 = 0
      req_1 = this%errors%required_to_load_record()
!
      call mem%batch_setup(batch_i, req_0, req_1, 'construct_new_G_elements')
!
      call this%errors%prepare_records([batch_i])
!
      do current_i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(current_i_batch)
!
         call this%errors%load(prev_e, batch_i)
!
         do i = batch_i%first, batch_i%get_last()
!
            i_index = this%record_to_index(i)
!
            G(element, i_index) = ddot(this%n_equations,                   &
                                       e,                                  &
                                       1,                                  &
                                       prev_e(:, i - batch_i%first + 1),   &
                                       1)
!
            G(i_index, element) = G(element, i_index)
!
         enddo
      enddo
!
      call this%errors%free_records()
!
      call mem%batch_finalize()
!
   end subroutine construct_new_G_elements
!
!
   subroutine extrapolate_parameters(this, x)
!!
!!    Extrapolate parameters
!!    Written by Eirik F. Kjønstad, 2018-2021
!!
!!    Extrapolates the parameter vector and stores the result in x.
!!
      implicit none
!
      class(diis_tool), intent(inout) :: this
!
      real(dp), dimension(this%n_parameters), intent(out) :: x
!
      call this%solve_diis_equation()
!
      call this%construct_extrapolated_vector(this%parameters, x, this%n_parameters)
!
   end subroutine extrapolate_parameters
!
!
   subroutine solve_diis_equation(this)
!!
!!    Solve DIIS equation
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018-2021
!!
!!    Solve the DIIS equation G w = H for the DIIS weights w
!!
      implicit none
!
      class(diis_tool), intent(inout) :: this
!
      integer, dimension(:), allocatable :: ipiv ! Pivot integers (see dgesv routine)
!
      real(dp), dimension(:), allocatable   :: padded_H
      real(dp), dimension(:,:), allocatable :: padded_G
!
      integer :: info, dim_
!
      dim_ = this%dim_G + 1
!
      call mem%alloc(padded_H, dim_)
      call mem%alloc(padded_G, dim_, dim_)
!
      call this%construct_padded_H(padded_H)
      call this%construct_padded_G(padded_G)
!
      call mem%alloc(ipiv, dim_)
      ipiv = 0
!
      call dgesv(dim_,        &
                  1,          &
                  padded_G,   &
                  dim_,       &
                  ipiv,       &
                  padded_H,   & ! = solution to equation on exit
                  dim_,       &
                  info)
!
      call mem%dealloc(ipiv, dim_)
!
      if (info .ne. 0) &
         call output%error_msg('could not solve DIIS matrix equation! (i0)', ints=[info])
!
      this%weights = padded_H(1 : this%dim_G)
!
      call mem%dealloc(padded_H, dim_)
      call mem%dealloc(padded_G, dim_, dim_)
!
   end subroutine solve_diis_equation
!
!
   subroutine construct_padded_H(this, padded_H)
!!
!!    Construct padded H
!!    Written by Eirik F. Kjønstad, 2021
!!
!!    padded_H = (0 .... 0 -1)
!!
      implicit none
!
      class(diis_tool), intent(in) :: this
!
      real(dp), dimension(this%dim_G + 1), intent(out) :: padded_H
!
      padded_H(1 : this%dim_G) = zero
      padded_H(this%dim_G + 1) = -one
!
   end subroutine construct_padded_H
!
!
   subroutine construct_padded_G(this, padded_G)
!!
!!    Construct padded G
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
!!    Equal to G but padded as follows:
!!
!!       G(dim_G + 1, j) = -1,   j = 1, 2, 3, ..., dim_G
!!       G(i, dim_G + 1) = -1,   i = 1, 2, 3, ..., dim_G
!!
!!       G(dim_G + 1, dim_G + 1) = 0.
!!
      implicit none
!
      class(diis_tool), intent(inout) :: this
!
      real(dp), dimension(this%dim_G + 1, this%dim_G + 1), intent(out) :: padded_G
!
      integer :: i
!
      padded_G(1 : this%dim_G, 1 : this%dim_G) = this%G(1 : this%dim_G, 1 : this%dim_G)
!
      do i = 1, this%dim_G
!
         padded_G(this%dim_G + 1, i) = -one
         padded_G(i, this%dim_G + 1) = -one
!
      enddo
!
      padded_G(this%dim_G + 1, this%dim_G + 1) = zero
!
   end subroutine construct_padded_G
!
!
   subroutine construct_extrapolated_vector(this, y_vectors, y, dim_y)
!!
!!    Construct extrapolated vector
!!    Written by Eirik F. Kjønstad, Mar 2020
!!
!!    Constructs the extrapolated vector
!!
!!       y = sum_i w_i y_i
!!
!!    Here
!!
!!       y:         vector to calculate
!!       y_vectors: storer for y_i
!!
      implicit none
!
      class(diis_tool), intent(in) :: this
!
      class(record_storer), intent(inout) :: y_vectors
!
      integer, intent(in) :: dim_y
!
      real(dp), dimension(dim_y), intent(out) :: y
!
      real(dp), dimension(:,:), pointer, contiguous :: y_i
!
      integer :: req_0, req_1, i, current_i_batch
!
      type(batching_index), allocatable :: batch_i
!
      call zero_array(y, dim_y)
!
      req_0 = 0
      req_1 = y_vectors%required_to_load_record()
!
      batch_i = batching_index(this%dim_G)
!
      call mem%batch_setup(batch_i, req_0, req_1, 'construct_extrapolated_vector')
!
      call y_vectors%prepare_records([batch_i])
!
      do current_i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(current_i_batch)
!
         call y_vectors%load(y_i, batch_i)
!
         do i = batch_i%first, batch_i%get_last()
!
!           Add w_i y_i to y
!
            call daxpy(dim_y,                                     &
                       this%weights(this%record_to_index(i)),     &
                       y_i(:, i - batch_i%first + 1),             &
                       1,  y, 1)
!
         enddo
!
      enddo
!
      call y_vectors%free_records()
!
      call mem%batch_finalize()
!
   end subroutine construct_extrapolated_vector
!
!
   subroutine update_crop_error_and_parameter(this, x)
!!
!!    Update CROP error and parameter
!!    Written by Eirik F. Kjønstad, Nov-Des 2019
!!
!!    Replaces the current e and x to the extrapolated < e > and < x >
!!
!!    On entry, 'x' is assumed to contain the DIIS-extrapolated parameters.
!!
      implicit none
!
      class(diis_tool), intent(inout) :: this
!
      real(dp), dimension(this%n_parameters), intent(in) :: x
!
      real(dp), dimension(:), allocatable :: e ! Extrapolated error vector < e >
!
      integer :: record
!
!     Overwrite x with <x>
!
      record = this%index_to_record(this%dim_G)
!
      call this%parameters%copy_record_in(x, record)
!
!     Compute <e> and overwrite current e with <e>
!
      call mem%alloc(e, this%n_equations)
!
      call this%construct_extrapolated_vector(this%errors, e, this%n_equations)
!
      call this%errors%copy_record_in(e, record)
!
!     Use <e> to update DIIS matrix
!
      call this%construct_new_G_elements(this%G,         &
                                         e,              &
                                         this%dim_G,     &
                                         this%dimension_)
!
      call mem%dealloc(e, this%n_equations)
!
   end subroutine update_crop_error_and_parameter
!
!
end module diis_tool_class
