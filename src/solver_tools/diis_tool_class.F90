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
module diis_tool_class
!
!!
!!    DIIS tool class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2018
!!
!!    Objects of this class can be instantiated in order to solve a non-linear set of equations,
!!    e(x) = 0, where e and x are vectors of length n_equations and n_parameters.
!!
!!    The squared norm of the averaged error vector
!!
!!       < e(x) > = sum_k w_k e(x)_k
!!
!!    is minimized with respect to the condition sum_k w_k = 1. The k index
!!    runs over previous iterations, by default limited to the last 8 iterations
!!    (determined by dimension_).
!!
!!    The standard update of the parameters is
!!
!!       x = sum_k=1^n w_k x_k     (*)
!!
!!    Note that the "x_k" sent to the routine may be quasi-Newton update estimates,
!!    so that the DIIS-extrapolated x in (*) is the optimal combination of previous 
!!    quasi-Newton estimates (this is, e.g., the case for CC ground and 
!!    excited states).
!!
!!    See Pulay, P. Convergence acceleration of iterative sequences. 
!!    The case of SCF iteration. Chem. Phys. Lett. 1980, 73, 393−398.
!!
!!    Typical usage: 
!!
!!    do while (.not. converged) 
!!
!!       -> Calculate error vector e given the current parameters x
!!
!!       -> Overwrite x with the x_n you want to use in the extrapolation 
!!          step (*). This can, e.g., be a quasi-Newton update estimate based  
!!          on the current e and x. 
!!
!!       call diis%update(e, x)
!!
!!       -> now, x contains the extrapolated guess (*)
!!
!!    enddo 
!!
!!    In the tool, the DIIS matrix is denoted as G:
!!
!!       G_ij = e_i^T e_j,     
!!
!!    where i and j denote the record numbers.
!!
!
   use global_out, only : output
!
   use parameters
   use array_utilities, only : zero_array
!
   use sequential_file_class, only : sequential_file
   use record_storer_class, only : record_storer  
!
   use memory_manager_class, only : mem
   use batching_index_class, only: batching_index
   use timings_class, only: timings
!
!
   type :: diis_tool
!
      character(len=40) :: name_          ! determines the prefix of all DIIS files
!
      class(record_storer), allocatable :: e_vectors, x_vectors 
!
      integer, private :: iteration       ! Variable keeping track of the current DIIS iteration
                                          ! Note: defined to increment by +1 each time 'update' 
                                          ! is called.
!
      integer, private :: dimension_      ! Standard is 8, though it might be useful to 
                                          ! change this value
!
      integer, private :: n_parameters    ! The length of the x vectors
      integer, private :: n_equations     ! The length of the e vectors
!
      logical :: accumulate               ! Variable keeping track of wether we are  
!                                         ! cumulatively building the DIIS matrix or not
!                                         ! Note: default is to accumulate
!
      logical :: erase_history            ! (DIIS space) 1,2,3,4,...,8,1,2,3,... if true 
                                          ! (DIIS space) 1,2,3,4,...,8,8,8,8,... if false 
!
      logical :: crop                     ! Conjugate residual with optimal trial vectors 
                                          ! J. Chem. Theory Comput. 2015, 11, 4, 1518-1524
!
      real(dp), dimension(:,:), allocatable :: G ! DIIS matrix 
!
!     Index arrays 
!
!        - index:    index in the DIIS matrix according to DIIS tool 
!        - record:   index in record storer, i.e. storage position 
!
!        These begin to differ when old records are discarded from DIIS space. 
!        See cycle_left routine.
!
      integer, dimension(:), allocatable :: index_to_record
      integer, dimension(:), allocatable :: record_to_index 
!
   contains
!
      procedure :: update                          => update_diis_tool   
!
      procedure :: read_x                          => read_x_diis_tool
      procedure :: read_e                          => read_e_diis_tool
!
      procedure :: write_x                         => write_x_diis_tool
      procedure :: write_e                         => write_e_diis_tool
!
      procedure :: initialize_storers              => initialize_storers_diis_tool
      procedure :: finalize_storers                => finalize_storers_diis_tool
!
      procedure :: get_dim_G                       => get_dim_G_diis_tool
      procedure :: get_dimension                   => get_dimension_diis_tool
!
      procedure :: print_settings                  => print_settings_diis_tool
!
      procedure, private :: construct_padded_G              &
                         => construct_padded_G_diis_tool
!
      procedure, private :: construct_new_G_elements        &
                         => construct_new_G_elements_diis_tool
!
      procedure, private :: construct_full_G                &
                         => construct_full_G_diis_tool
!
      procedure, private :: write_e_and_x                   &
                         => write_e_and_x_diis_tool
!
      procedure, private :: update_crop_e_and_x             &
                         => update_crop_e_and_x_diis_tool 
!
      procedure, private :: cycle_left                      &
                         => cycle_left_diis_tool
!
      procedure, private :: construct_extrapolated_vector   &
                         => construct_extrapolated_vector_diis_tool
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
                  dimension_, accumulate, erase_history, crop, &
                  records_in_memory) result(diis)
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
!!    accumulate:          If true, the parts of the DIIS matrix constructed 
!!                         on "update" is the parts that have not been calculated
!!                         previously. If false, the entire DIIS matrix is re-
!!                         constructed when "update" is called. This is sometimes
!!                         neccessary because the previous records have to be 
!!                         updated to a new basis. Default is true.
!!
!!    erase_history:       Instead of keeping the previous 8 records 
!!                         (if DIIS dimension = 8), the history is erased 
!!                         every 8th call to "update". This usually leads to 
!!                         slower convergene (or oscillating behavior in the 
!!                         residuals). Default is false.
!!
!!    crop:                Use conjugate residual with optimal trial vectors 
!!                         (see J. Chem. Theory Comput. 2015, 11, 4, 1518-1524). Optional.
!!                         Default is false. 
!!
      type(diis_tool) :: diis
!
      character(len=*), intent(in) :: name_
!
      integer, intent(in) :: n_parameters
      integer, intent(in) :: n_equations
!
      integer, intent(in), optional :: dimension_
      logical, intent(in), optional :: accumulate
      logical, intent(in), optional :: erase_history
      logical, intent(in), optional :: crop
      logical, intent(in), optional :: records_in_memory
!
      diis%name_           = trim(name_)
      diis%n_parameters    = n_parameters
      diis%n_equations     = n_equations
 !
      diis%iteration       = 1 
      diis%dimension_      = 8
!
      diis%accumulate      = .true. 
      diis%erase_history   = .false. 
      diis%crop            = .false. 
!
      if (present(accumulate)) then 
!
         diis%accumulate = accumulate
!
      endif
!
      if (present(erase_history)) then
!
         diis%erase_history = erase_history
!
      endif
!
      if (present(dimension_)) then
!
         diis%dimension_ = dimension_
!
      endif 
!
      if (present(crop)) then 
!
         diis%crop = crop
!
      endif 
!
      if (diis%crop .and. diis%erase_history) then 
!
         call output%warning_msg("Asked for both 'crop' and 'erase history' in DIIS tool. " // &
                                 "Convergence will not improve by using CROP in this case.")
!
      endif
!
      call diis%print_settings()
!
      diis%e_vectors = record_storer(trim(diis%name_) // '_e',    &
                                     diis%n_equations,            &
                                     diis%dimension_,             &
                                     records_in_memory,           &
                                     delete=.true.)
!
      diis%x_vectors = record_storer(trim(diis%name_) // '_x',    &
                                     diis%n_parameters,           &
                                     diis%dimension_,             &
                                     records_in_memory,           &
                                     delete=.true.)
!
   end function new_diis_tool
!
!
   subroutine initialize_storers_diis_tool(diis)
!!
!!    Prepare storers 
!!    Written by Eirik F. Kjønstad, Nov 2019 
!!
!!    Prepares the storer objects for e and x vectors.
!!
!!    records_in_memory:   If true, keep DIIS records (previous x and e vectors) 
!!                         in memory. Otherwise, they are stored on disk.
!!
      implicit none 
!
      class(diis_tool) :: diis
!
      integer :: I
!
      call diis%e_vectors%initialize_storer()
      call diis%x_vectors%initialize_storer()
!
      call mem%alloc(diis%G, diis%dimension_, diis%dimension_)
!
      call mem%alloc(diis%index_to_record, diis%dimension_)
      call mem%alloc(diis%record_to_index, diis%dimension_)
!
      do I = 1, diis%dimension_
!
         diis%index_to_record(I) = I 
         diis%record_to_index(I) = I 
!
      enddo
!
   end subroutine initialize_storers_diis_tool
!
!
   subroutine finalize_storers_diis_tool(diis)
!!
!!    Cleanup storers 
!!    Written by Eirik F. Kjønstad, Nov 2019 
!!
!!    Performs cleanup for the different storers. 
!!    Cannot be done in the constructor. 
!!
      implicit none 
!
      class(diis_tool) :: diis
!
      call diis%e_vectors%finalize_storer()
      call diis%x_vectors%finalize_storer()
!
      call mem%dealloc(diis%G, diis%dimension_, diis%dimension_)
!
      call mem%dealloc(diis%index_to_record, diis%dimension_)
      call mem%dealloc(diis%record_to_index, diis%dimension_)
!
   end subroutine finalize_storers_diis_tool
!
!
    subroutine update_diis_tool(diis, e, x)
!!
!!    Update 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2018
!!
!!    Routine to call in each iteration of a DIIS-based solver. The 
!!    routine stores the current vectors (e and x) sent to the routine. 
!!    Then it constructs the DIIS matrix G and solves the DIIS equation 
!!
!!       G w = H,   G_ij = e_i^T e_j,    H_i = 0 
!!
!!    for the DIIS weights w. In practice, G and H are padded 
!!    (one extra column and row) to enforce 
!!
!!       sum_i w_i = 1.
!!
!!    On exit, 
!!
!!       x = sum_i w_i x_i.
!!
!!    Based on the CCSD specific routine written by Sarai D. Folkestad
!!    and Eirik F. Kjønstad, May 2017.
!!
      implicit none
!
      class(diis_tool) :: diis
!
      real(dp), dimension(diis%n_equations)  :: e
      real(dp), dimension(diis%n_parameters) :: x
!
      integer :: info
      integer :: dim_G
!
      real(dp), dimension(:), allocatable   :: padded_H
      real(dp), dimension(:,:), allocatable :: padded_G
!
      integer, dimension(:), allocatable :: ipiv ! Pivot integers (see dgesv routine)
!
      type(timings), allocatable :: timer 
!
      timer = timings('DIIS: total time in DIIS-update', 'v')
      call timer%turn_on()
!
!     Compute the current dimensionality of the problem 
!     (1, 2,..., 7, 8, 8, 8, 8,...) for the standard dimension_ = 8 without erasing history
!
      dim_G = diis%get_dim_G()
!
!     Write current e and x vector to file
!
      call diis%write_e_and_x(e, x, dim_G)
!
!     :: Solve the least squares problem, G w = H
!
!        G : DIIS matrix, G_ij = e_i^T e_j,
!        H : DIIS vector,  H_i = 0,
!
!     where i, j = 1, 2, ..., dim_G. To enforce normality
!     of the solution, G is extended with a row & column of -1's
!     and H with a -1 at the end.
!
!     Set the DIIS vector H
!
      call mem%alloc(padded_H, dim_G + 1)
!
      padded_H = zero
      padded_H(dim_G + 1) = -one
!
!     Set the padded DIIS matrix G
!
      call mem%alloc(padded_G, dim_G + 1, dim_G + 1)
!
      call diis%construct_padded_G(dim_G, padded_G, e)
!
!     Solve the DIIS equation G w = H for the DIIS weights w 
!
!     Note: on exit, the solution is in the padded_H,
!     provided info = 0 (see LAPACK documentation for more)
!
      call mem%alloc(ipiv, dim_G + 1)
      ipiv = 0
!
      call dgesv(dim_G + 1,         &
                  1,                &
                  padded_G,         &
                  dim_G + 1,        &
                  ipiv,             &
                  padded_H,         &
                  dim_G + 1,        &
                  info)
!
      call mem%dealloc(ipiv, dim_G + 1)
      call mem%dealloc(padded_G, dim_G + 1, dim_G + 1)
!
      if (info .ne. 0) &
         call output%error_msg('could not solve DIIS matrix equation! (i0)', ints=[info])
!
!     Construct the extrapolated x vector, i.e. the DIIS guess for the parameters 
!
      call diis%construct_extrapolated_vector(diis%x_vectors, x, padded_H, dim_G)
!
!     If CROP is enabled, we replace the current e and x vectors 
!     with the extrapolated e and x 
!
      if (diis%crop) call diis%update_crop_e_and_x(x,             &
                                                   padded_H,      &
                                                   dim_G)
!
      call mem%dealloc(padded_H, dim_G + 1)
!
      diis%iteration = diis%iteration + 1
!
      call timer%turn_off()
!
   end subroutine update_diis_tool
!
!
   subroutine update_crop_e_and_x_diis_tool(diis, x, H, dim_G)
!!
!!    Update CROP e and x 
!!    Written by Eirik F. Kjønstad, Nov-Des 2019 
!!
!!    Sets the current e and x to the extrapolated < e > and < x >
!!
!!    x:       extrapolated x, i.e. < x > 
!!    dim_G:   the current dimension of the DIIS matrix
!!    H:       current DIIS vector
!!
      implicit none 
!
      class(diis_tool) :: diis 
!
      real(dp), dimension(diis%n_parameters), intent(in) :: x 
!
      integer, intent(in) :: dim_G
!
      real(dp), dimension(dim_G), intent(in) :: H
!
      real(dp), dimension(:), allocatable :: e ! Extrapolated error vector < e >
!
!     Overwrite current x with the DIIS extrapolated < x > 
! 
      call diis%x_vectors%copy_record_in(x, diis%index_to_record(dim_G))
!
!     Compute extrapolated < e > and overwrite current e record with it 
!
      call mem%alloc(e, diis%n_equations)
!
      call diis%construct_extrapolated_vector(diis%e_vectors, e, H, dim_G)
      call diis%e_vectors%copy_record_in(e, diis%index_to_record(dim_G)) 
!
      if (diis%accumulate) then 
!
!        Recompute the new elements of the DIIS matrix (given the new e vector)
!        If not accumulate, these will be computed in the next iteration anyway.
!
         call diis%construct_new_G_elements(diis%G,         &
                                            e,              &
                                            dim_G,          &
                                            diis%dimension_)
!
      endif
!
      call mem%dealloc(e, diis%n_equations)
!
   end subroutine update_crop_e_and_x_diis_tool
!
!
   subroutine construct_extrapolated_vector_diis_tool(diis, y_vectors, y, w, dim_w)
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
!!       y :         extrapolated vector
!!       y_vectors : record storer for keeping y_i 
!!       w:          vector with DIIS weights w_i   
!!
      implicit none 
!
      class(diis_tool) :: diis 
!
      class(record_storer) :: y_vectors
!
      real(dp), dimension(y_vectors%record_dim), intent(out) :: y 
!
      integer, intent(in) :: dim_w 
!
      real(dp), dimension(dim_w), intent(in) :: w 
!
      real(dp), dimension(:,:), pointer, contiguous :: y_i 
!
      integer :: req_0, req_1, i, current_i_batch
!
      type(batching_index), allocatable :: batch_i 
!
      call zero_array(y, y_vectors%record_dim)
!
!     Prepare to batch over records i 
!
      req_0 = 0
      req_1 = y_vectors%required_to_load_record()
!
      batch_i = batching_index(dim_w)
!
      call mem%batch_setup(batch_i, req_0, req_1)
!
!     Construct extrapolated vector in batches 
!
      call y_vectors%prepare_records([batch_i])
!
      do current_i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(current_i_batch)
!
         call y_vectors%load(y_i, batch_i)
!
         do i = batch_i%first, batch_i%last 
!
!           Add w_i y_i to y 
!
            call daxpy(y_vectors%record_dim,          &
                       w(diis%record_to_index(i)),    &
                       y_i(:, i - batch_i%first + 1), &
                       1,                             &
                       y,                             &
                       1)
!
         enddo
!
      enddo ! end of i batches 
!
      call y_vectors%free_records()      
!
   end subroutine construct_extrapolated_vector_diis_tool
!
!
   function get_dim_G_diis_tool(diis)
!!
!!    Get current dimension 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
      implicit none
!
      class(diis_tool) :: diis
!
      integer :: get_dim_G_diis_tool
!
      if (diis%iteration .gt. diis%dimension_ .and. &
            .not. diis%erase_history) then 
!
         get_dim_G_diis_tool = diis%dimension_
!
      else
!
         get_dim_G_diis_tool = diis%iteration - &
                           diis%dimension_*((diis%iteration-1)/diis%dimension_) 
!
      endif
!
   end function get_dim_G_diis_tool
!
!
   function get_dimension_diis_tool(diis)
!!
!!    Get dimension 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
      implicit none
!
      class(diis_tool) :: diis
!
      integer :: get_dimension_diis_tool
!
      get_dimension_diis_tool = diis%dimension_
!
   end function get_dimension_diis_tool
!
!
   subroutine write_e_and_x_diis_tool(diis, e, x, dim_G)
!!
!!    Write current e and x 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
!!    Stores the current e and x vectors. Each index 1, 2, ..., diis_dim
!!    has its own record. If more iterations have passed than we keep on file
!!    (i.e., diis%iteration > diis%dimension_), then we do a first-in/
!!    last-out replacement unless erase_history is true. 
!!
      implicit none 
!
      class(diis_tool) :: diis
!
      real(dp), dimension(diis%n_equations), intent(in)  :: e
      real(dp), dimension(diis%n_parameters), intent(in) :: x
!
      integer, intent(in) :: dim_G
!
      logical :: cycle_left   
!
      cycle_left = (.not. diis%erase_history) .and. &
                  (dim_G .eq. diis%dimension_) .and. &
                 (dim_G .ne. diis%iteration)
!
      if (cycle_left) call diis%cycle_left()
!
      call diis%e_vectors%copy_record_in(e, diis%index_to_record(dim_G))
      call diis%x_vectors%copy_record_in(x, diis%index_to_record(dim_G))
!
   end subroutine write_e_and_x_diis_tool
!
!
   subroutine construct_padded_G_diis_tool(diis, dim_G, padded_G, e)
!!
!!    Construct padded G 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
!!    Constructs the padded DIIS matrix G. For i, j = 1, 2, .., dim_G, we have
!!
!!       G(i,j) = e_i^T e_j.
!!
!!    The padding consists of:
!!
!!       G(dim_G + 1, j) = -1,   j = 1, 2, 3, ..., dim_G 
!!       G(i, dim_G + 1) = -1,   i = 1, 2, 3, ..., dim_G 
!!
!!       G(dim_G + 1, dim_G + 1) = 0.
!!
      implicit none
!
      class(diis_tool) :: diis
!
      integer, intent(in) :: dim_G
!
      real(dp), dimension(dim_G + 1, dim_G + 1), intent(inout) :: padded_G
!
      real(dp), dimension(diis%n_equations), intent(in) :: e
!
      integer :: i
!
      padded_G = zero
!
      if (.not. diis%accumulate) then 
!
!        - Construct the entire DIIS matrix G_ij = e_i^T e_j 
!
         call diis%construct_full_G(padded_G, dim_G)
!
      else ! accumulate 
!
!        - Transfer the parts of the DIIS matrix already calculated 
!
         if (diis%iteration .gt. diis%dimension_ .and. .not. diis%erase_history) then 
!
!           Old matrix is full size; cut out the parts to keep 
!
            padded_G(1 : dim_G - 1, &
                     1 : dim_G - 1) = diis%G(2 : dim_G, &
                                             2 : dim_G)   

!
         elseif (dim_G .gt. 1) then 
!
!           Just copy the old matrix to the new 
!
            padded_G(1 : dim_G - 1, &
                     1 : dim_G - 1) = diis%G(1 : dim_G - 1, &
                                             1 : dim_G - 1)   
!
         endif 
!
!        - Compute the new elements of the DIIS matrix 
!
         call diis%construct_new_G_elements(padded_G,   &
                                            e,          &
                                            dim_G,      &
                                            dim_G + 1)
!
!        Keep a copy of the DIIS matrix for the next iteration 
!
         diis%G(1 : dim_G, 1 : dim_G) = padded_G(1 : dim_G, 1 : dim_G)
!
      endif
!
!     Add the (-1)-padding in the last row and column  
!
      do i = 1, dim_G 
!
         padded_G(dim_G + 1, i) = -one
         padded_G(i, dim_G + 1) = -one
!
      enddo
!
   end subroutine construct_padded_G_diis_tool
!
!
   subroutine construct_full_G_diis_tool(diis, G, dim_G)
!!
!!    Construct full G  
!!    Written by Eirik F. Kjønstad, Mar 2020 
!! 
!!    Computes the DIIS matrix
!!
!!       G(i,j) = e_i^T e_j = G(j,i)          i,j = 1,2, ..., dim_G
!!
!!    G:     G matrix 
!!    dim_G: dimensionality of G matrix 
!!
      implicit none 
!
      class(diis_tool) :: diis 
!
      integer, intent(in) :: dim_G 
!
      real(dp), dimension(dim_G + 1, dim_G + 1), intent(inout) :: G 
!
      real(dp), dimension(:,:), pointer, contiguous :: prev_e_i, prev_e_j  
!
      integer :: i, j, current_i_batch, current_j_batch, req_0, req_1_i, req_1_j, req_2
!
      integer :: i_index, j_index 
!
      type(batching_index), allocatable :: batch_i, batch_j 
!
      real(dp) :: ddot
!
!     Prepare for batching over i and j 
!
      req_0 = 0
!
      req_1_i = diis%e_vectors%required_to_load_record()
      req_1_j = diis%e_vectors%required_to_load_record()
!
      req_2 = 0
!
      batch_i = batching_index(dim_G)
      batch_j = batching_index(dim_G)
!
      call mem%batch_setup(batch_i, batch_j, req_0, req_1_i, req_1_j, req_2)
!
!     Allocate arrays to store e_i and e_j
!
      call diis%e_vectors%prepare_records([batch_i, batch_j])
!
      do current_i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(current_i_batch)
!
         call diis%e_vectors%load(prev_e_i, batch_i, 1)
!
         do current_j_batch = 1, current_i_batch
!
            call batch_j%determine_limits(current_j_batch)
!
            call diis%e_vectors%load(prev_e_j, batch_j, 2)
!
!           Compute lower triangle contributions to the DIIS matrix 
!
            do i = batch_i%first, batch_i%last 
               do j = batch_j%first, min(batch_j%last, i) 
!
                  i_index = diis%record_to_index(i)
                  j_index = diis%record_to_index(j)
!
                  G(i_index, j_index) = ddot(diis%n_equations,                        &
                                             prev_e_i(:, i - batch_i%first + 1),      &
                                             1,                                       &
                                             prev_e_j(:, j - batch_j%first + 1),      &
                                             1)
!
               enddo
            enddo 
!
         enddo ! end of j batches 
      enddo ! end of i batches 
!
!     Set upper triangle of DIIS matrix 
!
      do i = 1, dim_G
         do j = 1, i
!
            G(j, i) = G(i, j)
!
         enddo
      enddo 
!
      call diis%e_vectors%free_records()
!
   end subroutine construct_full_G_diis_tool
!
!
   subroutine construct_new_G_elements_diis_tool(diis,      &
                                                 G,         &
                                                 e,         &
                                                 element,   &
                                                 dim_G)
!!
!!    Construct new G elements 
!!    Written by Eirik F. Kjønstad, Mar 2020 
!! 
!!    Computes new elements of the DIIS matrix:
!!
!!       G(i,element) = e_i^T e = G(element,i)          i = 1,2, ..., element
!!
!!    e :  current DIIS error vector 
!!    G:   G matrix 
!!
      implicit none 
!
      class(diis_tool) :: diis 
!
      real(dp), dimension(diis%n_equations), intent(in) :: e 
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
!     Prepare to batch over records 
!
      batch_i = batching_index(element)
!
      req_0 = 0
      req_1 = diis%e_vectors%required_to_load_record()
!
      call mem%batch_setup(batch_i, req_0, req_1)
!
!     Construct new elements in DIIS matrix 
!
      call diis%e_vectors%prepare_records([batch_i])
!
      do current_i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(current_i_batch)
!
         call diis%e_vectors%load(prev_e, batch_i)
!
         do i = batch_i%first, batch_i%last
!  
            i_index = diis%record_to_index(i)
!
            G(element, i_index) = ddot(diis%n_equations,                      &
                                       e,                                     &
                                       1,                                     &
                                       prev_e(:, i - batch_i%first + 1),      &
                                       1)
!
            G(i_index, element) = G(element, i_index)
!
         enddo
!
      enddo
!
      call diis%e_vectors%free_records()
!
   end subroutine construct_new_G_elements_diis_tool
!
!
   subroutine read_x_diis_tool(diis, x_i, i)
!!
!!    Read x_i 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Mar 2019
!!
      implicit none
!
      class(diis_tool) :: diis
!
      real(dp), dimension(diis%n_parameters), intent(out) :: x_i
!
      integer, intent(in) :: i
!
      integer :: dim_G
!
      dim_G = diis%get_dim_G()
!
      if (i .gt. dim_G) call output%error_msg('Asked for an x not in use.')
!
      call diis%x_vectors%copy_record_out(x_i, diis%index_to_record(i))
!
   end subroutine read_x_diis_tool
!
!
   subroutine read_e_diis_tool(diis, e_i, i)
!!
!!    Read e_i
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Mar 2019
!!
      implicit none
!
      class(diis_tool) :: diis
!
      real(dp), dimension(diis%n_equations), intent(out) :: e_i
!
      integer, intent(in) :: i
!
      integer :: dim_G
!
      dim_G = diis%get_dim_G()
!
      if (i .gt. dim_G) call output%error_msg('Asked for an e not in use.')
!
      call diis%e_vectors%copy_record_out(e_i, diis%index_to_record(i))
!
   end subroutine read_e_diis_tool
!
!
   subroutine write_x_diis_tool(diis, x_i, i)
!!
!!    Write x_i 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Mar 2019
!!
      implicit none
!
      class(diis_tool) :: diis
!
      real(dp), dimension(diis%n_parameters), intent(in) :: x_i
!
      integer, intent(in) :: i
!
      integer :: dim_G
!
      dim_G = diis%get_dim_G()
!
      if (diis%accumulate) call output%error_msg('Can not set x for cumulative diis.')
      if (i .gt. dim_G) call output%error_msg('Asked for an x not in use.')
!
      call diis%x_vectors%copy_record_in(x_i, diis%index_to_record(i))
!
   end subroutine write_x_diis_tool
!
!
   subroutine write_e_diis_tool(diis, e_i, i)
!!
!!    Write e_i 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Mar 2019
!!
      implicit none
!
      class(diis_tool) :: diis
!
      real(dp), dimension(diis%n_equations), intent(in) :: e_i
!
      integer, intent(in) :: i
!
      integer :: dim_G
!
      dim_G = diis%get_dim_G()
!
      if (diis%accumulate) call output%error_msg('Can not set e for cumulative diis.')
      if (i .gt. dim_G) call output%error_msg('Asked for an e not in use.')
!
      call diis%e_vectors%copy_record_in(e_i, diis%index_to_record(i))
!
   end subroutine write_e_diis_tool
!
!
   subroutine cycle_left_diis_tool(diis)
!!
!!    Cycle left 
!!    Written by Eirik F. Kjønstad, 2019 
!!
!!    If we have a list of records A, B, C, D, this routine 
!!    reorders the list to B, C, D, A. Then, when overwriting 
!!    the last record (i.e. A), we remove the originally oldest 
!!    record and insert the newest record at the end. 
!!
      implicit none
!
      class(diis_tool) :: diis
!
      integer :: I, first_record 
!
!     Update record to storage position index array 
!
      first_record = diis%index_to_record(1)
!
      do I = 1, diis%dimension_ - 1
!
         diis%index_to_record(I) =  diis%index_to_record(I + 1)
!
      enddo 
!
      diis%index_to_record(diis%dimension_) = first_record
!
!     Update storage position to record index array
!
      do I = 1, diis%dimension_
!
         diis%record_to_index(diis%index_to_record(I)) = I
!
      enddo 
!
   end subroutine cycle_left_diis_tool
!
!
   subroutine print_settings_diis_tool(diis)
!!
!!    Print settings
!!    Written by Sarai D. Folkestad, 2020
!!
      implicit none
!
      class(diis_tool), intent(inout)  :: diis
!
      call output%printf('n', '- DIIS tool settings:', fs='(/t3,a)')
!
      call output%printf('n', 'DIIS dimension: (i3)', &
                         ints=[diis%dimension_], fs='(/t6,a)')
!
      if (diis%crop) then 
!
        call output%printf('n', 'Enabled CROP in the DIIS algorithm.', fs='(/t6,a)')
!
      endif
!
      if (diis%erase_history) then 
!
        call output%printf('v', 'DIIS history deleted when &
              &full DIIS dimension is reached.', fs='(/t6,a)')
!
      endif
!
      call output%newline('n')
!
   end subroutine print_settings_diis_tool
!
!
end module diis_tool_class
