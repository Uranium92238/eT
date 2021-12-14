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
module davidson_tool_class
!
!!
!!    Abstract Davidson davidson class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018-2019
!!
!
   use parameters
!
   use range_class
   use record_storer_class, only: record_storer
   use memory_manager_class, only: mem
   use global_out, only: output
   use array_utilities, only: get_l2_norm, copy_and_scale, zero_array
   use precondition_tool_class, only: precondition_tool
   use batching_index_class, only: batching_index
   use timings_class, only: timings
!
   type, abstract :: davidson_tool
!
      character(len=40) :: name_
!
      real(dp), dimension(:), allocatable :: omega_re ! omega = frequency or energy;
      real(dp), dimension(:), allocatable :: omega_im ! see eigen Davidson and linear Davidson tools
!
      real(dp), dimension(:,:), allocatable :: A_red
      real(dp), dimension(:,:), allocatable :: X_red
!
      class(record_storer), allocatable :: trials, transforms
!
      integer :: dim_red
      integer :: max_dim_red
      integer :: n_new_trials
!
      integer :: n_parameters
      integer :: n_solutions
!
      real(dp) :: lindep_threshold
!
      logical :: do_precondition
      class(precondition_tool), allocatable :: preconditioner
!
   contains
!
!     Procedures a user of the tool may need to use
!
      procedure :: set_trial            => set_trial_davidson_tool
      procedure :: get_trial            => get_trial_davidson_tool
      procedure :: add_new_trial        => add_new_trial_davidson_tool
      procedure :: set_transform        => set_transform_davidson_tool
      procedure :: get_transform        => get_transform_davidson_tool
      procedure :: first_new_trial      => first_new_trial_davidson_tool
      procedure :: last_new_trial       => last_new_trial_davidson_tool
!
      procedure :: set_preconditioner   => set_preconditioner_davidson_tool
!
      procedure :: red_dim_exceeds_max  => red_dim_exceeds_max_davidson_tool
      procedure :: update_reduced_dim   => update_reduced_dim_davidson_tool
!
      procedure :: update_reduced_space &
                => update_reduced_space_davidson_tool
!
      procedure :: print_settings       => print_settings_davidson_tool
!
!     Other routines
!
      procedure :: construct_solution                          &
                => construct_solution_davidson_tool
!
      procedure :: construct_AX                                &
                => construct_AX_davidson_tool
!
      procedure :: construct_full_space_vector                 &
                => construct_full_space_vector_davidson_tool
!
      procedure :: construct_reduced_matrix                    &
                => construct_reduced_matrix_davidson_tool
!
      procedure :: construct_reduced_submatrix                 &
                => construct_reduced_submatrix_davidson_tool
!
      procedure :: orthonormalize_trial_vecs                   &
                => orthonormalize_trial_vecs_davidson_tool
!
      procedure :: set_trials_to_solutions                     &
                => set_trials_to_solutions_davidson_tool
!
      procedure :: initialize                                  &
                => initialize_davidson_tool
!
      procedure :: cleanup                                     &
                => cleanup_davidson_tool
!
      procedure(destruct_reduced_space_quantities), deferred :: destruct_reduced_space_quantities
!
      procedure :: set_lindep_threshold                        &
                => set_lindep_threshold_davidson_tool
!
      procedure :: reset_reduced_space                         &
                => reset_reduced_space_davidson_tool
!
   end type davidson_tool
!
!
   abstract interface
!
      subroutine destruct_reduced_space_quantities(davidson)
!
         import :: davidson_tool
!
         implicit none
!
         class(davidson_tool) :: davidson
!
      end subroutine destruct_reduced_space_quantities
!
   end interface
!
!
contains
!
!
   subroutine initialize_davidson_tool(davidson)
!!
!!    Initialize
!!    Written by Eirik F. Kjønstad, Nov 2019
!!
!!    Initializes the storers for trials and transforms and allocates the frequencies/energies
!!
!!    records_in_memory: if true,  trials and transforms are stored in memory
!!                       if false, trials and transforms are stored on disk
!!
!!    Trial vectors c define the subspace. Transforms refer to the transformed
!!    trial vectors, i.e. rho = A c, where A is the linear transformation.
!!
      implicit none
!
      class(davidson_tool) :: davidson
!
      call davidson%reset_reduced_space()
!
      call davidson%trials%initialize()
      call davidson%transforms%initialize()
!
      call mem%alloc(davidson%omega_re, davidson%max_dim_red)
      call mem%alloc(davidson%omega_im, davidson%max_dim_red)
!
      call zero_array(davidson%omega_re, davidson%max_dim_red)
      call zero_array(davidson%omega_im, davidson%max_dim_red)
!
   end subroutine initialize_davidson_tool
!
!
   subroutine cleanup_davidson_tool(davidson)
!!
!!    Cleanup
!!    Written by Eirik F. Kjønstad, Nov 2019
!!
!!    Cleanup of variables (and associated memory and/or files)
!!    initialized in the 'initialize' routine.
!!
      implicit none
!
      class(davidson_tool) :: davidson
!
      call davidson%trials%finalize()
      call davidson%transforms%finalize()
!
      call mem%dealloc(davidson%omega_re, davidson%max_dim_red)
      call mem%dealloc(davidson%omega_im, davidson%max_dim_red)
!
      call davidson%destruct_reduced_space_quantities()
      if (davidson%do_precondition) call davidson%preconditioner%destruct_precondition_vector()
!
   end subroutine cleanup_davidson_tool
!
!
   subroutine set_trial_davidson_tool(davidson, c, n)
!!
!!    Set trial
!!    Written by Eirik F. Kjønstad, 2019
!!
!!    Stores the nth trial.
!!
      implicit none
!
      class(davidson_tool) :: davidson
!
      integer, intent(in) :: n
!
      real(dp), dimension(davidson%n_parameters), intent(in) :: c
!
      call davidson%trials%copy_record_in(c, n)
!
   end subroutine set_trial_davidson_tool
!
!
   subroutine get_trial_davidson_tool(davidson, c, n)
!!
!!    Get trial
!!    Written by Eirik F. Kjønstad, 2019
!!
!!    Retrieves the nth trial.
!!
      implicit none
!
      class(davidson_tool) :: davidson
!
      integer, intent(in) :: n
!
      real(dp), dimension(davidson%n_parameters), intent(out) :: c
!
      call davidson%trials%copy_record_out(c, n)
!
   end subroutine get_trial_davidson_tool
!
!
   subroutine set_transform_davidson_tool(davidson, rho, n)
!!
!!    Set transform
!!    Written by Eirik F. Kjønstad, 2019
!!
!!    Stores the nth transform.
!!
      implicit none
!
      class(davidson_tool) :: davidson
!
      integer, intent(in) :: n
!
      real(dp), dimension(davidson%n_parameters), intent(in) :: rho
!
      call davidson%transforms%copy_record_in(rho, n)
!
   end subroutine set_transform_davidson_tool
!
!
   subroutine get_transform_davidson_tool(davidson, rho, n)
!!
!!    Get transform
!!    Written by Eirik F. Kjønstad, 2019
!!
!!    Retrieves the nth transform.
!!
      implicit none
!
      class(davidson_tool) :: davidson
!
      integer, intent(in) :: n
!
      real(dp), dimension(davidson%n_parameters), intent(out) :: rho
!
      call davidson%transforms%copy_record_out(rho, n)
!
   end subroutine get_transform_davidson_tool
!
!
   function first_new_trial_davidson_tool(davidson) result(first)
!!
!!    First trial
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2019
!!
!!    Returns index of first trial vector.
!!
      implicit none
!
      class(davidson_tool), intent(in) :: davidson
!
      integer :: first
!
      first = davidson%dim_red - davidson%n_new_trials + 1
!
   end function first_new_trial_davidson_tool
!
!
   function last_new_trial_davidson_tool(davidson) result(last)
!!
!!    Last trial
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2019
!!
!!    Returns index of last trial vector.
!!
      implicit none
!
      class(davidson_tool), intent(in) :: davidson
!
      integer :: last
!
      last = davidson%dim_red
!
   end function last_new_trial_davidson_tool
!
!
   subroutine add_new_trial_davidson_tool(davidson, R, n)
!!
!!    Add new trial
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018-2019
!!
!!    R is the nth residual
!!
!!    Preconditions R, normalizes it, and adds it to the trial space
!!
      implicit none
!
      class(davidson_tool) :: davidson
!
      real(dp), dimension(davidson%n_parameters), intent(in) :: R
!
      integer, intent(in) :: n
!
      real(dp) :: norm_trial
      real(dp), dimension(:), allocatable :: trial
!
!     Precondition & normalize R
!
      call mem%alloc(trial, davidson%n_parameters)
      call dcopy(davidson%n_parameters, R, 1, trial, 1)
!
      if (davidson%do_precondition) &
         call davidson%preconditioner%do_(trial,                        &
                                          shift=davidson%omega_re(n),   &
                                          prefactor=-one)
!
      norm_trial = get_l2_norm(trial, davidson%n_parameters)
      call dscal(davidson%n_parameters, one/norm_trial, trial, 1)
!
!     Add it to the trial space
!
      davidson%n_new_trials = davidson%n_new_trials + 1
      call davidson%set_trial(trial, davidson%dim_red + davidson%n_new_trials)
!
      call mem%dealloc(trial, davidson%n_parameters)
!
   end subroutine add_new_trial_davidson_tool
!
!
   subroutine read_trial_davidson_tool(davidson, c, n)
!!
!!    Read trial vector
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2018
!!
!!    Reads the nth trial vector from file and places it in c.
!!
!!    If n is not passed, it reads the trial at wherever the cursor
!!    in the file currently is.
!!
      implicit none
!
      class(davidson_tool), intent(inout) :: davidson
!
      real(dp), dimension(davidson%n_parameters), intent(out) :: c
!
      integer, intent(in) :: n
!
      call davidson%trials%copy_record_out(c, n)
!
   end subroutine read_trial_davidson_tool
!
!
   subroutine write_trial_davidson_tool(davidson, c, n)
!!
!!    Write trial
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2018
!!
!!    Write nth trial vector c to file.
!!
      implicit none
!
      class(davidson_tool), intent(inout) :: davidson
!
      real(dp), dimension(davidson%n_parameters), intent(in) :: c
!
      integer, intent(in) :: n
!
      call davidson%trials%copy_record_in(c, n)
!
   end subroutine write_trial_davidson_tool
!
!
   subroutine orthonormalize_trial_vecs_davidson_tool(davidson)
!!
!!    Orthonormalize trial vecs
!!    Written by Eirik F. Kjønstad, Oct 2019 and Mar 2020
!!
!!    Orthonormalizes the new trial vectors against the existing
!!    trial vectors using the modified Gram-Schmidt (MGS) procedure.
!!
      implicit none
!
      class(davidson_tool), intent(inout) :: davidson
!
      real(dp), dimension(:,:), pointer, contiguous :: c_i, c_j
!
      real(dp) :: norm_c_j, r_ji, ddot
!
      integer :: i, j, k, req_0, req_1_i, req_1_j, req_2
!
      integer :: n_trials, n_new_trials, n_old_trials
!
      integer :: current_i_batch, current_j_batch
!
      type(batching_index), allocatable :: batch_i, batch_j
!
      type(range_), allocatable :: i_interval
!
      integer, parameter :: n_orthonormalizations = 2
!
      type(timings), allocatable :: timer
!
      logical :: perform_batching
!
      timer = timings('Davidson: orthonormalization of new trial vectors', 'v')
      call timer%turn_on()
!
!     Batch setup over trial vectors
!
      perform_batching = .true. ! Set to false in special case of no batch-setup
!
      req_0   = 0
      req_1_i = davidson%trials%required_to_load_record()
      req_1_j = davidson%trials%required_to_load_record()
      req_2   = 0
!
      n_trials     = davidson%last_new_trial()
      n_new_trials = davidson%last_new_trial() - davidson%first_new_trial() + 1
      n_old_trials = n_trials - n_new_trials
!
      if ((n_new_trials + &
           min(n_old_trials, 1))*(davidson%n_parameters)*dp .lt. mem%get_available()) then
!
!        Can keep all new vectors and at least one old vector
!        => hold all new vectors (j) & batch over old vectors (i)
!
         call output%printf('debug', 'Davidson orthogonalization: &
                                     &batching over old vectors')
!
         batch_j = batching_index(n_new_trials,             &
                                  offset=n_old_trials)
!
         call batch_j%do_single_batch()
!
         batch_i = batching_index(dimension_=n_old_trials,  &
                                  offset=0)
!
         if (batch_i%index_dimension .eq. 0) then
!
!           No old vectors; no batches over i
!
            call batch_i%do_not_batch()
            perform_batching = .false.
!
         else
!
!           Set up batching over old vectors i
!
            req_0 = req_1_j * n_new_trials
            call mem%batch_setup(batch_i, req_0, req_1_i, 'orthonormalize_trial_vecs 1')
!
         endif
!
      else
!
!        Cannot hold all new and a single old vector
!        => batch over both new (j) and old vectors (i)
!
         call output%printf('debug', 'Davidson orthogonalization: &
                                     &batching over old and new vectors')
!
         batch_i = batching_index(dimension_=n_trials,   &
                                  offset=0)
!
         batch_j = batching_index(n_new_trials,          &
                                  offset=n_old_trials)
!
         call mem%batch_setup(batch_i, batch_j, req_0, req_1_i, req_1_j, req_2, &
                              'orthonormalize_trial_vecs 2')
!
      endif
!
      call davidson%trials%prepare_records([batch_j, batch_i])
!
!     Modified Gram-Schmidt loops
!
      do k = 1, n_orthonormalizations
!
         do current_j_batch = 1, batch_j%num_batches
!
            call batch_j%determine_limits(current_j_batch)
!
            call davidson%trials%load(c_j, batch_j, 1)
!
!           Remove components along vectors i < batch of j
!
            do current_i_batch = 1, batch_i%num_batches
!
               call batch_i%determine_limits(current_i_batch)
!
!              Load only i vectors that are needed for the j batch
!
               i_interval = range_(batch_i%first,                                &
                                   min(batch_i%get_last(), batch_j%first - 1) -  &
                                   batch_i%first + 1)
!
               if (i_interval%length .lt. 1) cycle ! No i vectors to load;
                                                         ! go to next i batch
!
               call davidson%trials%load(c_i, i_interval, 2)
!
!              Remove i components along j vectors
!
               do j = batch_j%first, batch_j%get_last()
                  do i = i_interval%first, i_interval%get_last()
!
                     r_ji = ddot(davidson%n_parameters,              &
                                 c_j(:, j - batch_j%first + 1),      &
                                 1,                                  &
                                 c_i(:, i - i_interval%first + 1),   &
                                 1)
!
                     call daxpy(davidson%n_parameters,            &
                                -r_ji,                            &
                                c_i(:, i - i_interval%first + 1), &
                                1,                                &
                                c_j(:, j - batch_j%first + 1),    &
                                1)
!
                  enddo
               enddo
            enddo ! end of i batches
!
!           Remove components within the j batch and normalize
!
            do j = 1, batch_j%length
!
               do i = 1, j - 1
!
                  r_ji = ddot(davidson%n_parameters,  &
                              c_j(:, j),              &
                              1,                      &
                              c_j(:, i),              &
                              1)
!
                  call daxpy(davidson%n_parameters,   &
                             -r_ji,                   &
                             c_j(:, i),               &
                             1,                       &
                             c_j(:, j),               &
                             1)
!
               enddo
!
               norm_c_j = get_l2_norm(c_j(:, j), davidson%n_parameters)
!
               if (norm_c_j < davidson%lindep_threshold) &
                  call output%error_msg('detected linear dependence &
                                           &in Davidson trial space.')
!
               call dscal(davidson%n_parameters, one/norm_c_j, c_j(:, j), 1)
!
            enddo
!
!           Make sure the normalized batch of j vectors are stored
!
            call davidson%trials%store(c_j, batch_j)
!
         enddo ! end of j batches
      enddo ! end of normalization loops
!
      if (perform_batching) call mem%batch_finalize()
!
      call davidson%trials%free_records()
!
      call timer%turn_off()
!
   end subroutine orthonormalize_trial_vecs_davidson_tool
!
!
   subroutine construct_reduced_matrix_davidson_tool(davidson)
!!
!!    Construct reduced matrix
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2018 - Mar 2020
!!
!!    Constructs the reduced matrix
!!
!!       A_ij = c_i^T A c_j
!!            = c_i^T rho_j.
!!
!!    Elements computed in previous iterations are not recomputed.
!!
      implicit none
!
      class(davidson_tool), intent(inout) :: davidson
!
      real(dp), dimension(:,:), allocatable :: A_red_copy
!
      integer :: req_0, req_1_i, req_1_j, req_2
!
      type(batching_index), allocatable :: batch_i, batch_j
!
      integer :: prev_dim_red
!
      type(timings), allocatable :: timer
!
      timer = timings('Davidson: time to construct reduced matrix', 'v')
      call timer%turn_on()
!
!     Will batch over reduced space indices, i for trials and j for transforms
!     Required is always the same:
!
      req_0   = 0
      req_1_i = davidson%trials%required_to_load_record()
      req_1_j = davidson%transforms%required_to_load_record()
      req_2   = 0
!
      if (davidson%dim_red .eq. davidson%n_solutions) then
!
!        First iteration or reset of space: calculate all elements
!
         call mem%alloc(davidson%A_red, davidson%dim_red, davidson%dim_red)
!
         batch_i = batching_index(davidson%dim_red)
         batch_j = batching_index(davidson%dim_red)
!
         call mem%batch_setup(batch_i, batch_j, req_0, req_1_i, req_1_j, req_2, &
                              'construct_reduced_matrix 1')
!
         call davidson%construct_reduced_submatrix(batch_i, batch_j)
!
         call mem%batch_finalize()
!
      else
!
!        Not first iteration: calculate new elements only
!
!        Transfer previously calculated elements and
!        extend the dimensionality of the reduced matrix
!
         prev_dim_red = davidson%dim_red - davidson%n_new_trials
!
         call mem%alloc(A_red_copy, prev_dim_red, prev_dim_red)
!
         call dcopy(prev_dim_red**2, davidson%A_red, 1, A_red_copy, 1)
!
         call mem%dealloc(davidson%A_red, prev_dim_red, prev_dim_red)
!
         call mem%alloc(davidson%A_red, davidson%dim_red, davidson%dim_red)
!
         davidson%A_red(1:prev_dim_red, 1:prev_dim_red) = A_red_copy
!
         call mem%dealloc(A_red_copy, prev_dim_red, prev_dim_red)
!
!        Compute the new blocks in the reduced matrix
!
!        i index new, j index all
!
         batch_i = batching_index(dimension_=davidson%n_new_trials, &
                                  offset=prev_dim_red)
!
         batch_j = batching_index(dimension_=davidson%dim_red, &
                                  offset=0)
!
         call mem%batch_setup(batch_i, batch_j, req_0, req_1_i, req_1_j, req_2, &
                              'construct_reduced_matrix 2')
!
         call davidson%construct_reduced_submatrix(batch_i, batch_j)
!
         call mem%batch_finalize()
!
!        i index old, j index new
!
         batch_i = batching_index(dimension_=prev_dim_red, &
                                  offset=0)
!
         batch_j = batching_index(dimension_=davidson%n_new_trials, &
                                  offset=prev_dim_red)
!
         call mem%batch_setup(batch_i, batch_j, req_0, req_1_i, req_1_j, req_2, &
                              'construct_reduced_matrix 3')
!
         call davidson%construct_reduced_submatrix(batch_i, batch_j)
!
         call mem%batch_finalize()
!
      endif
!
      call timer%turn_off()
!
   end subroutine construct_reduced_matrix_davidson_tool
!
!
   subroutine construct_reduced_submatrix_davidson_tool(davidson, batch_i, batch_j)
!!
!!    Construct reduced submatrix
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018-2020
!!
!!    Constructs parts of the reduced matrix A_ij based on the prepared
!!    batching indices batch_i and batch_j.
!!
      implicit none
!
      class(davidson_tool), intent(inout) :: davidson
!
      type(batching_index), intent(inout) :: batch_i, batch_j
!
      integer :: current_i_batch, current_j_batch
!
      real(dp), dimension(:,:), pointer, contiguous :: c, rho
!
      call davidson%trials%prepare_records([batch_i])
      call davidson%transforms%prepare_records([batch_j])
!
      do current_i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(current_i_batch)
!
         call davidson%trials%load(c, batch_i)
!
         do current_j_batch = 1, batch_j%num_batches
!
            call batch_j%determine_limits(current_j_batch)
!
            call davidson%transforms%load(rho, batch_j)
!
            call dgemm('T', 'N',                                        &
                        batch_i%length,                                 &
                        batch_j%length,                                 &
                        davidson%n_parameters,                          &
                        one,                                            &
                        c,                                              &
                        davidson%n_parameters,                          &
                        rho,                                            &
                        davidson%n_parameters,                          &
                        zero,                                           &
                        davidson%A_red(batch_i%first, batch_j%first),   &
                        davidson%dim_red)
!
         enddo ! end of j batches
      enddo ! end of i batches
!
      call davidson%trials%free_records()
      call davidson%transforms%free_records()
!
   end subroutine construct_reduced_submatrix_davidson_tool
!
!
   subroutine construct_solution_davidson_tool(davidson, X, n)
!!
!!    Construct solution
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2018
!!
!!    Constructs
!!
!!       X_n = sum_i c_i (X_red_n)_i,
!!
!!    i.e. the current nth full space solution.
!!
      implicit none
!
      class(davidson_tool), intent(inout) :: davidson
!
      integer, intent(in) :: n
!
      real(dp), dimension(davidson%n_parameters), intent(out) :: X
!
      call davidson%construct_full_space_vector(davidson%trials, X, n)
!
   end subroutine construct_solution_davidson_tool
!
!
   subroutine construct_full_space_vector_davidson_tool(davidson, y_vectors, y, n)
!!
!!    Construct full space vector
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2018
!!
!!    Constructs
!!
!!       X_n = sum_i y_i (X_red_n)_i,
!!
!!    i.e. the current nth full space vector.
!!
!!    y_vectors is a storer containing the full space y vectors.
!!
!!    This is typically the trials or transforms, in which case this
!!    constructs the full space X_n and A X_n.
!!
      implicit none
!
      class(davidson_tool), intent(inout) :: davidson
!
      class(record_storer) :: y_vectors
!
      integer, intent(in) :: n
!
      real(dp), dimension(davidson%n_parameters), intent(out) :: y
!
      real(dp), dimension(:,:), pointer, contiguous :: y_i
!
      integer :: current_i_batch, req_0, req_1
!
      type(batching_index), allocatable :: batch_i
!
      if (n .gt. davidson%dim_red) then
!
        call output%error_msg('Tried to construct a full space vector (n = (i0)) that does not &
                              &correspond to a reduced space solution (dim_red = (i0)).', &
                              ints=[n, davidson%dim_red])
!
      endif
!
      req_0 = 0
      req_1 = y_vectors%required_to_load_record()
!
      batch_i = batching_index(davidson%dim_red)
!
      call mem%batch_setup(batch_i, req_0, req_1, 'construct_full_space_vector')
!
      call zero_array(y, davidson%n_parameters)
!
      call y_vectors%prepare_records([batch_i])
!
      do current_i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(current_i_batch)
!
         call y_vectors%load(y_i, batch_i)
!
         call dgemm('N', 'N',                         &
                     davidson%n_parameters,           &
                     1,                               &
                     batch_i%length,                  &
                     one,                             &
                     y_i,                             & ! y_alpha,i
                     davidson%n_parameters,           &
                     davidson%X_red(batch_i%first,n), & ! Xred_i,n
                     davidson%dim_red,                &
                     one,                             &
                     y,                               & ! y_alpha,n
                     davidson%n_parameters)
!
      enddo
!
      call mem%batch_finalize()
!
      call y_vectors%free_records()
!
   end subroutine construct_full_space_vector_davidson_tool
!
!
   subroutine construct_AX_davidson_tool(davidson, AX, n)
!!
!!    Construct AX
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2018
!!
!!    Constructs
!!
!!       AX_n = sum_i rho_i (X_red_n)_i,
!!
!!    where X_n is the current nth full space solution.
!!
      implicit none
!
      class(davidson_tool), intent(inout) :: davidson
!
      integer, intent(in) :: n
!
      real(dp), dimension(davidson%n_parameters), intent(out) :: AX
!
      call davidson%construct_full_space_vector(davidson%transforms, AX, n)
!
   end subroutine construct_AX_davidson_tool
!
!
   subroutine set_preconditioner_davidson_tool(davidson, preconditioner)
!!
!!    Set preconditioner
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    This routine saves the diagonal preconditioner to file. The
!!    assumption being that
!!
!!       preconditioner(i) ~ A(i,i),
!!
!!    where A is the coefficient matrix. The inverse of this diagonal
!!    matrix is a good approximation of A if A is diagonally dominant
!!    to some extent.
!!
      implicit none
!
      class(davidson_tool), intent(inout) :: davidson
!
      real(dp), dimension(davidson%n_parameters), intent(in) :: preconditioner
!
      davidson%do_precondition = .true.
      davidson%preconditioner = precondition_tool(davidson%n_parameters)
      call davidson%preconditioner%initialize_and_set_precondition_vector(preconditioner)
!
   end subroutine set_preconditioner_davidson_tool
!
!
   subroutine set_trials_to_solutions_davidson_tool(davidson)
!!
!!    Set trials to solutions
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    Set trials equal to solution on file.
!!    This should be used when max_dim_red.
!!
      implicit none
!
      class(davidson_tool), intent(inout) :: davidson
!
      real(dp), dimension(:,:), allocatable :: X
!
      integer :: solution
!
      call mem%alloc(X, davidson%n_parameters, davidson%n_solutions)
!
      do solution = 1, davidson%n_solutions
!
         call davidson%construct_solution(X(:,solution), solution)
!
      enddo
!
      call davidson%trials%copy_record_in(X, range_(1, davidson%n_solutions))
!
      call mem%dealloc(X, davidson%n_parameters, davidson%n_solutions)
!
!     Deallocate reduced space quantities (if they are allocated)
!
      call davidson%destruct_reduced_space_quantities()
!
      davidson%dim_red = 0
      davidson%n_new_trials = davidson%n_solutions
!
   end subroutine set_trials_to_solutions_davidson_tool
!
!
   subroutine update_reduced_dim_davidson_tool(davidson)
!!
!!    Update reduced dim
!!    Written by Eirik F. Kjønstad, Mar 2020
!!
!!    Updates the reduced dimensionality:
!!
!!       dim_red = dim_red + n_new_trials
!!
      implicit none
!
      class(davidson_tool), intent(inout) :: davidson
!
      davidson%dim_red = davidson%dim_red + davidson%n_new_trials
!
   end subroutine update_reduced_dim_davidson_tool
!
!
   subroutine update_reduced_space_davidson_tool(davidson)
!!
!!    Update reduced space
!!    Written by Sarai D. Folkestad, May 2021
!!
      implicit none
!
      class(davidson_tool), intent(inout) :: davidson
!
      if (davidson%red_dim_exceeds_max()) call davidson%set_trials_to_solutions()
      call davidson%update_reduced_dim()
      call davidson%orthonormalize_trial_vecs()
!
   end subroutine update_reduced_space_davidson_tool
!
!
   function red_dim_exceeds_max_davidson_tool(davidson) result(exceeds_max)
!!
!!    Reduced dimension exceeds max
!!    Written by Eirik F. Kjønstad, Mar 2020
!!
!!    Is the next reduced space dimensionality greater than or
!!    equal than the maximum? Then return true. Otherwise false.
!!
      implicit none
!
      class(davidson_tool), intent(in) :: davidson
!
      logical :: exceeds_max
!
      exceeds_max = (davidson%dim_red + davidson%n_new_trials .gt. davidson%max_dim_red)
!
   end function red_dim_exceeds_max_davidson_tool
!
!
   subroutine print_settings_davidson_tool(davidson)
!!
!!    Print settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
      implicit none
!
      class(davidson_tool), intent(inout) :: davidson
!
      call output%printf('n', ' - Davidson tool settings:' , fs='(/t3,a)')
!
      call output%printf('n', 'Number of parameters:          (i11)', &
                         ints=[davidson%n_parameters], fs='(/t6,a)')
!
      call output%printf('n', 'Number of requested solutions: (i11)', &
                         ints=[davidson%n_solutions], fs='(t6,a)')
!
      call output%printf('n', 'Max reduced space dimension:   (i11)', &
                         ints=[davidson%max_dim_red], fs='(t6,a)')
!
      call output%newline('n')
!
    end subroutine print_settings_davidson_tool
!
!
    subroutine set_lindep_threshold_davidson_tool(davidson, lindep_threshold)
!!
!!    Set lindep threshold
!!    Written by Sarai D. Folkestad, 2021
!!
      implicit none
!
      class(davidson_tool), intent(inout) :: davidson
!
      real(dp), intent(in) :: lindep_threshold
!
      davidson%lindep_threshold = lindep_threshold
!
    end subroutine set_lindep_threshold_davidson_tool
!
!
   subroutine reset_reduced_space_davidson_tool(davidson)
!!
!!    Reset reduced space
!!    Written by Sarai D. Folkestad, 2021
!!
      implicit none
!
      class(davidson_tool), intent(inout) :: davidson
!
      davidson%dim_red      = 0
      davidson%n_new_trials = davidson%n_solutions
!
   end subroutine reset_reduced_space_davidson_tool
!
!
end module davidson_tool_class
