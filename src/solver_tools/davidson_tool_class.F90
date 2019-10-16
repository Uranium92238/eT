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
module davidson_tool_class
!
!!
!!    Abstract Davidson davidson class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018-2019
!!
!
   use parameters
!
   use sequential_file_class, only: sequential_file
   use memory_manager_class, only: mem 
   use global_out, only: output
   use array_utilities, only: get_l2_norm, copy_and_scale, zero_array
!
!
   type, abstract :: davidson_tool
!
      character(len=40) :: name_ 
!
      real(dp), dimension(:,:), allocatable :: A_red
      real(dp), dimension(:,:), allocatable :: X_red
!
      type(sequential_file) :: trials, transforms, preconditioner
!
      integer :: dim_red
      integer :: max_dim_red
      integer :: n_new_trials
!
      integer :: n_parameters
      integer :: n_solutions
!
      real(dp) :: add_trial_threshold 
!
      logical :: do_precondition
!
   contains
!
!     Procedures a user of the tool may need to use 
!
      procedure, non_overridable :: iterate                    => iterate_davidson_tool
      procedure, non_overridable :: read_trial                 => read_trial_davidson_tool
      procedure, non_overridable :: write_trial                => write_trial_davidson_tool      
      procedure, non_overridable :: write_transform            => write_transform_davidson_tool  
      procedure, non_overridable :: construct_solution         => construct_solution_davidson_tool
!
      procedure, non_overridable :: first_trial                => first_trial_davidson_tool
      procedure, non_overridable :: last_trial                 => last_trial_davidson_tool
!
      procedure :: set_preconditioner                          => set_preconditioner_davidson_tool
!
      procedure(solve_reduced_problem), deferred :: solve_reduced_problem
!
!     Other routines 
!
      procedure, non_overridable :: read_transform             => read_transform_davidson_tool
!
      procedure, non_overridable :: construct_AX               => construct_AX_davidson_tool
      procedure, non_overridable :: construct_reduced_matrix   => construct_reduced_matrix_davidson_tool
!
      procedure :: precondition                                => precondition_davidson_tool
!
      procedure :: orthogonalize_against_trial_vecs            => orthogonalize_against_trial_vecs_davidson_tool
      procedure :: orthonormalize_trial_vecs                   => orthonormalize_trial_vecs_davidson_tool
!
      procedure :: set_trials_to_solutions                     => set_trials_to_solutions_davidson_tool
!
   end type davidson_tool
!
!
   abstract interface
!
!
      subroutine solve_reduced_problem(davidson)
!
         import :: davidson_tool
!
         implicit none 
!
         class(davidson_tool) :: davidson 
!
      end subroutine solve_reduced_problem
!
!
   end interface
!
contains
!
!
   subroutine iterate_davidson_tool(davidson)
!!
!!    Iterate 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2019 
!!
!!    Updates the reduced space dimension. It should be called at 
!!    the beginning of the iterative loop, after the initial set 
!!    of trial vectors have been written:
!!
!!    do while (.not. converged ...)
!!
!!       iteration = iteration + 1
!!       call davidson%iterate()
!!
!!       ... 
!!
!!    enddo
!!
!!    If the current reduced dimension exceeds the specified max_dim_red,
!!    the solutions are set as existing trials, with the current set of trials 
!!    obtained from residuals as additional new trials vectors.
!!
      implicit none 
!
      class(davidson_tool), intent(inout) :: davidson 
!
      if (davidson%dim_red >= davidson%max_dim_red) then
!
         call davidson%set_trials_to_solutions()
!
      else
!
         davidson%dim_red = davidson%dim_red + davidson%n_new_trials
!
      endif
!
   end subroutine iterate_davidson_tool
!
!
   function first_trial_davidson_tool(davidson) result(first)
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
   end function first_trial_davidson_tool
!
!
   function last_trial_davidson_tool(davidson) result(last)
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
   end function last_trial_davidson_tool
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
      integer, optional, intent(in) :: n
!
      call davidson%trials%open_('read')
!
      if (present(n)) then
! 
         call davidson%trials%rewind_()
         call davidson%trials%skip(n - 1)
!
      endif
!
      call davidson%trials%read_(c, davidson%n_parameters)
      call davidson%trials%close_()
!
   end subroutine read_trial_davidson_tool
!
!
   subroutine write_trial_davidson_tool(davidson, c, position_)
!!
!!    Write trial
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2018
!!
!!    Write trial vector c to file.
!!
!!    Optional argument "position_" must be either 'rewind' or 'append',
!!    and default is 'append'.
!!
      implicit none
!
      class(davidson_tool), intent(inout) :: davidson
!
      real(dp), dimension(davidson%n_parameters), intent(in) :: c
!
      character(len=*), optional, intent(in) :: position_
!
      character(len=:), allocatable :: local_position_
!
      local_position_ = 'append'
!
!     Was position_ passed?
!
      if (present(position_)) then
!
!        Sanity check on position_ variable
!
         if ((trim(position_) .ne. 'rewind') .and. (trim(position_) .ne. 'append')) then
!
            call output%error_msg('position_ specifier not recognized in write_trial_davidson_tool.')
!
         endif
!
         local_position_ = position_
!
      endif 
!
      call davidson%trials%open_('write', local_position_)
      call davidson%trials%write_(c, davidson%n_parameters)
      call davidson%trials%close_()
!
   end subroutine write_trial_davidson_tool
!
!
   subroutine orthonormalize_trial_vecs_davidson_tool(davidson)
!!
!!    Orthogonalize trial vecs  
!!    Written by Eirik F. Kjønstad, Mar 2019 
!!
!!    Orthonormalizes trial vectors to make sure that the 
!!    metric in the reduced problem, c_i^T c_j = delta_ij,
!!    which is always assumed. This is needed after adding 
!!    an initial set of trial vectors, upon restart from 
!!    previous solutions and when resetting the reduced 
!!    space using the last set of solutions.
!!
      implicit none 
!
      class(davidson_tool), intent(inout) :: davidson 
!
      real(dp), dimension(:,:), allocatable :: c
!
      real(dp), dimension(:), allocatable :: c_tmp
!
      real(dp) :: norm_c, ddot, projection 
!
      integer :: i, j, n_done
!
      call output%printf('Orthonormalizing trial vectors.', pl='n', fs='(/t3,a)')
!
!     Write the first trial vector to file 
!
      call mem%alloc(c, davidson%n_parameters, davidson%dim_red)
      call mem%alloc(c_tmp, davidson%n_parameters)
!
      n_done = 0 
!
      do i = 1, davidson%dim_red 
!
         call davidson%read_trial(c(:,i), i)
!
         c_tmp = c(:,i)
         do j = 1, n_done 
!
            projection = ddot(davidson%n_parameters, c_tmp, 1, c(1,j), 1)
            call daxpy(davidson%n_parameters, -projection, c(:,j), 1, c(:,i), 1)
!
         enddo
!
         norm_c = sqrt(ddot(davidson%n_parameters, c(1,i), 1, c(1,i), 1))
         c(:,i) = c(:,i)/norm_c 
!
         n_done = n_done + 1
!
      enddo
!
      call mem%dealloc(c_tmp, davidson%n_parameters)
!
      call davidson%trials%open_('write', 'rewind')
!
      do i = 1, davidson%dim_red
!
         call davidson%trials%write_(c(:,i), davidson%n_parameters)
!
      enddo 
!
      call davidson%trials%close_()
!
      call mem%dealloc(c, davidson%n_parameters, davidson%dim_red)
!
   end subroutine orthonormalize_trial_vecs_davidson_tool
!
!
   subroutine orthogonalize_against_trial_vecs_davidson_tool(davidson, R)
!!
!!    Orthogonalize against trial vectors 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    Orthogonalizes R against the existing trial vectors. Note 
!!    that this is usually done residual after residual, where 
!!    each new residual is orthogonalized against the previous 
!!    (where the number of previous changes with 'n_new_trials').
!!
!!    Note that it does not normalize R.
!!
      implicit none 
!
      class(davidson_tool), intent(inout) :: davidson 
!
      real(dp), dimension(davidson%n_parameters) :: R 
!
      real(dp) :: ddot, projection_of_R_on_c_i
!
      integer :: i 
!
      real(dp), dimension(:), allocatable :: c_i
!
      call mem%alloc(c_i, davidson%n_parameters)
!
      do i = 1, davidson%dim_red + davidson%n_new_trials
!
         call davidson%read_trial(c_i, i)
         projection_of_R_on_c_i = ddot(davidson%n_parameters, c_i, 1, R, 1)
!
         call daxpy(davidson%n_parameters, -projection_of_R_on_c_i, c_i, 1, R, 1)
!
      enddo 
!
      call mem%dealloc(c_i, davidson%n_parameters)
!
   end subroutine orthogonalize_against_trial_vecs_davidson_tool
!
!
   subroutine read_transform_davidson_tool(davidson, rho_n, n)
!!
!!    Read transformed vector
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2018
!!
!!    Read nth transformed vector rho_n from file.
!!
      implicit none
!
      class(davidson_tool), intent(inout) :: davidson
!
      real(dp), dimension(davidson%n_parameters), intent(out) :: rho_n
!
      integer, intent(in) :: n
!
      call davidson%transforms%open_('read', 'rewind')
      call davidson%transforms%skip(n - 1)
      call davidson%transforms%read_(rho_n, davidson%n_parameters)
      call davidson%transforms%close_()
!
   end subroutine read_transform_davidson_tool
!
!
   subroutine write_transform_davidson_tool(davidson, rho, position_)
!!
!!    Write transform
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2018
!!
!!    Write transformed vector rho to file.
!!
!!    Optional argument position_ must be either 'rewind' or 'append',
!!    and default is 'append'.
!!
      implicit none
!
      class(davidson_tool), intent(inout) :: davidson
!
      real(dp), dimension(davidson%n_parameters), intent(in) :: rho
!
      character(len=*), optional, intent(in) :: position_ 
!
      character(len=:), allocatable :: local_position_
!
      local_position_ = 'append'
!
!     Was position_ passed?
!
      if (present(position_)) then
!
!        Sanity check on position_ variable
!
         if ((trim(position_) .ne. 'rewind') .and. (trim(position_) .ne. 'append')) then
!
            call output%error_msg('position_ specifier not recognized in write_transform_davidson_tool.')
!
         endif
!
         local_position_ = position_
!
      endif 
!
      call davidson%transforms%open_('write', local_position_)
      call davidson%transforms%write_(rho, davidson%n_parameters)
      call davidson%transforms%close_()
!
   end subroutine write_transform_davidson_tool
!
!
   subroutine construct_reduced_matrix_davidson_tool(davidson, entire)
!!
!!    Construct reduced matrix
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2018
!!
!!    Constructs the reduced matrix 
!!
!!       A_red_ij = c_i^T * A c_j
!!                = c_i^T * rho_j,
!!      
!!    by reading the trial vectors c_i and transformed vectors rho_j from file and writes
!!    to file.
!!
      implicit none
!
      class(davidson_tool), intent(inout) :: davidson
!
      logical, optional, intent(in) :: entire ! Construct the entire matrix 
!
      real(dp) :: ddot
!
      logical :: entire_local
!
      real(dp), dimension(:,:), allocatable :: A_red_copy
!
      real(dp), dimension(:), allocatable :: c_i, rho_j
!
      integer :: i, j
!
      if (present(entire)) then
!
         entire_local = entire
!
      else
!
         entire_local = .false.
!
      endif
!
!     Open the files
!
      call davidson%trials%open_('read', 'rewind')
      call davidson%transforms%open_('read', 'rewind')
!
!     Allocate c and rho
!
      call mem%alloc(c_i, davidson%n_parameters)
      call mem%alloc(rho_j, davidson%n_parameters)
!
!     Construct reduced matrix: A_red_ij = c_i^T * A * c_j = c_i^T * rho_j
!
      if (davidson%dim_red .eq. davidson%n_solutions .or. entire_local) then ! First iteration
!
!        Make the entire reduced matrix
!
         call mem%alloc(davidson%A_red, davidson%dim_red, davidson%dim_red)
!
         do i = 1, davidson%dim_red
!
            call davidson%trials%read_(c_i, davidson%n_parameters)
!
            call davidson%transforms%rewind_()
!
            do j = 1, davidson%dim_red
!
               call davidson%transforms%read_(rho_j, davidson%n_parameters)
!
               davidson%A_red(i,j) = ddot(davidson%n_parameters, c_i, 1, rho_j, 1)
!
            enddo
         enddo
!
      else
!
!        Pad previous A_red
!
         call mem%alloc(A_red_copy, davidson%dim_red - davidson%n_new_trials, davidson%dim_red - davidson%n_new_trials)
         call copy_and_scale(one, davidson%A_red, A_red_copy, (davidson%dim_red - davidson%n_new_trials)**2)
!
         call mem%dealloc(davidson%A_red, davidson%dim_red - davidson%n_new_trials, davidson%dim_red - davidson%n_new_trials)
!
         call mem%alloc(davidson%A_red, davidson%dim_red, davidson%dim_red)
!
         davidson%A_red(1:davidson%dim_red - davidson%n_new_trials, 1:davidson%dim_red - davidson%n_new_trials) = A_red_copy
!
         call mem%dealloc(A_red_copy, davidson%dim_red - davidson%n_new_trials, davidson%dim_red - davidson%n_new_trials)
!
         call davidson%trials%rewind_()
!
         do i = 1, davidson%dim_red
!
            call davidson%trials%read_(c_i, davidson%n_parameters)
!           
            call davidson%transforms%rewind_()
!
            do j = 1, davidson%dim_red
!
               call davidson%transforms%read_(rho_j, davidson%n_parameters)
!
               if (j .le. davidson%dim_red - davidson%n_new_trials) then
!
                  if (i .gt. davidson%dim_red - davidson%n_new_trials) then 
!
                     davidson%A_red(i,j) = ddot(davidson%n_parameters, c_i, 1, rho_j, 1)
!
                  endif
!
               elseif (j .gt. davidson%dim_red - davidson%n_new_trials) then
!
                  davidson%A_red(i,j) = ddot(davidson%n_parameters, c_i, 1, rho_j, 1)
!
               endif
!
            enddo
!
         enddo
!
      endif
!
      call mem%dealloc(c_i, davidson%n_parameters)
      call mem%dealloc(rho_j, davidson%n_parameters)
!
!     Close files for trial vectors and transformed vectors
!
      call davidson%trials%close_()
      call davidson%transforms%close_()
!
   end subroutine construct_reduced_matrix_davidson_tool
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
      real(dp), dimension(:), allocatable :: c_i
!
      integer :: i
!
      call zero_array(X, davidson%n_parameters)
!
      call mem%alloc(c_i, davidson%n_parameters)
!
      call davidson%trials%open_('read', 'rewind')
!  
      do i = 1, davidson%dim_red
! 
         call davidson%trials%read_(c_i, davidson%n_parameters)
!
         call daxpy(davidson%n_parameters, davidson%X_red(i, n), c_i, 1, X, 1)
!
      enddo    
!
      call mem%dealloc(c_i, davidson%n_parameters)
      call davidson%trials%close_()
!
   end subroutine construct_solution_davidson_tool
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
      real(dp), dimension(:), allocatable :: rho_i
!
      integer :: i
!
      call zero_array(AX, davidson%n_parameters)
!
      call davidson%transforms%open_('read', 'rewind')
!
      call mem%alloc(rho_i, davidson%n_parameters)
!   
      do i = 1, davidson%dim_red
! 
         call davidson%transforms%read_(rho_i, davidson%n_parameters)
!
         call daxpy(davidson%n_parameters, davidson%X_red(i, n), rho_i, 1, AX, 1)
!
      enddo    
!
      call mem%dealloc(rho_i, davidson%n_parameters)
!
      call davidson%transforms%close_()
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
      call davidson%preconditioner%open_('write', 'rewind')
      call davidson%preconditioner%write_(preconditioner, davidson%n_parameters)
      call davidson%preconditioner%close_()
!
      davidson%do_precondition = .true.
!
   end subroutine set_preconditioner_davidson_tool
!
!
   subroutine precondition_davidson_tool(davidson, R, alpha)
!!
!!    Precondition 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
!!    Preconditions the vector R. 
!!
!!    Without the optional α:
!!
!!       R(i) <- R(i)/preconditioner(i).
!!
!!    With α:
!!
!!       R(i) <- R(i)/(preconditioner(i) - alpha).
!!
!!    α will typically be the current (approximated) eigenvalue
!!    if we are solving an eigenvalue problem.
!!
!!    However, if the user has not set any preconditioner, 
!!    this routine performs no action on R. 
!!
      implicit none 
!
      class(davidson_tool), intent(inout) :: davidson
!
      real(dp), dimension(davidson%n_parameters), intent(inout) :: R
!
      real(dp), intent(in), optional :: alpha
!
      real(dp), dimension(:), allocatable :: preconditioner
!
      integer :: i
!
      if (davidson%do_precondition) then 
!
         call mem%alloc(preconditioner, davidson%n_parameters)
!
         call davidson%preconditioner%open_('read', 'rewind')
         call davidson%preconditioner%read_(preconditioner, davidson%n_parameters)
         call davidson%preconditioner%close_()
!
         if (.not. present(alpha)) then
!
!$omp parallel do private(i)
            do i = 1, davidson%n_parameters
!
               R(i) = R(i)/preconditioner(i)
!
            enddo
!$omp end parallel do
!
         else
!
!$omp parallel do private(i)
            do i = 1, davidson%n_parameters
!
               R(i) = R(i)/(preconditioner(i)-alpha)
!
            enddo 
!$omp end parallel do
!
         endif
!
         call mem%dealloc(preconditioner, davidson%n_parameters)
!
      endif 
!
   end subroutine precondition_davidson_tool
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
      call davidson%write_trial(X(:,1), 'rewind')
!
      do solution = 2, davidson%n_solutions
!
         call davidson%write_trial(X(:,solution))
!
      enddo
!
      call mem%dealloc(X, davidson%n_parameters, davidson%n_solutions)
!
!     Delete transformed vectors file, if it is there
!  
      call davidson%transforms%delete_()
!
!     Delete A_red if it is there
!
      if (allocated(davidson%A_red)) call mem%dealloc(davidson%A_red, &
         davidson%dim_red - davidson%n_new_trials, davidson%dim_red - davidson%n_new_trials)
!
      davidson%dim_red = davidson%n_solutions
      davidson%n_new_trials = davidson%n_solutions
!
      call davidson%orthonormalize_trial_vecs()
!
   end subroutine set_trials_to_solutions_davidson_tool
!
!
end module davidson_tool_class
