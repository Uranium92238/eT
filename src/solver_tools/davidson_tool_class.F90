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
   use record_storer_class, only: record_storer
   use memory_storer_class, only: memory_storer
   use file_storer_class, only: file_storer
   use memory_manager_class, only: mem 
   use global_out, only: output
   use array_utilities, only: get_l2_norm, copy_and_scale, zero_array
   use precondition_tool_class, only: precondition_tool
!
!
   type, abstract :: davidson_tool
!
      character(len=40) :: name_ 
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
      procedure :: iterate                                     => iterate_davidson_tool 
      procedure :: construct_solution                          => construct_solution_davidson_tool
      procedure :: set_trial                                   => set_trial_davidson_tool 
      procedure :: get_trial                                   => get_trial_davidson_tool 
      procedure :: set_transform                               => set_transform_davidson_tool 
      procedure :: get_transform                               => get_transform_davidson_tool 
      procedure :: first_trial                                 => first_trial_davidson_tool
      procedure :: last_trial                                  => last_trial_davidson_tool
!
      procedure :: set_preconditioner                          => set_preconditioner_davidson_tool
!
      procedure(solve_reduced_problem), deferred :: solve_reduced_problem
!
!     Other routines 
!
      procedure, non_overridable :: construct_AX               => construct_AX_davidson_tool
      procedure, non_overridable :: construct_reduced_matrix   => construct_reduced_matrix_davidson_tool
!
      procedure :: orthonormalize_trial_vecs                   => orthonormalize_trial_vecs_davidson_tool
!
      procedure :: set_trials_to_solutions                     => set_trials_to_solutions_davidson_tool
!
      procedure :: initialize_trials_and_transforms            => initialize_trials_and_transforms_davidson_tool
      procedure :: finalize_trials_and_transforms              => finalize_trials_and_transforms_davidson_tool
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
   subroutine initialize_trials_and_transforms_davidson_tool(davidson, records_in_memory)
!!
!!    Initialize trials and transforms  
!!    Written by Eirik F. Kjønstad, Nov 2019 
!!
!!    Initializes the storers for trials and transforms.
!!
!!    records_in_memory: if true,  trials and transforms are stored in memory 
!!                       if false, trials and transforms are stored on disk 
!!
!!    Trial vectors c define the subspace. Transforms refer to the transformed 
!!    trial vectors, i.e. rho = A c, if A is the linear transformation of 
!!    the linear equation. 
!!
      implicit none 
!
      class(davidson_tool) :: davidson 
!
      logical, intent(in) :: records_in_memory
!
      if (records_in_memory) then 
!
         call output%printf('Reduced space basis and transforms are stored in memory.', &
                                    pl='m', fs='(/t6,a)')
!
         davidson%trials = memory_storer(trim(davidson%name_) // '_trials', &
                     davidson%n_parameters, davidson%max_dim_red + davidson%n_solutions)   
!
         davidson%transforms = memory_storer(trim(davidson%name_) // '_transforms', &
                     davidson%n_parameters, davidson%max_dim_red + davidson%n_solutions) 
!
      else
!
         call output%printf('Reduced space basis and transforms are stored on disk.', &
                                 pl='m', fs='(/t6,a)')
!
         davidson%trials = file_storer(trim(davidson%name_) // '_trials', &
                     davidson%n_parameters, davidson%max_dim_red + davidson%n_solutions, &
                     delete=.true.)   
!
         davidson%transforms = file_storer(trim(davidson%name_) // '_transforms', &
                     davidson%n_parameters, davidson%max_dim_red + davidson%n_solutions, &
                     delete=.true.)   
!
      endif 
!
      call davidson%trials%initialize_storer()
      call davidson%transforms%initialize_storer()
!
   end subroutine initialize_trials_and_transforms_davidson_tool
!
!
   subroutine finalize_trials_and_transforms_davidson_tool(davidson)
!!
!!    Finalize trials and transforms 
!!    Written by Eirik F. Kjønstad, Nov 2019 
!!
!!    Finalizes the storers. This means for files that they are closed 
!!    and deleted.
!!
      implicit none 
!
      class(davidson_tool) :: davidson 
!
      call davidson%trials%finalize_storer()
      call davidson%transforms%finalize_storer()
!
   end subroutine finalize_trials_and_transforms_davidson_tool
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
      call davidson%trials%set(c, n)
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
      call davidson%trials%get(c, n)
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
      call davidson%transforms%set(rho, n)
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
      call davidson%transforms%get(rho, n)
!
   end subroutine get_transform_davidson_tool
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
      if (davidson%dim_red + davidson%n_new_trials >= davidson%max_dim_red) then
!
         call davidson%set_trials_to_solutions()
!
      else 
!
         davidson%dim_red = davidson%dim_red + davidson%n_new_trials
!
      endif 
!
      call davidson%orthonormalize_trial_vecs()
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
      integer, intent(in) :: n
!
      call davidson%trials%get(c, n)
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
      call davidson%trials%set(c, n)
!
   end subroutine write_trial_davidson_tool
!
!
   subroutine orthonormalize_trial_vecs_davidson_tool(davidson)
!!
!!    Orthonormalize trial vecs  
!!    Written by Eirik F. Kjønstad, Oct 2019 
!!
!!    Orthonormalizes the full set of trial vectors using the 
!!    modified Gram-Schmidt procedure, see e.g. Iterative Methods 
!!    for Sparse Linear Systems, 2nd ed., Yousef Saad.
!!
      implicit none 
!
      class(davidson_tool), intent(inout) :: davidson 
!
      real(dp), dimension(:), allocatable :: c_i, c_j
!
      real(dp) :: norm_c_j, r_ji, ddot 
!
      integer :: i, j
!
      call output%printf('Orthonormalizing trial vectors.', pl='v', fs='(/t3,a)')
!
      call mem%alloc(c_i, davidson%n_parameters)
      call mem%alloc(c_j, davidson%n_parameters)
!
!     First orthonormalization of new trials 
!
      do j = davidson%first_trial(), davidson%last_trial() 
!
         call davidson%get_trial(c_j, j)
!
         do i = 1, j - 1
!
            call davidson%get_trial(c_i, i)
!
            r_ji = ddot(davidson%n_parameters, c_j, 1, c_i, 1)
            call daxpy(davidson%n_parameters, -r_ji, c_i, 1, c_j, 1)
!
         enddo
!
         norm_c_j = sqrt(ddot(davidson%n_parameters, c_j, 1, c_j, 1))
!
         if (norm_c_j < davidson%lindep_threshold) then 
!
            call output%error_msg('detected linear dependence in trial space in davidson_tool.')
!
         endif
!
         call dscal(davidson%n_parameters, one/norm_c_j, c_j, 1)
!
         call davidson%set_trial(c_j, j)
!
      enddo
!
!     Second orthonormalization to avoid accumulation of noise
!
      do j = davidson%first_trial(), davidson%last_trial() 
!
         call davidson%get_trial(c_j, j)
!
         do i = 1, j - 1
!
            call davidson%get_trial(c_i, i)
!
            r_ji = ddot(davidson%n_parameters, c_j, 1, c_i, 1)
            call daxpy(davidson%n_parameters, -r_ji, c_i, 1, c_j, 1)
!
         enddo
!
         norm_c_j = sqrt(ddot(davidson%n_parameters, c_j, 1, c_j, 1))
!
         if (norm_c_j < davidson%lindep_threshold) then 
!
             call output%error_msg('detected linear dependence in trial space in davidson_tool.')
!
         endif
!
         call dscal(davidson%n_parameters, one/norm_c_j, c_j, 1)
!
         call davidson%set_trial(c_j, j)
!
      enddo
!
      call mem%dealloc(c_i, davidson%n_parameters)
      call mem%dealloc(c_j, davidson%n_parameters)
!
   end subroutine orthonormalize_trial_vecs_davidson_tool
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
            call davidson%get_trial(c_i, i)
!
            do j = 1, davidson%dim_red
!
               call davidson%get_transform(rho_j, j)
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
         do i = 1, davidson%dim_red
!
            call davidson%get_trial(c_i, i)
!           
            do j = 1, davidson%dim_red
!
               call davidson%get_transform(rho_j, j)
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
      do i = 1, davidson%dim_red
! 
         call davidson%get_trial(c_i, i)
!
         call daxpy(davidson%n_parameters, davidson%X_red(i, n), c_i, 1, X, 1)
!
      enddo    
!
      call mem%dealloc(c_i, davidson%n_parameters)
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
      call mem%alloc(rho_i, davidson%n_parameters)
!   
      do i = 1, davidson%dim_red
! 
         call davidson%get_transform(rho_i, i)
!
         call daxpy(davidson%n_parameters, davidson%X_red(i, n), rho_i, 1, AX, 1)
!
      enddo    
!
      call mem%dealloc(rho_i, davidson%n_parameters)
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
      davidson%preconditioner = precondition_tool(preconditioner, davidson%n_parameters)
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
      do solution = 1, davidson%n_solutions
!
         call davidson%set_trial(X(:,solution), solution)
!
      enddo
!
      call mem%dealloc(X, davidson%n_parameters, davidson%n_solutions)
!
!     Delete A_red if it is there
!
      if (allocated(davidson%A_red)) call mem%dealloc(davidson%A_red, &
         davidson%dim_red - davidson%n_new_trials, davidson%dim_red - davidson%n_new_trials)
!
      davidson%dim_red = davidson%n_solutions
      davidson%n_new_trials = davidson%n_solutions
!
   end subroutine set_trials_to_solutions_davidson_tool
!
!
end module davidson_tool_class
