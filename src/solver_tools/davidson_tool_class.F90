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
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2018
!!
!
   use kinds
   use parameters
!
   use file_class
   use disk_manager_class
   use memory_manager_class
   use array_utilities
!
!
   type, abstract :: davidson_tool
!
      character(len=40) :: name 
!
      real(dp), dimension(:,:), allocatable :: A_red
      real(dp), dimension(:,:), allocatable :: X_red
!
      type(file) :: X, trials, transforms, preconditioner, projector
!
      integer :: dim_red
      integer :: max_dim_red
!
      integer :: n_parameters
      integer :: n_solutions
      integer :: n_new_trials
!
      real(dp) :: residual_threshold
!
      logical :: do_precondition
      logical :: do_projection
!
      integer :: current_n_trials
!
   contains
!
!     Read and write routines
!
      procedure, non_overridable :: read_trial    => read_trial_davidson_tool
      procedure, non_overridable :: write_trial   => write_trial_davidson_tool
!
      procedure, non_overridable :: read_transform  => read_transform_davidson_tool
      procedure, non_overridable :: write_transform => write_transform_davidson_tool  
!
      procedure :: read_solution => read_solution_davidson_tool
!
!     Other procedures
!
      procedure, non_overridable :: construct_reduced_matrix => construct_reduced_matrix_davidson_tool
      procedure, non_overridable :: construct_X              => construct_X_davidson_tool
      procedure, non_overridable :: construct_AX             => construct_AX_davidson_tool
!
      procedure :: set_preconditioner               => set_preconditioner_davidson_tool
      procedure :: set_projector                    => set_projector_davidson_tool
      procedure :: precondition                     => precondition_davidson_tool
      procedure :: projection                       => projection_davidson_tool
      procedure :: orthogonalize_against_trial_vecs => orthogonalize_against_trial_vecs_davidson_tool
      procedure :: orthonormalize_trial_vecs        => orthonormalize_trial_vecs_davidson_tool
!
      procedure :: set_A_red => set_A_red_davidson_tool
      procedure :: get_A_red => get_A_red_davidson_tool
!
      procedure :: set_trials_to_solutions => set_trials_to_solutions_davidson_tool
!
!     Deferred routines  
!
      procedure(solve_reduced_problem), deferred :: solve_reduced_problem
!
      procedure :: read_max_dim_red => read_max_dim_red_davidson_tool
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
   subroutine read_trial_davidson_tool(davidson, c_i, n)
!!
!!    Read trial vector
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2018
!!
!!    Read n'th trial vector from file
!!
      implicit none
!
      class(davidson_tool) :: davidson
!
      real(dp), dimension(davidson%n_parameters, 1) :: c_i
!
      integer, optional :: n
!
      integer :: ioerror
!
      call disk%open_file(davidson%trials, 'read')
!
      if (present(n)) call davidson%trials%prepare_to_read_line(n)
!
      read(davidson%trials%unit, iostat = ioerror) c_i
      if (ioerror .ne. 0) call output%error_msg('reading trial vectors file.')
!
      call disk%close_file(davidson%trials)
!
   end subroutine read_trial_davidson_tool
!
!
   subroutine write_trial_davidson_tool(davidson, c_i, position)
!!
!!    Write trial
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2018
!!
!!    Write trial to file
!!
!!    Optional argument position must be either 'rewind' or 'append'
!!
      implicit none
!
      class(davidson_tool) :: davidson
!
      real(dp), dimension(davidson%n_parameters, 1) :: c_i
!
      character(len=*), optional :: position
!
      integer :: ioerror 
!
!     Was position passed ?
!
      if (present(position)) then
!
!        Sanity check on position variable
!
         if ((trim(position) .ne. 'rewind') .and. (trim(position) .ne. 'append')) then
!
            call output%error_msg('Position specifier not recognized!')
!
         endif
!
         call disk%open_file(davidson%trials, 'write', position)
!
      else
!
         call disk%open_file(davidson%trials, 'write', 'append')
!
      endif 
!
      write(davidson%trials%unit, iostat = ioerror) c_i
      if (ioerror .ne. 0) call output%error_msg('writing trial vectors file.')
!
      call disk%close_file(davidson%trials)
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
      class(davidson_tool) :: davidson 
!
      real(dp), dimension(:,:), allocatable :: c
      real(dp), dimension(:), allocatable :: c_tmp
      real(dp) :: norm_c, ddot
      integer :: i,j,n_done
!
      write(output%unit, '(/t3,a,i0,a)') 'Orthonormalizing the trial vectors.'
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
            c(:,i) = c(:,i) - ddot(davidson%n_parameters, c_tmp, 1, c(1,j), 1)*c(:,j)
!
         enddo
!
         norm_c = sqrt(ddot(davidson%n_parameters, c(1,j), 1, c(1,j), 1))
         c(:,i) = c(:,i)/norm_c 
!
         n_done = n_done + 1
!
      enddo
!
      call disk%open_file(davidson%trials, 'readwrite')
      rewind(davidson%trials%unit)
      do i = 1, davidson%dim_red
!
         write(davidson%trials%unit) c(:,i)
!
      enddo 
!
      call disk%close_file(davidson%trials)
!
      call mem%dealloc(c, davidson%n_parameters, davidson%dim_red)
      call mem%dealloc(c_tmp, davidson%n_parameters)
!
   end subroutine orthonormalize_trial_vecs_davidson_tool
!
!
   subroutine read_solution_davidson_tool(davidson, solution, n)
!!
!!    Read solution 
!!    Written by Eirik F. Kjønstad, Dec 2018
!!
      implicit none 
!
      class(davidson_tool) :: davidson 
!
      real(dp), dimension(davidson%n_parameters, 1) :: solution 
!
      integer :: n
!
      integer :: ioerror
!
      call disk%open_file(davidson%X, 'read')
      call davidson%X%prepare_to_read_line(n)
!
      ioerror = 0
      read(davidson%X%unit, iostat=ioerror) solution 
      if (ioerror .ne. 0) call output%error_msg('could not read davidson solution.')
!
      call disk%close_file(davidson%X)
!
   end subroutine read_solution_davidson_tool
!
!
   subroutine read_transform_davidson_tool(davidson, rho_i, n)
!!
!!    Read transformed vector
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2018
!!
!!    Read n'th transformed vector from file
!!
      implicit none
!
      class(davidson_tool) :: davidson
!
      real(dp), dimension(davidson%n_parameters, 1) :: rho_i
!
      integer :: n
!
      integer :: ioerror
!
      call disk%open_file(davidson%transforms, 'read')
!
      call davidson%transforms%prepare_to_read_line(n)
!
      read(davidson%transforms%unit, iostat = ioerror) rho_i
      if (ioerror .ne. 0) call output%error_msg('reading transformed vectors file 2.')
!
      call disk%close_file(davidson%transforms)
!
   end subroutine read_transform_davidson_tool
!
!
   subroutine write_transform_davidson_tool(davidson, rho_i, position)
!!
!!    Write transform
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2018
!!
!!    Write transformed vector to file
!!
!!    Optional argument position must be either 'rewind' or 'append'
!!
      implicit none
!
      class(davidson_tool) :: davidson
!
      real(dp), dimension(davidson%n_parameters, 1) :: rho_i
!
      character(len=*), optional, intent(in) :: position 
!
      integer :: ioerror   
!
!     Was position passed ?
!
      if (present(position)) then
!
!        Sanity check on position variable
!
         if ((trim(position) .ne. 'rewind') .and. (trim(position) .ne. 'append')) then
!
            call output%error_msg('Position specifier not recognized!')
!
         endif
!
         call disk%open_file(davidson%transforms, 'write', position)
!
      else
!
         call disk%open_file(davidson%transforms, 'write', 'append')
!
      endif 
!
      write(davidson%transforms%unit, iostat = ioerror) rho_i
      if (ioerror .ne. 0) call output%error_msg('writing transformed vector to file.')
!
      call disk%close_file(davidson%transforms)
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
      class(davidson_tool) :: davidson
!
      logical, optional, intent(in) :: entire ! Construct the entire matrix 
!
      real(dp) :: ddot
!
      logical :: entire_local
!
      real(dp), dimension(:,:), allocatable :: A_red_copy, c_i, rho_j
!
      integer :: i, j, ioerror
!
      if (present(entire)) then
         entire_local = entire
      else
         entire_local = .false.
      endif
!
!     Ask disk manager to open the files
!
      call disk%open_file(davidson%trials, 'read')
      call disk%open_file(davidson%transforms, 'read')
!
!     Rewind files
!
      rewind(davidson%trials%unit)
      rewind(davidson%transforms%unit)
!
!     Allocate c and rho
!
      call mem%alloc(c_i, davidson%n_parameters, 1)
      call mem%alloc(rho_j, davidson%n_parameters, 1)
!
!     Construct reduced matrix: A_red_ij = c_i^T * A * c_j = c_i^T * rho_i
!
      if (davidson%dim_red .eq. davidson%n_solutions .or. entire_local) then ! First iteration
!
!        Make the entire reduced matrix
!
         call mem%alloc(davidson%A_red, davidson%dim_red, davidson%dim_red)
!
         do i = 1, davidson%dim_red
!
            read(davidson%trials%unit, iostat=ioerror) c_i
            if (ioerror .ne. 0) call output%error_msg('reading trial vector file.')
!
            rewind(davidson%transforms%unit)
            do j = 1, davidson%dim_red
!
               read(davidson%transforms%unit, iostat=ioerror) rho_j
               if (ioerror .ne. 0) call output%error_msg('reading transformed vector file.')
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
         A_red_copy = davidson%A_red
!
         call mem%dealloc(davidson%A_red, davidson%dim_red - davidson%n_new_trials, davidson%dim_red - davidson%n_new_trials)
!
         call mem%alloc(davidson%A_red, davidson%dim_red, davidson%dim_red)
!
         davidson%A_red(1:davidson%dim_red - davidson%n_new_trials, 1:davidson%dim_red - davidson%n_new_trials) = A_red_copy
!
         call mem%dealloc(A_red_copy, davidson%dim_red - davidson%n_new_trials, davidson%dim_red - davidson%n_new_trials)
!
         rewind(davidson%trials%unit)
!
         do i = 1, davidson%dim_red
!
            read(davidson%trials%unit, iostat=ioerror) c_i
            if (ioerror .ne. 0) call output%error_msg('reading trial vector file.')
!           
            rewind(davidson%transforms%unit)
!
            do j = 1, davidson%dim_red
!
               read(davidson%transforms%unit, iostat=ioerror) rho_j
               if (ioerror .ne. 0) call output%error_msg('reading transformed vector file 3.')
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
      call mem%dealloc(c_i, davidson%n_parameters, 1)
      call mem%dealloc(rho_j, davidson%n_parameters, 1)
!
!     Close files for trial vectors and transformed vectors
!
      call disk%close_file(davidson%trials)
      call disk%close_file(davidson%transforms)
!
   end subroutine construct_reduced_matrix_davidson_tool
!
!
   subroutine construct_X_davidson_tool(davidson, X, n)
!!
!!    Construct X
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2018
!!
!!       X_n = sum_i c_i (X_red_n)_i
!!
      implicit none
!
      class(davidson_tool) :: davidson
!
      integer :: n
!
      real(dp), dimension(davidson%n_parameters, 1) :: X
!
      real(dp), dimension(:,:), allocatable :: c_i
!
      integer :: i, ioerror
!
      X = zero
!
      call mem%alloc(c_i, davidson%n_parameters, 1)
!
      call disk%open_file(davidson%trials, 'read')
      rewind(davidson%trials%unit)
!   
      do i = 1, davidson%dim_red
! 
         read(davidson%trials%unit, iostat = ioerror) c_i
         if (ioerror .ne. 0) call output%error_msg('reading trial vector file.')
!
         call daxpy(davidson%n_parameters, davidson%X_red(i, n), c_i, 1, X, 1)
!
      enddo    
!
      call mem%dealloc(c_i, davidson%n_parameters, 1)
      call disk%close_file(davidson%trials)
!
   end subroutine construct_X_davidson_tool
!
!
   subroutine construct_AX_davidson_tool(davidson, AX, n)
!!
!!    Construct X
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2018
!!
!!       AX_n = sum_i rho_i (X_red_n)_i
!!
      implicit none
!
      class(davidson_tool) :: davidson
!
      integer :: n
!
      real(dp), dimension(davidson%n_parameters, 1) :: AX
!
      real(dp), dimension(:,:), allocatable :: rho_i
!
      integer :: i, ioerror
!
      AX = zero
!
      call mem%alloc(rho_i, davidson%n_parameters, 1)
!
      call disk%open_file(davidson%transforms, 'read')
      rewind(davidson%transforms%unit)
!   
      do i = 1, davidson%dim_red
! 
         read(davidson%transforms%unit, iostat = ioerror) rho_i
         if (ioerror .ne. 0) call output%error_msg('reading transformed vector file.')
!
         call daxpy(davidson%n_parameters, davidson%X_red(i, n), rho_i, 1, AX, 1)
!
      enddo    
!
      call mem%dealloc(rho_i, davidson%n_parameters, 1)
!
      call disk%close_file(davidson%transforms)
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
      class(davidson_tool) :: davidson
!
      real(dp), dimension(davidson%n_parameters, 1), intent(in) :: preconditioner 
!
      call disk%open_file(davidson%preconditioner, 'write', 'rewind')
      write(davidson%preconditioner%unit) preconditioner 
      call disk%close_file(davidson%preconditioner)
!
      davidson%do_precondition = .true.
!
   end subroutine set_preconditioner_davidson_tool
!
!
   subroutine set_projector_davidson_tool(davidson, projector)
!!
!!    Set projector 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
!!    This routine saves the projector to file. 
!!
      implicit none 
!
      class(davidson_tool) :: davidson
!
      real(dp), dimension(davidson%n_parameters, 1), intent(in) :: projector 
!
      call disk%open_file(davidson%projector, 'write', 'rewind')
      write(davidson%projector%unit) projector 
      call disk%close_file(davidson%projector)
!
      davidson%do_projection = .true.
!
   end subroutine set_projector_davidson_tool
!
!
   subroutine precondition_davidson_tool(davidson, R)
!!
!!    Precondition 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
!!    Preconditions the vector R, i.e. 
!!
!!       R(i) <- R(i)/preconditioner(i)
!!
!!    However, if the user has not set any preconditioner, 
!!    this routine performs no action on R. 
!!
      implicit none 
!
      class(davidson_tool) :: davidson
!
      real(dp), dimension(davidson%n_parameters, 1), intent(inout) :: R
!
      real(dp), dimension(:,:), allocatable :: preconditioner
!
      integer :: i 
!
!
      if (davidson%do_precondition) then 
!
         call mem%alloc(preconditioner, davidson%n_parameters, 1)
!
         call disk%open_file(davidson%preconditioner, 'read')
         rewind(davidson%preconditioner%unit)
         read(davidson%preconditioner%unit) preconditioner
         call disk%close_file(davidson%preconditioner)
!
         do i = 1, davidson%n_parameters
!
            R(i, 1) = R(i, 1)/preconditioner(i, 1)
!
         enddo 
!
         call mem%dealloc(preconditioner, davidson%n_parameters, 1)
!
      endif 
!
   end subroutine precondition_davidson_tool
!
!
   subroutine projection_davidson_tool(davidson, R)
!!
!!    Projection
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
!!    Project the vector R, i.e. 
!!
!!       R(i) <- R(i)*projector(i)
!!
!!    However, if the user has not set any preconditioner, 
!!    this routine performs no action on R. 
!!
      implicit none 
!
      class(davidson_tool) :: davidson
!
      real(dp), dimension(davidson%n_parameters, 1), intent(inout) :: R
!
      real(dp), dimension(:,:), allocatable :: projector
!
      integer :: i 
!
      if (davidson%do_projection) then 
!
         call mem%alloc(projector, davidson%n_parameters, 1)
!
         call disk%open_file(davidson%projector, 'read')
         rewind(davidson%projector%unit)
         read(davidson%projector%unit) projector
         call disk%close_file(davidson%projector)
!
         do i = 1, davidson%n_parameters
!
            R(i, 1) = R(i, 1)*projector(i, 1)
!
         enddo 
!
         call mem%dealloc(projector, davidson%n_parameters, 1)
!
      endif 
!
   end subroutine projection_davidson_tool
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
      class(davidson_tool) :: davidson 
!
      real(dp), dimension(davidson%n_parameters, 1) :: R 
!
      real(dp) :: ddot, projection_of_R_on_c_i
!
      integer :: i 
!
      real(dp), dimension(:,:), allocatable :: c_i
!
      call mem%alloc(c_i, davidson%n_parameters, 1)
!
      do i = 1, davidson%dim_red + davidson%n_new_trials
!
         call davidson%read_trial(c_i, i)
         projection_of_R_on_c_i = ddot(davidson%n_parameters, c_i, 1, R, 1)
!
         call daxpy(davidson%n_parameters, - projection_of_R_on_c_i, c_i, 1, R, 1)
!
      enddo 
!
      call mem%dealloc(c_i, davidson%n_parameters, 1)
!
   end subroutine orthogonalize_against_trial_vecs_davidson_tool
!
!
   subroutine set_A_red_davidson_tool(davidson, A_red)
!!
!!    Set A reduced 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
!!    Sets the reduced coefficient matrix A_red equal to the
!!    array sent to the routine. This is a non-standard routine 
!!    that should only be used if the matrix e.g. depends on 
!!    a parameter (scaling certain elements, for instance),
!!    that must be handled outside the Davidson object. 
!!
!!    Usually, construction of A_red is done solely by 
!!    the Davidson object, based on trial vectors and 
!!    transformed vectors on file. 
!!
      implicit none 
!
      class(davidson_tool) :: davidson 
!
      real(dp), dimension(davidson%dim_red, davidson%dim_red), intent(in) :: A_red 
!
      davidson%A_red = A_red 
!
   end subroutine set_A_red_davidson_tool
!
!
   subroutine get_A_red_davidson_tool(davidson, A_red)
!!
!!    Get A reduced 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(davidson_tool) :: davidson 
!
      real(dp), dimension(davidson%dim_red, davidson%dim_red), intent(inout) :: A_red 
!
      A_red = davidson%A_red 
!
   end subroutine get_A_red_davidson_tool
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
      class(davidson_tool) :: davidson 
!
      real(dp), dimension(:,:), allocatable :: X, c_i
!
      integer :: solution, i
!
      real(dp) :: projection_of_X_on_c_i, ddot, norm
!
      call disk%open_file(davidson%trials, 'readwrite')
      call disk%open_file(davidson%X, 'read')
!
      rewind(davidson%X%unit)
      rewind(davidson%trials%unit)
!
      call mem%alloc(X, davidson%n_parameters, 1)
      call mem%alloc(c_i, davidson%n_parameters, 1)
!
      do solution = 1, davidson%n_solutions
!
         read(davidson%X%unit) X
!
         rewind(davidson%trials%unit)
!
         do i = 1, solution - 1
!
            read(davidson%trials%unit) c_i
!
            projection_of_X_on_c_i = ddot(davidson%n_parameters, c_i, 1, X, 1)
!
            call daxpy(davidson%n_parameters, - projection_of_X_on_c_i, c_i, 1, X, 1)
!
         enddo
!
         norm = sqrt(ddot(davidson%n_parameters, X, 1, X, 1))
         call dscal(davidson%n_parameters, one/norm, X, 1)
!
         write(davidson%trials%unit) X
!
      enddo
!
      call mem%dealloc(X, davidson%n_parameters, 1)
      call mem%dealloc(c_i, davidson%n_parameters, 1)
!
      rewind(davidson%trials%unit)
!
      call disk%close_file(davidson%trials)
      call disk%close_file(davidson%X)
!
!     Delete transformed vectors file, if it is there
!
      call disk%delete(davidson%transforms)
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
   subroutine read_max_dim_red_davidson_tool(davidson)
!!
!!    Read max reduced dimension
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none
!
      class(davidson_tool) :: davidson 
!
      character(len=100) :: line
!
      rewind(input%unit)
!
      do 
!
         read(input%unit, '(a100)') line
         line = remove_preceding_blanks(line)
!
         if (trim(line) .eq. 'end geometry') then
            backspace(input%unit)
            return
         endif
!
         if (line(1:22) == 'max reduced dimension:') then
! 
            read(line(23:100), *) davidson%max_dim_red
            return
! 
         endif
!
      enddo
!
   end subroutine read_max_dim_red_davidson_tool
!
!
end module davidson_tool_class
