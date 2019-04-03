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
module eigen_davidson_tool_class
!
!!
!!    Eigenvalue davidson tool class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2018
!!
!!    A tool to help solve an eigenvalue equation A X_n = omega_n X_n
!!    for symmetric A using the Davidson algorithm. It is tailored to
!!    be usable also in cases where A cannot be stored, but where the 
!!    transformation X -> A X is implemented. A typical Davidson loop
!!    will use this tool as follows:
!!
!!       To-write.
!!
!
   use kinds
   use parameters
!
   use file_class
   use disk_manager_class
   use memory_manager_class
   use davidson_tool_class
   use array_utilities
   use array_analysis
!
   type, extends(davidson_tool) :: eigen_davidson_tool 
!
      real(dp) :: eigenvalue_threshold
!
      type(file) :: projector
!
      logical :: do_projection
!
      real(dp), dimension(:), allocatable :: omega_re 
      real(dp), dimension(:), allocatable :: omega_im
!
   contains 
!
!     Prepare and cleanup
!
      procedure :: prepare                            => prepare_eigen_davidson_tool 
      procedure :: cleanup                            => cleanup_eigen_davidson_tool
!
!     Get-Set routines
!
      procedure :: get_eigenvalue                     => get_eigenvalue_eigen_davidson_tool
!
!     Eigenvalue Davidson specific routines
!
      procedure :: construct_next_trial_vec           => construct_next_trial_vec_eigen_davidson_tool
      procedure :: solve_reduced_problem              => solve_reduced_problem_eigen_davidson_tool
      procedure :: construct_residual                 => construct_residual_eigen_davidson_tool
!
!     Projector routines
!
      procedure :: set_projector                      => set_projector_eigen_davidson_tool
      procedure :: projection                         => projection_eigen_davidson_tool
!
!     Private helper routines
!
      procedure, private :: construct_re_residual     => construct_re_residual_eigen_davidson_tool
      procedure, private :: construct_im_residual     => construct_im_residual_eigen_davidson_tool
!  
      procedure, private :: initialize_omega_re       => initialize_omega_re_eigen_davidson_tool
      procedure, private :: initialize_omega_im       => initialize_omega_im_eigen_davidson_tool
!  
      procedure, private :: destruct_omega_re         => destruct_omega_re_eigen_davidson_tool
      procedure, private :: destruct_omega_im         => destruct_omega_im_eigen_davidson_tool
!
   end type eigen_davidson_tool
!
contains
!
!
   subroutine prepare_eigen_davidson_tool(davidson, name, n_parameters, n_solutions, &
                                             residual_threshold, eigenvalue_threshold)
!!
!!    Initialize 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(eigen_davidson_tool) :: davidson 
!
      character(len=*), intent(in) :: name
!
      integer, intent(in)      :: n_parameters, n_solutions  
      real(dp), intent(in)     :: residual_threshold, eigenvalue_threshold  
!
      davidson%n_parameters = n_parameters
      davidson%n_solutions  = n_solutions
!
      davidson%residual_threshold   = residual_threshold
      davidson%eigenvalue_threshold = eigenvalue_threshold
!
      davidson%name = trim(name)
!
      call davidson%trials%init(trim(davidson%name) // '_trials', 'sequential', 'unformatted')
      call davidson%transforms%init(trim(davidson%name) // '_transforms', 'sequential', 'unformatted')
      call davidson%preconditioner%init(trim(davidson%name) // '_preconditioner', 'sequential', 'unformatted')
      call davidson%projector%init(trim(davidson%name) // '_projector', 'sequential', 'unformatted')
!
!     For safety, delete old files if they are on disk
!
       call disk%delete(davidson%trials)
       call disk%delete(davidson%transforms)
!
      davidson%do_precondition   = .false.         ! Switches to true if 'set_preconditioner' is called
      davidson%do_projection     = .false.         ! Switches to true if 'set_projection' is called
!
      davidson%dim_red           = n_solutions     ! Initial dimension equal to number of solutions
      davidson%n_new_trials      = n_solutions 
!
      davidson%max_dim_red = min(n_solutions*50, 150)   
!
      davidson%current_n_trials = 0
!
      call davidson%read_max_dim_red()     
!
   end subroutine prepare_eigen_davidson_tool
!
!
   subroutine cleanup_eigen_davidson_tool(davidson)
!!
!!    Finalize 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(eigen_davidson_tool) :: davidson 
!
      call disk%open_file(davidson%trials, 'write', 'rewind')
      call disk%open_file(davidson%projector, 'write', 'rewind')
      call disk%open_file(davidson%transforms, 'write', 'rewind')
      call disk%open_file(davidson%preconditioner, 'write', 'rewind')
!
      call disk%close_file(davidson%trials, 'delete')
      call disk%close_file(davidson%projector, 'delete')
      call disk%close_file(davidson%transforms, 'delete')
      call disk%close_file(davidson%preconditioner, 'delete')
!
   end subroutine cleanup_eigen_davidson_tool
!  
!
   subroutine initialize_omega_im_eigen_davidson_tool(davidson)
!!
!!
      implicit none 
!
      class(eigen_davidson_tool) :: davidson 
!
      if (.not. allocated(davidson%omega_im)) call mem%alloc(davidson%omega_im, davidson%n_solutions)
!
   end subroutine initialize_omega_im_eigen_davidson_tool
!  
!
   subroutine initialize_omega_re_eigen_davidson_tool(davidson)
!!
!!
      implicit none 
!
      class(eigen_davidson_tool) :: davidson 
!
      if (.not. allocated(davidson%omega_re)) call mem%alloc(davidson%omega_re, davidson%n_solutions)
!
   end subroutine initialize_omega_re_eigen_davidson_tool
!  
!
   subroutine destruct_omega_im_eigen_davidson_tool(davidson)
!!
!!
      implicit none 
!
      class(eigen_davidson_tool) :: davidson 
!
      if ( allocated(davidson%omega_im)) call mem%dealloc(davidson%omega_im, davidson%n_solutions)
!
   end subroutine destruct_omega_im_eigen_davidson_tool
!  
!
   subroutine destruct_omega_re_eigen_davidson_tool(davidson)
!!
!!
      implicit none 
!
      class(eigen_davidson_tool) :: davidson 
!
      if ( allocated(davidson%omega_re)) call mem%dealloc(davidson%omega_re, davidson%n_solutions)
!
   end subroutine destruct_omega_re_eigen_davidson_tool
!
!
   real(dp) function get_eigenvalue_eigen_davidson_tool(davidson,n)
!!
!!    Get eigenvalue
!!    Written by Eirik F. Kjønstad, Dec 2018 
!!
!!    Returns the nth eigenvalue.
!!
      implicit none 
!
      class(eigen_davidson_tool), intent(in) :: davidson 
!
      integer, intent(in) :: n 
!
      get_eigenvalue_eigen_davidson_tool = davidson%omega_re(n)
!
   end function get_eigenvalue_eigen_davidson_tool
!
!
   subroutine solve_reduced_problem_eigen_davidson_tool(davidson)
!!
!!    Solve reduced problem 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    For the current reduced matrix A, this routine solves the eigenvalue problem 
!!
!!       A X_n = omega_n X_n,    n = 1, 2, ..., dim_red.
!!
!!    On exit, the eigenvalues are stored in the vectors omega_re and omega_im,
!!    where the real and imaginary parts of eigenvalue n is stored in the nth 
!!    row of these vectors.
!!
      implicit none 
!
      class(eigen_davidson_tool) :: davidson 
!
      real(dp), dimension(:), allocatable :: work
!
      real(dp), dimension(:), allocatable :: omega_re
      real(dp), dimension(:), allocatable :: omega_im
!
      integer, dimension(:), allocatable :: index_list
!
      real(dp), dimension(:,:), allocatable :: X_red
      real(dp), dimension(:,:), allocatable :: A_red ! Safe copy to avoid BLAS overwrite
!
      integer :: info = 0, j = 0, i = 0, worksize
      real(dp)  :: dummy =0.0, optwork
!
!     Solve reduced eigenvalue problem
!
      info = 0
!
      call mem%alloc(X_red, davidson%dim_red, davidson%dim_red)
      call mem%alloc(A_red, davidson%dim_red, davidson%dim_red)
!
      X_red = zero
      call dcopy(davidson%dim_red**2, davidson%A_red, 1, A_red, 1)
!
      call mem%alloc(omega_re, davidson%dim_red)
      call mem%alloc(omega_im, davidson%dim_red)
!
      omega_re = zero
      omega_im = zero
!
!     Find optimal work size, lwork = -1
      call dgeev('N','V',             &
                  davidson%dim_red,   &
                  A_red,              &
                  davidson%dim_red,   &
                  omega_re,           &
                  omega_im,           &
                  dummy,              &              
                  1,                  &
                  X_red,              &
                  davidson%dim_red,   &
                  optwork,            &   
                  -1,                 &
                  info)
!
      worksize = int(optwork+0.01)
!
      call mem%alloc(work, worksize)
      work = zero
!      
      call dgeev('N','V',             &
                  davidson%dim_red,   &
                  A_red,              &
                  davidson%dim_red,   &
                  omega_re,           &
                  omega_im,           &
                  dummy,              &              
                  1,                  &
                  X_red,              &
                  davidson%dim_red,   &
                  work,               &   
                  worksize,           &
                  info)
!
      call mem%dealloc(A_red, davidson%dim_red, davidson%dim_red)
!
      if (info .ne. 0) call output%error_msg('could not solve reduced equation in Davidson davidson.')
!
      call mem%dealloc(work, worksize)
!
!     Find lowest n_solutions eigenvalues and sort them (the corresponding indices
!     are placed in the integer array index_list)
!
      call davidson%initialize_omega_im()
      call davidson%initialize_omega_re()
!
      davidson%omega_re = zero
      davidson%omega_im = zero
!
      call mem%alloc(index_list, davidson%n_solutions)
!
      call get_n_lowest(davidson%n_solutions, davidson%dim_red, &
                        omega_re, davidson%omega_re, index_list)
!
!     Pick out solution vectors and imaginary parts of eigenvalues according to index_list
!
      if (allocated(davidson%X_red)) &
         call mem%dealloc(davidson%X_red, davidson%dim_red - davidson%n_new_trials, davidson%n_solutions)
!
      call mem%alloc(davidson%X_red, davidson%dim_red, davidson%n_solutions)
      davidson%X_red = zero
!
      do j = 1, davidson%n_solutions
!
         do i = 1, davidson%dim_red
!
            davidson%X_red(i, j) = X_red(i, index_list(j))
!
         enddo
!
         davidson%omega_im(j) = omega_im(index_list(j))
!
      enddo
!
      call mem%dealloc(X_red, davidson%dim_red, davidson%dim_red)
      call mem%dealloc(omega_re, davidson%dim_red)
      call mem%dealloc(omega_im, davidson%dim_red)
      call mem%dealloc(index_list, davidson%n_solutions)
!
   end subroutine solve_reduced_problem_eigen_davidson_tool
!
!
   subroutine construct_residual_eigen_davidson_tool(davidson, R, X, norm_X, n)
!!
!!    Construct residual 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    Constructs the nth full space residual 
!!
!!       R = A X - omega X,
!!
!!    without any preconditioning (this is applied after a call
!!    to construct residual from the construct next trial vector
!!    routine). Here, X is the nth eigenvector and omega is the 
!!    eigenvalue. The residual is normalized by the norm of the 
!!    solution X.
!!
      implicit none 
!
      class(eigen_davidson_tool), intent(in) :: davidson 
!
      real(dp), dimension(davidson%n_parameters)                :: R 
      real(dp), dimension(davidson%n_parameters), intent(in) :: X 
!
      integer, intent(in)  :: n
      real(dp), intent(in) :: norm_X 
!
      if (davidson%omega_im(n) .eq. zero) then  ! standard case: the nth root is not part of a complex pair
!
         call davidson%construct_AX(R, n) ! set R = AX 
!
        call daxpy(davidson%n_parameters, - davidson%omega_re(n), X, 1,  R, 1)
        call dscal(davidson%n_parameters, one/norm_X, R, 1)
!
      else
!
!        If it's the first root of the complex pair, construct the real residual; if it's the second, 
!        construct instead the imaginary residual
!
         if (n .eq. 1) then
!
            if (davidson%n_solutions .lt. 2) call output%error_msg('add one more root to treat the complex pair.')
!
            call davidson%construct_re_residual(R, X, norm_X, n)
!
         elseif (n .eq. davidson%n_solutions) then 
!
            if (davidson%omega_re(n-1) .ne. davidson%omega_re(n)) &
                        call output%error_msg('add one more root to treat the complex pair.')
!
            call davidson%construct_im_residual(R, X, norm_X, n)
!
         else ! neither first or last, so it's safe to look at n + 1 and n - 1 
!
            if (davidson%omega_re(n) .eq. davidson%omega_re(n - 1)) then
!
               call davidson%construct_im_residual(R, X, norm_X, n)
!
            elseif (davidson%omega_re(n) .eq. davidson%omega_re(n + 1) ) then 
!
               call davidson%construct_re_residual(R, X, norm_X, n)
!
            else ! should never happen, but just in case, let the user know there's a bug
!
               call output%error_msg('something went very wrong when trying to construct imaginary residual.')
!
            endif 
!
         endif
!
      endif
!
   end subroutine construct_residual_eigen_davidson_tool
!
!
   subroutine construct_re_residual_eigen_davidson_tool(davidson, R, X_re, norm_X_re, n)
!!
!!    Construct real residual 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
!!    For a complex pair, 
!!
!!       X+ = X_re + i X_im 
!!       X- = X_re - i X_im,
!!
!!    we will find the real part X_re as the first of the two roots,
!!    and the imaginary part X_im as the second of the two roots 
!!    (in X_red). This routine constructs the residual associated
!!    with the real part. As such, the construct residual routine 
!!    makes sure the n-th root is X_re and (n + 1)-th is X_im (and
!!    exists, which might not be the case if too few roots are 
!!    requested). On exit, the real residual is placed in R: 
!!    
!!       R = (A X_re - omega_re X_re) + omega_im X_im 
!!
!!    The residual is divided by the norm of X+ (or, 
!!    equivalently, X-).
!!
      implicit none 
!
      class(eigen_davidson_tool), intent(in) :: davidson 
!
      real(dp), dimension(davidson%n_parameters)             :: R 
      real(dp), dimension(davidson%n_parameters), intent(in) :: X_re  
!
      integer, intent(in)  :: n 
      real(dp), intent(in) :: norm_X_re
!
      real(dp), dimension(:), allocatable :: X_im 
!
      real(dp) :: norm_X_im
!
      call mem%alloc(X_im, davidson%n_parameters)  
!
      call davidson%construct_X(X_im, n + 1) ! set X_im
      call davidson%construct_AX(R, n)       ! set R = A X_re 
!
      call daxpy(davidson%n_parameters, - davidson%omega_re(n), X_re, 1, R, 1) 
      call daxpy(davidson%n_parameters, davidson%omega_im(n), X_im, 1, R, 1) 
!
      norm_X_im = get_l2_norm(X_im, davidson%n_parameters)
!
      call dscal(davidson%n_parameters, one/(sqrt(norm_X_re**2 + norm_X_im**2)), R, 1)
!
      call mem%dealloc(X_im, davidson%n_parameters) 
!
   end subroutine construct_re_residual_eigen_davidson_tool
!
!
   subroutine construct_im_residual_eigen_davidson_tool(davidson, R, X_im, norm_X_im, n)
!!
!!    Construct imaginary residual 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
!!    For a complex pair, 
!!
!!       X+ = X_re + i X_im 
!!       X- = X_re - i X_im,
!!
!!    we will find the real part X_re as the first of the two roots,
!!    and the imaginary part X_im as the second of the two roots 
!!    (in X_red). This routine constructs the residual associated
!!    with the imaginary part. As such, the construct residual 
!!    routine makes sure the n-th root is X_im and (n - 1)-th 
!!    is X_re. On exit, the imaginary residual is placed in R: 
!!    
!!       R = (A X_im - omega_re X_im) - omega_im X_re 
!!
!!    The residual is divided by the norm of X+ (or, equivalently, 
!!    X-). Note that omega_im of root n is the negative of omega_im 
!!    of root n - 1. To avoid inconsistency, we use omega_im from 
!!    the previous root (otherwise, we would change the sign of the
!!    last constribution).
!!
      implicit none 
!
      class(eigen_davidson_tool), intent(in) :: davidson 
!
      real(dp), dimension(davidson%n_parameters)             :: R 
      real(dp), dimension(davidson%n_parameters), intent(in) :: X_im  
!
      integer, intent(in)  :: n 
      real(dp), intent(in) :: norm_X_im
!
      real(dp), dimension(:), allocatable :: X_re  
!
      real(dp) :: norm_X_re
!
      call mem%alloc(X_re, davidson%n_parameters)  
!
      call davidson%construct_X(X_re, n - 1) ! set X_re 
      call davidson%construct_AX(R, n)       ! set R = A X_im 
!
      call daxpy(davidson%n_parameters, - davidson%omega_re(n), X_im, 1, R, 1) 
      call daxpy(davidson%n_parameters, - davidson%omega_im(n - 1), X_re, 1, R, 1) 
!
      norm_X_re = get_l2_norm(X_re, davidson%n_parameters)
!
      call dscal(davidson%n_parameters, one/(sqrt(norm_X_re**2 + norm_X_im**2)), R, 1)
!
      call mem%dealloc(X_re, davidson%n_parameters)  
!
   end subroutine construct_im_residual_eigen_davidson_tool
!
!
   subroutine construct_next_trial_vec_eigen_davidson_tool(davidson, residual_norm, n)
!!
!!    Construct next trial vector  
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
!!    This routine decides whether to construct a new trial vector based on 
!!    the current residual norm. The new trial vector is a preconditioned 
!!    residual orthogonalized against the search space (the existing trial
!!    vectors). If it turns out that the root is not converged and the 
!!    new trial vector contains no new components (linear dependence),
!!    the routine terminates in an error. 
!!
!!    Note that the residual norm is set on exit and need not to be 
!!    calculated beforehand. 
!!    
      implicit none 
!
      class(eigen_davidson_tool) :: davidson 
!
      real(dp), intent(out) :: residual_norm 
!
      integer, optional, intent(in) :: n 
!
      integer :: k ! k = n, where k is set to 1 if n is not present 
!
      real(dp) :: norm_X, norm_new_trial, norm_residual, norm_precond_residual
!
      real(dp), dimension(:), allocatable :: R, X 
!
      if (present(n)) then
! 
         k = n
!
      else
!
         k = 1
!
      endif
!
!     Construct full space solution vector X, 
!     and the associated residual R 
!
      call mem%alloc(X, davidson%n_parameters)
      call mem%alloc(R, davidson%n_parameters)
!
      call davidson%construct_X(X, k) 
      norm_X = get_l2_norm(X, davidson%n_parameters) 
!
      call davidson%construct_residual(R, X, norm_X, k)
!
      call dscal(davidson%n_parameters, one/norm_X, X, 1)
!
      call mem%dealloc(X, davidson%n_parameters)
!
!     Calculate the norm of the residual and test for convergence. If not
!     converged, the residual is preconditioned & we remove components already 
!     in the search space by orthogonalizing it against existing trial vectors.
!
      norm_residual = get_l2_norm(R, davidson%n_parameters)
      residual_norm = norm_residual
!
      if (norm_residual .gt. davidson%residual_threshold) then 
!
         call davidson%projection(R)
         call davidson%precondition(R)
!
         norm_precond_residual = get_l2_norm(R, davidson%n_parameters)
         call dscal(davidson%n_parameters, one/norm_precond_residual, R, 1)
!
         call davidson%orthogonalize_against_trial_vecs(R)
         norm_new_trial = get_l2_norm(R, davidson%n_parameters)
!
         if (norm_new_trial .gt. davidson%residual_threshold) then
!
            davidson%n_new_trials = davidson%n_new_trials + 1
            call dscal(davidson%n_parameters, one/norm_new_trial, R, 1)
            call davidson%write_trial(R, 'append')
!
         endif
!
      endif 
!
      call mem%dealloc(R, davidson%n_parameters)
!
   end subroutine construct_next_trial_vec_eigen_davidson_tool
!
!
   subroutine set_projector_eigen_davidson_tool(davidson, projector)
!!
!!    Set projector 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
!!    This routine saves the projector to file. 
!!
      implicit none 
!
      class(eigen_davidson_tool) :: davidson
!
      real(dp), dimension(davidson%n_parameters), intent(in) :: projector 
!
      call disk%open_file(davidson%projector, 'write', 'rewind')
      write(davidson%projector%unit) projector 
      call disk%close_file(davidson%projector)
!
      davidson%do_projection = .true.
!
   end subroutine set_projector_eigen_davidson_tool
!
!
   subroutine projection_eigen_davidson_tool(davidson, R)
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
      class(eigen_davidson_tool) :: davidson
!
      real(dp), dimension(davidson%n_parameters), intent(inout) :: R
!
      real(dp), dimension(:), allocatable :: projector
!
      integer :: i 
!
      if (davidson%do_projection) then 
!
         call mem%alloc(projector, davidson%n_parameters)
!
         call disk%open_file(davidson%projector, 'read')
         rewind(davidson%projector%unit)
         read(davidson%projector%unit) projector
         call disk%close_file(davidson%projector)
!
         do i = 1, davidson%n_parameters
!
            R(i) = R(i)*projector(i)
!
         enddo 
!
         call mem%dealloc(projector, davidson%n_parameters)
!
      endif 
!
   end subroutine projection_eigen_davidson_tool
!
!
end module eigen_davidson_tool_class
