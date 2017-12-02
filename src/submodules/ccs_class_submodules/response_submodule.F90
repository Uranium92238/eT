submodule (ccs_class) response
!
!!
!!    Response submodule (CCS)
!!    Written by Eirik F. Kjønstad and Sarai Dery Folkestad, June 2017
!!
!!    Contains the following family of procedures of the CCS class: 
!!
!!    response_driver:                       directs the solution of molecular properties.
!!    response_solver:                       solves for a particular molecular property.
!!    initialize_response:                   finds a suitable start trial vector & requests the construction 
!!                                           of the so-called response gradient vector (F).
!!    solve_reduced_response_equation:       solves the response equation in the reduced space
!!                                           of trial vectors (in the given iteration).
!!    construct_gradient_vector:             constructs the gradient vector (F) & saves it to disk.
!!    construct_next_response_trial_vectors: finds the next trial vectors from the residual resulting
!!                                           from the reduced space solution. 
!!    construct_reduced_matrix:              constructs the reduced matrix (the reduced Jacobian, or Jacobian^T).
!!    construct_reduced_gradient:            constructs the reduced gradient vector.
!!    
!!    Note: in order to implement new properties, changes needs to be made in the "transform_trial_vectors"
!!    and "construct_gradient_vector" routines. The behavior of these routines are governed by the value of
!!    the "response_task" string, which must be set in the response driver. 
!!
!
   implicit none 
!
!  Some variables available to all routines of the module
!
   integer(i15) :: iteration = 1
   integer(i15) :: max_iterations = 75
!
!  Variables to handle convergence criterea
!
   logical :: converged = .false. ! True if the residual has converged 
!
contains 
!
!
   module subroutine response_driver_ccs(wf)
!!
!!    Response driver (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
!!    Directs the solution of molculear properties for CCS. The
!!    routine is inherited is to be inherited unaltered in the CC hierarchy. 
!!
      implicit none 
!
      class(ccs) :: wf 
!
!     Let the user know the response driver is running
!
      write(unit_output,'(t3,a)')  ':: Response solver (Davidson)'
      write(unit_output,'(t3,a)')  ':: E. F. Kjønstad, S. D. Folkestad, June 2017' 
      flush(unit_output)    

!     Run the general solver routine (file names are given
!     by the task, i.e., the file 'right_eigenvectors' contains
!     the right eigenvectors)
!
      call wf%response_solver
!
   end subroutine response_driver_ccs
!
!
   module subroutine response_preparations_ccs(wf)
!!
!!
!!
      implicit none
!
      class(ccs) :: wf 
!
      if (wf%tasks%multipliers) then
!
         wf%tasks%current = 'multipliers'
!
         wf%excited_state_specifications%right = .false.
         wf%excited_state_specifications%left  = .true.
!
         if (wf%tasks%core_excited_state) then
!
            wf%excited_state_specifications%solution_file = 'left_core'
!
         else
!
            wf%excited_state_specifications%solution_file = 'left_valence'
!
         endif
!
      endif
!
   end subroutine response_preparations_ccs
!
!
   module subroutine response_solver_ccs(wf)
!!
!!    Response Solver (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
!!    Solves the response equation
!!
!!       A X = F,
!!
!!    where F is the gradient vector and A is either the Jacobian (A) or 
!!    the transposed Jacobian (A^T) matrix. 
!!
!!    The equation is solved by constructing a set of trial vectors, c_i, and 
!!    solving the projected/reduced equations associated with the reduced A and
!!    reduced F vectors:
!!
!!       A_ij = c_i^T A c_j,    F_i = c_i^T F 
!!
!!    The equation is precondition with the orbital energy differences, and the space 
!!    of trial vectors is expanded according to a Davidson algorithm.
!!
      implicit none 
!
      class(ccs) :: wf
!
      integer(i15) :: reduced_dim  = 0  ! Current dimension of the reduced space (i.e., the number of trial vectors)
      integer(i15) :: n_new_trials = 0 ! Number of new trial vectors for the current iteration
!
!     Solution for the reduced eigenvalue problem
!
      real(dp), dimension(:,:), allocatable :: solution_vector_reduced
!
!     Initialize variables 
!
      reduced_dim  = 1
      n_new_trials = 1
!
!     Initialization: - construct the gradient vector and save to file
!                     - use the gradient vector and the orbital differences 
!                       to get the first trial vector 
!
      call wf%initialize_response
!
      write(unit_output,'(/t3,a)')   'Iteration      Residual norm'
      write(unit_output,'(t3,a)')    '----------------------------' 
      flush(unit_output)
!
      do while (.not. converged .and. iteration .le. max_iterations)
!
!        Allocate solution vector 
!
         call wf%mem%alloc(solution_vector_reduced, reduced_dim, 1)
         solution_vector_reduced = zero 
!
!        Transform new trial vector (rho_i = A c_i for new trials i)
!
         call wf%transform_trial_vectors(reduced_dim - n_new_trials + 1, reduced_dim)
!
!        Solve the reduced linear equation (A X = F)
!
         call wf%solve_reduced_response_equation(solution_vector_reduced, reduced_dim, n_new_trials)
!
!        Construct next trial vectors (& test for convergence)
!
         call wf%construct_next_response_trial_vectors(solution_vector_reduced, reduced_dim, n_new_trials)
!
         call wf%mem%dealloc(solution_vector_reduced, reduced_dim, 1)
!
         if (converged) then 
!
            write(unit_output,'(/t3,a,i2,a/)') 'Converged in ',iteration,' iterations!'
!
         else
!
            iteration = iteration + 1
!
         endif
!
      enddo
!
   end subroutine response_solver_ccs
!
!
   module subroutine initialize_response_ccs(wf)
!!
!!    Initialize response (CCS)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, June 2017
!!
!!    Performs two tasks:
!!
!!    1. Initializes start trial vector for response solver. We use the 
!!    the diagonal approximation D of A (and A^T) to form the trial vector:
!!
!!      A X = F => X ~ D^-1 F 
!!
!!    The diagonal approximation of A consists of the orbital energy differences.
!!
!!    2. Constructs the vector F and saves it to file.
!!
      implicit none 
!
      class(ccs) :: wf 
!
      real(dp), dimension(:,:), allocatable :: orbital_diff
      real(dp), dimension(:,:), allocatable :: gradient_vector ! F 
!
      integer :: I = 0
!
      integer(i15) :: ioerror = 0
      integer(i15) :: unit_resp_trial_vecs = 0
      integer(i15) :: unit_grad_vec = 0 
!
!     Open response trial vectors file    
!
      call generate_unit_identifier(unit_resp_trial_vecs)
      open(unit=unit_resp_trial_vecs, file='trial_vec', action='write', status='unknown', &
           access='direct', form='unformatted', recl=dp*(wf%n_parameters), iostat=ioerror)  
!
!     Construct the gradient vector 
!    
      call wf%construct_gradient_vector
!
!     Read gradient vector 
!
      call generate_unit_identifier(unit_grad_vec)
      open(unit=unit_grad_vec, file='gradient_vec', action='read', status='unknown', &
           access='direct', form='unformatted', recl=dp*(wf%n_parameters), iostat=ioerror) 
!
      call wf%mem%alloc(gradient_vector, wf%n_parameters, 1)
      gradient_vector = zero 
!
      read(unit_grad_vec, rec=1, iostat=ioerror) gradient_vector
      if (ioerror .ne. 0) write(unit_output,*) 'WARNING!!'
!
      close(unit_grad_vec)
!
!     Calculate the orbital differences
!
      call wf%mem%alloc(orbital_diff, wf%n_parameters, 1)
      orbital_diff = zero 
!
      call wf%calculate_orbital_differences(orbital_diff)
!
!     Transform by D^-1 - i.e., scale by orbital differences (F <- D^-1 F)
!
      do I = 1, wf%n_parameters
!
         gradient_vector(I, 1) = gradient_vector(I, 1)/orbital_diff(I, 1)
!
      enddo
!
!     Write the trial vector to file 
!
      write(unit_resp_trial_vecs, rec=1, iostat=ioerror) & 
                           (gradient_vector(I,1), I = 1, wf%n_parameters)
!
      call wf%mem%dealloc(gradient_vector, wf%n_parameters, 1)
      call wf%mem%dealloc(orbital_diff, wf%n_parameters, 1)
!
!     Close response trial vectors file
!
      close(unit_resp_trial_vecs)
!
   end subroutine initialize_response_ccs
!
!
   module subroutine solve_reduced_response_equation_ccs(wf, solution_vector_reduced, reduced_dim, n_new_trials)
!!
!!    Solve Reduced Response Equation
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, June 2017
!!
!!    Constructs the reduced A matrix and the reduced F vector,
!!    and solves the (reduced space) linear equation A X = F for X. 
!!
      implicit none
!
      class(ccs) :: wf
!
      integer(i15) :: reduced_dim, n_new_trials
! 
      real(dp), dimension(reduced_dim, 1) :: solution_vector_reduced
!
      real(dp), dimension(:,:), allocatable :: A_reduced ! Reduced Jacobi (or Jacobi^T)
      real(dp), dimension(:,:), allocatable :: F_reduced ! Reduced gradient vector 
!
      integer(i15), dimension(:,:), allocatable :: ipiv ! Pivot integers (see dgesv routine)
!
      integer :: info = -1 ! Error integer for dgesv routine (LU factorization)
!
!     Allocate the reduced Jacobi matrix & the reduced gradient vector
!
      call wf%mem%alloc(A_reduced, reduced_dim, reduced_dim)
      call wf%mem%alloc(F_reduced, reduced_dim, 1)
!
      A_reduced = zero 
      F_reduced = zero 
!
!     Construct the reduced Jacobi and reduced gradient vectors 
!
      call wf%construct_reduced_matrix(A_reduced, reduced_dim, n_new_trials)
      call wf%construct_reduced_gradient(F_reduced, reduced_dim, n_new_trials)
!
!     Solve the eigenvalue problem
!
!     Note: on exit, the solution is in the F_reduced vector,
!     provided info = 0 (see LAPACK documentation for more)
!
      call wf%mem%alloc_int(ipiv, reduced_dim, 1)
      ipiv = 0
!
      info = 0
!
      call dgesv(reduced_dim,  &
                  1,           &
                  A_reduced,   &
                  reduced_dim, &
                  ipiv,        &
                  F_reduced,   &
                  reduced_dim, &
                  info)
!
      if (info .ne. 0) then 
!
         write(unit_output,*) 'Error: could not solve reduced response equation.', info
         stop
!
      endif 
!
!     Save the solution in the solution vector 
!
      solution_vector_reduced = F_reduced
!
!     Deallocations the reduced Jacobi matrix & the reduced gradient vector
!
      call wf%mem%dealloc(A_reduced, reduced_dim, reduced_dim)
      call wf%mem%dealloc(F_reduced, reduced_dim, 1)
!
   end subroutine solve_reduced_response_equation_ccs
!
!
   module subroutine construct_gradient_vector_ccs(wf)
!!
!!    Construct Gradient Vector (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
!!    Constructs the gradient vector, given the current value of "response_task",
!!    and stores the vector on disk for use by the solver. 
!!
      implicit none 
!
      class(ccs) :: wf
!
      real(dp), dimension(:,:), allocatable :: gradient_vector ! F 
!
      integer(i15) :: unit_grad_vec = 0
      integer(i15) :: ioerror = 0
!
!     Allocate and construct the gradient vector 
!
      call wf%mem%alloc(gradient_vector, wf%n_parameters, 1)
      gradient_vector = zero 
!
      if (wf%response_task=='multipliers') then 
!
         write(unit_output,'(/t3,a)') 'Requested the solution of the multiplier equation (t-bar).'
         call wf%construct_eta(gradient_vector)
!
      else
!
         write(unit_output,*) 'Error: Requested gradient vector not recognized.'
         stop
!
      endif
!
!     Open gradient vector file 
!
      call generate_unit_identifier(unit_grad_vec)
      open(unit=unit_grad_vec, file='gradient_vec', action='write', status='unknown', &
           access='direct', form='unformatted', recl=dp*(wf%n_parameters), iostat=ioerror)
!
!     Write gradient vector to file and deallocate
!
      write(unit_grad_vec, rec=1, iostat=ioerror) gradient_vector
      call wf%mem%dealloc(gradient_vector, wf%n_parameters, 1)
!
!     Close file 
!
      close(unit_grad_vec)
!
   end subroutine construct_gradient_vector_ccs
!
!
   module subroutine construct_reduced_matrix_ccs(wf, A_reduced, reduced_dim, n_new_trials)
!!
!!    Construct Reduced Matrix (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
!!    Constructs A_ij = c_i^T A c_j by reading from file and constructing the missing elements.
!!
      implicit none
!
      class(ccs) :: wf
!
      integer(i15) :: reduced_dim, n_new_trials
!
      real(dp), dimension(reduced_dim, reduced_dim) :: A_reduced
!
      real(dp), dimension(:,:), allocatable :: c_i   ! Trial vector i 
      real(dp), dimension(:,:), allocatable :: rho_j ! Transformed trial vector j
!
      integer(i15) :: i = 0, j = 0
!
      real(dp) :: ddot
!
      integer(i15) :: unit_reduced_jacobi = 0
      integer(i15) :: unit_trial_vecs = 0
      integer(i15) :: unit_rho = 0
!
      integer(i15) :: ioerror = 0
!
!     Open files with trial vectors, transformed and nontransformed 
! 
      call generate_unit_identifier(unit_trial_vecs)
      open(unit=unit_trial_vecs, file='trial_vec', action='read', status='old', &
           access='direct', form='unformatted', recl=dp*(wf%n_parameters), iostat=ioerror)
!
      call generate_unit_identifier(unit_rho)
      open(unit=unit_rho, file='transformed_vec', action='read', status='old', &
           access='direct', form='unformatted', recl=dp*(wf%n_parameters), iostat=ioerror) 
!
!     If first iteration, create reduced matrix file. Otherwise,
!     read the part of the reduced matrix already stored on file. 
!
      if (iteration .eq. 1) then 
!
         call generate_unit_identifier(unit_reduced_jacobi)
         open(unit=unit_reduced_jacobi, file='reduced_jacobi', action='readwrite', status='unknown',&
               form='unformatted', iostat=ioerror)
!
      else
!
         call generate_unit_identifier(unit_reduced_jacobi)
         open(unit=unit_reduced_jacobi, file='reduced_jacobi', action='readwrite', status='old',&
          form='unformatted', iostat=ioerror)
!
         rewind(unit_reduced_jacobi)
         read(unit_reduced_jacobi) &
               ((A_reduced(i,j), i = 1, reduced_dim-n_new_trials), j = 1, reduced_dim - n_new_trials)     
!
      endif
!
!     Allocate trial (c) and transformed trial (rho) vectors 
!
      call wf%mem%alloc(c_i, wf%n_parameters, 1)
      call wf%mem%alloc(rho_j, wf%n_parameters, 1)
!
      c_i   = zero
      rho_j = zero
!
!     If first iteration, form the first block of the reduced matrix. Otherwise,
!     form the three new square blocks.
!
      if (iteration .eq. 1) then 
!
         do i = 1,reduced_dim
!
           read(unit_trial_vecs, rec=i, iostat=ioerror) c_i
!
            do j = 1,reduced_dim
!
               read(unit_rho, rec=j, iostat=ioerror) rho_j
               A_reduced(i,j) = ddot(wf%n_parameters, c_i, 1, rho_j, 1)
!
            enddo
!
         enddo
!
      else ! iteration ≠ 1
!
         do i = 1,reduced_dim
!
            read(unit_trial_vecs, rec=i, iostat=ioerror) c_i
!
            do j = reduced_dim - n_new_trials + 1, reduced_dim
!
               read(unit_rho, rec=j, iostat=ioerror) rho_j
               A_reduced(i,j) = ddot(wf%n_parameters, c_i, 1, rho_j, 1)
!
            enddo
!
         enddo
!
         do j = 1, reduced_dim - n_new_trials
!
            read(unit_rho, rec=j, iostat=ioerror) rho_j
!
            do i = reduced_dim - n_new_trials + 1, reduced_dim
!
               read(unit_trial_vecs, rec=i, iostat=ioerror) c_i
               A_reduced(i,j) = ddot(wf%n_parameters, c_i, 1, rho_j, 1)
!
            enddo
!
         enddo
!
      endif
!
!     Close files for trial vectors and transformed vectors
!
      close(unit_trial_vecs)
      close(unit_rho)
!
!     Write current reduced Jacobi matrix to file
!
      rewind(unit_reduced_jacobi)
      write(unit_reduced_jacobi) ((A_reduced(i,j), i = 1, reduced_dim), j = 1, reduced_dim)
!
!     Close reduced Jacobi file
!
      close(unit_reduced_jacobi)
!
   end subroutine construct_reduced_matrix_ccs
!
!
   module subroutine construct_reduced_gradient_ccs(wf, F_reduced, reduced_dim, n_new_trials)
!!
!!    Construct Reduced Gradient (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
!!    Constructs F_i = c_i^T F by reading from file and constructing the missing elements.
!!
      implicit none
!
      class(ccs) :: wf
!
      integer(i15) :: reduced_dim, n_new_trials
!
      real(dp), dimension(:,:), allocatable :: F         ! Gradient vector 
      real(dp), dimension(reduced_dim, 1)   :: F_reduced ! Reduced gradient vector
!
      real(dp), dimension(:,:), allocatable :: c_i ! Trial vector i 
!
      integer(i15) :: i = 0
!
      integer(i15) :: unit_trial_vecs = 0
      integer(i15) :: unit_reduced_gradient = 0
      integer(i15) :: unit_grad_vec = 0
!
      integer(i15) :: ioerror = 0
!
      real(dp) :: ddot
!
!     Open file with trial vectors
! 
      call generate_unit_identifier(unit_trial_vecs)
      open(unit=unit_trial_vecs, file='trial_vec', action='read', status='old', &
           access='direct', form='unformatted', recl=dp*(wf%n_parameters), iostat=ioerror)
!
!     If first iteration, create reduced gradient file. Otherwise,
!     read the part of the reduced gradient already stored on file. 
!
      if (iteration .eq. 1) then 
!
         call generate_unit_identifier(unit_reduced_gradient)
         open(unit=unit_reduced_gradient, file='reduced_gradient', action='readwrite', status='unknown',&
               form='unformatted', iostat=ioerror)
!
      else
!
         call generate_unit_identifier(unit_reduced_gradient)
         open(unit=unit_reduced_gradient, file='reduced_gradient', action='readwrite', status='old',&
          form='unformatted', iostat=ioerror)
!
         rewind(unit_reduced_gradient)
         read(unit_reduced_gradient) &
               (F_reduced(i,1), i = 1, reduced_dim - n_new_trials)     
!
      endif
!
!     Allocate trial (c) vector
!
      call wf%mem%alloc(c_i, wf%n_parameters, 1)
      c_i = zero
!
!     Allocate and read gradient vector (F)
!
      call wf%mem%alloc(F, wf%n_parameters, 1)
      F = zero
!
!     Read gradient vector 
!
      call generate_unit_identifier(unit_grad_vec)
      open(unit=unit_grad_vec, file='gradient_vec', action='read', status='unknown', &
           access='direct', form='unformatted', recl=dp*(wf%n_parameters), iostat=ioerror) 
!
      read(unit_grad_vec, rec=1, iostat=ioerror) F
!
      close(unit_grad_vec)
!
!     Form the new elements of the reduced gradient vector 
!
      do i = reduced_dim - n_new_trials + 1, reduced_dim
!
         read(unit_trial_vecs, rec=i, iostat=ioerror) c_i
         F_reduced(i, 1) = ddot(wf%n_parameters, c_i, 1, F, 1)
!
      enddo
!
!     Close files for trial vectors 
!
      close(unit_trial_vecs)
!
!     Write current reduced gradient vector to file
!
      rewind(unit_reduced_gradient)
      write(unit_reduced_gradient) (F_reduced(i,1), i = 1, reduced_dim)
!
!     Close reduced gradient file
!
      close(unit_reduced_gradient)
!
   end subroutine construct_reduced_gradient_ccs   
!
!
   module subroutine construct_next_response_trial_vectors_ccs(wf, solution_vector_reduced, reduced_dim, n_new_trials)
!!
!!    Construct Next Response Trial Vectors (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
!!    Constructs the next trial vector by constructing the residual vector
!!    
!!       R = (A*X - F)/|X|,
!!
!!    and orthogonalizing it against the previous trial vectors.
!!
!!    Residual vectors are preconditioned before orthogonalization.
!!    This is done by dividing by the orbital differences.
!!    
!!    If the norm of orthogonal vector is very small 
!!    (i.e. high degree of linear dependence on previous trial vectors)
!!    it is scrapped. If norm sufficiently large, the vector is normalized and
!!    stored in trial_vec file, to be used in the next iteration.
!!
!!    The routine also constructs full space solution vectors and stores them
!!    in file solution_vectors when the residuals are converged.
!! 
      implicit none
!
      class(ccs) :: wf
!
      integer(i15) :: reduced_dim, n_new_trials
!
      real(dp), dimension(reduced_dim, 1) :: solution_vector_reduced ! X_reduced 
!
      integer(i15) :: unit_rho = 0
      integer(i15) :: unit_trial_vecs = 0
      integer(i15) :: unit_solution = 0
      integer(i15) :: unit_grad_vec = 0
!
      integer(i15) :: ioerror = 0
!
      integer(i15) :: trial = 0 ! Loops over the trial vectors
!
      real(dp) :: norm_solution_vector = zero
      real(dp) :: norm_new_trial = zero
      real(dp) :: norm_residual = zero
      real(dp) :: conv_test ! Norm to test for convergence
!
      real(dp), parameter :: linear_dependence_threshold = 1.0D-18
!
      real(dp), dimension(:,:), allocatable :: solution_vector ! X 
      real(dp), dimension(:,:), allocatable :: gradient_vector ! F
      real(dp), dimension(:,:), allocatable :: residual        ! R
      real(dp), dimension(:,:), allocatable :: orbital_diff    ! Orbital differences
!
      real(dp), dimension(:,:), allocatable :: c_i 
      real(dp), dimension(:,:), allocatable :: rho_i
!
      integer(i15) :: i = 0, a = 0, ai = 0
!
      real(dp) :: ddot, dot_prod
!
!     Open trial vector files, transformed and untransformed (rho_i = A c_i), 
!     and the solution vector file 
!
      call generate_unit_identifier(unit_trial_vecs)
      open(unit=unit_trial_vecs, file='trial_vec', action='readwrite', status='unknown', &
      access='direct', form='unformatted', recl=dp*(wf%n_parameters), iostat=ioerror)
!
      call generate_unit_identifier(unit_rho)
      open(unit=unit_rho, file='transformed_vec', action='read', status='unknown', &
      access='direct', form='unformatted', recl=dp*(wf%n_parameters), iostat=ioerror) 
!
      call generate_unit_identifier(unit_solution)
      open(unit=unit_solution, file=wf%response_task, action='write', status='unknown', &
      access='direct', form='unformatted', recl=dp*(wf%n_parameters), iostat=ioerror) 
!
!     Construct full space solution vector X 
!
      call wf%mem%alloc(solution_vector, wf%n_parameters, 1)
      solution_vector = zero 
!
      call wf%mem%alloc(c_i, wf%n_parameters, 1)
!
      do trial = 1, reduced_dim
!
!        X = X + solution_vector_reduced_i * c_i 
! 
         c_i = zero 
         read(unit_trial_vecs, rec=trial, iostat=ioerror) c_i
         call daxpy(wf%n_parameters, solution_vector_reduced(trial, 1), c_i, 1, solution_vector, 1)
!
      enddo
!
      call wf%mem%dealloc(c_i, wf%n_parameters, 1)
!
!     Calculate the norm of the solution vector & deallocate 
!
      norm_solution_vector = sqrt(ddot(wf%n_parameters, solution_vector, 1, solution_vector, 1))
!
!     Construct A X and place into residual R 
!
      call wf%mem%alloc(residual, wf%n_parameters, 1)
      residual = zero 
!
      call wf%mem%alloc(rho_i, wf%n_parameters, 1)
!
      do trial = 1, reduced_dim
!
!        R = R + solution_vector_reduced_i * rho_i 
! 
         rho_i = zero 
         read(unit_rho, rec=trial, iostat=ioerror) rho_i
         call daxpy(wf%n_parameters, solution_vector_reduced(trial, 1), rho_i, 1, residual, 1)
!
      enddo    
!
!     Read the gradient vector from disk 
!
      call wf%mem%alloc(gradient_vector, wf%n_parameters, 1)
      gradient_vector = zero 
!
      call generate_unit_identifier(unit_grad_vec)
      open(unit=unit_grad_vec, file='gradient_vec', action='read', status='unknown', &
           access='direct', form='unformatted', recl=dp*(wf%n_parameters), iostat=ioerror) 
!
      read(unit_grad_vec, rec=1, iostat=ioerror) gradient_vector
!
      close(unit_grad_vec)
!
!     Subtract it from the residual (R = A X - F)
!
      call daxpy(wf%n_parameters, -one, gradient_vector, 1, residual, 1)
!
      call wf%mem%dealloc(gradient_vector, wf%n_parameters, 1)
!
!     Calculate the norm of the residual (|| A X - F || / || X ||)
!     and print to output 
!
      norm_residual = sqrt(ddot(wf%n_parameters, residual, 1, residual, 1))
!
      conv_test = norm_residual/norm_solution_vector
!
      if (conv_test .le. wf%response_specifications%residual_threshold) then 
!
         converged = .true.
!
!        Write solution vector to disk 
!
         write(unit_solution, rec=1, iostat=ioerror) solution_vector
!
      endif
!
!     Print information to output 
!
      write(unit_output,'(t3,i2,13x,e10.4)') iteration, conv_test 
      flush(unit_output)
!
!     Precondition the residual by inverse orbital energy differences
!
      call wf%mem%alloc(orbital_diff, wf%n_parameters, 1)
      orbital_diff = zero
!
      call wf%calculate_orbital_differences(orbital_diff)
!
      do i = 1, wf%n_parameters
!
         residual(i, 1) = residual(i,1)/orbital_diff(i,1)
!
      enddo
!
      call wf%mem%dealloc(orbital_diff, wf%n_parameters, 1)
!
!     Orthogonalize the residual against other trials vectors 
!
      call dscal(wf%n_parameters, one/norm_residual, residual, 1) ! Normalize residual 
!
      call wf%mem%alloc(c_i, wf%n_parameters, 1)
!
!     prod_i (I - c_i*c_i^T)*Res = prod_i (Res - c_i*c_i^T*Res)
!
      do trial = 1, reduced_dim
!
         c_i = zero
         read(unit_trial_vecs, rec=trial, iostat=ioerror) c_i
!
         dot_prod = ddot(wf%n_parameters, c_i, 1, residual, 1)
         call daxpy(wf%n_parameters, -dot_prod, c_i, 1, residual,1)
!
      enddo
!
      call wf%mem%dealloc(c_i, wf%n_parameters, 1)
!
      n_new_trials = 0
!
      norm_new_trial = sqrt(ddot(wf%n_parameters, residual, 1, residual, 1))
!
      if ((norm_new_trial .gt. linear_dependence_threshold) .and. .not. converged) then
!
         n_new_trials = n_new_trials + 1
         call dscal(wf%n_parameters, one/norm_new_trial, residual, 1)
         write(unit_trial_vecs, rec=n_new_trials+reduced_dim, iostat=ioerror) residual
!     
      endif      
!
      if (n_new_trials .eq. 0 .and. .not. converged) then
!
         write(unit_output,*) 'Error: linear dependency resulted in no new trial vectors.'
         stop
!
      endif
!
!     Close all files
!
      close(unit_trial_vecs)
      close(unit_rho)
      close(unit_solution)
!
!     Update dimension of reduced space 
!
      reduced_dim = reduced_dim + n_new_trials
!
      call wf%mem%dealloc(residual, wf%n_parameters, 1)
      call wf%mem%dealloc(solution_vector, wf%n_parameters, 1)
!
   end subroutine construct_next_response_trial_vectors_ccs
!
!
end submodule response