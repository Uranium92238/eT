submodule(ccs_class) excited_state
!
!!
!!    Excited state sub(CCS)
!!    Written by Eirik F. Kjønstad and Sarai Dery Folkestad, May 2017
!!
!!    Contains the following family of procedures of the CCS class:
!!
!!    excited_state_driver:                directs the solution of excited state problems.
!!    excited_state_solver:                solves the excited state problem.
!!    solve_reduced_eigenvalue_equation:   solves the excited state problem in the projected/reduced space 
!!                                         of trial vectors (in a given iteration).
!!    construct_next_trial_vectors:        finds the new trial vectors resulting from the residuals found 
!!                                         by solving the reduced eigenvalue equation.
!!    calculate_orbital_differences:       calculates the orbital differences (used for preconditioning and
!!                                         start vector guess).
!!    transform_trial_vectors:             transforms (by A) the new trial vectors and saves them to disk.
!!    find_start_trial_indices:            find the indices corresponding to the lowest orbital differences.
!!    trial_vectors_from_stored_solutions: finds suitable start trial vectors from stored solutions (for restart).
!!    
!
   implicit none 
!
!  Some variables available to all routines of the module
!
   integer(i15) :: iteration = 1
   integer(i15) :: max_iterations = 75 ! E: we move this to calculation settings later
!
!  Variables to handle convergence criterea
!
   logical :: converged = .false. ! True iff both the energy and the equations have converged 
!
   logical :: converged_energy   = .false.
   logical :: converged_residual = .false.
!
   logical :: timings = .true.
   logical :: print_vectors = .true.
!
!
contains
!
!
   module subroutine excited_state_driver_ccs(wf)
!!
!!    Excited state driver (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
!!    Directs the solution of the excited state problem for CCS. The
!!    routine is inherited is to be inherited unaltered in the CC hierarchy. 
!!
!!    Note: it is only necessary to alter this routine if the excited states are 
!!    solved for by a different algorithm (such as in similarity constrained CC, 
!!    where the excited states and ground state are determined simultaneously).
!!
      implicit none 
!
      class(ccs) :: wf 
!
!     Let the user know the excited state driver is running
!
      write(unit_output,'(t3,a)')    ':: Excited state solver (Davidson)'
      write(unit_output,'(t3,a/)')   ':: E. F. Kjønstad, S. D. Folkestad, May 2017'
      write(unit_output,'(t3,a,i3,a,a,a)') &
                                     'Requested ',wf%tasks%n_singlet_states,' ', trim(wf%name), ' singlet states.'
      write(unit_output,'(t3,a,i3,a,a,a)') &
                                     'Requested ',wf%tasks%n_triplet_states,' ', trim(wf%name), ' triplet states.'     
!
!     Preparations for excited state solver 
!
      call wf%excited_state_preparations
!
!     Set current task to excited state calculation 
! 
      wf%current_task = 'excited_state'
!
!     Run the general solver routine (file names are given
!     by the task, i.e., the file 'right_valence' contains
!     the right eigenvectors)
!
      call wf%excited_state_solver
!
!     Final work and preparations for other tasks (such as property calculations)
!
      call wf%excited_state_cleanup
!
   end subroutine excited_state_driver_ccs
!
!
   module subroutine excited_state_preparations_ccs(wf)
!!
!!    Excited State Preparations (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
!!    A routine for preparation tasks (if any). Can be overwritten
!!    in descendants if other preparations prove necessary.    
!!
      class(ccs) :: wf 
!
!     Store voov-electronic repulsion integrals to file if there is space
!
      call wf%store_t1_vo_ov_electronic_repulsion
!
   end subroutine excited_state_preparations_ccs
!
!
   module subroutine excited_state_cleanup_ccs(wf)
!!
!!    Excited State Cleanup (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
!!    A routine for cleanup tasks (if any). Can be overwritten
!!    in descendants if other cleanups prove necessary.    
!!
      class(ccs) :: wf 
!
!     Nothing yet!
!
   end subroutine excited_state_cleanup_ccs
!
!
   module subroutine excited_state_solver_ccs(wf)
!!
!!    Excited State Solver
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!    Directs the solution of the excited states using a Davidson algorithm.
!!    The routine aims to find the right eigenvectors of the matrix A:
!!
!!       A X = e X,  
!!
!!    and the eigenvalues which corresponds to the excitation energies. 
!!
!!    The matrix A above can be the Jacobian (the usual A) or the transposed Jacobian (A^T), or,
!!    in principle, any matrix. The transformation used is determined by the value of the wavefunction's 
!!    "response_task" variable (if "right_eigenvectors", use A; if "left_eigenvectors", use A^T; and so on).
!!    The selection of A is done in the routine transform_trial_vectors.
!!
!!    The problem is solved in a reduced space. To find n roots, n start trial vectors {c_i}_i=1, 
!!    n are generated according to the lowest orbital differences. Then a reduced space Jacobian 
!!    is constructed,
!! 
!!       A_red_ij = c_i^T * A c_j,
!!
!!    and the eigenvalues e and eigenvectors x of this matrix are found.
!!    The full space vectors {X_j}_j=1,n are then given by
!!
!!       X_j = sum_i x_j_i*c_i, 
!! 
!!    and the j'th residual vector is given by
!!
!!       R_j = (A*X_j - e*X_j)/|X_j|.
!!
!!    If the norm of this residual is sufficiently small (and the excitation energies 
!!    are converged within a given threshold), convergence is reached. If not, new trial 
!!    vectors will be generated by orthogonalizing the residual vector against the previous 
!!    trial vectors and then normalizing them, thereby expanding the dimension 
!!    of the reduced space for the next iteration.
!!
!!    The linear system (equivalently, the residual) is preconditioned with a diagonal 
!!    matrix with elements equal to the inverse orbital differences.
!!   
      implicit none
!
      class(ccs) :: wf
!
      integer(i15) :: reduced_dim  = 0 ! Current dimension of the reduced space (i.e., the number of trial vectors)
      integer(i15) :: n_new_trials = 0 ! Number of new trial vectors for the current iteration 
!
!     Solution for the reduced eigenvalue problem
!
      real(dp), dimension(:,:), allocatable :: eigenvalues_Re_new ! Eigenvalues of current iteration 
      real(dp), dimension(:,:), allocatable :: eigenvalues_Im_new
!
      real(dp), dimension(:,:), allocatable :: eigenvalues_Re_old ! Eigenvalues of previous iteration 
      real(dp), dimension(:,:), allocatable :: eigenvalues_Im_old
!
      real(dp), dimension(:,:), allocatable :: solution_vectors_reduced
      real(dp), dimension(:,:), allocatable :: solution
!
      integer(i15) :: state = 0, unit_solution = 0, ioerror = 0 ! For looping over the states
!
      real(dp) :: start_excited_state_solver, end_excited_state_solver
      real(dp) :: start_excited_state_iter, end_excited_state_iter
!
!     Start timings
!
      call cpu_time(start_excited_state_solver)
!
!     Test for triplet calculation, and stop if so - not yet implemented
!
      if (.not. wf%tasks%n_triplet_states .eq. 0) then
         write(unit_output,'(/t3,a/)') 'Triplet excitations not implemented.'
      endif
      flush(unit_output)
!
      converged = .false. 
!
      converged_energy   = .false.
      converged_residual = .false.
!
      iteration = 1
!
!     Initialize for excited state calculation
!
      reduced_dim  = wf%tasks%n_singlet_states
      n_new_trials = wf%tasks%n_singlet_states
!
!
      call wf%initialize_excited_states
!
!     Find start trial vectors and store them to the trial_vec file
!
      call wf%initialize_trial_vectors
!
!     Allocate and initialize eigenvalue arrays
!
      call allocator(eigenvalues_Im_old, wf%tasks%n_singlet_states, 1)
      call allocator(eigenvalues_Re_old, wf%tasks%n_singlet_states, 1)
      call allocator(eigenvalues_Im_new, wf%tasks%n_singlet_states, 1)
      call allocator(eigenvalues_Re_new, wf%tasks%n_singlet_states, 1)
!
      eigenvalues_Im_old = zero
      eigenvalues_Re_old = zero
      eigenvalues_Im_new = zero
      eigenvalues_Re_new = zero
!
!     Enter iterative loop
!
      do while (.not. converged .and. iteration .le. max_iterations) 
         if (timings) call cpu_time(start_excited_state_iter)
!
!        Prints 
!
         write(unit_output,'(/t3,a,i3)') 'Iteration:', iteration
         write(unit_output,'(t3,a,i3/)') 'Reduced space dimension:', reduced_dim
         flush(unit_output)
!
!        Transform new trial vectors  
!        rho_i = A * c_i
!
         call wf%transform_trial_vectors(reduced_dim - n_new_trials + 1, reduced_dim)
!
!        Allocate solution vectors for reduced problem
!
         call allocator(solution_vectors_reduced, reduced_dim, wf%tasks%n_singlet_states)
         solution_vectors_reduced = zero
!
!        Solve the reduced eigenvalue problem
!
         call wf%solve_reduced_eigenvalue_equation(eigenvalues_Re_new, eigenvalues_Im_new, &
                                                   solution_vectors_reduced, reduced_dim, n_new_trials)
!
!        Test for energy convergence of all states 
!
         converged_energy = .true.
!
         do state = 1, wf%tasks%n_singlet_states
!
            if (abs(eigenvalues_Re_new(state,1)-eigenvalues_Re_old(state,1)) .gt. wf%settings%energy_threshold) then
!
               converged_energy = .false.
!
            endif
!
         enddo
!
!        Save excitation energies for next iteration 
!
         call dcopy(wf%tasks%n_singlet_states, eigenvalues_Im_new, 1, eigenvalues_Im_old, 1)
         call dcopy(wf%tasks%n_singlet_states, eigenvalues_Re_new, 1, eigenvalues_Re_old, 1)
!
!        Get next trial vectors and test for convergence of residuals 
!
         call wf%construct_next_trial_vectors(eigenvalues_Re_new, eigenvalues_Im_new, &
                                                solution_vectors_reduced, reduced_dim, n_new_trials)
!
!        Test for convergence of the energies and residuals
!
         if (converged_energy .and. converged_residual) then
            converged = .true.
         else
            iteration = iteration + 1
         endif
!
         call deallocator(solution_vectors_reduced, reduced_dim, wf%tasks%n_singlet_states)
!
         if (timings) then
            call cpu_time(end_excited_state_iter)
            write(unit_output,'(t3,a35,i5,a5,f14.8/)') 'Total time (seconds) of iteration ',&
                                     iteration, ' : ' ,end_excited_state_iter - start_excited_state_iter
            flush(unit_output)
!
         endif
!
      enddo ! End of iterative loop 
!
!     Prints
!
      if ( converged ) then
         write(unit_output,'(/t3,a,i2,a/)')  'Converged in ', iteration, ' iterations!'
!
         if (print_vectors) then
!
            call wf%print_excited_state_info
!
         endif  
!
      else
         write(unit_output,'(/t3,a/)') 'Max number of iterations performed without convergence!'
      endif
!
!     Final deallocations
!
      call deallocator(eigenvalues_Im_old, reduced_dim, 1)
      call deallocator(eigenvalues_Re_old, reduced_dim, 1)
      call deallocator(eigenvalues_Im_new, reduced_dim, 1)
      call deallocator(eigenvalues_Re_new, reduced_dim, 1)
!
!     End and print timings
!
      call cpu_time(end_excited_state_solver)
!
      write(unit_output,'(t3,a27,f14.8/)') 'Total time (seconds):', end_excited_state_solver - start_excited_state_solver
      flush(unit_output)
!
   end subroutine excited_state_solver_ccs
!
!
   module subroutine solve_reduced_eigenvalue_equation_ccs(wf, eigenvalues_Re, eigenvalues_Im, &
                                                                solution_vectors_reduced, reduced_dim, n_new_trials)
!!
!!    Solve Reduced Eigenvalue Equation
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!    Constructs the reduced A matrix, solves its eigenvalue equation,
!!    and returns its first n eigenvalues and eigenvectors (reduced space
!!    solution vectors).
!!
      implicit none
!
      class(ccs) :: wf
!
      integer(i15) :: reduced_dim, n_new_trials
!
      real(dp), dimension(wf%tasks%n_singlet_states,1) :: eigenvalues_Re
      real(dp), dimension(wf%tasks%n_singlet_states,1) :: eigenvalues_Im
! 
      real(dp), dimension(reduced_dim, wf%tasks%n_singlet_states) :: solution_vectors_reduced
!
      real(dp), dimension(:,:), allocatable :: A_red
      real(dp), dimension(:,:), allocatable :: c_i
      real(dp), dimension(:,:), allocatable :: rho_j
      real(dp), dimension(:,:), allocatable :: solution_vectors_reduced_all 
      real(dp), dimension(:,:), allocatable :: eigenvalues_Re_all
      real(dp), dimension(:,:), allocatable :: eigenvalues_Im_all
      real(dp), dimension(:,:), allocatable :: work
!
      real(dp) :: ddot, dummy
!
      integer(i15), dimension(:,:), allocatable :: index_list
!
      integer(i15) :: unit_trial_vecs = 0, unit_rho = 0, unit_reduced_jacobi = 0, ioerror = 0 
!
      integer(i15) :: i = 0, j = 0
!
      integer :: info = -1 
!
      call allocator(A_red, reduced_dim, reduced_dim)
      A_red = zero
!
!     -::- Prepare to solve the eigenvalue problem -::-
!
!     Prepare files
!  
      call generate_unit_identifier(unit_trial_vecs)
         open(unit=unit_trial_vecs, file='trial_vec', action='read', status='unknown', &
           access='direct', form='unformatted', recl=dp*(wf%n_parameters), iostat=ioerror)
!
      call generate_unit_identifier(unit_rho)
         open(unit=unit_rho, file='transformed_vec', action='read', status='unknown', &
           access='direct', form='unformatted', recl=dp*(wf%n_parameters), iostat=ioerror)      
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
         read(unit_reduced_jacobi) ((A_red(i,j),i = 1, reduced_dim-n_new_trials), j=1, reduced_dim-n_new_trials)
!
      endif
!
!     Allocate c and rho
!
      call allocator(c_i, wf%n_parameters, 1)
      call allocator(rho_j, wf%n_parameters, 1)
      c_i   = zero
      rho_j = zero
!
!     Construct reduced Jacobi matrix 
!     A_red_ij = c_i^T * A * c_j = c_i^T * rho_i
!
      if (iteration .eq. 1) then
!
!        If first iteration, make the entire reduced matrix 
!
         do i = 1,reduced_dim
           read(unit_trial_vecs, rec=i, iostat=ioerror) c_i
            do j = 1,reduced_dim
               read(unit_rho, rec=j, iostat=ioerror) rho_j
               A_red(i,j) = ddot(wf%n_parameters, c_i, 1, rho_j, 1)
            enddo
         enddo
!
      else
!
!        If not, make the elements of the reduced matrix not constructed
!        in previous iterations 
!
         do i = 1,reduced_dim
!
            read(unit_trial_vecs, rec=i, iostat=ioerror) c_i
!
            do j = reduced_dim - n_new_trials + 1, reduced_dim
               read(unit_rho, rec=j, iostat=ioerror) rho_j
               A_red(i,j) = ddot(wf%n_parameters, c_i, 1, rho_j, 1)
            enddo
!
         enddo
!
         do j = 1, reduced_dim - n_new_trials
!
            read(unit_rho, rec=j, iostat=ioerror) rho_j
!
            do i = reduced_dim - n_new_trials + 1, reduced_dim
               read(unit_trial_vecs, rec=i, iostat=ioerror) c_i
               A_red(i,j) = ddot(wf%n_parameters, c_i, 1, rho_j, 1)
            enddo
!
         enddo
!
      endif
!
      call deallocator(c_i, wf%n_parameters, 1)
      call deallocator(rho_j, wf%n_parameters, 1)
!
!     Close files for trial vectors and transformed vectors
!
      close(unit_trial_vecs)
      close(unit_rho)
!
!     Write current reduced Jacobi matrix to file
!
      rewind(unit_reduced_jacobi)
      write(unit_reduced_jacobi) ((A_red(i,j), i = 1, reduced_dim), j = 1, reduced_dim)
!
!     Close reduced Jacobi file
!
      close(unit_reduced_jacobi)
!
!     -::- Solve reduced eigenvalue problem -::-
!
!     Allocate arrays for eigenvalues and eigenvectors
!
      call allocator(solution_vectors_reduced_all, reduced_dim, reduced_dim)
      solution_vectors_reduced_all = zero
!
      call allocator(eigenvalues_Re_all, reduced_dim, 1)
      call allocator(eigenvalues_Im_all, reduced_dim, 1)
      eigenvalues_Re_all = zero
      eigenvalues_Im_all = zero
!
      call allocator(work, 4*reduced_dim, 1)
      work = zero
!
!     Solve reduced eigenvalue problem
!
      info = 0
      call dgeev('N','V',                       &
                  reduced_dim,                  &
                  A_red,                        &
                  reduced_dim,                  &
                  eigenvalues_Re_all,           &
                  eigenvalues_Im_all,           &
                  dummy,                        &
                  1,                            &
                  solution_vectors_reduced_all, &
                  reduced_dim,                  &
                  work,                         &
                  4*reduced_dim,                &
                  info)
      if (info .ne. 0) then 
         write(unit_output,*)  'WARNING: Error while finding solution', info
         stop
      endif
!
      call deallocator(work, 4*reduced_dim, 1)
!
!     Deallocate reduced Jacobi matrix
!
      call deallocator(A_red, reduced_dim, reduced_dim)
!
!     Find lowest eigenvalues and sort them (the corresponding indices
!     are placed in the integer array index_list)
!
      call allocator_int(index_list, wf%tasks%n_singlet_states, 1)
      index_list = 0
!
      call get_n_lowest(wf%tasks%n_singlet_states, reduced_dim, eigenvalues_Re_all, eigenvalues_Re, index_list)
!
!     Pick out solution vectors and imaginary parts of eigenvalues according to index_list 
!
      do i = 1, reduced_dim
         do j = 1, wf%tasks%n_singlet_states
!
            solution_vectors_reduced(i,j) = solution_vectors_reduced_all(i,index_list(j,1))
            eigenvalues_Im = eigenvalues_Im_all(index_list(j,1), 1)
!
         enddo
      enddo
!
!     Final deallocatons
!
      call deallocator(solution_vectors_reduced_all, reduced_dim, reduced_dim)
      call deallocator(eigenvalues_Im_all, reduced_dim, 1)
      call deallocator(eigenvalues_Re_all, reduced_dim, 1)
!
      call deallocator_int(index_list,wf%tasks%n_singlet_states,1)
!
   end subroutine solve_reduced_eigenvalue_equation_ccs
!
!
   module subroutine construct_next_trial_vectors_ccs(wf, eigenvalues_Re, eigenvalues_Im, &
                                                solution_vectors_reduced, &
                                                reduced_dim, n_new_trials)
!!
!!    Construct Next Trial Vectors    
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!    Constructs the next eigenvectors by constructing the residual vectors
!!    
!!       R_j = (A*X_j - e*X_j)/|X_j|,
!!
!!    orthogonalizing them against the other trial vectors.
!!
!!    Residual vectors are preconditioned before orthogonalization.
!!    This is done by dividing by the orbital differences.
!!    
!!    If norm of orthogonal vector is very small 
!!    (i.e. high degree of linear dependence on previous trial vectors)
!!    it is scrapped. If norm sufficiently large, vector is normalized and
!!    stored in trial_vec file, to be used in the next iteration.
!!
!!    Routine also constructs full space solution vectors and stores them
!!    in file solution_vectors 
!! 
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%tasks%n_singlet_states,1) :: eigenvalues_Re
      real(dp), dimension(wf%tasks%n_singlet_states,1) :: eigenvalues_Im
!
      integer(i15) :: reduced_dim
      integer(i15) :: n_new_trials
!
      real(dp), dimension(reduced_dim, wf%tasks%n_singlet_states) :: solution_vectors_reduced

!
!     Local variables 
!
      integer(i15) :: root = 0, trial = 0, i = 0, j = 0
!
      real(dp), dimension(:,:), allocatable :: solution_vector
      real(dp), dimension(:,:), allocatable :: residual
      real(dp), dimension(:,:), allocatable :: orbital_diff
!
      real(dp), dimension(:,:), allocatable :: c_i
      real(dp), dimension(:,:), allocatable :: c_j
      real(dp), dimension(:,:), allocatable :: rho_i
!
      real(dp) :: norm_new_trial = zero
      real(dp) :: norm_residual = zero
      real(dp) :: norm_solution_vector = zero
!
      real(dp) :: conv_test = zero
      real(dp) :: dot_prod = zero  
!
      integer(i15) :: ioerror = 0
      integer(i15) :: unit_trial_vecs = 0, unit_rho = 0, unit_solution = 0 ! Unit identifiers for files 
!
      real(dp) :: ddot
      character(100) :: iostring
!
!     Prepare necessary files
!
      call generate_unit_identifier(unit_trial_vecs)
      open(unit=unit_trial_vecs, file='trial_vec', action='readwrite', status='unknown', &
        access='direct', form='unformatted', recl=dp*(wf%n_parameters), iostat=ioerror)
      if (ioerror .ne. 0) write(unit_output,*) 'Error while opening trial vecs file'
!
      call generate_unit_identifier(unit_rho)
      open(unit=unit_rho, file='transformed_vec', action='read', status='unknown', &
      access='direct', form='unformatted', recl=dp*(wf%n_parameters), iostat=ioerror) 
      if (ioerror .ne. 0) write(unit_output,*) 'Error while opening transformed vecs file'
!
      call generate_unit_identifier(unit_solution)
      open(unit=unit_solution, file=wf%excited_state_task, action='write', status='unknown', &
      access='direct', form='unformatted', recl=dp*(wf%n_parameters), iostat=ioerror) 
      if (ioerror .ne. 0) write(unit_output,*) 'Error while opening solution file'
!
      call allocator(residual, wf%n_parameters, 1)
!
      n_new_trials = 0
      converged_residual = .true.
!
      write(unit_output,'(t3,a)') 'Root       Eigenvalue (Re)      Eigenvalue (Im)      Residual norm'
      write(unit_output,'(t3,a)') '-------------------------------------------------------------------'
!
!     For each of the roots
!
      do root = 1, wf%tasks%n_singlet_states
!
         residual = zero
!
!        :: Create fullspace vector X and calculate norm ||X|| ::
!
         call allocator(solution_vector, wf%n_parameters, 1)
         solution_vector = zero
!
         call allocator(c_i, wf%n_parameters, 1)
!
!        Calculate X_j = sum_i x_j_i * c_i (for the root j)
!
         do trial = 1, reduced_dim
!
            c_i = zero
            read(unit_trial_vecs, rec=trial, iostat=ioerror, iomsg = iostring) c_i
            call daxpy(wf%n_parameters, solution_vectors_reduced(trial,root), c_i, 1, solution_vector, 1)
!
            if (ioerror .ne. 0) write(unit_output,*) 'Error reading trial vecs in get_next_trial_vectors', ioerror, iostring
!
         enddo
!
         call deallocator(c_i, wf%n_parameters, 1)
!
!        Calculate norm of solution vector
!
         norm_solution_vector = sqrt(ddot(wf%n_parameters, solution_vector, 1, solution_vector, 1))
!
!        Normalize and write full space solution vectors to file
!
         call dscal(wf%n_parameters, one/norm_solution_vector, solution_vector, 1)
         write(unit_solution, rec=root, iostat=ioerror) solution_vector
!
!        :: Calculate residual ::
!
         call dcopy(wf%n_parameters, solution_vector, 1, residual, 1) ! R = X 
!        
         call deallocator(solution_vector, wf%n_parameters, 1)
!
!        Residual = AX - eX
!
!        - eX : 
!
         call dscal(wf%n_parameters, -eigenvalues_Re(root, 1), residual, 1)
!
!        + AX : 
!
         call allocator(rho_i, wf%n_parameters, 1)
!
         do trial = 1, reduced_dim
!
            rho_i = zero
            read(unit_rho, rec=trial, iostat=ioerror) rho_i
            call daxpy(wf%n_parameters, solution_vectors_reduced(trial,root), rho_i, 1, residual, 1)
!
            if (ioerror .ne. 0) write(unit_output,*) 'Error reading tranf vecs in get_next_trial_vectors', ioerror
!
         enddo
!
         call deallocator(rho_i, wf%n_parameters, 1)
!
!        Calculate norm of residual || AX - eX ||
!
         norm_residual = sqrt(ddot(wf%n_parameters, residual, 1, residual, 1))
!
!        Calculate residual norm and check convergence criteria on residual norms
!        ||AX-eX||/||X||
!
         conv_test = norm_residual/norm_solution_vector
!
         if (conv_test .gt. wf%settings%equation_threshold) converged_residual = .false.
!
!        Prints
!
         write(unit_output,'(t3,i2,5x,f14.8,7x,f14.8,11x,e10.4)') root, eigenvalues_Re(root, 1), &
                                                                eigenvalues_Im(root, 1), norm_residual/norm_solution_vector
         flush(unit_output)
!
!        :: Precondition the residual by inverse orbital energy differences ::
!
         call wf%precondition_residual(residual)
!
!        :: Orthogonalize the residual to the other trial vectors ::
!
!        Calculate norm of the residual
!
         norm_residual = sqrt(ddot(wf%n_parameters, residual, 1, residual, 1))
!
!        Normalize the residual
!
         call dscal(wf%n_parameters, one/norm_residual, residual, 1)
!
!        Orthogonalize new trial vector (the residual) against old trial vectors
!
         call allocator(c_i, wf%n_parameters, 1)
!
!        prod_i (I - c_i*c_i^T)*Res = prod_i (Res - c_i*c_i^T*Res)
!
         do trial = 1, reduced_dim + n_new_trials
!
            c_i = zero
            read(unit_trial_vecs, rec=trial, iostat=ioerror) c_i
            if (ioerror .ne. 0) write(unit_output,*) 'Error reading trial vecs in get_next_trial_vectors'
!
            dot_prod = ddot(wf%n_parameters, c_i, 1, residual, 1)
            call daxpy(wf%n_parameters, -dot_prod, c_i, 1, residual,1)
!
         enddo
!
         call deallocator(c_i, wf%n_parameters, 1)
!
!        Calculate norm of residual
!
         norm_new_trial = sqrt(ddot(wf%n_parameters, residual, 1, residual, 1))
!
!        Test for linear dependency on old trial vectors
!        If norm sufficiently high new vector is normalized and written to file
!
         if ((norm_new_trial .gt. wf%settings%equation_threshold) .and. (conv_test .gt. wf%settings%equation_threshold)) then
!
            n_new_trials = n_new_trials + 1
            call dscal(wf%n_parameters, one/norm_new_trial, residual, 1)
            write(unit_trial_vecs, rec=n_new_trials+reduced_dim, iostat=ioerror) residual
!     
         endif
!  
      enddo
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
      call deallocator(residual, wf%n_parameters, 1)
!
      write(unit_output,'(t3,a)') '-------------------------------------------------------------------'
!
   end subroutine construct_next_trial_vectors_ccs
!
!
   module subroutine initialize_trial_vectors_ccs(wf)
!!
!!    Initialize trial vectors (CCS)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!    Wrapper for initialization of trial vectors:
!!
!!    If restart, then checks if old solution file exists, 
!!    then uses old solutions as new trial vectors
!!
!!    If not restart: 
!!    initialize_trial_vectors_valence is called for regular excited state calculation
!!    initialize_trial_vectors_core is called for cvs calculation
!!     
!!
      implicit none
!
      class(ccs) :: wf
!     
!
!     If restart use old solution vectors for first start vectors
!
      if (wf%settings%restart) then 
!
         write(unit_output,'(/t3,a)') 'Requested restart. Using old solution vectors as trial vectors.'
         call wf%trial_vectors_from_stored_solutions
         return
!
      endif
!
!     Test for ionization or excitation
!
      if (wf%tasks%excited_state .or. wf%tasks%core_excited_state) then ! Excitation
!
!        Test for valence or core
!
         if ((wf%excited_state_task .eq. 'right_valence') .or. &
            (wf%excited_state_task .eq. 'left_valence')) then
!
            call wf%initialize_trial_vectors_valence
!
         elseif (wf%excited_state_task .eq. 'right_core') then
!
            call wf%initialize_trial_vectors_core
!
         endif
!
      elseif (wf%tasks%ionized_state .or. wf%tasks%core_ionized_state) then ! Ionization
!
!        Test for valence or core
!
         if ((wf%excited_state_task .eq. 'right_valence') .or. &
            (wf%excited_state_task .eq. 'left_valence')) then
!
            call wf%initialize_trial_vectors_valence_ionization
!
         elseif (wf%excited_state_task .eq. 'right_core') then
!
            call wf%initialize_trial_vectors_core_ionization
!
         endif
      endif
!
   end subroutine initialize_trial_vectors_ccs
!
!
   module subroutine initialize_trial_vectors_valence_ccs(wf)
!!
!!    Initialize trial vectors valence
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad
!!
!!    Initializes start trial vectors for the calculation of 
!!    singlet excited states and writes them to file 'trial_vecs'.
!!    Initializes start trial vectors for the calculation of 
!!    singlet excited states and writes them to file 'trial_vecs'.
!!
!!    n start vectors are constructed by finding the n lowest orbital differences,      
!!    where n = n_singlet_states. Vector i has a 1.0D0 at the element corresponding to the i'th lowest
!!    orbital difference and 0.0d0 everywhere else
!!
      implicit none
!
      class(ccs) :: wf
!       
      integer(i15), dimension(:,:), allocatable :: index_lowest_obital_diff
!
      real(dp), dimension(:,:), allocatable :: c
! 
      integer(i15) :: i = 0, j = 0
!
      integer(i15) :: unit_trial_vecs = 0, unit_rho = 0, ioerror = 0
!
!     Allocate array for the indices of the lowest orbital differences
!
      call allocator_int( index_lowest_obital_diff, wf%tasks%n_singlet_states, 1)
      index_lowest_obital_diff = zero

!
!     Find indecies of lowest orbital differences
!
      call wf%find_start_trial_indices(index_lowest_obital_diff)
!
!     Generate start trial vectors c and write to file
!
      call allocator(c, wf%n_parameters, 1)
!
!     Prepare for writing trial vectors to file
!
      call generate_unit_identifier(unit_trial_vecs)
      open(unit=unit_trial_vecs, file='trial_vec', action='write', status='unknown', &
           access='direct', form='unformatted', recl=dp*(wf%n_parameters), iostat=ioerror)
!
      do i = 1, (wf%tasks%n_singlet_states)
         c = zero
         c(index_lowest_obital_diff(i,1),1) = one
         write(unit_trial_vecs, rec=i, iostat=ioerror) (c(j,1), j = 1, wf%n_parameters)
      enddo
!
!     Close file
!     
      close(unit_trial_vecs)
!
!     Deallocate c
!
      call deallocator(c, wf%n_parameters, 1)
!
!     Deallocate index_lowest_obital_diff
!
      call deallocator_int( index_lowest_obital_diff, wf%tasks%n_singlet_states, 1)
!
      end subroutine initialize_trial_vectors_valence_ccs
!
!
   module subroutine trial_vectors_from_stored_solutions_ccs(wf)
!!
!!    Trial vectors from old solutions,
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, June 2017
!!
!!    Restart: Use old solutions as trial vectors
!!    Trial Vectors from Stored Solutions (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!    Reads the solutions from file and uses them as the first trial
!!    vectors in the iterative loop.
!!
      implicit none
!
      class(ccs) :: wf
!
      logical      :: solution_exists = .false.
      logical      :: more_trials = .true.
!
      integer(i15) :: ioerror = 0, unit_solution = 0, unit_trial_vecs = 0
      integer(i15) :: number_of_solutions = 0
!
      integer(i15) :: i = 0, j = 0
!
      real(dp), dimension(:,:), allocatable :: c_i, c_j
!
      real(dp) :: ddot, dot_prod = zero, norm = zero 
!
!     Open solution vector file - if it does not exist return
!
      inquire(file=wf%excited_state_task, exist=solution_exists)
!
!     If no solution vector file, return and use orbital differences.
!
      if (.not. solution_exists) return
!
!     Open files
!
      call generate_unit_identifier(unit_trial_vecs)
      open(unit=unit_trial_vecs, file='trial_vec', action='readwrite', status='old', &
         access='direct', form='unformatted', recl=dp*(wf%n_parameters), iostat=ioerror)
!
      call generate_unit_identifier(unit_solution)
!
      open(unit=unit_solution, file=wf%excited_state_task, action='read', status='unknown', &
         access='direct', form='unformatted', recl=dp*(wf%n_parameters), iostat=ioerror) 
!
!     Allocate c_i
!
      call allocator(c_i, wf%n_parameters, 1)
      c_i = zero
!
      i = 1
!
      do while ((i .le. wf%tasks%n_singlet_states) .and. more_trials)
!
!        Read old solutions and count them
!
         read(unit_solution, rec=i, iostat=ioerror) c_i
         if (ioerror .ne. 0) write(unit_output,*) 'Error reading solution vecs'
!
         if (ioerror .eq. 0) then
!
            write(unit_trial_vecs, rec = i) c_i
!
         else
!
            more_trials = .false.
!
         endif
!
         i = i + 1
!
      enddo
!
!     Deallocate c_i 
!
      call deallocator(c_i, wf%n_parameters, 1)
!
!     Close solution file
!
      close(unit_solution)
!
!     Allocate c_i and c_j
!
      call allocator(c_i, wf%n_parameters, 1)
      call allocator(c_j, wf%n_parameters, 1)
      c_i = zero
      c_j = zero
!
!     Reorthonormalize trial vectors
!
      do i = 1, wf%tasks%n_singlet_states
!
         read(unit_trial_vecs, rec=i, iostat=ioerror) c_i
!
         do j = 1, i-1
!
            read(unit_trial_vecs, rec=j, iostat=ioerror) c_j
            dot_prod = ddot(wf%n_parameters, c_j, 1, c_i, 1)
            call daxpy(wf%n_parameters, -dot_prod, c_j, 1, c_i, 1)
!  
            norm = sqrt(ddot(wf%n_parameters, c_i, 1, c_i, 1))
            call dscal(wf%n_parameters, one/norm, c_i, 1)
!
         enddo
         write(unit_trial_vecs, rec = i)c_i
      enddo
!  
      call deallocator(c_i, wf%n_parameters, 1)
      call deallocator(c_j, wf%n_parameters, 1)  
!
!     Close trial vector file
!
      close(unit_trial_vecs)     
!
   end subroutine trial_vectors_from_stored_solutions_ccs
!
!
   module subroutine find_start_trial_indices_ccs(wf, index_list)
!!
!!    Find Start Trial Indices (CCS) 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
      implicit none
!
      class(ccs) :: wf
      integer(i15), dimension(wf%tasks%n_singlet_states,1), intent(inout) :: index_list
!
      real(dp), dimension(:,:), allocatable :: orbital_diff
      real(dp), dimension(:,:), allocatable :: lowest_orbital_diff
!
      integer(i15) :: a = 0, i = 0, j = 0
!
      integer(i15) :: ai = 0
!
      real(dp)     :: max
      integer(i15) :: max_pos 
!
      real(dp)     :: swap     = zero
      integer(i15) :: swap_int = 0
!
!     Test if there are user specified trial vectors
!
      if (wf%tasks%user_specified_start_vector) then
         index_list = wf%tasks%start_vectors
      else
!
!        Allocate orbital_diff
!
         call allocator(orbital_diff,wf%n_parameters,1)
         orbital_diff = zero
!
!        Calculate orbital differences
!
         call wf%calculate_orbital_differences(orbital_diff)
!
!        Finding lowest orbital differences
!
         call allocator(lowest_orbital_diff, wf%tasks%n_singlet_states, 1)
!        
         lowest_orbital_diff = zero
!
         call get_n_lowest(wf%tasks%n_singlet_states, wf%n_parameters, orbital_diff, lowest_orbital_diff, index_list)
!
         call deallocator(orbital_diff,wf%n_parameters,1)
!
         call deallocator(lowest_orbital_diff, wf%tasks%n_singlet_states, 1)
!
      endif
!
   end subroutine find_start_trial_indices_ccs
!
!
   module subroutine calculate_orbital_differences_ccs(wf,orbital_diff)

!!
!!    Calculate Orbital Differences (CCS)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad May 2017
!!
      implicit none
!
      class(ccs) :: wf
      real(dp), dimension(wf%n_parameters, 1) :: orbital_diff
!
      integer(i15) :: a = 0, i = 0
      integer(i15) :: ai = 0
!
      do i = 1, wf%n_o
         do a = 1, wf%n_v
            ai = index_two(a, i, wf%n_v)
            orbital_diff(ai, 1) = wf%fock_diagonal(a + wf%n_o, 1) - wf%fock_diagonal(i, 1)
         enddo
      enddo
!
   end subroutine calculate_orbital_differences_ccs
!
!
   module subroutine transform_trial_vectors_ccs(wf, first_trial, last_trial)
!!
!!    Transform trial vectors (CCS)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!    Each trial vector in first_trial to last_trial is read from file and
!!    transformed before the transformed vector is written to file.
!!
      implicit none
!
      class(ccs) :: wf
!
      integer(i15), intent(in) :: first_trial, last_trial ! Which trial_vectors we are to transform
!
      real(dp), dimension(:,:), allocatable :: c_a_i
!
      integer(i15) :: unit_trial_vecs = 0, unit_rho = 0, ioerror = 0
      integer(i15) :: trial = 0 
!
!     Allocate c_a_i
!
      call allocator(c_a_i, wf%n_v, wf%n_o)
      c_a_i = zero 
!
!     Open trial vector and transformed vector files
!
      call generate_unit_identifier(unit_trial_vecs)
      open(unit=unit_trial_vecs, file='trial_vec', action='read', status='unknown', &
        access='direct', form='unformatted', recl=dp*(wf%n_v)*(wf%n_o), iostat=ioerror)
!
      call generate_unit_identifier(unit_rho)
      open(unit=unit_rho, file='transformed_vec', action='write', status='unknown', &
        access='direct', form='unformatted', recl=dp*(wf%n_v)*(wf%n_o), iostat=ioerror)
!
!     For each trial vector: Read, transform and write
!               
      do trial = first_trial, last_trial
!
         read(unit_trial_vecs, rec=trial, iostat=ioerror) c_a_i
!
         if (wf%current_task == 'excited_state') then
!
            if (wf%excited_state_task == 'right_valence' .or. wf%excited_state_task == 'right_core') then
!
               call wf%jacobian_ccs_transformation(c_a_i)
!
            elseif (wf%excited_state_task == 'left_valence') then
!               
               call wf%jacobian_transpose_ccs_transformation(c_a_i)
!
            else
!
               write(unit_output,*) 'Error: Excited state task not recognized'
               stop
!
            endif
!
         elseif(wf%current_task == 'response') then
!
            if (wf%response_task =='left_eigenvectors') then
!
               call wf%jacobian_transpose_ccs_transformation(c_a_i)
!
            elseif (wf%response_task == 'multipliers') then 
!
               call wf%jacobian_transpose_ccs_transformation(c_a_i)
!
            else
!
               write(unit_output,*) 'Error: Response task not recognized'
               stop
!
            endif
!
         else
!
            write(unit_output,*) 'Error: Current task not recognized'
            stop
!
         endif
!
!        -::- Projections -::-
!
!        Test for core calculation 
!
         if (wf%tasks%core_excited_state .or. wf%tasks%core_ionized_state) then
!  
!           Project out contamination from valence contributions
!
            call wf%cvs_rho_a_i_projection(c_a_i)
!
         endif
!
!        Test for ionization calculation
!
         if (wf%tasks%ionized_state .or. wf%tasks%core_ionized_state) then
!
!           Project out contamination from regular excitations
!
            call wf%ionization_rho_a_i_projection(c_a_i)
!
         endif
!
!        Write transformed vector to file
!
         write(unit_rho, rec=trial, iostat=ioerror) c_a_i
          
      enddo
!
      close(unit_trial_vecs) 
      close(unit_rho)                                
!
!     Deallocate c_a_i
!
      call deallocator(c_a_i, wf%n_v, wf%n_o)
!
   end subroutine transform_trial_vectors_ccs
!
!
      module subroutine initialize_excited_states_ccs(wf)
!!
         implicit none 
!    
         class(ccs) :: wf
!
         call wf%initialize_amplitudes
!
      end subroutine initialize_excited_states_ccs
!
!
      module subroutine precondition_residual_ccs(wf, residual)
!!
!!       Precondition residual
!!       Written by Sarai D. Folkestad, Aug. 2017
!!
!!       Calls precondition_residual_valence for normal excited state calculation
!!       Calls precondition_residual_core for cvs calculation
!!
         implicit none
!
         class(ccs) :: wf
         real(dp), dimension(wf%n_parameters ,1) :: residual
!
         if (wf%tasks%excited_state .or. wf%tasks%core_excited_state) then       
!
            if ((wf%excited_state_task .eq. 'right_valence') .or. &
                (wf%excited_state_task .eq. 'left_valence')) then
!
               call wf%precondition_residual_valence(residual)
!
            elseif (wf%excited_state_task .eq. 'right_core') then
!
               call wf%precondition_residual_core(residual)
!
            endif
!
         elseif (wf%tasks%ionized_state) then
!
               if ((wf%excited_state_task .eq. 'right_valence') .or. &
                (wf%excited_state_task .eq. 'left_valence')) then
!
               call wf%precondition_residual_valence_ionization(residual)
!
            elseif (wf%excited_state_task .eq. 'right_core') then
!
               call wf%precondition_residual_core_ionization(residual)
!
            endif
!
         endif
!
      end subroutine precondition_residual_ccs
!
!
      module subroutine precondition_residual_valence_ccs(wf, residual)
!!
!!       Precondition residual valence
!!       Written by Sarai D. Folkestad, Aug. 2017
!!
!!       Divide elements of residual by orbital difference       
!!
         implicit none
!
         class(ccs) :: wf
         real(dp), dimension(wf%n_parameters ,1) :: residual
!
         integer(i15) :: i = 0
!
         real(dp), dimension(:,:), allocatable :: orbital_diff
!   
         call allocator(orbital_diff, wf%n_parameters, 1)
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
         call deallocator(orbital_diff, wf%n_parameters, 1)
!
      end subroutine precondition_residual_valence_ccs
!
!
     
      module subroutine print_excited_state_info_ccs(wf)
!!
!!
!!
         implicit none
!  
         class(ccs) :: wf
!
         integer(i15) :: unit_solution = -1, unit_es_info = -1, ioerror = 0
         integer(i15) :: state = 0
!
         real(dp), dimension(:,:), allocatable :: solution
!  
!        Read solution vectors 
!  
         call generate_unit_identifier(unit_solution)
!
         open(unit=unit_solution, file=wf%excited_state_task, action='read', status='unknown', &
         access='direct', form='unformatted', recl=dp*(wf%n_parameters), iostat=ioerror) 
!
         if (ioerror .ne. 0) write(unit_output,*) 'Error while opening solution file'
!
!        Open info file 
!  
         call generate_unit_identifier(unit_es_info)
!
         open(unit=unit_es_info, file='excited_state_information', action='write', status='unknown', &
         access='sequential', form='formatted', iostat=ioerror)
         rewind(unit_es_info) 
!
         if (ioerror .ne. 0) write(unit_output,*) 'Error while opening excited_state_information file'

         call allocator(solution, wf%n_parameters, 1)
         do state = 1, wf%tasks%n_singlet_states
!  
            solution = zero
            read(unit_solution, rec=state) solution
!  
            write(unit_es_info,'(/a33)')'----------------------------------'
            write(unit_es_info,'(a30,i2, a1)')'Components of solution vector', state, ':'
            write(unit_es_info,'(a33/)')'----------------------------------'
            call wf%print_excitation_vector(solution, unit_es_info)
!  
         enddo
!
         call deallocator(solution, wf%n_parameters,1) 
!
         close(unit_solution)
         close(unit_es_info)
!
      end subroutine print_excited_state_info_ccs
!
!
      module subroutine print_excitation_vector_ccs(wf, vec, unit_id)
!!
!!
!!
         implicit none
!  
         class(ccs) :: wf
!
         real(dp), dimension(wf%n_parameters, 1) :: vec
!
         integer(i15) :: unit_id     
!
         integer(i15) :: a = 0, i = 0, ai = 0
!
         write(unit_id,'(2a6,a12)')'a', 'i', 'coeff'
         write(unit_id,'(t3,a)')'-------------------------'
!
         do a = 1, wf%n_v
            do i = 1, wf%n_o
!  
               ai = index_two(a, i, wf%n_v)
               if (abs(vec(ai, 1)) .gt. 1.0D-03) then
                  write(unit_id,'(2i6,f12.4)') a, i, vec(ai, 1)
               endif
!
         enddo
      enddo
!
      end subroutine print_excitation_vector_ccs
!
!
end submodule excited_state
