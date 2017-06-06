submodule (ccs_class) excited_state
!
!!
!!    Excited state  submodule (CCS)
!!    Written by Eirik F. Kjønstad and Sarai Dery Folkestad, May 2017
!!
!!    Contains the following family of procedures of the CCS class:
!!
!
   implicit none 
!
!  Some variables available to all routines of the module
!
   integer(i15) :: iteration = 1, max_iterations = 50
!
!     Variables to handle convergence criterea
!
   logical :: converged = .false. ! True iff both the energy and the equations have converged 
!
   logical :: converged_energy   = .false.
   logical :: converged_residual = .false.
!
!
contains
!
!
   module subroutine excited_state_solver_ccs(wf)
!!
!!    Excited State Solver
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad May 2017
!!
!!    Directs the solution of the excited states using a Davidson algorithm.
!!    The routine aims to find the right eigenvectors of the Jacobian matrix
!!
!!       AX = eX,  
!!
!!    and the eigenvalues which corresponds to the excitation energies.
!!
!!    The problem is solved in reduced space. To find n roots, n start trial vectors {c_i}_i=1,n are 
!!    generated according to the lowest orbital differences. Then a reduced space Jacobian is 
!!    constructed 
!! 
!!       A_red_ij = c_i^T * A c_j
!!
!!    and the eigenvalues e and eigenvectors x of this matrix are found.
!!    The n full space vectors {X_j}_j=1,n are then given by
!!
!!       X_j = sum_i x_j_i*c_i, 
!! 
!!    and the j'th residual vector is given by
!!
!!       R_j = (A*X_j - e*X_j)/|X_j|,
!!
!!    and only if the norm of this residual is sufficiently small 
!!    (and the excitation energies are converged within a given threshold)
!!    is convergence reached. If not, then new trial vectors will be generated 
!!    by orthogonalizing the residual vector against the previous trial vectors
!!    and then normalizing them, and consequently expanding the dimension 
!!    of the reduced space for the next iteration.
!!   

      implicit none
!
      class(ccs) :: wf
!
      integer(i15) :: reduced_dim  = 0 
      integer(i15) :: n_new_trials = 0
!
!  Solution for the reduced eigenvalue problem
!
      real(dp), dimension(:,:), allocatable :: eigenvalues_Re_new
      real(dp), dimension(:,:), allocatable :: eigenvalues_Im_new
      real(dp), dimension(:,:), allocatable :: eigenvalues_Re_old
      real(dp), dimension(:,:), allocatable :: eigenvalues_Im_old
      real(dp), dimension(:,:), allocatable :: solution_vectors_reduced
!
      integer(i15) :: i = 0
!
      real(dp) :: start_excited_state_solver, end_excited_state_solver
!
!     Start timings
!
      call cpu_time(start_excited_state_solver)
!
!     Let the user know the excited state solver is running
!
      write(unit_output,'(/t3,a)')   ':: Excited state solver (Davidson)'
      write(unit_output,'(t3,a/)')   ':: E. F. Kjønstad, S. D. Folkestad, May 2017'
      write(unit_output,'(t3,a,i3,a,a,a)') &
                                     'Requested ',wf%tasks%n_singlet_states,' ', trim(wf%name), ' singlet states.'
      write(unit_output,'(t3,a,i3,a,a,a/)') &
                                     'Requested ',wf%tasks%n_triplet_states,' ', trim(wf%name), ' triplet states.'
!
!     Test for n_triplet_states - Not implemented
!
      if (.not. wf%tasks%n_triplet_states .eq. 0) then
         write(unit_output,'(t3,a/)') 'Triplet excitations not implemented.'
      endif
      flush(unit_output)
!
!     Initialize for excited state calculation
!
      reduced_dim        = wf%tasks%n_singlet_states
      n_new_trials = wf%tasks%n_singlet_states
!
!     Find start trial vectors and store them to file trial_vec
!
      call wf%initialize_trial_vectors
!
!     Allocate eigenvalue arrays
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
!     Start of iterative loop
!
      do while (.not. converged .and. iteration .le. max_iterations) 
!
!     Prints 
!
         write(unit_output,'(/t3,a,i3/)') 'Iter.', iteration
         write(unit_output,'(t3,a,i3)')'Reduced space dimension:', reduced_dim
         write(unit_output,'(t3,a/)')'----------------------------'
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
!        Test energy convergence criteria
!
         converged_energy = .true.
!
         do i = 1, wf%tasks%n_singlet_states
            if (abs(eigenvalues_Re_new(i,1)-eigenvalues_Re_old(i,1)) .gt. wf%settings%energy_threshold) then
               converged_energy = .false.
            endif
         enddo
!
!        Updating excitation energies
!
         call dcopy(wf%tasks%n_singlet_states, eigenvalues_Im_new, 1, eigenvalues_Im_old, 1)
         call dcopy(wf%tasks%n_singlet_states, eigenvalues_Re_new, 1, eigenvalues_Re_old, 1)
!
!        Get next trial vectors
!
         call wf%construct_next_trial_vectors(eigenvalues_Re_new, eigenvalues_Im_new, &
                                                solution_vectors_reduced, reduced_dim, n_new_trials)
!
!        Test for convergence
!
         if (converged_energy .and. converged_residual) then
            converged = .true.
         else
            iteration = iteration + 1
         endif
!
         call deallocator(solution_vectors_reduced, reduced_dim, wf%tasks%n_singlet_states)
!
      enddo
!
!     Prints
!
      if ( converged ) then
            write(unit_output,'(/t3,a,i2,a/)')  'Converged in ', iteration, ' iterations!'
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
   module subroutine solve_reduced_eigenvalue_equation_ccs(wf, eigenvalues_Re, eigenvalues_Im,&
                                                                solution_vectors_reduced, reduced_dim, n_new_trials)
!!
!!    Solve reduced eigenvalue problem
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad May 2017
!!
!!    Constructes the reduced A matrix, solves the reduced eigenvalue problem
!!    and returns the first n eigenvalues and eigenvectors.
!!
      implicit none
!
      class(ccs)                                                  :: wf
      integer(i15)                                                :: reduced_dim, n_new_trials
      real(dp), dimension(wf%tasks%n_singlet_states,1)            :: eigenvalues_Re
      real(dp), dimension(wf%tasks%n_singlet_states,1)            :: eigenvalues_Im 
      real(dp), dimension(reduced_dim, wf%tasks%n_singlet_states) :: solution_vectors_reduced
!
      real(dp), dimension(:,:), allocatable :: A_red
      real(dp), dimension(:,:), allocatable :: c_i
      real(dp), dimension(:,:), allocatable :: rho_j
      real(dp), dimension(:,:), allocatable :: eigenvectors  ! S: So maybe a new name for this? 
      real(dp), dimension(:,:), allocatable :: eigenvalues_full_Re
      real(dp), dimension(:,:), allocatable :: eigenvalues_full_Im
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
      integer      :: info = -1 
!
!
      call allocator(A_red, reduced_dim, reduced_dim)
      A_red = zero
!
!     -::- Prepare to solve the eigenvalue problem -::-
!
!     Prepare files
!  
      call generate_unit_identifier(unit_trial_vecs)
         open(unit=unit_trial_vecs, file='trial_vec', action='read', status='old', &
           access='direct', form='unformatted', recl=dp*(wf%n_parameters), iostat=ioerror)
!
      call generate_unit_identifier(unit_rho)
         open(unit=unit_rho, file='transformed_vec', action='read', status='old', &
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
!     Construct reduced jacobi matrix 
!     A_red_ij = c_i^T * A * c_j = c_i^T * rho_i
!
      if (iteration .eq. 1) then
         do i = 1,reduced_dim
           read(unit_trial_vecs, rec=i, iostat=ioerror) c_i
            do j = 1,reduced_dim
               read(unit_rho, rec=j, iostat=ioerror) rho_j
               A_red(i,j) = ddot(wf%n_parameters, c_i, 1, rho_j, 1)
            enddo
         enddo
      else
         do i = 1,reduced_dim
           read(unit_trial_vecs, rec=i, iostat=ioerror) c_i
            do j = reduced_dim - n_new_trials + 1, reduced_dim
             read(unit_rho, rec=j, iostat=ioerror) rho_j
               A_red(i,j) = ddot(wf%n_parameters, c_i, 1, rho_j, 1)
            enddo
         enddo
         do j = 1, reduced_dim - n_new_trials
           read(unit_rho, rec=j, iostat=ioerror) rho_j
           do i = reduced_dim - n_new_trials + 1, reduced_dim
              read(unit_trial_vecs, rec=i, iostat=ioerror) c_i
              A_red(i,j) = ddot(wf%n_parameters, c_i, 1, rho_j, 1)
            enddo
         enddo
      endif
!
      call deallocator(c_i, wf%n_parameters, 1)
      call deallocator(rho_j, wf%n_parameters, 1)
!
!
!     Close files for trial vectors and transformed vectors
!
      close(unit_trial_vecs)
      close(unit_rho)
!
!     Write new A_red to file
!
      rewind(unit_reduced_jacobi)
      write(unit_reduced_jacobi) ((A_red(i,j), i = 1, reduced_dim), j = 1, reduced_dim)
!
!     Close file
!
      close(unit_reduced_jacobi)
!
!     -::- Solve reduced eigenvalue problem -::-
!
!     Allocate arrays for eigenvalues and eigenvectors
!
      call allocator(eigenvectors, reduced_dim, reduced_dim)
      eigenvectors = zero
!
      call allocator(eigenvalues_full_Re, reduced_dim, 1)
      call allocator(eigenvalues_full_Im, reduced_dim, 1)
      eigenvalues_full_Re = zero
      eigenvalues_full_Im = zero
!
      call allocator(work, 4*reduced_dim, 1)
      work = zero
!
!     Solve eigenvalue problem
!
      call dgeev('N','V',              &
                  reduced_dim,         &
                  A_red,               &
                  reduced_dim,         &
                  eigenvalues_full_Re, &
                  eigenvalues_full_Im, &
                  dummy,               &
                  1,                   &
                  eigenvectors,        &
                  reduced_dim,         &
                  work,                &
                  4*reduced_dim,       &
                  info)
!
      call deallocator(work, 4*reduced_dim, 1)
!
!     Deallocate A_red
!
      call deallocator(A_red, reduced_dim, reduced_dim)
!
!     Find lowest eigenvalues and sort them
!
      call allocator_int(index_list,wf%tasks%n_singlet_states,1)
      index_list = 0
!
      call get_n_lowest(wf%tasks%n_singlet_states, reduced_dim, eigenvalues_full_Re, eigenvalues_Re, index_list)
!
!     Select solution vectors and imaginary parts of eigenvalues according to index_list 
!
      do i = 1, reduced_dim
         do j = 1, wf%tasks%n_singlet_states
!
            solution_vectors_reduced(i,j) = eigenvectors(i,index_list(j,1))
            eigenvalues_Im = eigenvalues_full_Im(index_list(j,1), 1)
!
         enddo
      enddo
!
!     Final deallocatons
!
      call deallocator(eigenvectors, reduced_dim, reduced_dim)
      call deallocator(eigenvalues_full_Im, reduced_dim, 1)
      call deallocator(eigenvalues_full_Re, reduced_dim, 1)
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
!!    Get next trial vectors.    
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
      real(dp), dimension(wf%tasks%n_singlet_states,1)            :: eigenvalues_Re
      real(dp), dimension(wf%tasks%n_singlet_states,1)            :: eigenvalues_Im
      real(dp), dimension(reduced_dim, wf%tasks%n_singlet_states) :: solution_vectors_reduced
      integer(i15) :: reduced_dim, n_new_trials
!
      integer(i15) :: root = 0, trial = 0, i = 0, j= 0
!
      real(dp), dimension(:,:), allocatable :: solution_vector
      real(dp), dimension(:,:), allocatable :: residual
      real(dp), dimension(:,:), allocatable :: orbital_diff
!
      real(dp), dimension(:,:), allocatable :: c_i
      real(dp), dimension(:,:), allocatable :: c_j
      real(dp), dimension(:,:), allocatable :: rho_i
!
      
      real(dp) :: norm_solution_vector = zero, norm_residual = zero, norm_new_trial = zero, dot_prod = zero  
!
      integer(i15) :: unit_trial_vecs = 0, unit_rho = 0, unit_solution = 0, ioerror = 0
!
      real(dp) :: ddot
!
!     Prepare necessary files
!
      call generate_unit_identifier(unit_trial_vecs)
      open(unit=unit_trial_vecs, file='trial_vec', action='readwrite', status='old', &
        access='direct', form='unformatted', recl=dp*(wf%n_parameters), iostat=ioerror)
!
      call generate_unit_identifier(unit_rho)
      open(unit=unit_rho, file='transformed_vec', action='read', status='old', &
      access='direct', form='unformatted', recl=dp*(wf%n_parameters), iostat=ioerror) 
!
      call generate_unit_identifier(unit_solution)
      open(unit=unit_solution, file='solution_vectors', action='write', status='unknown', &
      access='direct', form='unformatted', recl=dp*(wf%n_parameters), iostat=ioerror) 
!
      call allocator(residual, wf%n_parameters, 1)
!
      n_new_trials = 0
      converged_residual = .true.
!
      write(unit_output,'(/t3,a)') 'Root    Eigenvalue (Re)      Eigenvalue (Im)      Residual norm'
      write(unit_output,'(t3,a)') '----------------------------------------------------------------'
!
!     For each of the roots
!
      do root = 1, wf%tasks%n_singlet_states
         residual = zero
!
!       :: Create fullspace vector X and calculate norm ||X|| ::
!
!        Xj = sum_i x_j_i c_i
!
         call allocator(solution_vector, wf%n_parameters, 1)
         solution_vector = zero
!
         call allocator(c_i, wf%n_parameters, 1)
!
!        Xj = sum_i x_j_i * c_i
!
         do trial = 1, reduced_dim
!
            c_i = zero
            read(unit_trial_vecs, rec=trial, iostat=ioerror) c_i
            call daxpy(wf%n_parameters, solution_vectors_reduced(trial,root), c_i, 1, solution_vector, 1)
!
            if (ioerror .ne. 0) write(unit_output,*) 'Error reading trial vecs in get_next_trial_vectors'
!
         enddo
!
         call deallocator(c_i, wf%n_parameters, 1)
!
!        Calculate norm of solution vector
!
         norm_solution_vector = sqrt(ddot(wf%n_parameters, solution_vector, 1, solution_vector, 1))
!
!        :: Calculate residual ::
!
         call dcopy(wf%n_parameters, solution_vector, 1, residual, 1)
!
!        Normalize and write full space solution vectors to file
!
         call dscal(wf%n_parameters, one/norm_solution_vector, solution_vector, 1)
         write(unit_solution, rec=root, iostat=ioerror) solution_vector
!        
         call deallocator(solution_vector, wf%n_parameters, 1)
!
!        Residual = (AX - eX)
!
!        - eX : 
!
         call dscal(wf%n_parameters, -eigenvalues_Re(root, 1), residual, 1)
!         
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
            if (ioerror .ne. 0) write(unit_output,*) 'Error reading tranf vecs in get_next_trial_vectors'
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
         if (norm_residual/norm_solution_vector  .gt. wf%settings%ampeqs_threshold) then
            converged_residual = .false.
         endif
!
!        Prints
!
         write(unit_output,'(t3,i2,5x,f14.8,7x,f14.8,11x,e10.4)') root, eigenvalues_Re(root, 1), &
                                                                eigenvalues_Im(root, 1), norm_residual/norm_solution_vector
         flush(unit_output)
!
!        :: Preconditioning ::
!
         call allocator(orbital_diff, wf%n_parameters, 1)
         orbital_diff = zero
!
         call wf%calculate_orbital_differences(orbital_diff)
!
         do i = 1, wf%n_parameters
            residual(i, 1) = residual(i,1)/orbital_diff(i,1)
         enddo
!
         call deallocator(orbital_diff, wf%n_parameters, 1)
!
!        :: Orthogonalization ::
!
!        Calculate norm of residual
!
         norm_residual = sqrt(ddot(wf%n_parameters, residual, 1, residual, 1))
!
!        Normalize residual
!
         call dscal(wf%n_parameters, one/norm_residual, residual, 1)
!
!        Orthogonalize new trial vector against old trial vectors
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
         if ((norm_new_trial .gt. wf%settings%ampeqs_threshold) ) then
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
      write(unit_output,'(t3,a/)') '----------------------------------------------------------------'
!
      write(unit_output,*)'n_new_trials', n_new_trials
   end subroutine construct_next_trial_vectors_ccs
!
!
end submodule