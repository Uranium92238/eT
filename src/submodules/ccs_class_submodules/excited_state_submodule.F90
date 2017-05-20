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
   integer(i15) :: iteration = 1, max_iterations = 25
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
!!    Solver for excited states
!!
      implicit none
!
      class(ccs) :: wf
!
      integer(i15) :: n_red        = 0 
      integer(i15) :: n_new_trials = 0
!
!  Solution for the reduced eigenvalue problem
!
      real(dp), dimension(:,:), allocatable :: eigenvalues_Re_new
      real(dp), dimension(:,:), allocatable :: eigenvalues_Im_new
      real(dp), dimension(:,:), allocatable :: eigenvalues_Re_old
      real(dp), dimension(:,:), allocatable :: eigenvalues_Im_old
      real(dp), dimension(:,:), allocatable :: eigenvectors
!
      integer(i15) :: i = 0
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
!
!     Initialize for excited state calculation
!
      n_red        = wf%tasks%n_singlet_states
      n_new_trials = wf%tasks%n_singlet_states
!
      call wf%initialize_trial_vectors
!
!     Allocate eigenvalue arrays
!
      call allocator(eigenvalues_Im_old, wf%tasks%n_singlet_states, 1)
      call allocator(eigenvalues_Re_old, wf%tasks%n_singlet_states, 1)
      call allocator(eigenvalues_Im_new, wf%tasks%n_singlet_states, 1)
      call allocator(eigenvalues_Re_new, wf%tasks%n_singlet_states, 1)
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
         write(unit_output,'(t3,a,i3)')'Reduced space dimension:', n_red
         write(unit_output,'(t3,a/)')'----------------------------'
!
!        Transform new trial vectors
!
         call wf%transform_trial_vecs(n_red - n_new_trials + 1, n_red)
!
!        Allocate solution vectors for reduced problem
!
         call allocator(eigenvectors, n_red, wf%tasks%n_singlet_states)
         eigenvectors = zero
!
!        Solve the reduced eigenvalue problem
!
         call wf%solve_reduced_eigenvalue_problem(eigenvalues_Re_new, eigenvalues_Im_new, eigenvectors, n_red, n_new_trials)
!
!        Test energy convergence criteria
!
         converged_energy = .true.
!
         do i = 1, wf%tasks%n_singlet_states
            if ( abs(eigenvalues_Re_new(i,1)-eigenvalues_Re_old(i,1)) .gt. wf%settings%energy_threshold) then
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
         call wf%get_next_trial_vectors(eigenvalues_Re_new, eigenvalues_Im_new, eigenvectors, n_red, n_new_trials)
!
!        Test for convergence
!
         if (converged_energy .and. converged_residual) then
            converged = .true.
         else
            iteration = iteration + 1
         endif
!
         call deallocator(eigenvectors, n_red, wf%tasks%n_singlet_states)
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
      call deallocator(eigenvalues_Im_old, n_red, 1)
      call deallocator(eigenvalues_Re_old, n_red, 1)
      call deallocator(eigenvalues_Im_new, n_red, 1)
      call deallocator(eigenvalues_Re_new, n_red, 1)
!
   end subroutine excited_state_solver_ccs
!
!
   module subroutine solve_reduced_eigenvalue_problem_ccs(wf, eigenvalues_Re, eigenvalues_Im, solution_vectors, n_red, n_new_trials)
!!
!!    Solve reduced eigenvalue problem
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad May 2017
!!
!!    Constructes the reduced A matrix, solves the reduced eigenvalue problem
!!    and returns the first n eigenvalues and eigenvectors.
!!
      implicit none
!
      class(ccs)                                            :: wf
      integer(i15)                                          :: n_red, n_new_trials
      real(dp), dimension(wf%tasks%n_singlet_states,1)      :: eigenvalues_Re
      real(dp), dimension(wf%tasks%n_singlet_states,1)      :: eigenvalues_Im 
      real(dp), dimension(n_red, wf%tasks%n_singlet_states) :: solution_vectors
!
      real(dp), dimension(:,:), allocatable :: A_red
      real(dp), dimension(:,:), allocatable :: c_i
      real(dp), dimension(:,:), allocatable :: rho_j
      real(dp), dimension(:,:), allocatable :: eigenvectors
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
      call allocator(A_red, n_red, n_red)
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
          form='formatted', iostat=ioerror)
!
      else
!
         call generate_unit_identifier(unit_reduced_jacobi)
         open(unit=unit_reduced_jacobi, file='reduced_jacobi', action='readwrite', status='old',&
          form='formatted', iostat=ioerror)
!
         rewind(unit_reduced_jacobi)
         read(unit_reduced_jacobi,*) ((A_red(i,j),i = 1, n_red-n_new_trials), j=1, n_red-n_new_trials)
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
!
      if (iteration .eq. 1) then
         do i = 1,n_red
           read(unit_trial_vecs, rec=i, iostat=ioerror) c_i
            do j = 1,n_red
               read(unit_rho, rec=j, iostat=ioerror) rho_j
               A_red(i,j) = ddot(wf%n_parameters, c_i, 1, rho_j, 1)
            enddo
         enddo
      else
         do i = 1,n_red
           read(unit_trial_vecs, rec=i, iostat=ioerror) c_i
            do j = n_red - n_new_trials + 1, n_red
               read(unit_rho, rec=j, iostat=ioerror) rho_j
               A_red(i,j) = ddot(wf%n_parameters, c_i, 1, rho_j, 1)
            enddo
         enddo
         do j = 1, n_red - n_new_trials
           read(unit_rho, rec=j, iostat=ioerror) rho_j
           do i = n_red - n_new_trials + 1, n_red
              read(unit_trial_vecs, rec=i, iostat=ioerror) c_i
              A_red(i,j) = ddot(wf%n_parameters, c_i, 1, rho_j, 1)
            enddo
         enddo
      endif
!
      call deallocator(c_i, wf%n_parameters, 1)
      call deallocator(rho_j, wf%n_parameters, 1)
!
!      write(unit_output,*) 'Reduced A:'
!      do i = 1, n_red
!            write(unit_output,*)(A_red(i,j), j=1,n_red)
!      enddo
!
!     Close files for trial vectors and transformed vectors
!
      close(unit_trial_vecs)
      close(unit_rho)
!
!     Write new A_red to file
!
      rewind(unit_reduced_jacobi)
      write(unit_reduced_jacobi,*) ((A_red(i,j), i = 1, n_red), j = 1, n_red)
!
!     Close file
!
      close(unit_reduced_jacobi)
!
!     -::- Solve reduced eigenvalue problem -::-
!
!     Allocate arrays for eigenvalues and eigenvectors
!
      call allocator(eigenvectors, n_red, n_red)
      eigenvectors = zero
!
      call allocator(eigenvalues_full_Re, n_red, 1)
      call allocator(eigenvalues_full_Im, n_red, 1)
      eigenvalues_full_Re = zero
      eigenvalues_full_Im = zero
!
      call allocator(work, 4*n_red, 1)
      work = zero
!
!     Solve eigenvalue problem
!
      call dgeev('N','V',              &
                  n_red,               &
                  A_red,               &
                  n_red,               &
                  eigenvalues_full_Re, &
                  eigenvalues_full_Im, &
                  dummy,               &
                  1,                   &
                  eigenvectors,        &
                  n_red,               &
                  work,                &
                  4*n_red,             &
                  info)
!
!     Deallocate A_red
!
      call deallocator(A_red, n_red, n_red)
!
!     Sort eigenvalues, and crop
!
      call allocator_int(index_list,wf%tasks%n_singlet_states,1)
      index_list = 0
      call get_n_lowest(wf%tasks%n_singlet_states, n_red, eigenvalues_full_Re, eigenvalues_Re, index_list)
!
!     Select solution vectors and imaginary parts of eigenvalues according to index_list 
!
      do i = 1, n_red
         do j = 1, wf%tasks%n_singlet_states
            solution_vectors(i,j) = eigenvectors(i,index_list(j,1))
            eigenvalues_Im = eigenvalues_full_Im(index_list(j,1), 1)
         enddo
      enddo
!
!     Write eigenvalues_Re
!
!
!      write(unit_output,*)'Eigenvalues'
!      write(unit_output,*) (eigenvalues_Re(i,1), i=1,wf%tasks%n_singlet_states)
!
!      write(unit_output,*) (eigenvalues_full_Re(i,1), i=1,n_red)
!
!     Final deallocatons
!
      call deallocator(eigenvectors, n_red, n_red)
      call deallocator(eigenvalues_full_Im, n_red, 1)
      call deallocator(eigenvalues_full_Re, n_red, 1)
      call deallocator_int(index_list,wf%tasks%n_singlet_states,1)
!
   end subroutine solve_reduced_eigenvalue_problem_ccs
!
!
   module subroutine get_next_trial_vectors_ccs(wf, eigenvalues_Re_new, eigenvalues_Im_new, &
                                                solution_vectors_red, &
                                                n_red, n_new_trials)
!!
!!
!!
         implicit none
!
         class(ccs) :: wf
         real(dp), dimension(wf%tasks%n_singlet_states,1) :: eigenvalues_Re_new
         real(dp), dimension(wf%tasks%n_singlet_states,1) :: eigenvalues_Im_new
         real(dp), dimension(n_red, wf%tasks%n_singlet_states) :: solution_vectors_red
         integer(i15) :: n_red, n_new_trials
!
         integer(i15) :: root = 0, trial = 0, i = 0
!
         real(dp), dimension(:,:), allocatable :: full_space_solution_vector
         real(dp), dimension(:,:), allocatable :: residual
         real(dp) :: norm_solution_vector  
         real(dp), dimension(:,:), allocatable :: residual_norms
!
         real(dp), dimension(:,:), allocatable :: c_i
         real(dp), dimension(:,:), allocatable :: rho_i
         real(dp), dimension(:,:), allocatable :: project_out
         real(dp), dimension(:,:), allocatable :: temp
!
         integer(i15) :: unit_trial_vecs = 0, unit_rho = 0, unit_solution = 0, ioerror = 0
         real(dp) :: ddot
!
!        Prepare necessary files
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
         call allocator(residual_norms, wf%tasks%n_singlet_states, 1)
         residual_norms = zero
!
         n_new_trials = 0
         converged_residual = .true.
!
         write(unit_output,'(/t3,a)') 'Root    Eigenvalue (Re)      Eigenvalue (Im)      Residual norm'
         write(unit_output,'(t3,a)') '----------------------------------------------------------------'
!
!        For each of the roots
!
         do root = 1, wf%tasks%n_singlet_states
!
!           Create fullspace vector X and calculate norm ||X||
!           Xj = sum_i x_j_i c_i
!
            call allocator(full_space_solution_vector, wf%n_parameters, 1)
            full_space_solution_vector = zero
!
            call allocator(c_i, wf%n_parameters, 1)
!
!           Xj = sum_i x_j_i * c_i
!
            do trial = 1, n_red
!
               c_i = zero
               read(unit_trial_vecs, rec=trial, iostat=ioerror) c_i
               call daxpy(wf%n_parameters, solution_vectors_red(trial,root), c_i, 1, full_space_solution_vector, 1)
!
            enddo
!
            call deallocator(c_i, wf%n_parameters, 1)
!
!           Calculate norm of solution vector
!
            norm_solution_vector = sqrt(ddot(wf%n_parameters,full_space_solution_vector, 1, full_space_solution_vector, 1))
!
!           Start to calculate residual
!
            call dcopy(wf%n_parameters, full_space_solution_vector, 1, residual, 1)
!
!           Normalize and write full space solution vectors to file
!
            call dscal(wf%n_parameters, one/norm_solution_vector, full_space_solution_vector, 1)
            write(unit_solution, rec=root, iostat=ioerror) full_space_solution_vector
            call deallocator(full_space_solution_vector, wf%n_parameters, 1)
!
!           Residual = (AX - eX)/||X||
!
!           -eX : 
!
            call dscal(wf%n_parameters, -eigenvalues_Re_new(root, 1), residual, 1)
!
!           AX : 
!
            call allocator(rho_i, wf%n_parameters, 1)
!
            do trial = 1, n_red
!
               rho_i = zero
               read(unit_rho, rec=trial, iostat=ioerror) rho_i
               call daxpy(wf%n_parameters, solution_vectors_red(trial,root), rho_i, 1, residual, 1)
!
            enddo
!
            call deallocator(rho_i, wf%n_parameters, 1)
            call dscal(wf%n_parameters, norm_solution_vector, residual, 1)
!
!           Calculate residual norm and check convergence criteria on residual norms
!
            residual_norms(root, 1) = sqrt(ddot(wf%n_parameters, residual, 1, residual, 1)) 
            if (residual_norms(root,1) .gt. wf%settings%ampeqs_threshold) then
               converged_residual = .false.
            endif
!
!           Prints
!
            write(unit_output,'(t3,i2,5x,f14.8,7x,f14.8,11x,e10.4)') root, eigenvalues_Re_new(root, 1), &
                                                                   eigenvalues_Im_new(root, 1), residual_norms(root, 1)
!
!           Orthogonalize new trial vector against old trial vectors
!
            call allocator(c_i, wf%n_parameters, 1)
            call allocator(project_out,wf%n_parameters,wf%n_parameters)
            call allocator(temp, wf%n_parameters, 1)
!
!           prod_i (I - c_i*c_i^T)*Residual
!
            do trial = 1, n_red
               c_i = zero
               read(unit_trial_vecs, rec=trial, iostat=ioerror) c_i
!
!             Y = -c_i*c_i^T : 
!
               call dgemm('N', 'T', &
                           wf%n_parameters, &
                           wf%n_parameters, &
                           1,               &
                           -one,            &
                           c_i,             &
                           wf%n_parameters, &
                           c_i,             &
                           wf%n_parameters, &
                           zero,            &
                           project_out,     &
                           wf%n_parameters)
!
!              Y = I + Y : 
!
               do i = 1, wf%n_parameters
                  project_out(i,i) = project_out(i,i) + 1
               enddo
!
!              Residual*Y : 
!
               call dgemm('N','N', &
                  wf%n_parameters, &
                  1,               &
                  wf%n_parameters, &
                  one,             &
                  project_out,     &
                  wf%n_parameters, &
                  residual,        &
                  wf%n_parameters, &
                  zero,            &
                  temp,            &
                  wf%n_parameters  )
!
                  call dcopy(wf%n_parameters, temp, 1, residual, 1)
!
            enddo
!
            call deallocator(c_i, wf%n_parameters, 1)
            call deallocator(project_out,wf%n_parameters,wf%n_parameters)
            call deallocator(temp, wf%n_parameters, 1)
!
!           Calculate norm of new trial vector
!
            norm_solution_vector = sqrt(ddot(wf%n_parameters, residual, 1, residual, 1))
!
!           Test for linear dependency on old trial vectors
!           If norm sufficiently high new vector is normalized and written to file
!
            if (norm_solution_vector .gt. wf%settings%ampeqs_threshold) then
!
               n_new_trials = n_new_trials + 1
               call dscal(wf%n_parameters, one/norm_solution_vector, residual, 1)
               write(unit_trial_vecs, rec=n_new_trials+n_red, iostat=ioerror) residual
!
            endif
!
         enddo
!
!        Close all files
!
         close(unit_trial_vecs)
         close(unit_rho)
         close(unit_solution)
!
!        Update dimension of reduced space 
!
         n_red = n_red + n_new_trials
!
         call deallocator(residual, wf%n_parameters, 1)
         call deallocator(residual_norms, wf%tasks%n_singlet_states, 1)
!
         write(unit_output,'(t3,a/)') '----------------------------------------------------------------'
      end subroutine get_next_trial_vectors_ccs
!
!
end submodule