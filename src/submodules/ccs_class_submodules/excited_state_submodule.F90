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
   integer(i15) :: iteration = 1 
!
   integer(i15) :: n_red        = 0 
   integer(i15) :: n_new_trials = 0
!
!
contains
!
!
   module subroutine excited_state_solver_ccs(wf)
!
!
      implicit none
!
      class(ccs) :: wf 
!
!     Let the user know the ground state solver is running
!
      write(unit_output,'(/t3,a)')   ':: Excited state solver (Davidson)'
      write(unit_output,'(t3,a/)')   ':: E. F. Kjønstad, S. D. Folkestad, May 2017'
      write(unit_output,'(t3,a,i3,a,a,a/)') &
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
      call wf%initialize_trial_vectors
!     call other initialization routines.
!
      n_red        = wf%tasks%n_singlet_states
      n_new_trials = wf%tasks%n_singlet_states
!
      call wf%transform_trial_vecs(n_red - n_new_trials + 1, n_red)
!
      call wf%solve_reduced_eigenvalue_problem
!
   end subroutine excited_state_solver_ccs
!
!
   module subroutine solve_reduced_eigenvalue_problem_ccs(wf)
!
      implicit none
!
      class(ccs) :: wf 
!
      real(dp), dimension(:,:), allocatable :: A_red
      real(dp), dimension(:,:), allocatable :: c_i
      real(dp), dimension(:,:), allocatable :: rho_j
      real(dp), dimension(:,:), allocatable :: eig_Re
      real(dp), dimension(:,:), allocatable :: eig_Im
      real(dp), dimension(:,:), allocatable :: eig_vec
      real(dp), dimension(:,:), allocatable :: work
      real(dp) :: ddot, dummy
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
!     :: Prepare to solve the eigenvalue problem ::
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
         open(unit=unit_reduced_jacobi, file='reduced_jacobi', action='readwrite', status='new',&
          form='unformatted', iostat=ioerror)
!
      else
!
         call generate_unit_identifier(unit_reduced_jacobi)
         open(unit=unit_reduced_jacobi, file='reduced_jacobi', action='readwrite', status='old',&
          form='unformatted', iostat=ioerror)
!
         rewind(unit_reduced_jacobi)
         read(unit_reduced_jacobi) ((A_red(i,j),i = 1, n_red-n_new_trials), j=1, n_red-n_new_trials)
      endif
!
!     Allocate c and rho
!
      call allocator(c_i, wf%n_parameters, 1)
      call allocator(rho_j, wf%n_parameters, 1)
      c_i   = zero
      rho_j = zero
!
!     Construct reduced jacobi matrix ! S: there is probably a smarter way of making A_red. Here new c are read twice
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
      write(unit_output,*) 'Reduced A:'
      do i = 1, n_red
            write(unit_output,*)(A_red(i,j), j=1,n_red)
      enddo
      close(unit_trial_vecs)
      close(unit_rho)
!
!     Write new A_red to file
!
      rewind(unit_reduced_jacobi)
      write(unit_reduced_jacobi) A_red
!
!     Close file
!
      close(unit_reduced_jacobi)
!
!     :: Solve reduced eigenvalue problem ::
!
!     Allocate arrays for eigenvalues and eigenvectors
!
      call allocator(eig_Re, n_red, 1)
      call allocator(eig_Im, n_red, 1)
      eig_Re = zero
      eig_Im = zero
!
      call allocator(eig_vec, n_red, n_red)
      eig_vec = zero
!
      call allocator(work, 4*n_red, 1)
      work = zero
      call dgeev('N','V',  &
                  n_red,   &
                  A_red,   &
                  n_red,   &
                  eig_Re,  &
                  eig_Im,  &
                  dummy,   &
                  1,       &
                  eig_vec, &
                  n_red,   &
                  work,    &
                  4*n_red, &
                  info)
      call vec_print(eig_Re,n_red,1)
      call vec_print(eig_Im,n_red,1)
!
!     Deallocate A_red
!
      call deallocator(A_red, n_red, n_red)
      call deallocator(eig_Re, n_red, 1)
      call deallocator(eig_Im, n_red, 1)
      call deallocator(eig_vec, n_red, n_red)
!
   end subroutine solve_reduced_eigenvalue_problem_ccs
!
!
end submodule