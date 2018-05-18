submodule (ccs_class) ground_state 
!
!!
!!    Ground state submodule (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!    Consists of the following module subroutines of the CCS:
!!
!!    ground_state_solver:        Controls the iterative loop, calling in turn
!!                                the calculation of the energy, the amplitude equations 
!!                                (and its norm), and the new_amplitudes routine.
!!
!!    new_amplitudes:             Calculates the quasi-Newton estimate and passes the 
!!                                information needed by the DIIS routine.
!!
!!    diis_ccs:                   This routine saves the quasi-Newton estimate Δ t and
!!                                t + Δ t to file. It uses the previous estimates to
!!                                select the amplitudes t for the next iteration.
!!    
!!
!!    calc_ampeqs:                Updates the amplitude equations for the current amplitudes.
!!    calc_ampeqs_norm:           Calculates the norm of the amplitude equations.
!!    calc_quasi_Newton_singles:  Calculates the singles part of the quasi-Newton estimate.
!!
!!    Can be inherited by models of the same level (e.g. CC2) without modification.
!!
!!    When inherited by higher level models (e.g. CCSD), the new_amplitudes and calc_ampeqs_norm
!!    routines should be overridden to account for the doubles quasi-Newton estimate, amplitudes, 
!!    and projection vector.
!!
!
   implicit none
!
!  Some variables available to all routines of the module
!
   integer(i15) :: iteration = -1 
!
   integer(i15) :: unit_dt          = -1   ! Unit identifier for Δ t_i file 
   integer(i15) :: unit_t_dt        = -1   ! Unit identifier for t_i + Δ t_i file
   integer(i15) :: unit_diis_matrix = -1   ! Unit identifier for DIIS matrix file
!
   integer(i15), parameter :: diis_dim = 9 ! The maximum dimension of the DIIS matrix, plus 1 
!
!
contains
!
!
   module subroutine ground_state_driver_ccs(wf)
!!
!!    Ground state driver (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
!!    Directs the solution of the ground state problem for CCS. The
!!    routine is written so as to be inherited unaltered in the CC hierarchy. 
!!
      implicit none 
!
      class(ccs) :: wf  
!
!     Let the user know the ground state solver is running
!
      write(unit_output,'(/t3,a)')   ':: Ground state solver (DIIS)'
      write(unit_output,'(t3,a/)')   ':: S. D. Folkestad, E. F. Kjønstad, May 2017'
!
      write(unit_output,'(t3,a,a,a/)')  'Settings for ',trim(wf%name), ' ground state calculation:'
!
      write(unit_output,'(t6,a20,e9.2)') 'Energy threshold:   ', wf%ground_state_specifications%energy_threshold
      write(unit_output,'(t6,a20,e9.2)') 'Residual threshold: ', wf%ground_state_specifications%residual_threshold
      flush(unit_output)
!
!     Preparations for ground state solver 
!
      call wf%ground_state_preparations
!
!     Run the solver routine
!
      call wf%ground_state_solver
!
!     Final work and preparations for other tasks (such as excited state calculations)
!
      call wf%ground_state_cleanup
!
   end subroutine ground_state_driver_ccs
!
!
   module subroutine ground_state_preparations_ccs(wf)
!!
!!    Ground State Preparations (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
!!    A routine for preparation tasks (if any). Can be overwritten
!!    in descendants if other preparations prove necessary.    
!!
      class(ccs) :: wf 
!
      wf%tasks%current = 'ground_state'
!
!     Allocate amplitudes (if not allocated) and calculate number of amplitudes 
!
      call wf%initialize_amplitudes 
!
!     Allocate projection vector 
!
      call wf%initialize_omega   
!
   end subroutine ground_state_preparations_ccs
!
!
   module subroutine ground_state_cleanup_ccs(wf)
!!
!!    Ground State Cleanup (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
!!    A routine for cleanup tasks (if any). Can be overwritten
!!    in descendants if other cleanups prove necessary.    
!!
      class(ccs) :: wf 
!
      call wf%destruct_amplitudes
      call wf%destruct_omega
!
   end subroutine ground_state_cleanup_ccs
!
!
   module subroutine ground_state_solver_ccs(wf)
!!
!!    Ground State Solver 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!    Directs the solution of the ground state amplitude equations
!!    using a DIIS algorithm. The problem the routine solves is 
!!
!!       X_mu(t) = 0, where t = { t_mu }_mu 
!!
!!    For standard coupled cluster theories, the vector X is the
!!    projection vector (omega).
!!
      implicit none
!
      class(ccs) :: wf 
!
      type(diis) :: ground_state_diis 
!
      real(dp) :: prev_energy
      real(dp) :: ampeqs_norm
!
      real(dp) :: start_gs_solver, end_gs_solver
!
      logical :: converged_energy = .false.
      logical :: converged_ampeqs = .false.
!
      logical :: converged = .false. ! True iff both the energy and the equations have converged 
!
!     If restart, read amplitudes from disk 
!
      if (wf%ground_state_specifications%restart) then 
!
         write(unit_output,'(/t3,a)') 'Requested restart. Preparing for restart.'
         call wf%ground_state_restart
!
      endif
!
!     Initialize DIIS object  
!
      call ground_state_diis%init('ground_state', wf%n_parameters)
!
!     Open DIIS files 
!
!       call generate_unit_identifier(unit_dt)
!       open(unit=unit_dt,file='diis_dt',status='unknown',form='unformatted')
! !
!       call generate_unit_identifier(unit_t_dt)
!       open(unit=unit_t_dt,file='diis_t_dt',status='unknown',form='unformatted')
! !
!       call generate_unit_identifier(unit_diis_matrix)
!       open(unit=unit_diis_matrix,file='diis_matrix',status='unknown',form='unformatted')
!
!     Enter iterative loop
!
      write(unit_output,'(/t3,a)')   'Iter.    Energy               Norm of amplitude eq.'
      write(unit_output,'(t3,a)')    '---------------------------------------------------' 
      flush(unit_output)
!
!     Make sure the initial energy is up to date for first iteration
!  
      call wf%calc_energy
!
      iteration = 1
      converged_energy = .false.
      converged_ampeqs = .false.
      converged = .false.
!
      call cpu_time(start_gs_solver)
!
      do while ((.not. converged) .and. (iteration .le. wf%ground_state_specifications%max_iterations))
!
!        Save the previous energy 
!
         prev_energy = wf%energy 
!
!        Update the energy 
!
         call wf%calc_energy
!
!        Update the Fock matrix 
!
         call wf%construct_fock
!
!        Construct the current amplitude equations vector,
!        and calculate the norm of the amplitude equations
!
         call wf%calc_ampeqs
         call wf%calc_ampeqs_norm(ampeqs_norm)
!
!        Check for convergence of the energy and the amplitude equations
!
         converged_energy = abs(wf%energy-prev_energy) .lt. wf%ground_state_specifications%energy_threshold
         converged_ampeqs = ampeqs_norm                .lt. wf%ground_state_specifications%residual_threshold
!
!        Print information to output 
!
         write(unit_output,'(t3,i3,4x,f19.12,4x,e10.4)') iteration, wf%energy, ampeqs_norm 
         flush(unit_output) ! Flush so that the user can follow each iteration in real-time
!
!        Perform DIIS update if convergence hasn't been reached
!
         if (converged_energy .and. converged_ampeqs) then
!
            converged = .true.
!
            write(unit_output,'(t3,a)')    '---------------------------------------------------' 
            if (iteration .eq. 1 .and. wf%name .ne. 'CCS') write(unit_output,'(/t3,a,/t3,a)') &
                                                                  'Note: residual converged in first iteration.', &
                                                                    'Energy convergence therefore not tested in this calculation.'
!
            write(unit_output,'(/t3,a,i3,a/)')  'Converged in ', iteration, ' iterations!'
!
         else
!
!           Update the amplitudes for the next iteration 
!       
        !    call wf%new_amplitudes
            call wf%new_amplitudes(ground_state_diis)
            iteration = iteration + 1
!
         endif
!
      enddo
!
!     Close the DIIS files
!
      ! close(unit_dt)
      ! close(unit_t_dt)
      ! close(unit_diis_matrix)
!
      call cpu_time(end_gs_solver)
!
!     Print summary
!
      write(unit_output,'(/t3,a,a,a/)')'Summary of ', trim(wf%name), ' ground state calculation:'
      write(unit_output,'(t6,a25,f19.12)')  'Total energy (hartrees):  ', wf%energy
      write(unit_output,'(t6,a25,f19.12/)') 'Total time CPU (seconds): ', end_gs_solver - start_gs_solver
      flush(unit_output)
!
!     Save the amplitudes 
!
       call wf%save_amplitudes
!
!     Issue an error & stop if the equations did not converge 
!
      if (.not. converged) then 
!
         write(unit_output,'(/t3,a)') 'Error: Ground state equations did not converge.'
         write(unit_output,'(t3,a/)') 'Consider increasing the maximum number of iterations.'
         stop
!
      endif
!
   end subroutine ground_state_solver_ccs
!
!
   module subroutine ground_state_restart_ccs(wf)
!!
!!    Ground state restart (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Dec 2017
!!
!!    Called if restart of the ground state is requested.   
!!
      class(ccs) :: wf 
!
      call wf%read_amplitudes   
!
   end subroutine ground_state_restart_ccs
!
!
   module subroutine calc_ampeqs_ccs(wf)
!!
!!    Calculate amplitude e quations (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!    Constructs the amplitude equations vector (the projection vector 
!!    in CCS) for the amplitudes of the current iteration of the ground 
!!    state solver. It also calculates the norm of the amplitude equations, 
!!    which is zero when the equations are exactly solved.
!!
      implicit none 
!
      class(ccs) :: wf 
!
!     Update the projection vector 
!
      call wf%construct_omega
!
   end subroutine calc_ampeqs_ccs
!
!
    module subroutine calc_ampeqs_norm_ccs(wf, ampeqs_norm)
!!
!!     Calculate Amplitude Equations Norm (CCS)
!!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
      implicit none 
!
      class(ccs) :: wf 
!
      real(dp) :: ampeqs_norm 
!
      real(dp) :: ddot ! For dot product
!
      ampeqs_norm = zero
      ampeqs_norm = ddot(wf%n_t1am, wf%omega1, 1, wf%omega1, 1)
      ampeqs_norm = sqrt(ampeqs_norm)
!
   end subroutine calc_ampeqs_norm_ccs
!
!
    module subroutine new_amplitudes_ccs(wf, diis_ground_state)
!!
!!    New Amplitudes (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!    Directs the calculation of the quasi-Newton estimate Δ t_i, 
!!    and t_i + Δ t_i, and calls the DIIS routine to save & get 
!!    the amplitudes for the next iteration.
!!
      implicit none 
!
      class(ccs) :: wf 
!
      class(diis) :: diis_ground_state
!
      real(dp), dimension(:,:), allocatable :: dt   ! Δ t_i
      real(dp), dimension(:,:), allocatable :: t_dt ! t_i + Δ t_i
!
!     Allocate Δ t_i and t_i + Δ t_i vectors 
! 
      call wf%mem%alloc(dt, wf%n_parameters, 1)
      call wf%mem%alloc(t_dt, wf%n_parameters, 1)
!
      dt   = zero 
      t_dt = zero 
!
!     Calculate Δ t_i
!
      call wf%calc_quasi_Newton_singles(dt)
!
!     Set t_i + Δ t_i 
!
      call dcopy(wf%n_t1am, dt, 1, t_dt, 1)           ! t_dt = Δ t_i 
      call daxpy(wf%n_t1am, one, wf%t1am, 1, t_dt, 1) ! t_dt = t_i + Δ t_i
!
!     Save estimates to file and get the next amplitudes
!     (they are placed in dt on exit from diis) 
!
      call wf%diis(dt, t_dt)
!
!     Set the new amplitudes 
!
      call dcopy(wf%n_t1am, dt, 1, wf%t1am, 1)
!
!     Deallocate vectors 
!
      call wf%mem%dealloc(dt, wf%n_parameters, 1)
      call wf%mem%dealloc(t_dt, wf%n_parameters, 1)
!
   end subroutine new_amplitudes_ccs
!
!
    module subroutine calc_quasi_Newton_singles_ccs(wf,dt)
!!
!!    Calculate quasi-Newton estimate (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!    Calculates the quasi-Newton estimate Δ t_i (singles part)
!!    and places the contribution in the dt vector (of length n_parameters,
!!    with singles first, then doubles, etc. if inherited)
!!
      implicit none 
!
      class(ccs) :: wf 
!
      real(dp), dimension(wf%n_parameters, 1) :: dt
!
      integer(i15) :: a = 0, i = 0, ai = 0
!
!     Calculate the singles Δ t_i contribbution
!
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            ai = index_two(a, i, wf%n_v)
!
            dt(ai, 1) = - wf%omega1(a, i)/(wf%fock_diagonal(wf%n_o + a, 1) - &
                                           wf%fock_diagonal(i, 1))
!
         enddo
      enddo
!
   end subroutine calc_quasi_Newton_singles_ccs
!
!
    module subroutine diis_ccs(wf, dt, t_dt)
!!
!!    DIIS routine
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!    The next amplitudes are 
!!
!!       t_n+1 = sum_k w_k (t_k + dt_k), 
!! 
!!    where the weights w_k in front of the quasi-Newton estimate dt_k
!!    are determined so as to minimize 
!!
!!       f(w_k) = sum_k w_k dt_k, 
!!
!!    with the constraint that g(w_k) = sum_k w_k - 1 = 0.
!!
      implicit none 
!
      class(ccs), intent(in) :: wf 
!
      real(dp), dimension(wf%n_parameters, 1) :: dt 
      real(dp), dimension(wf%n_parameters, 1) :: t_dt 
!
      real(dp), dimension(:,:), allocatable :: dt_i ! To hold previous Δ t_i temporarily
!
      real(dp) :: ddot
!
      integer(i15) :: i = 0, j = 0
!
      integer      :: info = -1         ! Error integer for dgesv routine (LU factorization)
      integer(i15) :: current_index = 0 ! Progressing as follows: 1,2,...,7,8,1,2,...
!
      real(dp), dimension(:,:), allocatable :: diis_vector
      real(dp), dimension(:,:), allocatable :: diis_matrix
!
      integer(i15), dimension(diis_dim) :: ipiv = 0 ! Pivot integers (see dgesv routine)
!
!     Set the current index 
!
      current_index = iteration - (diis_dim-1)*((iteration-1)/(diis_dim-1)) ! 1,2,...,7,8,1,2,...
!
!     :: Save (Δ t_i) and (t_i + Δ t_i) to files ::
!
      if (current_index .eq. 1) then  
         rewind(unit_dt)
         rewind(unit_t_dt)
      endif
!
      write(unit_dt)   (dt(i,1), i = 1, wf%n_parameters)
      write(unit_t_dt) (t_dt(i,1), i = 1, wf%n_parameters)
!
!     :: Solve the least squares problem, G * w = H ::
!
!        G : DIIS matrix, G_ij = Δ t_i Δ t_j,
!        H : DIIS vector,  H_i = 0,
!
!     where i, j = 1, 2, ..., current_index. To enforce normality
!     of the solution, G is extended with a row & column of -1's 
!     and H with a -1 at the end.
!
!     First set the DIIS vector to one 
!
      call wf%mem%alloc(diis_vector,current_index+1,1)
      diis_vector = zero 
!
!     Allocate the DIIS matrix and read in previous matrix elements
!
      call wf%mem%alloc(diis_matrix, current_index+1, current_index+1)
      diis_matrix = zero 
!
      if (current_index .gt. 1) then 
!
         rewind(unit_diis_matrix)
!
         do j = 1, current_index - 1
            do i = 1, current_index - 1
!
               read(unit_diis_matrix) diis_matrix(i,j)
!
            enddo
         enddo
!
      endif
!
!     Get the parts of the DIIS matrix G not constructed in 
!     the previous iterations 
!
      call wf%mem%alloc(dt_i, wf%n_parameters, 1) ! Allocate temporary holder of quasi-Newton estimates
      dt_i = zero 
!
      rewind(unit_dt)
!
      do i = 1, current_index
!
         read(unit_dt) (dt_i(j,1), j = 1, wf%n_parameters) 
!
         diis_matrix(current_index,i) = ddot(wf%n_parameters, dt, 1, dt_i, 1) 
         diis_matrix(i,current_index) = diis_matrix(current_index,i)
!
         diis_matrix(current_index+1,i) = -one
         diis_matrix(i,current_index+1) = -one 
!
      enddo
!
      diis_vector(current_index+1,1) = -one
!
!     Write the current DIIS matrix to file 
!
      rewind(unit_diis_matrix)
!
      do j = 1, current_index
         do i = 1, current_index
!
            write(unit_diis_matrix) diis_matrix(i,j)
!
         enddo
      enddo
!
!     Solve the DIIS equation 
!
!     Note: on exit, the solution is in the diis_vector,
!     provided info = 0 (see LAPACK documentation for more)
!
      call dgesv(current_index+1,  &
                  1,               &
                  diis_matrix,     &
                  current_index+1, &
                  ipiv,            &
                  diis_vector,     &
                  current_index+1, &
                  info)
!
!     :: Update the amplitudes (placed in dt on exit) ::
!
      dt = zero
!
      rewind(unit_t_dt)
!
      do i = 1, current_index
!
!        Read the t_i + Δ t_i vector 
!
         t_dt = zero
         read(unit_t_dt) (t_dt(j, 1), j = 1, wf%n_parameters)
!
!        Add w_i (t_i + Δ t_i) to the amplitudes 
!
         call daxpy(wf%n_parameters, diis_vector(i, 1), t_dt, 1, dt, 1)
!
      enddo
!
!     Deallocations 
!
      call wf%mem%dealloc(dt_i, wf%n_parameters, 1)
      call wf%mem%dealloc(diis_vector, current_index + 1, 1)
      call wf%mem%dealloc(diis_matrix, current_index + 1, current_index+1)
!
   end subroutine diis_ccs
!
!
end submodule ground_state
