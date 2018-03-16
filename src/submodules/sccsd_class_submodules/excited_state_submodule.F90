submodule (sccsd_class) excited_state
!
!!
!!    Excited state submodule (SCCSD)
!!    Written by Eirik F. Kjønstad, June 2017
!!
!!    Contains the following subroutines of the SCCSD class:
!!
!!    excited_state_driver: controls the solution of the coupled set of equations 
!!                          for the ground and excited states of SCCSD. 
!!    
!
   implicit none 
!
   integer(i15) :: diis_dim = 7
!
contains
!
!
   module subroutine excited_state_driver_sccsd(wf)
!!
!!    Excited state driver (SCCSD)
!!    Written by Eirik F. Kjønstad, June 2017
!!
!!    Directs the solution of the excited state problem for SCCSD. 
!!
      implicit none 
!
      class(sccsd) :: wf 
!
      logical :: converged = .false.
!
      integer(i15) :: iteration = 1
      integer(i15) :: max_iterations = 100
   !   integer(i15) :: max_iterations = 1 ! To run CCSD & compute overlap & then quit.
!
!     File unit identifiers & error integer 
!
      integer(i15) :: unit_triples      = 0 ! Stores converged triple amplitude for restart
      integer(i15) :: unit_sccsd_cycles = 0 ! Stores information, for prints, of iterative loop
!
      integer(i15) :: ioerror = 0
!
      integer(i15) :: I = 0 ! Over iteration loops - for printing 
!
!     Damping parameters for DIIS solver 
!
      real(dp)            :: diis_damper = 0.80D0
      real(dp), parameter :: diis_damper_no_history = 0.5D0
!
!     Timings
!
      real(dp) :: begin_timer, end_timer
!
      real(dp) :: t_dt_local, dt_local
!
!     Holding triple amplitude & overlap during printout at the end 
!
      real(dp) :: temp_triples, temp_overlap
!
!     Let the user know the excited state driver is running,
!     and print some basic information 
!
      call cpu_time(begin_timer)
!
      write(unit_output,'(t3,a)') ':: Ground and excited state SCCSD driver'
      write(unit_output,'(t3,a)') ':: E. F. Kjønstad, June 2017'
!
      write(unit_output,'(/t3,a)') 'The driver will solve the ground and excited state'
      write(unit_output,'(t3,a)')  'equations in a series of cycles. The triple amplitude'
      write(unit_output,'(t3,a)')  'is altered in each cycle until the generalized overlap'
      write(unit_output,'(t3,a/)') 'between the constrained states is zero.'
!
      write(unit_output,'(t3,a/)') 'A single cycle consists of four steps:'
!
      write(unit_output,'(t6,a)')     '  I. Solve ground state equations'
      write(unit_output,'(t6,a)')     ' II. Solve multiplier equation*'
      write(unit_output,'(t6,a)')     'III. Solve excited state eigenvalue problem'
      write(unit_output,'(t6,a/)')    ' IV. Calculate overlap'
!
      write(unit_output,'(t6,a)')     '* This step is (currently) omitted for ground state intersections,'
      write(unit_output,'(t6,a/)')    '  and will in the future also be omitted for excited states.'
!
      write(unit_output,'(t3,a)')  'A DIIS-algorithm updates the triple amplitude to make'
      write(unit_output,'(t3,a/)') 'the generalized overlap vanish.'
!
      write(unit_output,'(t3,a,i2,a,a,a)') &
               'Requested ',wf%excited_state_specifications%n_singlet_states,' ', trim(wf%name), ' singlet states.'
      write(unit_output,'(t3,a,i2,a,a,a)') &
               'Requested ',wf%excited_state_specifications%n_triplet_states,' ', trim(wf%name), ' triplet states.' 
!
!     Read the triples amplitude (if it exists)
!
      call wf%read_triples
!
!     Open file for storing cycle convergence information 
! 
      call generate_unit_identifier(unit_sccsd_cycles)
      open(unit=unit_sccsd_cycles, file='sccsd_cycles', action='readwrite', status='unknown', &
           access='direct', form='unformatted', recl=dp, iostat=ioerror)
!
!     Enter DIIS-self-consistent SCC cycles
!
      do while (.not. converged .and. iteration .le. max_iterations)
!
         write(unit_output,'(/t3,a29,i3)')     'Similarity constrained cycle:', iteration 
         write(unit_output,'(t3,a53,f16.12)')  'In this cycle, we use the following triple amplitude:', wf%triples
!
         if (wf%state_A .eq. 0) then ! Ground state intersection
!
            call wf%ground_state_intersection_cycle(iteration)
!
         else
!
            call wf%excited_state_intersection_cycle(iteration)
!
         endif
!
         write(unit_output,'(/t3,a)') 'Summary of similarity constrained cycle:' 
!
         write(unit_output,'(/t6,a20,f16.12)') 'Triple amplitude:   ', wf%triples
         write(unit_output,'(t6,a20,f16.12)')  'Generalized overlap:', wf%overlap
!
!        Store to cycles file (triples, then overlap)
!
         write(unit_sccsd_cycles, rec=2*iteration-1, iostat=ioerror) wf%triples
         write(unit_sccsd_cycles, rec=2*iteration, iostat=ioerror)   wf%overlap
!
!        Check for convergence of the overlap 
!
         if (abs(wf%overlap) .le. wf%scc_settings%overlap_threshold) then 
!
            converged = .true.
            write(unit_output,'(/t3,a,i2,a/)')  'The similarity constrained equations converged in ', iteration, ' cycles!'
!
            call generate_unit_identifier(unit_triples)
            open(unit_triples, file='triples', status='unknown', form='unformatted')
!
            rewind(unit_triples)
            write(unit_triples) wf%triples 
            close(unit_triples)
!
         else
!
!           Update amplitude
!
            dt_local   = wf%overlap 
            t_dt_local = wf%overlap + wf%triples 
!
            call wf%sccsd_diis(dt_local, t_dt_local, iteration)
!
!           Update amplitudes according to DIIS, with dampings
!
            if (iteration .eq. 0) then
!
               wf%triples = diis_damper_no_history*dt_local + (one-diis_damper_no_history)*wf%triples
!
            else
!
               wf%triples = diis_damper*dt_local + (one-diis_damper)*wf%triples
!
            endif
!
            if (iteration .eq. 1) then
!
!              Set restarts to true for subsequent iterations 
!
               wf%ground_state_specifications%restart  = .true.
               wf%excited_state_specifications%restart = .true.
               wf%response_specifications%restart      = .true.
!
            endif
!
            iteration = iteration + 1
!
         endif
!
      enddo
!
      if (.not. converged) then ! Let the user know before quitting
!
         write(unit_output,'(/t3,a)') 'Error: reached maximum number of cycles without convergence.'
         stop
!
      else ! Display the final energies, etc. 
!
!        Display cycle information
!
         write(unit_output,'(/t6,a5,4x,a16,5x,a19)')'Cycle', 'triple amplitude', 'generalized overlap'
         write(unit_output,'(t6,a)')'--------------------------------------------------'
!
         do i = 1, iteration
!
            read(unit_sccsd_cycles, rec=2*i-1, iostat=ioerror) temp_triples
            read(unit_sccsd_cycles, rec=2*i, iostat=ioerror) temp_overlap
!
            write(unit_output,'(t6,i3,6x,f15.12,6x,f15.12)') i, temp_triples, temp_overlap
!
         enddo
!
         write(unit_output,'(t6,a)')'--------------------------------------------------'
!
!        Display converged amplitudes and energies
!
         write(unit_output,'(/t3,a)') 'Summary of SCCSD ground and excited state calculation:'
!
         write(unit_output,'(/t6,a27,f19.12)') 'Triple amplitude:          ', wf%triples
         write(unit_output,'(t6,a27,f19.12/)') 'Ground state energy [a.u.]:', wf%energy
!
         write(unit_output,'(t6,a10,4x,a13,11x,a11,9x,a14)')'Excitation', 'energy [a.u.]', 'energy [eV]', 'energy [cm^-1]'
         write(unit_output,'(t6,a)')'---------------------------------------------------------------------------------'
!
         do I = 1, wf%excited_state_specifications%n_singlet_states
!
!           Print energy of excitation in eV, hartree and cm^-1
!
            write(unit_output,'(t6,i3,6x,f19.12,5x,f19.12,5x,f19.12)') i, wf%excitation_energies(i,1),        &
                                                                         wf%excitation_energies(i,1)*27.211399, &
                                                                         wf%excitation_energies(i,1)*219474.63
         enddo
!
         write(unit_output,'(t6,a)')'---------------------------------------------------------------------------------'
         write(unit_output,'(t6,a)') '1 a.u. = 27.211399 eV'
         write(unit_output,'(t6,a)') '1 a.u. = 219474.63 cm^-1'
!
      endif
!
      call cpu_time(end_timer)
      write(unit_output,'(/t6,a25,f8.2)') 'Total CPU time (seconds):', end_timer-begin_timer
!
      close(unit_sccsd_cycles)
!
   end subroutine excited_state_driver_sccsd
!
!
   module subroutine excited_state_intersection_cycle_sccsd(wf, iteration)
!!
!!    Excited state intersection cycle (SCCSD)
!!    Written by Eirik F. Kjønstad, June 2017
!!
!!    Performs all the tasks needed in a single cycle for an 
!!    excited state intersection 
!!
      implicit none 
!
      class(sccsd) :: wf 
!
      integer(i15) :: iteration 
!
!     Part I. Solve for the ground state 
!
      write(unit_output,'(/t3,a,i2)') 'Part I. Solving the ground state equation'
      wf%tasks%current = 'ground_state'
      call wf%ground_state_driver
!
!     Part II. Solve for the multipliers
!
      write(unit_output,'(t3,a,i2)') 'Part II. Solving the multiplier equation'
!
      call wf%excited_state_preparations
!
      write(unit_output,'(/t3,a)')  ':: Response solver (Davidson)'
      write(unit_output,'(t3,a)')   ':: E. F. Kjønstad, S. D. Folkestad, June 2017' 
      flush(unit_output)  
!
      wf%tasks%multipliers = .true.
      wf%tasks%current = 'multipliers'
      call wf%response_solver
      wf%tasks%multipliers = .false.
!
!     Part III. Solve for the right excited states 
!
      write(unit_output,'(t3,a)')  'Part III. Solving the eigenvalue equation for the right eigenvectors'
!
      write(unit_output,'(/t3,a)')  ':: Excited state solver (Davidson)'
      write(unit_output,'(t3,a)')   ':: E. F. Kjønstad, S. D. Folkestad, May 2017'
!
      wf%excited_state_specifications%right = .true.  ! Use A transformation
      wf%excited_state_specifications%left  = .false. ! Not A^T transformation
      wf%tasks%current = 'excited_state'
!
      call wf%excited_state_solver
      call wf%excited_state_cleanup
!
!     Part IV. Calculate the overlap between the similarity constrained states 
!
      write(unit_output,'(/t3,a)') 'Part IV. Computing the generalized overlap between the states'
!
!     Test whether one or both of the eigenvectors have flipped vis-a-vis
!     the last iteration 
!
      call wf%eigenvector_controller(iteration) 
!
!     Compute the generalized overlap 
!
      call wf%calc_overlap
!
   end subroutine excited_state_intersection_cycle_sccsd
!
!
   module subroutine ground_state_intersection_cycle_sccsd(wf, iteration)
!!
!!    Ground state intersection cycle (SCCSD)
!!    Written by Eirik F. Kjønstad, June 2017
!!
!!    Performs all the tasks needed in a single cycle for a 
!!    ground state intersection 
!!
      implicit none 
!
      class(sccsd) :: wf 
!
      integer(i15) :: iteration 
!
!     Part I. Solve for the ground state 
!
      write(unit_output,*) 'Memory in this iteration:', wf%mem%available
!
      write(unit_output,'(/t3,a,i2)') 'Part I. Solving the ground state equation'
      wf%tasks%current = 'ground_state'
      call wf%ground_state_driver
!
!     Part II. Solve for the right excited states 
!
      write(unit_output,'(t3,a)')  'Part III. Solving the eigenvalue equation for the right eigenvectors'
!
      write(unit_output,'(/t3,a)')  ':: Excited state solver (Davidson)'
      write(unit_output,'(t3,a)')   ':: E. F. Kjønstad, S. D. Folkestad, May 2017'
!
      call wf%excited_state_preparations
!
      write(unit_output,*) 'Memory:', wf%mem%available
!
      wf%excited_state_specifications%right = .true.  ! Use A transformation
      wf%excited_state_specifications%left  = .false. ! Not A^T transformation
      wf%tasks%current = 'excited_state'
!
      call wf%excited_state_solver
      call wf%excited_state_cleanup
!
!     Part IV. Calculate the overlap between the similarity constrained states 
!
      write(unit_output,'(/t3,a)') 'Part IV. Computing the generalized overlap between the states'
      flush(unit_output)
!
!     Test whether the excited eigenvector has flipped, vis-a-vis
!     the last iteration 
!
      call wf%ground_state_eigenvector_controller(iteration) 
!
!     Compute the generalized overlap 
!
      call wf%calc_ground_state_overlap
!
   end subroutine ground_state_intersection_cycle_sccsd
!
!
   module subroutine sccsd_diis_sccsd(wf, dt, t_dt, iteration)
!!
!!    SCCSD DIIS routine
!!    Written by Eirik F. Kjønstad, 2017
!!
!!    Adapted from the CCS DIIS routine by Sarai D. Folkestad and Eirik F. Kjønstad.
!!
!!    The next amplitude is 
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
!!    In the SCCSD loop, the quasi-Newton estimate is the value of the overlap. 
!!
      implicit none 
!
      class(sccsd) :: wf 
!
      integer(i15), intent(in) :: iteration
!
      real(dp) :: dt 
      real(dp) :: t_dt 
!
      real(dp) :: dt_i ! To hold previous Δ t_i temporarily
!
      real(dp) :: ddot
!
      integer(i15) :: unit_dt = 0, unit_t_dt = 0, unit_diis_matrix = 0
!
      integer(i15) :: i = 0, j = 0, ij = 0, ioerror = 0
!
      integer      :: info = -1         ! Error integer for dgesv routine (LU factorization)
      integer(i15) :: current_index = 0 ! Progressing as follows: 1,2,...,7,8,1,2,...
!
      real(dp), dimension(:,:), allocatable :: diis_vector
      real(dp), dimension(:,:), allocatable :: diis_matrix
!
      integer(i15), dimension(9) :: ipiv = 0 ! Pivot integers (see dgesv routine)
!
!     Set the current index 
!
      current_index = iteration - (diis_dim-1)*((iteration-1)/(diis_dim-1)) ! 1,2,...,7,8,1,2,...
!
!     Open DIIS files 
!
      call generate_unit_identifier(unit_dt)
      open(unit=unit_dt, file='sccsd_dt', action='readwrite', status='unknown', &
           access='direct', form='unformatted', recl=dp, iostat=ioerror) 
!
      call generate_unit_identifier(unit_t_dt)
      open(unit=unit_t_dt, file='sccsd_t_dt', action='readwrite', status='unknown', &
           access='direct', form='unformatted', recl=dp, iostat=ioerror) 
!
      call generate_unit_identifier(unit_diis_matrix)
      open(unit=unit_diis_matrix, file='sccsd_diis_matrix', action='readwrite', status='unknown', &
           access='direct', form='unformatted', recl=dp, iostat=ioerror)   
!
      if (ioerror .ne. 0) then 
!
         write(unit_output,'(t3,a)') 'Error: could not open SCCSD DIIS files'
         stop 
!
      endif
!
!     :: Save (Δ t_i) and (t_i + Δ t_i) to files ::
!
      write(unit_dt, rec=current_index, iostat=ioerror) dt
!
      if (ioerror .ne. 0) then 
!
         write(unit_output,'(t3,a)') 'Error: could not write to SCCSD dt vector'
         stop 
!
      endif
!
      write(unit_t_dt, rec=current_index, iostat=ioerror) t_dt
!
      if (ioerror .ne. 0) then 
!
         write(unit_output,'(t3,a)') 'Error: could not write to SCCSD t_dt vector'
         stop 
!
      endif
!
!     :: Solve the least squares problem, G * w = H ::
!
!        G : DIIS matrix,  G_ij = Δ t_i Δ t_j,
!        H : DIIS vector,  H_i = 0,
!
!     where i, j = 1, 2, ..., current_index. To enforce normality
!     of the solution, G is extended with a row & column of -1's 
!     and H with a -1 at the end.
!
!     First set the DIIS vector to one 
!
      call wf%mem%alloc(diis_vector, current_index+1, 1)
      diis_vector = zero 
!
!     Allocate the DIIS matrix and read in previous matrix elements
!
      call wf%mem%alloc(diis_matrix, current_index+1, current_index+1)
      diis_matrix = zero 
!
      if (current_index .gt. 1) then 
!
         do j = 1, current_index - 1
            do i = 1, current_index - 1
!
               ij = index_packed(i,j)
               read(unit_diis_matrix, rec=ij, iostat=ioerror) diis_matrix(i,j)
!
            enddo
         enddo
!
      endif
!
      if (ioerror .ne. 0) then 
!
         write(unit_output,'(t3,a)') 'Error: could not read SCCSD DIIS matrix'
         stop 
!
      endif
!
!     Get the parts of the DIIS matrix G not constructed in 
!     the previous iterations 
!
      dt_i = zero 
!
      do i = 1, current_index
!
         read(unit_dt, rec=i, iostat=ioerror) dt_i
!
         diis_matrix(current_index,i) = dt*dt_i 
         diis_matrix(i,current_index) = diis_matrix(current_index,i)
!
         diis_matrix(current_index+1,i) = -one
         diis_matrix(i,current_index+1) = -one 
!
      enddo
!
      if (ioerror .ne. 0) then 
!
         write(unit_output,'(t3,a)') 'Error: could not read SCCSD dt vector'
         stop 
!
      endif
!
      diis_vector(current_index+1,1) = -one
!
!     Write the current DIIS matrix to file 
!
      do j = 1, current_index
         do i = 1, current_index
!
            ij = index_packed(i,j)
            write(unit_diis_matrix, rec=ij, iostat=ioerror) diis_matrix(i,j)
!
         enddo
      enddo
!
      if (ioerror .ne. 0) then 
!
         write(unit_output,'(t3,a)') 'Error: could not write SCCSD DIIS matrix'
         stop 
!
      endif
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
      do i = 1, current_index
!
!        Read the t_i + Δ t_i vector 
!
         t_dt = zero
         read(unit_t_dt, rec=i, iostat=ioerror) t_dt
!
!        Add w_i (t_i + Δ t_i) to the amplitudes 
!
         dt = dt + diis_vector(i, 1)*t_dt
!
      enddo
!
      if (ioerror .ne. 0) then 
!
         write(unit_output,'(t3,a)') 'Error: could not read SCCSD t_dt vector'
         stop 
!
      endif
!
!     Deallocations 
!
      call wf%mem%dealloc(diis_vector, current_index + 1, 1)
      call wf%mem%dealloc(diis_matrix, current_index + 1, current_index+1)
!
!
!     Close files 
!
      close(unit_dt)
      close(unit_t_dt)
      close(unit_diis_matrix)
!
   end subroutine sccsd_diis_sccsd
!
!
   module subroutine eigenvector_controller_sccsd(wf, iteration)
!!
!!    Eigenvector controller (SCCSD)
!!    Written by Eirik F. Kjønstad, Dec 2017
!!
      implicit none 
!
      class(sccsd) :: wf 
!
      integer(i15), intent(in) :: iteration 
!
      real(dp), dimension(:,:), allocatable :: rA 
      real(dp), dimension(:,:), allocatable :: rB
!
      real(dp), dimension(:,:), allocatable :: prev_rA 
      real(dp), dimension(:,:), allocatable :: prev_rB
!
      real(dp) :: ddot
!
      integer(i15) :: ioerror = 0
!
      real(dp) :: dot_product
!
      integer(i15) :: unit_prev_solution = -1
      integer(i15) :: unit_solution = -1
!
!     Open the solutions file 
!
      call generate_unit_identifier(unit_solution)
      open(unit=unit_solution, file='right_valence', action='readwrite', status='unknown', &
            access='direct', form='unformatted', recl=dp*(wf%n_parameters), iostat=ioerror) 
!
!     Read state A from disk 
!
      call wf%mem%alloc(rA, wf%n_parameters, 1)
      read(unit_solution, rec=wf%state_A, iostat=ioerror) rA 
!
!     Open previous solution file 
!
      call generate_unit_identifier(unit_prev_solution)
      open(unit=unit_prev_solution, file='prev_right_valence', action='readwrite', status='unknown', &
            access='direct', form='unformatted', recl=dp*(wf%n_parameters), iostat=ioerror) 
!
      if (iteration .eq. 1) then 
!
!        Simply write the current solutions to the previous solutions file 
!
         write(unit_prev_solution, rec=1, iostat=ioerror) rA
!
      else 
!
!        Read in previous solution 
!
         call wf%mem%alloc(prev_rA, wf%n_parameters, 1)
         read(unit_prev_solution, rec=1, iostat=ioerror) prev_rA 
!
!        Calculate dot product between current and previous solution 
!
         dot_product = zero
         dot_product = ddot(wf%n_parameters, prev_rA, 1, rA, 1)
!
!        If dot product is negative, then we are pretty sure the vector has flipped 
!
         if (dot_product .lt. zero) then 
!
            write(unit_output,'(/t6,a)')          'Eigencontroller: I think state A has flipped (r -> -r).'
            write(unit_output,'(t6,a42,f16.12)')  'Overlap with solution from last iteration:', dot_product
!
            write(unit_output,'(/t6,a)')          'Flipping back & saving consistent solution to file for the overlap calculation.'
!
            rA = -rA 
            write(unit_solution, rec=wf%state_A, iostat=ioerror) rA  
            write(unit_prev_solution, rec=1, iostat=ioerror) rA
!
         else
!
            write(unit_prev_solution, rec=1, iostat=ioerror) rA
!
         endif 
!
         call wf%mem%dealloc(prev_rA, wf%n_parameters, 1)
!
      endif 
!
      call wf%mem%dealloc(rA, wf%n_parameters, 1)
!
!     Read state B from disk 
!
      call wf%mem%alloc(rB, wf%n_parameters, 1)
      read(unit_solution, rec=wf%state_B, iostat=ioerror) rB 
!
!
      if (iteration .eq. 1) then 
!
!        Simply write the current solutions to the previous solutions file 
!
         write(unit_prev_solution, rec=2, iostat=ioerror) rB
!
      else 
!
!        Read in previous solution 
!
         call wf%mem%alloc(prev_rB, wf%n_parameters, 1)
         read(unit_prev_solution, rec=2, iostat=ioerror) prev_rB 
!
!        Calculate dot product between current and previous solution 
!
         dot_product = zero
         dot_product = ddot(wf%n_parameters, prev_rB, 1, rB, 1)
!
!        If dot product is negative, then we are pretty sure the vector has flipped 
!
         if (dot_product .lt. zero) then 
!
            write(unit_output,'(/t6,a)')          'Eigencontroller: I think state B has flipped (r -> -r).'
            write(unit_output,'(t6,a42,f16.12)')  'Overlap with solution from last iteration:', dot_product
!
            write(unit_output,'(/t6,a)')          'Flipping back & saving consistent solution to file for the overlap calculation.'
!
            rB = -rB 
            write(unit_solution, rec=wf%state_B, iostat=ioerror) rB 
            write(unit_prev_solution, rec=2, iostat=ioerror) rB
!
         else
!
            write(unit_prev_solution, rec=2, iostat=ioerror) rB
!
         endif 
!
         call wf%mem%dealloc(prev_rB, wf%n_parameters, 1)
!
      endif 
!
      call wf%mem%dealloc(rB, wf%n_parameters, 1)
!
!     Close files 
!
      close(unit_prev_solution)
      close(unit_solution)
!
!     Deallocations 
!
      call wf%mem%dealloc(rA, wf%n_parameters, 1)
      call wf%mem%dealloc(rB, wf%n_parameters, 1)
!
   end subroutine eigenvector_controller_sccsd
!
!
   module subroutine ground_state_eigenvector_controller_sccsd(wf, iteration)
!!
!!    Ground state eigenvector controller (SCCSD)
!!    Written by Eirik F. Kjønstad, Feb 2018
!!
      implicit none 
!
      class(sccsd) :: wf 
!
      integer(i15), intent(in) :: iteration 
!
      real(dp), dimension(:,:), allocatable :: rB
!
      real(dp), dimension(:,:), allocatable :: prev_rB
!
      real(dp) :: ddot
!
      integer(i15) :: ioerror = 0
!
      real(dp) :: dot_product
!
      integer(i15) :: unit_prev_solution = -1
      integer(i15) :: unit_solution = -1
!
!     Open the solutions file 
!
      call generate_unit_identifier(unit_solution)
      open(unit=unit_solution, file='right_valence', action='readwrite', status='unknown', &
            access='direct', form='unformatted', recl=dp*(wf%n_parameters), iostat=ioerror) 
!
!     Open previous solution file 
!
      call generate_unit_identifier(unit_prev_solution)
      open(unit=unit_prev_solution, file='prev_right_valence', action='readwrite', status='unknown', &
            access='direct', form='unformatted', recl=dp*(wf%n_parameters), iostat=ioerror)
!
!     Read state B from disk 
!
      call wf%mem%alloc(rB, wf%n_parameters, 1)
      read(unit_solution, rec=wf%state_B, iostat=ioerror) rB 
!
!
      if (iteration .eq. 1) then 
!
!        Simply write the current solutions to the previous solutions file 
!
         write(unit_prev_solution, rec=1, iostat=ioerror) rB
!
      else 
!
!        Read in previous solution 
!
         call wf%mem%alloc(prev_rB, wf%n_parameters, 1)
         read(unit_prev_solution, rec=1, iostat=ioerror) prev_rB 
!
!        Calculate dot product between current and previous solution 
!
         dot_product = zero
         dot_product = ddot(wf%n_parameters, prev_rB, 1, rB, 1)
!
!        If dot product is negative, then we are pretty sure the vector has flipped 
!
         if (dot_product .lt. zero) then 
!
            write(unit_output,'(/t6,a)')          'Eigencontroller: I think the excited state has flipped (r -> -r).'
            write(unit_output,'(t6,a42,f16.12)')  'Overlap with solution from last iteration:', dot_product
!
            write(unit_output,'(/t6,a)')          'Flipping back & saving consistent solution to file for the overlap calculation.'
!
            rB = -rB 
            write(unit_solution, rec=wf%state_B, iostat=ioerror) rB 
            write(unit_prev_solution, rec=1, iostat=ioerror) rB
!
         else
!
            write(unit_prev_solution, rec=1, iostat=ioerror) rB
!
         endif 
!
         call wf%mem%dealloc(prev_rB, wf%n_parameters, 1)
!
      endif 
!
      call wf%mem%dealloc(rB, wf%n_parameters, 1)
!
!     Close files 
!
      close(unit_prev_solution)
      close(unit_solution)
!
!     Deallocations 
!
      call wf%mem%dealloc(rB, wf%n_parameters, 1)
!
   end subroutine ground_state_eigenvector_controller_sccsd
!
!
end submodule excited_state
