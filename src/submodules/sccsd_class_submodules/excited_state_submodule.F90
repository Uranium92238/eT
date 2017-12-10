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
      logical :: differentiate_overlap = .true.
!
      real(dp) :: previous_overlap
      real(dp) :: derivative_overlap
!
      real(dp) :: displacement
!
      integer(i15) :: iteration = 1
      integer(i15) :: max_iterations = 50
!
      integer(i15) :: unit_triples
!
      integer(i15) :: i = 0
!
      real(dp) :: begin_timer, end_timer
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
      write(unit_output,'(t6,a)')  '  I. Solve ground state equations'
      write(unit_output,'(t6,a)')  ' II. Solve multiplier equation'
      write(unit_output,'(t6,a)')  'III. Solve excited state eigenvalue problem'
      write(unit_output,'(t6,a/)') ' IV. Calculate overlap'
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
!     Set displacement for numerical differentiation 
!
      displacement = wf%scc_settings%overlap_threshold/10D0
!
!     Enter self-consistent similarity constrained CC loop 
!
      do while (.not. converged .and. iteration .le. max_iterations)
!
         write(unit_output,'(/t3,a29,i3)')     'Similarity constrained cycle:', iteration 
         write(unit_output,'(t3,a53,f16.12)')  'In this cycle, we use the following triple amplitude:', wf%triples
!
!        Solve for the ground state 
!
         write(unit_output,'(/t3,a,i2)') 'Part I. Solving the ground state equation'
         wf%tasks%current = 'ground_state'
         call wf%ground_state_driver
!
!        Solve for the multipliers
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
!        Solve for the right excited states 
!
         write(unit_output,'(t3,a)')  'Part III. Solving the eigenvalue equation for the right eigenvectors'
!
         write(unit_output,'(/t3,a)')   ':: Excited state solver (Davidson)'
         write(unit_output,'(t3,a)')    ':: E. F. Kjønstad, S. D. Folkestad, May 2017'
!
         wf%excited_state_specifications%right = .true.  ! Use A transformation
         wf%excited_state_specifications%left  = .false. ! Not A^T transformation
         wf%tasks%current = 'excited_state'
!
         call wf%excited_state_solver
         call wf%excited_state_cleanup
!
!        Calculate the overlap between the similarity constrained states 
!
         write(unit_output,'(/t3,a)') 'Part IV. Computing the generalized overlap between the states...'
         call wf%calc_overlap
!
         if (iteration .ne. 1) then 
!
!           Control whether one of the vectors have flipped (r -> -r), 
!           causing a change in the sign of the overlap 
!
            if (abs(previous_overlap + wf%overlap) .lt. wf%scc_settings%overlap_threshold) then 
!
               write(unit_output,'(/t3,a)') 'It seems that the overlap has flipped sign. Reversing.'
               wf%overlap = - wf%overlap 
!
            endif
!
         endif
!
         write(unit_output,'(/t3,a)') 'Summary of similarity constrained cycle:' 
!
         write(unit_output,'(/t6,a20,f16.12)') 'Triple amplitude:   ', wf%triples
         write(unit_output,'(t6,a20,f16.12)')  'Generalized overlap:', wf%overlap 
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
            if (differentiate_overlap) then ! Numerically differentiate overlap for Newton-Raphson step 
!
               write(unit_output,'(/t3,a,a)') 'Doing a small displacement to ', & 
                                              'get the derivative of the overlap.'
!
               wf%triples = wf%triples + displacement
!
               differentiate_overlap = .false. ! Next time, do Newton-Raphson
!
            else ! Perform Newton-Raphson step
!
               derivative_overlap = (wf%overlap - previous_overlap)/displacement 
!
               wf%triples = wf%triples - previous_overlap/derivative_overlap
!
               differentiate_overlap = .true. ! Next time, differentiate overlap 
!
            endif
!
            iteration = iteration + 1
!
         endif
!
         previous_overlap = wf%overlap 
!
      enddo
!
      if (.not. converged) then ! Let the user know before quitting
!
         write(unit_output,'(/t3,a)') 'Reached maximum number of cycles without convergence.'
!
      else ! Display the final energies, etc. 
!
         write(unit_output,'(t3,a)') 'Summary of SCCSD ground and excited state calculation:'
!
         write(unit_output,'(/t6,a27,f19.12)') 'Triple amplitude:          ', wf%triples
         write(unit_output,'(t6,a27,f19.12/)') 'Ground state energy [a.u.]:', wf%energy
!
         write(unit_output,'(t6,a10,4x,a13,11x,a11,9x,a14)')'Excitation', 'energy [a.u.]', 'energy [eV]', 'energy [cm^-1]'
         write(unit_output,'(t6,a)')'---------------------------------------------------------------------------------'
!
         do i = 1, wf%excited_state_specifications%n_singlet_states
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
   end subroutine excited_state_driver_sccsd
!
!
end submodule excited_state