submodule (sccsd_class) excited_state
!
!!
!!    Excited state submodule (SCCSD)
!!    Written by Eirik F. Kjønstad, June 2017
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
!     Let the user know the excited state driver is running,
!     and print some basic information 
!
      write(unit_output,'(/t3,a)')    ':: Excited state solver (DIIS for the ground state, Davidson for the excited states)'
      write(unit_output,'(t3,a/)')    ':: E. F. Kjønstad, S. D. Folkestad, May-June 2017'
      write(unit_output,'(t3,a,i2,a,a,a)') &
               'Requested ',wf%excited_state_specifications%n_singlet_states,' ', trim(wf%name), ' singlet states.'
      write(unit_output,'(t3,a,i2,a,a,a)') &
               'Requested ',wf%excited_state_specifications%n_triplet_states,' ', trim(wf%name), ' triplet states.'     
!      
      write(unit_output, '(/t3,a,/t3,a,i1,a,i1,a)') &
                                     'Constraining the Jacobian matrix to be nondefective in', &
                                     'the subspace spanned by the states ', wf%state_A, ' and ', wf%state_B, '.'
!
      write(unit_output,'(/t3,a)') 'Triple excitation:'
!
      write(unit_output,'(/t6,a,3i2)') 'Occupied:', wf%I, wf%J, wf%K
      write(unit_output,'(t6,a,3i2/)') 'Virtual: ', wf%A, wf%B, wf%C
!
!     Set displacement for numerical differentiation 
!
      displacement = wf%scc_settings%overlap_threshold/10D0
!
!     Converge residuals and energies of the ground state to, at least, 
!     the accuracy of the overlap 
!
      if (wf%excited_state_specifications%energy_threshold .gt. wf%scc_settings%overlap_threshold) then 
! 
         wf%excited_state_specifications%energy_threshold = wf%scc_settings%overlap_threshold
!
      endif 
!
      if (wf%excited_state_specifications%residual_threshold .gt. wf%scc_settings%overlap_threshold) then 
! 
         wf%excited_state_specifications%residual_threshold = wf%scc_settings%overlap_threshold
!
      endif 
!
!     Converge residuals and energies of the excited states to, at least, 
!     the accuracy of the overlap 
!
      if (wf%excited_state_specifications%energy_threshold .gt. wf%scc_settings%overlap_threshold) then 
! 
         wf%excited_state_specifications%energy_threshold = wf%scc_settings%overlap_threshold
!
      endif 
!
      if (wf%excited_state_specifications%residual_threshold .gt. wf%scc_settings%overlap_threshold) then 
!
         wf%excited_state_specifications%residual_threshold = wf%scc_settings%overlap_threshold
!
      endif
!
!     Enter self-consistent similarity constrained CC loop 
!
      do while (.not. converged .and. iteration .le. max_iterations)
!
         write(unit_output,'(/t3,a,i2)') ':: Similarity constrained iteration ', iteration 
!
!        Solve for the ground state 
!
         write(unit_output,'(/t3,a,i2)') 'I. Solving the equations for the ground state:'
         wf%tasks%current = 'ground_state'
         call wf%ground_state_driver
!
!        Solve for the multipliers
!
         write(unit_output,'(/t3,a,i2)') 'II. Solving the equations for the multipliers:'
!
         call wf%excited_state_preparations
!
         wf%tasks%multipliers = .true.
         wf%tasks%current = 'multipliers'
         call wf%response_solver
         wf%tasks%multipliers = .false.
!
!        Solve for the right excited states 
!
         write(unit_output,'(/t3,a)') 'III. Solving the eigenvalue equation for the right eigenvectors:'
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
         write(unit_output,'(/t3,a)') 'IV. Computing the generalized overlap between the constrained states.'
         call wf%calc_overlap
!
         if (iteration .ne. 1) then 
!
!           Control whether one of the vectors have flipped (r -> -r), 
!           causing a change in the sign of the overlap 
!
            if (abs(previous_overlap + wf%overlap) .lt. wf%scc_settings%overlap_threshold) then 
!
               write(unit_output,'(/t3,a)') 'It seems that the overlap has switched sign. Reversing.'
               wf%overlap = - wf%overlap 
!
            endif
!
         endif
!
         write(unit_output,'(//t3,a,i2)') 'Similarity constrained CC summary, iteration', iteration 
!
         write(unit_output,'(/t3,a21,f16.12)') 'Triples amplitude:', wf%triples
         write(unit_output,'(t3,a21,f16.12)')  'Q*Q overlap:',       wf%overlap 
!
!        Check for convergence of the overlap 
!
         if (abs(wf%overlap) .le. wf%scc_settings%overlap_threshold) then 
!
            converged = .true.
            write(unit_output,'(/t3,a,i2,a/)')  'Converged in ', iteration, ' iterations!'
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
   end subroutine excited_state_driver_sccsd
!
!
end submodule excited_state