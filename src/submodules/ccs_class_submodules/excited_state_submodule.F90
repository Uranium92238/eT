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
   integer(i15) :: iteration = -1 
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
!
      n_red        = wf%tasks%n_singlet_states
      n_new_trials = wf%tasks%n_singlet_states
!
      call wf%transform_trial_vecs(n_red - n_new_trials + 1, n_red)
!
   end subroutine excited_state_solver_ccs
!
!
end submodule