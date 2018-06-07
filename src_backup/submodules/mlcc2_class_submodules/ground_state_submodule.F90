submodule (mlcc2_class) ground_state 
!
!!
!!    Ground state submodule (mlcc2)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Jun 2017
!!
!!
!
   implicit none
!
contains
!
   module subroutine ground_state_preparations_mlcc2(wf)
!!
!!    Ground State Preparations (mlcc2)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
!!    A routine for preparation tasks (if any). Can be overwritten
!!    in descendants if other preparations prove necessary.    
!!
      implicit none
!
      class(mlcc2) :: wf 
!
      wf%tasks%current = 'ground_state'
!
!     Allocate amplitudes (if not allocated) and calculate number of amplitudes 
!
      call wf%initialize_single_amplitudes 
!
!     Allocate projection vector 
!
      call wf%initialize_omega   
!
   end subroutine ground_state_preparations_mlcc2
!
!
   module subroutine ground_state_cleanup_mlcc2(wf)
!!
!!    Ground State Cleanup (mlcc2)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
!!    A routine for cleanup tasks (if any). Can be overwritten
!!    in descendants if other cleanups prove necessary.    
!!
      implicit none
!
      class(mlcc2) :: wf 
!
      call wf%save_amplitudes
!
      call wf%destruct_single_amplitudes
!
      call wf%destruct_omega
!
   end subroutine ground_state_cleanup_mlcc2
!
!
end submodule ground_state