submodule (ccsdpt_class) ground_state
!
!!    
!!     Ground state submodule (CCSD(T))
!!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!    
!!    
!!     Consists of the following subroutines of the CCSD(T) module:
!!     
!!     destruct_ground_state: before destroying the converged CCSD amplitudes,
!!                            the non-iterative CCSD(T) correction is to be added
!!                            to the energy (call wf%calc_energy_correction).
!!    
!
   implicit none 
!
!
contains
!
!
   module subroutine destruct_ground_state_ccsdpt(wf)
!!
!!    Destruct Ground State (CCSD(T))
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
!!    Calculates the CCSD(T) energy correction to the CCSD
!!    energy, thereafter deallocating the amplitudes and the 
!!    projection vector.
!!
      implicit none
!
      class(ccsdpt) :: wf
!  
      call wf%calc_energy_correction ! Add CCSD(T) energy correction 
      call wf%destruct_amplitudes
      call wf%destruct_omega
!
   end subroutine destruct_ground_state_ccsdpt
!
!
end submodule ground_state 
