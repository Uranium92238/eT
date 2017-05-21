submodule (ccsdpt_class) omega
!
!!
!!    Omega submodule (CCSD(T))
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!
!!    Contains the following procedures of the CCSD(T) class:
!!
!!    construct_omega:  constructs the projection vector (omega1, omega2).
!!                      Calls the CCSD construction of omega (not CC3, its 
!!                      immediate parent).
!!
!
   implicit none 
!
!
contains
!
!
   subroutine construct_omega_ccsdpt(wf) 
!
!     Construct Omega (CCSD(T))
!     Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!
!     Directs the construction of the projection vector < mu | exp(-T) H exp(T) | R >
!     for the current amplitudes of the object wfn.
!
      class(ccsdpt) :: wf
!
      call construct_omega_ccsd(wf)
!
   end subroutine construct_omega_ccsdpt
!
! 
end submodule omega