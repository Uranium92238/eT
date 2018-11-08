submodule (ccsd_class) omega_ccsd
!
!!
!!    Omega submodule (ccsd)
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 
!!    and Andreas Skeidsvoll, 2018
!!
!!    Routines to construct 
!!
!!    Ω =  < mu | exp(-T) H exp(T) | R >
!!
!
   implicit none
!
!
contains
!
!
   module subroutine construct_omega_cc3(wf, omega)
!!
!!    Construct omega (CC3)
!!    Written by Alex C. Paul and Rolf H. Myhre 2018
!!
!!    Directs the construction of the projection vector < mu | exp(-T) H exp(T) | R >
!!    for the current amplitudes of the object wfn
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_amplitudes, 1), intent(inout) :: omega
!
      real(dp), dimension(:,:), allocatable :: omega1, omega2
!
      call mem%alloc(omega1, wf%n_v, wf%n_o)
      call mem%alloc(omega2, (wf%n_v)*(wf%n_o)*((wf%n_v)*(wf%n_o) +1)/2, 1)
!
!     Set the omega vector to zero
!
      omega1 = zero
      omega2 = zero
!
!     Construct singles contributions
!
      call wf%omega_ccsd_a1(omega1)
      call wf%omega_ccsd_b1(omega1)
      call wf%omega_ccsd_c1(omega1)
!
      call wf%omega_ccs_a1(omega1)
!
!     Construct doubles contributions
!
      call wf%omega_ccsd_a2(omega2)
      call wf%omega_ccsd_b2(omega2)
      call wf%omega_ccsd_c2(omega2)
      call wf%omega_ccsd_d2(omega2)
      call wf%omega_ccsd_e2(omega2)
!
      call dcopy((wf%n_o)*(wf%n_v), omega1, 1, omega, 1)
      call dcopy((wf%n_o)*(wf%n_v)*((wf%n_o)*(wf%n_v)+1)/2, omega2, 1, omega((wf%n_o)*(wf%n_v)+1, 1), 1)
!
      call mem%dealloc(omega1, wf%n_v, wf%n_o)
      call mem%dealloc(omega2, (wf%n_v)*(wf%n_o)*((wf%n_v)*(wf%n_o) +1)/2, 1)
!
   end subroutine construct_omega_cc3
!
!
end submodule
