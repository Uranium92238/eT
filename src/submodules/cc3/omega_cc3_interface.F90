!
   module subroutine construct_omega_cc3(wf, omega)
!!
!!   Construct omega (CC3)
!!   Alex C. Paul and Rolf H. Myhre 2018
!!
     implicit none
!
     class(cc3), intent(in) :: wf
!
     real(dp), dimension(wf%n_amplitudes, 1), intent(inout) :: omega
!
   end subroutine construct_omega_cc3
!
!
   module subroutine omega_cc3_a_cc3(wf, omega1, omega2)
!!
!!    CC3 Omega terms
!!    Alex C. Paul and Rolf H. Myhre 2018
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: omega1
      real(dp), dimension((wf%n_v)*(wf%n_o)*((wf%n_v)*(wf%n_o) +1)/2, 1), intent(inout) :: omega2
!
   end subroutine omega_cc3_a_cc3
!
!
   module subroutine omega_cc3_integrals_cc3(wf)
!!
!!    Construct integrals need in CC3 Omega and store on disk
!!    Rolf H. Myhre, January 2019
!!
      implicit none
!
      class(cc3) :: wf
!
   end subroutine omega_cc3_integrals_cc3
!
