
   module subroutine omega_cc3_a_cc3(wf, omega)
!!
!!    CC3 Omega terms
!!    Alex C. Paul and Rolf H. Myhre 2018
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_amplitudes, 1), intent(inout) :: omega
!
   end subroutine omega_cc3_a_cc3
!
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

