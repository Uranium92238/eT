
   module subroutine omega_ccsd_a1_ccsd(wf, omega1)
!!
!!    Omega A1 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout):: omega1
!
   end subroutine omega_ccsd_a1_ccsd
!
!
   module subroutine omega_ccsd_b1_ccsd(wf, omega1)
!!
!!    Omega B1 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout):: omega1
!
   end subroutine omega_ccsd_b1_ccsd
!
!
   module subroutine omega_ccsd_c1_ccsd(wf, omega1)
!!
!!    Omega C1 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout):: omega1
!
   end subroutine omega_ccsd_c1_ccsd
!
!
   module subroutine omega_ccsd_a2_ccsd(wf, omega2)
!!
!!    Omega A2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_v)*(wf%n_o)*((wf%n_v)*(wf%n_o)+1)/2, 1), intent(inout):: omega2
!
   end subroutine omega_ccsd_a2_ccsd
!
!
   module subroutine omega_ccsd_b2_ccsd(wf, omega2)
!!
!!    Omega B2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_v)*(wf%n_o)*((wf%n_v)*(wf%n_o)+1)/2, 1), intent(inout):: omega2
!
   end subroutine omega_ccsd_b2_ccsd
!
!
   module subroutine omega_ccsd_c2_ccsd(wf, omega2)
!!
!!    Omega C2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_v)*(wf%n_o)*((wf%n_v)*(wf%n_o)+1)/2, 1), intent(inout):: omega2
!
   end subroutine omega_ccsd_c2_ccsd
!
!
   module subroutine omega_ccsd_d2_ccsd(wf, omega2)
!!
!!    Omega D2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_v)*(wf%n_o)*((wf%n_v)*(wf%n_o)+1)/2, 1), intent(inout):: omega2
!
   end subroutine omega_ccsd_d2_ccsd
!
!
   module subroutine omega_ccsd_e2_ccsd(wf, omega2)
!!
!!    Omega E2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_v)*(wf%n_o)*((wf%n_v)*(wf%n_o)+1)/2, 1), intent(inout):: omega2
!
   end subroutine omega_ccsd_e2_ccsd
!
!
   module subroutine construct_omega_ccsd(wf, omega)
!!
!!    Construct omega (CCSD)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
         class(ccsd), intent(inout) :: wf
!
         real(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: omega
!
   end subroutine construct_omega_ccsd
!
!
