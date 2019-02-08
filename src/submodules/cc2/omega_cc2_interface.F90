
   module subroutine construct_omega_cc2(wf, omega)
!!
!!    Construct omega 
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 
!!    and Linda Goletto, 2018
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_amplitudes, 1), intent(inout) :: omega
!
   end subroutine construct_omega_cc2
!
!
   module subroutine omega_cc2_a1_cc2(wf, omega, eps_o, eps_v)
!!
!!    Omega CC2 A1 term
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 
!!    and Linda Goletto, 2018
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_amplitudes, 1), intent(inout) :: omega
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
   end subroutine omega_cc2_a1_cc2
!
!
   module subroutine omega_cc2_b1_cc2(wf, omega, eps_o, eps_v)
!!
!!    Omega CC2 B1 term
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 
!!    and Linda Goletto, 2018
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_amplitudes, 1), intent(inout) :: omega
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
   end subroutine omega_cc2_b1_cc2
!
!
   module subroutine omega_cc2_c1_cc2(wf, omega, eps_o, eps_v)
!!
!!    Omega CC2 C1 term
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 
!!    and Linda Goletto, 2018
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_amplitudes, 1), intent(inout) :: omega
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
   end subroutine omega_cc2_c1_cc2