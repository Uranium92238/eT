submodule (cc2_class) omega_cc2
!
!!
!!    Omega submodule (cc2)
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 
!!    and Linda Goletto, 2018
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
   module subroutine construct_omega_cc2(wf, omega)
!!
!!    Construct omega 
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 
!!    and Linda Goletto, 2018
!!
!!    Directs the construction of the projection vector < mu | exp(-T) H exp(T) | R >
!!    for the current amplitudes of the object wfn
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_amplitudes, 1), intent(inout) :: omega
!
      omega = zero
!
      call wf%omega_ccs_a1(omega)
!
      call wf%omega_cc2_a1(omega)
      call wf%omega_cc2_b1(omega)
      call wf%omega_cc2_c1(omega)
!
   end subroutine construct_omega_cc2
!
!
   module subroutine omega_cc2_a1_cc2(wf, omega)
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
!
   end subroutine omega_cc2_a1_cc2
!
!
   module subroutine omega_cc2_b1_cc2(wf, omega)
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
!
   end subroutine omega_cc2_b1_cc2
!
!
   module subroutine omega_cc2_c1_cc2(wf, omega)
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
!
   end subroutine omega_cc2_c1_cc2
!
end submodule
