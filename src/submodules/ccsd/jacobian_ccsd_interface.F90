   module subroutine jacobian_transform_trial_vector_ccsd(wf, c_i)
!!
!!    Jacobian transform trial vector
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes) :: c_i
!
   end subroutine jacobian_transform_trial_vector_ccsd
!
!
   module subroutine jacobian_ccsd_transformation_ccsd(wf, c)
!!
!!    Jacobian transformation (CCSD)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes) :: c
!
   end subroutine jacobian_ccsd_transformation_ccsd
!
!
  module subroutine jacobian_ccsd_a1_ccsd(wf, rho_ai, c_ai)
!!
!!    Jacobian CCSD A1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o) :: rho_ai
!
   end subroutine jacobian_ccsd_a1_ccsd
!
!
   module subroutine jacobian_ccsd_b1_ccsd(wf, rho_ai, c_aibj)
!!
!!    Jacobian CCSD B1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o) :: rho_ai
!
   end subroutine jacobian_ccsd_b1_ccsd
!
!
   module subroutine jacobian_ccsd_c1_ccsd(wf, rho_ai, c_aibj)
!!
!!    Jacobian CCSD C1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o) :: rho_ai
!
   end subroutine jacobian_ccsd_c1_ccsd
!
!
    module subroutine jacobian_ccsd_d1_ccsd(wf, rho_ai, c_bicj)
!!
!!    Jacobian CCSD D1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_bicj
      real(dp), dimension(wf%n_v, wf%n_o) :: rho_ai
!
   end subroutine jacobian_ccsd_d1_ccsd
!
!
   module subroutine jacobian_ccsd_a2_ccsd(wf, rho_aibj, c_ai)
!!
!!    Jacobian CCSD A2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)           :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: rho_aibj
!
   end subroutine jacobian_ccsd_a2_ccsd
!
!
   module subroutine jacobian_ccsd_b2_ccsd(wf, rho_aibj, c_ai)
!!
!!    Jacobian CCSD B2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)           :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: rho_aibj
!
   end subroutine jacobian_ccsd_b2_ccsd
!
!
   module subroutine jacobian_ccsd_c2_ccsd(wf, rho_aibj, c_ai)
!!
!!    Jacobian CCSD C2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)           :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: rho_aibj
!
   end subroutine jacobian_ccsd_c2_ccsd
!
!
  module subroutine jacobian_ccsd_d2_ccsd(wf, rho_aibj, c_ai)
!!
!!    Jacobian CCSD D2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)           :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: rho_aibj
!
   end subroutine jacobian_ccsd_d2_ccsd
!
!
    module subroutine jacobian_ccsd_e2_ccsd(wf, rho_aibj, c_aick)
!!
!!    Jacobian CCSD E2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: rho_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_aick
!
   end subroutine jacobian_ccsd_e2_ccsd
!
!
   module subroutine jacobian_ccsd_f2_ccsd(wf, rho_aibj, c_aibj)
!!
!!    Jacobian CCSD F2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: rho_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_aibj
!
   end subroutine jacobian_ccsd_f2_ccsd
!
!
   module subroutine jacobian_ccsd_g2_ccsd(wf, rho_aibj, c_aibj)
!!
!!    Jacobian CCSD G2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: rho_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_aibj
!
   end subroutine jacobian_ccsd_g2_ccsd
!
!
   module subroutine jacobian_ccsd_h2_ccsd(wf, rho_aibj, c_aibj)
!!
!!    Jacobian CCSD H2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: rho_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_aibj
!
   end subroutine jacobian_ccsd_h2_ccsd
!
!
   module subroutine jacobian_ccsd_i2_ccsd(wf, rho_aibj, c_aibj)
!!
!!    Jacobian CCSD I2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o) :: rho_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_aibj
!
   end subroutine jacobian_ccsd_i2_ccsd
!
!
   module subroutine jacobian_ccsd_j2_ccsd(wf, rho_abij, c_abij)
!!
!!    Jacobian CCSD J2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o) :: rho_abij
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: c_abij
!
   end subroutine jacobian_ccsd_j2_ccsd
!
!
   module subroutine jacobian_ccsd_k2_ccsd(wf, rho_abij, c_abij)
!!
!!    Jacobian CCSD K2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o) :: rho_abij
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: c_abij
!
   end subroutine jacobian_ccsd_k2_ccsd
!
