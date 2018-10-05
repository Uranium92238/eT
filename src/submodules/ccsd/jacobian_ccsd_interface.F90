
   module subroutine jacobian_transform_trial_vector_ccsd(wf, c_i)
!!
!!    Jacobian transform trial vector 
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!    and Andreas Skeidsvoll, 2018
!!
      class(ccsd), intent(in) :: wf 
!
      real(dp), dimension(wf%n_amplitudes, 1) :: c_i
!
   end subroutine jacobian_transform_trial_vector_ccsd
!
!
   module subroutine jacobian_ccsd_transformation_ccsd(wf, c)
!!
!!    Jacobian transformation (CCSD)
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!    and Andreas Skeidsvoll, 2018
!!
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_amplitudes, 1) :: c
!
   end subroutine jacobian_ccsd_transformation_ccsd
!
!
  module subroutine jacobian_ccsd_a1_ccsd(wf, rho_a_i, c_a_i)
!!
!!    Jacobian CCSD A1
!!    Written by Sarai D. Folkestad Eirik F. Kjønstad
!!    and Andreas Skeidsvoll, 2018
!!
      implicit none
!
      class(ccsd) :: wf
!
!     Vectors sent to the routine
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_a_i
      real(dp), dimension(wf%n_v, wf%n_o) :: rho_a_i
!
   end subroutine jacobian_ccsd_a1_ccsd
!
!
   module subroutine jacobian_ccsd_b1_ccsd(wf, rho_a_i, c_ai_bj)
!!
!!    Jacobian CCSD B1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_o)*(wf%n_v),(wf%n_o)*(wf%n_v)), intent(in) :: c_ai_bj
      real(dp), dimension(wf%n_v, wf%n_o) :: rho_a_i
!
   end subroutine jacobian_ccsd_b1_ccsd
!
!
   module subroutine jacobian_ccsd_c1_ccsd(wf, rho_a_i, c_ai_bj)
!!
!!    Jacobian CCSD C1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!    and Andreas Skeidsvoll, 2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_o)*(wf%n_v),(wf%n_o)*(wf%n_v)), intent(in) :: c_ai_bj
      real(dp), dimension(wf%n_v, wf%n_o) :: rho_a_i
!
   end subroutine jacobian_ccsd_c1_ccsd
!
!
    module subroutine jacobian_ccsd_d1_ccsd(wf, rho_a_i, c_bi_cj)
!!
!!    Jacobian CCSD D1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!    and Andreas Skeidsvoll, 2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_o)*(wf%n_v),(wf%n_o)*(wf%n_v)), intent(in) :: c_bi_cj
      real(dp), dimension(wf%n_v, wf%n_o) :: rho_a_i
!
   end subroutine jacobian_ccsd_d1_ccsd
!
!
   module subroutine jacobian_ccsd_a2_ccsd(wf, rho_ai_bj, c_a_i)
!!
!!    Jacobian CCSD A2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)           :: c_a_i
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
!
   end subroutine jacobian_ccsd_a2_ccsd
!
!
   module subroutine jacobian_ccsd_b2_ccsd(wf, rho_ai_bj, c_a_i)
!!
!!    Jacobian CCSD B2
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad
!!    and Andreas Skeidsvoll, 2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)           :: c_a_i
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
!
   end subroutine jacobian_ccsd_b2_ccsd
!
!
   module subroutine jacobian_ccsd_c2_ccsd(wf, rho_ai_bj, c_a_i)
!!
!!    Jacobian CCSD C2
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad
!!    and Andreas Skeidsvoll, 2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)           :: c_a_i
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
!
   end subroutine jacobian_ccsd_c2_ccsd
!
!
  module subroutine jacobian_ccsd_d2_ccsd(wf, rho_ai_bj, c_a_i)
!!
!!    Jacobian CCSD D2
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad
!!    and Andreas Skeidsvoll, 2018
!!
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)           :: c_a_i
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
!
   end subroutine jacobian_ccsd_d2_ccsd
!
!
    module subroutine jacobian_ccsd_e2_ccsd(wf, rho_ai_bj, c_ai_ck)
!!
!!    Jacobian CCSD E2
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!    and Andreas Skeidsvoll, 2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)), intent(in) :: c_ai_ck
!
   end subroutine jacobian_ccsd_e2_ccsd
!
!
   module subroutine jacobian_ccsd_f2_ccsd(wf, rho_ai_bj, c_ai_bj)
!!
!!       Jacobian CCSD F2
!!       Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!       and Andreas Skeidsvoll, 2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)), intent(in) :: c_ai_bj
!
   end subroutine jacobian_ccsd_f2_ccsd
!
!
   module subroutine jacobian_ccsd_g2_ccsd(wf, rho_ai_bj, c_ai_bj)
!!
!!    Jacobian CCSD G2
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!    and Andreas Skeidsvoll, 2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)), intent(in) :: c_ai_bj
!
   end subroutine jacobian_ccsd_g2_ccsd
!
!
   module subroutine jacobian_ccsd_h2_ccsd(wf, rho_ai_bj, c_ai_bj)
!!
!!       Jacobian CCSD H2
!!       Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!       and Andreas Skeidsvoll
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
         real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)), intent(in) :: c_ai_bj
!
   end subroutine jacobian_ccsd_h2_ccsd
!
!
   module subroutine jacobian_ccsd_i2_ccsd(wf, rho_ai_bj, c_ai_bj)
!!
!!    Jacobian CCSD I2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad
!!    and Andreas Skeidsvoll, 2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)) :: rho_ai_bj
      real(dp), dimension((wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v)), intent(in) :: c_ai_bj
!
   end subroutine jacobian_ccsd_i2_ccsd
!
!
   module subroutine jacobian_ccsd_j2_ccsd(wf, rho_ab_ij, c_ab_ij)
!!
!!     Jacobian CCSD J2
!!     Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!     and Andreas Skeidsvoll, 2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_v)**2, (wf%n_o)**2) :: rho_ab_ij
      real(dp), dimension((wf%n_v)**2, (wf%n_o)**2), intent(in) :: c_ab_ij
!
   end subroutine jacobian_ccsd_j2_ccsd
!
!
   module subroutine jacobian_ccsd_k2_ccsd(wf, rho_ab_ij, c_ab_ij)
!!
!!    Jacobian CCSD K2
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!    and Andreas Skeidsvoll, 2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_v)**2, (wf%n_o)**2) :: rho_ab_ij
      real(dp), dimension((wf%n_v)**2, (wf%n_o)**2), intent(in) :: c_ab_ij
!
   end subroutine jacobian_ccsd_k2_ccsd
!