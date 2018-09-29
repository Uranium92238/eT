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