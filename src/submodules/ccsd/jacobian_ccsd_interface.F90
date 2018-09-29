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