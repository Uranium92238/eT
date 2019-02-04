!
!
   module subroutine F_transform_vector_ccsd(wf, c)
!!
!!    F transformation
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018     
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_amplitudes, 1), intent(inout) :: c
!
   end subroutine F_transform_vector_ccsd
!
!
   module subroutine F_ccsd_a1_1_ccsd(wf, c_aibj, rho_ai)
!!
!!    F transformation A1,1 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)      :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)                   :: rho_ai
!
   end subroutine F_ccsd_a1_1_ccsd
!
!
   module subroutine F_ccsd_a2_1_ccsd(wf, c_ai, rho_aibj)
!!
!!    F transformation A2,1 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)      :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout)   :: rho_aibj
!
   end subroutine F_ccsd_a2_1_ccsd
!
!
!
!
   module subroutine F_ccsd_a1_2_ccsd(wf, c_ai, rho_ai, tbar_aibj)
!!
!!    F transformation A1,2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                   :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)                :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)   :: tbar_aibj
!
   end subroutine F_ccsd_a1_2_ccsd
!
!
   module subroutine F_ccsd_b1_2_ccsd(wf, c_ai, rho_ai, tbar_aibj)
!!
!!    F transformation B1,2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                   :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)                :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)   :: tbar_aibj
!
   end subroutine F_ccsd_b1_2_ccsd
!
!
   module subroutine F_ccsd_c1_2_ccsd(wf, c_ai, rho_ai, tbar_aibj)
!!
!!    F transformation C1,2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                   :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)                :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)   :: tbar_aibj
!
   end subroutine F_ccsd_c1_2_ccsd
!
!
   module subroutine F_ccsd_d1_2_ccsd(wf, c_aibj, rho_ai, tbar_aibj)
!!
!!    F transformation D1,2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)   :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)                :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)   :: tbar_aibj
!
   end subroutine F_ccsd_d1_2_ccsd
!
!
   module subroutine F_ccsd_e1_2_ccsd(wf, c_aibj, rho_ai, tbar_aibj)
!!
!!    F transformation E1,2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)   :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)                :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)   :: tbar_aibj
!
   end subroutine F_ccsd_e1_2_ccsd
!
!
   module subroutine F_ccsd_f1_2_ccsd(wf, c_aibj, rho_ai, tbar_aibj)
!!
!!    F transformation F1,2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_F1,2 = -(F_ib tbar_ckaj + F_ja tbar_ckbi) c_bjck 
!!
!!    The first two terms of (67)
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)              :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: tbar_aibj
!
   end subroutine F_ccsd_f1_2_ccsd
!
!
   module subroutine F_ccsd_g1_2_ccsd(wf, c_aibj, rho_ai, tbar_aibj)
!!
!!    F transformation G1,2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_G1,2 = -(L_iljb tbar_ckal + L_jlia tbar_ckbl) c_bjck 
!!
!!    The third and fourth terms of (67)
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)              :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: tbar_aibj
!
   end subroutine F_ccsd_g1_2_ccsd
!
!
   module subroutine F_ccsd_h1_2_ccsd(wf, c_aibj, rho_ai, tbar_aibj)
!!
!!    F transformation H1,2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_H1,2 = (L_dajb tbar_ckdi + L_dbia tbar_ckdj) c_bjck 
!!
!!    The last two terms of (67)
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)              :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: tbar_aibj
!
   end subroutine F_ccsd_h1_2_ccsd
!
!
   module subroutine F_ccsd_i1_2_ccsd(wf, c_ai, rho_ai, tbar_aibj, t_aibj)
!!
!!    F transformation I1,2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_I1,2 = - L_ldic tbar_bjak (t_bjck c_dl + t_bjdl c_ck) 
!!               - L_ialc tbar_bjdk (t_bjck c_dl + t_bjdl c_ck)
!!               - L_ldka tbar_bjci (t_bjck c_dl + t_bjdl c_ck)
!!
!!    Equation (23) in Sarai's document
!!    Equation (69) in Eirik's document
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                 :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)              :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: tbar_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: t_aibj
!
   end subroutine F_ccsd_i1_2_ccsd
!
!
   module subroutine F_ccsd_j1_2_ccsd(wf, c_ai, rho_ai, tbar_aibj, t_aibj)
!!
!!    F transformation J1,2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_J1,2 = (g_kbid tbar_cjal + g_jcid tbar_bkal + g_kajd tbar_cibl 
!!                g_jakd tbar_bicl + g_ibkd tbar_ajcl + g_lakb tbar_dicj) t_ckdl c_bj
!!
!!    Equation (68)
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                 :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)              :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: tbar_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: t_aibj
!
   end subroutine F_ccsd_j1_2_ccsd
!
!
