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
!
!
   module subroutine F_ccsd_a2_2_ccsd(wf, c_ai, rho_aibj, tbar_aibj)
!!
!!    F transformation a2,2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                      :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout)   :: rho_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)      :: tbar_aibj
!
   end subroutine F_ccsd_a2_2_ccsd
!
!
   module subroutine F_ccsd_b2_2_ccsd(wf, c_ai, rho_aibj, tbar_aibj)
!!
!!    F transformation b2,2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                      :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout)   :: rho_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)      :: tbar_aibj
!
   end subroutine F_ccsd_b2_2_ccsd
!
!
   module subroutine F_ccsd_c2_2_ccsd(wf, c_ai, rho_aibj, tbar_aibj)
!!
!!    F transformation c2,2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                      :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout)   :: rho_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)      :: tbar_aibj
!
   end subroutine F_ccsd_c2_2_ccsd
!
!
   module subroutine F_ccsd_d2_2_ccsd(wf, c_ai, rho_aibj, tbar_aibj)
!!
!!    F transformation d2,2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                      :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout)   :: rho_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)      :: tbar_aibj
!
   end subroutine F_ccsd_d2_2_ccsd
!
!
   module subroutine F_ccsd_e2_2_ccsd(wf, c_ai, rho_aibj, tbar_aibj)
!!
!!    F transformation e2,2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                      :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout)   :: rho_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)      :: tbar_aibj
!
   end subroutine F_ccsd_e2_2_ccsd
!
!
   module subroutine F_ccsd_f2_2_ccsd(wf, c_aibj, rho_aibj, tbar_aibj)
!!
!!    F transformation f2,2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)      :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout)   :: rho_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)      :: tbar_aibj
!
   end subroutine F_ccsd_f2_2_ccsd
!
!
   module subroutine F_ccsd_g2_2_ccsd(wf, c_aibj, rho_aibj, tbar_aibj)
!!
!!    F transformation g2,2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)      :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout)   :: rho_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)      :: tbar_aibj
!
   end subroutine F_ccsd_g2_2_ccsd
