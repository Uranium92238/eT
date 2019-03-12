!
!
   module subroutine F_transform_vector_ccs(wf, c)
!!
!!    F transformation
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018     
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_amplitudes), intent(inout) :: c
!
   end subroutine F_transform_vector_ccs
!
!
   module subroutine F_ccs_a1_0_ccs(wf, c_ai, rho_ai)
!!
!!    F transformation A1,0 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_A1,0_ai = 2 L_iajb c_bj
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
!
   end subroutine F_ccs_a1_0_ccs
!
!
   module subroutine F_ccs_a1_1_ccs(wf, c_ai, rho_ai)
!!
!!    F transformation A1,1 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
!
   end subroutine F_ccs_a1_1_ccs
!
!
   module subroutine F_ccs_b1_1_ccs(wf, c_ai, rho_ai)
!!
!!    F transformation B1,1 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
!
   end subroutine F_ccs_b1_1_ccs
!
!
   module subroutine F_ccs_c1_1_ccs(wf, c_ai, rho_ai)
!!
!!    F transformation C1,1 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
!
   end subroutine F_ccs_c1_1_ccs
!
