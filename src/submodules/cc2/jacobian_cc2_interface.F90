
   module subroutine effective_jacobian_transformation_cc2(wf, c)
!!
!!    Effective jacobian transformation
!!    Written by Eirik F. Kjønstad and Sarai Dery Folkestad
!!    Linda Goletto, and Alexander Paul, Dec 2018
!!
      implicit none
!
      class(cc2) :: wf
!
      real(dp), dimension(wf%n_amplitudes, 1) :: c
!
   end subroutine effective_jacobian_transformation_cc2
!
!
   module subroutine jacobian_cc2_A1_cc2(wf, rho_ai, c_bj)
!!
!!    Jacobian CC2 A1
!!    Written by Eirik F. Kjønstad and Sarai Dery Folkestad
!!    Linda Goletto, and Alexander Paul, Dec 2018
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
!     Vectors sent to the routine
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)     :: c_bj
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)  :: rho_ai
!
    end subroutine jacobian_cc2_A1_cc2
!
!
    module subroutine jacobian_cc2_B1_cc2(wf, rho_ai, c_bj, eps_o, eps_v)
!!
!!    Jacobian CC2 B1
!!    Written by Eirik F. Kjønstad and Sarai Dery Folkestad
!!    Linda Goletto, and Alexander Paul, Dec 2018
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
!     Vectors sent to the routine
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_bj
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
   end subroutine jacobian_cc2_B1_cc2
!
!
   module subroutine effective_jacobian_cc2_a1_cc2(wf, omega, rho_ai, c_ai, eps_o, eps_v)
!!
!!    Effective Jacobian CC2 A1
!!    Written by Eirik F. Kjønstad and Sarai Dery Folkestad
!!    Linda Goletto, and Alexander Paul, Jan 2019
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
   end subroutine effective_jacobian_cc2_a1_cc2
!
!
   module subroutine effective_jacobian_cc2_b1_cc2(wf, omega, rho_ai, c_ai, eps_o, eps_v)
!!
!!    Effective jacobian B1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!    Linda Goletto, and Alexander Paul, Jan 2019
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
!
      real(dp), dimension(wf%n_o), intent(in)  :: eps_o
      real(dp), dimension(wf%n_v), intent(in)  :: eps_v
!
   end subroutine effective_jacobian_cc2_b1_cc2
!
!
   module subroutine effective_jacobian_cc2_c1_cc2(wf, omega, rho_ai, c_cj, eps_o, eps_v)
!!
!!    Jacobian CC2 E1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,
!!    Linda Goletto, and Alexander Paul, Jan 2019
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_cj
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
end subroutine effective_jacobian_cc2_c1_cc2
!
!
   module subroutine effective_jacobian_cc2_d1_cc2(wf, omega, rho_ai, c_bl, eps_o, eps_v)
!!
!!    Jacobian CC2 E1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,
!!    Linda Goletto, and Alexander Paul, Jan 2019
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_bl
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
end subroutine effective_jacobian_cc2_d1_cc2
!
!
   module subroutine effective_jacobian_cc2_e1_cc2(wf, omega, rho_ai, c_bl, eps_o, eps_v)
!!
!!    Jacobian CC2 F1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,
!!    Linda Goletto, and Alexander Paul, Dec 2018
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_bl
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
   end subroutine effective_jacobian_cc2_e1_cc2
!
!
   module subroutine effective_jacobian_cc2_f1_cc2(wf, omega, rho_ai, c_ai, eps_o, eps_v)
!!
!!    Jacobian CC2 f1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,
!!    Linda Goletto, and Alexander Paul, Dec 2018
!!
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
!     Sent to the routine
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
   end subroutine effective_jacobian_cc2_f1_cc2
!
!
