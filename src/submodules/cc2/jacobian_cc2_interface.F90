
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

   module subroutine construct_jacobian_cc2_A1_cc2(wf, rho_a_i, c_b_j)
 !!
 !!    Jacobian CC2 A1
 !!    Written by Eirik F. Kjønstad and Sarai Dery Folkestad
 !!    Linda Goletto, and Alexander Paul, Dec 2018
 !!
       implicit none
 !
       class(cc2) :: wf
 !
 !     Vectors sent to the routine
 !
       real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_b_j
       real(dp), dimension(wf%n_v, wf%n_o) :: rho_a_i
 !
    end subroutine construct_jacobian_cc2_A1_cc2
