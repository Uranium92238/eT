
   module subroutine jacobian_transpose_transform_trial_vector_ccsd(wf, c_i)
!!
!!    Jacobian transpose transform trial vector 
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!    and Andreas Skeidsvoll, 2018
!!
      class(ccsd), intent(in) :: wf 
!
      real(dp), dimension(wf%n_es_amplitudes, 1) :: c_i
!
   end subroutine jacobian_transpose_transform_trial_vector_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_transformation_ccsd(wf, b)
!!
!!    Jacobian transpose transformation (CCSD)
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!    and Andreas Skeidsvoll, 2018
!!
      implicit none 
!
      class(ccsd) :: wf 
      real(dp), dimension(wf%n_es_amplitudes, 1) :: b
!
   end subroutine jacobian_transpose_ccsd_transformation_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_a1_ccsd(wf, sigma_a_i, b_a_i)
!!
!!    Jacobian transpose CCSD A1 
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad
!!    and Andreas Skeidsvoll, 2018
!!
      implicit none 
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: b_a_i 
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_a_i 
!
   end subroutine jacobian_transpose_ccsd_a1_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_b1_ccsd(wf, sigma_a_i, b_a_i)
!!
!!    Jacobian transpose CCSD B1 
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad
!!    and Andreas Skeidsvoll, 2018
!!
      implicit none 
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: b_a_i 
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_a_i 
!
   end subroutine jacobian_transpose_ccsd_b1_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_c1_ccsd(wf, sigma_a_i, b_ai_bj)
!!
!!    Jacobian transpose CCSD C1 
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad
!!    and Andreas Skeidsvoll, 2018
!!
      implicit none 
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o)                       :: sigma_a_i 
      real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: b_ai_bj 
!
   end subroutine jacobian_transpose_ccsd_c1_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_d1_ccsd(wf, sigma_a_i, b_ai_bj)
!!
!!    Jacobian transpose CCSD D1 
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad
!!    and Andreas Skeidsvoll, 2018
!!
      implicit none 
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o)                       :: sigma_a_i 
      real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: b_ai_bj 
!
   end subroutine jacobian_transpose_ccsd_d1_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_e1_ccsd(wf, sigma_a_i, b_ai_bj)
!!
!!    Jacobian transpose CCSD E1 
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad
!!    and Andreas Skeidsvoll, 2018
!!
      implicit none 
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o)                       :: sigma_a_i 
      real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: b_ai_bj 
!
   end subroutine jacobian_transpose_ccsd_e1_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_f1_ccsd(wf, sigma_a_i, b_ai_bj)
!!
!!    Jacobian transpose CCSD F1 
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad
!!    and Andreas Skeidsvoll, 2018
!!
      implicit none 
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o)                       :: sigma_a_i 
      real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: b_ai_bj 
!
   end subroutine jacobian_transpose_ccsd_f1_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_g1_ccsd(wf, sigma_a_i, b_ai_bj)
!!
!!    Jacobian transpose CCSD G1 
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad
!!    and Andreas Skeidsvoll, 2018
!!
      implicit none 
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o)                       :: sigma_a_i 
      real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: b_ai_bj
!
   end subroutine jacobian_transpose_ccsd_g1_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_a2_ccsd(wf, sigma_ai_bj, b_a_i)
!!
!!    Jacobian transpose CCSD A2 
!!    Written by Sarai D. Folkestad Eirik F. Kjønstad
!!    and Andreas Skeidsvoll, 2018
!!
      implicit none 
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o)                       :: b_a_i  
      real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: sigma_ai_bj 
!
   end subroutine jacobian_transpose_ccsd_a2_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_b2_ccsd(wf, sigma_ai_bj, b_ai_bj)
!!
!!    Jacobian transpose CCSD B2 
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad
!!    and Andreas Skeidsvoll, 2018
!!
      implicit none 
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: b_ai_bj 
      real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: sigma_ai_bj 
!
   end subroutine jacobian_transpose_ccsd_b2_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_c2_ccsd(wf, sigma_ai_bj, b_ai_bj)
!!
!!    Jacobian transpose CCSD C2 
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad
!!    and Andreas Skeidsvoll
!!
      implicit none 
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: b_ai_bj 
      real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: sigma_ai_bj 
!
   end subroutine jacobian_transpose_ccsd_c2_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_d2_ccsd(wf, sigma_ai_bj, b_ai_bj)
!!
!!    Jacobian transpose CCSD D2 
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad
!!    and Andreas Skeidsvoll, 2018
!!
      implicit none 
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: b_ai_bj 
      real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: sigma_ai_bj
!
   end subroutine jacobian_transpose_ccsd_d2_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_e2_ccsd(wf, sigma_ai_bj, b_ai_bj)
!!
!!    Jacobian transpose CCSD E2 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
      implicit none 
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: b_ai_bj 
      real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: sigma_ai_bj
!
   end subroutine jacobian_transpose_ccsd_e2_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_f2_ccsd(wf, sigma_ai_bj, b_ai_bj)
!!
!!    Jacobian transpose CCSD F2 
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad
!!    and Andreas Skeidsvoll, 2018
!!
      implicit none 
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: b_ai_bj 
      real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: sigma_ai_bj
!
   end subroutine jacobian_transpose_ccsd_f2_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_g2_ccsd(wf, sigma_ai_bj, b_ai_bj)
!!
!!    Jacobian transpose CCSD G2 
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad
!!    and Andreas Skeidsvoll, 2018
!!
      implicit none 
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: b_ai_bj 
      real(dp), dimension((wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o)) :: sigma_ai_bj
!
   end subroutine jacobian_transpose_ccsd_g2_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_h2_ccsd(wf, sigma_ab_ij, b_ab_ij)
!!
!!    Jacobian transpose CCSD H2 
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad
!!    and Andreas Skeidsvoll, 2018
!!
      implicit none 
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_v)**2, (wf%n_o)**2) :: b_ab_ij
      real(dp), dimension((wf%n_v)**2, (wf%n_o)**2) :: sigma_ab_ij
!
   end subroutine jacobian_transpose_ccsd_h2_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_i2_ccsd(wf, sigma_ab_ij, b_ab_ij)
!!
!!    Jacobian transpose CCSD I2 
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad
!!    and Andreas Skeidsvoll, 2018
!!
      implicit none 
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_v)**2, (wf%n_o)**2) :: b_ab_ij
      real(dp), dimension((wf%n_v)**2, (wf%n_o)**2) :: sigma_ab_ij
!
   end subroutine jacobian_transpose_ccsd_i2_ccsd
