   module subroutine construct_multiplier_equation_cc2(wf, equation)
!!
!!    Construct multiplier equation 
!!    Written by Sarai D. Folkestad, Feb 2019
!!
      implicit none 
!
      class(cc2), intent(in) :: wf 
!
      real(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: equation 
!
   end subroutine construct_multiplier_equation_cc2