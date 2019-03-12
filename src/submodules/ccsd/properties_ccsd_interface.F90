
   module subroutine construct_etaX_ccsd(wf, Xoperator, etaX)
!!
!!    Construct etaX
!!    Written by Josefine H. Andersen, Februrary 2019
!!
      implicit none
!
      class(ccsd), intent(In) :: wf
!
      character(len=*), intent(in)            :: Xoperator
      real(dp), dimension(wf%n_es_amplitudes, 1), intent(inout) :: etaX
!
   end subroutine construct_etaX_ccsd
!
!
   module subroutine construct_etaX_ccs_singles_ccsd(wf, Xoperator, etaX)
!!
!!    Construct CCS term of CCSD etaX singles
!!    Written by Josefine H. Andersen, February 2019
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      character(len=*), intent(in) :: Xoperator
!
      real(dp), dimension(wf%n_es_amplitudes, 1), intent(inout) :: etaX
!
   end subroutine construct_etaX_ccs_singles_ccsd
!
!
   module subroutine construct_etaX_singles_q1_ccsd(wf, Xoperator, etaX)
!!
!!    Construct Q1 term in etaX singles contribution
!!    Written by Josefine H. Andersen, February 2019
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      character(len=*), intent(in) :: Xoperator
      !real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: etaX
      real(dp), dimension(wf%n_es_amplitudes, 1), intent(inout) :: etaX
!
   end subroutine construct_etaX_singles_q1_ccsd
!
!
   module subroutine construct_etaX_singles_q2_ccsd(wf, Xoperator, etaX)
!!
!!    Construct Q2 term in etaX singles contribution
!!    Written by Josefine H. Andersen, February 2019
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      character(len=*), intent(in) :: Xoperator
      !real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: etaX
      real(dp), dimension(wf%n_es_amplitudes, 1), intent(inout) :: etaX
!
   end subroutine construct_etaX_singles_q2_ccsd
!
!
   module subroutine construct_etaX_doubles_q1_ccsd(wf, Xoperator, etaX)
!!
!!    Construct Q1 term in etaX doubles contribution
!!    Written by Josefine H. Andersen, February 2019
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      character(len=*), intent(in)            :: Xoperator
      real(dp), dimension(wf%n_es_amplitudes, 1), intent(inout) :: etaX
!
   end subroutine construct_etaX_doubles_q1_ccsd 
!
!
   module subroutine construct_etaX_doubles_q2_ccsd(wf, Xoperator, etaX)
!!
!!    Construct Q2 term in etaX doubles contribution
!!    Written by Josefine H. Andersen, February 2019
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      character(len=*), intent(in)            :: Xoperator
      real(dp), dimension(wf%n_es_amplitudes, 1), intent(inout) :: etaX
!
   end subroutine construct_etaX_doubles_q2_ccsd 
!
!
   module subroutine construct_csiX_ccsd(wf, Xoperator, csiX)
!!
!!    Construct csiX
!!    Written by Josefine H. Andersen, February 2019
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      character(len=*), intent(in) :: Xoperator
      real(dp), dimension(wf%n_es_amplitudes, 1), intent(inout) :: csiX
!
   end subroutine construct_csiX_ccsd 
!
!
   module subroutine construct_csiX_singles_ccsd(wf, Xoperator, csiX)
!!
!!    Construct csiX singles contribution
!!    Written by Josefine H. Andersen, February 2019
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      character(len=*), intent(in) :: Xoperator
      real(dp), dimension(wf%n_es_amplitudes, 1), intent(inout) :: csiX
!
   end subroutine construct_csiX_singles_ccsd
!
!
   module subroutine construct_csiX_doubles_ccsd(wf, Xoperator, csiX)
!!
!!    Construct csiX doubles contribution
!!    Written by Josefine H. Andersen, February 2019
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      character(len=*), intent(in) :: Xoperator
      real(dp), dimension(wf%n_es_amplitudes, 1), intent(inout) :: csiX
!
   end subroutine construct_csiX_doubles_ccsd
!
!
   module subroutine get_eom_contribution_ccsd(wf, etaX, csiX, Xoperator)
!!
!!    Add EOM contribution to csiX vector, CCSD
!!    Written by Josefine H. Andersen
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      character(len=*), intent(inout), optional :: Xoperator
!
      real(dp), dimension(wf%n_es_amplitudes, 1), intent(inout) :: etaX
      real(dp), dimension(wf%n_es_amplitudes, 1), intent(in)    :: csiX
!
   end subroutine get_eom_contribution_ccsd
!
!
   module subroutine get_eom_doubles_contribution_ccsd(wf, Xoperator, etaX)
!!
!!    Build EOM contribution to etaX in CCSD
!!    Written by Josefine H. Andersen, February 2019
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      character(len=*), intent(in) :: Xoperator
!      
      real(dp), dimension(wf%n_es_amplitudes, 1), intent(inout) :: etaX
!
   end subroutine get_eom_doubles_contribution_ccsd
!
!

