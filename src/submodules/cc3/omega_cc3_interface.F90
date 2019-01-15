!
   module subroutine construct_omega_cc3(wf, omega)
!!
!!   Construct omega (CC3)
!!   Alex C. Paul and Rolf H. Myhre 2018
!!
     implicit none
!
     class(cc3), intent(in) :: wf
!
     real(dp), dimension(wf%n_amplitudes, 1), intent(inout) :: omega
!
   end subroutine construct_omega_cc3
!
!
   module subroutine omega_cc3_a_cc3(wf, omega1, omega2)
!!
!!    CC3 Omega terms
!!    Alex C. Paul and Rolf H. Myhre 2018
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: omega1
      real(dp), dimension((wf%n_v)*(wf%n_o)*((wf%n_v)*(wf%n_o) +1)/2), intent(inout) :: omega2
!
   end subroutine omega_cc3_a_cc3
!
!
   module subroutine omega_cc3_integrals_cc3(wf)
!!
!!    Construct integrals need in CC3 Omega and store on disk
!!    Rolf H. Myhre, January 2019
!!
      implicit none
!
      class(cc3) :: wf
!
   end subroutine omega_cc3_integrals_cc3
!
!
   module subroutine omega_cc3_vvv_reader_cc3(wf,batch_x,g_bdcx,g_dbxc)
!!
!!    Read in the intgrals needed in the current batches
!!
!!    Rolf H. Myhre, January 2019
!!
      implicit none
!
      real(dp), dimension(:,:,:,:), intent(inout) :: g_bdcx
      real(dp), dimension(:,:,:,:), intent(inout) :: g_dbxc
!
      type(batching_index), intent(in) :: batch_x
!
      class(cc3) :: wf
!
   end subroutine omega_cc3_vvv_reader_cc3
!
!
   module subroutine omega_cc3_ov_vv_reader_cc3(wf,batch_y,batch_x,g_lycx,g_ylxc,L_ybxc)
!!
!!    Read the ljck, jlkc and jbkc integrals needed in the current batches
!!
!!    Rolf H. Myhre, January 2019
!!
      implicit none
!
      real(dp), dimension(:,:,:,:), intent(inout) :: g_lycx
      real(dp), dimension(:,:,:,:), intent(inout) :: g_ylxc
      real(dp), dimension(:,:,:,:), intent(inout) :: L_ybxc
!
      type(batching_index), intent(in) :: batch_x, batch_y
!
      class(cc3) :: wf
!
   end subroutine omega_cc3_ov_vv_reader_cc3
!
!
