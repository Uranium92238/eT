   module subroutine jacobian_transform_trial_vector_cc3(wf, omega, c_i)
!!
!!   Jacobian transform trail vector (CC3)
!!   Alexander Paul and Rolf H. Myhre February 2019
!!
     implicit none
!
     class(cc3), intent(in) :: wf
!
     real(dp), dimension(wf%n_amplitudes, 1) :: c_i
!
   end subroutine jacobian_transform_trial_vector_cc3
!
!
   module subroutine jacobian_cc3_transformation_cc3(wf, omega, c)
!!
!!    Jacobian transformation (CC3)
!!    Alexander Paul and Rolf H. Myhre Feb 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), intent(in) :: omega
      real(dp), dimension(wf%n_amplitudes, 1), intent(inout) :: c
!
   end subroutine jacobian_cc3_transformation_cc3
!
!
   module subroutine jacobian_cc3_integrals_cc3(wf)
!!
!!    Construct integrals need in CC3 Jacobian and store on disk
!!    Alexander Paul and Rolf H. Myhre February 2019
!!
      implicit none
!
      class(cc3) :: wf
!
   end subroutine jacobian_cc3_integrals_cc3
!
!
   module subroutine jacobian_cc3_vvv_reader_cc3(wf,batch_x,g_bdcx,g_dbxc)
!!
!!    Read in the intgrals needed in the current batches
!!    Alexander Paul and Rolf H. Myhre February 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      type(batching_index), intent(in) :: batch_x
!
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_v,batch_x%length), intent(out) :: g_bdcx
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_v,batch_x%length), intent(out) :: g_dbxc
!
   end subroutine jacobian_cc3_vvv_reader_cc3
!
!
   module subroutine jacobian_cc3_ov_vv_reader_cc3(wf,batch_y,batch_x,g_lycx,g_ylxc,L_ybxc)
!!
!!    Read the ljck, jlkc and jbkc integrals needed in the current batches
!!    Alexander Paul and Rolf H. Myhre Feb 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      type(batching_index), intent(in) :: batch_x, batch_y
!
      real(dp), dimension(wf%n_o,wf%n_v,batch_y%length,batch_x%length), intent(out) :: g_lycx
      real(dp), dimension(wf%n_v,wf%n_o,batch_y%length,batch_x%length), intent(out) :: g_ylxc
      real(dp), dimension(wf%n_v,wf%n_v,batch_y%length,batch_x%length), intent(out) :: L_ybxc
!
   end subroutine jacobian_cc3_ov_vv_reader_cc3
!
!
   module subroutine jacobian_cc3_W_calc_cc3(wf, i, j, k, t_abc, u_abc, t_abji, &
                                                g_bdci, g_bdcj, g_bdck, &
                                                g_ljci, g_lkci, g_lkcj, g_licj, g_lick, g_ljck)
!!
!!    Calculate intermediates W^abc_ijk needed for the T^abc_ijk and C^abc_ijk amplitudes
!!    Alexander Paul and Rolf H. Myhre Feb 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      integer(i15), intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: t_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: u_abc
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in)   :: t_abji
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)           :: g_bdci
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)           :: g_bdcj
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)           :: g_bdck
!
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                   :: g_ljci
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                   :: g_lkci
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                   :: g_lkcj
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                   :: g_licj
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                   :: g_lick
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                   :: g_ljck
!
!
   end subroutine jacobian_cc3_W_calc_cc3
!
!
   module subroutine jacobian_cc3_eps_cc3(wf, omega, i, j, k, t_abc)
!!
!!    Divide W^abc_ijk with -epsilon^abc_ijk to obtain T^abc_ijk
!!    Alexander Paul and Rolf H. Myhre Feb 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      integer(i15), intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_v), intent(inout) :: t_abc
!
!
   end subroutine jacobian_cc3_eps_cc3
!
!
   module subroutine jacobian_cc3_rho1_cc3(wf, i, j, k, t_abc, u_abc, rho1, rho2, F_kc, &
                                          L_jbic, L_kbic, L_kbjc, L_ibjc, L_ibkc, L_jbkc)
!!
!!    Calculate the triples contribution to rho1 and
!!    the Fock contribution to rho2
!!    Alexander Paul and Rolf H. Myhre Feb 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      integer(i15), intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)              :: t_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)             :: u_abc
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)                   :: rho1
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout)   :: rho2
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                      :: F_kc
!
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_jbic
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_kbic
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_kbjc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_ibjc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_ibkc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_jbkc
!
   end subroutine jacobian_cc3_rho1_cc3
!
!
   module subroutine jacobian_cc3_rho2_cc3(wf, i, j, k, t_abc, u_abc, v_abc, rho2, &
                                          g_dbic, g_dbjc, g_dbkc, &
                                          g_jlic, g_klic, g_kljc, g_iljc, g_ilkc, g_jlkc)
!!
!!    Calculate the triples contribution to rho1 and
!!    the Fock contribution to rho2
!!    Alexander Paul and Rolf H. Myhre Feb 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      integer(i15), intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)              :: t_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)             :: u_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)             :: v_abc
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout)   :: rho2
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)              :: g_dbic
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)              :: g_dbjc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)              :: g_dbkc
!
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: g_jlic
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: g_klic
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: g_kljc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: g_iljc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: g_ilkc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: g_jlkc
!
   end subroutine jacobian_cc3_rho2_cc3
!
