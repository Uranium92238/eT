   module subroutine effective_jacobian_transformation_cc3(wf, omega, c)
!!
!!    Jacobian transformation (CC3)
!!    Alexander Paul and Rolf H. Myhre Feb 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: c
!
   end subroutine effective_jacobian_transformation_cc3
!
!
   module subroutine jacobian_cc3_A_cc3(wf, omega, c_ai, c_abji, rho_ai, rho_abij)
!!
!!    CC3 jacobian terms
!!    Alex C. Paul and Rolf H. Myhre, Feb 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: c_abji
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout) :: rho_abij
!
   end subroutine jacobian_cc3_A_cc3
!
!
   module subroutine jacobian_cc3_c1_integrals_cc3(wf, c_ai)
!!
!!    Construct c1 transformed integrals needed in CC3 Jacobian
!!    Alexander Paul and Rolf H. Myhre February 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
!
   end subroutine jacobian_cc3_c1_integrals_cc3
!
!
   module subroutine jacobian_cc3_construct_fock_ia_c1_cc3(wf, c_ai, F_ia_c1)
!!
!!    Calculates C1 transformed elements of the Fock matrix required for the CC3 jacobian
!!    Rolf H. Myhre and Alexander Paul, February 2019
!!
!!    F_ia_c1 = sum_j L_iajj' = sum_j 2 g_iajj' - g_ij'ja
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(out) :: F_ia_c1
!
   end subroutine jacobian_cc3_construct_fock_ia_c1_cc3
!
!
   module subroutine jacobian_cc3_c3_vvv_reader_cc3(wf, batch_x, g_bdcx, g_dbxc, g_bdcx_c1)
!!
!!    Read the bdck, dbkc and c1-transformed bdck integrals in the current batch
!!    needed for the c3-contributions
!!
!!    Rolf H. Myhre and Alexander Paul, Feb 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      type(batching_index), intent(in) :: batch_x
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_bdcx
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_dbxc
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_bdcx_c1
!
   end subroutine jacobian_cc3_c3_vvv_reader_cc3
!
!
   module subroutine jacobian_cc3_t3_vvv_reader_cc3(wf, batch_x, g_bdcx, g_dbxc_c1)
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
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_bdcx
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_dbxc_c1
!
   end subroutine jacobian_cc3_t3_vvv_reader_cc3
!
!
   module subroutine jacobian_cc3_dbic_reader_cc3(wf, g_dbxc_c1)
!!
!!    Read the c1-transformed dbkc integral needed for the t3-contribution (non batching)
!!
!!    Rolf H. Myhre and Alexander Paul, Feb 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_dbxc_c1
!
   end subroutine jacobian_cc3_dbic_reader_cc3
!
!
   module subroutine jacobian_cc3_c3_ov_vv_reader_cc3(wf, batch_y, batch_x, g_lycx, g_ylxc, L_ybxc, &
                                                      g_lycx_c1)
!!
!!    Read the ljck, jlkc, jbkc and c1-transformed ljck integrals in the current batches
!!    needed for the c3-contribution
!!
!!    Rolf H. Myhre and Alexander Paul, Feb 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      type(batching_index), intent(in) :: batch_x, batch_y
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_lycx
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_ylxc
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: L_ybxc
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_lycx_c1
!
   end subroutine jacobian_cc3_c3_ov_vv_reader_cc3
!
!
   module subroutine jacobian_cc3_t3_ov_vv_reader_cc3(wf, batch_y, batch_x, g_lycx, g_ylxc_c1)
!!
!!    Read the ljck and c1-transformed jlkc integrals in the current batch
!!    needed for the t3-contribution
!!
!!    Rolf H. Myhre and Alexander Paul, Feb 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      type(batching_index), intent(in) :: batch_x, batch_y
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_lycx
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_ylxc_c1
!
   end subroutine jacobian_cc3_t3_ov_vv_reader_cc3
!
!
   module subroutine jacobian_cc3_jlkc_reader_cc3(wf, g_ylxc_c1)
!!
!!    Read the c1-transformed jlkc  needed for the t3-contribution (non batching)
!!
!!    Rolf H. Myhre and Alexander Paul, Feb 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_ylxc_c1
!
   end subroutine jacobian_cc3_jlkc_reader_cc3
!
!
   module subroutine jacobian_cc3_c3_calc_cc3(wf, omega, i, j, k, c_abc, u_abc, t_abji, c_abji,    &
                                                g_bdci, g_bdcj, g_bdck, g_ljci, g_lkci,            &
                                                g_lkcj, g_licj, g_lick, g_ljck,                    &
                                                g_bdci_c1, g_bdcj_c1, g_bdck_c1, g_ljci_c1,        &
                                                g_lkci_c1, g_lkcj_c1, g_licj_c1, g_lick_c1, g_ljck_c1)
!!
!!    Construct c^abc_ijk amplitudes
!!
!!    c^abc = (omega - ε^abc_ijk)^-1 * P^abc_ijk (sum_d c^ad_ij g_ckbd - sum_l c^ab_il g_cklj
!!             + sum_d t^ad_ij g'_bdck - sum_l t^ab_il g'_cklj
!!
!!    Alexander Paul and Rolf H. Myhre Feb 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), intent(in) :: omega
!
      integer, intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: c_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: u_abc
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in)   :: t_abji
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in)   :: c_abji
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
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)           :: g_bdci_c1
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)           :: g_bdcj_c1
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)           :: g_bdck_c1
!
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                   :: g_ljci_c1
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                   :: g_lkci_c1
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                   :: g_lkcj_c1
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                   :: g_licj_c1
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                   :: g_lick_c1
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                   :: g_ljck_c1
!
!
   end subroutine jacobian_cc3_c3_calc_cc3
!
!
   module subroutine jacobian_cc3_t3_calc_cc3(wf, i, j, k, t_abc, u_abc, t_abji,   &
                                                g_bdci, g_bdcj, g_bdck, g_ljci, g_lkci,   &
                                                g_lkcj, g_licj, g_lick, g_ljck)
!!
!!    Construct t^abc_ijk amplitudes
!!
!!    t^abc = (omega - ε^abc_ijk)^-1 * P^abc_ijk (sum_d t^ad_ij (bd|ck) - sum_l t^ab_il (ck|lj))
!!
!!    Alexander Paul and Rolf H. Myhre Feb 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: i, j, k
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
   end subroutine jacobian_cc3_t3_calc_cc3
!
!
   module subroutine jacobian_cc3_fock_rho2_cc3(wf, i, j, k, t_abc, u_abc, rho_abij, F_kc)
!!
!!    Calculate the Fock contribution to rho2 for fixed i,j and k
!!
!!    rho_abji =+ sum_kc (C^abc_ijk - C^cba_ijk) F_kc
!!
!!    Alexander Paul and Rolf H. Myhre, Feb 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)              :: t_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)             :: u_abc
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout)   :: rho_abij
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                      :: F_kc
!
   end subroutine jacobian_cc3_fock_rho2_cc3
!
!
