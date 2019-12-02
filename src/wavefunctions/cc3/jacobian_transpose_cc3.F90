!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2019 the authors of eT
!
!  eT is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  eT is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program. If not, see <https://www.gnu.org/licenses/>.
!
!
submodule (cc3_class) jacobian_transpose
!
!!
!!    Jacobian transpose submodule (cc3)
!!    Written by Alexander C. Paul and Rolf H. Myhre, March 2019
!!
!!    Routines for the linear transform of trial
!!    vectors by the transpose of the Jacobian matrix
!!
!!    sigma_i = A^T * c_i,
!!
!!    where
!!
!!    A_mu,nu = < mu| exp(-T) [H, tau_nu] exp(T)|R >.
!!
!
   implicit none
!
!
contains
!
!
   module subroutine effective_jacobian_transpose_transformation_cc3(wf, omega, c, cvs)
!!
!!    Effective Jacobian transpose transformation (CC3)
!!    Alexander C. Paul and Rolf H. Myhre, March 2019
!!
!!    Directs the transformation by the transpose of the  CC3 Jacobi matrix,
!!
!!       A_mu,nu = < mu| exp(-T) [H, tau_nu] exp(T)|R >,
!!
!!    The transformation is performed as sigma^T = c^T A, where c is the vector
!!    sent to the routine. On exit, the vector c is equal to sigma (the transformed
!!    vector).
!!
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: c
!
      logical, intent(in) :: cvs
!
      real(dp), dimension(:,:), allocatable :: c_ai
      real(dp), dimension(:,:,:,:), allocatable :: c_aibj, c_abij
!
      real(dp), dimension(:,:), allocatable :: sigma_ai
      real(dp), dimension(:,:,:,:), allocatable :: sigma_aibj, sigma_abij
!
      type(timings) :: cc3_timer, cc3_timer_t3_a1, cc3_timer_t3_b1, cc3_timer_c3
      type(timings) :: ccsd_timer
!
      cc3_timer_t3_a1 = timings('Time in CC3 T3 a1')
      cc3_timer_t3_b1 = timings('Time in CC3 T3 b1')
      cc3_timer_c3   = timings('Time in CC3 C3')
      cc3_timer      = timings('Total CC3 contribution')
      ccsd_timer     = timings('Total CCSD contribution')
!
!     Allocate and zero the transformed singles vector
!
      call mem%alloc(sigma_ai, wf%n_v, wf%n_o)
      call zero_array(sigma_ai, wf%n_v*wf%n_o)
!
      call mem%alloc(c_ai, wf%n_v, wf%n_o)
!
      call dcopy(wf%n_t1, c, 1, c_ai, 1)
!
!     :: CCS contributions to the transformed singles vector ::
!
      call ccsd_timer%turn_on()
!
      call wf%jacobian_transpose_ccs_a1(sigma_ai, c_ai)
      call wf%jacobian_transpose_ccs_b1(sigma_ai, c_ai)
!
!     :: CCSD contributions to the transformed singles vector ::
!
      call wf%jacobian_transpose_doubles_a1(sigma_ai, c_ai, wf%u_aibj)
!
!     Allocate the incoming unpacked doubles vector
!
      call mem%alloc(c_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call squareup(c(wf%n_t1 + 1 : wf%n_es_amplitudes), c_aibj, wf%n_t1)
!
      call wf%jacobian_transpose_doubles_b1(sigma_ai, c_aibj)
      call wf%jacobian_transpose_ccsd_d1(sigma_ai, c_aibj)
      call wf%jacobian_transpose_ccsd_e1(sigma_ai, c_aibj)
      call wf%jacobian_transpose_ccsd_f1(sigma_ai, c_aibj)
      call wf%jacobian_transpose_ccsd_g1(sigma_ai, c_aibj)
!
!     :: CCSD contributions to the transformed doubles vector ::
!     Allocate unpacked transformed vector
!
      call mem%alloc(sigma_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call zero_array(sigma_aibj, (wf%n_v*wf%n_o)**2)
!
!     Contributions from singles vector c
!
      !call wf%jacobian_transpose_ccsd_a2(sigma_aibj, c_ai)
      call wf%jacobian_transpose_doubles_a2(sigma_aibj, c_ai)
!
!     Contributions from doubles vector c
!
      call wf%jacobian_transpose_ccsd_b2(sigma_aibj, c_aibj)
      call wf%jacobian_transpose_ccsd_c2(sigma_aibj, c_aibj)
      call wf%jacobian_transpose_ccsd_d2(sigma_aibj, c_aibj)
      call wf%jacobian_transpose_ccsd_e2(sigma_aibj, c_aibj)
      call wf%jacobian_transpose_ccsd_f2(sigma_aibj, c_aibj)
      call wf%jacobian_transpose_ccsd_g2(sigma_aibj, c_aibj)
!
!     Compute CC3 contributions to sigma_ai and sigma_aibj and symmetrize sigma_aibj
!     CCSD H2 and I2 are already symmetric and will be computed afterwards
!
      call ccsd_timer%freeze()
!
      call mem%alloc(sigma_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(c_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call sort_1234_to_1324(sigma_aibj, sigma_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_1324(c_aibj, c_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(sigma_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(c_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call cc3_timer%turn_on()
!
!     CC3-Contributions from the T3-amplitudes
      call cc3_timer_t3_a1%turn_on()
!
      call wf%jacobian_transpose_cc3_t3_a1(c_abij, sigma_ai)
!
      call cc3_timer_t3_a1%turn_off()
!
      call cc3_timer_t3_b1%turn_on()
!
      call wf%jacobian_transpose_cc3_t3_b1(c_abij, sigma_ai, cvs)
!
      call cc3_timer_t3_b1%turn_off()
!
!     CC3-Contributions from the C3-amplitudes
      call cc3_timer_c3%turn_on()
!
      call wf%jacobian_transpose_cc3_c3_a(omega, c_ai, c_abij, sigma_ai, sigma_abij, cvs)
!
      call cc3_timer_c3%turn_off()
!
      call cc3_timer%turn_off()
!
!     Done with singles vector c; Overwrite the incoming singles c vector for exit
!
      call ccsd_timer%turn_on()
!
      call mem%dealloc(c_ai, wf%n_v, wf%n_o)
!
      call dcopy(wf%n_t1, sigma_ai, 1, c, 1)
!
      call mem%dealloc(sigma_ai, wf%n_v, wf%n_o)
!
!     Last two CCSD-terms (H2, I2) are already symmetric.
!     Perform the symmetrization sigma_aibj = P_ij^ab sigma_aibj
!
      call symmetrize_12_and_34(sigma_abij, wf%n_v, wf%n_o)
!
!     Compute CCSD H2 and I2 contributions
!
      call wf%jacobian_transpose_ccsd_h2(sigma_abij, c_abij)
      call wf%jacobian_transpose_ccsd_i2(sigma_abij, c_abij)
!
      call mem%dealloc(c_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call ccsd_timer%turn_off()
!
!     overwrite the incoming, packed doubles c vector with sigma_abij
!     and pack in for exit
!
      call packin(c(wf%n_t1 + 1 : wf%n_es_amplitudes), sigma_abij, wf%n_v, wf%n_o)
!
      call mem%dealloc(sigma_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
   end subroutine effective_jacobian_transpose_transformation_cc3
!
!
   module subroutine jacobian_transpose_cc3_t3_a1_cc3(wf, c_abij, sigma_ai)
!!
!!    Jacobian transpose T3 A1 term
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!!    Computes the first contribution of the T3 amplitudes to sigma_1
!!
!!    Reads in the intermediates X_abid and X_ajil prepared in prepare_jacobian_transpose
!!    contracts with c_abij and adds to sigma_ai
!!
!!    sigma_dl +=  sum_abi X_abid * C_abil + sum_aik C_daji * X_ajil
!!
!!    where: X_abid = - sum_jck (2t^abc_ijk - t^cba_ijk - t^acb_ijk) * g_kcjd
!!           X_ajil = - sum_bck (2t^abc_ijk - t^cba_ijk - t^acb_ijk) * g_lbkc
!!    
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: c_abij
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(out) :: sigma_ai
!
      real(dp), dimension(:,:,:,:), allocatable :: X_abid
      real(dp), dimension(:,:,:,:), allocatable :: X_ajil
!
      type(batching_index) :: batch_d
      integer :: d_batch
      integer :: req_0, req_d
!
!     :: X_abid term ::
!
      batch_d = batching_index(wf%n_v)
!
      req_0 = 0
      req_d = wf%n_o * wf%n_v**2
!
      call mem%batch_setup(batch_d, req_0, req_d)
!
      call wf%X_abid%open_('read')
!
      call batch_d%determine_limits(1)
      call mem%alloc(X_abid, wf%n_v, wf%n_v, wf%n_o, batch_d%length)
!
      do d_batch = 1, batch_d%num_batches
!
         call batch_d%determine_limits(d_batch)
!
!        Read in X_abid written with compound index "id" as record
!
         call compound_record_reader(wf%n_o, batch_d, wf%X_abid, X_abid)
!
         call dgemm('T','N',                    & ! X is transposed
                     batch_d%length,            &
                     wf%n_o,                    &
                     wf%n_o*wf%n_v**2,          &
                     one,                       &
                     X_abid,                    & ! X_d_abi
                     wf%n_o*wf%n_v**2,          &
                     c_abij,                    & ! c_abi_l
                     wf%n_o*wf%n_v**2,          &
                     one,                       &
                     sigma_ai(batch_d%first,1), & ! sigma_dl
                     wf%n_v)
!
      enddo
!
      call batch_d%determine_limits(1)
      call mem%dealloc(X_abid, wf%n_v, wf%n_v, wf%n_o, batch_d%length)
!
      call wf%X_abid%close_()
!
!     :: X_ajil term ::
!
      call wf%X_ajil%open_('read')
!
      call mem%alloc(X_ajil, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call single_record_reader(wf%n_o, wf%X_ajil, X_ajil)
!
      call wf%X_ajil%close_()
!
      call dgemm('N','N',              &
                  wf%n_v,              &
                  wf%n_o,              &
                  wf%n_v * wf%n_o**2,  &
                  one,                 &
                  c_abij,              & ! C_d_aji
                  wf%n_v,              &
                  X_ajil,              & ! X_aji_l
                  wf%n_v * wf%n_o**2,  &
                  one,                 &
                  sigma_ai,            & ! sigma_dl
                  wf%n_v)
!
      call mem%dealloc(X_ajil, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
   end subroutine jacobian_transpose_cc3_t3_a1_cc3
!
!
   module subroutine jacobian_transpose_cc3_t3_b1_cc3(wf, c_abij, sigma_ai, cvs)
!!
!!    Jacobian transpose T3 B1 term
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!!    Constructs t^abc_ijk for fixed ijk and contracts with c_abij
!!    The intermediate X_ai is then contracted with L_iald
!!
!!    sigma_dl +=  sum_abcijk C^bc_jk (t^abc_ijk - t^bac_ijk) L_iald
!!             +=  sum_ai X_ai * L_iald
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: c_abij
!
      logical, intent(in) :: cvs
!
!     Arrays for triples amplitudes
      real(dp), dimension(:,:,:), allocatable :: t_abc
      real(dp), dimension(:,:,:), allocatable :: u_abc
!
!     Unpacked doubles amplitudes
      real(dp), dimension(:,:,:,:), allocatable :: t_abij
!
!     Integrals and Pointers
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdci
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdcj
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdck
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdci_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdcj_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdck_p => null()
!
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ljci
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_lkci
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_lkcj
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_licj
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_lick
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ljck
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_ljci_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_lkci_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_lkcj_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_licj_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_lick_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_ljck_p => null()
!
      real(dp), dimension(:,:,:,:), allocatable :: L_iald
      real(dp), dimension(:,:,:,:), allocatable :: L_jbkc
!
!     Intermediate
      real(dp), dimension(:,:), allocatable :: X_ai
!
      type(batching_index) :: batch_i, batch_j, batch_k, batch_l
      integer :: i, j, k, i_rel, j_rel, k_rel
      integer :: i_batch, j_batch, k_batch, l_batch ! used for the current batch
      integer :: req_0, req_1, req_2, req_3
!
      logical :: ijk_core
      integer :: i_cvs
!
!     :: Construct intermediate X_ai ::
!
      call mem%alloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call squareup_and_sort_1234_to_1324(wf%t2, t_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Setup and Batching loops
!
      req_0 = 2*wf%n_v**3 + wf%n_v*wf%n_o
      req_1 = wf%n_v**3
      req_2 = wf%n_o*wf%n_v
      req_3 = 0
!
      batch_i = batching_index(wf%n_o)
      batch_j = batching_index(wf%n_o)
      batch_k = batching_index(wf%n_o)
!
      call mem%batch_setup_ident(batch_i, batch_j, batch_k, &
                                 req_0, req_1, req_2, req_3, buffer_size = zero)
!
!     Allocate integral arrays
!
!     Split up so that the integral and amplitude arrays are closer in mem
!
      if (batch_i%num_batches .eq. 1) then ! no batching
!
         call mem%alloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%alloc(g_ljci, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      else ! batching
!
         call batch_i%determine_limits(1)
!
         call mem%alloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%alloc(g_bdcj, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%alloc(g_bdck, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
         call mem%alloc(g_ljci, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_lkci, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_lkcj, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_licj, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_lick, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_ljck, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
!
      endif
!
!     Arrays for the triples amplitudes
      call mem%alloc(t_abc, wf%n_v, wf%n_v, wf%n_v)
      call mem%alloc(u_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     Array for the whole intermediate X_ai
      call mem%alloc(X_ai, wf%n_v, wf%n_o)
      call zero_array(X_ai, wf%n_v*wf%n_o)
!
      call wf%g_bdck_t%open_('read')
      call wf%g_ljck_t%open_('read')
!
      do i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(i_batch)
!
         call single_record_reader(batch_i, wf%g_bdck_t, g_bdci)
         g_bdci_p => g_bdci
!
         do j_batch = 1, i_batch
!
            call batch_j%determine_limits(j_batch)
!
            call compound_record_reader(batch_j, batch_i, wf%g_ljck_t, g_ljci)
            g_ljci_p => g_ljci
!
            if (j_batch .ne. i_batch) then
!
               call single_record_reader(batch_j, wf%g_bdck_t, g_bdcj)
               g_bdcj_p => g_bdcj
!
               call compound_record_reader(batch_i, batch_j, wf%g_ljck_t, g_licj)
               g_licj_p => g_licj
!
            else
!
               g_bdcj_p => g_bdci
!
               g_licj_p => g_ljci
!
            endif
!
            do k_batch = 1, j_batch
!
               call batch_k%determine_limits(k_batch)
!
               if (k_batch .ne. j_batch) then ! k_batch != j_batch, k_batch != i_batch
!
                  call single_record_reader(batch_k, wf%g_bdck_t, g_bdck)
                  g_bdck_p => g_bdck
! 
                  call compound_record_reader(batch_k, batch_i, wf%g_ljck_t, g_lkci)
                  g_lkci_p => g_lkci
!
                  call compound_record_reader(batch_i, batch_k, wf%g_ljck_t, g_lick)
                  g_lick_p => g_lick
!
                  call compound_record_reader(batch_k, batch_j, wf%g_ljck_t, g_lkcj)
                  g_lkcj_p => g_lkcj
!
                  call compound_record_reader(batch_j, batch_k, wf%g_ljck_t, g_ljck)
                  g_ljck_p => g_ljck
!
!
               else if (k_batch .eq. i_batch) then ! k_batch == j_batch == i_batch
!
                  g_bdck_p => g_bdci
!
                  g_lkci_p => g_ljci
                  g_lick_p => g_ljci
!
                  g_lkcj_p => g_ljci
                  g_ljck_p => g_ljci
!
               else ! k_batch == j_batch != i_batch
!
                  g_bdck_p => g_bdcj
!
                  g_lkci_p => g_ljci
                  g_lick_p => g_licj
!
                  call compound_record_reader(batch_k, batch_j, wf%g_ljck_t, g_lkcj)
                  g_lkcj_p => g_lkcj
                  g_ljck_p => g_lkcj
!
               endif
!
!              sigma_dl +=  sum_abcijk C^bc_jk (t^abc_ijk - t^bac_ijk) L_iald
!              CVS: in principle check j,k and l but due to the symmetry in L_iald
!                   we can also check i,j,k
!
               do i = batch_i%first, batch_i%last
!
                  i_rel = i - batch_i%first + 1
!
                  do j = batch_j%first, min(batch_j%last, i)
!
                     j_rel = j - batch_j%first + 1
!
                     do k = batch_k%first, min(batch_k%last, j)
!
                        if (k .eq. i) then ! k == j == i
                           cycle
                        end if
!
!                       Check if at least one index i,j,k is a core orbital
                        if(cvs) then
!
                           ijk_core = .false.
!
                           do i_cvs = 1, wf%n_core_MOs
!
                              if(     i .eq. wf%core_MOs(i_cvs)   &
                                 .or. j .eq. wf%core_MOs(i_cvs)   &
                                 .or. k .eq. wf%core_MOs(i_cvs))  then
!
                                 ijk_core = .true.
!
                              end if
!
                           end do
!
!                          Cycle if i,j,k are not core orbitals
                           if (.not. ijk_core) cycle
!
                        end if
!
                        k_rel = k - batch_k%first + 1
!
!                       Construct t^{abc}_{ijk} for given i, j, k
!
                        call wf%construct_W(i, j, k, t_abc, u_abc, t_abij, &
                                            g_bdci_p(:,:,:,i_rel),         &
                                            g_bdcj_p(:,:,:,j_rel),         &
                                            g_bdck_p(:,:,:,k_rel),         &
                                            g_ljci_p(:,:,j_rel,i_rel),     &
                                            g_lkci_p(:,:,k_rel,i_rel),     &
                                            g_lkcj_p(:,:,k_rel,j_rel),     &
                                            g_licj_p(:,:,i_rel,j_rel),     &
                                            g_lick_p(:,:,i_rel,k_rel),     &
                                            g_ljck_p(:,:,j_rel,k_rel))
!
                        call wf%divide_by_orbital_differences(i, j, k, t_abc)
!
                        call wf%construct_x_ai_intermediate(i, j, k, t_abc, u_abc, X_ai, c_abij)
!
                     enddo ! loop over k
                  enddo ! loop over j
               enddo ! loop over i
            enddo ! batch_k
         enddo ! batch_j
      enddo ! batch_i
!
!     Close files
!
      call wf%g_bdck_t%close_()
      call wf%g_ljck_t%close_()
!
!     Deallocate the integral arrays
!
      if (batch_i%num_batches .eq. 1) then ! no batching
!
         call mem%dealloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%dealloc(g_ljci, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      else ! batching
!
         call batch_i%determine_limits(1)
!
         call mem%dealloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%dealloc(g_bdcj, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%dealloc(g_bdck, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
         call mem%dealloc(g_ljci, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_lkci, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_lkcj, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_licj, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_lick, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_ljck, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
!
      endif
!
!     Deallocate amplitude arrays
!
      call mem%dealloc(t_abc, wf%n_v, wf%n_v, wf%n_v)
      call mem%dealloc(u_abc, wf%n_v, wf%n_v, wf%n_v)
!
      call mem%dealloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
!
!     :: sigma_dl += sum_ai X_ai * L_iald ::
!
!
      req_0 = 0
      req_1 = wf%n_v**2 * wf%n_o
!
      call wf%L_jbkc_t%open_('read')
!
      batch_l = batching_index(wf%n_o)
      call mem%batch_setup(batch_l, req_0, req_1)
      call batch_l%determine_limits(1)
!
      call mem%alloc(L_jbkc, wf%n_v, wf%n_v, wf%n_o, batch_l%length) ! ordered adil
      call mem%alloc(L_iald, wf%n_v, wf%n_o, wf%n_v, batch_l%length) ! ordered aidl
!
      do l_batch = 1, batch_l%num_batches
!
         call batch_l%determine_limits(l_batch)
!
         call compound_record_reader(wf%n_o, batch_l, wf%L_jbkc_t, L_jbkc)
!
         call sort_1234_to_1324(L_jbkc, L_iald, wf%n_v, wf%n_v, wf%n_o, batch_l%length)
!
!        sigma_dl += sum_ck X_ai * L_iald
!
         call dgemv('T',                        &
                     wf%n_v * wf%n_o,           &
                     wf%n_v * batch_l%length,   &
                     one,                       &
                     L_iald,                    & ! L_dl_ai
                     wf%n_v * wf%n_o,           &
                     X_ai,                      & ! X_ai
                     1,                         &
                     one,                       &
                     sigma_ai(1,batch_l%first), & ! sigma_dl
                     1)
!
      enddo
!
      call batch_l%determine_limits(1)
      call mem%dealloc(L_jbkc, wf%n_v, wf%n_v, wf%n_o, batch_l%length)
      call mem%dealloc(L_iald, wf%n_v, wf%n_o, wf%n_v, batch_l%length)
!
      call wf%L_jbkc_t%close_()
!
      call mem%dealloc(X_ai, wf%n_v, wf%n_o)
!
   end subroutine jacobian_transpose_cc3_t3_b1_cc3
!
!
   module subroutine construct_x_ai_intermediate_cc3(wf, i, j, k, t_abc, u_abc, x_ai, c_bcjk)
!!
!!    Constructs the intermediate X_ai
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!!    X_ai = sum_abcijk (t^abc_ijk - t^bac_ijk) C^bc_jk 
!!
!!    All permutations for i,j,k have to be considered due to the restrictions in the i,j,k loops
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)           :: t_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: u_abc
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(out)                  :: X_ai
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in)   :: c_bcjk
!
      
!     Construct u_acb = (t_acb - t_bca)
!     u_acb zero if i == k, but this is never true
!
      call construct_132_minus_231(t_abc, u_abc, wf%n_v)
!
!     X_ai += sum_bc (t^acb - t^bca) * C_bckj
!
      call dgemv('N',               &
                  wf%n_v,           & 
                  wf%n_v**2,        & 
                  one,              &
                  u_abc,            & ! u_a_bc
                  wf%n_v,           &
                  c_bcjk(:,:,k,j),  & ! c_bc,kj
                  1,                &
                  one,              &
                  X_ai(:,i),        & ! X_a,i
                  1)
!
!     X_ak += -sum_bc (t^acb - t^bca) * C_bcij
!
      call dgemv('N',               &
                  wf%n_v,           & 
                  wf%n_v**2,        & 
                  -one,             &
                  u_abc,            & ! u_a_bc
                  wf%n_v,           &
                  c_bcjk(:,:,i,j),  & ! c_bc,ij
                  1,                &
                  one,              &
                  X_ai(:,k),        & ! X_a,k
                  1)
!
!     Only this condition needed. If 2 indices would be equal u would be zero 
!     or the contribution would be equal to one of the first two terms
      if (k .ne. j .and. j .ne. i) then
!
!        Construct u_abc = t_abc - t_bac
!
         call construct_123_minus_213(t_abc, u_abc, wf%n_v)
!
!        X_ai += sum_bc (t^abc - t^bac) * C_bcjk
!
         call dgemv('N',               &
                     wf%n_v,           & 
                     wf%n_v**2,        & 
                     one,              &
                     u_abc,            & ! u_a_bc
                     wf%n_v,           &
                     c_bcjk(:,:,j,k),  & ! c_bc,jk
                     1,                &
                     one,              &
                     X_ai(:,i),        & ! X_a,i
                     1)
!
!        X_aj += -sum_bc (t^abc - t^bac) * C_bcik
!
         call dgemv('N',               &
                     wf%n_v,           & 
                     wf%n_v**2,        & 
                     -one,             &
                     u_abc,            & ! u_a_bc
                     wf%n_v,           &
                     c_bcjk(:,:,i,k),  & ! c_bc,ik
                     1,                &
                     one,              &
                     X_ai(:,j),        & ! X_a,j
                     1)
!
!        Construct u_cba = t_cba - t_cab
!
         call construct_321_minus_312(t_abc, u_abc, wf%n_v)
!
!        X_ak += sum_bc (t^cba - t^cab) * C_bcji
!
         call dgemv('N',               &
                     wf%n_v,           & 
                     wf%n_v**2,        & 
                     one,              &
                     u_abc,            & ! u_a_bc
                     wf%n_v,           &
                     c_bcjk(:,:,j,i),  & ! c_bc,ji
                     1,                &
                     one,              &
                     X_ai(:,k),        & ! X_a,k
                     1)
!
!        X_aj += -sum_bc (t^cba - t^cab) * C_bcki
!
         call dgemv('N',               &
                     wf%n_v,           & 
                     wf%n_v**2,        & 
                     -one,             &
                     u_abc,            & ! u_a_bc
                     wf%n_v,           &
                     c_bcjk(:,:,k,i),  & ! c_bc,ki
                     1,                &
                     one,              &
                     X_ai(:,j),        & ! X_a,j
                     1)
!
      end if
!
   end subroutine construct_x_ai_intermediate_cc3


   module subroutine jacobian_transpose_cc3_c3_a_cc3(wf, omega,            &
                                                     c_ai, c_abij,         &
                                                     sigma_ai, sigma_abij, &
                                                     cvs)
!!
!!    Contributions of the c3/L3 terms
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!!
!!    Construct C^abc_ijk in single batches of ijk and compute the contributions
!!    to the singles and doubles part of the outgoing vector
!!
!!    The construction of C3 is split into contributions 
!!    from outer products and matrix multiplications
!!
!!    1 array for each Permutation of C_abc will be used 
!!    to reduce the amount of N^7-contractions and sorting
!!
!!    c_mu3 = (omega - epsilon^abc_ijk)^-1 (c_mu1 < mu1| [H,tau_nu3] |R > 
!!                                        + c_mu2 < mu2| [H,tau_nu3] |R >)
!!
!!    sigma_1 += c_mu3 < mu3| [[H,T_2],tau_nu1] |R >
!!    sigma_2 += c_mu3 < mu3| [H,tau_ nu2] |R >
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: c_abij
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout) :: sigma_abij
!
      logical, intent(in) :: cvs
!
!     Arrays for triples amplitudes
      real(dp), dimension(:,:,:), allocatable :: c_abc
      real(dp), dimension(:,:,:), allocatable :: u_abc
      real(dp), dimension(:,:,:), allocatable :: v_abc
!
!     Unpacked doubles amplitudes
      real(dp), dimension(:,:,:,:), allocatable :: t_abij
!
      real(dp), dimension(:,:), allocatable :: F_ov_ck
!
!     Arrays for intermediates
!     cannot hold the whole Y_bcek array
      real(dp), dimension(:,:,:,:), allocatable, target :: Y_bcei
      real(dp), dimension(:,:,:,:), allocatable, target :: Y_bcej
      real(dp), dimension(:,:,:,:), allocatable, target :: Y_bcek
      real(dp), dimension(:,:,:,:), contiguous, pointer :: Y_bcei_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: Y_bcej_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: Y_bcek_p => null()
!
      real(dp), dimension(:,:,:,:), allocatable :: Y_cmjk
!
!     Integrals and Pointers
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdci
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdcj
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdck
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdci_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdcj_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdck_p => null()
!
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_dbic
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_dbjc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_dbkc
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_dbic_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_dbjc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_dbkc_p => null()
!
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ljci
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_lkci
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_lkcj
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_licj
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_lick
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ljck
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_ljci_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_lkci_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_lkcj_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_licj_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_lick_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_ljck_p => null()
!
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_jlic
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_klic
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_kljc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_iljc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ilkc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_jlkc
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_jlic_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_klic_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_kljc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_iljc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_ilkc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_jlkc_p => null()
!
!     L_kbjc and L_jbkc each other's transpose,
!     but utilising this makes the code more complicated and
!     error prone without any huge advantages
!
      real(dp), dimension(:,:,:,:), allocatable, target  :: L_ibjc
      real(dp), dimension(:,:,:,:), allocatable, target  :: L_ibkc
      real(dp), dimension(:,:,:,:), allocatable, target  :: L_jbkc
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: L_ibjc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: L_ibkc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: L_jbkc_p => null()
!
      integer              :: i, j, k, i_rel, j_rel, k_rel
      type(batching_index) :: batch_i, batch_j, batch_k
      integer              :: i_batch, j_batch, k_batch
      integer              :: req_0, req_1, req_2, req_3
!
      logical :: ijk_core
      integer :: i_cvs
!
!
!     Set up arrays for amplitudes
      call mem%alloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call squareup_and_sort_1234_to_1324(wf%t2, t_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Setup and Batching loops
!
      req_0 = 3*wf%n_v**3 + wf%n_v*wf%n_o**3 + wf%n_o*wf%n_v
      req_1 = 3*(wf%n_v)**3
      req_2 = 2*(wf%n_o)*(wf%n_v) + (wf%n_v)**2
      req_3 = 0
!
      batch_i = batching_index(wf%n_o)
      batch_j = batching_index(wf%n_o)
      batch_k = batching_index(wf%n_o)
!
      call mem%batch_setup_ident(batch_i, batch_j, batch_k, &
                           req_0, req_1, req_2, req_3, buffer_size = zero)
!
!     Allocate integral arrays
!
!     Split up so that the integral and amplitude arrays are closer in mem
!
      if (batch_i%num_batches .eq. 1) then ! no batching
!
         call mem%alloc(g_dbic, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%alloc(g_jlic, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         call mem%alloc(L_ibjc, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      else ! batching
!
         call batch_i%determine_limits(1)
!
!        Ordered such that batching indices are at the end
!
         call mem%alloc(g_dbic, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%alloc(g_dbjc, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%alloc(g_dbkc, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
         call mem%alloc(g_jlic, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%alloc(g_klic, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%alloc(g_kljc, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%alloc(g_iljc, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%alloc(g_ilkc, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%alloc(g_jlkc, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
!
         call mem%alloc(L_ibjc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(L_ibkc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(L_jbkc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
!
      endif
!
!     Fock matrix subblock: Resorting for easier contractions later
      call mem%alloc(F_ov_ck, wf%n_v, wf%n_o)
      call sort_12_to_21(wf%fock_ia, F_ov_ck, wf%n_o, wf%n_v)
!
!     Arrays for the triples amplitudes and intermediates
      call mem%alloc(c_abc, wf%n_v, wf%n_v, wf%n_v)
      call mem%alloc(u_abc, wf%n_v, wf%n_v, wf%n_v)
      call mem%alloc(v_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     Remaining integral arrays
!
      if (batch_i%num_batches .eq. 1) then ! no batching
!
         call mem%alloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%alloc(g_ljci, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
         call mem%alloc(Y_bcei, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
!
      else ! batching
!
         call batch_i%determine_limits(1)
!
         call mem%alloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%alloc(g_bdcj, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%alloc(g_bdck, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
         call mem%alloc(g_ljci, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_lkci, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_lkcj, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_licj, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_lick, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_ljck, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
!
         call mem%alloc(Y_bcei, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%alloc(Y_bcej, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%alloc(Y_bcek, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
      end if
!
!     vooo Intermediate
      call mem%alloc(Y_cmjk, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call zero_array(Y_cmjk, wf%n_v*wf%n_o**3)
!
      call wf%g_bdck_t%open_('read')
      call wf%g_ljck_t%open_('read')
      call wf%g_dbkc_t%open_('read')
      call wf%g_jlkc_t%open_('read')
      call wf%L_jbkc_t%open_('read')
!
      wf%Y_bcek = direct_file('Y_bcek', wf%n_v**3)
      call wf%Y_bcek%open_()
!
      do i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(i_batch)
!
         call single_record_reader(batch_i, wf%g_bdck_t, g_bdci, wf%g_dbkc_t, g_dbic)
         g_bdci_p => g_bdci
         g_dbic_p => g_dbic
!
         call zero_array(Y_bcei, batch_i%length*wf%n_v**3)
         Y_bcei_p => Y_bcei
!
         do j_batch = 1, i_batch
!
            call batch_j%determine_limits(j_batch)
!
            call compound_record_reader(batch_j, batch_i, wf%g_ljck_t, g_ljci, &
                                        wf%g_jlkc_t, g_jlic)
            call compound_record_reader(batch_i, batch_j, wf%L_jbkc_t, L_ibjc)
!
            g_ljci_p => g_ljci
            g_jlic_p => g_jlic
            L_ibjc_p => L_ibjc
!
            if (j_batch .ne. i_batch) then
!
               call single_record_reader(batch_j, wf%g_bdck_t, g_bdcj, wf%g_dbkc_t, g_dbjc)
               g_bdcj_p => g_bdcj
               g_dbjc_p => g_dbjc
!
               call single_record_reader(batch_j, wf%Y_bcek, Y_bcej)
               Y_bcej_p => Y_bcej
!
               call compound_record_reader(batch_i, batch_j, wf%g_ljck_t, g_licj, wf%g_jlkc_t, g_iljc)
               g_licj_p => g_licj
               g_iljc_p => g_iljc
!
            else
!
               g_bdcj_p => g_bdci
               g_dbjc_p => g_dbic
!
               Y_bcej_p => Y_bcei
!
               g_licj_p => g_ljci
               g_iljc_p => g_jlic
!
            endif
!
            do k_batch = 1, j_batch
!
               call batch_k%determine_limits(k_batch)
!
               if (k_batch .ne. j_batch) then ! k_batch != j_batch, k_batch != i_batch
!
                  call single_record_reader(batch_k, wf%g_bdck_t, g_bdck, wf%g_dbkc_t, g_dbkc)
                  g_bdck_p => g_bdck
                  g_dbkc_p => g_dbkc
!
                  call single_record_reader(batch_k, wf%Y_bcek, Y_bcek)
                  Y_bcek_p => Y_bcek
! 
                  call compound_record_reader(batch_k, batch_i, wf%g_ljck_t, g_lkci, wf%g_jlkc_t, g_klic)
!
                  g_lkci_p => g_lkci
                  g_klic_p => g_klic
!
                  call compound_record_reader(batch_i, batch_k, wf%g_ljck_t, g_lick, wf%g_jlkc_t, g_ilkc, & 
                                              wf%L_jbkc_t, L_ibkc)
!
                  g_lick_p => g_lick
                  g_ilkc_p => g_ilkc
                  L_ibkc_p => L_ibkc
!
                  call compound_record_reader(batch_k, batch_j, wf%g_ljck_t, g_lkcj, wf%g_jlkc_t, g_kljc)
                  g_lkcj_p => g_lkcj
                  g_kljc_p => g_kljc
!
                  call compound_record_reader(batch_j, batch_k, wf%g_ljck_t, g_ljck, wf%g_jlkc_t, g_jlkc, & 
                                              wf%L_jbkc_t, L_jbkc)
!
                  g_ljck_p => g_ljck
                  g_jlkc_p => g_jlkc
                  L_jbkc_p => L_jbkc
!
               else if (k_batch .eq. i_batch) then ! k_batch == j_batch == i_batch
!
                  g_bdck_p => g_bdci
                  g_dbkc_p => g_dbic
!
                  Y_bcek_p => Y_bcei
!
                  g_lkci_p => g_ljci
                  g_klic_p => g_jlic
!
                  g_lick_p => g_ljci
                  g_ilkc_p => g_jlic
                  L_ibkc_p => L_ibjc
!
                  g_lkcj_p => g_ljci
                  g_kljc_p => g_jlic
!
                  g_ljck_p => g_ljci
                  g_jlkc_p => g_jlic
                  L_jbkc_p => L_ibjc
!
               else ! k_batch == j_batch != i_batch
!
                  g_bdck_p => g_bdcj
                  g_dbkc_p => g_dbjc
!
                  Y_bcek_p => Y_bcej
!
                  g_lkci_p => g_ljci
                  g_klic_p => g_jlic
!
                  g_lick_p => g_licj
                  g_ilkc_p => g_iljc
                  L_ibkc_p => L_ibjc
!
                  call compound_record_reader(batch_k, batch_j, wf%g_ljck_t, g_lkcj, wf%g_jlkc_t,  &
                                             g_kljc, wf%L_jbkc_t, L_jbkc)
                  g_lkcj_p => g_lkcj
                  g_kljc_p => g_kljc
!
                  g_ljck_p => g_lkcj
                  g_jlkc_p => g_kljc
                  L_jbkc_p => L_jbkc
!
               endif
!
               do i = batch_i%first, batch_i%last
!
                  i_rel = i - batch_i%first + 1
!
                  do j = batch_j%first, min(batch_j%last, i)
!
                     j_rel = j - batch_j%first + 1
!
                     do k = batch_k%first, min(batch_k%last, j)
!
                        if (i .eq. j .and. i .eq. k) then
                           cycle
                        end if
!
!                       c3_calc does not zero out the array
                        call zero_array(c_abc, wf%n_v**3)
!
!                       Check if at least one index i,j,k is a core orbital
                        if(cvs) then
!
                           ijk_core = .false.
!
                           do i_cvs = 1, wf%n_core_MOs
!
                              if(     i .eq. wf%core_MOs(i_cvs)   &
                                 .or. j .eq. wf%core_MOs(i_cvs)   &
                                 .or. k .eq. wf%core_MOs(i_cvs))  then
!
                                 ijk_core = .true.
!
                              end if
!
                           end do
!
!                          Cycle if i,j,k are not core orbitals
                           if (.not. ijk_core) cycle
!
                        end if
!
                        k_rel = k - batch_k%first + 1
!
!                       Construct C^abc_ijk for given i, j, k
!                       and construct the intermediates Y_bcek, Y_cmjk 
!                       and calculate contributions to sigma2
!
                        call wf%jacobian_transpose_cc3_c3_calc(i, j ,k, c_ai, c_abij,       &
                                                               c_abc, u_abc,                &
                                                               v_abc, F_ov_ck,              &
                                                               L_ibjc_p(:,:,i_rel,j_rel),   &
                                                               L_ibkc_p(:,:,i_rel,k_rel),   &
                                                               L_jbkc_p(:,:,j_rel,k_rel),   &
                                                               g_dbic_p(:,:,:,i_rel),       &
                                                               g_dbjc_p(:,:,:,j_rel),       &
                                                               g_dbkc_p(:,:,:,k_rel),       &
                                                               g_jlic_p(:,:,j_rel,i_rel),   &
                                                               g_klic_p(:,:,k_rel,i_rel),   &
                                                               g_kljc_p(:,:,k_rel,j_rel),   &
                                                               g_iljc_p(:,:,i_rel,j_rel),   &
                                                               g_ilkc_p(:,:,i_rel,k_rel),   &
                                                               g_jlkc_p(:,:,j_rel,k_rel))
!
                        call wf%divide_by_orbital_differences(i, j, k, c_abc, omega)
!
                        call wf%jacobian_transpose_cc3_a_n7(i, j, k, c_abc, u_abc, sigma_abij, &
                                                            g_bdci_p(:,:,:,i_rel),             &
                                                            g_bdcj_p(:,:,:,j_rel),             &
                                                            g_bdck_p(:,:,:,k_rel),             &
                                                            g_ljci_p(:,:,j_rel,i_rel),         &
                                                            g_lkci_p(:,:,k_rel,i_rel),         &
                                                            g_lkcj_p(:,:,k_rel,j_rel),         &
                                                            g_licj_p(:,:,i_rel,j_rel),         &
                                                            g_lick_p(:,:,i_rel,k_rel),         &
                                                            g_ljck_p(:,:,j_rel,k_rel))
!
                        call wf%construct_y_intermediates(i, j, k, c_abc, u_abc, t_abij,  &
                                                          Y_cmjk,                         &
                                                          Y_bcei_p(:,:,:,i_rel),          &
                                                          Y_bcej_p(:,:,:,j_rel),          &
                                                          Y_bcek_p(:,:,:,k_rel))
!
                     enddo ! loop over k
                  enddo ! loop over j
               enddo ! loop over i
!
!              write the intermediate Y_bcek to file. 
!              Will be read in after the loops for the contractions to sigma_ai
!
               if (k_batch .ne. j_batch) then !k_batch != j_batch, k_batch != i_batch
                  call single_record_writer(batch_k, wf%Y_bcek, Y_bcek)
               endif
!
            enddo ! batch_k
!
            if (j_batch .ne. i_batch) then
               call single_record_writer(batch_j, wf%Y_bcek, Y_bcej)
            endif
!
         enddo ! batch_j
!
         call single_record_writer(batch_i, wf%Y_bcek, Y_bcei)
!
      enddo ! batch_i
!
      call wf%g_bdck_t%close_()
      call wf%g_ljck_t%close_()
      call wf%g_dbkc_t%close_()
      call wf%g_jlkc_t%close_()
      call wf%L_jbkc_t%close_()
!
!     Deallocate the integral arrays
!
      if (batch_i%num_batches .eq. 1) then
!
         call mem%dealloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%dealloc(g_dbic, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%dealloc(g_ljci, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
         call mem%dealloc(g_jlic, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         call mem%dealloc(L_ibjc, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
         call mem%dealloc(Y_bcei, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
!
      else
         call batch_i%determine_limits(1)
!
         call mem%dealloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%dealloc(g_bdcj, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%dealloc(g_bdck, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
         call mem%dealloc(g_dbic, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%dealloc(g_dbjc, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%dealloc(g_dbkc, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
         call mem%dealloc(g_ljci, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_lkci, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_lkcj, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_licj, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_lick, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_ljck, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
!
         call mem%dealloc(g_jlic, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%dealloc(g_klic, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%dealloc(g_kljc, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%dealloc(g_iljc, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%dealloc(g_ilkc, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%dealloc(g_jlkc, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
!
         call mem%dealloc(L_ibjc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(L_ibkc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(L_jbkc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
!
         call mem%dealloc(Y_bcei, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%dealloc(Y_bcej, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%dealloc(Y_bcek, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
      endif
!
!     Deallocate amplitudes arrays and Fock matrix
!
      call mem%dealloc(c_abc, wf%n_v, wf%n_v, wf%n_v)
      call mem%dealloc(u_abc, wf%n_v, wf%n_v, wf%n_v)
      call mem%dealloc(v_abc, wf%n_v, wf%n_v, wf%n_v)
!
      call mem%dealloc(F_ov_ck, wf%n_v, wf%n_o)
!
      call mem%dealloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
!     Contribution of the Y_cmjk to sigma1
!
      call wf%jacobian_transpose_cc3_c3_a1_y_o(sigma_ai, Y_cmjk)
!
      call mem%dealloc(Y_cmjk, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
!     Contribution of the Y_bcek to sigma1
!
      call wf%jacobian_transpose_cc3_c3_b1_y_v(sigma_ai)
!
      call wf%Y_bcek%close_()
!
   end subroutine jacobian_transpose_cc3_c3_a_cc3
!
!
   module subroutine jacobian_transpose_cc3_c3_calc_cc3(wf, i, j ,k, c_ai, c_abij, &
                                                        c_abc, u_abc,              &
                                                        v_abc, F_ov_ck,               &
                                                        L_ibjc, L_ibkc, L_jbkc,    &
                                                        g_dbic, g_dbjc, g_dbkc,    &
                                                        g_jlic, g_klic, g_kljc,    &
                                                        g_iljc, g_ilkc, g_jlkc)
!!
!!    Construct c3 amplitudes for fixed i,j,k
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!!    C^abc_ijk 
!!    = (omega - epsilon^abc_ijk)^-1 P^abc_ijk (C__ai*L_jbkc - C__ak*L_jbic 
!!                                  + C_abij*F_kc - C__abik*F_jc
!!                                  + sum_l (C__ablk g_iljc - C__abil L_jlkc) 
!!                                  - sum_d (C__adjk g_ibdc - C__adij L_dbkc))
!!
!!    Contributions from outer products:
!!    P^abc_ijk (C_ai*L_jbkc - C_ak*L_jbic + Cabij*F_kc - C_abik*F_jc)
!!
!!    Contributions from matrix multiplication:
!!      sum_l P^abc_ijk (C_ablk g_iljc + C_abil g_jckl - 2 C_abil g_jlkc) 
!!    - sum_d P^abc_ijk (C_adjk g_ibdc + C_adij g_dckb - 2 C_adij g_dbkc)
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: c_abij
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: c_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: u_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: v_abc
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: F_ov_ck
!
!     L_ibjc ordered bc,ij
      real(dp), dimension(wf%n_v, wf%n_v), intent(in) :: L_ibjc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in) :: L_ibkc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in) :: L_jbkc
!
!     g_dbkc ordered bcd,k
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in) :: g_dbic
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in) :: g_dbjc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in) :: g_dbkc
!
!     g_jlkc ordered cl,jk
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_jlic
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_klic
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_kljc
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_iljc
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_ilkc
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_jlkc
!
!
!     :: Contribution 1 ::
!
!     c_abc <- u_abc = -sum_l (2* c_ablj g_ilkc - c_aclj g_ilkb - c_cblj g_ilka)
!
      call dgemm('N', 'T',          &
                  wf%n_v**2,        &
                  wf%n_v,           &
                  wf%n_o,           &
                  -one,             &
                  c_abij(:,:,:,j),  & ! c_ab_l,j c_ac_l,j c_cb_l,j
                  wf%n_v**2,        &
                  g_ilkc,           & ! g_c_l,ik g_b_l,ik g_a_l,ik
                  wf%n_v,           &
                  zero,             &
                  u_abc,            &
                  wf%n_v**2)
!
!     c_abc <- u_abc = sum_d (2* c_adij g_dbkc - c_adij g_dckb - c_cdij g_dbka)
!
      call dgemm('N', 'T',          &
                  wf%n_v,           &
                  wf%n_v**2,        &
                  wf%n_v,           &
                  one,              &
                  c_abij(:,:,i,j),  & ! c_a_d,ij c_c_d,ij
                  wf%n_v,           &
                  g_dbkc,           & ! g_bc_d,k g_cb_d,k g_ba_d,k
                  wf%n_v**2,        &
                  one,              &
                  u_abc,            &
                  wf%n_v)
!
!     c_abc <- u_abc = 2*L_iajb*c_ck - L_iajc*c_bk - L_icjb*c_ak
!
      call dger(wf%n_v**2, &
               wf%n_v,     &
               one,        &
               L_ibjc,     & ! L_ab,ij L_ac,ij L_cb,ij
               1,          &
               c_ai(:,k),  & ! c_c,k c_b,k c_a,k
               1,          &
               u_abc,      &
               wf%n_v**2)
!
!     c_abc <- u_abc = 2*c^ab_ij*F_kc - c^ac_ij*F_kb - c^cb_ij*F_ka
!
      call dger(wf%n_v**2,       &
               wf%n_v,           &
               one,              &
               c_abij(:,:,i,j),  & ! c_ab,ij c_ac,ij c_cb,ij
               1,                &
               F_ov_ck(:,k),     & ! F_c,k F_b,k F_a,k
               1,                &
               u_abc,            &
               wf%n_v**2)
!
!     c_abc <- u_abc <- v_abc = - sum_l (2*c_bali g_jlkc - c_bcli g_jlka - c_cali g_jlkb)
!
      call dgemm('N', 'T',          &
                  wf%n_v**2,        &
                  wf%n_v,           &
                  wf%n_o,           &
                  -one,             &
                  c_abij(:,:,:,i),  & ! c_ba_l,i c_bc_l,i c_ca_l,i
                  wf%n_v**2,        &
                  g_jlkc,           & ! g_c_l,jk g_b_l,jk g_a_l,jk
                  wf%n_v,           &
                  zero,             &
                  v_abc,            &
                  wf%n_v**2)
!
!     c_abc <- u_abc <- v_abc = sum_d (2*c_bdji g_dakc - c_bdji g_dcka - c_cdji g_dakb)
!
      call dgemm('N', 'T',          &
                  wf%n_v,           &
                  wf%n_v**2,        &
                  wf%n_v,           &
                  one,              &
                  c_abij(:,:,j,i),  & ! c_b_d,ji c_c_d,ji
                  wf%n_v,           &
                  g_dbkc,           & ! g_ac_d,k g_ca_d,k g_ab_d,k
                  wf%n_v**2,        &
                  one,              &
                  v_abc,            &
                  wf%n_v)
!
      call sort_123_to_213_and_add(v_abc, u_abc, wf%n_v, wf%n_v, wf%n_v)
      call add_two_123_min_132_min_321(u_abc, c_abc, wf%n_v)
!
!
!     :: Contribution 2 ::
!
!     c_abc <- u_abc = - sum_l (2* c_aclk g_iljb - c_ablk g_iljc - c_bclk g_ilja)
!
      call dgemm('N', 'T',          &
                  wf%n_v**2,        &
                  wf%n_v,           &
                  wf%n_o,           &
                  -one,             &
                  c_abij(:,:,:,k),  & ! c_ac_l,k c_ab_l,k c_bc_l,k
                  wf%n_v**2,        &
                  g_iljc,           & ! g_b_l,ij g_c_l,ij g_a_l,ij
                  wf%n_v,           &
                  zero,             &
                  u_abc,            &
                  wf%n_v**2)
!
!     c_abc <- u_abc = sum_d (2* c_adik g_dcjb - c_adik g_dbjc - c_bdik g_dcja)
!
      call dgemm('N', 'T',          &
                  wf%n_v,           &
                  wf%n_v**2,        &
                  wf%n_v,           &
                  one,              &
                  c_abij(:,:,i,k),  & ! c_a_d,ik c_b_d,ik
                  wf%n_v,           &
                  g_dbjc,           & ! g_cb_d,j g_bc_d,j g_ca_d,j
                  wf%n_v**2,        &
                  one,              &
                  u_abc,            &
                  wf%n_v)
!
!     c_abc <- u_abc = 2*L_iakc*c_bj - L_iakb*c_cj - L_ibkc*c_aj
!
      call dger(wf%n_v**2, &
               wf%n_v,     &
               one,        &
               L_ibkc,     & ! L_ac,ik L_ab,ik L_bc,ik
               1,          &
               c_ai(:,j),  & ! c_b,j c_c,j c_a,j
               1,          &
               u_abc,      &
               wf%n_v**2)
!
!     c_abc <- u_abc = 2*c^ac_ik*F_jb - c^ab_ik*F_jc - c^bc_ik*F_ja
!
      call dger(wf%n_v**2,       &
               wf%n_v,           &
               one,              &
               c_abij(:,:,i,k),  & ! c_ac,ik c_ab,ik c_bc,ik
               1,                &
               F_ov_ck(:,j),     & ! F_b,j F_c,j F_a,j
               1,                &
               u_abc,            &
               wf%n_v**2)
!
!     c_cab <- u_abc <- v_abc = - sum_l (2*c_cali g_kljb - c_cbli g_klja - c_bali g_kljc)
!
      call dgemm('N', 'T',          &
                  wf%n_v**2,        &
                  wf%n_v,           &
                  wf%n_o,           &
                  -one,             &
                  c_abij(:,:,:,i),  & ! c_ca_l,i c_cb_l,i c_ba_l,i
                  wf%n_v**2,        &
                  g_kljc,           & ! g_b_l,kj g_a_l,kj g_c_l,kj
                  wf%n_v,           &
                  zero,             &
                  v_abc,            &
                  wf%n_v**2)
!
!     c_abc <- u_abc <- v_abc = sum_d (2*c_cdki g_dajb - c_cdki g_dbja - c_bdki g_dajc)
!
      call dgemm('N', 'T',          &
                  wf%n_v,           &
                  wf%n_v**2,        &
                  wf%n_v,           &
                  one,              &
                  c_abij(:,:,k,i),  & ! c_c_d,ki c_b_d,ki
                  wf%n_v,           &
                  g_dbjc,           & ! g_ab_d,j g_ba_d,j g_ac_d,j
                  wf%n_v**2,        &
                  one,              &
                  v_abc,            &
                  wf%n_v)
!
      call sort_123_to_213_and_add(v_abc, u_abc, wf%n_v, wf%n_v, wf%n_v)
      call add_two_132_min_231_min_123(u_abc, c_abc, wf%n_v)
!
!
!     :: Contribution 3 ::
!
!     c_abc <- u_abc = - sum_l (2*c_bclk g_jlia - c_balk g_jlic - c_aclk g_jlib)
!
      call dgemm('N', 'T',          &
                  wf%n_v**2,        &
                  wf%n_v,           &
                  wf%n_o,           &
                  -one,             &
                  c_abij(:,:,:,k),  & ! c_bc_l,k c_ba_l,k c_ac_l,k
                  wf%n_v**2,        &
                  g_jlic,           & ! g_a_l,ji g_c_l,ji g_b_l,ji
                  wf%n_v,           &
                  zero,             &
                  u_abc,            &
                  wf%n_v**2)
!
!     c_abc <- u_abc = sum_d (2* c_bdjk g_dcia - c_bdjk g_daic - c_adjk g_dcib)
!
      call dgemm('N', 'T',          &
                  wf%n_v,           &
                  wf%n_v**2,        &
                  wf%n_v,           &
                  one,              &
                  c_abij(:,:,j,k),  & ! c_b_d,jk c_a_d,jk
                  wf%n_v,           &
                  g_dbic,           & ! g_ca_d,i g_ac_d,i g_cb_d,i
                  wf%n_v**2,        &
                  one,              &
                  u_abc,            &
                  wf%n_v)
!
!     c_abc <- u_abc = 2*L_jbkc*c_ai - L_jbka*c_ci - L_jakc*c_bi
!
      call dger(wf%n_v**2, &
               wf%n_v,     &
               one,        &
               L_jbkc,     & ! L_bc,jk L_ba,jk L_ac,jk
               1,          &
               c_ai(:,i),  & ! c_a,i c_c,i c_b,i
               1,          &
               u_abc,      &
               wf%n_v**2)
!
!     c_abc <- u_abc = 2*c^bc_jk*F_ia - c^ba_jk*F_ic - c^ac_jk*F_ib
!
      call dger(wf%n_v**2,       &
               wf%n_v,           &
               one,              &
               c_abij(:,:,j,k),  & ! c_bc,jk c_ba,jk c_ac,jk
               1,                &
               F_ov_ck(:,i),     & ! F_a,i F_c,i F_b,i
               1,                &
               u_abc,            &
               wf%n_v**2)
!
!     c_abc <- u_abc <- v_abc = - sum_l (2*c_cblj g_klia - c_calj g_klib - c_ablj g_klic)
!
      call dgemm('N', 'T',          &
                  wf%n_v**2,        &
                  wf%n_v,           &
                  wf%n_o,           &
                  -one,             &
                  c_abij(:,:,:,j),  & ! c_cb_l,j c_ca_l,j c_ab_l,j
                  wf%n_v**2,        &
                  g_klic,           & ! g_a_l,ki g_b_l,ki g_c_l,ki
                  wf%n_v,           &
                  zero,             &
                  v_abc,            &
                  wf%n_v**2)
! 
!     c_abc <- u_abc <- v_abc = sum_d (2*c_cdkj g_dbia - c_cdkj g_daib - c_adkj g_dbic)
!
      call dgemm('N', 'T',          &
                  wf%n_v,           &
                  wf%n_v**2,        &
                  wf%n_v,           &
                  one,              &
                  c_abij(:,:,k,j),  & ! c_c_d,kj c_a_d,kj
                  wf%n_v,           &
                  g_dbic,           & ! g_ba_d,i g_ab_d,i g_bc_d,i
                  wf%n_v**2,        &
                  one,              &
                  v_abc,            &
                  wf%n_v)
!
      call sort_123_to_213_and_add(v_abc, u_abc, wf%n_v, wf%n_v, wf%n_v)
      call add_two_231_min_213_min_132(u_abc, c_abc, wf%n_v)
!
!
   end subroutine jacobian_transpose_cc3_c3_calc_cc3
!
!
   module subroutine jacobian_transpose_cc3_a_n7_cc3(wf, i, j, k, c_abc, u_abc, sigma_abij,  &
                                                     g_bdci, g_bdcj, g_bdck, g_ljci, g_lkci, &
                                                     g_lkcj, g_licj, g_lick, g_ljck)
!!
!!    Jacobian transpose A2 term
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!!    sigma_adij =   sum_ckd c^abc_ijk g_bdck
!!    sigma_abil = - sum_cki c^bac_ijk g_lick
!!
!!    All permutations for i,j,k have to be considered 
!!    due to the restrictions in the i,j,k loops
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)           :: c_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: u_abc
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(out)  :: sigma_abij
!
!     g_bdck ordered as dbc,k
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)           :: g_bdci
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)           :: g_bdcj
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)           :: g_bdck
!
!     g_ljck ordered as lc,jk
      real(dp), dimension(wf%n_o, wf%n_v), intent(in)                   :: g_ljci
      real(dp), dimension(wf%n_o, wf%n_v), intent(in)                   :: g_lkci
      real(dp), dimension(wf%n_o, wf%n_v), intent(in)                   :: g_lkcj
      real(dp), dimension(wf%n_o, wf%n_v), intent(in)                   :: g_licj
      real(dp), dimension(wf%n_o, wf%n_v), intent(in)                   :: g_lick
      real(dp), dimension(wf%n_o, wf%n_v), intent(in)                   :: g_ljck
!
!     sigma_adij += sum_bc,k c^abc_ijk g_bdck
!
      call dgemm('N', 'T',             &
                  wf%n_v,              &
                  wf%n_v,              &
                  wf%n_v**2,           &
                  one,                 &
                  c_abc,               & ! c_a_bc
                  wf%n_v,              &
                  g_bdck,              & ! g_d_bc,k
                  wf%n_v,              &
                  one,                 &
                  sigma_abij(:,:,i,j), & ! sigma_a_d,ij
                  wf%n_v)
!
!     sigma_ablj += - sum_c c^abc_ijk g_lick
!
      call dgemm('N', 'T',             &
                  wf%n_v**2,           &
                  wf%n_o,              &
                  wf%n_v,              &
                  -one,                &
                  c_abc,               & ! c_ab_c
                  wf%n_v**2,           &
                  g_lick,              & ! g_l_c,ik
                  wf%n_o,              &
                  one,                 &
                  sigma_abij(:,:,:,j), & ! sigma_ab_l,j
                  wf%n_v**2)
!
!
!     sigma_adki += sum_bc,j c^bca_ijk g_bdcj
!
      call dgemm('T', 'T',             &
                  wf%n_v,              &
                  wf%n_v,              &
                  wf%n_v**2,           &
                  one,                 &
                  c_abc,               & ! c_bc_a
                  wf%n_v**2,           &
                  g_bdcj,              & ! g_d_bc,j
                  wf%n_v,              &
                  one,                 &
                  sigma_abij(:,:,k,i), & ! sigma_a_d,ki
                  wf%n_v)
!
!     sigma_ablk += -sum_c c^cab_ijk g_ljci
!
      call dgemm('T', 'T',             &
                  wf%n_v**2,           &
                  wf%n_o,              &
                  wf%n_v,              &
                  -one,                &
                  c_abc,               & ! c_c_ab
                  wf%n_v,              &
                  g_ljci,              & ! g_l_c,ji
                  wf%n_o,              &
                  one,                 &
                  sigma_abij(:,:,:,k), & ! sigma_ab_l,k
                  wf%n_v**2)
!
!
      call sort_123_to_231(c_abc, u_abc, wf%n_v, wf%n_v, wf%n_v)
!     123_to_231(cab) -> abc
!     123_to_231(bca) -> cab
!
!     sigma_adjk += sum_bc,i c^cab_ijk g_bdci
!
      call dgemm('N', 'T',             &
                  wf%n_v,              &
                  wf%n_v,              &
                  wf%n_v**2,           &
                  one,                 &
                  u_abc,               & ! c_a_bc
                  wf%n_v,              &
                  g_bdci,              & ! g_d_bc,i
                  wf%n_v,              &
                  one,                 &
                  sigma_abij(:,:,j,k), & ! sigma_a_d,jk
                  wf%n_v)
!
!     sigma_abli += -sum_c c^bca_ijk g_lkcj
!
      call dgemm('T', 'T',             &
                  wf%n_v**2,           &
                  wf%n_o,              &
                  wf%n_v,              &
                  -one,                &
                  u_abc,               & ! c_c_ab
                  wf%n_v,              &
                  g_lkcj,              & ! g_l_c,kj
                  wf%n_o,              &
                  one,                 &
                  sigma_abij(:,:,:,i), & ! sigma_ab_l,i
                  wf%n_v**2)
!
!
      if (k .ne. j .and. j .ne. i) then
!
         call sort_123_to_213(c_abc, u_abc, wf%n_v, wf%n_v, wf%n_v)
!        123_to_213(bac) -> abc
!        123_to_213(acb) -> cab
!        123_to_213(cba) -> bca
!
!        sigma_adji += sum_bc,k c^bac_ijk g_bdck
!
         call dgemm('N', 'T',             &
                     wf%n_v,              &
                     wf%n_v,              &
                     wf%n_v**2,           &
                     one,                 &
                     u_abc,               & ! c_a_bc
                     wf%n_v,              &
                     g_bdck,              & ! g_d_bc,k
                     wf%n_v,              &
                     one,                 &
                     sigma_abij(:,:,j,i), & ! sigma_a_d,ji
                     wf%n_v)
!
!        sigma_abli += - sum_c c^bac_ijk g_ljck
!
         call dgemm('N', 'T',             &
                     wf%n_v**2,           &
                     wf%n_o,              &
                     wf%n_v,              &
                     -one,                &
                     u_abc,               & ! c_ab_c
                     wf%n_v**2,           &
                     g_ljck,              & ! g_l_c,jk
                     wf%n_o,              &
                     one,                 &
                     sigma_abij(:,:,:,i), & ! sigma_ab_l,i
                     wf%n_v**2)
!
!
!        sigma_adkj += sum_bc c^cba_ijk g_bdci
!
         call dgemm('T', 'T',             &
                     wf%n_v,              &
                     wf%n_v,              &
                     wf%n_v**2,           &
                     one,                 &
                     u_abc,               & ! c_bc_a
                     wf%n_v**2,           &
                     g_bdci,              & ! g_d_bc,i
                     wf%n_v,              &
                     one,                 &
                     sigma_abij(:,:,k,j), & ! sigma_a_d,kj
                     wf%n_v)
!
!        sigma_ablk += - sum_c c^acb_ijk g_licj
!
         call dgemm('T', 'T',             &
                     wf%n_v**2,           &
                     wf%n_o,              &
                     wf%n_v,              &
                     -one,                &
                     u_abc,               & ! c_c_ab
                     wf%n_v,              &
                     g_licj,              & ! g_l_c,ij
                     wf%n_o,              &
                     one,                 &
                     sigma_abij(:,:,:,k), & ! sigma_ab_l,k
                     wf%n_v**2)
!
!
         call sort_123_to_321(c_abc, u_abc, wf%n_v, wf%n_v, wf%n_v)
!        123_to_321(acb) -> bca
!        123_to_321(cba) -> abc
!
!        sigma_adik += sum_bc c^acb_ijk g_bdcj
!
         call dgemm('T', 'T',             &
                     wf%n_v,              &
                     wf%n_v,              &
                     wf%n_v**2,           &
                     one,                 &
                     u_abc,               & ! c_bc_a
                     wf%n_v**2,           &
                     g_bdcj,              & ! g_d_bc,j
                     wf%n_v,              &
                     one,                 &
                     sigma_abij(:,:,i,k), & ! sigma_ad,ik
                     wf%n_v)
!
!        sigma_ablj += - sum_c c^cba_ijk g_lkci
!
         call dgemm('N', 'T',             &
                     wf%n_v**2,           &
                     wf%n_o,              &
                     wf%n_v,              &
                     -one,                &
                     u_abc,               & ! c_ab_c
                     wf%n_v**2,           &
                     g_lkci,              & ! g_l_c,ki
                     wf%n_o,              &
                     one,                 &
                     sigma_abij(:,:,:,j), & ! sigma_ab_l,j
                     wf%n_v**2)
!
      end if
!
!
   end subroutine jacobian_transpose_cc3_a_n7_cc3
!
!
   module subroutine construct_y_intermediates_cc3(wf, i, j, k, c_abc, u_abc, t_abij, y_cmjk,   &
                                                   Y_bcei, Y_bcej, Y_bcek)
!!
!!    Construct Y intermediates
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!!    Y_bcei = sum_aij c^abc_ijk * t^ae_ij
!!    Y_cmjk = sum_abj c^bac_ijk * t^ba_mj
!!
!!    All permutations for i,j,k have to be considered due to the restrictions in the i,j,k loops
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)           :: c_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: u_abc
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in)   :: t_abij
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_o, wf%n_o), intent(out)  :: Y_cmjk
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: Y_bcei
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: Y_bcej
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: Y_bcek
!
!     Y_cmik = sum_ab,j c^abc_ijk t^ab_mj
!
      call dgemm('T','N',           &
                  wf%n_v,           &
                  wf%n_o,           &
                  wf%n_v**2,        &
                  one,              &
                  c_abc,            & ! c_ab_c
                  wf%n_v**2,        &
                  t_abij(:,:,:,j),  & ! t_ab_m,j
                  wf%n_v**2,        &
                  one,              &
                  Y_cmjk(:,:,i,k),  & ! Y_c_m,ik
                  wf%n_v)
!
!     Y_bcek = sum_a,ij c^abc_ijk t^ae_ij
!
      call dgemm('T','N',           &
                  wf%n_v**2,        &
                  wf%n_v,           &
                  wf%n_v,           &
                  one,              &
                  c_abc,            & ! c_a_bc
                  wf%n_v,           &
                  t_abij(:,:,i,j),  & ! t_a_e,ij
                  wf%n_v,           &
                  one,              &
                  Y_bcek,           & ! Y_bc_e,k
                  wf%n_v**2)
!
!
!     Y_cmji = sum_ab,k c^cab_ijk t^ab_mk
!
      call dgemm('N','N',           &
                  wf%n_v,           &
                  wf%n_o,           &
                  wf%n_v**2,        &
                  one,              &
                  c_abc,            & ! c_c_ab
                  wf%n_v,           &
                  t_abij(:,:,:,k),  & ! t_ab_m,k
                  wf%n_v**2,        &
                  one,              &
                  Y_cmjk(:,:,j,i),  & ! Y_c_m,ji
                  wf%n_v)
!
!     Y_bcej = sum_a,ki c^bca_ijk t^ae_ki
!
      call dgemm('N','N',           &
                  wf%n_v**2,        &
                  wf%n_v,           &
                  wf%n_v,           &
                  one,              &
                  c_abc,            & ! c_bc_a
                  wf%n_v**2,        &
                  t_abij(:,:,k,i),  & ! t_a_e_ki
                  wf%n_v,           &
                  one,              &
                  Y_bcej,           & ! Y_bc_e,j
                  wf%n_v**2)
!
!
      call sort_123_to_231(c_abc, u_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     123->231(bca) = cab
!     123->231(cab) = abc
!
!     Y_cmkj = sum_ab,i c^bca_ijk t^ab_mi
!
      call dgemm('N','N',           &
                  wf%n_v,           &
                  wf%n_o,           &
                  wf%n_v**2,        &
                  one,              &
                  u_abc,            & ! c_c_ab
                  wf%n_v,           &
                  t_abij(:,:,:,i),  & ! t_ab_m,i
                  wf%n_v**2,        &
                  one,              &
                  Y_cmjk(:,:,k,j),  & ! Y_c_m,kj
                  wf%n_v)
!
!
!     Y_bcei = sum_a,jk c^cab_ijk t^ae_jk
!
      call dgemm('T','N',           &
                  wf%n_v**2,        &
                  wf%n_v,           &
                  wf%n_v,           &
                  one,              &
                  u_abc,            & ! c_a_bc
                  wf%n_v,           &
                  t_abij(:,:,j,k),  & ! t_a_e,jk
                  wf%n_v,           &
                  one,              &
                  Y_bcei,           & ! Y_bc_e,i
                  wf%n_v**2)
!
      if (k .ne. j .and. j .ne. i) then
!
         call sort_123_to_213(c_abc, u_abc, wf%n_v, wf%n_v, wf%n_v)
!
!        123->213(bac) = abc
!        123->213(acb) = cab
!        123->213(cba) = bca
!
!        Y_cmjk = sum_ab,i c^bac_ijk t^ab_mi
!
         call dgemm('T','N',           &
                     wf%n_v,           &
                     wf%n_o,           &
                     wf%n_v**2,        &
                     one,              &
                     u_abc,            & ! c_ab_c
                     wf%n_v**2,        &
                     t_abij(:,:,:,i),  & ! t_ab_m,i
                     wf%n_v**2,        &
                     one,              &
                     Y_cmjk(:,:,j,k),  & ! Y_c_m,jk
                     wf%n_v)
!
!        Y_bcek = sum_a,ji c^bac_ijk t^ae_ji
!
         call dgemm('T','N',           &
                     wf%n_v**2,        &
                     wf%n_v,           &
                     wf%n_v,           &
                     one,              &
                     u_abc,            & ! c_a_bc
                     wf%n_v,           &
                     t_abij(:,:,j,i),  & ! t_a_e,ji
                     wf%n_v,           &
                     one,              &
                     Y_bcek,           & ! Y_bc_e,k
                     wf%n_v**2)
!
!
!        Y_cmij = sum_ab,k c^acb_ijk t^ab_mk
!
         call dgemm('N','N',           &
                     wf%n_v,           &
                     wf%n_o,           &
                     wf%n_v**2,        &
                     one,              &
                     u_abc,            & ! c_c_ab
                     wf%n_v,           &
                     t_abij(:,:,:,k),  & ! t_ab_m,k
                     wf%n_v**2,        &
                     one,              &
                     Y_cmjk(:,:,i,j),  & ! Y_c_m,ij
                     wf%n_v)
!
!        Y_bcei = sum_a,kj c^cba t^ae_kj
!
         call dgemm('N','N',           &
                     wf%n_v**2,        &
                     wf%n_v,           &
                     wf%n_v,           &
                     one,              &
                     u_abc,            & ! c_bc_a
                     wf%n_v**2,        &
                     t_abij(:,:,k,j),  & ! t_a_e,kj
                     wf%n_v,           &
                     one,              &
                     Y_bcei,           & ! Y_bc_e,i
                     wf%n_v**2)
!
!
         call sort_123_to_132(c_abc, u_abc, wf%n_v, wf%n_v, wf%n_v)
!
!        123->132(cba) = cab
!        123->132(acb) = abc
!
!        Y_cmki = sum_ab,j c^cba_ijk t^ab_mj
!
         call dgemm('N','N',           &
                     wf%n_v,           &
                     wf%n_o,           &
                     wf%n_v**2,        &
                     one,              &
                     u_abc,            & ! c_c_ab
                     wf%n_v,           &
                     t_abij(:,:,:,j),  & ! t_ab_m,j
                     wf%n_v**2,        &
                     one,              &
                     Y_cmjk(:,:,k,i),  & ! Y_c_m,ki
                     wf%n_v)
!
!
!        Y_bcej = sum_a,ik c^acb_ijk t^ae_ik
!
         call dgemm('T','N',           &
                     wf%n_v**2,        &
                     wf%n_v,           &
                     wf%n_v,           &
                     one,              &
                     u_abc,            & ! c_a_bc
                     wf%n_v,           &
                     t_abij(:,:,i,k),  & ! t_a_e,ik
                     wf%n_v,           &
                     one,              &
                     Y_bcej,           & ! Y_bc_e,j
                     wf%n_v**2)
!
      end if
!
   end subroutine construct_y_intermediates_cc3
!
!
   module subroutine jacobian_transpose_cc3_c3_a1_y_o_cc3(wf, sigma_ai, Y_cmjk)
!!
!!    Jacobian tranpose contribution of Y_vooo
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!!    sigma_1 += sum_mjk Y_cmjk * g_mjlk
!!    sigma_1 += sum_cmj g_mjcd * Y_cmjk
!!    sigma_1 += sum_cmk g_leck * Y_cmjk
!!   
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_o, wf%n_o), intent(in) :: Y_cmjk
!
      real(dp), dimension(:,:,:,:), allocatable :: g_mjlk ! ordered g_mjkl
      real(dp), dimension(:,:,:,:), allocatable :: g_cdmj ! ordered g_cmjd
      real(dp), dimension(:,:,:,:), allocatable :: g_cjmd ! ordered g_cmjd
!
!     arrays for resorting
      real(dp), dimension(:,:,:,:), allocatable :: Y_cmkj
!
      type(batching_index) :: batch_l, batch_d
      integer :: l_batch, d_batch
      integer :: req_0, req_1
!
!     :: Term 1: sigma_cl += sum_mjk Y_cmjk * g_mjlk ::
!
      call wf%g_mjlk_t%open_('read')
!
      batch_l = batching_index(wf%n_o)
!
      req_0 = 0
      req_1 = wf%n_o**3
!
      call mem%batch_setup(batch_l, req_0, req_1)
!
      call batch_l%determine_limits(1)
      call mem%alloc(g_mjlk, wf%n_o, wf%n_o, wf%n_o, batch_l%length) ! g_mjk#l
!
      do l_batch = 1, batch_l%num_batches
!
         call batch_l%determine_limits(l_batch)
!
         call single_record_reader(batch_l, wf%g_mjlk_t, g_mjlk)
!
         call dgemm('N','N',                    &
                     wf%n_v,                    &
                     batch_l%length,            &
                     wf%n_o**3,                 &
                     one,                       &
                     Y_cmjk,                    & ! Y_c_mjk
                     wf%n_v,                    &
                     g_mjlk,                    & ! g_mjk_l
                     wf%n_o**3,                 &
                     one,                       &
                     sigma_ai(1,batch_l%first), &
                     wf%n_v)
!
      enddo ! l_batch
!
      call batch_l%determine_limits(1)
      call mem%dealloc(g_mjlk, wf%n_o, wf%n_o, wf%n_o, batch_l%length)
!
      call wf%g_mjlk_t%close_()
!
!
!     :: Term 2: sigma_dk += - sum_cmj g_cdmj * Y_cmjk ::
!
      call wf%g_cdlk_t%open_('read')
!
      batch_d = batching_index(wf%n_v)
!
      req_0 = 0
      req_1 = 2*wf%n_v*wf%n_o**2
!
      call mem%batch_setup(batch_d, req_0, req_1)
!
      call batch_d%determine_limits(1)
      call mem%alloc(g_cdmj, wf%n_v, wf%n_o, wf%n_o, batch_d%length) ! cmj#d
!
      do d_batch = 1, batch_d%num_batches
!
         call batch_d%determine_limits(d_batch)
!
!        read g_cdmj stored as cmj#d
         call compound_record_reader(wf%n_o, batch_d, wf%g_cdlk_t, g_cdmj)
!
         call dgemm('T','N',                       &
                     batch_d%length,               &
                     wf%n_o,                       &
                     wf%n_v*wf%n_o**2,             &
                     -one,                         &
                     g_cdmj,                       & ! g_cmj_#d
                     wf%n_v*wf%n_o**2,             &
                     Y_cmjk,                       & ! Y_cmj_k
                     wf%n_v*wf%n_o**2,             &
                     one,                          &
                     sigma_ai(batch_d%first,1),    &
                     wf%n_v)
!
      enddo
!
      call batch_d%determine_limits(1)
      call mem%dealloc(g_cdmj, wf%n_v, wf%n_o, wf%n_o, batch_d%length) ! cmj#d
!
      call wf%g_cdlk_t%close_()
!
!
!     :: Term 3: sigma_dk += - sum_cmj g_cjmd * Y_cmkj ::
!
      call mem%alloc(Y_cmkj, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call sort_1234_to_1243(Y_cmjk, Y_cmkj, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call wf%g_ckld_t%open_('read')
!
      batch_d = batching_index(wf%n_v)
!
      req_0 = 0
      req_1 = wf%n_v*wf%n_o**2
!
      call mem%batch_setup(batch_d, req_0, req_1)
!
      call batch_d%determine_limits(1)
!
      call mem%alloc(g_cjmd, wf%n_v, wf%n_o, wf%n_o, batch_d%length) ! cmj#d
!
      do d_batch = 1, batch_d%num_batches
!
         call batch_d%determine_limits(d_batch)
!
         call compound_record_reader(wf%n_o, batch_d, wf%g_ckld_t, g_cjmd) ! cmj#d
!
         call dgemm('T','N',                       &
                     batch_d%length,               &
                     wf%n_o,                       &
                     wf%n_v*wf%n_o**2,             &
                     -one,                         &
                     g_cjmd,                       & ! g_cmj_#d
                     wf%n_v*wf%n_o**2,             &
                     Y_cmkj,                       & ! Y_cmj_k
                     wf%n_v*wf%n_o**2,             &
                     one,                          &
                     sigma_ai(batch_d%first,1),    &
                     wf%n_v)
!
      enddo
!
      call batch_d%determine_limits(1)
      call mem%dealloc(g_cjmd, wf%n_v, wf%n_o, wf%n_o, batch_d%length)
!
      call wf%g_ckld_t%close_()
!
      call mem%dealloc(Y_cmkj, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
   end subroutine jacobian_transpose_cc3_c3_a1_y_o_cc3
!
!
   module subroutine jacobian_transpose_cc3_c3_b1_y_v_cc3(wf, sigma_ai)
!!
!!    Jacobian transpose contribution Y_vvvo
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!!    sigma_1 += sum_bec g_becd * X_bcek
!!    sigma_1 += sum_cek X_bcek * g_leck
!!    sigma_1 += sum_bek X_bcek * g_lkbe
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
!
      real(dp), dimension(:,:,:,:), allocatable :: g_becd ! vvvv sorted as bce#d
      real(dp), dimension(:,:,:,:), allocatable :: g_ckle ! ovvo sorted as cl#ke
      real(dp), dimension(:,:,:,:), allocatable :: g_celk ! vvoo sorted as cle#k
!
!     arrays for reordering
      real(dp), dimension(:,:,:,:), allocatable :: g_cekl
!
      real(dp), dimension(:,:,:,:), allocatable :: Y_bcek
      real(dp), dimension(:,:,:,:), allocatable :: Y_cbek
!
      type(batching_index) :: batch_k, batch_d
      integer :: k_batch, d_batch
      integer :: req_0, req_k, req_d, req_2
!
!     :: Term 1: sigma_dk += sum_bec g_becd * Y_bcek ::
!
      call wf%g_becd_t%open_('read')
!
      batch_k = batching_index(wf%n_o)
      batch_d = batching_index(wf%n_v)
!
      req_0 = 0
      req_k = wf%n_v**3
      req_d = wf%n_v**3
      req_2 = 0
!
      call mem%batch_setup(batch_k, batch_d, req_0, req_k, req_d, req_2)
!
      call batch_k%determine_limits(1)
      call mem%alloc(Y_bcek, wf%n_v, wf%n_v, wf%n_v, batch_k%length)
!
      call batch_d%determine_limits(1)
      call mem%alloc(g_becd, wf%n_v, wf%n_v, wf%n_v, batch_d%length) 
!
      do k_batch = 1, batch_k%num_batches
!
         call batch_k%determine_limits(k_batch)
!
         call single_record_reader(batch_k, wf%Y_bcek, Y_bcek)
!
         do d_batch = 1, batch_d%num_batches
!
            call batch_d%determine_limits(d_batch)
            call single_record_reader(batch_d, wf%g_becd_t, g_becd)
!
            call dgemm('T','N',                                &
                        batch_d%length,                        &
                        batch_k%length,                        &
                        wf%n_v**3,                             &
                        one,                                   &
                        g_becd,                                & ! g_bce_d
                        wf%n_v**3,                             &
                        Y_bcek,                                & ! Y_bce_k
                        wf%n_v**3,                             &
                        one,                                   &
                        sigma_ai(batch_d%first,batch_k%first), & ! sigma_a_i
                        wf%n_v)
!
         enddo ! d_batch
!
      enddo ! k_batch
!
      call batch_k%determine_limits(1)
      call mem%dealloc(Y_bcek, wf%n_v, wf%n_v, wf%n_v, batch_k%length)
!
      call batch_d%determine_limits(1)
      call mem%dealloc(g_becd, wf%n_v, wf%n_v, wf%n_v, batch_d%length)
!
      call wf%g_becd_t%close_()
!
!
!     :: Term 2: sigma_bl += - sum_cek Y_bcek * g_ckle ::
!
      batch_k = batching_index(wf%n_o)
!
      req_0 = 0
      req_k = wf%n_v**3 + 2*wf%n_o*wf%n_v**2
!
      call mem%batch_setup(batch_k, req_0, req_k)
!
      call batch_k%determine_limits(1)
      call mem%alloc(Y_bcek, wf%n_v, wf%n_v, wf%n_v, batch_k%length)
!
      call wf%g_ckld_t%open_('read')
!
      do k_batch = 1, batch_k%num_batches
!
         call batch_k%determine_limits(k_batch)
!
         call mem%alloc(g_ckle, wf%n_v, wf%n_o, batch_k%length, wf%n_v) ! read in cl#ke
         call mem%alloc(g_cekl, wf%n_v, wf%n_v, batch_k%length, wf%n_o) ! sort to ce#kl
!
         call compound_record_reader(wf%n_v, batch_k, wf%g_ckld_t, g_ckle, .true.) ! stored as cl#k#e
         call sort_1234_to_1432(g_ckle, g_cekl, wf%n_v, wf%n_o, batch_k%length, wf%n_v) ! cl#ke -> ce#kl
!
         call single_record_reader(batch_k, wf%Y_bcek, Y_bcek)
!
         call dgemm('N','N',                    &
                     wf%n_v,                    &
                     wf%n_o,                    &
                     batch_k%length*wf%n_v**2,  &
                     -one,                      &
                     Y_bcek,                    & ! Y_b_cek
                     wf%n_v,                    &
                     g_cekl,                    & ! g_cek_l
                     batch_k%length*wf%n_v**2,  &
                     one,                       &
                     sigma_ai,                  &
                     wf%n_v)
!
         call mem%dealloc(g_ckle, wf%n_v, wf%n_o, batch_k%length, wf%n_v) ! cl#ke
         call mem%dealloc(g_cekl, wf%n_v, wf%n_v, batch_k%length, wf%n_o) ! ce#kl
!
      enddo
!
      call batch_k%determine_limits(1)
      call mem%dealloc(Y_bcek, wf%n_v, wf%n_v, wf%n_v, batch_k%length)
!
      call wf%g_ckld_t%close_()
!
!     
!     :: Term 3: sigma_bl += - sum_cek Y_cbek * g_celk ::
!
      call wf%g_cdlk_t%open_('read')
!
      req_0 = 0
      req_k = 2*wf%n_v**3 + 2*wf%n_o*wf%n_v**2
!
      call mem%batch_setup(batch_k, req_0, req_k)
!
      call batch_k%determine_limits(1)
      call mem%alloc(Y_bcek, wf%n_v, wf%n_v, wf%n_v, batch_k%length)
      call mem%alloc(Y_cbek, wf%n_v, wf%n_v, wf%n_v, batch_k%length)
!
      do k_batch = 1, batch_k%num_batches
!
         call batch_k%determine_limits(k_batch)
!
         call mem%alloc(g_celk, wf%n_v, wf%n_o, batch_k%length, wf%n_v) ! read in cl#ke
         call mem%alloc(g_cekl, wf%n_v, wf%n_v, batch_k%length, wf%n_o) ! sort to ce#kl
!
         call single_record_reader(batch_k, wf%Y_bcek, Y_bcek)
         call sort_1234_to_2134(Y_bcek, Y_cbek, wf%n_v, wf%n_v, wf%n_v, batch_k%length)
!
         call compound_record_reader(wf%n_v, batch_k, wf%g_cdlk_t, g_celk, .true.) ! stored as cl#k#e
         call sort_1234_to_1432(g_celk, g_cekl, wf%n_v, wf%n_o, batch_k%length, wf%n_v) ! cl#ke -> ce#kl
!
         call dgemm('N','N',                    &
                     wf%n_v,                    &
                     wf%n_o,                    &
                     batch_k%length*wf%n_v**2,  &
                     -one,                      &
                     Y_cbek,                    & ! Y_b_cek
                     wf%n_v,                    &
                     g_cekl,                    & ! g_cek_l
                     batch_k%length*wf%n_v**2,  &
                     one,                       &
                     sigma_ai,                  &
                     wf%n_v)
!
         call mem%dealloc(g_celk, wf%n_v, wf%n_o, batch_k%length, wf%n_v) ! cl#ke
         call mem%dealloc(g_cekl, wf%n_v, wf%n_v, batch_k%length, wf%n_o) ! ce#kl
!
      enddo
!
      call batch_k%determine_limits(1)
      call mem%dealloc(Y_bcek, wf%n_v, wf%n_v, wf%n_v, batch_k%length)
      call mem%dealloc(Y_cbek, wf%n_v, wf%n_v, wf%n_v, batch_k%length)
!
      call wf%g_cdlk_t%close_()
!
   end subroutine jacobian_transpose_cc3_c3_b1_y_v_cc3
!
!
end submodule jacobian_transpose
