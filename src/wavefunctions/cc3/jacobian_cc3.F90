!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2020 the authors of eT
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
submodule (cc3_class) jacobian
!
!!
!!    Jacobian submodule
!!
!!    Routines for the linear transform of trial
!!    vectors by the Jacobian matrix
!!
!!    rho_i = A * c_i,
!!
!!    where
!!
!!    A_mu,nu = < mu| exp(-T) [H, tau_nu] exp(T) |R >.
!!
!
   implicit none
!
!
contains
!
!
   module subroutine effective_jacobian_transformation_cc3(wf, omega, c)
!!
!!    Effective Jacobian transformation (CC3)
!!    Written by Alexander C. Paul and Rolf H. Myhre, Feb 2019
!!
!!    Directs the transformation by the CC3 Jacobi matrix,
!!
!!       A_mu,nu = < mu| exp(-T) [H, tau_nu] exp(T) | R >,
!!
!!    where the basis employed for the brackets is biorthonormal.
!!    The transformation is rho = A c, i.e.,
!!
!!       rho_mu = (A c)_mu 
!!       = sum_ck A_mu,ck c_ck + 1/2 sum_ckdl A_mu,ckdl c_ckdl (1 + delta_ck,dl)
!!
!!    On exit, c is overwritten by rho. That is, c(ai) = rho_a_i,
!!    and c(aibj) = rho_aibj.
!!
      use array_utilities, only: scale_diagonal
!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), intent(in) :: omega
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: c
!
      real(dp), dimension(:,:), allocatable :: c_ai
      real(dp), dimension(:,:,:,:), allocatable :: c_aibj, c_abij
!
      real(dp), dimension(:,:), allocatable :: rho_ai
      real(dp), dimension(:,:,:,:), allocatable :: rho_aibj, rho_abij
!
      type(timings), allocatable :: cc3_timer, ccsd_timer, timer
!
      timer           = timings('Jacobian CC3')
      cc3_timer       = timings('Jacobian CC3 (CC3 contribution)', pl='normal')
      ccsd_timer      = timings('Jacobian CC3 (CCSD contribution)', pl='normal')
!
      call timer%turn_on()
!
!     Allocate and zero the transformed vector (singles part)
!
      call mem%alloc(rho_ai, wf%n_v, wf%n_o)
      call zero_array(rho_ai, wf%n_v*wf%n_o)
!
      call mem%alloc(c_ai, wf%n_v, wf%n_o)
!
      call dcopy(wf%n_t1, c(1:wf%n_t1), 1, c_ai, 1)
!
!     :: CCS contributions to the singles c vector ::
!
      call ccsd_timer%turn_on()
!
      call wf%jacobian_ccs_a1(rho_ai, c_ai)
      call wf%jacobian_ccs_b1(rho_ai, c_ai)
!
!     :: CCSD contributions to the transformed singles vector ::
!
      call wf%jacobian_doubles_a1(rho_ai, c_ai)
!
!     Allocate the incoming unpacked doubles vector
!
      call mem%alloc(c_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call squareup(c(wf%n_t1 + 1 : wf%n_es_amplitudes), c_aibj, wf%n_t1)
!
!     Scale the doubles vector by 1 + delta_ai,bj, i.e.
!     redefine to c_ckdl = c_ckdl (1 + delta_ck,dl)
!
      call scale_diagonal(two, c_aibj, wf%n_t1)
!
      call wf%jacobian_doubles_b1(rho_ai, c_aibj)
      call wf%jacobian_doubles_c1(rho_ai, c_aibj)
      call wf%jacobian_doubles_d1(rho_ai, c_aibj)
!
!     :: CCSD contributions to the transformed doubles vector ::
!     Allocate unpacked transformed vector
!
      call mem%alloc(rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call zero_array(rho_aibj, (wf%n_v*wf%n_o)**2)
!
!     Contributions from singles vector c
!
      call wf%jacobian_doubles_a2(rho_aibj, c_ai)
      call wf%jacobian_ccsd_b2(rho_aibj, c_ai)
      call wf%jacobian_ccsd_c2(rho_aibj, c_ai)
      call wf%jacobian_ccsd_d2(rho_aibj, c_ai)
!
!     Contributions from doubles vector c
!
      call wf%jacobian_ccsd_e2(rho_aibj, c_aibj)
      call wf%jacobian_ccsd_f2(rho_aibj, c_aibj)
      call wf%jacobian_ccsd_g2(rho_aibj, c_aibj)
      call wf%jacobian_ccsd_h2(rho_aibj, c_aibj)
      call wf%jacobian_ccsd_i2(rho_aibj, c_aibj)
!
!     Compute CC3 contributions to rho_ai and rho_aibj and symmetrise rho_aibj
!     CCSD J2 and K2 are already symmetric and will be computed afterwards. The CCSD L2 
!     term is also already symmetric, but is computed as aibj and is computed before 
!     reordering, where the half*c2 cancels out in later symmetrization of ai-bj.
!
      call mem%alloc(rho_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(c_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call sort_1234_to_1324(c_aibj, c_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call sort_1234_to_1324(rho_aibj, rho_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(c_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call ccsd_timer%freeze()
!
      call cc3_timer%turn_on()
!
!     CC3-Contributions from the T3-amplitudes
      call wf%jacobian_cc3_t3_a2(c_ai, rho_abij)
!
      call wf%jacobian_cc3_t3_b2(c_ai, rho_abij)
!
!     CC3-Contributions from the C3-amplitudes
!
      call wf%jacobian_cc3_c3_a(omega, c_ai, c_abij, rho_ai, rho_abij)
!
      call cc3_timer%freeze()
!
!     Done with singles vector c; Overwrite the incoming singles c vector for exit
!
      call mem%dealloc(c_ai, wf%n_v, wf%n_o)
!
      call dcopy(wf%n_t1, rho_ai, 1, c, 1)
!
      call mem%dealloc(rho_ai, wf%n_v, wf%n_o)
!
!     Last two CCSD-terms (J2, K2) are already symmetric.
!     Perform the symmetrization rho_ai_bj = P_ij^ab rho_ai_bj
!
      call ccsd_timer%turn_on()
!
      call symmetrize_12_and_34(rho_abij, wf%n_v, wf%n_o)
!
!     Compute CCSD J2 and K2 contributions
!
      call wf%jacobian_ccsd_j2(rho_abij, c_abij)
      call wf%jacobian_ccsd_k2(rho_abij, c_abij)
      call wf%omega_ccsd_a2(rho_abij, c_abij)
!
      call mem%dealloc(c_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call ccsd_timer%freeze()
!
!     divide by the biorthonormal factor 1 + delta_ai,bj
! 
      call scale_diagonal(half, rho_abij, wf%n_v, wf%n_o)
!
!     overwrite the incoming, packed doubles c vector with rho_abij
!     and pack in for exit
!
      call packin(c(wf%n_t1 + 1 : wf%n_es_amplitudes), rho_abij, wf%n_v, wf%n_o)
!
      call mem%dealloc(rho_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call timer%turn_off()
      call cc3_timer%turn_off()
      call ccsd_timer%turn_off()
!
   end subroutine effective_jacobian_transformation_cc3
!
!
   module subroutine jacobian_cc3_t3_a2_cc3(wf, c_ai, rho_abij)
!!
!!    Jacobian CC3 A2
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!!    Reads in the intermediates X_abid and X_ajil prepared in 
!!    prepare_jacobian_transform contracts with c_ai and adds to rho_abij
!!
!!    rho_abil += sum_abi X_abid * C_dl
!!    rho_daji += sum_aik C_dl * X_ajil
!!
!!    where: X_abid = - sum_jck (2 t^abc_ijk - t^cba_ijk - t^acb_ijk) * g_kcjd
!!           X_ajil = - sum_bck (2 t^abc_ijk - t^cba_ijk - t^acb_ijk) * g_lbkc
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(out) :: rho_abij
!
      real(dp), dimension(:,:,:,:), allocatable :: X_abid
      real(dp), dimension(:,:,:,:), allocatable :: X_ajil
!
      type(batching_index) :: batch_d
      integer :: d_batch
      integer :: req_0, req_d
!
      type(timings), allocatable :: cc3_timer_t3_a2
!
      cc3_timer_t3_a2 = timings('Time in CC3 T3 a2', pl='verbose') 
      call cc3_timer_t3_a2%turn_on()
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
         call wf%X_abid%read_compound_full_batch(X_abid, wf%n_o, batch_d)
!
         call dgemm('N','N',                 &
                     wf%n_o*wf%n_v**2,       &
                     wf%n_o,                 &
                     batch_d%length,         &
                     one,                    &
                     X_abid,                 & ! X_abi_d
                     wf%n_o*wf%n_v**2,       &
                     c_ai(batch_d%first, 1), & ! c_d_l
                     wf%n_v,                 &
                     one,                    &
                     rho_abij,               & ! rho_abil
                     wf%n_o*wf%n_v**2)
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
      call wf%X_ajil%read_4(X_ajil, 1, wf%n_o)
!
      call wf%X_ajil%close_()
!
      call dgemm('N','T',           &
                  wf%n_v,           &
                  wf%n_v*wf%n_o**2, &
                  wf%n_o,           &
                  one,              &
                  c_ai,             & ! C_d_l
                  wf%n_v,           &
                  X_ajil,           & ! X_aji_l
                  wf%n_v*wf%n_o**2, &
                  one,              &
                  rho_abij,         & ! rho_daji
                  wf%n_v)
!
      call mem%dealloc(X_ajil, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call cc3_timer_t3_a2%turn_off()
!
   end subroutine jacobian_cc3_t3_a2_cc3
!
!
   module subroutine jacobian_cc3_t3_b2_cc3(wf, c_ai, rho_abij)
!!
!!    Jacobian CC3 B2
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!!    sigma_abij +=  sum_ckdl C^d_l L_kcld (t^abc_ijk - t^bac_ijk)
!!               +=  sum_ck F_kc_c1 * (t^abc_ijk - t^bac_ijk)
!!
!!    Constructs t^abc_ijk for fixed ijk contracts with
!!    the c1-transformed Fock matrix
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout) :: rho_abij
!
!     Arrays for triples amplitudes
      real(dp), dimension(:,:,:), allocatable :: t_abc
      real(dp), dimension(:,:,:), allocatable :: u_abc
!
!     Unpacked doubles amplitudes
      real(dp), dimension(:,:,:,:), allocatable :: t_abij
!
      real(dp), dimension(:,:), allocatable :: F_ov_ck_c1
!
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
      integer              :: i, j, k, i_rel, j_rel, k_rel
      type(batching_index) :: batch_i, batch_j, batch_k
      integer              :: i_batch, j_batch, k_batch
      integer              :: req_0, req_1, req_2, req_3
!
      logical :: ijk_core
!
      type(timings), allocatable :: cc3_timer_t3_b2
!
      cc3_timer_t3_b2 = timings('Time in CC3 T3 b2', pl='verbose')
      call cc3_timer_t3_b2%turn_on()
!
      call mem%alloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call squareup_and_sort_1234_to_1324(wf%t2, t_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     C1 transformed Fock matrix
!
      call mem%alloc(F_ov_ck_c1, wf%n_v, wf%n_o)
      call wf%construct_c1_fock(c_ai, F_ov_ck_c1)
!
!     Setup and Batching loops for the T3-contributions to rho2
!
      req_0 = 2*wf%n_v**3
      req_1 = (wf%n_v)**3
      req_2 = (wf%n_o)*(wf%n_v)
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
         call mem%alloc(g_ljci, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
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
      call wf%g_bdck_t%open_('read')
      call wf%g_ljck_t%open_('read')
!
      do i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(i_batch)
!
         call wf%g_bdck_t%read_interval(g_bdci, batch_i)
         g_bdci_p => g_bdci
!
         do j_batch = 1, i_batch
!
            call batch_j%determine_limits(j_batch)
!
            call wf%g_ljck_t%read_compound(g_ljci, batch_j, batch_i)
            g_ljci_p => g_ljci
!
            if (j_batch .ne. i_batch) then
!
               call wf%g_bdck_t%read_interval(g_bdcj, batch_j)
               g_bdcj_p => g_bdcj
!
               call wf%g_ljck_t%read_compound(g_licj, batch_i, batch_j)
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
               if (k_batch .ne. j_batch) then !k_batch != j_batch, k_batch != i_batch
!
                  call wf%g_bdck_t%read_interval(g_bdck, batch_k)
                  g_bdck_p => g_bdck
! 
                  call wf%g_ljck_t%read_compound(g_lkci, batch_k, batch_i)
                  g_lkci_p => g_lkci
!
                  call wf%g_ljck_t%read_compound(g_lick, batch_i, batch_k)
                  g_lick_p => g_lick
!
                  call wf%g_ljck_t%read_compound(g_lkcj, batch_k, batch_j)
                  g_lkcj_p => g_lkcj
!
                  call wf%g_ljck_t%read_compound(g_ljck, batch_j, batch_k)
                  g_ljck_p => g_ljck
!
!
               else if (k_batch .eq. i_batch) then !k_batch == j_batch == i_batch
!
                  g_bdck_p => g_bdci
!
                  g_lkci_p => g_ljci
                  g_lick_p => g_ljci
!
                  g_lkcj_p => g_ljci
                  g_ljck_p => g_ljci
!
               else !k_batch == j_batch != i_batch
!
                  g_bdck_p => g_bdcj
!
                  g_lkci_p => g_ljci
                  g_lick_p => g_licj
!
                  call wf%g_ljck_t%read_compound(g_lkcj, batch_k, batch_j)
                  g_lkcj_p => g_lkcj
                  g_ljck_p => g_lkcj
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
                        if (k .eq. i) then ! k == j == i
                           cycle
                        end if
!
!                       Check if at least one index i,j,k is a core orbital
!                       Here t3 contributes to rho2 and can, thus, be restricted as well
!
                        if(wf%cvs) then
!
                           ijk_core = .false.
!
                           if(     any(wf%core_MOs .eq. i) &
                              .or. any(wf%core_MOs .eq. j) &
                              .or. any(wf%core_MOs .eq. k)) then
!
                              ijk_core = .true.
!
                           end if
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
                        call wf%jacobian_cc3_b2_fock(i, j, k, t_abc, u_abc, rho_abij, F_ov_ck_c1)
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
      if (batch_i%num_batches .eq. 1) then !no batching
!
         call mem%dealloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%dealloc(g_ljci, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
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
      call mem%dealloc(F_ov_ck_c1, wf%n_v, wf%n_o)
!
      call mem%dealloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call cc3_timer_t3_b2%turn_off()
!
   end subroutine jacobian_cc3_t3_b2_cc3
!
!
   module subroutine jacobian_cc3_c3_a_cc3(wf, omega, c_ai, c_abij, rho_ai, rho_abij)
!!
!!    Contribution of the C3/R3 terms
!!    Written by Alexander C. Paul and Rolf H. Myhre, Feb 2019
!!
!!    Construct C^abc_ijk in single batches of ijk and compute the contributions
!!    to the singles and doubles part of the outgoing vector
!!
!!    The triples amplitudes are expressed in terms of doubles amplitudes:
!!    C_3 = (omega - epsilon_mu3)^-1 (< mu3| [H,C_2] | HF > 
!!                                  + < mu3| [[H,C_1],T_2] |HF >)
!!
!!    c^abc = (omega - epsilon^abc_ijk)^-1 * P^abc_ijk 
!!             (sum_d c^ad_ij g_ckbd - sum_l c^ab_il g_cklj
!!            + sum_d t^ad_ij g'_bdck - sum_l t^ab_il g'_cklj
!!
!!    rho1 += < mu1| [H,C_3] |R >
!!    rho2 += < mu2| [H,C_3] |R >
!!
!!    Based on omega_cc3_a_cc3 written by Rolf H. Myhre
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
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout) :: rho_abij
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
!     T1 transformed integrals
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
      real(dp), dimension(:,:,:,:), allocatable, target  :: L_jbic
      real(dp), dimension(:,:,:,:), allocatable, target  :: L_kbic
      real(dp), dimension(:,:,:,:), allocatable, target  :: L_kbjc
      real(dp), dimension(:,:,:,:), allocatable, target  :: L_ibjc
      real(dp), dimension(:,:,:,:), allocatable, target  :: L_ibkc
      real(dp), dimension(:,:,:,:), allocatable, target  :: L_jbkc
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: L_jbic_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: L_kbic_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: L_kbjc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: L_ibjc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: L_ibkc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: L_jbkc_p => null()
!
!     C1 transformed integrals
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdci_c1
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdcj_c1
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdck_c1
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdci_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdcj_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdck_c1_p => null()
!
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ljci_c1
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_lkci_c1
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_lkcj_c1
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_licj_c1
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_lick_c1
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ljck_c1
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_ljci_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_lkci_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_lkcj_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_licj_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_lick_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_ljck_c1_p => null()
!
      integer              :: i, j, k, i_rel, j_rel, k_rel
      type(batching_index) :: batch_i, batch_j, batch_k
      integer              :: i_batch, j_batch, k_batch
      integer              :: req_0, req_1, req_2, req_3
!
      logical :: ijk_core
!
      type(timings), allocatable :: cc3_timer_c3
!
      cc3_timer_c3    = timings('Time in CC3 C3', pl='verbose')
      call cc3_timer_c3%turn_on()
!
!     Set up required c1-transformed integrals
      call wf%construct_c1_integrals(c_ai)
!
!     Set up arrays for amplitudes
      call mem%alloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call squareup_and_sort_1234_to_1324(wf%t2, t_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      batch_i = batching_index(wf%n_o)
      batch_j = batching_index(wf%n_o)
      batch_k = batching_index(wf%n_o)
!
!     Setup and Batching loops for the C3-contributions to rho1 and rho2
!
      req_0 = 3*wf%n_v**3 + wf%n_o*wf%n_v
      req_1 = 3*(wf%n_v)**3
      req_2 = 3*(wf%n_o)*(wf%n_v) + (wf%n_v)**2
      req_3 = 0
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
         call mem%alloc(g_ljci, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
         call mem%alloc(g_bdci_c1, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%alloc(g_ljci_c1, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      else ! batching
!
         call batch_i%determine_limits(1)
!
!        Ordered such that batching indices are at the end
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
         call mem%alloc(g_bdci_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%alloc(g_bdcj_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%alloc(g_bdck_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
         call mem%alloc(g_ljci_c1, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_lkci_c1, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_lkcj_c1, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_licj_c1, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_lick_c1, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_ljck_c1, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
!
      endif
!
!     Resorting of the Fock-Matrix for easier contractions later
      call mem%alloc(F_ov_ck, wf%n_v, wf%n_o)
      call sort_12_to_21(wf%fock_ia, F_ov_ck, wf%n_o, wf%n_v)
!
!     Arrays for the triples amplitudes and intermediates
      call mem%alloc(c_abc, wf%n_v, wf%n_v, wf%n_v)
      call mem%alloc(u_abc, wf%n_v, wf%n_v, wf%n_v)
      call mem%alloc(v_abc, wf%n_v, wf%n_v, wf%n_v)
!
      call zero_array(c_abc, wf%n_v**3)
!
!     Remaining integral arrays
!
      if (batch_i%num_batches .eq. 1) then ! no batching
!
         call mem%alloc(g_dbic, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%alloc(g_jlic, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         call mem%alloc(L_jbic, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
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
         call mem%alloc(L_jbic, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(L_kbic, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(L_kbjc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(L_ibjc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(L_ibkc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(L_jbkc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
!
      endif
!
      call wf%g_bdck_t%open_('read')
      call wf%g_ljck_t%open_('read')
      call wf%g_dbkc_t%open_('read')
      call wf%g_jlkc_t%open_('read')
      call wf%L_jbkc_t%open_('read')
!
      call wf%g_bdck_c%open_('read')
      call wf%g_ljck_c%open_('read')
!
      do i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(i_batch)
!
         call wf%g_bdck_t%read_interval(g_bdci, batch_i)
         call wf%g_dbkc_t%read_interval(g_dbic, batch_i)
         call wf%g_bdck_c%read_interval(g_bdci_c1, batch_i)
!
         g_bdci_p => g_bdci
         g_dbic_p => g_dbic
         g_bdci_c1_p => g_bdci_c1
!
         do j_batch = 1, i_batch
!
            call batch_j%determine_limits(j_batch)
!
            call wf%g_ljck_t%read_compound(g_ljci, batch_j, batch_i)
            call wf%g_jlkc_t%read_compound(g_jlic, batch_j, batch_i)
            call wf%L_jbkc_t%read_compound(L_jbic, batch_j, batch_i)
            call wf%g_ljck_c%read_compound(g_ljci_c1, batch_j, batch_i)
!
            g_ljci_p => g_ljci
            g_jlic_p => g_jlic
            L_jbic_p => L_jbic
!
            g_ljci_c1_p => g_ljci_c1
!
            if (j_batch .ne. i_batch) then
!
               call wf%g_bdck_t%read_interval(g_bdcj, batch_j)
               call wf%g_dbkc_t%read_interval(g_dbjc, batch_j)
               call wf%g_bdck_c%read_interval(g_bdcj_c1, batch_j)
!
               g_bdcj_p => g_bdcj
               g_dbjc_p => g_dbjc
               g_bdcj_c1_p => g_bdcj_c1
!
               call wf%g_ljck_t%read_compound(g_licj, batch_i, batch_j)
               call wf%g_jlkc_t%read_compound(g_iljc, batch_i, batch_j)
               call wf%L_jbkc_t%read_compound(L_ibjc, batch_i, batch_j)
               call wf%g_ljck_c%read_compound(g_licj_c1, batch_i, batch_j)
!
               g_licj_p => g_licj
               g_iljc_p => g_iljc
               L_ibjc_p => L_ibjc
!
               g_licj_c1_p => g_licj_c1
!
            else
!
               g_bdcj_p => g_bdci
               g_dbjc_p => g_dbic
!
               g_licj_p => g_ljci
               g_iljc_p => g_jlic
               L_ibjc_p => L_jbic
!
               g_bdcj_c1_p => g_bdci_c1
!
               g_licj_c1_p => g_ljci_c1
!
            endif
!
            do k_batch = 1, j_batch
!
               call batch_k%determine_limits(k_batch)
!
               if (k_batch .ne. j_batch) then ! k_batch != j_batch, k_batch != i_batch
!
                  call wf%g_bdck_t%read_interval(g_bdck, batch_k)
                  call wf%g_dbkc_t%read_interval(g_dbkc, batch_k)
                  call wf%g_bdck_c%read_interval(g_bdck_c1, batch_k)
!
                  g_bdck_p => g_bdck
                  g_dbkc_p => g_dbkc
                  g_bdck_c1_p => g_bdck_c1
! 
                  call wf%g_ljck_t%read_compound(g_lkci, batch_k, batch_i)
                  call wf%g_jlkc_t%read_compound(g_klic, batch_k, batch_i)
                  call wf%L_jbkc_t%read_compound(L_kbic, batch_k, batch_i)
                  call wf%g_ljck_c%read_compound(g_lkci_c1, batch_k, batch_i)
!
                  g_lkci_p => g_lkci
                  g_klic_p => g_klic
                  L_kbic_p => L_kbic
                  g_lkci_c1_p => g_lkci_c1
!
                  call wf%g_ljck_t%read_compound(g_lick, batch_i, batch_k)
                  call wf%g_jlkc_t%read_compound(g_ilkc, batch_i, batch_k)
                  call wf%L_jbkc_t%read_compound(L_ibkc, batch_i, batch_k)
                  call wf%g_ljck_c%read_compound(g_lick_c1, batch_i, batch_k)
!
                  g_lick_p => g_lick
                  g_ilkc_p => g_ilkc
                  L_ibkc_p => L_ibkc
                  g_lick_c1_p => g_lick_c1
!
                  call wf%g_ljck_t%read_compound(g_lkcj, batch_k, batch_j)
                  call wf%g_jlkc_t%read_compound(g_kljc, batch_k, batch_j)
                  call wf%L_jbkc_t%read_compound(L_kbjc, batch_k, batch_j)
                  call wf%g_ljck_c%read_compound(g_lkcj_c1, batch_k, batch_j)
!
                  g_lkcj_p => g_lkcj
                  g_kljc_p => g_kljc
                  L_kbjc_p => L_kbjc
                  g_lkcj_c1_p => g_lkcj_c1
!
                  call wf%g_ljck_t%read_compound(g_ljck, batch_j, batch_k)
                  call wf%g_jlkc_t%read_compound(g_jlkc, batch_j, batch_k)
                  call wf%L_jbkc_t%read_compound(L_jbkc, batch_j, batch_k)
                  call wf%g_ljck_c%read_compound(g_ljck_c1, batch_j, batch_k)
!
                  g_ljck_p => g_ljck
                  g_jlkc_p => g_jlkc
                  L_jbkc_p => L_jbkc
                  g_ljck_c1_p => g_ljck_c1
!
               else if (k_batch .eq. i_batch) then ! k_batch == j_batch == i_batch
!
                  g_bdck_p => g_bdci
                  g_dbkc_p => g_dbic
                  g_bdck_c1_p => g_bdci_c1
!
                  g_lkci_p => g_ljci
                  g_klic_p => g_jlic
                  L_kbic_p => L_jbic
                  g_lkci_c1_p => g_ljci_c1
!
                  g_lick_p => g_ljci
                  g_ilkc_p => g_jlic
                  L_ibkc_p => L_jbic
                  g_lick_c1_p => g_ljci_c1
!
                  g_lkcj_p => g_ljci
                  g_kljc_p => g_jlic
                  L_kbjc_p => L_jbic
                  g_lkcj_c1_p => g_ljci_c1
!
                  g_ljck_p => g_ljci
                  g_jlkc_p => g_jlic
                  L_jbkc_p => L_jbic
                  g_ljck_c1_p => g_ljci_c1
!
               else ! k_batch == j_batch != i_batch
!
                  g_bdck_p => g_bdcj
                  g_dbkc_p => g_dbjc
                  g_bdck_c1_p => g_bdcj_c1
!
                  g_lkci_p => g_ljci
                  g_klic_p => g_jlic
                  L_kbic_p => L_jbic
                  g_lkci_c1_p => g_ljci_c1
!
                  g_lick_p => g_licj
                  g_ilkc_p => g_iljc
                  L_ibkc_p => L_ibjc
                  g_lick_c1_p => g_licj_c1
!
                  call wf%g_ljck_t%read_compound(g_lkcj, batch_k, batch_j)
                  call wf%g_jlkc_t%read_compound(g_kljc, batch_k, batch_j)
                  call wf%L_jbkc_t%read_compound(L_kbjc, batch_k, batch_j)
                  call wf%g_ljck_c%read_compound(g_lkcj_c1, batch_k, batch_j)
!
                  g_lkcj_p => g_lkcj
                  g_kljc_p => g_kljc
                  L_kbjc_p => L_kbjc
                  g_lkcj_c1_p => g_lkcj_c1
!
                  g_ljck_p => g_lkcj
                  g_jlkc_p => g_kljc
                  L_jbkc_p => L_kbjc
                  g_ljck_c1_p => g_lkcj_c1
!
               endif
!
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
!                       Check if at least one index i,j,k is a core orbital
!
                        if(wf%cvs) then
!
                           ijk_core = .false.
!
                           if(     any(wf%core_MOs .eq. i) &
                              .or. any(wf%core_MOs .eq. j) &
                              .or. any(wf%core_MOs .eq. k)) then
!
                              ijk_core = .true.
!
                           end if
!
!                          Cycle if i,j,k are not core orbitals
                           if (.not. ijk_core) cycle
!
                        end if
!
                        k_rel = k - batch_k%first + 1
!
!                       Construct C^{abc}_{ijk} for given i, j, k
!                       and calculate contributions to rho1 and rho2
!                       Using c1-transformed integrals the terms have the same form 
!                       as the omega terms (where t_abc = c_abc)
!
!                       Therefore the contributions to the c3-amplitudes can be computed 
!                       using the same routine once for t1-transformed and once for 
!                       c1-transformed integrals
!
                        call wf%construct_W(i, j, k, c_abc, u_abc, c_abij,  &
                                            g_bdci_p(:,:,:,i_rel),          &
                                            g_bdcj_p(:,:,:,j_rel),          &
                                            g_bdck_p(:,:,:,k_rel),          &
                                            g_ljci_p(:,:,j_rel,i_rel),      &
                                            g_lkci_p(:,:,k_rel,i_rel),      &
                                            g_lkcj_p(:,:,k_rel,j_rel),      &
                                            g_licj_p(:,:,i_rel,j_rel),      &
                                            g_lick_p(:,:,i_rel,k_rel),      &
                                            g_ljck_p(:,:,j_rel,k_rel))
!
                        call wf%construct_W(i, j, k, c_abc, u_abc, t_abij,  &
                                            g_bdci_c1_p(:,:,:,i_rel),       &
                                            g_bdcj_c1_p(:,:,:,j_rel),       &
                                            g_bdck_c1_p(:,:,:,k_rel),       &
                                            g_ljci_c1_p(:,:,j_rel,i_rel),   &
                                            g_lkci_c1_p(:,:,k_rel,i_rel),   &
                                            g_lkcj_c1_p(:,:,k_rel,j_rel),   &
                                            g_licj_c1_p(:,:,i_rel,j_rel),   &
                                            g_lick_c1_p(:,:,i_rel,k_rel),   &
                                            g_ljck_c1_p(:,:,j_rel,k_rel),   &
                                            overwrite = .false.) ! overwrite R_abc
!
                        call wf%divide_by_orbital_differences(i, j, k, c_abc, omega)
!
                        call wf%omega_cc3_a_n7(i, j, k, c_abc, u_abc, v_abc,  &
                                               rho_abij,                      &
                                               g_dbic_p(:,:,:,i_rel),         &
                                               g_dbjc_p(:,:,:,j_rel),         &
                                               g_dbkc_p(:,:,:,k_rel),         &
                                               g_jlic_p(:,:,j_rel,i_rel),     &
                                               g_klic_p(:,:,k_rel,i_rel),     &
                                               g_kljc_p(:,:,k_rel,j_rel),     &
                                               g_iljc_p(:,:,i_rel,j_rel),     &
                                               g_ilkc_p(:,:,i_rel,k_rel),     &
                                               g_jlkc_p(:,:,j_rel,k_rel))
!
                        call wf%omega_cc3_a_n6(i, j, k, c_abc, u_abc,      &
                                               rho_ai, rho_abij, F_ov_ck,  &
                                               L_jbic_p(:,:,j_rel,i_rel),  &
                                               L_kbic_p(:,:,k_rel,i_rel),  &
                                               L_kbjc_p(:,:,k_rel,j_rel),  &
                                               L_ibjc_p(:,:,i_rel,j_rel),  &
                                               L_ibkc_p(:,:,i_rel,k_rel),  &
                                               L_jbkc_p(:,:,j_rel,k_rel))
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
      call wf%g_dbkc_t%close_()
      call wf%g_jlkc_t%close_()
      call wf%L_jbkc_t%close_()
!
      call wf%g_bdck_c%close_()
      call wf%g_ljck_c%close_()
!
!     Deallocate the integral arrays
!
      if (batch_i%num_batches .eq. 1) then ! no batching
!
         call mem%dealloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%dealloc(g_dbic, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%dealloc(g_ljci, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
         call mem%dealloc(g_jlic, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         call mem%dealloc(L_jbic, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
         call mem%dealloc(g_bdci_c1, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%dealloc(g_ljci_c1, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      else ! batching
!
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
         call mem%dealloc(L_jbic, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(L_kbic, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(L_kbjc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(L_ibjc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(L_ibkc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(L_jbkc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
!
         call mem%dealloc(g_bdci_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%dealloc(g_bdcj_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%dealloc(g_bdck_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
         call mem%dealloc(g_ljci_c1, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_lkci_c1, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_lkcj_c1, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_licj_c1, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_lick_c1, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_ljck_c1, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
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
      call cc3_timer_c3%turn_off()
!
   end subroutine jacobian_cc3_c3_a_cc3
!
!
   module subroutine construct_c1_integrals_cc3(wf, c_ai)
!!
!!    Construct c1 transformed integrals
!!    Alexander C. Paul and Rolf H. Myhre Feb 2019
!!
!!    g'_bdck = (b'd|ck) + (bd|c'k) + (bd|ck')   ordered as dbc,k
!!    g'_ljck = (lj'|ck) + (lj|ck') + (lj|c'k)   ordered as lc,jk
!!
!!    Based on omega_cc3_integrals_cc3 written by Rolf H. Myhre
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
!
      real(dp), dimension(:,:,:,:), allocatable :: g_pqrs ! Array for constructed integrals
      real(dp), dimension(:,:,:,:), allocatable :: h_pqrs ! Array for sorted integrals
!
      real(dp), dimension(:,:,:), allocatable :: L_J_ck_c1              ! c1 transformed Cholesky vector
      real(dp), dimension(:,:,:), allocatable :: L_J_ck, L_J_bd, L_J_lj ! Cholesky vectors
!
      type(batching_index) :: batch_k
!
      integer :: req_0, req_k
      integer :: k_batch
!
      batch_k = batching_index(wf%n_o)
!
!     g'_bdck = (b'd|ck) + (bd|c'k) + (bd|ck')
!
      req_0 = 2*(wf%n_v)**2*(wf%integrals%n_J) + (wf%n_v)*(wf%n_o)*(wf%integrals%n_J)
      req_k = max((wf%integrals%n_J)*(wf%n_v) + (wf%n_v)**3, 2*(wf%n_v)**3)
!
      call mem%alloc(L_J_bd, wf%integrals%n_J, wf%n_v, wf%n_v)
!
      call mem%batch_setup(batch_k, req_0, req_k)
!
      wf%g_bdck_c = direct_stream_file('g_bdck_c', wf%n_v**3)
      call wf%g_bdck_c%open_('write')
!
      do k_batch = 1, batch_k%num_batches
!
         call batch_k%determine_limits(k_batch)
!
!        :: Term 1: g_b'dck += sum_J L_J_b'd L_J_ck ::
!
!        here L_J_bd is used for the c1-transformed cholesky vector
         call wf%integrals%construct_cholesky_ab_c1(L_J_bd, c_ai, 1, wf%n_v, 1, wf%n_v)
!
         call mem%alloc(L_J_ck, wf%integrals%n_J, wf%n_v, batch_k%length)
         call wf%integrals%get_cholesky_t1(L_J_ck, wf%n_o + 1, wf%n_o + wf%n_v, batch_k%first, batch_k%last)
!
         call mem%alloc(g_pqrs, wf%n_v, wf%n_v, wf%n_v, batch_k%length)
!
         call dgemm('T', 'N',                   &
                     (wf%n_v)**2,               &
                     (wf%n_v)*(batch_k%length), &
                     wf%integrals%n_J,          &
                     one,                       &
                     L_J_bd,                    & ! L_bd_J  b is c1-transformed
                     wf%integrals%n_J,          &
                     L_J_ck,                    & ! L_J_ck
                     wf%integrals%n_J,          &
                     zero,                      &
                     g_pqrs,                    & ! (b'd|ck)
                     (wf%n_v)**2)
!
         call mem%dealloc(L_J_ck, wf%integrals%n_J, wf%n_v, batch_k%length)
!
!        :: Term 2: g_bdc'k += sum_J L_J_bd L_J_c'k_c1 ::
!
         call mem%alloc(L_J_ck_c1, wf%integrals%n_J, wf%n_v, batch_k%length)
         call wf%integrals%construct_cholesky_ai_a_c1(L_J_ck_c1, c_ai, 1, wf%n_v, batch_k%first, batch_k%last)
!
         call wf%integrals%get_cholesky_t1(L_J_bd, wf%n_o + 1, wf%n_o + wf%n_v, wf%n_o + 1, wf%n_o + wf%n_v)
!
         call dgemm('T', 'N',                   &
                     (wf%n_v)**2,               &
                     (wf%n_v)*(batch_k%length), &
                     wf%integrals%n_J,          &
                     one,                       &
                     L_J_bd,                    & ! L_bd_J
                     wf%integrals%n_J,          &
                     L_J_ck_c1,                 & ! L_J_ck  c is c1-transformed
                     wf%integrals%n_J,          &
                     one,                       &
                     g_pqrs,                    & ! (bd|c'k)
                     (wf%n_v)**2)
!
!        :: Term 3: g_bdck' += sum_J L_bd_J L_ck_J_c1 ::
!
         call wf%integrals%construct_cholesky_ai_i_c1(L_J_ck_c1, c_ai, 1, wf%n_v, batch_k%first, batch_k%last)
!
         call dgemm('T', 'N',                   &
                     (wf%n_v)**2,               &
                     (wf%n_v)*(batch_k%length), &
                     wf%integrals%n_J,          &
                     one,                       &
                     L_J_bd,                    & ! L_bd_J
                     wf%integrals%n_J,          &
                     L_J_ck_c1,                 & ! L_J_ck  k is c1-transformed
                     wf%integrals%n_J,          &
                     one,                       &
                     g_pqrs,                    & ! (bd|ck')
                     (wf%n_v)**2)
!
         call mem%dealloc(L_J_ck_c1, wf%integrals%n_J, wf%n_v, batch_k%length)
!
!        Sort from g_pqrs = (b'd|ck) + (bd|c'k) + (bd|ck') to h_pqrs ordered as dbck
!
         call mem%alloc(h_pqrs, wf%n_v, wf%n_v, wf%n_v, batch_k%length) ! order dbck
!
         call sort_1234_to_2134(g_pqrs, h_pqrs, wf%n_v, wf%n_v, wf%n_v, batch_k%length)
!
         call mem%dealloc(g_pqrs, wf%n_v, wf%n_v, wf%n_v, batch_k%length)
!
!        Write to file
!
         call wf%g_bdck_c%write_interval(h_pqrs, batch_k)
!
         call mem%dealloc(h_pqrs, wf%n_v, wf%n_v, wf%n_v, batch_k%length)
!
      enddo ! batch_k
!
      call mem%dealloc(L_J_bd, wf%integrals%n_J, wf%n_v, wf%n_v)
!
      call wf%g_bdck_c%close_()
!
!
!     g'_ljck = (lj'|ck) + (lj|ck') + (lj|c'k) ordered as lc,jk
!
!
      req_0 = max((wf%integrals%n_J)*(wf%n_v)**2, &
                  2*(wf%integrals%n_J)*(wf%n_o)**2 + (wf%integrals%n_J)*(wf%n_v)*(wf%n_o))
      req_k = max(2*(wf%n_v)*(wf%n_o)**2, (wf%n_v)*(wf%n_o)**2 + (wf%integrals%n_J)*(wf%n_v))
!
      call mem%alloc(L_J_lj, wf%integrals%n_J, wf%n_o, wf%n_o)
!
      call mem%batch_setup(batch_k,req_0,req_k)
!
      wf%g_ljck_c = direct_stream_file('g_ljck_c', wf%n_v*wf%n_o)
      call wf%g_ljck_c%open_('write')
!
      do k_batch = 1, batch_k%num_batches
!
         call batch_k%determine_limits(k_batch)
!
!        :: Term 1: g_lj'ck = sum_J L_J_jl_c1 L_J_ck ::
!
         call wf%integrals%construct_cholesky_ij_c1(L_J_lj, c_ai, 1, wf%n_o, 1, wf%n_o)
!
         call mem%alloc(L_J_ck, wf%integrals%n_J, wf%n_v, batch_k%length)
         call wf%integrals%get_cholesky_t1(L_J_ck, wf%n_o + 1, wf%n_o + wf%n_v, batch_k%first, batch_k%last)
!
         call mem%alloc(g_pqrs, wf%n_o, wf%n_o, wf%n_v, batch_k%length)
!
         call dgemm('T', 'N',                   &
                     (wf%n_o)**2,               &
                     (wf%n_v)*(batch_k%length), &
                     wf%integrals%n_J,          &
                     one,                       &
                     L_J_lj,                    & ! L_lj_J  j is c1-transformed
                     wf%integrals%n_J,          &
                     L_J_ck,                    & ! L_J_ck
                     wf%integrals%n_J,          &
                     zero,                      &
                     g_pqrs,                    & ! (lj'|ck)
                     (wf%n_o)**2)
!
         call mem%dealloc(L_J_ck, wf%integrals%n_J, wf%n_v, batch_k%length)
!
!        :: Term 2: g_ljck' = sum_J L_J_lj L_J_ck_c1 ::
!
         call mem%alloc(L_J_ck_c1, wf%integrals%n_J, wf%n_v, batch_k%length)
         call wf%integrals%construct_cholesky_ai_i_c1(L_J_ck_c1, c_ai, 1, wf%n_v, batch_k%first, batch_k%last)
!
         call wf%integrals%get_cholesky_t1(L_J_lj, 1, wf%n_o, 1, wf%n_o)
!
         call dgemm('T', 'N',                   &
                     (wf%n_o)**2,               &
                     (wf%n_v)*(batch_k%length), &
                     wf%integrals%n_J,          &
                     one,                       &
                     L_J_lj,                    & ! L_lj_J
                     wf%integrals%n_J,          &
                     L_J_ck_c1,                 & ! L_J_ck  k is c1_transformed
                     wf%integrals%n_J,          &
                     one,                       &
                     g_pqrs,                    & ! (lj|ck')
                     (wf%n_o)**2)
!
!        :: Term 3: g_bdck' = sum_J L_bd_J L_ck_J_c1 ::
!
         call wf%integrals%construct_cholesky_ai_a_c1(L_J_ck_c1, c_ai, 1, wf%n_v, batch_k%first, batch_k%last)
!
         call dgemm('T', 'N',                   &
                     (wf%n_o)**2,               &
                     (wf%n_v)*(batch_k%length), &
                     wf%integrals%n_J,          &
                     one,                       &
                     L_J_lj,                    & ! L_lj_J
                     wf%integrals%n_J,          &
                     L_J_ck_c1,                 & ! L_J_ck  c is c1_transformed
                     wf%integrals%n_J,          &
                     one,                       &
                     g_pqrs,                    & ! (lj|c'k)
                     (wf%n_o)**2)
!
         call mem%dealloc(L_J_ck_c1, wf%integrals%n_J, wf%n_v, batch_k%length)
!
!        Sort from g_pqrs = (lj'|ck) + (lj|ck') + (lj|c'k) to h_pqrs
!
         call mem%alloc(h_pqrs, wf%n_o, wf%n_v, wf%n_o, batch_k%length) ! lcjk
!
         call sort_1234_to_1324(g_pqrs,h_pqrs, wf%n_o, wf%n_o, wf%n_v, batch_k%length)
!
         call mem%dealloc(g_pqrs, wf%n_o, wf%n_o, wf%n_v, batch_k%length)
!
         call wf%g_ljck_c%write_compound_full_batch(h_pqrs, wf%n_o, batch_k)
!
         call mem%dealloc(h_pqrs, wf%n_o, wf%n_v, wf%n_o, batch_k%length)
!
      enddo ! batch_k
!
      call mem%dealloc(L_J_lj, wf%integrals%n_J, wf%n_o, wf%n_o)
!
      call wf%g_ljck_c%close_()
!
   end subroutine construct_c1_integrals_cc3
!
!
   module subroutine construct_c1_fock_cc3(wf, c_ai, F_ia_c1)
!!
!!    Construct C1-transformed fock matrix ov-block
!!    Written by Alexander C. Paul, Feb 2019
!!
!!    Calculates C1-transformed occupied-virtual elements of the Fock matrix
!!    required for the CC3 jacobian and returns it ordered as n_v, n_o
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
      real(dp), dimension(:,:,:), allocatable :: L_J_ia
      real(dp), dimension(:,:,:), allocatable :: L_J_jk_c1
!
      real(dp), dimension(:,:,:,:), allocatable :: g_iajk
!
      integer :: i, a, j
!
      call zero_array(F_ia_c1, wf%n_v*wf%n_o)
!
!     Construct the integrals from the Cholesky Vectors
!
      call mem%alloc(L_J_ia, wf%integrals%n_J, wf%n_o, wf%n_v)
!
      call wf%integrals%get_cholesky_t1(L_J_ia, 1, wf%n_o, wf%n_o + 1, wf%n_o + wf%n_v)
!
      call mem%alloc(L_J_jk_c1, wf%integrals%n_J, wf%n_o, wf%n_o)
!
      call wf%integrals%construct_cholesky_ij_c1(L_J_jk_c1, c_ai, 1, wf%n_o, 1, wf%n_o)
!
      call mem%alloc(g_iajk, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('T', 'N',             &
                  (wf%n_v)*(wf%n_o),   &
                  (wf%n_o)**2,         &
                  wf%integrals%n_J,    &
                  one,                 &
                  L_J_ia,              & ! L_ia_J
                  wf%integrals%n_J,    &
                  L_J_jk_c1,           & ! L_J_jk'
                  wf%integrals%n_J,    &
                  zero,                &
                  g_iajk,              & ! (ia|jk')
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(L_J_jk_c1, wf%integrals%n_J, wf%n_o, wf%n_o)
      call mem%dealloc(L_J_ia, wf%integrals%n_J, wf%n_o, wf%n_v)
!
!     Add contributions and resort to F_ia_c1(a,i)
!
!$omp parallel do private(a,i,j)
      do i = 1, wf%n_o
         do a = 1, wf%n_v
            do j = 1, wf%n_o
!
               F_ia_c1(a,i) = F_ia_c1(a,i) + two*g_iajk(i,a,j,j) - g_iajk(j,a,i,j)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(g_iajk, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
   end subroutine construct_c1_fock_cc3
!
!
   module subroutine jacobian_cc3_b2_fock_cc3(wf, i, j, k, t_abc, u_abc, rho_abij, F_ov_ck)
!!
!!    Jacobian CC3 contribution c1-transformed fock matrix
!!    Written by Alexander C. Paul and Rolf H. Myhre, Feb 2019
!!
!!    rho_2 =+ P^{ab}_{ij} sum_kc (t^abc_ijk - t^cba_ijk) F_kc
!!
!!    The permutations of i,j,k are necessary 
!!    due to the index restrictions in the batching loops
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
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                      :: F_ov_ck
!
!
!     Construct u_abc = t_abc - t_cba
!
      call construct_123_minus_321(t_abc, u_abc, wf%n_v)
!
!     rho_abij += sum_c (t^abc - t^cba)*F_kc
!
      call dgemv('N',               &
                 wf%n_v**2,         &
                 wf%n_v,            &
                 one,               &
                 u_abc,             &
                 wf%n_v**2,         &
                 F_ov_ck(:,k),      &
                 1,                 &
                 one,               &
                 rho_abij(:,:,i,j), &
                 1)
!
!     rho_abkj += sum_c (t^cba - t^abc)*F_ic
!     if i == k, but this is never true
!
      call dgemv('N',               &
                 wf%n_v**2,         &
                 wf%n_v,            &
                 -one,              &
                 u_abc,             &
                 wf%n_v**2,         &
                 F_ov_ck(:,i),      &
                 1,                 &
                 one,               &
                 rho_abij(:,:,k,j), &
                 1)
!
!
      if (i .ne. j .and. j .ne. k) then
!
!        Construct u_abc = t_acb - t_cab
!
         call construct_132_minus_312(t_abc, u_abc, wf%n_v)
!
!        rho_abik += sum_c (t^acb - t^cab)*F_jc
!        Term only needed for j .ne. k and u_abc is 0 if i .eq. j
!
         call dgemv('N',               &
                    wf%n_v**2,         &
                    wf%n_v,            &
                    one,               &
                    u_abc,             &
                    wf%n_v**2,         &
                    F_ov_ck(:,j),      &
                    1,                 &
                    one,               &
                    rho_abij(:,:,i,k), &
                    1)
!
!        rho_abjk += sum_c (t^cab - t^acb)*F_ic
!        Term only needed for j .ne. k and i .ne. j
!
         call dgemv('N',                  &
                     wf%n_v**2,           &
                     wf%n_v,              &
                     -one,                &
                     u_abc,               &
                     wf%n_v**2,           &
                     F_ov_ck(:,i),        &
                     1,                   &
                     one,                 &
                     rho_abij(:,:,j,k),   &
                     1)
!
!        Construct u_abc = t_bac - t_bca
!
         call construct_213_minus_231(t_abc, u_abc, wf%n_v)
!
!
!        rho_abij += sum_c (t^bac - t^bca)*F_kc
!        Term only needed for i .ne. j and u_abc is 0 if j .eq. k
!
         call dgemv('N',               &
                    wf%n_v**2,         &
                    wf%n_v,            &
                    one,               &
                    u_abc,             &
                    wf%n_v**2,         &
                    F_ov_ck(:,k),      &
                    1,                 &
                    one,               &
                    rho_abij(:,:,j,i), &
                    1)
!
!        rho_abki += sum_c (t^bca - t^bac)*F_jc
!        Term only needed for i .ne. j and u_abc is 0 if j .eq. k
!
         call dgemv('N',               &
                    wf%n_v**2,         &
                    wf%n_v,            &
                    -one,              &
                    u_abc,             &
                    wf%n_v**2,         &
                    F_ov_ck(:,j),      &
                    1,                 &
                    one,               &
                    rho_abij(:,:,k,i), &
                    1)
!
      end if
!
   end subroutine jacobian_cc3_b2_fock_cc3
!
!
end submodule jacobian
