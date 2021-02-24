!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2021 the authors of eT
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
submodule (ccsd_class) omega_ccsd
!
!!
!!    Omega submodule
!!
!!    Routines to construct 
!!
!!    Omega_mu =  < mu | exp(-T) H exp(T) | R >
!!
!!    Transfered to the current eT program from the first version 
!!    of eT by Andreas Skeidsvoll and Sarai D. Folkestad, 2018.
!!
!
   implicit none
!
!
contains
!
!
   module subroutine construct_omega_ccsd(wf, omega)
!!
!!    Construct omega (CCSD)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Directs the construction of the projection vector < mu | exp(-T) H exp(T) | R >
!!    for the current amplitudes of the object wfn
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
!
      real(dp), dimension(wf%n_t1 + wf%n_t2), intent(out) :: omega
!
      real(dp), dimension(:,:,:,:), allocatable :: t_aibj, t_abij
      real(dp), dimension(:,:,:,:), allocatable :: omega_aibj, omega_abij
!
      type(timings), allocatable :: timer 
!
      timer = timings('Construct CCSD Omega', pl='normal')
      call timer%turn_on()
!
      call zero_array(omega, wf%n_t1 + wf%n_t2)
!
!     Construct singles contributions
!
      call wf%ccs%construct_omega(omega(1 : wf%n_t1))
!
      call wf%construct_u_aibj()
!
      call wf%omega_doubles_a1(omega(1 : wf%n_t1), wf%u_aibj)
      call wf%omega_doubles_b1(omega(1 : wf%n_t1), wf%u_aibj)
      call wf%omega_doubles_c1(omega(1 : wf%n_t1), wf%u_aibj)
!
!     Construct doubles contributions
!
      call mem%alloc(omega_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call zero_array(omega_aibj, wf%n_t1**2)
!
      call mem%alloc(t_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_aibj, wf%n_t1)
!
      call wf%omega_ccsd_c2(omega_aibj, t_aibj)
      call wf%omega_ccsd_d2(omega_aibj, t_aibj)
      call wf%omega_ccsd_e2(omega_aibj, t_aibj)
!
      call symmetric_sum(omega_aibj, wf%n_t1)
!
      call mem%alloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_1234_to_1324(t_aibj, t_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(t_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(omega_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_1234_to_1324(omega_aibj, omega_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(omega_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%omega_ccsd_b2(omega_abij, t_abij)
      call mem%dealloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call scale_diagonal(half, omega_abij, wf%n_v, wf%n_o)
!
      call packin(omega(wf%n_t1+1:), omega_abij, wf%n_v, wf%n_o)
      call mem%dealloc(omega_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call wf%omega_ccsd_a2(omega(wf%n_t1+1:), wf%t2, right=.true.)
!
      call timer%turn_off()
!
   end subroutine construct_omega_ccsd
!
!
   module subroutine omega_ccsd_a2_ccsd(wf, omega2, t2, right, diagonal_factor)
!!
!!    Omega A2 term
!!    Written by Rolf H. Myhre, Eirik F. Kjønstad, and Sarai D. Folkestad, Jan 2020
!!      
!!    Computes the most expensive CCSD terms omega^ab_ij = sum_cd g_acbd t^cd_ij
!!
!!    The algorithm exploits the symmetry of T2, Omega and g by computing the 
!!    symmetric (+) and antisymmtric (-) combinations of all three entities, 
!!    resulting in an asymptotic scaling of 1/4 * n_v^4 * n_o^2
!!    In addition, the vvvv integrals must typically be computed from the Cholesky vector,
!!    scaling as 1/2 * n_v^4 * n_J
!!
!!    omega^ab_ij += omega+^ab_ij + omega-^ab_ij
!!    omega^ab_ji += omega+^ab_ij - omega-^ab_ij
!!
!!    t+^cd_ij = (1-1/2*delta_cd)*(1-1/2*delta_ij)*(t^cd_ij + t^cd_ji) (c>=d, i>=j)
!!    t-^cd_ij = t^cd_ij - t^cd_ji (c>d, i>j)
!!
!!    g+_cdab = g_acbd + g_adbc (c>=d, a>=b)
!!    g-_cdab = g_acbd - g_adbc (c>d, a>b)
!!
!!    omega+^ab_ij = (1-1/2*delta_ab)*sum_cd g+_acbd t+^cd_ij (a>=b, i>=j, c>=d)
!!    omega-^ab_ij = sum_cd g-_acbd t-^cd_ij (a>b, i>j, c>d)
!!
!!    Note that the diagonal elements of the minus combinations are zero and not computed
!!
!!    right : logical, if true, g+- will be ordered as 2413, else, they will be ordered as 1324
!!
!!    diagonal_factor: real, optional, biorthonormal factor for excited state calculations
!!                                     if present and right is true:
!!                                        diagonal elements of t+ will be scaled, 
!!                                     elsif present and right is false:
!!                                        diagonal elements of omega+ will be scaled
!!
!!    In general, _p dimensions are given by dim*(dim+1)/2, 
!!          while _m dimensions are given by dim*(dim-1)/2
!!
      use packed_array_utilities_r
!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v*wf%n_o*(wf%n_v*wf%n_o+1)/2), intent(inout) :: omega2
      real(dp), dimension(wf%n_v*wf%n_o*(wf%n_v*wf%n_o+1)/2), intent(in)    :: t2
      logical, intent(in)  :: right
      real(dp), optional, intent(in) :: diagonal_factor
!
!     +- T2
!
      real(dp), dimension(:,:), allocatable :: t_p_ijcd
      real(dp), dimension(:,:), allocatable :: t_m_ijcd
!
!     +- integrals
!
      real(dp), dimension(:,:), allocatable :: g_p_cdab
      real(dp), dimension(:,:), allocatable :: g_m_cdab
!
!     Work and pointers for integrals and omega2 +-
!
      real(dp), dimension(:), allocatable, target :: work
      real(dp), dimension(:), pointer, contiguous :: g_vvvv        => null()
      real(dp), dimension(:), pointer, contiguous :: omega2_p_ijab => null()
      real(dp), dimension(:), pointer, contiguous :: omega2_m_ijab => null()
!
!     Batching and memory handling variables
!
      integer :: req0, req_a, req_b, req_ab, max_req
!
      integer :: current_a_batch
      integer :: current_b_batch    
      integer :: first_p, last_p, first_q, last_q
      integer :: first_r, last_r, first_s, last_s
!
      integer :: ab_p_dim, ab_m_dim, work_dim
      integer :: n_v_p, n_v_m, n_o_p, n_o_m
      integer :: ab_p, ab_m
      logical :: dim_gt_one
!
      type(batching_index) :: batch_a
      type(batching_index) :: batch_b
!
      type(timings) :: ccsd_a2_timer, ccsd_a2_integral_timer
!
      ccsd_a2_timer = timings('omega ccsd a2', pl='verbose')
      ccsd_a2_integral_timer = timings('omega ccsd a2 g_abcd', pl='verbose')
!
      call ccsd_a2_timer%turn_on()
!
!     Compute the symmetric and antisymmetric packed n_v and n_o
!
      n_v_p = wf%n_v*(wf%n_v+1)/2
      n_o_p = wf%n_o*(wf%n_o+1)/2
      n_v_m = wf%n_v*(wf%n_v-1)/2
      n_o_m = wf%n_o*(wf%n_o-1)/2
      dim_gt_one = wf%n_v .gt. 1 .and. wf%n_o .gt. 1
!
!     Set up T2 +-
!
      call mem%alloc(t_p_ijcd, n_o_p, n_v_p)
      call mem%alloc(t_m_ijcd, n_o_m, n_v_m)
!
      call construct_plus_minus_2413_from_packed(t2, t_p_ijcd, t_m_ijcd, wf%n_o, wf%n_v)
!
      if(present(diagonal_factor) .and. right) then
         call scale_double_packed_diagonal(t_p_ijcd, wf%n_o, wf%n_v, diagonal_factor)
      endif
!
      call scale_double_packed_blocks(t_p_ijcd, wf%n_o, wf%n_v, p_factor=half, q_factor=half)
!
!     Set up batching, first see if can get away with a single batch
!
      batch_a = batching_index(wf%n_v)
      batch_b = batching_index(wf%n_v)
!
      max_req = n_v_p*wf%n_v + wf%n_v**2*(wf%n_v**2+1)/2 + wf%n_v*wf%n_o*(wf%n_v*wf%n_o+1)/2
      call wf%eri%get_eri_t1_packed_mem('vv', max_req, 1, wf%n_v, qp=right)
!
      req0 = 0
      req_a = 0
      req_b = 0
      call wf%eri%get_eri_t1_mem('vvvv', req_a, req_b, 1, wf%n_v, 1, wf%n_v, qp=right, sr=right)
      req_ab = wf%n_v**2 + (max(wf%n_v, wf%n_o))**2
!
      call mem%batch_setup(batch_a, batch_b, req0, req_a, req_b, req_ab, req_single_batch=max_req)
!
!     Figure out the batch dependent dimensions and allocate g+- and work
!
      if (batch_a%num_batches .eq. 1) then
         ab_p_dim   = wf%n_v*(wf%n_v+1)/2
         ab_m_dim   = wf%n_v*(wf%n_v-1)/2
      else
         ab_p_dim   = batch_a%max_length**2
         ab_m_dim   = batch_a%max_length**2
      endif
      work_dim = max(ab_p_dim*wf%n_v**2, ab_p_dim*n_o_p + ab_m_dim*n_o_m)
!
      call mem%alloc(g_p_cdab, n_v_p, ab_p_dim)
      call mem%alloc(g_m_cdab, n_v_m, ab_m_dim)
      call mem%alloc(work,work_dim)
!
      g_vvvv(1:wf%n_v**2*ab_p_dim)    => work
      omega2_p_ijab(1:n_o_p*ab_p_dim) => work
      omega2_m_ijab(1:n_o_m*ab_m_dim) => work(n_o_p*ab_p_dim+1:)
!
!     Set batch independent integral dimensions depending on left or right calculation
!
      if(right) then
         first_q = 1
         last_q  = wf%n_v
         first_s = 1
         last_s  = wf%n_v
      else
         first_p = 1
         last_p  = wf%n_v
         first_r = 1
         last_r  = wf%n_v
      endif
!
      do current_a_batch = 1, batch_a%num_batches
!
         call batch_a%determine_limits(current_a_batch)
!
         do current_b_batch = 1, current_a_batch
!
            call batch_b%determine_limits(current_b_batch)
!
!           Set batch dependent integral dimensions depending on left or right calculation
!
            if(right) then
               first_p = batch_a%first
               last_p  = batch_a%last
               first_r = batch_b%first
               last_r  = batch_b%last
            else
               first_q = batch_a%first
               last_q  = batch_a%last
               first_s = batch_b%first
               last_s  = batch_b%last
            endif
!
!           If batches are the same, we can compute the packed integral 
!           and double packed g and omega, else we need the full integral block
!
            if (current_b_batch .eq. current_a_batch) then
!
!              Packed ab dimensions for current batch set
!
               ab_p = batch_a%length*(batch_a%length+1)/2
               ab_m = batch_a%length*(batch_a%length-1)/2
!
               call ccsd_a2_integral_timer%turn_on()
!
               call wf%eri%get_eri_t1_packed('vv', g_vvvv, &
                                             first_p, last_p, first_q, last_q, qp=right)
!
               call ccsd_a2_integral_timer%freeze()
!
               call construct_plus_minus_1324_from_RFP(work, g_p_cdab, g_m_cdab, &
                                                       wf%n_v, batch_a%length)
!
               call dgemm('N','N',       &
                          n_o_p,         &
                          ab_p,          &
                          n_v_p,         &
                          half,          &
                          t_p_ijcd,      & !t+_ij_cd
                          n_o_p,         &
                          g_p_cdab,      & !g+_cd_ab
                          n_v_p,         &
                          zero,          &
                          omega2_p_ijab, & !omega+_ij_ab
                          n_o_p)
!
               if (dim_gt_one .and. ab_m .gt. 0) then
                  call dgemm('N','N',       &
                             n_o_m,         &
                             ab_m,          &
                             n_v_m,         &
                             half,          &
                             t_m_ijcd,      & !t-_ij_cd
                             n_o_m,         &
                             g_m_cdab,      & !g-_cd_ab
                             n_v_m,         &
                             zero,          &
                             omega2_m_ijab, & !omega-_ij_ab
                             n_o_m)
               endif
!
               call scale_double_packed_blocks(omega2_p_ijab, wf%n_o, batch_a%length, q_factor=half)
!
               if(present(diagonal_factor) .and. .not. right) then
                  call scale_double_packed_diagonal(omega2_p_ijab, wf%n_o, batch_a%length, &
                                                    diagonal_factor)
               endif
!
               call add_double_packed_plus_minus(omega2, omega2_p_ijab, omega2_m_ijab, &
                                                 wf%n_o, batch_a%length, wf%n_v, batch_a%first-1)
!
            else
!
!              Unpacked ab dimensions for current batch set
!
               ab_p = batch_a%length*batch_b%length
               ab_m = batch_a%length*batch_b%length
!
               call ccsd_a2_integral_timer%turn_on()
!
               call wf%eri%get_eri_t1('vvvv', g_vvvv,                   &
                                      first_p, last_p, first_q, last_q, &
                                      first_r, last_r, first_s, last_s, &
                                      qp=right, sr=right)
!
               call ccsd_a2_integral_timer%freeze()
!
               call construct_plus_minus_1324_from_full(g_vvvv, g_p_cdab, g_m_cdab, &
                                                        wf%n_v, batch_a%length, batch_b%length)
!
               call dgemm('N','N',        &
                           n_o_p,         &
                           ab_p,          &
                           n_v_p,         &
                           half,          &
                           t_p_ijcd,      & !t+_ij_cd
                           n_o_p,         &
                           g_p_cdab,      & !g+_cd_ab
                           n_v_p,         &
                           zero,          &
                           omega2_p_ijab, & !omega+_ij_ab
                           n_o_p)
!
               if (dim_gt_one) then
                  call dgemm('N','N',        &
                              n_o_m,         &
                              ab_m,          &
                              n_v_m,         &
                              half,          &
                              t_m_ijcd,      & !t-_ij_cd
                              n_o_m,         &
                              g_m_cdab,      & !g-_cd_ab
                              n_v_m,         &
                              zero,          &
                              omega2_m_ijab, & !omega-_ij_ab
                              n_o_m)
               endif
!
               call add_single_packed_plus_minus(omega2, omega2_p_ijab, omega2_m_ijab,   &
                                                 wf%n_o, batch_a%length, batch_b%length, &
                                                 wf%n_v, batch_a%first-1, batch_b%first-1)
!
            endif
!
         enddo ! End batching over b
      enddo ! End batching over a
!
      g_vvvv        => null()
      omega2_p_ijab => null()
      omega2_m_ijab => null()
!
      call mem%dealloc(g_p_cdab, n_v_p, ab_p_dim)
      call mem%dealloc(g_m_cdab, n_v_m, ab_m_dim)
      call mem%dealloc(work,work_dim)
!
      call mem%dealloc(t_p_ijcd, n_o_p, n_v_p)
      call mem%dealloc(t_m_ijcd, n_o_m, n_v_m)
!
      call ccsd_a2_timer%turn_off()
      call ccsd_a2_integral_timer%turn_off()
!
   end subroutine omega_ccsd_a2_ccsd
!
!
   module subroutine omega_ccsd_b2_ccsd(wf, omega_abij, t_abij)
!!
!!    Omega B2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Omega B2 = g_aibj + sum_(kl) t_akbl*(g_kilj + sum_(cd) t_cidj * g_kcld)
!!
!!    Structure: g_kilj is constructed first and reordered as g_klij.
!!    Then the contraction over cd is performed, and the results added to g_klij.
!!    t_ak_bl is then reordered as t_ab_kl and the contraction over kl is performed.
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_o,wf%n_o), intent(inout):: omega_abij
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_o,wf%n_o), intent(in):: t_abij
!
!     Integrals
!
      real(dp), dimension(:,:,:,:), allocatable :: g_aibj
      real(dp), dimension(:,:,:,:), allocatable :: g_kcld
      real(dp), dimension(:,:,:,:), allocatable :: g_klcd
      real(dp), dimension(:,:,:,:), allocatable :: g_klij
      real(dp), dimension(:,:,:,:), allocatable :: g_kilj
!
      type(timings) :: ccsd_b2_timer
!
      ccsd_b2_timer = timings('omega ccsd b2', pl='verbose')
      call ccsd_b2_timer%turn_on()
!
!     Construct g_aibj and add to omega2 
!
      call mem%alloc(g_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call wf%eri%get_eri_t1('vovo', g_aibj)
      call add_1324_to_1234(one, g_aibj, omega_abij,  wf%n_v,  wf%n_v,  wf%n_o,  wf%n_o)
      call mem%dealloc(g_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Allocate and construct g_kilj
!
      call mem%alloc(g_kilj, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call wf%eri%get_eri_t1('oooo', g_kilj)
!
      call mem%alloc(g_klij,  wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call sort_1234_to_1324(g_kilj, g_klij, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call mem%dealloc(g_kilj, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
!     Allocate and construct g_kcld
!
      call mem%alloc(g_kcld, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%eri%get_eri_t1('ovov', g_kcld)
!
!     Reorder g_kcld as g_klcd
!
      call mem%alloc(g_klcd,  wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
      call sort_1234_to_1324(g_kcld, g_klcd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%dealloc(g_kcld, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call dgemm('N','N',      &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  (wf%n_v)**2, &
                  one,         &
                  g_klcd,      &
                  (wf%n_o)**2, &
                  t_abij,      &
                  (wf%n_v)**2, &
                  one,         &
                  g_klij,      &
                  (wf%n_o)**2)
!
      call mem%dealloc(g_klcd,  wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
!     omega_abij = sum_(kl) t_ab_kl*X_kl_ij
!
      call dgemm('N','N',      &
                  (wf%n_v)**2, &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  one,         &
                  t_abij,      & ! t_ab_kl
                  (wf%n_v)**2, &
                  g_klij,      &
                  (wf%n_o)**2, &
                  one,         &
                  omega_abij,  &
                  (wf%n_v)**2)
!
      call mem%dealloc(g_klij,  wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call ccsd_b2_timer%turn_off()
!
   end subroutine omega_ccsd_b2_ccsd
!
!
   module subroutine omega_ccsd_c2_ccsd(wf, omega_aibj, t_aibj)
!!
!!    Omega C2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Omega C2 = -1/2 * sum_(ck) t_bk_cj*(g_kiac -1/2 sum_(dl)t_al_di * g_kdlc)
!!                    - sum_(ck) t_bk_ci*(g_kj_ac - sum_(dl)t_al_dj * g_kdlc)
!!                    - 1/2 * sum_ck u_jk^bc g_acki
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(inout):: omega_aibj
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(in):: t_aibj
!
!     Integrals
!
      real(dp), dimension(:,:,:,:), allocatable :: g_kdlc
      real(dp), dimension(:,:,:,:), allocatable :: g_dlck
      real(dp), dimension(:,:,:,:), allocatable :: g_kiac
      real(dp), dimension(:,:,:,:), allocatable :: g_aick
!
!     Reordered T2 amplitudes
!
      real(dp), dimension(:,:,:,:), allocatable :: t_aidl
      real(dp), dimension(:,:,:,:), allocatable :: t_ckbj
!
!    Intermediates for matrix multiplication
!
      real(dp), dimension(:,:,:,:), allocatable :: X_aick
      real(dp), dimension(:,:,:,:), allocatable :: Y_aibj
!
!     Reordered U2 amplitudes
!
      real(dp), dimension(:,:,:,:), allocatable :: omega_a_batch
!
!     Indices
!
      integer :: a, b, c
      integer :: i, j, k
!
!     Batching and memory handling
!
      integer :: req0, req1
      integer :: current_a_batch = 0
!
      type(batching_index) :: batch_a
!
      type(timings) :: ccsd_c2_timer
!
      ccsd_c2_timer = timings('omega ccsd c2', pl='verbose')
      call ccsd_c2_timer%turn_on()
!
!     Sort t_al_di = t_li^ad as t_aidl (1234 to 1432)
!
      call mem%alloc(t_aidl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call sort_1234_to_1432(t_aibj, t_aidl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Get g_kdlc
!
      call mem%alloc(g_kdlc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%eri%get_eri_t1('ovov', g_kdlc)
!
!     Sort g_kdlc to g_dlck (1234 to 2341)
!
      call mem%alloc(g_dlck,  wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call sort_1234_to_2341(g_kdlc, g_dlck, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%dealloc(g_kdlc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     -1/2*sum_(dl) t_aidl*g_dlck = X_aick
!
      call mem%alloc(X_aick,  wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  -half,             &
                  t_aidl,            &
                  (wf%n_o)*(wf%n_v), &
                  g_dlck,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_aick,            &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(g_dlck,  wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(t_aidl,  wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Allocate a holder for - 1/2 * sum_ck u_jk^bc g_acki,
!     constructed in batches over the a index below
!
!     Constructing g_kiac
!
!     Prepare for batching
!
      req0 = wf%n_o**2*wf%eri%n_J
!
      req1 = wf%n_v*wf%eri%n_J + (wf%n_o)*(wf%n_v**2)
!
      batch_a = batching_index(wf%n_v)
!
      call mem%batch_setup(batch_a, req0, req1)
!
!     Loop over the number of a batches
!
      do current_a_batch = 1, batch_a%num_batches
!
         call batch_a%determine_limits(current_a_batch)
!
!        Allocate and construct g_kiac
!
         call mem%alloc(g_kiac, wf%n_o, wf%n_o, batch_a%length, wf%n_v)
!
         call wf%eri%get_eri_t1('oovv', g_kiac, first_r=batch_a%first, last_r=batch_a%last)
!
!        X_aick = X_aick + g_kiac
!
!$omp parallel do private(a, i, k, c)
         do k = 1, wf%n_o
            do c = 1, wf%n_v
               do i = 1, wf%n_o
                  do a = 1, batch_a%length
!
                     X_aick(a+batch_a%first-1, i, c, k) = X_aick(a+batch_a%first-1, i, c, k) &
                                                        + g_kiac(k, i, a, c)
!
                  enddo
               enddo
            enddo
         enddo
!$omp end parallel do
!
!        Calculate the contribution to the term
!
!        omega_aibj = - 1/2 * sum_ck u_jk^bc g_acki
!
         call mem%alloc(g_aick, batch_a%length, wf%n_o, wf%n_v, wf%n_o)
!
         call sort_1234_to_3241(g_kiac, g_aick, wf%n_o, wf%n_o, batch_a%length, wf%n_v)
!
         call mem%dealloc(g_kiac, wf%n_o, wf%n_o, batch_a%length, wf%n_v)
!
!        - 1/2 * sum_ck u_jk^bc g_acki = -1/2 * sum_ck g_aick u_ckbj
!
         call mem%alloc(omega_a_batch, batch_a%length, wf%n_o, wf%n_v, wf%n_o)
!
         call dgemm('N','N',                 &
                  (wf%n_o)*(batch_a%length), &
                  (wf%n_o)*(wf%n_v),         &
                  (wf%n_o)*(wf%n_v),         &
                  -one/two,                  &
                  g_aick,                    &
                  (wf%n_o)*(batch_a%length), &
                  wf%u_aibj,                 & ! u_ck_bj
                  (wf%n_o)*(wf%n_v),         &
                  zero,                      &
                  omega_a_batch,             &
                  (wf%n_o)*batch_a%length)
!
         call mem%dealloc(g_aick, batch_a%length, wf%n_o, wf%n_v, wf%n_o)
!
!$omp parallel do private(a, i, b, j)
         do j = 1, wf%n_o
            do b = 1, wf%n_v
               do i = 1, wf%n_o
                  do a = 1, batch_a%length
!
                     omega_aibj(a + batch_a%first - 1, i, b, j) = omega_aibj(a + batch_a%first - 1, i, b, j)&
                                                                + omega_a_batch(a, i, b, j)               
!
                  enddo
               enddo
            enddo
         enddo
!$omp end parallel do
!
         call mem%dealloc(omega_a_batch, batch_a%length, wf%n_o, wf%n_v, wf%n_o)
!
      enddo ! End of batching
!
!     Add the - 1/2 * sum_ck u_jk^bc g_acki term to omega
!
!
!     Reorder t_bkcj = t_kj^bc as t_ckbj
!
      call mem%alloc(t_ckbj,  wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_3214(t_aibj, t_ckbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Form intermediate Y_aibj = - sum_(ck) X_aick*t_ckbj
!
      call mem%alloc(Y_aibj,  wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  -one,              &
                  X_aick,            &
                  (wf%n_o)*(wf%n_v), &
                  t_ckbj,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  Y_aibj,            &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(X_aick,  wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(t_ckbj,  wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Omega_aibj,1 = P_ai_bj ( 1/2*Y_aibj + Y_aj_bi )
!
!$omp parallel do private(i, a, j, b)
         do j = 1, wf%n_o
            do b = 1, wf%n_v
               do i = 1, wf%n_o
                  do a = 1, wf%n_v
!
                     omega_aibj(a, i, b, j) = omega_aibj(a, i, b, j) + half*Y_aibj(a, i, b, j) + Y_aibj(a, j, b, i) 
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(Y_aibj,  wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call ccsd_c2_timer%turn_off()
!
   end subroutine omega_ccsd_c2_ccsd
!
!
   module subroutine omega_ccsd_d2_ccsd(wf, omega_aibj, t_aibj)
!!
!!    Omega D2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Calculates the D2 term,
!!
!!      D2: sum_ck u_jk^bc g_aikc
!!        + 1/4 * sum_ck u_jk^bc sum_dl L_ldkc u_il^ad,
!!
!!    where
!!
!!        u_jk^bc = 2 * t_jk^bc - t_kj^bc,
!!        L_ldkc  = 2 * g_ldkc  - g_lckd.
!!
!!    The first and second terms are referred to as D2.1 and D2.2.
!!
!!    The routine adds the terms in the following order: D2.2, D2.1
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(inout):: omega_aibj
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(in):: t_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ldkc         ! g_ldkc
      real(dp), dimension(:,:,:,:), allocatable :: L_ldkc         ! L_ldkc = 2 * g_ldkc - g_lckd
      real(dp), dimension(:,:,:,:), allocatable :: u_aild         ! u_il^ad = 2 * t_il^ad - t_li^ad
      real(dp), dimension(:,:,:,:), allocatable :: Z_aikc         ! An intermediate, see below
      real(dp), dimension(:,:,:), allocatable   :: L_J_kc, L_J_ai ! Cholesky vectors, term 1
      real(dp), dimension(:,:,:), allocatable   :: X_bj_J         ! Intermediate, term 1
!
      type(timings) :: ccsd_d2_timer
!
      ccsd_d2_timer = timings('omega ccsd d2', pl='verbose')
      call ccsd_d2_timer%turn_on()
!
!     :: Calculate the D2.2 term of omega ::
!
!     Form L_ld_kc = L_ldkc = 2*g_ldkc(ld,kc) - g_ldkc(lc,kd)
!
      call mem%alloc(g_ldkc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call wf%eri%get_eri_t1('ovov', g_ldkc)
!
      call mem%alloc(L_ldkc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call copy_and_scale(two, g_ldkc, L_ldkc, (wf%n_o)**2*(wf%n_v)**2)
      call add_1432_to_1234(-one, g_ldkc, L_ldkc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%dealloc(g_ldkc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%alloc(u_aild, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      
      call zero_array(u_aild, wf%n_t1**2)
      call add_1243_to_1234(two, t_aibj, u_aild, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call add_1342_to_1234(-one, t_aibj, u_aild, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
!     Form the intermediate Z_aikc = sum_dl u_aild L_ldkc and set it to zero
!
      call mem%alloc(Z_aikc,  wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  u_aild,            &
                  (wf%n_o)*(wf%n_v), &
                  L_ldkc,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  Z_aikc,            &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(L_ldkc,  wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     Form the D2.2 term, 1/4 sum_kc Z_aikc u_kc_bj = 1/4 sum_kc Z_aikc(ai,kc) u_aild(bj,kc)
!
      call dgemm('N','T',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one/four,          &
                  Z_aikc,            &
                  (wf%n_o)*(wf%n_v), &
                  u_aild,            &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  omega_aibj,        &
                  (wf%n_o)*(wf%n_v))
!
!     Some justification for the above matrix multiplication. We have
!
!           1/4 * sum_ck (sum_dl u_il^ad L_ldkc) u_jk^bc = 1/4 * sum_ck Z_ai,kc u_kc,bj,
!
!     where Z_ai,kc = sum_dl u_ai,ld L_ld,kc. Note that u_aild(ai,ld) = u_il^ad,
!     which means that u_aild(bj,kc)^T = u_aild(kc,bj) = u_kj^cb = u_jk^bc.
!
      call mem%dealloc(Z_aikc, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
!     :: Calculate the D2.1 term of omega ::
!
!     u_jk^bc g_aikc = L_ai^J (u_bjkc L_kc^J) = L_ai^J X_bj^J
!
!     X_bj_J = u_bjkc L_J_kc
!
      call mem%alloc(L_J_kc, wf%eri%n_J, wf%n_o, wf%n_v)
      call wf%eri%get_cholesky_t1(L_J_kc, 1, wf%n_o, wf%n_o + 1, wf%n_mo)
!
      call mem%alloc(X_bj_J, wf%n_v, wf%n_o, wf%eri%n_J)
!
      call dgemm('N', 'T',           &
                  (wf%n_o)*(wf%n_v), &
                  wf%eri%n_J,        &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  u_aild,            & ! u_bj,kc
                  (wf%n_o)*(wf%n_v), &
                  L_J_kc,            &
                  wf%eri%n_J,        &
                  zero,              &
                  X_bj_J,            &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(L_J_kc, wf%eri%n_J, wf%n_o, wf%n_v)
!
      call mem%alloc(L_J_ai, wf%eri%n_J, wf%n_v, wf%n_o)
      call wf%eri%get_cholesky_t1(L_J_ai, wf%n_o + 1, wf%n_mo, 1, wf%n_o)
!
!     omega_aibj =+ L_J_ai X_bj_J
!
      call dgemm('T', 'T',           &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  wf%eri%n_J,        &
                  one,               &
                  L_J_ai,            &
                  wf%eri%n_J,        &
                  X_bj_J,            &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  omega_aibj,        &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(X_bj_J, wf%n_v, wf%n_o, wf%eri%n_J)
      call mem%dealloc(L_J_ai, wf%eri%n_J, wf%n_v, wf%n_o)
      call mem%dealloc(u_aild, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call ccsd_d2_timer%turn_off()
!
   end subroutine omega_ccsd_d2_ccsd
!
!
   module subroutine omega_ccsd_e2_ccsd(wf, omega_aibj, t_aibj)
!!
!!    Omega E2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Calculates the E2 term,
!!
!!      E2: sum_c t_ij^ac (F_bc - sum_dkl g_ldkc u_kl^bd)
!!        - sum_k t_ik^ab (F_kj + sum_cdl g_ldkc u_lj^dc),
!!
!!    where
!!
!!        u_kl^bc = 2 * t_kl^bc - t_lk^bc.
!!
!!    The first term is referred to as the E2.1 term, and comes out ordered as (b,jai).
!!    The second term is referred to as the E2.2 term, and comes out ordered as (aib,j).
!!
!!    Both are permuted added to the projection vector element omega2(ai,bj) of
!!    the wavefunction object wf.
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(inout):: omega_aibj
      real(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(in):: t_aibj
!
!     Vectors for E2.1 term
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ldkc      ! g_ldkc
      real(dp), dimension(:,:,:,:), allocatable :: u_bldk      ! u_kl^bd
      real(dp), dimension(:,:), allocatable :: X_b_c           ! An intermediate, see below for definition
!
!     Vectors for E2.2 term
!
      real(dp), dimension(:,:), allocatable :: Y_k_j        ! An intermediate, see below for definition
!
      type(timings) :: ccsd_e2_timer 
!
      ccsd_e2_timer = timings('omega ccsd e2', pl='verbose')
      call ccsd_e2_timer%turn_on()
!
!     :: Calculate the E2.1 term of omega ::
!
!     Form u_bldk = u_kl^bd = u_bkdl
!                  = 2 * t_kl^bd - t_lk^bd = 2 * t_bkdl(bk,dl) - t_bkdl(bl, dk)
!
      call mem%alloc(u_bldk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      
      call copy_and_scale(-one, t_aibj, u_bldk, (wf%n_v)**2*(wf%n_o)**2)
      call add_1432_to_1234(two, t_aibj, u_bldk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Form g_ldkc = g_ldkc
!
      call mem%alloc(g_ldkc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call wf%eri%get_eri_t1('ovov', g_ldkc)
!
!     Make the intermediate X_b_c = F_bc - sum_dkl g_ldkc u_kl^bd and set to zero
!
      call mem%alloc(X_b_c, wf%n_v, wf%n_v)
!
!     First, copy the virtual-virtual Fock matrix into the intermediate
!
      call dcopy((wf%n_v)**2, wf%fock_ab, 1, X_b_c, 1) ! X_b_c = F_bc
!
!     Then, add the second contribution,
!     - sum_dkl g_ldkc u_kl^bd = - sum_dkl u_b_kdl * g_kdl_c, to X_b_c
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  -one,                 &
                  u_bldk,               & ! u_b_ldk
                  wf%n_v,               &
                  g_ldkc,               & ! g_ldk_c
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  X_b_c,                &
                  wf%n_v)
!
!     Form t_cjai = t_ij^ac = t_ji^ca
!
!     Form the E2.1 term
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_v,               &
                  one,                  &
                  X_b_c,                &
                  wf%n_v,               &
                  t_aibj,               & ! t_c_jai
                  wf%n_v,               &
                  one,                  &
                  omega_aibj,           & ! omega_bjai But we will symmetrize later
                  wf%n_v)
!
!     Add the E2.1 term to the omega vector and deallocations
!
      call mem%dealloc(X_b_c, wf%n_v, wf%n_v)
!
!     :: Calculate E2.2 term of omega ::
!
!     Make the intermediate Y_k_j = F_kj  + sum_cdl u_lj^dc g_ldkc
!                                 = F_kj  + sum_cdl u_lj^dc g_kcld
!                                 = F_k_j + sum_cdl g_k_cld * u_cld_j
!
      call mem%alloc(Y_k_j, wf%n_o, wf%n_o)
!
!     Copy the occupied-occupied Fock matrix, such that Y_k_j = F_kj
!
      call dcopy((wf%n_o)**2, wf%fock_ij, 1, Y_k_j, 1)
!
!     Add sum_cdl g_k_dlc u_dlc_j to Y_k_j, such that
!     Y_k_j = F_k_j + sum_cdl g_k_cld u_cld_j
!
!     Note:
!
!     g_ldkc(kc,ld) = g_kcld              -> pretend that this is g_k_cld
!     u_b_ldk(c,ldj) = u_jl^cd (= u_lj^dc) -> pretend that this is u_cld_j
!
      call dgemm('N','N',              &
                 wf%n_o,               &
                 wf%n_o,               &
                 (wf%n_o)*(wf%n_v)**2, &
                 one,                  &
                 g_ldkc,               & ! g_k_cld
                 wf%n_o,               &
                 u_bldk,               & ! u_cld_j
                 (wf%n_o)*(wf%n_v)**2, &
                 one,                  &
                 Y_k_j,                &
                 wf%n_o)
!
      call mem%dealloc(u_bldk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(g_ldkc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     Calculate the E2.2 term,
!     - sum_k t_aib_k Y_k_j = - sum_k t_ik^ab (F_kj + sum_cdl g_ldkc u_lj^dc)
!
!     Note: t_cjai = t_ji^ca => t_cjai(ai,bk) = t_ik^ab;
!     thus, we can treat t_cjai as t_aib_k = t_ik^ab.
!
      call dgemm('N','N',              &
                 (wf%n_o)*(wf%n_v)**2, &
                 wf%n_o,               &
                 wf%n_o,               &
                 -one,                 &
                 t_aibj,               & ! t_aib_k
                 (wf%n_o)*(wf%n_v)**2, &
                 Y_k_j,                &
                 wf%n_o,               &
                 one,                  &
                 omega_aibj,           & ! omega2_aib_j
                 (wf%n_o)*(wf%n_v)**2)
!
!     Deallocate Y_k_j and the amplitudes
!
      call mem%dealloc(Y_k_j, wf%n_o, wf%n_o)
!
      call ccsd_e2_timer%turn_off()
!
   end subroutine omega_ccsd_e2_ccsd
!
!
   module subroutine construct_u_aibj_ccsd(wf)
!!
!!    Construct u_aibj
!!    Written by Tor S. Haugland, Nov 2019
!!
!!    Constructs 
!!
!!       u_aibj = 2t_aibj - t_ajbi
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
!
      real(dp), dimension(:,:,:,:), allocatable :: t_aibj
!
      type(timings) :: timer
!
      timer = timings('Construct u_aibj', pl='verbose')
      call timer%turn_on()
!
      call mem%alloc(t_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_aibj, wf%n_t1)
!
      call copy_and_scale(two,    t_aibj, wf%u_aibj, wf%n_v**2 * wf%n_o**2)
      call add_1432_to_1234(-one, t_aibj, wf%u_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(t_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine construct_u_aibj_ccsd
!
!
end submodule omega_ccsd 
