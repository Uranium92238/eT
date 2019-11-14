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
submodule (ccsd_class) omega_ccsd_complex
!
!!
!!    Omega submodule (CCSD)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Routines to construct 
!!
!!    Ω =  < mu | exp(-T) H exp(T) | R >
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
   module subroutine construct_omega_ccsd_complex(wf, omega)
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
      complex(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: omega
!
      complex(dp), dimension(:,:), allocatable :: omega1
      complex(dp), dimension(:,:,:,:), allocatable :: t_aibj, t_abij
      complex(dp), dimension(:,:,:,:), allocatable :: omega_aibj, omega_abij
!
!     Construct singles contributions
!
      call mem%alloc(omega1, wf%n_v, wf%n_o)
      call zero_array_complex(omega1, wf%n_t1)
!
      call wf%omega_ccs_a1_complex(omega1)
!
      call wf%construct_u_aibj_complex()
!
      call wf%omega_doubles_a1_complex(omega1, wf%u_aibj_complex)
      call wf%omega_doubles_b1_complex(omega1, wf%u_aibj_complex)
      call wf%omega_doubles_c1_complex(omega1, wf%u_aibj_complex)
!
      call zcopy(wf%n_t1, omega1, 1, omega, 1)
!
      call mem%dealloc(omega1, wf%n_v, wf%n_o)
!
!     Construct doubles contributions
!
      call mem%alloc(omega_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call zero_array_complex(omega_aibj, wf%n_t1**2)
!
      call mem%alloc(t_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2_complex, t_aibj, wf%n_t1)
!
      call wf%omega_ccsd_c2_complex(omega_aibj, t_aibj)
      call wf%omega_ccsd_d2_complex(omega_aibj, t_aibj)
      call wf%omega_ccsd_e2_complex(omega_aibj, t_aibj)
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
      call wf%omega_ccsd_a2_complex(omega_abij, t_abij)
      call wf%omega_ccsd_b2_complex(omega_abij, t_abij)
      call mem%dealloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call scale_diagonal(half_complex, omega_abij, wf%n_v, wf%n_o)
!
      call packin(omega(wf%n_t1+1 : wf%n_gs_amplitudes), omega_abij, wf%n_v, wf%n_o)
      call mem%dealloc(omega_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
   end subroutine construct_omega_ccsd_complex
!
!
   module subroutine omega_ccsd_a2_ccsd_complex(wf, omega_abij, t_abij)
!!
!!    Omega A2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!      
!!    A2 = sum_(cd)g_acbd * t_cidj 
!!
!!    Structure: Batching over both a and b 
!!                t^+_ci_dj = t_cidj + t_di_cj
!!                t^-_ci_dj = t_cidj - t_di_cj
!!                g^+_ac_bd = g_acbd + g_bc_ad
!!                g^-_ac_bd = g_acbd - g_bc_ad
!!
!!                omega_A2_ai_bj = 1/4*(g^+_ac_bd*t^+_ci_dj + g^-_ac_bd*t^-_ci_dj) = omega_A2_bj_ai
!!                omega_A2_aj_bi = 1/4*(g^+_ac_bd*t^+_ci_dj - g^-_ac_bd*t^-_ci_dj) = omega_A2_bi_aj
!!
      implicit none
!
      class(ccsd) :: wf
!
      complex(dp), dimension(wf%n_v,wf%n_v,wf%n_o,wf%n_o), intent(inout) :: omega_abij
      complex(dp), dimension(wf%n_v,wf%n_v,wf%n_o,wf%n_o), intent(in):: t_abij
!
!     Integrals
!
      complex(dp), dimension(:,:,:,:), allocatable :: g_acbd
      complex(dp), dimension(:,:), allocatable :: g_p_abcd
      complex(dp), dimension(:,:), allocatable :: g_m_abcd
!
!     Reordered T2 amplitudes
!
      complex(dp), dimension(:,:), allocatable :: t_p_cdij
      complex(dp), dimension(:,:), allocatable :: t_m_cdij
!
!     Reordered omega 2
!
      complex(dp), dimension(:,:), allocatable :: omega2_p_abij
      complex(dp), dimension(:,:), allocatable :: omega2_m_abij
!
!     Indices
!
      integer :: a, b, c, d, a_full, b_full
      integer :: i, j
!
      integer :: ab, cd
      integer :: ai, aj, bj, bi, ci, cj, dj, di
      integer :: ij
!
      integer :: aibj, biaj, cidj, dicj 
!
      complex(dp) :: diag_factor
!
!     Batching and memory handling variables
!
      integer :: req0, req1_a, req1_b, rec2
!
      integer :: current_a_batch
      integer :: current_b_batch    
!
      integer :: n_v_packed, n_o_packed, batch_a_packed
!
      type(batching_index) :: batch_a
      type(batching_index) :: batch_b
!
      type(timings) :: ccsd_a2_timer, ccsd_a2_integral_timer
!
      ccsd_a2_timer = timings('omega ccsd a2')
      ccsd_a2_integral_timer = timings('omega ccsd a2 g_abcd')
!      
      call ccsd_a2_timer%turn_on()
!
!     Some helpful integers
!
      n_v_packed = wf%n_v*(wf%n_v+1)/2
      n_o_packed = wf%n_o*(wf%n_o+1)/2
!
      req0 = 2*(n_v_packed)*(n_o_packed)
!
      req1_a = wf%integrals%n_J*wf%n_v 
      req1_b = wf%integrals%n_J*wf%n_v 
!
      rec2 = wf%n_v**2 + 2*(n_o_packed) + 2*(n_v_packed)
!
!     Initialize batching variables
!
      batch_a = batching_index(wf%n_v)
      batch_b = batching_index(wf%n_v)
!
      call mem%batch_setup(batch_a, batch_b, req0, req1_a, req1_b, rec2)
!
      if (batch_a%num_batches /= batch_b%num_batches) then ! Should not happen, but just to be safe 
!
        call output%error_msg('Expected same-sized batches in omega_ccsd_a2_complex but something went wrong in mem%batch_setup.')
!
      endif 
!
!     Start looping over a-batches
!
      do current_a_batch = 1, batch_a%num_batches
!
         call batch_a%determine_limits(current_a_batch)
         batch_a_packed = batch_a%length*(batch_a%length+1)/2
!
         do current_b_batch = 1, current_a_batch
!
            call batch_b%determine_limits(current_b_batch)
!
            call mem%alloc(g_acbd, batch_a%length, wf%n_v, batch_b%length, wf%n_v)
!
            call ccsd_a2_integral_timer%turn_on()
!
            call wf%get_vvvv_complex(g_acbd,         &
                              batch_a%first, &
                              batch_a%last,  &
                              1,             &
                              wf%n_v,        &
                              batch_b%first, &
                              batch_b%last,  &
                              1,             &
                              wf%n_v)
!
            call ccsd_a2_integral_timer%freeze()
!
            if (current_b_batch .eq. current_a_batch) then
!
!              Allocate for +-g, +-t
!
               call mem%alloc(g_p_abcd, batch_a_packed, n_v_packed)
               call mem%alloc(g_m_abcd, batch_a_packed, n_v_packed)
               call mem%alloc(t_p_cdij, n_v_packed, n_o_packed)
               call mem%alloc(t_m_cdij, n_v_packed, n_o_packed)
!
!              Reorder g_ca_db to g_abcd and t_cidj to t_cdij
!
!$omp parallel do private(a,b,c,d,ab,cd,diag_factor)
               do c = 1, wf%n_v
                  do d = 1, c
!
                     cd = (c*(c-3)/2) + c + d
!
                     if (c .ne. d) then
                        diag_factor = two_complex
                     else
                        diag_factor = one_complex
                     endif
!
                     do b = 1, batch_b%length
                        do a = 1, b
!
                           ab = (b*(b-3)/2) + a + b
!
                           g_p_abcd(ab, cd) = diag_factor*(g_acbd(a, c, b, d) + g_acbd(a, d, b, c))
                           g_m_abcd(ab, cd) = diag_factor*(g_acbd(a, d, b, c) - g_acbd(a, c, b, d)) !a and b and c and d switched
!
                        enddo
                     enddo
                  enddo
               enddo
!$omp end parallel do
!
!$omp parallel do private(i,j,c,d,ij)
              do i = 1, wf%n_o
                 do j = 1, i
!
                     ij = (i*(i-3)/2) + i + j 
!
                     do c = 1, wf%n_v
                        do d = 1, c
!
                           cd = (c*(c-3)/2) + c + d
!
                           t_p_cdij(cd, ij) = t_abij(c, d, i, j) + t_abij(d, c, i, j)
                           t_m_cdij(cd, ij) = t_abij(c, d, i, j) - t_abij(d, c, i, j)
!
                       enddo
                    enddo
                 enddo
              enddo
!$omp end parallel do
!
!              Dellocate g_acbd
!
               call mem%dealloc(g_acbd, batch_a%length, wf%n_v, batch_b%length, wf%n_v)
!
!              Allocate omega +-
!
               call mem%alloc(omega2_p_abij, batch_a_packed, n_o_packed)
               call mem%alloc(omega2_m_abij, batch_a_packed, n_o_packed)
!
!              omega2_ab_ij = sum_(cd) g_abcd*t_cdij
!
               call zgemm('N','N',        &
                          batch_a_packed, &
                          n_o_packed,     &
                          n_v_packed,     &
                          one_complex/four_complex,       &
                          g_p_abcd,       &
                          batch_a_packed, &
                          t_p_cdij,       &
                          n_v_packed,     &
                          zero_complex,           &
                          omega2_p_abij,  &
                          batch_a_packed)
!
               call zgemm('N','N',        &
                          batch_a_packed, &
                          n_o_packed,     &
                          n_v_packed,     &
                          one_complex/four_complex,       &
                          g_m_abcd,       &
                          batch_a_packed, &
                          t_m_cdij,       &
                          n_v_packed,     &
                          zero_complex,           &
                          omega2_m_abij,  &
                          batch_a_packed )
!
!             Deallocate +-g, +-t
!
              call mem%dealloc(g_p_abcd, batch_a_packed, n_v_packed)
              call mem%dealloc(g_m_abcd, batch_a_packed, n_v_packed)
              call mem%dealloc(t_p_cdij, n_v_packed, n_o_packed)
              call mem%dealloc(t_m_cdij, n_v_packed, n_o_packed)
!
!$omp parallel do private(i, j, a, b, ij, ai, aj, bj, bi, ab, aibj, biaj)
               do i = 1, wf%n_o
                  do j = 1, i
!
                     ij = (i*(i-3)/2) + i + j
!
                     do a = 1, batch_a%length
!
                        
                        a_full = a + batch_a%first - 1
!
                        do b = 1, a
!
                           ab = (a*(a-3)/2) + a + b 
                           b_full = b + batch_b%first - 1
!
                           if (a .ne. b .and. i .ne. j) then
!
                              omega_abij(a_full, b_full, i, j) = omega_abij(a_full, b_full, i, j) &
                                          + omega2_p_abij(ab, ij) + omega2_m_abij(ab, ij)
!
                              omega_abij(b_full, a_full, i, j) = omega_abij(b_full, a_full, i, j) &
                                           + omega2_p_abij(ab, ij) - omega2_m_abij(ab, ij)
!
                              omega_abij(b_full, a_full, j, i) = omega_abij(b_full, a_full, j, i) &
                                           + omega2_p_abij(ab, ij) + omega2_m_abij(ab, ij)
!
                              omega_abij(a_full, b_full, j, i) = omega_abij(a_full, b_full, j, i) &
                                          + omega2_p_abij(ab, ij) - omega2_m_abij(ab, ij)
!
                           elseif (a == b .and. i .ne. j) then
!
                              omega_abij(a_full, b_full, i, j) = omega_abij(a_full, b_full, i, j) &
                                          + omega2_p_abij(ab, ij)
!
                              omega_abij(a_full, b_full, j, i) = omega_abij(a_full, b_full, j, i) &
                                          + omega2_p_abij(ab, ij)
!
                           elseif (a .ne. b .and. i == j) then
!
                              omega_abij(a_full, b_full, i, j) = omega_abij(a_full, b_full, i, j) &
                                          + omega2_p_abij(ab, ij)
!
                              omega_abij(b_full, a_full, i, j) = omega_abij(b_full, a_full, i, j) &
                                           + omega2_p_abij(ab, ij)
!
                           elseif (a == b .and. i == j) then
!
                              omega_abij(a_full, b_full, i, j) = omega_abij(a_full, b_full, i, j) &
                                          + omega2_p_abij(ab, ij) + omega2_m_abij(ab, ij)
!
                           endif
!
                        enddo
                     enddo
                  enddo
               enddo
!$omp end parallel do
!
!              Deallocate omega +-
!
               call mem%dealloc(omega2_p_abij, batch_a_packed, n_o_packed)
               call mem%dealloc(omega2_m_abij, batch_a_packed, n_o_packed)
!
            else
!
!              Allocate for +-g, +-t
!
               call mem%alloc(g_p_abcd, (batch_a%length)*(batch_b%length), n_v_packed)
               call mem%alloc(g_m_abcd, (batch_a%length)*(batch_b%length), n_v_packed)
               call mem%alloc(t_p_cdij, n_v_packed, n_o_packed)
               call mem%alloc(t_m_cdij, n_v_packed, n_o_packed)
!
!$omp parallel do private(a,b,c,d,ab,cd,diag_factor)
               do c = 1, wf%n_v
                  do d = 1, c
!
                     cd = (c*(c-3)/2) + c + d
!
                     if (c .ne. d) then
                        diag_factor = two_complex
                     else
                        diag_factor = one_complex
                     endif
!
                     do  b = 1, batch_b%length
                        do a = 1, batch_a%length
!
                           ab = (b-1)*batch_a%length + a
!
                           g_p_abcd(ab, cd) = diag_factor*(g_acbd(a, c, b, d) + g_acbd(a, d, b, c))
                           g_m_abcd(ab, cd) = diag_factor*(g_acbd(a, c, b, d) - g_acbd(a, d, b, c))
!
                        enddo
                     enddo
                  enddo
               enddo
!$omp end parallel do
!
!$omp parallel do schedule(static) private(c,d,i,j,cd,ij,ci,cj,di,dj,cidj,dicj)
               do i = 1, wf%n_o
                  do j = 1, i
!
                     ij = (i*(i-3)/2) + i + j 
!
                     do c = 1, wf%n_v
                        do d = 1, c
!
                           cd = (c*(c-3)/2) + c + d
!
                           t_p_cdij(cd, ij) = t_abij(c, d, i, j) + t_abij(d, c, i, j)
                           t_m_cdij(cd, ij) = t_abij(c, d, i, j) - t_abij(d, c, i, j)
!
                        enddo
                     enddo
                  enddo
               enddo
!$omp end parallel do
!
!              Dellocate g_acbd
!
               call mem%dealloc(g_acbd, batch_a%length, wf%n_v, batch_b%length, wf%n_v)
!
!              Allocate omega +-
!
               call mem%alloc(omega2_p_abij, (batch_a%length)*(batch_b%length), n_o_packed)
               call mem%alloc(omega2_m_abij, (batch_a%length)*(batch_b%length), n_o_packed)
!
!              omega2_ab_ij = sum_(cd) g_abcd*t_cdij
!
               call zgemm('N','N',                            &
                           (batch_a%length)*(batch_b%length), &
                           n_o_packed,                        &
                           n_v_packed,                        &
                           one_complex/four_complex,                          &
                           g_p_abcd,                          &
                           (batch_a%length)*(batch_b%length), &
                           t_p_cdij,                          &
                           n_v_packed,                        &
                           zero_complex,                              &
                           omega2_p_abij,                     &
                           (batch_a%length)*(batch_b%length))
!
               call zgemm('N','N',                            &
                           (batch_a%length)*(batch_b%length), &
                           n_o_packed,                        &
                           n_v_packed,                        &
                           one_complex/four_complex,                          &
                           g_m_abcd,                          &
                           (batch_a%length)*(batch_b%length), &
                           t_m_cdij,                          &
                           n_v_packed,                        &
                           zero_complex,                              &
                           omega2_m_abij,                     &
                           (batch_a%length)*(batch_b%length))
!
!              Deallocate +-g, +-t
!
               call mem%dealloc(g_p_abcd, (batch_a%length)*(batch_b%length), n_v_packed)
               call mem%dealloc(g_m_abcd, (batch_a%length)*(batch_b%length), n_v_packed)
               call mem%dealloc(t_p_cdij, n_v_packed, n_o_packed)
               call mem%dealloc(t_m_cdij, n_v_packed, n_o_packed)
!
!$omp parallel do private(i, j, a, b, ij, ai, aj, bj, bi, ab, aibj, biaj)
               do i = 1, wf%n_o
                  do j = 1, i
!
                     ij = (i*(i-3)/2) + i + j
!
                     do a = 1, batch_a%length

                        a_full = a + batch_a%first - 1
!
                        do b = 1, batch_b%length
!
                           ab = batch_a%length*(b - 1) + a
                           b_full = b + batch_b%first - 1
!
                           if (i .ne. j) then
!
                              omega_abij(a_full, b_full, i, j) = omega_abij(a_full, b_full, i, j) &
                                          + omega2_p_abij(ab, ij) + omega2_m_abij(ab, ij)
!
                              omega_abij(b_full, a_full, i, j) = omega_abij(b_full, a_full, i, j) &
                                           + omega2_p_abij(ab, ij) - omega2_m_abij(ab, ij)
!
                              omega_abij(b_full, a_full, j, i) = omega_abij(b_full, a_full, j, i) &
                                           + omega2_p_abij(ab, ij) + omega2_m_abij(ab, ij)
!
                              omega_abij(a_full, b_full, j, i) = omega_abij(a_full, b_full, j, i) &
                                          + omega2_p_abij(ab, ij) - omega2_m_abij(ab, ij)
!
                           elseif (i == j) then
!
                              omega_abij(a_full, b_full, i, j) = omega_abij(a_full, b_full, i, j) &
                                          + omega2_p_abij(ab, ij)
!
                              omega_abij(b_full, a_full, i, j) = omega_abij(b_full, a_full, i, j) &
                                           + omega2_p_abij(ab, ij)
!
                           endif
!
                        enddo
                     enddo
                  enddo
               enddo
!$omp end parallel do
!
!              Deallocate omega +-
!
               call mem%dealloc(omega2_p_abij, (batch_a%length)*(batch_b%length), n_o_packed)
               call mem%dealloc(omega2_m_abij, (batch_a%length)*(batch_b%length), n_o_packed)
!
            endif
!
         enddo ! End batching over b
      enddo ! End batching over a
!
      call ccsd_a2_timer%turn_off()
      call ccsd_a2_integral_timer%turn_off()
!
   end subroutine omega_ccsd_a2_ccsd_complex
!
!
   module subroutine omega_ccsd_b2_ccsd_complex(wf, omega_abij, t_abij)
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
      complex(dp), dimension(wf%n_v,wf%n_v,wf%n_o,wf%n_o), intent(inout):: omega_abij
      complex(dp), dimension(wf%n_v,wf%n_v,wf%n_o,wf%n_o), intent(in):: t_abij
!
!     Integrals
!
      complex(dp), dimension(:,:,:,:), allocatable :: g_aibj
      complex(dp), dimension(:,:,:,:), allocatable :: g_kcld
      complex(dp), dimension(:,:,:,:), allocatable :: g_klcd
      complex(dp), dimension(:,:,:,:), allocatable :: g_klij
      complex(dp), dimension(:,:,:,:), allocatable :: g_kilj
!
!     Reordered T2 apmlitudes
!
     ! complex(dp), dimension(:,:,:,:), allocatable :: t_cdij
!
!     Reordered omega
!
   !   complex(dp), dimension(:,:,:,:), allocatable :: omega_abij
   !   complex(dp), dimension(:,:,:,:), allocatable :: omega_aibj
!
      type(timings) :: ccsd_b2_timer
!
      ccsd_b2_timer = timings('omega ccsd b2')
      call ccsd_b2_timer%turn_on()
!
!     Construct g_aibj and add to omega2 
!
      call mem%alloc(g_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call wf%get_vovo_complex(g_aibj)
      call add_1324_to_1234(one_complex, g_aibj, omega_abij,  wf%n_v,  wf%n_v,  wf%n_o,  wf%n_o)
      call mem%dealloc(g_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Allocate and construct g_kilj
!
      call mem%alloc(g_kilj, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call wf%get_oooo_complex(g_kilj)
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
      call wf%get_ovov_complex(g_kcld)
!
!     Reorder g_kcld as g_klcd
!
      call mem%alloc(g_klcd,  wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
      call sort_1234_to_1324(g_kcld, g_klcd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%dealloc(g_kcld, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     Reorder t_cidj as t_cdij
!
      call zgemm('N','N',      &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  (wf%n_v)**2, &
                  one_complex,         &
                  g_klcd,      &
                  (wf%n_o)**2, &
                  t_abij,      &
                  (wf%n_v)**2, &
                  one_complex,         &
                  g_klij,      &
                  (wf%n_o)**2)
!
!     Deallocate t_cdij and g_klcd
!
      call mem%dealloc(g_klcd,  wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
!     omega_abij = sum_(kl) t_ab_kl*X_kl_ij
!
    !  call mem%alloc(omega_abij,  wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call zgemm('N','N',      &
                  (wf%n_v)**2, &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  one_complex,         &
                  t_abij,      & ! t_ab_kl
                  (wf%n_v)**2, &
                  g_klij,      &
                  (wf%n_o)**2, &
                  one_complex,         &
                  omega_abij,  &
                  (wf%n_v)**2)
!
      call mem%dealloc(g_klij,  wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call ccsd_b2_timer%turn_off()
!
   end subroutine omega_ccsd_b2_ccsd_complex
!
!
   module subroutine omega_ccsd_c2_ccsd_complex(wf, omega_aibj, t_aibj)
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
      complex(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(inout):: omega_aibj
      complex(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(in):: t_aibj
!
!     Integrals
!
      complex(dp), dimension(:,:,:,:), allocatable :: g_kdlc
      complex(dp), dimension(:,:,:,:), allocatable :: g_dlck
      complex(dp), dimension(:,:,:,:), allocatable :: g_kiac
      complex(dp), dimension(:,:,:,:), allocatable :: g_aick
!
!     Reordered T2 amplitudes
!
      complex(dp), dimension(:,:,:,:), allocatable :: t_aidl
      complex(dp), dimension(:,:,:,:), allocatable :: t_ckbj
!
!    Intermediates for matrix multiplication
!
      complex(dp), dimension(:,:,:,:), allocatable :: X_aick
      complex(dp), dimension(:,:,:,:), allocatable :: Y_aibj
!
!     Reordered U2 amplitudes
!
      complex(dp), dimension(:,:,:,:), allocatable :: omega_a_batch
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
      ccsd_c2_timer = timings('omega ccsd c2')
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
      call wf%get_ovov_complex(g_kdlc)
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
      call zgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  -half_complex,             &
                  t_aidl,            &
                  (wf%n_o)*(wf%n_v), &
                  g_dlck,            &
                  (wf%n_o)*(wf%n_v), &
                  zero_complex,              &
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
      req0 = wf%n_o**2*wf%integrals%n_J
!
      req1 = wf%n_v*wf%integrals%n_J + (wf%n_o)*(wf%n_v**2)
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
         call wf%get_oovv_complex(g_kiac,                        &
                           1, wf%n_o,                    &
                           1, wf%n_o,                    &
                           batch_a%first, batch_a%last,  &
                           1, wf%n_v)
!
!        X_aick = X_aick + g_kiac
!
!$omp parallel do private(a, i, k, c)
         do k = 1, wf%n_o
            do c = 1, wf%n_v
               do i = 1, wf%n_o
                  do a = 1, batch_a%length
!
                     X_aick(a+batch_a%first-1, i, c, k) = X_aick(a+batch_a%first-1, i, c, k) + g_kiac(k, i, a, c)
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
         call zgemm('N','N',                 &
                  (wf%n_o)*(batch_a%length), &
                  (wf%n_o)*(wf%n_v),         &
                  (wf%n_o)*(wf%n_v),         &
                  -one_complex/two_complex,                  &
                  g_aick,                    &
                  (wf%n_o)*(batch_a%length), &
                  wf%u_aibj_complex,                 & ! u_ck_bj
                  (wf%n_o)*(wf%n_v),         &
                  zero_complex,                      &
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
      call zgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  -one_complex,              &
                  X_aick,            &
                  (wf%n_o)*(wf%n_v), &
                  t_ckbj,            &
                  (wf%n_o)*(wf%n_v), &
                  zero_complex,              &
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
                     omega_aibj(a, i, b, j) = omega_aibj(a, i, b, j) + half_complex*Y_aibj(a, i, b, j) + Y_aibj(a, j, b, i) 
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
   end subroutine omega_ccsd_c2_ccsd_complex
!
!
   module subroutine omega_ccsd_d2_ccsd_complex(wf, omega_aibj, t_aibj)
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
      complex(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(inout):: omega_aibj
      complex(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(in):: t_aibj
!
      complex(dp), dimension(:,:,:,:), allocatable :: g_ldkc ! g_ldkc
      complex(dp), dimension(:,:,:,:), allocatable :: L_ldkc ! L_ldkc = 2 * g_ldkc - g_lckd
      complex(dp), dimension(:,:,:,:), allocatable :: u_aild ! u_il^ad = 2 * t_il^ad - t_li^ad
      complex(dp), dimension(:,:,:,:), allocatable :: Z_aikc ! An intermediate, see below
      complex(dp), dimension(:,:,:,:), allocatable :: g_aikc ! g_aikc
!
      type(timings) :: ccsd_d2_timer
!
      ccsd_d2_timer = timings('omega ccsd d2')
      call ccsd_d2_timer%turn_on()
!
!     :: Calculate the D2.2 term of omega ::
!
!     Form L_ld_kc = L_ldkc = 2*g_ldkc(ld,kc) - g_ldkc(lc,kd)
!
      call mem%alloc(g_ldkc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call wf%get_ovov_complex(g_ldkc)
!
      call mem%alloc(L_ldkc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call copy_and_scale_complex(two_complex, g_ldkc, L_ldkc, (wf%n_o)**2*(wf%n_v)**2)
      call add_1432_to_1234(-one_complex, g_ldkc, L_ldkc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%dealloc(g_ldkc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%alloc(u_aild, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      
      call zero_array_complex(u_aild, wf%n_t1**2)
      call add_1243_to_1234(two_complex, t_aibj, u_aild, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call add_1342_to_1234(-one_complex, t_aibj, u_aild, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
!     Form the intermediate Z_aikc = sum_dl u_aild L_ldkc and set it to zero
!
      call mem%alloc(Z_aikc,  wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call zgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one_complex,               &
                  u_aild,            &
                  (wf%n_o)*(wf%n_v), &
                  L_ldkc,            &
                  (wf%n_o)*(wf%n_v), &
                  zero_complex,              &
                  Z_aikc,            &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(L_ldkc,  wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     Form the D2.2 term, 1/4 sum_kc Z_aikc u_kc_bj = 1/4 sum_kc Z_aikc(ai,kc) u_aild(bj,kc)
!
      call zgemm('N','T',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one_complex/four_complex,          &
                  Z_aikc,            &
                  (wf%n_o)*(wf%n_v), &
                  u_aild,            &
                  (wf%n_o)*(wf%n_v), &
                  one_complex,               &
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
      call mem%alloc(g_aikc, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call wf%get_voov_complex(g_aikc)
!
!     Calculate the D2.1 term, sum_ck u_jk^bc g_aikc = sum_ck g_aikc(ai,kc) u_aild(bj,kc)
!
      call zgemm('N','T',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one_complex,               &
                  g_aikc,            &
                  (wf%n_o)*(wf%n_v), &
                  u_aild,            &
                  (wf%n_o)*(wf%n_v), &
                  one_complex,               &
                  omega_aibj,        &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(u_aild, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call mem%dealloc(g_aikc, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call ccsd_d2_timer%turn_off()
!
   end subroutine omega_ccsd_d2_ccsd_complex
!
!
   module subroutine omega_ccsd_e2_ccsd_complex(wf, omega_aibj, t_aibj)
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
      complex(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(inout):: omega_aibj
      complex(dp), dimension(wf%n_v,wf%n_o,wf%n_v,wf%n_o), intent(in):: t_aibj
!
!     Vectors for E2.1 term
!
      complex(dp), dimension(:,:,:,:), allocatable :: g_ldkc      ! g_ldkc
      complex(dp), dimension(:,:,:,:), allocatable :: u_bldk      ! u_kl^bd
      complex(dp), dimension(:,:), allocatable :: X_b_c           ! An intermediate, see below for definition
!
!     Vectors for E2.2 term
!
      complex(dp), dimension(:,:), allocatable :: Y_k_j        ! An intermediate, see below for definition
!
      type(timings) :: ccsd_e2_timer 
!
      ccsd_e2_timer = timings('omega ccsd e2')
      call ccsd_e2_timer%turn_on()
!
!     :: Calculate the E2.1 term of omega ::
!
!     Form u_bldk = u_kl^bd = u_bkdl
!                  = 2 * t_kl^bd - t_lk^bd = 2 * t_bkdl(bk,dl) - t_bkdl(bl, dk)
!
      call mem%alloc(u_bldk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      
      call copy_and_scale_complex(-one_complex, t_aibj, u_bldk, (wf%n_v)**2*(wf%n_o)**2)
      call add_1432_to_1234(two_complex, t_aibj, u_bldk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Form g_ldkc = g_ldkc
!
      call mem%alloc(g_ldkc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call wf%get_ovov_complex(g_ldkc)
!
!     Make the intermediate X_b_c = F_bc - sum_dkl g_ldkc u_kl^bd and set to zero
!
      call mem%alloc(X_b_c, wf%n_v, wf%n_v)
!
!     First, copy the virtual-virtual Fock matrix into the intermediate
!
      call zcopy((wf%n_v)**2, wf%fock_ab_complex, 1, X_b_c, 1) ! X_b_c = F_bc
!
!     Then, add the second contribution,
!     - sum_dkl g_ldkc u_kl^bd = - sum_dkl u_b_kdl * g_kdl_c, to X_b_c
!
      call zgemm('N','N',               &
                  wf%n_v,               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  -one_complex,                 &
                  u_bldk,               & ! u_b_ldk
                  wf%n_v,               &
                  g_ldkc,               & ! g_ldk_c
                  (wf%n_v)*(wf%n_o)**2, &
                  one_complex,                  &
                  X_b_c,                &
                  wf%n_v)
!
!     Form t_cjai = t_ij^ac = t_ji^ca
!
!     Form the E2.1 term
!
      call zgemm('N','N',               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_v,               &
                  one_complex,                  &
                  X_b_c,                &
                  wf%n_v,               &
                  t_aibj,               & ! t_c_jai
                  wf%n_v,               &
                  one_complex,                  &
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
      call zcopy((wf%n_o)**2, wf%fock_ij_complex, 1, Y_k_j, 1)
!
!     Add sum_cdl g_k_dlc u_dlc_j to Y_k_j, such that
!     Y_k_j = F_k_j + sum_cdl g_k_cld u_cld_j
!
!     Note:
!
!     g_ldkc(kc,ld) = g_kcld              -> pretend that this is g_k_cld
!     u_b_ldk(c,ldj) = u_jl^cd (= u_lj^dc) -> pretend that this is u_cld_j
!
      call zgemm('N','N',              &
                 wf%n_o,               &
                 wf%n_o,               &
                 (wf%n_o)*(wf%n_v)**2, &
                 one_complex,                  &
                 g_ldkc,               & ! g_k_cld
                 wf%n_o,               &
                 u_bldk,               & ! u_cld_j
                 (wf%n_o)*(wf%n_v)**2, &
                 one_complex,                  &
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
      call zgemm('N','N',              &
                 (wf%n_o)*(wf%n_v)**2, &
                 wf%n_o,               &
                 wf%n_o,               &
                 -one_complex,                 &
                 t_aibj,               & ! t_aib_k
                 (wf%n_o)*(wf%n_v)**2, &
                 Y_k_j,                &
                 wf%n_o,               &
                 one_complex,                  &
                 omega_aibj,           & ! omega2_aib_j
                 (wf%n_o)*(wf%n_v)**2)
!
!     Deallocate Y_k_j and the amplitudes
!
      call mem%dealloc(Y_k_j, wf%n_o, wf%n_o)
!
      call ccsd_e2_timer%turn_off()
!
   end subroutine omega_ccsd_e2_ccsd_complex
!
!
   module subroutine construct_u_aibj_ccsd_complex(wf)
!!
!!    Construct u_aibj
!!    Written by Tor S. Haugland, Nov 2019
!!
!!    Construct
!!       u_aibj = 2t_aibj - t_ajbi
!!
      implicit none
!
      class(ccsd),                  intent(inout)          :: wf
!
      complex(dp), dimension(:,:,:,:), allocatable :: t_aibj
!
      type(timings) :: timer
!
      timer = timings('Construct u_aibj', pl='debug')
      call timer%turn_on()
!
      call mem%alloc(t_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2_complex, t_aibj, wf%n_t1)
!
      call copy_and_scale_complex(two_complex,    t_aibj, wf%u_aibj_complex, wf%n_v**2 * wf%n_o**2)
      call add_1432_to_1234(-one_complex, t_aibj, wf%u_aibj_complex, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(t_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine construct_u_aibj_ccsd_complex
!
!
end submodule omega_ccsd_complex
