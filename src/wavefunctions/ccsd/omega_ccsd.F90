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
submodule (ccsd_class) omega_ccsd
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
      real(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: omega
!
      real(dp), dimension(:,:), allocatable :: omega1
      real(dp), dimension(:), allocatable :: omega2
!
      call mem%alloc(omega1, wf%n_v, wf%n_o)
      call mem%alloc(omega2, (wf%n_v)*(wf%n_o)*((wf%n_v)*(wf%n_o) +1)/2)
!
!     Set the omega vector to zero
!
      omega1 = zero
      omega2 = zero
!
!     Construct singles contributions
!
      call wf%omega_ccsd_a1(omega1)
      call wf%omega_ccsd_b1(omega1)
      call wf%omega_ccsd_c1(omega1)
!
      call wf%omega_ccs_a1(omega1)
!
!     Construct doubles contributions
!
      call wf%omega_ccsd_a2(omega2)
      call wf%omega_ccsd_b2(omega2)
      call wf%omega_ccsd_c2(omega2)
      call wf%omega_ccsd_d2(omega2)
      call wf%omega_ccsd_e2(omega2)
!
      call wf%from_biorthogonal_to_biorthonormal(omega2)
!
      call dcopy((wf%n_o)*(wf%n_v), omega1, 1, omega, 1)
      call dcopy((wf%n_o)*(wf%n_v)*((wf%n_o)*(wf%n_v)+1)/2, omega2, 1, omega((wf%n_o)*(wf%n_v)+1), 1)
!
      call mem%dealloc(omega1, wf%n_v, wf%n_o)
      call mem%dealloc(omega2, (wf%n_v)*(wf%n_o)*((wf%n_v)*(wf%n_o) +1)/2)
!
   end subroutine construct_omega_ccsd
!
!
   module subroutine omega_ccsd_a1_ccsd(wf, omega1)
!!
!!    Omega A1 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Calculates the A1 term,
!!
!!       A1: sum_ckd g_adkc * u_ki^cd,
!!
!!    and adds it to the singles projection vector (omega1) of
!!    the wavefunction object wf.
!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout):: omega1
!
      integer :: current_a_batch = 0
!
      type(batching_index) :: batch_a
!
      real(dp), dimension(:,:,:,:), allocatable :: u_dkci, t_dkci
      real(dp), dimension(:,:,:,:), allocatable :: g_adkc
!
      integer :: req0, req1
!
      type(timings) :: ccsd_a1_timer
!
      call ccsd_a1_timer%init('omega ccsd a1')
      call ccsd_a1_timer%start()
!
!     u_ki^cd = 2*t_ki^cd - t_ik^cd (ordered as u_dkci)
!
      call mem%alloc(t_dkci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_dkci, (wf%n_v)*(wf%n_o))
!
      call mem%alloc(u_dkci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      u_dkci = -t_dkci
      call add_1432_to_1234(two, t_dkci, u_dkci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(t_dkci,  wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Batch over a to hold g_adkc
!
      req0 = wf%n_o*wf%integrals%n_J*wf%n_v
!
      req1 = wf%n_v*wf%integrals%n_J + wf%n_v**2*(wf%n_o)
!
      call batch_a%init(wf%n_v)
!
      call mem%batch_setup(batch_a, req0, req1)
!
      do current_a_batch = 1, batch_a%num_batches
!
         call batch_a%determine_limits(current_a_batch)
!
         call mem%alloc(g_adkc, batch_a%length, wf%n_v, wf%n_o, wf%n_v)
!
         call wf%get_vvov(g_adkc,                        &
                           batch_a%first, batch_a%last,  &
                           1, wf%n_v,                    &
                           1, wf%n_o,                    &
                           1, wf%n_v)
!
         call dgemm('N','N',                     &
                     batch_a%length,             &
                     wf%n_o,                     &
                     (wf%n_o)*(wf%n_v)**2,       &
                     one,                        &
                     g_adkc,                     & ! g_a_dkc
                     batch_a%length,             &
                     u_dkci,                     & ! u_dkc_i
                     (wf%n_o)*(wf%n_v)**2,       &
                     one,                        &
                     omega1(batch_a%first, 1),   &
                     wf%n_v)
!
         call mem%dealloc(g_adkc, batch_a%length, wf%n_v, wf%n_o, wf%n_v)
!
      enddo ! End of batches of the index a
!
      call mem%dealloc(u_dkci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call ccsd_a1_timer%freeze()
      call ccsd_a1_timer%switch_off()
!
   end subroutine omega_ccsd_a1_ccsd
!
!
   module subroutine omega_ccsd_b1_ccsd(wf, omega1)
!!
!!    Omega B1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Calculates the B1 term,
!!
!!       B1:   - sum_ckl u_kl^ac * g_kilc,
!!
!!    and adds it to the singles projection vector (omega1) of
!!    the wavefunction object wf
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout):: omega1
!
      real(dp), dimension(:,:,:,:), allocatable :: g_lcki ! g_kilc
      real(dp), dimension(:,:,:,:), allocatable :: t_alck ! t_kl^ac
      real(dp), dimension(:,:,:,:), allocatable :: u_alck ! u_kl^ac = 2 t_kl^ac - t_lk^ac
!
      type(timings) :: ccsd_b1_timer 
!  
      call ccsd_b1_timer%init('omega ccsd b1')
      call ccsd_b1_timer%start()
!
!     Form u_alck = u_kl^ac = 2 * t_kl^ac - t_lk^ac
!     Square up amplitudes and reorder: t_akcl to t_alck
!
      call mem%alloc(g_lcki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call wf%get_ovoo(g_lcki)
!
!     u_alck = 2 * t_kl^ac - t_lk^ac = 2 * t_alck(alck) - t_alck(akcl)
!
      call mem%alloc(t_alck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_alck, (wf%n_v)*(wf%n_o))
!
      call mem%alloc(u_alck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      u_alck = zero
      call daxpy(((wf%n_v)*(wf%n_o))**2, -one, t_alck, 1, u_alck, 1)
!
      call add_1432_to_1234(two, t_alck, u_alck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(t_alck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',                 &
                  wf%n_v,                 &
                  wf%n_o,                 &
                  (wf%n_v)*((wf%n_o)**2), &
                  -one,                   &
                  u_alck,                 & ! u_a_lck
                  wf%n_v,                 &
                  g_lcki,                 & ! g_lck_i
                  (wf%n_v)*((wf%n_o)**2), &
                  one,                    &
                  omega1,                 &
                  wf%n_v)
!
      call mem%dealloc(u_alck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(g_lcki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call ccsd_b1_timer%freeze()
      call ccsd_b1_timer%switch_off()
!
   end subroutine omega_ccsd_b1_ccsd
!
!
   module subroutine omega_ccsd_c1_ccsd(wf, omega1)
!!
!!    Omega C1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Calculates the C1 term of omega,
!!
!!       C1: sum_ck F_kc*u_aick,
!!
!!    and adds it to the projection vector (omega1) of
!!    the wavefunction object wf
!!
!!    u_aikc = 2*t_ckai - t_ciak
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout):: omega1
!
      real(dp), dimension(:,:), allocatable :: F_c_k ! F_kc
!
      real(dp), dimension(:,:,:,:), allocatable :: u_aick
      real(dp), dimension(:,:,:,:), allocatable :: t_aick
!
      type(timings) :: ccsd_c1_timer 
!
      call ccsd_c1_timer%init('omega ccsd c1')
      call ccsd_c1_timer%start()
!
!     Form u_aick = u_ik^ac = 2*t_ik^ac - t_ki^ac
!
      call mem%alloc(t_aick,  wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_aick, (wf%n_o)*(wf%n_v))
!
      call mem%alloc(u_aick,  wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      u_aick = zero
!
      call add_1432_to_1234(-one, t_aick, u_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call daxpy((wf%n_o)**2*(wf%n_v)**2, two, t_aick, 1, u_aick, 1)
!
      call mem%dealloc(t_aick,  wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Reorder the Fock matrix F_ck = F_kc
!
      call mem%alloc(F_c_k, wf%n_v, wf%n_o)
!
      call sort_12_to_21(wf%fock_ia, F_c_k, wf%n_o, wf%n_v)
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  1,                 &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  u_aick,            &
                  (wf%n_o)*(wf%n_v), &
                  F_c_k,             & ! Equal to F_kc
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  omega1,            &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(u_aick,  wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(F_c_k, wf%n_v, wf%n_o)
!
      call ccsd_c1_timer%freeze()
      call ccsd_c1_timer%switch_off()
!
   end subroutine omega_ccsd_c1_ccsd
!
!
   module subroutine omega_ccsd_a2_ccsd(wf, omega2)
!!
!!    Omega A2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!      
!!    A2 = g_aibj + sum_(cd)g_acbd * t_cidj = A2.1 + A.2.2
!!
!!    Structure: Batching over both a and b for A2.2.
!!                t^+_ci_dj = t_cidj + t_di_cj
!!                t^-_ci_dj = t_cidj - t_di_cj
!!                g^+_ac_bd = g_acbd + g_bc_ad
!!                g^-_ac_bd = g_acbd - g_bc_ad
!!
!!                omega_A2.2_ai_bj = 1/4*(g^+_ac_bd*t^+_ci_dj + g^-_ac_bd*t^-_ci_dj) = omega_A2.2_bj_ai
!!                omega_A2.2_aj_bi = 1/4*(g^+_ac_bd*t^+_ci_dj - g^-_ac_bd*t^-_ci_dj) = omega_A2.2_bi_aj
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_v)*(wf%n_o)*((wf%n_v)*(wf%n_o)+1)/2), intent(inout):: omega2
!
!     Integrals
!
      real(dp), dimension(:,:,:,:), allocatable :: g_aibj
      real(dp), dimension(:,:,:,:), allocatable :: g_acbd
      real(dp), dimension(:,:), allocatable :: g_p_abcd
      real(dp), dimension(:,:), allocatable :: g_m_abcd
!
!     Reordered T2 amplitudes
!
      real(dp), dimension(:,:), allocatable :: t_p_cdij
      real(dp), dimension(:,:), allocatable :: t_m_cdij
!
!     Reordered omega 2
!
      real(dp), dimension(:,:), allocatable :: omega2_p_abij
      real(dp), dimension(:,:), allocatable :: omega2_m_abij
!
!     Indices
!
      integer :: a, b, c, d
      integer :: i, j
!
      integer :: ab, cd
      integer :: ai, aj, bj, bi, ci, cj, dj, di
      integer :: ij
!
      integer :: aibj, biaj, cidj, dicj 
!
      real(dp) :: diag_factor
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
      call ccsd_a2_timer%init('omega ccsd a2')
      call ccsd_a2_integral_timer%init('omega ccsd a2 g_abcd')
!      
      call ccsd_a2_timer%start()
!
!     Some helpful integers
!
      n_v_packed = wf%n_v*(wf%n_v+1)/2
      n_o_packed = wf%n_o*(wf%n_o+1)/2
!
!     :: Calculate the A2.1 term of omega ::
!
!     Create g_aibj
!
      call mem%alloc(g_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%get_vovo(g_aibj)
!
!     Add A2.1 to Omega 2
!
      call add_to_packed(omega2, g_aibj, (wf%n_o)*(wf%n_v))

      call mem%dealloc(g_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!    ::  Calculate the A2.2 term  of omega ::
!
      req0 = 2*(n_v_packed)*(n_o_packed)
!
      req1_a = wf%integrals%n_J*wf%n_v 
      req1_b = wf%integrals%n_J*wf%n_v 
!
      rec2 = 2*wf%n_v**2 + 2*(n_o_packed)
!
!     Initialize batching variables
!
      call batch_a%init(wf%n_v)
      call batch_b%init(wf%n_v)
!
      call mem%batch_setup(batch_a, batch_b, req0, req1_a, req1_b, rec2)
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
            call ccsd_a2_integral_timer%start()
!
            call wf%get_vvvv(g_acbd,         &
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
                        diag_factor = two
                     else
                        diag_factor = one
                     endif
!
                     do  b = 1, batch_b%length
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
!$omp parallel do private(i,j,c,d,ij,cd,ci,cj,di,dj,cidj,dicj)
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
                           ci = (i-1)*wf%n_v + c
                           dj = (j-1)*wf%n_v + d
                           cj = (j-1)*wf%n_v + c
                           di = (i-1)*wf%n_v + d
!
                           cidj = (ci*(ci-3)/2) + ci + dj
                           dicj = (max(cj,di)*(max(cj,di)-3)/2) + cj + di 
!
                           t_p_cdij(cd, ij) = wf%t2(cidj) + wf%t2(dicj)
                           t_m_cdij(cd, ij) = wf%t2(cidj) - wf%t2(dicj)
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
               call dgemm('N','N',        &
                          batch_a_packed, &
                          n_o_packed,     &
                          n_v_packed,     &
                          one/four,       &
                          g_p_abcd,       &
                          batch_a_packed, &
                          t_p_cdij,       &
                          n_v_packed,     &
                          zero,           &
                          omega2_p_abij,  &
                          batch_a_packed)
!
               call dgemm('N','N',        &
                          batch_a_packed, &
                          n_o_packed,     &
                          n_v_packed,     &
                          one/four,       &
                          g_m_abcd,       &
                          batch_a_packed, &
                          t_m_cdij,       &
                          n_v_packed,     &
                          zero,           &
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
                        ai = (i-1)*wf%n_v + (a + batch_a%first - 1)
                        aj = (j-1)*wf%n_v + (a + batch_a%first - 1)
!
                        do b = 1, a
!
!
                           bi = (i-1)*wf%n_v + (b + batch_b%first - 1)
                           bj = (j-1)*wf%n_v + (b + batch_b%first - 1)
!
                           ab = (a*(a-3)/2) + a + b 
!
                           aibj = (max(ai,bj)*(max(ai,bj)-3)/2) + ai + bj
!
!                          Reorder into omega2_aibj
!
                           omega2(aibj) = omega2(aibj) &
                                                + omega2_p_abij(ab, ij) + omega2_m_abij(ab, ij)
!
                           if (a .ne. b .and. i .ne. j) then
!
                              biaj = (max(bi,aj)*(max(bi,aj)-3)/2) + bi + aj 
                              omega2(biaj) = omega2(biaj) &
                                           + omega2_p_abij(ab, ij) - omega2_m_abij(ab, ij)
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
                        diag_factor = two
                     else
                        diag_factor = one
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
                           ci = (i-1)*wf%n_v + c
                           dj = (j-1)*wf%n_v + d
                           cj = (j-1)*wf%n_v + c
                           di = (i-1)*wf%n_v + d
!
                           cidj = (max(ci,dj)*(max(ci,dj)-3)/2) + ci + dj 
                           dicj = (max(cj,di)*(max(cj,di)-3)/2) + cj + di 
!
                           t_p_cdij(cd, ij) = wf%t2(cidj) + wf%t2(dicj)
                           t_m_cdij(cd, ij) = wf%t2(cidj) - wf%t2(dicj)
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
               call dgemm('N','N',                            &
                           (batch_a%length)*(batch_b%length), &
                           n_o_packed,                        &
                           n_v_packed,                        &
                           one/four,                          &
                           g_p_abcd,                          &
                           (batch_a%length)*(batch_b%length), &
                           t_p_cdij,                          &
                           n_v_packed,                        &
                           zero,                              &
                           omega2_p_abij,                     &
                           (batch_a%length)*(batch_b%length))
!
               call dgemm('N','N',                            &
                           (batch_a%length)*(batch_b%length), &
                           n_o_packed,                        &
                           n_v_packed,                        &
                           one/four,                          &
                           g_m_abcd,                          &
                           (batch_a%length)*(batch_b%length), &
                           t_m_cdij,                          &
                           n_v_packed,                        &
                           zero,                              &
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
!
                        ai = wf%n_v*(i - 1) + a + batch_a%first - 1 ! A is full-space a index
                        aj = wf%n_v*(j - 1) + a + batch_a%first - 1 ! A is full-space a index
!
                        do b = 1, batch_b%length
!
                           bj = wf%n_v*(j - 1) + b + batch_b%first - 1 ! B is full-space b index
                           bi = wf%n_v*(i - 1) + b + batch_b%first - 1 ! B is full-space b index
!
                           ab = batch_a%length*(b - 1) + a
!
                           aibj = max(ai, bj)*(max(ai, bj)-3)/2 + ai + bj
!
!                          Reorder into omega2_aibj
!
                           omega2(aibj) = omega2(aibj) &
                                        + omega2_p_abij(ab, ij) + omega2_m_abij(ab, ij)
!
                           if (i .ne. j) then
!
                              biaj = max(bi, aj)*(max(bi, aj)-3)/2 + bi + aj
                              omega2(biaj) = omega2(biaj) &
                                           + omega2_p_abij(ab, ij) - omega2_m_abij(ab, ij)
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
      call ccsd_a2_timer%freeze()
!
      call ccsd_a2_timer%switch_off()
      call ccsd_a2_integral_timer%switch_off()
!
   end subroutine omega_ccsd_a2_ccsd
!
!
   module subroutine omega_ccsd_b2_ccsd(wf, omega2)
!!
!!    Omega B2
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!!    Omega B2 = sum_(kl) t_ak_bl*(g_kilj + sum_(cd) t_cidj * g_kcld)
!!
!!    Structure: g_kilj is constructed first and reordered as g_klij.
!!    Then the contraction over cd is performed, and the results added to g_klij.
!!    t_ak_bl is then reordered as t_ab_kl and the contraction over kl is performed.
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_v)*(wf%n_o)*((wf%n_v)*(wf%n_o)+1)/2), intent(inout):: omega2
!
!     Integrals
!
      real(dp), dimension(:,:,:,:), allocatable :: g_kcld
      real(dp), dimension(:,:,:,:), allocatable :: g_klcd
      real(dp), dimension(:,:,:,:), allocatable :: g_klij
      real(dp), dimension(:,:,:,:), allocatable :: g_kilj
!
!     Reordered T2 apmlitudes
!
      real(dp), dimension(:,:,:,:), allocatable :: t_cdij
!
!     Reordered omega
!
      real(dp), dimension(:,:,:,:), allocatable :: omega_abij
      real(dp), dimension(:,:,:,:), allocatable :: omega_aibj
!
      type(timings) :: ccsd_b2_timer
!
      call ccsd_b2_timer%init('omega ccsd b2')
      call ccsd_b2_timer%start()
!
!     Allocate and construct g_kilj
!
      call mem%alloc(g_kilj, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call wf%get_oooo(g_kilj)
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
      call wf%get_ovov(g_kcld)
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
      call mem%alloc(t_cdij,  wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call squareup_and_sort_1234_to_1324(wf%t2, t_cdij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',      &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  (wf%n_v)**2, &
                  one,         &
                  g_klcd,      &
                  (wf%n_o)**2, &
                  t_cdij,      &
                  (wf%n_v)**2, &
                  one,         &
                  g_klij,      &
                  (wf%n_o)**2)
!
!     Deallocate t_cdij and g_klcd
!
      call mem%dealloc(g_klcd,  wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
!     omega_abij = sum_(kl) t_ab_kl*X_kl_ij
!
      call mem%alloc(omega_abij,  wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N','N',      &
                  (wf%n_v)**2, &
                  (wf%n_o)**2, &
                  (wf%n_o)**2, &
                  one,         &
                  t_cdij,      & ! t_ab_kl
                  (wf%n_v)**2, &
                  g_klij,      &
                  (wf%n_o)**2, &
                  zero,        &
                  omega_abij,  &
                  (wf%n_v)**2)
!
      call mem%dealloc(t_cdij,  wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(g_klij,  wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
!     Reorder into omega2
!
      call mem%alloc(omega_aibj,  wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_1324(omega_abij, omega_aibj, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call mem%dealloc(omega_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call add_to_packed(omega2, omega_aibj, (wf%n_o)*(wf%n_v))
      call mem%dealloc(omega_aibj,  wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call ccsd_b2_timer%freeze()
      call ccsd_b2_timer%switch_off()
!
   end subroutine omega_ccsd_b2_ccsd
!
!
   module subroutine omega_ccsd_c2_ccsd(wf, omega2)
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
      real(dp), dimension((wf%n_v)*(wf%n_o)*((wf%n_v)*(wf%n_o)+1)/2), intent(inout):: omega2
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
      real(dp), dimension(:,:,:,:), allocatable :: t_bkcj
!
!     Intermediates for matrix multiplication
!
      real(dp), dimension(:,:,:,:), allocatable :: X_aick
      real(dp), dimension(:,:,:,:), allocatable :: Y_aibj
!
!     Reordered U2 amplitudes
!
      real(dp), dimension(:,:,:,:), allocatable :: u_ckbj
      real(dp), dimension(:,:,:,:), allocatable :: omega2_aibj ! Holds term temporarily
      real(dp), dimension(:,:,:,:), allocatable :: omega_a_batch
!
!     Indices
!
      integer :: a, b, c
      integer :: i, j, k
!
      integer :: ai, aj, bi, bj
!
      integer :: aibj
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
      call ccsd_c2_timer%init('omega ccsd c2')
      call ccsd_c2_timer%start()
!
!     Sort t_al_di = t_li^ad as t_aidl (1234 to 1432)
!
      call mem%alloc(t_aidl,  wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call squareup_and_sort_1234_to_1432(wf%t2, t_aidl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Get g_kdlc
!
      call mem%alloc(g_kdlc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_ovov(g_kdlc)
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
!     Deallocate L_ia_J and g_dlck,
!
      call mem%dealloc(g_dlck,  wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(t_aidl,  wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Form u_ckbj = u_jk^bc = u_kj^cb = 2 * t_jk^bc - t_kj^bc
!
      call mem%alloc(t_ckbj,  wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_ckbj, (wf%n_o)*(wf%n_v))
!
      call mem%alloc(u_ckbj,  wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      u_ckbj = zero
!
      call add_1432_to_1234(-one, t_ckbj, u_ckbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call daxpy((wf%n_o)**2*(wf%n_v)**2, two, t_ckbj, 1, u_ckbj, 1)
!
      call mem%dealloc(t_ckbj,  wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Allocate a holder for - 1/2 * sum_ck u_jk^bc g_acki,
!     constructed in batches over the a index below
!
      call mem%alloc(omega2_aibj,  wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      omega2_aibj = zero
!
!     Constructing g_kiac
!
!     Prepare for batching
!
      req0 = wf%n_o**2*wf%integrals%n_J
!
      req1 = wf%n_v*wf%integrals%n_J + (wf%n_o)*(wf%n_v**2)
!
      call batch_a%init(wf%n_v)
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
         call wf%get_oovv(g_kiac,                        &
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
         call dgemm('N','N',                 &
                  (wf%n_o)*(batch_a%length), &
                  (wf%n_o)*(wf%n_v),         &
                  (wf%n_o)*(wf%n_v),         &
                  -one/two,                  &
                  g_aick,                    &
                  (wf%n_o)*(batch_a%length), &
                  u_ckbj,                    &
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
                     omega2_aibj(a + batch_a%first - 1, i, b, j) = omega_a_batch(a, i, b, j)               
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
!     Deallocate reordered u_ckbj vector
!
      call mem%dealloc(u_ckbj,  wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Add the - 1/2 * sum_ck u_jk^bc g_acki term to omega
!
!     Omega2(aibj, 1) =+ omega2_aibj(ai,bj) + omega2_bj_ai(bj,ai)
!
      call symmetric_sum(omega2_aibj, (wf%n_o)*(wf%n_v))         ! symmetrize
      call add_to_packed(omega2, omega2_aibj, (wf%n_o)*(wf%n_v)) ! add to packed
!
      call mem%dealloc(omega2_aibj,  wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Reorder t_bkcj = t_kj^bc as t_ckbj
!
      call mem%alloc(t_bkcj,  wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_bkcj, (wf%n_o)*(wf%n_v))
!
      call mem%alloc(t_ckbj,  wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_3214(t_bkcj, t_ckbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(t_bkcj,  wf%n_v, wf%n_o, wf%n_v, wf%n_o)
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
!$omp parallel do private(i, a, j, b, ai, bj, aj, bi, aibj)
         do i = 1, wf%n_o
            do a = 1, wf%n_v
!
               ai = wf%n_v*(i - 1) + a
!
               do j = 1, wf%n_o
                  do b = 1, wf%n_v
!
                  bj = wf%n_v*(j - 1) + b
!
                  if (ai .ge. bj) then
!
                     aj = wf%n_v*(j - 1) + a
                     bi = wf%n_v*(i - 1) + b
!
                     aibj = max(ai, bj)*(max(ai, bj)-3)/2 + ai + bj
!
                     omega2(aibj) = omega2(aibj) + half*Y_aibj(a, i, b, j) + Y_aibj(a, j, b, i) &
                                                        + half*Y_aibj(b, j, a, i) + Y_aibj(b, i, a, j)
!
                  endif
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(Y_aibj,  wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call ccsd_c2_timer%freeze()
      call ccsd_c2_timer%switch_off()
!
   end subroutine omega_ccsd_c2_ccsd
!
!
   module subroutine omega_ccsd_d2_ccsd(wf, omega2)
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
!!    All terms are added to the omega vector of the wavefunction object wf.
!!
!!    The routine adds the terms in the following order: D2.2, D2.1
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension((wf%n_v)*(wf%n_o)*((wf%n_v)*(wf%n_o)+1)/2), intent(inout):: omega2
!
      real(dp), dimension(:,:), allocatable :: omega2_aibj ! For storing D2.2 & D2.1
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ldkc ! g_ldkc
      real(dp), dimension(:,:,:,:), allocatable :: L_ldkc ! L_ldkc = 2 * g_ldkc - g_lckd
      real(dp), dimension(:,:,:,:), allocatable :: t_aidl ! t_il^ad
      real(dp), dimension(:,:,:,:), allocatable :: u_aidl ! u_il^ad = 2 * t_il^ad - t_li^ad
      real(dp), dimension(:,:,:,:), allocatable :: u_aild ! u_il^ad = 2 * t_il^ad - t_li^ad
      real(dp), dimension(:,:,:,:), allocatable :: Z_aikc ! An intermediate, see below
      real(dp), dimension(:,:,:,:), allocatable :: g_aikc ! g_aikc
!
      type(timings) :: ccsd_d2_timer
!
      call ccsd_d2_timer%init('omega ccsd d2')
      call ccsd_d2_timer%start()
!
!     :: Calculate the D2.2 term of omega ::
!
!     Form L_ld_kc = L_ldkc = 2*g_ldkc(ld,kc) - g_ldkc(lc,kd)
!
      call mem%alloc(g_ldkc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call wf%get_ovov(g_ldkc)
!
      call mem%alloc(L_ldkc,  wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      L_ldkc = zero
!
      call add_1432_to_1234(-one, g_ldkc, L_ldkc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call daxpy((wf%n_o)**2*(wf%n_v)**2, two, g_ldkc, 1, L_ldkc, 1)
!
      call mem%dealloc(g_ldkc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     Form u_aild = u_il^ad = 2 * t_il^ad - t_li^ad
!
      call mem%alloc(t_aidl,  wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_aidl, (wf%n_o)*(wf%n_v))
!
      call mem%alloc(u_aidl,  wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      u_aidl = zero
!
      call add_1432_to_1234(-one, t_aidl, u_aidl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call daxpy((wf%n_o)**2*(wf%n_v)**2, two, t_aidl, 1, u_aidl, 1)
!
      call mem%dealloc(t_aidl,  wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(u_aild,  wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_1243(u_aidl, u_aild, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(u_aidl,  wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Form the intermediate Z_aikc = sum_dl u_aild L_ldkc and set it to zero
!
      call mem%alloc(Z_aikc,  wf%n_v, wf%n_o, wf%n_v, wf%n_o)
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
      call mem%dealloc(L_ldkc,  wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Form the D2.2 term, 1/4 sum_kc Z_aikc u_kc_bj = 1/4 sum_kc Z_aikc(ai,kc) u_aild(bj,kc)
!
      call mem%alloc(omega2_aibj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
                  zero,              &
                  omega2_aibj,       &
                  (wf%n_o)*(wf%n_v))
!
!     Some justification for the above matrix multiplication. We have
!
!           1/4 * sum_ck (sum_dl u_il^ad L_ldkc) u_jk^bc = 1/4 * sum_ck Z_ai,kc u_kc,bj,
!
!     where Z_ai,kc = sum_dl u_ai,ld L_ld,kc. Note that u_aild(ai,ld) = u_il^ad,
!     which means that u_aild(bj,kc)^T = u_aild(kc,bj) = u_kj^cb = u_jk^bc.
!
      call mem%dealloc(Z_aikc,  wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     :: Calculate the D2.1 term of omega ::
!
!     Form g_aikc = g_aikc
!
      call mem%alloc(g_aikc, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call wf%get_voov(g_aikc)
!
!     Calculate the D2.1 term, sum_ck u_jk^bc g_aikc = sum_ck g_aikc(ai,kc) u_aild(bj,kc)
!
      call dgemm('N','T',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  g_aikc,            &
                  (wf%n_o)*(wf%n_v), &
                  u_aild,            &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  omega2_aibj,       &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(u_aild,  wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(g_aikc, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
!     Add the D2.1 term to the omega vector
!
      call symmetrize_and_add_to_packed(omega2, omega2_aibj, (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(omega2_aibj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call ccsd_d2_timer%freeze()
      call ccsd_d2_timer%switch_off()
!
   end subroutine omega_ccsd_d2_ccsd
!
!
   module subroutine omega_ccsd_e2_ccsd(wf, omega2)
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
      real(dp), dimension((wf%n_v)*(wf%n_o)*((wf%n_v)*(wf%n_o)+1)/2), intent(inout):: omega2
!
!     Vectors for E2.1 term
!
      real(dp), dimension(:,:), allocatable :: omega2_bjai ! For storing the E2.1 term temporarily
      real(dp), dimension(:,:,:,:), allocatable :: g_ldkc      ! g_ldkc
      real(dp), dimension(:,:,:,:), allocatable :: u_bldk      ! u_kl^bd
      real(dp), dimension(:,:,:,:), allocatable :: t_bkdl      ! t_kl^bd
      real(dp), dimension(:,:,:,:), allocatable :: u_bkdl      ! u_kl^bd
      real(dp), dimension(:,:), allocatable :: X_b_c           ! An intermediate, see below for definition
      real(dp), dimension(:,:,:,:), allocatable :: t_cjai      ! t_ij^ac
!
!     Vectors for E2.2 term
!
      real(dp), dimension(:,:), allocatable :: omega2_aibj ! For storing the E2.2 term temporarily
      real(dp), dimension(:,:), allocatable :: Y_k_j        ! An intermediate, see below for definition
!
      type(timings) :: ccsd_e2_timer 
!
      call ccsd_e2_timer%init('omega ccsd e2')
      call ccsd_e2_timer%start()
!
!     :: Calculate the E2.1 term of omega ::
!
!     Form u_bldk = u_kl^bd = u_bkdl
!                  = 2 * t_kl^bd - t_lk^bd = 2 * t_bkdl(bk,dl) - t_bkdl(bl, dk)
!
      call mem%alloc(t_bkdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_bkdl, (wf%n_v)*(wf%n_o))
!
      call mem%alloc(u_bkdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      u_bkdl = zero
!
      call add_1432_to_1234(-one, t_bkdl, u_bkdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call daxpy((wf%n_v)**2*(wf%n_o)**2, two, t_bkdl, 1, u_bkdl, 1)
!
      call mem%dealloc(t_bkdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(u_bldk,  wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_1432(u_bkdl, u_bldk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(u_bkdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Form g_ldkc = g_ldkc
!
      call mem%alloc(g_ldkc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call wf%get_ovov(g_ldkc)
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
      call mem%alloc(t_cjai, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_cjai, (wf%n_v)*(wf%n_o))
!
!     Form the E2.1 term
!
      call mem%alloc(omega2_bjai, (wf%n_o)*(wf%n_v), (wf%n_v)*(wf%n_o))
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_v,               &
                  one,                  &
                  X_b_c,                &
                  wf%n_v,               &
                  t_cjai,               & ! t_c_jai
                  wf%n_v,               &
                  zero,                 &
                  omega2_bjai,          & ! omega2_b_jai
                  wf%n_v)
!
!     Add the E2.1 term to the omega vector and deallocations
!
      call symmetrize_and_add_to_packed(omega2, omega2_bjai, (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(omega2_bjai, (wf%n_o)*(wf%n_v), (wf%n_v)*(wf%n_o))
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
      call mem%alloc(omega2_aibj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call dgemm('N','N',              &
                 (wf%n_o)*(wf%n_v)**2, &
                 wf%n_o,               &
                 wf%n_o,               &
                 -one,                 &
                 t_cjai,               & ! t_aib_k
                 (wf%n_o)*(wf%n_v)**2, &
                 Y_k_j,                &
                 wf%n_o,               &
                 zero,                 &
                 omega2_aibj,          & ! omega2_aib_j
                 (wf%n_o)*(wf%n_v)**2)
!
!     Deallocate Y_k_j and the amplitudes
!
      call mem%dealloc(Y_k_j, wf%n_o, wf%n_o)
      call mem%dealloc(t_cjai, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Add the E2.2 term to the omega vector
!
      call symmetrize_and_add_to_packed(omega2, omega2_aibj, (wf%n_o)*(wf%n_v))
!
!     Deallocate the E2.2 term
!
      call mem%dealloc(omega2_aibj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call ccsd_e2_timer%freeze()
      call ccsd_e2_timer%switch_off()
!
   end subroutine omega_ccsd_e2_ccsd
!
!
end submodule omega_ccsd 
