!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2022 the authors of eT
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
submodule (mlccsd_class) omega_mlccsd
!
!!
!!    Omega submodule (MLCCSD)
!!
!!    Routines to construct
!!
!!    Omega =  < mu | exp(-T) H exp(T) | R >
!!
!!    Based on omega_ccsd.F90 written by Sarai D. Folkestad
!!    and Eirik F. Kjønstad
!!
!!    Modified for MLCCSD by restricting indices to the active space.
!!
!
   implicit none
!
contains
!
!
   module subroutine construct_omega_mlccsd(wf, omega)
!!
!!    Construct omega (MLCCSD)
!!    Written by Sarai D. Folkestad, 2017-2019
!!
!!    Directs the construction of the projection vector < mu | exp(-T) H exp(T) | R >
!!    for the current amplitudes of the object wf
!!
      use array_utilities, only : scale_diagonal, zero_array
      use reordering, only: squareup, packin
!
      implicit none
!
      class(mlccsd), intent(inout) :: wf
!
      real(dp), dimension(wf%n_gs_amplitudes), intent(out) :: omega
!
      real(dp), dimension(:,:,:,:), allocatable :: omega_aibj
      real(dp), dimension(:,:,:,:), allocatable :: x_aibj
!
      integer :: n_a_o, n_a_v
!
      type(timings) :: timer
!
      timer = timings('Construct MLCCSD Omega', pl='n')
      call timer%turn_on()
!
!     Set the omega vector to zero
!
      call zero_array(omega, wf%n_gs_amplitudes)
!
      n_a_o = wf%n_cc2_o + wf%n_ccsd_o
      n_a_v = wf%n_cc2_v + wf%n_ccsd_v
!
      call wf%ccs%construct_omega(omega(1 : wf%n_t1))
!
!     Construct singles contributions
!
      call wf%construct_u_aibj()
!
      call wf%omega_cc2_a1(omega(1 : wf%n_t1), n_a_o, n_a_v, 1, 1, n_a_o, n_a_v)
      call wf%omega_cc2_b1(omega(1 : wf%n_t1), n_a_o, n_a_v, 1, 1, n_a_o, n_a_v)
      call wf%omega_cc2_c1(omega(1 : wf%n_t1), n_a_o, n_a_v, 1, 1)
!
      call mem%alloc(x_aibj, n_a_v, n_a_o, n_a_v, n_a_o)
      call squareup(wf%x2, x_aibj, n_a_o*n_a_v)
!
!     Construct doubles contributions
!
      call mem%alloc(omega_aibj, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
      call zero_array(omega_aibj, (wf%n_ccsd_v**2)*(wf%n_ccsd_o**2))
!
      call wf%omega_ccsd_a2(omega_aibj)
      call wf%omega_ccsd_b2(omega_aibj, x_aibj)
      call wf%omega_ccsd_c2(omega_aibj, x_aibj)
      call wf%omega_ccsd_d2(omega_aibj, x_aibj)
      call wf%omega_ccsd_e2(omega_aibj, x_aibj)
      call wf%omega_ccsd_f2(omega_aibj, x_aibj)
!
      call scale_diagonal(half, omega_aibj, wf%n_ccsd_v*wf%n_ccsd_o)
!
      call packin(omega(wf%n_t1 + 1:wf%n_gs_amplitudes), omega_aibj, wf%n_ccsd_v*wf%n_ccsd_o)
!
      call mem%dealloc(omega_aibj, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
      call mem%dealloc(x_aibj, n_a_v, n_a_o, n_a_v, n_a_o)
!
      call timer%turn_off()
!
   end subroutine construct_omega_mlccsd
!
!
   module subroutine omega_ccsd_a2_mlccsd(wf, omega_aibj)
!!
!!    Omega A2 term
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!    Modified for MLCCSD by Sarai D. Folkestad 2017-2019
!!
!!    A2 = g_aibj
!!
!!    INDEX RESTRICTIONS
!!
!!    a, b, i, j are CCSD indices
!!
      implicit none
!
      class(mlccsd) :: wf
!
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o), &
               intent(inout):: omega_aibj
!
!     Integrals
!
      real(dp), dimension(:,:,:,:), allocatable :: g_aibj
!
      type(timings) :: timer
!
      integer :: n_a_o, n_a_v
!
      timer = timings('Omega MLCCSD a2', pl='v')
!
      call timer%turn_on()
!
      n_a_o = wf%n_cc2_o + wf%n_ccsd_o
      n_a_v = wf%n_cc2_v + wf%n_ccsd_v
!
      call mem%alloc(g_aibj, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
      call wf%eri_t1%get('vovo', g_aibj, 1, wf%n_ccsd_v, 1, wf%n_ccsd_o, &
                                             1, wf%n_ccsd_v, 1, wf%n_ccsd_o)
!
      call daxpy((wf%n_ccsd_v**2)*(wf%n_ccsd_o**2), one, g_aibj, 1, omega_aibj, 1)

      call mem%dealloc(g_aibj, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
      call timer%turn_off()
!
   end subroutine omega_ccsd_a2_mlccsd
!
!
   module subroutine omega_ccsd_b2_mlccsd(wf, omega_aibj, x_aibj)
!!
!!    Omega B2 term
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!    Modified for MLCCSD by Sarai D. Folkestad 2017-2019
!!
!!    B2 =  sum_(cd)g_acbd * x_cidj = B2.1 + B.2.2
!!
!!    Structure: Batching over both a and b for B2.2.
!!                t^+_ci_dj = t_cidj + t_di_cj
!!                t^-_ci_dj = t_cidj - t_di_cj
!!                g^+_ac_bd = g_acbd + g_bc_ad
!!                g^-_ac_bd = g_acbd - g_bc_ad
!!
!!       omega_B2.2_ai_bj = 1/4*(g^+_ac_bd*t^+_ci_dj + g^-_ac_bd*t^-_ci_dj) = omega_A2.2_bj_ai
!!       omega_B2.2_aj_bi = 1/4*(g^+_ac_bd*t^+_ci_dj - g^-_ac_bd*t^-_ci_dj) = omega_A2.2_bi_aj
!!
!!
!!    INDEX RESTRICTIONS
!!
!!    a, b, i, j are CCSD indices
!!
!!    c, d, are CC2 indices
!!
!
      implicit none
!
      class(mlccsd) :: wf
!
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o), &
               intent(inout):: omega_aibj
!
      real(dp), dimension(wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o, &
                          wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o), intent(in) :: x_aibj
!
!
!     Integrals
!
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
      integer :: ij
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
      type(timings) :: timer, integral_timer
!
      integer :: n_a_o, n_a_v
!
      timer = timings('Omega MLCCSD B2', pl='v')
      call timer%turn_on()
      integral_timer = timings('Omega MLCCSD B2 g_abcd', pl='v')
!
      n_a_o = wf%n_cc2_o + wf%n_ccsd_o
      n_a_v = wf%n_cc2_v + wf%n_ccsd_v
!
!     Some helpful integers
!
      n_v_packed = n_a_v*(n_a_v+1)/2
      n_o_packed = wf%n_ccsd_o*(wf%n_ccsd_o+1)/2
!
      req0 = 2*(n_v_packed)*(n_o_packed)
!
      req1_a = 0
      req1_b = 0
!
      rec2 = 2*(n_a_v)**2 + 2*(n_o_packed)
!
!     Initialize batching variables
!
      batch_a = batching_index(wf%n_ccsd_v)
      batch_b = batching_index(wf%n_ccsd_v)
!
      call mem%batch_setup(batch_a, batch_b, req0, req1_a, req1_b, rec2, &
                           tag='omega_ccsd_b2_mlccsd')
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
            call mem%alloc(g_acbd, batch_a%length, n_a_v, batch_b%length, n_a_v)
!
            call integral_timer%turn_on()
!
            call wf%eri_t1%get('vvvv', g_acbd,                      &
                                   batch_a%first, batch_a%get_last(),   &
                                   1, n_a_v,                            &
                                   batch_b%first, batch_b%get_last(),   &
                                   1, n_a_v)
!
            call integral_timer%freeze()
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
!$omp parallel do private(a, b, c, d, ab, cd, diag_factor)
               do c = 1, n_a_v
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
!
                           g_m_abcd(ab, cd) = diag_factor*(g_acbd(a, d, b, c) - g_acbd(a, c, b, d))
!                          a and b and c and d switched
!
                        enddo
                     enddo
                  enddo
               enddo
!$omp end parallel do
!
!$omp parallel do private(i, j, c, d, ij, cd)
              do i = 1, wf%n_ccsd_o
                 do j = 1, i
!
                     ij = (i*(i-3)/2) + i + j
!
                     do c = 1, n_a_v
                        do d = 1, c
!
                           cd = (c*(c-3)/2) + c + d
!
                           t_p_cdij(cd, ij) = x_aibj(c, i, d, j) &
                                            + x_aibj(d, i, c, j)
                           t_m_cdij(cd, ij) = x_aibj(c, i, d, j) &
                                            - x_aibj(d, i, c, j)
!
                       enddo
                    enddo
                 enddo
              enddo
!$omp end parallel do
!
!              Dellocate g_acbd
!
               call mem%dealloc(g_acbd, batch_a%length, n_a_v, batch_b%length, n_a_v)
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
!$omp parallel do private(i, j, a, b, ij, ab)
               do i = 1, wf%n_ccsd_o
                  do j = 1, i
!
                     ij = (i*(i-3)/2) + i + j
!
                     do a = 1, batch_a%length
                        do b = 1, a
!
                           ab = (a*(a-3)/2) + a + b
!
!                          Reorder into omega2_aibj
!
                           omega_aibj(a + batch_a%first - 1, i, b + batch_b%first - 1, j) =  &
                                              omega_aibj(a + batch_a%first - 1, i,  &
                                                         b + batch_b%first - 1, j)  &
                                             + omega2_p_abij(ab, ij) + omega2_m_abij(ab, ij)
!
                           if (a .ne. b .and. i .ne. j) then
!
                              omega_aibj(b + batch_b%first - 1, i, a + batch_a%first - 1, j) = &
                                            omega_aibj(b + batch_b%first - 1, i, &
                                                       a + batch_a%first - 1, j) &
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
!$omp parallel do private(a, b, c, d, ab, cd, diag_factor)
               do c = 1, n_a_v
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
!$omp parallel do schedule(static) private(c, d, i, j, cd, ij)
               do i = 1, wf%n_ccsd_o
                  do j = 1, i
!
                     ij = (i*(i-3)/2) + i + j
!
                     do c = 1, n_a_v
                        do d = 1, c
!
                           cd = (c*(c-3)/2) + c + d
!
                           t_p_cdij(cd, ij) = x_aibj(c, i, d, j) &
                                            + x_aibj(d, i, c, j)
!
                           t_m_cdij(cd, ij) = x_aibj(c, i,  d, j) &
                                            - x_aibj(d, i,  c, j)
!
                        enddo
                     enddo
                  enddo
               enddo
!$omp end parallel do
!
!              Dellocate g_acbd
!
               call mem%dealloc(g_acbd, batch_a%length, n_a_v, batch_b%length, n_a_v)
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
!$omp parallel do private(i, j, a, b, ij, ab)
               do i = 1, wf%n_ccsd_o
                  do j = 1, i
!
                     ij = (i*(i-3)/2) + i + j
!
                     do a = 1, batch_a%length
!
                        do b = 1, batch_b%length
!
                           ab = batch_a%length*(b - 1) + a
!
!                          Reorder into omega2_aibj
!
                           omega_aibj(a + batch_a%first - 1, i, b + batch_b%first - 1, j) =  &
                                 omega_aibj(a + batch_a%first - 1, i, b + batch_b%first - 1, j)&
                                 + omega2_p_abij(ab, ij) + omega2_m_abij(ab, ij)
!
                           if (i .ne. j) then
!
                              omega_aibj(b + batch_b%first - 1, i, a + batch_a%first - 1, j) = &
                                    omega_aibj(b + batch_b%first - 1, i, a + batch_a%first - 1, j) &
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
      call mem%batch_finalize()
!
      call timer%turn_off()
      call integral_timer%turn_off()
!
   end subroutine omega_ccsd_b2_mlccsd
!
!
   module subroutine omega_ccsd_c2_mlccsd(wf, omega_aibj, x_aibj)
!!
!!    Omega C2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!    Modified for MLCCSD by Sarai D. Folkestad 2017-2019
!!
!!    C2 = sum_(kl) x_ak_bl*(g_kilj + sum_(cd) x_cidj * g_kcld)
!!
!!    INDEX RESTRICTIONS:
!!
!!    a, i, b, j are CCSD orbitals
!!
!!    k, l, c, d are CC2 + CCSD orbitals
!!
      use reordering, only: sort_1234_to_1324, add_1324_to_1234
!
      implicit none
!
      class(mlccsd) :: wf
!
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o), &
               intent(inout):: omega_aibj
!
      real(dp), dimension(wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o, &
                          wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o), intent(in) :: x_aibj
!
!
!     Integrals
!
      real(dp), dimension(:,:,:,:), allocatable :: g_kcld
      real(dp), dimension(:,:,:,:), allocatable :: g_klcd
      real(dp), dimension(:,:,:,:), allocatable :: g_klij
      real(dp), dimension(:,:,:,:), allocatable :: g_kilj
!
!     Reordered X2 amplitudes
!
      real(dp), dimension(:,:,:,:), allocatable :: x_cdij
      real(dp), dimension(:,:,:,:), allocatable :: x_abkl
!
!     Reordered omega
!
      real(dp), dimension(:,:,:,:), allocatable :: omega_abij
!
      type(timings) :: timer
!
      integer :: c, d, i, j, a, b, k, l
!
      integer :: n_a_o, n_a_v
!
      timer = timings('Omega MLCCSD C2', pl='v')
      call timer%turn_on()
!
      n_a_o = wf%n_cc2_o + wf%n_ccsd_o
      n_a_v = wf%n_cc2_v + wf%n_ccsd_v
!
!     Allocate and construct g_kilj
!
      call mem%alloc(g_kilj, n_a_o, wf%n_ccsd_o, n_a_o, wf%n_ccsd_o)
!
      call wf%eri_t1%get('oooo', g_kilj, &
                             1, n_a_o,       &
                             1, wf%n_ccsd_o, &
                             1, n_a_o,       &
                             1, wf%n_ccsd_o)
!
      call mem%alloc(g_klij, n_a_o, n_a_o, wf%n_ccsd_o, wf%n_ccsd_o)
!
      call sort_1234_to_1324(g_kilj, g_klij, n_a_o, wf%n_ccsd_o, n_a_o, wf%n_ccsd_o)
!
      call mem%dealloc(g_kilj, n_a_o, wf%n_ccsd_o, n_a_o, wf%n_ccsd_o)
!
!     Allocate and construct g_kcld
!
      call mem%alloc(g_kcld, n_a_o, n_a_v, n_a_o, n_a_v)
!
      call wf%eri_t1%get('ovov', g_kcld, 1, n_a_o, 1, n_a_v, 1, n_a_o, 1, n_a_v)
!
!     Reorder g_kcld as g_klcd
!
      call mem%alloc(g_klcd, n_a_o, n_a_o, n_a_v, n_a_v)
!
      call sort_1234_to_1324(g_kcld, g_klcd, n_a_o, n_a_v, n_a_o, n_a_v)
!
      call mem%dealloc(g_kcld, n_a_o, n_a_v, n_a_o, n_a_v)
!
!     Reorder x_cidj as x_cdij
!
      call mem%alloc(x_cdij, n_a_v, n_a_v, wf%n_ccsd_o, wf%n_ccsd_o)
!
!$omp parallel do private (j, i, c, d) collapse(2)
      do j = 1, wf%n_ccsd_o
         do i = 1, wf%n_ccsd_o
            do c = 1, n_a_v
               do d = 1, n_a_v
!
                  x_cdij(c, d, i, j) = x_aibj(c, i, d, j)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call dgemm('N','N',           &
                  (n_a_o)**2,       &
                  (wf%n_ccsd_o)**2, &
                  (n_a_v)**2,       &
                  one,              &
                  g_klcd,           &
                  (n_a_o)**2,       &
                  x_cdij,           &
                  (n_a_v)**2,       &
                  one,              &
                  g_klij,           &
                  (n_a_o)**2)
!
!     Deallocate t_cdij and g_klcd
!
      call mem%dealloc(g_klcd, n_a_o, n_a_o, n_a_v, n_a_v)!
      call mem%dealloc(x_cdij, n_a_v, n_a_v, wf%n_ccsd_o, wf%n_ccsd_o)
!
!     omega_abij = sum_(kl) x_ab_kl*X_kl_ij
!
      call mem%alloc(x_abkl, wf%n_ccsd_v, wf%n_ccsd_v, n_a_o, n_a_o)
!
!$omp parallel do private (a, b, k, l) collapse(2)
      do l = 1, n_a_o
         do k = 1, n_a_o
            do b = 1, wf%n_ccsd_v
              do a = 1, wf%n_ccsd_v
!
                  x_abkl(a, b, k, l) = x_aibj(a, k, b, l)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%alloc(omega_abij, wf%n_ccsd_v, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_o)
!
      call dgemm('N','N',           &
                  (wf%n_ccsd_v)**2, &
                  (wf%n_ccsd_o)**2, &
                  (n_a_o)**2,       &
                  one,              &
                  x_abkl,           &
                  (wf%n_ccsd_v)**2, &
                  g_klij,           &
                  (n_a_o)**2,       &
                  zero,             &
                  omega_abij,       &
                  (wf%n_ccsd_v)**2)
!
      call mem%dealloc(g_klij, n_a_o, n_a_o, wf%n_ccsd_o, wf%n_ccsd_o)
      call mem%dealloc(x_abkl, wf%n_ccsd_v, wf%n_ccsd_v, n_a_o, n_a_o)
!
!     Reorder into omega_aibj and symmetrize
!
      call add_1324_to_1234(one, omega_abij, omega_aibj, wf%n_ccsd_v, &
                           wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
      call mem%dealloc(omega_abij, wf%n_ccsd_v, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_o)
!
      call timer%turn_off()
!
   end subroutine omega_ccsd_c2_mlccsd
!
!
   module subroutine omega_ccsd_d2_mlccsd(wf, omega_aibj, x_aibj)
!!
!!    Omega D2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!    Modified for MLCCSD by Sarai D. Folkestad 2017-2019
!!
!!    D2 = -1/2 * sum_(ck) x_bk_cj*(g_kiac - 1/2 sum_(dl)x_al_di * g_kdlc)
!!                    - sum_(ck) x_bk_ci*(g_kj_ac - sum_(dl)x_al_dj * g_kdlc)
!!                    - 1/2 * sum_ck u_jk^bc g_acki
!!
!!    INDEX RESTRICTIONS:
!!
!!    a, i, b, j are CCSD indices
!!
!!    c, k, d, l are CCSD + CC2 indices
!!
      use reordering, only: sort_1234_to_2341, sort_1234_to_3241, symmetric_sum
      use array_utilities, only: zero_array
!
      implicit none

      class(mlccsd) :: wf
!
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o), &
                  intent(inout):: omega_aibj
!
      real(dp), dimension(wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o, &
                          wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o), intent(in) :: x_aibj
!
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
      real(dp), dimension(:,:,:,:), allocatable :: x_aidl
      real(dp), dimension(:,:,:,:), allocatable :: x_ckbj
!
!     Intermediates for matrix multiplication
!
      real(dp), dimension(:,:,:,:), allocatable :: Y_aick
      real(dp), dimension(:,:,:,:), allocatable :: I_aibj
!
!     Reordered U2 amplitudes
!
      real(dp), dimension(:,:,:,:), allocatable :: u_ckbj
      real(dp), dimension(:,:,:,:), allocatable :: omega2_aibj ! Holds term temporarily
      real(dp), dimension(:,:,:,:), allocatable :: omega_a_batch
!
!     Indices
!
      integer :: a, b, c, d
      integer :: i, j, k, l
!
!     Batching and memory handling
!
      integer :: req0, req1
      integer :: current_a_batch = 0
!
      type(batching_index) :: batch_a
!
      type(timings) :: timer
!
      integer :: n_a_o, n_a_v
!
      timer = timings('Omega MLCCSD D2', pl='v')
      call timer%turn_on()
!
      n_a_o = wf%n_cc2_o + wf%n_ccsd_o
      n_a_v = wf%n_cc2_v + wf%n_ccsd_v
!
!     Sort x_aldi = x_li^ad as x_aidl
!
      call mem%alloc(x_aidl, wf%n_ccsd_v, wf%n_ccsd_o, n_a_v, n_a_o)
!
!$omp parallel do private (a, i, d, l) collapse(2)
      do l = 1, n_a_o
         do d = 1, n_a_v
            do i = 1, wf%n_ccsd_o
               do a = 1, wf%n_ccsd_v
!
                  x_aidl(a, i, d, l) =x_aibj(a, l, d, i)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%alloc(g_kdlc, n_a_o, n_a_v, n_a_o, n_a_v)
!
      call wf%eri_t1%get('ovov', g_kdlc, 1, n_a_o, 1, n_a_v, 1, n_a_o, 1, n_a_v)
!
!     Sort g_kdlc to g_dlck
!
      call mem%alloc(g_dlck, n_a_v, n_a_o, n_a_v, n_a_o)
!
      call sort_1234_to_2341(g_kdlc, g_dlck, n_a_o, n_a_v, n_a_o, n_a_v)
!
      call mem%dealloc(g_kdlc, n_a_o, n_a_v, n_a_o, n_a_v)
!
!     -1/2*sum_(dl) x_aidl*g_dlck = Y_aick
!
      call mem%alloc(Y_aick, wf%n_ccsd_v, wf%n_ccsd_o, n_a_v, n_a_o)
!
      call dgemm('N','N',                       &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),  &
                  (n_a_v)*(n_a_o),              &
                  (n_a_v)*(n_a_o),              &
                  -half,                        &
                  x_aidl,                       &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),  &
                  g_dlck,                       &
                  (n_a_v)*(n_a_o),              &
                  zero,                         &
                  Y_aick,                       &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v))
!
      call mem%dealloc(g_dlck, n_a_v, n_a_o, n_a_v, n_a_o)
      call mem%dealloc(x_aidl, wf%n_ccsd_v, wf%n_ccsd_o, n_a_v, n_a_o)
!
!     Form u_ckbj = u_jk^bc  = 2 * x_jk^bc - x_kj^bc
!
      call mem%alloc(u_ckbj, n_a_v, n_a_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
!$omp parallel do private (j, b, k, c) collapse(2)
      do j = 1, wf%n_ccsd_o
         do b = 1, wf%n_ccsd_v
            do k = 1, n_a_o
               do c = 1, n_a_v
!
                  u_ckbj(c, k, b, j) = &
                              two*x_aibj(c, k, b, j) &
                              - one*x_aibj(b, k, c, j)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
!     Allocate a holder for - 1/2 * sum_ck u_jk^bc g_acki,
!     constructed in batches over the a index below
!
      call mem%alloc(omega2_aibj,  wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
      call zero_array(omega2_aibj, (wf%n_ccsd_o**2)*(wf%n_ccsd_v**2))
!
!     Constructing g_kiac
!
!     Prepare for batching
!
      req0 = (wf%n_ccsd_o)*(n_a_o)*wf%eri_t1%n_J
!
      req1 = (n_a_v)*wf%eri_t1%n_J + (wf%n_ccsd_o)*(n_a_o)*(n_a_v)
!
      batch_a = batching_index(wf%n_ccsd_v)
!
      call mem%batch_setup(batch_a, req0, req1, tag='omega_ccsd_d2_mlccsd')
!
!     Loop over the number of a batches
!
      do current_a_batch = 1, batch_a%num_batches
!
         call batch_a%determine_limits(current_a_batch)
!
!        Allocate and construct g_kiac
!
         call mem%alloc(g_kiac, n_a_o, wf%n_ccsd_o, batch_a%length, n_a_v)
!
         call wf%eri_t1%get('oovv', g_kiac, 1, n_a_o, 1, wf%n_ccsd_o, &
                                                batch_a%first, batch_a%get_last(), 1, n_a_v)
!
!        X_aick = X_aick + g_kiac
!
!$omp parallel do private(a, i, k, c) collapse(2)
         do k = 1, n_a_o
            do c = 1, n_a_v
               do i = 1, wf%n_ccsd_o
                  do a = 1, batch_a%length
!
                     Y_aick(a+batch_a%first-1, i, c, k) = Y_aick(a+batch_a%first-1, i, c, k) &
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
         call mem%alloc(g_aick, batch_a%length, wf%n_ccsd_o, n_a_v, n_a_o)
!
         call sort_1234_to_3241(g_kiac, g_aick, n_a_o, wf%n_ccsd_o, batch_a%length, n_a_v)
!
         call mem%dealloc(g_kiac, n_a_o, wf%n_ccsd_o, batch_a%length, n_a_v)
!
!        - 1/2 * sum_ck u_jk^bc g_acki = -1/2 * sum_ck g_aick u_ckbj
!
         call mem%alloc(omega_a_batch, batch_a%length, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
         call dgemm('N','N',                       &
                  (wf%n_ccsd_o)*(batch_a%length),  &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),     &
                  (n_a_o)*(n_a_v),                 &
                  -half,                           &
                  g_aick,                          &
                  (wf%n_ccsd_o)*(batch_a%length),  &
                  u_ckbj,                          &
                  (n_a_o)*(n_a_v),                 &
                  zero,                            &
                  omega_a_batch,                   &
                  (wf%n_ccsd_o)*batch_a%length)
!
         call mem%dealloc(g_aick, batch_a%length, wf%n_ccsd_o, n_a_v, n_a_o)
!
!$omp parallel do private(a, i, b, j) collapse(2)
         do j = 1, wf%n_ccsd_o
            do b = 1, wf%n_ccsd_v
               do i = 1, wf%n_ccsd_o
                  do a = 1, batch_a%length
!
                     omega2_aibj(a + batch_a%first - 1, i, b, j) = &
                                    omega2_aibj(a + batch_a%first - 1, i, b, j) &
                                    + omega_a_batch(a, i, b, j)
!
                  enddo
               enddo
            enddo
         enddo
!$omp end parallel do
!
         call mem%dealloc(omega_a_batch, batch_a%length, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
      enddo ! End of batching
!
      call mem%batch_finalize()
!
      call mem%dealloc(u_ckbj, n_a_v, n_a_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
!      =+ omega2_aibj(ai,bj) + omega2_bj_ai(bj,ai)
!
      call symmetric_sum(omega2_aibj, (wf%n_ccsd_o)*(wf%n_ccsd_v))   ! symmetrize
!
      call daxpy(((wf%n_ccsd_o)*(wf%n_ccsd_v))**2, one, omega2_aibj, 1, omega_aibj, 1)
!
      call mem%dealloc(omega2_aibj, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
!     Reorder x_bkcj = x_kj^bc as x_ckbj
!
      call mem%alloc(x_ckbj,  n_a_v, n_a_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
!$omp parallel do private (b, j, k, c) collapse(2)
      do j = 1, wf%n_ccsd_o
         do b = 1, wf%n_ccsd_v
            do k = 1, n_a_o
               do c = 1, n_a_v
!
                  x_ckbj(c, k, b, j) = x_aibj(b, k, c, j)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
!     Form intermediate Y_aibj = - sum_(ck) X_aick x_ckbj
!
      call mem%alloc(I_aibj,  wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
      call dgemm('N','N',                       &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),  &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),  &
                  (n_a_o)*(n_a_v),              &
                  -one,                         &
                  Y_aick,                       &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),  &
                  x_ckbj,                       &
                  (n_a_o)*(n_a_v),              &
                  zero,                         &
                  I_aibj,                       &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v))
!
      call mem%dealloc(Y_aick,  wf%n_ccsd_v, wf%n_ccsd_o, n_a_v, n_a_o)
      call mem%dealloc(x_ckbj,  n_a_v, n_a_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
!     Omega_aibj,1 = P_ai_bj ( 1/2*I_aibj + I_aj_bi )
!
      call mem%alloc(omega2_aibj, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
!$omp parallel do private(i, a, j, b) collapse(2)
         do i = 1, wf%n_ccsd_o
            do a = 1, wf%n_ccsd_v
               do j = 1, wf%n_ccsd_o
                  do b = 1, wf%n_ccsd_v
!
                  omega2_aibj(a, i, b, j) = half*I_aibj(a, i, b, j) + I_aibj(a, j, b, i)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(I_aibj,  wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
      call symmetric_sum(omega2_aibj, wf%n_ccsd_v*wf%n_ccsd_o)
!
      call daxpy(((wf%n_ccsd_o)*(wf%n_ccsd_v))**2, one, omega2_aibj, 1, omega_aibj, 1)
!
      call mem%dealloc(omega2_aibj, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
      call timer%turn_off()
!
   end subroutine omega_ccsd_d2_mlccsd
!
!
   module subroutine omega_ccsd_e2_mlccsd(wf, omega_aibj, x_aibj)
!!
!!    Omega E2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!    Modified for MLCCSD by Sarai D. Folkestad 2017-2019
!!
!!    Calculates the D2 term,
!!
!!      E2: sum_ck u_jk^bc g_aikc
!!        + 1/4 * sum_ck u_jk^bc sum_dl L_ldkc u_il^ad,
!!
!!    where
!!
!!        u_jk^bc = 2 * t_jk^bc - t_kj^bc,
!!        L_ldkc  = 2 * g_ldkc  - g_lckd.
!!
!!    The first and second terms are referred to as E2.1 and E2.2.
!!    All terms are added to the omega vector of the wavefunction object wf.
!!
!!    The routine adds the terms in the following order: E2.2, E2.1
!!
!!
!!    INDEX RESTRICTIONS:
!!
!!    a, i, b, j are CCSD indices
!!
!!    c, d, k, l are CC2 + CCSD indices
!!
      use reordering, only: add_1432_to_1234, symmetric_sum
      use array_utilities, only: copy_and_scale
!
      implicit none
!
      class(mlccsd) :: wf
!
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o), &
                  intent(inout):: omega_aibj
!
      real(dp), dimension(wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o, &
                          wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o), intent(in) :: x_aibj
!
!
      real(dp), dimension(:,:,:,:), allocatable :: omega2_aibj ! For storing E2.2 & E2.1
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ldkc
      real(dp), dimension(:,:,:,:), allocatable :: L_ldkc
      real(dp), dimension(:,:,:,:), allocatable :: u_aild
      real(dp), dimension(:,:,:,:), allocatable :: Z_aikc
      real(dp), dimension(:,:,:), allocatable   :: L_J_kc, L_J_ai ! Cholesky vectors, term 1
      real(dp), dimension(:,:,:), allocatable   :: X_bj_J         ! Intermediate, term 1
!
      type(timings) :: timer
!
      integer :: a, i, d, l
!
      integer :: n_a_o, n_a_v
!
      timer = timings('Omega MLCCSD E2', pl='v')
      call timer%turn_on()
!
      n_a_o = wf%n_cc2_o + wf%n_ccsd_o
      n_a_v = wf%n_cc2_v + wf%n_ccsd_v
!
!     Form L_ld_kc = L_ldkc = 2*g_ldkc(ld,kc) - g_ldkc(lc,kd)
!
      call mem%alloc(g_ldkc, n_a_o, n_a_v, n_a_o, n_a_v)
      call wf%eri_t1%get('ovov', g_ldkc, 1, n_a_o, 1, n_a_v, 1, n_a_o, 1, n_a_v)
!
      call mem%alloc(L_ldkc, n_a_o, n_a_v, n_a_o, n_a_v)
!
      call copy_and_scale(two, g_ldkc, L_ldkc, (n_a_o**2)*(n_a_v**2))
      call add_1432_to_1234(-one, g_ldkc, L_ldkc, n_a_o, n_a_v, n_a_o, n_a_v)
!
      call mem%dealloc(g_ldkc, n_a_o, n_a_v, n_a_o, n_a_v)
!
!     Form u_aild = u_il^ad = 2 * x_il^ad - x_li^ad ordered as u_aild
!
      call mem%alloc(u_aild,  wf%n_ccsd_v, wf%n_ccsd_o, n_a_o, n_a_v)
!
!$omp parallel do private (a, i, d, l)
      do d = 1, n_a_v
         do l = 1, n_a_o
            do i = 1, wf%n_ccsd_o
               do a = 1, wf%n_ccsd_v
!
                  u_aild(a, i, l, d) =&
                              two*x_aibj(a, i, d, l)&
                                 - x_aibj(a, l, d, i)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
!     Form the intermediate Z_aikc = sum_dl u_aild L_ldkc and set it to zero
!
      call mem%alloc(Z_aikc, wf%n_ccsd_v, wf%n_ccsd_o, n_a_o, n_a_v)
!
      call dgemm('N','N',                       &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),  &
                  (n_a_o)*(n_a_v),              &
                  (n_a_o)*(n_a_v),              &
                  one,                          &
                  u_aild,                       &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),  &
                  L_ldkc,                       &
                  (n_a_o)*(n_a_v),              &
                  zero,                         &
                  Z_aikc,                       &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v))
!
      call mem%dealloc(L_ldkc, n_a_o, n_a_v, n_a_o, n_a_v)
!
!      1/4 sum_kc Z_aikc u_kc_bj = 1/4 sum_kc Z_aikc(ai,kc) u_aild(bj,kc)
!
      call mem%alloc(omega2_aibj, (wf%n_ccsd_v), (wf%n_ccsd_o), (wf%n_ccsd_v), (wf%n_ccsd_o))
!
      call dgemm('N','T',                       &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),  &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),  &
                  (n_a_o)*(n_a_v),              &
                  one/four,                     &
                  Z_aikc,                       &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),  &
                  u_aild,                       & ! u_bjkc
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),  &
                  zero,                         &
                  omega2_aibj,                  &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v))
!
      call mem%dealloc(Z_aikc,  wf%n_ccsd_v, wf%n_ccsd_o, n_a_o, n_a_v)
!
!     u_jk^bc g_aikc = L_ai^J (u_bjkc L_kc^J) = L_ai^J X_bj^J
!
!     X_bj_J = u_bjkc L_J_kc
!
      call mem%alloc(L_J_kc, wf%eri_t1%n_J, n_a_o, n_a_v)
      call wf%L_t1%get(L_J_kc, 1, n_a_o, wf%n_o + 1, wf%n_o + n_a_v)
!
      call mem%alloc(X_bj_J, wf%n_ccsd_v, wf%n_ccsd_o, wf%eri_t1%n_J)
!
      call dgemm('N', 'T',                     &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v), &
                  wf%eri_t1%n_J,                  &
                  (n_a_o)*(n_a_v),             &
                  one,                         &
                  u_aild,                      & ! u_bj,kc
                  (wf%n_ccsd_o)*(wf%n_ccsd_v), &
                  L_J_kc,                      &
                  wf%eri_t1%n_J,                  &
                  zero,                        &
                  X_bj_J,                      &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v))
!
      call mem%dealloc(L_J_kc, wf%eri_t1%n_J, n_a_o, n_a_v)
!
      call mem%alloc(L_J_ai, wf%eri_t1%n_J, wf%n_ccsd_v, wf%n_ccsd_o)
      call wf%L_t1%get(L_J_ai, wf%n_o + 1, wf%n_o + wf%n_ccsd_v, 1, wf%n_ccsd_o)
!
!     omega2_aibj =+ L_J_ai X_bj_J
!
      call dgemm('T', 'T',                     &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v), &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v), &
                  wf%eri_t1%n_J,                  &
                  one,                         &
                  L_J_ai,                      &
                  wf%eri_t1%n_J,                  &
                  X_bj_J,                      &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v), &
                  one,                         &
                  omega2_aibj,                 &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v))
!
      call mem%dealloc(X_bj_J, wf%n_ccsd_v, wf%n_ccsd_o, wf%eri_t1%n_J)
      call mem%dealloc(L_J_ai, wf%eri_t1%n_J, wf%n_ccsd_v, wf%n_ccsd_o)
      call mem%dealloc(u_aild,  wf%n_ccsd_v, wf%n_ccsd_o, n_a_v, n_a_o)
!
!     Add the D2.1 term to the omega vector
!
      call symmetric_sum(omega2_aibj, (wf%n_ccsd_o)*(wf%n_ccsd_v))
!
      call daxpy(((wf%n_ccsd_o)*(wf%n_ccsd_v))**2, one, omega2_aibj, 1, omega_aibj, 1)
!
      call mem%dealloc(omega2_aibj, (wf%n_ccsd_v), (wf%n_ccsd_o), (wf%n_ccsd_v), (wf%n_ccsd_o))
!
      call timer%turn_off()
!
   end subroutine omega_ccsd_e2_mlccsd
!
!
   module subroutine omega_ccsd_f2_mlccsd(wf, omega_aibj, x_aibj)
!!
!!    Omega F2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!    Modified for MLCCSD by Sarai D. Folkestad 2017-2019
!!
!!    Calculates the F2 term,
!!
!!      F2: sum_c x_ij^ac (F_bc - sum_dkl g_ldkc u_kl^bd)
!!        - sum_k x_ik^ab (F_kj + sum_cdl g_ldkc u_lj^dc),
!!
!!    where
!!
!!        u_kl^bc = 2 * x_kl^bc - x_lk^bc.
!!
!!    INDEX RESTRICTIONS:
!!
!!    a, i, b, j are CCSD indices
!!
!!    c, d, k, l are CC2 + CCSD indices
!!
      use reordering, only: symmetric_sum
!
      implicit none
!
      class(mlccsd) :: wf
!
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o), &
                  intent(inout):: omega_aibj
!
      real(dp), dimension(wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o, &
                          wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o), intent(in) :: x_aibj
!
!
      real(dp), dimension(:,:,:,:), allocatable :: omega2_bjai
      real(dp), dimension(:,:,:,:), allocatable :: g_ldkc
      real(dp), dimension(:,:,:,:), allocatable :: u_bldk
      real(dp), dimension(:,:,:,:), allocatable :: u_cldj
      real(dp), dimension(:,:,:,:), allocatable :: x_cjai
      real(dp), dimension(:,:,:,:), allocatable :: x_aibk
!
      real(dp), dimension(:,:), allocatable :: X_bc
      real(dp), dimension(:,:), allocatable :: Y_kj
!
      type(timings) :: timer
!
      integer :: k, l, d, b, c, a, i, j
!
      integer :: n_a_o, n_a_v
!
      timer = timings('Omega MLCCSD F2', pl='v')
      call timer%turn_on()
!
      n_a_o = wf%n_cc2_o + wf%n_ccsd_o
      n_a_v = wf%n_cc2_v + wf%n_ccsd_v
!
!     Form u_bkdl = 2 * x_kl^bd - x_lk^bd  ordered as u_bldk
!
      call mem%alloc(u_bldk,  wf%n_ccsd_v, n_a_o, n_a_v, n_a_o)
!
!$omp parallel do private(k, d, l, b) collapse(2)
      do k = 1, n_a_o
         do d = 1, n_a_v
            do l = 1, n_a_o
               do b = 1, wf%n_ccsd_v
!
                  u_bldk(b, l, d, k) = two*x_aibj(b,k,d,l) &
                                       - x_aibj(b,l,d,k)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
!
!     Form g_ldkc = g_ldkc
!
      call mem%alloc(g_ldkc, n_a_o, n_a_v, n_a_o, n_a_v)
      call wf%eri_t1%get('ovov', g_ldkc, 1, n_a_o, 1, n_a_v, 1, n_a_o, 1, n_a_v)
!
!     Make the intermediate X_b_c = F_bc - sum_dkl g_ldkc u_kl^bd and set to zero
!
      call mem%alloc(X_bc, wf%n_ccsd_v, n_a_v)
!
!     First, copy the virtual-virtual Fock matrix into the intermediate
!
!$omp parallel do private(c, b) collapse(2)
      do c = 1, n_a_v
         do b = 1, wf%n_ccsd_v
!
            X_bc(b,c) = wf%fock_ab(b, c)
!
         enddo
      enddo
!$omp end parallel do
!
!     Then, add the second contribution,
!     - sum_dkl g_ldkc u_kl^bd = - sum_dkl u_b_kdl * g_kdl_c, to X_b_c
!
      call dgemm('N','N',                          &
                  wf%n_ccsd_v,                     &
                  n_a_v,                           &
                  (n_a_v)*(n_a_o**2),              &
                  -one,                            &
                  u_bldk,                          & ! u_b_ldk
                  wf%n_ccsd_v,                     &
                  g_ldkc,                          & ! g_ldk_c
                  (n_a_v)*(n_a_o**2),              &
                  one,                             &
                  X_bc,                            &
                  wf%n_ccsd_v)
!
      call mem%dealloc(u_bldk, wf%n_ccsd_v, n_a_o, n_a_v, n_a_o)
!
      call mem%alloc(omega2_bjai, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
      call mem%alloc(x_cjai, n_a_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
!$omp parallel do private(i, j, a, c) collapse(2)
      do i = 1, wf%n_ccsd_o
         do a = 1, wf%n_ccsd_v
            do j = 1, wf%n_ccsd_o
               do c = 1, n_a_v
!
                  x_cjai(c, j, a, i) = x_aibj(c, j, a, i)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call dgemm('N','N',                          &
                  wf%n_ccsd_v,                     &
                  (wf%n_ccsd_v)*(wf%n_ccsd_o)**2,  &
                  n_a_v,                           &
                  one,                             &
                  X_bc,                            &
                  wf%n_ccsd_v,                     &
                  x_cjai,                          &
                  n_a_v,                           &
                  zero,                            &
                  omega2_bjai,                     &
                  wf%n_ccsd_v)
!
      call mem%dealloc(X_bc, wf%n_ccsd_v, n_a_v)
      call mem%dealloc(x_cjai, n_a_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
!     Make the intermediate Y_k_j = F_kj  + sum_cdl u_lj^dc g_ldkc
!
      call mem%alloc(Y_kj, n_a_o, wf%n_ccsd_o)
!
!     Copy the occupied-occupied Fock matrix, such that Y_k_j = F_kj
!
!$omp parallel do private (j,k) collapse(2)
      do j = 1, wf%n_ccsd_o
         do k = 1, n_a_o
!
            Y_kj(k,j) = wf%fock_ij(k, j)
!
         enddo
      enddo
!$omp end parallel do
!
!     Add sum_cdl g_k_dlc u_dlc_j to Y_k_j, such that
!     Y_k_j = F_k_j + sum_cdl g_k_cld u_cld_j
!
!     Note:
!
!     g_ldkc(kc,ld) = g_kcld -> pretend that this is g_k_cld
!
      call mem%alloc(u_cldj,  n_a_v, n_a_o, n_a_v, wf%n_ccsd_o)
!
!$omp parallel do private (j, d, l, c) collapse(2)
      do j = 1, wf%n_ccsd_o
         do d = 1, n_a_v
            do l = 1, n_a_o
               do c = 1, n_a_v
!
                  u_cldj(c, l, d, j) = two*x_aibj(c, j, d, l) &
                                       - x_aibj(c, l, d, j)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call dgemm('N','N',             &
                 n_a_o,               &
                 wf%n_ccsd_o,         &
                 (n_a_o)*(n_a_v**2),  &
                 one,                 &
                 g_ldkc,              & ! g_k_cld
                 n_a_o,               &
                 u_cldj,              &
                 (n_a_o)*(n_a_v**2),  &
                 one,                 &
                 Y_kj,                &
                 n_a_o)
!
      call mem%dealloc(u_cldj, n_a_v, n_a_o, n_a_v, wf%n_ccsd_o)
      call mem%dealloc(g_ldkc, n_a_o, n_a_v, n_a_o, n_a_v)
!
!     - sum_k x_aib_k Y_k_j = - sum_k x_ik^ab (F_kj + sum_cdl g_ldkc u_lj^dc)
!
      call mem%alloc(x_aibk, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, n_a_o)
!
!$omp parallel do private (k, b, i, a) collapse(2)
      do k = 1, n_a_o
         do b = 1, wf%n_ccsd_v
            do i = 1, wf%n_ccsd_o
               do a = 1, wf%n_ccsd_v
!
                  x_aibk(a, i, b, k) = x_aibj(a, i, b, k)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call dgemm('N','N',                          &
                 (wf%n_ccsd_o)*(wf%n_ccsd_v)**2,   &
                 wf%n_ccsd_o,                      &
                 n_a_o,                            &
                 -one,                             &
                 x_aibk,                           &
                 (wf%n_ccsd_o)*(wf%n_ccsd_v)**2,   &
                 Y_kj,                             &
                 n_a_o,                            &
                 one,                              &
                 omega2_bjai,                      & ! omega2_aibj but we will symmetrize anyway
                 (wf%n_ccsd_o)*(wf%n_ccsd_v)**2)
!
!
      call mem%dealloc(Y_kj, n_a_o, wf%n_ccsd_o)
      call mem%dealloc(x_aibk, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, n_a_o)
!
      call symmetric_sum(omega2_bjai, (wf%n_ccsd_o)*(wf%n_ccsd_v))
      call daxpy(((wf%n_ccsd_o)*(wf%n_ccsd_v))**2, one, omega2_bjai, 1, omega_aibj, 1)
!
      call mem%dealloc(omega2_bjai, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
      call timer%turn_off()
!
   end subroutine omega_ccsd_f2_mlccsd
!
!
end submodule omega_mlccsd
