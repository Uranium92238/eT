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
submodule (lowmem_cc2_class) omega_lowmem_cc2
!
!!
!!    Omega submodule 
!!
!!    Routines to construct 
!!
!!    Ω =  < mu | exp(-T) H exp(T) | R >
!!
!
   implicit none
!
!
contains
!
   module subroutine construct_omega_lowmem_cc2(wf, omega)
!!
!!    Construct omega 
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 
!!    Linda Goletto, and Alexander C. Paul, Dec 2018
!!
!!    Direqts the construction of the projection vector < mu | exp(-T) H exp(T) | R >
!!    for the current wavefunction amplitudes.
!!
      implicit none
!
      class(lowmem_cc2), intent(inout) :: wf
!
      real(dp), dimension(wf%n_t1), intent(out) :: omega
!
      type(timings), allocatable :: timer 
!
      timer = timings('Construct omega lowmem-cc2', pl='normal')
      call timer%turn_on()
!
      call zero_array(omega, wf%n_t1)
!
      call wf%omega_ccs_a1(omega)
!
      call wf%omega_cc2_a1(omega, wf%orbital_energies(1:wf%n_o), wf%orbital_energies(wf%n_o + 1 : wf%n_mo))
      call wf%omega_cc2_b1(omega, wf%orbital_energies(1:wf%n_o), wf%orbital_energies(wf%n_o + 1 : wf%n_mo))
      call wf%omega_cc2_c1(omega, wf%orbital_energies(1:wf%n_o), wf%orbital_energies(wf%n_o + 1 : wf%n_mo))
!
      call timer%turn_off()
!
   end subroutine construct_omega_lowmem_cc2
!
!
   module subroutine omega_cc2_a1_lowmem_cc2(wf, omega, eps_o, eps_v)
!!
!!    Omega CC2 A1 term
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 
!!    Linda Goletto, and Alexander C. Paul, Dec 2018
!!
!!    Calculates the A1 term,
!!
!!       A1: sum_cjb u_bj_ci * g_abjc,
!!
!!    with 
!!       
!!       u_bj_ci = 2*t_bj_ci - t_bi_cj
!!
!!    and
!!
!!       t_bj_ci = - g_bjci/ε^{bc}_{ji}
!!
!!    and adds it to the projection vector omega
!!
      implicit none
!
      class(lowmem_cc2), intent(inout) :: wf
!
      real(dp), dimension(wf%n_t1), intent(inout) :: omega
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
      real(dp), dimension(:,:,:,:), allocatable :: g_bicj
      real(dp), dimension(:,:,:,:), allocatable :: u_bjci
      real(dp), dimension(:,:,:,:), allocatable :: g_abjc
!
      real(dp) :: eps_ci
!
      integer :: b, j, c, i
!
      integer :: req0, req1_b, req1_c, req2
!
      integer :: current_b_batch, current_c_batch
!
      type(batching_index) :: batch_b, batch_c
!
      type(timings), allocatable :: timer 
!
      timer = timings('omega cc2 a1 lowmem', pl='verbose')
      call timer%turn_on()
!
      req0 = 0
!
      req1_b = (wf%eri%n_J)*(wf%n_v)
      req1_c = (wf%eri%n_J)*(wf%n_o)
!
      req2 =  2*(wf%n_o**2) + (wf%n_o)*(wf%n_v)
!
      batch_b = batching_index(wf%n_v)
      batch_c = batching_index(wf%n_v)
!
      call mem%batch_setup(batch_b, batch_c, req0, req1_b, req1_c, req2)
!
      do current_b_batch = 1, batch_b%num_batches
!
         call batch_b%determine_limits(current_b_batch)
!
         do current_c_batch = 1, batch_c%num_batches
!
            call batch_c%determine_limits(current_c_batch)
!
!           Construct u_bicj = -2g_bicj + g_bjci/(e_bjci), ordered as u_bjci
!
            call mem%alloc(g_bicj, batch_b%length,wf%n_o, batch_c%length,wf%n_o)
            call mem%alloc(u_bjci, batch_b%length,wf%n_o, batch_c%length,wf%n_o)
!
            call wf%eri%get_eri_t1('vovo', g_bicj, batch_b%first, batch_b%get_last(), 1, wf%n_o, &
                                                   batch_c%first, batch_c%get_last(), 1, wf%n_o)
!
!$omp parallel do schedule(static) private(i, j, c, b, eps_ci) collapse(2)
            do i = 1, wf%n_o
               do c = 1, (batch_c%length)
!
                  eps_ci = eps_v(c + batch_c%first - 1) - eps_o(i)
!
                  do  j = 1, wf%n_o
                     do b = 1, (batch_b%length)
!                 
                        u_bjci(b,j,c,i) = -(two*g_bicj(b,i,c,j) - g_bicj(b,j,c,i))&
                                                                /(eps_ci +  &
                                                                  eps_v(b + batch_b%first - 1)&
                                                                - eps_o(j))
                      enddo
                   enddo
               enddo
            enddo
!$omp end parallel do
!
            call mem%dealloc(g_bicj, batch_b%length,wf%n_o, batch_c%length,wf%n_o)
!
!           Omega_ai += sum_bjc u_bicj g_abjc       
!
            call mem%alloc(g_abjc, batch_b%length,wf%n_v, batch_c%length,wf%n_o)
!
            call wf%eri%get_eri_t1('vvov', g_abjc, 1, wf%n_v, batch_b%first, batch_b%get_last(), &
                                                   1, wf%n_o, batch_c%first, batch_c%get_last())
!
            call dgemm('N','N',                                   &
                        wf%n_v,                                   &
                        wf%n_o,                                   &
                        (batch_b%length)*(batch_c%length)*wf%n_o, &
                        one,                                      &
                        g_abjc,                                   &
                        wf%n_v,                                   &
                        u_bjci,                                   &
                        (batch_b%length)*(batch_c%length)*wf%n_o, &
                        one,                                      &
                        omega,                                    &
                        wf%n_v)
!
            call mem%dealloc(g_abjc, batch_b%length, wf%n_v, batch_c%length, wf%n_o)
            call mem%dealloc(u_bjci, batch_b%length, wf%n_o, batch_c%length, wf%n_o)
!
         enddo
      enddo
!
      call mem%batch_finalize()
!
      call timer%turn_off()
!
   end subroutine omega_cc2_a1_lowmem_cc2
!
!
   module subroutine omega_cc2_b1_lowmem_cc2(wf, omega, eps_o, eps_v)
!!
!!    Omega CC2 B1 term
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 
!!    Linda Goletto, and Alexander C. Paul, Dec 2018
!!
!!    Calculates the B1 term,
!!
!!       B1: - sum_ckl (2g_kb_ji - g_jb_ki) * t_aj_bk,
!!
!!    with
!!
!!       t_aj_bk = - g_ajbk/ε^{ab}_{jk}
!!
!!    and adds it to the projection vector (omega) of
!!    the wavefunction object wf.
!!
      implicit none
!
      class(lowmem_cc2), intent(inout) :: wf
!
      real(dp), dimension(wf%n_t1), intent(inout) :: omega
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ajbk
      real(dp), dimension(:,:,:,:), allocatable :: g_jbki
      real(dp), dimension(:,:,:,:), allocatable :: g_kbji
!
      integer :: a, b, j, k
!
      integer :: req0, req1_b, req1_j, req1_k, req2_bj, req2_bk, req2_jk, req3
!
      integer :: current_b_batch, current_j_batch, current_k_batch
!
      type(batching_index) :: batch_b, batch_j, batch_k
!
      type(timings), allocatable :: timer 
!
      timer = timings('omega cc2 b1 lowmem', pl='verbose')
      call timer%turn_on()
!
      req0 = 0
!
      req1_b = wf%eri%n_J
      req1_j = wf%eri%n_J*(wf%n_v)
      req1_k = wf%eri%n_J*(wf%n_o)
!
      req2_bj = wf%eri%n_J
      req2_bk = wf%eri%n_J
      req2_jk = 0
!
      req3 = (wf%n_v) + 2*(wf%n_o) 
!
      batch_b = batching_index(wf%n_v)
      batch_j = batching_index(wf%n_o)
      batch_k = batching_index(wf%n_o)
!
      call mem%batch_setup(batch_b, batch_j, batch_k, req0, req1_b, &
                    req1_j, req1_k, req2_bj, req2_bk, req2_jk, req3)
!
      do current_b_batch = 1, batch_b%num_batches
!
         call batch_b%determine_limits(current_b_batch)
!
         do current_j_batch = 1, batch_j%num_batches
!
            call batch_j%determine_limits(current_j_batch)
!
            do current_k_batch = 1, batch_k%num_batches
!
               call batch_k%determine_limits(current_k_batch)
!
!               Construct t_ajbk = - g_ajbk/e_ajbk, 
!               store in g_ajbk to avoid double memory requirement
!
               call mem%alloc(g_ajbk, wf%n_v, batch_j%length,&
                              batch_b%length, batch_k%length)
!
               call wf%eri%get_eri_t1('vovo', g_ajbk, 1, wf%n_v, &
                                                      batch_j%first, batch_j%get_last(), &
                                                      batch_b%first, batch_b%get_last(), &
                                                      batch_k%first, batch_k%get_last())
!
!$omp parallel do schedule(static) private(k, j, a, b)
               do a = 1, wf%n_v
                  do j = 1, (batch_j%length)
                     do b = 1, (batch_b%length)
                        do k = 1, (batch_k%length)
!
                           g_ajbk(a,j,b,k) = - g_ajbk(a,j,b,k)/(eps_v(a) &
                                                   + eps_v(b + batch_b%first - 1)&
                                                   - eps_o(j + batch_j%first - 1) &
                                                   - eps_o(k + batch_k%first - 1))
!
                        enddo
                     enddo
                  enddo
               enddo
!$omp end parallel do
!
               call mem%alloc(g_jbki, batch_j%length, batch_b%length, &
                              batch_k%length, wf%n_o)
!
               call wf%eri%get_eri_t1('ovoo', g_jbki,              &
                                      batch_j%first, batch_j%get_last(), &
                                      batch_b%first, batch_b%get_last(), &
                                      batch_k%first, batch_k%get_last(), & 
                                      1, wf%n_o)
!
!              Omega_ai += sum_jbk t_ajbk g_jbki       
!
               call dgemm('N','N',                                             &
                           wf%n_v,                                             &
                           wf%n_o,                                             &
                           (batch_j%length)*(batch_k%length)*(batch_b%length), &
                           one,                                                &
                           g_ajbk,                                             &
                           wf%n_v,                                             &
                           g_jbki,                                             &
                           (batch_j%length)*(batch_k%length)*(batch_b%length), &
                           one,                                                &
                           omega,                                              &
                           wf%n_v)
!
               call mem%dealloc(g_jbki, batch_j%length, batch_b%length, &
                                batch_k%length, wf%n_o)
!
!              Construct g_kbij and reorder as g_jbki
!
               call mem%alloc(g_kbji, batch_k%length, batch_b%length, &
                                      batch_j%length, wf%n_o)
!
               call wf%eri%get_eri_t1('ovoo', g_kbji,              &
                                      batch_k%first, batch_k%get_last(), &
                                      batch_b%first, batch_b%get_last(), &
                                      batch_j%first, batch_j%get_last(), & 
                                      1, wf%n_o)
!
               call mem%alloc(g_jbki, batch_j%length, batch_b%length, &
                              batch_k%length, wf%n_o)
!
               call sort_1234_to_3214(g_kbji, g_jbki, (batch_k%length), &
                     (batch_b%length), (batch_j%length), wf%n_o)
!
               call mem%dealloc(g_kbji, batch_k%length, batch_b%length, &
                                        batch_j%length, wf%n_o)
!
!              Omega_ai += sum_bjk t_ajbk g_kbji
!
               call dgemm('N','N',                                             &
                           wf%n_v,                                             &
                           wf%n_o,                                             &
                           (batch_j%length)*(batch_k%length)*(batch_b%length), &
                           -two,                                               &
                           g_ajbk,                                             &
                           wf%n_v,                                             &
                           g_jbki,                                             &
                           (batch_j%length)*(batch_k%length)*(batch_b%length), &
                           one,                                                &
                           omega,                                              &
                           wf%n_v)
!
               call mem%dealloc(g_jbki, batch_j%length, batch_b%length, &
                                        batch_k%length, wf%n_o)
               call mem%dealloc(g_ajbk, wf%n_v, batch_j%length,&
                                      batch_b%length, batch_k%length)
!
            enddo
         enddo
      enddo
!
      call mem%batch_finalize()
!
      call timer%turn_off()
!
   end subroutine omega_cc2_b1_lowmem_cc2
!
!
   module subroutine omega_cc2_c1_lowmem_cc2(wf, omega, eps_o, eps_v)
!!
!!    Omega CC2 C1 term
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 
!!    Linda Goletto, and Alexander C. Paul, Dec 2018
!!
!!    Calculates the C1 term,
!!
!!       C1: sum_bj u_ai_bj * F_{jb},
!!
!!    with 
!!       
!!       u_ai_bj = 2*t_ai_bj - t_aj_bi
!!
!!    and
!!
!!       t_ai_bj = - g_aibj/ε^{ab}_{ij}
!!
!!    and adds it to the projection vector (omega) of
!!    the wavefunction object wf.
!!
      implicit none
!
      class(lowmem_cc2), intent(inout) :: wf
!
      real(dp), dimension(wf%n_t1), intent(inout) :: omega
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
      real(dp), dimension(:,:,:,:), allocatable :: g_aibj
      real(dp), dimension(:,:,:,:), allocatable :: u_aibj
      real(dp), dimension(:,:), allocatable :: F_bj
!
      integer :: i, j, a, b
!
      integer :: req0, req1_j, req1_i, req2, omega_offset
!
      integer :: current_j_batch, current_i_batch
!
      type(batching_index) :: batch_j, batch_i
!
      type(timings), allocatable :: timer 
!
      timer = timings('omega cc2 c1 lowmem', pl='verbose')
      call timer%turn_on()
!
      req0 = 0
!
      req1_j = (wf%n_v)*(wf%eri%n_J)
      req1_i = (wf%n_v)*(wf%eri%n_J)
!
      req2 =  2*(wf%n_v)**2
!
      batch_i = batching_index(wf%n_o)
      batch_j = batching_index(wf%n_o)
!
      call mem%batch_setup(batch_i, batch_j, req0, req1_i, req1_j, req2)
!
      do current_i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(current_i_batch)
!
         do current_j_batch = 1, batch_j%num_batches
!
            call batch_j%determine_limits(current_j_batch)
!
!           Construct u_aibj = -(2g_aibj - g_ajbi)/e_aibj
!
            call mem%alloc(g_aibj, wf%n_v, batch_i%length, wf%n_v, batch_j%length)
!
            call wf%eri%get_eri_t1('vovo', g_aibj, 1, wf%n_v, batch_i%first, batch_i%get_last(),  &
                                                   1, wf%n_v, batch_j%first, batch_j%get_last())
!
            call mem%alloc(u_aibj, wf%n_v, batch_i%length, wf%n_v, batch_j%length)
!
!           Construct u_aibj
!
!$omp parallel do schedule(static) private(i, j, a, b)
            do b = 1, wf%n_v 
               do j = 1, batch_j%length
                  do i = 1, batch_i%length 
                     do a = 1, wf%n_v
!
                        u_aibj(a,i,b,j) = (-two*g_aibj(a,i,b,j)+g_aibj(b,i,a,j))/(eps_v(a) &
                                                             + eps_v(b) &
                                                             - eps_o(i + batch_i%first - 1) &
                                                             - eps_o(j + batch_j%first - 1)) 
!
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
            call mem%dealloc(g_aibj, wf%n_v, batch_i%length, wf%n_v, batch_j%length)
!
!           Reorder F_jb as F_bj
!
            call mem%alloc(F_bj, wf%n_v, batch_j%length)
!
!$omp parallel do schedule(static) private(b,j)
            do b = 1, wf%n_v 
               do j = 1, batch_j%length
!
                     F_bj(b, j) = wf%fock_ia(j + batch_j%first - 1, b)
!
               enddo
            enddo
!$omp end parallel do
!
!           Omega_ai += u_aibj F_jb
!
            omega_offset = wf%n_v*(batch_i%first - 1) + 1
!
            call dgemm('N', 'N',                     &
                       (batch_i%length)*wf%n_v,      &
                       1,                            &
                       (batch_j%length)*wf%n_v,      &
                       one,                          &
                       u_aibj,                       &
                       wf%n_v*(batch_i%length),      &
                       F_bj,                         &
                       wf%n_v*(batch_j%length),      &
                       one,                          &
                       omega(omega_offset),          &
                       wf%n_v*wf%n_o)
!        
            call mem%dealloc(u_aibj, wf%n_v, batch_i%length, wf%n_v, batch_j%length)
            call mem%dealloc(F_bj, wf%n_v, batch_j%length)
!
         enddo
      enddo
!
      call mem%batch_finalize()
!
      call timer%turn_off()
!
   end subroutine omega_cc2_c1_lowmem_cc2
!
!
end submodule omega_lowmem_cc2
