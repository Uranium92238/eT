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
submodule (doubles_class) omega_doubles
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
!
   module subroutine omega_doubles_a1_doubles(wf, omega, u)
!!
!!    Omega doubles A1 term
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, Jan 2019
!!
!!    Calculates the A1 term,
!!
!!       A1: sum_ckd u_bicj g_abjc = sum_ckd u_bjc_i * g_a_bjc,
!!
!!    and adds it to the projection vector omega
!!
      implicit none
!
      class(doubles), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: omega
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: u
!
      real(dp), dimension(:,:,:,:), allocatable :: g_abjc, u_bjci
!
      real(dp), dimension(:,:), allocatable :: omega_ai
!
      type(batching_index) :: batch_a
!
      integer :: req0, req1
!
      integer :: a, i, current_a_batch
!
      type(timings) :: timer 
!  
      timer = timings('omega doubles a1', pl='verbose')
      call timer%turn_on()
!
      call mem%alloc(u_bjci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Reorder u_bicj to u_bjci
!
      call sort_1234_to_1432(u, u_bjci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      req0 = (wf%n_o)*(wf%n_v)*(wf%integrals%n_J)
      req1 = (wf%n_v)**2*(wf%n_o) + (wf%n_v)*(wf%integrals%n_J)
!
      batch_a = batching_index(wf%n_v)
!
      call mem%batch_setup(batch_a, req0, req1)
!
      do current_a_batch = 1, batch_a%num_batches
!
         call batch_a%determine_limits(current_a_batch)
!
         call mem%alloc(g_abjc, batch_a%length, wf%n_v, wf%n_o, wf%n_v)
!
         call wf%get_vvov(g_abjc,                       &
                           batch_a%first, batch_a%last, &
                           1, wf%n_v,                   &
                           1, wf%n_o,                   &
                           1, wf%n_v)
!
         call mem%alloc(omega_ai, batch_a%length, wf%n_o)
!
         call dgemm('N','N',               &
                     batch_a%length,       &
                     wf%n_o,               &
                     (wf%n_o)*(wf%n_v)**2, &
                     one,                  &
                     g_abjc,               & ! g_a_bjc
                     batch_a%length,       &
                     u_bjci,               & ! u_bjc_i
                     (wf%n_o)*(wf%n_v)**2, &
                     zero,                 &
                     omega_ai,             &
                     batch_a%length)
!
         call mem%dealloc(g_abjc, batch_a%length, wf%n_v, wf%n_o, wf%n_v)
!
!$omp parallel do private(i, a)
         do i = 1, wf%n_o
            do a = 1, batch_a%length
!
               omega(a + batch_a%first - 1, i) = omega(a + batch_a%first - 1, i) &
                                                + omega_ai(a, i)
!
            enddo
         enddo
!$omp end parallel do
!
         call mem%dealloc(omega_ai, batch_a%length, wf%n_o)
!
      enddo ! batch_a
!
      call mem%dealloc(u_bjci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine omega_doubles_a1_doubles
!
!
   module subroutine omega_doubles_b1_doubles(wf, omega, u)
!!
!!    Omega doubles B1 term
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, Jan 2019
!!
!!    Calculates the B1 term,
!!
!!       B1: - sum_ckl g_kb,ji * u_aj,bk,
!!
!!    with
!!
!!      u_aj_bk = 2t_aj,bk - t_ak,bj
!!
!!    and adds it to the projection vector (omega)
!!
      implicit none
!
      class(doubles), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: omega
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: u
!
      real(dp), dimension(:,:,:,:), allocatable :: g_kbji
      real(dp), dimension(:,:,:,:), allocatable :: g_jbki
!
      type(timings) :: timer 
!  
      timer = timings('omega doubles b1', pl='verbose')
      call timer%turn_on()
!
!     g_kbji ordered as g_jbki
!
      call mem%alloc(g_kbji, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call wf%get_ovoo(g_kbji)
!
      call mem%alloc(g_jbki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call sort_1234_to_3214(g_kbji, g_jbki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call mem%dealloc(g_kbji, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
!     omega_ai += - sum_ckl g_kb,ji * u_aj,bk
!
      call dgemm('N', 'N',             &
                  wf%n_v,              &
                  wf%n_o,              &
                  (wf%n_o**2)*wf%n_v,  &
                  -one,                &
                  u,                   & ! u_a_jbk
                  wf%n_v,              &
                  g_jbki,              & ! g_jbk_i
                  (wf%n_o**2)*wf%n_v,  &
                  one,                 &
                  omega,               & ! omega_a_i
                  wf%n_v)
!
      call mem%dealloc(g_jbki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine omega_doubles_b1_doubles
!
!
   module subroutine omega_doubles_c1_doubles(wf, omega, u)
!!
!!    Omega doubles C1 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2019
!!
!!    Calculates the C1 term,
!!
!!       C1: sum_bj u_ai,bj * F_jb,
!!
!!    with
!!
!!       u_ai_bj = 2*t_ai_bj - t_aj_bi
!!
!!
      implicit none
!
      class(doubles), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: omega
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: u
!
      real(dp), dimension(:,:), allocatable :: F_bj
!
      type(timings) :: timer
!  
      timer = timings('omega doubles c1', pl='verbose')
      call timer%turn_on()
!
      call mem%alloc(F_bj, wf%n_o, wf%n_v)
      call sort_12_to_21(wf%fock_ia, F_bj, wf%n_o, wf%n_v)
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  1,                 &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  u,                 & ! u_ai_bj
                  (wf%n_o)*(wf%n_v), &
                  F_bj,              &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  omega,             &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(F_bj, wf%n_o, wf%n_v)
!
      call timer%turn_off()
!
    end subroutine omega_doubles_c1_doubles
!
!
end submodule omega_doubles
