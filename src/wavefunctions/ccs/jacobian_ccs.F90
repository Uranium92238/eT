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
submodule (ccs_class) jacobian_ccs
!
!!
!!    Jacobian submodule
!!
!!    Routines for the linear transform of trial
!!    vectors by the Jacobian matrix
!!
!!    ρ_i = A * c_i,
!!
!!    where
!!
!!    A_μ,ν = < μ | exp(-T) [H, τ_ν] exp(T) | R >.
!!
!
   implicit none
!
!
contains
!
!
   module subroutine prepare_for_jacobian_ccs(wf)
!!
!!    Prepare for jacobian
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2019
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
!     For now, do nothing.
!
      call output%printf('v', '- No preparations for the ' // trim(wf%name_) // &
                         ' excited state equation.', fs='(/t3,a)')
!
   end subroutine prepare_for_jacobian_ccs
!
!
   module subroutine jacobian_transformation_ccs(wf, c)
!!
!!    Jacobian transformation
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!    Directs the transformation by the CCSD Jacobi matrix,
!!
!!       A_μ,ν = < μ | exp(-T) [H, τ_ν] exp(T) | R >.
!!
!!    In particular,
!!
!!       rho_mu = (A c)_mu = sum_ck A_mu,ck c_ck.
!!
!!    On exit, c is overwritten by rho.
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: c
!
      real(dp), dimension(:,:), allocatable :: rho
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian CCS transformation', pl='normal')
      call timer%turn_on()
!
!     Allocate the transformed vector & add the terms to it
!
      call mem%alloc(rho, wf%n_v, wf%n_o)
      call zero_array(rho, wf%n_t1)
!
      call wf%jacobian_ccs_a1(rho, c)
      call wf%jacobian_ccs_b1(rho, c)
!
!     Then overwrite the c vector with the transformed vector
!
      call dcopy((wf%n_o)*(wf%n_v), rho, 1, c, 1)
      call mem%dealloc(rho, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_transformation_ccs
!
!
   module subroutine jacobian_ccs_a1_ccs(wf, rho1, c1)
!!
!!    Jacobian CCS A1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!    Calculates the A1 term,
!!
!!       sum_b F_ab c_bi - sum_j F_ji c_aj,
!!
!!    and adds it to the rho vector.
!!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c1
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho1
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian CCS A1', pl='verbose')
      call timer%turn_on()
!
!     sum_b F_a_b c_b_i 
!
      call dgemm('N', 'N',     &
                  wf%n_v,      &
                  wf%n_o,      &
                  wf%n_v,      &
                  one,         &
                  wf%fock_ab,  &
                  wf%n_v,      &
                  c1,          &
                  wf%n_v,      &
                  one,         &
                  rho1,        &
                  wf%n_v)
!
!     - sum_j c_a_j F_j_i
!
      call dgemm('N','N',      &
                  wf%n_v,      &
                  wf%n_o,      &
                  wf%n_o,      &
                  -one,        &
                  c1,          &
                  wf%n_v,      &
                  wf%fock_ij,  &
                  wf%n_o,      &
                  one,         &
                  rho1,        &
                  wf%n_v)
!
      
      call timer%turn_off()
!
   end subroutine jacobian_ccs_a1_ccs
!
!
   module subroutine jacobian_ccs_b1_ccs(wf, rho_ai, c_bj)
!!
!!    Jacobian CCS B1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,
!!    Linda Goletto, and Alexander C. Paul, Dec 2018
!!
!!    Calculates the A1 term
!!
!!       A1: sum_bj (2 g_aijb - g_abji) * c_bj
!!
!!    and adds it to rho_ai.
!!
!!    Separate calculation of both terms due to batching.
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_bj
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
!
      real(dp), dimension(:,:), allocatable :: c_jb
!
      real(dp), dimension(:,:,:,:), allocatable :: g_aijb, g_abji
!
      integer :: req0, req1_i, req1_j, req1_b, req2
!
      integer :: current_i_batch, current_j_batch, current_b_batch
!
      integer :: rho_offset, j, b
!
      type(batching_index) :: batch_i, batch_j, batch_b
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian CCS B1', pl='verbose')
      call timer%turn_on()
!
!     :: Term 1: rho_ai = sum_bj 2 g_aijb * c_bj ::
!
      req0 = 0
!
      req1_i = (wf%integrals%n_J)*(wf%n_v)
      req1_j = (wf%integrals%n_J)*(wf%n_v)
!
      req2 = (wf%n_v)**2
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
            call mem%alloc(g_aijb, wf%n_v, (batch_i%length), (batch_j%length), wf%n_v)
!
            call wf%get_voov(g_aijb,                        &
                              1, wf%n_v,                    &
                              batch_i%first, batch_i%last,  &
                              batch_j%first, batch_j%last,  &
                              1, wf%n_v)
!
            call mem%alloc(c_jb, (batch_j%length), wf%n_v)
!
!$omp parallel do private(j, b)
            do b = 1, wf%n_v
               do j = 1, batch_j%length
!
                  c_jb(j, b) = c_bj(b, j + batch_j%first - 1)
!
               enddo
            enddo
!$omp end parallel do
!
!           rho_a_i = rho_a_i + sum_bj 2 g_aijb * c_bj
!
            rho_offset = wf%n_v*(batch_i%first - 1) + 1
!
            call dgemm('N', 'N',                   &
                        (wf%n_v)*(batch_i%length), &
                        1,                         &
                        (batch_j%length)*(wf%n_v), &
                        two,                       &
                        g_aijb,                    & ! g_ai_jb
                        (wf%n_v)*(batch_i%length), &
                        c_jb,                      & ! c_jb
                        (batch_j%length)*(wf%n_v), &
                        one,                       &
                        rho_ai(1, batch_i%first),  & 
                        (wf%n_v)*(wf%n_o))
!
            call mem%dealloc(c_jb, (batch_j%length), wf%n_v)
            call mem%dealloc(g_aijb, wf%n_v,(batch_i%length), (batch_j%length), wf%n_v)
!
         enddo ! batch_j
      enddo !batch_i
!
!     :: Term 2 rho_ai = - g_abji * c_bj::
!
      req1_i = max((wf%integrals%n_J)*(wf%n_v), (wf%integrals%n_J)*(wf%n_o))
      req1_b = max((wf%integrals%n_J)*(wf%n_v), (wf%integrals%n_J)*(wf%n_o))
!
      req2 = 2*(wf%n_o)*(wf%n_v)
!
      batch_i = batching_index(wf%n_o)
      batch_b = batching_index(wf%n_v)
!
      call mem%batch_setup(batch_i, batch_b, req0, req1_i, req1_b, req2)
!
      do current_i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(current_i_batch)
!
         do current_b_batch = 1, batch_b%num_batches
!
            call batch_b%determine_limits(current_b_batch)
!
            call mem%alloc(g_abji, wf%n_v, (batch_b%length), wf%n_o, (batch_i%length))
!
            call wf%get_vvoo(g_abji,                        &
                              1, wf%n_v,                    &
                              batch_b%first, batch_b%last,  &
                              1, wf%n_o,                    &
                              batch_i%first, batch_i%last)
!
!           Sort g_abji(a,b,j,i) as g_abji(a,i,j,b)
!
            call mem%alloc(g_aijb, wf%n_v, (batch_i%length), wf%n_o, (batch_b%length))
            call sort_1234_to_1432(g_abji, g_aijb, wf%n_v, (batch_b%length), wf%n_o, (batch_i%length))
!
            call mem%dealloc(g_abji, wf%n_v, (batch_b%length), wf%n_o, (batch_i%length))
!
            call mem%alloc(c_jb, wf%n_o, (batch_b%length))
!
!$omp parallel do private(j, b)
            do j = 1, wf%n_o
               do b = 1, batch_b%length
!
                  c_jb(j, b) = c_bj(b + batch_b%first - 1, j)
!
               enddo
            enddo
!$omp end parallel do
!
!           rho_a_i = rho_a_i - sum_bj g_aijb * c_jb
!
            rho_offset = wf%n_v*(batch_i%first - 1) + 1
!
            call dgemm('N', 'N',                   &
                        (wf%n_v)*(batch_i%length), &
                        1,                         &
                        (wf%n_o)*(batch_b%length), &
                        -one,                      &
                        g_aijb,                    & ! g_ai_jb
                        (wf%n_v)*(batch_i%length), &
                        c_jb,                      & ! c_jb
                        (wf%n_o)*(batch_b%length), &
                        one,                       & 
                        rho_ai(1, batch_i%first),  & ! HALLO
                        (wf%n_v)*(wf%n_o))
!
            call mem%dealloc(g_aijb, wf%n_v, (batch_i%length), wf%n_o,(batch_b%length))
            call mem%dealloc(c_jb, wf%n_o, (batch_b%length))
!
         enddo ! batch_b
      enddo ! batch_i
!
      call timer%turn_off()
!
   end subroutine jacobian_ccs_b1_ccs
!
!
end submodule jacobian_ccs
