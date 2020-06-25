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
submodule (ccs_class) jacobian_transpose_ccs_complex
!
!!
!!    Jacobian transpose submodule
!!
!!    Routines for the linear transform of trial
!!    vectors by the transpose of the Jacobian matrix
!!
!!    σ_i = A^T * b_i,
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
   module subroutine prepare_for_jacobian_transpose_ccs_complex(wf)
!!
!!    Prepare for jacobian transpose
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
   end subroutine prepare_for_jacobian_transpose_ccs_complex
!
!
   module subroutine jacobian_transpose_transformation_ccs_complex(wf, b)
!!
!!    Jacobian transpose transformation 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
!!    Calculates the transpose Jacobian transformation, i.e., the transformation
!!    by the transpose of the Jacobian matrix
!!
!!       A_mu,nu = < mu | exp(-T) [H, τ_nu] exp(T) | R >.
!!
!!    In particular,
!!
!!       sigma_mu = (b^T A)_mu = sum_ck b_ck A_ck,mu.
!!
!!    On exit, b is overwritten by sigma.
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      complex(dp), dimension(wf%n_es_amplitudes), intent(inout) :: b
!
      complex(dp), dimension(:,:), allocatable :: sigma_ai
!
      type(timings), allocatable :: timer 
!
      timer = timings('Jacobian transpose CCS', pl='normal')
      call timer%turn_on()
!
!     Allocate the transformed vector & add the terms to it
!
      call mem%alloc(sigma_ai, wf%n_v, wf%n_o)
      call zero_array_complex(sigma_ai, wf%n_t1)
!
      call wf%jacobian_transpose_ccs_a1_complex(sigma_ai, b)
      call wf%jacobian_transpose_ccs_b1_complex(sigma_ai, b)
!
!     Then overwrite the b vector with the transformed vector
!
      call zcopy((wf%n_o)*(wf%n_v), sigma_ai, 1, b, 1)
      call mem%dealloc(sigma_ai, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_transpose_transformation_ccs_complex
!
!
   module subroutine jacobian_transpose_ccs_a1_ccs_complex(wf, sigma_ai, b_ai)
!!
!!    Jacobian transpose A1 (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
!!    Calculates the A1 term,
!!
!!       sum_c b_ci F_ca - sum_k b_ak F_ik,
!!
!!    and adds it to the sigma-vector (b^T -> sigma^T = b^T A).
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      complex(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      complex(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: b_ai
!
      type(timings), allocatable :: timer 
!
      timer = timings('Jacobian transpose CCS A1', pl='verbose')
      call timer%turn_on()
!
!     Add sum_c F_ca b_ci = sum_c F_ac^T b_ci
!
      call zgemm('T','N',     &
                  wf%n_v,     &
                  wf%n_o,     &
                  wf%n_v,     &
                  one_complex,        &
                  wf%fock_ab_complex, &
                  wf%n_v,     &
                  b_ai,       &
                  wf%n_v,     &
                  one_complex,        &
                  sigma_ai,   &
                  wf%n_v)
!
!     Add - sum_k b_ak F_ik = - sum_k b_ak F_ki^T
!
      call zgemm('N','T',     &
                  wf%n_v,     &
                  wf%n_o,     &
                  wf%n_o,     &
                  -one_complex,       &
                  b_ai,       &
                  wf%n_v,     &
                  wf%fock_ij_complex, &
                  wf%n_o,     &
                  one_complex,        &
                  sigma_ai,   &
                  wf%n_v)
!
      call timer%turn_off()
!
   end subroutine jacobian_transpose_ccs_a1_ccs_complex
!
!
   module subroutine jacobian_transpose_ccs_b1_ccs_complex(wf, sigma_ai, b_ai)
!!
!!    Jacobian transpose B1 (CCS)
!!    Written by Sarai D. Folkestad, Jun 2019
!!
!!    Calculates the (CCS) B1 term of the Jacobian transpose 
!!    transfromation. 
!!
!!       B1 = sum_bj L_bjia b_bj
!!          = sum_bj (2 g_bjia b_bj - g_baij b_bj)
!!    
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      complex(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      complex(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: b_ai
!
      integer :: req0, req1_j, req1_a, req2
      integer :: current_a_batch, current_j_batch
!
      type(batching_index) batch_a, batch_j
!
      complex(dp), dimension(:,:), allocatable :: sigma_ia
      complex(dp), dimension(:,:,:,:), allocatable :: L_bjia, g_baij
!
      integer :: i, a
!
      type(timings), allocatable :: timer 
!
      timer = timings('Jacobian transpose CCS B1', pl='verbose')
      call timer%turn_on()
!
      call mem%alloc(sigma_ia, wf%n_o, wf%n_v)
      call zero_array_complex(sigma_ia, wf%n_t1)
!
      req0 = 0
!
      req1_a = max((wf%integrals%n_J)*(wf%n_v),(wf%integrals%n_J)*(wf%n_o))
      req1_j = max((wf%integrals%n_J)*(wf%n_v),(wf%integrals%n_J)*(wf%n_o))
!
      req2 = (wf%n_o)*(wf%n_v)*2
!
      batch_a = batching_index(wf%n_v)
      batch_j = batching_index(wf%n_o)
!
      call mem%batch_setup(batch_a, batch_j, req0, req1_a, req1_j, req2)
!
      do current_a_batch = 1, batch_a%num_batches
!
         call batch_a%determine_limits(current_a_batch)
!
         do current_j_batch = 1, batch_j%num_batches
!
            call batch_j%determine_limits(current_j_batch)
!
            call mem%alloc(L_bjia, wf%n_v, batch_j%length, wf%n_o, batch_a%length)
!
            call wf%get_voov_complex(L_bjia,                        &
                              1, wf%n_v,                    &
                              batch_j%first, batch_j%last,  &
                              1, wf%n_o,                    &
                              batch_a%first, batch_a%last)
!
            call zscal((wf%n_v)*(wf%n_o)*(batch_a%length)*(batch_j%length), two_complex, L_bjia, 1)
!
            call mem%alloc(g_baij, wf%n_v, batch_a%length, wf%n_o, batch_j%length)
!
            call wf%get_vvoo_complex(g_baij,                     &
                           1, wf%n_v,                    &
                           batch_a%first, batch_a%last,  &
                           1, wf%n_o,                    &
                           batch_j%first, batch_j%last)
!
            call add_1432_to_1234(-one_complex, g_baij, L_bjia, wf%n_v, batch_j%length, wf%n_o, batch_a%length)
!
            call mem%dealloc(g_baij, wf%n_v, batch_a%length, wf%n_o, batch_j%length)
!
            call zgemm('N', 'N',                   &
                        1,                         &
                        wf%n_o*(batch_a%length),   &
                        (batch_j%length)*(wf%n_v), &
                        one_complex,                       &
                        b_ai(1, batch_j%first),    & ! b_bj
                        1,                         &
                        L_bjia,                    &
                        (batch_j%length)*(wf%n_v), &
                        one_complex,                       &
                        sigma_ia(1, batch_a%first),&
                        1)
!
            call mem%dealloc(L_bjia, wf%n_v, batch_j%length, wf%n_o, batch_a%length)
!
         enddo ! batch j
      enddo ! batch a
!
!$omp parallel do private(a, i)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            sigma_ai(a,i) = sigma_ai(a,i) + sigma_ia(i,a)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(sigma_ia, wf%n_o, wf%n_v)
!
      call timer%turn_off()
!
   end subroutine jacobian_transpose_ccs_b1_ccs_complex
!
!
end submodule jacobian_transpose_ccs_complex
