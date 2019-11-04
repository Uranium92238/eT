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
submodule (ccs_class) jacobian_transpose_ccs
!
!!
!!    Jacobian transpose submodule (CCS)
!!    Set up by Andreas Skeidsvoll, Aug 2019
!!
!
   implicit none
!
!
contains
!
!
   module subroutine prepare_for_jacobian_transpose_ccs(wf)
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
      call output%printf('- No preparations for the ' // trim(wf%name_) // &
                         ' excited state equation.', pl='v', fs='(/t3,a)')
!
   end subroutine prepare_for_jacobian_transpose_ccs
!
!
   module subroutine jacobian_transpose_transformation_ccs(wf, b)
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
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: b
!
      real(dp), dimension(:,:), allocatable :: sigma_ai
!
!     Allocate the transformed vector & add the terms to it
!
      call mem%alloc(sigma_ai, wf%n_v, wf%n_o)
      call zero_array(sigma_ai, wf%n_t1)
!
      call wf%jacobian_transpose_ccs_a1(sigma_ai, b)
      call wf%jacobian_transpose_ccs_b1(sigma_ai, b)
!
!     Then overwrite the b vector with the transformed vector
!
      call dcopy((wf%n_o)*(wf%n_v), sigma_ai, 1, b, 1)
      call mem%dealloc(sigma_ai, wf%n_v, wf%n_o)
!
   end subroutine jacobian_transpose_transformation_ccs
!
!
   module subroutine jacobian_transpose_ccs_a1_ccs(wf, sigma_ai, b_ai)
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
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: b_ai
!
!     Add sum_c F_ca b_ci = sum_c F_ac^T b_ci
!
      call dgemm('T','N',     &
                  wf%n_v,     &
                  wf%n_o,     &
                  wf%n_v,     &
                  one,        &
                  wf%fock_ab, &
                  wf%n_v,     &
                  b_ai,       &
                  wf%n_v,     &
                  one,        &
                  sigma_ai,   &
                  wf%n_v)
!
!     Add - sum_k b_ak F_ik = - sum_k b_ak F_ki^T
!
      call dgemm('N','T',     &
                  wf%n_v,     &
                  wf%n_o,     &
                  wf%n_o,     &
                  -one,       &
                  b_ai,       &
                  wf%n_v,     &
                  wf%fock_ij, &
                  wf%n_o,     &
                  one,        &
                  sigma_ai,   &
                  wf%n_v)
!
   end subroutine jacobian_transpose_ccs_a1_ccs
!
!
   module subroutine jacobian_transpose_ccs_b1_ccs(wf, sigma_ai, b_ai)
!!
!!    Jacobian transpose B1 (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
!!    Calculates the B1 term,
!!
!!       sum_ck L_ckia b_ck
!!
!!    and adds it to the sigma-vector (b^T -> sigma^T = b^T A).
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: b_ai
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ckia ! g_ckia
      real(dp), dimension(:,:,:,:), allocatable :: g_caik ! g_caik
!
      real(dp), dimension(:,:,:,:), allocatable :: L_aick ! L_ckia = 2 * g_ckia - g_caik
!
      integer :: k, c, i, a
!
      integer              :: req0, req1, current_a_batch
      type(batching_index) :: batch_a
!
!     :: Construct L_aick = L_ckia
!
      call mem%alloc(g_ckia, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call wf%get_voov(g_ckia)
!
      call mem%alloc(L_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call zero_array(L_aick, (wf%n_o)**2*(wf%n_v)**2)
!
      batch_a = batching_index(wf%n_v)
!
      req0 = wf%integrals%n_J*wf%n_o**2 ! L_ik^J
      req1 = wf%n_v*wf%n_o**2 + wf%integrals%n_J*wf%n_v ! g_caik
!
      call mem%batch_setup(batch_a, req0, req1)
!
      do current_a_batch = 1, batch_a%num_batches
!
!        Set part of L_aick = L_ckia = 2 * g_ckia - g_caik for current a batch
!
         call batch_a%determine_limits(current_a_batch)
!
         call mem%alloc(g_caik, wf%n_v, batch_a%length, wf%n_o, wf%n_o)
!
         call wf%get_vvoo(g_caik,         &
                           1,             &
                           wf%n_v,        &
                           batch_a%first, &
                           batch_a%last,  &
                           1,             &
                           wf%n_o,        &
                           1,             &
                           wf%n_o)
!
!$omp parallel do private(k,c,i,a)
         do k = 1, wf%n_o
            do c = 1, wf%n_v
               do i = 1, wf%n_o
                  do a = 1, batch_a%length
!
                     L_aick(a + batch_a%first - 1,i,c,k) = two*g_ckia(c,k,i,a + batch_a%first - 1) - g_caik(c,a,i,k)
!
                  enddo
               enddo
            enddo
         enddo
!$omp end parallel do
!
         call mem%dealloc(g_caik, wf%n_v, batch_a%length, wf%n_o, wf%n_o)
!
      enddo ! End of batches over a
!
      call mem%dealloc(g_ckia, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
!     :: Add sum_ck L_ckia b_ck = sum_ck L_aick b_ck to sigma
!
      call dgemm('N','N',            &
                  (wf%n_v)*(wf%n_o), &
                  1,                 &
                  (wf%n_v)*(wf%n_o), &
                  one,               &
                  L_aick,            &
                  (wf%n_v)*(wf%n_o), &
                  b_ai,              & ! "b_ai"
                  (wf%n_v)*(wf%n_o), &
                  one,               &
                  sigma_ai,          & ! "sigma_ai"
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(L_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine jacobian_transpose_ccs_b1_ccs
!
!
end submodule jacobian_transpose_ccs