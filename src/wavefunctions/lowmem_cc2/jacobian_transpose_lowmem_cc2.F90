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
submodule (lowmem_cc2_class) jacobian_transpose
!
!!
!!    Jacobian transpose submodule
!!
!!    Routines for the linear transform of trial
!!    vectors by the Jacobian transpose matrix
!!
!!    σ_i = A^T * c_i,
!!
!!    where
!!
!!    A_μ,ν = < μ |exp(-T) [H, τ_ν] exp(T) | R >.
!!
!
   implicit none
!
!
contains
!
!
   module subroutine effective_jacobian_transpose_transformation_lowmem_cc2(wf, omega, b, sigma, &
                                                                            cvs, rm_core)
!!
!!    Effective Jacobian transpose transformation
!!    Written by Sarai Dery Folkestad, Jun 2019
!!
      use array_initialization, only: zero_array
!
      implicit none
!
      class(lowmem_cc2) :: wf
!
      real(dp), intent(in) :: omega
      real(dp), dimension(wf%n_t1), intent(in)  :: b
      real(dp), dimension(wf%n_t1), intent(out) :: sigma
!
      logical, intent(in) :: cvs, rm_core
!
      real(dp), dimension(:), allocatable :: eps_o
      real(dp), dimension(:), allocatable :: eps_v
!
      type(timings), allocatable :: timer
!
      timer = timings('Effective Jacobian transpose lowmem CC2', pl='normal')
      call timer%turn_on()
!
      if (cvs) call output%error_msg('CVS not yet implemented for lowmem CC2.')
      if (rm_core) call output%error_msg('"Remove core" not yet implemented for lowmem CC2.')
!
!     Zero the transformed vector
!
      call zero_array(sigma, wf%n_t1)
!
!     :: CCS contributions to the singles c vector ::
!
      call wf%ccs%jacobian_transpose_transformation(b, sigma)
!
!     :: CC2 contributions to the transformed singles vector ::
!
      call mem%alloc(eps_o, wf%n_o)
      call mem%alloc(eps_v, wf%n_v)
!
      call dcopy(wf%n_o, wf%orbital_energies, 1, eps_o, 1)
      call dcopy(wf%n_v, wf%orbital_energies(wf%n_o + 1), 1, eps_v, 1)
!
      call wf%jacobian_transpose_cc2_a1(sigma, b, eps_o, eps_v)
      call wf%jacobian_transpose_cc2_b1(sigma, b, eps_o, eps_v)
!
      call wf%effective_jacobian_transpose_cc2_a1(omega, sigma, b, eps_o, eps_v)
      call wf%effective_jacobian_transpose_cc2_b1(omega, sigma, b, eps_o, eps_v)
      call wf%effective_jacobian_transpose_cc2_c1(omega, sigma, b, eps_o, eps_v)
      call wf%effective_jacobian_transpose_cc2_d1(omega, sigma, b, eps_o, eps_v)
      call wf%effective_jacobian_transpose_cc2_e1(omega, sigma, b, eps_o, eps_v)
      call wf%effective_jacobian_transpose_cc2_f1(omega, sigma, b, eps_o, eps_v)
!
      call mem%dealloc(eps_o, wf%n_o)
      call mem%dealloc(eps_v, wf%n_v)
!
      call timer%turn_off()
!
   end subroutine effective_jacobian_transpose_transformation_lowmem_cc2
!
!
   module subroutine jacobian_transpose_cc2_a1_lowmem_cc2(wf, sigma_ai, b_ai, eps_o, eps_v)
!!
!!    Jacobian transpose A1 (CC2)
!!    Written by Sarai D. Folkestad, Jun 2019
!!
!!    Constructs the Jacobian transpose A1 term
!!
!!       A1 = sum_ckbj u_bjck L_iakc b_bj - t_bjck L_kcja b_bi
!!
!!    where
!!
!!       u_bjck = 2t_bjck - t_bkcj
!!
!!    and
!!
!!       t_bjck = - g_bjck/ε_bjck
!!
!!    Batching over j and k, we will construct the intermediates
!!
!!       X_ck = sum_bj u_bjck b_bj = - 2 sum_bj g_bjck b_bj / ε_bjck
!!                                       + sum_bj g_cjbk b_bj / ε_bjck
!!
!!       Y_ba = sum_bj t_bjck L_kcja = - sum_bj g_bjck L_kcja / ε_bjck
!!
      use reordering, only: add_3214_to_1234, add_3412_to_1234
      use reordering, only: add_1432_to_1234, add_2143_to_1234, add_2341_to_1234
!
      implicit none
!
      class(lowmem_cc2), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: b_ai
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
      real(dp), dimension(:,:), allocatable :: X_ck, Y_ba
      real(dp), dimension(:,:,:,:), allocatable :: g_bjck, g_iakc, u_ckbj
      real(dp), dimension(:,:,:,:), allocatable :: L_aick, L_jcka, g_kcja
!
      integer :: b, j, c, k
!
      integer :: req0, req1_i, req1_k, req2
      integer, dimension(2) :: req1
!
      integer :: current_i_batch, current_j_batch, current_k_batch
!
      type(batching_index) :: batch_i, batch_j, batch_k
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian transpose CC2 A1', pl='verbose')
      call timer%turn_on()
!
      req0 = 0
!
      req1 = wf%eri_t1%get_memory_estimate('ovov', 1, wf%n_v, 1, wf%n_v)
      req2 = 2*wf%n_v**2
!
      batch_j = batching_index(wf%n_o)
      batch_k = batching_index(wf%n_o)
!
      call mem%alloc(X_ck, wf%n_v, wf%n_o, set_zero=.true.)
!
      call mem%alloc(Y_ba, wf%n_v, wf%n_v, set_zero=.true.)
!
      call mem%batch_setup(batch_j, batch_k, req0, req1(1), req1(2), req2, &
                           tag='jacobian_transpose_cc2_a1_lowmem_cc2 1')
!
      do current_j_batch = 1, batch_j%num_batches
!
         call batch_j%determine_limits(current_j_batch)
!
         do current_k_batch = 1, batch_k%num_batches
!
            call batch_k%determine_limits(current_k_batch)
!
!           L_kcja ordered as L_jcka
!
            call mem%alloc(g_kcja, batch_k%length, wf%n_v, batch_j%length, wf%n_v)
            call wf%eri_t1%get('ovov', g_kcja,               &
                                   batch_k%first, batch_k%get_last(),  &
                                   1, wf%n_v,                    &
                                   batch_j%first, batch_j%get_last(),  &
                                   1, wf%n_v)
!
            call mem%alloc(L_jcka, batch_j%length, wf%n_v, batch_k%length, wf%n_v, set_zero=.true.)
!
            call add_3214_to_1234(two, g_kcja, L_jcka, batch_j%length, wf%n_v, batch_k%length, wf%n_v)
            call add_3412_to_1234(-one, g_kcja, L_jcka, batch_j%length, wf%n_v, batch_k%length, wf%n_v)
!
            call mem%dealloc(g_kcja, batch_k%length, wf%n_v, batch_j%length, wf%n_v)
!
            call mem%alloc(g_bjck, wf%n_v, batch_j%length, wf%n_v, batch_k%length)
!
            call wf%eri_t1%get('vovo', g_bjck, 1, wf%n_v, batch_j%first, batch_j%get_last(),  &
                                                   1, wf%n_v, batch_k%first, batch_k%get_last())
!
!           t_bjck = - g_bjck/ε_bjck
!
            do b = 1, wf%n_v
               do j = 1, batch_j%length
                  do c = 1, wf%n_v
                     do k = 1, batch_k%length
!
                        g_bjck(b, j, c, k) = -g_bjck(b, j, c, k)/(- eps_o(j + batch_j%first - 1) &
                                                                  - eps_o(k + batch_k%first - 1) &
                                                                  + eps_v(b) + eps_v(c))
!
                     enddo
                  enddo
               enddo
            enddo
!
!           Y_ba = - sum_ckj L_kcja * t_bjck
!
            call dgemm('N', 'N', &
                        wf%n_v,  &
                        wf%n_v,  &
                        (batch_k%length)*(batch_j%length)*(wf%n_v), &
                        -one,    &
                        g_bjck,  &
                        wf%n_v,  &
                        L_jcka,  &
                        (batch_k%length)*(batch_j%length)*(wf%n_v), &
                        one,     &
                        Y_ba,    &
                        wf%n_v)
!
            call mem%dealloc(L_jcka, batch_j%length, wf%n_v, batch_k%length, wf%n_v)
!
!           X_ck = sum_bj u_ckbj b_bj
!
            call mem%alloc(u_ckbj, wf%n_v, batch_k%length, wf%n_v, batch_j%length, set_zero=.true.)
!
            call add_3412_to_1234(two, g_bjck, u_ckbj, wf%n_v, batch_k%length, wf%n_v, batch_j%length)
            call add_1432_to_1234(-one, g_bjck, u_ckbj, wf%n_v, batch_k%length, wf%n_v, batch_j%length)
!
            call mem%dealloc(g_bjck, wf%n_v, batch_j%length, wf%n_v, batch_k%length)
!
            call dgemm('N', 'N',                   &
                        (batch_k%length)*(wf%n_v), &
                        1,                         &
                        (batch_j%length)*(wf%n_v), &
                        one,                       &
                        u_ckbj,                    &
                        (batch_k%length)*(wf%n_v), &
                        b_ai(1, batch_j%first),    & ! b_bj
                        (batch_j%length)*(wf%n_v), &
                        one,                       &
                        X_ck(1, batch_k%first),    &
                        wf%n_v*wf%n_o)
!
            call mem%dealloc(u_ckbj, wf%n_v, batch_k%length, wf%n_v, batch_j%length)
!
         enddo
      enddo
!
      call mem%batch_finalize()
!
!     Term 2: sigma_ai += sum_b Y_ba b_bi
!
      call dgemm('T', 'N',    &
                  wf%n_v,     &
                  wf%n_o,     &
                  wf%n_v,     &
                  one,        &
                  Y_ba,       &
                  wf%n_v,     &
                  b_ai,       & ! b_bi
                  wf%n_v,     &
                  one,        &
                  sigma_ai,   &
                  wf%n_v)
!
      call mem%dealloc(Y_ba, wf%n_v, wf%n_v)
!
!     sigma_ai + = sum_ck X_ck L_iakc
!
      req0 = 0
!
      req1_i = wf%eri_t1%n_J*wf%n_v
      req1_k = wf%eri_t1%n_J*wf%n_v
!
      req2 =  2*wf%n_v**2
!
      batch_i = batching_index(wf%n_o)
      batch_k = batching_index(wf%n_o)
!
      call mem%batch_setup(batch_i, batch_k, req0, req1_i, req1_k, req2, &
                           tag='jacobian_transpose_cc2_a1_lowmem_cc2 2')
!
      do current_i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(current_i_batch)
!
         do current_k_batch = 1, batch_k%num_batches
!
            call batch_k%determine_limits(current_k_batch)
!
!           L_iakc ordered as L_aick
!
            call mem%alloc(g_iakc, batch_i%length, wf%n_v, batch_k%length, wf%n_v)
!
            call wf%eri_t1%get('ovov', g_iakc,               &
                                   batch_i%first, batch_i%get_last(),  &
                                   1, wf%n_v,                    &
                                   batch_k%first, batch_k%get_last(),  &
                                   1, wf%n_v)
!
            call mem%alloc(L_aick, wf%n_v, batch_i%length, wf%n_v, batch_k%length, set_zero=.true.)
!
            call add_2143_to_1234(two, g_iakc, L_aick, wf%n_v, batch_i%length, wf%n_v, batch_k%length)
            call add_2341_to_1234(-one, g_iakc, L_aick, wf%n_v, batch_i%length, wf%n_v, batch_k%length)
!
            call mem%dealloc(g_iakc, batch_i%length, wf%n_v, batch_k%length, wf%n_v)
!
            call dgemm('N', 'N',                      &
                        wf%n_v*batch_i%length,        &
                        1,                            &
                        wf%n_v*batch_k%length,        &
                        one,                          &
                        L_aick,                       &
                        wf%n_v*batch_i%length,        &
                        X_ck(1,batch_k%first),        &
                        wf%n_o*wf%n_v,                &
                        one,                          &
                        sigma_ai(1, batch_i%first),   &
                        wf%n_v*wf%n_o)
!
            call mem%dealloc(L_aick, wf%n_v, batch_i%length, wf%n_v, batch_k%length)
!
         enddo
      enddo
!
      call mem%batch_finalize()
!
      call mem%dealloc(X_ck, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_transpose_cc2_a1_lowmem_cc2
!
!
   module subroutine jacobian_transpose_cc2_b1_lowmem_cc2(wf, sigma_ai, b_ai, eps_o, eps_v)
!!
!!    Jacobian transpose B1 (CC2)
!!    Written by Sarai D. Folkestad, Jun 2019
!!
!!    Constructs the Jacobian transpose A1 term
!!
!!       B1 =  - sum_ckbj t_bjck L_kcib b_aj
!!
!!    where
!!
!!       t_bjck = - g_bjck/ε_bjck
!!
!!    Batching over b and c, we will construct the intermediate
!!
!!       X_ij = sum_bck t_ckbj L_ibkc
!!
      use reordering, only: add_3412_to_1234, add_1432_to_1234
!
      implicit none
!
      class(lowmem_cc2), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: b_ai
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
      real(dp), dimension(:,:), allocatable :: X_ij
      real(dp), dimension(:,:,:,:), allocatable :: g_ckbj, L_ickb, g_ibkc
!
      integer :: b, j, c, k
!
      integer :: req0, req1_b, req1_c, req2
!
      integer :: current_b_batch, current_c_batch
!
      type(batching_index) :: batch_b, batch_c
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian transpose CC2 B1', pl='verbose')
      call timer%turn_on()
!
      req0 = 0
!
      req1_b = wf%n_o*(wf%eri_t1%n_J)
      req1_c = wf%n_o*(wf%eri_t1%n_J)
!
      req2 = 2*wf%n_o**2
!
      call mem%alloc(X_ij, wf%n_o, wf%n_o, set_zero=.true.)
!
      batch_b = batching_index(wf%n_v)
      batch_c = batching_index(wf%n_v)
!
      call mem%batch_setup(batch_b, batch_c, req0, req1_b, req1_c, req2, &
                           tag='jacobian_transpose_cc2_b1_lowmem_cc2')
!
      do current_b_batch = 1, batch_b%num_batches
!
         call batch_b%determine_limits(current_b_batch)
!
         do current_c_batch = 1, batch_c%num_batches
!
            call batch_c%determine_limits(current_c_batch)
!
!           L_ibkc ordered as L_ickb
!
            call mem%alloc(g_ibkc, wf%n_o, batch_b%length, wf%n_o, batch_c%length)
!
            call wf%eri_t1%get('ovov', g_ibkc,               &
                                   1, wf%n_o,                    &
                                   batch_b%first, batch_b%get_last(),  &
                                   1, wf%n_o,                    &
                                   batch_c%first, batch_c%get_last())
!
            call mem%alloc(L_ickb, wf%n_o, batch_c%length, wf%n_o, batch_b%length, set_zero=.true.)
!
            call add_3412_to_1234(-one, g_ibkc, L_ickb, wf%n_o, batch_c%length, wf%n_o, batch_b%length)
!
            call add_1432_to_1234(two, g_ibkc, L_ickb, wf%n_o, batch_c%length, wf%n_o, batch_b%length)
!
            call mem%dealloc(g_ibkc, wf%n_o, batch_b%length, wf%n_o, batch_c%length)
!
            call mem%alloc(g_ckbj, batch_c%length, wf%n_o, batch_b%length, wf%n_o)
!
            call wf%eri_t1%get('vovo', g_ckbj, batch_c%first, batch_c%get_last(), 1, wf%n_o, &
                                                   batch_b%first, batch_b%get_last(), 1, wf%n_o)
!
!           t_ckbj = - g_ckbj/ε_ckbj
!
            do b = 1, batch_b%length
               do j = 1, wf%n_o
                  do c = 1, batch_c%length
                     do k = 1, wf%n_o
!
                        g_ckbj(c, k, b, j) = -g_ckbj(c, k, b, j)/(  eps_v(b + batch_b%first - 1) &
                                                                  + eps_v(c + batch_c%first - 1) &
                                                                  - eps_o(j) - eps_o(k))
!
                     enddo
                  enddo
               enddo
            enddo
!
            call dgemm('N', 'N', &
                        wf%n_o,  &
                        wf%n_o,  &
                        wf%n_o*(batch_b%length)*(batch_c%length), &
                        one,     &
                        L_ickb,  &
                        wf%n_o,  &
                        g_ckbj,  &
                        wf%n_o*(batch_b%length)*(batch_c%length), &
                        one,     &
                        X_ij,    &
                        wf%n_o)
!
            call mem%dealloc(L_ickb, wf%n_o, batch_c%length, wf%n_o, batch_b%length)
            call mem%dealloc(g_ckbj, batch_c%length, wf%n_o, batch_b%length, wf%n_o)
!
         enddo
      enddo
!
      call mem%batch_finalize()
!
      call dgemm('N', 'T',    &
                  wf%n_v,     &
                  wf%n_o,     &
                  wf%n_o,     &
                  -one,       &
                  b_ai,       & ! b_aj
                  wf%n_v,     &
                  X_ij,       &
                  wf%n_o,     &
                  one,        &
                  sigma_ai,   &
                  wf%n_v)
!
      call mem%dealloc(X_ij, wf%n_o, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_transpose_cc2_b1_lowmem_cc2
!
!
   module subroutine effective_jacobian_transpose_cc2_a1_lowmem_cc2(wf, omega, sigma_ai, b_ai, eps_o, eps_v)
!!
!!    Effective jacobian transpose CC2 A1
!!    Written by Sarai D. Folkestad, Jul 2019
!!
!!    sigma_ai =+ g_bjca [P_bj,ci( 2 F_jb b_ci - F_ib b_cj)] 1/(ω - ε_cibj)
!!
!!    The term is calculated in batches over the b and c.
!!
!
      implicit none
!
      class(lowmem_cc2), intent(inout) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: b_ai
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
      integer :: req0, req1_b, req1_c, req2
!
      integer :: current_b_batch, current_c_batch
!
      type(batching_index) :: batch_b, batch_c
!
      integer :: i, j, b, c
!
      real(dp), dimension(:,:,:,:), allocatable :: X_bjci, g_bjca
!
      type(timings), allocatable :: timer
!
      timer = timings('Effective jacobian transpose CC2 A1', pl='verbose')
      call timer%turn_on()
!
      req0 = 0
!
      req1_b = (wf%eri_t1%n_J)*wf%n_o
      req1_c = (wf%eri_t1%n_J)*wf%n_v
!
      req2 = 2*(wf%n_v)*(wf%n_o)
!
      batch_b = batching_index(wf%n_v)
      batch_c = batching_index(wf%n_v)
!
      call mem%batch_setup(batch_b, batch_c, req0, req1_b, req1_c, req2, &
                           tag='effective_jacobian_transpose_cc2_a1_lowmem_cc2')
!
      do current_b_batch = 1, batch_b%num_batches
!
         call batch_b%determine_limits(current_b_batch)
!
         do current_c_batch = 1, batch_c%num_batches
!
            call batch_c%determine_limits(current_c_batch)
!
            call mem%alloc(g_bjca, batch_b%length, wf%n_o, batch_c%length, wf%n_v)
!
            call wf%eri_t1%get('vovv', g_bjca, batch_b%first, batch_b%get_last(), 1, wf%n_o, &
                                                   batch_c%first, batch_c%get_last(), 1, wf%n_v)
!
            call mem%alloc(X_bjci, batch_b%length, wf%n_o, batch_c%length, wf%n_o)
!
!$omp parallel do private(i, c, j, b)
            do i = 1, wf%n_o
               do c = 1, batch_c%length
                  do j = 1, wf%n_o
                     do b = 1, batch_b%length
!
                        X_bjci(b, j, c, i) = (two*wf%fock_ia(j, b + batch_b%first - 1)*b_ai(c + batch_c%first - 1, i) &
                                           +  two*wf%fock_ia(i, c + batch_c%first - 1)*b_ai(b + batch_b%first - 1, j) &
                                           -  wf%fock_ia(i, b + batch_b%first - 1)*b_ai(c + batch_c%first - 1, j) &
                                           -  wf%fock_ia(j, c + batch_c%first - 1)*b_ai(b + batch_b%first - 1, i))/ &
                                             (omega - eps_v(b + batch_b%first - 1) - eps_v(c + batch_c%first - 1) &
                                             + eps_o(i) + eps_o(j))
!
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
            call dgemm('T', 'N',                                     &
                        wf%n_v,                                      &
                        wf%n_o,                                      &
                        (batch_b%length)*(batch_c%length)*(wf%n_o),  &
                        one,                                         &
                        g_bjca,                                      &
                        (batch_b%length)*(batch_c%length)*(wf%n_o),  &
                        X_bjci,                                      &
                        (batch_b%length)*(batch_c%length)*(wf%n_o),  &
                        one,                                         &
                        sigma_ai,                                    &
                        wf%n_v)
!
            call mem%dealloc(X_bjci, batch_b%length, wf%n_o, batch_c%length, wf%n_o)
            call mem%dealloc(g_bjca, batch_b%length, wf%n_o, batch_c%length, wf%n_v)
!
         enddo
      enddo
!
      call mem%batch_finalize()
!
      call timer%turn_off()
!
   end subroutine effective_jacobian_transpose_cc2_a1_lowmem_cc2
!
!
   module subroutine effective_jacobian_transpose_cc2_b1_lowmem_cc2(wf, omega, sigma_ai, b_ai, eps_o, eps_v)
!!
!!    Effective jacobian transpose CC2 B1
!!    Written by Sarai D. Folkestad, Jul 2019
!!
!!    sigma_ai =+ g_bjca [P_bj,ci L_dcjb b_di] 1/(ω - ε_cibj)
!!
!!    The term is calculated in batches over the b and c.
!!
      use reordering, only: sort_1234_to_4321
      use reordering, only: add_2341_to_1234, add_2143_to_1234, add_4123_to_1234
!
      implicit none
!
      class(lowmem_cc2), intent(inout) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: b_ai
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
      integer :: req0, req1_b, req1_c, req2
!
      integer :: current_b_batch, current_c_batch
!
      type(batching_index) :: batch_b, batch_c
!
      integer :: i, j, b, c
!
      real(dp), dimension(:,:,:,:), allocatable :: g_bjca, Y_bjci, X_icjb, X_jbic, g_dcjb, g_dbic
!
      type(timings), allocatable :: timer
!
      timer = timings('Effective jacobian transpose CC2 B1', pl='verbose')
      call timer%turn_on()
!
      req0 = 0
!
      req1_b = max((wf%eri_t1%n_J)*wf%n_v, (wf%eri_t1%n_J)*wf%n_o)
      req1_c = max((wf%eri_t1%n_J)*wf%n_v, (wf%eri_t1%n_J)*wf%n_o)
!
      req2 = 3*(wf%n_v)*(wf%n_o)
!
      batch_b = batching_index(wf%n_v)
      batch_c = batching_index(wf%n_v)
!
      call mem%batch_setup(batch_b, batch_c, req0, req1_b, req1_c, req2, &
                           tag='effective_jacobian_transpose_cc2_b1_lowmem_cc2')
!
      do current_b_batch = 1, batch_b%num_batches
!
         call batch_b%determine_limits(current_b_batch)
!
         do current_c_batch = 1, batch_c%num_batches
!
            call batch_c%determine_limits(current_c_batch)
!
            call mem%alloc(g_dcjb, wf%n_v, batch_c%length, wf%n_o, batch_b%length)
            call wf%eri_t1%get('vvov', g_dcjb, 1, wf%n_v, batch_c%first, batch_c%get_last(), &
                                                   1, wf%n_o, batch_b%first, batch_b%get_last())
!
            call mem%alloc(X_icjb, wf%n_o, batch_c%length, wf%n_o, batch_b%length)
!
            call dgemm('T', 'N',                                     &
                        wf%n_o,                                      &
                        (batch_c%length)*(batch_b%length)*(wf%n_o),  &
                        wf%n_v,                                      &
                        one,                                         &
                        b_ai,                                        & ! b_di
                        wf%n_v,                                      &
                        g_dcjb,                                      &
                        wf%n_v,                                      &
                        zero,                                        &
                        X_icjb,                                      &
                        wf%n_o)
!
            call mem%dealloc(g_dcjb, wf%n_v, batch_c%length, wf%n_o, batch_b%length)
!
            call mem%alloc(Y_bjci, batch_b%length, wf%n_o, batch_c%length, wf%n_o)
!
            call sort_1234_to_4321(X_icjb, Y_bjci, wf%n_o, batch_c%length, wf%n_o, batch_b%length)
!
            call dscal((batch_c%length)*(wf%n_o**2)*(batch_b%length), two, Y_bjci, 1)
!
            call add_2341_to_1234(-one, X_icjb, Y_bjci, batch_b%length, wf%n_o, batch_c%length, wf%n_o)
!
            call mem%dealloc(X_icjb, wf%n_o, batch_c%length, wf%n_o, batch_b%length)
!
            call mem%alloc(g_dbic, wf%n_v, batch_b%length, wf%n_o, batch_c%length)
            call wf%eri_t1%get('vvov', g_dbic, 1, wf%n_v, batch_b%first, batch_b%get_last(), &
                                                   1, wf%n_o, batch_c%first, batch_c%get_last())
!
            call mem%alloc(X_jbic, wf%n_o, batch_b%length, wf%n_o, batch_c%length)
!
            call dgemm('T', 'N',                                     &
                        wf%n_o,                                      &
                        (batch_b%length)*(wf%n_o)*(batch_c%length),  &
                        wf%n_v,                                      &
                        one,                                         &
                        b_ai,                                        & ! b_dj
                        wf%n_v,                                      &
                        g_dbic,                                      &
                        wf%n_v,                                      &
                        zero,                                        &
                        X_jbic,                                      &
                        wf%n_o)
!
            call mem%dealloc(g_dbic, wf%n_v, batch_b%length, wf%n_o, batch_c%length)
!
            call add_2143_to_1234(two, X_jbic, Y_bjci, batch_b%length, wf%n_o, batch_c%length, wf%n_o)
            call add_4123_to_1234(-one, X_jbic, Y_bjci, batch_b%length, wf%n_o, batch_c%length, wf%n_o)
!
            call mem%dealloc(X_jbic, wf%n_o, batch_b%length, wf%n_o, batch_c%length)
!
!
!$omp parallel do private(i, c, j, b)
            do i = 1, wf%n_o
               do c = 1, batch_c%length
                  do j = 1, wf%n_o
                     do b = 1, batch_b%length
!
                        Y_bjci(b, j, c, i) = Y_bjci(b, j, c, i)/ &
                                             (omega - (eps_v(b + batch_b%first - 1) + eps_v(c + batch_c%first - 1) &
                                             - eps_o(i) - eps_o(j)))
!
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
            call mem%alloc(g_bjca, batch_b%length, wf%n_o, batch_c%length, wf%n_v)
!
            call wf%eri_t1%get('vovv', g_bjca, batch_b%first, batch_b%get_last(), 1, wf%n_o, &
                                                   batch_c%first, batch_c%get_last(), 1, wf%n_v)
!
            call dgemm('T', 'N',                                    &
                       wf%n_v,                                      &
                       wf%n_o,                                      &
                       (batch_b%length)*(batch_c%length)*(wf%n_o),  &
                       one,                                         &
                       g_bjca,                                      &
                       (batch_b%length)*(batch_c%length)*(wf%n_o),  &
                       Y_bjci,                                      &
                       (batch_b%length)*(batch_c%length)*(wf%n_o),  &
                       one,                                         &
                       sigma_ai,                                    &
                       wf%n_v)
!
            call mem%dealloc(g_bjca, batch_b%length, wf%n_o, batch_c%length, wf%n_v)
            call mem%dealloc(Y_bjci, batch_b%length, wf%n_o, batch_c%length, wf%n_o)
!
         enddo
      enddo
!
      call mem%batch_finalize()
!
      call timer%turn_off()
!
   end subroutine effective_jacobian_transpose_cc2_b1_lowmem_cc2
!
!
   module subroutine effective_jacobian_transpose_cc2_c1_lowmem_cc2(wf, omega, sigma_ai, b_ai, eps_o, eps_v)
!!
!!    Effective jacobian transpose CC2 C1
!!    Written by Sarai D. Folkestad, Jul 2019
!!
!!    sigma_ai =- g_bjca [P_bj,ci L_ikjb b_ck] 1/(ω - ε_cibj)
!!
!!    The term is calculated in batches over the k, b, c
!!
      use array_initialization, only: copy_and_scale
      use reordering, only: add_3214_to_1234, add_2143_to_1234, add_4321_to_1234
!
      implicit none
!
      class(lowmem_cc2), intent(inout) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: b_ai
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
      real(dp), dimension(:,:,:,:), allocatable :: g_bjca, X_bjci, Y_jbic, Y_icjb
      real(dp), dimension(:,:,:,:), allocatable :: g_icjk, g_jbik, L_icjk, L_jbik
!
      type(batching_index) :: batch_k, batch_b, batch_c
!
      integer :: current_k_batch, current_c_batch, current_b_batch
      integer :: req0, req1_k, req1_c, req1_b, req2_kb, req2_kc, req2_bc, req3
!
      integer :: i, j, c, b
!
      type(timings), allocatable :: timer
!
      timer = timings('Effective jacobian transpose CC2 C1', pl='verbose')
      call timer%turn_on()
!
      req0 = 0
!
      req1_k = (wf%eri_t1%n_J)*(wf%n_o)
      req1_c = max((wf%eri_t1%n_J)*(wf%n_v),(wf%eri_t1%n_J)*(wf%n_o))
      req1_b = max((wf%eri_t1%n_J)*(wf%n_v),(wf%eri_t1%n_J)*(wf%n_o))
!
      req2_kc = 2*(wf%n_o)*(wf%n_v)
      req2_kb = 2*(wf%n_o)*(wf%n_v)
      req2_bc = max(2*(wf%n_o**2),(wf%n_o**2) + (wf%n_o)*(wf%n_v))
!
      req3 = 0
!
      batch_k = batching_index(wf%n_o)
      batch_b = batching_index(wf%n_v)
      batch_c = batching_index(wf%n_v)
!
      call mem%batch_setup(batch_k, batch_b, batch_c, req0, req1_k, req1_b, req1_c, &
                           req2_kb, req2_kc, req2_bc, req3, &
                           tag='effective_jacobian_transpose_cc2_c1_lowmem_cc2')
!
      do current_b_batch = 1, batch_b%num_batches
!
         call batch_b%determine_limits(current_b_batch)
!
         do current_c_batch = 1, batch_c%num_batches
!
            call batch_c%determine_limits(current_c_batch)
!
            call mem%alloc(X_bjci, batch_b%length, wf%n_o, batch_c%length, wf%n_o, set_zero=.true.)
!
            do current_k_batch = 1, batch_k%num_batches
!
               call batch_k%determine_limits(current_k_batch)
!
!              Construct L_jbik and contract with b_ck
!
               call mem%alloc(g_jbik, wf%n_o, batch_b%length, wf%n_o, batch_k%length)
!
               call wf%eri_t1%get('ovoo', g_jbik,                         &
                                      1, wf%n_o, batch_b%first, batch_b%get_last(), &
                                      1, wf%n_o, batch_k%first, batch_k%get_last())
!
               call mem%alloc(L_jbik, wf%n_o, batch_b%length, wf%n_o, batch_k%length)
!
               call copy_and_scale(two, g_jbik, L_jbik, (wf%n_o**2)*(batch_b%length)*(batch_k%length))
               call add_3214_to_1234(-one, g_jbik, L_jbik, wf%n_o, batch_b%length, wf%n_o, batch_k%length)
!
               call mem%dealloc(g_jbik, wf%n_o, batch_b%length, wf%n_o, batch_k%length)
!
               call mem%alloc(Y_jbic, wf%n_o, batch_b%length, wf%n_o, batch_c%length)
!
               call dgemm('N', 'T',                            &
                           (wf%n_o**2)*(batch_b%length),       &
                           batch_c%length,                     &
                           batch_k%length,                     &
                           one,                                &
                           L_jbik,                             &
                           (wf%n_o**2)*(batch_b%length),       &
                           b_ai(batch_c%first, batch_k%first), &
                           wf%n_v,                             &
                           zero,                               &
                           Y_jbic,                             &
                           (wf%n_o**2)*(batch_b%length))
!
               call mem%dealloc(L_jbik, wf%n_o, batch_b%length, wf%n_o, batch_k%length)
!
               call add_2143_to_1234(one, Y_jbic, X_bjci, batch_b%length, wf%n_o, batch_c%length, wf%n_o)
!
               call mem%dealloc(Y_jbic, wf%n_o, batch_b%length, wf%n_o, batch_c%length)
!
!              Construct L_icjk and contract with b_bk
!
               call mem%alloc(g_icjk, wf%n_o, batch_c%length, wf%n_o, batch_k%length)
!
               call wf%eri_t1%get('ovoo', g_icjk,                         &
                                      1, wf%n_o, batch_c%first, batch_c%get_last(), &
                                      1, wf%n_o, batch_k%first, batch_k%get_last())
!
               call mem%alloc(L_icjk, wf%n_o, batch_c%length, wf%n_o, batch_k%length)
!
               call copy_and_scale(two, g_icjk, L_icjk, (wf%n_o**2)*(batch_c%length)*(batch_k%length))
               call add_3214_to_1234(-one, g_icjk, L_icjk, wf%n_o, batch_c%length, wf%n_o, batch_k%length)
!
               call mem%dealloc(g_icjk, wf%n_o, batch_c%length, wf%n_o, batch_k%length)
!
               call mem%alloc(Y_icjb, wf%n_o, batch_c%length, wf%n_o, batch_b%length)
!
               call dgemm('N', 'T',                            &
                           (wf%n_o**2)*(batch_c%length),       &
                           batch_b%length,                     &
                           batch_k%length,                     &
                           one,                                &
                           L_icjk,                             &
                           (wf%n_o**2)*(batch_c%length),       &
                           b_ai(batch_b%first, batch_k%first), &
                           wf%n_v,                             &
                           zero,                               &
                           Y_icjb,                             &
                           (wf%n_o**2)*(batch_c%length))
!
               call mem%dealloc(L_icjk, wf%n_o, batch_c%length, wf%n_o, batch_k%length)
!
               call add_4321_to_1234(one, Y_icjb, X_bjci, batch_b%length, wf%n_o, batch_c%length, wf%n_o)
!
               call mem%dealloc(Y_icjb, wf%n_o, batch_c%length, wf%n_o, batch_b%length)
!
            enddo
!
!$omp parallel do private(i, c, j, b)
            do i = 1, wf%n_o
               do c = 1, batch_c%length
                  do j = 1, wf%n_o
                     do b = 1, batch_b%length
!
                        X_bjci(b, j, c, i) = X_bjci(b, j, c, i)/ &
                                             (omega - (eps_v(b + batch_b%first - 1) + eps_v(c + batch_c%first - 1) &
                                             - eps_o(i) - eps_o(j)))
!
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
!           Construct g_bjca and contract with X_bjci
!
            call mem%alloc(g_bjca, batch_b%length, wf%n_o, batch_c%length, wf%n_v)
!
            call wf%eri_t1%get('vovv', g_bjca, batch_b%first, batch_b%get_last(), 1, wf%n_o, &
                                                   batch_c%first, batch_c%get_last(), 1, wf%n_v)
!
            call dgemm('T', 'N',                                     &
                        wf%n_v,                                      &
                        wf%n_o,                                      &
                        (batch_b%length)*(batch_c%length)*(wf%n_o),  &
                        -one,                                        &
                        g_bjca,                                      &
                        (batch_b%length)*(batch_c%length)*(wf%n_o),  &
                        X_bjci,                                      &
                        (batch_b%length)*(batch_c%length)*(wf%n_o),  &
                        one,                                         &
                        sigma_ai,                                    &
                        wf%n_v)
!
            call mem%dealloc(X_bjci, batch_b%length, wf%n_o, batch_c%length, wf%n_o)
            call mem%dealloc(g_bjca, batch_b%length, wf%n_o, batch_c%length, wf%n_v)
!
         enddo
      enddo
!
      call mem%batch_finalize()
!
      call timer%turn_off()
!
   end subroutine effective_jacobian_transpose_cc2_c1_lowmem_cc2
!
!
   module subroutine effective_jacobian_transpose_cc2_d1_lowmem_cc2(wf, omega, sigma_ai, b_ai, eps_o, eps_v)
!!
!!    Effective jacobian transpose CC2 D1
!!    Written by Sarai D. Folkestad, Jul 2019
!!
!!    sigma_ai =- g_bjik [P_bj,ak( 2 F_jb b_ak - F_kb b_aj)] 1/(ω - ε_bjak)
!!
!!    The term is calculated in batches over the k and j.
!!
!
      implicit none
!
      class(lowmem_cc2), intent(inout) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: b_ai
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
      integer :: req0, req1_j, req1_k, req2
!
      integer :: current_j_batch, current_k_batch
!
      type(batching_index) :: batch_k, batch_j
!
      integer :: k, j, b, a
!
      real(dp), dimension(:,:,:,:), allocatable :: X_akbj, g_ikbj
!
      type(timings), allocatable :: timer
!
      timer = timings('Effective jacobian transpose CC2 D1', pl='verbose')
      call timer%turn_on()
!
      req0 = 0
!
      req1_j = wf%eri_t1%n_J*wf%n_v
      req1_k = wf%eri_t1%n_J*wf%n_o
!
      req2 = 2*(wf%n_v)*(wf%n_o)
!
      batch_j = batching_index(wf%n_o)
      batch_k = batching_index(wf%n_o)
!
      call mem%batch_setup(batch_j, batch_k, req0, req1_j, req1_k, req2, &
                           tag='effective_jacobian_transpose_cc2_d1_lowmem_cc2')
!
      do current_j_batch = 1, batch_j%num_batches
!
         call batch_j%determine_limits(current_j_batch)
!
         do current_k_batch = 1, batch_k%num_batches
!
            call batch_k%determine_limits(current_k_batch)
!
            call mem%alloc(g_ikbj, wf%n_o, batch_k%length, wf%n_v, batch_j%length)
!
            call wf%eri_t1%get('oovo', g_ikbj, 1, wf%n_o, batch_k%first, batch_k%get_last(), &
                                                   1, wf%n_v, batch_j%first, batch_j%get_last())
!
            call mem%alloc(X_akbj, wf%n_v, batch_k%length, wf%n_v, batch_j%length)
!
!$omp parallel do private(a, k, b, j)
            do j = 1, batch_j%length
               do b = 1, wf%n_v
                  do k = 1, batch_k%length
                     do a = 1, wf%n_v
!
                        X_akbj(a, k, b, j) = (two*wf%fock_ia(j + batch_j%first - 1, b)*b_ai(a, k + batch_k%first - 1) &
                                           +  two*wf%fock_ia(k + batch_k%first - 1, a)*b_ai(b, j + batch_j%first - 1) &
                                           -  wf%fock_ia(k + batch_k%first - 1, b)*b_ai(a, j + batch_j%first - 1) &
                                           -  wf%fock_ia(j + batch_j%first - 1, a)*b_ai(b, k + batch_k%first - 1))/ &
                                             (omega - eps_v(b) - eps_v(a) &
                                             + eps_o(k + batch_k%first - 1) + eps_o(j+ batch_j%first - 1))
!
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
            call dgemm('N', 'T',                                     &
                        wf%n_v,                                      &
                        wf%n_o,                                      &
                        (batch_k%length)*(batch_j%length)*(wf%n_v),  &
                        -one,                                        &
                        X_akbj,                                      &
                        wf%n_v,                                      &
                        g_ikbj,                                      &
                        wf%n_o,                                      &
                        one,                                         &
                        sigma_ai,                                    &
                        wf%n_v)
!
            call mem%dealloc(X_akbj, wf%n_v, batch_k%length, wf%n_v, batch_j%length)
            call mem%dealloc(g_ikbj, wf%n_o, batch_k%length, wf%n_v, batch_j%length)
!
         enddo
      enddo
!
      call mem%batch_finalize()
!
      call timer%turn_off()
!
   end subroutine effective_jacobian_transpose_cc2_d1_lowmem_cc2
!
!
   module subroutine effective_jacobian_transpose_cc2_e1_lowmem_cc2(wf, omega, sigma_ai, b_ai, eps_o, eps_v)
!!
!!    Effective jacobian transpose CC2 E1
!!    Written by Sarai D. Folkestad, Jul 2019
!!
!!    sigma_ai =+ g_bjik [P_bj,ak( L_jbkl b_al)] 1/(ω - ε_bjak)
!!
!!    The term is calculated in batches over the k and j.
!!
      use reordering, only: add_4321_to_1234, add_4123_to_1234
      use reordering, only: add_2143_to_1234, add_2341_to_1234
!
      implicit none
!
      class(lowmem_cc2), intent(inout) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: b_ai
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
      integer :: req0, req1_j, req1_k, req2
!
      integer :: current_j_batch, current_k_batch
!
      type(batching_index) :: batch_k, batch_j
!
      integer :: k, j, b, a
!
      real(dp), dimension(:,:,:,:), allocatable :: Y_akbj, X_jbka, X_kajb, g_jbkl, g_kajl, g_ikbj
!
      type(timings), allocatable :: timer
!
      timer = timings('Effective jacobian transpose CC2 E1', pl='verbose')
      call timer%turn_on()
!
      req0 = 0
!
      req1_j = (wf%eri_t1%n_J)*wf%n_v
      req1_k = (wf%eri_t1%n_J)*wf%n_v
!
      req2 = 2*(wf%n_v)*(wf%n_o)
!
      batch_j = batching_index(wf%n_o)
      batch_k = batching_index(wf%n_o)
!
      call mem%batch_setup(batch_j, batch_k, req0, req1_j, req1_k, req2, &
                           tag='effective_jacobian_transpose_cc2_e1_lowmem_cc2')
!
      do current_j_batch = 1, batch_j%num_batches
!
         call batch_j%determine_limits(current_j_batch)
!
         do current_k_batch = 1, batch_k%num_batches
!
            call batch_k%determine_limits(current_k_batch)
!
            call mem%alloc(g_jbkl, batch_j%length, wf%n_v, batch_k%length, wf%n_o)
!
            call wf%eri_t1%get('ovoo', g_jbkl, batch_j%first, batch_j%get_last(), 1, wf%n_v, &
                                                   batch_k%first, batch_k%get_last(), 1, wf%n_o)
!
            call mem%alloc(X_jbka, batch_j%length, wf%n_v, batch_k%length, wf%n_v)
!
            call dgemm('N', 'T',                                  &
                        wf%n_v*(batch_k%length)*(batch_j%length), &
                        wf%n_v,                                   &
                        wf%n_o,                                   &
                        one,                                      &
                        g_jbkl,                                   &
                        wf%n_v*(batch_k%length)*(batch_j%length), &
                        b_ai,                                     & ! b_al
                        wf%n_v,                                   &
                        zero,                                     &
                        X_jbka,                                   &
                        wf%n_v*(batch_k%length)*(batch_j%length))
!
            call mem%dealloc(g_jbkl, batch_j%length, wf%n_v, batch_k%length, wf%n_o)
!
            call mem%alloc(Y_akbj, wf%n_v, batch_k%length, wf%n_v, batch_j%length, set_zero=.true.)
!
            call add_4321_to_1234(two, X_jbka, Y_akbj, wf%n_v, batch_k%length, wf%n_v, batch_j%length)
            call add_4123_to_1234(-one, X_jbka, Y_akbj, wf%n_v, batch_k%length, wf%n_v, batch_j%length)
!
            call mem%dealloc(X_jbka, batch_j%length, wf%n_v, batch_k%length, wf%n_v)
!
            call mem%alloc(g_kajl, batch_k%length, wf%n_v, batch_j%length, wf%n_o)
!
            call wf%eri_t1%get('ovoo', g_kajl, batch_k%first, batch_k%get_last(), 1, wf%n_v, &
                                                   batch_j%first, batch_j%get_last(), 1, wf%n_o)
!
            call mem%alloc(X_kajb, batch_k%length, wf%n_v, batch_j%length, wf%n_v)
!
            call dgemm('N', 'T',                                  &
                        wf%n_v*(batch_k%length)*(batch_j%length), &
                        wf%n_v,                                   &
                        wf%n_o,                                   &
                        one,                                      &
                        g_kajl,                                   &
                        wf%n_v*(batch_k%length)*(batch_j%length), &
                        b_ai,                                     & ! b_bl
                        wf%n_v,                                   &
                        zero,                                     &
                        X_kajb,                                   &
                        wf%n_v*(batch_k%length)*(batch_j%length))
!
            call mem%dealloc(g_kajl, batch_k%length, wf%n_v, batch_j%length, wf%n_o)
!
            call add_2143_to_1234(two, X_kajb, Y_akbj, wf%n_v, batch_k%length, wf%n_v, batch_j%length)
            call add_2341_to_1234(-one, X_kajb, Y_akbj, wf%n_v, batch_k%length, wf%n_v, batch_j%length)
!
            call mem%dealloc(X_kajb, batch_k%length, wf%n_v, batch_j%length, wf%n_v)
!
!$omp parallel do private(a, k, b, j)
            do j = 1, batch_j%length
               do b = 1, wf%n_v
                  do k = 1, batch_k%length
                     do a = 1, wf%n_v
!
                        Y_akbj(a, k, b, j) = Y_akbj(a, k, b, j)/ (omega - eps_v(b) - eps_v(a) &
                                             + eps_o(k + batch_k%first - 1) + eps_o(j+ batch_j%first - 1))
!
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
            call mem%alloc(g_ikbj, wf%n_o, batch_k%length, wf%n_v, batch_j%length)
!
            call wf%eri_t1%get('oovo', g_ikbj, 1, wf%n_o, batch_k%first, batch_k%get_last(),  &
                                                   1, wf%n_v, batch_j%first, batch_j%get_last())
!
            call dgemm('N', 'T',                                     &
                        wf%n_v,                                      &
                        wf%n_o,                                      &
                        (batch_k%length)*(batch_j%length)*(wf%n_v),  &
                        one,                                         &
                        Y_akbj,                                      &
                        wf%n_v,                                      &
                        g_ikbj,                                      &
                        wf%n_o,                                      &
                        one,                                         &
                        sigma_ai,                                    &
                        wf%n_v)
!
            call mem%dealloc(Y_akbj, wf%n_v, batch_k%length, wf%n_v, batch_j%length)
            call mem%dealloc(g_ikbj, wf%n_o, batch_k%length, wf%n_v, batch_j%length)
!
         enddo
      enddo
!
      call mem%batch_finalize()
!
      call timer%turn_off()
!
   end subroutine effective_jacobian_transpose_cc2_e1_lowmem_cc2
!
!
   module subroutine effective_jacobian_transpose_cc2_f1_lowmem_cc2(wf, omega, sigma_ai, b_ai, eps_o, eps_v)
!!
!!    Effective jacobian transpose CC2 F1
!!    Written by Sarai D. Folkestad, Jul 2019
!!
!!    sigma_ai =- g_bjik [P_bj,ak( L_cajb b_ck)] 1/(ω - ε_bjak)
!!
!!    The term is calculated in batches over the k, j and c indices.
!!
      use array_initialization, only: copy_and_scale
      use reordering, only: add_1432_to_1234, add_2143_to_1234, add_4321_to_1234
!
      implicit none
!
      class(lowmem_cc2), intent(inout) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: b_ai
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
      type(batching_index) :: batch_k, batch_j, batch_c
!
      integer :: current_k_batch, current_c_batch, current_j_batch
      integer :: req0, req1_k, req1_c, req1_j, req2_kj, req2_kc, req2_jc, req3
!
      integer :: a, k, b, j
!
      real(dp), dimension(:,:,:,:), allocatable :: X_akbj, Y_kajb, Y_jbka
      real(dp), dimension(:,:,:,:), allocatable :: L_cajb, L_cbka, g_cajb, g_cbka, g_ikbj
!
      type(timings), allocatable :: timer
!
      timer = timings('Effective jacobian transpose CC2 F1', pl='verbose')
      call timer%turn_on()
!
      req0 = 0
!
      req1_k = wf%eri_t1%n_J*wf%n_v
      req1_c = wf%eri_t1%n_J*wf%n_v
      req1_j = wf%eri_t1%n_J*wf%n_v
!
      req2_kc = 2*(wf%n_v**2)
      req2_kj = 2*(wf%n_v**2)
      req2_jc = 2*(wf%n_v**2)
!
      req3 = 0
!
      batch_k = batching_index(wf%n_o)
      batch_j = batching_index(wf%n_o)
      batch_c = batching_index(wf%n_v)
!
      call mem%batch_setup(batch_k, batch_j, batch_c, req0, req1_k, req1_j, req1_c, &
                           req2_kj, req2_kj, req2_jc, req3, &
                           tag='effective_jacobian_transpose_cc2_f1_lowmem_cc2')
!
      do current_k_batch = 1, batch_k%num_batches
!
         call batch_k%determine_limits(current_k_batch)
!
         do current_j_batch = 1, batch_j%num_batches
!
            call batch_j%determine_limits(current_j_batch)
!
            call mem%alloc(X_akbj, wf%n_v, batch_k%length, wf%n_v, batch_j%length, set_zero=.true.)
!
            do current_c_batch = 1, batch_c%num_batches
!
               call batch_c%determine_limits(current_c_batch)
!
!              Construct L_cajb and contract with b_ck
!
               call mem%alloc(g_cajb, batch_c%length, wf%n_v, batch_j%length, wf%n_v)
               call wf%eri_t1%get('vvov', g_cajb, batch_c%first, batch_c%get_last(), 1, wf%n_v, &
                                                      batch_j%first, batch_j%get_last(), 1, wf%n_v)
!
               call mem%alloc(L_cajb, batch_c%length, wf%n_v, batch_j%length, wf%n_v)
!
               call copy_and_scale(two, g_cajb, L_cajb, (wf%n_v**2)*batch_c%length*batch_j%length)
               call add_1432_to_1234(-one, g_cajb, L_cajb, batch_c%length, wf%n_v, batch_j%length, wf%n_v)
!
               call mem%dealloc(g_cajb, batch_c%length, wf%n_v, batch_j%length, wf%n_v)
!
               call mem%alloc(Y_kajb, batch_k%length, wf%n_v, batch_j%length, wf%n_v)
!
               call dgemm('T', 'N', &
                           batch_k%length,                     &
                           (wf%n_v**2)*batch_j%length,         &
                           batch_c%length,                     &
                           one,                                &
                           b_ai(batch_c%first, batch_k%first), & ! b_ck
                           wf%n_v,                             &
                           L_cajb,                             &
                           batch_c%length,                     &
                           zero,                               &
                           Y_kajb,                             &
                           batch_k%length)
!
               call mem%dealloc(L_cajb, batch_c%length, wf%n_v, batch_j%length, wf%n_v)
!
               call add_2143_to_1234(one, Y_kajb, X_akbj, wf%n_v, batch_k%length, wf%n_v, batch_j%length)
!
               call mem%dealloc(Y_kajb, batch_k%length, wf%n_v, batch_j%length, wf%n_v)
!
!              Construct L_cbka and contract with c_cj
!
               call mem%alloc(g_cbka, batch_c%length, wf%n_v, batch_k%length, wf%n_v)
!
               call wf%eri_t1%get('vvov', g_cbka, batch_c%first, batch_c%get_last(), 1, wf%n_v, &
                                                      batch_k%first, batch_k%get_last(), 1, wf%n_v)
!
               call mem%alloc(L_cbka, batch_c%length, wf%n_v, batch_k%length, wf%n_v)
!
               call copy_and_scale(two, g_cbka, L_cbka, (wf%n_v**2)*batch_c%length*batch_k%length)
               call add_1432_to_1234(-one, g_cbka, L_cbka, batch_c%length, wf%n_v, batch_k%length, wf%n_v)
!
               call mem%dealloc(g_cbka, batch_c%length, wf%n_v, batch_k%length, wf%n_v)
!
               call mem%alloc(Y_jbka, batch_j%length, wf%n_v, batch_k%length, wf%n_v)
!
               call dgemm('T', 'N', &
                           batch_j%length,                     &
                           (wf%n_v**2)*batch_k%length,         &
                           batch_c%length,                     &
                           one,                                &
                           b_ai(batch_c%first, batch_j%first), & ! b_cj
                           wf%n_v,                             &
                           L_cbka,                             &
                           batch_c%length,                     &
                           zero,                               &
                           Y_jbka,                             &
                           batch_j%length)
!
               call mem%dealloc(L_cbka, batch_c%length, wf%n_v, batch_k%length, wf%n_v)
!
               call add_4321_to_1234(one, Y_jbka, X_akbj, wf%n_v, batch_k%length, wf%n_v, batch_j%length)

               call mem%dealloc(Y_jbka, batch_j%length, wf%n_v, batch_k%length, wf%n_v)
!
            enddo
!
!$omp parallel do private(a, k, b, j)
            do j = 1, batch_j%length
               do b = 1, wf%n_v
                  do k = 1, batch_k%length
                     do a = 1, wf%n_v
!
                        X_akbj(a, k, b, j) = X_akbj(a, k, b, j)/ (omega - eps_v(b) - eps_v(a) &
                                             + eps_o(k + batch_k%first - 1) + eps_o(j+ batch_j%first - 1))
!
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
            call mem%alloc(g_ikbj, wf%n_o, batch_k%length, wf%n_v, batch_j%length)
!
            call wf%eri_t1%get('oovo', g_ikbj, 1, wf%n_o, batch_k%first, batch_k%get_last(), &
                                                   1, wf%n_v, batch_j%first, batch_j%get_last())
!
            call dgemm('N', 'T',                                     &
                        wf%n_v,                                      &
                        wf%n_o,                                      &
                        (batch_k%length)*(batch_j%length)*(wf%n_v),  &
                        -one,                                         &
                        X_akbj,                                      &
                        wf%n_v,                                      &
                        g_ikbj,                                      &
                        wf%n_o,                                      &
                        one,                                         &
                        sigma_ai,                                    &
                        wf%n_v)
!
            call mem%dealloc(g_ikbj, wf%n_o, batch_k%length, wf%n_v, batch_j%length)
            call mem%dealloc(X_akbj, wf%n_v, batch_k%length, wf%n_v, batch_j%length)
!
         enddo
      enddo
!
      call mem%batch_finalize()
!
      call timer%turn_off()
!
   end subroutine effective_jacobian_transpose_cc2_f1_lowmem_cc2
!
!
end submodule jacobian_transpose
