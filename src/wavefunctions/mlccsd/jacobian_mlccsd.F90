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
submodule (mlccsd_class) jacobian_mlccsd
!
!!
!!    Jacobian submodule (MLCCSD)
!!    Written by Sarai D. Folkestad, 2017-2019
!!
!!    Routines for the linear transform of trial
!!    vectors by the Jacobian matrix
!!
!!    rho_i = A * c_i,
!!
!!    where
!!
!!    A_mu,nu = < mu | exp(-T) [H, τ_nu] exp(T) | R >.
!!
!!
!
   implicit none
!
!
contains
!
!
   module subroutine prepare_for_jacobian_mlccsd(wf)
!!
!!    Prepare for Jacobian
!!    Written by Sarai D. Folkestad, 2019
!!
      implicit none 
!
      class(mlccsd), intent(inout) :: wf
!
      real(dp), dimension(:,:,:,:), allocatable :: L_ckdl, g_kcld, g_ckdl, g_klcd
!
      integer :: n_a_o, n_a_v 
!
      type(timings), allocatable :: timer
!
      timer = timings('Prepare for jacobian MLCCSD', pl='n')
      call timer%turn_on()
!
      n_a_o = wf%n_ccsd_o + wf%n_cc2_o
      n_a_v = wf%n_ccsd_v + wf%n_cc2_v
!
      call wf%read_amplitudes()
!
      call wf%construct_x2()
!
      call wf%construct_u_aibj()
!
      if (wf%do_cc2) then
!
!        MO transform AO-fock
!
         call wf%mo_transform(wf%ao_fock, wf%mo_fock)
!
      endif
!
      call wf%save_jacobian_a1_intermediates(n_a_o, n_a_v, 1, 1)
!
      call wf%save_jacobian_c2_intermediates()
      call wf%save_jacobian_d2_intermediate()
!
      call mem%alloc(g_kcld, n_a_o, n_a_v, n_a_o, n_a_v)
      call wf%eri%get_eri_t1('ovov', g_kcld, 1, n_a_o, 1, n_a_v, 1, n_a_o, 1, n_a_v) 
!
      call mem%alloc(L_ckdl, n_a_v, n_a_o, n_a_v, n_a_o)
      call zero_array(L_ckdl, n_a_v**2 * n_a_o**2)
      call add_2143_to_1234(two, g_kcld, L_ckdl, n_a_v, n_a_o, n_a_v, n_a_o)
      call add_2341_to_1234(-one, g_kcld, L_ckdl, n_a_v, n_a_o, n_a_v, n_a_o)
!
      call wf%save_jacobian_e2_intermediate(L_ckdl)
      call wf%save_jacobian_g2_intermediates(L_ckdl)
!
      call mem%dealloc(L_ckdl, n_a_v, n_a_o, n_a_v, n_a_o)
!
      call mem%alloc(g_ckdl, n_a_v, n_a_o, n_a_v, n_a_o)
      call sort_1234_to_2143(g_kcld, g_ckdl, n_a_o, n_a_v, n_a_o, n_a_v)
      call mem%dealloc(g_kcld, n_a_o, n_a_v, n_a_o, n_a_v)
!
      call wf%save_jacobian_h2_intermediates(g_ckdl)
!
!     Reorder g_kcld to g_klcd
!
      call mem%alloc(g_klcd, n_a_o, n_a_o, n_a_v, n_a_v)
      call sort_1234_to_2413(g_ckdl, g_klcd, n_a_v, n_a_o, n_a_v, n_a_o)
      call mem%dealloc(g_ckdl, n_a_v, n_a_o, n_a_v, n_a_o)
!
      call wf%save_jacobian_j2_intermediates(g_klcd)
!
      call mem%dealloc(g_klcd, n_a_o, n_a_o, n_a_v, n_a_v)
!
      call timer%turn_off()
!
   end subroutine prepare_for_jacobian_mlccsd
!
!
   module subroutine jacobian_transformation_mlccsd(wf, c)
!!
!!    Jacobian transformation (MLCCSD)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    Directs the transformation by the CCSD Jacobi matrix,
!!
!!       A_μ,ν = < μ | exp(-T) [H, τ_ν] exp(T) | R >,
!!
!!    where the basis employed for the brackets is biorthonormal.
!!    The transformation is rho = A c, i.e.,
!!
!!       rho_mu = (A c)_mu = sum_ck A_mu,ck c_ck
!!                  + 1/2 sum_ckdl A_mu,ckdl c_ckdl (1 + delta_ck,dl).
!!
!!    On exit, c is overwritten by rho. That is, c_ai = rho_ai,
!!    and c_aibj = rho_aibj.
!!
      implicit none
!
      class(mlccsd), intent(inout) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: c
!
      real(dp), dimension(:,:,:,:), allocatable :: c_aibj, c_abij
!
      real(dp), dimension(:,:), allocatable :: rho_ai
      real(dp), dimension(:,:,:,:), allocatable :: rho_aibj, rho_aibj_ccsd, rho_abij
!
      integer :: i, j, a, b, ai, bj, aibj ! Index
!
      integer :: n_a_o, n_a_v 
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian MLCCSD', pl='n')
      call timer%turn_on()
!
      n_a_o = wf%n_ccsd_o + wf%n_cc2_o
      n_a_v = wf%n_ccsd_v + wf%n_cc2_v
!
!     Allocate and zero the transformed vector (singles part)
!
      call mem%alloc(rho_ai, wf%n_v, wf%n_o)
      call zero_array(rho_ai, wf%n_t1)
!
!     :: CCS contributions to the singles c vector ::
!
      call wf%jacobian_ccs_a1(rho_ai, c(1:wf%n_t1))
      call wf%jacobian_ccs_b1(rho_ai, c(1:wf%n_t1))
!
!     :: CC2/CCSD contributions to the transformed singles vector ::
!
      call wf%jacobian_cc2_a1(rho_ai, c(1:wf%n_t1), n_a_o, n_a_v, 1, 1)
!
      call mem%alloc(c_aibj, n_a_v, n_a_o, n_a_v, n_a_o)
!
!$omp parallel do schedule(static) private(a, i, b, j, ai, bj, aibj) collapse(2)
      do i = 1, n_a_o
         do a = 1, n_a_v
!
            ai = n_a_v*(i - 1) + a
!
            do j = 1, n_a_o
               do b = 1, n_a_v
!
                  bj = n_a_v*(j - 1) + b
!
                  if (ai .ge. bj) then
!
                     aibj = ai*(ai-3)/2 + ai + bj
!
                     c_aibj(a,i,b,j) = c(wf%n_t1 + aibj)
                     c_aibj(b,j,a,i) = c(wf%n_t1 + aibj)
!
                  endif
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
!     Scale the doubles vector by 1 + delta_ai,bj, i.e.
!     redefine to c_ckdl = c_ckdl (1 + delta_ck,dl)
!
!$omp parallel do schedule(static) private(a, i) collapse(2)
      do i = 1, n_a_o
         do a = 1, n_a_v
!
         c_aibj(a,i,a,i) = two*c_aibj(a,i,a,i)
!
         enddo
      enddo
!$omp end parallel do
!
      call wf%jacobian_cc2_b1(rho_ai, c_aibj, n_a_o, n_a_v, &
                             1, 1, n_a_o, n_a_v)
!
!     :: CCSD/CC2 contributions to the transformed doubles vector ::
!
!     Allocate unpacked transformed vector
!
      call mem%alloc(rho_aibj, n_a_v, n_a_o, n_a_v, n_a_o)
      call zero_array(rho_aibj, n_a_v**2 * n_a_o**2)
!
      call wf%jacobian_cc2_a2(rho_aibj, c(1:wf%n_t1), n_a_o, n_a_v,    &
                             1, 1, n_a_o, n_a_v)   
!
      if (wf%do_cc2) call wf%jacobian_cc2_b2(rho_aibj, c_aibj)
!
      call symmetric_sum(rho_aibj, (n_a_v)*(n_a_o))
!
      call mem%alloc(rho_aibj_ccsd, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
      call zero_array(rho_aibj_ccsd, (wf%n_ccsd_v**2)*(wf%n_ccsd_o**2))
!
      call wf%jacobian_ccsd_b2(rho_aibj_ccsd, c(1:wf%n_t1))
      call wf%jacobian_ccsd_c2(rho_aibj_ccsd, c(1:wf%n_t1))
      call wf%jacobian_ccsd_d2(rho_aibj_ccsd, c(1:wf%n_t1))
!
!     Done with singles vector c; overwrite it with
!     transformed vector for exit
!
      call dcopy(wf%n_t1, rho_ai, 1, c, 1)
      call mem%dealloc(rho_ai, wf%n_v, wf%n_o)
!
!     Contributions from doubles vector c
!
      call wf%jacobian_ccsd_e2(rho_aibj_ccsd, c_aibj)
      call wf%jacobian_ccsd_f2(rho_aibj_ccsd, c_aibj)
      call wf%jacobian_ccsd_g2(rho_aibj_ccsd, c_aibj)
      call wf%jacobian_ccsd_h2(rho_aibj_ccsd, c_aibj)
      call wf%jacobian_ccsd_i2(rho_aibj_ccsd, c_aibj)
!
!     Symmetrize rho_aibj before the last two terms are added.
!
!     rho_aibj <- P^{ab}_{ij} rho_aibj = rho_aibj + rho_bjai
!
      call symmetric_sum(rho_aibj_ccsd, (wf%n_ccsd_v)*(wf%n_ccsd_o))
!
!     O^2V^4 term is already coded in omega_mlccsd.F90
!
      call wf%omega_ccsd_b2(rho_aibj_ccsd, c_aibj)
!
      call mem%alloc(c_abij, n_a_v, n_a_v, n_a_o, n_a_o)
      call mem%alloc(rho_abij, wf%n_ccsd_v, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_o)
!
      call zero_array(rho_abij, (wf%n_ccsd_v**2)*(wf%n_ccsd_o**2))
!
      call sort_1234_to_1324(c_aibj, c_abij, n_a_v, n_a_o, n_a_v, n_a_o)
      call mem%dealloc(c_aibj, n_a_v, n_a_o, n_a_v, n_a_o)
!
      call wf%jacobian_ccsd_j2(rho_abij, c_abij)
      call wf%jacobian_ccsd_k2(rho_abij, c_abij)
!
      call mem%dealloc(c_abij, n_a_v, n_a_v, n_a_o, n_a_o)
!
      call add_1324_to_1234(one, rho_abij, rho_aibj_ccsd, &
                              wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
      call mem%dealloc(rho_abij, wf%n_ccsd_v, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_o)
!
!     Order rho_abij back into rho_aibj & divide by
!     the biorthonormal factor 1 + delta_ai,bj
!
!$omp parallel do private(j, b, i, a) collapse(2)
      do j = 1, wf%n_ccsd_o
         do b = 1, wf%n_ccsd_v
            do i = 1, wf%n_ccsd_o
               do a = 1, wf%n_ccsd_v
!
                  rho_aibj(a, i, b, j) = rho_aibj(a, i, b, j) + rho_aibj_ccsd(a,i,b,j)

!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(rho_aibj_ccsd, wf%n_ccsd_v, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_o) 
!
!$omp parallel do schedule(static) private(a,i,b,j) collapse(2)
         do i = 1, n_a_o
            do a = 1, n_a_v
!
            rho_aibj(a,i,a,i) = half*rho_aibj(a,i,a,i)
!
         enddo
      enddo
!$omp end parallel do
!
!     Overwrite the incoming doubles c vector & pack in
!
!$omp parallel do schedule(static) private(a, i, b, j, ai, bj, aibj) collapse(2)
         do j = 1, n_a_o
            do b = 1, n_a_v
!
            bj = n_a_v*(j - 1) + b
!
            do i = 1, n_a_o
               do a = 1, n_a_v
!
                  ai = n_a_v*(i - 1) + a
!
                  if (ai .ge. bj) then
!
                     aibj = ai*(ai-3)/2 + ai + bj
!
                     c((wf%n_o)*(wf%n_v) + aibj) = rho_aibj(a,i,b,j)
!
                  endif
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(rho_aibj, n_a_v, n_a_o, n_a_v, n_a_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_transformation_mlccsd
!
!
   module subroutine jacobian_cc2_b2_mlccsd(wf, rho_aibj, c_aibj)
!!
!!    Jacobian CC2 B2
!!    Written by Sarai Dery Folkestad, 2019
!!
!!    The doubles-doubles part of the CC2 Jacobian in non-canonical basis  
!!
!!       rho_aibj += F_bc c_cjai - F_kj c_aibk
!!
!!    INDEX RESTRICTIONS:
!!
!!    All indices are CC2 indices
!!
      implicit none
!
      class(mlccsd) :: wf
!
      real(dp), dimension(wf%n_cc2_v + wf%n_ccsd_v, wf%n_cc2_o + wf%n_ccsd_o,&
                     wf%n_cc2_v + wf%n_ccsd_v, wf%n_cc2_o + wf%n_ccsd_o), intent(inout) :: c_aibj
      real(dp), dimension(wf%n_cc2_v + wf%n_ccsd_v, wf%n_cc2_o + wf%n_ccsd_o,&
                     wf%n_cc2_v + wf%n_ccsd_v, wf%n_cc2_o + wf%n_ccsd_o), intent(inout) :: rho_aibj   
!
      integer :: n_a_o, n_a_v 
!
      real(dp), dimension(:,:,:,:), allocatable :: rho_aibj_cc2, c_aibj_cc2
!
      integer :: a, i, b, j
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian B2 MLCCSD (CC2)', pl='v')
      call timer%turn_on()
!
      n_a_o = wf%n_ccsd_o + wf%n_cc2_o
      n_a_v = wf%n_ccsd_v + wf%n_cc2_v
!
      call mem%alloc(c_aibj_cc2, n_a_v, n_a_o, n_a_v, n_a_o)
      call mem%alloc(rho_aibj_cc2, n_a_v, n_a_o, n_a_v, n_a_o)
!
      call copy_and_scale(one, c_aibj, c_aibj_cc2, (n_a_o**2)*(n_a_v**2))
!
      call dgemm('N', 'N',                            &
                  n_a_v,                              &
                  (n_a_o**2)*(n_a_v),                 &
                  n_a_v,                              &
                  one,                                &
                  wf%mo_fock(wf%n_o + 1, wf%n_o + 1), & ! F_bc
                  wf%n_mo,                            &
                  c_aibj_cc2,                         & ! c_cjai
                  n_a_v,                              &
                  zero,                               &
                  rho_aibj_cc2,                       & ! We will symmetrize after
                  n_a_v)
!
      call dgemm('N', 'N',             &
                  (n_a_v**2)*(n_a_o),  &
                  n_a_o,               &
                  n_a_o,               &
                  -one,                &
                  c_aibj_cc2,          & ! c_aibk
                  (n_a_v**2)*(n_a_o),  &
                  wf%mo_fock,          & ! F_kj 
                  wf%n_mo,             &
                  one,                 &
                  rho_aibj_cc2,        &
                  (n_a_v**2)*(n_a_o))
!
!     Zero out CCSD part 
!
!$omp parallel do private(j, b, i, a) collapse(2)
      do j = 1, wf%n_ccsd_o
         do b = 1, wf%n_ccsd_v
            do i = 1, wf%n_ccsd_o
               do a = 1, wf%n_ccsd_v
!
                  rho_aibj_cc2(a, i, b, j) = zero
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call daxpy((n_a_o**2)*(n_a_v**2), one, rho_aibj_cc2, 1, rho_aibj, 1)
!
      call mem%dealloc(c_aibj_cc2, n_a_v, n_a_o, n_a_v, n_a_o)
      call mem%dealloc(rho_aibj_cc2, n_a_v, n_a_o, n_a_v, n_a_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_cc2_b2_mlccsd
!
!
   module subroutine jacobian_ccsd_b2_mlccsd(wf, rho_aibj, c_ai)
!!
!!    Jacobian MLCCSD B2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    rho_aibj^B2 = - sum_kc (F_kc x_ij^ac c_bk + F_kc x_ik^ab c_cj)
!!
!!    INDEX RESTRICTIONS:
!!
!!    a, i, b, j are CCSD indices
!!
!!    k is full index in first term and CCSD + CC2 in second term 
!!    c is full index in second term and CCSD + CC2 in first term 
!!
      implicit none
!
      class(mlccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                         :: c_ai
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o) :: rho_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: x_cjai   
      real(dp), dimension(:,:,:,:), allocatable :: x_aibk   
      real(dp), dimension(:,:,:,:), allocatable :: Y_kjai  
!
      real(dp), dimension(:,:), allocatable :: Y_kj         
!
      type(timings), allocatable :: timer
!
      integer :: n_a_o, n_a_v 
!
      integer :: j, i, a, c, ai, cj, aicj, k, bk, aibk, b
!
      timer = timings('Jacobian MLCCSD B2', pl='v')
      call timer%turn_on()
!
      n_a_o = wf%n_ccsd_o + wf%n_cc2_o
      n_a_v = wf%n_ccsd_v + wf%n_cc2_v
!
!     :: Term 1. - sum_kc F_kc x_ij^ac c_bk ::
!
!     Order the amplitudes as x_ca_ij = x_ij^ac
!
      call mem%alloc(x_cjai, n_a_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
!$omp parallel do private (a, i, c, j, ai, cj, aicj) collapse(2)
      do i = 1, wf%n_ccsd_o
         do a = 1, wf%n_ccsd_v
!
            ai = n_a_v*(i - 1) + a
!
            do j = 1, wf%n_ccsd_o
               do c = 1, n_a_v
!
                  cj = n_a_v*(j - 1) + c
!
                  aicj = max(ai, cj)*(max(ai,cj) - 3)/2 + ai + cj
!
                  x_cjai(c, j, a, i) = wf%x2(aicj)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
!     Form the intermediate Y_k_aij = sum_c F_k_c x_c_aij
!
      call mem%alloc(Y_kjai, wf%n_o, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
      call dgemm('N', 'N',                            &
                  wf%n_o,                             &
                  (wf%n_ccsd_v)*(wf%n_ccsd_o**2),     &
                  n_a_v,                              &
                  one,                                &
                  wf%fock_ia,                         & ! F_k,c
                  wf%n_o,                             &
                  x_cjai,                             & ! x_c_jai
                  n_a_v,                              &
                  zero,                               &
                  Y_kjai,                             &
                  wf%n_o)
!
      call mem%dealloc(x_cjai, n_a_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
!     Form rho_b_aij = sum_k c_ai(b,k) X_k_aij(k,aij)
!
      call dgemm('N', 'N',                         &
                  wf%n_ccsd_v,                     &
                  (wf%n_ccsd_v)*(wf%n_ccsd_o)**2,  &
                  wf%n_o,                          &
                  -one,                            &
                  c_ai,                            & ! c_bk
                  wf%n_v,                          &
                  Y_kjai,                          & 
                  wf%n_o,                          &
                  zero,                            &
                  rho_aibj,                        & ! Will symmetrize later
                  wf%n_ccsd_v)  
!
      call mem%dealloc(Y_kjai, wf%n_o, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_o)
!
!     :: Term 2. - sum_kc F_kc x_ik^ab c_cj ::
!
!     Form Y_kj = sum_c F_kc c_cj 
!
      call mem%alloc(Y_kj, n_a_o, wf%n_ccsd_o)
!
      call dgemm('N','N',                          &
                  n_a_o,                           &
                  wf%n_ccsd_o,                     &
                  wf%n_v,                          &
                  one,                             &
                  wf%fock_ia,                      & ! F_k,c
                  wf%n_o,                          &
                  c_ai,                            & ! c_c,j
                  wf%n_v,                          &
                  zero,                            &
                  Y_kj,                            &
                  n_a_o)
!
!     Form rho_aib_j = - sum_k x_aib_k Y_k_j
!
      call mem%alloc(x_aibk, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, n_a_o)
!
!$omp parallel do private(k, b, i, a, bk, ai, aibk) collapse(2) 
      do k = 1, n_a_o
         do b = 1, wf%n_ccsd_v
!
            bk = n_a_v*(k - 1) + b
!
            do i = 1, wf%n_ccsd_o
               do a = 1, wf%n_ccsd_v
!
                  ai = n_a_v*(i-1) + a
                  aibk = max(ai,bk)*(max(ai,bk)-3)/2 + ai + bk
!
                  x_aibk(a, i, b, k) = wf%x2(aibk)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
       call dgemm('N','N',                          &
                   (wf%n_ccsd_o)*(wf%n_ccsd_v)**2,  &
                   wf%n_ccsd_o,                     &
                   n_a_o,                           &
                   -one,                            &
                   x_aibk,                          & 
                   (wf%n_ccsd_o)*(wf%n_ccsd_v)**2,  &
                   Y_kj,                            &
                   n_a_o,                           &
                   one,                             &
                   rho_aibj,                        & ! rho_aib,j
                   (wf%n_ccsd_o)*(wf%n_ccsd_v)**2)
!
      call mem%dealloc(Y_kj, n_a_o, wf%n_ccsd_o)
      call mem%dealloc(x_aibk, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, n_a_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_ccsd_b2_mlccsd
!
!
   module subroutine jacobian_ccsd_c2_1_mlccsd(wf, rho_aibj, c_ai)
!!
!!    Jacobian MLCCSD C2-1
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    rho_aibj^C2 += sum_kcl g_ljkc x_ki^ac c_bl
!!
!!    INDEX RESTRICTIONS:
!!
!!     a, i, b, j are CCSD indices
!!
!!     c, k are CCSD + CC2, l is unrestricted
!!
      implicit none
!
      class(mlccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o), &
                                           intent(inout) :: rho_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: X_ljai
!
      integer :: n_a_o, n_a_v 
!
      n_a_o = wf%n_ccsd_o + wf%n_cc2_o
      n_a_v = wf%n_ccsd_v + wf%n_cc2_v
!
!     :: Term 1. sum_kcl g_ljkc x_ki^ac c_bl ::
!
      call mem%alloc(X_ljai, wf%n_o, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
!     Form X_ljai = sum_ck g_ljkc x_ki^ac
!
      call wf%jacobian_c2_intermediate_oovo_1%open_('read', 'rewind')
      call wf%jacobian_c2_intermediate_oovo_1%read_(X_ljai, (wf%n_o)*(wf%n_ccsd_v)*(wf%n_ccsd_o**2))
      call wf%jacobian_c2_intermediate_oovo_1%close_('keep')
!
!     Calculate rho_b_jai = sum_l c_bl X_ljai
!
      call dgemm('N', 'N',                         &
                  wf%n_ccsd_v,                     &
                  (wf%n_ccsd_v)*(wf%n_ccsd_o)**2,  &
                  wf%n_o,                          &
                  one,                             &
                  c_ai,                            & ! c_bl
                  wf%n_v,                          &
                  X_ljai,                          & 
                  wf%n_o,                          &
                  one,                             &
                  rho_aibj,                        & ! We will symmetrize after
                  wf%n_ccsd_v)
!
      call mem%dealloc(X_ljai, wf%n_o, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
   end subroutine jacobian_ccsd_c2_1_mlccsd
!
!
   module subroutine jacobian_ccsd_c2_2_mlccsd(wf, rho_aibj, c_ai)
!!
!!    Jacobian MLCCSD C2-2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    rho_aibj^C2 += sum_kcl g_ljkc x_li^bc c_ak 
!!
!!    INDEX RESTRICTIONS:
!!
!!    a, i, b, j are CCSD indices
!!
!!    l and c are CCSD + CC2, k is unrestricted
!!
      implicit none
!
      class(mlccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
!
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, &
                              wf%n_ccsd_v, wf%n_ccsd_o), intent(inout) :: rho_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: X_kjbi
      real(dp), dimension(:,:,:,:), allocatable :: rho_ajbi
!
      integer :: n_a_o, n_a_v 
!
      n_a_o = wf%n_ccsd_o + wf%n_cc2_o
      n_a_v = wf%n_ccsd_v + wf%n_cc2_v
!
!     :: Term 1: sum_kcl g_ljkc x_li^bc c_ak ::
!
!     X_kjbi = g_ljkc x_li^bc
!
      call mem%alloc(X_kjbi, wf%n_o, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
      wf%jacobian_c2_intermediate_oovo_2 = sequential_file('jacobian_c2_intermediate_oovo_2_ccsd')
      call wf%jacobian_c2_intermediate_oovo_2%open_('read', 'rewind')
      call wf%jacobian_c2_intermediate_oovo_2%read_(X_kjbi, (wf%n_o)*(wf%n_ccsd_v)*(wf%n_ccsd_o**2))
      call wf%jacobian_c2_intermediate_oovo_2%close_('keep')
!
!     Calculate rho_a_jbi = sum_k c_ak X_kjbi
!
      call mem%alloc(rho_ajbi, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
      call dgemm('N', 'N',                         &
                  wf%n_ccsd_v,                     &
                  (wf%n_ccsd_v)*(wf%n_ccsd_o)**2,  &
                  wf%n_o,                          &
                  one,                             &
                  c_ai,                            & ! c_ak
                  wf%n_v,                          &
                  X_kjbi,                          &
                  wf%n_o,                          &
                  zero,                            &
                  rho_ajbi,                        & 
                  wf%n_ccsd_v)
!
      call mem%dealloc(X_kjbi, wf%n_o, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
      call add_1432_to_1234(one, rho_ajbi, rho_aibj, &
                           wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
      call mem%dealloc(rho_ajbi, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
     end subroutine jacobian_ccsd_c2_2_mlccsd
!
!
   module subroutine jacobian_ccsd_c2_3_mlccsd(wf, rho_aibj, c_ai, g_ljkc)
!!
!!    Jacobian MLCCSD C2-3
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    rho_aibj^C2 += sum_kcl g_ljkc x_lk^ba c_ci
!!
!!    INDEX RESTRICTIONS:
!!
!!    a, i, b, j are CCSD indices
!!
!!    k and l are CCSD + CC2, c is unrestricted
!!
      implicit none
!
      class(mlccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
!
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, &
                  wf%n_ccsd_v, wf%n_ccsd_o), intent(inout) :: rho_aibj
!
      real(dp), dimension(wf%n_ccsd_o + wf%n_cc2_o, wf%n_ccsd_o, &
                  wf%n_ccsd_o + wf%n_cc2_o, wf%n_v), intent(in) :: g_ljkc
!
      integer :: n_a_o, n_a_v 
!
      real(dp), dimension(:,:,:,:), allocatable :: Y_ljki
      real(dp), dimension(:,:,:,:), allocatable :: Y_klij
      real(dp), dimension(:,:,:,:), allocatable :: x_bakl
      real(dp), dimension(:,:,:,:), allocatable :: rho_baij
!
      integer :: b, a, k, l, bl, ak, blak
!
      n_a_o = wf%n_ccsd_o + wf%n_cc2_o
      n_a_v = wf%n_ccsd_v + wf%n_cc2_v
!
!     :: Term 1. sum_kcl g_ljkc x_lk^ba c_ci ::
!
!     Form the intermediate Y_ljki = sum_c g_ljkc c_ci
!
      call mem%alloc(Y_ljki, n_a_o, wf%n_ccsd_o, n_a_o, wf%n_ccsd_o)
!
      call dgemm('N','N',                    &
                  (n_a_o**2)*wf%n_ccsd_o,    &
                  wf%n_ccsd_o,               &
                  wf%n_v,                    &
                  one,                       &
                  g_ljkc,                    & 
                  (n_a_o**2)*wf%n_ccsd_o,    &
                  c_ai,                      & ! c_c_i
                  wf%n_v,                    &
                  zero,                      &
                  Y_ljki,                    & 
                  (n_a_o**2)*wf%n_ccsd_o)
!
!     Reorder to Y_ljki as Y_kjli
!
      call mem%alloc(Y_klij, n_a_o, n_a_o, wf%n_ccsd_o, wf%n_ccsd_o)
!
      call sort_1234_to_3142(Y_ljki, Y_klij, n_a_o, wf%n_ccsd_o, n_a_o, wf%n_ccsd_o)
!
      call mem%dealloc(Y_ljki, n_a_o, wf%n_ccsd_o, n_a_o, wf%n_ccsd_o)
!
!     Order amplitudes as x_ba_kl = x_lk^ba
!
      call mem%alloc(x_bakl, wf%n_ccsd_v, wf%n_ccsd_v, n_a_o, n_a_o)
!
!$omp parallel do private (l, k, a, b, ak, bl, blak) collapse(2)
      do l = 1, n_a_o
         do k = 1, n_a_o
            do a = 1, wf%n_ccsd_v
!
               ak = n_a_v*(k - 1) + a
!
               do b = 1, wf%n_ccsd_v
!
                  bl = n_a_v*(l - 1) + b
                  blak = max(bl,ak)*(max(bl,ak) - 3)/2 + bl + ak
!
                  x_bakl(b, a, k, l) = wf%x2(blak)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
!     Calculate rho_baij = sum_kcl g_ljkc x_lk^ba c_ci
!                         = sum_kl x_bakl Y_klij
!
      call mem%alloc(rho_baij, wf%n_ccsd_v, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_o)
!
      call dgemm('N','N',           &
                  (wf%n_ccsd_v)**2, &
                  (wf%n_ccsd_o)**2, &
                  (n_a_o)**2,       &
                  one,              &
                  x_bakl,           & 
                  (wf%n_ccsd_v)**2, &
                  Y_klij,           & 
                  (n_a_o)**2,       &
                  zero,             &
                  rho_baij,         & ! rho_ba_ij
                  (wf%n_ccsd_v)**2)
!
      call mem%dealloc(Y_klij, n_a_o, n_a_o, wf%n_ccsd_o, wf%n_ccsd_o)
      call mem%dealloc(x_bakl, wf%n_ccsd_v, wf%n_ccsd_v, n_a_o, n_a_o)
!
      call add_3124_to_1234(one, rho_baij, rho_aibj, &
                           wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
      call mem%dealloc(rho_baij, wf%n_ccsd_v, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_o)
!
     end subroutine jacobian_ccsd_c2_3_mlccsd
!
!
   module subroutine jacobian_ccsd_c2_4_mlccsd(wf, rho_aibj, c_ai, L_ljck)
!!
!!    Jacobian MLCCSD C2-4
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    rho_aibj^C2 += - sum_kcl L_ljkc x_il^ab c_ck  
!!
!!    INDEX RESTRICTIONS:
!!
!!    a, i, b, j are CCSD indices
!!
!!    l is CCSD + CC2, c and k are unrestricted
!!
      implicit none
!
      class(mlccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
!
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, &
                  wf%n_ccsd_v, wf%n_ccsd_o), intent(inout) :: rho_aibj
!
      real(dp), dimension(wf%n_ccsd_o + wf%n_cc2_o, &
                  wf%n_ccsd_o, wf%n_v, wf%n_o), intent(in) :: L_ljck
!
      integer :: n_a_o, n_a_v 
!
      real(dp), dimension(:,:), allocatable :: Y_lj
      real(dp), dimension(:,:,:,:), allocatable :: x_aibl
!
      integer :: a, i, b, l, ai, bl, aibl
!
      n_a_o = wf%n_ccsd_o + wf%n_cc2_o
      n_a_v = wf%n_ccsd_v + wf%n_cc2_v
!
!    :: Term 1. - sum_kcl L_ljkc x_il^ab c_ck ::
!
!     Calculate the intermediate Y_lj = sum_ck L_ljck c_ck
!
      call mem%alloc(Y_lj, n_a_o, wf%n_ccsd_o)
!
      call dgemm('N', 'N',                &
                  (n_a_o)*(wf%n_ccsd_o),  &
                  1,                      &
                  (wf%n_o)*(wf%n_v),      &
                  one,                    &
                  L_ljck,                 & 
                  (n_a_o)*(wf%n_ccsd_o),  &
                  c_ai,                   & ! c_ck
                  (wf%n_o)*(wf%n_v),      &
                  zero,                   &
                  Y_lj,                   & 
                  (n_a_o)*(wf%n_ccsd_o))
!
!     Order the amplitudes as x_ai_bl = x_il^ab
!
      call mem%alloc(x_aibl, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, n_a_o)
!
!$omp parallel do private(l, b, i, a, bl, ai, aibl) collapse(2)
      do l = 1, n_a_o
         do b = 1, wf%n_ccsd_v
!
            bl = n_a_v*(l - 1) + b
!
            do i = 1, wf%n_ccsd_o
               do a = 1, wf%n_ccsd_v
!
                  ai = n_a_v*(i - 1) + a
                  aibl = max(ai,bl)*(max(ai,bl) - 3)/2 + ai + bl
!
                  x_aibl(a, i, b, l) = wf%x2(aibl)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
!     Form rho_aibj =+ - sum_l x_il^ab Y_lj = - sum_l x_aib_l Y_lj
!
      call dgemm('N','N',                          &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v)**2,  &
                  wf%n_ccsd_o,                     &
                  n_a_o,                           &
                  -one,                            &
                  x_aibl,                          & 
                  (wf%n_ccsd_o)*(wf%n_ccsd_v)**2,  &
                  Y_lj,                            & 
                  n_a_o,                           &
                  one,                             &
                  rho_aibj,                        & 
                  (wf%n_ccsd_o)*(wf%n_ccsd_v)**2)
!
      call mem%dealloc(Y_lj, n_a_o, wf%n_ccsd_o)
      call mem%dealloc(x_aibl, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, n_a_o)
!
     end subroutine jacobian_ccsd_c2_4_mlccsd
!
!
   module subroutine jacobian_ccsd_c2_mlccsd(wf, rho_aibj, c_ai)
!!
!!    Jacobian CCSD C2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
      implicit none
!
      class(mlccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o), &
                                           intent(inout) :: rho_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ljkc
      real(dp), dimension(:,:,:,:), allocatable :: g_ljkc_3
!
      real(dp), dimension(:,:,:,:), allocatable :: L_ljck
!
      integer :: n_a_o, n_a_v, l, j, c, k
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian MLCCSD C2', pl='v')
      call timer%turn_on()
!
      n_a_o = wf%n_ccsd_o + wf%n_cc2_o
      n_a_v = wf%n_ccsd_v + wf%n_cc2_v               
!
      call wf%jacobian_ccsd_c2_1(rho_aibj, c_ai) 
      call wf%jacobian_ccsd_c2_2(rho_aibj, c_ai)
!
      call mem%alloc(g_ljkc, wf%n_o, wf%n_ccsd_o, wf%n_o, wf%n_v)
!
      call wf%eri%get_eri_t1('ooov', g_ljkc, 1, wf%n_o, 1, wf%n_ccsd_o, 1, wf%n_o, 1, wf%n_v)
!
      call mem%alloc(g_ljkc_3, n_a_o, wf%n_ccsd_o, n_a_o, wf%n_v)
!
!$omp parallel do private(c, k, l, j)
      do c = 1, wf%n_v
         do k = 1, n_a_o
            do j = 1, wf%n_ccsd_o
               do l = 1, n_a_o
!
                  g_ljkc_3(l, j, k, c) = g_ljkc(l, j, k, c)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call wf%jacobian_ccsd_c2_3(rho_aibj, c_ai, g_ljkc_3)
!
      call mem%dealloc(g_ljkc_3, n_a_o, wf%n_ccsd_o, n_a_o, wf%n_v)
!
      call mem%alloc(L_ljck, n_a_o, wf%n_ccsd_o, wf%n_v, wf%n_o)
!
!$omp parallel do private(c, k, l, j)
      do k = 1, wf%n_o
         do c = 1, wf%n_v
            do j = 1, wf%n_ccsd_o
               do l = 1, n_a_o
!
                  L_ljck(l, j, c, k) = two*g_ljkc(l, j, k, c) - g_ljkc(k, j, l, c) 
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!    
      call mem%dealloc(g_ljkc, wf%n_o, wf%n_ccsd_o, wf%n_o, wf%n_v)
!
      call wf%jacobian_ccsd_c2_4(rho_aibj, c_ai, L_ljck)
!
      call mem%dealloc(L_ljck, n_a_o, wf%n_ccsd_o, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_ccsd_c2_mlccsd
!
!
  module subroutine jacobian_ccsd_d2_1_mlccsd(wf, rho_aibj, c_ai)
!!
!!    Jacobian CCSD D2-1
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    rho_aibj^D2 += - sum_kcd g_kcbd x_ij^cd c_ak
!!
!     INDEX RESTRICTIONS 
!
!     k unrestricted, c, d CC2 + CCSD
!
      implicit none
!
      class(mlccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o), &
                  intent(inout) :: rho_aibj
!
      integer :: n_a_o, n_a_v 
!
      real(dp), dimension(:,:,:,:), allocatable :: X_kibj
!
      n_a_o = wf%n_ccsd_o + wf%n_cc2_o
      n_a_v = wf%n_ccsd_v + wf%n_cc2_v
!
!     :: Term 1. - sum_kcd g_kcbd x_ij^cd c_ak ::
!
!     X_kibj = g_kcbd x_ij^cd
!     
      call mem%alloc(X_kibj, wf%n_o, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
      wf%jacobian_d2_intermediate = sequential_file('jacobian_d2_intermediate_ccsd')
      call wf%jacobian_d2_intermediate%open_('read', 'rewind')
!
      call wf%jacobian_d2_intermediate%read_(X_kibj, wf%n_o*(wf%n_ccsd_o**2)*wf%n_ccsd_v)
!
      call wf%jacobian_d2_intermediate%close_('keep')
!
!     Form rho_aibj = - sum_k c_ak X_kibj 
!
      call dgemm('N', 'N',                            &
                  wf%n_ccsd_v,                        &
                  (wf%n_ccsd_v)*(wf%n_ccsd_o)**2,     &
                  wf%n_o,                             &
                  -one,                               &
                  c_ai,                               & ! c_a_k
                  wf%n_v,                             &
                  X_kibj,                             & ! X_k_ibj
                  wf%n_o,                             &
                  one,                                &
                  rho_aibj,                           & ! rho_a_ibj
                  wf%n_ccsd_v)
!
      call mem%dealloc(X_kibj, wf%n_o, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
   end subroutine jacobian_ccsd_d2_1_mlccsd
!
!
  module subroutine jacobian_ccsd_d2_2_mlccsd(wf, rho_aibj, c_ai, batch_b, g_dkbc)
!!
!!    Jacobian CCSD D2-2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    rho_aibj^D2 += - sum_kcd g_kcbd x_kj^ad c_ci
!!
!     INDEX RESTRICTIONS 
!
!     c unrestricted, k, d CC2 + CCSD
!
      implicit none
!
      class(mlccsd) :: wf
!
      type(batching_index), intent(in) :: batch_b
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, batch_b%length, wf%n_ccsd_o), &
                     intent(inout) :: rho_aibj
!
      real(dp), dimension(wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o, &
                           batch_b%length, wf%n_v), intent(in) :: g_dkbc
!
      integer :: n_a_o, n_a_v 
!  
      real(dp), dimension(:,:,:,:), allocatable :: Y_dkbi
      real(dp), dimension(:,:,:,:), allocatable :: x_ajdk
      real(dp), dimension(:,:,:,:), allocatable :: rho_ajbi
!
      integer :: a, j, d, k, ak, dj, akdj
!
      n_a_o = wf%n_ccsd_o + wf%n_cc2_o
      n_a_v = wf%n_ccsd_v + wf%n_cc2_v
!
!     :: Term 1. - sum_kcd g_kcbd x_kj^ad c_ci ::
!
!     Form the intermediate Y_dkbi = sum_c g_kcbd c_ci
!
      call mem%alloc(Y_dkbi, n_a_v, n_a_o, batch_b%length, wf%n_ccsd_o)
!
      call dgemm('N','N',                                   &
                  n_a_v*(n_a_o)*(batch_b%length),           &
                  (wf%n_ccsd_o),                            &
                  wf%n_v,                                   &
                  one,                                      &
                  g_dkbc,                                   &
                  n_a_v*(n_a_o)*(batch_b%length),           &
                  c_ai,                                     & ! c_ci
                  wf%n_v,                                   &
                  zero,                                     &
                  Y_dkbi,                                   &
                  n_a_v*(n_a_o)*(batch_b%length))
!
!     Order x_akdj as x_ajdk
!
      call mem%alloc(x_ajdk, wf%n_ccsd_v, wf%n_ccsd_o, n_a_v, n_a_o)
!
!$omp parallel do private(k, d, j, a, ak, dj, akdj) collapse(2)
      do k = 1, n_a_o 
         do d = 1, n_a_v
            do j = 1, wf%n_ccsd_o
!
               dj = n_a_v*(j-1) + d
!
               do a = 1, wf%n_ccsd_v
!
                  ak = n_a_v*(k - 1) + a
                  akdj = max(ak,dj)*(max(ak,dj)-3)/2 + ak + dj
!
                  x_ajdk(a, j, d, k) = wf%x2(akdj)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do 
!
!    Calculate rho_ib_aj = - sum_dk Y_dkbi x_ajdk
!
     call mem%alloc(rho_ajbi, wf%n_ccsd_v, wf%n_ccsd_o, batch_b%length, wf%n_ccsd_o)
!
     call dgemm('N','N',                           &
                 (wf%n_ccsd_o)*(wf%n_ccsd_v),      &
                 (wf%n_ccsd_o)*(batch_b%length),   &
                 (n_a_o)*(n_a_v),                  &
                 -one,                             &
                 x_ajdk,                           &
                 (wf%n_ccsd_o)*(wf%n_ccsd_v),      &
                 Y_dkbi,                           & 
                 (n_a_o)*(n_a_v),                  &
                 zero,                             &
                 rho_ajbi,                         & 
                 (wf%n_ccsd_o)*(wf%n_ccsd_v))
!
      call mem%dealloc(Y_dkbi, n_a_v, n_a_o, batch_b%length, wf%n_ccsd_o)
      call mem%dealloc(x_ajdk, wf%n_ccsd_v, wf%n_ccsd_o, n_a_v, n_a_o)
!
      call add_1432_to_1234(one, rho_ajbi, rho_aibj, &
            wf%n_ccsd_v, wf%n_ccsd_o, batch_b%length, wf%n_ccsd_o)
!
      call mem%dealloc(rho_ajbi, wf%n_ccsd_v, wf%n_ccsd_o, batch_b%length, wf%n_ccsd_o)
!
   end subroutine jacobian_ccsd_d2_2_mlccsd
!
!
   module subroutine jacobian_ccsd_d2_3_mlccsd(wf, rho_aibj, c_ai, batch_b, g_kcbd)
!!
!!    Jacobian CCSD D2-3
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    rho_aibj^D2 += - sum_kcd g_kcbd x_ik^ca c_dj
!!
!     INDEX RESTRICTIONS 
!
!     d unrestricted, k, c CC2 + CCSD
!
      implicit none
!
      class(mlccsd) :: wf
!
      type(batching_index), intent(in) :: batch_b
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
!
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, &
                           batch_b%length, wf%n_ccsd_o), intent(inout) :: rho_aibj
!
      real(dp), dimension(wf%n_ccsd_o + wf%n_cc2_o, &
                  wf%n_ccsd_v + wf%n_cc2_v, batch_b%length, wf%n_v), intent(in) :: g_kcbd
!
      integer :: n_a_o, n_a_v 
!
      real(dp), dimension(:,:,:,:), allocatable :: Y_kcbj
      real(dp), dimension(:,:,:,:), allocatable :: x_aikc
!
      integer :: a, i, c, k, ci, ak, ciak
!
      n_a_o = wf%n_ccsd_o + wf%n_cc2_o
      n_a_v = wf%n_ccsd_v + wf%n_cc2_v
!
!     :: Term 1. - sum_kcd g_kcbd x_ik^ca c_dj ::
!
!     Form the intermediate Y_ckbj = sum_d g_kcbd c_dj
!
      call mem%alloc(Y_kcbj, n_a_o, n_a_v, batch_b%length, wf%n_ccsd_o)
!
      call dgemm('N','N',                             &
                  (n_a_v)*(n_a_o)*(batch_b%length),   &
                  wf%n_ccsd_o,                        &
                  wf%n_v,                             &
                  one,                                &
                  g_kcbd,                             & 
                  (n_a_v)*(n_a_o)*(batch_b%length),   &
                  c_ai,                               & ! c_dj
                  wf%n_v,                             &
                  zero,                               &
                  Y_kcbj,                             &
                  (n_a_v)*(n_a_o)*(batch_b%length))
!
!     x_ciak ordered as x_aikc
!
      call mem%alloc(x_aikc, wf%n_ccsd_v, wf%n_ccsd_o, n_a_o, n_a_v)
!
!$omp parallel do private(c, k, i, a, ci, ak, ciak) collapse(2)
      do c = 1, n_a_v
         do k = 1, n_a_o
            do i = 1, wf%n_ccsd_o
!
               ci = n_a_v*(i - 1) + c
!
               do a = 1, wf%n_ccsd_v
!
                  ak = n_a_v*(k - 1) + a
                  ciak = max(ci, ak)*(max(ci, ak) - 3)/2 + ci + ak
!
                  x_aikc(a, i, k, c) = wf%x2(ciak)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call dgemm('N','N',                             &
                     (wf%n_ccsd_v)*(wf%n_ccsd_o),     &
                     (wf%n_ccsd_o)*(batch_b%length),  &
                     (n_a_v)*(n_a_o),                 &
                     -one,                            &
                     x_aikc,                          & 
                     (wf%n_ccsd_v)*(wf%n_ccsd_o),     &
                     Y_kcbj,                          & 
                     (n_a_v)*(n_a_o),                 &
                     one,                             &
                     rho_aibj,                        & 
                     (wf%n_ccsd_v)*(wf%n_ccsd_o))
!
      call mem%dealloc(Y_kcbj, n_a_o, n_a_v, batch_b%length, wf%n_ccsd_o)
      call mem%dealloc(x_aikc, wf%n_ccsd_v, wf%n_ccsd_o, n_a_o, n_a_v)
!
   end subroutine jacobian_ccsd_d2_3_mlccsd
!
!
  module subroutine jacobian_ccsd_d2_4_mlccsd(wf, rho_aibj, c_ai, batch_b, L_ckbd)
!!
!!    Jacobian CCSD D2-4
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    rho_aibj^D2 += sum_kcd L_kcbd x_ik^ac c_dj
!!
!     INDEX RESTRICTIONS 
!
!    d unrestricted, k, c CC2 + CCSD
!
      implicit none
!
      class(mlccsd) :: wf
!
      type(batching_index), intent(in) :: batch_b
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
!
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, &
                  batch_b%length, wf%n_ccsd_o), intent(inout) :: rho_aibj
!
      real(dp), dimension(wf%n_ccsd_v + wf%n_cc2_v, &
                  wf%n_ccsd_o + wf%n_cc2_o, batch_b%length, wf%n_v), intent(in) :: L_ckbd
!
      integer :: n_a_o, n_a_v 
!
      real(dp), dimension(:,:,:,:), allocatable :: Y_ckbj
      real(dp), dimension(:,:,:,:), allocatable :: x_aick
!
      integer :: a, i, c, k, ai, ck, aick
!
      n_a_o = wf%n_ccsd_o + wf%n_cc2_o
      n_a_v = wf%n_ccsd_v + wf%n_cc2_v
!
!     :: Term 1.  sum_kcd L_kcbd x_ik^ac c_dj ::
!
!     Y_ckbj = sum_d L_kcbd c_dj
!
      call mem%alloc(Y_ckbj, n_a_v, n_a_o, batch_b%length, wf%n_ccsd_o)
!
      call dgemm('N','N',                             &
                  (n_a_v)*(n_a_o)*(batch_b%length),   &
                  wf%n_ccsd_o,                        &
                  wf%n_v,                             &
                  one,                                &
                  L_ckbd,                             & 
                  (n_a_v)*(n_a_o)*(batch_b%length),   &
                  c_ai,                               & ! c_dj
                  wf%n_v,                             &
                  zero,                               &
                  Y_ckbj,                             & 
                  (n_a_v)*(n_a_o)*(batch_b%length))
!
      call mem%alloc(x_aick, wf%n_ccsd_v, wf%n_ccsd_o, n_a_v, n_a_o)
!
!$omp parallel do private (k, c, i, a, ck, ai, aick) collapse(2)
      do k = 1, n_a_o
         do c = 1, n_a_v
!
            ck = n_a_v*(k - 1) + c
!
            do i = 1, wf%n_ccsd_o
               do a = 1, wf%n_ccsd_v
!
                  ai = n_a_v*(i - 1) + a
                  aick = max(ai,ck)*(max(ai,ck) - 3)/2 + ai + ck
!
                  x_aick(a, i, c, k) = wf%x2(aick)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
!     Form rho_aibj =  sum_ck x_aick Y_ckbj
!

      call dgemm('N','N',                          &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),     &
                  (wf%n_ccsd_o)*(batch_b%length),  &
                  (n_a_o)*(n_a_v),                 &
                  one,                             &
                  x_aick,                          & 
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),     &
                  Y_ckbj,                          & 
                  (n_a_o)*(n_a_v),                 &
                  one,                             &
                  rho_aibj,                        & ! rho_ai_bj
                  (wf%n_ccsd_o)*(wf%n_ccsd_v))
!
      call mem%dealloc(Y_ckbj, n_a_v, n_a_o, batch_b%length, wf%n_ccsd_o)
      call mem%dealloc(x_aick, wf%n_ccsd_v, wf%n_ccsd_o, n_a_v, n_a_o)
!
   end subroutine jacobian_ccsd_d2_4_mlccsd
!
!
  module subroutine jacobian_ccsd_d2_5_mlccsd(wf, rho_aibj, c_ai, batch_b, L_ckbd)
!!
!!    Jacobian CCSD D2-4
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    rho_aibj^D2 +=  sum_kcd L_kcbd x_ij^ad c_ck
!!
!     INDEX RESTRICTIONS 
!
!     c, k unrestricted, d CC2 + CCSD
!
      implicit none
!
      class(mlccsd) :: wf
!
      type(batching_index), intent(in) :: batch_b
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
!
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, &
                  batch_b%length, wf%n_ccsd_o), intent(inout) :: rho_aibj
!
      real(dp), dimension(wf%n_v, wf%n_o, batch_b%length, &
                  wf%n_ccsd_v + wf%n_cc2_v), intent(in) :: L_ckbd
!
      integer :: n_a_o, n_a_v 
!
      real(dp), dimension(:,:), allocatable :: X_bd
!
      real(dp), dimension(:,:,:,:), allocatable :: x_aijd
      real(dp), dimension(:,:,:,:), allocatable :: rho_baij
!
      integer :: d, j, i, a, dj, ai, aidj
!
      n_a_o = wf%n_ccsd_o + wf%n_cc2_o
      n_a_v = wf%n_ccsd_v + wf%n_cc2_v
!
!     :: Term 1.  sum_kcd L_kcbd x_ij^ad c_ck ::
!
     call mem%alloc(X_bd, (batch_b%length), (n_a_v))
!
      call dgemm('N','N',                    &
                  1,                         &
                  (batch_b%length)*(n_a_v),  &
                  (wf%n_o)*(wf%n_v),         &
                  one,                       &
                  c_ai,                      & ! c_ck
                  1,                         &
                  L_ckbd,                    & ! L_ck_bd
                  (wf%n_o)*(wf%n_v),         &
                  zero,                      &
                  X_bd,                      & ! X_bd
                  1)
!
      call mem%alloc(x_aijd, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_o, n_a_v)
!
!$omp parallel do private(d, j, i, a, dj, ai, aidj) collapse(2)
      do d = 1, n_a_v
         do j = 1, wf%n_ccsd_o
!
            dj = n_a_v*(j - 1) + d
!
            do i = 1, wf%n_ccsd_o
               do a = 1, wf%n_ccsd_v
!
                  ai = n_a_v*(i - 1) + a
                  aidj = max(ai,dj)*(max(ai,dj) - 3)/2 + ai + dj
!
                  x_aijd(a, i, j, d) = wf%x2(aidj)
!
               enddo
            enddo
         enddo
      enddo 
!$omp end parallel do
!
!     Form rho_b_aij =  =  sum_d X_bd t_d_aij
!
      call mem%alloc(rho_baij, batch_b%length, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_o)
!
      call dgemm('N','T',                          &
                  batch_b%length,                  &
                  (wf%n_ccsd_v)*(wf%n_ccsd_o)**2,  &
                  n_a_v,                           &
                  one,                             &
                  X_bd,                            & 
                  batch_b%length,                  &
                  x_aijd,                          &
                  (wf%n_ccsd_v)*(wf%n_ccsd_o)**2,  &
                  zero,                            &
                  rho_baij,                        & ! rho_b_aij
                  batch_b%length)
!
         call mem%dealloc(X_bd, (batch_b%length), (n_a_v))
         call mem%dealloc(x_aijd, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_o, n_a_v)
!
         call add_3124_to_1234(one, rho_baij, rho_aibj,&
                  wf%n_ccsd_v, wf%n_ccsd_o, batch_b%length, wf%n_ccsd_o)
!
         call mem%dealloc(rho_baij, batch_b%length, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_o)
!
   end subroutine jacobian_ccsd_d2_5_mlccsd
!
!
   module subroutine jacobian_ccsd_d2_mlccsd(wf, rho_aibj, c_ai)
!!
!!    Jacobian CCSD D2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
      implicit none
!
      class(mlccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_ai
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o), &
                        intent(inout) :: rho_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: rho_aibj_batch
!
!
      real(dp), dimension(:,:,:,:), allocatable :: g_kcbd
      real(dp), dimension(:,:,:,:), allocatable :: g_dkbc
      real(dp), dimension(:,:,:,:), allocatable :: g_kcbd_3
!
      real(dp), dimension(:,:,:,:), allocatable :: L_ckbd
      real(dp), dimension(:,:,:,:), allocatable :: L_ckbd_4
      real(dp), dimension(:,:,:,:), allocatable :: L_ckbd_5
      integer :: n_a_o, n_a_v 
!
      integer :: rec1, rec0
      integer :: current_b_batch
!
      type(batching_index) :: batch_b

      integer :: b, i, j, a, d, k, c
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian MLCCSD D2', pl='v')
      call timer%turn_on()
!
      n_a_o = wf%n_ccsd_o + wf%n_cc2_o
      n_a_v = wf%n_ccsd_v + wf%n_cc2_v
!
      call wf%jacobian_ccsd_d2_1(rho_aibj, c_ai)
!
!     Initialize batching variable
!
      rec0 = max((wf%eri%n_J)*(wf%n_o)*(wf%n_v), &
                 (wf%n_ccsd_o**2)*(n_a_v), &
                 (n_a_o)*(n_a_v)*(wf%n_ccsd_o)*(wf%n_ccsd_v), &
                 (wf%n_ccsd_o**2)*(wf%n_ccsd_v)*(n_a_v))

      rec1 = (wf%n_o)*(wf%n_v**2) + max((wf%n_o)*(wf%n_v**2), &
             2*(wf%n_o)*(wf%n_ccsd_o**2) + (wf%n_ccsd_o**2)*(wf%n_ccsd_v) + (n_a_v**2)*(wf%n_o), &
             (n_a_o)*(n_a_v)*(wf%n_v) + (wf%n_ccsd_o)*(n_a_o)*(n_a_v)&
              + (wf%n_ccsd_o**2)*(wf%n_ccsd_v), &
             n_a_v + (wf%n_ccsd_o**2)*(wf%n_ccsd_v))

!
      batch_b = batching_index(wf%n_ccsd_v)
      call mem%batch_setup(batch_b, rec0, rec1)
!
!     Start looping over b-batches
!
      do current_b_batch = 1, batch_b%num_batches
!
!        Get batching limits for current b-batch
!
         call batch_b%determine_limits(current_b_batch)
!
         call mem%alloc(rho_aibj_batch, wf%n_ccsd_v, wf%n_ccsd_o, batch_b%length, wf%n_ccsd_o)
         call zero_array(rho_aibj_batch, (wf%n_ccsd_o**2)*(batch_b%length)*(wf%n_ccsd_v))
!
         call mem%alloc(g_kcbd, wf%n_o, wf%n_v, batch_b%length, wf%n_v)
!
         call wf%eri%get_eri_t1('ovvv', g_kcbd, 1, wf%n_o, 1, wf%n_v, &
                                                batch_b%first, batch_b%last, 1, wf%n_v)
!
         call mem%alloc(g_dkbc, n_a_v, n_a_o, batch_b%length, wf%n_v)
!
!$omp parallel do private (d, k, b, c)
         do d = 1, n_a_v
            do k = 1, n_a_o
               do b = 1, batch_b%length
                  do c = 1, wf%n_v
!
                     g_dkbc(d, k, b, c) = g_kcbd(k , c, b, d)
!
                  enddo
               enddo
            enddo
         enddo
!$omp end parallel do

!
         call wf%jacobian_ccsd_d2_2(rho_aibj_batch, c_ai, batch_b, g_dkbc)
!
         call mem%dealloc(g_dkbc, n_a_v, n_a_o, batch_b%length, wf%n_v)
!
         call mem%alloc(g_kcbd_3, n_a_o, n_a_v, batch_b%length, wf%n_v)
!
!$omp parallel do private (d, k, b, c)
         do d = 1, wf%n_v
            do b = 1, batch_b%length
               do c = 1, n_a_v
                  do k = 1, n_a_o
!
                     g_kcbd_3(k, c, b, d) = g_kcbd(k, c, b, d)
!
                  enddo
               enddo
            enddo
         enddo
!$omp end parallel do
!        
         call wf%jacobian_ccsd_d2_3(rho_aibj_batch, c_ai, batch_b, g_kcbd_3)
!
         call mem%dealloc(g_kcbd_3, n_a_o, n_a_v, batch_b%length, wf%n_v)
!
         call mem%alloc(L_ckbd, wf%n_v, wf%n_o, batch_b%length, wf%n_v)
         call zero_array(L_ckbd, (wf%n_v**2)*(wf%n_o)*(batch_b%length))
!
         call add_2134_to_1234(two, g_kcbd, L_ckbd, wf%n_v, wf%n_o, batch_b%length, wf%n_v)
         call add_2431_to_1234(-one, g_kcbd, L_ckbd, wf%n_v, wf%n_o, batch_b%length, wf%n_v)
!
         call mem%dealloc(g_kcbd, wf%n_o, wf%n_v, batch_b%length, wf%n_v)
!
         call mem%alloc(L_ckbd_4, n_a_v, n_a_o, batch_b%length, wf%n_v)
!
!$omp parallel do private (d, k, b, c)
         do d = 1, wf%n_v
            do b = 1, batch_b%length
               do k = 1, n_a_o
                  do c = 1, n_a_v
!
                     L_ckbd_4(c, k, b, d) = L_ckbd(c, k, b, d)
!
                  enddo
               enddo
            enddo
         enddo
!$omp end parallel do
!
         call wf%jacobian_ccsd_d2_4(rho_aibj_batch, c_ai, batch_b, L_ckbd_4)
!
         call mem%dealloc(L_ckbd_4, n_a_v, n_a_o, batch_b%length, wf%n_v)
!
         call mem%alloc(L_ckbd_5, wf%n_v, wf%n_o, batch_b%length, n_a_v)
!
!$omp parallel do private (d, k, b, c)
         do d = 1, n_a_v
            do b = 1, batch_b%length
               do k = 1, wf%n_o
                  do c = 1, wf%n_v
!
                     L_ckbd_5(c, k, b, d) = L_ckbd(c, k, b, d)
!
                  enddo
               enddo
            enddo
         enddo
!$omp end parallel do
!
         call mem%dealloc(L_ckbd, wf%n_v, wf%n_o, batch_b%length, wf%n_v)
!
         call wf%jacobian_ccsd_d2_5(rho_aibj_batch, c_ai, batch_b, L_ckbd_5)
!
         call mem%dealloc(L_ckbd_5, wf%n_v, wf%n_o, batch_b%length, n_a_v) 
!
!$omp parallel do private(j, b, i, a) collapse(2)
         do j = 1, wf%n_ccsd_o
            do b = 1, batch_b%length
               do i = 1, wf%n_ccsd_o
                  do a = 1, wf%n_ccsd_v
!
                     rho_aibj(a, i, b + batch_b%first - 1, j) = &
                              rho_aibj(a, i, b + batch_b%first - 1, j) &
                                 + rho_aibj_batch(a, i, b, j)
!
                  enddo
               enddo
            enddo
         enddo
!$omp end parallel do
!
         call mem%dealloc(rho_aibj_batch, wf%n_ccsd_v, wf%n_ccsd_o, batch_b%length, wf%n_ccsd_o)
!   
      enddo
!
      call timer%turn_off()
!
   end subroutine jacobian_ccsd_d2_mlccsd
!
!
    module subroutine jacobian_ccsd_e2_mlccsd(wf, rho_aibj, c_aibj)
!!
!!    Jacobian CCSD E2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    rho_aibj^E2 = 2 sum_dlck x_bj,dl * L_kc,ld * c_ai,ck
!!                  - sum_dlck x_bj,dl * L_kc,ld * c_ak,ci
!!
!!                = 2 Y_ck,bj * c_ai,ck - Y_ck,bj * c_ak,ci
!!
!!                = Y_ck,bj  (2 c_ai,ck - c_ak,ci)
!!
!!
!!    INDEX RESTRICTIONS:
!!
!!       a, i, b, j are CCSD indices
!!
!!       d, l, c, k are CC2 + CCSD indices
!!
      implicit none
!
      class(mlccsd) :: wf
!
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o) :: rho_aibj
      real(dp), dimension(wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o, &
                         wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o), intent(in) :: c_aibj
!
      type(timings), allocatable :: timer
!
      integer :: n_a_o, n_a_v 
!
      real(dp), dimension(:,:,:,:), allocatable :: c_aick
      real(dp), dimension(:,:,:,:), allocatable :: X_ckbj
!
      integer :: a, i, c, k
!
      timer = timings('Jacobian MLCCSD E2', pl='v')
      call timer%turn_on()
!
      n_a_o = wf%n_ccsd_o + wf%n_cc2_o
      n_a_v = wf%n_ccsd_v + wf%n_cc2_v
!
      call mem%alloc(X_ckbj, n_a_v, n_a_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
      wf%jacobian_e2_intermediate = sequential_file('jacobian_e2_intermediate_ccsd')
      call wf%jacobian_e2_intermediate%open_('read', 'rewind')
!
      call wf%jacobian_e2_intermediate%read_(X_ckbj, n_a_o*n_a_v*(wf%n_ccsd_o)*(wf%n_ccsd_v))
!
      call wf%jacobian_e2_intermediate%close_('keep')
!
!     c_aick = (2 c_ai,ck - c_ak,ci)
!
      call mem%alloc(c_aick, wf%n_ccsd_v, wf%n_ccsd_o, n_a_v, n_a_o)
!
!$omp parallel do private (k, c, i, a) collapse(2)
      do k = 1, n_a_o
         do c = 1, n_a_v
            do i = 1, wf%n_ccsd_o
               do a = 1, wf%n_ccsd_v
!
                  c_aick(a, i, c, k) =- c_aibj(a, k, c, i) &
                                       + two*c_aibj(a, i, c, k)
                                       
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call dgemm('N', 'N',                      &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),  &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),  &
                  (n_a_o)*(n_a_v),              &
                  one,                          &
                  c_aick,                       & 
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),  &
                  X_ckbj,                       &
                  (n_a_o)*(n_a_v),              &
                  one,                          &
                  rho_aibj,                     &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v))
!
      call mem%dealloc(X_ckbj, n_a_v, n_a_o, wf%n_ccsd_v, wf%n_ccsd_o)
      call mem%dealloc(c_aick, wf%n_ccsd_v, wf%n_ccsd_o, n_a_v, n_a_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_ccsd_e2_mlccsd
!
!
   module subroutine jacobian_ccsd_f2_mlccsd(wf, rho_aibj, c_aibj)
!!
!!    Jacobian MLCCSD F2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!       rho_aibj^F2 = - sum_ckdl x_ai,dj * L_kc,ld * c_bl,ck
!!                     - sum_ckdl x_ai_bl * L_kc,ld * c_ck,dj
!!
!!    L_kc,ld = 2*g_kc,ld - g_kd,lc = 2*g_kcld(kc,ld) - 2*g_kcld(kd,lc)
!!
!!    INDEX RESTRICTIONS:
!!
!!    a, i, b, j are CCSD indices
!!
!!    c, k, d, l are CCSD + CC2 indices
!!
      implicit none
!
      class(mlccsd) :: wf
!
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o) :: rho_aibj
      real(dp), dimension(wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o, &
                         wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o), intent(in) :: c_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: g_kcld
      real(dp), dimension(:,:,:,:), allocatable :: L_ckdl
      real(dp), dimension(:,:,:,:), allocatable :: x_aijd
      real(dp), dimension(:,:,:,:), allocatable :: x_aibl
      real(dp), dimension(:,:,:,:), allocatable :: rho_aijb
      real(dp), dimension(:,:), allocatable :: Y_db
      real(dp), dimension(:,:), allocatable :: Z_jl
!
      integer :: n_a_o, n_a_v 
!
      type(timings), allocatable :: timer
!
      integer :: d, j, i, a, b, l, ai, dj, aidj, bl, aibl
!
      timer = timings('Jacobian MLCCSD F2', pl='v')
      call timer%turn_on()
!
      n_a_o = wf%n_ccsd_o + wf%n_cc2_o
      n_a_v = wf%n_ccsd_v + wf%n_cc2_v
!
!     Construct L_kcld ordered as L_ckdl
!
      call mem%alloc(g_kcld, n_a_o, n_a_v, n_a_o, n_a_v)
!
      call wf%eri%get_eri_t1('ovov', g_kcld, 1, n_a_o, 1, n_a_v, 1, n_a_o, 1, n_a_v)
!
      call mem%alloc(L_ckdl, n_a_v, n_a_o, n_a_v, n_a_o)
      call zero_array(L_ckdl, (n_a_v**2)*(n_a_o**2))
!
      call add_2341_to_1234(-one, g_kcld, L_ckdl, n_a_v, n_a_o, n_a_v, n_a_o)
      call add_2143_to_1234(two, g_kcld, L_ckdl, n_a_v, n_a_o, n_a_v, n_a_o)
!
      call mem%dealloc(g_kcld, n_a_o, n_a_v, n_a_o, n_a_v)
!
!     :: Term 1: - sum_ckdl x_aidj * L_kcld * c_blck
!
!     Y_db = sum_clk L_dlck * c_blck
!
      call mem%alloc(Y_db, n_a_v, wf%n_ccsd_v)
!
      call dgemm('N', 'T',             &
                  n_a_v,               &
                  wf%n_ccsd_v,         &
                  (n_a_o**2)*(n_a_v),  &
                  one,                 &
                  L_ckdl,              & ! L_dlck
                  n_a_v,               &
                  c_aibj,              & ! c_blck
                  n_a_v,               &
                  zero,                &
                  Y_db,                & 
                  n_a_v)
!
      call mem%alloc(x_aijd, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_o, n_a_v)
!
!$omp parallel do private(d, j, i, a, dj, ai, aidj) collapse(2)
      do d = 1, n_a_v
         do j = 1, wf%n_ccsd_o
!
            dj = n_a_v*(j - 1) + d
!
            do i = 1, wf%n_ccsd_o 
               do a = 1, wf%n_ccsd_v
!
                  ai = n_a_v*(i - 1) + a
!
                  aidj = max(ai,dj)*(max(ai,dj)-3)/2 + ai + dj
                  x_aijd(a, i, j, d) = wf%x2(aidj)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%alloc(rho_aijb, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_o, wf%n_ccsd_v)
!
!     rho_aij_b = sum_d t_aijd*Y_d_b
!
      call dgemm('N','N',                             &
                  ((wf%n_ccsd_o)**2)*(wf%n_ccsd_v),   &
                  wf%n_ccsd_v,                        &
                  n_a_v,                              &
                  -one,                               &
                  x_aijd,                             &
                  ((wf%n_ccsd_o)**2)*(wf%n_ccsd_v),   &
                  Y_db,                               &
                  n_a_v,                              &
                  zero,                               &
                  rho_aijb,                           &
                  ((wf%n_ccsd_o)**2)*(wf%n_ccsd_v))
!
      call mem%dealloc(x_aijd, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_o, n_a_v)
      call mem%dealloc(Y_db, n_a_v, wf%n_ccsd_v)
!
      call add_1243_to_1234(one, rho_aijb, rho_aibj, &
                  wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
      call mem%dealloc(rho_aijb, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_o, wf%n_ccsd_v)
!
!     :: Term 2: - sum_ckdl x_aibl * L_kcld * c_ckdj ::
!
      call mem%alloc(Z_jl, wf%n_ccsd_o, n_a_o)
!
!     Z_jl = sum_ckd L_ckdl c_ckdj
!
      call dgemm('T', 'N',             &
                  wf%n_ccsd_o,         &
                  n_a_o,               &
                  (n_a_v**2)*(n_a_o),  &
                  one,                 &
                  c_aibj,              &  ! c_ckdj
                  (n_a_v**2)*(n_a_o),  &
                  L_ckdl,              & 
                  (n_a_v**2)*(n_a_o),  &
                  zero,                &
                  Z_jl,                &
                  wf%n_ccsd_o)
!
      call mem%dealloc(L_ckdl, n_a_v, n_a_o, n_a_v, n_a_o)
!
      call mem%alloc(x_aibl, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, n_a_o)
!
!$omp parallel do private (l, b, i, a, bl, ai, aibl) collapse(2)
      do l = 1, n_a_o
         do b = 1, wf%n_ccsd_v
!
            bl = n_a_v*(l-1) + b
!
            do i = 1, wf%n_ccsd_o
               do a = 1, wf%n_ccsd_v
!
                  ai = n_a_v*(i-1) + a
                  aibl = max(ai,bl)*(max(ai,bl)-3)/2 + ai + bl
!
                  x_aibl(a,i,b,l) = wf%x2(aibl)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call dgemm('N','T',                             &
                  ((wf%n_ccsd_v)**2)*(wf%n_ccsd_o),   &
                  wf%n_ccsd_o,                        &
                  n_a_o,                              &
                  -one,                               &
                  x_aibl,                             & 
                  ((wf%n_ccsd_v)**2)*(wf%n_ccsd_o),   &
                  Z_jl,                               & 
                  wf%n_ccsd_o,                        &
                  one,                                &
                  rho_aibj,                           & 
                  ((wf%n_ccsd_v)**2)*(wf%n_ccsd_o))
!
      call mem%dealloc(x_aibl, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, n_a_o)
      call mem%dealloc(Z_jl, wf%n_ccsd_o, n_a_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_ccsd_f2_mlccsd
!
!
   module subroutine jacobian_ccsd_g2_mlccsd(wf, rho_aibj, c_aibj)
!!
!!    Jacobian MLCCSD G2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    rho_aibj^G2 =  - sum_ckdl x_bl,dj * L_kc,ld * c_ai,ck
!!                   - sum_ckdl x_ck_bl * L_kc,ld * c_ai,dj
!!                   - sum_ckld x_ck,dj * L_kc,ld * c_ai,bl
!!
!!    L_kc,ld = 2*g_kc,ld - g_kd,lc = 2*g_kcld(kc,ld) - 2*g_kcld(kd,lc)
!!
!!    INDEX RESTRICTIONS:
!!
!!    a, i, b, j are CCSD indices
!!
!!    c, k, d, l are CCSD + CC2 indices
!!
      implicit none
!
      class(mlccsd) :: wf
!
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o) :: rho_aibj
      real(dp), dimension(wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o, &
                         wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o), intent(in) :: c_aibj
!
!
      integer :: n_a_o, n_a_v 
!
      real(dp), dimension(:,:,:,:), allocatable :: X_ckbj
      real(dp), dimension(:,:,:,:), allocatable :: c_aick
      real(dp), dimension(:,:,:,:), allocatable :: c_djai
      real(dp), dimension(:,:,:,:), allocatable :: c_aibl
!
      real(dp), dimension(:,:), allocatable :: X_db
      real(dp), dimension(:,:), allocatable :: X_lj
!
      type(timings), allocatable :: timer
!
      integer :: a, i, b, j, c, k, d, l
!
      timer = timings('Jacobian MLCCSD G2', pl='v')
      call timer%turn_on()
!
      n_a_o = wf%n_ccsd_o + wf%n_cc2_o
      n_a_v = wf%n_ccsd_v + wf%n_cc2_v
!
      call mem%alloc(X_ckbj, n_a_v, n_a_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
      call wf%jacobian_g2_intermediate_vovo%open_('read', 'rewind')
!
      call wf%jacobian_g2_intermediate_vovo%read_(X_ckbj, n_a_o*n_a_v*(wf%n_ccsd_o)*(wf%n_ccsd_v))
!
      call wf%jacobian_g2_intermediate_vovo%close_('keep')
!
!     :: Term 1: - sum_ckdl x_bl,dj * L_kc,ld * c_ai,ck  ::
!
      call mem%alloc(c_aick, wf%n_ccsd_v, wf%n_ccsd_o, n_a_v, n_a_o)
!
!$omp parallel do private (k, c, i, a) collapse(2)
      do k = 1, n_a_o
         do c = 1, n_a_v
            do i = 1, wf%n_ccsd_o
               do a = 1, wf%n_ccsd_v
!
                  c_aick(a, i, c, k) = c_aibj(a, i, c, k)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call dgemm('N', 'N',                      &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),  &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),  &
                  (n_a_o)*(n_a_v),              &
                  -one,                         &
                  c_aick,                       & 
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),  &
                  X_ckbj,                       & 
                  (n_a_o)*(n_a_v),              &
                  one,                          &
                  rho_aibj,                     &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v))
!
      call mem%dealloc(X_ckbj, n_a_v, n_a_o, wf%n_ccsd_v, wf%n_ccsd_o)
      call mem%dealloc(c_aick, wf%n_ccsd_v, wf%n_ccsd_o, n_a_v, n_a_o)
!
!     :: Term 2: - sum_ckdl x_blck * L_kcld * c_aidj = - sum_d X_db c_aidj
!
      call mem%alloc(X_db, n_a_v, wf%n_ccsd_v)
!
      call wf%jacobian_g2_intermediate_vv%open_('read', 'rewind')
!
      call wf%jacobian_g2_intermediate_vv%read_(X_db, n_a_v*(wf%n_ccsd_v))
!
      call wf%jacobian_g2_intermediate_vv%close_('keep')
!
      call mem%alloc(c_djai, n_a_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
!$omp parallel do private(i, a, j, d) collapse(2)
      do i = 1, wf%n_ccsd_o
         do a = 1, wf%n_ccsd_v
            do j = 1, wf%n_ccsd_o
               do d = 1, n_a_v
!
                  c_djai(d, j, a, i) = c_aibj(d, j, a, i)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call dgemm('T','N',                             &
                  wf%n_ccsd_v,                        &
                  ((wf%n_ccsd_o)**2)*(wf%n_ccsd_v),   &
                  n_a_v,                              &
                  -one,                               &
                  X_db,                               & 
                  n_a_v,                              &
                  c_djai,                             & 
                  n_a_v,                              &
                  one,                                &
                  rho_aibj,                           & ! we will symmetrize after
                  (wf%n_ccsd_v))
!
      call mem%dealloc(c_djai, n_a_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
      call mem%dealloc(X_db, n_a_v, wf%n_ccsd_v)
!
!     :: Term 3: - sum_ckld x_ckdj * L_kcld * c_aibl = -sum_l X_lj c_aibl
!
      call mem%alloc(X_lj, n_a_o, wf%n_ccsd_o)
      call wf%jacobian_g2_intermediate_oo%open_('read', 'rewind')
!
      call wf%jacobian_g2_intermediate_oo%read_(X_lj, n_a_o*(wf%n_ccsd_o))
!
      call wf%jacobian_g2_intermediate_oo%close_('keep')
!
      call mem%alloc(c_aibl, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, n_a_o)
!
!$omp parallel do private (l, b, i, a) collapse(2)
      do l = 1, n_a_o
         do b = 1, wf%n_ccsd_v
            do i = 1, wf%n_ccsd_o
               do a = 1, wf%n_ccsd_v
!
                  c_aibl(a, i, b, l) = c_aibj(a, i, b, l)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call dgemm('N','N',                             &
                  ((wf%n_ccsd_v)**2)*(wf%n_ccsd_o),   &
                  wf%n_ccsd_o,                        &
                  n_a_o,                              &
                  -one,                               &
                  c_aibl,                             &
                  ((wf%n_ccsd_v)**2)*(wf%n_ccsd_o),   &
                  X_lj,                               & 
                  n_a_o,                              &
                  one,                                &
                  rho_aibj,                           &
                  ((wf%n_ccsd_v)**2)*(wf%n_ccsd_o))
!
      call mem%dealloc(X_lj, n_a_o, wf%n_ccsd_o)
      call mem%dealloc(c_aibl, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, n_a_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_ccsd_g2_mlccsd
!
!
   module subroutine jacobian_ccsd_h2_mlccsd(wf, rho_aibj, c_aibj)
!!
!!    Jacobian CCSD H2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!       rho_aibj^H2 =   sum_ckdl x_ci,ak * g_kc,ld * c_bl,dj
!!                     + sum_ckdl x_cj,al * g_kc,ld * c_bk,di
!!
!!    INDEX RESTRICTIONS
!!
!!    a, i, b, j are CCSD indices
!!    c, k, d, l are CC2 + CCSD indices
!!
!!
      implicit none
!
      class(mlccsd) :: wf

!
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o) :: rho_aibj
      real(dp), dimension(wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o, &
                         wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o), intent(in) :: c_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: X_ajdk
      real(dp), dimension(:,:,:,:), allocatable :: X_aidl
      real(dp), dimension(:,:,:,:), allocatable :: c_dlbj
      real(dp), dimension(:,:,:,:), allocatable :: c_dkbi
      real(dp), dimension(:,:,:,:), allocatable :: rho_ajbi
!
      integer :: i, b, j, k, d, l
!
      type(timings), allocatable :: timer
!
      integer :: n_a_o, n_a_v
!
      timer = timings('Jacobian MLCCSD H2', pl='v')
      call timer%turn_on()
!
      n_a_o = wf%n_ccsd_o + wf%n_cc2_o
      n_a_v = wf%n_ccsd_v + wf%n_cc2_v
!
!     :: Term 1: sum_ckld x_ci,ak * g_kc,ld * c_bl,dj ::
!
      call mem%alloc(X_aidl, wf%n_ccsd_v, wf%n_ccsd_o, n_a_v, n_a_o)
!
      call wf%jacobian_h2_intermediate_vovo_1%open_('read', 'rewind')
!
      call wf%jacobian_h2_intermediate_vovo_1%read_(X_aidl, n_a_v*n_a_o*(wf%n_ccsd_o)*(wf%n_ccsd_v))
!
      call wf%jacobian_h2_intermediate_vovo_1%close_('keep')
!
      call mem%alloc(c_dlbj, n_a_v, n_a_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
!$omp parallel do private(j, b, d, l) collapse(2)
      do j = 1, wf%n_ccsd_o
         do b = 1, wf%n_ccsd_v
            do d = 1, n_a_v
               do l = 1, n_a_o
!
                  c_dlbj(d, l, b, j) = c_aibj(d, j, b, l)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call dgemm('N', 'N',                      &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),  &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),  &
                  (n_a_o)*(n_a_v),              &
                  one,                          &
                  X_aidl,                       & 
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),  &
                  c_dlbj,                       &
                  (n_a_o)*(n_a_v),              &
                  one,                          &
                  rho_aibj,                     &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v))
!
      call mem%dealloc(c_dlbj, n_a_v, n_a_o, wf%n_ccsd_v, wf%n_ccsd_o)
      call mem%dealloc(X_aidl, wf%n_ccsd_v, wf%n_ccsd_o, n_a_v, n_a_o)
!
!     :: Term 2: sum_ckdl t_cjal * g_kcld * c_bkdi
!
      call mem%alloc(X_ajdk, wf%n_ccsd_v, wf%n_ccsd_o, n_a_v, n_a_o)
!
      call wf%jacobian_h2_intermediate_vovo_2%open_('read', 'rewind')
!
      call wf%jacobian_h2_intermediate_vovo_2%read_(X_ajdk, n_a_v*n_a_o*(wf%n_ccsd_o)*(wf%n_ccsd_v))
!
      call wf%jacobian_h2_intermediate_vovo_2%close_('keep')
!
!     Reorder c_bkdi as c_dkbi
!
      call mem%alloc(c_dkbi, n_a_v, n_a_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
!$omp parallel do private(i, b, d, k) collapse(2)
      do i = 1, wf%n_ccsd_o
         do b = 1, wf%n_ccsd_v
            do d = 1, n_a_v
               do k = 1, n_a_o
!
                  c_dkbi(d, k, b, i) = c_aibj(b, k, d, i)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%alloc(rho_ajbi, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
!     rho_ajbi = sum_kd  Y_ajkd * c_kdbi
!
      call dgemm('N', 'N',                      &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),  &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),  &
                  (n_a_o)*(n_a_v),              &
                  one,                          &
                  X_ajdk,                       & 
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),  &
                  c_dkbi,                       &
                  (n_a_o)*(n_a_v),              &
                  zero,                         &
                  rho_ajbi,                     &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v))
!
      call mem%dealloc(c_dkbi, n_a_v, n_a_o, wf%n_ccsd_v, wf%n_ccsd_o)
      call mem%dealloc(X_ajdk, wf%n_ccsd_v, wf%n_ccsd_o, n_a_v, n_a_o)
!
!     Reorder into rho_aibj
!
      call add_1432_to_1234(one, rho_ajbi, rho_aibj, &
                        wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
      call mem%dealloc(rho_ajbi, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_ccsd_h2_mlccsd
!
!
   module subroutine jacobian_ccsd_i2_1_mlccsd(wf, rho_aibj, c_aibj)
!!
!!    Jacobian CCSD I2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    rho_aibj^I2 =  sum_c F_bc * c_ai,cj - sum_k F_kj * c_ai,bk
!!
!!    INDEX RESTRICTIONS:
!!
!!       a, i, b, j are CCSD indices
!!
!!       c is CC2 index 
!!
      implicit none
!
      class(mlccsd) :: wf

!
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o) :: rho_aibj
      real(dp), dimension(wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o, &
                         wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o), intent(in) :: c_aibj
!
      integer :: n_a_o, n_a_v
!
      real(dp), dimension(:,:,:,:), allocatable :: c_cjai, c_aibk
!
      integer :: a, i, j, c, k, b
!
      n_a_o = wf%n_ccsd_o + wf%n_cc2_o
      n_a_v = wf%n_ccsd_v + wf%n_cc2_v
!
      call mem%alloc(c_cjai, n_a_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
!$omp parallel do private (i, a, j, c) collapse(2)
      do i = 1, wf%n_ccsd_o
         do a = 1, wf%n_ccsd_v
            do j = 1, wf%n_ccsd_o
               do c = 1, n_a_v
!
                  c_cjai(c, j, a, i) = c_aibj(c, j, a, i)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call dgemm('N', 'N',                                     &
                  wf%n_ccsd_v,                                 &
                  (wf%n_ccsd_o**2)*(wf%n_ccsd_v),              &
                  n_a_v,                                       &
                  one,                                         &
                  wf%fock_ab,                                  & ! F_bc
                  wf%n_v,                                      &
                  c_cjai,                                      & 
                  n_a_v,                                       &
                  one,                                         &
                  rho_aibj,                                    & ! We will symmetrize after
                  wf%n_ccsd_v)
!
      call mem%dealloc(c_cjai, n_a_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
      call mem%alloc(c_aibk, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, n_a_o)
!
!$omp parallel do private (k, b, i, a) collapse(2)
      do k = 1, n_a_o
         do b = 1, wf%n_ccsd_v
            do i = 1, wf%n_ccsd_o
               do a = 1, wf%n_ccsd_v
!
                  c_aibk(a, i, b, k) = c_aibj(a, i, b, k)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call dgemm('N', 'N',                         &
                  (wf%n_ccsd_v**2)*(wf%n_ccsd_o),  &
                  wf%n_ccsd_o,                     &
                  n_a_o,                           &
                  -one,                            &
                  c_aibk,                          & ! c_aibk
                  (wf%n_ccsd_v**2)*(wf%n_ccsd_o),  &
                  wf%fock_ij,                      &
                  wf%n_o,                          & ! F_kj
                  one,                             &
                  rho_aibj,                        &
                  (wf%n_ccsd_v**2)*(wf%n_ccsd_o))      
!
      call mem%dealloc(c_aibk, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, n_a_o)
!
   end subroutine jacobian_ccsd_i2_1_mlccsd
!
!
   module subroutine jacobian_ccsd_i2_2_mlccsd(wf, rho_aibj, c_aibj)
!!
!!    Jacobian CCSD I2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    rho_aibj^I2 =  sum_ck L_bj,kc*c_ai,ck - sum_ck ( g_kc,bj*c_ak,ci + g_ki,bc*c_ak,cj )
!!
!!    INDEX RESTRICTIONS:
!!
!!       a, i, b, j are CCSD indices
!!
!!       c, k is CC2 + CCSD index 
!!
      implicit none
!
      class(mlccsd) :: wf
!
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o) :: rho_aibj
      real(dp), dimension(wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o, &
                         wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o), intent(in) :: c_aibj
!
      integer :: n_a_o, n_a_v
!
      integer :: a, i, c, k, j
!
      real(dp), dimension(:,:,:,:), allocatable :: g_bjkc
      real(dp), dimension(:,:,:,:), allocatable :: Y_kcai
      real(dp), dimension(:,:,:,:), allocatable :: g_bckj
      real(dp), dimension(:,:,:,:), allocatable :: g_ckbj
      real(dp), dimension(:,:,:,:), allocatable :: c_aick
      real(dp), dimension(:,:,:,:), allocatable :: c_ajck
      real(dp), dimension(:,:,:,:), allocatable :: rho_ajbi
!
      n_a_o = wf%n_ccsd_o + wf%n_cc2_o
      n_a_v = wf%n_ccsd_v + wf%n_cc2_v
!
!     sum_ck L_bj,kc*c_ai,ck - sum_ck ( g_kc,bj*c_ak,ci + g_ki,bc*c_ak,cj ) ::
!
!     sum_ck ( g_bj,kc*(2*c_ai,ck - c_ak,ci) - g_bc,kj*c_ai,ck - g_ki,bc*c_ak,cj )
!     sum_ck ( g_bj,kc*Y_kcai - g_bc,kj*c_ai,ck - g_ki,bc*c_ak,cj )   
!
!     Construct g_bj,kc
!
      call mem%alloc(g_bjkc, wf%n_ccsd_v, wf%n_ccsd_o, n_a_o, n_a_v)
!
      call wf%eri%get_eri_t1('voov', g_bjkc, 1, wf%n_ccsd_v, 1, wf%n_ccsd_o, 1, n_a_o, 1, n_a_v)
!
      call mem%alloc(Y_kcai, n_a_o, n_a_v, wf%n_ccsd_v, wf%n_ccsd_o)
!
!     Y_kcai = (2*c_ai,ck - c_ak,ci) 
!
!$omp parallel do private (i, a, c, k) collapse(2)
      do i = 1, wf%n_ccsd_o
         do a = 1, wf%n_ccsd_v
            do c = 1, n_a_v
               do k = 1, n_a_o
!
                  Y_kcai(k, c, a, i) = two*c_aibj(a , i, c, k) &
                                         - c_aibj(a, k, c, i)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call dgemm('N', 'N',                      &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),  &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),  &
                  (n_a_o)*(n_a_v),              &
                  one,                          &
                  g_bjkc,                       & 
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),  &
                  Y_kcai,                       & 
                  (n_a_o)*(n_a_v),              &
                  one,                          & 
                  rho_aibj,                     & ! We will symmeterize after
                  (wf%n_ccsd_o)*(wf%n_ccsd_v))
!
      call mem%dealloc(Y_kcai, n_a_o, n_a_v, wf%n_ccsd_v, wf%n_ccsd_o)
      call mem%dealloc(g_bjkc, wf%n_ccsd_v, wf%n_ccsd_o, n_a_o, n_a_v)

      call mem%alloc(g_bckj, wf%n_ccsd_v, n_a_v, n_a_o, wf%n_ccsd_o)
!
      call wf%eri%get_eri_t1('vvoo', g_bckj, 1, wf%n_ccsd_v, 1, n_a_v, 1, n_a_o, 1, wf%n_ccsd_o) 
!
      call mem%alloc(g_ckbj, n_a_v, n_a_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
      call sort_1234_to_2314(g_bckj, g_ckbj, wf%n_ccsd_v, n_a_v, n_a_o, wf%n_ccsd_o)
!
      call mem%dealloc(g_bckj, wf%n_ccsd_v, n_a_v, n_a_o, wf%n_ccsd_o)
!
!     rho_aibj += - sum_ck c_aick * g_ckbj
!
      call mem%alloc(c_aick, wf%n_ccsd_v, wf%n_ccsd_o, n_a_v, n_a_o)
!
!$omp parallel do private (k, c, i, a) collapse(2)
      do k = 1, n_a_o
         do c = 1, n_a_v
            do i = 1, wf%n_ccsd_o 
               do a = 1, wf%n_ccsd_v
!
                  c_aick(a, i, c, k) = c_aibj(a, i, c, k)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call dgemm('N', 'N',                      &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),  &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),  &
                  (n_a_o)*(n_a_v),              &
                  -one,                         &
                  c_aick,                       & 
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),  &
                  g_ckbj,                       & 
                  (n_a_o)*(n_a_v),              &
                  one,                          &
                  rho_aibj,                     & 
                  (wf%n_ccsd_o)*(wf%n_ccsd_v))
!
      call mem%dealloc(c_aick, wf%n_ccsd_v, wf%n_ccsd_o, n_a_v, n_a_o)
!
!     Reorder  c_ak,cj to c_aj_ck
!
      call mem%alloc(c_ajck, wf%n_ccsd_v, wf%n_ccsd_o, n_a_v, n_a_o)
!
!$omp parallel do private(k, c, j, a) collapse(2)
      do k = 1, n_a_o
         do c = 1, n_a_v
            do j = 1, wf%n_ccsd_o 
               do a = 1, wf%n_ccsd_v
!
                  c_ajck(a, j, c, k) = c_aibj(a, k, c, j)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%alloc(rho_ajbi, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
      call dgemm('N', 'N',                      &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),  &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),  &
                  (n_a_o)*(n_a_v),              &
                  -one,                         &
                  c_ajck,                       & ! c_aj_ck
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),  &
                  g_ckbj,                       &  ! g_ck_bi
                  (n_a_o)*(n_a_v),              &
                  zero,                         &
                  rho_ajbi,                     & ! rho_aj_bi
                  (wf%n_ccsd_o)*(wf%n_ccsd_v))
!
      call mem%dealloc(g_ckbj, n_a_v, n_a_o, wf%n_ccsd_v, wf%n_ccsd_o)
      call mem%dealloc(c_ajck, wf%n_ccsd_v, wf%n_ccsd_o, n_a_v, n_a_o)
!
      call add_1432_to_1234(one, rho_ajbi, rho_aibj, &
               wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
      call mem%dealloc(rho_ajbi, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
   end subroutine jacobian_ccsd_i2_2_mlccsd
!
!      
   module subroutine jacobian_ccsd_i2_mlccsd(wf, rho_aibj, c_aibj)
!!
!!    Jacobian CCSD I2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    rho_aibj^I2 =  sum_c F_bc * c_ai,cj - sum_k F_kj * c_ai,bk
!!                   + sum_ck L_bj,kc * c_ai,ck
!!                   - sum_ck ( g_kc,bj * c_ak,ci + g_ki,bc * c_ak,cj )
!!
!!
      implicit none
!
      class(mlccsd) :: wf
!
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o) :: rho_aibj
      real(dp), dimension(wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o, &
                          wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o), intent(in) :: c_aibj
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian MLCCSD I2', pl='v')

      call timer%turn_on()
!
      call wf%jacobian_ccsd_i2_1(rho_aibj, c_aibj)
      call wf%jacobian_ccsd_i2_2(rho_aibj, c_aibj)
!
      call timer%turn_off()
!
   end subroutine jacobian_ccsd_i2_mlccsd
!
!
   module subroutine jacobian_ccsd_j2_mlccsd(wf, rho_abij, c_abij)
!!
!!    Jacobian MLCCSD J2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!       rho_abij^J2 =    sum_ckld x_ci,dj * g_kc,ld * c_ak,bl
!!                      + sum_ckdl x_ak,bl * g_kc,ld * c_ci,dj
!!
      implicit none
!
      class(mlccsd) :: wf
!
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_o) :: rho_abij
      real(dp), dimension(wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_v + wf%n_cc2_v, &
                          wf%n_ccsd_o + wf%n_cc2_o, wf%n_ccsd_o + wf%n_cc2_o), intent(in) :: c_abij
!
      type(timings), allocatable :: timer
!
      real(dp), dimension(:,:,:,:), allocatable :: g_klcd
      real(dp), dimension(:,:,:,:), allocatable :: x_abkl
      real(dp), dimension(:,:,:,:), allocatable :: c_cdij, c_abkl
      real(dp), dimension(:,:,:,:), allocatable :: X_klij
!
      integer :: i, j, d, c, a, b, k, l, ak, bl, akbl
!
      integer :: n_a_o, n_a_v
!
      timer = timings('Jacobian MLCCSD J2', pl='v')
      call timer%turn_on()
!
      n_a_o = wf%n_ccsd_o + wf%n_cc2_o
      n_a_v = wf%n_ccsd_v + wf%n_cc2_v
!
      call mem%alloc(X_klij, n_a_o, n_a_o, wf%n_ccsd_o, wf%n_ccsd_o)
!
      call wf%jacobian_j2_intermediate_oooo%open_('read', 'rewind')
      call wf%jacobian_j2_intermediate_oooo%read_(X_klij, n_a_o**2 * (wf%n_ccsd_o**2))
      call wf%jacobian_j2_intermediate_oooo%close_('keep')
!
      call mem%alloc(c_abkl, wf%n_ccsd_v, wf%n_ccsd_v, n_a_o, n_a_o)
!
!$omp parallel do private(l, k, b, a) collapse(2)
      do l = 1, n_a_o
         do k = 1, n_a_o
            do b = 1,wf%n_ccsd_v
               do a = 1, wf%n_ccsd_v
!
                  c_abkl(a, b, k, l) = c_abij(a, b, k, l)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
!     rho_abij += c_ab_kl * X_kl_ij
!
      call dgemm('N', 'N',          &
                  (wf%n_ccsd_v)**2, &
                  (wf%n_ccsd_o)**2, &
                  (n_a_o)**2,       &
                  one,              &
                  c_abkl,           & 
                  (wf%n_ccsd_v)**2, &
                  X_klij,           & ! X_kl_ij
                  (n_a_o)**2,       &
                  one,              &
                  rho_abij,         & 
                  (wf%n_ccsd_v)**2)
!
      call mem%dealloc(c_abkl, wf%n_ccsd_v, wf%n_ccsd_v, n_a_o, n_a_o)
!
      call mem%alloc(g_klcd, n_a_o, n_a_o, n_a_v, n_a_v)
!
      call wf%jacobian_j2_intermediate_oovv%open_('read', 'rewind')
      call wf%jacobian_j2_intermediate_oovv%read_(g_klcd, n_a_v**2 * n_a_o**2)
      call wf%jacobian_j2_intermediate_oovv%close_('keep')
!
!     X_kl_ij = g_kl_cd * c_cd_ij
!
      call mem%alloc(c_cdij, n_a_v, n_a_v, wf%n_ccsd_o, wf%n_ccsd_o)
!
!$omp parallel do private(j, i, c, d) collapse(2)
      do j = 1, wf%n_ccsd_o
         do i = 1, wf%n_ccsd_o
            do d = 1, n_a_v
               do c = 1, n_a_v
!
                  c_cdij(c, d, i, j) = c_abij(c, d, i, j)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call dgemm('N', 'N',          &
                  (n_a_o)**2,       &
                  (wf%n_ccsd_o)**2, &
                  (n_a_v)**2,       &
                  one,              &
                  g_klcd,           & ! g_kl_cd
                  (n_a_o)**2,       &
                  c_cdij,           & ! c_cd_ij
                  (n_a_v)**2,       &
                  zero,             &
                  X_klij,           & ! X_kl_ij
                  (n_a_o)**2)
!
      call mem%dealloc(c_cdij, n_a_v, n_a_v, wf%n_ccsd_o, wf%n_ccsd_o)
      call mem%dealloc(g_klcd, n_a_o, n_a_o, n_a_v, n_a_v)
!
!     rho_abij += x_abkl * X_klij
!
      call mem%alloc(x_abkl, wf%n_ccsd_v, wf%n_ccsd_v, n_a_o, n_a_o)
!
!$omp parallel do private(l,k,b, a, ak, bl, akbl) collapse(2)
      do l = 1, n_a_o
         do k = 1, n_a_o
            do b = 1, wf%n_ccsd_v
!
               bl = n_a_v*(l - 1) + b
!
               do a = 1, wf%n_ccsd_v
!
                  ak = n_a_v*(k - 1) + a
                  akbl = max(ak,bl)*(max(ak,bl)-3)/2 + ak + bl
!
                  x_abkl(a, b, k, l) = wf%x2(akbl)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call dgemm('N', 'N',             &
                  (wf%n_ccsd_v)**2,    &
                  (wf%n_ccsd_o)**2,    &
                  (n_a_o)**2,          &
                  one,                 &
                  x_abkl,              & 
                  (wf%n_ccsd_v)**2,    &
                  X_klij,              & ! X_kl_ij
                  (n_a_o)**2,          &
                  one,                 &
                  rho_abij,            & ! rho_ab_ij
                  (wf%n_ccsd_v)**2)
!
      call mem%dealloc(X_klij, n_a_o, n_a_o, wf%n_ccsd_o, wf%n_ccsd_o)
      call mem%dealloc(x_abkl, wf%n_ccsd_v, wf%n_ccsd_v, n_a_o, n_a_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_ccsd_j2_mlccsd
!
!
   module subroutine jacobian_ccsd_k2_mlccsd(wf, rho_abij, c_abij)
!!
!!    Jacobian MLCCSD K2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    rho_abij^K2 =    sum_kl g_kilj * c_akbl
!!                       + sum_cd g_acbd * c_cidj (in omega term)
!!
!!    For the last term we batch over a and b and
!!    add each batch to rho_aibj
!!
      implicit none
!
      class(mlccsd) :: wf
!
      real(dp), dimension(wf%n_ccsd_v, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_o) :: rho_abij
      real(dp), dimension(wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_v + wf%n_cc2_v, &
                         wf%n_ccsd_o + wf%n_cc2_o, wf%n_ccsd_o + wf%n_cc2_o), intent(in) :: c_abij
!
      real(dp), dimension(:,:,:,:), allocatable :: g_kilj
      real(dp), dimension(:,:,:,:), allocatable :: g_klij
      real(dp), dimension(:,:,:,:), allocatable :: c_abkl
!
      integer :: a, b, k, l
!
      type(timings), allocatable :: timer
!
      integer :: n_a_o, n_a_v
!
      timer = timings('Jacobian MLCCSD K2', pl='v')
      call timer%turn_on()
!
      n_a_o = wf%n_ccsd_o + wf%n_cc2_o
      n_a_v = wf%n_ccsd_v + wf%n_cc2_v
!
      call mem%alloc(g_kilj, n_a_o, wf%n_ccsd_o, n_a_o, wf%n_ccsd_o)
!
      call wf%eri%get_eri_t1('oooo', g_kilj, &
                             1, n_a_o,       &
                             1, wf%n_ccsd_o, &
                             1, n_a_o,       &
                             1, wf%n_ccsd_o)
!
!     Reorder g_kilj to g_klij
!
      call mem%alloc(g_klij, n_a_o, n_a_o, wf%n_ccsd_o, wf%n_ccsd_o)
!
      call sort_1234_to_1324(g_kilj, g_klij, n_a_o, wf%n_ccsd_o, n_a_o, wf%n_ccsd_o)
!
      call mem%dealloc(g_kilj, n_a_o, wf%n_ccsd_o, n_a_o, wf%n_ccsd_o)
!
!     rho_abij += sum_kl g_kilj * c_akbl = sum_kl c_abij(a,b,k,l) g_klij
!
      call mem%alloc(c_abkl, wf%n_ccsd_v, wf%n_ccsd_v, n_a_o, n_a_o)
!
!$omp parallel do private(l, k, a, b) collapse(2)
      do l = 1, n_a_o
         do k = 1, n_a_o
            do b = 1, wf%n_ccsd_v
               do a = 1, wf%n_ccsd_v
!
                  c_abkl(a, b, k, l) = c_abij(a, b, k, l)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call dgemm('N', 'N',          &
                  (wf%n_ccsd_v)**2, &
                  (wf%n_ccsd_o)**2, &
                  (n_a_o)**2,       &
                  one,              &
                  c_abkl,           & 
                  (wf%n_ccsd_v)**2, &
                  g_klij,           & ! g_kl_ij
                  (n_a_o)**2,       &
                  one,             &
                  rho_abij,         & ! rho_ab_ij
                  (wf%n_ccsd_v)**2)
!
      call mem%dealloc(c_abkl, wf%n_ccsd_v, wf%n_ccsd_v, n_a_o, n_a_o)
!
      call mem%dealloc(g_klij, n_a_o, n_a_o, wf%n_ccsd_o, wf%n_ccsd_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_ccsd_k2_mlccsd
!
!
   module subroutine save_jacobian_c2_intermediates_mlccsd(wf)
!!
!!    Save jacobian c2 intermediates
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2019
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    Constructs the intermediates for Jacobian C2: 
!!
!!       X_ljai = sum_ck g_ljkc x_ki^ac - sum_kc L_ljkc x_ik^ac
!!
!!          Index restrictions: a, i, j are CCSD indices, 
!!          c and k are CCSD + CC2 indices, and l is unrestricted 
!!
!!       X_kjbi = g_ljkc x_li^bc
!!
!!          Index restrictions: b, i, j are CCSD indices, 
!!          c and l are CCSD + CC2 indices, and k is unrestricted 
!!
!!    used in the c2-term. This is done only once in prepare_for_jacobian
!!    and the intermediate is stored in the file jacobian_c2_intermediate_oovo_1
!!    and jacobian_c2_intermediate_oovo_2 which are wf variables.
!!
      implicit none
!
      class(mlccsd) :: wf
!
      type(timings), allocatable :: timer
!
      real(dp), dimension(:,:,:,:), allocatable :: X_ljai 
      real(dp), dimension(:,:,:,:), allocatable :: X_kjbi 
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ljck, g_kjcl, L_ljck
!
      integer :: n_a_o, n_a_v
      integer :: j, l
!
      real(dp), dimension(:,:,:,:), allocatable :: x_ckai, g_ljkc
!
      integer :: i, a, k, c, ci, ak, ciak, ai, ck, ckai
!
      n_a_o = wf%n_ccsd_o + wf%n_cc2_o
      n_a_v = wf%n_ccsd_v + wf%n_cc2_v
!
      timer = timings('Jacobian MLCCSD C2 intermediates', pl='v')
      call timer%turn_on()
!
      call mem%alloc(g_ljkc, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call wf%eri%get_eri_t1('ooov', g_ljkc, 1, wf%n_o, 1, wf%n_o, 1, wf%n_o, 1, wf%n_v)
!
      call mem%alloc(x_ckai, n_a_v, n_a_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
!     Order x_ciak as x_ckai and restrict active indices (a, i) to CCSD space
!
!$omp parallel do private (i, a, k, c, ak, ci, ciak) collapse(2)
      do i = 1, wf%n_ccsd_o
         do a = 1, wf%n_ccsd_v
            do k = 1, n_a_o
!
               ak = n_a_v*(k-1) + a
!
               do c = 1, n_a_v
!
                  ci = n_a_v*(i - 1) + c
                  ciak = max(ci, ak)*(max(ci,ak)-3)/2 + ci + ak
!
                  x_ckai(c, k, a, i) = wf%x2(ciak)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
!     X_ljai = sum_ck g_ljkc x_ki^ac
!
!     Index restrictions: a, i, j are CCSD indices, 
!     c and k are CCSD + CC2 indices, and l is unrestricted 
!
!     Reorder g_ljkc as g_ljck and restrict indices (j, c, k) to active space
!
      call mem%alloc(g_ljck, wf%n_o, wf%n_ccsd_o, n_a_v, n_a_o)
!
!$omp parallel do private (l, j, k, c)
      do k = 1, n_a_o
         do c = 1, n_a_v
            do j = 1, wf%n_ccsd_o
               do l = 1, wf%n_o
!
                  g_ljck(l, j, c, k) = g_ljkc(l, j, k, c)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%alloc(X_ljai, wf%n_o, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
      call dgemm('N', 'N',                   &
                  wf%n_o*wf%n_ccsd_o,        &
                  wf%n_ccsd_v*wf%n_ccsd_o,   &
                  n_a_v*n_a_o,               &
                  one,                       &
                  g_ljck,                    &
                  wf%n_o*wf%n_ccsd_o,        &
                  x_ckai,                    &
                  n_a_v*n_a_o,               &
                  zero,                      &
                  X_ljai,                    &
                  wf%n_o*wf%n_ccsd_o)
!
      call mem%dealloc(g_ljck, wf%n_o, wf%n_ccsd_o, n_a_v, n_a_o)
!
!       X_kjbi = g_ljkc x_li^bc
!
!          Index restrictions: b, i, j are CCSD indices, 
!          c and l are CCSD + CC2 indices, and k is unrestricted 
!
!     Reorder g_ljkc as g_kjcl and restrict indices (j, c, l) to active space
!
      call mem%alloc(g_kjcl, wf%n_o, wf%n_ccsd_o, n_a_v, n_a_o)
!
!$omp parallel do private (l, j, k, c)
      do l = 1, n_a_o
         do c = 1, n_a_v
            do j = 1, wf%n_ccsd_o
               do k = 1, wf%n_o
!
                  g_kjcl(k, j, c, l) = g_ljkc(l, j, k, c)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%alloc(X_kjbi, wf%n_o, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
      call dgemm('N', 'N',                      &
                  (wf%n_o)*(wf%n_ccsd_o),       &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),  &
                  (n_a_o)*(n_a_v),              &
                  one,                          &
                  g_kjcl,                       &  ! g_kj_cl 
                  (wf%n_o)*(wf%n_ccsd_o),       &
                  x_ckai,                       &  ! x_cl_bi (x^cb_il)
                  (n_a_o)*(n_a_v),              &
                  zero,                         &
                  X_kjbi,                       & 
                  (wf%n_o)*(wf%n_ccsd_o))
!
      call mem%dealloc(g_kjcl, wf%n_o, wf%n_ccsd_o, n_a_v, n_a_o)
      call mem%dealloc(x_ckai, n_a_v, n_a_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
      wf%jacobian_c2_intermediate_oovo_2 = sequential_file('jacobian_c2_intermediate_oovo_2_ccsd')
      call wf%jacobian_c2_intermediate_oovo_2%open_('write', 'rewind')
      call wf%jacobian_c2_intermediate_oovo_2%write_(X_kjbi,(wf%n_o)*(wf%n_ccsd_v)*(wf%n_ccsd_o**2))
      call wf%jacobian_c2_intermediate_oovo_2%close_('keep')
!
      call mem%dealloc(X_kjbi, wf%n_o, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
!     X_ljai += - sum_kc L_ljkc x_ik^ac
!
!        Index restrictions: a, i, j are CCSD indices, 
!        c and k are CCSD + CC2 indices, and l is unrestricted 
!
!     Construct L_ljck = 2 g_ljkc - g_kjlc
!
      call mem%alloc(L_ljck, wf%n_o, wf%n_ccsd_o, n_a_v, n_a_o)
!
!$omp parallel do private (k, c, j, l)
      do k = 1, n_a_o
         do c = 1, n_a_v
            do j = 1, wf%n_ccsd_o
               do l = 1, wf%n_o
!
                  L_ljck(l, j, c, k) = two*g_ljkc(l, j, k, c) &
                                         - g_ljkc(k, j, l, c) 
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(g_ljkc, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
!     Restrict active indices (a, i) to the CCSD space
!
      call mem%alloc(x_ckai, n_a_v, n_a_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
!$omp parallel do private (i, a, k, c, ai, ck, ckai) collapse(2)
      do i = 1, wf%n_ccsd_o
         do a = 1, wf%n_ccsd_v
!
            ai = n_a_v*(i - 1) + a 
!
            do k = 1, n_a_o
               do c = 1, n_a_v
!
                  ck = n_a_v*(k - 1) + c
                  ckai = max(ck, ai)*(max(ck,ai)-3)/2 + ck + ai
!
                  x_ckai(c, k, a, i) = wf%x2(ckai)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call dgemm('N','N',                       &
                  (wf%n_o)*(wf%n_ccsd_o),       &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),  &
                  (n_a_o)*(n_a_v),              &
                  -one,                         &
                  L_ljck,                       & 
                  (wf%n_o)*(wf%n_ccsd_o),       &
                  x_ckai,                       & 
                  (n_a_o)*(n_a_v),              &
                  one,                          &
                  X_ljai,                       & 
                  (wf%n_o)*(wf%n_ccsd_o))
!
      call mem%dealloc(x_ckai, n_a_v, n_a_o, wf%n_ccsd_v, wf%n_ccsd_o)
      call mem%dealloc(L_ljck, wf%n_o, wf%n_ccsd_o, n_a_v, n_a_o)
!
      wf%jacobian_c2_intermediate_oovo_1 = sequential_file('jacobian_c2_intermediate_oovo_1_ccsd')
      call wf%jacobian_c2_intermediate_oovo_1%open_('write', 'rewind')
      call wf%jacobian_c2_intermediate_oovo_1%write_(X_ljai, (wf%n_o)*(wf%n_ccsd_v)*(wf%n_ccsd_o**2))
      call wf%jacobian_c2_intermediate_oovo_1%close_('keep')
!
      call mem%dealloc(X_ljai, wf%n_o, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
      call timer%turn_off()
!
   end subroutine save_jacobian_c2_intermediates_mlccsd
!
!
   module subroutine save_jacobian_d2_intermediate_mlccsd(wf)
!!
!!    Save jacobian d2 intermediate
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2019
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    Constructs the intermediate 
!!
!!       X_kibj = sum_dl g_kcbd x_ij^cd 
!!
!!       Index restrictions: b, i, j are CCSD indices,  
!!       c and d are CC2 + CCSD indices, and k is unrestricted
!!
!!    used in the d2-term. This is done only once in prepare_for_jacobian
!!    and the intermediate is stored in the file jacobian_d2_intermediate
!!    which is a wf variable.
!!
      implicit none
!
      class(mlccsd) :: wf
!
      type(timings), allocatable :: timer
!
      integer :: req1, req0
      integer :: current_b_batch
!
      type(batching_index) :: batch_b
!
      real(dp), dimension(:,:,:,:), allocatable :: g_bdkc, g_cdkb, X_ijkb, X_kibj, x_ijcd
!
      integer :: n_a_o, n_a_v
!
      integer :: i, k, c, ci, dj, cidj, d, j, b
!
      n_a_o = wf%n_ccsd_o + wf%n_cc2_o
      n_a_v = wf%n_ccsd_v + wf%n_cc2_v
!
      timer = timings('Jacobian CCSD D2 intermediate', pl='v')
      call timer%turn_on()
!
      call mem%alloc(X_kibj, wf%n_o, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
      call zero_array(X_kibj, wf%n_o*(wf%n_ccsd_o**2)*wf%n_ccsd_v)
!
!     Order amplitudes as x_ij_cd = x_ij^cd = x_ci_dj
!
      call mem%alloc(x_ijcd, wf%n_ccsd_o, wf%n_ccsd_o, n_a_v, n_a_v)
!
!$omp parallel do private (d, c, j, i, ci, dj, cidj)
      do d = 1, n_a_v
         do c = 1, n_a_v
            do j = 1, wf%n_ccsd_o
!
               dj = n_a_v*(j - 1) + d
!
               do i = 1, wf%n_ccsd_o
!
                  ci = n_a_v*(i - 1) + c
                  cidj = max(ci, dj)*(max(ci,dj)-3)/2 + ci + dj
!
                  x_ijcd(i, j, c, d) = wf%x2(cidj)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
!     Initialize batching variable
!
      req0 = n_a_v*(wf%n_o)*(wf%eri%n_J)
      req1 = n_a_v*(wf%eri%n_J) + 2*(wf%n_o)*(n_a_v)**2 + (wf%n_o)*(wf%n_ccsd_o**2)
!
      batch_b = batching_index(wf%n_ccsd_v)
      call mem%batch_setup(batch_b, req0, req1)
!
!     Start looping over b-batches
!
      do current_b_batch = 1, batch_b%num_batches
!
!        Get batching limits for current b-batch
!
         call batch_b%determine_limits(current_b_batch)
!
!        Form g_kc_db = g_kcbd
!
         call mem%alloc(g_bdkc, batch_b%length, n_a_v, wf%n_o, n_a_v)
!
         call wf%eri%get_eri_t1('vvov', g_bdkc, batch_b%first, batch_b%last, 1, n_a_v, &
                                                1, wf%n_o, 1, n_a_v)
!
!        Reorder g_bd_kc to g_cd_kb (= g_kcbd), i.e. 1234 to 4231
!
         call mem%alloc(g_cdkb, n_a_v, n_a_v, wf%n_o, batch_b%length)
!
         call sort_1234_to_4231(g_bdkc, g_cdkb, batch_b%length, n_a_v, wf%n_o, n_a_v)
!
         call mem%dealloc(g_bdkc, batch_b%length, n_a_v, wf%n_o, n_a_v)
!
!        Form intermediate X_ij_kb = sum_cd g_kcdb t_ij^cd
!                                  = sum_cd t_ij_cd g_cd_kb
!
         call mem%alloc(X_ijkb, wf%n_ccsd_o, wf%n_ccsd_o, wf%n_o, batch_b%length)
!
         call dgemm('N', 'N',                   &
                     (wf%n_ccsd_o)**2,          &
                     (wf%n_o)*(batch_b%length), &
                     (n_a_v)**2,                &
                     one,                       &
                     x_ijcd,                    & ! x_ij_cd
                     (wf%n_ccsd_o)**2,          &
                     g_cdkb,                    & ! g_cd_kb
                     (n_a_v)**2,                &
                     zero,                      &
                     X_ijkb,                    & ! X_ij_kb
                     (wf%n_ccsd_o)**2)
!
         call mem%dealloc(g_cdkb, n_a_v, n_a_v, wf%n_o, batch_b%length)
!
!        Reorder and add to full space X_kibj = X_ij_kb
!
!$omp parallel do private(j, b, i, k)
         do j = 1, wf%n_ccsd_o
            do b = 1, batch_b%length
               do i = 1, wf%n_ccsd_o
                  do k = 1, wf%n_o
!
                     X_kibj(k, i, b + batch_b%first - 1, j) = X_ijkb(i, j, k, b)
!
                  enddo
               enddo
            enddo
         enddo
!$omp end parallel do
!
         call mem%dealloc(X_ijkb, wf%n_ccsd_o, wf%n_ccsd_o, wf%n_o, batch_b%length)
!
      enddo
!
!     Store intermediate to file 
!
      wf%jacobian_d2_intermediate = sequential_file('jacobian_d2_intermediate_ccsd')
      call wf%jacobian_d2_intermediate%open_('write', 'rewind')
      call wf%jacobian_d2_intermediate%write_(X_kibj, wf%n_o*(wf%n_ccsd_o**2)*wf%n_ccsd_v)
      call wf%jacobian_d2_intermediate%close_('keep')
!
      call mem%dealloc(X_kibj, wf%n_o, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
      call mem%dealloc(x_ijcd, wf%n_ccsd_o, wf%n_ccsd_o, n_a_v, n_a_v)
!
      call timer%turn_off()
!
   end subroutine save_jacobian_d2_intermediate_mlccsd
!
!
   module subroutine save_jacobian_e2_intermediate_mlccsd(wf, L_ckdl)
!!
!!    Save jacobian e2 intermediate
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2019
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    Constructs the intermediate 
!!
!!       X_ckbj = sum_dl L_kcld x_jl^bd 
!!
!!       Index restrictions: b, j are CCSD indices,  
!!       l, d, k and c are CC2 + CCSD indices.
!!
!!    used in the e2-term. This is done only once in prepare_for_jacobian
!!    and the intermediate is stored in the file jacobian_e2_intermediate
!!    which is a wf variable.
!!
      implicit none
!
      class(mlccsd) :: wf
!
      real(dp), dimension(wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o, &
               wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o), intent(in) :: L_ckdl
!
      real(dp), dimension(:,:,:,:), allocatable :: x_dlbj
      real(dp), dimension(:,:,:,:), allocatable :: x_ckbj
!
      type(timings), allocatable :: timer
!
      integer :: n_a_o, n_a_v
!
      integer :: l, d, j, b, bj, dl, bjdl
!
      timer = timings('Jacobian CCSD E2 intermediate', pl='v')
      call timer%turn_on()
!
      n_a_o = wf%n_ccsd_o + wf%n_cc2_o
      n_a_v = wf%n_ccsd_v + wf%n_cc2_v
!
!     X_ckbj = sum_dl L_kcld x_jl^bd 
!
!     Index restrictions: b, j are CCSD indices,  
!     l, d, k and c are CC2 + CCSD indices.
!
      call mem%alloc(x_dlbj, n_a_v, n_a_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
!$omp parallel do private(j, b, l, d, bj, dl, bjdl) collapse(2)
      do j = 1, wf%n_ccsd_o
         do b = 1, wf%n_ccsd_v
!
            bj = n_a_v*(j - 1) + b
!
            do l = 1, n_a_o
               do d = 1, n_a_v
!
                  dl = n_a_v*(l - 1) + d
                  bjdl = max(bj,dl)*(max(bj,dl) - 3)/2 + bj + dl
!
                  x_dlbj(d,l,b,j) = wf%x2(bjdl)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%alloc(X_ckbj, n_a_v, n_a_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
     call dgemm('N', 'N',                       &
                 (n_a_o)*(n_a_v),               &
                 (wf%n_ccsd_o)*(wf%n_ccsd_v),   &
                 (n_a_o)*(n_a_v),               &
                 one,                           &
                 L_ckdl,                        & 
                 (n_a_o)*(n_a_v),               &
                 x_dlbj,                        & 
                 (n_a_o)*(n_a_v),               &
                 zero,                          &
                 X_ckbj,                        & 
                 (n_a_o)*(n_a_v))
!
      call mem%dealloc(x_dlbj, n_a_v, n_a_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
!     Store intermediate to file 
!
      wf%jacobian_e2_intermediate = sequential_file('jacobian_e2_intermediate_ccsd')
      call wf%jacobian_e2_intermediate%open_('write', 'rewind')
      call wf%jacobian_e2_intermediate%write_(X_ckbj, n_a_o*n_a_v*(wf%n_ccsd_o)*(wf%n_ccsd_v))
      call wf%jacobian_e2_intermediate%close_('keep')
!
      call mem%dealloc(X_ckbj, n_a_v, n_a_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
      call timer%turn_off()
!
   end subroutine save_jacobian_e2_intermediate_mlccsd
!
!
   module subroutine save_jacobian_g2_intermediates_mlccsd(wf, L_ckdl)
!!
!!    Save jacobian g2 intermediates
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2019
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    Constructs the intermediate 
!!
!!       X_ckbj = sum_dl L_kcld x_lj^bd 
!!
!!       X_db = sum_ckl L_kcld x_lk^bc 
!!
!!       X_lj = sum_cdl L_kcld x_kj^cd
!!
!!       Index restrictions: b, j are CCSD indices,  
!!       l, d, k and c are CC2 + CCSD indices.
!!
!!    used in the g2-term. This is done only once in prepare_for_jacobian
!!    and the intermediate is stored in the files jacobian_G2_intermediate_vovo
!!    jacobian_G2_intermediate_vv, and jacobian_G2_intermediate_oo which
!!    are wf variables.
!!
      implicit none
!
      class(mlccsd) :: wf
!
      real(dp), dimension(wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o, &
               wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o), intent(in) :: L_ckdl
!
      real(dp), dimension(:,:,:,:), allocatable :: x_blck, x_dlbj, X_ckbj
      real(dp), dimension(:,:), allocatable :: X_db, X_lj
!
      type(timings), allocatable :: timer
!
      integer :: j, b, l, d, bl, dj, bldj
      integer :: n_a_o, n_a_v
!
      timer = timings('Jacobian CCSD G2 intermediate', pl='v')
      call timer%turn_on()
!
      n_a_o = wf%n_ccsd_o + wf%n_cc2_o
      n_a_v = wf%n_ccsd_v + wf%n_cc2_v
!
!     X_ckbj = sum_dl L_kcld x_lj^bd 
!
!     Index restrictions: b, j are CCSD indices,  
!     l, d, k and c are CC2 + CCSD indices.
!
!     x_bldj ordered as x_dlbj
!
      call mem%alloc(x_dlbj, n_a_v, n_a_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
!$omp parallel do private(l,d,j,b, dj, bl, bldj) collapse(2)
      do j = 1, wf%n_ccsd_o
         do b = 1, wf%n_ccsd_v
            do l = 1, n_a_o
!
               bl = n_a_v*(l - 1) + b
!
               do d = 1, n_a_v

                  dj = n_a_v*(j - 1) + d
                  bldj = max(bl,dj)*(max(bl,dj)-3)/2 + bl + dj
!
                  x_dlbj(d, l, b, j) = wf%x2(bldj)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%alloc(X_ckbj, n_a_v, n_a_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
!     X_ckbj = sum_dl x_bldj * L_kcld 
!
      call dgemm('N', 'N',                      &
                  (n_a_o)*(n_a_v),              &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),  &
                  (n_a_o)*(n_a_v),              &
                  one,                          &
                  L_ckdl,                       & 
                  (n_a_o)*(n_a_v),              &
                  x_dlbj,                       &
                  (n_a_o)*(n_a_v),              &
                  zero,                         &
                  X_ckbj,                       & 
                  (n_a_o)*(n_a_v))
!
      call mem%dealloc(x_dlbj, n_a_v, n_a_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
!     Store intermediate to file 
!
      wf%jacobian_g2_intermediate_vovo = sequential_file('jacobian_g2_intermediate_vovo_ccsd')
      call wf%jacobian_g2_intermediate_vovo%open_('write', 'rewind')
      call wf%jacobian_g2_intermediate_vovo%write_(X_ckbj, n_a_o*n_a_v*(wf%n_ccsd_o)*(wf%n_ccsd_v))
      call wf%jacobian_g2_intermediate_vovo%close_('keep')
!
      call mem%dealloc(X_ckbj, n_a_v, n_a_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
!     X_db = sum_dl L_kcld x_lk^bc 
!
!     Index restrictions: b, j are CCSD indices,  
!     l, d, k and c are CC2 + CCSD indices.
!
      call mem%alloc(x_blck, n_a_v, n_a_o, n_a_v, n_a_o)
      call squareup(wf%x2, x_blck, (n_a_v)*(n_a_o))
!
      call mem%alloc(X_db, n_a_v, wf%n_ccsd_v)
!
      call dgemm('N', 'T',            &
                  n_a_v,              &
                  wf%n_ccsd_v,        &
                  (n_a_o**2)*(n_a_v), &
                  one,                &
                  L_ckdl,             & ! L_d_lck
                  n_a_v,              &
                  x_blck,             & ! x_b_lck
                  n_a_v,              &
                  zero,               &
                  X_db,               &
                  n_a_v)
!
!     Store intermediate to file 
!
      wf%jacobian_g2_intermediate_vv = sequential_file('jacobian_g2_intermediate_vv_ccsd')
      call wf%jacobian_g2_intermediate_vv%open_('write', 'rewind')
      call wf%jacobian_g2_intermediate_vv%write_(X_db, n_a_v*(wf%n_ccsd_v))
      call wf%jacobian_g2_intermediate_vv%close_('keep')
!
      call mem%dealloc(X_db, n_a_v, wf%n_ccsd_v)
!
!     X_lj = sum_cdl L_kcld x_kj^ck
!
!     Index restrictions: b, j are CCSD indices,  
!     l, d, k and c are CC2 + CCSD indices.
!
!     NOTE: treat x_blck as x_ckdj 
!
      call mem%alloc(X_lj, n_a_o, wf%n_ccsd_o)
!
!     X_lj = sum_ckd L_lckd * x_ckdj
!
      call dgemm('T', 'N',             &
                  n_a_o,               &
                  wf%n_ccsd_o,         &
                  (n_a_v**2)*(n_a_o),  &
                  one,                 &
                  L_ckdl,              & 
                  (n_a_v**2)*(n_a_o),  &
                  x_blck,              & ! x_ckdj
                  (n_a_v**2)*(n_a_o),  &
                  zero,                &
                  X_lj,                & ! X_l_j
                  n_a_o)
!
      call mem%dealloc(x_blck, n_a_v, n_a_o, n_a_v, n_a_o)
!
!     Store intermediate to file 
!
      wf%jacobian_g2_intermediate_oo = sequential_file('jacobian_g2_intermediate_oo_ccsd')
      call wf%jacobian_g2_intermediate_oo%open_('write', 'rewind')
      call wf%jacobian_g2_intermediate_oo%write_(X_lj, n_a_o*(wf%n_ccsd_o))
      call wf%jacobian_g2_intermediate_oo%close_('keep')
!
      call mem%dealloc(X_lj, n_a_o, wf%n_ccsd_o)
!
      call timer%turn_off()
!
   end subroutine save_jacobian_g2_intermediates_mlccsd
!
!
   module subroutine save_jacobian_h2_intermediates_mlccsd(wf, g_ckdl)
!!
!!    Save jacobian h2 intermediates
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2019
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    Constructs the intermediate 
!!
!!       X_aidl = x^ca_ik g_kcld
!!
!!       X_ajdk = x^ca_jl g_kcld
!!
!!       Index restrictions: a and i are CCSD indices,  
!!       l, d, k and c are CC2 + CCSD indices.
!!
!!    used in the h2-term. This is done only once in prepare_for_jacobian
!!    and the intermediate is stored in the file jacobian_h2_intermediate_ovov_1
!!    and jacobian_h2_intermediate_ovov_2 which are wf variables.
!!
      implicit none
!
      class(mlccsd) :: wf
!
      real(dp), dimension(wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o, &
               wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_o + wf%n_cc2_o), intent(in) :: g_ckdl
!
      real(dp), dimension(:,:,:,:), allocatable :: x_aick, X_aidl, x_ajcl, g_cldk, X_ajdk
!
      type(timings), allocatable :: timer
!
      integer :: n_a_o, n_a_v
      integer :: k, c, i, a, ci, ak, akci, j, cj, al, alcj, l
!
      timer = timings('Jacobian CCSD H2 intermediate', pl='v')
      call timer%turn_on()
!
      n_a_o = wf%n_ccsd_o + wf%n_cc2_o
      n_a_v = wf%n_ccsd_v + wf%n_cc2_v
!
!     X_aidl = x^ca_ik g_kcld
!
!     Index restrictions: a and i are CCSD indices,  
!     l, d, k and c are CC2 + CCSD indices.
!
!     x_ak,ci ordered as x_aikc
!
      call mem%alloc(x_aick, wf%n_ccsd_v, wf%n_ccsd_o, n_a_v, n_a_o)
!
!$omp parallel do private (c, k, i, a, ci, ak, akci) collapse(2)
      do k = 1, n_a_o
         do c = 1, n_a_v
            do i = 1, wf%n_ccsd_o
!
               ci = n_a_v*(i - 1) + c
!
               do a = 1, wf%n_ccsd_v
!
                  ak = n_a_v*(k - 1) + a
                  akci = max(ak,ci)*(max(ak,ci)-3)/2 + ak + ci
!
                  x_aick(a, i, c, k) = wf%x2(akci)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
! 
      call mem%alloc(X_aidl, wf%n_ccsd_v, wf%n_ccsd_o, n_a_v, n_a_o)
!
      call dgemm('N', 'N',                      &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),  &
                  (n_a_o)*(n_a_v),              &
                  (n_a_o)*(n_a_v),              &
                  one,                          &
                  x_aick,                       & 
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),  &
                  g_ckdl,                       & 
                  (n_a_o)*(n_a_v),              &
                  zero,                         &
                  X_aidl,                       & 
                  (wf%n_ccsd_o)*(wf%n_ccsd_v))
!
      call mem%dealloc(x_aick, wf%n_ccsd_v, wf%n_ccsd_o, n_a_v, n_a_o)
!
!     Store intermediate to file 
!
      wf%jacobian_h2_intermediate_vovo_1 = sequential_file('jacobian_h2_intermediate_vovo_1_ccsd')
      call wf%jacobian_h2_intermediate_vovo_1%open_('write', 'rewind')
      call wf%jacobian_h2_intermediate_vovo_1%write_(X_aidl, &
                  n_a_v*n_a_o*(wf%n_ccsd_o)*(wf%n_ccsd_v))
      call wf%jacobian_h2_intermediate_vovo_1%close_('keep')
!
      call mem%dealloc(X_aidl, wf%n_ccsd_v, wf%n_ccsd_o, n_a_v, n_a_o)
!
!     X_ajdk = x^ca_jl g_kcld
!
!     Index restrictions: a i, and j are CCSD indices,  
!     l, d, k and c are CC2 + CCSD indices.
!
!     x_alcj ordered as x_ajlc
!
      call mem%alloc(x_ajcl, wf%n_ccsd_v, wf%n_ccsd_o, n_a_v, n_a_o)
!
!$omp parallel do private (c, l, j, a, cj, al, alcj) collapse(2)
      do c = 1, n_a_v
         do l = 1, n_a_o
            do j = 1, wf%n_ccsd_o
!
               cj = n_a_v*(j - 1) + c
!
               do a = 1, wf%n_ccsd_v
!
                  al = n_a_v*(l - 1) + a
                  alcj = max(al,cj)*(max(al,cj)-3)/2 + al + cj
!
                  x_ajcl(a, j, c, l) = wf%x2(alcj)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%alloc(g_cldk, n_a_v, n_a_o, n_a_v, n_a_o)
      call sort_1234_to_1432(g_ckdl, g_cldk, n_a_v, n_a_o, n_a_v, n_a_o)
!
      call mem%alloc(X_ajdk, wf%n_ccsd_v, wf%n_ccsd_o, n_a_v, n_a_o)
!
!     X_ajkd = sum_lc t_ajlc * g_lckd
!
      call dgemm('N', 'N',                      &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),  &
                  (n_a_o)*(n_a_v),              &
                  (n_a_o)*(n_a_v),              &
                  one,                          &
                  x_ajcl,                       &
                  (wf%n_ccsd_o)*(wf%n_ccsd_v),  &
                  g_cldk,                       & ! g_lckd
                  (n_a_o)*(n_a_v),              &
                  zero,                         &
                  X_ajdk,                       & ! X_aj_kd
                  (wf%n_ccsd_o)*(wf%n_ccsd_v))
!
      call mem%dealloc(g_cldk, n_a_v, n_a_o, n_a_v, n_a_o)
      call mem%dealloc(x_ajcl, wf%n_ccsd_v, wf%n_ccsd_o, n_a_v, n_a_o)
!
!     Store intermediate to file 
!
      wf%jacobian_h2_intermediate_vovo_2 = sequential_file('jacobian_h2_intermediate_vovo_2_ccsd')
      call wf%jacobian_h2_intermediate_vovo_2%open_('write', 'rewind')
      call wf%jacobian_h2_intermediate_vovo_2%write_(X_ajdk, &
               n_a_v*n_a_o*(wf%n_ccsd_o)*(wf%n_ccsd_v))
      call wf%jacobian_h2_intermediate_vovo_2%close_('keep')
!
      call mem%dealloc(X_ajdk, wf%n_ccsd_v, wf%n_ccsd_o, n_a_v, n_a_o)
!
      call timer%turn_off()
!
   end subroutine save_jacobian_h2_intermediates_mlccsd
!
!
   module subroutine save_jacobian_j2_intermediates_mlccsd(wf, g_klcd)
!!
!!    Save jacobian j2 intermediates
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2019
!!
!!    Adapted for MLCCSD by Sarai D. Folkestad, 2019
!!
!!    Constructs the intermediate 
!!
!!       X_klij = g_kcld x^cd_ij
!!
!!       Index restrictions: i and j are CCSD indices,  
!!       l, d, k and c are CC2 + CCSD indices.
!!
!!    Also saves the integral g_kcld ordered as g_klcd
!!
!!    used in the j2-term. This is done only once in prepare_for_jacobian
!!    and the intermediate is stored in the file jacobian_j2_intermediate_oooo
!!    and jacobian_ji_intermediate_oovv which are wf variables.
!!
      implicit none
!
      class(mlccsd) :: wf
!
      real(dp), dimension(wf%n_ccsd_o + wf%n_cc2_o, wf%n_ccsd_o + wf%n_cc2_o, &
               wf%n_ccsd_v + wf%n_cc2_v, wf%n_ccsd_v + wf%n_cc2_v), intent(in) :: g_klcd
!  
      real(dp), dimension(:,:,:,:), allocatable :: x_cdij, X_klij
!
      type(timings), allocatable :: timer
!
      integer :: n_a_o, n_a_v
      integer :: j, i, d, c, dj, ci, cidj
!
!     X_klij = g_kcld x^cd_ij
!
!     Index restrictions: i and j are CCSD indices,  
!     l, d, k and c are CC2 + CCSD indices.
!
      timer = timings('Jacobian MLCCSD J2 intermediate', pl='v')
      call timer%turn_on()
!
      n_a_o = wf%n_ccsd_o + wf%n_cc2_o
      n_a_v = wf%n_ccsd_v + wf%n_cc2_v
!
!     Reorder amplitudes x_ij^cd as x_ijcd
!
      call mem%alloc(x_cdij, n_a_v, n_a_v, wf%n_ccsd_o, wf%n_ccsd_o)
!
!$omp parallel do private(i, j, d, c, ci, dj, cidj) collapse(2)
      do j = 1, wf%n_ccsd_o
         do i = 1, wf%n_ccsd_o
            do d = 1, n_a_v
!
               dj = n_a_v*(j - 1) + d
!
               do c = 1, n_a_v
!
                  ci = n_a_v*(i - 1) + c
                  cidj = max(ci,dj)*(max(ci,dj)-3)/2 + ci + dj
!
                  x_cdij(c, d, i, j) = wf%x2(cidj)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%alloc(X_klij, n_a_o, n_a_o, wf%n_ccsd_o, wf%n_ccsd_o)
!
      call dgemm('N', 'N',          &
                  (n_a_o)**2,       &
                  (wf%n_ccsd_o)**2, &
                  (n_a_v)**2,       &
                  one,              &
                  g_klcd,           & 
                  (n_a_o)**2,       &
                  x_cdij,           &
                  (n_a_v)**2,       &
                  zero,             &
                  X_klij,           & 
                  (n_a_o)**2)
!
      call mem%dealloc(x_cdij, n_a_v, n_a_v, wf%n_ccsd_o, wf%n_ccsd_o)
!
!     Store intermediate to file 
!
      wf%jacobian_j2_intermediate_oooo = sequential_file('jacobian_j2_intermediate_oooo_ccsd')
      call wf%jacobian_j2_intermediate_oooo%open_('write', 'rewind')
      call wf%jacobian_j2_intermediate_oooo%write_(X_klij, n_a_o**2 * (wf%n_ccsd_o**2))
      call wf%jacobian_j2_intermediate_oooo%close_('keep')
!
      call mem%dealloc(X_klij, n_a_o, n_a_o, wf%n_ccsd_o, wf%n_ccsd_o)
!
!     Storing the integrals g_kcld ordered as g_klcd
!
      wf%jacobian_j2_intermediate_oovv = sequential_file('jacobian_j2_intermediate_oovv_ccsd')
      call wf%jacobian_j2_intermediate_oovv%open_('write', 'rewind')
      call wf%jacobian_j2_intermediate_oovv%write_(g_klcd, n_a_v**2 * n_a_o**2)
      call wf%jacobian_j2_intermediate_oovv%close_('keep')
!
      call timer%turn_off()
!
   end subroutine save_jacobian_j2_intermediates_mlccsd
!
!
end submodule jacobian_mlccsd
