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
submodule (ccs_class) F_ccs
!
!!
!!    F submodule (CCS)
!!
!!    Routines for the linear transform of
!!    vectors by the F matrix
!!
!!    ρ = F * c,
!!
!!    where
!!
!!    F_μ,ν = < Λ' | [[ exp(T) H_0 exp(T), τ_μ ], τ_ν ] | R >,
!!
!!    Where < Λ' | = < R | + sum_μ tbar_μ < μ |
!!
!!    F transformation routines for CCS written by Sarai D. Folkestad and
!!    Eirik F. Kjønstad, Feb 2019, and debugged by Josefine H. Andersen,
!!    spring 2019.
!!
!
   implicit none
!
!
contains
!
!
   module subroutine F_transformation_ccs(wf, c, rho)
!!
!!    F transformation
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2019
!!
!!    Directs the transformation by the F matrix.
!!
!!     F_mu,nu = < tbar | [[H-bar,tau_mu],tau_nu] | HF >
!!
!!    Modified for ccsd F transformation by
!!    (A. K. Schnack-Petersen and) Eirik F. Kjønstad Sep 2021
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: c
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: rho
!
      real(dp), dimension(:), allocatable :: tbar
!
      call mem%alloc(tbar, wf%n_es_amplitudes)
      call wf%get_full_multipliers(tbar)
!
      call wf%F_x_transformation(c, rho, tbar)
!
      call mem%dealloc(tbar, wf%n_es_amplitudes)
!
   end subroutine F_transformation_ccs
!
!
   module subroutine F_x_transformation_ccs(wf, c, rho, x)
!!
!!    F(x) transformation
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2019
!!
!!    Directs the transformation by the F(X) matrix, defined by
!!
!!     F(X)_mu,nu = < X | [[H-bar,tau_mu],tau_nu] | HF >
!!
!!    where < X | = < HF | + < mu | X_mu.
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: c
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: rho
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: x
!
      call wf%F_x_mu_transformation(c, rho, x)
!
      call wf%F_ccs_a1_0(c, rho)
!
   end subroutine F_x_transformation_ccs
!
!
   module subroutine F_x_mu_transformation_ccs(wf, c, rho, x)
!!
!!    F(X) mu transformation
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2019
!!
!!    Directs the transformation by the excited configuration contribution (mu)
!!    of the F matrix:
!!
!!       F'(X)_mu,nu = < X | [[H-bar,tau_mu],tau_nu] | HF >
!!
!!    with < X / = < mu / X_mu. For the full F transformation, see F_transformation
!!    routine.
!!
!!    Edited by A. K. Schnack-Petersen og E. F. Kjønstad, Sep 2021
!!
      use array_initialization, only: zero_array
!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(wf%n_t1), intent(in) :: c
      real(dp), dimension(wf%n_t1), intent(out) :: rho
!
      real(dp), dimension(wf%n_t1), intent(in) :: x
!
      call zero_array(rho, wf%n_t1)
!
      call wf%F_ccs_a1_1(c, rho, x)
      call wf%F_ccs_b1_1(c, rho, x)
      call wf%F_ccs_c1_1(c, rho, x)
!
   end subroutine F_x_mu_transformation_ccs
!
!
   module subroutine F_ccs_a1_0_ccs(wf, c_ai, rho_ai)
!!
!!    F transformation A1,0 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2019
!!
!!    rho_A1,0_ai = 2 * L_iajb * c_bj
!!
      use reordering, only: add_2143_to_1234, add_2341_to_1234
!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
!
!     Local variables
!
      real(dp), dimension(:,:,:,:), allocatable :: g_iajb
      real(dp), dimension(:,:,:,:), allocatable :: L_aibj
!
      call mem%alloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call wf%eri_t1%get('ovov', g_iajb)
!
!     L_iajb = 2 g_iajb - g_ibja (ordered as L_aibj)
!
      call mem%alloc(L_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o, set_zero=.true.)
      call add_2143_to_1234(two, g_iajb, L_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call add_2341_to_1234(-one, g_iajb, L_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     rho_ai += 2 * L_iajb c_bj
!
      call dgemm('N', 'N',             &
                  (wf%n_v)*(wf%n_o),   &
                  1,                   &
                  (wf%n_v)*(wf%n_o),   &
                  two,                 &
                  L_aibj,              &
                  (wf%n_v)*(wf%n_o),   &
                  c_ai,                & ! c_bj
                  (wf%n_v)*(wf%n_o),   &
                  one,                 &
                  rho_ai,              &
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(L_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine F_ccs_a1_0_ccs
!
!
   module subroutine F_ccs_a1_1_ccs(wf, c_ai, rho_ai, tbar_ai)
!!
!!    F transformation A1,1 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2019
!!
!!    rho_A1,1_ai = - (F_ib * tbar_aj + F_ja * tbar_bi) * c_bj
!!
!!    Modified for ccsd F transformation by
!!    (A. K. Schnack-Petersen and) Eirik F. Kjønstad Sep 2021
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: tbar_ai
!
!     Local variables
!
      real(dp), dimension(:,:), allocatable :: X_ij
      real(dp), dimension(:,:), allocatable :: X_ji
!
!     Term 1 : - F_ib * tbar_aj * c_bj
!
!     X_ij = F_ib * c_bj
!
      call mem%alloc(X_ij, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N',    &
                  wf%n_o,     &
                  wf%n_o,     &
                  wf%n_v,     &
                  one,        &
                  wf%fock_ia, & ! F_i,b
                  wf%n_o,     &
                  c_ai,       & ! c_b,j
                  wf%n_v,     &
                  zero,       &
                  X_ij,       &
                  wf%n_o)
!
!     rho_ai -= tbar_aj * X_ij
!
      call dgemm('N', 'T',    &
                  wf%n_v,     &
                  wf%n_o,     &
                  wf%n_o,     &
                  -one,       &
                  tbar_ai,    & ! tbar_a,j
                  wf%n_v,     &
                  X_ij,       & ! X_j,i
                  wf%n_o,     &
                  one,        &
                  rho_ai,     &
                  wf%n_v)
!
      call mem%dealloc(X_ij, wf%n_o, wf%n_o)
!
!     Term 2 : - F_ja * tbar_bi * c_bj
!
!     X_ji = c_bj * tbar_bi
!
      call mem%alloc(X_ji, wf%n_o, wf%n_o)
!
      call dgemm('T', 'N',    &
                  wf%n_o,     &
                  wf%n_o,     &
                  wf%n_v,     &
                  one,        &
                  c_ai,       & ! c_bj
                  wf%n_v,     &
                  tbar_ai,    & ! tbar_b,i
                  wf%n_v,     &
                  zero,       &
                  X_ji,       &
                  wf%n_o)
!
!     rho_ai -= F_ja * X_ji
!
      call dgemm('T', 'N',    &
                  wf%n_v,     &
                  wf%n_o,     &
                  wf%n_o,     &
                  -one,       &
                  wf%fock_ia, & ! F_j,a
                  wf%n_o,     &
                  X_ji,       &
                  wf%n_o,     &
                  one,        &
                  rho_ai,     &
                  wf%n_v)
!
      call mem%dealloc(X_ji, wf%n_o, wf%n_o)
!
   end subroutine F_ccs_a1_1_ccs
!
!
   module subroutine F_ccs_b1_1_ccs(wf, c_ai, rho_ai, tbar_ai)
!!
!!    F transformation B1,1 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2019
!!
!!    rho_B1,1_ai = - (L_ikjb * tbar_ak + L_jkia * tbar_bk) * c_bj
!!
!!    Modified for ccsd F transformation by
!!    (A. K. Schnack-Petersen and) Eirik F. Kjønstad Sep 2021
!!
      use reordering, only: add_2143_to_1234, add_4123_to_1234
!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: tbar_ai
!
!     Local variables
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ikjb
      real(dp), dimension(:,:,:,:), allocatable :: L_kibj
!
      real(dp), dimension(:,:), allocatable :: X_ki
      real(dp), dimension(:,:), allocatable :: X_kj
!
!     Term 1: - L_ikjb * tbar_ak * c_bj
!
!     L_ikjb = 2 g_ikjb - g_jkib (ordered as L_kibj)
!
      call mem%alloc(g_ikjb, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call wf%eri_t1%get('ooov', g_ikjb)
!
      call mem%alloc(L_kibj, wf%n_o, wf%n_o, wf%n_v, wf%n_o, set_zero=.true.)
      call add_2143_to_1234(two, g_ikjb, L_kibj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call add_4123_to_1234(-one, g_ikjb, L_kibj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_ikjb, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
!     X_ki = L_kibj * c_bj
!
      call mem%alloc(X_ki, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N',             &
                  wf%n_o**2,           &
                  1,                   &
                  (wf%n_o)*(wf%n_v),   &
                  one,                 &
                  L_kibj,              & ! L_ki_jb
                  wf%n_o**2,           &
                  c_ai,                & ! c_bj
                  (wf%n_o)*(wf%n_v),   &
                  zero,                &
                  X_ki,                &
                  wf%n_o**2)
!
!      rho_ai += tbar_ak * X_ki
!
      call dgemm('N', 'N',             &
                  wf%n_v,              &
                  wf%n_o,              &
                  wf%n_o,              &
                  -one,                &
                  tbar_ai,             & ! tbar_ak
                  wf%n_v,              &
                  X_ki,                &
                  wf%n_o,              &
                  one,                 &
                  rho_ai,              &
                  wf%n_v)
!
      call mem%dealloc(X_ki, wf%n_o, wf%n_o)
!
!     Term 2: - L_jkia * tbar_bk * c_bj
!
!     X_kj = tbar_bk c_bj
!
      call mem%alloc(X_kj, wf%n_o, wf%n_o)
!
      call dgemm('T', 'N',             &
                  wf%n_o,              &
                  wf%n_o,              &
                  wf%n_v,              &
                  one,                 &
                  tbar_ai,             & ! tbar_bk
                  wf%n_v,              &
                  c_ai,                & ! c_bj
                  wf%n_v,              &
                  zero,                &
                  X_kj,                &
                  wf%n_o)
!
!     rho_ai -= L_kjai * X_kj
!
      call dgemm('T', 'N',             &
                  (wf%n_v)*(wf%n_o),   &
                  1,                   &
                  (wf%n_o)**2,         &
                  -one,                &
                  L_kibj,              & ! L_kj_ai
                  wf%n_o**2,           &
                  X_kj,                &
                  wf%n_o**2,           &
                  one,                 &
                  rho_ai,              &
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(X_kj, wf%n_o, wf%n_o)
      call mem%dealloc(L_kibj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine F_ccs_b1_1_ccs
!
!
   module subroutine F_ccs_c1_1_ccs(wf, c_ai, rho_ai, tbar_ai)
!!
!!    F transformation C1,1 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2019
!!
!!    rho_C1,1_ai = (L_cajb * tbar_ci + L_cbia * tbar_cj) * c_bj
!!                = (X_iajb + X_jbia) * c_bj
!!
!!    Modified for ccsd F transformation by
!!    (A. K. Schnack-Petersen and) Eirik F. Kjønstad Sep 2021
!!
      use batching_index_class, only: batching_index
      use reordering, only: add_1432_to_1234, sort_1234_to_2143
!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_v,wf%n_o), intent(in)     :: tbar_ai
!
!     Local variables
!
      real(dp), dimension(:,:,:,:), allocatable :: X_iajb
      real(dp), dimension(:,:,:,:), allocatable :: X_aibj
      real(dp), dimension(:,:,:,:), allocatable :: g_cajb
      real(dp), dimension(:,:,:,:), allocatable :: L_cajb
!
      type(batching_index) :: batch_c
!
      integer :: current_c_batch, req0, req1
!
      call mem%alloc(X_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v, set_zero=.true.)
!
!     Prepare for batching over index c
!
      req0 = (wf%eri_t1%n_J)*(wf%n_o)*(wf%n_v)
!
      req1 = max(2*(wf%n_o)*(wf%n_v**2), (wf%n_o)*(wf%n_v**2) + (wf%eri_t1%n_J)*(wf%n_v))
!
      batch_c = batching_index(wf%n_v)
!
      call mem%batch_setup(batch_c, req0, req1, tag='F_ccs_c1_1_ccs')
!
      do current_c_batch = 1, batch_c%num_batches
!
!        Determine the limits for the current c-batch
!
         call batch_c%determine_limits(current_c_batch)
!
!        L_cajb = 2 g_cajb - g_cbja
!
         call mem%alloc(g_cajb, batch_c%length, wf%n_v, wf%n_o, wf%n_v)
!
         call wf%eri_t1%get('vvov', g_cajb, first_p = batch_c%first, &
                                                last_p = batch_c%get_last())
!
         call mem%alloc(L_cajb, batch_c%length, wf%n_v, wf%n_o, wf%n_v)
!
         call dcopy(batch_c%length*wf%n_v**2*(wf%n_o), g_cajb, 1, L_cajb, 1)
         call dscal(batch_c%length*wf%n_v**2*(wf%n_o), two, L_cajb, 1)
!
         call add_1432_to_1234(-one, g_cajb, L_cajb, batch_c%length, wf%n_v, wf%n_o, wf%n_v)
!
         call mem%dealloc(g_cajb, batch_c%length, wf%n_v, wf%n_o, wf%n_v)
!
!        X_iajb = tbar_ci * L_cajb
!
         call dgemm('T', 'N',                   &
                     wf%n_o,                    &
                     (wf%n_v**2)*wf%n_o,        &
                     batch_c%length,            &
                     one,                       &
                     tbar_ai(batch_c%first, 1), & ! tbar_c_i
                     wf%n_v,                    &
                     L_cajb,                    & ! L_c_ajb
                     batch_c%length,            &
                     one,                       &
                     X_iajb,                    & ! X_i_ajb
                     wf%n_o)
!
         call mem%dealloc(L_cajb, batch_c%length, wf%n_v, wf%n_o, wf%n_v)
!
      enddo ! batches of c
!
      call mem%batch_finalize()
!
      call mem%alloc(X_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call sort_1234_to_2143(X_iajb, X_aibj, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%dealloc(X_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     rho_ai += X_aibj * c_bj
!
      call dgemm('N', 'N',                      &
                  (wf%n_v)*(wf%n_o),            &
                  1,                            &
                  (wf%n_v)*(wf%n_o),            &
                  one,                          &
                  X_aibj,                       & ! X_ai_bj
                  (wf%n_v)*(wf%n_o),            &
                  c_ai,                         & ! c_bj
                  (wf%n_v)*(wf%n_o),            &
                  one,                          &
                  rho_ai,                       &
                  (wf%n_v)*(wf%n_o))
!
!     rho_ai += X_bjai * c_bj
!
      call dgemm('T', 'N',                      &
                  (wf%n_v)*(wf%n_o),            &
                  1,                            &
                  (wf%n_v)*(wf%n_o),            &
                  one,                          &
                  X_aibj,                       & ! X_bj_ai
                  (wf%n_v)*(wf%n_o),            &
                  c_ai,                         & ! c_bj
                  (wf%n_v)*(wf%n_o),            &
                  one,                          &
                  rho_ai,                       &
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(X_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine F_ccs_c1_1_ccs
!
end submodule F_ccs
