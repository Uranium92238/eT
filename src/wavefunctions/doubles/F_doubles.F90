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
submodule (doubles_class) F_doubles
!
!!
!!    F submodule (doubles)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
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
!
   implicit none
!
!
contains
!
!
   module subroutine F_x_mu_transformation_doubles(wf, c, rho, x)
!!
!!    F(X) mu transformation
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
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
      use array_utilities, only: scale_diagonal, zero_array
      use reordering, only: squareup, symmetric_sum, packin
!
      implicit none
!
      class(doubles), intent(inout) :: wf
!
      real(dp), dimension(wf%n_t1 + wf%n_t2), intent(in) :: c
      real(dp), dimension(wf%n_t1 + wf%n_t2), intent(out) :: rho
!
      real(dp), dimension(wf%n_t1 + wf%n_t2), intent(in) :: x
!
      real(dp), dimension(:,:), allocatable     :: c_ai, x_ai, rho_ai
      real(dp), dimension(:,:,:,:), allocatable :: c_aibj, x_aibj, rho_aibj 
!
      type(timings), allocatable :: timer
!
      timer = timings('F transformation doubles', 'm')
      call timer%turn_on()
!
      call mem%alloc(c_ai,   wf%n_v, wf%n_o)
      call mem%alloc(rho_ai, wf%n_v, wf%n_o)
      call mem%alloc(x_ai,   wf%n_v, wf%n_o)
!
      call dcopy(wf%n_t1, x, 1, x_ai, 1)
!
      call zero_array(rho_ai, wf%n_v*wf%n_o)
!
      call dcopy(wf%n_t1, c, 1, c_ai, 1)
!
      call wf%F_ccs_a1_1(c_ai, rho_ai, x_ai)
      call wf%F_ccs_b1_1(c_ai, rho_ai, x_ai)
      call wf%F_ccs_c1_1(c_ai, rho_ai, x_ai)
!
      call mem%alloc(x_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(x(wf%n_t1 + 1:), x_aibj, wf%n_t1)
!
      call wf%F_doubles_a1_2(c_ai, rho_ai, x_aibj)
      call wf%F_doubles_b1_2(c_ai, rho_ai, x_aibj)
      call wf%F_doubles_c1_2(c_ai, rho_ai, x_aibj)
!
      call mem%dealloc(x_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(c_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call squareup(c(wf%n_t1 + 1:), c_aibj, wf%n_t1)
!
      call scale_diagonal(two, c_aibj, wf%n_t1)
!
      call wf%F_doubles_a1_1(c_aibj, rho_ai, x_ai)
!
      call mem%dealloc(c_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dcopy(wf%n_t1, rho_ai, 1, rho, 1)
      call mem%dealloc(rho_ai, wf%n_v, wf%n_o)
!
      call mem%alloc(rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call zero_array(rho_aibj, (wf%n_v**2)*(wf%n_o**2))
!
      call wf%F_doubles_a2_1(c_ai, rho_aibj, x_ai)
!
      call mem%dealloc(x_ai, wf%n_v, wf%n_o)
      call mem%dealloc(c_ai, wf%n_v, wf%n_o)
!
      call symmetric_sum(rho_aibj, wf%n_t1)
!
      call packin(rho(wf%n_t1+1:), rho_aibj, wf%n_t1)
!
      call mem%dealloc(rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine F_x_mu_transformation_doubles
!
!
   module subroutine F_doubles_a1_1_doubles(wf, c_aibj, rho_ai, tbar_ai)
!!
!!    F transformation A1,1 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_A1,1_ai = 2 L_iajb * c_bjck * tbar_ck
!!                - (L_jbic * tbar_ak + L_iajc * tbar_bk + L_jbka * tbar_ci) * c_bjck
!!
!!    Edited by A. K. Schnack-Petersen og E. F. Kjønstad, Sep 2021
!!
      use array_utilities, only: zero_array
      use reordering, only: add_2143_to_1234, add_4123_to_1234, sort_1234_to_3214
!
      implicit none
!
      class(doubles), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)              :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                 :: tbar_ai
!
!     Local variables
!
      real(dp), dimension(:,:,:,:), allocatable :: g_iajb
      real(dp), dimension(:,:,:,:), allocatable :: L_aibj
      real(dp), dimension(:,:,:,:), allocatable :: c_cjbk
!
      real(dp), dimension(:,:), allocatable :: X_bj
      real(dp), dimension(:,:), allocatable :: X_ac
      real(dp), dimension(:,:), allocatable :: X_ki
      real(dp), dimension(:,:), allocatable :: X_cj
!
!     L_iajb = 2 g_iajb - g_jaib (ordered as L_aibj)
!
      call mem%alloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%eri_t1%get('ovov',g_iajb)
!
      call mem%alloc(L_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call zero_array(L_aibj, (wf%n_v**2)*(wf%n_o**2))
      call add_2143_to_1234(two, g_iajb, L_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call add_4123_to_1234(-one, g_iajb, L_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     Term 1: 2 * L_iajb * c_bjck * tbar_ck
!
!     X_bj = c_bjck * tbar_ck
!
      call mem%alloc(X_bj, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',             &
                  (wf%n_v)*(wf%n_o),   &
                  1,                   &
                  (wf%n_v)*(wf%n_o),   &
                  one,                 &
                  c_aibj,              & ! c_bj,ck
                  (wf%n_v)*(wf%n_o),   &
                  tbar_ai,             &
                  (wf%n_v)*(wf%n_o),   &
                  zero,                &
                  X_bj,                &
                  (wf%n_v)*(wf%n_o))
!
!     rho_ai += 2 * L_aibj * X_bj
!
      call dgemm('N', 'N',             &
                  (wf%n_v)*(wf%n_o),   &
                  1,                   &
                  (wf%n_v)*(wf%n_o),   &
                  two,                 &
                  L_aibj,              & ! L_ai,bj
                  (wf%n_v)*(wf%n_o),   &
                  X_bj,                &
                  (wf%n_v)*(wf%n_o),   &
                  one,                 &
                  rho_ai,              &
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(X_bj, wf%n_v, wf%n_o)
!
!     Term 2: - L_jbic * tbar_ak * c_bjck
!
!     X_ki = c_bjck L_bjci
!
      call mem%alloc(X_ki, wf%n_o, wf%n_o)
!
      call dgemm('T', 'N',             &
                  wf%n_o,              &
                  wf%n_o,              &
                  (wf%n_v**2)*wf%n_o,  &
                  one,                 &
                  c_aibj,              & ! c_bjc,k
                  (wf%n_v**2)*wf%n_o,  &
                  L_aibj,              & ! L_bjc,i
                  (wf%n_v**2)*wf%n_o,  &
                  zero,                &
                  X_ki,                &
                  wf%n_o)
!
!     rho_ai -=  tbar_ak * X_ki
!
      call dgemm('N', 'N',             &
                  wf%n_v,              &
                  wf%n_o,              &
                  wf%n_o,              &
                  -one,                &
                  tbar_ai,             & ! tbar_a,k
                  wf%n_v,              &
                  X_ki,                &
                  wf%n_o,              &
                  one,                 &
                  rho_ai,              &
                  wf%n_v)

!
      call mem%dealloc(X_ki, wf%n_o, wf%n_o)
!
!     Term 3: -  L_iajc * tbar_bk * c_bjck
!
!     Reorder c_bjck to c_cjbk
!
      call mem%alloc(c_cjbk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_3214(c_aibj, c_cjbk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(X_cj, wf%n_v, wf%n_o)
!
!     X_cj = c_cjbk t_bk
!
      call dgemm('N', 'N',             &
                  (wf%n_v)*(wf%n_o),   &
                  1,                   &
                  (wf%n_v)*(wf%n_o),   &
                  one,                 &
                  c_cjbk,              & ! c_cj,bk
                  (wf%n_v)*(wf%n_o),   &
                  tbar_ai,             & ! t_bk
                  (wf%n_v)*(wf%n_o),   &
                  zero,                &
                  X_cj,                &
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(c_cjbk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     rho_ai += L_aicj * X_cj
!
      call dgemm('N', 'N',             &
                  (wf%n_v)*(wf%n_o),   &
                  1,                   &
                  (wf%n_v)*(wf%n_o),   &
                  -one,                &
                  L_aibj,              & ! L_ai,ck
                  (wf%n_v)*(wf%n_o),   &
                  X_cj,                &
                  (wf%n_v)*(wf%n_o),   &
                  one,                 &
                  rho_ai,              &
                  (wf%n_v)*(wf%n_o))

!
      call mem%dealloc(X_cj, wf%n_v, wf%n_o)

!
!     Term 4: - L_jbka * tbar_ci * c_bjck
!
!     X_ac = L_jbak * c_bjck
!
      call mem%alloc(X_ac, wf%n_v, wf%n_v)
!
      call dgemm('N', 'T',             &
                  wf%n_v,              &
                  wf%n_v,              &
                  (wf%n_o**2)*wf%n_v,  &
                  one,                 &
                  L_aibj,              &  ! L_a,kbj
                  wf%n_v,              &
                  c_aibj,              &  ! c_c,kbj
                  wf%n_v,              &
                  zero,                &
                  X_ac,                &
                  wf%n_v)
!
      call mem%dealloc(L_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     rho_ai += X_ac * tbar_ci
!
      call dgemm('N', 'N',             &
                  wf%n_v,              &
                  wf%n_o,              &
                  wf%n_v,              &
                  -one,                &
                  X_ac,                &
                  wf%n_v,              &
                  tbar_ai,             & ! tbar,ci
                  wf%n_v,              &
                  one,                 &
                  rho_ai,              &
                  wf%n_v)
!
      call mem%dealloc(X_ac, wf%n_v, wf%n_v)
!
   end subroutine F_doubles_a1_1_doubles
!
!
   module subroutine F_doubles_a2_1_doubles(wf, c_ai, rho_aibj, tbar_ai)
!!
!!    F transformation A2,1 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_A2,1_aibj = 2 * L_iakc * tbar_bj * c_ck
!!                   - (L_jbic * tbar_ak + L_jbka * tbar_ci + L_kcib * tbar_aj) c_ck
!!
!!    Edited by A. K. Schnack-Petersen og E. F. Kjønstad, Sep 2021
!!
      use array_utilities, only: copy_and_scale
      use reordering, only: add_1432_to_1234, sort_12_to_21
      use reordering, only: add_4321_to_1234, add_2143_to_1234
!
      implicit none
!
      class(doubles), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)      :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout)   :: rho_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)      :: tbar_ai
!
      real(dp), dimension(:,:,:,:), allocatable :: L_iakc, g_iakc
      real(dp), dimension(:,:,:,:), allocatable :: rho_jbia, rho_iajb
      real(dp), dimension(:,:,:,:), allocatable :: X_jbik
!
      real(dp), dimension(:,:), allocatable :: X_ia, X_ib, X_ik, c_kc
!
      integer :: a, i, b, j
!
!     Construct L_iakc = 2 g_iakc - g_icka
!
      call mem%alloc(L_iakc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%alloc(g_iakc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%eri_t1%get('ovov',g_iakc)
!
      call copy_and_scale(two, g_iakc, L_iakc, wf%n_o**2*wf%n_v**2)
      call add_1432_to_1234(-one, g_iakc, L_iakc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%dealloc(g_iakc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     :: Term 1: 2 * L_iakc * tbar_bj * c_ck
!
!     X_ia = 2 * L_iakc * c_ck
!
      call mem%alloc(c_kc, wf%n_o, wf%n_v)
      call sort_12_to_21(c_ai, c_kc, wf%n_v, wf%n_o)
!
      call mem%alloc(X_ia, wf%n_o, wf%n_v)
!
      call dgemm('N','N',               &
                  (wf%n_o)*(wf%n_v),    &
                  1,                    &
                  (wf%n_o)*(wf%n_v),    &
                  two,                  &
                  L_iakc,               &
                  (wf%n_o)*(wf%n_v),    &
                  c_kc,                 &
                  (wf%n_o)*(wf%n_v),    &
                  zero,                 &
                  X_ia,                 &
                  (wf%n_o)*(wf%n_v))
!
!     rho_aibj =+ X_ia * tbar_bj
!
!     Collapse(2) is a workaround for Intel compiler bug
!     that somehow messes up the loop with -O3 when there
!     is not collapse statement
!$omp parallel do private(j,b,i,a) collapse(2)
      do j = 1, wf%n_o
         do b = 1, wf%n_v
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  rho_aibj(a,i,b,j) = rho_aibj(a,i,b,j) + X_ia(i,a)*tbar_ai(b,j)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(X_ia, wf%n_o, wf%n_v)
!
!     :: Term 2: - L_jbic * tbar_ak * c_ck
!
!     X_jbik = - L_jbic * c_ck
!
      call mem%alloc(X_jbik, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N','N',               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_o,               &
                  wf%n_v,               &
                  -one,                 &
                  L_iakc,               & ! L_jbi,c
                  (wf%n_v)*(wf%n_o)**2, &
                  c_ai,                 & ! c_c,k
                  wf%n_v,               &
                  zero,                 &
                  X_jbik,               & ! X_jbi,k
                  (wf%n_v)*(wf%n_o)**2)
!
!     rho_jbia = X_jbik * tbar_ak
!
      call mem%alloc(rho_jbia, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call dgemm('N','T',               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_v,               &
                  wf%n_o,               &
                  one,                  &
                  X_jbik,               & ! X_jbi,k
                  (wf%n_v)*(wf%n_o)**2, &
                  tbar_ai,              & ! tbar_a,k
                  wf%n_v,               &
                  zero,                 &
                  rho_jbia,             & ! rho_jbi,a
                  (wf%n_v)*(wf%n_o)**2)
!
      call mem%dealloc(X_jbik, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
!     rho_aibj =+ rho_jbia
!
      call add_4321_to_1234(one, rho_jbia, rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(rho_jbia, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     :: Term 3: - L_kajb * tbar_ci * c_ck
!
!     X_ik = -tbar_ci * c_ck
!
      call mem%alloc(X_ik, wf%n_o, wf%n_o)
!
      call dgemm('T','N',               &
                  wf%n_o,               &
                  wf%n_o,               &
                  wf%n_v,               &
                  -one,                 &
                  tbar_ai,              & ! tbar_c,i
                  wf%n_v,               &
                  c_ai,                 & ! c_c,k
                  wf%n_v,               &
                  zero,                 &
                  X_ik,                 &
                  wf%n_o)
!
!     rho_iajb = X_ik * L_kajb
!
      call mem%alloc(rho_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call dgemm('N','N',               &
                  wf%n_o,               &
                  (wf%n_o)*(wf%n_v)**2, &
                  wf%n_o,               &
                  one,                  &
                  X_ik,                 & ! X_i,k
                  wf%n_o,               &
                  L_iakc,               & ! L_k,ajb
                  wf%n_o,               &
                  zero,                 &
                  rho_iajb,             & ! rho_i,ajb
                  wf%n_o)
!
      call mem%dealloc(X_ik, wf%n_o, wf%n_o)
!
!     rho_aibj =+ rho_iajb
!
      call add_2143_to_1234(one, rho_iajb, rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(rho_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     :: Term 4: - L_kcib * tbar_aj * c_ck
!
!     X_ib = - c_ck * L_kcib
!
      call mem%alloc(X_ib, wf%n_o, wf%n_v)
!
      call dgemm('N','N',               &
                  1,                    &
                  (wf%n_o)*(wf%n_v),    &
                  (wf%n_o)*(wf%n_v),    &
                  -one,                 &
                  c_kc,                 & ! c_1,kc
                  1,                    &
                  L_iakc,               & ! L_kc,ib
                  (wf%n_o)*(wf%n_v),    &
                  zero,                 &
                  X_ib,                 & ! X_1,ib
                  1)
!
!     rho_aibj =+ X_ib * tbar_aj
!
!     Collapse(2) is a workaround for Intel (2021) compiler bug
!     that somehow messes up the loop with -O3 when there
!     is not collapse statement
!$omp parallel do private(j,b,i,a) collapse(2)
      do j = 1, wf%n_o
         do b = 1, wf%n_v
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  rho_aibj(a,i,b,j) = rho_aibj(a,i,b,j) + X_ib(i,b)*tbar_ai(a,j)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(X_ib, wf%n_o, wf%n_v)
      call mem%dealloc(c_kc, wf%n_o, wf%n_v)
      call mem%dealloc(L_iakc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
   end subroutine F_doubles_a2_1_doubles
!
!
   module subroutine F_doubles_a1_2_doubles(wf, c_ai, rho_ai, tbar_aibj)
!!
!!    F transformation A1,2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_A1,1 = - (g_ibck * tbar_ajck + g_jack * tbar_bick) * c_bj
!!
!!    Edited by A. K. Schnack-Petersen og E. F. Kjønstad, Sep 2021
!!
      use reordering, only: sort_1234_to_4123, sort_1234_to_2341, add_21_to_12
!
      implicit none
!
      class(doubles), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                   :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)                :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)   :: tbar_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ckib
      real(dp), dimension(:,:,:,:), allocatable :: X_ckij, X_jcki, X_jick, X_ickj
!
      real(dp), dimension(:,:), allocatable :: rho_ia
!
!     :: Term 1: - g_ckib * tbar_ajck * c_bj
!
      call mem%alloc(g_ckib, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call wf%eri_t1%get('voov',g_ckib)
!
!     X_ckij = - g_ckib * c_bj
!
      call mem%alloc(X_ckij, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call dgemm('N','N',               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_o,               &
                  wf%n_v,               &
                  -one,                 &
                  g_ckib,               & ! g_cki,b
                  (wf%n_v)*(wf%n_o)**2, &
                  c_ai,                 & ! c_b,j
                  wf%n_v,               &
                  zero,                 &
                  X_ckij,               & ! X_cki,j
                  (wf%n_v)*(wf%n_o)**2)
!
!     rho_ai = tbar_ajck * X_ckij
!
      call mem%alloc(X_jcki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call sort_1234_to_4123(X_ckij, X_jcki, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call mem%dealloc(X_ckij, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  wf%n_o,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  tbar_aibj,            & ! tbar_a,jck
                  wf%n_v,               &
                  X_jcki,               & ! X_jck,i
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  rho_ai,               & ! rho_a,i
                  wf%n_v)
!
      call mem%dealloc(X_jcki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
!     :: Term 2: - g_ckja * tbar_bick * c_bj
!
!     X_jick = - c_bj * tbar_bick
!
      call mem%alloc(X_jick, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T','N',               &
                  wf%n_o,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_v,               &
                  -one,                 &
                  c_ai,                 & ! c_b,j
                  wf%n_v,               &
                  tbar_aibj,            & ! tbar_b,ick
                  wf%n_v,               &
                  zero,                 &
                  X_jick,               & ! X_j,ick
                  wf%n_o)
!
!     rho_ia = X_jick * g_ckja
!
      call mem%alloc(X_ickj, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call sort_1234_to_2341(X_jick, X_ickj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(X_jick, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(rho_ia, wf%n_o, wf%n_v)
!
      call dgemm('N','N',               &
                  wf%n_o,               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  X_ickj,               & ! X_i,ckj
                  wf%n_o,               &
                  g_ckib,               & ! g_ckj,a
                  (wf%n_v)*(wf%n_o)**2, &
                  zero,                 &
                  rho_ia,               & ! rho_i,a
                  wf%n_o)
!
      call add_21_to_12(one, rho_ia, rho_ai, wf%n_v, wf%n_o)
!
      call mem%dealloc(X_ickj, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(g_ckib, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call mem%dealloc(rho_ia, wf%n_o, wf%n_v)
!
   end subroutine F_doubles_a1_2_doubles
!
!
   module subroutine F_doubles_b1_2_doubles(wf, c_ai, rho_ai, tbar_aibj)
!!
!!    F transformation B1,2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_B1,1 = - (g_ikcb * tbar_akcj + g_jkca * tbar_bkci) * c_bj
!!
!!    Edited by A. K. Schnack-Petersen og E. F. Kjønstad, Sep 2021
!!
      implicit none
!
      class(doubles), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                   :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)                :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)   :: tbar_aibj
!
!     Local variables
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ikcb
      real(dp), dimension(:,:,:,:), allocatable :: X_ikcj
      real(dp), dimension(:,:,:,:), allocatable :: X_jkci
!
!     Term 1: - g_ikcb * tbar_akcj * c_bj
!
      call mem%alloc(g_ikcb, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call wf%eri_t1%get('oovv',g_ikcb)
!
!     X_ikcj = g_ikcb * c_bj
!
      call mem%alloc(X_ikcj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',              &
                  (wf%n_o**2)*(wf%n_v), &
                  wf%n_o,               &
                  wf%n_v,               &
                  one,                  &
                  g_ikcb,               & ! g_ikc,b
                  (wf%n_o**2)*(wf%n_v), &
                  c_ai,                 & ! c_b,j
                  wf%n_v,               &
                  zero,                 &
                  X_ikcj,               &
                  (wf%n_o**2)*(wf%n_v))
!
!     rho_ai -= tbar_akcj * X_ikcj^T
!
      call dgemm('N', 'T',              &
                  wf%n_v,               &
                  wf%n_o,               &
                  (wf%n_o**2)*(wf%n_v), &
                  -one,                 &
                  tbar_aibj,            & ! tbar_a,kcj
                  wf%n_v,               &
                  X_ikcj,               & ! X_i,kcj
                  wf%n_o,               &
                  one,                  &
                  rho_ai,               &
                  wf%n_v)
!
      call mem%dealloc(X_ikcj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     Term 2:  - g_jkca * tbar_bkci * c_bj
!
!     X_jkci = c_bj^T * tbar_bkci
!
      call mem%alloc(X_jkci, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'N',              &
                  wf%n_o,               &
                  (wf%n_o**2)*wf%n_v,   &
                  wf%n_v,               &
                  one,                  &
                  c_ai,                 & ! c_b,j
                  wf%n_v,               &
                  tbar_aibj,            & ! t_b,kci
                  wf%n_v,               &
                  zero,                 &
                  X_jkci,               & ! X_j,kci
                  wf%n_o)
!
!     rho_ai -= g_jkca^T * X_jkci
!
      call dgemm('T', 'N',              &
                  wf%n_v,               &
                  wf%n_o,               &
                  (wf%n_o**2)*wf%n_v,   &
                  -one,                 &
                  g_ikcb,               & ! g_jkc,a
                  (wf%n_o**2)*wf%n_v,   &
                  X_jkci,               & ! X_jkc,i
                  (wf%n_o**2)*wf%n_v,   &
                  one,                  &
                  rho_ai,               &
                  wf%n_v)
!
      call mem%dealloc(g_ikcb, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call mem%dealloc(X_jkci, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine F_doubles_b1_2_doubles
!
!
   module subroutine F_doubles_c1_2_doubles(wf, c_ai, rho_ai, tbar_aibj)
!!
!!    F transformation C1,2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_C1,1 = (g_ikjl * tbar_akbl - g_cadb * tbar_cidj) * c_bj
!!
!!    Edited by A. K. Schnack-Petersen og E. F. Kjønstad, Sep 2021
!!
      use batching_index_class, only: batching_index
      use reordering, only: sort_1234_to_3412, sort_1234_to_3124
!
      implicit none
!
      class(doubles), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                   :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)                :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)   :: tbar_aibj
!
!     Local variables
!
      real(dp), dimension(:,:,:,:), allocatable :: X_jlak
      real(dp), dimension(:,:,:,:), allocatable :: X_akjl
      real(dp), dimension(:,:,:,:), allocatable :: g_ikjl
      real(dp), dimension(:,:,:,:), allocatable :: X_bcdi
      real(dp), dimension(:,:,:,:), allocatable :: X_dbci
      real(dp), dimension(:,:,:,:), allocatable :: tbar_jcdi
      real(dp), dimension(:,:,:,:), allocatable :: g_dbca
!
      type(batching_index) :: batch_c, batch_d
!
      integer :: req0, req1_c, req1_d, req2, current_c_batch, current_d_batch
      integer :: c, d, i, j
!
!     Term 1:  g_ikjl * tbar_blak * c_bj
!
!     X_jlak = c_bj^T * tbar_blak
!
      call mem%alloc(X_jlak, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'N',             &
                  wf%n_o,              &
                  (wf%n_o**2)*wf%n_v,  &
                  wf%n_v,              &
                  one,                 &
                  c_ai,                & ! c_b,j
                  wf%n_v,              &
                  tbar_aibj,           & ! tbar_b,lak
                  wf%n_v,              &
                  zero,                &
                  X_jlak,              & ! X_j,lak
                  wf%n_o)
!
!     Reorder X_jlak to X_akjl
!
      call mem%alloc(X_akjl, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call sort_1234_to_3412(X_jlak, X_akjl, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(X_jlak, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(g_ikjl, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
      call wf%eri_t1%get('oooo',g_ikjl)
!
!     rho_ai += X_akjl * g_ikjl^T
!
      call dgemm('N', 'T',             &
                  wf%n_v,              &
                  wf%n_o,              &
                  wf%n_o**3,           &
                  one,                 &
                  X_akjl,              &
                  wf%n_v,              &
                  g_ikjl,              &
                  wf%n_o,              &
                  one,                 &
                  rho_ai,              &
                  wf%n_v)
!
      call mem%dealloc(X_akjl, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%dealloc(g_ikjl, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
!     Term 2: - g_cadb * tbar_cidj * c_bj
!
      req0 = 0
!
      req1_c =  wf%eri_t1%n_J*wf%n_v
      req1_d =  wf%eri_t1%n_J*wf%n_v
!
      req2 = (wf%n_o)*(wf%n_v) + wf%n_v**2
!
      batch_c = batching_index(wf%n_v)
      batch_d = batching_index(wf%n_v)
!
      call mem%batch_setup(batch_c, batch_d, req0, req1_c, req1_d, req2, &
                           tag='F_doubles_c1_2_doubles')
!
      do current_c_batch = 1, batch_c%num_batches
!
         call batch_c%determine_limits(current_c_batch)
!
         do current_d_batch = 1, batch_d%num_batches
!
            call batch_d%determine_limits(current_d_batch)
!
            call mem%alloc(tbar_jcdi, wf%n_o, batch_c%length, batch_d%length, wf%n_o)
!
!           Reorder tbar_cidj to tbar_jcdi
!
!$omp parallel do private(i, d, c, j)
            do i = 1, wf%n_o
               do d = 1, batch_d%length
                  do c = 1, batch_c%length
                     do j = 1, wf%n_o
!
                        tbar_jcdi(j, c, d, i) = tbar_aibj(c + batch_c%first - 1, i, d + batch_d%first - 1, j)
!
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
!           X_bcdi = c_bj * tbar_jcdi ( = tbar_cidj * c_bj )
!
            call mem%alloc(X_bcdi, wf%n_v, batch_c%length, batch_d%length, wf%n_o)
!
            call dgemm('N', 'N',                                     &
                        wf%n_v,                                      &
                        (batch_d%length)*(batch_c%length)*(wf%n_o),  &
                        wf%n_o,                                      &
                        one,                                         &
                        c_ai,                                        & ! c_bj
                        wf%n_v,                                      &
                        tbar_jcdi,                                   & ! tbar_j,cdi
                        wf%n_o,                                      &
                        zero,                                        &
                        X_bcdi,                                      & ! X_b,cdi
                        wf%n_v)
!
            call mem%dealloc(tbar_jcdi, wf%n_o, batch_c%length, batch_d%length, wf%n_o)
!
!           Reorder X_bcdi to X_dbci
!
            call mem%alloc(X_dbci, batch_d%length, wf%n_v, batch_c%length, wf%n_o)
            call sort_1234_to_3124(X_bcdi, X_dbci, wf%n_v, batch_c%length, batch_d%length, wf%n_o)
            call mem%dealloc(X_bcdi, wf%n_v, batch_c%length, batch_d%length, wf%n_o)
!
            call mem%alloc(g_dbca, batch_d%length, wf%n_v, batch_c%length, wf%n_v)
            call wf%eri_t1%get('vvvv',g_dbca,                               &
                                   batch_d%first, batch_d%get_last(),           &
                                   1, wf%n_v,                                   &
                                   batch_c%first, batch_c%get_last(),           &
                                   1, wf%n_v)
!
            call dgemm('T', 'N',                                     &
                        wf%n_v,                                      &
                        wf%n_o,                                      &
                        (batch_d%length)*(batch_c%length)*(wf%n_v),  &
                        one,                                         &
                        g_dbca,                                      & ! g_dbc,a
                        (batch_d%length)*(batch_c%length)*(wf%n_v),  &
                        X_dbci,                                      & ! X_dbc,i
                        (batch_d%length)*(batch_c%length)*(wf%n_v),  &
                        one,                                         &
                        rho_ai,                                      &
                        wf%n_v)
!
            call mem%dealloc(g_dbca, batch_d%length, wf%n_v, batch_c%length, wf%n_v)
            call mem%dealloc(X_dbci, batch_d%length, wf%n_v, batch_c%length, wf%n_o)
!
         enddo ! batch_d
      enddo ! batch_c
!
      call mem%batch_finalize()
!
   end subroutine F_doubles_c1_2_doubles
!
!
end submodule F_doubles
