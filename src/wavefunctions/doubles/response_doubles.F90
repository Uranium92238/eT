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
submodule (doubles_class) response_doubles
!
!!
!!    Response properties submodule
!!
!!    Routines for construction of the right-hand-side, eta^X,
!!    and left-hand-side, xi^X vectors and the left-hand-side (D^L)
!!    and right-hand-side (D^R) transition densities for transition moments.
!!
!!    Equation-of-motion (EOM):
!!
!!    (Following Koch, H., Kobayashi, R., Sanches de Merás, A., and Jørgensen, P.,
!!    J. Chem. Phys. 100, 4393 (1994))
!!
!!          eta_mu^X,EOM = < Lambda| [X, tau_mu] |CC >
!!                         + (< Lambda| tau_mu X |CC > - tbar_mu < Lambda| X |CC >)
!!                       = eta^{X,0} + eta^{X,corr}
!!
!!    Where the last two terms are called the EOM-corrections and the first term also
!!    appears in LR-CC.
!!
!!    The left-hand-side vector is the same in EOM-CC and LR-CC:
!!
!!          xi^X_mu = < mu| exp(-T) X exp(T)|R >
!!
!!    Density Matrices:
!!
!!    In general a CC density matrix can be written as:
!!
!!          D_pq = < X| e^(-T) E_pq e^T |Y >
!!
!!    where X and Y are left and right state vectors with contributions
!!    from a reference determinant and excited determinants (< mu|, |nu >):
!!
!!          D_pq =             X_ref < HF| e^(-T) E_pq e^T |HF >  Y_ref
!!                 + sum_mu    X_mu  < mu| e^(-T) E_pq e^T |HF >  Y_ref
!!                 + sum_mu    X_ref < HF| e^(-T) E_pq e^T |mu >  Y_mu
!!                 + sum_mu,nu X_mu  < mu| e^(-T) E_pq e^T |nu >  Y_nu
!!
!!    Depending on the type of density matrix (Ground state, transition ,
!!    excited state, interstate transition) different states and thus different
!!    amplitudes X_ref, X_mu, Y_ref and Y_mu will contribute.
!!
!!    In EOM theory the states can be written as the following vectors:
!!
!!          |CC >     = R_0 = (1, 0)
!!          |Lambda > = L_0 = (1, tbar_mu)
!!          |R_k >    = R_k = (-sum_mu(tbar_mu*R_mu), R_mu)
!!          |L_k >    = L_k = (0, L_mu)
!!
!!    The routine names derive from the contribution of the vectors:
!!
!!       ref_ref: first component of the vector for the left and right state
!!
!!       mu_ref:  second component of the vector for the left and
!!                first component of the vector for the right state
!!
!!       ref_mu:  first component of the vector for the left and
!!                second component of the vector for the right state
!!
!!       mu_nu:   second component of the vector for the left and right state
!!
!!
!!    The EOM transition density matrices are constructed as follows:
!!
!!          D^L_pq = < k| E_pq |CC >
!!          D^R_pq = < Lambda| E_pq |k >
!!
!!    where |k > and < k| are the eigenvectors of the Jacobian
!!    with the amplitudes R_mu, L_mu
!!
!!          |k > = - tbar R_k |CC > + sum_mu (tau_mu |CC > R_{k,mu})
!!          < k| = sum_mu L_{k,mu} < mu| e^-T
!!
!!    For the left transition density all the ground state terms can be reused,
!!    if tbar is replaced by L_k and the ref_ref term is neglected.
!!
!
   implicit none
!
!
contains
!
!
   module subroutine construct_right_transition_density_doubles(wf, density, state, R)
!!
!!    Construct right transition density (EOM)
!!    Written by Alexander C. Paul, June 2019
!!
!!          D^R_pq = < Lambda| E_pq |R_n >
!!
!!    where |R_n > is the right eigenvector of the Jacobian
!!    with amplitudes R_mu
!!
!!          |R_n > = sum_mu (tau_mu R_{n,mu} - tbar_mu R_{n,mu}) |CC >
!!                 = (r0 + sum_mu tau_mu R_{n,mu}) |CC >
!!
!!    Contributions to the right transition density are split as follows:
!!
!!    D^R_pq = sum_mu < HF| E_pq |mu > R_mu + sum_mu tbar_mu < mu| E_pq |HF >
!!           + sum_mu,nu tbar_mu < mu| E_pq |nu > R_nu
!!
!!    The last term is separated in "mu_nu_density_terms" as it is identical
!!    for the right transition density and excited state densities
!!
      use array_utilities, only: get_trace
!
      implicit none
!
      class(doubles) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(out) :: density
!
      integer, intent(in) :: state
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: R
!
      type(timings)     :: timer
      character(len=25) :: timer_name
!
      write(timer_name, '(a,i0,a)') 'Density <0|E_pq|',state,'>'
      timer = timings(trim(timer_name), pl='m')
      call timer%turn_on()
!
      call wf%mu_nu_density_terms(density, 0, [wf%t1bar, wf%t2bar], &
                               state, wf%r0(state), R)
!
      call wf%density_ccs_ref_mu_ov(density, R)
      call wf%density_mu_mu_oo(density, -wf%r0(state))
      call wf%density_mu_ref(density, wf%density, wf%r0(state))
!
      call timer%turn_off()
!
      call output%printf('debug', 'Trace (a0): (f15.12)', chars=[trim(timer_name)], &
                          reals=[get_trace(density, wf%n_mo)], fs='(/t6,a)')
!
   end subroutine construct_right_transition_density_doubles
!
!
   module subroutine mu_nu_density_terms_doubles(wf, density, m, L, n, r0, R)
!!
!!    density mu nu terms
!!    Written by Alexander C. Paul, May 2021
!!
!!    Constructs terms of the form:
!!       sum_mu,nu L_mu < mu| E_pq |nu > R_nu
!!
!!    corresponding to terms of the right transition density
!!    and excited state densities.
!!
!!    NB: Terms where mu == nu are separated out in construct_es_density
!!        and construct_right_transition_density
!!
      use array_utilities, only: scale_diagonal, zero_array
      use reordering, only: squareup
!
      implicit none
!
      class(doubles) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(out) :: density
!
      integer, intent(in) :: m, n
!
!     L might only be contiguous in the ranges (1:n_t1) and (1+t1: n_t1+n_t2)
!     as L can also be a combined array of t1bar + t2bar
      real(dp), dimension(wf%n_t1+wf%n_t2), intent(in) :: L
!
      real(dp), intent(out) :: r0
      real(dp), dimension(wf%n_t1+wf%n_t2), intent(in) :: R
!
      real(dp), dimension(:,:,:,:), allocatable :: L2, R2
!
      real(dp) :: ddot
!
      type(timings)     :: timer
      character(len=40) :: timer_name
!
      write(timer_name, '(a,i0,a,i0,a)') 'Doubles contribution to <', m, '|E_pq|',n,'>'
      timer = timings(trim(timer_name), pl='v')
!
      call zero_array(density, wf%n_mo**2)
      r0 = zero
!
      call wf%ccs%mu_nu_density_terms(density, m, L, n, r0, R)
!
      call timer%turn_on() ! Only doubles contribution
!
      call mem%alloc(L2, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(L(wf%n_t1+1:), L2, wf%n_v*wf%n_o)
!
      call wf%density_doubles_mu_nu_ov(density, L2, R(1:wf%n_t1))
      call wf%density_doubles_mu_nu_vo(density, L2, R(1:wf%n_t1))
!
      call mem%alloc(R2, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(R(wf%n_t1 + 1 : wf%n_es_amplitudes), R2, wf%n_t1)
      call scale_diagonal(two, R2, wf%n_t1)
!
      r0 = r0 - half*ddot(wf%n_t1**2, R2, 1, L2, 1)
!
      call wf%density_doubles_mu_ref_oo(density, L2, R2)
      call wf%density_doubles_mu_ref_vv(density, L2, R2)
!
      call mem%dealloc(L2, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%density_doubles_mu_ref_ov(density, L(1:wf%n_t1), R2)
!
      call mem%dealloc(R2, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine mu_nu_density_terms_doubles
!
!
   module subroutine density_doubles_mu_nu_ov_doubles(wf, density, tbar_aibj, R_ai)
!!
!!    One electron density (EOM) excited-determinant/excited-determinant ov-term
!!    Written by Alexander C. Paul, June 2019
!!
!!    Computes terms of the form:
!!
!!          D_pq += sum_mu,nu X_mu < mu| e^(-T) E_pq e^T |nu > Y_nu
!!
!!    explicit term in this routine:
!!          D^R_kc += sum_abij R^a_i tbar^ab_ij (2t^bc_jk - t^bc_kj)
!!                   -sum_abij tbar^ab_ij (R^b_k t^ac_ij + R^c_j t^ab_ik)
!!
      use reordering, only: squareup, add_1432_to_1234
!
      implicit none
!
      class(doubles) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: density
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: tbar_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: R_ai
!
      real(dp), dimension(:,:,:,:), allocatable :: u_bjck, t_aick, t_aicj
      real(dp), dimension(:,:), allocatable :: X_bj, X_ck, X_jk, X_bc
!
      integer :: k, c
!
!     :: Term 1: sum_abij R^a_i tbar^ab_ij (2t^bc_jk - t^bc_kj) ::
!
      call mem%alloc(X_bj, wf%n_v, wf%n_o)
!
      call dgemv('T',            &
                  wf%n_o*wf%n_v, &
                  wf%n_o*wf%n_v, &
                  one,           &
                  tbar_aibj,     & ! tbar_bj_ai
                  wf%n_o*wf%n_v, &
                  R_ai,          & ! R_ai
                  1,             &
                  zero,          &
                  X_bj,          &
                  1)
!
      call mem%alloc(t_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_aick, wf%n_v*wf%n_o)
!
      call mem%alloc(u_bjck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dcopy((wf%n_v*wf%n_o)**2, t_aick, 1, u_bjck, 1)
      call dscal((wf%n_v*wf%n_o)**2, two, u_bjck, 1)
      call add_1432_to_1234(-one, t_aick, u_bjck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(t_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(X_ck, wf%n_v, wf%n_o)
!
      call dgemv('T',            &
                  wf%n_o*wf%n_v, &
                  wf%n_o*wf%n_v, &
                  one,           &
                  u_bjck,        & ! u_ck_bj
                  wf%n_o*wf%n_v, &
                  X_bj,          & ! X_bj
                  1,             &
                  zero,          &
                  X_ck,          &
                  1)
!
      call mem%dealloc(X_bj, wf%n_v, wf%n_o)
      call mem%dealloc(u_bjck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     :: Term 2: -sum_abij tbar^ab_ij (R^b_k t^ac_ij + R^c_j t^ab_ik) ::
!
      call mem%alloc(t_aicj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_aicj, wf%n_v*wf%n_o)
!
      call mem%alloc(X_bc, wf%n_v, wf%n_v)
!
!     X_bc = sum_aij tbar_aibj t_aicj = sum_abi tbar_bjai t_cjai
!
      call dgemm('N', 'T',          &
                  wf%n_v,           &
                  wf%n_v,           &
                  wf%n_v*wf%n_o**2, &
                  one,              &
                  tbar_aibj,        & ! tbar_b_jai
                  wf%n_v,           &
                  t_aicj,           & ! t_c_jai
                  wf%n_v,           &
                  zero,             &
                  X_bc,             &
                  wf%n_v)
!
!     X_jk = sum_aib tbar_aibj t_aibk
!
      call mem%alloc(X_jk, wf%n_o, wf%n_o)
!
      call dgemm('T', 'N',          &
                  wf%n_o,           &
                  wf%n_o,           &
                  wf%n_o*wf%n_v**2, &
                  one,              &
                  tbar_aibj,        & ! tbar_aib_j
                  wf%n_o*wf%n_v**2, &
                  t_aicj,           & ! t_aib_k
                  wf%n_o*wf%n_v**2, &
                  zero,             &
                  X_jk,             &
                  wf%n_o)
!
      call mem%dealloc(t_aicj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T','N',  &
                  wf%n_v,  &
                  wf%n_o,  &
                  wf%n_v,  &
                  -one,    &
                  X_bc,    & ! X_c_b
                  wf%n_v,  &
                  R_ai,    & ! R_b_k
                  wf%n_v,  &
                  one,     &
                  X_ck,    &
                  wf%n_v)
!
      call mem%dealloc(X_bc, wf%n_v, wf%n_v)
!
      call dgemm('N','N',  &
                  wf%n_v,  &
                  wf%n_o,  &
                  wf%n_o,  &
                  -one,    &
                  R_ai,    & ! R_c_j
                  wf%n_v,  &
                  X_jk,    & ! X_j_k
                  wf%n_o,  &
                  one,     &
                  X_ck,    &
                  wf%n_v)
!
      call mem%dealloc(X_jk, wf%n_o, wf%n_o)
!
!$omp parallel do private(c, k)
      do c = 1, wf%n_v
         do k = 1, wf%n_o
!
            density(k, wf%n_o + c) = X_ck(c, k) + density(k, wf%n_o + c)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(X_ck, wf%n_v, wf%n_o)
!
   end subroutine density_doubles_mu_nu_ov_doubles
!
!
   module subroutine density_doubles_mu_nu_vo_doubles(wf, density, tbar_aibj, R_ai)
!!
!!    One electron density (EOM) excited-determinant/excited-determinant vo-term
!!    Written by Alexander C. Paul, June 2019
!!
!!    Computes terms of the form:
!!
!!          D_pq += sum_mu,nu X_mu < mu| e^(-T) E_pq e^T |nu > Y_nu
!!
!!    explicit term in this routine:
!!          D^R_bj += sum_ai R^a_i tbar^ab_ij
!!
      implicit none
!
      class(doubles) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: density
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: tbar_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: R_ai
!
      real(dp), dimension(:,:), allocatable :: rho_vo
!
      integer :: i, a
!
      call mem%alloc(rho_vo, wf%n_v, wf%n_o)
!
      call dgemv('T',            &
                  wf%n_o*wf%n_v, &
                  wf%n_o*wf%n_v, &
                  one,           &
                  tbar_aibj,     & ! tbar_bj_ai
                  wf%n_o*wf%n_v, &
                  R_ai,          & ! R_ai
                  1,             &
                  zero,          &
                  rho_vo,        &
                  1)
!
!$omp parallel do private(a, i)
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            density(wf%n_o + a, i) = rho_vo(a, i) + density(wf%n_o + a, i)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(rho_vo, wf%n_v, wf%n_o)
!
   end subroutine density_doubles_mu_nu_vo_doubles
!
   module subroutine construct_eom_etaX_doubles(wf, X, xiX, etaX)
!!
!!    Construct EOM etaX
!!    Written by Sarai D. Folkestad, May 2019
!!
!!    Constructs the EOM effective etaX vector, adding the EOM
!!    correction to etaX.
!!
      implicit none
!
      class(doubles), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: xiX
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: etaX
!
      call wf%construct_etaX(X, etaX)
!
      call wf%etaX_eom_a(etaX, xiX)
!
      call wf%etaX_eom_doubles_a1(X, etaX(1:wf%n_t1))
!
   end subroutine construct_eom_etaX_doubles
!
!
   module subroutine construct_etaX_doubles(wf, X, etaX)
!!
!!    Construct eta^X
!!    Written by Josefine H. Andersen, 2019
!!
!!    Adapted by Sarai D. Folekstad, Apr 2019
!!
!!    Constructs left-hand-side vector etaX:
!!
!!       eta^X_mu = < Lambda| [X, tau_mu] |CC >
!!
      use array_utilities, only: zero_array
      use reordering, only: symmetrize_and_add_to_packed
!
      implicit none
!
      class(doubles), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: etaX
!
      real(dp), dimension(:,:), allocatable :: etaX_ai
      real(dp), dimension(:,:,:,:), allocatable :: etaX_aibj
!
      integer :: a, i, ai
!
      call zero_array(etaX, wf%n_es_amplitudes)
!
!     etaX_ai:
!
      call mem%alloc(etaX_ai, wf%n_v, wf%n_o)
      call zero_array(etaX_ai, (wf%n_o*wf%n_v))
!
      call wf%etaX_ccs_a1(X, etaX_ai)
      call wf%etaX_ccs_b1(X, etaX_ai)
!
      call wf%etaX_doubles_a1(X, etaX_ai)
!
!$omp parallel do private (a, i, ai)
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            ai = wf%n_v*(i - 1) + a
!
            etaX(ai) = etaX_ai(a,i)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(etaX_ai, wf%n_v, wf%n_o)
!
!     etaX_aibj:
!
      call mem%alloc(etaX_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call zero_array(etaX_aibj, (wf%n_o*wf%n_v)**2)
!
      call wf%etaX_doubles_a2(X, etaX_aibj)
      call wf%etaX_doubles_b2(X, etaX_aibj)
!
      call symmetrize_and_add_to_packed(etaX(wf%n_t1 + 1 : wf%n_es_amplitudes), &
                                        etaX_aibj, wf%n_v*wf%n_o)
!
      call mem%dealloc(etaX_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine construct_etaX_doubles
!
!
   module subroutine etaX_doubles_a1_doubles(wf, X, etaX_ai)
!!
!!    etaX CCSD A1
!!    Written by Josefine H. Andersen, Feb 2019
!!
!!    Adapted by Sarai D. Folkestad, Apr 2019
!!
!!       A1 = - sum_ckdl (tb_ckal X_id t_ckdl + tb_ckdi X_la t_ckdl)
!!
      use reordering, only: squareup
!
      implicit none
!
      class(doubles), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: etaX_ai
!
      real(dp), dimension(:,:), allocatable :: X_id ! X_la
!
      real(dp), dimension(:,:,:,:), allocatable :: tb_ckal ! tb_ckdi
      real(dp), dimension(:,:,:,:), allocatable :: t_ckdl
!
      real(dp), dimension(:,:), allocatable :: I_ad    ! intermediate, first term
      real(dp), dimension(:,:), allocatable :: I_li    ! intermediate, second term
!
      integer :: d, i
!
      call mem%alloc(X_id, wf%n_o, wf%n_v)
!
!$omp parallel do private(d, i)
      do d = 1, wf%n_v
         do i = 1,wf%n_o
!
            X_id(i, d) = X(i, d + wf%n_o)
!
         enddo
      enddo
!$omp end parallel do
!
!     Squareup multipliers
!
      call mem%alloc(tb_ckal, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2bar, tb_ckal, wf%n_t1)
!
!     Read amplitudes and order as t_lck_d = t_kl^cd
!
      call mem%alloc(t_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_ckdl, wf%n_t1)
!
!     :: First term: - sum_ckdl tb_ckal X_id t_ckdl
!
!     I_a_d = sum_ckl tb_a_lck t_lck_d = sum_ckl tb_ckal t_kl^cd
!
      call mem%alloc(I_ad, wf%n_v, wf%n_v)
!
      call dgemm('N','T',           &
                  wf%n_v,           &
                  wf%n_v,           &
                  wf%n_v*wf%n_o**2, &
                  one,              &
                  tb_ckal,          & ! tbar_a_lck
                  wf%n_v,           &
                  t_ckdl,           & ! t_d_lck
                  wf%n_v,           &
                  zero,             &
                  I_ad,             &
                  wf%n_v)
!
!     Add   - sum_ckdl tb_ckal X_id t_kl^cd
!         = - sum_d I_a_d X_id
!         = - sum_d I_a_d X_i_a^T(d,i)
!
      call dgemm('N','T',     &
                  wf%n_v,     &
                  wf%n_o,     &
                  wf%n_v,     &
                  -one,       &
                  I_ad,       &
                  wf%n_v,     &
                  X_id,       &
                  wf%n_o,     &
                  one,        &
                  etaX_ai,    &
                  wf%n_v)
!
      call mem%dealloc(I_ad, wf%n_v, wf%n_v)
!
!     :: Second term: sum_ckdl tb_ckdi X_la t_ckdl
!
!     X_l_i = sum_ckd t_l_ckd tb_ckd_i  = sum_ckd tb_ckdi t_kl^cd
!
      call mem%alloc(I_li, wf%n_o, wf%n_o)
!
      call dgemm('T','N',           &
                  wf%n_o,           &
                  wf%n_o,           &
                  wf%n_v**2*wf%n_o, &
                  one,              &
                  t_ckdl,           & ! t_ckd_l
                  wf%n_v**2*wf%n_o, &
                  tb_ckal,          & ! tbar_ckd_i
                  wf%n_v**2*wf%n_o, &
                  zero,             &
                  I_li,             &
                  wf%n_o)
!
      call mem%dealloc(t_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(tb_ckal, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Add - sum_ckdl b_ckdi X_la t_kl^cd = - sum_l X_la I_l_i = - sum_l X_i_a^T(a,l) I_l_i(l,i)
!
      call dgemm('T','N',  &
                  wf%n_v,  &
                  wf%n_o,  &
                  wf%n_o,  &
                  -one,    &
                  X_id,    &
                  wf%n_o,  &
                  I_li,    &
                  wf%n_o,  &
                  one,     &
                  etaX_ai, &
                  wf%n_v)
!
      call mem%dealloc(I_li, wf%n_o, wf%n_o)
      call mem%dealloc(X_id, wf%n_o, wf%n_v)
!
   end subroutine etaX_doubles_a1_doubles
!
   module subroutine etaX_doubles_a2_doubles(wf, X, etaX_aibj)
!!
!!    etaX CCSD A2
!!    Written by Josefine H. Andersen, Feb 2019
!!
!!    Adapted by Sarai D. Folkestad, Apr 2019
!!
!!    Constructs the A2 term of the etaX vector
!!
!!       A2 = 2 X_jb tb_ai - X_ib tb_aj
!!
      implicit none
!
      class(doubles), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout) :: etaX_aibj
!
      integer :: i, a, j, b
!
!$omp parallel do private(a, i, b, j)
       do j = 1, wf%n_o
         do b = 1, wf%n_v
!
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  etaX_aibj(a, i, b, j) = etaX_aibj(a, i, b, j) &
                                        + two*X(j, b + wf%n_o)*wf%t1bar(a, i) &
                                        - X(i, b + wf%n_o)*wf%t1bar(a, j)
!
               enddo
            enddo
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine etaX_doubles_a2_doubles
!
!
   module subroutine etaX_doubles_b2_doubles(wf, X, etaX_aibj)
!!
!!    etaX CCSD B2
!!    Written by Josefine H. Andersen, Feb 2019
!!
!!    Adapted by Sarai D. Folkestad, Apr 2019
!!
!!    Constructs the B2 term of the etaX vector
!!
!!       B2 = sum_c tb_aicj X_cb - sum_k tb_aibk X_jk
!!
      use reordering, only: squareup, sort_1234_to_1243, add_1243_to_1234
!
      implicit none
!
      class(doubles), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout) :: etaX_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: etaX_aijb
!
      real(dp), dimension(:,:,:,:), allocatable :: tb_aibj
      real(dp), dimension(:,:,:,:), allocatable :: tb_aijc
!
      real(dp), dimension(:,:), allocatable :: X_cb
      real(dp), dimension(:,:), allocatable :: X_jk
!
      integer :: b, c, k, j
!
!     Get and squareup multipliers
!
      call mem%alloc(tb_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2bar, tb_aibj, wf%n_v*wf%n_o)
!
!     Reorder multipiers to tb_aijc
!
      call mem%alloc(tb_aijc, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call sort_1234_to_1243(tb_aibj, tb_aijc, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     :: First term: sum_c tb_aicj X_cb
!
      call mem%alloc(X_cb, wf%n_v, wf%n_v)
!
!$omp parallel do private(b, c)
      do b = 1, wf%n_v
         do c = 1, wf%n_v
!
            X_cb(c,b) = X(c + wf%n_o, b + wf%n_o)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%alloc(etaX_aijb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N','N',           &
                  wf%n_v*wf%n_o**2, &
                  wf%n_v,           &
                  wf%n_v,           &
                  one,              &
                  tb_aijc,          &
                  wf%n_v*wf%n_o**2, &
                  X_cb,             &
                  wf%n_v,           &
                  zero,             &
                  etaX_aijb,        &
                  wf%n_v*wf%n_o**2)
!
      call mem%dealloc(tb_aijc, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call mem%dealloc(X_cb, wf%n_v, wf%n_v)
!
      call add_1243_to_1234(one, etaX_aijb, etaX_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(etaX_aijb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
!     :: Second term: -sum_k tb_aick X_jk
!
      call mem%alloc(X_jk, wf%n_o, wf%n_o)
!
!$omp parallel do private(j, k)
      do k = 1, wf%n_o
         do j = 1, wf%n_o
!
            X_jk(j,k) = X(j,k)
!
         enddo
      enddo
!$omp end parallel do
!
      call dgemm('N','T',           &
                  wf%n_o*wf%n_v**2, &
                  wf%n_o,           &
                  wf%n_o,           &
                  -one,             &
                  tb_aibj,          &
                  wf%n_o*wf%n_v**2, &
                  X_jk,             &
                  wf%n_o,           &
                  one,              &
                  etaX_aibj,        &
                  wf%n_o*wf%n_v**2)
!
      call mem%dealloc(tb_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X_jk, wf%n_o, wf%n_o)
!
   end subroutine etaX_doubles_b2_doubles
!
!
   module subroutine construct_xiX_doubles(wf, X, xiX)
!!
!!    Construct xiX
!!    Written by Josefine H. Andersen, 2019
!!
!!    Adapted by Sarai D. Folkestad
!!
!!    Constructs xi^X_mu :
!!
!!       xi^X_mu = < mu| exp(-T) X exp(T)|R >
!!
      use array_utilities, only: zero_array
      use reordering, only: symmetrize_and_add_to_packed
!
      implicit none
!
      class(doubles), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: xiX
!
      real(dp), dimension(:,:), allocatable :: xiX_ai
      real(dp), dimension(:,:,:,:), allocatable :: xiX_aibj
!
      integer :: a, i, ai
!
      call zero_array(xiX, wf%n_es_amplitudes)
!
!     xiX_ai
!
      call mem%alloc(xiX_ai, wf%n_v, wf%n_o)
      call zero_array(xiX_ai, (wf%n_o*wf%n_v))
!
      call wf%xiX_ccs_a1(X, xiX_ai)
      call wf%xiX_doubles_a1(X, xiX_ai)
!
!$omp parallel do private (a, i, ai)
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            ai = wf%n_v*(i - 1) + a
!
            xiX(ai) = xiX_ai(a,i)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(xiX_ai, wf%n_v, wf%n_o)
!
!     xiX_aibj
!
      call mem%alloc(xiX_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call zero_array(xiX_aibj, (wf%n_o*wf%n_v)**2)
!
      call wf%xiX_doubles_a2(X, xiX_aibj)
!
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            xiX_aibj(a, i, a, i) = half*xiX_aibj(a, i, a, i)
!
         enddo
      enddo
!
      call symmetrize_and_add_to_packed(xiX(wf%n_t1 + 1 : wf%n_es_amplitudes), &
                                        xiX_aibj, wf%n_v*wf%n_o)
!
      call mem%dealloc(xiX_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine construct_xiX_doubles
!
!
   module subroutine xiX_doubles_a1_doubles(wf, X, xiX_ai)
!!
!!    xiX CCSD A1
!!    Written by Josefine H. Andersen, Feb 2019
!!
!!    Adapted by Sarai D. Folkestad, Apr 2010
!!
!!    Constructs the A1 term of xiX
!!
!!       A1 = sum_ck u_aick X_kc,
!!
!!    where u_aick = 2t_ckai - t_ciak
!!
      use array_utilities, only: zero_array
      use reordering, only: squareup, add_1432_to_1234
!
      implicit none
!
      class(doubles), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: xiX_ai
!
      real(dp), dimension(:,:,:,:), allocatable   :: u_aick
      real(dp), dimension(:,:,:,:), allocatable   :: t_aick
!
      real(dp), dimension(:,:), allocatable   :: X_ck
!
      integer :: k, c
!
!     X_kc ordered as X_ck
!
      call mem%alloc(X_ck, wf%n_v, wf%n_o)
!
!$omp parallel do private(k, c)
      do k = 1, wf%n_o
         do c = 1, wf%n_v
!
            X_ck(c, k) = X(k, c + wf%n_o)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%alloc(t_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_aick,  wf%n_o*wf%n_v)
!
!     Form u_aick = 2 t_ai_ck - t_ak_ci
!
      call mem%alloc(u_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call zero_array(u_aick, (wf%n_o*wf%n_v)**2)
!
      call add_1432_to_1234(-one, t_aick, u_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call daxpy(wf%n_o**2 * wf%n_v**2, two, t_aick, 1, u_aick, 1)
!
      call mem%dealloc(t_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     sum_ck u_ai_ck X_kc
!
      call dgemm('N','N',        &
                  wf%n_o*wf%n_v, &
                  1,             &
                  wf%n_o*wf%n_v, &
                  one,           &
                  u_aick,        &
                  wf%n_o*wf%n_v, &
                  X_ck,          &
                  wf%n_o*wf%n_v, &
                  one,           &
                  xiX_ai,        &
                  wf%n_o*wf%n_v)
!
      call mem%dealloc(u_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X_ck, wf%n_v, wf%n_o)
!
   end subroutine xiX_doubles_a1_doubles
!
!
   module subroutine xiX_doubles_a2_doubles(wf, X, xiX_aibj)
!!
!!    xiX CCSD A2
!!    Written by Josefine H. Andersen, Feb 2019
!!
!!    Adapted by Sarai D. Folekstad
!!
!!    Construct xiX A2 contribution:
!!
!!       A2 = sum_c t_aicj X_bc - sum_k t_aibk X_kj
!!
      use reordering, only: squareup
!
      implicit none
!
      class(doubles), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout) :: xiX_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: t_cjai
!
      real(dp), dimension(:,:), allocatable :: X_bc
      real(dp), dimension(:,:), allocatable :: X_kj
!
      integer :: b, c, j, k
!
      call mem%alloc(t_cjai, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_cjai, wf%n_v*wf%n_o)
!
!     :: First term: sum_c t_aicj X_bc
!
      call mem%alloc(X_bc, wf%n_v, wf%n_v)
!
!$omp parallel do private(b, c)
      do c = 1, wf%n_v
         do b = 1, wf%n_v
!
            X_bc(b, c) = X(b + wf%n_o, c + wf%n_o)
!
         enddo
      enddo
!$omp end parallel do
!
      call dgemm('N','N',          &
                 wf%n_v,           &
                 wf%n_v*wf%n_o**2, &
                 wf%n_v,           &
                 one,              &
                 X_bc,             &
                 wf%n_v,           &
                 t_cjai,           &
                 wf%n_v,           &
                 one,              &
                 xiX_aibj,         & ! xiX_bjai, will symmetrize anyhow
                 wf%n_v)
!
      call mem%dealloc(X_bc, wf%n_v, wf%n_v)
!
!     :: Second term: -sum_k t_ai_bk X_kj
!
      call mem%alloc(X_kj, wf%n_o, wf%n_o)
!
!$omp parallel do private(b, c)
      do j = 1, wf%n_o
         do k = 1, wf%n_o
!
            X_kj(k, j) = X(k, j)
!
         enddo
      enddo
!$omp end parallel do
!
      call dgemm('N','N',          &
                 wf%n_o*wf%n_v**2, &
                 wf%n_o,           &
                 wf%n_o,           &
                 -one,             &
                 t_cjai,           &
                 wf%n_o*wf%n_v**2, &
                 X_kj,             &
                 wf%n_o,           &
                 one,              &
                 xiX_aibj,         &
                 wf%n_o*wf%n_v**2)
!
      call mem%dealloc(X_kj, wf%n_o, wf%n_o)
      call mem%dealloc(t_cjai, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine xiX_doubles_a2_doubles
!
!
   module subroutine etaX_eom_doubles_a1_doubles(wf, X, etaX_ai)
!!
!!    etaX EOM CCSD A1
!!    Written by Josefine H. Andersen, Feb 2019
!!
!!    Adapted by Sarai D. Folkestad
!!
!!    Constructs the A1 correction term to eta^X for EOM
!!
!!       A1 = sum_ck tb_aick X_ck + sum_ckdl tb_aick u_ckdl X_ld,
!!
!!    where u_ckdl = 2*t_ckdl - t_cldk
!!
      use array_utilities, only: zero_array
      use reordering, only: squareup, add_1432_to_1234
!
      implicit none
!
      class(doubles), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: etaX_ai
!
      real(dp), dimension(:,:,:,:), allocatable :: tb_aibj
      real(dp), dimension(:,:,:,:), allocatable :: t_ckdl
      real(dp), dimension(:,:,:,:), allocatable :: u_ckdl
!
      real(dp), dimension(:,:,:,:), allocatable :: I_aidl
!
      real(dp), dimension(:,:), allocatable :: X_ck
      real(dp), dimension(:,:), allocatable :: X_dl
!
      integer :: c, k, d, l
!
!     :: First term: sum_ck tb_aick X_ck
!
      call mem%alloc(tb_aibj, wf%n_v, wf%n_o,  wf%n_v, wf%n_o)
      call squareup(wf%t2bar, tb_aibj, wf%n_v*wf%n_o)
!
      call mem%alloc(X_ck, wf%n_v, wf%n_o)
!
!$omp parallel do private(k, c)
      do k = 1, wf%n_o
         do c = 1, wf%n_v
!
            X_ck(c, k) = X(c + wf%n_o, k)
!
         enddo
      enddo
!$omp end parallel do
!
      call dgemm('N','N',       &
                 wf%n_v*wf%n_o, &
                 1,             &
                 wf%n_v*wf%n_o, &
                 one,           &
                 tb_aibj,       &
                 wf%n_v*wf%n_o, &
                 X_ck,          &
                 wf%n_v*wf%n_o, &
                 one,           &
                 etaX_ai,       &
                 wf%n_v*wf%n_o)
!
      call mem%dealloc(X_ck, wf%n_v, wf%n_o)
!
!     :: Second term: sum_ckdl tb_aick u_ckdl X_ld
!
!     X_ld ordered as X_dl
!
      call mem%alloc(X_dl, wf%n_v, wf%n_o)
!
!$omp parallel do private(l, d)
      do l = 1, wf%n_o
         do d = 1, wf%n_v
!
            X_dl(d, l) = X(l, d + wf%n_o)
!
         enddo
      enddo
!$omp end parallel do
!
!     Form u_aick
!
      call mem%alloc(t_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_ckdl, wf%n_v*wf%n_o)
!
      call mem%alloc(u_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call zero_array(u_ckdl, (wf%n_o*wf%n_v)**2)
!
      call add_1432_to_1234(-one, t_ckdl, u_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call daxpy((wf%n_o*wf%n_v)**2, two, t_ckdl, 1, u_ckdl, 1)
!
      call mem%dealloc(t_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Form the intermediate I_ai_dl = sum_ck tb_ai_ck u_ck_dl
!
      call mem%alloc(I_aidl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',       &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 wf%n_v*wf%n_o, &
                 one,           &
                 tb_aibj,       &
                 wf%n_v*wf%n_o, &
                 u_ckdl,        &
                 wf%n_v*wf%n_o, &
                 zero,          &
                 I_aidl,        &
                 wf%n_v*wf%n_o)
!
      call mem%dealloc(tb_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(u_ckdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Form sum_dl I_ai_dl X_ld^T
!
      call dgemm('N','N',       &
                 wf%n_v*wf%n_o, &
                 1,             &
                 wf%n_v*wf%n_o, &
                 one,           &
                 I_aidl,        &
                 wf%n_v*wf%n_o, &
                 X_dl,          &
                 wf%n_v*wf%n_o, &
                 one,           &
                 etaX_ai,       &
                 wf%n_v*wf%n_o)
!
            call mem%dealloc(I_aidl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
            call mem%dealloc(X_dl, wf%n_v, wf%n_o)
!
   end subroutine etaX_eom_doubles_a1_doubles
!
!
   module subroutine etaX_eom_a_doubles(wf, etaX, xiX)
!!
!!    Get eom contribution
!!    Written by Josefine H. Andersen, Feb 2019
!!
!!    Add EOM contribution to etaX vector
!!
!!    EOM correction:  eta^X,corr_mu += tbar_mu (xi * tbar)
!!
      implicit none
!
      class(doubles), intent(in) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: etaX
      real(dp), dimension(wf%n_es_amplitudes), intent(in)    :: xiX
!
      real(dp) :: X_cc
      real(dp) :: ddot
!
      real(dp), dimension(:), allocatable :: multipliers
!
      call mem%alloc(multipliers, wf%n_es_amplitudes)
!
      call dcopy(wf%n_t1, wf%t1bar, 1, multipliers(1:wf%n_t1), 1)
      call dcopy(wf%n_t2, wf%t2bar, 1, multipliers(wf%n_t1 + 1: wf%n_es_amplitudes), 1)
!
      X_cc = ddot(wf%n_es_amplitudes, multipliers, 1, xiX, 1)
!
      call daxpy(wf%n_es_amplitudes, -X_cc, multipliers, 1, etaX, 1)
!
      call mem%dealloc(multipliers, wf%n_es_amplitudes)
!
   end subroutine etaX_eom_a_doubles
!
!
end submodule response_doubles
