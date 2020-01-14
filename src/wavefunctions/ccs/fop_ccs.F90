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
submodule (ccs_class) fop_ccs
!
!!
!!    First order properties submodule
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
   module subroutine construct_right_transition_density_ccs(wf, state)
!!
!!    Construct right one-electron transition density (EOM)
!!    Written by Alexander C. Paul, June 2019
!!
!!          D^R_pq = < Lambda| E_pq |k >
!!
!!    where |k > is the right eigenvector of the Jacobian
!!    with amplitudes R_mu
!!
!!          |k > = sum_mu (tau_mu R_{k,mu} - tbar_mu R_{k,mu}) |CC >
!!
!!    Contributions to the right transition density are split as follows:
!!
!!          D^R_pq = sum_mu D_pq(ref-mu) R_mu + sum_mu tbar_mu D_pq(mu-ref)
!!                   + sum_mu,nu tbar_mu D_pq(mu-nu) R_nu
!!
      implicit none
!
      class(ccs) :: wf
!
      integer, intent(in) :: state
!
      real(dp), dimension(:), allocatable :: R_k
!
      real(dp) :: ddot, tbar_R_overlap
!
      type(timings) :: R_TDM_timer
!
      R_TDM_timer = timings('Right transition density')
!
      call R_TDM_timer%turn_on()
!
      call zero_array(wf%right_transition_density, (wf%n_mo)**2)
!
      call mem%alloc(R_k, wf%n_es_amplitudes)
      call wf%read_excited_state(R_k, state, 'right')
!
      call wf%density_ccs_mu_nu_oo(wf%right_transition_density, wf%t1bar, R_k)
      call wf%density_ccs_ref_mu_ov(wf%right_transition_density, R_k)
      call wf%density_ccs_mu_nu_vv(wf%right_transition_density, wf%t1bar, R_k)
!
      tbar_R_overlap = ddot(wf%n_t1, wf%t1bar, 1, R_k, 1)
!
      call wf%density_mu_mu_oo(wf%right_transition_density, tbar_R_overlap)
!
      call wf%density_mu_ref(wf%right_transition_density, &
                               wf%density,                &
                               tbar_R_overlap)
!
      call mem%dealloc(R_k, wf%n_es_amplitudes)
!
      call R_TDM_timer%turn_off()
!
   end subroutine construct_right_transition_density_ccs
!
!
   module subroutine density_ccs_mu_nu_oo_ccs(wf, density, tbar_ai, R_ai)
!!
!!    One electron density (EOM) excited-determinant/excited-determinant oo-term 
!!    Written by Alexander C. Paul, June 2019
!!
!!    Computes terms of the form:
!!
!!          D_pq += sum_mu X_mu * < mu| e^(-T) E_pq e^T |nu > Y_nu
!!
!!    explicit term in this routine:
!!          D^R_ij = -sum_a R_ai tbar_aj
!!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: density
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: tbar_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: R_ai
!
      call dgemm('T', 'N', &
                  wf%n_o,  &
                  wf%n_o,  &
                  wf%n_v,  &
                  -one,    &
                  R_ai,    & ! R_a_i
                  wf%n_v,  &
                  tbar_ai, & ! tbar_a_j
                  wf%n_v,  &
                  one,     &
                  density, &
                  wf%n_mo)
!
   end subroutine density_ccs_mu_nu_oo_ccs
!
!
   module subroutine density_ccs_ref_mu_ov_ccs(wf, density, R_ai)
!!
!!    One electron density (EOM) reference/excited-determinant ov-term 
!!    Written by Alexander C. Paul, June 2019
!!
!!    Computes terms of the form:
!!
!!          D_pq += sum_mu < HF| e^(-T) E_pq e^T |nu > Y_mu
!!
!!    explicit term in this routine:
!!          D^R_ia = 2*R^a_i 
!!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: density
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: R_ai
!
      integer :: i, a
!
!$omp parallel do private(a, i)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            density(i, wf%n_o + a) = two*R_ai(a, i) + density(i, wf%n_o + a)
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine density_ccs_ref_mu_ov_ccs
!
!
   module subroutine density_ccs_mu_nu_vv_ccs(wf, density, tbar_ai, R_ai)
!!
!!    One electron density (EOM) excited-determinant/excited-determinant vv-term 
!!    Written by Alexander C. Paul, June 2019
!!
!!    Computes terms of the form:
!!
!!          D_pq += sum_mu,nu X_mu < mu| e^(-T) E_pq e^T |nu > Y_nu
!!
!!    explicit term in this routine:
!!          D^R_ab = sum_i R_bi tbar_ai
!!      
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: density
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: tbar_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: R_ai
!
      call dgemm('N', 'T', &
                  wf%n_v,  &
                  wf%n_v,  &
                  wf%n_o,  &
                  one,     &
                  tbar_ai, & ! tbar_a_i
                  wf%n_v,  &
                  R_ai,    & ! R_b_i
                  wf%n_v,  &
                  one,     &
                  density(wf%n_o+1, wf%n_o+1),  &
                  wf%n_mo)
!
   end subroutine density_ccs_mu_nu_vv_ccs
!
!
   module subroutine density_mu_mu_oo_ccs(wf, density, scaling_factor)
!!
!!    One electron density (EOM) excited determinant/excited determinant
!!    Written by Alexander C. Paul, June 2019
!!
!!    Should be reusable in the whole CC hierachy
!!
!!    Computes the following term containing the same excited contributions (mu,mu) 
!!    of the left and right vector
!!
!!          D_pq += sum_mu X_mu < mu| tau_mu e^(-T) E_pq e^T |HF > Y_mu
!!               += sum_mu X_mu *Y_mu < HF| e^(-T) E_pq e^T |HF >
!!
!!    This is the Hartree fock density scaled by the overlap X_mu*Y_mu
!!
!!    For the right transition density of any level in CC:
!!
!!       D^R_pq +=  sum_mu R_{k,mu} tbar_mu * < HF| e^-T E_pq e^T |HF >
!!              +=  sum_mu R_{k,mu} tbar_mu * 2 delta_pq delta_p,occ 
!!
!!       scaling_factor: overlap of tbar- and R-vector
!!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: density
!
      real(dp), intent(in) :: scaling_factor
!
      integer :: i
!
!        D^R_ii += sum_mu 2 * tbar_mu * R_mu
!
!$omp parallel do private(i)
      do i = 1, wf%n_o
!
         density(i, i) = density(i,i) + two*scaling_factor
!
      enddo
!$omp end parallel do
!
   end subroutine density_mu_mu_oo_ccs
!
!
   module subroutine density_mu_ref_ccs(wf, density_out, density_in, scaling_factor)
!!
!!    One electron density (EOM) excited determinant/reference term
!!    Written by Alexander C. Paul, June 2019
!!
!!    Should be reusable in the whole CC hierachy
!!
!!    Computes term that is a density matrix scaled by Y_ref
!!
!!          D_pq += sum_mu X_mu < mu| e^(-T) E_pq e^T |HF > Y_ref
!!
!!    For the right transition density of any level in CC:
!!
!!       D^R_pq -= sum_mu,nu R_{k,mu} tbar_mu* tbar_nu < nu| e^-T E_pq e^T |HF >
!!              -= sum_mu R_{k,mu} tbar_mu D^GS_pq
!!              -= scaling_factor D^GS_pq
!!
!!       density_out   : right_transition density
!!       density_in    : ground state density 
!!       scaling_factor: overlap of tbar- and R-vector
!!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: density_out
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in)    :: density_in
!
      real(dp), intent(in) :: scaling_factor
!
!     D^R_pq -= sum_mu tbar_mu * R_mu D^GS_pq
!
      call daxpy(wf%n_mo**2,      &
                 -scaling_factor, &
                 density_in,      &
                 1,               &
                 density_out,     &
                 1)
!
   end subroutine density_mu_ref_ccs
!
!
   module subroutine construct_left_transition_density_ccs(wf, state)
!!
!!    Construct left one-electron transition density (EOM)
!!    Written by Alexander C. Paul, June 2019
!!
!!          D^L_pq = < k| E_pq |CC >
!!
!!    where < k| is the left eigenvector of the Jacobian
!!    with amplitudes L_mu
!!
!!          < k| = sum_mu L_{k,mu} < mu| e^-T
!!
      implicit none
!
      class(ccs) :: wf
!
      integer, intent(in) :: state
!
      real(dp), dimension(:), allocatable :: L_k
!
      type(timings) :: L_TDM_timer
!
      L_TDM_timer = timings('Left transition density')
!
      call L_TDM_timer%turn_on()
!
      call mem%alloc(L_k, wf%n_es_amplitudes)
      call wf%read_excited_state(L_k, state, 'left')
!
      call zero_array(wf%left_transition_density, (wf%n_mo)**2)
!
      call wf%density_ccs_mu_ref_vo(wf%left_transition_density, L_k)
!
      call mem%dealloc(L_k, wf%n_es_amplitudes)
!
      call L_TDM_timer%turn_off()
!
   end subroutine construct_left_transition_density_ccs
!
!
   module subroutine construct_eom_etaX_ccs(wf, X, csiX, etaX)
!!
!!    Construct EOM etaX
!!    Written by Sarai D. Folkestad, May 2019
!!
!!    Constructs the EOM effective etaX vector, adding the EOM
!!    correction to etaX. 
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in)    :: csiX
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: etaX
!

      call wf%construct_etaX(X, etaX)
!
      call wf%etaX_eom_a(etaX, csiX)
!
   end subroutine construct_eom_etaX_ccs
!
!
   module subroutine construct_etaX_ccs(wf, X, etaX)
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
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: etaX
!
      call zero_array(etaX, wf%n_es_amplitudes)
!
      call wf%etaX_ccs_a1(X, etaX)
      call wf%etaX_ccs_b1(X, etaX)
!
   end subroutine construct_etaX_ccs
!
!
   module subroutine etaX_ccs_a1_ccs(wf, X, etaX_ai)
!!
!!    Construct etaX A1 
!!    Written by Josefine H. Andersen, 2019
!!
!!    Adapted by Sarai D. Folekstad, Apr 2019
!!
!!    Adds the A1 term of eta_ai^X:
!!
!!       A1 = 2X_ia
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: etaX_ai
!
      integer :: a, i
!
!$omp parallel do private(a, i)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            etaX_ai(a, i) = etaX_ai(a, i) + two*X(i, wf%n_o + a)
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine etaX_ccs_a1_ccs
!
!
   module subroutine etaX_ccs_b1_ccs(wf, X, etaX_ai)
!!
!!    Construct etaX B1
!!    Written by Josefine H. Andersen
!!
!!    Adapted by Sarai D. Folkestad, Apr 2019
!!
!!    Adds the B1 term of eta_ai^X:
!!
!!       B1 = sum_c tb_ci X_ca - sum_k tb_ak X_ik
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!      
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)   :: etaX_ai
!      
!
      real(dp), dimension(:,:), allocatable :: X_ca
      real(dp), dimension(:,:), allocatable :: X_ik
!
      integer :: i, k, a, c
!
!     :: First term  sum_c tb_ci X_ca
!
      call mem%alloc(X_ca, wf%n_v, wf%n_v)
!
!$omp parallel do private(a, c)
      do a = 1, wf%n_v
         do c = 1, wf%n_v
!
            X_ca(c, a) = X(c + wf%n_o, a + wf%n_o)
!
         enddo
      enddo
!$omp end parallel do
!
!
      call dgemm('T','N',  &
                 wf%n_v,   &
                 wf%n_o,   &
                 wf%n_v,   &
                 one,      &
                 X_ca,     &
                 wf%n_v,   &
                 wf%t1bar, & !tbar_ci
                 wf%n_v,   &
                 one,      &
                 etaX_ai,  &
                 wf%n_v)
!         
      call mem%dealloc(X_ca, wf%n_v, wf%n_v)
!      
!     :: Second term  - sum_k tb_ak X_ik
!
      call mem%alloc(X_ik, wf%n_o, wf%n_o)
!
!$omp parallel do private(k, i)
      do k = 1, wf%n_o
         do i = 1, wf%n_o
!
            X_ik(i, k) = X(i, k)
!
         enddo
      enddo
!$omp end parallel do
!      
      call dgemm('N','T',   &
                 wf%n_v,    &
                 wf%n_o,    &
                 wf%n_o,    &
                 -one,      &
                 wf%t1bar,  & ! tbar_ak
                 wf%n_v,    &
                 X_ik,      &
                 wf%n_o,    &
                 one,       &
                 etaX_ai,   &
                 wf%n_v)
!         
      call mem%dealloc(X_ik, wf%n_o, wf%n_o)
!
   end subroutine etaX_ccs_b1_ccs
!
!
   module subroutine construct_csiX_ccs(wf, X, csiX)
!!
!!    Construct csiX
!!    Written by Josefine H. Andersen, 2019
!!
!!    Adapted by Sarai D. Folkestad
!!
!!    Constructs xi^X_mu :
!!
!!       xi^X_mu = < mu| exp(-T) X exp(T)|R >
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!      
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: csiX
!
      call zero_array(csiX, wf%n_es_amplitudes)
!
      call wf%csiX_ccs_a1(X, csiX)
!
   end subroutine construct_csiX_ccs
!
!
   module subroutine csiX_ccs_a1_ccs(wf, X, csiX_ai)
!!
!!    Construct csiX A1 
!!    Written by Josefine H. Andersen, Feb 2019
!!
!!    Adds the A1 term to csiX:
!! 
!!       xi^X_ai += X_ai
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!      
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: csiX_ai
!
      integer :: a, i
!
!$omp parallel do private(a, i)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            csiX_ai(a, i) = csiX_ai(a, i) + X(wf%n_o + a, i)
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine csiX_ccs_a1_ccs
!
!
   module subroutine etaX_eom_a_ccs(wf, etaX, csiX)
!!
!!    EtaX EOM A
!!    Written by Josefine H. Andersen, Feb 2019
!!
!!    Add EOM A correction to etaX vector:
!!
!!       A:  eta^X,corr_mu += tbar_mu (xi * tbar) 
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: etaX
      real(dp), dimension(wf%n_es_amplitudes), intent(in)    :: csiX
!
      real(dp) :: X_cc
      real(dp) :: ddot
!
      X_cc = ddot(wf%n_es_amplitudes, wf%t1bar, 1, csiX, 1)
!
      call daxpy(wf%n_es_amplitudes, -X_cc, wf%t1bar, 1, etaX, 1)
!
   end subroutine etaX_eom_a_ccs
!
!
   module subroutine calculate_lr_transition_strength_ccs(wf, S, etaX, csiX, state, T_l, T_r, M)
!!
!!    Calculate LR transition strength
!!    Written by Josefine H. Andersen, February 2019
!!
!!    Given etaX and csiX, this routine calculates the left and right transition 
!!    moments T_l and T_r for the state number "state" and the transition strength 
!!    S = T_l * T_r.
!! 
!!    The left and right states L and R are read from file and made binormal by the routine.
!!
!!    Modified by Eirik F. Kjønstad, Nov 2019:
!!    Changed to calculate the LR transition strength instead of the EOM value.
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), intent(inout) :: S
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: etaX
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: csiX
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: M
!
      real(dp), intent(out) :: T_l, T_r
      integer, intent(in)   :: state
!
      real(dp), dimension(:), allocatable :: L, R
!
      real(dp) :: ddot
!
      call mem%alloc(L, wf%n_es_amplitudes)
      call mem%alloc(R, wf%n_es_amplitudes)
!
      call wf%read_excited_state(L, state, 'left')
      call wf%read_excited_state(R, state, 'right')
!
!     Left and right transition moments
!
      T_r = ddot(wf%n_es_amplitudes, etaX, 1, R, 1) + ddot(wf%n_es_amplitudes, M, 1, csiX, 1)
      T_l = ddot(wf%n_es_amplitudes, L, 1, csiX, 1)
!
!     Transition strength
!
      S  = T_l * T_r
!
      call mem%dealloc(L, wf%n_es_amplitudes)
      call mem%dealloc(R, wf%n_es_amplitudes)
!
   end subroutine calculate_lr_transition_strength_ccs
!
!
end submodule fop_ccs
