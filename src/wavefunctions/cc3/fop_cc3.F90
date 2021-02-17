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
submodule (cc3_class) fop_cc3
!
!!
!!    First order properties submodule
!!
!!    Routines for construction of the CC3 transition densities.
!!
!!    Equation-of-motion (EOM):
!!
!!    (Following Koch, H., Kobayashi, R., Sanches de Merás, A., and Jørgensen, P.,
!!    J. Chem. Phys. 100, 4393 (1994))
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
   module subroutine construct_left_transition_density_cc3(wf, state)
!!
!!    Construct left one-electron transition density
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
      class(cc3) :: wf
!
      integer, intent(in) :: state
!
      real(dp), dimension(:), allocatable :: L_k
!
      real(dp), dimension(:,:), allocatable :: L_ai, L_ia
!
      real(dp), dimension(:,:,:,:), allocatable :: L_aibj, t_aibj
      real(dp), dimension(:,:,:,:), allocatable :: L_abij, t_abij
      real(dp), dimension(:,:,:,:), allocatable :: L_ijab, t_ijab
!
      real(dp), dimension(:,:), allocatable :: density_oo, density_ov, density_vv
!
      integer :: i, j, a, b
!
      type(timings) :: L_TDM_timer
      type(timings) :: cc3_timer, ccsd_timer
      type(timings) :: cc3_ijk_timer, cc3_abc_timer
!
      L_TDM_timer    = timings('Left transition density', pl='m')
      cc3_ijk_timer  = timings('CC3 left TDM ijk batching', pl='v')
      cc3_abc_timer  = timings('CC3 left TDM abc batching', pl='v')
      cc3_timer      = timings('Total CC3 contribution left TDM', pl='n')
      ccsd_timer     = timings('Total CCSD contribution left TDM', pl='n')
!
      call L_TDM_timer%turn_on()
!
      call mem%alloc(L_k, wf%n_es_amplitudes)
      call wf%read_excited_state(L_k, state, state, 'left')
!
      call ccsd_timer%turn_on()
!
      call zero_array(wf%left_transition_density, wf%n_mo**2)
!
!     Allocate the singles part of the excitation vector
!
      call mem%alloc(L_ai, wf%n_v, wf%n_o)
!
      call dcopy(wf%n_t1, L_k, 1, L_ai, 1)
!
      call wf%density_ccs_mu_ref_vo(wf%left_transition_density, L_ai)
!
      call mem%alloc(t_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_aibj, (wf%n_v)*(wf%n_o))
!
      call wf%density_doubles_mu_ref_ov(wf%left_transition_density, L_ai, t_aibj)
!
!     Allocate and unpack doubles part of the excitation vector
!
      call mem%alloc(L_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call squareup(L_k(wf%n_t1+1 : wf%n_es_amplitudes), L_aibj, wf%n_v*wf%n_o)
!
      call mem%dealloc(L_k, wf%n_es_amplitudes)
!
      call wf%density_doubles_mu_ref_oo(wf%left_transition_density, L_aibj, t_aibj)
      call wf%density_doubles_mu_ref_vv(wf%left_transition_density, L_aibj, t_aibj)
!
      call ccsd_timer%turn_off()
!
!     :: CC3 contributions ::
!     -----------------------
!
      call cc3_timer%turn_on()
!
      call mem%alloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_1234_to_1324(t_aibj, t_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(t_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Construct covariant L2 = 1/3 (2L^ab_ij + L^ba_ij)
      call mem%alloc(L_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call sort_1234_to_1324(L_aibj, L_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call dscal(wf%n_t1**2, two*third, L_abij, 1)
      call add_2314_to_1234(third, L_aibj, L_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(L_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     :: CC3 contribution to ov- and vv-part ::   
!     ::         in batches of i,j,k         ::
!
      call mem%alloc(density_ov, wf%n_o, wf%n_v)
      call mem%alloc(density_vv, wf%n_v, wf%n_v)
      call zero_array(density_ov, wf%n_v*wf%n_o)
      call zero_array(density_vv, wf%n_v**2)
!
      call cc3_ijk_timer%turn_on()
      call wf%density_cc3_mu_ref_ijk(density_ov, density_vv,             &
                                     wf%left_excitation_energies(state), &
                                     L_ai, L_abij, t_abij,               &
                                     cvs=wf%cvs, keep_Y=.false.)
      call cc3_ijk_timer%turn_off()
!
!     Add CC3 ov and vv contributions to the density matrix
!
!$omp parallel do private(a, i)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            wf%left_transition_density(i, wf%n_o+a) = density_ov(i, a)  & 
                              + wf%left_transition_density(i, wf%n_o+a)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(density_ov, wf%n_o, wf%n_v)
!
!$omp parallel do private(a, b)
      do b = 1, wf%n_v
         do a = 1, wf%n_v
!
            wf%left_transition_density(wf%n_o+a, wf%n_o+b) = density_vv(a, b) &
                              + wf%left_transition_density(wf%n_o+a, wf%n_o+b)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(density_vv, wf%n_v, wf%n_v)
!
!     :: CC3 contribution to oo-part ::
!     ::     in batches of a,b,c     ::
!
!     Need t2, L2 in ijab ordering for the oo term
      call mem%alloc(t_ijab, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call sort_1234_to_3412(t_abij, t_ijab, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call mem%alloc(L_ijab, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call sort_1234_to_3412(L_abij, L_ijab, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(L_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call mem%alloc(L_ia, wf%n_o, wf%n_v)
      call sort_12_to_21(L_ai, L_ia, wf%n_v, wf%n_o)
      call mem%dealloc(L_ai, wf%n_v, wf%n_o)
!
      call mem%alloc(density_oo, wf%n_o, wf%n_o)
      call zero_array(density_oo, wf%n_o**2)
!
      call cc3_abc_timer%turn_on()
      call wf%density_cc3_mu_ref_abc(density_oo,                         &
                                     wf%left_excitation_energies(state), &
                                     L_ia, L_ijab, t_ijab, cvs=wf%cvs)
      call cc3_abc_timer%turn_off()
!
      call mem%dealloc(L_ia, wf%n_o, wf%n_v)
      call mem%dealloc(t_ijab, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call mem%dealloc(L_ijab, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
!$omp parallel do private(i, j)
      do j = 1, wf%n_o
         do i = 1, wf%n_o
!
            wf%left_transition_density(i, j) = density_oo(i, j)   &
                              + wf%left_transition_density(i, j)
                                                
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(density_oo, wf%n_o, wf%n_o)
      call cc3_timer%turn_off()
!
      call L_TDM_timer%turn_off()
!
   end subroutine construct_left_transition_density_cc3
!
!
   module subroutine construct_right_transition_density_cc3(wf, state)
!!
!!    Construct right one-electron transition density
!!    Written by Alexander C. Paul, June 2019
!!
!!          rho^R_pq = < Lambda| E_pq |k >
!!
!!    where |k > is the right eigenvector of the Jacobian
!!    with amplitudes R_mu
!!
!!          |k > = sum_mu (tau_mu |CC > R_{k,mu} - tbar_mu |CC > R_{k,mu}) 
!!
      use array_utilities, only: scale_diagonal
!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: state
!
      real(dp), dimension(:), allocatable :: R_k
!
      real(dp), dimension(:,:), allocatable :: R_ai
      real(dp), dimension(:,:), allocatable :: tbar_ia
!
      real(dp), dimension(:,:,:,:), allocatable :: tbar_aibj, R_aibj
      real(dp), dimension(:,:,:,:), allocatable :: tbar_abij, R_abij
      real(dp), dimension(:,:,:,:), allocatable :: tbar_ijab, R_ijab
!
      real(dp), dimension(:,:), allocatable :: density_oo, density_ov 
      real(dp), dimension(:,:), allocatable :: density_vo, density_vv
!
      real(dp) :: tbar_R_overlap, ddot
!
      integer :: i, j, a, b
!
      type(timings) :: R_TDM_timer, cc3_timer, ccsd_timer
      type(timings) :: cc3_ijk_timer, cc3_abc_timer, cc3_int_timer
!
      R_TDM_timer    = timings('Right transition density', pl='m')
      cc3_ijk_timer  = timings('CC3 right TDM ijk batching', pl='v')
      cc3_abc_timer  = timings('CC3 right TDM abc batching', pl='v')
      cc3_int_timer  = timings('CC3 right TDM contributions from intermediates', pl='v')
      cc3_timer      = timings('Total CC3 contribution right TDM', pl='n')
      ccsd_timer     = timings('Total CCSD contribution right TDM', pl='n')
!
      call R_TDM_timer%turn_on()
!
      call ccsd_timer%turn_on()
!
      call zero_array(wf%right_transition_density, (wf%n_mo)**2)
!
      call mem%alloc(R_k, wf%n_es_amplitudes)
      call wf%read_excited_state(R_k, state, state, 'right')
!
      call mem%alloc(R_ai, wf%n_v, wf%n_o)
!
      call dcopy(wf%n_t1, R_k, 1, R_ai, 1)
!
      call wf%density_ccs_mu_nu_oo(wf%right_transition_density, wf%t1bar, R_ai)
      call wf%density_ccs_ref_mu_ov(wf%right_transition_density, R_ai)
      call wf%density_ccs_mu_nu_vv(wf%right_transition_density, wf%t1bar, R_ai)
!
      tbar_R_overlap = ddot(wf%n_t1, wf%t1bar, 1, R_ai, 1)
!
      call mem%alloc(tbar_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2bar, tbar_aibj, wf%n_v*wf%n_o)
!
      call wf%density_doubles_mu_nu_ov(wf%right_transition_density, tbar_aibj, R_ai)
      call wf%density_doubles_mu_nu_vo(wf%right_transition_density, tbar_aibj, R_ai)
!
!     Allocate and unpack doubles part of the excitation vector
!
      call mem%alloc(R_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call squareup(R_k(wf%n_t1 + 1 : wf%n_es_amplitudes), R_aibj, wf%n_t1)
!
      call mem%dealloc(R_k, wf%n_es_amplitudes)
!
!     Scale the doubles vector by 1 + δ_ai,bj
!
      call scale_diagonal(two, R_aibj, wf%n_t1)
!
      tbar_R_overlap = tbar_R_overlap + half*ddot((wf%n_v*wf%n_o)**2, &
                                                   R_aibj,            &
                                                   1,                 &
                                                   tbar_aibj,         &
                                                   1)
!
!
      call wf%density_doubles_mu_ref_ov(wf%right_transition_density, &
                                        wf%t1bar, R_aibj)
!
      call wf%density_doubles_mu_ref_oo(wf%right_transition_density, &
                                        tbar_aibj, R_aibj)
      call wf%density_doubles_mu_ref_vv(wf%right_transition_density, &
                                        tbar_aibj, R_aibj)
!
      call ccsd_timer%turn_off()
!
!     :: CC3 contributions ::
!     -----------------------
!
      call cc3_timer%turn_on()
!
!     :: CC3-Contribution including intermediates of the GS-density matrix ::
!
      call mem%alloc(density_ov, wf%n_o, wf%n_v)
!
      call cc3_int_timer%turn_on()
      call wf%density_cc3_mu_nu_ov(density_ov, R_ai)
!
      call mem%alloc(density_oo, wf%n_o, wf%n_o)
      call zero_array(density_oo, wf%n_o**2)
!
      call mem%alloc(density_vv, wf%n_v, wf%n_v)
      call zero_array(density_vv, wf%n_v**2)
!
      call wf%density_cc3_mu_nu_oo_ov_vv(density_oo, density_ov, &
                                         density_vv, R_ai, R_aibj)
      call cc3_int_timer%freeze()
!
!     :: CC3 contribution in batches of i,j,k ::
!
      call mem%alloc(R_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_1234_to_1324(R_aibj, R_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(R_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Construct covariant tbar2 = 1/3 (2tbar^ab_ij + tbar^ba_ij)
      call mem%alloc(tbar_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call sort_1234_to_1324(tbar_aibj, tbar_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call dscal(wf%n_t1**2, two*third, tbar_abij, 1)
      call add_2314_to_1234(third, tbar_aibj, tbar_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(tbar_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(density_vo, wf%n_v, wf%n_o)
      call zero_array(density_vo, wf%n_v*wf%n_o)
!
      call cc3_ijk_timer%turn_on
      call wf%density_cc3_mu_nu_ijk(density_ov, density_vo, density_vv,  &
                                    wf%right_excitation_energies(state), &
                                    wf%t1bar, tbar_abij,                 &
                                    R_ai, R_abij, tbar_R_overlap)
      call cc3_ijk_timer%turn_off
!
!     :: Contribution of the vo block to the ov block ::
!
      call cc3_int_timer%turn_on
      call wf%density_cc3_mu3_nu2_ov(density_ov, density_vo)
      call cc3_int_timer%freeze()
!
!     Done with vo-block - add to density matrix
!$omp parallel do private(a, i)
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            wf%right_transition_density(wf%n_o+a, i) = density_vo(a, i)  & 
                                       + wf%right_transition_density(wf%n_o+a, i)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(density_vo, wf%n_v, wf%n_o)
!
!     Done with vv-block - add to density matrix
!$omp parallel do private(a, b)
      do b = 1, wf%n_v
         do a = 1, wf%n_v
!
            wf%right_transition_density(wf%n_o+a, wf%n_o+b) = density_vv(a, b) &
                             + wf%right_transition_density(wf%n_o+a, wf%n_o+b)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(density_vv, wf%n_v, wf%n_v)
!
!$omp parallel do private(a, i)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            wf%right_transition_density(i, wf%n_o+a) = density_ov(i, a)  & 
                                       + wf%right_transition_density(i, wf%n_o+a)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(density_ov, wf%n_o, wf%n_v)
!
!     :: CC3 contribution in batches of a,b,c ::
!
!     Need tbar2, R2 in ijab ordering for the oo term
      call mem%alloc(R_ijab, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call sort_1234_to_3412(R_abij, R_ijab, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(R_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call mem%alloc(tbar_ijab, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call sort_1234_to_3412(tbar_abij, tbar_ijab, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(tbar_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call mem%alloc(tbar_ia, wf%n_o, wf%n_v)
      call sort_12_to_21(wf%t1bar, tbar_ia, wf%n_v, wf%n_o)
!
      call cc3_abc_timer%turn_on
      call wf%density_cc3_mu_nu_abc(density_oo, wf%right_excitation_energies(state), &
                                    tbar_ia, tbar_ijab, R_ai, R_ijab)
      call cc3_abc_timer%turn_off
!
      call mem%dealloc(R_ai, wf%n_v, wf%n_o) 
      call mem%dealloc(tbar_ia, wf%n_o, wf%n_v)
      call mem%dealloc(tbar_ijab, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call mem%dealloc(R_ijab, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
!     Done with oo-block - add to density matrix
!
!$omp parallel do private(i, j)
      do j = 1, wf%n_o
         do i = 1, wf%n_o
!
            wf%right_transition_density(i, j) =  wf%right_transition_density(i, j) &
                                                + density_oo(i, j)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(density_oo, wf%n_o, wf%n_o)
!
      call cc3_timer%turn_off()
!
!     Contribution of the ground state density scaled by 
!     the right-hand side reference term (- sum_mu tbar_mu*R_mu)
!
      call wf%density_mu_mu_oo(wf%right_transition_density, tbar_R_overlap)
!
      call wf%density_mu_ref(wf%right_transition_density, &
                             wf%density,                  &
                             tbar_R_overlap)
!
      call R_TDM_timer%turn_off()
!
   end subroutine construct_right_transition_density_cc3
!
!
   module subroutine density_cc3_mu_nu_ov_cc3(wf, density_ov, R_ai)
!!
!!    One electron density (EOM) excited-determinant/excited-determinant ov-term 
!!    Written by Alexander C. Paul, August 2019
!!
!!    Computes terms of the form:
!!
!!          D_pq += sum_mu,nu X_mu < mu| e^(-T) E_pq e^T |nu > Y_nu
!!
!!    explicit term in this routine:
!!
!!       D^R_ld = -1/2 sum_{abcijk} tbar^abc_ijk(R^c_l t^abd_ijk + R^d_k t^abc_ijl)
!!              = sum_{abcijk}( -1/2 tbar^abc_ijk t^abc_ijl R^d_k 
!!                              -1/2 tbar^abc_ijk t^abd_ijk R^c_l)
!!              = sum_k D_lk R^d_k - sum_c D_cd R^c_l
!!
!!    Calculates the contribution of the oo- and vv-blocks
!!    of the GS density to the ov-block of the right TDM
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_v), intent(inout) :: density_ov
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: R_ai
!
!     rho^R_ld += sum_k D_lk R^d_k
!
      call dgemm('N', 'T',                &
                  wf%n_o,                 &
                  wf%n_v,                 &
                  wf%n_o,                 &
                  one,                    &
                  wf%GS_cc3_density_oo,   & ! D_l_k
                  wf%n_o,                 &
                  R_ai,                   & ! R_d_k
                  wf%n_v,                 &
                  zero,                   &
                  density_ov,             &
                  wf%n_o)
!
!     rho^R_ld -= sum_c R^c_l D_cd
!
      call dgemm('T', 'N',                &
                  wf%n_o,                 &
                  wf%n_v,                 &
                  wf%n_v,                 &
                  -one,                   &
                  R_ai,                   & ! R_c_l
                  wf%n_v,                 &
                  wf%GS_cc3_density_vv,   & ! D_c_d
                  wf%n_v,                 &
                  one,                    &
                  density_ov,             &
                  wf%n_o)
!
   end subroutine density_cc3_mu_nu_ov_cc3
!
!
   module subroutine density_cc3_mu_nu_oo_ov_vv_cc3(wf, density_oo, density_ov, &
                                                    density_vv, R_ai, R_aibj)
!!
!!    One electron density (EOM) excited-determinant/excited-determinant 
!!    oo-, ov- and vv-term 
!!    Written by Alexander C. Paul, August 2019
!!
!!    Computes terms of the form:
!!
!!          D_pq += sum_mu,nu X_mu < mu| e^(-T) E_pq e^T |nu > Y_nu
!!
!!    explicit terms in this routine:
!!
!!          D^R_lk -= sum{abcij}  tbar^abc_ijk t^bc_jl R^a_i
!!                 -= sum{ai}  Y_alki R^a_i
!!
!!          D^R_ld -= sum{abcijk} tbar^abc_ijk (t^bd_jk R^ac_il + t^bc_jl R^ad_ik)
!!                 -= sum{aci} Y_cadi R^ac_il + sum{aik} Y_alki R^ad_ik
!!
!!          D^R_cd += sum{abijk}  tbar^abc_ijk t^bd_jk R^a_i
!!                 += sum{ai}  Y_cadi R^a_i
!!
!!    Calculates the contribution of the Y_clik- and Y_dbck-intermediates
!!    to the oo-, ov- and vv-block of the right TDM
!!
!!    The intermediates are constructed and written to file during
!!    the construction of the GS-density and while solving for the multipliers
!!
!!          Y_ebck = sum_aij tbar^abc_ijk * t^ae_ij
!!          Y_clik = sum_abj tbar^abc_ijk * t^ab_lj
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_o), intent(inout) :: density_oo
      real(dp), dimension(wf%n_o, wf%n_v), intent(inout) :: density_ov
      real(dp), dimension(wf%n_v, wf%n_v), intent(inout) :: density_vv
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: R_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: R_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: R_jbai
!
      real(dp), dimension(:,:,:,:), allocatable :: Y_dcai
      real(dp), dimension(:,:,:,:), allocatable :: Y_alki, Y_lkai
!
      real(dp), dimension(:,:), allocatable :: rho_dc
!
!     File copied in save_tbar_intermediates
      type(direct_stream_file) :: Y_ebck_tbar
!
      type(batching_index) :: batch_i
      integer :: i_batch
      integer :: req_0 , req_i
!
!     :: Terms involving Y_dcai ::
!
      call mem%alloc(rho_dc, wf%n_v , wf%n_v)
      call zero_array(rho_dc, wf%n_v**2)
!
      call mem%alloc(R_jbai, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call sort_1234_to_4312(R_aibj, R_jbai, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      Y_ebck_tbar = direct_stream_file('Y_ebck_tbar', wf%n_v**3)
      call Y_ebck_tbar%open_('read')
!
      batch_i = batching_index(wf%n_o)
!
      req_0 = 0
      req_i = wf%n_v**3
!
      call mem%batch_setup(batch_i, req_0, req_i)
!
      do i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(i_batch)
!
         call mem%alloc(Y_dcai, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
         call Y_ebck_tbar%read_interval(Y_dcai, batch_i)
!
!        rho^R_cd += sum{ai} Y_dcai R^a_i
!
         call dgemv('N',                     &
                     wf%n_v**2,              &
                     wf%n_v*batch_i%length,  &
                     one,                    &
                     Y_dcai,                 & ! Y_dc_ai
                     wf%n_v**2,              &
                     R_ai(:,batch_i%first),  & ! R_ai
                     1,                      &
                     one,                    &
                     rho_dc,                 &
                     1)
!
!        rho^R_ld -= sum{aci} Y_dcai R^ac_il
!
         call dgemm('N', 'T',                      &
                     wf%n_o,                       &
                     wf%n_v,                       &
                     batch_i%length*wf%n_v**2,     &
                     -one,                         &
                     R_jbai(:,:,:,batch_i%first),  & ! R_l_cai
                     wf%n_o,                       &
                     Y_dcai,                       & ! Y_d_cai
                     wf%n_v,                       &
                     one,                          &
                     density_ov,                   &
                     wf%n_o)
!
         call mem%dealloc(Y_dcai, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
      enddo
!
      call mem%dealloc(R_jbai, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      call add_21_to_12(one, rho_dc, density_vv, wf%n_v, wf%n_v)
!
      call mem%dealloc(rho_dc, wf%n_v ,wf%n_v)
!
      call Y_ebck_tbar%close_()
!
!     :: Terms involving Y_alki ::
!
      call wf%Y_cmjk_tbar%open_('read')
!
      call mem%alloc(Y_alki, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call wf%Y_cmjk_tbar%read_(Y_alki, 1, wf%n_o)
!
      call wf%Y_cmjk_tbar%close_()
!
      call mem%alloc(Y_lkai, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_2314(Y_alki, Y_lkai, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%dealloc(Y_alki, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
!     rho^R_kl -= sum{ai} Y_alki R^a_i     
!
      call dgemv('N',            &
                  wf%n_o**2,     &
                  wf%n_v*wf%n_o, &
                  -one,          &
                  Y_lkai,        & ! Y_lk_ai
                  wf%n_o**2,     &
                  R_ai,          & ! R_ai
                  1,             &
                  one,           &
                  density_oo,    &
                  1)
!
!     rho^R_ld -= sum{aik} Y_alki R^ad_ik
!
      call dgemm('N', 'T',          &
                  wf%n_o,           &
                  wf%n_v,           &
                  wf%n_o**2*wf%n_v, &
                  -one,             &
                  Y_lkai,           & ! Y_l_kai
                  wf%n_o,           &
                  R_aibj,           & ! R_d_kai
                  wf%n_v,           &
                  one,              &
                  density_ov,       &
                  wf%n_o)
!
      call mem%dealloc(Y_lkai, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine density_cc3_mu_nu_oo_ov_vv_cc3
!
!
   module subroutine density_cc3_mu_nu_ijk_cc3(wf, density_ov, density_vo, &
                                               density_vv, omega,          &
                                               tbar_ai, tbar_abij,         &
                                               R_ai, R_abij, tbar_R_overlap)
!!
!!    One electron density excited-determinant/excited-determinant term 
!!    in batches of the occupied orbitals i,j,k
!!    Written by Alexander C. Paul, July 2019
!!
!!    Computes terms of the form:
!!
!!          D_pq += sum_mu,nu X_mu * < mu| e^(-T) E_pq e^T |nu > Y_nu
!!
!!    where X_mu and Y_nu are general amplitude (tbar or L)
!!
!!
!!    Construct R^abc_ijk and tbar^abc_ijk in batches of i,j,k and compute
!!    the contributions to the vo- and vv-block of the right TDM.
!!
!!    explicit terms in this routine
!!
!!    R_mu3 = (omega - eps_mu3)^-1 (< mu3| [H,R_2] |HF > 
!!                                + < mu3| [[H,R_1],T_2] |HF >)
!!    tbar_mu3 = (- eps_mu3)^-1 (tbar_mu1 < mu1| [H,tau_nu3] |R > 
!!                             + tbar_mu2 < mu2| [H,tau_nu3] |R >
!!
!!    vo-part:
!!          D^R_ck += 1/2 sum{abij} tbar^abc_ijk R^ab_ij
!!
!!    vv-part:
!!          D^R_cd += 1/2 sum_{abijk} tbar^abc_ijk R^abd_ijk
!!
!!    ov-part:
!!          D^R_kc += tbar^ab_ij (R^abc_ijk - R^abc_ikj)
!!
!!    Also construct the intermediate Z_bcjk needed for the ov-block
!!          Z_bcjk = sum{ai} tbar^abc_ijk R^a_i
!!
      use array_utilities, only: copy_and_scale
!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), intent(inout) :: tbar_R_overlap
!
      real(dp), dimension(wf%n_o, wf%n_v), intent(inout) :: density_ov
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: density_vo
      real(dp), dimension(wf%n_v, wf%n_v), intent(inout) :: density_vv
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: tbar_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: tbar_abij
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: R_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: R_abij
!
      real(dp), dimension(:,:,:,:), allocatable :: t_abij
!
      real(dp), dimension(:,:,:), allocatable :: R_abc
      real(dp), dimension(:,:,:), allocatable :: tbar_abc
      real(dp), dimension(:,:,:), allocatable :: u_abc
!
!     Help array used for sorting integrals and amplitudes
      real(dp), dimension(:,:,:,:), allocatable :: sorting
!
!     Integrals for the construction of the T3-amplitudes
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdci
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdcj
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdck
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdci_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdcj_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdck_p => null()
!
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ljci
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_lkci
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_lkcj
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_licj
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_lick
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ljck
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_ljci_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_lkci_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_lkcj_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_licj_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_lick_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_ljck_p => null()
!
!     Integrals for the construction of the triples multipliers
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_dbic
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_dbjc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_dbkc
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_dbic_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_dbjc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_dbkc_p => null()
!
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_jlic
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_klic
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_kljc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_iljc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ilkc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_jlkc
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_jlic_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_klic_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_kljc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_iljc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_ilkc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_jlkc_p => null()
!
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ibjc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ibkc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_jbkc
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_ibjc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_ibkc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_jbkc_p => null()
!
!     C1 transformed integrals
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdci_c1
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdcj_c1
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdck_c1
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdci_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdcj_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdck_c1_p => null()
!
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ljci_c1
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_lkci_c1
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_lkcj_c1
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_licj_c1
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_lick_c1
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ljck_c1
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_ljci_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_lkci_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_lkcj_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_licj_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_lick_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_ljck_c1_p => null()
!
!     Array for the intermediate
      real(dp), dimension(:,:,:,:), allocatable :: Z_bcjk, covariant_Z
      real(dp), dimension(:,:), allocatable     :: density_ai
!
      type(batching_index) :: batch_i, batch_j, batch_k
      integer :: i_batch, j_batch, k_batch
      integer :: i, j, k, i_rel, j_rel, k_rel
      integer :: req_0, req_1, req_2, req_3, req_i, req_1_eri
      integer :: req_single_batch
      real(dp) :: ddot
!
      call mem%alloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call squareup_and_sort_1234_to_1324(wf%t2, t_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      batch_i = batching_index(wf%n_o)
      batch_j = batching_index(wf%n_o)
      batch_k = batching_index(wf%n_o)
!
!     Memory for sorting array and getting the integrals
      call wf%estimate_mem_c1_integral_setup(req_0, req_1_eri)
      req_0 = req_0 + 3*wf%n_v**3 + wf%n_v*wf%n_o + (wf%n_v*wf%n_o)**2
      req_1_eri = req_1_eri + max(wf%n_v**3, wf%n_o**2*wf%n_v)
!
!     Need less memory if we don't need to batch, so we overwrite the maximum 
!     required memory in batch_setup
!
      req_single_batch = req_0 + req_1_eri*wf%n_o + 3*wf%n_v**3*wf%n_o &
                       + 3*wf%n_v*wf%n_o**3 + (wf%n_v*wf%n_o)**2
!
      req_1 = 3*wf%n_v**3
      req_i = req_1 + req_1_eri ! Mem for integral setup only needed for 1 index.
      req_2 = 6*wf%n_o*wf%n_v + wf%n_v**2
      req_3 = 0
!
      call mem%batch_setup(batch_i, batch_j, batch_k,  &
                           req_0, req_i, req_1, req_1, &
                           req_2, req_2, req_2, req_3, &
                           req_single_batch=req_single_batch)
!
      call mem%alloc(R_abc, wf%n_v, wf%n_v, wf%n_v)
      call mem%alloc(u_abc, wf%n_v, wf%n_v, wf%n_v)
      call mem%alloc(tbar_abc, wf%n_v, wf%n_v, wf%n_v)
!
      call mem%alloc(density_ai, wf%n_v, wf%n_o)
      call zero_array(density_ai, wf%n_t1)
!
      call mem%alloc(Z_bcjk, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call zero_array(Z_bcjk, wf%n_t1**2)
!
      if (batch_i%num_batches .eq. 1) then ! no batching
!
         call mem%alloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%alloc(g_bdci_c1, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%alloc(g_dbic, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
!
         call mem%alloc(g_ljci, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         call mem%alloc(g_ljci_c1, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
         call mem%alloc(g_jlic, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
         call mem%alloc(g_ibjc, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
         if (wf%n_o .le. wf%n_v) then
            call mem%alloc(sorting, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         else
            call mem%alloc(sorting, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         end if
!
      else ! batching
!
         call mem%alloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%alloc(g_bdcj, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%alloc(g_bdck, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
!
         call mem%alloc(g_ljci, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_lkci, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_lkcj, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_licj, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_lick, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_ljck, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
!
         call mem%alloc(g_bdci_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%alloc(g_bdcj_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%alloc(g_bdck_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
!
         call mem%alloc(g_ljci_c1, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_lkci_c1, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_lkcj_c1, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_licj_c1, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_lick_c1, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_ljck_c1, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
!
         call mem%alloc(g_dbic, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%alloc(g_dbjc, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%alloc(g_dbkc, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
!
         call mem%alloc(g_jlic, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_klic, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_kljc, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_iljc, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_ilkc, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_jlkc, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
!
         call mem%alloc(g_ibjc, wf%n_v, wf%n_v, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_ibkc, wf%n_v, wf%n_v, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_jbkc, wf%n_v, wf%n_v, batch_i%max_length, batch_i%max_length)
!
         if (wf%n_o .le. wf%n_v) then
            call mem%alloc(sorting, wf%n_v, wf%n_v, wf%n_v,batch_i%max_length)
         else
            call mem%alloc(sorting, wf%n_v, wf%n_o, wf%n_o,batch_i%max_length)
         end if
!
      endif

!
!     Loop over the batches in i,j,k
!     Read integrals and assign pointers
!     Without pointers we'll have to use three times as much
!     memory for the non-batching case
!
      do i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(i_batch)
!
         call wf%setup_vvvo(g_bdci, g_bdci_p, sorting, batch_i)
         call wf%setup_vvvo(g_bdci_c1, g_bdci_c1_p, sorting, batch_i, c_ai=R_ai)
!
         call wf%setup_vvov(g_dbic, g_dbic_p, sorting, batch_i, left=.true.)
!
         do j_batch = 1, i_batch
!
            call batch_j%determine_limits(j_batch)
!
            call wf%setup_oovo(g_ljci, g_ljci_p, sorting, batch_j, batch_i)
            call wf%setup_oovo(g_ljci_c1, g_ljci_c1_p, sorting, batch_j, batch_i, c_ai=R_ai)
!
            call wf%setup_ooov(g_jlic, g_jlic_p, sorting, batch_j, batch_i)
!
            call wf%setup_ovov(g_ibjc, g_ibjc_p, sorting, batch_i, batch_j)
!
            if (j_batch .ne. i_batch) then
!
               call wf%setup_vvvo(g_bdcj, g_bdcj_p, sorting, batch_j)
               call wf%setup_vvvo(g_bdcj_c1, g_bdcj_c1_p, sorting, batch_j, c_ai=R_ai)
!
               call wf%setup_vvov(g_dbjc, g_dbjc_p, sorting, batch_j, left=.true.)
!
               call wf%setup_oovo(g_licj, g_licj_p, sorting, batch_i, batch_j)
               call wf%setup_oovo(g_licj_c1, g_licj_c1_p, sorting, batch_i, batch_j, c_ai=R_ai)
!
               call wf%setup_ooov(g_iljc, g_iljc_p, sorting, batch_i, batch_j)
!
            else
!
               call wf%point_vvvo(g_bdcj_p, g_bdci, batch_j%length)
               call wf%point_vvvo(g_bdcj_c1_p, g_bdci_c1, batch_j%length)
!
               call wf%point_vvvo(g_dbjc_p, g_dbic, batch_j%length)
!
               call wf%point_vooo(g_licj_p, g_ljci, batch_i%length, batch_j%length)
               call wf%point_vooo(g_licj_c1_p, g_ljci_c1, batch_i%length, batch_j%length)
!
               call wf%point_vooo(g_iljc_p, g_jlic, batch_i%length, batch_j%length)
!
            endif
!
            do k_batch = 1, j_batch
!
               call batch_k%determine_limits(k_batch)
!
               if (k_batch .ne. j_batch) then ! k_batch != j_batch, k_batch != i_batch
!
                  call wf%setup_vvvo(g_bdck, g_bdck_p, sorting, batch_k)
                  call wf%setup_vvvo(g_bdck_c1, g_bdck_c1_p, sorting, batch_k, c_ai=R_ai)
!
                  call wf%setup_vvov(g_dbkc, g_dbkc_p, sorting, batch_k, left=.true.)
!
                  call wf%setup_oovo(g_lick, g_lick_p, sorting, batch_i, batch_k)
                  call wf%setup_oovo(g_ljck, g_ljck_p, sorting, batch_j, batch_k)
                  call wf%setup_oovo(g_lkci, g_lkci_p, sorting, batch_k, batch_i)
                  call wf%setup_oovo(g_lkcj, g_lkcj_p, sorting, batch_k, batch_j)
!
                  call wf%setup_oovo(g_lick_c1, g_lick_c1_p, sorting, batch_i, batch_k, c_ai=R_ai)
                  call wf%setup_oovo(g_ljck_c1, g_ljck_c1_p, sorting, batch_j, batch_k, c_ai=R_ai)
                  call wf%setup_oovo(g_lkci_c1, g_lkci_c1_p, sorting, batch_k, batch_i, c_ai=R_ai)
                  call wf%setup_oovo(g_lkcj_c1, g_lkcj_c1_p, sorting, batch_k, batch_j, c_ai=R_ai)
!
                  call wf%setup_ooov(g_ilkc, g_ilkc_p, sorting, batch_i, batch_k)
                  call wf%setup_ooov(g_jlkc, g_jlkc_p, sorting, batch_j, batch_k)
                  call wf%setup_ooov(g_klic, g_klic_p, sorting, batch_k, batch_i)
                  call wf%setup_ooov(g_kljc, g_kljc_p, sorting, batch_k, batch_j)
!
                  call wf%setup_ovov(g_ibkc, g_ibkc_p, sorting, batch_i, batch_k)
                  call wf%setup_ovov(g_jbkc, g_jbkc_p, sorting, batch_j, batch_k)
!
               else if (k_batch .eq. i_batch) then ! k_batch == j_batch == i_batch
!
                  call wf%point_vvvo(g_bdck_p, g_bdci, batch_k%length)
                  call wf%point_vvvo(g_bdck_c1_p, g_bdci_c1, batch_k%length)
!
                  call wf%point_vvvo(g_dbkc_p, g_dbic, batch_k%length)
!
                  call wf%point_vooo(g_lick_p, g_ljci, batch_i%length, batch_k%length)
                  call wf%point_vooo(g_ljck_p, g_ljci, batch_j%length, batch_k%length)
                  call wf%point_vooo(g_lkci_p, g_ljci, batch_k%length, batch_i%length)
                  call wf%point_vooo(g_lkcj_p, g_ljci, batch_k%length, batch_j%length)
!
                  call wf%point_vooo(g_lick_c1_p, g_ljci_c1, batch_i%length, batch_k%length)
                  call wf%point_vooo(g_ljck_c1_p, g_ljci_c1, batch_j%length, batch_k%length)
                  call wf%point_vooo(g_lkci_c1_p, g_ljci_c1, batch_k%length, batch_i%length)
                  call wf%point_vooo(g_lkcj_c1_p, g_ljci_c1, batch_k%length, batch_j%length)
!
                  call wf%point_vooo(g_ilkc_p, g_jlic, batch_i%length, batch_k%length)
                  call wf%point_vooo(g_jlkc_p, g_jlic, batch_j%length, batch_k%length)
                  call wf%point_vooo(g_klic_p, g_jlic, batch_k%length, batch_i%length)
                  call wf%point_vooo(g_kljc_p, g_jlic, batch_k%length, batch_j%length)
!
                  call wf%point_vvoo(g_ibkc_p, g_ibjc, batch_i%length, batch_k%length)
                  call wf%point_vvoo(g_jbkc_p, g_ibjc, batch_j%length, batch_k%length)
!
               else ! k_batch == j_batch != i_batch
!
                  call wf%point_vvvo(g_bdck_p, g_bdcj, batch_k%length)
                  call wf%point_vvvo(g_bdck_c1_p, g_bdcj_c1, batch_k%length)
!
                  call wf%point_vvvo(g_dbkc_p, g_dbjc, batch_k%length)
!
                  call wf%setup_oovo(g_lkcj, g_lkcj_p, sorting, batch_k, batch_j)
                  call wf%point_vooo(g_lick_p, g_licj, batch_i%length, batch_k%length)
                  call wf%point_vooo(g_ljck_p, g_lkcj, batch_j%length, batch_k%length)
                  call wf%point_vooo(g_lkci_p, g_ljci, batch_k%length, batch_i%length)
!
                  call wf%setup_oovo(g_lkcj_c1, g_lkcj_c1_p, sorting, batch_k, batch_j, c_ai=R_ai)
                  call wf%point_vooo(g_lick_c1_p, g_licj_c1, batch_i%length, batch_k%length)
                  call wf%point_vooo(g_ljck_c1_p, g_lkcj_c1, batch_j%length, batch_k%length)
                  call wf%point_vooo(g_lkci_c1_p, g_ljci_c1, batch_k%length, batch_i%length)
!
                  call wf%setup_ooov(g_kljc, g_kljc_p, sorting, batch_k, batch_j)
                  call wf%point_vooo(g_ilkc_p, g_iljc, batch_i%length, batch_k%length)
                  call wf%point_vooo(g_jlkc_p, g_kljc, batch_j%length, batch_k%length)
                  call wf%point_vooo(g_klic_p, g_jlic, batch_k%length, batch_i%length)
!
                  call wf%setup_ovov(g_jbkc, g_jbkc_p, sorting, batch_j, batch_k)
                  call wf%point_vvoo(g_ibkc_p, g_ibjc, batch_i%length, batch_k%length)
!
               endif
!
               do i = batch_i%first, batch_i%last
!
                  i_rel = i - batch_i%first + 1
!
                  do j = batch_j%first, min(batch_j%last, i)
!
                     j_rel = j - batch_j%first + 1
!
                     do k = batch_k%first, min(batch_k%last, j)
!
                        if (i .eq. j .and. i .eq. k) then
                           cycle
                        end if
!
                        k_rel = k - batch_k%first + 1
!
!                       Check if at least one index i,j,k is a core orbital
                        if(wf%cvs) then
!
                           if(.not. (any(wf%core_MOs .eq. i) &
                              .or.   any(wf%core_MOs .eq. j) &
                              .or.   any(wf%core_MOs .eq. k))) cycle
!
                        end if
!
!                       Also terms that only include the triples multipliers 
!                       can be cycled if we use CVS
!                       Because they are contracted with R^a_i or R^ab_ij
!                       which have to contain at least one core orbital.
!
!                       Construct covariant W_abc for given i,j,k
!                       L_abc is obtained by dividing the linear combination
!                       4W_abc - 2W_bac - 2W_cba - 2W_acb + W_bca + W_cab
!                       by omega - eps^abc_ijk
!
                        call wf%construct_W(i, j, k, sorting, tbar_abc,  &
                                            tbar_abij,                  &
                                            g_dbic_p(:,:,:,i_rel),      &
                                            g_dbjc_p(:,:,:,j_rel),      &
                                            g_dbkc_p(:,:,:,k_rel),      &
                                            g_jlic_p(:,:,j_rel,i_rel),  &
                                            g_klic_p(:,:,k_rel,i_rel),  &
                                            g_kljc_p(:,:,k_rel,j_rel),  &
                                            g_iljc_p(:,:,i_rel,j_rel),  &
                                            g_ilkc_p(:,:,i_rel,k_rel),  &
                                            g_jlkc_p(:,:,j_rel,k_rel))
!
                        call wf%outer_product_terms_l3(i, j, k, tbar_ai, tbar_abij, &
                                                       tbar_abc, wf%fock_ia,        &
                                                       g_ibjc_p(:,:,i_rel,j_rel),   &
                                                       g_ibkc_p(:,:,i_rel,k_rel),   &
                                                       g_jbkc_p(:,:,j_rel,k_rel))
!
                        call construct_contravariant_t3(tbar_abc, u_abc, wf%n_v)
!
                        call wf%divide_by_orbital_differences(i, j, k, tbar_abc, zero)
!
                        call wf%density_cc3_mu_nu_vo(i, j, k, tbar_abc, u_abc,   &
                                                               density_vo, R_abij)
!
!                       Construct intermediate used for the ov-block:
!                          Z_bcjk = sum_ai tbar^abc_ijk R^a_i
!
                        call wf%construct_Z_intermediate(i, j, k, tbar_abc, u_abc,  &
                                                         Z_bcjk, R_ai)
!
!                       Construct R^{abc}_{ijk} for given i, j, k
!                       Using c1-transformed integrals the terms have the same form 
!                       as the omega terms (where t_abc = R_abc)
!
                        call wf%construct_V(i, j, k, u_abc, R_abc,        &
                                            t_abij, R_abij,               &
                                            g_bdci_p(:,:,:,i_rel),        &
                                            g_bdcj_p(:,:,:,j_rel),        &
                                            g_bdck_p(:,:,:,k_rel),        &
                                            g_bdci_c1_p(:,:,:,i_rel),     &
                                            g_bdcj_c1_p(:,:,:,j_rel),     &
                                            g_bdck_c1_p(:,:,:,k_rel),     &
                                            g_ljci_p(:,:,j_rel,i_rel),    &
                                            g_lkci_p(:,:,k_rel,i_rel),    &
                                            g_lkcj_p(:,:,k_rel,j_rel),    &
                                            g_licj_p(:,:,i_rel,j_rel),    &
                                            g_lick_p(:,:,i_rel,k_rel),    &
                                            g_ljck_p(:,:,j_rel,k_rel),    &
                                            g_ljci_c1_p(:,:,j_rel,i_rel), &
                                            g_lkci_c1_p(:,:,k_rel,i_rel), &
                                            g_lkcj_c1_p(:,:,k_rel,j_rel), &
                                            g_licj_c1_p(:,:,i_rel,j_rel), &
                                            g_lick_c1_p(:,:,i_rel,k_rel), &
                                            g_ljck_c1_p(:,:,j_rel,k_rel))
!
                        call wf%divide_by_orbital_differences(i, j, k, R_abc, omega)
!
                        call wf%density_cc3_mu_ref_vv(i, j, k, density_vv, R_abc, &
                                                      u_abc, tbar_abc, sorting)
!
!                       Overlap:
!                          tbar_R_overlap = sum_{ai >= bj >= ck} L^abc_ijk R^abc_ijk
!                                         = 1/6 sum_abcijk L^abc_ijk R^abc_ijk
!
!                       Due to the restrictions on the loops the factor of 1/6 
!                       vanishes, but we need to account for double counting 
!                       if 2 indices are the same e.g.:
!                             tbar^abc_ijk*R^abc_ijk = tbar^abc_jik*R^abc_jik
!
                        if (i .ne. j .and. j .ne. k) then
                           tbar_R_overlap = tbar_R_overlap &
                                          + ddot(wf%n_v**3, tbar_abc, 1 , R_abc, 1)
                        else
                           tbar_R_overlap = tbar_R_overlap &
                                          + half*ddot(wf%n_v**3, tbar_abc, 1 , R_abc, 1)
                        end if
!
!                       Need the contravariant R to use construct_x_ai
                        call construct_contravariant_t3(R_abc, u_abc, wf%n_v)
!
                        call wf%construct_x_ai_intermediate(i, j, k, R_abc, u_abc, &
                                                            tbar_abij, density_ai)
!
                     enddo ! loop over k
                  enddo ! loop over j
               enddo ! loop over i
            enddo ! batch_k
         enddo ! batch_j
      enddo ! batch_i
!
      call mem%dealloc(R_abc, wf%n_v, wf%n_v, wf%n_v)
      call mem%dealloc(u_abc, wf%n_v, wf%n_v, wf%n_v)
      call mem%dealloc(tbar_abc, wf%n_v, wf%n_v, wf%n_v)
!
      if (batch_i%num_batches .eq. 1) then ! no batching
!
         call mem%dealloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%dealloc(g_dbic, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%dealloc(g_ljci, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         call mem%dealloc(g_jlic, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         call mem%dealloc(g_ibjc, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
         call mem%dealloc(g_bdci_c1, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%dealloc(g_ljci_c1, wf%n_v , wf%n_o, wf%n_o, wf%n_o)
!
         if (wf%n_o .le. wf%n_v) then
            call mem%dealloc(sorting, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         else
            call mem%dealloc(sorting, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         end if
!
      else ! batching
!
         call mem%dealloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%dealloc(g_bdcj, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%dealloc(g_bdck, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
!
         call mem%dealloc(g_dbic, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%dealloc(g_dbjc, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%dealloc(g_dbkc, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
!
         call mem%dealloc(g_ljci, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_lkci, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_lkcj, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_licj, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_lick, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_ljck, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
!
         call mem%dealloc(g_jlic, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_klic, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_kljc, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_iljc, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_ilkc, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_jlkc, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
!
         call mem%dealloc(g_ibjc, wf%n_v, wf%n_v, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_ibkc, wf%n_v, wf%n_v, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_jbkc, wf%n_v, wf%n_v, batch_i%max_length, batch_i%max_length)
!
         call mem%dealloc(g_bdci_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%dealloc(g_bdcj_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%dealloc(g_bdck_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
!
         call mem%dealloc(g_ljci_c1, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_lkci_c1, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_lkcj_c1, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_licj_c1, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_lick_c1, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_ljck_c1, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
!
         if (wf%n_o .le. wf%n_v) then
            call mem%dealloc(sorting, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         else
            call mem%dealloc(sorting, wf%n_v, wf%n_o, wf%n_o, batch_i%max_length)
         end if
!
      endif
!
!     Add contribution of to the ov-block
      call add_21_to_12(one, density_ai, density_ov, wf%n_o, wf%n_v)
!
      call mem%dealloc(density_ai, wf%n_v, wf%n_o)
!
!     Contribution of the Z_bcjk intermediate
!
!     Construct covariant_Z = 1/3 (2 Z_bcjk + Z_bckj)
!
      call mem%alloc(covariant_Z, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call copy_and_scale(two*third, Z_bcjk, covariant_Z, wf%n_t1**2)
      call add_2134_to_1234(third, Z_bcjk, covariant_Z, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call mem%dealloc(Z_bcjk, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call wf%density_cc3_mu_nu_ov_Z_term(density_ov, covariant_Z, t_abij)
!
      call mem%dealloc(covariant_Z, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
   end subroutine density_cc3_mu_nu_ijk_cc3
!
!
   module subroutine density_cc3_mu_nu_vo_cc3(wf, i, j, k, tbar_abc,   &
                                              v_abc, density_vo, R_abij)
!!
!!    One electron density (EOM) excited-determinant/excited-determinant vo-term 
!!    Written by Alexander C. Paul, August 2019
!!
!!    Computes terms of the form:
!!
!!          D_pq += sum_mu,nu X_mu < mu| e^(-T) E_pq e^T |nu > Y_nu
!!
!!    explicit term in this routine:
!!
!!          D^R_ck += sum_{abij} tbar^abc_ijk R^ab_ij
!!
!!    All permutations for i,j,k have to be considered 
!!    due to the restrictions in the i,j,k loops
!!
!!    based on construct_x_ai_intermediate_cc3
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)         :: tbar_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)        :: v_abc
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)              :: density_vo
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: R_abij
!
!     rho^R_ck += 1/2 sum_{abij} tbar^abc_ijk R^ab_ij
      call wf%omega1_cc3_permutation(k, tbar_abc, R_abij(:,:,i,j), density_vo, half)
!
!     rho^R_ci += 1/2 sum_{abjk} tbar^cab_ijk R^ab_jk
      call wf%construct_X_vo_permutation(i, tbar_abc, R_abij(:,:,j,k), density_vo, half)
!
!     bca -> cab
!
      call sort_123_to_231(tbar_abc, v_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     rho^R_cj += 1/2 sum_{abik} tbar^bca_ijk R^ab_ki
      call wf%construct_X_vo_permutation(j, v_abc, R_abij(:,:,k,i), density_vo, half)
!
      if (k .ne. j .and. j .ne. i) then
!
!        bac -> abc
!        acb -> cab
!
         call sort_123_to_213(tbar_abc, v_abc, wf%n_v, wf%n_v, wf%n_v)
!
!        rho^R_ck += 1/2 sum_{abij} tbar^bac_ijk R^ab_ji
         call wf%omega1_cc3_permutation(k, v_abc, R_abij(:,:,j,i), density_vo, half)
!
!        rho^R_cj += 1/2 sum_{abik} tbar^acb_ijk R^ab_ik
         call wf%construct_X_vo_permutation(j, v_abc, R_abij(:,:,i,k), density_vo, half)
!
!        cba -> cab
!
         call sort_123_to_132(tbar_abc, v_abc, wf%n_v, wf%n_v, wf%n_v)
!
!        rho^R_ci += 1/2 sum_{abik} tbar^cba_ijk R^ab_kj
         call wf%construct_X_vo_permutation(i, v_abc, R_abij(:,:,k,j), density_vo, half)
!
      end if
!
   end subroutine density_cc3_mu_nu_vo_cc3
!
!
   module subroutine construct_Z_intermediate_cc3(wf, i, j, k, tbar_abc, &
                                                  v_abc, Z_bcjk, R_ai)
!!
!!    Constructs Z-intermediate 
!!    Written by Alexander C. Paul, August 2019
!!
!!    based on jacobian_cc3_b2_fock_cc3
!!    Written by Alexander C. Paul and Rolf H. Myhre, Feb 2019
!!
!!    needed for the ov-block of the right TDM
!!
!!       Z_bcjk += sum_{ai} tbar^abc_ijk R^a_i
!!
!!    All permutations for i,j,k have to be considered 
!!    due to the restrictions in the i,j,k loops
!!      
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)           :: tbar_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: v_abc
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(out)  :: Z_bcjk
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)   :: R_ai
!
!     Z_bcjk += sum_{ai} tbar^abc_ijk R^a_i
      call wf%omega2_fock_cc3_permutation(j, k, tbar_abc, R_ai(:,i), Z_bcjk, one)
!
!     Z_bcij += sum_{ak} tbar^bca_ijk R^a_k
      call wf%rho2_fock_cc3_permutation(i, j, tbar_abc, R_ai(:,k), Z_bcjk, one)
!
!     cab -> abc
      call sort_123_to_231(tbar_abc, v_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     Z_bcki += sum_{aj} tbar^cab_ijk R^a_j
      call wf%omega2_fock_cc3_permutation(k, i, v_abc, R_ai(:,j), Z_bcjk, one)
!
      if (k .ne. j .and. j .ne. i) then
!
!        bac -> abc
!        cba -> bca
         call sort_123_to_213(tbar_abc, v_abc, wf%n_v, wf%n_v, wf%n_v)
!
!        Z_bcik += sum_{aj} tbar^bac_ijk R^a_j
         call wf%omega2_fock_cc3_permutation(i, k, v_abc, R_ai(:,j), Z_bcjk, one)
!
!        Z_bcji += sum_{ak} tbar^cba_ijk R^a_k
         call wf%rho2_fock_cc3_permutation(j, i, v_abc, R_ai(:,k), Z_bcjk, one)
!
!        acb -> bca
         call sort_123_to_321(tbar_abc, v_abc, wf%n_v, wf%n_v, wf%n_v)
!
!        Z_bckj += sum_{ak} tbar^acb_ijk R^a_i
         call wf%rho2_fock_cc3_permutation(k, j, v_abc, R_ai(:,i), Z_bcjk, one)
!
      end if
!
   end subroutine construct_Z_intermediate_cc3
!
!
   module subroutine density_cc3_mu3_nu2_ov_cc3(wf, density_ov, density_vo)
!!
!!    One electron density (EOM) triples-determinant/doubles-determinant ov-term 
!!    Written by Alexander C. Paul, August 2019
!!
!!    Computes terms of the form:
!!
!!          D_pq += sum_mu,nu X_mu < mu| e^(-T) E_pq e^T |nu > Y_nu
!!
!!    explicit term in this routine:
!!
!!          D^R_ld += sum{abcijk} 1/2 tbar^abc_ijk R^ab_ij(2t^cd_kl-t^cd_lk)
!!                 += sum{ck} D^R_ck (2t^cd_kl-t^cd_lk)
!!
!!    Calculates the contribution of the vo-block of the right TDM to
!!    the ov-block of the right TDM
!!
      use array_utilities, only: copy_and_scale
!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_v), intent(inout) :: density_ov
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: density_vo
!
      real(dp), dimension(:,:,:,:), allocatable :: t_ckld, u_ckld
!
      call mem%alloc(t_ckld, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call squareup_and_sort_1234_to_1243(wf%t2, t_ckld, wf%n_v, wf%n_o,  wf%n_v, wf%n_o)
!
      call mem%alloc(u_ckld, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call copy_and_scale(two, t_ckld, u_ckld, wf%n_t1**2)
      call add_1324_to_1234(-one, t_ckld, u_ckld, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call mem%dealloc(t_ckld, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemv('T',         &
                  wf%n_t1,    &
                  wf%n_t1,    &
                  one,        &
                  u_ckld,     & ! u_ck_ld
                  wf%n_t1,    &
                  density_vo, & ! rho_ck
                  1,          &
                  one,        &
                  density_ov, &
                  1)
!
      call mem%dealloc(u_ckld, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
   end subroutine density_cc3_mu3_nu2_ov_cc3
!
!
   module subroutine density_cc3_mu_nu_abc_cc3(wf, density_oo, omega, &
                                               tbar_ia, tbar_ijab, R_ai, R_ijab)
!!
!!    One electron density excited-determinant/excited-determinant term 
!!    in batches of the virtual orbitals a,b,c
!!    Written by Alexander C. Paul, August 2019
!!
!!    Computes terms of the form:
!!
!!          D_pq += sum_mu,nu X_mu * < mu| e^(-T) E_pq e^T |nu > Y_nu
!!
!!    where X_mu and Y_nu are general amplitude (tbar or L)
!!
!!    Construct R^abc_ijk and tbar^abc_ijk in batches of a,b,c
!!    and compute the contribution to the oo-block of the right TDM.
!!
!!    explicit terms in this routine
!!
!!    R_mu3 = (omega - eps_mu3)^-1 (< mu3| [H,R_2] |HF > 
!!                                + < mu3| [[H,R_1],T_2] |HF >)
!!
!!    tbar_mu3 = (- eps_mu3)^-1 (tbar_mu1 < mu1| [H,tau_nu3] |R > 
!!                             + tbar_mu2 < mu2| [H,tau_nu3] |R >
!!
!!    oo-part:
!!          D^R_kl -= 1/2 sum_{abcij} tbar^abc_ijl R^abc_ijk
!!
      use omp_lib
! 
      implicit none
!
      class(cc3) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_o, wf%n_o), intent(inout) :: density_oo
!
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: tbar_ia
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_v, wf%n_v), intent(in) :: tbar_ijab
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: R_ai
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_v, wf%n_v), intent(in) :: R_ijab
!
      real(dp), dimension(:,:,:,:), allocatable :: t_ijab
!
      real(dp), dimension(:,:,:,:), allocatable :: R_ijk
      real(dp), dimension(:,:,:,:), allocatable :: tbar_ijk
      real(dp), dimension(:,:,:,:), allocatable :: u_ijk
      real(dp), dimension(:,:,:,:), allocatable :: v_ijk
!
!     Help array used for sorting integrals and amplitudes
      real(dp), dimension(:,:,:,:), allocatable :: sorting
!
!     Integrals for the construction of the R3-amplitudes
      real(dp), dimension(:,:,:,:), allocatable, target :: g_ljak
      real(dp), dimension(:,:,:,:), allocatable, target :: g_ljbk
      real(dp), dimension(:,:,:,:), allocatable, target :: g_ljck
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_ljak_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_ljbk_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_ljck_p => null()
!
      real(dp), dimension(:,:,:,:), allocatable, target :: g_bdak
      real(dp), dimension(:,:,:,:), allocatable, target :: g_cdak
      real(dp), dimension(:,:,:,:), allocatable, target :: g_cdbk
      real(dp), dimension(:,:,:,:), allocatable, target :: g_adbk
      real(dp), dimension(:,:,:,:), allocatable, target :: g_adck
      real(dp), dimension(:,:,:,:), allocatable, target :: g_bdck
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_bdak_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_cdak_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_cdbk_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_adbk_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_adck_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_bdck_p => null()
!
!     c1-transformed integrals 
      real(dp), dimension(:,:,:,:), allocatable, target :: g_ljak_c1
      real(dp), dimension(:,:,:,:), allocatable, target :: g_ljbk_c1
      real(dp), dimension(:,:,:,:), allocatable, target :: g_ljck_c1
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_ljak_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_ljbk_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_ljck_c1_p => null()
!
      real(dp), dimension(:,:,:,:), allocatable, target :: g_bdak_c1
      real(dp), dimension(:,:,:,:), allocatable, target :: g_cdak_c1
      real(dp), dimension(:,:,:,:), allocatable, target :: g_cdbk_c1
      real(dp), dimension(:,:,:,:), allocatable, target :: g_adbk_c1
      real(dp), dimension(:,:,:,:), allocatable, target :: g_adck_c1
      real(dp), dimension(:,:,:,:), allocatable, target :: g_bdck_c1
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_bdak_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_cdak_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_cdbk_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_adbk_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_adck_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_bdck_c1_p => null()
!
!     Integrals for the construction of the multipliers
      real(dp), dimension(:,:,:,:), allocatable, target :: g_jlka
      real(dp), dimension(:,:,:,:), allocatable, target :: g_jlkb
      real(dp), dimension(:,:,:,:), allocatable, target :: g_jlkc
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_jlka_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_jlkb_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_jlkc_p => null()
!
      real(dp), dimension(:,:,:,:), allocatable, target :: g_dbka
      real(dp), dimension(:,:,:,:), allocatable, target :: g_dcka
      real(dp), dimension(:,:,:,:), allocatable, target :: g_dckb
      real(dp), dimension(:,:,:,:), allocatable, target :: g_dakb
      real(dp), dimension(:,:,:,:), allocatable, target :: g_dakc
      real(dp), dimension(:,:,:,:), allocatable, target :: g_dbkc
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_dbka_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_dcka_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_dckb_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_dakb_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_dakc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_dbkc_p => null()
!
      real(dp), dimension(:,:,:,:), allocatable, target :: g_jakb
      real(dp), dimension(:,:,:,:), allocatable, target :: g_jakc
      real(dp), dimension(:,:,:,:), allocatable, target :: g_jbkc
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_jakb_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_jakc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_jbkc_p => null()
!
!     Temporary array for each thread
      real(dp), dimension(:,:,:), allocatable :: density_oo_thread
!
      integer              :: a, b, c, a_rel, b_rel, c_rel
      type(batching_index) :: batch_a, batch_b, batch_c
      integer              :: a_batch, b_batch, c_batch
      integer              :: req_0, req_1, req_2, req_3, req_1_eri, req_a
      integer              :: n_threads, thread_n
      integer              :: req_single_batch
!
      thread_n  = 1
      n_threads = 1
!
!$    n_threads = omp_get_max_threads()
      call mem%alloc(density_oo_thread, wf%n_o, wf%n_o, n_threads)
      call zero_array(density_oo_thread, wf%n_o**2*n_threads)
!
!     Set up arrays for amplitudes
      call mem%alloc(t_ijab, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call squareup_and_sort_1234_to_2413(wf%t2, t_ijab, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      batch_a = batching_index(wf%n_v)
      batch_b = batching_index(wf%n_v)
      batch_c = batching_index(wf%n_v)
!
!     Memory for sorting array and getting the integrals
      call wf%estimate_mem_integral_setup_abc(req_0, req_1_eri)
      req_1_eri = req_1_eri + max(wf%n_v**2*wf%n_o, wf%n_o**2*wf%n_v)
      req_0 = req_0 + 4*wf%n_o**3*n_threads
!
!     Need less memory if we don't need to batch, so we overwrite the maximum 
!     required memory in batch_setup
!
      req_single_batch = req_0 + req_1_eri*wf%n_v + 3*wf%n_v**3*wf%n_o &
                       + 3*wf%n_v*wf%n_o**3 + (wf%n_v*wf%n_o)**2
!
      req_1 = 3*wf%n_o**3
      req_a = req_1 + req_1_eri ! Mem for integral setup only needed for 1 index.
      req_2 = 6*wf%n_o*wf%n_v + wf%n_o**2
      req_3 = 0
!
      call mem%batch_setup(batch_a, batch_b, batch_c,  &
                           req_0, req_a, req_1, req_1, &
                           req_2, req_2, req_2, req_3, &
                           req_single_batch=req_single_batch)
!
      call mem%alloc(R_ijk, wf%n_o, wf%n_o, wf%n_o, n_threads)
      call mem%alloc(u_ijk, wf%n_o, wf%n_o, wf%n_o, n_threads)
      call mem%alloc(v_ijk, wf%n_o, wf%n_o, wf%n_o, n_threads)
      call mem%alloc(tbar_ijk, wf%n_o, wf%n_o, wf%n_o, n_threads)
!
      if (batch_a%num_batches .eq. 1) then ! no batching
!
         call mem%alloc(g_ljak, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
         call mem%alloc(g_bdak, wf%n_v, wf%n_o, wf%n_v, wf%n_v)
!
         call mem%alloc(g_ljak_c1, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
         call mem%alloc(g_bdak_c1, wf%n_v, wf%n_o, wf%n_v, wf%n_v)
!
         call mem%alloc(g_jlka, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
         call mem%alloc(g_dbka, wf%n_o, wf%n_v, wf%n_v, wf%n_v)
         call mem%alloc(g_jakb, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
         if (wf%n_o .lt. wf%n_v) then
            call mem%alloc(sorting, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         else
            call mem%alloc(sorting, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         end if
!
      else ! batching
!
         call mem%alloc(g_ljak, wf%n_o, wf%n_o, wf%n_o, batch_a%max_length)
         call mem%alloc(g_ljbk, wf%n_o, wf%n_o, wf%n_o, batch_a%max_length)
         call mem%alloc(g_ljck, wf%n_o, wf%n_o, wf%n_o, batch_a%max_length)
!
         call mem%alloc(g_bdak, wf%n_v, wf%n_o, batch_a%max_length, batch_a%max_length)
         call mem%alloc(g_cdak, wf%n_v, wf%n_o, batch_a%max_length, batch_a%max_length)
         call mem%alloc(g_cdbk, wf%n_v, wf%n_o, batch_a%max_length, batch_a%max_length)
         call mem%alloc(g_adbk, wf%n_v, wf%n_o, batch_a%max_length, batch_a%max_length)
         call mem%alloc(g_adck, wf%n_v, wf%n_o, batch_a%max_length, batch_a%max_length)
         call mem%alloc(g_bdck, wf%n_v, wf%n_o, batch_a%max_length, batch_a%max_length)
!
         call mem%alloc(g_ljak_c1, wf%n_o, wf%n_o, wf%n_o, batch_a%max_length)
         call mem%alloc(g_ljbk_c1, wf%n_o, wf%n_o, wf%n_o, batch_a%max_length)
         call mem%alloc(g_ljck_c1, wf%n_o, wf%n_o, wf%n_o, batch_a%max_length)
!
         call mem%alloc(g_bdak_c1, wf%n_v, wf%n_o, batch_a%max_length, batch_a%max_length)
         call mem%alloc(g_cdak_c1, wf%n_v, wf%n_o, batch_a%max_length, batch_a%max_length)
         call mem%alloc(g_cdbk_c1, wf%n_v, wf%n_o, batch_a%max_length, batch_a%max_length)
         call mem%alloc(g_adbk_c1, wf%n_v, wf%n_o, batch_a%max_length, batch_a%max_length)
         call mem%alloc(g_adck_c1, wf%n_v, wf%n_o, batch_a%max_length, batch_a%max_length)
         call mem%alloc(g_bdck_c1, wf%n_v, wf%n_o, batch_a%max_length, batch_a%max_length)
!
         call mem%alloc(g_jlka, wf%n_o, wf%n_o, wf%n_o, batch_a%max_length)
         call mem%alloc(g_jlkb, wf%n_o, wf%n_o, wf%n_o, batch_a%max_length)
         call mem%alloc(g_jlkc, wf%n_o, wf%n_o, wf%n_o, batch_a%max_length)
!
         call mem%alloc(g_dbka, wf%n_o, wf%n_v, batch_a%max_length, batch_a%max_length)
         call mem%alloc(g_dcka, wf%n_o, wf%n_v, batch_a%max_length, batch_a%max_length)
         call mem%alloc(g_dckb, wf%n_o, wf%n_v, batch_a%max_length, batch_a%max_length)
         call mem%alloc(g_dakb, wf%n_o, wf%n_v, batch_a%max_length, batch_a%max_length)
         call mem%alloc(g_dakc, wf%n_o, wf%n_v, batch_a%max_length, batch_a%max_length)
         call mem%alloc(g_dbkc, wf%n_o, wf%n_v, batch_a%max_length, batch_a%max_length)
!
         call mem%alloc(g_jakb, wf%n_o, wf%n_o, batch_a%max_length, batch_a%max_length)
         call mem%alloc(g_jakc, wf%n_o, wf%n_o, batch_a%max_length, batch_a%max_length)
         call mem%alloc(g_jbkc, wf%n_o, wf%n_o, batch_a%max_length, batch_a%max_length)
!
         if (wf%n_o .le. wf%n_v) then
            call mem%alloc(sorting, wf%n_o, wf%n_v, wf%n_v, batch_a%max_length)
         else
            call mem%alloc(sorting, wf%n_o, wf%n_o, wf%n_o, batch_a%max_length)
         end if
!
      endif
!
!     Loop over the batches in a,b,c
!     Read integrals and assign pointers
!     Without pointers we'll have to use three times as much
!     memory for the non-batching case
!
      do a_batch = 1, batch_a%num_batches
!
         call batch_a%determine_limits(a_batch)
!
         call wf%setup_oovo_abc(g_ljak, g_ljak_p, sorting, batch_a)
         call wf%setup_ooov_abc(g_jlka, g_jlka_p, sorting, batch_a)
         call wf%setup_oovo_abc(g_ljak_c1, g_ljak_c1_p, sorting, batch_a, c_ai=R_ai)
!
         do b_batch = 1, a_batch
!
            call batch_b%determine_limits(b_batch)
!
            call wf%setup_vvvo_abc(g_bdak, g_bdak_p, sorting, batch_b, batch_a)
            call wf%setup_vvov_abc(g_dbka, g_dbka_p, sorting, batch_b, batch_a)
            call wf%setup_ovov_abc(g_jakb, g_jakb_p, sorting, batch_a, batch_b)
            call wf%setup_vvvo_abc(g_bdak_c1, g_bdak_c1_p, sorting, batch_b, batch_a, c_ai=R_ai)
!
            if (b_batch .ne. a_batch) then
!
               call wf%setup_oovo_abc(g_ljbk, g_ljbk_p, sorting, batch_b)
               call wf%setup_ooov_abc(g_jlkb, g_jlkb_p, sorting, batch_b)
               call wf%setup_vvvo_abc(g_adbk, g_adbk_p, sorting, batch_a, batch_b)
               call wf%setup_vvov_abc(g_dakb, g_dakb_p, sorting, batch_a, batch_b)
               call wf%setup_oovo_abc(g_ljbk_c1, g_ljbk_c1_p, sorting, batch_b, c_ai=R_ai)
               call wf%setup_vvvo_abc(g_adbk_c1, g_adbk_c1_p, sorting, batch_a, batch_b, c_ai=R_ai)
!
            else
!
               call wf%point_oovo_abc(g_ljbk_p, g_ljak, batch_b%length)
               call wf%point_ooov_abc(g_jlkb_p, g_jlka, batch_b%length)
               call wf%point_vvvo_abc(g_adbk_p, g_bdak, batch_a%length, batch_b%length)
               call wf%point_vvov_abc(g_dakb_p, g_dbka, batch_a%length, batch_b%length)
               call wf%point_oovo_abc(g_ljbk_c1_p, g_ljak_c1, batch_b%length)
               call wf%point_vvvo_abc(g_adbk_c1_p, g_bdak_c1, batch_a%length, batch_b%length)
!
            endif
!
            do c_batch = 1, b_batch
!
               call batch_c%determine_limits(c_batch)
!
               if (c_batch .ne. b_batch) then ! c_batch != b_batch and c_batch != a_batch
!
                  call wf%setup_oovo_abc(g_ljck, g_ljck_p, sorting, batch_c)
!
                  call wf%setup_ooov_abc(g_jlkc, g_jlkc_p, sorting, batch_c)
!
                  call wf%setup_vvvo_abc(g_cdak, g_cdak_p, sorting, batch_c, batch_a)
                  call wf%setup_vvvo_abc(g_adck, g_adck_p, sorting, batch_a, batch_c)
                  call wf%setup_vvvo_abc(g_cdbk, g_cdbk_p, sorting, batch_c, batch_b)
                  call wf%setup_vvvo_abc(g_bdck, g_bdck_p, sorting, batch_b, batch_c)
!
                  call wf%setup_vvov_abc(g_dcka, g_dcka_p, sorting, batch_c, batch_a)
                  call wf%setup_vvov_abc(g_dakc, g_dakc_p, sorting, batch_a, batch_c)
                  call wf%setup_vvov_abc(g_dbkc, g_dbkc_p, sorting, batch_b, batch_c)
                  call wf%setup_vvov_abc(g_dckb, g_dckb_p, sorting, batch_c, batch_b)
!
                  call wf%setup_ovov_abc(g_jakc, g_jakc_p, sorting, batch_a, batch_c)
                  call wf%setup_ovov_abc(g_jbkc, g_jbkc_p, sorting, batch_b, batch_c)
!
                  call wf%setup_oovo_abc(g_ljck_c1, g_ljck_c1_p, sorting, batch_c, c_ai=R_ai)
!
                  call wf%setup_vvvo_abc(g_cdak_c1, g_cdak_c1_p, sorting, batch_c, batch_a, c_ai=R_ai)
                  call wf%setup_vvvo_abc(g_adck_c1, g_adck_c1_p, sorting, batch_a, batch_c, c_ai=R_ai)
                  call wf%setup_vvvo_abc(g_cdbk_c1, g_cdbk_c1_p, sorting, batch_c, batch_b, c_ai=R_ai)
                  call wf%setup_vvvo_abc(g_bdck_c1, g_bdck_c1_p, sorting, batch_b, batch_c, c_ai=R_ai)
!
               else if (c_batch .eq. a_batch) then ! c_batch = b_batch = a_batch
!
                  call wf%point_oovo_abc(g_ljck_p, g_ljak, batch_c%length)
                  call wf%point_ooov_abc(g_jlkc_p, g_jlka, batch_c%length)
!
                  call wf%point_vvvo_abc(g_cdak_p, g_bdak, batch_c%length, batch_a%length)
                  call wf%point_vvvo_abc(g_adck_p, g_bdak, batch_a%length, batch_c%length)
                  call wf%point_vvvo_abc(g_cdbk_p, g_bdak, batch_c%length, batch_b%length)
                  call wf%point_vvvo_abc(g_bdck_p, g_bdak, batch_b%length, batch_c%length)
!
                  call wf%point_vvov_abc(g_dcka_p, g_dbka, batch_c%length, batch_a%length)
                  call wf%point_vvov_abc(g_dakc_p, g_dbka, batch_a%length, batch_c%length)
                  call wf%point_vvov_abc(g_dckb_p, g_dbka, batch_c%length, batch_b%length)
                  call wf%point_vvov_abc(g_dbkc_p, g_dbka, batch_b%length, batch_c%length)
!
                  call wf%point_ovov_abc(g_jakc_p, g_jakb, batch_a%length, batch_c%length)
                  call wf%point_ovov_abc(g_jbkc_p, g_jakb, batch_b%length, batch_c%length)
!
                  call wf%point_oovo_abc(g_ljck_c1_p, g_ljak_c1, batch_c%length)
!
                  call wf%point_vvvo_abc(g_cdak_c1_p, g_bdak_c1, batch_c%length, batch_a%length)
                  call wf%point_vvvo_abc(g_adck_c1_p, g_bdak_c1, batch_a%length, batch_c%length)
                  call wf%point_vvvo_abc(g_cdbk_c1_p, g_bdak_c1, batch_c%length, batch_b%length)
                  call wf%point_vvvo_abc(g_bdck_c1_p, g_bdak_c1, batch_b%length, batch_c%length)
!
               else ! c_batch == b_batch != a_batch
!
                  call wf%point_oovo_abc(g_ljck_p, g_ljbk, batch_c%length)
                  call wf%point_ooov_abc(g_jlkc_p, g_jlkb, batch_c%length)
!
                  call wf%setup_vvvo_abc(g_cdbk, g_cdbk_p, sorting, batch_c, batch_b)
                  call wf%point_vvvo_abc(g_bdck_p, g_cdbk, batch_b%length, batch_c%length)
                  call wf%point_vvvo_abc(g_cdak_p, g_bdak, batch_c%length, batch_a%length)
                  call wf%point_vvvo_abc(g_adck_p, g_adbk, batch_a%length, batch_c%length)
!
                  call wf%setup_vvov_abc(g_dckb, g_dckb_p, sorting, batch_c, batch_b)
                  call wf%point_vvov_abc(g_dcka_p, g_dbka, batch_c%length, batch_a%length)
                  call wf%point_vvov_abc(g_dakc_p, g_dakb, batch_a%length, batch_c%length)
                  call wf%point_vvov_abc(g_dbkc_p, g_dckb, batch_b%length, batch_c%length)
!
                  call wf%setup_ovov_abc(g_jbkc, g_jbkc_p, sorting, batch_c, batch_b)
                  call wf%point_ovov_abc(g_jakc_p, g_jakb, batch_a%length, batch_c%length)
!
                  call wf%point_oovo_abc(g_ljck_c1_p, g_ljbk_c1, batch_c%length)
!
                  call wf%setup_vvvo_abc(g_cdbk_c1, g_cdbk_c1_p, sorting, batch_c, batch_b, c_ai=R_ai)
                  call wf%point_vvvo_abc(g_bdck_c1_p, g_cdbk_c1, batch_b%length, batch_c%length)
                  call wf%point_vvvo_abc(g_cdak_c1_p, g_bdak_c1, batch_c%length, batch_a%length)
                  call wf%point_vvvo_abc(g_adck_c1_p, g_adbk_c1, batch_a%length, batch_c%length)
!
               endif
!
!$omp parallel do private(a, a_rel, b, b_rel, c, c_rel, thread_n) & 
!$omp schedule(guided, 1)
               do a = batch_a%first, batch_a%last
!
                  a_rel = a - batch_a%first + 1
!
!$                thread_n = omp_get_thread_num() + 1
!
                  do b = batch_b%first, min(batch_b%last, a)
!
                     b_rel = b - batch_b%first + 1
!
                     do c = batch_c%first, min(batch_c%last, b)
!
                        if (a .eq. b .and. a .eq. c) then
                           cycle
                        end if
!
                        c_rel = c - batch_c%first + 1
!
!                       Construct R^abc_ijk for given a,b,c
!                       Using c1-transformed integrals the terms have the same form 
!                       as the omega terms (where t_ijk = R_ijk)
!
!                       Therefore the contributions to the c3-amplitudes can be computed 
!                       using the same routine once for t1-transformed and once for 
!                       c1-transformed integrals
!
                        call wf%construct_W_abc(a, b, c,                   &
                                                R_ijk(:,:,:,thread_n),     &
                                                u_ijk(:,:,:,thread_n),     &
                                                R_ijab,                    &
                                                g_ljak_p(:,:,:,a_rel),     &
                                                g_ljbk_p(:,:,:,b_rel),     &
                                                g_ljck_p(:,:,:,c_rel),     &
                                                g_bdak_p(:,:,b_rel,a_rel), &
                                                g_cdak_p(:,:,c_rel,a_rel), &
                                                g_cdbk_p(:,:,c_rel,b_rel), &
                                                g_adbk_p(:,:,a_rel,b_rel), &
                                                g_adck_p(:,:,a_rel,c_rel), &
                                                g_bdck_p(:,:,b_rel,c_rel))
!
                        call wf%construct_W_abc(a, b, c,                      &
                                                R_ijk(:,:,:,thread_n),        &
                                                u_ijk(:,:,:,thread_n),        &
                                                t_ijab,                       &
                                                g_ljak_c1_p(:,:,:,a_rel),     &
                                                g_ljbk_c1_p(:,:,:,b_rel),     &
                                                g_ljck_c1_p(:,:,:,c_rel),     &
                                                g_bdak_c1_p(:,:,b_rel,a_rel), &
                                                g_cdak_c1_p(:,:,c_rel,a_rel), &
                                                g_cdbk_c1_p(:,:,c_rel,b_rel), &
                                                g_adbk_c1_p(:,:,a_rel,b_rel), &
                                                g_adck_c1_p(:,:,a_rel,c_rel), &
                                                g_bdck_c1_p(:,:,b_rel,c_rel), &
                                                overwrite = .false.) ! Do not overwrite R_ijk
!
                        call wf%divide_by_orbital_differences_abc(a, b, c, &
                                                                  R_ijk(:,:,:,thread_n), &
                                                                  omega, wf%cvs)
!
                        call wf%construct_W_abc(a, b, c,                   &
                                                tbar_ijk(:,:,:,thread_n),  &
                                                v_ijk(:,:,:,thread_n),     &
                                                tbar_ijab,                 &
                                                g_jlka_p(:,:,:,a_rel),     &
                                                g_jlkb_p(:,:,:,b_rel),     &
                                                g_jlkc_p(:,:,:,c_rel),     &
                                                g_dbka_p(:,:,b_rel,a_rel), &
                                                g_dcka_p(:,:,c_rel,a_rel), &
                                                g_dckb_p(:,:,c_rel,b_rel), &
                                                g_dakb_p(:,:,a_rel,b_rel), &
                                                g_dakc_p(:,:,a_rel,c_rel), &
                                                g_dbkc_p(:,:,b_rel,c_rel))
!
                        call wf%outer_product_terms_l3_abc(a, b, c, tbar_ia,          &
                                                           tbar_ijab,                 &
                                                           tbar_ijk(:,:,:,thread_n),  & 
                                                           wf%fock_ia,                &
                                                           g_jakb_p(:,:,a_rel,b_rel), &
                                                           g_jakc_p(:,:,a_rel,c_rel), &
                                                           g_jbkc_p(:,:,b_rel,c_rel))
!
                        call construct_contravariant_t3(tbar_ijk(:,:,:,thread_n), &
                                                               v_ijk(:,:,:,thread_n), wf%n_o)
!
                        call wf%divide_by_orbital_differences_abc(a, b, c, tbar_ijk(:,:,:,thread_n))
!
                        call wf%density_cc3_mu_ref_oo(a, b, c, density_oo_thread(:,:,thread_n), &
                                                      R_ijk(:,:,:,thread_n),                    &
                                                      u_ijk(:,:,:,thread_n),                    &
                                                      tbar_ijk(:,:,:,thread_n),                 &
                                                      v_ijk(:,:,:,thread_n))
!
                     enddo ! loop over c
                  enddo ! loop over b
               enddo ! loop over a
!$omp end parallel do
            enddo ! batch_c
         enddo ! batch_b
      enddo ! batch_a
!
!     Collect contributions to the oo-block from all threads
!
      do a = 1, n_threads
         call daxpy(wf%n_o**2, one, density_oo_thread(:,:,a), 1, density_oo, 1)
      enddo
!
      call mem%dealloc(density_oo_thread, wf%n_o, wf%n_o, n_threads)
!
      call mem%dealloc(R_ijk, wf%n_o, wf%n_o, wf%n_o, n_threads)
      call mem%dealloc(u_ijk, wf%n_o, wf%n_o, wf%n_o, n_threads)
      call mem%dealloc(v_ijk, wf%n_o, wf%n_o, wf%n_o, n_threads)
      call mem%dealloc(tbar_ijk, wf%n_o, wf%n_o, wf%n_o, n_threads)
!
      call mem%dealloc(t_ijab, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
      if (batch_a%num_batches .eq. 1) then ! no batching
!
         call mem%dealloc(g_ljak, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
         call mem%dealloc(g_bdak, wf%n_v, wf%n_o, wf%n_v, wf%n_v)
!
         call mem%dealloc(g_jlka, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
         call mem%dealloc(g_dbka, wf%n_o, wf%n_v, wf%n_v, wf%n_v)
         call mem%dealloc(g_jakb, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
         call mem%dealloc(g_ljak_c1, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
         call mem%dealloc(g_bdak_c1, wf%n_v, wf%n_o, wf%n_v, wf%n_v)
!
         if (wf%n_o .lt. wf%n_v) then
            call mem%dealloc(sorting, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         else
            call mem%dealloc(sorting, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         end if      
!
      else ! batching
!
         call mem%dealloc(g_ljak, wf%n_o, wf%n_o, wf%n_o, batch_a%max_length)
         call mem%dealloc(g_ljbk, wf%n_o, wf%n_o, wf%n_o, batch_a%max_length)
         call mem%dealloc(g_ljck, wf%n_o, wf%n_o, wf%n_o, batch_a%max_length)
!
         call mem%dealloc(g_bdak, wf%n_v, wf%n_o, batch_a%max_length, batch_a%max_length)
         call mem%dealloc(g_cdak, wf%n_v, wf%n_o, batch_a%max_length, batch_a%max_length)
         call mem%dealloc(g_cdbk, wf%n_v, wf%n_o, batch_a%max_length, batch_a%max_length)
         call mem%dealloc(g_adbk, wf%n_v, wf%n_o, batch_a%max_length, batch_a%max_length)
         call mem%dealloc(g_adck, wf%n_v, wf%n_o, batch_a%max_length, batch_a%max_length)
         call mem%dealloc(g_bdck, wf%n_v, wf%n_o, batch_a%max_length, batch_a%max_length)
!
         call mem%dealloc(g_ljak_c1, wf%n_o, wf%n_o, wf%n_o, batch_a%max_length)
         call mem%dealloc(g_ljbk_c1, wf%n_o, wf%n_o, wf%n_o, batch_a%max_length)
         call mem%dealloc(g_ljck_c1, wf%n_o, wf%n_o, wf%n_o, batch_a%max_length)
!
         call mem%dealloc(g_bdak_c1, wf%n_v, wf%n_o, batch_a%max_length, batch_a%max_length)
         call mem%dealloc(g_cdak_c1, wf%n_v, wf%n_o, batch_a%max_length, batch_a%max_length)
         call mem%dealloc(g_cdbk_c1, wf%n_v, wf%n_o, batch_a%max_length, batch_a%max_length)
         call mem%dealloc(g_adbk_c1, wf%n_v, wf%n_o, batch_a%max_length, batch_a%max_length)
         call mem%dealloc(g_adck_c1, wf%n_v, wf%n_o, batch_a%max_length, batch_a%max_length)
         call mem%dealloc(g_bdck_c1, wf%n_v, wf%n_o, batch_a%max_length, batch_a%max_length)
!
         call mem%dealloc(g_jlka, wf%n_o, wf%n_o, wf%n_o, batch_a%max_length)
         call mem%dealloc(g_jlkb, wf%n_o, wf%n_o, wf%n_o, batch_a%max_length)
         call mem%dealloc(g_jlkc, wf%n_o, wf%n_o, wf%n_o, batch_a%max_length)
!
         call mem%dealloc(g_dbka, wf%n_o, wf%n_v, batch_a%max_length, batch_a%max_length)
         call mem%dealloc(g_dcka, wf%n_o, wf%n_v, batch_a%max_length, batch_a%max_length)
         call mem%dealloc(g_dckb, wf%n_o, wf%n_v, batch_a%max_length, batch_a%max_length)
         call mem%dealloc(g_dakb, wf%n_o, wf%n_v, batch_a%max_length, batch_a%max_length)
         call mem%dealloc(g_dakc, wf%n_o, wf%n_v, batch_a%max_length, batch_a%max_length)
         call mem%dealloc(g_dbkc, wf%n_o, wf%n_v, batch_a%max_length, batch_a%max_length)
!
         call mem%dealloc(g_jakb, wf%n_o, wf%n_o, batch_a%max_length, batch_a%max_length)
         call mem%dealloc(g_jakc, wf%n_o, wf%n_o, batch_a%max_length, batch_a%max_length)
         call mem%dealloc(g_jbkc, wf%n_o, wf%n_o, batch_a%max_length, batch_a%max_length)
!
         if (wf%n_o .le. wf%n_v) then
            call mem%dealloc(sorting, wf%n_o, wf%n_v, wf%n_v, batch_a%max_length)
         else
            call mem%dealloc(sorting, wf%n_o, wf%n_o, wf%n_o, batch_a%max_length)
         end if
!
      endif
!
   end subroutine density_cc3_mu_nu_abc_cc3
!
!
   module subroutine density_cc3_mu_nu_ov_Z_term_cc3(wf, density_ov, Z_bcjk, t_abij)
!!
!!    One electron density (EOM) excited-determinant/excited-determinant 
!!    ov-term from Z_intermediate
!!    Written by Alexander C. Paul, August 2019
!!
!!    Computes terms of the form:
!!
!!          D_pq += sum_mu,nu X_mu < mu| e^(-T) E_pq e^T |nu > Y_nu
!!
!!    explicit term in this routine:
!!
!!       D^R_ld += sum_{abcijk} tbar^abc_ijk R^a_i (t^bcd_jkl - t^bcd_jlk)
!!              += sum_{bcjk} Z_bcjk (t^bcd_jkl - t^bcd_jlk)
!!              += sum_{bcjk} ~Z_bcjk t^bcd_jkl
!!
!!    Construct t3-amplitudes in batches of i,j,k and calculate the
!!    contribution to the ov-block involving the Z-intermediate
!!
!!       t_mu3 = -< mu3| [U,T2] |HF > (eps_mu3)^-1
!!
!!    NB: the covariant Z intermediate (~Z = 1/3 (2 Z_bcjk + Z_bckj))
!!        is used in this routine.
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_v), intent(inout) :: density_ov
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: Z_bcjk
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) ::t_abij
!
      real(dp), dimension(:,:,:), allocatable :: t_abc
!
!     Help array used for sorting integrals and amplitudes
      real(dp), dimension(:,:,:,:), allocatable :: sorting
!
!     Integrals for the construction of the T3-amplitudes
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdci
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdcj
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdck
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdci_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdcj_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdck_p => null()
!
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ljci
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_lkci
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_lkcj
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_licj
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_lick
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ljck
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_ljci_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_lkci_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_lkcj_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_licj_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_lick_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_ljck_p => null()
!
      real(dp), dimension(:,:), allocatable :: density_ai
!
      integer :: i, j, k, i_rel, j_rel, k_rel
      type(batching_index) :: batch_i, batch_j, batch_k
      integer :: i_batch, j_batch, k_batch
      integer :: req_0, req_1, req_2, req_3, req_i, req_1_eri
      integer :: req_single_batch
!
      call mem%alloc(density_ai, wf%n_v, wf%n_o)
      call zero_array(density_ai, wf%n_t1)
!
      batch_i = batching_index(wf%n_o)
      batch_j = batching_index(wf%n_o)
      batch_k = batching_index(wf%n_o)
!
!     Memory for sorting array and getting the integrals
      call wf%estimate_mem_integral_setup(req_0, req_1_eri)
      req_0 = req_0 + wf%n_v**3
      req_1_eri = req_1_eri + max(wf%n_v**3, wf%n_o**2*wf%n_v)
!
!     Need less memory if we don't need to batch, so we overwrite the maximum 
!     required memory in batch_setup
!
      req_single_batch = req_0 + req_1_eri*wf%n_o + wf%n_v**3*wf%n_o &
                       + wf%n_v*wf%n_o**3
!
      req_1 = wf%n_v**3
      req_i = req_1 + req_1_eri ! Mem for integral setup only needed for 1 index.
      req_2 = 2*wf%n_o*wf%n_v
      req_3 = 0
!
      call mem%batch_setup(batch_i, batch_j, batch_k,  &
                           req_0, req_i, req_1, req_1, &
                           req_2, req_2, req_2, req_3, &
                           req_single_batch=req_single_batch)
!
      call mem%alloc(t_abc, wf%n_v, wf%n_v, wf%n_v)
!
      if (batch_i%num_batches .eq. 1) then ! no batching
!
         call mem%alloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%alloc(g_ljci, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
         if (wf%n_o .le. wf%n_v) then
            call mem%alloc(sorting, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         else
            call mem%alloc(sorting, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         end if
!
      else ! batching
!
         call mem%alloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%alloc(g_bdcj, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%alloc(g_bdck, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
!
         call mem%alloc(g_ljci, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_lkci, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_lkcj, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_licj, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_lick, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_ljck, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
!
         if (wf%n_o .le. wf%n_v) then
            call mem%alloc(sorting, wf%n_v, wf%n_v, wf%n_v,batch_i%max_length)
         else
            call mem%alloc(sorting, wf%n_v, wf%n_o, wf%n_o,batch_i%max_length)
         end if
!
      endif
!
!     Loop over the batches in i,j,k
!     Read integrals and assign pointers
!     Without pointers we'll have to use three times as much
!     memory for the non-batching case
!
      do i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(i_batch)
!
         call wf%setup_vvvo(g_bdci, g_bdci_p, sorting, batch_i)
!
         do j_batch = 1, i_batch
!
            call batch_j%determine_limits(j_batch)
!
            call wf%setup_oovo(g_ljci, g_ljci_p, sorting, batch_j, batch_i)
!
            if (j_batch .ne. i_batch) then
!
               call wf%setup_vvvo(g_bdcj, g_bdcj_p, sorting, batch_j)
!
               call wf%setup_oovo(g_licj, g_licj_p, sorting, batch_i, batch_j)
!
            else
!
               call wf%point_vvvo(g_bdcj_p, g_bdci, batch_j%length)
!
               call wf%point_vooo(g_licj_p, g_ljci, batch_i%length, batch_j%length)
!
            endif
!
            do k_batch = 1, j_batch
!
               call batch_k%determine_limits(k_batch)
!
               if (k_batch .ne. j_batch) then ! k_batch != j_batch, k_batch != i_batch
!
                  call wf%setup_vvvo(g_bdck, g_bdck_p, sorting, batch_k)
!
                  call wf%setup_oovo(g_lick, g_lick_p, sorting, batch_i, batch_k)
                  call wf%setup_oovo(g_ljck, g_ljck_p, sorting, batch_j, batch_k)
                  call wf%setup_oovo(g_lkci, g_lkci_p, sorting, batch_k, batch_i)
                  call wf%setup_oovo(g_lkcj, g_lkcj_p, sorting, batch_k, batch_j)
!
!
               else if (k_batch .eq. i_batch) then !k_batch == j_batch == i_batch
!
                  call wf%point_vvvo(g_bdck_p, g_bdci, batch_k%length)
!
                  call wf%point_vooo(g_lick_p, g_ljci, batch_i%length, batch_k%length)
                  call wf%point_vooo(g_ljck_p, g_ljci, batch_j%length, batch_k%length)
                  call wf%point_vooo(g_lkci_p, g_ljci, batch_k%length, batch_i%length)
                  call wf%point_vooo(g_lkcj_p, g_ljci, batch_k%length, batch_j%length)
!
               else ! k_batch == j_batch != i_batch
!
                  call wf%point_vvvo(g_bdck_p, g_bdcj, batch_k%length)
!
                  call wf%setup_oovo(g_lkcj, g_lkcj_p, sorting, batch_k, batch_j)
                  call wf%point_vooo(g_lick_p, g_licj, batch_i%length, batch_k%length)
                  call wf%point_vooo(g_ljck_p, g_lkcj, batch_j%length, batch_k%length)
                  call wf%point_vooo(g_lkci_p, g_ljci, batch_k%length, batch_i%length)
!
               endif
!
               do i = batch_i%first, batch_i%last
!
                  i_rel = i - batch_i%first + 1
!
                  do j = batch_j%first, min(batch_j%last, i)
!
                     j_rel = j - batch_j%first + 1
!
                     do k = batch_k%first, min(batch_k%last, j)
!
                        if (i .eq. j .and. i .eq. k) then
                           cycle
                        end if
!
                        k_rel = k - batch_k%first + 1
!
!                       construct t3 for fixed i,j,k
                        call wf%construct_W(i, j, k, sorting, t_abc, t_abij, &
                                            g_bdci_p(:,:,:,i_rel),          &
                                            g_bdcj_p(:,:,:,j_rel),          &
                                            g_bdck_p(:,:,:,k_rel),          &
                                            g_ljci_p(:,:,j_rel,i_rel),      &
                                            g_lkci_p(:,:,k_rel,i_rel),      &
                                            g_lkcj_p(:,:,k_rel,j_rel),      &
                                            g_licj_p(:,:,i_rel,j_rel),      &
                                            g_lick_p(:,:,i_rel,k_rel),      &
                                            g_ljck_p(:,:,j_rel,k_rel))
!
                        call construct_contravariant_t3(t_abc, sorting, wf%n_v)
                        call wf%divide_by_orbital_differences(i, j, k, t_abc)
!
                        call wf%construct_x_ai_intermediate(i, j, k, t_abc, sorting, &
                                                            Z_bcjk, density_ai)
!
                     enddo ! loop over k
                  enddo ! loop over j
               enddo ! loop over i
            enddo ! batch_k
         enddo ! batch_j
      enddo ! batch_i
!
      call add_21_to_12(one, density_ai, density_ov, wf%n_o, wf%n_v)
!
      call mem%dealloc(density_ai, wf%n_v, wf%n_o)
!
      call mem%dealloc(t_abc, wf%n_v, wf%n_v, wf%n_v)
!
      if (batch_i%num_batches .eq. 1) then ! no batching
!
         call mem%dealloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%dealloc(g_ljci, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
         if (wf%n_o .le. wf%n_v) then
            call mem%dealloc(sorting, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         else
            call mem%dealloc(sorting, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         end if
!
      else ! batching
!
         call mem%dealloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%dealloc(g_bdcj, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%dealloc(g_bdck, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
!
         call mem%dealloc(g_ljci, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_lkci, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_lkcj, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_licj, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_lick, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_ljck, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
!
         if (wf%n_o .le. wf%n_v) then
            call mem%dealloc(sorting, wf%n_v, wf%n_v, wf%n_v,batch_i%max_length)
         else
            call mem%dealloc(sorting, wf%n_v, wf%n_o, wf%n_o,batch_i%max_length)
         end if
!
      endif
!
   end subroutine density_cc3_mu_nu_ov_Z_term_cc3
!
!
end submodule
