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
submodule (ccs_class) fop_ccs
!
!!
!!    First order properties submodule (CCS)
!!    Written by Josefine H. Andersen, Mar 2019
!!
!!    Adapted by Sarai D. Folkestad
!!
!!    Routines for construction of the right-hand-side, η^X, and left-hand-side, ξ^X
!!    vectors and the left-hand-side (ρ^L) and right-hand-side (ρ^R) transition densities
!!    for transition moments.
!!
!!    Equation-of-motion (EOM):
!!
!!    (Following Koch, H., Kobayashi, R., Sanches de Merás, A., and Jørgensen, P.,
!!    J. Chem. Phys. 100, 4393 (1994))
!!
!!       η_μ^X,EOM =  < Λ | [X, τ_μ] | CC > + (< Λ | τ_μ X | CC >  - tbar_μ < Λ | X | CC > )
!!                 = η^{X,0} + η^{X,corr}
!!
!!    Where the last two terms are called the EOM-corrections and the first term also 
!!    appears in LR-CC.
!!
!!    The left-hand-side vector is the same in EOM-CC and LR-CC:   
!!
!!       ξ^X_μ = < μ | exp(-T) X exp(T)| R >
!!
!!    The transition density matrices are construct as follows:
!!
!!       ρ^L_pq = < k | E_pq | CC >
!!       ρ^R_pq = < Λ | E_pq | k >
!!
!!    where |k> and <k| are the eigenvectors of the Jacobian with amplitudes R_μ, L_μ
!!
!!       | k > = sum_μ (τ_μ | CC > R_{k,μ} - tbar_μ | CC > R_{k,μ})
!!       < k | = sum_μ L_{k,μ} < μ | e^-T
!!
!
   implicit none
!
!
contains
!
!
   module subroutine construct_right_transition_density_ccs(wf, R_k)
!!
!!    Construct right one-electron transition density for the state k
!!    Written by Alexander Paul, June 2019
!!
!!          ρ^R_pq = < Λ | E_pq | k >
!!
!!    where |k> is the right eigenvector of the Jacobian
!!    with amplitudes R_μ
!!
!!          | k > = sum_μ (τ_μ | CC > R_{k,μ} - tbar_μ | CC > R_{k,μ}) 
!!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: R_k
!
      call zero_array(wf%right_transition_density, (wf%n_mo)**2)
!
      call wf%right_transition_density_ccs_oo(wf%t1bar, R_k)
      call wf%right_transition_density_ccs_ov(R_k)
      call wf%right_transition_density_ccs_vv(wf%t1bar, R_k)
      call wf%right_transition_density_ccs_gs_contr(wf%t1bar, R_k)
!
   end subroutine construct_right_transition_density_ccs
!
!
   module subroutine right_transition_density_ccs_oo_ccs(wf, tbar_ai, R_ai)
!!
!!    Right transition density oo contribution
!!    Written by Alexander Paul, June 2019
!!
!!    ρ^R_ij = -sum_a R_ai tbar_aj
!!      
      implicit none
!
      class(ccs) :: wf
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
                  wf%right_transition_density,  &
                  wf%n_mo)
!
   end subroutine right_transition_density_ccs_oo_ccs
!
!
   module subroutine right_transition_density_ccs_ov_ccs(wf, R_ai)
!!
!!    Right transition density ov
!!    Written by Alexander Paul, June 2019
!!
!!    ρ^R_pq = 2*R_ai 
!!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: R_ai
!
      integer :: i, a
!
!$omp parallel do private(a, i)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            wf%right_transition_density(i, wf%n_o + a) = two*R_ai(a, i) &
                           + wf%right_transition_density(i, wf%n_o + a)
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine right_transition_density_ccs_ov_ccs
!
!
   module subroutine right_transition_density_ccs_vv_ccs(wf, tbar_ai, R_ai)
!!
!!    Right transition density vv contribution
!!    Written by Alexander Paul, June 2019
!!
!!    ρ^R_ab = sum_i R_bi tbar_ai
!!      
      implicit none
!
      class(ccs) :: wf
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
                  wf%right_transition_density(wf%n_o+1, wf%n_o+1),  &
                  wf%n_mo)
!
   end subroutine right_transition_density_ccs_vv_ccs
!
!
   module subroutine right_transition_density_ccs_gs_contr_ccs(wf, tbar_ai, R_ai)
!!
!!    Right transition density, contribution from the ground state density
!!    ρ^R_pq -= sum_μν R_{k,μ}tbar_μ tbar_ν < ν |e^-T E_pq e^T| HF >
!!           -= sum_μ R_{k,μ}tbar_μ (D_GS - D_HF)
!!    CCS:
!!    ρ^R_pq = ρ^R_ai = sum_bj R^b_j tbar^b_j tbar^a_i
!!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: tbar_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: R_ai
!
      real(dp) :: ddot, scaling_factor
!
      integer :: i, a
!
      scaling_factor = ddot(wf%n_v*wf%n_o, tbar_ai, 1, R_ai, 1)
!
!$omp parallel do private(a, i)
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            wf%right_transition_density(wf%n_o + a, i) = scaling_factor*tbar_ai(a, i)  &
                     + wf%right_transition_density(wf%n_o + a, i)
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine right_transition_density_ccs_gs_contr_ccs
!
!
   module subroutine construct_left_transition_density_ccs(wf, L_k)
!!
!!    Construct left one-electron transition density for the state k
!!    Written by Alexander Paul, June 2019
!!
!!          ρ^L_pq = < k | E_pq | CC >
!!
!!    where <k| is the left eigenvector of the Jacobian
!!    with amplitudes L_μ
!!
!!          < k | = sum_μ L_{k,μ} < μ | e^-T
!!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: L_k
!
      call zero_array(wf%left_transition_density, (wf%n_mo)**2)
!
      call wf%gs_one_el_density_ccs_vo(wf%left_transition_density, L_k)
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
!!    Construct η^X
!!    Written by Josefine H. Andersen, 2019
!!
!!    Adapted by Sarai D. Folekstad, Apr 2019
!!
!!    Constructs left-hand-side vector etaX:
!!
!!       η^X_μ = < Λ | [X, τ_μ] | CC >
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
!!    Adds the A1 term of η_ai^X:
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
!!    Adds the B1 term of η_ai^X:
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
!!    Constructs ξ^X_μ :
!!
!!       ξ^X_μ = < μ | exp(-T) X exp(T)| R >
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
!!       ξ^X_ai =+ X_ai
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
!!       A:  η^X,corr_μ += tbar_μ (ξ * tbar) 
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
      real(dp), dimension(:), allocatable :: multipliers
!
      call mem%alloc(multipliers, wf%n_es_amplitudes)
!
      call dcopy(wf%n_t1, wf%t1bar, 1, multipliers, 1)
!
      X_cc = ddot(wf%n_es_amplitudes, multipliers, 1, csiX, 1)
!
      call daxpy(wf%n_es_amplitudes, -X_cc, multipliers, 1, etaX, 1)
!
      call mem%dealloc(multipliers, wf%n_es_amplitudes)
!
   end subroutine etaX_eom_a_ccs
!
!
   module subroutine calculate_transition_strength_ccs(wf, S, etaX, csiX, state, T_l, T_r)
!!
!!    Calculate transition strength
!!    Written by Josefine H. Andersen, February 2019
!!
!!    Given etaX and csiX, this routine calculates the left and right transition 
!!    moments T_l and T_r for the state number "state" and the transition strength 
!!    S = T_l * T_r.
!! 
!!    The left and right states L and R are read from file and made binormal by the routine.
!!
!!    Consistency/sanity check: Eirik F. Kjønstad, Aug 2019
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
      real(dp), intent(out) :: T_l, T_r
      integer, intent(in)   :: state
!
      real(dp), dimension(:), allocatable :: L, R
!
      real(dp) :: ddot
!
!     Read states and make them binormal by scaling the left vector 
!
      call mem%alloc(L, wf%n_es_amplitudes)
      call mem%alloc(R, wf%n_es_amplitudes)
!
      call wf%binormalize_L_wrt_R(L, R, state)
!
!     Left and right transition moments
!
      T_r = ddot(wf%n_es_amplitudes, etaX, 1, R, 1)
      T_l = ddot(wf%n_es_amplitudes, L, 1, csiX, 1)
!
!     Transition strength
!
      S  = T_l * T_r
!
      call mem%dealloc(L, wf%n_es_amplitudes)
      call mem%dealloc(R, wf%n_es_amplitudes)
!
   end subroutine calculate_transition_strength_ccs
!
!
end submodule fop_ccs
