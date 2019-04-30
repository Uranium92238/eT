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
!!    vectors for transition moments.
!!
!!    Equation-of-motion (EOM):
!!
!!   (Followinng Koch, K., Kobayashi, R., Sanches de Merás, A., and Jørgensen, P.,
!!    J. Chem. Phys. 100, 4393 (1994))
!!
!!       η_μ^X,EOM =  < Λ | [X, τ_μ] | CC > + (< Λ | τ_μ X | R >  - tbar_μ < Λ | X | R > )
!!
!!
!!    Where the last two terms are called the EOM-corrections and the first term also 
!!    appears in LR-CC.
!!
!!    The left-hand-side vector is the same in EOM-CC and LR-CC:   
!!
!!       ξ^X_μ = < μ | exp(-T) X exp(T)| R >
!!
!
   implicit none
!
!
contains
!
   module subroutine construct_etaX_ccs(wf, X, etaX)
!!
!!    Construct η^X
!!    Written by Josefine H. Andersen, 2019
!!
!!    Adapted by Sarai D. Folekstad, Apr 2019
!!
!!    Handling the construction of left-hand-side vector etaX 
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
      etaX = zero
!
      call wf%etaX_ccs_a1(X, etaX)
      call wf%etaX_ccs_b1(X, etaX)
!
   end subroutine construct_etaX_ccs
!
!
   module subroutine etaX_ccs_a1_ccs(wf, X, etaX_ai)
!!
!!    Construct etaX A1 (CCS)
!!    Written by Josefine H. Andersen, 2019
!!
!!    Adapted by Sarai D. Folekstad, Apr 2019
!!
!!    Constructs the A1 term of η_ai^X:
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
!!    Constructs the B1 term of η_ai^X:
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
!!    Handling the construction of ξ^X_μ :
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
      csiX = zero
!
      call wf%csiX_ccs_a1(X, csiX)
!
   end subroutine construct_csiX_ccs
!
!
   module subroutine csiX_ccs_a1_ccs(wf, X, csiX_ai)
!!
!!    Construct csiX 
!!    Written by Josefine H. Andersen, Feb 2019
!!
!!       ξ^X_ai = X_ai
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
   module subroutine get_eom_contribution_ccs(wf, etaX, csiX, X)
!!
!!    Add EOM contribution to etaX vector
!!    Written by Josefine H. Andersen, Feb 2019
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: X
!
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: etaX
      real(dp), dimension(wf%n_es_amplitudes), intent(in)    :: csiX
!
      call wf%get_eom_xcc_contribution(etaX, csiX)
!
   end subroutine get_eom_contribution_ccs
!
!
   module subroutine get_eom_xcc_contribution_ccs(wf, etaX, csiX)
!!
!!    Get eom contribution
!!    Written by Josefine H. Andersen, Feb 2019
!!
!!    Add EOM contribution to etaX vector
!!
!!       EOM correction:  η^X,corr_μ += tbar_μ (ξ * tbar) 
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
      call wf%get_multipliers(multipliers)
!
      X_cc = ddot(wf%n_es_amplitudes, multipliers, 1, csiX, 1)
!
      call daxpy(wf%n_es_amplitudes, -X_cc, multipliers, 1, etaX, 1)
!
      call mem%dealloc(multipliers, wf%n_es_amplitudes)
!
   end subroutine get_eom_xcc_contribution_ccs
!
   module subroutine scale_left_excitation_vector_ccs(wf, L, R)
!!
!!    Make left and right excitation vectors biorthogonal by scaling left vector
!!    Written by Josefine H. Andersen, Feb 2019
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: L
      real(dp), dimension(wf%n_es_amplitudes), intent(in)    :: R
!
      real(dp) :: norm, ddot
!
      norm = ddot(wf%n_es_amplitudes, L, 1, R, 1)
!
      call dscal(wf%n_es_amplitudes, one/norm, L, 1)
!
   end subroutine scale_left_excitation_vector_ccs
!
!
   module subroutine calculate_transition_strength_ccs(wf, S, etaX, csiX, state, T_l, T_r)
!!
!!    Calculate transition strength for spectra
!!    Written by Josefine H. Andersen, February 2019
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
      real(dp), dimension(:), allocatable :: L_n, R_n
!
      real(dp) :: ddot
!
      call mem%alloc(L_n, wf%n_es_amplitudes)
      call mem%alloc(R_n, wf%n_es_amplitudes)
!
      call wf%read_excited_state(L_n, state, 'left')
      call wf%read_excited_state(R_n, state, 'right')
!
      call wf%scale_left_excitation_vector(L_n, R_n)
!
!     Left and right transition moments
!
      T_r = ddot(wf%n_es_amplitudes, etaX, 1, R_n, 1)
      T_l = ddot(wf%n_es_amplitudes, L_n, 1, csiX, 1)
!
!     Transition strength
!
      S  = T_l * T_r
!
      call mem%dealloc(L_n, wf%n_es_amplitudes)
      call mem%dealloc(R_n, wf%n_es_amplitudes)
!
   end subroutine calculate_transition_strength_ccs
!
end submodule fop_ccs
