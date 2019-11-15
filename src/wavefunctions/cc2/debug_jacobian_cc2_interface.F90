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
   module subroutine omega_for_jacobian_debug_cc2(wf, omega, t)
!!
!!    Omega for Jacobian debug
!!    Written by Sarai D. Folkestad, Sep. 2019
!!
!!    Constructs Ω for CC2, both singles and
!!    doubles part. 
!!
!!    The construct omega routine in CC2
!!    calculates u from the integrals, hence 
!!    it is not called directly (as is done 
!!    for CCS and CCSD), but u is constructed and 
!!    the individual terms are called (ccs_a1, 
!!    doubles_a1-doubles_c1). Additionally, the 
!!    doubles part of Ω is constructed. This routine 
!!    is only used in this debug procedure.
!!
!!
      implicit none
!
      class(cc2), intent(inout) :: wf
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: omega
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: t
!
   end subroutine omega_for_jacobian_debug_cc2
!
!
   module subroutine amplitudes_for_jacobian_debug_cc2(wf, t)
!!
!!    Amplitudes for Jacobian debug
!!    Written by Sarai D. Folkestad, Sep 2019
!!
!!    Calculates the amplitudes 
!!    with dimension n_es_amplitudes
!!
!!    Gets the t1-amplitudes 
!!    Constructs the t2-amplitudes
!!
      implicit none
!
      class(cc2), intent(inout) :: wf
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: t
!
   end subroutine amplitudes_for_jacobian_debug_cc2
!
!
   module subroutine normalization_for_jacobian_debug_cc2(wf, A_numerical_mu_nu, nu)
!!
!!    Normalization for Jacobian debug
!!    Written by Sarai D. Folkestad, Sep 2019
!!
!!    Differentiation wrt. doubles amplitudes (nu > wf%n_t1)
!!    will yield factor one half on the diagonal nu = aiai
!!
!!    This routine scales these diagonal elements 
!!    by factor two
!!
!!    Differentiation of Ω_μ2 wrt. doubles amplitudes
!!    gives factor two on the off-diagonal aibj, ai .ne. bj
!!
!!    A_numerical_μ2_aibj = dΩ_μ2/dt_aibj = 2 ε_μ2,aibj δ_μ2,aibj
!!
!!    For CC2, the doubles-doubles block of A is defined as
!!
!!    A_μ2,ν2 = delta_μ2,ν2 ε_ν2
!!
!!    This routine scales these off-diagonal elements 
!!    by factor half
!!
      implicit none
!
      class(cc2), intent(in) :: wf
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: A_numerical_mu_nu
      integer, intent(in) :: nu
!
   end subroutine normalization_for_jacobian_debug_cc2
