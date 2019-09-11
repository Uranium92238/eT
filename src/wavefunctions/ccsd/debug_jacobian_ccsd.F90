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
submodule (ccsd_class) debug_jacobian_ccsd
!
!!
!!    Debug Jacobian (CCSD)
!!
!!    Routines debug analytical Jacobian by comparing to 
!!    Jacobian computed by numerical differentiation of 
!!    Ω 
!! 
!!       A_μν = dΩ_μ/dt_ν =  d(< μ | exp(-T) H exp(T) | R >)/dt_ν
!!
!!    The numerical differentiation is done using
!!    central differences.
!!
!
   implicit none
!
!
contains
!
!
!
   module subroutine normalization_for_jacobian_debug_ccsd(wf, A_numerical_mu_nu, nu)
!!
!!    Normalization for jacobian debug
!!    Written by Sarai D. Folkestad and Tor S. Haugland, Sep 2019
!!
!!    Differentiation wrt. doubles amplitudes (nu > wf%n_t1)
!!    will yield factor one half on the diagonal nu = aiai
!!
!!    A_numerical_μ_aiai = dΩ_μ/dt_aiai = 1/2 <μ|[H^T, tau_aiai]|R>
!!
!!    This routine scales these diagonal elements by factor two
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: A_numerical_mu_nu
!
      integer, intent(in) :: nu
!
      logical :: scale_
!
      integer :: ai, aibj, aiai
!
      if (nu .le. wf%n_t1) return
!
      scale_ = .false. 
!
      aibj = nu - wf%n_t1
!
      do ai = 1, wf%n_t1
!
         aiai = ai*(ai-3)/2 + 2*ai
         if (aiai == aibj) scale_ = .true.
!
      enddo
!
      if (scale_) call dscal(wf%n_es_amplitudes, two, A_numerical_mu_nu, 1)
!
   end subroutine normalization_for_jacobian_debug_ccsd
!

!
end submodule debug_jacobian_ccsd