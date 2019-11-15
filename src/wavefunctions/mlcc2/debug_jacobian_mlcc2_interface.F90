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
   module subroutine omega_for_jacobian_debug_mlcc2(wf, omega, t)
!!
!!    Omega for Jacobian debug
!!    Written by Sarai D. Folkestad, Sep. 2019
!!
!!    Calculates omega 
!!    with dimension n_es_amplitudes for
!!    t given on input. For methods where 
!!
!!       n_es_amplitudes = n_gs_amplitudes
!!
!!    this simply entails calling construct omega.
!!    For e.g. CC2 this routine must be overwritten
!!    to obtain Ω_μ2
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: omega
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: t
!
   end subroutine omega_for_jacobian_debug_mlcc2
!
!
   module subroutine amplitudes_for_jacobian_debug_mlcc2(wf, t)
!!
!!    Amplitudes for Jacobian debug
!!    Written by Sarai D. Folkestad, Sep 2019
!!
!!    Calculates the amplitudes 
!!    with dimension n_es_amplitudes
!!
!!    For methods where 
!!
!!       n_es_amplitudes = n_gs_amplitudes
!!
!!    this simply entails calling get_amplitudes.
!!    For e.g. CC2 this routine must be overwritten
!!    to obtain t_μ2
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: t
!
   end subroutine amplitudes_for_jacobian_debug_mlcc2
!
!
   module subroutine normalization_for_jacobian_debug_mlcc2(wf, A_numerical_mu_nu, nu)
!!
!!    Normalization for Jacobian debug
!!    Written by Sarai D. Folkestad, Sep 2019
!!
!!    Routine does nothing for CCS
!!
      implicit none
!
      class(mlcc2), intent(in) :: wf
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: A_numerical_mu_nu
      integer, intent(in) :: nu
!
   end subroutine normalization_for_jacobian_debug_mlcc2
