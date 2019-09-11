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
   module subroutine omega_for_jacobian_debug_ccs(wf, omega, t)
!!
!!    Omega for jacoboan debug
!!    Written by Sarai D. Folkestad
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: omega
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: t
!
   end subroutine omega_for_jacobian_debug_ccs
!
!
   module subroutine amplitudes_for_jacobian_debug_ccs(wf, t)
!!
!!    Amplitudes for jacobian debug
!!    Written by Sarai D. Folkestad
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: t
!
   end subroutine amplitudes_for_jacobian_debug_ccs
!
!
   module subroutine numerical_test_jacobian_ccs(wf, store_on_file)
!!
!!    Numerical test for Jacobian
!!    Written by Sarai D. Folkestad, Sep 2019
!!
      implicit none
!
      class(ccs) :: wf
!
      logical, optional :: store_on_file
!
   end subroutine numerical_test_jacobian_ccs
!
!
   module subroutine normalization_for_jacobian_debug_ccs(wf, A_numerical_mu_nu, nu)
!!
!!    Normalization for jacobian debug
!!    Written by Sarai D. Folkestad
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: A_numerical_mu_nu
!
      integer, intent(in) :: nu
!
   end subroutine normalization_for_jacobian_debug_ccs
!
!