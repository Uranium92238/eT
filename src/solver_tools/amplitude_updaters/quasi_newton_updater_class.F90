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
module quasi_newton_updater_class
!!
!!    Quasi-Newton updater class module
!!    Written by Eirik F. Kjønstad, 2020
!! 
!!    Updates the amplitudes according to a (zeroth-order) quasi-Newton estimate,
!!
!!       amplitudes_mu = amplitudes_mu - epsilon_mu^-1 residual_mu,
!!
!!    where epsilon_mu denotes orbital energy differences.
!!
   use kinds 
   use amplitude_updater_class,  only: amplitude_updater
   use ccs_class,                only: ccs
   use memory_manager_class,     only: mem
!
   implicit none
!
   type, extends(amplitude_updater) :: quasi_newton_updater
!
   contains
!
      procedure :: calculate_update => calculate_update_quasi_newton_updater
!
   end type quasi_newton_updater
!
!
   interface quasi_newton_updater
!
      procedure :: new_quasi_newton_updater
!
   end interface quasi_newton_updater
!
!
contains 
!
!
   pure function new_quasi_newton_updater(n_amplitudes,     &
                                          scale_amplitudes, &
                                          scale_residual) result(this)
!!
!!    New Quasi-Newton updater
!!    Written by Eirik F. Kjønstad, 2020
!!
!!    scale_amplitudes / scale_residual: if true, scales the double amplitudes by a factor of two
!!
      implicit none 
!
      integer, intent(in) :: n_amplitudes 
      logical, intent(in) :: scale_amplitudes 
      logical, intent(in) :: scale_residual
!
      type(quasi_newton_updater) :: this 
!
      this%n_amplitudes     = n_amplitudes
      this%scale_amplitudes = scale_amplitudes
      this%scale_residual   = scale_residual
!
   end function new_quasi_newton_updater
!
!
   subroutine calculate_update_quasi_newton_updater(this, wf, residual, update)
!!
!!    Calculate update
!!    Written by Eirik F. Kjønstad, 2020
!!
!!    Calculates the amplitude update based on the 
!!    standard quasi-Newton orbital-differences estimate:
!!
!!       dt_mu = - residual_mu/epsilon_mu
!!
!!    In order to be consistent with earlier implementation, this routine
!!    also preconditions the residual - which is not strictly necessary but
!!    has been kept for now. Disabling it will require regeneration of tests.
!!
      class(quasi_newton_updater), intent(in) :: this  
!
      class(ccs), intent(inout) :: wf 
!
      real(dp), dimension(this%n_amplitudes), intent(inout) :: residual
      real(dp), dimension(this%n_amplitudes), intent(inout) :: update
!
      integer :: k
!
      call wf%get_orbital_differences(update, this%n_amplitudes)
!
!$omp parallel do private(k)
      do k = 1, this%n_amplitudes
!
         residual(k) = -residual(k)/update(k)
         update(k) = residual(k)
!
      enddo
!$omp end parallel do
!
   end subroutine calculate_update_quasi_newton_updater
!
!
end module quasi_newton_updater_class
