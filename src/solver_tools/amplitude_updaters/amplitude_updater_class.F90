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
module amplitude_updater_class
!!
!!    Amplitude updator class module
!!    Written by Eirik F. Kjønstad, 2020
!! 
!!    Abstract class defining the interface for the task 
!!    of updating amplitudes in various solvers.
!!
   use parameters
   use ccs_class, only: ccs
   use memory_manager_class, only: mem
!
   implicit none
!
   type, abstract :: amplitude_updater
!
      integer :: n_amplitudes
!
      logical :: scale_amplitudes ! Scales double amplitudes to be consistent with the scaling 
                                  ! applied in the wavefunction. This is needed in the Newton-
                                  ! Raphson updater (a descendant), where the amplitude update (dt)
                                  ! is non-scaled due to the basis in the Jacobian and omega.
                                  ! It must then be scaled when added to the current amplitudes 
                                  ! (t = t + dt).
!
      logical :: scale_residual   ! Residual also needs to be scaled in some cases. This is not 
                                  ! strictly necessary, but is kept as an option in order to not
                                  ! change the original eT implementation of CC ground state DIIS.
!
   contains
!
      procedure :: update => update_amplitude_updater
!
      procedure(calculate_update), deferred :: calculate_update
!
   end type amplitude_updater 
!
!
   abstract interface
!
      subroutine calculate_update(this, wf, residual, update)
!!
!!       Calculate update
!!
!!       Interface for routine that calculates the update of the amplitudes 
!!       (update) based on the wave function and the given residual vector.
!!
         import :: ccs, dp, amplitude_updater
!
         implicit none
!
         class(amplitude_updater), intent(inout) :: this 
!
         class(ccs), intent(inout) :: wf 
!
         real(dp), dimension(this%n_amplitudes), intent(inout) :: residual 
         real(dp), dimension(this%n_amplitudes), intent(inout) :: update 
!
      end subroutine calculate_update
!
   end interface
!
!
contains 
!
!
   subroutine update_amplitude_updater(this, wf, amplitudes, residual)
!!
!!    Update
!!    Written by Eirik F. Kjønstad, 2020
!!
!!    Calculates the amplitude update and adds it to 'amplitudes'.
!!
!!    Might also modify the residual (e.g. precondition it) depending on 
!!    the implementation of 'calculate_update' and the value of 'scale_residual'.
!!
      implicit none 
!
      class(amplitude_updater), intent(inout) :: this  
!
      class(ccs), intent(inout) :: wf 
!
      real(dp), dimension(this%n_amplitudes), intent(inout) :: amplitudes 
      real(dp), dimension(this%n_amplitudes), intent(inout) :: residual 
!
      real(dp), dimension(:), allocatable :: update 
!
      call mem%alloc(update, this%n_amplitudes)
!
      call this%calculate_update(wf, residual, update)
!
      if (this%scale_amplitudes) call wf%scale_amplitudes(update) ! 'Update' and 'amplitude' needs
                                                                  ! to have the same diagonal 
                                                                  ! scaling in the daxpy below.
!
      call daxpy(this%n_amplitudes, one, update, 1, amplitudes, 1)
!
      call mem%dealloc(update, this%n_amplitudes)
!
      if (this%scale_residual) call wf%scale_amplitudes(residual)
!
   end subroutine update_amplitude_updater
!
!
end module amplitude_updater_class
