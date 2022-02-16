!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2022 the authors of eT
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
module cc_lr_eta_xi_calculator_class
!
!!
!! CC LR eta and xi calculator class
!! Written by Eirik F. Kjønstad, Jan 2022
!!
!
   use parameters
   use ccs_class, only: ccs
   use cc_eta_xi_calculator_class, only: cc_eta_xi_calculator
   use memory_manager_class, only: mem
!
   implicit none
!
   type, extends(cc_eta_xi_calculator) :: cc_lr_eta_xi_calculator
!
   contains
!
      procedure, public :: calculate &
                        => calculate_cc_lr_eta_xi_calculator
!
   end type cc_lr_eta_xi_calculator
!
!
contains
!
!
   function new_cc_lr_eta_xi_calculator(wf) result(this)
!!
!!    New cc lr eta xi calculator
!!    Written by Eirik F. Kjønstad, Jan 2022
!!
      implicit none
!
      class(ccs), target, intent(in) :: wf
!
      type(cc_lr_eta_xi_calculator) :: this
!
      this%wf => wf
!
   end function new_cc_lr_eta_xi_calculator
!
!
   subroutine calculate_cc_lr_eta_xi_calculator(this, xiX, etaX)
!!
!!    Calculate
!!    Written by Eirik F. Kjønstad, Jan 2022
!!
      implicit none
!
      class(cc_lr_eta_xi_calculator) :: this
!
      real(dp), dimension(this%wf%n_es_amplitudes, 3), intent(out) :: xiX
      real(dp), dimension(this%wf%n_es_amplitudes, 3), intent(out) :: etaX
!
      integer :: q
!
      real(dp), dimension(:,:,:), allocatable :: X
!
      call mem%alloc(X, this%wf%n_mo, this%wf%n_mo, 3)
!
      call this%wf%get_t1_oei('dipole', X)
!
      do q = 1, 3
!
         call this%wf%construct_xiX(X(:,:,q), xiX(:,q))
         call this%wf%construct_etaX(X(:,:,q), etaX(:,q))
!
      end do
!
      call mem%dealloc(X, this%wf%n_mo, this%wf%n_mo, 3)
!
   end subroutine calculate_cc_lr_eta_xi_calculator
!
!
end module cc_lr_eta_xi_calculator_class
