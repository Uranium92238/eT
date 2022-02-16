!
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
module cc_engine_class
!
!!
!! Coupled cluster engine class
!! Written by Tor S. Haugland and Eirik F. Kjønstad, 2019-2021
!!
!
   use ccs_class, only: ccs
   use cc_wavefunctions_class, only: cc_wavefunctions
   use global_out, only: output
!
   implicit none
!
   type, abstract :: cc_engine
!
      class(cc_wavefunctions), allocatable :: allowed_cc_wfs
!
   contains
!
      procedure(ignite), deferred, public :: ignite
      procedure(set_allowed_wfs), deferred, public :: set_allowed_wfs
!
      procedure, public :: check_wavefunctions
!
   end type cc_engine
!
!
   abstract interface
!
      subroutine ignite(this, wf)
!
         import :: ccs, cc_engine
!
         implicit none
!
         class(cc_engine), intent(inout) :: this
!
         class(ccs), intent(inout) :: wf
!
      end subroutine ignite
!
!
      subroutine set_allowed_wfs(this)
!
         import :: cc_engine
!
         implicit none
!
         class(cc_engine), intent(inout) :: this
!
      end subroutine set_allowed_wfs
!
   end interface
!
contains
!
   subroutine check_wavefunctions(this, wf)
!!
!!    Check wavefunction
!!    Written by Sarai D. Folkestad, Feb 2022
!!
      implicit none
!
      class(cc_engine), intent(in) :: this
      class(ccs), intent(in) :: wf
!
      logical :: allowed
!
      allowed = this%allowed_cc_wfs%is_allowed(wf%name_)
!
      if (.not. allowed) &
         call output%error_msg('wave function ' // trim(wf%name_) // &
                               ' is not allowed for the selected engine!')
!
   end subroutine check_wavefunctions
!
end module cc_engine_class
