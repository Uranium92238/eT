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
module cc_wavefunctions_class
!
!!
!!    CC wavefunctions struct
!!    Written by  Sarai D. Folkestad, Feb 2022
!!
!
   use global_out, only: output
!
   implicit none
!
   type :: cc_wavefunctions
!
      character(len=50), dimension(9) :: names
      logical, dimension(9) :: allowed
!
   contains
!
      procedure, public :: set
      procedure, public :: is_allowed
!
   end type cc_wavefunctions
!
   interface cc_wavefunctions
!
      procedure :: new_cc_wavefunctions
!
   end interface cc_wavefunctions
!
contains
!
   function new_cc_wavefunctions() result(this)
!!
!!    New
!!    Written by  Sarai D. Folkestad, Feb 2022
!!
      implicit none
!
      type(cc_wavefunctions) :: this
!
      this%names = [character(len=50) :: &
                  'mp2', &
                  'ccs', &
                  'cc2', &
                  'low memory cc2', &
                  'ccsd', &
                  'ccsd(t)', &
                  'cc3', &
                  'mlcc2', &
                  'mlccsd']
!
      this%allowed = .false.
!
   end function new_cc_wavefunctions
!
!
   subroutine set(this, name_, allowed)
!!
!!    Set
!!    Written by  Sarai D. Folkestad, Feb 2022
!!
      implicit none
!
      class(cc_wavefunctions), intent(inout) :: this
      character(len=*),        intent(in)    :: name_
      logical,                 intent(in)    :: allowed
!
      integer :: i
!
      do i = 1, len(this%names)
!
         if (trim(name_) == trim(this%names(i))) then
!
            this%allowed(i) = allowed
            return
!
         endif
!
      enddo
!
      call output%error_msg('did not recognize the CC wave function!')
!
   end subroutine set
!
   function is_allowed(this, name_) result(allowed)
!!
!!    Is allowed
!!    Written by  Sarai D. Folkestad, Feb 2022
!!
      implicit none

      class(cc_wavefunctions), intent(in) :: this
      character(len=*),        intent(in) :: name_
!
      logical :: allowed
!
      integer :: i
!
!
      allowed = .false.
      do i = 1, len(this%names)
!
         if (trim(name_) == trim(this%names(i))) then
!
            allowed = this%allowed(i)
            return
!
         endif
!
      enddo
!
      call output%error_msg('did not recognize the CC wave function!')
!
   end function is_allowed
!
end module cc_wavefunctions_class
