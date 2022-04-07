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
module cc_F_transformation_class
!
!!
!! CC F transformation (abstract)
!! Written by Eirik F. Kjønstad, Jan 2022
!!
!
   use parameters
   use ccs_class, only: ccs
!
   implicit none
!
   type, abstract :: cc_F_transformation
!
      class(ccs), pointer :: wf
!
   contains
!
      procedure(transform), public, deferred :: transform
!
   end type cc_F_transformation
!
!
   abstract interface
!
      subroutine transform(this, t, Ft)
!
         import :: cc_F_transformation, dp
!
         implicit none
!
         class(cc_F_transformation) :: this
!
         real(dp), dimension(this%wf%n_es_amplitudes), intent(in)  :: t
         real(dp), dimension(this%wf%n_es_amplitudes), intent(out) :: Ft
!
      end subroutine transform
!
   end interface
!
!
end module cc_F_transformation_class
