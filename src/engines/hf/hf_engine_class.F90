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
module hf_engine_class
!
!!
!! Hartree-Fock engine class
!! Written by Sarai D. Folkestad, Eirik F. Kjønstad, and Alexander C. Paul, 2018-2022
!!
!
   use hf_class, only: hf
!
   implicit none
!
   type, abstract :: hf_engine
!
   contains
!
      procedure(ignite), deferred, public :: ignite
!
   end type hf_engine
!
!
   abstract interface
!
      subroutine ignite(this, wf)
!
         import :: hf, hf_engine
!
         implicit none
!
         class(hf_engine), intent(inout) :: this
!
         class(hf), intent(inout) :: wf
!
      end subroutine ignite
!
   end interface
!
!
end module hf_engine_class
