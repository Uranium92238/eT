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
module mm_atom_class
!
!!
!!    MM atom class module
!!    Written by Tommaso Giovannini, Eirik F. Kjønstad and Sarai D. Folkestad, 2018 and 2020
!!
!
   use parameters
   use point_charge_class, only: point_charge
!
   implicit none
!
   type, extends(point_charge) :: mm_atom
!
!     Electronegativity, chemical hardness and periodic table symbol
!
      real(dp), private :: chi   
      real(dp), private :: eta
      character(len=2)  :: symbol
!
   contains
!
      procedure :: set_chi &
                => set_chi_mm_atom
!
      procedure :: set_eta &
                => set_eta_mm_atom
!
      procedure :: get_chi &
                => get_chi_mm_atom
!
      procedure :: get_eta &
                => get_eta_mm_atom
!
   end type  mm_atom
!
   interface  mm_atom
!
      procedure :: new_mm_atom
!
   end interface  mm_atom
!
!
contains
!
!
   pure function new_mm_atom(r, symbol) result(atom)
!!
!!    New mm_atom 
!!    Written by Sarai D. Folkestad, 2020
!!
      implicit none 
!
      type(mm_atom) :: atom
!
      real(dp), dimension(3), intent(in) :: r
      character(len=2),       intent(in) :: symbol 
!
      call atom%set_r(r)
      atom%symbol = symbol
!
   end function new_mm_atom
!
!
   pure subroutine set_chi_mm_atom(this, chi)
!!
!!    Set chi
!!    Written by Sarai D. Folkestad, 2020
!!
      implicit none
!
      class(mm_atom), intent(inout)  :: this
      real(dp),       intent(in)     :: chi
!
      this%chi = chi
!
   end subroutine set_chi_mm_atom
!
!
   pure subroutine set_eta_mm_atom(this, eta)
!!
!!    Set eta
!!    Written by Sarai D. Folkestad, 2020
!!
      implicit none
!
      class(mm_atom), intent(inout)  :: this
      real(dp),       intent(in)     :: eta
!
      this%eta = eta
!
   end subroutine set_eta_mm_atom
!
!
   pure function get_chi_mm_atom(this) result(chi)
!!
!!    Get chi
!!    Written by Sarai D. Fochilkestad, 2020
!!
      implicit none
!
      class(mm_atom), intent(in) :: this
      real(dp)                   :: chi
!
      chi = this%chi
!
   end function get_chi_mm_atom
!
!
   pure function get_eta_mm_atom(this) result(eta)
!!
!!    Get eta
!!    Written by Sarai D. Folkestad, 2020
!!
      implicit none
!
      class(mm_atom), intent(in)  :: this
      real(dp)                    :: eta
!
      eta = this%eta
!
   end function get_eta_mm_atom
!
end module mm_atom_class
