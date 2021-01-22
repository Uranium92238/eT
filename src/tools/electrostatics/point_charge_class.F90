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
!  but WiTHOUT ANY WArrANTY; without even the implied warranty of
!  MErCHANTABiLiTY or FiTNESS FOr A PArTiCULAr PUrPOSE. See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program. if not, see <https://www.gnu.org/licenses/>.
!
!
module point_charge_class
!
!!
!!    Point charge class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Oct 2020
!!
!
   use kinds
   use parameters
   use memory_manager_class, only: mem
!
   implicit none 
!
   type :: point_charge
!
      real(dp), dimension(3), private :: r
      real(dp),               private :: q
!
   contains 
!
      procedure :: set_r &
                => set_r_point_charge
!
      procedure :: set_q &
                => set_q_point_charge
!
      procedure :: get_r &
                => get_r_point_charge
!
      procedure :: get_q &
                => get_q_point_charge
!
   end type point_charge
!
!
   interface point_charge
!
      procedure :: new_point_charge 
!
   end interface point_charge
!
contains  
!
!
   pure function new_point_charge(r, q) result(this)
!!
!!    New point charge
!!    Written by Sarai D. Folkestad, 2020
!!
      implicit none
!
      real(dp), dimension(3), intent(in)  :: r
      real(dp),               intent(in)  :: q
      type(point_charge)                  :: this
!
      this%r = r
      this%q = q
!
   end function new_point_charge
!
!
   pure subroutine set_r_point_charge(this, r)
!!
!!    Set r
!!    Written by Sarai D. Folkestad, 2020
!!
      implicit none
!
      class(point_charge), intent(inout)  :: this
      real(dp), dimension(3), intent(in)  :: r
!
      this%r = r
!
   end subroutine set_r_point_charge
!
!
   pure subroutine set_q_point_charge(this, q)
!!
!!    Set q
!!    Written by Sarai D. Folkestad, 2020
!!
      implicit none
!
      class(point_charge), intent(inout)  :: this
      real(dp),            intent(in)     :: q
!
      this%q = q
!
   end subroutine set_q_point_charge
!
!
   pure function get_r_point_charge(this) result(r)
!!
!!    Get r
!!    Written by Sarai D. Folkestad, 2020
!!
      implicit none
!
      class(point_charge), intent(in)  :: this
      real(dp), dimension(3)           :: r
!
      r = this%r
!
   end function get_r_point_charge
!
!
   pure function get_q_point_charge(this) result(q)
!!
!!    Get q
!!    Written by Sarai D. Folkestad, 2020
!!
      implicit none
!
      class(point_charge), intent(in)  :: this
      real(dp)                         :: q
!
      q = this%q
!
   end function get_q_point_charge
!
!
end module point_charge_class
