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
module J_screener_class
!!
!! J screener class
!! Written by Sarai D. Folkestad, 2021
!!
   use kinds
   use ao_tool_class, only: ao_tool
   use abstract_G_screener_class, only: abstract_G_screener
!
   implicit none
!
   type, extends(abstract_G_screener) :: J_screener
!
      real(dp) :: J_threshold
!
   contains
!
      procedure :: s1s2s3_cycle &
                => s1s2s3_cycle_J_screener
!
      procedure :: s1s2s3s4_cycle &
                => s1s2s3s4_cycle_J_screener
!
      procedure :: cleanup &
                => cleanup_J_screener
!
      procedure :: get_D_max &
                => get_D_max_J_screener
!
   end type J_screener
!
   interface  J_screener
!
      procedure :: new_J_screener
!
   end interface  J_screener
!
contains
!
!
   pure function new_J_screener(J_threshold) result(this)
!!
!!    New J screening
!!    Written by Sarai D. Folkestad, Jun 2021
!!
      implicit none
!
      type(J_screener) :: this
      real(dp), intent(in) :: J_threshold
!
      this%J_threshold = J_threshold
!
   end function new_J_screener
!
   pure function s1s2s3_cycle_J_screener(this, ao, s1, s2, s3, s1s2) result(cycle_)
!!
!!    s1s2s3 cycle
!!    Written by Eirik F. Kjønstad, Linda Goletto, and Sarai D. Folkestad 2020-2021
!!
      implicit none
!
      class(J_screener),   intent(in)  :: this
      class(ao_tool),      intent(in)  :: ao
      integer,             intent(in) :: s1, s2, s3, s1s2
!
      logical :: cycle_
!
      real(dp) :: density
!
      cycle_ = .false.
!
      density = max(this%sh_max_D(s3), this%shp_max_D(s1, s2))
!
      if (ao%cs_eri_max(s1s2, 1) * ao%cs_eri_max_sh(s3) * density < this%J_threshold) &
         cycle_ = .true.
!
   end function s1s2s3_cycle_J_screener
!
!
   pure function s1s2s3s4_cycle_J_screener(this, ao, s1, s2, s3, s4, s1s2, s3s4) result(cycle_)
!!
!!    s1s2s3s4 cycle
!!    Written by Eirik F. Kjønstad, Linda Goletto, and Sarai D. Folkestad 2020-2021
!!
      implicit none
!
      class(J_screener),   intent(in) :: this
      class(ao_tool),      intent(in) :: ao
      integer,             intent(in) :: s1, s2, s3, s4, s1s2, s3s4
!
      logical :: cycle_
!
      real(dp) :: density, cs_eri_max_s1s2s3s4
!
      cycle_ = .false.
!
      cs_eri_max_s1s2s3s4 = ao%cs_eri_max(s1s2, 1) * ao%cs_eri_max(s3s4, 2)

      if (cs_eri_max_s1s2s3s4 * this%max_D < this%J_threshold) then
         cycle_ = .true.
         return
      endif
!
      density = max(this%shp_max_D(s3,s4), &
                      this%shp_max_D(s1,s2))
!
      if (density * cs_eri_max_s1s2s3s4 < this%J_threshold) cycle_ = .true.
!
   end function s1s2s3s4_cycle_J_screener
!
!
   subroutine cleanup_J_screener(this)
!!
!!    Cleanup
!!    Written by Eirik F. Kjønstad, Linda Goletto, and Sarai D. Folkestad 2020-2021
!!
      implicit none
      class(J_screener), intent(inout)  :: this
!
      call this%cleanup_D_lists()
!
   end subroutine cleanup_J_screener
!
!
   pure function get_D_max_J_screener(this, s1, s2, s3, s4) result(D_max)
!
      implicit none
!
      class(J_screener), intent(in)  :: this
      real(dp) :: D_max
      integer, intent(in) :: s1, s2, s3, s4
!
!
      D_max = max(this%shp_max_D(s3,s4), &
                      this%shp_max_D(s1,s2))
!
   end function get_D_max_J_screener
!
end module J_screener_class
