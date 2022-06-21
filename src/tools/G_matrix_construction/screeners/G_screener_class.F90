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
module G_screener_class
!!
!! G screener class
!! Written by Sarai D. Folkestad, 2021
!!
   use kinds
   use ao_tool_class,             only: ao_tool
   use abstract_G_screener_class, only: abstract_G_screener
!
   implicit none
!
   type, extends(abstract_G_screener) :: G_screener
!
      real(dp) :: K_threshold, J_threshold
!
   contains
!
      procedure :: s1s2s3_cycle &
                => s1s2s3_cycle_G_screener
!
      procedure :: s1s2s3s4_cycle &
                => s1s2s3s4_cycle_G_screener
!
      procedure :: cleanup &
                => cleanup_G_screener
!
      procedure :: get_D_max &
                => get_D_max_G_screener
!
   end type G_screener
!
   interface  G_screener
!
      procedure :: new_G_screener
!
   end interface  G_screener
!
contains
!
!
   pure function new_G_screener(J_threshold, K_threshold) result(this)
!!
!!    New
!!    Written by Sarai D. Folkestad, Jun 2021
!!
      implicit none
!
      type(G_screener) :: this
      real(dp), intent(in) :: J_threshold, K_threshold
!
      this%J_threshold = J_threshold
      this%K_threshold = K_threshold
!
   end function new_G_screener
!
!
   pure function s1s2s3_cycle_G_screener(this, ao, s1, s2, s3, s1s2) result(cycle_)
!!
!!    s1s2s3 cycle
!!    Written by Eirik F. Kjønstad, Linda Goletto, and Sarai D. Folkestad 2020-2021
!!
      implicit none
!
      class(G_screener),   intent(in) :: this
      class(ao_tool),      intent(in) :: ao
      integer,             intent(in) :: s1, s2, s3, s1s2
!
      logical :: cycle_
!
      real(dp) :: density_J, density_K
!
      cycle_ = .false.
!
      density_J = max(this%sh_max_D(s3), this%shp_max_D(s1, s2))
!
      density_K = max(this%shp_max_D(s3, s2),  &
                      this%shp_max_D(s3, s1),  &
                      this%sh_max_D(s2),       &
                      this%sh_max_D(s1))
!
      if (ao%cs_eri_max(s1s2, 1) * ao%cs_eri_max_sh(s3) * density_J < this%J_threshold .and. &
          ao%cs_eri_max(s1s2, 1) * ao%cs_eri_max_sh(s3) * density_K < this%K_threshold) &
         cycle_ = .true.
!
   end function s1s2s3_cycle_G_screener
!
!
   pure function s1s2s3s4_cycle_G_screener(this, ao, s1, s2, s3, s4, s1s2, s3s4) result(cycle_)
!!
!!    s1s2s3s4 cycle
!!    Written by Eirik F. Kjønstad, Linda Goletto, and Sarai D. Folkestad 2020-2021
!!
      implicit none
!
      class(G_screener),   intent(in) :: this
      class(ao_tool),      intent(in) :: ao
      integer,             intent(in) :: s1, s2, s3, s4, s1s2, s3s4
!
      logical :: cycle_
!
      real(dp) :: density_J, density_K, cs_eri_max_s1s2s3s4
!
      real(dp) :: skip_threshold
!
      cycle_ = .false.
!
      skip_threshold = min(this%K_threshold, this%J_threshold)
!
      cs_eri_max_s1s2s3s4 = ao%cs_eri_max(s1s2, 1) * ao%cs_eri_max(s3s4, 2)

      if (cs_eri_max_s1s2s3s4 * this%max_D < skip_threshold) then
         cycle_ = .true.
         return
      endif
!
      density_J = max(this%shp_max_D(s3,s4), &
                      this%shp_max_D(s1,s2))
!
      density_K = max(this%shp_max_D(s3,s2), &
                      this%shp_max_D(s3,s1), &
                      this%shp_max_D(s4,s2), &
                      this%shp_max_D(s1,s4))
!
      if (density_K * cs_eri_max_s1s2s3s4 < this%K_threshold .and. &
          density_J * cs_eri_max_s1s2s3s4 < this%J_threshold) cycle_ = .true.
!
   end function s1s2s3s4_cycle_G_screener
!
!
   subroutine cleanup_G_screener(this)
!!
!!    Cleanup
!!    Written by Eirik F. Kjønstad, Linda Goletto, and Sarai D. Folkestad 2020-2021
!!
      implicit none
      class(G_screener), intent(inout)  :: this
!
      call this%cleanup_D_lists()
!
   end subroutine cleanup_G_screener
!
!
   pure function get_D_max_G_screener(this, s1, s2, s3, s4) result(D_max)
!
      implicit none
!
      class(G_screener), intent(in)  :: this
      real(dp) :: D_max
      integer, intent(in) :: s1, s2, s3, s4
!
      real(dp) :: density_J, density_K
!
!
      density_J = max(this%shp_max_D(s3,s4), &
                      this%shp_max_D(s1,s2))
!
      density_K = max(this%shp_max_D(s3,s2), &
                   this%shp_max_D(s3,s1), &
                   this%shp_max_D(s4,s2), &
                   this%shp_max_D(s1,s4))
!
      D_max = max(density_J, density_K)
!
   end function get_D_max_G_screener
!
!
end module G_screener_class
