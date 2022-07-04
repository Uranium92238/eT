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
module K_MO_screener_class
!!
!! K MO screener class
!! Written by Sarai D. Folkestad, 2021
!!
   use kinds
   use ao_tool_class,             only: ao_tool
   use abstract_G_screener_class, only: abstract_G_screener
   use memory_manager_class, only: mem
!
   implicit none
!
   type, extends(abstract_G_screener) :: K_MO_screener
!
      real(dp) :: K_threshold
!
      real(dp), dimension(:,:), pointer      :: C
      real(dp), dimension(:),   allocatable  :: sh_max_C
      real(dp)                               :: max_C
!
      integer, dimension(:), allocatable :: sig_s1s2
      integer, dimension(:), allocatable :: sig_s3
      logical, dimension(:), allocatable :: s1_is_sig
!
      integer :: n_sig_s3, n_sig_shp
      integer :: n_mo, n_ao
!
   contains
!
!     Deferred routines
!
      procedure :: s1s2s3_cycle &
                => s1s2s3_cycle_K_MO_screener
!
      procedure :: s1s2s3s4_cycle &
                => s1s2s3s4_cycle_K_MO_screener
!
      procedure :: cleanup &
                => cleanup_K_MO_screener
!
      procedure :: get_D_max &
                => get_D_max_K_MO_screener
!
!     Overwritten routines
!
      procedure :: s1s2s3_exit &
                => s1s2s3_exit_K_MO_screener
!
      procedure :: prescreening &
                => prescreening_K_MO_screener
!
      procedure :: get_sig_s1s2 &
                => get_sig_s1s2_K_MO_screener
!
      procedure :: get_sig_s3 &
                => get_sig_s3_K_MO_screener
!
      procedure :: get_n_sig_s3 &
                => get_n_sig_s3_K_MO_screener
!
!     Private routines
!
      procedure, private :: cleanup_sh_lists
      procedure, private :: cleanup_C_lists
      procedure, private :: make_C_lists
      procedure, private :: determine_sig_s1
      procedure, private :: determine_sig_s3

!
   end type K_MO_screener
!
   interface  K_MO_screener
      procedure :: new_K_MO_screener
   end interface  K_MO_screener
!
contains
!
!
   function new_K_MO_screener(K_threshold, n_ao, n_mo, C) result(this)
!!
!!    New
!!    Written by Sarai D. Folkestad, Jun 2021
!!
      implicit none
!
      type(K_MO_screener) :: this
!
      real(dp),                        intent(in)         :: K_threshold
      integer,                         intent(in)         :: n_ao, n_mo
      real(dp), dimension(n_ao, n_mo), intent(in), target :: C
!
      this%K_threshold = K_threshold
!
      this%n_mo = n_mo
      this%n_ao = n_ao
!
      this%C(1:this%n_ao, 1:this%n_mo) => C(:,:)
!
   end function new_K_MO_screener
!
!
   pure function s1s2s3_cycle_K_MO_screener(this, ao, s1, s2, s3, s1s2) result(cycle_)
!!
!!    s1s2s3 cycle
!!    Written by Eirik F. Kjønstad, Linda Goletto, and Sarai D. Folkestad 2020-2021
!!
      implicit none
!
      class(K_MO_screener),intent(in) :: this
      class(ao_tool),      intent(in) :: ao
      integer,             intent(in) :: s1, s2, s3, s1s2
!
      logical :: cycle_
      real(dp) :: D_max_shp_s3s2, D_max_shp_s3s1, cs_eri_max_s1s2s3
!
      real(dp) :: density_K
!
      cycle_ = .false.
!
      D_max_shp_s3s2 = this%shp_max_D(s3, s2)
      D_max_shp_s3s1 = this%shp_max_D(s3, s1)
!
      cs_eri_max_s1s2s3 = ao%cs_eri_max(s1s2, 1) * ao%cs_eri_max_sh(s3)
!
      density_K = max(D_max_shp_s3s2 * this%sh_max_C(s1) * this%max_C,  &
                      D_max_shp_s3s1 * this%sh_max_C(s2) * this%max_C,  &
                      this%sh_max_D(s2) * this%sh_max_C(s1) * this%sh_max_C(s3), &
                      this%sh_max_D(s1) * this%sh_max_C(s2) * this%sh_max_C(s3))
!
      if (density_K * cs_eri_max_s1s2s3 < this%K_threshold) cycle_ = .true.
!
   end function s1s2s3_cycle_K_MO_screener
!
!
   pure function s1s2s3s4_cycle_K_MO_screener(this, ao, s1, s2, s3, s4, s1s2, s3s4) result(cycle_)
!!
!!    s1s2s3s4 cycle
!!    Written by Eirik F. Kjønstad, Linda Goletto, and Sarai D. Folkestad 2020-2021
!!
      implicit none
!
      class(K_MO_screener), intent(in) :: this
      class(ao_tool),       intent(in) :: ao
      integer,              intent(in) :: s1, s2, s3, s4, s1s2, s3s4
!
      logical :: cycle_
!
      real(dp) :: density_K, cs_eri_max_s1s2s3s4, D_max_shp_s3s2, D_max_shp_s3s1
!
      cycle_ = .false.
!
      D_max_shp_s3s2 = this%shp_max_D(s3, s2)
      D_max_shp_s3s1 = this%shp_max_D(s3, s1)
!
      cs_eri_max_s1s2s3s4 = ao%cs_eri_max(s1s2, 1) * ao%cs_eri_max(s3s4, 2)
!
      density_K = max(D_max_shp_s3s2 * this%sh_max_C(s1) * this%sh_max_C(s4),         &
                      D_max_shp_s3s1 * this%sh_max_C(s2) * this%sh_max_C(s4),         &
                      this%shp_max_D(s4,s2) * this%sh_max_C(s1) * this%sh_max_C(s3),  &
                      this%shp_max_D(s1,s4) * this%sh_max_C(s2) * this%sh_max_C(s3))
!
      if (density_K * cs_eri_max_s1s2s3s4 < this%K_threshold) cycle_ = .true.
!
   end function s1s2s3s4_cycle_K_MO_screener
!
!
   pure function get_D_max_K_MO_screener(this, s1, s2, s3, s4) result(D_max)
!
      implicit none
!
      class(K_MO_screener), intent(in)  :: this
      real(dp) :: D_max
      integer, intent(in) :: s1, s2, s3, s4
!
      real(dp) :: D_max_shp_s3s2, D_max_shp_s3s1
!
      D_max_shp_s3s2 = this%shp_max_D(s3, s2)
      D_max_shp_s3s1 = this%shp_max_D(s3, s1)
!
      D_max = max(D_max_shp_s3s2 * this%sh_max_C(s1) * this%sh_max_C(s4),         &
                  D_max_shp_s3s1 * this%sh_max_C(s2) * this%sh_max_C(s4),         &
                  this%shp_max_D(s4,s2) * this%sh_max_C(s1) * this%sh_max_C(s3),  &
                  this%shp_max_D(s1,s4) * this%sh_max_C(s2) * this%sh_max_C(s3))
!
   end function get_D_max_K_MO_screener
!
!
   subroutine cleanup_K_MO_screener(this)
!!
!!    Cleanup
!!    Written by Sarai D. Folkestad, 2021
!!
      implicit none
      class(K_MO_screener), intent(inout)  :: this
!
      call this%cleanup_D_lists()
      call this%cleanup_C_lists()
      call this%cleanup_sh_lists()
!
      this%C => null()
!
   end subroutine cleanup_K_MO_screener
!
!
   subroutine prescreening_K_MO_screener(this, ao, D)
!!
!!    Prescreening
!!    Written by Sarai D. Folkestad 2021
!!
      implicit none
!
      class(K_MO_screener),            intent(inout)  :: this
      class(ao_tool),                  intent(in)     :: ao
      real(dp), dimension(ao%n, ao%n), intent(in)     :: D
!
      this%n_sh = ao%n_sh
      this%n_sig_shp = ao%n_sig_eri_shp
!
      call this%make_D_lists(ao, D)
      call this%make_C_lists(ao)
      call this%determine_sig_s1(ao)
      call this%determine_sig_s3(ao)
!
   end subroutine prescreening_K_MO_screener
!
!
   pure function get_sig_s3_K_MO_screener(this, s3_tilde, s1) result(s3)
!!
!!    Get significant s3
!!    Written by Sarai D. Folkestad 2021
!!
      use warning_suppressor, only: do_nothing
      implicit none

      class(K_MO_screener), intent(in)  :: this
      integer,              intent(in) :: s3_tilde
      integer,              intent(in) :: s1
!
      integer :: s3
!
      call do_nothing(s1)
      s3 = this%sig_s3(s3_tilde)
!
   end function get_sig_s3_K_MO_screener
!
!
   pure function get_sig_s1s2_K_MO_screener(this, s1s2_tilde) result(s1s2)
!!
!!    Get significant s1s2
!!    Written by Sarai D. Folkestad 2021
!!
      implicit none

      class(K_MO_screener), intent(in)  :: this
      integer,              intent(in) :: s1s2_tilde
!
      integer :: s1s2
!
      s1s2 = this%sig_s1s2(s1s2_tilde)
!
   end function get_sig_s1s2_K_MO_screener
!
!
   pure function get_n_sig_s3_K_MO_screener(this, s1) result(n_sig_s3)
!!
!!    Get n significant s3
!!    Written by Sarai D. Folkestad 2021
!!
      use warning_suppressor, only: do_nothing
      implicit none
!
      class(K_MO_screener), intent(in) :: this
      integer,              intent(in) :: s1
!
      integer :: n_sig_s3
!
      call do_nothing(s1)
      n_sig_s3 = this%n_sig_s3
!
   end function get_n_sig_s3_K_MO_screener
!
!
   subroutine determine_sig_s3(this, ao)
!!
!!    Determine significant s3
!!    Written by Sarai D. Folkestad 2021
!
      use array_utilities, only: quicksort_ascending_int
!
      implicit none
!
      class(K_MO_screener), intent(inout) :: this
      class(ao_tool),             intent(in)    :: ao
!
      integer  :: s3, count
!
!     Count the number of significant s1 (i.e. also the significant s3)
!
      this%n_sig_s3 = 0
      do s3 = 1, ao%n_sh
!
         if (this%s1_is_sig(s3)) this%n_sig_s3 = this%n_sig_s3 + 1
!
      enddo
!
      call mem%alloc(this%sig_s3, this%n_sig_s3)
!
!     Sort significant s3 according to max_p C_w,p * ao%cs_eri_max_sh(s3), for AO w in s3.
!
      count = 0
      do s3 = 1, ao%n_sh
         if (this%s1_is_sig(s3)) then
!
            count = count + 1
            this%sig_s3(count) = s3
!
         endif
      enddo
!
      call quicksort_ascending_int(this%sig_s3, this%n_sig_s3)
!
   end subroutine determine_sig_s3
!
!
   subroutine determine_sig_s1(this, ao)
!!
!!    Determine significant s3
!!    Written by Eirik F. Kjønstad and Linda Goletto, Aug 2020
!!
      implicit none
!
      class(K_MO_screener), intent(inout) :: this
      class(ao_tool),       intent(in)    :: ao
!
      integer :: s1, s2, s1s2, s1s2_packed
      real(dp) :: cs_eri_max_s1s2, density_K
!
!     Compute the number of significant shell pairs given the current density
!     and save their indices
!
      call mem%alloc(this%s1_is_sig, this%n_sh)
      this%s1_is_sig = .false.
!
      call mem%alloc(this%sig_s1s2, this%n_sig_shp)
      this%sig_s1s2 = 0
!
      this%n_sig_s1s2 = 0
!
      do s1s2 = 1, ao%n_sig_eri_shp
!
         s1s2_packed = ao%cs_eri_max_indices(s1s2, 3)
!
         s1 = ao%cs_eri_max_indices(s1s2_packed, 1)
         s2 = ao%cs_eri_max_indices(s1s2_packed, 2)
!
         cs_eri_max_s1s2 = ao%cs_eri_max(s1s2, 1)*ao%g_max
!
         density_K = max(this%sh_max_D(s2)*this%sh_max_C(s1)*this%max_C, &
                         this%sh_max_D(s1)*this%sh_max_C(s2)*this%max_C)
!
         if (density_K * cs_eri_max_s1s2 < this%K_threshold) cycle
!
         this%s1_is_sig(s1) = .true.
!
         this%n_sig_s1s2 = this%n_sig_s1s2 + 1
         this%sig_s1s2(this%n_sig_s1s2) = s1s2
!
      enddo
!
   end subroutine determine_sig_s1
!
!
   pure function s1s2s3_exit_K_MO_screener(this, ao, s1, s2, s3, s1s2) result(exit_)
!!
!!    s1s2s3 exit
!!    Written by Eirik F. Kjønstad, Linda Goletto, and Sarai D. Folkestad 2020-2021
!!
!
      use warning_suppressor, only: do_nothing
!
      implicit none
!
      class(K_MO_screener),   intent(in) :: this
      class(ao_tool),         intent(in) :: ao
      integer,                intent(in) :: s1, s2, s3, s1s2
!
      logical :: exit_
!
      call do_nothing(this)
      call do_nothing(ao)
      call do_nothing(s2)
      call do_nothing(s1s2)
!
      exit_ = (s3 .gt. s1)
!
   end function s1s2s3_exit_K_MO_screener
!
!
   subroutine make_C_lists(this, ao)
!!
!!    Make C lists
!!    Written by Eirik F. Kjønstad, Linda Goletto, and Sarai D. Folkestad, 2020-2021
!!
      implicit none
!
      class(K_MO_screener), intent(inout) :: this
      class(ao_tool),       intent(in)    :: ao
!
      call mem%alloc(this%sh_max_C, this%n_sh)
!
      call ao%construct_sh_X_max_index_1(this%n_mo, this%C, this%sh_max_C)
      this%max_C = maxval(this%sh_max_C)
!
   end subroutine make_C_lists
!
!
   subroutine cleanup_sh_lists(this)
!!
!!    Cleanup shell lists
!!    Written by Sarai D. Folkestad, 2021
!!
      implicit none
!
      class(K_MO_screener), intent(inout) :: this
!
      if (allocated(this%s1_is_sig))call mem%dealloc(this%s1_is_sig, this%n_sh)
      if (allocated(this%sig_s1s2))call mem%dealloc(this%sig_s1s2,this%n_sig_shp)
      if (allocated(this%sig_s3))call mem%dealloc(this%sig_s3, this%n_sig_s3)
!
   end subroutine cleanup_sh_lists
!
!
   subroutine cleanup_C_lists(this)
!!
!!    Cleanup C lists
!!    Written by Sarai D. Folkestad, 2021
!!
      implicit none
!
      class(K_MO_screener), intent(inout) :: this
!
      if (allocated(this%sh_max_C)) call mem%dealloc(this%sh_max_C, this%n_sh)
!
   end subroutine cleanup_C_lists
!
!
end module K_MO_screener_class
