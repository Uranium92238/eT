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
module K_screener_class
!!
!! K screener class
!! Written by Sarai D. Folkestad, 2021
!!
   use kinds
   use ao_tool_class,               only: ao_tool
   use memory_manager_class,        only: mem
   use abstract_G_screener_class,   only: abstract_G_screener
!
   implicit none
!
   type, extends(abstract_G_screener) :: K_screener
!
      real(dp) :: K_threshold
!
      integer :: n_s3_for_all_s1
!
      integer, dimension(:), allocatable  :: sig_s3_offsets
      integer, dimension(:), allocatable  :: sig_s3
      integer, dimension(:), allocatable  :: n_sig_s3
!
   contains
!
!     Deferred routines
!
      procedure :: s1s2s3_cycle &
                => s1s2s3_cycle_K_screener
!
      procedure :: s1s2s3s4_cycle &
                => s1s2s3s4_cycle_K_screener
!
      procedure :: cleanup &
                => cleanup_K_screener
!
      procedure :: get_D_max &
                => get_D_max_K_screener
!
!     Overwritten routines
!
      procedure :: get_sig_s3 &
                => get_sig_s3_K_screener
!
      procedure :: get_n_sig_s3 &
                => get_n_sig_s3_K_screener
!
      procedure :: prescreening &
                => prescreening_K_screener
!
!     Private routines
!
      procedure, private :: determine_sig_s3
      procedure, private :: cleanup_s3_lists
      procedure, private :: set_n_sig_s3
      procedure, private :: is_s3_sig_for_s1
      procedure, private :: set_sig_s3
!
   end type K_screener
!
   interface  K_screener
      procedure :: new_K_screener
   end interface  K_screener
!
contains
!
!
   pure function new_K_screener(K_threshold) result(this)
!!
!!    New
!!    Written by Sarai D. Folkestad, Jun 2021
!!
      implicit none
!
      type(K_screener) :: this
      real(dp), intent(in) :: K_threshold
!
      this%K_threshold = K_threshold
!
   end function new_K_screener
!
!
   pure function s1s2s3_cycle_K_screener(this, ao, s1, s2, s3, s1s2) result(cycle_)
!!
!!    s1s2s3 cycle
!!    Written by Eirik F. Kjønstad, Linda Goletto, and Sarai D. Folkestad 2020-2021
!!
      implicit none
!
      class(K_screener), intent(in)  :: this
      class(ao_tool), intent(in)  :: ao
      logical :: cycle_
      integer, intent(in) :: s1, s2, s3, s1s2
!
      real(dp) :: density_K
!
      cycle_ = .false.
!
      density_K = max(this%shp_max_D(s3, s2),  &
                      this%shp_max_D(s3, s1),  &
                      this%sh_max_D(s2),       &
                      this%sh_max_D(s1))
!
      if (ao%cs_eri_max(s1s2, 1) * ao%cs_eri_max_sh(s3) * density_K < this%K_threshold) &
         cycle_ = .true.
!
   end function s1s2s3_cycle_K_screener
!
!
   pure function s1s2s3s4_cycle_K_screener(this, ao, s1, s2, s3, s4, s1s2, s3s4) result(cycle_)
!!
!!    s1s2s3s4 cycle
!!    Written by Eirik F. Kjønstad, Linda Goletto, and Sarai D. Folkestad 2020-2021
!!
      implicit none
!
      class(K_screener), intent(in)  :: this
      class(ao_tool), intent(in)  :: ao
      logical :: cycle_
      integer, intent(in) :: s1, s2, s3, s4, s1s2, s3s4
!
      real(dp) :: density_K, cs_eri_max_s1s2s3s4
!
      cycle_ = .false.
!
      cs_eri_max_s1s2s3s4 = ao%cs_eri_max(s1s2, 1) * ao%cs_eri_max(s3s4, 2)

      if (cs_eri_max_s1s2s3s4 * this%max_D < this%K_threshold) then
         cycle_ = .true.
         return
      endif
!
      density_K = max(this%shp_max_D(s3,s2), &
                      this%shp_max_D(s3,s1), &
                      this%shp_max_D(s4,s2), &
                      this%shp_max_D(s1,s4))
!
      if (density_K * cs_eri_max_s1s2s3s4 < this%K_threshold) cycle_ = .true.
!
   end function s1s2s3s4_cycle_K_screener
!
!
   subroutine cleanup_K_screener(this)
!!
!!    Cleanup
!!    Written by Sarai D. Folkestad, 2021
!!
      implicit none
!
      class(K_screener), intent(inout)  :: this
!
      call this%cleanup_D_lists()
      call this%cleanup_s3_lists()
!
   end subroutine cleanup_K_screener
!
!
   pure function get_sig_s3_K_screener(this, s3_tilde, s1) result(s3)
!!
!!    Get significant s3
!!    Written by Sarai D. Folkestad 2021
!!
      implicit none

      class(K_screener), intent(in) :: this
      integer,           intent(in) :: s3_tilde
      integer,           intent(in) :: s1
!
      integer :: s3
!
      s3 = this%sig_s3(this%sig_s3_offsets(s1) - 1 + s3_tilde)
!
   end function get_sig_s3_K_screener
!
!
   pure function get_n_sig_s3_K_screener(this, s1) result(n_sig_s3)
!!
!!    Get n significant s3
!!    Written by Sarai D. Folkestad 2021
!!
      implicit none
!
      class(K_screener),   intent(in) :: this
      integer,             intent(in) :: s1
!
      integer :: n_sig_s3
!
      n_sig_s3 = this%n_sig_s3(s1)
!
   end function get_n_sig_s3_K_screener
!
!
   subroutine prescreening_K_screener(this, ao, D)
!!
!!    Prescreening
!!    Written by Sarai D. Folkestad 2021
!!
      implicit none
!
      class(K_screener),               intent(inout)  :: this
      class(ao_tool),                  intent(in)     :: ao
      real(dp), dimension(ao%n, ao%n), intent(in)     :: D
!
      this%n_sh = ao%n_sh
      this%n_sig_s1s2 = ao%n_sig_eri_shp
!
      call this%make_D_lists(ao, D)
      call this%determine_sig_s3(ao)
!
   end subroutine prescreening_K_screener
!
!
   subroutine cleanup_s3_lists(this)
!!
!!    Cleanup s3 lists
!!    Written by Sarai D. Folkestad, 2021
!!
      implicit none
!
      class(K_screener), intent(inout) :: this
!
      if (allocated(this%sig_s3_offsets)) &
            call mem%dealloc(this%sig_s3_offsets, this%n_sh)
!
      if (allocated(this%n_sig_s3)) &
            call mem%dealloc(this%n_sig_s3, this%n_sh)
!
      if (allocated(this%sig_s3)) &
            call mem%dealloc(this%sig_s3, this%n_s3_for_all_s1)
!
   end subroutine cleanup_s3_lists
!
!
   subroutine determine_sig_s3(this, ao)
!!
!!    Determine significant s3
!!    Written by Sarai D. Folkestad 2021
!!
      implicit none
!
      class(K_screener),   intent(inout)  :: this
      class(ao_tool),      intent(in)     :: ao
!
      logical, dimension(:,:), allocatable :: s3_is_significant
!
      call mem%alloc(this%sig_s3_offsets, ao%n_sh)
      call mem%alloc(this%n_sig_s3, ao%n_sh)
!
      call mem%alloc(s3_is_significant, ao%n_sh, ao%n_sh)
      call this%is_s3_sig_for_s1(ao, s3_is_significant)
!
      call this%set_n_sig_s3(ao, s3_is_significant)
!
      call mem%alloc(this%sig_s3, this%n_s3_for_all_s1)
!
      call this%set_sig_s3(ao, s3_is_significant)
      call mem%dealloc(s3_is_significant, ao%n_sh, ao%n_sh)
!
   end subroutine determine_sig_s3
!
!
   subroutine is_s3_sig_for_s1(this, ao, s3_is_significant)
!!
!!    Is s3 significant for s1
!!    Written by Eirik F. Kjønstad, Linda Goletto, and Sarai D. Folkestad, 2020-2021
!!
!!    Returns a logical array with elements s3_is_significant(s3, s1)
!!    which is true if s3 is significant for s1
!!
      use array_initialization, only: set_logicals
!
      implicit none
!
      class(K_screener),                     intent(in)  :: this
      class(ao_tool),                        intent(in)  :: ao
      logical, dimension(ao%n_sh, ao%n_sh),  intent(out) :: s3_is_significant
!
      logical, dimension(:,:), allocatable :: significant_s3_or_s4_for_s1
!
      integer :: s1, s, s2, s2_tilde, s3, s4, s4_tilde, s1s2
      real(dp) :: cs_eri_max_s1s2s, density_K
!
      call set_logicals(s3_is_significant, ao%n_sh**2, .false.)
      call mem%alloc(significant_s3_or_s4_for_s1, ao%n_sh, ao%n_sh, set_to=.false.)
!
!$omp parallel do private(s1, s, s2_tilde, s2, s1s2, cs_eri_max_s1s2s, density_K) schedule(dynamic)
      do s1 = 1, ao%n_sh
         do s = 1, ao%n_sh
            do s2_tilde = 1, ao%n_sig_s2_for_s1(s1)
!
               s2 = ao%sig_s2_for_s1(s2_tilde, s1)
!
               if (s2 .gt. s1) exit
!
               s1s2 = (max(s1,s2)*(max(s1,s2)-3)/2) + s1 + s2
!
               cs_eri_max_s1s2s = ao%cs_eri_max(s1s2, 2)*ao%cs_eri_max_sh(s)
!
               density_K = max(this%shp_max_D(s,s1), this%shp_max_D(s2,s))
!
               if (cs_eri_max_s1s2s * density_K > this%K_threshold) then
!
                  significant_s3_or_s4_for_s1(s, s1) = .true.
                  exit
!
               endif
            enddo
         enddo
      enddo
!$omp end parallel do
!
!
!$omp parallel do private(s1, s3, s4_tilde, s4) schedule(dynamic)
      do s1 = 1, ao%n_sh
         do s3 = 1, s1
!
            s3_is_significant(s3, s1) = significant_s3_or_s4_for_s1(s3, s1)
!
            if (s3_is_significant(s3, s1)) cycle
!
            do s4_tilde = 1, ao%n_sig_s2_for_s1(s3)
!
               s4 = ao%sig_s2_for_s1(s4_tilde, s3)
!
               if (significant_s3_or_s4_for_s1(s4, s1)) then
!
                  s3_is_significant(s3, s1) = .true.
                  exit
!
               endif
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(significant_s3_or_s4_for_s1, ao%n_sh, ao%n_sh)
!
   end subroutine is_s3_sig_for_s1
!
!
   subroutine set_sig_s3(this, ao, s3_is_significant)
!!
!!    Set significant s3
!!    Written by Eirik F. Kjønstad, Linda Goletto, and Sarai D. Folkestad, 2020-2021
!!
      implicit none
!
      class(K_screener),                    intent(inout) :: this
      class(ao_tool),                       intent(in)    :: ao
      logical, dimension(ao%n_sh, ao%n_sh), intent(in)    :: s3_is_significant
!
      integer :: s1, s3, counter
!
      do s1 = 1, ao%n_sh
!
         counter = 0
         do s3 = 1, s1
!
            if (s3_is_significant(s3,s1)) then
!
               this%sig_s3(this%sig_s3_offsets(s1) + counter) = s3
               counter = counter + 1
!
            endif
!
         enddo
      enddo
   end subroutine set_sig_s3
!
!
   subroutine set_n_sig_s3(this, ao, s3_is_significant)
!!
!!    Set n sig s3 for s1
!!    Written by Eirik F. Kjønstad, Linda Goletto, and Sarai D. Folkestad, 2020-2021
!!
      implicit none
!
      class(K_screener),                    intent(inout) :: this
      class(ao_tool),                       intent(in)    :: ao
      logical, dimension(ao%n_sh, ao%n_sh), intent(in)    :: s3_is_significant
!
      integer :: s1, s3
!
      this%n_s3_for_all_s1 = 0
      this%n_sig_s3 = 0
!
      do s1 = 1, ao%n_sh
         do s3 = 1, s1
!
            if (s3_is_significant(s3,s1)) then
!
               this%n_s3_for_all_s1 = this%n_s3_for_all_s1 + 1
               this%n_sig_s3(s1) = this%n_sig_s3(s1) + 1
!
            endif
!
         enddo
      enddo
!
      this%sig_s3_offsets = 0
      this%sig_s3_offsets(1) = 1
!
      do s1 = 2, ao%n_sh
!
         this%sig_s3_offsets(s1) = this%sig_s3_offsets(s1 - 1) + this%n_sig_s3(s1 - 1)
!
      enddo
!
   end subroutine set_n_sig_s3
!
!
   pure function get_D_max_K_screener(this, s1, s2, s3, s4) result(D_max)
!
      implicit none
!
      class(K_screener), intent(in)  :: this
      real(dp) :: D_max
      integer, intent(in) :: s1, s2, s3, s4
!
      D_max = max(this%shp_max_D(s3,s2), &
                   this%shp_max_D(s3,s1), &
                   this%shp_max_D(s4,s2), &
                   this%shp_max_D(s1,s4))
!
   end function get_D_max_K_screener
!
end module K_screener_class
