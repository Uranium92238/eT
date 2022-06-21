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
module abstract_G_screener_class
!!
!! Abstract G screener class
!! Written by Sarai D. Folkestad, 2021
!!
!! Handles the screening for construction of the two-electron
!! contribution to the Fock matrix (G)
!!
   use kinds
   use ao_tool_class, only: ao_tool
   use memory_manager_class, only: mem
!
   implicit none
!
   type, abstract :: abstract_G_screener
!
      real(dp), dimension(:,:), allocatable :: shp_max_D
      real(dp), dimension(:),   allocatable :: sh_max_D
!
      real(dp) :: max_D
!
      integer :: n_sh, n_sig_s1s2
!
   contains
!
      procedure (s1s2s3_cycle_abstract),     deferred :: s1s2s3_cycle
      procedure (s1s2s3s4_cycle_abstract),   deferred :: s1s2s3s4_cycle
      procedure (cleanup_abstract),          deferred :: cleanup
      procedure (get_D_max_abstract),        deferred :: get_D_max
!
      procedure, public :: s1s2s3_exit &
                        => s1s2s3_exit_abstract
!
      procedure, public :: get_sig_s3 &
                        => get_sig_s3_abstract
!
      procedure, public :: get_sig_s1s2 &
                        => get_sig_s1s2_abstract
!
      procedure, public :: get_n_sig_s3 &
                        => get_n_sig_s3_abstract
!
      procedure, public :: prescreening &
                        => prescreening_abstract
!
      procedure, public, non_overridable :: cleanup_D_lists
      procedure, public, non_overridable :: make_D_lists
!
   end type abstract_G_screener
!
   abstract interface
!
      pure function s1s2s3_cycle_abstract(this, ao, s1, s2, s3, s1s2) result(cycle_)
!
         import abstract_G_screener
         import ao_tool
!
         implicit none
!
         class(abstract_G_screener), intent(in)  :: this
         class(ao_tool), intent(in)  :: ao
         logical :: cycle_
         integer, intent(in) :: s1, s2, s3, s1s2
!
      end function s1s2s3_cycle_abstract
!
      pure function s1s2s3s4_cycle_abstract(this, ao, s1, s2, s3, s4, s1s2, s3s4) result(cycle_)
!
         import abstract_G_screener
         import ao_tool
!
         implicit none
!
         class(abstract_G_screener), intent(in)  :: this
         class(ao_tool), intent(in)  :: ao
         logical :: cycle_
         integer, intent(in) :: s1, s2, s3, s4, s1s2, s3s4
!
      end function s1s2s3s4_cycle_abstract
!
      subroutine cleanup_abstract(this)
!
         import abstract_G_screener
!
         implicit none
!
         class(abstract_G_screener), intent(inout)  :: this
!
      end subroutine cleanup_abstract
!
      pure function get_D_max_abstract(this, s1, s2, s3, s4) result(D_max)
!
         use parameters
         import abstract_G_screener
!
         implicit none
!
         class(abstract_G_screener), intent(in)  :: this
         real(dp) :: D_max
         integer, intent(in) :: s1, s2, s3, s4
!
      end function get_D_max_abstract
!
   end interface
!
contains
!
!
   subroutine make_D_lists(this, ao, D)
!!
!!    Make D lists
!!    Written by Sarai D. Folkestad, 2021
!!
!
      use array_utilities, only: get_abs_max
!
      implicit none
!
      class(abstract_G_screener),      intent(inout) :: this
      class(ao_tool),                  intent(in)    :: ao
      real(dp), dimension(ao%n, ao%n), intent(in)    :: D
!
      call mem%alloc(this%shp_max_D, this%n_sh, this%n_sh)
      call mem%alloc(this%sh_max_D, this%n_sh)
!
      call ao%construct_shp_X(D, this%shp_max_D)
      call ao%construct_sh_X_from_shp_X(this%shp_max_D, this%sh_max_D)
      this%max_D = get_abs_max(this%sh_max_D, this%n_sh)
!
   end subroutine make_D_lists
!
!
   subroutine cleanup_D_lists(this)
!!
!!    Cleanup D lists
!!    Written by Sarai D. Folkestad, 2021
!!
      implicit none
!
      class(abstract_G_screener), intent(inout)  :: this
!
      if (allocated(this%shp_max_D)) &
            call mem%dealloc(this%shp_max_D, this%n_sh, this%n_sh)
!
      if (allocated(this%sh_max_D)) &
            call mem%dealloc(this%sh_max_D, this%n_sh)
!
   end subroutine cleanup_D_lists
!
!
   subroutine prescreening_abstract(this, ao, D)
!!
!!    Prescreening
!!    Written by Sarai D. Folkestad 2021
!!
      implicit none
!
      class(abstract_G_screener),      intent(inout)  :: this
      class(ao_tool),                  intent(in)     :: ao
      real(dp), dimension(ao%n, ao%n), intent(in)     :: D
!
      this%n_sh = ao%n_sh
      this%n_sig_s1s2 = ao%n_sig_eri_shp
!
      call this%make_D_lists(ao, D)
!
   end subroutine prescreening_abstract
!
!
   pure function get_sig_s3_abstract(this, s3_tilde, s1) result(s3)
!!
!!    Get significant s3
!!    Written by Sarai D. Folkestad 2021
!!
      use warning_suppressor, only: do_nothing
      implicit none

      class(abstract_G_screener), intent(in)  :: this
      integer,                    intent(in) :: s3_tilde
      integer,                    intent(in) :: s1
!
      integer :: s3
!
      call do_nothing(s1)
      call do_nothing(this)
      s3 = s3_tilde
!
   end function get_sig_s3_abstract
!
!
   pure function get_sig_s1s2_abstract(this, s1s2_tilde) result(s1s2)
!!
!!    Get significant s1s2
!!    Written by Sarai D. Folkestad 2021
!!
      use warning_suppressor, only: do_nothing
      implicit none

      class(abstract_G_screener), intent(in)  :: this
      integer,                    intent(in) :: s1s2_tilde
!
      integer :: s1s2
!
      call do_nothing(this)
      s1s2 = s1s2_tilde
!
   end function get_sig_s1s2_abstract
!
!
   pure function get_n_sig_s3_abstract(this, s1) result(n_sig_s3)
!!
!!    Get n significant s3
!!    Written by Sarai D. Folkestad 2021
!!
      use warning_suppressor, only: do_nothing
      implicit none
!
      class(abstract_G_screener), intent(in) :: this
      integer,                    intent(in) :: s1
!
      integer :: n_sig_s3
!
      call do_nothing(this)
      n_sig_s3 = s1
!
   end function get_n_sig_s3_abstract
!
!
   pure function s1s2s3_exit_abstract(this, ao, s1, s2, s3, s1s2) result(exit_)
!!
!!    s1s2s3 exit
!!    Written by Eirik F. Kjønstad, Linda Goletto, and Sarai D. Folkestad 2020-2021
!!
!
      use warning_suppressor, only: do_nothing
!
      implicit none
!
      class(abstract_G_screener),   intent(in) :: this
      class(ao_tool),               intent(in) :: ao
      integer,                      intent(in) :: s1, s2, s3, s1s2
!
      logical :: exit_
!
      call do_nothing(this)
      call do_nothing(ao)
      call do_nothing(s1)
      call do_nothing(s2)
      call do_nothing(s3)
      call do_nothing(s1s2)
!
      exit_ = .false.
!
   end function s1s2s3_exit_abstract
!
!
end module abstract_G_screener_class
