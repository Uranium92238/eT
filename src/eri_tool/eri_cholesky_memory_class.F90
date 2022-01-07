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
!
module eri_cholesky_memory_class
!
!!
!! ERI Cholesky memory class
!! Written by Sarai D. Folkestad, Sep 2021
!!
!! Handles the memory storage and delivery of the
!! Cholesky vectors
!!
   use parameters
!
   use global_out,   only : output
!
   use array_3D_list_class,      only: array_3D_list
   use range_class,              only: range_
   use memory_manager_class,     only: mem
!
   use abstract_eri_cholesky_class, only: abstract_eri_cholesky
!
   implicit none
!
   type, extends(abstract_eri_cholesky) :: eri_cholesky_memory
!
      type(array_3D_list), allocatable, private ::  L_Jpq
!
   contains
!
      procedure, private :: copy_in_block
      procedure, private :: copy_out_block
      procedure, private :: allocate_L_Jpq
!
      procedure, public :: initialize &
                        => initialize_eri_cholesky_memory
!
      procedure, public :: get &
                        => get_eri_cholesky_memory
!
      procedure, public :: set &
                        => set_eri_cholesky_memory
!
      final :: cleanup
!
      procedure, public :: get_memory_estimate &
                        => get_memory_estimate_eri_cholesky_memory
!
   end type eri_cholesky_memory
!
   interface eri_cholesky_memory
!
      procedure :: new_eri_cholesky_memory
!
   end interface eri_cholesky_memory
!
contains
!
   function new_eri_cholesky_memory() result(this)
!!
!!    New
!!    Written by Sarai D. Folkestad, Sep 2021
!!
      implicit none
!
      type(eri_cholesky_memory) :: this
!
      this%n_J       = -1
      this%n_ranges  = -1
      this%dim_      = -1
!
      this%L_Jpq = array_3D_list()
!
   end function new_eri_cholesky_memory
!
!
   subroutine initialize_eri_cholesky_memory(this, n_J, n_ranges, range_lengths)
!!
!!    Initialize
!!    Written by Sarai D. Folkestad, Sep 2021
!!
      implicit none

      class(eri_cholesky_memory), intent(inout) :: this
!
      integer, intent(in) :: n_ranges, n_J
      integer, dimension(n_ranges), intent(in) :: range_lengths
!
      call this%set_dimensions(n_J, n_ranges, range_lengths)
      call this%set_index_ranges()
!
      call this%allocate_L_Jpq()
!
   end subroutine initialize_eri_cholesky_memory
!
!
   subroutine allocate_L_Jpq(this)
!!
!!    Allocate L_Jpq
!!    Written by Sarai D. Folkestad, Sep 2021
!!
      implicit none
!
      class(eri_cholesky_memory), intent(inout) :: this
!
      integer :: p, q
!
      this%n_blocks = 0
!
      do p = 1, this%n_ranges
         do q = 1, this%n_ranges
!
            call this%L_Jpq%push_back(this%n_J, this%range_lengths(p), this%range_lengths(q))
!
            this%n_blocks = this%n_blocks + 1
!
         enddo
      enddo
!
   end subroutine allocate_L_Jpq
!
!
   subroutine set_eri_cholesky_memory(this, L_Jpq, first_p, last_p, first_q, last_q)
!!
!!    Set
!!    Written by Sarai D. Folkestad, Sep 2021
!!
      implicit none
!
      class(eri_cholesky_memory),                            intent(inout):: this
      integer,                                                       intent(in)   :: first_p, last_p
      integer,                                                       intent(in)   :: first_q, last_q
      real(dp), dimension(this%n_J, first_p:last_p, first_q:last_q), intent(in)   :: L_Jpq
!
      type(range_) :: subrange_p, subrange_q, p_range, q_range
!
      integer :: block_, p, q
!
      p_range = range_(first_p, last_p - first_p + 1)
      q_range = range_(first_q, last_q - first_q + 1)
!
      block_ = 0
!
      do p = 1, this%n_ranges
         do q = 1, this%n_ranges
!
            block_ = block_ + 1
!
            if (.not. (p_range%overlaps(this%index_ranges(p)) .and. &
                       q_range%overlaps(this%index_ranges(q)))) cycle
!
!           Get (sub)range of p and q for the current block.
            call this%get_block_range_overlap(block_,            &
                                             p_range, q_range,   &
                                             subrange_p,         &
                                             subrange_q)
!
            call this%copy_in_block(block_,               &
                                    subrange_p,           &
                                    subrange_q, L_Jpq,    &
                                    p_range, q_range)
!
         enddo
      enddo
!
   end subroutine set_eri_cholesky_memory
!
!
   subroutine get_eri_cholesky_memory(this, L_Jpq, first_p, last_p, first_q, last_q)
!!
!!    Get
!!    Written by Sarai D. Folkestad, Sep 2021
!!
      implicit none
!
      class(eri_cholesky_memory),                                    intent(inout):: this
      integer,                                                       intent(in)   :: first_p, last_p
      integer,                                                       intent(in)   :: first_q, last_q
      real(dp), dimension(this%n_J, first_p:last_p, first_q:last_q), intent(out)  :: L_Jpq
!
      type(range_) :: subrange_p, subrange_q, p_range, q_range
!
      integer :: block_, p, q
!
      p_range = range_(first_p, last_p - first_p + 1)
      q_range = range_(first_q, last_q - first_q + 1)
!
      block_ = 0
!
      do p = 1, this%n_ranges
         do q = 1, this%n_ranges
!
            block_ = block_ + 1
!
            if (.not. (p_range%overlaps(this%index_ranges(p)) .and. &
                       q_range%overlaps(this%index_ranges(q)))) cycle
!
!           Get (sub)range of p and q for the current block.
            call this%get_block_range_overlap(block_,            &
                                             p_range, q_range,   &
                                             subrange_p,         &
                                             subrange_q)
!
            call this%copy_out_block(block_,     &
                              subrange_p,        &
                              subrange_q, L_Jpq, &
                              p_range, q_range)
!
         enddo
      enddo
!
   end subroutine get_eri_cholesky_memory
!
!
   subroutine copy_in_block(this, block_, set_range_p, set_range_q, L_Jpq, p_range, q_range)
!!
!!    Copy in block_
!!    Written by Sarai D. Folkestad, Sep 2021
!!
      implicit none
!
      class(eri_cholesky_memory), intent(inout) :: this
!
      type(range_), intent(in) :: set_range_p, set_range_q
      type(range_), intent(in) :: p_range, q_range
!
      real(dp), dimension(this%n_J, p_range%length, q_range%length), intent(in) :: L_Jpq
!
      integer, intent(in) :: block_
!
      real(dp), dimension(:,:,:), pointer :: block_pointer
!
      integer :: p, q, k, r, s
      integer, dimension(2) :: r_and_s
!
      r_and_s = this%get_range_indices_from_block(block_)
      r = r_and_s(1)
      s = r_and_s(2)
!
      call this%L_Jpq%get_array_n(block_pointer, block_)
!
!$omp parallel do private(p, q, k)
      do q = set_range_q%first, set_range_q%get_last()
         do p = set_range_p%first, set_range_p%get_last()
            do k = 1, this%n_J

               block_pointer(k, p - this%index_ranges(r)%first + 1, &
                                q - this%index_ranges(s)%first + 1) = &
                       L_Jpq(k, p - p_range%first + 1, q - q_range%first + 1)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine copy_in_block
!
!
   subroutine copy_out_block(this, block_, set_range_p, set_range_q, L_Jpq, p_range, q_range)
!!
!!    Copy out block_
!!    Written by Sarai D. Folkestad, Sep 2021
!!
      implicit none
!
      class(eri_cholesky_memory), intent(inout) :: this
!
      type(range_), intent(in) :: set_range_p, set_range_q
      type(range_), intent(in) :: p_range, q_range
!
      real(dp), dimension(this%n_J, p_range%length, q_range%length), intent(out) :: L_Jpq
!
      integer, intent(in) :: block_
!
      integer :: r, s, p, q, k
      integer, dimension(2) :: r_and_s
!
      real(dp), dimension(:,:,:), pointer :: block_pointer
!
      r_and_s = this%get_range_indices_from_block(block_)
      r = r_and_s(1)
      s = r_and_s(2)
!
      call this%L_Jpq%get_array_n(block_pointer, block_)
!
!$omp parallel do private(p, q, k)
      do q = set_range_q%first, set_range_q%get_last()
         do p = set_range_p%first, set_range_p%get_last()
            do k = 1, this%n_J
!
               L_Jpq(k, p - p_range%first + 1, q - q_range%first + 1) = &
                     block_pointer(k, p - this%index_ranges(r)%first + 1, &
                                      q - this%index_ranges(s)%first + 1)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine copy_out_block
!
!
   subroutine cleanup(this)
!!
!!    Cleanup
!!    Written by Sarai D. Folkestad, Sep 2021
!!
      implicit none
!
      type(eri_cholesky_memory),  intent(inout)  :: this
!
      if (allocated(this%range_lengths)) call mem%dealloc(this%range_lengths, this%n_ranges)
!
      if (allocated(this%L_Jpq)) call this%L_Jpq%finalize()
!
   end subroutine cleanup
!
!
   function get_memory_estimate_eri_cholesky_memory(this, first_p, last_p, first_q, last_q) result(memory)
!!
!!    Memory estimate
!!    Written by Sarai D. Folkestad, Sep 2021
!!
!!    Returns the memory required by the Cholesky tool
!!    to handle the indicated Cholesky vector block_
!!
!
      use warning_suppressor, only: do_nothing
!
      implicit none
!
      class(eri_cholesky_memory), intent(in) :: this
!
      integer, intent(in) :: first_p, last_p
      integer, intent(in) :: first_q, last_q
!
      integer :: memory
!
      call do_nothing(this)
      call do_nothing([first_p, last_p, first_q, last_q])
!
      memory = 0
!
   end function get_memory_estimate_eri_cholesky_memory
!
!
end module eri_cholesky_memory_class
