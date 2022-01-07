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
module eri_cholesky_disk_class
!
!!
!! ERI Cholesky disk class
!! Written by Sarai D. Folkestad, Sep 2021
!!
!! Handles the disk storage and delivery of the
!! Cholesky vectors
!!
!
   use parameters
!
   use global_out,   only : output
!
   use direct_stream_file_class, only: direct_stream_file
   use range_class,              only: range_
   use memory_manager_class,     only: mem
!
   use abstract_eri_cholesky_class, only: abstract_eri_cholesky
!
   implicit none
!
   type, extends(abstract_eri_cholesky) :: eri_cholesky_disk
!
      type(direct_stream_file), dimension(:), allocatable, private :: file_
      character(len=40), private :: name_
!
   contains
!
      procedure, private :: write_to_block
      procedure, private :: read_from_block
!
      procedure, public :: initialize &
                        => initialize_eri_cholesky
!
      procedure, public :: get &
                        => get_eri_cholesky
!
      procedure, public :: set &
                        => set_eri_cholesky
!
      procedure, public :: get_memory_estimate &
                        => get_memory_estimate_eri_cholesky
!
      final :: cleanup
!
   end type eri_cholesky_disk
!
   interface eri_cholesky_disk
!
      procedure :: new_eri_cholesky
!
   end interface eri_cholesky_disk
!
!
contains
!
!
   pure function new_eri_cholesky(name_) result(this)
!!
!!    New
!!    Written by Sarai D. Foklkestad, Sep 2021
!!
      implicit none
!
      type(eri_cholesky_disk) :: this
      character(len=*), intent(in) :: name_
!
      this%n_J       = -1
      this%n_ranges  = -1
      this%dim_      = -1
      this%n_blocks  = -1
      this%name_     = name_
!
   end function new_eri_cholesky
!
!
   subroutine initialize_eri_cholesky(this, n_J, n_ranges, range_lengths)
!!
!!    Initialize
!!    Written by Sarai D. Foklkestad, Sep 2021
!!
      implicit none

      class(eri_cholesky_disk), intent(inout) :: this
!
      integer,                      intent(in) :: n_ranges, n_J
      integer, dimension(n_ranges), intent(in) :: range_lengths
!
      character (len=4) :: block_number
!
      integer :: i
!
      call this%set_dimensions(n_J, n_ranges, range_lengths)
      call this%set_index_ranges()
!
      allocate(this%file_(this%n_blocks))
!
      do i = 1, this%n_blocks
!
         write(block_number, '(i4.4)') i
         this%file_(i) = direct_stream_file('cholesky_'//trim(this%name_)//'_block_'//block_number,&
                                            this%n_J, w_size=dp)
      enddo
!
   end subroutine initialize_eri_cholesky
!
!
   subroutine set_eri_cholesky(this, L_Jpq, first_p, last_p, first_q, last_q)
!!
!!    Set
!!    Written by Sarai D. Foklkestad, Sep 2021
!!
      implicit none
!
      class(eri_cholesky_disk),                                      intent(inout):: this
      integer,                                                       intent(in)   :: first_p, last_p
      integer,                                                       intent(in)   :: first_q, last_q
      real(dp), dimension(this%n_J, first_p:last_p, first_q:last_q), intent(in)   :: L_Jpq
!
      type(range_) :: subrange_p, subrange_q, p_range, q_range
!
      integer :: block_, p, q, J, k, l
!
      real(dp), dimension(:,:,:), allocatable :: L_Jpq_copy
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
            call this%get_block_range_overlap(block_,           &
                                             p_range, q_range,  &
                                             subrange_p,       &
                                             subrange_q)
!
            if (subrange_p == p_range) then
!
               call this%write_to_block(block_, &
                                       subrange_p, &
                                       subrange_q, L_Jpq(:,:,subrange_q%first:subrange_q%get_last()))
            else
!
               call mem%alloc(L_Jpq_copy, this%n_J, subrange_p%length, subrange_q%length)
!
!$omp parallel do private(k, l, J)
               do k = 1, subrange_q%length
                  do l = 1, subrange_p%length
                     do J = 1, this%n_J
!
                        L_Jpq_copy(J, l, k) &
                           = L_Jpq(J, l + subrange_p%first - 1,&
                                      k + subrange_q%first - 1)
!
                     enddo
                  enddo
               enddo
!$omp end parallel do
!
               call this%write_to_block(block_, &
                                       subrange_p, &
                                       subrange_q, L_Jpq_copy)
!
               call mem%dealloc(L_Jpq_copy, this%n_J, subrange_p%length, subrange_q%length)
!
            endif
!
         enddo
      enddo
!
   end subroutine set_eri_cholesky
!
!
   subroutine get_eri_cholesky(this, L_Jpq, first_p, last_p, first_q, last_q)
!!
!!    Get
!!    Written by Sarai D. Foklkestad, Sep 2021
!!
      implicit none
!
      class(eri_cholesky_disk),                                      intent(inout):: this
      integer,                                                       intent(in)   :: first_p, last_p
      integer,                                                       intent(in)   :: first_q, last_q
      real(dp), dimension(this%n_J, first_p:last_p, first_q:last_q), intent(out)  :: L_Jpq
!
      type(range_) :: subrange_p, subrange_q, p_range, q_range
!
      integer :: block_, p, q, k, l, J
!
      real(dp), dimension(:,:,:), allocatable :: L_Jpq_copy
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
            if (subrange_p == p_range) then
!
               call this%read_from_block(block_, &
                                 subrange_p, &
                                 subrange_q, L_Jpq(:,:,subrange_q%first:subrange_q%get_last()))
            else
!
               call mem%alloc(L_Jpq_copy, this%n_J, subrange_p%length, subrange_q%length)
               call this%read_from_block(block_,     &
                                         subrange_p, &
                                         subrange_q, &
                                         L_Jpq_copy)
!
!$omp parallel do private(k, l, J)
               do k = 1, subrange_q%length
                  do l = 1, subrange_p%length
                     do J = 1, this%n_J
!
                        L_Jpq(J, l + subrange_p%first - 1, k + subrange_q%first - 1) = L_Jpq_copy(J, l, k)
!
                     enddo
                  enddo
               enddo
!$omp end parallel do
!
               call mem%dealloc(L_Jpq_copy, this%n_J, subrange_p%length, subrange_q%length)
!
            endif
!
         enddo
      enddo
!
   end subroutine get_eri_cholesky
!
!
   subroutine write_to_block(this, block_, write_p, write_q, L_Jpq)
!!
!!    Write to block_
!!    Written by Sarai D. Folkestad, Sep 2021
!!
      implicit none
!
      class(eri_cholesky_disk),                                      intent(inout) :: this
      type(range_),                                                  intent(in) :: write_p, write_q
      real(dp), dimension(this%n_J, write_p%length, write_q%length), intent(in) :: L_Jpq
!
      integer :: block_
      integer :: r, s, first_record, last_record, q
      integer :: p_reduced_first, q_reduced_first, p_reduced_last, q_reduced_last
!
      integer, dimension(2) :: r_and_s
!
      call this%file_(block_)%open_('write')
!
      r_and_s = this%get_range_indices_from_block(block_)
      r = r_and_s(1)
      s = r_and_s(2)
!
      p_reduced_first = write_p%first - this%index_ranges(r)%first + 1
      p_reduced_last = write_p%get_last() - this%index_ranges(r)%first + 1
!
      if (this%index_ranges(r) == write_p) then
!
         q_reduced_first = write_q%first - this%index_ranges(s)%first + 1
         q_reduced_last = write_q%get_last() - this%index_ranges(s)%first + 1
!
         first_record = this%index_ranges(r)%length*(q_reduced_first - 1) + p_reduced_first
         last_record = this%index_ranges(r)%length*(q_reduced_last - 1) + p_reduced_last
!
         call this%file_(block_)%write_(L_Jpq, first_record, last_record)
!
      else
!
         do q = write_q%first, write_q%get_last()
!
            q_reduced_first = q - this%index_ranges(s)%first + 1
            q_reduced_last = q - this%index_ranges(s)%first + 1
!
            first_record = this%index_ranges(r)%length*(q_reduced_first - 1) + p_reduced_first
            last_record = this%index_ranges(r)%length*(q_reduced_last - 1) + p_reduced_last
!
            call this%file_(block_)%write_(L_Jpq(:,:,q), first_record, last_record)
!
         enddo
!
      endif
!
      call this%file_(block_)%close_()
!
   end subroutine write_to_block
!
!
   subroutine read_from_block(this, block_, read_p, read_q, L_Jpq)
!!
!!    Read from block_
!!    Written by Sarai D. Folkestad, Sep 2021
!!
      implicit none
!
      class(eri_cholesky_disk),                                    intent(inout)  :: this
      type(range_),                                                intent(in)     :: read_p, read_q
      real(dp), dimension(this%n_J, read_p%length, read_q%length), intent(out)    :: L_Jpq
      integer :: block_
!
      integer, dimension(2) :: r_and_s
!
      integer :: r, s, first_record, last_record, q
      integer :: p_reduced_first, q_reduced_first, p_reduced_last, q_reduced_last, q_reduced
!
      call this%file_(block_)%open_('read')
!
      r_and_s = this%get_range_indices_from_block(block_)
      r = r_and_s(1)
      s = r_and_s(2)
!
      p_reduced_first = read_p%first - this%index_ranges(r)%first + 1
      p_reduced_last = read_p%get_last() - this%index_ranges(r)%first + 1
!
      if (this%index_ranges(r) == read_p) then
!
         q_reduced_first = read_q%first - this%index_ranges(s)%first + 1
         q_reduced_last = read_q%get_last() - this%index_ranges(s)%first + 1
!
         first_record = this%index_ranges(r)%length*(q_reduced_first - 1) + p_reduced_first
         last_record = this%index_ranges(r)%length*(q_reduced_last - 1) + p_reduced_last
!
         call this%file_(block_)%read_(L_Jpq, first_record, last_record)
!
      else
!
         do q = read_q%first, read_q%get_last()
!
            q_reduced = q - this%index_ranges(s)%first + 1
!
            first_record = this%index_ranges(r)%length*(q_reduced - 1) + p_reduced_first
            last_record  = this%index_ranges(r)%length*(q_reduced - 1) + p_reduced_last
!
            call this%file_(block_)%read_(L_Jpq(:,:, q - read_q%first + 1), first_record, last_record)
!
         enddo
!
      endif
!
      call this%file_(block_)%close_()
!
   end subroutine read_from_block
!
!
   function get_memory_estimate_eri_cholesky(this, first_p, last_p, first_q, last_q) result(memory)
!!
!!    Memory estimate
!!    Written by Sarai D. Folkestad, Sep 2021
!!
!!    Returns the memory required by the Cholesky tool
!!    to handle the indicated Cholesky vector block_
!!
      implicit none
!
      class(eri_cholesky_disk), intent(in) :: this
!
      integer, intent(in) :: first_p, last_p
      integer, intent(in) :: first_q, last_q
!
      integer :: memory
!
      type(range_) :: p_range, q_range
      logical :: p_full
!
      p_range = range_(first_p, last_p - first_p + 1)
      q_range = range_(first_q, last_q - first_q + 1)
!
      p_full = .false.
      if (any(p_range == this%index_ranges)) p_full = .true.
!
      memory = 0
      if (.not. p_full) memory = p_range%length*q_range%length*this%n_J
!
   end function get_memory_estimate_eri_cholesky
!
!
   subroutine cleanup(this)
!!
!!    Cleanup
!!    Written by Sarai D. Folkestad, Sep 2021
!!
      implicit none
!
      type(eri_cholesky_disk),  intent(inout)  :: this
!
      if (allocated(this%range_lengths)) call mem%dealloc(this%range_lengths, this%n_ranges)
!
   end subroutine cleanup
!
end module eri_cholesky_disk_class
