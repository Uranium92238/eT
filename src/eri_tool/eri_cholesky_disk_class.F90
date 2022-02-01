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
   use global_out,                  only : output
   use direct_stream_file_class,    only: direct_stream_file
   use range_class,                 only: range_
   use memory_manager_class,        only: mem
   use abstract_eri_cholesky_class, only: abstract_eri_cholesky
!
   implicit none
!
   type, extends(abstract_eri_cholesky) :: eri_cholesky_disk
!
      type(direct_stream_file), dimension(:), allocatable, private :: file_
!
      character(len=40), private :: name_
!
   contains
!
      procedure, private :: write_to_block
      procedure, private :: read_from_block
!
      procedure, public :: initialize &
                        => initialize_eri_cholesky_disk
!
      procedure, public :: get_block &
                        => get_block_eri_cholesky_disk
!
      procedure, public :: set_block &
                        => set_block_eri_cholesky_disk
!
      procedure, public :: get_memory_estimate &
                        => get_memory_estimate_eri_cholesky_disk
!
      procedure, public :: load_memory_estimate &
                        => load_memory_estimate_eri_cholesky_disk
!
      procedure, public :: load_block &
                        => load_block_eri_cholesky_disk
!
      procedure, public :: offload_block &
                        => offload_block_eri_cholesky_disk
!
      final :: cleanup
!
   end type eri_cholesky_disk
!
   interface eri_cholesky_disk
!
      procedure :: new_eri_cholesky_disk
!
   end interface eri_cholesky_disk
!
!
contains
!
!
   pure function new_eri_cholesky_disk(name_) result(this)
!!
!!    New
!!    Written by Sarai D. Foklkestad, Sep 2021
!!
!
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
   end function new_eri_cholesky_disk
!
!
   subroutine initialize_eri_cholesky_disk(this, n_J, n_ranges, range_lengths)
!!
!!    Initialize
!!    Written by Sarai D. Foklkestad, Sep 2021
!!
      use cholesky_block_list_class,      only: cholesky_block_list
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
      this%L_Jpq_loaded = cholesky_block_list(this%n_J)
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
   end subroutine initialize_eri_cholesky_disk
!
!
   subroutine set_block_eri_cholesky_disk(this, block_, overlap_p, overlap_q, &
                                    L_Jpq, range_p, range_q)
!!
!!    Set block
!!    Written by Sarai D. Folkestad, Sep 2021
!!
!!    Sets a subblock of the given block, given by overlap_p and overlap_q
!!
      implicit none
!
      class(eri_cholesky_disk), intent(inout):: this
!
      type(range_),  intent(in) :: overlap_p, overlap_q
         type(range_), intent(in) :: range_p, range_q
      integer,       intent(in) :: block_
!
      real(dp), dimension(this%n_J, range_p%length, range_q%length), intent(in)   :: L_Jpq
!
      real(dp), dimension(:,:,:), allocatable :: L_Jpq_copy
!
      integer :: p, q, J
!
      if (overlap_p == range_p) then
!
         call this%write_to_block(block_, &
                                 overlap_p, &
                                 overlap_q, L_Jpq(:,:,1:overlap_q%length))
      else
!
         call mem%alloc(L_Jpq_copy, this%n_J, overlap_p%length, overlap_q%length)
!
!$omp parallel do private(p, q, J)
      do q = overlap_q%first, overlap_q%get_last()
         do p = overlap_p%first, overlap_p%get_last()
               do J = 1, this%n_J
!
                  L_Jpq_copy(J, p - overlap_p%first + 1, q - overlap_q%first + 1) &
                     = L_Jpq(J, p - range_p%first + 1, q - range_q%first + 1)
!
               enddo
            enddo
         enddo
!$omp end parallel do
!
         call this%write_to_block(block_, &
                                 overlap_p, &
                                 overlap_q, L_Jpq_copy)
!
         call mem%dealloc(L_Jpq_copy, this%n_J, overlap_p%length, overlap_q%length)
!
      endif
!
   end subroutine set_block_eri_cholesky_disk
!
!
   subroutine get_block_eri_cholesky_disk(this, block_, overlap_p, overlap_q, &
                                    L_Jpq, range_p, range_q)
!!
!!    Get block
!!    Written by Sarai D. Folkestad, Sep 2021
!!
!!    Gets a subblock of the given block, given by overlap_p and overlap_q
!!
      implicit none
!
      class(eri_cholesky_disk), intent(inout):: this
!
      type(range_), intent(in) :: overlap_p, overlap_q
      type(range_), intent(in) :: range_p, range_q
      integer,      intent(in) :: block_
!
      real(dp), dimension(this%n_J, range_p%length, range_q%length), intent(out)   :: L_Jpq
!
      real(dp), dimension(:,:,:), allocatable :: L_Jpq_copy
!
      integer ::  p, q, J
!
      type(range_) :: range_r, range_s
!
      call this%get_ranges_from_block(block_, range_r, range_s)
!
      if (overlap_p == range_p) then
!
         call this%read_from_block(block_, &
                           overlap_p, &
                           overlap_q, L_Jpq(:,:,1:overlap_q%length))
      else
!
         call mem%alloc(L_Jpq_copy, this%n_J, overlap_p%length, overlap_q%length)
         call this%read_from_block(block_,    &
                                   overlap_p, &
                                   overlap_q, &
                                   L_Jpq_copy)
!
!$omp parallel do private(q, p, J)
         do q = overlap_q%first, overlap_q%get_last()
            do p = overlap_p%first, overlap_p%get_last()
               do J = 1, this%n_J
!
                  L_Jpq(J, p - range_p%first + 1, q - range_q%first + 1) &
                        = L_Jpq_copy(J, p - overlap_p%first + 1, &
                                        q - overlap_q%first + 1)
!
               enddo
            enddo
         enddo
!$omp end parallel do
!
         call mem%dealloc(L_Jpq_copy, this%n_J, overlap_p%length, overlap_q%length)
!
      endif
!
   end subroutine get_block_eri_cholesky_disk
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
      integer :: first_record, last_record, q
      integer :: p_reduced_first, q_reduced_first, p_reduced_last, q_reduced_last
!
      type(range_) :: range_r, range_s
!
      call this%file_(block_)%open_('write')
!
      call this%get_ranges_from_block(block_, range_r, range_s)
!
      p_reduced_first = write_p%first - range_r%first + 1
      p_reduced_last  = write_p%get_last() - range_r%first + 1
!
      if (range_r == write_p) then
!
         q_reduced_first = write_q%first - range_s%first + 1
         q_reduced_last  = write_q%get_last() - range_s%first + 1
!
         first_record = range_r%length*(q_reduced_first - 1) + p_reduced_first
         last_record  = range_r%length*(q_reduced_last  - 1) + p_reduced_last
!
         call this%file_(block_)%write_(L_Jpq, first_record, last_record)
!
      else
!
         do q = write_q%first, write_q%get_last()
!
            q_reduced_first = q - range_s%first + 1
            q_reduced_last  = q - range_s%first + 1
!
            first_record = range_r%length*(q_reduced_first - 1) + p_reduced_first
            last_record  = range_r%length*(q_reduced_last  - 1) + p_reduced_last
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
      type(range_) :: range_r, range_s
!
      integer :: first_record, last_record, q
      integer :: p_reduced_first, q_reduced_first, p_reduced_last, q_reduced_last, q_reduced
!
      call this%file_(block_)%open_('read')
!
      call this%get_ranges_from_block(block_, range_r, range_s)
!
      p_reduced_first = read_p%first - range_r%first + 1
      p_reduced_last = read_p%get_last() - range_r%first + 1
!
      if (range_r == read_p) then
!
         q_reduced_first = read_q%first - range_s%first + 1
         q_reduced_last = read_q%get_last() - range_s%first + 1
!
         first_record = range_r%length*(q_reduced_first - 1) + p_reduced_first
         last_record = range_r%length*(q_reduced_last - 1) + p_reduced_last
!
         call this%file_(block_)%read_(L_Jpq, first_record, last_record)
!
      else
!
         do q = read_q%first, read_q%get_last()
!
            q_reduced = q - range_s%first + 1
!
            first_record = range_r%length*(q_reduced - 1) + p_reduced_first
            last_record  = range_r%length*(q_reduced - 1) + p_reduced_last
!
            call this%file_(block_)%read_(L_Jpq(:,:, q - read_q%first + 1), &
                                          first_record, last_record)
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
   function get_memory_estimate_eri_cholesky_disk(this, first_p, last_p, first_q, last_q) result(memory)
!!
!!    Get memory estimate
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
      type(range_) :: range_p, range_q
!
      memory = 0
!
      if (any(range_p == this%index_ranges)) return
!
      range_p = range_(first_p, last_p - first_p + 1)
      range_q = range_(first_q, last_q - first_q + 1)
!
      memory = range_p%length*range_q%length*this%n_J
!
   end function get_memory_estimate_eri_cholesky_disk
!
!
   function load_memory_estimate_eri_cholesky_disk(this, first_p, last_p, first_q, last_q) result(memory)
!!
!!    Load memory estimate
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
      type(range_) :: range_p, range_q
!
      range_p = range_(first_p, last_p - first_p + 1)
      range_q = range_(first_q, last_q - first_q + 1)
!
      memory = range_p%length*range_q%length*this%n_J
!
   end function load_memory_estimate_eri_cholesky_disk
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
!
   subroutine load_block_eri_cholesky_disk(this,  L_Jpq, first_p, last_p, first_q, last_q)
!!
!!    Load block
!!    Written by Sarai D. Folkestad, Jan 2022
!!
      implicit none
!
      class(eri_cholesky_disk), intent(inout):: this
      integer,                     intent(in)   :: first_p, last_p
      integer,                     intent(in)   :: first_q, last_q
!
      real(dp), dimension(:,:,:), pointer,   intent(out)   :: L_Jpq
!
      type(range_) :: range_p, range_q
!
      range_p = range_(first_p, last_p - first_p + 1)
      range_q = range_(first_q, last_q - first_q + 1)
!
      call this%load_block_to_linked_list(range_p, range_q)
!
      call this%L_Jpq_loaded%get_array(L_Jpq, range_p, range_q)
!
   end subroutine load_block_eri_cholesky_disk
!
!
   subroutine offload_block_eri_cholesky_disk(this, first_p, last_p, first_q, last_q)
!!
!!    Offload block
!!    Written by Sarai D. Folkestad, Jan 2022
!!
      implicit none
!
      class(eri_cholesky_disk), intent(inout):: this
      integer,                      intent(in)    :: first_p, last_p
      integer,                      intent(in)    :: first_q, last_q
!
      type(range_)      :: range_p, range_q
!
      range_p = range_(first_p, last_p - first_p + 1)
      range_q = range_(first_q, last_q - first_q + 1)
!
      if (this%L_Jpq_loaded%node_in_list(range_p, range_q)) then
!
         call this%offload_block_from_linked_list(range_p, range_q)
!
      else
!
         call output%error_msg('cannot find Cholesky block!')
!
      endif
!
   end subroutine offload_block_eri_cholesky_disk
!
!
end module eri_cholesky_disk_class
