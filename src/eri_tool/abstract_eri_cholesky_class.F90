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
module abstract_eri_cholesky_class
!
!!
!! Abstract ERI Cholesky class
!! Written by Sarai D. Folkestad, Sep 2021
!!
!! Defines the interface for the ERI Cholesky tools.
!!
!! These classes handles the storage and delivery of
!! the Cholesky vectors.
!!
!! The vectors may be stored with a given block structure:
!!
!! Upon calling the 'initialize' routine
!! a blocks structure may be requested by specifying a
!! set of indices.
!!
!! Cholesky vectors are of size (dim_**2 * n_J)
!!
!! Example 1:
!!
!!    call L%initialize(n_J, 2, [n_o, n_v])
!!
!! results in each index p and q of L_pq^J
!! being partitioned into 2 ranges of lengths n_o
!! and n_v. Resulting in 4 Cholesky vector blocks
!!
!! Example 2:
!!
!!    call L%initialize(n_J, 1, [n_ao])
!!
!! results in the indices w and x of L_wx^J
!! to not be partitioned.
!!
!!
!! The Cholesky vectors can be set using the 'set' routine
!! and retrieved in two ways:
!!
!! 1. 'get' - Requires allocation of an array for the vectors,
!!            this may lead to higher memory usage than necessary,
!!            if the vectors are stored in memory.
!!
!! 2. 'load_block'         - Loads the requested block of the Cholesky factor into memory
!!                           and returns a pointer to the block (or a contiguous part of it)
!!    'offload_block'      - Frees the memory used for the loaded block, if any.
!!
   use parameters
!
   use global_out,                only : output
   use range_class,               only: range_
   use memory_manager_class,      only: mem
   use observable_class,          only: observable
   use cholesky_block_list_class, only: cholesky_block_list
!
   implicit none
!
   type, extends(observable), abstract :: abstract_eri_cholesky
!
      integer :: dim_, n_J
!
      integer :: n_ranges, n_blocks
!
      type(range_), dimension(:), allocatable :: index_ranges
      integer,      dimension(:), allocatable :: range_lengths
!
      type(cholesky_block_list), allocatable :: L_Jpq_loaded
!
   contains
!
!     Routines that should only be used/overwritten in descendants
!
      procedure :: set_dimensions
      procedure :: set_index_ranges
      procedure :: pq_in_block
      procedure :: get_block_overlaps
      procedure :: is_contiguous_in_block
      procedure :: get_contiguous_block_index
      procedure :: load_block_to_linked_list
      procedure :: offload_block_from_linked_list
      procedure :: get_ranges_from_block
!
      procedure (get_block_abstract), deferred :: get_block
      procedure (set_block_abstract), deferred :: set_block
!
!     Public routines
!
      procedure, public :: set
      procedure, public :: get
      procedure, public :: basis_transformation
      procedure, public :: set_equal_to
!
      procedure (initialize_abstract),           public, deferred :: initialize
      procedure (load_block_abstract),           public, deferred :: load_block
      procedure (offload_block_abstract),        public, deferred :: offload_block
      procedure (get_memory_estimate_abstract),  public, deferred :: get_memory_estimate
      procedure (load_memory_estimate_abstract), public, deferred :: load_memory_estimate
!
      procedure, private :: overlaps_with_block
!
   end type abstract_eri_cholesky
!
   abstract interface
!
      subroutine initialize_abstract(this, n_J, n_ranges, range_lengths)
!
         use parameters
         import abstract_eri_cholesky
!
         implicit none

         class(abstract_eri_cholesky), intent(inout) :: this
!
         integer,                      intent(in) :: n_ranges, n_J
         integer, dimension(n_ranges), intent(in) :: range_lengths
!
      end subroutine initialize_abstract
!
      subroutine set_block_abstract(this, block_, overlap_p, overlap_q, &
                                    L_Jpq, range_p, range_q)
!
         use parameters
         import abstract_eri_cholesky
         import range_
!!
         implicit none
!
         class(abstract_eri_cholesky), intent(inout):: this
!
         type(range_), intent(in) :: overlap_p, overlap_q
         type(range_), intent(in) :: range_p, range_q
         integer,      intent(in) :: block_
!
         real(dp), dimension(this%n_J, range_p%length, range_q%length), intent(in)   :: L_Jpq
!
      end subroutine set_block_abstract
!
      subroutine get_block_abstract(this, block_, overlap_p, overlap_q, &
                                    L_Jpq, range_p, range_q)
!
         use parameters
         import abstract_eri_cholesky
         import range_
!
         implicit none
!
         class(abstract_eri_cholesky), intent(inout):: this
!
         type(range_), intent(in) :: overlap_p, overlap_q
         type(range_), intent(in) :: range_p, range_q

         integer,      intent(in) :: block_
!
         real(dp), dimension(this%n_J, range_p%length, range_q%length), intent(out)   :: L_Jpq
!
      end subroutine get_block_abstract
!
!
      subroutine load_block_abstract(this, L_Jpq, first_p, last_p, first_q, last_q)
!
         use parameters
         import abstract_eri_cholesky
!
         implicit none
!
         class(abstract_eri_cholesky), intent(inout) :: this
         integer,                      intent(in)    :: first_p, last_p
         integer,                      intent(in)    :: first_q, last_q
!
         real(dp), dimension(:,:,:), pointer,   intent(out)   :: L_Jpq
!
      end subroutine load_block_abstract
!
!
      subroutine offload_block_abstract(this, first_p, last_p, first_q, last_q)
!
         import abstract_eri_cholesky
!
         implicit none
!
         class(abstract_eri_cholesky), intent(inout) :: this
         integer,                      intent(in)    :: first_p, last_p
         integer,                      intent(in)    :: first_q, last_q
!
      end subroutine offload_block_abstract
!
      function get_memory_estimate_abstract(this, first_p, last_p, first_q, last_q) result(memory)
!
         use parameters
         import abstract_eri_cholesky
!
         implicit none
!
         class(abstract_eri_cholesky), intent(in) :: this
!
         integer, intent(in) :: first_p, last_p
         integer, intent(in) :: first_q, last_q
!
         integer :: memory
!
      end function get_memory_estimate_abstract
!
!
      function load_memory_estimate_abstract(this, first_p, last_p, first_q, last_q) result(memory)
!
         use parameters
         import abstract_eri_cholesky
!
         implicit none
!
         class(abstract_eri_cholesky), intent(in) :: this
!
         integer, intent(in) :: first_p, last_p
         integer, intent(in) :: first_q, last_q
!
         integer :: memory
!
      end function load_memory_estimate_abstract
!
   end interface
!
contains
!
   subroutine set_dimensions(this, n_J, n_ranges, range_lengths)
!!
!!    Set dimensions
!!    Written by Sarai D. Folkestad, Sep 2021
!!
      implicit none
!
      class(abstract_eri_cholesky), intent(inout) :: this
!
      integer, intent(in) :: n_ranges, n_J
      integer, dimension(n_ranges) :: range_lengths
!
      this%n_J      = n_J
      this%n_ranges = n_ranges
      this%n_blocks = n_ranges**2
!
      call mem%alloc(this%range_lengths, this%n_ranges)
      this%range_lengths = range_lengths
!
      this%dim_ = sum(this%range_lengths)
!
   end subroutine set_dimensions
!
   subroutine set_index_ranges(this)
!!
!!    Set index ranges
!!    Written by Sarai D. Folkestad, Sep 2021
!!
      implicit none
!
      class(abstract_eri_cholesky), intent(inout) :: this

      integer :: current_first, block_
!
      allocate(this%index_ranges(this%n_ranges))
!
      current_first = 1
!
      do block_ = 1, this%n_ranges
!
         this%index_ranges(block_) = range_(current_first, this%range_lengths(block_))
         current_first = current_first + this%range_lengths(block_)
!
      enddo
!
   end subroutine set_index_ranges
!
!
   subroutine get_ranges_from_block(this, block_, range_p, range_q)
!!
!!    Get ranges indices from block
!!    Written by Sarai D. Folkestad, Sep 2021
!!
!!    Returns the range indices for the given block_
!!
!!    E.g. for n_ranges = 2 -> n_blocks = 4
!!
!!    block_ = 3 -> p_range_index = 2, q_range_index = 1
!!
      implicit none
!
      class(abstract_eri_cholesky), intent(in) :: this
      integer,                      intent(in) :: block_
!
      type(range_), intent(out) :: range_p, range_q
!
      integer :: p_range_index, q_range_index
      integer :: pq_block
!
      pq_block = 0
!
      do p_range_index = 1, this%n_ranges
         do q_range_index = 1, this%n_ranges
!
            pq_block = pq_block + 1
            if (pq_block == block_) then
!
               range_p = this%index_ranges(p_range_index)
               range_q = this%index_ranges(q_range_index)
!
            endif
!
         enddo
      enddo
!
   end subroutine get_ranges_from_block
!
!
   function pq_in_block(this, p, q, block_) result(is_in_block)
!!
!!    pq in block
!!    Written by Sarai D. Folkestad, Sep 2021
!!
!!    Returns .true. if the indices p and q are
!!    in the block_
!!
      implicit none
!
      class(abstract_eri_cholesky), intent(in) :: this
      integer, intent(in) :: p, q, block_
!
      logical :: is_in_block
!
      type(range_) :: range_r, range_s
!
      call this%get_ranges_from_block(block_, range_r, range_s)
!
      is_in_block = range_r%contains_(p) .and. range_s%contains_(q)
!
   end function pq_in_block
!
!
   subroutine get_block_overlaps(this, block_, range_p, range_q, &
                                 overlap_p_range, overlap_q_range)
!!
!!    Get block range overlap
!!    Written by Sarai D. Folkestad, Sep 2021
!!
!!    Returns the overlap between the
!!    ranges in the given block_ and the provided
!!    ranges (range_p, range_q)
!!
!
      implicit none
!
      class(abstract_eri_cholesky), intent(in)  :: this
      integer,                      intent(in)  :: block_
      type(range_),                 intent(in)  :: range_p, range_q
      type(range_),                 intent(out) :: overlap_p_range, overlap_q_range
!
      type(range_) :: range_r, range_s
!
      call this%get_ranges_from_block(block_, range_r, range_s)
!
      overlap_p_range = range_r%get_overlap(range_p)
      overlap_q_range = range_s%get_overlap(range_q)
!
   end subroutine get_block_overlaps
!
!
   subroutine basis_transformation(this, T)
!!
!!    Basis transformation
!!    Written by Rolf H. Myhre, May 2020
!!
!!    Based on routine written by Sarai D. Folkestad
!!
!!    Updates the Cholesky vectors using the
!!    transformation matrix T
!!
!!       L_J -> L'_J = T L_J T^T,
!!
!!    In cases where we can't hold two full vectors in memory,
!!    we have to construct all the intermediate vectors in batches and write them to temp_file
!!    before reading them in and constructing the final vectors.
!!
!
      use batching_index_class, only: batching_index
      use direct_stream_file_class, only: direct_stream_file
      use reordering, only: sort_123_to_132
!
      implicit none
!
      class(abstract_eri_cholesky), intent(inout) :: this
!
      real(dp), dimension(this%dim_, this%dim_), intent(in) :: T
!
      real(dp), dimension(:), allocatable, target :: L_J_1,   L_J_2
      real(dp), dimension(:,:,:), pointer         :: L_J_1_p, L_J_2_p
!
      type(batching_index) :: batcher
!
      logical :: all_in_mem
!
      integer :: batch, r, req_1
!
      type(direct_stream_file) :: temp_file
!
      batcher = batching_index(this%dim_)
!
      req_1 = 2*this%n_J*this%dim_ + this%get_memory_estimate(1, 1, 1, this%dim_)
!
      call mem%batch_setup(batcher, 0, req_1, tag='Cholesky basis transformation')
      all_in_mem = (batcher%num_batches .eq. 1)
!
      call mem%alloc(L_J_1, this%n_J * batcher%max_length * this%dim_)
      call mem%alloc(L_J_2, this%n_J * batcher%max_length * this%dim_)
!
      if (.not. all_in_mem) then
         temp_file = direct_stream_file('temp_file', this%n_J, dp, 'new')
         call temp_file%open_('readwrite')
      endif
!
!     Construct intermediate vectors L''_J_rq = sum_s L_J_rs T_qs
!
      do batch = 1, batcher%num_batches
!
         call batcher%determine_limits(batch)
!
         L_J_1_p(1:this%n_J, 1:batcher%length, 1:this%dim_) => L_J_1
         L_J_2_p(1:this%n_J, 1:batcher%length, 1:this%dim_) => L_J_2
!
         call this%get(L_J_1, batcher%first, batcher%get_last(), 1, this%dim_)
!
         call dgemm('N', 'T',                                     &
                    this%n_J*batcher%length, this%dim_, this%dim_,&
                    one,                                          &
                    L_J_1_p, this%n_J*batcher%length,             &
                    T, this%dim_,                                 &
                    zero,                                         &
                    L_J_2_p, this%n_J*batcher%length)
!
         L_J_1_p(1:this%n_J, 1:this%dim_, 1:batcher%length) => L_J_1
!
         call sort_123_to_132(L_J_2_p, L_J_1_p, this%n_J, batcher%length, this%dim_)
!
         if (.not. all_in_mem) then
            call temp_file%write_(L_J_1, (batcher%first-1)*this%dim_+1, &
                                  batcher%get_last()*this%dim_)
         endif
!
      enddo
!
!     Construct final vectors L'_J_pq = sum_r T_pr L''_J_rq
!
      do batch = 1, batcher%num_batches
!
         call batcher%determine_limits(batch)
!
         L_J_1_p(1:this%n_J, 1:batcher%length, 1:this%dim_) => L_J_1
         L_J_2_p(1:this%n_J, 1:batcher%length, 1:this%dim_) => L_J_2
!
         if (.not. all_in_mem) then
            do r = 1, this%dim_
               call temp_file%read_(L_J_1_p(:, :, r), &
                                    (r-1)*this%dim_ + batcher%first, &
                                    (r-1)*this%dim_ + batcher%get_last())
            enddo
         endif
!
         call dgemm('N', 'T',                                       &
                    this%n_J*batcher%length, this%dim_, this%dim_,  &
                    one,                                            &
                    L_J_1_p, this%n_J*batcher%length,               &
                    T, this%dim_,                                   &
                    zero,                                           &
                    L_J_2_p, this%n_J*batcher%length)

         L_J_1_p(1:this%n_J, 1:this%dim_, 1:batcher%length) => L_J_1
!
         call sort_123_to_132(L_J_2_p, L_J_1_p, this%n_J, batcher%length, this%dim_)
!
         call this%set(L_J_1_p, 1, this%dim_, batcher%first, batcher%get_last())
!
      enddo

      if (.not. all_in_mem) then
         call temp_file%close_('delete')
      endif
!
      call mem%dealloc(L_J_1, this%n_J * batcher%max_length * this%dim_)
      call mem%dealloc(L_J_2, this%n_J * batcher%max_length * this%dim_)
!
      call mem%batch_finalize()
      call this%notify_observers()
!
   end subroutine basis_transformation
!
   subroutine set_equal_to(this, that)
!!
!!    Set equal to
!!    Written by Sarai D. Folkestad, Sep 2021
!!
!
      use batching_index_class, only: batching_index
!
      implicit none
!
      class(abstract_eri_cholesky), intent(inout) :: this
      class(abstract_eri_cholesky), intent(inout) :: that

      type(batching_index), allocatable :: batch_q
!
      integer :: req_0, req_1, batch
!
      real(dp), dimension(:,:,:), allocatable :: L_J_pq
!
      batch_q = batching_index(this%dim_)
!
      req_0 = 0
      req_1 = this%n_J*this%dim_ &
            + that%get_memory_estimate(1, this%dim_, 1, 1) &
            + this%get_memory_estimate(1, this%dim_, 1, 1)
!
      call mem%batch_setup(batch_q, req_0, req_1, tag='Copy Cholesky vector')
!
      call mem%alloc(L_J_pq, this%n_J, this%dim_, batch_q%max_length)
!
      do batch = 1, batch_q%num_batches
!
         call batch_q%determine_limits(batch)
!
         call that%get(L_J_pq, 1, this%dim_, batch_q%first, batch_q%get_last())
         call this%set(L_J_pq, 1, this%dim_, batch_q%first, batch_q%get_last())
!
      enddo
!
      call mem%dealloc(L_J_pq, this%n_J, this%dim_, batch_q%max_length)
!
      call mem%batch_finalize()
!
      call this%notify_observers()
!
   end subroutine set_equal_to
!
!
   function is_contiguous_in_block(this, range_p, range_q) result(contiguous_)
!!
!!    Is contiguous in block
!!    Written by Sarai Dery Folkestad, Jan 2022
!!
!!    Returns true if the given ranges correspond to contiguous chunk
!!    in a block
!!
      implicit none
!
      class(abstract_eri_cholesky), intent(in) :: this
      type(range_), intent(in) :: range_p, range_q
!
      logical :: contiguous_
!
      integer :: p, q
!
      contiguous_ = .false.
!
      do p = 1, this%n_ranges
!
         if (.not. range_p == this%index_ranges(p)) cycle
!
         do q = 1, this%n_ranges
!
            if (this%index_ranges(q)%contains_(range_q)) then
               contiguous_ = .true.
               return
            endif
!
         enddo
      enddo
!
   end function is_contiguous_in_block
!
!
   function get_contiguous_block_index(this, range_p, range_q) result(index)
!!
!!    Get contiguous block index
!!    Written by Sarai Dery Folkestad, Jan 2022
!!
!!    Returns the block where the given ranges correspond to a contiguous chunk,
!!    returns -1 if no such block exists
!!
      implicit none
!
      class(abstract_eri_cholesky), intent(inout):: this
      type(range_), intent(in) :: range_p, range_q
!
      integer :: index
      integer ::  block_, p, q
!
      index = -1
      block_ = 0
!
      do p = 1, this%n_ranges
!
         do q = 1, this%n_ranges
!
            block_ = block_ + 1
!
            if (.not. range_p == this%index_ranges(p)) cycle
!
            if (this%index_ranges(q)%contains_(range_q)) then
               index = block_
               return
            endif
!
         enddo
      enddo
!
   end function get_contiguous_block_index
!
!
   subroutine load_block_to_linked_list(this, range_p, range_q)
!!
!!    Load block to linked list
!!    Written by Sarai Dery Folkestad, Jan 2022
!!
!!    Loads a block (given by the ranges) into the
!!    L_pq_loaded list
!!
      implicit none
!
      class(abstract_eri_cholesky), intent(inout):: this
      type(range_),                 intent(in)   :: range_p, range_q
!
      real(dp), dimension(:,:,:), pointer :: L_ptr
!
      call this%L_Jpq_loaded%push_back(range_p, range_q)
      call this%L_Jpq_loaded%get_array(L_ptr, range_p, range_q)
!
      call this%get(L_ptr, range_p%first, range_p%get_last(), range_q%first, range_q%get_last())
!
   end subroutine load_block_to_linked_list
!
!
   subroutine offload_block_from_linked_list(this, range_p, range_q)
!!
!!    Offload block from linked list
!!    Written by Sarai Dery Folkestad, Jan 2022
!!
!!    Offloads a block (given by the ranges) into the
!!    L_pq_loaded list
!!
      implicit none
!
      class(abstract_eri_cholesky), intent(inout):: this
      type(range_),                 intent(in)   :: range_p, range_q
!
      call this%L_Jpq_loaded%remove(range_p, range_q)
!
   end subroutine offload_block_from_linked_list
!
!
   subroutine set(this, L_Jpq, first_p, last_p, first_q, last_q)
!!
!!    Set
!!    Written by Sarai D. Folkestad, Sep 2021
!!
      implicit none
!
      class(abstract_eri_cholesky), intent(inout) :: this
      integer,                      intent(in)    :: first_p, last_p
      integer,                      intent(in)    :: first_q, last_q
!
      real(dp), dimension(this%n_J, first_p:last_p, first_q:last_q), intent(in)   :: L_Jpq
!
      type(range_) :: overlap_p, overlap_q, range_p, range_q
!
      integer :: block_
!
      range_p = range_(first_p, last_p - first_p + 1)
      range_q = range_(first_q, last_q - first_q + 1)
!
      do block_ = 1, this%n_blocks
!
         if (.not. this%overlaps_with_block(range_p, range_q, block_)) cycle
         call this%get_block_overlaps(block_, range_p, range_q, overlap_p, overlap_q)
         call this%set_block(block_,  overlap_p, overlap_q, L_Jpq, range_p, range_q)
!
      enddo
!
   end subroutine set
!
!
   subroutine get(this, L_Jpq, first_p, last_p, first_q, last_q)
!!
!!    Get
!!    Written by Sarai D. Folkestad, Sep 2021
!!
      implicit none
!
      class(abstract_eri_cholesky), intent(inout):: this
      integer,                      intent(in)   :: first_p, last_p
      integer,                      intent(in)   :: first_q, last_q
!
      real(dp), dimension(this%n_J, first_p:last_p, first_q:last_q), intent(out)   :: L_Jpq
!
      type(range_) :: overlap_p, overlap_q, range_p, range_q
!
      integer :: block_
!
      range_p = range_(first_p, last_p - first_p + 1)
      range_q = range_(first_q, last_q - first_q + 1)
!
      do block_ = 1, this%n_blocks
!
         if (.not. this%overlaps_with_block(range_p, range_q, block_)) cycle
         call this%get_block_overlaps(block_, range_p, range_q, overlap_p, overlap_q)
         call this%get_block(block_, overlap_p, overlap_q, L_Jpq, range_p, range_q)
!
      enddo
!
   end subroutine get
!
!
   function overlaps_with_block(this, range_p, range_q, block_) result(overlaps)
!!
!!    Overlaps with block
!!    Written by Sarai D. Folkestad, Sep 2021
!!
      implicit none
!
      class(abstract_eri_cholesky), intent(in):: this
!
      type(range_), intent(in) :: range_p, range_q
!
      integer, intent(in) :: block_
!
      logical :: overlaps
!
      type(range_) :: range_r, range_s
!
      call this%get_ranges_from_block(block_, range_r, range_s)
      overlaps = (range_p%overlaps(range_r)) .and. (range_q%overlaps(range_s))
!
   end function overlaps_with_block
!
!
end module abstract_eri_cholesky_class
