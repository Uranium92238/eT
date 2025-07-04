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
module batching_index_class
!
!!
!!    Batching index class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Dec 2017
!!
!!    The batching index class represents a single index (e.g. a in g_abcd)
!!    that is being batched over. During a batching process, it contains information
!!    relevant to the restriction of that index, such as the first and last indices,
!!    the length of the batching region as well as the total number of batches.
!!
!!    A typical use of the batching index is as follows:
!!
!!       type(batching_index) :: batch_a
!!
!!       batch_a = batching_index(wf%n_v) -> initializes batching object for an index of dimension
!!                                            equal to the number of virtual orbitals
!!
!!       call mem%batch_setup(batch_a, req0, req1) -> determines the number of batches and
!!                                                   and the largest batching length & saves
!!                                                   that information in the batch_a object
!!
!!       do current_a_batch = 1, batch_a%num_batches
!!
!!          call batch_a%determine_limits(current_a_batch) -> determines the first and last index values,
!!                                                            and the length, of the batching range
!!                                                            for the specific batch 'current_a_batch',
!!                                                            saving that info. in batch_a
!!
!!          Do some stuff. Useful & available info:
!!
!!             batch_a%first  -> first value of index a in the current batch
!!             batch_a%last   -> last value of index a in the current batch
!!             batch_a%length -> length of current batch, last - first + 1
!!
!!       enddo
!!
!!    The case is similar for a batching over several indices simultaneously. See
!!    the memory manager "batch_setup" routines for details.
!!
!!
!
   use kinds
   use parameters
   use global_out, only : output
   use range_class
!
   type, extends(range_) :: batching_index
!
!     Values determined by the memory manager
!
      integer :: max_length      ! Maximum length of a batch
      integer :: num_batches     ! The number of batches in total for the index
!
!     Value that must be initialized by user
!
      integer :: index_dimension ! Full length of index (e.g., typically n_vir for virtual index)
!
!     Offset to add when calculating batch limits
!
      integer :: offset
!
!     Logical for sanity check
!
      logical :: initialized = .false.
!
   contains
!
!     Routine that sets the batch dependent variables,
!     first, last and length, based on which batch it is
!
      procedure :: determine_limits    => determine_limits_batching_index
!
!     Debug option:
!     Forced batching routine called by memory manager to ensure
!     batching regardless of available memory. Batch size is randomly generated.
!
      procedure :: force_batch         => force_batch_batching_index
!
!     Make sure that the batching index does not batch (num batches: 0)
!
      procedure :: do_not_batch        => do_not_batch_batching_index
!
!     Make sure that the batching index does a single batch (num batches: 1)
!
      procedure :: do_single_batch     => do_single_batch_batching_index
!
   end type batching_index
!
!
   interface batching_index
!
      procedure :: new_batching_index
!
   end interface batching_index
!
!
contains
!
!
   function new_batching_index(dimension_, offset) result(batch_p)
!!
!!    New batching index
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Dec 2017
!!    Modified by Eirik F. Kjønstad, Mar 2020
!!
!!    Note: every batching index must be initialized!
!!
!!    dimension_: the total length of the index, e.g. the number of virtual
!!                orbitals "n_v" for a virtual index "a".
!!
!!    offset:     (optional) offset for the index; e.g., if we have an
!!                MO index p that is restricted to virtual orbitals, we can
!!                set dimension_ to "n_v" and the offset to "n_o" (so that the
!!                index can take the values n_o + 1, n_o + 2, ..., n_o + n_v).
!!                Default is offset=0.
!!
!!    Eirik F. Kjønstad, Mar 2020: added offset.
!!
      implicit none
!
      type(batching_index) :: batch_p
!
      integer, intent(in)           :: dimension_
      integer, intent(in), optional :: offset
!
      batch_p%index_dimension = dimension_
!
      batch_p%initialized = .true.
!
      batch_p%max_length  = 0
      batch_p%num_batches = 0
!
      batch_p%offset = 0
      if (present(offset)) batch_p%offset = offset
!
   end function new_batching_index
!
!
   subroutine determine_limits_batching_index(this, batch_number)
!!
!!    Determine limits
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Dec 2017
!!    Modified by Eirik F. Kjønstad, Mar 2020
!!
!!    Given the batch number, this routine determine the first and
!!    and last values, for the index, as well as the length of the
!!    current batching interval.
!!
!!    Eirik F. Kjønstad, Mar 2020: added offset.
!!
      implicit none
!
      class(batching_index) :: this
!
      integer, intent(in) :: batch_number ! The current batch
!
!     Sanity check
!
      if (.not. this%initialized) then
!
         call output%error_msg('a non-initialized batching variable was used.')
!
      endif
!
!     Determine limits of batch, q = first, first + 1, ..., last
!
      if (batch_number .lt. this%num_batches) then
         this%range_ = range_(this%offset + 1 + (batch_number-1)*this%max_length, &
                              this%max_length)
      else
         this%range_ = range_(this%offset + 1 + (batch_number-1)*this%max_length, &
                              this%index_dimension - (batch_number-1)*this%max_length)
      endif
!
   end subroutine determine_limits_batching_index
!
!
   subroutine force_batch_batching_index(this)
!!
!!    Force batch
!!    Written by Eirik F. Kjønstad, Sep 2019
!!
      implicit none
!
      class(batching_index) :: this
!
      real(dp) :: some_number_between_0_and_1 ! [0, 1)
!
      call random_number(some_number_between_0_and_1)
!
!     1, 2, 3, ..., index_dimension - 1
      this%max_length = 1 + floor((this%index_dimension - 1)*some_number_between_0_and_1)
!
      this%num_batches = (this%index_dimension-1)/(this%max_length)+1
!
      call output%printf('debug', 'Forced batch of index with dimension: (i0). &
                         &Number of batches: (i0), Max length of batch: (i0)', &
                         ints=[this%index_dimension, this%num_batches,   &
                         this%max_length], ll=90)
!
   end subroutine force_batch_batching_index
!
!
   subroutine do_not_batch_batching_index(this)
!!
!!    Do not batch
!!    Written by Eirik F. Kjønstad, Apr 2020
!!
!!    Does initialization needed for performing no batching,
!!    i.e. setting the number of batches to 0.
!!
      implicit none
!
      class(batching_index), intent(inout) :: this
!
      this%max_length  = 0
      this%num_batches = 0
!
   end subroutine do_not_batch_batching_index
!
!
   subroutine do_single_batch_batching_index(this)
!!
!!    Do single batch
!!    Written by Eirik F. Kjønstad, Apr 2020
!!
!!    Does initialization needed to do one batch with entire range,
!!    i.e. setting the number of batches to 1.
!!
      implicit none
!
      class(batching_index), intent(inout) :: this
!
      this%max_length  = this%index_dimension
      this%num_batches = 1
!
      call this%determine_limits(1)
!
   end subroutine do_single_batch_batching_index
!
!
end module batching_index_class
