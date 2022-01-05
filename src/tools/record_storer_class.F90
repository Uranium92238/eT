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
module record_storer_class
!
!!
!!    Record storer class module
!!    Written by Eirik F. Kjønstad, 2020
!!
!!    Object for storing a set of records (real, double precision). Provides
!!    a single interface for records that are either stored on file or in memory,
!!    depending on user settings or the memory available to the calculation.
!!
!
   use kinds
   use range_class
   use memory_manager_class, only: mem
   use global_out, only: output
   use batching_index_class, only: batching_index
   use array_list_class, only: array_list
   use direct_stream_file_class, only: direct_stream_file
!
   implicit none
!
   type :: record_storer
!
!     Name of storer - determines file name for stored records
!
      character(len=255), private :: name_ = 'no_name'
!
!     Maximum number of records & dimensionality of each record
!
      integer, private :: n_records
      integer, private :: record_dim
!
!     Linked list of arrays for storing records in memory
!
      type(array_list), allocatable, private :: linked_arrays
!
!     Pointer for array (if in_memory = true)
!
      real(dp), dimension(:,:), pointer, contiguous, private :: array_p
!
!     Are the records stored in memory?
!
      logical, private :: in_memory
!
!     File for storing records
!
      type(direct_stream_file), allocatable, private :: file_
!
!     Storing ranges for each of the blocks loaded into memory
!
      integer, private :: n_blocks
      type(range_), dimension(:), allocatable, private :: blocks
!
   contains
!
!     Three main routines to access records in storer
!
!     - Prepare to hold a set of records in memory
!     - Load sets a pointer to array where the records are stored
!     - Free up space allocated with prepare
!
      procedure, public :: prepare_records &
                        => prepare_records_record_storer
!
      generic, public   :: load &
                        => load_records, &
                           load_range
!
      procedure, public :: free_records &
                        => free_records_record_storer
!
!     Initializations and finalizations
!
      procedure, public :: initialize &
                        => initialize_record_storer
!
      procedure, public :: finalize &
                        => finalize_record_storer
!
!     # elements that will be allocated on load of one record
!
      procedure, public :: required_to_load_record &
                        => required_to_load_record_record_storer
!
!     Set/get records by copying out/into storer
!
      generic, public :: copy_record_in &
                      => copy_record_in_single_record, &
                         copy_record_in_range
!
      generic, public :: copy_record_out &
                      => copy_record_out_single_record, &
                         copy_record_out_range
!
!     Store changes made to records temporarily held in memory
!
      procedure, public :: store &
                        => store_record_storer
!
      procedure, private :: update_loaded_blocks ! If file storage, some records may be loaded into memory
                                                 ! If so, this routine updates these loaded blocks when the
                                                 ! records are changed
!
      procedure, private :: load_records
      procedure, private :: load_range
!
      procedure, private :: copy_record_out_single_record
      procedure, private :: copy_record_out_range
!
      procedure, private :: copy_record_in_single_record
      procedure, private :: copy_record_in_range
!
   end type record_storer
!
!
   interface record_storer
!
      procedure :: new_record_storer
!
   end interface record_storer
!
!
contains
!
!
   function new_record_storer(name_,            &
                              record_dim,       &
                              n_records,        &
                              records_in_mem) result(storer)
!!
!!    Record storer constructor
!!    Writen by Eirik F. Kjønstad, 2019
!!
!!    name_:            Name of storer.
!!
!!    record_dim:       The length of each record (i.e. the dimensionality of the array
!!                      you want to store in each record)
!!
!!    n_records:        The number of records (you can store arrays in
!!                      records 1,2,3,..., n_records)
!!
!!    records_in_mem:   Store records in memory if true.
!!                      Stores records on direct stream file if false.
!!
      implicit none
!
      type(record_storer) :: storer
!
      character(len=*), intent(in) :: name_
!
      integer, intent(in) :: record_dim
      integer, intent(in) :: n_records
      logical, intent(in) :: records_in_mem
!
      storer%name_       = trim(name_)
      storer%record_dim  = record_dim
      storer%n_records   = n_records
!
      storer%in_memory = records_in_mem
      storer%n_blocks = 0
!
      if (storer%in_memory) then
!
         call output%printf('n', 'Storage ('// trim(storer%name_)//'): memory', fs='(t6,a)')
!
      else
!
         call output%printf('n', 'Storage ('//trim(storer%name_)//'): file', fs='(t6,a)')
!
      endif
!
   end function new_record_storer
!
!
   function required_to_load_record_record_storer(storer) result(required)
!!
!!    Required to load record
!!    Written by Eirik F. Kjønstad, Mar 2020
!!
!!    Returns the number of elements that will be allocated
!!    if a single record is loaded.
!!
!!    Equals the record dimension if records are read from file into memory.
!!    Equals zero if records are already in memory.
!!
      implicit none
!
      class(record_storer) :: storer
!
      integer :: required
!
      if (storer%in_memory) then
!
         required = 0
!
      else
!
         required = storer%record_dim
!
      endif
!
   end function required_to_load_record_record_storer
!
!
   subroutine prepare_records_record_storer(storer, batches)
!!
!!    Prepare records
!!    Written by Eirik F. Kjønstad, 2020
!!
!!    Makes preparations needed to load records for one or more
!!    batching indices. No preparations if in memory. Allocates arrays
!!    for each batching index if records are on file.
!!
!!    batches: Array of batching indices. The number of batching indices
!!             you need is equal to the number of records you need to hold
!!             simultaneously. E.g., if you have a loop over i and j and
!!             possibly need to hold different sets of records for i and j.
!!
!!             This is to ensure that a file storer can temporarily load
!!             i and j records independently into memory.
!!
!!    The routine allocates enough memory to hold any batch, i.e.
!!    it allocates the maximum batch length for the index. The allocated
!!    number of elements is therefore
!!
!!       # elements allocated = record_dim * ( batches(1)%max_length
!!                                           + batches(2)%max_length
!!                                           + ... )
!!
      implicit none
!
      class(record_storer) :: storer
!
      type(batching_index), dimension(:) :: batches
!
      integer :: k
!
      if (.not. storer%in_memory) then
!
!        Open file with read-write access
!
         call storer%file_%open_('readwrite')
!
!        Set up an array for each batching index
!        by adding arrays to the linked list
!
!        Records will be read into these arrays
!
         storer%n_blocks = size(batches)
!
         do k = 1, storer%n_blocks
!
           call storer%linked_arrays%push_back(storer%record_dim, &
                                               batches(k)%max_length)
!
         enddo
!
!        Allocate blocks for holding range information
!
         allocate(storer%blocks(storer%n_blocks))
!
      endif
!
   end subroutine prepare_records_record_storer
!
!
   subroutine free_records_record_storer(storer)
!!
!!    Free records
!!    Written by Eirik F. Kjønstad, 2020
!!
!!    Frees memory temporarily used to keep records.
!!
!!    No action if in memory. Deallocates arrays allocated in
!!    the prepare_records routine if on file.
!!
      implicit none
!
      class(record_storer) :: storer
!
      if (.not. storer%in_memory) then
!
!        Deallocate all arrays temporarily allocated in the linked list
!
         call storer%linked_arrays%finalize()
!
!        Deallocate blocks for holding range information
!
         deallocate(storer%blocks)
!
         storer%n_blocks = 0
!
!        Close file
!
         call storer%file_%close_()
!
      endif
!
   end subroutine free_records_record_storer
!
!
   subroutine copy_record_out_single_record(storer, x, n)
!!
!!    Copy record out
!!    Written by Eirik F. Kjønstad, 2020
!!
!!    Copies the nth record into x.
!!
      implicit none
!
      class(record_storer) :: storer
!
      real(dp), dimension(storer%record_dim), intent(out) :: x
!
      integer, intent(in) :: n
!
      if (.not. storer%in_memory) then
!
!        If record is on file, read record from file
!
         call storer%file_%open_('read')
!
         call storer%file_%read_(x, n, n)
!
         call storer%file_%close_()
!
      else
!
!        If records are in memory, copy the record
!
         call dcopy(storer%record_dim,    &
                    storer%array_p(:, n), &
                    1,                    &
                    x,                    &
                    1)
!
      endif
!
   end subroutine copy_record_out_single_record
!
!
   subroutine copy_record_out_range(storer, x, a_range)
!!
!!    Copy record out
!!    Written by Eirik F. Kjønstad, 2020
!!
!!    Copies an range of records into x.
!!
      implicit none
!
      class(record_storer) :: storer
!
      class(range_), intent(in) :: a_range
!
      real(dp), dimension(1), intent(out) :: x
!
      if (.not. storer%in_memory) then
!
!        If on file, read from file
!
         call storer%file_%open_('read')
!
         call storer%file_%read_(x, a_range%first, a_range%get_last())
!
         call storer%file_%close_()
!
      else
!
!        If records are in memory, copy the record
!
         call dcopy(storer%record_dim*a_range%length, &
                    storer%array_p(:, a_range%first), &
                    1, x, 1)
!
      endif
!
   end subroutine copy_record_out_range
!
!
   subroutine load_records(storer, x, first, last, block_)
!!
!!    Load records
!!    Written by Eirik F. Kjønstad, 2020
!!
!!    Specific procedure for the generic "load".
!!
!!    x:         pointer to load records into
!!
!!    first:     first record to load into x
!!    last:      last record to load into x
!!
!!    block_:    (optional, integer) x is pointed to this
!!               block in memory. Default is 1.
!!
!!               The blocks are ordered in the same order as "batches" in the
!!               prepare_records routine. E.g., if [batch_i, batch_j] was sent to
!!               prepare_records, then the i records should be loaded with block=1
!!               and the j records with block=2.
!!
!!               When loading records from file, we will read from file only if the
!!               requested block is not already loaded into an earlier block. E.g.,
!!               for a loop over i and j (corresponding to blocks 1 and 2), then we
!!               check whether j records are already loaded into the i block. If they
!!               are, then x is set to point to the i block. Similarly, for i, j, and
!!               k (blocks 1, 2, and 3), j can point to i blocks and k can point to
!!               i and j blocks. This is done to avoid incorrectly overwriting
!!               loaded pointers.
!!
!!               If the records are in memory, the value of "block" will not
!!               be used. In this case x is just set to the specified range.
!!
      implicit none
!
      class(record_storer) :: storer
!
      integer, intent(in) :: first, last
!
      real(dp), dimension(:,:), pointer, contiguous :: x
!
      integer, intent(in), optional :: block_
!
      integer :: block_local, block_length, length
!
      integer :: stored_in_block, k, first_offset, last_offset
!
      block_local = 1
      if (present(block_)) block_local = block_
!
!     Sanity check: is range valid?
!
      length = last - first + 1
!
      if (length .lt. 0) then
!
         call output%error_msg('Tried to load invalid range [(i0),(i0)]', ints=[first, last])
!
      endif
!
      if (.not. storer%in_memory) then
!
!        Check if the records are already stored in another block
!
         stored_in_block = -1
!
         do k = 1, storer%n_blocks
!
            if (storer%blocks(k)%contains_(range_(first, last))) stored_in_block = k
!
         enddo
!
         if (stored_in_block .ne. -1 .and. stored_in_block .le. block_local) then
!
!           If records are already loaded into an earlier block,
!           we set the pointer to that block
!
            first_offset = first - storer%blocks(stored_in_block)%first + 1
            last_offset  = last  - storer%blocks(stored_in_block)%first + 1
!
            call storer%linked_arrays%get_array(x, stored_in_block, first_offset, last_offset)
!
         else
!
!           Otherwise, we load into the specified block
!           Get the pointer to the block and read records into it
!
            call storer%linked_arrays%get_array(x, block_local, n_columns=block_length)
!
!           Sanity check: does the specified length exceed the block length?
!
            if (length .gt. block_length) then
!
               call output%error_msg('Tried to read (i0) elements into block nr. (i0), &
                                     &but that block can only hold (i0) elements.', &
                                     ints=[length, block_local, block_length])
!
            endif
!
!           Read records
!
            call storer%file_%read_(x, first, last)
!
!           Store the block range information in the storer
!
            storer%blocks(block_local) = range_(first, last)
!
         endif
!
      else
!
!        If records are in memory, set x to the requested subset of records
!
         x => storer%array_p(:, first : last)
!
      endif
!
   end subroutine load_records
!
!
   subroutine load_range(storer, x, a_range, block_)
!!
!!    Load range
!!    Written by Eirik F. Kjønstad, 2020
!!
!!    Specific procedure for the generic "load".
!!
!!    x:       pointer to load records into
!!
!!    a_range: range of records to load into x
!!
!!    block_:  (optional, integer) x is pointed to this
!!             block in memory. Default is 1.
!!
!!             The blocks are ordered in the same order as "batches" in the
!!             prepare_records routine. E.g., if [batch_i, batch_j] was sent to
!!             prepare_records, then the i records should be loaded with block=1
!!             and the j records with block=2.
!!
!!             If the records are in memory, the value of "block" will not
!!             be used. In this case x is just set to the specified range.
!!
      implicit none
!
      class(record_storer) :: storer
!
      class(range_), intent(in) :: a_range
!
      real(dp), dimension(:,:), pointer, contiguous :: x
!
      integer, intent(in), optional :: block_
!
      call storer%load_records(x, a_range%first, a_range%get_last(), block_)
!
   end subroutine load_range
!
!
   subroutine store_record_storer(storer, x, a_range)
!!
!!    Store
!!    Written by Eirik F. Kjønstad, 2020
!!
!!    Stores changes made to records temporarily held in memory.
!!
!!    x:       pointer to records
!!    a_range: range of records associated with x
!!
!!    If records are kept in memory, the changes made to
!!    the loaded pointer is permanent and this routine does not
!!    have any effect.
!!
!!    If records are kept on file, the changes made to
!!    the loaded pointer is lost when a new set of records
!!    is loaded. This routine stores the records to file
!!    so that they are not lost when a new set of records
!!    are loaded into memory.
!!
      implicit none
!
      class(record_storer) :: storer
!
      class(range_), intent(in) :: a_range
!
      real(dp), dimension(:,:), pointer, contiguous :: x
!
      if (.not. storer%in_memory) then
!
!        Write to file
!
         call storer%file_%write_(x,                  &
                                  a_range%first,      &
                                  a_range%get_last())
!
!        Make sure blocks in memory are updated
!
         call storer%update_loaded_blocks(x, a_range%first, a_range%get_last())
!
!
      endif
!
   end subroutine store_record_storer
!
!
   subroutine store_single_record(storer, x, n)
!!
!!    Store single record
!!    Written by Eirik F. Kjønstad, 2020
!!
!!    Stores changes made to a given record.
!!
!!    If records are kept in memory, the changes made to
!!    the loaded pointer is permanent and this routine does not
!!    have any effect.
!!
!!    If records are kept on file, the changes made to
!!    the loaded pointer is lost when a new set of records
!!    is loaded. This routine stores the records to file
!!    so that they are not lost when a new set of records
!!    are loaded into memory.
!!
      implicit none
!
      class(record_storer) :: storer
!
      integer, intent(in) :: n
!
      real(dp), dimension(:,:) :: x
!
      if (.not. storer%in_memory) then
!
!        Write record to file
!
         call storer%file_%write_(x, n, n)
!
!        Make sure blocks in memory are updated
!
         call storer%update_loaded_blocks(x, n, n)
!
      endif
!
   end subroutine store_single_record
!
!
   subroutine copy_record_in_range(storer, x, a_range)
!!
!!    Set range
!!    Written by Eirik F. Kjønstad, 2020
!!
!!    Saves array x in the storer by writing to file or copying the array.
!!
      implicit none
!
      class(record_storer) :: storer
!
      class(range_), intent(in) :: a_range
!
      real(dp), dimension(1,1), intent(in) :: x
!
      if (.not. storer%in_memory) then
!
!        Save records to file
!
         call storer%file_%open_('write')
!
         call storer%file_%write_(x, a_range%first, a_range%get_last())
!
         call storer%file_%close_()
!
!        If range is already loaded into memory, overwrite block as well
!
         call storer%update_loaded_blocks(x, a_range%first, a_range%get_last())
!
      else
!
!        If records are in memory, copy into memory
!
         call dcopy(storer%record_dim*a_range%length,                      &
                    x,                                                     &
                    1,                                                     &
                    storer%array_p(:, a_range%first : a_range%get_last()), &
                    1)
!
      endif
!
   end subroutine copy_record_in_range
!
!
   subroutine copy_record_in_single_record(storer, x, n)
!!
!!    Set copy single record
!!    Written by Eirik F. Kjønstad, 2020
!!
!!    Saves array x in the storer by writing to file or copying the array.
!!
      implicit none
!
      class(record_storer) :: storer
!
      integer, intent(in) :: n
!
      real(dp), dimension(storer%record_dim), intent(in) :: x
!
      if (.not. storer%in_memory) then
!
!        If records are on file, save records on file
!
         call storer%file_%open_('write')
!
         call storer%file_%write_(x, n, n)
!
         call storer%file_%close_()
!
!        If record is loaded in memory, overwrite record
!
         call storer%update_loaded_blocks(x, n, n)
!
      else
!
!        If records are in memory, copy into memory
!
         call dcopy(storer%record_dim,       &
                    x,                       &
                    1,                       &
                    storer%array_p(:, n),    &
                    1)
!
      endif
!
   end subroutine copy_record_in_single_record
!
!
   subroutine update_loaded_blocks(storer, x, first, last)
!!
!!    Update loaded blocks
!!    Written by Eirik F. Kjønstad, Apr 2020
!!
!!    Copies x into loaded blocks corresponding to the range given by first and last.
!!
!!    Assumes that storer%in_memory is false and should only be called in this case.
!!
      implicit none
!
      class(record_storer) :: storer
!
      integer, intent(in) :: first, last
!
      real(dp), dimension(storer%record_dim, last - first + 1), intent(in) :: x
!
      real(dp), dimension(:,:), pointer, contiguous :: array_p
!
      integer :: k, n, length, first_offset, last_offset
!
      do k = 1, storer%n_blocks
!
!        Skip to next block if [first, last] has no elements in common with the block
!
         if (.not. storer%blocks(k)%contains_(range_(first, last))) cycle
!
         if (storer%blocks(k)%contains_(range_(first, last))) then
!
!           range is subset of the kth block: copy entire block
!
            length = last - first + 1
!
            first_offset = first - storer%blocks(k)%first + 1
            last_offset  = last  - storer%blocks(k)%first + 1
!
            call storer%linked_arrays%get_array(array_p, k, first_offset, last_offset)
!
            call dcopy(storer%record_dim*length,   &
                       x,                          &
                       1,                          &
                       array_p,                    &
                       1)
!
         else
!
!           range is not subset: copy elements one-by-one
!
            call storer%linked_arrays%get_array(array_p, k)
!
            do n = first, last
!
               if (storer%blocks(k)%contains_(n)) then
!
                  call dcopy(storer%record_dim,                            &
                             x(:, n - first + 1),                          &
                             1,                                            &
                             array_p(:, n - storer%blocks(k)%first + 1),   &
                             1)
!
               endif
!
            enddo
!
         endif
!
      enddo
!
   end subroutine update_loaded_blocks
!
!
   subroutine initialize_record_storer(storer)
!!
!!    Initialize
!!    Written by Eirik F. Kjønstad, Nov 2019
!!
      implicit none
!
      class(record_storer) :: storer
!
      call output%printf('debug', 'Doing preparations for storer (a0)', &
                         chars=[storer%name_], fs='(/t3,a)')
!
!     Allocate file and linked-list object for storing records in memory
!
      storer%file_         = direct_stream_file(storer%name_, storer%record_dim)
      storer%linked_arrays = array_list()
!
      if (storer%in_memory) then
!
!        If in memory, allocate a single array in the linked list
!        to hold all records, and set pointer to this array
!
         call storer%linked_arrays%push_back(storer%record_dim, &
                                             storer%n_records)
!
         call storer%linked_arrays%get_array(storer%array_p, 1)
!
      endif
!
   end subroutine initialize_record_storer
!
!
   subroutine finalize_record_storer(storer)
!!
!!    Finalize
!!    Written by Eirik F. Kjønstad, Nov 2019
!!
      implicit none
!
      class(record_storer) :: storer
!
      call output%printf('debug', 'Doing finalizations for storer (a0)', &
                         chars=[storer%name_], fs='(/t3,a)')
!
      if (storer%in_memory) then
!
         call storer%linked_arrays%finalize()
!
      else
!
         call storer%file_%delete_()
!
      endif
!
   end subroutine finalize_record_storer
!
!
end module record_storer_class
