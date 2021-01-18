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
   use memory_manager_class, only: mem   
   use global_out, only: output  
   use interval_class, only: interval
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
      character(len=255) :: name_ = 'no_name'
!
!     Maximum number of records & dimensionality of each record
!
      integer :: n_records 
      integer :: record_dim
!
!     Linked list of arrays for storing records in memory 
!
      type(array_list), allocatable :: linked_arrays
!
!     Pointer for array (if in_memory = true)
!
      real(dp), dimension(:,:), pointer, contiguous :: array_p
!
!     Are the records stored in memory?
!
      logical :: in_memory 
!
!     Delete records on finalization?
!     Records are dumped to file if records are stored in memory and delete is false
!
      logical :: delete 
!
!     File for storing records
!
      type(direct_stream_file), allocatable :: file_
!
!     Storing intervals for each of the blocks loaded into memory 
!
      integer :: n_blocks 
      type(interval), dimension(:), allocatable :: blocks
!
   contains
!
!     Prepare to hold a set of records in memory
!
      procedure :: prepare_records           => prepare_records_record_storer
!
!     Sets pointer to array where the records are stored 
!
      procedure :: load_records_record_storer
      procedure :: load_interval_record_storer
!
      generic :: load                        => load_records_record_storer, &
                                                load_interval_record_storer
!
!     Free up space allocated to temporarily keep records in memory 
!
      procedure :: free_records              => free_records_record_storer
!
!     Store changes made to records temporarily held in memory 
!
      procedure :: store                     => store_record_storer 
!
!     Set copy: set record(s) by copying into storer 
!
      procedure :: copy_record_in_single_record_record_storer
      procedure :: copy_record_in_interval_record_storer
!
      generic   :: copy_record_in            => copy_record_in_single_record_record_storer, &
                                                copy_record_in_interval_record_storer
!
!     Get copy: get record(s) by copying out of storer  
!
      procedure :: copy_record_out_single_record_record_storer
      procedure :: copy_record_out_interval_record_storer
!
      generic   :: copy_record_out           => copy_record_out_single_record_record_storer, &
                                                copy_record_out_interval_record_storer
!
!     # elements that will be allocated on load of one record 
!
      procedure :: required_to_load_record   => required_to_load_record_record_storer
!
!     Dumps all records on file if they are kept in memory 
!
      procedure :: write_all                 => write_all_record_storer
!
!     When storing/setting, make sure blocks loaded into memory are updated 
!
      procedure :: update_block              => update_block_record_storer
!
!     Initializations and finalizations
!
!     For memory storer, allocates/deallocates full array. 
!     Opens/closes file for file storer.
!
      procedure :: initialize_storer         => initialize_storer_record_storer 
      procedure :: finalize_storer           => finalize_storer_record_storer 
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
                              records_in_mem,   &
                              delete) result(storer)
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
!!    delete:           Delete file if records are stored on file.
!!                      If records are stored in memory, the records are dumped 
!!                      on file when the storer is finalized.
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
      logical, intent(in) :: delete
!
      storer%name_       = trim(name_)
      storer%record_dim  = record_dim 
      storer%n_records   = n_records 
!
      storer%in_memory = records_in_mem
      storer%delete    = delete 
      storer%n_blocks = 0
!
      if (storer%in_memory) then 
!
         call output%printf('v', trim(storer%name_) // &
                                 ' is stored in memory.', fs='(/t3,a)')
!
      else 
!
         call output%printf('v', trim(storer%name_) // &
                                 ' is stored on file.', fs='(/t3,a)')
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
!        Allocate blocks for holding interval information
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
!        Deallocate blocks for holding interval information
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
   subroutine copy_record_out_single_record_record_storer(storer, x, n)
!!
!!    Get copy 
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
   end subroutine copy_record_out_single_record_record_storer
!
!
   subroutine copy_record_out_interval_record_storer(storer, x, interval_)
!!
!!    Get copy 
!!    Written by Eirik F. Kjønstad, 2020 
!!
!!    Copies an interval of records into x.
!!
      implicit none 
!
      class(record_storer) :: storer 
!
      class(interval), intent(in) :: interval_ 
!
      real(dp), dimension(storer%record_dim, interval_%length), intent(out) :: x 
!
      if (.not. storer%in_memory) then 
!
!        If on file, read from file 
!
         call storer%file_%open_('read')
!
         call storer%file_%read_(x, interval_%first, interval_%last)
!
         call storer%file_%close_()
!
      else 
!
!        If records are in memory, copy the record 
!
         call dcopy(storer%record_dim*interval_%length,     &
                    storer%array_p(:, interval_%first),     &
                    1,                                      &
                    x,                                      &
                    1)
!
      endif       
!
   end subroutine copy_record_out_interval_record_storer
!  
!
   subroutine load_records_record_storer(storer, x, first, last, block_)
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
!!               be used. In this case x is just set to the specified interval.
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
!     Sanity check: is interval valid?
!
      length = last - first + 1
!
      if (length .lt. 0) then 
!
         call output%error_msg('Tried to load invalid interval [(i0),(i0)]', ints=[first, last])
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
            if (storer%blocks(k)%is_subset(first, last)) stored_in_block = k
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
!           Store the block interval information in the storer 
!
            storer%blocks(block_local) = interval(first, last)
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
   end subroutine load_records_record_storer
!
!
   subroutine load_interval_record_storer(storer, x, interval_, block_)
!!
!!    Load interval 
!!    Written by Eirik F. Kjønstad, 2020
!!
!!    Specific procedure for the generic "load".
!!
!!    x:         pointer to load records into  
!!
!!    interval_: interval of records to load into x
!!
!!    block_:    (optional, integer) x is pointed to this 
!!               block in memory. Default is 1. 
!!
!!               The blocks are ordered in the same order as "batches" in the 
!!               prepare_records routine. E.g., if [batch_i, batch_j] was sent to 
!!               prepare_records, then the i records should be loaded with block=1  
!!               and the j records with block=2.
!!
!!               If the records are in memory, the value of "block" will not 
!!               be used. In this case x is just set to the specified interval.
!!
      implicit none 
!
      class(record_storer) :: storer
!
      class(interval), intent(in) :: interval_
!
      real(dp), dimension(:,:), pointer, contiguous :: x 
!
      integer, intent(in), optional :: block_
!
      call storer%load(x, interval_%first, interval_%last, block_)
!
   end subroutine load_interval_record_storer
!
!
   subroutine store_record_storer(storer, x, interval_)
!!
!!    Store
!!    Written by Eirik F. Kjønstad, 2020 
!!
!!    Stores changes made to records temporarily held in memory. 
!!
!!    x:         pointer to records 
!!    interval_: interval of records associated with x 
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
      class(interval), intent(in) :: interval_
!
      real(dp), dimension(:,:), pointer, contiguous :: x 
!
      if (.not. storer%in_memory) then 
!
!        Write to file 
!
         call storer%file_%write_(x,                  &
                                  interval_%first,    &
                                  interval_%last)
!
!        Make sure blocks in memory are updated 
!
         call storer%update_block(x, interval_%first, interval_%last)
!
!
      endif
!
   end subroutine store_record_storer
!
!
   subroutine store_single_record_record_storer(storer, x, n)
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
         call storer%update_block(x, n, n)
!
      endif
!
   end subroutine store_single_record_record_storer
!
!
   subroutine copy_record_in_interval_record_storer(storer, x, interval_)
!!
!!    Set interval 
!!    Written by Eirik F. Kjønstad, 2020
!!
!!    Saves array x in the storer by writing to file or copying the array.
!!
      implicit none 
!
      class(record_storer) :: storer
!
      class(interval), intent(in) :: interval_
!
      real(dp), dimension(storer%record_dim, interval_%length), intent(in) :: x 
!
      if (.not. storer%in_memory) then 
!
!        Save records to file 
!
         call storer%file_%open_('write')
!
         call storer%file_%write_(x, interval_%first, interval_%last)
!
         call storer%file_%close_()
!
!        If interval is already loaded into memory, overwrite block as well
!
         call storer%update_block(x, interval_%first, interval_%last)
!
      else 
!
!        If records are in memory, copy into memory 
!
         call dcopy(storer%record_dim*interval_%length,                    &
                    x,                                                     &
                    1,                                                     &
                    storer%array_p(:, interval_%first : interval_%last),   &
                    1)
!
      endif 
!
   end subroutine copy_record_in_interval_record_storer
!
!
   subroutine copy_record_in_single_record_record_storer(storer, x, n)
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
         call storer%update_block(x, n, n)
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
   end subroutine copy_record_in_single_record_record_storer
!
!
   subroutine update_block_record_storer(storer, x, first, last)
!!
!!    Update block  
!!    Written by Eirik F. Kjønstad, Apr 2020
!!
!!    Copies x into loaded blocks corresponding to the interval given by first and last.
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
         if (storer%blocks(k)%empty_intersection(first, last)) cycle
!
         if (storer%blocks(k)%is_subset(first, last)) then 
!
!           Interval is subset of the kth block: copy entire block 
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
!           Interval is not subset: copy elements one-by-one
!
            call storer%linked_arrays%get_array(array_p, k)
!
            do n = first, last 
!
               if (storer%blocks(k)%element_is_member(n)) then 
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
   end subroutine update_block_record_storer
!
!
   subroutine write_all_record_storer(storer)
!!
!!    Write all
!!    Written by Eirik F. Kjønstad, Mar 2020 
!!
!!    Dumps records stored in memory to file. 
!!    Assumes that records are in memory and file is closed.
!!
      implicit none 
!
      class(record_storer) :: storer 
!
      call storer%file_%open_('write')
!
      call storer%file_%write_(storer%array_p, 1, storer%n_records)
!
      call storer%file_%close_()
!
   end subroutine write_all_record_storer
!
!
   subroutine initialize_storer_record_storer(storer)
!!
!!    Initialize storer  
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
   end subroutine initialize_storer_record_storer
!
!
   subroutine finalize_storer_record_storer(storer)
!!
!!    Finalize storer  
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
!        Dump records on file if requested 
!
         if (.not. storer%delete) &
            call storer%write_all()
!
!        Deallocate array to keep records 
!
         call storer%linked_arrays%finalize()
!
      else
!
!        Delete file if requested 
!
         if (storer%delete .and. storer%file_%exists()) call storer%file_%delete_()
!
      endif 
!
   end subroutine finalize_storer_record_storer
!
!
end module record_storer_class
