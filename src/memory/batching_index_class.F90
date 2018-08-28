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
!!       call batch_a%init(wf%n_v) -> initializes batching object for an index of dimension
!!                                    equal to the number of virtual orbitals
!!
!!       call wf%mem%num_batch(batch_a, required) -> determines the number of batches and
!!                                                   and the largest batching length, saving
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
!!    The case is similar for a two-batching process, which can also be handled
!!    by the wavefunction's memory manager object 'mem'.
!!
!!
!
   use kinds
   use parameters
   use file_class
!
!  ::::::::::::::::::::::::::::::::::::::::::::::::
!  -::- Definition of the batching_index class -::-
!  ::::::::::::::::::::::::::::::::::::::::::::::::
!
   type :: batching_index
!
!     Values relating the limits and length of the current batch
!     (set by determine_limits procedure for a given batch)
!
      integer(i15) :: first  = 0 ! Current first value of index
      integer(i15) :: last   = 0 ! Current last value of index
      integer(i15) :: length = 0 ! Current length of batch (last - first + 1)
!
!     Values relating the the size of the batch and the total number of batches
!     (set by memory manager routines)
!
      integer(i15) :: max_length = 0  ! Maximum length of batch (most batches will be of this size, but typically not all)
      integer(i15) :: num_batches = 0 ! The number of batches in total for the index
!
!     Value that must be initialized by user
!
      integer(i15) :: index_dimension = 0 ! Full length of index (e.g., typically n_vir for virtual index)
!
!     Logical for initialization (for sanity check)
!
      logical :: initialized = .false.
!
   contains
!
!     Initialization routine (sets the index dimension)
!
      procedure :: init => init_batching_index
!
!     Routine that sets the batch dependent variables,
!     first, last and length, based on which batch it is
!
      procedure :: determine_limits => determine_limits_batching_index
!
   end type batching_index
!
!
contains
!
!
   subroutine init_batching_index(batch_p, dimension)
!!
!!    Init (batching index)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Dec 2017
!!
!!    Note: every batching index must be initialized!
!!    The 'dimension' variable specifies the total length of the
!!    batching index, e.g. the number of virtuals for a virtual index.
!!
      implicit none
!
      class(batching_index) :: batch_p
!
      integer(i15), intent(in) :: dimension
!
      batch_p%index_dimension = dimension
      batch_p%initialized = .true.
!
   end subroutine init_batching_index
!
!
   subroutine determine_limits_batching_index(batch_p, batch_number)
!!
!!    Determine limits
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Dec 2017
!!
!!    Given the batch number, this routine determine the first and
!!    and last values, for the index, as well as the length of the
!!    current batching interval.
!!
      implicit none
!
      class(batching_index) :: batch_p ! p is general index (can be virtual or occupied or other)
!
      integer(i15), intent(in) :: batch_number ! The current batch
!
!     Sanity check
!
      if (.not. batch_p%initialized) then
!
         write(output%unit,'(t3,a)') 'Error: a non-initialized batching variable was used.'
         stop
!
      endif
!
!     Determine limits of batch, q = first, first + 1, ..., last
!
      batch_p%first = 1 + (batch_number-1)*(batch_p%max_length)
      batch_p%last  = min((batch_p%max_length)+(batch_number-1)*(batch_p%max_length), batch_p%index_dimension)
!
!     Calculate the length of the batch
!
      batch_p%length = batch_p%last - batch_p%first + 1
!
   end subroutine determine_limits_batching_index
!
!
end module batching_index_class
