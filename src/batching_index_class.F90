module batching_index_class
!
!!
!!                    Batching index class module                                 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Dec 2017         
!!
!!    The batching index class represents a single index (e.g. a in g_abcd)
!!    that is being batched over. It contains information relevant to the 
!!    restriction of that index.
!! 
!
!  :::::::::::::::::::::::::::::::::::
!  -::- Modules used by the class -::-
!  :::::::::::::::::::::::::::::::::::
!
   use types
   use input_output
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
      integer(i15) :: first      = 0 ! Current first value of index 
      integer(i15) :: last       = 0 ! Current last value of index
      integer(i15) :: length     = 0 ! Current length of batch (last - first + 1)
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
!!    Note: every batching index has to be initialized - otherwise it
!!          will not work as expected!
!!
      implicit none 
!
      class(batching_index) :: batch_p
!
      integer(i15), intent(in) :: dimension 
!
      batch_p%index_dimension = dimension 
!
   end subroutine init_batching_index
!
!
   subroutine determine_limits_batching_index(batch_p, batch_number)
!!
!!    Determine limits 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Dec 2017
!!
      implicit none 
!
      class(batching_index) :: batch_p ! p is general index (can be virtual or occupied or other)
!
      integer(i15), intent(in) :: batch_number ! The current batch 
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