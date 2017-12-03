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
      integer(i15) :: first       ! Current first value of index 
      integer(i15) :: last        ! Current last value of index
      integer(i15) :: length      ! Current length of batch (last - first + 1)
      integer(i15) :: max_length  ! Maximum length of batch (most batches will be of this size, but typically not all)
!
      integer(i15) :: full_length ! Full length of index (e.g., typically n_vir for virtual index)
!
   contains 
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
      batch_p%last  = min((batch_p%max_length)+(batch_number-1)*(batch_p%max_length), batch_p%full_length)
!
!     Calculate the length of the batch 
!
      batch_p%length = batch_p%last - batch_p%first + 1
!
   end subroutine determine_limits_batching_index
!
!
end module batching_index_class