module batching
!!
!!    Batching module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2018
!!
!!    Note, June 2018, EFK. TEMPORARY MODULE. TO BE REMOVED ONCE ALL BATCHING
!!    IS DONE THROUGH THE MEMORY MANAGER USING BATCHING INDEX OBJECTS.
!!
!
   use input_output
   use types
!
   implicit none
!
contains
!
!
   subroutine num_batch(required,available,max_batch_length,n_batch,batch_dimension)
!!
!!    Number of batches
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2017
!!
!!    Calculates number of batches
!!
!!    Batching structure will be:
!!    With rest:     (n_batch-1)*(max_batch_length) + rest = required
!!    Without rest:  (n_batch)*(max_batch_length) = required
!!
      implicit none
!
!
      integer(i15), intent(in)           :: available, batch_dimension
      integer(i15)                       :: max_batch_length,n_batch
      integer(i15)                       :: required
      integer(i15)                       :: buffer
!
!     Adding buffer for required
!
      buffer = required/10
!
      required = required + buffer
!
      if (required .lt. available) then
         n_batch = 1
         max_batch_length = batch_dimension
         return
      endif
!
!     Max batch size
!
      max_batch_length = available/(required/batch_dimension)
!
!     Number of full batches
!
      n_batch=batch_dimension/max_batch_length
!
!     Test for rest
!
      if (n_batch*max_batch_length .lt. batch_dimension) then
         n_batch = n_batch+1
      endif
!
   end subroutine num_batch
!
!
   subroutine num_two_batch(required,available,max_batch_length,n_batch,batch_dimension)
!!
!!    Number of batches
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2017
!!
!!    Calculates number of batches when two batching variables are needed
!!
!!    Batching structure will be:
!!    With rest:     (n_batch-1)*(max_batch_length) + rest = required
!!    Without rest:  (n_batch)*(max_batch_length) = required
!!
      implicit none
!
      integer(i15), intent(in) :: available, batch_dimension
      integer(i15)             :: max_batch_length,n_batch,i,buffer,required
!
      buffer = required/10
!
      required = required + buffer
!
      n_batch = 1
!
      if (required .lt. available) then
            n_batch = 1
            max_batch_length = batch_dimension
            return
      endif
!
      do i = 1, batch_dimension
         if (available .gt. required/i**2) then
!
            n_batch = i
            max_batch_length = batch_dimension/n_batch
!
!           Test for rest
!
            if (n_batch*max_batch_length .lt. batch_dimension) then
               n_batch = n_batch + 1
            endif
!
            return
         endif
      enddo
!
   end subroutine num_two_batch
!
!
   subroutine batch_limits(first,last,batch_number,max_batch_length,batch_dimension)
!!
!!     Batch limits
!!     Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2017
!!
!!     Find batch limits (first and last)
!!
!!     batch_number: the current batch (1,2,...,n_batch)
!!     max_batch_length: the length of each batch (except the last, which may be a rest, see n_one_batch routine)
!!     batch_dimension: the dimensionality of the batching variable (e.g., n_vir for a virtual index)
!!
      implicit none
!
      integer(i15) :: first,last
      integer(i15), intent(in) :: batch_number,max_batch_length,batch_dimension
!
      first = 1 + (batch_number-1)*max_batch_length
      last  = min(max_batch_length+(batch_number-1)*max_batch_length,batch_dimension)
!
   end subroutine batch_limits
!
!
end module batching
