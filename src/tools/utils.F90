module utils
!
!!
!!    Utilities module 
!!    Written by Sarai D. Folkstad and Eirik F. Kjønstad, 28 Feb 2017
!!
!!    Contains:
!!    
!!    Index functions:
!!       index_packed: Calculates the packed index of symetric matrix
!!       index_two:    Calculates the compound index of two indices
!!       index_three   Calculates the compound index given by three indices 
!!
!!    Matrix utilities:
!!       packin:      packs in symetric matrix
!!       packed_size: Returns size of packed matrix
!!       squeareup:   squares up symmetric matrix
!!
!!    Batching subroutines:
!!        num_batch:     Calculates the number of batches needed.
!!        num_two_batch: Calculates the number of batches needed 
!!                       for to batching variables with equal number of batches. 
!!        batch_limits:  Returns batch start index and batch end index.
!!
! 
   use input_output
   use types
!
contains
!
!
   integer(i15) function index_packed(i,j)
!!
!!    Packed index    
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2017
!!
!!    Calculates the packed index of symetric matrix.
!!
      implicit none
!
      integer(i15), intent(in) :: i,j
!
      index_packed = (max(i,j)*(max(i,j)-3)/2) + i + j
!
   end function index_packed
!
!
   integer(i15) function index_three(p,q,r,dim_p,dim_q)
!!
!!    Three index compound
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2017
!!
!!    Returns the compound index (pqr)
!!
      implicit none
!
      integer(i15), intent(in) :: p, q, r, dim_p, dim_q
!
      index_three = dim_p*(dim_q*(r-1)+q-1)+p
!
   end function index_three
!
!
   integer(i15) function index_two(p,q,dim_p)
!!
!!    Two index compound
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2017
!!
!!    Returns the compound index (pq)
!!
      implicit none
!
      integer(i15), intent(in) :: p,q,dim_p
!
      index_two = dim_p*(q-1)+p
!
   end function index_two
!
!
   integer(i15) function packed_size(N)
!!
!!    Packed size    
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2017
!!
!!    Returns size of packed symmetric matrices
!!    of dimension N x N (triangular elements) 
!!   
      implicit none
!
      integer(i15), intent(in) :: N
!
      packed_size = N*(N+1)/2
!
   end function packed_size
!
!
   subroutine squareup(packed,unpacked,N)
!!
!!    Square up packed symmetric matrix
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2017
!!
!!    Squares up to full dimension (N x N) of packed matrices.
!!
      implicit none
!
      integer(i15), intent(in) :: N
!
      real(dp), dimension(:,:), intent(in) :: packed
      real(dp), dimension(:,:)             :: unpacked
!
      integer(i15) :: i = 0, j = 0
!
      do i = 1, N
         do j = 1, N
            unpacked(i, j) = packed(index_packed(i,j), 1)
         enddo
      enddo
!
   end subroutine
!
!
   subroutine packin(packed,unpacked,N)
!!
!!    Pack in symmetric matrix
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2017
!!
!!    Pack down full square matrix of dimension N x N.
!!
      implicit none
!
      integer(i15), intent(in) :: N
!
      real(dp), dimension(:,:) :: packed
      real(dp), dimension(:,:),intent(in) :: unpacked
!
      integer(i15) :: i = 0, j = 0
!
      do i = 1, N
         do j = 1, N
            packed(index_packed(i, j), 1) = unpacked(i, j)
         enddo
      enddo
!
   end subroutine
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
!  Max batch size
!
      max_batch_length = available/(required/batch_dimension)
!
!  Number of full batches
!
      n_batch=batch_dimension/max_batch_length
!
!  Test for rest
!
      if (n_batch*max_batch_length .lt. batch_dimension) then
         n_batch = n_batch+1
      endif
!
   end subroutine num_batch
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
!     
      integer(i15), intent(in)           :: available, batch_dimension
      integer(i15)                       :: max_batch_length,n_batch,i,buffer,required
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
      if (available .gt. buffer/i**2) then ! E: insert logical for success!
!
         n_batch = i
         max_batch_length = batch_dimension/n_batch
!
!        Test for rest
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
end module utils
