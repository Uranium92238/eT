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
!  :::::::::::::::::::::::::
!  -::- Index functions -::-
!  :::::::::::::::::::::::::
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
!     Debug sanity check 
!
!       if (p .eq. 0 .or. q .eq. 0 .or. r .eq. 0) write(unit_output,*) 'WARNING: one of the indices in index_three is zero!',p,q,r
! !
!       if (p .gt. dim_p) write(unit_output,*) 'WARNING: first index exceeds its dimension', p, dim_p
!       if (q .gt. dim_q) write(unit_output,*) 'WARNING: first index exceeds its dimension', q, dim_q
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
! !     Debug sanity check 
! !
!       if (p .eq. 0 .or. q .eq. 0) write(unit_output,*) 'WARNING: one of the indices in index_two is zero!',p,q
! !
!       if (p .gt. dim_p) write(unit_output,*) 'WARNING: first index exceeds its dimension', p, dim_p
!
   end function index_two
!
! ::::::::::::::::::::::::
! -:- Matrix utilities -:-
! ::::::::::::::::::::::::
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
   subroutine squareup_to_compound(packed,unpacked,N,M)
!!
!!    Square up packed symmetric matrix
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2017
!!
!!    Squares up to full dimension ((N x N), M) of packed matrices.
!!
      implicit none
!
      integer(i15), intent(in) :: N,M
!
      real(dp), dimension(N*(N+1)/2,M), intent(in) :: packed
      real(dp), dimension(N*N,M)                   :: unpacked
!
      integer(i15) :: i = 0, j = 0
!
      do i = 1, N
         do j = 1, N
            unpacked(index_two(i,j,N), 1:M) = packed(index_packed(i,j), 1:M)
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
!
            ! if (abs(unpacked(i, j) - unpacked(j, i)) .gt. 10D-8) then 
            !    write(unit_output,*) 'WARNING: Attempting to pack non-symmetric matrix'
            !    write(unit_output,*) 'Make sure code is bug-free. Information will be lost.'
            !    write(unit_output,*) unpacked(i, j), unpacked(j, i)
            ! endif
!
            packed(index_packed(i, j), 1) = unpacked(i, j)
!
         enddo
      enddo
!
   end subroutine
!
! :::::::::::::::::::::::::::::::
!  -:- Batching functionality -:-
! :::::::::::::::::::::::::::::::
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
   subroutine get_n_lowest(n, size, vec, sorted_short_vec, index_list)
!!
!!    Get n lowest elements
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!    Finds the n lowest values of vec,
!!    sorts them, and returns them in sorted_short_vec 
!!    together with an index list refering to the indices of the 
!!    lowest elements in the original vector.
!!
      implicit none
!
      integer(i15) :: n    ! Number of elements wanted
      integer(i15) :: size ! Size of original vector
!
      real(dp), dimension(size, 1) :: vec
      real(dp), dimension(n, 1)    :: sorted_short_vec
!
      integer(i15), dimension(n, 1) ::index_list
!
!     Variables for sorting
!
      real(dp)     :: max
      integer(i15) :: max_pos
!
      real(dp)     :: swap     = zero
      integer(i15) :: swap_int = 0
!
      integer(i15) :: i = 0, j = 0
!
!        Placing the n first elements of vec into sorted_short_vec
!
         sorted_short_vec(1,1) = vec(1,1)
         index_list(1,1) = 1
!
         max = sorted_short_vec(1,1)
         max_pos = 1
!
         do i = 2, n
!
            sorted_short_vec(i,1) = vec(i,1)
            index_list(i,1) = i
!
            if (sorted_short_vec(i,1) .ge. max) then
!
               max = sorted_short_vec(i,1)
               max_pos = i
!
            endif
         enddo
!
!        Looping through the rest of vec to find lowest values
!
         do i = n + 1, size
            if (vec(i,1) .lt. max) then
!
               sorted_short_vec(max_pos,1) = vec(i,1)
               index_list(max_pos,1) = i
               max = vec(i,1)
!
               do j = 1, n
                  if (sorted_short_vec(j, 1) .gt. max) then
!
                     max = sorted_short_vec(j, 1)
                     max_pos = j
!
                  endif
               enddo
            endif
         enddo
!
!        Sorting sorted_short_vec
!
         do i = 1, n
            do j = 1, n - 1
               if (sorted_short_vec(j,1) .gt. sorted_short_vec(j+1, 1)) then
!
                  swap = sorted_short_vec(j,1)
                  sorted_short_vec(j,1) = sorted_short_vec(j+1, 1)
                  sorted_short_vec(j+1, 1) = swap
!
                  swap_int = index_list(j, 1)
                  index_list(j,1) = index_list(j + 1,1)
                  index_list(j + 1,1) = swap_int
!
               endif
            enddo
         enddo     
!
   end subroutine get_n_lowest
!
   function check_orthogonality(A, M, N)
!!
!!   Check orthogonality
!!   Written by Eirik F. Kjønstad and Sarai D. Folkestad, June 2017
!!
!!    Check if columns of A are orthogonal. A is (M x N) matrix.
!!    Returns logical.
!!
      use workspace
!
      implicit none
!  
      integer(i15)             :: M
      integer(i15)             :: N
      real(dp), dimension(M,N) :: A
      logical                  :: check_orthogonality
!
      integer(i15) :: i = 0, j = 0
      real(dp), dimension(:,:), allocatable :: a_i, a_j
      real(dp) :: ddot
!
      check_orthogonality = .true.
!
      call allocator(a_i, M, 1)
      call allocator(a_j, M, 1)
!
      do i = 1, N
         a_i(:,1) = A(:,i)
         do j = 1, i-1
            a_j(:,1) = A(:,j)
            if (abs(ddot(M,a_i, 1, a_j, 1)) .gt. 1.0d-07) then
               check_orthogonality = .false.
               return
            endif
         enddo
      enddo
!
      call deallocator(a_i, M, 1)
      call deallocator(a_j, M, 1)
!
   end function check_orthogonality
! 
end module utils
