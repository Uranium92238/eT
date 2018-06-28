module array_analysis
!
!!
!!    Array analysis module
!!    Written by Sarai D. Folkstad and Eirik F. Kjønstad, May 2018
!!
!!    This module contains routines that are used to analyze excitation
!!    vectors (i.e., not change them). In particular, routines to print the
!!    dominant excitation vector elements (two-index for X_ai; four-index
!!    for X_aibj), though indices can be general (X_pq, X_pqrs).
!!
!
   use file_class
   use index
   use kinds
   use parameters
   use memory_manager_class
!
   implicit none
!
contains
!
!
   subroutine print_dominant_two_index(x_pq, dim_p, dim_q, label_p, label_q)
!!
!!    Print dominant elements of two index vector
!!    Written by Eirik F. Kjønstad, May 2018
!!
      implicit none
!
      real(dp), dimension(:,:) :: x_pq
!
      integer(i15) :: dim_p
      integer(i15) :: dim_q
!
      character(len=1), optional :: label_p
      character(len=1), optional :: label_q
!
      integer(i15), dimension(20, 1) :: index_list
!
      integer(i15) :: I = 0, p = 0, q = 0
!
!     Print banner
!
      write(output%unit,'(t6,a32)')'------------------------------------------------------'
!
      if (present(label_p) .and. present(label_q)) then
!
         write(output%unit,'(t6,a3, 8x, a3, 8x, a10)') label_p, label_q, 'amplitude'
!
      elseif (.not. present(label_p) .and. .not. present(label_q)) then
!
         write(output%unit,'(t6,a3, 8x, a3, 8x, a10)') 'p', 'q', 'amplitude'
!
      else
!
         write(output%unit,'(/t6,a)') 'Error: when printing vector, either specify all or no index labels (p, q)'
         stop
!
      endif
!
      write(output%unit,'(t6,a32)')'------------------------------------------------------'
!
!     Determine the 20 largest elements of the vector
!
      index_list = 0
      call determine_dominant_elements(x_pq, dim_p*dim_q, index_list, 20)
!
      do I = 1, 20
!
         if (index_list(I,1) .eq. 0) then ! The vector x has less than 20 elements. Ignore.
!
!           Do nothing
!
         else
!
!           Print vector element and index if above a certain threshold
!
            if (abs(x_pq(index_list(I, 1), 1)) .lt. 1.0D-03) then
!
               exit ! Never print very, very small amplitudes
!
            else
!
!              Invert the compound index
!
               call invert_compound_index(index_list(I,1), p, q, dim_p, dim_q)
!
!              Print record, i.e. the Ith largest element with indices
!
               write(output%unit,'(t6,i3, 8x,i3, 10x, f8.5)') p, q, x_pq(index_list(I,1), 1)
!
            endif
!
         endif
!
      enddo
!
      write(output%unit,'(t6,a32)')'------------------------------------------------------'
!
   end subroutine print_dominant_two_index
!
!
   subroutine print_dominant_four_index(x_pqrs, dim_p, dim_q, dim_r, dim_s, &
                                                label_p, label_q, label_r, label_s)
!!
!!    Print dominant elements of packed four index vector
!!    Written by Eirik F. Kjønstad, May 2018
!!
      implicit none
!
      real(dp), dimension(:,:) :: x_pqrs
!
      integer(i15) :: dim_p
      integer(i15) :: dim_q
      integer(i15) :: dim_r
      integer(i15) :: dim_s
!
      character(len=1), optional :: label_p
      character(len=1), optional :: label_q
      character(len=1), optional :: label_r
      character(len=1), optional :: label_s
!
      integer(i15), dimension(20, 1) :: index_list
!
      integer(i15) :: I = 0, p = 0, q = 0, r = 0, s = 0
!
      integer(i15) :: pq = 0, rs = 0
!
!     Sanity check
!
      if (dim_p*dim_q .ne. dim_r*dim_s) then
!
         write(output%unit,'(/t6, a)') 'Error: four-index printing of dominant elements only supports'
         write(output%unit,'(t6, a)')  'square symmetric matrices'
         stop
!
      endif
!
!     Print banner
!
      write(output%unit,'(t6,a54)')'------------------------------------------------------'
!
      if (present(label_p) .and. present(label_q) .and. present(label_r) .and. present(label_s)) then
!
         write(output%unit,'(t6, a3, 8x, a3, 8x, a3, 8x, a3, 8x, a10)')'a','i','b','j', 'amplitude'
!
      elseif (.not. present(label_p) .and. .not. present(label_q) .and. .not. present(label_r) .and. .not. present(label_s)) then
!
         write(output%unit,'(t6, a3, 8x, a3, 8x, a3, 8x, a3, 8x, a10)') label_p, label_q, label_r, label_s, 'amplitude'
!
      else
!
         write(output%unit,'(/t6,a)') 'Error: when printing vector, either specify all or no index labels (p, q, r, s)'
         stop
!
      endif
!
      write(output%unit,'(t6,a54)')'------------------------------------------------------'
!
!     Determine the 20 largest elements of the vector
!
      index_list = 0
      call determine_dominant_elements(x_pqrs, dim_p*dim_q*(dim_p*dim_q+1)/2, index_list, 20)
!
      do I = 1, 20
!
         if (index_list(I,1) .eq. 0) then ! The vector x has less than 20 elements. Ignore.
!
!           Do nothing
!
         else
!
!           Print vector element and index if above a certain threshold
!
            if (abs(x_pqrs(index_list(I, 1), 1)) .lt. 1.0D-03) then
!
               exit ! Never print very, very small amplitudes
!
            else
!
!              Invert the packed index
!
               call invert_packed_index(index_list(I,1), pq, rs, dim_p*dim_q)
!
               call invert_compound_index(pq, p, q, dim_p, dim_q)
               call invert_compound_index(rs, r, s, dim_r, dim_s)
!
!              Print record, i.e. the Ith largest element with indices
!
                write(output%unit,'(t6,i3, 8x,i3, 8x,i3, 8x, i3, 10x, f8.5)') &
                                    p, q, r, s, x_pqrs(index_list(I,1), 1)
!
            endif
!
         endif
!
      enddo
!
      write(output%unit,'(t6,a54)')'------------------------------------------------------'
!
   end subroutine print_dominant_four_index
!
!
   subroutine determine_dominant_elements(x, dim_x, index_list, n_to_sort)
!!
!!    Determine dominant elements in a vector
!!    Written by Eirik F. Kjønstad, May 2018
!!
!!    On exit from this routine, index_list(K, 1) is the index for the Kth
!!    largest (in absolute value) of the x array, K = 1, 2, ..., n_to_sort.
!!
!!    If dim_x < n_to_sort, all elements of X are sorted. Remember to zero
!!    the index_list array before sending it to the routine. If  dim_x < n_to_sort,
!!    index_list(K,1) will remain 0 for K > dim_x.
!!
!!    Note, EFK June 2018. There is a bug in this routine if two elements
!!    are identical. Use get_n_lowest routine instead, if possible.
!!
      implicit none
!
      real(dp), dimension(:,:) :: x
      integer(i15)             :: dim_x
!
      integer(i15) :: n_to_sort ! Number of elements of x to sort
!
      integer(i15), dimension(n_to_sort,1) :: index_list ! List of the sorted indices (p, q)
!
      integer(i15) :: counter = 0, I = 0
      real(dp)     :: largest
!
      do counter = 1, min(n_to_sort, dim_x)
!
         largest = abs(x(1, 1))
!
         do I = counter + 1, dim_x
!
            if (counter .eq. 1) then
!
                if (abs(x(I,1)) .gt. largest) then
!
                   index_list(counter, 1) = I ! Set new maximum index
                   largest = abs(x(I,1))      ! Set new maximum value
!
               endif
!
            else
!
                if (abs(x(I,1)) .gt. largest                            .and. &
                    abs(x(I,1)) .le. abs(x(index_list(counter-1, 1),1)) .and. &
                              I .ne. index_list(counter-1, 1)) then
!
                   index_list(counter, 1) = I ! Set new maximum index
                   largest = abs(x(I,1))      ! Set new maximum value
!
               endif
!
            endif
!
         enddo
!
      enddo
!
   end subroutine determine_dominant_elements
!
!
   subroutine get_n_lowest(n, size, vec, sorted_short_vec, index_list)
!!
!!    Get n lowest elements
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!    Finds the n lowest values of vec, sorts them, and returns them
!!    in sorted_short_vec, together with an index list refering to the
!!    indices of the lowest elements in the original vector.
!!
      implicit none
!
      integer(i15) :: n    ! Number of elements wanted
      integer(i15) :: size ! Size of original vector
!
      real(dp), dimension(size, 1) :: vec
      real(dp), dimension(n, 1)    :: sorted_short_vec
!
      integer(i15), dimension(n, 1) :: index_list
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
!
   subroutine get_n_highest(n, size, vec, sorted_short_vec, index_list)
!!
!!    Get n highest elements
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!    Finds the n highest values of vec, sorts them, and returns them
!!    in sorted_short_vec, together with an index list refering to the
!!    indices of the highest elements in the original vector.
!!
      implicit none
!
      integer(i15) :: n    ! Number of elements wanted
      integer(i15) :: size ! Size of original vector
!
      real(dp), dimension(size, 1) :: vec
      real(dp), dimension(n, 1)    :: sorted_short_vec
!
      integer(i15), dimension(n, 1) :: index_list
!
!     Variables for sorting
!
      real(dp)     :: min
      integer(i15) :: min_pos
!
      real(dp)     :: swap     = zero
      integer(i15) :: swap_int = 0
!
      integer(i15) :: i = 0, j = 0
!
!        Placing the n first elements of vec into sorted_short_vec
!
         sorted_short_vec = zero
         sorted_short_vec(1,1) = vec(1,1)
         index_list(1,1) = 1
!
         min = sorted_short_vec(1,1)
         min_pos = 1
!
         do i = 2, n
!
            sorted_short_vec(i,1) = vec(i,1)
            index_list(i,1) = i
!
            if (sorted_short_vec(i,1) .le. min) then
!
               min = sorted_short_vec(i,1)
               min_pos = i
!
            endif
         enddo
!
!        Looping through the rest of vec to find lowest values
!
         do i = n + 1, size
            if (vec(i,1) .gt. min) then
!
               sorted_short_vec(min_pos,1) = vec(i,1)
               index_list(min_pos,1) = i
               min = vec(i,1)
!
               do j = 1, n
                  if (sorted_short_vec(j, 1) .lt. min) then
!
                     min = sorted_short_vec(j, 1)
                     min_pos = j
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
               if (sorted_short_vec(j,1) .lt. sorted_short_vec(j+1, 1)) then
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
   end subroutine get_n_highest
!
!
   function check_orthogonality(A, M, N)
!!
!!    Check orthogonality
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, June 2017
!!
!!    Check if columns of A are orthogonal. A is (M x N) matrix.
!!    Returns logical.
!!
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
      call mem%alloc(a_i, M, 1)
      call mem%alloc(a_j, M, 1)
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
      call mem%dealloc(a_i, M, 1)
      call mem%dealloc(a_j, M, 1)
!
   end function check_orthogonality
!
!
end module array_analysis
