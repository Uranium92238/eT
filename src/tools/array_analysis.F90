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
   use kinds
   use parameters
   use memory_manager_class
!
   implicit none
!
contains
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
!!    index_list(K) will remain 0 for K > dim_x.
!!
!!    Note, EFK June 2018. There is a bug in this routine if two elements
!!    are identical. Use get_n_lowest routine instead, if possible.
!!
      implicit none
!
      real(dp), dimension(:) :: x
      integer             :: dim_x
!
      integer :: n_to_sort ! Number of elements of x to sort
!
      integer, dimension(n_to_sort) :: index_list ! List of the sorted indices (p, q)
!
      integer :: counter = 0, I = 0
      real(dp)     :: largest
!
      do counter = 1, min(n_to_sort, dim_x)
!
         largest = abs(x(1))
!
         do I = counter + 1, dim_x
!
            if (counter .eq. 1) then
!
                if (abs(x(I)) .gt. largest) then
!
                   index_list(counter) = I ! Set new maximum index
                   largest = abs(x(I))      ! Set new maximum value
!
               endif
!
            else
!
                if (abs(x(I)) .gt. largest                            .and. &
                    abs(x(I)) .le. abs(x(index_list(counter-1))) .and. &
                              I .ne. index_list(counter-1)) then
!
                   index_list(counter) = I ! Set new maximum index
                   largest = abs(x(I))      ! Set new maximum value
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
      integer :: n    ! Number of elements wanted
      integer :: size ! Size of original vector
!
      real(dp), dimension(size) :: vec
      real(dp), dimension(n)    :: sorted_short_vec
!
      integer, dimension(n) :: index_list
!
!     Variables for sorting
!
      real(dp)     :: max
      integer :: max_pos
!
      real(dp)     :: swap     = zero
      integer :: swap_int = 0
!
      integer :: i = 0, j = 0
!
!        Placing the n first elements of vec into sorted_short_vec
!
         sorted_short_vec(1) = vec(1)
         index_list(1) = 1
!
         max = sorted_short_vec(1)
         max_pos = 1
!
         do i = 2, n
!
            sorted_short_vec(i) = vec(i)
            index_list(i) = i
!
            if (sorted_short_vec(i) .ge. max) then
!
               max = sorted_short_vec(i)
               max_pos = i
!
            endif
         enddo
!
!        Looping through the rest of vec to find lowest values
!
         do i = n + 1, size
            if (vec(i) .lt. max) then
!
               sorted_short_vec(max_pos) = vec(i)
               index_list(max_pos) = i
               max = vec(i)
!
               do j = 1, n
                  if (sorted_short_vec(j) .gt. max) then
!
                     max = sorted_short_vec(j)
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
               if (sorted_short_vec(j) .gt. sorted_short_vec(j+1)) then
!
                  swap = sorted_short_vec(j)
                  sorted_short_vec(j) = sorted_short_vec(j+1)
                  sorted_short_vec(j+1) = swap
!
                  swap_int = index_list(j)
                  index_list(j) = index_list(j + 1)
                  index_list(j + 1) = swap_int
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
      integer :: n    ! Number of elements wanted
      integer :: size ! Size of original vector
!
      real(dp), dimension(size) :: vec
      real(dp), dimension(n)    :: sorted_short_vec
!
      integer, dimension(n) :: index_list
!
!     Variables for sorting
!
      real(dp)     :: min
      integer :: min_pos
!
      real(dp)     :: swap     = zero
      integer :: swap_int = 0
!
      integer :: i = 0, j = 0
!
!        Placing the n first elements of vec into sorted_short_vec
!
         sorted_short_vec = zero
         sorted_short_vec(1) = vec(1)
         index_list(1) = 1
!
         min = sorted_short_vec(1)
         min_pos = 1
!
         do i = 2, n
!
            sorted_short_vec(i) = vec(i)
            index_list(i) = i
!
            if (sorted_short_vec(i) .le. min) then
!
               min = sorted_short_vec(i)
               min_pos = i
!
            endif
         enddo
!
!        Looping through the rest of vec to find lowest values
!
         do i = n + 1, size
            if (vec(i) .gt. min) then
!
               sorted_short_vec(min_pos) = vec(i)
               index_list(min_pos) = i
               min = vec(i)
!
               do j = 1, n
                  if (sorted_short_vec(j) .lt. min) then
!
                     min = sorted_short_vec(j)
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
               if (sorted_short_vec(j) .lt. sorted_short_vec(j+1)) then
!
                  swap = sorted_short_vec(j)
                  sorted_short_vec(j) = sorted_short_vec(j+1)
                  sorted_short_vec(j+1) = swap
!
                  swap_int = index_list(j)
                  index_list(j) = index_list(j + 1)
                  index_list(j + 1) = swap_int
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
      integer             :: M
      integer             :: N
      real(dp), dimension(M,N) :: A
      logical                  :: check_orthogonality
!
      integer :: i = 0, j = 0
      real(dp), dimension(:), allocatable :: a_i, a_j
      real(dp) :: ddot
!
      check_orthogonality = .true.
!
      call mem%alloc(a_i, M)
      call mem%alloc(a_j, M)
!
      do i = 1, N
         a_i(:) = A(:,i)
         do j = 1, i-1
            a_j(:) = A(:,j)
            if (abs(ddot(M,a_i, 1, a_j, 1)) .gt. 1.0d-07) then
               check_orthogonality = .false.
               return
            endif
         enddo
      enddo
!
      call mem%dealloc(a_i, M)
      call mem%dealloc(a_j, M)
!
   end function check_orthogonality
!
!
   recursive subroutine quicksort_recursive(vec, first, last)
!!
!!    Recursive implementation of quicksort (descending order)
!!
!!    Stolen from https://gist.github.com/t-nissie/479f0f16966925fa29ea
!!
      implicit none
!
      real(dp), dimension(:), intent(inout) :: vec
      integer, intent(in) :: first, last
!
      real(dp) :: pivot, temp
      integer :: i, j

      pivot = vec((first+last)/2)
!
      i = first
      j = last
!
      do
!
         do while (vec(i) > pivot)
!
            i = i + 1
!
         end do
!
         do while (pivot > vec(j))
!
            j = j - 1
!
         end do
!
         if (i >= j) exit
!
         temp = vec(i)
         vec(i) = vec(j)
         vec(j) = temp
!
         i = i + 1
         j = j - 1
!
      end do
!
      if (first < i - 1) call quicksort_recursive(vec, first, i-1)
      if (j + 1 < last)  call quicksort_recursive(vec, j+1, last)
!
   end subroutine quicksort_recursive
!
   recursive subroutine quicksort_with_index_recursive(vec, index_list, first, last)
!!
!!    Recursive implementation of quicksort with index list (descending order)
!!
!!    Adapted from https://gist.github.com/t-nissie/479f0f16966925fa29ea
!!    index_list stuff added by Sarai D. Folkestad
!!
!!
      implicit none
!
      real(dp), dimension(:), intent(inout) :: vec
      integer, dimension(:), intent(inout) :: index_list
      integer, intent(in) :: first, last
!
      real(dp) :: pivot, temp
      integer :: temp_index
      integer :: i, j

      pivot = vec((first+last)/2)
!
      i = first
      j = last
!
      do
!
         do while (vec(i) > pivot)
!
            i = i + 1
!
         end do
!
         do while (pivot > vec(j))
!
            j = j - 1
!
         end do
!
         if (i >= j) exit
!
         temp = vec(i)
         temp_index = index_list(i)
!
         vec(i) = vec(j)
         vec(j) = temp
!
         index_list(i) = index_list(j)
         index_list(j) = temp_index
!
         i = i + 1
         j = j - 1
!
      end do
!
      if (first < i - 1) call quicksort_with_index_recursive(vec, index_list, first, i-1)
      if (j + 1 < last)  call quicksort_with_index_recursive(vec, index_list, j+1, last)
!
   end subroutine quicksort_with_index_recursive
!
   subroutine quicksort_descending(vec, dim)
!!
!!    Wrapper for recursive quicksort routine
!!
      implicit none
!
      integer, intent(in) :: dim
      real(dp), dimension(dim), intent(inout) :: vec
!
      call quicksort_recursive(vec, 1, dim)
!
   end subroutine quicksort_descending
!
   subroutine quicksort_with_index_descending(vec, index_list, dim)
!!
!!    Wrapper for recursive quicksort with index list
!!
      implicit none
!
      integer, intent(in) :: dim
      real(dp), dimension(dim), intent(inout) :: vec
      integer, dimension(dim), intent(inout) :: index_list
!
      integer :: i
!
      do i = 1, dim
         index_list(i) = i
      enddo
!
      call quicksort_with_index_recursive(vec, index_list, 1, dim)
!
   end subroutine quicksort_with_index_descending
!
   subroutine quicksort_with_index_ascending(vec, index_list, dim)
!!
!!    Wrapper for recursive quicksort with index list ascending order
!!
      implicit none
!
      integer, intent(in) :: dim
      real(dp), dimension(dim), intent(inout) :: vec
      integer, dimension(dim), intent(inout) :: index_list
!
      call dscal(dim, -one, vec, 1)
!
      call quicksort_with_index_descending(vec, index_list, dim)
!
      call dscal(dim, -one, vec, 1)
!
   end subroutine quicksort_with_index_ascending
!
!
   recursive subroutine quicksort_recursive_int(vec, first, last)
!!
!!    Recursive implementation of quicksort (descending order)
!!
!!    Stolen from https://gist.github.com/t-nissie/479f0f16966925fa29ea
!!
      implicit none
!
      integer, dimension(:), intent(inout) :: vec
      integer, intent(in) :: first, last
!
      integer :: pivot, temp
      integer :: i, j

      pivot = vec((first+last)/2)
!
      i = first
      j = last
!
      do
!
         do while (vec(i) > pivot)
!
            i = i + 1
!
         end do
!
         do while (pivot > vec(j))
!
            j = j - 1
!
         end do
!
         if (i >= j) exit
!
         temp = vec(i)
         vec(i) = vec(j)
         vec(j) = temp
!
         i = i + 1
         j = j - 1
!
      end do
!
      if (first < i - 1) call quicksort_recursive_int(vec, first, i-1)
      if (j + 1 < last)  call quicksort_recursive_int(vec, j+1, last)
!
   end subroutine quicksort_recursive_int
!
   recursive subroutine quicksort_with_index_recursive_int(vec, index_list, first, last)
!!
!!    Recursive implementation of quicksort with index list (descending order)
!!
!!    Adapted from https://gist.github.com/t-nissie/479f0f16966925fa29ea
!!    index_list stuff added by Sarai D. Folkestad
!!
!!
      implicit none
!
      integer, dimension(:), intent(inout) :: vec
      integer, dimension(:), intent(inout) :: index_list
      integer, intent(in) :: first, last
!
      integer :: pivot, temp
      integer :: temp_index
      integer :: i, j

      pivot = vec((first+last)/2)
!
      i = first
      j = last
!
      do
!
         do while (vec(i) > pivot)
!
            i = i + 1
!
         end do
!
         do while (pivot > vec(j))
!
            j = j - 1
!
         end do
!
         if (i >= j) exit
!
         temp = vec(i)
         temp_index = index_list(i)
!
         vec(i) = vec(j)
         vec(j) = temp
!
         index_list(i) = index_list(j)
         index_list(j) = temp_index
!
         i = i + 1
         j = j - 1
!
      end do
!
      if (first < i - 1) call quicksort_with_index_recursive_int(vec, index_list, first, i-1)
      if (j + 1 < last)  call quicksort_with_index_recursive_int(vec, index_list, j+1, last)
!
   end subroutine quicksort_with_index_recursive_int
!
   subroutine quicksort_descending_int(vec, dim)
!!
!!    Wrapper for recursive quicksort routine
!!
      implicit none
!
      integer, intent(in) :: dim
      integer, dimension(dim), intent(inout) :: vec
!
      call quicksort_recursive_int(vec, 1, dim)
!
   end subroutine quicksort_descending_int
!
   subroutine quicksort_with_index_descending_int(vec, index_list, dim)
!!
!!    Wrapper for recursive quicksort with index list
!!
      implicit none
!
      integer, intent(in) :: dim
      integer, dimension(dim), intent(inout) :: vec
      integer, dimension(dim), intent(inout) :: index_list
!
      integer :: i
!
      do i = 1, dim
         index_list(i) = i
      enddo
!
      call quicksort_with_index_recursive_int(vec, index_list, 1, dim)
!
   end subroutine quicksort_with_index_descending_int
!
   subroutine quicksort_with_index_ascending_int(vec, index_list, dim)
!!
!!    Wrapper for recursive quicksort with index list ascending order
!!
      implicit none
!
      integer, intent(in) :: dim
      integer, dimension(dim), intent(inout) :: vec
      integer, dimension(dim), intent(inout) :: index_list
!
      vec = -1*vec
!
      call quicksort_with_index_descending_int(vec, index_list, dim)
!
      vec = -1*vec
!
   end subroutine quicksort_with_index_ascending_int
!
end module array_analysis
