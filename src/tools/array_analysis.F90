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
      real(dp), dimension(:) :: x_pq
!
      integer :: dim_p
      integer :: dim_q
!
      character(len=1), optional :: label_p
      character(len=1), optional :: label_q
!
      integer, dimension(20) :: index_list
!
      integer :: I = 0, p = 0, q = 0
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
         if (index_list(I) .eq. 0) then ! The vector x has less than 20 elements. Ignore.
!
!           Do nothing
!
         else
!
!           Print vector element and index if above a certain threshold
!
            if (abs(x_pq(index_list(I))) .lt. 1.0D-03) then
!
               exit ! Never print very, very small amplitudes
!
            else
!
!              Invert the compound index
!
               call invert_compound_index(index_list(I), p, q, dim_p, dim_q)
!
!              Print record, i.e. the Ith largest element with indices
!
               write(output%unit,'(t6,i3, 8x,i3, 10x, f8.5)') p, q, x_pq(index_list(I))
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
      real(dp), dimension(:) :: x_pqrs
!
      integer :: dim_p
      integer :: dim_q
      integer :: dim_r
      integer :: dim_s
!
      character(len=1), optional :: label_p
      character(len=1), optional :: label_q
      character(len=1), optional :: label_r
      character(len=1), optional :: label_s
!
      integer, dimension(20) :: index_list
!
      integer :: I = 0, p = 0, q = 0, r = 0, s = 0
!
      integer :: pq = 0, rs = 0
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
         if (index_list(I) .eq. 0) then ! The vector x has less than 20 elements. Ignore.
!
!           Do nothing
!
         else
!
!           Print vector element and index if above a certain threshold
!
            if (abs(x_pqrs(index_list(I))) .lt. 1.0D-03) then
!
               exit ! Never print very, very small amplitudes
!
            else
!
!              Invert the packed index
!
               call invert_packed_index(index_list(I), pq, rs, dim_p*dim_q)
!
               call invert_compound_index(pq, p, q, dim_p, dim_q)
               call invert_compound_index(rs, r, s, dim_r, dim_s)
!
!              Print record, i.e. the Ith largest element with indices
!
                write(output%unit,'(t6,i3, 8x,i3, 8x,i3, 8x, i3, 10x, f8.5)') &
                                    p, q, r, s, x_pqrs(index_list(I))
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
