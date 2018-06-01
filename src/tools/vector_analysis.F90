module vector_analysis
!
!!
!!    Vector analysis module
!!    Written by Sarai D. Folkstad and Eirik F. Kjønstad, May 2018
!!
!
   use input_output
   use utils
   use types
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
      write(unit_output,'(t6,a32)')'------------------------------------------------------'
!
      if (present(label_p) .and. present(label_q)) then
!
         write(unit_output,'(t6,a3, 8x, a3, 8x, a10)') label_p, label_q, 'amplitude'
!
      elseif (.not. present(label_p) .and. .not. present(label_q)) then
!
         write(unit_output,'(t6,a3, 8x, a3, 8x, a10)') 'p', 'q', 'amplitude'
!
      else
!
         write(unit_output,'(/t6,a)') 'Error: when printing vector, either specify all or no index labels (p, q)'
         stop
!
      endif
!
      write(unit_output,'(t6,a32)')'------------------------------------------------------'
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
               write(unit_output,'(t6,i3, 8x,i3, 10x, f8.5)') p, q, x_pq(index_list(I,1), 1)
!
            endif
!
         endif
!
      enddo
!
      write(unit_output,'(t6,a32)')'------------------------------------------------------'
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
         write(unit_output,'(/t6, a)') 'Error: four-index printing of dominant elements only supports'
         write(unit_output,'(t6, a)')  'square symmetric matrices'
         stop
!
      endif
!
!     Print banner
!
      write(unit_output,'(t6,a54)')'------------------------------------------------------'
!
      if (present(label_p) .and. present(label_q) .and. present(label_r) .and. present(label_s)) then
!
         write(unit_output,'(t6, a3, 8x, a3, 8x, a3, 8x, a3, 8x, a10)')'a','i','b','j', 'amplitude'
!
      elseif (.not. present(label_p) .and. .not. present(label_q) .and. .not. present(label_r) .and. .not. present(label_s)) then
!
         write(unit_output,'(t6, a3, 8x, a3, 8x, a3, 8x, a3, 8x, a10)') label_p, label_q, label_r, label_s, 'amplitude'
!
      else
!
         write(unit_output,'(/t6,a)') 'Error: when printing vector, either specify all or no index labels (p, q, r, s)'
         stop
!
      endif
!
      write(unit_output,'(t6,a54)')'------------------------------------------------------'
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
                write(unit_output,'(t6,i3, 8x,i3, 8x,i3, 8x, i3, 10x, f8.5)') &
                                    p, q, r, s, x_pqrs(index_list(I,1), 1)
!
            endif
!
         endif
!
      enddo
!
      write(unit_output,'(t6,a54)')'------------------------------------------------------'
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
   subroutine invert_compound_index(pq, p, q, dim_p, dim_q)
!!
!!    Invert compound index
!!    Written by Eirik F. Kjønstad, May 2018
!!
!!    Given the compound index pq = dim_p*(q-1) + p, this routine
!!    determines p and q.
!!
      implicit none
!
      integer(i15), intent(in) :: pq
!
      integer(i15) :: p, q
      integer(i15), intent(in) :: dim_p, dim_q
!
      integer(i15) :: I = 0
!
!     Since dim_p*q >= pq > dim_p*(q-1) by construction, we can determine q
!
      do I = 1, dim_q
!
         if  (pq .gt. dim_p*(I-1) .and. pq .le. dim_p*I) then
!
            q = I
!
         endif
!
      enddo
!
!     From q, it is straight-forward to determine p
!
      p = pq - dim_p*(q-1)
!
!     Sanity test
!
      if (index_two(p, q, dim_p) .ne. pq) then
!
         write(unit_output,'(/t3,a)') 'Error: inversion of compound index unsuccessful'
         stop
!
      endif
!
   end subroutine invert_compound_index
!
!
   subroutine invert_packed_index(pq, p, q, dim)
!!
!!    Invert packed index
!!    Written by Eirik F. Kjønstad, June 2018
!!
!!    Returns the indices (p, q) in the upper triangular part of the symmetric matrix.
!!
      implicit none
!
      integer(i15) :: pq
!
      integer(i15) :: p, q
!
      integer(i15) :: dim
!
      integer(i15) :: I = 0, J = 0
!
!     Loop through upper triangular part until the index matches
!
      do I = 1, dim
         do J = I, dim
!
            if (index_packed(I, J) .eq. pq) then
!
               p = I
               q = J
!
            endif
!
         enddo
      enddo
!
!     Sanity check
!
      if (index_packed(p, q) .ne. pq) then
!
         write(unit_output,'(/t6,a)') 'Error: inversion of packed index unsuccessful'
         stop
!
      endif
!
   end subroutine invert_packed_index
!
!
   real(dp) function dot_product(x, y, n)
!!
!!    Calculate dot product
!!    Written by Eirik F. Kjønstad, June 2018
!!
!!    Returns the dot product of x and y, two vectors of length n
!!
      implicit none
!
      integer(i15), intent(in) :: n
!
      real(dp), dimension(:,:), intent(in) :: x
      real(dp), dimension(:,:), intent(in) :: y
!
      real(dp) :: ddot
!
      dot_product = ddot(n, x, 1, y, 1)
!
   end function dot_product
!
!
end module vector_analysis
