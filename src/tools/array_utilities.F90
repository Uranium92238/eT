!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2021 the authors of eT
!
!  eT is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  eT is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program. If not, see <https://www.gnu.org/licenses/>.
!
!
module array_utilities
!
!!
!!    Array utilities module
!!
!!    Routines that perform various operations on arrays.
!!    Reordering routines are gathered in a separate module (see reordering.F90).
!!
!
   use parameters
   use global_out, only : output
   use memory_manager_class, only : mem
!
   implicit none
!
   interface sandwich
      procedure :: sandwich_real, &
                   sandwich_complex
   end interface sandwich
!
   interface symmetric_sandwich_right_transposition
      procedure :: symmetric_sandwich_right_transposition_real, &
                   symmetric_sandwich_right_transposition_complex, &
                   symmetric_sandwich_right_transposition_complex_between_real, &
                   symmetric_sandwich_right_transposition_replace
   end interface symmetric_sandwich_right_transposition
!
   interface scale_diagonal
      procedure :: scale_real_diagonal_by_real, &
                   scale_complex_diagonal_by_real, &
                   scale_complex_diagonal_by_complex, &
                   scale_real_4_diagonal_by_real, &
                   scale_complex_4_diagonal_by_real, &
                   scale_complex_4_diagonal_by_complex, &
                   scale_real_4_diagonal_by_real_1324, &
                   scale_complex_4_diagonal_by_real_1324, &
                   scale_complex_4_diagonal_by_complex_1324, &
                   scale_real_packed_4_diagonal_by_real, &
                   scale_complex_packed_4_diagonal_by_real, &
                   scale_complex_packed_4_diagonal_by_complex
   end interface scale_diagonal
!
   interface entrywise_product
      procedure :: entrywise_product_in_place, &
                   entrywise_product_in_place_2dim, &
                   entrywise_product_in_place_3dim, &
                   entrywise_product_in_place_4dim, &
                   entrywise_product_
   end interface entrywise_product
!
contains
!
!
   subroutine zero_array(x, n)
!!
!!    Zero array 
!!    Written by Eirik F. Kjønstad, 2018 
!!
!!    Sets the array x of length n to zero.
!!
      implicit none 
!
      integer, intent(in) :: n 
!
      real(dp), dimension(n), intent(out) :: x 
!
      integer :: I 
!
!$omp parallel do private(I) schedule(static)
      do I = 1, n 
!
         x(I) = zero 
!
      enddo
!$omp end parallel do
!
   end subroutine zero_array
!
!
   subroutine zero_array_complex(X, n)
!!
!!    Zero array 
!!    Written by Eirik F. Kjønstad, 2018
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
      implicit none 
!
      integer, intent(in) :: n 
!
      complex(dp), dimension(n), intent(out) :: X 
!
      integer :: I 
!
!$omp parallel do private(I) schedule(static)
      do I = 1, n 
!
         X(I) = zero_complex
!
      enddo
!$omp end parallel do
!
   end subroutine zero_array_complex
!
!
   subroutine zero_array_int(x, n)
!!
!!    Zero array int
!!    Written by Eirik F. Kjønstad, 2020 
!!
!!    Sets the integer array x of length n to zero.
!!
      implicit none 
!
      integer, intent(in) :: n 
!
      integer, dimension(n), intent(out) :: x 
!
      integer :: I 
!
!$omp parallel do private(I) schedule(static)
      do I = 1, n 
!
         x(I) = 0 
!
      enddo
!$omp end parallel do
!
   end subroutine zero_array_int
!
!
   subroutine identity_array(x, n)
!!
!!    Identity array 
!!    Written by Ida-Marie Hoyvik, Oct 2019 
!!
!!    Sets the array x of dimension nxn to be the identity matrix.
!!
      implicit none 
!
      integer, intent(in) :: n 
!
      real(dp), dimension(n,n), intent(out) :: x 
!
      integer :: i 
!
      call zero_array(x,n*n)
!
!$omp parallel do private(i) schedule(static)
      do i = 1, n 
!
         x(i,i) = one 
!
      enddo
!$omp end parallel do
!
   end subroutine identity_array
!
!
   pure function is_significant(x, n, threshold, screening) result(is_significant_)
!!
!!    Is significant
!!    Written by Eirik F. Kjønstad and Sarai D. Folkstad, June 2018
!!
!!    Returns true if one element of x (of length n) is greater, in absolute value,
!!    than "threshold": abs(x(i)) > threshold.
!!
!!    screening (optional): if a screening vector is passed to the routine,
!!    then the testing criterion is abs(x(i)*screening(i)) > threshold.
!!
      implicit none
!
      integer, intent(in) :: n
!
      real(dp), dimension(n), intent(in)  :: x
      real(dp), dimension(n), intent(in), optional  :: screening
!
      real(dp), intent(in)  :: threshold
!
      integer :: i
!
      logical :: is_significant_
!
      is_significant_ = .false.
!
      if (present(screening)) then
!
         do i = 1, n
!
            if (abs(x(i)*screening(i)) .gt. threshold) then
!
               is_significant_ = .true.
               return
!
            endif
!
         enddo
!
      else
!
         do i = 1, n
!
            if (abs(x(i)) .gt. threshold) then
!
               is_significant_ = .true.
               return
!
            endif
!
         enddo
      endif
!
   end function is_significant
!
!
   function n_significant(x, n, threshold) result(n_significant_)
!!
!!    Number of significant
!!    Written by Eirik F. Kjønstad and Sarai D. Folkstad, June 2018
!!
!!    Returns the number of elements in a vector x (of dimension n) larger, 
!!    in absolute value, than "threshold": #i such that abs(x(i)) > threshold.
!!
      implicit none
!
      integer, intent(in) :: n
!
      real(dp), dimension(n), intent(in)  :: x
!
      real(dp), intent(in)  :: threshold
!
      integer :: n_significant_, i = 0
!
      n_significant_ = 0
!
      do i = 1, n
!
         if (abs(x(i)) .gt. threshold) then
!
            n_significant_ = n_significant_ + 1
!
         endif
!
      enddo
!
   end function n_significant
!
!
   subroutine reduce_vector(vec, vec_reduced, block_firsts, block_significant, n_blocks, dim_, dim_reduced)
!!
!!    Reduce vector
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Cuts the significant blocks out of a vector and places them in a reduced size vector
!!
!!    vec:  vector to cut out blocks from
!!    dim_: dimension of vec 
!!
!!    vec_reduced: vector to place the cut-out blocks into (intent: out)
!!    dim_reduced: dimension of vec_reduced 
!!
!!    n_blocks:          total number of blocks 
!!    block_significant: a logical array, block_significant(i) == .true. means that i is a significant block 
!!    block_firsts:      an integer array, of dimension n_blocks + 1, containing the first index of each significant block;
!!                       the (n_blocks)+1-th element equals (last element of the final block) + 1.
!!
      implicit none
!
      integer, intent(in) :: dim_, dim_reduced, n_blocks
!
      logical, dimension(n_blocks), intent(in) :: block_significant
      integer, dimension(n_blocks + 1), intent(in) :: block_firsts
!
      real(dp), dimension(dim_), intent(in) :: vec
      real(dp), dimension(dim_reduced), intent(out) :: vec_reduced
!
      integer :: block_, current_pos, first, last, size_
!
      current_pos = 1
!
      do block_ = 1, n_blocks
!
         if (block_significant(block_)) then
!
            first = block_firsts(block_)
            last  = block_firsts(block_ + 1) - 1
            size_ = last - first + 1
!
            vec_reduced(current_pos : current_pos + size_ - 1) = vec(first : last)
            current_pos = current_pos + size_
!
         endif
!
      enddo
!
   end subroutine reduce_vector
!
!
   subroutine reduce_vector_int(vec, vec_reduced, block_firsts, block_significant, n_blocks, dim_, dim_reduced)
!!
!!    Reduce vector int
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Cuts the significant blocks out of a vector and places them in a reduced size vector
!!
!!    vec:  vector to cut out blocks from
!!    dim_: dimension of vec 
!!
!!    vec_reduced: vector to place the cut-out blocks into (intent: out)
!!    dim_reduced: dimension of vec_reduced 
!!
!!    n_blocks:          total number of blocks 
!!    block_significant: a logical array, block_significant(i) == .true. means that i is a significant block 
!!    block_firsts:      an integer array, of dimension n_blocks + 1, containing the first index of each significant block;
!!                       the (n_blocks)+1-th element equals (last element of the final block) + 1.
!!
      implicit none
!
      integer, intent(in) :: dim_, dim_reduced, n_blocks
!
      logical, dimension(n_blocks), intent(in) :: block_significant
      integer, dimension(n_blocks + 1), intent(in) :: block_firsts
!
      integer, dimension(dim_), intent(in) :: vec
      integer, dimension(dim_reduced), intent(out) :: vec_reduced
!
      integer :: block_, current_pos, first, last, size_
!
      current_pos = 1
!
      do block_ = 1, n_blocks
!
         if (block_significant(block_)) then
!
            first = block_firsts(block_)
            last  = block_firsts(block_ + 1) - 1
            size_ = last - first + 1
!
            vec_reduced(current_pos : current_pos + size_ - 1) = vec(first : last)
            current_pos = current_pos + size_
!
         endif
!
      enddo
!
   end subroutine reduce_vector_int
!
!
   subroutine reduce_array(array, array_reduced, block_firsts, block_significant, n_blocks, dim_, dim_reduced, columns)
!!
!!    Reduce array
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Cuts the significant row blocks out of a array and places them in a reduced size array 
!!
!!    array:   array to cut out blocks from
!!    dim_:    number of rows in array 
!!
!!    array_reduced: array to place the cut-out blocks into (intent: out)
!!    dim_reduced:   number of rows in array_reduced
!!
!!    columns: number of columns in array and array_reduced 
!!
!!    n_blocks:          total number of blocks 
!!    block_significant: a logical array, block_significant(i) == .true. means that i is a significant block 
!!    block_firsts:      an integer array, of dimension n_blocks + 1, containing the first index of each significant block;
!!                       the (n_blocks)+1-th element equals (last element of the final block) + 1.
!!
      implicit none
!
      integer, intent(in) :: dim_, dim_reduced, n_blocks, columns
!
      logical, dimension(n_blocks), intent(in) :: block_significant
      integer, dimension(n_blocks + 1), intent(in) :: block_firsts
!
      real(dp), dimension(dim_, columns), intent(in) :: array
      real(dp), dimension(dim_reduced, columns), intent(out) :: array_reduced
!
      integer :: block_, current_pos, first, last, size_, I
!
!$omp parallel do schedule(static) private(I, current_pos, block_, first, last, size_)
      do I = 1, columns
!
         current_pos = 1
!
         do block_ = 1, n_blocks
!
            if (block_significant(block_)) then
!
               first = block_firsts(block_)
               last  = block_firsts(block_ + 1) - 1
               size_ = last - first + 1
!
               array_reduced(current_pos : current_pos + size_ - 1, I) = array(first : last, I)
!
               current_pos = current_pos + size_
!
            endif
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine reduce_array
!
!
   subroutine reduce_array_int(array, array_reduced, block_firsts, block_significant, n_blocks, dim_, dim_reduced, columns)
!!
!!    Reduce array
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Cuts the significant row blocks out of a array and places them in a reduced size array
!!
!!    array:   array to cut out blocks from
!!    dim_:    number of rows in array 
!!
!!    array_reduced: array to place the cut-out blocks into (intent: out)
!!    dim_reduced:   number of rows in array_reduced
!!
!!    columns: number of columns in array and array_reduced 
!!
!!    n_blocks:          total number of blocks 
!!    block_significant: a logical array, block_significant(i) == .true. means that i is a significant block 
!!    block_firsts:      an integer array, of dimension n_blocks + 1, containing the first index of each significant block;
!!                       the (n_blocks)+1-th element equals (last element of the final block) + 1.
!!
      implicit none
!
      integer, intent(in) :: dim_, dim_reduced, n_blocks, columns
!
      logical, dimension(n_blocks), intent(in) :: block_significant
      integer, dimension(n_blocks + 1), intent(in) :: block_firsts
!
      integer, dimension(dim_, columns), intent(in) :: array
      integer, dimension(dim_reduced, columns), intent(out) :: array_reduced
!
      integer :: block_, current_pos, first, last, size_, I
!
!$omp parallel do schedule(static) private(I, current_pos, block_, first, last, size_)
      do I = 1, columns
!
         current_pos = 1
!
         do block_ = 1, n_blocks
!
            if (block_significant(block_)) then
!
               first = block_firsts(block_)
               last  = block_firsts(block_ + 1) - 1
               size_ = last - first + 1
!
               array_reduced(current_pos : (current_pos + size_ - 1), I) = array(first : last, I)
!
               current_pos = current_pos + size_
!
            endif
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine reduce_array_int
!
!
   subroutine cholesky_decomposition_limited_diagonal(matrix, cholesky_vectors, dim_, &
                                                     n_vectors, threshold, n_included_diagonals, &
                                                     included_diagonals, n_vectors_requested)
!!
!!    Cholesky decomposition limited diagonal,
!!    Written by Sarai Dery Folkestad, June 2017.
!!
!!    Cholesky decomposition with pivots selected from a subset of the diagonals.
!!
!!    Routine is used for decomposition of density to construct active orbitals.
!!
!!    The number of pivots may specified through the optional
!!    argument n_vectors_requested  
!!
!!    On exit matrix_xy = matrix_xy - sum_J L_xJ*LyJ
!!
      implicit none
!
      integer, intent(in) :: dim_, n_included_diagonals
      integer, intent(out) :: n_vectors
!
      real(dp), intent(in) :: threshold
!
      integer, intent(in), optional :: n_vectors_requested
!
      real(dp), dimension(dim_, dim_), intent(inout) :: matrix
      real(dp), dimension(dim_, n_included_diagonals), intent(out) :: cholesky_vectors
!
      integer, dimension(n_included_diagonals), intent(in) :: included_diagonals
!
      integer, dimension(:), allocatable :: used_diag
!
      integer :: i, j, index_max, n_max_pivots
      real(dp) :: max_diagonal
!
      real(dp), dimension(:), allocatable :: diagonal
      real(dp), dimension(:), allocatable :: temp_cholesky_vector
!
      real(dp), parameter :: tolerance = 1.0d-10
!
      n_max_pivots = n_included_diagonals
!
      if (present(n_vectors_requested)) n_max_pivots = n_vectors_requested
!
      call mem%alloc(diagonal, dim_)
!
      do i = 1, dim_
!
         diagonal(i) = matrix(i, i)
!
      enddo
!
      cholesky_vectors = zero
!
      call mem%alloc(used_diag, dim_)
      used_diag = 0
!
      do i = 1, n_max_pivots
!
         n_vectors = i
!
!        Find the maximum diagonal
!
         index_max = 0
         max_diagonal = 0.0d0
!
         do j = 1, n_included_diagonals
!
            if (abs(diagonal(included_diagonals(j))) .gt. abs(max_diagonal)) then
!
               max_diagonal = diagonal(included_diagonals(j))
               index_max    = included_diagonals(j)
!
            endif
!
         enddo
!
!        Check against threshold and whether diagonal is negative
!
         if (max_diagonal .lt. 0.0d0) then
            if (abs(max_diagonal) .gt. tolerance) then
!
               call output%error_msg('Found negative diagonal in cholesky decomposition.')
!
            endif
         endif
!
         if (abs(max_diagonal) .lt. threshold) then
!
            n_vectors = n_vectors - 1
!
            call mem%dealloc(used_diag, dim_)
            call mem%dealloc(diagonal, dim_)
!
!           On exit, cholesky vectors subtracted from matrix
!
            call dgemm('N', 'T',          &
                        dim_,             &
                        dim_,             &
                        n_vectors,        &
                        -one,             &
                        cholesky_vectors, &
                        dim_,             &
                        cholesky_vectors, &
                        dim_,             &
                        one,              &
                        matrix,           &
                        dim_)
!
            return
!
         else
!
            used_diag(n_vectors) = index_max
!
         endif
!
!        Cholesky vectors
!
         cholesky_vectors(:, n_vectors) = matrix(:, index_max)
!
         if (n_vectors .gt. 1) then
!
            call mem%alloc(temp_cholesky_vector, n_vectors - 1)
            temp_cholesky_vector(:) = cholesky_vectors(index_max, 1 : n_vectors - 1)
!
            call dgemm('N', 'T',                         &
                        dim_,                            &
                        1,                               &
                        n_vectors - 1,                   &
                        -one,                            &
                        cholesky_vectors,                &
                        dim_,                            &
                        temp_cholesky_vector,            &
                        1,                               &
                        one,                             &
                        cholesky_vectors(1, n_vectors),  &
                        dim_)
!
            call mem%dealloc(temp_cholesky_vector, n_vectors - 1)
!
         endif
!
         do j = 1, n_vectors - 1
!
            cholesky_vectors(used_diag(j), n_vectors) = zero
!
         enddo
!
         call dscal(dim_, one/sqrt(max_diagonal), cholesky_vectors(1, n_vectors), 1)
!
         do j = 1, dim_
!
            diagonal(j) = diagonal(j) - cholesky_vectors(j, n_vectors)**2
!
         enddo
!
         diagonal(index_max) = zero
!
      enddo
!
      call mem%dealloc(used_diag, dim_)
      call mem%dealloc(diagonal, dim_)
!
!     On exit, cholesky vectors subtracted from matrix
!
      call dgemm('N', 'T',          &
                  dim_,             &
                  dim_,             &
                  n_vectors,        &
                  -one,             &
                  cholesky_vectors, &
                  dim_,             &
                  cholesky_vectors, &
                  dim_,             &
                  one,              &
                  matrix,           &
                  dim_)
!
   end subroutine cholesky_decomposition_limited_diagonal
!
!
   subroutine full_cholesky_decomposition_system(matrix, cholesky_vectors, dim_, &
                                                 n_vectors, threshold, used_diag)
!!
!!    Cholesky decomposition
!!    Written by Sarai Dery Folkestad, June 2017.
!!
!!    Wrapper for Lapack decomposition routine dpstrf.
!!
!!    matrix:                          matrix M we want to decompose, P^T M P = L L^T 
!!    cholesky_vectors (intent: out):  L 
!!    dim_:                            M and L are (dim) x (dim) arrays
!!
!!    used_diag (intent: out):         Vector (of dimension dim) containing the information in P:
!!                                     P(used_diag(j), j) = one. Other elements of P are zero.
!!                      
!!    threshold:                       Threshold to use in the Cholesky decomposition.
!!    
!!    n_vectors:                       Number of Cholesky vectors
!! 
      implicit none
!
      integer, intent(in) :: dim_
      integer, intent(out) :: n_vectors
!
      real(dp), intent(in) :: threshold
!
      real(dp), dimension(dim_, dim_), intent(in)  :: matrix
      real(dp), dimension(dim_, dim_), intent(out) :: cholesky_vectors
!
      integer, dimension(dim_) :: used_diag
!
      real(dp), dimension(:), allocatable :: work  ! work array for LAPACK
!
      integer :: info
      integer :: I, J
!
      cholesky_vectors = matrix
!
      call mem%alloc(work, (2*dim_))
!
!     DPSTRF computes the Cholesky factorization with complete pivoting
!     of a real symmetric positive semidefinite matrix.
!
      call dpstrf('L',        &
            dim_,             &
            cholesky_vectors, &
            dim_,             &
            used_diag,        &
            n_vectors,        &
            threshold,        &
            work,             &
            info)
!
      call mem%dealloc(work, (2*dim_))
!
      do I = 1, dim_ ! Zero upper unreferenced triangle
         do J = 1, I - 1
!
            cholesky_vectors(J, I) = zero
!
         enddo
      enddo
!
      if (info .lt. 0) then
         call output%error_msg('Cholesky decomposition failed! Something wrong in call to dpstrf')
      end if
!
   end subroutine full_cholesky_decomposition_system
!
!
   subroutine full_cholesky_decomposition_effective(matrix, cholesky_vectors, dim_, &
                                                    n_vectors, threshold, used_diag)
!!
!!    Cholesky decomposition,
!!    Written by Sarai Dery Folkestad, June 2017.
!!
!!
      implicit none
!
      integer, intent(in) :: dim_
      integer, intent(out) :: n_vectors
!
      real(dp), intent(in) :: threshold
!
      real(dp), dimension(dim_, dim_), intent(inout) :: matrix
      real(dp), dimension(dim_, dim_), intent(out) :: cholesky_vectors
!
      integer, dimension(dim_), intent(out) :: used_diag
!
      integer :: i, j, index_max
      real(dp) :: max_diagonal, min_diagonal
!
      real(dp), dimension(:), allocatable :: diagonal
      real(dp), dimension(:), allocatable :: temp_cholesky_vector
!
      real(dp), parameter :: tolerance = 1.0d-10
!
      used_diag = 0
!
!     Looping over the number of cholesky vectors
!
      call mem%alloc(diagonal, dim_)
!
      do i = 1, dim_
!
         diagonal(i) = matrix(i, i)
!
      enddo
!
      do i = 1, dim_
!
         n_vectors = i
!
!        Find the maximum diagonal
!
         index_max = 0
         max_diagonal = 0.0d0
!
         do j = 1, dim_
!
            if (abs(diagonal(j)) .gt. abs(max_diagonal)) then
!
               max_diagonal = diagonal(j)
               index_max    = j
!
            endif
!
         enddo
!
!        Check against threshold and whether diagonal is negative
!
         if (max_diagonal .lt. 0.0d0) then
            if (abs(max_diagonal) .gt. tolerance) then
!
               call output%error_msg('Found negative diagonal in cholesky decomposition.')
!
            endif
         endif
!
         if (abs(max_diagonal) .lt. threshold) then
!
            n_vectors = n_vectors - 1
!
            min_diagonal = 1.0D10
!
            do j = 1, dim_
!
            call dgemm('N', 'T',          &
                        dim_,             &
                        dim_,             &
                        n_vectors,        &
                        -one,             &
                        cholesky_vectors, &
                        dim_,             &
                        cholesky_vectors, &
                        dim_,             &
                        one,              &
                        matrix,           &
                        dim_)
               if (diagonal(j) .lt. min_diagonal) min_diagonal = diagonal(j)

!
            enddo
!
            call output%printf('n', 'The smallest diagonal after decomposition &
                               &is: (e12.4)', reals=[min_diagonal], fs='(/t6,a)')
            call mem%dealloc(diagonal, dim_)
!
            return
!
         else
!
            used_diag(n_vectors) = index_max
!
         endif
!
!        Cholesky vectors
!
         cholesky_vectors(:, n_vectors) = matrix(:, index_max)
!
         if (n_vectors .gt. 1) then
!
            call mem%alloc(temp_cholesky_vector, n_vectors - 1)
            temp_cholesky_vector(:) = cholesky_vectors(index_max, 1 : n_vectors - 1)
!
            call dgemm('N', 'T',                         &
                        dim_,                            &
                        1,                               &
                        n_vectors - 1,                   &
                        -one,                            &
                        cholesky_vectors,                &
                        dim_,                            &
                        temp_cholesky_vector,            &
                        1,                               &
                        one,                             &
                        cholesky_vectors(1, n_vectors),  &
                        dim_)
!
            call mem%dealloc(temp_cholesky_vector, n_vectors - 1)
!
         endif
!
         do j = 1, n_vectors - 1
!
            cholesky_vectors(used_diag(j), n_vectors) = zero
!
         enddo
!
         call dscal(dim_, one/sqrt(max_diagonal), cholesky_vectors(1, n_vectors), 1)
!
         do j = 1, dim_
!
            diagonal(j) = diagonal(j) - cholesky_vectors(j, n_vectors)**2
!
         enddo
!
         diagonal(index_max) = zero
!
         do j = 1, dim_
!
            matrix(j,index_max) = 0.0D0
            matrix(index_max,j) = 0.0D0
!
         enddo
!
      enddo
!
      min_diagonal = 1.0D10
!
      do j = 1, dim_
!
            if (diagonal(j) .lt. min_diagonal) min_diagonal = diagonal(j)

!
      enddo
!
      call output%printf('n', 'The smallest diagonal after decomposition is: (e12.4)', &
                         reals=[min_diagonal], fs='(/t6,a)')
!
      call mem%dealloc(diagonal, dim_)
!
   end subroutine full_cholesky_decomposition_effective
!
!
   subroutine invert(Ainv, A, n)
!!
!!    Invert matrix A
!!    Written by Sarai D. Folkestad, 2018
!!
!!    Inverts n x n - matrix A and places the result in Ainv.
!!
      implicit none
!
      integer, intent(in) :: n
!
      real(dp), dimension(n,n), intent(in) :: A
      real(dp), dimension(n,n), intent(inout) :: Ainv
!
!     Store A in Ainv to prevent it from being overwritten by LAPACK
      call dcopy(n**2, A, 1, Ainv, 1)
!
      call invert_in_place(Ainv, n)
!
   end subroutine invert
!
!
   subroutine invert_in_place(A, n)
!!
!!    Invert in place
!!    Written by Sarai D. Folkestad, 2018
!!
!!    Inverts n x n - matrix A and places the result in A.
!!
      implicit none
!
      integer, intent(in) :: n
!
      real(dp), dimension(n,n), intent(inout) :: A
!
      real(dp), dimension(n) :: work  ! work array for LAPACK
      integer, dimension(n) :: ipiv   ! pivot indices
      integer :: info
!
!     DGETRF computes an LU factorization of a general M-by-N matrix A
!     using partial pivoting with row interchanges.
!
      call dgetrf(n, n, A, n, ipiv, info)
!
      if (info /= 0) &
         call output%error_msg('Matrix is numerically singular! dgetrf error integer: (i0)', ints=[info])
!
!     DGETRI computes the inverse of a matrix using the LU factorization
!     computed by DGETRF.
!
      call dgetri(n, A, n, ipiv, work, n, info)
!
      if (info /= 0) &
         call output%error_msg('Matrix inversion failed! dgetri error integer: (i0)', ints=[info])
!
   end subroutine invert_in_place
!
!
   subroutine invert_lower_triangular(Ainv, A, n)
!!
!!    Invert lower triagonal
!!    Written by Sarai D. Folkestad, 2018
!!
!!    Inverts lower triangular n x n - matrix A and places the result in Ainv.
!!
      implicit none
!
      integer, intent(in) :: n
!
      real(dp), dimension(n, n), intent(in) :: A
      real(dp), dimension(n, n), intent(out) :: Ainv
!
      integer :: info
!
!     Store A in Ainv to prevent it from being overwritten by LAPACK
!
      Ainv = A
!
!     dtrtri computes the inverse of a real upper or lower triangular
!     matrix A.
!
      call dtrtri('l','n', n, Ainv, n, info)
!
      if (info /= 0) then 
         call output%error_msg('Matrix inversion failed.' // &
                               ' "Dtrtri" finished with info: (i0)', ints=[info])
      end if
!
   end subroutine invert_lower_triangular
!
!
   pure real(dp) function get_abs_max(X, n)
!!
!!    Get absolute maximum value of vector
!!    Written by Eirik F. Kjønstad, 2018
!!
!!    Returns the largest element, in absolute valute value, of X (length: n)
!!
      implicit none 
!
      integer, intent(in) :: n
!
      real(dp), dimension(n), intent(in) :: X
!
      integer :: i
!
      get_abs_max = zero
!
      do i = 1, n
!
         if (abs(X(i)) .gt. get_abs_max) get_abs_max = abs(X(i))
!
      enddo
!
   end function get_abs_max
!
!
   subroutine get_abs_max_w_index(X, n, abs_max, index_)
!!
!!    Get absolute maximum value and index of vector
!!    Written by Sarai D. Folkestad, Sep 2019
!!
!!    Computes the largest element, in absolute valute value, of X (length: n),
!!    and the index of X associated with that value. 
!!
!!    The routine is a copy of the function get_abs_max 
!!    written by Eirik F. Kjønstad. S.D.F added the
!!    index_ stuff.
!!
      implicit none 
!
      integer, intent(in) :: n
!
      real(dp), dimension(n), intent(in) :: X
!
      real(dp) :: abs_max
!
      integer :: index_
!
      integer :: i
!
      abs_max = zero
      index_ = 0
!
      do i = 1, n
!
         if (abs(X(i)) .gt. abs_max) then
            abs_max = abs(X(i))
            index_ = i
         endif
!
      enddo
!
   end subroutine get_abs_max_w_index
!
!
   subroutine sandwich_real(X, A, B, n, left)
!!
!!    Sandwich real
!!    Written by Eirik F. Kjønstad, 2018
!!
!!    Overwrites the n-by-n matrix X with a sandwich-ed X:
!!
!!       X <- A^T X B
!!
!!    All matrices are n x n.
!!
!!    left (optional): if true, transpose the left factor, A (standard); 
!!    if false, tranpose the right factor, B. 
!!
      implicit none
!
      integer, intent(in) :: n
!
      real(dp), dimension(n, n), intent(in) :: A
      real(dp), dimension(n, n), intent(in) :: B
!
      real(dp), dimension(n, n) :: X
!
      logical, optional, intent(in) :: left 
!
      real(dp), dimension(:, :), allocatable :: tmp
!
      call mem%alloc(tmp, n, n)
!
      if (present(left)) then 
!
         if (left) then ! Transpose the left factor, X <- A^T X B
!
            call dgemm('N', 'N', &
                        n,       &
                        n,       &
                        n,       &
                        one,     &
                        X,       &
                        n,       &
                        B,       &
                        n,       &
                        zero,    &
                        tmp,     & ! tmp = X B
                        n)
!
            call dgemm('T', 'N', &
                        n,       &
                        n,       &
                        n,       &
                        one,     &
                        A,       &
                        n,       &
                        tmp,     &
                        n,       &
                        zero,    &
                        X,       & ! X = A^T tmp = A^T X B
                        n)
!
         else ! Transpose the right factor, X <- A X B^T
!
            call dgemm('N', 'T', &
                        n,       &
                        n,       &
                        n,       &
                        one,     &
                        X,       &
                        n,       &
                        B,       &
                        n,       &
                        zero,    &
                        tmp,     & ! tmp = X B^T
                        n)
!
            call dgemm('N', 'N', &
                        n,       &
                        n,       &
                        n,       &
                        one,     &
                        A,       &
                        n,       &
                        tmp,     &
                        n,       &
                        zero,    &
                        X,       & ! X = A tmp = A X B^T
                        n)
!
         endif
!
      else ! Transpose left factor (standard)
!  
            call dgemm('N', 'N', &
                        n,       &
                        n,       &
                        n,       &
                        one,     &
                        X,       &
                        n,       &
                        B,       &
                        n,       &
                        zero,    &
                        tmp,     & ! tmp = X B
                        n)
!
            call dgemm('T', 'N', &
                        n,       &
                        n,       &
                        n,       &
                        one,     &
                        A,       &
                        n,       &
                        tmp,     &
                        n,       &
                        zero,    &
                        X,       & ! X = A^T tmp = A^T X B
                        n)
!
      endif
!
      call mem%dealloc(tmp, n, n)
!
   end subroutine sandwich_real
!
!
   subroutine sandwich_complex(X, A, B, n, left)
!!
!!    Sandwich complex
!!    Written by Andreas Skeidsvoll, Jan 2020
!!
!!    Overwrites the n-by-n matrix X with a sandwich-ed X:
!!
!!       X <- A^T X B
!!
!!    All matrices are n x n.
!!
!!    left (optional): if true, transpose the left factor, A (standard); 
!!    if false, tranpose the right factor, B. 
!!
!!    Based on sandwich_real by Eirik F. Kjønstad, 2018
!!
      implicit none
!
      integer, intent(in) :: n
!
      complex(dp), dimension(n, n), intent(in) :: A
      complex(dp), dimension(n, n), intent(in) :: B
!
      complex(dp), dimension(n, n) :: X
!
      logical, optional, intent(in) :: left 
!
      complex(dp), dimension(:, :), allocatable :: tmp
!
      call mem%alloc(tmp, n, n)
!
      if (present(left)) then 
!
         if (left) then ! Transpose the left factor, X <- A^T X B
!
            call zgemm('N', 'N',      &
                        n,            &
                        n,            &
                        n,            &
                        one_complex,  &
                        X,            &
                        n,            &
                        B,            &
                        n,            &
                        zero_complex, &
                        tmp,          & ! tmp = X B
                        n)
!
            call zgemm('T', 'N',      &
                        n,            &
                        n,            &
                        n,            &
                        one_complex,  &
                        A,            &
                        n,            &
                        tmp,          &
                        n,            &
                        zero_complex, &
                        X,            & ! X = A^T tmp = A^T X B
                        n)
!
         else ! Transpose the right factor, X <- A X B^T
!
            call zgemm('N', 'T',      &
                        n,            &
                        n,            &
                        n,            &
                        one_complex,  &
                        X,            &
                        n,            &
                        B,            &
                        n,            &
                        zero_complex, &
                        tmp,          & ! tmp = X B^T
                        n)
!
            call zgemm('N', 'N',      &
                        n,            &
                        n,            &
                        n,            &
                        one_complex,  &
                        A,            &
                        n,            &
                        tmp,          &
                        n,            &
                        zero_complex, &
                        X,            & ! X = A tmp = A X B^T
                        n)
!
            endif
!
         else ! Transpose left factor (standard)
!  
            call zgemm('N', 'N',      &
                        n,            &
                        n,            &
                        n,            &
                        one_complex,  &
                        X,            &
                        n,            &
                        B,            &
                        n,            &
                        zero_complex, &
                        tmp,          & ! tmp = X B
                        n)
!
            call zgemm('T', 'N',      &
                        n,            &
                        n,            &
                        n,            &
                        one_complex,  &
                        A,            &
                        n,            &
                        tmp,          &
                        n,            &
                        zero_complex, &
                        X,            & ! X = A^T tmp = A^T X B
                        n)
!
         endif
!
      call mem%dealloc(tmp, n, n)
!
   end subroutine sandwich_complex
!
!
   subroutine symmetric_sandwich(Xr, X, A, m, n)
!!
!!    Symmetric sandwich
!!    Written by Eirik F. Kjønstad, 2018
!!
!!    Constructs the similary transformed matrix X_r
!!
!!       Xr = A^T X A
!!
      implicit none
!
      integer, intent(in) :: m, n
!
      real(dp), dimension(m, n), intent(in) :: A
!
      real(dp), dimension(m, m) :: X
      real(dp), dimension(n, n) :: Xr
!
      real(dp), dimension(:, :), allocatable :: tmp
!
      call mem%alloc(tmp, m, n)
!
      call dgemm('N', 'N', &
                  m,       &
                  n,       &
                  m,       &
                  one,     &
                  X,       &
                  m,       &
                  A,       &
                  m,       &
                  zero,    &
                  tmp,     & ! tmp = X A
                  m)
!
      call dgemm('T', 'N', &
                  n,       &
                  n,       &
                  m,       &
                  one,     &
                  A,       &
                  m,       &
                  tmp,     &
                  m,       &
                  zero,    &
                  Xr,      & ! X = A^T tmp = A^T X B
                  n)
!
      call mem%dealloc(tmp, m, n)
!
   end subroutine symmetric_sandwich
!
!
   subroutine symmetric_sandwich_right_transposition_real(Xr, X, A, m, n)
!!
!!    Symmetric sandwich right transposition
!!    Written by Eirik F. Kjønstad, 2018
!!
!!    Constructs the similarity transformation
!!    wrt the matrix A^T
!!
!!       Xr = A X A^T
!!
      implicit none
!
      integer, intent(in) :: m, n
!
      real(dp), dimension(m, n), intent(in)  :: A
      real(dp), dimension(n, n), intent(in)  :: X
      real(dp), dimension(m, m), intent(out) :: Xr
!
      real(dp), dimension(:, :), allocatable :: tmp
!
      call mem%alloc(tmp, n, m)
!
      call dgemm('N', 'T', &
                  n,       &
                  m,       &
                  n,       &
                  one,     &
                  X,       &
                  n,       &
                  A,       &
                  m,       &
                  zero,    &
                  tmp,     & ! tmp = X A^T
                  n)
!
      call dgemm('N', 'N', &
                  m,       &
                  m,       &
                  n,       &
                  one,     &
                  A,       &
                  m,       &
                  tmp,     &
                  n,       &
                  zero,    &
                  Xr,      & ! X = A tmp = A X A^T 
                  m)
!
      call mem%dealloc(tmp, n, m)
!
   end subroutine symmetric_sandwich_right_transposition_real
!
!
   subroutine symmetric_sandwich_right_transposition_complex_between_real(Xr, X, A, m, n)
!!
!!    Symmetric sandwich right transposition
!!    Written by Eirik F. Kjønstad, 2018
!!
!!    Constructs the similarity transformation
!!    wrt the matrix A^T
!!
!!       Xr = A X A^T
!!
!!    Modified by Andreas S. Skeidsvoll, Jul 2020
!!    To be used with real A and complex X matrices.
!!
      implicit none
!
      integer, intent(in) :: m, n
!
      real(dp), dimension(m, n), intent(in)  :: A
      complex(dp), dimension(n, n), intent(in)  :: X
      complex(dp), dimension(m, m), intent(out) :: Xr
!
      complex(dp), dimension(:, :), allocatable :: tmp
!
      call mem%alloc(tmp, n, m)
!
      call zgemm('N', 'T',            &
                  n,                  &
                  m,                  &
                  n,                  &
                  one_complex,        &
                  X,                  &
                  n,                  &
                  cmplx(A, zero, dp), &
                  m,                  &
                  zero_complex,       &
                  tmp,                & ! tmp = X A^T
                  n)
!
      call zgemm('N', 'N',            &
                  m,                  &
                  m,                  &
                  n,                  &
                  one_complex,        &
                  cmplx(A, zero, dp), &
                  m,                  &
                  tmp,                &
                  n,                  &
                  zero_complex,       &
                  Xr,                 & ! X = A tmp = A X A^T 
                  m)
!
      call mem%dealloc(tmp, n, m)
!
   end subroutine symmetric_sandwich_right_transposition_complex_between_real
!
!
   subroutine symmetric_sandwich_right_transposition_complex(Xr, X, A, m, n)
!!
!!    Symmetric sandwich right transposition
!!    Written by Eirik F. Kjønstad, 2018
!!
!!    Constructs the similarity transformation
!!    wrt the matrix A^T
!!
!!       Xr = A X A^T
!!
      implicit none
!
      integer, intent(in) :: m, n
!
      complex(dp), dimension(m, n), intent(in)  :: A
      complex(dp), dimension(n, n), intent(in)  :: X
      complex(dp), dimension(m, m), intent(out) :: Xr
!
      complex(dp), dimension(:, :), allocatable :: tmp
!
      call mem%alloc(tmp, n, m)
!
      call zgemm('N', 'T',      &
                  n,            &
                  m,            &
                  n,            &
                  one_complex,  &
                  X,            &
                  n,            &
                  A,            &
                  m,            &
                  zero_complex, &
                  tmp,          & ! tmp = X A^T
                  n)
!
      call zgemm('N', 'N',      &
                  m,            &
                  m,            &
                  n,            &
                  one_complex,  &
                  A,            &
                  m,            &
                  tmp,          &
                  n,            &
                  zero_complex, &
                  Xr,           & ! X = A tmp = A X A^T 
                  m)
!
      call mem%dealloc(tmp, n, m)
!
   end subroutine symmetric_sandwich_right_transposition_complex
!
!
   subroutine symmetric_sandwich_right_transposition_replace(X, A, m)
!!
!!    Symmetric sandwich right transposition replace
!!    Written by Eirik F. Kjønstad, 2018
!!
!!    Overwrites the n-by-n matrix X with a sandwich-ed X:
!!
!!       X <- A X A^T
!!
      implicit none
!
      integer, intent(in) :: m
!
      real(dp), dimension(m, m), intent(in) :: A
!
      real(dp), dimension(m, m), intent(inout) :: X
!
      real(dp), dimension(:, :), allocatable :: tmp
!
      call mem%alloc(tmp, m, m)
!
      call dgemm('N', 'T', &
                  m,       &
                  m,       &
                  m,       &
                  one,     &
                  X,       &
                  m,       &
                  A,       &
                  m,       &
                  zero,    &
                  tmp,     & ! tmp = X A^T
                  m)
!
      call dgemm('N', 'N', &
                  m,       &
                  m,       &
                  m,       &
                  one,     &
                  A,       &
                  m,       &
                  tmp,     &
                  m,       &
                  zero,    &
                  X,       & ! X = A tmp = A X A^T 
                  m)
!
      call mem%dealloc(tmp, m, m)
!
   end subroutine symmetric_sandwich_right_transposition_replace
!
!
   real(dp) function get_l2_norm(X, n)
!!
!!    Get L^2 norm 
!!    Written by Eirik F. Kjønstad, Aug 2018 
!!
!!    Returns the L^2 norm of the n-dimensional vector X: 
!!
!!       sqrt( sum_i=1^n X_i^2 )
!!
      implicit none 
!
      integer, intent(in) :: n
!
      real(dp), dimension(n), intent(in) :: X 
!
      real(dp) :: ddot 
!
      get_l2_norm = sqrt(ddot(n, X, 1, X, 1))
!
   end function get_l2_norm
!
!
   subroutine transpose_(A, A_trans, dim_)
!!
!!    Transpose 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
!!    Saves the transpose of A in A_trans. The matrices are assumed to 
!!    be dim_ x dim_.
!!
      implicit none
!
      integer :: dim_, n, m
!
      real(dp), dimension(dim_, dim_), intent(in) :: A
      real(dp), dimension(dim_, dim_), intent(out) :: A_trans
!
!$omp parallel do private(m, n) shared(A, A_trans)
      do m = 1, dim_ 
         do n = 1, dim_ 
!
            A_trans(n, m) = A(m, n)
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine transpose_
!
!
   subroutine copy_and_scale(alpha, X, Y, n)
!!
!!    Copy and scale
!!    Written by Sarai D. Folkestad, May 2019
!!
!!    Sets Y as:
!!
!!       Y = alpha*X
!!
!!    X and Y are vectors of length n, and alpha is a real number.
!!
      implicit none
!
      integer, intent(in) :: n
!
      real(dp), dimension(n), intent(out) :: Y
      real(dp), dimension(n), intent(in) :: X
!
      real(dp), intent(in) :: alpha
!
      integer :: i
!
!$omp parallel do private(i)
      do i = 1, n
!
         Y(i) = alpha*X(i)
!
      enddo
!$omp end parallel do 
!
   end subroutine copy_and_scale
!
!
   subroutine copy_and_scale_complex(alpha, X, Y, n)
!!
!!    Copy and scale
!!    Written by Sarai D. Folkestad, May 2019
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Y = alpha*X
!!
      implicit none
!
      integer, intent(in) :: n
!
      complex(dp), dimension(n), intent(out) :: Y
      complex(dp), dimension(n), intent(in) :: X
!
      complex(dp), intent(in) :: alpha
!
      integer :: i
!
!$omp parallel do private(i)
      do i = 1, n
!
         Y(i) = alpha*X(i)
!
      enddo
!$omp end parallel do 
!
   end subroutine copy_and_scale_complex
!
!
   subroutine get_n_lowest(n, size, vec, sorted_short_vec, index_list)
!!
!!    Get n lowest elements
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!    Finds the n lowest values of "vec" (of length "size"), sorts them, and returns them
!!    in "sorted_short_vec", together with an index list (of length "n") refering to the 
!!    indices of the lowest elements in the original vector "vec".
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
!!    Finds the n highest values of "vec" (of length "size"), sorts them, and returns them
!!    in "sorted_short_vec", together with an index list (of length "n") refering to the 
!!    indices of the lowest elements in the original vector "vec".
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
   recursive subroutine quicksort_recursive(vec, first, last)
!!
!!    Recursive implementation of quicksort (descending order)
!!
!!    Adapted from quicksort.f by "t-nissie" (https://gist.github.com/t-nissie/479f0f16966925fa29ea)
!!    licensed under GPLv3.
!!
!!    Modified by Sarai D. Folkestad, 6. Aug. 2018
!!
!!    - Type real*8 changed to real(dp)       
!!    - Variable x renamed to pivot
!!    - Variable t renamed to temp
!!
!     quicksort.f -*-f90-*-
!     Author: t-nissie
!     License: GPLv3
!     Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
!
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
!
   recursive subroutine quicksort_with_index_recursive(vec, index_list, first, last)
!!
!!    Recursive implementation of quicksort with index list (descending order)
!!
!!    Adapted from quicksort.f by "t-nissie" (https://gist.github.com/t-nissie/479f0f16966925fa29ea)
!!    licensed under GPLv3.
!!
!!    "index_list" : integer array which stores the original index of ordered elements
!!                   e.g., in index_list(1) is the original index of the largest element
!!
!!
!!    Modified by Sarai D. Folkestad, 11. Nov. 2018
!!
!!    - Type real*8 changed to real(dp)       
!!    - Variable x renamed to pivot
!!    - Variable t renamed to temp
!!    - Added index_list
!!
!     quicksort.f -*-f90-*-
!     Author: t-nissie
!     License: GPLv3
!     Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
!
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
!
   subroutine quicksort_descending(vec, dim_)
!!
!!    Wrapper for recursive quicksort routine
!!
      implicit none
!
      integer, intent(in) :: dim_
      real(dp), dimension(dim_), intent(inout) :: vec
!
      call quicksort_recursive(vec, 1, dim_)
!
   end subroutine quicksort_descending
!
!
   subroutine quicksort_with_index_descending(vec, index_list, dim_)
!!
!!    Wrapper for recursive quicksort with index list
!!
      implicit none
!
      integer, intent(in) :: dim_
      real(dp), dimension(dim_), intent(inout) :: vec
      integer, dimension(dim_), intent(out) :: index_list
!
      integer :: i
!
      do i = 1, dim_
         index_list(i) = i
      enddo
!
      call quicksort_with_index_recursive(vec, index_list, 1, dim_)
!
   end subroutine quicksort_with_index_descending
!
!
   subroutine quicksort_with_index_ascending(vec, index_list, dim_)
!!
!!    Wrapper for recursive quicksort with index list ascending order
!!
      implicit none
!
      integer, intent(in) :: dim_
      real(dp), dimension(dim_), intent(inout) :: vec
      integer, dimension(dim_), intent(out) :: index_list
!
      call dscal(dim_, -one, vec, 1)
!
      call quicksort_with_index_descending(vec, index_list, dim_)
!
      call dscal(dim_, -one, vec, 1)
!
   end subroutine quicksort_with_index_ascending
!
!
   recursive subroutine quicksort_recursive_int(vec, first, last)
!!
!!    Recursive implementation of quicksort (descending order)
!!
!!    Adapted from quicksort.f by "t-nissie" (https://gist.github.com/t-nissie/479f0f16966925fa29ea)
!!    licensed under GPLv3.
!!
!!    Modified by Sarai D. Folkestad, 6. Aug. 2018
!!   
!!    - Variable x renamed to pivot
!!    - Variable t renamed to temp
!!    - Added index_list
!!
!!    Modified by Andreas Skeidsvoll, 26. Feb. 2019
!!
!!    - Type real*8 changed to integer    
!!
!     quicksort.f -*-f90-*-
!     Author: t-nissie
!     License: GPLv3
!     Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
!
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
!
   recursive subroutine quicksort_with_index_recursive_int(vec, index_list, first, last)
!!
!!    Recursive implementation of quicksort with index list (descending order)
!!
!!    Adapted from quicksort.f by "t-nissie" (https://gist.github.com/t-nissie/479f0f16966925fa29ea)
!!    licensed under GPLv3.
!!
!!    "index_list" : integer array which stores the original index of ordered elements
!!                   e.g., in index_list(1) is the original index of the largest element
!!
!!
!!    Modified by Sarai D. Folkestad, 11. Nov. 2018
!!   
!!    - Variable x renamed to pivot
!!    - Variable t renamed to temp
!!    - Added index_list
!!
!!    Modified by Andreas Skeidsvoll, 26. Feb. 2019
!!
!!    - Type real*8 changed to integer    
!!
!     quicksort.f -*-f90-*-
!     Author: t-nissie
!     License: GPLv3
!     Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
!
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
!
   subroutine quicksort_descending_int(vec, dim_)
!!
!!    Wrapper for recursive quicksort routine
!!
      implicit none
!
      integer, intent(in) :: dim_
      integer, dimension(dim_), intent(inout) :: vec
!
      call quicksort_recursive_int(vec, 1, dim_)
!
   end subroutine quicksort_descending_int
!
!
   subroutine quicksort_with_index_descending_int(vec, index_list, dim_)
!!
!!    Wrapper for recursive quicksort with index list
!!
      implicit none
!
      integer, intent(in) :: dim_
      integer, dimension(dim_), intent(inout) :: vec
      integer, dimension(dim_), intent(inout) :: index_list
!
      integer :: i
!
      do i = 1, dim_
         index_list(i) = i
      enddo
!
      call quicksort_with_index_recursive_int(vec, index_list, 1, dim_)
!
   end subroutine quicksort_with_index_descending_int
!
!
   subroutine quicksort_with_index_ascending_int(vec, index_list, dim_)
!!
!!    Wrapper for recursive quicksort with index list ascending order
!!
      implicit none
!
      integer, intent(in) :: dim_
      integer, dimension(dim_), intent(inout) :: vec
      integer, dimension(dim_), intent(inout) :: index_list
!
      vec = -1*vec
!
      call quicksort_with_index_descending_int(vec, index_list, dim_)
!
      vec = -1*vec
!
   end subroutine quicksort_with_index_ascending_int
!
!
   subroutine generalized_diagonalization(dim_, A, M, V, e)
!!
!!    Generalized diagonalization
!!    Written by Sarai D. Folkestad, Oct 2019
!!
!!    Solves the generalized eigenvalue
!!    problem
!!
!!       A V = M X e
!!
!!    'dim_' : dimension of matrix A
!!
!!    'A' : Matrix to be diagonalized
!!
!!    'M' : Metric
!!
!!    'V' : Eigenvectors
!!
!!    'e' : Eigenvalues (optional)
!!
!!    Wrapper for dsygv (LAPACK)
!!
      implicit none
!
      integer, intent(in) :: dim_
!
      real(dp), dimension(dim_, dim_), intent(in) :: A
      real(dp), dimension(dim_, dim_), intent(in) :: M
!
      real(dp), dimension(dim_, dim_), intent(out) :: V
!
      real(dp), dimension(dim_), optional, intent(out) :: e
!
      real(dp), dimension(:), allocatable :: work, eigenvalues
!
      integer :: info, lwork
!
!     Work space inquiry
!
      call dcopy(dim_**2, A, 1, V, 1)
!
      call mem%alloc(work, 1)
      lwork = -1
!
      call mem%alloc(eigenvalues, dim_)
!
      call dsygv(1, 'V', 'L',       &
                  dim_,             &
                  V,                &
                  dim_,             &
                  M,                &
                  dim_,             &
                  eigenvalues,      &
                  work,             &
                  lwork,            &
                  info)
!
      lwork = ceiling(work(1))
      call mem%dealloc(work, 1)
!
!     Diagonalization
!
      call mem%alloc(work, lwork)
!
      call dsygv(1, 'V', 'L',       &
                  dim_,             &
                  V,                &
                  dim_,             &
                  M,                &
                  dim_,             &
                  eigenvalues,      &
                  work,             &
                  lwork,            &
                  info)
!
      if (info .ne. 0) call output%error_msg('Could not solve generalized eigenvalue problem')
!
      call mem%dealloc(work, lwork)
!
      if (present(e)) then
!
         call dcopy(dim_, eigenvalues, 1, e, 1)
!
      endif
!
      call mem%dealloc(eigenvalues, dim_)
!
   end subroutine generalized_diagonalization
!
!
   function count_n_true(dim_, logical_list) result(n_true)
!!
!!    Count n true
!!    Written by Sarai D. Folkestad, Oct 2018
!!
!!    Counts the number of .true.
!!    in the list of logicals ('logical_list')
!!    with dimension 'dim_'
!!
      implicit none
!
      integer, intent(in) :: dim_
      logical, dimension(dim_), intent(in) :: logical_list
!
      integer :: n_true
!
      integer :: I
!
      n_true = 0
!
      do I = 1, dim_
!
         if (logical_list(I)) n_true = n_true + 1
!
      enddo
!
   end function count_n_true
!
!
   subroutine extract_columns_of_matrix(n_extract, n_rows, n_columns, M, M_extracted, extract)
!!
!!    Extract columns of matrix
!!    Written by Ida-Marie Høyvik and Sarai D. Folkestad, Oct 2019
!!
!!    Extracts columns of a matrix, based on the logical list
!!    'extract' where extract(I)=.true. defines the instruction 
!!    to extract. The number of elements to extract
!!    'n_extract' must already be determined and passed to the routine.
!!
      implicit none
!
      integer, intent(in) :: n_extract, n_rows, n_columns
!
      real(dp), dimension(n_rows, n_columns), intent(in) :: M
      real(dp), dimension(n_rows, n_extract), intent(out) :: M_extracted
!
      logical, dimension(n_columns) :: extract
!
      integer :: counter, I
!
      if (n_extract .gt. n_columns) &
         call output%error_msg('Can not extract more columns than the original matrix contains: ' &
                               //'mismatch beween n_extract and number of true in the extract list')
!
      counter = 0
!
      do I = 1, n_columns
!  
         if (extract(I)) then
!
            counter = counter + 1
!
            if (counter .gt. n_extract) call output%error_msg('Tried to extract non-existent column from matrix')
!
            call dcopy(n_rows, M(1,I), 1, M_extracted(1,counter), 1)
!     
         endif
!
      enddo
!
   end subroutine extract_columns_of_matrix
!
!
   function are_vectors_parallel(x, y, dim_, threshold) result(parallel)
!!
!!    Are vectors parallel
!!    Written by Alexander C. Paul and Rolf H. Myhre, Oct 2019
!!
!!    Checks if 2 vectors are parallel within a given threshold.
!!    For 2 parallel vectors x and y = a*x
!!          sqrt(x*x * y*y) = |x*y|
!!
      implicit none
!
      logical :: parallel
!
      integer, intent(in) :: dim_
!
      real(dp), dimension(dim_), intent(in) :: x, y
!
      real(dp), intent(in) :: threshold
!
      real(dp) :: ddot, x_dot_x, y_dot_y, x_dot_y
!
      x_dot_x = ddot(dim_, x, 1, x, 1)
!
      y_dot_y = ddot(dim_, y, 1, y, 1)
!
      x_dot_y = ddot(dim_, x, 1, y, 1)
!
      if(abs(sqrt(x_dot_x*y_dot_y) - abs(x_dot_y)) .gt. threshold) then
!
         parallel = .false.
!
      else
!
         parallel = .true.
!
      end if
!
   end function are_vectors_parallel
!
!
   subroutine check_for_parallel_vectors(vectors, dim_, n_vectors,    &
                                         n_non_parallel_vectors,     &
                                         threshold, parallel)
!!
!!    Check for parallel vectors
!!    written by Alexander Paul, Oct 2019
!!
!!    Checks for parallel vectors in a given array of vectors
!!    returns the number of not parallel vectors and which vectors parallel
!!
!!    vectors:  array of vectors to be checked
!!    parallel: parallel(p) will be true on out if vector p is parallel to a
!!              previous vector for example, if there are three vectors and they
!!              are all parallel, parallel will be [.false., .true., .true.] on output
!!
      implicit none
!
      integer, intent(in)  :: n_vectors
      integer, intent(out) :: n_non_parallel_vectors
!
      integer, intent(in) :: dim_
!
      real(dp), dimension(dim_, n_vectors), intent(in) :: vectors
!
      real(dp), intent(in) :: threshold
!
      logical, dimension(n_vectors), intent(out), optional :: parallel
!
      logical, dimension(n_vectors) :: parallel_temp
!
      integer :: p, q
!
      n_non_parallel_vectors = n_vectors
!
      do p = 1, n_vectors
!
         parallel_temp(p) = .false.
!
         do q = 1, p-1
!
!        Skip if the vector q was already parallel to another p
!
            if (parallel_temp(q)) cycle
!
               if(are_vectors_parallel(vectors(:,p), vectors(:,q), dim_, threshold)) then
!
                  parallel_temp(p) = .true.
!
                  n_non_parallel_vectors = n_non_parallel_vectors - 1
!
               end if
!
         end do
      end do
!
      if (present(parallel)) then
         parallel = parallel_temp
      end if
!
   end subroutine check_for_parallel_vectors
!
!
   subroutine entrywise_product_in_place(dim_, X, Y)
!!
!!    entrywise product in place
!!    Written by Alexander C. Paul, Sep 2019
!!
!!    Also called Hadamard or Schur product
!!
!!    Scales each element of a vector 
!!    by the corresponding element of another vector
!!             X(p) = X(p) * Y(p)
!!
!!    NB: Vector X contains product on exit
!!
!!    Used eg. for applying a projector
!!
      implicit none
!
      integer, intent(in) :: dim_
!
      real(dp), dimension(dim_), intent(inout)  :: X
      real(dp), dimension(dim_), intent(in)     :: Y
!
      integer :: p
!
!$omp parallel do private(p)
      do p = 1, dim_
!
         X(p) = X(p)*Y(p)
!
      enddo  
!$omp end parallel do
!
   end subroutine entrywise_product_in_place
!
!
   subroutine entrywise_product_in_place_2dim(dim_, X, Y)
!!
!!    entrywise product in place of 2 dimensional vectors
!!    Written by Alexander C. Paul, Sep 2019
!!
!!    Also called Hadamard or Schur product
!!
!!    Scales each element of a 2 dimensional vector 
!!    by the corresponding element of another vector
!!             X(p) = X(p) * Y(p)
!!
!!    NB: Vector X contains product on exit
!!
!!    dim_ : full dimension of the vector, dim_1*dim_2
!!
!!    Used eg. for applying a projector
!!
      implicit none
!
      integer, intent(in) :: dim_
!
      real(dp), dimension(:,:), intent(inout)  :: X
      real(dp), dimension(:,:), intent(in)     :: Y
!
      call entrywise_product_in_place(dim_, X, Y)
!
   end subroutine entrywise_product_in_place_2dim
!
!
   subroutine entrywise_product_in_place_3dim(dim_, X, Y)
!!
!!    entrywise product in place 3-dimensional vectors
!!    Written by Alexander C. Paul, Sep 2019
!!
!!    Also called Hadamard or Schur product
!!
!!    Scales each element of a 3 dimensional vector 
!!    by the corresponding element of another vector
!!             X(p) = X(p) * Y(p)
!!
!!    NB: Vector X contains product on exit
!!
!!    dim_ : full dimension of the vector, dim_1*dim_2*dim_3
!!
!!    Used eg. for applying a projector
!!
      implicit none
!
      integer, intent(in) :: dim_
!
      real(dp), dimension(:,:,:), intent(inout)  :: X
      real(dp), dimension(:,:,:), intent(in)     :: Y
!
      call entrywise_product_in_place(dim_, X, Y)
!
   end subroutine entrywise_product_in_place_3dim
!
!
   subroutine entrywise_product_in_place_4dim(dim_, X, Y)
!!
!!    entrywise product in place 4-dimensional vectors
!!    Written by Alexander C. Paul, Sep 2019
!!
!!    Also called Hadamard or Schur product
!!
!!    Scales each element of a 4 dimensional vector 
!!    by the corresponding element of another vector
!!             X(p) = X(p) * Y(p)
!!
!!    NB: Vector X contains product on exit
!!
!!    dim_ : full dimension of the vector, dim_1*dim_2*dim_3*dim_4
!!
!!    Used eg. for applying a projector
!!
      implicit none
!
      integer, intent(in) :: dim_
!
      real(dp), dimension(:,:,:,:), intent(inout)  :: X
      real(dp), dimension(:,:,:,:), intent(in)     :: Y
!
      call entrywise_product_in_place(dim_, X, Y)
!
   end subroutine entrywise_product_in_place_4dim
!
!
   subroutine entrywise_product_(dim_, X, Y, Z)
!!
!!    entrywise product
!!    Written by Alexander C. Paul, Sep 2019
!!
!!    Also called Hadamard or Schur product
!!
!!    Scales each element of a vector 
!!    by the corresponding element of another vector
!!             Z(p) = X(p) * Y(p)
!!
!!    Result will be returned as the Z-vector
!!
      implicit none
!
      integer, intent(in) :: dim_
!
      real(dp), dimension(dim_), intent(in)  :: X
      real(dp), dimension(dim_), intent(in)  :: Y
      real(dp), dimension(dim_), intent(out) :: Z
!
      integer :: p
!
!$omp parallel do private(p)
      do p = 1, dim_
!
         Z(p) = X(p)*Y(p)
!
      enddo  
!$omp end parallel do
!
   end subroutine entrywise_product_
!
!
   subroutine scale_real_diagonal_by_real(alpha, X, dim_)
!!
!!    Scale diagonal of real array by real  
!!    Written by Anders Hutcheson, Oct 2019 
!!
!!    X is of dimension dim_^2 and we scale every 
!!    dim_+1 element(diagonal) by alpha,
!!    because we only scale the diagonal elements 
!!    the first argument is dim_    
!!

      implicit none
!
      real(dp), intent(in) :: alpha
!
      integer, intent(in) :: dim_
!
      real(dp), dimension(dim_, dim_), intent(inout) :: X 
!
      call dscal(dim_, alpha, X, dim_+1)
!
   end subroutine scale_real_diagonal_by_real
!
!
   subroutine scale_complex_diagonal_by_real(alpha, X, dim_)
!!
!!    Scale diagonal of complex array by real
!!    Written by Anders Hutcheson, Oct 2019 
!!
!!    X is of dimension dim_^2 and we scale every 
!!    dim_+1 element(diagonal) by alpha,
!!    because we only scale the diagonal elements 
!!    the first argument is dim_ 
!!
      implicit none
!
      real(dp), intent(in) :: alpha
!
      integer, intent(in) :: dim_
!
      complex(dp), dimension(dim_, dim_), intent(inout) :: X    
!
      call zdscal(dim_, alpha, X, dim_+1)
!
   end subroutine scale_complex_diagonal_by_real
!
!
   subroutine scale_complex_diagonal_by_complex(alpha, X, dim_)
!!
!!    Scale diagonal of complex array by complex
!!    Written by Anders Hutcheson, Oct 2019
!!
!!    X is of dimension dim_^2 and we scale every 
!!    dim_+1 element(diagonal) by alpha,
!!    because we only scale the diagonal elements 
!!    the first argument is dim_  
!! 
      implicit none
!
      complex(dp), intent(in) :: alpha
!
      integer, intent(in) :: dim_
!
      complex(dp), dimension(dim_, dim_), intent(inout) :: X    
!
      call zscal(dim_, alpha, X, dim_+1)
!
   end subroutine scale_complex_diagonal_by_complex
!
!
   subroutine scale_real_4_diagonal_by_real(alpha, X, dim_)
!!
!!    Scale diagonal of 4 dimensional real array by real 
!!    Written by Anders Hutcheson, Oct 2019 
!!
!!    The matrix has to have the following dimensions,
!!    X(dim_1,dim_2,dim_1,dim_2) and 
!!    dim_ is then given by dim_1*dim_2  
!!
      implicit none
!
      real(dp), intent(in) :: alpha
!
      real(dp), dimension(:,:,:,:), intent(inout) :: X
!
      integer, intent(in) :: dim_
!
      call scale_real_diagonal_by_real(alpha, X, dim_)
!
   end subroutine scale_real_4_diagonal_by_real
!
!
   subroutine scale_complex_4_diagonal_by_real(alpha, X, dim_)
!!
!!    Scale diagonal of 4 dimensional complex array by real  
!!    Written by Anders Hutcheson, Oct 2019 
!! 
!!    The matrix has to have the following dimensions,
!!    X(dim_1,dim_2,dim_1,dim_2) and 
!!    dim_ is then given by dim_1*dim_2  
!!
      implicit none
!
      real(dp), intent(in) :: alpha
!
      complex(dp), dimension(:,:,:,:), intent(inout) :: X
!
      integer, intent(in) :: dim_
!
      call scale_complex_diagonal_by_real(alpha, X, dim_)
!
   end subroutine scale_complex_4_diagonal_by_real
!
!
   subroutine scale_complex_4_diagonal_by_complex(alpha, X, dim_)
!!
!!    Scale diagonal of 4 dimensional complex array by complex  
!!    Written by Anders Hutcheson, Oct 2019 
!! 
!!    The matrix has to have the following dimensions,
!!    X(dim_1,dim_2,dim_1,dim_2) and 
!!    dim_ is then given by dim_1*dim_2   
!!
      implicit none
!
      complex(dp), intent(in) :: alpha
!
      complex(dp), dimension(:,:,:,:), intent(inout) :: X
!
      integer, intent(in) :: dim_
!
      call scale_complex_diagonal_by_complex(alpha, X, dim_)
!
   end subroutine scale_complex_4_diagonal_by_complex
!
!
   subroutine scale_real_4_diagonal_by_real_1324(alpha, X, dim_p, dim_q)
!!
!!    Scale diagonal of 4 dimensional real array sorted 1324 by real 
!!    Written by Anders Hutcheson, Oct 2019 
!!
!!    The matrix has to have the following dimensions,
!!    X(dim_p, dim_p, dim_q, dim_q) 
!!
      implicit none
!
      real(dp), intent(in) :: alpha
!
      integer, intent(in) :: dim_p, dim_q
!
      real(dp), dimension(dim_p,dim_p,dim_q,dim_q), intent(inout) :: X
!      
      integer :: p, q
!
!$omp parallel do schedule(static) private(p, q)
      do q = 1, dim_q
         do p = 1, dim_p
!
            X(p,p,q,q) = alpha*X(p,p,q,q)
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine scale_real_4_diagonal_by_real_1324
!
!
   subroutine scale_complex_4_diagonal_by_real_1324(alpha, X, dim_p, dim_q)
!!
!!    Scale diagonal of 4 dimensional complex array sorted 1324 by real 
!!    Written by Anders Hutcheson, Oct 2019 
!!
!!    The matrix has to have the following dimensions,
!!    X(dim_p, dim_p, dim_q, dim_q) 
!! 
      implicit none
!
      real(dp), intent(in) :: alpha
!
      integer, intent(in) :: dim_p, dim_q
!
      complex(dp), dimension(dim_p,dim_p,dim_q,dim_q), intent(inout) :: X
!      
      integer :: p, q
!
!$omp parallel do schedule(static) private(p, q)
      do q = 1, dim_q
         do p = 1, dim_p
!
            X(p,p,q,q) = alpha*X(p,p,q,q)
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine scale_complex_4_diagonal_by_real_1324
!
!
   subroutine scale_complex_4_diagonal_by_complex_1324(alpha, X, dim_p, dim_q)
!!
!!    Scale diagonal of 4 dimensional complex array sorted 1324 by complex 
!!    Written by Anders Hutcheson, Oct 2019 
!!
!!    The matrix has to have the following dimensions,
!!    X(dim_p, dim_p, dim_q, dim_q) 
!! 
      implicit none
!
      complex(dp), intent(in) :: alpha
!
      integer, intent(in) :: dim_p, dim_q
!
      complex(dp), dimension(dim_p,dim_p,dim_q,dim_q), intent(inout) :: X
!      
      integer :: p, q
!
!$omp parallel do schedule(static) private(p, q)
      do q = 1, dim_q
         do p = 1, dim_p
!
            X(p,p,q,q) = alpha*X(p,p,q,q)
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine scale_complex_4_diagonal_by_complex_1324
!
!
   subroutine scale_real_packed_4_diagonal_by_real(alpha, X, dim_pq)
!!
!!    Scale real packed 4 diagonal by real
!!    Written by Andreas Skeidsvoll, Nov 2019
!!
!!    Scales the diagonal of a real packed 4-dimensional matrix by a real scalar alpha. The
!!    unpacked matrix is symmetric under the exchange of the two first and last dimensions, and has
!!    been packed as
!!
!!    X(dim_p, dim_q, dim_p, dim_q) -> X(dim_pq*(dim_pq+1)/2)
!!
!!    where dim_pq = dim_p*dim_q
!!
      implicit none 
!
      real(dp), intent(in) :: alpha
!
      integer :: dim_pq
!
      real(dp), dimension(dim_pq*(dim_pq+1)/2), intent(inout) :: X 
!
      integer :: pq, pqpq
!
!$omp parallel do private(pq,pqpq)
      do pq = 1, dim_pq
!
         pqpq = pq*(pq+1)/2
!
         X(pqpq) = alpha*X(pqpq)
!
      enddo
!$omp end parallel do
!
   end subroutine scale_real_packed_4_diagonal_by_real
!
!
   subroutine scale_complex_packed_4_diagonal_by_real(alpha, X, dim_pq)
!!
!!    Scale complex packed 4 diagonal by real
!!    Written by Andreas Skeidsvoll, Nov 2019
!!
!!    Scales the diagonal of a complex packed 4-dimensional matrix by a real scalar alpha. The
!!    unpacked matrix is symmetric under the exchange of the two first and last dimensions, and has
!!    been packed as
!!
!!    X(dim_p, dim_q, dim_p, dim_q) -> X(dim_pq*(dim_pq+1)/2)
!!
!!    where dim_pq = dim_p*dim_q
!!
      implicit none 
!
      real(dp), intent(in) :: alpha
!
      integer :: dim_pq
!
      complex(dp), dimension(dim_pq*(dim_pq+1)/2), intent(inout) :: X 
!
      integer :: pq, pqpq
!
      complex(dp) :: alpha_complex
!
      alpha_complex = cmplx(alpha, zero, dp)
!
!$omp parallel do private(pq,pqpq)
      do pq = 1, dim_pq
!
         pqpq = pq*(pq+1)/2
!
         X(pqpq) = alpha_complex*X(pqpq)
!
      enddo
!$omp end parallel do
!
   end subroutine scale_complex_packed_4_diagonal_by_real
!
!
   subroutine scale_complex_packed_4_diagonal_by_complex(alpha, X, dim_pq)
!!
!!    Scale complex packed 4 diagonal by complex
!!    Written by Andreas Skeidsvoll, Nov 2019
!!
!!    Scales the diagonal of a complex packed 4-dimensional matrix by a complex scalar alpha. The
!!    unpacked matrix is symmetric under the exchange of the two first and last dimensions, and has
!!    been packed as
!!
!!    X(dim_p, dim_q, dim_p, dim_q) -> X(dim_pq*(dim_pq+1)/2)
!!
!!    where dim_pq = dim_p*dim_q
!!
      implicit none 
!
      complex(dp), intent(in) :: alpha
!
      integer :: dim_pq
!
      complex(dp), dimension(dim_pq*(dim_pq+1)/2), intent(inout) :: X 
!
      integer :: pq, pqpq
!
!$omp parallel do private(pq,pqpq)
      do pq = 1, dim_pq
!
         pqpq = pq*(pq+1)/2
!
         X(pqpq) = alpha*X(pqpq)
!
      enddo
!$omp end parallel do
!
   end subroutine scale_complex_packed_4_diagonal_by_complex
!
!
   complex(dp) function our_zdotu(n, zx, incx, zy, incy)
!!
!!    Our zdotu
!!    Written by Andreas Skeidsvoll, Nov 2019
!!
!!    A custom zdotu routine, to make zdotu work on Macs
!!
      implicit none
!
      integer :: incx, incy, n
!
      complex(dp) :: zx(*), zy(*)
!
      integer :: i
!
      if ((incx .ne. 0) .and. (incy .ne. 0)) continue
!
      our_zdotu = zero_complex
!
!$omp parallel do schedule(static) private(i) reduction(+:our_zdotu)
      do i = 1, n
!
         our_zdotu = our_zdotu + zx(i)*zy(i)
!
      enddo
!$omp end parallel do
!
   end function our_zdotu
!
!
   subroutine block_diagonalize_symmetric(A,  n_total, n_blocks, block_dim, diagonal, flip_vectors)
!!
!!    Block diagonalize symmetric
!!    Written by Sarai D. Folkestad, Mar 2020
!!
!!    Block diagonalizes the symmetric matrix A (an n_total x n_total matrix)
!!
!!    n_blocks : The number of blocks
!!
!!    block_dim : array of dimensions for the blocks
!!                Elements sum up to n_total
!!
!!    doagonal : Contains the diagonal elements of the 
!!               matrix A after block diagonalization
!!
!!    This routine uses the convention that eigenvectors are flipped
!!    such that the largest element of the eigenvectors are positive.
!!
!!
      implicit none
!
      integer, intent(in) :: n_total, n_blocks
      integer, dimension(n_blocks), intent(in) :: block_dim
!
      real(dp), dimension(n_total, n_total), intent(inout) :: A
      real(dp), dimension(n_total), intent(out) :: diagonal
!
      logical, intent(in), optional :: flip_vectors
!
      real(dp), dimension(:), allocatable :: work
!
      real(dp), dimension(:), allocatable :: eigenvalues
!
      integer  :: info, i, block_, offset, index_max
      real(dp), dimension(1) :: worksize
      logical  :: local_flip_vectors
!
      local_flip_vectors = .true.
      if (present(flip_vectors)) local_flip_vectors = flip_vectors
!
      offset = 1
!
      do block_ = 1, n_blocks
!
         if (block_dim(block_) .gt. 0) then
!
            call mem%alloc(eigenvalues, block_dim(block_))
!  
!           Work size query
!  
            call dsyev('V','U',              &
                        block_dim(block_),   &
                        A(offset, offset),   &
                        n_total,             &
                        eigenvalues,         &
                        worksize,            &
                        -1,                  &
                        info)
!  
            if (info .ne. 0) call output%error_msg('could not perform DSYEV worksize query')
!
!           Diagonalize block
!
            call mem%alloc(work, ceiling(worksize(1)))
!
            call dsyev('V','U',               &
                        block_dim(block_),    &
                        A(offset, offset),    &
                        n_total,              &
                        eigenvalues,          &
                        work,                 &
                        ceiling(worksize(1)), &
                        info)
!
            if (info .ne. 0) then
               call output%error_msg('Diagonalization of block failed.' // &
                                    ' "Dsyev" finished with info: (i0)', ints=[info])
            end if
!
            call mem%dealloc(work, ceiling(worksize(1)))
!
!           Convention for eigenvectors of block, maximum element is positive
!
            if (local_flip_vectors) then
!
               do i = offset, offset + block_dim(block_) - 1
!
                  index_max = maxloc(abs(A(offset : offset + block_dim(block_) - 1, i)), dim=1)
                  index_max = index_max + offset - 1
!
                  if (A(index_max, i) .lt. zero) &
                     call dscal(block_dim(block_), &
                               -one,               &
                               A(offset : offset + block_dim(block_) - 1, i), 1)
!
               enddo
!
            endif
!
!           Setting diagonal
!
            do i = 1, block_dim(block_)
!
               diagonal(i + offset - 1) = eigenvalues(i)
!
            enddo
!
            call mem%dealloc(eigenvalues, block_dim(block_))
!
         elseif (block_dim(block_) .lt. 0) then
!
            call output%error_msg('block diagonalization attempted with negative dimension.')
!
         endif
!
         offset = offset + block_dim(block_)
!
      enddo
!
      if ((offset - 1) .ne. n_total) &
         call output%error_msg('dimensions do not add to n_total in block diagonalization')
!
   end subroutine block_diagonalize_symmetric
!
!
   subroutine gram_schmidt_biorthonormalization(L, R, dim_, n_vectors, threshold)
!!
!!    Gram-Schmidt Biorthonormalization
!!    Written by Alexander C. Paul, Oct 2019
!!
!!    Routine to biorthonormalize two linear independent sets of vectors L and
!!    R following Kohaupt, L., Rocky Mountain J. Math., 44, 1265, (2014)
!!
!!    NB: Order of the R vectors can change from input to output
!!        The R vectors are ordered such that the corresponding L vector 
!!        with non-zero overlap are at the same position
!!
!!    The p-th vector is determined in terms of the previously biorthogonalized
!!    vectors q
!!
!!    Biorthogonal:  L'_p = L_p - sum_q < L_p, R"_q>  L"_q
!!    Biorthonormal: L"_p = < L'_p, R_p >^(-1)  L'_p = < L'_p, R"_p >^(-1)  L'_p
!!
!!    Biorthonormal: R"_p = R_p - sum_q < R_p, L"_q >  R"_q
!!
!!    Modified by Andreas S. Skeidsvoll, Nov 2020
!!    Biorthogonalization against previous vectors is now done irrespectively
!!    of the number of zero overlaps with the following vectors (unless all
!!    overlaps are zero).
!!
      implicit none
!
      integer, intent(in) :: dim_, n_vectors
      real(dp), dimension(dim_, n_vectors), intent(inout)  :: L, R
      real(dp), intent(in) :: threshold
!
      real(dp) :: overlap, max_overlap, temp
!
      integer :: p, q, index_max_overlap
!
      real(dp), external :: ddot
!
      do p = 1, n_vectors ! looping through the left vectors
!
!        :: Biorthogonalize left vector against previous right vectors ::
!        ----------------------------------------------------------------
!
         do q = 1, p-1
!
!           Construct biorthogonal left vector:
!              L'_p = L_p - sum_q < L_p, R"_q>  L"_q
!
            overlap = ddot(dim_,      &
                           L(:,p), 1, &
                           R(:,q), 1)
!
            call daxpy(dim_,      &
                       -overlap,  &
                       L(:,q), 1, &
                       L(:,p), 1)
!
         enddo
!
!        Look for the right vector with the maximum overlap with L_p. Start at
!        q = p, as the overlaps with the first (p-1) vectors are zero.
!
         index_max_overlap = p
         max_overlap = zero
!
         do q = p, n_vectors ! right vectors
!
            overlap = ddot(dim_,      &
                           L(:,p), 1, &
                           R(:,q), 1)
!
            if (abs(overlap) .gt. abs(max_overlap)) then
!
               index_max_overlap = q
               max_overlap = overlap
!
            end if
!
         end do
!
         if (abs(max_overlap) .lt. threshold) then
!
            call output%printf('m', 'Overlaps between left vector with &
                                     &index (i0) and right vectors are &
                                     &less than threshold: (e8.3).', &
                               reals=[threshold], ints=[p])
!
            call output%error_msg('Trying to binormalize nonoverlapping vectors.')
!
         end if
!
!        :: Biorthonormalization ::
!        --------------------------
!
!        L"_p = < L'_p, R_p >^(-1)  L'_p
!
         call dscal(dim_, one/max_overlap, L(:,p), 1)
!
         do q = 1, p-1
!
!           Construct biorthonormal right vector:
!              R"_p = R_p - sum_q < R_p, L"_q >  R"_q
!
            overlap = ddot(dim_,                  &
                           R(:,index_max_overlap), 1, &
                           L(:,q), 1)
!
            call daxpy(dim_,      &
                       -overlap,  &
                       R(:,q), 1, &
                       R(:,index_max_overlap), 1)
!
         end do
!
!        Exchange R_p and R_index_max_overlap to achieve the correct ordering
!
         if(p .ne. index_max_overlap) then
!
!$omp parallel do private(temp, q)
            do q = 1, dim_
!
               temp   = R(q,p)
               R(q,p) = R(q,index_max_overlap)
               R(q,index_max_overlap) = temp
!
            end do
!$omp end parallel do
!
         end if
!
      end do ! Loop over p (left vectors)
!
   end subroutine gram_schmidt_biorthonormalization
!
!
   subroutine diagonalize_symmetric(A, dim_, diagonal)
!!
!!    Diagonalize symmetric
!!    Written by Sarai D. Folkestad, Jun 2020
!!
!!    Diagonalizes the symmetric matrix A (an dim_ x dim_ matrix)
!!
!!    diagonal : Contains the diagonal elements of the 
!!               matrix A after block diagonalization
!!
!!    This routine uses the convention that eigenvectors are flipped
!!    such that the largest element of the eigenvectors are positive.
!!
!!    This is a wrapper using the block_diagonalization routine
!!    which in turn uses the lapack dsyeev routine
!!
      implicit none
!
      integer, intent(in) :: dim_
!
      real(dp), dimension(dim_, dim_), intent(inout)   :: A
      real(dp), dimension(dim_), optional, intent(out) :: diagonal
!
      real(dp), dimension(:), allocatable :: eigenvalues
!
      if (present(diagonal)) then
!
         call block_diagonalize_symmetric(A,  dim_, 1, [dim_], diagonal)
!
      else
!
         call mem%alloc(eigenvalues, dim_)
         call block_diagonalize_symmetric(A,  dim_, 1, [dim_], eigenvalues)
         call mem%dealloc(eigenvalues, dim_)
!
      endif
!
   end subroutine diagonalize_symmetric
!
!
   subroutine constant_array(x, n, const)
!!
!!    Constant array 
!!    Written by Sarai D. Folkestad, 2021
!!
!!    Sets the array x of length n to const.
!!
      implicit none 
!
      integer, intent(in) :: n 
!
      real(dp), dimension(n), intent(out) :: x 
!
      real(dp), intent(in) :: const
!
      integer :: I 
!
!$omp parallel do private(I) schedule(static)
      do I = 1, n 
!
         x(I) = const
!
      enddo
!$omp end parallel do
!
   end subroutine constant_array
!
!
   subroutine generalized_diagonalization_symmetric(A, S, dim_, e)
!!
!!    Generalized diagonalization symmetric matrix
!!    Written by Sarai D. Folkestad,  2020
!!
!!    Solves the linear equation
!!
!!       A X = S X e.
!!    
!!    On exit the vectors (X) are stored in A
!!    On exit S is overwritten
!! 
      implicit none
!
      integer, intent(in) :: dim_
!
      real(dp), dimension(dim_, dim_), intent(inout)   :: A
      real(dp), dimension(dim_, dim_), intent(inout)   :: S
      real(dp), dimension(dim_), optional, intent(out) :: e
!
      real(dp), dimension(:), allocatable :: work
!
      integer :: info
!
      info = 0
!
      call mem%alloc(work, 4*dim_)
!
      call dsygv(1, 'V', 'L', &
                  dim_,       &
                  A,          & 
                  dim_,       &
                  S,          &
                  dim_,       &
                  e,          &
                  work,       &
                  4*(dim_),   &
                  info)
!
      call mem%dealloc(work, 4*dim_)   
!
      if (info .ne. 0) call output%error_msg('in generalized diagonalization (array_utilities.F90)')   
!
   end subroutine generalized_diagonalization_symmetric
!
!
   subroutine copy_integer(X, Y, n, alpha, beta)
!!
!!    Copy integer
!!    Written by Alexander C. Paul, Feb 2021
!!
!!    Sets Y as:
!!
!!       Y = X
!!
!!    X and Y are vectors of length n
!!
      implicit none
!
      integer, intent(in) :: n
!
      integer, dimension(n), intent(out) :: Y
      integer, dimension(n), intent(in) :: X
!
      integer, optional, intent(in) :: alpha, beta
      integer :: alpha_, beta_
!
      integer :: i
!
      alpha_ = 1
      beta_  = 0
      if (present(alpha)) alpha_ = alpha
      if (present(beta))  beta_  = beta
!
      if (beta_ .eq. 0) then
!
!$omp parallel do private(i)
         do i = 1, n
!
            Y(i) = alpha_*X(i)
!
         enddo
!$omp end parallel do 
!
      else
!
!$omp parallel do private(i)
         do i = 1, n
!
            Y(i) = beta_*Y(i) + alpha_*X(i)
!
         enddo
!$omp end parallel do 
!
      end if
!
   end subroutine copy_integer
!
end module array_utilities
