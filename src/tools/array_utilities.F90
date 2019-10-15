!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2019 the authors of eT
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
!!    Written by Sarai D. Folkstad and Eirik F. Kjønstad, May 2018
!!
!!    This module contains routines that perform various operations on arrays
!!    and that do not belong elsewhere (such as in the reordering module,
!!    where all such routines are gathered for convenience).
!!
!
   use parameters
   use global_out, only : output
   use memory_manager_class, only : mem
!
   implicit none
!
contains
!
!
   integer function get_max_index(x, n)
!!
!!    Get max index
!!    Written by Eirik F. Kjønstad, 2018
!!
      implicit none
!
      integer, intent(in) :: n
      real(dp), dimension(n), intent(in) :: X
!
      integer :: I
      real(dp)     :: max_val
!
      get_max_index = 1
      max_val = X(1)
      do I = 2, n
!
         if (X(I) .gt. max_val) then
!
            get_max_index = I
            max_val = X(I)
!
         endif
!
      enddo
!
   end function get_max_index
!
!
   subroutine zero_array(X, n)
!!
!!    Zero array 
!!    Written by Eirik F. Kjønstad, 2018 
!!
      implicit none 
!
      integer, intent(in) :: n 
!
      real(dp), dimension(n), intent(out) :: X 
!
      integer :: I 
!
!$omp parallel do private(I) schedule(static)
      do I = 1, n 
!
         X(I) = zero 
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
   logical function is_significant(vec, n, threshold, screening)
!!
!!    Is vector significant ?
!!    Written by Eirik F. Kjønstad and Sarai D. Folkstad, June 2018
!!
!!    Returns true if all elements are below threshold
!!
      implicit none
!
      integer, intent(in) :: n
!
      real(dp), dimension(n), intent(in)  :: vec
      real(dp), dimension(n), intent(in), optional  :: screening
!
      real(dp), intent(in)  :: threshold
!
      integer :: i = 0
!
      is_significant = .false.
!
      if (present(screening)) then
         do i = 1, n
!
            if (abs(vec(i)*screening(i)) .gt. threshold) then
!
               is_significant = .true.
               return
!
            endif
!
         enddo
      else
!
         do i = 1, n
!
            if (abs(vec(i)) .gt. threshold) then
!
               is_significant = .true.
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
   integer function n_significant(vec, n, threshold)
!!
!!    Number of significant in vector
!!    Written by Eirik F. Kjønstad and Sarai D. Folkstad, June 2018
!!
!!    Returns the number of elements in vector larger than threshold
!!
      implicit none
!
      integer, intent(in) :: n
!
      real(dp), dimension(n), intent(in)  :: vec
!
      real(dp), intent(in)  :: threshold
!
      integer :: i = 0
!
      n_significant = 0
!
      do i = 1, n
!
         if (abs(vec(i)) .gt. threshold) then
!
            n_significant = n_significant + 1
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
!!    Cuts the significant blocks out of a vector and places them in a reduced size
!!    vector
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
!!    Reduce vector
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Cuts the significant blocks out of a vector and places them in a reduced size
!!    vector
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
!!    Cuts the significant row blocks out of a array and places them in a reduced size
!!    array
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
   subroutine reduce_array_column(array, array_reduced, block_firsts, block_significant, n_blocks, dim_, dim_reduced, rows)
!!
!!    Reduce array column
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Cuts the significant column blocks out of an array and places them in a reduced size array
!!
      implicit none
!
      integer, intent(in) :: dim_, dim_reduced, n_blocks, rows
!
      logical, dimension(n_blocks), intent(in) :: block_significant
      integer, dimension(n_blocks + 1), intent(in) :: block_firsts
!
      real(dp), dimension(rows, dim_), intent(in) :: array
      real(dp), dimension(rows, dim_reduced), intent(out) :: array_reduced
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
            array_reduced(:, current_pos : current_pos + size_ - 1) = array(:, first : last)
!
            current_pos = current_pos + size_
!
         endif
!
      enddo
!
   end subroutine reduce_array_column
!
!
   subroutine reduce_array_int(array, array_reduced, block_firsts, block_significant, n_blocks, dim_, dim_reduced, columns)
!!
!!    Reduce array
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Cuts the significant row blocks out of a array and places them in a reduced size
!!    array
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
            array_reduced(current_pos : (current_pos + size_ - 1), :) = array(first : last, :)
!
            current_pos = current_pos + size_
!
         endif
!
      enddo
!
   end subroutine reduce_array_int
!
!
   subroutine full_cholesky_decomposition(matrix, cholesky_vectors, dim_, n_vectors,&
                                        threshold, used_diag)
!!
!!    Cholesky decomposition,
!!    Written by Sarai Dery Folkestad, June 2017
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
      integer, dimension(dim_), optional, intent(out) :: used_diag
!
      integer :: i, j, k, index_max
      real(dp) :: max_diagonal
!
      real(dp), parameter :: tolerance = 1.0d-10
!
      if (present(used_diag)) used_diag = 0
!
!     Looping over the number of cholesky vectors
!
      do i = 1, dim_
         n_vectors = i
!
!        Find the maximum diagonal
!
         index_max = 0
         max_diagonal = 0.0d0
!
         do j = 1, dim_
!
            if (abs(matrix(j, j)) .gt. abs(max_diagonal)) then
!
               max_diagonal = matrix(j,j)
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
               write(output%unit,*)'Error: Found negative diagonal in cholesky decomposition.'
               stop
!
            endif
         endif
!
         if (abs(max_diagonal) .lt. threshold) then
!
            n_vectors = n_vectors - 1
            return
         else
            if (present(used_diag)) used_diag(n_vectors) = index_max
         endif
!
!        Cholesky vectors
!
         do j = 1, dim_
!
            cholesky_vectors(j,i) = matrix(j, index_max)/sqrt(max_diagonal)
!
         enddo
!
!        Subtract from matrix
!
         do j = 1, dim_
            do k = 1, dim_
!
               matrix(k,j) = matrix(k,j) - cholesky_vectors(k,i)*cholesky_vectors(j,i)
!
            enddo
         enddo
!
         do j = 1, dim_
            matrix(j,index_max) = 0.0D0
            matrix(index_max,j) = 0.0D0
         enddo
!
      enddo

   end subroutine full_cholesky_decomposition
!
!
   subroutine full_cholesky_decomposition_effective(matrix, cholesky_vectors, dim_, n_vectors,&
                                        threshold, used_diag)
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
               write(output%unit,*)'Error: Found negative diagonal in cholesky decomposition.'
               stop
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
               if (diagonal(j) .lt. min_diagonal) min_diagonal = diagonal(j)

!
            enddo
!
            write(output%unit, '(t3, a46, e12.4)') 'The smallest diagonal after decomposition is: ', min_diagonal
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
      write(output%unit, '(t3, a46, e12.4)') 'The smallest diagonal after decomposition is: ', min_diagonal
!
      call mem%dealloc(diagonal, dim_)
!
   end subroutine full_cholesky_decomposition_effective
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
!!    Routine is used for decomposition of density to construct active
!!    orbitals.
!!
!!    The number of pivots may specified through the optional
!!    argument n_vectors_requested  
!!
!!    On exit matrix_xy = matrix_xy - sum_J L_xJ*LyJ
!!
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
               write(output%unit,*)'Error: Found negative diagonal in cholesky decomposition.'
               stop
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
                        one,              &
                        cholesky_vectors, &
                        dim_,             &
                        cholesky_vectors, &
                        dim_,             &
                        -one,             &
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
 subroutine full_cholesky_decomposition_system(matrix, cholesky_vectors, dim_, n_vectors, &
                                                   threshold, used_diag)
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
      allocate(work(2*dim_))
!
!     DPSTRF computes the Cholesky factorization with complete pivoting
!     of a real symmetric positive semidefinite matrix.
!
      call dpstrf('L',      &
            dim_,              &
            cholesky_vectors, &
            dim_,              &
            used_diag,        &
            n_vectors,        &
            threshold,        &
            work,             &
            info)
!
      deallocate(work)
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
         write(*,*) info
         stop 'Cholesky decomposition failed! Something wrong in call to dpstrf'
      end if
!
   end subroutine full_cholesky_decomposition_system
!
!
   subroutine inv(Ainv, A, n)
!!
!!    Invert matrix A
!!
      implicit none
!
      integer, intent(in) :: n
!
      real(dp), dimension(n,n), intent(inout) :: A
      real(dp), dimension(n,n), intent(inout) :: Ainv

      real(dp), dimension(n) :: work  ! work array for LAPACK
      integer, dimension(n) :: ipiv   ! pivot indices
      integer :: info
!
!     Store A in Ainv to prevent it from being overwritten by LAPACK
!
      Ainv = A
!
!     DGETRF computes an LU factorization of a general M-by-N matrix A
!     using partial pivoting with row interchanges.
!
      call DGETRF(n, n, Ainv, n, ipiv, info)
!
      if (info /= 0) then
         stop 'Matrix is numerically singular!'
      end if
!
!     DGETRI computes the inverse of a matrix using the LU factorization
!     computed by DGETRF.
!
      call DGETRI(n, Ainv, n, ipiv, work, n, info)
!
      if (info /= 0) then
         stop 'Matrix inversion failed!'
      end if
!
   end subroutine inv
!
!
   subroutine inv_lower_tri(Ainv, A, n)
!!
!!    Invert lower triagonal matrix A
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
!     DTRTRI computes the inverse of a real upper or lower triangular
!     matrix A.
!
      call DTRTRI('l','n', n, Ainv, n, info)
!
      if (info /= 0) then
         write(output%unit, *) 'Error: matrix inversion failed!', info
         stop
      end if
!
   end subroutine inv_lower_tri
!
!
   subroutine symmetrize(M, n)
!!
!!    Symmetrize matrix
!!    Written by Eirik F. Kjønstad, 2018
!!
!!     M <- 1/2 * (M + M^T)
!!
      implicit none
!
      integer, intent(in) :: n
!
      integer :: i, j
!
      real(dp), dimension(n, n) :: M
      real(dp), dimension(:, :), allocatable :: MT
!
      call mem%alloc(MT, n, n)
!
      MT = transpose(M)
!
      do j = 1, n
         do i = 1, n
!
            M(i, j) = half*(M(i, j) + MT(i, j))
!
         enddo
      enddo
!
      call mem%dealloc(MT, n, n)
!
   end subroutine symmetrize
!
!
   subroutine anti_symmetrize(M, n)
!!
!!    Antisymmetrize matrix
!!    Written by Eirik F. Kjønstad, 2018
!!
!!     M <- 1/2 * (M - M^T)
!!
      implicit none
!
      integer, intent(in) :: n
!
      integer :: i, j
!
      real(dp), dimension(n, n) :: M
      real(dp), dimension(:, :), allocatable :: MT
!
      call mem%alloc(MT, n, n)
!
      MT = transpose(M)
!
      do j = 1, n
         do i = 1, n
!
            M(i, j) = half*(M(i, j) - MT(i, j))
!
         enddo
      enddo
!
      call mem%dealloc(MT, n, n)
!
   end subroutine anti_symmetrize
!
!
   real(dp) function get_abs_max(X, n)
!!
!!    Get absolute maximum value of vector
!!    Written by Eirik F. Kjønstad, 2018
!!
!!    Returns the largest element of |X|, i.e. the largest 
!!    element in terms of absolute value
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
!!    Returns the largest element of |X|, i.e. the largest 
!!    element in terms of absolute value
!!
!!    The routine is a copy of the function get_abs_max 
!!    written by Eirik F. Kjønstad. S.D.F added the
!!    index_ stuff
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
   subroutine sandwich(X, A, B, n, left)
!!
!!    Sandwich
!!    Written by Eirik F. Kjønstad, 2018
!!
!!    Overwrites the n-by-n matrix X with a sandwich-ed X:
!!
!!       X <- A^T X B
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

!
      call mem%dealloc(tmp, n, n)
!
   end subroutine sandwich
!
!
   subroutine symmetric_sandwich(Xr, X, A, m, n)
!!
!!    Symmetric sandwich
!!    Written by Eirik F. Kjønstad, 2018
!!
!!    Overwrites the n-by-n matrix X with a sandwich-ed X:
!!
!!       X <- A^T X A
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
   subroutine symmetric_sandwich_right(Xr, X, A, m, n)
!!
!!    Symmetric sandwich using right transposition
!!    Written by Eirik F. Kjønstad, 2018
!!
!!    Overwrites the n-by-n matrix X with a sandwich-ed X:
!!
!!       X <- A X A^T
!!
      implicit none
!
      integer, intent(in) :: m, n
!
      real(dp), dimension(m, n), intent(in) :: A
!
      real(dp), dimension(n, n) :: X
      real(dp), dimension(m, m) :: Xr
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
   end subroutine symmetric_sandwich_right
!
!
   subroutine commute(A, B, AcB, n)
!!
!!    Commute 
!!    Written by Eirik F. Kjønstad, 2018
!!
!!    Calculates the commutator of A and B, two square matrices 
!!    of dimension n, and places the result in AB. 
!!
      implicit none 
!
      integer, intent(in) :: n 
!
      real(dp), dimension(n, n), intent(in) :: A 
      real(dp), dimension(n, n), intent(in) :: B
!
      real(dp), dimension(n, n) :: AcB ! [A, B] = AB - BA on exit 
!
      call dgemm('N', 'N', &
                  n,       &
                  n,       &
                  n,       &
                  one,     &   
                  A,       &
                  n,       &
                  B,       &
                  n,       &
                  zero,    &
                  AcB,     & ! AcB = AB 
                  n)
!
      call dgemm('N', 'N', &
                  n,       &
                  n,       &
                  n,       &
                  -one,    &   
                  B,       &
                  n,       &
                  A,       &
                  n,       &
                  one,     &
                  AcB,     & ! AcB = AcB - BA = AB - BA 
                  n)      
!
   end subroutine commute
!
!
   real(dp) function get_l2_norm(X, n)
!!
!!    Get L^2 norm 
!!    Written by Eirik F. Kjønstad, Aug 2018 
!!
!!    Returns the L^2 norm of the n-dimensional X vector, 
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
   subroutine print_vector(A, n, indent)
!!
!!    Print vector 
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Suitable to print vector of size < 1000.
!!
      implicit none 
!
      integer, intent(in) :: n
!
      real(dp), dimension(n), intent(in) :: A 
!
      integer :: I 
!
      character(len=*)   :: indent ! indentation
      character(len=255) :: adv    ! advance
      character(len=255) :: frmt   ! format
      character(len=255) :: sep    ! column separation
!
      write(output%unit, *)
      sep = trim(indent)
      adv = 'no'
      do I = 1, n 
!
         frmt = '(t' // trim(sep) // ', i3, f18.12)'
         write(output%unit, frmt, advance=trim(adv)) I, A(I) 
!
         if (mod(I, 3) .eq. 2) then 
!
            adv = 'yes'
            sep = '5'
!
         elseif (mod(I, 3) .eq. 0) then
!
            adv = 'no'
            sep = indent
!
         else
!
            adv = 'no'
            sep = '5'
!
         endif
!
      enddo
      write(output%unit, *)
!
   end subroutine print_vector
!
!
   subroutine trans(A, A_trans, dim_)
!!
!!    Transpose 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
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
   end subroutine trans
!
!
   subroutine copy_and_scale(alpha, X, Y, n)
!!
!!    Copy and scale
!!    Written by Sarai D. Folkestad, May 2019
!!
!!    Y = alpha*X
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
   subroutine print_matrix(name_, matrix, idim1, idim2)
!!    
!!    Print square matrix 
!!    Written by tommaso giovannini, march 2019
!!    
      implicit none
!
      character(len=*) :: name_
!      
      integer :: idim1, idim2, i, j, k, nelmfinal, nvec, start
!      
      real(kind=dp), dimension(idim1,idim2)    :: matrix
      real(kind=dp), dimension(:,:,:), allocatable :: dummy
      real(kind=dp), dimension(:,:,:), allocatable :: dummy2
!     
      character(len=2)  :: start_string
      character(len=12) :: frmt0
      character(len=13) :: frmt1
      character(len=10) :: frmt2
!
      call output%warning_msg('Printing all elements below 1.0d-6 as 0.0d0')
!
      frmt0="(t5,65('='))"
      frmt1="(3x,5(7x,i5))"
!
!     count how many vectors are there
!
      nvec      = ceiling(float(idim2)/5.0d0)
!         
      if(mod(idim2,5).eq.0) then 
!         
         nelmfinal = 5
!         
      else 
!         
         nelmfinal = mod(idim2,5)
!         
      endif
!         
      call mem%alloc(dummy, idim1, nvec-1, 5)
!      
      call mem%alloc(dummy2, idim1, 1, nelmfinal)
!      
      if(nvec.ne.1) then
!$omp parallel do collapse(3) 
         do i = 1, nvec-1
!         
            do j = 1, idim1
!         
               do k = 1, 5
!         
                  if(abs(matrix(j,k+(i-1)*5)).gt.1.0d-6) then
!         
                     dummy(j,i,k) = matrix(j,k+(i-1)*5)
!         
                  else
!         
                     dummy(j,i,k) = 0.0d0
!         
                  endif
!         
               enddo
!         
            enddo
!         
         enddo
!$omp end parallel do
      endif
!         
      do j = 1, idim1
!         
         do k = 1, nelmfinal
!         
            if(abs(matrix(j,k+(nvec-1)*5)).gt.1.0d-6) then
!         
               dummy2(j,1,k) = matrix(j,k+(nvec-1)*5)
!         
            else
!         
               dummy2(j,1,k) = 0.0d0
!         
            endif
!         
         enddo
!         
      enddo
!         
!     place the name in the middle of the string
!
      start = 32-len(name_)/2
      write(start_string,'(i2)') start
      frmt2  ="(t5,"// start_string //"x,a)"
!      
      write(output%unit,*)
      write(output%unit,frmt0 )
      write(output%unit,frmt2) name_
      write(output%unit,frmt0 )
!      
      do i = 1, nvec
!      
         if(i.ne.nvec) write(output%unit,frmt1) (k+(i-1)*5,k=1,5)
         if(i.eq.nvec) write(output%unit,frmt1) (k+(i-1)*5,k=1,nelmfinal)
!         
         do j = 1, idim1
!         
            if(i.ne.nvec) write(output%unit,'(3x,i4,2x,5(e11.4,1x))') j,(dummy(j,i,k), k = 1, 5 )
!         
            if(i.eq.nvec) then 
!         
               write(output%unit,'(3x,i4,2x,5(e11.4,1x))') j,(dummy2(j,1,k), k = 1, nelmfinal )
!         
            endif   
!         
         enddo
!         
      enddo
!         
      write(output%unit,frmt0 )
!         
      call mem%dealloc(dummy, idim1, nvec-1, 5)
!      
      call mem%dealloc(dummy2, idim1, 1, nelmfinal)
!
!
   end subroutine print_matrix
!
!
end module array_utilities
