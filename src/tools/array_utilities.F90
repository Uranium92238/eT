!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2022 the authors of eT
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
   interface get_trace
      procedure :: get_trace_r, &
                   get_trace_c
   end interface get_trace
!
   interface dsfrk_
      procedure :: dsfrk_r, &
                   dsfrk_c
   end interface dsfrk_
!
contains
!
   subroutine dsfrk_r(n, k, alpha, A, beta, C)
!!
!!    dsfrk real
!!    Written by Sarai D. Folkestad
!!
!!    Wrapper for lapack dsfrk
!!
!!       C = alpha*A**T*A + beta*C,
!!
!!    where C is returned in upper rectangular
!!    full packed format.
!!
      implicit none
!
      integer :: n, k
      real(dp) :: alpha, beta
      real(dp), dimension(k,n) :: A
      real(dp), dimension((n + mod(n+1,2)), n*(n+1)/2) :: C
!
      call dsfrk('n', 'U', 'T', n, k, alpha, A, k, beta, C)
!
   end subroutine dsfrk_r
!
   subroutine dsfrk_c(n, k, alpha, A, beta, C)
!!
!!    dsfrk complex
!!    Written by Rolf H. Myhre
!!
!!    lapack dsfrk
!!
!!       C = alpha*A**T*A + beta*C,
!!
!!    does not exist for complex A and C. In stead
!!    must make two calls to dsyrk and one call to dgemm.
!!
!!    C is returned in upper rectangular
!!    full packed format.
!!
      implicit none
!
      integer :: n, k
      complex(dp) :: alpha, beta
      complex(dp), dimension(k,n) :: A
      complex(dp), dimension((n + mod(n+1,2)), n*(n+1)/2) :: C
!
      integer :: n1, n2, off
!
      n1 = n/2
      n2 = n - n1
      off = mod(n+1,2)
!
      call zsyrk('L', 'T', n1, k, alpha, A, k, &
                 beta, C(n2+1+off,1), n+off)
!
      call zsyrk('U', 'T', n2, k, alpha, A(1, n2+off), k, &
                 beta, C(n1+1,1), n+off)
!
      call zgemm('T', 'n', n1, n2, k, alpha, A, &
                 k, A(1, n2+off), k, beta, C, n+off)
!
   end subroutine dsfrk_c
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
      call zero_array(x,n**2)
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
   subroutine full_cholesky_decomposition(matrix, cholesky_vectors, dim_, &
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
      call dcopy(dim_**2, matrix, 1, cholesky_vectors, 1)
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
   end subroutine full_cholesky_decomposition
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
   function get_root_mean_square(X, n) result(rms)
!!
!!    Get root mean square
!!    Written by Eirik F. Kjønstad, 2021
!!
!!    Get root-mean-square, RMS = sqrt(X^T X / n)
!!
      implicit none
!
      integer, intent(in) :: n
      real(dp), dimension(n), intent(in) :: X
!
      real(dp) :: ddot, rms
!
      rms = ddot(n, X, 1, X, 1)/real(n, kind=dp)
!
      rms = sqrt(rms)
!
   end function get_root_mean_square
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
   function zdot(n, zx, incx, zy, incy) result(dot)
!!
!!    zdot
!!    Written by Andreas Skeidsvoll, Nov 2019
!!
!!    A custom zdotu routine, to make zdotu work on Macs
!!
      implicit none
!
      integer :: incx, incy, n
!
      complex(dp) :: zx(*), zy(*)
      complex(dp) :: dot
!
      integer :: i
!
      if ((incx .ne. 0) .and. (incy .ne. 0)) continue
!
      dot = zero_complex
!
!$omp parallel do schedule(static) private(i) reduction(+:dot)
      do i = 1, n
!
         dot = dot + zx(i)*zy(i)
!
      enddo
!$omp end parallel do
!
   end function zdot
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
   function get_trace_r(x, n) result(trace)
!!
!!    Compute trace real
!!    Written by Alexander C. Paul, Apr 2020
!!
!!    Computes trace of a quadratic (n x n) array x.
!!
      implicit none
!
      integer, intent(in) :: n
!
      real(dp), dimension(n,n), intent(in) :: x
!
      real(dp) :: trace
      integer  :: p
!
      trace = zero
!
!$omp parallel do private(p) schedule(static) reduction(+:trace)
      do p = 1, n
!
         trace = trace + x(p,p)
!
      enddo
!$omp end parallel do
!
   end function get_trace_r
!
!
   function get_trace_c(x, n) result(trace)
!!
!!    Compute trace complex
!!    Written by Alexander C. Paul, Apr 2020
!!
!!    Computes trace of a quadratic (n x n) array x.
!!
      implicit none
!
      integer, intent(in) :: n
!
      complex(dp), dimension(n,n), intent(in) :: x
!
      complex(dp) :: trace
      integer :: p
!
      trace = zero
!
!$omp parallel do private(p) schedule(static) reduction(+:trace)
      do p = 1, n
!
         trace = trace + x(p,p)
!
      enddo
!$omp end parallel do
!
   end function get_trace_c
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
!
   subroutine add_to_subblock(scalar, x, y, p_range, q_range)
!!
!!    Add block
!!    Written by Alexander C. Paul, Apr 2021
!!
!!    Add elements of x to the block of y defined by first/last p and q.
!!    y is assumed to
!!
      use range_class
!
      implicit none
!
      real(dp), intent(in) :: scalar
      class(range_), intent(in) :: p_range, q_range
!
      real(dp), dimension(:,:), intent(inout) :: y
      real(dp), intent(in), dimension(p_range%length, q_range%length) :: x
!
      integer :: p, q, pp, qq
!
!$omp parallel do schedule(static) private(p,q, pp, qq)
      do q = 1, q_range%length
         qq = q_range%first + (q-1)*q_range%step
!
         do p = 1, p_range%length
            pp = p_range%first + (p-1)*p_range%step
!
            y(pp,qq) = y(pp,qq) + scalar*x(p,q)
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_to_subblock
!
!
   subroutine calculate_pseudoinverse(Ainv, A, m, n, tau)
!!
!!    Calculate pseudoinverse
!!    Written by Eirik F. Kjønstad, 2021
!!
!!    The pseudoinverse is computed via an SVD:
!!
!!       A    = U S V^T
!!       Ainv = V S-1 U^T
!!
!!    In S-1, diagonals with singular values (S) below tau*max(S) are set to zero.
!!
      implicit none
!
      integer, intent(in) :: m, n
!
      real(dp), dimension(m,n), intent(in)  :: A
      real(dp), dimension(n,m), intent(out) :: Ainv
!
      real(dp), intent(in) :: tau
!
      real(dp), dimension(:), allocatable   :: S, Sinv
      real(dp), dimension(:,:), allocatable :: U, VT, Acopy
!
      real(dp), dimension(:), allocatable :: work
      integer :: lwork
!
      integer :: error_integer, k, l, j, min_mn
!
      real(dp) :: max_singular
!
      min_mn = min(m,n)
!
      call mem%alloc(Acopy, m, n)
      call copy_and_scale(one, A, Acopy, m*n)
!
      lwork = 2000 * min_mn
      call mem%alloc(work, lwork)
!
      call mem%alloc(S, min_mn)
      call mem%alloc(U, m, min_mn)
      call mem%alloc(VT, min_mn, n)
!
      error_integer = 0
!
      call dgesvd('S', 'S', &
                   m, &
                   n, &
                   Acopy, &
                   m, &
                   S, &
                   U, &
                   m, &
                   VT, &
                   min_mn, &
                   work, &
                   lwork, &
                   error_integer)
!
      call mem%dealloc(work, lwork)
      call mem%dealloc(Acopy, m, n)
!
      if (error_integer .ne. 0) &
         call output%error_msg('SVD failed in pseudoinverse: (i0)', ints=[error_integer])
!
      max_singular = maxval(S, min_mn)
!
      call mem%alloc(Sinv, min_mn)
!
      do k = 1, min_mn
!
         if (S(k) .lt. tau*max_singular) then
!
            Sinv(k) = zero
!
         else
!
            Sinv(k) = one/S(k)
!
         endif
!
      enddo
!
      call mem%dealloc(S, min_mn)
!
!     (Ainv)_kl = V_kj S-1_j U^T_jl
!
      call zero_array(Ainv, m*n)
!
!$omp parallel do schedule(static) private(k,l,j)
      do k = 1, n
         do l = 1, m
            do j = 1, min_mn
!
               Ainv(k,l) = Ainv(k,l) + VT(j,k)*Sinv(j)*U(l,j)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(U, m, min_mn)
      call mem%dealloc(Sinv, min_mn)
      call mem%dealloc(VT, min_mn, n)
!
   end subroutine calculate_pseudoinverse
!
!
   function get_euclidean_distance(point1, point2, N) result(distance)
!!
!!    Get euclidean distance
!!    Written by Alexander C. Paul, Nov 2021
!!
      implicit none
!
      integer, intent(in) :: n
      real(dp), dimension(n), intent(in) :: point1, point2
!
      real(dp) :: distance
!
      integer  :: i
!
      distance = zero
!
      do i = 1, n
         distance = distance + (point1(i) - point2(i))**2
      end do
!
      distance = sqrt(distance)
!
   end function get_euclidean_distance
!
!
   function find_location(array, value_, n) result(index)
!!
!!    Find location
!!    Written by Alexander C. Paul, Jan 2022
!!
!!    Returns index of element of array equal to value_
!!    If no identical element is found the result will be 0.
!!
!!    Required because GCC 8 does not know findloc
!!
      implicit none
!
      integer, intent(in) :: n
      character(len=*), dimension(n), intent(in) :: array
      character(len=*), intent(in) :: value_
!
      integer :: index
!
      do index = 1, n
         if (value_ == array(index)) return
      end do
!
      index = 0
!
   end function find_location
!
!
end module array_utilities
