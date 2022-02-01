!
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
module reordering
!
!!
!!    Reordering module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2018
!!
!!    Here are all routines that are used to re-sort arrays and (e.g., sort_123_to_132)
!!    performing re-sort-addition (add_132_to_123), along with routines that both squareup
!!    (e.g., t2am -> t_ai_bj) and re-sort (t_ai_bj -> t_aj_bi) in a single step; also included
!!    are routines that squareup (t_aibj -> t_ai_bj) and packin (t_ai_bj -> t_aibj) and
!!    related routines.
!!
!
   use parameters
!
   implicit none
!
   interface sort_12_to_21
      procedure :: sort_12_to_21_r, &
                   sort_12_to_21_c, &
                   sort_12_compound_to_21
   end interface sort_12_to_21
!
   interface add_21_to_12
      procedure :: add_21_to_12, &
                   add_21_to_12_complex
   end interface add_21_to_12
!
   interface symmetric_sum
      procedure :: symmetric_sum_real, &
                   symmetric_sum_complex, &
                   symmetric_sum_4, &
                   symmetric_sum_4_complex
   end interface
!
   interface sort_123_to_312
      procedure :: sort_123_to_312, &
                   sort_123_to_312_complex
   end interface sort_123_to_312
!
   interface sort_123_to_132
      procedure :: sort_123_to_132, &
                   sort_123_to_132_complex
   end interface sort_123_to_132
!
   interface sort_123_to_213
      procedure :: sort_123_to_213, &
                   sort_123_to_213_complex
   end interface sort_123_to_213
!
   interface add_132_to_123
      procedure :: add_132_to_123, &
                   add_132_to_123_complex
   end interface add_132_to_123
!
   interface sort_1234_to_3412
      procedure :: sort_1234_to_3412, &
                   sort_1234_to_3412_complex
   end interface sort_1234_to_3412
!
   interface sort_1234_to_4132
      procedure :: sort_1234_to_4132, &
                   sort_1234_to_4132_complex
   end interface sort_1234_to_4132
!
   interface sort_1234_to_4123
      procedure :: sort_1234_to_4123, &
                   sort_1234_to_4123_complex
   end interface sort_1234_to_4123
!
   interface sort_1234_to_3142
      procedure :: sort_1234_to_3142, &
                   sort_1234_to_3142_complex
   end interface sort_1234_to_3142
!
   interface squareup_and_sort_1234_to_2413
      procedure :: squareup_and_sort_1234_to_2413, &
                   squareup_and_sort_1234_to_2413_complex
   end interface squareup_and_sort_1234_to_2413
!
   interface squareup_and_sort_1234_to_2314
      procedure :: squareup_and_sort_1234_to_2314, &
                   squareup_and_sort_1234_to_2314_complex
   end interface squareup_and_sort_1234_to_2314
!
   interface sort_1234_to_2314
      procedure :: sort_1234_to_2314, &
                   sort_1234_to_2314_complex
   end interface sort_1234_to_2314
!
   interface sort_1234_to_2134
      procedure :: sort_1234_to_2134, &
                   sort_1234_to_2134_complex
   end interface sort_1234_to_2134
!
   interface sort_1234_to_2413
      procedure :: sort_1234_to_2413, &
                   sort_1234_to_2413_complex
   end interface sort_1234_to_2413
!
   interface sort_1234_to_2431
      procedure :: sort_1234_to_2431, &
                   sort_1234_to_2431_complex
   end interface sort_1234_to_2431
!
   interface sort_1234_to_1324
      procedure :: sort_1234_to_1324, &
                   sort_1234_to_1324_complex
   end interface sort_1234_to_1324
!
   interface sort_1234_to_2143
      procedure :: sort_1234_to_2143, &
                   sort_1234_to_2143_complex
   end interface
!
   interface squareup_and_sort_1234_to_1324
      procedure :: squareup_and_sort_1234_to_1324, &
                   squareup_and_sort_1234_to_1324_complex
   end interface squareup_and_sort_1234_to_1324
!
   interface sort_1234_to_1432
      procedure :: sort_1234_to_1432, &
                   sort_1234_to_1432_complex
   end interface sort_1234_to_1432
!
   interface sort_1234_to_4312
      procedure :: sort_1234_to_4312, &
                   sort_1234_to_4312_complex
   end interface sort_1234_to_4312
!
   interface sort_1234_to_1423
      procedure :: sort_1234_to_1423, &
                   sort_1234_to_1423_complex
   end interface sort_1234_to_1423
!
   interface add_1423_to_1234
      procedure :: add_1423_to_1234, &
                   add_1423_to_1234_complex
   end interface add_1423_to_1234
!
   interface add_1432_to_1234
      procedure :: add_1432_to_1234, &
                   add_1432_to_1234_complex
   end interface add_1432_to_1234
!
   interface add_1324_to_1234
      procedure :: add_1324_to_1234, &
                   add_1324_to_1234_complex
   end interface add_1324_to_1234
!
   interface add_1243_to_1234
      procedure :: add_1243_to_1234, &
                   add_1243_to_1234_complex
   end interface add_1243_to_1234
!
   interface add_3421_to_1234
      procedure :: add_3421_to_1234, &
                   add_3421_to_1234_complex
   end interface add_3421_to_1234
!
   interface add_2341_to_1234
      procedure :: add_2341_to_1234, &
                   add_2341_to_1234_complex
   end interface add_2341_to_1234
!
   interface add_2143_to_1234
      procedure :: add_2143_to_1234, &
                   add_2143_to_1234_complex
   end interface add_2143_to_1234
!
   interface add_3214_to_1234
      procedure :: add_3214_to_1234, &
                   add_3214_to_1234_complex
   end interface add_3214_to_1234
!
   interface add_2431_to_1234
      procedure :: add_2431_to_1234, &
                   add_2431_to_1234_complex
   end interface add_2431_to_1234
!
   interface add_4321_to_1234
      procedure :: add_4321_to_1234, &
                   add_4321_to_1234_complex
   end interface add_4321_to_1234
!
   interface add_4123_to_1234
      procedure :: add_4123_to_1234, &
                   add_4123_to_1234_complex
   end interface add_4123_to_1234
!
   interface squareup_and_sort_1234_to_1432
      procedure :: squareup_and_sort_1234_to_1432, &
                   squareup_and_sort_1234_to_1432_complex
   end interface squareup_and_sort_1234_to_1432
!
   interface sort_1234_to_3214
      procedure :: sort_1234_to_3214, &
                   sort_1234_to_3214_complex
   end interface sort_1234_to_3214
!
   interface sort_1234_to_3241
      procedure :: sort_1234_to_3241, &
                   sort_1234_to_3241_complex
   end interface sort_1234_to_3241
!
   interface sort_1234_to_1243
      procedure :: sort_1234_to_1243, &
                   sort_1234_to_1243_complex
   end interface sort_1234_to_1243
!
   interface sort_1234_to_2341
      procedure :: sort_1234_to_2341, &
                   sort_1234_to_2341_complex
   end interface sort_1234_to_2341
!
   interface symmetrize_and_add_to_packed
      procedure :: symmetrize_and_add_to_packed_real, &
                   symmetrize_and_add_4_to_packed
   end interface symmetrize_and_add_to_packed
!
   interface squareup_anti
      procedure :: squareup_anti, &
                   squareup_anti_complex
   end interface squareup_anti
!
   interface squareup
      procedure :: squareup_real, &
                   squareup_complex, &
                   squareup_to_4, &
                   squareup_to_4_complex
   end interface squareup
!
   interface packin
      procedure :: packin_real,      &
                   packin_complex,   &
                   packin_4_real,    &
                   packin_4_complex, &
                   packin_4_from_1324_order_real, &
                   packin_4_from_1324_order_complex
   end interface packin
!
!
contains
!
!     -::- Two-index re-sort and re-sort-add routines -::-
!     ----------------------------------------------------
!
   subroutine sort_12_to_21_r(x_pq, x_qp, dim_p, dim_q)
!!
!!    Sort 12 to 21 real
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_p_q to x_q_p (i.e., 12 to 21).
!!
!!    The unordered array x_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s.
!!    The ordered array x_sp_rq is assumed allocated as dim_s*dim_p x dim_r*dim_q.
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q
!
      real(dp), dimension(dim_p, dim_q), intent(in)  :: x_pq
      real(dp), dimension(dim_q, dim_p), intent(out) :: x_qp
!
      integer :: p, q
!
!$omp parallel do schedule(static) private(p,q)
      do p = 1, dim_p
         do q = 1, dim_q
!
            x_qp(q, p) = x_pq(p, q)
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_12_to_21_r
!
   subroutine sort_12_compound_to_21(x_pq, x_qp, dim_p, dim_q)
!!
!!    Sort 12 compound to 21
!!    Written by Alexander C. Paul, May 2021
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q
!
      real(dp), dimension(dim_p*dim_q), intent(in)  :: x_pq
      real(dp), dimension(dim_q, dim_p), intent(out) :: x_qp
!
      call sort_12_to_21_r(x_pq, x_qp, dim_p, dim_q)
!
   end subroutine sort_12_compound_to_21
!
!
   subroutine sort_12_to_21_c(x_p_q, x_q_p, dim_p, dim_q)
!!
!!    Sort 12 to 21 complex
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Reorders the array x_p_q to x_q_p (i.e., 12 to 21).
!!
!!    The unordered array x_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s.
!!    The ordered array x_sp_rq is assumed allocated as dim_s*dim_p x dim_r*dim_q.
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q
!
      complex(dp), dimension(dim_p, dim_q), intent(in)  :: x_p_q
      complex(dp), dimension(dim_q, dim_p), intent(out) :: x_q_p
!
      integer :: p, q
!
!$omp parallel do schedule(static) private(p,q)
      do p = 1, dim_p
         do q = 1, dim_q
!
            x_q_p(q, p) = x_p_q(p, q)
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_12_to_21_c
!
!
   subroutine add_21_to_12(scalar, x, y_p_q, dim_p, dim_q)
!!
!!    Add 21 to 12
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_p_q(p,q) = y_p_q(p,q) + scalar*x(q,p)
!!
!!    The unordered array y_p_q is assumed allocated as dim_p x dim_q,
!!    and x accordingly.
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q
!
      real(dp), dimension(dim_p, dim_q), intent(inout) :: y_p_q
      real(dp), dimension(dim_q, dim_p), intent(in)    :: x
!
      integer :: p, q
!
!$omp parallel do schedule(static) private(p, q)
      do q = 1, dim_q
         do p = 1, dim_p
!
               y_p_q(p,q) = y_p_q(p,q) + scalar*x(q,p)
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_21_to_12
!
!
   subroutine add_21_to_12_complex(scalar, x, y_p_q, dim_p, dim_q)
!!
!!    Add 21 to 12
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Performs:
!!
!!       y_p_q(p,q) = y_p_q(p,q) + scalar*x(q,p)
!!
!!    The unordered array y_p_q is assumed allocated as dim_p x dim_q,
!!    and x accordingly.
!!
      implicit none
!
      complex(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q
!
      complex(dp), dimension(dim_p, dim_q) :: y_p_q
      complex(dp), dimension(dim_q, dim_p), intent(in) :: x
!
      integer :: p, q
!
!$omp parallel do schedule(static) private(p, q)
      do q = 1, dim_q
         do p = 1, dim_p
!
               y_p_q(p,q) = y_p_q(p,q) + scalar*x(q,p)
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_21_to_12_complex
!
!
   subroutine symmetric_sum_real(x, dim_)
!!
!!    Symmetric sum
!!    Written by Eirik F. Kjønstad, Dec 2017
!!
!!    Performs the action
!!
!!       x(p,q) = x(p,q) + x(q,p)
!!
!!    without making a separate copy of x.
!!
!!    Note: the effect is equal to the routine "add_21_to_12" with scalar = 1,
!!          where a second copy of x is not needed.
!!
      implicit none
!
      integer, intent(in) :: dim_
!
      real(dp), dimension(dim_, dim_), intent(inout) :: x
!
      integer :: p, q
!
!     Overwrite the lower triangular part of the matrix
!
!$omp parallel do private(p, q)
      do q = 1, dim_
         do p = q, dim_
!
            x(p,q) = x(p,q) + x(q,p)
!
         enddo
      enddo
!$omp end parallel do
!
!     Copy the lower triangular part to the upper triangular part
!
!$omp parallel do private(p, q)
      do p = 1, dim_
         do q = p + 1, dim_
!
            x(p,q) = x(q,p)
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine symmetric_sum_real
!
!
   subroutine symmetric_sum_complex(x, dim_)
!!
!!    Symmetric sum
!!    Written by Eirik F. Kjønstad, Dec 2017
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Performs the action
!!
!!       x(p,q) = x(p,q) + x(q,p)
!!
!!    without making a separate copy of x.
!!
!!    Note: the effect is equal to the routine "add_21_to_12" with scalar = 1,
!!          where a second copy of x is not needed.
!!
      implicit none
!
      integer, intent(in) :: dim_
!
      complex(dp), dimension(dim_, dim_), intent(inout) :: x
!
      integer :: p, q
!
!     Overwrite the lower triangular part of the matrix
!
!$omp parallel do private(p, q)
      do q = 1, dim_
         do p = q, dim_
!
            x(p,q) = x(p,q) + x(q,p)
!
         enddo
      enddo
!$omp end parallel do
!
!     Copy the lower triangular part to the upper triangular part
!
!$omp parallel do private(p, q)
      do p = 1, dim_
         do q = p + 1, dim_
!
            x(p,q) = x(q,p)
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine symmetric_sum_complex
!
!
   subroutine symmetric_sum_4(x, dim_)
!!
!!    Symmetric sum
!!    Written by Andreas Skeidsvoll, Sep 2019
!!
!!    Performs the action
!!
!!       x(p,q) = x(p,q) + x(q,p)
!!
!!
      implicit none
!
      integer, intent(in) :: dim_
!
      real(dp), dimension(:,:,:,:), intent(inout) :: x
!
      call symmetric_sum_real(x, dim_)
!
   end subroutine symmetric_sum_4
!
!
   subroutine symmetric_sum_4_complex(x, dim_)
!!
!!    Symmetric sum
!!    Written by Andreas Skeidsvoll, Sep 2019
!!
!!    Performs the action
!!
!!       x(p,q) = x(p,q) + x(q,p)
!!
!!
      implicit none
!
      integer, intent(in) :: dim_
!
      complex(dp), dimension(:,:,:,:), intent(inout) :: x
!
      call symmetric_sum_complex(x, dim_)
!
   end subroutine symmetric_sum_4_complex
!
!
   subroutine symmetrize_12_and_34(x, dim_p, dim_r)
!!
!!    Symmetric sum
!!    Written by Rolf H. Myhre, May 2019
!!
!!    Performs the action
!!
!!       x(p,q,r,s) = x(p,q,r,s) + x(q,p,s,r)
!!
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_r
!
      real(dp), dimension(dim_p, dim_p, dim_r, dim_r), intent(inout) :: x
!
      integer :: p, q, r, s, p_lim
!
!     Lower triangles first
!
!$omp parallel do private(p, q, r, s, p_lim)
      do s = 1, dim_r
         do r = s, dim_r
            do q = 1, dim_p
!
               if(r .ne. s) then
                  p_lim = 1
               else
                  p_lim = q
               endif
!
               do p = p_lim, dim_p
!
                  x(p,q,r,s) = x(p,q,r,s) + x(q,p,s,r)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
!  And the upper triangles
!
!$omp parallel do private(p, q, r, s, p_lim)
      do s = 1, dim_r
         do r = 1, s
            do q = 1, dim_p
!
               if(r .ne. s) then
                  p_lim = dim_p
               else
                  p_lim = q-1
               endif
!
               do p = 1, p_lim
!
                  x(p,q,r,s) = x(q,p,s,r)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine symmetrize_12_and_34
!
!
!     -::- Three-index re-sort and re-sort-add routines -::-
!     ------------------------------------------------------
!
   subroutine sort_123_to_312(x_pqr, x_rpq, dim_p, dim_q, dim_r)
!!
!!    Sort 123 to 312
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       x_rpq(r,p,q) = x_pqr(p,q,r)
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r
!
      real(dp), dimension(dim_p, dim_q, dim_r), intent(in)  :: x_pqr
      real(dp), dimension(dim_r, dim_p, dim_q), intent(out) :: x_rpq
!
      integer :: r, q, p
!
!$omp parallel do schedule(static) private(r,q,p)
      do q = 1, dim_q
         do p = 1, dim_p
            do r = 1, dim_r
!
               x_rpq(r,p,q) = x_pqr(p,q,r)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_123_to_312
!
!
   subroutine sort_123_to_312_complex(x_pqr, x_rpq, dim_p, dim_q, dim_r)
!!
!!    Sort 123 to 312
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Performs:
!!
!!       x_rpq(r,p,q) = x_pqr(p,q,r)
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r
!
      complex(dp), dimension(dim_p, dim_q, dim_r), intent(in)  :: x_pqr
      complex(dp), dimension(dim_r, dim_p, dim_q), intent(out) :: x_rpq
!
      integer :: r, q, p
!
!$omp parallel do schedule(static) private(r,q,p)
      do q = 1, dim_q
         do p = 1, dim_p
            do r = 1, dim_r
!
               x_rpq(r,p,q) = x_pqr(p,q,r)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_123_to_312_complex
!
!
   subroutine sort_123_to_231(x_pqr, x_qrp, dim_p, dim_q, dim_r)
!!
!!    Sort 123 to 312
!!    Written by Rolf H. Myhre, April 2019
!!
!!    Performs:
!!
!!       x_qrp(q,r,p) = x_pqr(p,q,r)
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r
!
      real(dp), dimension(dim_p, dim_q, dim_r), intent(in)  :: x_pqr
      real(dp), dimension(dim_q, dim_r, dim_p), intent(out) :: x_qrp
!
      integer :: r, q, p
!
!$omp parallel do schedule(static) private(r,q,p)
      do p = 1, dim_p
         do r = 1, dim_r
            do q = 1, dim_q
!
               x_qrp(q,r,p) = x_pqr(p,q,r)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_123_to_231
!
!
   subroutine add_312_to_123(scalar, x_rpq, x_pqr, dim_p, dim_q, dim_r)
!!
!!    Sort 123 to 312 and add
!!    Written by Andreas Skeidsvoll, Sep 2019
!!
!!    Performs:
!!
!!       x_pqr(p,q,r) = x_pqr(p,q,r) + scalar*x_rpq(r,p,q)
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r
!
      real(dp), intent(in) :: scalar
!
      real(dp), dimension(dim_r, dim_p, dim_q), intent(in)       :: x_rpq
      real(dp), dimension(dim_p, dim_q, dim_r), intent(inout)    :: x_pqr
!
      integer :: r, q, p
!
!$omp parallel do schedule(static) private(r,q,p)
      do q = 1, dim_q
         do p = 1, dim_p
            do r = 1, dim_r
!
               x_pqr(p,q,r) = x_pqr(p,q,r) +  scalar*x_rpq(r,p,q)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_312_to_123
!
!
   subroutine sort_123_to_132(x_pqr, x_prq, dim_p, dim_q, dim_r)
!!
!!    Sort 123 to 132
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       x_prq(p,r,q) = x_pqr(p,q,r)
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r
!
      real(dp), dimension(dim_p, dim_q, dim_r), intent(in)  :: x_pqr
      real(dp), dimension(dim_p, dim_r, dim_q), intent(out) :: x_prq
!
      integer :: r, q, p
!
!$omp parallel do schedule(static) private(r,q,p)
      do q = 1, dim_q
         do r = 1, dim_r
            do p = 1, dim_p
!
               x_prq(p,r,q) = x_pqr(p,q,r)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_123_to_132
!
!
   subroutine sort_123_to_132_complex(x_pqr, x_prq, dim_p, dim_q, dim_r)
!!
!!    Sort 123 to 132
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Performs:
!!
!!       x_prq(p,r,q) = x_pqr(p,q,r)
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r
!
      complex(dp), dimension(dim_p, dim_q, dim_r), intent(in)  :: x_pqr
      complex(dp), dimension(dim_p, dim_r, dim_q), intent(out) :: x_prq
!
      integer :: r, q, p
!
!$omp parallel do schedule(static) private(r,q,p)
      do q = 1, dim_q
         do r = 1, dim_r
            do p = 1, dim_p
!
               x_prq(p,r,q) = x_pqr(p,q,r)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_123_to_132_complex
!
!
   subroutine sort_123_to_213(x_pqr, x_qpr, dim_p, dim_q, dim_r)
!!
!!    Sort 123 to 132
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       x_qpr(q,p,r) = x_pqr(p,q,r)
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r
!
      real(dp), dimension(dim_p, dim_q, dim_r), intent(in)  :: x_pqr
      real(dp), dimension(dim_q, dim_p, dim_r), intent(out) :: x_qpr
!
      integer :: p, q, r
!
!$omp parallel do schedule(static) private(p, q, r)
      do r = 1, dim_r
         do p = 1, dim_p
            do q = 1, dim_q
!
               x_qpr(q,p,r) = x_pqr(p,q,r)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_123_to_213
!
!
   subroutine sort_123_to_213_complex(x_pqr, x_qpr, dim_p, dim_q, dim_r)
!!
!!    Sort 123 to 132
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Performs:
!!
!!       x_qpr(q,p,r) = x_pqr(p,q,r)
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r
!
      complex(dp), dimension(dim_p, dim_q, dim_r), intent(in)  :: x_pqr
      complex(dp), dimension(dim_q, dim_p, dim_r), intent(out) :: x_qpr
!
      integer :: p, q, r
!
!$omp parallel do schedule(static) private(p, q, r)
      do r = 1, dim_r
         do p = 1, dim_p
            do q = 1, dim_q
!
               x_qpr(q,p,r) = x_pqr(p,q,r)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_123_to_213_complex
!
!
   subroutine add_213_to_123(scalar, x_qpr, y_pqr, dim_p, dim_q, dim_r)
!!
!!    Add 213 to 123
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pqr(pqr) = y_pqr(pqr) + scalar*x(qpr)
!!
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r
!
      real(dp), dimension(dim_p, dim_q, dim_r) :: y_pqr
      real(dp), dimension(dim_q, dim_p, dim_r), intent(in) :: x_qpr
!
      integer :: r, q, p
!
!$omp parallel do schedule(static) private(r,q,p)
      do r = 1, dim_r
         do q = 1, dim_q
            do p = 1, dim_p
!
               y_pqr(p,q,r) = y_pqr(p,q,r) + scalar*x_qpr(q,p,r)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_213_to_123
!
!
   subroutine add_321_to_123(scalar, x_rqp, y_pqr, dim_p, dim_q, dim_r)
!!
!!    Add 321 to 123
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pqr(pqr) = y_pqr(pqr) + scalar*x(rqp)
!!
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r
!
      real(dp), dimension(dim_p, dim_q, dim_r) :: y_pqr
      real(dp), dimension(dim_q, dim_p, dim_r), intent(in) :: x_rqp
!
      integer :: r, q, p
!
!$omp parallel do schedule(static) private(r,q,p)
      do r = 1, dim_r
         do q = 1, dim_q
            do p = 1, dim_p
!
               y_pqr(p,q,r) = y_pqr(p,q,r) + scalar*x_rqp(r,q,p)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_321_to_123
!
!
   subroutine add_132_to_123(scalar, x_prq, y_pqr, dim_p, dim_q, dim_r)
!!
!!    Add 132 to 123
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pqr(p,q,r) = y_pqr(p,q,r) + scalar * x_prq(p,r,q)
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r
!
      real(dp), dimension(dim_p, dim_q, dim_r) :: y_pqr
      real(dp), dimension(dim_p, dim_r, dim_q), intent(in) :: x_prq
!
      integer :: p, q, r
!
!$omp parallel do schedule(static) private(r,q,p)
      do r = 1, dim_r
         do q = 1, dim_q
            do p = 1, dim_p
!
               y_pqr(p,q,r) = y_pqr(p,q,r) + scalar*x_prq(p,r,q)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_132_to_123
!
!
   subroutine add_132_to_123_complex(scalar, x_prq, y_pqr, dim_p, dim_q, dim_r)
!!
!!    Add 132 to 123
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Performs:
!!
!!       y_pqr(p,q,r) = y_pqr(p,q,r) + scalar * x_prq(p,r,q)
!!
      implicit none
!
      complex(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r
!
      complex(dp), dimension(dim_p, dim_q, dim_r) :: y_pqr
      complex(dp), dimension(dim_p, dim_r, dim_q), intent(in) :: x_prq
!
      integer :: p, q, r
!
!$omp parallel do schedule(static) private(r,q,p)
      do r = 1, dim_r
         do q = 1, dim_q
            do p = 1, dim_p
!
               y_pqr(p,q,r) = y_pqr(p,q,r) + scalar*x_prq(p,r,q)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_132_to_123_complex
!
!
!     -::- Four-index re-sort and re-sort-add routines -::-
!     -----------------------------------------------------
!
   subroutine add_3124_to_1234(scalar, x_rpqs, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 3124 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x_rpqs(r,p,q,s)
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      real(dp), dimension(dim_r, dim_p, dim_q, dim_s), intent(in) :: x_rpqs
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_rpqs(r,p,q,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_3124_to_1234
!
!
   subroutine sort_1234_to_3412(x_pqrs, x_rspq, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 3412
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pqrs to x_rspq (i.e., 1234 to 3412).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_r, dim_s, dim_p, dim_q), intent(out) :: x_rspq
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do q = 1, dim_q
         do p = 1, dim_p
            do s = 1, dim_s
               do r = 1, dim_r
!
                  x_rspq(r, s, p, q) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_3412
!
!
   subroutine sort_1234_to_3412_complex(x_pqrs, x_rspq, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 3412
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Reorders the array x_pqrs to x_rspq (i.e., 1234 to 3412).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      complex(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      complex(dp), dimension(dim_r, dim_s, dim_p, dim_q), intent(out) :: x_rspq
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do q = 1, dim_q
         do p = 1, dim_p
            do s = 1, dim_s
               do r = 1, dim_r
!
                  x_rspq(r, s, p, q) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_3412_complex
!
!
   subroutine sort_1234_to_4132(x_pqrs, x_sprq, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 4132
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pqrs to x_sprq (i.e., 1234 to 4132).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_s, dim_p, dim_r, dim_q), intent(out) :: x_sprq
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do q = 1, dim_q
         do r = 1, dim_r
            do p = 1, dim_p
               do s = 1, dim_s
!
                  x_sprq(s,p,r,q) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_4132
!
!
   subroutine sort_1234_to_4132_complex(x_pqrs, x_sprq, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 4132
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Reorders the array x_pqrs to x_sprq (i.e., 1234 to 4132).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      complex(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      complex(dp), dimension(dim_s, dim_p, dim_r, dim_q), intent(out) :: x_sprq
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do q = 1, dim_q
         do r = 1, dim_r
            do p = 1, dim_p
               do s = 1, dim_s
!
                  x_sprq(s,p,r,q) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_4132_complex
!
!
   subroutine sort_1234_to_4123(x_pqrs, x_spqr, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 4132
!!    Written by Eirik F. Kjønstad, Rolf H. Myhre and Andreas Skeidsvoll, 2018
!!
!!    Reorders the array x_pqrs to x_spqr (i.e., 1234 to 4123).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_s, dim_p, dim_q, dim_r), intent(out) :: x_spqr
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do q = 1, dim_q
         do r = 1, dim_r
            do p = 1, dim_p
               do s = 1, dim_s
!
                  x_spqr(s,p,q,r) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_4123
!
!
   subroutine sort_1234_to_4123_complex(x_pqrs, x_spqr, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 4132
!!    Written by Eirik F. Kjønstad, Rolf H. Myhre and Andreas Skeidsvoll, 2018
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Reorders the array x_pqrs to x_spqr (i.e., 1234 to 4123).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      complex(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      complex(dp), dimension(dim_s, dim_p, dim_q, dim_r), intent(out) :: x_spqr
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do q = 1, dim_q
         do r = 1, dim_r
            do p = 1, dim_p
               do s = 1, dim_s
!
                  x_spqr(s,p,q,r) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_4123_complex
!
!
   subroutine sort_1234_to_3124(x_pqrs, x_rpqs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 3124
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pqrs to x_rpqs (i.e., 1234 to 3124).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_r, dim_p, dim_q, dim_s), intent(out) :: x_rpqs
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do q = 1, dim_q
            do p = 1, dim_p
               do r = 1, dim_r
!
                  x_rpqs(r,p,q,s) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_3124
!
!
   subroutine sort_1234_to_3142(x_pqrs, x_rpsq, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 3142
!!    Written by Eirik F. Kjønstad, Rolf H. Myhre and
!!    Sarai D. Folkestad, 2018
!!
!!    Reorders the array x_pqrs to x_rpsq (i.e., 1234 to 3142).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_r, dim_p, dim_s, dim_q), intent(out) :: x_rpsq
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do q = 1, dim_q
         do s = 1, dim_s
            do p = 1, dim_p
               do r = 1, dim_r
!
                  x_rpsq(r,p,s,q) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_3142
!
!
   subroutine sort_1234_to_3142_complex(x_pqrs, x_rpsq, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 3142
!!    Written by Eirik F. Kjønstad, Rolf H. Myhre and
!!    Sarai D. Folkestad, 2018
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Reorders the array x_pqrs to x_rpsq (i.e., 1234 to 3142).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      complex(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      complex(dp), dimension(dim_r, dim_p, dim_s, dim_q), intent(out) :: x_rpsq
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do q = 1, dim_q
         do s = 1, dim_s
            do p = 1, dim_p
               do r = 1, dim_r
!
                  x_rpsq(r,p,s,q) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_3142_complex
!
!
   subroutine squareup_and_sort_1234_to_2413(x_pqrs, x_qspr, dim_p, dim_q, dim_r, dim_s)
!!
!!    Squareup and sort 1234 to 2413
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Square up and reorder the array x_pqrs to x_qspr (i.e., 1234 to 2413).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(((dim_p*dim_q+1)*dim_r*dim_s)/2), intent(in) :: x_pqrs
      real(dp), dimension(dim_q, dim_s, dim_p, dim_r), intent(out)     :: x_qspr
!
      integer :: p, q, r, s, rs, pq, pqrs
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,pqrs)
      do r = 1, dim_r
         do p = 1, dim_p
            do s = 1, dim_s
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
!
                  pq = dim_p*(q-1) + p
                  pqrs = (max(pq,rs)*(max(pq,rs)-3)/2) + pq + rs
!
                  x_qspr(q,s,p,r) = x_pqrs(pqrs)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine squareup_and_sort_1234_to_2413
!
!
   subroutine squareup_and_sort_1234_to_2413_complex(x_pqrs, x_qspr, dim_p, dim_q, dim_r, dim_s)
!!
!!    Squareup and sort 1234 to 2413
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Square up and reorder the array x_pqrs to x_qspr (i.e., 1234 to 2413).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      complex(dp), dimension(((dim_p*dim_q+1)*dim_r*dim_s)/2), intent(in) :: x_pqrs
      complex(dp), dimension(dim_q, dim_s, dim_p, dim_r), intent(out)     :: x_qspr
!
      integer :: p, q, r, s, rs, pq, pqrs
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,pqrs)
      do r = 1, dim_r
         do p = 1, dim_p
            do s = 1, dim_s
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
!
                  pq = dim_p*(q-1) + p
                  pqrs = (max(pq,rs)*(max(pq,rs)-3)/2) + pq + rs
!
                  x_qspr(q,s,p,r) = x_pqrs(pqrs)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine squareup_and_sort_1234_to_2413_complex
!
!
   subroutine squareup_and_sort_1234_to_2314(x_pqrs, x_qrps, dim_p, dim_q, dim_r, dim_s)
!!
!!    Squareup and sort 1234 to 2314
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Square up and reorder the array x_pqrs to x_qspr (i.e., 1234 to 2314).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(((dim_p*dim_q+1)*dim_r*dim_s)/2), intent(in) :: x_pqrs
      real(dp), dimension(dim_q, dim_r, dim_p, dim_s), intent(out)     :: x_qrps
!
      integer :: p, q, r, s, rs, pq, pqrs
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,pqrs)
      do s = 1, dim_s
         do p = 1, dim_p
            do r = 1, dim_r
!
               rs = dim_r*(s-1) + r
!
               do q = 1, dim_q
!
                  pq = dim_p*(q-1) + p
                  pqrs = (max(pq,rs)*(max(pq,rs)-3)/2) + pq + rs
!
                  x_qrps(q,r,p,s) = x_pqrs(pqrs)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine squareup_and_sort_1234_to_2314
!
!
   subroutine squareup_and_sort_1234_to_2314_complex(x_pqrs, x_qrps, dim_p, dim_q, dim_r, dim_s)
!!
!!    Squareup and sort 1234 to 2314
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Square up and reorder the array x_pqrs to x_qspr (i.e., 1234 to 2413).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      complex(dp), dimension(((dim_p*dim_q+1)*dim_r*dim_s)/2), intent(in) :: x_pqrs
      complex(dp), dimension(dim_q, dim_r, dim_p, dim_s), intent(out)     :: x_qrps
!
      integer :: p, q, r, s, rs, pq, pqrs
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,pqrs)
      do s = 1, dim_s
         do p = 1, dim_p
            do r = 1, dim_r
!
               rs = dim_r*(s-1) + r
!
               do q = 1, dim_q
!
                  pq = dim_p*(q-1) + p
                  pqrs = (max(pq,rs)*(max(pq,rs)-3)/2) + pq + rs
!
                  x_qrps(q,r,p,s) = x_pqrs(pqrs)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine squareup_and_sort_1234_to_2314_complex
!
!
   subroutine sort_1234_to_2314(x_pqrs, x_qrps, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 4132
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pqrs to x_qrps (i.e., 1234 to 2314).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_q, dim_r, dim_p, dim_s), intent(out) :: x_qrps
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do p = 1, dim_p
            do r = 1, dim_r
               do q = 1, dim_q
!
                  x_qrps(q,r,p,s) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_2314
!
!
   subroutine sort_1234_to_2314_complex(x_pqrs, x_qrps, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 4132
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Reorders the array x_pqrs to x_qrps (i.e., 1234 to 2314).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      complex(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      complex(dp), dimension(dim_q, dim_r, dim_p, dim_s), intent(out) :: x_qrps
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do p = 1, dim_p
            do r = 1, dim_r
               do q = 1, dim_q
!
                  x_qrps(q,r,p,s) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_2314_complex
!
!
   subroutine sort_1234_to_2134(x_pqrs, x_qprs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 2134
!!    Written by Sarai D. Folkestad,
!!    Eirik F. Kjønstad and Rolf H. Myhre, 2018
!!
!!    Reorders the array x_pqrs to x_qprs (i.e., 1234 to 2134).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_q, dim_p, dim_r, dim_s), intent(out) :: x_qprs
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do p = 1, dim_p
               do q = 1, dim_q
!
                  x_qprs(q, p, r, s) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_2134
!
!
   subroutine sort_1234_to_2134_complex(x_pqrs, x_qprs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 2134
!!    Written by Sarai D. Folkestad,
!!    Eirik F. Kjønstad and Rolf H. Myhre, 2018
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Reorders the array x_pqrs to x_qprs (i.e., 1234 to 2134).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      complex(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      complex(dp), dimension(dim_q, dim_p, dim_r, dim_s), intent(out) :: x_qprs
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do p = 1, dim_p
               do q = 1, dim_q
!
                  x_qprs(q, p, r, s) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_2134_complex
!
!
   subroutine sort_1234_to_2143(x_pqrs, x_qpsr, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 2143
!!    Written by Sarai D. Folkestad,
!!    Eirik F. Kjønstad and Rolf H. Myhre, 2018
!!
!!    Reorders the array x_pqrs to x_qpsr (i.e., 1234 to 2143).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_q, dim_p, dim_s, dim_r), intent(out) :: x_qpsr
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do r = 1, dim_r
         do s = 1, dim_s
            do p = 1, dim_p
               do q = 1, dim_q
!
                  x_qpsr(q, p, s, r) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_2143
!
!
   subroutine sort_1234_to_2143_complex(x_pqrs, x_qpsr, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 2143
!!    Written by Sarai D. Folkestad,
!!    Eirik F. Kjønstad and Rolf H. Myhre, 2018
!!
!!    Reorders the array x_pqrs to x_qpsr (i.e., 1234 to 2143).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      complex(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      complex(dp), dimension(dim_q, dim_p, dim_s, dim_r), intent(out) :: x_qpsr
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do r = 1, dim_r
         do s = 1, dim_s
            do p = 1, dim_p
               do q = 1, dim_q
!
                  x_qpsr(q, p, s, r) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_2143_complex
!
!
   subroutine sort_1234_to_2413(x_pqrs, x_qspr, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 2413
!!    Written by Sarai D. Folkestad,
!!    Eirik F. Kjønstad and Rolf H. Myhre, 2018
!!
!!    Reorders the array x_pqrs to x_qspr (i.e., 1234 to 2413).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_q, dim_s, dim_p, dim_r), intent(out) :: x_qspr
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do r = 1, dim_r
         do p = 1, dim_p
            do s = 1, dim_s
               do q = 1, dim_q
!
                  x_qspr(q, s, p, r) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_2413
!
!
   subroutine sort_1234_to_2413_complex(x_pqrs, x_qspr, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 2413
!!    Written by Sarai D. Folkestad,
!!    Eirik F. Kjønstad and Rolf H. Myhre, 2018
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Reorders the array x_pqrs to x_qspr (i.e., 1234 to 2413).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      complex(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      complex(dp), dimension(dim_q, dim_s, dim_p, dim_r), intent(out) :: x_qspr
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do r = 1, dim_r
         do p = 1, dim_p
            do s = 1, dim_s
               do q = 1, dim_q
!
                  x_qspr(q, s, p, r) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_2413_complex
!
!
   subroutine sort_1234_to_2431(x_pqrs, x_qsrp, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 2431
!!    Written by Sarai D. Folkestad,
!!    Eirik F. Kjønstad and Rolf H. Myhre, 2018
!!
!!    Reorders the array x_pqrs to x_qsrp (i.e., 1234 to 2431).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in) :: x_pqrs
      real(dp), dimension(dim_q, dim_s, dim_r, dim_p)             :: x_qsrp
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(p,q,r,s)
      do p = 1, dim_p
         do r = 1, dim_r
            do s = 1, dim_s
               do q = 1, dim_q
!
                  x_qsrp(q,s,r,p) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_2431
!
!
   subroutine sort_1234_to_2431_complex(x_pqrs, x_qsrp, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 2431
!!    Written by Sarai D. Folkestad,
!!    Eirik F. Kjønstad and Rolf H. Myhre, 2018
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Reorders the array x_pqrs to x_qsrp (i.e., 1234 to 2431).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      complex(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in) :: x_pqrs
      complex(dp), dimension(dim_q, dim_s, dim_r, dim_p)             :: x_qsrp
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(p,q,r,s)
      do p = 1, dim_p
         do r = 1, dim_r
            do s = 1, dim_s
               do q = 1, dim_q
!
                  x_qsrp(q,s,r,p) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_2431_complex
!
!
   subroutine sort_1234_to_3421(x_pqrs, x_rsqp, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 3421
!!    Written by Sarai D. Folkestad,
!!    Eirik F. Kjønstad and Rolf H. Myhre, 2018
!!
!!    Reorders the array x_pqrs to x_rsqp (i.e., 1234 to 3421).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in) :: x_pqrs
      real(dp), dimension(dim_r, dim_s, dim_q, dim_p)             :: x_rsqp
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(p,q,r,s)
      do p = 1, dim_p
         do q = 1, dim_q
            do s = 1, dim_s
               do r = 1, dim_r
!
                  x_rsqp(r,s,q,p) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_3421
!
!
   subroutine sort_1234_to_1324(x_pqrs, x_prqs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 1324
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pq_rs to x_pr_qs (i.e., 1234 to 1324).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_p, dim_r, dim_q, dim_s), intent(out) :: x_prqs
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do q = 1, dim_q
            do r = 1, dim_r
               do p = 1, dim_p
!
                  x_prqs(p, r, q, s) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_1324
!
!
   subroutine sort_1234_to_1324_complex(x_pqrs, x_prqs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 1324
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Reorders the array x_pq_rs to x_pr_qs (i.e., 1234 to 1324).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      complex(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      complex(dp), dimension(dim_p, dim_r, dim_q, dim_s), intent(out) :: x_prqs
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do q = 1, dim_q
            do r = 1, dim_r
               do p = 1, dim_p
!
                  x_prqs(p, r, q, s) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_1324_complex
!
!
   subroutine squareup_and_sort_1234_to_1324(x_pqrs, x_prqs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Squareup and sort 1234 to 1324
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pq_rs to x_pr_qs (i.e., 1234 to 1324).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(((dim_p*dim_q+1)*(dim_r*dim_s)/2)), intent(in) :: x_pqrs
      real(dp), dimension(dim_p, dim_r, dim_q, dim_s), intent(out)       :: x_prqs
!
      integer :: p, q, r, s, pq, rs, pqrs
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,pqrs)
      do s = 1, dim_s
         do q = 1, dim_q
            do r = 1, dim_r
!
               rs = dim_r*(s-1) + r
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  pqrs = (max(pq,rs)*(max(pq,rs)-3)/2) + pq + rs
!
                  x_prqs(p,r,q,s) = x_pqrs(pqrs)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine squareup_and_sort_1234_to_1324
!
!
   subroutine squareup_and_sort_1234_to_1324_complex(x_pqrs, x_prqs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Squareup and sort 1234 to 1324
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Reorders the array x_pq_rs to x_pr_qs (i.e., 1234 to 1324).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      complex(dp), dimension(((dim_p*dim_q+1)*(dim_r*dim_s)/2)), intent(in) :: x_pqrs
      complex(dp), dimension(dim_p, dim_r, dim_q, dim_s), intent(out)       :: x_prqs
!
      integer :: p, q, r, s, pq, rs, pqrs
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,pqrs)
      do s = 1, dim_s
         do q = 1, dim_q
            do r = 1, dim_r
!
               rs = dim_r*(s-1) + r
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  pqrs = (max(pq,rs)*(max(pq,rs)-3)/2) + pq + rs
!
                  x_prqs(p,r,q,s) = x_pqrs(pqrs)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine squareup_and_sort_1234_to_1324_complex
!
!
   subroutine sort_1234_to_2341(x_pqrs, x_qrsp, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 2341
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pqrs to x_qrsp (i.e., 1234 to 2341).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_q, dim_r, dim_s, dim_p), intent(out) :: x_qrsp
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  x_qrsp(q,r,s,p) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_2341
!
!
   subroutine sort_1234_to_2341_complex(x_pqrs, x_qrsp, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 2341
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pqrs to x_qrsp (i.e., 1234 to 2341).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      complex(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      complex(dp), dimension(dim_q, dim_r, dim_s, dim_p), intent(out) :: x_qrsp
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  x_qrsp(q,r,s,p) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_2341_complex
!
!
   subroutine sort_1234_to_1342(x_pqrs, x_prsq, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 1342
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pqrs to x_prsq (i.e., 1234 to 1342).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_p, dim_r, dim_s, dim_q), intent(out) :: x_prsq
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do q = 1, dim_q
         do s = 1, dim_s
            do r = 1, dim_r
               do p = 1, dim_p
!
                  x_prsq(p,r,s,q) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_1342
!
!
   subroutine sort_1234_to_1432(x_pqrs, x_psrq, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 1432
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pqrs to x_psrq (i.e., 1234 to 1432).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_p, dim_s, dim_r, dim_q), intent(out) :: x_psrq
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do q = 1, dim_q
         do r = 1, dim_r
            do s = 1, dim_s
               do p = 1, dim_p
!
                  x_psrq(p,s,r,q) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_1432
!
!
   subroutine sort_1234_to_1432_complex(x_pqrs, x_psrq, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 1432
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Reorders the array x_pqrs to x_psrq (i.e., 1234 to 1432).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      complex(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      complex(dp), dimension(dim_p, dim_s, dim_r, dim_q), intent(out) :: x_psrq
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do q = 1, dim_q
         do r = 1, dim_r
            do s = 1, dim_s
               do p = 1, dim_p
!
                  x_psrq(p,s,r,q) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_1432_complex
!
!
   subroutine sort_1234_to_4312(x_pqrs, x_srpq, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 4312
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pqrs to x_srpq (i.e., 1234 to 4312).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_s, dim_r, dim_p, dim_q), intent(out) :: x_srpq
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(p,q,r,s)
      do q = 1, dim_q
         do p = 1, dim_p
            do r = 1, dim_r
               do s = 1, dim_s
!
                  x_srpq(s,r,p,q) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_4312
!
!
   subroutine sort_1234_to_4312_complex(x_pqrs, x_srpq, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 4312
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Reorders the array x_pqrs to x_srpq (i.e., 1234 to 4312).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      complex(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      complex(dp), dimension(dim_s, dim_r, dim_p, dim_q), intent(out) :: x_srpq
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(p,q,r,s)
      do q = 1, dim_q
         do p = 1, dim_p
            do r = 1, dim_r
               do s = 1, dim_s
!
                  x_srpq(s,r,p,q) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_4312_complex
!
!
   subroutine sort_1234_to_1423(x_pqrs, x_psqr, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 1423
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pqrs to x_psqr (i.e., 1234 to 1423).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_p, dim_s, dim_q, dim_r), intent(out) :: x_psqr
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do r = 1, dim_r
         do q = 1, dim_q
            do s = 1, dim_s
               do p = 1, dim_p
!
                  x_psqr(p,s,q,r) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_1423
!
!
   subroutine sort_1234_to_1423_complex(x_pqrs, x_psqr, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 1423
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Reorders the array x_pqrs to x_psqr (i.e., 1234 to 1423).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      complex(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      complex(dp), dimension(dim_p, dim_s, dim_q, dim_r), intent(out) :: x_psqr
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do r = 1, dim_r
         do q = 1, dim_q
            do s = 1, dim_s
               do p = 1, dim_p
!
                  x_psqr(p,s,q,r) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_1423_complex
!
!
   subroutine add_1423_to_1234(scalar, x_psqr, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 1423 to 1234
!!    Written by Eirik F. Kjønstad, Rolf H. Myhre and Andreas Skeidsvoll, 2018
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x_psqr(p,s,q,r)
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      real(dp), dimension(dim_p, dim_s, dim_q, dim_r), intent(in)    :: x_psqr
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_psqr(p,s,q,r)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_1423_to_1234
!
!
   subroutine add_1423_to_1234_complex(scalar, x_psqr, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 1423 to 1234
!!    Written by Eirik F. Kjønstad, Rolf H. Myhre and Andreas Skeidsvoll, 2018
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x_psqr(p,s,q,r)
!!
      implicit none
!
      complex(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      complex(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      complex(dp), dimension(dim_p, dim_s, dim_q, dim_r), intent(in)    :: x_psqr
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_psqr(p,s,q,r)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_1423_to_1234_complex
!
!
   subroutine add_1432_to_1234(scalar, x_psrq, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 1432 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x(p,s,r,q)
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      real(dp), dimension(dim_p, dim_s, dim_r, dim_q), intent(in)    :: x_psrq
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_psrq(p,s,r,q)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_1432_to_1234
!
!
   subroutine add_1432_to_1234_complex(scalar, x_psrq, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 1432 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x(p,s,r,q)
!!
      implicit none
!
      complex(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      complex(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      complex(dp), dimension(dim_p, dim_s, dim_r, dim_q), intent(in)    :: x_psrq
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_psrq(p,s,r,q)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_1432_to_1234_complex
!
!
   subroutine add_1342_to_1234(scalar, x_prsq, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 1342 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x_prsq(p,r,s,q)
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      real(dp), dimension(dim_p, dim_r, dim_s, dim_q), intent(in)    :: x_prsq
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_prsq(p,r,s,q)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_1342_to_1234
!
!
   subroutine add_1324_to_1234(scalar, x_prqs, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 1342 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Andreas Skeidsvoll, Apr 2019: Created by modifying add_1342_to_1234
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x_prqs(p,r,q,s)
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      real(dp), dimension(dim_p, dim_r, dim_q, dim_s), intent(in)    :: x_prqs
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_prqs(p,r,q,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_1324_to_1234
!
!
   subroutine add_1324_to_1234_complex(scalar, x_prqs, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 1342 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!    Modified by Andreas Skeidsvoll, Apr 2019: Created by modifying add_1342_to_1234
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x_prqs(p,r,q,s)
!!
      implicit none
!
      complex(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      complex(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      complex(dp), dimension(dim_p, dim_r, dim_q, dim_s), intent(in)    :: x_prqs
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_prqs(p,r,q,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_1324_to_1234_complex
!
!
   subroutine add_1243_to_1234(scalar, x_pqsr, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 1243 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x_pqsr(p,q,s,r)
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      real(dp), dimension(dim_p, dim_q, dim_s, dim_r), intent(in)    :: x_pqsr
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_pqsr(p,q,s,r)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_1243_to_1234
!
!
   subroutine add_1243_to_1234_complex(scalar, x_pqsr, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 1243 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x_pqsr(p,q,s,r)
!!
      implicit none
!
      complex(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      complex(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      complex(dp), dimension(dim_p, dim_q, dim_s, dim_r), intent(in)    :: x_pqsr
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_pqsr(p,q,s,r)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_1243_to_1234_complex
!
!
   subroutine add_3412_to_1234(scalar, x_rspq, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 3412 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x(r,s,p,q)
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      real(dp), dimension(dim_r, dim_s, dim_p, dim_q), intent(in)    :: x_rspq
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(p,q,r,s)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_rspq(r,s,p,q)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_3412_to_1234
!
!
   subroutine add_3421_to_1234(scalar, x_rsqp, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 3421 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x_rsqp(r,s,q,p)
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      real(dp), dimension(dim_r, dim_s, dim_q, dim_p), intent(in)    :: x_rsqp
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(p,q,r,s)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_rsqp(r,s,q,p)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_3421_to_1234
!
!
   subroutine add_3421_to_1234_complex(scalar, x_rsqp, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 3421 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x_rsqp(r,s,q,p)
!!
      implicit none
!
      complex(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      complex(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      complex(dp), dimension(dim_r, dim_s, dim_q, dim_p), intent(in)    :: x_rsqp
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(p,q,r,s)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_rsqp(r,s,q,p)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_3421_to_1234_complex
!
!
   subroutine add_2341_to_1234(scalar, x_qrsp, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 2341 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x_qrsp(q,r,s,p)
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      real(dp), dimension(dim_q, dim_r, dim_s, dim_p), intent(in)    :: x_qrsp
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_qrsp(q,r,s,p)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_2341_to_1234
!
!
   subroutine add_2341_to_1234_complex(scalar, x_qrsp, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 2341 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x_qrsp(q,r,s,p)
!!
      implicit none
!
      complex(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      complex(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      complex(dp), dimension(dim_q, dim_r, dim_s, dim_p), intent(in)    :: x_qrsp
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_qrsp(q,r,s,p)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_2341_to_1234_complex
!
!
   subroutine add_2143_to_1234(scalar, x_qpsr, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 2143 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x(q,p,s,r)
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      real(dp), dimension(dim_q, dim_p, dim_s, dim_r), intent(in)    :: x_qpsr
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_qpsr(q,p,s,r)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_2143_to_1234
!
!
   subroutine add_2143_to_1234_complex(scalar, x_qpsr, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 2143 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x(q,p,s,r)
!!
      implicit none
!
      complex(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      complex(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      complex(dp), dimension(dim_q, dim_p, dim_s, dim_r), intent(in)    :: x_qpsr
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_qpsr(q,p,s,r)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_2143_to_1234_complex
!
!
   subroutine add_2134_to_1234(scalar, x_qprs, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 2143 to 1234
!!    Written by Sarai D. Folkestad,
!!    Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x_qprs(q,p,r,s)
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      real(dp), dimension(dim_q, dim_p, dim_r, dim_s), intent(in)    :: x_qprs
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_qprs(q,p,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_2134_to_1234
!
!
   subroutine add_3214_to_1234(scalar, x_rqps, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 3214 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x(r,q,p,s)
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      real(dp), dimension(dim_r, dim_q, dim_p, dim_s), intent(in)    :: x_rqps
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_rqps(r,q,p,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_3214_to_1234
!
!
   subroutine add_3214_to_1234_complex(scalar, x_rqps, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 3214 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x(r,q,p,s)
!!
      implicit none
!
      complex(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      complex(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      complex(dp), dimension(dim_r, dim_q, dim_p, dim_s), intent(in)    :: x_rqps
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_rqps(r,q,p,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_3214_to_1234_complex
!
!
   subroutine add_4231_to_1234(scalar, x_sqrp, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 4231 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x(s,q,r,p)
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      real(dp), dimension(dim_s, dim_q, dim_r, dim_p), intent(in)    :: x_sqrp
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_sqrp(s,q,r,p)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_4231_to_1234
!
!
   subroutine add_2413_to_1234(scalar, x_qspr, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 2431 to 1234
!!    Written by Sarai D. Folkestad,
!!    Eirik F. Kjønstad and Alexander C. Paul, Jan 2019
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x_qspr(q,s,p,r)
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      real(dp), dimension(dim_q, dim_s, dim_p, dim_r), intent(in)    :: x_qspr
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_qspr(q,s,p,r)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_2413_to_1234
!
!
   subroutine add_2431_to_1234(scalar, x_qsrp, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 2431 to 1234
!!    Written by Sarai D. Folkestad,
!!    Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x_qsrp(q,s,r,p)
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      real(dp), dimension(dim_q, dim_s, dim_r, dim_p), intent(in)    :: x_qsrp
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_qsrp(q,s,r,p)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_2431_to_1234
!
!
   subroutine add_2431_to_1234_complex(scalar, x_qsrp, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 2431 to 1234
!!    Written by Sarai D. Folkestad,
!!    Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x_qsrp(q,s,r,p)
!!
      implicit none
!
      complex(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      complex(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      complex(dp), dimension(dim_q, dim_s, dim_r, dim_p), intent(in)    :: x_qsrp
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_qsrp(q,s,r,p)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_2431_to_1234_complex
!
!
   subroutine add_4213_to_1234(scalar, x_sqpr, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 4213 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x_sqpr(s,q,p,r)
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      real(dp), dimension(dim_s, dim_q, dim_p, dim_r), intent(in)    :: x_sqpr
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_sqpr(s,q,p,r)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_4213_to_1234
!
!
   subroutine add_4321_to_1234(scalar, x_srqp, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 4321 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x_srqp(s,r,q,p)
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      real(dp), dimension(dim_s, dim_r, dim_q, dim_p), intent(in)    :: x_srqp
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_srqp(s,r,q,p)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_4321_to_1234
!
!
   subroutine add_4321_to_1234_complex(scalar, x_srqp, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 4321 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x_srqp(s,r,q,p)
!!
      implicit none
!
      complex(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      complex(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      complex(dp), dimension(dim_s, dim_r, dim_q, dim_p), intent(in)    :: x_srqp
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_srqp(s,r,q,p)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_4321_to_1234_complex
!
!
   subroutine sort_1234_to_4321(x_pqrs, x_srqp, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 4321
!!    Written by Eirik F. Kjønstad, Rolf H. Myhre and Andreas Skeidsvoll, 2018
!!
!!    Reorders the array x_pqrs to x_srqp (i.e., 1234 to 4321).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_s, dim_r, dim_q, dim_p), intent(out) :: x_srqp
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do p = 1, dim_p
         do q = 1, dim_q
            do r = 1, dim_r
               do s = 1, dim_s
   !
                  x_srqp(s,r,q,p) = x_pqrs(p,q,r,s)
   !
               enddo
            enddo
         enddo
      enddo
   !$omp end parallel do
!
   end subroutine sort_1234_to_4321
!
!
   subroutine add_4123_to_1234(scalar, x_spqr, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 4123 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x_spqr(s,p,q,r)
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      real(dp), dimension(dim_s, dim_p, dim_q, dim_r), intent(in)    :: x_spqr
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_spqr(s,p,q,r)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_4123_to_1234
!
!
   subroutine add_4123_to_1234_complex(scalar, x_spqr, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 4123 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x_spqr(s,p,q,r)
!!
      implicit none
!
      complex(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      complex(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      complex(dp), dimension(dim_s, dim_p, dim_q, dim_r), intent(in)    :: x_spqr
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_spqr(s,p,q,r)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_4123_to_1234_complex
!
!
   subroutine add_4132_to_1234(scalar, x_sprq, y_pqrs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 4132 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar * x_sprq(s,p,r,q)
!!
      implicit none
!
      real(dp), intent(in) :: scalar
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(inout) :: y_pqrs
      real(dp), dimension(dim_s, dim_p, dim_r, dim_q), intent(in)    :: x_sprq
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do r = 1, dim_r
            do q = 1, dim_q
               do p = 1, dim_p
!
                  y_pqrs(p,q,r,s) = y_pqrs(p,q,r,s) + scalar*x_sprq(s,p,r,q)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_4132_to_1234
!
!
   subroutine squareup_and_sort_1234_to_1432(x_pqrs, x_psrq, dim_p, dim_q, dim_r, dim_s)
!!
!!    Squareup and sort 1234 to 1432
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the packed array x_pqrs to unpacked x_psrq (i.e., 1234 to 1432).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(((dim_p*dim_q+1)*(dim_r*dim_s))/2), intent(in) :: x_pqrs
      real(dp), dimension(dim_p, dim_s, dim_r, dim_q), intent(out)       :: x_psrq
!
      integer :: p, q, r, s, pq, rs, pqrs
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,pqrs)
      do q = 1, dim_q
         do r = 1, dim_r
            do s = 1, dim_s
!
               rs = dim_r*(s-1) + r
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
!
                  pqrs = (max(pq,rs)*(max(pq,rs)-3)/2) + pq + rs
!
                  x_psrq(p,s,r,q) = x_pqrs(pqrs)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine squareup_and_sort_1234_to_1432
!
!
   subroutine squareup_and_sort_1234_to_1432_complex(x_pqrs, x_psrq, dim_p, dim_q, dim_r, dim_s)
!!
!!    Squareup and sort 1234 to 1432
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Reorders the packed array x_pqrs to unpacked x_psrq (i.e., 1234 to 1432).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      complex(dp), dimension(((dim_p*dim_q+1)*(dim_r*dim_s))/2), intent(in) :: x_pqrs
      complex(dp), dimension(dim_p, dim_s, dim_r, dim_q), intent(out)       :: x_psrq
!
      integer :: p, q, r, s, pq, rs, pqrs
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,pqrs)
      do q = 1, dim_q
         do r = 1, dim_r
            do s = 1, dim_s
!
               rs = dim_r*(s-1) + r
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
!
                  pqrs = (max(pq,rs)*(max(pq,rs)-3)/2) + pq + rs
!
                  x_psrq(p,s,r,q) = x_pqrs(pqrs)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine squareup_and_sort_1234_to_1432_complex
!
!
   subroutine squareup_and_sort_1234_to_1243(x_pqrs, x_pqsr, dim_p, dim_q, dim_r, dim_s)
!!
!!    Squareup and sort 1234 to 1243
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the packed array x_pqrs to unpacked x_pqsr (i.e., 1234 to 1432).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(((dim_p*dim_q+1)*(dim_r*dim_s))/2), intent(in) :: x_pqrs
      real(dp), dimension(dim_p, dim_q, dim_s, dim_r), intent(out)       :: x_pqsr
!
      integer :: p, q, r, s, pq, rs, pqrs
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,pqrs)
      do r = 1, dim_r
         do s = 1, dim_s
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  pqrs = (max(pq,rs)*(max(pq,rs)-3)/2) + pq + rs
!
                  x_pqsr(p,q,s,r) = x_pqrs(pqrs)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine squareup_and_sort_1234_to_1243
!
!
   subroutine squareup_and_sort_1234_to_1423(x_pqrs, x_psqr, dim_p, dim_q, dim_r, dim_s)
!!
!!    Squareup and sort 1234 to 1423
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the packed array x_pqrs to unpacked x_psqr (i.e., 1234 to 1423).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(((dim_p*dim_q+1)*(dim_r*dim_s))/2), intent(in) :: x_pqrs
      real(dp), dimension(dim_p, dim_s, dim_q, dim_r), intent(out)       :: x_psqr
!
      integer :: p, q, r, s, pq, rs, pqrs
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,pqrs)
      do r = 1, dim_r
         do q = 1, dim_q
            do s = 1, dim_s
!
               rs = dim_r*(s-1) + r
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  pqrs = (max(pq,rs)*(max(pq,rs)-3)/2) + pq + rs
!
                  x_psqr(p,s,q,r) = x_pqrs(pqrs)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine squareup_and_sort_1234_to_1423
!
!
   subroutine sort_1234_to_3214(x_pqrs, x_rqps, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 3214
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pq_rs to x_rq_ps (i.e., 1234 to 3214).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_r, dim_q, dim_p, dim_s), intent(out) :: x_rqps
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do p = 1, dim_p
            do q = 1, dim_q
               do r = 1, dim_r
!
                  x_rqps(r,q,p,s) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_3214
!
!
   subroutine sort_1234_to_3214_complex(x_pqrs, x_rqps, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 3214
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Reorders the array x_pq_rs to x_rq_ps (i.e., 1234 to 3214).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      complex(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      complex(dp), dimension(dim_r, dim_q, dim_p, dim_s), intent(out) :: x_rqps
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do s = 1, dim_s
         do p = 1, dim_p
            do q = 1, dim_q
               do r = 1, dim_r
!
                  x_rqps(r,q,p,s) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_3214_complex
!
!
   subroutine sort_1234_to_4231(x_pqrs, x_sqrp, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 4231
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pqrs to x_sqrp (i.e., 1234 to 4231).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_s, dim_q, dim_r, dim_p), intent(out) :: x_sqrp
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do p = 1, dim_p
         do r = 1, dim_r
            do q = 1, dim_q
               do s = 1, dim_s
!
                  x_sqrp(s,q,r,p) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_4231
!
!
   subroutine sort_1234_to_4213(x_pqrs, x_sqpr, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 4213
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pq_rs to x_sq_pr (i.e., 1234 to 4213).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_s, dim_q, dim_p, dim_r), intent(out) :: x_sqpr
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do r = 1, dim_r
         do p = 1, dim_p
            do q = 1, dim_q
               do s = 1, dim_s
!
                  x_sqpr(s,q,p,r) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_4213
!
!
   subroutine sort_1234_to_3241(x_pqrs, x_rqsp, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 3241
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pqrs to x_rqsp (i.e., 1234 to 3241).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_r, dim_q, dim_s, dim_p), intent(out) :: x_rqsp
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do p = 1, dim_p
         do s = 1, dim_s
            do q = 1, dim_q
               do r = 1, dim_r
!
                  x_rqsp(r,q,s,p) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_3241
!
!
   subroutine sort_1234_to_3241_complex(x_pqrs, x_rqsp, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 3241
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Reorders the array x_pqrs to x_rqsp (i.e., 1234 to 3241).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      complex(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      complex(dp), dimension(dim_r, dim_q, dim_s, dim_p), intent(out) :: x_rqsp
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(s,r,q,p)
      do p = 1, dim_p
         do s = 1, dim_s
            do q = 1, dim_q
               do r = 1, dim_r
!
                  x_rqsp(r,q,s,p) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_3241_complex
!
!
   subroutine sort_1234_to_1243(x_pqrs, x_pqsr, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 1243
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pqrs to x_pqsr (i.e., 1234 to 1243).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      real(dp), dimension(dim_p, dim_q, dim_s, dim_r), intent(out) :: x_pqsr
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(p, q, r, s)
      do r = 1,dim_r
         do s = 1,dim_s
            do q = 1, dim_q
               do p = 1, dim_p
!
                  x_pqsr(p,q,s,r) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_1243
!
!
   subroutine sort_1234_to_1243_complex(x_pqrs, x_pqsr, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 1243
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Reorders the array x_pqrs to x_pqsr (i.e., 1234 to 1243).
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      complex(dp), dimension(dim_p, dim_q, dim_r, dim_s), intent(in)  :: x_pqrs
      complex(dp), dimension(dim_p, dim_q, dim_s, dim_r), intent(out) :: x_pqsr
!
      integer :: p, q, r, s
!
!$omp parallel do schedule(static) private(p, q, r, s)
      do r = 1,dim_r
         do s = 1,dim_s
            do q = 1, dim_q
               do p = 1, dim_p
!
                  x_pqsr(p,q,s,r) = x_pqrs(p,q,r,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_1243_complex
!
!     -::- Squareup and packin and related routines -::-
!     --------------------------------------------------
!
   subroutine add_to_packed_real(packed, unpacked, N)
!!
!!    Add to packed
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2017
!!
!!    Adds a symmetric unpacked N x N matrix to a packed N x N
!!    matrix.
!!
      implicit none
!
      integer, intent(in) :: N
!
      real(dp), dimension(N*(N+1)/2), intent(inout) :: packed
      real(dp), dimension(N, N), intent(in)         :: unpacked
!
      integer :: i, j, ij
!
!$omp parallel do schedule(static) private(i,j, ij)
      do i = 1, N
         do j = i, N
            ij = (max(i,j)*(max(i,j)-3)/2) + i + j
            packed(ij) = packed(ij) + unpacked(i,j)
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_to_packed_real
!
!
   subroutine symmetrize_and_add_to_packed_real(packed, unpacked, N)
!!
!!    Symmetrize and add to packed
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2017
!!
!!    Symmetrizes unpacked N x N matrix and adds the symmtrized sum
!!    to a packed N x N matrix.
!!
      implicit none
!
      integer, intent(in) :: N
!
      real(dp), dimension(N*(N+1)/2), intent(out) :: packed
      real(dp), dimension(N,N), intent(in)        :: unpacked
!
      integer :: i, j, ij
!
!$omp parallel do schedule(static) private(i,j, ij)
      do i = 1, N
         do j = i, N
            ij = (max(i,j)*(max(i,j)-3)/2) + i + j
            packed(ij) = packed(ij) + unpacked(i,j) + unpacked(j,i)
         enddo
      enddo
!$omp end parallel do
!
   end subroutine symmetrize_and_add_to_packed_real
!
!
   subroutine symmetrize_and_add_4_to_packed(packed, unpacked, N)
!!
!!    Symmetrize and add 4 to packed
!!    Written by Andreas Skeidsvoll, Sep 2019
!!
!!    Symmetrizes 4-dimensional unpacked matrix and adds the symmetrized sum
!!    to a packed matrix.
!!
      implicit none
!
      integer, intent(in) :: N
!
      real(dp), dimension(N*(N+1)/2), intent(out) :: packed
      real(dp), dimension(:,:,:,:), intent(in)    :: unpacked
!
      call symmetrize_and_add_to_packed_real(packed, unpacked, N)
!
   end subroutine symmetrize_and_add_4_to_packed
!
!
   subroutine symmetrize_add_contra_to_packed(x, y, dim_p, dim_q)
!!
!!    Symmetrize and add contra to packed
!!    Written by Rolf H. Myhre and Alexander C. Paul, Aug 2020
!!
!!    Symmetrize and add an unsymmetric unpacked pq x rs matrix ordered prqs
!!       (2 X_pqrs - X_rqps)
!!    to a packed symmetric matrix Y_pqrs
!!
!!    Performs:
!!
!!       Y_pqrs = 1/3(2X_pqrs + X_rqps + 2X_rspq + X_psrq)
!!
!!    Note that this routine requires dim_p = dim_r and dim_q = dim_s
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q
!
      real(dp), dimension(dim_p, dim_p, dim_q, dim_q), intent(in)       :: x
      real(dp), dimension(dim_p*dim_q*(dim_p*dim_q+1)/2), intent(inout) :: y
!
      integer :: p, q, r, s, pq, rs, pqrs, p_end
!
!$omp parallel do schedule(static) collapse(2) private(s,r,q,p,pq,rs,pqrs,p_end)
      do s = 1, dim_q
         do r = 1, dim_p
            do q = 1, s
!
               if (s .ne. q) then
                  p_end = dim_p
               else
                  p_end = r
               end if
!
               do p = 1, p_end
!
                  pq = dim_p*(q-1)+p
                  rs = dim_p*(s-1)+r
                  pqrs = max(pq,rs)*(max(pq,rs)-3)/2 + pq + rs
!
                  y(pqrs) = y(pqrs) &
                          + third*(two*x(p,r,q,s) + x(r,p,q,s) &
                          +        two*x(r,p,s,q) + x(p,r,s,q))
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine symmetrize_add_contra_to_packed
!
!
   subroutine symmetrize_add_to_packed(x, y, dim_p, dim_q)
!!
!!    Symmetrize and add to packed
!!    Written by Rolf H. Myhre and Alexander C. Paul, Aug 2020
!!
!!    Symmetrize and add an unsymmetric unpacked pq x rs matrix X_prqs
!!    to a packed symmetric matrix Y_pqrs
!!
!!    Note that this routine requires dim_p = dim_r and dim_q = dim_s
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q
!
      real(dp), dimension(dim_p, dim_p, dim_q, dim_q), intent(in)       :: x
      real(dp), dimension(dim_p*dim_q*(dim_p*dim_q+1)/2), intent(inout) :: y
!
      integer :: p, q, r, s, pq, rs, pqrs, p_end
!
!$omp parallel do schedule(static) collapse(2) private(s,r,q,p,pq,rs,pqrs,p_end)
      do s = 1, dim_q
         do r = 1, dim_p
            do q = 1, s
!
               if (s .ne. q) then
                  p_end = dim_p
               else
                  p_end = r
               end if
!
               do p = 1, p_end
!
                  pq = dim_p*(q-1)+p
                  rs = dim_p*(s-1)+r
                  pqrs = max(pq,rs)*(max(pq,rs)-3)/2 + pq + rs
!
                  y(pqrs) = y(pqrs) + (x(p,r,q,s) + x(r,p,s,q))
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine symmetrize_add_to_packed
!
!
   subroutine squareup_anti(packed, unpacked, N)
!!
!!    Square up packed antisymmetric matrix
!!    Written by Eirik F. Kjønstad, 2018
!!
!!    Squares up to full dimension (N x N) of packed matrix. The packed
!!    antisymmetric matrix contains the strictly lower triangular part
!     of the full unpacked matrix.
!!
      implicit none
!
      integer, intent(in) :: N
!
      real(dp), dimension(N*(N-1)/2), intent(in)  :: packed
      real(dp), dimension(N,N), intent(out)       :: unpacked
!
      integer :: i, j
!
!     Set diagonal to zero
!
      do i = 1, N
!
         unpacked(i, i) = zero
!
      enddo
!
!     Set lower and upper strictly triangular parts
!
      do i = 2, N
         do j = 1, i - 1
!
            unpacked(i, j) = packed((max(i-1,j)*(max(i-1,j)-3)/2) + i-1 + j)
            unpacked(j, i) = -packed((max(i-1,j)*(max(i-1,j)-3)/2) + i-1 + j)
!
         enddo
      enddo
!
   end subroutine squareup_anti
!
!
   subroutine squareup_anti_complex(packed, unpacked, N)
!!
!!    Square up packed antisymmetric matrix
!!    Written by Eirik F. Kjønstad, 2018
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Squares up to full dimension (N x N) of packed matrix. The packed
!!    antisymmetric matrix contains the strictly lower triangular part
!     of the full unpacked matrix.
!!
      implicit none
!
      integer, intent(in) :: N
!
      complex(dp), dimension(N*(N-1)/2), intent(in)  :: packed
      complex(dp), dimension(N,N), intent(out)       :: unpacked
!
      integer :: i, j
!
!     Set diagonal to zero_complex
!
      do i = 1, N
!
         unpacked(i, i) = zero_complex
!
      enddo
!
!     Set lower and upper strictly triangular parts
!
      do i = 2, N
         do j = 1, i - 1
!
            unpacked(i, j) = packed((max(i-1,j)*(max(i-1,j)-3)/2) + i-1 + j)
            unpacked(j, i) = -packed((max(i-1,j)*(max(i-1,j)-3)/2) + i-1 + j)
!
         enddo
      enddo
!
   end subroutine squareup_anti_complex
!
!
   subroutine squareup_real(packed,unpacked,N)
!!
!!    Square up packed symmetric matrix
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2017
!!
!!    Squares up to full dimension (N x N) of packed matrices.
!!
      implicit none
!
      integer, intent(in) :: N
!
      real(dp), dimension(N*(N+1)/2), intent(in) :: packed
      real(dp), dimension(N,N), intent(out)      :: unpacked
!
      integer :: i, j
!
!$omp parallel do schedule(static) private(i,j)
      do j = 1, N
         do i = 1, N
            unpacked(i, j) = packed((max(i,j)*(max(i,j)-3)/2) + i + j)
         enddo
      enddo
!$omp end parallel do
!
   end subroutine squareup_real
!
!
   subroutine squareup_complex(packed,unpacked,N)
!!
!!    Square up packed symmetric matrix
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2017
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Squares up to full dimension (N x N) of packed matrices.
!!
      implicit none
!
      integer, intent(in) :: N
!
      complex(dp), dimension(N*(N+1)/2), intent(in) :: packed
      complex(dp), dimension(N,N), intent(out)      :: unpacked
!
      integer :: i, j
!
!$omp parallel do schedule(static) private(i,j)
      do j = 1, N
         do i = 1, N
            unpacked(i, j) = packed((max(i,j)*(max(i,j)-3)/2) + i + j)
         enddo
      enddo
!$omp end parallel do
!
   end subroutine squareup_complex
!
!
   subroutine squareup_to_4(packed,unpacked,N)
!!
!!    Square up packed symmetric matrix
!!    Written by Andreas Skeidsvoll, Sep 2019
!!
!!    Squares up to full dimension of four dimensional packed matrices.
!!
      implicit none
!
      integer, intent(in) :: N
!
      real(dp), dimension(N*(N+1)/2), intent(in) :: packed
      real(dp), dimension(:,:,:,:), intent(out) :: unpacked
!
      call squareup_real(packed, unpacked, N)
!
   end subroutine squareup_to_4
!
!
   subroutine squareup_to_4_complex(packed,unpacked,N)
!!
!!    Square up packed symmetric matrix
!!    Written by Andreas Skeidsvoll, Sep 2019
!!
!!    Squares up to full dimension of four dimensional packed matrices.
!!
      implicit none
!
      integer, intent(in) :: N
!
      complex(dp), dimension(N*(N+1)/2), intent(in) :: packed
      complex(dp), dimension(:,:,:,:), intent(out)  :: unpacked
!
      call squareup_complex(packed, unpacked, N)
!
   end subroutine squareup_to_4_complex
!
!
   subroutine packin_real(packed,unpacked,N)
!!
!!    Pack in symmetric matrix
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2017
!!
!!    Pack down full square matrix of dimension N x N.
!!
      implicit none
!
      integer, intent(in) :: N
!
      real(dp), dimension(N*(N+1)/2), intent(out) :: packed
      real(dp), dimension(n,N), intent(in)        :: unpacked
!
      integer :: i, j
!
!$omp parallel do schedule(static) private(i,j)
      do i = 1, N
         do j = 1, i
!
            packed((max(i,j)*(max(i,j)-3)/2) + i + j) = unpacked(i, j)
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine packin_real
!
!
   subroutine packin_complex(packed,unpacked,N)
!!
!!    Pack in symmetric matrix
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2017
!!    Modified by Andreas Skeidsvoll, Sep 2019: Changed real arrays to complex
!!
!!    Pack down full square matrix of dimension N x N.
!!
      implicit none
!
      integer, intent(in) :: N
!
      complex(dp), dimension(N*(N+1)/2), intent(out) :: packed
      complex(dp), dimension(n,N), intent(in)        :: unpacked
!
      integer :: i, j
!
!$omp parallel do schedule(static) private(i,j)
      do i = 1, N
         do j = 1, i
!
            packed((max(i,j)*(max(i,j)-3)/2) + i + j) = unpacked(i, j)
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine packin_complex
!
!
   subroutine packin_4_real(packed,unpacked,N)
!!
!!    Pack in 4-dimensional symmetric matrix
!!    Written by Andreas Skeidsvoll, Sep 2019
!!
!!    Pack down full 4-dimensional square matrix.
!!
      implicit none
!
      integer, intent(in) :: N
!
      real(dp), dimension(N*(N+1)/2), intent(out) :: packed
      real(dp), dimension(:,:,:,:), intent(in)    :: unpacked
!
      call packin_real(packed, unpacked, N)
!
   end subroutine packin_4_real
!
!
   subroutine packin_4_complex(packed,unpacked,N)
!!
!!    Pack in 4-dimensional symmetric matrix
!!    Written by Andreas Skeidsvoll, Sep 2019
!!
!!    Pack down full 4-dimensional square matrix.
!!
      implicit none
!
      integer, intent(in) :: N
!
      complex(dp), dimension(N*(N+1)/2), intent(out) :: packed
      complex(dp), dimension(:,:,:,:), intent(in)    :: unpacked
!
      call packin_complex(packed, unpacked, N)
!
   end subroutine packin_4_complex
!
!
   subroutine packin_4_from_1324_order_real(packed, unpacked, dim_p, dim_q)
!!
!!    Pack in symmetric matrix ordered 1324
!!    Written by Rolf H. Myhre, Oct 2019
!!
!!    Pack in unpacked array X_prqs to packed array Y_pqrs
!!    where dim_p = dim_r and dim_q = dim_s
!!
      implicit none
!
      integer, intent(in) :: dim_p
      integer, intent(in) :: dim_q
!
      real(dp), dimension(dim_p*dim_q*(dim_p*dim_q+1)/2), intent(out) :: packed
      real(dp), dimension(dim_p,dim_p,dim_q,dim_q), intent(in)        :: unpacked
!
      integer :: p, q, r, s, pq, rs, pqrs, r_end
!
!$omp parallel do schedule(static) private(p, q, r, s, pq, rs, pqrs, r_end)
      do q = 1, dim_q
         do p = 1, dim_p
!
            pq = dim_p*(q - 1) + p
!
            do s = 1, q
!
               if (s .ne. q) then
                  r_end = dim_p
               else
                  r_end = p
               endif
!
               do r = 1, r_end
!
                  rs = dim_p*(s - 1) + r
!
                  pqrs = pq*(pq-3)/2 + pq + rs
!
                  packed(pqrs) = unpacked(p,r,q,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine packin_4_from_1324_order_real
!
!
   subroutine packin_4_from_1324_order_complex(packed, unpacked, dim_p, dim_q)
!!
!!    Pack in symmetric complex matrix ordered 1324
!!    Written by Rolf H. Myhre, Oct 2019
!!
!!    Pack in unpacked array X_prqs to packed array Y_pqrs
!!    where dim_p = dim_r and dim_q = dim_s
!!
      implicit none
!
      integer, intent(in) :: dim_p
      integer, intent(in) :: dim_q
!
      complex(dp), dimension(dim_p*dim_q*(dim_p*dim_q+1)/2), intent(out) :: packed
      complex(dp), dimension(dim_p,dim_p,dim_q,dim_q), intent(in)        :: unpacked
!
      integer :: p, q, r, s, pq, rs, pqrs, r_end
!
!$omp parallel do schedule(static) private(p, q, r, s, pq, rs, pqrs, r_end)
      do q = 1, dim_q
         do p = 1, dim_p
!
            pq = dim_p*(q - 1) + p
!
            do s = 1, q
!
               if (s .ne. q) then
                  r_end = dim_p
               else
                  r_end = p
               endif
!
               do r = 1, r_end
!
                  rs = dim_p*(s - 1) + r
!
                  pqrs = (pq*(pq-1))/2 + rs
!
                  packed(pqrs) = unpacked(p,r,q,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine packin_4_from_1324_order_complex
!
!
   subroutine packin_and_add_from_1324_order(packed, unpacked, dim_p, dim_q)
!!
!!    Pack in and add symmetric matrix ordered 1324
!!    Written by Rolf H. Myhre, Oct 2019
!!
!!    Almost identical to non-adding routine by Rolf H. Myhre, Oct 2019
!!
!!    Pack in unpacked array X_prqs to packed array Y_pqrs
!!    where dim_p = dim_r and dim_q = dim_s
!!
      implicit none
!
      integer, intent(in) :: dim_p
      integer, intent(in) :: dim_q
!
      real(dp), dimension(dim_p*dim_q*(dim_p*dim_q+1)/2), intent(inout) :: packed
      real(dp), dimension(dim_p,dim_p,dim_q,dim_q), intent(in)          :: unpacked
!
      integer :: p, q, r, s, pq, rs, pqrs, r_end
!
!$omp parallel do schedule(static) private(p, q, r, s, pq, rs, pqrs, r_end)
      do q = 1, dim_q
         do p = 1, dim_p
!
            pq = dim_p*(q - 1) + p
!
            do s = 1, q
!
               if (s .ne. q) then
                  r_end = dim_p
               else
                  r_end = p
               endif
!
               do r = 1, r_end
!
                  rs = dim_p*(s - 1) + r
!
                  pqrs = pq*(pq-3)/2 + pq + rs
!
                  packed(pqrs) = packed(pqrs) + unpacked(p,r,q,s)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine packin_and_add_from_1324_order
!
!
   subroutine packin_anti(packed, unpacked, N)
!!
!!    Pack in anti-symmetric matrix
!!    Written by Eirik F. Kjønstad, 2018
!!
!!    Pack down full square anti-symmetric matrix of dimension N x N,
!!    where the strictly lower triangular part of the unpacked matrix
!!    is stored in packed form.
!!
      implicit none
!
      integer, intent(in) :: N
!
      real(dp), dimension(N*(N-1)/2), intent(out) :: packed
      real(dp), dimension(N,N), intent(in)        :: unpacked
!
      integer :: i, j
!
      do i = 2, N
         do j = 1, i - 1
!
            packed((max(i-1,j)*(max(i-1,j)-3)/2) + i-1 + j) = unpacked(i, j)
!
         enddo
      enddo
!
   end subroutine packin_anti
!
!
   subroutine construct_packed_contravariant(x, y, dim_p, dim_q)
!!
!!    Construct packed contravariant
!!    Written by Rolf H. Myhre and Alexander C. Paul, Aug 2020
!!
!!    Constructs
!!       Y_pqrs = (2 X_pqrs - X_rqps) for packed arrays X and Y
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q
!
      real(dp), dimension(dim_p*dim_q*(dim_p*dim_q+1)/2), intent(in)    :: x
      real(dp), dimension(dim_p*dim_q*(dim_p*dim_q+1)/2), intent(inout) :: y
!
      integer :: p, q, r, s, pq, rs, pqrs, p_end, rq, ps, rqps
!
!$omp parallel do schedule(static) collapse(2) &
!$omp private(p, q, r, s, rs, rq, pq, ps, pqrs, rqps, p_end)
      do s = 1, dim_q
         do r = 1, dim_p
!
            rs = dim_p*(s-1)+r
!
            do q = 1, s
!
               rq = dim_p*(q-1)+r
!
               if (s .ne. q) then
                  p_end = dim_p
               else
                  p_end = r
               end if
!
               do p = 1, p_end
!
                  pq = dim_p*(q-1)+p
                  ps = dim_p*(s-1)+p
                  pqrs = rs*(rs-3)/2 + pq + rs
                  rqps = max(rq,ps)*(max(rq,ps)-3)/2 + rq + ps
!
                  y(pqrs) = (two*x(pqrs) - x(rqps))
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine construct_packed_contravariant
!
!
   subroutine construct_packed_covariant(x, y, dim_p, dim_q)
!!
!!    Construct packed covariant
!!    Written by Rolf H. Myhre and Alexander C. Paul, Aug 2020
!!
!!    Constructs:
!!       Y_pqrs = 1/3(2X_pqrs + X_rqps)
!!
!!    Note that this routine requires dim_p = dim_r and dim_q = dim_s
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q
!
      real(dp), dimension(dim_p*dim_q*(dim_p*dim_q+1)/2), intent(in)  :: x
      real(dp), dimension(dim_p*dim_q*(dim_p*dim_q+1)/2), intent(out) :: y
!
      integer :: p, q, r, s, pq, rs, pqrs, p_end, rq, ps, rqps
!
!$omp parallel do schedule(static) collapse(2) &
!$omp private(p,q,r,s,rs, rq, pq, ps, pqrs, rqps, p_end)
      do s = 1, dim_q
         do r = 1, dim_p
!
            rs = dim_p*(s-1)+r
!
            do q = 1, s
!
               rq = dim_p*(q-1)+r
!
               if (s .ne. q) then
                  p_end = dim_p
               else
                  p_end = r
               end if
!
               do p = 1, p_end
!
                  pq = dim_p*(q-1)+p
                  ps = dim_p*(s-1)+p
                  pqrs = rs*(rs-3)/2 + pq + rs
                  rqps = max(rq,ps)*(max(rq,ps)-3)/2 + rq + ps
!
                  y(pqrs) = third*(two*x(pqrs) + x(rqps))
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine construct_packed_covariant
!
!
   subroutine construct_covariant_1324(x_pqrs, y_prqs, dim_p, dim_q)
!!
!!    Construct covariant 1324
!!    Written by Alexander C. Paul, Jan 2021
!!
!!    Constructs:
!!       Y_prqs = 1/3(2X_pqrs + X_rqps)
!!    from packed vector x_pqrs
!!
!!    Note that this routine requires dim_p = dim_r and dim_q = dim_s
!!
      implicit none
!
      integer, intent(in) :: dim_p, dim_q
!
      real(dp), dimension(((dim_p*dim_q+1)*(dim_p*dim_q)/2)), intent(in) :: x_pqrs
      real(dp), dimension(dim_p, dim_p, dim_q, dim_q), intent(out)       :: y_prqs
!
      integer :: p, q, r, s, pq, rs, ps, rq, pqrs, rqps
!
!$omp parallel do schedule(static) private(s,r,q,p, rs, rq, pq, ps, pqrs, rqps)
      do s = 1, dim_q
         do q = 1, dim_q
            do r = 1, dim_p
!
               rs = dim_p*(s-1) + r
               rq = dim_p*(q-1) + r
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  ps = dim_p*(s-1) + p
!
                  pqrs = max(pq,rs)*(max(pq,rs)-3)/2 + pq + rs
                  rqps = max(rq,ps)*(max(rq,ps)-3)/2 + rq + ps
!
                  y_prqs(p,r,q,s) = third*(two*x_pqrs(pqrs) + x_pqrs(rqps))
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine construct_covariant_1324
!
!
   subroutine construct_contravariant_t3(t, u, dim_)
!!
!!    Construct contravariant t
!!    Written by Rolf H. Myhre and Alexander C. Paul, August 2020
!!
!!    Construct:
!!       t_pqr = 4t_pqr - 2t_qpr - 2t_rqp - 2t_prq + t_qrp + t_rpq
!!
      use array_utilities, only: copy_and_scale
!
      implicit none
!
      integer, intent(in) :: dim_
!
      real(dp), dimension(dim_, dim_, dim_), intent(inout) :: t
      real(dp), dimension(dim_, dim_, dim_), intent(out)   :: u
!
      call copy_and_scale( two, t, u, dim_**3)
      call add_132_to_123(-one, t, u, dim_, dim_, dim_)
      call add_321_to_123(-one, t, u, dim_, dim_, dim_)
!
      call copy_and_scale( two, u, t, dim_**3)
!
      call add_213_to_123(-one, u, t, dim_, dim_, dim_)
!
   end subroutine construct_contravariant_t3
!
!
end module reordering
