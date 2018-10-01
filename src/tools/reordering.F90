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
   use kinds
   use parameters
   use index
   use file_class
!
   implicit none
!
contains
!
!     -::- Two-index re-sort and re-sort-add routines -::-
!     ----------------------------------------------------
!
   subroutine sort_12_to_21(x_p_q, x_q_p, dim_p, dim_q)
!!
!!    Sort 12 to 21
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_p_q to x_q_p (i.e., 12 to 21).
!!
!!    The unordered array x_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s.
!!    The ordered array x_sp_rq is assumed allocated as dim_s*dim_p x dim_r*dim_q.
!!
      implicit none
!
      integer(i15), intent(in) :: dim_p, dim_q
!
      real(dp), dimension(dim_p, dim_q), intent(in) :: x_p_q
      real(dp), dimension(dim_q, dim_p) :: x_q_p
!
      integer(i15) :: p, q
!
      do q = 1, dim_q
         do p = 1, dim_p
!
            x_q_p(q, p) = x_p_q(p, q)
!
         enddo
      enddo
!
   end subroutine sort_12_to_21
!
!
   subroutine add_21_to_12(gamma, x, y_p_q, dim_p, dim_q)
!!
!!    Add 21 to 12
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_p_q(p,q) = y_p_q(p,q) + gamma*x(q,p)
!!
!!    The unordered array y_p_q is assumed allocated as dim_p x dim_q,
!!    and x accordingly.
!!
      implicit none
!
      real(dp), intent(in) :: gamma
!
      integer(i15), intent(in) :: dim_p, dim_q
!
      real(dp), dimension(dim_p, dim_q) :: y_p_q
      real(dp), dimension(dim_q, dim_p), intent(in) :: x
!
      integer(i15) :: p, q
!
      do q = 1, dim_q
         do p = 1, dim_p
!
               y_p_q(p,q) = y_p_q(p,q) + gamma*x(q,p)
!
         enddo
      enddo
!
   end subroutine add_21_to_12
!
!
   subroutine symmetric_sum(x, dim)
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
!!    Note: the effect is equal to the routine "add_21_to_12" with gamma = 1,
!!          where a second copy of x is not needed.
!!
      implicit none
!
      integer(i15), intent(in) :: dim
!
      real(dp), dimension(dim, dim), intent(inout) :: x
!
      integer(i15) :: p, q
!
!     Overwrite the lower triangular part of the matrix
!
!$omp parallel do private(p, q)
      do q = 1, dim
         do p = q, dim
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
      do p = 1, dim
         do q = p + 1, dim
!
            x(p,q) = x(q,p)
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine symmetric_sum
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
!!       x_rpq(rpq, 1) = x_pqr(pqr, 1)
!!
      implicit none
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r
!
      real(dp), dimension(:, :), intent(in)    :: x_pqr
      real(dp), dimension(:, :), intent(inout) :: x_rpq
!
      integer(i15) :: pqr, rpq, r, q, p
!
!$omp parallel do schedule(static) private(r,q,p,pqr, rpq)
      do r = 1, dim_r
         do q = 1, dim_q
            do p = 1, dim_p
!
               pqr = dim_p*(dim_q*(r-1)+q-1)+p
               rpq = dim_r*(dim_p*(q-1)+p-1)+r
!
               x_rpq(rpq, 1) = x_pqr(pqr, 1)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_123_to_312
!
!
  subroutine sort_123_to_321(x_pqr, x_rqp, dim_p, dim_q, dim_r)
!!
!!    Sort 123 to 321
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       x_rqp(rqp, 1) = x_pqr(pqr, 1)
!!
     implicit none
!
     integer(i15), intent(in) :: dim_p, dim_q, dim_r
!
     real(dp), dimension(:, :), intent(in)    :: x_pqr
     real(dp), dimension(:, :), intent(inout) :: x_rqp
!
     integer(i15) :: pqr, rqp, r, q, p
!
!$omp parallel do schedule(static) private(r,q,p,pqr, rqp)
     do r = 1, dim_r
        do q = 1, dim_q
           do p = 1, dim_p
!
              pqr = dim_p*(dim_q*(r-1)+q-1)+p
              rqp = dim_r*(dim_q*(p-1)+q-1)+r
!
              x_rqp(rqp, 1) = x_pqr(pqr, 1)
!
           enddo
        enddo
     enddo
!$omp end parallel do
!
  end subroutine sort_123_to_321
!
!
   subroutine sort_123_to_132(x_pqr, x_prq, dim_p, dim_q, dim_r)
!!
!!    Sort 123 to 132
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       x_prq(prq, 1) = x_pqr(pqr, 1)
!!
      implicit none
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r
!
      real(dp), dimension(:, :), intent(in)    :: x_pqr
      real(dp), dimension(:, :), intent(inout) :: x_prq
!
      integer(i15) :: pqr, prq, r, q, p
!
!$omp parallel do schedule(static) private(r,q,p,pqr,prq)
      do r = 1, dim_r
         do q = 1, dim_q
            do p = 1, dim_p
!
               pqr = dim_p*(dim_q*(r-1)+q-1)+p
               prq = dim_p*(dim_r*(q-1)+r-1)+p
!
               x_prq(prq, 1) = x_pqr(pqr, 1)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_123_to_132
!
!
   subroutine sort_123_to_213(x_pqr, x_qpr, dim_p, dim_q, dim_r)
!!
!!    Sort 123 to 132
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       x_qpr(qpr, 1) = x_pqr(pqr, 1)
!!
      implicit none
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r
!
      real(dp), dimension(:, :), intent(in)    :: x_pqr
      real(dp), dimension(:, :), intent(inout) :: x_qpr
!
      integer(i15) :: pqr, qpr, r, q, p
!
!$omp parallel do schedule(static) private(r,q,p,pqr,qpr)
      do r = 1, dim_r
         do q = 1, dim_q
            do p = 1, dim_p
!
               pqr = dim_p*(dim_q*(r-1)+q-1)+p
               qpr = dim_q*(dim_p*(r-1)+p-1)+q
!
               x_qpr(qpr, 1) = x_pqr(pqr, 1)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_123_to_213
!
!
   subroutine add_321_to_123(gamma, x, y_pqr, dim_p, dim_q, dim_r)
!!
!!    Add 1243 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pqr(pqr,1) = y_pqr(pqr,1) + gamma*x(rqp,1)
!!
!!
      implicit none
!
      real(dp), intent(in) :: gamma
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r
!
      real(dp), dimension(:, :) :: y_pqr
      real(dp), dimension(:, :), intent(in) :: x
!
      integer(i15) :: pqr, rqp, r, q, p
!
!$omp parallel do schedule(static) private(r,q,p,pqr,rqp)
      do r = 1, dim_r
         do q = 1, dim_q
            do p = 1, dim_p
!
               pqr = dim_p*(dim_q*(r-1)+q-1)+p
               rqp = dim_r*(dim_q*(p-1)+q-1)+r
!
               y_pqr(pqr, 1) = y_pqr(pqr, 1) + gamma*x(rqp, 1)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_321_to_123
!
!
   subroutine add_132_to_123(gamma, x, y_pqr, dim_p, dim_q, dim_r)
!!
!!    Add 132 to 123
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pqr(pqr,1) = y_pqr(pqr, 1) + gamma * x(prq, 1)
!!
      implicit none
!
      real(dp), intent(in) :: gamma
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r
!
      real(dp), dimension(dim_p*dim_q*dim_r, 1) :: y_pqr
      real(dp), dimension(dim_r*dim_p*dim_q, 1), intent(in) :: x
!
      integer(i15) :: p, q, r, pqr, prq
!
!$omp parallel do schedule(static) private(r,q,p,pqr,prq)
      do r = 1, dim_r
         do q = 1, dim_q
            do p = 1, dim_p
!
               pqr = dim_p*(dim_q*(r-1)+q-1)+p
               prq = dim_p*(dim_r*(q-1)+r-1)+p
!
               y_pqr(pqr, 1) = y_pqr(pqr, 1) + gamma*x(prq, 1)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_132_to_123
!
!
   subroutine add_312_to_123(gamma, x, y_pqr, dim_p, dim_q, dim_r)
!!
!!    Add 312 to 123
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pqr(pqr,1) = y_pqr(pqr, 1) + gamma * x(rpq, 1)
!!
      implicit none
!
      real(dp), intent(in) :: gamma
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r
!
      real(dp), dimension(dim_p*dim_q*dim_r, 1) :: y_pqr
      real(dp), dimension(dim_r*dim_p*dim_q, 1), intent(in) :: x
!
      integer(i15) :: p, q, r, pqr, rpq
!
!$omp parallel do schedule(static) private(r,q,p,pqr,rpq)
      do r = 1, dim_r
         do q = 1, dim_q
            do p = 1, dim_p
!
               pqr = dim_p*(dim_q*(r-1)+q-1)+p
               rpq = dim_r*(dim_p*(q-1)+p-1)+r
!
               y_pqr(pqr, 1) = y_pqr(pqr, 1) + gamma*x(rpq, 1)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_312_to_123
!
!     -::- Four-index re-sort and re-sort-add routines -::-
!     -----------------------------------------------------
!
   subroutine add_3124_to_1234(gamma, x, y_pq_rs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 3124 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pq_rs(pq,rs) = y_pq_rs(pq,rs) + gamma * x(rp, qs)
!!
!!    The unordered array y_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s,
!!    and x accordingly.
!!
      implicit none
!
      real(dp), intent(in) :: gamma
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s) :: y_pq_rs
      real(dp), dimension(dim_r*dim_p, dim_q*dim_s), intent(in) :: x
!
      integer(i15) :: p, q, r, s, pq, rs, rp, qs
!
!$omp parallel do schedule(static) private(s,r,q,p,rs,qs,rp,pq)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
!
               qs = dim_q*(s-1) + q
!
               do p = 1, dim_p
!
                  rp = dim_r*(p-1) + r
                  pq = dim_p*(q-1) + p
!
                  y_pq_rs(pq, rs) = y_pq_rs(pq, rs) + gamma*x(rp, qs)
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
   subroutine sort_1234_to_3412(x_pq_rs, x_rs_pq, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 3412
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pq_rs to x_rs_pq (i.e., 1234 to 3412).
!!
!!    The unordered array x_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s.
!!    The ordered array x_rs_pq is assumed allocated as dim_r*dim_s x dim_p*dim_q.
!!
      implicit none
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s), intent(in)    :: x_pq_rs
      real(dp), dimension(dim_s*dim_p, dim_r*dim_q), intent(inout) :: x_rs_pq
!
      integer(i15) :: p, q, r, s, rs, pq
!
!$omp parallel do schedule(static) private(s,r,rs,q,p,pq)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
!
                  x_rs_pq(rs, pq) = x_pq_rs(pq, rs)
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
   subroutine sort_1234_to_4132(x_pq_rs, x_sp_rq, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 4132
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pq_rs to x_sp_rq (i.e., 1234 to 4132).
!!
!!    The unordered array x_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s.
!!    The ordered array x_sp_rq is assumed allocated as dim_s*dim_p x dim_r*dim_q.
!!
      implicit none
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s), intent(in) :: x_pq_rs
      real(dp), dimension(dim_s*dim_p, dim_r*dim_q), intent(inout) :: x_sp_rq
!
      integer(i15) :: p, q, r, s, rs, pq, sp, rq
!
!$omp parallel do schedule(static) private(s,r,rs,q,rq,p,pq,sp)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
!
               rq = dim_r*(q-1) + r
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  sp = dim_s*(p-1) + s
!
                  x_sp_rq(sp, rq) = x_pq_rs(pq, rs)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_4132
!
   subroutine squareup_and_sort_1234_to_4132(x_pqrs, x_sp_rq, dim_p, dim_q, dim_r, dim_s)
!!
!!    Square up and sort 1234 to 4132
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre
!!    and Andreas Skeidsvoll, 2018
!!
!!    Reorders the array x_pq_rs to x_sp_rq (i.e., 1234 to 4132).
!!
!!    The unordered array x_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s.
!!    The ordered array x_sp_rq is assumed allocated as dim_s*dim_p x dim_r*dim_q.
!!
      implicit none
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(:,:), intent(in) :: x_pqrs
      real(dp), dimension(dim_s*dim_p, dim_r*dim_q), intent(inout) :: x_sp_rq
!
      integer(i15) :: p, q, r, s, rs, pq, sp, rq, pqrs
!
!$omp parallel do schedule(static) private(s,r,rs,q,rq,p,pq,sp,pqrs)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
!
               rq = dim_r*(q-1) + r
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  sp = dim_s*(p-1) + s
                  pqrs = (max(pq,rs)*(max(pq,rs)-3)/2) + pq + rs
!
                  x_sp_rq(sp, rq) = x_pqrs(pqrs, 1)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine squareup_and_sort_1234_to_4132
!
!
   subroutine sort_1234_to_4123(x_pq_rs, x_sp_qr, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 4132
!!    Written by Eirik F. Kjønstad, Rolf H. Myhre and Andreas Skeidsvoll, 2018
!!
!!    Reorders the array x_pq_rs to x_sp_qr (i.e., 1234 to 4123).
!!
!!    The unordered array x_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s.
!!    The ordered array x_sp_qr is assumed allocated as dim_s*dim_p x dim_r*dim_q.
!!
      implicit none
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s), intent(in) :: x_pq_rs
      real(dp), dimension(dim_s*dim_p, dim_r*dim_q), intent(inout) :: x_sp_qr
!
      integer(i15) :: p, q, r, s, rs, pq, sp, qr
!
!$omp parallel do schedule(static) private(s,r,rs,q,qr,p,pq,sp)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
!
               qr = dim_q*(r-1) + q
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  sp = dim_s*(p-1) + s
!
                  x_sp_qr(sp, qr) = x_pq_rs(pq, rs)
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
   subroutine squareup_and_sort_1234_to_4123(x_pqrs, x_sp_qr, dim_p, dim_q, dim_r, dim_s)
!!
!!    Square up and srt 1234 to 4132
!!    Written by Eirik F. Kjønstad, Rolf H. Myhre and Andreas Skeidsvoll, 2018
!!
!!    Reorders the array x_pqrs to x_sp_qr (i.e., 1234 to 4123).
!!
!!    The unordered array x_pqrs is assumed allocated as dim_p*dim_q*dim_r*dim_s x 1.
!!    The ordered array x_sp_qr is assumed allocated as dim_s*dim_p x dim_r*dim_q.
!!
      implicit none
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(:,:), intent(in) :: x_pqrs
      real(dp), dimension(dim_s*dim_p, dim_r*dim_q), intent(inout) :: x_sp_qr
!
      integer(i15) :: p, q, r, s, rs, pq, sp, qr, pqrs
!
!$omp parallel do schedule(static) private(s,r,rs,q,qr,p,pq,sp,pqrs)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
!
               qr = dim_q*(r-1) + q
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  sp = dim_s*(p-1) + s
                  pqrs = (max(pq,rs)*(max(pq,rs)-3)/2) + pq + rs
!
                  x_sp_qr(sp, qr) = x_pqrs(pqrs, 1)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine squareup_and_sort_1234_to_4123
!
!
   subroutine sort_1234_to_3124(x_pq_rs, x_rp_qs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 3124
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pq_rs to x_rp_qs (i.e., 1234 to 3124).
!!
!!    The unordered array x_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s.
!!    The ordered array x_rp_qs is assumed allocated as dim_r*dim_p x dim_q*dim_s.
!!
      implicit none
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s), intent(in) :: x_pq_rs
      real(dp), dimension(dim_r*dim_p, dim_q*dim_s) :: x_rp_qs
!
      integer(i15) :: p, q, r, s, rs, pq, rp, qs
!
!$omp parallel do schedule(static) private(s,r,rs,q,rp,p,pq,qs)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
!
               qs = dim_q*(s-1) + q
!
               do p = 1, dim_p
!
                  rp = dim_r*(p-1) + r
                  pq = dim_p*(q-1) + p
!
                  x_rp_qs(rp, qs) = x_pq_rs(pq, rs)
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
   subroutine sort_1234_to_3142(x_pq_rs, x_rp_qs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 3142
!!    Written by Eirik F. Kjønstad, Rolf H. Myhre and
!!    Sarai D. Folkestad, 2018
!!
!!    Reorders the array x_pq_rs to x_rp_qs (i.e., 1234 to 3142).
!!
!!    The unordered array x_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s.
!!    The ordered array x_rp_qs is assumed allocated as dim_r*dim_p x dim_q*dim_s.
!!
      implicit none
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s), intent(in) :: x_pq_rs
      real(dp), dimension(dim_r*dim_p, dim_q*dim_s) :: x_rp_qs
!
      integer(i15) :: p, q, r, s, rs, pq, rp, sq
!
!$omp parallel do schedule(static) private(s,r,rs,q,rp,p,pq,sq)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
!
               sq = dim_s*(q-1) + s
!
               do p = 1, dim_p
!
                  rp = dim_r*(p-1) + r
                  pq = dim_p*(q-1) + p
!
                  x_rp_qs(rp, sq) = x_pq_rs(pq, rs)
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
   subroutine squareup_and_sort_1234_to_3124(x_pqrs, x_rp_qs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Squareup and sort 1234 to 3124
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pq_rs to x_rp_qs (i.e., 1234 to 3124).
!!
!!    The unordered array x_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s.
!!    The ordered array x_rp_qs is assumed allocated as dim_r*dim_p x dim_q*dim_s.
!!
      implicit none
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(:,:), intent(in) :: x_pqrs
      real(dp), dimension(dim_r*dim_p, dim_q*dim_s) :: x_rp_qs
!
      integer(i15) :: p, q, r, s, rs, pq, pqrs, rp, qs
!
!$omp parallel do schedule(static) private(s,r,rs,q,rp,p,pq,qs)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
!
               qs = dim_q*(s-1) + q
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  rp = dim_r*(p-1) + r
                  pqrs = (max(pq,rs)*(max(pq,rs)-3)/2) + pq + rs
!
                  x_rp_qs(rp, qs) = x_pqrs(pqrs, 1)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine squareup_and_sort_1234_to_3124
!
!
   subroutine squareup_and_sort_1234_to_2413(x_pqrs, x_qs_pr, dim_p, dim_q, dim_r, dim_s)
!!
!!    Squareup and sort 1234 to 2413
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pq_rs to x_qs_pr (i.e., 1234 to 2413).
!!
      implicit none
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(:,:), intent(in) :: x_pqrs
      real(dp), dimension(dim_q*dim_s, dim_p*dim_r) :: x_qs_pr
!
      integer(i15) :: p, q, r, s, rs, pq, pqrs, qs, pr
!
!$omp parallel do schedule(static) private(s,r,rs,q,pr,p,pq,qs)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
!
               qs = dim_q*(s-1) + q
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  pr = dim_p*(r-1) + p
                  pqrs = (max(pq,rs)*(max(pq,rs)-3)/2) + pq + rs
!
                  x_qs_pr(qs, pr) = x_pqrs(pqrs, 1)
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
   subroutine squareup_and_sort_1234_to_2341(x_pqrs, x_qs_pr, dim_p, dim_q, dim_r, dim_s)
!!
!!    Squareup and sort 1234 to 2341
!!    Written by Sarai D. Folkestad, 
!!    Eirik F. Kjønstad and Rolf H. Myhre, 2018
!!
!!    Reorders the array x_pq_rs to x_qr_sp (i.e., 1234 to 2413).
!!
      implicit none
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(:,:), intent(in) :: x_pqrs
      real(dp), dimension(dim_q*dim_s, dim_p*dim_r) :: x_qs_pr
!
      integer(i15) :: p, q, r, s, rs, pq, pqrs, qr, sp
!
!$omp parallel do schedule(static) private(s,r,rs,q,qr,p,pq,sp)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
!
               qr = dim_q*(r-1) + q
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  sp = dim_s*(p-1) + s
                  pqrs = (max(pq,rs)*(max(pq,rs)-3)/2) + pq + rs
!
                  x_qs_pr(qr, sp) = x_pqrs(pqrs, 1)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine squareup_and_sort_1234_to_2341
!
!
   subroutine sort_1234_to_2314(x_pq_rs, x_qr_ps, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 4132
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pq_rs to x_qr_ps (i.e., 1234 to 2314).
!!
!!    The unordered array x_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s.
!!    The ordered array x_qr_ps is assumed allocated as dim_q*dim_r x dim_p*dim_s.
!!
      implicit none
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s), intent(in) :: x_pq_rs
      real(dp), dimension(dim_q*dim_r, dim_p*dim_s) :: x_qr_ps
!
      integer(i15) :: p, q, r, s, rs, pq, qr, ps
!
!$omp parallel do schedule(static) private(s,r,rs,q,qr,p,pq,ps)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
!
               qr = dim_q*(r-1) + q
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  ps = dim_p*(s-1) + p
!
                  x_qr_ps(qr, ps) = x_pq_rs(pq, rs)
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
   subroutine sort_1234_to_2134(x_pq_rs, x_qr_ps, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 2134
!!    Written by Sarai D. Folkestad, 
!!    Eirik F. Kjønstad and Rolf H. Myhre, 2018
!!
!!    Reorders the array x_pq_rs to x_qp_rs (i.e., 1234 to 2314).
!!
!!    The unordered array x_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s.
!!    The ordered array x_qr_ps is assumed allocated as dim_q*dim_r x dim_p*dim_s.
!!
      implicit none
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s), intent(in) :: x_pq_rs
      real(dp), dimension(dim_q*dim_r, dim_p*dim_s) :: x_qr_ps
!
      integer(i15) :: p, q, r, s, rs, pq, qp
!
!$omp parallel do schedule(static) private(s,r,rs,q,qp,p,pq)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  qp = dim_q*(p-1) + q
!
                  x_qr_ps(qp, rs) = x_pq_rs(pq, rs)
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
   subroutine sort_1234_to_2413(x_pq_rs, x_qr_ps, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 2413
!!    Written by Sarai D. Folkestad, 
!!    Eirik F. Kjønstad and Rolf H. Myhre, 2018
!!
!!    Reorders the array x_pq_rs to x_qs_pr (i.e., 1234 to 2314).
!!
!!    The unordered array x_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s.
!!    The ordered array x_qr_ps is assumed allocated as dim_q*dim_r x dim_p*dim_s.
!!
      implicit none
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s), intent(in) :: x_pq_rs
      real(dp), dimension(dim_q*dim_r, dim_p*dim_s) :: x_qr_ps
!
      integer(i15) :: p, q, r, s, rs, pq, qs, pr
!
!$omp parallel do schedule(static) private(s,r,rs,q,qs,p,pq,pr)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
!
               qs = dim_q*(s-1) + q
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  pr = dim_p*(r-1) + p
!
                  x_qr_ps(qs, pr) = x_pq_rs(pq, rs)
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
   subroutine sort_1234_to_1324(x_pq_rs, x_pr_qs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 1324
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pq_rs to x_pr_qs (i.e., 1234 to 1324).
!!
!!    The unordered array x_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s.
!!    The ordered array x_pr_qs is assumed allocated as dim_p*dim_r x dim_q*dim_s.
!!
      implicit none
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s), intent(in) :: x_pq_rs
      real(dp), dimension(dim_p*dim_r, dim_q*dim_s) :: x_pr_qs
!
      integer(i15) :: p, q, r, s, rs, qs, pq, pr
!
!$omp parallel do schedule(static) private(s,r,rs,q,qs,p,pq,pr)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
!
               qs = dim_q*(s-1) + q
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  pr = dim_p*(r-1) + p
!
                  x_pr_qs(pr, qs) = x_pq_rs(pq, rs)
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
   subroutine squareup_and_sort_1234_to_1324(x_pqrs, x_pr_qs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Squareup and sort 1234 to 1324
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pq_rs to x_pr_qs (i.e., 1234 to 1324).
!!
!!    The unordered array x_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s.
!!    The ordered array x_pr_qs is assumed allocated as dim_p*dim_r x dim_q*dim_s.
!!
      implicit none
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(:,:), intent(in) :: x_pqrs
      real(dp), dimension(dim_p*dim_r, dim_q*dim_s) :: x_pr_qs
!
      integer(i15) :: p, q, r, s, rs, qs, pq, pr, pqrs
!
!$omp parallel do schedule(static) private(s,r,rs,q,qs,p,pq,pr,pqrs)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
!
               qs = dim_q*(s-1) + q
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  pr = dim_p*(r-1) + p
                  pqrs = (max(pq,rs)*(max(pq,rs)-3)/2) + pq + rs
!
                  x_pr_qs(pr, qs) = x_pqrs(pqrs, 1)
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
   subroutine sort_1234_to_2341(x_pq_rs, x_qr_sp, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 2341
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pq_rs to x_qr_sp (i.e., 1234 to 2341).
!!
!!    The unordered array x_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s.
!!    The ordered array x_qr_sp is assumed allocated as dim_q*dim_r x dim_s*dim_p.
!!
      implicit none
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s), intent(in) :: x_pq_rs
      real(dp), dimension(dim_q*dim_r, dim_s*dim_p) :: x_qr_sp
!
      integer(i15) :: p, q, r, s, pq, rs, qr, sp
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,qr,sp)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
!
               qr = dim_q*(r-1) + q
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  sp = dim_s*(p-1) + s
!
                  x_qr_sp(qr, sp) = x_pq_rs(pq, rs)
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
   subroutine sort_1234_to_1342(x_pq_rs, x_pr_sq, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 1342
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pq_rs to x_pr_sq (i.e., 1234 to 1342).
!!
!!    The unordered array x_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s.
!!    The ordered array x_pr_sq is assumed allocated as dim_p*dim_r x dim_s*dim_q.
!!
      implicit none
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s), intent(in) :: x_pq_rs
      real(dp), dimension(dim_p*dim_r, dim_s*dim_q) :: x_pr_sq
!
      integer(i15) :: p, q, r, s, pq, rs, pr, sq
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,pr,sq)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
!
               sq = dim_s*(q-1) + s
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  pr = dim_p*(r-1) + p
!
                  x_pr_sq(pr, sq) = x_pq_rs(pq, rs)
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
   subroutine sort_1234_to_1432(x_pq_rs, x_ps_rq, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 1432
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pq_rs to x_ps_rq (i.e., 1234 to 1432).
!!
!!    The unordered array x_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s.
!!    The ordered array x_ps_rq is assumed allocated as dim_p*dim_s x dim_r*dim_q.
!!
      implicit none
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s), intent(in)    :: x_pq_rs
      real(dp), dimension(dim_p*dim_s, dim_r*dim_q), intent(inout) :: x_ps_rq
!
      integer(i15) :: p, q, r, s, pq, rs, ps, rq
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,ps,rq)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
!
               rq = dim_r*(q-1) + r
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  ps = dim_p*(s-1) + p
!
                  x_ps_rq(ps, rq) = x_pq_rs(pq, rs)
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
   subroutine squareup_and_sort_1234_to_4312(x_pqrs, x_sr_pq, dim_p, dim_q, dim_r, dim_s)
!!
!!    Square up and sort 1234 to 4312
!!    Written by Eirik F. Kjønstad, Rolf H. Myhre
!!    and Andreas Skeidsvoll, 2018
!!
!!    Squares up and reorders the array x_pqrs to x_sr_pq (i.e., 1234 to 4312).
!!
!!    The unordered array x_pqrs is assumed allocated as dim_p*dim_q*dim_r*dim_s x 1.
!!    The ordered array x_sr_pq is assumed allocated as dim_s*dim_r x dim_p*dim_q.
!!
      implicit none
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(:,:), intent(in) :: x_pqrs
      real(dp), dimension(dim_s*dim_r, dim_p*dim_q), intent(inout) :: x_sr_pq
!
      integer(i15) :: p, q, r, s, pq, rs, sr, pqrs
!
!$omp parallel do schedule(static) private(s,r,rs,sr,q,p,pq,pqrs)
      do s = 1, dim_s
            do r = 1, dim_r
   !
               rs = dim_r*(s-1) + r
               sr = dim_s*(r-1) + s
   !
               do q = 1, dim_q
                  do p = 1, dim_p
   !
                     pq = dim_p*(q-1) + p
                     pqrs = (max(pq,rs)*(max(pq,rs)-3)/2) + pq + rs
   !
                     x_sr_pq(sr, pq) = x_pqrs(pqrs, 1)
   !
                  enddo
               enddo
            enddo
         enddo
!$omp end parallel do
!
   end subroutine squareup_and_sort_1234_to_4312
!
!
   subroutine sort_1234_to_4312(x_pq_rs, x_sr_pq, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 4312
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pq_rs to x_sr_pq (i.e., 1234 to 4312).
!!
!!    The unordered array x_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s.
!!    The ordered array x_sr_pq is assumed allocated as dim_s*dim_r x dim_p*dim_q.
!!
      implicit none
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s), intent(in)    :: x_pq_rs
      real(dp), dimension(dim_s*dim_r, dim_p*dim_q), intent(inout) :: x_sr_pq
!
      integer(i15) :: p, q, r, s, pq, rs, sr
!
!$omp parallel do schedule(static) private(s,r,rs,sr)
      do s = 1, dim_s
         do r = 1, dim_r
!
            sr = dim_s*(r-1) + s
            rs = dim_r*(s-1) + r
!
            x_sr_pq(sr, :) = x_pq_rs(:, rs)
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_4312
!
!
   subroutine sort_1234_to_1423(x_pq_rs, x_ps_qr, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 1423
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pq_rs to x_ps_qr (i.e., 1234 to 1423).
!!
!!    The unordered array x_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s.
!!    The ordered array x_ps_qr is assumed allocated as dim_p*dim_s x dim_q*dim_r.
!!
      implicit none
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s), intent(in)    :: x_pq_rs
      real(dp), dimension(dim_p*dim_s, dim_r*dim_q), intent(inout) :: x_ps_qr
!
      integer(i15) :: p, q, r, s, pq, rs, ps, qr
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,ps,qr)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
!
               qr = dim_q*(r-1) + q
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  ps = dim_p*(s-1) + p
!
                  x_ps_qr(ps, qr) = x_pq_rs(pq, rs)
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
   subroutine add_1423_to_1234(gamma, x, y_pq_rs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 1423 to 1234
!!    Written by Eirik F. Kjønstad, Rolf H. Myhre and Andreas Skeidsvoll, 2018
!!
!!    Performs:
!!
!!       y_pq_rs(pq,rs) = y_pq_rs(pq,rs) + gamma * x(ps, qr)
!!
!!    The unordered array y_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s,
!!    and x accordingly.
!!
      implicit none
!
      real(dp), intent(in) :: gamma
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s) :: y_pq_rs
      real(dp), dimension(dim_p*dim_s, dim_q*dim_r), intent(in) :: x
!
      integer(i15) :: p, q, r, s, pq, rs, ps, qr
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,ps,qr)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
!
               qr = dim_q*(r-1) + q
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  ps = dim_p*(s-1) + p
!
                  y_pq_rs(pq, rs) = y_pq_rs(pq, rs) + gamma*x(ps, qr)
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
   subroutine add_1432_to_1234(gamma, x, y_pq_rs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 1432 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pq_rs(pq,rs) = y_pq_rs(pq,rs) + gamma * x(ps, rq)
!!
!!    The unordered array y_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s,
!!    and x accordingly.
!!
      implicit none
!
      real(dp), intent(in) :: gamma
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s) :: y_pq_rs
      real(dp), dimension(dim_p*dim_s, dim_r*dim_q), intent(in) :: x
!
      integer(i15) :: p, q, r, s, pq, rs, ps, rq
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,ps,rq)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
!
               rq = dim_r*(q-1) + r
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  ps = dim_p*(s-1) + p
!
                  y_pq_rs(pq, rs) = y_pq_rs(pq, rs) + gamma*x(ps, rq)
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
   subroutine add_1342_to_1234(gamma, x, y_pq_rs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 1342 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pq_rs(pq,rs) = y_pq_rs(pq,rs) + gamma * x(pr, sq)
!!
!!    The unordered array y_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s,
!!    and x accordingly.
!!
      implicit none
!
      real(dp), intent(in) :: gamma
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s) :: y_pq_rs
      real(dp), dimension(dim_p*dim_r, dim_s*dim_q), intent(in) :: x
!
      integer(i15) :: p, q, r, s, pq, rs, pr, sq
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,pr,sq)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
!
               sq = dim_s*(q-1) + s
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  pr = dim_p*(r-1) + p
!
                  y_pq_rs(pq, rs) = y_pq_rs(pq, rs) + gamma*x(pr, sq)
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
   subroutine add_1243_to_1234(gamma, x, y_pq_rs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 1243 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pq_rs(pq,rs) = y_pq_rs(pq,rs) + gamma * x(pq, sr)
!!
!!    The unordered array y_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s,
!!    and x accordingly.
!!
      implicit none
!
      real(dp), intent(in) :: gamma
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s) :: y_pq_rs
      real(dp), dimension(dim_p*dim_q, dim_s*dim_r), intent(in) :: x
!
      integer(i15) :: p, q, r, s, pq, rs, sr
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,sr)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
            sr = dim_s*(r-1) + s
!
            do q = 1, dim_q
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
!
                  y_pq_rs(pq, rs) = y_pq_rs(pq, rs) + gamma*x(pq, sr)
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
   subroutine add_3412_to_1234(gamma, x, y_pq_rs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 3412 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pq_rs(pq,rs) = y_pq_rs(pq,rs) + gamma * x(rs, pq)
!!
!!    The unordered array y_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s,
!!    and x accordingly.
!!
      implicit none
!
      real(dp), intent(in) :: gamma
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s) :: y_pq_rs
      real(dp), dimension(dim_r*dim_s, dim_p*dim_q), intent(in) :: x
!
      integer(i15) :: p, q, r, s, pq, rs
!
      call add_21_to_12(gamma, x, y_pq_rs, dim_p*dim_q, dim_r*dim_s)
!
   end subroutine add_3412_to_1234
!
!
   subroutine add_3421_to_1234(gamma, x, y_pq_rs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 3421 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pq_rs(pq,rs) = y_pq_rs(pq,rs) + gamma * x(rs, qp)
!!
!!    The unordered array y_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s,
!!    and x accordingly.
!!
      implicit none
!
      real(dp), intent(in) :: gamma
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s) :: y_pq_rs
      real(dp), dimension(dim_r*dim_s, dim_p*dim_q), intent(in) :: x
!
      integer(i15) :: p, q, I, pq, qp
!
!$omp parallel do schedule(static) private(I,q,p,pq,qp)
      do I = 1, dim_r*dim_s
         do q = 1, dim_q
            do p = 1, dim_p
!
               pq = dim_p*(q-1) + p
               qp = dim_q*(p-1) + q
!
               y_pq_rs(pq,I) = y_pq_rs(pq,I) + gamma*x(I,qp)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_3421_to_1234
!
!
   subroutine add_2341_to_1234(gamma, x, y_pq_rs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 2341 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pq_rs(pq,rs) = y_pq_rs(pq,rs) + gamma * x(qr, sp)
!!
!!    The unordered array y_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s,
!!    and x accordingly.
!!
      implicit none
!
      real(dp), intent(in) :: gamma
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s) :: y_pq_rs
      real(dp), dimension(dim_q*dim_r, dim_s*dim_p), intent(in) :: x
!
      integer(i15) :: p, q, r, s, pq, rs, qr, sp
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,sp,qr)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
!
               qr = dim_q*(r-1) + q
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  sp = dim_s*(p-1) + s
!
                  y_pq_rs(pq, rs) = y_pq_rs(pq, rs) + gamma*x(qr, sp)
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
   subroutine add_2143_to_1234(gamma, x, y_pq_rs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 2143 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pq_rs(pq,rs) = y_pq_rs(pq,rs) + gamma * x(qp, sr)
!!
!!    The unordered array y_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s,
!!    and x accordingly.
!!
      implicit none
!
      real(dp), intent(in) :: gamma
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s) :: y_pq_rs
      real(dp), dimension(dim_q*dim_p, dim_s*dim_r), intent(in) :: x
!
      integer(i15) :: p, q, r, s, pq, rs, qp, sr
!
!$omp parallel do schedule(static) private(s,r,q,p,qp,rs,sr,pq)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
            sr = dim_s*(r-1) + s
!
            do q = 1, dim_q
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  qp = dim_q*(p-1) + q
!
                  y_pq_rs(pq, rs) = y_pq_rs(pq, rs) + gamma*x(qp, sr)
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
   subroutine add_3214_to_1234(gamma, x, y_pq_rs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 3214 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pq_rs(pq,rs) = y_pq_rs(pq,rs) + gamma * x(rq, ps)
!!
!!    The unordered array y_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s,
!!    and x accordingly.
!!
      implicit none
!
      real(dp), intent(in) :: gamma
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s) :: y_pq_rs
      real(dp), dimension(dim_r*dim_q, dim_p*dim_s), intent(in) :: x
!
      integer(i15) :: p, q, r, s, pq, rs, rq, ps
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,ps,rq)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
!
               rq = dim_r*(q-1) + r
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  ps = dim_p*(s-1) + p
!
                  y_pq_rs(pq, rs) = y_pq_rs(pq, rs) + gamma*x(rq, ps)
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
   subroutine add_4231_to_1234(gamma, x, y_pq_rs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 4231 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pq_rs(pq,rs) = y_pq_rs(pq,rs) + gamma * x(sq, rp)
!!
!!    The unordered array y_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s,
!!    and x accordingly.
!!
      implicit none
!
      real(dp), intent(in) :: gamma
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s) :: y_pq_rs
      real(dp), dimension(dim_s*dim_q, dim_r*dim_p), intent(in) :: x
!
      integer(i15) :: p, q, r, s, pq, rs, sq, rp
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,rp,sq)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
!
               sq = dim_s*(q-1) + s
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  rp = dim_r*(p-1) + r
!
                  y_pq_rs(pq, rs) = y_pq_rs(pq, rs) + gamma*x(sq, rp)
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
   subroutine add_4213_to_1234(gamma, x, y_pq_rs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 4213 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pq_rs(pq,rs) = y_pq_rs(pq,rs) + gamma * x(sq, pr)
!!
!!    The unordered array y_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s,
!!    and x accordingly.
!!
      implicit none
!
      real(dp), intent(in) :: gamma
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s) :: y_pq_rs
      real(dp), dimension(dim_s*dim_q, dim_p*dim_r), intent(in) :: x
!
      integer(i15) :: p, q, r, s, pq, rs, sq, pr
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,pr,sq)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
!
               sq = dim_s*(q-1) + s
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  pr = dim_p*(r-1) + p
!
                  y_pq_rs(pq, rs) = y_pq_rs(pq, rs) + gamma*x(sq, pr)
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
   subroutine add_4321_to_1234(gamma, x, y_pq_rs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 4321 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pq_rs(pq,rs) = y_pq_rs(pq,rs) + gamma * x(sr, qp)
!!
!!    The unordered array y_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s,
!!    and x accordingly.
!!
      implicit none
!
      real(dp), intent(in) :: gamma
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s) :: y_pq_rs
      real(dp), dimension(dim_s*dim_r, dim_q*dim_p), intent(in) :: x
!
      integer(i15) :: p, q, r, s, pq, rs, sr, qp
!
!$omp parallel do schedule(static) private(s,r,q,p,qp,rs,sr,pq)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
            sr = dim_s*(r-1) + s
!
            do q = 1, dim_q
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  qp = dim_q*(p-1) + q
!
                  y_pq_rs(pq, rs) = y_pq_rs(pq, rs) + gamma*x(sr, qp)
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
   subroutine sort_1234_to_4321(x_pq_rs, x_sr_qp, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 4321
!!    Written by Eirik F. Kjønstad, Rolf H. Myhre and Andreas Skeidsvoll, 2018
!!
!!    Reorders the array x_pq_rs to x_sr_qp (i.e., 1234 to 4321).
!!
!!    The unordered array x_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s.
!!    The ordered array x_sr_qp is assumed allocated as dim_s*dim_r x dim_q*dim_p.
!!
      implicit none
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s), intent(in)    :: x_pq_rs
      real(dp), dimension(dim_s*dim_r, dim_q*dim_p), intent(inout) :: x_sr_qp
!
      integer(i15) :: p, q, r, s, qp, rs, sr, pq
!
!$omp parallel do schedule(static) private(s,r,q,p,qp,rs,sr,pq)
      do s = 1, dim_s
            do r = 1, dim_r
   !
               rs = dim_r*(s-1) + r
               sr = dim_s*(r-1) + s
   !
               do q = 1, dim_q
                  do p = 1, dim_p
   !
                     pq = dim_p*(q-1) + p
                     qp = dim_q*(p-1) + q
   !
                     x_sr_qp(sr, qp) = x_pq_rs(pq, rs)
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
   subroutine add_4312_to_1234(gamma, x, y_pq_rs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 4312 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pq_rs(pq,rs) = y_pq_rs(pq,rs) + gamma * x(sr, pq)
!!
!!    The unordered array y_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s,
!!    and x accordingly.
!!
      implicit none
!
      real(dp), intent(in) :: gamma
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s) :: y_pq_rs
      real(dp), dimension(dim_s*dim_r, dim_p*dim_q), intent(in) :: x
!
      integer(i15) :: p, q, r, s, pq, rs, sr
!
!$omp parallel do schedule(static) private(s,r,q,p,sr,rs,pq)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
            sr = dim_s*(r-1) + s
!
            do q = 1, dim_q
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
!
                  y_pq_rs(pq, rs) = y_pq_rs(pq, rs) + gamma*x(sr, pq)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_4312_to_1234
!
!
   subroutine add_4123_to_1234(gamma, x, y_pq_rs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 4123 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pq_rs(pq,rs) = y_pq_rs(pq,rs) + gamma * x(sp, qr)
!!
!!    The unordered array y_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s,
!!    and x accordingly.
!!
      implicit none
!
      real(dp), intent(in) :: gamma
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s) :: y_pq_rs
      real(dp), dimension(dim_s*dim_p, dim_q*dim_r), intent(in) :: x
!
      integer(i15) :: p, q, r, s, pq, rs, sp, qr
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,qr,sp)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
!
               qr = dim_q*(r-1) + q
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  sp = dim_s*(p-1) + s
!
                  y_pq_rs(pq, rs) = y_pq_rs(pq, rs) + gamma*x(sp, qr)
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
   subroutine add_4132_to_1234(gamma, x, y_pq_rs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add 4132 to 1234
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pq_rs(pq,rs) = y_pq_rs(pq,rs) + gamma * x(sp, rq)
!!
!!    The unordered array y_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s,
!!    and x accordingly.
!!
      implicit none
!
      real(dp), intent(in) :: gamma
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s) :: y_pq_rs
      real(dp), dimension(dim_s*dim_p, dim_q*dim_r), intent(in) :: x
!
      integer(i15) :: p, q, r, s, pq, rs, sp, rq
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,rq,sp)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
!
               rq = dim_r*(q-1) + r
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  sp = dim_s*(p-1) + s
!
                  y_pq_rs(pq, rs) = y_pq_rs(pq, rs) + gamma*x(sp, rq)
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
   subroutine squareup_and_sort_1234_to_1432(x_pqrs, x_ps_rq, dim_p, dim_q, dim_r, dim_s)
!!
!!    Squareup and sort 1234 to 1432
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the packed array x_pq_rs to unpacked x_ps_rq (i.e., 1234 to 1432).
!!
!!    The unordered packed array x_pqrs is assumed allocated as an (pq,rs) packed x 1 array.
!!    This array will typically be the packed T2 amplitudes (wf%t2am).
!!
!!    The ordered array x_ps_rq is assumed allocated as dim_p*dim_s x dim_r*dim_q.
!!
      implicit none
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(:,:), intent(in) :: x_pqrs
      real(dp), dimension(dim_p*dim_s, dim_r*dim_q) :: x_ps_rq
!
      integer(i15) :: p, q, r, s, pq, rs, pqrs, ps, rq
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,ps,rq,pqrs)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
!
               rq = dim_r*(q-1) + r
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  ps = dim_p*(s-1) + p
                  pqrs = (max(pq,rs)*(max(pq,rs)-3)/2) + pq + rs
!
                  x_ps_rq(ps, rq) = x_pqrs(pqrs, 1)
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
   subroutine squareup_and_sort_1234_to_1243(x_pqrs, x_pq_sr, dim_p, dim_q, dim_r, dim_s)
!!
!!    Squareup and sort 1234 to 1243
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the packed array x_pq_rs to unpacked x_pq_sr (i.e., 1234 to 1432).
!!
!!    The unordered packed array x_pqrs is assumed allocated as an (pq,rs) packed x 1 array.
!!    This array will typically be the packed T2 amplitudes (wf%t2am).
!!
!!    The ordered array x_pq_sr is assumed allocated as dim_p*dim_q x dim_s*dim_r.
!!
      implicit none
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(:,:), intent(in) :: x_pqrs
      real(dp), dimension(dim_p*dim_q, dim_s*dim_r) :: x_pq_sr
!
      integer(i15) :: p, q, r, s, pq, rs, pqrs, sr
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,sr,pqrs)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
            sr = dim_s*(r-1) + s
!
            do q = 1, dim_q
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  pqrs = (max(pq,rs)*(max(pq,rs)-3)/2) + pq + rs
!
                  x_pq_sr(pq, sr) = x_pqrs(pqrs, 1)
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
   subroutine squareup_and_sort_1234_to_1423(x_pqrs, x_ps_qr, dim_p, dim_q, dim_r, dim_s)
!!
!!    Squareup and sort 1234 to 1423
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the packed array x_pq_rs to unpacked x_ps_qr (i.e., 1234 to 1423).
!!
!!    The unordered packed array x_pqrs is assumed allocated as an (pq,rs) packed x 1 array.
!!    This array will typically be the packed T2 amplitudes (wf%t2am).
!!
!!    The ordered array x_ps_qr is assumed allocated as dim_p*dim_s x dim_r*dim_q.
!!
      implicit none
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(:,:), intent(in) :: x_pqrs
      real(dp), dimension(dim_p*dim_s, dim_r*dim_q) :: x_ps_qr
!
      integer(i15) :: p, q, r, s, pq, rs, pqrs, ps, qr
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,ps,qr,pqrs)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
!
               qr = dim_q*(r-1) + q
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  ps = dim_p*(s-1) + p
                  pqrs = (max(pq,rs)*(max(pq,rs)-3)/2) + pq + rs
!
                  x_ps_qr(ps, qr) = x_pqrs(pqrs, 1)
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
   subroutine sort_1234_to_3214(x_pq_rs, x_rq_ps, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 3214
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pq_rs to x_rq_ps (i.e., 1234 to 3214).
!!
!!    The unordered array x_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s.
!!    The ordered array x_rq_ps is assumed allocated as dim_r*dim_q x dim_p*dim_s.
!!
      implicit none
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s), intent(in) :: x_pq_rs
      real(dp), dimension(dim_r*dim_q, dim_p*dim_s) :: x_rq_ps
!
      integer(i15) :: p, q, r, s, pq, rs, rq, ps
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,ps,rq)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
!
               rq = dim_r*(q-1) + r
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  ps = dim_p*(s-1) + p
!
                  x_rq_ps(rq, ps) = x_pq_rs(pq, rs)
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
   subroutine sort_1234_to_4231(x_pq_rs, x_sq_rp, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 4231
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pq_rs to x_sq_rp (i.e., 1234 to 4231).
!!
!!    The unordered array x_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s.
!!    The ordered array x_sq_rp is assumed allocated as dim_s*dim_q x dim_r*dim_p.
!!
      implicit none
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s), intent(in) :: x_pq_rs
      real(dp), dimension(dim_s*dim_q, dim_r*dim_p) :: x_sq_rp
!
      integer(i15) :: p, q, r, s, pq, rs, sq, rp
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,sq,rp)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
!
               sq = dim_s*(q-1) + s
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  rp = dim_r*(p-1) + r
!
                  x_sq_rp(sq, rp) = x_pq_rs(pq, rs)
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
   subroutine sort_1234_to_4213(x_pq_rs, x_sq_pr, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 4213
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pq_rs to x_sq_pr (i.e., 1234 to 4213).
!!
!!    The unordered array x_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s.
!!    The ordered array x_sq_pr is assumed allocated as dim_s*dim_q x dim_p*dim_r.
!!
      implicit none
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s), intent(in) :: x_pq_rs
      real(dp), dimension(dim_s*dim_q, dim_r*dim_p) :: x_sq_pr
!
      integer(i15) :: p, q, r, s, pq, rs, sq, pr
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,sq,pr)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
!
               sq = dim_s*(q-1) + s
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  pr = dim_p*(r-1) + p
!
                  x_sq_pr(sq, pr) = x_pq_rs(pq, rs)
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
   subroutine squareup_and_sort_1234_to_4213(x_pqrs, x_sq_pr, dim_p, dim_q, dim_r, dim_s)
!!
!!    Square up and sort 1234 to 4213
!!    Written by Eirik F. Kjønstad, Rolf H. Myhre
!!    and Andreas Skeidsvoll, 2018
!!
!!    Reorders the array x_pq_rs to x_sq_pr (i.e., 1234 to 4213).
!!
!!    The unordered array x_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s.
!!    The ordered array x_sq_pr is assumed allocated as dim_s*dim_q x dim_p*dim_r.
!!
      implicit none
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(:,:), intent(in) :: x_pqrs
      real(dp), dimension(dim_s*dim_q, dim_r*dim_p) :: x_sq_pr
!
      integer(i15) :: p, q, r, s, pq, rs, sq, pr, pqrs
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,sq,pr,pqrs)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
!
               sq = dim_s*(q-1) + s
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  pr = dim_p*(r-1) + p
                  pqrs = (max(pq,rs)*(max(pq,rs)-3)/2) + pq + rs
!
                  x_sq_pr(sq, pr) = x_pqrs(pqrs,1)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine squareup_and_sort_1234_to_4213
!
!
   subroutine squareup_and_sort_1234_to_3214(x_pqrs, x_rq_ps, dim_p, dim_q, dim_r, dim_s)
!!
!!    Squareup and sort 1234 to 3214
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the packed array x_pq_rs to unpacked x_rq_ps (i.e., 1234 to 3214).
!!
!!    The unordered packed array x_pqrs is assumed allocated as an (pq,rs) packed x 1 array.
!!    This array will typically be the packed T2 amplitudes (wf%t2am).
!!
!!    The ordered array x_rq_ps is assumed allocated as dim_r*dim_q x dim_p*dim_s.
!!
      implicit none
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(:, :), intent(in) :: x_pqrs
      real(dp), dimension(dim_r*dim_q, dim_p*dim_s) :: x_rq_ps
!
      integer(i15) :: p, q, r, s, pq, rs, pqrs, rq, ps
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,ps,rq,pqrs)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
!
               rq = dim_r*(q-1) + r
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  ps = dim_p*(s-1) + p
                  pqrs = (max(pq,rs)*(max(pq,rs)-3)/2) + pq + rs
!
                  x_rq_ps(rq, ps) = x_pqrs(pqrs, 1)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine squareup_and_sort_1234_to_3214
!
!
   subroutine sort_1234_to_3241(x_pq_rs, x_rq_sp, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 3241
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pq_rs to x_rq_sp (i.e., 1234 to 3241).
!!
!!    The unordered array x_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s.
!!    The ordered array x_rq_sp is assumed allocated as dim_r*dim_q x dim_s*dim_p.
!!
      implicit none
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s), intent(in) :: x_pq_rs
      real(dp), dimension(dim_r*dim_q, dim_s*dim_p) :: x_rq_sp
!
      integer(i15) :: p, q, r, s, pq, rs, rq, sp
!
!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,rq,sp)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
!
            do q = 1, dim_q
!
               rq = dim_r*(q-1) + r
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  sp = dim_s*(p-1) + s
!
                  x_rq_sp(rq, sp) = x_pq_rs(pq, rs)
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
   subroutine sort_1234_to_1243(x_pq_rs, x_pq_sr, dim_p, dim_q, dim_r, dim_s)
!!
!!    Sort 1234 to 1243
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Reorders the array x_pq_rs to x_pq_sr (i.e., 1234 to 1243).
!!
!!    The unordered array x_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s.
!!    The ordered array x_pq_sr is assumed allocated as dim_p*dim_q x dim_s*dim_r.
!!
      implicit none
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s
!
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s), intent(in) :: x_pq_rs
      real(dp), dimension(dim_p*dim_q, dim_s*dim_r) :: x_pq_sr
!
      integer(i15) :: I, r, s, sr, rs
!
!$omp parallel do schedule(static) private(I,r,s,sr,rs)
      do I = 1, dim_p*dim_q
!
         do s = 1, dim_s
            do r = 1, dim_r
!
               sr = dim_s*(r-1) + s
               rs = dim_r*(s-1) + r
!
               x_pq_sr(I, sr) = x_pq_rs(I, rs)
!
            enddo
         enddo
!
      enddo
!$omp end parallel do
!
   end subroutine sort_1234_to_1243
!
!
!     -::- Squareup and packin and related routines -::-
!     --------------------------------------------------
!
   subroutine add_to_packed(packed, unpacked, N)
!!
!!    Add to packed
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2017
!!
!!    Adds a symmetric unpacked N x N matrix to a packed N x N
!!    matrix.
!!
      implicit none
!
      integer(i15), intent(in) :: N
!
      real(dp), dimension(:,:), intent(inout) :: packed
      real(dp), dimension(:,:), intent(in) :: unpacked
!
      integer(i15) :: i, j, ij
!
!$omp parallel do schedule(static) private(i,j, ij)
      do i = 1, N
         do j = i, N
            ij = (max(i,j)*(max(i,j)-3)/2) + i + j
            packed(ij, 1) = packed(ij, 1) + unpacked(i,j)
         enddo
      enddo
!$omp end parallel do
!
   end subroutine add_to_packed
!
!
   subroutine symmetrize_and_add_to_packed(packed, unpacked, N)
!!
!!    Symmetrize and add to packed
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2017
!!
!!    Symmetrizes unpacked N x N matrix and adds the symmtrized sum
!!    to a packed N x N matrix.
!!
      implicit none
!
      integer(i15), intent(in) :: N
!
      real(dp), dimension(:,:) :: packed
      real(dp), dimension(:,:), intent(in) :: unpacked
!
      integer(i15) :: i, j, ij
!
!$omp parallel do schedule(static) private(i,j, ij)
      do i = 1, N
         do j = i, N
            ij = (max(i,j)*(max(i,j)-3)/2) + i + j
            packed(ij, 1) = packed(ij, 1) + unpacked(i,j) + unpacked(j,i)
         enddo
      enddo
!$omp end parallel do
!
   end subroutine symmetrize_and_add_to_packed
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
      integer(i15), intent(in) :: N
!
      real(dp), dimension(:,:), intent(in) :: packed
      real(dp), dimension(:,:)             :: unpacked
!
      integer(i15) :: i = 0, j = 0
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
            unpacked(i, j) = packed(index_packed(i - 1, j), 1)
            unpacked(j, i) = -unpacked(i, j)
!
         enddo
      enddo
!
   end subroutine squareup_anti
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
!$omp parallel do schedule(static) private(i,j)
      do i = 1, N
         do j = 1, N
            unpacked(i, j) = packed(index_packed(i,j), 1)
         enddo
      enddo
!$omp end parallel do
!
   end subroutine squareup
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
      integer(i15) :: i = 0, j = 0, ij = 0
!
!$omp parallel do schedule(static) private(i,j, ij)
      do i = 1, N
         do j = 1, N
            ij = N*(j - 1) + i
            unpacked(ij, :) = packed(index_packed(i,j), :)
         enddo
      enddo
!$omp end parallel do
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
!$omp parallel do schedule(static) private(i,j)
      do i = 1, N
         do j = 1, i
!
            packed(index_packed(i, j), 1) = unpacked(i, j)
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine packin
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
      integer(i15), intent(in) :: N
!
      real(dp), dimension(:,:) :: packed
      real(dp), dimension(:,:), intent(in) :: unpacked
!
      integer(i15) :: i = 0, j = 0
!
      do i = 2, N
         do j = 1, i - 1
!
            packed(index_packed(i - 1, j), 1) = unpacked(i, j)
!
         enddo
      enddo
!
   end subroutine packin_anti
!
!
end module reordering
