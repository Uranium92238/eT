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
!  ::::::::::::::::::::::::::
!  -::- Sorting routines -::-
!  ::::::::::::::::::::::::::
!
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
   pure subroutine sort_1234_to_4132(x_pq_rs, x_sp_rq, dim_p, dim_q, dim_r, dim_s)
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
!!$omp parallel do schedule(static) private(s,r,rs,q,rq,p,pq,sp)
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
!!$omp end parallel do
!
   end subroutine sort_1234_to_4132
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
!!$omp parallel do schedule(static) private(s,r,rs,q,rq,p,pq,sp)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r 
!
            do q = 1, dim_q
!
               rp = dim_r*(p-1) + r
!
               do p = 1, dim_p 
!
                  pq = dim_p*(q-1) + p 
                  qs = dim_q*(s-1) + q 
!
                  x_rp_qs(rp, qs) = x_pq_rs(pq, rs)
!
               enddo
            enddo
         enddo
      enddo
!!$omp end parallel do
!
   end subroutine sort_1234_to_3124
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
!!$omp parallel do schedule(static) private(s,r,rs,q,rq,p,pq,sp)
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
!!$omp end parallel do
!
   end subroutine squareup_and_sort_1234_to_3124
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
!!$omp parallel do schedule(static) private(s,r,rs,q,rq,p,pq,sp)
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
!!$omp end parallel do
!
   end subroutine sort_1234_to_2314
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
!!$omp parallel do schedule(static) private(s,r,rs,q,qs,p,pq,pr)
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
!!$omp end parallel do
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
!!$omp parallel do schedule(static) private(s,r,rs,q,qs,p,pq,pr,pqrs)
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
!!$omp end parallel do
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
!!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,qr,sp)
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
!!$omp end parallel do
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
!!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,qr,sp)
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
!!$omp end parallel do
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
!!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,ps,rq)
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
!!$omp end parallel do
!
   end subroutine sort_1234_to_1432
!
!
   subroutine add_sorted_1234_to_1432(delta, y_pq_rs, gamma, x_pq_rs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add sorted 1234 to 1432
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pq_rs(pq,rs) = delta*y_pq_rs(pq,rs) + gamma * x_pq_rs(ps, rq)
!!
!!    As it adds the exchange contribution, dim_q = dim_s by assumption.
!!
!!    The unordered array x_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s.
!!    The unordered array y_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s.
!!
      implicit none 
!
      real(dp), intent(in) :: gamma, delta 
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s 
!
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s) :: y_pq_rs 
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s), intent(in) :: x_pq_rs
!
      integer(i15) :: p, q, r, s, pq, rs, ps, rq
!
   !   call dscal(dim_p*dim_q*dim_r*dim_s, delta, y_pq_rs, 1)
!
!!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,ps,rq)
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
                  y_pq_rs(pq, rs) = delta*y_pq_rs(pq, rs) + gamma*x_pq_rs(ps, rq)
!
               enddo
            enddo
         enddo
      enddo
!!$omp end parallel do
!
   end subroutine add_sorted_1234_to_1432
!
!
   subroutine add_sorted_1234_to_4231(delta, y_pq_rs, gamma, x_pq_rs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Add sorted 1234 to 4231
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    Performs:
!!
!!       y_pq_rs(pq,rs) = delta*y_pq_rs(pq,rs) + gamma * x_pq_rs(sq, rp)
!!
!!    As it adds the exchange contribution, dim_q = dim_s by assumption.
!!
!!    The unordered array x_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s.
!!    The unordered array y_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s.
!!
      implicit none 
!
      real(dp), intent(in) :: gamma, delta 
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s 
!
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s) :: y_pq_rs 
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s), intent(in) :: x_pq_rs
!
      integer(i15) :: p, q, r, s, pq, rs, sq, rp 
!
   !   call dscal(dim_p*dim_q*dim_r*dim_s, delta, y_pq_rs, 1)
!
!!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,ps,rq)
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
                  y_pq_rs(pq, rs) = delta*y_pq_rs(pq, rs) + gamma*x_pq_rs(sq, rp)
!
               enddo
            enddo
         enddo
      enddo
!!$omp end parallel do
!
   end subroutine add_sorted_1234_to_4231
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
!!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,ps,rq,pqrs)
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
!!$omp end parallel do
!
   end subroutine squareup_and_sort_1234_to_1432
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
!!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,ps,rq)
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
!!$omp end parallel do
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
!!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,ps,rq)
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
!!$omp end parallel do
!
   end subroutine sort_1234_to_4231
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
!!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,ps,rq,pqrs)
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
!!$omp end parallel do
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
!!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,rq,sp)
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
!!$omp end parallel do
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
!!$omp parallel do schedule(static) private(I,r,s,sr,rs)
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
!!$omp end parallel do
!
   end subroutine sort_1234_to_1243
!
!
   subroutine set_two_1234_minus_1432(x_pq_rs, y_pq_rs, dim_p, dim_q, dim_r, dim_s)
!!
!!    Set two 1234 minus 1432
!!    Written by Eirik F. Kjønstad and Rolf H. Myhre, Dec 2017
!!
!!    A two-Coulomb minus exchange routine. More precisely,
!!    the routine sets 
!!
!!       y_pq_rs(pq,rs) = 2*x_pq_rs(pq,rs) - x_pq_rs(ps,rq),
!!
!!    where the array x_pq_rs is assumed allocated as dim_p*dim_q x dim_r*dim_s,
!!    with dim_q = dim_s 
!!
      implicit none 
!
      integer(i15), intent(in) :: dim_p, dim_q, dim_r, dim_s 
!
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s), intent(in) :: x_pq_rs 
      real(dp), dimension(dim_p*dim_q, dim_r*dim_s) :: y_pq_rs
!
      integer(i15) :: p, q, r, s, pq, rs, ps, rq
!
!!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,ps,rq)
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
                  y_pq_rs(pq, rs) = two*x_pq_rs(pq, rs)-x_pq_rs(ps, rq)
!
               enddo
            enddo
         enddo
      enddo
!!$omp end parallel do
!
   end subroutine set_two_1234_minus_1432
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
      real(dp), dimension(:,:) :: packed
      real(dp), dimension(:,:), intent(in) :: unpacked
!
      integer(i15) :: i, j, ij 
!
!!$omp parallel do schedule(static) private(i,j)
      do i = 1, N
         do j = i, N
            ij = (max(i,j)*(max(i,j)-3)/2) + i + j
            packed(ij, 1) = packed(ij, 1) + unpacked(i,j)
         enddo
      enddo
!!$omp end parallel do
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
!!$omp parallel do schedule(static) private(i,j)
      do i = 1, N
         do j = i, N
            ij = (max(i,j)*(max(i,j)-3)/2) + i + j
            packed(ij, 1) = packed(ij, 1) + unpacked(i,j) + unpacked(j,i)
         enddo
      enddo
!!$omp end parallel do
!
   end subroutine symmetrize_and_add_to_packed
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
!!$omp parallel do schedule(static) private(i,j)
      do i = 1, N
         do j = 1, N
            unpacked(i, j) = packed(index_packed(i,j), 1)
         enddo
      enddo
!!$omp end parallel do
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
!     Max batch size
!
      max_batch_length = available/(required/batch_dimension)
!
!     Number of full batches
!
      n_batch=batch_dimension/max_batch_length
!
!     Test for rest
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
!           Test for rest
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
