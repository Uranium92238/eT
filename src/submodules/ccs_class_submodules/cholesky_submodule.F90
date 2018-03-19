submodule(ccs_class) cholesky
!
!!
!!    Cholesky sub(CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!    Contains the following family of procedures of the CCS class:
!!
!!    get_cholesky_ij(L_ij_J):         returns the T1-transformed Cholesky vector L_ij^J 
!!    get_cholesky_ia(L_ia_J):         returns the T1-transformed Cholesky vector L_ia^J 
!!    get_cholesky_ai(L_ai_J):         returns the T1-transformed Cholesky vector L_ai^J
!!    get_cholesky_ab(L_ab_J, ...):    returns the T1-transformed Cholesky vector L_ab^J,
!!                                     but has options (...) for batching over the two 
!!                                     indices a and b
!!
!
   implicit none 
!
!
contains
!
!

   module subroutine get_cholesky_ij_ccs(wf, L_ij_J, i_first, i_last, j_first, j_last)
!!
!!    Get Cholesky IJ
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!    Reads and T1-transforms the IA Cholesky vectors:
!!     
!!       L_ij_J_T1 = L_ij_J + sum_a t_aj * L_ia_J
!!
!!    Memory required in routine:
!!
!!       2*n_J*(i_length)*n_v     -> for reading L_ia_J contribution and reordering
!!       i_length = i_last - i_first + 1
!!
!!    Optional arguments: i_first, i_last, j_first, j_last can be used in order to restrict indices
!!        
      implicit none 
!
      class(ccs)               :: wf
      integer(i15), optional   :: i_first, j_first     ! First index (can differ from 1 when batching or for mlcc)
      integer(i15), optional   :: i_last, j_last      ! Last index (can differ from n_o when batching or for mlcc)
      real(dp), dimension(:,:) :: L_ij_J
!
!     Local routine variables
!
      real(dp), dimension(:,:), allocatable :: L_ia_J ! L_ia^J
      real(dp), dimension(:,:), allocatable :: L_iJ_a ! L_ia^J reordered
      real(dp), dimension(:,:), allocatable :: L_iJ_k ! L_ik^J reordered
!
      integer(i15) :: i = 0, J = 0, a = 0, ij = 0, ia = 0, ik = 0, k = 0
      integer(i15) :: i_length, j_length
!
      if (present(i_first) .and. present(i_last) .and. present(j_first) .and. present(j_last)) then
!
         i_length = i_last - i_first + 1
         j_length = j_last - j_first + 1
!
!        Allocate
!
         call wf%mem%alloc(L_ia_J, i_length*(wf%n_v), wf%n_J)
         call wf%mem%alloc(L_iJ_a, i_length*(wf%n_J), wf%n_v)
!
!        Read the untransformed Cholesky vectors 
!
         call wf%read_cholesky_ij(L_ij_J, i_first, i_last, j_first, j_last)
         call wf%read_cholesky_ia(L_ia_J, i_first, i_last, 1, wf%n_v)
!
!        Reorder L_ia_J to L_iJ_a 
!
         do i = 1, i_length
            do J = 1, wf%n_J
               do a = 1, wf%n_v
!                 
                  iJ = index_two(i, J, i_length)
                  ia = index_two(i, a, i_length)
!
                  L_iJ_a(iJ, a) = L_ia_J(ia, J)
!
               enddo
            enddo
         enddo
!
!        Deallocate L_ia_J
!
         call wf%mem%dealloc(L_ia_J, i_length*(wf%n_v), wf%n_J)
!
!        Allocate L_iJ_k (L_ij^J reordered as L_iJ_j)
!
         call wf%mem%alloc(L_iJ_k, i_length*(wf%n_J), j_length)
!
!        T1-transformation
!
         call dgemm('N','N',            &
                     i_length*(wf%n_J), &
                     j_length,          &
                     wf%n_v,            &
                     one,               &
                     L_iJ_a,            &
                     i_length*(wf%n_J), &
                     wf%t1am(1, j_first),&
                     wf%n_v,            &
                     zero,              &
                     L_iJ_k,            &
                     i_length*(wf%n_J))
!
!        Place terms from L_iJ_k into L_ij_J
!
         do i = 1, i_length
            do k = 1, j_length
               do J = 1, wf%n_J
!                 
!                 Needed indices
!
                  iJ = index_two(i, J, i_length)
                  ik = index_two(i, k, i_length)
!
                  L_ij_J(ik, J) = L_ij_J(ik, J) + L_iJ_k(iJ, k)
!
               enddo
            enddo
         enddo
!
!        Deallocate L_iJ_k and L_iJ_a
!
         call wf%mem%dealloc(L_iJ_k, i_length*(wf%n_J), j_length)
         call wf%mem%dealloc(L_iJ_a, i_length*(wf%n_J), wf%n_v)
!
    elseif (.not. (present(i_first) .and. present(i_last) .and. present(j_first) .and. present(j_last))) then
!    Allocate
!
     call wf%mem%alloc(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
     call wf%mem%alloc(L_iJ_a, (wf%n_o)*(wf%n_J), wf%n_v)
!
!    Read the untransformed Cholesky vectors 
!
     call wf%read_cholesky_ij(L_ij_J)
     call wf%read_cholesky_ia(L_ia_J)
!
!    Reorder L_ia_J to L_iJ_a 
!
     do i = 1, wf%n_o
        do J = 1, wf%n_J
           do a = 1, wf%n_v
!             
              iJ = index_two(i, J, wf%n_o)
              ia = index_two(i, a, wf%n_o)
!
              L_iJ_a(iJ, a) = L_ia_J(ia, J)
!
           enddo
        enddo
     enddo
!
!    Deallocate L_ia_J
!
     call wf%mem%dealloc(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!    Allocate L_iJ_k
!
     call wf%mem%alloc(L_iJ_k, (wf%n_o)*(wf%n_J), wf%n_o)
!
!    T1-transformation
!
     call dgemm('N','N',            &
                 (wf%n_o)*(wf%n_J), &
                 wf%n_o,            &
                 wf%n_v,            &
                 one,               &
                 L_iJ_a,            &
                 (wf%n_o)*(wf%n_J), &
                 wf%t1am,           &
                 wf%n_v,            &
                 zero,              &
                 L_iJ_k,            &
                 (wf%n_o)*(wf%n_J))
!
!    Place terms from L_iJ_k into L_ij_J
!
     do i = 1, wf%n_o
        do k = 1, wf%n_o
           do J = 1, wf%n_J
!             
!             Needed indices
!
              iJ = index_two(i, J, wf%n_o)
              ik = index_two(i, k, wf%n_o)
!
              L_ij_J(ik, J) = L_ij_J(ik, J) + L_iJ_k(iJ, k)
!
           enddo
        enddo
     enddo
!
!    Deallocate L_iJ_k and L_iJ_a
!
     call wf%mem%dealloc(L_iJ_k, (wf%n_o)*(wf%n_J), wf%n_o)
     call wf%mem%dealloc(L_iJ_a, (wf%n_o)*(wf%n_J), wf%n_v)
    else
      write(unit_output, *) 'Error: in call to get_cholesky_ij'
      stop
   endif   
! 
   end subroutine get_cholesky_ij_ccs
!
!
   module subroutine get_cholesky_ia_ccs(wf, L_ia_J, i_first, i_last, a_first, a_last)
!!
!!    Get Cholesky IA
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!    Reads and T1-transforms IA Cholesky vectors
!!
!!       L_ia_J_T1 = L_ia_J (only reading necessary)
!!
!!    Memory required in routine:
!!
!!       No additional memory
!!
!!    Optional arguments: i_first, i_last, a_first, a_last can be used in order to restrict indices
!!
      implicit none 
!
      class(ccs)               :: wf
      integer(i15), optional   :: i_first, a_first   ! First index (can differ from 1 when batching or for mlcc)
      integer(i15), optional   :: i_last, a_last    ! Last index (can differ from n_v/n_o when batching or for mlcc)
      real(dp), dimension(:,:) :: L_ia_J
!
      call wf%read_cholesky_ia(L_ia_J,i_first, i_last, a_first, a_last)
!
   end subroutine get_cholesky_ia_ccs
!
!
   module subroutine get_cholesky_ai_ccs(wf, L_ai_J, a_first, a_last, i_first, i_last)

!!
!!    Get Cholesky AI
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!    Read and T1-transform Cholesky AI vectors:
!!    
!!       L_ai_J_T1 = L_ia_J - sum_j  t_aj*L_ji_J 
!!                          + sum_b  t_bi*L_ab_J
!!                          - sum_bj t_aj*t_bi*L_jb_J
!!
!!    Allocations in routine:
!!
!!      (1) n_J*(i_length)*(a_length) + 2*n_J*(a_length)*batch_length  ->  for L_ab_J contribution (batches of b)
!!      (2) n_J*(i_length)*n_v + 2*n_J*n_o*(i_length)                  ->  for L_ij_J contribution
!!      (3) 2*n_J*n_o*n_v                                              ->  for L_jb_J contribution
!!
!!      i_length = i_last - i_first + 1          
!!      a_length = a_last - a_first + 1          
!!
!!      (1) determines memory requirement. 
!!
!!    Optional arguments: i_first, i_last, a_first, a_last can be used in order to restrict indices
!!
      implicit none 
!
      class(ccs)               :: wf
      integer(i15), optional   :: i_first, a_first     ! First index (can differ from 1 when batching or for mlcc)
      integer(i15), optional   :: i_last, a_last      ! Last index (can differ from n_o/n_v when batching or for mlcc)
      real(dp), dimension(:,:) :: L_ai_J
!
!     Local routine variables
!
      logical :: reorder ! Reorder or not, when reading Cholesky AB 
!
!     Batch variables
!
      integer(i15) :: required = 0, available = 0, max_batch_length = 0, n_batch = 0, L_off = 0
      integer(i15) :: a_batch = 0, batch_start = 0, batch_end = 0, batch_length = 0
!
!     Indices
!
      integer(i15) :: a = 0, b = 0, J = 0, i = 0, ai = 0, Ja = 0
      integer(i15) :: ba = 0, k = 0, ik = 0, iJ = 0, kb = 0, kJ = 0, ab = 0
      integer(i15) :: a_length, i_length
!
!     Cholesky vectors (in many different orderings)
!
      real(dp), dimension(:,:), allocatable :: L_ba_J
      real(dp), dimension(:,:), allocatable :: L_ab_J
      real(dp), dimension(:,:), allocatable :: L_Ja_b
      real(dp), dimension(:,:), allocatable :: L_Ja_i
      real(dp), dimension(:,:), allocatable :: L_ik_J
      real(dp), dimension(:,:), allocatable :: L_k_iJ
      real(dp), dimension(:,:), allocatable :: L_a_iJ
      real(dp), dimension(:,:), allocatable :: L_kJ_b
      real(dp), dimension(:,:), allocatable :: L_kJ_i
      real(dp), dimension(:,:), allocatable :: L_kb_J
!
      if (present(i_first) .and. present(i_last) .and. present(a_first) .and. present(a_last)) then
         a_length = a_last - a_first + 1
         i_length = i_last - i_first + 1
!
!        Read L_ai^J from file 
!
         call wf%read_cholesky_ai(L_ai_J, a_first, a_last, i_first, i_last)
!                          
!        :: L_ab_J contributions ::
!
!        Allocate L_Ja_i
!
         call wf%mem%alloc(L_Ja_i, (wf%n_J)*a_length, i_length)
         L_Ja_i = zero
!
!        Set batching variables 
!
         required = 2*(wf%n_v)*(a_last - a_first + 1)*(wf%n_J)*4
         available = get_available()
         max_batch_length = 0
!
         n_batch = 0
         a_batch = 0
!
         batch_length = 0
         batch_start = 0
         batch_end = 0
!
!        Calculate the number of batches 
!
         call num_batch(required, available, max_batch_length, n_batch, a_length)
!
         do a_batch = 1, n_batch
!
!           Get start, end, and length of batch
!
            call batch_limits(batch_start, batch_end, a_batch, max_batch_length, a_length)
            batch_start  = batch_start + (a_first - 1)
            batch_end    = batch_end + (a_first - 1)
!
            if (batch_end .gt. a_last) batch_end = a_last
            batch_length = batch_end - batch_start + 1
!
!           Read L_ab_J 
! 
            call wf%mem%alloc(L_ab_J, (wf%n_v)*batch_length, wf%n_J) ! L_ab^J 
            L_ab_J = zero
!
             call wf%read_cholesky_ab(L_ab_J, batch_start, batch_end, 1, wf%n_v)           
!         
            call wf%mem%alloc(L_Ja_b, batch_length*(wf%n_J), wf%n_v)
            L_Ja_b = zero           
!
!           Reorder the Cholesky array L_ab_J
!   
            call sort_123_to_312(L_ab_J, L_Ja_b, batch_length, wf%n_v, wf%n_J)     
!
            call wf%mem%dealloc(L_ab_J, (wf%n_v)*batch_length, wf%n_J)   
!
!           Calculate sum_b L_Ja_b*t_b_i = L_Ja_i 
!           
            L_off = index_two(1, batch_start, wf%n_J)
!
            call dgemm('N','N',                &
                        batch_length*(wf%n_J), &
                        i_length,              &
                        wf%n_v,                &
                        one,                   &
                        L_Ja_b,                &
                        batch_length*(wf%n_J), &
                        wf%t1am(1,i_first),    &
                        wf%n_v,                &
                        one,                   &
                        L_Ja_i(L_off, 1),      &
                        a_length*(wf%n_J))
!
            call wf%mem%dealloc(L_Ja_b, batch_length*(wf%n_J), wf%n_v)
!
         enddo ! batching over a 
!
!        Add terms to T1-transformed Cholesky AI vector 
!
         do i = 1, i_length
            do a = 1, a_length
               do J = 1, wf%n_J
!
!                 Needed indices
!
                  Ja = index_two(J, a, wf%n_J)
                  ai = index_two(a, i, a_length)
!
                  L_ai_J(ai, J) = L_ai_J(ai, J) + L_Ja_i(Ja, i)
!
               enddo
            enddo
         enddo

!
!        Deallocate L_Ja_i
!
         call wf%mem%dealloc(L_Ja_i, (wf%n_J)*a_length, i_length)
!
!
!        :: L_ij_J contributions ::
!
!
!        Allocate L_a_iJ, L_ik_J, L_k_iJ
!
         call wf%mem%alloc(L_a_iJ, a_length, (wf%n_J)*i_length)
         call wf%mem%alloc(L_k_iJ, wf%n_o,  i_length*(wf%n_J))
!
         call wf%mem%alloc(L_ik_J, (wf%n_o)*i_length, wf%n_J)
!  
!        Read Cholesky IJ vectors
!
         call wf%read_cholesky_ij(L_ik_J, i_first, i_last, 1, wf%n_o) ! L_ik_J(ik,J) = L_ik^J
!
!        Reorder IJ Cholesky vectors
!
         do i = 1, i_length
            do k = 1, wf%n_o
               do J = 1, wf%n_J
!
!                 Needed indices
!
                  ik = index_two(i, k, i_length)
                  iJ = index_two(i, J, i_length)
!
                  L_k_iJ(k, iJ) = L_ik_J(ik, J) ! L_k_iJ(k,iJ) = L_ik^J
!
               enddo
            enddo
         enddo
!
!        Calculate -sum_k t_a_k*L_k_iJ = L_a_iJ  ! Here we assume L_ik^J = L_ki^J 
!
         call dgemm('N','N',            &
                     a_length,          &
                     i_length*(wf%n_J), &
                     wf%n_o,            &
                     -one,              &
                     wf%t1am(a_first,1),&
                     wf%n_v,            &
                     L_k_iJ,            &
                     wf%n_o,            &
                     zero,              &
                     L_a_iJ,            &
                     a_length)
!
!        Add terms to T1-transformation of L_ai_J
!
         do i = 1, i_length
            do a = 1, a_length
               do J = 1, wf%n_J
!
!                 Needed indices
!
                  ai = index_two(a, i, a_length)
                  iJ = index_two(i, J, i_length)
!
                  L_ai_J(ai, J) = L_ai_J(ai, J) + L_a_iJ(a, iJ)
!
               enddo
            enddo
         enddo
!
!        Deallocate L_a_iJ, L_ik_J, L_k_iJ
!
         call wf%mem%dealloc(L_a_iJ, a_length, (wf%n_J)*i_length)      
         call wf%mem%dealloc(L_k_iJ, wf%n_o, i_length*(wf%n_J))
!
         call wf%mem%dealloc(L_ik_J, (wf%n_o)*i_length, wf%n_J)
!
!
!        :: L_jb_J contributions ::    
!
!
         call wf%mem%alloc(L_kJ_b, (wf%n_o)*(wf%n_J), wf%n_v)
         call wf%mem%alloc(L_kb_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!        Read the Cholesky vector L_kb_J
!
         call wf%read_cholesky_ia(L_kb_J)
!
!        Reorder L_kb_J to L_kJ_b
!
         do k = 1, wf%n_o
            do b = 1, wf%n_v
               do J = 1, wf%n_J
!
                  kb = index_two(k, b, wf%n_o)
                  kJ = index_two(k, J, wf%n_o)
!
                  L_kJ_b(kJ, b) = L_kb_J(kb, J)
!
               enddo
            enddo
         enddo
!
!        Deallocate L_kb_J
!
         call wf%mem%dealloc(L_kb_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!        Allocate L_kJ_i 
!
         call wf%mem%alloc(L_kJ_i, (wf%n_o)*(wf%n_J), i_length)
!
!        Calculate sum_b L_kJ_b*t_b_i = L_kJ_i
!
         call dgemm('N','N',            &
                     (wf%n_o)*(wf%n_J), &
                     i_length,          &
                     wf%n_v,            &
                     one,               &
                     L_kJ_b,            &
                     (wf%n_o)*(wf%n_J), &
                     wf%t1am(1,i_first),&
                     wf%n_v,            &
                     zero,              &
                     L_kJ_i,            &
                     (wf%n_o)*(wf%n_J))
!
!        Deallocate L_kJ_b
!
         call wf%mem%dealloc(L_kJ_b, (wf%n_o)*(wf%n_J), wf%n_v)
!
!        Allocate L_k_iJ
!  
         call wf%mem%alloc(L_k_iJ, (wf%n_o), i_length*(wf%n_J))
!
!        Reorder L_kJ_i to L_k_iJ    
!
         do i = 1, i_length
            do k = 1, wf%n_o
               do J = 1, wf%n_J
!
                  kJ = index_two(k, J, wf%n_o)
                  iJ = index_two(i, J, i_length)
!
                  L_k_iJ(k, iJ) = L_kJ_i(kJ, i)
!
               enddo
            enddo
         enddo
!
!        Deallocate L_kJ_i
!
         call wf%mem%dealloc(L_kJ_i, (wf%n_o)*(wf%n_J), i_length)
!
!        Allocate L_a_iJ
!
         call wf%mem%alloc(L_a_iJ, a_length, i_length*(wf%n_J))
!         
!        Calculate sum_k t_a_k*L_k_iJ = L_a_iJ
!
         call dgemm('N','N',            &
                     a_length,          &
                     i_length*(wf%n_J), &
                     wf%n_o,            &
                     -one,              &
                     wf%t1am(a_first,1),&
                     wf%n_v,            &
                     L_k_iJ,            &
                     wf%n_o,            &
                     zero,              &
                     L_a_iJ,            &
                     a_length)
!
!        Add contribution to L ai_J
!
         do a = 1, a_length
            do i = 1, i_length
               do J = 1, wf%n_J
!
                  iJ = index_two(i, J, i_length)
                  ai = index_two(a, i, a_length)
!
                  L_ai_J(ai, J) = L_ai_J(ai, J) + L_a_iJ(a, iJ)
!
               enddo
            enddo
         enddo
!
!        Deallocations
!
         call wf%mem%dealloc(L_a_iJ, a_length, i_length*(wf%n_J))
         call wf%mem%dealloc(L_k_iJ, wf%n_o, i_length*(wf%n_J))
!
      elseif (.not.(present(i_first) .and. present(i_last) .and. present(a_first) .and. present(a_last))) then
!
!        Read L_ai^J from file 
!
         call wf%read_cholesky_ai(L_ai_J)
!                          
!        :: L_ab_J contributions ::
!
!        Allocate L_Ja_i
!
         call wf%mem%alloc(L_Ja_i, (wf%n_J)*(wf%n_v), wf%n_o)
         L_Ja_i = zero
!
!        Set batching variables 
!
         required = 2*(wf%n_v)**2*(wf%n_J)*4
         available = get_available()
         max_batch_length = 0
!
         n_batch = 0
         a_batch = 0
!
         batch_length = 0
         batch_start = 0
         batch_end = 0
!
!        Calculate the number of batches 
!
         call num_batch(required, available, max_batch_length, n_batch, wf%n_v)
!
         do a_batch = 1, n_batch
!
!           Get start, end, and length of batch
!
            call batch_limits(batch_start, batch_end, a_batch, max_batch_length, wf%n_v)
            batch_length = batch_end - batch_start + 1
!
            call wf%mem%alloc(L_ab_J, batch_length*(wf%n_v), wf%n_J) ! L_ab^J
            L_ab_J = zero
!
!           Read Cholesky AB vectors, batching over a
! 
            call wf%read_cholesky_ab(L_ab_J, batch_start, batch_end, 1, wf%n_v)
!
            call wf%mem%alloc(L_Ja_b, batch_length*(wf%n_J), wf%n_v)
            L_Ja_b = zero
!
!           Reorder the Cholesky array L_ab_J
!
            call sort_123_to_312(L_ab_J, L_Ja_b, batch_length, wf%n_v, wf%n_J)     
!
            call wf%mem%dealloc(L_ab_J, (wf%n_v)*batch_length, wf%n_J)    
!
!           Calculate sum_b L_Ja_b*t_b_i = L_Ja_i 
!           
            L_off = index_two(1, batch_start, wf%n_J)
!
            call dgemm('N','N',                &
                        batch_length*(wf%n_J), &
                        wf%n_o,                &
                        wf%n_v,                &
                        one,                   &
                        L_Ja_b,                &
                        batch_length*(wf%n_J), &
                        wf%t1am,               &
                        wf%n_v,                &
                        one,                   &
                        L_Ja_i(L_off, 1),      &
                        (wf%n_v)*(wf%n_J))
!
!           Deallocate  L_Ja_b
!
            call wf%mem%dealloc(L_Ja_b, batch_length*(wf%n_J), wf%n_v)
!
         enddo ! batching over a 
!
!        Add terms to T1-transformed Cholesky AI vector 
!
         do i = 1, wf%n_o
            do a = 1, wf%n_v
               do J = 1, wf%n_J
!
!                 Needed indices
!
                  Ja = index_two(J, a, wf%n_J)
                  ai = index_two(a, i, wf%n_v)
!
                  L_ai_J(ai, J) = L_ai_J(ai, J) + L_Ja_i(Ja, i)
!
               enddo
            enddo
         enddo
!
!        Deallocate L_Ja_i
!
         call wf%mem%dealloc(L_Ja_i, (wf%n_J)*(wf%n_v), wf%n_o)
!
!
!        :: L_ij_J contributions ::
!
!
!        Allocate L_a_iJ, L_ik_J, L_k_iJ
!
         call wf%mem%alloc(L_a_iJ, wf%n_v, (wf%n_J)*(wf%n_o))
         call wf%mem%alloc(L_k_iJ, wf%n_o, (wf%n_o)*(wf%n_J))
!
         call wf%mem%alloc(L_ik_J, (wf%n_o)**2, wf%n_J)
!  
!        Read Cholesky IJ vectors
!
         call wf%read_cholesky_ij(L_ik_J) ! L_ik_J(ik,J) = L_ik^J 
!
!        Reorder IJ Cholesky vectors
!
         do i = 1, wf%n_o
            do k = 1, wf%n_o
               do J = 1, wf%n_J
!
!                 Needed indices
!
                  ik = index_two(i, k, wf%n_o)
                  iJ = index_two(i, J, wf%n_o)
!
                  L_k_iJ(k, iJ) = L_ik_J(ik, J) ! L_k_iJ(k,iJ) = L_ik^J
!
               enddo
            enddo
         enddo
!
!        Calculate -sum_k t_a_k*L_k_iJ = L_a_iJ  ! Here we assume L_ik^J = L_ki^J 
!
         call dgemm('N','N',            &
                     wf%n_v,            &
                     (wf%n_o)*(wf%n_J), &
                     wf%n_o,            &
                     -one,              &
                     wf%t1am,           &
                     wf%n_v,            &
                     L_k_iJ,            &
                     wf%n_o,            &
                     zero,              &
                     L_a_iJ,            &
                     wf%n_v)
!
!        Add terms to T1-transformation of L_ai_J
!
         do i = 1, wf%n_o
            do a = 1, wf%n_v
               do J = 1, wf%n_J
!
!                 Needed indices
!
                  ai = index_two(a, i, wf%n_v)
                  iJ = index_two(i, J, wf%n_o)
!
                  L_ai_J(ai, J) = L_ai_J(ai, J) + L_a_iJ(a, iJ)
!
               enddo
            enddo
         enddo
!
!        Deallocate L_a_iJ, L_ik_J, L_k_iJ
!
         call wf%mem%dealloc(L_a_iJ, wf%n_v, (wf%n_J)*(wf%n_o))      
         call wf%mem%dealloc(L_k_iJ, wf%n_o, (wf%n_o)*(wf%n_J))
!
         call wf%mem%dealloc(L_ik_J, (wf%n_o)**2, wf%n_J)
!
!
!        :: L_jb_J contributions ::    
!
!
         call wf%mem%alloc(L_kJ_b, (wf%n_o)*(wf%n_J), wf%n_v)
         call wf%mem%alloc(L_kb_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!        Read the Cholesky vector L_kb_J
!
         call wf%read_cholesky_ia(L_kb_J)
!
!        Reorder L_kb_J to L_kJ_b
!
         do k = 1, wf%n_o
            do b = 1, wf%n_v
               do J = 1, wf%n_J
!
                  kb = index_two(k, b, wf%n_o)
                  kJ = index_two(k, J, wf%n_o)
!
                  L_kJ_b(kJ, b) = L_kb_J(kb, J)
!
               enddo
            enddo
         enddo
!
!        Deallocate L_kb_J
!
         call wf%mem%dealloc(L_kb_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!        Allocate L_kJ_i 
!
         call wf%mem%alloc(L_kJ_i, (wf%n_o)*(wf%n_J), wf%n_o)
!
!        Calculate sum_b L_kJ_b*t_b_i = L_kJ_i
!
         call dgemm('N','N',            &
                     (wf%n_o)*(wf%n_J), &
                     wf%n_o,            &
                     wf%n_v,            &
                     one,               &
                     L_kJ_b,            &
                     (wf%n_o)*(wf%n_J), &
                     wf%t1am,           &
                     wf%n_v,            &
                     zero,              &
                     L_kJ_i,            &
                     (wf%n_o)*(wf%n_J))
!
!        Deallocate L_kJ_b
!
         call wf%mem%dealloc(L_kJ_b, (wf%n_o)*(wf%n_J), wf%n_v)
!
!        Allocate L_k_iJ
!  
         call wf%mem%alloc(L_k_iJ, (wf%n_o), (wf%n_o)*(wf%n_J))
!
!        Reorder L_kJ_i to L_k_iJ    
!
         do i = 1, wf%n_o
            do k = 1, wf%n_o
               do J = 1, wf%n_J
!
                  kJ = index_two(k, J, wf%n_o)
                  iJ = index_two(i, J, wf%n_o)
!
                  L_k_iJ(k, iJ) = L_kJ_i(kJ, i)
!
               enddo
            enddo
         enddo
!
!        Deallocate L_kJ_i
!
         call wf%mem%dealloc(L_kJ_i, (wf%n_o)*(wf%n_J), wf%n_o)
!
!        Allocate L_a_iJ
!
         call wf%mem%alloc(L_a_iJ, wf%n_v, (wf%n_o)*(wf%n_J))
!         
!        Calculate sum_k t_a_k*L_k_iJ = L_a_iJ
!
         call dgemm('N','N',            &
                     wf%n_v,            &
                     (wf%n_o)*(wf%n_J), &
                     wf%n_o,            &
                     -one,              &
                     wf%t1am,           &
                     wf%n_v,            &
                     L_k_iJ,            &
                     wf%n_o,            &
                     zero,              &
                     L_a_iJ,            &
                     wf%n_v)
!
!        Add contribution to L ai_J
!
         do a = 1, wf%n_v
            do i = 1, wf%n_o
               do J = 1, wf%n_J
!
                  iJ = index_two(i, J, wf%n_o)
                  ai = index_two(a, i, wf%n_v)
!
                  L_ai_J(ai, J) = L_ai_J(ai, J) + L_a_iJ(a, iJ)
!
               enddo
            enddo
         enddo
!
!        Deallocations
!
         call wf%mem%dealloc(L_a_iJ, wf%n_v, (wf%n_o)*(wf%n_J))
         call wf%mem%dealloc(L_k_iJ, wf%n_o, (wf%n_o)*(wf%n_J))
!
      else
         write(unit_output, *) 'WARNING: Error in call to read_cholesky_ia'
         stop
      endif    
!
   end subroutine get_cholesky_ai_ccs
!
!
   module subroutine get_cholesky_ab_ccs(wf, L_ab_J, a_first, a_last, b_first, b_last)
!!
!!    Get Cholesky AB
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!    Reads and T1-transforms the IA Cholesky vectors:
!!
!!       L_ab_J_T1 = L_ab_J - sum_i t_ai*L_ib_J
!!
!!
!!    Required memory: 
!!
!!       n_J*b_length*a_length       ->   For reordering of L_ab_J / L_ba_J
!!       2*b_length*n_o*n_J          ->   For L_ib_J contribution
!!
!!      a_length = a_last - a_first + 1          
!!      b_length = b_last - b_first + 1 
!!
      implicit none
!
      class(ccs)               :: wf
      integer(i15), intent(in) :: a_first, b_first   ! First index (can differ from 1 when batching)
      integer(i15), intent(in) :: a_last, b_last    ! Last index  (can differ from n_v when batching)
      real(dp), dimension(((b_last - b_first + 1)*(a_last - a_first + 1)), wf%n_J) :: L_ab_J ! L_ab^J
!
!     Local routine variables
!
      integer(i15) :: memory_lef = 0
!
      integer :: unit_chol_ab = -1 ! Unit identifier for cholesky_ab file
!
      integer :: a = 0, b = 0, J = 0, i = 0, ia = 0, aJ = 0, ib = 0, Jb = 0, ab = 0, ba = 0
!
      real(dp), dimension(:,:), allocatable :: L_ib_J
      real(dp), dimension(:,:), allocatable :: L_Jb_i
      real(dp), dimension(:,:), allocatable :: L_Jb_a
      real(dp), dimension(:,:), allocatable :: L_a_Jb
      real(dp), dimension(:,:), allocatable :: L_i_Jb
!
      integer(i15) :: b_length, a_length
!
      a_length = a_last - a_first + 1
      b_length = b_last - b_first + 1
!
!     Allocate L_ib_J
!     
      call wf%mem%alloc(L_ib_J, (wf%n_o)*b_length, wf%n_J)
!
!     Read L_ia_J
!
!     Note: using L_ia_J instead of L_ai_J, here, to avoid two reorderings.
!           This is possible because of the symmetry L_ai_J(ai,J) == L_ia_J(ia,J).
!  
      call wf%read_cholesky_ia(L_ib_J, 1, wf%n_o, b_first, b_last)
!
!     Read L_ab_J for batch of b
!
      call wf%read_cholesky_ab(L_ab_J, a_first, a_last, b_first, b_last)
!
!     Allocate L_Jb,i for batch of b
!
      call wf%mem%alloc(L_Jb_i, (wf%n_J)*b_length, wf%n_o)
!
!     Reorder L_ib_J to L_Jb_i
!
      do i = 1, wf%n_o
         do b = 1, b_length
            do J = 1, wf%n_J
!
               ib = index_two(i, b, wf%n_o) 
               Jb = index_two(J, b, wf%n_J)
!
               L_Jb_i(Jb, i) = L_ib_J(ib, J)
!
            enddo
         enddo
      enddo
!
!     Dellocate L_ib_J
!  
      call wf%mem%dealloc(L_ib_J, (wf%n_o)*b_length, wf%n_J)
!
!     Allocate L_Jb_a for batch of b
!  
      call wf%mem%alloc(L_Jb_a, (wf%n_J)*b_length, a_length)
!
!     T1-transformation
!
      call dgemm('N','T',                &
                  (wf%n_J)*b_length,     & 
                  a_length,              &
                  wf%n_o,                &
                  -one,                  &
                  L_Jb_i,                &
                  (wf%n_J)*b_length,     &
                  wf%t1am(a_first, 1),   &
                  wf%n_v,                &
                  zero,                  &
                  L_Jb_a,                &
                  b_length*(wf%n_J))
!
!     Add terms of L_Jb_a to L_ab_J
!
      call add_321_to_123(one, L_Jb_a, L_ab_J, a_length, b_length, wf%n_J)
!
!     Dellocate L_Jb,i and L_Jb_a for batch of b
!
      call wf%mem%dealloc(L_Jb_a, (wf%n_J)*b_length, a_length)
      call wf%mem%dealloc(L_Jb_i, (wf%n_J)*b_length, wf%n_o)
!
   end subroutine get_cholesky_ab_ccs
!
end submodule cholesky
