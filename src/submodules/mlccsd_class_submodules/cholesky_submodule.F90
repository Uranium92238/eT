submodule (mlccsd_class) cholesky
!
!!
!!    Cholesky submodule (MLCCSD)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!

!
   implicit none 
!
!
contains
!
!
   module subroutine read_transform_cholesky_for_CC2_amplitude_mlccsd(wf)
!
      implicit none
!
      class(mlccsd) :: wf
!
!
      integer(i15) :: unit_chol_ao = 0, unit_chol_mo_ai = 0, ioerror = 0 
      integer(i15) :: J, i, a, ai 
!
      real(dp), dimension(:,:), allocatable :: C_o, C_v
      real(dp), dimension(:,:), allocatable :: chol_ao_sq, chol_ao
      real(dp), dimension(:,:), allocatable :: L_ai_J
      real(dp), dimension(:,:), allocatable :: X
!
      call allocator(C_o, wf%n_ao, wf%n_total_active_o)
      call allocator(C_v, wf%n_ao, wf%n_total_active_v)
!
      do i = 1, wf%n_total_active_o
!
         C_o(:,i) = wf%mo_coef_cc2_ccs(:,i)
!
      enddo
!
      do a = 1, wf%n_total_active_v
!
         C_v(:,a) = wf%mo_coef_cc2_ccs(: , wf%n_o + a)
!
      enddo
!
!     Read AO-cholesky
!
      call generate_unit_identifier(unit_chol_ao)
      open(unit=unit_chol_ao, file='mlcc_cholesky', status='old', form='formatted')
      rewind(unit_chol_ao)
!
!     Read the number of Cholesky vectors (n_J) and 
!     the number of atomic orbitals (n_ao)
!
      read(unit_chol_ao,*) wf%n_ao, wf%n_J
!
      call generate_unit_identifier(unit_chol_mo_ai)
!
      open(unit=unit_chol_mo_ai, file='cholesky_ai_cc2_amplitudes', action='write', status='unknown', &
           access='direct', form='unformatted', recl=dp*(wf%n_J), iostat=ioerror)
      if (ioerror .ne. 0) then
            write(unit_output,*)'WARNING: error while creating cholesky_ai_cc2_amplitudes'
            stop
         endif
!
      call allocator(L_ai_J, wf%n_total_active_o*wf%n_total_active_v, wf%n_J)
!
      do J = 1, wf%n_J
!
         call allocator(chol_ao, wf%n_ao*(wf%n_ao+1)/2, 1)
         call allocator(chol_ao_sq, wf%n_ao, wf%n_ao)
         chol_ao = zero
         chol_ao_sq = zero
!
         read(unit_chol_ao,*) (chol_ao(i,1), i = 1, wf%n_ao*(wf%n_ao+1)/2)
!         
         call squareup(chol_ao, chol_ao_sq, wf%n_ao)
!
         call deallocator(chol_ao, wf%n_ao*(wf%n_ao+1)/2, 1)
!
         call allocator(X, wf%n_ao, wf%n_total_active_o)
!
         call dgemm('N','N',     &
                     wf%n_ao,    &
                     wf%n_total_active_o, &
                     wf%n_ao,    &
                     one,        &
                     chol_ao_sq, &
                     wf%n_ao,    &
                     C_o,        &
                     wf%n_ao,    &
                     zero,       &
                     X,          &
                     wf%n_ao)
!
         call deallocator(chol_ao_sq, wf%n_ao, wf%n_ao)
!
         call dgemm('T','N',     &
                     wf%n_total_active_v, &
                     wf%n_total_active_o, &
                     wf%n_ao,    &
                     one,        &
                     C_v,        &
                     wf%n_ao,    &
                     X,          &
                     wf%n_ao,    &
                     zero,       &
                     L_ai_J(1,J),&
                     wf%n_total_active_v)
!
         call deallocator(X, wf%n_ao, wf%n_total_active_o)
      enddo
      
!
      close(unit_chol_ao)
      do ai = 1, (wf%n_total_active_v)*(wf%n_total_active_o)
         write(unit_chol_mo_ai, rec=ai) (L_ai_J(ai,j), j = 1, wf%n_J)
      enddo
      call deallocator(L_ai_J, (wf%n_total_active_o)*(wf%n_total_active_v), wf%n_J)
      close(unit_chol_mo_ai)
!
      call deallocator(C_o, wf%n_ao, wf%n_total_active_o)
      call deallocator(C_v, wf%n_ao, wf%n_total_active_v)
!
   end subroutine read_transform_cholesky_for_CC2_amplitude_mlccsd
!
!
   subroutine read_cholesky_ai_for_cc2_amplitudes_mlccsd(wf,L_ai_J, a_first, a_last, i_first, i_last)
!!
!!    Read Cholesky IA 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!    Reads the MO Cholesky IA (occ-vir) vectors from file and
!!    places them in the incoming L_ia_J matrix
!!
!!
!!    Optional arguments: i_first, i_last, a_first, a_last can be used in order to restrict indices
!!
         implicit none
!
         class(mlccsd)            :: wf    
         real(dp), dimension(:,:) :: L_ai_J ! L_ia^J
         integer(i15)             :: i_first, a_first     ! First index (can differ from 1 when batching or for mlcc)
         integer(i15)             :: i_last, a_last      ! Last index (can differ from n_o/n_v when batching or for mlcc)
!
!        Local routine variables
!
         integer(i15) :: unit_chol_mo_ai = -1 ! Unit identifier for cholesky_ia file
         integer(i15) :: ioerror = 0
!
         integer(i15) :: i = 0, j = 0, a = 0, ai = 0, ai_full = 0
!
         integer(i15) :: active_space
         integer(i15) :: n_active_o = 0, n_active_v = 0
!
         n_active_o = i_last - i_first + 1 
         n_active_v = a_last - a_first + 1 
!
!        Prepare for reading: generate unit idientifier, open, and rewind file
!
         call generate_unit_identifier(unit_chol_mo_ai)
         open(unit=unit_chol_mo_ai, file='cholesky_ai_cc2_amplitudes', action='read', status='unknown', &
              access='direct', form='unformatted', recl=dp*(wf%n_J), iostat=ioerror)
!
         if (ioerror .ne. 0) then
            write(unit_output,*)'WARNING: error while reading cholesky_ai_cc2_amplitudes'
            stop
         endif
!
!        Read Cholesky vectors into the L_ia_J matrix
!
         do i = 1, n_active_o
            do a = 1, n_active_v
               ai = index_two(a, i, n_active_v)
               ai_full = index_two(a + a_first - 1, i + i_first - 1, wf%n_total_active_v)
!
               read(unit_chol_mo_ai, rec=ai_full) (L_ai_J(ai,j), j = 1, wf%n_J)
            enddo
!
         enddo
!
!        Close file
!
         close(unit_chol_mo_ai)
!   
   end subroutine read_cholesky_ai_for_cc2_amplitudes_mlccsd
!
!
  module subroutine get_cholesky_ai_for_cc2_amplitudes_mlccsd(wf, L_ai_J, a_first, a_last, i_first, i_last)
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
!!
!!      (1) determines memory requirement. 
!!
      implicit none 
!
      class(mlccsd)            :: wf
      real(dp), dimension(:,:) :: L_ai_J
      integer(i15)             :: i_first, a_first     ! First index (can differ from 1 when batching or for mlcc)
      integer(i15)             :: i_last, a_last      ! Last index (can differ from n_o/n_v when batching or for mlcc)
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
      integer(i15) :: n_active_o = 0, n_active_v = 0
!
      n_active_o = i_last - i_first + 1
      n_active_v = a_last - a_first + 1
!
!     Read L_ai^J from file 
!
      call wf%read_cholesky_ai_for_cc2_amplitudes(L_ai_J, a_first, a_last, i_first, i_last)
!                       
!     :: L_ab_J contributions ::
!
!
!     Allocate L_Ja_i
!
      call allocator(L_Ja_i, (wf%n_J)*n_active_v, n_active_o)
      L_Ja_i = zero
!
!     Set batching variables 
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
!     Calculate the number of batches 
!
      call num_batch(required, available, max_batch_length, n_batch, n_active_v)
!
      do a_batch = 1, n_batch
!
!        Get start, end, and length of batch
!
         call batch_limits(batch_start, batch_end, a_batch, max_batch_length, n_active_v)
         batch_length = batch_end - batch_start + 1
!
         call allocator(L_ab_J, batch_length*(wf%n_v), wf%n_J) ! L_ab^J
         L_ab_J = zero
!
!        Read Cholesky AB vectors, batching over a
! 
         call wf%read_cholesky_ab(L_ab_J, batch_start, batch_end, 1, wf%n_v)
!
         call allocator(L_ba_J, (wf%n_v)*batch_length, wf%n_J) ! L_ab^J = L_ba_J(ba,J)
!
         L_ba_J = zero
!
         do b = 1, wf%n_v
           do a = 1, batch_length
!
             ba = index_two(b, a, wf%n_v)
             ab = index_two(a, b, batch_length)
!
             do J = 1, wf%n_J
!
               L_ba_J(ba, J) = L_ab_J(ab, J)
!
             enddo
           enddo
         enddo
!
         call deallocator(L_ab_J, batch_length*(wf%n_v), wf%n_J)
!
         call allocator(L_Ja_b, batch_length*(wf%n_J), wf%n_v)
         L_Ja_b = zero
!
!        Reorder the Cholesky array L_ba_J
!
         do a = 1, batch_length
            do b = 1, wf%n_v
               do J = 1, wf%n_J
!
!                 Needed indices
!
                  ba = index_two(b, a, wf%n_v)
                  Ja = index_two(J, a, wf%n_J)
!
                  L_Ja_b(Ja, b) = L_ba_J(ba, J) ! L_ab^J
!
               enddo
            enddo
         enddo
         call deallocator(L_ba_J, (wf%n_v)*batch_length, wf%n_J)    
!
!        Calculate sum_b L_Ja_b*t_b_i = L_Ja_i 
!        
         L_off = index_two(1, batch_start, wf%n_J)
!
!
         call dgemm('N','N',                &
                     batch_length*(wf%n_J), &
                     n_active_o,      &
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
!
!        Deallocate  L_Ja_b
!
         call deallocator(L_Ja_b, batch_length*(wf%n_J), wf%n_v)
!
      enddo ! batching over a 
!
!     Add terms to T1-transformed Cholesky AI vector 
!
      do i = 1, n_active_o
         do a = 1, n_active_v
            do J = 1, wf%n_J
!
!              Needed indices
!
               Ja = index_two(J, a, wf%n_J)
               ai = index_two(a, i, n_active_v)
!
               L_ai_J(ai, J) = L_ai_J(ai, J) + L_Ja_i(Ja, i)
!
            enddo
         enddo
      enddo
!
!     Deallocate L_Ja_i
!
      call deallocator(L_Ja_i, (wf%n_J)*n_active_v, n_active_o)
!
!
!     :: L_ij_J contributions ::
!
!
!     Allocate L_a_iJ, L_ik_J, L_k_iJ
!
      call allocator(L_a_iJ, n_active_v, (wf%n_J)*n_active_o)
      call allocator(L_k_iJ, wf%n_o, n_active_o*(wf%n_J))
!
      call allocator(L_ik_J, (wf%n_o)*n_active_o, wf%n_J)
!  
!     Read Cholesky IJ vectors
!
      call wf%read_cholesky_ij(L_ik_J, 1, n_active_o, 1, wf%n_o) ! L_ik_J(ik,J) = L_ik^J 
!
!     Reorder IJ Cholesky vectors
!
      do i = 1, n_active_o
         do k = 1, wf%n_o
            do J = 1, wf%n_J
!
!              Needed indices
!
               ik = index_two(i, k, n_active_o)
               iJ = index_two(i, J, n_active_o)
!
               L_k_iJ(k, iJ) = L_ik_J(ik, J) ! L_k_iJ(k,iJ) = L_ik^J
!
            enddo
         enddo
      enddo
!
!     Calculate -sum_k t_a_k*L_k_iJ = L_a_iJ  ! Here we assume L_ik^J = L_ki^J 
!
      call dgemm('N','N',            &
                  n_active_v,          &
                  n_active_o*(wf%n_J), &
                  wf%n_o,            &
                  -one,              &
                  wf%t1am,           &
                  wf%n_v,            &
                  L_k_iJ,            &
                  wf%n_o,            &
                  zero,              &
                  L_a_iJ,            &
                  n_active_v)
!
!     Add terms to T1-transformation of L_ai_J
!
      do i = 1, n_active_o
         do a = 1, n_active_v
            do J = 1, wf%n_J
!
!              Needed indices
!
               ai = index_two(a, i, n_active_v)
               iJ = index_two(i, J, n_active_o)
!
               L_ai_J(ai, J) = L_ai_J(ai, J) + L_a_iJ(a, iJ)
!
            enddo
         enddo
      enddo
!
!     Deallocate L_a_iJ, L_ik_J, L_k_iJ
!
      call deallocator(L_a_iJ, n_active_v, (wf%n_J)*n_active_o)      
      call deallocator(L_k_iJ, wf%n_o, n_active_o*(wf%n_J))
!
      call deallocator(L_ik_J, (wf%n_o)*n_active_o, wf%n_J)
!
!
!     :: L_jb_J contributions ::    
!
!
      call allocator(L_kJ_b, (wf%n_o)*(wf%n_J), wf%n_v)
      call allocator(L_kb_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!     Read the Cholesky vector L_kb_J
!
      call wf%read_cholesky_ia(L_kb_J)
!
!     Reorder L_kb_J to L_kJ_b
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
!     Deallocate L_kb_J
!
      call deallocator(L_kb_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!     Allocate L_kJ_i 
!
      call allocator(L_kJ_i, (wf%n_o)*(wf%n_J), n_active_o)
!
!     Calculate sum_b L_kJ_b*t_b_i = L_kJ_i
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_J), &
                  n_active_o,  &
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
!     Deallocate L_kJ_b
!
      call deallocator(L_kJ_b, (wf%n_o)*(wf%n_J), wf%n_v)
!
!     Allocate L_k_iJ
!  
      call allocator(L_k_iJ, (wf%n_o), n_active_o*(wf%n_J))
!
!     Reorder L_kJ_i to L_k_iJ    
!
      do i = 1, n_active_o
         do k = 1, wf%n_o
            do J = 1, wf%n_J
!
               kJ = index_two(k, J, wf%n_o)
               iJ = index_two(i, J, n_active_o)
!
               L_k_iJ(k, iJ) = L_kJ_i(kJ, i)
!
            enddo
         enddo
      enddo
!
!     Deallocate L_kJ_i
!
      call deallocator(L_kJ_i, (wf%n_o)*(wf%n_J), n_active_o)
!
!     Allocate L_a_iJ
!
      call allocator(L_a_iJ, n_active_v, n_active_o*(wf%n_J))
!      
!     Calculate sum_k t_a_k*L_k_iJ = L_a_iJ
!
      call dgemm('N','N',            &
                  n_active_v,          &
                  n_active_o*(wf%n_J), &
                  wf%n_o,            &
                  -one,              &
                  wf%t1am,           &
                  wf%n_v,            &
                  L_k_iJ,            &
                  wf%n_o,            &
                  zero,              &
                  L_a_iJ,            &
                  wf%n_total_active_v)
!
!     Add contribution to L ai_J
!
      do a = 1, n_active_v
         do i = 1, n_active_o
            do J = 1, wf%n_J
!
               iJ = index_two(i, J, n_active_o)
               ai = index_two(a, i, n_active_v)
!
               L_ai_J(ai, J) = L_ai_J(ai, J) + L_a_iJ(a, iJ)
!
            enddo
         enddo
      enddo
!
!     Deallocations
!
      call deallocator(L_a_iJ, n_active_v, n_active_o*(wf%n_J))
      call deallocator(L_k_iJ, wf%n_o, n_active_o*(wf%n_J))  
!
   end subroutine get_cholesky_ai_for_cc2_amplitudes_mlccsd
!
!
end submodule cholesky