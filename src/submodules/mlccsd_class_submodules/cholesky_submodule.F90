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
   subroutine read_transform_cholesky_for_CC2_amplitude_mlccsd(wf)
!
      implicit none
!
      class(mlccsd) :: wf
!
!
      integer(i15) :: unit_chol_ao = 0, ioerror = 0 
      integer(i15) :: unit_chol_mo_ai = -1
      integer(i15) :: unit_chol_mo_ia = -1
      integer(i15) :: unit_chol_mo_ij = -1
      integer(i15) :: unit_chol_mo_ab = -1
      integer(i15) :: unit_chol_mo_ab_tmp = -1
      integer(i15) :: J, i, a, ai, ia, k, ij, ij_rec, a, b, ab 
!
      real(dp), dimension(:,:), allocatable :: C_o, C_v
      real(dp), dimension(:,:), allocatable :: chol_ao_sq, chol_ao
      real(dp), dimension(:,:), allocatable :: L_ai_J
      real(dp), dimension(:,:), allocatable :: L_ij_J
      real(dp), dimension(:,:), allocatable :: L_ab_J
      real(dp), dimension(:,:), allocatable :: X
!
!     Batching variables
!
      integer(i15) :: b_batch = 0, b_first = 0, b_last = 0, b_length = 0
      integer(i15) :: required = 0, available = 0, n_batch = 0, batch_dimension = 0
      integer(i15) :: max_batch_length = 0
!
      integer(i15) :: throw_away_index = 0
      real(dp)     :: throw_away
!
      call allocator(C_o, wf%n_ao, wf%n_o)
      call allocator(C_v, wf%n_ao, wf%n_v)
!
      do i = 1, wf%n_o
!
         C_o(:,i) = wf%mo_coef_cc2_ccs(:,i)
!
      enddo
!
      do a = 1, wf%n_v
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
!     Open files for mo cholesky in CCS/CC2 block diagonal basis
!
      call generate_unit_identifier(unit_chol_mo_ai)
      call generate_unit_identifier(unit_chol_mo_ia)
      call generate_unit_identifier(unit_chol_mo_ij)
      call generate_unit_identifier(unit_chol_mo_ab)
!
!     ai-type
!
      open(unit=unit_chol_mo_ai, file='cholesky_ai_cc2_amplitudes', action='write', status='unknown', &
           access='direct', form='unformatted', recl=dp*(wf%n_J), iostat=ioerror)
!
!     ia-type
!
      open(unit=unit_chol_mo_ia, file='cholesky_ia_cc2_amplitudes', action='write', status='unknown', &
           access='direct', form='unformatted', recl=dp*(wf%n_J), iostat=ioerror)
!
!     ij-type
!
      open(unit=unit_chol_mo_ij, file='cholesky_ij_cc2_amplitudes', action='write', status='unknown', &
           access='direct', form='unformatted', recl=dp*(wf%n_J), iostat=ioerror)
!
!     ab-type
!
      open(unit=unit_chol_mo_ab, file='cholesky_ab_cc2_amplitudes', action='write', status='unknown', &
           access='direct', form='unformatted', recl=dp*(wf%n_J), iostat=ioerror)
!
      if (ioerror .ne. 0) then
            write(unit_output,*)'WARNING: error while creating cholesky_ai_cc2_amplitudes'
            stop
         endif
!
      call allocator(L_ai_J, (wf%n_o)*(wf%n_v), wf%n_J)
      L_ai_J = zero
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
         call allocator(X, wf%n_ao, wf%n_o)
!
         call dgemm('N','N',     &
                     wf%n_ao,    &
                     wf%n_o,     &
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
                     wf%n_v,     &
                     wf%n_o,     &
                     wf%n_ao,    &
                     one,        &
                     C_v,        &
                     wf%n_ao,    &
                     X,          &
                     wf%n_ao,    &
                     one,        &
                     L_ai_J(1,J),&
                     wf%n_v)
!
         call deallocator(X, wf%n_ao, wf%n_o)
      enddo
!      
!     ai-type
!
      do ai = 1, (wf%n_total_active_v)*(wf%n_total_active_o)
         write(unit_chol_mo_ai, rec=ai) (L_ai_J(ai,j), j = 1, wf%n_J)
      enddo
!
!     ia-type
!
      do a = 1, (wf%n_v)
         do i = 1, wf%n_o
            ia = index_two(i, a, wf%n_o)
            ai = index_two(a, i, wf%n_v)
            write(unit_chol_mo_ia, rec=ia) (L_ai_J(ai,j), j = 1, wf%n_J)
         enddo
      enddo

      call deallocator(L_ai_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
      call allocator(L_ij_J, (wf%n_o)*(wf%n_o), wf%n_J)
      L_ij_J = zero
!
!     Rewind ao cholesky file
!
      rewind(unit_chol_ao)
      read(unit_chol_ao,*) wf%n_ao, wf%n_J

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
         call allocator(X, wf%n_ao, wf%n_o)
!
         call dgemm('N','N',     &
                     wf%n_ao,    &
                     wf%n_o,     &
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
                     wf%n_o,     &
                     wf%n_o,     &
                     wf%n_ao,    &
                     one,        &
                     C_o,        &
                     wf%n_ao,    &
                     X,          &
                     wf%n_ao,    &
                     one,        &
                     L_ij_J(1,J),&
                     wf%n_o)
!
         call deallocator(X, wf%n_ao, wf%n_o)
      enddo
!
!     ij-type
!
      do ij = 1, (wf%n_o)*(wf%n_o)
            write(unit_chol_mo_ij, rec=ij) (L_ij_J(ij,j), j = 1, wf%n_J)
      enddo
!
      call deallocator(L_ij_J, (wf%n_o)*(wf%n_o), wf%n_J)
!
!     Cannot hold all vectors, must first write for each J, 
!     then read and write to direct access file in batches of b
!     
      call allocator(L_ab_J, wf%n_v, wf%n_v)
!
      rewind(unit_chol_ao) 
      read(unit_chol_ao,*) wf%n_ao, wf%n_J
!
!
      call generate_unit_identifier(unit_chol_mo_ab_tmp)
      open(unit_chol_mo_ab_tmp, file='cholesky_ab_tmp', status='unknown', form='unformatted')
      rewind(unit_chol_mo_ab_tmp)
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
         call allocator(X, wf%n_ao, wf%n_v)
!
         call dgemm('N','N',     &
                     wf%n_ao,    &
                     wf%n_v,     &
                     wf%n_ao,    &
                     one,        &
                     chol_ao_sq, &
                     wf%n_ao,    &
                     C_v,        &
                     wf%n_ao,    &
                     zero,       &
                     X,          &
                     wf%n_ao)
!
         call deallocator(chol_ao_sq, wf%n_ao, wf%n_ao)
!
         call dgemm('T','N',     &
                     wf%n_v,     &
                     wf%n_v,     &
                     wf%n_ao,    &
                     one,        &
                     C_v,        &
                     wf%n_ao,    &
                     X,          &
                     wf%n_ao,    &
                     zero,       &
                     L_ab_J,     &
                     wf%n_v)
!
         call deallocator(X, wf%n_ao, wf%n_v)
!
         write(unit_chol_mo_ab_tmp) ((L_ab_J(a,b), b = 1, a), a = 1, wf%n_v)
!
      enddo
!
      call deallocator(L_ab_J, (wf%n_v),(wf%n_v))
!
!     Read L_ab_J in batches over b
!
      required = ((wf%n_v)**2)*(wf%n_J)
!
      required = 4*required ! In words
      available = get_available()
!
      batch_dimension  = wf%n_v ! Batch over the virtual index b
      max_batch_length = 0      ! Initilization of unset variables 
      n_batch          = 0
!
      call num_batch(required, available, max_batch_length, n_batch, batch_dimension)           
!
!     Loop over the number of a batches 
!
      do b_batch = 1, n_batch
!  
         rewind(unit_chol_mo_ab_tmp)
!
         call batch_limits(b_first, b_last, b_batch, max_batch_length, batch_dimension)
         b_length = b_last - b_first + 1 
!
         call allocator(L_ab_J, (((b_length + 1)*b_length/2)+(wf%n_v - b_length - b_first + 1)*b_length), wf%n_J)
!
         if (b_first .ne. 1) then
!  
!           Calculate index of last element to throw away
!  
            throw_away_index = index_packed(wf%n_v, b_first - 1)
!  
!           Throw away all elements from 1 to throw_away_index, then read from batch start
!  
            do j = 1, wf%n_J
!
              read(unit_chol_mo_ab_tmp) (throw_away, i = 1, throw_away_index), &
                                    (L_ab_J(a,j), a = 1,(((b_length + 1)*b_length/2)+(wf%n_v - b_length - b_first + 1)*b_length))
!
            enddo
!
         else
!  
!           Read from the start of each entry
!  
            do j = 1, wf%n_J
!
              read(unit_chol_mo_ab_tmp) (L_ab_J(a,j), &
                        a = 1, (((b_length + 1)*b_length/2)+(wf%n_v - b_length - b_first + 1)*b_length))
!
            enddo
!
         endif
!
         do a = 1, wf%n_v
            do b = b_first, b_last
               ab = index_packed(a, b)
               write(unit_chol_mo_ab, rec=ab) (L_ab_J(ab, J), J = 1, wf%n_J)
            enddo
         enddo
!
         call deallocator(L_ab_J, (((b_length+1)*b_length/2)+(wf%n_v - b_length - b_first + 1)*b_length), wf%n_J)
!
      enddo
      close(unit_chol_mo_ab_tmp, status='delete')
!
      close(unit_chol_ao)
!
      close(unit_chol_mo_ai)
      close(unit_chol_mo_ia)
      close(unit_chol_mo_ij)
      close(unit_chol_mo_ab)
!
      call deallocator(C_o, wf%n_ao, wf%n_o)
      call deallocator(C_v, wf%n_ao, wf%n_v)
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
  subroutine read_cholesky_ia_for_cc2_amplitudes_mlccsd(wf,L_ia_J, i_first, i_last, a_first, a_last)
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
      class(mlccsd)                :: wf    
      integer(i15)             :: i_first, a_first     ! First index (can differ from 1 when batching or for mlcc) 
      integer(i15)             :: i_last, a_last      ! Last index (can differ from n_o when batching or for mlcc) 
      real(dp), dimension(:,:) :: L_ia_J ! L_ia^J
!
!     Local routine variables
!
      integer(i15) :: unit_chol_mo_ia = -1 ! Unit identifier for cholesky_ia file
      integer(i15) :: ioerror = 0
!
      integer(i15) :: i = 0, j = 0, a = 0, ia = 0, ia_full = 0
!
      integer(i15) :: i_length, a_length
!  
      i_length = i_last - i_first + 1
      a_length = a_last - a_first + 1
!
!     Prepare for reading: generate unit idientifier, open, and rewind file
!
      call generate_unit_identifier(unit_chol_mo_ia)
      open(unit=unit_chol_mo_ia, file='cholesky_ia_cc2_amplitudes', action='read', status='unknown', &
           access='direct', form='unformatted', recl=dp*(wf%n_J), iostat=ioerror)
      if (ioerror .ne. 0) then
         write(unit_output,*)'WARNING: error while reading cholesky_ij_direct'
         stop
      endif
!
!     Read Cholesky vectors into the L_ia_J matrix
!
      do i = 1, i_length
         do a = 1, a_length
!
            ia_full = index_two(i + i_first - 1, a + a_first - 1, wf%n_o)
            ia = index_two(i, a, i_length)
            read(unit_chol_mo_ia, rec=ia_full) (L_ia_J(ia,j), j = 1, wf%n_J)
!
         enddo
      enddo
!
!     Close file
!
      close(unit_chol_mo_ia)

!   
   end subroutine read_cholesky_ia_for_cc2_amplitudes_mlccsd
!
!
   subroutine read_cholesky_ij_for_cc2_amplitudes_mlccsd(wf,L_ij_J , i_first, i_last, j_first, j_last)
!!
!!    Read Cholesky IJ 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!    Reads the MO Cholesky IJ (occ-occ) vectors from file and 
!!    places them in the incoming L_ij_J matrix
!!
!!    Optional arguments: i_first, i_last, j_first, j_last can be used in order to restrict indices
!!
      implicit none
!
      class(mlccsd)            :: wf
      integer(i15)             :: i_first, j_first     ! First index (can differ from 1 when batching or for mlcc) 
      integer(i15)             :: i_last, j_last      ! Last index (can differ from n_o when batching or for mlcc)   
      real(dp), dimension(:,:) :: L_ij_J ! L_ij^J
!
!     Local routine variables 
!
      integer(i15) :: unit_chol_mo_ij = -1 ! Unit identifier for cholesky_ij file
      integer(i15) :: ioerror
!
      integer(i15) :: i = 0, j = 0, k = 0, ij = 0, ik = 0, ij_full
!
      integer(i15) :: i_length, j_length ! number of i and j elements
!
!
      i_length = i_last - i_first + 1
      j_length = j_last - j_first + 1
!
!     Prepare for reading: generate unit idientifier, open file, and rewind
!
      call generate_unit_identifier(unit_chol_mo_ij)
      open(unit=unit_chol_mo_ij, file='cholesky_ij_cc2_amplitudes', action='read', status='unknown', &
           access='direct', form='unformatted', recl=dp*(wf%n_J), iostat=ioerror)
      if (ioerror .ne. 0) then
         write(unit_output,*)'WARNING: error while reading cholesky_ij_direct'
         stop
      endif
!
!     Read the Cholesky vectors into the L_ij_J matrix
!
      do i = 1, i_length
         do k = 1, j_length
!
            ij_full = index_two((i + i_first - 1), (k + j_first - 1), wf%n_o)
!
            ij =  index_two(i, k, i_length)
!
            read(unit_chol_mo_ij, rec=ij_full) (L_ij_J(ij,J), J = 1, wf%n_J)
!
         enddo
      enddo
!
!     Close file
!
      close(unit_chol_mo_ij) 
!  
   end subroutine read_cholesky_ij_for_cc2_amplitudes_mlccsd
!
!
    subroutine read_cholesky_ab_for_cc2_amplitudes_mlccsd(wf, L_ab_J, a_first, a_last, b_first, b_last)
!!
!!    Read Cholesky AB 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2017
!!
!!    Reads the MO Cholesky AB (vir-vir) vectors from file and
!!    places them in the incoming L_ab_J matrix, with batching 
!!    if necessary
!!
!!    Optional arguments: b_first, b_last, a_first, a_last can be used in order to restrict indices
!!
      implicit none
!
      class(mlccsd)                :: wf
      integer(i15), intent(in) :: a_first, b_first   ! First index (can differ from 1 when batching  or for mlcc)
      integer(i15), intent(in) :: a_last, b_last    ! Last index  (can differ from n_v when batching or for mlcc)
      real(dp), dimension(((a_last - a_first + 1)*(b_last - b_first + 1)), wf%n_J) :: L_ab_J ! L_ab^J
!
!
      integer(i15) :: unit_chol_mo_ab_direct = -1 ! Unit identifier for cholesky_ab file
      integer(i15) :: ioerror = 0
!
      integer(i15) :: a = 0, b = 0, j = 0, i = 0, ab = 0, ab_full = 0
      integer(i15) :: a_length, b_length
!
      a_length = a_last - a_first + 1
      b_length = b_last - b_first + 1

!
!     Prepare for reading: generate unit identifier, open, and rewind file
!  
      call generate_unit_identifier(unit_chol_mo_ab_direct)
      open(unit=unit_chol_mo_ab_direct, file='cholesky_ab_cc2_amplitudes', action='read', status='unknown', &
            access='direct', form='unformatted', recl=dp*(wf%n_J), iostat=ioerror)

!
      if (ioerror .ne. 0) then
         write(unit_output,*)'WARNING: error while reading cholesky_ab_direct.', ioerror
         stop
      endif
!
         do a = 1, a_length
            do b = 1, b_length
               ab_full = index_packed(a + a_first - 1,b + b_first - 1)
               ab = index_two(a, b, a_length)
               read(unit_chol_mo_ab_direct, rec=ab_full) (L_ab_J(ab, J), J = 1, wf%n_J)
            enddo
         enddo
!
!        Close file
!        
      close(unit_chol_mo_ab_direct)
!
!
   end subroutine read_cholesky_ab_for_cc2_amplitudes_mlccsd
!
!
  subroutine get_cholesky_ai_for_cc2_amplitudes_mlccsd(wf, L_ai_J, a_first, a_last, i_first, i_last)
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
      integer(i15) :: a = 0, b = 0, J = 0, i = 0, ai = 0, Ja = 0, ia = 0
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
      real(dp), dimension(:,:), allocatable :: t1_a_i
      real(dp), dimension(:,:), allocatable :: X1_a_i
      integer(i15) :: n_active_o = 0, n_active_v = 0
!
      n_active_o = i_last - i_first + 1
      n_active_v = a_last - a_first + 1
!
      
      call allocator(X1_a_I, wf%n_v, wf%n_o)
!
      call dgemm('T','N',  &
                  wf%n_v,  &
                  wf%n_o,  &
                  wf%n_v,  &
                  one,     &
                  wf%T_v,  &
                  wf%n_v,  &
                  wf%t1am, &
                  wf%n_v,  &
                  zero,    &
                  X1_a_I,  &
                  wf%n_v)
!
      call allocator(t1_a_i, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N', &
                  wf%n_v,  &
                  wf%n_o,  &
                  wf%n_o,  &
                  one,     &
                  X1_a_I,  &
                  wf%n_v,  &
                  wf%T_o,  &
                  wf%n_o,  &
                  zero,    &
                  t1_a_i,  &
                  wf%n_v)
!
         call deallocator(X1_a_I, wf%n_v, wf%n_o)
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
         call wf%read_cholesky_ab_for_cc2_amplitudes(L_ab_J, batch_start, batch_end, 1, wf%n_v)
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
                     n_active_o,            &
                     wf%n_v,                &
                     one,                   &
                     L_Ja_b,                &
                     batch_length*(wf%n_J), &
                     t1_a_i,                &
                     wf%n_v,                &
                     one,                   &
                     L_Ja_i(L_off, 1),      &
                     (n_active_v)*(wf%n_J))
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
      call wf%read_cholesky_ij_for_cc2_amplitudes(L_ik_J, 1, n_active_o, 1, wf%n_o) ! L_ik_J(ik,J) = L_ik^J 
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
      call dgemm('N','N',              &
                  n_active_v,          &
                  n_active_o*(wf%n_J), &
                  wf%n_o,              &
                  -one,                &
                  t1_a_i,              &
                  wf%n_v,              &
                  L_k_iJ,              &
                  wf%n_o,              &
                  zero,                &
                  L_a_iJ,              &
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
      call wf%read_cholesky_ia_for_cc2_amplitudes(L_kb_J, 1, wf%n_o, 1, wf%n_v)
!
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
                  n_active_o,        &
                  wf%n_v,            &
                  one,               &
                  L_kJ_b,            &
                  (wf%n_o)*(wf%n_J), &
                  t1_a_i,            &
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
      call dgemm('N','N',              &
                  n_active_v,          &
                  n_active_o*(wf%n_J), &
                  wf%n_o,              &
                  -one,                &
                  t1_a_i,              &
                  wf%n_v,              &
                  L_k_iJ,              &
                  wf%n_o,              &
                  zero,                &
                  L_a_iJ,              &
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
      call deallocator(t1_a_i, wf%n_v, wf%n_o)
!
   end subroutine get_cholesky_ai_for_cc2_amplitudes_mlccsd
!
!
end submodule cholesky