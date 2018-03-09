submodule (sccsd_class) metric
!
!!
!!    Metric submodule (SCCSD)
!!    Written by Eirik F. Kjønstad, June 2017
!!
!!    Contains the following family of procedures of the SCCSD class:
!!
!!    calc_overlap:          calculates the generalized overlap r_A^T M r_B between 
!!                           the constrained eigenstates. For r_B -> M r_B, it uses
!!                           the metric_transformation routine.
!!
!!    metric_transformation: transforms a vector r by the SCCSD metric M: r -> M * r.
!!                           See the routine for the definition of M.
!!
!!    construct_q:           constructs the vector q_mu = < mu | e^T | R >.
!!
!!    Q_transformation:      transforms a vector r by the matrix Q_mu,nu = < mu | e^T | nu >,
!!                           or its transpose.
!!
!!    S_transformation:      transforms a vector r by the elementary overlap matrix, 
!!                           S_mu,nu = < mu | nu >. At the moment used solely for debug 
!!                           purposes.
!! 
!
   implicit none 
!
!
contains 
!
!
   module subroutine calc_overlap_sccsd(wf)
!!
!!    Calculate Overlap (SCCSD)
!!    Written by Eirik F. Kjønstad, June 2017
!!
!!    Calculates the overlap between the two constrained states,
!!    i.e., the value of r_A^T M r_B.
!!
      implicit none 
!
      class(sccsd) :: wf 
!
      real(dp), dimension(:,:), allocatable :: rA      ! (r_ai^A r_aibj^A)
      real(dp), dimension(:,:), allocatable :: rA_a_i  ! r_ai^A 
      real(dp), dimension(:,:), allocatable :: rA_aibj ! r_aibj^A
!
      real(dp), dimension(:,:), allocatable :: rB      ! (r_ai^B r_aibj^B)
      real(dp), dimension(:,:), allocatable :: rB_a_i  ! r_ai^B
      real(dp), dimension(:,:), allocatable :: rB_aibj ! r_aibj^B
!
      integer(i15) :: unit_solution = 0
      integer(i15) :: ioerror = 0
!
      real(dp) :: ddot 
!
!     Read eigenvectors from disk 
!
      call generate_unit_identifier(unit_solution)
      open(unit=unit_solution, file='right_valence', action='read', status='unknown', &
            access='direct', form='unformatted', recl=dp*(wf%n_parameters), iostat=ioerror) 
!
      if (ioerror .ne. 0) then 
!
         write(unit_output,'(t3,a)') 'Could not read excited state vectors in calc_overlap_sccsd'
         stop
!
      endif
!
!     Read state A from disk
!
      call allocator(rA, wf%n_parameters, 1)
      rA = zero
!
      read(unit_solution, rec=wf%state_A, iostat=ioerror) rA 
!
      call allocator(rA_a_i, wf%n_v, wf%n_o)
      call allocator(rA_aibj, wf%n_t2am, 1)
!
      rA_a_i  = zero
      rA_aibj = zero
!
      call dcopy(wf%n_t1am, rA, 1, rA_a_i, 1)
      call dcopy(wf%n_t2am, rA(wf%n_t1am + 1, 1), 1, rA_aibj, 1)
!
      call deallocator(rA, wf%n_parameters, 1)
!
!     Read state B from disk
!
      call allocator(rB, wf%n_parameters, 1)
      rB = zero
!
      read(unit_solution, rec=wf%state_B, iostat=ioerror) rB
!
      call allocator(rB_a_i, wf%n_v, wf%n_o)
      call allocator(rB_aibj, wf%n_t2am, 1)
!
      rB_a_i  = zero
      rB_aibj = zero
!
      call dcopy(wf%n_t1am, rB, 1, rB_a_i, 1)
      call dcopy(wf%n_t2am, rB(wf%n_t1am + 1, 1), 1, rB_aibj, 1)
!
      call deallocator(rB, wf%n_parameters, 1)
!
!     rB = M rB 
!
      call wf%metric_transformation(rB_a_i, rB_aibj)
!
!     overlap = rA^T M rB 
!
      wf%overlap = zero
      wf%overlap = ddot(wf%n_t1am, rA_a_i, 1, rB_a_i, 1) &
                 + ddot(wf%n_t2am, rA_aibj, 1, rB_aibj, 1)
!
!     Close file with eigenvectors 
!
      close(unit_solution)
!
!     Final deallocations
!
      call deallocator(rA_a_i, wf%n_v, wf%n_o)
      call deallocator(rA_aibj, wf%n_t2am, 1)   
!
      call deallocator(rB_a_i, wf%n_v, wf%n_o)
      call deallocator(rB_aibj, wf%n_t2am, 1)   
!
   end subroutine calc_overlap_sccsd
!
!
   module subroutine calc_ground_state_overlap_sccsd(wf)
!!
!!    Calculate ground state overlap (SCCSD)
!!    Written by Eirik F. Kjønstad, Feb 2018
!!
!!    Calculates the overlap between the two constrained states,
!!    i.e., the value of r_A^T M r_B, where r_A is the ground state,
!!    i.e., (1 0) in full space.
!!
!!    In this case, since r_A is known, it is more convenient to compute 
!!    the overlap directly in this routine. 
!!
!!       r_A^T M r_B = eta^T r_B + q^T S Q A r_B / (1 + q^T S q)
!!
      implicit none 
!
      class(sccsd) :: wf 
!
      real(dp), dimension(:,:), allocatable :: rB      ! (r_ai^B r_aibj^B)
      real(dp), dimension(:,:), allocatable :: rB_a_i  ! r_ai^B
      real(dp), dimension(:,:), allocatable :: rB_aibj ! r_aibj^B
!
      real(dp), dimension(:,:), allocatable :: q_a_i  ! q_ai
      real(dp), dimension(:,:), allocatable :: q_aibj ! q_aibj
!
      real(dp), dimension(:,:), allocatable :: tmp_q_a_i  ! q_ai
      real(dp), dimension(:,:), allocatable :: tmp_q_aibj ! q_aibj
!
      real(dp), dimension(:,:), allocatable :: eta ! (eta_ai eta_aibj)
!
      real(dp) :: etanorm, ccnorm 
!
      logical :: transpose
!
      integer(i15) :: unit_solution = 0
      integer(i15) :: ioerror = 0
!
      real(dp) :: ddot 
!
      wf%overlap = zero 
!
!     Construct the eta vector 
!
      call allocator(eta, wf%n_parameters, 1)
      eta = zero 
!
!     We need the amplitudes to use the Q and q routines,
!     and for constructing the eta vector 
!
      call wf%read_amplitudes
!
      call wf%construct_eta(eta)
!
!     Read eigenvectors from disk 
!
      call generate_unit_identifier(unit_solution)
      open(unit=unit_solution, file='right_valence', action='read', status='unknown', &
            access='direct', form='unformatted', recl=dp*(wf%n_parameters), iostat=ioerror) 
!
      if (ioerror .ne. 0) then 
!
         write(unit_output,'(t3,a)') 'Could not read excited state vectors in calc_ground_state_overlap_sccsd'
         stop
!
      endif
!
!     Read state B from disk
!
      call allocator(rB, wf%n_parameters, 1)
      rB = zero
!
      read(unit_solution, rec=wf%state_B, iostat=ioerror) rB
!
!     Compute the eta^T r_B contribution to the overlap 
!
      wf%overlap = wf%overlap + ddot(wf%n_parameters, eta, 1, rB, 1)
!
      etanorm = zero 
      etanorm = etanorm + ddot(wf%n_parameters, eta, 1, eta, 1)
      etanorm = sqrt(etanorm)
!
      write(unit_output,'(/t6,a26,f16.12)')  'Value of normalized eta*r:', wf%overlap/etanorm
!
      call deallocator(eta, wf%n_parameters, 1)
!
      call allocator(rB_a_i, wf%n_v, wf%n_o)
      call allocator(rB_aibj, wf%n_t2am, 1)
!
      rB_a_i  = zero
      rB_aibj = zero
!
      call dcopy(wf%n_t1am, rB, 1, rB_a_i, 1)
      call dcopy(wf%n_t2am, rB(wf%n_t1am + 1, 1), 1, rB_aibj, 1)
!
      call deallocator(rB, wf%n_parameters, 1)
!
!     r_B <- A r_B 
!
      call wf%jacobian_ccsd_transformation(rB_a_i, rB_aibj)
!
!     This transformation destructs the amplitudes,
!     therefore re-read to avoid problems 
!
      call wf%read_amplitudes
!
!     r_B <- Q A r_B 
!
      transpose = .false.
      call wf%Q_transformation(rB_a_i, rB_aibj, transpose)
!
!     r_B <- S Q A r_B 
!
      call wf%S_transformation(rB_a_i, rB_aibj)
!
!     Construct the q vector in two copies 
!
      call allocator(q_a_i, wf%n_v, wf%n_o)
      call allocator(q_aibj, wf%n_t2am, 1)
!
      q_a_i  = zero
      q_aibj = zero 
!
      call wf%construct_q(q_a_i, q_aibj)
!
      call allocator(tmp_q_a_i, wf%n_v, wf%n_o)
      call allocator(tmp_q_aibj, wf%n_t2am, 1)
!
      tmp_q_a_i  = zero
      tmp_q_aibj = zero 
!
      call wf%construct_q(tmp_q_a_i, tmp_q_aibj)
!
!     tmp_q <- S q 
!
      call wf%S_transformation(tmp_q_a_i, tmp_q_aibj)
!
!     Calculate coupled cluster norm, 1 + q^T S q 
!
      ccnorm = one + ddot(wf%n_t1am, q_a_i, 1, tmp_q_a_i, 1) &
                   + ddot(wf%n_t2am, q_aibj, 1, tmp_q_aibj, 1) ! 1 + q^T S q
!
      call deallocator(tmp_q_a_i, wf%n_v, wf%n_o)
      call deallocator(tmp_q_aibj, wf%n_t2am, 1)
!
!     Calculate the second overlap contribution, q^T S Q A r_B / (1 + q^T S q)
!
      wf%overlap = wf%overlap + (one/ccnorm)*(ddot(wf%n_t1am, q_a_i, 1, rB_a_i, 1) + &
                                              ddot(wf%n_t2am, q_aibj, 1, rB_aibj, 1))
!
!     Close file with eigenvectors 
!
      close(unit_solution)
!
!     Final deallocations 
!
      call deallocator(q_a_i, wf%n_v, wf%n_o)
      call deallocator(q_aibj, wf%n_t2am, 1)
!
      call deallocator(rB_a_i, wf%n_v, wf%n_o)
      call deallocator(rB_aibj, wf%n_t2am, 1)   
!
      call wf%destruct_amplitudes
!
   end subroutine calc_ground_state_overlap_sccsd
!
!
   module subroutine metric_transformation_sccsd(wf, r_a_i, r_aibj)
!!
!!    Metric Transformation (SCCSD)
!!    Written by Eirik F. Kjønstad, June 2017
!!
!!    Performs the SCCSD metric transformation,
!!
!!       r -> M r = (N t-bar t-bar^T + Q^T S Q 
!!             - Q^T S q t-bar^T - t-bar q^T S Q) r,
!!
!!    where t-bar is the multiplier vector, N the 
!!    coupled cluster norm, 1 + q^T S q, and
!!
!!       q_mu    = < mu | e^T | R >,
!!       Q_mu,nu = < mu | e^T | nu >,
!!       S_mu,nu = < mu | nu> (elementary overlap; that is, mu is not tilde basis)
!!
!!    Here, q and Q are in biorthonormal elementary basis. 
!!
      implicit none
!
      class(sccsd) :: wf 
!
      real(dp), dimension(wf%n_v, wf%n_o) :: r_a_i  ! r_ai 
      real(dp), dimension(wf%n_t2am, 1)   :: r_aibj ! r_aibj 
!
      real(dp), dimension(:,:), allocatable :: transformed_a_i  
      real(dp), dimension(:,:), allocatable :: transformed_aibj
!
      real(dp), dimension(:,:), allocatable :: temporary_a_i
      real(dp), dimension(:,:), allocatable :: temporary_aibj
!
      real(dp), dimension(:,:), allocatable :: temporary_two_a_i
      real(dp), dimension(:,:), allocatable :: temporary_two_aibj
!
      real(dp), dimension(:,:), allocatable :: multipliers
!
      integer(i15) :: unit_multipliers = 0
      integer(i15) :: ioerror = 0, i = 0, a = 0
!
      real(dp) :: ccnorm
      real(dp) :: dotproduct ! Holds various dot products
      real(dp) :: ddot
!
      logical :: transpose ! Q^T if true, Q if false 
!
!     We need the amplitudes to use the Q and q routines
!
      call wf%read_amplitudes
!
!     Allocate the array that holds the transformed vector,
!     and the temporary vector used for each individual contribution
!
      call allocator(temporary_a_i, wf%n_v, wf%n_o)
      call allocator(temporary_aibj, wf%n_t2am, 1)
!
      call allocator(transformed_a_i, wf%n_v, wf%n_o)
      call allocator(transformed_aibj, wf%n_t2am, 1)
!
      temporary_a_i  = zero 
      temporary_aibj = zero
!
      transformed_a_i  = zero 
      transformed_aibj = zero
!
!     :: Term 2. Q^T S Q r ::
!
!     temp = r
!     
      temporary_a_i  = r_a_i
      temporary_aibj = r_aibj
!
!     temp = Q r 
!
      transpose = .false.
      call wf%Q_transformation(temporary_a_i, temporary_aibj, transpose)
!
!     temp = S Q r 
!
      call wf%S_transformation(temporary_a_i, temporary_aibj)
!
!     temp = Q^T S Q r 
!
      transpose = .true.
      call wf%Q_transformation(temporary_a_i, temporary_aibj, transpose)
!
      call daxpy(wf%n_t1am, one, temporary_a_i, 1, transformed_a_i, 1)
      call daxpy(wf%n_t2am, one, temporary_aibj, 1, transformed_aibj, 1)
!
!     :: Term 4. - t-bar (q^T S Q r) ::
!
!     temp = r
!     
      temporary_a_i  = r_a_i
      temporary_aibj = r_aibj  
!
!     temp = Q r 
!
      transpose = .false.
      call wf%Q_transformation(temporary_a_i, temporary_aibj, transpose)
!
!     temp = S Q r 
!
      call wf%S_transformation(temporary_a_i, temporary_aibj)
!
!     dotproduct = q^T S Q r 
!
      call allocator(temporary_two_a_i, wf%n_v, wf%n_o)
      call allocator(temporary_two_aibj, wf%n_t2am, 1)
!
      temporary_two_a_i  = zero
      temporary_two_aibj = zero 
!
      call wf%construct_q(temporary_two_a_i, temporary_two_aibj)
!
      dotproduct = zero
      dotproduct = ddot(wf%n_t1am, temporary_two_a_i, 1, temporary_a_i, 1) &
                 + ddot(wf%n_t2am, temporary_two_aibj, 1, temporary_aibj, 1)
!
!     Read the multiplier vector 
!
      call generate_unit_identifier(unit_multipliers)
      open(unit=unit_multipliers, file='multipliers', action='read', status='unknown', &
      access='direct', form='unformatted', recl=dp*(wf%n_parameters), iostat=ioerror) 
!
      if (ioerror .ne. 0) then 
!
         write(unit_output,'(t3,a)') 'Could not open multipliers file in metric_transformation_sccsd'
         stop
!
      endif
!
      call allocator(multipliers, wf%n_parameters, 1) ! (t-bar_ai, t-bar_aibj)
      multipliers = zero
!
      read(unit_multipliers, rec=1, iostat=ioerror) multipliers
!
      if (ioerror .ne. 0) then 
!
         write(unit_output,'(t3,a)') 'Could not read multipliers in metric_transformation_sccsd'
         stop
!
      endif
!
      close(unit_multipliers)
!
!     Add - t-bar (q^T S Q r) to the transformed vector 
!
      call daxpy(wf%n_t1am, -dotproduct, multipliers, 1, transformed_a_i, 1)
      call daxpy(wf%n_t2am, -dotproduct, multipliers(wf%n_t1am + 1, 1), 1, transformed_aibj, 1)
!
!     :: Term 3. - Q^T S q (t-bar^T r) :: 
!
!     dotproduct = t-bar^T r 
!
      dotproduct = zero 
      dotproduct = ddot(wf%n_t1am, multipliers, 1, r_a_i, 1) &
                 + ddot(wf%n_t2am, multipliers(wf%n_t1am + 1, 1), 1, r_aibj, 1)
!
!     Recall that temporary_two holds the q-vector. 
!
!     temporary = q 
!
      temporary_a_i  = zero
      temporary_aibj = zero
!
      call dcopy(wf%n_t1am, temporary_two_a_i, 1, temporary_a_i, 1)
      call dcopy(wf%n_t2am, temporary_two_aibj, 1, temporary_aibj, 1)
!
!     temporary_two = S q 
!
      call wf%S_transformation(temporary_two_a_i, temporary_two_aibj)
!
!     Get the coupled cluster norm, N = 1 + q^T S q, for term 1
! 
      ccnorm = one + ddot(wf%n_t1am, temporary_a_i, 1, temporary_two_a_i, 1) &
                   + ddot(wf%n_t2am, temporary_aibj, 1, temporary_two_aibj, 1) ! 1 + q^T S q
!
!     temporary_two = Q^T S q 
!
      transpose = .true.
      call wf%Q_transformation(temporary_two_a_i, temporary_two_aibj, transpose)      
!
      call daxpy(wf%n_t1am, -dotproduct, temporary_two_a_i, 1, transformed_a_i, 1)
      call daxpy(wf%n_t2am, -dotproduct, temporary_two_aibj, 1, transformed_aibj, 1)
!
!     :: Term 1. N t-bar (t-bar^T r) ::
!
!     Recall: dotproduct = t-bar^T r
!
      call daxpy(wf%n_t1am, dotproduct*ccnorm, multipliers, 1, transformed_a_i, 1)
      call daxpy(wf%n_t2am, dotproduct*ccnorm, multipliers(wf%n_t1am + 1, 1), 1, transformed_aibj, 1)
!
      call deallocator(temporary_two_a_i, wf%n_v, wf%n_o)
      call deallocator(temporary_two_aibj, wf%n_t2am, 1)
!
      call deallocator(temporary_a_i, wf%n_v, wf%n_o)
      call deallocator(temporary_aibj, wf%n_t2am, 1)
!
      r_a_i  = transformed_a_i
      r_aibj = transformed_aibj
!
      call deallocator(transformed_a_i, wf%n_v, wf%n_o)
      call deallocator(transformed_aibj, wf%n_t2am, 1)
!
      call deallocator(multipliers, wf%n_parameters, 1)
!
      call wf%destruct_amplitudes
!
   end subroutine metric_transformation_sccsd
!
!
   module subroutine construct_q_sccsd(wf, q_a_i, q_aibj)
!!
!!    Construct q (SCCSD)
!!    Written by Eirik F. Kjønstad, June 2017
!!
!!    Constructs q_mu = < mu | e^T | R >.
!!
      implicit none 
!
      class(sccsd) :: wf 
!
      real(dp), dimension(wf%n_v, wf%n_o) :: q_a_i  ! q_ai 
      real(dp), dimension(wf%n_t2am, 1)   :: q_aibj ! q_aibj
!
      integer(i15) :: j = 0, b = 0, i = 0, a = 0, bj = 0, ai = 0, aibj = 0
!
      q_a_i  = zero
      q_aibj = zero 
!
!     q_ai = t_i^a 
!
      q_a_i = wf%t1am
!
!     q_aibj = (t_ij^ab + t_i^a t_j^b)/(1 + delta_ai,bj)
!
      do j = 1, wf%n_o
         do b = 1, wf%n_v
!
            bj = index_two(b, j, wf%n_v)
!
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  ai = index_two(a, i, wf%n_v)
!
                  aibj = index_packed(ai, bj)
!
                  if (a .eq. b .and. i .eq. j) then
!
                     q_aibj(aibj, 1) = half*(wf%t2am(aibj, 1) &
                                     + (wf%t1am(a,i))*(wf%t1am(b,j)))
!
                  else
!
                     q_aibj(aibj, 1) = wf%t2am(aibj, 1) &
                                     + (wf%t1am(a,i))*(wf%t1am(b,j))
!
                  endif                    
!
               enddo
            enddo
         enddo
      enddo
!
   end subroutine construct_q_sccsd    
!
!
   module subroutine Q_transformation_sccsd(wf, r_a_i, r_aibj, transpose)
!!
!!    Q transformation (SCCSD)
!!    Written by Eirik F. Kjønstad, June 2017
!!
!!    Transforms the vector r by Q_mu,nu = < mu | e^T | nu >
!!    or its transpose.
!!
      implicit none 
!
      class(sccsd) :: wf 
!
      real(dp), dimension(wf%n_v, wf%n_o) :: r_a_i  ! r_ai 
      real(dp), dimension(wf%n_t2am, 1)   :: r_aibj ! r_aibj
!
      logical :: transpose ! If true, transform by Q^T instead of Q  
!
      real(dp), dimension(:,:), allocatable :: transformed_a_i
      real(dp), dimension(:,:), allocatable :: transformed_aibj
!
      integer(i15) :: j = 0, b = 0, i = 0, a = 0, bj = 0, ai = 0, aibj = 0
      integer(i15) :: c = 0, k = 0, ck = 0, aick = 0
!
      call allocator(transformed_a_i, wf%n_v, wf%n_o)
      call allocator(transformed_aibj, wf%n_t2am, 1)
!
      transformed_a_i  = zero
      transformed_aibj = zero
!
      if (.not. transpose) then 
!
!        :: Q transformation :: 
!
         transformed_a_i = r_a_i
!
         do j = 1, wf%n_o
            do b = 1, wf%n_v
!
               bj = index_two(b, j, wf%n_v)
!
               do i = 1, wf%n_o
                  do a = 1, wf%n_v
!
                     ai = index_two(a, i, wf%n_v)
!
                     aibj = index_packed(ai, bj)
!
                     if (a .eq. b .and. i .eq. j) then
!
                        transformed_aibj(aibj, 1) = r_aibj(aibj, 1) &
                                                  + half*(r_a_i(a, i)*(wf%t1am(b, j)) &
                                                  + r_a_i(b, j)*(wf%t1am(a, i)))
!
                     else
!
                        transformed_aibj(aibj, 1) = r_aibj(aibj, 1) &
                                                  + r_a_i(a, i)*wf%t1am(b, j) &
                                                  + r_a_i(b, j)*wf%t1am(a, i) 
!
                     endif                  
!
                  enddo
               enddo
            enddo
         enddo
!
      else
!
!        :: Q^T transformation :: 
!
         transformed_aibj = r_aibj
!
         transformed_a_i = r_a_i
!
         do i = 1, wf%n_o
            do a = 1, wf%n_v
!
!              Do the sum over c & k: add sum_ck t_k^c r_aick
!
               ai = index_two(a, i, wf%n_v)
!
               do k = 1, wf%n_o
                  do c = 1, wf%n_v
!
                     ck = index_two(c, k, wf%n_v)
!
                     aick = index_packed(ai, ck)
!
                     transformed_a_i(a, i) = transformed_a_i(a, i) &
                                           + r_aibj(aick, 1)*(wf%t1am(c, k))
!
                  enddo
               enddo
!
            enddo
         enddo
!
      endif
!
!     Copy the transformed vector into the incoming vector
!
      r_a_i  = transformed_a_i
      r_aibj = transformed_aibj
!
      call deallocator(transformed_a_i, wf%n_v, wf%n_o)
      call deallocator(transformed_aibj, wf%n_t2am, 1)
!
   end subroutine Q_transformation_sccsd
! 
!
   module subroutine S_transformation_sccsd(wf, r_a_i, r_aibj)
!!
!!    S Transformation (SCCSD)
!!    Written by Eirik F. Kjønstad, June 2017
!!
!!    Calculates the transformation by the elementary overlap matrix S_mu,nu < mu | nu >,
!!    and places the result in the incoming vector r.
!!
!!    Note: used for debugging purposes (t-bar^T = q^T S Q / N in the FCI limit). 
!!
      class(sccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o) :: r_a_i 
      real(dp), dimension(wf%n_t2am, 1)   :: r_aibj
!
      real(dp), dimension(:,:), allocatable :: transformed_a_i
      real(dp), dimension(:,:), allocatable :: transformed_aibj
!
      integer(i15) :: i = 0, a = 0, j = 0, b = 0, ai = 0, bj = 0, aj = 0, bi = 0, aibj = 0, ajbi = 0
      real(dp) :: bnf
!
      call allocator(transformed_a_i, wf%n_v, wf%n_o)
      call allocator(transformed_aibj, wf%n_t2am, 1)
!
      transformed_a_i = zero 
      transformed_aibj = zero
!
!     Calculate the singles transformation
!
      call dcopy(wf%n_t1am,r_a_i,1,transformed_a_i,1)
      call dscal(wf%n_t1am,two,transformed_a_i,1)
!
!     Calculate the doubles transformation
!
      do i = 1,wf%n_o
         do a = 1,wf%n_v
!
            ai = index_two(a, i, wf%n_v)
!
            do j = 1,wf%n_o
               do b = 1,wf%n_v
!
                  bj = index_two(b, j, wf%n_v)
                  aj = index_two(a, j, wf%n_v)
                  bi = index_two(b, i, wf%n_v)
!
                  aibj = index_packed(ai,bj)
                  ajbi = index_packed(aj,bi)
!
!                 Compute the bi-orthonormal factor
!
                  if (i .eq. j .and. a .eq. b) then
                     bnf = two
                  else
                     bnf = one
                  endif
!
                  transformed_aibj(aibj,1) = two*bnf*(two*r_aibj(aibj,1)-r_aibj(ajbi,1))
!
               enddo
            enddo
         enddo
      enddo
!
!     Copy the transformed result into the vec: vec <- S * vec
!
      call dcopy(wf%n_t1am,transformed_a_i,1,r_a_i,1)
      call dcopy(wf%n_t2am,transformed_aibj,1,r_aibj,1)
!
   end subroutine S_transformation_sccsd
!
!
end submodule metric
