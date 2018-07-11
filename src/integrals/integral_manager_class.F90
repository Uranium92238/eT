module integral_manager_class
!
!!
!!    Integral_manager class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
!  Fortran interfaces to C++ routines
!
   use h_xy
   use s_xy
   use g_wxyz
   use L_wxyz
   use libint_initialization
!
!  Disk & memory class modules
!
   use index
   use reordering
   use array_utilities
   use array_analysis
!
   use file_class
   use memory_manager_class
!
   use molecular_system_class
!
   implicit none
!
!  Class definition
!
   type :: integral_manager
!
!     No attributes yet
!
   contains
!
      procedure :: cholesky_decompose => cholesky_decompose_integral_manager
!
      procedure :: get_ao_h_xy   => get_ao_h_xy_integral_manager   ! h_αβ
      procedure :: get_ao_s_xy   => get_ao_s_xy_integral_manager   ! s_αβ
      procedure :: get_ao_g_wxyz => get_ao_g_wxyz_integral_manager ! g_αβγδ
!
   end type integral_manager
!
!
!   interface
!
!
!      module subroutine cholesky_decompose_integral_manager(integrals, molecule)
!!
!!       Cholesky decompose
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!         implicit none
!
!         class(integral_manager) :: integrals
!         class(molecular_system) :: molecule
!
!      end subroutine cholesky_decompose_integral_manager
!
!
!   end interface
!
!
contains
!
!
   subroutine cholesky_decompose_integral_manager(integrals, molecule)
!!
!!    Cholesky decompose
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    In this routine, the electronic repulsion AO integrals g_αβ,γδ = (αβ|γδ)
!!    are Cholesky decomposed, yielding a set of Cholesky vectors L_αβ,J, J = 1,
!!    2, 3, ..., n_J. From these, the integrals can be reconstructed as
!!
!!       g_αβ,γδ = sum_J L_αβ,J L_αβ,J^T.
!!
!!    This reconstruction of g_αβ,γδ should not be performed directly by the developer.
!!    Specialized get-routines for integrals in MO basis is implemented, and AO integrals
!!    could be calculated directly using the libint-interface routine get_ao_g_wxyz.
!!
!!    The algorithm in the routine is based on the algorithm given on p. 333-334,
!!    Chapter 13, "Cholesky Decomposition Techniques in Electronic Structure Theory"
!!    (Aquilante et al.), in the monograph "Linear-Scaling Techniques in Computational
!!    Chemistry and Physics" (2011, R. Zalesny et al.)
!!
!!    To prepare, the diagonal terms, D_αβ = g_αβ,αβ, are screened in shell-pairs (AB),
!!    where shell-pairs for which all diagonals are less than a given threshold (e.g., 10-8)
!!    are ignored. Afterwards, a number (<= max_qual) of diagonals are qualified to be decomposed,
!!    Max_qual is 100 by default, but fewer AO pairs are qualified if there are less than 100
!!    diagonal elements which satisfy the qualification criterion, D_αβ > span * max(D_αβ). The
!!    ordering of the qualified diagonals are from the largest to the smallest diagonal within
!!    each shell pair, which themselves are ordered according to their maximum element (from
!!    largest to smallest). The span factor is by deafult 10-3. From the qualified AO pairs {γδ*},
!!    a set of Cholesky vectors are formed as follows:
!!
!!       L_αβ^J = (αβ|γδ*)/sqrt(D_γδ*),
!!
!!    where γδ* refers to the AO index pair that corresponds to a qualified AO pair. Each qualified
!!    AO pair gives rise to one Cholesky vector (labeled by J = 1, 2, 3, ...).
!!
!!    A qualified diagonal will only give rise to a Cholesky vector if the corresponding diagonal is
!!    above 10-8 after updating the diagonal elements using the larger diagonals among the qualified
!!    AO pairs. In particular, for a batch of qualified AO pairs, the diagonal is updated as
!!
!!       D_αβ <- D_αβ - (L_αβ^J)^2,
!!
!!    implying that certain qualified AO pairs γδ* might be ignored, since the associated diagonal elements
!!    D_γδ* could become less than the given threshold (e.g., 1.0D-8). The non-diagonal elements are similarly
!!    updated as
!!
!!       (αβ|γδ*) <- (αβ|γδ*) - L_αβ^J L_γδ*^J, (1)
!!
!!    where γδ* is restricted by the qualified AO pairs and αβ runs over all the significant AO pairs
!!    in the given cycle.
!!
!!    In each new cycle, the previously made Cholesky vectors are used to subtract from (αβ|γδ),
!!    according to (1), and the new qualified AO pairs are used as described above. The αβ index
!!    is restricted to significant diagonals (e.g., >= 1.0D-8), which is reduced from cycle to
!!    cycle. Only the significant parts {αβ} of the previous L_αβ^J vectors are held in memory,
!!    and the end result is a Cholesky basis: a set of AO pairs {γδ} that are to be used to construct
!!    the full set of Cholesky vectors.
!!
      implicit none
!
      class(integral_manager) :: integrals
      class(molecular_system) :: molecule
!
      type(file) :: screening_info, cholesky_ao_vectors, auxiliary, auxiliary_inverse
!
      integer(i15) :: n_s, n_sp
      integer(i15) :: n_cholesky, n_new_cholesky
      integer(i15) :: n_new_sig_sp, n_new_sig_aop
      integer(i15) :: sig_sp_counter
      integer(i15) :: n_sig_sp, n_sig_aop, n_qual_sp, sig_AB_sp
      integer(i15) :: n_old_qual_aop, n_qual_aop, n_qual_aop_in_sp,  n_qual_aop_in_prev_sps
      integer(i15) :: n_sp_in_basis, sp_in_basis, iteration
!
      integer(i15) :: A, B, C, D, AB_sp, CD_sp, AB
      integer(i15) :: I, K, J, L, KL, I_counter, J_counter, IJ
      integer(i15) :: w, x, y, z, xy, yz, wx, xy_packed, wx_packed, wxyz, yz_packed
      integer(i15) :: first, last, first_sig_aop, last_sig_aop, first_x, first_y
      integer(i15) :: sp, aop, current_qual, current_sig_sp, qual, qual_max, current_new_sig_sp
      integer(i15) :: current_aop_in_sp, n_left_out, n_vectors, size_AB, offset_yz, new_position, rec_offset
!
      type(interval) :: A_interval, B_interval, C_interval, D_interval
!
      logical, dimension(:, :), allocatable :: sig_sp, new_sig_sp, first_sig_sp
!
      real(dp), dimension(:,:), allocatable :: D_AB, D_xy, D_xy_new
      real(dp), dimension(:,:), allocatable :: g_AB_AB, g_AB_CD, g_CD_AB, g_wxyz, g_J_yz, L_K_yz
      real(dp), dimension(:,:), allocatable :: cholesky, cholesky_new
      real(dp), dimension(:,:), allocatable :: cholesky_copy
      real(dp), dimension(:,:), allocatable :: cholesky_tmp
      real(dp), dimension(:,:), allocatable :: sorted_max_in_sig_sp
      real(dp), dimension(:,:), allocatable :: sorted_qual_aop_in_sp
      real(dp), dimension(:,:), allocatable :: max_in_sig_sp
      real(dp), dimension(:,:), allocatable :: integrals_auxiliary
      real(dp), dimension(:,:), allocatable :: integrals_auxiliary_packed
      real(dp), dimension(:,:), allocatable :: integrals_auxiliary_new
      real(dp), dimension(:,:), allocatable :: auxiliary_basis_new
      real(dp), dimension(:,:), allocatable :: auxiliary_basis
      real(dp), dimension(:,:), allocatable :: auxiliary_basis_inverse
      real(dp), dimension(:,:), allocatable :: D_diff
!
      integer(i15), dimension(:,:), allocatable :: sig_sp_to_first_sig_aop
      integer(i15), dimension(:,:), allocatable :: new_sig_sp_to_first_sig_aop
      integer(i15), dimension(:,:), allocatable :: new_sig_aop_to_aos
      integer(i15), dimension(:,:), allocatable :: sig_aop_to_aos
      integer(i15), dimension(:,:), allocatable :: max_in_sig_sp_indices
      integer(i15), dimension(:,:), allocatable :: qual_aop, qual_sp
      integer(i15), dimension(:,:), allocatable :: qual_aop_copy, qual_sp_copy
      integer(i15), dimension(:,:), allocatable :: sorted_qual_aop_in_sp_indices
      integer(i15), dimension(:,:), allocatable :: sorted_max_sig_sp
      integer(i15), dimension(:,:), allocatable :: sig_sp_to_shells
      integer(i15), dimension(:,:), allocatable :: new_sig_sp_to_shells
      integer(i15), dimension(:,:), allocatable :: cholesky_basis, cholesky_basis_new
      integer(i15), dimension(:,:), allocatable :: sig_sp_to_first_sig_sp
      integer(i15), dimension(:,:), allocatable :: basis_shell_info_full
      integer(i15), dimension(:,:), allocatable :: basis_shell_info
      integer(i15), dimension(:,:), allocatable :: basis_aops_in_CD_sp, basis_aops_in_AB_sp
      integer(i15), dimension(:,:), allocatable :: keep_vectors, leave_out
!
      logical :: construct_more_choleskys, done, found
!
      real(dp) :: D_max, max_in_sp, D_max_full, max_diff, ddot
!
      real(dp), parameter :: threshold = 1.0D-8
      real(dp), parameter :: span      = 1.0D-2
!
      integer(i15), parameter :: max_qual = 100
!
      real(dp) :: s_prep_time, e_prep_time, s_select_basis_time, e_select_basis_time
      real(dp) :: s_build_basis_time, e_build_basis_time, s_decomp_time, e_decomp_time, full_decomp_time
      real(dp) :: e_invert_time, s_invert_time, e_build_vectors_time, s_build_vectors_time
      real(dp) :: s_integral_time, e_integral_time, full_integral_time
      real(dp) :: s_reduce_time, e_reduce_time, full_reduce_time
      real(dp) :: s_construct_time, e_construct_time, full_construct_time
!
      write(output%unit, '(/a51)') ':: Cholesky decomposition of two-electron integrals'
      flush(output%unit)
!
      write(output%unit, '(/a)') ' - Preparations'
      flush(output%unit)
!
      call cpu_time(s_prep_time)
!
      write(output%unit, '(a20, i10)')'Number of aos:      ', molecule%get_n_aos()
      write(output%unit, '(a20, i10)')'Number of ao pairs: ', &
                                          molecule%get_n_aos()*(molecule%get_n_aos()+1)/2
!
      n_s   = molecule%get_n_shells() ! Number of shells
      n_sp  = n_s*(n_s + 1)/2         ! Number of shell pairs packed
!
!     Pre-screening of full diagonal
!
      allocate(sig_sp(n_sp, 1))
      sig_sp = .false.
!
      sp = 1        ! Shell pair number
      n_sig_aop = 0 ! Number of significant AO pairs
      n_sig_sp  = 0 ! Number of significant shell pairs
!
      do B = 1, n_s
         do A = B, n_s
!
            A_interval = molecule%get_shell_limits(A)
            B_interval = molecule%get_shell_limits(B)
!
!           Construct diagonal D_AB for the given shell pair
!
            call mem%alloc(g_AB_AB, &
                     (A_interval%size)*(B_interval%size), &
                     (A_interval%size)*(B_interval%size))
!
            g_AB_AB = zero
            call integrals%get_ao_g_wxyz(g_AB_AB, A, B, A, B)
!
            call mem%alloc(D_AB, (A_interval%size)*(B_interval%size), 1)
!
            do xy = 1, (A_interval%size)*(B_interval%size)
!
               D_AB(xy, 1) = g_AB_AB(xy, xy)
!
            enddo
!
            call mem%dealloc(g_AB_AB, &
                     (A_interval%size)*(B_interval%size), &
                     (A_interval%size)*(B_interval%size))
!
!           Determine whether shell pair is significant
!
            sig_sp(sp, 1) = is_significant(D_AB, (A_interval%size)*(B_interval%size), threshold)
!
            call mem%dealloc(D_AB, (A_interval%size)*(B_interval%size), 1)
!
            if (sig_sp(sp, 1)) then
!
               n_sig_aop = n_sig_aop + &
                              get_size_sp(A_interval, B_interval)
!
               n_sig_sp = n_sig_sp + 1
!
            endif
!
            sp = sp + 1
!
         enddo
      enddo
!
      first_sig_sp = sig_sp
!
      write(output%unit, '(/a)')'Initial reduction of shell pairs:'
      write(output%unit, '(a33, 2x, i6)')'Total number of shell pairs:     ', n_sp
      write(output%unit, '(a33, 2x, i6)')'Significant shell pairs:         ', n_sig_sp
!
!     Construct significant diagonal
!
      call mem%alloc_int(sig_sp_to_first_sig_aop, n_sig_sp + 1, 1)
      sig_sp_to_first_sig_aop = 0
      sig_sp_to_first_sig_aop(n_sig_sp + 1, 1) = n_sig_aop + 1
!
      call mem%alloc_int(sig_sp_to_shells, n_sig_sp, 2) ! A B
      sig_sp_to_shells = 0
!
      call mem%alloc(D_xy, n_sig_aop, 1)
      D_xy = zero
!
      call mem%alloc_int(sig_aop_to_aos, n_sig_aop, 2)
      sig_aop_to_aos = 0
!
!     Note: allocated with length n_significant_sp + 1, last element is used for n_significant_aop
!     This is convenient because significant_sp_to_first_significant_aop will be used to calculate lengths.
!
      sp              = 1
      current_sig_sp  = 1
      first_sig_aop   = 1
!
      do B = 1, n_s
         do A = B, n_s
!
            if (sig_sp(sp, 1)) then
!
               sig_sp_to_first_sig_aop(current_sig_sp, 1) = first_sig_aop
!
               A_interval = molecule%get_shell_limits(A)
               B_interval = molecule%get_shell_limits(B)
!
               sig_sp_to_shells(current_sig_sp, 1) = A
               sig_sp_to_shells(current_sig_sp, 2) = B
!
               call mem%alloc(g_AB_AB, &
                     (A_interval%size)*(B_interval%size), &
                     (A_interval%size)*(B_interval%size))
!
               g_AB_AB = zero
               call integrals%get_ao_g_wxyz(g_AB_AB, A, B, A, B)
!
               if (A .eq. B) then
!
                  do x = 1, A_interval%size
                     do y = 1, B_interval%size
!
                        xy = A_interval%size*(y - 1) + x
                        xy_packed = (max(x,y)*(max(x,y)-3)/2) + x + y
!
                        D_xy(xy_packed + first_sig_aop - 1, 1) = g_AB_AB(xy, xy)
!
                        sig_aop_to_aos(xy_packed + first_sig_aop - 1, 1) = &
                                                      A_interval%first + x - 1
!
                        sig_aop_to_aos(xy_packed + first_sig_aop - 1, 2) = &
                                                      B_interval%first + y - 1
!
                     enddo
                  enddo
!
               else ! A ≠ B
!
                  do x = 1, (A_interval%size)
                     do y = 1, (B_interval%size)
!
                        xy = A_interval%size*(y - 1) + x
                        D_xy(xy + first_sig_aop - 1, 1) = g_AB_AB(xy,xy)
!
                        sig_aop_to_aos(xy + first_sig_aop - 1, 1) = A_interval%first + x - 1
                        sig_aop_to_aos(xy + first_sig_aop - 1, 2) = B_interval%first + y - 1
!
                     enddo
                  enddo
!
               endif
!
               call mem%dealloc(g_AB_AB, &
                     (A_interval%size)*(B_interval%size), &
                     (A_interval%size)*(B_interval%size))
!
               first_sig_aop = first_sig_aop + get_size_sp(A_interval, B_interval)
!
               current_sig_sp = current_sig_sp + 1
!
            endif ! End of if (significant)
!
            sp = sp + 1
!
         enddo
      enddo
!
!     Write screening_info_file containing 
!
!        1. number of significant shell pairs, number of significant ao pairs
!        2. sig_sp - vector of logicals to describe which shell pairs are significant 
!        3. D_xy = ( xy | xy ), the significant diagonal.
!
      call screening_info%init('screening_info', 'sequential', 'formatted')
      call disk%open_file(screening_info, 'write')
!
      write(screening_info%unit,*) n_sig_sp, n_sig_aop
      write(screening_info%unit,*) sig_sp
      write(screening_info%unit,*) D_xy
!
      call disk%close_file(screening_info)
!
!     Timings
!
      call cpu_time(e_prep_time)
      write(output%unit, '(/a17, e11.4, a5)')'Time to prepare: ', e_prep_time - s_prep_time, ' sec.'
      flush(output%unit)
!
!     Determining the basis
!
      write(output%unit, '(/a)') ' - Determinig the elements of the basis'
      flush(output%unit)
!
      call cpu_time(s_select_basis_time)
!
      n_cholesky = 0
      write(output%unit, '(/a)') 'Iter.    Max diagonal    #Qualified    #Cholesky    Cholesky array size'
      write(output%unit, '(a)') '----------------------------------------------------------------------'
      flush(output%unit)
!
      iteration = 0
!
      full_integral_time = 0
      full_reduce_time = 0
      full_construct_time = 0
!
      do while (.not. done)
!
      iteration = iteration + 1
!
!        Shell maximums and shell maximums indices vectors
!
         call mem%alloc(max_in_sig_sp, n_sig_sp, 1)
!
         max_in_sig_sp = zero
!
         do sp = 1, n_sig_sp
!
!           Get first and last indices of shell pair
!
            first = sig_sp_to_first_sig_aop(sp, 1)
            last  = sig_sp_to_first_sig_aop(sp + 1, 1) - 1
!
!           Determine the largest elements
!
            do I = first, last
!
               if (D_xy(I, 1) .gt. max_in_sig_sp(sp, 1)) then
!
                  max_in_sig_sp(sp, 1) = D_xy(I, 1)
!
               elseif ((D_xy(I, 1) .lt. 0.0d0) .and. (abs(D_xy(I, 1)) .gt. 1.0d-10)) then
!
                  write(output%unit, '(a44, e12.5)')'Error: found significant negative diagonal: ',D_xy(I, 1)
                  stop 
!
               endif
!
            enddo
!
         enddo
!
!        Sort from largest to smallest by determining an index array
!
         call mem%alloc_int(sorted_max_sig_sp, n_sig_sp,1)
         sorted_max_sig_sp = 0
!
         call mem%alloc(sorted_max_in_sig_sp, n_sig_sp, 1)
         sorted_max_in_sig_sp = zero
!
         call get_n_highest(n_sig_sp, n_sig_sp, max_in_sig_sp, sorted_max_in_sig_sp, sorted_max_sig_sp)
!
         D_max_full  = sorted_max_in_sig_sp(1, 1)
         n_qual_aop  = 0
         n_qual_sp   = 0
!
         call mem%dealloc(sorted_max_in_sig_sp, n_sig_sp, 1)
!
         call mem%alloc_int(qual_aop, max_qual, 3)
         call mem%alloc_int(qual_sp, n_sp, 3)
         qual_sp = 0
!
         do sp = 1, n_sig_sp
!
            current_sig_sp = sorted_max_sig_sp(sp, 1)
!
            first_sig_aop = sig_sp_to_first_sig_aop(current_sig_sp, 1)
            last_sig_aop  = sig_sp_to_first_sig_aop(current_sig_sp + 1, 1) - 1
!
            n_qual_aop_in_sp = 0
!
            do aop = first_sig_aop, last_sig_aop
!
               if ((D_xy(aop, 1) .ge. span*D_max_full ) .and. (n_qual_aop .lt. max_qual)) then
!
                  n_qual_aop_in_sp  = n_qual_aop_in_sp + 1
                  n_qual_aop        = n_qual_aop + 1
!
               endif
!
            enddo
!
            if (n_qual_aop_in_sp .ne. 0) then
!
               n_qual_sp = n_qual_sp + 1
!
               call mem%alloc_int(sorted_qual_aop_in_sp_indices, n_qual_aop_in_sp, 1)
               call mem%alloc(sorted_qual_aop_in_sp, n_qual_aop_in_sp, 1)
!
               call get_n_highest(n_qual_aop_in_sp, last_sig_aop - first_sig_aop + 1, &
                                 D_xy(first_sig_aop:last_sig_aop, 1), sorted_qual_aop_in_sp, &
                                 sorted_qual_aop_in_sp_indices)
!
               n_old_qual_aop = (n_qual_aop - n_qual_aop_in_sp)
!
               do aop = 1, n_qual_aop_in_sp
!
                  qual_aop(aop + n_old_qual_aop, 1) = sig_aop_to_aos(sorted_qual_aop_in_sp_indices(aop, 1) &
                                                               + first_sig_aop - 1, 1)
!
                  qual_aop(aop + n_old_qual_aop, 2) = sig_aop_to_aos(sorted_qual_aop_in_sp_indices(aop, 1) &
                                                               + first_sig_aop - 1, 2)
!
                  qual_aop(aop + n_old_qual_aop, 3) = sorted_qual_aop_in_sp_indices(aop, 1) &
                                                               + first_sig_aop - 1
!
               enddo
!
               first_x = sig_aop_to_aos(first_sig_aop, 1) ! alpha
               first_y = sig_aop_to_aos(first_sig_aop, 2) ! beta
!
               qual_sp(n_qual_sp, 1) = molecule%basis2shell(first_x)
               qual_sp(n_qual_sp, 2) = molecule%basis2shell(first_y)
               qual_sp(n_qual_sp, 3) = n_qual_aop_in_sp
!
               call mem%dealloc_int(sorted_qual_aop_in_sp_indices, n_qual_aop_in_sp, 1)
               call mem%dealloc(sorted_qual_aop_in_sp, n_qual_aop_in_sp, 1)
!
            endif
!
            if (n_qual_aop == max_qual) then
!
               exit
!
            endif
!
         enddo
!
         call mem%dealloc_int(sorted_max_sig_sp, n_sig_sp, 1)
!
!        Cut out the qualified parts of the aop and sp lists
!
         call mem%alloc_int(qual_aop_copy, n_qual_aop, 3)
         call mem%alloc_int(qual_sp_copy, n_qual_sp, 3)
!
         qual_aop_copy(:, :) = qual_aop(1 : n_qual_aop, :)
         qual_sp_copy(:, :)  = qual_sp(1 : n_qual_sp, :)
!
         call mem%dealloc_int(qual_aop, max_qual, 3)
         call mem%dealloc_int(qual_sp, n_sp, 3)
!
         call mem%alloc_int(qual_aop, n_qual_aop, 3)
         call mem%alloc_int(qual_sp, n_qual_sp, 3)
!
         qual_aop    = qual_aop_copy
         qual_sp     = qual_sp_copy
!
         call mem%dealloc_int(qual_aop_copy, n_qual_aop, 3)
         call mem%dealloc_int(qual_sp_copy, n_qual_sp, 3)
!
         n_qual_aop_in_prev_sps = 0
!
         call cpu_time(s_integral_time)
!
         call mem%alloc(g_wxyz, n_sig_aop, n_qual_aop)
!
         do CD_sp = 1, n_qual_sp
!
            C                = qual_sp(CD_sp, 1)
            D                = qual_sp(CD_sp, 2)
            n_qual_aop_in_sp = qual_sp(CD_sp, 3)
!
            C_interval = molecule%get_shell_limits(C)
            D_interval = molecule%get_shell_limits(D)
!
!           Calculate the ({wx} | J) integrals,
!           where {wx} is the screened list of integrals
!
            AB_sp                  = 1
            sig_AB_sp              = 1
!
            do B = 1, n_s
               do A = B, n_s
!
                  if (sig_sp(AB_sp, 1)) then
!
                     A_interval = molecule%get_shell_limits(A)
                     B_interval = molecule%get_shell_limits(B)
!
                     call mem%alloc(g_AB_CD, &
                                    (A_interval%size)*(B_interval%size), &
                                    (C_interval%size)*(D_interval%size))
!
                     call integrals%get_ao_g_wxyz(g_AB_CD, A, B, C, D)
!
                     do aop = 1, n_qual_aop_in_sp
!
                        y = qual_aop(aop + n_qual_aop_in_prev_sps, 1)
                        z = qual_aop(aop + n_qual_aop_in_prev_sps, 2)
!
                        yz = C_interval%size*(z - D_interval%first) + y - C_interval%first + 1
!
                        if (A == B) then
!
                           do w = 1, A_interval%size
                              do x = w, B_interval%size
!
                                 wx_packed = (max(w,x)*(max(w,x)-3)/2) + w + x
                                 wx = A_interval%size*(x-1) + w
!
                                 g_wxyz(sig_sp_to_first_sig_aop(sig_AB_sp, 1) + wx_packed - 1, aop + n_qual_aop_in_prev_sps) &
                                       = g_AB_CD(wx, yz)
!
                              enddo
                           enddo
!
                        else
!
                           do w = 1, A_interval%size
                              do x = 1, B_interval%size
!
                                 wx = A_interval%size*(x-1) + w
!
                                 g_wxyz(sig_sp_to_first_sig_aop(sig_AB_sp, 1) + wx - 1, aop + n_qual_aop_in_prev_sps) &
                                       = g_AB_CD(wx, yz)
!
                              enddo
                           enddo
!
                        endif
!
                     enddo
!
                     call mem%dealloc(g_AB_CD, &
                                       (A_interval%size)*(B_interval%size), &
                                       (C_interval%size)*(D_interval%size))
!
                     sig_AB_sp = sig_AB_sp + 1
!
                  endif
!
                  AB_sp = AB_sp + 1
!
               enddo ! A
            enddo ! B
!
            n_qual_aop_in_prev_sps = n_qual_aop_in_prev_sps + n_qual_aop_in_sp
!
         enddo ! cd_sp 
!
         call cpu_time(e_integral_time)
         full_integral_time = full_integral_time + e_integral_time - s_integral_time
!
!        Subtract old cholesky vectors
!
         call cpu_time(s_construct_time)
!
         if (n_cholesky .ne. 0) then
!
            call mem%alloc(cholesky_tmp, n_cholesky, n_qual_aop)
!
            do K = 1, n_qual_aop
!
               cholesky_tmp(:, K) = cholesky(qual_aop(K, 3), :)
!
            enddo
!
            call dgemm('N', 'N',          &
                        n_sig_aop,        &
                        n_qual_aop,       &
                        n_cholesky,       &
                        -one,             &
                        cholesky,         &
                        n_sig_aop,        &   
                        cholesky_tmp,     &
                        n_cholesky,       &
                        one,              &
                        g_wxyz,           &
                        n_sig_aop)

!
            call mem%dealloc(cholesky_tmp, n_cholesky, n_qual_aop)
!
            call mem%alloc_int(cholesky_basis, n_cholesky + n_qual_aop, 3)
            cholesky_basis(1 : n_cholesky, :) = cholesky_basis_new(:, :)
            call mem%dealloc_int(cholesky_basis_new, n_cholesky, 3)
!
         else
!
            call mem%alloc_int(cholesky_basis, n_qual_aop, 3)
            cholesky_basis = 0
!
         endif
!
         call mem%alloc(cholesky_new, n_sig_aop, n_qual_aop)
         cholesky_new = zero
!
         current_qual = 0
!
         construct_more_choleskys = .true.
!
         do while ((current_qual .lt. n_qual_aop) .and. construct_more_choleskys)
!
            current_qual = current_qual + 1
!
            D_max = zero
!
            do qual = 1, n_qual_aop
!
               xy = qual_aop(qual, 3)
!
               if (D_xy(xy, 1) .gt. D_max) then
!
                  qual_max = qual
                  D_max    = D_xy(xy, 1)
!
               endif
!
            enddo
!
            if ((D_max .gt. threshold) .and. (D_max .ge. span*D_max_full)) then
!
               cholesky_basis(n_cholesky + current_qual, 1) = qual_aop(qual_max, 1)
               cholesky_basis(n_cholesky + current_qual, 2) = qual_aop(qual_max, 2)
!
               A = molecule%basis2shell(qual_aop(qual_max, 1))
               B = molecule%basis2shell(qual_aop(qual_max, 2))
!
               cholesky_basis(n_cholesky + current_qual, 3) = get_sp_from_shells(A, B, n_s)
!
               cholesky_new(: , current_qual) = g_wxyz(:, qual_max)
!
               if (current_qual .gt. 1) then
!
                  call mem%alloc(cholesky_tmp, 1, current_qual - 1)
   !
                  cholesky_tmp(1, :) = cholesky_new(qual_aop(qual_max, 3), 1 : current_qual - 1)
!
                  call dgemm('N', 'T',                      &
                           n_sig_aop,                       &
                           1,                               &
                           current_qual - 1,                &
                           -one,                            &
                           cholesky_new,                    &
                           n_sig_aop,                       &
                           cholesky_tmp,                    &
                           1,                               &
                           one,                             &
                           cholesky_new(1, current_qual),   &
                           n_sig_aop) 
!
               call mem%dealloc(cholesky_tmp, 1, current_qual - 1)  
!
            endif
!
            g_wxyz(qual_aop(qual_max, 3), :) = zero
!
            call dscal(n_sig_aop, one/sqrt(D_max), cholesky_new(1, current_qual), 1)
!
               do xy = 1, n_sig_aop
!
                  D_xy(xy, 1) = D_xy(xy, 1) - cholesky_new(xy, current_qual)**2
!
               enddo
!
            else
!
               construct_more_choleskys = .false.
               current_qual = current_qual - 1
!
            endif
!
         enddo
!
         call cpu_time(e_construct_time)
         full_construct_time = full_construct_time + e_construct_time - s_construct_time
         call mem%dealloc(g_wxyz, n_sig_aop, n_qual_aop)
!
         n_new_cholesky = current_qual
!
!        Find new significant diagonals
!
         n_new_sig_sp = 0
         allocate(new_sig_sp(n_sig_sp,1))
!
         new_sig_sp = .false.
!
         sig_sp_counter = 0
!
         do sp = 1, n_sp
!
            if (sig_sp(sp, 1)) then
!
               sig_sp_counter = sig_sp_counter + 1
!
               first = sig_sp_to_first_sig_aop(sig_sp_counter, 1)
               last  = sig_sp_to_first_sig_aop(sig_sp_counter + 1, 1) - 1
!
               new_sig_sp(sig_sp_counter, 1) = is_significant(D_xy(first:last, 1), &
                                                last - first + 1, threshold)
!
               sig_sp(sp, 1) = new_sig_sp(sig_sp_counter, 1)
!
               if (new_sig_sp(sig_sp_counter, 1)) then
!
                  n_new_sig_sp = n_new_sig_sp + 1
!
               endif
!
            endif
!
         enddo
!
         call cpu_time(s_reduce_time)
!
         if (n_new_sig_sp .gt. 0) then
!
!           Update index lists: sps -> aops, aops -> aos, and sps -> full sps           
!
            call mem%alloc_int(new_sig_sp_to_first_sig_aop, n_new_sig_sp + 1, 1)
            new_sig_sp_to_first_sig_aop = 0
!
            call mem%alloc_int(sig_sp_to_first_sig_sp, n_sig_sp + 1, 1) ! 1 2 3 4 ... n_sig_sp, n_sig_sp + 1
            sig_sp_to_first_sig_sp(n_sig_sp + 1, 1) = n_sig_sp + 1
!
            current_new_sig_sp    = 1
            n_new_sig_aop = 0
            first_sig_aop = 1
!
            do sp = 1, n_sig_sp
!
               sig_sp_to_first_sig_sp(sp, 1) = sp
!
               if (new_sig_sp(sp, 1)) then
!
                  A = sig_sp_to_shells(sp, 1)
                  B = sig_sp_to_shells(sp, 2)
!
                  A_interval = molecule%get_shell_limits(A)
                  B_interval = molecule%get_shell_limits(B)
!
                  new_sig_sp_to_first_sig_aop(current_new_sig_sp, 1) = first_sig_aop
!
                  first_sig_aop = first_sig_aop + get_size_sp(A_interval, B_interval)
                  n_new_sig_aop = first_sig_aop - 1
!
                  current_new_sig_sp    = current_new_sig_sp + 1

               endif
!
            enddo
!
            new_sig_sp_to_first_sig_aop(current_new_sig_sp, 1) = n_new_sig_aop + 1
!
            call mem%alloc_int(new_sig_aop_to_aos, n_new_sig_aop, 2)
!
            call reduce_array_int(sig_aop_to_aos,       &
                               new_sig_aop_to_aos,      &
                               sig_sp_to_first_sig_aop, &
                               new_sig_sp,              &
                               n_sig_sp,                &
                               n_sig_aop,               &
                               n_new_sig_aop,           &
                               2)
!
            call mem%alloc_int(new_sig_sp_to_shells, n_new_sig_sp, 2)
            new_sig_sp_to_shells = 0
!
            call reduce_array_int(sig_sp_to_shells,        &
                                  new_sig_sp_to_shells,    &
                                  sig_sp_to_first_sig_sp,  &
                                  new_sig_sp,              &
                                  n_sig_sp,                &
                                  n_sig_sp,                &
                                  n_new_sig_sp,            &
                                  2)
!
            call mem%dealloc_int(sig_sp_to_first_sig_sp, n_sig_sp + 1, 1)
            call mem%dealloc_int(sig_sp_to_shells, n_sig_sp, 2)
            call mem%alloc_int(sig_sp_to_shells, n_new_sig_sp, 2)
!
            sig_sp_to_shells = new_sig_sp_to_shells
            call mem%dealloc_int(new_sig_sp_to_shells, n_new_sig_sp, 2)
!
            call mem%alloc(D_xy_new, n_new_sig_aop, 1)
!
           call reduce_vector(D_xy,                     &
                             D_xy_new,                  &
                             sig_sp_to_first_sig_aop,   &
                             new_sig_sp,                &
                             n_sig_sp,                  &
                             n_sig_aop,                 &
                             n_new_sig_aop)
!
            call mem%dealloc(D_xy, n_sig_aop, 1)
            call mem%alloc(D_xy, n_new_sig_aop, 1)
!
            call dcopy(n_new_sig_aop, D_xy_new, 1, D_xy, 1)
!
            call mem%dealloc(D_xy_new, n_new_sig_aop, 1)  
!
            if (n_cholesky == 0) then
!
               call mem%alloc(cholesky, n_new_sig_aop, n_new_cholesky)
!
               call reduce_array(cholesky_new,              &
                                  cholesky,                 &
                                  sig_sp_to_first_sig_aop,  &
                                  new_sig_sp,               &
                                  n_sig_sp,                 &
                                  n_sig_aop,                &
                                  n_new_sig_aop,            &
                                  n_new_cholesky)
!
               call mem%dealloc(cholesky_new, n_sig_aop, n_qual_aop)
!
            else
!
               call mem%alloc(cholesky_copy, n_new_sig_aop, n_new_cholesky + n_cholesky)
!
               call reduce_array(cholesky,                  &
                                  cholesky_copy,            &
                                  sig_sp_to_first_sig_aop,  &
                                  new_sig_sp,               &
                                  n_sig_sp,                 &
                                  n_sig_aop,                &
                                  n_new_sig_aop,            &
                                  n_cholesky)
!
               call reduce_array(cholesky_new,                       &
                                  cholesky_copy(1, n_cholesky + 1),  &
                                  sig_sp_to_first_sig_aop,           &
                                  new_sig_sp,                        &
                                  n_sig_sp,                          &
                                  n_sig_aop,                         &
                                  n_new_sig_aop,                     &
                                  n_new_cholesky)
!
               call mem%dealloc(cholesky, n_sig_aop, n_cholesky)
               call mem%alloc(cholesky, n_new_sig_aop, n_cholesky + n_new_cholesky)
               cholesky = cholesky_copy
               call mem%dealloc(cholesky_copy, n_new_sig_aop, n_cholesky + n_new_cholesky)
!
               call mem%dealloc(cholesky_new, n_sig_aop, n_qual_aop)
!
            endif
!
            call mem%alloc_int(cholesky_basis_new, n_cholesky + n_new_cholesky, 3)
            cholesky_basis_new(:, :) = cholesky_basis(1 : n_cholesky + n_new_cholesky, :)
            call mem%dealloc_int(cholesky_basis, n_cholesky + n_qual_aop, 3)
!
   !        Deallocate old lists & reallocate + copy over new lists
!
            deallocate(new_sig_sp)
!
            call mem%dealloc_int(sig_sp_to_first_sig_aop, n_sig_sp + 1, 1)
            call mem%alloc_int(sig_sp_to_first_sig_aop, n_new_sig_sp + 1, 1)
            sig_sp_to_first_sig_aop = new_sig_sp_to_first_sig_aop
            call mem%dealloc_int(new_sig_sp_to_first_sig_aop, n_new_sig_sp + 1, 1)
!
            call mem%dealloc_int(sig_aop_to_aos, n_sig_aop, 2)
            call mem%alloc_int(sig_aop_to_aos, n_new_sig_aop, 2)
            sig_aop_to_aos = new_sig_aop_to_aos
            call mem%dealloc_int(new_sig_aop_to_aos, n_new_sig_aop, 2)
!
            n_sig_sp = n_new_sig_sp
            n_sig_aop = n_new_sig_aop
!
            n_cholesky = n_cholesky + n_new_cholesky
!
            call mem%dealloc_int(qual_aop, n_qual_aop, 3)
            call mem%dealloc_int(qual_sp, n_qual_sp, 3)
!
            write(output%unit, '(i4, 5x,  e12.5, 4x, i4, 8x, i7, 6x, i10)') &
            iteration, D_max_full , n_qual_aop, n_cholesky, n_cholesky*n_sig_aop 
            flush(output%unit)
!
         else
!
            call mem%dealloc(cholesky_new, n_sig_aop, n_qual_aop)
            call mem%dealloc(cholesky, n_sig_aop, n_cholesky)
!
            call mem%alloc_int(cholesky_basis_new, n_cholesky + n_new_cholesky, 3)
            cholesky_basis_new(:, :) = cholesky_basis(1 : n_cholesky + n_new_cholesky, :)
            call mem%dealloc_int(cholesky_basis, n_cholesky + n_qual_aop, 3)
!
            n_cholesky = n_cholesky + n_new_cholesky
!
            done = .true.
!
            write(output%unit, '(i4, 5x,  e12.5, 4x, i4, 8x, i7, 6x, i10)') iteration, D_max_full, n_qual_aop, n_cholesky, 0 
            flush(output%unit)
!
         endif
         call cpu_time(e_reduce_time)
         full_reduce_time = full_reduce_time + e_reduce_time - s_reduce_time
!
      enddo
      write(output%unit, '(a)') '----------------------------------------------------------------------'
      flush(output%unit)
!
!     Timings
!
      call cpu_time(e_select_basis_time)
      write(output%unit, '(/a22, e11.4, a5)')'Time to select basis: ',&
                            e_select_basis_time - s_select_basis_time, ' sec.'
      write(output%unit, '(t6, a36, e11.4, a5)')'Time to construct integrals: ',&
                            full_integral_time, ' sec.'
      write(output%unit, '(t6, a36, e11.4, a5)')'Time to reduce arrays:       ',&
                            full_reduce_time, ' sec.'
      write(output%unit, '(t6, a36, e11.4, a5)')'Time to make vectors:        ',&
                            full_construct_time, ' sec.'

!
!     Building the auxiliary_basis
!
      write(output%unit, '(/a)') '- Building auxiliary basis'
      flush(output%unit)
      call cpu_time(s_build_basis_time)
!
!     Construct a list of all shell pairs (and shells) that contain elements of the basis 
!     and how many elements of the basis they contain 
!
      call mem%alloc_int(basis_shell_info_full, n_sp, 4) ! A, B, AB, n_basis_aops_in_sp
      basis_shell_info_full = 0
!
      n_sp_in_basis = 0
!
      do i = 1, n_cholesky
!
         A = molecule%basis2shell(cholesky_basis_new(i, 1))
         B = molecule%basis2shell(cholesky_basis_new(i, 2))
!
         AB = get_sp_from_shells(A, B, n_s)
!
         found = .false.
!
         do sp_in_basis = 1, n_sp_in_basis
!
            if (AB == basis_shell_info_full(sp_in_basis, 3)) then
               found = .true.
               basis_shell_info_full(sp_in_basis, 4) = basis_shell_info_full(sp_in_basis, 4) + 1 
               exit
            endif
!
         enddo
!
         if(.not. found) then
!
            n_sp_in_basis = n_sp_in_basis + 1
!
            basis_shell_info_full(n_sp_in_basis, 1) = A
            basis_shell_info_full(n_sp_in_basis, 2) = B
            basis_shell_info_full(n_sp_in_basis, 3) = AB         
            basis_shell_info_full(n_sp_in_basis, 4) = 1         
!
         endif
!
      enddo
!
      call mem%alloc_int(basis_shell_info, n_sp_in_basis, 4)
      basis_shell_info(:, :) = basis_shell_info_full(1:n_sp_in_basis, :)
      call mem%dealloc_int(basis_shell_info_full, n_sp, 4)
!
!     Construct integrals (J | J')
!     
      call mem%alloc(integrals_auxiliary_packed, n_cholesky*(n_cholesky+1)/2, 1)
!
      do CD_sp = 1, n_sp_in_basis
!
         C = basis_shell_info(CD_sp, 1)
         D = basis_shell_info(CD_sp, 2)
!
         C_interval = molecule%get_shell_limits(C)
         D_interval = molecule%get_shell_limits(D)
!
         call mem%alloc_int(basis_aops_in_CD_sp, basis_shell_info(CD_sp, 4), 3)
!
!        Determine which elements in the shell pair CD are elements of the basis
!
         current_aop_in_sp = 0
!
         do I = 1, n_cholesky
            if (cholesky_basis_new(I,3) == basis_shell_info(CD_sp, 3)) then
!
               current_aop_in_sp = current_aop_in_sp + 1
!
               basis_aops_in_CD_sp(current_aop_in_sp, 1) = cholesky_basis_new(I,1) - C_interval%first + 1
               basis_aops_in_CD_sp(current_aop_in_sp, 2) = cholesky_basis_new(I,2) - D_interval%first + 1
               basis_aops_in_CD_sp(current_aop_in_sp, 3) = I
!
            endif
         enddo
!
         do AB_sp = 1, n_sp_in_basis
!
            A = basis_shell_info(AB_sp, 1)
            B = basis_shell_info(AB_sp, 2)
!
            A_interval = molecule%get_shell_limits(A)
            B_interval = molecule%get_shell_limits(B)
!
            call mem%alloc_int(basis_aops_in_AB_sp, basis_shell_info(AB_sp, 4), 3)
!
!           Determine which elements in the shell pair AB are elements of the basis
!
            current_aop_in_sp = 0
!
            do I = 1, n_cholesky
               if (cholesky_basis_new(I,3) == basis_shell_info(AB_sp, 3)) then
!
                  current_aop_in_sp = current_aop_in_sp + 1
!
                  basis_aops_in_AB_sp(current_aop_in_sp, 1) = cholesky_basis_new(I,1) - A_interval%first + 1
                  basis_aops_in_AB_sp(current_aop_in_sp, 2) = cholesky_basis_new(I,2) - B_interval%first + 1
                  basis_aops_in_AB_sp(current_aop_in_sp, 3) = I
!
               endif
            enddo
!
!           Construct integrals
!
            call mem%alloc(g_AB_CD, &
                     (A_interval%size)*(B_interval%size), &
                     (C_interval%size)*(D_interval%size))
!
            call integrals%get_ao_g_wxyz(g_AB_CD, A, B, C, D)
!
!           Only keep those that correspond to elements of the basis
!
            do I = 1, basis_shell_info(AB_sp, 4)
               do J = 1, basis_shell_info(CD_sp, 4)
!
                  y = basis_aops_in_CD_sp(J, 1)
                  z = basis_aops_in_CD_sp(J, 2)
                  yz = C_interval%size*(z-1)+y
!
                  w = basis_aops_in_AB_sp(I, 1)
                  x = basis_aops_in_AB_sp(I, 2)
                  wx = A_interval%size*(x-1)+w
!
                  K = basis_aops_in_AB_sp(I, 3)
                  L = basis_aops_in_CD_sp(J, 3)
                  KL = (max(K,L)*(max(K,L)-3)/2) + K + L
!
                  integrals_auxiliary_packed(KL, 1) = g_AB_CD(wx, yz)
!
               enddo
            enddo
!
            call mem%dealloc(g_AB_CD, &
                     (A_interval%size)*(B_interval%size), &
                     (C_interval%size)*(D_interval%size))
!
            call mem%dealloc_int(basis_aops_in_AB_sp, basis_shell_info(AB_sp, 4), 3)       
!
         enddo ! AB
!
         call mem%dealloc_int(basis_aops_in_CD_sp, basis_shell_info(CD_sp, 4), 3)
!
      enddo ! CD
!
      call mem%dealloc_int(basis_shell_info, n_sp_in_basis, 4)
!
      call mem%alloc(auxiliary_basis, n_cholesky, n_cholesky)
!
      call mem%alloc(integrals_auxiliary, n_cholesky, n_cholesky)
      call squareup(integrals_auxiliary_packed, integrals_auxiliary, n_cholesky)
      call mem%dealloc(integrals_auxiliary_packed, n_cholesky*(n_cholesky + 1)/2, 1)
!
!     Cholesky decomposition of (J | J') to find auxiliary basis using full pivoting
!
      n_vectors = 0
      call mem%alloc_int(keep_vectors, n_cholesky, 1)
!
      call cpu_time(s_decomp_time)
!
      call full_cholesky_decomposition_effective(integrals_auxiliary, auxiliary_basis, &
                                          n_cholesky, n_vectors, threshold, keep_vectors)
!
      call cpu_time(e_decomp_time)
      full_decomp_time = e_decomp_time - s_decomp_time
!
      call mem%dealloc(integrals_auxiliary, n_cholesky, n_cholesky) 
!
!     There might be a reduction in the number of Cholesky vectors at this point (n_vectors < n_cholesky)
!     We only keep the rows (and columns) corresponding to the diagonals used for the decomposition,
!     also the cholesky vectors should be ordered as a lower triangular matrix for effective inversion.
!     The cholesky_basis array is updated.
!
      call mem%alloc(auxiliary_basis_new, n_vectors, n_vectors)
      call mem%alloc_int(cholesky_basis, n_vectors, 3)
!
      do I = 1, n_cholesky
!
         found = .false.
!
         do K = 1, n_vectors
            if (keep_vectors(K, 1) == I) then
!
               found = .true.
               exit
!
            endif
         enddo
!
         if (found) then
!
            cholesky_basis(K, :) = cholesky_basis_new(I, :)
            auxiliary_basis_new(K, :) = auxiliary_basis(I, 1 : n_vectors)
!
         endif

      enddo

      call mem%dealloc_int(cholesky_basis_new, n_cholesky, 3)
      call mem%dealloc_int(keep_vectors, n_cholesky, 1)
!
      call mem%dealloc(auxiliary_basis, n_cholesky, n_cholesky)
      call mem%alloc(auxiliary_basis, n_vectors, n_vectors)
      auxiliary_basis = auxiliary_basis_new
      call mem%dealloc(auxiliary_basis_new, n_vectors, n_vectors)
!
      n_cholesky = n_vectors
!
      write(output%unit, '(t6, a34, i7)')'Final number of cholesky vectors: ', n_cholesky
!
!     Update the basis_shell_info array which contains information of which shell pairs (and shells)
!     contain elements of the basis and how many elements of the basis they contain. 
!
      call mem%alloc_int(basis_shell_info_full, n_sp, 4) ! A, B, AB, n_basis_aops_in_sp
      basis_shell_info_full = 0
!
      n_sp_in_basis = 0
!
      do i = 1, n_cholesky
!
         A = molecule%basis2shell(cholesky_basis(i, 1))
         B = molecule%basis2shell(cholesky_basis(i, 2))
!
         AB = get_sp_from_shells(A, B, n_s)
!
         found = .false.
!
         do sp_in_basis = 1, n_sp_in_basis
!
              if (AB == basis_shell_info_full(sp_in_basis, 3)) then
                 found = .true.
                 basis_shell_info_full(sp_in_basis, 4) = basis_shell_info_full(sp_in_basis, 4) + 1 
                 exit
              endif
!
         enddo
!
         if(.not. found) then
!
            n_sp_in_basis = n_sp_in_basis + 1
!
            basis_shell_info_full(n_sp_in_basis, 1) = A
            basis_shell_info_full(n_sp_in_basis, 2) = B
            basis_shell_info_full(n_sp_in_basis, 3) = AB         
            basis_shell_info_full(n_sp_in_basis, 4) = 1         
!
         endif
!
      enddo
!
      call mem%alloc_int(basis_shell_info, n_sp_in_basis, 4)
      basis_shell_info(:, :) = basis_shell_info_full(1:n_sp_in_basis, :)
      call mem%dealloc_int(basis_shell_info_full, n_sp, 4)
!
      call cpu_time(e_build_basis_time)
      write(output%unit, '(/a21, e11.4, a5)')'Time to build basis: ',&
                            e_build_basis_time - s_build_basis_time, ' sec.'
      write(output%unit, '(t6, a25, e11.4, a5)')'Time to decompose (J|K): ',&
                            full_decomp_time, ' sec'  
!
!     Write auxiliary basis to file               
!
      call auxiliary%init('auxiliary_basis', 'sequential', 'unformatted')
      call disk%open_file(auxiliary, 'write')
!
      write(auxiliary%unit) auxiliary_basis
!
      call disk%close_file(auxiliary)
!
      write(output%unit, '(/a)') '- Inverting auxiliary basis'
      flush(output%unit)
      call cpu_time(s_invert_time)
!
      call mem%alloc(auxiliary_basis_inverse, n_cholesky, n_cholesky)
      call inv_lower_tri(auxiliary_basis_inverse, auxiliary_basis, n_cholesky)
      call mem%dealloc(auxiliary_basis, n_cholesky, n_cholesky)
!
!     Write inverse of auxiliary basis to file               
!
      call auxiliary_inverse%init('auxiliary_basis_inverse', 'sequential', 'unformatted')
      call disk%open_file(auxiliary_inverse, 'write')
!
      write(auxiliary_inverse%unit) auxiliary_basis_inverse
!
      call disk%close_file(auxiliary_inverse)
!
!     Timings
!
      call cpu_time(e_invert_time)
      write(output%unit, '(/a16, e11.4, a5)')'Time to invert: ',&
                            e_invert_time - s_invert_time, ' sec.' 
!
      write(output%unit, '(/a)') '- Construct Cholesky vectors'
      flush(output%unit)
      call cpu_time(s_build_vectors_time)
      full_integral_time = 0
      full_construct_time = 0
!
!     Prepare file for AO Cholesky vectors
!
      call cholesky_ao_vectors%init('cholesky_ao_xy', 'direct', 'unformatted', dp*n_cholesky)
!
      rec_offset = 0
!
!     Construct (K | yz) and do matrix multiplication 
!     sum_K (K | J)^-1 (J | yz) in batches of yz
!
      done = .false.
!
      do while (.not. done)
!
         call cpu_time(s_integral_time) 
!
!        Determine size of batch
!   
         sp = 0
         size_AB = 0
!
         do B = 1, n_s
            do A = B, n_s
!
               sp = sp + 1
!
               A_interval = molecule%get_shell_limits(A)
               B_interval = molecule%get_shell_limits(B)
!
               if (first_sig_sp(sp, 1)) then
!
                  size_AB = size_AB + get_size_sp(A_interval, B_interval)
!
                  if (2*size_AB*n_cholesky .gt. mem%available) then
!
                     size_AB = size_AB - get_size_sp(A_interval, B_interval)
                     sp = sp - 1
!
                  endif
!
               endif
!
            enddo
         enddo
!
!        Construct g_J_yz = (J | yz)
!
         call mem%alloc(g_J_yz, n_cholesky, size_AB)
         g_J_yz = zero
         offset_yz = 0
!
         do B = 1, n_s 
            do A = B, n_s
!
               AB_sp = get_sp_from_shells(A, B, n_s)
!

               if (first_sig_sp(AB_sp, 1) .and. AB_sp .le. sp) then
!
                  first_sig_sp(AB_sp, 1) = .false.
!
                  A_interval = molecule%get_shell_limits(A)
                  B_interval = molecule%get_shell_limits(B)
!
                  do CD_sp = 1, n_sp_in_basis

!
                     C = basis_shell_info(CD_sp, 1)
                     D = basis_shell_info(CD_sp, 2)
!
                     C_interval = molecule%get_shell_limits(C)
                     D_interval = molecule%get_shell_limits(D)
!
                     call mem%alloc_int(basis_aops_in_CD_sp, basis_shell_info(CD_sp, 4), 3)
!
                     current_aop_in_sp = 0
!
                     do I = 1, n_cholesky
                        if (cholesky_basis(I,3) == basis_shell_info(CD_sp, 3)) then
!
                           current_aop_in_sp = current_aop_in_sp + 1
!
                           basis_aops_in_CD_sp(current_aop_in_sp, 1) = cholesky_basis(I,1) - C_interval%first + 1
                           basis_aops_in_CD_sp(current_aop_in_sp, 2) = cholesky_basis(I,2) - D_interval%first + 1
                           basis_aops_in_CD_sp(current_aop_in_sp, 3) = I
!
                        endif
                     enddo
!
                     call mem%alloc(g_CD_AB, &
                              (C_interval%size)*(D_interval%size), &
                              (A_interval%size)*(B_interval%size))
!
                     call integrals%get_ao_g_wxyz(g_CD_AB, C, D, A, B)
!
                     if (A == B) then
!
                        do J = 1, basis_shell_info(CD_sp, 4)
                           w = basis_aops_in_CD_sp(J, 1)
                           x = basis_aops_in_CD_sp(J, 2)
                           L = basis_aops_in_CD_sp(J, 3)
!
                           wx = C_interval%size*(x-1)+w
!
                           do y = 1, A_interval%size
                              do z = 1, B_interval%size
                                 
!
                                    yz_packed = (max(y,z)*(max(y,z)-3)/2) + y + z
                                    yz = A_interval%size*(z-1) + y
!
                                    g_J_yz(L, yz_packed + offset_yz) = g_CD_AB(wx, yz)
!
                                 enddo
                              enddo
                           enddo
!
                        else
!
                           do J = 1, basis_shell_info(CD_sp, 4)
                              w = basis_aops_in_CD_sp(J, 1)
                              x = basis_aops_in_CD_sp(J, 2)
                              L = basis_aops_in_CD_sp(J, 3)
!
                              wx = C_interval%size*(x-1) + w
!
                              do y = 1, A_interval%size
                                 do z = 1, B_interval%size
!
                                    yz = A_interval%size*(z-1) + y
!
                                    g_J_yz(L, yz + offset_yz) = g_CD_AB(wx, yz)
!
                                 enddo
                              enddo
                           enddo
!
                        endif
!
                        call mem%dealloc(g_CD_AB,                 &
                           (C_interval%size)*(D_interval%size),   &
                           (A_interval%size)*(B_interval%size))
!
                     call mem%dealloc_int(basis_aops_in_CD_sp, basis_shell_info(CD_sp, 4), 3)
!
                  enddo ! CD
!
                  offset_yz = offset_yz + get_size_sp(A_interval, B_interval)
!
               endif
!
            enddo ! B
!        
         enddo ! A
!
         call cpu_time(e_integral_time) 
         full_integral_time = full_integral_time + e_integral_time - s_integral_time
!
         call cpu_time(s_construct_time)
!
!        L_K_yz = sum_K (K | J)^-1 (J | yz)
!
         call mem%alloc(L_K_yz, n_cholesky, size_AB)
!
         call dgemm('N', 'N',                     &
                       n_cholesky,                &
                       size_AB,                   &
                       n_cholesky,                &
                       one,                       &
                       auxiliary_basis_inverse,   &
                       n_cholesky,                &
                       g_J_yz,                    &
                       n_cholesky,                &
                       zero,                      &
                       L_K_yz,                    &
                       n_cholesky)
!
         call mem%dealloc(g_J_yz, n_cholesky, size_AB)
         call cpu_time(e_construct_time)
         full_construct_time = full_construct_time + e_construct_time - s_construct_time
!
!        Write vectors to file
!   
         call disk%open_file(cholesky_ao_vectors, 'write')
!
         do I = 1, size_AB
!
            write(cholesky_ao_vectors%unit, rec=I + rec_offset) (L_K_yz(J, I), J = 1, n_cholesky)
!
         enddo
!
         rec_offset = rec_offset + size_AB
!
         call disk%close_file(cholesky_ao_vectors)
!
         done = .true.
!
         do I = 1, n_sp
            if (first_sig_sp(I, 1)) then
               done = .false.
               exit
            endif
         enddo
!
      enddo ! done
!
!     Timings
!
      call cpu_time(e_build_vectors_time)
      write(output%unit, '(/a23, e11.4, a5)')'Time to build vectors: ',&
                            e_build_vectors_time - s_build_vectors_time, ' sec.' 
      write(output%unit, '(t6, a36, e11.4, a5)')'Time to construct integrals: ',&
                            full_integral_time, ' sec.'
      write(output%unit, '(t6, a36, e11.4, a5)')'Time to make vectors:        ',&
                            full_construct_time, ' sec.'

!
      call mem%dealloc(auxiliary_basis_inverse, n_cholesky, n_cholesky)
!
!     Test diagonal
!
      call disk%open_file(screening_info, 'read')
!
      read(screening_info%unit, *) n_sig_sp, n_sig_aop
!
      call mem%alloc(D_diff, n_sig_aop, 1)
!
      read(screening_info%unit, *) sig_sp
      read(screening_info%unit, *) D_diff
!
      call disk%close_file(screening_info)
!
      call mem%alloc(L_K_yz, n_cholesky, 1)
!
      call disk%open_file(cholesky_ao_vectors, 'read')
!
      do aop = 1, n_sig_aop
!
         read(cholesky_ao_vectors%unit, rec=aop) (L_K_yz(J, 1), J = 1, n_cholesky)
!
         D_diff(aop, 1) = D_diff(aop, 1) - ddot(n_cholesky, L_K_yz, 1, L_K_yz, 1)
!
      enddo
!
      call disk%close_file(cholesky_ao_vectors)
!
      call mem%dealloc(L_K_yz, n_cholesky, 1)
!
      max_diff = zero
!
      do aop = 1, n_sig_aop
         if (abs(D_diff(aop, 1)) .gt. max_diff) max_diff = abs(D_diff(aop, 1))
      enddo
!
      write(output%unit, '(a60, e12.4)')'Maximal difference between approximate and actual diagonal: ', max_diff
      
      call mem%dealloc(D_diff, n_sig_aop, 1)
!
   end subroutine cholesky_decompose_integral_manager
!
!
   integer(i15) function get_size_sp(A_interval, B_interval)
!!
!!    Get size shell pair
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Returns size of diagonal for a given shell pair.
!!
      implicit none
!
      type(interval) :: A_interval
      type(interval) :: B_interval
!
      if (A_interval%first == B_interval%first) then
!
         get_size_sp = A_interval%size*(A_interval%size + 1)/2
!
      else
!
         get_size_sp = (A_interval%size)*(B_interval%size)
!
      endif
!
   end function get_size_sp
!
!
   integer(i15) function get_sp_from_shells(s1, s2, n_s)
!!
!!    Get shell pair from shells,
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      integer(i15), intent(in) :: s1, s2, n_s
      integer(i15) :: A, B
!
      get_sp_from_shells = 1
!
      do B = 1, n_s
         do A = B, n_s
!
            if (s1 == A .and. s2 == B)  return
!
            get_sp_from_shells = get_sp_from_shells + 1
!
         enddo
      enddo

!
   end function get_sp_from_shells
!
!
   subroutine get_ao_h_xy_integral_manager(int, h)
!!
!!    Get h_αβ integrals
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Fortran wrapper for the C++ routine that calculates and
!!    saves the h_αβ integral_manager in the array h.
!!
      implicit none
!
      class(integral_manager) :: int
!
      real(kind=8), dimension(:,:), intent(inout) :: h
!
      call get_ao_h_xy(h)
!
   end subroutine get_ao_h_xy_integral_manager
!
!
   subroutine get_ao_s_xy_integral_manager(int, s)
!!
!!    Get s_αβ integrals
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Fortran wrapper for the C++ routine that calculates and
!!    saves the s_αβ integral_manager in the array s.
!!
      implicit none
!
      class(integral_manager) :: int
!
      real(kind=8), dimension(:,:), intent(inout) :: s
!
      call get_ao_s_xy(s)
!
   end subroutine get_ao_s_xy_integral_manager
!
!
   subroutine get_ao_g_wxyz_integral_manager(int, g, s1, s2, s3, s4)
!!
!!    Get g_αβγδ integrals
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Fortran wrapper for the C++ routine that calculates and
!!    saves the g_αβγδ integral_manager in the array g.
!!
!!    s1 is first shell index and so on.
!!
      implicit none
!
      class(integral_manager) :: int
!
      real(kind=8), dimension(:,:), intent(inout) :: g
!
      integer(kind=8), intent(in) :: s1, s2, s3, s4
!
      call get_ao_g_wxyz(g, s1, s2, s3, s4)
!
   end subroutine get_ao_g_wxyz_integral_manager
!
!
end module integral_manager_class
