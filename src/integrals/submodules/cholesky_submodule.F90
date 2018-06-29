submodule (integral_manager_class) cholesky
!
!!
!!    Cholesky submodule
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!
   use index
   use array_utilities
   use array_analysis
!
   implicit none
!
!
contains
!
!
   module subroutine cholesky_decompose_integral_manager(integrals, molecule)
!!
!!    Cholesky decompose
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    The routine has the following index mappings:
!!
!!    ....
!!
!!
      implicit none
!
      class(integral_manager) :: integrals
      class(molecular_system) :: molecule
!
      integer(i15) :: n_s, n_sp
!
      integer(i15) :: A, B
!
      real(dp), parameter :: threshold = 1.0D-8
      real(dp), parameter :: span      = 1.0D-3
!
      integer(i15), parameter :: max_qual = 100
!
      n_s   = molecule%get_n_shells() ! Number of shells
      n_sp  = n_s*(n_s + 1)/2         ! Number of shell pairs packed
!
!     Pre-screening of full diagonal
!
      allocate(sig_sp(n_sp, 1))
      sig_sp = .false.
!
      sp = 1                ! Shell pair number
      n_sig_aop = 0 ! Number of significant AO pairs
      n_sig_sp  = 0 ! Number of significant shell pairs
!
      do B = 1, n_s
         do A = B, n_s
!
            A_interval = get_shell_limits(A)
            B_interval = get_shell_limits(B)
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
            sig_sp(sp) = is_significant(D_AB, (A_interval%size)*(B_interval%size), threshold))
!
            call mem%dealloc(D_AB, (A_interval%size)*(B_interval%size), 1)
!
            if (sig_sp(sp)) then
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
!     Construct significant diagonal
!
      call mem%alloc_int(sig_sp_to_first_sig_aop, n_sig_sp + 1, 1)
      sig_sp_to_first_sig_aop = 0
      sig_sp_to_first_sig_aop(n_sig_sp + 1, 1) = n_sig_aop + 1
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
            if (sig_sp(sp)) then
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
                        D_xy(xy_packed + first_sig_aop - 1, 1) = g_wxyz(xy, xy)
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
                        D_xy(xy + first_sig_aop - 1, 1) = g_wxyz(xy,xy)
!
                        sig_aop_to_aos(xy + first_sig_aop - 1, 1) = A_interval%first + x - 1
                        sig_aop_to_aos(xy + first_sig_aop - 1, 2) = B_interval%first + y - 1
!
                     enddo
                  enddo
!
               endif
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
      n_cholesky = 0
!
      do while (.not. done)
!
!        Shell maximums and shell maximums indices vectors
!
         call mem%alloc_int(max_in_sig_sp_indices, n_sig_sp, 1)
         call mem%alloc(max_in_sig_sp, n_sig_sp, 1)
!
         max_in_sp = zero
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
               if (diag_xy(I, 1) .gt. max_in_sig_sp(sp, 1)) then
!
                  max_in_sig_sp(sp, 1)       = D_xy(I, 1)
                  max_in_sig_sp_indices(sp)  = I
!
               endif
!
            enddo
!
         enddo
!
!        Sort from largest to smallest by determining an index array
!
         call mem%alloc(sorted_max_sig_sp, n_sig_sp,1)
         sorted_max_sig_sp = 0
!
         call mem%alloc(sorted_max_in_sig_sp, n_sig_sp, 1)
         sorted_max_in_sig_sp = zero
!
         call get_n_highest(n_sig_sp, n_sig_sp, max_in_sig_sp, sorted_max_in_sig_sp, sorted_max_siq_sp)
!
         D_max       = sorted_max_in_siq_sp(1, 1)
         n_qual_aop  = 0
         n_qual_sp   = 0
!
         call mem%alloc_int(qual_aop, max_qual, 3)
         call mem%alloc_int(qual_sp_to_shells, n_sp, 3)
!
         do sp = 1, n_sig_sp
!
            current_sig_sp = sorted_max_sig_sp(sp, 1)
!
            first_sig_aop  = sig_sp_to_first_sig_aop(current_sig_sp, 1)
            last_sig_aop   = sig_sp_to_first_sig_aop(current_sig_sp + 1, 1) - 1
!
            n_qual_aop_in_sp = 0
!
            do aop = first_sig_aop, last_sig_aop
!
               if ((D_xy(aop, 1) .ge. span*D_max ) .and. (n_qual_aop .lt. max_qual)) then
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
               qual_sp(sp, 1) = molecule%basis2shell(first_x)
               qual_sp(sp, 2) = molecule%basis2shell(first_y)
               qual_sp(sp, 1) = n_qual_aop_in_sp
!
               call mem%dealloc_int(sorted_qual_aop_in_sp_indices, n_qual_in_sp, 1)
               call mem%dealloc(sorted_qual_aop_in_sp, n_qual_in_sp, 1)
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
!        Cut out the qualified parts of the aop and sp lists
!
         call mem%alloc_int(qual_aop_copy, n_qual_aop, 3)
         call mem%alloc_int(qual_sp_copy, n_qual_sp, 3)
!
         qual_aop_copy(:, :) = qual_aop(1 : n_qual_aop, : )
         qual_sp_copy(:, :)  = qual_sp(1 : n_qual_sp, : )
!
         call mem%dealloc_int(qual_aop, 100, 3)
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
         call mem%alloc(g_wxyz, n_sig_aop, n_qual_aop)
!
         do CD_sp = 1, n_qual_sp
!
            C                = qual_sp(sp, 1)
            D                = qual_sp(sp, 2)
            n_qual_aop_in_sp = qual_sp(sp, 3)
!
            C_interval = molecule%get_shell_limits(C)
            D_interval = molecule%get_shell_limits(D)
!
!           Calculate the ({wx} | J) integrals,
!           where {wx} is the screened list of integrals
!
            AB_sp                  = 1
            sig_AB_sp              = 1
            n_qual_aop_in_prev_sps = 0
!
            do B = 1, n_s
               do A = B, n_s
!
                  if (sig_sp(AB_sp)) then
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
                                 g_wxyz(sig_sp_to_first_sig_aop(sig_AB_sp) + wx_packed - 1, aop + n_qual_aop_in_prev_sps) &
                                       = g_ABCD(wx, yz)
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
                                 g_wxyz(sig_sp_to_first_sig_aop(sig_AB_sp) + wx - 1, aop + n_qual_aop_in_prev_sps) &
                                       = g_ABCD(wx, yz)
!
                              enddo
                           enddo
!
                        endif
!
                     enddo
!
                     call mem%dealloc(g_ABCD, &
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
!        Subtract old cholesky vectors
!
         if (n_cholesky .ne. 0) then
!
            call dgemm('N', 'T',    &
                        n_sig_aop,  &
                        n_qual_aop, &
                        n_cholesky, &
                        -one,       &
                        cholesky,   &
                        n_sig_aop,  &
                        cholesky,   &
                        n_sig_aop,  &
                        one,        &
                        g_wxyz,     &
                        n_sig_aop)
!
         endif
!
         call mem%alloc(cholesky_new, n_sig_aop, n_qual_aop)
         cholesky_new = zero
!
         current_qual = 0
!
         D_max = one
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
               if (D_xy(xy, 1) .ge. D_max) then
!
                  qual_max = qual
                  max_xy = xy
                  D_max    = D_xy(xy, 1)
!
               endif
!
            enddo
!
            if (D_max .gt. threshold) then
!
               cholesky_basis(n_cholesky + current_qual, 1) = qual_aop(qual_max, 1)
               cholesky_basis(n_cholesky + current_qual, 2) = qual_aop(qual_max, 2)
!
               cholesky_new(: , current_qual) = g_wxyz(:, qual_max)/sqrt(D_max)
!
               do xy = 1, n_sig_aop
!
                  D_xy(xy, 1) = D_xy(xy, 1) - cholesky_new(xy, current_qual)**2
!
                  do K = 1, n_qual_aop
!
                     g_wxyz(xy, K) = g_wxyz(xy, K) - cholesky_new(xy, current_qual)*cholesky_new(K, current_qual)
!
                  enddo
!
               enddo
!
            else
!
               construct_more_choleskys = .false.
!
            endif
!
         enddo
!
         n_new_cholesky = current_qual
!
!        Find new significant diagonals
!
         n_new_sig_sp = 0
         allocate(new_sig_sp(n_sig_sp,1))
         new_sig_sp = .false.
!
         do sp = 1, n_sig_sp
!
            first = sig_sp_to_first_sig_aop(sp, 1)
            last  = sig_sp_to_first_sig_aop(sp + 1, 1) - 1
!
            new_sig_sp(sp) = is_significant(D_xy(first:last, 1), &
                                            last - first + 1,    &
                                            threshold)
!
            if (new_sig_sp(sp)) n_new_sig_sp = n_new_sig_sp + 1
!
         enddo
!
!        Update index lists: sps -> aops and aops -> aos
!
         call mem%alloc_int(new_sig_sp_to_first_sig_aop, n_new_sig_sp, 1)
         new_sig_sp_to_first_sig_aop = 0
!
         new_sig_sp    = 1
         n_new_sig_aop = 0
         first_sig_aop = 1
!
         do sp = 1, n_sig_sp
!
            if (new_sig_sp(sp)) then
!
               A = sig_sp_to_shells(sp, 1)
               B = sig_sp_to_shells(sp, 2)
!
               A_interval = get_shell_limits(A)
               B_interval = get_shell_limits(B)
!
               new_sig_sp_to_first_sig_aop(new_sig_sp, 1) = first_sig_aop
!
               first_sig_aop = first_sig_aop + get_size_sp(A_interval, B_interval)
               n_new_sig_aop = first_sig_aop - 1
               new_sig_sp    = new_sig_sp + 1
            endif
!
         enddo
!
         call reduce_vector(sig_aop_to_aos,          &
                            new_sig_aop_to_aos,      &
                            sig_sp_to_first_sig_aop, &
                            new_sig_sp,              &
                            n_sig_sp,                &
                            n_sig_aop,               &
                            n_new_sig_aop)
!
         call mem%alloc(D_xy_new, n_new_sig_sp, 1)
!
         call reduce_vector(D_xy, D_xy_new, sig_sp_to_first_sig_aop, new_sig, n_sig_sp, n_sig_aop, n_new_sig_aop)
!
         call mem%dealloc(D_xy, n_sig_aop, 1)
         call mem%alloc(D_xy, n_new_sig_aop, 1)
!
         call dcopy(n_new_sig_aop, D_xy_new, 1, D_xy, 1)
!
         call mem%dealloc(D_xy_new, n_new_sig_aop)
!
         if (n_cholesky == 0) then
!
            call mem%alloc(cholesky, n_new_sig_aop, n_new_cholesky)
!
            call reduce_array(cholesky_new,             &
                               cholesky,                 &
                               sig_sp_to_first_sig_aop,  &
                               new_sig,                  &
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
!!
            call reduce_array(cholesky_new,                       &
                               cholesky_copy(1, n_cholesky + 1),  &
                               sig_sp_to_first_sig_aop,           &
                               new_sig_sp,                        &
                               n_sig_sp,                          &
                               n_sig_aop,                         &
                               n_new_sig_aop,                     &
                               n_new_cholesky)
!
            call mem%dealloc(cholesky_new, n_sig_aop, n_qual_aop)
!
         endif
!
!        Deallocate old lists & reallocate + copy over new lists
!
         deallocate(sig_sp)
         allocate(sig_sp(n_new_sig_sp,1))
         sig_sp = new_sig_sp
         deallocate(new_sig_sp)
!
         call mem%dealloc_int(sig_sp_to_first_sig_aop, n_sig_sp + 1, 1)
         call mem%alloc_int(sig_sp_to_first_sig_aop, n_new_sig_sp + 1, 1)
         sig_sp_to_first_sig_aop = new_sig_sp_to_first_sig_aop
         call mem%dealloc_int(new_sig_sp_to_first_sig_aop, n_new_sig_sp + 1, 1)
!
         call mem%dealloc_int(sig_sp_to_first_sig_aop, n_sig_sp, 1)
         call mem%alloc_int(sig_sp_to_first_sig_aop, n_new_sig_sp, 1)
         sig_sp_to_first_sig_aop = new_sig_sp_to_first_sig_aop
         call mem%dealloc_int(new_sig_sp_to_first_sig_aop, n_new_sig_sp, 1)
!
         n_sig_sp = n_new_sig_sp
!
         n_cholesky = n_cholesky + n_new_cholesky
!
         call mem%dealloc_int(qual_aop, n_qual_aop, 3)
         call mem%dealloc_int(qual_sp, n_qual_sp, 3)
!
      enddo
!
!
   end subroutine cholesky_decompose_integral_manager
!
!
   function get_size_sp(A_interval, B_interval)
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
end submodule cholesky
