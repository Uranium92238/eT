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
      implicit none
!
      class(integral_manager) :: integrals
      class(molecular_system) :: molecule
!
      integer(i15) :: n_s, n_sp
!
      integer(i15) :: A, B
!
      integer(i15) :: n_significant_aop
      integer(i15) :: n_significant_sp
!
      logical, dimension(:,:), allocatable :: significant
!
      real(dp), parameter :: threshold = 1.0D-8
!
      n_s   = molecule%get_n_shells() ! Number of shells
      n_sp  = n_s*(n_s + 1)/2         ! Number of shell pairs packed
!
!     Pre-screening of full diagonal
!
      allocate(significant(n_sp, 1))
      significant = .false.
!
      sp = 1                ! Shell pair number
      n_significant_aop = 0 ! Number of significant AO pairs
      n_significant_sp  = 0 ! Number of significant shell pairs
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
            significant(sp) = is_significant(D_AB, (A_interval%size)*(B_interval%size), threshold))
!
            call mem%dealloc(D_AB, (A_interval%size)*(B_interval%size), 1)
!
            if (significant(sp)) then
!
               n_significant_aop = n_significant_aop + &
                              get_size_sp(A_interval, B_interval)
!
               n_significant_sp = n_significant_sp + 1
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
      call mem%alloc_int(significant_sp_to_first_significant_aop(n_significant_sp))
      significant_sp_to_first_significant_aop = 0
      first_significant_aop = 1
!
      sp = 1
      significant_sp = 1
!
      do B = 1, n_s
         do A = B, n_s
!
            if (significant(sp)) then
!
               significant_sp_to_first_significant_aop(significant_sp, 1) = first_significant_aop
!
               A_interval = molecule%get_shell_limits(A)
               B_interval = molecule%get_shell_limits(B)
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
                        D_xy(xy_packed + first_significant_aop - 1, 1) = g_wxyz(xy, xy)
!
                        significant_aop_to_aos(xy_packed + first_significant_aop - 1, 1) = &
                                                                           A_interval%first + x - 1
!
                        significant_aop_to_aos(xy_packed + first_significant_aop - 1, 2) = &
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
                        D_xy(xy + first_significant_aop - 1, 1) = g_wxyz(xy,xy)
!
                        significant_aop_to_aos(xy + first_significant_aop - 1, 1) = A_interval%first + x - 1
                        significant_aop_to_aos(xy + first_significant_aop - 1, 2) = B_interval%first + y - 1
!
                     enddo
                  enddo
!
               endif
!
               first_significant_aop = first_significant_aop + get_size_sp(A_interval, B_interval)
!
               significant_sp = significant_sp + 1
!
            endif ! End of if (significant)
!
            sp = sp + 1
!
         enddo
      enddo
!
!
   end subroutine cholesky_decompose_integral_manager
!
      integer(i15) :: n_s, n_sp
      integer(i15) :: dim_screened, n_screened_sp
      integer(i15) :: first, last
      integer(i15) :: sp, screened_sp, ab_sp, ao, cd_sp, screened_ab_sp
      integer(i15) :: aop, current_screened_sp
      integer(i15) :: A, B, I, C, D, J, K
      integer(i15) :: x, y, xy, xy_packed, w, wx, wx_packed, z, yz
      integer(i15) :: offset
      integer(i15) :: first_screened_aop, last_screened_aop, first_x, first_y
      integer(i15) :: max_qualified, n_qualified_in_sp, n_qualified, n_old_qualified, n_qualified_sp
      integer(i15) :: current_qual, qual_offset, qual, qual_max
!
      real(dp) :: diag_max, span
!
      type(interval) :: A_interval, B_interval, C_interval, D_interval
!
      logical, dimension(:), allocatable        :: screened
!
      real(dp), dimension(:), allocatable       :: sp_offsets
!
      real(dp), dimension(:,:), allocatable     :: g_wxyz, g_ABCD
      real(dp), dimension(:,:), allocatable     :: diag_xy
      real(dp), dimension(:,:), allocatable     :: max_in_sp
      real(dp), dimension(:,:), allocatable     :: sorted_max_in_sp
      real(dp), dimension(:,:), allocatable     :: sorted_qualified_in_sp
      real(dp), dimension(:,:), allocatable     :: cholesky
!
      integer(i15), dimension(:), allocatable   :: max_in_sp_indices, screened_sp_offsets
!
      integer(i15), dimension(:,:), allocatable :: diag_to_aos
      integer(i15), dimension(:,:), allocatable :: sorted_max_sp
      integer(i15), dimension(:,:), allocatable :: qualified_aop
      integer(i15), dimension(:,:), allocatable :: qualified_aop_copy
      integer(i15), dimension(:,:), allocatable :: qualified_sp
      integer(i15), dimension(:,:), allocatable :: qualified_sp_copy
      integer(i15), dimension(:,:), allocatable :: sorted_qualified_in_sp_indices
!
      n_s   = molecule%get_n_shells()     ! number of shells
      n_sp  = n_s*(n_s + 1)/2             ! number of shell pairs packed
!
      write(output%unit, *) 'Number of shells: ', n_s
      write(output%unit, *) 'Number of shell pairs: ', n_sp
!
!     Determine the number of screened diagonals without storing
!
      allocate(screened(n_sp))
      allocate(sp_offsets(n_sp))
!
      n_screened_sp  = 0
      dim_screened   = 0
!
      sp_offsets   = zero
      screened     = .true.

      sp = 1
      offset = 1
!
      do B = 1, n_s
         do A = B, n_s
!
            A_interval = molecule%get_shell_limits(A)
            B_interval = molecule%get_shell_limits(B)
!
            call mem%alloc(g_wxyz, &
                     (A_interval%size)*(B_interval%size), &
                     (A_interval%size)*(B_interval%size))
!
            g_wxyz = zero
!
            call integrals%get_ao_g_wxyz(g_wxyz, A, B, A, B)
!
            sp_offsets(sp) = offset
!
!           Determine whether shell pair is screened out
!
            screened(sp) = .true.
!
            do I = 1, (A_interval%size)*(B_interval%size)
!
               if (g_wxyz(I,I) .gt. 1.0D-8) then
!
                  screened(sp) = .false.
!
               endif
!
            enddo
!
            call mem%dealloc(g_wxyz, &
                     (A_interval%size)*(B_interval%size), &
                     (A_interval%size)*(B_interval%size))
!
!           Update AO pair offsets for shell pair
!
            if (A .eq. B) then
!
               offset = offset + A_interval%size*(A_interval%size + 1)/2
!
            else
!
               offset = offset + (A_interval%size)*(B_interval%size)
!
            endif
!
!           Add contribution to screened diagonal dimension, if not screened,
!           and set +1 to the dimension of surviving shell pairs
!
            if (.not. screened(sp)) then
!
               n_screened_sp = n_screened_sp + 1
!
               if (A .eq. B) then
!
                  dim_screened = dim_screened + A_interval%size*(A_interval%size + 1)/2
!
               else
!
                  dim_screened = dim_screened + (A_interval%size)*(B_interval%size)
!
               endif
!
            endif
!
            sp = sp + 1
!
         enddo
      enddo
!
      write(output%unit, *) 'Dim of screened diagonal: ', dim_screened
      flush(output%unit)
!
      call mem%alloc(diag_xy, dim_screened, 1)
      diag_xy = zero
!
      call mem%alloc_int(diag_to_aos, dim_screened, 2)
      diag_to_aos = 0
!
!     Construct the screened diagonal and simultaneously get new offsets
!
      offset = 1
      sp = 1
      screened_sp = 1
!
      allocate(screened_sp_offsets(n_screened_sp))
!
      do B = 1, n_s
         do A = B, n_s
!
            write(output%unit, *) '1'
            flush(output%unit)
!
            A_interval = molecule%get_shell_limits(A)
            B_interval = molecule%get_shell_limits(B)
!
            if (.not. screened(sp)) then
!
               screened_sp_offsets(screened_sp) = offset
!
               call mem%alloc(g_wxyz, (A_interval%size)*(B_interval%size), (A_interval%size)*(B_interval%size))
               call integrals%get_ao_g_wxyz(g_wxyz, A, B, A, B)
!
               if (A .eq. B) then
!
                  do x = 1, A_interval%size
                     do y = 1, B_interval%size
!
                        xy = A_interval%size*(y - 1) + x
                        xy_packed = (max(x,y)*(max(x,y)-3)/2) + x + y
!
                        diag_xy(offset + xy_packed - 1, 1) = g_wxyz(xy, xy)
!
                        diag_to_aos(offset + xy_packed - 1, 1) = A_interval%first + x - 1
                        diag_to_aos(offset + xy_packed - 1, 2) = B_interval%first + y - 1
!
                     enddo
                  enddo
!
                  offset = offset + A_interval%size*(A_interval%size + 1)/2
!
               else
!
                  do x = 1, (A_interval%size)
                     do y = 1, (B_interval%size)
!
                        xy = A_interval%size*(y - 1) + x
                        diag_xy(offset + xy - 1, 1) = g_wxyz(xy,xy)
!
                        diag_to_aos(offset + xy - 1, 1) = A_interval%first + x - 1
                        diag_to_aos(offset + xy - 1, 2) = B_interval%first + y - 1
!
                     enddo
                  enddo
!
                  offset = offset + (A_interval%size)*(B_interval%size)
!
               endif
!
               call mem%dealloc(g_wxyz, (A_interval%size)*(B_interval%size), (A_interval%size)*(B_interval%size))
!
               screened_sp = screened_sp + 1
!
            endif
!
            sp = sp + 1
!
         enddo
      enddo
!
!     Shell maximums and shell maximums indices vectors
!
      allocate(max_in_sp_indices(n_screened_sp))
      call mem%alloc(max_in_sp, n_screened_sp, 1)
!
      max_in_sp = zero
!
      do sp = 1, n_screened_sp
!
!        Get first and last indices of shell pair
!
         first = screened_sp_offsets(sp)
!
         if (sp .eq. n_screened_sp) then
!
            last = dim_screened
!
         else
!
            last = screened_sp_offsets(sp + 1) - 1
!
         endif
!
!        Determine the largest elements
!
         do I = first, last
!
            if (diag_xy(I, 1) .gt. max_in_sp(sp, 1)) then
!
               max_in_sp(sp, 1) = diag_xy(I, 1)
               max_in_sp_indices(sp) = I
!
            endif
!
         enddo
!
      enddo
!
!     Sort from largest to smallest by determining an index array
!
      allocate(sorted_max_sp(n_screened_sp,1))
      sorted_max_sp = 0
!
      call mem%alloc(sorted_max_in_sp, n_screened_sp, 1)
      sorted_max_in_sp = zero
!
      call get_n_highest(n_screened_sp, n_screened_sp, &
                        max_in_sp, sorted_max_in_sp, sorted_max_sp)
!
      diag_max       = sorted_max_in_sp(1, 1)
      n_qualified    = 0
      n_qualified_sp = 0
      span           = 1.0D-3
      max_qualified  = 100
!
      call mem%alloc_int(qualified_aop, 100, 3)
      call mem%alloc_int(qualified_sp, n_sp, 3)
!
      do sp = 1, n_screened_sp
!
         current_screened_sp = sorted_max_sp(sp, 1)
!
         first_screened_aop = screened_sp_offsets(current_screened_sp)
!
         if (current_screened_sp .eq. n_screened_sp) then
!
            last_screened_aop = dim_screened
!
         else
!
            last_screened_aop = screened_sp_offsets(current_screened_sp + 1) - 1
!
         endif
!
         n_qualified_in_sp = 0
!
         do aop = first_screened_aop, last_screened_aop
!
            if ((diag_xy(aop, 1) .ge. span*diag_max ) .and. (n_qualified .lt. max_qualified)) then
!
               n_qualified_in_sp = n_qualified_in_sp + 1
               n_qualified       = n_qualified + 1
!
            endif
!
         enddo
!
         if (n_qualified_in_sp .ne. 0) then
!
            n_qualified_sp = n_qualified_sp + 1
!
            call mem%alloc_int(sorted_qualified_in_sp_indices, n_qualified_in_sp, 1)
            call mem%alloc(sorted_qualified_in_sp, n_qualified_in_sp, 1)
!
            call get_n_highest(n_qualified_in_sp, last_screened_aop - first_screened_aop + 1, &
                              diag_xy(first_screened_aop:last_screened_aop, 1), sorted_qualified_in_sp, &
                              sorted_qualified_in_sp_indices)
!
            n_old_qualified = (n_qualified - n_qualified_in_sp)
!
            do aop = 1, n_qualified_in_sp
!
               qualified_aop(aop + n_old_qualified, 1) = diag_to_aos(sorted_qualified_in_sp_indices(aop, 1) &
                                                            + first_screened_aop - 1, 1)
!
               qualified_aop(aop + n_old_qualified, 2) = diag_to_aos(sorted_qualified_in_sp_indices(aop, 1) &
                                                            + first_screened_aop - 1, 2)
!
               qualified_aop(aop + n_old_qualified, 3) = sorted_qualified_in_sp_indices(aop, 1) &
                                                            + first_screened_aop - 1
!
            enddo
!
            first_x = diag_to_aos(first_screened_aop, 1) ! alpha
            first_y = diag_to_aos(first_screened_aop, 2) ! beta
!
            qualified_sp(sp, 1) = molecule%basis2shell(first_x)
            qualified_sp(sp, 2) = molecule%basis2shell(first_y)
            qualified_sp(sp, 3) = n_qualified_in_sp
!
            call mem%dealloc_int(sorted_qualified_in_sp_indices, n_qualified_in_sp, 1)
            call mem%dealloc(sorted_qualified_in_sp, n_qualified_in_sp, 1)
!
         endif
!
         if (n_qualified == max_qualified) then
!
            exit
!
         endif
!
      enddo
!
!     Cut out the qualified parts of the aop and sp lists
!
      call mem%alloc_int(qualified_aop_copy, n_qualified, 3)
      call mem%alloc_int(qualified_sp_copy, n_qualified_sp, 3)
!
      qualified_aop_copy(:, :) = qualified_aop(1:n_qualified, :)
      qualified_sp_copy(:, :)  = qualified_sp(1:n_qualified_sp, :)
!
      call mem%dealloc_int(qualified_aop, 100, 3)
      call mem%dealloc_int(qualified_sp, n_sp, 3)
!
      call mem%alloc_int(qualified_aop, n_qualified, 3)
      call mem%alloc_int(qualified_sp, n_qualified_sp, 3)
!
      qualified_aop = qualified_aop_copy
      qualified_sp = qualified_sp_copy
!
      call mem%dealloc_int(qualified_aop_copy, n_qualified, 3)
      call mem%dealloc_int(qualified_sp_copy, n_qualified_sp, 3)
!
      qual_offset = 0
      offset = 0
!
      call mem%alloc(g_wxyz, &
                        dim_screened, &
                        n_qualified)
!
      do cd_sp = 1, n_qualified_sp
!
         C                 = qualified_sp(sp, 1)
         D                 = qualified_sp(sp, 2)
         n_qualified_in_sp = qualified_sp(sp, 3)
!
         C_interval = molecule%get_shell_limits(C)
         D_interval = molecule%get_shell_limits(D)
!
!        Calculate the ({wx} | J) integrals,
!        where {wx} is the screened list of integrals
!
         ab_sp = 1
         screened_ab_sp = 1
         offset = 0
!
         do B = 1, n_s
            do A = B, n_s
!
               if (.not. screened(ab_sp)) then
!
                  A_interval = molecule%get_shell_limits(A)
                  B_interval = molecule%get_shell_limits(B)
!
                  call mem%alloc(g_ABCD, &
                                 (A_interval%size)*(B_interval%size), &
                                 (C_interval%size)*(D_interval%size))
!
                  call integrals%get_ao_g_wxyz(g_ABCD, A, B, C, D)
!
                  do ao = 1, n_qualified_in_sp
!
                     y = qualified_aop(ao + offset, 1)
                     z = qualified_aop(ao + offset, 2)
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
                              g_wxyz(screened_sp_offsets(screened_ab_sp) + wx_packed - 1, ao + qual_offset) = g_ABCD(wx, yz)
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
                              g_wxyz(screened_sp_offsets(screened_ab_sp) + wx - 1, ao + qual_offset) = g_ABCD(wx, yz)
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
                  screened_ab_sp = screened_ab_sp + 1
!
               endif
!
               ab_sp = ab_sp + 1
!
            enddo ! A
         enddo ! B
!
         qual_offset = qual_offset + n_qualified_in_sp
         offset = offset + n_qualified_in_sp
!
      enddo ! cd_sp
!
      call mem%alloc(cholesky, dim_screened, n_qualified )
      cholesky = zero
!
      current_qual = 1
!
      diag_max = one
!
      do while ((current_qual .lt. n_qualified) .and. (diag_max .gt. 1.0d-8))
!
         diag_max = zero
!
         do qual = 1, n_qualified
!
            xy = qualified_aop(qual, 3)
!
            if (diag_xy(xy, 1) .ge. diag_max) then
!
               qual_max = qual
               diag_max = diag_xy(xy, 1)
!
            endif
!
         enddo
!
         cholesky(: , current_qual) = g_wxyz (:, qual_max)/sqrt(diag_max)
!
         do xy = 1, dim_screened
!
            diag_xy(xy, 1) = diag_xy(xy, 1) - cholesky(xy, current_qual)**2
!
            do K = 1, n_qualified
!
               g_wxyz(xy, K) = g_wxyz(xy, K) - cholesky(xy, current_qual)*cholesky(K, current_qual)
!
            enddo
!
         enddo
!
         current_qual = current_qual + 1
!
      enddo
!
      call mem%dealloc(cholesky, dim_screened, n_qualified)
!
      deallocate(screened_sp_offsets)
      call mem%dealloc_int(qualified_aop, n_qualified, 3)
      call mem%dealloc_int(qualified_sp, n_qualified_sp, 3)
!
   end subroutine cholesky_decompose_integral_manager
!
!
   module function get_size_sp(A_interval, B_interval)
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
   end function get_size_AB
!
!
end submodule cholesky
