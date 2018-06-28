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
      integer(i15) :: dim_screened, n_screened_sp
      integer(i15) :: first, last
      integer(i15) :: sp, screened_sp
      integer(i15) :: aop, current_screened_sp
      integer(i15) :: A, B, I
      integer(i15) :: x, y, xy, xy_packed
      integer(i15) :: offset
      integer(i15) :: first_screened_aop, last_screened_aop, first_x, first_y
      integer(i15) :: max_qualified, n_qualified_in_sp, n_qualified, n_old_qualified
!
      real(dp) :: diag_max, span
!
      type(interval) :: A_interval, B_interval
!
      logical, dimension(:), allocatable        :: screened
!
      real(dp), dimension(:), allocatable       :: sp_offsets, screened_sp_offsets
!
      real(dp), dimension(:,:), allocatable     :: g_wxyz
      real(dp), dimension(:,:), allocatable     :: diag_xy
      real(dp), dimension(:,:), allocatable     :: max_in_sp
      real(dp), dimension(:,:), allocatable     :: sorted_max_in_sp
      real(dp), dimension(:,:), allocatable     :: sorted_qualified_in_sp
!
      integer(i15), dimension(:), allocatable   :: max_in_sp_indices
!
      integer(i15), dimension(:,:), allocatable :: diag_to_aos
      integer(i15), dimension(:,:), allocatable :: sorted_max_sp
      integer(i15), dimension(:,:), allocatable :: qualified_aop
      integer(i15), dimension(:,:), allocatable :: qualified_sp
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
                        diag_to_aos(offset + xy_packed - 1, 1) = B_interval%first + y - 1
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
                        diag_to_aos(offset + xy - 1, 1) = B_interval%first + y - 1
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
      call mem%alloc_int(qualified_aop, 100, 2)
      call mem%alloc_int(qualified_sp, n_s, 2)
!
      do sp = 1, n_screened_sp
!
         current_screened_sp = sorted_max_sp(sp, 1)
!
         first_screened_aop = screened_sp_offsets(current_screened_sp)
!
         if (sp == n_screened_sp) then
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
            if ((diag_xy(aop, 1) .le. span*diag_max ) .and. (n_qualified .lt. max_qualified)) then
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
                              diag_xy(first:last, 1), sorted_qualified_in_sp, &
                              sorted_qualified_in_sp_indices)
!
            n_old_qualified = (n_qualified - n_qualified_in_sp)
!
            do aop = 1, n_qualified_in_sp
!
               qualified_aop(aop + n_old_qualified, 1) = diag_to_aos(sorted_qualified_in_sp_indices(aop, 1) &
                                                            + first_screened_aop - 1, 1)
               qualified_aop(aop + n_old_qualified, 2) = diag_to_aos(sorted_qualified_in_sp_indices(aop, 1) &
                                                            + first_screened_aop - 1, 2)
!
            enddo
!
            first_x = diag_to_aos(first_screened_aop, 1) ! alpha
            first_y = diag_to_aos(first_screened_aop, 2) ! beta
!
            qualified_sp(sp, 1) = molecule%basis2shell(first_x)
            qualified_sp(sp, 2) = molecule%basis2shell(first_y)
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
      call mem%alloc_int(qualified_aop_copy, n_qualified, 2)
      call mem%alloc_int(qualified_sp_copy, n_qualified_sp, 2)
!
      qualified_aop_copy(:, :) = qualified_aop(1:n_qualified, :)
      qualified_sp_copy(:, :)  = qualified_sp(1:n_qualified_sp, :)
!
      call mem%dealloc_int(qualified_aop, 100, 2)
      call mem%dealloc_int(qualified_sp, n_s, 2)
!
      call mem%alloc_int(qualified_aop, n_qualified, 2)
      call mem%alloc_int(qualified_sp, n_qualified_sp, 2)
!
      qualified_aop = qualified_aop_copy
      qualified_sp = qualified_sp_copy
!
      call mem%dealloc_int(qualified_aop_copy, n_qualified, 2)
      call mem%dealloc_int(qualified_sp_copy, n_qualified_sp, 2)
!




!
      J = 1
!
      do cd_sp = 1, n_qualified_sp
!
         C = qualified_sp(sp, 1)
         D = qualified_sp(sp, 2)
!
         C_interval = get_shell_limits(C)
         D_interval = get_shell_limits(D)
!
!        Calculate the ({wx} | J) integrals,
!        where {wx} is the screened list of integrals
!
         call mem%alloc(g_wxyz, &
                        dim_screened_diagonal, &
                        n_qualified)
!
         ab_sp = 1
         screened_ab_sp = 1
!
         do B = 1, n_shells
            do A = B, n_shells
!
               if (.not. screened(ab_sp)) then
!
                  A_interval = get_shell_limits(A)
                  B_interval = get_shell_limits(B)
!
                  call mem%alloc(g_ABCD, &
                                 (A_interval%size)*(B_interval%size),
                                 (C_interval%size)*(D_interval%size))
!
                  call integrals%get_ao_g_wxyz(g_ABCD, A, B, C, D)
!

!
                  call mem%dealloc(g_ABCD, &
                                    (A_interval%size)*(B_interval%size),
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
      enddo ! cd_sp
!
!
      deallocate(screened_sp_offsets)
!
   end subroutine cholesky_decompose_integral_manager
!
!
end submodule cholesky
