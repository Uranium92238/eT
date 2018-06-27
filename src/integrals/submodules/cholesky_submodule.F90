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
      integer(i15) :: n_shells
      integer(i15) :: n_shell_pairs
!
      integer(i15) :: offset, counter, counter_reduced
!
      logical, dimension(:), allocatable        :: screened
      integer(i15), dimension(:), allocatable   :: offsets, reduced_offsets
!
      real(dp), dimension(:,:), allocatable :: D_xy, g_wxyz
!
      integer(i15), dimension(:,:), allocatable :: index_list
!
      integer(i15)   :: A, B
      integer(i15)   :: AB = 0, shell_pair = 0
      type(interval) :: A_intval, B_intval
!
      integer(i15) :: x = 0, y = 0, xy = 0, xy_packed = 0, I = 0
!
      integer(i15) :: max_index = 0, first = 0, last = 0
!
      integer(i15) :: dim_screened_diagonal, dim_screened_shell_pairs
!
      integer(i15), dimension(:), allocatable :: shell_max_indices
      integer(i15), dimension(:,:), allocatable :: sorted_shell_max_indices
      real(dp), dimension(:,:), allocatable :: shell_max, sorted_shell_max
!
      n_shells = molecule%get_n_shells()
      n_shell_pairs = n_shells*(n_shells + 1)/2
!
      write(output%unit, *) 'Number of shells: ', n_shells
      write(output%unit, *) 'Number of shell pairs: ', n_shell_pairs
!
!     Determine the number of screened diagonals without storing
!
      allocate(screened(n_shell_pairs))
      allocate(offsets(n_shell_pairs))
!
      dim_screened_shell_pairs = 0
      dim_screened_diagonal    = 0
      offsets = 0
!
      screened = .true.
      counter = 1
      offset = 1
!
      do B = 1, n_shells
         do A = B, n_shells
!
            A_intval = molecule%get_shell_limits(A)
            B_intval = molecule%get_shell_limits(B)
!
            call mem%alloc(g_wxyz, (A_intval%size)*(B_intval%size), (A_intval%size)*(B_intval%size))
            g_wxyz = zero
!
            call integrals%get_ao_g_wxyz(g_wxyz, A, B, A, B)
!
            offsets(counter) = offset
!
!           Determine whether shell pair is screened out
!
            screened(counter) = .true.
!
            do I = 1, (A_intval%size)*(B_intval%size)
!
               if (g_wxyz(I,I) .gt. 1.0D-8) then
!
                  screened(counter) = .false.
!
               endif
!
            enddo
!
            call mem%dealloc(g_wxyz, (A_intval%size)*(B_intval%size), (A_intval%size)*(B_intval%size))
!
!           Update AO pair offsets for shell pair
!
            if (A .eq. B) then
!
               offset = offset + A_intval%size*(A_intval%size + 1)/2
!
            else
!
               offset = offset + (A_intval%size)*(B_intval%size)
!
            endif
!
!           Add contribution to screened diagonal dimension, if not screened,
!           and set +1 to the dimension of surviving shell pairs
!
            if (.not. screened(counter)) then
!
               dim_screened_shell_pairs = dim_screened_shell_pairs + 1
!
               if (A .eq. B) then
!
                  dim_screened_diagonal = dim_screened_diagonal + &
                              A_intval%size*(A_intval%size + 1)/2
!
               else
!
                  dim_screened_diagonal = dim_screened_diagonal + &
                                 (A_intval%size)*(B_intval%size)
!
               endif
!
            endif
!
            counter = counter + 1
!
         enddo
      enddo
!
      write(output%unit, *) 'Dim of screened diagonal: ', dim_screened_diagonal
      flush(output%unit)
!
      call mem%alloc(D_xy, dim_screened_diagonal, 1)
      D_xy = zero
!
      call mem%alloc_int(index_list, dim_screened_diagonal, 2)
      index_list = 0
!
!     Construct the screened diagonal and simultaneously get new offsets
!
      offset = 1
      counter = 1
      counter_reduced = 1
!
      allocate(reduced_offsets(dim_screened_shell_pairs))
!
      do B = 1, n_shells
         do A = B, n_shells
!
            A_intval = molecule%get_shell_limits(A)
            B_intval = molecule%get_shell_limits(B)
!
            if (.not. screened(counter)) then
!
               reduced_offsets(counter_reduced) = offset
!
               call mem%alloc(g_wxyz, (A_intval%size)*(B_intval%size), (A_intval%size)*(B_intval%size))
               call integrals%get_ao_g_wxyz(g_wxyz, A, B, A, B)
!
               if (A .eq. B) then
!
                  do x = 1, A_intval%size
                     do y = 1, B_intval%size
!
                        xy = A_intval%size*(y - 1) + x
                        xy_packed = (max(x,y)*(max(x,y)-3)/2) + x + y
!
                        D_xy(offset + xy_packed - 1, 1) = g_wxyz(xy, xy)
!
                        index_list(offset + xy_packed - 1, 1) = A_intval%first + x - 1
                        index_list(offset + xy_packed - 1, 1) = B_intval%first + y - 1
!
                     enddo
                  enddo
!
                  offset = offset + A_intval%size*(A_intval%size + 1)/2
!
               else
!
                  do x = 1, (A_intval%size)
                     do y = 1, (B_intval%size)
!
                        xy = A_intval%size*(y - 1) + x
                        D_xy(offset + xy - 1, 1) = g_wxyz(xy,xy)
!
                        index_list(offset + xy - 1, 1) = A_intval%first + x - 1
                        index_list(offset + xy - 1, 1) = B_intval%first + y - 1
!
                     enddo
                  enddo
!
                  offset = offset + (A_intval%size)*(B_intval%size)
!
               endif
!
               call mem%dealloc(g_wxyz, (A_intval%size)*(B_intval%size), (A_intval%size)*(B_intval%size))
!
               counter_reduced = counter_reduced + 1
!
            endif
!
            counter = counter + 1
!
         enddo
      enddo
!
!     Find the index of the maximum value of D_xy
!
      max_index = get_max_index(D_xy, dim_screened_diagonal)
!
!     Shell maximums and shell maximums indices vectors
!
      allocate(shell_max_indices(dim_screened_shell_pairs))
      call mem%alloc(shell_max, dim_screened_shell_pairs, 1)
!
      shell_max = zero
!
      do shell_pair = 1, dim_screened_shell_pairs
!
!        Get first and last indices of shell pair
!
         first = reduced_offsets(shell_pair)
!
         if (shell_pair .eq. dim_screened_shell_pairs) then
!
            last = dim_screened_diagonal
!
         else
!
            last = reduced_offsets(shell_pair + 1) - 1
!
         endif
!
!        Determine the largest elements
!
         do I = first, last
!
            if (D_xy(I, 1) .gt. shell_max(shell_pair, 1)) then
!
               shell_max(shell_pair, 1) = D_xy(I, 1)
               shell_max_indices(shell_pair) = I
!
            endif
!
         enddo
!
      enddo
!
!     Sort from largest to smallest by determining an index array
!
      allocate(sorted_shell_max_indices(dim_screened_shell_pairs,1))
      sorted_shell_max_indices = 0
!
      call mem%alloc(sorted_shell_max, dim_screened_shell_pairs, 1)
      sorted_shell_max = zero
!
      call get_n_lowest(dim_screened_shell_pairs, dim_screened_shell_pairs, &
                        shell_max, sorted_shell_max, sorted_shell_max_indices)
!
      do shell_pair = 1, dim_screened_shell_pairs
!
         write(output%unit, *) 'The ', shell_pair, ' smallest shell is ', sorted_shell_max_indices(shell_pair, 1), &
                                 'with the value ', shell_max(sorted_shell_max_indices(shell_pair, 1), 1)
!
      enddo
!
      deallocate(reduced_offsets)
!
   end subroutine cholesky_decompose_integral_manager
!
!
end submodule cholesky
