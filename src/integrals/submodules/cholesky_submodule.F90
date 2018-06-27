submodule (integral_manager_class) cholesky
!
!!
!!    Cholesky submodule
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
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
      integer(i15) :: offset, counter
!
      logical, dimension(:), allocatable        :: screened
      integer(i15), dimension(:), allocatable   :: offsets, reduced_offsets
!
      real(dp), dimension(:,:), allocatable :: D_xy, g_wxyz
!  
      integer(i15), dimension(:,:), allocatable :: index_list
!
      integer(kind = 4)   :: A, B
      integer(i15)        :: AB = 0
      type(interval) :: A_intval, B_intval
!
      integer(i15) :: x = 0, y = 0, xy = 0, xy_packed = 0, I = 0
!
      integer(i15) :: dim_screened_diagonal, dim_screened_shell_pairs
!
      n_shells = molecule%get_n_shells()
      n_shell_pairs = n_shells*(n_shells + 1)/2
!
!     Determine the number of screened diagonals without storing
!
      allocate(screened(n_shell_pairs))
      allocate(offsets(n_shell_pairs))
!
      dim_screened_shell_pairs = 0
      dim_screened_diagonal    = 0
!
      offsets = 1
!
      do B = 1, n_shells
         do A = B, n_shells
!
            A_intval = molecule%get_shell_limits(A)
            B_intval = molecule%get_shell_limits(B)
!
            call mem%alloc(g_wxyz, (A_intval%size)*(B_intval%size), (A_intval%size)*(B_intval%size))
            call integrals%get_ao_g_wxyz(g_wxyz, A, B, A, B)
!
            AB = index_packed(A, B)
            offsets(AB) = offset
!
!           Determine whether shell pair is screened out
!
            screened(AB) = .true.
!
            do I = 1, (A_intval%size)*(B_intval%size)
!
               if (g_wxyz(I,I) .gt. 1.0D-8) then
!
                  screened(AB) = .false.
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
            if (.not. screened(AB)) then
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
         enddo
      enddo
!
      call mem%alloc(D_xy, dim_screened_diagonal, 1)
      D_xy = zero
!
      call mem%alloc_int(index_list, dim_screened_diagonal, 2)
      index_list = 0
!
!     Construct the screened diagonal
!
      offset = 1
!
      do B = 1, n_shells
         do A = B, n_shells
!
            A_intval = molecule%get_shell_limits(A)
            B_intval = molecule%get_shell_limits(B)
!
            AB = index_packed(A, B)
!
            if (.not. screened(AB)) then
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
            endif
!
         enddo
      enddo
!
!     Make reduced offset list for the diagonals that are kept
!
      call mem%alloc(reduced_offsets, dim_screened_shell_pairs, 1)
!
      counter = 1
!
      do AB = 1, n_shell_pairs
!
         if (.not. screened(AB)) then
!
            reduced_offsets(counter, 1) = offsets(AB, 1)
            counter = counter + 1
!
         endif
!
      enddo



!     Old code below
!
!     Determine the number of diagonals and offsets for shell pairs
!
!
!       dim_D_xy = 0
!       offset = 1
! !
!       do B = 1, n_shells
!          do A = B, n_shells
! !
!             A_intval = molecule%get_shell_limits(A)
!             B_intval = molecule%get_shell_limits(B)
! !
!             AB = index_packed(A, B)
!             offsets(AB) = offset
! !
!             if (A .eq. B) then
! !
!                dim_D_xy = dim_D_xy + A_intval%size*(A_intval%size + 1)/2
!                offset = offset + A_intval%size*(A_intval%size + 1)/2
! !
!             else
! !
!                dim_D_xy = dim_D_xy + (A_intval%size)*(B_intval%size)
!                offset = offset + (A_intval%size)*(B_intval%size)
! !
!             endif
! !
!          enddo
!       enddo
! !
! !     Construct the diagonal and screen if below a threshold
! !
!       call mem%alloc(D_xy, dim_D_xy, 1)
!       D_xy = zero
! !
!       screened = .true.
! !
!       counter = 1
! !
!        do B = 1, n_shells
!          do A = B, n_shells
! !
!             A_intval = molecule%get_shell_limits(A)
!             B_intval = molecule%get_shell_limits(B)
! !
!             call mem%alloc(g_wxyz, (A_intval%size)*(B_intval%size), (A_intval%size)*(B_intval%size))
!             call integrals%get_g_wxyz(g_wxyz, A, B, A, B)
! !
!             if (A == B) then
! !
!                do x = 1, A_intval%size
!                   do y = 1, x
! !
!                      xy = A_intval%size*(y-1)+x
!                      xy_packed = (max(i,j)*(max(i,j)-3)/2) + i + j
! !
!                      D_xy(offsets(counter, 1) + xy_packed - 1, 1) = g_wxyz(xy, xy)
! !
!                      if (D_xy(offsets(counter) + xy_packed - 1, 1) > 1.0d-8) then
! !
!                         screened(counter) = .false.
! !
!                      endif
! !
!                   enddo
!                enddo
! !
!             else
! !
!                do x = 1, A_intval%size
!                   do y = 1, B_intval%size
! !
!                      xy = A_intval%size*(y-1)+x
! !
!                      D_xy(offsets(counter, 1) + xy - 1, 1) = g_wxyz(xy, xy)
! !
!                      if (D_xy(offsets(counter) + xy - 1, 1) > 1.0d-8) then
! !
!                         screened(counter) = .false.
! !
!                      endif
! !
!                   enddo
!                enddo
! !
!             endif
! !
!             call mem%dealloc(g_wxyz, (A_intval%size)*(B_intval%size), (A_intval%size)*(B_intval%size))
! !
!             counter = counter + 1
! !
!          enddo
!       enddo
! !
! !     Determine the dimension of screened diagonal array
! !
!       new_dim_D_xy = 0
! !
!       do AB = 1, n_shell_pairs
! !
!          if (.not. screened(AB)) then
! !
!             if (AB == n_shell_pairs) then
! !
!                new_dim_D_xy = new_dim_D_xy + (dim_D_xy - offset(AB) + 1)
! !
!             else
! !
!                new_dim_D_xy = new_dim_D_xy + (offset(AB + 1) - offset(AB) + 1)
! !
!             endif
! !
!          endif
! !
!       enddo
! !
! !     Set the screened diagonals from the full diagonal array
! !
!       mem%alloc(new_D_xy, new_dim_D_xy, 1)
!       new_D_xy = zero
! !
!       offset = 0
!       do AB = 1, n_shell_pairs
! !
!          if (.not. screened(AB)) then
! !
!             if (AB .eq. n_shell_pairs) then
! !
!                length = dim_D_xy - offsets(AB)
!                new_D_xy(offset:new_dim_D_xy) = D_xy(offsets(AB):dim_D_xy)
! !
!             else
! !
!                length = offsets(AB + 1) - offsets(AB)
!                new_D_xy(offset:(offset + length - 1)) = D_xy(offsets(AB):(offsets(AB + 1) - 1))
! !
!             endif
! !
!             offset = offset + length
! !
!          endif
! !
!       enddo
! !
!       call mem%dealloc(D_xy, dim_D_xy, 1)
!
   end subroutine cholesky_decompose_integral_manager
!
!
end submodule cholesky
