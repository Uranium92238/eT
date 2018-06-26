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
      integer(i15) :: offset
!
      logical, dimension(:), allocatable        :: screened
      integer(i15), dimension(:), allocatable   :: offsets
!
      real(dp), dimension(:,:), allocatable :: D_xy
!
      integer(i15)   :: A, B, AB = 0
      type(interval) :: A_intval, B_intval
!
      integer(i15) :: xy = 0, xy_packed = 0
!
      n_shells = molecule%get_n_shells()
      n_shell_pairs = n_shells*(n_shells + 1)/2
!
!     Determine the number of diagonals and offsets for shell pairs
!
      allocate(screened(n_shell_pairs))
      allocate(offsets(n_shell_pairs))
!
      dim_D_xy = 0
      offset = 1
!
      do B = 1, n_shells
         do A = B, n_shells
!
            A_intval = molecule%get_shell_limits(A)
            B_intval = molecule%get_shell_limits(B)
!
            AB = index_packed(A, B)
            offsets(AB) = offset
!
            if (A .eq. B) then
!
               dim_D_xy = dim_D_xy + A_intval%size*(A_intval%size + 1)/2
               offset = offset + A_intval%size*(A_intval%size + 1)/2
!
            else
!
               dim_D_xy = dim_D_xy + (A_intval%size)*(B_intval%size)
               offset = offset + (A_intval%size)*(B_intval%size)
!
            endif
!
         enddo
      enddo
!
!     Construct the diagonal and screen if below a threshold
!
      call mem%alloc(D_xy, dim_D_xy, 1)
      D_xy = zero
!
      screened = .true.
!
      counter = 1
!
       do B = 1, n_shells
         do A = B, n_shells
!
            A_intval = molecule%get_shell_limits(A)
            B_intval = molecule%get_shell_limits(B)
!
            call mem%alloc(g_wxyz, (A_intval%size)*(B_intval%size), (A_intval%size)*(B_intval%size))
            call integrals%get_g_wxyz(g_wxyz, A, B, A, B)
!
            if (A == B) then
!
               do x = 1, A_intval%size
                  do y = 1, x
!
                     xy = A_intval%size*(y-1)+x
                     xy_packed = (max(i,j)*(max(i,j)-3)/2) + i + j
!
                     D_xy(offsets(counter, 1) + xy_packed - 1, 1) = g_wxyz(xy, xy)
!
                     if (D_xy(offsets(counter) + xy_packed - 1, 1) > 1.0d-8) then
!
                        screened(counter) = .false.
!
                     endif
!
                  enddo
               enddo
!
            else
!
               do x = 1, A_intval%size
                  do y = 1, B_intval%size
!
                     xy = A_intval%size*(y-1)+x
!
                     D_xy(offsets(counter, 1) + xy - 1, 1) = g_wxyz(xy, xy)
!
                     if (D_xy(offsets(counter) + xy - 1, 1) > 1.0d-8) then
!
                        screened(counter) = .false.
!
                     endif
!
                  enddo
               enddo
!
            endif
!
            call mem%dealloc(g_wxyz, (A_intval%size)*(B_intval%size), (A_intval%size)*(B_intval%size))
!
            counter = counter + 1
!
         enddo
      enddo
!
!     Determine the dimension of screened diagonal array
!
      new_dim_D_xy = 0
!
      do AB = 1, n_shell_pairs
!
         if (.not. screened(AB)) then
!
            if (AB == n_shell_pairs) then
!
               new_dim_D_xy = new_dim_D_xy + (dim_D_xy - offset(AB) + 1)
!
            else
!
               new_dim_D_xy = new_dim_D_xy + (offset(AB + 1) - offset(AB) + 1)
!
            endif
!
         endif
!
      enddo
!
!     Set the screened diagonals from the full diagonal array
!
      mem%alloc(new_D_xy, new_dim_D_xy, 1)
      new_D_xy = zero
!
      offset = 0
      do AB = 1, n_shell_pairs
!
         if (.not. screened(AB)) then
!
            if (AB .eq. n_shell_pairs) then
!
               length = dim_D_xy - offsets(AB)
               new_D_xy(offset:new_dim_D_xy) = D_xy(offsets(AB):dim_D_xy)
!
            else
!
               length = offsets(AB + 1) - offsets(AB)
               new_D_xy(offset:(offset + length - 1)) = D_xy(offsets(AB):(offsets(AB + 1) - 1))
!
            endif
!
            offset = offset + length
!
         endif
!
      enddo
!
      call mem%dealloc(D_xy, dim_D_xy, 1)
!
   end subroutine cholesky_decompose_integral_manager
!
!
end submodule cholesky
