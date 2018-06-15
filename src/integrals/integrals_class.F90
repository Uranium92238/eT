module integrals_class
!
!!
!!    Integrals class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   use h_xy    ! For one electron integrals, h_xy
   use s_xy    ! For one electron overlaps, s_xy
   use g_xyzw  ! For two electron integrals, g_xyzw
!
   use index
   use reordering
!
   use file_class
   use memory_manager_class
!
   implicit none
!
   type :: integrals
!
   contains
!
      procedure :: get_ao_xy => get_ao_xy_integrals
!
   end type integrals
!
contains
!
   subroutine get_ao_xy_integrals(int, n_ao)
!!
!!    Get AO XY integrals
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(integrals) :: int
!
      integer(i15) :: n_ao
!
      real(kind=8), dimension(:,:), allocatable :: h, s, g, g_reord, c, d, f, u,x, gd
!
      real(kind=8), dimension(:,:), allocatable :: temp, L_inv, rho_i, rho_j, s_copy, c_occ
!
      real(kind=8), dimension(:,:), allocatable :: ci, cj, eig_re, eig_im, work
!
      logical ::converged
!
      real(dp) :: ddot, dot_ij, dot_ii, dummy, norm, scf_energy
!
      integer(i15) :: i = 0, j = 0, info = 0, iteration = 1
!
!     Get h_xy
!
      call mem%alloc(h, n_ao, n_ao)
      h = zero
!
      call get_ao_xy(h)
!
      do i = 1, n_ao
         do j = 1, n_ao
            write(output%unit,*) 'i j h_ij - h_ji', i, j, h(i,j)-h(j,i)
         enddo
      enddo
!
!     Get g_xyzw
!
      call mem%alloc(g, n_ao**2, n_ao**2)
      g = zero
!
      call get_ao_g_xyzw(g)
!
!     Initialize the molecular orbital coefficients
!
      call mem%alloc(c, n_ao, n_ao)
      c = zero
!
!     Compute the initial density matrix
!
      call mem%alloc(d, n_ao, n_ao)
      d = zero
!  Use soad guess?
      d(1,1) = one
      d(2,2) = one
!
   converged = .false.
   iteration = 1

   do while (.not. converged .and. iteration .lt. 500)
!
!     Get s_xy
!
      call mem%alloc(s, n_ao, n_ao)
      s = zero
!
      call get_ao_s_xy(s)
!
!     Calculate the initial Fock matrix
!
      call mem%alloc(f, n_ao, n_ao)
      f = zero
!
!     Term 1. f_xy = h_xy
!
      call dcopy(n_ao**2, h, 1, f, 1) ! 1E term
!
!     Term 2. f_xy =+ sum_wz g_xywz D_wz
!
      call dgemm('N','N',  &
                  n_ao**2, &
                  1,       &
                  n_ao**2, &
                  one,     &
                  g,       &
                  n_ao**2, &
                  d,       &
                  n_ao**2, &
                  one,     &
                  f,       &
                  n_ao**2)
!
!     Term 3. f_xy =+ - 1/2 sum_wz g_xzwy D_wz
!
!     g_reord(xywz) = g(xzwy)
!             1432      1234
!
      call mem%alloc(g_reord, n_ao**2, n_ao**2)
      g_reord = zero
!
      call sort_1234_to_1432(g, g_reord, n_ao, n_ao, n_ao, n_ao)
!
      call dgemm('N','N',  &
                  n_ao**2, &
                  1,       &
                  n_ao**2, &
                  -half,   &
                  g_reord, &
                  n_ao**2, &
                  d,       &
                  n_ao**2, &
                  one,     &
                  f,       &
                  n_ao**2)
      call mem%alloc(GD, n_ao, n_ao)
      GD = zero
!
!     Term 1. GD_xy =+ 2 sum_wz g_xywz D_wz
!
      call dgemm('N','N',  &
                  n_ao**2, &
                  1,       &
                  n_ao**2, &
                  two,     &
                  g,       &
                  n_ao**2, &
                  d,       &
                  n_ao**2, &
                  zero,     &
                  GD,      &
                  n_ao**2)
!
!     Term 2. f_xy =+ - sum_wz g_xzwy D_wz
!
!     g_reord(xywz) = g(xzwy)
!             1432      1234
!
      call dgemm('N','N',  &
                  n_ao**2, &
                  1,       &
                  n_ao**2, &
                  -one,    &
                  g_reord, &
                  n_ao**2, &
                  d,       &
                  n_ao**2, &
                  one,     &
                  GD,      &
                  n_ao**2)
!
      scf_energy = zero
      call calc_energy(scf_energy, h, d, GD, n_ao)
!
      call mem%dealloc(GD, n_ao, n_ao)
!
      write(output%unit,*) 'HF energy:', scf_energy
!
!     Solve the eigenvalue problem
!
      call mem%alloc(work, 4*n_ao, 1)
      work = zero
!
      call mem%alloc(eig_re, n_ao, 1)
      eig_re = zero
!
      call dsygv(1, 'V', &
                  'L', &
                  n_ao, &
                  f, & ! contains the eigenvectors on exit
                  n_ao, &
                  s, &
                  n_ao, &
                  eig_re, &
                  work, &
                  4*n_ao, &
                  info)
!
      c = f
!
      do i = 1, n_ao
         write(output%unit, *) 'i eig', i, eig_re(i,1)
      enddo
!
!     Update the density matrix with new c's
!
      call mem%alloc(c_occ, n_ao, 1) ! H2
!
      c_occ(:,1) = c(:,1)
!
      call dgemm('N','T', &
                  n_ao,   &
                  n_ao,   &
                  1,      &
                  two,    &
                  c_occ,  &
                  n_ao,   &
                  c_occ,  &
                  n_ao,   &
                  zero,   & ! Update AO density matrix!
                  d,      &
                  n_ao)
!
      call mem%dealloc(c_occ, n_ao, 1) ! H2
!
      do i = 1, n_ao
         write(output%unit, *) 'i d_ii', i, d(i,i)
      enddo
!
!     Calculate the HF energy
!
!     G(D) according to Helgaker and co
!
      call mem%alloc(GD, n_ao, n_ao)
      GD = zero
!
!     Term 1. GD_xy =+ 2 sum_wz g_xywz D_wz
!
      call dgemm('N','N',  &
                  n_ao**2, &
                  1,       &
                  n_ao**2, &
                  two,     &
                  g,       &
                  n_ao**2, &
                  d,       &
                  n_ao**2, &
                  zero,     &
                  GD,      &
                  n_ao**2)
!
!     Term 2. f_xy =+ - sum_wz g_xzwy D_wz
!
!     g_reord(xywz) = g(xzwy)
!             1432      1234
!
      call dgemm('N','N',  &
                  n_ao**2, &
                  1,       &
                  n_ao**2, &
                  -one,    &
                  g_reord, &
                  n_ao**2, &
                  d,       &
                  n_ao**2, &
                  one,     &
                  GD,      &
                  n_ao**2)
!
      scf_energy = zero
      call calc_energy(scf_energy, h, d, GD, n_ao)
!
      write(output%unit,*) 'SCF energy', iteration, scf_energy
!
      call mem%dealloc(f, n_ao, n_ao)
      call mem%dealloc(g_reord, n_ao**2, n_ao**2)
      call mem%dealloc(GD, n_ao, n_ao)
      call mem%dealloc(s_copy, n_ao, n_ao)
!
      iteration = iteration + 1
!
   enddo
!
   end subroutine get_ao_xy_integrals
!
!
   subroutine calc_energy(scf_energy, h, d, GD, n_ao)
!
      implicit none
!
      real(dp) :: scf_energy
!
      integer(i15) :: n_ao, i
!
      real(dp), dimension(:,:) :: h, d, GD
!
      real(dp), dimension(:,:), allocatable :: dh, dGD
!
!     E = Tr (d * h) + 1/4 * Tr (d * GD) + (hnuc)
!
      call mem%alloc(dh, n_ao, n_ao)
!
      call dgemm('N','N', &
                  n_ao, &
                  n_ao, &
                  n_ao, &
                  one, &
                  d, &
                  n_ao, &
                  h, &
                  n_ao, &
                  zero, &
                  dh, &
                  n_ao)
!
      call mem%alloc(dGD, n_ao, n_ao)
!
      call dgemm('N','N', &
                  n_ao, &
                  n_ao, &
                  n_ao, &
                  one, &
                  d, &
                  n_ao, &
                  GD, &
                  n_ao, &
                  zero, &
                  dGD, &
                  n_ao)
!
      scf_energy = zero
!
      do i = 1, n_ao
!
         scf_energy = scf_energy + dh(i,i) + (one/four)*dGD(i,i)
!
      enddo
!
   end subroutine calc_energy
!
! Returns the inverse of a matrix calculated by finding the LU
! decomposition.  Depends on LAPACK.
function inv(A) result(Ainv)
  real(dp), dimension(:,:), intent(in) :: A
  real(dp), dimension(size(A,1),size(A,2)) :: Ainv

  real(dp), dimension(size(A,1)) :: work  ! work array for LAPACK
  integer, dimension(size(A,1)) :: ipiv   ! pivot indices
  integer :: n, info

  ! External procedures defined in LAPACK
  external DGETRF
  external DGETRI

  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  Ainv = A
  n = size(A,1)

  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  call DGETRF(n, n, Ainv, n, ipiv, info)

  if (info /= 0) then
     stop 'Matrix is numerically singular!'
  end if

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  call DGETRI(n, Ainv, n, ipiv, work, n, info)

  if (info /= 0) then
     stop 'Matrix inversion failed!'
  end if
end function inv
!
end module integrals_class
