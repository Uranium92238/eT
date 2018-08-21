module eigen_davidson_tool_class
!
!!
!!    Abstract Davidson solver class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2018
!!
!!
!
   use kinds
   use parameters
!
   use file_class
   use disk_manager_class
   use memory_manager_class
!
   type, extends(davidson_tool) :: eigen_davidson_tool 
!
      real(dp), dimension(:,:), allocatable :: omega_re 
      real(dp), dimension(:,:), allocatable :: omega_im
!
   contains 
!
      procedure :: solve_reduced_problem => solve_reduced_problem_eigen_davidson_tool
!
   end type eigen_davidson_tool
!
contains
!
!
   subroutine solve_reduced_problem_eigen_davidson_tool(davidson)
!!
!!    Solve reduced problem 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    For the current reduced matrix A, this routine solves the eigenvalue problem 
!!
!!       A X_n = omega_n X_n,    n = 1, 2, ..., dim_red.
!!
!!    On exit, the eigenvalues are stored in the vectors omega_re and omega_im,
!!    where the real and imaginary parts of eigenvalue n is stored in the nth 
!!    row of these vectors.
!!
      implicit none 
!
      class(eigen_davidson_tool) :: davidson 
!
      real(dp), dimension(:,:), allocatable :: work
!
      real(dp), dimension(:,:), allocatable :: omega_re
      real(dp), dimension(:,:), allocatable :: omega_im
!
      integer(i15), dimension(:,:), allocatable :: index_list
!
      integer(i15) :: dummy = 0, info = 0, j = 0, i = 0
!
!     -::- Solve reduced eigenvalue problem -::-
!
      call mem%alloc(work, 4*(solver%dim_red), 1)
!      
      work = zero
      info = 0
!
      call mem%alloc(X_red, davidson%dim_red, davidson%dim_red)
      X_red = zero
!
      call mem%alloc(omega_re, solver%dim_red, 1)
      call mem%alloc(omega_im, solver%dim_red, 1)
!
      omega_re = zero
      omega_im = zero
!
      call dgeev('N','V',             &
                  davidson%dim_red,   &
                  A_red,              &
                  davidson%dim_red,   &
                  omega_re,           &
                  omega_im,           &
                  dummy,              &              
                  1,                  &
                  X_red,              &
                  davidson%dim_red,   &
                  work,               &
                  4*davidson%dim_red, &
                  info)
!
      if (info .ne. 0) call call output%error_msg('could not solve reduced equation in Davidson solver.')
!
      call mem%dealloc(work, 4*davidson%reduced_space, 1)
!
!     Find lowest n_solutions eigenvalues and sort them (the corresponding indices
!     are placed in the integer array index_list)
!
      call mem%alloc(davidson%omega_re, solver%n_solutions, 1)
      call mem%alloc(davidson%omega_im, solver%n_solutions, 1)
!
      davidson%omega_re = zero
      davidson%omega_im = zero
!
      call mem%alloc_int(index_list, davidson%n_solutions, 1)
!
      call get_n_lowest(davidson%n_solutions, davidson%dim_red, &
                        eigen_re, davidson%eigen_re, index_list)
!
!     Pick out solution vectors and imaginary parts of eigenvalues according to index_list
!
      call mem%alloc(davidson%X_red, davidson%dim_red, davidson%n_solutions)
      davidson%X_red = zero
!
      do j = 1, davidson%n_solutions
!
         do i = 1, davidson%dim_red
!
            davidson%X_red(i, j) = X_red(i, index_list(j,1))
!
         enddo
!
         davidson%omega_im(j, 1) = omega_im(index_list(j,1), 1)
!
      enddo 
!
      call mem%dealloc(X_red, davidson%reduced_space, solver%reduced_space)
      call mem%dealloc(omega_re, davidson%reduced_space, 1)
      call mem%dealloc(omega_im, davidson%reduced_space, 1)
      call mem%dealloc_int(index_list, davidson%n_solutions, 1)
!
   end subroutine solve_reduced_problem_eigen_davidson_tool
!
!
end module eigen_davidson_tool_class
