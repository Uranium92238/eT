module eigen_davidson_tool_class
!
!!
!!    Eigenvalue davidson tool class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2018
!!
!!    A tool to help solve an eigenvalue equation A X_n = omega_n X_n
!!    for symmetric A using the Davidson algorithm. It is tailored to
!!    be usable also in cases where A cannot be stored, but where the 
!!    transformation X -> A X is implemented. A typical Davidson loop
!!    will use this tool as follows:
!!
!!       To-write.
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
      procedure :: construct_residual    => construct_residual_eigen_davidson_tool
!
      procedure :: construct_re_residual => construct_re_residual_eigen_davidson_tool
      procedure :: construct_im_residual => construct_im_residual_eigen_davidson_tool
!
   end type eigen_davidson_tool
!
contains
!
!
   subroutine solve_reduced_problem_eigen_davidson_tool(davidson)
!!
!!    Solve reduced problem 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
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
   subroutine construct_residual_eigen_davidson_tool(davidson, R, n)
!!
!!    Construct residual 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    Constructs the nth full space residual 
!!
!!       R = A X - omega X,
!!
!!    without any preconditioning (this is applied after a call
!!    to construct residual from the construct next trial vector
!!    routine). Here, X is the nth eigenvector and omega is the 
!!    eigenvalue. 
!!
      implicit none 
!
      class(eigen_davidson_tool), intent(in) :: davidson 
!
      real(dp), dimension(davidson%n_parameters, 1) :: R 
!
      integer(i15), intent(in) :: n
!
      real(dp), dimension(:,:), allocatable :: X 
!
      call mem%alloc(X, davidson%n_parameters, 1)
!
      if (omega_im(n, 1) .eq. zero) then  ! standard case: the nth root is not part of a complex pair
!
         call davidson%construct_X(X, n)  ! set X
         call davidson%construct_AX(R, n) ! set R = AX 
!
         R = R - omega_re(n, 1)*X         ! R = AX - omega*X 
!
      else ! the nth root is part of a complex pair 
!
!        If it's the first root of the pair, construct the real residual; if it's the second, 
!        construct instead the imaginary residual  
!
         if (n .eq. 1) then
!
            if (davidson%n_solutions .lt. 2) call error_msg('add one more root to treat the complex pair.')
!
            call davidson%construct_re_residual(R, n)
!
         elseif (n .eq. davidson%n_solutions) then 
!
            if (omega_re(n-1, 1) .ne. omega_re(n, 1)) call error_msg('add one more root to treat the complex pair.')
!
            call davidson%construct_im_residual(R, n)
!
         else ! neither first or last, so it's safe to look at n + 1 and n - 1 
!
            if (omega_re(n, 1) .eq. omega_re(n - 1, 1)) then
!
               call davidson%construct_im_residual(R, n)
!
            elseif (omega_re(n, 1) .eq. omega_re(n + 1, 1) ) then 
!
               call davidson%construct_re_residual(R, n)
!
            else ! should never happen, but just in case, let the user know there's a bug
!
               call error_msg('something went very wrong when trying to construct imaginary residual.')
!
            endif 
!
         endif
!
      endif
!
      call mem%dealloc(X, davidson%n_parameters)

!
   end subroutine construct_residual_eigen_davidson_tool
!
!
   subroutine construct_re_residual_eigen_davidson_tool(davidson, R, n)
!!
!!    Construct real residual 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
!!    For a complex pair, 
!!
!!       X+ = X_re + i X_im 
!!       X- = X_re - i X_im,
!!
!!    we will find the real part X_re as the first of the two roots,
!!    and the imaginary part X_im as the second of the two roots 
!!    (in X_red). This routine constructs the residual associated
!!    with the real part. As such, the construct residual routine 
!!    makes sure the n-th root is X_re and (n + 1)-th is X_im (and
!!    exists, which might not be the case if too few roots are 
!!    requested). On exit, the real residual is placed in R: 
!!    
!!       R = (A X_re - omega_re X_re) + omega_im X_im 
!!
!!    The residual is divided by the norm of X+ (or, 
!!    equivalently, X-).
!!
      implicit none 
!
      class(eigen_davidson_tool), intent(in) :: davidson 
!
      real(dp), dimension(davidson%n_parameters, 1) :: R 
!
      integer(i15), intent(in) :: n 
!
      real(dp), dimension(:,:), allocatable :: X_re  
      real(dp), dimension(:,:), allocatable :: X_im 
!
      real(dp) :: norm_X_re
      real(dp) :: norm_X_im
!
      call mem%alloc(X_re, davidson%n_parameters, 1)  
      call mem%alloc(X_im, davidson%n_parameters, 1)  
!
      call davidson%construct_X(X_re, n)     ! set X_re 
      call davidson%construct_X(X_im, n + 1) ! set X_im
      call davidson%construct_AX(R, n)       ! set R = A X_re 
!
      R = R - omega_re(n, 1)*X_re + omega_im(n, 1)*X_im 
!
      norm_X_re = get_l2_norm(X_re, davidson%n_parameters)
      norm_X_im = get_l2_norm(X_im, davidson%n_parameters)
!
      R = R/(sqrt(norm_X_re**2 + norm_X_im**2))
!
      call mem%dealloc(X_re, davidson%n_parameters, 1)  
      call mem%dealloc(X_im, davidson%n_parameters, 1) 
!
   end subroutine construct_re_residual_eigen_davidson_tool
!
!
   subroutine construct_im_residual_eigen_davidson_tool(davidson, R, n)
!!
!!    Construct imaginary residual 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
!!    For a complex pair, 
!!
!!       X+ = X_re + i X_im 
!!       X- = X_re - i X_im,
!!
!!    we will find the real part X_re as the first of the two roots,
!!    and the imaginary part X_im as the second of the two roots 
!!    (in X_red). This routine constructs the residual associated
!!    with the imaginary part. As such, the construct residual 
!!    routine makes sure the n-th root is X_im and (n - 1)-th 
!!    is X_re. On exit, the imaginary residual is placed in R: 
!!    
!!       R = (A X_im - omega_re X_im) - omega_im X_re 
!!
!!    The residual is divided by the norm of X+ (or, equivalently, 
!!    X-). Note that omega_im of root n is the negative of omega_im 
!!    of root n - 1. To avoid inconsistency, we use omega_im from 
!!    the previous root (otherwise, we would change the sign of the
!!    last constribution).
!!
      implicit none 
!
      class(eigen_davidson_tool), intent(in) :: davidson 
!
      real(dp), dimension(davidson%n_parameters, 1) :: R 
!
      integer(i15), intent(in) :: n 
!
      real(dp), dimension(:,:), allocatable :: X_re  
      real(dp), dimension(:,:), allocatable :: X_im 
!
      real(dp) :: norm_X_re
      real(dp) :: norm_X_im
!
      call mem%alloc(X_re, davidson%n_parameters, 1)  
      call mem%alloc(X_im, davidson%n_parameters, 1)  
!
      call davidson%construct_X(X_re, n - 1) ! set X_re 
      call davidson%construct_X(X_im, n)     ! set X_im
      call davidson%construct_AX(R, n)       ! set R = A X_im 
!
      R = R - omega_re(n, 1)*X_im - omega_im(n - 1, 1)*X_re
!
      norm_X_re = get_l2_norm(X_re, davidson%n_parameters)
      norm_X_im = get_l2_norm(X_im, davidson%n_parameters)
!
      R = R/(sqrt(norm_X_re**2 + norm_X_im**2))
!
      call mem%dealloc(X_re, davidson%n_parameters, 1)  
      call mem%dealloc(X_im, davidson%n_parameters, 1) 
!
   end subroutine construct_im_residual_eigen_davidson_tool
!
!
end module eigen_davidson_tool_class
