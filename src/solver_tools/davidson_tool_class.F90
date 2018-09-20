module davidson_tool_class
!
!!
!!    Abstract Davidson davidson class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2018
!!
!
   use kinds
   use parameters
!
   use file_class
   use disk_manager_class
   use memory_manager_class
   use array_utilities
!
!
   type, abstract :: davidson_tool
!
      character(len=40) :: name 
!
      real(dp), dimension(:,:), allocatable :: A_red
      real(dp), dimension(:,:), allocatable :: X_red
!
      type(file) :: X, trials, transforms, preconditioner
!
      integer(i15) :: dim_red
      integer(i15) :: max_dim_red
!
      integer(i15) :: n_parameters
      integer(i15) :: n_solutions
      integer(i15) :: n_new_trials
!
      real(dp) :: residual_threshold
!
      logical :: do_precondition
!
   contains
!
!     Read and write routines
!
      procedure, non_overridable :: read_trial    => read_trial_davidson_tool
      procedure, non_overridable :: write_trial   => write_trial_davidson_tool
!
      procedure, non_overridable :: read_transform  => read_transform_davidson_tool
      procedure, non_overridable :: write_transform => write_transform_davidson_tool  
!
!     Other procedures
!
      procedure, non_overridable :: construct_reduced_matrix => construct_reduced_matrix_davidson_tool
      procedure, non_overridable :: construct_X              => construct_X_davidson_tool
      procedure, non_overridable :: construct_AX             => construct_AX_davidson_tool
!
      procedure :: set_preconditioner               => set_preconditioner_davidson_tool
      procedure :: precondition                     => precondition_davidson_tool
      procedure :: orthogonalize_against_trial_vecs => orthogonalize_against_trial_vecs_davidson_tool
!
      procedure :: set_A_red => set_A_red_davidson_tool
      procedure :: get_A_red => get_A_red_davidson_tool
!
      procedure :: add_initial_trial_vec => add_initial_trial_vec_davidson_tool
!
      procedure :: set_trials_to_solutions => set_trials_to_solutions_davidson_tool
!
!     Deferred routines  
!
      procedure(solve_reduced_problem), deferred :: solve_reduced_problem
      procedure(construct_residual), deferred    :: construct_residual
!
      procedure :: read_max_dim_red => read_max_dim_red_davidson_tool
!
   end type davidson_tool
!
!
   abstract interface
!
!
      subroutine solve_reduced_problem(davidson)
!
         import :: davidson_tool
!
         implicit none 
!
         class(davidson_tool) :: davidson 
!
      end subroutine solve_reduced_problem
!
!
      subroutine construct_residual(davidson, R, X, norm_X, n)
!
         import :: davidson_tool, dp, i15
!
         implicit none 
!
         class(davidson_tool), intent(in) :: davidson 
!
         real(dp), dimension(davidson%n_parameters, 1)             :: R 
         real(dp), dimension(davidson%n_parameters, 1), intent(in) :: X 
!
         real(dp), intent(in) :: norm_X 
!
         integer(i15), intent(in) :: n
!
      end subroutine construct_residual
!
   end interface
!
contains
!
!
   subroutine read_trial_davidson_tool(davidson, c_i, n)
!!
!!    Read trial vector
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2018
!!
!!    Read n'th trial vector from file
!!
      implicit none
!
      class(davidson_tool) :: davidson
!
      real(dp), dimension(davidson%n_parameters, 1) :: c_i
!
      integer(i15) :: n
!
      integer(i15) :: ioerror
!
      call disk%open_file(davidson%trials, 'read')
!
      call davidson%trials%prepare_to_read_line(n)
!
      read(davidson%trials%unit, iostat = ioerror) c_i
      if (ioerror .ne. 0) call output%error_msg('reading trial vectors file.')
!
      call disk%close_file(davidson%trials)
!
   end subroutine read_trial_davidson_tool
!
!
   subroutine write_trial_davidson_tool(davidson, c_i, position)
!!
!!    Write trial
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2018
!!
!!    Write trial to file
!!
!!    Optional argument position must be either 'rewind' or 'append'
!!
      implicit none
!
      class(davidson_tool) :: davidson
!
      real(dp), dimension(davidson%n_parameters, 1) :: c_i
!
      character(len=*), optional :: position
!
      integer(i15) :: ioerror 
!
!     Was position passed ?
!
      if (present(position)) then
!
!        Sanity check on position variable
!
         if ((trim(position) .ne. 'rewind') .and. (trim(position) .ne. 'append')) then
!
            call output%error_msg('Position specifier not recognized!')
!
         endif
!
         call disk%open_file(davidson%trials, 'write', position)
!
      else
!
         call disk%open_file(davidson%trials, 'write', 'append')
!
      endif 
!
      write(davidson%trials%unit, iostat = ioerror) c_i
      if (ioerror .ne. 0) call output%error_msg('writing trial vectors file.')
!
      call disk%close_file(davidson%trials)
!
   end subroutine write_trial_davidson_tool
!
!
   subroutine add_initial_trial_vec_davidson_tool(davidson, c)
!!
!!    Add initial trial vector 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    Adds a trial vector to the search space, though it is not
!!    assumed to be orthogonal to existing trial vectors on entry.
!!    Components along existing vectors are removed before the 
!!    vector is stored to file.
!!
      implicit none 
!
      class(davidson_tool) :: davidson 
!
      real(dp), dimension(davidson%n_parameters, 1) :: c 
!
      real(dp) :: ddot, projection_of_c_on_c_i, norm_c
!
      integer(i15) :: i 
!
      real(dp), dimension(:,:), allocatable :: c_i
!
      if (davidson%dim_red .eq. 0) then 
!
         norm_c = get_l2_norm(c, davidson%n_parameters)
         c = c/norm_c
!
         call disk%open_file(davidson%trials, 'write', 'rewind')
         write(davidson%trials%unit) c 
         call disk%close_file(davidson%trials)
!
      else 
!
         call mem%alloc(c_i, davidson%n_parameters, 1)
!
         do i = 1, davidson%dim_red
!
            call davidson%read_trial(c_i, i)
            projection_of_c_on_c_i = ddot(davidson%n_parameters, c_i, 1, c, 1)
!
            c = c - projection_of_c_on_c_i*c_i 
!
         enddo 
!
         norm_c = get_l2_norm(c, davidson%n_parameters)
         c = c/norm_c
!
         call disk%open_file(davidson%trials, 'write', 'append')
         write(davidson%trials%unit) c 
         call disk%close_file(davidson%trials)
!
         call mem%dealloc(c_i, davidson%n_parameters, 1)
!
      endif
!
      davidson%dim_red = davidson%dim_red + 1
!
   end subroutine add_initial_trial_vec_davidson_tool
!
!
   subroutine read_transform_davidson_tool(davidson, rho_i, n)
!!
!!    Read transformed vector
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2018
!!
!!    Read n'th transformed vector from file
!!
      implicit none
!
      class(davidson_tool) :: davidson
!
      real(dp), dimension(davidson%n_parameters, 1) :: rho_i
!
      integer(i15) :: n
!
      integer(i15) :: ioerror
!
      call disk%open_file(davidson%transforms, 'read')
!
      call davidson%transforms%prepare_to_read_line(n)
!
      read(davidson%transforms%unit, iostat = ioerror) rho_i
      if (ioerror .ne. 0) call output%error_msg('reading transformed vectors file 2.')
!
      call disk%close_file(davidson%transforms)
!
   end subroutine read_transform_davidson_tool
!
!
   subroutine write_transform_davidson_tool(davidson, rho_i, position)
!!
!!    Write transform
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2018
!!
!!    Write transformed vector to file
!!
!!    Optional argument position must be either 'rewind' or 'append'
!!
      implicit none
!
      class(davidson_tool) :: davidson
!
      real(dp), dimension(davidson%n_parameters, 1) :: rho_i
!
      character(len=*), optional, intent(in) :: position 
!
      integer(i15) :: ioerror   
!
!     Was position passed ?
!
      if (present(position)) then
!
!        Sanity check on position variable
!
         if ((trim(position) .ne. 'rewind') .and. (trim(position) .ne. 'append')) then
!
            call output%error_msg('Position specifier not recognized!')
!
         endif
!
         call disk%open_file(davidson%transforms, 'write', position)
!
      else
!
         call disk%open_file(davidson%transforms, 'write', 'append')
!
      endif 
!
      write(davidson%transforms%unit, iostat = ioerror) rho_i
      if (ioerror .ne. 0) call output%error_msg('writing transformed vector to file.')
!
      call disk%close_file(davidson%transforms)
!
   end subroutine write_transform_davidson_tool
!
!
   subroutine construct_reduced_matrix_davidson_tool(davidson, entire)
!!
!!    Construct reduced matrix
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2018
!!
!!    Constructs the reduced matrix 
!!
!!       A_red_ij = c_i^T * A c_j
!!                = c_i^T * rho_j,
!!      
!!    by reading the trial vectors c_i and transformed vectors rho_j from file and writes
!!    to file.
!!
      implicit none
!
      class(davidson_tool) :: davidson
!
      logical, optional, intent(in) :: entire ! Construct the entire matrix 
!
      real(dp) :: ddot
!
      logical :: entire_local
!
      real(dp), dimension(:,:), allocatable :: A_red_copy, c_i, rho_j
!
      integer(i15) :: i, j, ioerror
!
      if (present(entire)) then
         entire_local = entire
      else
         entire_local = .false.
      endif
!
!     Ask disk manager to open the files
!
      call disk%open_file(davidson%trials, 'read')
      call disk%open_file(davidson%transforms, 'read')
!
!     Rewind files
!
      rewind(davidson%trials%unit)
      rewind(davidson%transforms%unit)
!
!     Allocate c and rho
!
      call mem%alloc(c_i, davidson%n_parameters, 1)
      call mem%alloc(rho_j, davidson%n_parameters, 1)
!
!     Construct reduced matrix: A_red_ij = c_i^T * A * c_j = c_i^T * rho_i
!
      if (davidson%dim_red .eq. davidson%n_solutions .or. entire_local) then ! First iteration
!
!        Make the entire reduced matrix
!
         call mem%alloc(davidson%A_red, davidson%dim_red, davidson%dim_red)
!
         do i = 1, davidson%dim_red
!
            read(davidson%trials%unit, iostat=ioerror) c_i
            if (ioerror .ne. 0) call output%error_msg('reading trial vector file.')
!
            rewind(davidson%transforms%unit)
            do j = 1, davidson%dim_red
!
               read(davidson%transforms%unit, iostat=ioerror) rho_j
               if (ioerror .ne. 0) call output%error_msg('reading transformed vector file.')
!
               davidson%A_red(i,j) = ddot(davidson%n_parameters, c_i, 1, rho_j, 1)
!
            enddo
         enddo
!
      else
!
!        Pad previous A_red
!
         call mem%alloc(A_red_copy, davidson%dim_red - davidson%n_new_trials, davidson%dim_red - davidson%n_new_trials)
         A_red_copy = davidson%A_red
!
         call mem%dealloc(davidson%A_red, davidson%dim_red - davidson%n_new_trials, davidson%dim_red - davidson%n_new_trials)
!
         call mem%alloc(davidson%A_red, davidson%dim_red, davidson%dim_red)
!
         davidson%A_red(1:davidson%dim_red - davidson%n_new_trials, 1:davidson%dim_red - davidson%n_new_trials) = A_red_copy
!
         call mem%dealloc(A_red_copy, davidson%dim_red - davidson%n_new_trials, davidson%dim_red - davidson%n_new_trials)
!
         rewind(davidson%trials%unit)
!
         do i = 1, davidson%dim_red
!
            read(davidson%trials%unit, iostat=ioerror) c_i
            if (ioerror .ne. 0) call output%error_msg('reading trial vector file.')
!
            
            rewind(davidson%transforms%unit)
!
            do j = 1, davidson%dim_red
!
               read(davidson%transforms%unit, iostat=ioerror) rho_j
               if (ioerror .ne. 0) call output%error_msg('reading transformed vector file 3.')
!
               if (j .le. davidson%dim_red - davidson%n_new_trials) then
!
                  if (i .gt. davidson%dim_red - davidson%n_new_trials) then 
!
                     davidson%A_red(i,j) = ddot(davidson%n_parameters, c_i, 1, rho_j, 1)
!
                  endif
!
               elseif (j .gt. davidson%dim_red - davidson%n_new_trials) then
!
                  davidson%A_red(i,j) = ddot(davidson%n_parameters, c_i, 1, rho_j, 1)
!
               endif
!
            enddo
!
         enddo
!
      endif
!
      call mem%dealloc(c_i, davidson%n_parameters, 1)
      call mem%dealloc(rho_j, davidson%n_parameters, 1)
!
!     Close files for trial vectors and transformed vectors
!
      call disk%close_file(davidson%trials)
      call disk%close_file(davidson%transforms)
!
   end subroutine construct_reduced_matrix_davidson_tool
!
!
   subroutine construct_X_davidson_tool(davidson, X, n)
!!
!!    Construct X
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2018
!!
!!       X_n = sum_i c_i (X_red_n)_i
!!
      implicit none
!
      class(davidson_tool) :: davidson
!
      integer(i15) :: n
!
      real(dp), dimension(davidson%n_parameters, 1) :: X
!
      real(dp), dimension(:,:), allocatable :: c_i
!
      integer(i15) :: i, ioerror
!
      X = zero
!
      call mem%alloc(c_i, davidson%n_parameters, 1)
!
      call disk%open_file(davidson%trials, 'read')
      rewind(davidson%trials%unit)
!   
      do i = 1, davidson%dim_red
! 
         read(davidson%trials%unit, iostat = ioerror) c_i
         if (ioerror .ne. 0) call output%error_msg('reading trial vector file.')
!
         call daxpy(davidson%n_parameters, davidson%X_red(i, n), c_i, 1, X, 1)
!
      enddo    
!
      call mem%dealloc(c_i, davidson%n_parameters, 1)
      call disk%close_file(davidson%trials)
!
   end subroutine construct_X_davidson_tool
!
!
   subroutine construct_AX_davidson_tool(davidson, AX, n)
!!
!!    Construct X
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2018
!!
!!       AX_n = sum_i rho_i (X_red_n)_i
!!
      implicit none
!
      class(davidson_tool) :: davidson
!
      integer(i15) :: n
!
      real(dp), dimension(davidson%n_parameters, 1) :: AX
!
      real(dp), dimension(:,:), allocatable :: rho_i
!
      integer(i15) :: i, ioerror
!
      AX = zero
!
      call mem%alloc(rho_i, davidson%n_parameters, 1)
!
      call disk%open_file(davidson%transforms, 'read')
      rewind(davidson%transforms%unit)
!   
      do i = 1, davidson%dim_red
! 
         read(davidson%transforms%unit, iostat = ioerror) rho_i
         if (ioerror .ne. 0) call output%error_msg('reading transformed vector file.')
!
         call daxpy(davidson%n_parameters, davidson%X_red(i, n), rho_i, 1, AX, 1)
!
      enddo    
!
      call mem%dealloc(rho_i, davidson%n_parameters, 1)
!
      call disk%close_file(davidson%transforms)
!
   end subroutine construct_AX_davidson_tool
!
!
   subroutine set_preconditioner_davidson_tool(davidson, preconditioner)
!!
!!    Set preconditioner 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
!!    This routine saves the diagonal preconditioner to file. The 
!!    assumption being that 
!!
!!       preconditioner(i) ~ A(i,i),
!!
!!    where A is the coefficient matrix. The inverse of this diagonal 
!!    matrix is a good approximation of A if A is diagonally dominant
!!    to some extent.
!!
      implicit none 
!
      class(davidson_tool) :: davidson
!
      real(dp), dimension(davidson%n_parameters, 1), intent(in) :: preconditioner 
!
      call disk%open_file(davidson%preconditioner, 'write', 'rewind')
      write(davidson%preconditioner%unit) preconditioner 
      call disk%close_file(davidson%preconditioner)
!
      davidson%do_precondition = .true.
!
   end subroutine set_preconditioner_davidson_tool
!
!
   subroutine precondition_davidson_tool(davidson, R)
!!
!!    Precondition 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
!!    Preconditions the vector R, i.e. 
!!
!!       R(i) <- R(i)/preconditioner(i)
!!
!!    However, if the user has not set any preconditioner, 
!!    this routine performs no action on R. 
!!
      implicit none 
!
      class(davidson_tool) :: davidson
!
      real(dp), dimension(davidson%n_parameters, 1), intent(inout) :: R
!
      real(dp), dimension(:,:), allocatable :: preconditioner
!
      integer(i15) :: i 
!
!
      if (davidson%do_precondition) then 
!
         call mem%alloc(preconditioner, davidson%n_parameters, 1)
!
         call disk%open_file(davidson%preconditioner, 'read')
         rewind(davidson%preconditioner%unit)
         read(davidson%preconditioner%unit) preconditioner
         call disk%close_file(davidson%preconditioner)
!
         do i = 1, davidson%n_parameters
!
            R(i, 1) = R(i, 1)/preconditioner(i, 1)
!
         enddo 
!
         call mem%dealloc(preconditioner, davidson%n_parameters, 1)
!
      endif 
!
   end subroutine precondition_davidson_tool
!
!
   subroutine orthogonalize_against_trial_vecs_davidson_tool(davidson, R)
!!
!!    Orthogonalize against trial vectors 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    Orthogonalizes R against the existing trial vectors. Note 
!!    that this is usually done residual after residual, where 
!!    each new residual is orthogonalized against the previous 
!!    (where the number of previous changes with 'n_new_trials').
!!
!!    Note that it does not normalize R.
!!
      implicit none 
!
      class(davidson_tool) :: davidson 
!
      real(dp), dimension(davidson%n_parameters, 1) :: R 
!
      real(dp) :: ddot, projection_of_R_on_c_i
!
      integer(i15) :: i 
!
      real(dp), dimension(:,:), allocatable :: c_i
!
      call mem%alloc(c_i, davidson%n_parameters, 1)
!
      do i = 1, davidson%dim_red + davidson%n_new_trials
!
         call davidson%read_trial(c_i, i)
         projection_of_R_on_c_i = ddot(davidson%n_parameters, c_i, 1, R, 1)
!
         call daxpy(davidson%n_parameters, - projection_of_R_on_c_i, c_i, 1, R, 1)
!
      enddo 
!
      call mem%dealloc(c_i, davidson%n_parameters, 1)
!
   end subroutine orthogonalize_against_trial_vecs_davidson_tool
!
!
   subroutine set_A_red_davidson_tool(davidson, A_red)
!!
!!    Set A reduced 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
!!    Sets the reduced coefficient matrix A_red equal to the
!!    array sent to the routine. This is a non-standard routine 
!!    that should only be used if the matrix e.g. depends on 
!!    a parameter (scaling certain elements, for instance),
!!    that must be handled outside the Davidson object. 
!!
!!    Usually, construction of A_red is done solely by 
!!    the Davidson object, based on trial vectors and 
!!    transformed vectors on file. 
!!
      implicit none 
!
      class(davidson_tool) :: davidson 
!
      real(dp), dimension(davidson%dim_red, davidson%dim_red), intent(in) :: A_red 
!
      davidson%A_red = A_red 
!
   end subroutine set_A_red_davidson_tool
!
!
   subroutine get_A_red_davidson_tool(davidson, A_red)
!!
!!    Get A reduced 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(davidson_tool) :: davidson 
!
      real(dp), dimension(davidson%dim_red, davidson%dim_red), intent(inout) :: A_red 
!
      A_red = davidson%A_red 
!
   end subroutine get_A_red_davidson_tool
!  
!
   subroutine set_trials_to_solutions_davidson_tool(davidson)
!!
!!    Set trials to solutions
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
!!    Set trials equal to solution on file.
!!    This should be used when max_dim_red is reached and 
!!    if we restart from old solutions
!!
      implicit none
!
      class(davidson_tool) :: davidson 
!
      real(dp), dimension(:,:), allocatable :: X
!
      integer(i15) :: solution
!
      call disk%open_file(davidson%trials, 'write', 'rewind')
      call disk%open_file(davidson%X, 'read')
      rewind(davidson%X%unit)
!
      call mem%alloc(X, davidson%n_parameters, 1)
!
      do solution = 1, davidson%n_solutions
!
         read(davidson%X%unit) X 
         write(davidson%trials%unit) X
!
      enddo
!
      call mem%dealloc(X, davidson%n_parameters, 1)
!
      rewind(davidson%trials%unit)
!
      call disk%close_file(davidson%trials)
      call disk%close_file(davidson%X)
!
!     Delete transformed vectors file, if it is there
!
      call disk%open_file(davidson%transforms, 'write', 'rewind')
      call disk%close_file(davidson%transforms, 'delete')
!
!     Delete A_red if it is there
!
      if (allocated(davidson%A_red)) call mem%dealloc(davidson%A_red, &
         davidson%dim_red - davidson%n_new_trials, davidson%dim_red - davidson%n_new_trials)
!
      davidson%dim_red = davidson%n_solutions
      davidson%n_new_trials = davidson%n_solutions
!
   end subroutine set_trials_to_solutions_davidson_tool
!
!
   subroutine read_max_dim_red_davidson_tool(davidson)
!!
!!    Read max reduced dimension
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none
!
      class(davidson_tool) :: davidson 
!
      character(len=100) :: line
!
      rewind(input%unit)
!
      do 
!
         read(input%unit, '(a100)') line
         line = remove_preceding_blanks(line)
!
         if (trim(line) .eq. 'end geometry') then
            backspace(input%unit)
            return
         endif
!
         if (line(1:22) == 'max reduced dimension:') then
! 
            read(line(23:100), *) davidson%max_dim_red
            return
! 
         endif
!
      enddo
!
   end subroutine read_max_dim_red_davidson_tool
!
end module davidson_tool_class
