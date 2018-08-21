module davidson_tool_class
!
!!
!!    Abstract Davidson davidson class module
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
!
   type, abstract :: davidson_tool
!
      real(dp), dimension(:,:), allocatable :: A_red
      real(dp), dimension(:,:), allocatable :: X_red
!
      type(file) :: X, trials, transforms
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
   contains
!
!     Read and write routines
!
      procedure, non_overridable :: read    => read_davidson_tool
      procedure, non_overridable :: write   => write_davidson_tool
!
      procedure, non_overridable :: read_trials    => read_trials_davidson_tool
      procedure, non_overridable :: write_trials   => write_trials_davidson_tool
!
      procedure, non_overridable :: read_transforms  => read_transforms_davidson_tool
      procedure, non_overridable :: write_transforms => write_transforms_davidson_tool  
!
!     Other procedures
!
      procedure, non_overridable :: construct_reduced_matrix => construct_reduced_matrix_davidson_tool
!
!     Deferred routines
!
      procedure(essential_davidson), deferred :: initialize
      procedure(essential_davidson), deferred :: finalize   
!
      procedure(essential_davidson), deferred :: solve_reduced_problem
      procedure(essential_davidson), deferred :: construct_residual
      procedure(essential_davidson), deferred :: construct_solution
!
   end type davidson_tool
!
!
   abstract interface
!
      subroutine essential_davidson(davidson)
!
         import :: davidson_tool
!
         implicit none 
!
         class(davidson_tool) :: davidson 
!
      end subroutine essential_davidson
!
   end interface
!
contains
!
!
   subroutine read_davidson_tool(davidson)
!!
!!
!!
      implicit none
!
      class(davidson_tool) :: davidson
!
   end subroutine read_davidson_tool
!
!
   subroutine write_davidson_tool(davidson)
!!
!!
!!
      implicit none
!
      class(davidson_tool) :: davidson
!
   end subroutine write_davidson_tool
!
!
   subroutine read_trials_davidson_tool(davidson)
!!
!!
!!
      implicit none
!
      class(davidson_tool) :: davidson
!
   end subroutine read_trials_davidson_tool
!
!
   subroutine write_trials_davidson_tool(davidson)
!!
!!
!!
      implicit none
!
      class(davidson_tool) :: davidson
!
   end subroutine write_trials_davidson_tool
!
!
   subroutine read_transforms_davidson_tool(davidson)
!!
!!
!!
      implicit none
!
      class(davidson_tool) :: davidson
!
   end subroutine read_transforms_davidson_tool
!
!
   subroutine write_transforms_davidson_tool(davidson)
!!
!!
!!
      implicit none
!
      class(davidson_tool) :: davidson
!
   end subroutine write_transforms_davidson_tool
!
!
   subroutine construct_reduced_matrix_davidson_tool(davidson)
!!
!!    Construct reduced matrix
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, June 2018
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
      real(dp) :: ddot
!
      real(dp), dimension(:,:), allocatable :: A_red_copy, c_i, rho_j
!
      integer(i15) :: i, j, ioerror
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
      if (davidson%dim_red .eq. davidson%n_solutions) then ! First iteration
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
         do i = 1, davidson%dim_red
!
            read(davidson%trials%unit, iostat=ioerror) c_i
            if (ioerror .ne. 0) call output%error_msg('reading trial vector file.')
!
            do j = 1, davidson%dim_red
!
               read(davidson%transforms%unit, iostat=ioerror) rho_j
               if (ioerror .ne. 0) call output%error_msg('reading transformed vector file.')
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
end module davidson_tool_class
