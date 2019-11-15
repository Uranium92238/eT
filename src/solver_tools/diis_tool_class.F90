!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2019 the authors of eT
!
!  eT is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  eT is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program. If not, see <https://www.gnu.org/licenses/>.
!
!
module diis_tool_class
!
!!
!!    DIIS tool class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2018
!!
!!    Objects of this class can be instantiated in order to solve a non-linear set of equations,
!!    e(x) = 0, where e and x are vectors of length n_equations and n_parameters.
!!
!!    The squared norm of the averaged error vector
!!
!!       < e(x) > = sum_k w_k e(x)_k
!!
!!    is minimized with respect to the condition sum_k w_k = 1. The k index
!!    runs over previous iterations, by default limited to the last 8 iterations
!!    (determined by dimension_).
!!
!!    The standard update of the parameters is
!!
!!       x = sum_k=1^n w_k x_k     (*)
!!
!!    Note that the "x_k" sent to the routine may be quasi-Newton update estimates,
!!    so that the DIIS-extrapolated x in (*) is the optimal combination of previous 
!!    quasi-Newton estimates (this is, e.g., the case for CC ground and 
!!    excited states).
!!
!!    See Pulay, P. Convergence acceleration of iterative sequences. 
!!    The case of SCF iteration. Chem. Phys. Lett. 1980, 73, 393−398.
!!
!!    Typical usage: 
!!
!!    do while (.not. converged) 
!!
!!       -> Calculate error vector e given the current parameters x
!!
!!       -> Overwrite x with the x_n you want to use in the extrapolation 
!!          step (*). This can, e.g., be a quasi-Newton update estimate based  
!!          on the current e and x. 
!!
!!       call diis%update(e, x)
!!
!!       -> now, x contains the extrapolated guess (*)
!!
!!    enddo 
!!
!
   use global_out, only : output
!
   use parameters
   use array_utilities, only : zero_array
!
   use sequential_file_class, only : sequential_file
   use record_storer_class, only : record_storer  
   use file_storer_class, only : file_storer  
   use memory_storer_class, only : memory_storer  
!
   use memory_manager_class, only : mem
!
!
   type :: diis_tool
!
      character(len=40) :: name_          ! determines the prefix of all DIIS files
!
      type(sequential_file), private    :: diis_matrix 
      class(record_storer), allocatable :: e_vectors, x_vectors 
!
      integer, private :: iteration       ! Variable keeping track of the current DIIS iteration
!                                         ! Note: defined to increment by +1 each time 'update' is called.
!
      integer, private :: dimension_      ! Standard is 8, though it might be useful to change this value
!
      integer, private :: n_parameters    ! The length of the x vectors
      integer, private :: n_equations     ! The length of the e vectors
!
      logical :: accumulate               ! Variable keeping track of wether we are  
!                                         ! cumulatively building the DIIS matrix or not
!                                         ! Note: default is to accumulate
!
      logical :: erase_history            ! (DIIS space) 1,2,3,4,...,8,1,2,3,... if true 
                                          ! (DIIS space) 1,2,3,4,...,8,8,8,8,... if false 
!
   contains
!
      procedure :: update                          => update_diis_tool   
!
      procedure :: get_current_dim                 => get_current_dim_diis_tool
!
      procedure :: read_x                          => read_x_diis_tool
      procedure :: read_e                          => read_e_diis_tool
!
      procedure :: write_x                         => write_x_diis_tool
      procedure :: write_e                         => write_e_diis_tool
!
      procedure, private :: construct_diis_matrix  => construct_diis_matrix_diis_tool
      procedure, private :: write_current_vecs     => write_current_vecs_diis_tool
!
      final :: destructor 
!
      procedure :: initialize_storers              => initialize_storers_diis_tool
      procedure :: finalize_storers                => finalize_storers_diis_tool
!
   end type diis_tool
!
!
   interface diis_tool 
!
      procedure :: new_diis_tool
!
   end interface diis_tool 
!
!
contains
!
!
   function new_diis_tool(name_, n_parameters, n_equations, &
                        dimension_, accumulate, erase_history) result(diis)
!!
!!    New DIIS tool
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2018
!!
!!    name_:               Name of DIIS tool (used for temporary files)
!!
!!    n_parameters:        Number of parameters (dimensionality of the x vectors)
!!
!!    n_equations:         Number of equations (dimensionality of the e vectors)
!!
!!    dimension_:          Number of records to use in the DIIS extrapolation. 
!!                         Default is to use the previous 8 records. 
!!
!!    accumulate:          If true, the parts of the DIIS matrix constructed 
!!                         on "update" is the parts that have not been calculated
!!                         previously. If false, the entire DIIS matrix is re-
!!                         constructed when "update" is called. This is sometimes
!!                         neccessary because the previous records have to be 
!!                         updated to a new basis. Default is true.
!!
!!    erase_history:       Instead of keeping the previous 8 records 
!!                         (if DIIS dimension = 8), the history is erased 
!!                         every 8th call to "update". This usually leads to 
!!                         slower convergene (or oscillating behavior in the 
!!                         residuals). Default is false.
!!
      type(diis_tool) :: diis
!
      character(len=*), intent(in) :: name_
!
      integer, intent(in) :: n_parameters
      integer, intent(in) :: n_equations
!
      integer, intent(in), optional :: dimension_
      logical, intent(in), optional :: accumulate
      logical, intent(in), optional :: erase_history
!
      diis%name_            = trim(name_)
      diis%n_parameters     = n_parameters
      diis%n_equations      = n_equations
      diis%iteration        = 1 
      diis%dimension_       = 8
!
      diis%accumulate       = .true. 
      diis%erase_history    = .false. 
!
      if (present(accumulate)) then 
!
         diis%accumulate = accumulate
!
      endif
!
      if (present(erase_history)) then
!
         diis%erase_history = erase_history
!
      endif
!
      if (present(dimension_)) then
!
         diis%dimension_ = dimension_
!
      endif 
!
!     Initialize DIIS matrix 
!
      diis%diis_matrix = sequential_file(trim(diis%name_) // '_matrix')
!
   end function new_diis_tool
!
!
   subroutine destructor(diis)
!!
!!    Destructor 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, July 2018
!!
      implicit none
!
      type(diis_tool) :: diis
!
      if (diis%iteration .gt. 1 .and. diis%accumulate) then
         call diis%diis_matrix%delete_()
      endif
!
   end subroutine destructor
!
!
   subroutine finalize_storers_diis_tool(diis)
!!
!!    Cleanup storers 
!!    Written by Eirik F. Kjønstad, Nov 2019 
!!
!!    Performs cleanup for the different storers. 
!!    Cannot be done in the constructor. 
!!
      implicit none 
!
      class(diis_tool) :: diis
!
      call diis%e_vectors%finalize_storer()
      call diis%x_vectors%finalize_storer()
!
   end subroutine finalize_storers_diis_tool
!
!
   subroutine initialize_storers_diis_tool(diis, records_in_memory)
!!
!!    Prepare storers 
!!    Written by Eirik F. Kjønstad, Nov 2019 
!!
!!    Prepares the storer objects for e and x vectors.
!!
!!    records_in_memory:   If true, keep DIIS records (previous x and e vectors) 
!!                         in memory. Otherwise, they are stored on disk.
!!
      implicit none 
!
      class(diis_tool) :: diis
!
      logical, intent(in) :: records_in_memory
!
      if (records_in_memory) then 
!
         call output%printf('Storing DIIS records in memory.', pl='v', fs='(/t3,a)')
!
         diis%e_vectors = memory_storer(trim(diis%name_) // '_e', &
                           diis%n_equations, diis%dimension_)
!
         diis%x_vectors = memory_storer(trim(diis%name_) // '_x', &
                           diis%n_parameters, diis%dimension_)
!
      else 
!
         call output%printf('Storing DIIS records on file.', pl='v', fs='(/t3,a)')
!
         diis%e_vectors = file_storer(trim(diis%name_) // '_e', &
                           diis%n_equations, diis%dimension_, delete=.true.)
         diis%x_vectors = file_storer(trim(diis%name_) // '_x', &
                           diis%n_parameters, diis%dimension_, delete=.true.)
!
      endif 
!
      call diis%e_vectors%initialize_storer()
      call diis%x_vectors%initialize_storer()
!
   end subroutine initialize_storers_diis_tool
!
!
    subroutine update_diis_tool(diis, e, x)
!!
!!    Update 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2018
!!
!!    Routine to call in each iteration of a DIIS-based solver. The 
!!    routine stores the current vectors (e and x) sent to the routine. 
!!    Then it constructs the DIIS matrix G and solves the DIIS equation 
!!
!!       G w = H,   G_ij = e_i^T e_j,    H_i = 0 
!!
!!    for the DIIS weights w. In practice, G and H are padded 
!!    (one extra column and row) to enforce 
!!
!!       sum_i w_i = 1.
!!
!!    On exit, 
!!
!!       x = sum_i w_i x_i.
!!
!!    Based on the CCSD specific routine written by Sarai D. Folkestad
!!    and Eirik F. Kjønstad, May 2017.
!!
      implicit none
!
      class(diis_tool) :: diis
!
      real(dp), dimension(diis%n_equations)  :: e
      real(dp), dimension(diis%n_parameters) :: x
!
      real(dp), dimension(:), allocatable :: x_i ! To hold previous x_i temporarily
!
      integer :: i = 0
!
      integer :: info = -1         ! Error integer for dgesv routine (LU factorization)
      integer :: current_dim 
!
      real(dp), dimension(:), allocatable :: diis_vector
      real(dp), dimension(:,:), allocatable :: diis_matrix
!
      integer, dimension(:), allocatable :: ipiv ! Pivot integers (see dgesv routine)
!
!     Compute the current dimensionality of the problem 
!     (1, 2,..., 7, 8, 8, 8, 8,...) for the standard dimension_ = 8 without erasing history
!
      current_dim = diis%get_current_dim()
!
!     Write current e and x vector to file
!
      call diis%write_current_vecs(e, x, current_dim)
!
!     :: Solve the least squares problem, G w = H
!
!        G : DIIS matrix, G_ij = Δx_i Δx_j,
!        H : DIIS vector,  H_i = 0,
!
!     where i, j = 1, 2, ..., current_dim. To enforce normality
!     of the solution, G is extended with a row & column of -1's
!     and H with a -1 at the end.
!
!     Set the DIIS vector H
!
      call mem%alloc(diis_vector, current_dim + 1)
!
      diis_vector = zero
      diis_vector(current_dim + 1) = -one
!
!     Set the DIIS matrix G
!
      call mem%alloc(diis_matrix, current_dim + 1, current_dim + 1)
!
      call diis%construct_diis_matrix(current_dim, diis_matrix, e)
!
!     Solve the DIIS equation G w = H for the DIIS weights w 
!
!     Note: on exit, the solution is in the diis_vector,
!     provided info = 0 (see LAPACK documentation for more)
!
      call mem%alloc(ipiv, current_dim + 1)
      ipiv = 0
!
      call dgesv(current_dim + 1,  &
                  1,               &
                  diis_matrix,     &
                  current_dim + 1, &
                  ipiv,            &
                  diis_vector,     &
                  current_dim + 1, &
                  info)
!
      call mem%dealloc(ipiv, current_dim + 1)
      call mem%dealloc(diis_matrix, current_dim + 1, current_dim + 1)
!
      if (info .ne. 0) call output%error_msg('could not solve DIIS matrix equation! (i0)', ints=[info])
!
!     Update the parameters x
!
      call zero_array(x, diis%n_parameters)
!
      call mem%alloc(x_i, diis%n_parameters)
!
      do i = 1, current_dim
!
!        Read the x_i vector
!
         call diis%x_vectors%get(x_i, i)
!
!        Add w_i x_i to x 
!
         call daxpy(diis%n_parameters, diis_vector(i), x_i, 1, x, 1)
!
      enddo
!
      call mem%dealloc(x_i, diis%n_parameters)
!
      call mem%dealloc(diis_vector, current_dim + 1)
!
      diis%iteration = diis%iteration + 1
!
   end subroutine update_diis_tool
!
!
   function get_current_dim_diis_tool(diis)
!!
!!    Get current dimension 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
      implicit none
!
      class(diis_tool) :: diis
!
      integer :: get_current_dim_diis_tool
!
      if (diis%iteration .gt. diis%dimension_ .and. &
            .not. diis%erase_history) then 
!
         get_current_dim_diis_tool = diis%dimension_
!
      else
!
         get_current_dim_diis_tool = diis%iteration - &
                           diis%dimension_*((diis%iteration-1)/diis%dimension_) 
!
      endif
!
   end function get_current_dim_diis_tool
!
!
   subroutine write_current_vecs_diis_tool(diis, e, x, current_dim)
!!
!!    Write current vectors 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
!!    Writes the current e and x to file. Each index 1, 2, ..., diis_dim
!!    has its own record. If more iterations have passed than we keep on file
!!    (i.e., diis%iteration > diis%dimension_), then we do a first-in/
!!    last-out replacement unless erase_history is true. 
!!
      implicit none 
!
      class(diis_tool) :: diis
!
      real(dp), dimension(diis%n_equations), intent(in)  :: e
      real(dp), dimension(diis%n_parameters), intent(in) :: x
!
      integer, intent(in) :: current_dim
!
      logical :: cycle_left   
!
      cycle_left = (.not. diis%erase_history) .and. &
                  (current_dim .eq. diis%dimension_) .and. &
                  (current_dim .ne. diis%iteration)
!
      if (cycle_left) then 
!
         call diis%e_vectors%cycle_left() 
         call diis%x_vectors%cycle_left() 
!
      endif 
!
      call diis%e_vectors%set(e, current_dim)
      call diis%x_vectors%set(x, current_dim)
!
   end subroutine write_current_vecs_diis_tool
!
!
   subroutine construct_diis_matrix_diis_tool(diis, current_dim, diis_matrix, e)
!!
!!    Construct DIIS matrix
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
      implicit none
!
      class(diis_tool) :: diis
!
      integer, intent(in) :: current_dim
!
      real(dp), dimension(current_dim + 1, current_dim + 1), intent(inout) :: diis_matrix
!
      real(dp), dimension(diis%n_equations), intent(in) :: e
!
      real(dp), dimension(:), allocatable :: e_i  ! To hold previous e_i temporarily
      real(dp), dimension(:), allocatable :: e_j  ! To hold previous e_j temporarily
!
      integer :: i, j
!
      real(dp), dimension(:,:), allocatable :: diis_matrix_tmp 
!
      real(dp) :: ddot
!
      diis_matrix = zero
!
      if (diis%accumulate) then
!
         call diis%diis_matrix%open_()
!
         if (diis%iteration .gt. diis%dimension_ .and. .not. diis%erase_history) then
!
!           First entry goes out, eliminating the (old) 1st row & column.
!
               call mem%alloc(diis_matrix_tmp, current_dim, current_dim)
!
               call diis%diis_matrix%rewind_()
!
               do j = 1, current_dim
                  do i = 1, current_dim
!
                     call diis%diis_matrix%read_(diis_matrix_tmp(i,j))
!
                  enddo
               enddo
!
               diis_matrix(1 : (current_dim - 1), 1 : (current_dim - 1)) = &
               diis_matrix_tmp(2 : current_dim, 2 : current_dim)      
!
               call mem%dealloc(diis_matrix_tmp, current_dim, current_dim)
!
         elseif (current_dim .gt. 1) then
!
!           Still doing the first set of vectors: just read the previous matrix 
!
            call diis%diis_matrix%rewind_()
!
            do j = 1, current_dim - 1
               do i = 1, current_dim - 1
!
                  call diis%diis_matrix%read_(diis_matrix(i,j))
!
               enddo
            enddo
!
         endif
!
!        Get the parts of the DIIS matrix G not constructed in
!        the previous iterations
!
         call mem%alloc(e_i, diis%n_equations)
!
         do i = 1, current_dim
!
            call diis%e_vectors%get(e_i, i)
!
            diis_matrix(current_dim, i) = ddot(diis%n_equations, e, 1, e_i, 1)
            diis_matrix(i, current_dim) = diis_matrix(current_dim, i)
!
            diis_matrix(current_dim + 1, i) = -one
            diis_matrix(i, current_dim + 1) = -one
!
         enddo
!
         call mem%dealloc(e_i, diis%n_equations)
!
!        Write the current DIIS matrix to file
!
         call diis%diis_matrix%rewind_()
!
         do j = 1, current_dim
            do i = 1, current_dim
!
               call diis%diis_matrix%write_(diis_matrix(i,j))
!
            enddo
         enddo
!
         call diis%diis_matrix%close_() 
!
      else ! Non-cumulative diis -> build diis matrix from scratch
!         
         call mem%alloc(e_i, diis%n_equations)
         call mem%alloc(e_j, diis%n_equations)
!
         do i = 1, current_dim
!
            call diis%e_vectors%get(e_i, i)
!
            do j = 1, i
!
               call diis%e_vectors%get(e_j, j)
!
               diis_matrix(i, j) = ddot(diis%n_equations, e_i, 1, e_j, 1)
               diis_matrix(j, i) = diis_matrix(i, j)
!
            enddo
!
            diis_matrix(current_dim + 1, i) = -one
            diis_matrix(i, current_dim + 1) = -one
!
         enddo 
!
         call mem%dealloc(e_i, diis%n_equations)
         call mem%dealloc(e_j, diis%n_equations)
!
      endif
!
   end subroutine construct_diis_matrix_diis_tool
!
!
   subroutine read_x_diis_tool(diis, x_i, i)
!!
!!    Get x_i 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Mar 2019
!!
      implicit none
!
      class(diis_tool) :: diis
!
      real(dp), dimension(diis%n_parameters), intent(out) :: x_i
!
      integer, intent(in) :: i
!
      integer :: current_dim
!
      current_dim = diis%get_current_dim()
!
      if (i .gt. current_dim) call output%error_msg('Asked for an x not in use.')
!
      call diis%x_vectors%get(x_i, i)
!
   end subroutine read_x_diis_tool
!
!
   subroutine read_e_diis_tool(diis, e_i, i)
!!
!!    Get e_i
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Mar 2019
!!
      implicit none
!
      class(diis_tool) :: diis
!
      real(dp), dimension(diis%n_equations), intent(out) :: e_i
!
      integer, intent(in) :: i
!
      integer :: current_dim
!
      current_dim = diis%get_current_dim()
!
      if (i .gt. current_dim) call output%error_msg('Asked for an e not in use.')
!
      call diis%e_vectors%get(e_i, i)
!
   end subroutine read_e_diis_tool
!
!
   subroutine write_x_diis_tool(diis, x_i, i)
!!
!!    Set x_i 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Mar 2019
!!
      implicit none
!
      class(diis_tool) :: diis
!
      real(dp), dimension(diis%n_parameters), intent(in) :: x_i
!
      integer, intent(in) :: i
!
      integer :: current_dim
!
      current_dim = diis%get_current_dim()
!
      if (diis%accumulate) call output%error_msg('Can not set x for cumulative diis.')
      if (i .gt. current_dim) call output%error_msg('Asked for an x not in use.')
!
      call diis%x_vectors%set(x_i, i)
!
   end subroutine write_x_diis_tool
!
!
   subroutine write_e_diis_tool(diis, e_i, i)
!!
!!    Set e_i 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Mar 2019
!!
      implicit none
!
      class(diis_tool) :: diis
!
      real(dp), dimension(diis%n_equations), intent(in) :: e_i
!
      integer, intent(in) :: i
!
      integer :: current_dim
!
      current_dim = diis%get_current_dim()
!
      if (diis%accumulate) call output%error_msg('Can not set e for cumulative diis.')
      if (i .gt. current_dim) call output%error_msg('Asked for an e not in use.')
!
      call diis%e_vectors%set(e_i, i)
!
   end subroutine write_e_diis_tool
!
!
end module diis_tool_class
