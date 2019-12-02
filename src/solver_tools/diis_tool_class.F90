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
!!    In the tool, the DIIS matrix is denoted as G:
!!
!!       G_ij = e_i^T e_j,     
!!
!!    where i and j denote the record numbers.
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
      class(record_storer), allocatable :: e_vectors, x_vectors 
!
      integer, private :: iteration       ! Variable keeping track of the current DIIS iteration
                                          ! Note: defined to increment by +1 each time 'update' 
                                          ! is called.
!
      integer, private :: dimension_      ! Standard is 8, though it might be useful to 
                                          ! change this value
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
      logical :: crop                     ! Conjugate residual with optimal trial vectors 
                                          ! J. Chem. Theory Comput. 2015, 11, 4, 1518-1524
!
      real(dp), dimension(:,:), allocatable :: G ! DIIS matrix 
!
   contains
!
      procedure :: update                          => update_diis_tool   
!
      procedure :: read_x                          => read_x_diis_tool
      procedure :: read_e                          => read_e_diis_tool
!
      procedure :: write_x                         => write_x_diis_tool
      procedure :: write_e                         => write_e_diis_tool
!
      procedure :: initialize_storers              => initialize_storers_diis_tool
      procedure :: finalize_storers                => finalize_storers_diis_tool
!
      procedure :: get_dim_G                       => get_dim_G_diis_tool
!
      procedure, private :: construct_padded_G     => construct_padded_G_diis_tool
      procedure, private :: write_e_and_x          => write_e_and_x_diis_tool
      procedure, private :: update_crop_e_and_x    => update_crop_e_and_x_diis_tool 
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
                  dimension_, accumulate, erase_history, crop) result(diis)
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
!!    crop:                Use conjugate residual with optimal trial vectors 
!!                         (see J. Chem. Theory Comput. 2015, 11, 4, 1518-1524). Optional.
!!                         Default is false. 
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
      logical, intent(in), optional :: crop
!
      diis%name_            = trim(name_)
      diis%n_parameters     = n_parameters
      diis%n_equations      = n_equations
 !
      diis%iteration        = 1 
      diis%dimension_       = 8
!
      diis%accumulate       = .true. 
      diis%erase_history    = .false. 
      diis%crop             = .false. 
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
      if (present(crop)) then 
!
         diis%crop = crop
!
      endif 
!
      if (diis%crop .and. diis%erase_history) then 
!
         call output%warning_msg("Asked for both 'crop' and 'erase history' in DIIS tool. " // &
                                 "Convergence will not improve by using CROP in this case.")
!
      endif
!
   end function new_diis_tool
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
      call mem%dealloc(diis%G, diis%dimension_, diis%dimension_)
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
      call mem%alloc(diis%G, diis%dimension_, diis%dimension_)
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
      real(dp), dimension(:), allocatable :: x_i   ! ith parameter vector x_i 
!
      integer :: i = 0
!
      integer :: info = -1 
      integer :: dim_G 
!
      real(dp), dimension(:), allocatable   :: padded_H
      real(dp), dimension(:,:), allocatable :: padded_G
!
      integer, dimension(:), allocatable :: ipiv ! Pivot integers (see dgesv routine)
!
!     Compute the current dimensionality of the problem 
!     (1, 2,..., 7, 8, 8, 8, 8,...) for the standard dimension_ = 8 without erasing history
!
      dim_G = diis%get_dim_G()
!
!     Write current e and x vector to file
!
      call diis%write_e_and_x(e, x, dim_G)
!
!     :: Solve the least squares problem, G w = H
!
!        G : DIIS matrix, G_ij = e_i^T e_j,
!        H : DIIS vector,  H_i = 0,
!
!     where i, j = 1, 2, ..., dim_G. To enforce normality
!     of the solution, G is extended with a row & column of -1's
!     and H with a -1 at the end.
!
!     Set the DIIS vector H
!
      call mem%alloc(padded_H, dim_G + 1)
!
      padded_H = zero
      padded_H(dim_G + 1) = -one
!
!     Set the padded DIIS matrix G
!
      call mem%alloc(padded_G, dim_G + 1, dim_G + 1)
!
      call diis%construct_padded_G(dim_G, padded_G, e)
!
!     Solve the DIIS equation G w = H for the DIIS weights w 
!
!     Note: on exit, the solution is in the padded_H,
!     provided info = 0 (see LAPACK documentation for more)
!
      call mem%alloc(ipiv, dim_G + 1)
      ipiv = 0
!
      call dgesv(dim_G + 1,         &
                  1,                &
                  padded_G,      &
                  dim_G + 1,        &
                  ipiv,             &
                  padded_H,         &
                  dim_G + 1,        &
                  info)
!
      call mem%dealloc(ipiv, dim_G + 1)
      call mem%dealloc(padded_G, dim_G + 1, dim_G + 1)
!
      if (info .ne. 0) &
         call output%error_msg('could not solve DIIS matrix equation! (i0)', ints=[info])
!
!     Update the parameters x
!
      call zero_array(x, diis%n_parameters)
!
      call mem%alloc(x_i, diis%n_parameters)
!
      do i = 1, dim_G
!
!        Read the x_i vector
!
         call diis%x_vectors%get(x_i, i)
!
!        Add w_i x_i to x 
!
         call daxpy(diis%n_parameters, padded_H(i), x_i, 1, x, 1)
!
      enddo
!
      call mem%dealloc(x_i, diis%n_parameters)
!
!     If CROP is enabled, we replace the current e and x vectors 
!     with the extrapolated e and x 
!
      if (diis%crop) call diis%update_crop_e_and_x(x,             &
                                                   padded_H,      &
                                                   dim_G)
!
      call mem%dealloc(padded_H, dim_G + 1)
!
      diis%iteration = diis%iteration + 1
!
   end subroutine update_diis_tool
!
!
   subroutine update_crop_e_and_x_diis_tool(diis, x, H, dim_G)
!!
!!    Update CROP e and x 
!!    Written by Eirik F. Kjønstad, Nov-Des 2019 
!!
!!    Sets the current e and x to the extrapolated < e > and < x >
!!
!!    x:             extrapolated x, i.e. < x > 
!!    dim_G:   the current dimension of the DIIS matrix
!!    H:             current DIIS vector
!!
      implicit none 
!
      class(diis_tool) :: diis 
!
      real(dp), dimension(diis%n_parameters), intent(in) :: x 
!
      integer, intent(in) :: dim_G
!
      real(dp), dimension(dim_G), intent(in) :: H
!
      real(dp), dimension(:), allocatable :: e_i   ! ith error vector 
      real(dp), dimension(:), allocatable :: e     ! Extrapolated error vector < e >
!
      integer :: i
!
      real(dp) :: ddot
!
!     Overwrite current x with the DIIS extrapolated < x > 
! 
      call diis%x_vectors%set(x, dim_G)
!
!     Compute extrapolated < e > and overwrite current e record with it 
!
      call mem%alloc(e, diis%n_equations)
      call mem%alloc(e_i, diis%n_equations)
!
      call zero_array(e, diis%n_equations)
!
      do i = 1, dim_G
!
!        Read the e_i vector
!
         call diis%e_vectors%get(e_i, i)
!
!        Add w_i e_i to e 
!
         call daxpy(diis%n_equations, H(i), e_i, 1, e, 1)
!
      enddo
!
      call diis%e_vectors%set(e, dim_G) 
!
      if (diis%accumulate) then 
!
!        Recompute the relevant elements of the DIIS matrix 
!
         do i = 1, dim_G
!
            call diis%e_vectors%get(e_i, i)
!
            diis%G(dim_G, i) = ddot(diis%n_equations, e, 1, e_i, 1)
            diis%G(i, dim_G) = diis%G(dim_G, i)
!
         enddo
!
      endif
!
      call mem%dealloc(e, diis%n_equations)
      call mem%dealloc(e_i, diis%n_equations)
!
   end subroutine update_crop_e_and_x_diis_tool
!
!
   function get_dim_G_diis_tool(diis)
!!
!!    Get current dimension 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
      implicit none
!
      class(diis_tool) :: diis
!
      integer :: get_dim_G_diis_tool
!
      if (diis%iteration .gt. diis%dimension_ .and. &
            .not. diis%erase_history) then 
!
         get_dim_G_diis_tool = diis%dimension_
!
      else
!
         get_dim_G_diis_tool = diis%iteration - &
                           diis%dimension_*((diis%iteration-1)/diis%dimension_) 
!
      endif
!
   end function get_dim_G_diis_tool
!
!
   subroutine write_e_and_x_diis_tool(diis, e, x, dim_G)
!!
!!    Write current e and x 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
!!    Stores the current e and x vectors. Each index 1, 2, ..., diis_dim
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
      integer, intent(in) :: dim_G
!
      logical :: cycle_left   
!
      cycle_left = (.not. diis%erase_history) .and. &
                  (dim_G .eq. diis%dimension_) .and. &
                  (dim_G .ne. diis%iteration)
!
      if (cycle_left) then 
!
         call diis%e_vectors%cycle_left() 
         call diis%x_vectors%cycle_left() 
!
      endif 
!
      call diis%e_vectors%set(e, dim_G)
      call diis%x_vectors%set(x, dim_G)
!
   end subroutine write_e_and_x_diis_tool
!
!
   subroutine construct_padded_G_diis_tool(diis, dim_G, padded_G, e)
!!
!!    Construct padded G 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
!!    Constructs the padded DIIS matrix G. For i, j = 1, 2, .., dim_G, we have
!!
!!       G(i,j) = e_i^T e_j.
!!
!!    The padding consists of:
!!
!!       G(dim_G + 1, j) = -1,   j = 1, 2, 3, ..., dim_G 
!!       G(i, dim_G + 1) = -1,   i = 1, 2, 3, ..., dim_G 
!!
!!       G(dim_G + 1, dim_G + 1) = 0.
!!
      implicit none
!
      class(diis_tool) :: diis
!
      integer, intent(in) :: dim_G
!
      real(dp), dimension(dim_G + 1, dim_G + 1), intent(inout) :: padded_G
!
      real(dp), dimension(diis%n_equations), intent(in) :: e
!
      real(dp), dimension(:), allocatable :: e_i  ! To hold previous e_i temporarily
      real(dp), dimension(:), allocatable :: e_j  ! To hold previous e_j temporarily
!
      integer :: i, j
!
      real(dp) :: ddot
!
      padded_G = zero
!
      if (.not. diis%accumulate) then 
!
!        Construct the entire DIIS matrix 
!
         call mem%alloc(e_i, diis%n_equations)
         call mem%alloc(e_j, diis%n_equations)
!
         do i = 1, dim_G
!
            call diis%e_vectors%get(e_i, i)
!
            do j = 1, i
!
               call diis%e_vectors%get(e_j, j)
!
               padded_G(i, j) = ddot(diis%n_equations, e_i, 1, e_j, 1)
               padded_G(j, i) = padded_G(i, j)
!
            enddo
!
         enddo 
!
         call mem%dealloc(e_i, diis%n_equations)
         call mem%dealloc(e_j, diis%n_equations)
!
      else ! accumulate 
!
!        Transfer the parts of the DIIS matrix already calculated 
!
         if (diis%iteration .gt. diis%dimension_ .and. .not. diis%erase_history) then 
!
!           Old matrix is full size; cut out the parts to keep 
!
            padded_G(1 : dim_G - 1, &
                     1 : dim_G - 1) = diis%G(2 : dim_G, &
                                             2 : dim_G)   

!
         elseif (dim_G .gt. 1) then 
!
!           Just copy the old matrix to the new 
!
            padded_G(1 : dim_G - 1, &
                     1 : dim_G - 1) = diis%G(1 : dim_G - 1, &
                                             1 : dim_G - 1)   
!
         endif 
!
!        Compute the new elements of the DIIS matrix 
!
         call mem%alloc(e_i, diis%n_equations)
!
         do i = 1, dim_G
!
            call diis%e_vectors%get(e_i, i)
!
            padded_G(dim_G, i) = ddot(diis%n_equations, e, 1, e_i, 1)
            padded_G(i, dim_G) = padded_G(dim_G, i)
!
         enddo
!
         call mem%dealloc(e_i, diis%n_equations)
!
!        Keep a copy of the DIIS matrix for the next iteration 
!
         diis%G(1 : dim_G, 1 : dim_G) = padded_G(1 : dim_G, 1 : dim_G)
!
      endif
!
!     Add the (-1)-padding in the last row and column  
!
      do i = 1, dim_G 
!
         padded_G(dim_G + 1, i) = -one
         padded_G(i, dim_G + 1) = -one
!
      enddo
!
   end subroutine construct_padded_G_diis_tool
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
      integer :: dim_G
!
      dim_G = diis%get_dim_G()
!
      if (i .gt. dim_G) call output%error_msg('Asked for an x not in use.')
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
      integer :: dim_G
!
      dim_G = diis%get_dim_G()
!
      if (i .gt. dim_G) call output%error_msg('Asked for an e not in use.')
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
      integer :: dim_G
!
      dim_G = diis%get_dim_G()
!
      if (diis%accumulate) call output%error_msg('Can not set x for cumulative diis.')
      if (i .gt. dim_G) call output%error_msg('Asked for an x not in use.')
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
      integer :: dim_G
!
      dim_G = diis%get_dim_G()
!
      if (diis%accumulate) call output%error_msg('Can not set e for cumulative diis.')
      if (i .gt. dim_G) call output%error_msg('Asked for an e not in use.')
!
      call diis%e_vectors%set(e_i, i)
!
   end subroutine write_e_diis_tool
!
!
end module diis_tool_class
