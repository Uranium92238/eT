module memory_manager_class
!
!!
!!    Memory manager class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Dec 2017
!!
!!    The memory manager class handles the memory used by the model calculation,
!!    and there is an object called 'mem' in the wavefunction object of this class.
!!
!!    To account for the available memory, all large arrays must be allocated and deallocated
!!    using the memory manager. A typical usage of the 'mem' object is as follows:
!!
!!       real(dp), dimension(:,:), allocatable :: array -> declares an allocatable array
!!
!!       call wf%mem%alloc(array, M, N)   -> allocates the array of dimension M x N
!!
!!       ... Do stuff with the array
!!
!!       call wf%mem%dealloc(array, M, N) -> deallocates the array of dimension M x N
!!
!!
!!    NB! Large arrays MUST always be allocated using the memory manager object. Small arrays,
!!    integers, strings, etc., which use a negligible amount of memory, do not need to pass through
!!    the memory manager.
!!
!!
!!    The 'alloc' and 'dealloc' routines allow the memory manager keep track of the
!!    the memory available at a given time. From the specified total memory, the class
!!    can then set the batching information in a batching index (see the batching
!!    index class). See the num_batch and num_two_batches procedures for more details.
!!
!
   use kinds
   use parameters
   use file_class
   use batching_index_class
   use io_utilities
!
!  ::::::::::::::::::::::::::::::::::::::::::::::::
!  -::- Definition of the memory_manager class -::-
!  ::::::::::::::::::::::::::::::::::::::::::::::::
!
   type :: memory_manager
!
!     The total amount of memory specified by user (standard: 8 GB)
!
      integer(i15) :: total = 8000000000
!
!     The amount of memory currently available, based on the arrays currently allocated
!     (memory used by objects and local variables are not included in this estimate)
!
      integer(i15) :: available = 8000000000
!
!     Buffer for handling batches (standard: 10%). This means in practice that 'required
!     memory' estimates are increased by 10% in case they miss they slightly underestimate
!     the correct memory requirements
!
      integer(i15) :: buffer = 10 ! 10%
!
   contains
!
!     Initialization routine (used if user specifies a memory different from standard)
!
      procedure :: prepare => prepare_memory_manager
!
!     Allocation and deallocation routines for double precision arrays
!
      procedure :: alloc       => alloc_memory_manager
      procedure :: dealloc     => dealloc_memory_manager
!
!     Allocation and deallocation routines for integer arrays
!
      procedure :: alloc_int   => alloc_int_memory_manager
      procedure :: dealloc_int => dealloc_int_memory_manager
!
!     Routines for determining the number of batches
!
      procedure :: num_batch     => num_batch_memory_manager     ! For one-index batch
      procedure :: num_two_batch => num_two_batch_memory_manager ! For two-index batches
!
      procedure :: read_settings  => read_settings_memory_manager
      procedure :: print_settings => print_settings_memory_manager
!
!     Routines that ask for information
!
      procedure :: room_for_n_arrays_of_size => room_for_n_arrays_of_size_memory_manager
!
      procedure :: batch_setup_1_memory_manager
      procedure :: batch_setup_2_memory_manager
      procedure :: batch_setup_3_memory_manager
      generic   :: batch_setup => batch_setup_1_memory_manager, batch_setup_2_memory_manager, &
                                  batch_setup_3_memory_manager
!
   end type memory_manager
!
!  Main memory object
!
   type(memory_manager) :: mem
!
contains
!
!
   subroutine prepare_memory_manager(mem)
!!
!!    Prepare (memory manager)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Dec 2017
!!
!!    Prepares the memory manager object by setting the
!!    total and initial available memory. 
!!
      implicit none
!
      class(memory_manager) :: mem
!
      if (requested_section('memory')) then
!
         call mem%read_settings()
!
      else
!
!        Set default value 
!
         mem%total = 8*1000000000
!
      endif
!
      mem%available = mem%total
!
      call mem%print_settings()
!
   end subroutine prepare_memory_manager
!
!
   subroutine alloc_memory_manager(mem, array, M, N)
!!
!!    Alloc (memory manager)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Dec 2017
!!
!!    Allocates a double precision array and updates the available
!!    memory accordingly.
!!
      implicit none
!
      class(memory_manager) :: mem
!
      real(dp), dimension(:,:), allocatable :: array
!
      integer(i15), intent(in) :: M, N ! First and second dimension of array that is being allocated
!
      integer(i15) :: size_array ! Total size of array (M*N)
      integer(i15) :: stat = 0
      integer(i15) :: error = 0
!
      size_array = M*N
!
!     Allocate array and check whether allocation was successful
!
      allocate(array(M,N), stat = error)
!
      if (stat .ne. 0) then
!
         call output%error_msg('could not allocate array with #elements =', size_array)
!
      endif
!
!     Update the available memory
!
!     The 'double precision' type (see types.F90) is typically 8 bytes,
!     though it might differ due to its definition in terms of precision.
!
      mem%available = mem%available - dp*size_array
!
!     Check if there is no more memory (defined as being no more memory
!     left of what was specified by user as available)
!
      if (mem%available .lt. 0) then
!
         call output%error_msg('user-specified memory insufficient.')
!
      endif
!
   end subroutine alloc_memory_manager
!
!
   subroutine dealloc_memory_manager(mem, array, M, N)
!!
!!    Dealloc (memory manager)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Dec 2017
!!
!!    Deallocates a double precision array and updates the available
!!    memory accordingly.
!!
      implicit none
!
      class(memory_manager) :: mem
!
      real(dp), dimension(:,:), allocatable :: array
!
      integer(i15), intent(in) :: M, N ! First and second dimension of array that is being allocated
!
      integer(i15) :: size_array ! Total size of array (M*N)
      integer(i15) :: stat = 0
      integer(i15) :: error = 0
!
      size_array = M*N
!
!     Deallocate array and check whether deallocation was successful
!
      deallocate(array, stat = error)
!
      if (stat .ne. 0) then
!
         call output%error_msg('could not deallocate array with #elements =', size_array)
!
      endif
!
!     Update the available memory
!
!     The 'double precision' type (see types.F90) is typically 8 bytes,
!     though it might differ due to its definition in terms of precision.
!
      mem%available = mem%available + dp*size_array
!
   end subroutine dealloc_memory_manager
!
!
   subroutine alloc_int_memory_manager(mem, array, M, N)
!!
!!    Alloc int (memory manager)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Dec 2017
!!
!!    Allocates an integer array and updates the available
!!    memory accordingly.
!!
      implicit none
!
      class(memory_manager) :: mem
!
      integer(i15), dimension(:,:), allocatable :: array
!
      integer(i15), intent(in) :: M, N ! First and second dimension of array that is being allocated
!
      integer(i15) :: size_array ! Total size of array (M*N)
      integer(i15) :: stat = 0
      integer(i15) :: error = 0
!
      size_array = M*N
!
!     Allocate array and check whether allocation was successful
!
      allocate(array(M,N), stat = error)
!
      if (stat .ne. 0) then
!
         call output%error_msg('Error: could not allocate array with #elements =', size_array)
!
      endif
!
!     Update the available memory
!
!     The 'integer 15', or i15, type (see types.F90) is typically 8 bytes,
!     though it might differ due to its definition in terms of precision.
!
      mem%available = mem%available - i15*size_array
!
!     Check if there is no more memory (defined as being no more memory
!     left of what was specified by user as available)
!
      if (mem%available .lt. 0) then
!
         call output%error_msg('user-specified memory insufficient.')
!
      endif
!
   end subroutine alloc_int_memory_manager
!
!
   subroutine dealloc_int_memory_manager(mem, array, M, N)
!!
!!    Dealloc int (memory manager)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Dec 2017
!!
!!    Deallocates an integer array and updates the available
!!    memory accordingly.
!!
      implicit none
!
      class(memory_manager) :: mem
!
      integer(i15), dimension(:,:), allocatable :: array
!
      integer(i15), intent(in) :: M, N ! First and second dimension of array that is being allocated
!
      integer(i15) :: size_array ! Total size of array (M*N)
      integer(i15) :: stat = 0
      integer(i15) :: error = 0
!
      size_array = M*N
!
!     Deallocate array and check whether deallocation was successful
!
      deallocate(array, stat = error)
!
      if (stat .ne. 0) then
!
         call output%error_msg('could not deallocate array with #elements =', size_array)
!
      endif
!
!     Update the available memory
!
!     The 'integer 15', or i15, type (see types.F90) is typically 4 bytes,
!     though it might differ due to its definition in terms of precision.
!
      mem%available = mem%available + i15*size_array
!
   end subroutine dealloc_int_memory_manager
!
!
   subroutine num_batch_memory_manager(mem, batch_p, required)
!!
!!    Number of batches
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Dec 2017
!!
!!    Given the required memory, this routine determines, for a one-index batching,
!!    the maximum batch length and the number of batches in total.
!!
      implicit none
!
      class(memory_manager) :: mem
!
      class(batching_index) :: batch_p ! The index being batched over
!
      integer(i15) :: required
!
!     Add buffer to required estimate
!
      required = required + required/(mem%buffer)
!
      if (required .lt. mem%available) then
!
!        No need to batch
!
         batch_p%num_batches = 1
         batch_p%max_length = batch_p%index_dimension
!
         return
!
      endif
!
!     We need to batch
!
!     Determine maximum batch length
!
      batch_p%max_length = (mem%available)/(required/(batch_p%index_dimension))
!
!     Number of full batches
!
      batch_p%num_batches = (batch_p%index_dimension)/(batch_p%max_length)
!
!     Test for rest not included in the preceding integer division
!
      if ((batch_p%num_batches)*(batch_p%max_length) .lt. batch_p%index_dimension) then
!
         batch_p%num_batches = batch_p%num_batches + 1
!
      endif
!
   end subroutine num_batch_memory_manager
!
!
   subroutine num_two_batch_memory_manager(mem, batch_p, batch_q, required)
!!
!!    Number of two-batches
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Dec 2017
!!
!!    Given the required memory, this routine determines, for two-index batching,
!!    the maximum batch length and the number of batches in total of both batching
!!    indices.
!!
      implicit none
!
      class(memory_manager) :: mem
!
      class(batching_index) :: batch_p
      class(batching_index) :: batch_q
!
      integer(i15) :: required
!
      integer(i15) :: i = 0
!
      required = required + required/(mem%buffer)
!
      batch_p%num_batches = 1
      batch_q%num_batches = 1
!
      if (required .lt. mem%available) then
!
!        No need to batch
!
         batch_p%num_batches = 1
         batch_q%num_batches = 1
!
         batch_p%max_length = batch_p%index_dimension
         batch_q%max_length = batch_q%index_dimension
!
         return
!
      endif
!
!     We need to batch
!
!     Test whether two different-length indices are requested,
!     because this feature is not yet implemented
!
      if (batch_p%index_dimension .ne. batch_q%index_dimension) then
!
         call output%error_msg('batching over indices of different lengths is not yet implemented')
!
      endif
!
!     Determine number of batches
!
      do i = 1, batch_p%index_dimension
!
         if (mem%available .gt. required/i**2) then
!
            batch_p%num_batches = i
            batch_q%num_batches = i
!
            batch_p%max_length = (batch_p%index_dimension)/(batch_p%num_batches)
            batch_q%max_length = (batch_q%index_dimension)/(batch_q%num_batches)
!
!           Test for rest
!
            if ((batch_p%num_batches)*(batch_p%max_length) .lt. batch_p%index_dimension) then
!
               batch_p%num_batches = batch_p%num_batches + 1
               batch_p%num_batches = batch_p%num_batches + 1
!
            endif
!
            return
!
         endif
!
      enddo
!
   end subroutine num_two_batch_memory_manager
!
!
   subroutine read_settings_memory_manager(mem)
!!
!!    Read settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!  
      class(memory_manager) :: mem
!
      integer(i15) :: n_specs, i
!
      character(len=100) :: line
!
      mem%total = 8.0d0
!
      call move_to_section('memory', n_specs)
!
      do i = 1, n_specs
!
         read(input%unit, '(a100)') line
         line = remove_preceding_blanks(line)
!
         if (line(1:10) == 'available:' ) then
!
            read(line(11:100), *) mem%total
!
         endif
!
      enddo
!
      mem%total = mem%total*1000000000
!
   end subroutine read_settings_memory_manager
!
!
   subroutine print_settings_memory_manager(mem)
!!
!!    Print settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!  
      class(memory_manager) :: mem
!
      write(output%unit, '(t3, a38, i5, a)') 'Memory available for calculation:     ', mem%total/1000000000, ' GB'
!
   end subroutine print_settings_memory_manager
!
!
   integer(i15) function room_for_n_arrays_of_size_memory_manager(mem, M)
!!
!!    Room for number of arrays of size
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
!!    Returns the number of double precision arrays of size M that can currently be held 
!!    in memory, based on what is - according to the manager - available.
!!
      implicit none 
!
      class(memory_manager), intent(in) :: mem 
!
      integer(i15), intent(in) :: M
!
      room_for_n_arrays_of_size_memory_manager = (mem%available)/(dp*M)
!
   end function room_for_n_arrays_of_size_memory_manager
!
!
   subroutine batch_setup_1_memory_manager(mem, batch_p, req0, req1, element_size)
!!
!!    Setup batching 
!!    Written by Rolf H. Myhre and Eirik F. Kjønstad, December 2018
!!
!!    Batching setup for a single index.
!!
!!    batch_p:  Initialized batching object.
!!
!!    req0:     Memory required that does not change with the index dimension.
!!              E.g., n_o**2*n_v**2 for (vo|vo) if none of the indices 
!!              in the integral is batched over.
!!
!!    req1:     Memory required per batching index (linear with batch size). 
!!              E.g., n_v**3 for (vv|vo) when batching over the 
!!              occupied index.
!!
      implicit none
!
      class(memory_manager) :: mem
!
      class(batching_index) :: batch_p ! The index being batched over
!
      integer(i15), intent(in) :: req0
      integer(i15), intent(in) :: req1
!
      integer(i15), intent(in), optional :: element_size
!
      integer(i15) :: req0_tot
      integer(i15) :: req1_min
      integer(i15) :: req_min
!
      integer(i15) :: req_tot
!
      integer(i15) :: e_size
!
      e_size = dp
      if(present(element_size)) then
         e_size = element_size
      endif
!
      req0_tot = (req0 + req0/(mem%buffer))*e_size
      req1_min = (req1 + req1/(mem%buffer))*e_size 
!
      req_min = req0_tot + req1_min 
      req_tot = req0_tot + req1_min*batch_p%index_dimension
!
      if (req_tot .lt. mem%available) then
!
!        No need to batch
!
         batch_p%num_batches = 1
         batch_p%max_length  = batch_p%index_dimension
!
      else if (req_min .gt. mem%available) then
!
!        Not enough memory for a batch
!
         write(output%unit,'(t3,a,i14,a,i14)') 'Need at least', req_min, 'but only have ', mem%available
         call output%error_msg('Not enough memory for a batch.')
!
      else
!
!        We need to batch
!
!        Determine maximum batch length
!
         batch_p%max_length = (mem%available - req0_tot)/(req1_min)
!
!        Number of full batches
!
         batch_p%num_batches = (batch_p%index_dimension-1)/(batch_p%max_length)+1
!
      endif
!
   end subroutine batch_setup_1_memory_manager
!
!
   subroutine batch_setup_2_memory_manager(mem, batch_p, batch_q, req0, req1_p, req1_q, req2, element_size)
!!
!!    Setup batching 
!!    Written by Rolf H. Myhre and Eirik F. Kjønstad, Dec 2018
!!
!!    Batching setup for two batching indices.
!!
!!    batch_p: Initialized batching object
!!    batch_q: Initialized batching object
!!
!!    req0: required memory that does not scale with batch size 
!! 
!!    req1: required memory that scales linearly with p batch size
!!    req1: required memory that scales linearly with q batch size
!!
!!    req2: required memory that scales quadratically with batch size
!!
!!    element_size: memory per element, default is double precision
!!
!!    If you are batching over i and j and need to keep g_abij, g_abci and g_abcj in memory, 
!!    req1 = 2*n_v**3 (both for i and j) and req2 = n_v**2. Memory per batch is then 
!!    batch_size*req1 + batch_size**2*req2
!!
!!    If you are batching over a and j and need to keep g_abij, g_abci and g_abcj in memory,
!!    req1_a = n_o*n_v**2, req1_j = 0, and req2 = n_o*n_v + n_v**2. Note that one integral (g_abci)
!!    scales linearly with the a-index but that there are no such integrals for the j-index.
!!
!!    Be careful with symmetries and permutations!
!!
      implicit none
!
      class(memory_manager) :: mem
!
      class(batching_index) :: batch_p ! An index being batched over
      class(batching_index) :: batch_q ! An index being batched over
!
      integer(i15), intent(in) :: req0
      integer(i15), intent(in) :: req1_p
      integer(i15), intent(in) :: req1_q
      integer(i15), intent(in) :: req2
!
      integer(i15), intent(in), optional :: element_size
!
      logical :: figgered_out
!
      integer(i15) :: req0_tot
      integer(i15) :: req1_p_min
      integer(i15) :: req1_q_min 
      integer(i15) :: req2_min
      integer(i15) :: req_min
      integer(i15) :: req_tot 
!
      integer(i15) :: p_elements, q_elements
!
      integer(i15) :: e_size
!
      e_size = dp
      if(present(element_size)) then
         e_size = element_size
      endif
!
      req0_tot   = (req0 + req0/(mem%buffer))*e_size 
      req1_p_min = (req1_p + req1_p/(mem%buffer))*e_size
      req1_q_min = (req1_q + req1_q/(mem%buffer))*e_size
      req2_min = (req2 + req2/(mem%buffer))*e_size
!
      req_min = req0_tot + req1_p_min + req1_q_min + req2_min 
!
      req_tot = req0_tot + req1_p_min*(batch_p%index_dimension) &
                         + req1_q_min*(batch_q%index_dimension) &
                         + req2_min*(batch_p%index_dimension)*(batch_q%index_dimension)
!
      if (req_tot .lt. mem%available) then
!
!        No need to batch
!
         batch_p%num_batches = 1
         batch_p%max_length  = batch_p%index_dimension
!
         batch_q%num_batches = 1
         batch_q%max_length  = batch_q%index_dimension
!
      else if (req_min .gt. mem%available) then
!
!        Not enough memory for a batch
!
         write(output%unit,'(t3,a,i14,a,i14)') 'Need ', req_min, 'but only have ', mem%available
         call output%error_msg('Not enough memory for a batch.')
!
      else
!
!        We need to batch
!
!        Figure out how many we have room for
!
!        I. First, try to increment both indices simultaneously
!
         p_elements = 1
         q_elements = 1
!
         figgered_out = .false.
         do while (.not. figgered_out                              &
                     .and. p_elements .lt. batch_p%index_dimension &
                     .and. q_elements .lt. batch_q%index_dimension)
!
            if ((p_elements+1)*(q_elements+1)*req2_min &
                  + (p_elements+1)*req1_p_min          &
                  + (q_elements+1)*req1_q_min          &
                  + req0 .lt. mem%available) then 
!
               p_elements = p_elements + 1 ! can hold +1 batch size 
               q_elements = q_elements + 1
!
            else
!
               figgered_out = .true.       ! cannot hold +1 batch size 
!
            endif
!
         enddo
!
!        II. If simultaneous incrementation was not sufficient,
!            then try to increment the largest index further. This is 
!            guaranteed to work, so let's just go ahead and increment
!            with no safeguards in place.
! 
         if (.not. figgered_out) then 
!
            if (batch_p%index_dimension .gt. batch_q%index_dimension) then 
!
!              Increment p 
!
               do while (((p_elements+1)*q_elements*req2_min &
                           + (p_elements+1)*req1_p_min       &
                           + q_elements*req1_q_min           &
                           + req0) .lt. mem%available)
!
                  p_elements = p_elements + 1
!
               enddo
!
            elseif (batch_p%index_dimension .lt. batch_q%index_dimension) then 
!
!              Increment q 
!
               do while ((p_elements*(q_elements+1)*req2_min &
                           + p_elements*req1_p_min           &
                           + (q_elements+1)*req1_q_min       &
                           + req0) .lt. mem%available)
!
                  q_elements = q_elements + 1
!
               enddo               
!
            else
!
               call output%error_msg('Something went very wrong! Expected different-sized' // &
                                      'indices, but got same-sized indices (in batching setup).')
!
            endif
!
            figgered_out = .true. 
!
         endif
!
         batch_p%max_length = p_elements         
         batch_q%max_length = q_elements         
!
!        Figure out how many batches
!
         batch_p%num_batches = (batch_p%index_dimension-1)/(batch_p%max_length)+1
         batch_q%num_batches = (batch_q%index_dimension-1)/(batch_q%max_length)+1
!
      endif
!
   end subroutine batch_setup_2_memory_manager
!
!
   subroutine batch_setup_3_memory_manager(mem, batch_p, batch_q, batch_r, req0, req1_p, req1_q, &
                                       req1_r, req2_pq, req2_pr, req2_qr, req3, element_size)
!!
!!    Setup batching 
!!    This is setup for two batch indices
!!    Written by Rolf H. Myhre December 2018
!!
!!    Batching setup for three batching indices.
!!
!!    batch_p: Initialized batching object
!!    batch_q: Initialized batching object
!!    batch_r: Initialized batching object
!!
!!    req0: required memory that does not scale with batch size 
!! 
!!    req1_p: required memory that scales linearly with p batch size
!!    req1_q: required memory that scales linearly with q batch size
!!    req1_r: required memory that scales linearly with r batch size
!!
!!    req2_pq: required memory that scales quadratically with pq batch size
!!    req2_pr: required memory that scales quadratically with pr batch size
!!    req2_qr: required memory that scales quadratically with qr batch size
!!
!!    req3: required memory that scales cubically with batch indices pqr
!!
!!    element_size: memory per element, default is double precision
!!
!!    Be careful with symmetries and permutations!
!!
      implicit none
!
      class(memory_manager) :: mem
!
      class(batching_index) :: batch_p ! An index being batched over
      class(batching_index) :: batch_q ! An index being batched over
      class(batching_index) :: batch_r ! An index being batched over
!
      integer(i15), intent(in) :: req0
! 
      integer(i15), intent(in) :: req1_p 
      integer(i15), intent(in) :: req1_q 
      integer(i15), intent(in) :: req1_r 
!
      integer(i15), intent(in) :: req2_pq
      integer(i15), intent(in) :: req2_pr
      integer(i15), intent(in) :: req2_qr 
!
      integer(i15), intent(in) :: req3
!
      integer(i15), intent(in), optional :: element_size
!
      integer(i15) :: req0_tot
!
      integer(i15) :: req1_p_min
      integer(i15) :: req1_q_min 
      integer(i15) :: req1_r_min
! 
      integer(i15) :: req2_pq_min
      integer(i15) :: req2_pr_min
      integer(i15) :: req2_qr_min
!
      integer(i15) :: req3_min
!
      integer(i15) :: req_min
      integer(i15) :: req_tot 
!
      integer(i15) :: p_elements, q_elements, r_elements
!
      logical :: found_batch_size, p_incremented, q_incremented, r_incremented
!
      integer(i15) :: e_size
!
      e_size = dp
      if(present(element_size)) then
         e_size = element_size
      endif
!
      req0_tot   = (req0 + req0/(mem%buffer))*e_size
! 
      req1_p_min = (req1_p + req1_p/(mem%buffer))*e_size
      req1_q_min = (req1_q + req1_q/(mem%buffer))*e_size
      req1_r_min = (req1_r + req1_r/(mem%buffer))*e_size
!
      req2_pq_min = (req2_pq + req2_pq/(mem%buffer))*e_size
      req2_pr_min = (req2_pr + req2_pr/(mem%buffer))*e_size
      req2_qr_min = (req2_qr + req2_qr/(mem%buffer))*e_size
!
      req3_min = (req3 + req3/(mem%buffer))*e_size
!
      req_min = req0_tot + req1_p_min + req1_q_min + req1_r_min &
                           + req2_pq_min + req2_pr_min + req2_qr_min + req3_min  
!
      req_tot = req0_tot + req1_p_min*(batch_p%index_dimension) &
                         + req1_q_min*(batch_q%index_dimension) &
                         + req1_r_min*(batch_r%index_dimension) &
                         + req2_pq_min*(batch_p%index_dimension)*(batch_q%index_dimension) &
                         + req2_pr_min*(batch_p%index_dimension)*(batch_r%index_dimension) &
                         + req2_qr_min*(batch_q%index_dimension)*(batch_r%index_dimension) &
                         + req3_min*(batch_p%index_dimension)*(batch_q%index_dimension)*(batch_r%index_dimension)
!
      if (req_tot .lt. mem%available) then
!
!        No need to batch
!
         batch_p%num_batches = 1
         batch_p%max_length  = batch_p%index_dimension
!
         batch_q%num_batches = 1
         batch_q%max_length  = batch_q%index_dimension
!
         batch_r%num_batches = 1
         batch_r%max_length  = batch_r%index_dimension
!
      else if (req_min .gt. mem%available) then
!
!        Not enough memory for a batch
!
         write(output%unit,'(t3,a,i14,a,i14)') 'Need ', (req_min), 'but only have ', &
                                               mem%available
         call output%error_msg('Not enough memory for a batch')
!
      else
!
!        First, try to increment both indices simultaneously
!
         p_elements = 1
         q_elements = 1
         r_elements = 1

!
         found_batch_size = .false.
         p_incremented = .true.
         q_incremented = .true.
         r_incremented = .true.
!
         do while (.not. found_batch_size .and. (p_incremented .or. q_incremented .or. r_incremented))
!
            if ((p_elements) .lt. batch_p%index_dimension) then
               p_elements = p_elements + 1
               p_incremented = .true.
            else
               p_incremented = .false.
            endif
!
            if ((q_elements) .lt. batch_q%index_dimension) then
               q_elements = q_elements + 1
               q_incremented = .true.
            else
               q_incremented = .false.
            endif
!
            if ((r_elements) .lt. batch_r%index_dimension) then
               r_elements = r_elements + 1
               r_incremented = .true.
            else
               r_incremented = .false.
            endif
!
            if ( (p_elements)*(q_elements)*(r_elements)*req3_min &
                  + (p_elements)*(q_elements)*req2_pq &
                  + (p_elements)*(r_elements)*req2_pr &
                  + (q_elements)*(r_elements)*req2_qr &
                  + (p_elements)*req1_p_min          &
                  + (q_elements)*req1_q_min          &
                  + (r_elements)*req1_r_min          &
                  + req0 .ge. mem%available) then 
!
                  found_batch_size = .true.       ! cannot hold +1 batch size 
                  if (p_incremented) p_elements = p_elements - 1
                  if (q_incremented) q_elements = q_elements - 1
                  if (r_incremented) r_elements = r_elements - 1
!
            endif
!
         enddo
!
         batch_p%max_length = p_elements         
         batch_q%max_length = q_elements        
         batch_r%max_length = r_elements         
!
!        Figure out how many batches
!
         batch_p%num_batches = (batch_p%index_dimension-1)/(batch_p%max_length)+1
         batch_q%num_batches = (batch_q%index_dimension-1)/(batch_q%max_length)+1
         batch_r%num_batches = (batch_r%index_dimension-1)/(batch_r%max_length)+1
!
      endif
!
   end subroutine batch_setup_3_memory_manager
!
!
end module memory_manager_class
