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
      procedure :: alloc_2_memory_manager
      generic   :: alloc       => alloc_2_memory_manager
!
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
   subroutine alloc_2_memory_manager(mem, array, M, N)
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
   end subroutine alloc_2_memory_manager
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
end module memory_manager_class
