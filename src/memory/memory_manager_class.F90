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
!!       call mem%alloc(array, M, N)   -> allocates the array of dimension M x N
!!
!!       ... Do stuff with the array
!!
!!       call mem%dealloc(array, M, N) -> deallocates the array of dimension M x N
!!
!!    Analogous calls are made to make one, three and four dimensional tensors as well,
!!    e.g. call mem%alloc(X, M, N, K) for an M x N x K tensor.
!!
!!    Note: Large arrays MUST always be allocated using the memory manager object. Small arrays,
!!    integers, strings, etc., which use a negligible amount of memory, do not have to pass 
!!    through the memory manager.
!!
!!    The 'alloc' and 'dealloc' routines allow the memory manager keep track of the
!!    the memory available at a given time. From the specified total memory, the class
!!    can then set the batching information for a batching index (see the batching
!!    index class) or set of such indices (see the batch_setup routines below).
!!
!
   use kinds
   use parameters
   use global_out, only : output
   use global_in, only : input
   use batching_index_class, only : batching_index
   use io_utilities
!
!     Debug option:
!     Require that batch setup always gives batching 
!
#ifdef _FORCED_BATCHING
      logical, parameter, private :: force_batch = .true.
#else
      logical, parameter, private :: force_batch = .false.    
#endif 
!
!  Class definition 
!
   type :: memory_manager
!
!     The total amount of memory specified by user (standard: 8 GB)
!
      integer(i15) :: total
!
!     The amount of memory currently available, based on the arrays currently allocated
!     (memory used by objects and local variables are not included in this estimate)
!
      integer(i15) :: available
!
   contains
!
!     Initialization routine (used if user specifies a memory different from standard)
!
      procedure :: prepare          => prepare_memory_manager
!
      procedure :: check_for_leak   => check_for_leak_memory_manager
!
!     Allocation and deallocation routines for arrays
!
      procedure :: alloc_1_memory_manager
      procedure :: alloc_2_memory_manager
      procedure :: alloc_3_memory_manager
      procedure :: alloc_4_memory_manager
      procedure :: alloc_5_memory_manager
      procedure :: alloc_int_1_memory_manager
      procedure :: alloc_int_2_memory_manager
      procedure :: alloc_int_3_memory_manager
      procedure :: alloc_int_4_memory_manager
      generic   :: alloc => alloc_1_memory_manager, &
                            alloc_2_memory_manager, &
                            alloc_3_memory_manager, &
                            alloc_4_memory_manager, &
                            alloc_5_memory_manager, &
                            alloc_int_1_memory_manager, &
                            alloc_int_2_memory_manager, &
                            alloc_int_3_memory_manager, &
                            alloc_int_4_memory_manager
!
      procedure :: dealloc_1_memory_manager
      procedure :: dealloc_2_memory_manager
      procedure :: dealloc_3_memory_manager
      procedure :: dealloc_4_memory_manager
      procedure :: dealloc_5_memory_manager
      procedure :: dealloc_int_1_memory_manager
      procedure :: dealloc_int_2_memory_manager
      procedure :: dealloc_int_3_memory_manager
      procedure :: dealloc_int_4_memory_manager
      generic   :: dealloc => dealloc_1_memory_manager, &
                              dealloc_2_memory_manager, &
                              dealloc_3_memory_manager, &
                              dealloc_4_memory_manager, &
                              dealloc_5_memory_manager, &
                              dealloc_int_1_memory_manager, &
                              dealloc_int_2_memory_manager, &
                              dealloc_int_3_memory_manager, &
                              dealloc_int_4_memory_manager
!
!     Routines for determining the number of batches
!
      procedure :: batch_setup_1_memory_manager
      procedure :: batch_setup_2_memory_manager
      procedure :: batch_setup_3_memory_manager
      generic   :: batch_setup         => batch_setup_1_memory_manager, batch_setup_2_memory_manager, &
                                          batch_setup_3_memory_manager
!
      procedure :: batch_setup_3_ident_memory_manager
      generic   :: batch_setup_ident   => batch_setup_3_ident_memory_manager
!
      procedure :: read_settings       => read_settings_memory_manager
      procedure :: print_settings      => print_settings_memory_manager
!
      procedure :: get_available       => get_available_memory_manager
      procedure :: print_available     => print_available_memory_manager
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
!     Set standard and read settings 
!
      mem%total = 8 ! GB
!
      call mem%read_settings()
!
!     Convert from GB to B
!
      mem%total = mem%total*1000000000
!
      mem%available = mem%total
!
      call mem%print_settings()
!
   end subroutine prepare_memory_manager
!
!
   subroutine check_for_leak_memory_manager(mem)
!!
!!    Check for leak 
!!    Written by Eirik F. Kjønstad, Apr 2019 
!!
!!    Issues a warning if there has been a leak since the 
!!    the memory manager was prepared. Should only be called 
!!    at the end of the program, when all arrays that were
!!    allocated since mem%prepare() should have been deallocated.
!!
      implicit none 
!
      class(memory_manager), intent(in) :: mem 
!
      if (mem%available .ne. mem%total) then 
!
         write(output%unit, '(/t3,a)')  'Mismatch in memory according to eT and specified on input:'
!
         write(output%unit, '(/t6, a27, f11.4, a)') 'Memory available (eT):     ', &
                         real(mem%available)/real(1000000000), ' GB'
         write(output%unit, '(t6, a27, f11.4, a)') 'Memory available (input):   ', &
                         real(mem%total)/real(1000000000), ' GB'
!
         call output%warning_msg('Deallocations may be missing or specified with incorrect dimensionalities.')
!
      endif 
!
   end subroutine check_for_leak_memory_manager
!
!
   integer(i15) function get_available_memory_manager(mem)
!!
!!    Get available  
!!    Written by Eirik F. Kjønstad, Jan 2019 
!!
      implicit none 
!
      class(memory_manager), intent(in) :: mem 
!
      get_available_memory_manager = mem%available
!
   end function get_available_memory_manager
!
!
   subroutine print_available_memory_manager(mem)
!!
!!    Get available  
!!    Written by Eirik F. Kjønstad, Jan 2019 
!!
      implicit none 
!
      class(memory_manager), intent(in) :: mem 
!
      write(output%unit, '(t3, a38, i5, a)') 'Currently available memory:     ', &
                         mem%available/1000000, ' MB'
!
   end subroutine print_available_memory_manager
!
!
   subroutine alloc_1_memory_manager(mem, array, M)
!!
!!    Alloc (memory manager)
!!    Written by Rolf H. Myhre, January 2019
!!
!!    Allocates a one dimensional double precision array and updates the available
!!    memory accordingly.
!!
      implicit none
!
      class(memory_manager) :: mem
!
      real(dp), dimension(:), allocatable :: array
!
      integer, intent(in) :: M ! Dimension of array that is being allocated
!
      integer :: size_array ! Total size of array (M)
      integer :: error = 0
!
      size_array = M
!
!     Allocate array and check whether allocation was successful
!
      allocate(array(M), stat = error)
!
      if (error .ne. 0) then
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
   end subroutine alloc_1_memory_manager
!
!
   subroutine alloc_2_memory_manager(mem, array, M, N)
!!
!!    Alloc (memory manager)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Dec 2017
!!
!!    Allocates a two dimensional double precision array and updates the available
!!    memory accordingly.
!!
      implicit none
!
      class(memory_manager) :: mem
!
      real(dp), dimension(:,:), allocatable :: array
!
      integer, intent(in) :: M, N ! First and second dimension of array that is being allocated
!
      integer :: size_array ! Total size of array (M*N)
      integer :: error = 0
!
      size_array = M*N
!
!     Allocate array and check whether allocation was successful
!
      allocate(array(M,N), stat = error)
!
      if (error .ne. 0) then
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
   subroutine alloc_3_memory_manager(mem, array, M, N, O)
!!
!!    Alloc (memory manager)
!!    Written by Rolf H. Myhre, January 2019
!!
!!    Allocates a three dimensional double precision array and updates the available
!!    memory accordingly.
!!
      implicit none
!
      class(memory_manager) :: mem
!
      real(dp), dimension(:,:,:), allocatable :: array
!
      integer, intent(in) :: M, N, O ! First, second and third dimension of array 
!
      integer :: size_array ! Total size of array (M*N*O)
      integer :: error = 0
!
      size_array = M*N*O
!
!     Allocate array and check whether allocation was successful
!
      allocate(array(M,N,O), stat = error)
!
      if (error .ne. 0) then
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
   end subroutine alloc_3_memory_manager
!
!
   subroutine alloc_4_memory_manager(mem, array, M, N, O, P)
!!
!!    Alloc (memory manager)
!!    Written by Rolf H. Myhre, January 2019
!!
!!    Allocates a four dimensional double precision array and updates the available
!!    memory accordingly.
!!
      implicit none
!
      class(memory_manager) :: mem
!
      real(dp), dimension(:,:,:,:), allocatable :: array
!
      integer, intent(in) :: M, N, O, P ! First, second, third and fourth dimension of array 
!
      integer :: size_array ! Total size of array (M*N*O*P)
      integer :: error = 0
!
      size_array = M*N*O*P
!
!     Allocate array and check whether allocation was successful
!
      allocate(array(M,N,O,P), stat = error)
!
      if (error .ne. 0) then
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
   end subroutine alloc_4_memory_manager
!
!
   subroutine alloc_5_memory_manager(mem, array, M, N, O, P, Q)
!!
!!    Alloc (memory manager)
!!    Written by Rolf H. Myhre, January 2019
!!
!!    Allocates a five dimensional double precision array and updates the available
!!    memory accordingly.
!!
      implicit none
!
      class(memory_manager) :: mem
!
      real(dp), dimension(:,:,:,:,:), allocatable :: array
!
      integer, intent(in) :: M, N, O, P, Q ! First, second, third, fourth, fifth dimension of array 
!
      integer :: size_array ! Total size of array (M*N*O*P*Q)
      integer :: error = 0
!
      size_array = M*N*O*P*Q
!
!     Allocate array and check whether allocation was successful
!
      allocate(array(M,N,O,P,Q), stat = error)
!
      if (error .ne. 0) then
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
   end subroutine alloc_5_memory_manager
!
!
   subroutine dealloc_1_memory_manager(mem, array, M)
!!
!!    Dealloc (memory manager)
!!    Written by Rolf H. Myhre, January 2019
!!
!!    Deallocates a one dimensional double precision array and updates the available
!!    memory accordingly.
!!
      implicit none
!
      class(memory_manager) :: mem
!
      real(dp), dimension(:), allocatable :: array
!
      integer, intent(in) :: M ! Dimension of array 
!
      integer :: size_array ! Total size of array (M)
      integer :: error = 0
!
      size_array = M
!
!     Deallocate array and check whether deallocation was successful
!
      deallocate(array, stat = error)
!
      if (error .ne. 0) then
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
   end subroutine dealloc_1_memory_manager
!
!
   subroutine dealloc_2_memory_manager(mem, array, M, N)
!!
!!    Dealloc (memory manager)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Dec 2017
!!
!!    Deallocates a two dimensional double precision array and updates the available
!!    memory accordingly.
!!
      implicit none
!
      class(memory_manager) :: mem
!
      real(dp), dimension(:,:), allocatable :: array
!
      integer, intent(in) :: M, N ! First and second dimension of array 
!
      integer :: size_array ! Total size of array (M*N)
      integer :: error = 0
!
      size_array = M*N
!
!     Deallocate array and check whether deallocation was successful
!
      deallocate(array, stat = error)
!
      if (error .ne. 0) then
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
   end subroutine dealloc_2_memory_manager
!
!
   subroutine dealloc_3_memory_manager(mem, array, M, N, O)
!!
!!    Dealloc (memory manager)
!!    Written by Rolf H. Myhre, January 2019
!!
!!    Deallocates a three dimensional double precision array and updates the available
!!    memory accordingly.
!!
      implicit none
!
      class(memory_manager) :: mem
!
      real(dp), dimension(:,:,:), allocatable :: array
!
      integer, intent(in) :: M, N, O ! First, second and third dimension of array
!
      integer :: size_array ! Total size of array (M*N*O)
      integer :: error = 0
!
      size_array = M*N*O
!
!     Deallocate array and check whether deallocation was successful
!
      deallocate(array, stat = error)
!
      if (error .ne. 0) then
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
   end subroutine dealloc_3_memory_manager
!
!
   subroutine dealloc_4_memory_manager(mem, array, M, N, O, P)
!!
!!    Dealloc (memory manager)
!!    Written by Rolf H. Myhre, January 2019
!!
!!    Deallocates a four dimensional double precision array and updates the available
!!    memory accordingly.
!!
      implicit none
!
      class(memory_manager) :: mem
!
      real(dp), dimension(:,:,:,:), allocatable :: array
!
      integer, intent(in) :: M, N, O, P ! First, second, third and fourth dimension of array
!
      integer :: size_array ! Total size of array (M*N*O*P)
      integer :: error = 0
!
      size_array = M*N*O*P
!
!     Deallocate array and check whether deallocation was successful
!
      deallocate(array, stat = error)
!
      if (error .ne. 0) then
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
   end subroutine dealloc_4_memory_manager
!
!
   subroutine dealloc_5_memory_manager(mem, array, M, N, O, P, Q)
!!
!!    Dealloc (memory manager)
!!    Written by Rolf H. Myhre, January 2019
!!
!!    Deallocates a five dimensional double precision array and updates the available
!!    memory accordingly.
!!
      implicit none
!
      class(memory_manager) :: mem
!
      real(dp), dimension(:,:,:,:,:), allocatable :: array
!
      integer, intent(in) :: M, N, O, P, Q ! First, second, third, fourth, fifth dimension of array
!
      integer :: size_array ! Total size of array (M*N*O*P)
      integer :: error = 0
!
      size_array = M*N*O*P*Q
!
!     Deallocate array and check whether deallocation was successful
!
      deallocate(array, stat = error)
!
      if (error .ne. 0) then
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
   end subroutine dealloc_5_memory_manager
!
!
   subroutine alloc_int_1_memory_manager(mem, array, M)
!!
!!    Alloc int (memory manager)
!!    Written by Rolf H. Myhre, January 2019
!!
!!    Allocates a one dimensional integer array and updates the available
!!    memory accordingly.
!!
      implicit none
!
      class(memory_manager) :: mem
!
      integer, dimension(:), allocatable :: array
!
      integer, intent(in) :: M ! Dimension of array 
!
      integer :: size_array ! Total size of array (M)
      integer :: error = 0
      integer :: int_size
!
      size_array = M
!
!     Allocate array and check whether allocation was successful
!
      allocate(array(M), stat = error)
!
      if (error .ne. 0) then
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
      int_size = storage_size(array(1))/8
      mem%available = mem%available - int_size*size_array
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
   end subroutine alloc_int_1_memory_manager
!
!
   subroutine alloc_int_2_memory_manager(mem, array, M, N)
!!
!!    Alloc int (memory manager)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Dec 2017
!!
!!    Allocates a two dimensional integer array and updates the available
!!    memory accordingly.
!!
      implicit none
!
      class(memory_manager) :: mem
!
      integer, dimension(:,:), allocatable :: array
!
      integer, intent(in) :: M, N ! First and second dimension of array 
!
      integer :: size_array ! Total size of array (M*N)
      integer :: error = 0
      integer :: int_size
!
      size_array = M*N
!
!     Allocate array and check whether allocation was successful
!
      allocate(array(M,N), stat = error)
!
      if (error .ne. 0) then
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
      int_size = storage_size(array(1,1))/8
      mem%available = mem%available - int_size*size_array
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
   end subroutine alloc_int_2_memory_manager
!
!
   subroutine alloc_int_3_memory_manager(mem, array, M, N, O)
!!
!!    Alloc int (memory manager)
!!    Written by Rolf H. Myhre, January 2019
!!
!!    Allocates a three dimensional integer array and updates the available
!!    memory accordingly.
!!
      implicit none
!
      class(memory_manager) :: mem
!
      integer, dimension(:,:,:), allocatable :: array
!
      integer, intent(in) :: M, N, O ! First, second and third dimension of array 
!
      integer :: size_array ! Total size of array (M*N*O)
      integer :: error = 0
      integer :: int_size
!
      size_array = M*N*O
!
!     Allocate array and check whether allocation was successful
!
      allocate(array(M,N,O), stat = error)
!
      if (error .ne. 0) then
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
      int_size = storage_size(array(1,1,1))/8
      mem%available = mem%available - int_size*size_array
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
   end subroutine alloc_int_3_memory_manager
!
!
   subroutine alloc_int_4_memory_manager(mem, array, M, N, O, P)
!!
!!    Alloc int (memory manager)
!!    Written by Rolf H. Myhre, January 2019
!!
!!    Allocates a four dimensional integer array and updates the available
!!    memory accordingly.
!!
      implicit none
!
      class(memory_manager) :: mem
!
      integer, dimension(:,:,:,:), allocatable :: array
!
      integer, intent(in) :: M, N, O, P ! First, second, third and fourth dimension of array 
!
      integer :: size_array ! Total size of array (M*N*O*P)
      integer :: error = 0
      integer :: int_size
!
      size_array = M*N*O*P
!
!     Allocate array and check whether allocation was successful
!
      allocate(array(M,N,O,P), stat = error)
!
      if (error .ne. 0) then
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
      int_size = storage_size(array(1,1,1,1))/8
      mem%available = mem%available - int_size*size_array
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
   end subroutine alloc_int_4_memory_manager
!
!
   subroutine dealloc_int_1_memory_manager(mem, array, M)
!!
!!    Dealloc int (memory manager)
!!    Written by Rolf H. Myhre, January 2019
!!
!!    Deallocates a two dimensional integer array and updates the available
!!    memory accordingly.
!!
      implicit none
!
      class(memory_manager) :: mem
!
      integer, dimension(:), allocatable :: array
!
      integer, intent(in) :: M ! Dimension of array 
!
      integer :: size_array ! Total size of array (M*N)
      integer :: error = 0
      integer :: int_size
!
      size_array = M
!
!     Deallocate array and check whether deallocation was successful
!
      deallocate(array, stat = error)
!
      if (error .ne. 0) then
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
      int_size = storage_size(array(1))/8
      mem%available = mem%available + int_size*size_array
!
   end subroutine dealloc_int_1_memory_manager
!
!
   subroutine dealloc_int_2_memory_manager(mem, array, M, N)
!!
!!    Dealloc int (memory manager)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Dec 2017
!!
!!    Deallocates a two dimensional integer array and updates the available
!!    memory accordingly.
!!
      implicit none
!
      class(memory_manager) :: mem
!
      integer, dimension(:,:), allocatable :: array
!
      integer, intent(in) :: M, N ! First and second dimension of array 
!
      integer :: size_array ! Total size of array (M*N)
      integer :: error = 0
      integer :: int_size
!
      size_array = M*N
!
!     Deallocate array and check whether deallocation was successful
!
      deallocate(array, stat = error)
!
      if (error .ne. 0) then
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
      int_size = storage_size(array(1,1))/8
      mem%available = mem%available + int_size*size_array
!
   end subroutine dealloc_int_2_memory_manager
!
!
   subroutine dealloc_int_3_memory_manager(mem, array, M, N, O)
!!
!!    Dealloc int (memory manager)
!!    Written by Rolf H. Myhre, January 2019
!!
!!    Deallocates a three dimensional integer array and updates the available
!!    memory accordingly.
!!
      implicit none
!
      class(memory_manager) :: mem
!
      integer, dimension(:,:,:), allocatable :: array
!
      integer, intent(in) :: M, N, O ! First, second and third dimension of array
!
      integer :: size_array ! Total size of array (M*N*O)
      integer :: error = 0
      integer :: int_size
!
      size_array = M*N*O
!
!     Deallocate array and check whether deallocation was successful
!
      deallocate(array, stat = error)
!
      if (error .ne. 0) then
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
      int_size = storage_size(array(1,1,1))/8
      mem%available = mem%available + int_size*size_array
!
   end subroutine dealloc_int_3_memory_manager
!
!
   subroutine dealloc_int_4_memory_manager(mem, array, M, N, O, P)
!!
!!    Dealloc int (memory manager)
!!    Written by Rolf H. Myhre, January 2019
!!
!!    Deallocates a four dimensional integer array and updates the available
!!    memory accordingly.
!!
      implicit none
!
      class(memory_manager) :: mem
!
      integer, dimension(:,:,:,:), allocatable :: array
!
      integer, intent(in) :: M, N, O, P ! First, second, third and fourth dimension of array
!
      integer :: size_array ! Total size of array (M*N*O*P)
      integer :: error = 0
      integer :: int_size
!
      size_array = M*N*O*P
!
!     Deallocate array and check whether deallocation was successful
!
      deallocate(array, stat = error)
!
      if (error .ne. 0) then
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
      int_size = storage_size(array(1,1,1,1))/8
      mem%available = mem%available + int_size*size_array
!
   end subroutine dealloc_int_4_memory_manager
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
      call input%get_keyword_in_section('available', 'memory', mem%total)
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
      write(output%unit, '(t3, a38, i5, a)') 'Memory available for calculation:     ', &
                         mem%total/1000000000, ' GB'
!
   end subroutine print_settings_memory_manager
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
      integer, intent(in) :: req0
      integer, intent(in) :: req1
!
      integer, intent(in), optional :: element_size
!
      integer :: req0_tot
      integer :: req1_min
      integer :: req_min
!
      integer :: req_tot
!
      integer :: e_size
!
      if (.not. batch_p%initialized) then
!
         call output%error_msg('batch_setup_1 called on uninitialized batch')
!
      endif
!
      e_size = dp
      if(present(element_size)) then
         e_size = element_size
      endif
!
      req0_tot = req0*e_size
      req1_min = req1*e_size
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
         batch_p%max_length = int((mem%available - req0_tot)/(req1_min))
!
!        Number of full batches
!
         batch_p%num_batches = (batch_p%index_dimension-1)/(batch_p%max_length)+1
!
      endif
!
      if (force_batch) call batch_p%force_batch()
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
!!    req1_p: required memory that scales linearly with p batch size
!!    req1_q: required memory that scales linearly with q batch size
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
      integer, intent(in) :: req0
      integer, intent(in) :: req1_p
      integer, intent(in) :: req1_q
      integer, intent(in) :: req2
!
      integer, intent(in), optional :: element_size
!
      logical :: figured_out
!
      integer :: req0_tot
      integer :: req1_p_min
      integer :: req1_q_min 
      integer :: req2_min
      integer :: req_min
      integer :: req_tot 
!
      integer :: p_elements, q_elements
!
      integer :: e_size
!
      if ((.not. batch_p%initialized) .or. (.not. batch_q%initialized)) then
!
         call output%error_msg('batch_setup_2 called on uninitialized batch')
!
      endif
!
      e_size = dp
      if(present(element_size)) then
         e_size = element_size
      endif
!
      req0_tot   = req0*e_size
      req1_p_min = req1_p*e_size
      req1_q_min = req1_q*e_size
      req2_min = req2*e_size
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
         figured_out = .false.
         do while (.not. figured_out                              &
                     .and. p_elements .lt. batch_p%index_dimension &
                     .and. q_elements .lt. batch_q%index_dimension)
!
            if (((p_elements+1)*(q_elements+1)*req2_min &
                  + (p_elements+1)*req1_p_min          &
                  + (q_elements+1)*req1_q_min          &
                  + req0) .lt. mem%available) then
!
               p_elements = p_elements + 1 ! can hold +1 batch size
               q_elements = q_elements + 1
!
            else
!
               figured_out = .true.       ! cannot hold +1 batch size
!
            endif
!
         enddo
!
!
!        II. If simultaneous incrementation was not sufficient,
!             then try to increment the largest index further. This is
!             guaranteed to work, so let's just go ahead and increment
!             with no safeguards in place.
!
         if (.not. figured_out) then
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
            figured_out = .true.
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
!     Debug feature: enforced random batching 
!
      if (force_batch) then 
!
         if (batch_p%index_dimension == batch_q%index_dimension) then 
!
            call batch_p%force_batch()
!
            batch_q%max_length  = batch_p%max_length
            batch_q%num_batches = batch_p%num_batches
!
         else 
!
            call batch_p%force_batch()
            call batch_q%force_batch()
!
         endif 
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
      integer, intent(in) :: req0
! 
      integer, intent(in) :: req1_p 
      integer, intent(in) :: req1_q 
      integer, intent(in) :: req1_r 
!
      integer, intent(in) :: req2_pq
      integer, intent(in) :: req2_pr
      integer, intent(in) :: req2_qr 
!
      integer, intent(in) :: req3
!
      integer, intent(in), optional :: element_size
!
      integer :: req0_tot
!
      integer :: req1_p_min
      integer :: req1_q_min 
      integer :: req1_r_min
! 
      integer :: req2_pq_min
      integer :: req2_pr_min
      integer :: req2_qr_min
!
      integer :: req3_min
!
      integer :: req_min
      integer :: req_tot 
!
      integer :: p_elements, q_elements, r_elements
!
      logical :: found_batch_size, p_incremented, q_incremented, r_incremented
!
      integer :: e_size
!
      if ((.not. batch_p%initialized) .or. (.not. batch_q%initialized) .or. (.not. batch_r%initialized)) then
!
         call output%error_msg('batch_setup_3 called on uninitialized batch')
!
      endif
!
      e_size = dp
      if(present(element_size)) then
         e_size = element_size
      endif
!
      req0_tot   = req0*e_size
!
      req1_p_min = req1_p*e_size
      req1_q_min = req1_q*e_size
      req1_r_min = req1_r*e_size
!
      req2_pq_min = req2_pq*e_size
      req2_pr_min = req2_pr*e_size
      req2_qr_min = req2_qr*e_size
!
      req3_min = req3*e_size
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
                  + (p_elements)*(q_elements)*req2_pq_min &
                  + (p_elements)*(r_elements)*req2_pr_min &
                  + (q_elements)*(r_elements)*req2_qr_min &
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
      if (force_batch) then 
!
        call batch_p%force_batch()
        call batch_q%force_batch()
        call batch_r%force_batch()
!
      endif 
!
   end subroutine batch_setup_3_memory_manager
!
!
   subroutine batch_setup_3_ident_memory_manager(mem, batch_p, batch_q, batch_r, &
                                                 req0, req1, req2, req3, buffer_size, element_size)
!!
!!    Setup batching 
!!    This is setup for three batch indices
!!    with identical memory requirements
!!    Written by Rolf H. Myhre January 2019
!!
!!    Batching setup for three batching indices.
!!
!!    batch_p: Initialized batching object
!!    batch_q: Initialized batching object
!!    batch_r: Initialized batching object
!!
!!    req0: required memory that does not scale with batch size 
!!    req1: required memory for 1 of the 3 identical indices 
!!          scales linearly with batch size
!!          total req1 is 3*req1 for batches in p,q and r
!!    req2: required memory for 1 pair (e.g. pr) of the 3 identical indices
!!          scales quadratically with batch size
!!          total req2 is 6*req2 due to permutations of the indices (pq,qp,pr,rp,qr,rq)
!!    req3: required memory for 1 triple (e.g. rqp) of the 3 idenctical indices
!!          scales cubically with batch size
!!          total req3 is 6*req3 due to permutations of the indices (pqr,qpr,rqp,prq,qrp,rpq)
!!
!!    buffer_size: overwrite the default buffer size of 10%
!!
!!    element_size: memory per element, default is double precision
!!
!!    Be careful with symmetries and permutations!
!!
!!    example: CC3 omega/jacobian/jacobian_transpose
!!             use restricted loops over i≥j≥k to reduce the cost
!!             In case we use the integral g_ljci:
!!             1 array if we are not batching
!!             6 smaller arrays g_licj, g_lick, g_ljck, g_ljci, g_lkci, g_lkcj
!!             needed in case we are batching but:
!!             req2 = wf%n_o*wf%n_v 
!!             for this example because the factor of 6 is included in this routine
!!
      implicit none
!
      class(memory_manager) :: mem
!
      class(batching_index) :: batch_p ! An index being batched over
      class(batching_index) :: batch_q ! An index being batched over
      class(batching_index) :: batch_r ! An index being batched over
!
      integer, intent(in) :: req0
      integer, intent(in) :: req1 
      integer, intent(in) :: req2
      integer, intent(in) :: req3
!
      integer, intent(in), optional :: element_size
      real(dp), intent(in), optional :: buffer_size
!
      integer :: req0_tot
!
      integer :: req1_min
      integer :: req2_min
      integer :: req3_min
!
      integer :: req_min
      integer :: req_tot 
!
      integer :: elements
!
      logical :: found_batch_size, incremented
!
      integer :: e_size
      real(dp) :: buff
!
      if ((.not. batch_p%initialized) .or. (.not. batch_q%initialized) .or. (.not. batch_r%initialized)) then
!
         call output%error_msg('batch_setup_3 called on uninitialized batch')
!
      endif
!
      e_size = dp
      if(present(element_size)) then
         e_size = element_size
      endif
!
      buff = zero
      if(present(buffer_size)) then
         buff = buffer_size
      endif
!
      if (batch_p%index_dimension .ne. batch_q%index_dimension .or. &
          batch_p%index_dimension .ne. batch_r%index_dimension) then
!
         call output%error_msg('Index dimensions not identical in batch_setup_3_ident')
!
      endif
!
      req0_tot   = (req0 + int(req0*buff))*e_size
      req1_min   = (req1 + int(req1*buff))*e_size
      req2_min   = (req2 + int(req2*buff))*e_size
      req3_min   = (req3 + int(req3*buff))*e_size
!
!     Minimal required memory for batches of size 1
      req_min = req0_tot + 3*req1_min + 6*req2_min + 6*req3_min
!
!     Required memory to hold complete arrays in mem (no batching)
      req_tot = req0_tot + req1_min*(batch_p%index_dimension)     &
                         + req2_min*(batch_p%index_dimension)**2  &
                         + req3_min*(batch_p%index_dimension)**3
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
         elements = 1
         found_batch_size = .false.
         incremented = .true.
!
         do while (.not. found_batch_size .and. incremented)
!
            if ((elements) .lt. batch_p%index_dimension) then
               elements = elements + 1
               incremented = .true.
            else
               incremented = .false.
            endif
!
            if (  6*elements**3*req3_min  &
                + 6*elements**2*req2_min  &
                + 3*elements*req1_min     &
                + req0_tot .ge. mem%available) then 
!
               found_batch_size = .true.       ! cannot hold +1 batch size 
!
                  if (incremented) elements = elements - 1
!
            endif
!
         enddo
!
         batch_p%max_length = elements         
         batch_q%max_length = elements        
         batch_r%max_length = elements         
!
!        Figure out how many batches
!
         batch_p%num_batches = (batch_p%index_dimension-1)/(batch_p%max_length)+1
         batch_q%num_batches = (batch_q%index_dimension-1)/(batch_q%max_length)+1
         batch_r%num_batches = (batch_r%index_dimension-1)/(batch_r%max_length)+1
!
      endif
!
      if (force_batch) then 
!
         call batch_p%force_batch()
!
         batch_q%max_length  = batch_p%max_length
         batch_q%num_batches = batch_p%num_batches
!
         batch_r%max_length  = batch_p%max_length
         batch_r%num_batches = batch_p%num_batches
!
      endif 
!
   end subroutine batch_setup_3_ident_memory_manager
!
!
end module memory_manager_class
