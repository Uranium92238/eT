module memory_manager_class
!
!!
!!                   Memory manager class module                                 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Dec 2017         
!!
!!    The memory manager class handles the memory used by the model calculation,
!!    and there is an object called 'mem' in the wavefunction object of this class. 
!!
!!    Note: large arrays must always be allocated using the memory manager object. Small arrays,
!!          integers, strings, etc., which use a negligible amount of memory, are not considered 
!!          by the memory manager.
!! 
!
!  :::::::::::::::::::::::::::::::::::
!  -::- Modules used by the class -::-
!  :::::::::::::::::::::::::::::::::::
!
   use types
   use input_output
!
!  ::::::::::::::::::::::::::::::::::::::::::::::::
!  -::- Definition of the memory_manager class -::-
!  ::::::::::::::::::::::::::::::::::::::::::::::::
!
   type :: memory_manager 
!
!     The total amount of memory specified by user (standard: 30 GB)
!
      integer(i15) :: total = 30000000000 
!
!     The amount of memory currently available, based on the arrays currently allocated 
!     (memory used by objects and local variables are not included in this estimate)
!
      integer(i15) :: available
!
   contains 
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
   end type memory_manager                                                                            
!
!
contains
!
!
   module subroutine alloc_memory_manager(mem, array, M, N)
!!
!!    Alloc (Memory Manager)
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
      integer(i15) :: size ! Total size of array (M*N)
      integer(i15) :: stat = 0
      integer(i15) :: error = 0
!
      size = M*N 
!
!     Allocate array and check whether allocation was successful
!
      allocate(array(M,N), stat = error)
!
      if (stat .ne. 0) then 
!
         write(unit_output,'(t3,a,i15)') 'Error: could not allocate array with #elements =', size
         stop
!
      endif
!
!     Update the available memory 
!
!     The 'double precision' type (see types.F90) is typically 8 bytes,
!     though it might differ due to its definition in terms of precision.
!
      mem%available = mem%available - dp*size 
!
!     Check if there is no more memory (defined as being no more memory
!     left of what was specified by user as available)
!
      if (mem%available .lt. 0) then 
!
         write(unit_output,'(t3,a)') "Error: user-specified memory insufficient."
         stop
!
      endif         
!
   end subroutine alloc_memory_manager
!
!
   module subroutine dealloc_memory_manager(mem, array, M, N)
!!
!!    Dealloc (Memory Manager)
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
      integer(i15) :: size ! Total size of array (M*N)
      integer(i15) :: stat = 0
      integer(i15) :: error = 0
!
      size = M*N 
!
!     Deallocate array and check whether deallocation was successful
!
      deallocate(array, stat = error)
!
      if (stat .ne. 0) then 
!
         write(unit_output,'(t3,a,i15)') 'Error: could not deallocate array with #elements =', size
         stop
!
      endif
!
!     Update the available memory 
!
!     The 'double precision' type (see types.F90) is typically 8 bytes,
!     though it might differ due to its definition in terms of precision.
!
      mem%available = mem%available + dp*size       
!
   end subroutine dealloc_memory_manager
!
!
   module subroutine alloc_int_memory_manager(mem, array, M, N)
!!
!!    Alloc Int (Memory Manager)
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
      integer(i15) :: size ! Total size of array (M*N)
      integer(i15) :: stat = 0
      integer(i15) :: error = 0
!
      size = M*N 
!
!     Allocate array and check whether allocation was successful
!
      allocate(array(M,N), stat = error)
!
      if (stat .ne. 0) then 
!
         write(unit_output,'(t3,a,i15)') 'Error: could not allocate array with #elements =', size
         stop
!
      endif
!
!     Update the available memory 
!
!     The 'integer 15', or i15, type (see types.F90) is typically 4 bytes,
!     though it might differ due to its definition in terms of precision.
!
      mem%available = mem%available - i15*size 
!
!     Check if there is no more memory (defined as being no more memory
!     left of what was specified by user as available)
!
      if (mem%available .lt. 0) then 
!
         write(unit_output,'(t3,a)') "Error: user-specified memory insufficient."
         stop
!
      endif         
!
   end subroutine alloc_int_memory_manager
!
!
   module subroutine dealloc_int_memory_manager(mem, array, M, N)
!!
!!    Dealloc Int (Memory Manager)
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
      integer(i15) :: size ! Total size of array (M*N)
      integer(i15) :: stat = 0
      integer(i15) :: error = 0
!
      size = M*N 
!
!     Deallocate array and check whether deallocation was successful
!
      deallocate(array, stat = error)
!
      if (stat .ne. 0) then 
!
         write(unit_output,'(t3,a,i15)') 'Error: could not deallocate array with #elements =', size
         stop
!
      endif
!
!     Update the available memory 
!
!     The 'integer 15', or i15, type (see types.F90) is typically 4 bytes,
!     though it might differ due to its definition in terms of precision.
!
      mem%available = mem%available + i15*size       
!
   end subroutine dealloc_int_memory_manager
!
!
end module memory_manager_class