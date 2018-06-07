module workspace
!
!!
!!    Workspace module
!!    Written by Henrik Koch, Rolf H. Myhre, Eirik Kjønstad and Sarai Folkestad, Jan 2017
!!
!!    THIS IS BEING PHASED OUT & HAS BEEN REPLACED BY THE MEMORY MANAGER OBJECT.
!!    THE MANY CALLS TO "GET_AVAILABLE" MEANS THAT WE CANNOT REMOVE IT JUST YET.
!!
!!    Manages program memory usage and contains:
!!
!!    work_init:       Initializes the memory management variables.
!!
!!    allocator:       Allocation of double precission array of two dimmensions (M,N).
!!                     Updates memory management variables.
!!
!!    deallocator:     Deallocation of double precission array of two dimmensions (M,N).
!!                     Updates memory management variables.
!!
!!    allocator_int:   Allocation of integer array of two dimmensions (M,N).
!!                     Updates memory management variables.
!!
!!    deallocator_int: Deallocation of integer array of two dimmensions (M,N)
!!                     Updates memory management variables.
!!
!!    get_available:   Returns available memory.
!!
!
   use types
   use input_output
!
   implicit none
!
   integer, private :: work_length  = 0
   integer, private :: work_remains = 0
   integer, private :: work_used    = 0
!
   integer(i15)     :: mem = 4000000000 ! ca. 30 gb
!
!
contains
!
!
   subroutine work_init
!
!     Work Initilization
!
!     Initializes memory management variables
!
      implicit none
!
      work_length  = mem
      work_remains = mem
      work_used    = 0
!
   end subroutine work_init
!
!
   subroutine allocator(elm, M, N)
!
!     Allocator
!     Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2017
!
!     Allocates array and updates memory variables
!
      implicit none
!
      integer, intent(in)                    :: M, N
      real(dp), dimension(:,:), allocatable  :: elm
      integer                                :: size
      integer                                :: stat = 0, error = 0
!
      logical :: debug = .false.
!
      size = M*N
!
      allocate(elm(M,N),stat = error)
!
      if (stat .ne. 0) then
         write(unit_output,'(t3,a,i15)') 'Allocation error! Could not allocate array of size (M*N):', size
         stop
      endif
!
      if(debug) write(unit_output,*) work_remains, 4*size
!
      work_remains = work_remains - 4*size
      work_used    = work_used    + 4*size
!
      if (work_remains .lt. 0) then
         write(unit_output,'(t3,a)') "Error: user-specified memory insufficient."
         stop
      endif
!
   end subroutine allocator
!
!
   subroutine deallocator(elm, M, N)
!
!     Deallocator
!     Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2017
!
!     Deallocation and update of memory information
!
      implicit none
!
      real(dp), dimension(:,:), allocatable   :: elm
      integer                                 :: stat = 0, error = 0
      integer, intent(in)                     :: M, N
      integer                                 :: size
!
      size = M*N
!
      deallocate(elm,stat = error)
!
      if (stat .ne. 0) then
!
         write(unit_output,'(t3,a,i15)') 'Deallocation error! Could not deallocate array of size (M*N):', size
         stop
!
      endif
!
      work_remains = work_remains + 4*size
      work_used    = work_used    - 4*size
!
   end subroutine deallocator
!
!
   subroutine allocator_int(elm, M, N)
!
!     Allocator Integer Arrays
!     Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2017
!
!     Allocates array and updates memory information
!
      implicit none
!
      integer, intent(in)                  :: M, N
      integer, dimension(:,:), allocatable :: elm
      integer                              :: size
      integer                              :: stat = 0, error = 0
!
      size = M*N
!
      allocate(elm(M,N),stat = error)
!
      if (stat .ne. 0) then
!
         write(unit_output,'(t3,a,i15)') 'Allocation error! Could not allocate array of size (M*N):', size
         stop
!
      endif
!
      work_remains = work_remains - 2*size
      work_used    = work_used    + 2*size
!
      if (work_remains .lt. 0) then
!
         write(unit_output,'(t3,a)') "Error: user-specified memory insufficient."
         stop
!
      endif
!
   end subroutine allocator_int
!
!
   subroutine deallocator_int(elm, M, N)
!
!     Deallocator Integer Arrays
!     Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2017
!
!     Deallocates array and updates memory information
!
      implicit none
!
      integer, dimension(:,:), allocatable :: elm
      integer                              :: stat = 0, error = 0
      integer, intent(in)                  :: M, N
      integer                              :: size
!
      size = M*N
!
      deallocate(elm,stat = error)
!
      if (stat .ne. 0) then
!
         write(unit_output,'(t3,a,i15)') 'Deallocation error! Could not deallocate array of size (M*N):', size
         stop
!
      endif
!
      work_remains = work_remains + 2*size
      work_used    = work_used    - 2*size
!
   end subroutine deallocator_int
!
!
   integer function get_available()
!
!     Get Available
!     Written by Eirik F. Kjønstad and Sarai D. Folkestad, Mar 2017
!
!     Returns the available memory
!
      get_available = work_remains
!
   end function get_available
!
!
end module workspace
