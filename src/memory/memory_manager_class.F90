!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2022 the authors of eT
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
   use parameters
   use global_out, only : output
   use batching_index_class, only : batching_index
   use memory_tracker_class, only : memory_tracker
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
      integer(i64), private :: total
!
!     The amount of memory currently available, based on the arrays currently allocated
!     (memory used by objects and local variables are not included in this estimate)
!
      integer(i64), private :: available
!
!     Maximum amount of memory used at the same time
!
      integer(i64), private :: max_used
!
!     Unit for memory, default is GB
!
      character(len=2), private :: units
!
!     Batch memory tracking
!
      logical, private :: batching_on
      type(memory_tracker), allocatable :: batch_mem_tracker
!
   contains
!
!     Check if there are memory leaks - on exit of program
!
      procedure :: check_for_leak   => check_for_leak_memory_manager
!
!     Allocation and deallocation routines for arrays
!
      procedure :: alloc_r_1_memory_manager
      procedure :: alloc_r_2_memory_manager
      procedure :: alloc_r_3_memory_manager
      procedure :: alloc_r_4_memory_manager
      procedure :: alloc_r_5_memory_manager
      procedure :: alloc_c_1_memory_manager
      procedure :: alloc_c_2_memory_manager
      procedure :: alloc_c_3_memory_manager
      procedure :: alloc_c_4_memory_manager
      procedure :: alloc_i_1_memory_manager
      procedure :: alloc_i_2_memory_manager
      procedure :: alloc_l_1_memory_manager
      generic   :: alloc             => alloc_r_1_memory_manager,       &
                                        alloc_r_2_memory_manager,       &
                                        alloc_r_3_memory_manager,       &
                                        alloc_r_4_memory_manager,       &
                                        alloc_r_5_memory_manager,       &
                                        alloc_c_1_memory_manager,       &
                                        alloc_c_2_memory_manager,       &
                                        alloc_c_3_memory_manager,       &
                                        alloc_c_4_memory_manager,       &
                                        alloc_i_1_memory_manager,       &
                                        alloc_i_2_memory_manager,       &
                                        alloc_l_1_memory_manager
!
      procedure :: dealloc_r_1_memory_manager
      procedure :: dealloc_r_2_memory_manager
      procedure :: dealloc_r_3_memory_manager
      procedure :: dealloc_r_4_memory_manager
      procedure :: dealloc_r_5_memory_manager
      procedure :: dealloc_c_1_memory_manager
      procedure :: dealloc_c_2_memory_manager
      procedure :: dealloc_c_3_memory_manager
      procedure :: dealloc_c_4_memory_manager
      procedure :: dealloc_i_1_memory_manager
      procedure :: dealloc_i_2_memory_manager
      procedure :: dealloc_l_1_memory_manager
      generic   :: dealloc           => dealloc_r_1_memory_manager,     &
                                        dealloc_r_2_memory_manager,     &
                                        dealloc_r_3_memory_manager,     &
                                        dealloc_r_4_memory_manager,     &
                                        dealloc_r_5_memory_manager,     &
                                        dealloc_c_1_memory_manager,     &
                                        dealloc_c_2_memory_manager,     &
                                        dealloc_c_3_memory_manager,     &
                                        dealloc_c_4_memory_manager,     &
                                        dealloc_i_1_memory_manager,     &
                                        dealloc_i_2_memory_manager,     &
                                        dealloc_l_1_memory_manager
!
!     Routines for determining the number of batches
!
      procedure :: batch_setup_1_memory_manager
      procedure :: batch_setup_2_memory_manager
      procedure :: batch_setup_3_memory_manager
      generic   :: batch_setup       => batch_setup_1_memory_manager, &
                                        batch_setup_2_memory_manager, &
                                        batch_setup_3_memory_manager
!
      procedure :: batch_finalize &
                => batch_finalize_memory_manager
!
      procedure, private :: initialize_batching_tracker
!
      procedure :: update_memory_after_alloc        => update_memory_after_alloc_memory_manager
      procedure :: update_memory_after_dealloc      => update_memory_after_dealloc_memory_manager
!
      procedure, nopass :: print_allocation_error   => print_allocation_error_memory_manager
      procedure, nopass :: print_deallocation_error => print_deallocation_error_memory_manager
!
!     Print of settings
!
      procedure :: print_settings    => print_settings_memory_manager
!
!     Get and print of memory
!
      procedure :: print_max_used                  => print_max_used_memory_manager
      procedure :: print_available                 => print_available_memory_manager
      procedure :: get_available                   => get_available_memory_manager
      procedure, nopass :: get_memory_as_character => get_memory_as_character_memory_manager
!
   end type memory_manager
!
!
   interface memory_manager
!
      procedure :: new_memory_manager
!
   end interface memory_manager
!
!  Main memory object
!
   type(memory_manager) :: mem
!
contains
!
!
   function new_memory_manager(total, units) result(mem)
!!
!!    New memory manager
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Dec 2017
!!
!!    Creates the memory manager object and sets the
!!    total and initial available memory.
!!
      implicit none
!
      type(memory_manager) :: mem
!
      integer(i64), intent(in) :: total
!
      character(len=*), intent(in) :: units
!
!     Set standard and read settings
!
!     Default is 8 GB
!
      mem%total = total
      mem%units = trim(units)
!
!     Convert from current unit to B
!
      if (mem%units == 'gb') then
!
         mem%total =  mem%total*1000000000
!
      elseif (mem%units == 'mb') then
!
         mem%total =  mem%total*1000000
!
      elseif (mem%units == 'kb') then
!
         mem%total =  mem%total*1000
!
      elseif (trim(mem%units) == 'b') then
!
!       Do nothing
!
      else
!
         call output%error_msg('did not recognize the memory unit specified in input')
!
      endif
!
      mem%available = mem%total
      mem%max_used = mem%total - mem%available
!
      mem%batching_on = .false.
!
      call mem%print_settings()
!
   end function new_memory_manager
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
      character(len=200) :: difference_string
!
      if (mem%available .ne. mem%total) then
!
         call output%printf('m', 'Mismatch in memory according to eT and &
                            &specified on input:', fs='(/t3,a)')
!
         call output%printf('m', 'Memory available (eT):    (a0)', &
                            chars=[mem%get_memory_as_character(mem%available,.true.)], fs='(/t6,a)')
!
         call output%printf('m', 'Memory available (input): (a0)', &
                            chars=[mem%get_memory_as_character(mem%total,.true.)], fs='(t6,a)')
!
!
         difference_string = mem%get_memory_as_character(mem%total-mem%available,.true.)
         call output%printf('m', 'Difference:               (a0)', &
                            chars=[trim(difference_string)], fs='(t6,a)')
!
         call output%error_msg('Deallocations are missing or specified with &
                                 &incorrect dimensionalities.')
!
      endif
!
   end subroutine check_for_leak_memory_manager
!
!
   pure function get_available_memory_manager(mem) result(memory)
!!
!!    Get available
!!    Written by Eirik F. Kjønstad, Jan 2019
!!
      implicit none
!
      class(memory_manager), intent(in) :: mem
!
      integer(i64) :: memory
!
      memory = mem%available
!
   end function get_available_memory_manager
!
!
   pure function get_memory_as_character_memory_manager(input_mem, all_digits) result(memory)
!!
!!    Get available memory as character
!!    Written by Alexander C. Paul, Oct 2019
!!
!!    Receives an 8 byte integer containing the memory in byte.
!!    Returns character containing the same the number with a reasonable unit
!!
!!    If all_digits is .true. the full memory is returned in bytes
!!
      implicit none
!
      integer(i64), intent(in) :: input_mem
!
      logical, intent(in), optional :: all_digits
!
      character(len=17) :: memory
!
      logical :: all_digits_local
!
!     Print all digits? (i.e. give memory in B)
!
      all_digits_local = .false.
      if (present(all_digits)) all_digits_local = all_digits
!
      if (all_digits_local) then
!
            write(memory,'(i0, a)') input_mem, ' B'
            memory = trim(adjustl(memory))
!
      else if (abs(input_mem) .lt. 1d6) then
!
         write(memory,'(f10.3, a)') dble(input_mem)/1.0d3, ' KB'
         memory = trim(adjustl(memory))
!
      else if (abs(input_mem) .lt. 1d9) then
!
         write(memory,'(f10.6, a)') dble(input_mem)/1.0d6, ' MB'
         memory = trim(adjustl(memory))
!
      else if (abs(input_mem) .lt. 1d12) then
!
         write(memory,'(f10.6, a)') dble(input_mem)/1.0d9, ' GB'
         memory = trim(adjustl(memory))
!
      else if (abs(input_mem) .lt. 1d15) then
!
         write(memory,'(f13.6, a)') dble(input_mem)/1.0d12, ' TB'
         memory = trim(adjustl(memory))
!
      end if
!
   end function get_memory_as_character_memory_manager
!
!
   subroutine print_available_memory_manager(mem)
!!
!!    Print available
!!    Written by Eirik F. Kjønstad, Jan 2019
!!
      implicit none
!
      class(memory_manager), intent(in) :: mem
!

      call output%printf('m', 'Currently available memory: (a0)', &
                         chars=[mem%get_memory_as_character(mem%available, .true.)])
!
   end subroutine print_available_memory_manager
!
!
   subroutine print_max_used_memory_manager(mem)
!!
!!    Print maximum used memory
!!    Written by Alexander C. Paul, May 2020
!!
      implicit none
!
      class(memory_manager), intent(in) :: mem
!

      call output%printf('n', 'Peak memory usage during the execution of eT: (a0)', &
                         chars=[mem%get_memory_as_character(mem%max_used)], fs='(/t3,a)')
!
   end subroutine print_max_used_memory_manager
!
!
   subroutine alloc_r_1_memory_manager(mem, array, M)
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
      character(len=100) :: error_msg
!
      size_array = M
!
!     Allocate array and check whether allocation was successful
!
      allocate(array(M), stat = error, errmsg = error_msg)
!
      if (error .ne. 0) then
         call mem%print_allocation_error(size_array, error_msg)
      endif
!
!     Update the available memory
!
!     The 'double precision' type (see types.F90) is typically 8 bytes,
!     though it might differ due to its definition in terms of precision.
!
      call mem%update_memory_after_alloc(size_array, dp)
!
   end subroutine alloc_r_1_memory_manager
!
!
   subroutine alloc_r_2_memory_manager(mem, array, M, N)
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
      character(len=100) :: error_msg
!
      size_array = M*N
!
!     Allocate array and check whether allocation was successful
!
      allocate(array(M,N), stat = error, errmsg = error_msg)
!
      if (error .ne. 0) then
         call mem%print_allocation_error(size_array, error_msg)
      endif
!
!     Update the available memory
!
!     The 'double precision' type (see types.F90) is typically 8 bytes,
!     though it might differ due to its definition in terms of precision.
!
      call mem%update_memory_after_alloc(size_array, dp)
!
   end subroutine alloc_r_2_memory_manager
!
!
   subroutine alloc_r_3_memory_manager(mem, array, M, N, O)
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
      character(len=100) :: error_msg
!
      size_array = M*N*O
!
!     Allocate array and check whether allocation was successful
!
      allocate(array(M,N,O), stat = error, errmsg = error_msg)
!
      if (error .ne. 0) then
         call mem%print_allocation_error(size_array, error_msg)
      endif
!
!     Update the available memory
!
!     The 'double precision' type (see types.F90) is typically 8 bytes,
!     though it might differ due to its definition in terms of precision.
!
      call mem%update_memory_after_alloc(size_array, dp)
!
   end subroutine alloc_r_3_memory_manager
!
!
   subroutine alloc_r_4_memory_manager(mem, array, M, N, O, P)
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
      character(len=100) :: error_msg
!
      size_array = M*N*O*P
!
!     Allocate array and check whether allocation was successful
!
      allocate(array(M,N,O,P), stat = error, errmsg = error_msg)
!
      if (error .ne. 0) then
         call mem%print_allocation_error(size_array, error_msg)
      endif
!
!     Update the available memory
!
!     The 'double precision' type (see types.F90) is typically 8 bytes,
!     though it might differ due to its definition in terms of precision.
!
      call mem%update_memory_after_alloc(size_array, dp)
!
   end subroutine alloc_r_4_memory_manager
!
!
   subroutine alloc_r_5_memory_manager(mem, array, M, N, O, P, Q)
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
      character(len=100) :: error_msg
!
      size_array = M*N*O*P*Q
!
!     Allocate array and check whether allocation was successful
!
      allocate(array(M,N,O,P,Q), stat = error, errmsg = error_msg)
!
      if (error .ne. 0) then
         call mem%print_allocation_error(size_array, error_msg)
      endif
!
!     Update the available memory
!
!     The 'double precision' type (see types.F90) is typically 8 bytes,
!     though it might differ due to its definition in terms of precision.
!
      call mem%update_memory_after_alloc(size_array, dp)
!
   end subroutine alloc_r_5_memory_manager
!
!
   subroutine alloc_c_1_memory_manager(mem, array, M)
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
      complex(dp), dimension(:), allocatable :: array
!
      integer, intent(in) :: M ! Dimension of array that is being allocated
!
      integer :: size_array ! Total size of array (M)
      integer :: error = 0
!
      character(len=100) :: error_msg
!
      size_array = M
!
!     Allocate array and check whether allocation was successful
!
      allocate(array(M), stat = error, errmsg = error_msg)
!
      if (error .ne. 0) then
         call mem%print_allocation_error(size_array, error_msg)
      endif
!
!     Update the available memory
!
!     The 'double precision' type (see types.F90) is typically 8 bytes,
!     though it might differ due to its definition in terms of precision.
!
      call mem%update_memory_after_alloc(size_array, 2*dp)
!
   end subroutine alloc_c_1_memory_manager
!
!
   subroutine alloc_c_2_memory_manager(mem, array, M, N)
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
      complex(dp), dimension(:,:), allocatable :: array
!
      integer, intent(in) :: M, N ! First and second dimension of array that is being allocated
!
      integer :: size_array ! Total size of array (M*N)
      integer :: error = 0
!
      character(len=100) :: error_msg
!
      size_array = M*N
!
!     Allocate array and check whether allocation was successful
!
      allocate(array(M,N), stat = error, errmsg = error_msg)
!
      if (error .ne. 0) then
         call mem%print_allocation_error(size_array, error_msg)
      endif
!
!     Update the available memory
!
!     The 'double precision' type (see types.F90) is typically 8 bytes,
!     though it might differ due to its definition in terms of precision.
!
      call mem%update_memory_after_alloc(size_array, 2*dp)
!
   end subroutine alloc_c_2_memory_manager
!
!
   subroutine alloc_c_3_memory_manager(mem, array, M, N, O)
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
      complex(dp), dimension(:,:,:), allocatable :: array
!
      integer, intent(in) :: M, N, O ! First, second and third dimension of array
!
      integer :: size_array ! Total size of array (M*N*O)
      integer :: error = 0
!
      character(len=100) :: error_msg
!
      size_array = M*N*O
!
!     Allocate array and check whether allocation was successful
!
      allocate(array(M,N,O), stat = error, errmsg = error_msg)
!
      if (error .ne. 0) then
         call mem%print_allocation_error(size_array, error_msg)
      endif
!
!     Update the available memory
!
!     The 'double precision' type (see types.F90) is typically 8 bytes,
!     though it might differ due to its definition in terms of precision.
!
      call mem%update_memory_after_alloc(size_array, 2*dp)
!
   end subroutine alloc_c_3_memory_manager
!
!
   subroutine alloc_c_4_memory_manager(mem, array, M, N, O, P)
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
      complex(dp), dimension(:,:,:,:), allocatable :: array
!
      integer, intent(in) :: M, N, O, P ! First, second, third and fourth dimension of array
!
      integer :: size_array ! Total size of array (M*N*O*P)
      integer :: error = 0
!
      character(len=100) :: error_msg
!
      size_array = M*N*O*P
!
!     Allocate array and check whether allocation was successful
!
      allocate(array(M,N,O,P), stat = error, errmsg = error_msg)
!
      if (error .ne. 0) then
         call mem%print_allocation_error(size_array, error_msg)
      endif
!
!     Update the available memory
!
!     The 'double precision' type (see types.F90) is typically 8 bytes,
!     though it might differ due to its definition in terms of precision.
!
      call mem%update_memory_after_alloc(size_array, 2*dp)
!
   end subroutine alloc_c_4_memory_manager
!
!
   subroutine dealloc_r_1_memory_manager(mem, array, M)
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
      character(len=100) :: error_msg
!
      size_array = M
!
!     Deallocate array and check whether deallocation was successful
!
      deallocate(array, stat = error, errmsg = error_msg)
!
      if (error .ne. 0) then
         call mem%print_deallocation_error(size_array, error_msg)
      endif
!
!     Update the available memory
!
!     The 'double precision' type (see types.F90) is typically 8 bytes,
!     though it might differ due to its definition in terms of precision.
!
      call mem%update_memory_after_dealloc(size_array, dp)
!
   end subroutine dealloc_r_1_memory_manager
!
!
   subroutine dealloc_r_2_memory_manager(mem, array, M, N)
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
      character(len=100) :: error_msg
!
      size_array = M*N
!
!     Deallocate array and check whether deallocation was successful
!
      deallocate(array, stat = error, errmsg = error_msg)
!
      if (error .ne. 0) then
         call mem%print_deallocation_error(size_array, error_msg)
      endif
!
!     Update the available memory
!
!     The 'double precision' type (see types.F90) is typically 8 bytes,
!     though it might differ due to its definition in terms of precision.
!
      call mem%update_memory_after_dealloc(size_array, dp)
!
   end subroutine dealloc_r_2_memory_manager
!
!
   subroutine dealloc_r_3_memory_manager(mem, array, M, N, O)
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
      character(len=100) :: error_msg
!
      size_array = M*N*O
!
!     Deallocate array and check whether deallocation was successful
!
      deallocate(array, stat = error, errmsg = error_msg)
!
      if (error .ne. 0) then
         call mem%print_deallocation_error(size_array, error_msg)
      endif
!
!     Update the available memory
!
!     The 'double precision' type (see types.F90) is typically 8 bytes,
!     though it might differ due to its definition in terms of precision.
!
      call mem%update_memory_after_dealloc(size_array, dp)
!
   end subroutine dealloc_r_3_memory_manager
!
!
   subroutine dealloc_r_4_memory_manager(mem, array, M, N, O, P)
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
      character(len=100) :: error_msg
!
      size_array = M*N*O*P
!
!     Deallocate array and check whether deallocation was successful
!
      deallocate(array, stat = error, errmsg = error_msg)
!
      if (error .ne. 0) then
         call mem%print_deallocation_error(size_array, error_msg)
      endif
!
!     Update the available memory
!
!     The 'double precision' type (see types.F90) is typically 8 bytes,
!     though it might differ due to its definition in terms of precision.
!
      call mem%update_memory_after_dealloc(size_array, dp)
!
   end subroutine dealloc_r_4_memory_manager
!
!
   subroutine dealloc_r_5_memory_manager(mem, array, M, N, O, P, Q)
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
      integer :: size_array ! Total size of array (M*N*O*P*Q)
      integer :: error = 0
!
      character(len=100) :: error_msg
!
      size_array = M*N*O*P*Q
!
!     Deallocate array and check whether deallocation was successful
!
      deallocate(array, stat = error, errmsg = error_msg)
!
      if (error .ne. 0) then
         call mem%print_deallocation_error(size_array, error_msg)
      endif
!
!     Update the available memory
!
!     The 'double precision' type (see types.F90) is typically 8 bytes,
!     though it might differ due to its definition in terms of precision.
!
      call mem%update_memory_after_dealloc(size_array, dp)
!
   end subroutine dealloc_r_5_memory_manager
!
!
   subroutine dealloc_c_1_memory_manager(mem, array, M)
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
      complex(dp), dimension(:), allocatable :: array
!
      integer, intent(in) :: M ! Dimension of array
!
      integer :: size_array ! Total size of array (M)
      integer :: error = 0
!
      character(len=100) :: error_msg
!
      size_array = M
!
!     Deallocate array and check whether deallocation was successful
!
      deallocate(array, stat = error, errmsg = error_msg)
!
      if (error .ne. 0) then
         call mem%print_deallocation_error(size_array, error_msg)
      endif
!
!     Update the available memory
!
!     The 'double precision' type (see types.F90) is typically 8 bytes,
!     though it might differ due to its definition in terms of precision.
!
      call mem%update_memory_after_dealloc(size_array, 2*dp)
!
   end subroutine dealloc_c_1_memory_manager
!
!
   subroutine dealloc_c_2_memory_manager(mem, array, M, N)
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
      complex(dp), dimension(:,:), allocatable :: array
!
      integer, intent(in) :: M, N ! First and second dimension of array
!
      integer :: size_array ! Total size of array (M*N)
      integer :: error = 0
!
      character(len=100) :: error_msg
!
      size_array = M*N
!
!     Deallocate array and check whether deallocation was successful
!
      deallocate(array, stat = error, errmsg = error_msg)
!
      if (error .ne. 0) then
         call mem%print_deallocation_error(size_array, error_msg)
      endif
!
!     Update the available memory
!
!     The 'double precision' type (see types.F90) is typically 8 bytes,
!     though it might differ due to its definition in terms of precision.
!
      call mem%update_memory_after_dealloc(size_array, 2*dp)
!
   end subroutine dealloc_c_2_memory_manager
!
!
   subroutine dealloc_c_3_memory_manager(mem, array, M, N, O)
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
      complex(dp), dimension(:,:,:), allocatable :: array
!
      integer, intent(in) :: M, N, O ! First, second and third dimension of array
!
      integer :: size_array ! Total size of array (M*N*O)
      integer :: error = 0
!
      character(len=100) :: error_msg
!
      size_array = M*N*O
!
!     Deallocate array and check whether deallocation was successful
!
      deallocate(array, stat = error, errmsg = error_msg)
!
      if (error .ne. 0) then
         call mem%print_deallocation_error(size_array, error_msg)
      endif
!
!     Update the available memory
!
!     The 'double precision' type (see types.F90) is typically 8 bytes,
!     though it might differ due to its definition in terms of precision.
!
      call mem%update_memory_after_dealloc(size_array, 2*dp)
!
   end subroutine dealloc_c_3_memory_manager
!
!
   subroutine dealloc_c_4_memory_manager(mem, array, M, N, O, P)
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
      complex(dp), dimension(:,:,:,:), allocatable :: array
!
      integer, intent(in) :: M, N, O, P ! First, second, third and fourth dimension of array
!
      integer :: size_array ! Total size of array (M*N*O*P)
      integer :: error = 0
!
      character(len=100) :: error_msg
!
      size_array = M*N*O*P
!
!     Deallocate array and check whether deallocation was successful
!
      deallocate(array, stat = error, errmsg = error_msg)
!
      if (error .ne. 0) then
         call mem%print_deallocation_error(size_array, error_msg)
      endif
!
!     Update the available memory
!
!     The 'double precision' type (see types.F90) is typically 8 bytes,
!     though it might differ due to its definition in terms of precision.
!
      call mem%update_memory_after_dealloc(size_array, 2*dp)
!
   end subroutine dealloc_c_4_memory_manager
!
!
   subroutine alloc_i_1_memory_manager(mem, array, M)
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
!
      character(len=100) :: error_msg
!
      size_array = M
!
!     Allocate array and check whether allocation was successful
!
      allocate(array(M), stat = error, errmsg = error_msg)
!
      if (error .ne. 0) then
         call mem%print_allocation_error(size_array, error_msg)
      endif
!
!     Update the available memory
!     Check integer size
!
      call mem%update_memory_after_alloc(size_array, int_size)
!
   end subroutine alloc_i_1_memory_manager
!
!
   subroutine alloc_i_2_memory_manager(mem, array, M, N)
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
!
      character(len=100) :: error_msg
!
      size_array = M*N
!
!     Allocate array and check whether allocation was successful
!
      allocate(array(M,N), stat = error, errmsg = error_msg)
!
      if (error .ne. 0) then
         call mem%print_allocation_error(size_array, error_msg)
      endif
!
!     Update the available memory
!     Check integer size
!
      call mem%update_memory_after_alloc(size_array, int_size)
!
   end subroutine alloc_i_2_memory_manager
!
!
   subroutine dealloc_i_1_memory_manager(mem, array, M)
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
!
      character(len=100) :: error_msg
!
      size_array = M
!
!     Deallocate array and check whether deallocation was successful
!
      deallocate(array, stat = error, errmsg = error_msg)
!
      if (error .ne. 0) then
         call mem%print_deallocation_error(size_array, error_msg)
      endif
!
!     Update the available memory
!     Check integer size
!
      mem%available = mem%available + int_size*size_array
!
   end subroutine dealloc_i_1_memory_manager
!
!
   subroutine dealloc_i_2_memory_manager(mem, array, M, N)
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
!
      character(len=100) :: error_msg
!
      size_array = M*N
!
!     Deallocate array and check whether deallocation was successful
!
      deallocate(array, stat = error, errmsg = error_msg)
!
      if (error .ne. 0) then
         call mem%print_deallocation_error(size_array, error_msg)
      endif
!
!     Update the available memory
!     Check integer size
!
      mem%available = mem%available + int_size*size_array
!
   end subroutine dealloc_i_2_memory_manager
!
!
   subroutine alloc_l_1_memory_manager(mem, array, M)
!!
!!    Alloc log (memory manager)
!!    Written by Rolf H. Myhre, September 2019
!!
!!    Allocates a one dimensional logical array and updates the available
!!    memory accordingly.
!!
      implicit none
!
      class(memory_manager) :: mem
!
      logical, dimension(:), allocatable :: array
!
      integer, intent(in) :: M ! Dimension of array
!
      integer :: size_array ! Total size of array (M)
      integer :: error = 0
      integer :: log_size
!
      character(len=100) :: error_msg
!
      size_array = M
!
!     Allocate array and check whether allocation was successful
!
      allocate(array(M), stat = error, errmsg = error_msg)
!
      if (error .ne. 0) then
         call mem%print_allocation_error(size_array, error_msg)
      endif
!
!     Update the available memory
!     Figure out how big a logical is.
!
      log_size = storage_size(array(1))/8
      call mem%update_memory_after_alloc(size_array, log_size)
!
   end subroutine alloc_l_1_memory_manager
!
!
   subroutine dealloc_l_1_memory_manager(mem, array, M)
!!
!!    Dealloc log (memory manager)
!!    Written by Rolf H. Myhre, September 2019
!!
!!    Deallocates a one dimensional logical array and updates the available
!!    memory accordingly.
!!
      implicit none
!
      class(memory_manager) :: mem
!
      logical, dimension(:), allocatable :: array
!
      integer, intent(in) :: M ! Dimension of array
!
      integer :: size_array ! Total size of array (M*N)
      integer :: error = 0
      integer :: log_size
!
      character(len=100) :: error_msg
!
      size_array = M
!
!     Deallocate array and check whether deallocation was successful
!
      deallocate(array, stat = error, errmsg = error_msg)
!
      if (error .ne. 0) then
         call mem%print_deallocation_error(size_array, error_msg)
      endif
!
!     Update the available memory
!
      log_size = storage_size(array(1))/8
      mem%available = mem%available + log_size*size_array
!
   end subroutine dealloc_l_1_memory_manager
!
!
   subroutine print_allocation_error_memory_manager(size_array, error_msg)
!!
!!    Check allocation error
!!    Written by Alexander C. Paul, March 2020
!!
      implicit none
!
      integer, intent(in) :: size_array
      character (len=*), intent(in) :: error_msg
!
      call output%printf('m', error_msg, fs='(/t3,a)')
      call output%printf('m', 'Note: Error message from gfortran might not be accurate.', &
                          fs='(t3,a)')
      call output%error_msg('Could not allocate array with #elements = (i0).', &
                             ints=[size_array], ffs='(t3,a)')
!
   end subroutine print_allocation_error_memory_manager
!
!
   subroutine print_deallocation_error_memory_manager(size_array, error_msg)
!!
!!    Check deallocation error
!!    Written by Alexander C. Paul, March 2020
!!
      implicit none
!
      integer, intent(in) :: size_array
      character (len=*), intent(in) :: error_msg
!
      call output%printf('m', error_msg)
      call output%printf('m', 'Note: Error message from gfortran might not be accurate.')
      call output%error_msg('Could not deallocate array with #elements = (i0).', &
                             ints=[size_array])
!
   end subroutine print_deallocation_error_memory_manager
!
!
   subroutine update_memory_after_alloc_memory_manager(mem, size_array, size_type)
!!
!!    Update memory after allocation
!!    Written by Alexander C. Paul, May 2020
!!
!!    size_array: total size of the array allocated
!!    size_type : storage size of one element of the array in Byte
!!
      implicit none
!
      class(memory_manager) :: mem
!
      integer, intent(in) :: size_array, size_type
!
      integer(i64) :: bytes
!
      bytes = int(size_array*size_type, kind=i64)
!
      mem%available = mem%available - bytes
!
      if (mem%available .lt. 0) then
!
         call output%error_msg('User-specified memory insufficient in mem%alloc. &
                               &Tried to allocate array with (i0) elements.', &
                                ints=[size_array], ll=50)
!
      endif
!
!     Update max used memory if needed
      if (mem%max_used < (mem%total - mem%available)) &
          mem%max_used =  mem%total - mem%available
!
      if (mem%batching_on) call mem%batch_mem_tracker%update(bytes)
!
   end subroutine update_memory_after_alloc_memory_manager
!
!
   subroutine update_memory_after_dealloc_memory_manager(mem, size_array, size_type)
!!
!!    Update memory after deallocation
!!    Written by Alexander C. Paul, May 2020
!!
!!    size_array: total size of the array allocated
!!    size_type : storage size of one element of the array in Byte
!!
      implicit none
!
      class(memory_manager) :: mem
!
      integer, intent(in) :: size_array, size_type
!
      integer(i64) :: bytes
!
      bytes = int(size_array*size_type, kind=i64)
!
      mem%available = mem%available + bytes
!
      if (mem%batching_on) call mem%batch_mem_tracker%update(-bytes)
!
   end subroutine update_memory_after_dealloc_memory_manager
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
      call output%printf('m', 'Memory available for calculation: ' //  &
                         mem%get_memory_as_character(mem%total))
!
   end subroutine print_settings_memory_manager
!
!
   subroutine batch_finalize_memory_manager(mem)
!!
!!    Batch finalize
!!    Written by Eirik F. Kjønstad, June 2021
!!
!!    Must be called after a batching loop is finished.
!!
!!    The routine turns of batching mode and deallocates the
!!    memory tracker for the batching procedure.
!!
      implicit none
!
      class(memory_manager), intent(inout) :: mem
!
      mem%batching_on = .false.
!
      if (allocated(mem%batch_mem_tracker)) then
!
         deallocate(mem%batch_mem_tracker)
!
      else
!
         call output%error_msg('Asked to finalize batch, but batching tracker &
                               &not allocated! Was batch_finalize already called &
                               &for the current batching loop?')
!
      endif
!
   end subroutine batch_finalize_memory_manager
!
!
   subroutine initialize_batching_tracker(mem, max_memory_usage, tag)
!!
!!    Initialize batching tracker
!!    Written by Eirik F. Kjønstad, June 2021
!!
!!    To be called when batching has been determined.
!!
!!    Makes sure memory usage is tracked during the batching loops.
!!
      implicit none
!
      class(memory_manager), intent(inout) :: mem
!
      integer(i64), intent(in) :: max_memory_usage
      character(len=*), intent(in) :: tag
!
      if (mem%batching_on) then
!
         call output%error_msg('Tried to initialize memory tracker for batching loop, &
                               &but the memory manager is already in batching mode! &
                               &Have you forgotten to finalize the previous batching loop?')
!
      endif
!
      mem%batching_on = .true.
      mem%batch_mem_tracker = memory_tracker(max_memory_usage, tag)
!
   end subroutine initialize_batching_tracker
!
!
   subroutine batch_setup_1_memory_manager(mem, batch_p, req0, req1, tag, element_size)
!!
!!    Batch setup 1
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
      character(len=*), intent(in) :: tag
!
      integer, intent(in), optional :: element_size
!
!
      integer(i64):: req0_tot
      integer(i64):: req1_min
      integer(i64):: req_min
!
      integer(i64):: req_tot
!
      integer :: e_size
      character(len=17), allocatable :: reqChar
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
      req0_tot = int(req0*e_size, kind=i64)
      req1_min = int(req1*e_size, kind=i64)
!
      req_min = req0_tot + req1_min
      req_tot = req0_tot + req1_min*int(batch_p%index_dimension, kind=i64)
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
!        Hack because intel flips out if we put two functions in chars=[]
!
         reqChar = mem%get_memory_as_character(req_min, .true.)
         call output%printf('m', 'Need at least (a0) but only have (a0)', &
                            chars=[reqChar, mem%get_memory_as_character(mem%available, .true.)])
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
      if (batch_p%num_batches > 1) call output%printf('v', 'Batching in (a0)', chars=[tag])
!
      call mem%initialize_batching_tracker(req0_tot + req1_min*int(batch_p%max_length, kind=i64), &
                                           tag)
!
   end subroutine batch_setup_1_memory_manager
!
!
   subroutine batch_setup_2_memory_manager(mem, batch_p, batch_q, req0, req1_p, req1_q, &
                                           req2, tag, element_size, req_single_batch)
!!
!!    Batch setup 2
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
!!    req_single_batch: optional specifying the minimal memory needed to not batch
!!
!!    If you are batching over i and j and need to keep g_abij, g_abci and g_abcj in memory,
!!    req1_i = n_v**3, req1_j = n_v**3 and req2 = n_v**2.
!!    Memory per batch is then batch_size*(req1_i + req1_j) + batch_size**2*req2
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
      character(len=*), intent(in) :: tag
!
      integer, intent(in), optional :: element_size
      integer, intent(in), optional :: req_single_batch
!
!
      logical :: figured_out
!
      integer(i64):: req0_tot
      integer(i64):: req1_p_min
      integer(i64):: req1_q_min
      integer(i64):: req2_min
      integer(i64):: req_min
      integer(i64):: req_tot
!
      integer(i64):: p_elements, q_elements
!
      integer :: e_size
      character(len=17), allocatable :: reqChar
!
      integer(i64) :: max_memory_usage
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
      req0_tot   = int(req0*e_size, kind=i64)
      req1_p_min = int(req1_p*e_size, kind=i64)
      req1_q_min = int(req1_q*e_size, kind=i64)
      req2_min   = int(req2*e_size, kind=i64)
!
      req_min = req0_tot + req1_p_min + req1_q_min + req2_min
!
!     Determine or copy the memory needed to not batch
!
      if (present(req_single_batch)) then
         req_tot = int(req_single_batch*e_size, kind=i64)
      else
!
         req_tot = req0_tot + req1_p_min*int(batch_p%index_dimension, kind=i64)  &
                  + req1_q_min*int(batch_q%index_dimension, kind=i64)  &
                  + req2_min*int(batch_p%index_dimension, kind=i64)    &
                     *int(batch_q%index_dimension, kind=i64)
!
      end if
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
!        Hack because intel flips out if we put two functions in chars=[]
         reqChar = mem%get_memory_as_character(req_min, .true.)
         call output%printf('m', 'Need at least (a0) but only have (a0)', &
                            chars=[reqChar, mem%get_memory_as_character(mem%available, .true.)])
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
         do while (.not. figured_out                                    &
                     .and. int(p_elements) .lt. batch_p%index_dimension &
                     .and. int(q_elements) .lt. batch_q%index_dimension)
!
            if (((p_elements+1)*(q_elements+1)*req2_min &
                  + (p_elements+1)*req1_p_min          &
                  + (q_elements+1)*req1_q_min          &
                  + req0_tot) .lt. mem%available) then
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
                           + req0_tot) .lt. mem%available)
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
                           + req0_tot) .lt. mem%available)
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
         batch_p%max_length = int(p_elements)
         batch_q%max_length = int(q_elements)
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
      max_memory_usage = req0_tot + &
                         req1_p_min*int(batch_p%max_length, kind=i64) + &
                         req1_q_min*int(batch_q%max_length, kind=i64) + &
                         req2_min*int(batch_q%max_length, kind=i64)     &
                                 *int(batch_q%max_length, kind=i64)
!
      if (batch_p%num_batches .eq. 1 .and. &
          batch_q%num_batches .eq. 1) then
!
         if (present(req_single_batch)) then
!
            max_memory_usage = int(req_single_batch*e_size, kind=i64)
!
         endif
!
      endif
!
      if (any([batch_p%num_batches, batch_q%num_batches] > 1)) then
         call output%printf('v', 'Batching in (a0)', chars=[tag])
      end if
!
      call mem%initialize_batching_tracker(max_memory_usage, tag)
!
   end subroutine batch_setup_2_memory_manager
!
!
   subroutine batch_setup_3_memory_manager(mem, batch_p, batch_q, batch_r, req0, &
                                           req1_p, req1_q, req1_r, req2_pq,      &
                                           req2_pr, req2_qr, req3, tag,          &
                                           element_size, req_single_batch)
!!
!!    Batch setup 3
!!    This is setup for three batch indices
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
!!    req_single_batch: optional specifying the minimal memory needed to not batch
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
      character(len=*), intent(in) :: tag
!
      integer, intent(in), optional :: element_size
      integer, intent(in), optional :: req_single_batch
!
!
      integer(i64):: req0_tot
!
      integer(i64):: req1_p_min
      integer(i64):: req1_q_min
      integer(i64):: req1_r_min
!
      integer(i64):: req2_pq_min
      integer(i64):: req2_pr_min
      integer(i64):: req2_qr_min
!
      integer(i64):: req3_min
!
      integer(i64):: req_min
      integer(i64):: req_tot
!
      integer(i64):: p_elements, q_elements, r_elements
!
      logical :: found_batch_size, p_incremented, q_incremented, r_incremented
!
      integer :: e_size
      character(len=17), allocatable :: reqChar
!
      integer(i64) :: max_memory_usage
!
      if ((.not. batch_p%initialized)        &
            .or. (.not. batch_q%initialized) &
            .or. (.not. batch_r%initialized)) then
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
      req0_tot   = int(req0*e_size, kind=i64)
!
      req1_p_min = int(req1_p*e_size, kind=i64)
      req1_q_min = int(req1_q*e_size, kind=i64)
      req1_r_min = int(req1_r*e_size, kind=i64)
!
      req2_pq_min = int(req2_pq*e_size, kind=i64)
      req2_pr_min = int(req2_pr*e_size, kind=i64)
      req2_qr_min = int(req2_qr*e_size, kind=i64)
!
      req3_min = int(req3*e_size, kind=i64)
!
      req_min = req0_tot + req1_p_min + req1_q_min + req1_r_min &
              + req2_pq_min + req2_pr_min + req2_qr_min + req3_min
!
!     Determine or copy the memory needed to not batch
!
      if (present(req_single_batch)) then
         req_tot = int(req_single_batch*e_size, kind=i64)
      else
!
         req_tot = req0_tot + req1_p_min                              &
                              *int(batch_p%index_dimension, kind=i64) &
                            + req1_q_min                              &
                              *int(batch_q%index_dimension, kind=i64) &
                            + req1_r_min                              &
                              *int(batch_r%index_dimension, kind=i64) &
                            + req2_pq_min                             &
                              *int(batch_p%index_dimension, kind=i64) &
                              *int(batch_q%index_dimension, kind=i64) &
                            + req2_pr_min                             &
                              *int(batch_p%index_dimension, kind=i64) &
                              *int(batch_r%index_dimension, kind=i64) &
                            + req2_qr_min                             &
                              *int(batch_q%index_dimension, kind=i64) &
                              *int(batch_r%index_dimension, kind=i64) &
                            + req3_min                                &
                              *int(batch_p%index_dimension, kind=i64) &
                              *int(batch_q%index_dimension, kind=i64) &
                              *int(batch_r%index_dimension, kind=i64)
!
      end if
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
!        Hack because intel flips out if we put two functions in chars=[]
         reqChar = mem%get_memory_as_character(req_min, .true.)
         call output%printf('m', 'Need at least (a0) but only have (a0)', &
                            chars=[reqChar, mem%get_memory_as_character(mem%available, .true.)])
         call output%error_msg('Not enough memory for a batch.')
!
      else
!
!        First, try to increment all indices simultaneously
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
         do while (.not. found_batch_size &
                   .and. (p_incremented .or. q_incremented .or. r_incremented))
!
            if (int(p_elements) .lt. batch_p%index_dimension) then
               p_elements = p_elements + 1
               p_incremented = .true.
            else
               p_incremented = .false.
            endif
!
            if (int(q_elements) .lt. batch_q%index_dimension) then
               q_elements = q_elements + 1
               q_incremented = .true.
            else
               q_incremented = .false.
            endif
!
            if (int(r_elements) .lt. batch_r%index_dimension) then
               r_elements = r_elements + 1
               r_incremented = .true.
            else
               r_incremented = .false.
            endif
!
            if (    p_elements*q_elements*r_elements*req3_min    &
                  + p_elements*q_elements           *req2_pq_min &
                  + p_elements*r_elements           *req2_pr_min &
                  + q_elements*r_elements           *req2_qr_min &
                  + p_elements                      *req1_p_min  &
                  + q_elements                      *req1_q_min  &
                  + r_elements                      *req1_r_min  &
                  + req0_tot .ge. mem%available) then
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
         batch_p%max_length = int(p_elements)
         batch_q%max_length = int(q_elements)
         batch_r%max_length = int(r_elements)
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
         if (batch_p%index_dimension .eq. batch_q%index_dimension .and. &
             batch_p%index_dimension .eq. batch_r%index_dimension) then
!
            batch_q%max_length  = batch_p%max_length
            batch_q%num_batches = batch_p%num_batches
!
            batch_r%max_length  = batch_p%max_length
            batch_r%num_batches = batch_p%num_batches
!
         else
!
            call batch_q%force_batch()
            call batch_r%force_batch()
!
         end if
!
      endif
!
      max_memory_usage = req0_tot + &
                         req1_p_min*int(batch_p%max_length, kind=i64)  +   &
                         req1_q_min*int(batch_q%max_length, kind=i64)  +   &
                         req1_r_min*int(batch_r%max_length, kind=i64)  +   &
                         req2_pq_min*int(batch_p%max_length, kind=i64)     &
                                    *int(batch_q%max_length, kind=i64) +   &
                         req2_pr_min*int(batch_p%max_length, kind=i64)     &
                                    *int(batch_r%max_length, kind=i64) +   &
                         req2_qr_min*int(batch_q%max_length, kind=i64)     &
                                    *int(batch_r%max_length, kind=i64) +   &
                         req3_min   *int(batch_p%max_length, kind=i64)     &
                                    *int(batch_q%max_length, kind=i64)     &
                                    *int(batch_r%max_length, kind=i64)
!
      if (batch_p%num_batches .eq. 1 .and. &
          batch_q%num_batches .eq. 1 .and. &
          batch_r%num_batches .eq. 1) then
!
         if (present(req_single_batch)) then
!
            max_memory_usage = int(req_single_batch*e_size, kind=i64)
!
         endif
!
      endif
!
      if (any([batch_p%num_batches, batch_q%num_batches, batch_r%num_batches] > 1)) then
         call output%printf('v', 'Batching in (a0)', chars=[tag])
      end if
!
      call mem%initialize_batching_tracker(max_memory_usage, tag)
!
   end subroutine batch_setup_3_memory_manager
!
!
end module memory_manager_class
