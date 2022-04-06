!
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
module eri_cholesky_memory_class
!
!!
!! ERI Cholesky memory class
!! Written by Sarai D. Folkestad, Sep 2021
!!
!! Handles the memory storage and delivery of the
!! Cholesky vectors
!!
   use parameters
!
   use global_out,                  only : output
   use cholesky_block_list_class,   only: cholesky_block_list
   use range_class,                 only: range_
   use memory_manager_class,        only: mem
   use abstract_eri_cholesky_class, only: abstract_eri_cholesky
!
   implicit none
!
   type, extends(abstract_eri_cholesky) :: eri_cholesky_memory
!
      type(cholesky_block_list), allocatable, private ::  L_Jpq
!
   contains
!
      procedure, private :: allocate_L_Jpq
!
      procedure, public :: initialize &
                        => initialize_eri_cholesky_memory
!
      procedure, public :: get_block &
                        => get_block_eri_cholesky_memory
!
      procedure, public :: set_block &
                        => set_block_eri_cholesky_memory
!
      procedure, public :: get_memory_estimate &
                        => get_memory_estimate_eri_cholesky_memory
!
      procedure, public :: load_memory_estimate &
                        => load_memory_estimate_eri_cholesky_memory
!
      procedure, public :: load_block &
                        => load_block_eri_cholesky_memory
!
      procedure, public :: offload_block &
                        => offload_block_eri_cholesky_memory
!
      final :: cleanup
!
   end type eri_cholesky_memory
!
   interface eri_cholesky_memory
!
      procedure :: new_eri_cholesky_memory
!
   end interface eri_cholesky_memory
!
contains
!
   function new_eri_cholesky_memory() result(this)
!!
!!    New
!!    Written by Sarai D. Folkestad, Sep 2021
!!
      implicit none
!
      type(eri_cholesky_memory) :: this
!
      this%n_J       = -1
      this%n_ranges  = -1
      this%dim_      = -1
!
   end function new_eri_cholesky_memory
!
!
   subroutine initialize_eri_cholesky_memory(this, n_J, n_ranges, range_lengths)
!!
!!    Initialize
!!    Written by Sarai D. Folkestad, Sep 2021
!!
      implicit none

      class(eri_cholesky_memory), intent(inout) :: this
!
      integer, intent(in) :: n_ranges, n_J
      integer, dimension(n_ranges), intent(in) :: range_lengths
!
      call this%set_dimensions(n_J, n_ranges, range_lengths)
      call this%set_index_ranges()
!
      this%L_Jpq = cholesky_block_list(this%n_J)
      this%L_Jpq_loaded = cholesky_block_list(this%n_J)
!
      call this%allocate_L_Jpq()
!
   end subroutine initialize_eri_cholesky_memory
!
!
   subroutine allocate_L_Jpq(this)
!!
!!    Allocate L_Jpq
!!    Written by Sarai D. Folkestad, Sep 2021
!!
      implicit none
!
      class(eri_cholesky_memory), intent(inout) :: this
!
      integer :: p, q
!
      this%n_blocks = 0
!
      do p = 1, this%n_ranges
         do q = 1, this%n_ranges
!
            this%n_blocks = this%n_blocks + 1
            call this%L_Jpq%push_back(this%index_ranges(p), this%index_ranges(q))
!
         enddo
      enddo
!
   end subroutine allocate_L_Jpq
!
!
   subroutine set_block_eri_cholesky_memory(this, block_, overlap_p, overlap_q,  &
                                    L_Jpq, range_p, range_q)
!!
!!    Set block
!!    Written by Sarai D. Folkestad, Sep 2021
!!
!!    Sets a subblock of the given block, given by overlap_p and overlap_q
!!
      implicit none
!
      class(eri_cholesky_memory),                            intent(inout):: this
!
      type(range_), intent(in) :: overlap_p, overlap_q
      type(range_), intent(in) :: range_p, range_q
      integer, intent(in) :: block_
!
      real(dp), dimension(this%n_J, range_p%length, range_q%length), intent(in)   :: L_Jpq
!
      real(dp), dimension(:,:,:), pointer :: block_pointer
!
      integer :: p, q, J
!
      type(range_) :: range_r, range_s
!
      call this%get_ranges_from_block(block_, range_r, range_s)
      call this%L_Jpq%get_array(block_pointer, range_r, range_s)
!
!$omp parallel do private(p, q, J)
      do q = overlap_q%first, overlap_q%get_last()
         do p = overlap_p%first, overlap_p%get_last()
            do J = 1, this%n_J

               block_pointer(J, p - range_r%first + 1, &
                                q - range_s%first + 1) = &
                       L_Jpq(J, p - range_p%first + 1, q - range_q%first + 1)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine set_block_eri_cholesky_memory
!
!
   subroutine get_block_eri_cholesky_memory(this, block_, overlap_p, overlap_q, &
                                    L_Jpq, range_p, range_q)
!!
!!    Get block
!!    Written by Sarai D. Folkestad, Sep 2021
!!
!!    Gets a subblock of the given block, given by overlap_p and overlap_q
!!
      implicit none
!
      class(eri_cholesky_memory), intent(inout):: this
!
      type(range_), intent(in) :: overlap_p, overlap_q
      type(range_), intent(in) :: range_p, range_q
      integer,      intent(in) :: block_
!
      real(dp), dimension(this%n_J, range_p%length, range_q%length), intent(out)   :: L_Jpq
!
      integer :: p, q, J
!
      real(dp), dimension(:,:,:), pointer :: block_pointer
!
      type(range_) :: range_r, range_s
!
      call this%get_ranges_from_block(block_, range_r, range_s)
!
      call this%L_Jpq%get_array(block_pointer, range_r, range_s)
!
!$omp parallel do private(p, q, J)
      do q = overlap_q%first, overlap_q%get_last()
         do p = overlap_p%first, overlap_p%get_last()
            do J = 1, this%n_J
!
               L_Jpq(J, p - range_p%first + 1, q - range_q%first + 1) = &
                     block_pointer(J, p - range_r%first + 1, q - range_s%first + 1)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine get_block_eri_cholesky_memory
!
!
   subroutine cleanup(this)
!!
!!    Cleanup
!!    Written by Sarai D. Folkestad, Sep 2021
!!
      implicit none
!
      type(eri_cholesky_memory),  intent(inout)  :: this
!
      if (allocated(this%range_lengths)) call mem%dealloc(this%range_lengths, this%n_ranges)
      if (allocated(this%L_Jpq)) call this%L_Jpq%finalize()
 !
   end subroutine cleanup
!
!
   function get_memory_estimate_eri_cholesky_memory(this, first_p, last_p, first_q, last_q) result(memory)
!!
!!    Memory estimate
!!    Written by Sarai D. Folkestad, Sep 2021
!!
!!    Returns the memory required by the Cholesky tool
!!    to handle the indicated Cholesky vector block_
!!
!
      use warning_suppressor, only: do_nothing
!
      implicit none
!
      class(eri_cholesky_memory), intent(in) :: this
!
      integer, intent(in) :: first_p, last_p
      integer, intent(in) :: first_q, last_q
!
      integer :: memory
!
      call do_nothing(this)
      call do_nothing([first_p, last_p, first_q, last_q])
!
      memory = 0
!
   end function get_memory_estimate_eri_cholesky_memory
!
!
   function load_memory_estimate_eri_cholesky_memory(this, first_p, last_p, first_q, last_q) result(memory)
!!
!!    Memory estimate
!!    Written by Sarai D. Folkestad, Sep 2021
!!
!!    Returns the memory required by the Cholesky tool
!!    to handle the indicated Cholesky vector block_
!!
!
      use warning_suppressor, only: do_nothing
!
      implicit none
!
      class(eri_cholesky_memory), intent(in) :: this
!
      integer, intent(in) :: first_p, last_p
      integer, intent(in) :: first_q, last_q
!
      integer :: memory
!
      type(range_) :: range_p, range_q
!
      range_p = range_(first_p, last_p - first_p + 1)
      range_q = range_(first_q, last_q - first_q + 1)
!
      memory = 0
!
      if (this%is_contiguous_in_block(range_p, range_q)) return
!
      memory = range_p%length*range_q%length*this%n_J
!
   end function load_memory_estimate_eri_cholesky_memory
!
!
   subroutine load_block_eri_cholesky_memory(this, L_Jpq, first_p, last_p, first_q, last_q)
!!
!!    Load block
!!    Written by Sarai D. Folkestad, Jan 2022
!!
      implicit none
!
      class(eri_cholesky_memory), intent(inout):: this
      integer,                     intent(in)   :: first_p, last_p
      integer,                     intent(in)   :: first_q, last_q
!
         real(dp), dimension(:,:,:), pointer,   intent(out)   :: L_Jpq
!
      type(range_) :: range_p, range_q

      integer           :: block_
      type(range_)      :: range_r, range_s
!
      range_p = range_(first_p, last_p - first_p + 1)
      range_q = range_(first_q, last_q - first_q + 1)
!
      if (this%is_contiguous_in_block(range_p, range_q)) then
!
         block_ = this%get_contiguous_block_index(range_p, range_q)
         call this%get_ranges_from_block(block_, range_r, range_s)
!
         call this%L_Jpq%get_array(L_Jpq,                             &
                                  range_r, range_s,                   &
                                  range_q%first - range_s%first + 1,  &
                                  range_q%get_last() - range_s%first + 1)
!
      else
!

         call this%load_block_to_linked_list(range_p, range_q)
         call this%L_Jpq_loaded%get_array(L_Jpq, range_p, range_q)
!
      endif
!
   end subroutine load_block_eri_cholesky_memory
!
!
   subroutine offload_block_eri_cholesky_memory(this, first_p, last_p, first_q, last_q)
!!
!!    Offload block
!!    Written by Sarai D. Folkestad, Jan 2022
!!
!!    Note: in the case of the block being stored in this%L_pq_J
!!          nothing should happen when this routine is called.
!!
      implicit none
!
      class(eri_cholesky_memory), intent(inout):: this
      integer,                    intent(in)    :: first_p, last_p
      integer,                    intent(in)    :: first_q, last_q
!
      type(range_)      :: range_p, range_q
!
      range_p = range_(first_p, last_p - first_p + 1)
      range_q = range_(first_q, last_q - first_q + 1)
!
      if (this%L_Jpq_loaded%node_in_list(range_p, range_q)) then
!
         call this%offload_block_from_linked_list(range_p, range_q)
!
!     block not in either of the lists
      else if (.not. this%is_contiguous_in_block(range_p, range_q)) then
!
         call output%error_msg('cannot offload Cholesky block!')
!
      endif
!
   end subroutine offload_block_eri_cholesky_memory
!
!
end module eri_cholesky_memory_class
