!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2021 the authors of eT
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
module precondition_tool_class
!
!!
!!    Precondition tool class module
!!    Written by Eirik F. Kjønstad, 2019
!! 
!!    Vectors can be preconditioned by X as follows:
!!
!!       Initialization:
!!
!!          precondition_tool = precondition_tool(X, dim_X)
!!
!!       Then you can use
!!
!!          call precondition_tool%do(R) -> R(I) = R(I)/X(I)
!!
!!       or 
!!           
!!          call precondition_tool%do(R, shift=alpha) -> R(I) = R(I)/(X(I) - alpha).
!!
!!    Class was made to replace, and is based on, precondition functionality coded 
!!    in several solvers and solver tools originally written by Sarai D. Folkestad 
!!    and Eirik F. Kjønstad. 
!!
!
   use kinds 
   use parameters
   use memory_manager_class, only: mem
   use global_out, only: output
!
   implicit none
!
   type :: precondition_tool
!
      integer :: dim_
      real(dp), dimension(:), allocatable :: preconditioner   
!
   contains
!
      procedure :: do_ => do_precondition_tool 
!
      procedure :: destruct_precondition_vector => destruct_precondition_vector_precondition_tool
!
   end type precondition_tool
!
!
   interface precondition_tool 
!
      procedure :: new_precondition_tool
!
   end interface precondition_tool 
!
!
contains 
!
!
   function new_precondition_tool(preconditioner, dim_) result(tool)
!!
!!    New precondition tool 
!!    Written by Eirik F. Kjønstad, 2019
!!
!!    preconditioner: vector X that will precondition vectors as 
!!                    R(I) -> R(I)/(X(I) - shift).
!!
!!    dim_: length of the preconditioner vector X  
!!
      implicit none 
!
      type(precondition_tool) :: tool 
!
      integer, intent(in) :: dim_ 
      real(dp), dimension(dim_), intent(in) :: preconditioner 
!
      tool%dim_ = dim_ 
!
      call mem%alloc(tool%preconditioner, tool%dim_)
      call dcopy(tool%dim_, preconditioner, 1, tool%preconditioner, 1)
!
   end function new_precondition_tool
!
!
   subroutine do_precondition_tool(tool, R, shift, prefactor)
!!
!!    Do 
!!    Written by Eirik F. Kjønstad, 2019
!!
!!    Does the precondition operation:
!!
!!       R(i) = prefactor*R(i)/(preconditioner(i) - shift).
!!
!!    The preconditioner is set by the constructor.
!!
!!    shift:     Optional real number. Default is zero. 
!!    prefactor: Optional real number. Default is one. 
!!
      implicit none
!
      class(precondition_tool), intent(in) :: tool
!
      real(dp), dimension(tool%dim_), intent(inout) :: R 
!
      real(dp), intent(in), optional :: shift 
      real(dp), intent(in), optional :: prefactor 
!
      real(dp) :: shift_local
      real(dp) :: prefactor_local
!
      integer :: i
!
!     Set optionals 
!
      if (present(shift)) then 
!
         shift_local = shift 
!
      else 
!
         shift_local = zero 
!
      endif
!
      if (present(prefactor)) then 
!
         prefactor_local = prefactor 
!
      else 
!
         prefactor_local = one
!
      endif
!
!     Precondition R 
!
!$omp parallel do private(i)
      do i = 1, tool%dim_
!
         R(i) = prefactor_local*R(i)/(tool%preconditioner(i) - shift_local)
!
      enddo  
!$omp end parallel do
!
   end subroutine do_precondition_tool
!
!
   subroutine destruct_precondition_vector_precondition_tool(tool)
!!
!!    Destruct precondition vector 
!!    Written by Eirik F. Kjønstad, 2019 
!!
      implicit none 
!
      class(precondition_tool), intent(inout) :: tool 
!
      if (allocated(tool%preconditioner)) call mem%dealloc(tool%preconditioner, tool%dim_)
!
   end subroutine destruct_precondition_vector_precondition_tool
!
!
end module precondition_tool_class
