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
module conjugate_gradient_tool_class
!
!!
!!    Conjugate gradient solver tool class module
!!    Written by Åsmund H. Tveten and Eirik F. Kjønstad, 2019
!!
!
   use kinds
   use parameters
!
   use file_class
   use disk_manager_class
   use memory_manager_class
!
!
   type :: conjugate_gradient_tool
!
      private 
      integer  :: iteration 
      integer  :: n_parameters
      real(dp) :: beta ! Fletcher-Reeves CG constant
!
      real(dp), dimension(:), allocatable, private :: previous_gradient
!
   contains
!
      procedure, public :: get_next_direction => get_next_direction_conjugate_gradient_tool
!
      procedure, private :: compute_beta      => compute_beta_conjugate_gradient_tool
!
   end type conjugate_gradient_tool
!
!
   interface conjugate_gradient_tool 
!
      procedure :: new_conjugate_gradient_tool
!
   end interface conjugate_gradient_tool 
!
!
contains
!
!
   function new_conjugate_gradient_tool(n_parameters) result(cg)
!!
!!    New conjugate gradient tool
!!    Written by Åsmund H. Tveten and Eirik F. Kjønstad, 2019
!!
      type(conjugate_gradient_tool) :: cg
!
      integer, intent(in) :: n_parameters
!
      cg%n_parameters = n_parameters
      cg%iteration    = 0
!
      call mem%alloc(cg%previous_gradient, cg%n_parameters)
!
   end function new_conjugate_gradient_tool
!
!
   subroutine get_next_direction_conjugate_gradient_tool(cg, gradient, descent_direction)
!!
!!    Get next direction 
!!    Written by Åsmund H. Tveten and Eirik F. Kjønstad, 2019
!!
!!    Gradient:          current gradient 
!!    Descent_direction: next direction according to the CG algorithm   
!!
      implicit none 
!
      class(conjugate_gradient_tool), intent(inout) :: cg 
!
      real(dp), dimension(cg%n_parameters), intent(in)      :: gradient  
      real(dp), dimension(cg%n_parameters), intent(inout)   :: descent_direction  
!
      cg%iteration = cg%iteration + 1
!
      if (cg%iteration .eq. 1) then 
!
         descent_direction = -gradient
!
      else
!
         call cg%compute_beta(gradient)
         descent_direction = cg%beta*descent_direction - gradient 
         write(output%unit, *) 'beta is : ', cg%beta
!
      endif 
!
      cg%previous_gradient = gradient
!
   end subroutine get_next_direction_conjugate_gradient_tool
!
!
   subroutine compute_beta_conjugate_gradient_tool(cg, gradient)
!!
!!    Compute beta 
!!    Written by Åsmund H. Tveten and Eirik F. Kjønstad, 2019
!!
!!    Computes beta as the Flecter-Reeves constant.
!!
      implicit none 
!
      class(conjugate_gradient_tool), intent(inout) :: cg 
!
      real(dp), dimension(cg%n_parameters), intent(in) :: gradient 
!
      real(dp) :: gradient2, previous_gradient2, ddot  
!
      real(dp), dimension(:), allocatable :: temp 
!
      call mem%alloc(temp, cg%n_parameters)
      temp = gradient - cg%previous_gradient
!
      gradient2 = ddot(cg%n_parameters, gradient, 1, temp, 1)
      previous_gradient2 = ddot(cg%n_parameters, cg%previous_gradient, 1, cg%previous_gradient, 1)
!
      cg%beta = max(gradient2/previous_gradient2, zero)
!
      call mem%dealloc(temp, cg%n_parameters)  
!
   end subroutine compute_beta_conjugate_gradient_tool
!
!
end module conjugate_gradient_tool_class
