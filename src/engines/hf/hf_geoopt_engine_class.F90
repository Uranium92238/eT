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
module hf_geoopt_engine_class
!
!!
!! HF geometry optimization engine class
!! Written by Sarai D. Folkestad, Eirik F. Kjønstad, and Alexander C. Paul, 2018-2022
!!
!
   use hf_engine_class, only: hf_engine
   use hf_class,        only: hf
!
   use hf_geoopt_task_class, only: hf_geoopt_task
!
   implicit none
!
   type, extends(hf_engine) :: hf_geoopt_engine
!
      type(hf_geoopt_task), allocatable, private :: geometry_optimization
!
   contains
!
      procedure, public :: ignite => ignite_hf_geoopt_engine
!
   end type hf_geoopt_engine
!
contains
!
!
   subroutine ignite_hf_geoopt_engine(this, wf)
!!
!!    Ignite
!!    Written by Alexander C. Paul, May 2022
!!
      implicit none
!
      class(hf_geoopt_engine), intent(inout) :: this
      class(hf), intent(inout) :: wf
!
      this%geometry_optimization = hf_geoopt_task()
      call this%geometry_optimization%execute(wf)
!
   end subroutine ignite_hf_geoopt_engine
!
!
   end module hf_geoopt_engine_class
   