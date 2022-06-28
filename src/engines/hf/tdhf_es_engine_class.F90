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
module tdhf_es_engine_class
!
!!
!! TDHF excited state engine class
!! Written by Sarai D. Folkestad, Alexander C. Paul, 2018-2022
!!
!
   use hf_engine_class, only: hf_engine
   use hf_class,        only: hf
!
   use tdhf_es_task_class, only: tdhf_es_task
!
   implicit none
!
   type, extends(hf_engine) :: tdhf_es_engine
!
      type(tdhf_es_task), allocatable, private :: tdhf
!
   contains
!
      procedure, public :: ignite => ignite_hf_es_engine
!
   end type tdhf_es_engine
!
contains
!
!
   subroutine ignite_hf_es_engine(this, wf)
!!
!!    Ignite
!!    Written by Alexander C. Paul, May 2022
!!
      use hf_gs_engine_class, only: hf_gs_engine
!
      implicit none
!
      class(tdhf_es_engine), intent(inout) :: this
      class(hf),           intent(inout) :: wf
!
      type(hf_gs_engine), allocatable :: gs_engine
!
      gs_engine = hf_gs_engine()
      call gs_engine%ignite(wf)
!
      this%tdhf = tdhf_es_task()
      call this%tdhf%execute(wf)
!
   end subroutine ignite_hf_es_engine
!
!
   end module tdhf_es_engine_class
