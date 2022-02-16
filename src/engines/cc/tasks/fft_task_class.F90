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
module fft_task_class
!
!!
!! Fast fourier transform (FFT) task class
!! Written by Eirik F. Kjønstad and Alexander C. Paul, Feb 2022
!!
!
   use ccs_class,     only: ccs
   use cc_task_class, only: cc_task
!
   implicit none
!
   type, extends(cc_task) :: fft_task
!
      character(len=:), allocatable :: f, file_tag
!
   contains
!
      procedure, public :: execute &
                        => execute_fft_task
!
   end type fft_task
!
!
   interface fft_task
!
      procedure :: new_fft_task
!
   end interface fft_task
!
!
contains
!
!
   function new_fft_task(f, file_tag) result(this)
!!
!!    New fft task
!!    Written by Alexander C. Paul, Feb 2022
!!
      implicit none
!
      type(fft_task) :: this
!
      character(len=*), intent(in) :: f, file_tag
!
      this%f = f
      this%file_tag = file_tag
!
      this%name_ = 'Fast fourier transform of ' // this%f
!
   end function new_fft_task
!
!
   subroutine execute_fft_task(this, wf)
!!
!!    Execute
!!    Written by Alexander C. Paul, Feb 2022
!!
      use complex_fft_class, only: complex_fft
      use warning_suppressor, only: do_nothing
!
      implicit none
!
      class(fft_task), intent(inout) :: this
!
      class(ccs), intent(inout), target :: wf
!
      type(complex_fft), allocatable :: fft_solver
!
      call this%print_header()
      call this%start_timer()
!
      call do_nothing(wf)
!
      fft_solver = complex_fft(this%f, this%file_tag)
!
      call fft_solver%run()
      call fft_solver%cleanup()
!
      deallocate(fft_solver)
!
      call this%end_timer()
!
   end subroutine execute_fft_task
!
!
end module fft_task_class
