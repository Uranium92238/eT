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
module timings_class
!
!!
!!    Timings class module
!!    Written by Eirik F. Kjønstad, Dec 2018
!!
!!    Class to handle a timing in the program. Usage:
!!
!!       type(timings) :: A1_timer
!!
!!       A1_timer = new_timer('name')
!!       call A1_timer%turn_on()
!!
!!       ... do stuff ...
!!
!!       call A1_timer%turn_off()
!!
!!    You can freeze and turn_on the clock repeatedly, e.g. in a loop:
!!
!!    do (...)
!!
!!       call A1_timer%turn_on()
!!
!!       ... stuff that we want to time 
!!
!!       call A1_timer%freeze()
!!
!!       ... stuff that we don't want to time 
!!
!!    enddo
!!
!!    call A1_timer%turn_off()
!!
!!    The object will automatically print the time elapsed to
!!    the timing file in eT when turn_off is called, but you can 
!!    also request the elapsed time after the timer has been 
!!    turned off (e.g. for prints to the main output file):
!!
!!       call A1_timer%turn_off() 
!!
!!       wall_time = A1_timer%get_elapsed_time('wall')
!!       cpu_time  = A1_timer%get_elapsed_time('cpu')
!!
!!    A timer that has been turned off may be reused,
!!    though the tag will be the same (indistinguishable
!!    in output). Iteration timers in solvers are an 
!!    example where this is useful. In that case, do a 
!!    reset to zero out the timings of a turned off timer:
!!
!!    do while (.not. converged)
!!
!!       call iteration_timer%turn_on()
!!
!!       ... do stuff ...
!!
!!       call iteration_timer%turn_off()
!!       call iteration_timer%reset()
!!
!!    enddo
!!
!
   use parameters
   use global_files
   use output_file_class
!
   implicit none
!
   type timings 
!
      private 
      character(len=100) :: tag 
!
      logical :: on
!
      real(dp) :: elapsed_wall_time
      real(dp) :: elapsed_cpu_time 
!
      real(dp) :: wall_time_start
      real(dp) :: cpu_time_start
!
      real(dp) :: wall_time_end
      real(dp) :: cpu_time_end
!
   contains 
!
      procedure :: turn_on                => turn_on_timings
      procedure :: freeze                 => freeze_timings
      procedure :: turn_off               => turn_off_timings
!
      procedure :: reset                  => reset_timings
!
      procedure :: get_elapsed_time       => get_elapsed_time_timings
!
      procedure, private :: print_times   => print_times_timings
!
   end type timings 
!
!
   interface timings
!
      procedure :: new_timer
!
   end interface timings 
!
!
contains
!
!
   function new_timer(tag) result(timer)
!!
!!    Timer constructor routine 
!!    Written by Eirik F. Kjønstad, Dec 2018 
!!
!!    Initializes timer. Tag is the name of the timer,
!!    as shown in the timing output file when turn_off()
!!    is called.
!!
      type(timings) :: timer 
!
      character(len=*), intent(in) :: tag 
!
!     Set name & then set all times to zero 
!
      timer%tag = tag
      call timer%reset()
!
   end function new_timer
!
!
   subroutine turn_on_timings(timer)
!!
!!    Turn on  
!!    Written by Eirik F. Kjønstad, Dec 2018 
!!
!!    Turns on the wall and CPU clocks.
!!
      implicit none 
!
      class(timings) :: timer 
      integer :: counter, c_rate
!
      call cpu_time(timer%cpu_time_start)
      call system_clock(count=counter, count_rate=c_rate)
      timer%wall_time_start = real(counter,dp)/c_rate
!
      timer%on = .true.
!
   end subroutine turn_on_timings
!
!
   subroutine freeze_timings(timer)
!!
!!    Freeze 
!!    Written by Eirik F. Kjønstad, Dec 2018 
!!
!!    Freezes the wall and CPU clocks and updates 
!!    the accumulated elapsed times. 
!!
      implicit none 
!
      class(timings) :: timer 
      integer :: counter, c_rate
!
      call cpu_time(timer%cpu_time_end)
      call system_clock(count=counter, count_rate=c_rate)
      timer%wall_time_end = real(counter,dp)/c_rate
!
      timer%elapsed_cpu_time  = timer%elapsed_cpu_time + (timer%cpu_time_end - timer%cpu_time_start)
      timer%elapsed_wall_time = timer%elapsed_wall_time + (timer%wall_time_end - timer%wall_time_start)
!
      timer%on = .false.
!
   end subroutine freeze_timings
!
!
   subroutine reset_timings(timer)
!!
!!    Reset 
!!    Written by Eirik F. Kjønstad, Dec 2018 
!!
!!    Sets all times to zero.
!!
      implicit none 
!
      class(timings), intent(inout) :: timer 
!
      timer%elapsed_wall_time    = zero 
      timer%elapsed_cpu_time     = zero 
      timer%wall_time_start      = zero 
      timer%wall_time_end        = zero 
      timer%cpu_time_start       = zero 
      timer%cpu_time_end         = zero 
!
   end subroutine reset_timings
!
!
   subroutine turn_off_timings(timer)
!!
!!    Turn off 
!!    Written by Eirik F. Kjønstad, Dec 2018 
!!
!!    Prints the wall and CPU times to the timing file,
!!    stopping the clock if it was not frozen beforehand.
!!
      implicit none 
!
      class(timings), intent(inout) :: timer 
!
      if (timer%on) call timer%freeze() 
      call timer%print_times()
!
   end subroutine turn_off_timings
!
!
   subroutine print_times_timings(timer)
!!
!!    Print times 
!!    Written by Eirik F. Kjønstad, Dec 2018 
!!
!!    Prints the wall and CPU times to the timing file.
!!
      implicit none 
!
      class(timings), intent(in) :: timer 
!
      write(timing%unit, '(/t3,a)') timer%tag
      write(timing%unit, '(t3,a17,f20.2)')  'wall time (sec): ', timer%elapsed_wall_time
      write(timing%unit, '(t3,a17,f20.2)')  'cpu time (sec):  ', timer%elapsed_cpu_time
!
   end subroutine print_times_timings
!
!
   function get_elapsed_time_timings(timer, what) result(elapsed)
!!
!!    Get elapsed time 
!!    Written by Eirik F. Kjønstad, Dec 2018 
!!
!!    Returns the elapsed time.
!!
!!       what == 'wall' => returns elapsed wall time 
!!       what == 'cpu'  => returns elapsed CPU time 
!!
      implicit none 
!
      class(timings), intent(in) :: timer 
!
      character(len=*), intent(in) :: what
!
      real(dp) :: elapsed 
!
      if (timer%on) call output%error_msg("Don't ask for elapsed time when the clock is on.")
!
      elapsed = zero
      if (trim(what) == 'wall') then 
!
         elapsed = timer%elapsed_wall_time
!
      elseif (trim(what) == 'cpu') then 
!
         elapsed = timer%elapsed_cpu_time
!
      else
!
         call output%error_msg('Did not recognize requested time in timings object ' // timer%tag)
!
      endif 
!
   end function get_elapsed_time_timings
!
!
end module timings_class
