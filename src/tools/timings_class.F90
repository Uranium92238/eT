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
!!       call A1_timer%init('name')
!!       call A1_timer%start()
!!
!!       ... do stuff ...
!!
!!       call A1_timer%freeze()
!!       call A1_timer%switch_off()
!!
!!    You can freeze and start the clock repeatedly, e.g. in a loop:
!!
!!    do (...)
!!
!!       call A1_timer%start()
!!
!!       ... stuff that we want to time 
!!
!!       call A1_timer%freeze()
!!
!!       ... stuff that we don't want to time 
!!
!!    enddo
!!
!!    call A1_timer%switch_off()
!!
!!    The object will automatically print the time elapsed to
!!    the timing file in eT when swtich_off is called, but you can 
!!    also request the elapsed time by calling 
!! 
!!       wall_time = A1_timer%get_elapsed_time('wall')
!!       cpu_time  = A1_timer%get_elapsed_time('cpu')
!!
!!    Note that a switched off clock is reset! This means that 
!!    you must collect the timings before you switch it off.
!!    Otherwise they will just equal zero.
!!
!!    A timer that has been switched off may be reused,
!!    though the tag will be the same (indistinguishable
!!    in output).  
!!
!
   use file_class
   use parameters
!
   implicit none
!
   type timings 
!
      character(len=100) :: tag 
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
      procedure :: init             => init_timings
!
      procedure :: start            => start_timings
      procedure :: freeze           => freeze_timings
      procedure :: reset            => reset_timings
!
      procedure :: switch_off       => switch_off_timings
!
      procedure :: print_times      => print_times_timings
      procedure :: get_elapsed_time => get_elapsed_time_timings
!
   end type timings 
!
contains
!
!
   subroutine init_timings(timer,tag)
!!
!!    Init 
!!    Written by Eirik F. Kjønstad, Dec 2018 
!!
!!    Initializes timer. Tag is the name of the timer,
!!    as shown in the timing output file when switch_off()
!!    is called.
!!
      class(timings) :: timer 
!
      character(len=*) :: tag 
!
!     Set name & then set all times to zero 
!
      timer%tag = tag
      call timer%reset()
!
   end subroutine init_timings
!
!
   subroutine start_timings(timer)
!!
!!    Start 
!!    Written by Eirik F. Kjønstad, Dec 2018 
!!
!!    Starts the wall and CPU clocks.
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
   end subroutine start_timings
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
      timer%elapsed_wall_time = zero 
      timer%elapsed_cpu_time  = zero 
      timer%wall_time_start   = zero 
      timer%wall_time_end     = zero 
      timer%cpu_time_start    = zero 
      timer%cpu_time_end      = zero 
!
   end subroutine reset_timings
!
!
   subroutine switch_off_timings(timer)
!!
!!    Switch off 
!!    Written by Eirik F. Kjønstad, Dec 2018 
!!
!!    Prints the wall and CPU times to the timing file,
!!    then resets the timer.
!!
      implicit none 
!
      class(timings), intent(inout) :: timer 
!
!     Print & reset 
!
      call timer%print_times()
      call timer%reset()
!
   end subroutine switch_off_timings
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
      write(timing%unit, '(t3,a17,f12.5)')  'wall time (sec): ', timer%elapsed_wall_time
      write(timing%unit, '(t3,a17,f12.5)')  'cpu time (sec):  ', timer%elapsed_cpu_time
!
   end subroutine print_times_timings
!
!
   real(dp) function get_elapsed_time_timings(timer, what)
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
      character(len=*) :: what
!
      get_elapsed_time_timings = zero
      if (what == 'wall') then 
!
         get_elapsed_time_timings = timer%elapsed_wall_time
!
      elseif (what == 'cpu') then 
!
         get_elapsed_time_timings = timer%elapsed_cpu_time
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
