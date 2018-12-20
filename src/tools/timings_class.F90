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
!!    If you wish to print times repeatedly without switching 
!!    off the clock, you can do so by calling print_times().
!!    This call should only be used when the clock is freezed. 
!!    It will give you the currently accumulated time. 
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
!!    as shown in the timing output file.
!!
      class(timings) :: timer 
!
      character(len=*) :: tag 
!
      timer%tag = tag
!
      timer%elapsed_wall_time = zero 
      timer%elapsed_cpu_time  = zero 
      timer%wall_time_start   = zero 
      timer%wall_time_end     = zero 
      timer%cpu_time_start    = zero 
      timer%cpu_time_end      = zero 
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
!
      real(dp) :: omp_get_wtime
!
      call cpu_time(timer%cpu_time_start)
      timer%wall_time_start = omp_get_wtime()
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
!
      real(dp) :: omp_get_wtime
!
      call cpu_time(timer%cpu_time_end)
      timer%wall_time_end = omp_get_wtime()
!
      timer%elapsed_cpu_time  = timer%elapsed_cpu_time + (timer%cpu_time_end - timer%cpu_time_start)
      timer%elapsed_wall_time = timer%elapsed_wall_time + (timer%wall_time_end - timer%wall_time_start)
!
   end subroutine freeze_timings
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
      call timer%init(timer%tag)
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
      write(timing%unit, '(t3,a17,f19.12)')  'wall time (sec): ', timer%elapsed_wall_time
      write(timing%unit, '(t3,a17,f19.12)')  'cpu time (sec):  ', timer%elapsed_cpu_time
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
