module cc_response_solver_settings_class
!!
!!    CC response solver settings
!!    Written by Eirik Kjøntad, Oct 2022
!!
   use parameters
!
   implicit none 
!
   type :: cc_response_solver_settings
!
      integer  :: max_iterations
      logical  :: records_in_memory
      real(dp) :: threshold  
!
   end type cc_response_solver_settings
!
   interface cc_response_solver_settings
!
      procedure :: new_cc_response_solver_settings
!
   end interface cc_response_solver_settings
!
contains 
!
!
   function new_cc_response_solver_settings() result(this)
!!
!!    New
!!    Written by Eirik Kjøntad, Oct 2022
!!
      use global_in, only: input 
      use global_out, only: output 
!
      implicit none 
!
      type(cc_response_solver_settings) :: this 
!
      character(len=200) :: storage
!
      this%threshold = 1.0d-3 
      call input%get_keyword('threshold', 'solver cc response', this%threshold)
!
      this%max_iterations = 100
      call input%get_keyword('max iterations', 'solver cc response', this%max_iterations)
!
      storage = 'disk'
      call input%get_keyword('storage', 'solver cc response', storage)
!
      if (trim(storage) == 'memory') then
!
         this%records_in_memory = .true.
!
      elseif (trim(storage) == 'disk') then
!
         this%records_in_memory = .false.
!
      else
!
         call output%error_msg('Could not recognize keyword storage in solver cc response: ' // &
                                 trim(storage))
!
      endif
!
   end function new_cc_response_solver_settings
!
!
end module cc_response_solver_settings_class