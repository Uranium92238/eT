module io_eT_program
!
!!
!!    IO eT-program module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkstad, August 2018
!!
!
   use kinds
   use file_class
   use io_utilities
!
contains
!
   integer(i15) function get_n_methods()
!!
!!
!!
      implicit none
!
      character(len=100) :: line
!
      rewind(input%unit)
!
      get_n_methods = 0
!
      do 
!
         read(input%unit, '(a100)') line
         line = remove_preceding_blanks(line)
!
         if (trim(line) .eq. 'end geometry') call output%error_msg('method section not detected.')
!
         if (trim(line) == 'method') then
!
            do
!
               read(input%unit, '(a100)') line
               line = remove_preceding_blanks(line)
!
               if (trim(line) .eq. 'end method')  return
!
               get_n_methods = get_n_methods + 1
!
            enddo
!
         endif
!
      enddo
!
   end function get_n_methods
!
!
   logical function requested_method(method)
!!
!!
!!
      implicit none
!
      character(len=*) :: method
!
      character(len=100) :: line
!
      rewind(input%unit)
!
      requested_method = .false.
!
      do 
!
         read(input%unit, '(a100)') line
         line = remove_preceding_blanks(line)
!
         if (trim(line) .eq. 'end geometry') call output%error_msg('method section not detected.')
!
         if (trim(line) == 'method') then
!
            do
!
               read(input%unit, '(a100)') line
               line = remove_preceding_blanks(line)
!
               if (trim(line) .eq. 'end method') return
!
               if (trim(line) .eq. method) then
                  requested_method = .true.
                  return
               endif
!
            enddo
!
         endif
!
      enddo
!
   end function requested_method
!
!
   logical function requested_task(task)
!!
!!
!!
      implicit none
!
      character(len=*) :: task
!
      character(len=100) :: line
!
      rewind(input%unit)
!
      requested_task = .false.
!
      do 
!
         read(input%unit, '(a100)') line
         line = remove_preceding_blanks(line)
!
         if (trim(line) .eq. 'end geometry') call output%error_msg('task section not detected.')
!
         if (trim(line) == 'do') then
!
            do
!
               read(input%unit, '(a100)') line
               line = remove_preceding_blanks(line)
!
               if (trim(line) .eq. 'end do') return
!
               if (trim(line) .eq. task) then
                  requested_task = .true.
                  return
               endif
!
            enddo
!
         endif
!
      enddo
!
   end function requested_task
!
end module io_eT_program