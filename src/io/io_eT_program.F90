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
   integer function get_n_methods()
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
!
   subroutine read_cc_methods(n_methods, cc_methods)
!!
!!    Read cc methods, 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkstad, August 2018
!!
!!    Reads cc methods (and mp2 for now) and orders them such that 
!!    cc calculations will be performed in correct order.
!!
!!
      implicit none
!
      integer :: n_methods
!
      character(len=40), dimension(n_methods) :: cc_methods
!
      character(len=100) :: line
!
      integer :: n_methods_total, count_cc_methods, i
!
      count_cc_methods = 0
!
!     Read for ccs
!
      call move_to_section('method', n_methods_total)
!
      do i = 1, n_methods_total
!
         read(input%unit, '(a100)') line
         line = remove_preceding_blanks(line)   
!
         if (trim(line) == 'ccs' .or. trim(line) == 'CCS' .or. trim(line) == 'Ccs')  then 
!
            count_cc_methods = count_cc_methods + 1
            cc_methods(count_cc_methods) = 'ccs'
            exit
!
         endif
!
      enddo
!
!     Read for mp2
!
      call move_to_section('method', n_methods_total)
!
      do i = 1, n_methods_total
!
         read(input%unit, '(a100)') line
         line = remove_preceding_blanks(line)   
!
         if (trim(line) == 'mp2' .or. trim(line) == 'MP2' .or. trim(line) == 'Mp2')  then 
!
            count_cc_methods = count_cc_methods + 1
            cc_methods(count_cc_methods) = 'mp2'
            exit
!
         endif
!
      enddo
!
!     Read for cc2
!
      call move_to_section('method', n_methods_total)
!
      do i = 1, n_methods_total
!
         read(input%unit, '(a100)') line
         line = remove_preceding_blanks(line)   
!
         if (trim(line) == 'cc2' .or. trim(line) == 'CC2' .or. trim(line) == 'Cc2')  then 
!
            count_cc_methods = count_cc_methods + 1
            cc_methods(count_cc_methods) = 'cc2'
            exit
!
         endif
!
      enddo
!
!     Read for lowmem-cc2
!
      call move_to_section('method', n_methods_total)
!
      do i = 1, n_methods_total
!
         read(input%unit, '(a100)') line
         line = remove_preceding_blanks(line)   
!
         if (trim(line) == 'lowmem-cc2')  then 
!
            count_cc_methods = count_cc_methods + 1
            cc_methods(count_cc_methods) = 'lowmem-cc2'
            exit
!
         endif
!
      enddo
!
!     Read for ccsd
!
      call move_to_section('method', n_methods_total)
!
      do i = 1, n_methods_total
!
         read(input%unit, '(a100)') line
         line = remove_preceding_blanks(line)   
!
         if (trim(line) == 'ccsd' .or. trim(line) == 'CCSD' .or. trim(line) == 'Ccsd')  then 
!
            count_cc_methods = count_cc_methods + 1
            cc_methods(count_cc_methods) = 'ccsd'
            exit
!
         endif
!
      enddo
!
      if (count_cc_methods .ne. n_methods) call output%error_msg('reading cc methods.')
!
   end subroutine read_cc_methods
!
!
   subroutine select_engine(engine)
!!
!!    Select engine, 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkstad, August 2018
!!
!!    Select cc engine
!!
      implicit none
!
      character(len=40) :: engine
!
!     Read for multipliers
!
      if (requested_task('multipliers') .or. requested_task('Multipliers')) then
!
         engine = 'multipliers'
         return
!
      endif
!
!     Read for excited state
!
      if (requested_task('excited state') .or. requested_task('Excited state')) then
!
         engine = 'excited state'
         return
!
      endif
!
!     Read for ground state
!
      if (requested_task('ground state') .or. requested_task('Ground state')) then
!
         engine = 'ground state'
         return
!
      endif 
!
      call output%error_msg('No cc tasks recognized.')
!
   end subroutine select_engine
!
end module io_eT_program
