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
program eT_program
!
!!
!!  eT - a coupled cluster program
!!  Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!
   use kinds
   use file_class
   use disk_manager_class
   use memory_manager_class
   use libint_initialization
   use timings_class
!
   use omp_lib
!
   implicit none
!
!  Other variables
!
   integer :: n_threads = 1
!
!  Timer object
!
   type(timings) :: eT_timer
!
!  Prepare input, output and timing file
!
   call output%init('eT.out')
   call disk%open_file(output, 'write', 'rewind')
!
   call input%init('eT.inp')
   call disk%open_file(input, 'read')
!
   call timing%init('timing.out')
   call disk%open_file(timing, 'write', 'rewind')
!
!  Prepare memory manager and disk manager
!
   call mem%prepare()
   call disk%prepare()
!
   call eT_timer%init("Total time in eT")
   call eT_timer%start()
!
!  Print program banner
!
   write(output%unit,'(///t24,a)')                       'eT - a coupled cluster program '
   write(output%unit,'(t8,a)')         'Original authors: Sarai D. Folkestad, Eirik F. Kjønstad, and Henrik Koch'
   flush(output%unit)
!
   write(output%unit,'(/t3, a)')     '----------------------------------------------------------------------------------'
   write(output%unit,'(t4, a, a)')    'Author:                ','Contribution(s):'
   write(output%unit,'(t3, a)')      '----------------------------------------------------------------------------------'
   write(output%unit,'(t4, a, a)')    'Sarai D. Folkestad     ','Program design, HF, CCS, CC2, CCSD, Libint-interface,'
   write(output%unit,'(t4, a, a)')    '                       ','Cholesky decomposition, Davidson-tool, CVS'
   write(output%unit,'(t4, a, a)')    'Linda Goletto          ','CC2'
   write(output%unit,'(t4, a, a)')    'Eirik F. Kjønstad      ','Program design, HF, UHF, CCS, CC2, CCSD, DIIS-tool,'
   write(output%unit,'(t4, a, a)')    '                       ','Cholesky decomposition, Libint-interface, Davidson-tool'
   write(output%unit,'(t4, a, a)')    'Rolf H. Myhre          ','CC3, Runtest-interface, Launch script'
   write(output%unit,'(t4, a, a)')    'Alexander Paul         ','CC2, CC3'
   write(output%unit,'(t4, a, a)')    'Andreas Skeidsvoll     ','MP2'
   write(output%unit,'(t3,a)')       '----------------------------------------------------------------------------------'
   write(output%unit,'(t4,a/)')       'Other contributors: A. Balbi, M. Scavino'
   flush(output%unit)
!
!$   n_threads = omp_get_max_threads()
!
   if (n_threads .eq. 1) then
!
      write(output%unit,'(t3,a,i0,a/)')   'Running on ',n_threads, ' OMP thread'
!
   else
!
      write(output%unit,'(t3,a,i0,a/)')   'Running on ',n_threads, ' OMP threads'
!
   endif
!
   call input%check_for_errors()
!
   call initialize_libint()
!
!  Hartree-Fock calculation
!
   if (input%requested_reference_calculation()) call reference_calculation()
!!
!!  Coupled cluster calculation
!!
!   if (input%requested_cc_calculation()) call cc_calculation()
!!
   call finalize_libint()
!
   call eT_timer%freeze()
   call eT_timer%switch_off()
!
   write(output%unit, '(/t3,a)') 'eT terminated successfully!'
!
   call disk%close_file(output)
   call disk%close_file(input)
   call disk%close_file(timing)
!
   !call mem%cleanup()
   !call disk%cleanup()
!
end program eT_program

subroutine reference_calculation()
!!
!! Reference calculation
!! Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
!! Directs the reference state calculation for eT
!!
   use hf_class
   use uhf_class
   use hf_engine_class
!
   implicit none
!
!  Possible reference wavefunctions   
!
   type(hf), allocatable  :: hf_wf
   type(uhf), allocatable :: uhf_wf
!
!  Engine
!
   type(hf_engine)   :: ref_engine
!
!  Other variables
!
   character(len=21) :: reference_wf
!
   reference_wf = input%get_reference_wf()
!
   if (trim(reference_wf) == 'hf') then
!
      allocate(hf_wf)
!
      call ref_engine%ignite(hf_wf)    
!
      deallocate(hf_wf)
!
   elseif (trim(reference_wf) == 'uhf') then
!
      allocate(uhf_wf)
!
      call ref_engine%ignite(uhf_wf)  
!
      deallocate(uhf_wf)
!
   else
!
      call output%error_msg('did not recognize the reference wavefunction ' // trim(reference_wf) //'.')
!
   endif
!
end subroutine reference_calculation