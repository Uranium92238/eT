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
!!  Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2019
!!
!
   use kinds
   use file_class
   use disk_manager_class
   use memory_manager_class
   use libint_initialization
   use molecular_system_class
   use timings_class
!
   use omp_lib
!
   implicit none
!
   integer :: io_error
   integer :: n_threads
!
!  Molecular system object 
!
   type(molecular_system) :: system 
!
!  Timer object
!
   type(timings) :: eT_timer
!
!  Prepare input, output and timing file
!
   call output%init('eT.out')
   open(newunit=output%unit, file=output%name, access=output%access, &
      action='write', status='unknown', form=output%format, iostat=io_error)
!
   call input%init('eT.inp')
   open(newunit=input%unit, file=input%name, access=input%access, &
      action='read', status='unknown', form=input%format, iostat=io_error)
!
   call timing%init('timing.out')
   open(newunit=timing%unit, file=timing%name, access=timing%access, &
      action='write', status='unknown', form=timing%format, iostat=io_error)
!
   if (io_error /= 0) stop 'Error: could not open eT files (.inp/.out)'
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
   n_threads = 1
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
!  Prepare memory manager and disk manager
!
   call mem%prepare()
   call disk%prepare()
!
   call initialize_libint()
!
!  Prepare molecular system 
!
   call system%prepare()
!
!  Hartree-Fock calculation
!
   if (input%requested_reference_calculation()) call reference_calculation(system)
!
!  Coupled cluster calculation
!
   if (input%requested_cc_calculation()) call cc_calculation(system)
!
   call finalize_libint()
!
   call eT_timer%freeze()
   call eT_timer%switch_off()
!
   write(output%unit, '(/t3,a)') 'eT terminated successfully!'
!
   close(output%unit)
   close(input%unit)
   close(timing%unit)
!
end program eT_program
!
!
subroutine reference_calculation(system)
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
   type(molecular_system) :: system
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
      call hf_wf%prepare(system)
      call ref_engine%ignite(hf_wf)
      call hf_wf%cleanup() 
!
      deallocate(hf_wf)
!
   elseif (trim(reference_wf) == 'uhf') then
!
      allocate(uhf_wf)
!
      call uhf_wf%prepare(system)
      call ref_engine%ignite(uhf_wf)
      call uhf_wf%cleanup()
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
!
!
subroutine cc_calculation(system)
!!
!! Coupled cluster calculation
!! Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
!! Directs the coupled cluster calculation for eT
!!
   use ccs_class
   use cc2_class
   use lowmem_cc2_class
   use cc3_class
   use mp2_class
!
   use gs_engine_class
   use es_engine_class
   use zop_engine_class
!
   implicit none
!
   type(molecular_system) :: system
!
!  Possible coupled cluster wavefunctions   
!
   type(ccs), target          :: ccs_wf
   type(cc2), target          :: cc2_wf
   type(lowmem_cc2), target   :: lowmem_cc2_wf
   type(ccsd),target          :: ccsd_wf
   type(cc3), target          :: cc3_wf
   type(mp2), target          :: mp2_wf
!
   class(ccs), pointer :: cc_wf
!
!  Possible engines
!
   type(gs_engine)   :: gs_cc_engine
   type(es_engine)   :: es_cc_engine
   type(zop_engine)  :: zop_cc_engine
!
!  Other variables
!
   character(len=21) :: cc_wf_name
!
   cc_wf_name = input%get_cc_wf()
!
   select case (trim(cc_wf_name))
!
      case ('ccs')
!
         cc_wf => ccs_wf
!
      case ('cc2')
!
         cc_wf => cc2_wf
!
      case ('lowmem-cc2')
!
         cc_wf => lowmem_cc2_wf
!
      case ('ccsd')
!
         cc_wf => ccsd_wf
!
      case ('cc3')
!
         cc_wf => cc3_wf
!
      case ('mp2')
!
         cc_wf => mp2_wf
!
      case default
!
         call output%error_msg('could not recognize CC method ' // trim(cc_wf_name) // '.')
!
   end select
!
   if (input%requested_keyword_in_section('excited state', 'do')) then
!
      call cc_wf%prepare(system)
      call es_cc_engine%ignite(cc_wf)
      call cc_wf%cleanup()   
!
   elseif (input%requested_keyword_in_section('zop','do')) then 
!
      call cc_wf%prepare(system)
      call zop_cc_engine%ignite(cc_wf)
      call cc_wf%cleanup()
!
   elseif (input%requested_keyword_in_section('ground state', 'do')) then
!
      call cc_wf%prepare(system)
      call gs_cc_engine%ignite(cc_wf) 
      call cc_wf%cleanup()  
!
   else
!
      call output%error_msg('could not recognize coupled cluster task.')
!
   endif
!
end subroutine cc_calculation
