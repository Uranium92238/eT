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
   use direct_file_class, only : direct_file
   use sequential_file_class, only : sequential_file
   use global_files, only : input
   use timings_class, only : timings
   use disk_manager_class, only : disk
   use memory_manager_class
   use libint_initialization
   use molecular_system_class
!
   use omp_lib
!
   implicit none
!
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
   output = output_file('eT.out')
   call output%open_()
!
   input = input_file('eT.inp')
   call input%open_()
!
   timing = output_file('timing.out')
   call timing%open_()
!
   eT_timer = timings("Total time in eT")
   call eT_timer%turn_on()
!
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
   write(output%unit,'(t4, a, a)')    'Josefine H. Andersen   ','First order properties'
   write(output%unit,'(t4, a, a)')    'Sarai D. Folkestad     ','Program design, HF, CCS, CC2, CCSD, Libint-interface,'
   write(output%unit,'(t4, a, a)')    '                       ','Cholesky decomposition, Davidson-tool, CVS, DIIS-tool'
   write(output%unit,'(t4, a, a)')    '                       ','Zeroth order properties, First order properties'
   write(output%unit,'(t4, a, a)')    'Linda Goletto          ','CC2'
   write(output%unit,'(t4, a, a)')    'Eirik F. Kjønstad      ','Program design, HF, UHF, CCS, CC2, CCSD, DIIS-tool,'
   write(output%unit,'(t4, a, a)')    '                       ','Cholesky decomposition, Libint-interface, Davidson-tool'
   write(output%unit,'(t4, a, a)')    '                       ','Zeroth order properties, First order properties,       '
   write(output%unit,'(t4, a, a)')    '                       ','BFGS-tool                                              '
   write(output%unit,'(t4, a, a)')    'Rolf H. Myhre          ','CC3, Runtest-interface, Launch script'
   write(output%unit,'(t4, a, a)')    'Alexander Paul         ','CC2, CC3, First order properties'
   write(output%unit,'(t4, a, a)')    'Andreas Skeidsvoll     ','MP2'
   write(output%unit,'(t4, a, a)')    'Åsmund H. Tveten       ','HF'
   write(output%unit,'(t3,a)')       '----------------------------------------------------------------------------------'
   write(output%unit,'(t4,a/)')       'Other contributors: A. Balbi, M. Scavino'
   call output%flush_()   
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
   call eT_timer%turn_off()
!
   call mem%check_for_leak()
!
   write(output%unit, '(/t3,a)') 'eT terminated successfully!'
!
   call output%close_()
   call input%close_()
   call timing%close_()
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
   use molecular_system_class, only: molecular_system
!
   use global_files, only: input 
   use output_file_class, only: output
!
   use hf_class, only: hf 
   use uhf_class, only: uhf 
   use hf_engine_class, only: hf_engine 
!
   use hf_geoopt_engine_class, only: hf_geoopt_engine 
!
   implicit none
!
   type(molecular_system) :: system
!
   class(hf), allocatable  :: ref_wf
!
   type(hf_engine)         :: ref_engine
   type(hf_geoopt_engine)  :: ref_geoopt_engine 
!
   character(len=21) :: ref_wf_name
!
   ref_wf_name = input%get_reference_wf()
!
   if (trim(ref_wf_name) == 'hf') then
!
      ref_wf = hf(system)
!
   elseif (trim(ref_wf_name) == 'uhf') then
!
      ref_wf = uhf(system)
!
   else
!
      call output%error_msg('did not recognize the reference wavefunction ' &
                                    // trim(ref_wf_name) //'.')
!
   endif
!
   if (input%requested_keyword_in_section('ground state geoopt', 'do')) then 
!
      ref_geoopt_engine = hf_geoopt_engine()
      call ref_geoopt_engine%ignite(ref_wf)
!
   else 
!
      ref_engine = hf_engine()
      call ref_engine%ignite(ref_wf)
!
   endif 
!
   call ref_wf%cleanup()
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
   use global_files, only: input
   use output_file_class, only: output 
!
   use molecular_system_class, only: molecular_system
!
   use ccs_class, only: ccs 
   use cc2_class, only: cc2 
   use lowmem_cc2_class, only: lowmem_cc2
   use ccsd_class, only: ccsd
   use cc3_class, only: cc3
   use mp2_class, only: mp2 
!
   use abstract_engine_class, only: abstract_engine
   use gs_engine_class, only: gs_engine
   use es_engine_class, only: es_engine
   use zop_engine_class, only: zop_engine 
   use fop_engine_class, only: fop_engine 
!
   implicit none
!
   type(molecular_system) :: system
!
   class(ccs), allocatable :: cc_wf
   class(abstract_engine), allocatable :: cc_engine 
!
   character(len=21) :: cc_wf_name
!
   if (.not. input%requested_reference_calculation()) &
      call output%error_msg('to run CC calculation reference wavefunction must be specified.')
!
   cc_wf_name = input%get_cc_wf()
!
   select case (trim(cc_wf_name))
!
      case ('ccs')
!
         cc_wf = ccs(system)
!
      case ('cc2')
!
         cc_wf = cc2(system)
!
      case ('lowmem-cc2')
!
         cc_wf = lowmem_cc2(system)
!
      case ('ccsd')
!
         cc_wf = ccsd(system)
!
      case ('cc3')
!
         cc_wf = cc3(system)
!
      case ('mp2')
!
         cc_wf = mp2(system)
!
      case default
!
         call output%error_msg('could not recognize CC method ' // trim(cc_wf_name) // '.')
!
   end select
!
   if (input%requested_keyword_in_section('fop', 'do')) then
!
      cc_engine = fop_engine()
!
   elseif (input%requested_keyword_in_section('excited state', 'do')) then
!
      cc_engine = es_engine()
!
   elseif (input%requested_keyword_in_section('zop', 'do')) then 
!
      cc_engine = zop_engine()
!
   elseif (input%requested_keyword_in_section('ground state', 'do')) then
!
      cc_engine = gs_engine()
!
   else
!
      call output%error_msg('could not recognize coupled cluster task.')
!
   endif
!
   call cc_engine%ignite(cc_wf)
   call cc_wf%cleanup()
!
end subroutine cc_calculation
