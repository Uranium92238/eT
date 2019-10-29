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
   use global_in 
   use global_out
   use timings_class, only : timings, timing
   use memory_manager_class, only : mem
   use libint_initialization, only : initialize_libint, finalize_libint
   use molecular_system_class, only : molecular_system
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
   eT_timer = timings("Total time in eT", pl='minimal')
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
   write(output%unit,'(t4, a, a)')    'Josefine H. Andersen   ','first order properties'
   write(output%unit,'(t4, a, a)')    'Sarai D. Folkestad     ','program design, HF, CCS, CC2, CCSD, libint-interface,'
   write(output%unit,'(t4, a, a)')    '                       ','Cholesky decomposition, Davidson-tool, CVS, DIIS-tool'
   write(output%unit,'(t4, a, a)')    '                       ','zeroth order properties, first order properties, IP'
   write(output%unit,'(t4, a, a)')    '                       ','frozen core, MLCC2, MLHF'
   write(output%unit,'(t4, a, a)')    'Tommaso Giovannini     ','QM/MM, Basis Sets'
   write(output%unit,'(t4, a, a)')    'Linda Goletto          ','CC2, MLHF'
   write(output%unit,'(t4, a, a)')    'Tor S. Haugland        ','SAD'
   write(output%unit,'(t4, a, a)')    'Anders Hutcheson       ','frozen HF orbitals'
   write(output%unit,'(t4, a, a)')    'Ida-Marie Høyvik       ','MLHF, frozen HF orbitals'
   write(output%unit,'(t4, a, a)')    'Eirik F. Kjønstad      ','program design, HF, UHF, CCS, CC2, CCSD, DIIS-tool,'
   write(output%unit,'(t4, a, a)')    '                       ','Cholesky decomposition, Libint-interface, Davidson-tool'
   write(output%unit,'(t4, a, a)')    '                       ','zeroth order properties, first order properties, SAD   '
   write(output%unit,'(t4, a, a)')    '                       ','BFGS-tool                                              '
   write(output%unit,'(t4, a, a)')    'Rolf H. Myhre          ','CC3, Runtest-interface, launch script, file system'
   write(output%unit,'(t4, a, a)')    'Alexander Paul         ','CC2, CC3, first order properties'
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
!  Set print level in output and timing files
   call set_global_print_levels()
!
!  Prepare memory manager and disk manager
!
   call mem%prepare()
!
   call initialize_libint()
!
!  Prepare molecular system 
!
   system = molecular_system()
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
   call system%cleanup()
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
   use global_in,  only: input 
   use global_out, only: output
!
   use hf_class, only: hf 
   use uhf_class, only: uhf 
   use mlhf_class, only: mlhf 
!
   use reference_engine_class, only: reference_engine 
   use hf_geoopt_engine_class, only: hf_geoopt_engine 
!
   implicit none
!
   type(molecular_system) :: system
!
   class(hf), allocatable  :: ref_wf
!
   class(reference_engine), allocatable :: ref_engine
!
   character(len=25) :: ref_wf_name
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
   elseif (trim(ref_wf_name) == 'mlhf') then
!
      ref_wf = mlhf(system)
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
      ref_engine = hf_geoopt_engine()
!
   else 
!
      ref_engine = reference_engine()
!
   endif 
!
   call ref_engine%ignite(ref_wf)
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
   use global_in,  only: input
   use global_out, only: output 
!
   use molecular_system_class, only: molecular_system
!
   use ccs_class, only: ccs 
   use cc2_class, only: cc2 
   use lowmem_cc2_class, only: lowmem_cc2
   use ccsd_class, only: ccsd
   use cc3_class, only: cc3
   use mp2_class, only: mp2 
   use mlcc2_class, only: mlcc2
!
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
   class(gs_engine), allocatable :: cc_engine 
!
   character(len=25) :: cc_wf_name
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
      case ('mlcc2')
!
         cc_wf = mlcc2(system)
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
!
!
subroutine set_global_print_levels()
!!
!! Set global print levels
!! Written by Rolf H. Myhre, Oct. 2019
!!
!! Reads and sets the global print levels for the output file
!! and the timing file from input.
!!
   use global_out, only: output 
   use global_in, only: input
   use timings_class, only : timing
!
   character(len=200) :: print_level
!
!  Set default
   print_level = 'normal' 
!
!  Overwrite print_level if keyword is present
   call input%get_keyword_in_section('output print level', 'print', print_level)
!
!  This is the only place this routine is allowed to be called
   call output%set_global_print_level(print_level)
!
!  Repeat for timing file
!  Set default
   print_level = 'normal' 
!
!  Overwrite print_level if keyword is present
   call input%get_keyword_in_section('timing print level', 'print', print_level)
!
!  This is the only place this routine is allowed to be called
   call timing%set_global_print_level(print_level)
!
end subroutine set_global_print_levels
!
