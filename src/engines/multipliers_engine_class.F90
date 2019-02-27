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
module multipliers_engine_class
!!
!!    Coupled cluster multiplier engine class module 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
   use abstract_engine_class
   use ccs_class
   use eri_cd_solver_class
   use davidson_cc_es_solver_class
   use davidson_cc_ip_solver_class
   use davidson_cvs_cc_es_solver_class
   use diis_cc_gs_solver_class
   use davidson_cc_multipliers_class
   use diis_cc_multipliers_class
!
   type, extends(abstract_engine) :: multipliers_engine 
!
      character(len=100) :: algorithm
!
   contains 
!
      procedure :: prepare          => prepare_multipliers_engine
      procedure :: run              => run_multipliers_engine
      procedure :: cleanup          => cleanup_multipliers_engine
      procedure :: read_algorithm   => read_algorithm_multipliers_engine
!
   end type multipliers_engine 
!
contains
!
   subroutine prepare_multipliers_engine(engine)
!!
!!    Prepare 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
      implicit none 
!
      class(multipliers_engine) :: engine 
!
      engine%name_ = 'Multipliers engine'
      engine%algorithm = 'davidson'
!
      call engine%read_algorithm()
!
   end subroutine prepare_multipliers_engine
!
!
   subroutine run_multipliers_engine(engine, wf)
!!
!!    Run 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
      implicit none 
!
      class(multipliers_engine) :: engine 
!
      class(ccs) :: wf
!
      type(eri_cd_solver), allocatable                   :: eri_chol_solver
      type(diis_cc_gs_solver), allocatable               :: cc_gs_solver
!
      type(davidson_cc_multipliers), allocatable   :: cc_multipliers_davidson
      type(diis_cc_multipliers), allocatable       :: cc_multipliers_diis
!
      write(output%unit, '(/t3,a,a)') '- Running ', trim(engine%name_)
!
!     Cholesky decomposition 
!
      allocate(eri_chol_solver)
!
      call eri_chol_solver%prepare(wf%system)
      call eri_chol_solver%run(wf%system)
!
      call eri_chol_solver%cholesky_vecs_diagonal_test(wf%system)
      call eri_chol_solver%full_test_cholesky_vecs(wf%system)
!
      call eri_chol_solver%construct_mo_cholesky_vecs(wf%system, wf%n_mo, wf%orbital_coefficients)
!
      call wf%integrals%prepare(eri_chol_solver%n_cholesky, wf%n_o, wf%n_v)
!
      call eri_chol_solver%cleanup()
      deallocate(eri_chol_solver)
!
!     Ground state solution 
!
      allocate(cc_gs_solver)
!
      call cc_gs_solver%prepare(wf)
      call cc_gs_solver%run(wf)
      call cc_gs_solver%cleanup(wf)
!
      deallocate(cc_gs_solver)
!
      call wf%integrals%write_t1_cholesky(wf%t1)
!
!     Multiplier equation
!
      if (engine%algorithm .ne. 'davidson' .and. engine%algorithm .ne. 'diis') then
!
         call output%error_msg('Could not recognize algorithm for multiplier equation.')
!
      elseif (wf%name_ == 'cc2' .or. engine%algorithm == 'diis') then
!
         allocate(cc_multipliers_diis)
!
         call cc_multipliers_diis%prepare(wf)
         call cc_multipliers_diis%run(wf)
         call cc_multipliers_diis%cleanup(wf)
!
         deallocate(cc_multipliers_diis)
!
      elseif (engine%algorithm == 'davidson') then
!
         allocate(cc_multipliers_davidson)
!
         call cc_multipliers_davidson%prepare(wf)
         call cc_multipliers_davidson%run(wf)
         call cc_multipliers_davidson%cleanup(wf)
!
         deallocate(cc_multipliers_davidson)
!
      endif
!
   end subroutine run_multipliers_engine
!
!
   subroutine cleanup_multipliers_engine(engine)
!!
!!    Cleanup 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
      implicit none 
!
      class(multipliers_engine) :: engine 
!
      write(output%unit, '(/t3,a,a)') '- Cleaning up ', trim(engine%name_)
!
   end subroutine cleanup_multipliers_engine
!
!
   subroutine read_algorithm_multipliers_engine(engine)
!!
!!    Read algorithm
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
      implicit none
!
      class(multipliers_engine), intent(inout) :: engine 
!
      character(len=100) :: line
!
      integer :: i, n_records
!
      if (requested_section('multipliers')) then
         call move_to_section('multipliers', n_records)
!
         do i = 1, n_records
!
            read(input%unit, '(a100)') line
            line = remove_preceding_blanks(line)
!
            if (line(1:10) == 'algorithm:') then
!
               engine%algorithm = line(11:100)
               engine%algorithm = remove_preceding_blanks(engine%algorithm)
               return
!
            endif
!
         enddo
      endif
!
   end subroutine read_algorithm_multipliers_engine
!
!
end module multipliers_engine_class
