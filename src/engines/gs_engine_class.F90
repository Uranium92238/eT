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
module gs_engine_class
!!
!!    Coupled cluster ground state engine class module 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
   use abstract_engine_class
   use ccs_class
   use eri_cd_class
!
   use diis_cc_gs_class
!
   type, extends(abstract_engine) :: gs_engine 
!
   contains 
!
      procedure :: prepare => prepare_gs_engine
      procedure :: run     => run_gs_engine
      procedure :: cleanup => cleanup_gs_engine
!
!
   end type gs_engine 
!
contains
!
   subroutine prepare_gs_engine(engine)
!!
!!    Prepare 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
      implicit none 
!
      class(gs_engine) :: engine 
!
      engine%name_ = 'Ground state engine'
!
   end subroutine prepare_gs_engine
!
!
   subroutine run_gs_engine(engine, wf)
!!
!!    Run 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
      implicit none 
!
      class(gs_engine) :: engine 
!
      class(ccs) :: wf
!
      type(eri_cd), allocatable     :: eri_chol_solver
      type(diis_cc_gs), allocatable :: cc_gs_solver 
!
      write(output%unit, '(/t3,a,a)') '- Running ', trim(engine%name_)
!
!     Cholesky decoposition 
!
      allocate(eri_chol_solver)
!
      call eri_chol_solver%prepare(wf%system)
      call eri_chol_solver%run(wf%system)
!
      call eri_chol_solver%cholesky_vecs_diagonal_test(wf%system)
!
      call eri_chol_solver%construct_mo_cholesky_vecs(wf%system, wf%n_mo, wf%orbital_coefficients)
!
      call wf%integrals%prepare(eri_chol_solver%n_cholesky, wf%n_o, wf%n_v)
!
      call eri_chol_solver%cleanup()
      deallocate(eri_chol_solver)
!
!     Ground state solution (avoid starting solver if there are no equations to solve)
!
      if (trim(wf%name_) == 'MP2') then 
!
         call wf%calculate_energy()
         call wf%print_wavefunction_summary()
!
      else 
!
         allocate(cc_gs_solver)
!
         call cc_gs_solver%prepare(wf)
         call cc_gs_solver%run(wf)
         call cc_gs_solver%cleanup(wf)
!
         deallocate(cc_gs_solver)
!
      endif
!
   end subroutine run_gs_engine
!
!
   subroutine cleanup_gs_engine(engine)
!!
!!    Cleanup 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018 
!!
      implicit none 
!
      class(gs_engine) :: engine 
!
      write(output%unit, '(/t3,a,a)') '- Cleaning up ', trim(engine%name_)
!
   end subroutine cleanup_gs_engine
!
!
end module gs_engine_class
