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
   use ccs_class
   use abstract_engine_class
!
   type, extends(abstract_engine) :: gs_engine 
!
      character(len=200) :: multipliers_algorithm
!
      character(len=200) :: gs_algorithm  
!
   contains 
!
      procedure :: prepare                               => prepare_gs_engine
      procedure :: run                                   => run_gs_engine
!
      procedure :: do_ground_state                       => do_ground_state_gs_engine 
!
      procedure :: do_multipliers                        => do_multipliers_gs_engine 
!
      procedure, nopass :: calculate_dipole_moment       => calculate_dipole_moment_gs_engine
!
      procedure, nopass :: calculate_quadrupole_moment   => calculate_quadrupole_moment_gs_engine
      procedure, nopass :: remove_trace                  => remove_trace_gs_engine
!
      procedure :: read_settings                         => read_settings_gs_engine
! 
      procedure :: read_gs_settings                      => read_gs_settings_gs_engine 
!
   end type gs_engine 
!
!
contains
!
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
      engine%name_                 = 'ground state coupled cluster'
!
      engine%timer = timings(trim(engine%name_) // ' engine')
      call engine%timer%turn_on()
!
      engine%multipliers_algorithm = 'davidson'
      engine%gs_algorithm          = 'diis'
!
      call engine%read_settings()
!
      engine%tasks = [character(len=150) ::                                                                          &
            'Cholesky decomposition of the ERI-matrix',                                                              &
            'Calculation of the ground state amplitudes ('//trim(engine%gs_algorithm)//'-algorithm)', &
            'Calculation of the ground state energy']
!
      engine%description = 'Calculates the ground state CC wavefunction | CC > = exp(T) | R >'
!
   end subroutine prepare_gs_engine
!
!
   subroutine read_settings_gs_engine(engine)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019 
!!
      implicit none 
!
      class(gs_engine) :: engine 
!
      call engine%read_gs_settings()
!
   end subroutine read_settings_gs_engine
!
!
   subroutine read_gs_settings_gs_engine(engine)
!!
!!    Read GS engine settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019 
!!
      implicit none 
!
      class(gs_engine) :: engine 
!
      call input%get_keyword_in_section('algorithm', 'solver cc multipliers', engine%multipliers_algorithm)
      call input%get_keyword_in_section('algorithm', 'solver cc gs', engine%gs_algorithm )
!
   end subroutine read_gs_settings_gs_engine
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
      class(ccs)       :: wf
!
      call engine%print_banner(wf)
!
!     Cholesky decoposition of the electron repulsion integrals 
!
      call engine%do_cholesky(wf, wf%orbital_coefficients)
!
!     Determine ground state
!
      call engine%do_ground_state(wf)
!
   end subroutine run_gs_engine
!
!
   subroutine do_ground_state_gs_engine(engine, wf)
!!
!!    Do ground state   
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019 
!!
!!    Performs a ground state coupled cluster calculation. Determines 
!!    the amplitudes, i.e. the (right) ground state | CC > = e^T | R >.
!!
!!    After this routine, the wavefunction has the cluster amplitudes 
!!    stored in memory.     
!!
      use diis_cc_gs_class
      use newton_raphson_cc_gs_class
!
      implicit none 
!
      class(gs_engine), intent(in) :: engine 
!
      class(ccs), intent(inout) :: wf 
!
      type(diis_cc_gs), allocatable :: diis_solver 
      type(newton_raphson_cc_gs), allocatable :: newton_raphson_solver 
!
      if (trim(wf%name_) == 'mp2') then 
!
         call wf%integrals%write_t1_cholesky(wf%t1)
         call wf%calculate_energy()
!
         call wf%print_wavefunction_summary()
!
      elseif (trim(engine%gs_algorithm) == 'diis') then 
!
         allocate(diis_solver)
!
         call diis_solver%prepare(wf)
         call diis_solver%run(wf)
         call diis_solver%cleanup(wf)
!
         deallocate(diis_solver)
!
      elseif (trim(engine%gs_algorithm) == 'newton-raphson') then 
!
         allocate(newton_raphson_solver)
!
         call newton_raphson_solver%prepare(wf)
         call newton_raphson_solver%run(wf)
         call newton_raphson_solver%cleanup(wf)
!
         deallocate(newton_raphson_solver)
!
      endif 
!
   end subroutine do_ground_state_gs_engine
!
!
   subroutine do_multipliers_gs_engine(engine, wf)
!!
!!    Do multipliers  
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019 
!!
!!    Performs a (left) ground state coupled cluster calculation. Determines 
!!    the amplitudes, i.e. the (left) ground state < Λ |  = (< R | + sum_mu tbar_mu < mu |) e^-T.
!!
!!    After this routine, the wavefunction has the cluster amplitudes 
!!    and the multipliers stored in memory.   
!!
      use diis_cc_multipliers_class
      use davidson_cc_multipliers_class
!
      implicit none 
!
      class(gs_engine), intent(in) :: engine 
!
      class(ccs), intent(inout) :: wf 
!
      type(diis_cc_multipliers), allocatable     :: diis_solver
      type(davidson_cc_multipliers), allocatable :: davidson_solver
!
      if (trim(engine%multipliers_algorithm) == 'davidson') then 
!
         allocate(davidson_solver)
!
         call davidson_solver%prepare(wf)
         call davidson_solver%run(wf)
         call davidson_solver%cleanup(wf)
!
         deallocate(davidson_solver)
!
      elseif (trim(engine%multipliers_algorithm) == 'diis') then 
!
         allocate(diis_solver)
!
         call diis_solver%prepare(wf)
         call diis_solver%run(wf)
         call diis_solver%cleanup(wf)
!
         deallocate(diis_solver)
!
      else
!
         call output%error_msg('Could not recognize multipliers algorithm ' // trim(engine%multipliers_algorithm))
!
      endif 
!
   end subroutine do_multipliers_gs_engine
!
!
   subroutine calculate_dipole_moment_gs_engine(wf, electronic, nuclear, total)
!!
!!    Calculate dipole moment 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019 
!!
!!    Determines < mu > = < Λ | mu | CC > for the one-electron operator mu. 
!!    The routine assumes that the multipliers and amplitudes have been determined 
!!    and that the one-electron density is calculated and in memory. 
!!
      implicit none 
!
      class(ccs), intent(in) :: wf 
!
      real(dp), dimension(3), intent(out) :: electronic
      real(dp), dimension(3), intent(out) :: nuclear
      real(dp), dimension(3), intent(out) :: total 
!
      integer :: k 
!
      real(dp), dimension(:,:,:), allocatable :: mu_pqk
!
!     Get the integrals mu_pqk for components k = 1, 2, 3 in the T1-transformed basis  
!  
      call mem%alloc(mu_pqk, wf%n_mo, wf%n_mo, 3)
      call wf%construct_mu(mu_pqk)
!
!     Get electronic expectation value contribution 
!
      do k = 1, 3 
!
         electronic(k) = wf%calculate_expectation_value(mu_pqk(:,:,k))
!
      enddo 
!
      call mem%dealloc(mu_pqk, wf%n_mo, wf%n_mo, 3)
!
!     Get nuclear expectation value contribution, then sum the two
!
      call wf%system%get_nuclear_dipole(nuclear)
!
      total = electronic + nuclear 
!
   end subroutine calculate_dipole_moment_gs_engine
!
!
   subroutine calculate_quadrupole_moment_gs_engine(wf, electronic, nuclear, total)
!!
!!    Calculate quadrupole moment 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019 
!!
!!    Determines < q > = < Λ | q | CC > for the one-electron operator q. 
!!    The routine assumes that the multipliers and amplitudes have been determined 
!!    and that the one-electron density is calculated and in memory. 
!!
      implicit none 
!
      class(ccs), intent(in) :: wf 
!
      real(dp), dimension(6), intent(out) :: electronic
      real(dp), dimension(6), intent(out) :: nuclear
      real(dp), dimension(6), intent(out) :: total 
!
      integer :: k 
!
      real(dp), dimension(:,:,:), allocatable :: q_pqk
!
!     Get the integrals q_pqk for components k = 1, 2, ..., 6 in the T1-transformed basis  
!  
      call mem%alloc(q_pqk, wf%n_mo, wf%n_mo, 6)
      call wf%construct_q(q_pqk)
!
!     Get electronic expectation value contribution 
!
      do k = 1, 6
!
         electronic(k) = wf%calculate_expectation_value(q_pqk(:,:,k))
!
      enddo 
!
      call mem%dealloc(q_pqk, wf%n_mo, wf%n_mo, 6)
!
!     Get nuclear expectation value contribution, then sum the two
!
      call wf%system%get_nuclear_quadrupole(nuclear)
!
      total = electronic + nuclear 
!
   end subroutine calculate_quadrupole_moment_gs_engine
!
!
   subroutine remove_trace_gs_engine(M)
!!
!!    Remove trace 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
!!    The assumption here is that M is a 2-tensor ordered as xx, xy, xz, yy, yz, and zz,
!!    where the other elements of the tensor are given by symmetry, such as for the quadrupole 
!!    moment. Thus, this routine can be called after a call to "calculate quadrupole moment"
!!    to make the moment trace-free.
!!
      implicit none 
!
      real(dp), dimension(6), intent(inout) :: M
!
      real(dp) :: trace_
!
      trace_ = M(1) + M(4) + M(6)
!
      M(1) = (three*M(1) - trace_)/two
      M(4) = (three*M(4) - trace_)/two
      M(6) = (three*M(6) - trace_)/two
!
      M(2) = (three*M(2))/two
      M(3) = (three*M(3))/two
      M(5) = (three*M(5))/two
!
   end subroutine remove_trace_gs_engine
!
!
end module gs_engine_class
