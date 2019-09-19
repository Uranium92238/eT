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
module lowmem_cc2_class
!
!!
!!    Low-memory coupled cluster singles and perturbative doubles (CC2) class module
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,
!!    Linda Goletto and Alexander Paul, 2018
!!
!
   use ccs_class
!
   implicit none
!
   type, extends(ccs) :: lowmem_cc2
!
   contains
!
      procedure :: construct_omega  => construct_omega_lowmem_cc2
!
      procedure :: omega_cc2_a1     => omega_cc2_a1_lowmem_cc2
      procedure :: omega_cc2_b1     => omega_cc2_b1_lowmem_cc2
      procedure :: omega_cc2_c1     => omega_cc2_c1_lowmem_cc2
!
      procedure :: calculate_energy => calculate_energy_lowmem_cc2
!
      procedure :: construct_excited_state_equation => construct_excited_state_equation_lowmem_cc2
!
      procedure :: effective_jacobian_transformation => effective_jacobian_transformation_lowmem_cc2
!
      procedure :: jacobian_cc2_a1 => jacobian_cc2_a1_lowmem_cc2
      procedure :: jacobian_cc2_b1 => jacobian_cc2_b1_lowmem_cc2
!
      procedure :: effective_jacobian_cc2_a1 => effective_jacobian_cc2_a1_lowmem_cc2
      procedure :: effective_jacobian_cc2_b1 => effective_jacobian_cc2_b1_lowmem_cc2
      procedure :: effective_jacobian_cc2_c1 => effective_jacobian_cc2_c1_lowmem_cc2
      procedure :: effective_jacobian_cc2_d1 => effective_jacobian_cc2_d1_lowmem_cc2
      procedure :: effective_jacobian_cc2_e1 => effective_jacobian_cc2_e1_lowmem_cc2
      procedure :: effective_jacobian_cc2_f1 => effective_jacobian_cc2_f1_lowmem_cc2
!
      procedure :: effective_jacobian_transpose_transformation => effective_jacobian_transpose_transformation_lowmem_cc2
      procedure :: jacobian_transpose_ccs_b1 => jacobian_transpose_ccs_b1_lowmem_cc2
      procedure :: jacobian_transpose_cc2_a1 => jacobian_transpose_cc2_a1_lowmem_cc2
      procedure :: jacobian_transpose_cc2_b1 => jacobian_transpose_cc2_b1_lowmem_cc2
      procedure :: effective_jacobian_transpose_cc2_a1   => effective_jacobian_transpose_cc2_a1_lowmem_cc2
      procedure :: effective_jacobian_transpose_cc2_b1   => effective_jacobian_transpose_cc2_b1_lowmem_cc2
      procedure :: effective_jacobian_transpose_cc2_c1   => effective_jacobian_transpose_cc2_c1_lowmem_cc2
      procedure :: effective_jacobian_transpose_cc2_d1   => effective_jacobian_transpose_cc2_d1_lowmem_cc2
      procedure :: effective_jacobian_transpose_cc2_e1   => effective_jacobian_transpose_cc2_e1_lowmem_cc2
      procedure :: effective_jacobian_transpose_cc2_f1   => effective_jacobian_transpose_cc2_f1_lowmem_cc2
!
   end type lowmem_cc2
!
!
   interface
!
      include "omega_lowmem_cc2_interface.F90"
      include "jacobian_lowmem_cc2_interface.F90"
      include "jacobian_transpose_lowmem_cc2_interface.F90"
      include "zop_lowmem_cc2_interface.F90"
!
   end interface
!
!
   interface lowmem_cc2
!
      procedure :: new_lowmem_cc2
!
   end interface lowmem_cc2
!
!
contains
!
!
   function new_lowmem_cc2(system) result(wf)
!!
!!    New lowmem CC2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      type(lowmem_cc2) :: wf
!
      class(molecular_system), target, intent(in) :: system 
!
      wf%name_ = 'low memory cc2'
!
      wf%system => system
!
      call wf%read_hf()
!
      call wf%initialize_files()
!
      call wf%initialize_orbital_coefficients()
      call wf%initialize_orbital_energies()
!
      call wf%read_orbital_coefficients()
      call wf%read_orbital_energies()
!
      wf%bath_orbital = .false.
      wf%frozen_core = .false.
      wf%cvs = .false.
!
      call wf%read_settings()
!
      if (wf%bath_orbital) call wf%make_bath_orbital()
      if (wf%frozen_core) call wf%remove_core_orbitals()
!
      wf%n_t1            = (wf%n_o)*(wf%n_v)
      wf%n_gs_amplitudes = wf%n_t1
      wf%n_es_amplitudes = wf%n_t1
!
      call wf%write_cc_restart()
!
      call wf%initialize_fock_ij()
      call wf%initialize_fock_ia()
      call wf%initialize_fock_ai()
      call wf%initialize_fock_ab()
!
   end function new_lowmem_cc2
!
!
   subroutine construct_excited_state_equation_lowmem_cc2(wf, X, R, w, r_or_l)
!!
!!    Construct excited state equation
!!    Written by Eirik F. Kjønstad, Dec 2018
!!
!!    Constructs R = AX - wX, where w = X^T A X and norm(X) = sqrt(X^T X) = 1
!!
!!    Note I: we assume that X is normalized. If it is not,
!!    please normalize before calling the routine.
!!
!!    Note II: this routine constructs the excited state equation
!!    for standard CC models and the effective (!) excited state
!!    equation in perturbative models. In the lowmem_CC2 routine, for
!!    instance, X and R will be n_o*n_v vectors and A(w) will
!!    depend on the excitation energy w. See, e.g., Weigend and
!!    Hättig's RI-lowmem_CC2 paper for more on this topic. This means
!!    that w should be the previous w-value when entering the
!!    routine (so that A(w)X may be constructed approximately)
!!    in perturbative models.
!!
!!    Note III: the routine is used by the DIIS excited state solver.
!!
      implicit none
!
      class(lowmem_cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in)    :: X
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: R
!
      character(len=*), intent(in) :: r_or_l
!
      real(dp), intent(inout) :: w
!
      real(dp), dimension(:), allocatable :: X_copy
!
      real(dp) :: ddot
!
!     Construct residual based on previous excitation energy w
!
      call mem%alloc(X_copy, wf%n_es_amplitudes)
!
      call dcopy(wf%n_es_amplitudes, X, 1, X_copy, 1)
!
      if (r_or_l .eq. "right") then
         call wf%effective_jacobian_transformation(w, X_copy) ! X_copy <- AX
      elseif (r_or_l .eq. "left") then
         call wf%effective_jacobian_transpose_transformation(w, X_copy) ! X_copy <- A^TX
      endif

!
      call dcopy(wf%n_es_amplitudes, X_copy, 1, R, 1)
      call daxpy(wf%n_es_amplitudes, -w, X, 1, R, 1)
!
!     Update excitation energy w
!
      w = ddot(wf%n_es_amplitudes, X, 1, X_copy, 1)
      call mem%dealloc(X_copy, wf%n_es_amplitudes)
!
   end subroutine construct_excited_state_equation_lowmem_cc2
!
!
end module lowmem_cc2_class
