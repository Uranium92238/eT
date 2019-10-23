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
module wavefunction_class
!
!!
!!    Wavefunction class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   use parameters
   use array_utilities, only : zero_array
   use global_out, only : output
   use sequential_file_class, only : sequential_file
   use memory_manager_class, only : mem
   use molecular_system_class, only : molecular_system
   use interval_class, only : interval
   use libint_initialization, only : initialize_potential_c
   use global_in, only: input
!
   implicit none
!
   type, abstract :: wavefunction
!
      character(len=40) :: name_
!
      real(dp) :: energy
!
      integer :: n_ao
      integer :: n_mo
      integer :: n_o
      integer :: n_v
!
      type(molecular_system), pointer :: system
!
      real(dp), dimension(:,:), allocatable :: orbital_coefficients
      real(dp), dimension(:), allocatable :: orbital_energies
!
      type(sequential_file) :: orbital_coefficients_file
      type(sequential_file) :: orbital_energies_file
!
!     Frozen orbital variables. Frozen orbitals are typically frozen core or frozen HF orbitals.
!
      real(dp), dimension(:,:), allocatable :: mo_fock_fc_contribution 
      real(dp), dimension(:,:), allocatable :: mo_fock_frozen_hf_contribution 
!
      type(sequential_file) :: mo_fock_fc_file, mo_fock_frozen_hf_file
!
      logical :: frozen_core
      logical :: frozen_hf_mos
!
      real(dp) :: cholesky_orbital_threshold = 1.0D-2
!
   contains
!
      procedure :: initialize_orbital_coefficients => initialize_orbital_coefficients_wavefunction
      procedure :: initialize_orbital_energies     => initialize_orbital_energies_wavefunction
!
      procedure :: destruct_orbital_coefficients   => destruct_orbital_coefficients_wavefunction
      procedure :: destruct_orbital_energies       => destruct_orbital_energies_wavefunction
!
      procedure :: get_ao_x_wx                     => get_ao_x_wx_wavefunction
      procedure :: get_ao_x_wx_1der                => get_ao_x_wx_1der_wavefunction
!
      procedure :: get_ao_h_wx                     => get_ao_h_wx_wavefunction
      procedure :: get_ao_h_wx_1der                => get_ao_h_wx_1der_wavefunction
!
      procedure :: get_ao_s_wx                     => get_ao_s_wx_wavefunction
      procedure :: get_ao_s_wx_1der                => get_ao_s_wx_1der_wavefunction
!
      procedure :: get_ao_mu_wx                    => get_ao_mu_wx_wavefunction
      procedure :: get_ao_q_wx                     => get_ao_q_wx_wavefunction
      procedure :: get_ao_v_wx                     => get_ao_v_wx_wavefunction
!
      procedure :: get_mo_mu                       => get_mo_mu_wavefunction
      procedure :: get_mo_h                        => get_mo_h_wavefunction
      procedure :: get_mo_q                        => get_mo_q_wavefunction
!
      procedure :: mo_transform                    => mo_transform_wavefunction
!
      procedure :: initialize_wavefunction_files   => initialize_wavefunction_files_wavefunction
!
      procedure :: read_orbital_coefficients       => read_orbital_coefficients_wavefunction
      procedure :: save_orbital_coefficients       => save_orbital_coefficients_wavefunction
      procedure :: read_orbital_energies           => read_orbital_energies_wavefunction
      procedure :: save_orbital_energies           => save_orbital_energies_wavefunction
!
      procedure :: is_restart_safe                 => is_restart_safe_wavefunction 
!
      procedure(gradient_function), deferred :: &
               construct_molecular_gradient        
!
      procedure :: projected_atomic_orbitals => projected_atomic_orbitals_wavefunction
      procedure :: get_orbital_overlap       => get_orbital_overlap_wavefunction
      procedure :: lovdin_orthonormalization => lovdin_orthonormalization_wavefunction
!
      procedure :: construct_ao_electrostatics                 => construct_ao_electrostatics_wavefunction       ! V_wx, E_wx, V(D), E(D)
      procedure :: update_h_wx_mm                              => update_h_wx_mm_hf
!  
      procedure :: construct_and_write_mo_cholesky             => construct_and_write_mo_cholesky_wavefunction      
!
      procedure :: construct_orbital_block_by_density_cd       => construct_orbital_block_by_density_cd_wavefunction
!
      procedure :: initialize_mo_fock_fc_contribution          => initialize_mo_fock_fc_contribution_wavefunction                
      procedure :: destruct_mo_fock_fc_contribution            => destruct_mo_fock_fc_contribution_wavefunction  
      procedure :: initialize_mo_fock_frozen_hf_contribution   => initialize_mo_fock_frozen_hf_contribution_wavefunction                
      procedure :: destruct_mo_fock_frozen_hf_contribution     => destruct_mo_fock_frozen_hf_contribution_wavefunction  
!
      procedure :: read_frozen_orbitals_settings   => read_frozen_orbitals_settings_wavefunction 
!
   end type wavefunction 
!
!
   abstract interface 
!
      subroutine gradient_function(wf, E_qk)
!!
!!       Gradient function
!!       Written by Eirik F. Kjønstad, 2019
!!
!!       Deferred routine to compute molecular derivative E_qk of the energy E.
!!
!!       q refers to x, y, z (q = 1,2,3)
!!       k refers to atom number (k = 1,2,3,...,n_atoms) 
!!
         import :: wavefunction, dp
!
         implicit none 
!
         class(wavefunction), intent(in) :: wf
!
         real(dp), dimension(3, wf%system%n_atoms), intent(inout) :: E_qk 
!
      end subroutine gradient_function
!
!
      subroutine construct_ao_x_wx(molecule, X, A, B)
!!
!!       Construct AO x_wx 
!!       Written by Eirik F. Kjønstad, 2019
!!
!!       Abstract interface for routines, belonging to molecular system,
!!       that calculates integrals X for the shell pair (A, B).
!!
         use kinds, only: dp 
         use molecular_system_class, only: molecular_system
!
         implicit none 
!
         class(molecular_system), intent(in) :: molecule
!
         integer, intent(in) :: A, B 
!
         real(dp), dimension(molecule%shell_limits(A)%length, &
                     molecule%shell_limits(B)%length), intent(out) :: X
!
      end subroutine construct_ao_x_wx
!
!
      subroutine construct_ao_x_wx_1der(molecule, X_Ax, X_Ay, X_Az, X_Bx, X_By, X_Bz, A, B)
!!
!!       Construct AO x_wx_1der 
!!       Written by Eirik F. Kjønstad, 2019
!!
!!       Abstract interface for routines, belonging to molecular system,
!!       that calculates first-derivative integral contributions X for the shell pair (A, B).
!!
!!       Ax: derivative with respect to the x coordinate of the nuclei that shell A is centered on
!!       Ay: derivative with respect to the y coordinate of the nuclei that shell A is centered on
!!       Az: derivative with respect to the z coordinate of the nuclei that shell A is centered on
!!       Bx: derivative with respect to the x coordinate of the nuclei that shell B is centered on
!!       By: derivative with respect to the y coordinate of the nuclei that shell B is centered on
!!       Bz: derivative with respect to the z coordinate of the nuclei that shell B is centered on
!!
         use kinds, only: dp 
         use molecular_system_class, only: molecular_system 
!
         implicit none 
!
         class(molecular_system), intent(in) :: molecule
!
         integer, intent(in) :: A, B 
!
         real(dp), dimension(molecule%shell_limits(A)%length, molecule%shell_limits(B)%length), intent(out) :: X_Ax
         real(dp), dimension(molecule%shell_limits(A)%length, molecule%shell_limits(B)%length), intent(out) :: X_Ay
         real(dp), dimension(molecule%shell_limits(A)%length, molecule%shell_limits(B)%length), intent(out) :: X_Az
         real(dp), dimension(molecule%shell_limits(A)%length, molecule%shell_limits(B)%length), intent(out) :: X_Bx
         real(dp), dimension(molecule%shell_limits(A)%length, molecule%shell_limits(B)%length), intent(out) :: X_By
         real(dp), dimension(molecule%shell_limits(A)%length, molecule%shell_limits(B)%length), intent(out) :: X_Bz
!
      end subroutine construct_ao_x_wx_1der
!
   end interface 
!
!
contains
!
!
   subroutine get_ao_x_wx_wavefunction(wf, construct_ao_x_AB, x_wx)
!!
!!    Get AO X 
!!    Written by Eirik F. Kjønstad, 2019 
!!
!!    Computes the AO integral matrix x_wx using  
!!    the construction routine construct_ao_x_AB. This routine
!!    constructs the x_AB contributions to x_wx for the shells A and B
!!
      implicit none 
!
      class(wavefunction), intent(in) :: wf 
!
      procedure(construct_ao_x_wx) :: construct_ao_x_AB
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(out) :: x_wx  
!
      type(interval) :: A_interval, B_interval
!
      integer :: x, y, A, B 
!
      real(dp), dimension(:,:), pointer                        :: x_AB_p 
      real(dp), dimension(wf%system%max_shell_size**2), target :: x_AB
!
!$omp parallel do private(A, B, x_AB, x_AB_p, A_interval, B_interval, x, y) schedule(static)
      do A = 1, wf%system%n_s
!
         A_interval = wf%system%shell_limits(A)
!
         do B = 1, A
!
            B_interval = wf%system%shell_limits(B)
!
            call construct_ao_x_AB(wf%system, x_AB, A, B)
!
            x_AB_p(1 : A_interval%length, 1 : B_interval%length) => x_AB(1 : A_interval%length*B_interval%length)
!
            do x = 1, A_interval%length
               do y = 1, B_interval%length
!
                  x_wx(A_interval%first - 1 + x, B_interval%first - 1 + y) = x_AB_p(x, y)
                  x_wx(B_interval%first - 1 + y, A_interval%first - 1 + x) = x_AB_p(x, y)
!
               enddo
            enddo
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine get_ao_x_wx_wavefunction
!
!
   subroutine get_ao_x_wx_1der_wavefunction(wf, construct_ao_x_AB_1der, x_wxqk)
!!
!!    Get AO X 1der 
!!    Written by Eirik F. Kjønstad, 2019 
!!
!!    Computes the geometrically differentiated AO integral matrix x_wxqk using  
!!    the construction routine construct_ao_x_AB_1der. This routine
!!    constructs the x_ABqk contributions to x_wxqk for the shells A and B.
!!
!!    The derivative indices for the nuclei are:
!! 
!!    q = 1, 2, 3 (x, y, z)
!!    k = 1, 2, 3, ..., n_atoms 
!!
!!    w and x are AO indices.
!!
      implicit none 
!
      class(wavefunction), intent(in) :: wf 
!
      procedure(construct_ao_x_wx_1der) :: construct_ao_x_AB_1der
!
      real(dp), dimension(wf%n_ao, wf%n_ao, 3, wf%system%n_atoms) :: x_wxqk 
!
      integer :: A, B, A_atom, B_atom, w, x, q, w_f, x_f 
!
      real(dp), dimension((wf%system%max_shell_size**2)*6), target :: x_ABqk 
!
      real(dp), dimension(:,:,:,:), pointer, contiguous :: x_ABqk_p 
!
      type(interval) :: A_interval, B_interval 
!     
      call zero_array(x_wxqk, 3*(wf%system%n_atoms)*(wf%n_ao)**2)
!
      do A = 1, wf%system%n_s
!
         A_atom     = wf%system%shell2atom(A)
         A_interval = wf%system%shell_limits(A)
!
         do B = 1, A
!
            B_atom     = wf%system%shell2atom(B)
            B_interval = wf%system%shell_limits(B)
!
            x_ABqk_p(1 : A_interval%length, 1 : B_interval%length, 1 : 3, 1 : 2) &
                                 => x_ABqk(1 : (A_interval%length)*(B_interval%length)*6)
!
            call construct_ao_x_AB_1der(wf%system,            &
                                       x_ABqk_p(:,:,1,1),   &
                                       x_ABqk_p(:,:,2,1),   &
                                       x_ABqk_p(:,:,3,1),   &
                                       x_ABqk_p(:,:,1,2),   &
                                       x_ABqk_p(:,:,2,2),   &
                                       x_ABqk_p(:,:,3,2),   &
                                       A, B)
!
            do q = 1, 3
               do w = 1, A_interval%length
                  do x = 1, B_interval%length
!
                     w_f = A_interval%first - 1 + w
                     x_f = B_interval%first - 1 + x
!
                     x_wxqk(w_f, x_f, q, A_atom) = x_wxqk(w_f, x_f, q, A_atom) + x_ABqk_p(w, x, q, 1)
                     x_wxqk(x_f, w_f, q, A_atom) = x_wxqk(x_f, w_f, q, A_atom) + x_ABqk_p(w, x, q, 1)
!
                     x_wxqk(w_f, x_f, q, B_atom) = x_wxqk(w_f, x_f, q, B_atom) + x_ABqk_p(w, x, q, 2)
                     x_wxqk(x_f, w_f, q, B_atom) = x_wxqk(x_f, w_f, q, B_atom) + x_ABqk_p(w, x, q, 2)
!
                  enddo
               enddo
            enddo
!
            nullify(x_ABqk_p)
!
         enddo
      enddo
!
   end subroutine get_ao_x_wx_1der_wavefunction
!
!
   subroutine initialize_orbital_coefficients_wavefunction(wf)
!!
!!    Initialize MO coefficients
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(wavefunction) :: wf
!
      if (.not. allocated(wf%orbital_coefficients)) call mem%alloc(wf%orbital_coefficients, wf%n_ao, wf%n_mo)
!
   end subroutine initialize_orbital_coefficients_wavefunction
!
!
   subroutine initialize_wavefunction_files_wavefunction(wf)
!!
!!    Initialize wavefunction files 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
      implicit none 
!
      class(wavefunction) :: wf 
!
      wf%orbital_coefficients_file = sequential_file('orbital_coefficients')
      wf%orbital_energies_file = sequential_file('orbital_energies')
!
   end subroutine initialize_wavefunction_files_wavefunction
!
!
   subroutine destruct_orbital_coefficients_wavefunction(wf)
!!
!!    Destruct MO coefficients
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(wavefunction) :: wf
!
      if (allocated(wf%orbital_coefficients)) call mem%dealloc(wf%orbital_coefficients, wf%n_ao, wf%n_mo)
!
   end subroutine destruct_orbital_coefficients_wavefunction
!
!
   subroutine initialize_orbital_energies_wavefunction(wf)
!!
!!    Initialize orbital energies
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(wavefunction) :: wf
!
      if (.not. allocated(wf%orbital_energies)) call mem%alloc(wf%orbital_energies, wf%n_mo)
!
   end subroutine initialize_orbital_energies_wavefunction
!
!
   subroutine destruct_orbital_energies_wavefunction(wf)
!!
!!    Destruct AO overlap
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(wavefunction) :: wf
!
      if (allocated(wf%orbital_energies)) call mem%dealloc(wf%orbital_energies, wf%n_mo)
!
   end subroutine destruct_orbital_energies_wavefunction
!
!
   subroutine mo_transform_wavefunction(wf, X_wx, Y_pq)
!!
!!    MO transform 
!!    Written by Eirik F. Kjønstad, Sep 2018 
!!
!!    Performs MO transformation of X and saves the result in Y:
!!
!!       Y_pq = sum_wx C_wp X_wx C_xq
!!
      implicit none 
!
      class(wavefunction), intent(in) :: wf 
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in)    :: X_wx 
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: Y_pq  
!
      real(dp), dimension(:,:), allocatable :: Z_wq ! = sum_x X_wx C_xq
!
      call mem%alloc(Z_wq, wf%n_ao, wf%n_mo)
!
      call dgemm('N', 'N',                 &
                  wf%n_ao,                 &
                  wf%n_mo,                 &
                  wf%n_ao,                 &
                  one,                     &
                  X_wx,                    &
                  wf%n_ao,                 &
                  wf%orbital_coefficients, & ! C_xq
                  wf%n_ao,                 &
                  zero,                    &
                  Z_wq,                    &
                  wf%n_ao)
!
      call dgemm('T', 'N',                 &
                  wf%n_mo,                 &
                  wf%n_mo,                 &
                  wf%n_ao,                 &
                  one,                     &
                  wf%orbital_coefficients, & ! C_wp 
                  wf%n_ao,                 &
                  Z_wq,                    &
                  wf%n_ao,                 &
                  zero,                    &
                  Y_pq,                    &
                  wf%n_mo)
!
      call mem%dealloc(Z_wq, wf%n_ao, wf%n_mo)
!
   end subroutine mo_transform_wavefunction
!
!
   subroutine get_ao_h_wx_wavefunction(wf, h)
!!
!!    Get AO h 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
!!    Constructs the one-electron Hamiltonian integrals h_pq 
!!    in the AO basis. 
!!
      use molecular_system_class, only: construct_ao_h_wx_molecular_system
!
      implicit none 
!
      class(wavefunction), intent(in) :: wf 
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(out) :: h 
!
      call wf%get_ao_x_wx(construct_ao_h_wx_molecular_system, h)
!
   end subroutine get_ao_h_wx_wavefunction
!
!
   subroutine get_ao_s_wx_wavefunction(wf, s)
!!
!!    Get AO s
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
!!    Constructs the one-electron overlap matrix s
!!    in the AO basis. 
!!
      use molecular_system_class, only: construct_ao_s_wx_molecular_system
!
      implicit none 
!
      class(wavefunction), intent(in) :: wf 
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(out) :: s 
!
      call wf%get_ao_x_wx(construct_ao_s_wx_molecular_system, s)
!
   end subroutine get_ao_s_wx_wavefunction
!
!
   subroutine get_ao_s_wx_1der_wavefunction(wf, s_wxqk)
!!
!!    Get AO s 1st derivative
!!    Written by Eirik F. Kjønstad, June 2019 
!!
!!    Constructs the derivative of the AO overlap matrix.
!!
      use molecular_system_class, only: construct_ao_s_wx_1der_molecular_system
!     
      implicit none 
!
      class(wavefunction), intent(in) :: wf 
!
      real(dp), dimension(wf%n_ao, wf%n_ao, 3, wf%system%n_atoms), intent(out) :: s_wxqk
! 
      call wf%get_ao_x_wx_1der(construct_ao_s_wx_1der_molecular_system, s_wxqk)
!
   end subroutine get_ao_s_wx_1der_wavefunction
!
!
   subroutine get_ao_h_wx_1der_wavefunction(wf, h_wxqk)
!!
!!    Get AO h 1st derivative
!!    Written by Eirik F. Kjønstad, June 2019 
!!
!!    Constructs the 1st derivative of the nuclear-electron attraction integrals.
!!
      use molecular_system_class, only: construct_ao_h_wx_kinetic_1der_molecular_system
!
      implicit none 
!
      class(wavefunction), intent(in) :: wf 
!
      real(dp), dimension(wf%n_ao, wf%n_ao, 3, wf%system%n_atoms), intent(out) :: h_wxqk
!
      real(dp), dimension(:,:,:,:), allocatable :: h_wxqk_nuclear 
!
      call wf%get_ao_x_wx_1der(construct_ao_h_wx_kinetic_1der_molecular_system, h_wxqk)
!
      call mem%alloc(h_wxqk_nuclear, wf%n_ao, wf%n_ao, 3, wf%system%n_atoms)
!
      call wf%system%construct_ao_h_wx_nuclear_1der(h_wxqk_nuclear, wf%n_ao)
!
      call daxpy(3*(wf%system%n_atoms)*wf%n_ao**2, one, h_wxqk_nuclear, 1, h_wxqk, 1)
!
      call mem%dealloc(h_wxqk_nuclear, wf%n_ao, wf%n_ao, 3, wf%system%n_atoms)
!
   end subroutine get_ao_h_wx_1der_wavefunction
!
!
   subroutine get_ao_mu_wx_wavefunction(wf, mu_X, mu_Y, mu_Z)
!!
!!    Get AO mu
!!    Written by Eirik F. Kjønstad, Sep 2018 
!!
!!    Constructs the full dipole integrals mu
!!    for the X, Y, and Z components.
!!
      implicit none 
!
      class(wavefunction), intent(in) :: wf 
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(out) :: mu_X
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(out) :: mu_Y
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(out) :: mu_Z
!
      type(interval) :: A_interval, B_interval
!
      integer :: x, y, A, B
!
      real(dp), dimension(:,:), pointer :: mu_AB_X_p
      real(dp), dimension(:,:), pointer :: mu_AB_Y_p
      real(dp), dimension(:,:), pointer :: mu_AB_Z_p
!
      real(dp), dimension(wf%system%max_shell_size**2), target :: mu_AB_X 
      real(dp), dimension(wf%system%max_shell_size**2), target :: mu_AB_Y 
      real(dp), dimension(wf%system%max_shell_size**2), target :: mu_AB_Z 
!
!$omp parallel do private(A, B, A_interval, B_interval, x, y, mu_AB_X, mu_AB_Y, mu_AB_Z, &
!$omp                       &   mu_AB_X_p, mu_AB_Y_p, mu_AB_Z_p)
      do A = 1, wf%system%n_s
!
         A_interval = wf%system%shell_limits(A)
!
         do B = 1, A
!
            B_interval = wf%system%shell_limits(B)
!
            call wf%system%construct_ao_mu_wx(mu_AB_X, mu_AB_Y, mu_AB_Z, A, B)
!
            mu_AB_X_p(1 : A_interval%length, 1 : B_interval%length) => mu_AB_X(1 : A_interval%length*B_interval%length)
            mu_AB_Y_p(1 : A_interval%length, 1 : B_interval%length) => mu_AB_Y(1 : A_interval%length*B_interval%length)
            mu_AB_Z_p(1 : A_interval%length, 1 : B_interval%length) => mu_AB_Z(1 : A_interval%length*B_interval%length)
!
            do x = 1, A_interval%length
               do y = 1, B_interval%length
!
                  mu_X(A_interval%first - 1 + x, B_interval%first - 1 + y) = mu_AB_X_p(x, y)
                  mu_X(B_interval%first - 1 + y, A_interval%first - 1 + x) = mu_AB_X_p(x, y)
!
                  mu_Y(A_interval%first - 1 + x, B_interval%first - 1 + y) = mu_AB_Y_p(x, y)
                  mu_Y(B_interval%first - 1 + y, A_interval%first - 1 + x) = mu_AB_Y_p(x, y)
!
                  mu_Z(A_interval%first - 1 + x, B_interval%first - 1 + y) = mu_AB_Z_p(x, y)
                  mu_Z(B_interval%first - 1 + y, A_interval%first - 1 + x) = mu_AB_Z_p(x, y)
!
               enddo
            enddo
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine get_ao_mu_wx_wavefunction
!
!
   subroutine get_ao_q_wx_wavefunction(wf, q_xx, q_xy, q_xz, q_yy, q_yz, q_zz)
!!
!!    Get AO q
!!    Written by Eirik F. Kjønstad, Sep 2018 
!!
!!    Uses the integral tool to construct the quadrupole integrals.
!!
      implicit none 
!
      class(wavefunction), intent(in) :: wf 
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(out) :: q_xx
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(out) :: q_xy
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(out) :: q_xz
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(out) :: q_yy
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(out) :: q_yz
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(out) :: q_zz
!
      type(interval) :: A_interval, B_interval
!
      integer :: x, y, A, B
!
      real(dp), dimension(:,:), pointer :: q_AB_xx_p
      real(dp), dimension(:,:), pointer :: q_AB_xy_p
      real(dp), dimension(:,:), pointer :: q_AB_xz_p
      real(dp), dimension(:,:), pointer :: q_AB_yy_p
      real(dp), dimension(:,:), pointer :: q_AB_yz_p
      real(dp), dimension(:,:), pointer :: q_AB_zz_p
!
      real(dp), dimension(wf%system%max_shell_size**2), target :: q_AB_xx 
      real(dp), dimension(wf%system%max_shell_size**2), target :: q_AB_xy 
      real(dp), dimension(wf%system%max_shell_size**2), target :: q_AB_xz 
      real(dp), dimension(wf%system%max_shell_size**2), target :: q_AB_yy 
      real(dp), dimension(wf%system%max_shell_size**2), target :: q_AB_yz 
      real(dp), dimension(wf%system%max_shell_size**2), target :: q_AB_zz 
!
!$omp parallel do private (A, B, x, y, A_interval, B_interval, q_AB_xx, q_AB_xy, q_AB_xz, q_AB_yy, &
!$omp   & q_AB_yz, q_AB_zz, q_AB_xx_p, q_AB_xy_p, q_AB_xz_p, q_AB_yy_p, q_AB_yz_p, q_AB_zz_p)
      do A = 1, wf%system%n_s
!
         A_interval = wf%system%shell_limits(A)
!
         do B = 1, A
!
            B_interval = wf%system%shell_limits(B)
!
            call wf%system%construct_ao_q_wx(q_AB_xx, q_AB_xy, q_AB_xz, q_AB_yy, q_AB_yz, q_AB_zz, A, B)
!
            q_AB_xx_p(1 : A_interval%length, 1 : B_interval%length) => q_AB_xx(1 : A_interval%length*B_interval%length)
            q_AB_xy_p(1 : A_interval%length, 1 : B_interval%length) => q_AB_xy(1 : A_interval%length*B_interval%length)
            q_AB_xz_p(1 : A_interval%length, 1 : B_interval%length) => q_AB_xz(1 : A_interval%length*B_interval%length)
            q_AB_yy_p(1 : A_interval%length, 1 : B_interval%length) => q_AB_yy(1 : A_interval%length*B_interval%length)
            q_AB_yz_p(1 : A_interval%length, 1 : B_interval%length) => q_AB_yz(1 : A_interval%length*B_interval%length)
            q_AB_zz_p(1 : A_interval%length, 1 : B_interval%length) => q_AB_zz(1 : A_interval%length*B_interval%length)
!
             do x = 1, A_interval%length
                do y = 1, B_interval%length
!
                   q_xx(A_interval%first - 1 + x, B_interval%first - 1 + y) = q_AB_xx_p(x, y)
                   q_xx(B_interval%first - 1 + y, A_interval%first - 1 + x) = q_AB_xx_p(x, y)
!
                   q_xy(A_interval%first - 1 + x, B_interval%first - 1 + y) = q_AB_xy_p(x, y)
                   q_xy(B_interval%first - 1 + y, A_interval%first - 1 + x) = q_AB_xy_p(x, y)
!
                   q_xz(A_interval%first - 1 + x, B_interval%first - 1 + y) = q_AB_xz_p(x, y)
                   q_xz(B_interval%first - 1 + y, A_interval%first - 1 + x) = q_AB_xz_p(x, y)
!
                   q_yy(A_interval%first - 1 + x, B_interval%first - 1 + y) = q_AB_yy_p(x, y)
                   q_yy(B_interval%first - 1 + y, A_interval%first - 1 + x) = q_AB_yy_p(x, y)
!
                   q_yz(A_interval%first - 1 + x, B_interval%first - 1 + y) = q_AB_yz_p(x, y)
                   q_yz(B_interval%first - 1 + y, A_interval%first - 1 + x) = q_AB_yz_p(x, y)
!
                   q_zz(A_interval%first - 1 + x, B_interval%first - 1 + y) = q_AB_zz_p(x, y)
                   q_zz(B_interval%first - 1 + y, A_interval%first - 1 + x) = q_AB_zz_p(x, y)
!
                enddo
             enddo
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine get_ao_q_wx_wavefunction
!
!
   subroutine save_orbital_coefficients_wavefunction(wf)
!!
!!    Save orbital coefficients 
!!    Written by Eirik F. Kjønstad, Oct 2018 
!!
      implicit none 
!
      class(wavefunction), intent(inout) :: wf 
!
      call wf%orbital_coefficients_file%open_('write', 'rewind')
!
      call wf%orbital_coefficients_file%write_(wf%orbital_coefficients, wf%n_ao*wf%n_mo)
!
      call wf%orbital_coefficients_file%close_
!
   end subroutine save_orbital_coefficients_wavefunction
!
!
   subroutine read_orbital_coefficients_wavefunction(wf)
!!
!!    Save orbital coefficients 
!!    Written by Eirik F. Kjønstad, Oct 2018 
!!
      implicit none 
!
      class(wavefunction), intent(inout) :: wf 
!
      call wf%orbital_coefficients_file%open_('read', 'rewind')
!
      call wf%orbital_coefficients_file%read_(wf%orbital_coefficients, wf%n_ao*wf%n_mo)
!
      call wf%orbital_coefficients_file%close_
!
   end subroutine read_orbital_coefficients_wavefunction
!
!
   subroutine save_orbital_energies_wavefunction(wf)
!!
!!    Save orbital energies 
!!    Written by Eirik F. Kjønstad, Oct 2018 
!!
      implicit none 
!
      class(wavefunction), intent(inout) :: wf 
!
      call wf%orbital_energies_file%open_('write', 'rewind')
!
      call wf%orbital_energies_file%write_(wf%orbital_energies, wf%n_mo)
!
      call wf%orbital_energies_file%close_
!
   end subroutine save_orbital_energies_wavefunction
!
!
   subroutine read_orbital_energies_wavefunction(wf)
!!
!!    Save orbital energies 
!!    Written by Eirik F. Kjønstad, Oct 2018 
!!
      implicit none 
!
      class(wavefunction), intent(inout) :: wf 
!
      call wf%orbital_energies_file%open_('read', 'rewind')
!
      call wf%orbital_energies_file%read_(wf%orbital_energies, wf%n_mo)
!
      call wf%orbital_energies_file%close_
!
   end subroutine read_orbital_energies_wavefunction
!
!
   subroutine is_restart_safe_wavefunction(wf, task)
!!
!!    Is restart safe?
!!    Written by Eirik F. Kjønstad, Mar 2019 
!!
      implicit none 
!
      class(wavefunction) :: wf 
!
      character(len=*), intent(in) :: task 
!
      call output%error_msg('Cannot restart for task ' // trim(task) // ' from abstract wavefunction ' // trim(wf%name_))
!
   end subroutine is_restart_safe_wavefunction
!
!
   subroutine get_mo_mu_wavefunction(wf, mu_pqk)
!!
!!    Get MO dipole operator
!!    Written by Sarai D. Folekstad, Apr 2019
!!
!!    Constructs 
!!
!!       mu_pqk, k = 1, 2, 3 (x, y, z)
!!
!!    in the MO basis.
!!
      implicit none
!      
      class(wavefunction), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo, 3), intent(out) :: mu_pqk
!
      real(dp), dimension(:,:,:), allocatable :: mu_wxk
!
      call mem%alloc(mu_wxk, wf%n_ao, wf%n_ao, 3)
!
      call wf%get_ao_mu_wx(mu_wxk(:,:,1), mu_wxk(:,:,2), mu_wxk(:,:,3))
!
      call wf%mo_transform(mu_wxk(:,:,1), mu_pqk(:,:,1))
      call wf%mo_transform(mu_wxk(:,:,2), mu_pqk(:,:,2))
      call wf%mo_transform(mu_wxk(:,:,3), mu_pqk(:,:,3))
!
      call mem%dealloc(mu_wxk, wf%n_ao, wf%n_ao, 3)
!
   end subroutine get_mo_mu_wavefunction
!
!
   subroutine get_mo_q_wavefunction(wf, q_pqk)
!!
!!    Get MO quadrupole operator
!!    Written by Eirik F. Kjønstad, 2019
!!
!!    Constructs 
!!
!!       q_pqk, k = 1, 2, 3, 4, 5, 6 (xx, xy, xz, yy, yz, and zz)
!!
!!    in the MO basis.
!!
      implicit none
!      
      class(wavefunction), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo, 6) :: q_pqk
!
      real(dp), dimension(:,:,:), allocatable :: q_wxk
!
      call mem%alloc(q_wxk, wf%n_ao, wf%n_ao, 6)
!
      call wf%get_ao_q_wx(q_wxk(:,:,1), q_wxk(:,:,2), q_wxk(:,:,3), q_wxk(:,:,4), q_wxk(:,:,5), q_wxk(:,:,6))
!
      call wf%mo_transform(q_wxk(:,:,1), q_pqk(:,:,1))
      call wf%mo_transform(q_wxk(:,:,2), q_pqk(:,:,2))
      call wf%mo_transform(q_wxk(:,:,3), q_pqk(:,:,3))
      call wf%mo_transform(q_wxk(:,:,4), q_pqk(:,:,4))
      call wf%mo_transform(q_wxk(:,:,5), q_pqk(:,:,5))
      call wf%mo_transform(q_wxk(:,:,6), q_pqk(:,:,6))
!
      call mem%dealloc(q_wxk, wf%n_ao, wf%n_ao, 6)
!
   end subroutine get_mo_q_wavefunction
!
!
   subroutine get_mo_h_wavefunction(wf, h)
!!
!!    Get MO h 
!!    Written by Sarai D. Folekstad, Apr 2019
!!
      implicit none
!      
      class(wavefunction), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo) :: h
      real(dp), dimension(:,:), allocatable :: h_wx
!
      call mem%alloc(h_wx, wf%n_ao, wf%n_ao)
!
      call wf%get_ao_h_wx(h_wx)
!
      call wf%mo_transform(h_wx, h)
!
      call mem%dealloc(h_wx, wf%n_ao, wf%n_ao)
!
   end subroutine get_mo_h_wavefunction
!
!
   subroutine projected_atomic_orbitals_wavefunction(wf, D, PAO_coeff, n_orbitals, first_ao)
!!
!!    Projected atomic orbitals
!!    Written by Linda Goletto and Sarai D. Folkestad, Jun 2019
!!
!!    Constructs the projected atomic orbitals (PAOs)
!!
!!       C^PAO = I - DS
!!
!!    NOTE: n_orbitals must be equal to n_ao unless
!!    a restricted orbital space is requested.
!! 
!!    For restricted space optional argument first_ao   
!!    may be used.
!!
      use array_utilities, only: zero_array
!
      implicit none
!
      class(wavefunction), intent(in) :: wf
!
      integer, intent(in) :: n_orbitals
!  
      real(dp), dimension(wf%n_ao,wf%n_ao), intent(in)      :: D
      real(dp), dimension(wf%n_ao,n_orbitals), intent(out)  :: PAO_coeff
!
      integer, intent(in), optional :: first_ao
! 
      integer :: local_first, local_last
!
      integer :: i
!
      real(dp), dimension(:,:), allocatable :: S
!
!     Set offsets
!
      local_first = 1
      if (present(first_ao)) local_first = first_ao  
!
      local_last = local_first + n_orbitals - 1
!
!     Sanity checks and safety guards
!
      if (n_orbitals .gt. wf%n_ao) call output%error_msg('number of orbitals for PAOs exceeds number of AOs')
      if ((local_first .lt. 1) .or. (local_first .gt. wf%n_ao)) call output%error_msg('First PAO exceeds number of AOs')
      if ((local_last .lt. 1) .or. (local_last .gt. wf%n_ao)) call output%error_msg('Last PAO exceeds number of AOs')
!
!     Construct PAO coefficients
!
      call zero_array(PAO_coeff, (wf%n_ao)*(n_orbitals))
!
!$omp parallel do private(i)
      do i = 1, n_orbitals
!
         PAO_coeff(i + local_first - 1,i) = one
!
      enddo
!$omp end parallel do
!  
      call mem%alloc(S, wf%n_ao, wf%n_ao)
      S = zero
!
      call wf%get_ao_s_wx(S)
!
      call dgemm('N', 'N',             &
                  wf%n_ao,             &
                  n_orbitals,          &
                  wf%n_ao,             &
                  -one,                &
                  D,                   &
                  wf%n_ao,             &
                  S(1,local_first),    &
                  wf%n_ao,             &
                  one,                 &
                  PAO_coeff,           &
                  wf%n_ao)
!
      call mem%dealloc(S, wf%n_ao, wf%n_ao)
!
   end subroutine projected_atomic_orbitals_wavefunction
!
!
   subroutine get_orbital_overlap_wavefunction(wf, orbital_coeff, n_orbitals, S)
!!
!!    Get orbital overlap
!!    Written by Sarai D. Folkestad and Linda Goletto, Jun 2019
!!
!!    Construct the orbital overlap
!!
!!       S = C^T S_AO C
!!
!!    for some set of orbital coefficients.
!!
      implicit none
!
      class(wavefunction) :: wf
!
      integer, intent(in) :: n_orbitals
!
      real(dp), dimension(wf%n_ao, n_orbitals), intent(in) :: orbital_coeff
!
      real(dp), dimension(n_orbitals, n_orbitals), intent(out) :: S
!
      real(dp), dimension(:,:), allocatable :: S_ao, X
!
      call mem%alloc(S_ao, wf%n_ao, wf%n_ao)
      call wf%get_ao_s_wx(S_ao)
!
      call mem%alloc(X, n_orbitals, wf%n_ao)
!
      call dgemm('T', 'N',       &
                  n_orbitals,    &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  one,           &
                  orbital_coeff, &
                  wf%n_ao,       &
                  S_ao,          &
                  wf%n_ao,       &
                  zero,          &
                  X,             &
                  n_orbitals)
!
      call mem%dealloc(S_ao, wf%n_ao, wf%n_ao)
!
      call dgemm('N', 'N',       &
                  n_orbitals,    &
                  n_orbitals,    &
                  wf%n_ao,       &
                  one,           &
                  X,             &
                  n_orbitals,    &
                  orbital_coeff, &
                  wf%n_ao,       &
                  zero,          &
                  S,             &
                  n_orbitals)
!
      call mem%dealloc(X, n_orbitals, wf%n_ao)
!
   end subroutine get_orbital_overlap_wavefunction
!
!
   subroutine lovdin_orthonormalization_wavefunction(wf, orbital_coeff, S, n_orbitals, rank)
!!
!!    Lövdin orthonormalization 
!!    Written by Linda Goletto and sarai D. Folkestad, Jun 2019
!!
!!    Orthonormalizes the orbital_coeff using Lövdin 
!!    orthonormalization
!!
!!    The orbital overlap
!!
!!       S = C^T S_AO C
!!
!!    is diagonalized
!!
!!       S = U λ U^T
!!
!!    and the orbitals are updated according to
!!
!!       C = C U λ^(-1/2).
!!
!!    Linear dependence is removed by screening on the eigenvalues
!!    with a threshold of 1.0 * 10^-6
!!
      use array_utilities, only: copy_and_scale
!
      implicit none
!
      class(wavefunction) :: wf
!
      integer, intent(in)  :: n_orbitals
      integer, intent(out) :: rank
!
      real(dp), dimension(wf%n_ao, n_orbitals), intent(in)     :: orbital_coeff
      real(dp), dimension(n_orbitals, n_orbitals), intent(inout)  :: S
!
      real(dp), dimension(:), allocatable :: eigenvalues, work, inv_sqrt_eig
      real(dp), dimension(:,:), allocatable :: orbital_coeff_copy
!
      integer :: i, info
!
      real(dp), parameter :: threshold = 1.0d-6
!
!     Diagonalize S
!
      call mem%alloc(eigenvalues, n_orbitals)
!
      call mem%alloc(work, 3*n_orbitals-1)
!
      call dscal(n_orbitals**2, -one, S, 1)
!
      call dsyev('V', 'L',       &
               n_orbitals,       &
               S,                &
               n_orbitals,       &
               eigenvalues,      &
               work,             &
               3*n_orbitals-1,   &
               info)
!
      call mem%dealloc(work, 3*n_orbitals-1)
!
!     Find the rank of S
!
      rank = 0
!
      do i = 1, n_orbitals
!
         if (abs(eigenvalues(i)) .lt. threshold) exit
!
         rank = rank + 1
!
      enddo
!
!     Invert and square-root the eigenvalues (screening done here)
!
      call mem%alloc(inv_sqrt_eig, rank)
!
      do i = 1, rank
!
         inv_sqrt_eig(i) = 1/(sqrt(abs(eigenvalues(i))))
!
      enddo
!
      call mem%dealloc(eigenvalues, n_orbitals)
!
!     S = U λ^(-1/2)
!
      do i = 1, rank
!
         S(:, i) = S(:, i)*inv_sqrt_eig(i)
!
      enddo
!
      call mem%dealloc(inv_sqrt_eig, rank)
!
      do i = rank + 1, n_orbitals
!
         S(:,i) = zero
!
      enddo
!
!     Transform orbital coefficients
!
      call mem%alloc(orbital_coeff_copy, wf%n_ao, n_orbitals)
!
      call copy_and_scale(one, orbital_coeff, orbital_coeff_copy, (n_orbitals)*(wf%n_ao))
!
      call dgemm('N', 'N',                   &
                  wf%n_ao,                   &
                  rank,                      &
                  n_orbitals,                &
                  one,                       &
                  orbital_coeff_copy,        &
                  wf%n_ao,                   &
                  S,                         &
                  n_orbitals,                &
                  zero,                      &
                  orbital_coeff,             &
                  wf%n_ao)
!
      call mem%dealloc(orbital_coeff_copy, wf%n_ao, n_orbitals)
!
   end subroutine lovdin_orthonormalization_wavefunction
!
!
   subroutine get_ao_v_wx_wavefunction(wf, V)
!!
!!    Get AO v
!!    Written by Tommaso Giovannini, May 2019
!!
!!    Constructs the elecrostatic potential integral matrix v.
!!
      use molecular_system_class
!
      implicit none 
!
      class(wavefunction), intent(in) :: wf 
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(out) :: V
!
      call wf%get_ao_x_wx(construct_ao_v_wx_molecular_system, V)
      call dscal(wf%n_ao**2, two, V, 1)
!
   end subroutine get_ao_v_wx_wavefunction
!
!
   subroutine construct_ao_electrostatics_wavefunction(wf, multipole, elec_nucl, what, elec_fock, &
                                                     & property_points, ao_density)
!!
!!    Construct electrostatic properties
!!    Written by Tommaso Giovannini, April 2019 
!!
!!-------------------------------------------------------------------------------
!!
!!    multipole       = 0          potential 
!!                      1          electric field
!!                                 
!!    elec_nucl       = 0          do only electronic contribution to property
!!                      1          do both electronic and nuclear contribution 
!!                    
!!    what            = 'fock'     fock contribution, i.e. q*V_αβ or mu*E_αβ
!!                      'prop'     requested property calculated at MM points V(D) or E(D)
!!
!!-------------------------------------------------------------------------------
!!    optional arguments
!!
!!    elec_fock       = matrix for fock contribution (dimension wf%n_ao, wf%n_ao)
!!                      MANDATORY if what.eq.'fock'
!!
!!    property_points = array for property contribution 
!!                      (dimension wf%system%mm%n_atoms if multipole = 0 or 
!!                       dimension 3*wf%system%mm%n_atoms if multipole = 1)
!!                      MANDATORY if what.eq.'prop'
!!
!!-------------------------------------------------------------------------------
!!
!!    Implementation through LibInt interface
!!
      implicit none 
!       
!     input variables
! 
      class(wavefunction) :: wf
      integer :: multipole
      integer :: elec_nucl
      character(len=4) :: what
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(inout), optional :: elec_fock
      real(dp), dimension(:), intent(inout), optional :: property_points
      real(dp), dimension(:,:), intent(inout), optional :: ao_density
! 
!     stuff for libint
! 
      integer(i6) :: n_variables_i6
      real(dp), dimension(1) :: fake_charge
! 
!     internal variables
! 
      integer  :: mm_atom
      integer  :: qm_atom
!
      real(dp) :: ddot
      real(dp) :: distQMMM
      real(dp), dimension(wf%n_ao,wf%n_ao) :: fock_internal_ao
!     
!     Sanity check for Fock creation 
!       
      if(elec_nucl.eq.1.and.what.eq.'fock') &
         call output%error_msg('Do both electronic and nuclear and what = fock')
! 
!     checking informations on elec_fock and property_points
! 
      if(what.eq.'fock') then
!      
         if(.not.present(elec_fock)) &
            call output%error_msg('Error in call electrostatic fock creation')
!      
      else if(what.eq.'prop') then
!      
         if(.not.present(property_points)) then 
!         
            call output%error_msg('Error in call electrostatic property points creation')
!            
         else 
!         
            if (.not.present(ao_density)) &
               call output%error_msg('Error in call electrostatic property points creation: no ao density')
!               
         endif
!         
         if(multipole.eq.0) then
!         
            if(size(property_points).ne.wf%system%mm%n_atoms) &
               call output%error_msg('Wrong dimension of property points in electrostatics integral')
!               
         else if(multipole.ge.1) then
!         
            call output%error_msg('Multipoles integral greater than 1 NYI')
!           
         endif
!         
      endif
!               
! 
!     Interface with LibInt for calculation of Fock
! 
      if(what.eq.'fock') then
! 
         if(multipole.eq.0) then 
!
            n_variables_i6 = int(size(wf%system%mm%charge),kind=i6)
!            
            call initialize_potential_c(wf%system%mm%charge, &
                                        wf%system%mm%coordinates*angstrom_to_bohr, &
                                        n_variables_i6)
            call wf%get_ao_v_wx(elec_fock)
!            
         else if(multipole.ge.1) then 
!            
            call output%error_msg('Multipole greater than 1 NYI')
!
         endif
! 
!     Calculation of properties @ points
! 
!     Electronic part: interface to LibInt
!     Nuclear part   : explicit calculation
!      
      else if(what.eq.'prop') then
!      
         if(multipole.eq.0) then
!      
            call zero_array(property_points,wf%system%mm%n_atoms)
!      
            do mm_atom = 1, wf%system%mm%n_atoms
!      
               call zero_array(fock_internal_ao,wf%n_ao*wf%n_ao)
!               
               n_variables_i6 = int(1,kind=i6)
               fake_charge(1) = 1.0d0
!               
               call initialize_potential_c(fake_charge, &
                                           wf%system%mm%coordinates(:,mm_atom)*angstrom_to_bohr, &
                                           n_variables_i6)
               call wf%get_ao_v_wx(fock_internal_ao)
!      
               property_points(mm_atom) = -one/two * ddot((wf%n_ao)**2,ao_density,1,fock_internal_ao,1)
!      
               if(elec_nucl.eq.1) then
!      
                  do qm_atom = 1, wf%system%n_atoms
!      
                     distQMMM =  dsqrt((wf%system%atoms(qm_atom)%x - wf%system%mm%coordinates(1,mm_atom))**2 + &
                                       (wf%system%atoms(qm_atom)%y - wf%system%mm%coordinates(2,mm_atom))**2 + &
                                       (wf%system%atoms(qm_atom)%z - wf%system%mm%coordinates(3,mm_atom))**2)
!      
                     distQMMM = angstrom_to_bohr*distQMMM
!      
                     property_points(mm_atom) = property_points(mm_atom) &
                                              - (wf%system%atoms(qm_atom)%number_)/distQMMM
!      
                  enddo
!      
               endif
!      
            enddo
!      
         else if (multipole.ge.1) then 
!         
            call output%error_msg('Multipole greater than 1 NYI')
!            
         endif
!      
      endif
!
!
   end subroutine construct_ao_electrostatics_wavefunction
!
!
   subroutine update_h_wx_mm_hf(wf, h_wx_eff)
!!
!!    Update one-electron Hamiltonian with fixed charges
!!    for non-polarizable QM/MM
!!    Written by Tommaso Giovannini, July 2019 for QM/MM
!!
      implicit none
!
      class(wavefunction) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(inout) :: h_wx_eff
!
      if(.not.allocated(wf%system%mm%nopol_h_wx)) then 
!      
         call mem%alloc(wf%system%mm%nopol_h_wx, wf%n_ao, wf%n_ao)
!         
         call zero_array(wf%system%mm%nopol_h_wx, wf%n_ao*wf%n_ao)
!
         call wf%construct_ao_electrostatics(0,0,'fock',elec_fock=wf%system%mm%nopol_h_wx)
!         
         call output%print_matrix('debug', 'Electrostatic Embedding h:', &
                                  wf%system%mm%nopol_h_wx, wf%n_ao, wf%n_ao)
      endif
!         
      h_wx_eff = h_wx_eff + half * wf%system%mm%nopol_h_wx
!         
!
      call output%print_matrix('debug', 'h_eff (QM + Electrostatic Embedding)', & 
                               h_wx_eff, wf%n_ao, wf%n_ao)
!
!
   end subroutine update_h_wx_mm_hf
!
!
   subroutine construct_and_write_mo_cholesky_wavefunction(wf, n_mo, mo_coeff, mo_cholesky_file)
!!
!!    Construct and write MO Cholesky 
!!    Written by Sarai D. Folkestad, Sep 2019
!!
!!    Reads the AO Cholesky vectors, transforms them
!!    by the MO coefficients and writes them to file.
!!
!!    
!!    n_mo:             Number of MOs in mo_coeff
!!
!!    mo_coeff:         MO coefficients used to transform to MO
!!                      basis. Dimension (wf%n_ao, n_mo)
!!
!!    mo_cholesky_file: File where the MO Cholesky vectors
!!                      are stored. This file must be initialized
!!                      before being passed to the routine
!!
!
      use direct_file_class, only: direct_file
      use batching_index_class, only: batching_index
      use reordering, only: sort_123_to_132
      use timings_class, only: timings
!
      implicit none
!
      class(wavefunction), intent(in) :: wf
!
      integer, intent(in) :: n_mo
!
      real(dp), dimension(wf%n_ao, n_mo), intent(in) :: mo_coeff
!
      type(direct_file), intent(inout) :: mo_cholesky_file
!
      real(dp), dimension(:,:,:), allocatable :: L_Jxy, L_Jxp, L_Jpx, L_Jpq
!
      integer :: p, q, x, y, rec_xy, rec_pq
!
      type(batching_index) :: batch_p, batch_y
!
      integer :: req0, req1_p, req1_y, req2, current_p_batch, current_y_batch
!
      type(timings) :: timer
!
      timer = timings('MO transform and write Cholesky vectors')
      call timer%turn_on()
!
      call wf%system%ao_cholesky_file%open_('read')
      call mo_cholesky_file%open_('write')
!
      req0 = 0
!
      req1_p =  max(2*wf%system%n_J*wf%n_ao, wf%system%n_J*wf%n_ao + wf%system%n_J*n_mo)
      req1_y = wf%system%n_J*wf%n_ao
!
      req2 = 0
!
!     Initialize batching variables
!
      batch_p = batching_index(n_mo)
      batch_y = batching_index(wf%n_ao)
!
      call mem%batch_setup(batch_p, batch_y, req0, req1_p, req1_y, req2)
!
      do current_p_batch = 1, batch_p%num_batches
!
         call batch_p%determine_limits(current_p_batch)
!
         call mem%alloc(L_Jxp, wf%system%n_J, wf%n_ao, batch_p%length)
         call zero_array(L_Jxp, (wf%system%n_J)*(wf%n_ao)*(batch_p%length))
!
         do current_y_batch = 1, batch_y%num_batches
!
            call batch_y%determine_limits(current_y_batch)
!
            call mem%alloc(L_Jxy, wf%system%n_J, wf%n_ao, batch_y%length)
!
            do x = 1, wf%n_ao
               do y = 1, batch_y%length
!
                  rec_xy = max(x,y + batch_y%first - 1)*(max(x,y + batch_y%first - 1)-3)/2 + x + y + batch_y%first - 1
!
                  call wf%system%ao_cholesky_file%read_(L_Jxy(:,x,y), rec_xy)
!
               enddo
            enddo
!
            call dgemm('N', 'N',                         &
                  wf%n_ao*wf%system%n_J,                 &
                  batch_p%length,                        &
                  batch_y%length,                        &
                  one,                                   &
                  L_Jxy,                                 &
                  wf%n_ao*wf%system%n_J,                 &
                  mo_coeff(batch_y%first,batch_p%first), &
                  wf%n_ao,                               &
                  one,                                   &
                  L_Jxp,                                 &
                  wf%n_ao*wf%system%n_J)
!
            call mem%dealloc(L_Jxy, wf%system%n_J, wf%n_ao, batch_y%length)
!
         enddo ! y batches
!
         call mem%alloc(L_Jpx, wf%system%n_J, batch_p%length, wf%n_ao)
         call sort_123_to_132(L_Jxp, L_Jpx, wf%system%n_J, wf%n_ao, batch_p%length)
         call mem%dealloc(L_Jxp, wf%system%n_J, wf%n_ao, batch_p%length)
!
         call mem%alloc(L_Jpq, wf%system%n_J, batch_p%length, n_mo)
!
         call dgemm('N', 'N',                         &
                     batch_p%length*wf%system%n_J,    &
                     n_mo,                            &
                     wf%n_ao,                         &
                     one,                             &
                     L_Jpx,                           &
                     batch_p%length*wf%system%n_J,    &
                     mo_coeff,                        &
                     wf%n_ao,                         &
                     zero,                            &
                     L_Jpq,                           &
                     batch_p%length*wf%system%n_J)
!
         call mem%dealloc(L_Jpx, wf%system%n_J, batch_p%length, wf%n_ao)
!
         do p = batch_p%first, batch_p%last
            do q = 1, p
!
               rec_pq = max(p,q)*(max(p,q)-3)/2 + p + q
!
               call mo_cholesky_file%write_(L_Jpq(:,p - batch_p%first + 1,q), rec_pq)
!
            enddo
         enddo

!
         call mem%dealloc(L_Jpq, wf%system%n_J, batch_p%length, n_mo)
!
      enddo ! p batches
!
      call wf%system%ao_cholesky_file%close_('keep')
      call mo_cholesky_file%close_('keep')
!
      call timer%turn_off()
!
   end subroutine construct_and_write_mo_cholesky_wavefunction
!
!
   subroutine construct_orbital_block_by_density_cd_wavefunction(wf, D, n_vectors, threshold, mo_offset, active_aos)
!!
!!    Construct orbital block  by Cholesky decomposition for density
!!    Written by Sarai D. Folkestad, Feb 2019
!!
!!    Cholesky decomposition of density D plus
!!    update of corresponding wavefunction MOs
!!
!!    See A. M. J. Sánchez de Merás, H. Koch, 
!!    I. G. Cuesta, and L. Boman (J. Chem. Phys. 132, 204105 (2010))
!!    for more information on active space generation
!!    using Cholesky decomposition
!!
!!    'D' : Density to be decomposed 
!!          Given in the AO basis
!!
!!    'n_vectors' : The number of Cholesky 
!!                  vectors (orbitals) constructed
!!
!!    'threshold' : The Cholesky decomposition threshold
!!
!!    'mo_offset' : Offset used to set the new MOs
!!
!!    'active_aos' : List of AOs on active atoms (optional)
!!                   This is used in case we want to partially
!!                   decompose density to get active orbitals
!!
!!
!
      use array_utilities, only: cholesky_decomposition_limited_diagonal, full_cholesky_decomposition_effective
!
      implicit none
!
      class(wavefunction), intent(inout) :: wf
!
      real(dp), dimension(wf%n_ao,wf%n_ao), intent(inout) :: D
!
      integer, intent(out) :: n_vectors
!
      real(dp), intent(in) :: threshold
!
      integer, intent(in) :: mo_offset
!
      integer, dimension(:), optional :: active_aos
!
      integer :: n_active_aos, ao, mo
!
      real(dp), dimension(:,:), allocatable :: cholesky_vec
      integer, dimension(:), allocatable  :: keep_vectors
!
      call mem%alloc(cholesky_vec, wf%n_ao, wf%n_ao)
!
      if (present(active_aos)) then
!
!        Active space generation by CD choosing pivots 
!        only on active atoms
!
         n_active_aos = size(active_aos)
!
         if (n_active_aos .gt. wf%n_ao) call output%error_msg('More active AOs than total AOs')
!      
         call cholesky_decomposition_limited_diagonal(D, cholesky_vec, wf%n_ao, &
                                                      n_vectors, threshold, &
                                                      n_active_aos, active_aos)
! 
!        Set the  MOs to be the ones to be frozen for CC
!
!$omp parallel do private(ao, mo)
         do mo = 1, n_vectors
            do ao = 1, wf%n_ao
!
               wf%orbital_coefficients(ao, mo_offset + mo) = cholesky_vec(ao, mo)
!
            enddo
         enddo
!$omp end parallel do
!
      else
!
         call mem%alloc(keep_vectors, wf%n_ao)
!
!        Full CD of density (used for inactive densities)
!
         call full_cholesky_decomposition_effective(D, cholesky_vec, &
                                             wf%n_ao, n_vectors, &
                                             threshold, keep_vectors)
!
!
         call mem%dealloc(keep_vectors, wf%n_ao)
!
!
!$omp parallel do private(ao, mo)
         do mo = 1, n_vectors
            do ao = 1, wf%n_ao
!
               wf%orbital_coefficients(ao, mo_offset + mo) = cholesky_vec(ao, mo)
!
            enddo
         enddo
!$omp end parallel do
!
      endif
!
      call mem%dealloc(cholesky_vec, wf%n_ao, wf%n_ao)
!
   end subroutine construct_orbital_block_by_density_cd_wavefunction
!
!
   subroutine read_frozen_orbitals_settings_wavefunction(wf)
!!
!!    Read frozen orbitals 
!!    Written by Sarai D. Folkestad, Oct 2019
!!
!!    Reads the frozen orbitals section of the input
!!
!!     - Frozen core 
!!
!!    - Frozen hf orbitals
!!
!!    This routine is used at HF level to prepare mos and 
!!    frozen fock contributions.
!!
!!    This routine is read at cc level to figure out if there
!!    should be frozen fock contributions
!!
      implicit none
!
      class(wavefunction) :: wf
!
      wf%frozen_core    = .false.
      wf%frozen_hf_mos  = .false.
!
      if (input%requested_keyword_in_section('core', 'frozen orbitals')) wf%frozen_core = .true.
      if (input%requested_keyword_in_section('hf', 'frozen orbitals')) wf%frozen_hf_mos = .true.
!
   end subroutine read_frozen_orbitals_settings_wavefunction
!
!
   subroutine initialize_mo_fock_frozen_hf_contribution_wavefunction(wf)
!!
!!    Initialize Fock frozen HF orbitals contributions
!!    Written by Sarai D. Folkestad, Sep. 2019
!!
      implicit none
!
      class(wavefunction) :: wf
!
      if (.not. allocated(wf%mo_fock_frozen_hf_contribution)) call mem%alloc(wf%mo_fock_frozen_hf_contribution, wf%n_mo, wf%n_mo)
!
   end subroutine initialize_mo_fock_frozen_hf_contribution_wavefunction
!
!
   subroutine destruct_mo_fock_frozen_hf_contribution_wavefunction(wf)
!!
!!    Destruct Fock frozen HF orbitals contributions
!!    Written by Sarai D. Folkestad, Sep. 2019
!!
      implicit none
!
      class(wavefunction) :: wf
!
      if (allocated(wf%mo_fock_frozen_hf_contribution)) call mem%dealloc(wf%mo_fock_frozen_hf_contribution, wf%n_mo, wf%n_mo)
!
   end subroutine destruct_mo_fock_frozen_hf_contribution_wavefunction
!
!
   subroutine initialize_mo_fock_fc_contribution_wavefunction(wf)
!!
!!    Initialize Fock frozen core
!!    Written by Sarai D. Folkestad, Sep. 2019
!!
      implicit none
!
      class(wavefunction) :: wf
!
      if (.not. allocated(wf%mo_fock_fc_contribution)) call mem%alloc(wf%mo_fock_fc_contribution, wf%n_mo, wf%n_mo)
!
   end subroutine initialize_mo_fock_fc_contribution_wavefunction
!
!
   subroutine destruct_mo_fock_fc_contribution_wavefunction(wf)
!!
!!    Destruct Fock frozen core
!!    Written by Sarai D. Folkestad, Sep. 2019
!!
      implicit none
!
      class(wavefunction) :: wf
!
      if (allocated(wf%mo_fock_fc_contribution)) call mem%dealloc(wf%mo_fock_fc_contribution, wf%n_mo, wf%n_mo)
!
   end subroutine destruct_mo_fock_fc_contribution_wavefunction
!
!
end module wavefunction_class
