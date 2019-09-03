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
   use array_utilities, only : zero_array, print_matrix
   use global_out, only : output
   use sequential_file_class, only : sequential_file
   use memory_manager_class, only : mem
   use molecular_system_class, only : molecular_system
   use interval_class, only : interval
   use libint_initialization
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
!     stuff for QM/MM
!  
!
      type(sequential_file) :: orbital_coefficients_file
      type(sequential_file) :: orbital_energies_file
!
   contains
!
      procedure :: initialize_orbital_coefficients => initialize_orbital_coefficients_wavefunction
      procedure :: initialize_orbital_energies     => initialize_orbital_energies_wavefunction
!
      procedure :: destruct_orbital_coefficients   => destruct_orbital_coefficients_wavefunction
      procedure :: destruct_orbital_energies       => destruct_orbital_energies_wavefunction
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
      procedure :: construct_ao_electrostatics              => construct_ao_electrostatics_wavefunction       ! V_αβ, E_αβ, V(D), E(D)
      procedure :: update_h_wx_mm                           => update_h_wx_mm_hf
!      
   end type wavefunction 
!
!
   abstract interface 
!
      subroutine gradient_function(wf, E_qk)
!
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
   end interface 
!
!
contains
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
!!    Constructs the full one-electron h matrix.
!!
      implicit none 
!
      class(wavefunction), intent(in) :: wf 
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(out) :: h 
!
      type(interval) :: A_interval, B_interval
!
      integer :: x, y, A, B 
!
      real(dp), dimension(:,:), pointer                        :: h_AB_p 
      real(dp), dimension(wf%system%max_shell_size**2), target :: h_AB
!
!$omp parallel do private(A, B, h_AB, h_AB_p, A_interval, B_interval, x, y) schedule(static)
      do A = 1, wf%system%n_s
!
         A_interval = wf%system%shell_limits(A)
!
         do B = 1, A
!
            B_interval = wf%system%shell_limits(B)
!
            call wf%system%construct_ao_h_wx(h_AB, A, B)
!
            h_AB_p(1 : A_interval%size, 1 : B_interval%size) => h_AB(1 : A_interval%size*B_interval%size)
!
            do x = 1, A_interval%size
               do y = 1, B_interval%size
!
                  h(A_interval%first - 1 + x, B_interval%first - 1 + y) = h_AB_p(x, y)
                  h(B_interval%first - 1 + y, A_interval%first - 1 + x) = h_AB_p(x, y)
!
               enddo
            enddo
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine get_ao_h_wx_wavefunction
!
!
   subroutine get_ao_s_wx_wavefunction(wf, s)
!!
!!    Get AO s
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
!!    Constructs the full one-electron overlap matrix s.
!!
      implicit none 
!
      class(wavefunction), intent(in) :: wf 
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(out) :: s 
!
      type(interval) :: A_interval, B_interval
!
      integer :: x, y, A, B
!
      real(dp), dimension(:,:), pointer                        :: s_AB_p 
      real(dp), dimension(wf%system%max_shell_size**2), target :: s_AB 
!
!$omp parallel do private(A, B, A_interval, B_interval, s_AB, s_AB_p, x, y)
      do A = 1, wf%system%n_s
!
         A_interval = wf%system%shell_limits(A)
!
         do B = 1, A
!
            B_interval = wf%system%shell_limits(B)
!
            call wf%system%construct_ao_s_wx(s_AB, A, B)
!
            s_AB_p(1 : A_interval%size, 1 : B_interval%size) => s_AB(1 : A_interval%size*B_interval%size)
!
            do x = 1, A_interval%size
               do y = 1, B_interval%size
!
                  s(A_interval%first - 1 + x, B_interval%first - 1 + y) = s_AB_p(x, y)
                  s(B_interval%first - 1 + y, A_interval%first - 1 + x) = s_AB_p(x, y)
!
               enddo
            enddo
!
            nullify(s_AB_p)
!
         enddo
      enddo
!$omp end parallel do
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
      implicit none 
!
      class(wavefunction), intent(in) :: wf 
!
      real(dp), dimension(wf%n_ao, wf%n_ao, 3, wf%system%n_atoms) :: s_wxqk
!
      integer :: A, B, A_atom, B_atom, w, x, q, w_f, x_f 
!
      real(dp), dimension((wf%system%max_shell_size**2)*6), target :: s_ABqk 
!
      real(dp), dimension(:,:,:,:), pointer, contiguous :: s_ABqk_p 
!
      type(interval) :: A_interval, B_interval 
!     
      s_wxqk = zero
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
            s_ABqk_p(1 : A_interval%size, 1 : B_interval%size, 1 : 3, 1 : 2) &
                                 => s_ABqk(1 : (A_interval%size)*(B_interval%size)*6)
!
            call wf%system%construct_ao_s_wx_1der(s_ABqk_p(:,:,1,1),    &
                                                   s_ABqk_p(:,:,2,1),   &
                                                   s_ABqk_p(:,:,3,1),   &
                                                   s_ABqk_p(:,:,1,2),   &
                                                   s_ABqk_p(:,:,2,2),   &
                                                   s_ABqk_p(:,:,3,2),   &
                                                   A, B)
!
            do q = 1, 3
               do w = 1, A_interval%size
                  do x = 1, B_interval%size
!
                     w_f = A_interval%first - 1 + w
                     x_f = B_interval%first - 1 + x
!
                     s_wxqk(w_f, x_f, q, A_atom) = s_wxqk(w_f, x_f, q, A_atom) + s_ABqk_p(w, x, q, 1)
                     s_wxqk(x_f, w_f, q, A_atom) = s_wxqk(x_f, w_f, q, A_atom) + s_ABqk_p(w, x, q, 1)
!
                     s_wxqk(w_f, x_f, q, B_atom) = s_wxqk(w_f, x_f, q, B_atom) + s_ABqk_p(w, x, q, 2)
                     s_wxqk(x_f, w_f, q, B_atom) = s_wxqk(x_f, w_f, q, B_atom) + s_ABqk_p(w, x, q, 2)
!
                  enddo
               enddo
            enddo
!
            nullify(s_ABqk_p)
!
         enddo
      enddo
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
      implicit none 
!
      class(wavefunction), intent(in) :: wf 
!
      real(dp), dimension(wf%n_ao, wf%n_ao, 3, wf%system%n_atoms), intent(out) :: h_wxqk
!
      integer :: A, B, A_atom, B_atom, w, x, q, w_f, x_f 
!
      real(dp), dimension((wf%system%max_shell_size**2)*3*2), target :: h_ABqk 
!
      real(dp), dimension(:,:,:,:), pointer, contiguous :: h_ABqk_p 
!
      type(interval) :: A_interval, B_interval 
!
      h_wxqk = zero
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
            h_ABqk_p(1 : A_interval%size, 1 : B_interval%size, 1 : 3, 1 : 2) &
                                 => h_ABqk(1 : A_interval%size*B_interval%size*3*2)
!
            call wf%system%construct_ao_h_wx_kinetic_1der(h_ABqk_p(:,:,1,1),     &
                                                            h_ABqk_p(:,:,2,1),   &
                                                            h_ABqk_p(:,:,3,1),   &
                                                            h_ABqk_p(:,:,1,2),   &
                                                            h_ABqk_p(:,:,2,2),   &
                                                            h_ABqk_p(:,:,3,2),   &
                                                            A, B)
!
            do q = 1, 3
               do w = 1, A_interval%size
                  do x = 1, B_interval%size
!
                     w_f = A_interval%first - 1 + w
                     x_f = B_interval%first - 1 + x
!
                     h_wxqk(w_f, x_f, q, A_atom) = h_wxqk(w_f, x_f, q, A_atom) + h_ABqk_p(w, x, q, 1)
                     h_wxqk(w_f, x_f, q, B_atom) = h_wxqk(w_f, x_f, q, B_atom) + h_ABqk_p(w, x, q, 2)
!
                  enddo
               enddo
            enddo
!
            if (A .ne. B) then 
!
               do q = 1, 3
                  do w = 1, A_interval%size
                     do x = 1, B_interval%size
   !
                        w_f = A_interval%first - 1 + w
                        x_f = B_interval%first - 1 + x
   !
                        h_wxqk(x_f, w_f, q, A_atom) = h_wxqk(x_f, w_f, q, A_atom) + h_ABqk_p(w, x, q, 1)
                        h_wxqk(x_f, w_f, q, B_atom) = h_wxqk(x_f, w_f, q, B_atom) + h_ABqk_p(w, x, q, 2)
   !
                     enddo
                  enddo
               enddo
!
            endif 
!
            nullify(h_ABqk_p)
!
         enddo
      enddo
!
      do A = 1, wf%system%n_s 
         do B = 1, A
!
            call wf%system%construct_and_add_ao_h_wx_nuclear_1der(h_wxqk, A, B, wf%n_ao)
!
         enddo
      enddo
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
            mu_AB_X_p(1 : A_interval%size, 1 : B_interval%size) => mu_AB_X(1 : A_interval%size*B_interval%size)
            mu_AB_Y_p(1 : A_interval%size, 1 : B_interval%size) => mu_AB_Y(1 : A_interval%size*B_interval%size)
            mu_AB_Z_p(1 : A_interval%size, 1 : B_interval%size) => mu_AB_Z(1 : A_interval%size*B_interval%size)
!
            do x = 1, A_interval%size
               do y = 1, B_interval%size
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
            q_AB_xx_p(1 : A_interval%size, 1 : B_interval%size) => q_AB_xx(1 : A_interval%size*B_interval%size)
            q_AB_xy_p(1 : A_interval%size, 1 : B_interval%size) => q_AB_xy(1 : A_interval%size*B_interval%size)
            q_AB_xz_p(1 : A_interval%size, 1 : B_interval%size) => q_AB_xz(1 : A_interval%size*B_interval%size)
            q_AB_yy_p(1 : A_interval%size, 1 : B_interval%size) => q_AB_yy(1 : A_interval%size*B_interval%size)
            q_AB_yz_p(1 : A_interval%size, 1 : B_interval%size) => q_AB_yz(1 : A_interval%size*B_interval%size)
            q_AB_zz_p(1 : A_interval%size, 1 : B_interval%size) => q_AB_zz(1 : A_interval%size*B_interval%size)
!
             do x = 1, A_interval%size
                do y = 1, B_interval%size
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
   subroutine get_mo_mu_wavefunction(wf, mu)
!!
!!    Get MO dipole operator
!!    Written by Sarai D. Folekstad, Apr 2019
!!
      implicit none
!      
      class(wavefunction), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo, 3) :: mu
      real(dp), dimension(:,:), allocatable :: mu_X_wx, mu_Y_wx, mu_Z_wx
!
      call mem%alloc(mu_X_wx, wf%n_ao, wf%n_ao)
      call mem%alloc(mu_Y_wx, wf%n_ao, wf%n_ao)
      call mem%alloc(mu_Z_wx, wf%n_ao, wf%n_ao)
!
      call wf%get_ao_mu_wx(mu_X_wx, mu_Y_wx, mu_Z_wx)
!
      call wf%mo_transform(mu_X_wx, mu(1,1,1))
      call wf%mo_transform(mu_Y_wx, mu(1,1,2))
      call wf%mo_transform(mu_Z_wx, mu(1,1,3))
!
      call mem%dealloc(mu_X_wx, wf%n_ao, wf%n_ao)
      call mem%dealloc(mu_Y_wx, wf%n_ao, wf%n_ao)
      call mem%dealloc(mu_Z_wx, wf%n_ao, wf%n_ao)
!
   end subroutine get_mo_mu_wavefunction
!
!
   subroutine get_mo_h_wavefunction(wf, h)
!!
!!    Get MO h (also for QM/MM calculations)
!!    Written by Sarai D. Folekstad, Apr 2019
!!    Modified by Tommaso Giovannini, July 2019 
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
      if(wf%system%mm_calculation.and.wf%system%mm%forcefield.eq.'non-polarizable') &
         h_wx = h_wx + half*wf%system%mm%nopol_h_wx
!         
      if(wf%system%mm_calculation.and.wf%system%mm%forcefield.ne.'non-polarizable') &
         h_wx = h_wx + half*wf%system%mm%pol_emb_fock
!
      call wf%mo_transform(h_wx, h)
!
      call mem%dealloc(h_wx, wf%n_ao, wf%n_ao)
!
   end subroutine get_mo_h_wavefunction
!
!
   subroutine get_ao_v_wx_wavefunction(wf, V)
!!
!!    Get AO h 
!!    Written by Tommaso Giovannini, May 2019
!!
!!    Uses the integral tool to construct the elecrostatic potential
!!
      implicit none 
!
      class(wavefunction), intent(in) :: wf 
!
      real(dp), dimension(wf%n_ao, wf%n_ao) :: V
!
      type(interval) :: A_interval, B_interval
!
      integer :: x, y, A, B
!
      real(dp), dimension(:,:), allocatable :: V_AB 
!
!$omp parallel do &
!$omp private(A, B, V_AB, A_interval, B_interval, x, y) schedule(static)
      do A = 1, wf%system%n_s
!
         A_interval = wf%system%shell_limits(A)
!
         do B = 1, A
!
            B_interval = wf%system%shell_limits(B)
!
            call mem%alloc(V_AB, A_interval%size, B_interval%size)
            call wf%system%construct_ao_v_wx(V_AB, A, B)
!
             do x = 1, A_interval%size
                do y = 1, B_interval%size
!
                   V(A_interval%first - 1 + x, B_interval%first - 1 + y) = two * V_AB(x, y)
                   V(B_interval%first - 1 + y, A_interval%first - 1 + x) = two * V_AB(x, y)
!
                enddo
             enddo
!
            call mem%dealloc(V_AB, A_interval%size, B_interval%size)
!
         enddo
      enddo
!$omp end parallel do
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
         if(wf%system%mm%verbose.ge.3) then 
!          
            call print_matrix('Electrostatic Embedding h:',wf%system%mm%nopol_h_wx,wf%n_ao,wf%n_ao)
!         
            flush(output%unit)
!
         endif
!         
      endif
!         
      h_wx_eff = h_wx_eff + half * wf%system%mm%nopol_h_wx
!         
      if(wf%system%mm%verbose.ge.3) then 
!      
        call print_matrix('h_eff (QM + Electrostatic Embedding)',h_wx_eff,wf%n_ao,wf%n_ao) 
!         
        flush(output%unit)
!         
      endif
!
!
   end subroutine update_h_wx_mm_hf
!
!
end module wavefunction_class