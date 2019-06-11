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
   use kinds
   use file_class
   use disk_manager_class
   use molecular_system_class
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
      type(file) :: orbital_coefficients_file
      type(file) :: orbital_energies_file
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
      procedure :: get_ao_h_wx_1der_numerical      => get_ao_h_wx_1der_numerical_wavefunction
!
      procedure :: get_ao_s_wx                     => get_ao_s_wx_wavefunction
      procedure :: get_ao_s_wx_1der                => get_ao_s_wx_1der_wavefunction
      procedure :: get_ao_s_wx_1der_numerical      => get_ao_s_wx_1der_numerical_wavefunction
!
      procedure :: get_ao_mu_wx                    => get_ao_mu_wx_wavefunction
      procedure :: get_ao_q_wx                     => get_ao_q_wx_wavefunction
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
      call wf%orbital_coefficients_file%init('orbital_coefficients', 'sequential', 'unformatted')
      call wf%orbital_energies_file%init('orbital_energies', 'sequential', 'unformatted')
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
      implicit none 
!
      class(wavefunction), intent(in) :: wf 
!
      real(dp), dimension(wf%n_ao, wf%n_ao, 3, wf%system%n_atoms) :: s_wxqk
!
      integer :: A, B, A_atom, B_atom, w, x, q
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
         A_atom     = wf%system%shell_to_atom(A)
         A_interval = wf%system%shell_limits(A)
!
         do B = 1, A
!
            B_atom     = wf%system%shell_to_atom(B)
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
                     s_wxqk(A_interval%first - 1 + w, B_interval%first - 1 + x, q, A_atom) = s_ABqk_p(w, x, q, 1)
                     s_wxqk(B_interval%first - 1 + x, A_interval%first - 1 + w, q, A_atom) = s_ABqk_p(w, x, q, 1)
!
                     s_wxqk(A_interval%first - 1 + w, B_interval%first - 1 + x, q, B_atom) = s_ABqk_p(w, x, q, 2)
                     s_wxqk(B_interval%first - 1 + x, A_interval%first - 1 + w, q, B_atom) = s_ABqk_p(w, x, q, 2)
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
   subroutine get_ao_s_wx_1der_numerical_wavefunction(wf, s_wxqk, dx)
!!
!!    Get AO s 1st derivative numerically
!!    Written by Eirik F. Kjønstad, June 2019 
!!
!!    Uses forward differences to calculate the derivative.
!!
      implicit none 
!
      class(wavefunction), intent(in) :: wf 
!
      real(dp), dimension(wf%n_ao, wf%n_ao, 3, wf%system%n_atoms) :: s_wxqk
!
      real(dp), intent(in) :: dx ! Displacement length in bohr
!
      real(dp), dimension(:,:), allocatable :: s_wx, s_wx_displaced, R_qk, R_qk_displaced
!
      integer :: k, q, w, x 
!
!     Get s at the reference geometry 
!
      call mem%alloc(s_wx, wf%n_ao, wf%n_ao)
      call wf%get_ao_s_wx(s_wx)
!
      call mem%alloc(R_qk, 3, wf%system%n_atoms)
      call mem%alloc(R_qk_displaced, 3, wf%system%n_atoms)
!
      R_qk = wf%system%get_geometry()
      call mem%alloc(s_wx_displaced, wf%n_ao, wf%n_ao)
!
      do k = 1, wf%system%n_atoms
         do q = 1, 3
!
!           Get s at the displaced geometry 
!
            R_qk_displaced = R_qk 
            R_qk_displaced(q,k) = R_qk_displaced(q,k) + dx 
            call wf%system%set_geometry(R_qk_displaced)
!
            call wf%get_ao_s_wx(s_wx_displaced)
!
!           Use difference from reference to compute derivative
!
            do w = 1, wf%n_ao
               do x = 1, wf%n_ao 
!
                  s_wxqk(w,x,q,k) = (s_wx_displaced(w,x) - s_wx(w,x))/dx
!
               enddo
            enddo
!
         enddo
      enddo
!    
      call wf%system%set_geometry(R_qk)
!
      call mem%dealloc(s_wx, wf%n_ao, wf%n_ao)
      call mem%dealloc(R_qk, 3, wf%system%n_atoms)
      call mem%dealloc(R_qk_displaced, 3, wf%system%n_atoms)
      call mem%dealloc(s_wx_displaced, wf%n_ao, wf%n_ao)
!
   end subroutine get_ao_s_wx_1der_numerical_wavefunction
!
!
   subroutine get_ao_h_wx_1der_wavefunction(wf, h_wxqk)
!!
!!    Get AO h 1st derivative
!!    Written by Eirik F. Kjønstad, June 2019 
!!
      implicit none 
!
      class(wavefunction), intent(in) :: wf 
!
      real(dp), dimension(wf%n_ao, wf%n_ao, 3, wf%system%n_atoms), intent(inout) :: h_wxqk
!
      integer :: A, B, A_atom, B_atom, w, x, q, k
!
      real(dp), dimension((wf%system%max_shell_size**2)*3*(wf%system%n_atoms)), target :: h_ABqk 
!
      real(dp), dimension(:,:,:,:), pointer, contiguous :: h_ABqk_p 
!
      type(interval) :: A_interval, B_interval 
!
      h_wxqk = zero
!
      do A = 1, wf%system%n_s
!
         A_atom     = wf%system%shell_to_atom(A)
         A_interval = wf%system%shell_limits(A)
!
         do B = 1, A
!
            B_atom     = wf%system%shell_to_atom(B)
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
                  h_wxqk(A_interval%first - 1 + w, B_interval%first - 1 + x, q, A_atom) = h_ABqk_p(w, x, q, 1)
                  h_wxqk(B_interval%first - 1 + x, A_interval%first - 1 + w, q, A_atom) = h_ABqk_p(w, x, q, 1)
                  h_wxqk(A_interval%first - 1 + w, B_interval%first - 1 + x, q, B_atom) = h_ABqk_p(w, x, q, 2)
                  h_wxqk(B_interval%first - 1 + x, A_interval%first - 1 + w, q, B_atom) = h_ABqk_p(w, x, q, 2)
!
                  enddo
               enddo
            enddo
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
   subroutine get_ao_h_wx_1der_numerical_wavefunction(wf, s_wxqk, dx)
!!
!!    Get AO h 1st derivative numerically
!!    Written by Eirik F. Kjønstad, June 2019 
!!
!!    Uses forward differences to calculate the derivative.
!!
      implicit none 
!
      class(wavefunction), intent(in) :: wf 
!
      real(dp), dimension(wf%n_ao, wf%n_ao, 3, wf%system%n_atoms) :: s_wxqk
!
      real(dp), intent(in) :: dx ! Displacement length in bohr
!
      real(dp), dimension(:,:), allocatable :: s_wx, s_wx_displaced, R_qk, R_qk_displaced
!
      integer :: k, q, w, x 
!
!     Get s at the reference geometry 
!
      call mem%alloc(s_wx, wf%n_ao, wf%n_ao)
      call wf%get_ao_h_wx(s_wx)
!
      call mem%alloc(R_qk, 3, wf%system%n_atoms)
      call mem%alloc(R_qk_displaced, 3, wf%system%n_atoms)
!
      R_qk = wf%system%get_geometry()
      call mem%alloc(s_wx_displaced, wf%n_ao, wf%n_ao)
!
      do k = 1, wf%system%n_atoms
         do q = 1, 3
!
!           Get s at the displaced geometry 
!
            R_qk_displaced = R_qk 
            R_qk_displaced(q,k) = R_qk_displaced(q,k) + dx 
            call wf%system%set_geometry(R_qk_displaced)
!
            call wf%get_ao_h_wx(s_wx_displaced)
!
!           Use difference from reference to compute derivative
!
            do w = 1, wf%n_ao
               do x = 1, wf%n_ao 
!
                  s_wxqk(w,x,q,k) = (s_wx_displaced(w,x) - s_wx(w,x))/dx
!
               enddo
            enddo
!
         enddo
      enddo
!    
      call wf%system%set_geometry(R_qk)
!
      call mem%dealloc(s_wx, wf%n_ao, wf%n_ao)
      call mem%dealloc(R_qk, 3, wf%system%n_atoms)
      call mem%dealloc(R_qk_displaced, 3, wf%system%n_atoms)
      call mem%dealloc(s_wx_displaced, wf%n_ao, wf%n_ao)
!
   end subroutine get_ao_h_wx_1der_numerical_wavefunction
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
      call disk%open_file(wf%orbital_coefficients_file, 'write', 'rewind')
!
      write(wf%orbital_coefficients_file%unit) wf%orbital_coefficients
!
      call disk%close_file(wf%orbital_coefficients_file)
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
      call wf%is_restart_safe('ground state')
!
      call disk%open_file(wf%orbital_coefficients_file, 'read', 'rewind')
!
      read(wf%orbital_coefficients_file%unit) wf%orbital_coefficients
!
      call disk%close_file(wf%orbital_coefficients_file)
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
      call disk%open_file(wf%orbital_energies_file, 'write', 'rewind')
!
      write(wf%orbital_energies_file%unit) wf%orbital_energies
!
      call disk%close_file(wf%orbital_energies_file)
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
      call wf%is_restart_safe('ground state')
!
      call disk%open_file(wf%orbital_energies_file, 'read', 'rewind')
!
      read(wf%orbital_energies_file%unit) wf%orbital_energies
!
      call disk%close_file(wf%orbital_energies_file)
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
end module wavefunction_class
