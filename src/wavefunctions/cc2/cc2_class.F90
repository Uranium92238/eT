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
module cc2_class
!
!!
!!    Coupled cluster singles and perturbative doubles (CC2) class module
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 2018
!!
!
   use ccs_class
!
   implicit none
!
   type, extends(ccs) :: cc2
!
      integer :: n_t2
!
      real(dp), dimension(:,:,:,:), allocatable :: u
!
      type(file) :: r2_file, l2_file
!
   contains
!
      procedure :: construct_u                                 => construct_u_cc2
!
      procedure :: construct_omega                             => construct_omega_cc2
!
      procedure :: omega_cc2_a1                                => omega_cc2_a1_cc2
      procedure :: omega_cc2_b1                                => omega_cc2_b1_cc2
      procedure :: omega_cc2_c1                                => omega_cc2_c1_cc2
!
      procedure :: calculate_energy                            => calculate_energy_cc2
!
      procedure :: jacobian_transform_trial_vector             => jacobian_transform_trial_vector_cc2
      procedure :: jacobian_cc2_transformation                 => jacobian_cc2_transformation_cc2
!
      procedure :: jacobian_cc2_a1                             => jacobian_cc2_a1_cc2
      procedure :: jacobian_cc2_b1                             => jacobian_cc2_b1_cc2
      procedure :: jacobian_cc2_a2                             => jacobian_cc2_a2_cc2
      procedure :: jacobian_cc2_b2                             => jacobian_cc2_b2_cc2
!
      procedure :: prepare_for_jacobian_transpose              => prepare_for_jacobian_transpose_cc2
!
      procedure :: jacobian_transpose_transform_trial_vector   => jacobian_transpose_transform_trial_vector_cc2
      procedure :: jacobian_transpose_cc2_transformation       => jacobian_transpose_cc2_transformation_cc2
!
      procedure :: jacobian_transpose_cc2_a1                   => jacobian_transpose_cc2_a1_cc2
      procedure :: jacobian_transpose_cc2_b1                   => jacobian_transpose_cc2_b1_cc2
      procedure :: jacobian_transpose_cc2_a2                   => jacobian_transpose_cc2_a2_cc2
      procedure :: jacobian_transpose_cc2_b2                   => jacobian_transpose_cc2_b2_cc2
!
      procedure :: initialize_u                                => initialize_u_cc2 
      procedure :: destruct_u                                  => destruct_u_cc2 
!
      procedure :: initialize_amplitudes                       => initialize_amplitudes_cc2 
      procedure :: destruct_amplitudes                         => destruct_amplitudes_cc2 
!
      procedure :: get_es_orbital_differences                  => get_es_orbital_differences_cc2
!
      procedure :: construct_multiplier_equation               => construct_multiplier_equation_cc2
!
      procedure :: get_cvs_projector                           => get_cvs_projector_cc2
!
      procedure :: initialize_files                            => initialize_files_cc2
      procedure :: initialize_doubles_files                    => initialize_doubles_files_cc2
!
      procedure :: save_doubles_vector                         => save_doubles_vector_cc2
      procedure :: read_doubles_vector                         => read_doubles_vector_cc2
!
      procedure :: restart_excited_state                       => restart_excited_state_cc2
      procedure :: save_excited_state                          => save_excited_state_cc2
      procedure :: is_restart_safe                             => is_restart_safe_cc2
!
   end type cc2
!
   interface
!
      include "omega_cc2_interface.F90"
      include "jacobian_cc2_interface.F90"
      include "jacobian_transpose_cc2_interface.F90"
!
   end interface 
!
!
   interface cc2
!
      procedure :: new_cc2 
!
   end interface cc2
!
!
contains
!
!
   function new_cc2(system) result(wf)
!!
!!    New CC2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      type(cc2) :: wf
!
      class(molecular_system), target, intent(in) :: system 
!
      type(file) :: hf_restart_file 
!
      wf%name_ = 'cc2'
!
      wf%system => system
!
      call hf_restart_file%init('hf_restart_file', 'sequential', 'unformatted')
!
      call disk%open_file(hf_restart_file, 'read', 'rewind')
!
      read(hf_restart_file%unit) wf%n_ao 
      read(hf_restart_file%unit) wf%n_mo 
      read(hf_restart_file%unit) 
      read(hf_restart_file%unit) wf%n_o  
      read(hf_restart_file%unit) wf%n_v  
      read(hf_restart_file%unit) wf%hf_energy  
!
      call disk%close_file(hf_restart_file)
!
      wf%n_t1            = (wf%n_o)*(wf%n_v)
      wf%n_t2            = wf%n_t1*(wf%n_t1+1)/2
      wf%n_gs_amplitudes = wf%n_t1
      wf%n_es_amplitudes = wf%n_t1 + wf%n_t2
!
      call wf%initialize_files()
!
      call wf%initialize_orbital_coefficients()
      call wf%initialize_orbital_energies()
!
      call wf%read_orbital_coefficients()
      call wf%read_orbital_energies()
!
      call wf%initialize_fock_ij()
      call wf%initialize_fock_ia()
      call wf%initialize_fock_ai()
      call wf%initialize_fock_ab()
!
   end function new_cc2
!
!
   subroutine calculate_energy_cc2(wf)
!!
!!    Calculate energy 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Jan 2019
!!
!!    E = E_HF + sum_aibj (t_i^a*t_j^b + t_ij^ab) L_iajb
!!
      class(cc2), intent(inout) :: wf 
!
      real(dp), dimension(:,:,:,:), allocatable :: g_aibj, g_iajb 
!
      real(dp) :: correlation_energy
!
      integer :: a, i, b, j
!
      call mem%alloc(g_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_vovo(g_aibj)
      call wf%get_ovov(g_iajb)
!
      correlation_energy = zero 
!
!$omp parallel do private(a,i,b,j) reduction(+:correlation_energy)
      do b = 1, wf%n_v
         do i = 1, wf%n_o 
            do j = 1, wf%n_o 
               do a = 1, wf%n_v
!
                  correlation_energy = correlation_energy +                                &
                                       (wf%t1(a, i)*wf%t1(b, j) -                          &
                                       (g_aibj(a,i,b,j))/(wf%orbital_energies(wf%n_o + a)  &
                                                         + wf%orbital_energies(wf%n_o + b) &
                                                         - wf%orbital_energies(i)           &
                                                         - wf%orbital_energies(j)))         &
                                       *(two*g_iajb(i,a,j,b)-g_iajb(i,b,j,a))
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do 
!
      call mem%dealloc(g_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      wf%energy = wf%hf_energy + correlation_energy
!
   end subroutine calculate_energy_cc2
!
!
   subroutine construct_u_cc2(wf)
!!
!!    Construct U 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Jan 2019
!!
!!    Construct
!!
!!       u_aibj = 2t_aibj - t_ajbi
!!
!!    with
!!
!!       t_aibj = - g_aibj/ε_aibj
!!
!!    where
!!
!!       ε_aibj = ε_a - ε_i + ε_b - ε_j 
!!
!!    and ε_r is the r'th orbital energy.
!!
      implicit none
!
      class(cc2), intent(inout) :: wf
!
      real(dp), dimension(:,:,:,:), allocatable :: g_aibj
!
      integer :: a, i, b, j
!
      call mem%alloc(g_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)  
!
      call wf%get_vovo(g_aibj)
!
      do b = 1, wf%n_v 
         do j = 1, wf%n_o 
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  wf%u(a, i, b, j) = (two*g_aibj(a, i, b, j) - g_aibj(a, j, b, i))/ &
                                       (  wf%orbital_energies(i) + wf%orbital_energies(j) &
                                        - wf%orbital_energies(wf%n_o + a) &
                                        - wf%orbital_energies(wf%n_o + b))
!
               enddo
            enddo
         enddo 
      enddo
!    
      call mem%dealloc(g_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)      
!
   end subroutine construct_u_cc2
!
!
   subroutine initialize_u_cc2(wf)
!!
!!    Initialize u 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Jan 2019
!!
      implicit none
!
      class(cc2) :: wf
!
      if (.not. allocated(wf%u)) call mem%alloc(wf%u, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine initialize_u_cc2
!
!
   subroutine destruct_u_cc2(wf)
!!
!!    Initialize u 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Jan 2019
!!
      implicit none
!
      class(cc2) :: wf
!
      if (allocated(wf%u)) call mem%dealloc(wf%u, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine destruct_u_cc2
!
!
   subroutine initialize_amplitudes_cc2(wf)
!!
!!    Initialize amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Jan 2019
!!
      implicit none
!
      class(cc2) :: wf
!
      call wf%initialize_t1()
      call wf%initialize_u()
!
   end subroutine initialize_amplitudes_cc2
!
!
   subroutine destruct_amplitudes_cc2(wf)
!!
!!    Destruct amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Jan 2019
!!
      implicit none
!
      class(cc2) :: wf
!
      call wf%destruct_t1()
      call wf%destruct_u()
!
   end subroutine destruct_amplitudes_cc2
!
!
   subroutine get_es_orbital_differences_cc2(wf, orbital_differences, N)
!!
!!    Get orbital differences 
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!    and Andreas Skeidsvoll, 2018
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      integer, intent(in) :: N 
      real(dp), dimension(N), intent(inout) :: orbital_differences
!
      integer :: a, i, ai, b, j, bj, aibj
!
!$omp parallel do schedule(static) private(a, i, b, j, ai, bj, aibj) 
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            ai = wf%n_v*(i - 1) + a
!
            orbital_differences(ai) = wf%orbital_energies(a + wf%n_o) - wf%orbital_energies(i)
!
            do j = 1, wf%n_o 
               do b = 1, wf%n_v
!
                  bj = wf%n_v*(j-1) + b 
!
                  if (ai .ge. bj) then
!
                     aibj = (ai*(ai-3)/2) + ai + bj
!
                     orbital_differences(aibj + (wf%n_o)*(wf%n_v)) = wf%orbital_energies(a + wf%n_o)     &
                                                                     - wf%orbital_energies(i) &
                                                                     +  wf%orbital_energies(b + wf%n_o)  &
                                                                     - wf%orbital_energies(j)
!
                  endif
!
               enddo
            enddo  
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine get_es_orbital_differences_cc2
!
!
   subroutine construct_multiplier_equation_cc2(wf, equation)
!!
!!    Construct multiplier equation 
!!    Written by Sarai D. Folkestad, Feb 2019
!!
!!    Constructs 
!!
!!       t-bar^T A + eta,
!!
!!    and places the result in 'equation'.
!!
!!    Solves analytically for tbar_aibj
!!
!!       tbar_aibj = - (η_aibj + sum_ai tbar_ai A_ai,aibj)/ε_aibj
!!
!!    where
!!
!!       η_aibj = 2 L_iajb       
!!
!!    and uses this to set up 'equation'
!!
!!       η_ai + sum_bj tbar_bj A_bj,ai + sum_bjck tbar_bjck A_{bjck,ai}
!!
      implicit none 
!
      class(cc2), intent(in) :: wf 
!
      real(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: equation 
!
      real(dp), dimension(:), allocatable :: eta 
      real(dp), dimension(:,:,:,:), allocatable :: t2bar
      real(dp), dimension(:,:,:,:), allocatable :: g_iajb
!
      integer :: a, b, i, j
!
!     Construct t2bar
!
      call mem%alloc(t2bar, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      t2bar = zero
!
!     t2bar = sum_ai tbar_ai A_ai,aibj
!
      call wf%jacobian_transpose_cc2_a2(t2bar, wf%t1bar)
!
      call mem%alloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call wf%get_ovov(g_iajb)
!
!     t2bar += η_aibj
!
      call add_2143_to_1234(four, g_iajb, t2bar, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call add_2341_to_1234(-two, g_iajb, t2bar, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     t2bar = t2bar/(-ε_aibj)
!
!$omp parallel do private(a, b, i, j)
      do b = 1, wf%n_v
         do j = 1, wf%n_o
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  t2bar(a, i, b, j) = t2bar(a, i, b, j)/(- wf%orbital_energies(a + wf%n_o) &
                                                         -  wf%orbital_energies(b + wf%n_o) &
                                                         +  wf%orbital_energies(i) &
                                                         +  wf%orbital_energies(j))
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
!     Set up the multipliers equation
!
!     equation = sum_bj tbar_bj A_bj,ai
!
      equation = zero
!  
      call wf%jacobian_transpose_ccs_a1(equation, wf%t1bar)
      call wf%jacobian_transpose_ccs_b1(equation, wf%t1bar)
      call wf%jacobian_transpose_cc2_a1(equation, wf%t1bar)
!
!     equation += sum_bjck tbar_bjck A_{bjck,ai}
!
      call wf%jacobian_transpose_cc2_b1(equation, t2bar)
!
      call mem%dealloc(t2bar, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Add eta, equation = t-bar^T A + eta 
!
      call mem%alloc(eta, wf%n_gs_amplitudes)
      call wf%construct_eta(eta)
!
      call daxpy(wf%n_gs_amplitudes, one, eta, 1, equation, 1)
!
      call mem%dealloc(eta, wf%n_gs_amplitudes)
!
   end subroutine construct_multiplier_equation_cc2
!
!
   subroutine get_cvs_projector_cc2(wf, projector, n_cores, core_MOs)
!!
!!    Get CVS projector
!!    Written by Sarai D. Folkestad, Oct 2018
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: projector
!
      integer, intent(in) :: n_cores
!
      integer, dimension(n_cores), intent(in) :: core_MOs
!
      integer :: core, i, a, ai, j, b, bj, aibj
!
      projector = zero
!
      do core = 1, n_cores
!
        i = core_MOs(core)
!
!$omp parallel do private (a, ai, j, b, bj, aibj)
        do a = 1, wf%n_v
!
           ai = wf%n_v*(i - 1) + a
           projector(ai) = one
!
            do j = 1, wf%n_o 
               do b = 1, wf%n_v
!
                  bj = wf%n_v*(j - 1) + b
                  aibj = max(ai, bj)*(max(ai, bj) - 3)/2 + ai + bj
!                  
                  projector(aibj + (wf%n_o)*(wf%n_v)) = one
!
               enddo
            enddo
        enddo
!$omp end parallel do
!
     enddo
!
   end subroutine get_cvs_projector_cc2
!
!
   subroutine save_doubles_vector_cc2(wf, X, n, file_)
!!
!!    Save doubles vector state 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Mar 2019 
!!
!!    Writes doubles vector "X" to the sequential
!!    and unformatted file "file_".
!!    
!!    NB! If n = 1, then the routine WILL REWIND the file before writing,
!!    thus DELETING every record in the file. For n >=2, we just append to
!!    the file. The purpose of this setup is that the files should be saved in 
!!    the correct order, from n = 1 to n = # states.
!!
      implicit none 
!
      class(cc2), intent(inout) :: wf 
!
      real(dp), dimension(wf%n_t2), intent(in) :: X 
!
      integer, intent(in) :: n ! state number 
!
      type(file) :: file_
!
      call disk%open_file(file_, 'write', 'append')
!
      if (n .eq. 1) rewind(file_%unit)
!
      write(file_%unit) X
!
      call disk%close_file(file_, 'keep')
!
   end subroutine save_doubles_vector_cc2
!
!
   subroutine read_doubles_vector_cc2(wf, X, n, file_)
!!
!!    Read doubles vector state 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Mar 2019 
!!
!!    Reads doubles vector "X" from the "n"'th line
!!    of the sequential and unformatted file "file_".
!!
      implicit none 
!
      class(cc2), intent(inout) :: wf 
!
      real(dp), dimension(wf%n_t2), intent(out) :: X 
!
      integer, intent(in) :: n ! state number 
!
      type(file) :: file_
!
      call disk%open_file(file_, 'read')
!
      call file_%prepare_to_read_line(n)
!
      read(file_%unit) X
!
      call disk%close_file(file_, 'keep')
!
   end subroutine read_doubles_vector_cc2
!
!
   subroutine save_excited_state_cc2(wf, X, n, side)
!!
!!    Save excited state 
!!    Written by Eirik F. Kjønstad, Mar 2019 
!!
!!    Saves an excited state to disk. Since the solvers 
!!    keep these vectors in full length, we receive a vector 
!!    in full length (n_es_amplitudes), and then distribute 
!!    the different parts of that vector to singles, doubles, etc.,
!!    files (if there are doubles, etc.).
!!
!!    NB! If n = 1, then the routine WILL REWIND the files before writing,
!!    thus DELETING every record in the file. For n >=2, we just append to
!!    the file. The purpose of this setup is that the files should be saved in 
!!    the correct order, from n = 1 to n = # states. 
!!
      implicit none 
!
      class(cc2), intent(inout) :: wf 
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: X 
!
      integer, intent(in) :: n ! state number 
!
      character(len=*), intent(in) :: side ! 'left' or 'right' 
!
      if (trim(side) == 'right') then 
!
         call wf%save_singles_vector(X(1 : wf%n_t1), n, wf%r1_file)
         call wf%save_doubles_vector(X(wf%n_t1 + 1 : wf%n_es_amplitudes), n, wf%r2_file)
!
      elseif (trim(side) == 'left') then 
!
         call wf%save_singles_vector(X(1 : wf%n_t1), n, wf%l1_file)
         call wf%save_doubles_vector(X(wf%n_t1 + 1 : wf%n_es_amplitudes), n, wf%l2_file)
!
      else
!
         call output%error_msg('Tried to save an excited state, but argument side not recognized: ' // side)
!
      endif
!
   end subroutine save_excited_state_cc2
!
!
   subroutine restart_excited_state_cc2(wf, X, n, side)
!!
!!    Restart excited state 
!!    Written by Sarai D. Fokestad, Mar 2019 
!!
!!    Wrapper for setting trial vectors to excited states on file
!!
      implicit none 
!
      class(cc2), intent(inout) :: wf 
!
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: X 
!
      integer, intent(in) :: n ! state number 
!
      character(len=*), intent(in) :: side ! 'left' or 'right' 
!
      integer :: a, i, b, j, ai, bj, aibj
      integer :: n_excited_states
!
      real(dp), dimension(:,:,:,:), allocatable :: r2_aibj, l2_aibj
!
      real(dp), dimension(:), allocatable :: omega
!
      call wf%is_restart_safe('excited state')
!
!     Check if we have read doubles vectors.
!     If not, set up doubles.
!
      if (trim(side) == 'right') then 
!
!        Read singles vector
!
         call wf%read_singles_vector(X(1 : wf%n_t1), n, wf%r1_file)
!
!        Read or construct doubles vector
!
         if (wf%r2_file%exists()) then
!
            call wf%read_doubles_vector(X(wf%n_t1 + 1 : wf%n_es_amplitudes), n, wf%r2_file)
!
         else
!
!           Construct r2 from r1
!
!           r2_aibj = (A_aibj,ck r1_ck)/(- ε_aibj + omega)
!
            call mem%alloc(r2_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
            r2_aibj = zero
!
            call wf%jacobian_cc2_a2(r2_aibj, X(1:(wf%n_o)*(wf%n_v)))
!
            n_excited_states = wf%get_n_excitation_energies_on_file()
!
            call mem%alloc(omega, n_excited_states)
!
            call wf%read_excitation_energies(n_excited_states, omega)
!
!$omp parallel do private(a, i, b, j, aibj, ai, bj)
            do a = 1, wf%n_v
               do i = 1, wf%n_o 
!
                  ai = wf%n_v*(i - 1) + a
!
                  do b = 1, wf%n_v
                     do j = 1, wf%n_o
!
                        bj = wf%n_v*(j - 1) + b
!
                        if (ai .ge. bj) then
!
                           aibj = ai*(ai - 3)/2 + ai + bj
!
                           X(aibj + (wf%n_v)*(wf%n_o)) = r2_aibj(a, i, b, j)&
                                                     /(omega(n) - wf%orbital_energies(a + wf%n_o) &
                                                                - wf%orbital_energies(b + wf%n_o) &
                                                                + wf%orbital_energies(i) &
                                                                + wf%orbital_energies(j) )
!
                        endif
!
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
            call mem%dealloc(omega, n_excited_states)
            call mem%dealloc(r2_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
         endif
!
      elseif (trim(side) == 'left') then
!
!        Read singles vector
!
         call wf%read_singles_vector(X(1 : wf%n_t1), n, wf%l1_file)
!
!        Read or construct doubles vector
!
         if (wf%l2_file%exists()) then
!
            call wf%read_doubles_vector(X(wf%n_t1 + 1 : wf%n_es_amplitudes), n, wf%l2_file)
!
         else
!
!           Construct l2 from l1
!
!           r2_aibj = (A^T_aibj,ck r1_ck)/(- ε_aibj + omega)
!
            call mem%alloc(l2_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
            l2_aibj = zero
!
            call wf%jacobian_transpose_cc2_a2(l2_aibj, X(1:(wf%n_o)*(wf%n_v)))
!
            n_excited_states = wf%get_n_excitation_energies_on_file()
!
            call mem%alloc(omega, n_excited_states)
!
            call wf%read_excitation_energies(n_excited_states, omega)
!
!$omp parallel do private(a, i, b, j, aibj, ai, bj)
            do b = 1, wf%n_v
               do j = 1, wf%n_o 
!
                  bj = wf%n_v*(j - 1) + b
!
                  do i = 1, wf%n_o
                     do a = 1, wf%n_v
!
                        ai = wf%n_v*(i - 1) + a
!
                        if (ai .ge. bj) then
!
                           aibj = ai*(ai - 3)/2 + ai + bj
!
                           X(aibj + (wf%n_v)*(wf%n_o)) = l2_aibj(a, i, b, j)&
                                                     /(omega(n) - wf%orbital_energies(a + wf%n_o) &
                                                                - wf%orbital_energies(b + wf%n_o) &
                                                                + wf%orbital_energies(i) &
                                                                + wf%orbital_energies(j) )
!
                        endif
!
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
            call mem%dealloc(omega, n_excited_states)
            call mem%dealloc(l2_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
         endif
!
      endif
!
   end subroutine restart_excited_state_cc2
!
!
   subroutine is_restart_safe_cc2(wf, task)
!!
!!    Is restart safe?
!!    Written by Eirik F. Kjønstad, Mar 2019 
!!
      implicit none 
!
      class(cc2) :: wf 
!
      character(len=*), intent(in) :: task 
!
      integer :: n_o, n_v, n_gs_amplitudes, n_es_amplitudes
!
      call disk%open_file(wf%restart_file, 'read', 'rewind')
!
      read(wf%restart_file%unit) n_o
      read(wf%restart_file%unit) n_v
      read(wf%restart_file%unit) n_gs_amplitudes
      read(wf%restart_file%unit) n_es_amplitudes
!
      call disk%close_file(wf%restart_file)
!
      if (n_o .ne. wf%n_o) call output%error_msg('attempted to restart from inconsistent number ' // &
                                                   'of occupied orbitals.')
!
      if (n_v .ne. wf%n_v) call output%error_msg('attempted to restart from inconsistent number ' // &
                                                   'of virtual orbitals.')
!
      if (trim(task) == 'ground state') then 
!
         if (n_gs_amplitudes .ne. wf%n_gs_amplitudes) &
            call output%error_msg('attempted to restart from inconsistent number ' // &
                                    'of ground state amplitudes.')    
!
      elseif (trim(task) == 'excited state') then    
!
         if (n_es_amplitudes .eq. wf%n_t1) then 
!
!           OK! We allow restart from CCS-like models (e.g. CC2 lowmem or CCS itself) in (highmem) CC2.
!           
         elseif (n_es_amplitudes .ne. wf%n_es_amplitudes) then
!
            call output%error_msg('attempted to restart from inconsistent number ' // &
                                    'of excited state amplitudes.')     
!
         endif
!
      else
!
         call output%error_msg('attempted to restart, but the task was not recognized: ' // task)
!
      endif   
!
   end subroutine is_restart_safe_cc2
!
!
   subroutine initialize_files_cc2(wf)
!!
!!    Initialize files 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
!!    Initializes the wavefucntion files for wavefunction parameters.
!!
      class(cc2) :: wf 
!
      call wf%initialize_wavefunction_files()
      call wf%initialize_cc_files()
      call wf%initialize_singles_files()
      call wf%initialize_doubles_files()
!
   end subroutine initialize_files_cc2
!
!
   subroutine initialize_doubles_files_cc2(wf)
!!
!!    Initialize doubles files 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
      class(cc2) :: wf 
!
      call wf%l2_file%init('l2', 'sequential', 'unformatted')
!
      call wf%r2_file%init('r2', 'sequential', 'unformatted')
!
   end subroutine initialize_doubles_files_cc2
!
!
!
end module cc2_class
