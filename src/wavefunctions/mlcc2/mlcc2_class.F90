!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2020 the authors of eT
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
module mlcc2_class
!
!!
!!    Multilevel Coupled cluster singles and perturbative doubles (MLCC2) class module  
!!    Written by Sarai D. Folkestad, June 2017 and spring 2019
!!
!!    This MLCC2 wavefunction is given by
!!
!!       | MLCC2 > = exp(T_1 + S_2) | HF >
!!
!!    where T_1 is the single excitation operator from standard CC
!!    and S_2 is a double excitation operator which only contains
!!    excitations within an active orbital space (the CC2 orbitals)
!!
!!    The S_2 operator is determined to first order in the 
!!    fluctuation potential U (H = F + U).
!!
!!    This class handles the ground and excited states for 
!!    MLCC2.
!!
!!    For further references on MLCC see:
!!
!!       Myhre, R. H., & Koch, H., JCP, 145(4), 044111 (2016)
!!
!
   use ccs_class
!
   implicit none
!
   type, extends(ccs) :: mlcc2
!
!     Requested levels
!
      logical :: do_ccs
      logical :: do_cc2
!
!     Orbital indices for levels
!
      integer :: n_ccs_o
      integer :: n_ccs_v
!
      integer :: n_cc2_o
      integer :: n_cc2_v
!
      integer :: first_ccs_o
      integer :: first_ccs_v
      integer :: last_ccs_o
      integer :: last_ccs_v
!
      integer :: first_cc2_o
      integer :: first_cc2_v
      integer :: last_cc2_o
      integer :: last_cc2_v
!
!     Orbital type string
!
      character(len=200) :: cc2_orbital_type
!
!     cc2 variables
!
      integer :: n_x2 ! n_s2
!
      real(dp), dimension(:,:,:,:), allocatable :: u_aibj
!
!     For properties
!
      real(dp), dimension(:), allocatable :: x2
      real(dp), dimension(:), allocatable :: t2bar
!
      integer :: n_nto_states
      integer :: n_cnto_states
!
      integer, dimension(:), allocatable :: cnto_states
      integer, dimension(:), allocatable :: nto_states
!
      logical :: cnto_restart

      type(sequential_file) :: jacobian_a1_intermediate_vv
      type(sequential_file) :: jacobian_a1_intermediate_oo
!
   contains
!
      procedure :: print_amplitude_info                              => print_amplitude_info_mlcc2
!
      procedure :: cleanup                                           => cleanup_mlcc2
!
!     Initializations and destructions
!
      procedure :: initialize_u_aibj                                 => initialize_u_aibj_mlcc2
      procedure :: destruct_u_aibj                                   => destruct_u_aibj_mlcc2
!
      procedure :: initialize_amplitudes                             => initialize_amplitudes_mlcc2
      procedure :: destruct_amplitudes                               => destruct_amplitudes_mlcc2
!
      procedure :: initialize_x2                                     => initialize_x2_mlcc2
      procedure :: destruct_x2                                       => destruct_x2_mlcc2
!
      procedure :: initialize_t2bar                                  => initialize_t2bar_mlcc2
      procedure :: destruct_t2bar                                    => destruct_t2bar_mlcc2
      procedure :: destruct_multipliers                              => destruct_multipliers_mlcc2
!
      procedure :: initialize_nto_states                             => initialize_nto_states_mlcc2
      procedure :: destruct_nto_states                               => destruct_nto_states_mlcc2
!
      procedure :: initialize_cnto_states                            => initialize_cnto_states_mlcc2
      procedure :: destruct_cnto_states                              => destruct_cnto_states_mlcc2
!
!     Set and get
!
      procedure :: get_es_orbital_differences                        => get_es_orbital_differences_mlcc2
!
      procedure :: determine_n_x2_amplitudes                         => determine_n_x2_amplitudes_mlcc2
      procedure :: determine_n_es_amplitudes                         => determine_n_es_amplitudes_mlcc2
      procedure :: determine_n_gs_amplitudes                         => determine_n_gs_amplitudes_mlcc2
!
      procedure :: get_cvs_projector                                 => get_cvs_projector_mlcc2
      procedure :: set_cvs_start_indices                             => set_cvs_start_indices_mlcc2
!
!     File handling
!
      procedure :: read_mlcc_settings                                => read_mlcc_settings_mlcc2
      procedure, non_overridable :: read_orbital_settings            => read_orbital_settings_mlcc2
      procedure, non_overridable :: read_cc2_orbital_settings        => read_cc2_orbital_settings_mlcc2
!
!     Orbital routines
!
      procedure :: mo_preparations                                   => mo_preparations_mlcc2
!
      procedure :: print_orbital_space                               => print_orbital_space_mlcc2
      procedure :: check_orbital_space                               => check_orbital_space_mlcc2
      procedure :: check_orthonormality_of_MOs                       => check_orthonormality_of_MOs_mlcc2
!
      procedure :: read_cnto_transformation_matrices                 => read_cnto_transformation_matrices_mlcc2
      procedure :: construct_ccs_cnto_transformation_matrices        => construct_ccs_cnto_transformation_matrices_mlcc2
!
      procedure :: construct_cholesky_orbitals                       => construct_cholesky_orbitals_mlcc2
!
      procedure :: construct_block_diagonal_fock_orbitals            => construct_block_diagonal_fock_orbitals_mlcc2
!
      procedure :: construct_block_diagonal_fock_mos_2_level         => construct_block_diagonal_fock_mos_2_level_mlcc2
!
      procedure :: construct_cntos                                   => construct_cntos_mlcc2
      procedure :: construct_M_and_N_cnto                            => construct_M_and_N_cnto_mlcc2
      procedure :: ccs_calculation_for_cntos                         => ccs_calculation_for_cntos_mlcc2
!
      procedure :: construct_M_nto                                   => construct_M_nto_mlcc2
      procedure :: construct_ccs_nto_transformation_matrix           => construct_ccs_nto_transformation_matrix_mlcc2
      procedure :: construct_mixed_nto_canonical_orbitals            => construct_mixed_nto_canonical_orbitals_mlcc2
      procedure :: construct_paos                                    => construct_paos_mlcc2
!
      procedure :: orbital_partitioning                              => orbital_partitioning_mlcc2
      procedure :: diagonalize_M_and_N                               => diagonalize_M_and_N_mlcc2
!
      procedure :: contruct_mo_basis_transformation                  => contruct_mo_basis_transformation_mlcc2
      procedure :: update_MO_fock_contributions                      => update_MO_fock_contributions_mlcc2
!
!     Ground state routines
!
      procedure :: calculate_energy                                  => calculate_energy_mlcc2
      procedure :: construct_omega                                   => construct_omega_mlcc2
!
      procedure :: omega_cc2_a1                                      => omega_cc2_a1_mlcc2
      procedure :: omega_cc2_b1                                      => omega_cc2_b1_mlcc2
      procedure :: omega_cc2_c1                                      => omega_cc2_c1_mlcc2
!
      procedure :: construct_multiplier_equation                     => construct_multiplier_equation_mlcc2
!
!     Jacobian transformation routines
!
      procedure :: prepare_for_jacobian                              => prepare_for_jacobian_mlcc2
!
      procedure :: jacobian_transformation                           => jacobian_transformation_mlcc2
!
      procedure :: jacobian_cc2_a1                                   => jacobian_cc2_a1_mlcc2
      procedure :: jacobian_cc2_b1                                   => jacobian_cc2_b1_mlcc2
      procedure :: jacobian_cc2_a2                                   => jacobian_cc2_a2_mlcc2
      procedure :: jacobian_cc2_b2                                   => jacobian_cc2_b2_mlcc2
!
      procedure :: save_jacobian_a1_intermediates                    => save_jacobian_a1_intermediates_mlcc2
!
!     Jacobian transpose routines
!
      procedure :: jacobian_transpose_transformation                 => jacobian_transpose_transformation_mlcc2
!
      procedure :: jacobian_transpose_cc2_a1                         => jacobian_transpose_cc2_a1_mlcc2
      procedure :: jacobian_transpose_cc2_b1                         => jacobian_transpose_cc2_b1_mlcc2
      procedure :: jacobian_transpose_cc2_a2                         => jacobian_transpose_cc2_a2_mlcc2
      procedure :: jacobian_transpose_cc2_b2                         => jacobian_transpose_cc2_b2_mlcc2
!
!     Amplitudes
!
      procedure :: construct_x2                                      => construct_x2_mlcc2
      procedure :: construct_t2bar                                   => construct_t2bar_mlcc2
      procedure :: construct_u_aibj                                  => construct_u_aibj_mlcc2
!
!     Debug 
!
      procedure :: omega_for_jacobian_debug                          => omega_for_jacobian_debug_mlcc2
      procedure :: amplitudes_for_jacobian_debug                     => amplitudes_for_jacobian_debug_mlcc2
      procedure :: normalization_for_jacobian_debug                  => normalization_for_jacobian_debug_mlcc2
      procedure :: construct_omega_doubles                           => construct_omega_doubles_mlcc2
!
!     Restart
!
      procedure :: is_restart_safe                                   => is_restart_safe_mlcc2
      procedure :: write_cc_restart                                  => write_cc_restart_mlcc2
!
   end type mlcc2
!
   interface mlcc2
!
      procedure :: new_mlcc2
!
   end interface mlcc2
!
   interface
!
      include "./orbitals_mlcc2_interface.F90"
      include "./omega_mlcc2_interface.F90"
      include "./jacobian_mlcc2_interface.F90"
      include "./jacobian_transpose_mlcc2_interface.F90"
      include "./initialize_destruct_mlcc2_interface.F90"
      include "./debug_jacobian_mlcc2_interface.F90"
!
   end interface 
!
contains
!
!
   function new_mlcc2(system, template_wf) result(wf)
!!
!!    New mlcc2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Adapted by Sarai D. Folkestad from CCS constructer, 2019
!!
      implicit none
!
      type(mlcc2) :: wf
!
      class(molecular_system), target, intent(in) :: system 
!
      class(wavefunction) :: template_wf
!
      wf%name_ = 'mlcc2'
!
      wf%n_cc2_o = 0
      wf%n_cc2_v = 0
      wf%n_ccs_o = 0
      wf%n_ccs_v = 0
!
      wf%need_g_abcd = .false.
!
      wf%cholesky_orbital_threshold = 1.0d-2
!
      call wf%general_cc_preparations(system)
      call wf%set_variables_from_template_wf(template_wf)
      call wf%print_banner()
!
      if (wf%bath_orbital) call output%error_msg('Bath orbitals not yet implemented for MLCC2')
!
      call wf%read_mlcc_settings()
!
      wf%n_t1 = (wf%n_o)*(wf%n_v)
      wf%n_gs_amplitudes = wf%n_t1
!
      call wf%initialize_fock()
!
      call wf%print_amplitude_info()
!
   end function new_mlcc2
!
!
   subroutine print_amplitude_info_mlcc2(wf)
!!
!!    Print amplitude info
!!    Written by Sarai D. Folkestad, Dec 2019
!!
!!
      implicit none
!
      class(mlcc2), intent(in) :: wf
!
      call wf%ccs%print_amplitude_info()  
!
   end subroutine print_amplitude_info_mlcc2
!
!
!
   subroutine print_orbital_space_mlcc2(wf)
!!
!!    Print orbital space
!!    Written by Sarai D. Folkestad, Jan 2019
!!
      implicit none
!
      class(mlcc2) :: wf
!
      call output%printf('m', '- MLCC2 orbital partitioning:',fs='(/t3,a)')
!
      call output%printf('m', 'Orbital type: ' // trim(wf%cc2_orbital_type), fs='(/t6,a)')
      call output%printf('m', 'Number occupied cc2 orbitals: (i4)', &
                         ints=[wf%n_cc2_o], fs='(/t6,a)')
      call output%printf('m', 'Number virtual cc2 orbitals:  (i4)', &
                         ints=[wf%n_cc2_v], fs='(t6,a)')
      call output%printf('m', 'Number occupied ccs orbitals: (i4)', &
                         ints=[wf%n_ccs_o], fs='(/t6,a)')
      call output%printf('m', 'Number virtual ccs orbitals:  (i4)', &
                         ints=[wf%n_ccs_v], fs='(t6,a)')
!
   end subroutine print_orbital_space_mlcc2
!
!
   subroutine read_mlcc_settings_mlcc2(wf)
!!
!!    Read MLCC settings
!!    Written by Sarai D. Folkestad, Apr 2019
!!
!!    Reads the cc and mlcc sections of the
!!    input file.
!!
!!    Calls two routines that handle 
!!    reading of orbital type etc. The routine will be
!!    overwritten for MLCCSD.
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      if (.not. input%requested_section('mlcc')) &
         call output%error_msg('cannot do mlcc calculation without mlcc section in eT.inp')
!
      call input%get_required_keyword_in_section('cc2 orbitals', 'mlcc', wf%cc2_orbital_type)
!
!     Get general orbital specifications
!
      call wf%read_orbital_settings(wf%cc2_orbital_type)
!
!     Read CC2 orbital specifications
!
      call wf%read_cc2_orbital_settings()
!
   end subroutine read_mlcc_settings_mlcc2
!
!
   subroutine read_cc2_orbital_settings_mlcc2(wf)
!!
!!    Read CC2 orbital settings
!!    Written by Sarai D. Folkestad, Apr 2019
!!
!!    Reads the CC2 specific orbital settings:
!!
!!       - reads the 'cc2 orbitals' keyword (requested keyword)
!!
!!    If CC2 orbitals are CNTOs or NTOs/Canonical:
!!
!!       - reads the number of occupied cntos and ntos (requested keyword for CNTO/NTO)
!!       - reads the number of virtual cntos and canonical
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      if (trim(wf%cc2_orbital_type) == 'cnto-approx') then
!
         call input%get_required_keyword_in_section('cnto occupied cc2', 'mlcc', wf%n_cc2_o)
!
         wf%n_cc2_v = wf%n_cc2_o*(wf%n_v/wf%n_o)
!
         call input%get_keyword_in_section('cnto virtual cc2', 'mlcc', wf%n_cc2_v)
!
      elseif (trim(wf%cc2_orbital_type) == 'nto-canonical') then
!
         call input%get_required_keyword_in_section('nto occupied cc2', 'mlcc', wf%n_cc2_o)
!
         wf%n_cc2_v = wf%n_cc2_o*(wf%n_v/wf%n_o)
!
         call input%get_keyword_in_section('canonical virtual cc2', 'mlcc', wf%n_cc2_v)
!
      endif
!
   end subroutine read_cc2_orbital_settings_mlcc2
!
!
   subroutine read_orbital_settings_mlcc2(wf, orbital_type)
!!
!!    Read orbital settings
!!    Written by Sarai D. Folkestad, Apr 2019
!!
!!    Reads the general, not the cc level specific,
!!    settings for the given 'orbital_type'
!!
!!    The routine will also be used for mlccsd
!!
!!    CNTO settings:
!!
!!       - 'cnto restart'
!!       - 'cnto states' (requested)
!!    
!!     NTO settings:
!!
!!       - 'nto states' (requested)
!!   
!!    Cholesky and Cholesky/PAO settings:
!!
!!       - 'cholesky threshold'
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      character(len=*), intent(in) :: orbital_type
!
      if (trim(orbital_type) == 'cnto' .or. trim(orbital_type) == 'cnto-approx') then
!
         wf%cnto_restart = .false.
!
         if (input%requested_keyword_in_section('cnto restart', 'mlcc')) wf%cnto_restart = .true.
!
         wf%n_cnto_states = input%get_n_elements_for_keyword_in_section('cnto states', 'mlcc')
!
         if (wf%n_cnto_states == 0) &
               call output%error_msg('to construct CNTOs excitation vectors must be specified.')
!
         call wf%initialize_cnto_states()
         call input%get_array_for_keyword_in_section('cnto states', 'mlcc', &
                  wf%n_cnto_states, wf%cnto_states)
!
      elseif (trim(orbital_type) == 'cholesky') then
!
         call input%get_keyword_in_section('cholesky threshold', 'mlcc', &
                  wf%cholesky_orbital_threshold)
!
      elseif (trim(orbital_type) == 'nto-canonical') then
!
         wf%n_nto_states = input%get_n_elements_for_keyword_in_section('nto states', 'mlcc')
!
         if (wf%n_nto_states == 0) &
               call output%error_msg('to construct NTOs excitation vectors must be specified.')
!
         call wf%initialize_nto_states()
         call input%get_array_for_keyword_in_section('nto states', 'mlcc', &
               wf%n_nto_states, wf%nto_states)
!
      elseif (trim(orbital_type) == 'cholesky-pao') then
!
         call input%get_keyword_in_section('cholesky threshold', 'mlcc', &
               wf%cholesky_orbital_threshold)
!
      else 
!
         call output%error_msg('could not recognize the orbital type')
!
      endif
!
   end subroutine read_orbital_settings_mlcc2
!
!
   subroutine construct_u_aibj_mlcc2(wf)
!!
!!    Construct u_aibj
!!    Written by Sarai D. Folkestad, Jan 2019
!!
!!    Construct
!!
!!       u_aibj = 2s_aibj - s_ajbi
!!
!!    with
!!
!!       s_aibj = - g_aibj/ε_aibj
!!
!!    where
!!
!!       ε_aibj = ε_a - ε_i + ε_b - ε_j 
!!
!!    and ε_r is the r'th orbital energy.
!!
!!    Assumes that x2 is already allocated 
!!    and constructed.
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      real(dp), dimension(:,:,:,:), allocatable :: s_aibj
!
      call mem%alloc(s_aibj, wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o)  
!
      call squareup(wf%x2, s_aibj, wf%n_cc2_v*wf%n_cc2_o)
!
      call copy_and_scale(two, s_aibj, wf%u_aibj, (wf%n_cc2_v**2)*(wf%n_cc2_o**2))
      call add_1432_to_1234(-one, s_aibj, wf%u_aibj, wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o)
!    
      call mem%dealloc(s_aibj, wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o)     
!
   end subroutine construct_u_aibj_mlcc2
!
!
   subroutine calculate_energy_mlcc2(wf)
!!
!!    Calculate energy 
!!    Written by Sarai D. Folkestad, Feb 2019
!!
!!       E = E_HF + sum_aibj t_i^a*t_j^b L_iajb + sum_aibj s_ij^ab L_iajb
!!
!!    with
!!
!!       s_aibj = - g_aibj/ε_aibj, {a,i,b,j} are CC2 orbitals
!! 
      class(mlcc2), intent(inout) :: wf 
!
      real(dp), dimension(:,:,:,:), allocatable :: g_aibj, g_iajb 
!
      real(dp) :: correlation_energy
!
      integer :: a, i, b, j
!
      call mem%alloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_ovov(g_iajb)
!
      correlation_energy = zero 
!
!     t1-contribution
!
!$omp parallel do private(a,i,b,j) reduction(+:correlation_energy)
      do b = 1, wf%n_v
         do i = 1, wf%n_o 
            do j = 1, wf%n_o 
               do a = 1, wf%n_v
!
                  correlation_energy = correlation_energy + &
                        (wf%t1(a, i)*wf%t1(b, j)*(two*g_iajb(i,a,j,b)-g_iajb(i,b,j,a)))
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%alloc(g_iajb, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v)
      call mem%alloc(g_aibj, wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o)
!
      call wf%get_ovov(g_iajb, &
                        wf%first_cc2_o, wf%last_cc2_o, &
                        wf%first_cc2_v, wf%last_cc2_v, &
                        wf%first_cc2_o, wf%last_cc2_o, &
                        wf%first_cc2_v, wf%last_cc2_v)
!
      call wf%get_vovo(g_aibj, &
                        wf%first_cc2_v, wf%last_cc2_v, &
                        wf%first_cc2_o, wf%last_cc2_o, &
                        wf%first_cc2_v, wf%last_cc2_v, &
                        wf%first_cc2_o, wf%last_cc2_o)
!
!$omp parallel do private(a,i,b,j) reduction(+:correlation_energy)
      do b = 1, wf%n_cc2_v
         do i = 1, wf%n_cc2_o 
            do j = 1, wf%n_cc2_o 
               do a = 1, wf%n_cc2_v
!
                  correlation_energy = correlation_energy - (g_aibj(a,i,b,j)/              &
                                    (wf%orbital_energies(wf%n_o + a + wf%first_cc2_v - 1)  &
                                   + wf%orbital_energies(wf%n_o + b + wf%first_cc2_v - 1)  &
                                   - wf%orbital_energies(i + wf%first_cc2_o - 1)           &
                                   - wf%orbital_energies(j + wf%first_cc2_o - 1)))         &
                                       *(two*g_iajb(i,a,j,b)-g_iajb(i,b,j,a))
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do 
!
      call mem%dealloc(g_iajb, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v)
      call mem%dealloc(g_aibj, wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o)
!
      wf%correlation_energy = correlation_energy
!
      wf%energy = wf%hf_energy + correlation_energy
!
   end subroutine calculate_energy_mlcc2
!
!
   subroutine check_orbital_space_mlcc2(wf)
!!
!!    Check orbital space
!!    Written by Sarai D. Folkestad, Jan 2019
!!
!!    Some sanity checks to ensure that nothing is very wrong
!!    with the orbital partitioning
!!
      implicit none
!
      class(mlcc2) :: wf
!
      call wf%check_orthonormality_of_MOs()
!
      if (wf%n_cc2_o .eq. 0) &
               call output%error_msg('no occupied cc2 orbitals in mlcc2 calculation.')
!
      if (wf%n_cc2_v .eq. 0) &
               call output%error_msg('no virtual cc2 orbitals in mlcc2 calculation.')
!
      if (wf%n_ccs_o + wf%n_ccs_v .eq. 0)  &
            call output%warning_msg('no ccs orbitals in mlcc2 calculation, ' //&
               'recomended to run standard cc2 code.')
!
      if (wf%n_cc2_o .lt. 0 .or. wf%n_ccs_o .lt. 0 .or. &
          wf%n_cc2_v .lt. 0 .or. wf%n_ccs_v .lt. 0 ) &
            call output%error_msg('negative orbital space size.')
!
      if (wf%n_cc2_o + wf%n_ccs_o .ne. wf%n_o) &
            call output%error_msg('occupied spaces do not add to n_o.')
!
      if (wf%n_cc2_v + wf%n_ccs_v .ne. wf%n_v) &
            call output%error_msg('virtual spaces do not add to n_v.')
!
   end subroutine check_orbital_space_mlcc2
!
!
   subroutine get_es_orbital_differences_mlcc2(wf, orbital_differences, N)
!!
!!    Get orbital differences 
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 2019
!!
      implicit none
!
      class(mlcc2), intent(in) :: wf
!
      integer, intent(in) :: N 
      real(dp), dimension(N), intent(inout) :: orbital_differences
!
      integer :: a, i, ai, b, j, bj, aibj
!
!$omp parallel do schedule(static) private(a, i, ai) 
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            ai = wf%n_v*(i - 1) + a
!
            orbital_differences(ai) = wf%orbital_energies(a + wf%n_o) - wf%orbital_energies(i)
!
         enddo
      enddo
!$omp end parallel do
!
!$omp parallel do schedule(static) private(a, i, b, j, ai, bj, aibj) 
      do a = 1, wf%n_cc2_v
         do i = 1, wf%n_cc2_o
!
            ai = wf%n_cc2_v*(i - 1) + a
!
            do j = 1, wf%n_cc2_o 
               do b = 1, wf%n_cc2_v
!
                  bj = wf%n_cc2_v*(j-1) + b 
!
                  if (ai .ge. bj) then
!
                     aibj = (ai*(ai-3)/2) + ai + bj
!
                     orbital_differences(aibj + (wf%n_o)*(wf%n_v)) =  &
                                 wf%orbital_energies(a + wf%n_o + wf%first_cc2_v - 1) &
                                 - wf%orbital_energies(i + wf%first_cc2_o - 1) &
                                 + wf%orbital_energies(b + wf%n_o + wf%first_cc2_v - 1) &
                                 - wf%orbital_energies(j + wf%first_cc2_o - 1)
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
   end subroutine get_es_orbital_differences_mlcc2
!
!
   subroutine construct_multiplier_equation_mlcc2(wf, equation)
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
      class(mlcc2), intent(in) :: wf 
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
      call mem%alloc(t2bar, wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o)
      call zero_array(t2bar, (wf%n_cc2_v**2)*(wf%n_cc2_o)**2)
!
!     t2bar = sum_ai tbar_ai A_ai,aibj
!
      call wf%jacobian_transpose_cc2_a2(t2bar, wf%t1bar, wf%n_cc2_o, wf%n_cc2_v, &
                                       wf%first_cc2_o, wf%first_cc2_v, wf%last_cc2_o, wf%last_cc2_v)
!
      call symmetric_sum(t2bar, wf%n_cc2_o*wf%n_cc2_v)
!
      call mem%alloc(g_iajb, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v)
!      
      call wf%get_ovov(g_iajb,                           &
                        wf%first_cc2_o, wf%last_cc2_o,   &
                        wf%first_cc2_v, wf%last_cc2_v,   &
                        wf%first_cc2_o, wf%last_cc2_o,   &
                        wf%first_cc2_v, wf%last_cc2_v)
!
!     t2bar += η_aibj
!
      call add_2143_to_1234(four, g_iajb, t2bar, wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o)
      call add_2341_to_1234(-two, g_iajb, t2bar, wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o)
!
      call mem%dealloc(g_iajb, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v)
!
!     t2bar = t2bar/(-ε_aibj)
!
!$omp parallel do private(a, b, i, j)
      do b = 1, wf%n_cc2_v
         do j = 1, wf%n_cc2_o
            do i = 1, wf%n_cc2_o
               do a = 1, wf%n_cc2_v
!
                  t2bar(a, i, b, j) = t2bar(a, i, b, j)/&
                     (-  wf%orbital_energies(a + wf%n_o + wf%first_cc2_v - 1) &
                      -  wf%orbital_energies(b + wf%n_o + wf%first_cc2_v - 1) &
                      +  wf%orbital_energies(i + wf%first_cc2_o - 1) &
                      +  wf%orbital_energies(j + wf%first_cc2_o - 1))
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
      call zero_array(equation, wf%n_gs_amplitudes)
!  
      call wf%jacobian_transpose_ccs_a1(equation, wf%t1bar)
      call wf%jacobian_transpose_ccs_b1(equation, wf%t1bar)
      call wf%jacobian_transpose_cc2_a1(equation, wf%t1bar, wf%n_cc2_o, wf%n_cc2_v, &
                                 wf%first_cc2_o, wf%first_cc2_v)
!
!
!     equation += sum_bjck tbar_bjck A_{bjck,ai}
!
      call wf%jacobian_transpose_cc2_b1(equation, t2bar, wf%n_cc2_o, wf%n_cc2_v, &
                           wf%first_cc2_o, wf%first_cc2_v, wf%last_cc2_o, wf%last_cc2_v)
!
      call mem%dealloc(t2bar, wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o)
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
   end subroutine construct_multiplier_equation_mlcc2
!
!
   subroutine construct_x2_mlcc2(wf)
!!
!!    Construct x2 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Jan 2019
!!
!!    Construct
!!
!!       s_aibj = - g_aibj/ε_aibj
!!
!!    where
!!
!!       ε_aibj = ε_a - ε_i + ε_b - ε_j 
!!
!!    and ε_r is the r'th orbital energy.
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      real(dp), dimension(:,:,:,:), allocatable :: g_aibj, s_aibj
!
      integer :: a, i, b, j
!
      call mem%alloc(g_aibj, wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o)  
      call mem%alloc(s_aibj, wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o)  
!
      call wf%get_vovo(g_aibj,                           &
                        wf%first_cc2_v, wf%last_cc2_v,   &
                        wf%first_cc2_o, wf%last_cc2_o,   &
                        wf%first_cc2_v, wf%last_cc2_v,   &
                        wf%first_cc2_o, wf%last_cc2_o)
!
!$omp parallel do private(a, i, b, j)
      do b = 1, wf%n_cc2_v 
         do j = 1, wf%n_cc2_o 
            do i = 1, wf%n_cc2_o
               do a = 1, wf%n_cc2_v
!
                  s_aibj(a, i, b, j) = (g_aibj(a, i, b, j))/ &
                           (-  wf%orbital_energies(a + wf%n_o + wf%first_cc2_v - 1) &
                            -  wf%orbital_energies(b + wf%n_o + wf%first_cc2_v - 1) &
                            +  wf%orbital_energies(i + wf%first_cc2_o - 1) &
                            +  wf%orbital_energies(j + wf%first_cc2_o - 1))

!
               enddo
            enddo
         enddo 
      enddo
!$omp end parallel do
!
      call packin(wf%x2, s_aibj, wf%n_cc2_o*wf%n_cc2_v)    
!
      call mem%dealloc(g_aibj, wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o)      
      call mem%dealloc(s_aibj, wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o)      
!
   end subroutine construct_x2_mlcc2
!
!
   subroutine construct_t2bar_mlcc2(wf)
!!
!!    Construct t2bar
!!    Written by Sarai D. Folkestad, May, 2019
!!
      implicit none
!
      class(mlcc2) :: wf
!
      real(dp), dimension(:,:,:,:), allocatable :: t2bar, g_iajb
!
      integer :: a, i, b, j
!
      call mem%alloc(t2bar, wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o)
      call zero_array(t2bar, (wf%n_cc2_v**2)*(wf%n_cc2_o**2))
!
!     t2bar = sum_ai tbar_ai A_ai,aibj
!
      call wf%jacobian_transpose_cc2_a2(t2bar, wf%t1bar, wf%n_cc2_o, wf%n_cc2_v, &
                                       wf%first_cc2_o, wf%first_cc2_v, wf%last_cc2_o, wf%last_cc2_v)
!
      call symmetric_sum(t2bar, wf%n_cc2_o*wf%n_cc2_v)
!
      call mem%alloc(g_iajb, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v)
!
      call wf%get_ovov(g_iajb,                           &
                        wf%first_cc2_o, wf%last_cc2_o,   &
                        wf%first_cc2_v, wf%last_cc2_v,   &
                        wf%first_cc2_o, wf%last_cc2_o,   &
                        wf%first_cc2_v, wf%last_cc2_v)
!
!     t2bar += η_aibj
!
      call add_2143_to_1234(four, g_iajb, t2bar, wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o)
      call add_2341_to_1234(-two, g_iajb, t2bar, wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o)
!
      call mem%dealloc(g_iajb, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v)
!
!     t2bar = t2bar/(-ε_aibj)
!
!$omp parallel do private(a, b, i, j)
      do b = 1, wf%n_cc2_v
         do j = 1, wf%n_cc2_o
            do i = 1, wf%n_cc2_o
               do a = 1, wf%n_cc2_v
!
                  t2bar(a, i, b, j) = t2bar(a, i, b, j)/&
                           (-  wf%orbital_energies(a + wf%n_o + wf%first_cc2_v - 1) &
                            -  wf%orbital_energies(b + wf%n_o + wf%first_cc2_v - 1) &
                            +  wf%orbital_energies(i + wf%first_cc2_o - 1) &
                            +  wf%orbital_energies(j + wf%first_cc2_o - 1))
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call packin(wf%t2bar, t2bar, (wf%n_cc2_v)*(wf%n_cc2_o))
!
      call mem%dealloc(t2bar, wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o)
!
   end subroutine construct_t2bar_mlcc2
!
!
   subroutine get_cvs_projector_mlcc2(wf, projector, n_cores, core_MOs)
!!
!!    Get CVS projector
!!    Written by Sarai D. Folkestad, Oct 2018
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: projector
!
      integer, intent(in) :: n_cores
!
      integer, dimension(n_cores), intent(in) :: core_MOs
!
      integer :: core, i, a, ai, j, b, bj, aibj
!
      call zero_array(projector, wf%n_es_amplitudes)
!
      do core = 1, n_cores
!
         i = core_MOs(core)
!
         if ((i .lt. wf%first_cc2_o) .or. (i  .gt. wf%last_cc2_o)) then 
            call output%error_msg('Core orbital (i0) is not CC2 orbital', ints=[i])
         end if
!
!$omp parallel do private (a, ai)
         do a = 1, wf%n_v
!
           ai = wf%n_v*(i - 1) + a
           projector(ai) = one
!
         enddo
!$omp end parallel do
!
!$omp parallel do private (a, ai, j, b, bj, aibj)
         do a = 1, wf%n_cc2_v
!
            ai = wf%n_cc2_v*(i - 1) + a
!
            do j = 1, wf%n_cc2_o 
               do b = 1, wf%n_cc2_v
!
                  bj = wf%n_cc2_v*(j - 1) + b
!
                  aibj = max(ai, bj)*(max(ai, bj) - 3)/2 + ai + bj
!                  
                  projector(aibj + wf%n_o*wf%n_v) = one
!
               enddo
            enddo
        enddo
!$omp end parallel do
!
     enddo
!
   end subroutine get_cvs_projector_mlcc2
!
!
   subroutine set_cvs_start_indices_mlcc2(wf, start_indices)
!!
!!    Set CVS start indices
!!    Written by Sarai D. Folkestad
!!
      implicit none
!
      class(mlcc2), intent(in) :: wf
!
      integer, dimension(wf%n_singlet_states)  :: start_indices
!
!     Local variables
!
      integer ::  mo
!
!
!     Sanity check on core MOs -> Are they in active space?
!
      do mo = 1, wf%n_core_MOs
!
         if ((wf%core_MOs(mo) .lt. wf%first_cc2_o) &
         .or. (wf%core_MOs(mo) .gt. wf%last_cc2_o)) &
            call output%error_msg('Active MO not in active space for MLCC2 calculation.')
!
      enddo
!
      call wf%ccs%set_cvs_start_indices(start_indices)
!
   end subroutine set_cvs_start_indices_mlcc2
!
!
   subroutine mo_preparations_mlcc2(wf)
!!
!!    MO preparations
!!    Written by Sarai D. Folkestad, Sep 2019
!!
!!    Partitions the orbitals and
!!    determines the number of amplitudes.
!!
!!    Determines the MLCC basis (occupied-occupied and virtual-virtual 
!!    Fock matrices are block diagonal).
!!
!!    Prepares the MO Cholesky vectors
!!
!!    Transforms all frozen constributions to the Fock matrix
!!    from the old (canonical) MO basis to the MLCC basis.
!!    This update is done twice, once after orbital partitioning 
!!    and once after occupied-occupied and virtual-virtual 
!!    Fock matrices are block diagonalized.
!!
!!
      implicit none
!
      class(mlcc2) :: wf
!
      real(dp), dimension(:,:), allocatable :: canonical_orbitals
      real(dp), dimension(:,:), allocatable :: partitioning_orbitals
!
!     Keep canonical orbitals (for transformation of frozen MO fock terms)
!
      call mem%alloc(canonical_orbitals, wf%n_ao, wf%n_mo)
      call dcopy(wf%n_ao*wf%n_mo, wf%orbital_coefficients, 1, canonical_orbitals, 1)
!
!     Construct partitioning orbital basis, and determine active spaces
!
      call wf%orbital_partitioning()
!
      call wf%determine_n_x2_amplitudes()
      call wf%determine_n_gs_amplitudes()
      call wf%determine_n_es_amplitudes()
!
      wf%integrals = mo_integral_tool(wf%n_o, wf%n_v, wf%system%n_J, wf%need_g_abcd)
!
      call wf%integrals%initialize_storage()
!
      call wf%construct_and_save_mo_cholesky(wf%n_mo, wf%orbital_coefficients)
!
!     Frozen fock terms transformed from the canonical MO basis to 
!     the basis of orbital partitioning
!
      if (wf%exists_frozen_fock_terms) &
         call wf%update_MO_fock_contributions(canonical_orbitals)
!
      call mem%dealloc(canonical_orbitals, wf%n_ao, wf%n_mo)
!
      call wf%initialize_t1()
      call zero_array(wf%t1, wf%n_t1)
!
      call wf%integrals%update_t1_integrals(wf%t1)
!
      call wf%construct_fock()
      call wf%destruct_t1()
!
!     Keep partitioning orbital basis (for transformation of frozen MO fock terms)
!
      call mem%alloc(partitioning_orbitals, wf%n_ao, wf%n_mo)
      call dcopy(wf%n_ao*wf%n_mo, wf%orbital_coefficients, 1, partitioning_orbitals, 1)
!
!     Construct MLCC orbital basis
!
      call wf%construct_block_diagonal_fock_orbitals()
!
      call wf%construct_and_save_mo_cholesky(wf%n_mo, wf%orbital_coefficients)
!
!     Frozen fock terms transformed from the basis of orbital partitioning to 
!     the MLCC basis
!
      if (wf%exists_frozen_fock_terms) &
         call wf%update_MO_fock_contributions(partitioning_orbitals)
!
      call mem%dealloc(partitioning_orbitals, wf%n_ao, wf%n_mo)
!
      call wf%construct_and_save_mo_cholesky(wf%n_mo, wf%orbital_coefficients)
!
      call wf%check_orbital_space()
      call wf%print_orbital_space()
!
   end subroutine mo_preparations_mlcc2
!
!
   subroutine cleanup_mlcc2(wf)
!!
!!    Cleanup
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad and
!!    Alexander C. Paul , 2018
!!
      implicit none
!
      class(mlcc2) :: wf
!
      call wf%destruct_amplitudes()
      call wf%destruct_multipliers()

      call wf%destruct_orbital_coefficients()
      call wf%destruct_orbital_energies()
!
      call wf%destruct_right_excitation_energies()
      call wf%destruct_left_excitation_energies()
!
      call wf%destruct_core_MOs()
!
      call wf%destruct_fock()
      call wf%destruct_mo_fock_frozen()
!
      call wf%integrals%cleanup()
!
      call wf%destruct_nto_states()
      call wf%destruct_cnto_states()
!
      if (allocated(wf%l_files)) call wf%l_files%finalize_storer()
      if (allocated(wf%r_files)) call wf%r_files%finalize_storer()
!
      call output%printf('v', '- Cleaning up (a0) wavefunction', &
                         chars=[trim(convert_to_uppercase(wf%name_))], fs='(/t3, a)')
!
   end subroutine cleanup_mlcc2
!
   subroutine determine_n_x2_amplitudes_mlcc2(wf)
!!
!!    Determine number of x2 amplitudes
!!    Written by Sarai D. Folkestad, Apr 2019
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      wf%n_x2 = (wf%n_cc2_v)*(wf%n_cc2_o)*((wf%n_cc2_v)*(wf%n_cc2_o) + 1)/2
!
   end subroutine determine_n_x2_amplitudes_mlcc2
!
!
   subroutine determine_n_es_amplitudes_mlcc2(wf)
!!
!!    Determine number of es amplitudes
!!    Written by Sarai D. Folkestad, Apr 2019
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      wf%n_es_amplitudes = wf%n_t1 + wf%n_x2
!
   end subroutine determine_n_es_amplitudes_mlcc2
!
!
   subroutine determine_n_gs_amplitudes_mlcc2(wf)
!!
!!    Determine number of gs amplitudes
!!    Written by Sarai D. Folkestad, Apr 2019
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      wf%n_gs_amplitudes = wf%n_t1
!
   end subroutine determine_n_gs_amplitudes_mlcc2
!
!
   subroutine is_restart_safe_mlcc2(wf, task)
!!
!!    Is restart safe?
!!    Written by Eirik F. Kjønstad, Mar 2019 
!!
!!    'task' : Which type of restart we are considering.
!!             can be either 'ground state' or 'excited state'
!!
!!    Modified by Sarai D. Folkestad, Nov 2019
!!    
!!    Modified for MLCC
!!
      implicit none 
!
      class(mlcc2) :: wf 
!
      character(len=*), intent(in) :: task 
!
      integer :: n_o, n_v, n_gs_amplitudes, n_es_amplitudes
      integer :: n_cc2_o, n_cc2_v, n_ccs_o, n_ccs_v
!
      character(len=200) :: cc2_orbital_type
!
      call wf%restart_file%open_('read', 'rewind')
!
      call wf%restart_file%read_(n_o)
      call wf%restart_file%read_(n_v)
      call wf%restart_file%read_(n_gs_amplitudes)
      call wf%restart_file%read_(n_es_amplitudes)
!
      if (n_o .ne. wf%n_o) &
        call output%error_msg('attempted to restart from inconsistent number of occupied orbitals.')
!
      if (n_v .ne. wf%n_v) &
         call output%error_msg('attempted to restart from inconsistent number of virtual orbitals.')
!
      call wf%restart_file%read_(cc2_orbital_type) 
!
      if (cc2_orbital_type .ne. wf%cc2_orbital_type) &
         call output%error_msg('attempted MLCC restart ' // &
         'with inconsistent orbital type.')
!
      call wf%restart_file%read_(n_ccs_o)
      call wf%restart_file%read_(n_ccs_v)
      call wf%restart_file%read_(n_cc2_o)
      call wf%restart_file%read_(n_cc2_v)
!
      if (n_ccs_o .ne. wf%n_ccs_o) &
         call output%error_msg('attempted to restart from inconsistent ' // &
                                                   'number of ccs occupied orbitals.')
!
      if (n_ccs_v .ne. wf%n_ccs_v) &
         call output%error_msg('attempted to restart from inconsistent ' // &
                                                   'number of ccs virtual orbitals.')
!
      if (n_cc2_o .ne. wf%n_cc2_o) &
         call output%error_msg('attempted to restart from inconsistent ' // &
                                                   'number of cc2 occupied orbitals.')
!
      if (n_cc2_v .ne. wf%n_cc2_v) &
         call output%error_msg('attempted to restart from inconsistent ' // &
                                                   'number of cc2 virtual orbitals.')
!
      call wf%restart_file%close_()
!
      if (trim(task) == 'ground state') then 
!
         if (n_gs_amplitudes .ne. wf%n_gs_amplitudes) &
            call output%error_msg('attempted to restart from inconsistent number ' // &
                                    'of ground state amplitudes.')    
!
      elseif (trim(task) == 'excited state') then    
!
         if (n_es_amplitudes .ne. wf%n_es_amplitudes) &
            call output%error_msg('attempted to restart from inconsistent number ' // &
                                    'of excited state amplitudes.')     
!
      else
!
         call output%error_msg('attempted to restart, but the task was not recognized: ' // task)
!
      endif   
!
   end subroutine is_restart_safe_mlcc2
!
!
   subroutine write_cc_restart_mlcc2(wf)
!!
!!    Write CC restart file
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
!!    Modified by Sarai D. Folkestad, Nov 2019
!!    
!!    Modified for MLCC
!!
      implicit none
!
      class(mlcc2) :: wf 
!
!     Write information to restart file 
!
      call wf%restart_file%open_('write', 'rewind')
!
      call wf%restart_file%write_(wf%n_o)
      call wf%restart_file%write_(wf%n_v)
      call wf%restart_file%write_(wf%n_gs_amplitudes)
      call wf%restart_file%write_(wf%n_es_amplitudes)
      call wf%restart_file%write_(wf%cc2_orbital_type)
      call wf%restart_file%write_(wf%n_ccs_o)
      call wf%restart_file%write_(wf%n_ccs_v)
      call wf%restart_file%write_(wf%n_cc2_o)
      call wf%restart_file%write_(wf%n_cc2_v)
!
      call wf%restart_file%close_()
!
   end subroutine write_cc_restart_mlcc2
!
!
   subroutine contruct_mo_basis_transformation_mlcc2(wf, C1, C2, T)
!!
!!    Construct MO basis transformation 
!!    Written by Sarai D. Folekstad, Nov 2019
!!
!!    Constructs a transformation matrix 'T' which
!!    takes a matrix from one molecular orbital basis
!!    to another. 
!!
!!    'C1' : coefficients of the MO basis we end up in
!!
!!    'C2' : coefficients of the MO basis we start out with 
!!
!!    The transformation matrix is defined as
!!
!!       T = C1^T S C2
!!    
!!    where S is the AO overlap matrix.
!!
      implicit none
!
      class(mlcc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_mo), intent(in) :: C1, C2
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(out) :: T
!
      real(dp), dimension(:,:), allocatable :: S, X
!
      call mem%alloc(S, wf%n_ao, wf%n_ao)
      call wf%get_ao_s_wx(S)
!
!     X = C1^T S
!
      call mem%alloc(X, wf%n_mo, wf%n_ao)
!
      call dgemm('T', 'N', &
                  wf%n_mo, &
                  wf%n_ao, &
                  wf%n_ao, &
                  one,     &
                  C1,      &
                  wf%n_ao, &
                  S,       &
                  wf%n_ao, &
                  zero,    &
                  X,       &
                  wf%n_mo)
!
      call mem%dealloc(S, wf%n_ao, wf%n_ao)
!
!     T = X C2
!
      call dgemm('N', 'N', &
                  wf%n_mo, &
                  wf%n_mo, &
                  wf%n_ao, &
                  one,     &
                  X,       &
                  wf%n_mo, &
                  C2,      &
                  wf%n_ao, &
                  zero,    &
                  T,       &
                  wf%n_mo)
!
      call mem%dealloc(X, wf%n_mo, wf%n_ao)
!
   end subroutine contruct_mo_basis_transformation_mlcc2
!
!
   subroutine update_MO_fock_contributions_mlcc2(wf, C_old)
!!
!!    Updates MO Fock contributions 
!!    Written by Sarai D. Folkestad, Nov 2019
!!
!!    Updates the frozen contributions to the fock matrix
!!    from the old MO basis (C_old) to the current
!!    (wf%orbital_coefficients)
!!
!
      use array_utilities, only : symmetric_sandwich_right_transposition_replace
!
      implicit none
!
      class(mlcc2) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_mo), intent(in) :: C_old
!
      real(dp), dimension(:,:), allocatable :: T
!
      call mem%alloc(T, wf%n_mo, wf%n_mo)
!
      call wf%contruct_mo_basis_transformation(wf%orbital_coefficients, C_old, T)
!
      call symmetric_sandwich_right_transposition_replace(wf%mo_fock_frozen, T, wf%n_mo)
!
      call mem%dealloc(T, wf%n_mo, wf%n_mo)
!
   end subroutine update_MO_fock_contributions_mlcc2
!
!
end module mlcc2_class
! 
