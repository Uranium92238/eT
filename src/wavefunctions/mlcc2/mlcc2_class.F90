!construct_c1_integrals
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2021 the authors of eT
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
   use stream_file_class, only: stream_file
!
   implicit none
!
   type, extends(ccs) :: mlcc2
!
!     Orbital indices for levels
!
      integer :: n_ccs_o
      integer :: n_ccs_v
!
      integer :: n_cc2_o
      integer :: n_cc2_v
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
      logical :: nto_restart
!
      logical :: restart_orbitals
!
      type(sequential_file) :: jacobian_a1_intermediate_vv
      type(sequential_file) :: jacobian_a1_intermediate_oo
!
      type(stream_file) :: orbital_coefficients_mlcc_file
      type(stream_file) :: orbital_energies_mlcc_file
!
      type(stream_file) :: T_nto_o_file
      type(stream_file) :: T_cnto_o_file
      type(stream_file) :: T_cnto_v_file
!
   contains
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
      procedure :: get_orbital_differences                           => get_orbital_differences_mlcc2
!
      procedure :: determine_n_x2_amplitudes                         => determine_n_x2_amplitudes_mlcc2
      procedure :: determine_n_es_amplitudes                         => determine_n_es_amplitudes_mlcc2
      procedure :: determine_n_gs_amplitudes                         => determine_n_gs_amplitudes_mlcc2
!
      procedure :: get_cvs_projector                                 => get_cvs_projector_mlcc2
      procedure :: set_cvs_start_indices                             => set_cvs_start_indices_mlcc2
!
!     Read input
!
      procedure :: read_mlcc_settings                                => read_mlcc_settings_mlcc2
      procedure :: read_orbital_settings                             => read_orbital_settings_mlcc2
!
!     Orbital routines
!
      procedure :: mo_preparations                                   => mo_preparations_mlcc2
!
      procedure :: general_mlcc_mo_preparations &
                  => general_mlcc_mo_preparations_mlcc2
!
      procedure :: mo_preparations_from_restart &
                  => mo_preparations_from_restart_mlcc2
!
      procedure :: print_orbital_space                               => print_orbital_space_mlcc2
      procedure :: check_orbital_space                               => check_orbital_space_mlcc2
      procedure :: check_orthonormality_of_MOs                       => check_orthonormality_of_MOs_mlcc2
!
      procedure :: read_cnto_transformation_matrices &
                     => read_cnto_transformation_matrices_mlcc2
      procedure :: write_cnto_transformation_matrices &
                     => write_cnto_transformation_matrices_mlcc2
!
!
      procedure :: read_nto_transformation_matrix &
                     => read_nto_transformation_matrix_mlcc2
      procedure :: write_nto_transformation_matrix &
                     => write_nto_transformation_matrix_mlcc2
!
      procedure :: construct_ccs_cnto_transformation_matrices        => construct_ccs_cnto_transformation_matrices_mlcc2
!
      procedure :: construct_cholesky_orbitals                       => construct_cholesky_orbitals_mlcc2
!
      procedure :: construct_block_diagonal_fock_orbitals            => construct_block_diagonal_fock_orbitals_mlcc2
!
      procedure :: construct_cntos                                   => construct_cntos_mlcc2
      procedure :: construct_M_and_N_cnto                            => construct_M_and_N_cnto_mlcc2
      procedure :: construct_M_and_N_singles_cnto                    => construct_M_and_N_singles_cnto_mlcc2
      procedure :: ccs_calculation_for_cntos                         => ccs_calculation_for_cntos_mlcc2
      procedure :: add_doubles_M_and_N_cnto                          => add_doubles_M_and_N_cnto_mlcc2
!
      procedure :: construct_M_nto                                   => construct_M_nto_mlcc2
      procedure :: construct_ccs_nto_transformation_matrix           => construct_ccs_nto_transformation_matrix_mlcc2
      procedure :: construct_mixed_nto_canonical_orbitals            => construct_mixed_nto_canonical_orbitals_mlcc2
      procedure :: construct_paos                                    => construct_paos_mlcc2
!
      procedure :: orbital_partitioning                              => orbital_partitioning_mlcc2
      procedure :: diagonalize_M_and_N                               => diagonalize_M_and_N_mlcc2
!
      procedure :: update_MO_fock_contributions                      => update_MO_fock_contributions_mlcc2
!
      procedure :: construct_semicanonical_mlcc_orbitals &
                => construct_semicanonical_mlcc_orbitals_mlcc2
!
!     Ground state routines
!
      procedure :: construct_fock                                    => construct_fock_mlcc2
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
!     Restart
!
      procedure :: read_doubles_vector                               => read_doubles_vector_mlcc2
      procedure :: save_doubles_vector                               => save_doubles_vector_mlcc2
      procedure :: read_excitation_vector_file                       => read_excitation_vector_file_mlcc2
      procedure :: save_excitation_vector_on_file                    => save_excitation_vector_on_file_mlcc2
      procedure :: get_restart_vector                                => get_restart_vector_mlcc2
!
!     Summary
!
      procedure :: print_X1_diagnostics                              => print_X1_diagnostics_mlcc2
!
!     Initialize wavefunction
!
      procedure :: initialize                                        => initialize_mlcc2
!
!     File handling
!
      procedure :: save_mlcc_orbitals                                => save_mlcc_orbitals_mlcc2
      procedure :: read_mlcc_orbitals                                => read_mlcc_orbitals_mlcc2
!
   end type mlcc2
!
   interface
!
      include "./orbitals_mlcc2_interface.F90"
      include "./file_handling_mlcc2_interface.F90"
      include "./omega_mlcc2_interface.F90"
      include "./jacobian_mlcc2_interface.F90"
      include "./jacobian_transpose_mlcc2_interface.F90"
      include "./initialize_destruct_mlcc2_interface.F90"
      include "./fock_mlcc2_interface.F90"
!
   end interface 
!
contains
!
!
   subroutine initialize_mlcc2(wf, template_wf)
!!
!!    Initialize
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Adapted by Sarai D. Folkestad from CCS constructer, 2019
!!
      use citation_class,           only : citation
      use citation_printer_class,   only : eT_citations
!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      class(wavefunction), intent(in) :: template_wf
!
      type(citation), allocatable :: reference 
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
      wf%restart_orbitals = .false.
      wf%cnto_restart = .false.
      wf%nto_restart = .false.
!
      wf%cholesky_orbital_threshold = 1.0d-2
!
      wf%cc2_orbital_type = 'none'
!
      call wf%general_cc_preparations()
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
      wf%orbital_coefficients_mlcc_file = stream_file('orbitals_mlcc2')
      wf%orbital_energies_mlcc_file = stream_file('orbital_energies_mlcc2')
!
      if (trim(wf%cc2_orbital_type) == 'cnto' .or. trim(wf%cc2_orbital_type) == 'cnto-approx') then
!
         wf%T_cnto_o_file = stream_file('cnto_M_transformation')
         wf%T_cnto_v_file = stream_file('cnto_N_transformation')
!
      elseif (trim(wf%cc2_orbital_type) == 'nto-canonical') then
!
         wf%T_nto_o_file = stream_file('nto_M_transformation')
!
      endif
!
      reference = citation(implementation = 'MLCC2 and MLCCSD',                             &
                           journal        = 'J. Chem. Theory Comput.',                      &
                           title_         = 'Multilevel CC2 and CCSD Methods with &
                                             &Correlated Natural Transition Orbitals',      &
                           volume         = '16',                                           &
                           issue          = '1',                                            &
                           pages          = '179–189',                                      &
                           year           = '2019',                                         &
                           doi            = '10.1021/acs.jctc.9b00701',                     &
                           authors        = [character(len=25) :: 'Sarai Dery Folkestad',   &
                                                                  'Henrik Koch'])
!
      call eT_citations%add(reference)
!
   end subroutine initialize_mlcc2
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
      if (.not. input%is_section_present('mlcc')) &
         call output%error_msg('cannot do mlcc calculation without mlcc section in eT.inp')
!
      wf%restart_orbitals = input%is_keyword_present('orbital restart', 'mlcc') .or. &
                            input%is_keyword_present('restart', 'do')
!
      call input%get_required_keyword('cc2 orbitals', 'mlcc', wf%cc2_orbital_type)
!
!     Read orbital settings
!
      call wf%read_orbital_settings(wf%cc2_orbital_type, 'cc2', wf%n_cc2_o, wf%n_cc2_v)
!
   end subroutine read_mlcc_settings_mlcc2
!
!
   subroutine read_orbital_settings_mlcc2(wf, orbital_type, level_string, n_level_o, n_level_v)
!!
!!    Read orbital settings
!!    Written by Sarai D. Folkestad, Apr 2019
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
!!    If orbitals are CNTOs or NTOs/Canonical:
!!
!!       - reads the number of occupied cntos and ntos (requested keyword for CNTO/NTO)
!!       - reads the number of virtual cntos and canonical
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      character(len=*), intent(in) :: orbital_type
      character(len=*), intent(in) :: level_string
!
      integer, intent(out) :: n_level_o, n_level_v
!
      if (trim(orbital_type) == 'cnto' .or. trim(orbital_type) == 'cnto-approx') then
!
         wf%cnto_restart = .false.
!
         if (input%is_keyword_present('cnto restart', 'mlcc') .or. &
            input%is_keyword_present('restart', 'do')) wf%cnto_restart = .true.
!
         wf%n_cnto_states = input%get_n_elements_for_keyword('cnto states', 'mlcc')
!
         if (wf%n_cnto_states == 0) &
               call output%error_msg('to construct CNTOs excitation vectors must be specified.')
!
         call wf%initialize_cnto_states()
!
         call input%get_array_for_keyword('cnto states', 'mlcc', &
               wf%n_cnto_states, wf%cnto_states)
!
!
         call input%get_required_keyword('cnto occupied ' // trim(level_string), &
               'mlcc', n_level_o)
!
         n_level_v = n_level_o*(wf%n_v/wf%n_o)
!
         call input%get_keyword('cnto virtual ' // trim(level_string), 'mlcc', n_level_v)
!
!
      elseif (trim(orbital_type) == 'cholesky') then
!
         call input%get_keyword('cholesky threshold', 'mlcc', &
                  wf%cholesky_orbital_threshold)
!
      elseif (trim(orbital_type) == 'nto-canonical') then
!
         wf%nto_restart = .false.
!
         if (input%is_keyword_present('nto restart', 'mlcc') .or. &
            input%is_keyword_present('restart', 'do')) wf%nto_restart = .true.
!
         wf%n_nto_states = input%get_n_elements_for_keyword('nto states', 'mlcc')
!
         if (wf%n_nto_states == 0) &
               call output%error_msg('to construct NTOs excitation vectors must be specified.')
!
         call wf%initialize_nto_states()
         call input%get_array_for_keyword('nto states', 'mlcc', &
               wf%n_nto_states, wf%nto_states)
!
         call input%get_required_keyword('nto occupied ' // trim(level_string), &
               'mlcc', n_level_o)
!
         n_level_v = n_level_o*(wf%n_v/wf%n_o)
!
         call input%get_keyword('canonical virtual ' // trim(level_string), &
               'mlcc', n_level_v)
!
      elseif (trim(orbital_type) == 'cholesky-pao') then
!
         call input%get_keyword('cholesky threshold', 'mlcc', &
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
      call wf%eri%get_eri_t1('ovov', g_iajb, 1, wf%n_o, 1, wf%n_v, 1, wf%n_o, 1, wf%n_v)
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
      call wf%eri%get_eri_t1('ovov', g_iajb,                &
                             1, wf%n_cc2_o, &
                             1, wf%n_cc2_v, &
                             1, wf%n_cc2_o, &
                             1, wf%n_cc2_v)
!
      call wf%eri%get_eri_t1('vovo', g_aibj,                &
                             1, wf%n_cc2_v, &
                             1, wf%n_cc2_o, &
                             1, wf%n_cc2_v, &
                             1, wf%n_cc2_o)
!
!$omp parallel do private(a,i,b,j) reduction(+:correlation_energy)
      do b = 1, wf%n_cc2_v
         do i = 1, wf%n_cc2_o 
            do j = 1, wf%n_cc2_o 
               do a = 1, wf%n_cc2_v
!
                  correlation_energy = correlation_energy - (g_aibj(a,i,b,j)/&
                                    (wf%orbital_energies(wf%n_o + a)  &
                                   + wf%orbital_energies(wf%n_o + b)  &
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
      call wf%check_orthonormality_of_MOs()
!
   end subroutine check_orbital_space_mlcc2
!
!
   subroutine get_orbital_differences_mlcc2(wf, orbital_differences, N)
!!
!!    Get orbital differences 
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 2019
!!
      implicit none
!
      class(mlcc2), intent(in) :: wf
!
      integer, intent(in) :: N 
      real(dp), dimension(N), intent(out) :: orbital_differences
!
      integer :: a, i, ai, b, j, bj, aibj
!
      call wf%ccs%get_orbital_differences(orbital_differences, wf%n_t1)
!
      if (N .eq. wf%n_t1) return ! Requested only singles orbital differences
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
                                 wf%orbital_energies(a + wf%n_o) &
                                 - wf%orbital_energies(i) &
                                 + wf%orbital_energies(b + wf%n_o) &
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
   end subroutine get_orbital_differences_mlcc2
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
      class(mlcc2), intent(inout) :: wf 
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
                                       1, 1, wf%n_cc2_o, wf%n_cc2_v)
!
      call symmetric_sum(t2bar, wf%n_cc2_o*wf%n_cc2_v)
!
      call mem%alloc(g_iajb, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v)
!      
      call wf%eri%get_eri_t1('ovov', g_iajb,                &
                             1, wf%n_cc2_o, &
                             1, wf%n_cc2_v, &
                             1, wf%n_cc2_o, &
                             1, wf%n_cc2_v)
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
                     (-  wf%orbital_energies(a + wf%n_o) &
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
      call zero_array(equation, wf%n_gs_amplitudes)
!  
      call wf%jacobian_transpose_ccs_a1(equation, wf%t1bar)
      call wf%jacobian_transpose_ccs_b1(equation, wf%t1bar)
      call wf%jacobian_transpose_cc2_a1(equation, wf%t1bar, wf%n_cc2_o, wf%n_cc2_v, &
                                 1, 1)
!
!
!     equation += sum_bjck tbar_bjck A_{bjck,ai}
!
      call wf%jacobian_transpose_cc2_b1(equation, t2bar, wf%n_cc2_o, wf%n_cc2_v, &
                           1, 1, wf%n_cc2_o, wf%n_cc2_v)
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
      call wf%eri%get_eri_t1('vovo', g_aibj,1, wf%n_cc2_v, &
                                            1, wf%n_cc2_o, &
                                            1, wf%n_cc2_v, &
                                            1, wf%n_cc2_o)
!
!$omp parallel do private(a, i, b, j)
      do b = 1, wf%n_cc2_v 
         do j = 1, wf%n_cc2_o 
            do i = 1, wf%n_cc2_o
               do a = 1, wf%n_cc2_v
!
                  s_aibj(a, i, b, j) = (g_aibj(a, i, b, j))/ &
                           (-  wf%orbital_energies(a + wf%n_o) &
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
                                       1, 1, wf%n_cc2_o, wf%n_cc2_v)
!
      call symmetric_sum(t2bar, wf%n_cc2_o*wf%n_cc2_v)
!
      call mem%alloc(g_iajb, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v)
!
      call wf%eri%get_eri_t1('ovov', g_iajb,                &
                             1, wf%n_cc2_o, &
                             1, wf%n_cc2_v, &
                             1, wf%n_cc2_o, &
                             1, wf%n_cc2_v)
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
                           (-  wf%orbital_energies(a + wf%n_o) &
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
         if ((i .lt. 1) .or. (i  .gt. wf%n_cc2_o)) then 
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
!     Sanity check on core MOs -> Are they in active space?
!
      do mo = 1, wf%n_core_MOs
!
         if ((wf%core_MOs(mo) .lt. 1) &
         .or. (wf%core_MOs(mo) .gt. wf%n_cc2_o)) &
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
      implicit none
!
      class(mlcc2) :: wf
!
      logical :: has_restart_files
!
      has_restart_files = wf%orbital_coefficients_mlcc_file%exists() .and. &
                          wf%orbital_energies_mlcc_file%exists()
      if (.not. has_restart_files) wf%restart_orbitals = .false.
!
      has_restart_files = wf%T_cnto_o_file%exists() .and. wf%T_cnto_v_file%exists()
      if (.not. has_restart_files) wf%cnto_restart = .false.
!
      has_restart_files = wf%T_nto_o_file%exists()
      if (.not. has_restart_files) wf%nto_restart = .false.
!
      if (wf%restart_orbitals) then
!
         call output%printf('m', 'Requested orbital restart, &
                           &reading orbitals and orbital energies')
!
         call wf%mo_preparations_from_restart()
!
      else
!
         call wf%general_mlcc_mo_preparations()
!
      endif
!
      call wf%print_orbital_space()
      call wf%check_orbital_space()
!
   end subroutine mo_preparations_mlcc2
!
!
   subroutine general_mlcc_mo_preparations_mlcc2(wf)
!!
!!    General MLCC MO prepatations
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
      implicit none
!
      class(mlcc2) :: wf
!
      real(dp), dimension(:,:), allocatable :: canonical_orbitals
      real(dp), dimension(:,:), allocatable :: partitioning_orbitals
      real(dp), dimension(:,:), allocatable :: T
!
!      Keep canonical orbitals (for transformation of frozen MO fock terms)
!
      call mem%alloc(canonical_orbitals, wf%ao%n, wf%n_mo)
      call dcopy(wf%ao%n*wf%n_mo, wf%orbital_coefficients, 1, canonical_orbitals, 1)
!
!     Construct partitioning orbital basis, and determine active spaces
!
      call wf%orbital_partitioning()
!
      call mem%alloc(T, wf%n_mo, wf%n_mo)
      call wf%contruct_mo_basis_transformation(wf%orbital_coefficients, canonical_orbitals, T)
      call wf%eri%update_cholesky_mo(T)
!
      call wf%determine_n_x2_amplitudes()
      call wf%determine_n_gs_amplitudes()
      call wf%determine_n_es_amplitudes()
!
!     Frozen fock terms transformed from the canonical MO basis to 
!     the basis of orbital partitioning
!
      if (wf%exists_frozen_fock_terms) &
         call wf%update_MO_fock_contributions(canonical_orbitals)
!
      call mem%dealloc(canonical_orbitals, wf%ao%n, wf%n_mo)
!
      call wf%initialize_t1()
      call zero_array(wf%t1, wf%n_t1)
      call wf%eri%set_t1_to_mo()
!
      call wf%construct_fock(task = 'block diagonalize fock')
      call wf%destruct_t1()
!
!     Keep partitioning orbital basis (for transformation of frozen MO fock terms)
!
      call mem%alloc(partitioning_orbitals, wf%ao%n, wf%n_mo)
      call dcopy(wf%ao%n*wf%n_mo, wf%orbital_coefficients, 1, partitioning_orbitals, 1)
!
!     Construct MLCC orbital basis
!
      call wf%construct_semicanonical_mlcc_orbitals()
!
      call wf%contruct_mo_basis_transformation(wf%orbital_coefficients, partitioning_orbitals, T)
      call wf%eri%update_cholesky_mo(T)
      call mem%dealloc(T, wf%n_mo, wf%n_mo)
!
!     Frozen fock terms transformed from the basis of orbital partitioning to 
!     the MLCC basis
!
      if (wf%exists_frozen_fock_terms) &
         call wf%update_MO_fock_contributions(partitioning_orbitals)
!
      call mem%dealloc(partitioning_orbitals, wf%ao%n, wf%n_mo)
!
!     Print MLCC orbital coefficients to file
!
      call wf%save_mlcc_orbitals()
!
   end subroutine general_mlcc_mo_preparations_mlcc2
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
      call output%printf('v', '- Cleaning up (a0) wavefunction', &
                         chars=[trim(convert_to_uppercase(wf%name_))], fs='(/t3, a)')
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
      call wf%destruct_frozen_CCT()
!
      call wf%eri%cleanup()
!
      call wf%destruct_nto_states()
      call wf%destruct_cnto_states()
!
      deallocate(wf%ao)
      if (wf%embedded) deallocate(wf%embedding)
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
      real(dp), dimension(wf%ao%n, wf%n_mo), intent(in) :: C_old
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
   subroutine print_X1_diagnostics_mlcc2(wf, X, label)
!!
!!    Print X1 diagnostics
!!    Written by Sarai D. Folkestad, Nov 2019       
!!
      implicit none
!
      class(mlcc2), intent(in) :: wf
!     
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: X
!
      character(len=1), intent(in) :: label
!
      real(dp), dimension(:), allocatable :: X_internal
!
      real(dp) :: internal_fraction
!
      integer :: a, i, ai, ai_full
!
      call wf%ccs%print_X1_diagnostics(X, label)
!
      call mem%alloc(X_internal, wf%n_cc2_v*wf%n_cc2_o)
!
      do i = 1, wf%n_cc2_o
         do a = 1, wf%n_cc2_v
!
            ai = wf%n_cc2_v*(i-1) + a 
            ai_full = wf%n_v*(i-1) + a
!
            X_internal(ai) = X(ai_full)
!
         enddo
      enddo
!
      internal_fraction = get_l2_norm(X_internal, wf%n_cc2_v*wf%n_cc2_o)&
                           /get_l2_norm(X,wf%n_es_amplitudes)
!
      call output%printf('n', 'MLCC diagnostics:', fs='(/t6,a)')
!
      call output%printf('n', '|(a0)1^internal|/|(a0)| =  (f19.12)', &
            reals=[internal_fraction], chars=[label, label], fs='(/t6,a)')
!
      internal_fraction = get_l2_norm(X_internal, wf%n_cc2_v*wf%n_cc2_o)&
                           /get_l2_norm(X(1:wf%n_t1), wf%n_t1)
!
      call output%printf('n', '|(a0)1^internal|/|(a0)1| = (f19.12)', &
            reals=[internal_fraction], chars=[label, label], fs='(t6,a)')
!
      call mem%dealloc(X_internal, wf%n_cc2_v*wf%n_cc2_o)
!
   end subroutine print_X1_diagnostics_mlcc2
!
!
   subroutine mo_preparations_from_restart_mlcc2(wf)
!!
!!    General MO prepatations from restart
!!    Written by Sarai D. Folkestad, May 2020
!!
!!    Reads MLCC orbitals and partitionings from file
!!    and transforms frozen Fock matrices and Cholesky vectors to 
!!    the MLCC basis
!!
      implicit none
!
      class(mlcc2) :: wf
!
      real(dp), dimension(:,:), allocatable :: canonical_orbitals
      real(dp), dimension(:,:), allocatable :: T
!
      call output%warning_msg('Number of orbitals in active and inactive spaces read &
                              &from orbital_coefficients_mlcc. Make sure restarted &
                              &calculation and restart files are consistent!')
!
!     Keep canonical orbitals (for transformation of frozen MO fock terms)
!
      call mem%alloc(canonical_orbitals, wf%ao%n, wf%n_mo)
      call dcopy(wf%ao%n*wf%n_mo, wf%orbital_coefficients, 1, canonical_orbitals, 1)
!
!     Read partitionings from restart file
!
      call wf%read_mlcc_orbitals()
!
      call mem%alloc(T, wf%n_mo, wf%n_mo)
      call wf%contruct_mo_basis_transformation(wf%orbital_coefficients, canonical_orbitals, T)
      call wf%eri%update_cholesky_mo(T)
      call mem%dealloc(T, wf%n_mo, wf%n_mo)
!
      call wf%determine_n_x2_amplitudes()
      call wf%determine_n_gs_amplitudes()
      call wf%determine_n_es_amplitudes()
!
!     Frozen fock terms transformed from the canonical MO basis to 
!     the basis of orbital partitioning
!
      if (wf%exists_frozen_fock_terms) &
         call wf%update_MO_fock_contributions(canonical_orbitals)
!
      call mem%dealloc(canonical_orbitals, wf%ao%n, wf%n_mo)
!
   end subroutine mo_preparations_from_restart_mlcc2
!
end module mlcc2_class
