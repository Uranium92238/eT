!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2022 the authors of eT
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
module mlccsd_class
!
!!
!!    Multilevel Coupled cluster singles and doubles (MLCCSD) class module
!!    Written by Sarai D. Folkestad, June 2017 and spring 2019
!!
!!    This MLCCSD wavefunction is given by
!!
!!       | MLCCSD > = exp(T_1 + S_2 + T_2) | HF >
!!
!!    where T_1 is the single excitation operator from standard CC,
!!    S_2 is a double excitation operator which only contains
!!    excitations within an active orbital space (the CC2 and CCSD orbitals)
!!    and T_2 is a double excitation operator which only contains
!!    excitations within another active orbital space (the CCSD orbitals)
!!
!!    The S_2 operator is determined to first order in the
!!    fluctuation potential U (H = F + U). T_2 is treated to infinite order
!!    in U.
!!
!!    This class handles the ground and excited states for
!!    MLCCSD.
!!
!!    For further references on MLCC see:
!!
!!       Myhre, R. H., & Koch, H., JCP, 145(4), 044111 (2016)
!!       Folkestad, S. D., & Koch, H.,  JCTC 16, 1, 179-189 (2020)
!!
!
   use mlcc2_class, only: mlcc2
!
   use parameters
   use global_out, only: output
   use timings_class, only: timings
   use memory_manager_class, only: mem
   use stream_file_class, only: stream_file
   use sequential_file_class, only: sequential_file
   use batching_index_class, only: batching_index
!
   implicit none
!
   type, extends(mlcc2) :: mlccsd
!
!     Requested levels
!
      logical :: do_ccs
      logical :: do_cc2
!
!     Orbital indices for levels
!
      integer :: n_ccsd_o
      integer :: n_ccsd_v
!
!     Orbital type string
!
      character(len=200) :: ccsd_orbital_type
!
      real(dp), dimension(:,:), allocatable  :: orbital_coefficients_cc2
      real(dp), dimension(:), allocatable    :: orbital_energies_cc2
!
      integer :: n_t2
!
      real(dp), dimension(:), allocatable :: t2
      real(dp), dimension(:,:), allocatable :: O_o, O_v ! Orthogonal transformation matrix
                                                        ! which transforms between MLCCSD basis
                                                        ! and basis for CC2 amplitude determination
!
!     Intermediates for Jacobian transformation
!
      type(sequential_file) :: jacobian_c2_intermediate_oovo_1
      type(sequential_file) :: jacobian_c2_intermediate_oovo_2
      type(sequential_file) :: jacobian_d2_intermediate
      type(sequential_file) :: jacobian_e2_intermediate
      type(sequential_file) :: jacobian_g2_intermediate_vovo
      type(sequential_file) :: jacobian_g2_intermediate_vv
      type(sequential_file) :: jacobian_g2_intermediate_oo
      type(sequential_file) :: jacobian_h2_intermediate_vovo_1
      type(sequential_file) :: jacobian_h2_intermediate_vovo_2
      type(sequential_file) :: jacobian_j2_intermediate_oooo
      type(sequential_file) :: jacobian_j2_intermediate_oovv
!
   contains
!
!     Read settings
!
      procedure :: read_mlcc_settings           => read_mlcc_settings_mlccsd
!
!     MO preparations
!
      procedure :: mo_preparations              => mo_preparations_mlccsd
      procedure :: determine_n_x2_amplitudes    => determine_n_x2_amplitudes_mlccsd
      procedure :: determine_n_gs_amplitudes    => determine_n_gs_amplitudes_mlccsd
      procedure :: check_orbital_space          => check_orbital_space_mlccsd
      procedure :: print_orbital_space          => print_orbital_space_mlccsd
!
      procedure :: construct_mlccsd_basis_transformation_matrix &
                  => construct_mlccsd_basis_transformation_matrix_mlccsd
!
!     Orbital partitioning
!
      procedure :: orbital_partitioning         => orbital_partitioning_mlccsd
!
      procedure :: construct_cholesky_orbitals  => construct_cholesky_orbitals_mlccsd
      procedure :: construct_paos               => construct_paos_mlccsd
!
      procedure :: construct_cc2_cnto_transformation_matrices &
                  => construct_cc2_cnto_transformation_matrices_mlccsd
!
      procedure :: cc2_calculation_for_cntos    => cc2_calculation_for_cntos_mlccsd
!
      procedure :: construct_semicanonical_mlcc_orbitals &
                => construct_semicanonical_mlcc_orbitals_mlccsd
!
      procedure :: construct_orbitals_cc2 => construct_orbitals_cc2_mlccsd
!
!     Omega
!
      procedure :: construct_fock               => construct_fock_mlccsd
      procedure :: construct_omega              => construct_omega_mlccsd
!
      procedure :: omega_ccsd_a2                => omega_ccsd_a2_mlccsd
      procedure :: omega_ccsd_b2                => omega_ccsd_b2_mlccsd
      procedure :: omega_ccsd_c2                => omega_ccsd_c2_mlccsd
      procedure :: omega_ccsd_d2                => omega_ccsd_d2_mlccsd
      procedure :: omega_ccsd_e2                => omega_ccsd_e2_mlccsd
      procedure :: omega_ccsd_f2                => omega_ccsd_f2_mlccsd
!
      procedure :: construct_x2                 => construct_x2_mlccsd
      procedure :: construct_u_aibj             => construct_u_aibj_mlccsd
      procedure :: get_orbital_differences      => get_orbital_differences_mlccsd
      procedure :: construct_cc2_amplitudes     => construct_cc2_amplitudes_mlccsd
!
      procedure :: calculate_energy             => calculate_energy_mlccsd
      procedure :: set_initial_amplitudes_guess => set_initial_amplitudes_guess_mlccsd
      procedure :: set_t2_to_cc2_guess          => set_t2_to_cc2_guess_mlccsd
!
      procedure :: get_amplitudes               => get_amplitudes_mlccsd
      procedure :: set_amplitudes               => set_amplitudes_mlccsd
!
      procedure :: form_newton_raphson_t_estimate &
                  => form_newton_raphson_t_estimate_mlccsd
!
!     Jacobian transformation
!
      procedure :: jacobian_transformation    => jacobian_transformation_mlccsd
!
      procedure :: jacobian_ccsd_d2_1  => jacobian_ccsd_d2_1_mlccsd
      procedure :: jacobian_ccsd_d2_2  => jacobian_ccsd_d2_2_mlccsd
      procedure :: jacobian_ccsd_d2_3  => jacobian_ccsd_d2_3_mlccsd
      procedure :: jacobian_ccsd_d2_4  => jacobian_ccsd_d2_4_mlccsd
      procedure :: jacobian_ccsd_d2_5  => jacobian_ccsd_d2_5_mlccsd
      procedure :: jacobian_ccsd_d2    => jacobian_ccsd_d2_mlccsd
!
      procedure :: jacobian_ccsd_c2_1  => jacobian_ccsd_c2_1_mlccsd
      procedure :: jacobian_ccsd_c2_2  => jacobian_ccsd_c2_2_mlccsd
      procedure :: jacobian_ccsd_c2_3  => jacobian_ccsd_c2_3_mlccsd
      procedure :: jacobian_ccsd_c2_4  => jacobian_ccsd_c2_4_mlccsd
      procedure :: jacobian_ccsd_c2    => jacobian_ccsd_c2_mlccsd
!
      procedure :: jacobian_ccsd_b2    => jacobian_ccsd_b2_mlccsd
      procedure :: jacobian_cc2_b2     => jacobian_cc2_b2_mlccsd
!
      procedure :: jacobian_ccsd_e2    => jacobian_ccsd_e2_mlccsd
      procedure :: jacobian_ccsd_f2    => jacobian_ccsd_f2_mlccsd
      procedure :: jacobian_ccsd_g2    => jacobian_ccsd_g2_mlccsd
      procedure :: jacobian_ccsd_h2    => jacobian_ccsd_h2_mlccsd
!
      procedure :: jacobian_ccsd_i2    => jacobian_ccsd_i2_mlccsd
      procedure :: jacobian_ccsd_i2_1  => jacobian_ccsd_i2_1_mlccsd
      procedure :: jacobian_ccsd_i2_2  => jacobian_ccsd_i2_2_mlccsd
!
      procedure :: jacobian_ccsd_k2    => jacobian_ccsd_k2_mlccsd
      procedure :: jacobian_ccsd_j2    => jacobian_ccsd_j2_mlccsd
!
      procedure :: prepare_for_jacobian => prepare_for_jacobian_mlccsd
!
      procedure :: save_jacobian_c2_intermediates  => save_jacobian_c2_intermediates_mlccsd
      procedure :: save_jacobian_d2_intermediate   => save_jacobian_d2_intermediate_mlccsd
      procedure :: save_jacobian_e2_intermediate   => save_jacobian_e2_intermediate_mlccsd
      procedure :: save_jacobian_g2_intermediates  => save_jacobian_g2_intermediates_mlccsd
      procedure :: save_jacobian_h2_intermediates  => save_jacobian_h2_intermediates_mlccsd
      procedure :: save_jacobian_j2_intermediates  => save_jacobian_j2_intermediates_mlccsd
!
!     Excited state
!
      procedure :: print_X1_diagnostics       =>  print_X1_diagnostics_mlccsd
!
!     CVS
!
      procedure :: get_cvs_projector         => get_cvs_projector_mlccsd
      procedure :: set_cvs_start_indices     => set_cvs_start_indices_mlccsd
!
!     Initialize/destruct
!
      procedure :: initialize_u_aibj            => initialize_u_aibj_mlccsd
      procedure :: destruct_u_aibj              => destruct_u_aibj_mlccsd
!
      procedure :: initialize_t2                => initialize_t2_mlccsd
      procedure :: destruct_t2                  => destruct_t2_mlccsd
!
      procedure :: initialize_amplitudes        => initialize_amplitudes_mlccsd
      procedure :: destruct_amplitudes          => destruct_amplitudes_mlccsd
!
      procedure :: initialize_orbital_coefficients_cc2 &
                  => initialize_orbital_coefficients_cc2_mlccsd
!
      procedure :: destruct_orbital_coefficients_cc2 &
                  => destruct_orbital_coefficients_cc2_mlccsd
!
      procedure :: initialize_orbital_energies_cc2 &
                  => initialize_orbital_energies_cc2_mlccsd
!
      procedure :: destruct_orbital_energies_cc2 &
                  => destruct_orbital_energies_cc2_mlccsd
!
!
      procedure :: initialize_O_o                 => initialize_O_o_mlccsd
      procedure :: destruct_O_o                   => destruct_O_o_mlccsd
      procedure :: initialize_O_v                 => initialize_O_v_mlccsd
      procedure :: destruct_O_v                   => destruct_O_v_mlccsd
!
!     Read/save
!
      procedure :: read_amplitudes                => read_amplitudes_mlccsd
      procedure :: save_amplitudes                => save_amplitudes_mlccsd
      procedure :: read_excitation_vector_file    => read_excitation_vector_file_mlccsd
      procedure :: save_excitation_vector_on_file => save_excitation_vector_on_file_mlccsd
      procedure :: get_restart_vector             => get_restart_vector_mlccsd
!
      procedure :: save_mlcc_orbitals             => save_mlcc_orbitals_mlccsd
      procedure :: read_mlcc_orbitals             => read_mlcc_orbitals_mlccsd
!
!     Cleanup
!
      procedure :: cleanup                        => cleanup_mlccsd
!
!     Initialize
!
      procedure :: initialize                     => initialize_mlccsd
!
!     Restart
!
      procedure :: mo_preparations_from_restart   => mo_preparations_from_restart_mlccsd
!
      procedure :: scale_amplitudes => scale_amplitudes_mlccsd
!
   end type mlccsd
!
   interface
!
      include "./orbitals_mlccsd_interface.F90"
      include "./omega_mlccsd_interface.F90"
      include "./initialize_destruct_mlccsd_interface.F90"
      include "./file_handling_mlccsd_interface.F90"
      include "./set_get_mlccsd_interface.F90"
      include "./jacobian_mlccsd_interface.F90"
      include "./fock_mlccsd_interface.F90"
!
   end interface
!
contains
!
!
   subroutine initialize_mlccsd(wf, template_wf)
!!
!!    Initialize mlccsd
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Adapted by Sarai D. Folkestad from CCS constructer, 2020
!!
      use citation_printer_class, only: eT_citations
      use wavefunction_class,     only: wavefunction
      use stream_file_class,      only: stream_file
!
      implicit none
!
      class(mlccsd), intent(inout) :: wf
!
      class(wavefunction), intent(in) :: template_wf
!
!     If we have a CC2 level, we will set the AO fock matrix from the template wavefunction.
!     The AO fock is currently only constructed for the reference wavefunctions.
!     Therefore we will stop if template_wf is not MLHF or RHF.
!
      if (trim(template_wf%name_) .ne. 'rhf' .and. trim(template_wf%name_) .ne. 'mlhf') &
         call output%error_msg('in initialization of the MLCCSD wavefunction.')
!
      wf%name_ = 'mlccsd'
!
      wf%n_ccsd_o = 0
      wf%n_ccsd_v = 0
      wf%n_cc2_o  = 0
      wf%n_cc2_v  = 0
      wf%n_ccs_o  = 0
      wf%n_ccs_v  = 0
!
      wf%do_ccs = .false.
      wf%do_cc2 = .false.
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
      wf%ccsd_orbital_type = 'none'
!
      call wf%general_cc_preparations()
      call wf%set_variables_from_template_wf(template_wf)
      call wf%print_banner()
!
      if (wf%ao%has_ghost_atoms()) &
         call output%warning_msg("Ghosts are experimental in multilevel.")
!
      if (wf%bath_orbital) call output%error_msg('Bath orbitals not yet implemented for MLCCSD')
!
      call wf%read_mlcc_settings()
!
      wf%n_t1 = (wf%n_o)*(wf%n_v)
!
      call wf%initialize_fock()
!
      if (wf%do_cc2) then
!
         call wf%initialize_orbital_energies_cc2()
         call wf%initialize_orbital_coefficients_cc2()
!
         call wf%initialize_ao_fock()
!
         call dcopy(wf%ao%n**2, template_wf%ao_fock, 1, wf%ao_fock, 1)
!
         call wf%initialize_mo_fock()
!
      endif
!
      wf%orbital_coefficients_mlcc_file = stream_file('orbitals_mlccsd')
      wf%orbital_energies_mlcc_file = stream_file('orbital_energies_mlccsd')
!
      if (trim(wf%ccsd_orbital_type) == 'cnto' .or. trim(wf%ccsd_orbital_type) == 'cnto-approx') then
!
         wf%T_cnto_o_file = stream_file('cnto_M_transformation')
         wf%T_cnto_v_file = stream_file('cnto_N_transformation')
!
      endif
!
      call eT_citations%add('MLCC2 and MLCCSD')
!
   end subroutine initialize_mlccsd
!
!
   subroutine read_mlcc_settings_mlccsd(wf)
!!
!!    Read MLCC settings
!!    Written by Sarai D. Folkestad, Apr 2019
!!
!!    Reads the mlcc sections of the
!!    input file.
!!
!!    Calls two routines that handle
!!    reading of orbital type etc.
!!
      use global_in, only: input
!
      implicit none
!
      class(mlccsd), intent(inout) :: wf
!
      logical :: do_ccsd
!
      if (.not. input%is_section_present('mlcc')) &
         call output%error_msg('cannot do mlcc calculation without mlcc section in eT.inp')
!
      wf%do_ccs  = input%is_string_in_cs_list('levels', 'mlcc', 'ccs')
      wf%do_cc2  = input%is_string_in_cs_list('levels', 'mlcc', 'cc2')
      do_ccsd = input%is_string_in_cs_list('levels', 'mlcc', 'ccsd')
!
      if (.not. do_ccsd) then
!
         call output%error_msg('CCSD level was not requested in the MLCCSD calculation.')
!
      endif
!
      wf%restart_orbitals = input%is_keyword_present('orbital restart', 'mlcc')
!
!     Get orbital types
!
      if (wf%do_cc2 .and. wf%do_ccs) &
            call input%get_required_keyword('cc2 orbitals', 'mlcc', wf%cc2_orbital_type)
!
      call input%get_required_keyword('ccsd orbitals', 'mlcc', wf%ccsd_orbital_type)
!
!     Get specific CC2 and CCSD orbital settings
!
      if (wf%do_cc2 .and. wf%do_ccs) &
         call wf%read_orbital_settings(wf%cc2_orbital_type, 'cc2', wf%n_cc2_o, wf%n_cc2_v)
!
      call wf%read_orbital_settings(wf%ccsd_orbital_type, 'ccsd', wf%n_ccsd_o, wf%n_ccsd_v)
!
   end subroutine read_mlcc_settings_mlccsd
!
!
   subroutine mo_preparations_mlccsd(wf)
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
!!    If there is a CC2 space, we must construct the transformation matrix which transforms
!!    between the MLCCSD basis and the basis where s amplitudes are constructed.
!!
      implicit none
!
      class(mlccsd) :: wf
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
      if (wf%do_cc2) call wf%construct_orbitals_cc2()
!
      call wf%print_orbital_space()
      call wf%check_orbital_space()
!
   end subroutine mo_preparations_mlccsd
!
!
   subroutine determine_n_x2_amplitudes_mlccsd(wf)
!!
!!    Determine number of x2 amplitudes
!!    Written by Sarai D. Folkestad, Apr 2019
!!
      implicit none
!
      class(mlccsd), intent(inout) :: wf
!
      integer :: n_v, n_o
!
      n_v = wf%n_cc2_v + wf%n_ccsd_v
      n_o = wf%n_cc2_o + wf%n_ccsd_o
!
      wf%n_x2 = n_v*n_o*(n_v*n_o + 1)/2
      wf%n_t2 = (wf%n_ccsd_v)*(wf%n_ccsd_o)*((wf%n_ccsd_v)*(wf%n_ccsd_o) + 1)/2
!
   end subroutine determine_n_x2_amplitudes_mlccsd
!
!
   subroutine determine_n_gs_amplitudes_mlccsd(wf)
!!
!!    Determine number of gs amplitudes
!!    Written by Sarai D. Folkestad, Apr 2019
!!
      implicit none
!
      class(mlccsd), intent(inout) :: wf
!
      wf%n_gs_amplitudes = wf%n_t1 + wf%n_t2
!
   end subroutine determine_n_gs_amplitudes_mlccsd
!
!
   subroutine check_orbital_space_mlccsd(wf)
!!
!!    Check orbital space
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Jan 2019
!!
      implicit none
!
      class(mlccsd) :: wf
!
      if (wf%n_ccsd_v == 0) &
         call output%error_msg('no virtual ccsd orbitals in mlccsd calulation.')
!
      if (wf%n_ccsd_o == 0) &
         call output%error_msg('no occupied ccsd orbitals in mlccsd calulation.')
!
      if (wf%n_ccs_o + wf%n_ccs_v + wf%n_cc2_o + wf%n_cc2_v  == 0)  &
         call output%warning_msg('no inactive orbitals in mlccsd calulation, '// &
                                 'recomended to run standard ccsd code.')
!
      if (wf%n_cc2_o < 0 .or. wf%n_ccs_o < 0 .or. &
          wf%n_cc2_v < 0 .or. wf%n_ccs_v < 0 ) call output%error_msg('negative orbital space size.')
!
      if (wf%n_ccsd_o +wf%n_cc2_o + wf%n_ccs_o .ne. wf%n_o) &
         call output%error_msg('occupied spaces do not add to n_o.')
!
      if (wf%n_ccsd_v + wf%n_cc2_v + wf%n_ccs_v .ne. wf%n_v) &
         call output%error_msg('virtual spaces do not add to n_v.')
!
      call wf%check_orthonormality_of_MOs()
!
   end subroutine check_orbital_space_mlccsd
!
!
   subroutine print_orbital_space_mlccsd(wf)
!!
!!    Print orbital space
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Jan 2019
!!
      implicit none
!
      class(mlccsd) :: wf
!
      call output%printf('m', '- MLCCSD orbital partitioning:',fs='(/t3,a)')
!
      call output%printf('m', 'Orbital type: ' // trim(wf%ccsd_orbital_type), fs='(/t6,a)')
!
      call output%printf('m', 'Number occupied ccsd orbitals: (i4)', &
                        ints=[wf%n_ccsd_o], fs='(/t6, a)')
      call output%printf('m', 'Number virtual ccsd orbitals:  (i4)', &
                        ints=[wf%n_ccsd_v], fs='(t6, a)')
      call output%printf('m', 'Number occupied cc2 orbitals:  (i4)', &
                        ints=[wf%n_cc2_o], fs='(/t6, a)')
      call output%printf('m', 'Number virtual cc2 orbitals:   (i4)', &
                        ints=[wf%n_cc2_v], fs='(t6, a)')
      call output%printf('m', 'Number occupied ccs orbitals:  (i4)', &
                        ints=[wf%n_ccs_o], fs='(/t6, a)')
      call output%printf('m', 'Number virtual ccs orbitals:   (i4)', &
                        ints=[wf%n_ccs_v], fs='(t6, a)')
!
   end subroutine print_orbital_space_mlccsd
!
!
   subroutine construct_x2_mlccsd(wf)
!!
!!    Construct X2
!!    Written by Sarai D. Folkestad, Jul 2019
!!
!!    Construct X2 in packed form and store in wf%x2
!!
!!
      use array_utilities, only: zero_array
      use reordering, only: packin
!
      implicit none
!
      class(mlccsd), intent(inout) :: wf
!
      real(dp), dimension(:,:,:,:), allocatable :: x_aibj
!
      integer :: a, i, b, j, ai, bj, aibj
!
      call mem%alloc(x_aibj, wf%n_ccsd_v + wf%n_cc2_v, &
                             wf%n_ccsd_o + wf%n_cc2_o, &
                             wf%n_ccsd_v + wf%n_cc2_v, &
                             wf%n_ccsd_o + wf%n_cc2_o)
!
      call zero_array(x_aibj, ((wf%n_ccsd_v + wf%n_cc2_v)**2)*((wf%n_ccsd_o + wf%n_cc2_o)**2))
!
      if (wf%do_cc2) call wf%construct_cc2_amplitudes(x_aibj)
!
!$omp parallel do private (j, b, a, i, bj, ai, aibj) collapse(2)
      do j = 1, wf%n_ccsd_o
         do b = 1, wf%n_ccsd_v
!
            bj = wf%n_ccsd_v*(j - 1) + b
!
            do i = 1, wf%n_ccsd_o
               do a = 1, wf%n_ccsd_v
!
                  ai = wf%n_ccsd_v*(i - 1) + a
!
                  aibj = max(ai,bj)*(max(ai,bj)-3)/2 + ai + bj
!
                  x_aibj(a, i, b, j) = x_aibj(a, i, b, j) + wf%t2(aibj)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call packin(wf%x2, x_aibj, (wf%n_ccsd_v + wf%n_cc2_v)*(wf%n_ccsd_o + wf%n_cc2_o))
!
      call mem%dealloc(x_aibj, wf%n_ccsd_v + wf%n_cc2_v, &
                               wf%n_ccsd_o + wf%n_cc2_o, &
                               wf%n_ccsd_v + wf%n_cc2_v, &
                               wf%n_ccsd_o + wf%n_cc2_o)
!
   end subroutine construct_x2_mlccsd
!
!
   subroutine construct_u_aibj_mlccsd(wf)
!!
!!    Construct u
!!    Written by Sarai D. Folkestad, Jul 2019
!!
!!    Note: assumes wf%x2 in memory
!!
      use reordering, only: squareup, add_1432_to_1234
      use array_utilities, only: copy_and_scale
!
      implicit none
!
      class(mlccsd), intent(inout) :: wf
!
      real(dp), dimension(:,:,:,:), allocatable :: x_aibj
!
      call mem%alloc(x_aibj,wf%n_ccsd_v + wf%n_cc2_v, &
                            wf%n_ccsd_o + wf%n_cc2_o, &
                            wf%n_ccsd_v + wf%n_cc2_v, &
                            wf%n_ccsd_o + wf%n_cc2_o)
!
      call squareup(wf%x2, x_aibj, (wf%n_ccsd_v + wf%n_cc2_v)*(wf%n_ccsd_o + wf%n_cc2_o))
!
      call copy_and_scale(two, x_aibj, wf%u_aibj, &
                     ((wf%n_ccsd_v + wf%n_cc2_v)**2)*((wf%n_ccsd_o + wf%n_cc2_o)**2))
!
      call add_1432_to_1234(-one, x_aibj, wf%u_aibj, wf%n_ccsd_v + wf%n_cc2_v, &
                                                wf%n_ccsd_o + wf%n_cc2_o, &
                                                wf%n_ccsd_v + wf%n_cc2_v, &
                                                wf%n_ccsd_o + wf%n_cc2_o)
!
      call mem%dealloc(x_aibj,wf%n_ccsd_v + wf%n_cc2_v, &
                              wf%n_ccsd_o + wf%n_cc2_o, &
                              wf%n_ccsd_v + wf%n_cc2_v, &
                              wf%n_ccsd_o + wf%n_cc2_o)
!
   end subroutine construct_u_aibj_mlccsd
!
!
   subroutine get_orbital_differences_mlccsd(wf, orbital_differences, N)
!!
!!    Get orbital differences
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 2019
!!
      implicit none
!
      class(mlccsd), intent(in) :: wf
!
      integer, intent(in) :: N
      real(dp), dimension(N), intent(out) :: orbital_differences
!
      integer :: a, i, ai, b, j, bj, aibj, n_a_o, n_a_v
!
      if ((N .ne. wf%n_gs_amplitudes) .and. (N .ne. wf%n_es_amplitudes)) &
         call output%error_msg('Could not recognize length of orbital differences vector')
!
      call wf%ccs%get_orbital_differences(orbital_differences(1:wf%n_t1), wf%n_t1)
!
      if (N .eq. wf%n_es_amplitudes) then
!
         n_a_o = wf%n_ccsd_o + wf%n_cc2_o
         n_a_v = wf%n_ccsd_v + wf%n_cc2_v
!
      else ! do_cc2 = .true. and N = wf%n_gs_amplitudes
!
         n_a_o = wf%n_ccsd_o
         n_a_v = wf%n_ccsd_v
!
      endif
!
!$omp parallel do schedule(static) private(a, i, b, j, ai, bj, aibj)
      do a = 1, n_a_v
         do i = 1, n_a_o
!
            ai = n_a_v*(i - 1) + a
!
            do j = 1, n_a_o
               do b = 1, n_a_v
!
                  bj = n_a_v*(j - 1) + b
!
                  if (ai .ge. bj) then
!
                     aibj = (ai*(ai-3)/2) + ai + bj
!
                     orbital_differences(aibj + (wf%n_o)*(wf%n_v)) = &
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
   end subroutine get_orbital_differences_mlccsd
!
!
   subroutine construct_cc2_amplitudes_mlccsd(wf, s_aibj)
!!
!!    Construct CC2 amplitudes
!!    Written by Sarai D. Folkestad, 2019
!!
!!    If wf%do_cc2 = .true., the CC2 amplitudes will
!!    be calculated in the basis where the occupied-occupied
!!    and virtual-virtual Fock matrices are diagonal in the
!!    CC2-CCSD blocks.
!!
!!       s_AIBJ = -g_AIBJ/e_AIBJ
!!
!!    The amplitudes are afterwards
!!    transformed to the MLCCSD basis
!!
!!       s_aibj = sum_AIBJ U_Aa U_Ii U_Bb U_Jj s_AIBJ.
!!
      use reordering, only: sort_123_to_132, sort_1234_to_2143
!
      implicit none
!
      class(mlccsd) :: wf
!
      real(dp), dimension(wf%n_ccsd_v + wf%n_cc2_v,wf%n_ccsd_o + wf%n_cc2_o, &
                          wf%n_ccsd_v + wf%n_cc2_v,wf%n_ccsd_o + wf%n_cc2_o), intent(out) :: s_aibj
!
      real(dp), dimension(:,:,:), allocatable :: L_J_ai
      real(dp), dimension(:,:,:), allocatable :: X_J_aj
      real(dp), dimension(:,:,:), allocatable :: X_J_ja
      real(dp), dimension(:,:,:), allocatable :: X_J_jb
!
      real(dp), dimension(:,:,:,:), allocatable :: s_iajb
!
      real(dp), dimension(:,:,:,:), allocatable :: X_kajb
      real(dp), dimension(:,:,:,:), allocatable :: X_kajd
      real(dp), dimension(:,:,:,:), allocatable :: X_akdj
      real(dp), dimension(:,:,:,:), allocatable :: X_ckdj
!
      integer :: n_a_o, n_a_v, last_a_o, last_a_v ! Total and last  active indices
!
      integer :: a, i, b, j
!
      n_a_o = wf%n_ccsd_o + wf%n_cc2_o
      n_a_v = wf%n_ccsd_v + wf%n_cc2_v
!
      last_a_o  = n_a_o
      last_a_v  = n_a_v
!
!     Read Cholesky vectors L^J_ai and transform them from the MLCCSD
!     basis to the basis for amplitude construction
!
      call mem%alloc(L_J_ai, wf%eri_t1%n_J, n_a_v, n_a_o)
!
      call wf%L_t1%get(L_J_ai, &
                           1 + wf%n_o,          &
                           last_a_v + wf%n_o,   &
                           1, last_a_o)
!
      call mem%alloc(X_J_aj, wf%eri_t1%n_J, n_a_v, n_a_o)
!
      call dgemm('N', 'T', &
                  wf%eri_t1%n_J*(n_a_v), &
                  n_a_o,                    &
                  n_a_o,                    &
                  one,                      &
                  L_J_ai,                   &
                  wf%eri_t1%n_J*(n_a_v), &
                  wf%O_o,                   &
                  n_a_o,                    &
                  zero,                     &
                  X_J_aj,                   &
                  wf%eri_t1%n_J*(n_a_v))
!
      call mem%dealloc(L_J_ai, wf%eri_t1%n_J, n_a_v, n_a_o)
!
      call mem%alloc(X_J_ja, wf%eri_t1%n_J, n_a_o, n_a_v)
!
      call sort_123_to_132(X_J_aj, X_J_ja, wf%eri_t1%n_J, n_a_v, n_a_o)
!
      call mem%dealloc(X_J_aj, wf%eri_t1%n_J, n_a_v, n_a_o)
!
      call mem%alloc(X_J_jb, wf%eri_t1%n_J, n_a_o, n_a_v)
!
      call dgemm('N', 'T', &
                  wf%eri_t1%n_J*(n_a_o), &
                  n_a_v,                    &
                  n_a_v,                    &
                  one,                      &
                  X_J_ja,                   &
                  wf%eri_t1%n_J*(n_a_o), &
                  wf%O_v,                   &
                  n_a_v,                    &
                  zero,                     &
                  X_J_jb,                   &
                  wf%eri_t1%n_J*(n_a_o))
!
      call mem%dealloc(X_J_ja, wf%eri_t1%n_J, n_a_o, n_a_v)
!
!     Construct the amplitudes in this basis
!     s_AIBJ = - g_AIBJ/ε_AIBJ
!
      call mem%alloc(s_iajb, n_a_o, n_a_v, n_a_o, n_a_v)

!
      call dgemm('T', 'N', &
                  (n_a_o)*(n_a_v),   &
                  (n_a_o)*(n_a_v),   &
                  wf%eri_t1%n_J,  &
                  -one,              &
                  X_J_jb,            &
                  wf%eri_t1%n_J,  &
                  X_J_jb,            &
                  wf%eri_t1%n_J,  &
                  zero,              &
                  s_iajb,            &
                  (n_a_o)*(n_a_v))
!
      call mem%dealloc(X_J_jb, wf%eri_t1%n_J, n_a_o, n_a_v)
!
!$omp parallel do private(b, j, a, i) collapse(2)
      do b = 1, n_a_v
         do j = 1, n_a_o
            do a = 1, n_a_v
               do i = 1, n_a_o
!
                  s_iajb(i, a, j, b) = s_iajb(i, a, j, b)/&
                        (wf%orbital_energies_cc2(wf%n_o + a) &
                        +  wf%orbital_energies_cc2(wf%n_o + b) &
                        -  wf%orbital_energies_cc2(i) &
                        -  wf%orbital_energies_cc2(j) )
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
!     Transform back to the MLCCSD basis
!
      call mem%alloc(X_kajb, n_a_o, n_a_v, n_a_o, n_a_v)
!
      call dgemm('T', 'N',          &
                  n_a_o,            &
                  n_a_o*(n_a_v**2), &
                  n_a_o,            &
                  one,              &
                  wf%O_o,           &
                  n_a_o,            &
                  s_iajb,           &
                  n_a_o,            &
                  zero,             &
                  X_kajb,           &
                  n_a_o)
!
      call mem%dealloc(s_iajb, n_a_o, n_a_v, n_a_o, n_a_v)
!
      call mem%alloc(X_kajd, n_a_o, n_a_v, n_a_o, n_a_v)
!
      call dgemm('N', 'N',          &
                  (n_a_o**2)*n_a_v, &
                  n_a_v,            &
                  n_a_v,            &
                  one,              &
                  X_kajb,           &
                  (n_a_o**2)*n_a_v, &
                  wf%O_v,           &
                  n_a_v,            &
                  zero,             &
                  X_kajd,           &
                  (n_a_o**2)*n_a_v)
!
      call mem%dealloc(X_kajb, n_a_o, n_a_v, n_a_o, n_a_v)
!
      call mem%alloc(X_akdj, n_a_v, n_a_o, n_a_v, n_a_o)
!
      call sort_1234_to_2143(X_kajd, X_akdj, n_a_o, n_a_v, n_a_o, n_a_v)
!
      call mem%dealloc(X_kajd, n_a_o, n_a_v, n_a_o, n_a_v)
!
      call mem%alloc(X_ckdj, n_a_v, n_a_o, n_a_v, n_a_o)
!
      call dgemm('T', 'N',          &
                  n_a_v,            &
                  n_a_v*(n_a_o**2), &
                  n_a_v,            &
                  one,              &
                  wf%O_v,           &
                  n_a_v,            &
                  X_akdj,           &
                  n_a_v,            &
                  zero,             &
                  X_ckdj,           &
                  n_a_v)
!
      call mem%dealloc(X_akdj, n_a_v, n_a_o, n_a_v, n_a_o)
!
      call dgemm('N', 'N',          &
                  (n_a_v**2)*n_a_o, &
                  n_a_o,            &
                  n_a_o,            &
                  one,              &
                  X_ckdj,           &
                  (n_a_v**2)*n_a_o, &
                  wf%O_o,           &
                  n_a_o,            &
                  zero,             &
                  s_aibj,           &
                  (n_a_v**2)*n_a_o)
!
      call mem%dealloc(X_ckdj, n_a_v, n_a_o, n_a_v, n_a_o)
!
   end subroutine construct_cc2_amplitudes_mlccsd
!
!
   subroutine calculate_energy_mlccsd(wf)
!!
!!    Calculate energy (MLCCSD)
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad,
!!    Andreas Skeidsvoll, 2018
!!
!!    Calculates the CCSD energy. This is only equal to the actual
!!    energy when the ground state equations are solved, of course.
!!
!!       E = E_hf + sum_aibj (t_ij^ab + t_i^a t_j^b) L_iajb
!!
!!    Modified for MLCCSD by Sarai D. Folkestad, 2019
!!
      use reordering, only: squareup
!
      implicit none
!
      class(mlccsd), intent(inout) :: wf
!
      real(dp), dimension(:,:,:,:), allocatable :: g_iajb, x_aibj
!
      real(dp) :: correlation_energy
!
      integer :: a, i, b, j
      integer :: n_doubles_o, n_doubles_v
!
      type(timings) :: timer
!
      timer = timings('Calculate energy', pl='n')
      call timer%turn_on()
!
      n_doubles_v = wf%n_ccsd_v + wf%n_cc2_v
      n_doubles_o = wf%n_ccsd_o + wf%n_cc2_o
!
      call wf%ccs%calculate_energy()
!
      correlation_energy = zero
!
      call mem%alloc(x_aibj, n_doubles_v, n_doubles_o, n_doubles_v, n_doubles_o)
      call wf%construct_x2()
      call squareup(wf%x2, x_aibj, n_doubles_v*n_doubles_o)
!
      call mem%alloc(g_iajb, n_doubles_o, n_doubles_v, n_doubles_o, n_doubles_v)
      call wf%eri_t1%get('ovov', g_iajb, &
                         1, n_doubles_o, &
                         1, n_doubles_v, &
                         1, n_doubles_o, &
                         1, n_doubles_v)
!
!$omp parallel do private(a,i,j,b) reduction(+:correlation_energy)
      do a = 1, n_doubles_v
         do i = 1, n_doubles_o
            do j = 1, n_doubles_o
               do b = 1, n_doubles_v
!
                  correlation_energy = correlation_energy + &
                                    x_aibj(a,i,b,j)*(two*g_iajb(i, a, j, b) &
                                    - g_iajb(i, b, j, a))
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(g_iajb, n_doubles_o, n_doubles_v, n_doubles_o, n_doubles_v)
      call mem%dealloc(x_aibj, n_doubles_v, n_doubles_o, n_doubles_v, n_doubles_o)
!
      wf%energy = wf%energy + correlation_energy
      wf%correlation_energy = wf%correlation_energy + correlation_energy
!
      call timer%turn_off()
!
   end subroutine calculate_energy_mlccsd
!
!
   subroutine set_initial_amplitudes_guess_mlccsd(wf, restart)
!!
!!    Set initial amplitudes guess
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!    Adapted by Alexander C. Paul to use the restart logical, Oct 2020
!!
      use array_utilities, only: zero_array
!
      implicit none
!
      class(mlccsd), intent(inout) :: wf
!
      logical, intent(in)        :: restart
!
      integer :: n_amplitudes_read
!
      if (.not. restart) then
!
         call zero_array(wf%t1, wf%n_t1)
         call wf%set_t2_to_cc2_guess()
!
      else
!
         if (wf%t_file%exists()) then
!
            call output%printf('m', 'Requested restart. Reading in solution from file.', &
                          fs='(/t3,a)')
!
            call wf%read_amplitudes(n_amplitudes_read)
!
            if(n_amplitudes_read == wf%n_gs_amplitudes) then
!
               call wf%construct_t1_cholesky(wf%t1, wf%L_mo, wf%L_t1)
!
            else if (n_amplitudes_read == wf%n_t1) then
!
               call wf%construct_t1_cholesky(wf%t1, wf%L_mo, wf%L_t1)
!
               call wf%set_t2_to_cc2_guess()
!
            else
!
               call output%error_msg('Did not recognize number of t-amplitudes on file &
                                     &expected (i0) or (i0) found (i0)', &
                                     ints=[wf%n_gs_amplitudes, wf%n_t1, n_amplitudes_read])
!
            end if
!
         else
!
            call zero_array(wf%t1, wf%n_t1)
!
            call wf%set_t2_to_cc2_guess()
!
         end if
!
      end if
!
   end subroutine set_initial_amplitudes_guess_mlccsd
!
!
   subroutine set_t2_to_cc2_guess_mlccsd(wf)
!!
!!    Set t2 amplitudes guess
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
!!    t_aibj = - g_aibj/ε_aibj
!!
!!    Modified for MLCCSD by Sarai D. Folkestad, Nov 2019
!!
      implicit none
!
      class(mlccsd) :: wf
!
      real(dp), dimension(:,:,:,:), allocatable :: g_aibj
!
      integer :: a, b, i, j, ai, bj, aibj
!
      call mem%alloc(g_aibj, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
      call wf%eri_t1%get('vovo', g_aibj, 1, wf%n_ccsd_v, 1, wf%n_ccsd_o, &
                                             1, wf%n_ccsd_v, 1, wf%n_ccsd_o)
!
!$omp parallel do schedule(static) private(a, i, b, j, ai, bj, aibj)
      do a = 1, wf%n_ccsd_v
         do i = 1, wf%n_ccsd_o
!
            ai = wf%n_ccsd_v*(i-1) + a
!
            do j = 1, wf%n_ccsd_o
               do b = 1, wf%n_ccsd_v
!
                  bj = wf%n_ccsd_v*(j-1) + b
!
                  if (ai .ge. bj) then
!
                     aibj = (ai*(ai-3)/2) + ai + bj
!
                     wf%t2(aibj) = g_aibj(a,i,b,j)/&
                                          (wf%orbital_energies(i) +          &
                                          wf%orbital_energies(j)  -          &
                                          wf%orbital_energies(a + wf%n_o) -  &
                                          wf%orbital_energies(b + wf%n_o))
!
                  endif
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(g_aibj, wf%n_ccsd_v, wf%n_ccsd_o, wf%n_ccsd_v, wf%n_ccsd_o)
!
   end subroutine set_t2_to_cc2_guess_mlccsd
!
!
   subroutine form_newton_raphson_t_estimate_mlccsd(wf, t, dt)
!!
!!    Form Newton-Raphson t estimate
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
!!    Here, t is the full amplitude vector and dt is the correction to the amplitude vector.
!!
!!    The correction is assumed to be obtained from either
!!    solving the Newton-Raphson equation
!!
!!       A dt = -omega,
!!
!!    where A and omega are given in the biorthonormal basis,
!!    or from the quasi-Newton equation (A ~ diagonal with diagonal = epsilon)
!!
!!        dt = -omega/epsilon
!!
!!    Epsilon is the vector of orbital differences.
!!
!!    On exit, t = t + dt, where the appropriate basis change has been accounted
!!    for (in particular for the double amplitudes in CCSD wavefunctions). Also,
!!    dt is expressed in the basis compatible with t.
!!
      implicit none
!
      class(mlccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: dt
      real(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: t
!
      integer :: ai, aiai
!
!     Change dt doubles diagonal to match the definition of the
!     double amplitudes
!
!$omp parallel do private(ai, aiai)
      do ai = 1, wf%n_ccsd_o*wf%n_ccsd_v
!
         aiai = ai*(ai - 3)/2 + 2*ai
         dt(wf%n_t1 + aiai) = two*dt(wf%n_t1 + aiai)
!
      enddo
!$omp end parallel do
!
!     Add the dt vector to the t vector
!
      call daxpy(wf%n_gs_amplitudes, one, dt, 1, t, 1)
!
   end subroutine form_newton_raphson_t_estimate_mlccsd
!
!
   subroutine cleanup_mlccsd(wf)
!!
!!    Cleanup
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad and
!!    Alexander C. Paul , 2018
!!
      use string_utilities, only: convert_to_uppercase
!
      implicit none
!
      class(mlccsd) :: wf
!
      call output%printf('v', '- Cleaning up (a0) wavefunction', &
                         chars=[trim(convert_to_uppercase(wf%name_))], fs='(/t3, a)')
!
      call wf%destruct_amplitudes()
      call wf%destruct_multipliers()
!
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
      call wf%destruct_nto_states()
      call wf%destruct_cnto_states()
!
      call wf%destruct_orbital_energies_cc2()
      call wf%destruct_orbital_coefficients_cc2()
!
      call wf%destruct_O_o()
      call wf%destruct_O_v()
!
      call wf%destruct_mo_fock()
      call wf%destruct_ao_fock()
!
      deallocate(wf%ao)
      if (wf%embedded) deallocate(wf%embedding)
!
!     To avoid memory leaks with intel, explicit deallocations
      deallocate(wf%eri_t1)
      deallocate(wf%L_mo)
      deallocate(wf%L_t1)
!
   end subroutine cleanup_mlccsd
!
!
   subroutine print_X1_diagnostics_mlccsd(wf, X, label)
!!
!!    Print X1 diagnostics
!!    Written by Sarai D. Folkestad, Nov 2019
!!
      use array_utilities, only: get_l2_norm
!
      implicit none
!
      class(mlccsd), intent(in) :: wf
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
      call mem%alloc(X_internal, wf%n_ccsd_v*wf%n_ccsd_o)
!
      do i = 1, wf%n_ccsd_o
         do a = 1, wf%n_ccsd_v
!
            ai = wf%n_ccsd_v*(i - 1) + a
            ai_full = wf%n_v*(i - 1) + a
!
            X_internal(ai) = X(ai_full)
!
         enddo
      enddo
!
      internal_fraction = get_l2_norm(X_internal, wf%n_ccsd_v*wf%n_ccsd_o)&
                           /get_l2_norm(X,wf%n_es_amplitudes)
!
      call output%printf('n', 'MLCC diagnostics:', fs='(/t6,a)')
!
      call output%printf('n', '|(a0)1^internal|/|(a0)| =  (f19.12)', &
            reals=[internal_fraction], chars=[label, label], fs='(/t6,a)')
!
      internal_fraction = get_l2_norm(X_internal, wf%n_ccsd_v*wf%n_ccsd_o)&
                           /get_l2_norm(X(1:wf%n_t1), wf%n_t1)
!
      call output%printf('n', '|(a0)1^internal|/|(a0)1| = (f19.12)', &
            reals=[internal_fraction], chars=[label, label], fs='(t6,a)')
!
      call mem%dealloc(X_internal, wf%n_ccsd_v*wf%n_ccsd_o)
!
   end subroutine print_X1_diagnostics_mlccsd
!
!
   subroutine get_cvs_projector_mlccsd(wf, projector, n_cores, core_MOs)
!!
!!    Get CVS projector
!!    Written by Sarai D. Folkestad, Oct 2018
!!
      use array_utilities, only: zero_array
!
      implicit none
!
      class(mlccsd), intent(inout) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: projector
!
      integer, intent(in) :: n_cores
!
      integer, dimension(n_cores), intent(in) :: core_MOs
!
      integer :: core, i, a, ai, j, b, bj, aibj, n_a_o, n_a_v
!
      call zero_array(projector, wf%n_es_amplitudes)
!
      do core = 1, n_cores
!
         i = core_MOs(core)
!
         if (i  .gt. wf%n_ccsd_o) then
            call output%error_msg('Core orbital (i0) is not CCSD orbital', ints=[i])
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
      n_a_o = wf%n_ccsd_o + wf%n_cc2_o
      n_a_v = wf%n_ccsd_v + wf%n_cc2_v
!
!$omp parallel do private (a, ai, j, b, bj, aibj)
         do a = 1, n_a_v
!
            ai = n_a_v*(i - 1) + a
!
            do j = 1, n_a_o
               do b = 1, n_a_v
!
                  bj = n_a_v*(j - 1) + b
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
   end subroutine get_cvs_projector_mlccsd
!
!
   subroutine set_cvs_start_indices_mlccsd(wf, start_indices)
!!
!!    Set CVS start indices
!!    Written by Sarai D. Folkestad
!!
      implicit none
!
      class(mlccsd), intent(in) :: wf
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
         if (wf%core_MOs(mo) .gt. wf%n_ccsd_o) &
            call output%error_msg('Active MO not in active space for MLCCSD calculation.')
!
      enddo
!
      call wf%ccs%set_cvs_start_indices(start_indices)
!
   end subroutine set_cvs_start_indices_mlccsd
!
!
   subroutine mo_preparations_from_restart_mlccsd(wf)
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
      class(mlccsd) :: wf
!
      real(dp), dimension(:,:), allocatable :: canonical_orbitals
      real(dp), dimension(:,:), allocatable :: T
!
      call output%warning_msg('Number of orbitals in active and inactive spaces read &
                              &from orbital_coefficients_mlcc. Make sure restarted &
                              &calculation and restart files are consistent!')
!
      call mem%alloc(canonical_orbitals, wf%ao%n, wf%n_mo)
      call dcopy(wf%ao%n*wf%n_mo, wf%orbital_coefficients, 1, canonical_orbitals, 1)
!
      call wf%read_mlcc_orbitals()
!
      call mem%alloc(T, wf%n_mo, wf%n_mo)
      call wf%construct_mo_basis_transformation(wf%orbital_coefficients, canonical_orbitals, T)
      call wf%L_mo%basis_transformation(T)
      call mem%dealloc(T, wf%n_mo, wf%n_mo)
!
      call wf%L_t1%set_equal_to(wf%L_mo)
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
   end subroutine mo_preparations_from_restart_mlccsd
!
!
   subroutine scale_amplitudes_mlccsd(wf, t)
!!
!!    Scale amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
!!    Scales t to conform with the convention used in the wavefunction:
!!
!!       t1 <- t1
!!       t2_aiai <- two * t2_aiai
!!       ...
!!
      use array_utilities, only: scale_diagonal
!
      implicit none
!
      class(mlccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: t
!
      call scale_diagonal(two,                                    &
                          t(wf%n_t1 + 1 : wf%n_gs_amplitudes),    &
                          wf%n_ccsd_o*wf%n_ccsd_v)
!
   end subroutine scale_amplitudes_mlccsd
!
!
end module mlccsd_class
