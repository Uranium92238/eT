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
   use mlcc2_class
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
      procedure :: construct_block_diagonal_fock_orbitals &
                  => construct_block_diagonal_fock_orbitals_mlccsd
!
      procedure :: construct_cholesky_orbitals  => construct_cholesky_orbitals_mlccsd
      procedure :: construct_paos               => construct_paos_mlccsd
!      
      procedure :: construct_cc2_cnto_transformation_matrices &
                  => construct_cc2_cnto_transformation_matrices_mlccsd
!
      procedure :: cc2_calculation_for_cntos    => cc2_calculation_for_cntos_mlccsd
!
!     Omega
!
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
      procedure :: get_gs_orbital_differences   => get_gs_orbital_differences_mlccsd
      procedure :: construct_cc2_amplitudes     => construct_cc2_amplitudes_mlccsd
!
      procedure :: calculate_energy             => calculate_energy_mlccsd
      procedure :: set_initial_amplitudes_guess => set_initial_amplitudes_guess_mlccsd
      procedure :: set_t2_to_mp2_guess          => set_t2_to_mp2_guess_mlccsd
!
      procedure :: get_amplitudes               => get_amplitudes_mlccsd
      procedure :: set_amplitudes               => set_amplitudes_mlccsd
!
      procedure :: form_newton_raphson_t_estimate &
                  => form_newton_raphson_t_estimate_mlccsd
!
!     Initialize/destruct
!
      procedure :: initialize_u_aibj            => initialize_u_aibj_mlccsd
      procedure :: destruct_u_aibj              => destruct_u_aibj_mlccsd
      procedure :: initialize_t2                => initialize_t2_mlccsd
      procedure :: destruct_t2                  => destruct_t2_mlccsd
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
      procedure :: initialize_O_o               => initialize_O_o_mlccsd
      procedure :: destruct_O_o                 => destruct_O_o_mlccsd
      procedure :: initialize_O_v               => initialize_O_v_mlccsd
      procedure :: destruct_O_v                 => destruct_O_v_mlccsd
!
!
!     Read/save
!
      procedure :: read_amplitudes              => read_amplitudes_mlccsd
      procedure :: save_amplitudes              => save_amplitudes_mlccsd
!
!     Cleanup 
!
      procedure :: cleanup                      => cleanup_mlccsd
!
!     Initialize
!
      procedure :: initialize                   => initialize_mlccsd
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
!
   end interface 
!
contains
!
!
   subroutine initialize_mlccsd(wf, template_wf)
!!
!!    New mlccsd
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Adapted by Sarai D. Folkestad from CCS constructer, 2020
!!
      implicit none
!
      class(mlccsd), intent(inout) :: wf 
!
      class(wavefunction), intent(in) :: template_wf
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
      wf%cholesky_orbital_threshold = 1.0d-2
!
      call wf%general_cc_preparations()
      call wf%set_variables_from_template_wf(template_wf)
      call wf%print_banner()
!
      if (wf%bath_orbital) call output%error_msg('Bath orbitals not yet implemented for MLCCSD')
!      
      call wf%read_mlcc_settings()
!
      wf%n_t1 = (wf%n_o)*(wf%n_v)
!
      if (wf%do_cc2) then
!
         call wf%initialize_orbital_energies_cc2()
         call wf%initialize_orbital_coefficients_cc2()
!
      endif
!
      call wf%initialize_fock()
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
      implicit none
!
      class(mlccsd), intent(inout) :: wf
!
      logical :: do_ccsd
!
      if (.not. input%requested_section('mlcc')) &
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
!     Get orbital types
!
      if (wf%do_cc2 .and. wf%do_ccs) &
            call input%get_required_keyword_in_section('cc2 orbitals', 'mlcc', wf%cc2_orbital_type)
!     
      call input%get_required_keyword_in_section('ccsd orbitals', 'mlcc', wf%ccsd_orbital_type)
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
      call wf%general_mlcc_mo_preparations()
!
      if (wf%do_cc2) then
!
         call wf%initialize_O_o()
         call wf%initialize_O_v()
!
         call wf%construct_mlccsd_basis_transformation_matrix()
!
      endif
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
   subroutine get_gs_orbital_differences_mlccsd(wf, orbital_differences, N)
!!
!!    Get ground state orbital differences 
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 2019
!!
!!    Calculate the orbital differences for the ground state calculation 
!!
      implicit none
!
      class(mlccsd), intent(in) :: wf
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
      do a = 1, wf%n_ccsd_v
         do i = 1, wf%n_ccsd_o
!
            ai = wf%n_ccsd_v*(i - 1) + a
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
   end subroutine get_gs_orbital_differences_mlccsd
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
      call mem%alloc(L_J_ai, wf%integrals%n_J, n_a_v, n_a_o)
!
      call wf%integrals%get_cholesky_t1(L_J_ai, &
                           1 + wf%n_o,          &
                           last_a_v + wf%n_o,   &
                           1, last_a_o)
!
      call mem%alloc(X_J_aj, wf%integrals%n_J, n_a_v, n_a_o)
!
      call dgemm('N', 'T', &
                  wf%integrals%n_J*(n_a_v), &
                  n_a_o,                    &
                  n_a_o,                    &
                  one,                      &
                  L_J_ai,                   &
                  wf%integrals%n_J*(n_a_v), &
                  wf%O_o,                   &
                  n_a_o,                    &
                  zero,                     &
                  X_J_aj,                   &
                  wf%integrals%n_J*(n_a_v))
!
      call mem%dealloc(L_J_ai, wf%integrals%n_J, n_a_v, n_a_o)
!
      call mem%alloc(X_J_ja, wf%integrals%n_J, n_a_o, n_a_v)
!
      call sort_123_to_132(X_J_aj, X_J_ja, wf%integrals%n_J, n_a_v, n_a_o)
!
      call mem%dealloc(X_J_aj, wf%integrals%n_J, n_a_v, n_a_o)
!
      call mem%alloc(X_J_jb, wf%integrals%n_J, n_a_o, n_a_v)
!
      call dgemm('N', 'T', &
                  wf%integrals%n_J*(n_a_o), &
                  n_a_v,                    &
                  n_a_v,                    &
                  one,                      &
                  X_J_ja,                   &
                  wf%integrals%n_J*(n_a_o), &
                  wf%O_v,                   &
                  n_a_v,                    &
                  zero,                     &
                  X_J_jb,                   &
                  wf%integrals%n_J*(n_a_o))
!
      call mem%dealloc(X_J_ja, wf%integrals%n_J, n_a_o, n_a_v)
!
!     Construct the amplitudes in this basis
!     s_AIBJ = - g_AIBJ/ε_AIBJ
!
      call mem%alloc(s_iajb, n_a_o, n_a_v, n_a_o, n_a_v)

!
      call dgemm('T', 'N', &
                  (n_a_o)*(n_a_v),   &
                  (n_a_o)*(n_a_v),   &
                  wf%integrals%n_J,  &
                  -one,              &
                  X_J_jb,            &
                  wf%integrals%n_J,  &
                  X_J_jb,            &
                  wf%integrals%n_J,  &
                  zero,              &
                  s_iajb,            &
                  (n_a_o)*(n_a_v))
!
      call mem%dealloc(X_J_jb, wf%integrals%n_J, n_a_o, n_a_v)
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
      n_doubles_v = wf%n_ccsd_v + wf%n_cc2_v
      n_doubles_o = wf%n_ccsd_o + wf%n_cc2_o
!
      call mem%alloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_ovov(g_iajb)
!
      correlation_energy = zero
!
!$omp parallel do private(a,i,b,j) reduction(+:correlation_energy)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
            do j = 1, wf%n_o
               do b = 1, wf%n_v
!
                  correlation_energy = correlation_energy + &
                                 (wf%t1(a,i))*(wf%t1(b,j))* &
                                 (two*g_iajb(i,a,j,b) - g_iajb(i,b,j,a))
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%alloc(x_aibj, n_doubles_v, n_doubles_o, n_doubles_v, n_doubles_o)
      call wf%construct_x2()
      call squareup(wf%x2, x_aibj, n_doubles_v*n_doubles_o)
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
      call mem%dealloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%dealloc(x_aibj, n_doubles_v, n_doubles_o, n_doubles_v, n_doubles_o)
!
      wf%energy = wf%hf_energy + correlation_energy
      wf%correlation_energy = correlation_energy
!
   end subroutine calculate_energy_mlccsd
!
!
   subroutine set_initial_amplitudes_guess_mlccsd(wf)
!!
!!    Set initial amplitudes guess
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(mlccsd) :: wf
!
      call zero_array(wf%t1, wf%n_t1)
!
      call wf%set_t2_to_mp2_guess()
!
   end subroutine set_initial_amplitudes_guess_mlccsd
!
!
   subroutine set_t2_to_mp2_guess_mlccsd(wf)
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
      call wf%get_vovo(g_aibj, &
                        1, wf%n_ccsd_v, &  
                        1, wf%n_ccsd_o, &  
                        1, wf%n_ccsd_v, &  
                        1, wf%n_ccsd_o)  
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
   end subroutine set_t2_to_mp2_guess_mlccsd
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
      call wf%integrals%cleanup()
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
      if (allocated(wf%l_files)) call wf%l_files%finalize_storer()
      if (allocated(wf%r_files)) call wf%r_files%finalize_storer()
!
   end subroutine cleanup_mlccsd
!
!
end module mlccsd_class
