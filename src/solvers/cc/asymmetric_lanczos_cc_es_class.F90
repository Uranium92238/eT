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
module asymmetric_lanczos_cc_es_class
!
!!
!!    Asymmetric Lanczos excited state solver class
!!    Written by Torsha Moitra and Sarai D. Folkestad
!!    and Sonia Coriani, Sep-Nov 2019
!!
!!    Calculates the excitation energies and oscillator strengths
!!    by use of an asymmetric Lanczos procedure.
!!
!!    A tridiagonal matrix T, of dimension n (= 'chain length'),
!!    is constructed through an asymmetric Lanczos procedure
!!    such that
!!
!!       T = P^T A Q,
!!
!!    where A is the Jacobian matrix and the vectors P and Q
!!    satisfy
!!
!!       P^T Q = I.
!!
!!    Note that the dimension of the A matrix is greater or equal to
!!    n.
!!
!!    The eigenvalues of A may then be approximated by the
!!    eigenvalues of T
!!
!!       L^T_n T R_n = e
!!
!!    and the oscillator strengths may be calculated through
!!
!!       f(j)^q = (2/3) * omega(j) norm_eta^q norm_xi^q L^T_n(j1) R_n(1j)
!!
!!    See for further details:
!!       S. Coriani et. al., J. Chem. Theory Comput. 2012, 8, 5, 1616-1628
!!
   use parameters
!
   use abstract_solver_class, only: abstract_solver
!
   use ccs_class,             only: ccs
   use timings_class,         only: timings
   use global_out,            only: output
   use memory_manager_class,  only: mem
!
   use abstract_projection_tool_class, only: abstract_projection_tool
   use es_cvs_projection_tool_class,   only: es_cvs_projection_tool
   use null_projection_tool_class,     only: null_projection_tool

   implicit none

   type, extends(abstract_solver) :: asymmetric_lanczos_cc_es
!
      character (len=100) :: name_
      character (len=100) :: tag
      character (len=100) :: author
      character (len=300) :: description1
      character (len=200) :: description2
      character (len=200) :: normalization
!
      integer :: chain_length
!
      character (len=40) :: es_type
!
      real(dp), dimension(:), allocatable :: energies
!
      class(abstract_projection_tool), allocatable :: projector
!
      class(ccs), pointer :: wf
!
      logical :: restart
!
   contains
!
      procedure :: run     => run_asymmetric_lanczos_cc_es
      procedure :: cleanup => cleanup_asymmetric_lanczos_cc_es
!
      procedure, private :: initialize_energies
      procedure, private :: destruct_energies
      procedure, private :: read_settings
      procedure, private :: print_settings
      procedure, private :: calculate_oscillator_strength
      procedure, private :: prepare_wf_for_excited_state
      procedure, private :: print_summary
      procedure, private :: print_banner
!
   end type asymmetric_lanczos_cc_es
!
   interface asymmetric_lanczos_cc_es
!
      procedure :: new_asymmetric_lanczos_cc_es
!
   end interface asymmetric_lanczos_cc_es
!
contains
!
!
   function new_asymmetric_lanczos_cc_es(wf) result(this)
!!
!!    New asymmetric lanczos CC ES
!!    Written by Torsha Moitra and Sarai D. Folkestad, Sep 2019
!!
!!    Constructor for the asymmetric lanczos excited state solver.
!!
!!       - Sets printables and prints the banner
!!       - Sets defaults, reads settings amd prints the settings
!!       - Prepares the wavefunction for Jacobian
!!          and Jacobian transposetransformation
!!
      use string_utilities, only: convert_to_uppercase
!
      implicit none
!
      class(ccs), intent(inout), target :: wf
!
      type(asymmetric_lanczos_cc_es) :: this
!
      this%total_timer = timings('Asymmetric Lanczos ES solver')
      call this%total_timer%turn_on()
!
!     Set printables
!
      this%name_  = 'Asymmetric Lanczos excited state solver'
      this%tag    = 'asymmetric lanczos'
!
      this%description1 = 'An asymmetric Lanczos solver builds a reduced space&
               & tridiagonal representation of the CC Jacobian of dimension defined&
               & by the chain length.&
               & Diagonalization of this representation gives the eigenvalues and eigenvectors.'
!
      this%description2 = 'A complete description of the algorithm can be found in &
               & S.Coriani et al., J. Chem. Theory Comput. 2012, 8, 5, 1616-1628.'
!
      call this%print_banner()
!
      this%wf => wf
!
!     Set defaults
!
      this%normalization = 'asymmetric'
      this%es_type       = 'valence'
!
      call this%read_settings()
!
!     Chain length cannot exceed the number of excited state amplitudes
!     (actual dimension of A in the basis)
!
      if(this%chain_length .gt. wf%n_es_amplitudes) then
!
         this%chain_length = wf%n_es_amplitudes
!
         call output%warning_msg('Chain length asked for is reset to: (i0)', &
                                 ints=[this%chain_length])
!
      endif
!
      call this%print_settings()
!
      call this%initialize_energies()
!
!     Prepare wavefunction for Jacobian transformations
!
      call this%prepare_wf_for_excited_state()
!
      if (trim(this%es_type) == 'valence') then
!
         this%projector = null_projection_tool(wf%n_es_amplitudes)
!
      elseif (trim(this%es_type) == 'core') then
!
         call wf%read_cvs_settings()
         this%projector = es_cvs_projection_tool(wf)
!
      else
!
         call output%error_msg('did not recognize the excited state type: ' //trim(this%es_type))
!
      endif
!
   end function new_asymmetric_lanczos_cc_es
!
!
   subroutine run_asymmetric_lanczos_cc_es(this)
!!
!!    Run
!!    Written by Torsha Moitra, Sarai D. Folkestad and Sonia Coriani 2019
!!
!!    Calculate excitation energies and oscillator strengths:
!!
!!       - The initial vectors p and q (seeds) are set
!!          (NOTE: these are now always eta and xi)
!!
!!       - The Krylov subspace is built using the asymmetric Lanczos
!!          procedure. The tridiagonal matrix T is built.
!!
!!       - T is diagonalized and the oscillator strengths are calculated.
!!
!!
!!    Modified by Andreas Skeidvoll and Torsha Moitra to add restart.
!!
      use asymmetric_lanczos_tool_class,  only : asymmetric_lanczos_tool
!
      implicit none

      class(asymmetric_lanczos_cc_es), intent(inout) :: this
!
      type(asymmetric_lanczos_tool), allocatable :: lanczos
!
      real(dp), dimension(:,:,:), allocatable   :: dipole_length
      real(dp), dimension(:,:), allocatable     :: reduced_left, reduced_right
      real(dp), dimension(:), allocatable       :: etaX, xiX
      real(dp), dimension(:), allocatable       :: p_and_q
      real(dp), dimension(:), allocatable       :: Aq, pA
      real(dp), dimension(:), allocatable       :: eigenvalues_Re, eigenvalues_Im
      real(dp), dimension(:), allocatable       :: oscillator_strength
!
      real(dp) :: norm_eta_times_norm_xi, ddot, overlap, factor, norm_p1q1
!
      integer :: cartesian, iteration, start_chain_length
!
      character(len=1), dimension(3), parameter :: cartesian_string = ['X', 'Y', 'Z']
!
!     Construct the dipole integrals
!
      call mem%alloc(dipole_length, this%wf%n_mo, this%wf%n_mo, 3)
      call this%wf%get_t1_oei('dipole', dipole_length)
!
      call mem%alloc(xiX, this%wf%n_es_amplitudes)
      call mem%alloc(etaX, this%wf%n_es_amplitudes)
!
      do cartesian=1, 3
!
!        Construct the eta and xi vectors
!
         call this%wf%construct_xiX(dipole_length(:,:,cartesian), xiX)
         call this%wf%construct_eom_etaX(dipole_length(:,:,cartesian), xiX, etaX)
!
         call this%projector%project(xiX)
         call this%projector%project(etaX)
!
!        Prepare initial seeds (p1,q1) as binormalized etaX and xiX vectors
!
!        The binormalization procedure is either 'symmetric' or 'asymmetric'
!
         overlap = ddot(this%wf%n_es_amplitudes, etaX, 1, xiX, 1)
!
         if(this%normalization=="symmetric")then
!
            norm_p1q1 = sqrt(abs(overlap))
!
            factor = one/norm_p1q1
!
            call dscal(this%wf%n_es_amplitudes, factor, xiX, 1)
            call dscal(this%wf%n_es_amplitudes, sign(factor,overlap), etaX, 1)
!
         else if (this%normalization == "asymmetric") then
!
            factor = one/overlap
!
            call dscal(this%wf%n_es_amplitudes, factor, etaX, 1)
!
         else
!
            call output%error_msg('did not recognize the binormalization procedure')
!
         endif
!
!        Product of the norm of eta and the norm of xi, to be used for oscillator strengths
!
         norm_eta_times_norm_xi = overlap
!
         lanczos = asymmetric_lanczos_tool(this%wf%n_es_amplitudes, this%chain_length,&
                                           this%normalization, this%restart, &
                                           cartesian_string(cartesian))
!
         call lanczos%initialize()
!
         if (.not. this%restart) then
!
            call lanczos%save_p(etaX, 1)
            call lanczos%save_q(xiX, 1)
!
         endif
!
         call mem%alloc(p_and_q, this%wf%n_es_amplitudes)
!
         call mem%alloc(Aq, this%wf%n_es_amplitudes)
         call mem%alloc(pA, this%wf%n_es_amplitudes)
!
!        Build the Krylov subspace and the tridiagonal matrix T
!        Loops over the chain-length
!
         start_chain_length = lanczos%get_start_chain_length()
         do iteration = start_chain_length, this%chain_length
!
!           Transforms p -> pA and q -> Aq
!
!           q is stored in p_and_q
            call lanczos%get_q(p_and_q, iteration)
            call this%wf%jacobian_transformation(p_and_q, Aq)
!
!           p is stored in p_and_q
            call lanczos%get_p(p_and_q, iteration)
            call this%wf%jacobian_transpose_transformation(p_and_q, pA)
!
            call this%projector%project(Aq)
            call this%projector%project(pA)
!
            call lanczos%expand_subspace(pA, Aq, iteration)
!
!           Checks if the value of beta is too small then stops iteration
!           Resets chain length to current iteration.
            if (lanczos%chain_length .lt. this%chain_length) then
!
               call output%warning_msg('Chain length value reset = (i4)', &
                  & ints=[lanczos%chain_length], fs='(/t6,a)')
!
               exit
!
            endif
!
         enddo
!
         this%chain_length = lanczos%chain_length
!
         call mem%dealloc(Aq, this%wf%n_es_amplitudes)
         call mem%dealloc(pA, this%wf%n_es_amplitudes)
         call mem%dealloc(p_and_q, this%wf%n_es_amplitudes)
!
         call mem%alloc(reduced_left, this%chain_length, this%chain_length)
         call mem%alloc(reduced_right, this%chain_length, this%chain_length)
         call mem%alloc(eigenvalues_Re, this%chain_length)
         call mem%alloc(eigenvalues_Im, this%chain_length)
!
         call lanczos%diagonalize_T(reduced_left, reduced_right, eigenvalues_Re, eigenvalues_Im)
!
         call lanczos%cleanup()
!
         call mem%alloc(oscillator_strength, this%chain_length)
!
         call this%calculate_oscillator_strength(reduced_left, reduced_right, eigenvalues_Re, &
                                                 eigenvalues_Im, norm_eta_times_norm_xi, &
                                                 oscillator_strength)
!
         call mem%dealloc(reduced_left, this%chain_length, this%chain_length)
         call mem%dealloc(reduced_right, this%chain_length, this%chain_length)
!
!        print results from Asymmetric Lanczos solver
!
         call this%print_summary(cartesian, eigenvalues_Re, eigenvalues_Im, &
                                    oscillator_strength)
!
         call mem%dealloc(eigenvalues_Re, this%chain_length)
         call mem%dealloc(eigenvalues_Im, this%chain_length)
         call mem%dealloc(oscillator_strength, this%chain_length)
!
      enddo
!
      call mem%dealloc(xiX,this%wf%n_es_amplitudes)
      call mem%dealloc(etaX,this%wf%n_es_amplitudes)
      call mem%dealloc(dipole_length,this%wf%n_mo,this%wf%n_mo,3)
!
   end subroutine run_asymmetric_lanczos_cc_es
!
!
   subroutine print_settings(this)
!!
!!    Print settings
!!    Written by Torsha Moitra, Nov 2019
!!    Added lanczos solver print by
!!
      implicit none
!
      class(asymmetric_lanczos_cc_es) :: this
!
      call output%printf('m', '- Settings for coupled cluster excited state &
                         &solver (' //trim(this%tag) // '):', fs='(/t3,a)')
!
      call output%printf('m', 'Chain length: (i6)', ints=[this%chain_length], &
                         fs='(/t6,a)')
!
      call output%printf('m', 'Biorthonormalization procedure: (a0)', &
                         chars=[trim(this%normalization)], fs='(t6,a)')
!
   end subroutine print_settings
!
!
   subroutine read_settings(this)
!!
!!    Read settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!    Adapted by Torsha Moitra for Asymmetric lanczos solver
!!
      use global_in, only: input
!
      implicit none
!
      class(asymmetric_lanczos_cc_es), intent(inout) :: this
!
      call input%get_keyword('chain length', 'solver cc es', this%chain_length)
!
      call input%get_keyword('lanczos normalization', 'solver cc es', &
                                                   this%normalization)
!
      if (input%is_keyword_present('core excitation', 'solver cc es')) &
            this%es_type = 'core'
!
      if (input%is_keyword_present('ionization', 'solver cc es')) &
            call output%error_msg('Ionization not possible with the asymmetric Lanczos solver')
!
      this%restart = input%is_keyword_present('restart', 'solver cc es') .or. &
                     input%is_keyword_present('restart', 'do')
!
   end subroutine read_settings
!
!
   subroutine initialize_energies(this)
!!
!!    Initialize energies
!!    Written by Torsha Moitra and Sarai D. Folkestad, Sep 2019
!!
      implicit none
!
      class(asymmetric_lanczos_cc_es) :: this
!
      call mem%alloc(this%energies, this%chain_length)
!
   end subroutine initialize_energies
!
!
   subroutine destruct_energies(this)
!!
!!    Destruct energies
!!    Written by Torsha Moitra and Sarai D. Folkestad, Sep 2019
!!
      implicit none
!
      class(asymmetric_lanczos_cc_es) :: this
!
      if (allocated(this%energies)) call mem%dealloc(this%energies, this%chain_length)
!
   end subroutine destruct_energies
!
!
   subroutine calculate_oscillator_strength(this, L, R, eigenvalues_Re, &
                                            eigenvalues_Im, norm_eta_times_norm_xi, &
                                            oscillator_strength)
!!
!!    Calculate oscillator strength
!!    Written by Sonia Coriani, Torsha Moitra and Sarai D. Folkestad, 2019
!!
!!    Constructs the oscillator strength
!!
!!       f(j)^q = (2/3) * omega(j) * norm_eta^q * norm_xi^q * L^T(j1) * R(1j)
!!
!!    for the cartesian component q
!!
!!    L and R are reduced space eigenvectors of the order of chain length.
!!
!
      use array_initialization, only: zero_array
!
      implicit none
!
      class(asymmetric_lanczos_cc_es), intent(in) :: this
!
      real(dp), dimension(this%chain_length, this%chain_length), intent(in) :: L
      real(dp), dimension(this%chain_length, this%chain_length), intent(in) :: R
!
      real(dp), dimension(this%chain_length), intent(in) :: eigenvalues_Re
      real(dp), dimension(this%chain_length), intent(in) :: eigenvalues_Im
!
      real(dp), intent(in) :: norm_eta_times_norm_xi
!
      real(dp), dimension(this%chain_length), intent(out) :: oscillator_strength

      integer :: j
!
      call zero_array(oscillator_strength, this%chain_length)
!
!$omp parallel do private(j)
      do j = 1, this%chain_length
!
!        If Eigenvalue is complex, skip oscillator strength
         if (abs(eigenvalues_Im(j)) .gt. 1.0d-12) cycle
!
         oscillator_strength(j) = L(1,j) * R(1,j)
!
      enddo
!$omp end parallel do
!
!$omp parallel do private(j)
      do j = 1, this%chain_length
!
         oscillator_strength(j) = eigenvalues_Re(j) * oscillator_strength(j)
!
      enddo
!$omp end parallel do
!
      call dscal(this%chain_length, (two/three)*norm_eta_times_norm_xi, oscillator_strength, 1)
!
   end subroutine calculate_oscillator_strength
!
!
   subroutine prepare_wf_for_excited_state(this)
!!
!!    Prepare wf for excited state
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2019
!!
      implicit none
!
      class(asymmetric_lanczos_cc_es) :: this
!
      call this%wf%prepare_for_jacobian()
      call this%wf%prepare_for_jacobian_transpose()
!
   end subroutine prepare_wf_for_excited_state
!
!
   subroutine print_summary(this, component, eigenvalues_Re, &
               eigenvalues_Im, oscillator_strength)
!!
!!    Print summary
!!    Written by Sarai D. Folkestad and Torsha Moitra, Nov 2019
!!
      use output_file_class, only: output_file
      use array_utilities,   only: quicksort_with_index_ascending
!
      implicit none
!
      class(asymmetric_lanczos_cc_es), intent(in) :: this
!
      integer, intent(in) :: component
!
      real(dp), dimension(this%chain_length), intent(inout)  :: eigenvalues_Re
      real(dp), dimension(this%chain_length), intent(in)     :: eigenvalues_Im
      real(dp), dimension(this%chain_length), intent(in)     :: oscillator_strength
!
      integer, dimension(:), allocatable :: index_list
!
      character(len=1), dimension(3) :: component_string = ['X', 'Y', 'Z']
!
      integer :: j
!
      integer, parameter :: print_in_output = 10
!
      character(len=6) :: chain_length_string
!
      type(output_file) :: spectrum_file
!
!     Order the real eigenvalues, we print from lowest to highest
!
      call mem%alloc(index_list, this%chain_length)
      call quicksort_with_index_ascending(eigenvalues_Re, index_list, &
                                          this%chain_length)
!
!     Initialize file for entire spectrum for the given component
!
      write(chain_length_string, '(i6)') this%chain_length
!
      spectrum_file = output_file('eT.lanczos' // trim(adjustl(chain_length_string)) // '_' &
                                   //component_string(component))
!
      call spectrum_file%open_()
!
!     Print entire spectrum
!
      call spectrum_file%printf('m', 'Nbr.      energy (Re) [eV]      energy (Im) [eV]    &
                                & energy (Re) [a.u]    energy (Im) [a.u.]   Osc. strength  ',   &
                                fs='(a)', ll=500)
!
      do j = 1, this%chain_length
!
         call spectrum_file%printf('m', '(i4) (f20.12)  (f20.12) (f20.12)' //    &
                                   ' (f20.12)     (f16.12)', ints=[j],           &
                                    reals=[eigenvalues_Re(j)*Hartree_to_eV,      &
                                    eigenvalues_Im(index_list(j))*Hartree_to_eV, &
                                    eigenvalues_Re(j),                           &
                                    eigenvalues_Im(index_list(j)),               &
                                    oscillator_strength(index_list(j))],         &
                                    fs='(a)', ll=500)
!
      enddo
!
      call spectrum_file%close_()
!
!     Print summary for eT.out
!
      if (component == 1) then
!
         call output%printf('m', '- Summary of the asymmetric Lanczos solver &
                            &for excited states', fs='(/t3,a)')
!
         call output%printf('m', 'Printing the (i0) lowest excited states for &
                            &each Cartesian component of the electric dipole moment', &
                            ints = [min(print_in_output, this%chain_length)], ffs='(/t6,a)', fs='(t6,a)')
!
      endif
!
      call output%printf('m', 'Component: (a1) ', &
                         chars=[component_string(component)], fs='(/t6,a)')
!
      call output%printf('m', 'State.      energy [a.u]         energy [eV]     &
                         &    Osc. strength  ', fs='(/t6, a)')
!
      call output%print_separator('minimal', 70, '-', fs='(t6, a)')
!

      do j = 1, min(print_in_output, this%chain_length)
!
         call output%printf('m', '(i4)  (f20.12) (f20.12)    (f16.12)', &
                            ints=[j], reals=[eigenvalues_Re(j), &
                            eigenvalues_Re(j)*Hartree_to_eV, &
                            oscillator_strength(index_list(j))], fs='(t6,a)')
!
      enddo
!
      call output%print_separator('minimal', 70, '-', fs='(t6, a)')
!
      call output%printf('m', 'For full spectrum see file: (a0)', &
                        chars=[spectrum_file%get_name()], fs='(t6,a)')
!
      call mem%dealloc(index_list, this%chain_length)
!
   end subroutine print_summary
!
!
   subroutine print_banner(this)
!!
!!    Print banner
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(asymmetric_lanczos_cc_es) :: this
!
      call output%printf('m', ' - ' // trim(this%name_), fs='(/t3,a)')
      call output%print_separator('m', len(trim(this%name_)) + 6, '-')
!
      call output%printf('n', this%description1, ffs='(/t3,a)', fs='(t3,a)')
      call output%printf('n', this%description2, ffs='(/t3,a)', fs='(t3,a)')
!
   end subroutine print_banner
!
!
   subroutine cleanup_asymmetric_lanczos_cc_es(this)
!!
!!    Cleanup
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      use string_utilities, only: convert_to_uppercase
!
      implicit none
!
      class(asymmetric_lanczos_cc_es), intent(inout) :: this
!
      call this%destruct_energies()
!
      call this%total_timer%turn_off()
!
      call output%printf('m', '- Finished solving the ' //  &
                         trim(convert_to_uppercase(this%wf%name_)) // ' excited &
                         &state equations.', fs='(/t3,a)')
!
      call output%printf('m', 'Total wall time (sec): (f20.5)', &
                         reals=[this%total_timer%get_elapsed_time('wall')], fs='(/t6,a)')
      call output%printf('m', 'Total cpu time (sec):  (f20.5)', &
                         reals=[this%total_timer%get_elapsed_time('cpu')], fs='(t6,a)')
!
   end subroutine cleanup_asymmetric_lanczos_cc_es
!
!
end module asymmetric_lanczos_cc_es_class
