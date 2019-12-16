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
!!       f(j)^q = (2/3) * omega(j) norm_eta^q norm_csi^q L^T_n(j1) R_n(1j)
!!
!!    See for further details: 
!!       S. Coriani et. al., J. Chem. Theory Comput. 2012, 8, 5, 1616-1628
!! 
   use parameters
!
   use ccs_class,             only : ccs
   use timings_class,         only : timings 
   use global_out,            only : output
   use global_in,             only : input
   use memory_manager_class,  only : mem
   use string_utilities,      only : convert_to_uppercase
!
   use es_projection_tool_class,          only : es_projection_tool
   use es_cvs_projection_tool_class,      only : es_cvs_projection_tool
   use es_valence_projection_tool_class,  only : es_valence_projection_tool

   implicit none

   type :: asymmetric_lanczos_cc_es

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
      type(timings) :: timer
!
      class(es_projection_tool), allocatable :: projector
!
   contains
!
      procedure :: read_settings                   =>  read_settings_asymmetric_lanczos_cc_es
      procedure :: print_settings                  =>  print_settings_asymmetric_lanczos_cc_es
!
      procedure :: initialize_energies             => initialize_energies_asymmetric_lanczos_cc_es
      procedure :: destruct_energies               => destruct_energies_asymmetric_lanczos_cc_es
!
      procedure :: calculate_oscillator_strength   &
                  => calculate_oscillator_strength_asymmetric_lanczos_cc_es
!
      procedure :: run                             => run_asymmetric_lanczos_cc_es
      procedure :: cleanup                         => cleanup_asymmetric_lanczos_cc_es
!
      procedure, nopass :: prepare_wf_for_excited_state  &
                           => prepare_wf_for_excited_state_asymmetric_lanczos_cc_es
!
      procedure :: print_summary => print_summary_asymmetric_lanczos_cc_es
      procedure :: print_banner  => print_banner_asymmetric_lanczos_cc_es
!
   end type asymmetric_lanczos_cc_es
!
!
   interface asymmetric_lanczos_cc_es
!           
      procedure :: new_asymmetric_lanczos_cc_es
!
   end interface asymmetric_lanczos_cc_es
!
!
contains
!
!   
   function new_asymmetric_lanczos_cc_es(wf) result(solver)
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
!
      use array_utilities, only: zero_array  
      implicit none
!
      type(asymmetric_lanczos_cc_es) :: solver
!
      class(ccs), intent(inout) :: wf
!
      solver%timer = timings(trim(convert_to_uppercase(wf%name_)) // ' excited state') 
      call solver%timer%turn_on()
!
!     Set printables
!
      solver%name_  = 'Asymmetric Lanczos excited state solver'
      solver%tag    = 'asymmetric lanczos'
!
      solver%description1 = 'An asymmetric Lanczos solver builds a reduced space&
               & tridiagonal representation of the CC Jacobian of dimension defined&
               & by the chain length.&
               & Diagonalization of this representation gives the eigenvalues and eigenvectors.'
!
      solver%description2 = 'A complete description of the algorithm can be found in &
               & S.Coriani et al., J. Chem. Theory Comput. 2012, 8, 5, 1616-1628.'
!
      call solver%print_banner()
!
!     Set defaults
!
      solver%normalization = 'asymmetric'
      solver%es_type       = 'valence'
!   
      call solver%read_settings()
!
!     Chain length cannot exceed the number of excited state amplitudes
!     (actual dimension of A in the basis)
!
      if(solver%chain_length .gt. wf%n_es_amplitudes) then
!
          solver%chain_length = wf%n_es_amplitudes
!
          call output%warning_msg('Chain length asked for is reset to: (i0)', &
              & ints=[solver%chain_length])
!
      endif
!
      call solver%print_settings()
!
      call solver%initialize_energies()
!
!     Prepare wavefunction for Jacobian transformations
!
      call solver%prepare_wf_for_excited_state(wf)
!
      if (trim(solver%es_type) == 'valence') then
!
         solver%projector = es_valence_projection_tool()
!
      elseif (trim(solver%es_type) == 'core') then
!
         call wf%read_cvs_settings()
         solver%projector = es_cvs_projection_tool(wf)
!
      else
!
         call output%error_msg('did not recognize the excited state type: ' //trim(solver%es_type))
!
      endif
!
   end function new_asymmetric_lanczos_cc_es
!
!
   subroutine run_asymmetric_lanczos_cc_es (solver,wf)
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
!
      use asymmetric_lanczos_tool_class,  only : asymmetric_lanczos_tool
!
      implicit none

      class(asymmetric_lanczos_cc_es) :: solver
!
      class(ccs) :: wf
!
      type(asymmetric_lanczos_tool), allocatable :: lanczos  
!
      real(dp), dimension(:,:,:), allocatable   :: dipole_length
      real(dp), dimension(:,:), allocatable     :: reduced_left, reduced_right
      real(dp), dimension(:), allocatable       :: etaX, csiX
      real(dp), dimension(:), allocatable       :: Aq, pA
      real(dp), dimension(:), allocatable       :: eigenvalues_Re, eigenvalues_Im
      real(dp), dimension(:), allocatable       :: oscillator_strength
!
      real(dp) :: norm_eta_times_norm_csi, ddot, overlap, factor, norm_p1q1
!
      integer :: cartesian, iteration
!
!     Construct the dipole integrals
!
      call mem%alloc(dipole_length, wf%n_mo, wf%n_mo, 3)
      call wf%construct_mu(dipole_length)
!
      call mem%alloc(csiX, wf%n_es_amplitudes)     
      call mem%alloc(etaX, wf%n_es_amplitudes)     
!
      do cartesian=1, 3
!
!        Construct the eta and csi vectors
!
         call wf%construct_csiX(dipole_length(:,:,cartesian), csiX)
         call wf%construct_eom_etaX(dipole_length(:,:,cartesian), csiX, etaX)
!
         if (solver%projector%active) then
!
            call solver%projector%do_(csiX)
            call solver%projector%do_(etaX)
!
         endif
!
!        Prepare initial seeds (p1,q1) as binormalized etaX and csiX vectors
!        
!        The binormalization procedure is either 'symmetric' or 'asymmetric'
!
         overlap = ddot(wf%n_es_amplitudes, etaX, 1, csiX, 1)
!
         if(solver%normalization=="symmetric")then
!
            norm_p1q1 = dsqrt(dabs(overlap))
!
            factor = one/norm_p1q1
!
            call dscal(wf%n_es_amplitudes, factor, csiX, 1)
            call dscal(wf%n_es_amplitudes, dsign(factor,overlap), etaX, 1)
!         
         else if (solver%normalization=="asymmetric") then
!
            factor = one/overlap
!
            call dscal(wf%n_es_amplitudes, factor, etaX, 1)
!
         else 
!
            call output%error_msg('did not recognize the binormalization procedure')
!
         endif
!
!        Product of the norm of eta and the norm of csi, to be used for oscillator strengths
!
         norm_eta_times_norm_csi = overlap
!
!        Initialize the asymmetric Lanczos tool
!
         lanczos = asymmetric_lanczos_tool(wf%n_es_amplitudes, solver%chain_length, &
                                             etaX, csiX, solver%normalization)
!
         call mem%alloc(Aq, wf%n_es_amplitudes)
         call mem%alloc(pA, wf%n_es_amplitudes)
!
!        Build the Krylov subspace and the tridiagonal matrix T
!        Loops over the chain-length 
!
         do iteration = 1, solver%chain_length 
!
!           p_(i) and q_(i) are stored in pA and Aq
!             
            call lanczos%read_p(pA, iteration)
            call lanczos%read_q(Aq, iteration)
!
!           Transforms p -> pA and q -> Aq
!
            call wf%jacobian_transformation(Aq) 
            call wf%jacobian_transpose_transformation(pA)
!
            if (solver%projector%active) then
!
               call solver%projector%do_(Aq)
               call solver%projector%do_(pA)
!
            endif
!            
!           Calculates alpha_(i) 
! 
            call lanczos%calculate_alpha(Aq, iteration)
!
!           Calculates beta_(i), gamme_(i), p_(i+1) and q_(i+1)
!
            if (iteration .lt. solver%chain_length) then
!
              call lanczos%calculate_beta_gamma_p_q(iteration, Aq, pA)  
!
            end if
!
!           Checks if the value of beta is too small then stops iteration
!           Resets chain length to current iteration. 
!
            if (lanczos%chain_length .lt. solver%chain_length) then 
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
         solver%chain_length = lanczos%chain_length  
!
         call mem%dealloc (Aq, wf%n_es_amplitudes)
         call mem%dealloc (pA, wf%n_es_amplitudes)
!
!        Diagonalize T to get the eigenvalues and eigenvectors 
!
         call mem%alloc(reduced_left, solver%chain_length, solver%chain_length)
         call mem%alloc(reduced_right, solver%chain_length, solver%chain_length)
!
         call mem%alloc(eigenvalues_Re, solver%chain_length)
         call mem%alloc(eigenvalues_Im, solver%chain_length)
!
         call lanczos%diagonalize_T(reduced_left, reduced_right, eigenvalues_Re, eigenvalues_Im)
!
!        Done with the tool, clean it up
!
         call lanczos%cleanup()
!
         call mem%alloc(oscillator_strength, solver%chain_length) 
!
!        Calculate oscillator strength
!
         call solver%calculate_oscillator_strength(reduced_left, reduced_right, eigenvalues_Re, &
                                                   norm_eta_times_norm_csi, oscillator_strength) 
!
         call mem%dealloc(reduced_left, solver%chain_length, solver%chain_length)
         call mem%dealloc(reduced_right, solver%chain_length, solver%chain_length)
!
!        print results from Asymmetric Lanczos solver
!
         call solver%print_summary(wf, cartesian, eigenvalues_Re, eigenvalues_Im, &
                                    oscillator_strength)
!
         call mem%dealloc(eigenvalues_Re, solver%chain_length)
         call mem%dealloc(eigenvalues_Im, solver%chain_length)
         call mem%dealloc(oscillator_strength, solver%chain_length)     
!
      enddo
!
      call mem%dealloc(csiX,wf%n_es_amplitudes)
      call mem%dealloc(etaX,wf%n_es_amplitudes)
      call mem%dealloc(dipole_length,wf%n_mo,wf%n_mo,3)
!
   end subroutine run_asymmetric_lanczos_cc_es
!
!
   subroutine print_settings_asymmetric_lanczos_cc_es(solver)
!!
!!    Print settings 
!!    Written by Torsha Moitra, Nov 2019
!!    Added lanczos solver print by 
!!
      implicit none 
!
      class(asymmetric_lanczos_cc_es) :: solver 
!
      call output%printf('m', '- Settings for coupled cluster excited state &
                         &solver (' //trim(solver%tag) // '):', fs='(/t3,a)')
!
      call output%printf('m', 'Chain length: (i6)', ints=[solver%chain_length], &
                         fs='(/t6,a)')
!
     call output%printf('m', 'Biorthonormalization procedure: (a0)', &
                        chars=[trim(solver%normalization)], fs='(t6,a)')
!
   end subroutine print_settings_asymmetric_lanczos_cc_es
!
!
   subroutine read_settings_asymmetric_lanczos_cc_es(solver)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!    Adapted by Torsha Moitra for Asymmetric lanczos solver
!!
      implicit none 
!
      class(asymmetric_lanczos_cc_es) :: solver 
!
      call input%get_keyword_in_section('chain length', 'solver cc es', solver%chain_length)
!
      call input%get_keyword_in_section('lanczos normalization', 'solver cc es', &
                                                   solver%normalization)
!
      if (input%requested_keyword_in_section('core excitation', 'solver cc es')) &
            solver%es_type = 'core'
!
      if (input%requested_keyword_in_section('ionization', 'solver cc es')) &
            call output%error_msg('Ionization not possible with the asymmetric Lanczos solver')
!
   end subroutine read_settings_asymmetric_lanczos_cc_es
!
!
   subroutine initialize_energies_asymmetric_lanczos_cc_es(solver)
!!
!!    Initialise energies
!!    Written by Torsha Moitra and Sarai D. Folkestad, Sep 2019
!!
      implicit none
!
      class(asymmetric_lanczos_cc_es) :: solver
!
      call mem%alloc(solver%energies, solver%chain_length)
!
   end subroutine initialize_energies_asymmetric_lanczos_cc_es
!
!
   subroutine destruct_energies_asymmetric_lanczos_cc_es(solver)
!!
!!    Destruct energies
!!    Written by Torsha Moitra and Sarai D. Folkestad, Sep 2019
!!
      implicit none
!
      class(asymmetric_lanczos_cc_es) :: solver
!
      if (allocated(solver%energies)) call mem%dealloc(solver%energies, solver%chain_length)
!
   end subroutine destruct_energies_asymmetric_lanczos_cc_es
!
!
   subroutine calculate_oscillator_strength_asymmetric_lanczos_cc_es(solver, L, R, eigenvalues, &
            norm_eta_times_norm_csi, oscillator_strength)
!!
!!    Calculate oscillator strength
!!    Written by Sonia Coriani, Torsha Moitra and Sarai D. Folkestad, 2019
!!
!!    Constructs the oscillator strength
!!
!!       f(j)^q = (2/3) * omega(j) * norm_eta^q * norm_csi^q * L^T(j1) * R(1j)
!!
!!    for the cartesian component q
!!
!!    L and R are reduced space eigenvectors of the order of chain length.
!!
      implicit none
!
      class(asymmetric_lanczos_cc_es), intent(in) :: solver
!      
      real(dp), dimension(solver%chain_length, solver%chain_length), intent(in) :: L
      real(dp), dimension(solver%chain_length, solver%chain_length), intent(in) :: R
!
      real(dp), dimension(solver%chain_length), intent(in) :: eigenvalues
!
      real(dp), intent(in) :: norm_eta_times_norm_csi
!
      real(dp), dimension(solver%chain_length), intent(out) :: oscillator_strength

      integer :: j
!
!$omp parallel do private(j)
      do j = 1, solver%chain_length
!
         oscillator_strength(j) = L(1,j) * R(1,j)
!
      enddo
!$omp end parallel do
!
!$omp parallel do private(j)
      do j = 1, solver%chain_length
!
         oscillator_strength(j) = eigenvalues(j) * oscillator_strength(j)
!
      enddo
!$omp end parallel do
!
      call dscal(solver%chain_length, (two/three)*norm_eta_times_norm_csi, oscillator_strength, 1)
!
   end subroutine calculate_oscillator_strength_asymmetric_lanczos_cc_es
!
!
   subroutine prepare_wf_for_excited_state_asymmetric_lanczos_cc_es(wf)
!!
!!    Prepare wf for excited state
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2019
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      call wf%prepare_for_jacobian()
      call wf%prepare_for_jacobian_transpose()
!
   end subroutine prepare_wf_for_excited_state_asymmetric_lanczos_cc_es
!
!
   subroutine print_summary_asymmetric_lanczos_cc_es(solver, wf, component, eigenvalues_Re, &
               eigenvalues_Im, oscillator_strength)
!!
!!    Print summary 
!!    Written by Sarai D. Folkestad and Torsha Moitra, Nov 2019
!!
!
      use output_file_class,              only: output_file
      use array_utilities,                only: quicksort_with_index_ascending
!
      implicit none
!
      class(asymmetric_lanczos_cc_es), intent(in) :: solver
!
      class(ccs), intent(in) :: wf
!
      integer, intent(in) :: component
!
      real(dp), dimension(solver%chain_length), intent(inout)  :: eigenvalues_Re
      real(dp), dimension(solver%chain_length), intent(in)     :: eigenvalues_Im
      real(dp), dimension(solver%chain_length), intent(in)     :: oscillator_strength
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
      call mem%alloc(index_list, solver%chain_length)
      call quicksort_with_index_ascending(eigenvalues_Re, index_list, &
                                          solver%chain_length)
!
!     Initialize file for entire spectrum for the given component
!
      write(chain_length_string, '(i6)') solver%chain_length
!
      spectrum_file = output_file(trim(wf%system%name_)//'_'// &
               trim(adjustl(chain_length_string))//'_'//component_string(component))
!
      call spectrum_file%open_
!
!     Print entire spectrum
!
      call spectrum_file%printf('m', 'Nbr.      energy (Re) [eV]      energy (Im) [eV]    &
                                & energy (Re) [a.u]    energy (Im) [a.u.]   Osc. strength  ',   &
                                fs='(a)', ll=500)
!
      do j = 1, solver%chain_length
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
!     Print summary for eT.out
!
      if (component == 1) then
!
         call output%printf('m', '- Summary of the asymmetric Lanczos solver &
                            &for excited states', fs='(/t3,a)')
!
         call output%printf('m', 'Printing the (i0) lowest excited states for &
                            &each Cartesian component of the electric dipole moment', &
                            ints = [min(print_in_output, solver%chain_length)], ffs='(/t6,a)', fs='(t6,a)')
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

      do j = 1, min(print_in_output, solver%chain_length)
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
      call output%printf('m', 'For full spectrum see file: ' //  &
                        trim(wf%system%name_) // '_' //  &
                        trim(adjustl(chain_length_string)) // '_' //  &
                        component_string(component), fs='(t6,a)')
!
      call mem%dealloc(index_list, solver%chain_length)
!
   end subroutine print_summary_asymmetric_lanczos_cc_es
!
!
   subroutine print_banner_asymmetric_lanczos_cc_es(solver)
!!
!!    Print banner
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none 
!
      class(asymmetric_lanczos_cc_es) :: solver 
!
      call output%printf('m', ' - ' // trim(solver%name_), fs='(/t3,a)')
      call output%print_separator('m', len(trim(solver%name_)) + 6, '-')
!
      call output%printf('n', solver%description1, ffs='(/t3,a)', fs='(t3,a)')
      call output%printf('n', solver%description2, ffs='(/t3,a)', fs='(t3,a)')
!
   end subroutine print_banner_asymmetric_lanczos_cc_es
!
!
   subroutine cleanup_asymmetric_lanczos_cc_es(solver, wf)
!!
!!    Cleanup 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(asymmetric_lanczos_cc_es)  :: solver
      class(ccs), intent(in)           :: wf
!
      call solver%destruct_energies()
!
      call solver%timer%turn_off()
!
      call output%printf('m', '- Finished solving the ' //  &
                         trim(convert_to_uppercase(wf%name_)) // ' excited &
                         &state equations.', fs='(/t3,a)')
!
      call output%printf('m', 'Total wall time (sec): (f20.5)', &
                         reals=[solver%timer%get_elapsed_time('wall')], fs='(/t6,a)')
      call output%printf('m', 'Total cpu time (sec):  (f20.5)', &
                         reals=[solver%timer%get_elapsed_time('cpu')], fs='(t6,a)')
!
   end subroutine cleanup_asymmetric_lanczos_cc_es
!
!
end module asymmetric_lanczos_cc_es_class
