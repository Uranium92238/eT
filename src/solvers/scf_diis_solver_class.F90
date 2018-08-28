module scf_diis_solver_class
!
!!
!!		Self-consistent field DIIS solver class module
!!		Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!!    This solver employs a direct inversion of the iterative
!!    subspace (DIIS) algorithm to get an effective density, based
!!    on previous AO densities and gradients, and from this density
!!    it diagonalizes the Roothan-Hall equation, giving a new density
!!    and gradient, Fock matrix and energy. The process is repeated
!!    until convergence.
!!
!
   use kinds
   use hf_solver_class
   use diis_tool_class
   use file_class
   use hf_class
   use disk_manager_class
   use libint_initialization
!
   implicit none
!
   type, extends(hf_solver) :: scf_diis_solver
!
      integer(i15) :: diis_dimension = 8
!
   contains
!
      procedure :: prepare   => prepare_scf_diis_solver
      procedure :: run       => run_scf_diis_solver
      procedure :: cleanup   => cleanup_scf_diis_solver
!
      procedure :: print_banner => print_banner_scf_diis_solver
!
      procedure :: read_settings          => read_settings_scf_diis_solver 
      procedure :: read_scf_diis_settings => read_scf_diis_settings_scf_diis_solver
!
   end type scf_diis_solver
!
!
contains
!
!
   subroutine prepare_scf_diis_solver(solver, wf)
!!
!!    Prepare 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(scf_diis_solver) :: solver
!
      class(hf) :: wf
!
!     Read settings (thresholds, etc.)
!
      call solver%read_settings()
!
      write(output%unit, *) 'Residual threshold: ', solver%residual_threshold
!
!     Set AO density to superposition of atomic densities (SAD)
!
      call wf%initialize_ao_density()
      call wf%set_ao_density_to_sad()
!
   end subroutine prepare_scf_diis_solver
!
!
   subroutine run_scf_diis_solver(solver, wf)
!!
!!    Run 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    This routine solves the Roothan-Hall equations in each macro-iteration,
!!    using the resulting density and Fock matrix together with previous densities
!!    and errors to make a new effective density matrix.
!!
      implicit none
!
      class(scf_diis_solver) :: solver
!
      class(hf) :: wf
!
      type(diis_tool) :: diis_manager
!
      logical :: converged
      logical :: converged_energy
      logical :: converged_residual
!
      real(dp) :: max_grad, energy, prev_energy
!
      real(dp) :: ddot
!
      integer(i15) :: iteration
!
      real(dp), dimension(:,:), allocatable :: D
      real(dp), dimension(:,:), allocatable :: F 
      real(dp), dimension(:,:), allocatable :: Po 
      real(dp), dimension(:,:), allocatable :: Pv 
      real(dp), dimension(:,:), allocatable :: G 
      real(dp), dimension(:,:), allocatable :: ao_fock 
      real(dp), dimension(:,:), allocatable :: h_wx 
      real(dp), dimension(:,:), allocatable :: prev_ao_density 
!
      integer(i15) :: n_s 
!
      real(dp), dimension(:,:), allocatable     :: sp_eri_schwarz
      integer(i15), dimension(:,:), allocatable :: sp_eri_schwarz_list
!
!     Print solver banner
!
      call solver%print_banner()
!
!     :: Construct screening vectors, as well as a degeneracy vector, 
!     used to construct AO Fock efficiently
!
      n_s = wf%system%get_n_shells()
!
      call mem%alloc(sp_eri_schwarz, n_s*(n_s + 1)/2, 2)
      call mem%alloc_int(sp_eri_schwarz_list, n_s*(n_s + 1)/2, 3)
!
      call wf%construct_sp_eri_schwarz(sp_eri_schwarz, sp_eri_schwarz_list, n_s)
!
!     :: Initialize the DIIS object,
!     which is used to get new effective densities from previous densities and gradients
!
      call diis_manager%init('hf_diis', (wf%n_ao)*(wf%n_ao + 1)/2, wf%n_ao**2, solver%diis_dimension)
!
!     :: Get an initial AO density and initial MO coefficients
!     by solving the Roothan-Hall equations for the SAD density
!
      call wf%initialize_ao_fock()
      call wf%construct_ao_fock_SAD(solver%coulomb_thr, solver%exchange_thr, solver%coulomb_precision)
!
      call wf%initialize_mo_coefficients()
!
      call wf%do_roothan_hall()
!
      call wf%construct_ao_density() 
      call wf%construct_ao_fock(sp_eri_schwarz, sp_eri_schwarz_list, n_s, &
                           solver%coulomb_thr, solver%exchange_thr, solver%coulomb_precision)  
!
      call mem%alloc(Po, wf%n_ao, wf%n_ao)
      call mem%alloc(Pv, wf%n_ao, wf%n_ao)
!
      call mem%alloc(ao_fock, wf%n_ao, wf%n_ao)          ! Holds Fock matrix, when wf%ao_fock 
                                                         ! becomes the DIIS Fock matrix 
!
      call mem%alloc(prev_ao_density, wf%n_ao, wf%n_ao)
!
      call mem%alloc(G, wf%n_ao, wf%n_ao)                ! Gradient 
      call mem%alloc(F, wf%n_ao*(wf%n_ao + 1)/2, 1)      ! Fock matrix packed 
!
      call wf%construct_projection_matrices(Po, Pv)
      call wf%construct_roothan_hall_gradient(G, Po, Pv)
!
      max_grad = get_abs_max(G, (wf%n_ao)**2)
!
      call packin(F, wf%ao_fock, wf%n_ao)
      call diis_manager%update(G, F)
!
      iteration = 1
      converged = .false.
!
      converged_energy   = .false.
      converged_residual = .false.
!
      prev_energy = zero
!
      call mem%alloc(h_wx, wf%n_ao, wf%n_ao)
      call get_ao_h_xy(h_wx)
!
      write(output%unit, '(t3,a)') 'Iteration    Energy (a.u.)        Max(grad.)    ΔE (a.u.)'
      write(output%unit, '(t3,a)') '----------------------------------------------------------'
!
      do while (.not. converged .and. iteration .le. solver%max_iterations)         
!
!        Set current energy
!
         energy = wf%hf_energy
!
!        Print current iteration information
!
         write(output%unit, '(t3,i3,10x,f17.12,4x,e10.4,4x,e10.4)') iteration, wf%hf_energy, &
                                          max_grad, abs(wf%hf_energy-prev_energy)
         flush(output%unit)
!
!        Test for convergence:
!
         converged_energy   = abs(energy-prev_energy) .lt. solver%energy_threshold
         converged_residual = max_grad                .lt. solver%residual_threshold
!
         converged = converged_residual .and. converged_energy
!
         if (converged) then
!
            write(output%unit, '(t3,a)') '---------------------------------------------------'
            write(output%unit, '(/t3,a13,i3,a12/)') 'Converged in ', iteration, ' iterations!'
!
         else
!
!           Solve the Roothan-Hall equation, then update the AO density 
!           and Fock matrix using the solution
!
            call wf%do_roothan_hall()
!
            prev_energy     = wf%hf_energy
            prev_ao_density = wf%ao_density
!
            if (iteration .ne. 1) wf%ao_fock = ao_fock 
!
            call wf%construct_ao_density()
            call wf%construct_ao_fock_cumulative(sp_eri_schwarz, sp_eri_schwarz_list,              &
                                                   n_s, prev_ao_density, h_wx, solver%coulomb_thr, &
                                                   solver%exchange_thr, solver%coulomb_precision)
!
!           Calculate the gradient from the current density
!
            call wf%construct_projection_matrices(Po, Pv)
            call wf%construct_roothan_hall_gradient(G, Po, Pv)
!
            max_grad = get_abs_max(G, (wf%n_ao)**2)
!
!           Update AO Fock matrix by DIIS,
!           but keep a copy of the actual Fock matrix for 
!           cumulative construction 
!
            call packin(F, wf%ao_fock, wf%n_ao)
            call diis_manager%update(G, F)
!
            ao_fock = wf%ao_fock  
            call wf%set_ao_fock(F)
!
         endif
!
         iteration = iteration + 1
!
      enddo
!
      call mem%dealloc(sp_eri_schwarz, n_s*(n_s + 1)/2, 2)
      call mem%alloc_int(sp_eri_schwarz_list, n_s*(n_s + 1)/2, 3)
!
      call mem%dealloc(Po, wf%n_ao, wf%n_ao)
      call mem%dealloc(Pv, wf%n_ao, wf%n_ao)
!
      call mem%dealloc(G, wf%n_ao, wf%n_ao)           
      call mem%dealloc(F, wf%n_ao*(wf%n_ao + 1)/2, 1) 
!
      call mem%dealloc(h_wx, wf%n_ao, wf%n_ao)
!
      call wf%destruct_ao_overlap()
!
!     Initialize engine (make final deallocations, and other stuff)
!
      call diis_manager%finalize()
!
      if (.not. converged) then 
!
         write(output%unit, '(t3,a)')   '---------------------------------------------------'
         write(output%unit, '(/t3,a)')  'Was not able to converge the equations in the given'
         write(output%unit, '(t3,a/)')  'number of maximum iterations.'
         stop
!
      endif 
!
   end subroutine run_scf_diis_solver
!
!
   subroutine cleanup_scf_diis_solver(solver, wf)
!!
!! 	Cleanup 
!! 	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(scf_diis_solver) :: solver
!
      class(hf) :: wf
!
      logical :: do_mo_transformation
!
      integer(i15) :: i
!
!     When finished, we do a final Roothan-Hall step,
!     to express the Fock matrix in the canonical MO basis 
!
      do_mo_transformation = .true.
      call wf%do_roothan_hall(do_mo_transformation)
!
   end subroutine cleanup_scf_diis_solver
!
!
   subroutine print_banner_scf_diis_solver(solver)
!!
!!    Print banner
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none 
!
      class(scf_diis_solver) :: solver 
!
      write(output%unit, '(/t3,a)') ':: Direct-integral self-consistent field DIIS Hartree-Fock solver'
      write(output%unit, '(t3,a/)') ':: E. F. Kjønstad, S. D. Folkestad, 2018'
!
      write(output%unit, '(t3,a)')  'This solver employs a direct inversion of the iterative'
      write(output%unit, '(t3,a)')  'subspace (DIIS) algorithm to get an effective density, based'
      write(output%unit, '(t3,a)')  'on previous AO densities and gradients, and from this density'
      write(output%unit, '(t3,a)')  'it diagonalizes the Roothan-Hall equation, giving a new density'
      write(output%unit, '(t3,a)')  'and gradient, Fock matrix and energy. The process is repeated'
      write(output%unit, '(t3,a/)') 'until convergence.'
      flush(output%unit)
!
   end subroutine print_banner_scf_diis_solver
!
!
   subroutine read_settings_scf_diis_solver(solver)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(scf_diis_solver) :: solver 
!
      call solver%read_hf_solver_settings()
      call solver%read_scf_diis_settings()
!
   end subroutine read_settings_scf_diis_solver
!
!
   subroutine read_scf_diis_settings_scf_diis_solver(solver)
!!
!!    Read SCF DIIS settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
!!    Reads settings specific to the class. 
!!
      implicit none 
!
      class(scf_diis_solver) :: solver 
!
      integer(i15) :: n_records, i 
!
      character(len=100) :: line, value 
!
      if (requested_section('hf')) then ! User has requested something 
!
         call move_to_section('hf', n_records)
!
         do i = 1, n_records
!
            read(input%unit, '(a100)') line
            line = remove_preceding_blanks(line)
!
            write(output%unit, *) trim(line)
!
            if (line(1:15) == 'diis_dimension:') then
!
               value = line(16:100)
               value = remove_preceding_blanks(value)
               read(value, *) solver%diis_dimension
               return
!
            endif
!
         enddo
!
      endif 
!
   end subroutine read_scf_diis_settings_scf_diis_solver
!
!
end module scf_diis_solver_class
