module hf_class
!
!!
!!    Hartree-Fock (HF) class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   use kinds
   use file_class
   use atom_class
   use reordering
   use index
   use integrals_class
   use molecule_class
   use disk_manager_class
!
   implicit none
!
   type :: hf
!
      character(len=40) :: name = 'HF'
!
      integer(i15) :: n_ao
      integer(i15) :: n_mo
!
      integer(i15) :: n_o
      integer(i15) :: n_v
!
      real(dp) :: hf_energy
!
      real(dp), dimension(:,:), allocatable :: ao_density
      real(dp), dimension(:,:), allocatable :: ao_overlap
      real(dp), dimension(:,:), allocatable :: ao_fock
!
      real(dp), dimension(:,:), allocatable :: mo_coefficients
      real(dp), dimension(:,:), allocatable :: orbital_energies
!
      type(molecule)  :: molecule
      type(integrals) :: integrals
!
	contains
!
      procedure :: initialize => initialize_hf
      procedure :: finalize   => finalize_hf
!
      procedure :: construct_ao_density => construct_ao_density_hf
      procedure :: construct_ao_fock    => construct_ao_fock_hf
!
      procedure :: calculate_hf_energy  => calculate_hf_energy_hf
!
      procedure :: get_n_hf_parameters => get_n_hf_parameters_hf
      procedure :: get_n_hf_equations  => get_n_hf_equations_hf
!
      procedure :: set_initial_hf_parameters => set_initial_hf_parameters_hf
      procedure :: calculate_hf_equations    => calculate_hf_equations_hf
!
      procedure :: get_ao_density_from_mo_coefficients => get_ao_density_from_mo_coefficients_hf
!
   end type hf
!
!
contains
!
!
   subroutine initialize_hf(wf)
!!
!! 	Initialize
!!  	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp) :: repulsion
!
      real(dp), dimension(:,:), allocatable :: h
!
      integer(i15) :: i = 0
!
!     Initialize molecule
!
      call wf%molecule%initialize ! Read atoms
      call wf%molecule%write      ! Write an xyz file for the read geometry
!
!     Determine the number of atomic orbitals and (standard) number
!     of molecular orbitals
!
      wf%n_ao = 0
      call get_n_aos(wf%n_ao)
!
      wf%n_mo = wf%n_ao
!
!     Determine the number of occupied and virtual molecular orbitals
!
      wf%n_o = (wf%molecule%get_n_electrons())/2
      wf%n_v = wf%n_mo - wf%n_o
!
!     Allocate member arrays
!
      call mem%alloc(wf%mo_coefficients, wf%n_ao, wf%n_ao)
      call mem%alloc(wf%ao_density, wf%n_ao, wf%n_ao)
      call mem%alloc(wf%ao_fock, wf%n_ao, wf%n_ao)
      call mem%alloc(wf%ao_overlap, wf%n_ao, wf%n_ao)
!
      wf%mo_coefficients = zero
      wf%ao_density = zero
      wf%ao_fock = zero
      wf%ao_overlap = zero
!
   end subroutine initialize_hf
!
!
   subroutine finalize_hf(wf)
!!
!!    Finalize
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
   end subroutine finalize_hf
!
!
   subroutine construct_ao_density_hf(wf)
!!
!!    Construct AO density
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Routine which calculates D_αβ = sum_i C_αi C_βi,
!!    where C are the MO coefficients.
!!
      implicit none
!
      class(hf) :: wf
!
      call dgemm('N', 'T',            &
                  wf%n_ao,            &
                  wf%n_ao,            &
                  wf%n_o,             &
                  two,                &
                  wf%mo_coefficients, &
                  wf%n_ao,            &
                  wf%mo_coefficients, &
                  wf%n_ao,            &
                  zero,               &
                  wf%ao_density,      &
                  wf%n_ao)
!
   end subroutine construct_ao_density_hf
!
!
   subroutine construct_ao_fock_hf(wf)
!!
!!    Construct AO Fock matrix
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Routine which calculates
!!
!!       F_αβ = h_αβ + sum_γδ g_αβγδ D_γδ - 1/2 * sum_γδ g_αδγβ D_γδ,
!!
!!    where D is the AO density.
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(:,:), allocatable :: g_xyzw
      real(dp), dimension(:,:), allocatable :: g_xwzy
!
!     F_αβ = h_αβ
!
      wf%ao_fock = zero
      call get_ao_h_xy(wf%ao_fock)
!
!     F_αβ =+ sum_γδ g_αβγδ D_γδ
!
      call mem%alloc(g_xyzw, (wf%n_ao)**2, (wf%n_ao)**2)
      g_xyzw = zero
!
      call get_ao_g_xyzw(g_xyzw)
!
      call dgemm('N', 'N',       &
                  (wf%n_ao)**2,  &
                  1,             &
                  (wf%n_ao)**2,  &
                  one,           &
                  g_xyzw,        & ! g_αβ_γδ
                  (wf%n_ao)**2,  &
                  wf%ao_density, & ! D_γδ
                  (wf%n_ao)**2,  &
                  one,           &
                  wf%ao_fock,    & ! F_αβ
                  (wf%n_ao)**2)
!
      call mem%alloc(g_xwzy, (wf%n_ao)**2, (wf%n_ao)**2)
      g_xwzy = zero
!
!     F_αβ =+ (-1/2)*sum_γδ g_αδγβ D_γδ
!
      call sort_1234_to_1432(g_xyzw, g_xwzy, wf%n_ao, wf%n_ao, wf%n_ao, wf%n_ao)
!
      call mem%dealloc(g_xyzw, (wf%n_ao)**2, (wf%n_ao)**2)
!
      call dgemm('N', 'N',       &
                  (wf%n_ao)**2,  &
                  1,             &
                  (wf%n_ao)**2,  &
                  -half,         &
                  g_xwzy,        & ! g_αβ_γδ = g_αδγβ
                  (wf%n_ao)**2,  &
                  wf%ao_density, & ! D_γδ
                  (wf%n_ao)**2,  &
                  one,           &
                  wf%ao_fock,    & ! F_αβ
                  (wf%n_ao)**2)
!
      call mem%dealloc(g_xwzy, (wf%n_ao)**2, (wf%n_ao)**2)
!
   end subroutine construct_ao_fock_hf
!
!
   subroutine calculate_hf_energy_hf(wf)
!!
!!    Calculate HF energy
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Calculates the Hartree-Fock energy,
!!
!!       E = Tr(h D) + 1/4 * Tr(D G(D)),
!!
!!    where D is the AO density and
!!
!!       G(D)_αβ = 2 sum_γδ g_αβγδ D_γδ - sum_γδ g_αδγβ D_γδ
!!
!!    The traces are calculated as dot products (since B is symmetric here):
!!
!!       Tr(AB) = sum_x (AB)_xx = sum_xy A_xy B_yx = sum_xy A_xy B_xy.
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp) :: ddot
!
      real(dp), dimension(:,:), allocatable :: GD_xy
!
      real(dp), dimension(:,:), allocatable :: h_xy
      real(dp), dimension(:,:), allocatable :: g_xyzw
      real(dp), dimension(:,:), allocatable :: g_xwzy
!
!     Construct G(D)
!
      call mem%alloc(g_xyzw, (wf%n_ao)**2, (wf%n_ao)**2)
      g_xyzw = zero
!
      call get_ao_g_xyzw(g_xyzw)
!
!     G(D)_αβ = 2 sum_γδ g_αβγδ D_γδ
!
      call mem%alloc(GD_xy, wf%n_ao, wf%n_ao)
!
      call dgemm('N', 'N',       &
                  (wf%n_ao)**2,  &
                  1,             &
                  (wf%n_ao)**2,  &
                  two,           &
                  g_xyzw,        & ! g_αβ_γδ
                  (wf%n_ao)**2,  &
                  wf%ao_density, & ! D_γδ
                  (wf%n_ao)**2,  &
                  zero,          &
                  GD_xy,         & ! G(D)_αβ
                  (wf%n_ao)**2)
!
!     G(D)_αβ =+ (-1) sum_γδ g_αδγβ D_γδ
!
      call mem%alloc(g_xwzy, (wf%n_ao)**2, (wf%n_ao)**2)
      g_xwzy = zero
!
      call sort_1234_to_1432(g_xyzw, g_xwzy, wf%n_ao, wf%n_ao, wf%n_ao, wf%n_ao)
!
      call mem%dealloc(g_xyzw, (wf%n_ao)**2, (wf%n_ao)**2)
!
      call dgemm('N', 'N',       &
                  (wf%n_ao)**2,  &
                  1,             &
                  (wf%n_ao)**2,  &
                  -one,          &
                  g_xwzy,        & ! g_αβ_γδ = g_αδγβ
                  (wf%n_ao)**2,  &
                  wf%ao_density, & ! D_γδ
                  (wf%n_ao)**2,  &
                  one,           &
                  GD_xy,         & ! G(D)_αβ
                  (wf%n_ao)**2)
!
      call mem%dealloc(g_xwzy, (wf%n_ao)**2, (wf%n_ao)**2)
!
      call mem%alloc(h_xy, wf%n_ao, wf%n_ao)
      h_xy = zero
!
      call get_ao_h_xy(h_xy)
!
      wf%hf_energy = wf%molecule%get_nuclear_repulsion()
!
      wf%hf_energy = wf%hf_energy + ddot((wf%n_ao)**2, h_xy, 1, wf%ao_density, 1)
      wf%hf_energy = wf%hf_energy + (one/four)*ddot((wf%n_ao)**2, wf%ao_density, 1, GD_xy, 1)
!
      call mem%dealloc(h_xy, wf%n_ao, wf%n_ao)
      call mem%dealloc(GD_xy, wf%n_ao, wf%n_ao)
!
   end subroutine calculate_hf_energy_hf
!
!
   function get_n_hf_parameters_hf(wf)
!!
!!    Get number of HF parameters
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    When solving the equations, we change the AO density matrix,
!!    which we can change in packed form, giving n_ao*(n_ao + 1)/2
!!    number of parameters.
!!
      implicit none
!
      class(hf) :: wf
!
      integer(i15) :: get_n_hf_parameters_hf
!
      get_n_hf_parameters_hf = (wf%n_ao)*(wf%n_ao + 1)/2
!
   end function get_n_hf_parameters_hf
!
!
   function get_n_hf_equations_hf(wf)
!!
!!    Get number of HF equations
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    The equations to solve are F_ia = 0, where F is the MO Fock
!!    matrix. The number of equations is n_o*n_v.
!!
      implicit none
!
      class(hf) :: wf
!
      integer(i15) :: get_n_hf_equations_hf
!
      get_n_hf_equations_hf = (wf%n_o)*(wf%n_v)
!
   end function get_n_hf_equations_hf
!
!
   subroutine set_initial_hf_parameters_hf(wf, D_xy)
!!
!!    Set initial HF parameters
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Sets the initial AO density by using the superposition of
!!    atomic densities (SAD) guess
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(:,:) :: D_xy ! AO density, D_αβ, packed
!
      D_xy(index_packed(1,1), 1) = one
      D_xy(index_packed(6,6), 1) = one
!
   end subroutine set_initial_hf_parameters_hf
!
!
   subroutine calculate_hf_equations_hf(wf, F_ia, D_xy)
!!
!!    Calculate HF equations
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(:,:) :: D_xy ! AO density, D_αβ, packed
      real(dp), dimension(wf%n_o, wf%n_v) :: F_ia ! Occ-vir block of MO Fock matrix
!
      real(dp), dimension(:,:), allocatable :: work
      real(dp), dimension(:,:), allocatable :: ao_fock_copy
!
      real(dp), dimension(:,:), allocatable :: X_mat, Y_mat
!
      real(dp) :: ddot, norm
!
      integer(i15) :: info = 0, i = 0, a = 0, ia = 0
!
!     Square up the AO density and save it
!
      wf%ao_density = zero
      call squareup(D_xy, wf%ao_density, wf%n_ao)
!
!     Calculate the energy from the AO density
!
      call wf%calculate_hf_energy()
!
      write(output%unit, *) 'Energy: ', wf%hf_energy
      flush(output%unit)
!
!     Calculate the AO Fock matrix from the AO density
!
      wf%ao_fock = zero
      call wf%construct_ao_fock()
!
!     Solve the Roothan-Hall equations
!
      call mem%alloc(wf%orbital_energies, wf%n_ao, 1)
      wf%orbital_energies = zero
!
      call get_ao_s_xy(wf%ao_overlap)
!
      call mem%alloc(work, 4*wf%n_ao, 1)
      work = zero
!
      call mem%alloc(ao_fock_copy, wf%n_ao, wf%n_ao)
!
      ao_fock_copy = zero
      ao_fock_copy = wf%ao_fock
!
      call dsygv(1, 'V',               &
                  'L',                 &
                  wf%n_ao,             &
                  wf%ao_fock,          & ! contains the orbital coefficients on exit
                  wf%n_ao,             &
                  wf%ao_overlap,       & ! gets scrambled
                  wf%n_ao,             &
                  wf%orbital_energies, &
                  work,                &
                  4*(wf%n_ao),         &
                  info)
!
      write(output%unit,*) 'info:', info
!
      call mem%dealloc(work, 4*wf%n_ao, 1)
!
      wf%mo_coefficients = zero
      wf%mo_coefficients = wf%ao_fock
!
     ! wf%ao_fock = ao_fock_copy
      call wf%construct_ao_density()
      call wf%construct_ao_fock()
!
      call mem%dealloc(ao_fock_copy, wf%n_ao, wf%n_ao)
!
!     Transform AO Fock matrix to MO Fock matrix and set F_ia
!
      call mem%alloc(X_mat, wf%n_ao, wf%n_ao)
      call mem%alloc(Y_mat, wf%n_ao, wf%n_ao)
!
      call dgemm('N', 'N',            &
                  wf%n_ao,            &
                  wf%n_ao,            &
                  wf%n_ao,            &
                  one,                &
                  wf%ao_fock,         &
                  wf%n_ao,            &
                  wf%mo_coefficients, &
                  wf%n_ao,            &
                  zero,               &
                  X_mat,              & ! X = F^ao C
                  wf%n_ao)
!
      call dgemm('T', 'N',            &
                  wf%n_ao,            &
                  wf%n_ao,            &
                  wf%n_ao,            &
                  one,                &
                  wf%mo_coefficients, &
                  wf%n_ao,            &
                  X_mat,              &
                  wf%n_ao,            &
                  zero,               &
                  Y_mat,              & ! Y = C^T F^ao C
                  wf%n_ao)
!
      call mem%dealloc(X_mat, wf%n_ao, wf%n_ao)
!
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            ia = index_two(i,a,wf%n_o)
            F_ia(ia, 1) = Y_mat(i, wf%n_o + a)
!
            write(output%unit,*) 'f_ia', F_ia(ia, 1), Y_mat(i, wf%n_o + a)
!
         enddo
      enddo
!
      norm = ddot((wf%n_o)*(wf%n_v), F_ia, 1, F_ia, 1)
!
      write(output%unit, *) 'Error norm: ', sqrt(norm)
      flush(output%unit)
!
      call mem%dealloc(Y_mat, wf%n_ao, wf%n_ao)
      call mem%dealloc(wf%orbital_energies, wf%n_ao, 1)
!
   end subroutine calculate_hf_equations_hf
!
!
   subroutine get_ao_density_from_mo_coefficients_hf(wf, D_xy)
!!
!!    Get AO density from MO coefficients
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(:,:) :: D_xy
!
      call wf%construct_ao_density()
      call packin(D_xy, wf%ao_density, wf%n_ao)
!
   end subroutine get_ao_density_from_mo_coefficients_hf
!
end module hf_class
