module hf_class
!
!!
!!    Hartree-Fock (HF) class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   use wavefunction_class
!
   use reordering
   use array_utilities
   use interval_class
   use index
!
   implicit none
!
   type, extends(wavefunction):: hf
!
      real(dp) :: hf_energy
!
      real(dp), dimension(:,:), allocatable :: ao_density
      real(dp), dimension(:,:), allocatable :: ao_overlap
      real(dp), dimension(:,:), allocatable :: ao_fock
      real(dp), dimension(:,:), allocatable :: orbital_energies
!
      type(integral_manager) :: integrals
!
	contains
!
!     Initialize and finalize wavefunction
!
      procedure :: initialize => initialize_hf
      procedure :: finalize   => finalize_hf
!
!     Construction routines for various wavefunction variables
!
      procedure :: construct_ao_density => construct_ao_density_hf
      procedure :: construct_ao_fock    => construct_ao_fock_hf
      procedure :: construct_mo_fock    => construct_mo_fock_hf
      procedure :: construct_ao_overlap => construct_ao_overlap_hf
      procedure :: calculate_hf_energy  => calculate_hf_energy_hf
!
!     Rotate and purification routines for the AO density
!
      procedure :: rotate_ao_density    => rotate_ao_density_hf
      procedure :: purify_ao_density    => purify_ao_density_hf
!
!     Cholesky decomposition of AO density and overlap
!
      procedure :: decompose_ao_density => decompose_ao_density_hf
      procedure :: decompose_ao_overlap => decompose_ao_overlap_hf
!
!     Solve the Roothan Hall equation FC = SCe by diagonalization
!
      procedure :: solve_roothan_hall => solve_roothan_hall_hf
!
!     Get and set routines for wavefunction variables
!
      procedure :: get_ao_density => get_ao_density_hf
      procedure :: get_fock_ov    => get_fock_ov_hf
!
      procedure :: set_ao_density   => set_ao_density_hf
!
      procedure :: set_SAD => set_SAD_hf
!
!     Initialize and destruct routines for wavefunction variables
!
      procedure :: initialize_ao_density       => initialize_ao_density_hf
      procedure :: initialize_ao_fock          => initialize_ao_fock_hf
      procedure :: initialize_mo_coefficients  => initialize_mo_coefficients_hf
      procedure :: initialize_ao_overlap       => initialize_ao_overlap_hf
      procedure :: initialize_orbital_energies => initialize_orbital_energies_hf
!
      procedure :: destruct_ao_density         => destruct_ao_density_hf
      procedure :: destruct_ao_fock            => destruct_ao_fock_hf
      procedure :: destruct_mo_coefficients    => destruct_mo_coefficients_hf
      procedure :: destruct_ao_overlap         => destruct_ao_overlap_hf
      procedure :: destruct_orbital_energies   => destruct_orbital_energies_hf
!
!     Routines that construct different components of the Roothan-Hall 1st order Newton equations
!
      procedure :: construct_projection_matrices               => construct_projection_matrices_hf
      procedure :: project_redundant_rotations                 => project_redundant_rotations_hf
      procedure :: construct_roothan_hall_hessian              => construct_roothan_hall_hessian_hf
      procedure :: construct_roothan_hall_gradient             => construct_roothan_hall_gradient_hf
      procedure :: construct_stationary_roothan_hall_condition => construct_stationary_roothan_hall_condition_hf
!
      procedure :: construct_sp_eri_schwarz =>  construct_sp_eri_schwarz_hf
!
   end type hf
!
!
contains
!
!
   subroutine initialize_hf(wf)
!!
!!    Initialize
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Sets low-memory member variables that are needed throughout - in contrast
!!    to larger arrays, which are allocated when needed - and calls initialization
!!    routines for the molecular system object and the Libint integral library.
!!
      implicit none
!
      class(hf) :: wf
!
      wf%name = 'HF'
!
!
      call wf%system%initialize() ! Initialize molecular system -> Should include SOAD
!
      wf%n_ao = 0
      call get_n_aos(wf%n_ao)
!
      wf%n_mo = wf%n_ao
!
      wf%n_o = (wf%system%get_n_electrons())/2
      wf%n_v = wf%n_mo - wf%n_o
!
!     Initialize Libint engines to be used
!
      call initialize_coulomb()
      call initialize_kinetic()
      call initialize_nuclear()
      call initialize_overlap() ! SHOULD THESE BE INITIALIZED IN THE ENGINE ?
      call initialize_overlap()
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
!     Nothing here yet
!
   end subroutine finalize_hf
!
!
   subroutine initialize_ao_density_hf(wf)
!!
!!    Initialize AO density
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      call mem%alloc(wf%ao_density, wf%n_ao, wf%n_ao)
!
   end subroutine initialize_ao_density_hf
!
!
   subroutine initialize_orbital_energies_hf(wf)
!!
!!    Initialize orbital energies
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
       call mem%alloc(wf%orbital_energies, wf%n_mo, 1)
       wf%orbital_energies = zero
!
   end subroutine initialize_orbital_energies_hf
!
!
   subroutine initialize_ao_fock_hf(wf)
!!
!!    Initialize AO Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      call mem%alloc(wf%ao_fock, wf%n_ao, wf%n_ao)
      wf%ao_fock = zero
!
   end subroutine initialize_ao_fock_hf
!
!
   subroutine initialize_mo_coefficients_hf(wf)
!!
!!    Initialize MO coefficients
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      call mem%alloc(wf%orbital_coefficients, wf%n_ao, wf%n_mo)
      wf%orbital_coefficients = zero
!
   end subroutine initialize_mo_coefficients_hf
!
!
   subroutine initialize_ao_overlap_hf(wf)
!!
!!    Initialize AO overlap
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      call mem%alloc(wf%ao_overlap, wf%n_ao, wf%n_ao)
      wf%ao_overlap = zero
!
   end subroutine initialize_ao_overlap_hf
!
!
   subroutine destruct_ao_overlap_hf(wf)
!!
!!    Destruct AO overlap
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      call mem%dealloc(wf%ao_overlap, wf%n_ao, wf%n_ao)
!
   end subroutine destruct_ao_overlap_hf
!
!
   subroutine destruct_orbital_energies_hf(wf)
!!
!!    Destruct AO overlap
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      call mem%dealloc(wf%orbital_energies, wf%n_mo, 1)
!
   end subroutine destruct_orbital_energies_hf
!
!
   subroutine destruct_ao_density_hf(wf)
!!
!!    Destruct AO density
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      call mem%dealloc(wf%ao_density, wf%n_ao, wf%n_ao)
!
   end subroutine destruct_ao_density_hf
!
!
   subroutine destruct_ao_fock_hf(wf)
!!
!!    Destruct AO Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      call mem%dealloc(wf%ao_fock, wf%n_ao, wf%n_ao)
!
   end subroutine destruct_ao_fock_hf
!
!
   subroutine destruct_mo_coefficients_hf(wf)
!!
!!    Destruct MO coefficients
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      call mem%dealloc(wf%orbital_coefficients, wf%n_ao, wf%n_mo)
!
   end subroutine destruct_mo_coefficients_hf
!
!
   subroutine set_SAD_hf(wf)
!!
!!    Set SAD
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      integer(i15) :: ao
!
      real(dp), dimension(:,:), allocatable :: density_diagonal
!
      call wf%initialize_ao_density()
!
!     Set initial density to superposition of atomic densities (SAD) guess
!
      call mem%alloc(density_diagonal, wf%n_ao, 1)
      call wf%system%SAD(wf%n_ao, density_diagonal)
!
      do ao = 1, wf%n_ao
!
         wf%ao_density(ao, ao) = density_diagonal(ao, 1)
!
      enddo
!
      call mem%dealloc(density_diagonal, wf%n_ao, 1)
!
   end subroutine set_SAD_hf
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
!!    Density is packed
!!
      implicit none
!
      class(hf) :: wf
!
      integer(i15) :: i, x, y, xy
!
      wf%ao_density = zero
!
      call dgemm('N', 'T',                   &
                  wf%n_ao,                   &
                  wf%n_ao,                   &
                  wf%n_o,                    &
                  two,                       &
                  wf%orbital_coefficients,   &
                  wf%n_ao,                   &
                  wf%orbital_coefficients,   &
                  wf%n_ao,                   &
                  zero,                      &
                  wf%ao_density,             &
                  wf%n_ao)
!
   end subroutine construct_ao_density_hf
!
!
   subroutine construct_sp_eri_schwarz_hf(wf, sp_eri_schwarz, n_s)
!!
!!
!!
      implicit none
!
      class(hf) :: wf 
!
      integer(i15) :: n_s
!
      real(dp), dimension(n_s, n_s) :: sp_eri_schwarz
!
!     Local variables
!
      integer(i15) :: s1, s2
!
      real(dp) :: max
!
      real(dp), dimension(:,:), allocatable :: g
!
      type(interval) :: A_interval, B_interval
!
!$omp parallel do private(s1, s2, A_interval, B_interval, g, max) schedule(dynamic)
      do s1 = 1, n_s
         do s2 = 1, s1
!
            A_interval = wf%system%get_shell_limits(s1)
            B_interval = wf%system%get_shell_limits(s2)
!
            call mem%alloc(g, (A_interval%size)*(B_interval%size), &
                              (A_interval%size)*(B_interval%size))
!
            call wf%system%ao_integrals%get_ao_g_wxyz(g, s1, s2, s1, s2)
!
            max = get_abs_max(g, (A_interval%size)*(B_interval%size)**2)
!
            call mem%dealloc(g, (A_interval%size)*(B_interval%size), &
                                (A_interval%size)*(B_interval%size))
!
            sp_eri_schwarz(s1, s2) = max
            sp_eri_schwarz(s2, s1) = max
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine construct_sp_eri_schwarz_hf
!
!
   subroutine determine_degeneracy_hf(wf, degeneracy, n_s)
!!
!!
!!
      implicit none
!
      class(hf) :: wf
!
      integer(i15) :: n_s
!
      real(dp), dimension(n_s**2, n_s**2), allocatable :: degeneracy
!
      integer(i15) :: s1, s2, s3, s4, s1s2, s3s4, deg_12, deg_34, deg_12_34
!
!$omp parallel do private(s1, s2, s3, s4, s1s2, s3s4, deg_12, deg_34, deg_12_34)
      do s1 = 1, n_s
         do s2 = 1, s1
!
            s1s2 = n_s*(s2 - 1) + s1
!
            if (s1 .eq. s2) then
!
               deg_12 = one
!
            else
!
               deg_12 = two
!
            endif
!
            do s3 = 1, s1
!
               if (s3 .eq. s1) then
!
                  s4_max = s2
!
               else
!
                  s4_max = s3
!
               endif
!
               do s4 = 1, s4_max
!
                  if (s3 .eq. s4) then
!
                     deg_34 = one
!
                  else
!
                     deg_34 = two
!
                  endif
!
                  if (s3 .eq. s1) then
!
                     if (s2 .eq. s4) then
!
                        deg_12_34 = one
!
                     else
!
                        deg_12_34 = two
!
                     endif
!
                  else
!
                     deg_12_34 = two
!
                  endif
!
                  s3s4 = n_s*(s4 - 1) + s3
!
                  degeneracy(s1s2, s3s4) = deg_12*deg_34*deg_12_34
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine determine_degeneracy_hf
!
!
   subroutine construct_ao_fock_hf(wf, sp_eri_schwarz, n_s)
!!
!!    Construct AO Fock matrix
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Calculates
!!
!!       F_αβ = h_αβ + sum_γδ g_αβγδ D_γδ - 1/2 * sum_γδ g_αδγβ D_γδ,
!!
!!    where D is the AO density. This routine is integral direct, and
!!    it calculates the Hartree-Fock energy by default.
!!
      implicit none
!
      class(hf) :: wf
!
      integer(i15) :: n_s
!
      real(dp), dimension(n_s, n_s) :: sp_eri_schwarz
!
      real(dp), dimension(:,:), allocatable :: ao_fock_packed
      real(dp), dimension(:,:), allocatable :: X_wz, h_wx, h_wx_square
      integer(i15) :: w, x, y, z, wx, yz, w_red, x_red, y_red, z_red
!
      real(dp), dimension(:,:), allocatable :: F1
      real(dp), dimension(:,:), allocatable :: F2
      real(dp), dimension(:,:), allocatable :: F3
      real(dp), dimension(:,:), allocatable :: F4
      real(dp), dimension(:,:), allocatable :: F5
      real(dp), dimension(:,:), allocatable :: F6
!
      integer(i15) :: s1, s2, s3, s4, s4_max, s1s2, s3s4
!
      type(interval) :: A_interval
      type(interval) :: B_interval
      type(interval) :: C_interval
      type(interval) :: D_interval
!
      real(dp) :: deg_12, deg_34, deg_12_34, deg, ddot, norm
      real(dp) :: temp, temp1, temp2, temp3, temp4, temp5, temp6
      real(dp) :: max, max_D_schwarz, max_eri_schwarz
!
      real(dp), dimension(:,:), allocatable :: degeneracy
      real(dp), dimension(:,:), allocatable :: g, D, sp_density_schwarz
!
      real(dp) :: start_timer, end_timer, omp_get_wtime
!
      start_timer = omp_get_wtime()
!
      call mem%alloc(sp_density_schwarz, n_s, n_s)
!
!$omp parallel do private(s1, s2, A_interval, B_interval, D, max)
      do s1 = 1, n_s
         do s2 = 1, s1
!
            A_interval = wf%system%get_shell_limits(s1)
            B_interval = wf%system%get_shell_limits(s2)
!
            call mem%alloc(D, (A_interval%size), (B_interval%size))
!
            D = wf%ao_density(A_interval%first : A_interval%last, B_interval%first : B_interval%last)
!
            max = get_abs_max(D, (A_interval%size)*(B_interval%size))
!
            call mem%dealloc(D, (A_interval%size), (B_interval%size))
!
            sp_density_schwarz(s1, s2) = max
            sp_density_schwarz(s2, s1) = max
!
         enddo
      enddo
!$omp end parallel do
!
!     Calculate max' prescreening
!
      max_D_schwarz     = get_abs_max(sp_density_schwarz, n_s**2)
      max_eri_schwarz   = get_abs_max(sp_eri_schwarz, n_s**2)
!
      end_timer = omp_get_wtime()
   !   write(output%unit, '(t3,a32,f9.1)') 'Schwarz screening array (sec.): ', end_timer - start_timer
!
      start_timer = omp_get_wtime()
!
      wf%ao_fock = zero
!
!$omp parallel do &
!$omp private(s1, s2, s3, s4, deg, s4_max, temp, s1s2, s3s4, &
!$omp A_interval, B_interval, C_interval, D_interval, w, x, y, z, wx, yz, temp1, temp2, temp3, &
!$omp temp4, temp5, temp6, F1, F2, F3, F4, F5, F6, w_red, x_red, y_red, z_red, g) schedule(dynamic)
      do s1 = 1, n_s
         do s2 = 1, s1
!
            if (sp_eri_schwarz(s1, s2)*(max_D_schwarz)*(max_eri_schwarz) .lt. 1.0d-12) continue
!
            s1s2 = n_s*(s2 - 1) + s1
!
            do s3 = 1, s1
!
               if (s3 .eq. s1) then
!
                  s4_max = s2
!
               else
!
                  s4_max = s3
!
               endif
!
               do s4 = 1, s4_max
!
                  if (sp_eri_schwarz(s1, s2)*(max_D_schwarz)*sp_eri_schwarz(s3, s4) .lt. 1.0d-12) continue
!
                  s3s4 = n_s*(s4 - 1) + s3
                  deg = degeneracy(s1s2, s3s4) ! Shell degeneracy
!
                  A_interval = wf%system%get_shell_limits(s1)
                  B_interval = wf%system%get_shell_limits(s2)
                  C_interval = wf%system%get_shell_limits(s3)
                  D_interval = wf%system%get_shell_limits(s4)
!
                  call mem%alloc(g, (A_interval%size)*(B_interval%size), &
                                    (C_interval%size)*(D_interval%size))
                  call wf%system%ao_integrals%get_ao_g_wxyz(g, s1, s2, s3, s4)
!
                  call mem%alloc(F1, A_interval%size, B_interval%size) ! F_wx
                  call mem%alloc(F6, C_interval%size, B_interval%size) ! F_yx
                  call mem%alloc(F2, C_interval%size, D_interval%size) ! F_yz
                  call mem%alloc(F3, A_interval%size, D_interval%size) ! F_wz
                  call mem%alloc(F4, B_interval%size, D_interval%size) ! F_xz
                  call mem%alloc(F5, A_interval%size, C_interval%size) ! F_wy
!
                  F1 = zero
                  F2 = zero
                  F3 = zero
                  F4 = zero
                  F5 = zero
                  F6 = zero
!
!                 Add Fock matrix contributions
!
                  do y = C_interval%first, C_interval%last
                     do z = D_interval%first, D_interval%last
!
                        y_red = y - C_interval%first + 1
                        z_red = z - D_interval%first + 1
!
                        yz = (C_interval%size)*(z_red - 1) + y_red
!
                        do w = A_interval%first, A_interval%last
                           do x = B_interval%first, B_interval%last
!
                              w_red = w - A_interval%first + 1
                              x_red = x - B_interval%first + 1
!
                              wx = (A_interval%size)*(x_red - 1) + w_red
!
                              temp = deg*g(wx, yz)
!
                              temp1 = (one/two)*temp*wf%ao_density(y, z)
                              temp2 = (one/two)*temp*wf%ao_density(w, x)
                              temp3 = (one/eight)*temp*wf%ao_density(y, x)
                              temp4 = (one/eight)*temp*wf%ao_density(y, w)
                              temp5 = (one/eight)*temp*wf%ao_density(z, x)
                              temp6 = (one/eight)*temp*wf%ao_density(w, z)
!
                              F1(w_red, x_red) = F1(w_red, x_red) + temp1
                              F2(y_red, z_red) = F2(y_red, z_red) + temp2
                              F3(w_red, z_red) = F3(w_red, z_red) - temp3
                              F4(x_red, z_red) = F4(x_red, z_red) - temp4
                              F5(w_red, y_red) = F5(w_red, y_red) - temp5
                              F6(y_red, x_red) = F6(y_red, x_red) - temp6
!
                           enddo
                        enddo
                     enddo
                  enddo
!
                  call mem%dealloc(g, (A_interval%size)*(B_interval%size), &
                                    (C_interval%size)*(D_interval%size))
!
                  do x = B_interval%first, B_interval%last
                     do w = A_interval%first, A_interval%last
!
                        w_red = w - A_interval%first + 1
                        x_red = x - B_interval%first + 1
!$omp atomic update
                        wf%ao_fock(w, x) = wf%ao_fock(w, x) + F1(w_red, x_red)
!
                     enddo
                  enddo
!
                  do x = B_interval%first, B_interval%last
                     do y = C_interval%first, C_interval%last
!
                        x_red = x - B_interval%first + 1
                        y_red = y - C_interval%first + 1
!$omp atomic update
                        wf%ao_fock(y, x) = wf%ao_fock(y, x) + F6(y_red, x_red)
!
                     enddo
                  enddo
!
                  do z = D_interval%first, D_interval%last
                     do y = C_interval%first, C_interval%last
!
                        y_red = y - C_interval%first + 1
                        z_red = z - D_interval%first + 1
!$omp atomic update
                        wf%ao_fock(y, z) = wf%ao_fock(y, z) + F2(y_red, z_red)
!
                     enddo
                  enddo
!
                  do z = D_interval%first, D_interval%last
                     do w = A_interval%first, A_interval%last
!
                        w_red = w - A_interval%first + 1
                        z_red = z - D_interval%first + 1
!$omp atomic update
                        wf%ao_fock(w, z) = wf%ao_fock(w, z) + F3(w_red, z_red)
!
                     enddo
                  enddo
!
                  do z = D_interval%first, D_interval%last
                     do x = B_interval%first, B_interval%last
!
                        x_red = x - B_interval%first + 1
                        z_red = z - D_interval%first + 1
!$omp atomic update
                        wf%ao_fock(x, z) = wf%ao_fock(x, z) + F4(x_red, z_red)
!
                     enddo
                  enddo
!
                  do y = C_interval%first, C_interval%last
                     do w = A_interval%first, A_interval%last
!
                        w_red = w - A_interval%first + 1
                        y_red = y - C_interval%first + 1
!$omp atomic update
                        wf%ao_fock(w, y) = wf%ao_fock(w, y) + F5(w_red, y_red)
!
                     enddo
                  enddo
!
                  call mem%dealloc(F1, A_interval%size, B_interval%size) ! F_wx
                  call mem%dealloc(F6, C_interval%size, B_interval%size) ! F_yx
                  call mem%dealloc(F2, C_interval%size, D_interval%size) ! F_yz
                  call mem%dealloc(F3, A_interval%size, D_interval%size) ! F_wz
                  call mem%dealloc(F4, B_interval%size, D_interval%size) ! F_xz
                  call mem%dealloc(F5, A_interval%size, C_interval%size) ! F_wy
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(sp_density_schwarz, n_s, n_s)
!
      call mem%dealloc(degeneracy, n_s**2, n_s**2)
!
      call symmetric_sum(wf%ao_fock, wf%n_ao)
!
      wf%ao_fock = wf%ao_fock*half
!
      call mem%alloc(h_wx, wf%n_ao, wf%n_ao)
      call get_ao_h_xy(h_wx)
!
      call wf%calculate_hf_energy(wf%ao_fock, h_wx)
!
      wf%ao_fock = wf%ao_fock + h_wx
!
      call mem%dealloc(h_wx, wf%n_ao*(wf%n_ao+1)/2, 1)
!
      end_timer = omp_get_wtime()
!
   end subroutine construct_ao_fock_hf
!
!
   subroutine calculate_hf_energy_hf(wf, half_GD_wx, h_wx)
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
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: half_GD_wx
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: h_wx
!
      wf%hf_energy = wf%system%get_nuclear_repulsion()
!
      wf%hf_energy = wf%hf_energy + ddot((wf%n_ao)**2, h_wx, 1, wf%ao_density, 1)
      wf%hf_energy = wf%hf_energy + two*(one/four)*ddot((wf%n_ao)**2, wf%ao_density, 1, half_GD_wx, 1)

   end subroutine calculate_hf_energy_hf
!
!
   subroutine set_ao_density_hf(wf, D)
!!
!!    Set AO density
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Sets the AO density from input
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(:,:) :: D ! Packed
!
      call squareup(D, wf%ao_density, wf%n_ao)
!
   end subroutine set_ao_density_hf
!
!
   subroutine solve_roothan_hall_hf(wf)
!!
!!    Solve Roothan Hall
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(:,:), allocatable :: work
      real(dp), dimension(:,:), allocatable :: ao_overlap_copy
!
      real(dp) :: ddot, norm
!
      integer(i15) :: info = 0
!
      wf%orbital_energies = zero
!
      call mem%alloc(work, 4*wf%n_ao, 1)
      work = zero
!
      call mem%alloc(ao_overlap_copy, wf%n_ao, wf%n_ao)
!
      wf%orbital_coefficients = wf%ao_fock
!
      ao_overlap_copy = wf%ao_overlap
!
      call dsygv(1, 'V', 'L',                &
                  wf%n_ao,                   &
                  wf%orbital_coefficients,   & ! ao_fock on entry orbital coefficients on exit
                  wf%n_ao,                   &
                  ao_overlap_copy,           &
                  wf%n_ao,                   &
                  wf%orbital_energies,       &
                  work,                      &
                  4*(wf%n_ao),               &
                  info)
!
      call mem%dealloc(work, 4*wf%n_ao, 1)
      call mem%dealloc(ao_overlap_copy, wf%n_ao, wf%n_ao)
!
   end subroutine solve_roothan_hall_hf
!
!
   subroutine construct_ao_overlap_hf(wf)
!!
!!    Construct AO overlap
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      call get_ao_s_xy(wf%ao_overlap)
!
   end subroutine construct_ao_overlap_hf
!
!
   subroutine construct_mo_fock_hf(wf, F_pq)
!!
!!    Construct MO Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(:,:) :: F_pq
!
      real(dp), dimension(:,:), allocatable :: X
!
      call mem%alloc(X, wf%n_ao, wf%n_ao)
!
      call dgemm('N', 'N',                   &
                  wf%n_ao,                   &
                  wf%n_ao,                   &
                  wf%n_ao,                   &
                  one,                       &
                  wf%ao_fock,                &
                  wf%n_ao,                   &
                  wf%orbital_coefficients,   &
                  wf%n_ao,                   &
                  zero,                      &
                  X,                         & ! X = F^ao C
                  wf%n_ao)
!
      call dgemm('T', 'N',                   &
                  wf%n_ao,                   &
                  wf%n_ao,                   &
                  wf%n_ao,                   &
                  one,                       &
                  wf%orbital_coefficients,   &
                  wf%n_ao,                   &
                  X,                         &
                  wf%n_ao,                   &
                  zero,                      &
                  F_pq,                      & ! F = C^T F^ao C
                  wf%n_ao)
!
      call mem%dealloc(X, wf%n_ao, wf%n_ao)
!
   end subroutine construct_mo_fock_hf
!
!
   subroutine get_fock_ov_hf(wf, F)
!!
!!    Get HF equations
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Constructs the occupied-virtual block of the Fock MO matrix,
!!    and returns the result in the array F.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_v)   :: F ! F_ia
      real(dp), dimension(:,:), allocatable :: X
!
      call mem%alloc(X, wf%n_ao, wf%n_v)
!
      call dgemm('N', 'N',                                  &
                  wf%n_ao,                                  &
                  wf%n_v,                                   &
                  wf%n_ao,                                  &
                  one,                                      &
                  wf%ao_fock,                               &
                  wf%n_ao,                                  &
                  wf%orbital_coefficients(1, wf%n_o + 1),   &
                  wf%n_ao,                                  &
                  zero,                                     &
                  X,                                        &
                  wf%n_ao)
!
      call dgemm('T', 'N',                   &
                  wf%n_o,                    &
                  wf%n_v,                    &
                  wf%n_ao,                   &
                  one,                       &
                  wf%orbital_coefficients,   &
                  wf%n_ao,                   &
                  X,                         &
                  wf%n_ao,                   &
                  zero,                      &
                  F,                         &
                  wf%n_o)
!
      call mem%dealloc(X, wf%n_ao, wf%n_v)
!
   end subroutine get_fock_ov_hf
!
!
   subroutine get_ao_density_hf(wf, D)
!!
!!    Get AO density
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Packs the wavefunction's AO density into D.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(:,:) :: D
!
      call packin(D, wf%ao_density, wf%n_ao)
!
   end subroutine get_ao_density_hf
!
!
   subroutine decompose_ao_density_hf(wf)
!!
!!    Decompose AO density
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Does a Cholesky decomposition of the AO density matrix,
!!
!!       D^AO_wx = sum_J L_w,J L_J,x^T,
!!
!!    and sets the MO coefficients accordingly.
!!
      implicit none
!
      class(hf) :: wf
!
      integer(i15), dimension(:, :), allocatable :: used_diag
!
      real(dp), dimension(:,:), allocatable :: perm_matrix
!
      real(dp), dimension(:,:), allocatable :: csc
      real(dp), dimension(:,:), allocatable :: tmp
!
      integer(i15) :: rank, i, j
!
      allocate(used_diag(wf%n_ao, 1))
!
      wf%ao_density = half*wf%ao_density
      call full_cholesky_decomposition_system(wf%ao_density, wf%orbital_coefficients, wf%n_ao, rank,&
                                                      1.0D-12, used_diag)
      wf%ao_density = two*wf%ao_density
!
!     Make permutation matrix P
!
      call mem%alloc(perm_matrix, wf%n_ao, wf%n_ao)
!
      perm_matrix = zero
!
      do j = 1, wf%n_ao
!
         perm_matrix(used_diag(j,1), j) = one
!
      enddo
!
      deallocate(used_diag)
!
!     Sanity check
!
      call mem%alloc(tmp, wf%n_ao, wf%n_ao)
!
      call dgemm('N','N', &
                  wf%n_ao, &
                  wf%n_ao, &
                  wf%n_ao, &
                  one, &
                  perm_matrix, &
                  wf%n_ao, &
                  wf%orbital_coefficients, &
                  wf%n_ao, &
                  zero, &
                  tmp, &
                  wf%n_ao)
!
      wf%orbital_coefficients = tmp
!
      call mem%dealloc(tmp, wf%n_ao, wf%n_ao)
!
   end subroutine decompose_ao_density_hf
!
!
   subroutine decompose_ao_overlap_hf(wf, L)
!!
!!    Decompose AO overlap
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Does a Cholesky decomposition of the AO density matrix,
!!
!!       S^AO_wx = sum_J L_w,J L_J,x^T,
!!
!!    and sets the MO coefficients accordingly.
!!
      implicit none
!
      class(hf) :: wf
!
      integer(i15), dimension(:), allocatable :: used_diag
!
      real(dp), dimension(wf%n_ao, wf%n_ao) :: L
!
      real(dp), dimension(:,:), allocatable :: perm_matrix
!
      real(dp), dimension(:,:), allocatable :: csc
      real(dp), dimension(:,:), allocatable :: tmp
!
      integer(i15) :: rank, i, j
!
      allocate(used_diag(wf%n_ao))
!
      L = zero
!
      call full_cholesky_decomposition_system(wf%ao_overlap, L, wf%n_ao, rank, &
                                                      1.0D-32, used_diag)
!
      if (rank .lt. wf%n_ao) write(output%unit, *) 'Warning: rank lower than full dim for S = L L^T'
!
!     Make permutation matrix P
!
      call mem%alloc(perm_matrix, wf%n_ao, wf%n_ao)
!
      perm_matrix = zero
!
      do j = 1, wf%n_ao
!
         perm_matrix(used_diag(j), j) = one
!
      enddo
!
      deallocate(used_diag)
!
!     Make matrix L
!
      call mem%alloc(tmp, wf%n_ao, wf%n_ao)
!
      call dgemm('N','N', &
                  wf%n_ao, &
                  wf%n_ao, &
                  wf%n_ao, &
                  one, &
                  perm_matrix, &
                  wf%n_ao, &
                  L, &
                  wf%n_ao, &
                  zero, &
                  tmp, &
                  wf%n_ao)
!
      L = tmp
!
      call mem%dealloc(tmp, wf%n_ao, wf%n_ao)
      call mem%dealloc(perm_matrix, wf%n_ao, wf%n_ao)
!
   end subroutine decompose_ao_overlap_hf
!
!
   subroutine rotate_ao_density_hf(wf, X)
!!
!!    Rotate AO density
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Performs an update of the AO density according to a first-order
!!    truncation of the BCH expansion:
!!
!!       D^AO <- exp(-X S) D^AO exp(S X) ~ D^AO + [D^AO, X]_S + 1/2 [[D^AO, X]_S, X]_S
!!                                       ~ D_AO + C + D.
!!
!!    This routine does not - currently - use the second order term, D, although the
!!    code is present (commented out).
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: X
!
      real(dp), dimension(:,:), allocatable :: M
      real(dp), dimension(:,:), allocatable :: C
      real(dp), dimension(:,:), allocatable :: prev_D
!
   !   real(dp), dimension(:,:), allocatable :: D
!
!     Construct C = [D^AO, X]_S = D^AO S X - X S D^AO
!
      call mem%alloc(M, wf%n_ao, wf%n_ao)
!
      call dgemm('N','N',        &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  one,           &
                  wf%ao_overlap, &
                  wf%n_ao,       &
                  X,             &
                  wf%n_ao,       &
                  zero,          &
                  M,             & ! M = S X
                  wf%n_ao)
!
   !   call mem%alloc(D, wf%n_ao, wf%n_ao)
      call mem%alloc(C, wf%n_ao, wf%n_ao)
!
      call dgemm('N','N',        &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  one,           &
                  wf%ao_density, &
                  wf%n_ao,       &
                  M,             &
                  wf%n_ao,       &
                  zero,          &
                  C,             & ! C = D^AO M = D^AO S X
                  wf%n_ao)
!
      call dgemm('N','N',        &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  one,           &
                  wf%ao_overlap, &
                  wf%n_ao,       &
                  wf%ao_density, &
                  wf%n_ao,       &
                  zero,          &
                  M,             & ! M = S D^AO
                  wf%n_ao)
!
      call dgemm('N','N',        &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  -one,          &
                  X,             &
                  wf%n_ao,       &
                  M,             &
                  wf%n_ao,       &
                  one,           &
                  C,             & ! C = C - X M = C - X S D^AO = D^AO S X - X S D^AO
                  wf%n_ao)
!
!     Calculate D = 1/2 [[D^AO, X]_S, X]_S = 1/2 [C, X]_S = 1/2 (C S X - X S C)
!
!       call dgemm('N','N',        &
!                   wf%n_ao,       &
!                   wf%n_ao,       &
!                   wf%n_ao,       &
!                   one,           &
!                   wf%ao_overlap, &
!                   wf%n_ao,       &
!                   X,             &
!                   wf%n_ao,       &
!                   zero,          &
!                   M,             & ! M = S X
!                   wf%n_ao)
! !
!       call dgemm('N','N',        &
!                   wf%n_ao,       &
!                   wf%n_ao,       &
!                   wf%n_ao,       &
!                   half,          &
!                   C,             &
!                   wf%n_ao,       &
!                   M,             &
!                   wf%n_ao,       &
!                   zero,          &
!                   D,             & ! D = 1/2 * C S X
!                   wf%n_ao)
! !
!       call dgemm('N','N',        &
!                   wf%n_ao,       &
!                   wf%n_ao,       &
!                   wf%n_ao,       &
!                   one,           &
!                   X,             &
!                   wf%n_ao,       &
!                   wf%ao_overlap, &
!                   wf%n_ao,       &
!                   zero,          &
!                   M,             & ! M = X S
!                   wf%n_ao)
! !
!       call dgemm('N','N',        &
!                   wf%n_ao,       &
!                   wf%n_ao,       &
!                   wf%n_ao,       &
!                   -half,         &
!                   M,             &
!                   wf%n_ao,       &
!                   C,             &
!                   wf%n_ao,       &
!                   one,           &
!                   D,             & ! D = D - 1/2 * X S C
!                   wf%n_ao)
!
      wf%ao_density = wf%ao_density + C
   !   wf%ao_density = wf%ao_density + C + D
!
      call mem%dealloc(M, wf%n_ao, wf%n_ao)
    !  call mem%dealloc(D, wf%n_ao, wf%n_ao)
      call mem%dealloc(C, wf%n_ao, wf%n_ao)
!
   end subroutine rotate_ao_density_hf
!
!
   subroutine purify_ao_density_hf(wf, threshold)
!!
!!    Purify AO density
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Purifies a non-idempotent AO density matrix - typically arising from
!!    the non-exact rotation of the density - by the following fixed point algorithm:
!!
!!       D^AO <- 3/4 D^AO S D^AO - 1/2 D^AO S D^AO S D^AO
!!
!!    To check whether the final result is consistent, it is possible to verify that
!!    1/2 D^AO S is idempotent (which was done during debug of routine).
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), intent(in) :: threshold
!
      real(dp), dimension(:,:), allocatable :: M ! Arrays to temporarily hold matrix products
      real(dp), dimension(:,:), allocatable :: N ! Arrays to temporarily hold matrix products
!
      real(dp), dimension(:,:), allocatable :: prev_ao_density ! Holds previous density matrix
!
      real(dp) :: ddot, error
!
      logical :: pure = .false.
!
      integer(i15) :: iteration, p, q
      integer(i15), parameter :: max_iterations = 50
!
      iteration = 1
!
      call mem%alloc(M, wf%n_ao, wf%n_ao)
      call mem%alloc(N, wf%n_ao, wf%n_ao)
      call mem%alloc(prev_ao_density, wf%n_ao, wf%n_ao)
!
      pure = .false.
      error = zero
!
      do while (.not. pure .and. iteration .le. max_iterations)
!
         prev_ao_density = wf%ao_density
!
         call dgemm('N','N',        &
                     wf%n_ao,       &
                     wf%n_ao,       &
                     wf%n_ao,       &
                     one,           &
                     wf%ao_overlap, &
                     wf%n_ao,       &
                     wf%ao_density, &
                     wf%n_ao,       &
                     zero,          &
                     M,             & ! M = S D^AO
                     wf%n_ao)
!
         call dgemm('N','N',        &
                     wf%n_ao,       &
                     wf%n_ao,       &
                     wf%n_ao,       &
                     one,           &
                     wf%ao_density, &
                     wf%n_ao,       &
                     M,             &
                     wf%n_ao,       &
                     zero,          &
                     N,             & ! N = D^AO M = D^AO S D^AO
                     wf%n_ao)
!
         wf%ao_density = (three/two)*N ! D^AO = 3 D^AO S D^AO
!
         call dgemm('N','N',        &
                     wf%n_ao,       &
                     wf%n_ao,       &
                     wf%n_ao,       &
                     one,           &
                     wf%ao_overlap, &
                     wf%n_ao,       &
                     N,             &
                     wf%n_ao,       &
                     zero,          &
                     M,             & ! M = S N = S D^AO S D^AO
                     wf%n_ao)
!
         call dgemm('N','N',          &
                     wf%n_ao,         &
                     wf%n_ao,         &
                     wf%n_ao,         &
                     one,             &
                     prev_ao_density, &
                     wf%n_ao,         &
                     M,               &
                     wf%n_ao,         &
                     zero,            &
                     N,               & ! N = D^AO M = D^AO S D^AO S D^AO
                     wf%n_ao)
!
         wf%ao_density = wf%ao_density - (one/two)*N
!
         M = wf%ao_density - prev_ao_density
!
         error = sqrt(ddot((wf%n_ao)**2, M, 1, M, 1))
!
         if (error .lt. threshold) then
!
            pure = .true.
!
         endif
!
         iteration = iteration + 1
!
      enddo
!
      if (.not. pure) then
!
         write(output%unit, '(t3,a49,f16.12)') 'Error: could not purify AO density. Final error: ', error
         stop
!
      endif
!
      call mem%dealloc(M, wf%n_ao, wf%n_ao)
      call mem%dealloc(N, wf%n_ao, wf%n_ao)
      call mem%dealloc(prev_ao_density, wf%n_ao, wf%n_ao)
!
   end subroutine purify_ao_density_hf
!
!
   subroutine construct_projection_matrices_hf(wf, Po, Pv)
!!
!!    Construct projection matrices
!!    Written by Eirik F. Kjønstad, 2018
!!
!!    Constructs Po = D S
!!               Pv = 1 - Po
!!
!!    D: AO density matrix, S: AO overlap matrix. Both are assumed
!!    to be allocated and properly set.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao) :: Po
      real(dp), dimension(wf%n_ao, wf%n_ao) :: Pv
!
      real(dp), dimension(:,:), allocatable :: tmp
!
      integer(i15) :: w, x
!
!     Po = D S
!
      call dgemm('N','N',        &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  half,          &
                  wf%ao_density, &
                  wf%n_ao,       &
                  wf%ao_overlap, &
                  wf%n_ao,       &
                  zero,          &
                  Po,            &
                  wf%n_ao)
!
!     Pv = I
!
      Pv = zero
!
      do x = 1, wf%n_ao
!
         Pv(x, x) = one
!
      enddo
!
!     Pv = I - Po
!
      do x = 1, wf%n_ao
         do w = 1, wf%n_ao
!
            Pv(w, x) = Pv(w, x) - Po(w, x)
!
         enddo
      enddo
!
   end subroutine construct_projection_matrices_hf
!
!
   subroutine project_redundant_rotations_hf(wf, X, Po, Pv)
!!
!!    Project redundant rotations
!!    Written by Eirik F. Kjønstad, 2018
!!
!!    Here, X is an antisymmetric rotations matrix on entering the routine,
!!    where some parameters are redundant for rotating the AO density. On exit,
!!    the redundant parameters have been projected out of X.
!!
!!    To achieve this, we set
!!
!!       X <- Po X Pv^T + Pv X Po^T,
!!
!!    where Po = D S and Pv = 1 - Po. In Po, D is the AO density and S
!!    is the AO overlap matrix.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao) :: X
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: Po
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: Pv
!
      real(dp), dimension(:, :), allocatable :: tmp
!
      integer(i15) :: p, q
!
!     Construct
!
!        tmp = X Pv^T   =>   tmp^T = Pv X^T = - Pv X
!
      call mem%alloc(tmp, wf%n_ao, wf%n_ao)
!
      call dgemm('N', 'T', &
                  wf%n_ao, &
                  wf%n_ao, &
                  wf%n_ao, &
                  one,     &
                  X,       &
                  wf%n_ao, &
                  Pv,      &
                  wf%n_ao, &
                  zero,    &
                  tmp,     &
                  wf%n_ao)
!
!     X = Po X Pv^T = Po tmp
!
      call dgemm('N', 'N', &
                  wf%n_ao, &
                  wf%n_ao, &
                  wf%n_ao, &
                  one,     &
                  Po,      &
                  wf%n_ao, &
                  tmp,     &
                  wf%n_ao, &
                  zero,    &
                  X,       &
                  wf%n_ao)
!
!     X = X + Pv X Po^T = X - (-Pv X) Po^T = X - tmp^T Po^T
!
      call dgemm('T', 'T', &
                  wf%n_ao, &
                  wf%n_ao, &
                  wf%n_ao, &
                  -one,    &
                  tmp,     &
                  wf%n_ao, &
                  Po,      &
                  wf%n_ao, &
                  one,     &
                  X,       &
                  wf%n_ao)
!
      call mem%dealloc(tmp, wf%n_ao, wf%n_ao)
!
   end subroutine project_redundant_rotations_hf
!
!
   subroutine construct_roothan_hall_hessian_hf(wf, H, Po, Pv)
!!
!!    Construct Roothan-Hall Hessian
!!    Written by Eirik F. Kjønstad, 2018
!!
!!    Constructs the Roothan-Hall Hessian,
!!
!!       H = Fvv - Foo = Pv^T F Pv - Po^T F Po,
!!
!!    where Po = D S and Pv = 1 - Po. In Po, D is the AO density and S
!!    is the AO overlap matrix.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: Po
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: Pv
!
      real(dp), dimension(wf%n_ao, wf%n_ao) :: H
!
      real(dp), dimension(:, :), allocatable :: tmp
!
!     Construct tmp = Fvv = Pv^T F Pv and set H = tmp = Fvv
!
      call mem%alloc(tmp, wf%n_ao, wf%n_ao)
!
      tmp = wf%ao_fock
      call sandwich(tmp, Pv, Pv, wf%n_ao)
!
      H = tmp
!
!     Construct tmp = Foo = Po^T F Po and set H = H - tmp = Fvv - Foo
!
      tmp = wf%ao_fock
      call sandwich(tmp, Po, Po, wf%n_ao)
!
      H = H - tmp
!
      call mem%dealloc(tmp, wf%n_ao, wf%n_ao)
!
   end subroutine construct_roothan_hall_hessian_hf
!
!
   subroutine construct_roothan_hall_gradient_hf(wf, G, Po, Pv)
!!
!!    Construct Roothan-Hall gradient
!!    Written by Eirik F. Kjønstad, 2018
!!
!!    Constructs the Roothan-Hall gradient,
!!
!!       G = Fov - Fvo = Po^T F Pv - Pv^T F Po,
!!
!!    where Po = D S and Pv = 1 - Po. In Po, D is the AO density and S
!!    is the AO overlap matrix.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: Po
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: Pv
!
      real(dp), dimension(wf%n_ao, wf%n_ao) :: G
!
      real(dp), dimension(:, :), allocatable :: tmp
!
!     Construct tmp = Fov = Po^T F Pv and set G = tmp = Fov
!
      call mem%alloc(tmp, wf%n_ao, wf%n_ao)
!
      tmp = wf%ao_fock
      call sandwich(tmp, Po, Pv, wf%n_ao)
!
      G = tmp
!
!     Construct tmp = Fvo = Pv^T F Po and set H = H - tmp = Fov - Fvo
!
      tmp = wf%ao_fock
      call sandwich(tmp, Pv, Po, wf%n_ao)
!
      G = G - tmp
!
      call mem%dealloc(tmp, wf%n_ao, wf%n_ao)
!
   end subroutine construct_roothan_hall_gradient_hf
!
!
   subroutine construct_stationary_roothan_hall_condition_hf(wf, RHC, H, X, G, S, level_shift)
!!
!!    Construct stationary Roothan-Hall condition
!!    Written by Eirik F. Kjønstad, 2018
!!
!!    Sets
!!
!!       RHC = H X S + S X H + G - level_shift S,
!!
!!    which equals zero on convergence of the Roothan-Hall Newton equations.
!!    Note that if similarity transformed H, S, and G are used (Y <- V-1 Y V-T),
!!    then the iterated solution X' = V^T X V, from which the actual X is easily
!!    extractable.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao) :: RHC
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: H
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: X
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: G
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: S
!
      real(dp), intent(in), optional :: level_shift
!
      real(dp), dimension(:,:), allocatable :: tmp
!
!     Construct tmp = X S => tmp^T = S^T X^T = S X^T = - S X
!
      call mem%alloc(tmp, wf%n_ao, wf%n_ao)
!
      call dgemm('N', 'N',       &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  one,           &
                  X,             &
                  wf%n_ao,       &
                  S,             &
                  wf%n_ao,       &
                  zero,          &
                  tmp,           &
                  wf%n_ao)
!
!     RHC = H X S = H tmp
!
      call dgemm('N', 'N',       &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  one,           &
                  H,             &
                  wf%n_ao,       &
                  tmp,           &
                  wf%n_ao,       &
                  zero,          &
                  RHC,           &
                  wf%n_ao)
!
!     RHC = RHC - (- S X) H = RHC - tmp^T H = H X S + S X H
!
      call dgemm('T', 'N',       &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  -one,          &
                  tmp,           &
                  wf%n_ao,       &
                  H,             &
                  wf%n_ao,       &
                  one,           &
                  RHC,           &
                  wf%n_ao)
!
      call mem%dealloc(tmp, wf%n_ao, wf%n_ao)
!
!     RHC = RHC + G = H X S + S X H + G
!
      RHC = RHC + G
!
      if (present(level_shift)) then 
!
!        RHC = RHC - level_shift S = H X S + S X H + G - level_shift S
!
         RHC = RHC - level_shift*S
!
      endif
!
   end subroutine construct_stationary_roothan_hall_condition_hf
!
!
end module hf_class
