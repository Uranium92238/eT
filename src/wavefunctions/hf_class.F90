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
   use interval_class
   use index
   use integral_manager_class
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
!
      real(dp), dimension(:,:), allocatable :: orbital_energies
!
      real(dp), dimension(:,:), allocatable :: g_wxyz
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
      procedure :: construct_gradient   => construct_gradient_hf
!
      procedure :: calculate_hf_energy  => calculate_hf_energy_hf
!
!     Solve the Roothan Hall equation FC = SCe by diagonalization
!
      procedure :: solve_roothan_hall => solve_roothan_hall_hf
!
!     Get and set routines for wavefunction variables
!
      procedure :: get_ao_density      => get_ao_density_hf
      procedure :: get_n_hf_parameters => get_n_hf_parameters_hf
      procedure :: get_n_hf_equations  => get_n_hf_equations_hf
      procedure :: get_hf_equations    => get_hf_equations_hf
!
      procedure :: set_ao_density      => set_ao_density_hf
!
!     Initialize and destruct routines for wavefunction variables
!
      procedure :: initialize_ao_density       => initialize_ao_density_hf
      procedure :: initialize_ao_fock          => initialize_ao_fock_hf
      procedure :: initialize_mo_coefficients  => initialize_mo_coefficients_hf
      procedure :: initialize_ao_overlap       => initialize_ao_overlap_hf
      procedure :: initialize_orbital_energies => initialize_orbital_energies_hf
!
      procedure :: destruct_ao_density       => destruct_ao_density_hf
      procedure :: destruct_ao_fock          => destruct_ao_fock_hf
      procedure :: destruct_mo_coefficients  => destruct_mo_coefficients_hf
      procedure :: destruct_ao_overlap       => destruct_ao_overlap_hf
      procedure :: destruct_orbital_energies => destruct_orbital_energies_hf
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
      call wf%system%initialize()
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
      wf%ao_density = zero
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
   subroutine construct_ao_fock_hf(wf)
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
      integer(i15) :: w, x, y, z, wx, yz, w_red, x_red, y_red, z_red
!
      real(dp), dimension(:,:), allocatable :: h_wx
!
      real(dp), dimension(:,:), allocatable :: F1
      real(dp), dimension(:,:), allocatable :: F2
      real(dp), dimension(:,:), allocatable :: F3
      real(dp), dimension(:,:), allocatable :: F4
      real(dp), dimension(:,:), allocatable :: F5
      real(dp), dimension(:,:), allocatable :: F6
!
      integer(i15) :: s1, s2, s3, s4, s4_max, n_s, s1s2, s3s4
!
      type(interval) :: A_interval
      type(interval) :: B_interval
      type(interval) :: C_interval
      type(interval) :: D_interval
!
      real(dp) :: deg_12, deg_34, deg_12_34, deg, ddot, norm
      real(dp) :: temp, temp1, temp2, temp3, temp4, temp5, temp6
!
      real(dp), dimension(:,:), allocatable :: degeneracy
      real(dp), dimension(:,:), allocatable :: g
!
      logical, dimension(:, :), allocatable :: schwarz
!
      real(dp) :: start_timer, end_timer, omp_get_wtime
!
      n_s = wf%system%get_n_shells()
!
      start_timer = omp_get_wtime()
!
      call mem%alloc(degeneracy, n_s**2, n_s**2)
      degeneracy = zero
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
      end_timer = omp_get_wtime()
   !   write(output%unit, '(t3,a44,f9.1)') 'Construct integral degeneracy array (sec.): ', end_timer - start_timer
!
      start_timer = omp_get_wtime()
!
      allocate(schwarz(n_s, n_s))
      schwarz = .false.
!
!$omp parallel do private(s1, s2, A_interval, B_interval, g, norm)
      do s1 = 1, n_s
         do s2 = 1, s1
!
            A_interval = wf%system%get_shell_limits(s1)
            B_interval = wf%system%get_shell_limits(s2)
!
            call mem%alloc(g, (A_interval%size)*(B_interval%size), &
                              (A_interval%size)*(B_interval%size))
!
            call wf%integrals%get_ao_g_wxyz(g, s1, s2, s1, s2)
!
            norm = sqrt(ddot((A_interval%size)**2*(B_interval%size)**2, g, 1, g, 1))
!
            if (norm .lt. 1.0D-12) then
!
               schwarz(s1,s2) = .true.
               schwarz(s2,s1) = .true.
!
            endif
!
            call mem%dealloc(g, (A_interval%size)*(B_interval%size), &
                                (A_interval%size)*(B_interval%size))
!
         enddo
      enddo
!$omp end parallel do
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
            if (schwarz(s1, s2)) continue
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
                  call wf%integrals%get_ao_g_wxyz(g, s1, s2, s3, s4)
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
      deallocate(schwarz)
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
      call mem%dealloc(h_wx, wf%n_ao, wf%n_ao)
!
      end_timer = omp_get_wtime()
   !   write(output%unit, '(t3,a59,f9.1)') 'Actual Fock construction and calculation of energy (sec.): ', end_timer - start_timer
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
      class(hf), intent(in) :: wf
!
      integer(i15) :: get_n_hf_equations_hf
!
      get_n_hf_equations_hf = (wf%n_o)*(wf%n_v)
!
   end function get_n_hf_equations_hf
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
      wf%ao_density = zero
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
   subroutine get_hf_equations_hf(wf, F)
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
   end subroutine get_hf_equations_hf
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
   subroutine construct_gradient_hf(wf, gradient_norm)
!!
!!    Construct energy gradient
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Here, ∇E is understood as the derivative of the HF energy with
!!    with respect to X, where the update
!!
!!       D^AO <- exp(-X S) D^AO exp(S X)
!!
!!    is used in density-based Hartree-Fock theory. The full matrix
!!    is calculated as
!!
!!       ∇E = 8 ( M - M^T ),
!!
!!    where M = S D^AO F^AO.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      integer(i15) :: w, x
!
      real(dp) :: gradient_norm
!
      real(dp), dimension(:,:), allocatable :: M
      real(dp), dimension(:,:), allocatable :: gradient ! ∇E
!
!     Calculate M matrix
!
      call mem%alloc(gradient, wf%n_ao, wf%n_ao)
      call mem%alloc(M, wf%n_ao, wf%n_ao)
!
      call dgemm('N','N',        &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  one,           &
                  wf%ao_density, &
                  wf%n_ao,       &
                  wf%ao_fock,    &
                  wf%n_ao,       &
                  zero,          &
                  gradient,      & ! Used temporarily
                  wf%n_ao)
!
      call dgemm('N','N',        &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  one,           &
                  wf%ao_overlap, &
                  wf%n_ao,       &
                  gradient,      &
                  wf%n_ao,       &
                  zero,          &
                  M,             &
                  wf%n_ao)
!
!     Set gradient = M and subtract the transpose, with prefactor 8
!
      do x = 1, wf%n_ao
         do w = 1, wf%n_ao
!
            gradient(w, x) = eight*(M(w, x) - M(x, w))
!
         enddo
      enddo
!
      call mem%dealloc(M, wf%n_ao, wf%n_ao)
!
      gradient_norm = zero
!
      do x = 1, wf%n_ao
         do w = 1, x
!
            gradient_norm = gradient_norm + gradient(w, x)**2
!
         enddo
      enddo
!
      call mem%dealloc(gradient, wf%n_ao, wf%n_ao)
!
      gradient_norm = sqrt(gradient_norm)
!
   end subroutine construct_gradient_hf
!
!
end module hf_class
