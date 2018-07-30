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
      procedure :: initialize => initialize_hf
      procedure :: finalize   => finalize_hf
!
!     Construct various HF arrays (density, AO Fock, ...)
!
      procedure :: construct_ao_density => construct_ao_density_hf
      procedure :: construct_ao_fock    => construct_ao_fock_hf
      procedure :: construct_mo_fock    => construct_mo_fock_hf
      procedure :: construct_ao_overlap => construct_ao_overlap_hf
!
      procedure :: calculate_hf_energy  => calculate_hf_energy_hf
!
!     Solve Roothan Hall equations
!
      procedure :: solve_roothan_hall => solve_roothan_hall_hf
!
!     Get and set routines
!
      procedure :: get_hf_equations    => get_hf_equations_hf
      procedure :: get_ao_density      => get_ao_density_hf
      procedure :: get_n_hf_parameters => get_n_hf_parameters_hf
      procedure :: get_n_hf_equations  => get_n_hf_equations_hf
!
      procedure :: set_ao_density      => set_ao_density_hf
!
!     Initialize and destruct routines
!
      procedure :: initialize_ao_density       => initialize_ao_density_hf
      procedure :: initialize_ao_fock          => initialize_ao_fock_hf
      procedure :: initialize_mo_coefficients  => initialize_mo_coefficients_hf
      procedure :: initialize_ao_overlap       => initialize_ao_overlap_hf
      procedure :: initialize_orbital_energies => initialize_orbital_energies_hf
      procedure :: initialize_g_wxyz           => initialize_g_wxyz_hf
!
      procedure :: destruct_ao_density       => destruct_ao_density_hf
      procedure :: destruct_ao_fock          => destruct_ao_fock_hf
      procedure :: destruct_mo_coefficients  => destruct_mo_coefficients_hf
      procedure :: destruct_ao_overlap       => destruct_ao_overlap_hf
      procedure :: destruct_orbital_energies => destruct_orbital_energies_hf
      procedure :: destruct_g_wxyz           => destruct_g_wxyz_hf
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
      implicit none
!
      class(hf) :: wf
!
      wf%name = 'HF'
!
      call wf%system%initialize() ! Initialize molecular system -> Should include SOAD
!
      wf%n_ao = 0
      call get_n_aos(wf%n_ao) ! Should this be a molecular system routine?
!
      wf%n_mo = wf%n_ao
!
      wf%n_o = (wf%system%get_n_electrons())/2 ! We only treat closed shell systems
      wf%n_v = wf%n_mo - wf%n_o
!
!     Initialize libint engines
!
      call initialize_coulomb()
      call initialize_kinetic()
      call initialize_nuclear()
      call initialize_overlap() ! SHOULD THESE BE INITIALIZED IN THE ENGINE ?
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
   subroutine initialize_g_wxyz_hf(wf)
!!
!!    Initialize g_wxyz
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      integer(i15) :: n_s
!
      real(dp), dimension(:,:), allocatable :: g_AB_CD
!
      type(interval) :: A_interval
      type(interval) :: B_interval
      type(interval) :: C_interval
      type(interval) :: D_interval
!
      integer(i15) :: A_s, B_s, C_s, D_s
!
      integer(i15) :: a, b, c, d, ab_full, cd_full, ab_reduced, cd_reduced
!
      real(dp) :: s_constr_g, e_constr_g
      real(dp) :: s_timer, e_timer
      real(dp) :: integral_time, copy_time, alloc_time
!
      call cpu_time(s_constr_g)
!
      call mem%alloc(wf%g_wxyz, (wf%n_ao)**2, (wf%n_ao)**2)
      wf%g_wxyz = zero
!
      n_s = wf%system%get_n_shells()
!
      write(output%unit, *) 'Starting to construct g_xyzw...'
      flush(output%unit)
!
      integral_time = zero
      copy_time = zero
      alloc_time = zero
!
!$omp parallel do &
!$omp private(A_s, B_s, C_s, D_s, g_AB_CD, a, b, c, d, ab_reduced, cd_reduced, ab_full, cd_full, &
!$omp A_interval, B_interval, C_interval, D_interval, s_timer, e_timer)
      do D_s = 1, n_s
         do C_s = 1, n_s
            do B_s = 1, n_s
               do A_s = 1, n_s
!
                  A_interval = wf%system%get_shell_limits(A_s)
                  B_interval = wf%system%get_shell_limits(B_s)
                  C_interval = wf%system%get_shell_limits(C_s)
                  D_interval = wf%system%get_shell_limits(D_s)
!
                  call cpu_time(s_timer)
!
                  call mem%alloc(g_AB_CD, &
                     (A_interval%size)*(B_interval%size), &
                     (C_interval%size)*(D_interval%size))
!
                  call cpu_time(e_timer); alloc_time = alloc_time + e_timer - s_timer
!
                  call cpu_time(s_timer)
!
                  call wf%integrals%get_ao_g_wxyz(g_AB_CD, A_s, B_s, C_s, D_s)
!
                  call cpu_time(e_timer); integral_time = integral_time + e_timer - s_timer
!
                  call cpu_time(s_timer)
!
                  do d = 1, D_interval%size
                     do c = 1, C_interval%size
                        do b = 1, B_interval%size
                           do a = 1, A_interval%size
!
                              ab_reduced = index_two(a, b, A_interval%size)
                              cd_reduced = index_two(c, d, C_interval%size)
!
                              ab_full = index_two(a + A_interval%first - 1, b + B_interval%first - 1, wf%n_ao)
                              cd_full = index_two(c + C_interval%first - 1, d + D_interval%first - 1, wf%n_ao)
!
                              wf%g_wxyz(ab_full, cd_full) = g_AB_CD(ab_reduced, cd_reduced)
!
                           enddo
                        enddo
                     enddo
                  enddo
!
                  call cpu_time(e_timer); copy_time = copy_time + e_timer - s_timer
!
                  call cpu_time(s_timer)
!
                  call mem%dealloc(g_AB_CD, &
                     (A_interval%size)*(B_interval%size), &
                     (C_interval%size)*(D_interval%size))
!
                  call cpu_time(e_timer); alloc_time = alloc_time + e_timer - s_timer
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call cpu_time(e_constr_g)
      write(output%unit, *) 'CPU time to construct g (sec) : ', e_constr_g - s_constr_g
!
      write(output%unit, *) 'Breakdown of CPU time: '
!
      write(output%unit, *) 'Integral time: ', integral_time
      write(output%unit, *) 'Alloc time: ', alloc_time
      write(output%unit, *) 'Copy time: ', copy_time
!
      write(output%unit, *) 'Done with constructing g_xyzw'
      flush(output%unit)
!
   end subroutine initialize_g_wxyz_hf
!
!
   subroutine destruct_g_wxyz_hf(wf)
!!
!!    Initialize g_wxyz
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      call mem%dealloc(wf%g_wxyz, (wf%n_ao)**2, (wf%n_ao)**2)
!
   end subroutine destruct_g_wxyz_hf
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
      integer(i15) :: w, x, y, z, wx, yz, w_red, x_red, y_red, z_red
!
      real(dp), dimension(:,:), allocatable :: h_wx
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
      real(dp) :: deg_12, deg_34, deg_12_34, deg
      real(dp) :: temp, temp1, temp2, temp3, temp4, temp5, temp6
!
      real(dp), dimension(:,:), allocatable :: degeneracy
!
      real(dp) :: start_timer, end_timer, omp_get_wtime
!
      wf%ao_fock = zero
!
!     Loop over permutationally unique shell sets
!
      n_s = wf%system%get_n_shells()
!
      start_timer = omp_get_wtime()
!
      call mem%alloc(degeneracy, n_s**2, n_s**2)
      degeneracy = zero
!
!!$omp parallel do private(s1, s2, s3, s4, s1s2, s3s4, deg_12, deg_34, deg_12_34)
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
!!$omp end parallel do
!
      end_timer = omp_get_wtime()
      write(output%unit, *) 'Wall time used to make deg array (sec.): ', end_timer - start_timer
!
      start_timer = omp_get_wtime()
!
!$omp parallel do &
!$omp private(s1, s2, s3, s4, deg_12, deg_34, deg_12_34, deg, s4_max, temp, s1s2, s3s4, &
!$omp A_interval, B_interval, C_interval, D_interval, w, x, y, z, wx, yz, temp1, temp2, temp3, &
!$omp temp4, temp5, temp6, F1, F2, F3, F4, F5, F6, w_red, x_red, y_red, z_red)
      do s1 = 1, n_s
         do s2 = 1, s1
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
                  A_interval = wf%system%get_shell_limits(s4)
                  B_interval = wf%system%get_shell_limits(s3)
                  C_interval = wf%system%get_shell_limits(s2)
                  D_interval = wf%system%get_shell_limits(s1)
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
                        yz = (wf%n_ao)*(z - 1) + y
!
                        do w = A_interval%first, A_interval%last
                           do x = B_interval%first, B_interval%last
!
                              wx = (wf%n_ao)*(x - 1) + w
!
                              temp = deg*wf%g_wxyz(wx, yz)
!
                              temp1 = (one/two)*temp*wf%ao_density(y, z)
                              temp2 = (one/two)*temp*wf%ao_density(w, x)
                              temp3 = (one/eight)*temp*wf%ao_density(y, x)
                              temp4 = (one/eight)*temp*wf%ao_density(y, w)
                              temp5 = (one/eight)*temp*wf%ao_density(z, x)
                              temp6 = (one/eight)*temp*wf%ao_density(w, z)
!!$omp atomic update
!                               wf%ao_fock(w, x) = wf%ao_fock(w, x) + temp1
! !!$omp atomic update
!                               wf%ao_fock(y, x) = wf%ao_fock(y, x) - temp6
! !!$omp atomic update
!                               wf%ao_fock(y, z) = wf%ao_fock(y, z) + temp2
! !!$omp atomic update
!                               wf%ao_fock(w, z) = wf%ao_fock(w, z) - temp3
! !!$omp atomic update
!                               wf%ao_fock(x, z) = wf%ao_fock(x, z) - temp4
! !!$omp atomic update
!                               wf%ao_fock(w, y) = wf%ao_fock(w, y) - temp5
!
                              w_red = w - A_interval%first + 1
                              x_red = x - B_interval%first + 1
                              y_red = y - C_interval%first + 1
                              z_red = z - D_interval%first + 1
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
      call mem%dealloc(degeneracy, n_s**2, n_s**2)
!
      call symmetric_sum(wf%ao_fock, wf%n_ao)
!
      wf%ao_fock = wf%ao_fock*half
!
      end_timer = omp_get_wtime()
      write(output%unit, *) 'Wall time used for two-electron part of Fock matrix (sec.): ', end_timer - start_timer
!
      call wf%calculate_hf_energy(wf%ao_fock)
!
      call mem%alloc(h_wx, wf%n_ao, wf%n_ao)
!
      call get_ao_h_xy(h_wx)
!
      wf%ao_fock = wf%ao_fock + h_wx
!
      call mem%dealloc(h_wx, wf%n_ao, wf%n_ao)
!
   end subroutine construct_ao_fock_hf
!
!
   subroutine calculate_hf_energy_hf(wf, half_GD_wx)
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
   !  real(dp), dimension(:,:), allocatable :: GD_wx
      real(dp), dimension(wf%n_ao, wf%n_ao) :: half_GD_wx
!
      real(dp), dimension(:,:), allocatable :: h_wx
      real(dp), dimension(:,:), allocatable :: g_wyzx
!
!     Construct G(D)
!
!     G(D)_αβ = 2 sum_γδ g_αβγδ D_γδ
!
!GD_wx = two*GD_wx
!       call mem%alloc(GD_wx, wf%n_ao, wf%n_ao)
! !
!       call dgemm('N', 'N',       &
!                   (wf%n_ao)**2,  &
!                   1,             &
!                   (wf%n_ao)**2,  &
!                   two,           &
!                   wf%g_wxyz,     & ! g_αβ_γδ
!                   (wf%n_ao)**2,  &
!                   wf%ao_density, & ! D_γδ
!                   (wf%n_ao)**2,  &
!                   zero,          &
!                   GD_wx,         & ! G(D)_αβ
!                   (wf%n_ao)**2)
! !
! !     G(D)_αβ =+ (-1) sum_γδ g_αδγβ D_γδ
! !
!       call mem%alloc(g_wyzx, (wf%n_ao)**2, (wf%n_ao)**2)
!       g_wyzx = zero
! !
!       call sort_1234_to_1432(wf%g_wxyz, g_wyzx, wf%n_ao, wf%n_ao, wf%n_ao, wf%n_ao)
! !
!       call dgemm('N', 'N',       &
!                   (wf%n_ao)**2,  &
!                   1,             &
!                   (wf%n_ao)**2,  &
!                   -one,          &
!                   g_wyzx,        & ! g_αβ_γδ = g_αδγβ
!                   (wf%n_ao)**2,  &
!                   wf%ao_density, & ! D_γδ
!                   (wf%n_ao)**2,  &
!                   one,           &
!                   GD_wx,         & ! G(D)_αβ
!                   (wf%n_ao)**2)
! !
!       call mem%dealloc(g_wyzx, (wf%n_ao)**2, (wf%n_ao)**2)
!
      call mem%alloc(h_wx, wf%n_ao, wf%n_ao)
      h_wx = zero
!
      call get_ao_h_xy(h_wx)
!
      wf%hf_energy = wf%system%get_nuclear_repulsion()
!
      wf%hf_energy = wf%hf_energy + ddot((wf%n_ao)**2, h_wx, 1, wf%ao_density, 1)
      wf%hf_energy = wf%hf_energy + two*(one/four)*ddot((wf%n_ao)**2, wf%ao_density, 1, half_GD_wx, 1)
!
      call mem%dealloc(h_wx, wf%n_ao, wf%n_ao)
   !   call mem%dealloc(GD_wx, wf%n_ao, wf%n_ao)
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
      wf%ao_overlap = zero
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
      class(hf) :: wf
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
!!    Constructs the MO Fock matrix and returns the occupied-virtual
!!    block of the full matrix.
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_v) :: F ! F_ia
!
      real(dp), dimension(:,:), allocatable :: F_pq
!
      integer(i15) :: i = 0, a = 0
!
      call mem%alloc(F_pq, wf%n_ao, wf%n_ao)
      F_pq = zero
!
      call wf%construct_mo_fock(F_pq)
!
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            F(i,a) = F_pq(i, wf%n_o + a)
!
         enddo
      enddo
!
      call mem%dealloc(F_pq, wf%n_ao, wf%n_ao)
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
      class(hf) :: wf
!
      real(dp), dimension(:,:) :: D
!
      call packin(D, wf%ao_density, wf%n_ao)
!
   end subroutine get_ao_density_hf
!
!
end module hf_class
