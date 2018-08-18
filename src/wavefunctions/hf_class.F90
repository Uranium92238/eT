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
   use array_analysis
   use interval_class
   use index
!
   implicit none
!
   type, extends(wavefunction):: hf
!
      integer(i15) :: n_so ! Number of linearly independent AOs
!
      real(dp) :: hf_energy
!
      real(dp), dimension(:,:), allocatable :: ao_density
      real(dp), dimension(:,:), allocatable :: ao_overlap
      real(dp), dimension(:,:), allocatable :: ao_fock
      real(dp), dimension(:,:), allocatable :: orbital_energies
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
!
      procedure :: construct_and_add => construct_and_add_hf
!
      procedure :: construct_ao_fock_SAD=> construct_ao_fock_SAD_hf
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
      procedure :: solve_roothan_hall => solve_roothan_hall_hf ! remove soon, does not belong to wavefunction 
!
!     Get and set routines for wavefunction variables
!
      procedure :: get_ao_density => get_ao_density_hf
      procedure :: get_fock_ov    => get_fock_ov_hf
!
      procedure :: set_ao_density        => set_ao_density_hf
      procedure :: set_ao_density_to_sad => set_ao_density_to_sad_hf
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
!
      procedure :: construct_roothan_hall_hessian              => construct_roothan_hall_hessian_hf
      procedure :: construct_roothan_hall_gradient             => construct_roothan_hall_gradient_hf
!
      procedure :: construct_sp_eri_schwarz => construct_sp_eri_schwarz_hf
      procedure :: determine_degeneracy     => determine_degeneracy_hf
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
!!    Density is packed
!!
      implicit none
!
      class(hf) :: wf
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
   subroutine construct_sp_eri_schwarz_hf(wf, sp_eri_schwarz, sp_eri_schwarz_list, n_s)
!!
!!    Construct shell-pair electronic-repulsion-integral Schwarz vector
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Computes a vector that contains the largest value (in absolute terms)
!!    of g_wxwx^1/2 for each shell pair (A,B), where w and x is in A and B, 
!!    respectively. These values are used to construct the AO Fock matrix
!!    without calculating integrals that are not needed.
!!
      implicit none
!
      class(hf) :: wf 
!
      integer(i15) :: n_s
!
      real(dp), dimension(n_s*(n_s + 1)/2, 1) :: sp_eri_schwarz
!
      integer(i15), dimension(n_s*(n_s + 1)/2, 3) :: sp_eri_schwarz_list
      integer(i15), dimension(:, :), allocatable  :: sp_eri_schwarz_list_copy
!
      integer(i15), dimension(:,:), allocatable :: sp_eri_schwarz_index_list
      real(dp), dimension(:,:),     allocatable :: sorted_sp_eri_schwarz
!
!     Local variables
!
      integer(i15) :: s1, s2, s1s2
!
      real(dp) :: maximum
!
      real(dp), dimension(:,:), allocatable :: g
!
      type(interval) :: A_interval, B_interval
!
!     Set the maximum element in each shell pair 
!
!$omp parallel do private(s1, s2, s1s2, A_interval, B_interval, g, maximum) schedule(dynamic)
      do s1 = 1, n_s
         do s2 = 1, s1
!
            s1s2 = (max(s1,s2)*(max(s1,s2)-3)/2) + s1 + s2
!
            A_interval = wf%system%shell_limits(s1)
            B_interval = wf%system%shell_limits(s2)
!
            call mem%alloc(g, (A_interval%size)*(B_interval%size), &
                              (A_interval%size)*(B_interval%size))
!
            call wf%system%ao_integrals%get_ao_g_wxyz(g, s1, s2, s1, s2)
!
            maximum = get_abs_max(g, (A_interval%size)*(B_interval%size)**2)
!
            call mem%dealloc(g, (A_interval%size)*(B_interval%size), &
                                (A_interval%size)*(B_interval%size))
!
            sp_eri_schwarz(s1s2, 1) = sqrt(maximum)
!
            sp_eri_schwarz_list(s1s2, 1) = s1 
            sp_eri_schwarz_list(s1s2, 2) = s2
!
         enddo
      enddo
!$omp end parallel do
!
!     Sort the sp_eri_schwarz vector and use the resulting index list 
!     to resort the sp_eri_schwarz_list matrix 
!
      call mem%alloc_int(sp_eri_schwarz_index_list, n_s*(n_s + 1)/2, 1)
      sp_eri_schwarz_index_list = 0
!
      call mem%alloc(sorted_sp_eri_schwarz, n_s*(n_s + 1)/2, 1)
      sorted_sp_eri_schwarz = sp_eri_schwarz
!
      call get_n_highest(n_s*(n_s + 1)/2, n_s*(n_s + 1)/2, sp_eri_schwarz, sorted_sp_eri_schwarz, sp_eri_schwarz_index_list)
!
      sp_eri_schwarz = sorted_sp_eri_schwarz
      call mem%dealloc(sorted_sp_eri_schwarz, n_s*(n_s + 1)/2, 1)
!
      call mem%alloc_int(sp_eri_schwarz_list_copy, n_s*(n_s + 1)/2, 2)
      sp_eri_schwarz_list_copy = sp_eri_schwarz_list
!
!$omp parallel do private(s1s2)
      do s1s2 = 1, n_s*(n_s + 1)/2
!
         sp_eri_schwarz_list(s1s2, :) = sp_eri_schwarz_list_copy(sp_eri_schwarz_index_list(s1s2, 1), :)
!
      enddo
!$omp end parallel do
!
      sp_eri_schwarz_list(:,3) = sp_eri_schwarz_index_list(:,1)
      call mem%dealloc_int(sp_eri_schwarz_index_list, n_s*(n_s + 1)/2, 1)
!
      call mem%dealloc_int(sp_eri_schwarz_list_copy, n_s*(n_s + 1)/2, 2)
!
   end subroutine construct_sp_eri_schwarz_hf
!
!
   subroutine determine_degeneracy_hf(wf, degeneracy, n_s)
!!
!!    Determine degeneracy vector
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    In the AO Fock construction, each shell quadruple has a degeneracy
!!    given by how many times the same integral appears in the actual
!!    unrestricted sum. The value degeneracy(s1s2,s3s4) gives the 
!!    degeneracy of g_wxyz, where w, x, y, z belong to shells s1, s2,
!!    s3 and s4, respectively.
!!
      implicit none
!
      class(hf) :: wf
!
      integer(i15) :: n_s
!
      real(dp), dimension(n_s**2, n_s**2) :: degeneracy
!
      integer(i15) :: s1, s2, s3, s4, s1s2, s3s4, deg_12, deg_34, deg_12_34, s4_max
!
!$omp parallel do private(s1, s2, s3, s4, s1s2, s3s4, deg_12, deg_34, deg_12_34) schedule(dynamic)
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
   subroutine construct_ao_fock_SAD_hf(wf)
!!
!!    Construct AO Fock matrix
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Calculates
!!
!!       F_αβ = h_αβ + sum_γ g_αβγγ D_γγ - 1/2 * sum_γ g_αγγβ D_γγ,
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
      real(dp), dimension(:,:), allocatable :: sp_eri_schwarz, sp_density_schwarz
!
      real(dp), dimension(:,:), allocatable :: h_wx
      integer(i15) :: x, y, z, xy, yz, xz, zz
!
      integer(i15) :: A, B, C
!
      type(interval) :: A_interval
      type(interval) :: B_interval
      type(interval) :: C_interval
!
      real(dp) :: maximum
!
      real(dp), dimension(:,:), allocatable :: g, g_C, g_K, D
!
      n_s = wf%system%get_n_shells()
!
     call mem%alloc(sp_density_schwarz, n_s, 1)
!
!$omp parallel do private(A, A_interval, D, maximum) schedule(dynamic)
      do A = 1, n_s
!
            A_interval = wf%system%shell_limits(A)
!
            call mem%alloc(D, (A_interval%size), (A_interval%size))
!
            D = wf%ao_density(A_interval%first : A_interval%last, A_interval%first : A_interval%last)
!
            maximum = get_abs_max(D, (A_interval%size)*(A_interval%size))
!
            call mem%dealloc(D, (A_interval%size), (A_interval%size))
!
            sp_density_schwarz(A, 1) = maximum
!
      enddo
!$omp end parallel do
!
   call mem%alloc(sp_eri_schwarz, n_s, n_s)
!
!$omp parallel do private(A, B, A_interval, B_interval, g, maximum) schedule(dynamic)
      do A = 1, n_s
         do B = 1, A
!
            A_interval = wf%system%shell_limits(A)
            B_interval = wf%system%shell_limits(B)
!
            call mem%alloc(g, (A_interval%size)*(B_interval%size), &
                              (A_interval%size)*(B_interval%size))
!
            call wf%system%ao_integrals%get_ao_g_wxyz(g, A, B, A, B)
!
            maximum = get_abs_max(g, (A_interval%size)*(B_interval%size)**2)
!
            call mem%dealloc(g, (A_interval%size)*(B_interval%size), &
                                (A_interval%size)*(B_interval%size))
!
            sp_eri_schwarz(A, B) = sqrt(maximum)
            sp_eri_schwarz(B, A) = sqrt(maximum)
!
         enddo
      enddo
!$omp end parallel do
!
      wf%ao_fock = zero
!
!$omp parallel do &
!$omp private(A, B, C, A_interval, B_interval, C_interval, x, y, z, xy, zz, xz, yz, &
!$omp g_C, g_K) schedule(dynamic)
      do A = 1, n_s
!
         A_interval = wf%system%shell_limits(A)
!
         do B = 1, A
!
            B_interval = wf%system%shell_limits(B)
!           
            do C = 1, n_s
!
               if (sp_eri_schwarz(A, B)*sp_eri_schwarz(C, C)*sp_density_schwarz(C, 1) .lt. 1.0d-10 .and. &
                  sp_eri_schwarz(A, C)*sp_eri_schwarz(B, C)*sp_density_schwarz(C, 1) .lt. 1.0d-8) continue
!
               C_interval = wf%system%shell_limits(C)
!
               call mem%alloc(g_C, (A_interval%size)*(B_interval%size), &
                                 (C_interval%size)*(C_interval%size))
               call mem%alloc(g_K, (A_interval%size)*(C_interval%size), &
                                 (B_interval%size)*(C_interval%size))
!
               call wf%system%ao_integrals%get_ao_g_wxyz(g_C, A, B, C, C)
               call wf%system%ao_integrals%get_ao_g_wxyz(g_K, A, C, B, C)
!
!              Add Fock matrix contributions
!
               if (A .ne. B) then
!
                  do x = A_interval%first, A_interval%last
                     do y = B_interval%first, B_interval%last
!
                        xy = A_interval%size*(y - B_interval%first) + x - A_interval%first + 1
!
                        do z = C_interval%first, C_interval%last
!
                           zz = C_interval%size*(z - C_interval%first) + z  - C_interval%first + 1
                           xz = A_interval%size*(z - C_interval%first) + x  - A_interval%first + 1
                           yz = B_interval%size*(z - C_interval%first) + y  - B_interval%first + 1
!
                           wf%ao_fock(x, y) = wf%ao_fock(x, y) + (g_C(xy, zz)- half*g_K(xz, yz))*wf%ao_density(z, z)
!
                        enddo
                     enddo
                  enddo
!
               else
!
                  do x = A_interval%first, A_interval%last
                     do y = A_interval%first, x
!
                        xy = A_interval%size*(y - A_interval%first) + x - A_interval%first + 1
!
                        do z = C_interval%first, C_interval%last
!
                           zz = C_interval%size*(z - C_interval%first) + z  - C_interval%first + 1
                           xz = A_interval%size*(z - C_interval%first) + x  - A_interval%first + 1
                           yz = B_interval%size*(z - C_interval%first) + y  - B_interval%first + 1
!
                           wf%ao_fock(x, y) = wf%ao_fock(x, y) + (g_C(xy, zz)- half*g_K(xz, yz))*wf%ao_density(z, z)
!
                        enddo
                     enddo
                  enddo
               endif
!                  
               call mem%dealloc(g_C, (A_interval%size)*(B_interval%size), &
                                 (C_interval%size)*(C_interval%size))
               call mem%dealloc(g_K, (A_interval%size)*(C_interval%size), &
                                 (B_interval%size)*(C_interval%size))
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
!$omp parallel do private(x) schedule(static)
      do x = 1, wf%n_ao
!
         wf%ao_fock(x, x) = half*wf%ao_fock(x, x)
!
      enddo
!$omp end parallel do
!
      call mem%dealloc(sp_density_schwarz, n_s, 1)
      call mem%dealloc(sp_eri_schwarz, n_s, n_s)
!
       call symmetric_sum(wf%ao_fock, wf%n_ao)
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
   end subroutine construct_ao_fock_SAD_hf
!
!
   subroutine construct_ao_fock_hf(wf, sp_eri_schwarz, sp_eri_schwarz_list, n_s, coulomb, exchange)
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
!!    Try this today: sort the integral screening vector (s1s2|s1s2) 
!!    according to size, together with an index list. Ask only to sort
!!    the elements above the threshold and ignore the rest! Then loop 
!!    over these s1 and s2. We should make the vector packed first, then 
!!    calculate the appropriate packed index in the loop s1 >= s2 when 
!!    actually constructing the Fock matrix.
!!
      implicit none
!
      class(hf) :: wf
!
      integer(i15) :: n_s
!
      real(dp), dimension(n_s*(n_s + 1)/2, 1)     :: sp_eri_schwarz
      integer(i15), dimension(n_s*(n_s + 1)/2, 3) :: sp_eri_schwarz_list ! list(s1s2, 1) = s1, list(s1s2, 2) = s2, list(s1s2, 3) = s1s2_sorted
!
      real(dp), optional :: coulomb, exchange ! Non-standard thresholds
!
      real(dp) :: coulomb_thr, exchange_thr ! Actual thresholds 
!
      real(dp), dimension(:,:), allocatable :: ao_fock_packed
      real(dp), dimension(:,:), allocatable :: X_wz, h_wx, h_wx_square
      integer(i15) :: w, x, y, z, wx, yz, w_red, x_red, y_red, z_red
!
      integer(i15) :: n_threads, omp_get_max_threads, thread_offset, thread, omp_get_thread_num
      real(dp), dimension(:,:), allocatable :: F 
!
      integer(i15) :: s1, s2, s3, s4, s4_max, s1s2, s3s4, s3s4_sorted
!
      type(interval) :: A_interval
      type(interval) :: B_interval
      type(interval) :: C_interval
      type(interval) :: D_interval
!
      logical :: skip
!
      integer(i15) :: n_sig_sp
!
      real(dp) :: deg_12, deg_34, deg_12_34, deg, ddot, norm
      real(dp) :: temp, temp1, temp2, temp3, temp4, temp5, temp6
      real(dp) :: maximum, max_D_schwarz, max_eri_schwarz
!
      real(dp), dimension(:,:), allocatable :: g, D, sp_density_schwarz
!
      type(interval), dimension(:), allocatable :: shell_limits 
!
!     Set thresholds to ignore Coulomb and exchange terms 
!
      if (present(coulomb)) then 
!
         coulomb_thr = coulomb 
!
      else
!
         coulomb_thr = 1.0D-10 
!
      endif 
!
      if (present(exchange)) then 
!
         exchange_thr = exchange 
!
      else
!
         exchange_thr = 1.0D-8
!
      endif 
!
!     Allocate and set shell limits vector 
!
   !   allocate(shell_limits(n_s))
   !   do s1 = 1, n_s 
!  
   !      shell_limits(s1) = wf%system%shell_limits(s1)
!
    !  enddo
!
      call mem%alloc(sp_density_schwarz, n_s, n_s)
!
!     Compute number of significant shell pairs (pre-screening)
!
      n_sig_sp = 0
      do s1s2 = 1, n_s*(n_s + 1)/2
!
         if (sp_eri_schwarz(s1s2, 1)**2 .lt. coulomb_thr) then
!
            exit
!
         else
!
            n_sig_sp = n_sig_sp + 1
!
         endif
!
      enddo
!
      write(output%unit, *) 'Number of shell pairs:', n_s*(n_s + 1)/2
      write(output%unit, *) 'Number of significant shell pairs:', n_sig_sp
!
!$omp parallel do private(s1, s2, A_interval, B_interval, D, maximum) schedule(dynamic)
      do s1 = 1, n_s
         do s2 = 1, s1
!
            A_interval = wf%system%shell_limits(s1)
            B_interval = wf%system%shell_limits(s2)
!
            call mem%alloc(D, (A_interval%size), (B_interval%size))
!
            D = wf%ao_density(A_interval%first : A_interval%last, B_interval%first : B_interval%last)
!
            maximum = get_abs_max(D, (A_interval%size)*(B_interval%size))
!
            call mem%dealloc(D, (A_interval%size), (B_interval%size))
!
            sp_density_schwarz(s1, s2) = maximum
            sp_density_schwarz(s2, s1) = maximum
!
         enddo
      enddo
!$omp end parallel do
!
!     Calculate maximum of all the shell pair maximums prescreening
!
      max_D_schwarz     = get_abs_max(sp_density_schwarz, n_s**2)
      max_eri_schwarz   = get_abs_max(sp_eri_schwarz, n_s*(n_s + 1)/2)
!
      n_threads = omp_get_max_threads()
      call mem%alloc(F, wf%n_ao, wf%n_ao*n_threads) ! [F(thr1) F(thr2) ...]
      F = zero 
!
      call wf%construct_and_add(F, n_threads, max_D_schwarz, max_eri_schwarz, & 
                                 sp_density_schwarz, sp_eri_schwarz, sp_eri_schwarz_list, &
                                 n_s, n_sig_sp, coulomb_thr, exchange_thr, wf%system%shell_limits)
!
      write(output%unit, *) 'Number of threads:', n_threads
!
! !$omp parallel do &
! !$omp private(s1, s2, s3, s4, deg, s4_max, temp, s1s2, s3s4, s3s4_sorted, deg_12, deg_34, deg_12_34, thread_offset, &
! !$omp A_interval, B_interval, C_interval, D_interval, w, x, y, z, wx, yz, temp1, temp2, temp3, &
! !$omp temp4, temp5, temp6, w_red, x_red, y_red, z_red, g, skip) schedule(dynamic)
!       do s1s2 = 1, n_sig_sp
! !
!          thread_offset = omp_get_thread_num()*wf%n_ao ! Start column of thread's Fock matrix 
! !
!          s1 = sp_eri_schwarz_list(s1s2, 1)
!          s2 = sp_eri_schwarz_list(s1s2, 2)
! !
!          A_interval = wf%system%shell_limits(s1)
!          B_interval = wf%system%shell_limits(s2)
! !
!          if (sp_eri_schwarz(s1s2, 1)*(max_D_schwarz)*(max_eri_schwarz) .lt. coulomb_thr) continue
! !
! !        s1 => s2, so s2/s1 = 0 when unequal, 1 when equal 
! !
! !        deg_12 = real(2-s1/s2)
! !
!          deg_12 = real(2-s2/s1, kind=dp)
! !
!          do s3 = 1, s1
! !
!             C_interval = wf%system%shell_limits(s3)
! !
! !           s3 <= s1, so s3/s1 = 0 if unequal, 1 if not 
! !
! !           s4_max = (s3/s1)*s3 + (1-s3/s1)*s2 
! !
!             s4_max = (s3/s1)*s2 + (1-s3/s1)*s3
! !
!             do s4 = 1, s4_max
! !
!                s3s4 = (max(s3,s4)*(max(s3,s4)-3)/2) + s3 + s4 
!                s3s4_sorted = sp_eri_schwarz_list(s3s4, 3)
!                temp = sp_eri_schwarz(s1s2, 1)*sp_eri_schwarz(s3s4_sorted, 1)
!                temp = sp_eri_schwarz(s1s2, 1)*max_eri_schwarz ! should be refined 
! !
!                if (temp*(max_D_schwarz) .lt. coulomb_thr) continue
! !
!                if (temp*sp_density_schwarz(s3,s4) .lt. coulomb_thr   .and. & ! F1
!                    temp*sp_density_schwarz(s1,s2) .lt. coulomb_thr   .and. & ! F2
!                    temp*sp_density_schwarz(s3,s2) .lt. exchange_thr  .and. & ! F3
!                    temp*sp_density_schwarz(s3,s1) .lt. exchange_thr  .and. & ! F4
!                    temp*sp_density_schwarz(s4,s2) .lt. exchange_thr  .and. & ! F5
!                    temp*sp_density_schwarz(s1,s4) .lt. exchange_thr) continue ! F6
! !
!                deg_34    = real(2-s4/s3, kind=dp)
!                deg_12_34 = min(1-s3/s1+2-min(s4/s2,s2/s4), 2)

!                deg = deg_12*deg_34*deg_12_34 ! Shell degeneracy
! !
!                D_interval = wf%system%shell_limits(s4)
! !
!                call mem%alloc(g, (A_interval%size)*(B_interval%size), &
!                                     (C_interval%size)*(D_interval%size))
! !
!                call wf%system%ao_integrals%get_ao_g_wxyz(g, s1, s2, s3, s4)
!                g = deg*g
! !
! !              Add Fock matrix contributions
! !
!                do x = B_interval%first, B_interval%last
!                   do w = A_interval%first, A_interval%last
!                      do z = D_interval%first, D_interval%last
!                         do y = C_interval%first, C_interval%last
! !
!                            y_red = y - C_interval%first + 1
!                            z_red = z - D_interval%first + 1
! !
!                            yz = (C_interval%size)*(z_red - 1) + y_red
! !
!                            w_red = w - A_interval%first + 1
!                            x_red = x - B_interval%first + 1
! !
!                            wx = (A_interval%size)*(x_red - 1) + w_red
! !
!                            temp = g(wx, yz)
! !
!                            temp1 = half*temp*wf%ao_density(y, z)
!                            temp2 = half*temp*wf%ao_density(w, x)
!                            temp3 = one_over_eight*temp*wf%ao_density(y, x)
!                            temp4 = one_over_eight*temp*wf%ao_density(y, w)
!                            temp5 = one_over_eight*temp*wf%ao_density(z, x)
!                            temp6 = one_over_eight*temp*wf%ao_density(w, z)
! !
!                            F(w, thread_offset + x) = F(w, thread_offset + x) + temp1
!                            F(y, thread_offset + z) = F(y, thread_offset + z) + temp2
!                            F(w, thread_offset + z) = F(w, thread_offset + z) - temp3
!                            F(x, thread_offset + z) = F(x, thread_offset + z) - temp4
!                            F(w, thread_offset + y) = F(w, thread_offset + y) - temp5
!                            F(y, thread_offset + x) = F(y, thread_offset + x) - temp6
! !
!                         enddo
!                      enddo
!                   enddo
!                enddo
! !
!                call mem%dealloc(g, (A_interval%size)*(B_interval%size), &
!                                     (C_interval%size)*(D_interval%size))
! !
!             enddo
!          enddo
!       enddo
! !$omp end parallel do
!
!
      call mem%dealloc(sp_density_schwarz, n_s, n_s)
    !  deallocate(shell_limits)
!
!     Put the accumulated Fock matrices from each thread into the Fock matrix 
!
      wf%ao_fock = zero
!
      do thread = 1, n_threads
!
         call daxpy(wf%n_ao**2, one, F(1, (thread-1)*wf%n_ao + 1), 1, wf%ao_fock, 1)
!
      enddo
!
      call mem%dealloc(F, wf%n_ao, wf%n_ao*n_threads) ! [F(thr1) F(thr2) ...]
!
      call symmetric_sum(wf%ao_fock, wf%n_ao)
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
   end subroutine construct_ao_fock_hf
!
!
   subroutine construct_and_add_hf(wf, F, n_threads, max_D_schwarz, max_eri_schwarz, & 
                                 sp_density_schwarz, sp_eri_schwarz, sp_eri_schwarz_list, &
                                 n_s, n_sig_sp, coulomb_thr, exchange_thr, shells)
!
      implicit none 
!
      class(hf), intent(in) :: wf 
!
      integer(i15), intent(in) :: n_threads 
      integer(i15), intent(in) :: n_s 
      integer(i15), intent(in) :: n_sig_sp
!
      type(interval), dimension(n_s), intent(in) :: shells
!
      real(dp), intent(in) :: max_D_schwarz, max_eri_schwarz
!
      real(dp), intent(in) :: coulomb_thr, exchange_thr
!
      real(dp), dimension(wf%n_ao, wf%n_ao*n_threads) :: F 
!
      real(dp), dimension(n_s*(n_s + 1)/2, 1), intent(in)     :: sp_eri_schwarz
      integer(i15), dimension(n_s*(n_s + 1)/2, 3), intent(in) :: sp_eri_schwarz_list
!
      real(dp), dimension(n_s, n_s), intent(in)               :: sp_density_schwarz
!
      real(dp) :: temp, temp1, temp2, temp3, temp4, temp5, temp6
      real(dp) :: deg, deg_12, deg_34, deg_12_34
!
      integer(i15) :: w, x, y, omp_get_thread_num, z, wx, yz, s1s2, s1, s2, s3, s4, s4_max 
      integer(i15) :: s3s4, s3s4_sorted, w_red, x_red, y_red, z_red, thread_offset, wxyz 
!
      type(interval) :: A_interval, B_interval, C_interval, D_interval
!
      real(dp) :: d1, d2, d3, d4, d5, d6 
!
      real(dp), dimension(:,:), allocatable :: g 
!
      integer(i15) :: max_shell_size, dim_s1, dim_s2, dim_s3, dim_s4 
      integer(i15) :: s1_first, s2_first, s3_first, s4_first 
      integer(i15) :: s1_last, s2_last, s3_last, s4_last
!
!     Determine largest shell size 
!
      max_shell_size = 0
      do s1 = 1, n_s 
!
         if (max_shell_size .lt. shells(s1)%size) max_shell_size = shells(s1)%size
!
      enddo
!
      call mem%alloc(g, max_shell_size**4, 1)
!
!
!$omp parallel do &
!$omp private(s1, s2, s3, s4, deg, s4_max, temp, s1s2, s3s4, s3s4_sorted, deg_12, deg_34, deg_12_34, thread_offset, &
!$omp w, x, y, z, wx, yz, temp1, temp2, temp3, d1, d2, d3, d4, d5, d6, s1_first, s2_first, s3_first, s4_first, &
!$omp temp4, temp5, temp6, w_red, x_red, dim_s1, dim_s2, dim_s3, dim_s4, y_red, z_red, wxyz, g, &
!$omp s1_last, s2_last, s3_last, s4_last) schedule(dynamic)
      do s1s2 = 1, n_sig_sp
!
         thread_offset = omp_get_thread_num()*wf%n_ao ! Start column of thread's Fock matrix 
!
         s1 = sp_eri_schwarz_list(s1s2, 1)
         s2 = sp_eri_schwarz_list(s1s2, 2)
!
      !   A_interval = wf%system%shell_limits(s1)
      !   B_interval = wf%system%shell_limits(s2)
!
         if (sp_eri_schwarz(s1s2, 1)*(max_D_schwarz)*(max_eri_schwarz) .lt. coulomb_thr) continue
!
!        s1 => s2, so s2/s1 = 0 when unequal, 1 when equal 
!
!        deg_12 = real(2-s1/s2)
!
       !  deg_12 = real(2-s2/s1, kind=dp)
         deg_12 = 2-s2/s1
!
               dim_s1 = shells(s1)%size
               dim_s2 = shells(s2)%size

!
               s1_first = shells(s1)%first
               s2_first = shells(s2)%first
!
               s1_last = shells(s1)%last
               s2_last = shells(s2)%last
!
         do s3 = 1, s1
!

               dim_s3 = shells(s3)%size

               s3_first = shells(s3)%first

               s3_last = shells(s3)%last
!
        !    C_interval = wf%system%shell_limits(s3)
!
!           s3 <= s1, so s3/s1 = 0 if unequal, 1 if not 
!
!           s4_max = (s3/s1)*s3 + (1-s3/s1)*s2 
!
            s4_max = (s3/s1)*s2 + (1-s3/s1)*s3
!
            do s4 = 1, s4_max
!
               s3s4 = (max(s3,s4)*(max(s3,s4)-3)/2) + s3 + s4 
               s3s4_sorted = sp_eri_schwarz_list(s3s4, 3)
               temp = sp_eri_schwarz(s1s2, 1)*sp_eri_schwarz(s3s4_sorted, 1)
               temp = sp_eri_schwarz(s1s2, 1)*max_eri_schwarz ! should be refined 
!
               if (temp*(max_D_schwarz) .lt. coulomb_thr) continue
!
               if (temp*sp_density_schwarz(s3,s4) .lt. coulomb_thr   .and. & ! F1
                   temp*sp_density_schwarz(s1,s2) .lt. coulomb_thr   .and. & ! F2
                   temp*sp_density_schwarz(s3,s2) .lt. exchange_thr  .and. & ! F3
                   temp*sp_density_schwarz(s3,s1) .lt. exchange_thr  .and. & ! F4
                   temp*sp_density_schwarz(s4,s2) .lt. exchange_thr  .and. & ! F5
                   temp*sp_density_schwarz(s1,s4) .lt. exchange_thr) continue ! F6
!
            !   deg_34    = real(2-s4/s3, kind=dp)
               deg_34    = 2-s4/s3
               deg_12_34 = min(1-s3/s1+2-min(s4/s2,s2/s4), 2)

               deg = deg_12*deg_34*deg_12_34 ! Shell degeneracy
!
           !    D_interval = wf%system%shell_limits(s4)
!
             !  call mem%alloc(g, (shells(s1)%size)*(shells(s2)%size), &
             !                       (shells(s3)%size)*(shells(s4)%size))
!
               dim_s4 = shells(s4)%size
!
               s4_first = shells(s4)%first
!
               s4_last = shells(s4)%last
!
               call wf%system%ao_integrals%get_ao_g_wxyz(g, s1, s2, s3, s4)
               g(1:dim_s1*dim_s2*dim_s3*dim_s4,1) = deg*g(1:dim_s1*dim_s2*dim_s3*dim_s4,1)
! w x y z 
!
!              Add Fock matrix contributions
!
               do z = s4_first, s4_last
!
                  z_red = z - s4_first + 1
!
                  do y = s3_first, s3_last 
!
                     y_red = y - s3_first + 1
!
                     d1 = wf%ao_density(y, z)
!
                     do x = s2_first, s2_last 
!
                        x_red = x - s2_first + 1
!
                        ! d3 = wf%ao_density(y, x)
                        ! d5 = wf%ao_density(z, x)
                         d3 = wf%ao_density(x, y)
                         d5 = wf%ao_density(x, z)
!
                        do w = s1_first, s1_last 
!
                           d2 = wf%ao_density(w, x)
                         !  d4 = wf%ao_density(y, w)
                           d4 = wf%ao_density(w, y)
                           d6 = wf%ao_density(w, z)
!
                           w_red = w - s1_first + 1
!
                       !    wxyz = dim_wx*(yz-1)+wx
                           !xyz = dim_s2*(dim_s3*(z_red-1)+y_red-1)+x_red
                           wxyz = dim_s1*(dim_s2*(dim_s3*(z_red-1)+y_red-1)+x_red-1)+w_red
!
                         !  temp = g(wx, yz)
                           temp = g(wxyz, 1)
!
                           temp1 = half*temp*d1
                           temp2 = half*temp*d2
                           temp3 = one_over_eight*temp*d3
                           temp4 = one_over_eight*temp*d4
                           temp5 = one_over_eight*temp*d5
                           temp6 = one_over_eight*temp*d6
!
                           F(w, thread_offset + x) = F(w, thread_offset + x) + temp1
                           F(y, thread_offset + z) = F(y, thread_offset + z) + temp2
                           F(w, thread_offset + z) = F(w, thread_offset + z) - temp3
                           F(x, thread_offset + z) = F(x, thread_offset + z) - temp4
                           F(w, thread_offset + y) = F(w, thread_offset + y) - temp5
                           F(y, thread_offset + x) = F(y, thread_offset + x) - temp6
!
                        enddo
                     enddo
                  enddo
               enddo
!
          !     call mem%dealloc(g, (shells(s1)%size)*(shells(s2)%size), &
          !                          (shells(s3)%size)*(shells(s4)%size))
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(g, max_shell_size**4, 1)
!
   end subroutine construct_and_add_hf
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
   subroutine set_ao_density_to_sad_hf(wf)
!!
!!    Set AO density to superposition of atomic densities (SAD) guess
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Careful! The AO density is allocated and zeroed by the initialize
!!    AO density routine and to avoid zeroing twice, this routine assumes
!!    that the non-diagonal terms of the density are already zero.
!!
!!    If called at a time when the density differs from zero, remember
!!    to zero it before calling this routine.
!!
      implicit none 
!
      class(hf) :: wf 
!
      real(dp), dimension(:,:), allocatable :: density_diagonal
!
      integer(i15) :: ao
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
   end subroutine set_ao_density_to_sad_hf
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
!!    Packs the AO density into D.
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
      integer(kind=4), dimension(:, :), allocatable :: used_diag
!
      real(dp), dimension(:,:), allocatable :: perm_matrix
!
      real(dp), dimension(:,:), allocatable :: csc
      real(dp), dimension(:,:), allocatable :: tmp
!
      integer(i15) :: rank
      integer(i15) :: i, j
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
!!    There is something seriously wrong with the way we handle 
!!    integers. For H2O/cc-pVDZ, used_diag becomes nonsense (large
!!    integers instead of pivots, thus giving a seg. fault)!
!!
!!    First case where eT version of BLAS routine does not work...
!!
      implicit none
!
      class(hf) :: wf
!
      integer(kind=4), dimension(:, :), allocatable :: used_diag
!
      real(dp), dimension(wf%n_ao, wf%n_ao) :: L
!
      real(dp), dimension(:,:), allocatable :: perm_matrix
!
      real(dp), dimension(:,:), allocatable :: csc
      real(dp), dimension(:,:), allocatable :: tmp
!
      integer(i15) :: rank
      integer(i15) :: i, j
!
      allocate(used_diag(wf%n_ao, 1))
      used_diag = 0
!
      L = zero
!
      call full_cholesky_decomposition_system(wf%ao_overlap, L, wf%n_ao, rank, &
                                                      1.0D-32, used_diag)
!
      ! write(output%unit, *) 'n_ao:', wf%n_ao 
      ! write(output%unit, *) 'rank:', rank
      ! write(output%unit, *) 'used_diag:', used_diag
      ! flush(output%unit)
!
!     Make permutation matrix P
!
      call mem%alloc(perm_matrix, wf%n_ao, wf%n_ao)
!
      perm_matrix = zero
!
      do j = 1, wf%n_ao
!
         perm_matrix(used_diag(j, 1), j) = one
!
      enddo
!
      deallocate(used_diag)
!
!     Make matrix L
!
      call mem%alloc(tmp, wf%n_ao, wf%n_ao)
!
      call dgemm('N','N',      &
                  wf%n_ao,     &
                  wf%n_ao,     &
                  wf%n_ao,     &
                  one,         &
                  perm_matrix, &
                  wf%n_ao,     &
                  L,           &
                  wf%n_ao,     &
                  zero,        &
                  tmp,         &
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
      wf%ao_density = wf%ao_density + C
!
      call mem%dealloc(M, wf%n_ao, wf%n_ao)
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
!!    Constructs Po = 1/2 D S
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
!     Po = 1/2 D S
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
!!    where Po = 1/2 D S and Pv = 1 - Po. In Po, D is the AO density and S
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
end module hf_class
