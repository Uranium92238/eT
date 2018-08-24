module mlhf_class
!
!!
!!    Multilevel Hartree-Fock (MLHF) class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   use hf_class
!
   use reordering
   use interval_class
   use index
   use active_atoms_info_class
!
   use eri_cd_solver_class
!
   use array_utilities
!
   implicit none
!
   type, extends(hf):: mlhf
!
   contains
!
      procedure :: initialize => initialize_mlhf
      procedure :: finalize => finalize_mlhf
!
      procedure :: construct_virtual_density => construct_virtual_density_mlhf
!
      procedure :: eri_decomp_test_w_active_dens => eri_decomp_test_w_active_dens_mlhf
!
   end type mlhf
!
!
contains
!
!
   subroutine initialize_mlhf(wf)
!!
!!    Initialize
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(mlhf) :: wf
!
      wf%name = 'MLHF'
!
      write(output%unit, '(a)')':: SAD-ML test'
      flush(output%unit)
!
      call wf%system%initialize() ! Initialize molecular system -> Should include SAD
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
      call wf%eri_decomp_test_w_active_dens()
!
   end subroutine initialize_mlhf
!
!
   subroutine finalize_mlhf(wf)
!!
!!    Finalize
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(mlhf) :: wf
!
!     Nothing here yet
!
   end subroutine finalize_mlhf
!
!
   subroutine eri_decomp_test_w_active_dens_mlhf(wf)
!
      implicit none
!
      class(mlhf) :: wf
!
      real(dp), dimension(:,:), allocatable :: cholesky_vectors_occ, cholesky_vectors_virt, V, ao_density_v
!
      integer(i15):: i, j, k, n_active_aos, ao_offset, active_ao_counter, n_vectors_occ, n_vectors_virt
      integer(i15):: a, x
!
      real(dp) :: max, e_construct_fock, s_construct_fock, omp_get_wtime
!
      integer(i15), dimension(:,:), allocatable :: active_aos, sp_eri_schwarz_list
!
      integer(i15) :: ao, n_active_occ, n_active_vir, n_s
!
      real(dp), dimension(:,:), allocatable :: density_diagonal, eri_deg, sp_eri_schwarz
!
      type(eri_cd_solver)  :: chol_solver
!
      call wf%initialize_ao_density()
!
!     Set initial density to superposition of atomic densities (SOAD) guess
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
!     Construct initial AO Fock from the SOAD density
!
      n_s = wf%system%get_n_shells()
!
!     Construct initial AO Fock from the SOAD density
!
      call wf%initialize_ao_fock()
!
      s_construct_fock = omp_get_wtime()
!
      write(output%unit,*)'set up fock'
      flush(output%unit)
!
      call wf%construct_ao_fock_SAD(1.0d-13, 1.0d-13, 1.0d-25)
      !write(output%unit,*)wf%hf_energy
!
      e_construct_fock = omp_get_wtime()
      write(output%unit, '(/a49, f11.2)')'Wall time to construct AO fock from SAD density: ', &
                                  e_construct_fock - s_construct_fock
      flush(output%unit)
!
!     Construct AO overlap matrix, Cholesky decompose it,
!     followed by preconditioning (making it the identity matrix
!     for this particular preconditioner - V)
!
      call wf%initialize_ao_overlap()
      call wf%construct_ao_overlap()
      call wf%decompose_ao_overlap_2() 
!
      write(output%unit, *)'Removed ', wf%n_ao - wf%n_so, 'AOs.'
!
!     Solve Roothan Hall once - using the SOAD guess - to get a decent AO density
!     on  which to start the preconditioned conjugate gradient (PCG) algorithm
!
      call wf%initialize_orbital_energies()
      call wf%initialize_mo_coefficients()
      call wf%do_roothan_hall() ! F^AO C = S C e to get new MOs C
!
!     Update the AO density
!
      call wf%construct_ao_density() ! Construct AO density from C
!
      n_active_aos = 0
!
      do i = 1, wf%system%n_active_atoms
!
         n_active_aos = n_active_aos + wf%system%atoms(i)%n_ao
!
      enddo 
!
      call mem%alloc_int(active_aos, n_active_aos, 1)
!
      do i = 1, n_active_aos
!
           active_aos(i, 1) = i 
!
      enddo
!
!     Determine number of active occupied
!
      n_active_occ = 0
!
      do i = 1, wf%system%n_active_atoms
!
        n_active_occ = n_active_occ + wf%system%atoms(i)%number
!
      enddo
!
      n_active_occ = n_active_occ/2
!
!     Determine number of active virtuals
!
      n_active_vir = n_active_occ * (wf%n_v/wf%n_o)
!
      call mem%alloc(ao_density_v, wf%n_ao, wf%n_ao)
      call wf%construct_virtual_density(ao_density_v)
!
      call mem%alloc(cholesky_vectors_occ, wf%n_ao, n_active_aos)
!
      call dscal(wf%n_ao**2, half, wf%ao_density, 1)
!
      call cholesky_decomposition_limited_diagonal(wf%ao_density, cholesky_vectors_occ, wf%n_ao, &
                                                     n_vectors_occ, 1.0d-2, n_active_aos, active_aos, n_active_occ)
!
      call mem%alloc(cholesky_vectors_virt, wf%n_ao, n_active_aos)
      cholesky_vectors_virt = zero
!
      call cholesky_decomposition_limited_diagonal(ao_density_v, cholesky_vectors_virt, wf%n_ao, &
                                                     n_vectors_virt, 1.0d-2, n_active_aos, active_aos, n_active_vir)
!
      call mem%dealloc_int(active_aos, n_active_aos, 1)
!
      call mem%alloc(V, wf%n_ao, 1)
!
      do x = 1, wf%n_ao
!
         max = 0.0d0
!
         do i = 1, n_vectors_occ
!
            if (cholesky_vectors_occ(x, i)**2 .gt. max) max = cholesky_vectors_occ(x, i)**2     
!
         enddo
!
         do a = 1, n_vectors_virt
!
            if (cholesky_vectors_virt(x, a)**2 .gt. max) max = cholesky_vectors_virt(x, a)**2 
!
         enddo
!
         V(x, 1) = max
!
      enddo
!
      call mem%dealloc(cholesky_vectors_virt, wf%n_ao, n_active_aos)
      call mem%dealloc(cholesky_vectors_occ, wf%n_ao, n_active_aos)
!
      write(output%unit, '(a27, i4)')'Number of occupied active: ', n_vectors_occ
      write(output%unit, '(a27, i4)')'Number of virtual active:  ', n_vectors_virt
      flush(output%unit)
!
!
      call chol_solver%initialize(wf%system)
      call chol_solver%solve(wf%system, V)
      call chol_solver%finalize()
!
!     Cholesky decomposition
!
      call mem%dealloc(V, wf%n_ao, 1)
!
   end subroutine eri_decomp_test_w_active_dens_mlhf
!
!
  subroutine construct_virtual_density_mlhf(wf, D_v)
!!
!!    Construct virtual density
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    D^V = S^-1 - D
!!        = P * L^-T * L^-1 * P^T - D
!!        = (L^-1 * P^T)^T * (L^-1 * P^T) - D
!!
!!    P = pivot matrix
!!    L = cholesky vectors of overlap
!!
      implicit none 
!
      class(mlhf) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao) :: D_v
!
      integer(i15) :: rank, i, x, y
!
      integer(kind=4), dimension(:), allocatable :: piv
!
      real(dp), dimension(:,:), allocatable :: L, L_inv, P,  L_inv_P_trans
!
      rank = 0
!
     allocate(piv(wf%n_ao))
     call mem%alloc(L, wf%n_ao, wf%n_ao)
!
     call full_cholesky_decomposition_system(wf%ao_overlap, L, wf%n_ao, rank, &
                                               1.0d-16, piv)
!
     call mem%alloc(L_inv, wf%n_ao, wf%n_ao)
!
     if (rank .eq. wf%n_ao) then
!
        call inv_lower_tri(L_inv, L, wf%n_ao)
!
     else
!
        write(output%unit) 'Error: does not yet work if S has lin dep'
        stop
!
     endif
!
     call mem%dealloc(L, wf%n_ao, wf%n_ao)
!
     call mem%alloc(P, wf%n_ao, wf%n_ao)
     P = zero
!
     do i = 1, wf%n_ao 
!
        P(piv(i), i) = one
!
     enddo
!
     deallocate(piv)
!
     call mem%alloc(L_inv_P_trans, wf%n_ao, wf%n_ao)
!
     call dgemm('N', 'T',       &
                 wf%n_ao,       &
                 wf%n_ao,       &
                 wf%n_ao,       &
                 one,           &
                 L_inv,         &
                 wf%n_ao,       &
                 P,             &
                 wf%n_ao,       &
                 zero,          &
                 L_inv_P_trans, &
                 wf%n_ao)
!
     call mem%dealloc(L_inv, wf%n_ao, wf%n_ao)
     call mem%dealloc(P, wf%n_ao, wf%n_ao)
!
     call dgemm('T', 'N',       &
                 wf%n_ao,       &
                 wf%n_ao,       &
                 wf%n_ao,       &
                 one,           &
                 L_inv_P_trans, &
                 wf%n_ao,       &
                 L_inv_P_trans, &
                 wf%n_ao,       &
                 zero,          &
                 D_v,           &
                 wf%n_ao)
!
     call mem%dealloc(L_inv_P_trans, wf%n_ao, wf%n_ao)
!
      call daxpy(wf%n_ao**2, -half, wf%ao_density, 1, D_v, 1)
!
  end subroutine construct_virtual_density_mlhf
!
end module mlhf_class
