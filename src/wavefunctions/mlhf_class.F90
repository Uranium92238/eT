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
      procedure :: prepare => prepare_mlhf
      procedure :: finalize => finalize_mlhf
!
      procedure :: construct_virtual_density => construct_virtual_density_mlhf
      procedure :: construct_virtual_density_from_MO => construct_virtual_density_from_MO_mlhf
!
      procedure :: eri_decomp_test_w_active_dens => eri_decomp_test_w_active_dens_mlhf
!
   end type mlhf
!
!
contains
!
!
   subroutine prepare_mlhf(wf)
!!
!!    Initialize
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(mlhf) :: wf
!
      integer(i15) :: ao
!
      real(dp), dimension(:,:), allocatable :: density_diagonal
!
      wf%name = 'MLHF'
!
      write(output%unit, '(a/)')':: SAD-ML test'
      flush(output%unit)
!
      call wf%system%prepare()
!
      call wf%system%print_geometry()
!
!     Construct AO overlap matrix, Cholesky decompose it,
!     followed by preconditioning (making it the identity matrix
!     for this particular preconditioner - V)
!
      wf%n_ao = wf%system%get_n_aos()
!
      call initialize_coulomb()
      call initialize_kinetic()
      call initialize_nuclear()
      call initialize_overlap()
!
      call wf%initialize_ao_overlap()
      call wf%construct_ao_overlap()
      call wf%decompose_ao_overlap() 
!
      wf%n_o = (wf%system%get_n_electrons())/2
      wf%n_v = wf%n_mo - wf%n_o
!
      write(output%unit, '(/t6, a8, i3, a23)')'Removed ', wf%n_ao - wf%n_mo, ' AOs due to linear dep.'
!
      call wf%initialize_ao_density()
!
!     Set initial density to superposition of atomic densities (SOAD) guess
!
      wf%ao_density = zero
!
      call wf%set_ao_density_to_sad()
!
      call wf%eri_decomp_test_w_active_dens()
!
   end subroutine prepare_mlhf
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
      integer(i15):: a
!
      real(dp) :: max_val, e_construct_fock, s_construct_fock, omp_get_wtime, x, y, z
!
      integer(i15), dimension(:,:), allocatable :: active_aos, sp_eri_schwarz_list
!
      integer(i15) :: n_active_occ, n_active_vir, n_s
!
      real(dp), dimension(:,:), allocatable :: eri_deg, sp_eri_schwarz, h_wx
!
      type(eri_cd_solver)  :: chol_solver
!
      call wf%initialize_ao_fock()
!
      s_construct_fock = omp_get_wtime()
!
      n_s = wf%system%get_n_shells()
!
      call wf%construct_ao_fock_SAD(1.0d-13, 1.0d-13, 1.0d-25)
!
    ! write(output%unit,*)wf%energy
!
    ! call mem%alloc(sp_eri_schwarz, n_s*(n_s + 1)/2, 2)
    ! call mem%alloc_int(sp_eri_schwarz_list, n_s*(n_s + 1)/2, 3)
!
    ! call wf%construct_sp_eri_schwarz(sp_eri_schwarz, sp_eri_schwarz_list, n_s)
!
!   ! Set the initial density guess and Fock matrix 
!
    ! call mem%alloc(h_wx, wf%n_ao, wf%n_ao)
    ! call wf%get_ao_h_wx(h_wx)
!
    ! call wf%update_fock_and_energy(sp_eri_schwarz, sp_eri_schwarz_list, n_s, h_wx)
    ! write(output%unit,*)wf%energy
    ! flush(output%unit)
    ! stop
!
      e_construct_fock = omp_get_wtime()
      write(output%unit, '(/t6, a49, f11.2)')'Wall time to construct AO fock from SAD density: ', &
                                  e_construct_fock - s_construct_fock
      flush(output%unit)
!
!     Solve Roothan Hall once - using the SOAD guess - to get a decent AO density
!     on  which to start the preconditioned conjugate gradient (PCG) algorithm
!
      call wf%initialize_orbital_energies()
      call wf%initialize_orbital_coefficients()
      call wf%do_roothan_hall(wf%ao_fock, wf%orbital_coefficients, wf%orbital_energies) ! F^AO C = S C e to get new MOs C
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
!     Add orbitals if bonds are capped. Assumes 2.5Å for covalent bonds.
!
      do i = wf%system%n_active_atoms + 1, wf%system%n_atoms
!
        do j = 1, wf%system%n_active_atoms
!
          x = (wf%system%atoms(j)%x - wf%system%atoms(i)%x)
          y = (wf%system%atoms(j)%y - wf%system%atoms(i)%y)
          z = (wf%system%atoms(j)%z - wf%system%atoms(i)%z)
!
          if (sqrt(x**2 + y**2 + z**2) .le. 1.5d0) n_active_occ = n_active_occ + 1 
!
        enddo
!
      enddo
!
!     Determine number of active virtuals
!
      n_active_vir = n_active_occ * (wf%n_v/wf%n_o)
!
      call mem%alloc(ao_density_v, wf%n_ao, wf%n_ao)
      call wf%construct_virtual_density_from_MO(ao_density_v)
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
      do j = 1, wf%n_ao
!
         max_val = 0.0d0
!
         do i = 1, n_vectors_occ
!
            if (cholesky_vectors_occ(j, i)**2 .gt. max_val) max_val = cholesky_vectors_occ(j, i)**2     
!
         enddo
!
         do a = 1, n_vectors_virt
!
            if (cholesky_vectors_virt(j, a)**2 .gt. max_val) max_val = cholesky_vectors_virt(j, a)**2 
!
         enddo
!
         V(j, 1) = max_val
!
      enddo
!
      call mem%dealloc(cholesky_vectors_virt, wf%n_ao, n_active_aos)
      call mem%dealloc(cholesky_vectors_occ, wf%n_ao, n_active_aos)
!
      write(output%unit, '(/t6, a27, i4)')'Number of occupied active: ', n_vectors_occ
      write(output%unit, '(t6, a27, i4/)')'Number of virtual active:  ', n_vectors_virt
      flush(output%unit)
!
!
      call chol_solver%prepare(wf%system)
      call chol_solver%run(wf%system, V)
      call chol_solver%cleanup()
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
                                               1.0d-6, piv)
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
!
  subroutine construct_virtual_density_from_MO_mlhf(wf, D_v)
!!
!!    Construct virtual density
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    D^V = sum_a C_alpha,a * C_beta,a 
!!
      implicit none 
!
      class(mlhf) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao) :: D_v

!
      call dgemm('N', 'T',                              &
                wf%n_ao,                                &
                wf%n_ao,                                &
                wf%n_mo - wf%n_o,                       &
                one,                                    &
                wf%orbital_coefficients(1, wf%n_o + 1), &
                wf%n_ao,                                &
                wf%orbital_coefficients(1, wf%n_o + 1), &
                wf%n_ao,                                &
                zero,                                   &  
                D_v,                                    &
                wf%n_ao)
!
  end subroutine construct_virtual_density_from_MO_mlhf
!
!
end module mlhf_class
