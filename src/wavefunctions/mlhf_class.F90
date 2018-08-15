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
   use active_atoms_class
!
   use eri_cd_solver_class
!
   use array_utilities
!
   implicit none
!
   type, extends(hf):: mlhf
!
      type(active_atoms) :: active_space
!
   contains
!
      procedure :: initialize => initialize_mlhf
      procedure :: finalize => finalize_mlhf
!
      procedure :: read_info => read_info_mlhf
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
      call wf%read_info
!
      write(output%unit, *)'Doing SAD/ml test'
      write(output%unit, *)'Active atoms', wf%active_space%atoms
      flush(output%unit)
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
      real(dp) :: max
!
      integer(i15), dimension(:,:), allocatable :: active_aos
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
      write(output%unit,*)'SAD'
      flush(output%unit)
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
      write(output%unit,*)'Fock'
      flush(output%unit)
!
      call mem%alloc(sp_eri_schwarz, n_s, n_s)
      call wf%construct_sp_eri_schwarz(sp_eri_schwarz, n_s)
!  
      call mem%alloc(eri_deg, n_s**2, n_s**2)
      call wf%determine_degeneracy(eri_deg, n_s)
!
!     Construct initial AO Fock from the SOAD density
!
      call wf%initialize_ao_fock()
      call wf%construct_ao_fock(sp_eri_schwarz, eri_deg, n_s)
!
      call mem%dealloc(eri_deg, n_s**2, n_s**2)
      call mem%dealloc(sp_eri_schwarz, n_s, n_s)
!
!     Construct AO overlap matrix, Cholesky decompose it,
!     followed by preconditioning (making it the identity matrix
!     for this particular preconditioner - V)
!
      call wf%initialize_ao_overlap()
      call wf%construct_ao_overlap()
!
      write(output%unit,*)'Fock diag'
      flush(output%unit)
!
!     Solve Roothan Hall once - using the SOAD guess - to get a decent AO density
!     on  which to start the preconditioned conjugate gradient (PCG) algorithm
!
      call wf%initialize_orbital_energies()
      call wf%initialize_mo_coefficients()
      call wf%solve_roothan_hall() ! F^AO C = S C e to get new MOs C
!
!     Update the AO density
!
      write(output%unit,*)'Decompose densities'
      flush(output%unit)
!
      call wf%construct_ao_density() ! Construct AO density from C
!
      n_active_aos = 0
!
      do i = 1, wf%active_space%n_active_atoms
!
         n_active_aos = n_active_aos + wf%system%atoms(wf%active_space%atoms(i, 1))%n_ao
!
      enddo 
!
      call mem%alloc_int(active_aos, n_active_aos, 1)
!
      ao_offset = 0
      active_ao_counter = 0
!
      do i = 1, wf%system%n_atoms
         do j = 1, wf%active_space%n_active_atoms
!
            if (i == wf%active_space%atoms(j, 1)) then
!
               do k = 1, wf%system%atoms(i)%n_ao
!
                  active_ao_counter = active_ao_counter + 1
!
                  active_aos(active_ao_counter, 1) = k + ao_offset
!   
               enddo
!
            endif
!
         enddo
!
         ao_offset = ao_offset + wf%system%atoms(i)%n_ao
!
      enddo
!
!     Determine number of active occupied
!
      n_active_occ = 0
!
      do i = 1, wf%active_space%n_active_atoms
!
        n_active_occ = n_active_occ + wf%system%atoms(wf%active_space%atoms(i, 1))%number
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
                                                     n_vectors_occ, 1.0d-9, n_active_aos, active_aos, n_active_occ)
!
      call mem%alloc(cholesky_vectors_virt, wf%n_ao, n_active_aos)
      cholesky_vectors_virt = zero
!
      call cholesky_decomposition_limited_diagonal(ao_density_v, cholesky_vectors_virt, wf%n_ao, &
                                                     n_vectors_virt, 1.0d-9, n_active_aos, active_aos, n_active_vir)
!
      write(output%unit, *)n_vectors_occ, n_active_occ, n_vectors_virt, n_active_vir
      flush(output%unit)
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
   subroutine read_info_mlhf(wf)
!
!!
!!    Read information
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!     mlhf
!!        n_active_atoms
!!        active_atoms
!!     end mlhf
!!
      implicit none
!
      class(mlhf) :: wf
!
      character(len=100) :: line
      character(len=100) :: current_basis
!
      integer(i15) :: i = 0,j = 0, ioerror, first, last
!
      call input%init('eT.inp', 'sequential', 'formatted')
      call disk%open_file(input, 'read')
      rewind(input%unit)
!
      read(input%unit,'(a)', iostat=ioerror) line
      line = remove_preceding_blanks(line)
!
      do while ((trim(line) .ne. 'end mlhf') .and. (line(1:2) .ne. 'do'))
!
         if (trim(line) == 'mlhf') then ! found cholesky section in input
!
            read(input%unit, *) wf%active_space%n_active_atoms
!
            call mem%alloc_int(wf%active_space%atoms, wf%active_space%n_active_atoms, 1)
!
            read(input%unit,'(a)', iostat=ioerror) line
            line = remove_preceding_blanks(line)
!
            if (line(1:1)=='[') then ! range given
!
               do i = 2, 100
!
                  if (line(i:i) == ',') exit
!
               enddo
!
               read(line(2:i-1), *) first
!
               do j = i, 100
!
                  if (line(j:j) == ']') exit
!
               enddo
!
               read(line(i+1:j-1), *) last
!
               if ((last - first + 1) .ne. wf%active_space%n_active_atoms) then
!
                  write(output%unit, *) 'Error: actual number of active atoms does not match number provided in input.'
                  stop
!
               endif
!
               do i = first, last
!
                  wf%active_space%atoms(i - first + 1, 1) = i
!
               enddo
!
            else ! specific atoms are given
!
               read(line, *) wf%active_space%atoms
!
            endif
!
         endif
!
         read(input%unit,'(a)') line
         line = remove_preceding_blanks(line)
!
      enddo
      backspace(input%unit)
!
      call disk%close_file(input)
!
   end subroutine read_info_mlhf
!
!
  subroutine construct_virtual_density_mlhf(wf, D_v)
!!
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
!  1. Hvordan skal jeg få rutiner fra HF (midlertidig arve derfra?)
!  2. Skrive inputleser
!  3. Teste dekomponering
!  4. Lage V_xy
!  5. ML
!
!
end module mlhf_class
