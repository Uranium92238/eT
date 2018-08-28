module hf_solver_class
!!
!!    Abstract HF solver class module 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
   use kinds 
   use parameters
   use hf_class
!
   implicit none 
!
   type, abstract :: hf_solver 
!
!     Thresholds and other general settings
!
      real(dp) :: energy_threshold   = 1.0D-6
      real(dp) :: residual_threshold = 1.0D-6
!
      integer(i15) :: max_iterations = 100
!
!     Information for the transformation from linearly dependent to independent AO set
!
      real(dp), dimension(:,:), allocatable :: cholesky_ao_overlap
      real(dp), dimension(:,:), allocatable :: permutation_matrix  
!
      real(dp) :: linear_dependence_threshold = 1.0D-6
!
      real(dp) :: coulomb_thr       = 1.0D-11 ! screening 
      real(dp) :: coulomb_precision = 1.0D-14 ! integral accuracy
      real(dp) :: exchange_thr      = 1.0D-11  ! screening 
!
   contains 
!
!     Main routines (called in order initialize, run, finalize)
!
      procedure(essential), deferred :: initialize 
      procedure(essential), deferred :: run 
      procedure(essential), deferred :: finalize
!
      procedure :: do_roothan_hall      => do_roothan_hall_hf_solver 
      procedure :: decompose_ao_overlap => decompose_ao_overlap_hf_solver
!
      procedure :: read_settings           => read_settings_hf_solver
      procedure :: read_hf_solver_settings => read_hf_solver_settings_hf_solver
!
   end type hf_solver
!
   abstract interface
!
      subroutine essential(solver, wf)
!
         import :: hf, hf_solver
!
         implicit none 
!
         class(hf_solver) :: solver 
!
         class(hf) :: wf 
!
      end subroutine essential
!
   end interface
!
contains 
!
   subroutine do_roothan_hall_hf_solver(solver, wf)
!!
!!    Do Roothan-Hall
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Solves the equation F C = S C e for the orbital coefficients C.
!!    More precisely, it solves the equation in a linearly independent
!!    subspace,
!!
!!       P^T F P (P^T C) = P^T S P P^T C e = L L^T (P^T C) e,
!!
!!    which differs from the AO basis if linear dependence is present 
!!    to within some given threshold. The components of the equation is 
!!    given by the Cholesky decomposition
!!
!!       P^T S P = L L^T,
!!
!!    where P is referred to as the 'permutation matrix' and L the 'cholesky
!!    ao overlap' (note that these are member variables of the solver which 
!!    must be set by a call to solver%decompose_ao_overlap(wf)). The number 
!!    of linearly independent orbitals is wf%n_mo, whereas the full number 
!!    is wf%n_ao.   
!!
      implicit none
!
      class(hf_solver) :: solver 
!
      class(hf) :: wf
!
      real(dp), dimension(:,:), allocatable :: work
      real(dp), dimension(:,:), allocatable :: metric 
      real(dp), dimension(:,:), allocatable :: ao_fock 
      real(dp), dimension(:,:), allocatable :: orbital_energies 
      real(dp), dimension(:,:), allocatable :: FP 
!
      real(dp) :: ddot, norm
!
      integer(i15) :: info = 0
!
      call mem%alloc(metric, wf%n_mo, wf%n_mo)
!
      call dgemm('N','T',                     &
                  wf%n_mo,                    &
                  wf%n_mo,                    &
                  wf%n_mo,                    &
                  one,                        &
                  solver%cholesky_ao_overlap, & 
                  wf%n_mo,                    &
                  solver%cholesky_ao_overlap, &
                  wf%n_mo,                    &
                  zero,                       &
                  metric,                     & ! metric = L L^T
                  wf%n_mo)
!
!     Allocate reduced space matrices 
!
      call mem%alloc(ao_fock, wf%n_mo, wf%n_mo)
      call mem%alloc(orbital_energies, wf%n_mo, 1)
!
!     Construct reduced space Fock matrix, F' = P^T F P,
!     which is to be diagonalized over the metric L L^T 
!
      call mem%alloc(FP, wf%n_ao, wf%n_mo)
!
      call dgemm('N','N',                    &
                  wf%n_ao,                   &
                  wf%n_mo,                   &
                  wf%n_ao,                   &
                  one,                       &
                  wf%ao_fock,                &
                  wf%n_ao,                   &
                  solver%permutation_matrix, &
                  wf%n_ao,                   &
                  zero,                      &  
                  FP,                        &
                  wf%n_ao)
!
      call dgemm('T','N',                    &
                  wf%n_mo,                   &
                  wf%n_mo,                   &
                  wf%n_ao,                   &
                  one,                       &
                  solver%permutation_matrix, &
                  wf%n_ao,                   &
                  FP,                        &
                  wf%n_ao,                   &
                  zero,                      &
                  ao_fock,                   & ! F' = P^T F P
                  wf%n_mo)   
!
      call mem%dealloc(FP, wf%n_mo, wf%n_ao)   
!
!     Solve F'C' = L L^T C' e
!
      info = 0
!
      call mem%alloc(work, 4*wf%n_mo, 1)
      work = zero
!
      call dsygv(1, 'V', 'L',       &
                  wf%n_mo,          &
                  ao_fock,          & ! ao_fock on entry, orbital coefficients on exit
                  wf%n_mo,          &
                  metric,           &
                  wf%n_mo,          &
                  orbital_energies, &
                  work,             &
                  4*(wf%n_mo),      &
                  info)
!
      call mem%dealloc(metric, wf%n_mo, wf%n_mo)
      call mem%dealloc(work, 4*wf%n_mo, 1)
      call mem%dealloc(orbital_energies, wf%n_mo, 1)
!
      if (info .ne. 0) then 
!
         write(output%unit, '(/t3,a/)') 'Error: could not solve Roothan-Hall equations.'
         stop
!
      endif
!
!     Transform back the solutions to original basis, C = P (P^T C) = P C'
!
      wf%orbital_coefficients = zero
!
      call dgemm('N','N',                    &
                  wf%n_ao,                   &
                  wf%n_mo,                   &
                  wf%n_mo,                   &
                  one,                       &
                  solver%permutation_matrix, &
                  wf%n_ao,                   &
                  ao_fock,                   & ! orbital coefficients 
                  wf%n_mo,                   &
                  zero,                      &
                  wf%orbital_coefficients,   &
                  wf%n_ao)
!
      call mem%dealloc(ao_fock, wf%n_mo, wf%n_mo)
!
   end subroutine do_roothan_hall_hf_solver
!
!
   subroutine decompose_ao_overlap_hf_solver(solver, wf)
!!
!!    Decompose AO overlap
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Performs a Cholesky decomposition of the AO overlap matrix S,
!!    to within a given threshold: 
!!
!!       P^T S P = L L^T 
!!
!!    The routine allocates and sets P and L, the 'permutation matrix'
!!    and the 'cholesky ao overlap'. Moreover, it sets the number of 
!!    linearly independent AOs, wf%n_mo, which is sometimes less than 
!!    wf%n_ao. From P and L, we can transform equations to the linearly 
!!    independent basis and back.
!!
      implicit none
!
      class(hf_solver) :: solver
!
      class(hf) :: wf
!
      integer(kind=4), dimension(:, :), allocatable :: used_diag
!
      real(dp), dimension(:, :), allocatable :: L
!
      integer(i15) :: i, j
!
      allocate(used_diag(wf%n_ao, 1))
      used_diag = 0
!
      call mem%alloc(L, wf%n_ao, wf%n_ao) ! Full Cholesky vector
      L = zero
!
      call full_cholesky_decomposition_system(wf%ao_overlap, L, wf%n_ao, wf%n_mo, &
                                          solver%linear_dependence_threshold, used_diag)
!
      call mem%alloc(solver%cholesky_ao_overlap, wf%n_mo, wf%n_mo) 
      solver%cholesky_ao_overlap(:,:) = L(1:wf%n_mo, 1:wf%n_mo)
!
      call mem%dealloc(L, wf%n_ao, wf%n_ao)
!
!     Make permutation matrix P
!
      call mem%alloc(solver%permutation_matrix, wf%n_ao, wf%n_mo)
!
      solver%permutation_matrix = zero
!
      do j = 1, wf%n_mo
!
         solver%permutation_matrix(used_diag(j, 1), j) = one
!
      enddo
!
      deallocate(used_diag)
!
   end subroutine decompose_ao_overlap_hf_solver
!
!
   subroutine read_settings_hf_solver(solver)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
!!    Reads the settings. This routine is to be overwritten by 
!!    descendants if more settings need to be set. 
!!
      implicit none 
!
      class(hf_solver) :: solver 
!
      call solver%read_hf_solver_settings()
!
   end subroutine read_settings_hf_solver
!
!
   subroutine read_hf_solver_settings_hf_solver(solver)
!!
!!    Read HF solver settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
!!    Reads the settings specific to this class.
!!
      implicit none 
!
      class(hf_solver) :: solver 
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
            if (line(1:17) == 'energy_threshold:') then
!
               value = line(18:100)
               value = remove_preceding_blanks(value)
               read(value, *) solver%energy_threshold
               return
!
            elseif (line(1:19) == 'residual_threshold:') then 
!
               value = line(20:100)
               value = remove_preceding_blanks(value)
               read(value, *) solver%residual_threshold
               return
!
            elseif (line(1:15) == 'max_iterations:') then 
!
               value = line(16:100)
               value = remove_preceding_blanks(value)
               read(value, *) solver%max_iterations
               return
!
            elseif (line(1:18) == 'coulomb_threshold:') then 
!
               value = line(19:100)
               value = remove_preceding_blanks(value)
               read(value, *) solver%coulomb_thr
               return
!
            elseif (line(1:19) == 'exchange_threshold:') then 
!
               value = line(20:100)
               value = remove_preceding_blanks(value)
               read(value, *) solver%exchange_thr
               return
!
            elseif (line(1:18) == 'coulomb_precision:') then 
!
               value = line(19:100)
               value = remove_preceding_blanks(value)
               read(value, *) solver%coulomb_precision
               return
!
            elseif (line(1:21) == 'linear_dep_threshold:') then 
!
               value = line(22:100)
               value = remove_preceding_blanks(value)
               read(value, *) solver%linear_dependence_threshold
               return
!
            endif
!
         enddo
!
      endif 
!
   end subroutine read_hf_solver_settings_hf_solver
!
!
end module hf_solver_class
