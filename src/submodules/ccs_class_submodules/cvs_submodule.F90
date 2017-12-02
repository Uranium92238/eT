submodule(ccs_class) cvs
!
!!
!!    CVS submodule(CCS)
!!    Written by Eirik F. Kjønstad and Sarai Dery Folkestad, May 2017
!!
!!    Contains the CVS routines for the CCS class and routines for core excitations.
!!    Note that this submodule contains both excited state routines and jacobian transformation routines.
!!
!!    -::- Subroutines in this submodule -::-
!!
!!    initialize_trial_vectors_core    - Driver for initialization of trial vectors
!!    find_start_trial_indices_core    - Finds indices of single core excitations of lowest energy to be used for start vectors
!!    precondition_residual_core       - Projects out contamination (valence contributions) and preconditions uncontaminated residual
!!    cvs_residual_projection          - Projection routine for residual
!!    core_jacobian_ccs_transformation - Jacobian transformation where contamination (valence contributions) is projected out
!!    cvs_rho_a_i_projection           - Projection routine for rho_a_i
!!    
!!    Helper subroutines:
!!
!!    find_core_mo - identifies the mo-index of the core orbital of interest.
!!
!!
!
!
contains
!
!
   module subroutine initialize_trial_vectors_core_ccs(wf)
!!
!!    Initialize trial vectors, for core excitation calculation 
!!    Written by Sarai D. Folkestad, Aug. 2017
!!
!!    Finds correct core MO, and selects the n_singlet_state lowest 
!!    orbital differences where one of the occupied indices corresponds to the 
!!    core MO
!!
      implicit none
!
      class(ccs) :: wf
!
      integer(i15), dimension(:,:), allocatable :: index_core_obital
!
      real(dp), dimension(:,:), allocatable :: c
! 
      integer(i15) :: i = 0, j = 0
!
      integer(i15) :: unit_trial_vecs = 0, unit_rho = 0, ioerror = 0
!
!     Allocate array for the indices of the lowest orbital differences
!
      call allocator_int( index_core_obital, wf%excited_state_specifications%n_singlet_states, 1)
      index_core_obital = zero
!
!     Find indecies of lowest orbital differences
!
      call wf%find_start_trial_indices_core(index_core_obital)
!
!     Generate start trial vectors c and write to file
!
      call allocator(c, wf%n_parameters, 1)
!
!     Prepare for writing trial vectors to file
!
      call generate_unit_identifier(unit_trial_vecs)
      open(unit=unit_trial_vecs, file='trial_vec', action='write', status='unknown', &
        access='direct', form='unformatted', recl=dp*(wf%n_parameters), iostat=ioerror)
!
      do i = 1, (wf%excited_state_specifications%n_singlet_states)
         c = zero
         c(index_core_obital(i,1),1) = one
         write(unit_trial_vecs, rec=i, iostat=ioerror) (c(j,1), j = 1, wf%n_parameters)
      enddo
!
!     Close file
!  
      close(unit_trial_vecs)
!
!     Deallocate c
!
      call deallocator(c, wf%n_parameters, 1)
!
!     Deallocate index_lowest_obital_diff
!
      call deallocator_int(index_core_obital, wf%excited_state_specifications%n_singlet_states, 1)
!
   end subroutine initialize_trial_vectors_core_ccs
!
!
   module subroutine find_start_trial_indices_core_ccs(wf, index_list)
!!
!!    Find indices for lowest orbital differences
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad
!!
!!
      implicit none
!
      class(ccs) :: wf
      integer(i15), dimension(wf%excited_state_specifications%n_singlet_states,1), intent(inout) :: index_list
!
      real(dp), dimension(:,:), allocatable     ::  sorted_short_vec
!
      integer(i15) :: a, i, counter
!
!     Find core mo(s)
!
      call allocator_int(wf%core_excited_state_specifications%index_core_mo, &
                           wf%core_excited_state_specifications%n_equivalent_cores, 1)
!
      call wf%find_core_mo
!
      counter = 1
      do a = 1, wf%n_v
         if ( counter .le.  wf%excited_state_specifications%n_singlet_states)  then
            do i = 1, wf%core_excited_state_specifications%n_equivalent_cores
!
               index_list(counter, 1) = index_two(a, wf%core_excited_state_specifications%index_core_mo(i, 1), wf%n_v)
!
               counter = counter + 1
!
               if ( counter .gt.  wf%excited_state_specifications%n_singlet_states) exit
!
            enddo
         endif
      enddo
!
   end subroutine find_start_trial_indices_core_ccs
!
   module subroutine find_core_mo_ccs(wf)
!!
!!    Find which mo are core mos
!!    Written by Sarai D. Folkestad, Aug. 2017
!!
      implicit none
!
      class(ccs) :: wf
!
      integer(i15) :: n_nuclei = 0, core = 0, ao = 0, mo = 0, ao_mo_index = 0
      integer(i15) :: first_ao_on_core = 0, counter = 0, n_aos_on_atoms = 0, i = 0, j = 0
      integer(i15), dimension(:,:), allocatable :: n_ao_on_center, ao_center_info, aos_on_atoms
!
!     :: Get center info ::
!
      call read_atom_info(n_nuclei, wf%n_ao)
!
!     n_ao_on_center contains number of aos on each atom
!     ao_center_info contains ao index (first column) belonging to each of the aos (second column)
!     
      call allocator_int(n_ao_on_center, n_nuclei, 1)      
      call allocator_int(ao_center_info, wf%n_ao, 2)
!
      call read_center_info(n_nuclei, wf%n_ao, n_ao_on_center, ao_center_info)
!
!     Find number of aos on atoms that we are interested in
!   
      n_aos_on_atoms = 0
      do i = 1, wf%core_excited_state_specifications%n_equivalent_cores
         n_aos_on_atoms = n_aos_on_atoms + n_ao_on_center(wf%core_excited_state_specifications%cores(i,1),1)
      enddo
!
      call allocator_int(aos_on_atoms, n_aos_on_atoms, 1)
!
      write(unit_output,*)wf%core_excited_state_specifications%n_equivalent_cores, wf%core_excited_state_specifications%cores(1,1)
      flush(unit_output)
!
      counter = 1    
      do i = 1, wf%core_excited_state_specifications%n_equivalent_cores 
         do j = 1, wf%n_ao
            if (ao_center_info(j,1) == wf%core_excited_state_specifications%cores(i,1)) then
!
               aos_on_atoms(counter,1) = ao_center_info(j,2)
               counter = counter + 1  
!         
             endif
         enddo
      enddo
!
      call deallocator_int(ao_center_info, wf%n_ao, 2)
!
!     :: Find core mo that has large ao component on the atom in question ::
!
      wf%core_excited_state_specifications%index_core_mo = zero
      counter = 0
!
      do ao = 1, n_aos_on_atoms
!
         do  mo = 1, wf%n_o
!
!           Determine wether orbital in core mo
!
            if(wf%fock_diagonal(mo, 1) .lt. -5.0d0) then

!
!              Determine wether the core orbital sits on the corect atom(s)
!
               ao_mo_index = index_two(aos_on_atoms(ao, 1), mo, wf%n_ao)
!
               if(abs(wf%mo_coef(ao_mo_index, 1)) .gt. 0.5d0) then
                  counter = counter + 1
!
                  if (counter .le. wf%core_excited_state_specifications%n_equivalent_cores) then
                     wf%core_excited_state_specifications%index_core_mo(counter, 1) = mo
                  endif
!
                     
               endif
!
            endif
!
         enddo
!
      enddo
!
!
      call deallocator_int(n_ao_on_center, n_nuclei, 1) 
      call deallocator_int(aos_on_atoms, n_aos_on_atoms, 1) 
!
!     Sanity check
!
      do core = 1, wf%core_excited_state_specifications%n_equivalent_cores
         if (wf%core_excited_state_specifications%index_core_mo(core, 1) .eq. 0) then
            write(unit_output,*)'WARNING: Found no core orbitals for core', wf%core_excited_state_specifications%cores(core, 1)
            stop
         endif
      enddo
!
   end subroutine find_core_mo_ccs
!
!
   module subroutine precondition_residual_core_ccs(wf, residual)
!!
!!    Precondition residual for core excited state calculation
!!    Written by Sarai D. Folkestad, Aug. 2017
!!
!!    Project out elements not corresponding to the core excitation
!!    Divide elements of residual by orbital difference
!!
      implicit none
!
      class(ccs) :: wf
      real(dp), dimension(wf%n_parameters ,1) :: residual
! 
      real(dp), dimension(:,:), allocatable :: orbital_diff
!
      integer(i15) :: i = 0, a = 0, ai = 0
!   
      call wf%cvs_residual_projection(residual)
!
      call allocator(orbital_diff, wf%n_parameters, 1)
      orbital_diff = zero
!
      call wf%calculate_orbital_differences(orbital_diff)
!
      do i = 1, wf%n_parameters
!
         residual(i, 1) = residual(i,1)/orbital_diff(i,1)
!
      enddo
!
      call deallocator(orbital_diff, wf%n_parameters, 1)
!
   end subroutine precondition_residual_core_ccs
!
!
   module subroutine cvs_residual_projection_ccs(wf, residual)
!!
!!    Residual projection for core excited state calculation (CCS)
!!    Written by Sarai D. Folkestad, Aug. 2017
!!
      implicit none
!
      class(ccs) :: wf
      real(dp), dimension(wf%n_parameters, 1) :: residual
!
      integer(i15) :: i = 0, a = 0, core = 0, ai = 0
!
      logical :: core_orbital
!
      do i = 1, wf%n_o
!
         core_orbital = .false.
         do core = 1, wf%core_excited_state_specifications%n_equivalent_cores
!
            if (i .eq. wf%core_excited_state_specifications%index_core_mo(core, 1)) core_orbital = .true.
!
         enddo
!
         if (.not. core_orbital) then
            do a = 1, wf%n_v
               ai = index_two(a, i, wf%n_v)
               residual(ai, 1) = zero
            enddo
         endif
!
      enddo
!
   end subroutine cvs_residual_projection_ccs
!
!
   module subroutine cvs_rho_a_i_projection_ccs(wf, vec_a_i)
!!
!!    CVS projection of rho_a_i, 
!!    Written by Sarai D. Folkestad, Aug. 2017
!!
!!    Projects out elements of rho that do not correspond to the core excitation.
!!
      implicit none
!
      class(ccs) :: wf
      real(dp), dimension(wf%n_v ,wf%n_o) :: vec_a_i
!
      integer(i15) :: i = 0, a = 0, core = 0
!
      logical :: core_orbital
!
      do i = 1, wf%n_o
!
         core_orbital = .false.
!
         do core = 1, wf%core_excited_state_specifications%n_equivalent_cores
!
            if (i .eq. wf%core_excited_state_specifications%index_core_mo(core, 1)) core_orbital = .true.
!
         enddo
!
         if (.not. core_orbital) then
            vec_a_i(:,i) = zero
         endif
!
      enddo
!
   end subroutine cvs_rho_a_i_projection_ccs
!
!
end submodule cvs
