submodule(ccs_class) ionized_state
!
!!
!!    Ionized state submodule(CCS)
!!    Written by Eirik F. Kjønstad and Sarai Dery Folkestad, October 2017
!!
!!    Contains the CCS class routines for ionization, both core and valence.
!!    Note that this submodule contains both excited state routines and jacobian transformation routines.
!!
!!    -::- Subroutines in this submodule -::-
!!
!!    ionized_state_driver                         - Driver for ionized state calculation
!!    initialize_trial_vectors_valence_ionization  - Initializes trial vectors for valence ionization calculation 
!!    initialize_trial_vectors_core_ionization     - Initializes trial vectors for core ionization calculation
!!    precondition_residual_valence_ionization     - Projects out contamination (regular excitations) 
!!                                                   and preconditions uncontaminated residual 
!!    ionization_residual_projection               - Projects out contamination (regular excitations) from residual
!!    ionization_rho_a_i_projection                - Projection routine for rho_a_i
!!                                                   (valence contributions and regular excitations) are projected out
!!    precondition_residual_core_ionization        - Projects out contaminations (valence contributiona and regular excitations) 
!!                                                   and preconditions uncontaminated residual 
!!    
!
   implicit none 
!
!
contains
!
!
   module subroutine ionized_state_driver_ccs(wf)
!!
!!    Ionizedstate driver (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
!!    Directs the solution of the ionizedstate problem for CCS. The
!!    routine is inherited is to be inherited unaltered in the CC hierarchy. 
!!
!!    Note: it is only necessary to alter this routine if the ionized states are 
!!    solved for by a different algorithm (such as in similarity constrained CC, 
!!    where the excited states and ground state are determined simultaneously).
!!
      implicit none 
!
      class(ccs) :: wf 
!
!     Let the user know the excited state driver is running
!
      write(unit_output,'(t3,a)')    ':: Ionized state solver (Davidson)'
      write(unit_output,'(t3,a/)')   ':: E. F. Kjønstad, S. D. Folkestad, May 2017'
      write(unit_output,'(t3,a,i2,a,a,a)') &
                                     'Requested ',wf%excited_state_specifications%n_singlet_states,&
                                       ' ', trim(wf%name), ' singlet states.'
!
!     Run the general solver routine (file names are given
!     by the task, i.e., the file 'right_valence' contains
!     the right eigenvectors)
!
      call wf%excited_state_solver
!
   end subroutine ionized_state_driver_ccs
!
!
   module subroutine initialize_trial_vectors_valence_ionization_ccs(wf)
!!
!!    Initialize trial vectors for valence ionized state
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad
!!
!!    Initializes start trial vectors for the calculation of 
!!    singlet excited states and writes them to file 'trial_vecs'.
!!
!!    n start vectors are constructed by finding the n lowest orbital differences,      
!!    where n = n_singlet_states. Vector i has a 1.0D0 at the element corresponding to the i'th lowest
!!    orbital difference and 0.0d0 everywhere else
!!
      implicit none
!
      class(ccs) :: wf
!
      integer(i15) :: mo = 0, i = 0, j = 0, diffuse_mo = 0
!       
      integer(i15), dimension(:,:), allocatable :: index_lowest_obital_diff 
!
      real(dp), dimension(:,:), allocatable :: c
!
      integer(i15) :: unit_trial_vecs = 0, ioerror = 0
!
!		Find MO corresponding to super diffuse orbital
!     
!     Super diffuse AO is the last AO (for dalton integrals and HF this corresponds to 
!     the ghost atom being declared at the end of the MOLECULE.INP file). Must locate the MO
!     with weight 1.0D0 on this AO.
!
      do mo = 1, wf%n_mo
!
         i = index_two(wf%n_ao, mo, wf%n_ao)
         if (abs(wf%mo_coef(i,1) - 1.0D0) .lt. 1.0D-7) diffuse_mo = mo
!
      enddo
!
!     Sanity check - Did we find diffuse mo?
!
      if (diffuse_mo .eq. 0) then
         write(unit_output,*)'Error: Super diffuse mo not found.'
         stop
      endif
!
!     Allocate array for the indices of the lowest orbital differences
!
      call wf%mem%alloc_int( index_lowest_obital_diff, wf%excited_state_specifications%n_singlet_states, 1)
      index_lowest_obital_diff = zero
!
!     Select start vectors corresponding to excitation from HOMO, HOMO - 1, ..., (HOMO - n_singlet_states + 1)
!     to super diffuse MO
!
      do i = 1, wf%excited_state_specifications%n_singlet_states
         j = index_two(diffuse_mo - wf%n_o, wf%n_o + 1 - i, wf%n_v)
         index_lowest_obital_diff(i, 1) = j
      enddo

      do i = 1, wf%excited_state_specifications%n_singlet_states
         write(unit_output,*)index_lowest_obital_diff(i, 1)
         write(unit_output,*)index_two(1,8, wf%n_v),index_two(1,7, wf%n_v)  
      enddo
!
!     Generate start trial vectors c and write to file
!
      call wf%mem%alloc(c, wf%n_parameters, 1)
!
!     Prepare for writing trial vectors to file
!
      call generate_unit_identifier(unit_trial_vecs)
      open(unit=unit_trial_vecs, file='trial_vec', action='write', status='unknown', &
           access='direct', form='unformatted', recl=dp*(wf%n_parameters), iostat=ioerror)
!
      do i = 1, (wf%excited_state_specifications%n_singlet_states)
         c = zero
         c(index_lowest_obital_diff(i,1),1) = one
         write(unit_trial_vecs, rec=i, iostat=ioerror) (c(j,1), j = 1, wf%n_parameters)
      enddo
!
!     Close file
!     
      close(unit_trial_vecs)
!
!     Deallocate c
!
      call wf%mem%dealloc(c, wf%n_parameters, 1)
!
!     Deallocate index_lowest_obital_diff
!
      call wf%mem%dealloc_int(index_lowest_obital_diff, wf%excited_state_specifications%n_singlet_states, 1)
!
   end subroutine initialize_trial_vectors_valence_ionization_ccs
!
!
   module subroutine initialize_trial_vectors_core_ionization_ccs(wf)
!!
!!    Initialize trial vectors for core ionized state
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad
!!
!!    Initializes start trial vectors for the calculation of 
!!    singlet excited states and writes them to file 'trial_vecs'.
!!
!!    n start vectors are constructed by finding the n lowest orbital differences,      
!!    where n = n_singlet_states. Vector i has a 1.0D0 at the element corresponding to the i'th lowest
!!    orbital difference and 0.0d0 everywhere else
!!
      implicit none
!
      class(ccs) :: wf
!
      integer(i15) :: mo = 0, i = 0, j = 0, k = 0, diffuse_mo = 0
!
      real(dp), dimension(:,:), allocatable :: c
!
      integer(i15) :: unit_trial_vecs = 0, ioerror = 0
!
!     Sanity check - does number of states correspond to number of equivalent cores ?
!
      if (wf%core_excited_state_specifications%n_equivalent_cores .ne. wf%excited_state_specifications%n_singlet_states) then
         write(unit_output,*)'Error: Using super diffuse orbital for XPS calc only allows for one root per equivalent core.'
         stop
      endif
!
!     Find MO corresponding to super diffuse orbital
!     
!     Super diffuse AO is the last AO (for dalton integrals and HF this corresponds to 
!     the ghost atom being declared at the end of the MOLECULE.INP file). Must locate the MO
!     with weight 1.0D0 on this AO.
!
      do mo = 1, wf%n_mo
!
         i = index_two(wf%n_ao, mo, wf%n_ao)
         if (abs(wf%mo_coef(i,1) - 1.0D0) .lt. 1.0D-7) diffuse_mo = mo
!
      enddo
!
!     Sanity check - Did we find diffuse mo?
!
      if (diffuse_mo .eq. 0) then
         write(unit_output,*)'Error: Super diffuse mo not found.'
         stop
      endif

!
!     Generate start trial vectors c and write to file
!
      call wf%mem%alloc(c, wf%n_parameters, 1)
!
!     Prepare for writing trial vectors to file
!
      call generate_unit_identifier(unit_trial_vecs)
      open(unit=unit_trial_vecs, file='trial_vec', action='write', status='unknown', &
           access='direct', form='unformatted', recl=dp*(wf%n_parameters), iostat=ioerror)
!
!     Find core mo(s)
!
      call wf%mem%alloc_int(wf%core_excited_state_specifications%index_core_mo,&
                         wf%core_excited_state_specifications%n_equivalent_cores, 1)
!
      call wf%find_core_mo

!
!     Select start vectors corresponding to excitation from CORE to super diffuse MO
!
      do i = 1, wf%core_excited_state_specifications%n_equivalent_cores
         c = zero
         k = index_two(diffuse_mo - wf%n_o, wf%core_excited_state_specifications%index_core_mo(i,1), wf%n_v)
         c(k,1) = one
         write(unit_trial_vecs, rec=i, iostat=ioerror) (c(j,1), j = 1, wf%n_parameters)
      enddo
!
!     Close file
!     
      close(unit_trial_vecs)
!
!     Deallocate c
!
      call wf%mem%dealloc(c, wf%n_parameters, 1)
!
   end subroutine initialize_trial_vectors_core_ionization_ccs
!
!
   module subroutine precondition_residual_valence_ionization_ccs(wf, residual)
!!
!!    Precondition residual valence
!!    Written by Sarai D. Folkestad, Aug. 2017
!!
!!    Projects out contamination from regular excitations
!!    Divide elements of residual by orbital difference       
!!
      implicit none
!
      class(ccs) :: wf
      real(dp), dimension(wf%n_parameters ,1) :: residual
!
      integer(i15) :: i = 0
!
      real(dp), dimension(:,:), allocatable :: orbital_diff      
!  
      call wf%ionization_residual_projection(residual)
!
      call wf%mem%alloc(orbital_diff, wf%n_parameters, 1)
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
      call wf%mem%dealloc(orbital_diff, wf%n_parameters, 1)
!
   end subroutine precondition_residual_valence_ionization_ccs
!
!
   module subroutine ionization_residual_projection_ccs(wf, residual)
!!
!!    Residual projection for core ionization (CCS)
!!    Written by Sarai D. Folkestad, Aug. 2017
!!
!!    Projects out contaminations from regular excitations
!!
      implicit none
!
      class(ccs) :: wf
      real(dp), dimension(wf%n_parameters, 1) :: residual
!
      integer(i15) :: i = 0, a = 0, diffuse_mo = 0, ai = 0, mo = 0
!
!     Find MO corresponding to super diffuse orbital
!  
!     Super diffuse AO is the last AO (for dalton integrals and HF this corresponds to 
!     the ghost atom being declared at the end of the MOLECULE.INP file). Must locate the MO
!     with weight 1.0D0 on this AO.
!
      do mo = 1, wf%n_mo
!  
         i = index_two(wf%n_ao, mo, wf%n_ao)
         if (abs(wf%mo_coef(i,1) - 1.0D0) .lt. 1.0D-7) diffuse_mo = mo
!
      enddo
!
!     Project out elements not corresponding to ionization
!
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            if (a .ne. (diffuse_mo - wf%n_o)) then
               ai = index_two(a, i, wf%n_v)
               residual(ai, 1) = 0.0d0
!
            endif
!
         enddo
      enddo
!
   end subroutine ionization_residual_projection_ccs
!
!
   module subroutine ionization_rho_a_i_projection_ccs(wf, rho_a_i)
!!
!!    Ionization rho_a_i projection (CCS).
!!    Written by Sarai D. Folkestad, Aug. 2017
!!
!!    Transformed vector projection for ionization, contamination from 
!!    regular excitations are projected out.
!!
      implicit none
!
      class(ccs) :: wf
      real(dp), dimension(wf%n_v, wf%n_o) :: rho_a_i
!
      integer(i15) :: i = 0, a = 0, diffuse_mo = 0, ai = 0, mo = 0
!
!     Find MO corresponding to super diffuse orbital
!  
!     Super diffuse AO is the last AO (for dalton integrals and HF this corresponds to 
!     the ghost atom being declared at the end of the MOLECULE.INP file). Must locate the MO
!     with weight 1.0D0 on this AO.
!
      do mo = 1, wf%n_mo
!  
         i = index_two(wf%n_ao, mo, wf%n_ao)
         if (abs(wf%mo_coef(i,1) - 1.0D0) .lt. 1.0D-7) diffuse_mo = mo
!
      enddo
!
!     Project out elements not corresponding to ionization
!
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            if (a .ne. (diffuse_mo - wf%n_o)) then
               rho_a_i(a, i) = 0.0d0
!
            endif
!
         enddo
      enddo
!
   end subroutine ionization_rho_a_i_projection_ccs
!
!
   module subroutine precondition_residual_core_ionization_ccs(wf, residual)
!  
!     Precondition residual core ionization
!     Written by Sarai D. Folkestad, Aug. 2017
!  
!     Projects out contaminations from valence contributions and all regular. excitations.
!     Divide elements of residual by corresponding orbital difference.       
!  
      implicit none
!
      class(ccs) :: wf
      real(dp), dimension(wf%n_parameters ,1) :: residual
!
      integer(i15) :: i = 0
!
      real(dp), dimension(:,:), allocatable :: orbital_diff      
!  
!     :: Projection ::
!
!     Project out valence
!
      call wf%cvs_residual_projection(residual)
!
!     Project out excitations
!
      call wf%ionization_residual_projection(residual)
!
!     Preconditioning by dividing with orbital differences
!
      call wf%mem%alloc(orbital_diff, wf%n_parameters, 1)
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
      call wf%mem%dealloc(orbital_diff, wf%n_parameters, 1)
!
   end subroutine precondition_residual_core_ionization_ccs
!
!
end submodule ionized_state