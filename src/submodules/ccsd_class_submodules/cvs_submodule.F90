submodule(ccsd_class) cvs
!
!!
!!    CVS submodule (CCSD)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!    Contains the CVS routines for the CCSD class.
!!    Note that this submodule contains both excited state routines and jacobian transformation routines.
!!
!
contains
!
!
   module subroutine cvs_residual_projection_ccsd(wf, residual)
!!
!!    Residual projection (CCSD), 
!!    Written by Sarai D. Folkestad Aug. 2017    
!!
      implicit none
!
      class(ccsd) :: wf
      real(dp), dimension(wf%n_parameters, 1) :: residual
!
      integer(i15) :: i = 0, a = 0, j = 0, b = 0, core = 0, ai = 0, bj = 0, aibj = 0
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
      do i = 1, wf%n_o
         do j = 1, wf%n_o
!
            core_orbital = .false.
            do core = 1, wf%core_excited_state_specifications%n_equivalent_cores
!
               if ((i .eq. wf%core_excited_state_specifications%index_core_mo(core, 1)) .or. &
                  (j .eq. wf%core_excited_state_specifications%index_core_mo(core, 1))) core_orbital = .true.
!
            enddo
!
            if (.not. core_orbital) then
               do a = 1,  wf%n_v
                  do b = 1,  wf%n_v
                     ai = index_two(a, i,  wf%n_v)
                     bj = index_two(b, j,  wf%n_v)
                     aibj = index_packed(ai, bj)

                     residual(wf%n_t1am + aibj, 1) = zero
                  enddo
               enddo
            endif
         enddo
      enddo
!
   end subroutine cvs_residual_projection_ccsd
!
!
    module subroutine cvs_rho_aibj_projection_ccsd(wf, vec_aibj)
!!
!!    Rho projection for CVS (CCSD),
!!    Written by Sarai D. Folkestad, Aug. 2017
!!
      implicit none
!
      class(ccsd) :: wf
      real(dp), dimension(:, :) :: vec_aibj
!
      integer(i15) :: i = 0, a = 0, j = 0, b = 0, core = 0, ai = 0, bj = 0, aibj = 0
!
      logical :: core_orbital      
!
      do i = 1, wf%n_o
       do j = 1, wf%n_o
!
          core_orbital = .false.
          do core = 1, wf%core_excited_state_specifications%n_equivalent_cores
!
             if ((i .eq. wf%core_excited_state_specifications%index_core_mo(core, 1)) .or. &
                (j .eq. wf%core_excited_state_specifications%index_core_mo(core, 1))) core_orbital = .true.
!
          enddo
!
          if (.not. core_orbital) then
             do a = 1, wf%n_v
                do b = 1, wf%n_v
!
                  ai = index_two(a, i, wf%n_v)
                  bj = index_two(b, j, wf%n_v)
!
                  if (ai .ge. bj) then
                    aibj = index_packed(ai, bj)
!
                    vec_aibj(aibj, 1) = zero
!
                  endif
                enddo
             enddo
          endif
       enddo
    enddo
!
  end subroutine cvs_rho_aibj_projection_ccsd
!
!
end submodule cvs
