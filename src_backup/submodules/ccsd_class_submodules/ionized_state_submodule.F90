submodule(ccsd_class) ionized_state
!
!!
!!    Ionized state submodule (CCSD)
!!    Written by Eirik F. Kjønstad and Sarai Dery Folkestad, October 2017
!!
!!    Contains the following family of procedures of the CCS class:
!!
!!    ionization_residual_projection: description missing
!!    ionization_rho_aibj_projection: -/-
!1
!
   implicit none 
!
!
contains
!
!
   module subroutine ionization_residual_projection_ccsd(wf, residual)
!!
!!    Residual projection for CVS
!!    Written by Sarai D. Folkestad, Aug. 2017
!!
      implicit none
!
      class(ccsd) :: wf
      real(dp), dimension(wf%n_parameters, 1) :: residual
!
      integer(i15) :: i = 0, a = 0, diffuse_mo = 0, ai = 0, mo = 0, b = 0, j = 0, bj = 0, aibj = 0
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
!
               ai = index_two(a, i, wf%n_v)
               residual(ai, 1) = 0.0d0
!
            endif
          enddo
       enddo
!
      do a = 1, wf%n_v
         do b = 1, wf%n_v
            do i = 1, wf%n_o
               do j = 1, wf%n_o
!
               if (.not.((a .eq. (diffuse_mo - wf%n_o) ) .or. (b .eq. (diffuse_mo - wf%n_o))) ) then
!
                  bj = index_two(b, j, wf%n_v)
                  aibj = index_packed(ai, bj)
                  residual((wf%n_o*wf%n_v) + aibj, 1) = 0.0d0
!
               endif
!
               enddo
            enddo
!
         enddo
      enddo
!
   end subroutine ionization_residual_projection_ccsd
!
!
   module subroutine ionization_rho_aibj_projection_ccsd(wf, rho_aibj)
!!
!!    Residual projection for CVS
!!    Written by Sarai D. Folkestad, Aug 2017
!!
      implicit none
!
      class(ccsd) :: wf
      real(dp), dimension(:,:) :: rho_aibj
!
      integer(i15) :: i = 0, a = 0, diffuse_mo = 0, ai = 0, mo = 0, b = 0, j = 0, bj = 0, aibj = 0
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
            do j = 1, wf%n_o
               do b = 1, wf%n_v
!
                  if ((a .ne. (diffuse_mo - wf%n_o) ).and. (b .ne. (diffuse_mo - wf%n_o)) ) then
!
                     bj = index_two(b, j, wf%n_v)
                     ai = index_two(a, i, wf%n_v)
!
                     if (ai .ge. bj) then
!
                        aibj=index_packed(ai,bj)
                        rho_aibj(aibj, 1) = 0.0d0
!
                     endif
!
                  endif
!
               enddo
            enddo
         enddo
      enddo
!
   end subroutine ionization_rho_aibj_projection_ccsd
!
!
end submodule ionized_state