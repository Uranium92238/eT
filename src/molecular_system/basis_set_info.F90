module basis_set_info
!
!!
!!    Basis set info module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!!
!
   use kinds
   use disk_manager_class
!
contains
!
   subroutine get_shells_to_fill_H_to_He(basis, info)
!!
!!
!!
      implicit none
!
      character(len=100) :: basis
!
      integer(i15), dimension(1, 2) :: info ! 1s
!
      info(1, 2) = 1 ! fill 1s first
!
      if (trim(basis) == 'cc-pVDZ' .or. trim(basis) == 'aug-cc-pVDZ') then ! DZ -> 2*(1s)
!
         info(1, 1) = 2
!
      elseif (trim(basis) == 'cc-pVTZ' .or. trim(basis) == 'aug-cc-pVTZ') then ! TZ -> 3*(1s)
!
         info(1, 1) = 3
!
      elseif (trim(basis) == 'cc-pVQZ' .or. trim(basis) == 'aug-cc-pVQZ') then ! QZ -> 4*(1s)
!
         info(1, 1) = 4
!
      elseif (trim(basis) == 'cc-pV5Z' .or. trim(basis) == 'aug-cc-pV5Z') then ! 5Z -> 5*(1s)
!
         info(1, 1) = 5
!
      elseif (trim(basis) == 'cc-pV6Z' .or. trim(basis) == 'aug-cc-pV6Z') then ! 6Z -> 6*(1s)
!
         info(1, 1) = 6
!
      elseif (trim(basis) == 'sto-3g') then ! TZ -> 1*(1s)
!
         info(1, 1) = 1
!
      else
!
         write(output%unit, *)'Error: the selected basis is not yet supported by eT.'
         stop
!
      endif
!
   end subroutine get_shells_to_fill_H_to_He
!
!
   subroutine get_shells_to_fill_Li_to_Ne(basis, info)
!!
!!
!!
      implicit none
!
      character(len=100) :: basis
!
      integer(i15), dimension(3, 2) :: info ! 1s, 2s, 2p
!
      info(1, 1) = 1
!
      info(1, 2) = 1 ! fill 1s first
      info(2, 2) = 2 ! 2s second
      info(3, 2) = 3 ! 2p last
!
      if (trim(basis) == 'cc-pVDZ' .or. trim(basis) == 'aug-cc-pVDZ') then ! DZ -> 1*(1s), 2*(2s,2p)
!
         info(2, 1) = 2
         info(3, 1) = 2
!
      elseif (trim(basis) == 'cc-pVTZ' .or. trim(basis) == 'aug-cc-pVTZ') then ! TZ -> 1*(1s), 3*(2s,2p)
!
         info(2, 1) = 3
         info(3, 1) = 3
!
      elseif (trim(basis) == 'cc-pVQZ' .or. trim(basis) == 'aug-cc-pVQZ') then ! QZ -> 1*(1s), 4*(2s,2p)
!
         info(2, 1) = 4
         info(3, 1) = 4
!
      elseif (trim(basis) == 'cc-pV5Z' .or. trim(basis) == 'aug-cc-pV5Z') then ! 5Z -> 1*(1s), 5*(2s,2p)
!
         info(2, 1) = 5
         info(3, 1) = 5
!
      elseif (trim(basis) == 'cc-pV6Z' .or. trim(basis) == 'aug-cc-pV6Z') then ! 6Z -> 1*(1s), 6*(2s,2p)
!
         info(2, 1) = 6
         info(3, 1) = 6
!
      elseif (trim(basis) == 'sto-3g') then ! SZ -> 1*(1s)
!
         info(2, 1) = 1
         info(3, 1) = 1
!
      else
!
         write(output%unit, *)'Error: the selected basis is not yet supported by eT.'
         stop
!
      endif
!
   end subroutine get_shells_to_fill_Li_to_Ne
!
!
   subroutine get_shells_to_fill_Na_to_Ar(basis, info)
!!
!!
!!
      implicit none
!
      character(len=100) :: basis
!
      integer(i15), dimension(5, 2) :: info !1s, 2s, 3s, 2p, 3p
!
!     Fill order
!
      info(1, 2) = 1 ! 1s
      info(2, 2) = 2 ! 2s
      info(4, 2) = 3 ! 2p 
      info(3, 2) = 4 ! 3s
      info(5, 2) = 5 ! 3p
!
      if (trim(basis) == 'cc-pVDZ' .or. trim(basis) == 'aug-cc-pVDZ') then ! DZ -> 1*(1s,2s,2p), 2*(3s,3p)
!
         info(1,1) = 1
         info(2,1) = 1
         info(4,1) = 1
!
         info(3,1) = 2
         info(5,1) = 2
!
      elseif (trim(basis) == 'cc-pVTZ' .or. trim(basis) == 'aug-cc-pVTZ') then ! TZ -> 1*(1s,2s,2p), 3*(3s,3p)
!
         info(1,1) = 1
         info(2,1) = 1
         info(4,1) = 1
!
         info(3,1) = 3
         info(5,1) = 3
!
      elseif (trim(basis) == 'cc-pVQZ' .or. trim(basis) == 'aug-cc-pVQZ') then ! QZ -> 1*(1s,2s,2p), 4*(3s,3p)
!
         info(1,1) = 1 ! 1s
         info(2,1) = 1 ! 2s
         info(4,1) = 1 ! 2p
!
         info(3,1) =4
         info(5,1) =4
!
      elseif (trim(basis) == 'cc-pV5Z' .or. trim(basis) == 'aug-cc-pV5Z') then ! 5Z -> 1*(1s,2s,2p), 5*(3s,3p)
!
         info(1,1) = 1 
         info(2,1) = 1 
         info(4,1) = 1 
!
         info(3,1) = 5
         info(5,1) = 5
!
      elseif (trim(basis) == 'cc-pV6Z' .or. trim(basis) == 'aug-cc-pV6Z') then ! 6Z -> 1*(1s,2s,2p), 6*(3s,3p)
!
         info(1,1) = 1
         info(2,1) = 1
         info(4,1) = 1
!
         info(3,1) = 6
         info(5,1) = 6
!
      elseif (trim(basis) == 'sto-3g') then ! SZ -> 1*(1s,2s,2p,3s,3p)
!
         info(1, 1) = 1
         info(2, 1) = 1
         info(3, 1) = 1
         info(4, 1) = 1
         info(5, 1) = 1
!
      else
!
         write(output%unit, *)'Error: the selected basis is not yet supported by eT.'
         stop
!
      endif
!
   end subroutine get_shells_to_fill_Na_to_Ar
!
!
!
!
   subroutine get_shells_to_fill_K_to_Ca(basis, info)
!!
!!
!!
      implicit none
!
      character(len=100) :: basis
!
      integer(i15), dimension(6, 2) :: info !1s, 2s, 3s, 4s, 2p, 3p
!
!     Fill order
!
      info(1, 2) = 1 ! 1s
      info(2, 2) = 2 ! 2s
      info(5, 2) = 3 ! 2p 
      info(3, 2) = 4 ! 3s
      info(6, 2) = 5 ! 3p
      info(4, 2) = 6 ! 4s
!
      if (trim(basis) == 'cc-pVDZ' .or. trim(basis) == 'aug-cc-pVDZ') then ! DZ -> 1*(1s,2s,2p,3s,3p), 2*(4s)
!
         info(1,1) = 1
         info(2,1) = 1
         info(3,1) = 1
         info(5,1) = 1
         info(6,1) = 1
!
         info(4,1) = 2
!
      elseif (trim(basis) == 'cc-pVTZ' .or. trim(basis) == 'aug-cc-pVTZ') then ! TZ -> 1*(1s,2s,2p,3s,3p), 3*(4s)
!
         info(1,1) = 1
         info(2,1) = 1
         info(3,1) = 1
         info(5,1) = 1
         info(6,1) = 1
!
         info(4,1) = 3
!
      elseif (trim(basis) == 'cc-pVQZ' .or. trim(basis) == 'aug-cc-pVQZ') then ! QZ -> 1*(1s,2s,2p,3s,3p), 4*(4s)
!
         info(1,1) = 1
         info(2,1) = 1
         info(3,1) = 1
         info(5,1) = 1
         info(6,1) = 1
!
         info(4,1) = 4
!
      elseif (trim(basis) == 'cc-pV5Z' .or. trim(basis) == 'aug-cc-pV5Z') then ! 5Z -> 1*(1s,2s,2p,3s,3p), 5*(4s)
!
         info(1,1) = 1
         info(2,1) = 1
         info(3,1) = 1
         info(5,1) = 1
         info(6,1) = 1
!
         info(4,1) = 5
!
      elseif (trim(basis) == 'cc-pV6Z' .or. trim(basis) == 'aug-cc-pV6Z') then ! 6Z -> 1*(1s,2s,2p,3s,3p), 6*(4s)
!
         info(1,1) = 1
         info(2,1) = 1
         info(3,1) = 1
         info(5,1) = 1
         info(6,1) = 1
!
         info(4,1) = 6
!
      elseif (trim(basis) == 'sto-3g') then ! SZ -> 1*(1s,2s,2p,3s,3p,4s)
!
         info(1, 1) = 1
         info(2, 1) = 1
         info(3, 1) = 1
         info(4, 1) = 1
         info(5, 1) = 1
         info(6, 1) = 1
!
      else
!
         write(output%unit, *)'Error: the selected basis is not yet supported by eT.'
         stop
!
      endif
!
   end subroutine get_shells_to_fill_K_to_Ca
end module basis_set_info
