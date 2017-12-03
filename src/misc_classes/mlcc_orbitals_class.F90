module mlcc_orbitals_class
!
!! 
!
!  :::::::::::::::::::::::::::::::::::
!  -::- Modules used by the class -::-
!  :::::::::::::::::::::::::::::::::::
!
   use types
   use input_reader
!
!  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  -::- Definition of the mlcc_orbitals class -::-
!  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
   type :: mlcc_orbitals
!
!  Orbitals used
!
   logical :: cholesky = .false.
   logical :: cnto     = .true.
!
   real(dp) :: delta_o = 1.0D-06
   real(dp) :: delta_v = 1.0D-06
!
!  Cholesky specific variables
!
   integer(i15)                              :: n_active_atoms
   integer(i15), dimension(:,:), allocatable :: active_atoms
!
   contains
!
!  -::- Procedure pointers -::-
!
   procedure :: orbital_reader => orbital_reader_mlcc_orb
!
   end type mlcc_orbitals                                                                             
!
contains
! 
!  -::- Procedure definitions -::-
! 
   subroutine orbital_reader_mlcc_orb(orbital_info, unit_input)
!!
!!
      implicit none
!
!
      integer(i15)         :: unit_input
!
      class(mlcc_orbitals) :: orbital_info
!
      character(len=40)    :: line
!
      read(unit_input,'(a40)') line 
      line = remove_preceding_blanks(line)
!
      if (trim(line) == '{') then
!
         do
!
            read(unit_input,'(a40)') line 
            line = remove_preceding_blanks(line)
!
            if (trim(line) == 'cnto') then ! CNTOs
!
               read(unit_input,'(a40)') line 
               line = remove_preceding_blanks(line)
!
               if (trim(line) == '{') then
!
                  do
!
                     read(unit_input,'(a40)') line 
                     line = remove_preceding_blanks(line)
!
                     if (trim(line) == 'delta_o:') then
                        read(unit_input,*) orbital_info%delta_o 
                        cycle
                     elseif (trim(line) == 'delta_v:') then
                        read(unit_input,*) orbital_info%delta_v
                        cycle
                     elseif (trim(line) == '}') then
                        return
                     elseif (trim(line) == 'end of eT input') then
                        write(unit_output,*)'Error: end of CNTO section of eT.inp was not found.'
                        stop
                     endif
!
                  enddo
!
               endif ! Use defaults (cntos with both thresholds at 1.0D-6)
               return
!
            elseif (trim(line) == 'cholesky') then ! cholesky orbitals
!
               orbital_info%cnto = .false.
               orbital_info%cholesky = .true.
!
               read(unit_input,'(a40)') line 
               line = remove_preceding_blanks(line)
!
               if (trim(line) == '{') then
!
                  do
!
                     read(unit_input,'(a40)') line 
                     line = remove_preceding_blanks(line)
!
                     if (trim(line) == 'n_active_atoms:') then
!
                        read(unit_input,*) orbital_info%n_active_atoms 
                        cycle
!
                     elseif (trim(line) == 'active_atoms:') then
!
                        if (orbital_info%n_active_atoms == 0) then
                           write(unit_output,*)&
                           'Error: n_active_atoms must be specified in eT.inp before active atoms are selected'
                           stop
                        endif
!
                        allocate(orbital_info%active_atoms(orbital_info%n_active_atoms, 1))
                        read(unit_input,*) orbital_info%active_atoms
                        cycle
!
                     elseif (trim(line) == '}') then
!
!                       Sanity check for cholesky orbitals 
!
                           if (orbital_info%n_active_atoms == 0 .or. .not. allocated(orbital_info%active_atoms)) then
                              write(unit_output,*)'Error: information needed for Cholesky orbitals was not found in eT.inp'
                              stop
                           endif
!
                        return
                     elseif (trim(line) == 'end of eT input') then
                        write(unit_output,*)'Error: end of cholesky section of eT.inp was not found.'
                        stop
                     endif
                  enddo
!
                  return
               else 
                  write(unit_output,*)'Error: Cholesky orbitals requires specification of the active space!'
                  stop
               endif
               
            elseif (trim(line) == '}') then
               return
            elseif (trim(line) == 'end of eT input.') then
               write(unit_output,*)'Error: Could not find end of orbital section of eT.inp!'
               stop
            endif

         enddo
!
      else ! Use defaults (cntos with both thresholds at 1.0D-6)
         return
      endif
!
   end subroutine orbital_reader_mlcc_orb
!
!
end module mlcc_orbitals_class