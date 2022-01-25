!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2022 the authors of eT
!
!  eT is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  eT is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program. If not, see <https://www.gnu.org/licenses/>.
!
!
module eT_references_class
!
!!
!!    eT references class
!!    Written by Alexander C. Paul, Jan 2022
!!
!
   use citation_class, only: citation
   use global_out, only: output
!
   implicit none
!
   type :: eT_references
!
      integer :: n_references
!
      type(citation),    dimension(:), allocatable :: references
      character(len=40), dimension(:), allocatable :: labels
!
   contains
!
      procedure :: get_doi_citation
      procedure :: get_apa_citation
!
   end type eT_references
!
   interface eT_references
!
      procedure :: new_eT_references
!
   end interface eT_references
!
contains
!
!
   function new_eT_references() result(this)
!!
!!    New eT references
!!    Written by Alexander C. Paul, Jan 2022
!!
      implicit none
!
      type(eT_references) :: this
      integer :: i
!
      this%n_references = 7
      allocate(this%references(this%n_references))
      allocate(this%labels(this%n_references))
!
      this%references(1) = citation(implementation = 'eT', &
                                    journal        = 'J. Chem. Phys.', &
                                    title_         = 'eT 1.0: An open source electronic structure &
                                                      &program with emphasis on coupled cluster &
                                                      &and multilevel methods', &
                                    volume         = '152', &
                                    issue          = '18', &
                                    pages          = '184103', &
                                    year           = '2020', &
                                    doi            = '10.1063/5.0004713', &
                                    authors        = [character(len=25) :: 'Sarai D. Folkestad', &
                                                                            'Eirik F. Kjønstad', &
                                                                           'Rolf H. Myhre', &
                                                                           'Josefine H. Andersen', &
                                                                           'Alice Balbi', &
                                                                           'Sonia Coriani', &
                                                                           'Tommaso Giovannini', &
                                                                           'Linda Goletto', &
                                                                           'Tor S. Haugland', &
                                                                           'Anders Hutcheson', &
                                                                           'Ida-Marie Høyvik', &
                                                                           'Torsha Moitra', &
                                                                           'Alexander C. Paul', &
                                                                           'Marco Scavino', &
                                                                           'Andreas S. Skeidsvoll', &
                                                                           'Åsmund H. Tveten', &
                                                                           'Henrik Koch'])
!
!
      this%references(2) = citation(implementation = 'Multimodel Newton algorithm', &
                                    journal        = 'J. Chem. Phys.', &
                                    title_         = 'Accelerated multimodel Newton-type algorithms &
                                                     &for faster convergence of ground and excited &
                                                     &state coupled cluster equations', &
                                    volume         = '153', &
                                    issue          = '1', &
                                    pages          = '014104', &
                                    year           = '2020', &
                                    doi            = '10.1063/5.0010989', &
                                    authors        = [character(len=25) :: 'Eirik F. Kjønstad', &
                                                                           'Sarai D. Folkestad', &
                                                                           'Henrik Koch'])
!
!
      this%references(3) = citation(implementation = 'Cholesky decomposition of ERIs', &
                                    journal        = 'J. Chem. Phys.', &
                                    title_         = 'An efficient algorithm for Cholesky &
                                                     &decomposition of electron repulsion integrals', &
                                    volume         = '150', &
                                    issue          = '19', &
                                    pages          = '194112', &
                                    year           = '2019', &
                                    doi            = '10.1063/1.5083802', &
                                    authors        = [character(len=25) :: 'Sarai D. Folkestad', &
                                                                           'Eirik F. Kjønstad', &
                                                                           'Henrik Koch'])
!
!
      this%references(4) = citation(implementation = 'Time-dependent CC', &
                                    journal        = 'Phys. Rev. A', &
                                    title_         = 'Time-dependent coupled-cluster theory for &
                                                     &ultrafast transient-absorption spectroscopy', &
                                    volume         = '102', &
                                    issue          = '2', &
                                    pages          = '023115', &
                                    year           = '2020', &
                                    doi            = '10.1103/PhysRevA.102.023115', &
                                    authors        = [character(len=25) :: 'Andreas S. Skeidsvoll', &
                                                                           'Alice Balbi', &
                                                                           'Henrik Koch'])
!
      this%references(5) = citation(implementation = 'CC3', &
                                    journal        = 'J. Chem. Theory Comput.', &
                                    title_         = 'New and Efficient Implementation of CC3', &
                                    volume         = '17', &
                                    issue          = '1', &
                                    pages          = '117–126', &
                                    year           = '2021', &
                                    doi            = '10.1021/acs.jctc.0c00686', &
                                    authors        = [character(len=25) :: 'Alexander C. Paul', &
                                                                           'Rolf H. Myhre', &
                                                                           'Henrik Koch'])
!
      this%references(6) = citation(implementation = 'MLCC2 and MLCCSD', &
                                    journal        = 'J. Chem. Theory Comput.', &
                                    title_         = 'Multilevel CC2 and CCSD Methods with &
                                                     &Correlated Natural Transition Orbitals', &
                                    volume         = '16', &
                                    issue          = '1', &
                                    pages          = '179–189', &
                                    year           = '2019', &
                                    doi            = '10.1021/acs.jctc.9b00701', &
                                    authors        = [character(len=25) :: 'Sarai Dery Folkestad', &
                                                                           'Henrik Koch'])
!
      this%references(7) = citation(implementation = 'QED-HF', &
                                    journal        = 'Phys. Rev. X', &
                                    title_         = 'Coupled Cluster Theory for Molecular Polaritons: &
                                                     &Changing ground and Excited States', &
                                    volume         = '10', &
                                    issue          = '4', &
                                    pages          = '041043', &
                                    year           = '2020', &
                                    doi            = '10.1103/PhysRevX.10.041043', &
                                    authors        = [character(len=25) :: 'Tor S. Haugland', &
                                                                           'Enrico Ronca', &
                                                                           'Eirik F. Kjønstad', &
                                                                           'Angel Rubio', &
                                                                           'Henrik Koch'])
!
      do i = 1, this%n_references
         this%labels(i) = this%references(i)%get_implementation()
      end do
!
   end function new_eT_references
!
!
   function get_doi_citation(this, label) result(doi_citation)
!!
!!    Get DOI citation
!!    Written by Alexander C. Paul, Jan 2022
!!
      use array_utilities, only: find_location
!
      implicit none
!
      class(eT_references), intent(in) :: this
!
      character(len=:), allocatable :: doi_citation
!
      character(len=40), intent(in)  :: label
!
      integer :: index
!
      index = find_location(this%labels, label, this%n_references)
!
      if (index == 0) then
         call output%error_msg("Citation with label (a0) not found in eT_references.F90", &
                               chars=[trim(label)])
      end if
!
      doi_citation = this%references(index)%get_doi_citation_string()
!
   end function get_doi_citation
!
!
   function get_apa_citation(this, label) result(apa_citation)
!!
!!    Get apa citation
!!    Written by Alexander C. Paul, Jan 2022
!!
      use array_utilities, only: find_location
!
      implicit none
!
      class(eT_references), intent(in) :: this
!
      character(len=:), allocatable :: apa_citation
!
      character(len=40), intent(in)  :: label
!
      integer :: index
!
      index = find_location(this%labels, label, this%n_references)
!
      if (index == 0) then
         call output%error_msg("Citation with label (a0) not found in eT_references.F90", &
                               chars=[trim(label)])
      end if
!
      apa_citation = this%references(index)%get_apa_citation_string()
!
   end function get_apa_citation
!
!
end module eT_references_class
