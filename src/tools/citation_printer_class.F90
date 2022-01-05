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
module citation_printer_class
!
!!
!!    Citation printer
!!    Written by Eirik F. Kjønstad, 2021 
!!
!!    Class for a citation printer.
!!
!!    Usage:
!!
!!       citation = new_citation(...)
!!       call citation_printer%add(citation)
!!
!!       Add more citations if you wish...
!!
!!       call citation_printer%print_(file)
!!        
!
   use citation_class,     only: citation
   use output_file_class,  only: output_file
   use global_in,          only: input
!
   implicit none
!
   type :: citation_printer
!
      integer, private :: n_citations 
      type(citation), dimension(100), private :: citations ! Assumes that no calculation will ask  
                                                           ! the user to cite more than 100 paper
!
      logical :: print_apa ! If true, prints full APA-style references;
                           ! Otherwise, only the DOI is printed.
!
   contains
!
      procedure, public :: add & 
                        => add_citation_printer
!
      procedure, public :: print_ &
                        => print_citation_printer
!
   end type citation_printer
!
   interface citation_printer
!
      procedure :: new_citation_printer
!
   end interface citation_printer
!
!  Main citation object
!
   type(citation_printer) :: eT_citations
!
contains
!
!
   function new_citation_printer() result(this)
!!
!!    New citation printer
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none
!
      type(citation_printer) :: this
!
      class(citation), allocatable :: eT_reference
!
      eT_reference = citation(implementation = 'eT',                                        &
                              journal        = 'J. Chem. Phys.',                                &
                              title_         = 'eT 1.0: An open source electronic structure &
                                                &program with emphasis on coupled cluster &
                                                &and multilevel methods',                       &
                              volume         = '152',                                           &
                              issue          = '18',                                            &
                              pages          = '184103',                                        &
                              year           = '2020',                                          &
                              doi            = '10.1063/5.0004713',                             &
                              authors        = [character(len=25) :: 'Sarai D. Folkestad',      &
                                                                     'Eirik F. Kjønstad',       &
                                                                     'Rolf H. Myhre',           &
                                                                     'Josefine H. Andersen',    &
                                                                     'Alice Balbi',             &
                                                                     'Sonia Coriani',           &
                                                                     'Tommaso Giovannini',      &
                                                                     'Linda Goletto',           &
                                                                     'Tor S. Haugland',         & 
                                                                     'Anders Hutcheson',        &
                                                                     'Ida-Marie Høyvik',        &
                                                                     'Torsha Moitra',           &
                                                                     'Alexander C. Paul',       &
                                                                     'Marco Scavino',           &
                                                                     'Andreas S. Skeidsvoll',   &
                                                                     'Åsmund H. Tveten',        &
                                                                     'Henrik Koch'])
!
      call this%add(eT_reference)
!
      this%n_citations = 1 ! Initially, there is just the eT reference
!
      this%print_apa = input%is_keyword_present('full references', 'print')
!
   end function new_citation_printer
!
!
   subroutine add_citation_printer(this, citation_)
!!
!!    Add
!!    Written by Eirik F. Kjønstad, 2021
!!
!!    Adds a citation.
!!
      implicit none 
!
      class(citation_printer), intent(inout) :: this 
!
      class(citation), intent(in) :: citation_ 
!
      this%n_citations = this%n_citations + 1
!
      this%citations(this%n_citations) = citation_
!
   end subroutine add_citation_printer
!
!
   subroutine print_citation_printer(this, file)
!!
!!    Print 
!!    Written by Eirik F. Kjønstad, 2021
!!
!!    Prints all citations to the output file named 'file'
!!
      implicit none 
!
      class(citation_printer), intent(inout) :: this 
!
      class(output_file), intent(inout) :: file 
!
      integer :: k 
!
      character(len=600) :: apa_citation, doi, implementation, doi_citation
!
      call file%printf('m', '- Implementation references:', fs='(/t3,a/)')
!
      do k = 1, this%n_citations
!
         doi            = this%citations(k)%get_doi()
         implementation = this%citations(k)%get_implementation()
!
         doi_citation = trim(implementation) // ": https://doi.org/" // trim(doi)
!
         call file%printf('m', doi_citation, fs='(t6,a)')
!
         if (this%print_apa) then 
!
            apa_citation = this%citations(k)%get_apa_citation_string()
!
            call file%printf('m', apa_citation, ffs='(/t9,a)', fs='(t9,a)')
!
            if (k .lt. this%n_citations) call file%newline('m')
!
         endif
!
      enddo
!
   end subroutine print_citation_printer
!
!
end module citation_printer_class
