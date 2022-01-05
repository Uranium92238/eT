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
module abstract_jacobian_transformer_class
!!
!!    Abstract Jacobian transformer
!!    Written by Eirik F. Kjønstad, 2021
!! 
!!    Wrapper classes for calling various Jacobian transformations in a single wave function.
!!    (e.g. lower-level transformations, such as CCSD transformation in a CC3 wavefunction)
!!    
   use kinds
   use ccs_class, only: ccs 
!
   implicit none
!
   type, abstract :: abstract_jacobian_transformer
!
      character(len=200) :: side ! 'right' or 'left' 
!
   contains
!
      procedure(prepare_for_jacobian), deferred :: prepare
      procedure(jacobian_transform), deferred   :: transform
!
   end type abstract_jacobian_transformer
!
   abstract interface
!
      subroutine prepare_for_jacobian(this, wf)
!!
!!       Prepare for Jacobian
!!       Written by Eirik F. Kjønstad, 2021
!!
!!       Interface for preparation of Jacobian transformation in descendants.
!!
!!       c: untransformed vector
!!       x: transformed vector
!!       w: excitation energy
!!
         import ccs
         import abstract_jacobian_transformer
!
         implicit none
!
         class(abstract_jacobian_transformer), intent(in) :: this
!
         class(ccs), intent(inout) :: wf
!
      end subroutine prepare_for_jacobian
!
      subroutine jacobian_transform(this, wf, c, x, w)
!!
!!       Jacobian transform
!!       Written by Eirik F. Kjønstad, 2021
!!
!!       Interface for Jacobian transformation in descendants.
!!
!!       c: untransformed vector
!!       x: transformed vector
!!       w: excitation energy
!!
         import ccs
         import abstract_jacobian_transformer
         import dp 
!
         implicit none
!
         class(abstract_jacobian_transformer), intent(in) :: this
!
         class(ccs), intent(inout) :: wf
!
         real(dp), dimension(wf%n_es_amplitudes), intent(in)  :: c 
         real(dp), dimension(wf%n_es_amplitudes), intent(out) :: x 
!
         real(dp), intent(in) :: w 
!
      end subroutine jacobian_transform
!
   end interface 
!
end module abstract_jacobian_transformer_class
