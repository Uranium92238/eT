!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2020 the authors of eT
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
module approximate_jacobian_transformer_class
!!
!!    Approximate Jacobian transformer
!!    Written by Eirik F. Kjønstad, 2021
!! 
!!    Wrapper object for calling the most accurate 
!!    approximation of the full Jacobian transformation
!!    that scales less steeply than the full transformation 
!!    (e.g., the CCSD transformation in the case of CC3).
!!
   use kinds
   use ccs_class,                            only: ccs 
   use abstract_jacobian_transformer_class,  only: abstract_jacobian_transformer
!
   implicit none
!
   type, extends(abstract_jacobian_transformer) :: approximate_jacobian_transformer
!
   contains
!
      procedure :: transform => transform_approximate_jacobian_transformer
      procedure :: prepare   => prepare_approximate_jacobian_transformer
!
   end type approximate_jacobian_transformer
!
!
   interface approximate_jacobian_transformer
!
      procedure :: new_approximate_jacobian_transformer
!
   end interface approximate_jacobian_transformer
!
!
contains
!
!
   pure function new_approximate_jacobian_transformer(side) result(this)
!!
!!    New approximate Jacobian transformer
!!    Written by Eirik F. Kjønstad, 2021
!!
!!    side: 'left' or 'right' (call A^T or A transformation, respectively)
!!
      implicit none 
!
      character(len=*), intent(in) :: side
!
      type(approximate_jacobian_transformer) :: this 
!
      this%side = trim(side)
!
   end function new_approximate_jacobian_transformer
!
!
   subroutine prepare_approximate_jacobian_transformer(this, wf)
!!
!!    Prepare
!!    Written by Eirik F. Kjønstad, 2021
!!
!!    Wrapper for preparation of approximate Jacobian transformation.
!!
!!    Calls the required routine in the wavefunction.
!!
      implicit none
!
      class(approximate_jacobian_transformer), intent(in) :: this
!
      class(ccs), intent(inout) :: wf
!
      call wf%prepare_for_approximate_Jacobians(this%side)
!
   end subroutine prepare_approximate_jacobian_transformer
!
!
   subroutine transform_approximate_jacobian_transformer(this, wf, c, x, w)
!!
!!    Transform
!!    Written by Eirik F. Kjønstad, 2021
!!
!!    Wrapper for approximate Jacobian transformation.
!!
!!    c: untransformed vector
!!    x: transformed vector
!!    w: excitation energy
!!
!!    Calls the required routine in the wavefunction.
!!
      implicit none
!
      class(approximate_jacobian_transformer), intent(in) :: this
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in)  :: c 
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: x 
!
      real(dp), intent(in) :: w 
!
      call wf%approximate_Jacobian_transform(this%side, c, x, w)
!
   end subroutine transform_approximate_jacobian_transformer
!
!
end module approximate_jacobian_transformer_class
