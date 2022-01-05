//
//
//  eT - a coupled cluster program
//  Copyright (C) 2016-2022 the authors of eT
//
//  eT is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  eT is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program. If not, see <https://www.gnu.org/licenses/>.
//
// -----------------------------------------------------------------------
/*

   Atom initialization routines
   Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018

*/
#include "ri_basis.h"

void get_n_ri_shells(int *n_shells){

   n_shells[0] = int(ri_basis.size());

}

void get_n_ri_ao(int *n_ao){

   n_ao[0] = int(ri_basis.nbf());

}

void get_ri_shell_size(int *shell, int *size){

   auto n = ri_basis[*shell - 1].size();
   size[0] = int(n);

}
