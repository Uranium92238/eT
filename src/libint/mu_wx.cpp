//
//
//  eT - a coupled cluster program
//  Copyright (C) 2016-2019 the authors of eT
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
 
  	Construct dipole integrals mu
  	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
 
*/
#include <iostream>

using namespace std;

#include <fstream>
#include <string>
#include <vector>
#include "mu_wx.h"
#include "extract_integrals.h"

#include <libint2.hpp>

#include "globals.h"

using namespace libint2;

void construct_ao_mu_wx(double *mu_X, double *mu_Y, double *mu_Z, int *s1, int *s2){

  const auto& buf_vec = dipole.results(); // will point to computed shell sets

  auto n1 = basis[*s1 - 1].size();        // Number of basis functions in shell 1
  auto n2 = basis[*s2 - 1].size();        // number of basis functions in shell 2

  dipole.compute(basis[*s1 - 1], basis[*s2 - 1]);

  // I don't know why, but Libint computes the overlap integrals for some reason (they are in buf_vec[0])!
  // => We treat this as junk & never use it. - Eirik

  auto ints_shellset_X = buf_vec[1];      // location of the computed mu_X integrals
  auto ints_shellset_Y = buf_vec[2];      // location of the computed mu_Y integrals
  auto ints_shellset_Z = buf_vec[3];      // location of the computed mu_Z integrals

  extract_integrals(mu_X, ints_shellset_X, n1, n2, -1.0e0);
  extract_integrals(mu_Y, ints_shellset_Y, n1, n2, -1.0e0);
  extract_integrals(mu_Z, ints_shellset_Z, n1, n2, -1.0e0);

  return;

}
