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
  
   	Construct h
   	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
  
*/

#include <iostream>

using namespace std;

#include <fstream>
#include <string>
#include <vector>
#include <omp.h>

#include <libint2.hpp>

#include "h_wx.h"

#include "globals.h"

using namespace libint2;

void construct_ao_h_wx(double *h, int *s1, int *s2){
/*
/   Add kinetic contribution
*/
  int thread = omp_get_thread_num();

  const auto& buf_vec = kinetic[thread].results(); // will point to computed shell sets

  kinetic[thread].compute(basis[*s1 - 1], basis[*s2 - 1]);

  auto ints_shellset = buf_vec[0];                // location of the computed integrals

  auto n1 = basis[*s1 - 1].size();                // Number of basis functions in shell 1
  auto n2 = basis[*s2 - 1].size();                // number of basis functions in shell 2

  if (ints_shellset == nullptr) {
    for(auto f1=0; f1!=n1; ++f1){
      for(auto f2=0; f2!=n2; ++f2){

      *(h + n1*f2+f1) = 0.0e0;

      }
    }
  }
  else{
    for(auto f1=0; f1!=n1; ++f1){
      for(auto f2=0; f2!=n2; ++f2){

      *(h + n1*f2+f1) = ints_shellset[f1*n2+f2]; 

      }
    }
  }
/*
/   Add nuclear attraction contribution
*/
  const auto& buf_vec_n = nuclear[thread].results();      // will point to computed shell sets

  nuclear[thread].compute(basis[*s1 - 1], basis[*s2 - 1]);
  auto ints_shellset_n = buf_vec_n[0];                    // location of the computed integrals

  if (ints_shellset_n != nullptr){
    for(auto f1=0; f1!=n1; ++f1){
      for(auto f2=0; f2!=n2; ++f2){

      *(h + n1*f2+f1) = *(h + n1*f2+f1) + ints_shellset_n[f1*n2+f2];

      }
    }
  }

  return;
}
