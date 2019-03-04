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
  
  	Construct electronic repulsion integrals g
  	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
  
*/
#include <iostream>

using namespace std;

#include <fstream>
#include <string>
#include <ctime>

#include "g_wxyz.h"

#include <libint2.hpp>
#include <vector>

#include "globals.h"

using namespace libint2;

void construct_ao_g_wxyz(double *g, int *s1, int *s2, int *s3, int *s4){

	int thread = omp_get_thread_num();

	const auto& buf_vec = electronic_repulsion_engines[thread].results(); // will point to computed shell sets

	auto n1 = basis[*s1 - 1].size(); // Number of basis functions in shell 1
	auto n2 = basis[*s2 - 1].size(); // Number of basis functions in shell 2
	auto n3 = basis[*s3 - 1].size(); // Number of basis functions in shell 3
	auto n4 = basis[*s4 - 1].size(); // Number of basis functions in shell 4

	electronic_repulsion_engines[thread].compute(basis[*s1 - 1], basis[*s2 - 1], basis[*s3 - 1], basis[*s4 - 1]);

	auto ints_1234 = buf_vec[0];     // Location of computed integrals

   if (ints_1234 == nullptr)
   {
      for(auto f1=0, f1234=0; f1!=n1; ++f1){

  	      for(auto f2=0; f2!=n2; ++f2){

            for(auto f3=0; f3!=n3; ++f3){

               for(auto f4=0; f4!=n4; ++f4, ++f1234){
                int ind_offset = n1*(n2*(n3*f4+f3)+f2)+f1;
	  				     *(g + ind_offset) = 0.0e0;

               }
            }
         }
      }
   }
   else
   {
      for(auto f1=0, f1234=0; f1!=n1; ++f1){

         for(auto f2=0; f2!=n2; ++f2){

            for(auto f3=0; f3!=n3; ++f3){

               for(auto f4=0; f4!=n4; ++f4, ++f1234){

                int ind_offset = n1*(n2*(n3*f4+f3)+f2)+f1;
                 *(g + ind_offset) = ints_1234[f1234];

               }
            }
         }
      }
   }

	return;
}

void construct_ao_g_wxyz_epsilon(double *g, int *s1, int *s2, int *s3, int *s4, double *epsilon, 
                                 int *thread, int *skip, int *n1, int *n2, int *n3, int *n4){

  electronic_repulsion_engines[*thread].set_precision(*epsilon);

  const auto& buf_vec = electronic_repulsion_engines[*thread].results(); // will point to computed shell sets

  electronic_repulsion_engines[*thread].compute(basis[*s1 - 1], basis[*s2 - 1], basis[*s3 - 1], basis[*s4 - 1]);

  auto ints_1234 = buf_vec[0]; // Location of computed integrals

   if (ints_1234 == nullptr)
   {
      *skip = 1;
   }
   else
   {
      *skip = 0;
      for(auto f1=0, f1234=0; f1!=*n1; ++f1){

         for(auto f2=0; f2!=*n2; ++f2){

            for(auto f3=0; f3!=*n3; ++f3){

               for(auto f4=0; f4!=*n4; ++f4, ++f1234){

                  int ind_offset = (*n1)*((*n2)*((*n3)*f4+f3)+f2)+f1;
                  *(g + ind_offset) = ints_1234[f1234];

               }
            }
         }
      }
   }

  return;
}
