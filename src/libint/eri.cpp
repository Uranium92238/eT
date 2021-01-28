//
//
//  eT - a coupled cluster program
//  Copyright (C) 2016-2021 the authors of eT
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
 
   Electron repulsion integral (eri) routines 
   Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2020
 
*/

using namespace std;

#include "eri.h"
#include "globals.h"
#include "omp_control.h"

void get_eri(double *g, 
             int *s1, int *s2, int *s3, int *s4, 
             double *epsilon, int *skip, 
             int *n1, int *n2, int *n3, int *n4){
/* 

    Get electron repulsion integrals (eri)
    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018-2020

    Calculates the electron repulsion integrals for a shell quartet (s1,s2,s3,s4) and stores the result in g.
    
    s1, s2, s3, s4  Shell quartet to calculate integrals for
    n1, n2, n3, n4  Number of AOs in the shells s1, s2, s3, s4
    epsilon         Libint precision value to use 
    skip            On exit, skip = 1 if Libint decided not to calculate any integrals; skip = 0 if 
                    Libint did calculate integrals. If skip = 1, the array g is not set to zero. If
                    you wish to 'calculate' the negligible integrals, you have to zero the elements of
                    g after the routine.

*/
  int thread = omp_get_thread_num();

  electronic_repulsion[thread].set_precision(*epsilon);

  const auto& buf_vec = electronic_repulsion[thread].results(); // will reference computed shell sets

  electronic_repulsion[thread].compute(basis[*s1 - 1], basis[*s2 - 1], basis[*s3 - 1], basis[*s4 - 1]);

  auto ints_1234 = buf_vec[0]; // Location of computed integrals

   if (ints_1234 == nullptr)
   {
      *skip = 1;
   }
   else
   {
      *skip = 0;

      for(auto f4=0; f4!=*n4; ++f4){
        for(auto f3=0; f3!=*n3; ++f3){
            for(auto f2=0; f2!=*n2; ++f2){
              for(auto f1=0; f1!=*n1; ++f1){

                  auto f1234 = (*n4)*((*n3)*((*n2)*f1+f2)+f3)+f4;
                  auto ind_offset = (*n1)*((*n2)*((*n3)*f4+f3)+f2)+f1;

                  *(g + ind_offset) = ints_1234[f1234];

               }
            }
         }
      }
   }

  return;
}

void get_eri_1der(double *g, int *s1, int *s2, int *s3, int *s4){
/* 

    Get electron repulsion integrals (eri)
    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018-2020

    Calculates the first derivative of the electron repulsion integrals 
    for a shell quartet (s1,s2,s3,s4) and stores the result in g.
    The ordering of g is wxyzqk, where w--z refers to AO indices, 
    q = 1,2,3 refers to x,y,z, and k = 1,2,3,4 refers to the centers associated
    with s1, s2, s3, and s4.
    
    s1, s2, s3, s4  Shell quartet to calculate integrals for

*/

  int thread = omp_get_thread_num();

  const auto& buf_vec = electronic_repulsion_1der[thread].results();

  std::size_t n1 = basis[*s1 - 1].size(); 
  std::size_t n2 = basis[*s2 - 1].size(); 
  std::size_t n3 = basis[*s3 - 1].size(); 
  std::size_t n4 = basis[*s4 - 1].size(); 

  electronic_repulsion_1der[thread].compute(basis[*s1 - 1], basis[*s2 - 1], 
                                    basis[*s3 - 1], basis[*s4 - 1]);

  auto offset = 0;

  for (auto k = 0, shell_set = 0; k != 4; ++k){ // Loop over shell centers/atoms
    for (auto q = 0; q != 3; ++q, ++shell_set){ // Loop over xyz on the given shell

        auto ints = buf_vec[shell_set];

        if (ints == nullptr)
        {
          for(std::size_t f1=0; f1!=n1; ++f1){

            for(std::size_t f2=0; f2!=n2; ++f2){

              for(std::size_t f3=0; f3!=n3; ++f3){

                for(std::size_t f4=0; f4!=n4; ++f4){

                  g[offset + n1*(n2*(n3*f4+f3)+f2)+f1] = 0.0e0;

                }
              }
            }  
          }
        }
        else
        {
          for(std::size_t f4=0; f4!=n4; ++f4){
            for(std::size_t f3=0; f3!=n3; ++f3){
              for(std::size_t f2=0; f2!=n2; ++f2){
                for(std::size_t f1=0; f1!=n1; ++f1){

                  auto f1234 = (n4)*((n3)*((n2)*f1+f2)+f3)+f4;
                  g[offset + n1*(n2*(n3*f4+f3)+f2)+f1] = ints[f1234];
  
                }
              }
            }
          }
        }

        offset = offset + n1*n2*n3*n4;

    }
  }
  return;
}
