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
 
  	Extract integrals
  	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
 
 	To extract integrals from Libint and place them into array in Fortran order.

*/
using namespace std;
#include "extract_integrals.h"

void extract_integrals(double *ints_fortran_order, const double *ints_cxx_order, int n1, int n2, double prefac){

   if (ints_cxx_order == nullptr){

      for(auto f1=0; f1!=n1; ++f1){
         for(auto f2=0; f2!=n2; ++f2){

            ints_fortran_order[n1*f2 + f1] = 0.0e0;

         }
      }
   }
   else{
      for(auto f1=0; f1!=n1; ++f1){
         for(auto f2=0; f2!=n2; ++f2){

            ints_fortran_order[n1*f2 + f1] = prefac*ints_cxx_order[f1*n2 + f2];

         }
      }
   }

}

