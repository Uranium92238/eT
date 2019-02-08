/*
/
/ 	Construct dipole integrals mu
/ 	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
/
*/
#include <iostream>

using namespace std;

#include <fstream>
#include <string>
#include <vector>
#include "utils.h"
#include "mu_wx.h"

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

  // mu_X

  if (ints_shellset_X == nullptr) {

    for(auto f1=0; f1!=n1; ++f1){
      for(auto f2=0; f2!=n2; ++f2){

        *(mu_X + n1*f2 + f1) = 0.0e0;

      }
    }
  }
  else{

    for(auto f1=0; f1!=n1; ++f1){
      for(auto f2=0; f2!=n2; ++f2){

        *(mu_X + n1*f2 + f1) = ints_shellset_X[f1*n2 + f2]; 

      }
    }
  }

  // mu_Y

  if (ints_shellset_Y == nullptr) {

    for(auto f1=0; f1!=n1; ++f1){
      for(auto f2=0; f2!=n2; ++f2){

        *(mu_Y + n1*f2+f1) = 0.0e0;

      }
    }
  }
  else{

    for(auto f1=0; f1!=n1; ++f1){
      for(auto f2=0; f2!=n2; ++f2){

        *(mu_Y + n1*f2 + f1) = ints_shellset_Y[f1*n2 + f2]; 

      }
    }
  }

  // mu_Z

  if (ints_shellset_Z == nullptr) {

    for(auto f1=0; f1!=n1; ++f1){
      for(auto f2=0; f2!=n2; ++f2){

        *(mu_Z + n1*f2 + f1) = 0.0e0;

      }
    }
  }
  else{

    for(auto f1=0; f1!=n1; ++f1){
      for(auto f2=0; f2!=n2; ++f2){

        *(mu_Z + n1*f2 + f1) = ints_shellset_Z[f1*n2 + f2]; 

      }
    }
  }

return;

}
