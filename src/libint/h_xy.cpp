/*
/
/ 	One-electron integral routines (for h_xy = h_αβ)
/ 	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
/
*/

#include <iostream>

using namespace std;

#include <fstream>
#include <string>
#include <vector>

#include <libint2.hpp>

#include "utils.h"
#include "h_xy.h"

#include "globals.h"

using namespace libint2;

void get_n_aos(long *n_ao){

	//initialize();

	*n_ao = basis.nbf();
  //cout << *n_ao << endl;

	//finalize();

	return;

}

void get_ao_h_xy(double *h){
//
// Initialize libint calculator
//
	int thread = omp_get_thread_num();
//
// Calculate number of basis functions from the basis object
//
	int num_aos = 0;
	num_aos = basis.nbf();
//
// ** Compute the kinetic energy part of one-electron integrals **
//
	auto shell2bf = basis.shell2bf(); // maps shell index to basis function index
                                   // shell2bf[0] = index of the first basis function in shell 0
                                   // shell2bf[1] = index of the first basis function in shell 1
                                   // ...
//
	const auto& buf_vec = kinetic[thread].results(); // will point to computed shell sets
//
	for(auto s1=0; s1!=basis.nshells(); ++s1) {
  		for(auto s2=0; s2!=basis.nshells(); ++s2) {
//
    		kinetic[thread].compute(basis[s1], basis[s2]);
    		auto ints_shellset = buf_vec[0];  // location of the computed integrals
    		if (ints_shellset == nullptr) {
      		continue;  // nullptr returned if the entire shell-set was screened out
        }

    		auto bf1 = shell2bf[s1];  // first basis function in first shell
    		auto n1 = basis[s1].size(); // number of basis functions in first shell
    		auto bf2 = shell2bf[s2];  // first basis function in second shell
    		auto n2 = basis[s2].size(); // number of basis functions in second shell

    		// integrals are packed into ints_shellset in row-major (C) form
    		// this iterates over integrals in this order
    		for(auto f1=0; f1!=n1; ++f1){
      		for(auto f2=0; f2!=n2; ++f2){
					*(h - 1 + index_two(bf1+1+f1, bf2+1+f2, num_aos)) = ints_shellset[f1*n2+f2];
				}
			}
  		}
	}
//
//	Compute the nuclear attraction energy part of one-electron integrals
//
// Make point charges (R_I, Z_I)
//
	nuclear[thread].set_params(make_point_charges(atoms));  // convert `atoms` to point charges
//
	const auto& buf_vec_n = nuclear[thread].results(); // will point to computed shell sets
                                          	  // const auto& is very important!
//
	for(auto s1=0; s1!=basis.nshells(); ++s1) {
  		for(auto s2=0; s2!=basis.nshells(); ++s2) {

    		nuclear[thread].compute(basis[s1], basis[s2]);
    		auto ints_shellset_n = buf_vec_n[0];  // location of the computed integrals
    		if (ints_shellset_n == nullptr)
      		continue;  // nullptr returned if the entire shell-set was screened out

    		auto bf1 = shell2bf[s1];  // first basis function in first shell
    		auto n1 = basis[s1].size(); // number of basis functions in first shell
    		auto bf2 = shell2bf[s2];  // first basis function in second shell
    		auto n2 = basis[s2].size(); // number of basis functions in second shell

    		// integrals are packed into ints_shellset in row-major (C) form
    		// this iterates over integrals in this order
    		for(auto f1=0; f1!=n1; ++f1){
      		for(auto f2=0; f2!=n2; ++f2){
					*(h - 1 + index_two(bf1+1+f1, bf2+1+f2, num_aos)) = *(h - 1 + index_two(bf1+1+f1, bf2+1+f2, num_aos)) + ints_shellset_n[f1*n2+f2];
				}
			}
//
  		}
	}
//
	// finalize();
//
	return;
}

void get_ao_h_xy_sp(double *h, long *s1, long *s2){
//
// Initialize libint calculator
//
// Calculate number of basis functions from the basis object
  int thread = omp_get_thread_num();
//
// ** Compute the kinetic energy part of one-electron integrals **
// 
  const auto& buf_vec = kinetic[thread].results(); // will point to computed shell sets
//
  kinetic[thread].compute(basis[*s1 - 1], basis[*s2 - 1]);
//
  auto ints_shellset = buf_vec[0];  // location of the computed integrals
//
  auto n1 = basis[*s1 - 1].size(); // Number of basis functions in shell 1
  auto n2 = basis[*s2 - 1].size(); // number of basis functions in shell 2

  if (ints_shellset == nullptr) {
    // integrals are packed into ints_shellset in row-major (C) form
    // this iterates over integrals in this order
    for(auto f1=0; f1!=n1; ++f1){
      for(auto f2=0; f2!=n2; ++f2){

      *(h + n1*f2+f1) = 0.0e0;

      }
    }
  }
  else{
    // integrals are packed into ints_shellset in row-major (C) form
    // this iterates over integrals in this order
    for(auto f1=0; f1!=n1; ++f1){
      for(auto f2=0; f2!=n2; ++f2){

      *(h + n1*f2+f1) = ints_shellset[f1*n2+f2];  //index_two(bf1+1+f1, bf2+1+f2, num_aos))

      }
    }
  }
//
//  Compute the nuclear attraction energy part of one-electron integrals
//
// Make point charges (R_I, Z_I)
//
  nuclear[thread].set_params(make_point_charges(atoms));  // convert `atoms` to point charges
//
  const auto& buf_vec_n = nuclear[thread].results(); // will point to computed shell sets
                                              // const auto& is very important!
//
  nuclear[thread].compute(basis[*s1 - 1], basis[*s2 - 1]);
  auto ints_shellset_n = buf_vec_n[0];  // location of the computed integrals
  if (ints_shellset_n != nullptr){
    // integrals are packed into ints_shellset in row-major (C) form
    // this iterates over integrals in this order
    for(auto f1=0; f1!=n1; ++f1){
      for(auto f2=0; f2!=n2; ++f2){

      *(h + n1*f2+f1) = *(h + n1*f2+f1) + ints_shellset_n[f1*n2+f2];

      }
    }
  }
//
//
  return;
}
