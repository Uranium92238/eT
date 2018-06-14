/*
/
/ 	Two-electron integral routines (for g_xyzw = g_αβγδ)
/ 	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
/
*/
#include "utils.h"
#include "g_xyzw.h"

#include <libint2.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace libint2;
using namespace std;

void get_ao_g_xyzw(double *g){

	initialize();

	string xyzfilename = "Water.xyz"; // see XYZ format description at http://en.wikipedia.org/wiki/XYZ_file_format
	ifstream input_file(xyzfilename);
	vector<Atom> atoms = read_dotxyz(input_file);

	BasisSet obs("cc-pVDZ", atoms);

	int num_aos = 0;

	for(auto s1=0; s1!=obs.size(); ++s1){

		num_aos = num_aos + obs[s1].size(); // Add number of basis functions in the given shell

	}

	Engine eri_engine(Operator::coulomb, obs.max_nprim(), obs.max_l());

	auto shell2bf = obs.shell2bf(); // maps shell index to basis function index
                                   // shell2bf[0] = index of the first basis function in shell 0
                                   // shell2bf[1] = index of the first basis function in shell 1
                                   // ...

	const auto& buf_vec = eri_engine.results(); // will point to computed shell sets

	for(auto s1=0; s1!=obs.size(); ++s1) {

		auto bf1 = shell2bf[s1];  // First basis function in shell 1
		auto n1 = obs[s1].size(); // Number of basis functions in shell 1

  		for(auto s2=0; s2!=obs.size(); ++s2) {

			auto bf2 = shell2bf[s2];  // First basis function in shell 2
			auto n2 = obs[s2].size(); // Number of basis functions in shell 2

			for (auto s3=0; s3!=obs.size(); ++s3) {

				auto bf3 = shell2bf[s3];  // First basis function in shell 3
				auto n3 = obs[s3].size(); // Number of basis functions in shell 3

				for (auto s4=0; s4!=obs.size(); ++s4) {

					auto bf4 = shell2bf[s4];  // First basis function in shell 4
					auto n4 = obs[s4].size(); // Number of basis functions in shell 4

					eri_engine.compute(obs[s1], obs[s2], obs[s3], obs[s4]);

					auto ints_1234 = buf_vec[0]; // Location of computed integrals
    				if (ints_1234 == nullptr)
      					continue;  // nullptr returned if all integrals were screened out

    				for(auto f1=0, f1234=0; f1!=n1; ++f1){

						auto bf1_ind = bf1 + f1;

      				for(auto f2=0; f2!=n2; ++f2){

							auto bf2_ind = bf2 + f2;

      					for(auto f3=0; f3!=n3; ++f3){

								auto bf3_ind = bf3 + f3;

      						for(auto f4=0; f4!=n4; ++f4, ++f1234){

									auto bf4_ind = bf4 + f4;

									*(g - 1 + index_four(bf1_ind+1,bf2_ind+1,bf3_ind+1,bf4_ind+1,num_aos,num_aos,num_aos)) = ints_1234[f1234];
								}
							}
						}
					}
				}
			}
  		}
	}

	finalize();

	return;
}
