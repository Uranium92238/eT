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

    Libint initialization of engines and basis set
    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018

*/
#include "libint_initialization.h"
#include "omp_control.h"

#include "eT_basis.h"



// Global variables declared in this file

// Basis set and atoms
eTBasis basis;
vector<Atom> atoms;

/*
   Vectors of Libint engines.
   Length of vectors is equal to number of threads
*/
vector<Engine> electronic_repulsion(omp_get_max_threads());
vector<Engine> electronic_repulsion_1der(omp_get_max_threads());
vector<Engine> kinetic(omp_get_max_threads());
vector<Engine> kinetic_1der(omp_get_max_threads());
vector<Engine> nuclear(omp_get_max_threads());
vector<Engine> nuclear_1der(omp_get_max_threads());
vector<Engine> overlap(omp_get_max_threads());
vector<Engine> overlap_1der(omp_get_max_threads());
vector<Engine> dipole(omp_get_max_threads());
vector<Engine> quadrupole(omp_get_max_threads());
vector<Engine> coulomb_external_charges(omp_get_max_threads());
vector<Engine> coulomb_external_unit_charges(omp_get_max_threads());


void export_geometry_and_basis_to_libint(const int nAtoms,
                                         const int *atomicNumbers,
                                         const int *atomicCharges,
                                         const double *atomicCoordinates,
                                         const char *basisSets,
                                         const int maxLength,
                                         const int *cartesians){

    /*
    Export geometry and basis set
    Written by Rolf H. Myhre, Mar. 2020
    Modified by Tor S. Haugland, May 2021. Add atomic charges.

    Exports the geometry and basis sets used from the Fortran side to Libint

    nAtoms is an int with number of atoms

    atomicNumbers is an int array of dimension n_atoms passed by reference

    atomicCharges is an int array of dimension n_atoms passed by reference. Usually
    equal to atomicNumbers, but zero for ghost atoms.

    atomicCordinates is a double array with dimension 3*n_atoms with atomic coordinates
    passed by reference

    basisSets is an array of c strings passed by reference

    maxLength is an integer indicating the storage length of the strings in basisSet

    cartesians is an array of ints passed by reference
    */

    vector<Atom> temporaryAtoms(nAtoms, Atom{0, 0.0, 0.0, 0.0});

    vector<string> basisList;
    vector<int> cartesiansList;
    vector<vector<int>> atomsWithBasis;

    bool newBasis;
    int nBasis;

    string basisString = string(&basisSets[0]);
    basisList.push_back(basisString);

    cartesiansList.push_back(cartesians[0]);
    atomsWithBasis.push_back({});

    nBasis = 1;
    int basisIndex = 0;

    for (int i = 0; i < nAtoms; i++)
    {
        basisIndex = nBasis - 1;

        temporaryAtoms[i].atomic_number = atomicNumbers[i];
        temporaryAtoms[i].x = atomicCoordinates[3*i];
        temporaryAtoms[i].y = atomicCoordinates[3*i+1];
        temporaryAtoms[i].z = atomicCoordinates[3*i+2];

        //Get current basis set
        basisString = string(&basisSets[i*maxLength]);

        newBasis = true;
        if(basisList[basisIndex] == basisString &&
           cartesiansList[basisIndex] == cartesians[i]){

            newBasis = false;
        }

        //If we found a new basis, add it to the list,
        //else just add the index to the list for the current basis
        if(newBasis){

            nBasis++;
            basisIndex = nBasis - 1;
            basisList.push_back(basisString);
            cartesiansList.push_back(cartesians[i]);
            atomsWithBasis.push_back({i});
        }
        else{
            atomsWithBasis[basisIndex].push_back(i);
        }
    }

    //Set global atoms variable
    atoms = temporaryAtoms;

    //Reset global basis variable
    eTBasis dummyBasis;
    basis = dummyBasis;

    //Temporary basis, one for each basis set in basisList
    BasisSet temporaryBasis;

    for(int i = 0; i < nBasis; i++){

        temporaryAtoms.clear();
        for (int const& atom: atomsWithBasis[i]){
            temporaryAtoms.push_back(atoms[atom]);
        }

        temporaryBasis = BasisSet(basisList[i], temporaryAtoms, true);

        //Set pure (spherical) or cartesian
        for(auto& shell: temporaryBasis) {
            for(auto& contraction: shell.contr) {

                if (cartesiansList[i] != 0) {
                   if (contraction.l >= 2 ) contraction.pure=false ;
                }
                else {
                   if (contraction.l >= 2 ) contraction.pure=true;
                }
            }
        }

        basis.add(temporaryBasis);
    }

    for (int i = 0; i < nAtoms; i++)
    {
        atoms[i].atomic_number = atomicCharges[i];
    }
}



void initialize_libint(){

    initialize();

}

void finalize_libint(){

    finalize();

}

void initialize_eri(const double eri_precision){

    electronic_repulsion[0] = Engine(Operator::coulomb, basis.max_nprim(), basis.max_l());
    electronic_repulsion_1der[0] = Engine(Operator::coulomb, basis.max_nprim(), basis.max_l(), 1);

    electronic_repulsion[0].set_precision(eri_precision);
    electronic_repulsion_1der[0].set_precision(eri_precision);

    for (int i = 1; i < omp_get_max_threads(); i++){

        electronic_repulsion[i] = electronic_repulsion[0];
        electronic_repulsion_1der[i] = electronic_repulsion_1der[0];

    }

}

void set_eri_precision(const double eri_precision){

    for (int i = 0; i < omp_get_max_threads(); i++){

        electronic_repulsion[i].set_precision(eri_precision);
        electronic_repulsion_1der[i].set_precision(eri_precision);

    }

}

void initialize_kinetic(){


    kinetic[0] = Engine(Operator::kinetic, basis.max_nprim(), basis.max_l());
    kinetic_1der[0] = Engine(Operator::kinetic, basis.max_nprim(), basis.max_l(),1);

    for (int i = 1; i < omp_get_max_threads(); i++){

        kinetic[i] = kinetic[0];
        kinetic_1der[i] = kinetic_1der[0];

    }

}

void initialize_nuclear(){

    vector<pair<double, array<double, 3>>> q;
    for (const auto& atom : atoms) {
      q.push_back({static_cast<double>(atom.atomic_number),
                   {{atom.x, atom.y, atom.z}}});
    }

    nuclear[0] = Engine(Operator::nuclear, basis.max_nprim(), basis.max_l());
    nuclear_1der[0] = Engine(Operator::nuclear, basis.max_nprim(), basis.max_l(), 1);

    nuclear[0].set_params(q);
    nuclear_1der[0].set_params(q);

    for (int i = 1; i < omp_get_max_threads(); i++){

        nuclear[i] = nuclear[0];
        nuclear_1der[i] = nuclear_1der[0];

    }

}

void initialize_overlap(){

    overlap[0] = Engine(Operator::overlap, basis.max_nprim(), basis.max_l());
    overlap_1der[0] = Engine(Operator::overlap, basis.max_nprim(), basis.max_l(),1);

    for (int i = 1; i < omp_get_max_threads(); i++){

        overlap[i] = overlap[0];
        overlap_1der[i] = overlap_1der[0];

    }

}

void initialize_dipole(){

    dipole[0] = Engine(Operator::emultipole1, basis.max_nprim(), basis.max_l()); // Note that expansion point = (0,0,0) by default

    for (int i = 1; i < omp_get_max_threads(); i++){
        dipole[i] = dipole[0];
    }

}

void initialize_quadrupole(){

    quadrupole[0] = Engine(Operator::emultipole2, basis.max_nprim(), basis.max_l()); // Note that expansion point = (0,0,0) by default

    for (int i = 0; i < omp_get_max_threads(); i++){
        quadrupole[i] = quadrupole[0];
    }

}

void initialize_coulomb_external_charges(const double *charges, const double *coordinates,
                                         const int n_points){

    coulomb_external_charges[0] = Engine(Operator::nuclear, basis.max_nprim(), basis.max_l());

    vector<pair<double, array<double, 3>>> points;

    for (int i = 0; i < n_points; i++) {

      int offset = (i-1)*3 + 3;
      points.push_back({(charges[i]),
                   {{*(coordinates+offset),*(coordinates+offset+1),
                    *(coordinates+offset+2)}}});

    }

    coulomb_external_charges[0].set_params(points);

    for (int i = 1; i < omp_get_max_threads(); i++){

        coulomb_external_charges[i] = coulomb_external_charges[0];

    }

}

void initialize_coulomb_external_unit_charges(const double *coordinates, const int n_points){

    coulomb_external_unit_charges[0] = Engine(Operator::nuclear, basis.max_nprim(), basis.max_l());

    vector<pair<double, array<double, 3>>> points;

    for (int i = 0; i < n_points; i++) {

      int offset = (i-1)*3 + 3;
      points.push_back({static_cast<double>(1.0e0),
                   {{*(coordinates+offset),*(coordinates+offset+1),
                    *(coordinates+offset+2)}}});
    }

    coulomb_external_unit_charges[0].set_params(points);

    for (int i = 1; i < omp_get_max_threads(); i++){

        coulomb_external_unit_charges[i] = coulomb_external_unit_charges[0];

    }

}
