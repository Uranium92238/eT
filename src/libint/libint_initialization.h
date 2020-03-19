//
//
//  eT - a coupled cluster program
//  Copyright (C) 2016-2020 the authors of eT
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

// check if the correct shell ordering is set-up in libint
// however throw an error

// define the LIBINT_CGSHELL_ORDERING, which contains information about
// Gaussian shell ordering
#include "libint2/config.h"
#if LIBINT_CGSHELL_ORDERING != LIBINT_CGSHELL_ORDERING_STANDARD
#error "Libint compiled with an cg shell ordering different from cartesian"
#endif

#ifdef __cplusplus
// Are we compiling this with a C++ compiler? Add extern "C" { ... }
extern "C" {
#endif

void initialize_libint();

void finalize_libint();

void initialize_coulomb();

void export_geometry_and_basis_to_libint(const int nAtoms,
                                         const int *atomicNumbers,
                                         const double *atomicCoordinates,
                                         const char *basisSets,
                                         const int maxLength,
                                         const int *cartesians);

void initialize_kinetic();

void initialize_nuclear();

void initialize_overlap();

void initialize_dipole();

void initialize_quadrupole();

void set_coulomb_precision(double *prec);

void initialize_potential(double *charges, double *coordinates, int *n_points);

void initialize_shell2atom();

#ifdef __cplusplus
}
#endif
