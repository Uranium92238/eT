//
//
//  eT - a coupled cluster program
//  Copyright (C) 2016-2022 the authors of eT
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

#ifdef __cplusplus
// Are we compiling this with a C++ compiler? Add extern "C" { ... }
extern "C" {
#else
#endif

void get_eri(double *g, 
             const int s1, const int s2, const int s3, const int s4, 
             const double epsilon_, int *skip, 
             const int n1, const int n2, const int n3, const int n4);

void get_eri_2c(double *g,
             const int J, const int K,
             const double epsilon_, int *skip,
             const int nJ, const int nK);

void get_eri_3c(double *g,
             const int J, const int s3, const int s4,
             const double epsilon_, int *skip,
             const int nJ, const int n3, const int n4);


void get_eri_1der(double *g, const int s1, const int s2, const int s3, const int s4);

#ifdef __cplusplus
}
#endif
