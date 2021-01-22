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

#ifdef __cplusplus
// Are we compiling this with a C++ compiler? Add extern "C" { ... }
extern "C" {
#else
#endif

void get_eri(double *g, int *s1, int *s2, int *s3, int *s4, double *epsilon, 
             int *skip, int *n1, int *n2, int *n3, int *n4);

void get_eri_1der(double *g, int *s1, int *s2, int *s3, int *s4);

#ifdef __cplusplus
}
#endif
