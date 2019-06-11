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

#ifdef __cplusplus
// Are we compiling this with a C++ compiler? Add extern "C" { ... }
extern "C" {
#else
#endif

void construct_ao_h_wx(double *h, int *s1, int *s2);

void construct_ao_h_wx_1der(double *h_1x, double *h_1y, double *h_1z, 
                  double *h_2x, double *h_2y, double *h_2z, int *s1, int *s2);

void construct_ao_h_wx_kinetic_1der(double *h_1x, double *h_1y, double *h_1z, 
               double *h_2x, double *h_2y, double *h_2z, int *s1, int *s2);

void construct_and_add_ao_h_wx_nuclear_1der(double *h_wxqk, int *s1, int *s2, int *n_ao);

#ifdef __cplusplus
}
#endif
