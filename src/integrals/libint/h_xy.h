#ifdef __cplusplus
// Are we compiling this with a C++ compiler? Add extern "C" { ... }
extern "C" {
#else
#endif

void get_ao_xy_kinetic(double *h);
void get_n_aos(int *n_ao);

#ifdef __cplusplus
}
#endif
