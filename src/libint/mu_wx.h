#ifdef __cplusplus
// Are we compiling this with a C++ compiler? Add extern "C" { ... }
extern "C" {
#else
#endif

void construct_ao_mu_wx(double *mu_X, double *mu_Y, double *mu_Z, int *s1, int *s2);

#ifdef __cplusplus
}
#endif
