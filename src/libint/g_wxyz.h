#ifdef __cplusplus
// Are we compiling this with a C++ compiler? Add extern "C" { ... }
extern "C" {
#else
#endif

void construct_ao_g_wxyz_epsilon(double *g, int *s1, int *s2, int *s3, int *s4, double *epsilon, 
                                 int *thread, int *skip, int *n1, int *n2, int *n3, int *n4);

void construct_ao_g_wxyz(double *g, int *s1, int *s2, int *s3, int *s4);

#ifdef __cplusplus
}
#endif
