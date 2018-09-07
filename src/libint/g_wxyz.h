#ifdef __cplusplus
// Are we compiling this with a C++ compiler? Add extern "C" { ... }
extern "C" {
#else
#endif

void construct_ao_g_wxyz_epsilon(double *g, long *s1, long *s2, long *s3, long *s4, double *epsilon, long *thread, long *skip,
							long *n1, long *n2, long *n3, long *n4);
void construct_ao_g_wxyz(double *g, long *s1, long *s2, long *s3, long *s4);

#ifdef __cplusplus
}
#endif
