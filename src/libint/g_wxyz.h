#ifdef __cplusplus
// Are we compiling this with a C++ compiler? Add extern "C" { ... }
extern "C" {
#else
#endif

void get_ao_g_wxyz_epsilon(double *g, long *s1, long *s2, long *s3, long *s4, double *epsilon);
void get_ao_g_wxyz(double *g, long *s1, long *s2, long *s3, long *s4);

#ifdef __cplusplus
}
#endif
