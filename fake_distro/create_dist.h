#include<string>
#include<vector>

using namespace std;

#define N_CALC 100


#define PI acos(-1.0)
#define GGG 1.536e6  /* Newton's Gravitational Constant in (km/s)^3 days Msun^(-1) */


#define MSEED 161803397
#define MBIG 1.0e9
#define MZ 0
#define FAC (1.0/MBIG)
#define MAXIT 100

void find_M2(double M2, double* f, double* df, double* model);
double get_inc(long *seed);

double ran3(long *);
double gauss_ran(double val, double val_err, long* seed);
double rtsafe(void (*funcd)(double,double *,double *,double *),double x1,double x2,double xacc,double* model);
