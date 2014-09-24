#include<string>
#include<vector>

using namespace std;

#define N_CALC 200

// Model Parameters
#define N_GAUSS 1    /* Number of Gaussians for WD distribution */
#define WD_DIST 2    /* 1-Gaussians, 2-Uniform */
#define WD_M_MIN 0.3 /* Minimum WD mass for flat distribution */
#define WD_M_MAX 1.1 /* Maximum WD mass for flat distribution */
#define KMAX 25.0    /* Minimum detectable radial velocity */

// NS Gaussian Parameters
#define ADD_NS 1     /* 0 - No NS, 1 - Include NS */
#define NS_RATE 0.1  /* Percent of stars that are NS */
#define NS_M 1.35    /* NS mean mass adopted  */
#define NS_SD 0.02   /* NS standard deviation adopted  */

// Constants
#define PI acos(-1.0)
#define GGG 1.536e6  /* Newton's Gravitational Constant in (km/s)^3 days Msun^(-1) */


#define MSEED 161803397
#define MBIG 1.0e9
#define MZ 0
#define FAC (1.0/MBIG)
#define MAXIT 100

double func_non(double M2, double* model);
double prob_M2_wrapper(double M2_min, double* model);
double prob_M2_model(double* model);
double eval_M2(double M2, double* model);
double prob_NS(double frac_NS, double M1, double M2_min);
double prob_NS_all(double frac_NS, double M1);
double gauss_prob(double val, double w, double val_err, double val_in);

void create_gauss_dist(double* M,double* sd,double* w);
void find_M2(double M2, double* f, double* df, double* model);
double get_inc(long *seed);

double ran3(long *);
double gauss_ran(double val, double val_err, long* seed);
double rtsafe(void (*funcd)(double,double *,double *,double *),double x1,double x2,double xacc,double* model);

