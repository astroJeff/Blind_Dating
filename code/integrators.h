#include<vector>

using namespace std;

#define MSEED 161803397
#define MBIG 1.0e9
#define MZ 0
#define FAC (1.0/MBIG)
#define MAXIT 100

#define JMAX 20   /* Max number of iterations in integrators */

void trapzd(double (*func)(double, double *), double x1, double x2, double* y, int nsteps, double* model);
void qtrap(double (*func)(double, double *), double x1, double x2, double eps, double* y, double* model);
void qromb(double (*func)(double, double *), double x1, double x2, double eps, double* y, double* model);
void polint2(double* xa,double* ya,int n,int order,double x,double* y,double* dy);
void polint(double* xa,double* ya,int n,double x,double* y,double* dy);


double ran3(long *);
void sort(vector<double>& x,long n_objs);
double rtsafe(void (*funcd)(double,double *,double *,double *),double x1,double x2,double xacc,double* model);
void rkdumb(double* vstart, int nvar, double x1, double x2, int nstep, double* x_cur, double* model,
            void (*derivs)(double, double *, double *, double *));
void rk4(double *y, double *dydx, int n, double x, double h, double *yout, double* model,
         void (*derivs)(double, double *, double *, double *));
