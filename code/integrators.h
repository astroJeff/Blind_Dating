#define JMAX 20   /* Max number of iterations in integrators */

void trapzd(double (*func)(double, double *), double x1, double x2, double* y, int nsteps, double* model);
void qtrap(double (*func)(double, double *), double x1, double x2, double eps, double* y, double* model);
void qromb(double (*func)(double, double *), double x1, double x2, double eps, double* y, double* model);
void polint2(double* xa,double* ya,int n,int order,double x,double* y,double* dy);
void polint(double* xa,double* ya,int n,double x,double* y,double* dy);
