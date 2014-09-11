#include<string>
#include<vector>

using namespace std;

#define N_CALC 100
#define BURN_IN 0
// #define M_DELTA 0.0002
// #define C_DELTA 0.0005
// #define SD_DELTA 0.0005
// #define SD_M2 0.00005


#define MSEED 161803397
#define MBIG 1.0e9
#define MZ 0
#define FAC (1.0/MBIG)
#define MAXIT 100

#define PI acos(-1.0)
#define GGG 1.536e6  /* Newton's Gravitational Constant in (km/s)^3 days Msun^(-1) */


void read_data(vector<string>& Names,vector<double>& M1,vector<double>& M1_err,vector<double>& K,vector<double>& K_err,vector<double>& Porb,vector<double>& Porb_err);
void initial_guess(double* mu,double* sd,double* w,vector<double>& inc,vector<double>& M2,vector<double>& M2_min,vector<string>& Names,vector<double>& M1,vector<double>& Porb,vector<double>& K,long* seed);
void mod_mix_gaus_gibbs(double mu[4],double sd[4],double w[4],vector<double>& inc,vector<double>& M2,vector<int>& C,vector<double>& M2_min,vector<string>& Names,vector<double>& M1,vector<double>& K,vector<double>& Porb,int print_flag,int l,ofstream* OUT,long* seed);

double prob_M2_model(double mu[4],double sd[4],double w[4],string Name,double M1,double M2_min,double K,double Porb,long* seed);
int prob_C_model(double M2,double mu[4],double sd[4],double w[4],long* seed);
void prob_mix_gauss(double mu[4],double sd[4],double w[4],vector<double>& M2,vector<int>& C);



double get_inc(long *seed);

double strtof(string temp);
int strtoi(string temp);
void split(string line, vector<string>& str_out);
double ran3(long *);
double gauss_ran(double val, double val_err, long* seed);
double gauss_prob(double val, double C, double val_err, double val_in);

void sort(vector<double>& x);
double rtsafe(void (*funcd)(double,double *,double *,double *),double x1,double x2,double xacc,double* model);
void rkdumb(double* vstart, int nvar, double x1, double x2, int nstep, double* x_cur, double* model,
            void (*derivs)(double, double *, double *, double *));
void rk4(double *y, double *dydx, int n, double x, double h, double *yout, double* model,
         void (*derivs)(double, double *, double *, double *));

void find_M2(double M2, double* f, double* df, double* model);
double eval_M2(double M2, double* model);
void integrate_M2(double M2, double* var, double* derivs, double* model);
void find_range(vector<double>& in_array,double* x_out,double* x_l,double* x_u);



// void print_results(double best_M[4],double best_C[4],double best_sd[4]);
// void iterate(double M0[4],double C[4],double sd[4],vector<double>& inc_m,vector<double>& M2_m,vector<double>& M2_min,vector<double>& M1,vector<double>& M1_err,vector<double>& M1_m,vector<double>& K,vector<double>& K_err,vector<double>& K_m,vector<double>& Porb,vector<double>& Porb_err,vector<double>& Porb_m,int print_flag,int k,ofstream* OUT,long* seed);
