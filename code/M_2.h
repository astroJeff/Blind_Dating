#include<string>
#include<vector>

using namespace std;

#define N_CALC 10000
#define BURN_IN 50

#define M2_MODEL 1 /* 1-Gaussians  */
#define MODEL_FIX_GAUSS 0   /* 0-Allow Gaussians to move 1-Fix Gaussians  */
#define N_GAUSS 1



#define PI acos(-1.0)
#define GGG 1.536e6  /* Newton's Gravitational Constant in (km/s)^3 days Msun^(-1) */


void read_data(vector<string>& Names,vector<double>& M1,vector<double>& M1_err,vector<double>& M2_min,vector<double>& K,vector<double>& K_err,vector<double>& Porb,vector<double>& Porb_err);
void initial_guess(double* mu,double* sd,double* w,vector<double>& inc,vector<double>& M2,vector<double>& M2_min,vector<string>& Names,vector<double>& M1,vector<double>& Porb,vector<double>& K,long* seed);
void mod_mix_gaus_gibbs(double mu[4],double sd[4],double w[4],vector<double>& inc,vector<double>& M2,vector<int>& C,vector<double>& M2_min,vector<string>& Names,vector<double>& M1,vector<double>& K,vector<double>& Porb,int print_flag,int l,ofstream* OUT,long* seed);
void iterate(double mu[4],double sd[4],double w[4],vector<double>& inc,vector<double>& M2,vector<int>& C,vector<double>& M2_min,vector<string>& Names,vector<double>& M1,vector<double>& K,vector<double>& Porb,int print_flag,int l,ofstream* OUT,long* seed);

double prob_M2_wrapper(double M2_min, double* model);
double prob_M2_model(double model[9]);
int prob_C_model(double M2,double mu[4],double sd[4],double w[4],long* seed);
void prob_mix_gauss(double mu[4],double sd[4],double w[4],vector<double>& M2,vector<int>& C,long* seed);

void find_gauss(double* mean,double* sd,vector<double> vals);



double get_inc(long *seed);

double strtof(string temp);
int strtoi(string temp);
void split(string line, vector<string>& str_out);

double gauss_ran(double val, double val_err, long* seed);
double gauss_prob(double val, double C, double val_err, double val_in);
double poisson_prob(double lambda, double n);
double poisson_ran(double lambda,long* seed);
double factorial(double n);

void find_M2(double M2, double* f, double* df, double* model);
double eval_M2(double M2, double* model);
void integrate_M2(double M2, double* f, double* df, double* model);
void find_range(vector<double>& in_array,double* x_out,double* x_l,double* x_u);

