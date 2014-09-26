#include<string>
#include<vector>

using namespace std;

#define N_CALC 100000
#define BURN_IN 500

#define M2_MODEL 1 /* 1-Gaussians  */
#define MODEL_FIX_GAUSS 0   /* 0-Allow Gaussians to move 1-Fix Gaussians  */
#define N_GAUSS 1
#define ADD_NS 1 /* Include a NS component */ 
#define USE_INC 1  /* Allow constraints from non-detection of eclipses */
                   /* Values taken from Hermes et al. (2014)  */
#define M2_GT_M1 1  /* Only allow solutions in which M2 > M1 */
                    /* This is because the observed WD must be lighter */

#define MOVE_MU 0.1
#define MOVE_SD 0.08
#define MOVE_W 0.01
#define MOVE_NS_FRAC 0.01

#define NS_MASS 1.35 /* NS Mass */

#define PI acos(-1.0)
#define GGG 1.536e6  /* Newton's Gravitational Constant in (km/s)^3 days Msun^(-1) */


void read_data(vector<string>& Names,vector<double>& M1,vector<double>& K,vector<double>& Porb,vector<double>& inc_max,vector<string>& inc_name);
void initial_guess(double* mu,double* sd,double* w,vector<double>& p_NS,double* frac_NS,vector<double>& inc_max,vector<string>& inc_name,vector<double>& M2,vector<double>& M2_min,vector<string>& Names,vector<double>& M1,vector<double>& Porb,vector<double>& K,long* seed);

void next_point(double mu[N_GAUSS],double sd[N_GAUSS],double w[N_GAUSS],double* frac_NS,double mu_new[N_GAUSS],double sd_new[N_GAUSS],double w_new[N_GAUSS],double* frac_NS_new,long* seed);
void move_to_point(double mu[N_GAUSS],double sd[N_GAUSS],double w[N_GAUSS],vector<double>& p_NS,double* frac_NS, double* P1,double mu_new[N_GAUSS],double sd_new[N_GAUSS],double w_new[N_GAUSS],vector<double>& p_NS_new,double* frac_NS_new,double* P2);

double P_model(double mu[N_GAUSS],double sd[N_GAUSS],double w[N_GAUSS],vector<double>& p_NS,double frac_NS,vector<double>& M2,vector<double>& M2_min,vector<string>& Names,vector<double>& M1,vector<double>& K,vector<double>& Porb,long* seed);

double prob_M2_model(double model[6+3*N_GAUSS]);
double prob_M2_wrapper(double M2_min, double* model);
double prob_NS(double frac_NS, double M1, double M2_min);
double prob_NS_all(double frac_NS, double M1);




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

