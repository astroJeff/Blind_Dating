#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<vector>
#include<cmath>
#include"non-detections.h"
#include"integrators.h"

using namespace std;

long N_accept=0;
long N_reject=0;

int main(int argc, char* argv[]){

  int j,index;
  long i,seed;

  double inc;
  double M1, Porb, M2, mu, sigma, frac_NS;
  double P_non;

  double model[5];
  double NS_model[2];

  double *w,*sd,*M0;


  
  seed = time(NULL);


  M0 = new double[N_GAUSS];
  sd = new double[N_GAUSS];
  w = new double[N_GAUSS];
  create_gauss_dist(M0,sd,w);
 
 
  // System parameters
  M1 = 0.25;
  Porb = 0.5;  // Assume a half-day orbit
  mu = 0.60;
  sigma = 0.30;
  frac_NS = 0.03;


  model[0] = M1;
  model[1] = Porb;
  model[2] = mu;   // mu
  model[3] = sigma;   // sigma
  model[4] = frac_NS;   // NS fraction


  P_non = 0.0;
  // Calculate integral
  qromb(func_non,0.0,0.4,1.0e-6,&P_non,model);

  // NS component
  inc = asin( pow(Porb/(2.0*PI*GGG)*(M1+NS_M)*(M1+NS_M),1.0/3.0) * KMAX / NS_M );
  P_non += (1.0 - cos(inc));

  cout << P_non << endl;
}


double func_non(double M2, double* model){

  double sini,inc;
  double Porb,M1;
  double mu, sigma, frac_NS;
  double prob;

  M1 = model[0];
  Porb = model[1];
  mu = model[2];
  sigma = model[3];
  frac_NS = model[4];


  sini = pow(Porb/(2.0*PI*GGG)*(M1+M2)*(M1+M2),1.0/3.0) * KMAX / M2;

  if(sini > 1.0) return gauss_prob(mu,1.0-frac_NS,sigma,M2);


  inc = asin(sini);

  prob = gauss_prob(mu,1.0-frac_NS,sigma,M2) * (1.0 - cos(inc));

  return prob;
}




double prob_M2_wrapper(double M2_min, double* model){
  double M1, P_M2;

  model[2] = M2_min;  
  M1 = model[0];

  P_M2 = prob_M2_model(model);

  return P_M2;
}


double prob_M2_model(double* model){
  double M2_min,M2_max,M2_prob;

  M2_prob = 0.0;

  if(WD_DIST == 1){
    M2_min = model[2];
    M2_max = 3.0;
  }else if(WD_DIST == 2){
    M2_min = WD_M_MIN;
    M2_max = WD_M_MAX;
  }


  // For all stars
  qromb(eval_M2,M2_min,M2_max,1.0e-5,&M2_prob,model);
  
  // For NS's
  //  qromb(eval_M2,NS_M-5*NS_SD,NS_M+5*NS_SD,1.0e-5,&M2_prob,model);

  return M2_prob;

}

double eval_M2(double M2, double* model){
  int i;
  
  double mu[N_GAUSS];
  double sd[N_GAUSS];
  double w[N_GAUSS];
  double M1,Porb,M2_min;
  double sini;
  double temp;



  M1 = model[0];     // Observed WD mass
  Porb = model[1];   // Orbital period
  M2_min = model[2];  


  for(i=0;i<N_GAUSS;i++){
    mu[i] = model[3*i+3];
    sd[i] = model[3*i+4];
    w[i] = model[3*i+5];
  }

  temp = 0.0;

  sini = (M2_min <= M2) ? pow((M1+M2)/(M1+M2_min) , 2.0/3.0) * M2_min/M2 : 0.0;

  if(WD_DIST == 1){

    for(i=0;i<N_GAUSS;i++) temp += gauss_prob(mu[i],w[i],sd[i],M2) * sini;
    
  }else if(WD_DIST == 2){

    if(WD_M_MIN < M2 && M2 < WD_M_MAX)  temp += sini;
    if(ADD_NS) temp += gauss_prob(NS_M,NS_RATE,NS_SD,M2) * sini;

  }


  return temp;
}




double prob_NS(double frac_NS, double M1, double M2_min){
  // P(NS) = sin(i)
  double M2 = NS_M;
  double prob;

  prob = frac_NS * pow((M1+M2)/(M1+M2_min), 2.0/3.0) * M2_min / M2;

  return (M2_min < M2) ? prob : 0.0;
}


double prob_NS_all(double frac_NS, double M1){
  // This is (should be) the analytic integral over all M2_min's for a given M1

  double xmin,xmax;
  double M2;
  //  double a1,a2;
  double a1,a2,a3,a4;

  xmin = 0.0;
  xmax = NS_M;  
  M2 = NS_M;

  a1 = 3.0 * xmax / M2 * pow((M1+M2)*(M1+M2)*(M1+xmax),1.0/3.0);
  a2 = -3.0 * xmin / M2 * pow((M1+M2)*(M1+M2)*(M1+xmin),1.0/3.0);

  a3 = -9.0/4.0 * pow(M1+M2,2.0/3.0) / M2 * pow(M1+xmax,4.0/3.0);
  a4 = 9.0/4.0 * pow(M1+M2,2.0/3.0) / M2 * pow(M1+xmin,4.0/3.0);
 
  
  //  a1 = -3.0/(4.0*M2) * (3.0*M1-xmax) * (M1+xmax) * pow((M1+M2)/(M1+xmax),2.0/3.0);
  //  a2 = -3.0/(4.0*M2) * (3.0*M1-xmin) * (M1+xmin) * pow((M1+M2)/(M1+xmin),2.0/3.0);
  
  //  return frac_NS*(a1-a2);
    return frac_NS*(a1+a2+a3+a4);
}



double gauss_prob(double val, double w, double val_err, double val_in){
  return (w/(val_err*sqrt(2.0*PI)) * exp(-(val-val_in)*(val-val_in)/(2.0*val_err*val_err)));
}



void create_gauss_dist(double* M,double* sd,double* w){
  int i;
  double weight;
  
  // Gaussian mixture model
  for(i=0;i<N_GAUSS;i++){
      
      if(i==0){
	M[i] = 0.7;
	w[i] = 0.9;
	sd[i] = 0.2;
      }
      
      if(i==1){
	M[i] = NS_M;
	w[i] = NS_RATE;
	sd[i] = NS_SD;
      }

      if(i==2){
	M[i] = 0.5;
	w[i] = 0.2;
	sd[i] = 0.1;
      }
      
    }
	
	// To normalize Gaussians
  weight = 0.0;
  for(i=0;i<N_GAUSS;i++) weight += w[i];
  for(i=0;i<N_GAUSS;i++) w[i] /= weight;
  
  if(ADD_NS) for(i=0;i<N_GAUSS;i++) w[i] *= (1.0-NS_RATE);
      
}
  


void find_M2(double M2, double* f, double* df, double* model){
  double M_tot,C;
  double M1, Mf, inc;
  
  M1 = model[0];
  Mf = model[1];
  inc = model[2];
  
  M_tot = M1+M2;
  C = sin(inc)*sin(inc)*sin(inc);
  
  *f = C*M2*M2*M2/(M_tot*M_tot) - Mf;
  *df = C*(3.0*M2*M2*M_tot*M_tot - 2.0*M2*M2*M2*M_tot) / (M_tot*M_tot*M_tot*M_tot);
}


double get_inc(long *seed){
  return acos(1.0 - ran3(seed));
}


double ran3(long *idum){ 
  static int inext, inextp;
  static long ma[56];
  static int iff=0;
  long mj, mk;
  int i, ii, k;
  
  if(*idum < 0 || iff == 0){
    iff=1;
    mj = abs(MSEED-abs(*idum));
    mj = mj % (long)MBIG;
    ma[55]=mj;
    mk=1;
    for(i=1;i<=54;i++){
      ii = (21*i) % 55;
      ma[ii] = mk;
      mk = mj - mk;
      if(mk < MZ) mk += (long)MBIG;
      mj = ma[ii];
    }
    for(k=1;k<=4;k++)
      for(i=1;i<=55;i++){
        ma[i] -= ma[1+(i+30) % 55];
        if(ma[i] < MZ) ma[i] += (long)MBIG;
      }
    inext=0;
    inextp=31;
    *idum=1;
  }
  if(++inext == 56) inext = 1;
  if(++inextp == 56) inextp = 1;
  mj=ma[inext]-ma[inextp];
  if(mj < MZ) mj += (long)MBIG;
  ma[inext]=mj;
  return mj*FAC;
}


double gauss_ran(double val, double val_err, long* seed){

  double ran_num;

  static int iset=0;
  static long myseed = 0;
  static float gset;

  double fac,rsq,v1,v2;

  if (*seed != -1)  myseed = *seed;

  if (myseed < 0) iset=0;
  if (iset == 0) {
    do {
      v1=2.0*ran3(&myseed)-1.0;
      v2=2.0*ran3(&myseed)-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;
    ran_num = v2*fac;
  } else {
    iset=0;
    ran_num = gset;
  }

  // Adjust the gaussian to the particular mu and sigma

  return ran_num*val_err+val;

}


double rtsafe(void (*funcd)(double,double *,double *,double *),double x1,double x2,double xacc,double* model){
  int j;
  double df,dx,dxold,f,fh,fl;
  double temp,xh,xl,rts;
  
  (*funcd)(x1,&fl,&df,model);
  (*funcd)(x2,&fh,&df,model);
  if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)){
//    cerr << "Root must be bracketed in rtsafe" << endl;
//    cerr << "M1 = " << M1 << " inc = " << inc << endl;
  }
  if (fl == 0.0) return x1;
  if (fh == 0.0) return x2;
  if (fl < 0.0) {
    xl=x1;
    xh=x2;
  } else {
    xh=x1;
    xl=x2;
  }
  rts=0.5*(x1+x2);
  dxold=fabs(x2-x1);
  dx=dxold;
  (*funcd)(rts,&f,&df,model);
  for (j=1;j<=MAXIT;j++) {
    if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0)
        || (fabs(2.0*f) > fabs(dxold*df))) {
      dxold=dx;
      dx=0.5*(xh-xl);
      rts=xl+dx;
      if (xl == rts) return rts;
    } else {
      dxold=dx;
      dx=f/df;
      temp=rts;
      rts -= dx;
      if (temp == rts) return rts;
    }
    if (fabs(dx) < xacc) return rts;
    (*funcd)(rts,&f,&df,model);
    if (f < 0.0)
      xl=rts;
    else
      xh=rts;
  }
  cerr << "Maximum number of iterations exceeded in rtsafe" << endl;
  return 0.0;
}

