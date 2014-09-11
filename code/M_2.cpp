#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<vector>
#include<cmath>
#include<iomanip>
#include"M_2.h"
#include"integrators.h"

using namespace std;


int main(int argc, char* argv[]){

  int k;
  int print_flag,n_objs;

  long i,seed;


  double mu[N_GAUSS];
  double sd[N_GAUSS];
  double w[N_GAUSS];

  ofstream OUT;
  
  // Observations
  vector<string> Names;
  vector<double> M1;
  vector<double> M1_err;
  vector<double> K;
  vector<double> K_err;
  vector<double> Porb;
  vector<double> Porb_err;
  
  


  
  seed = time(NULL);
  
  OUT.open("dist_params.dat");


  // Take in data
  read_data(Names,M1,M1_err,K,K_err,Porb,Porb_err);

  
  // Declare and initialize model variables
  n_objs = M1.size();


  vector<double> M2_min(n_objs,0.0);
  vector<double> M2(n_objs,0.0);
  vector<double> inc(n_objs,0.0);
  vector<int> C(n_objs,0);


  // Give initial values
  initial_guess(mu,sd,w,inc,M2,M2_min,Names,M1,Porb,K,&seed);



  // k number of walkers
  for(k=0;k<1;k++){



    // Iterate for BURN_IN times for "burn in"
    print_flag = 0;  // Don't print out for "burn in"
    
    for(i=0;i<BURN_IN;i++) mod_mix_gaus_gibbs(mu,sd,w,inc,M2,C,M2_min,Names,M1,K,Porb,print_flag,k,&OUT,&seed);



    print_flag = 1;  // Print out for iterations

    // Iterate for N_CALC times
    for(i=0;i<N_CALC;i++){

      // MCMC algorithm in here    
      mod_mix_gaus_gibbs(mu,sd,w,inc,M2,C,M2_min,Names,M1,K,Porb,print_flag,k,&OUT,&seed);
    
    }
  }


  // Print results
//  print_results(best_M,best_C,best_sd);

  


  OUT.close();

}


void initial_guess(double* mu,double* sd,double* w,vector<double>& inc,vector<double>& M2,vector<double>& M2_min,vector<string>& Names,vector<double>& M1,vector<double>& Porb,vector<double>& K,long* seed){
  int i;
  
  double Mf;
  double model[3];

  for(i=0;i<N_GAUSS;i++){
    if(i==0){
      mu[0] = 0.9;    // High mass WD
      sd[0] = 0.4;
      w[0] = 0.25;
    }

    if(i==1){
      mu[1] = 0.6;    // Standard WD
      sd[1] = 0.4;
      w[1] = 0.25;
    }

    if(i==2){
      mu[2] = 0.2;    // ELM WD
      sd[2] = 0.4;
      w[2] = 0.25;
    }

    if(i==3){
      mu[3] = 1.4;    // NS
      sd[3] = 0.4;
      w[3] = 0.25;
    }
  
  }

  
  for(i=0;i<M2.size();i++){

    inc[i] = get_inc(seed);

    Mf = Porb[i]*K[i]*K[i]*K[i]/(2.0*PI*GGG);

    model[0] = M1[i];
    model[1] = Mf;
    model[2] = PI/2.0;
    M2_min[i] = rtsafe(find_M2,0.0,10.0,1.0e-5,model);

    model[2] = inc[i];
    M2[i] = rtsafe(find_M2,M2_min[i],10.0,1.0e-5,model);

    if(Names[i] == "J0106-1000") M2[i] = 0.43;
    if(Names[i] == "J0651+2844") M2[i] = 0.50;
    if(Names[i] == "NLTT11748") M2[i] = 0.76;    

  }

  
}









void mod_mix_gaus_gibbs(double mu[N_GAUSS],double sd[N_GAUSS],double w[N_GAUSS],vector<double>& inc,vector<double>& M2,vector<int>& C,vector<double>& M2_min,vector<string>& Names,vector<double>& M1,vector<double>& K,vector<double>& Porb,int print_flag,int l,ofstream* OUT,long* seed){
  int i,j,k;
  int n_objs,nsteps;


  // Number of objects
  n_objs = K.size();
  

  // First, find P( M2, C | mu, sd, w )
  for(j=0;j<n_objs;j++){


    //  P( M2[j] | mu, sd, w )
    M2[j] = prob_M2_model(mu,sd,w,Names[j],M1[j],M2_min[j],K[j],Porb[j],seed);


    //  P( C[j] | M2[j], mu, sd, w )  
    C[j] = prob_C_model(M2[j],mu,sd,w,seed);
  
  }



// Second, find P( mu, sd, w | M2, C )
  prob_mix_gauss(mu,sd,w,M2,C);  
  
  
  
  if(print_flag){
    *OUT << setprecision(4) << l << " ";
    for(i=0;i<N_GAUSS;i++)  *OUT << mu[i] << " " << sd[i] << " " << w[i] << " ";    
    *OUT << endl;
  }

}





double prob_M2_model(double mu[N_GAUSS],double sd[N_GAUSS],double w[N_GAUSS],string Name,double M1,double M2_min,double K,double Porb,long* seed){
  //  P( M2[j] | mu, sd, w, M2_min[j] )

  int i,j;
  int nsteps;
  
  double M2;
  double A;
  double model[18];
  double M2_prob;
  

  model[0] = M1;
  model[1] = K;
  model[2] = Porb;
  model[3] = 1.0;   // This is where "A" goes
  model[4] = M2_min;
  model[5] = 0.0;

  if(M2_MODEL==1){
    // Copy variables to array for rkdumb
    for(i=0;i<N_GAUSS;i++){
      model[3*i+6] = mu[i];
      model[3*i+7] = sd[i];
      model[3*i+8] = w[i];
    }
  }


  if(Name == "J0106-1000") return 0.43;
  if(Name == "J0651+2844") return 0.50;
  if(Name == "NLTT11748") return 0.76;



// First, normalize probability -> Find A

  // Since non-zero probability at sin(i) = 1
  M2_prob = 0.0;
//  for(i=0;i<4;i++) M2_prob += gauss_prob(mu[i],w[i],sd[i],M2_min);
  
  // Integrate over full range [M2_min,5.0] to get total probability
  qromb(eval_M2,M2_min,5.0,1.0e-8,&M2_prob,model);
  
  // Normalize to 1
  A = 1.0/M2_prob;
  model[3] = A;

    
    
    
// Second, draw random number, then integrate to that point

  model[5] = ran3(seed);

  // Integrate over range [M2_min,M2], until we reach the target integrated area
  // From Maria's suggestion, use a root finder on a function that calls the integrator
  // Thanks for the idea, Maria!
  

  M2 = rtsafe(integrate_M2,M2_min,5.0,1.0e-8,model);
  
  


  return M2;


}


/*
void mix_gauss_mass(double M2, double* f, double* df, double* model){
  int i;

  double mu[N_GAUSS];
  double sd[N_GAUSS];
  double w[N_GAUSS];
  double area;
  double M2_min;


  area = model[0];

  if(M2_MODEL == 1){
    for(i=0;i<N_GAUSS;i++){
      mu[i] = model[3*i+1];
      sd[i] = model[3*i+2];
      w[i] = model[3*i+3];
    }
  }


  *f = -2.0*area;
  for(i=0;i<4;i++) *f += w[i] + w[i]*erf( (M2-mu[i])/(sd[i]*sqrt(2.0)) );

  
  *df = 0.0;
  for(i=0;i<4;i++) *df += w[i] * 2.0/sqrt(PI) * exp( -(M2-mu[i])*(M2-mu[i])/(2.0*sd[i]*sd[i]) );

}
*/


int prob_C_model(double M2,double mu[N_GAUSS],double sd[N_GAUSS],double w[N_GAUSS],long* seed){
  // Don't have to worry about sin(i) here - Affects all models at same M2 the same

  int i;

  double temp;
  double prob_temp[N_GAUSS];
  double prob_tot;
  double prob_ran;
  
  prob_tot = 0.0;
  for(i=0;i<N_GAUSS;i++){
    prob_temp[i] = gauss_prob(mu[i],w[i],sd[i],M2);
    prob_tot += prob_temp[i];
  }
  for(i=0;i<N_GAUSS;i++) prob_temp[i] = prob_temp[i]/prob_tot;

  prob_ran = ran3(seed);

  temp = 0.0;
  for(i=0;i<N_GAUSS;i++){
    temp += prob_temp[i];
    if(prob_ran < temp) return i;
  }
  
}


void prob_mix_gauss(double mu[N_GAUSS],double sd[N_GAUSS],double w[N_GAUSS],vector<double>& M2,vector<int>& C){  
  // To determine P( mu, sigma, w | M2, C )

  int i,j;
  int n_objs,n_C;
  
  double mean;

  n_objs = M2.size();

  for(i=0;i<N_GAUSS;i++){

    n_C = 0;
    mean = 0.0;
    for(j=0;j<n_objs;j++){

      if(C[j] == i){
        n_C ++;
    
        mean += M2[j];    
    
      }

    }  
  
    // If there are one or no data points in Gaussian keep same values as before
    if( n_C < 2 ){
      w[i] = 1.0/(double) n_objs;
    } else {
      // Weights are the number of data points in each Gaussian
      w[i] = (double)n_C / (double) n_objs;


      // Means are unweighted averages.
      // Can keep gaussians from moving in model
      // Change MODEL_FIX_GAUSS in M_2.h
      if(!MODEL_FIX_GAUSS) mu[i] = mean / (double) n_C;


      // Standard deviations are the calculated standard deviations
      sd[i] = 0.0;
      for(j=0;j<n_objs;j++){
        if(C[j] == i) sd[i] += (mu[i] - M2[j]) * (mu[i] - M2[j]);
      } 
      sd[i] = sqrt( 1.0 / (double(n_C) - 1.0) ) * sqrt(sd[i]);

    }    
  
  }


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


double eval_M2(double M2, double* model){
  int i;
  
  double mu[N_GAUSS];
  double sd[N_GAUSS];
  double w[N_GAUSS];
  double M1,K,Porb,A;
  double sini;
  double temp;


  M1 = model[0];     // Observed WD mass
  K = model[1];      // Line of sight orbital velocity
  Porb = model[2];   // Orbital period
  A = model[3];      // Normalization constant

  for(i=0;i<N_GAUSS;i++){
    mu[i] = model[3*i+6];
    sd[i] = model[3*i+7];
    w[i] = model[3*i+8];
  }

  temp = 0.0;
  for(i=0;i<N_GAUSS;i++){  

    sini = pow(Porb/(2.0*PI*GGG) * (M1+M2)*(M1+M2), 1.0/3.0) * K / M2;

    temp += A * gauss_prob(mu[i],w[i],sd[i],M2) * sini;

  }

  return temp;
}

void integrate_M2(double M2, double* f, double* df, double* model){
  int i;
  
  double mu[N_GAUSS];
  double sd[N_GAUSS];
  double w[N_GAUSS];
  
  double M1,K,Porb,A,M2_min;
  double area;


  M1 = model[0];     // Observed WD mass
  K = model[1];      // Line of sight orbital velocity
  Porb = model[2];   // Orbital period
  A = model[3];      // Normalization constant
  M2_min = model[4]; // Minimum M2
  area = model[5];   // Value to integrate to

  for(i=0;i<N_GAUSS;i++){
    mu[i] = model[3*i+6];
    sd[i] = model[3*i+7];
    w[i] = model[3*i+8];
  }


  // Integrate [M2_min,M2] to get area, update (*f)
  if(abs(M2-M2_min) < 1.0e-10){
    *f = 0.0;
  } else {
    qromb(eval_M2,M2_min,M2,1.0e-8,f,model);
  }
  
  // Subtract area we are trying to get to
  *f = (*f) - area;
  
  // Derivative of integral is just the function
  *df = eval_M2(*f,model);

}




void read_data(vector<string>& Names,vector<double>& M1,vector<double>& M1_err,vector<double>& K,vector<double>& K_err,vector<double>& Porb,vector<double>& Porb_err){
  fstream IN;
  
  double M1_temp, M1_err_temp;
  double K_temp, K_err_temp;
  double P_temp, P_err_temp;
  
  string Name_temp,line;
  
  vector<string> data;
  



  // Open data file
  IN.open("ELM_WD.dat");
  //  IN.open("dist_1_gauss.dat");


  // Read in data, add M2,min to vector M2
  while(getline(IN,line)){

    split(line,data);

    
    Name_temp = data[0];
    P_temp = strtof(data[1]);
    P_err_temp = strtof(data[2]);
    K_temp = strtof(data[3]);
    K_err_temp = strtof(data[4]);
    M1_temp = strtof(data[7]);
    M1_err_temp = 0.02;  

    // Don't worry about non-detections yet
    if (P_temp > 0.001){
      Names.push_back(Name_temp);
      Porb.push_back(P_temp);
      Porb_err.push_back(P_err_temp);
      K.push_back(K_temp);
      K_err.push_back(K_err_temp);      
      M1.push_back(M1_temp);
      M1_err.push_back(M1_err_temp);
    }
    

    /*
    Names.push_back(data[3]);
    Porb.push_back(strtof(data[1]));
    M1.push_back(strtof(data[0]));
    K.push_back(strtof(data[2]));
    */

  }
  
  IN.close();
}

double get_inc(long *seed){
  return acos(1.0 - ran3(seed));
}


double strtof(string temp){
  double out;
  istringstream str_in(temp);
  str_in >> out;
  return out;
}

int strtoi(string temp){
  int out;
  istringstream str_in(temp);
  str_in >> out;
  return out;
}

void split(string line, vector<string>& str_out){
  string tmp;
  size_t spot;
  str_out.clear();
  while(!line.empty()){
    spot = line.find(" ");
    if (spot>1000) spot = line.size();  // For the last element
    tmp = line.substr(0,spot);
    str_out.push_back(tmp);
    tmp.clear();
    line.erase(0,spot+1);
  }
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



double gauss_prob(double val, double w, double val_err, double val_in){
  return (w/(val_err*sqrt(2.0*PI)) * exp(-(val-val_in)*(val-val_in)/(2.0*val_err*val_err)));
}




void sort(vector<double>& x){
  // Heapsort algorithm

  long e,i,j,k,l;

  double tmp;

  e = N_CALC-1;

  // Build heap
  i = (int)(e/2);
  while(i>=0){
    j = i;
    while(j*2+1<=e){
      k = j*2+1;
      l = j;
      if(x[l] < x[k]) l = k;
      if(k+1 <= e && x[l] < x[k+1]) l = k+1;
      if(l!=j){
        tmp = x[j];
        x[j] = x[l];
        x[l] = tmp;
        j = l;
      } else break;
    }
    i--;
  }

  // Sort heap
  while(e>0){
    tmp = x[e];
    x[e] = x[0];
    x[0] = tmp;
    e--;
    
    j = 0;
    while(j*2+1<=e){
      k = j*2+1;
      l = j;
      if(x[l] < x[k]) l = k;
      if(k+1 <= e && x[l] < x[k+1]) l = k+1;
      if(l != j){
        tmp = x[j];
        x[j] = x[l];
        x[l] = tmp;
        j = l;
      }else break;
    }
  }
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






void rkdumb(double *vstart, int nvar, double x1, double x2, int nstep, double* x_cur,double* model,
            void (*derivs)(double, double *, double *, double *))
{
  int i,k;
  double x,h;
  double *v,*vout,*dv;

  v = new double[nvar];
  vout = new double[nvar];
  dv = new double[nvar];

  for (i=0;i<nvar;i++){
    v[i]=vstart[i];
  }
  x=x1;
  h=(x2-x1)/nstep;

  for (k=0;k<nstep;k++){
    (*derivs)(x,v,dv,model);
    
//    cout << x << " " << v[0] << " " << dv[0] << " ";
    
    rk4(v,dv,nvar,x,h,vout,model,derivs);
//    cout << v[0] << " " << dv[0] << endl;
    
    if((double)(x+h) == x) cerr << "Step size too small in routine rkdumb" << endl;
    x += h;
    (*x_cur) = x;
    for (i=0;i<nvar;i++){
      v[i]=vout[i];
    }
  }

  if (k==nstep)  for (i=0;i<nvar;i++) vstart[i]=v[i];
      
  delete [] dv;
  delete [] vout;
  delete [] v;
}


void rk4(double *y, double *dydx, int n, double x, double h, double *yout, double* model,
         void (*derivs)(double, double *, double *, double *))
{
  int i;
  double xh,hh,h6,*dym,*dyt,*yt;

  dym = new double[n];
  dyt = new double[n];
  yt = new double[n];

  hh=h*0.5;
  h6=h/6.0;
  xh=x+hh;
  for (i=0;i<n;i++) yt[i]=y[i]+hh*dydx[i];
  (*derivs)(xh,yt,dyt,model);
  for(i=0;i<n;i++) yt[i]=y[i]+hh*dyt[i];
  (*derivs)(xh,yt,dym,model);
  for(i=0;i<n;i++){
    yt[i]=y[i]+h*dym[i];
    dym[i] += dyt[i];
  }
  (*derivs)(x+h,yt,dyt,model);
  for(i=0;i<n;i++)
    yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
    
  delete [] yt;
  delete [] dyt;
  delete [] dym;
 }


void find_range(vector<double>& in_array,double* x_out,double* x_l,double* x_u){

  long i;

  sort(in_array);

  // Get the value one sigma away from the center
  *x_l = in_array[(long)(N_CALC*0.159)];
  *x_u = in_array[(long)(N_CALC*0.841)];
  *x_out = in_array[(long)(N_CALC*0.5)];
}







/*



void iterate(double M0[4],double C[4],double sd[4],vector<double>& inc_m,vector<double>& M2_m,vector<double>& M2_min,vector<double>& M1,vector<double>& M1_err,vector<double>& M1_m,vector<double>& K,vector<double>& K_err,vector<double>& K_m,vector<double>& Porb,vector<double>& Porb_err,vector<double>& Porb_m,int print_flag,int l,ofstream* OUT,long* seed){
  int i,j,k;
  int n_objs;

  double p1,p2,prob_part;
  double Mf_temp,M2_temp;
  double ran_val;
  double M0_new[4],C_new[4],sd_new[4];
  double C_tot;
  double sd_inc,sd_M2;
  
  // Number of objects
  n_objs = K.size();
  
  vector<double> inc_new(n_objs,0.0);
  vector<double> M1_new(n_objs,0.0);
  vector<double> M2_new(n_objs,0.0);
  vector<double> Porb_new(n_objs,0.0);
  vector<double> K_new(n_objs,0.0);




  // Determine probability for current point


  p1 = 0.0;
  for(i=0;i<n_objs;i++){
    prob_part = 0.0;

//    prob_part += log10(gauss_prob(M1[i],1.0,M1_err[i],M1_m[i]));
//    prob_part += log10(gauss_prob(Porb[i],1.0,Porb_err[i],Porb_m[i]));
//    prob_part += log10(gauss_prob(K[i],1.0,K_err[i],K_m[i]));

    // To make gaussians scale invariant
//    for(j=0;j<4;j++) prob_part -= sd_m[j];

    // M2 probability from model
    for(j=0;j<4;j++) prob_part += gauss_prob(M0[j],C[j],sd[j],M2_m[i]);

    // To account for inclination preference
    prob_part *= sin(inc_m[i]);

    // Determine total probability
    p1 += log10(prob_part);
  }

  
  
  
  // Find new points
  
  for(i=0;i<n_objs;i++){
  
    M1_new[i] = M1[i];
    Porb_new[i] = Porb[i];
    K_new[i] = K[i];

    // Select a new inclination angle
//    inc_new[i] = get_inc(seed);
    inc_new[i] = gauss_ran(inc_m[i],0.01,seed);
  
    // Select new M2 by root finding    
    Mf_temp = Porb[i]*K[i]*K[i]*K[i]/(2.0*PI*GGG);
    M2_new[i] = rtsafe(find_M2,M2_min[i],10.0,1.0e-5,M1[i],Mf_temp,inc_new[i]);
    
    // inc_new[i] = asin(K_new[i]/M2_new[i] * pow((M1_new[i]+M2_new[i])*(M1_new[i]+M2_new[i]) * Porb_new[i]/(2.0*PI*GGG),1.0/3.0));



  }


  // Select new coefficient distribution parameters
  
  for(i=0;i<4;i++){
    do{sd_new[i] = gauss_ran(sd[i],SD_DELTA,seed); }while(sd_new[i] < 0.0);
  }

  // Select new masses
  // Fix M0[3],M0_new[3] at 1.4 Msun
//  M0_new[3] = M0[3];
  do{
    for(i=0;i<3;i++){
      do{ M0_new[i] = gauss_ran(M0[i],M_DELTA,seed);}while(M0_new[i] < 0.0);
    }
    do{ M0_new[3] = gauss_ran(M0[3],M_DELTA,seed); }while(M0_new[3]< 1.3);
  }while(M0_new[0]>M0_new[1] || M0_new[1]>M0_new[2] || M0_new[2]>M0_new[3]);
  
  // Select coefficients
  for(i=0;i<4;i++){
    do{C_new[i] = gauss_ran(C[i],C_DELTA,seed);} while(C_new[i]<0.0);
  }    
  C_tot = C_new[0] + C_new[1] + C_new[2] + C_new[3];
  for(i=0;i<4;i++) C_new[i] = C_new[i] / C_tot;
  
    



  // Determine probability for next point

  p2 = 0.0;
  for(i=0;i<n_objs;i++){
    prob_part = 0.0;

//    prob_part += log10(gauss_prob(M1[i],1.0,M1_err[i],M1_new[i]));
//    prob_part += log10(gauss_prob(Porb[i],1.0,Porb_err[i],Porb_new[i]));
//    prob_part += log10(gauss_prob(K[i],1.0,K_err[i],K_new[i]));

    // To make gaussians scale invariant
//    for(j=0;j<4;j++) prob_part -= sd_new[j];

    // M2 probability from model
    for(j=0;j<4;j++) prob_part += gauss_prob(M0_new[j],C_new[j],sd_new[j],M2_new[i]);
    
    

    // To account for inclination preference
    prob_part *= sin(inc_new[i]);

    
    // Determine total probability
    p2 += log10(prob_part);


  }
 


  // Determine if move onto the new point
  // a logarithmic version
  
  ran_val = ran3(seed);
  if(ran_val < pow(10.0, p2 - p1)){

    for(i=0;i<4;i++){
      M0[i] = M0_new[i];
      C[i] = C_new[i];
      sd[i] = sd_new[i];
    }
    
    inc_m = inc_new;
    M1_m = M1_new;
    M2_m = M2_new;
    Porb_m = Porb_new;
    K_m = K_new;

    if(print_flag) N_accept ++;
  }else {
    if(print_flag) N_reject ++;
  }

  if(print_flag){
    *OUT << setprecision(4) << l << " " << p1 << " " << p2 << " ";
    *OUT << M0[0] << " " << C[0] << " " << sd[0] << " ";
    *OUT << M0[1] << " " << C[1] << " " << sd[1] << " ";
    *OUT << M0[2] << " " << C[2] << " " << sd[2] << " ";
    *OUT << M0[3] << " " << C[3] << " " << sd[3] << " ";

    *OUT << M0_new[0] << " " << C_new[0] << " " << sd_new[0] << " ";
    *OUT << M0_new[1] << " " << C_new[1] << " " << sd_new[1] << " ";
    *OUT << M0_new[2] << " " << C_new[2] << " " << sd_new[2] << " ";
    *OUT << M0_new[3] << " " << C_new[3] << " " << sd_new[3] << " ";
    
    
    *OUT << endl;
  }


}





void print_results(double best_M[4],double best_C[4],double best_sd[4]){
  int i;
  
  ofstream OUT;
  OUT.open("dist_Mf_inc.dat");

  cout << "Best fit is: " << endl;
  
  for (i=0;i<4;i++){
    cout << "Gaussian " << i << ": Centered at M=" << best_M[i];
    cout << " with a standard deviation, sigma=" << best_sd[i];
    cout << " and a normalized height of c=" << best_C[i] << endl;
  }


}


*/

