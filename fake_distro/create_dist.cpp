#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<vector>
#include<cmath>
#include"create_dist.h"

using namespace std;

long N_accept=0;
long N_reject=0;

int main(int argc, char* argv[]){

  int j,index;
  long i,seed;

  double temp_ran,temp_prob;
  double M1, Porb, vel, inc, Mf, M2, M2_min;
  
  double model[3];

  double *w,*sd,*M0;


  ofstream OUT;
  
  OUT.open("dist_gaus.dat");
  
  seed = time(NULL);




  if(WD_DIST == 1){
    M0 = new double[N_GAUSS];
    sd = new double[N_GAUSS];
    w = new double[N_GAUSS];
    create_gauss_dist(M0,sd,w);
  }
  
  OUT << "M1  Porb K  M2(actual) M2_min" << endl;


  // System parameters
  M1 = 0.25;
  Porb = 0.5;  // Assume a half-day orbit
  
  
  // Iterate for N_CALC times
  for(i=0;i<N_CALC;i++){
    
    temp_ran = ran3(&seed);
    
    
    if (temp_ran < NS_RATE){    // NS Distribution

      M2 = gauss_ran(NS_M,NS_SD,&seed);
      
    }else {                     // WD Distribution

      if(WD_DIST == 1){
	index = 0;
	temp_prob = 0.0;     
	temp_ran = ran3(&seed);
	
	// To pick which gaussian to choose from
	for(j=0;j<N_GAUSS;j++){
	  temp_prob += w[j];
	  if(temp_ran>temp_prob)  index = j+1;
	}    
	
	M2 = gauss_ran(M0[index],sd[index],&seed);
	
      }else {
	M2 = ran3(&seed) * (WD_M_MAX-WD_M_MIN) + WD_M_MIN;
      }
      
    }
    
    // Check with M2_min distribution??
    // Get a random inclination angle
    inc = get_inc(&seed);
    
    // Determine what the relative orbital velocity will be
    vel = M2 * sin(inc) * pow( 1.0 / ((M1+M2)*(M1+M2)) * 2.0 * PI * GGG / Porb, 1.0/3.0);
    
    // Then get the mass function
    Mf = Porb * vel*vel*vel / (2.0*PI*GGG);
    
    // Then solve for the apparent M2_min
    model[0] = M1;
    model[1] = Mf;
    model[2] = PI/2.0;  // inclination set to 90 degrees for M2_min
    M2_min = rtsafe(find_M2,0.0,10.0,1.0e-5,model);
    
    OUT << M1 << " " << Porb << " " << vel << " " << M2 << " " << M2_min << endl;
    
  }
  
  OUT.close();
  
  if(WD_DIST == 1){
    delete[] M0;
    delete[] sd;
    delete[] w;
  }

}


void create_gauss_dist(double* M,double* sd,double* w){
  int i;
  double weight;
  
  // Gaussian mixture model
  for(i=0;i<N_GAUSS;i++){
      
      if(i==0){
	M[i] = 0.7;
	w[i] = 0.4;
	sd[i] = 0.2;
      }
      
      if(i==1){
	M[i] = 0.9;
	w[i] = 0.2;
	sd[i] = 0.05;
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

