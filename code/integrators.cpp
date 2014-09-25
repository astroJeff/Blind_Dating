#include<iostream>
#include<cmath>
#include<vector>
#include"integrators.h"

using namespace std;

void trapzd(double (*func)(double, double *), double x1, double x2, double* y, int nsteps, double* model){

  int it,j;
  double dx,sum,tnm,x;
  
  if(nsteps==1){
    *y = 0.5 * (x2-x1) * (func(x1,model) + func(x2,model));
  }else {

    it = pow(2,nsteps-2);
    tnm = it;
    dx = (x2-x1) / tnm;   // Width of points
    x = x1 + 0.5*dx;
    
    sum = 0.0;
    for(j=0;j<it;j++){
      sum = sum + func(x,model);
      x = x + dx;
    }
    
    *y = 0.5 * (*y + (x2-x1)*sum/tnm);
  
  }

}



void qtrap(double (*func)(double, double *), double x1, double x2, double eps, double* y, double* model){

  int j,n,jmax;
  
  double olds;
  
  olds = -1.0e30;
  jmax = 20;
  
  for(j=0;j<jmax;j++){
  
    // Call trapzd with increasingly more steps taken
    n = j+1;
    trapzd(func,x1,x2,y,n,model);

    // Check if error target has been reached
    if (j > 4) {
      if ((abs(*y-olds) < eps*abs(olds)) || (*y==0 && olds==0)) return;
    }

    olds = *y;  
  }
  
  cerr << "Too many steps in qtrap" << endl;
}



void qromb(double (*func)(double, double *), double x1, double x2, double eps, double* y, double* model){

  int i,j,n,jmax;

  double y_err;
  double h[21];   // The 21 corresponds to jmax(=20)+1
  double s[21];  

  h[0] = 1.0;  
  jmax = 20;
  y_err = 0.0;

  for(j=0;j<jmax;j++){
    n = j+1;
    trapzd(func,x1,x2,&s[j],n,model);

 
    if (j > 3){
      polint2(h,s,n,n,0.0,y,&y_err);
      if (abs(y_err) < eps*abs(*y)) return;
    }
  
    s[j+1] = s[j];
    h[j+1] = 0.25*h[j];  // The 0.25 factor is important

  }

//  cerr << "Too many steps in qromb" << endl;
  cerr << "Too many steps in qromb " << model[6]<< " " << model[7] << " " << model[8] << endl;

}





void polint2(double* xa,double* ya,int n,int pts,double x,double* y,double* dy){
  /* 
    xa = array of x's, ya = array of y's, n = length of arrays, 
    pts = data pts to interpolate (probably more than 4), 
    x = x-value to interpret, y = output value, dy = error estimate
  */

  int i,m,ns;
  int index;
  
  double diff,diff_best;
  double den,dif,dift,ho,hp,w;
  double *xp,*yp,*c,*d;

  xp = new double[pts];
  yp = new double[pts];
  c = new double[pts];
  d = new double[pts];
  
  ns = 0;
  diff = abs(x-xa[0]);

  // Looking for the closest table entry
  // The index is stored in ns
  for(i=0;i<n;i++){
    diff_best = abs(x-xa[i]);
    if (diff_best < diff){
      ns = i;
      diff = diff_best;
    }
  }
  
  // NEED TO MAKE THIS FLEXIBLE FOR EXTRAPOLATIONS
  if( ns - (int)((double)pts/2.0) <= 0 ){
    index = 0;
  }else if( ns - (int)((double)pts/2.0) + pts >= n ){
    index = n - pts;
  }else {
    index = ns - (int)((double)pts/2.0);
  }

  // Create new array around closest table entry  
  for(i=0;i<pts;i++){
    xp[i] = xa[index+i];
    yp[i] = ya[index+i];
    c[i] = ya[index+i];
    d[i] = ya[index+i];
  }
  
  
  // Adjust ns to index around closest entry in new array 
  for(i=0;i<pts;i++) if ( (int)((double)pts/2.0) == i ) ns = i;
  
  // Starting value;
  *y = yp[ns];
  ns = ns-1;

  for(m=1;m<pts;m++){

    for(i=0;i<pts-m-1;i++){
      ho = xp[i] - x;
      hp = xp[i+m] - x;
      w = c[i+1] - d[i];
      den = ho-hp;
      if (den == 0) cerr << " Failure in polint" << endl;  
  
      den = w/den;
      d[i] = hp*den;
      c[i] = ho*den;
    }
    
    if (2*ns < pts-m) {
      *dy = c[ns+1];
    }else {
      *dy = d[ns];
      ns = ns-1;
    }
    
//    cout << m << " " << *dy << endl;
    
    *y = *y + *dy;
  }
  
  delete [] xp;
  delete [] yp;
  delete [] c;
  delete [] d;
}







void polint(double* xa,double* ya,int n,double x,double* y,double* dy){
  /* 
    xa = array of x's, ya = array of y's, n = length of arrays, 
    x = x-value to interpret, y = output value, dy = error estimate
  */

  int i,m,ns;
  
  double diff,diff_best;
  double den,dif,dift,ho,hp,w;
  double *c,*d;

  c = new double[n];
  d = new double[n];
  
  ns = 0;
  diff = abs(x-xa[0]);

  // Looking for the closest table entry
  // The index is stored in ns
  for(i=0;i<n;i++){
    diff_best = abs(x-xa[i]);
    if (diff_best < diff){
      ns = i;
      diff = diff_best;
    }
    
    c[i] = ya[i];
    d[i] = ya[i];
  }
  
  *y = ya[ns];
  ns = ns-1;

  for(m=1;m<n;m++){

    for(i=0;i<n-m-1;i++){
      ho = xa[i] - x;
      hp = xa[i+m] - x;
      w = c[i+1] - d[i];
      den = ho-hp;
      if (den == 0) cerr << " Failure in polint" << endl;  
  
      den = w/den;
      d[i] = hp*den;
      c[i] = ho*den;
    }
    
    if (2*ns < n-m) {
      *dy = c[ns+1];
    }else {
      *dy = d[ns];
      ns = ns-1;
    }
    
//    cout << m << " " << *dy << endl;
    
    *y = *y + *dy;
  }
  
  delete [] c;
  delete [] d;
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
        
    rk4(v,dv,nvar,x,h,vout,model,derivs);
    
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

double rtsafe(void (*funcd)(double,double *,double *,double *),double x1,double x2,double xacc,double* model){
  int j;
  double df,dx,dxold,f,fh,fl;
  double temp,xh,xl,rts;
  
  (*funcd)(x1,&fl,&df,model);
  (*funcd)(x2,&fh,&df,model);
  if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) cerr << "Root must be bracketed in rtsafe" << endl;
  
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


void sort(vector<double>& x,long n_objs){
  // Heapsort algorithm

  long e,i,j,k,l;

  double tmp;

  e = n_objs-1;

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


double ran3(long *idum){ 
  static int inext, inextp;
  static long ma[56];
  static int iff=0;
  long mj, mk;
  int i, ii, k;
  
  if(*idum < 0 || iff == 0){
    iff=1;
    mj = (long)abs(MSEED-abs((double)*idum));
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
