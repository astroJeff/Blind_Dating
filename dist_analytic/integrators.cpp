#include<iostream>
#include<cmath>
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

  cerr << "Too many steps in qromb" << endl;

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
