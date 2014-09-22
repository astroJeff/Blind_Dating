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

  int j,k;
  int n_objs;

  long i,seed;

  double P1,P2,Prob;
  double frac_NS,frac_NS_new;

  double mu[N_GAUSS];
  double sd[N_GAUSS];
  double w[N_GAUSS];
  double mu_new[N_GAUSS];
  double sd_new[N_GAUSS];
  double w_new[N_GAUSS];

  ofstream OUT;
  
  // Observations
  vector<string> Names;
  vector<double> M1;
  vector<double> M1_err;
  vector<double> K;
  vector<double> K_err;
  vector<double> Porb;
  vector<double> Porb_err;
  vector<double> M2_min;


  
  seed = time(NULL);
  
  OUT.open("dist_params.dat");


  // Take in data
  read_data(Names,M1,M1_err,M2_min,K,K_err,Porb,Porb_err);

  
  // Declare and initialize model variables
  n_objs = M1.size();
  vector<double> p_NS(n_objs,0.0);
  vector<double> p_NS_new(n_objs,0.0);
  vector<double> M2(n_objs,0.0);
  vector<double> inc(n_objs,0.0);
  vector<int> C(n_objs,0);


  // Give initial values
  initial_guess(mu,sd,w,p_NS,&frac_NS,inc,M2,M2_min,Names,M1,Porb,K,&seed);

  OUT << "Chain iteration mu_1 sd_1 w_1 w_NS Prob P(NS)_j" << endl;

  // k number of chains
  for(k=0;k<10;k++){


    // Metropolis-Hastings MCMC algorithm in here    
    // Burn in period
    P1 = P_model(mu,sd,w,p_NS,frac_NS,inc,M2,C,M2_min,Names,M1,K,Porb,&seed);
    for(i=0;i<BURN_IN;i++){
      next_point(mu,sd,w,&frac_NS,mu_new,sd_new,w_new,&frac_NS_new,&seed);
      P2 = P_model(mu_new,sd_new,w_new,p_NS_new,frac_NS_new,inc,M2,C,M2_min,Names,M1,K,Porb,&seed);
            
      if(ran3(&seed) < pow(10.0,P2-P1)) move_to_point(mu,sd,w,p_NS,&frac_NS,&P1,mu_new,sd_new,w_new,p_NS_new,&frac_NS_new,&P2);
    }
    
    
    // Iterate for N_CALC times
    P1 = P_model(mu,sd,w,p_NS,frac_NS,inc,M2,C,M2_min,Names,M1,K,Porb,&seed);
    for(i=0;i<N_CALC;i++){
      next_point(mu,sd,w,&frac_NS,mu_new,sd_new,w_new,&frac_NS_new,&seed);
      P2 = P_model(mu_new,sd_new,w_new,p_NS_new,frac_NS_new,inc,M2,C,M2_min,Names,M1,K,Porb,&seed);

      OUT << k << " " << i << " ";
      for(j=0;j<N_GAUSS;j++) OUT << setprecision(5) << mu[j] << " " << sd[j] << " " << w[j] << " ";
      OUT << frac_NS << " " << P1 << " ";
      for(j=0;j<p_NS.size();j++) OUT << setprecision(3) << p_NS[j] << " ";
      OUT << endl;
      
      if(ran3(&seed) < pow(10.0,P2-P1)) move_to_point(mu,sd,w,p_NS,&frac_NS,&P1,mu_new,sd_new,w_new,p_NS_new,&frac_NS_new,&P2);
            
    }
 
  }

  OUT.close();

}


void initial_guess(double mu[N_GAUSS],double sd[N_GAUSS],double w[N_GAUSS],vector<double>& p_NS,double* frac_NS,vector<double>& inc,vector<double>& M2,vector<double>& M2_min,vector<string>& Names,vector<double>& M1,vector<double>& Porb,vector<double>& K,long* seed){
  int i;
  
  double Mf;
  double model[3];

  for(i=0;i<N_GAUSS;i++){
    if(i==0){
      mu[0] = 0.90;    // High mass WD
      sd[0] = 0.20;
      w[0] = 0.90;
    }

    if(i==1){
      mu[1] = 0.6;    // Standard WD
      sd[1] = 0.2;
      w[1] = 0.46;
    }


    if(i==2){
      mu[2] = NS_MASS;    // NS
      sd[2] = 0.05;
      w[2] = 0.08;
    }
    
  }
  
  *frac_NS = 0.1;


// To set M2_min for each object

  for(i=0;i<M2.size();i++){
    
    Mf = Porb[i]*K[i]*K[i]*K[i]/(2.0*PI*GGG);

    model[0] = M1[i];
    model[1] = Mf;
    model[2] = PI/2.0;
    
    M2_min.push_back(rtsafe(find_M2,0.001,3.0,1.0e-6,model));
    
  }


/*

    inc[i] = get_inc(seed);

    cout << i << " " << K[i] << " " << M2_min[i] << endl;

    model[2] = inc[i];
    M2[i] = rtsafe(find_M2,M2_min[i],10.0,1.0e-5,model);

    if(Names[i] == "J0106-1000") M2[i] = 0.43;
    if(Names[i] == "J0651+2844") M2[i] = 0.50;
    if(Names[i] == "NLTT11748") M2[i] = 0.76;    

  }
*/
  
}


void next_point(double mu[N_GAUSS], double sd[N_GAUSS], double w[N_GAUSS], double* frac_NS,double mu_new[N_GAUSS],double sd_new[N_GAUSS],double w_new[N_GAUSS],double* frac_NS_new,long* seed){
  int i;
  double w_tot;
  
  for(i=0;i<N_GAUSS;i++){
    mu_new[i] = gauss_ran(mu[i],MOVE_MU,seed);
    sd_new[i] = gauss_ran(sd[i],MOVE_SD,seed);
    w_new[i] = gauss_ran(w[i],MOVE_W,seed);
    while(mu_new[i] < 0.0)  mu_new[i] = gauss_ran(mu[i],MOVE_MU,seed);
    while(sd_new[i] < 0.05)  sd_new[i] = gauss_ran(sd[i],MOVE_SD,seed);
    while(w_new[i] < 0.0)  w_new[i] = gauss_ran(w[i],MOVE_W,seed);
  }
  
  *frac_NS_new = gauss_ran(*frac_NS,MOVE_NS_FRAC,seed);
  while(*frac_NS_new < 0.0) *frac_NS_new = gauss_ran(*frac_NS,MOVE_NS_FRAC,seed);
  

  // Renormalize so weights add to unity
  w_tot = 0.0;
  for(i=0;i<N_GAUSS;i++) w_tot += w_new[i];
  
  for(i=0;i<N_GAUSS;i++) w_new[i] /= w_tot;  
  if(ADD_NS) for(i=0;i<N_GAUSS;i++) w_new[i] *= (1.0 - *frac_NS_new);
  
  
}


void move_to_point(double mu[N_GAUSS],double sd[N_GAUSS],double w[N_GAUSS],vector<double>& p_NS, double* frac_NS,double* P1,double mu_new[N_GAUSS],double sd_new[N_GAUSS],double w_new[N_GAUSS],vector<double>& p_NS_new,double* frac_NS_new,double* P2){
  int i;

  for(i=0;i<N_GAUSS;i++){
    mu[i] = mu_new[i];
    sd[i] = sd_new[i];
    w[i] = w_new[i];
  }
  
  (*frac_NS) = (*frac_NS_new);
  (*P1) = (*P2);


  // Move individual P(NS)'s - vector = is overloaded
  p_NS = p_NS_new;

}


double P_model(double mu[N_GAUSS],double sd[N_GAUSS],double w[N_GAUSS],vector<double>& p_NS,double frac_NS,vector<double>& inc,vector<double>& M2,vector<int>& C,vector<double>& M2_min,vector<string>& Names,vector<double>& M1,vector<double>& K,vector<double>& Porb,long* seed){
  int i,j,k;
  int n_objs;
  
  double P1,P2,P_NS_1,P_NS_2;
  double prob;
  double likelihood;
  double model[6+3*N_GAUSS];
  double NS_model[2];
  
  n_objs = K.size();


  likelihood = 0.0;
  
  for(i=0;i<N_GAUSS;i++){
    model[3*i+6] = mu[i];
    model[3*i+7] = sd[i];
    model[3*i+8] = w[i];
  }
  

  for(j=0;j<n_objs;j++){
  
    // Eclipsing systems
    // TEST THIS: Should the weights be unity????
    if(Names[j] == "J0106-1000"){
      if (ADD_NS) p_NS[j] = 0.0;
      likelihood += log10(gauss_prob(mu[0],w[0],sd[0],0.43));
      continue;
    }
    if(Names[j] == "J0651+2844"){
      if (ADD_NS) p_NS[j] = 0.0;
      likelihood += log10(gauss_prob(mu[0],w[0],sd[0],0.50));
      continue;
    }
    if(Names[j] == "J0345+1748"){
      if (ADD_NS) p_NS[j] = 0.0;
      likelihood += log10(gauss_prob(mu[0],w[0],sd[0],0.76));
      continue;
    }

  
  
    model[0] = M1[j];
    model[1] = K[j];
    model[2] = Porb[j];
    model[3] = 1.0;   // This is where "A" goes
    model[4] = M2_min[j];
    model[5] = frac_NS;  // Fraction of NS

    if(ADD_NS){
      NS_model[0] = frac_NS;
      NS_model[1] = M1[j];
    }
    
    // We want P(M2_min|mu,sd)
    P1 = prob_M2_model(model);
    if (ADD_NS) P_NS_1 = prob_NS(frac_NS,M1[j],M2_min[j]); 
     
    // Need to integrate over all M2_min
    qromb(prob_M2_wrapper,0.001,1.4,1.0e-4,&P2,model);
    if (ADD_NS) P_NS_2 = prob_NS_all(frac_NS,M1[j]);

        
    likelihood += log10( (P1+P_NS_1) / (P2 + P_NS_2) );

    if (ADD_NS) p_NS[j] = P_NS_1 / (P_NS_1 + P1);
  }
  
//  cout << P1 << " " << P2 << " " << likelihood << endl;

  return likelihood;

}



double prob_M2_model(double model[6+3*N_GAUSS]){
  //  P( M2[j] | mu, sd, w, M2_min[j] )

  int i,j;
  int nsteps;
  
  double M1,M2,A,K,Porb,M2_min;
  double mu,sd,w;
  double M2_prob;

/*
  model[0] = M1;
  model[1] = K;
  model[2] = Porb;
  model[3] = 1.0;   // This is where "A" goes
  model[4] = M2_min;
  model[5] = 0.0;


  model[6] = mu;
  model[7] = sd;
  model[8] = w;

*/

//  if(Name == "J0106-1000") return 0.43;
//  if(Name == "J0651+2844") return 0.50;
//  if(Name == "NLTT11748") return 0.76;

  // Normalize the function
  M2_prob = 0.0;
  M2_min = model[4];
  
  // Integrate over full range [M2_min,5.0] to get total probability
  qromb(eval_M2,M2_min,3.0,1.0e-5,&M2_prob,model);
  
  return M2_prob;

}


double prob_M2_wrapper(double M2_min, double* model){ 
  
  model[4] = M2_min;  

  return prob_M2_model(model);  
}



double prob_NS(double frac_NS, double M1, double M2_min){
  // P(NS) = sin(i)
  double M2 = NS_MASS;
  double prob;

  prob = frac_NS * pow((M1+M2)/(M1+M2_min), 2.0/3.0) * M2_min / M2;
 
  return (M2_min < M2) ? prob : 0.0;
}


double prob_NS_all(double frac_NS, double M1){
  // This is (should be) the analytic integral over all M2_min's for a given M1

  double xmin,xmax;
  double M2;
  double a1,a2,a3,a4;

  xmin = 0.0;
  xmax = NS_MASS;  
  M2 = NS_MASS;

//  a1 = 3.0 * xmax / M2 * pow((M1+M2)*(M1+M2)*(M1+xmax),1.0/3.0);
//  a2 = -3.0 * xmin / M2 * pow((M1+M2)*(M1+M2)*(M1+xmin),1.0/3.0);

//  a3 = -9.0/4.0 * pow(M1+M2,2.0/3.0) / M2 * pow(M1+xmax,4.0/3.0);
//  a4 = 9.0/4.0 * pow(M1+M2,2.0/3.0) / M2 * pow(M1+xmin,4.0/3.0);
  
  a1 = -3.0/(4.0*M2) * (3.0*M1-xmax) * (M1+xmax) * pow((M1+M2)/(M1+xmax),2.0/3.0);
  a2 = -3.0/(4.0*M2) * (3.0*M1-xmin) * (M1+xmin) * pow((M1+M2)/(M1+xmin),2.0/3.0);
  
  return frac_NS*(a1-a2);
//  return frac_NS*(a1+a2+a3+a4);
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
  double M1,K,Porb,A,M2_min;
  double sini;
  double temp;
  double frac_NS;


  M1 = model[0];     // Observed WD mass
  K = model[1];      // Line of sight orbital velocity
  Porb = model[2];   // Orbital period
  A = model[3];      // Normalization constant
  M2_min = model[4];  


  for(i=0;i<N_GAUSS;i++){
    mu[i] = model[3*i+6];
    sd[i] = model[3*i+7];
    w[i] = model[3*i+8];
  }

  temp = 0.0;
  for(i=0;i<N_GAUSS;i++){  


//    sini = pow(Porb/(2.0*PI*GGG) * (M1+M2)*(M1+M2), 1.0/3.0) * K / M2;
    sini = pow((M1+M2)/(M1+M2_min) , 2.0/3.0) * M2_min/M2;

    temp += A * gauss_prob(mu[i],w[i],sd[i],M2) * sini;

//    cout << i << " " << sini << " " << mu[i] << " " << sd[i] << " " << w[i] << " " << temp << endl;

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




void read_data(vector<string>& Names,vector<double>& M1,vector<double>& M1_err,vector<double>& M2_min,vector<double>& K,vector<double>& K_err,vector<double>& Porb,vector<double>& Porb_err){
  fstream IN;
  
  double M1_temp, M1_err_temp;
  double K_temp, K_err_temp;
  double P_temp, P_err_temp;
  double M2_min_temp;
  
  string Name_temp,line;
  
  vector<string> data;
  



  // Open data file
  IN.open("ELM_WD.dat");
//  IN.open("dist_1_gauss.dat");
//  IN.open("dist_07_NS.dat");
//  IN.open("dist_flat_NS.dat");

//  getline(IN,line);

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
    
//    Names.push_back(data[3]);
//    Porb.push_back(strtof(data[1]));
//    M1.push_back(strtof(data[0]));
//    K.push_back(strtof(data[2]));
//    M2_min.push_back(strtof(data[4]));
    

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




double poisson_prob(double lambda, double n){
  return pow(lambda,n)*exp(-lambda)/factorial(n);
}

double poisson_ran(double lambda,long* seed){
  int k;
  
  double ran_val,prob;
  
  k = 0;
  prob = 0.0;
  ran_val = ran3(seed);

  while(prob < ran_val){
    prob += exp(-lambda) * pow(lambda,k) / factorial(k);
    k ++;
  }
  k --;

  return k;
}



double factorial(double n){
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}






void find_range(vector<double>& in_array,double* x_out,double* x_l,double* x_u){

  long i;

  sort(in_array,N_CALC);

  // Get the value one sigma away from the center
  *x_l = in_array[(long)(N_CALC*0.159)];
  *x_u = in_array[(long)(N_CALC*0.841)];
  *x_out = in_array[(long)(N_CALC*0.5)];
}





