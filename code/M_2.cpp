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
  int n_objs, accept;

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
  vector<double> K;
  vector<double> Porb;
  vector<double> M2_min;
  vector<double> inc_max;
  vector<string> inc_name;


  
  seed = time(NULL);
  
  OUT.open("post_09_NS.dat");


  // Take in data
  read_data(Names,M1,K,Porb,inc_max,inc_name);

  
  // Declare and initialize model variables
  n_objs = M1.size();
  vector<double> p_NS(n_objs,0.0);
  vector<double> p_NS_new(n_objs,0.0);
  vector<double> M2(n_objs,0.0);


  // Give initial values
  initial_guess(mu,sd,w,p_NS,&frac_NS,inc_max,inc_name,M2,M2_min,Names,M1,Porb,K,&seed);

  OUT << "Chain iteration mu_1 sd_1 w_1 w_NS Prob P(NS)_j" << endl;

  // k number of chains
  for(k=0;k<10;k++){


    // Metropolis-Hastings MCMC algorithm in here    
    // Burn in period
    P1 = P_model(mu,sd,w,p_NS,frac_NS,M2,M2_min,Names,M1,K,Porb,&seed);
    for(i=0;i<BURN_IN;i++){
      next_point(mu,sd,w,&frac_NS,mu_new,sd_new,w_new,&frac_NS_new,&seed);
      P2 = P_model(mu_new,sd_new,w_new,p_NS_new,frac_NS_new,M2,M2_min,Names,M1,K,Porb,&seed);
            
      if(ran3(&seed) < pow(10.0,P2-P1)) move_to_point(mu,sd,w,p_NS,&frac_NS,&P1,mu_new,sd_new,w_new,p_NS_new,&frac_NS_new,&P2);
    }
    
    
    // Iterate for N_CALC times
    accept = 0;
    P1 = P_model(mu,sd,w,p_NS,frac_NS,M2,M2_min,Names,M1,K,Porb,&seed);
    for(i=0;i<N_CALC;i++){
      next_point(mu,sd,w,&frac_NS,mu_new,sd_new,w_new,&frac_NS_new,&seed);
      P2 = P_model(mu_new,sd_new,w_new,p_NS_new,frac_NS_new,M2,M2_min,Names,M1,K,Porb,&seed);

      OUT << k+1 << " " << i+1 << " ";
      for(j=0;j<N_GAUSS;j++) OUT << setprecision(5) << mu[j] << " " << sd[j] << " " << w[j] << " ";
      OUT << frac_NS << " " << P1 << " " << accept << " ";
      for(j=0;j<p_NS.size();j++) OUT << setprecision(3) << p_NS[j] << " ";
      OUT << endl;

      accept = 0; 
      if(ran3(&seed) < pow(10.0,P2-P1)){
        move_to_point(mu,sd,w,p_NS,&frac_NS,&P1,mu_new,sd_new,w_new,p_NS_new,&frac_NS_new,&P2);
        accept = 1;      
      }    
 
    }
 
  }

  OUT.close();

}


void initial_guess(double mu[N_GAUSS],double sd[N_GAUSS],double w[N_GAUSS],vector<double>& p_NS,double* frac_NS,vector<double>& inc_max,vector<string>& inc_name,vector<double>& M2,vector<double>& M2_min,vector<string>& Names,vector<double>& M1,vector<double>& Porb,vector<double>& K,long* seed){

  int i,j;
  
  double Mf;
  double model[3];

  for(i=0;i<N_GAUSS;i++){
    if(i==0){
      mu[0] = 0.90;    // High mass WD
      sd[0] = 0.30;
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
    

    // To use constraints based on non-detection of eclipses
    if(USE_INC){
      for(j=0;j<inc_max.size();j++){
	if(!Names[i].compare(inc_name[j])){
	  model[2] = inc_max[j] * PI / 180.0;
	  M2_min[i] = rtsafe(find_M2,0.001,3.0,1.0e-6,model);
	}
      }
    }

    // To require that M2 > M1 because the observed WD must
    // be the less massive one (inverse M-R relation)
    if(M2_GT_M1){
      if(M2_min[i] < M1[i]) M2_min[i] = M1[i];
    }
    
  }

  
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


double P_model(double mu[N_GAUSS],double sd[N_GAUSS],double w[N_GAUSS],vector<double>& p_NS,double frac_NS,vector<double>& M2,vector<double>& M2_min,vector<string>& Names,vector<double>& M1,vector<double>& K,vector<double>& Porb,long* seed){
  int i,j,k;
  int n_objs;
  
  double mf,mf_max;
  double mf1_3,Mtot4_3;
  double P1,P2,P_NS_1,P_NS_2;
  double prob;
  double likelihood;
  double model[7+3*N_GAUSS];
  
  n_objs = K.size();


  likelihood = 0.0;
  
  for(i=0;i<N_GAUSS;i++){
    model[3*i+7] = mu[i];
    model[3*i+8] = sd[i];
    model[3*i+9] = w[i];
  }
  

  for(j=0;j<n_objs;j++){
  
    // Eclipsing systems
    // TEST THIS: Should the weights be unity????
    if(Names[j] == "J0751-0141"){
      if (ADD_NS) p_NS[j] = 0.0;
      likelihood += log10(gauss_prob(mu[0],w[0],sd[0],0.97));
      continue;
    }
    if(Names[j] == "J0651+2844"){
      if (ADD_NS) p_NS[j] = 0.0;
      likelihood += log10(gauss_prob(mu[0],w[0],sd[0],0.50));
      continue;
    }
    if(Names[j] == "J0345+1748"){
      if (ADD_NS) p_NS[j] = 0.0;
      likelihood += log10(gauss_prob(mu[0],w[0],sd[0],0.72));
      continue;
    }


    mf = Porb[j] * K[j]*K[j]*K[j] / (2.0*PI*GGG);
  
    model[0] = M1[j];
    model[1] = K[j];
    model[2] = Porb[j];
    model[3] = mf;      // Mass function
    model[4] = M2_min[j];
    model[5] = frac_NS;  // Fraction of NS
    model[6] = pow(mf,1.0/3.0); // Mass function to the 1/3 power


    // We want P(M2_min|mu,sd)
    P1 = prob_M2_model(model);
    if (ADD_NS) P_NS_1 = prob_NS(model); 


    /*
    // Need to integrate over all mf
    mf_max = 2.0*2.0*2.0 / ((M1[j]+2.0) * (M1[j]+2.0));
    qromb(prob_Mf_wrapper,0.001,mf_max,1.0e-3,&P2,model);
    if (ADD_NS) P_NS_2 = prob_NS_all(frac_NS,M1[j]);
                 
    likelihood += log10( (P1+P_NS_1) / (P2 + P_NS_2) );
*/

    likelihood += log10( P1 + P_NS_1 );

    if (ADD_NS) p_NS[j] = P_NS_1 / (P_NS_1 + P1);
    
  }
  

  return likelihood;

}



double prob_M2_model(double* model){
  //  P( M2[j] | mu, sd, w, M2_min[j] )

  int i,j;
  int nsteps;
  
  double M1,M2,A,K,Porb,M2_min;
  double mu,sd,w;
  double M2_prob;


  // Normalize the function
  M2_prob = 0.0;
  M2_min = model[4];

  // Integrate over full range [M2_min,5.0] to get total probability
  qromb(eval_M2,M2_min+0.001,3.0,1.0e-3,&M2_prob,model);
  
  return M2_prob;

}


double prob_Mf_wrapper(double Mf, double* model){ 
  double M2_min,M1;
  double model_temp[3];
  double temp_prob;
  
  M1 = model[0];
  
  // To solve for M2_min for each Mf
  model_temp[0] = M1;
  model_temp[1] = Mf;
  model_temp[2] = PI/2.0;

  M2_min = rtsafe(find_M2,0.001,5.0,1.0e-4,model_temp);

  // Add each Mf, M2_min to the model
  model[3] = Mf;
  model[4] = M2_min;
  model[6] = pow(Mf,1.0/3.0);

  temp_prob = prob_M2_model(model);

  return temp_prob;  
}



double prob_NS(double* model){
  // P(NS) = sin(i)

  int i;
  
  double M1,M2,frac_NS,mf;
  double prob, M2_min;
  double NS_model[7+3*N_GAUSS];
  
  M1 = model[0];
  mf = model[3];
  M2_min = model[4];
  M2 = NS_MASS;
  frac_NS = model[5];



  if(NS_MODEL){
    // Modeling the NS distribution as a Gaussian

    for(i=0;i<7+3*N_GAUSS;i++) NS_model[i] = model[i];

    NS_model[7] = NS_MASS;
    NS_model[8] = 0.1;
    NS_model[9] = frac_NS;

    qromb(eval_M2,M2_min+0.001,3.0,1.0e-3,&prob,NS_model);

  }else {
    // Modeling the NS distribution as a delta function

    prob = frac_NS * pow(M1+M2,4.0/3.0) / (3.0*pow(mf,1.0/3.0) * M2 * sqrt(M2*M2 - pow(mf*(M1+M2)*(M1+M2),2.0/3.0)));
 
    return (M2_min < M2) ? prob : 0.0;
  }  



}


double prob_NS_all(double frac_NS, double M1){
  // This is (should be) the analytic integral over all M2_min's for a given M1

  double xmin,xmax;
  double M2;
  double a1,a2,a3,a4;

//  xmin = 0.0;
//  xmax = NS_MASS;  
//  M2 = NS_MASS;

//  a1 = 3.0 * xmax / M2 * pow((M1+M2)*(M1+M2)*(M1+xmax),1.0/3.0);
//  a2 = -3.0 * xmin / M2 * pow((M1+M2)*(M1+M2)*(M1+xmin),1.0/3.0);

//  a3 = -9.0/4.0 * pow(M1+M2,2.0/3.0) / M2 * pow(M1+xmax,4.0/3.0);
//  a4 = 9.0/4.0 * pow(M1+M2,2.0/3.0) / M2 * pow(M1+xmin,4.0/3.0);
  
//  a1 = -3.0/(4.0*M2) * (3.0*M1-xmax) * (M1+xmax) * pow((M1+M2)/(M1+xmax),2.0/3.0);
//  a2 = -3.0/(4.0*M2) * (3.0*M1-xmin) * (M1+xmin) * pow((M1+M2)/(M1+xmin),2.0/3.0);
  


  xmin = 0.0;
  xmax = NS_MASS*NS_MASS*NS_MASS/((M1+NS_MASS) * (M1+NS_MASS));
  M2 = NS_MASS;

  a1 = - 1.0/M2 * sqrt(M2*M2 - pow(xmax,2.0/3.0)*pow(M1+M2,4.0/3.0));
  a2 = - 1.0/M2 * sqrt(M2*M2 - pow(xmin,2.0/3.0)*pow(M1+M2,4.0/3.0));

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
  double M1,K,Porb,A,M2_min,mf;
  double Mtot4_3,mf1_3;
  double sini;
  double temp;
  double frac_NS;


  M1 = model[0];     // Observed WD mass
  K = model[1];      // Line of sight orbital velocity
  Porb = model[2];   // Orbital period
  mf = model[3];     // Mass function
  M2_min = model[4];  
  frac_NS = model[5];
  mf1_3 = model[6];

  Mtot4_3 = pow(M1+M2,4.0/3.0);



  for(i=0;i<N_GAUSS;i++){
    mu[i] = model[3*i+7];
    sd[i] = model[3*i+8];
    w[i] = model[3*i+9];
  }


  temp = 0.0;
  for(i=0;i<N_GAUSS;i++){  

    if(frac_NS == w[i]){
      temp += trunc_gauss(mu[i],w[i],sd[i],M2,1.3,2.0) * Mtot4_3 / ((3.0*mf1_3) * M2 * sqrt(M2*M2 - mf1_3*mf1_3 * Mtot4_3));
    }else {
      temp += trunc_gauss(mu[i],w[i],sd[i],M2,0.0,1.44) * Mtot4_3 / ((3.0*mf1_3) * M2 * sqrt(M2*M2 - mf1_3*mf1_3 * Mtot4_3));
    }
    //    temp += gauss_prob(mu[i],w[i],sd[i],M2) * Mtot4_3 / ((3.0*mf1_3) * M2 * sqrt(M2*M2 - mf1_3*mf1_3 * Mtot4_3));

  }

  return temp;
}



void read_data(vector<string>& Names,vector<double>& M1,vector<double>& K,vector<double>& Porb,vector<double>& inc_max,vector<string>& inc_name){

  int i, j;

  fstream IN;
  
  double M1_temp;
  double K_temp;
  double P_temp;

  string line;
  string Name_temp;
  
  vector<string> data;
  



  // Open data file
//  IN.open("ELM_WD.dat");
  IN.open("dist_09_NS.dat");

//  getline(IN,line);

  // Read in data, add M2,min to vector M2
  while(getline(IN,line)){

    split(line,data);

    
    Name_temp = data[2];
    P_temp = strtof(data[1]);
    K_temp = strtof(data[2]);
    M1_temp = strtof(data[0]);

    // Don't worry about non-detections yet
    if (P_temp > 0.001){
      Names.push_back(Name_temp);
      Porb.push_back(P_temp);
      K.push_back(K_temp);
      M1.push_back(M1_temp);
    }
    
//    Names.push_back(data[3]);
//    Porb.push_back(strtof(data[1]));
//    M1.push_back(strtof(data[0]));
//    K.push_back(strtof(data[2]));
//    M2_min.push_back(strtof(data[4]));
    

  }
  IN.close();


/*
  IN.open("ELM_inc.dat");
  
  // Get rid of header
  getline(IN,line);

  while(getline(IN,line)){

    split(line,data);

    inc_name.push_back(data[0]);
    inc_max.push_back(strtof(data[1]));
  }

  IN.close();
*/

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


double trunc_gauss(double val, double w, double val_err, double val_in, double g_min, double g_max){

  double prob, norm;

  prob = gauss_prob(val,w,val_err,val_in);
  
  norm = 0.5 * ( erf((g_max-val)/(sqrt(2.0)*val_err)) - erf((g_min-val)/(sqrt(2.0)*val_err)) );

  prob /= norm;
  
  return prob;

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





