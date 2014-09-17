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
  vector<double> M2_min;
  


  
  seed = time(NULL);
  
  OUT.open("dist_params.dat");


  // Take in data
  read_data(Names,M1,M1_err,M2_min,K,K_err,Porb,Porb_err);

  
  // Declare and initialize model variables
  n_objs = M1.size();


  vector<double> M2(n_objs,0.0);
  vector<double> inc(n_objs,0.0);
  vector<int> C(n_objs,0);


  // Give initial values
  initial_guess(mu,sd,w,inc,M2,M2_min,Names,M1,Porb,K,&seed);



  // k number of walkers
  for(k=0;k<1;k++){



    // Iterate for BURN_IN times for "burn in"
    print_flag = 0;  // Don't print out for "burn in"
    
 //   for(i=0;i<BURN_IN;i++) iterate(mu,sd,w,inc,M2,C,M2_min,Names,M1,K,Porb,print_flag,k,&OUT,&seed);


    print_flag = 1;  // Print out for iterations


    // Iterate for N_CALC times
    mu[0] = 0.60;
    for(i=0;i<10;i++){
      mu[0] += 0.02;

        
      sd[0] = 0.10;
      for(j=0;j<10;j++){
//    for(i=0;i<N_CALC;i++){
        sd[0] += 0.02;
    
    
      // MCMC algorithm in here    
//      mod_mix_gaus_gibbs(mu,sd,w,inc,M2,C,M2_min,Names,M1,K,Porb,print_flag,k,&OUT,&seed);
        iterate(mu,sd,w,inc,M2,C,M2_min,Names,M1,K,Porb,print_flag,k,&OUT,&seed);
      }
    
    }
  }



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
//      w[0] = 0.25;
      w[0] = 1.0;
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

/*  
  for(i=0;i<M2.size();i++){

    inc[i] = get_inc(seed);

    Mf = Porb[i]*K[i]*K[i]*K[i]/(2.0*PI*GGG);

    model[0] = M1[i];
    model[1] = Mf;
    model[2] = PI/2.0;
    M2_min[i] = rtsafe(find_M2,0.1,10.0,1.0e-5,model);

    cout << i << " " << K[i] << " " << M2_min[i] << endl;

    model[2] = inc[i];
    M2[i] = rtsafe(find_M2,M2_min[i],10.0,1.0e-5,model);

    if(Names[i] == "J0106-1000") M2[i] = 0.43;
    if(Names[i] == "J0651+2844") M2[i] = 0.50;
    if(Names[i] == "NLTT11748") M2[i] = 0.76;    

  }
*/
  
}



void iterate(double mu[N_GAUSS],double sd[N_GAUSS],double w[N_GAUSS],vector<double>& inc,vector<double>& M2,vector<int>& C,vector<double>& M2_min,vector<string>& Names,vector<double>& M1,vector<double>& K,vector<double>& Porb,int print_flag,int l,ofstream* OUT,long* seed){
  int i,j,k;
  int n_objs;
  
  double P1,P2;
  double prob;
  double likelihood;
  double model[9];
  
  n_objs = K.size();


  likelihood = 0.0;

  for(j=0;j<n_objs;j++){
    model[0] = M1[j];
    model[1] = K[j];
    model[2] = Porb[j];
    model[3] = 1.0;   // This is where "A" goes
    model[4] = M2_min[j];
    model[5] = 0.0;


    for(i=0;i<N_GAUSS;i++){
  
      model[6] = mu[i];
      model[7] = sd[i];
      model[8] = w[i];
    
      // We want P(M2_min|mu,sd)
      P1 = prob_M2_model(model);
    
    
      // Need to integrate over all M2_min
      qromb(prob_M2_wrapper,0.001,1.5,1.0e-4,&P2,model);
        
      likelihood += log10(P1/P2);

    }    

  }


  cout << mu[0] << " " << sd[0] << " " << likelihood << endl;

}



double prob_M2_wrapper(double M2_min, double* model){
  
  model[4] = M2_min;  

  return prob_M2_model(model);

}


double prob_M2_model(double model[9]){
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








void prob_mix_gauss(double mu[N_GAUSS],double sd[N_GAUSS],double w[N_GAUSS],vector<double>& M2,vector<int>& C,long* seed){  
  // To determine P( mu, sigma, w | M2, C )

  int i,j;
  int n_objs,n_C;
  
  double w_tot;
  double mean;
  
  vector<double> obj_M2;  
  
  n_objs = M2.size();
  w_tot = 0.0;


  for(i=0;i<N_GAUSS;i++){
  
  
    // Create vector of values
    obj_M2.clear();
    for(j=0;j<n_objs;j++) if(C[j] == i) obj_M2.push_back(M2[j]);
    
    
    
    // Determine mean and standard deviation
    find_gauss(&mu[i],&sd[i],obj_M2);
    
    
    
    // Determine weight - draw randomly from Poisson Distribution
    w[i] = poisson_ran(obj_M2.size(),seed);

    w_tot += w[i];
  }

  // Normalize so weights all add to unity
  for(i=0;i<N_GAUSS;i++) w[i] /= w_tot;
    

  cout << mu[0] << " " << sd[0] << " " << w[0] << endl;
  



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
  //  IN.open("ELM_WD.dat");
  IN.open("dist_1_gauss.dat");
  getline(IN,line);

  // Read in data, add M2,min to vector M2
  while(getline(IN,line)){

    split(line,data);

    /*
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
    */

    
    Names.push_back(data[3]);
    Porb.push_back(strtof(data[1]));
    M1.push_back(strtof(data[0]));
    K.push_back(strtof(data[2]));
    M2_min.push_back(strtof(data[4]));
    

  }
  
  IN.close();
}

double get_inc(long *seed){
  return acos(1.0 - ran3(seed));
}




void find_gauss(double* mean,double* sd,vector<double> vals){
  int i, n_objs;
	
  n_objs = vals.size();

  *mean = 0.0;
  for(i=0;i<n_objs;i++){
    *mean += vals[i];
  }
  *mean = *mean/(double)vals.size();


  *sd = 0.0;
  for(i=0;i<n_objs;i++) *sd += (*mean - vals[i]) * (*mean - vals[i]);
  
  *sd = sqrt( 1.0 / (double(n_objs) - 1.0) ) * sqrt(*sd);

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





