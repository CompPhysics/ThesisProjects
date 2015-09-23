/*


	mpicxx Roothan_VMC.cpp nrutil.c methods.cpp lib.cpp dfpmin.c lnsrch.c -o program.x   
	mpiexec -n 4 program.x    
 	



  CALCULATION OF THE SLATER-DETERMINANT FOR NEON AND BERYLLIUM
  WITH JASTROW-FACTOR
  To change between Be and Ne, just change the variable: no_of_particles.

  Variational Monte-Carlo with importance sampling and conjugate 
  gradient method CGM, includes Blocking.
  To activate blocking set the global variable blocking=true.

  Additional methods wich is called by dfpmin(), using only
  alpha and beta as parameters.

  Comment out cgm() in main if you dont want it to run.
  CGM method requires 1e7 cycles.

  If you change number of particles remember to update
  the delete_matrix(#,number of particles).

  Create an atom object f.ex atom helium.
 
  For Beryllium I fix electron 1,2 to have spin up,
  and 3,4 to have spin down.
  For Neon 1 to 5 have spin up, the rest spin down.

  I use analytical expressions for the gradient and 
  laplacian for 1s,2s and 2p.states.

  To use CGM I need the wavefunction, the whole (brute force) SD,
  when differentiating wrt alpha, beta. I have added an extra check
  for negative beta, which returns(if true) a fixed derivative and energy.


  JASTROW-FACTOR
  Create a upper triangular matrix with distances between particles.
  Takes into account the factor a which depends on two and two particles
  have opposit or equal spin-direction.
  This is general enough to handle any number of particles.

*/

#include <mpi.h>//parallellisering
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "lib.h"
#include "methods.hpp"
using namespace std;

//DEFINE GLOBAL VARIABLES
int numprocs;
int myrank;
long idum;
bool blocking = false; //activate blocking
double alphaalpha,betabeta;
bool deriv_bool;
double deriv_beta;



//separate filenames for blocking
ofstream blockofile;


ofstream outfile;


class atom{
public:
  int dimension, no_of_particles, charge;
  double **r_old, **r_new,**cumulative_E,**cumulative_E2, **mean_r_12,
    **qm_force_old, **qm_force_new;
  double* delta_E_array;
  int no_of_cycles, cycle,thermalization,max_variations,
    variate_alpha, variate_beta, accept;
  double step_length;//////////****step length and variational parameters
  
//SPIN MATRICES
  //SLATER DETERMINANT MATRICES AND VARIABLES
  double **D_up, **D_down, **temp_up, **temp_down;



  double energy, energy2, main_e, main_e2, variance;

  //Jastrow
  double **distance_old, **distance_new;

  //Roothan-Hartree

  double **root_helium, **root_beryllium, **root_neon, **root_magnesium, **root_silisium;
  double *rooth_factors;


  static const double h = 0.001;
  static const double h2 = 1e6;
  static const double D = 0.5;//Diffusion constant
  static const double time_step = 0.01;//time step

  atom();
  ~atom();
  void initialise();
  void mc_sampling(double,double,bool);
  double wave_function(double**,double, double);
  void output(int);
  
  //Slater
  double phi(double**r,double alpha,double beta, int m, int n);
  void set_up_slater(double **D_up, double **D_down, double, double);
  void update(double **D_up, double **D_down, int n, double ratio,double,double);
  double getRatio(double **D_up, double **D_down, int i, double alpha, double beta);
  void quantum_force(double**r,double**qm_force, double **D_up, double **D_down, double alpha, double beta);
  double local_energy(double **r,double **D_up, double **D_down, double alpha, double beta);
  double laplace_phi(double**r, double alpha,double beta, int n, int m);
  double gradient_phi(double**r,double alpha, int, int n, int m);

  //Jastrow
  double getJastrowRatio(double** r_old, double** r_new, double beta);
  void getDistance(double** distance, double** r_old, double beta);
  double diffR(double** r, int k, int l);
  double gradient_jastrow(double**r,int,int,int,double beta);
  //neon  
  void set_up_slater2(double**r ,double **D_up, double **D_down, double alpha, double beta);
  double det(double** A, int dim);

  //Roothan
  void roothan_initialise();
  double roothan_factor(int n_p, double chi);



  //sjekk
  void quantum_force2(double** r, double **qm_force, double wf_old,double alpha, double beta);
  double local_energy2(double **r, double wf_old,double alpha, double beta);
  /*************    SLATER DETERMINANT 1-PARTICLE WAVEFUNCTIONS    ******************************/

  //BERYLLIUM
  // inline function for 1s single-particle wave function
  inline double psi1s(double rad, double alpha) { 
    return exp(-alpha*rad);    
  }
  // inline function for 2s single-particle wave function
  inline double psi2s(double rad, double alpha) { 
    return rad*exp(-alpha*rad);
  }
  //derivative of 1s single-particle wave function
  
  inline double gradient_psi1s(double rad, double x_i, double alpha) { 
    return -alpha*x_i*exp(-alpha*rad)/rad;    
  }
  //derivative of 2s single-particle wave function
  inline double gradient_psi2s(double rad, double x_i, double alpha) { 
    return x_i*exp(-alpha*rad)*(1-alpha*rad)/rad;
  }
  //Laplace of 1s single-particle wave function
  inline double laplace_psi1s(double rad, double alpha) { 
    return alpha*exp(-alpha*rad)*(alpha*rad-2)/rad;    
  }
  //Laplace of 2s single-particle wave function
  inline double laplace_psi2s(double rad, double alpha) { 
    return exp(-rad*alpha)*(2.0/rad + alpha*alpha*rad - 4*alpha); 
  }
  //NEON
  //inline function for 2p single-particle wave function
  inline double psi2p(double rad, double alpha, double x_i) { 
    return x_i*exp(-alpha*rad);
  }
  //derivative of 2p single-particle wave function
  inline double gradient_psi2p_x(double rad, double alpha, double x_i){
    return exp(-alpha*rad)*(1-x_i*x_i*alpha/rad);
  }
  inline double gradient_psi2p_y(double rad, double alpha, double x_i, double y_i){
    return -x_i*y_i*alpha*exp(-alpha*rad)/rad;
  }
  inline double laplace_psi2p_x(double rad, double alpha, double x_i){
    return alpha*x_i*exp(-rad*alpha)*(-3 + x_i*x_i/rad/rad + x_i*x_i*alpha/rad)/rad;
  }
  inline double laplace_psi2p_y(double rad, double alpha, double x_i, double y_i){
    return x_i*alpha*exp(-alpha*rad)*(y_i*y_i/rad/rad - 1 + alpha*y_i*y_i/rad)/rad;
  }

  //MAGNESIUM
  inline double psi3s(double rad, double alpha){
    return rad*rad*exp(-alpha*rad);
  }	
  inline double gradient_psi3s(double rad, double x_i, double alpha){
    return exp(-alpha*rad)*(2*x_i - rad*alpha*x_i);
  }	
  inline double laplace_psi3s(double rad, double alpha){
    return exp(-alpha*rad)*(6 - 6*rad*alpha + alpha*alpha*rad*rad);	
  }	

  //SILISIUM

   inline double psi4p(double rad, double alpha, double x_i){
    return rad*rad*x_i*exp(-rad*alpha);
  }


  inline double gradient_psi4p_x(double rad, double alpha, double x_i){
    return (2*x_i*x_i + rad*rad - rad*x_i*x_i*alpha)*exp(-rad*alpha);
  }
 
  inline double gradient_psi4p_y(double rad, double alpha, double x_i, double y_i){
    return (2*y_i*x_i - rad*x_i*y_i*alpha)*exp(-rad*alpha);
  }

  inline double laplace_psi4p_x(double rad, double alpha, double x_i){
    return (6 - 3*x_i*x_i*alpha/rad - 3*rad*alpha + x_i*x_i*alpha*alpha)*exp(-rad*alpha)*x_i;
  }

  inline double laplace_psi4p_y(double rad, double alpha, double x_i, double y_i){
    return x_i*exp(-rad*alpha)*(2 - 3*y_i*y_i*alpha/rad - rad*alpha + alpha*alpha*y_i*y_i);
  }  //*/

};//end class atom

//constructor
atom::atom(){
  initialise();
  r_old = create_matrix(no_of_particles, dimension);  
  r_new = create_matrix(no_of_particles, dimension);  
  qm_force_old = create_matrix(no_of_particles, dimension);  
  qm_force_new = create_matrix(no_of_particles, dimension);  
  cumulative_E = create_matrix(max_variations+1,max_variations+1);
  cumulative_E2 = create_matrix(max_variations+1,max_variations+1);
  mean_r_12 = create_matrix(max_variations+1,max_variations+1);
  
  D_up = create_matrix(no_of_particles/2, no_of_particles/2);
  D_down = create_matrix(no_of_particles/2, no_of_particles/2);
  temp_up = create_matrix(no_of_particles/2, no_of_particles/2);
  temp_down = create_matrix(no_of_particles/2, no_of_particles/2);

  //Jastrow
  distance_old = create_matrix(no_of_particles,no_of_particles);  
  distance_new = create_matrix(no_of_particles,no_of_particles);  

  //Roothan-Hartree
  root_helium = create_matrix(5,2);
  root_beryllium = create_matrix(6,3);
  root_neon = create_matrix(6,5);
  root_magnesium = create_matrix(8,6);
  root_silisium = create_matrix(8,7);

  rooth_factors = new double[16];

}//end constructor

atom::~atom(){
  delete_matrix(D_up,no_of_particles/2);
  delete_matrix(D_down,no_of_particles/2);
  delete_matrix(temp_up,no_of_particles/2);
  delete_matrix(temp_down,no_of_particles/2);
  delete_matrix(root_helium, 5);
  delete_matrix(root_beryllium, 6);
  delete_matrix(root_neon, 6);
  delete_matrix(root_magnesium, 8);
  delete_matrix(root_silisium, 8);

 delete[] rooth_factors;
}
//Global helium object access from all methods
atom helium;

void atom::initialise(){

  no_of_particles=12;
  
  dimension=3;
  charge=no_of_particles;
  max_variations =0;
  step_length = 1.1;//This gives accept = 0.5
  no_of_cycles = 1.0e8;
  thermalization = 1.0e6;//0.2*no_of_cycles;
}// initialise()

void atom::mc_sampling(double alpha, double beta, bool deriv_bool){



	roothan_initialise();

	

  if(myrank==0){cout<<"\n/*Monte Carlo med alpha= "<< alpha << " og beta= " << beta <<endl;/**********************************/}
  //local variables
  double delta_E, green_ratio;
  
  

  double ratio=0;

 
  

  //JASTROW
  double jastrowRatio=0;

  //local variables for CGM
  double wf_ratio_beta, d_dbeta,exp_wf_energy_beta;
  wf_ratio_beta=0;
  deriv_beta=0;
  exp_wf_energy_beta=0;

  if(blocking==true){delta_E_array = new double[no_of_cycles];}
  energy=energy2=0; accept=0; delta_E=0;
 
  //random initial position, run over one particle at the time
  for(int i=0; i<no_of_particles; i++){
    for(int j=0; j<dimension; j++){
      r_old[i][j]=gaussian_deviate(&idum);
    }
  }
  
  //Initializng roothan-functions
  

  //SETTING UP THE SLATER MATRICES
  set_up_slater(D_up, D_down, alpha, beta);

  //SETTING UP THE DISTANCE MATRIX 
  getDistance(distance_old, r_old, beta);

  //temporary SD for the new quantum force
  set_up_slater(temp_up, temp_down, alpha, beta);

  //QUANTUM FORCE
  quantum_force(r_old, qm_force_old, D_up, D_down, alpha, beta);
 
  //LOOP OVER MONTE CARLO CYCLES
  for(cycle=0; cycle<no_of_cycles+thermalization; cycle++){

    //New position moving particle i
    for(int i=0; i<no_of_particles; i++){
      for(int j=0; j<dimension; j++){
	r_new[i][j]=r_old[i][j]+D*time_step*qm_force_old[i][j]
	  +gaussian_deviate(&idum)*sqrt(time_step); //y=x+D*timestep*F(x)+xsi  //*/

         
      }
      for (int k = 0; k < no_of_particles; k++) {
	if ( k != i) {
	  for ( int j=0; j < dimension; j++) {
	    r_new[k][j] = r_old[k][j];
	  }
	} 
      } 
      
      //SETTING UP THE DISTANCE MATRIX 
      getDistance(distance_new, r_new, beta);
      
      //SLATER DETERMINANT RATIO
      ratio = getRatio(D_up, D_down, i, alpha, beta);

      jastrowRatio = getJastrowRatio(distance_old, distance_new, beta);
      jastrowRatio = exp(jastrowRatio);
      
      //updating matrices
      for(int m=0; m<no_of_particles/2; m++){
	for(int n=0; n<no_of_particles/2; n++){
	 temp_up[m][n]= D_up[m][n];
	 temp_down[m][n]=D_down[m][n];
	}
      }
      //temporary update for quantum force
      update(temp_up, temp_down, i, ratio, alpha, beta);
      quantum_force(r_new, qm_force_new, temp_up, temp_down, alpha, beta);
      
      //IMPORTANCE SAMPLING
      double  greensfunction = 0.0;            
      for(int k=0; k<no_of_particles; k++){
	for(int j=0; j < dimension; j++) {
	  greensfunction += 0.5*(qm_force_old[k][j]+qm_force_new[k][j])*
	    (D*time_step*0.5*(qm_force_old[k][j]-qm_force_new[k][j])-r_new[k][j]+r_old[k][j]);
	}
      }
      greensfunction = exp(greensfunction);

      
///

      //METROPOLIS TEST
      if(ran1(&idum)<= greensfunction*ratio*ratio*jastrowRatio*jastrowRatio){
	for(int k=0; k<no_of_particles; k++){
	  for(int j=0; j<dimension; j++){
	    r_old[k][j]=r_new[k][j];
	    qm_force_old[k][j] = qm_force_new[k][j];
	  }
	}
	accept +=1;

	//UPDATE INVERSE SLATER MATRICES
	for(int m=0; m<no_of_particles/2; m++){
	  for(int n=0; n<no_of_particles/2; n++){
	    D_up[m][n]=temp_up[m][n];
	    D_down[m][n]=temp_down[m][n];
	  }
	}
	//UPDATE THE DISTANCE MATRIX
	for(int m=0; m<no_of_particles; m++){
	  for(int n=0; n<no_of_particles; n++){
	    distance_old[m][n]=distance_new[m][n];
	  }
	}
	
      }//end metropolis test
          
    }//end loop over particles

	

    //LOCAL ENERGY
    if(cycle>thermalization){
      
      delta_E = local_energy(r_old, D_up, D_down, alpha, beta);
      energy +=delta_E;
      energy2 += delta_E*delta_E;



      if(blocking==true){delta_E_array[cycle-thermalization]=delta_E;}
      
      //CGM DERIVATIVES USES BRUTE FORCE WAVEFUNCTIONS
      //compute derivative for dE_energy()
      if(deriv_bool==true){
	//d/d(alpha)
		
	//d/d(beta)
	d_dbeta  = (wave_function(r_new,alpha,beta+h)-
		    wave_function(r_new,alpha,beta-h))/(2*h); 
	wf_ratio_beta += d_dbeta/wave_function(r_new,alpha,beta);
		
	exp_wf_energy_beta +=(d_dbeta/wave_function(r_new,alpha,beta))*delta_E;
      }
    }//end energy sampling
  }//end loop over Monte Carlo cycles
  
  if(blocking==true){
    //RESULTS writing to separate file on each node, binary
    ostringstream f;
    f<<"blocks_rank"<<myrank<<".dat";
    blockofile.open(f.str().c_str(), ios::out | ios::binary);
    
    blockofile.write((char*)(delta_E_array+1), 
		     no_of_cycles*sizeof(double));
    
    blockofile.close();
    delete[] delta_E_array;
  }
  
  energy=energy/no_of_cycles;
  energy2=energy2/no_of_cycles;


  if(deriv_bool==true){
    
    wf_ratio_beta = wf_ratio_beta/no_of_cycles;

    deriv_beta  = 2*exp_wf_energy_beta/no_of_cycles  - 2*wf_ratio_beta*energy;
  }
  
  cout<<"Myrank ="<<myrank<<" MC(): "<<"mcs  energy = "<<energy<< " alpha = "<<alpha<<" beta = "<<beta<<endl;/*****************************/

  //The program stops for after one vmc round if blocking
  if(blocking==true){
    MPI_Barrier(MPI_COMM_WORLD);
  }
  //free memory


  MPI_Reduce(&energy, &main_e, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&energy2, &main_e2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


}//end mc_sampling()

/*************************    SLATER DETERMINANT WAVEFUNCTION    ********************
 *
 * Input: Position-matrix r
 *        Particle m 
 *        State n
 *        +variational parameters
 *
 *
 */

void atom::roothan_initialise(){

	if(charge==2){

		//Exponent in S-wave-functions
		root_helium[0][0]=1.41714;	///1s
		root_helium[1][0]=2.37682;	///1s
		root_helium[2][0]=4.39628;	///1s
		root_helium[3][0]=6.52699;	///1s
		root_helium[4][0]=7.94252;	///1s
		
		//1s expansion coefficients
		root_helium[0][1]=0.76838;
		root_helium[1][1]=0.22346;
		root_helium[2][1]=0.04082;
		root_helium[3][1]=-0.00994;
		root_helium[4][1]=0.00230;
	
		rooth_factors[0]=roothan_factor(1,root_helium[0][0]);
		rooth_factors[1]=roothan_factor(1,root_helium[1][0]);
		rooth_factors[2]=roothan_factor(1,root_helium[2][0]);
		rooth_factors[3]=roothan_factor(1,root_helium[3][0]);
		rooth_factors[4]=roothan_factor(1,root_helium[4][0]);

	}

	

	if(charge==4){

		//Exponent in S-wave-functions
		


		root_beryllium[0][0]=3.47116;	///1s
		root_beryllium[1][0]=6.36861;	///1s
		root_beryllium[2][0]=0.77820;	///2s
		root_beryllium[3][0]=0.94067;	///2s
		root_beryllium[4][0]=1.48725;	///2s
		root_beryllium[5][0]=2.71830;	///2s
		
		//1s expansion coefficients
		root_beryllium[0][1]=0.91796;
		root_beryllium[1][1]=0.08724;
		root_beryllium[2][1]=0.00108;
		root_beryllium[3][1]=-0.00199;
		root_beryllium[4][1]=0.00176;
		root_beryllium[5][1]=0.00628;
		
		//2s expansion coefficients
		root_beryllium[0][2]=-0.17092;
		root_beryllium[1][2]=-0.01455;
		root_beryllium[2][2]=0.21186;
		root_beryllium[3][2]=0.62499;
		root_beryllium[4][2]=0.26662;
		root_beryllium[5][2]=-0.09919;

		rooth_factors[0]=roothan_factor(1,root_beryllium[0][0]);
		rooth_factors[1]=roothan_factor(1,root_beryllium[1][0]);
		rooth_factors[2]=roothan_factor(2,root_beryllium[2][0]);
		rooth_factors[3]=roothan_factor(2,root_beryllium[3][0]);
		rooth_factors[4]=roothan_factor(2,root_beryllium[4][0]);
		rooth_factors[5]=roothan_factor(2,root_beryllium[5][0]);
		


	}
	
	if(charge==10){

		//Exponent in S-wave-functions
		root_neon[0][0]=9.48486;	///1s
		root_neon[1][0]=15.56590;	///1s
		root_neon[2][0]=1.96184;	///2s
		root_neon[3][0]=2.86423;	///2s
		root_neon[4][0]=4.82530;	///2s
		root_neon[5][0]=7.79242;	///2s
	
		//1s expansion coefficients
		root_neon[0][1]=0.93717;
		root_neon[1][1]=0.04899;
		root_neon[2][1]=0.00058;
		root_neon[3][1]=-0.00064;
		root_neon[4][1]=0.00551;
		root_neon[5][1]=0.01999;
	
		//2s expansion coefficients
		root_neon[0][2]=-0.23093;
		root_neon[1][2]=-0.00635;
		root_neon[2][2]=0.18620;
		root_neon[3][2]=0.66899;
		root_neon[4][2]=0.30910;
		root_neon[5][2]=-0.13871;
	
		rooth_factors[0]=roothan_factor(1,root_neon[0][0]);
		rooth_factors[1]=roothan_factor(1,root_neon[1][0]);
		rooth_factors[2]=roothan_factor(2,root_neon[2][0]);
		rooth_factors[3]=roothan_factor(2,root_neon[3][0]);
		rooth_factors[4]=roothan_factor(2,root_neon[4][0]);
		rooth_factors[5]=roothan_factor(2,root_neon[5][0]);

		//Exponent in P-wave-functions
		root_neon[0][3]=1.45208;	///2p
		root_neon[1][3]=2.38168;	///2p
		root_neon[2][3]=4.48489;	///2p
		root_neon[3][3]=9.13464;	///2p
	
		//2p expansion coefficients
		root_neon[0][4]=0.21799;
		root_neon[1][4]=0.53338;
		root_neon[2][4]=0.32933;
		root_neon[3][4]=0.01872;

		rooth_factors[6]=roothan_factor(2,root_neon[0][3]);
		rooth_factors[7]=roothan_factor(2,root_neon[1][3]);
		rooth_factors[8]=roothan_factor(2,root_neon[2][3]);
		rooth_factors[9]=roothan_factor(2,root_neon[3][3]);


	}

	if(charge==12){
	
		//Exponent in S-wave-functions
		root_magnesium[0][0]=12.01140;	///1s
		root_magnesium[1][0]=13.91620;	///3s
		root_magnesium[2][0]=9.48612;	///3s
		root_magnesium[3][0]=6.72188;	///3s
		root_magnesium[4][0]=4.24466;	///3s
		root_magnesium[5][0]=2.53466;	///3s
		root_magnesium[6][0]=1.46920;	///3s
		root_magnesium[7][0]=0.89084;	///3s
	
		//1s expansion coefficients
		root_magnesium[0][1]=0.96430;
		root_magnesium[1][1]=0.03548;
		root_magnesium[2][1]=0.02033;
		root_magnesium[3][1]=-0.00252;
		root_magnesium[4][1]=0.00162;
		root_magnesium[5][1]=-0.00038;
		root_magnesium[6][1]=0.00015;
		root_magnesium[7][1]=-0.00004;
		
		//2s expansion coefficients
		root_magnesium[0][2]=-0.24357;
		root_magnesium[1][2]=-0.00485;
		root_magnesium[2][2]=0.08002;
		root_magnesium[3][2]=0.39902;
		root_magnesium[4][2]=0.57358;
		root_magnesium[5][2]=0.05156;
		root_magnesium[6][2]=-0.00703;
		root_magnesium[7][2]=0.00161;
	
		//3s expansion coefficients
		root_magnesium[0][3]=0.04691;
		root_magnesium[1][3]=0.00144;
		root_magnesium[2][3]=-0.01850;
		root_magnesium[3][3]=-0.07964;
		root_magnesium[4][3]=-0.13478;
		root_magnesium[5][3]=-0.01906;
		root_magnesium[6][3]=0.48239;
		root_magnesium[7][3]=0.60221;

		rooth_factors[0]=roothan_factor(1,root_magnesium[0][0]);
		rooth_factors[1]=roothan_factor(3,root_magnesium[1][0]);
		rooth_factors[2]=roothan_factor(3,root_magnesium[2][0]);
		rooth_factors[3]=roothan_factor(3,root_magnesium[3][0]);
		rooth_factors[4]=roothan_factor(3,root_magnesium[4][0]);
		rooth_factors[5]=roothan_factor(3,root_magnesium[5][0]);
		rooth_factors[6]=roothan_factor(3,root_magnesium[6][0]);
		rooth_factors[7]=roothan_factor(3,root_magnesium[7][0]);
	



		//Exponent in P-wave-functions
		root_magnesium[0][4]=5.92580;	///2p
		root_magnesium[1][4]=7.98979;	///4p
		root_magnesium[2][4]=5.32964;	///4p
		root_magnesium[3][4]=3.71678;	///4p
		root_magnesium[4][4]=2.59986;	///4p
		
		//2p expansion coefficients
		root_magnesium[0][5]=0.52391;
		root_magnesium[1][5]=0.07012;
		root_magnesium[2][5]=0.31965;
		root_magnesium[3][5]=0.20860;
		root_magnesium[4][5]=0.03888;

		rooth_factors[8]=roothan_factor(2,root_magnesium[0][4]);
		rooth_factors[9]=roothan_factor(4,root_magnesium[1][4]);
		rooth_factors[10]=roothan_factor(4,root_magnesium[2][4]);
		rooth_factors[11]=roothan_factor(4,root_magnesium[3][4]);
		rooth_factors[12]=roothan_factor(4,root_magnesium[4][4]);



	}

        if(charge==14){

		//Exponent in S-wave-functions
		root_silisium[0][0]=14.01420;	///1S
		root_silisium[1][0]=16.39320;	///3S
		root_silisium[2][0]=10.87950;	///3S
		root_silisium[3][0]=7.72709;	///3S
		root_silisium[4][0]=5.16500;	///3S
		root_silisium[5][0]=2.97451;	///3S
		root_silisium[6][0]=2.14316;	///3S
		root_silisium[7][0]=1.09108;	///3S

		//1s expansion coefficients
		root_silisium[0][1]=0.96800;
		root_silisium[1][1]=0.03033;
		root_silisium[2][1]=0.02248;
		root_silisium[3][1]=-0.00617;
		root_silisium[4][1]=0.00326;
		root_silisium[5][1]=-0.00143;
		root_silisium[6][1]=0.00081;
		root_silisium[7][1]=-0.00016;
		
		//2s expansion coefficients
		root_silisium[0][2]=-0.25755;
		root_silisium[1][2]=-0.00446;
		root_silisium[2][2]=0.11153;
		root_silisium[3][2]=0.40339;
		root_silisium[4][2]=0.55032;
		root_silisium[5][2]=0.03381;
		root_silisium[6][2]=-0.00815;
		root_silisium[7][2]=0.00126;

		//3s expansion coefficients
		root_silisium[0][3]=0.06595;
		root_silisium[1][3]=0.00185;
		root_silisium[2][3]=-0.03461;
		root_silisium[3][3]=-0.10378;
		root_silisium[4][3]=-0.19229;
		root_silisium[5][3]=-0.06561;
		root_silisium[6][3]=0.59732;
		root_silisium[7][3]=0.55390;

		rooth_factors[0]=roothan_factor(1,root_silisium[0][0]);
		rooth_factors[1]=roothan_factor(3,root_silisium[1][0]);
		rooth_factors[2]=roothan_factor(3,root_silisium[2][0]);
		rooth_factors[3]=roothan_factor(3,root_silisium[3][0]);
		rooth_factors[4]=roothan_factor(3,root_silisium[4][0]);
		rooth_factors[5]=roothan_factor(3,root_silisium[5][0]);
		rooth_factors[6]=roothan_factor(3,root_silisium[6][0]);
		rooth_factors[7]=roothan_factor(3,root_silisium[7][0]);



		//Exponent in P-wave-functions
		root_silisium[0][4]=7.14360;	///2P
		root_silisium[1][4]=16.25720;	///4P
		root_silisium[2][4]=10.79720;	///4P
		root_silisium[3][4]=6.89724;	///4P
		root_silisium[4][4]=4.66598;	///4P
		root_silisium[5][4]=2.32046;	///4P
		root_silisium[6][4]=1.33470;	///4P
		root_silisium[7][4]=0.79318;	///4P

		//2p expansion coefficients
		root_silisium[0][5]=0.54290;
		root_silisium[1][5]=0.00234;
		root_silisium[2][5]=0.04228;
		root_silisium[3][5]=0.32155;
		root_silisium[4][5]=0.22474;
		root_silisium[5][5]=0.00732;
		root_silisium[6][5]=-0.00105;
		root_silisium[7][5]=0.00041;

		//3p expansion coefficients
		root_silisium[0][6]=-0.11535;
		root_silisium[1][6]=-0.00189;
		root_silisium[2][6]=-0.00473;
		root_silisium[3][6]=-0.07552;
		root_silisium[4][6]=0.01041;
		root_silisium[5][6]=0.46075;
		root_silisium[6][6]=0.57665;
		root_silisium[7][6]=0.06274;

		rooth_factors[8]=roothan_factor(2,root_silisium[0][4]);
		rooth_factors[9]=roothan_factor(4,root_silisium[1][4]);
		rooth_factors[10]=roothan_factor(4,root_silisium[2][4]);
		rooth_factors[11]=roothan_factor(4,root_silisium[3][4]);
		rooth_factors[12]=roothan_factor(4,root_silisium[4][4]);
		rooth_factors[13]=roothan_factor(4,root_silisium[5][4]);
		rooth_factors[14]=roothan_factor(4,root_silisium[6][4]);
		rooth_factors[15]=roothan_factor(4,root_silisium[7][4]);
		

	}
}

double atom::roothan_factor(int n_p, double chi){

  double factor1,factor2;
  factor1=1;
  for(int i=2*n_p; i>1; i--){
	factor1=factor1*i;
  }
  factor1=sqrt(factor1);
  factor1=1.0/factor1;

  factor2=2*chi;
  factor2=pow(factor2,(double)(n_p+0.5));


  return factor1*factor2;
}
double atom::phi(double**r,double alpha,double beta, int n, int m){



  double r_single_particle;
  double *argument1,*argument2;
  double temp =0;
  
  argument1 = new double[no_of_particles];
  argument2 = new double[no_of_particles];
  //zero arguments
  for (int i=0; i<no_of_particles; i++){ 
    argument1[i] = 0.0;
    argument2[i] = 0.0;
  }   
  
  for(int i=0; i<no_of_particles; i++){
    r_single_particle=0;
    for(int j=0; j<dimension; j++){
      r_single_particle += r[i][j]*r[i][j];//r^2=x^2+y^2+z^2
    }
    argument1[i] += sqrt(r_single_particle);
  }
  //Check which state 			
  if(n==0){
	//Check which atom
	


    	if(charge==2){
		temp = root_helium[0][1]*rooth_factors[0]*psi1s(argument1[m], root_helium[0][0]) + root_helium[1][1]*rooth_factors[1]*psi1s(argument1[m], root_helium[1][0]) + root_helium[2][1]*rooth_factors[2]*psi1s(argument1[m], root_helium[2][0]) + root_helium[3][1]*rooth_factors[3]*psi1s(argument1[m], root_helium[3][0]) + root_helium[4][1]*rooth_factors[4]*psi1s(argument1[m], root_helium[4][0]);
		}

	

	if(charge==4){
	
		temp = root_beryllium[0][1]*rooth_factors[0]*psi1s(argument1[m], root_beryllium[0][0]) + root_beryllium[1][1]*rooth_factors[1]*psi1s(argument1[m], root_beryllium[1][0]) + root_beryllium[2][1]*rooth_factors[2]*psi2s(argument1[m], root_beryllium[2][0]) + root_beryllium[3][1]*rooth_factors[3]*psi2s(argument1[m], root_beryllium[3][0]) +  root_beryllium[4][1]*rooth_factors[4]*psi2s(argument1[m], root_beryllium[4][0]) + root_beryllium[5][1]*rooth_factors[5]*psi2s(argument1[m], root_beryllium[5][0]);
		}

	


	if(charge==10){
		temp = root_neon[0][1]*rooth_factors[0]*psi1s(argument1[m], root_neon[0][0]) + root_neon[1][1]*rooth_factors[1]*psi1s(argument1[m], root_neon[1][0]) + root_neon[2][1]*rooth_factors[2]*psi2s(argument1[m], root_neon[2][0]) + root_neon[3][1]*rooth_factors[3]*psi2s(argument1[m], root_neon[3][0]) + root_neon[4][1]*rooth_factors[4]*psi2s(argument1[m], root_neon[4][0]) + root_neon[5][1]*rooth_factors[5]*psi2s(argument1[m], root_neon[5][0]);
		}
	if(charge==12){
		temp = root_magnesium[0][1]*rooth_factors[0]*psi1s(argument1[m], root_magnesium[0][0]) + 
		root_magnesium[1][1]*rooth_factors[1]*psi3s(argument1[m], root_magnesium[1][0]) + 
		root_magnesium[2][1]*rooth_factors[2]*psi3s(argument1[m], root_magnesium[2][0]) + 
		root_magnesium[3][1]*rooth_factors[3]*psi3s(argument1[m], root_magnesium[3][0]) + 
		root_magnesium[4][1]*rooth_factors[4]*psi3s(argument1[m], root_magnesium[4][0]) + 
		root_magnesium[5][1]*rooth_factors[5]*psi3s(argument1[m], root_magnesium[5][0]) + 
		root_magnesium[6][1]*rooth_factors[6]*psi3s(argument1[m], root_magnesium[6][0]) + 
		root_magnesium[7][1]*rooth_factors[7]*psi3s(argument1[m], root_magnesium[7][0]);
		}
	if(charge==14){
		temp = root_silisium[0][1]*rooth_factors[0]*psi1s(argument1[m], root_silisium[0][0]) + root_silisium[1][1]*rooth_factors[1]*psi3s(argument1[m], root_silisium[1][0]) + root_silisium[2][1]*rooth_factors[2]*psi3s(argument1[m], root_silisium[2][0]) + root_silisium[3][1]*rooth_factors[3]*psi3s(argument1[m], root_silisium[3][0]) + root_silisium[4][1]*rooth_factors[4]*psi3s(argument1[m], root_silisium[4][0]) + root_silisium[5][1]*rooth_factors[5]*psi3s(argument1[m], root_silisium[5][0]) + root_silisium[6][1]*rooth_factors[6]*psi3s(argument1[m], root_silisium[6][0]) + root_silisium[7][1]*rooth_factors[7]*psi3s(argument1[m], root_silisium[7][0]);

		}
  } 
  if(n==1){
  	if(charge==4){
		temp = root_beryllium[0][2]*rooth_factors[0]*psi1s(argument1[m], root_beryllium[0][0]) + root_beryllium[1][2]*rooth_factors[1]*psi1s(argument1[m], root_beryllium[1][0]) + root_beryllium[2][2]*rooth_factors[2]*psi2s(argument1[m], root_beryllium[2][0]) + 
		root_beryllium[3][2]*rooth_factors[3]*psi2s(argument1[m], root_beryllium[3][0]) +  root_beryllium[4][2]*rooth_factors[4]*psi2s(argument1[m], root_beryllium[4][0]) +  root_beryllium[5][2]*rooth_factors[5]*psi2s(argument1[m], root_beryllium[5][0]);
		}
	if(charge==10){
		temp = root_neon[0][2]*rooth_factors[0]*psi1s(argument1[m], root_neon[0][0]) + root_neon[1][2]*rooth_factors[1]*psi1s(argument1[m], root_neon[1][0]) + root_neon[2][2]*rooth_factors[2]*psi2s(argument1[m], root_neon[2][0]) + root_neon[3][2]*rooth_factors[3]*psi2s(argument1[m], root_neon[3][0]) + root_neon[4][2]*rooth_factors[4]*psi2s(argument1[m], root_neon[4][0]) + root_neon[5][2]*rooth_factors[5]*psi2s(argument1[m], root_neon[5][0]);
		} 
	if(charge==12){
		temp = 	root_magnesium[0][2]*rooth_factors[0]*psi1s(argument1[m], root_magnesium[0][0]) + root_magnesium[1][2]*rooth_factors[1]*psi3s(argument1[m], root_magnesium[1][0]) + root_magnesium[2][2]*rooth_factors[2]*psi3s(argument1[m], root_magnesium[2][0]) + root_magnesium[3][2]*rooth_factors[3]*psi3s(argument1[m], root_magnesium[3][0]) + root_magnesium[4][2]*rooth_factors[4]*psi3s(argument1[m], root_magnesium[4][0]) + root_magnesium[5][2]*rooth_factors[5]*psi3s(argument1[m], root_magnesium[5][0]) + root_magnesium[6][2]*rooth_factors[6]*psi3s(argument1[m], root_magnesium[6][0]) + root_magnesium[7][2]*rooth_factors[7]*psi3s(argument1[m], root_magnesium[7][0]) ;
		}
	if(charge==14){
		temp = root_silisium[0][2]*rooth_factors[0]*psi1s(argument1[m], root_silisium[0][0]) + root_silisium[1][2]*rooth_factors[1]*psi3s(argument1[m], root_silisium[1][0]) + root_silisium[2][2]*rooth_factors[2]*psi3s(argument1[m], root_silisium[2][0]) + root_silisium[3][2]*rooth_factors[3]*psi3s(argument1[m], root_silisium[3][0]) + root_silisium[4][2]*rooth_factors[4]*psi3s(argument1[m], root_silisium[4][0]) + root_silisium[5][2]*rooth_factors[5]*psi3s(argument1[m], root_silisium[5][0]) + root_silisium[6][2]*rooth_factors[6]*psi3s(argument1[m], root_silisium[6][0]) + root_silisium[7][2]*rooth_factors[7]*psi3s(argument1[m], root_silisium[7][0]);
 

		}
  } 
  if(n==2){
	if(charge==10){
		temp= root_neon[0][4]*rooth_factors[6]*psi2p(argument1[m], root_neon[0][3], r[m][0]) + root_neon[1][4]*rooth_factors[7]*psi2p(argument1[m], root_neon[1][3], r[m][0]) +  root_neon[2][4]*rooth_factors[8]*psi2p(argument1[m], root_neon[2][3], r[m][0]) + 
		root_neon[3][4]*rooth_factors[9]*psi2p(argument1[m], root_neon[3][3], r[m][0]);
	}
	if(charge==12){
		temp = root_magnesium[0][5]*rooth_factors[8]*psi2p(argument1[m], root_magnesium[0][4], r[m][0]) + root_magnesium[1][5]*rooth_factors[9]*psi4p(argument1[m], root_magnesium[1][4], r[m][0]) + root_magnesium[2][5]*rooth_factors[10]*psi4p(argument1[m], root_magnesium[2][4], r[m][0]) + root_magnesium[3][5]*rooth_factors[11]*psi4p(argument1[m], root_magnesium[3][4], r[m][0]) + root_magnesium[4][5]*rooth_factors[12]*psi4p(argument1[m], root_magnesium[4][4], r[m][0]);
		}

	if(charge==14){
		temp = root_silisium[0][5]*rooth_factors[8]*psi2p(argument1[m], root_silisium[0][4], r[m][0]) + root_silisium[1][5]*rooth_factors[9]*psi4p(argument1[m], root_silisium[1][4], r[m][0]) + root_silisium[2][5]*rooth_factors[10]*psi4p(argument1[m], root_silisium[2][4], r[m][0]) + root_silisium[3][5]*rooth_factors[11]*psi4p(argument1[m], root_silisium[3][4], r[m][0]) + root_silisium[4][5]*rooth_factors[12]*psi4p(argument1[m], root_silisium[4][4], r[m][0]) + root_silisium[5][5]*rooth_factors[13]*psi4p(argument1[m], root_silisium[5][4], r[m][0]) + root_silisium[6][5]*rooth_factors[14]*psi4p(argument1[m], root_silisium[6][4], r[m][0]) + root_silisium[7][5]*rooth_factors[15]*psi4p(argument1[m], root_silisium[7][4], r[m][0]);
		}
  }
  if(n==3){
    if(charge==10){
		temp= root_neon[0][4]*rooth_factors[6]*psi2p(argument1[m], root_neon[0][3], r[m][1]) + root_neon[1][4]*rooth_factors[7]*psi2p(argument1[m], root_neon[1][3], r[m][1]) +  root_neon[2][4]*rooth_factors[8]*psi2p(argument1[m], root_neon[2][3], r[m][1]) + 
		root_neon[3][4]*rooth_factors[9]*psi2p(argument1[m], root_neon[3][3], r[m][1]);
	}
	if(charge==12){
		temp = root_magnesium[0][5]*rooth_factors[8]*psi2p(argument1[m], root_magnesium[0][4], r[m][1]) + root_magnesium[1][5]*rooth_factors[9]*psi4p(argument1[m], root_magnesium[1][4], r[m][1]) + root_magnesium[2][5]*rooth_factors[10]*psi4p(argument1[m], root_magnesium[2][4], r[m][1]) + root_magnesium[3][5]*rooth_factors[11]*psi4p(argument1[m], root_magnesium[3][4], r[m][1]) + root_magnesium[4][5]*rooth_factors[12]*psi4p(argument1[m], root_magnesium[4][4], r[m][1]);
		}

	if(charge==14){
		temp = root_silisium[0][5]*rooth_factors[8]*psi2p(argument1[m], root_silisium[0][4], r[m][1]) + root_silisium[1][5]*rooth_factors[9]*psi4p(argument1[m], root_silisium[1][4], r[m][1]) + root_silisium[2][5]*rooth_factors[10]*psi4p(argument1[m], root_silisium[2][4], r[m][1]) + root_silisium[3][5]*rooth_factors[11]*psi4p(argument1[m], root_silisium[3][4], r[m][1]) + root_silisium[4][5]*rooth_factors[12]*psi4p(argument1[m], root_silisium[4][4], r[m][1]) + root_silisium[5][5]*rooth_factors[13]*psi4p(argument1[m], root_silisium[5][4], r[m][1]) + root_silisium[6][5]*rooth_factors[14]*psi4p(argument1[m], root_silisium[6][4], r[m][1]) + root_silisium[7][5]*rooth_factors[15]*psi4p(argument1[m], root_silisium[7][4], r[m][1]);
		}
  }
  if(n==4){
    if(charge==10){
		temp= root_neon[0][4]*rooth_factors[6]*psi2p(argument1[m], root_neon[0][3], r[m][2]) + root_neon[1][4]*rooth_factors[7]*psi2p(argument1[m], root_neon[1][3], r[m][2]) +  root_neon[2][4]*rooth_factors[8]*psi2p(argument1[m], root_neon[2][3], r[m][2]) + 
		root_neon[3][4]*rooth_factors[9]*psi2p(argument1[m], root_neon[3][3], r[m][2]);
	}
	if(charge==12){
		temp = root_magnesium[0][5]*rooth_factors[8]*psi2p(argument1[m], root_magnesium[0][4], r[m][2]) + root_magnesium[1][5]*rooth_factors[9]*psi4p(argument1[m], root_magnesium[1][4], r[m][2]) + root_magnesium[2][5]*rooth_factors[10]*psi4p(argument1[m], root_magnesium[2][4], r[m][2]) + root_magnesium[3][5]*rooth_factors[11]*psi4p(argument1[m], root_magnesium[3][4], r[m][2]) + root_magnesium[4][5]*rooth_factors[12]*psi4p(argument1[m], root_magnesium[4][4], r[m][2]);
		}

	if(charge==14){
		temp = root_silisium[0][5]*rooth_factors[8]*psi2p(argument1[m], root_silisium[0][4], r[m][2]) + root_silisium[1][5]*rooth_factors[9]*psi4p(argument1[m], root_silisium[1][4], r[m][2]) + root_silisium[2][5]*rooth_factors[10]*psi4p(argument1[m], root_silisium[2][4], r[m][2]) + root_silisium[3][5]*rooth_factors[11]*psi4p(argument1[m], root_silisium[3][4], r[m][2]) + root_silisium[4][5]*rooth_factors[12]*psi4p(argument1[m], root_silisium[4][4], r[m][2]) + root_silisium[5][5]*rooth_factors[13]*psi4p(argument1[m], root_silisium[5][4], r[m][2]) + root_silisium[6][5]*rooth_factors[14]*psi4p(argument1[m], root_silisium[6][4], r[m][2]) + root_silisium[7][5]*rooth_factors[15]*psi4p(argument1[m], root_silisium[7][4], r[m][2]);
		}
  }
   if(n==5){
   	if(charge==12){
		temp= root_magnesium[0][3]*rooth_factors[0]*psi1s(argument1[m], root_magnesium[0][0]) + root_magnesium[1][3]*rooth_factors[1]*psi3s(argument1[m], root_magnesium[1][0]) + root_magnesium[2][3]*rooth_factors[2]*psi3s(argument1[m], root_magnesium[2][0]) + root_magnesium[3][3]*rooth_factors[3]*psi3s(argument1[m], root_magnesium[3][0]) + root_magnesium[4][3]*rooth_factors[4]*psi3s(argument1[m], root_magnesium[4][0]) + root_magnesium[5][3]*rooth_factors[5]*psi3s(argument1[m], root_magnesium[5][0]) + root_magnesium[6][3]*rooth_factors[6]*psi3s(argument1[m], root_magnesium[6][0]) + root_magnesium[7][3]*rooth_factors[7]*psi3s(argument1[m], root_magnesium[7][0]); 

		}
	if(charge==14){
		temp = root_silisium[0][3]*rooth_factors[0]*psi1s(argument1[m], root_silisium[0][0]) + root_silisium[1][3]*rooth_factors[1]*psi3s(argument1[m], root_silisium[1][0]) + root_silisium[2][3]*rooth_factors[2]*psi3s(argument1[m], root_silisium[2][0]) + root_silisium[3][3]*rooth_factors[3]*psi3s(argument1[m], root_silisium[3][0]) + root_silisium[4][3]*rooth_factors[4]*psi3s(argument1[m], root_silisium[4][0]) + root_silisium[5][3]*rooth_factors[5]*psi3s(argument1[m], root_silisium[5][0]) + root_silisium[6][3]*rooth_factors[6]*psi3s(argument1[m], root_silisium[6][0]) + root_silisium[7][3]*rooth_factors[7]*psi3s(argument1[m], root_silisium[7][0]); 

		}

  } //*/
   if(n==6){
    		temp = root_silisium[0][6]*rooth_factors[8]*psi2p(argument1[m], root_silisium[0][4], r[m][1]) +
		root_silisium[1][6]*rooth_factors[9]*psi4p(argument1[m], root_silisium[1][4], r[m][1]) +
		root_silisium[2][6]*rooth_factors[10]*psi4p(argument1[m], root_silisium[2][4], r[m][1]) +
		root_silisium[3][6]*rooth_factors[11]*psi4p(argument1[m], root_silisium[3][4], r[m][1]) +
		root_silisium[4][6]*rooth_factors[12]*psi4p(argument1[m], root_silisium[4][4], r[m][1]) +
		root_silisium[5][6]*rooth_factors[13]*psi4p(argument1[m], root_silisium[5][4], r[m][1]) +
		root_silisium[6][6]*rooth_factors[14]*psi4p(argument1[m], root_silisium[6][4], r[m][1]) +
		root_silisium[7][6]*rooth_factors[15]*psi4p(argument1[m], root_silisium[7][4], r[m][1]); 
  }
//VET IKKEOM DET SKAL VÆRE 0,1 eller 2
  delete [] argument1;
  delete [] argument2;
  return temp;
}

/**********************  SETTING UP THE SLATER MATRICES   ********************  
 *
 *  If you write it as a 4x4-matrix, the determinant will be zero.
 *  Instead write it as a product of two 2x2-matrices. One for spin up
 *  and one for spin down. Valid as long as the Hamiltonian is independent
 *  of spin. Each matrix has an additional index.
 *  Spin up   = 0
 *  Spin down = 1
 *
 *  for spin opp og ned har jeg byttet indeksene for phi(), motsatt fra det
 *  jeg tror er riktig, men dette gir riktig svar.
 *
 */
void atom::set_up_slater(double **D_up, double **D_down, double alpha, double beta){
  //SPIN UP
  for(int i=0; i<no_of_particles/2; i++){
    for(int j=0; j<no_of_particles/2; j++){
      D_up[i][j] = phi(r_old,alpha,beta,j,i);/********/
    }
  }
  //SPIN DOWN
  for(int i=0; i<no_of_particles/2; i++){
    for(int j=0; j<no_of_particles/2; j++){
      D_down[i][j]=phi(r_old,alpha,beta,j,i+no_of_particles/2); 
    }
  }
  //inverting matrices
  inverse(D_up, no_of_particles/2);
  inverse(D_down, no_of_particles/2);
}//end set_up_slater()

/********************    UPDATE INVERSE SLATER MATRICES    ************************
 *
 *  If the new position is updated we 
 *  need to update the inverse Slater Matrix.
 *  I dont use D_up/down_-new because i never use it 
 *  again after it is overwritten.
 *
 */
void atom::update(double **D_up, double **D_down, int i, double ratio, double alpha, double beta){
  double *s_j;
  s_j = new double[no_of_particles/2]; 
  
  for(int l=0; l<no_of_particles/2; l++){
    s_j[l]=0;
  }
  
  if(i<no_of_particles/2){//(i==0 || i==1){
    //i != j
    for(int j=0; j<no_of_particles/2; j++){
      if(j!=i){
	for(int l=0; l<no_of_particles/2; l++){
	  s_j[j]+=phi(r_new,alpha,beta,l,i)*D_up[l][j];
	}
      } 
    }
    for(int j=0; j<no_of_particles/2; j++){
      if(j!=i){
	for(int k=0; k<no_of_particles/2; k++){
	  D_up[k][j]=D_up[k][j]-s_j[j]*D_up[k][i]/ratio;	    
	}
      }
    }
    //i'th column, i=j
    for(int k=0; k<no_of_particles/2; k++){
      D_up[k][i]=D_up[k][i]/ratio;	    
    }
  }
  else{//if(i==2 || i==3){
    //i != j
    i=i-no_of_particles/2;//for å være konsistent med løkkene over j!=i
    for(int j=0; j<no_of_particles/2; j++){
      if(j!=i){
	for(int l=0; l<no_of_particles/2; l++){
	  s_j[j]+=phi(r_new,alpha,beta,l,i+no_of_particles/2)*D_down[l][j];
	}
      } 
    }
    for(int j=0; j<no_of_particles/2; j++){
      if(j!=i){
	for(int k=0; k<no_of_particles/2; k++){
	  D_down[k][j]=D_down[k][j]-s_j[j]*D_down[k][i]/ratio;	    
	}
      }
    }
    //i'th column, i=j
    for(int k=0; k<no_of_particles/2; k++){
      D_down[k][i]=D_down[k][i]/ratio;	    
    }
    i=i+no_of_particles/2;
  }
  
  delete [] s_j;
}//end update()

/*************************     RATIO     ********************************
 *
 *       For example, moving particle 2:
 *       ratio = det(1,2new)*det(3,4)/(det(1,2)*det(3,4))
 *             = det(1,2new)/det(1,2)
 *
 *       Particle i has been moved
 *
 */
double atom::getRatio(double **D_up, double **D_down, int i, double alpha, double beta){
  double ratio=0;
  if(i<no_of_particles/2){//(i==0 || i==1){
    for(int j=0; j<no_of_particles/2; j++){
      ratio += phi(r_new,alpha,beta,j,i)*D_up[j][i];
    }
    return ratio;
  }
  else{//if(i==2 || i==3){
    for(int j=0; j<no_of_particles/2; j++){
      ratio += phi(r_new,alpha,beta,j,i)*D_down[j][i-no_of_particles/2];
    }
    return ratio;
  }
}//end getRatio()

/*******************     QUANTUM FORCE SLATER DETERMINANT   ******************************
 * QUANTUM FORCE F = 2/(psi)*del(psi)
 * fill the qm_force matrix wth the derivative of all
 * the particles, wrt dimension
 * 
 * del(psi)/psi = del(SD)/SD + del(J)/J 
 *
 *
 */
void atom::quantum_force(double**r,double**qm_force, double **D_up, double **D_down, double alpha, double beta){
  //zero
  for(int p=0; p<no_of_particles; p++){
    for(int q=0; q<dimension; q++){
      qm_force[p][q]=0;
    }
  }   
  //SLATER
  for(int p=0; p<no_of_particles; p++){
    if(p<no_of_particles/2){
      for(int q=0; q<dimension; q++){
	for(int l=0; l<no_of_particles/2; l++){
	  qm_force[p][q] += 2*gradient_phi(r,alpha,q,l,p)*D_up[l][p];
	}
      }
    }
    else{
      for(int q=0; q<dimension; q++){
	for(int l=0; l<no_of_particles/2; l++){
	  qm_force[p][q] += 2*gradient_phi(r,alpha,q,l,p)*D_down[l][p-no_of_particles/2];
	}
      }
    }
  }

  //JASTROW
 for(int p=0; p<no_of_particles; p++){
    for(int q=0; q<dimension; q++){
      for(int l=p+1; l<no_of_particles; l++){
	qm_force[p][q] += 2*gradient_jastrow(r,p,l,q,beta);
      }
      for(int l=0; l<p; l++){
	qm_force[p][q] += 2*gradient_jastrow(r,p,l,q,beta);
      }
    }
  }  //*/   

}//end quantum_force()

/*************************    GRADIENT SLATER DETERMINANT WAVEFUNCTION    ********************
 *
 *  This is a vector.
 *  Input: Position-matrix r
 *        coordinate x_i
 *        State n
 *        Particle m 
 *        +variational parameters
 *        
 *        gradient_phi(radius,alpha,coordinate l,state j,particle i)
 *
 */
double atom::gradient_phi(double**r,double alpha, int x_i, int n, int m){
double r_single_particle;
  double *argument1,*argument2;
  double temp=0;

  argument1 = new double[no_of_particles];
  argument2 = new double[no_of_particles];
  //zero arguments
  for (int i=0; i<no_of_particles; i++){ 
    argument1[i] = 0.0;
    argument2[i] = 0.0;
  }   
  
  for(int i=0; i<no_of_particles; i++){
    r_single_particle=0;
    for(int j=0; j<dimension; j++){
      r_single_particle += r[i][j]*r[i][j];//r^2=x^2+y^2+z^2
    }
    argument1[i] += sqrt(r_single_particle);
  }
 		
  if(n==0){
	if(charge==2){
		
		temp = root_helium[0][1]*rooth_factors[0]*gradient_psi1s(argument1[m], r[m][x_i], root_helium[0][0]) + root_helium[1][1]*rooth_factors[1]*gradient_psi1s(argument1[m], r[m][x_i], root_helium[1][0]) + root_helium[2][1]*rooth_factors[2]*gradient_psi1s(argument1[m], r[m][x_i], root_helium[2][0]) + root_helium[3][1]*rooth_factors[3]*gradient_psi1s(argument1[m], r[m][x_i], root_helium[3][0]) + root_helium[4][1]*rooth_factors[4]*gradient_psi1s(argument1[m], r[m][x_i], root_helium[4][0]);
		}
	if(charge==4){
		temp = root_beryllium[0][1]*rooth_factors[0]*gradient_psi1s(argument1[m], r[m][x_i], root_beryllium[0][0]) + root_beryllium[1][1]*rooth_factors[1]*gradient_psi1s(argument1[m], r[m][x_i], root_beryllium[1][0]) + root_beryllium[2][1]*rooth_factors[2]*gradient_psi2s(argument1[m], r[m][x_i], root_beryllium[2][0]) + root_beryllium[3][1]*rooth_factors[3]*gradient_psi2s(argument1[m], r[m][x_i], root_beryllium[3][0]) + root_beryllium[4][1]*rooth_factors[4]*gradient_psi2s(argument1[m], r[m][x_i], root_beryllium[4][0])+ root_beryllium[5][1]*rooth_factors[5]*gradient_psi2s(argument1[m], r[m][x_i], root_beryllium[5][0]);
		}
	if(charge==10){
		temp = root_neon[0][1]*rooth_factors[0]*gradient_psi1s(argument1[m], r[m][x_i], root_neon[0][0]) + root_neon[1][1]*rooth_factors[1]*gradient_psi1s(argument1[m], r[m][x_i], root_neon[1][0]) + root_neon[2][1]*rooth_factors[2]*gradient_psi2s(argument1[m], r[m][x_i], root_neon[2][0]) + root_neon[3][1]*rooth_factors[3]*gradient_psi2s(argument1[m], r[m][x_i], root_neon[3][0]) + root_neon[4][1]*rooth_factors[4]*gradient_psi2s(argument1[m], r[m][x_i], root_neon[4][0])+ root_neon[5][1]*rooth_factors[5]*gradient_psi2s(argument1[m], r[m][x_i], root_neon[5][0]); 
		
		}

	if(charge==12){
		root_magnesium[0][1]*rooth_factors[0]*gradient_psi1s(argument1[m], r[m][x_i], root_magnesium[0][0]) + root_magnesium[1][1]*rooth_factors[1]*gradient_psi3s(argument1[m], r[m][x_i], root_magnesium[1][0]) + root_magnesium[2][1]*rooth_factors[2]*gradient_psi3s(argument1[m], r[m][x_i], root_magnesium[2][0]) + root_magnesium[3][1]*rooth_factors[3]*gradient_psi3s(argument1[m], r[m][x_i], root_magnesium[3][0]) + root_magnesium[4][1]*rooth_factors[4]*gradient_psi3s(argument1[m], r[m][x_i], root_magnesium[4][0]) + root_magnesium[5][1]*rooth_factors[5]*gradient_psi3s(argument1[m], r[m][x_i], root_magnesium[5][0]) + root_magnesium[6][1]*rooth_factors[6]*gradient_psi3s(argument1[m], r[m][x_i], root_magnesium[6][0]) + root_magnesium[7][1]*rooth_factors[7]*gradient_psi3s(argument1[m], r[m][x_i], root_magnesium[7][0]);
		} 
	if(charge==14){
		root_silisium[0][1]*rooth_factors[0]*gradient_psi1s(argument1[m], r[m][x_i], root_magnesium[0][0]) + root_silisium[1][1]*rooth_factors[1]*gradient_psi3s(argument1[m], r[m][x_i], root_silisium[1][0]) + root_silisium[2][1]*rooth_factors[2]*gradient_psi3s(argument1[m], r[m][x_i], root_silisium[2][0]) + root_silisium[3][1]*rooth_factors[3]*gradient_psi3s(argument1[m], r[m][x_i], root_silisium[3][0]) + root_silisium[4][1]*rooth_factors[4]*gradient_psi3s(argument1[m], r[m][x_i], root_silisium[4][0]) + root_silisium[5][1]*rooth_factors[5]*gradient_psi3s(argument1[m], r[m][x_i], root_silisium[5][0]) + root_silisium[6][1]*rooth_factors[6]*gradient_psi3s(argument1[m], r[m][x_i], root_silisium[6][0]) + root_silisium[7][1]*rooth_factors[7]*gradient_psi3s(argument1[m], r[m][x_i], root_silisium[7][0]);
		} 

    
  } 
  if(n==1){
	if(charge==4){
		
		temp = root_beryllium[0][2]*rooth_factors[0]*gradient_psi1s(argument1[m], r[m][x_i], root_beryllium[0][0]) + root_beryllium[1][2]*rooth_factors[1]*gradient_psi1s(argument1[m], r[m][x_i], root_beryllium[1][0]) + root_beryllium[2][2]*rooth_factors[2]*gradient_psi2s(argument1[m], r[m][x_i], root_beryllium[2][0]) + root_beryllium[3][2]*rooth_factors[3]*gradient_psi2s(argument1[m], r[m][x_i], root_beryllium[3][0]) + root_beryllium[4][2]*rooth_factors[4]*gradient_psi2s(argument1[m], r[m][x_i], root_beryllium[4][0]) + root_beryllium[5][2]*rooth_factors[5]*gradient_psi2s(argument1[m], r[m][x_i], root_beryllium[5][0]);
		}
	
	if(charge==10){
		
		temp = root_neon[0][2]*rooth_factors[0]*gradient_psi1s(argument1[m], r[m][x_i], root_neon[0][0]) + root_neon[1][2]*rooth_factors[1]*gradient_psi1s(argument1[m], r[m][x_i], root_neon[1][0]) + root_neon[2][2]*rooth_factors[2]*gradient_psi2s(argument1[m], r[m][x_i], root_neon[2][0]) + root_neon[3][2]*rooth_factors[3]*gradient_psi2s(argument1[m], r[m][x_i], root_neon[3][0]) + root_neon[4][2]*rooth_factors[4]*gradient_psi2s(argument1[m], r[m][x_i], root_neon[4][0]) + root_neon[5][2]*rooth_factors[5]*gradient_psi2s(argument1[m], r[m][x_i], root_neon[5][0]);
		}
	
	if(charge==12){

		temp = root_magnesium[0][2]*rooth_factors[0]*gradient_psi1s(argument1[m], r[m][x_i], root_magnesium[0][0]) + root_magnesium[1][2]*rooth_factors[1]*gradient_psi3s(argument1[m], r[m][x_i], root_magnesium[1][0]) + root_magnesium[2][2]*rooth_factors[2]*gradient_psi3s(argument1[m], r[m][x_i], root_magnesium[2][0]) + root_magnesium[3][2]*rooth_factors[3]*gradient_psi3s(argument1[m], r[m][x_i], root_magnesium[3][0]) + root_magnesium[4][2]*rooth_factors[4]*gradient_psi3s(argument1[m], r[m][x_i], root_magnesium[4][0]) + root_magnesium[5][2]*rooth_factors[5]*gradient_psi3s(argument1[m], r[m][x_i], root_magnesium[5][0]) + root_magnesium[6][2]*rooth_factors[6]*gradient_psi3s(argument1[m], r[m][x_i], root_magnesium[6][0]) + root_magnesium[7][2]*rooth_factors[7]*gradient_psi3s(argument1[m], r[m][x_i], root_magnesium[7][0]); 


		}
	if(charge==14){

		temp = root_silisium[0][2]*rooth_factors[0]*gradient_psi1s(argument1[m], r[m][x_i], root_silisium[0][0]) + root_silisium[1][2]*rooth_factors[1]*gradient_psi3s(argument1[m], r[m][x_i], root_silisium[1][0]) + root_silisium[2][2]*rooth_factors[2]*gradient_psi3s(argument1[m], r[m][x_i], root_silisium[2][0]) + root_silisium[3][2]*rooth_factors[3]*gradient_psi3s(argument1[m], r[m][x_i], root_silisium[3][0]) + root_silisium[4][2]*rooth_factors[4]*gradient_psi3s(argument1[m], r[m][x_i], root_silisium[4][0]) + root_silisium[5][2]*rooth_factors[5]*gradient_psi3s(argument1[m], r[m][x_i], root_silisium[5][0]) + root_silisium[6][2]*rooth_factors[6]*gradient_psi3s(argument1[m], r[m][x_i], root_silisium[6][0]) + root_silisium[7][2]*rooth_factors[7]*gradient_psi3s(argument1[m], r[m][x_i], root_silisium[7][0]); 


		}
  } 
  if(n==2){
    if(x_i==0){
	if(charge==10){
		temp = root_neon[0][4]*rooth_factors[6]*gradient_psi2p_x(argument1[m], root_neon[0][3], r[m][x_i]) + root_neon[1][4]*rooth_factors[7]*gradient_psi2p_x(argument1[m], root_neon[1][3], r[m][x_i]) + root_neon[2][4]*rooth_factors[8]*gradient_psi2p_x(argument1[m], root_neon[2][3], r[m][x_i]) + root_neon[3][4]*rooth_factors[9]*gradient_psi2p_x(argument1[m], root_neon[3][3], r[m][x_i]); 

		}
	if(charge == 12){
		temp = 	root_magnesium[0][5]*rooth_factors[8]*gradient_psi2p_x(argument1[m], root_magnesium[0][4], r[m][x_i]) + root_magnesium[1][5]*rooth_factors[9]*gradient_psi4p_x(argument1[m], root_magnesium[1][4], r[m][x_i]) + root_magnesium[2][5]*rooth_factors[10]*gradient_psi4p_x(argument1[m], root_magnesium[2][4], r[m][x_i]) + root_magnesium[3][5]*rooth_factors[11]*gradient_psi4p_x(argument1[m], root_magnesium[3][4], r[m][x_i]) + root_magnesium[4][5]*rooth_factors[12]*gradient_psi4p_x(argument1[m], root_magnesium[4][4], r[m][x_i]);	
		}
	if(charge == 14){
		temp = root_silisium[0][5]*rooth_factors[8]*gradient_psi2p_x(argument1[m], root_silisium[0][4], r[m][x_i]) + root_silisium[1][5]*rooth_factors[9]*gradient_psi4p_x(argument1[m], root_silisium[1][4], r[m][x_i]) + root_silisium[2][5]*rooth_factors[10]*gradient_psi4p_x(argument1[m], root_silisium[2][4], r[m][x_i]) + root_silisium[3][5]*rooth_factors[11]*gradient_psi4p_x(argument1[m], root_silisium[3][4], r[m][x_i]) + root_silisium[4][5]*rooth_factors[12]*gradient_psi4p_x(argument1[m], root_silisium[4][4], r[m][x_i]) + root_silisium[5][5]*rooth_factors[13]*gradient_psi4p_x(argument1[m], root_silisium[5][4], r[m][x_i]) + root_silisium[6][5]*rooth_factors[14]*gradient_psi4p_x(argument1[m], root_silisium[6][4], r[m][x_i]) + root_silisium[7][5]*rooth_factors[15]*gradient_psi4p_x(argument1[m], root_silisium[7][4], r[m][x_i]);

		}
     
    }
    else{
	if(charge==10){
		temp = root_neon[0][4]*rooth_factors[6]*gradient_psi2p_y(argument1[m], root_neon[0][3],r[m][0],  r[m][x_i]) + root_neon[1][4]*rooth_factors[7]*gradient_psi2p_y(argument1[m], root_neon[1][3], r[m][0], r[m][x_i]) + root_neon[2][4]*rooth_factors[8]*gradient_psi2p_y(argument1[m], root_neon[2][3], r[m][0], r[m][x_i]) + root_neon[3][4]*rooth_factors[9]*gradient_psi2p_y(argument1[m], root_neon[3][3], r[m][0], r[m][x_i]); 

		}
	if(charge == 12){
		temp = 	root_magnesium[0][5]*rooth_factors[8]*gradient_psi2p_y(argument1[m], root_magnesium[0][4], r[m][0], r[m][x_i]) + root_magnesium[1][5]*rooth_factors[9]*gradient_psi4p_y(argument1[m], root_magnesium[1][4], r[m][0], r[m][x_i]) + root_magnesium[2][5]*rooth_factors[10]*gradient_psi4p_y(argument1[m], root_magnesium[2][4], r[m][0], r[m][x_i]) + root_magnesium[3][5]*rooth_factors[11]*gradient_psi4p_y(argument1[m], root_magnesium[3][4], r[m][0], r[m][x_i]) + root_magnesium[4][5]*rooth_factors[12]*gradient_psi4p_y(argument1[m], root_magnesium[4][4], r[m][0], r[m][x_i]);	
		}
	if(charge == 14){
		temp = root_silisium[0][5]*rooth_factors[8]*gradient_psi2p_y(argument1[m], root_silisium[0][4], r[m][0],r[m][x_i]) + root_silisium[1][5]*rooth_factors[9]*gradient_psi4p_y(argument1[m], root_silisium[1][4], r[m][0],r[m][x_i]) + root_silisium[2][5]*rooth_factors[10]*gradient_psi4p_y(argument1[m], root_silisium[2][4], r[m][0],r[m][x_i]) + root_silisium[3][5]*rooth_factors[11]*gradient_psi4p_y(argument1[m], root_silisium[3][4], r[m][0],r[m][x_i]) + root_silisium[4][5]*rooth_factors[12]*gradient_psi4p_y(argument1[m], root_silisium[4][4], r[m][0],r[m][x_i]) + root_silisium[5][5]*rooth_factors[13]*gradient_psi4p_y(argument1[m], root_silisium[5][4], r[m][0],r[m][x_i]) + root_silisium[6][5]*rooth_factors[14]*gradient_psi4p_y(argument1[m], root_silisium[6][4], r[m][0],r[m][x_i]) + root_silisium[7][5]*rooth_factors[15]*gradient_psi4p_y(argument1[m], root_silisium[7][4], r[m][0],r[m][x_i]);

		}
    }
  }
  if(n==3){
    if(x_i==1){
      if(charge==10){
		temp = root_neon[0][4]*rooth_factors[6]*gradient_psi2p_x(argument1[m], root_neon[0][3], r[m][x_i]) + root_neon[1][4]*rooth_factors[7]*gradient_psi2p_x(argument1[m], root_neon[1][3], r[m][x_i]) + root_neon[2][4]*rooth_factors[8]*gradient_psi2p_x(argument1[m], root_neon[2][3], r[m][x_i]) + root_neon[3][4]*rooth_factors[9]*gradient_psi2p_x(argument1[m], root_neon[3][3], r[m][x_i]); 

		}
	if(charge == 12){
		temp = 	root_magnesium[0][5]*rooth_factors[8]*gradient_psi2p_x(argument1[m], root_magnesium[0][4], r[m][x_i]) + root_magnesium[1][5]*rooth_factors[9]*gradient_psi4p_x(argument1[m], root_magnesium[1][4], r[m][x_i]) + root_magnesium[2][5]*rooth_factors[10]*gradient_psi4p_x(argument1[m], root_magnesium[2][4], r[m][x_i]) + root_magnesium[3][5]*rooth_factors[11]*gradient_psi4p_x(argument1[m], root_magnesium[3][4], r[m][x_i]) + root_magnesium[4][5]*rooth_factors[12]*gradient_psi4p_x(argument1[m], root_magnesium[4][4], r[m][x_i]);	
		}
	if(charge == 14){
		temp = root_silisium[0][5]*rooth_factors[8]*gradient_psi2p_x(argument1[m], root_silisium[0][4], r[m][x_i]) + root_silisium[1][5]*rooth_factors[9]*gradient_psi4p_x(argument1[m], root_silisium[1][4], r[m][x_i]) + root_silisium[2][5]*rooth_factors[10]*gradient_psi4p_x(argument1[m], root_silisium[2][4], r[m][x_i]) + root_silisium[3][5]*rooth_factors[11]*gradient_psi4p_x(argument1[m], root_silisium[3][4], r[m][x_i]) + root_silisium[4][5]*rooth_factors[12]*gradient_psi4p_x(argument1[m], root_silisium[4][4], r[m][x_i]) + root_silisium[5][5]*rooth_factors[13]*gradient_psi4p_x(argument1[m], root_silisium[5][4], r[m][x_i]) + root_silisium[6][5]*rooth_factors[14]*gradient_psi4p_x(argument1[m], root_silisium[6][4], r[m][x_i]) + root_silisium[7][5]*rooth_factors[15]*gradient_psi4p_x(argument1[m], root_silisium[7][4], r[m][x_i]);

		}
     
    }
    else{
      if(charge==10){
		temp = root_neon[0][4]*rooth_factors[6]*gradient_psi2p_y(argument1[m], root_neon[0][3],r[m][1],  r[m][x_i]) + root_neon[1][4]*rooth_factors[7]*gradient_psi2p_y(argument1[m], root_neon[1][3], r[m][1], r[m][x_i]) + root_neon[2][4]*rooth_factors[8]*gradient_psi2p_y(argument1[m], root_neon[2][3], r[m][1], r[m][x_i]) + root_neon[3][4]*rooth_factors[9]*gradient_psi2p_y(argument1[m], root_neon[3][3], r[m][1], r[m][x_i]); 

		}
	if(charge == 12){
		temp = 	root_magnesium[0][5]*rooth_factors[8]*gradient_psi2p_y(argument1[m], root_magnesium[0][4], r[m][1], r[m][x_i]) + root_magnesium[1][5]*rooth_factors[9]*gradient_psi4p_y(argument1[m], root_magnesium[1][4], r[m][1], r[m][x_i]) + root_magnesium[2][5]*rooth_factors[10]*gradient_psi4p_y(argument1[m], root_magnesium[2][4], r[m][1], r[m][x_i]) + root_magnesium[3][5]*rooth_factors[11]*gradient_psi4p_y(argument1[m], root_magnesium[3][4], r[m][1], r[m][x_i]) + root_magnesium[4][5]*rooth_factors[12]*gradient_psi4p_y(argument1[m], root_magnesium[4][4], r[m][1], r[m][x_i]);	
		}
	if(charge == 14){
		temp = root_silisium[0][5]*rooth_factors[8]*gradient_psi2p_y(argument1[m], root_silisium[0][4], r[m][1],r[m][x_i]) + root_silisium[1][5]*rooth_factors[9]*gradient_psi4p_y(argument1[m], root_silisium[1][4], r[m][1],r[m][x_i]) + root_silisium[2][5]*rooth_factors[10]*gradient_psi4p_y(argument1[m], root_silisium[2][4], r[m][1],r[m][x_i]) + root_silisium[3][5]*rooth_factors[11]*gradient_psi4p_y(argument1[m], root_silisium[3][4], r[m][1],r[m][x_i]) + root_silisium[4][5]*rooth_factors[12]*gradient_psi4p_y(argument1[m], root_silisium[4][4], r[m][1],r[m][x_i]) + root_silisium[5][5]*rooth_factors[13]*gradient_psi4p_y(argument1[m], root_silisium[5][4], r[m][1],r[m][x_i]) + root_silisium[6][5]*rooth_factors[14]*gradient_psi4p_y(argument1[m], root_silisium[6][4], r[m][1],r[m][x_i]) + root_silisium[7][5]*rooth_factors[15]*gradient_psi4p_y(argument1[m], root_silisium[7][4], r[m][1],r[m][x_i]);

		}
    }
  }
  if(n==4){
   if(x_i==2){
      if(charge==10){
		temp = root_neon[0][4]*rooth_factors[6]*gradient_psi2p_x(argument1[m], root_neon[0][3], r[m][x_i]) + root_neon[1][4]*rooth_factors[7]*gradient_psi2p_x(argument1[m], root_neon[1][3], r[m][x_i]) + root_neon[2][4]*rooth_factors[8]*gradient_psi2p_x(argument1[m], root_neon[2][3], r[m][x_i]) + root_neon[3][4]*rooth_factors[9]*gradient_psi2p_x(argument1[m], root_neon[3][3], r[m][x_i]); 

		}
	if(charge == 12){
		temp = 	root_magnesium[0][5]*rooth_factors[8]*gradient_psi2p_x(argument1[m], root_magnesium[0][4], r[m][x_i]) + root_magnesium[1][5]*rooth_factors[9]*gradient_psi4p_x(argument1[m], root_magnesium[1][4], r[m][x_i]) + root_magnesium[2][5]*rooth_factors[10]*gradient_psi4p_x(argument1[m], root_magnesium[2][4], r[m][x_i]) + root_magnesium[3][5]*rooth_factors[11]*gradient_psi4p_x(argument1[m], root_magnesium[3][4], r[m][x_i]) + root_magnesium[4][5]*rooth_factors[12]*gradient_psi4p_x(argument1[m], root_magnesium[4][4], r[m][x_i]);	
		}
	if(charge == 14){
		temp = root_silisium[0][5]*rooth_factors[8]*gradient_psi2p_x(argument1[m], root_silisium[0][4], r[m][x_i]) + root_silisium[1][5]*rooth_factors[9]*gradient_psi4p_x(argument1[m], root_silisium[1][4], r[m][x_i]) + root_silisium[2][5]*rooth_factors[10]*gradient_psi4p_x(argument1[m], root_silisium[2][4], r[m][x_i]) + root_silisium[3][5]*rooth_factors[11]*gradient_psi4p_x(argument1[m], root_silisium[3][4], r[m][x_i]) + root_silisium[4][5]*rooth_factors[12]*gradient_psi4p_x(argument1[m], root_silisium[4][4], r[m][x_i]) + root_silisium[5][5]*rooth_factors[13]*gradient_psi4p_x(argument1[m], root_silisium[5][4], r[m][x_i]) + root_silisium[6][5]*rooth_factors[14]*gradient_psi4p_x(argument1[m], root_silisium[6][4], r[m][x_i]) + root_silisium[7][5]*rooth_factors[15]*gradient_psi4p_x(argument1[m], root_silisium[7][4], r[m][x_i]);

		}
    }
    else{
      if(charge==10){
		temp = root_neon[0][4]*rooth_factors[6]*gradient_psi2p_y(argument1[m], root_neon[0][3],r[m][2],  r[m][x_i]) + root_neon[1][4]*rooth_factors[7]*gradient_psi2p_y(argument1[m], root_neon[1][3], r[m][2], r[m][x_i]) + root_neon[2][4]*rooth_factors[8]*gradient_psi2p_y(argument1[m], root_neon[2][3], r[m][2], r[m][x_i]) + root_neon[3][4]*rooth_factors[9]*gradient_psi2p_y(argument1[m], root_neon[3][3], r[m][2], r[m][x_i]); 

		}
	if(charge == 12){
		temp = 	root_magnesium[0][5]*rooth_factors[8]*gradient_psi2p_y(argument1[m], root_magnesium[0][4], r[m][2], r[m][x_i]) + root_magnesium[1][5]*rooth_factors[9]*gradient_psi4p_y(argument1[m], root_magnesium[1][4], r[m][2], r[m][x_i]) + root_magnesium[2][5]*rooth_factors[10]*gradient_psi4p_y(argument1[m], root_magnesium[2][4], r[m][2], r[m][x_i]) + root_magnesium[3][5]*rooth_factors[11]*gradient_psi4p_y(argument1[m], root_magnesium[3][4], r[m][2], r[m][x_i]) + root_magnesium[4][5]*rooth_factors[12]*gradient_psi4p_y(argument1[m], root_magnesium[4][4], r[m][2], r[m][x_i]);	
		}
	if(charge == 14){
		temp = root_silisium[0][5]*rooth_factors[8]*gradient_psi2p_y(argument1[m], root_silisium[0][4], r[m][2],r[m][x_i]) + root_silisium[1][5]*rooth_factors[9]*gradient_psi4p_y(argument1[m], root_silisium[1][4], r[m][2],r[m][x_i]) + root_silisium[2][5]*rooth_factors[10]*gradient_psi4p_y(argument1[m], root_silisium[2][4], r[m][2],r[m][x_i]) + root_silisium[3][5]*rooth_factors[11]*gradient_psi4p_y(argument1[m], root_silisium[3][4], r[m][2],r[m][x_i]) + root_silisium[4][5]*rooth_factors[12]*gradient_psi4p_y(argument1[m], root_silisium[4][4], r[m][2],r[m][x_i]) + root_silisium[5][5]*rooth_factors[13]*gradient_psi4p_y(argument1[m], root_silisium[5][4], r[m][2],r[m][x_i]) + root_silisium[6][5]*rooth_factors[14]*gradient_psi4p_y(argument1[m], root_silisium[6][4], r[m][2],r[m][x_i]) + root_silisium[7][5]*rooth_factors[15]*gradient_psi4p_y(argument1[m], root_silisium[7][4], r[m][2],r[m][x_i]);

		}
    }
  }
  if(n==5){
	if(charge==12){
		temp= root_magnesium[0][3]*rooth_factors[0]*gradient_psi1s(argument1[m], r[m][x_i], root_magnesium[0][0]) + root_magnesium[1][3]*rooth_factors[1]*gradient_psi3s(argument1[m], r[m][x_i], root_magnesium[1][0]) + root_magnesium[2][3]*rooth_factors[2]*gradient_psi3s(argument1[m], r[m][x_i], root_magnesium[2][0]) + root_magnesium[3][3]*rooth_factors[3]*gradient_psi3s(argument1[m], r[m][x_i], root_magnesium[3][0]) + root_magnesium[4][3]*rooth_factors[4]*gradient_psi3s(argument1[m], r[m][x_i], root_magnesium[4][0]) + root_magnesium[5][3]*rooth_factors[5]*gradient_psi3s(argument1[m], r[m][x_i], root_magnesium[5][0]) + root_magnesium[6][3]*rooth_factors[6]*gradient_psi3s(argument1[m], r[m][x_i], root_magnesium[6][0]) + root_magnesium[7][3]*rooth_factors[7]*gradient_psi3s(argument1[m], r[m][x_i], root_magnesium[7][0]); 

		}
	if(charge==14){
		temp = root_silisium[0][3]*rooth_factors[0]*gradient_psi1s(argument1[m], r[m][x_i], root_silisium[0][0]) +  root_silisium[1][3]*rooth_factors[1]*gradient_psi3s(argument1[m], r[m][x_i], root_silisium[1][0]) +  root_silisium[2][3]*rooth_factors[2]*gradient_psi3s(argument1[m], r[m][x_i], root_silisium[2][0]) +  root_silisium[3][3]*rooth_factors[3]*gradient_psi3s(argument1[m], r[m][x_i], root_silisium[3][0]) +  root_silisium[4][3]*rooth_factors[4]*gradient_psi3s(argument1[m], r[m][x_i], root_silisium[4][0]) +  root_silisium[5][3]*rooth_factors[5]*gradient_psi3s(argument1[m], r[m][x_i], root_silisium[5][0]) +  root_silisium[6][3]*rooth_factors[6]*gradient_psi3s(argument1[m], r[m][x_i], root_silisium[6][0]) +  root_silisium[7][3]*rooth_factors[7]*gradient_psi3s(argument1[m], r[m][x_i], root_silisium[7][0]);


		}
 } 

 //*/
 //VELGER GRUNNTILSTANDEN TIL Å VÆRE m_l=0, dvs solid harmonic= y.
  if(n==6){
      if(x_i==1){
	 temp = root_silisium[0][6]*rooth_factors[8]*gradient_psi2p_x(argument1[m], root_silisium[0][4], r[m][x_i]) + root_silisium[1][6]*rooth_factors[9]*gradient_psi4p_x(argument1[m], root_silisium[1][4], r[m][x_i]) + root_silisium[2][6]*rooth_factors[10]*gradient_psi4p_x(argument1[m], root_silisium[2][4], r[m][x_i]) + root_silisium[3][6]*rooth_factors[11]*gradient_psi4p_x(argument1[m], root_silisium[3][4], r[m][x_i]) + root_silisium[4][6]*rooth_factors[12]*gradient_psi4p_x(argument1[m], root_silisium[4][4], r[m][x_i]) + root_silisium[5][6]*rooth_factors[13]*gradient_psi4p_x(argument1[m], root_silisium[5][4], r[m][x_i]) + root_silisium[6][6]*rooth_factors[14]*gradient_psi4p_x(argument1[m], root_silisium[6][4], r[m][x_i]) + root_silisium[7][6]*rooth_factors[15]*gradient_psi4p_x(argument1[m], root_silisium[7][4], r[m][x_i]); 


         
         }
      else{
	  temp = root_silisium[0][6]*rooth_factors[8]*gradient_psi2p_y(argument1[m], root_silisium[0][4], r[m][1], r[m][x_i]) + root_silisium[1][6]*rooth_factors[9]*gradient_psi4p_y(argument1[m], root_silisium[1][4], r[m][1], r[m][x_i]) + root_silisium[2][6]*rooth_factors[10]*gradient_psi4p_y(argument1[m], root_silisium[2][4], r[m][1], r[m][x_i]) + root_silisium[3][6]*rooth_factors[11]*gradient_psi4p_y(argument1[m], root_silisium[3][4], r[m][1], r[m][x_i]) + root_silisium[4][6]*rooth_factors[12]*gradient_psi4p_y(argument1[m], root_silisium[4][4], r[m][1], r[m][x_i]) + root_silisium[5][6]*rooth_factors[13]*gradient_psi4p_y(argument1[m], root_silisium[5][4], r[m][1], r[m][x_i]) + root_silisium[6][6]*rooth_factors[14]*gradient_psi4p_y(argument1[m], root_silisium[6][4], r[m][1], r[m][x_i]) + root_silisium[7][6]*rooth_factors[15]*gradient_psi4p_y(argument1[m], root_silisium[7][4], r[m][1], r[m][x_i]); 

        
	 }
  }

  delete [] argument1;
  delete [] argument2;
  return temp;

}//end gradient_phi()

/*******************    LOCAL ENERGY SLATER DETERMINANT     ******************************
 *
 *
 */
double atom::local_energy(double **r,double **D_up, double **D_down, double alpha, double beta){
  //local variables
  double E_local, E_kinetic, E_potential, r_single_particle, r_12;
  //reset values
  E_local    =0;
  E_kinetic  =0;
  E_potential=0;

  //  KINETIC ENERGY FOR SLATER-DETERMINANT
  //  E_kinetic(local) = Laplace(det)/det.
  //  Sum over all particles and add the kinetic energy for each

  //LAPLACE SLATER-DETERMINANT ANALYTIC
  for(int i=0; i<no_of_particles; i++){
    if(i<no_of_particles/2){
      for(int j=0; j<no_of_particles/2; j++){
	E_kinetic -= laplace_phi(r,alpha,beta,j,i)*D_up[j][i];
      }
    }
    else{
      for(int j=0; j<no_of_particles/2; j++){
	E_kinetic -= laplace_phi(r,alpha,beta,j,i)*D_down[j][i-no_of_particles/2];
      }
    }
  }
  
 //CROSS TERM GRADIENT JASTROW AND SLATER ANALYTIC
  double tempE = 0;
  double **temp1 = create_matrix(no_of_particles, dimension);  
  double **temp2 = create_matrix(no_of_particles, dimension);  
  //zero
  for(int p=0; p<no_of_particles; p++){
    for(int q=0; q<dimension; q++){
      temp1[p][q]=0;
      temp2[p][q]=0;
    }
  }   
  //slater
  for(int p=0; p<no_of_particles; p++){
    if(p<no_of_particles/2){
      for(int q=0; q<dimension; q++){
	for(int l=0; l<no_of_particles/2; l++){
	  temp1[p][q] += gradient_phi(r,alpha,q,l,p)*D_up[l][p];
	}
      }
    }
    else{
      for(int q=0; q<dimension; q++){
	for(int l=0; l<no_of_particles/2; l++){
	  temp1[p][q] += gradient_phi(r,alpha,q,l,p)*D_down[l][p-no_of_particles/2];
	}
      }
    }
  }
  //jastrow
  for(int p=0; p<no_of_particles; p++){
    for(int q=0; q<dimension; q++){
      for(int l=p+1; l<no_of_particles; l++){
	temp2[p][q] += gradient_jastrow(r,p,l,q,beta);
      }
      for(int l=0; l<p; l++){
	temp2[p][q] += gradient_jastrow(r,p,l,q,beta);
      }
    }
  }   
  for(int p=0; p<no_of_particles; p++){
    for(int q=0; q<dimension; q++){
      tempE += temp1[p][q]*temp2[p][q];
    }
  }   
  //LAPLACE JASTROW ANALYTIC
  //(del J/J)^2
  double tempE2=0;
  for(int p=0; p<no_of_particles; p++){
    for(int q=0; q<dimension; q++){
      tempE2 += temp2[p][q]*temp2[p][q];
    }
  }   
  double tempE3=0;
  double a=0;
  //det andre leddet
  for(int p=0; p<no_of_particles; p++){
    for (int k = 0; k < no_of_particles; k++) {
      if ( k != p) {
	if(((p < no_of_particles/2) && (k <no_of_particles/2)) || ((p>=no_of_particles/2 && k>=no_of_particles/2))){
	   a=0.25;
	   tempE3 += 2*a/diffR(r,p,k)/pow((1+beta*diffR(r,p,k)),3); 
	 }
	 else{
	   a=0.5;
	   tempE3 += 2*a/diffR(r,p,k)/pow((1+beta*diffR(r,p,k)),3); 
	 }
      } 
    }
  }
  
  delete_matrix(temp1,no_of_particles);
  delete_matrix(temp2,no_of_particles);
 
  E_kinetic -= 2*tempE;
  E_kinetic -= tempE2;
  E_kinetic -= tempE3;   //*/

  E_kinetic *= 0.5;

  //POTENTIAL ENERGY
  //Coulomb potential
  for(int i=0; i<no_of_particles; i++){
    r_single_particle=0;
    for(int j=0; j<dimension; j++){
      r_single_particle += r[i][j]*r[i][j];//r^2=x^2+y^2+z^2
    }
    E_potential -=charge/sqrt(r_single_particle); 
  }
  //ee-interaction
  for(int i=0; i<no_of_particles-1; i++){
    for(int j=i+1; j<no_of_particles; j++){
      r_12=0;
      for(int k=0; k<dimension; k++){
	r_12 +=(r[i][k]-r[j][k])*(r[i][k]-r[j][k]);
      }
      E_potential +=1/sqrt(r_12);
    }
  }//*/
  E_local = E_potential + E_kinetic;
  
  return E_local;
}//end local_energy()

/*************************    LAPLACE SLATER DETERMINANT WAVEFUNCTION    ********************
 *
 * This is a scalar.
 * Input: Position-matrix r
 *        State n
 *        Particle m 
 *        +variational parameters
 *        
 *        laplace_phi(r,alpha,beta,j,i)
 *
 */
double atom::laplace_phi(double**r,double alpha,double beta, int n, int m){
  double r_single_particle;
  double *argument1,*argument2;
  double temp=0;

  argument1 = new double[no_of_particles];
  argument2 = new double[no_of_particles];
  //zero arguments
  for (int i=0; i<no_of_particles; i++){ 
    argument1[i] = 0.0;
    argument2[i] = 0.0;
  }   
  
  for(int i=0; i<no_of_particles; i++){
    r_single_particle=0;
    for(int j=0; j<dimension; j++){
      r_single_particle += r[i][j]*r[i][j];//r^2=x^2+y^2+z^2
    }
    argument1[i] += sqrt(r_single_particle);
  }
 			
  if(n==0){
    	if(charge==2){
		temp = root_helium[0][1]*rooth_factors[0]*laplace_psi1s(argument1[m], root_helium[0][0]) + root_helium[1][1]*rooth_factors[1]*laplace_psi1s(argument1[m], root_helium[1][0])  + root_helium[2][1]*rooth_factors[2]*laplace_psi1s(argument1[m], root_helium[2][0]) + root_helium[3][1]*rooth_factors[3]*laplace_psi1s(argument1[m], root_helium[3][0]) + root_helium[4][1]*rooth_factors[4]*laplace_psi1s(argument1[m], root_helium[4][0]); 
		}
	if(charge==4){
		temp = root_beryllium[0][1]*rooth_factors[0]*laplace_psi1s(argument1[m], root_beryllium[0][0]) + root_beryllium[1][1]*rooth_factors[1]*laplace_psi1s(argument1[m], root_beryllium[1][0]) + root_beryllium[2][1]*rooth_factors[2]*laplace_psi2s(argument1[m], root_beryllium[2][0]) + root_beryllium[3][1]*rooth_factors[3]*laplace_psi2s(argument1[m], root_beryllium[3][0]) + root_beryllium[4][1]*rooth_factors[4]*laplace_psi2s(argument1[m], root_beryllium[4][0]) + root_beryllium[5][1]*rooth_factors[5]*laplace_psi2s(argument1[m], root_beryllium[5][0]);

		}	
	if(charge==10){
		temp = root_neon[0][1]*rooth_factors[0]*laplace_psi1s(argument1[m], root_neon[0][0]) + root_neon[1][1]*rooth_factors[1]*laplace_psi1s(argument1[m], root_neon[1][0]) + root_neon[2][1]*rooth_factors[2]*laplace_psi2s(argument1[m], root_neon[2][0]) + root_neon[3][1]*rooth_factors[3]*laplace_psi2s(argument1[m], root_neon[3][0]) + root_neon[4][1]*rooth_factors[4]*laplace_psi2s(argument1[m], root_neon[4][0]) + root_neon[5][1]*rooth_factors[5]*laplace_psi2s(argument1[m], root_neon[5][0]); 
		}
	if(charge==12){
		temp = root_magnesium[0][1]*rooth_factors[0]*laplace_psi1s(argument1[m], root_magnesium[0][0]) + root_magnesium[1][1]*rooth_factors[1]*laplace_psi3s(argument1[m], root_magnesium[1][0]) + root_magnesium[2][1]*rooth_factors[2]*laplace_psi3s(argument1[m], root_magnesium[2][0]) + root_magnesium[3][1]*rooth_factors[3]*laplace_psi3s(argument1[m], root_magnesium[3][0]) + root_magnesium[4][1]*rooth_factors[4]*laplace_psi3s(argument1[m], root_magnesium[4][0]) + root_magnesium[5][1]*rooth_factors[5]*laplace_psi3s(argument1[m], root_magnesium[5][0]) + root_magnesium[6][1]*rooth_factors[6]*laplace_psi3s(argument1[m], root_magnesium[6][0]) + root_magnesium[7][1]*rooth_factors[7]*laplace_psi3s(argument1[m], root_magnesium[7][0]);



		}
	if(charge==14){
		temp = root_silisium[0][1]*rooth_factors[0]*laplace_psi1s(argument1[m], root_silisium[0][0]) + root_silisium[1][1]*rooth_factors[1]*laplace_psi3s(argument1[m], root_silisium[1][0]) + root_silisium[2][1]*rooth_factors[2]*laplace_psi3s(argument1[m], root_silisium[2][0]) + root_silisium[3][1]*rooth_factors[3]*laplace_psi3s(argument1[m], root_silisium[3][0]) + root_silisium[4][1]*rooth_factors[4]*laplace_psi3s(argument1[m], root_silisium[4][0]) + root_silisium[5][1]*rooth_factors[5]*laplace_psi3s(argument1[m], root_silisium[5][0]) + root_silisium[6][1]*rooth_factors[6]*laplace_psi3s(argument1[m], root_silisium[6][0]) + root_silisium[7][1]*rooth_factors[7]*laplace_psi3s(argument1[m], root_silisium[7][0]);


		}


	
  }  
  if(n==1){
	if(charge == 4){
		temp = root_beryllium[0][2]*rooth_factors[0]*laplace_psi1s(argument1[m], root_beryllium[0][0]) + root_beryllium[1][2]*rooth_factors[1]*laplace_psi1s(argument1[m], root_beryllium[1][0]) + root_beryllium[2][2]*rooth_factors[2]*laplace_psi2s(argument1[m], root_beryllium[2][0]) + root_beryllium[3][2]*rooth_factors[3]*laplace_psi2s(argument1[m], root_beryllium[3][0]) + root_beryllium[4][2]*rooth_factors[4]*laplace_psi2s(argument1[m], root_beryllium[4][0]) + root_beryllium[5][2]*rooth_factors[5]*laplace_psi2s(argument1[m], root_beryllium[5][0]); 

		}
	if(charge==10){
		temp = root_neon[0][2]*rooth_factors[0]*laplace_psi1s(argument1[m], root_neon[0][0]) + root_neon[1][2]*rooth_factors[1]*laplace_psi1s(argument1[m], root_neon[1][0]) + root_neon[2][2]*rooth_factors[2]*laplace_psi2s(argument1[m], root_neon[2][0]) + root_neon[3][2]*rooth_factors[3]*laplace_psi2s(argument1[m], root_neon[3][0]) + root_neon[4][2]*rooth_factors[4]*laplace_psi2s(argument1[m], root_neon[4][0]) + root_neon[5][2]*rooth_factors[5]*laplace_psi2s(argument1[m], root_neon[5][0]);  
		
		}

	if(charge==12){
		temp = root_magnesium[0][2]*rooth_factors[0]*laplace_psi1s(argument1[m], root_magnesium[0][0]) + root_magnesium[1][2]*rooth_factors[1]*laplace_psi3s(argument1[m], root_magnesium[1][0]) + root_magnesium[2][2]*rooth_factors[2]*laplace_psi3s(argument1[m], root_magnesium[2][0]) + root_magnesium[3][2]*rooth_factors[3]*laplace_psi3s(argument1[m], root_magnesium[3][0]) + root_magnesium[4][2]*rooth_factors[4]*laplace_psi3s(argument1[m], root_magnesium[4][0]) + root_magnesium[5][2]*rooth_factors[5]*laplace_psi3s(argument1[m], root_magnesium[5][0]) + root_magnesium[6][2]*rooth_factors[6]*laplace_psi3s(argument1[m], root_magnesium[6][0]) + root_magnesium[7][2]*rooth_factors[7]*laplace_psi3s(argument1[m], root_magnesium[7][0]); 


		}

	if(charge==14){
		temp = root_silisium[0][2]*rooth_factors[0]*laplace_psi1s(argument1[m], root_silisium[0][0]) + root_silisium[1][2]*rooth_factors[1]*laplace_psi3s(argument1[m], root_silisium[1][0]) + root_silisium[2][2]*rooth_factors[2]*laplace_psi3s(argument1[m], root_silisium[2][0]) + root_silisium[3][2]*rooth_factors[3]*laplace_psi3s(argument1[m], root_silisium[3][0]) + root_silisium[4][2]*rooth_factors[4]*laplace_psi3s(argument1[m], root_silisium[4][0]) + root_silisium[5][2]*rooth_factors[5]*laplace_psi3s(argument1[m], root_silisium[5][0]) + root_silisium[6][2]*rooth_factors[6]*laplace_psi3s(argument1[m], root_silisium[6][0]) + root_silisium[7][2]*rooth_factors[7]*laplace_psi3s(argument1[m], root_silisium[7][0]);
		}

   
  } 
  if(n==2){

	if(charge==10){

		temp = root_neon[0][4]*rooth_factors[6]*(laplace_psi2p_x(argument1[m],root_neon[0][3],r[m][0])
      +laplace_psi2p_y(argument1[m], root_neon[0][3], r[m][0], r[m][1])+laplace_psi2p_y(argument1[m], root_neon[0][3], r[m][0],r[m][2])) + root_neon[1][4]*rooth_factors[7]*(laplace_psi2p_x(argument1[m],root_neon[1][3],r[m][0])
      +laplace_psi2p_y(argument1[m], root_neon[1][3], r[m][0], r[m][1])+laplace_psi2p_y(argument1[m], root_neon[1][3], r[m][0],r[m][2])) +  root_neon[2][4]*rooth_factors[8]*(laplace_psi2p_x(argument1[m],root_neon[2][3],r[m][0])
      +laplace_psi2p_y(argument1[m], root_neon[2][3], r[m][0], r[m][1])+laplace_psi2p_y(argument1[m], root_neon[2][3], r[m][0],r[m][2])) +  root_neon[3][4]*rooth_factors[9]*(laplace_psi2p_x(argument1[m],root_neon[3][3],r[m][0])
      +laplace_psi2p_y(argument1[m], root_neon[3][3], r[m][0], r[m][1])+laplace_psi2p_y(argument1[m], root_neon[3][3], r[m][0],r[m][2])); 

		}

	if(charge==12){
		temp =  root_magnesium[0][5]*rooth_factors[8]*(laplace_psi2p_x(argument1[m],root_magnesium[0][4],r[m][0])
      +laplace_psi2p_y(argument1[m], root_magnesium[0][4], r[m][0], r[m][1])+laplace_psi2p_y(argument1[m], root_magnesium[0][4], r[m][0],r[m][2])) + root_magnesium[1][5]*rooth_factors[9]*(laplace_psi4p_x(argument1[m],root_magnesium[1][4],r[m][0])
      +laplace_psi4p_y(argument1[m], root_magnesium[1][4], r[m][0], r[m][1])+laplace_psi4p_y(argument1[m], root_magnesium[1][4], r[m][0],r[m][2])) + root_magnesium[2][5]*rooth_factors[10]*(laplace_psi4p_x(argument1[m],root_magnesium[2][4],r[m][0])
      +laplace_psi4p_y(argument1[m], root_magnesium[2][4], r[m][0], r[m][1])+laplace_psi4p_y(argument1[m], root_magnesium[2][4], r[m][0],r[m][2])) + root_magnesium[3][5]*rooth_factors[11]*(laplace_psi4p_x(argument1[m],root_magnesium[3][4],r[m][0])
      +laplace_psi4p_y(argument1[m], root_magnesium[3][4], r[m][0], r[m][1])+laplace_psi4p_y(argument1[m], root_magnesium[3][4], r[m][0],r[m][2]))  + root_magnesium[4][5]*rooth_factors[12]*(laplace_psi4p_x(argument1[m],root_magnesium[4][4],r[m][0])
      +laplace_psi4p_y(argument1[m], root_magnesium[4][4], r[m][0], r[m][1])+laplace_psi4p_y(argument1[m], root_magnesium[4][4], r[m][0],r[m][2])); 


	}

	if(charge==14){

		temp =  root_silisium[0][5]*rooth_factors[8]*(laplace_psi2p_x(argument1[m],root_silisium[0][4],r[m][0])
      +laplace_psi2p_y(argument1[m], root_silisium[0][4], r[m][0], r[m][1])+laplace_psi2p_y(argument1[m], root_silisium[0][4], r[m][0],r[m][2])) + root_silisium[1][5]*rooth_factors[9]*(laplace_psi4p_x(argument1[m],root_silisium[1][4],r[m][0])
      +laplace_psi4p_y(argument1[m], root_silisium[1][4], r[m][0], r[m][1])+laplace_psi4p_y(argument1[m], root_silisium[1][4], r[m][0],r[m][2])) + root_silisium[2][5]*rooth_factors[10]*(laplace_psi4p_x(argument1[m],root_silisium[2][4],r[m][0])
      +laplace_psi4p_y(argument1[m], root_silisium[2][4], r[m][0], r[m][1])+laplace_psi4p_y(argument1[m], root_silisium[2][4], r[m][0],r[m][2])) + root_silisium[3][5]*rooth_factors[11]*(laplace_psi4p_x(argument1[m],root_silisium[3][4],r[m][0])
      +laplace_psi4p_y(argument1[m], root_silisium[3][4], r[m][0], r[m][1])+laplace_psi4p_y(argument1[m], root_silisium[3][4], r[m][0],r[m][2])) + root_silisium[4][5]*rooth_factors[12]*(laplace_psi4p_x(argument1[m],root_silisium[4][4],r[m][0])
      +laplace_psi4p_y(argument1[m], root_silisium[4][4], r[m][0], r[m][1])+laplace_psi4p_y(argument1[m], root_silisium[4][4], r[m][0],r[m][2])) + root_silisium[5][5]*rooth_factors[13]*(laplace_psi4p_x(argument1[m],root_silisium[5][4],r[m][0])
      +laplace_psi4p_y(argument1[m], root_silisium[5][4], r[m][0], r[m][1])+laplace_psi4p_y(argument1[m], root_silisium[5][4], r[m][0],r[m][2])) + root_silisium[6][5]*rooth_factors[14]*(laplace_psi4p_x(argument1[m],root_silisium[6][4],r[m][0])
      +laplace_psi4p_y(argument1[m], root_silisium[6][4], r[m][0], r[m][1])+laplace_psi4p_y(argument1[m], root_silisium[6][4], r[m][0],r[m][2])) + root_silisium[7][5]*rooth_factors[15]*(laplace_psi4p_x(argument1[m],root_silisium[7][4],r[m][0])
      +laplace_psi4p_y(argument1[m], root_silisium[7][4], r[m][0], r[m][1])+laplace_psi4p_y(argument1[m], root_silisium[7][4], r[m][0],r[m][2]));

	}
  }
  if(n==3){
	if(charge==10){

		temp = root_neon[0][4]*rooth_factors[6]*(laplace_psi2p_x(argument1[m],root_neon[0][3],r[m][1])
      +laplace_psi2p_y(argument1[m], root_neon[0][3], r[m][1], r[m][0])+laplace_psi2p_y(argument1[m], root_neon[0][3], r[m][1],r[m][2])) + root_neon[1][4]*rooth_factors[7]*(laplace_psi2p_x(argument1[m],root_neon[1][3],r[m][1])
      +laplace_psi2p_y(argument1[m], root_neon[1][3], r[m][1], r[m][0])+laplace_psi2p_y(argument1[m], root_neon[1][3], r[m][1],r[m][2])) +  root_neon[2][4]*rooth_factors[8]*(laplace_psi2p_x(argument1[m],root_neon[2][3],r[m][1])
      +laplace_psi2p_y(argument1[m], root_neon[2][3], r[m][1], r[m][0])+laplace_psi2p_y(argument1[m], root_neon[2][3], r[m][1],r[m][2])) +  root_neon[3][4]*rooth_factors[9]*(laplace_psi2p_x(argument1[m],root_neon[3][3],r[m][1])
      +laplace_psi2p_y(argument1[m], root_neon[3][3], r[m][1], r[m][0])+laplace_psi2p_y(argument1[m], root_neon[3][3], r[m][1],r[m][2])); 

		}
	if(charge==12){
		temp =  root_magnesium[0][5]*rooth_factors[8]*(laplace_psi2p_x(argument1[m],root_magnesium[0][4],r[m][1])
	+laplace_psi2p_y(argument1[m], root_magnesium[0][4], r[m][1], r[m][0])+laplace_psi2p_y(argument1[m], root_magnesium[0][4], r[m][1],r[m][2])) + root_magnesium[1][5]*rooth_factors[9]*(laplace_psi4p_x(argument1[m],root_magnesium[1][4],r[m][1])
	+laplace_psi4p_y(argument1[m], root_magnesium[1][4], r[m][1], r[m][0])+laplace_psi4p_y(argument1[m], root_magnesium[1][4], r[m][1],r[m][2])) + root_magnesium[2][5]*rooth_factors[10]*(laplace_psi4p_x(argument1[m],root_magnesium[2][4],r[m][1])
	+laplace_psi4p_y(argument1[m], root_magnesium[2][4], r[m][1], r[m][0])+laplace_psi4p_y(argument1[m], root_magnesium[2][4], r[m][1],r[m][2])) + root_magnesium[3][5]*rooth_factors[11]*(laplace_psi4p_x(argument1[m],root_magnesium[3][4],r[m][1])
	+laplace_psi4p_y(argument1[m], root_magnesium[3][4], r[m][1], r[m][0])+laplace_psi4p_y(argument1[m], root_magnesium[3][4], r[m][1],r[m][2]))  + root_magnesium[4][5]*rooth_factors[12]*(laplace_psi4p_x(argument1[m],root_magnesium[4][4],r[m][1])
	+laplace_psi4p_y(argument1[m], root_magnesium[4][4], r[m][1], r[m][0])+laplace_psi4p_y(argument1[m], root_magnesium[4][4], r[m][1],r[m][2])); 
	
	
		}

	if(charge==14){

		temp =  root_silisium[0][5]*rooth_factors[8]*(laplace_psi2p_x(argument1[m],root_silisium[0][4],r[m][1])
      +laplace_psi2p_y(argument1[m], root_silisium[0][4], r[m][1], r[m][0])+laplace_psi2p_y(argument1[m], root_silisium[0][4], r[m][1],r[m][2])) + root_silisium[1][5]*rooth_factors[9]*(laplace_psi4p_x(argument1[m],root_silisium[1][4],r[m][1])
      +laplace_psi4p_y(argument1[m], root_silisium[1][4], r[m][1], r[m][0])+laplace_psi4p_y(argument1[m], root_silisium[1][4], r[m][1],r[m][2])) + root_silisium[2][5]*rooth_factors[10]*(laplace_psi4p_x(argument1[m],root_silisium[2][4],r[m][1])
      +laplace_psi4p_y(argument1[m], root_silisium[2][4], r[m][1], r[m][0])+laplace_psi4p_y(argument1[m], root_silisium[2][4], r[m][1],r[m][2])) + root_silisium[3][5]*rooth_factors[11]*(laplace_psi4p_x(argument1[m],root_silisium[3][4],r[m][1])
      +laplace_psi4p_y(argument1[m], root_silisium[3][4], r[m][1], r[m][0])+laplace_psi4p_y(argument1[m], root_silisium[3][4], r[m][1],r[m][2])) + root_silisium[4][5]*rooth_factors[12]*(laplace_psi4p_x(argument1[m],root_silisium[4][4],r[m][1])
      +laplace_psi4p_y(argument1[m], root_silisium[4][4], r[m][1], r[m][0])+laplace_psi4p_y(argument1[m], root_silisium[4][4], r[m][1],r[m][2])) + root_silisium[5][5]*rooth_factors[13]*(laplace_psi4p_x(argument1[m],root_silisium[5][4],r[m][1])
      +laplace_psi4p_y(argument1[m], root_silisium[5][4], r[m][1], r[m][0])+laplace_psi4p_y(argument1[m], root_silisium[5][4], r[m][1],r[m][2])) + root_silisium[6][5]*rooth_factors[14]*(laplace_psi4p_x(argument1[m],root_silisium[6][4],r[m][1])
      +laplace_psi4p_y(argument1[m], root_silisium[6][4], r[m][1], r[m][0])+laplace_psi4p_y(argument1[m], root_silisium[6][4], r[m][1],r[m][2])) + root_silisium[7][5]*rooth_factors[15]*(laplace_psi4p_x(argument1[m],root_silisium[7][4],r[m][1])
      +laplace_psi4p_y(argument1[m], root_silisium[7][4], r[m][1], r[m][0])+laplace_psi4p_y(argument1[m], root_silisium[7][4], r[m][1],r[m][2]));

	}

   
  }
  if(n==4){
	if(charge==10){

		temp = root_neon[0][4]*rooth_factors[6]*(laplace_psi2p_x(argument1[m],root_neon[0][3],r[m][2])
      +laplace_psi2p_y(argument1[m], root_neon[0][3], r[m][2], r[m][0])+laplace_psi2p_y(argument1[m], root_neon[0][3], r[m][2],r[m][1])) + root_neon[1][4]*rooth_factors[7]*(laplace_psi2p_x(argument1[m],root_neon[1][3],r[m][2])
      +laplace_psi2p_y(argument1[m], root_neon[1][3], r[m][2], r[m][0])+laplace_psi2p_y(argument1[m], root_neon[1][3], r[m][2],r[m][1])) +  root_neon[2][4]*rooth_factors[8]*(laplace_psi2p_x(argument1[m],root_neon[2][3],r[m][2])
      +laplace_psi2p_y(argument1[m], root_neon[2][3], r[m][2], r[m][0])+laplace_psi2p_y(argument1[m], root_neon[2][3], r[m][2],r[m][1])) +  root_neon[3][4]*rooth_factors[9]*(laplace_psi2p_x(argument1[m],root_neon[3][3],r[m][2])
      +laplace_psi2p_y(argument1[m], root_neon[3][3], r[m][2], r[m][0])+laplace_psi2p_y(argument1[m], root_neon[3][3], r[m][2],r[m][1])); 

		}
	if(charge==12){
		temp =  root_magnesium[0][5]*rooth_factors[8]*(laplace_psi2p_x(argument1[m],root_magnesium[0][4],r[m][2])
	+laplace_psi2p_y(argument1[m], root_magnesium[0][4], r[m][2], r[m][0])+laplace_psi2p_y(argument1[m], root_magnesium[0][4], r[m][2],r[m][1])) + root_magnesium[1][5]*rooth_factors[9]*(laplace_psi4p_x(argument1[m],root_magnesium[1][4],r[m][2])
	+laplace_psi4p_y(argument1[m], root_magnesium[1][4], r[m][2], r[m][0])+laplace_psi4p_y(argument1[m], root_magnesium[1][4], r[m][2],r[m][1])) + root_magnesium[2][5]*rooth_factors[10]*(laplace_psi4p_x(argument1[m],root_magnesium[2][4],r[m][2])
	+laplace_psi4p_y(argument1[m], root_magnesium[2][4], r[m][2], r[m][0])+laplace_psi4p_y(argument1[m], root_magnesium[2][4], r[m][2],r[m][1])) + root_magnesium[3][5]*rooth_factors[11]*(laplace_psi4p_x(argument1[m],root_magnesium[3][4],r[m][2])
	+laplace_psi4p_y(argument1[m], root_magnesium[3][4], r[m][2], r[m][0])+laplace_psi4p_y(argument1[m], root_magnesium[3][4], r[m][2],r[m][1]))  + root_magnesium[4][5]*rooth_factors[12]*(laplace_psi4p_x(argument1[m],root_magnesium[4][4],r[m][2])
	+laplace_psi4p_y(argument1[m], root_magnesium[4][4], r[m][2], r[m][0])+laplace_psi4p_y(argument1[m], root_magnesium[4][4], r[m][2],r[m][1])); 
	
	
		}

	if(charge==14){

		temp =  root_silisium[0][5]*rooth_factors[8]*(laplace_psi2p_x(argument1[m],root_silisium[0][4],r[m][2])
      +laplace_psi2p_y(argument1[m], root_silisium[0][4], r[m][2], r[m][0])+laplace_psi2p_y(argument1[m], root_silisium[0][4], r[m][2],r[m][1])) + root_silisium[1][5]*rooth_factors[9]*(laplace_psi4p_x(argument1[m],root_silisium[1][4],r[m][2])
      +laplace_psi4p_y(argument1[m], root_silisium[1][4], r[m][2], r[m][0])+laplace_psi4p_y(argument1[m], root_silisium[1][4], r[m][2],r[m][1])) + root_silisium[2][5]*rooth_factors[10]*(laplace_psi4p_x(argument1[m],root_silisium[2][4],r[m][2])
      +laplace_psi4p_y(argument1[m], root_silisium[2][4], r[m][2], r[m][0])+laplace_psi4p_y(argument1[m], root_silisium[2][4], r[m][2],r[m][1])) + root_silisium[3][5]*rooth_factors[11]*(laplace_psi4p_x(argument1[m],root_silisium[3][4],r[m][2])
      +laplace_psi4p_y(argument1[m], root_silisium[3][4], r[m][2], r[m][0])+laplace_psi4p_y(argument1[m], root_silisium[3][4], r[m][2],r[m][1])) + root_silisium[4][5]*rooth_factors[12]*(laplace_psi4p_x(argument1[m],root_silisium[4][4],r[m][2])
      +laplace_psi4p_y(argument1[m], root_silisium[4][4], r[m][2], r[m][0])+laplace_psi4p_y(argument1[m], root_silisium[4][4], r[m][2],r[m][1])) + root_silisium[5][5]*rooth_factors[13]*(laplace_psi4p_x(argument1[m],root_silisium[5][4],r[m][2])
      +laplace_psi4p_y(argument1[m], root_silisium[5][4], r[m][2], r[m][0])+laplace_psi4p_y(argument1[m], root_silisium[5][4], r[m][2],r[m][1])) + root_silisium[6][5]*rooth_factors[14]*(laplace_psi4p_x(argument1[m],root_silisium[6][4],r[m][2])
      +laplace_psi4p_y(argument1[m], root_silisium[6][4], r[m][2], r[m][0])+laplace_psi4p_y(argument1[m], root_silisium[6][4], r[m][2],r[m][1])) + root_silisium[7][5]*rooth_factors[15]*(laplace_psi4p_x(argument1[m],root_silisium[7][4],r[m][2])
      +laplace_psi4p_y(argument1[m], root_silisium[7][4], r[m][2], r[m][0])+laplace_psi4p_y(argument1[m], root_silisium[7][4], r[m][2],r[m][1]));

	}

    
  }
  if(n==5){
	if(charge==12){
		temp = root_magnesium[0][3]*rooth_factors[0]*laplace_psi1s(argument1[m], root_magnesium[0][0]) + root_magnesium[1][3]*rooth_factors[1]*laplace_psi3s(argument1[m], root_magnesium[1][0]) + root_magnesium[2][3]*rooth_factors[2]*laplace_psi3s(argument1[m], root_magnesium[2][0]) + root_magnesium[3][3]*rooth_factors[3]*laplace_psi3s(argument1[m], root_magnesium[3][0]) + root_magnesium[4][3]*rooth_factors[4]*laplace_psi3s(argument1[m], root_magnesium[4][0]) + root_magnesium[5][3]*rooth_factors[5]*laplace_psi3s(argument1[m], root_magnesium[5][0]) + root_magnesium[6][3]*rooth_factors[6]*laplace_psi3s(argument1[m], root_magnesium[6][0]) + root_magnesium[7][3]*rooth_factors[7]*laplace_psi3s(argument1[m], root_magnesium[7][0]);

	}
	if(charge==14){
		temp = root_silisium[0][3]*rooth_factors[0]*laplace_psi1s(argument1[m], root_silisium[0][0]) + root_silisium[1][3]*rooth_factors[1]*laplace_psi3s(argument1[m], root_silisium[1][0]) + root_silisium[2][3]*rooth_factors[2]*laplace_psi3s(argument1[m], root_silisium[2][0]) + root_silisium[3][3]*rooth_factors[3]*laplace_psi3s(argument1[m], root_silisium[3][0]) + root_silisium[4][3]*rooth_factors[4]*laplace_psi3s(argument1[m], root_silisium[4][0]) + root_silisium[5][3]*rooth_factors[5]*laplace_psi3s(argument1[m], root_silisium[5][0]) + root_silisium[6][3]*rooth_factors[6]*laplace_psi3s(argument1[m], root_silisium[6][0]) + root_silisium[7][3]*rooth_factors[7]*laplace_psi3s(argument1[m], root_silisium[7][0]);

	}


   
  } //Her velger jeg m_l=0, dvs solid harmonic = y. aLtså x_i=1.
  //*/

  if(n==6){
    temp = root_silisium[0][6]*rooth_factors[8]*(laplace_psi2p_x(argument1[m], root_silisium[0][4], r[m][1]) 
           +laplace_psi2p_y(argument1[m], root_silisium[0][4], r[m][1], r[m][0]) + laplace_psi2p_y(argument1[m], root_silisium[0][4], r[m][1], r[m][2])) + root_silisium[1][6]*rooth_factors[9]*(laplace_psi4p_x(argument1[m], root_silisium[1][4], r[m][1]) 
           +laplace_psi4p_y(argument1[m], root_silisium[1][4], r[m][1], r[m][0]) + laplace_psi4p_y(argument1[m], root_silisium[1][4], r[m][1], r[m][2])) + root_silisium[2][6]*rooth_factors[10]*(laplace_psi4p_x(argument1[m], root_silisium[2][4], r[m][1]) 
           +laplace_psi4p_y(argument1[m], root_silisium[2][4], r[m][1], r[m][0]) + laplace_psi4p_y(argument1[m], root_silisium[2][4], r[m][1], r[m][2])) + root_silisium[3][6]*rooth_factors[11]*(laplace_psi4p_x(argument1[m], root_silisium[3][4], r[m][1]) 
           +laplace_psi4p_y(argument1[m], root_silisium[3][4], r[m][1], r[m][0]) + laplace_psi4p_y(argument1[m], root_silisium[3][4], r[m][1], r[m][2])) + root_silisium[4][6]*rooth_factors[12]*(laplace_psi4p_x(argument1[m], root_silisium[4][4], r[m][1]) 
           +laplace_psi4p_y(argument1[m], root_silisium[4][4], r[m][1], r[m][0]) + laplace_psi4p_y(argument1[m], root_silisium[4][4], r[m][1], r[m][2])) + root_silisium[5][6]*rooth_factors[13]*(laplace_psi4p_x(argument1[m], root_silisium[5][4], r[m][1]) 
           +laplace_psi4p_y(argument1[m], root_silisium[5][4], r[m][1], r[m][0]) + laplace_psi4p_y(argument1[m], root_silisium[5][4], r[m][1], r[m][2])) + root_silisium[6][6]*rooth_factors[14]*(laplace_psi4p_x(argument1[m], root_silisium[6][4], r[m][1]) 
           +laplace_psi4p_y(argument1[m], root_silisium[6][4], r[m][1], r[m][0]) + laplace_psi4p_y(argument1[m], root_silisium[6][4], r[m][1], r[m][2])) + root_silisium[7][6]*rooth_factors[15]*(laplace_psi4p_x(argument1[m], root_silisium[7][4], r[m][1]) 
           +laplace_psi4p_y(argument1[m], root_silisium[7][4], r[m][1], r[m][0]) + laplace_psi4p_y(argument1[m], root_silisium[7][4], r[m][1], r[m][2]));
		
	


  }
  delete [] argument1;
  delete [] argument2;
  return temp;
}  //end laplace_phi()

//////////////////////////////////////////////////
//
//JASTROW-FACTOR METHODS
//
//////////////////////////////////////////////////

//get jastrow-ratio in metropolis test
double atom::getJastrowRatio(double** distance_old, double** distance_new,double beta){
  double jastrowRatio=0;

  for(int k=0; k<no_of_particles; k++){
    for(int l=0; l<k; l++){
      jastrowRatio += distance_new[l][k]-distance_old[l][k];
    }
  }
  for(int k=0; k<no_of_particles; k++){
    for(int l=k+1; l<no_of_particles; l++){
      jastrowRatio += distance_new[l][k]-distance_old[l][k];
    }
  }

  return jastrowRatio;
}//end getRatio()

//Upper triangular matrix with distances
void atom::getDistance(double** distance, double** r_old, double beta){
  double a=0;
  double temp; 

  for(int k=0; k<no_of_particles; k++){
    for(int l=k+1; l<no_of_particles; l++){
      temp = diffR(r_old, k,l);
      //spin up
      if(((k < no_of_particles/2) && (l <no_of_particles/2)) || ((k>=no_of_particles/2 && l>=no_of_particles/2))){
	a=0.25;
	distance[k][l] = a*temp/(1+beta*temp);
      }
      //spin down
      else{
	a=0.5;
	distance[k][l] = a*temp/(1+beta*temp);	
      }
    }
  }
}//end getDistance()

//distance between two particles
double atom:: diffR(double** r, int k, int l){
  double r12;
  r12=0;
  for(int j=0; j<dimension; j++){
    r12 += (r[k][j]-r[l][j])*(r[k][j]-r[l][j]);
  }
  r12 = sqrt(r12);
  return r12;
}//end diffR()

//gradient for all the particles
double atom::gradient_jastrow(double** r,int p, int q, int d, double beta){
  double temp1,temp2,temp3;
  double a=0;
  if(((p < no_of_particles/2) && (q <no_of_particles/2)) || ((p>=no_of_particles/2 && q>=no_of_particles/2))){
    a=0.25;
    temp1 = diffR(r,p,q);
    temp2 = 1+beta*temp1;
    temp3 = r[p][d]-r[q][d];
    return a*temp3/temp1/temp2/temp2; 
  }
  else{
    a=0.5;
    temp1 = diffR(r,p,q);
    temp2 = 1+beta*temp1;
    temp3 = r[p][d]-r[q][d];
    return a*temp3/temp1/temp2/temp2; 
  }
}//end gradien_jastrow()




/*************************    OLD WAVEFUNCTION    ********************
 *
 *Need this to find the derivative of the wavefunction wrt 
 *variational parameters
 *
 */
double atom::wave_function(double**r,double alpha,double beta){/*********************/
  double wf, r_single_particle, r_ij;
  double *argument1,*argument2;
  wf=0;

  argument1 = new double[no_of_particles];
  argument2 = new double[no_of_particles];
  //  zero arguments
  for (int i=0; i<no_of_particles; i++){ 
    argument1[i] = 0.0;
    argument2[i] = 0.0;
  }   
  
  for(int i=0; i<no_of_particles; i++){
    r_single_particle=0;
    for(int j=0; j<dimension; j++){
      r_single_particle += r[i][j]*r[i][j];//r^2=x^2+y^2+z^2
    }
    argument1[i] += sqrt(r_single_particle);
  }
  
  //SLATER DETERMINANT MATRICES AND VARIABLES
  double **D_up, **D_down;
  D_up = create_matrix(no_of_particles/2, no_of_particles/2);
  D_down = create_matrix(no_of_particles/2, no_of_particles/2);

  /***** BRUTE FORCE SLATER DETERMINANT,***********/
  set_up_slater2(r, D_up, D_down, alpha, beta);
 
  wf = det(D_up,no_of_particles/2)*det(D_down,no_of_particles/2);//5x5 det
  
  for(int i=0; i<no_of_particles-1; i++){
    for(int j=i+1; j<no_of_particles; j++){
      r_ij=0;
      for(int k=0; k<dimension; k++){
	r_ij +=(r[i][k]-r[j][k])*(r[i][k]-r[j][k]);
      }
      r_ij=sqrt(r_ij);
      if(((i < 5) && (j < 5)) || ((i>=5 && j>=5))){
	wf *= exp(r_ij/(4*(1+beta*r_ij)));
      }
      else{
	wf *= exp(r_ij/(2*(1+beta*r_ij)));
      }
    }
  }
 
  delete_matrix(D_up,no_of_particles/2);
  delete_matrix(D_down,no_of_particles/2);

  delete [] argument1;
  delete [] argument2;
  return wf;
}//end wave_function()

//brute force neon wavefunction
void atom::set_up_slater2(double**r, double **D_up, double **D_down, double alpha, double beta){
  //SPIN UP
  for(int i=0; i<no_of_particles/2; i++){
    for(int j=0; j<no_of_particles/2; j++){
      D_up[i][j] = phi(r,alpha,beta,j,i);/********/
    }
  }
  //SPIN DOWN
  for(int i=0; i<no_of_particles/2; i++){
    for(int j=0; j<no_of_particles/2; j++){
      D_down[i][j]=phi(r,alpha,beta,j,i+no_of_particles/2); 
    }
  }
  //inverting matrices
 //  inverse(D_up, no_of_particles/2);
//   inverse(D_down, no_of_particles/2);
}//end set_up_slater2()

//calculates the determinant to matrix A
double atom::det(double** A, int dim) {

  if (dim == 2)
    return A[0][0]*A[1][1] - A[0][1]*A[1][0];

  double sum = 0;

  for (int i = 0; i < dim; i++) {

    double** sub = new double*[dim-1];
    for (int j = 0; j < i; j++)
      sub[j] = &A[j][1];
    for (int j = i+1; j < dim; j++)
      sub[j-1] = &A[j][1];

    if(i % 2 == 0)
      sum += A[i][0] * det(sub,dim-1);
    else
      sum -= A[i][0] * det(sub,dim-1);

    delete[] sub;
  }

  return sum;
}

// //sjekk
// void atom::quantum_force2(double** r, double **qm_force, double wf_old,double alpha, double beta){
//   double **r_plus,**r_minus;
//   double wf_minus, wf_plus;
//   r_plus =  create_matrix(no_of_particles, dimension);  
//   r_minus = create_matrix(no_of_particles, dimension);  
//   for(int i=0; i<no_of_particles; i++){
//     for(int j=0; j<dimension; j++){
//       r_plus[i][j]=r_minus[i][j]=r[i][j];
//     }
//   }

//   //QUANTUM FORCE F = 2/(psi)*del(psi)
//   for(int p=0; p<no_of_particles; p++){
//     for(int q=0; q<dimension; q++){
//       r_plus[p][q] =r[p][q]+h;
//       r_minus[p][q]=r[p][q]-h;
//       wf_plus  = wave_function(r_plus,alpha,beta);
//       wf_minus = wave_function(r_minus,alpha,beta);
//       qm_force[p][q] =(wf_plus-wf_minus)/(h*wf_old); 
//       r_plus[p][q] = r[p][q];
//       r_minus[p][q] = r[p][q];
//     }
//   }
//   //free memory
//   delete_matrix(r_plus,no_of_particles);//4 = number of particles
//   delete_matrix(r_minus,no_of_particles);
// }//end quantum_force()

// double atom::local_energy2(double **r, double wf_old,double alpha, double beta){
//   //local variables
//   double E_local, wf_minus, wf_plus, E_kinetic, 
//     E_potential, r_single_particle, r_12;
  
//   //matrices with positions
//   double **r_plus,**r_minus;
//   r_plus =  create_matrix(no_of_particles, dimension);  
//   r_minus = create_matrix(no_of_particles, dimension);  
  
//   for(int i=0; i<no_of_particles; i++){
//     for(int j=0; j<dimension; j++){
//       r_plus[i][j]=r_minus[i][j]=r[i][j];
//     }
//   }

//   //reset values
//   E_local    =0;
//   E_kinetic  =0;
//   E_potential=0;

//   /*
//     KINETIC ENERGY
//     E_kinetic(local) = -d^2(wf_old)/dx^2/wf_old
    
//   */

//   for(int p=0; p<no_of_particles; p++){
//     for(int q=0; q<dimension; q++){
//       r_plus[p][q] =r[p][q]+h;
//       r_minus[p][q]=r[p][q]-h;
//       wf_plus  = wave_function(r_plus,alpha,beta);
//       wf_minus = wave_function(r_minus,alpha,beta);
//       E_kinetic -=(wf_minus+wf_plus-2*wf_old); 
//       r_plus[p][q] = r[p][q];
//       r_minus[p][q] = r[p][q];
//     }
//   }
//   E_kinetic = 0.5*h2*E_kinetic/wf_old;

//   //POTENTIAL ENERGY
    
//   //Coulomb potential
//   for(int i=0; i<no_of_particles; i++){
//     r_single_particle=0;
//     for(int j=0; j<dimension; j++){
//       r_single_particle += r[i][j]*r[i][j];//r^2=x^2+y^2+z^2
//     }
//     E_potential -=charge/sqrt(r_single_particle); 
//   }
//   //ee-interaction
//   for(int i=0; i<no_of_particles-1; i++){
//     for(int j=i+1; j<no_of_particles; j++){
//       r_12=0;
//       for(int k=0; k<dimension; k++){
// 	r_12 +=(r[i][k]-r[j][k])*(r[i][k]-r[j][k]);
//       }
//       E_potential +=1/sqrt(r_12);
//     }
//   }
//   E_local = E_potential + E_kinetic;
  
//   //free memory
//   delete_matrix(r_plus,no_of_particles);//remember to update these if you change number of particles
//   delete_matrix(r_minus,no_of_particles);

//   return E_local;
// }//end local_energy()






void input(){
  alphaalpha=1;
  betabeta = 1.379;
}

int main(int argc, char* argv[]){
 
  //MPI-INITIALISE
 
outfile.open("Roothaan_results.dat");

  MPI_Init(NULL,NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  
  idum = -time(NULL) - myrank;
  if(myrank==0){
    input();
  }
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&alphaalpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&betabeta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);



	helium.mc_sampling(alphaalpha, betabeta, false);

	if(myrank==0){
		helium.main_e=helium.main_e/numprocs;
		helium.main_e2=helium.main_e2/numprocs;
		helium.variance =helium.main_e2 - helium.main_e*helium.main_e; 
		cout<<"!ENERGY! : " << helium.main_e << endl;
		cout<<"!VARIANCE! : " << helium.variance << endl;
	}
if(myrank==0){
	outfile << setiosflags(ios::showpoint | ios::uppercase);
	outfile << setw(15) << setprecision(8) << helium.main_e;
	outfile << setw(15) << setprecision(8) << helium.variance;

		
	outfile << endl;
	}
	



//*/

outfile.close();

  //cgm();
  //END MPI
  MPI_Finalize();
  
  return 0;
}