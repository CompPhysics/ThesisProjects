/* Variational Monte Carlo for the hydrogen and helium atoms */
/* Note that numerical derivation is employed here           */
/* The cusp condition is not implemented                     */
/* For the hydrogen atom choose 1 for number of particles    */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/* declaraton of functions */

/* Function to calculate random numbers (from lib.c) */
double ran1(long *);

/* Function to declear matrix (from lib.c) */
void **matrix(int, int, int);

/* Function to read in data from screen  */
void initialise(/*int *,*/ int *, int *, int *, int *, int *, double *, 
		double *);

/* The Mc sampling for the variational Monte Carlo */
void mc_sampling(int, int, int, int, int, int, double, double, 
		  double *, double *, double *);

/* The variational wave function */
double wave_function(double **, double, double, int, int);

/* Function to find distance between particle i and particle j */
double distance(double **, double *, double, int, int, int);

/* The local energy */
double local_energy(double **, double, double, double, double, int, int);

/* prints to screen the results of the calculations  */
void output(int, int, double, double, double *, double *);

/* Begin of main program   */

int main(){

    int number_blocks, number_cycles, max_variations, thermalization;
    int dimension, number_particles; 
    double initial_correlation, delta_correlation, lowest_correlation;
    double *cumulative_e, *cumulative_e2;
    dimension = 3;
    /*   Read in data */
    initialise(/*&dimension,*/ &number_blocks, &number_particles, &max_variations, 
	       &number_cycles, &thermalization, &initial_correlation, 
	       &delta_correlation) ;

    cumulative_e = (double *) malloc((max_variations+1)*sizeof(double));
    cumulative_e2 = (double *) malloc((max_variations+1)*sizeof(double));

    /*  Do the mc sampling  */
    mc_sampling(dimension, number_blocks, number_particles, max_variations, 
		thermalization, number_cycles, initial_correlation,
		delta_correlation, &lowest_correlation, cumulative_e, cumulative_e2);


    /* Print out results  */
    output(max_variations, number_cycles, delta_correlation, lowest_correlation,  
	   cumulative_e, cumulative_e2);
    return 0;
}


/* Monte Carlo sampling with the Metropolis algorithm  */

void mc_sampling(int dimension, int number_blocks, int number_particles, 
		 int max_variations, int thermalization, int number_cycles, 
		 double initial_correlation, 
		 double delta_correlation, double *lowest_correlation, 
		 double *cumulative_e, double *cumulative_e2){

  int blocks, cycles, variate, accept, i, j;
  long idum;
  double wfnew, wfold, correlation5, energy, energy2, delta_e, lowest_energy;
  double correlation, **r_old, **r_new;

  /* constants needed in calculations: */
  double l, kinetic_factor, density, scale, four_epsilon;
  density = 0.3648*0.9;  //Schiff-Verlet in sigma units
  l = pow((number_particles/density),0.3333333);  //Simulation cube
  scale = l/4;
  kinetic_factor = 9.27647;  //5*(hbar**2/m)/sigma**2
  four_epsilon = 40.88;

  correlation = initial_correlation-delta_correlation;
  idum=-1;
  /* allocate matrices which contain the position of the particles */ 
  r_old = (double **) matrix( number_particles, dimension, sizeof(**r_old));
  r_new = (double **) matrix( number_particles, dimension, sizeof(**r_new));
  for (i = 0; i < number_particles; i++) { 
    for ( j=0; j < dimension; j++) {
      r_old[i][j] = r_new[i][j] = 0;
    }
  }
  for (i = 0; i < max_variations+1; i++){
    cumulative_e[i] = cumulative_e2[i] = 0;
  }
  lowest_energy = 0.;
  /* loop over variational parameters */ 
  for (variate=1; variate <= max_variations; variate++){
    /* initialisations of variational parameters and energies */
    correlation += delta_correlation;
    correlation5 = pow(correlation,5);
    energy = energy2 = 0; accept =0; delta_e=0;
    /*  initial trial position */
    for (i = 0; i < number_particles; i++) { 
      for ( j=0; j < dimension; j++) {
	r_old[i][j] = l*ran1(&idum);
      }
    }
    wfold = wave_function(r_old, correlation5, l, dimension, number_particles);
    /* loop over blocks */
    for (blocks = 0; blocks <= number_blocks; blocks++){
      /* loop over monte carlo cycles */
      energy=0;
      for (cycles = 1; cycles <= number_cycles+thermalization; cycles++){ 
	/* new position */
	for (i = 0; i < number_particles; i++) { 
	  for ( j=0; j < dimension; j++) {
	    /* attempted move */
	    r_new[i][j] = r_old[i][j]+scale*(ran1(&idum)-0.5);
	    /* move particle to th fundamental cell */
	    if (r_new[i][j]<0.)
	      r_new[i][j] += l;
	    if (r_new[i][j]>0.)
	      r_new[i][j] -= l;
	  }
	}
	wfnew = wave_function(r_new, correlation5, l, dimension, number_particles); 
	/* Metropolis test */
        printf("wfnew/wfold %12.5E%12.5E\n",  wfnew, wfold);
	if(ran1(&idum) <= wfnew/wfold) { 
	  for (i = 0; i < number_particles; i++) { 
	    for ( j=0; j < dimension; j++) {
	      r_old[i][j]=r_new[i][j];
	    }
	  }
	  wfold = wfnew;
	  accept = accept+1;
	}
	
	if ( cycles > thermalization ) {
	  /* compute local energy */ 
//	  delta_e = local_energy(r_old, correlation5, l, kinetic_factor, four_epsilon, dimension,number_particles);
	  /* update energies */ 
	  energy += delta_e;
	  //energy2 += delta_e*delta_e;
	} /* end if() */
      }   /* end of loop over MC trials  */ 
      energy /= number_cycles;
      printf("Block %d energy %12.5E\n", blocks, energy);
      if(blocks!=0){ // Neglect first block
	cumulative_e[variate] += energy;
	cumulative_e2[variate] += energy*energy;
      }
    } /* end of loop over blocks */
      printf("variational parameter %12.5E accepted steps %d\n", 
	     correlation, accept);
      /* update the energy average and its squared */
      cumulative_e[variate] /= number_blocks*number_particles;
      cumulative_e2[variate] /= number_blocks*number_particles*number_particles;
    
  }    /* end of loop over variational  steps */
}   /* end mc_sampling function  */


/* Function to compute the squared wave function  */

double  wave_function(double **r, double correlation5, double l, int dimension, int number_particles){

  int i, j;
  double wf, dist, lh, argument, *r_distance;
  
  r_distance = (double *) malloc((dimension)*sizeof(double));
  lh = l/2;
  argument=0;
  for(i=0; i<number_particles; i++){
    for(j=0; j<number_particles; j++){
      if(i != j){ // never test same
	dist = distance(r, r_distance, l, i, j, dimension);
	if(dist > lh) {dist = lh;} // cutoff at l/2

	argument += 1./pow(dist,5.);      

      } // end if

    } // end for(j)
  } // end for(i)
  argument *= -correlation5;

  wf = exp(argument);

  return wf; // wf**2, that is
}

double distance(double **r, double *r_distance, double l, int i, int j, int dimension){
  int dim;
  double dist;
  dist = 0;
  for(dim=0; dim < dimension; dim++){
    r_distance[dim] = abs(r[i][dim]-r[j][dim]); // always positive
    if(r_distance[dim]>l/2) {r_distance[dim] = l-r_distance[dim];}

    dist += r_distance[dim]*r_distance[dim];
        printf("wfnew/wfold %12.5E%12.5E\n",  dist, r_distance[dim]);
  }
 return sqrt(dist);
}

/* Function to calculate the local energy */

double local_energy(double **r, double correlation5, double l, double kinetic_factor, 
		    double four_epsilon, int dimension, int number_particles){

  int particle_i, particle_j, i, j, k;
  double local_energy, dist, dist6, lh, z1, z2, z3, *r_distance;
  local_energy = dist = 0;
  r_distance = (double *) malloc((dimension)*sizeof(double));
  lh = l/2;
  for(particle_i=0; particle_i<number_particles-1; particle_i++){
    for(particle_j=particle_i+1; particle_j<number_particles; particle_j++){
      dist = distance(r, r_distance, l, particle_i, particle_j, dimension);
      if(dist<=lh) // cutoff correlation
	local_energy += kinetic_factor/pow(dist,7);
      /* check image of x, y and z axis */
//      for(i=-1; i<=0; i++){
//	z1=(r_distance[0]+i*l)*(r_distance[0]+i*l);
//	for(j=-1; j<=0; j++){ 
//	  z2=(r_distance[1]+j*l)*(r_distance[1]+j*l);
//	  for(k=-1; k<=0; k++){
//	    z3=(r_distance[2]+k*l)*(r_distance[2]+k*l);
//	    dist=z1+z2+z3;
//	    dist6=dist*dist*dist;
//	    if(dist<=l*l) // cutoff potential at L
//	      local_energy += four_epsilon*(1./dist6/dist6-1./dist6);
//	  }
//	}
//      }
    }
  }
  return local_energy;
}

void initialise(/*int *dimension, */int *number_blocks, int *number_particles,
                int *max_variations, int *number_cycles, 
                int *thermalization, double *initial_correlation, 
		double *delta_correlation) {

  //  printf("dimensionality = ");
  //  scanf("%d",dimension);
  printf("number of blocks = ");
  scanf("%d",number_blocks);
  printf("number of particles = ");
  scanf("%d",number_particles);
  printf("maximum variational parameters = ");
  scanf("%d",max_variations);
  printf("# MC steps= ");
  scanf("%d", number_cycles);
  printf("# Thermalization  steps= ");
  scanf("%d", thermalization);
  printf("initial correlation= ");
  scanf("%lf", initial_correlation);
  printf("delta correlation= ");
  scanf("%lf", delta_correlation);

}  /* end of function initialise   */



void output(int max_variations, int number_cycles, double delta_correlation, 
	    double lowest_correlation, double *cumulative_e, double *cumulative_e2){

  int i;
  double correlation, variance, error;
  FILE *output_file1;
  /* open output file */
  output_file1 = fopen("out_test.dat", "w") ;
  correlation = 0;
  for( i=1; i <= max_variations; i++){
    correlation += delta_correlation;  
    variance = cumulative_e2[i]-cumulative_e[i]*cumulative_e[i];
    error = sqrt(variance/number_cycles);
    fprintf(output_file1, "%12.5E %12.5E %12.5E %12.5E \n", correlation,
	    cumulative_e[i],variance, error );
  }
  
  fclose (output_file1);
  
}  /* end of function output */        


