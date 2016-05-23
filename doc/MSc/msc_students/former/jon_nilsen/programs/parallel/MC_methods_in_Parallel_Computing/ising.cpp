//====================================================================
// PURPOSE:  2-D Ising Model
// LANGUAGE: C++
// COMPILER: xlC for IBM RS6000 PowerStations
// AUTHOR:   Chuanyi Ding  ding@arc.unm.edu
// DATE:     November 30, 1995
// ADDRESS:  Department of Chemical and Nuclear Engineering
//           The University of New Mexico
//           Albuquerque, NM 87131
//====================================================================

#include <iostream.h>
#include <math.h>
#include "assert.h"
#include "random_numbers.h"
#include "mpi.h"

void ising( RandomStream* rand1, RandomStream* rand2, int numProcs,
            int myRank, int num_random);
void makeMove( int l, double h, double v, double ran1, double ran2, 
	       RandomStream* rand1, double tau);
int Np( int n, int l );
int Nm( int n, int l );
void initialize( int l, double h, double v );

    int lat[11][11];
    double m, e;

int main( int argc, char* argv[] ) {
    int myRank, numProcs;
    MPI_Status status;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank( MPI_COMM_WORLD,&myRank );
    MPI_Comm_size( MPI_COMM_WORLD, &numProcs );

    // Gather the information needed to parameterize the simulation.
/*
    cout <<
    "Specify positive random number seed (or zero for random starting place): ";
    int32 seed;
    cin >> seed;

    cout << "Specify how many random walks : ";
    int num_random;
    cin >> num_random;
*/
    int seed1 = 1;
    int seed2 = 100;
    int num_random;

    if (myRank == 0)
       num_random = 4000;

    // sends num_random to all processors
    MPI_Bcast(&num_random, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Define two random number generators, rand1 handling which processor
    // for a random walk, rand2 handling the history of a random walk
    RandomStream *rand1, *rand2;
    if ( seed1 == 0 )       // random start
       rand1 = new RandomStream(RANDOMIZE);
    else
       rand1 = new RandomStream( seed1 );

    if ( seed2 == 0 )       // random start
       rand2 = new RandomStream(RANDOMIZE);
    else
       rand2 = new RandomStream( seed2 );

    float startTime, endTime;
    startTime = MPI_Wtime();
    ising( rand1, rand2, numProcs, myRank, num_random);
    endTime = MPI_Wtime();

    cout << "The running time (seconds) is " << endTime - startTime
         << " on processor " << myRank << endl;

    MPI_Finalize();

    return 0;
}


//
void ising( RandomStream* rand1, RandomStream* rand2, int numProcs,
            int myRank, int num_random ) {
   
    const int l = 10;
    const double h = 0.1;
    const double v = 0.5;

    const double tauMin = 1;
    const double tauMax = 10; 
    const int num_runs = 10;
    const double delTau = (tauMax - tauMin) / (num_runs - 1);

    // initialize
    initialize( l, h, v );

    double sumM, sumM2, sumE, sumE2;
    double ran1, ran2;

    for (int i=myRank; i<num_runs; i += numProcs) {

        double tau = tauMin + i * delTau;

        for (int j=0; j<100*l*l; j++ ) {
	    ran1 = rand1->next();
	    ran2 = rand1->next();
	    makeMove( l, h, v, ran1, ran2, rand1, tau );
        }

        sumM = 0;
	sumM2 = 0;
	sumE = 0;
	sumE2 = 0;

	for ( int i1=0; i1<num_random; i1++ ) {
	    ran1 = rand1->next();
	    ran2 = rand1->next();

	    makeMove( l, h, v, ran1, ran2, rand1, tau );
	    sumM += abs(m);
	    sumM2 += m * m;
	    sumE += e;
	    sumE2 += e * e;
     	}

	double aveM = sumM / num_random;
        double aveM2 = sumM2 / num_random;
        double aveE = sumE / num_random;
        double aveE2 = sumE2 / num_random;

	double delSqM = aveM2 - aveM * aveM;
        double delSqE = aveE2 - aveE * aveE;

        aveM /= l * l;
        aveE /= l * l;

        double temp = tau * tau * l * l; 
        double chi = delSqM / temp;
        double specHeat = delSqE / temp;

	cout << tau << "  " << aveM << "  " << aveE << "  "
 	     << chi << "  " << specHeat << "  " 
	     << sqrt( delSqM/num_random ) / temp << "  " 
	     << sqrt( delSqE/num_random ) / temp 
	     << endl;
    }

}



//
void makeMove( int l, double h, double v, double ran1, double ran2, 
		RandomStream* rand1, double tau ) {
   // selecting a point at i = 1,l; j = 1,l
   int i = int(1 + l * ran1);
   int j = int(1 + l * ran2);

   int nm = Nm(j, l);
   int np = Np(j, l);
 
   double delE = 2 * lat[i][j] * ( h + v * (lat[i][nm] + lat[i][np]) );

   if ( (delE < 0) || (rand1->next() < exp(-delE/tau)) ) {
      lat[i][j] *= -1;
      m += 2 * lat[i][j];
      e += delE;
   } 
} 

int Np( int n, int l ) {

   int np;
   if ( n == l )  
      np = 1;
   else
      np = n + 1;

   return np;
}
 

//
int Nm( int n, int l ) {

   int nm;
   if ( n == 1 )
      nm = l;
   else
      nm = n - 1;

   return nm;

}


// 
void initialize( int l, double h, double v ) {
   
   int iran = 99991;

   for ( int i=1; i<=l; i++ ) {
       for ( int j=1; j<=l; j++ ) {
   	   lat[i][j] = 1;
	}
   }

   m = l * l;
   e = -(h + 2 * v) * m;

}
	  
 
