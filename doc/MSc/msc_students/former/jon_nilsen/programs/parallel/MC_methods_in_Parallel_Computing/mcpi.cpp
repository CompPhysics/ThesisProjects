//====================================================================
// PURPOSE:  Evaluates the expected value
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

void mcpi( int numProcs, int myRank, int num_random, 
                   int num_dimension);
int inBoundary( double[], int );
double f( double[], int );
double pdf( double[], int );

int main( int argc, char* argv[] ) {
    int myRank, numProcs;
    MPI_Status status;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank( MPI_COMM_WORLD,&myRank );
    MPI_Comm_size( MPI_COMM_WORLD, &numProcs );

    // Gather the information needed to parameterize the simulation.
    
    int num_random;
    int num_dimension = 2;

    if (myRank == 0)  
       num_random = 100000;

    // sends num_random to all processors
    MPI_Bcast(&num_random, 1, MPI_INT, 0, MPI_COMM_WORLD);


    float startTime, endTime;
    startTime = MPI_Wtime();
    mcpi( numProcs, myRank, num_random, num_dimension );
    endTime = MPI_Wtime();

    cout << "The running time (seconds) is " << endTime - startTime 
         << " on processor " << myRank << endl;

    MPI_Finalize();

    return 0;
}

// mcpi 
void mcpi( int numProcs, int myRank, int num_random,
                   int num_dimension) {

    const double PI25DT = 3.141592653589793238462643;
    double pi = 0;
    double x, myPi = 0;      
    int count = 0;   // counts the number for each processor
    int totalCount = 0;
    int rankOfProc;
    int seed;
    MPI_Status status;

    RandomStream* rand[10];
    double ran[10];

    for ( int i=0; i<num_dimension; i++) {	// create r.n. generator
       // the seeds are different for the processors
       seed = (myRank + 1) * 100 + i * 10 + 1;	
       if ( seed == 0 )       			
          rand[i] = new RandomStream(RANDOMIZE);
       else
          rand[i] = new RandomStream( seed );
    }
 
    for ( i=myRank; i<num_random; i+=numProcs ) { // N random walks

//       if (myRank == 0) {

//          for ( int i0=0; i0<numProcs; i0++ ) {

             for ( int i1=0; i1<num_dimension; i1++) {
		ran[i1] = rand[i1]->next();
	     }

//           if ( i0 != (numProcs-1) ) 	// also produce ran[i] for self
//              MPI_Send(&ran,num_dimension,MPI_DOUBLE,i0+1,99,MPI_COMM_WORLD);
//          }
//       }
//       else {
//          MPI_Recv(&ran,num_dimension,MPI_DOUBLE,0,99,MPI_COMM_WORLD,&status);
//       } 
                         
       // each processor calculates a part of integration
       if ( inBoundary( ran, num_dimension ) ) {
          myPi += f( ran, num_dimension ) * pdf( ran, num_dimension );
          count++;
       }
    }
    
    myPi = 4 * myPi / num_random;

    MPI_Reduce(&myPi, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (myRank == 0) {
       cout << "pi is = " << pi 
            << " sigma = " << 4 * sqrt( (pi/4 - (pi/4)*(pi/4)) / num_random) 
            << " Error = " << abs(pi - PI25DT) << endl;
    }

    cout << " count of this processor = " << count << endl;
}


// chech if the point is in the boubdaries
int inBoundary( double x[], int num_dimension ) {

   double dd = 0;
   for ( int i=0; i<num_dimension; i++ ) {
       dd += x[i] * x[i];
   }

   if ( dd < 1 ) 
      return 1;

   return 0;
}
            


// evaluate f(x)
double f( double x[], int num_dimension ) {
   return 1;
}


// evaluate pdf
double pdf( double x[], int num_dimension ) {
    return 1;
} 
