//====================================================================
// PURPOSE:  Solving Poisson Equation  
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

void poissonSolver( RandomStream* rand, int numProcs, int myRank,
                    int num_random );
double EvalPoint( double xx, double yy, int side );
double MinDist( double xx, double yy, int& side );
 
    const double pi = 4*atan(1);
 
int main(int argc, char* argv[]) {
    int myRank, numProcs;
    MPI_Status status;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank( MPI_COMM_WORLD,&myRank );
    MPI_Comm_size( MPI_COMM_WORLD, &numProcs );

    // Gather the information needed to parameterize the simulation.
/*    cout <<
    "Specify positive random number seed (or zero for random starting place): ";
    int32 seed;
    cin >> seed;

    cout << "Specify how many random walks for a point: ";
    int num_random;
    cin >> num_random;
*/
    // the seeds are different for the processors
    int seed = (myRank + 1) * 10 + 1;
    int num_random = 1000;
   
    RandomStream* rand;
// Define random number generator
    if ( seed == 0 )       // random start
       rand = new RandomStream(RANDOMIZE);
    else
       rand = new RandomStream( seed );

    float startTime, endTime;
    startTime = MPI_Wtime();
    poissonSolver( rand, numProcs, myRank, num_random );
    endTime = MPI_Wtime();

    cout << "The running time (seconds) is " << endTime - startTime 
         << " on processor " << myRank << endl;

    MPI_Finalize();

    return 0;
}


// solving the Poisson equation with Dirichlet boundary condition
void poissonSolver( RandomStream* rand, int numProcs, int myRank, 
						      int num_random ) {

    double x[11], y[11]; 	// The numbers of the points
    int numPoint = 1;		// only middle point

    for (int i=0; i<numPoint; i++) {
        x[i] = (i + 1.) / ( numPoint + 1. );
        y[i] = (i + 1.) / ( numPoint + 1. );
    }
 
//    int ip = myRank;		// initializing the picked point index  
    int side = 0; 		// to which 
    double xx, yy, myU, u = 0;
    double tempU, mysqU, sqU = 0;
    double minDist, acceptDist = 0.002;
   
    for ( i=0; i<numPoint; i++ ) {
        
        // each processor will pick up a point for evaluation
//        if ( ip < numPoint ) {

           myU = 0;	// initializing u
	   mysqU = 0;
           // N random walks for each point
           for ( int ir=myRank; ir<num_random; ir+=numProcs ) {
 
              // initilaizing the position of the point for each random walk
              xx = x[i];
              yy = y[i];
 
              while (1) {

                 // calculating the minimum distance of selected point
                 // to the boundary
                 minDist = MinDist(xx, yy, side);
                 
                 if ( minDist < acceptDist ) break;
    
                 // each processor generates a random number, and make a random
                 // move, the point moves randomly to any point at the circle
                 // with a radius of minDist, and centered at (xx, yy)
                 double ran = rand->next();
 
                 xx += minDist * sin(2*pi*ran);
                 yy += minDist * cos(2*pi*ran);
              }
 
              // sum
              tempU = EvalPoint( xx, yy, side );
	      myU += tempU;
	      mysqU += tempU * tempU;
            }  // end of a ramdom walk

            myU /= num_random;
	    mysqU /= num_random;

	    MPI_Reduce(&myU, &u, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
            MPI_Reduce(&mysqU, &sqU, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
 
           // print out
	   if ( myRank == 0 )
              cout << "At the point of " << x[0] << ",  " << y[0] 
		   << " , u is = " << u 
		   << " , sigma = " << sqrt( (sqU - u*u) / num_random)
		   << endl;
//        }
//	ip += numProcs;
    } // end for all points
}       


// evaluates the minimum distance from the boundary
// here boundary is a square from (0,0) to (1,1)
double MinDist(double xx, double yy, int& side ) {
    // side = 1 for left, 2 for right, 3 for up, 4 for down 
    double minDist = xx;
    side = 1;		// left

    if ( (1-xx) < minDist ) {
       minDist = 1 - xx;
       side = 2; 	// right 
    }

    if ( (1-yy) < minDist ) {
       minDist = 1 - yy;
       side = 3;     	// up 
    }

    if ( yy < minDist ) {
       minDist = yy;
       side = 4;     	// down 
    }
 
    return minDist;
}

// evaluate u at the boundary
// u at boundary is sin(2 pi x) at the up and down sides
// 		    sin(2 pi y) at the left and right sides
double EvalPoint( double xx, double yy, int side ) {
   double u;

   if ( (side == 1) || (side == 2) ) {
      u = sin( pi * yy );
   }
   else if ( (side == 3) || (side == 4) ) {
      u = sin( pi * xx );
   }
   else 
      cout << "error for side" << endl;

   return u;
}
 
