//====================================================================
// PURPOSE:  Solving Poisson Equation  
//           Geometry -- 1*1 square, a 0.5*0.5 hole
//           Outside to inside
//           Boundary propagate -- the points by all processors
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
                    int num_grid_x, int num_grid_y, int num_random );
double f(double x, double y);
void grid( int&, int&, int, int );

    const double pi = 4*atan(1);
    int n1, n2, n3, n4;   // boundary position 
    int point = 1;        // where the calculated point is, initializing
                          // 1 means calculation starts from dowm side
//    int ix, iy;           // index of the position

int main(int argc, char* argv[]) {
    int myRank, numProcs;
    MPI_Status status;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank( MPI_COMM_WORLD, &myRank );
    MPI_Comm_size( MPI_COMM_WORLD, &numProcs );

    // Gather the information needed to parameterize the simulation.
/*    cout <<
    "Specify positive random number seed (or zero for random starting place): ";
    int32 seed;
    cin >> seed;

    cout << "Specify how many random walks for a point: ";
    int num_random;
    cin >> num_random;

    cout << "Specify how many grid in x direction ";
    int num_grid_x;
    cin >> num_grid_x;

    cout << "Specify how many grid in y direction ";
    int num_grid_y;
    cin >> num_grid_y;
*/
    int seed = 1;
    int num_random = 1000;
    int num_grid_x = 100;
    int num_grid_y = 100;
   
    RandomStream* rand;
// Define random number generator
    if ( seed == 0 )       // random start
       rand = new RandomStream(RANDOMIZE);
    else
       rand = new RandomStream( seed );

    float startTime, endTime;
    startTime = MPI_Wtime();
    poissonSolver( rand, numProcs, myRank, num_grid_x, num_grid_y, num_random );
    endTime = MPI_Wtime();

    cout << "The running time (seconds) is " << endTime - startTime 
         << "on processor " << myRank << endl;

    MPI_Finalize();

    return 0;
}


// solving the Poisson equation with Dirichlet boundary condition
void poissonSolver( RandomStream* rand, int numProcs, int myRank,
                    int num_grid_x, int num_grid_y, int N ) {

    const int nx = num_grid_x;      // the size of grids
    const int ny = num_grid_y; 
    const double delta = 1. / nx;   // step of x and y

    // initializing the outer grids to be calculated 
    n1 = 1;
    n2 = nx - 1;
    n3 = ny - 1;
    n4 = 1;

    // statically allocated size of an array
    const int nxx = 100; 
    const int nyy = 100;
    double u[nxx+1][nyy+1], x[nxx+1], y[nyy+1];
    int onBoundary[nxx+1][nyy+1]; 

/*  dynamical allocation for an array
    double* u; 
    u = new double[nx+1][ny+1];
    assert(u, "Unable to allocate memory!\n");
    double* x;
    x = new double[nx+1];
    assert(x, "Unable to allocate memory!\n");
    double* y;
    y = new double[ny+1];
    assert(y, "Unable to allocate memory!\n");
    int* onBoundary;
    onBoundary = new int[nx+1][ny+1];
    assert(onBoundary, "Unable to allocate memory!\n");
*/

    int ix, iy;
    int tmpX, tmpY;
    for ( ix=0; ix<=nx; ix++ ) {
        for ( iy=0; iy<=ny; iy++ )
            onBoundary[ix][iy] = 0;
    }
 
// set domain x(0,1) y(0,1)
    for ( ix=0; ix<=nx; ix++ )
        x[ix] = (1. * ix) / nx;
    for ( iy=0; iy<=ny; iy++ )
        y[iy] = (1. * iy) / ny;

// set outerboundary condition
    for ( ix=0; ix<=nx; ix++ ) {
        u[ix][0] = 0;    // lower
        onBoundary[ix][0] = 1;
        u[ix][ny] = 0;   // upper
        onBoundary[ix][ny] = 1;
    }

    for ( iy=0; iy<=ny; iy++ ) {
        u[0][iy] = 0;    // left 
        onBoundary[0][iy] = 1;
        u[nx][iy] = 0;   // right
        onBoundary[nx][iy] = 1;
    } 

// set inner boundary condition
    int inx1 = nx / 4;
    int inx2 = 3 * nx / 4;
    int iny1 = ny / 4;
    int iny2 = 3 * ny / 4;

    for ( ix=inx1; ix<=inx2; ix++ ) {
        u[ix][iny1] = 0;    // lower
        onBoundary[ix][iny1] = 1;
        u[ix][iny2] = 0;   // upper
        onBoundary[ix][iny2] = 1;
    }

    for ( iy=iny1; iy<=iny2; iy++ ) {
        u[inx1][iy] = 0;    // left
        onBoundary[inx1][iy] = 1;
        u[inx2][iy] = 0;   // right
        onBoundary[inx2][iy] = 1;
    }

    // initialize ix, iy 
    ix = myRank + 1; 
    iy = 1;
  
    int iix, iiy; 
    int flag = 0;   // (= 1) means calculation done
    double myU;
    while( 1 ) {  

       // calculate u[ix][iy] for a point (ix,iy) 
       double uxy = 0, sqU = 0; 

       if (onBoundary[ix][iy] == 1) {   // u[x][y] evaluated before
          myU = -1000;
          flag = 1;
       }
       else {
         for ( int i=1; i<=N; i++ ) {   // N random walks
             iix = ix;
             iiy = iy;
             double fxy = f( x[iix], y[iiy] );

             while ( 1 ) {              // while loop for one random walk 
                float ran = rand->next();
                if ( ran >= 0.75 ) {
                   iix++;
                }
                else if ( ran >= 0.50) {
                   iiy++; 
                }
                else if ( ran >= 0.25) {
                   iix--;   
                }
                else {
                   iiy--;  
                }
                if ( onBoundary[iix][iiy] == 1 ) break;
                fxy += f( x[iix], y[iiy] );
             }
             double tempU = u[iix][iiy] + fxy * delta * delta / 4;
	     uxy += tempU;
	     sqU += tempU * tempU;
          }
          myU = uxy / N;
	  sqU /= N;
          u[ix][iy] = myU;
          onBoundary[ix][iy] = 1;   // set this node as boundary
          
	  // output
	  if ( ix == iy )
	     cout <<"ix = " << ix << " iy = " << iy << " u = " << u[ix][iy] 
                  << " sigma = " << sqrt( (sqU - myU * myU)/N ) << endl;
       }  
         
       // sends u[ix][iy] to the other processors
       for ( int i = 0; i < numProcs; i++) {
          if (i != myRank) 
             MPI_Send(&myU, 1, MPI_DOUBLE, i, 99, MPI_COMM_WORLD);
       }

       // receives u[iix][iiy] from the other processors
       MPI_Status status;
       double u0;
       for ( i = 0; i < numProcs; i++) {
          if (i != myRank) { 
             MPI_Recv(&u0, 1, MPI_DOUBLE, i, 99, MPI_COMM_WORLD,&status);
             if (u0 <= -1000) {
                flag = 1;
             }
             else {   
                tmpX = ix;
                tmpY = iy;  
                grid( tmpX, tmpY, i - myRank, 0 );
                iix = tmpX;
             	iiy = tmpY;
                  
                u[iix][iiy] = u0;
                onBoundary[iix][iiy] = 1;
             }
          }  
       }


       // Barrier synchronization across all group members
       MPI_Barrier( MPI_COMM_WORLD );

       // check if the whole calculation done, flag = 1
       if (flag == 1) break;

       // set ix, iy of next node for the processor
       tmpX = ix;
       tmpY = iy;
       grid( tmpX, tmpY, numProcs, 1 );
       ix = tmpX;
       iy = tmpY;

       assert( ((ix < nx) && (iy < ny)), " out of given domain");
    }
    cout << "========  done ========= " << myRank << endl;
 
/*
    // chech if the data being communicated
    if (myRank == 1) {
      for ( iy=0; iy<=ny; iy++ ) {
        for ( ix=0; ix<=nx; ix++ )
            cout << "ix " << ix << " iy " << iy << " u " << u[ix][iy] 
                 << endl;
      }
    }
*/

}


// evaluate f(x,y)
double f(double x, double y) {
    return ( 2. * pi * pi * sin(pi*x) * sin(pi*y) );
} 


//
void grid(int& ix, int& iy, int numProcs, int flag) {

//    enum tokens { down, right, up, right};

//     cout << "NN" << n1 << n2 <<n3 << n4 << endl;
     switch ( point ) {
           case 1 :
		ix += numProcs;
		if ( (n2-n4) <= (numProcs+1) ) {
                   if ( ix > n2 ) {
                      ix = (n4 - 1) + (ix - n2);
                      iy++;
                   }
		}
		else {
		   if ( ix > n2 ) {
                      iy = n1 + (ix - n2);
      		      ix = n2;
                      if (flag == 1) {
                         n1++;
		         point = 2;
                      }
                   }
 		}
		break;
           case 2 :
                iy += numProcs;
                if ( iy > n3 ) {
                   ix = n2 - (iy -n3);
                   iy = n3;
                   if (flag == 1) {
		      n2--;
                      point = 3;
                   }
                }
                break;
           case 3 :
                ix -= numProcs;
                if ( ix < n4 ) {
                   iy = n3 + (ix - n4);
                   ix = n4;
                   if (flag == 1) {
		      n3--;
		      point = 4;
		   }
                }
                break;
           case 4 :
                iy -= numProcs;
                if ( iy < n1 ) {
                   ix = n4 - (iy - n1);
                   iy = n1;
                   if (flag == 1) {
                      n4++;
		      point = 1;
		   }
                }
     }
}

