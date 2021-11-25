//====================================================================
// PURPOSE:  Solving Poisson Equation  
//           Geometry -- 1*1 square, it is divided into 4 sub-domain
//           Row by row
//           Boundary propagate -- the points by all processors
// LANGUAGE: C++
// COMPILER: xlC for IBM RS6000 PowerStations
// AUTHOR:   Chuanyi Ding  ding@arc.unm.edu
// DATE:     November 30, 1995
// ADDRESS:  Department of Chemical and Nuclear Engineering
//           The University of New Mexico
//           Albuquerque, NM 87131
//====================================================================

#include "poisson.h"

    const double pi = 4*atan(1);
    int N;   // number of random walks

    const int nxx = 200;     // statically allocated size of an array
    const int nyy = 200;
    double u[nxx+1][nyy+1], x[nxx+1], y[nyy+1];
    int onBoundary[nxx+1][nyy+1];

    int myRank, numProcs;
    double uSend, uRec;      // u for a cell
    int flag = 0;   // (= 1) means calculation done

// solving the Poisson equation with Dirichlet boundary condition
void poissonSolver(RandomStream* rand, int num_grid_x, int num_grid_y, 
                   int num_random) {

    N = num_random;
    const int nx = num_grid_x;      // the size of grids
    const int ny = num_grid_y;
    const double delta = 1. / nx;   // step of x and y

    // set rank for each processor
    MPI_Status status;

    MPI_Comm_rank( MPI_COMM_WORLD,&myRank );
    MPI_Comm_size( MPI_COMM_WORLD, &numProcs );

    // set time
    float startTime, initEndTime, innerEndTime, endTime;
    startTime = MPI_Wtime();

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
    
    // intialize flags,  0 = not calculated, 1 = inner boundary
    //                   2 = outer boundary
    for ( ix=0; ix<=nx; ix++ ) {
        for ( iy=0; iy<=ny; iy++ )
            onBoundary[ix][iy] = 0;
    }
 
    // set domain x(0,1) y(0,1)
    for ( ix=0; ix<=nx; ix++ )
        x[ix] = (1. * ix) / nx;
    for ( iy=0; iy<=ny; iy++ )
        y[iy] = (1. * iy) / ny;

    // set boundary condition
    boundary(nx, ny);

    initEndTime = MPI_Wtime();

    cout << "The initialization time (seconds) is " << initEndTime - startTime
         << " on processor " << myRank << endl;

    // we devide the initial region in to 4 or 9 sub-region if the initial
    // region (nx, ny) is large enough, then calculate the parallel and 
    // vertical boundary, finally evaluate on each sub-region one by one
    int numInner = 2;    // determine the number of sub-region, numInner**2

    int ixp, iyp;
    
    // evaluate the parallel inner boundary line starting from leftest
    // initialize the index of the cell
    ix = myRank + 1; 
    iy = ny / numInner;
    
    // calcalute u for (ix,iy)
    while ( 1 ) {
       // calculate u[ix][iy] for a point (ix,iy), and send it to others
       sendU( rand, ix, iy, delta ); 
       
       // receives u[iix][iiy] from the other processors
       for ( int i = 0; i < numProcs; i++) {
           if (i != myRank) {
              MPI_Recv(&uRec, 1, MPI_DOUBLE, i, 99, MPI_COMM_WORLD,&status);
              if (uRec <= -1000) flag = 1;
              else {
                 ixp = ix + (i - myRank); // the index of inner boundary, 
                 iyp = iy;                // u at this cell received from
                                       // another processor
                 u[ixp][iyp] = uRec;
                 onBoundary[ixp][iyp] = 1;
              }
           }
        }

        // Barrier synchronization across all group members
        MPI_Barrier( MPI_COMM_WORLD );

        // check if the whole calculation done, flag = 1
        if (flag == 1) break;

        // set ix, iy of next node for the processor
        ix += numProcs;
        assert( ((ix < nx+10) && (iy < ny+10)), " out of given domain");
    }

    flag = 0;   // (= 1) means calculation done
    // evaluate the vertical inner boundary line starting from lowest
    // initialize the index of the cell
    ix = nx / numInner;
    iy = myRank + 1;

    onBoundary[nx/2][ny/2] = 0;    // reset middle point for vertical
                                   // not stop at this point

    // calcalute u for (ix,iy)
    while ( 1 ) {

       // calculate u[ix][iy] for a point (ix,iy), and send it to others
       sendU( rand, ix, iy, delta );

       // receives u[iix][iiy] from the other processors
       for ( int i = 0; i < numProcs; i++) {
           if (i != myRank) {
              MPI_Recv(&uRec, 1, MPI_DOUBLE, i, 99, MPI_COMM_WORLD,&status);
              if (uRec <= -1000) flag = 1;
              else {
                 ixp = ix;                  // the index of inner boundary,
                 iyp = iy + (i - myRank);   // u at this cell received from
                                            // another processor
                 u[ixp][iyp] = uRec;
                 onBoundary[ixp][iyp] = 1;
              }
           }
       }

       // Barrier synchronization across all group members
       MPI_Barrier( MPI_COMM_WORLD );

       // check if the whole calculation done, flag = 1
       if (flag == 1) break;

       // set ix, iy of next node for the processor
       iy += numProcs;
       assert( ((ix < nx+10) && (iy < ny+10)), " out of given domain");
    }

    cout << " ********** inner boundary done ************ " << endl;

    innerEndTime = MPI_Wtime();
    cout << "The inner set time (seconds) is "<<innerEndTime-initEndTime
         << " on processor " << myRank << endl;
 
    // set sub-region
    int nx1[4], ny1[4], nx2[4], ny2[4];  // boundary of each region

    // initialize boundary for each region 
    nx1[0] = 0;        // x  at lowerleft 
    ny1[0] = 0;        // y  at lowerleft 
    nx2[0] = nx / 2;   // x  at upperright
    ny2[0] = ny / 2;   // y  at upperright 

    nx1[1] = nx / 2;
    ny1[1] = 0;
    nx2[1] = nx;
    ny2[1] = ny / 2;

    nx1[2] = 0;
    ny1[2] = ny / 2;
    nx2[2] = nx / 2;
    ny2[2] = ny;
    
    nx1[3] = nx / 2;
    ny1[3] = ny / 2;
    nx2[3] = nx;
    ny2[3] = ny;

    for (int iRegion = 0; iRegion < 4; iRegion++) {
       flag = 0;
       ix = nx1[iRegion] + 1 + myRank;
       iy = ny1[iRegion] + 1;
     
       // calcalute u for (ix,iy) of each region
       while ( 1 ) {

          // calculate u[ix][iy] for a point (ix,iy), and send it to others
          sendU( rand, ix, iy, delta );

          // receives u[iix][iiy] from the other processors
          for ( int i = 0; i < numProcs; i++) {
              if (i != myRank) {
                 MPI_Recv(&uRec, 1, MPI_DOUBLE, i, 99, MPI_COMM_WORLD,&status);
                 if (uRec <= -1000) flag = 1;
                 else {
                    ixp = ix + (i - myRank); // the index of inner boundary,
                    iyp = iy;                // u at this cell received from
                                             // another processor
                    if ( ixp > (nx2[iRegion] - 1) ) {
                       ixp = nx1[iRegion] + (ixp - nx2[iRegion]) + 1;
                       iyp = iy + 1;
                    }
                    u[ixp][iyp] = uRec;
                    onBoundary[ixp][iyp] = 1;
                 }
              }
          } 

          // Barrier synchronization across all group members
          MPI_Barrier( MPI_COMM_WORLD );

          // check if the whole calculation done, flag = 1
          if (flag == 1) break;

          // set ix, iy of next node for the processor
          ix += numProcs;
          if ( ix > (nx2[iRegion] - 1) ) {
             ix = nx1[iRegion] + (ix - nx2[iRegion]) + 1;
             iy++;
          }
          assert( ((ix < nx+10) && (iy < ny+10)), " out of given domain");
       }
       cout << " ====== region " << iRegion << " is done" << endl;
    }

    endTime = MPI_Wtime();

    cout << "The running time (seconds) is " << endTime - startTime
         << " on processor " << myRank << endl;

    cout << "========  done ========= " << myRank << endl;
 

/*    // chech if the data being communicated
    if ( myRank == 0 ) {
       for ( iy=1; iy<=ny; iy+=5 ) {
           for ( ix=1; ix<=nx; ix+=5 )
               cout << "ix " << ix << " iy " << iy << " u " << u[ix][iy] 
                    << endl;
       }
    }
*/

}



// set outer boundary condition
void boundary(int nx, int ny) {
    for ( int ix=0; ix<=nx; ix++ ) {
        u[ix][0] = 0;    // lower
        onBoundary[ix][0] = 2;
        u[ix][ny] = 0;   // upper
        onBoundary[ix][ny] = 2;
    }

    for ( int iy=0; iy<=ny; iy++ ) {
        u[0][iy] = 0;    // left
        onBoundary[0][iy] = 2;
        u[nx][iy] = 0;   // right
        onBoundary[nx][iy] = 2;
    }
}



// evaluate f(x,y)
double f(double x, double y) {
    return ( 2. * pi * pi * sin(pi*x) * sin(pi*y) );
} 


// evaluate u by Monte Carlo method
double solveU(RandomStream* rand, int ix, int iy, double delta) {

    double uxy = 0, sqU = 0;

    for ( int i=1; i<=N; i++ ) {   // N random walks
        int iix = ix;
        int iiy = iy;
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

            if ( onBoundary[iix][iiy] > 0 ) break;
            fxy += f( x[iix], y[iiy] );
         }
         double tempU = u[iix][iiy] + fxy * delta * delta / 4;
 	 uxy += tempU;
	 sqU += tempU * tempU;
    }
    uxy /= N;
    sqU /= N;
    if (ix == iy)
       cout << "ix = " << ix << "iy = " << iy << "u = " << uxy 
	    << " sigma = " << sqrt( (sqU - uxy * uxy)/N ) << endl;
    return uxy;
} 


// calculate u[ix][iy] for a point (ix,iy), and send it to others
void sendU( RandomStream* rand, int ix, int iy, double delta ) {
   if ( onBoundary[ix][iy] > 0 ) {   // u[x][y] evaluated before
      uSend = -1000;
      flag = 1;
   }
   else {
      uSend = solveU(rand, ix, iy, delta);
      u[ix][iy] = uSend;
      onBoundary[ix][iy] = 1;   // set this node as boundary
   }

   // sends u[ix][iy] to the other processor iSend
   for ( int i = 0; i < numProcs; i++) {
       if (i != myRank)
       MPI_Send(&uSend, 1, MPI_DOUBLE, i, 99, MPI_COMM_WORLD);
   }
}


