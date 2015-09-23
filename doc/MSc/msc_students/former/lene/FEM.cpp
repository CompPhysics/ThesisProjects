#include <mpi.h>
#include "FEM.h"
#include <lpp/lapack.hh>




FEM::FEM(int elements, int local, int int_N_in, double step_length, int comm, double x_min){
//sette opp flere muligheter her? -parallell/ikke parallel, annend input, eks rmax+punkter
  //cout << "FEM konstruktør" << endl;
  M=elements;
  n_e=2; //test med 3?, funksjon for varierende antall
  int_N=4;//tester, har bare satt opp noen eksakte
  N=M*(n_e-1)+1;
  //l=N;
  h=step_length; 
  this->x_min=x_min;

  MPI_Comm_rank ( comm, &my_rank );
  MPI_Comm_size( comm, &P);
  this->comm=comm;
  
  mat=matrix(N,N); //bruk annen type her
  nymat=band_matrix(N,n_e);
  vec=matrix(N,N); //annet navn

  set_integration_points();

  potential=matrix(N,N); //bytte til array?
  x=matrix(N,N); //unødvendig?

  //sendetest --> henter forrige sin og neste sin - ok
//   int nabo=my_rank;
//   MPI_Status status;
//   MPI_Send(&my_rank, 1, MPI_DOUBLE, (my_rank+1+P)%P, 100, MPI_COMM_WORLD);

//   MPI_Recv(&nabo, 1, MPI_DOUBLE, (my_rank-1+P)%P, 100, MPI_COMM_WORLD, &status);
//   cout << "sender til : " << (my_rank+1+P)%P << ", mottar fra " << (my_rank-1+P)%P << endl;
//   cout << my_rank << "-"<< nabo << endl;

//   MPI_Send(&my_rank, 1, MPI_DOUBLE, (my_rank-1+P)%P, 100, MPI_COMM_WORLD);

//   MPI_Recv(&nabo, 1, MPI_DOUBLE, (my_rank+1+P)%P, 100, MPI_COMM_WORLD, &status);
//   cout << "sender til : " << (my_rank-1+P)%P << ", mottar fra " << (my_rank+1+P)%P << endl;
//   cout << my_rank << "-"<< nabo << endl;  
}

void FEM::test1D(double r_max){ 
  double boundary=r_max;
  h=2.0*boundary/M;
  for(int i=0; i<M; i++){ //bruker senere fra 1!
    for(int j=0; j<int_N; j++){
      for(int s=0; s<n_e; s++){
	x[i][j]+=N_i(s, xi[j])*(-boundary+(s+i)*h);//nodene
      }
      potential[i][j]=x[i][j]*x[i][j];
    }
  }
  solve(); 

}

void FEM::set_potential(int dim, int qn, double omega2, double K){ 
  double x, x2, nabla;

  if(dim==2) nabla=qn*qn-0.25;
  else if(dim==3) nabla=qn*(qn+1);
  else{ nabla=0; cout << "feil potensial!"<<endl;}

  for(int i=0; i<M; i++){ //bruker senere fra 1!
    for(int j=0; j<int_N; j++){
      x=0;
      
      for(int s=0; s<n_e; s++){
	x+=N_i(s, xi[j])*(x_min+(s+i)*h);//nodene-->sett opp tidligere/flere muligheter, regner ut dette flere ganger med samme
	//if(j==0) cout <<x_min+ (s+i)*h << endl; //ser riktig ut, skriver ut alle til den ene først?, fordi lite testsett? gjøres mange ganger
      }
      //cout << x <<endl;
      x2=x*x;
      potential[i][j]=nabla/x2 + K/x+omega2*x2;
     }
  }
}





void FEM::bc(){ //forbedre, trenger bare at ta til samme element --> bare på ytre rand!
//   cout << "P="<< P << endl;
//   cout << "me="<< my_rank << endl;
  if(my_rank==0){
    //cout << "rand1" <<endl;
    for(int i=0; i<n_e; i++){
      mat[i][0]=0;
      vec[i][0]=0;
      mat[0][i]=0;
      vec[0][i]=0;
    }
    mat[0][0]=vec[0][0]=1;
  }
  if(my_rank==P-1){
    //cout << "rand2" <<endl;
    for(int i=N-1; i>=N-n_e; i--){
      mat[i][N-1]=0;
      vec[i][N-1]=0;
      mat[N-1][i]=0;
      vec[N-1][i]=0;
    }
    mat[N-1][N-1]=vec[N-1][N-1]=1;
  }

  

}


double FEM::N_i(int i, double x){
  //n_e=2: 
  if(i==0) return 0.5-0.5*x;
  else if(i==1) return 0.5+0.5*x;
  
}
double FEM::dN_i(int i, double x){ //ta med konstant her? pass på utenfor!
  //n_e=2:
  if(i==0) return -0.5;
  else if(i==1) return 0.5;
}


void FEM::make_system(){
  matrix element_left=matrix(n_e,n_e);
  matrix element_right=matrix(n_e,n_e);

  for(int e=1; e<=M; e++){
    calc_element(e, element_left, element_right);  
    for(int i=0; i<n_e; i++){ //fyller global --> mer generell?
      for(int j=0; j<n_e; j++){
	mat[(n_e-1)*(e-1)+i][(n_e-1)*(e-1)+j]+=element_left[i][j];
	nymat[(n_e-1)*(e-1)+i][(n_e-1)*(e-1)+j]+=element_left[i][j];
	vec[(n_e-1)*(e-1)+i][(n_e-1)*(e-1)+j]+=element_right[i][j];
      }
    }
  }
 
//     mat.print();
//     cout << endl;
//    vec.print();
  
  //bc();
  cout << endl;

  //skriv ut i riktig rekkefølge
//   int flag;
//   MPI_Status status;
//   if(my_rank>0)
//     MPI_Recv(&flag, 1, MPI_INT, my_rank-1, 100, MPI_COMM_WORLD, &status);

//   cout << "Linknigssystem for " << my_rank << endl;
//   mat.print();
//   cout << endl;
//   nymat.print();

//   if (my_rank<P-1)
//     MPI_Send(&my_rank, 1, MPI_INT, my_rank+1, 100, MPI_COMM_WORLD);


  mvtest();

   //delete her?
}
void FEM::solve_system(){
  //cout << "solve system" << my_rank<< endl;
  //bytte til arpack -> pass på egenvektorer
  sol=Eigen(N-2,2,h); //0->N!!!!
  double* testvec=new double[2];
  testvec[0]=0;
  testvec[1]=0;
  double* eig=new double[N];//n?
  for(int i=0; i<N; i++) 
    eig[i]=0;
  gen_eigen(mat, vec, eig); //bytte ut med arpack!
    for(int i=0; i<N-2; i++) {
      //cout << 0.5*eig[i+2] << endl;
      sol.eigenvalue[i]=0.5*eig[i+2]; 
  }
  delete[] testvec;
  //cout << "ferdigsolve system" << my_rank<< endl;
}

double FEM::integrate(double* func){ 
  double sum=0;
  for(int i=0; i<int_N; i++){
    sum+=weights[i]*func[i];
  }
  return sum;
}

void FEM::calc_element(int e, matrix &left, matrix &right){ 
  //initialisere annet sted? Denne funksjonen kalles mange ganger? husk å slette!
  double* int1=new double[int_N];
  double* int2=new double[int_N];
  double M_ij;

  for(int i=0; i<n_e; i++){
    for(int j=0; j<n_e; j++){
      for(int k=0; k<int_N; k++){
	M_ij=h/2.0*N_i(i,xi[k])*N_i(j,xi[k]);
	int1[k]=2.0/h*dN_i(i,xi[k])*dN_i(j,xi[k])+M_ij*potential[e-1][k];
	int2[k]=M_ij;
      }
      left[i][j]=integrate(int1);
      right[i][j]=integrate(int2);
    }
  }

  delete[] int1;
  delete[] int2;
}

void FEM::gen_eigen(matrix A, matrix B, double* d){

  long int N=A.rows;
  long int M=A.cols;
  long int info;
  char uplo = 'l';
  long int type=1;

  double* B_array=B.get_array();
  double* A_array=A.get_array();

  double* e=new double[N-1];
  double* tau=new double[N-1];

  lpp::potrf (&uplo, &N, B_array, &M, &info);
  lpp::sygst (&type, &uplo, &N, A_array, &M, B_array, &M, &info);
  lpp::sytrd (&uplo, &N, A_array, &M, d, e, tau, &info);
  lpp::sterf (&N, d, e, &info);

}


void FEM::set_integration_points(){
  weights=new double[int_N];
  xi=new double[int_N];
 
  //http://mathworld.wolfram.com/Legendre-GaussQuadrature.html
  if(int_N==2){
    xi[0] = -sqrt(1.0/3.0);
    weights[0]= 1;
    xi[1]=-xi[0];
    weights[1]=1;    
  }
  else if(int_N==4){ 
    //gauss-legendre for polynomial of 4th order
    xi[0]= -sqrt(525+70*sqrt(30))/35.0;
    weights[0]= (18-sqrt(30))/36.0;
    xi[1]=-sqrt(525-70*sqrt(30))/35.0;
    weights[1]=(18+sqrt(30))/36.0;
    xi[2]=-xi[1];
    weights[2]=weights[1];
    xi[3]=-xi[0];
    weights[3]=weights[0];
  }else if(int_N==5){
    //gauss-legendre for polynomial of 5th order
    xi[0]=-sqrt(245+14*sqrt(70))/21.0;
    weights[0]=(322.0-13*sqrt(70))/900.0;
    xi[1]=-sqrt(245-14*sqrt(70))/21.0;
    weights[1]=(322.0+13*sqrt(70))/900.0;
    xi[2]=0;
    weights[2]=128.0/225.0;
    xi[3]=-xi[1];
    weights[3]=weights[1];
    xi[4]=-xi[0];
    weights[4]=weights[0];

  }else if(int_N%2==1){ //oddetall(trengs for simpsons?)
    //simpsons
    double step=2.0/(int_N-1);
    for(int i=0; i<int_N; i++) {
      xi[i]=-1+step*i;
      weights[i]=(2+2*(i%2))*step/3.0;
    }
    weights[0]*=0.5;
    weights[int_N-1]*=0.5;
  }else{
    //trapezoidal
    double step=2.0/(int_N-1);
    for(int i=0; i<int_N; i++) {
      xi[i]=-1+step*i;
      weights[i]=step;
    }
    weights[0]*=0.5;
    weights[int_N-1]*=0.5;
  }

}

void FEM::mvtest(){
  //cout << "tester matvec:" << endl;
  double start=1;
  double end;
  double prev, next;
  MPI_Status status;
int flag;



//   if(my_rank>0)
//     MPI_Recv(&start, 1, MPI_DOUBLE, my_rank-1, 2, MPI_COMM_WORLD, &status);
//   //cout << my_rank<<" [" << start << endl;
 
//   end=start+N-1;
//   //cout << my_rank << " ]" << end << endl;


//   if (my_rank<P-1)
//     MPI_Send(&end, 1, MPI_DOUBLE, my_rank+1, 2, MPI_COMM_WORLD);

//   double* test1 = new double[N];
//   double* test2 = new double[N];
//   for(int i=0; i<N; i++) {
//     test1[i]=start+i;
//     test2[i]= 0.23;
//   }

//   matvec(nymat, test1, test2);

//   if(my_rank>0)
//     MPI_Recv(&flag, 1, MPI_INT, my_rank-1, 100, MPI_COMM_WORLD, &status);

//   for(int i=0; i<N-1; i++) cout << test2[i] << endl;
//   if(my_rank==P-1) cout << test2[N-1] << endl;

//   if (my_rank<P-1)
//     MPI_Send(&my_rank, 1, MPI_INT, my_rank+1, 100, MPI_COMM_WORLD);

  call_P_ARPACK();

}

void FEM::matvec (matrix matrix, double* in, double* out){
  matrix.matvec(in, out);
  matvec_comm(out);
}

void FEM::matvec (band_matrix matrix, double* in, double* out){
  matrix.matvec(in, out);
  matvec_comm(out);
}

void FEM::matvec_comm(double* out){

  double prev, next;
  MPI_Status status;
  if (my_rank>0)
    MPI_Send(&out[0], 1, MPI_DOUBLE, my_rank-1, my_rank-1, MPI_COMM_WORLD);

  if (my_rank<P-1){
    MPI_Recv(&next, 1, MPI_DOUBLE, my_rank+1, my_rank, MPI_COMM_WORLD, &status);
    out[N-1]+=next;
  }

  if (my_rank<P-1)
    MPI_Send(&out[N-1], 1, MPI_DOUBLE, my_rank+1, my_rank+1, MPI_COMM_WORLD);

  if(my_rank>0){
    MPI_Recv(&prev, 1, MPI_DOUBLE, my_rank-1, my_rank, MPI_COMM_WORLD, &status);
    out[0]=prev;
  }
}


extern "C" { 
  double test_( double& );
  void pssaupd_( int& , int& , char, int&, char*, int& ,
		  double&, double*, int&, double**, int&, int*,
		  int*, double*, double*, int&, int& );
  void psnaupd_( int& , int& , char, int&, char*, int& ,
		  double&, double*, int&, double**, int&, int*,
		  int*, double*, double*, int&, int& );
  //void testmore_( const char*, const int&, const char*, const int& ); 
  //void lene_( double&);
}

void FEM::call_P_ARPACK(){

   double d1, d2;
   
   d1 = 10.0;
   d2 = test_( d1 ); // test is an external Fortran function 
   cout << "fortrantest: " << d1 << " " << d2 << endl;

   //int comm = MPI_COMM_WORLD;
   //cout << "comm2" << comm << endl; // ikke samme må sende med!
  int ido=0;
  char bmat = 'G';
  int nmax=256;
  int nloc=20;
  //if(my_rank<P-1) nloc--;
  char* which= "LM";
  int nev=1;
  double tol=0;
  double* resid=new double[nmax+1];
  int ncv=4;
  int ldv=nmax;
  double** v= new double*[nmax+1];
  for(int i=1; i<=nmax; i++) v[i]=new double[nmax+1];
  int* iparam=new int[11+1];
  int* ipntr=new int[14+1];
  double* workd= new double[3*nmax+1];
  int lworkl=3*nmax*nmax+6*nmax;
  double* workl= new double[lworkl+1];
  int info;
  iparam[1]=1;
  iparam[3]=1;
  iparam[7]=2;
  cout << "kaller pssaupd, nloc=" << nloc <<endl;

  psnaupd_( comm, ido, bmat, nloc, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info );
  cout << "info = " << info << endl;
  cout << "ido = " << ido << endl;
  
}
