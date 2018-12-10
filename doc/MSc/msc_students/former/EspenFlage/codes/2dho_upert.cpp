#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std; 
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

// initialize
double genUnpertEigenfuncInteract(int n, int m, double r);
double genPerturbation(double x);
double genLaguerre(int n, int m, double x);
double numeric_integration( int n, int n_mark,int m,double left, double right, int points);
double normalization(int n, int m,double left, double right, int points);
void tqli(vector<double> &d, vector<double> &e, vector<vector<double> > &z);
void tred2(vector<vector<double> > &a, vector<double> &d, vector<double> &e);
void sort_tridiag(vector<double> &d);

// main
int main (int argc, const char* argv[]) {

  if (argc!=8) {
    printf("Error in syntax. Use 'program n n_mark m r_max sampling_points' where n,m,r_max and sampling_points are integers\n");
    exit(1);
  }
  else {

    // define parameters
    static int N=atoi(argv[1]);
    static int M=atoi(argv[2]);
    static int n=atoi(argv[3]);
    static int m=atoi(argv[4]);
    static int n_mark_max=atoi(argv[5]);
    static double r_max=atof(argv[6]);
    static int number_of_step=atoi(argv[7]);
    static double r_min=0.000000001;


    // calculate 2. order perturbation
    int pert_maks=10;
    double pert_second=0.0;
    for (int i=0;i<pert_maks;i++) {
      if (i!=n) {
	pert_second+=(fabs(numeric_integration(n,i,m,r_min,r_max,number_of_step))*fabs(numeric_integration(n,i,m,r_min,r_max,number_of_step)))/(2*(n-i));
      }
    }
    // print results 
    printf("---------------------------------------------------------------\n");
    printf("Eigenvalues using perturbation theory for N=%d, M=%d, n=%d and m=%d\n",N,M,n,m);
    printf("\n");
    printf("E_0 = 2(N+n+1)+|M|+|m| = %d\n",2*(N+n+1)+abs(M)+abs(m));
    printf("\n"); 
    printf("E_1 = <n|H_1|n> = <%d|H_1|%d> = %f\n",n,n,numeric_integration(n,n,m,r_min,r_max,number_of_step));
    printf("\n");
    printf("E_2 = sum_m((|<m|H_1|n>|^2)/E_0-E_m) = %f\n",pert_second);
    printf("\n");
    printf("E = E_0+E_1+E_2 = %f\n",2*(N+n+1)+abs(M)+abs(m)+numeric_integration(n,n,m,r_min,r_max,number_of_step)+pert_second);
    printf("---------------------------------------------------------------\n");

    // calculate eigenvalues and put them on file eigen in matrix form
    ofstream outFile("eigen");
    vector<vector<double> > eigen((n_mark_max+1),vector<double>(n_mark_max+1));
    vector<double> d(n_mark_max);
    vector<double> e(n_mark_max);
    for (int i=1; i<=n_mark_max;i++) {
      for (int j=1; j<=n_mark_max;j++) {
	if (i==j) {
	  eigen[i][j]=(2*(N+(i-1)+1)+abs(m)+abs(M))+numeric_integration(i-1,j-1,m,r_min,r_max,number_of_step);
	}
	else {
	  eigen[i][j]=numeric_integration(i-1,j-1,m,r_min,r_max,number_of_step);
	}
	outFile << eigen[i][j];
	outFile.put(' ');					     
      }
      outFile << endl;
    }
    outFile.close();
    
    
    tred2(eigen,d,e);
    tqli(d,e,eigen);
    sort_tridiag(d);
    
    // print eigenvalues
    
    for (int i=1;i<=n_mark_max;i++) {
      printf("eigenvalue [n=%d]=%f\n",i-1,d[i]);
    }
    
  }
}


// generate the unperturbed eigenfunctions

double genUnpertEigenfuncInteract(int n, int m, double r) {
  return(pow(r,abs(m))*exp(-(r*r)/2.0)*genLaguerre(n,abs(m),r*r));

}


// generate perturbation

double genPerturbation(double x) {
    return((0.1*sqrt(2.0))/x);
}

// generates the Laguerre polynomials used in mass-center therm

double genLaguerre(int n,int m, double x) {

  vector<double> laguerpol(n+1);

  // calculate laguerre polynomial terms using recursive iterations

  if (m<0) {
    printf("Error in M. Use M>=0, terminating\n");
    exit ( 1 );
  }

  laguerpol[0] = 1.0;
  laguerpol[1]=(double)(m+1)-x;

  for (int i=2;i<=n;i++) {
    laguerpol[i]=(((double)(2*i+m-1)-x)*laguerpol[i-1]+(double)(-i-m+1)*laguerpol[i-2])/(double)i;
  }
  
  return(laguerpol[n]);
}

// numerical integration using trapeziodal rule

double numeric_integration(int n, int n_mark,  int m, double left, double right,int points) {

  double left_value,right_value,x,step;
  double sum;
  double norm_n=0.0;
  double norm_n_mark=0.0;
  step=(right-left)/((double)points);

  // normalize wf for
  // n
  norm_n=normalization(n,m,left,right,points);
  
  // then n'
  norm_n_mark=normalization(n_mark,m,left,right,points);


  // do <n'|pert|n> integrals
  sum=0.0;
  left_value=(norm_n*norm_n_mark*genUnpertEigenfuncInteract(n_mark,m,left)*genUnpertEigenfuncInteract(n,m,left)*genPerturbation(left)*left)/2;
  right_value=(norm_n*norm_n_mark*genUnpertEigenfuncInteract(n_mark,m,right)*genUnpertEigenfuncInteract(n,m,right)*genPerturbation(right)*right)/2;
  for(int i=1;i<=points-1;i++) {
    x=i*step+left; 
    sum+=norm_n*norm_n_mark*genUnpertEigenfuncInteract(n_mark,m,x)*genUnpertEigenfuncInteract(n,m,x)*genPerturbation(x)*x;
  }
  sum=(sum+left_value+right_value)*step;
  return(sum);
}

double normalization(int n, int m,double left, double right, int points) {

  double step=(right-left)/((double)points);
  double sum=0.0;
  double x=0.0;

  double left_value=(genUnpertEigenfuncInteract(n,m,left)*genUnpertEigenfuncInteract(n,m,left)*left)/2;
  double right_value=(genUnpertEigenfuncInteract(n,m,right)*genUnpertEigenfuncInteract(n,m,right)*right)/2;

  for(int i=1;i<=points-1;i++) {
    x=i*step+left;
    sum+=genUnpertEigenfuncInteract(n,m,x)*genUnpertEigenfuncInteract(n,m,x)*x;
  }
  sum=(sum+left_value+right_value)*step;
  return(1.0/sqrt(sum));
}

// tqli method, taken from Numerical Recipes. Finds eigenvalues and eigenvectors of a tridiagonal matrix.

void tqli(vector<double> &d, vector<double> &e, vector<vector<double> > &z) {
  double pythag(double a, double b);
  int n=d.size();
  int m,l,iter,i,k;
  double s,r,p,g,f,dd,c,b;
  
  for (i=2;i<=n;i++) e[i-1]=e[i];
  e[n]=0.0;
  for (l=1;l<=n;l++) {
    iter=0;
    do {
      for (m=l;m<=n-1;m++) {
	dd=fabs(d[m])+fabs(d[m+1]);
	if ((double)(fabs(e[m])+dd) == dd) break;
      }
      if (m != l) {
	if (iter++ == 30) printf("Too many iterations in tqli");
	g=(d[l+1]-d[l])/(2.0*e[l]);
	r=pythag(g,1.0);
	g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
	s=c=1.0;
	p=0.0;
	for (i=m-1;i>=l;i--) {
	  f=s*e[i];
	  b=c*e[i];
	  e[i+1]=(r=pythag(f,g));
	  if (r == 0.0) {
	    d[i+1] -= p;
	    e[m]=0.0;
	    break;
	  }
	  s=f/r;
	  c=g/r;
	  g=d[i+1]-p;
	  r=(d[i]-g)*s+2.0*c*b;
	  d[i+1]=g+(p=s*r);
	  g=c*r-b;
	  /*
	  for (k=1;k<=n;k++) {
	    f=z[k][i+1];
	    z[k][i+1]=s*z[k][i]+c*f;
	    z[k][i]=c*z[k][i]-s*f;
	  }
	  */
	}
	if (r == 0.0 && i >= l) continue;
	d[l] -= p;
	e[l]=g;
	e[m]=0.0;
      }
    } while (m != l);
  }
}

// pythag method needed for tqli, taken from Numerical Recipes

double pythag(double a, double b) {
  double absa,absb;
  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

// tred2 method taken from Numerical Recipes. Reduces the matrix to tridiagonal form.

void tred2(vector<vector<double> > &a, vector<double> &d, vector<double> &e) {
  int n=d.size();
  int l,k,j,i;
  float scale,hh,h,g,f;
  
  for (i=n;i>=2;i--) {
    l=i-1;
    h=scale=0.0;
    if (l > 1) {
      for (k=1;k<=l;k++)
	scale += fabs(a[i][k]);
      if (scale == 0.0)
	e[i]=a[i][l];
      else {
	for (k=1;k<=l;k++) {
	  a[i][k] /= scale;
	  h += a[i][k]*a[i][k];
	}
	f=a[i][l];
	g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
	e[i]=scale*g;
	h -= f*g;
	a[i][l]=f-g;
	f=0.0;
	for (j=1;j<=l;j++) {
	  a[j][i]=a[i][j]/h;
	  g=0.0;
	  for (k=1;k<=j;k++)
	    g += a[j][k]*a[i][k];
	  for (k=j+1;k<=l;k++)
	    g += a[k][j]*a[i][k];
	  e[j]=g/h;
	  f += e[j]*a[i][j];
	}
	hh=f/(h+h);
	for (j=1;j<=l;j++) {
	  f=a[i][j];
	  e[j]=g=e[j]-hh*f;
	  for (k=1;k<=j;k++)
	    a[j][k] -= (f*e[k]+g*a[i][k]);
	}
      }
    } else
      e[i]=a[i][l];
    d[i]=h;
  }
  d[1]=0.0;
  e[1]=0.0;

  /* Contents of this loop can be omitted if eigenvectors not
     wanted except for statement d[i]=a[i][i]; */
  for (i=1;i<=n;i++) {
    /*
    l=i-1;
    if (d[i]) {
      for (j=1;j<=l;j++) {
	g=0.0;
	for (k=1;k<=l;k++)
	  g += a[i][k]*a[k][j];
	for (k=1;k<=l;k++)
	  a[k][j] -= g*a[k][i];
      }
    }
    */
    d[i]=a[i][i];
    /*
    a[i][i]=1.0;
    for (j=1;j<=l;j++) a[j][i]=a[i][j]=0.0;
    */
  }
}

// sorting eigenvalues, smallest first.
void sort_tridiag(vector<double> &d) {
  int n=d.size();
  for (int i = 1; i < n; i++) {
    int k = i;
      double temp = d[i];
      for (int j = i+1; j <=n; j++) {
	if (d[j] < temp) {
	  k = j;
	  temp = d[j];
	}
      }
      if (k != i) {
	d[k] = d[i];
	d[i] = temp;
      }
  }
}
