//#include <iostream>
#include "singleParticleOrbitals.h"
#include "CoulombMatrix.h"


//Morten includes in lib.h
// #include <new>
// #include <cstdio>
// #include <cstdlib>
// #include <cstring>


using namespace std;
void create_object(singleParticleOrbitals* Base){
  Base = new singleParticleOrbitals(2, 2);
}

int main()
{


    // parameters later to be read from files
  int Nlevel;        // Maximum level of energy taken into account
  int dim;           // spatial dimension of the system
  double lambda;     // dimensionless coeff for the interaction strength


  // Prompt the user to enter the parameter one by one
  cout << "Enter the dimension of the system {2 or 3}: \n";
  cin >> dim;
  cout << "Enter the number of the maximum energy level to consider: \n";
  cin >> Nlevel;
  cout << "Enter the dimensionless coeff for the interaction strength (lambda): \n";
  cin >> lambda;


  // CREATE A SECOND TABLE LISTING STATES WITH SAME {M,S}
  singleParticleOrbitals Base(dim,Nlevel);
    
  
  //create_object(Base);
  //if (true)
  //Base = new singleParticleOrbitals(dim, Nlevel);
  
    
  // CREATE A SECOND TABLE LISTING STATES WITH SAME {M,S}
  CoulombMatrix CM(&Base, lambda);

  
#if 0
  // TEST FINDING ELEMENTS IN ANY MATRIX (not related to the CoulombMatrix class)
  int nb_couples=5;
  double ** matrix = new double * [nb_couples];
  for(int i=0;i<nb_couples;i++)//allocate memory to this matrix and fill it with random numbers
    {
      matrix[i]=new double [nb_couples];
      for(int j=0;j<nb_couples;j++)
	matrix[i][j]=rand();
    }
  //print all the content of the matrix
  for(int i=0;i<nb_couples;i++)
    for(int j=0;j<nb_couples;j++)
      printf("matrix[%d][%d]=%lf\n",i,j,matrix[i][j]);
  int nb_elem=(nb_couples-1)*(nb_couples-1+1)/2;
  printf("nb_elem=%d\n",nb_elem);
  double * matrix2 = new double [nb_elem];

  int elem;
  for(int i=0;i<nb_couples-1;i++)
    for(int j=i+1;j<nb_couples;j++)
      {
	elem=i*(2*nb_couples-i-1)/2+j-(i+1);
	matrix2[elem]=matrix[i][j];
	printf("matrix2[%d]=matrix[%d][%d]=%lf\n",elem,i,j,matrix2[elem]);
      }
#endif

#if 0
  // TEST READING THE STATE COUPLE INDICES OF A COULOMB MATRIX ELEMENT
  int indexElement;
  int channel,Bra,Ket;
  while(true)
    {
      // Prompt the user to enter the parameter one by one
      cout << "Enter a channel \n";
      cin >> channel;
      cout << "Enter the index element of the Coulomb matrix\n";
      cin >> indexElement;
      
      Bra= CM.mappingCMelementToCoupleStates (channel,indexElement,0);
      Ket= CM.mappingCMelementToCoupleStates (channel,indexElement,1);
      printf("Coulomb Element [%d]-> <%d|V|%d>=<%d %d|V|%d %d>\n",indexElement,Bra,Ket,Base.mappingCoupleToSingleStates(channel,Bra,0),Base.mappingCoupleToSingleStates(channel,Bra,1),Base.mappingCoupleToSingleStates(channel,Ket,0),Base.mappingCoupleToSingleStates(channel,Ket,1));
	  
    }	
#endif



#if 0
  // TEST CHECKING FACTORIAL FUNCTIONS
  // Prompt the user to enter the parameter one by one
  int entier;
  double entierD;
   
  while(true)
    {
      cout << "Enter a number \n";
      cin >> entier;
      entierD=(double) entier +1;
//       printf("entierD=%lf\n",entierD);
//       printf("factorial(%d)=%d\n",entier,factorial(entier));
//       printf("fastFactorial(%d)=%lf\n",entier,FastFactorial(entier));
//       printf("gamma(%d+1)=(%d)!=%lf\n",entier,entier,gamma(entierD));
//       printf("(gamma-fasFact)=%lf\n\n",(gamma(entierD)-FastFactorial(entier)));
//       printf("log(%d)=%lf\n",entier, log(entier));
//       printf("exp(%d)=%lf\n",entier, exp(entier));
//       printf("logFac(%d)=%lf\n",entier,LogFac(entier));
//      printf("log(%d!)=%lf\n",entier,log(gamma(entier+1)));
      printf("log(gamma(%d))=%lf\n",entier,log(gamma(entier)));
      printf("lgamma(%d)=%lf\n",entier,lgamma(entier));
      
  
    }

#endif


#if 0
  // TEST OF COMVERTING DOUBLE TO INT, not working the other way round
  int m;
  int j=2;

  while(true)
    {
      cout << "Enter a number for m\n";
      cin >> m;
      double gammaD= (abs(m)+m)*0.5;
      int gammaI=(int) gammaD;
      printf("m=%d\n",m);
      printf("gammaD=(abs(m)+m)/2=%lf\n",gammaD);
      printf("gammaI=%d\n",gammaI);
      printf("/*(-1)^(%d % 2) = %d\n",m,abs(m%2)*-2+1);
      printf("/*(-1)^(%d % 2) = %d\n",m,minusPower(m));
    }
#endif
  

#if 1
  // TEST FOR READING THE VALUE OF A COULOMB MATRIX ELEMENT WITH 4 STATE INDICES
 int s1,s2,s3,s4;

  while(true)
    {
      cout << "Enter 4 indices of single particle states\n";
      cin >> s1;
      cin >> s2;
      cin >> s3;      
      cin >> s4;
      double element= CM.mappingSingleStatesToCMelement(s1,s2,s3,s4);
      printf("<%02d,%02d|V|%02d,%02d>as = %lf\n",s1,s2,s3,s4,element);
    }
#endif

  printf("test1.cpp finished\n");

  return 0;
}
