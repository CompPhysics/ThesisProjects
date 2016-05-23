#ifndef CoulombMatrix_CPP
#define CoulombMatrix_CPP

// file        : CoulombMatrix.cpp
// description : implentation of class CoulombMatrix, used to generate and to store the matrix of the Coulomb interaction term when the single particle wave functions are expanded in a given basis set.

#include <iostream>
#include <cmath>
#include <cstdlib>
#include "CoulombMatrix.h"
#include "singleParticleOrbitals.h"
#include <iomanip>

using namespace std;

// Contructors
CoulombMatrix:: CoulombMatrix(const singleParticleOrbitals* S, const double lbd):Basis(S)
{
#if DEBUGCONSTRUCTOR
  printf("Construction starts: CoulombMatrix\n");
#endif
  allocateMemory();
  setLambda(lbd);
  fillingCoulombMatrix();

#if DEBUGCONSTRUCTOR
  printf("Construction finished: CoulombMatrix \n");
#endif
}


// Destructor
CoulombMatrix:: ~CoulombMatrix ()
{
  delete [] Coulomb_Matrix;
}

// Allocate memory for the Coulomb matrix according to the size of the closed shell system
void CoulombMatrix:: allocateMemory ()
{
#if COMMENT
  printf("---Start allocating memory for the Coulomb Matrix\n");
#endif
  int NbChannels= Basis->readNbChannels();
  Coulomb_Matrix= new double *[NbChannels];
  
  int * NbCouplesPerCh=Basis->read_NbCouplesPerChannel();
  int NbCMelements= 0;
  for(int i=0; i<NbChannels; i++)
    {
      NbCMelements = (NbCouplesPerCh[i]+1)*NbCouplesPerCh[i]/2;
      Coulomb_Matrix[i]= new double [NbCMelements];
      for(int j=0; j<NbCouplesPerCh[i]; j++)
	Coulomb_Matrix[i][j]=0;
    }

#if COMMENT
  printf("---Memory allocated for the Coulomb Matrix\n");
#endif
}


// update the value of the interaction strength
void CoulombMatrix:: setLambda(const double lbd)
{
  if(lbd < 0)
    {
      printf("\n############ ERROR MESSAGE: The interaction strength Lambda must be positive ############\n");
      exit(1);
    }
  lambda = lbd;
}


// Return the number of element in the Coulomb Matrix
int CoulombMatrix:: sizeCoulombMatrix ()
{
  int sizeCM=0;
  int NbChannels= Basis->readNbChannels();
  int * NbCouplesPerCh=Basis->read_NbCouplesPerChannel();
  int NbCMelements= 0;
  for(int i=0; i<NbChannels; i++)
    {
      NbCMelements = (NbCouplesPerCh[i]+1)*NbCouplesPerCh[i]/2;
      sizeCM+=NbCMelements;
    }
  return sizeCM;
}


// Computes only the minimum number of Coulomb Matrix elements necessary to generate the full Coulomb matrix (use the fact that the full matrix is symmetric and slso the properties of the antisymmetric operator)
void CoulombMatrix:: fillingCoulombMatrix ()
{
  FILE * filename= fopen ("CMatrix.m","w");

#if COMMENT
  printf("---Start filling for the Coulomb Matrix\n");
#endif
  int NbChannels= Basis->readNbChannels();
  int * NbCouplesPerCh=Basis->read_NbCouplesPerChannel();
  int NbCMelements= 0;
  int NBCouples;
  int elem=0;
  int sizeCM=sizeCoulombMatrix();
  int countElement=0;
  int modulo=(int) sizeCM/5;

#if CHECK2
  printf("format:\n n1 m1 n2 m2 n3 m3 n4 m4 = <1 2|V|3 4>as = Coulomb Matrix element\n");
#endif
#if COMMENT
  fprintf(filename,"format:\n n1 m1 n2 m2 n3 m3 n4 m4 = <1 2|V|3 4>as = Coulomb Matrix element\n");
#endif

  for(int ch=0; ch<NbChannels; ch++)
    {
#if CHECK2
      printf("\n-------- channel #%d\n",ch);
#endif
#if COMMENT
      fprintf(filename,"\n-------- channel #%d\n",ch);
#endif
      NBCouples=NbCouplesPerCh[ch];// number of couples in the given channel
      for(int i=0;i<NBCouples;i++)
	for(int j=i;j<NBCouples;j++)
       {
	 elem=i*(2*NBCouples-i+1)/2+j-i;// ATTENTION, depends on how is build the matrix
#if COMMENT
#if modulo > 0
	 countElement++;
	 if(!(countElement%modulo))
	   printf("--- Element #%04d/%d\n",countElement,sizeCM,modulo);
#endif	 
#endif
	 Coulomb_Matrix[ch][elem]=compute_CoulombMatrixElement(ch,i,j, filename);// compute <i|V|j>as = <ab|V|cd>as
       }
    }
fclose(filename);

#if COMMENT
  printf("---Coulomb Matrix filled\n");
#endif
}// end of fillingCoulombMatrix()


// Write to file the content of the sparse Coulomb matrix
void CoulombMatrix:: save2file(ostream & mFile)
  {
    cout << "\nSaving Coulomb Matrix to file ..."; 
    int NbChannels= Basis->readNbChannels();
    int * NbCouplesPerCh=Basis->read_NbCouplesPerChannel();
    int NBCouples;
    int elem=0;
    mFile << endl << endl << "%%%%%%%%   Coulomb Matrix element classified by channel and elements in the channel (2 couples of states) in the channel %%%%%%%%" << endl;

    for(int ch=0; ch<NbChannels; ch++)
      {
	mFile << endl << "%-------- channel #" << ch << endl;
	NBCouples=NbCouplesPerCh[ch];// number of couples in the given channel
	mFile << "CM_ch" << setfill ('0') << setw (4) << ch << "=sparse("<< NBCouples <<","<< NBCouples <<");" << endl;
	
	for(int i=0;i<NBCouples;i++)
	  for(int j=i;j<NBCouples;j++)
	    {
	      elem=i*(2*NBCouples-i+1)/2+j-i;// ATTENTION, depends on how is build the matrix
	      if(Coulomb_Matrix[ch][elem] != 0)
		{
		  mFile << "CM_ch";
		  mFile << setfill ('0') << setw (4);
		  mFile << ch;
		  mFile << "(" << i+1 << "," << j+1 << ")=";
		  mFile.setf(ios::fixed,ios::floatfield);
		  mFile.precision(24);
		  mFile << Coulomb_Matrix[ch][elem]<<";"<<endl;
		}
	    }
      }
    cout << "  saved :)\n";
  }// end of save2file()


// Compute the coulomb matrix element <i|V|j>as
double CoulombMatrix:: compute_CoulombMatrixElement (const int channel, const int i, const int j, FILE * filename)
{

 
  double CMelement=0;
  double directT=0;
  double exchangeT=0;
  int dim=Basis->Dimension();  // Calculation might depend on the dimension of the system
  int index_ml=Basis->readIndex_ml();
  int index_ms=Basis->readIndex_ms();

  // read the 2 single particle states of the 2 couple of state
  int stateBra1, stateBra2, stateKet1, stateKet2;// compute <i|V|j>as = <bra1 bra2|V| ket1 ket2>as
  stateBra1=Basis->mappingCoupleToSingleStates(channel,i,0);
  stateBra2=Basis->mappingCoupleToSingleStates(channel,i,1);
  stateKet1=Basis->mappingCoupleToSingleStates(channel,j,0);
  stateKet2=Basis->mappingCoupleToSingleStates(channel,j,1);
  
  // read the quantum numbers of each single particle state
  int n1,l1,ml1,ms1,n2,l2,ml2,ms2,n3,l3,ml3,ms3,n4,l4,ml4,ms4;// quantum numbers of the 4 single particle states
  n1=Basis->readQuantumNumberOfState(stateBra1,0);
  n2=Basis->readQuantumNumberOfState(stateBra2,0);
  n3=Basis->readQuantumNumberOfState(stateKet1,0);
  n4=Basis->readQuantumNumberOfState(stateKet2,0);
  if(dim == 3)
    {
      l1=Basis->readQuantumNumberOfState(stateBra1,1);
      l2=Basis->readQuantumNumberOfState(stateBra2,1);
      l3=Basis->readQuantumNumberOfState(stateKet1,1);
      l4=Basis->readQuantumNumberOfState(stateKet2,1);
    }
  ml1=Basis->readQuantumNumberOfState(stateBra1,index_ml);
  ml2=Basis->readQuantumNumberOfState(stateBra2,index_ml);
  ml3=Basis->readQuantumNumberOfState(stateKet1,index_ml);
  ml4=Basis->readQuantumNumberOfState(stateKet2,index_ml);

  ms1=Basis->readQuantumNumberOfState(stateBra1,index_ms);
  ms2=Basis->readQuantumNumberOfState(stateBra2,index_ms);
  ms3=Basis->readQuantumNumberOfState(stateKet1,index_ms);
  ms4=Basis->readQuantumNumberOfState(stateKet2,index_ms);
  
  
   if(dim==3)
    {
      printf("COULOMB MATRIX NOT AVAILABLE YET IN 3D\n");
      //exit(1);
    }
   else // dim ==2
    {
      // compute <12|V|34>as = <12|V|34> - <12|V|43> = directTerm - exchangeTerm
      if(ms1==ms4 && ms2==ms3)
	exchangeT= anisimovas(n1,ml1,n2,ml2,n3,ml3,n4,ml4);
      
      if(ms1==ms3 && ms2==ms4)
	directT= anisimovas(n1,ml1,n2,ml2,n4,ml4,n3,ml3);
      // Note that Rontani and Anisimovas defined the exchange term with their function <a(1)b(2)|V|c(2)d(1)>, whereas people usually keep the order of the particles <a(1)b(2)|V|c(1)d(2)>
      
      CMelement = lambda*(directT - exchangeT);

#if CHECK2
      //printf("  --- (ms1,ms4)=(%2d,%2d)-----(ms2,ms3)=(%2d,%2d)------  <%2d,%2d|V|%2d,%2d> = %lf\n",ms1,ms4,ms2,ms3,stateBra1,stateBra2,stateKet1,stateKet2,CMelement);

      // Output is written one line at a time, each line has one matrix element. The line has the format: n1 m1 n2 m2 n3 m3 n4 m4 element
      if(true)// CMelement != 0)
	{
	  printf("%2d %2d %2d %2d %2d %2d %2d %2d =",n1, ml1, n2, ml2, n3, ml3, n4, ml4);
	  printf(" <%2d,%2d|V|%2d,%2d>as = %+4.16f\n",stateBra1,stateBra2,stateKet1,stateKet2,CMelement);
	}
#endif
#if COMMENT
	  if(true)// CMelement != 0)
	    {
	      fprintf(filename,"%2d %2d %2d %2d %2d %2d %2d %2d =",n1, ml1, n2, ml2, n3, ml3, n4, ml4);
	      fprintf(filename," <%2d,%2d|V|%2d,%2d>as = %+4.16f\n",stateBra1,stateBra2,stateKet1,stateKet2,CMelement);
	    }
#endif
    }
  return CMelement;
}// end ofcompute_CoulombMatrixElement()


// compute one non-antisymmetric part of the Coulomb matrix element: <12|V|34> [see Asinimovas and Matulis (1997) for analytical expression]
// does NOT compute <12|V|34>as
double anisimovas (const int n1, const int m1, const int n2, const int m2, const int n3, const int m3, const int n4, const int m4)
{
  double matulis=0.0;
  if(m1+m2 == m3+m4)
    {
      double temp;
      int lambda;

      int gamma1=0;
      int gamma2=0;
      int gamma3=0;
      int gamma4=0;
      int G=0;
      for(int j1=0;j1<=n1;j1++)
	for(int j2=0;j2<=n2;j2++)
	  for(int j3=0;j3<=n3;j3++)
	    for(int j4=0;j4<=n4;j4++)
	      {
		gamma1 = (int) (j1 +j4 + 0.5*(abs(m1)+m1) + 0.5*(abs(m4)-m4));
		gamma2 = (int) (j2 +j3 + 0.5*(abs(m2)+m2) + 0.5*(abs(m3)-m3));
		gamma3 = (int) (j2 +j3 + 0.5*(abs(m3)+m3) + 0.5*(abs(m2)-m2));
		gamma4 = (int) (j1 +j4 + 0.5*(abs(m4)+m4) + 0.5*(abs(m1)-m1));

		G=gamma1+gamma2+gamma3+gamma4;

		double Lgratio1 = LogRatio1(j1,j2,j3,j4);
		double Lgproduct2 = LogProduct2(n1,m1,n2,m2,n3,m3,n4,m4,j1,j2,j3,j4);
		double Lgratio2 = LogRatio2(G);

		temp=0.0;
		lambda=0;
		for(int l1=0;l1<=gamma1;l1++)
		  for(int l2=0;l2<=gamma2;l2++)
		    for(int l3=0;l3<=gamma3;l3++)
		      for(int l4=0;l4<=gamma4;l4++)
			{
			  lambda=l1+l2+l3+l4;
			  if( (l1+l2)==(l3+l4) )
			    temp += minusPower(gamma2+gamma3-l2-l3)*exp(LogProduct3(l1,l2,l3,l4,gamma1,gamma2,gamma3,gamma4)+lgamma(1+lambda*0.5)+lgamma((G-lambda+1)*0.5));
			}


		matulis += minusPower(j1+j2+j3+j4) *exp(Lgratio1 + Lgproduct2 + Lgratio2) *temp;

	      }
      matulis *= Product1(n1,m1,n2,m2,n3,m3,n4,m4);

#if CHECK2
      // print both part of the antisymmetrized CM element if non-zero
      if(matulis !=0)
	printf(" <%2d %2d,%2d %2d|V|%2d %2d,%2d %2d> = %4.16f (non anti-sym)\n",n1,m1,n2,m2,n3,m3,n4,m4,matulis);
#endif
    }
  return matulis;
}




// Computes the first ratio in the Asinimovas expression
double LogRatio1(const int j1,const int j2,const int j3,const int j4)
{
  double temp = -LogFac(j1)-LogFac(j2)-LogFac(j3)-LogFac(j4);
#if CHECK2
  printf("LogRatio1=%lf\n",temp);
#endif
  return temp;
} 

// Computes the 2nd ratio in the Asinimovas expression
double LogRatio2(const int G)
{
  double temp = -1*(G+1)*0.5*log(2);
#if CHECK2
  printf("LogRatio2=%lf\n",temp);
#endif
  return temp;
} 


// Computes the log of the 2nd product in the Asinimovas expression
double LogProduct2(const int n1,const int m1,const int n2,const int m2,const int n3,const int m3,const int n4,const int m4,const int j1,const int j2,const int j3,const int j4)
{
  double temp = LogFac(n1+abs(m1))+LogFac(n2+abs(m2))+LogFac(n3+abs(m3))+LogFac(n4+abs(m4))-LogFac(n1-j1)-LogFac(n2-j2)-LogFac(n3-j3)-LogFac(n4-j4)-LogFac(j1+abs(m1))-LogFac(j2+abs(m2))-LogFac(j3+abs(m3))-LogFac(j4+abs(m4));
#if CHECK2
  printf("LogProd2=%lf\n",temp);
#endif
  return temp;
}


// Computes the log of the 3rd product in the Asinimovas expression
double LogProduct3(const int l1,const int l2,const int l3,const int l4,const int gamma1,const int gamma2,const int gamma3,const int gamma4)
{
  double temp = LogFac(gamma1)+LogFac(gamma2)+LogFac(gamma3)+LogFac(gamma4)-LogFac(l1)-LogFac(l2)-LogFac(l3)-LogFac(l4)-LogFac(gamma1-l1)-LogFac(gamma2-l2)-LogFac(gamma3-l3)-LogFac(gamma4-l4);
#if CHECK2
  printf("LogProd3=%d\n",temp);
#endif
  return temp;
}


// computes (-1)^k
int minusPower(const int k)
{
  int temp=abs(k%2)*-2+1;// gives 1 if k is even, gives -1 if k is odd
#if CHECK2
  printf("MinusPower=%d\n",temp);
#endif
  return temp;
}

// given 4 states (i.e.their index), this function reads the corresponding Coulomb Matrix elements
double CoulombMatrix:: mappingSingleStatesToCMelement (const int s1,const int s2,const int s3,const int s4) const
{
  double CMelement=0.0;
  int chBra=Basis->read_channelOfCouple(s1,s2);
  int chKet=Basis->read_channelOfCouple(s3,s4);
  if(chBra != -1 || chKet != -1) // at least one couple detected as made of 2 identical single particle states, CMelement should be zero
    {
      int signChange=1;
      if(chBra == chKet) //conservation of the total angular momentum and (spin????)
	{
      
	  int coupleBra=Basis->mappingSingleStatesToCouple(s1,s2);
	  int coupleKet=Basis->mappingSingleStatesToCouple(s3,s4);
	  // use properties of the antisymmetric operator to find the Coulomb
	  // matrix element <bra|V|ket>as=<ab|V|cd>as so that a<b, c<d and bra<ket
	  if(s1>s2)
	    signChange*=-1;
	  if(s3>s4)
	    signChange*=-1;
	  if(coupleBra>coupleKet)
	    CMelement=read_CoulombMatrixElement(chBra,coupleKet,coupleBra);
	  else
	    CMelement=read_CoulombMatrixElement(chBra,coupleBra,coupleKet);
	}
      CMelement*=signChange;
    }
  return CMelement;
}

// given the indices for 2 couple of states |i> and |j>, assuming i<j, read the corresponding Coulomb Matrix element
double CoulombMatrix:: read_CoulombMatrixElement(const int channel, const int bra, const int ket) const
{
  double CMelement=0.0;
  int * NCp=Basis->read_NbCouplesPerChannel();
  int indexElement = bra * (2*NCp[channel]-bra+1)/2+ket-bra; ////////ATTENTION , depend on how is build the Coulomb Element matrix
  CMelement=Coulomb_Matrix[channel][indexElement];
  return CMelement;
}

//  // this function reads the index of the corresponding couple of states for a given channel=2M+S and the index of an element in this channel, 
// int CoulombMatrix:: mappingCMelementToCoupleStates (const int channel, const int indexElement, const int CoupleStates)
//  {
//    int * NbCouplesPerCh=Basis->read_NbCouplesPerChannel();
//    int NBCouples=NbCouplesPerCh[channel];// number of couples in the given channel
//    int indexCouple[2];
//    int elem;// index of the element in the array representing the coulomb (sub)matrix
// #if CHECK
//    int NbCMelements = (NbCouplesPerCh[channel]-1)*NbCouplesPerCh[channel]/2;
//    if(indexElement >= NbCMelements)
//      printf("#### in mappingCMelementToCoupleStates(), you asked for a too big index of element  ####\n");
// #endif
//    for(int a=0;a<NBCouples-1;a++)
//      for(int b=a+1;b<NBCouples;b++)
//        {
// 	 elem=a*(2*NBCouples-a-1)/2+b-(a+1);
// 	 if(elem == indexElement)
// 	   {
// 	     indexCouple[0]=a;//bra <i| (=<ab|)
// 	     indexCouple[1]=b;//ket |j> (=|cd>)
// 	     //printf("--- <%d|V|%d>\n",indexCouple[0],indexCouple[1]);
// 	   }
//        }
//    return indexCouple[CoupleStates];
//  }

// computes log(n!)
double LogFac (const int n)
{
#if CHECK
  if(n>170)
    {
      printf("#### Too big integer in LogFac(n)!!!!#### \n\n");
      exit(1);
    }
#endif
  double temp=0;
  for(int i=2; i<n+1;i++)
    temp+=log(i);
#if CHECK2
  printf("LogFac=%lf\n",temp);
#endif
  return temp;
}

// computes first product of indices in the Anisimovas/Matulis expression
double Product1 (const int n1, const int ml1, const int n2, const int ml2, const int n3, const int ml3, const int n4, const int ml4)
{
  double temp=0;
  temp = LogFac(n1)+LogFac(n2)+LogFac(n3)+LogFac(n4)-LogFac(n1+abs(ml1))-LogFac(n2+abs(ml2))-LogFac(n3+abs(ml3))-LogFac(n4+abs(ml4));
  temp*=0.5;
#if CHECK2
  printf("Prod1=%lf\n",exp(temp));
#endif
  return exp(temp);
}

// // computes the factorial of a positive integer
// int factorial(int number) {
// 	int temp=0;

// 	if(number <= 1) return 1;

// 	temp = number * factorial(number - 1);
// 	return temp;
// }
// // End: function factorial()


// computes the factorial of a positive integer using log functions
double FastFactorial(int number) {
	double temp=0;
	for(int i=1;i<number+1;i++)
	  temp+=log(i);
	return  exp(temp);
}
// End: function factorial()


// // The function gamma() computes the gamma function of a real x
// // // taken on the web [http://www.crbond.com/math.htm]
// // double gamma(double x)
// // {
// //     int i,k,m;
// //     double ga,gr,r,z;

// //     static double g[] = {
// //         1.0,
// //         0.5772156649015329,
// //        -0.6558780715202538,
// //        -0.420026350340952e-1,
// //         0.1665386113822915,
// //        -0.421977345555443e-1,
// //        -0.9621971527877e-2,
// //         0.7218943246663e-2,
// //        -0.11651675918591e-2,
// //        -0.2152416741149e-3,
// //         0.1280502823882e-3,
// //        -0.201348547807e-4,
// //        -0.12504934821e-5,
// //         0.1133027232e-5,
// //        -0.2056338417e-6,
// //         0.6116095e-8,
// //         0.50020075e-8,
// //        -0.11812746e-8,
// //         0.1043427e-9,
// //         0.77823e-11,
// //        -0.36968e-11,
// //         0.51e-12,
// //        -0.206e-13,
// //        -0.54e-14,
// //         0.14e-14};

// //     if (x > 171.0) return 1e308;    // This value is an overflow flag.
// //     if (x == (int)x) {
// //         if (x > 0.0) {
// //             ga = 1.0;               // use factorial
// //             for (i=2;i<x;i++) {
// // 	      ga *= i; // for positiv integer
// //             }
// //          }
// //          else
// //             ga = 1e308;
// //      }
// //      else {
// //         if (fabs(x) > 1.0) {
// //             z = fabs(x);
// //             m = (int)z;
// //             r = 1.0;
// //             for (k=1;k<=m;k++) {
// //                 r *= (z-k);
// //             }
// //             z -= m;
// //         }
// //         else
// //             z = x;
// //         gr = g[24];
// //         for (k=23;k>=0;k--) {
// //             gr = gr*z+g[k];
// //         }
// //         ga = 1.0/(gr*z);
// //         if (fabs(x) > 1.0) {
// //             ga *= r;
// //             if (x < 0.0) {
// //                 ga = -M_PI/(x*ga*sin(M_PI*x));
// //             }
// //         }
// //     }
// //     return ga;
// // }
// // // end of the Gamma function



// lgamma.cpp -- log gamma function of real argument.
//      Algorithms and coefficient values from "Computation of Special
//      Functions", Zhang and Jin, John Wiley and Sons, 1996.
//
//  (C) 2003, C. Bond. All rights reserved.
//
//  Returns log(gamma) of real argument.
//  NOTE: Returns 1e308 if argument is 0 or negative.
//  taken on the web [http://www.crbond.com/math.htm]
double lgamma(double x)
{
#if CHECK
  if(x>=171)
    {
    printf("#### Too big integer in lgamma(x) to give accurate result!!!!#### \n\n");
    exit(1);
    }
#endif  
  
  double x0,x2,xp,gl,gl0;
  int n,k;
  static double a[] = {
    8.333333333333333e-02,
       -2.777777777777778e-03,
        7.936507936507937e-04,
       -5.952380952380952e-04,
        8.417508417508418e-04,
       -1.917526917526918e-03,
        6.410256410256410e-03,
       -2.955065359477124e-02,
        1.796443723688307e-01,
       -1.39243221690590};
    
    x0 = x;
    if (x <= 0.0) return 1e308;
    else if ((x == 1.0) || (x == 2.0)) return 0.0;
    else if (x <= 7.0) {
        n = (int)(7-x);
        x0 = x+n;
    }
    x2 = 1.0/(x0*x0);
    xp = 2.0*M_PI;
    gl0 = a[9];
    for (k=8;k>=0;k--) {
        gl0 = gl0*x2 + a[k];
    }
    gl = gl0/x0+0.5*log(xp)+(x0-0.5)*log(x0)-x0;
    if (x <= 7.0) {
        for (k=1;k<=n;k++) {
            gl -= log(x0-1.0);
            x0 -= 1.0;
        }
    }
#if CHECK2
    printf("lgamma=%lf\n",gl);
#endif
    return gl;
}// end of lgamma function



#endif
