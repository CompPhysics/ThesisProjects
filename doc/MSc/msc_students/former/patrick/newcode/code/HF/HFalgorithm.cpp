#ifndef HFalgorithm_CPP
#define HFalgorithm_CPP

// file        : HFalgorithm.cpp
// description : implentation of class HFalgorithm

#include <cmath>
#include <cstdlib>
#include "HFalgorithm.h"
#include "lib.h"

using namespace std;

// Contructors
HFalgorithm :: HFalgorithm(const singleParticleOrbitals* S, const singleParticleWaveFunctions* B, const singleParticleEnergies * E, const CoulombMatrix * CM, const double epsi):States(S), Basis(B), BasisEnergies(E),CMatrix(CM)
{
#if DEBUGCONSTRUCTOR
  printf("\nConstruction starts: HFalgorithm\n");
#endif
 nbS = States->readTotalNbStates();
 allocateMemory();
 setEpsilon(epsi);

#if DEBUGCONSTRUCTOR
  printf("Construction finished: HFalgorithm\n");
#endif
}

// Allocate memory to the matrix of wave functions
void HFalgorithm:: allocateMemory ()
{
#if COMMENT
  printf("---Start allocating memory for the energies\n");
#endif
  int NbStates= nbS;
  Range r(0,NbStates-1);
  HForbitalsCoeff.resize(r,r);//  matrix of the coefficient for the HF orbitals
  HForbitalsCoeff = 0.0;
    
  EffectivePotential.resize(r,r);// effective potential matrix
  EffectivePotential = 0.0;

  ParticleEnergies.resize(r); // single particle energies (eigenvalues of the problem)
  ParticleEnergies = 0.0;

  // Memory allocation for the matrix of the eigenvalue problem
  A = new double*[NbStates];
  for(int s=0;s <NbStates;s++)
    A[s]= new double[NbStates];
  
  
#if COMMENT
  printf("---Memory allocated for the energies\n");
#endif
}


// Destructor
HFalgorithm:: ~HFalgorithm ()
{
  deallocate();
}

// Free dynamic memory
void HFalgorithm:: deallocate ()
{
  delete [] A;
}

// update the value of the precision used in the Hartree-Fock algorithm (Energy precision)
void HFalgorithm:: setEpsilon(const double epsil)
{
  epsilon = abs(epsil);
}

// Initialisation of the single particle energies/wave functions
void HFalgorithm::initialization ()
{
  // initialsize the HF orbitals as one particle in each pure state of the basis set
  // initialize the single particle energies as HO single particle energies
  for(int s=0;s <nbS;s++)
    {
      HForbitalsCoeff(s,s) = 1.0;
      ParticleEnergies(s) = BasisEnergies->read_Energy(s);
      for(int s2=0;s2 <nbS;s2++)
	A[s][s2]=0.0;
    }
  //HForbitalsCoeff(1,0)=1;// test


  // normalize the resulting HF orbitals
  normalize_HForbitals();
#if COMMENT
  cout << "Initial Coeff.(HOstate, particle): " << HForbitalsCoeff  << endl;
#endif

  // compute the effective potential according to the initial HF orbitals
  EffectivePotential = 0.0;
  compute_EffectivePotential();
}  


// update the effective potential with the optimized single particle wfs
void HFalgorithm::compute_EffectivePotential ()
{
  // zeroing the matrix of the effective potential
   for(int alpha=0;alpha <nbS;alpha++)
     for(int gamma=0; gamma<nbS;gamma++)
       EffectivePotential(alpha, gamma)= 0.0;

  double C1,C2,CMelement;
  // Compute U(alpha, beta)=sum over particle, of the sum over states beta and delta of: V(alpha,beta,gamma,delta) * C(particle,beta) * C(particle,delta)
  for(int alpha=0;alpha <nbS;alpha++)
    for(int gamma=0; gamma<nbS;gamma++)
      {
	for(int particle=0; particle<nbS;particle++)
	  for(int beta=0; beta<nbS;beta++)
	    for(int delta=0; delta<nbS;delta++)
	      {
		C1 = HForbitalsCoeff(beta,particle);
		C2 = HForbitalsCoeff(delta,particle);
		CMelement = CMatrix->mappingSingleStatesToCMelement(alpha,beta,gamma,delta);
		EffectivePotential(alpha, gamma) += C1*C2*CMelement;
#if CHECK2
		printf("%d %d %d %d, particle:%d C1=%12.4E, C2=%12.4E, CMelemt=%12.4E, Eff.Pot[%d,%d]=%12.4E\n",alpha, beta, gamma, delta,particle,C1,C2,CMelement,alpha, gamma, EffectivePotential(alpha, gamma));
		cin.get();
#endif
	
	      }
      }// end of loops over element of U (alpha and gamma)
 

#if CHECK
  cout << "Effective Potential U(alpha, gamma): " <<EffectivePotential << endl;
#endif
}   


// Solve the linear system given by HF equations
void HFalgorithm::solve_HFequations ()
{
  // Solve the following eigenvalue problem: A.X=eX  where
  //     - A= HO single particle energies + effective potential
  //     - X: (eigenvector) corresponds to the coefficients of the expansion of one particle in the HO basis
  //     - e: (eigenvalues) corresponds to the eigenEnergy of one particle in the system
  

  // Zeroing A
  for(int alpha=0; alpha<nbS; alpha++)
    for(int gamma=0; gamma<nbS; gamma++)
      A[alpha][gamma] = 0.0;  

  // Build A
  for(int alpha=0; alpha<nbS; alpha++)
    {
      A[alpha][alpha] = ParticleEnergies(alpha);
       for(int gamma=0; gamma<nbS; gamma++)
  	A[alpha][gamma] += EffectivePotential(alpha,gamma);
    }

#if CHECK
  printf("\nMatrix to block-diagonalize");
  printScreen_A();
#endif
  
  double *diag, *eigenE;
  diag= new double [nbS]; 
  eigenE = new double [nbS];

  // Zeroing diag and eigenE
  for(int s=0;s<nbS;s++)
    {
      diag[s]=0.0;
      eigenE[s]=0.0;
    }
  
  // solve an eigenvalue problem using Householder's method for tridiagonalization and solve the eigenValue problem of the tridiagonal matrix
  tred2(A, nbS, diag, eigenE);
  tqli(diag,eigenE,nbS, A);// A now contains the normalized eigenvectors and diag the eigenvalues
  

  // Eigenvectors become the new HF orbital coefficients, don't do anything of eigenValues (eigenEnergies)
  for(int particle=0; particle<nbS;particle++)
      for(int s=0; s<nbS;s++)
	HForbitalsCoeff(s,particle)=A[s][particle];
 
#if COMMENT
  cout.precision(15);  
  cout.setf(ios::scientific);
  cout << "Coeff= (after solving HF eq.)";
  saveCoeff2file(cout);
#endif


#if 1
  printf("\neigenvalues: E={");
  for(int s=0;s<nbS-1;s++)
    {
      cout.precision(6);  
      cout.setf(ios::scientific);
      cout << diag[s] << ",";
    }
  cout << diag[nbS-1] << "}\n";
#endif

 }// end of HFalgorithm::solve_HFequations ()


// Update the total energy of the system
double HFalgorithm::compute_TotalEnergy()
{
  double HOenergy, CoulombRepulsionEnergy;
  HOenergy = 0.0;
  CoulombRepulsionEnergy = 0.0;
  double CMelement;
  
  for(int i=0; i<nbS;i++) // loop over particles, as many as # of states in closed shell systems
    {
      for(int alpha=0;alpha <nbS;alpha++)
	HOenergy += pow(HForbitalsCoeff(alpha,i),2)*ParticleEnergies(alpha);

      for(int j=0;j <nbS;j++)// second loop over particles for the contribution of interactions
	for(int alpha=0;alpha <nbS;alpha++)
	  for(int beta=0;beta <nbS;beta++)
	    for(int gamma=0; gamma<nbS;gamma++)
	      for(int delta=0; delta<nbS;delta++)
		{
		  CMelement = CMatrix->mappingSingleStatesToCMelement(alpha, beta, gamma, delta);
		  CoulombRepulsionEnergy += HForbitalsCoeff(alpha,i)*HForbitalsCoeff(beta,j)*CMelement*HForbitalsCoeff(gamma,i)*HForbitalsCoeff(delta,j);
		}
    }// end of loop over particles
  CoulombRepulsionEnergy *= 0.5;
  TotalEnergy=HOenergy + CoulombRepulsionEnergy;
#if CHECK2
  printf("\nHOEnergy Eho=%12.8E\n",HOenergy);
  printf("Coulomb Repulsion energy Ec=%12.8E\n",CoulombRepulsionEnergy);
#endif 

#if true
  cout.precision(15);  
  cout.setf(ios::scientific);
  //   cout.setf(ios::fixed,ios::floatfield);
  //   cout.precision(25);
  cout << "Total energy E=" << TotalEnergy;
#endif 
  double TotalEnergy;
} // end of compute_TotalEnergy()



void HFalgorithm::check_selfConsistency(){}              // Check for self consistency

// normalization of the vectors of coefficient describing HF orbitals
void HFalgorithm::normalize_HForbitals()
{
  double sum;
  int nbS=States->readTotalNbStates();
  for(int particle=0; particle<nbS; particle++)
    {
      sum =0.0;
      for(int c=0; c<nbS; c++)
	sum+=pow(HForbitalsCoeff(c,particle),2);
      
      sum=sqrt(sum);
      for(int c=0; c<nbS; c++)
	HForbitalsCoeff(c,particle)/=sum;
    }
}


// Print the matrix of the eigenvalue problem to screen
void HFalgorithm:: printScreen_A() const
{
  int count=0;
  printf("\nA=\n");
  for(int alpha=0; alpha<nbS; alpha++)
    for(int gamma=0; gamma<nbS; gamma++)
      {

	if(count==nbS)
	  {
	    printf("\n%12.4E ",A[alpha][gamma]);
	    count =0;
	  }
	else
	  printf("%12.4E ",A[alpha][gamma]);

	count++;
      }
  cout << endl << endl;
}

// Write the eigenvalue problem to file in matlab format
void HFalgorithm:: save2file(ostream & mFile)
{
    saveCoeff2file(mFile);
    saveEffPot2file(mFile);
  }// end of save2file()

// Write only the HF orbitals coeff. to file in matlab format
void HFalgorithm:: saveCoeff2file(ostream & mFile)
{
  mFile.precision(7);
  mFile.setf(ios::fixed,ios::floatfield);
  cout << "Saving the coeff. of the HF orbitals to file ..."; 
  mFile << endl << endl << "%%%%%%%%   Coefficients of the HF orbitals (one column per particle, one row per single particle state) %%%%%%%%" << endl;
  int sizeMatrix = HForbitalsCoeff.size();
  mFile << "HForbitals = [";
  for(int i=1;i<=sqrt(sizeMatrix);i++)
    {
      for(int j=1;j<=sqrt(sizeMatrix);j++)
		  mFile << HForbitalsCoeff(i-1,j-1)<< " ";
	
      mFile << endl;
    }
  mFile << "];" << endl;
  cout << "  saved :)\n";
}// end of saveCoeff2file()

// Write only the Effective Potential matrix to file in matlab format
void HFalgorithm:: saveEffPot2file(ostream & mFile)
{
  cout << "Saving the effective potential to file ..."; 
  mFile << endl << endl << "%%%%%%%%   Effective potential (row and column corresponding to single particle states in the order given by singleParticleOrbitals, so sorted by identical angular momentum) %%%%%%%%" << endl;
  int sizeMatrix = EffectivePotential.size();
  mFile << "U = [";
  for(int i=1;i<=sqrt(sizeMatrix);i++)
    {
      for(int j=1;j<=sqrt(sizeMatrix);j++)
	mFile << EffectivePotential(i-1,j-1)<< " ";
      mFile << endl;
    }
  mFile << "];" << endl;
  cout << "  saved :)\n";
}// end of saveEffPot2file()

// Write only the Total Energy to file in matlab format
void HFalgorithm:: saveEnergy2file(ostream & mFile)
{
  mFile.precision(20);
  mFile.setf(ios::fixed,ios::floatfield);
  mFile << "Total_Energy= " << TotalEnergy;
}

#endif // HFalgorithm_CPP

