#ifndef singleParticleEnergies_CPP
#define singleParticleEnergies_CPP

// file        : singleParticleEnergies.cpp
// description : implentation of class singleParticleEnergies

#include <iostream>
#include <cmath>
#include <cstdlib>
#include "singleParticleEnergies.h"

using namespace std;

// Contructors
singleParticleEnergies :: singleParticleEnergies(const singleParticleOrbitals* S, const singleParticleWaveFunctions* B):States(S), Basis(B)
{
#if DEBUGCONSTRUCTOR
  printf("\nConstruction starts: singleParticleEnergies\n");
#endif
 
 allocateMemory();
 computeEnergies();

#if DEBUGCONSTRUCTOR
  printf("Construction finished: singleParticleEnergies\n\n");
#endif
}

// Allocate memory to the matrix of wave functions
void singleParticleEnergies:: allocateMemory ()
{
#if COMMENT
  printf("---Start allocating memory for the energies\n");
#endif
  int NbStates= States->readTotalNbStates();
  Range r(0,NbStates-1);
  spEnergies.resize(r);// as many single particle energies as single particle states

#if COMMENT
  printf("---Memory allocated for the energies\n");
#endif
}


// Destructor
singleParticleEnergies:: ~singleParticleEnergies ()
{
  deallocate();
}

// Free dynamic memory
void singleParticleEnergies:: deallocate ()
{

}

// Compute the single particle energies
void singleParticleEnergies:: computeEnergies ()
{
  int dim=States->Dimension();
  int NbStates= States->readTotalNbStates();
  // quantum numbers
  int n=0;
  int l=0;
  int ml=0;
  int ms=0;
  int index_ml=States->readIndex_ml();
  int index_ms=States->readIndex_ms();
  
  // 2D CASE
  if(dim==2)
    {
      // if basis_type = Harmonic Osc. basis set
      if(Basis->basis_type == 1) 
	{
	  for(int s=0; s<NbStates;s++)
	    {
	      n=States->readQuantumNumberOfState(s,0);
	      ml=States->readQuantumNumberOfState(s,index_ml);
	      ms=States->readQuantumNumberOfState(s,index_ms);
	      spEnergies(s)=2*n+abs(ml)+1;
	      }
	}
      // if basis_type = hydrogen like basis set 
      else if(Basis->basis_type == 2)
	{
	  printf("\n Hydrogen like basis set not supported yet\n");
	  exit(1);
	}
      else //  if no defined basis set
	{
	  printf("\nBasis set not correctly defined\n");
	  exit(1);
	}
  
    }// end of 2D case
  
  // 3D CASE
  else if(dim==3)
    {
      // if basis_type = Harmonic Osc. basis set
      if(Basis->basis_type == 1) 
	{
	  for(int s=0; s<NbStates;s++)
	    {
	      n=States->readQuantumNumberOfState(s,0);
	      l=States->readQuantumNumberOfState(s,1);
	      ml=States->readQuantumNumberOfState(s,index_ml);
	      ms=States->readQuantumNumberOfState(s,index_ms);
	      spEnergies(s)=2*n+l+1.5;
	      }
	}
      // if basis_type = hydrogen like basis set 
      else if(Basis->basis_type == 2)
	{
	  printf("\n Hydrogen like basis set not supported yet\n");
	  exit(1);
	}
      else //  if no defined basis set
	{
	  printf("\nBasis set not correctly defined\n");
	  exit(1);
	}
    }// end of 3D case

#if CHECK2
  cout << "\nSingle particle energies:\n" << spEnergies << endl;
#endif

}// end of  computeEnergies()



// Read the energy of one single particle in the given basis
double singleParticleEnergies:: read_Energy(const int indexState) const
{
  double energy = spEnergies(indexState);
  return energy;
}






#endif // singleParticleEnergies_CPP
