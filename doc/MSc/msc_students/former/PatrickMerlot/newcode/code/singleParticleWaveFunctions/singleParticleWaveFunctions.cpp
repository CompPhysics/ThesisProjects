#ifndef singleParticleWaveFunctions_CPP
#define singleParticleWaveFunctions_CPP

// file        : singleParticleWaveFunctions.cpp
// description : implentation of class singleParticleWaveFunctions

#include <iostream>
#include <cmath>
#include <cstdlib>
#include "singleParticleWaveFunctions.h"

using namespace std;

// Contructors
singleParticleWaveFunctions :: singleParticleWaveFunctions(const singleParticleOrbitals* S):basis_type(HarmonicOscillatorBasis), States(S)
{
#if DEBUGCONSTRUCTOR
  printf("\nConstruction starts: singleParticleWaveFunctions\n");
#endif
 
 allocateMemory();

#if DEBUGCONSTRUCTOR
  printf("Construction finished: singleParticleWaveFunctions\n");
#endif
}

// Allocate memory to the matrix of wave functions
void singleParticleWaveFunctions:: allocateMemory ()
{
#if COMMENT
  printf("---Start allocating memory for the Wave functions\n");
#endif
  int NbStates= States->readTotalNbStates();
  WFs.resize(2);// resize here normally as NbStates x MaxStep, but need MaxStep to be initialized

#if COMMENT
  printf("---Memory allocated for the wave functions\n");
#endif
}


// Destructor
singleParticleWaveFunctions:: ~singleParticleWaveFunctions ()
{
  deallocate();
}

// Free dynamic memory
void singleParticleWaveFunctions:: deallocate ()
{

}

// Set basis type. Possibility are: HarmonicOscillatorBasis or HydrogenLikeBasis
void singleParticleWaveFunctions:: setBasisType(const basis_list x)   
{
  basis_type = x;
}

// // Read the type of basis defined for basis set
// basis_list singleParticleWaveFunctions:: read_basisType()
// {
//   return basis_type;
// }


#endif // singleParticleWaveFunctions_CPP
