#ifndef singleParticleEnergies_h_IS_INCLUDED
#define singleParticleEnergies_h_IS_INCLUDED

#include <iostream>
#include "singleParticleOrbitals.h"
#include "singleParticleWaveFunctions.h"
#include <blitz/array.h>

// file        : singleParticleEnergies.h
// description : definition of class singleParticleEnergies

// FIRST, do not integrate anything but just compute the single particle energies in the Harmonic oscillator basis set, so in 2D E(n,ml)=2n+|ml| and in 3D E(n,l,ml)=2n+l

using namespace std;
using namespace blitz;

class singleParticleEnergies
{
  
 private:
  const singleParticleOrbitals * States;      // Pointer to the single particle states
  const singleParticleWaveFunctions * Basis;  // Pointer to the single particle wave functions
  Array<double,1> spEnergies;                 // Matrix of the wave functions

  

public: 
  singleParticleEnergies (const singleParticleOrbitals* States, const singleParticleWaveFunctions* Basis);             // CONSTRUCTOR
  ~singleParticleEnergies ();                // Destructor
  void allocateMemory ();                    // Allocate memory to the matrix of wave functions
  void deallocate ();                        // Free dynamic memory
  void computeEnergies();                    // integrate the hamiltonian for a given potential
  double read_Energy(const int indexState) const;  // Read the energy of one single particle in the given basis

}; // end of class singleParticleEnergies

#endif // singleParticleEnergies_h_IS_INCLUDED
