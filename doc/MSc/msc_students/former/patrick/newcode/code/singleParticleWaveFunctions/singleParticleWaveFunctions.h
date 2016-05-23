#ifndef singleParticleWaveFunctions_h_IS_INCLUDED
#define singleParticleWaveFunctions_h_IS_INCLUDED

#include <iostream>
#include "singleParticleOrbitals.h"
#include <blitz/array.h>

// file        : singleParticleWaveFunctions.h
// description : definition of class singleParticleWaveFunctions

// NOT FULLY IMPLEMENTED HERE SINCE IT DO NOT ACTUALLY COMPUTE THE POLYNOMIALS USED FOR THE WAVE FUNCTIONS, but... IT HAS BEEN BUILD TO MAKE A LINK TO THE SingleParticleEnergies Class, and so to be easy to modify to generalize it later.

using namespace std;
using namespace blitz;

class singleParticleWaveFunctions
{
 public:
  enum basis_list { HarmonicOscillatorBasis=1, HydrogenLikeBasis=2 }; // List of possible choice of basis sets of wave functions
  basis_list basis_type;                      // Which basis polynomials to use: HarmonicOscillatorBasis or HydrogenLikeBasis

 private:
  const singleParticleOrbitals * States;      // Pointer to the single particle states
  int MaxStep;                                // Maximum size of a wave function (corresponding to the radius of the QD if it had dimensions of a length)
  Array<double,1> WFs;                        // Matrix of the wave functions
  
  

public: 
  singleParticleWaveFunctions (const singleParticleOrbitals* S);             // singleParticleWaveFunctions S();
  ~singleParticleWaveFunctions ();           // Destructor
  void allocateMemory ();                    // Allocate memory to the matrix of wave functions
  void deallocate ();                        // Free dynamic memory
  void setBasisType(const basis_list x);     // Set basis type. Possibility are: HarmonicOscillatorBasis or HydrogenLikeBasis
  void setMaxStep();
  //basis_list read_basisType();               // Read the type of basis defined for basis set

}; // end of class singleParticleWaveFunctions

#endif // singleParticleWaveFunctions_h_IS_INCLUDED
