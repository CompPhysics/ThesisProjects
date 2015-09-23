#ifndef HFalgorithm_h_IS_INCLUDED
#define HFalgorithm_h_IS_INCLUDED

#include <iostream>
#include "singleParticleOrbitals.h"
#include "singleParticleWaveFunctions.h"
#include "singleParticleEnergies.h"
#include "CoulombMatrix.h"
#include <blitz/array.h>

// file        : HFalgorithm.h
// description : definition of class HFalgorithm

using namespace std;
using namespace blitz;

class HFalgorithm
{
  
 private:
  const singleParticleOrbitals * States;            // Pointer to the one-body HO states (quantum numbers)
  const singleParticleWaveFunctions * Basis;        // Pointer to the one-body HO energy eigenstates
  const singleParticleEnergies * BasisEnergies;     // Pointer to the one-body HO energy eigenvalues
  const CoulombMatrix * CMatrix;                    // Pointer to the Coulomb Matrix object
  int nbS;                                          // Number of states in the basis set
  Array<double,1> ParticleEnergies;                 // Single Particle energies
  Array<double,2> HForbitalsCoeff;                  // Coeff. of the HF orbitals expanded in a given basis
  Array<double,2> EffectivePotential;               // Effective potential for a given set of HF orbitals
  double** A;                                       // Matrix of the eigenvalue problem
  double TotalEnergy;                               // total energy of the system
  double epsilon;                                   // Precision used for self-consistency
 
public: 
  HFalgorithm (const singleParticleOrbitals* States, const singleParticleWaveFunctions* Basis, const singleParticleEnergies * energyEigenValues, const CoulombMatrix * CM, const double epsi);             // CONSTRUCTOR
  ~HFalgorithm ();                           // Destructor
  void allocateMemory ();                    // Allocate memory to the matrix of wave functions
  void deallocate ();                        // Free dynamic memory
  void initialization ();                    // Initialisation of the single particle energies/wave functions
  void compute_EffectivePotential ();        // update the effective potential with the optimized single particle wfs
  void solve_HFequations ();                 // Solve the linear system given by HF equations
  double compute_TotalEnergy();              // Update the total energy of the system
  void check_selfConsistency();              // Check for self consistency
  void normalize_HForbitals();               // normalization of the vectors of coefficient describing HF orbitals

  // update the value of the precision used in the Hartree-Fock algorithm (Energy precision)
  void setEpsilon(const double epsil);
  void printScreen_A() const;                // Print the matrix of the eigenvalue problem to screen
  void save2file(ostream & mFile);           // Write the eigenvalue problem to file in matlab format
  void saveCoeff2file(ostream & mFile);      // Write only the HF orbitals coeff. to file in matlab format
  void saveEffPot2file(ostream & mFile);     // Write only the Effective Potential matrix to file in matlab format
  void saveEnergy2file(ostream & mFile);     // Write only the Total Energy to file in matlab format

}; // end of class HFalgorithm

#endif // HFalgorithm_h_IS_INCLUDED
