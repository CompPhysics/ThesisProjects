#include <iostream>
#include "singleParticleOrbitals.h"
#include "singleParticleWaveFunctions.h"
#include "singleParticleEnergies.h"

using namespace std;

int main()
{
   // parameters later to be read from files
  int Nlevel=2;  // Number of energy level taken into account
  int dim=2;   // spatial dimension of the system


  // Prompt the user to enter the parameter one by one
  cout << "Enter the dimension of the system {2 or 3}: \n";
  cin >> dim;
  cout << "Enter the number of the maximum energy level to consider: \n";
  cin >> Nlevel;

  // generate all possible single particles states
  singleParticleOrbitals States(dim, Nlevel);
   
  printf("you choose a basis with dimension %d and the max energy level E=%d, so %d possible states \n", States.Dimension(), States.maxEnergylevel(), States.update_totalNbStates());

  // for a type of basis set, compute all single particle wave fucntions corresponding to each single particle states
  singleParticleWaveFunctions Basis(&States);
  
  // compute all single particle energies for a given type of basis set
  singleParticleEnergies sp_energies(&States, &Basis);




  return 0;
}
