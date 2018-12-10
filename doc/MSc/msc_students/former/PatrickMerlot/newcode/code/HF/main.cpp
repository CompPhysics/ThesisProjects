#include <iostream>
#include <fstream>
#include "singleParticleOrbitals.h"
#include "CoulombMatrix.h"
#include "singleParticleWaveFunctions.h"
#include "singleParticleEnergies.h"
#include "HFalgorithm.h"
#include "ConfigFile.h"

#include <cstring>
#include <cstdio>


// file        : main.cpp
// description : Definition and solving energy of a quantum dot using Hartree-Fock method

using namespace std;
using std::string;
using std::cout;
using std::endl;

// declaration of functions
void cmd(char **, int &, int &, double &,double &, string &);

// Open output file
void setMatlabOutputFile(const string, ofstream &);

int main(int argc, char* argv[])
{
  // parameters later to be read from files
  int Nlevel;                       // Maximum level of energy taken into account, which determines the number of electrons in the closed-shell system
  int dim;                          // spatial dimension of the system
  double lambda;                    // characterizes the strength of the confinement
  double epsilon;                   // precision parameter for the self-consistency of the total energy
  std::string filename, nomF;             // Filename for matlab output.
  std::ofstream matlab;             // Stream for matlab output.
  
  // read input parameters from file 
  ConfigFile config( "parameters.inp" ); // load the configuration file
  dim = config.read<int>("dimension", 2);  // read input parameters or initialize them if missing in the input file
  Nlevel= config.read<int>("R", 1);
  lambda= config.read<double>("lambda", 1.0);
  epsilon= config.read<double>("epsilon", 1e-15);
  filename= config.read<string>("filename", "Qdot.m");
  bool save_states = config.read<bool>("save_states", false);   
  bool save_CoulombMatrix = config.read<bool>("save_CoulombMatrix", false); 
  bool save_EigValPb = config.read<bool>("save_EigValPb", false); 

  //overwrite them using the command line
  cmd(argv, dim, Nlevel, lambda, epsilon, filename);
  
  // Open MATLAB output file.
  // filename  File name of MATLAB script.
  setMatlabOutputFile(filename, matlab);

  // Prompt the parameters
  printf("\nInput Parameters:\n");
  cout << "Dimension of the system: d=" << dim << endl;
  cout << "Maximum energy level to consider: R=" << Nlevel << endl;
  cout << "Strength of confinement: lambda=" << lambda << endl;
  cout << "Precision for the HF algorithm: epsilon=" << epsilon << endl;
  cout << endl;

  // Print paramters to file
  if ( matlab.is_open() )
    {
      matlab << "     %%%%%%%%% Input Parameters %%%%%%%%   " << endl;
      matlab << "d= "<< dim;
      matlab << "           %%% Dimension of the system" << endl;
      matlab << "R= " << Nlevel;
      matlab << "           %%% Maximum energy level to consider" << endl;
      matlab.precision(5);
      matlab.setf(ios::fixed,ios::floatfield); matlab << "lambda= " << lambda ;
      matlab.setf(ios::scientific); matlab << "           %%% Strength of confinement" << endl;
      matlab << "epsilon= " << epsilon;
matlab << "           %%% Precision of the HF algorithm" << endl << endl << endl;
    }

  // CREATE TABLES OF QUANTUM STATES IN THE CLOSED SHELL SYSTEM
  singleParticleOrbitals States(dim, Nlevel);
    
#if COMMENT
  printf("Dimension %d and %d energy levels, so %d possible states \n", States.Dimension(), States.maxEnergylevel(), States.update_totalNbStates());
#endif
  
  // for a type of basis set, compute all single particle wave fucntions corresponding to each single particle states
  singleParticleWaveFunctions Basis(&States);
  
  
  // compute all HO single particle energies for a given type of basis set
  singleParticleEnergies sp_energies(&States, &Basis);

  // GENERATES THE MATRIX OF THE COULOMB INTERACTIONS
  if(dim == 2) // 3D to be implemented otherwise if no analytical expression
    {
      // Compute the matrix of the Coulomb interactions for a given value lambda of the confinement strength
      CoulombMatrix CM(&States,lambda);
  
      // HARTREE-FOCK ALGORITHM
      // --- allocate memory for matrices in the algorithm
      HFalgorithm algo(&States, &Basis, &sp_energies , &CM, epsilon);
      // --- guess the value of the initial coefficients for the HF orbitals and compute the first effective potential
      algo.initialization();
      
      // --- Self-consistent method
      double TotalEnergy = 0.0;
      double old_TotalEnergy = -1;
      int count_iter = 0;
      cout << "\n\n ### Iterative calculation of the total energy ###\n";
      double deltaEnergy=abs(old_TotalEnergy-TotalEnergy);
      while( deltaEnergy >= epsilon)
      //while(true)
	{
	  count_iter++;
	  cout << "iter #" << setfill ('0') << setw (3) << count_iter << ":   ";
	  old_TotalEnergy = TotalEnergy;
	  // solve the HF equations by solving an eigenvalue problem and update the HF orbitals coeff.
	  algo.solve_HFequations();

	  // compute the total energy of the system
	  TotalEnergy = algo.compute_TotalEnergy();
	  deltaEnergy=abs(old_TotalEnergy-TotalEnergy);
	  cout << "  DeltaEnergy=" << deltaEnergy << endl;
	  
	  // update the Effective potential with the new coeff.
	  algo.compute_EffectivePotential();
	  
	  // print to file the energy at each iteration
	  algo.saveEnergy2file(matlab);
	  matlab << "% iter #";
	  matlab << setfill ('0') << setw (3) << count_iter << " ";
	  matlab << endl;
	}// end of iterations

    
      // Prompt the parameters
      cout.precision(5);
      printf("\n\n\nRecap. Input Parameters:\n");
      cout << "Dimension of the system: d=" << dim << endl;
      cout << "Maximum energy level to consider: R=" << Nlevel << endl;
      cout << "Strength of confinement: lambda=" << lambda << endl;
      cout.setf(ios::scientific); cout << "epsilon= " << epsilon << endl;
      cout << "# particles = # states =" << States.update_totalNbStates() << endl;
      cout << "output filename = " << filename << endl;
      cout << endl;

      // Dump results to file
      if ( matlab.is_open() )
	{
	  matlab.setf(ios::fixed,ios::floatfield);     // floatfield not set
	  matlab.precision(20);
	  matlab << endl << endl <<"E" << " = " << TotalEnergy << "     % Hartree-Fock energy of the quantum dot"<< endl;
	  
	  // print to file the coefficients and the effective potential
	  if(save_EigValPb)
	    algo.save2file(matlab);
	  
	  // print the Coulomb Matrix to file
	  if(save_CoulombMatrix)
	    CM.save2file(matlab);
	  
	  // print the tables of states to file  
	  if(save_states)
	    States.save2file(matlab);
	}
      else
	{
	  printf("\nERROR: Output file not opened\n\n");
	  exit(1);
	}

    }// end if(dim == 2)

  return 0;
}// END OF MAIN


// Function to read out from Command line the parameters not set to default values
void cmd(char ** args, int & dim, int & Nlevel, double & lambda, double & epsilon, string & filename)
{
  for(char ** k = args+1, ** v = args+2; *k && *v; (k += 2), (v += 2))
    {
      if(strcmp(*k, "d") == 0 || strcmp(*k, "dim") == 0)
	dim = atoi(*v);
      else if(strcmp(*k, "Nlevel") == 0 || strcmp(*k, "R") == 0)
	Nlevel = atoi(*v);
      else if(strcmp(*k, "lambda") == 0 || strcmp(*k, "l") == 0)
	lambda  = atof(*v);
      else if(strcmp(*k, "epsilon") == 0 || strcmp(*k, "e") == 0)
	epsilon  = atof(*v);
      else if(strcmp(*k, "filename") == 0 || strcmp(*k, "file") == 0 || strcmp(*k, "o") == 0 || strcmp(*k, "output") == 0)
	{
	  filename  = *v;
	  char * outputFile = new char [filename.size()+1];
	  strcpy (outputFile, filename.c_str()); // outputFile now contains a c-string copy of filename
	  filename.assign(outputFile);
	}
      else fprintf(stderr, "Unrecognized key %s with value %s ignored\n", *k, *v);
    }
}//end of cmd()

// Open output file
void setMatlabOutputFile(const string filename, ofstream & matlab)
{
  matlab.open(filename.c_str());//fopen (filename.c_str(),"w");
}
