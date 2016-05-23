#ifndef DMC_importance_hpp_IS_INCLUDED
#define DMC_importance_hpp_IS_INCLUDED

// class doing diffuse monte carlo
// Usage:
// make an object of DMC (with infile and outfil name , 
// call the function diffMC(), and make sure the object 
// is deleted.

// standard library headers
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <new>

// other headers
#include "QickArray/QickArray.hpp"
#include "random/random.hpp"
#include "walker.hpp"

class DMC{

private:

  // copy walker k to walker l
  void copy_walker(int k, int l);
  
  // make space if no space
  int make_room(int ni, int nj, int nk);

  // x squared
  inline double sqr(double x){return x*x;}

  void initialize(std::string filename);

  void output(std::string filename);

  void oneMonteCarloStep(Random &ran, int k);

  void oneTimeStep(Random &ran);

  int no_of_walkers, desired_walkers, max_walkers, particles, 
    dimensions, steps, termalization, step_length;

  long idum;

  double tau, D, e_trial, std_dev;

  double eav, energy, energy2;

  int accept, trials;

  QickArray params;

  std::string infile, outfile;

  Walker *walkers;

public:

  DMC(std::string infile_, std::string outfile_);
  // destructor, note: output is called here
  ~DMC(){ output(outfile); delete[] walkers; }

  void diffMC();

};


inline DMC::DMC(std::string infile_, std::string outfile_){
  infile=infile_;
  outfile=outfile_;
  initialize(infile);
  walkers = new Walker[no_of_walkers](particles);
}

#endif
