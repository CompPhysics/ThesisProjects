#ifndef FuncUpDown_IS_INCLUDED
#define FuncUpDown_IS_INCLUDED

#include <cmath>
#include <iostream>
#include <fstream>
#include "../Coor/Coor.h"
#include "../SingleParticleFuncs/SingleParticleFuncs.h"
#include "../Domain/Domain.h"
#include "../Func/Func.h"

// ****************************************************************
// *                           FUNCUP                             *
// ****************************************************************
template <class Param>
class FuncUp : public FuncSetMultivar<Param> {
 protected:
  char* fixedParamsUp;

 public:
  FuncUp() {};
  virtual ~FuncUp() 
    {  
      for (int i=0; i<len; i++)
	delete function[i];
      delete function;
      delete[] values;
    }

  FuncUp(Param& coordinate, int __len) : 
    FuncSetMultivar<Param>(coordinate, __len) {}

  virtual void init(Param& coordinate, int __len, Domain* domain, 
		    int alphaVar) 
    { 
      fixedParamsUp   = domain->getFixedParamsUp();
      FuncSetMultivar<Param>::init(coordinate, __len);
      init(domain, alphaVar);
    }

  virtual void init() {
    FuncSetMultivar<Param>::init();
  }

  virtual void init(Domain* domain, int alphaVar) {
    int index = 0;
    ifstream ifile;
    ifile.open(fixedParamsUp);
    if (domain->getOrbitalType() == "Hydrogen") {
      if (domain->getUp1s())  {function[index]=new Hydr1s<Param>;  index++;}
      if (domain->getUp2s())  {function[index]=new Hydr2s<Param>;  index++;}
      if (domain->getUp2px()) {function[index]=new Hydr2px<Param>; index++;}
      if (domain->getUp2py()) {function[index]=new Hydr2py<Param>; index++;}
      if (domain->getUp2pz()) {function[index]=new Hydr2pz<Param>; index++;}
      if (index != len ) 
	cerr << "Error creating FuncUp!\nNot matching dimensions of number " 
	     << "particles spin up (" << len 
	     << ") and number of input functions\n";
      for (int i=0; i<len; i++) {
	function[i]->attach(domain->getNumAlpha(), 
			    domain->getAlphaParam(alphaVar));
      }
    }
    else if (domain->getOrbitalType() == "HartreeFock") {
      for (int i=0; i<len; i++) {
	function[i] = new HF<Param>;
	function[i]->readHFparams(ifile);
      }
    }
    for (int i=0; i<len; i++)
      function[i]->setNumberDimensions(domain->getNumDimensions());
    ifile.close();
  }
};


// ****************************************************************
// *                          FUNCDOWN                            *
// ****************************************************************
template <class Param>
class FuncDown : public FuncSetMultivar<Param> {
 protected:
  char* fixedParamsDown;
  
 public:
  FuncDown() {};
  virtual ~FuncDown()
    {  
      for (int i=0; i<len; i++)
	delete[] function[i];
      delete[] values;
    }
  FuncDown(Param& coordinate, int __len)  : 
    FuncSetMultivar<Param>(coordinate, __len) {}

  virtual void init(Param& coordinate, int __len, 
		    Domain* domain, int alphaVar) 
    {
      fixedParamsDown = domain->getFixedParamsDown();
      FuncSetMultivar<Param>::init(coordinate, __len);
      init(domain, alphaVar);
    } 
  
  virtual void init() {
    FuncSetMultivar<Param>::init();
  }
  
  virtual void init(Domain* domain, int alphaVar) {
    int index = 0;
    ifstream ifile;
    ifile.open(fixedParamsDown);
    if (domain->getOrbitalType() == "Hydrogen") {
      if (domain->getDown1s())  {function[index]=new Hydr1s<Param>;  index++;}
      if (domain->getDown2s())  {function[index]=new Hydr2s<Param>;  index++;}
      if (domain->getDown2px()) {function[index]=new Hydr2px<Param>; index++;}
      if (domain->getDown2py()) {function[index]=new Hydr2py<Param>; index++;}
      if (domain->getDown2pz()) {function[index]=new Hydr2pz<Param>; index++;}
      if (index != len ) cerr << "Error creating FuncUp!\nNot matching"  
			      << " dimensions of number particles spin up (" 
			      << len << ") and number of input functions\n";
      for (int i=0; i<len; i++)
	function[i]->attach(domain->getNumAlpha(), 
			    domain->getAlphaParam(alphaVar));
    } 
    else if (domain->getOrbitalType() == "HartreeFock") {
      for (int i=0; i<len; i++) {
	function[i] = new HF<Param>;
	function[i]->readHFparams(ifile);
      }
    }
    for (int i=0; i<len; i++)
      function[i]->setNumberDimensions(domain->getNumDimensions());
    ifile.close();
  }
};

#endif
