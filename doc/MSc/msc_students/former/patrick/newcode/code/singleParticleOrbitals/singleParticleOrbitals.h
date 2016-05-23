#ifndef singleParticleOrbitals_h_IS_INCLUDED
#define singleParticleOrbitals_h_IS_INCLUDED

#include <iostream>
#include <iomanip>

// file        : singleParticleOrbitals.h
// description : definition of class singleParticleOrbitals: used to generate and store all possible states of the closed shell system

using namespace std;


class singleParticleOrbitals
{
 private:
  int dim;                 // space dimensions {2 or 3}, safer to let it private parameter, like that you can't modify it from outside the class
  int index_n;             // nodal or principal quantum number
  int index_l;             // azimuthal quantum number (angular momentum)
  int index_ml;            // magnetic quantum number, (projection of angular momentum)
  int index_ms;            // spin projection quantum number
  int Mmax;                // biggest total angular momentum considering 2 particles
  
  int totalNbStates;       // gives the total number of states {n,l,m} summing over all the energy levels up to Nmax
  int NbChannels;          // Number of channel={M,S}
  int TotNbCouples;        // Total number of states couples |ij> with requirement i<j
  int Nmax;                // level of the maximum energy shell considered for the HO oscillator basis set, also called "energy cut" [Simen]
  int ** TableOfStates;    // table with the list of all single particle states (each state defined by a set of quantum number {n,l,ml,ms}) inside the basis
  int * lengthTabInteractionStates;// number of couples per channel
  int *** TabInteractionStates;    // list all the combinations of state couple with identical M and S such that {M,S}=channel(2M+S)
  int nbBlocks;            // # of blocks of states with either same {l, ml} in 3D or {ml} in 2D
  int * sizeBlock;         // # of states per block with the same {l, ml} in 3D or {ml} in 2D 

 public:                // members visible also outside the class
  singleParticleOrbitals (int _dim, int _Nmax);// singleParticleOrbitals S(2,8);
  ~singleParticleOrbitals ();               // Destructor
  singleParticleOrbitals& operator= (const singleParticleOrbitals& Sb);       // Sa = Sb;
  void deallocate (); // Free dynamic memory

  int Dimension () const;                      // int basisDimension = S.Dimension();
  int readIndex_l () const;                    // int index_l=readIndex_l();
  int readIndex_ml () const;                   // int index_ml=readIndex_ml();
  int readIndex_ms () const;                   // int index_ms=readIndex_ms();
  int maxEnergylevel () const;                 // int Nlevel = S.maxEnergyLevel();

  int readQuantumNumberOfState(const int indexState, const int indexNumber) const; //ex.: n =Basis.readQuantumNumberOfState(index_state,0); ml=Basis.readQuantumNumberOfState(index_state,Basis.readIndex_ml());
  int readTotalNbStates () const;              // int totNbStates = S.totalNbStates();
  int readNbChannels () const;                 // int NbChannels = S.readNbChannels();
  int readTotNbCouples () const;               // int totNbCouples = S.readTotNbCouples();
  int mappingCoupleToSingleStates (const int channel, const int indexCouple, const int leftRightState) const; // given ch=2M+S, and the index of a couple in this channel, read each Single particle states out of the table.
  int * read_NbCouplesPerChannel () const;     // int * NbCouplesPerChannel = S.read_NbCouplesPerChannel();
  int update_totalNbStates ();                 // int totalNbStates = S.totalNbStates();
  void update_Tables () ;                      // update_Tables();

  // given 2 single particle states, compute the corresponding channel index ch={2M+S}
  int read_channelOfCouple (const int indexState1, const int indexState2) const;  
  
  // given 2 single particle states, find the corresponding index for the couple.
  int mappingSingleStatesToCouple (const int indexState1, const int indexState2) const; 
  
  // sort the states in the TableOfStates with increasing values of the chosen quantum number
  void sort_TableOfStates (const int indexQuantumNumber); 
  
  // sort the states in the TableOfStates with increasing values of the chosen quantum number 1 first, then for a given value of the qn1, sort with increasing values of the qn2
  void sort_TableOfStates (const int indexQn1, const int indexQn2);
  
  // find the number of blocks with the same QN and the nb of states inside each block
  void update_infoBlocks (const int indexQuantumNumber); 

  // find the number of blocks with the same QN1 and QN2 and the size of each block (the nb of states) 
  void update_infoBlocks (const int indexQuantumNumber1,const int indexQuantumNumber2); 
  
  // print to screen the list of possible single particle states 
  void print_TableOfStates ();               

  // Write to file the tables of states
  void save2file(ostream & mFile);
};

#endif
