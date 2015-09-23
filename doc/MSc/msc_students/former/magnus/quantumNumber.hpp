/*
 * Code written by Magnus Pedersen Lohne
 * m.p.lohne@fys.uio.no
 * 
 * University of Oslo, August 2009
 *
 */


#ifndef QUANTUMNUMBER_HPP
#define QUANTUMNUMBER_HPP
#include <iomanip>


/*
 * QUANTUM NUMBER CLASS
 * base class
 *
 * Class for handeling quantum numbers of closed shell systems where the interaction
 * does not affect the total spin and total angular momentum  
 */

class quantumNumber {
  
public:
  int dim;             // dimensions
  int nbPart;          // number of particles in the system
  int nbShells;        // number of included shells in the single-particle wavefunction expansion
  int nbBasis;         // number of basis functions in the single-particle wavefunction expansion
  int fermiShell;      // fermi-shell
  int fermiState;      // fermi-state
  int nbParticleStates;// # of particle states (unoccupied states)
  int nbHoleStates;    // # of hole states (particle states)
  int nbChannels;      // # of channels, i.e. number of different set of total angular and spin for two-particle product state |ab>
  /*
   * table of states: containg all basis functions, i.e. all distinct sets of quantum numbers.
   * in two dimensions:   q1 q2 q3 q4  
   *             |1> ->  [.][ ][.][.]
   *             |2> ->  [.][.][.][.]
   *             |.> ->  [.][.][.][.]
   *              ...      .......
   */
  int** tableStates;
  /*
   * table of couples: containing all distinct couples of product states |ab> organized into channels;
   * total angular momentum and spin [M,S]. for given product state |ab>, total angular momentum and spin
   * M = m a + m_b
   * S = s_a + s_b
   * caracterize what channel (total angular momentum and spin) couple |ab> belong to. 
   * NB: |ab> and |ba> have obviously the same M and S, and hence we store only couples with distinct 
   * total quantum numbers a and b.
   * 
   *     channel      
   *      [M,S] -> [couple1][couple1][...][...]
   *       .... -> [couple1][couple2][...][...]
   *       ...  ->  ..........................
   * 
   * 
   *     [channel#]->[couple#]->[state1][state2]
   */
  int*** tableCouples;
  /*
   * number of coupels per channel
   */
  int* nbCouplesChannel;
  /*
   * table of channels [M][N]: table mapping channel# to total M and S
   */
  int** tableChannel; 
  /* 
   * table that map state couple |ab> to channel-index and couple-index
   */
  int*** stateMapChannel;
  /*
   * occupied states (i.e. hole states)
   */
  int* occStates;
  /*
   * unoccupied states (i.e. particle states)
   */
  int* unoccStates;   
  /*
   * constructor
   */
  //  quantumNumber(int, int, int, int);
  /*
   * set up table of states
   */
  virtual void setUp_tableStates() = 0;
  /*
   * set up table of occupied states
   */
  virtual void setUp_occStates() = 0;
  /* 
   * function that re-arrange the table of quantum numbers
   * into blocks of same angular and spin, but with different
   * radial quantum number n
   */
  virtual void reArr_tableStates() = 0;
  /*
   * set up table of couples
   */
  virtual void setUp_tableCouples() = 0;
  /*
   * set up table of mapping from couple |ab> to channel# and couple# 
   */
  virtual void setUp_stateMapChannel() = 0;
  /*
   * function that returns the bra-channel and the ket-channel
   * for incoming <alpha beta| and |gamma delta>, respectively
   */
  virtual void getChannels(int, int, int, int, int&, int&) = 0;
  /*
   * function that returns the bra-couple # and ket-couple #
   * for incoming <alpha beta| and |gamma delta>
   */
  virtual void getCouples(int, int, int, int, int&, int&) = 0;
  /*
   * write table of states to screen
   */
  virtual void print_tableStates() = 0;
  /*
   * write table of occupied states to screen
   */
  virtual void print_occStates() = 0;
  /*
   * write table of unoccupied states to screen
   */
  virtual void print_unoccStates() = 0;
  /*
   * destructor
   */
  virtual ~quantumNumber(){
    delete[] tableStates;
    delete[] tableCouples;
    delete[] nbCouplesChannel;
    delete[] tableChannel;
    delete[] stateMapChannel;
    delete[] occStates;
  }
};



/*
 * QDOT QUANTUM NUMBER CLASS
 * subclass of quantumNumber
 *
 */

class qdotQuantumNumber : public quantumNumber {

private:
  int n_max;      // maximal value of n
  int m_max;      // maximal value of m
  int index_n;    // radial quantum number index
  int index_m;    // angular (z-projection) quantum number index
  int index_s;    // spin quantum number index

public: 
  /*
   * constructor  
   */
  qdotQuantumNumber(int, int, int, int);
  /*
   * set up table of states
   */
  void setUp_tableStates();
  /*
   * function writes included shells with corresponding quantum numbers to screen
   */
  void printShells();
  /*
   * set up table of occupied states
   */
  void setUp_occStates();
  /*
   * set up table of unoccupied states
   */
  void setUp_unoccStates();
  /* 
   * function that re-arrange the table of quantum numbers
   * into blocks of same angular and spin, but with different
   * radial quantum number n
   */
  void reArr_tableStates();
  /*
   * function that prints all basis functions to screen
   */
  void print_tableStates();
  /*
   * set up table of couples
   */
  void setUp_tableCouples();
  /*
   * function that writes table of couples to screen
   */
  void print_tableCouples();
  /*
   * write occupied states to screen
   */
  void print_occStates();
  /*
   * write unoccupied states to screen
   */
  void print_unoccStates();
  /*
   * set up table of mapping from couple |ab> to channel# and couple# 
   */
  void setUp_stateMapChannel();
  /* 
   * function that writes table of mapping info to screen
   */
  void print_stateMapChannel();
   /*
   * function that returns the bra-channel and the ket-channel
   * for incoming <alpha beta| and |gamma delta>, respectively
   */
  void getChannels(int, int, int, int, int&, int&);
  /*
   * function that returns the bra-couple # and ket-couple #
   * for incoming <alpha beta| and |gamma delta>
   */
  void getCouples(int, int, int, int, int&, int&);
  /*
   * destructor
   */
  ~qdotQuantumNumber(){}

};


#endif 
