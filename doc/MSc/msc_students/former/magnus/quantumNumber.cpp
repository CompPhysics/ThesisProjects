/*
 * Code written by Magnus Pedersen Lohne
 * m.p.lohne@fys.uio.no
 * 
 * University of Oslo, August 2009
 *
 */


#include <iomanip>
#include <iostream>
#include <fstream>
#include <cmath>
#include "quantumNumber.hpp"



/*************************************
 * class functions for quantumNumber
 ************************************/

qdotQuantumNumber::qdotQuantumNumber(int nbPart, int dim, int nbShells, int fermiShell){
  
  // initialize class variables
  index_n = 0;        // n quantum number index  
  index_m = 1;        // m quantum number index
  index_s = 2;        // s quantum number index
  this->dim = dim;
  this->nbPart = nbPart;
  this->nbShells = nbShells;
  this->fermiShell = fermiShell;
    
  // determine number of basis functions for given number of shells
  int temp,i;
  if(nbShells == 1){
    nbBasis = 2;
  } 
  else {
    nbBasis = 2; 
    temp = 2;
    for(i=2; i<=nbShells; i++){
      temp += 2;
      nbBasis += temp;
    }
  }

  // determine fermi state
  if(fermiShell == 1){
    fermiState = 1;
  }
  else {
    fermiState = 1;
    temp = 2;
    for(i=2; i<=fermiShell; i++){
      temp += 2;
      fermiState += temp;
    }
  }

  // determine number of particle states
  nbHoleStates = fermiState + 1;

  // determine number of hole states
  nbParticleStates = nbBasis - fermiState - 1;
  
  // checking if input values are valid
  if( (fermiState+1) != nbPart ){
    printf("Usage:\n");
    printf("The number of states below (fermi-state invluded) the fermi state must be equal to the number of particles\n");
    exit(1);
  }
  else if(nbBasis < nbPart){
    printf("Usage:\n");
    printf("The number of basis functions must be equal to or larger than the number of particles in the system\n");
    exit(1);
  }
  else if(fermiShell > nbShells){
    printf("Usage:");
    printf("The number of included shells in the basis must be larger than the fermi shell\n");
    exit(1);
  }
  else if(dim!=2){
    printf("Usage:");
    printf("The program is at moment just written for the two-dimensional case\n");
    exit(1);
  }
     
  // initialize table of states
  tableStates = new int*[nbBasis];
  for(int i=0; i<nbBasis; i++){
    tableStates[i] = new int[3];
  }

  // print all included shells (with quantum numbers) to screen
  // printShells();
  
  // set up table of states
  setUp_tableStates();

  // set up table of occupied states (i.e. hole states)
  setUp_occStates();

  // set up table of unoccupied states (i.e. particle states)
  setUp_unoccStates();

  // re-arrange table of quantum numbers into blocks of same angular momentum and spin
  reArr_tableStates();

  // set up table of couples 
  setUp_tableCouples();

  // set up table of mapping-info 
  setUp_stateMapChannel();

}



void qdotQuantumNumber::printShells(){
  printf("\n");
  printf("*** SHELLS ***\n");
  int N,n,m,s;
  for(N=1; N<=nbShells; N++){        
    printf("%i  --->  ",N);
    printf("n=");
    for(n=0;n<=(N-1)/2;n++){              
      printf("%i,",n);
    }
    printf(" m=");
    for(m=-N+1; m<=(N-1); m+=2){
      printf("%i,",m);
    }
    printf(" s=");
    for(s=-1; s<=1; s+=2){
      printf("%i,",s);
    }
    printf("\n");  
  }
} // end printShells



void qdotQuantumNumber::setUp_tableStates(){
  int N,n,m,s;
  int state = 0;  // state variable
  for(N=1; N<=nbShells; N++){              // loop over included shells
    for(n=0;n<=(N-1)/2;n++){               // loop over n quantum number 
      for(m=-N+1; m<=(N-1); m+=2){         // loop over m quantum number
	if( (2*n + abs(m)) == (N-1)){      // test whether the state n,m is in the given shell B
	  for(s=-1; s<=1; s+=2){           // loop over spin quantum number
	    tableStates[state][index_n] = n;
	    tableStates[state][index_m] = m;
	    tableStates[state][index_s] = s;
	    state++;  // update state variable
	  }
	}
      }
    }
  }
} // end setUp_tableStates



void qdotQuantumNumber::setUp_occStates(){
  occStates = new int[nbHoleStates];
  for(int i=0; i<nbHoleStates; i++){
    occStates[i] = i;
  }
} // end setUp_occStates



void qdotQuantumNumber::setUp_unoccStates(){
  unoccStates = new int[nbParticleStates];
  int state = fermiState + 1;
  for(int i=0; i<nbParticleStates; i++){
    unoccStates[i] = state;
    state++;
  }
}



void qdotQuantumNumber::print_tableStates(){
  printf("\n");
  printf("-- Table of states\n");
  for(int state=0; state<nbBasis; state++){
    printf("   |%i>  -->  n=%i m=%i s=%i \n",state,tableStates[state][index_n],tableStates[state][index_m],tableStates[state][index_s]);
  }
  printf("\n");
} // end print_tableStates
 


void qdotQuantumNumber::print_occStates(){
  printf("\n");
  printf("-- Table of occupied states (hole states)\n");
  for(int state=0; state<nbHoleStates; state++){
    printf("   |%i>\n",occStates[state]);
  }
  printf("\n");
}


void qdotQuantumNumber::print_unoccStates(){
  printf("\n");
  printf("-- Table of unoccupied states (particle states)\n");
  if(nbParticleStates == 0){
    printf("   none\n");
  } else {
    for(int state=0; state<nbParticleStates; state++){
      printf("   |%i>\n",unoccStates[state]);
    }
  }
  printf("\n");
}


void qdotQuantumNumber::reArr_tableStates(){
  
  // variable declaration
  int s,m,n,state;   

  // initial max values
  n_max = 0;
  m_max = 0;
  
  // dobble pointer to the old quantum number array
  int** old_tableStates = new int*[nbBasis];
  for(int i=0; i<nbBasis; i++){
    old_tableStates[i] = new int[3];
  }
  for(int i=0; i<nbBasis; i++){
    for(int j=0; j<3; j++){
      old_tableStates[i][j] = tableStates[i][j];
    }
  }

  // determine maximum value of n and m
  for(state=0; state<nbBasis; state++){
    if(old_tableStates[state][index_n]>n_max){
    n_max = old_tableStates[state][index_n];
    }
    if(old_tableStates[state][index_m]>m_max){
      m_max = old_tableStates[state][index_m];
    }
  }

  // sort the states in blocks with equal m and s
  int index = 0;
  int index2 = 0;
  int index3 = 0;
  for(s=-1; s<=1; s+=2){
    for(m=-m_max; m<=m_max; m++){
      for(n=0; n<=n_max; n++){
	for(state=0; state<nbBasis; state++){
	  if(old_tableStates[state][index_n] == n && old_tableStates[state][index_m] == m && old_tableStates[state][index_s] == s){
	    tableStates[index][index_n] = old_tableStates[state][index_n];
	    tableStates[index][index_m] = old_tableStates[state][index_m];
	    tableStates[index][index_s] = old_tableStates[state][index_s];
	    if(state<nbHoleStates){
	      occStates[index2] = index;
	      index2++;
	    }
	    else {
	      unoccStates[index3] = index;
	      index3++;
	    }
	    index++;
	  }
	}
      }
    }
  }
  
  // free memory
  delete[] old_tableStates;
} // end reArr_tableStates



void qdotQuantumNumber::setUp_tableCouples(){

  // declare variables
  int i,j,k,l,m,n;
  int nbCouples, tot_nbCouples;
  int M_tot, S_tot;
  int indexChannel, indexCouple;
  
  // maximum value of total angular momentum (projection) M = m1 + m2
  int M_max = 2*m_max;
  // maximum value of total spin S = s1 + s2
  int S_max = 2;         
  
  // total number of "channels", that is: total number of possibles distinct sets [M,S]
  nbChannels = (S_max+1)*(2*M_max + 1);    

  // declare a channel-array
  tableChannel = new int*[nbChannels];
  for(i=0; i<nbChannels; i++){
    tableChannel[i] = new int[2];
  }  

  // declare a array containing numbers of couples for each channel 
  nbCouplesChannel = new int[nbChannels];
  
  indexChannel = 0;
  tot_nbCouples = 0;
  for(int M=-M_max; M<=M_max; M++){
    for(int S=-1; S<=1; S++){
      nbCouples = 0;
      tableChannel[indexChannel][0] = M;
      tableChannel[indexChannel][1] = S;      
      // count how many couples with total angular momentum and total
      // spin value M and S
      for(i=0; i<nbBasis; i++){
	for(j=i; j<nbBasis; j++){
	  M_tot = tableStates[i][index_m] + tableStates[j][index_m];
	  S_tot = (tableStates[i][index_s] + tableStates[j][index_s])/2; // devided ny 2 since s=+- 1 in program
	  if(M_tot == M && S_tot == S){
	    nbCouples++;
	  }
	}
      }
      
      // store numbers of states within channel 
      nbCouplesChannel[indexChannel] = nbCouples;
      tot_nbCouples += nbCouples;
      // update index variable, i.e. next channel
      indexChannel++;
    }
  }
  
  // initialize table of couples
  tableCouples = new int**[nbChannels];
  for(i=0; i<nbChannels; i++){
    tableCouples[i] = new int*[nbCouplesChannel[i]]; 
  }
  for(i=0; i<nbChannels; i++){
    for(j=0; j<nbCouplesChannel[i]; j++){
      tableCouples[i][j] = new int[2];
    }
  }

  // set up table of couples
  indexChannel = 0;
  for(int M=-M_max; M<=M_max; M++){
    for(int S=-1; S<=1; S++){
      indexCouple = 0;
      nbCouples = nbCouplesChannel[indexChannel];
      for(i=0; i<nbBasis; i++){
	for(j=i; j<nbBasis; j++){
	  M_tot = tableStates[i][index_m] + tableStates[j][index_m];
	  S_tot = (tableStates[i][index_s] + tableStates[j][index_s])/2;
	  if(M_tot == M && S_tot == S){
	    tableCouples[indexChannel][indexCouple][0] = i;
	    tableCouples[indexChannel][indexCouple][1] = j;
	    indexCouple++;
	  }
	}
      }
      indexChannel++;
    }
  }
} // end setUp_coupleStatesTable



void qdotQuantumNumber::print_tableCouples(){

  // variable decleration
  int i,j,k;
  printf("\n");
  printf("*** TABLE OG COUPLES ***\n");
  for(i=0; i<nbChannels; i++){    
    printf("Channel %i - [M=%i,S=%i]   -->  ",i,tableChannel[i][0],tableChannel[i][1]); 
    for(j=0; j<nbCouplesChannel[i]; j++){
      printf("|%i,%i> ",tableCouples[i][j][0],tableCouples[i][j][1]);
    }
    printf("\n");
  }
} // end print_tableCouples



void qdotQuantumNumber::setUp_stateMapChannel(){
  
  // variable decleration
  int i,j,k,l,indexChannel,indexCouple,nbCouples;
  
  // initialize table
  stateMapChannel = new int**[nbBasis];
  for(i=0; i<nbBasis; i++){
    stateMapChannel[i] = new int*[nbBasis];
  }
  for(i=0; i<nbBasis; i++){
    for(j=0; j<nbBasis; j++){
      stateMapChannel[i][j] = new int[2];
    }
  }  
  for(i=0; i<nbBasis; i++){
    for(j=0; j<nbBasis; j++){
      for(indexChannel = 0; indexChannel<nbChannels; indexChannel++){
	nbCouples = nbCouplesChannel[indexChannel];
	for(indexCouple = 0; indexCouple<nbCouples; indexCouple++){
	  k = tableCouples[indexChannel][indexCouple][0];
	  l = tableCouples[indexChannel][indexCouple][1];
	  if(k==i && l==j){
	    stateMapChannel[i][j][0] = indexChannel;
	    stateMapChannel[i][j][1] = indexCouple;
	  }
	}
      }
    }
  }
} // end setUp_coupleMapChannel



void qdotQuantumNumber::print_stateMapChannel(){
  
  // variable decleration
  int i,j;

  printf("\n");
  printf("*** TABLE OF MAPPING INFO STATE --> CHANNEL#,COUPLE#\n");
  for(i=0; i<nbBasis; i++){
    printf("|%i>   ",i);
  }
  for(i=0; i<nbBasis; i++){
    printf("\n|%i>   ",i);
    for(j=0; j<nbBasis; j++){
      printf("[%i,%i]   ",stateMapChannel[i][j][0], stateMapChannel[i][j][1]);
    }
  }
  printf("\n");
} // end print_stateMapChannel



void qdotQuantumNumber::getChannels(int alpha, int beta, int gamma, int delta, int& chanBra, int& chanKet){
  if(alpha<=beta && gamma<=delta){
    chanBra = stateMapChannel[alpha][beta][0];
    chanKet = stateMapChannel[gamma][delta][0];
  }
  else if(alpha<=beta && gamma>delta){
    chanBra = stateMapChannel[alpha][beta][0];
    chanKet = stateMapChannel[delta][gamma][0];
  }
  else if(alpha>beta && gamma<=delta){
    chanBra = stateMapChannel[beta][alpha][0];
    chanKet = stateMapChannel[gamma][delta][0];
  } 
  else {
    chanBra = stateMapChannel[beta][alpha][0];
    chanKet = stateMapChannel[delta][gamma][0];
  } 
}  



void qdotQuantumNumber::getCouples(int alpha, int beta, int gamma, int delta, int& coupBra, int& coupKet){
  if(alpha<=beta && gamma<=delta){
    coupBra = stateMapChannel[alpha][beta][1];
    coupKet = stateMapChannel[gamma][delta][1];
  }
  else if(alpha<=beta && gamma>delta){
    coupBra = stateMapChannel[alpha][beta][1];
    coupKet = stateMapChannel[delta][gamma][1];
  }
  else if(alpha>beta && gamma<=delta){
    coupBra = stateMapChannel[beta][alpha][1];
    coupKet = stateMapChannel[gamma][delta][1];
  } 
  else {
    coupBra = stateMapChannel[beta][alpha][1];
    coupKet = stateMapChannel[delta][gamma][1];
  }   
}
