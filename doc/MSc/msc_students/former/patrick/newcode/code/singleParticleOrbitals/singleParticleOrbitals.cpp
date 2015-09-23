#ifndef singleParticleOrbitals_CPP
#define singleParticleOrbitals_CPP

// file        : singleParticleOrbitals.cpp
// description : implentation of class singleParticleOrbitals

#include <cmath>
#include <cstdlib>
#include "singleParticleOrbitals.h"
#include <algorithm>

using namespace std;

// Contructors
singleParticleOrbitals:: singleParticleOrbitals(int _dim, int _Nmax):dim(_dim), Nmax(_Nmax)
{
#if DEBUGCONSTRUCTOR
  printf("\nConstruction starts: singleParticleOrbitals\n");
#endif
  if( !( _dim==2 || _dim==3) )
    {
      printf("ERROR MESSAGE: specify a dimension of 2 or 3 \n");
      exit(1);
    }
  index_n=0;
  if(dim == 2) // in 2 DIMENSIONS
    {
      index_ml=1;
      index_ms=2;
    }
  else if(dim == 3) // in 3 DIMENSIONS
    {
      index_l=1;
      index_ml=2;
      index_ms=3;
    }
  nbBlocks = 0;
  update_totalNbStates();
  update_Tables();

#if DEBUGCONSTRUCTOR
  printf("Construction finished: singleParticleOrbitals\n");
#endif
}

// singleParticleOrbitals Basis = S; (assigment operator, S and Basis are singleParticleOrbitals objects)
singleParticleOrbitals& singleParticleOrbitals:: operator= (const singleParticleOrbitals& S)
{
  dim=S.Dimension();
  index_n=0;
  index_l=S.readIndex_l(); 
  index_ml=S.readIndex_ml();
  index_ms=S.readIndex_ms();    
  NbChannels=S.readNbChannels();   
  TotNbCouples=S.readTotNbCouples();    
  Nmax=S.maxEnergylevel();        
  totalNbStates=S.readTotalNbStates();
  
  deallocate();
  // re-allocate memory
  TableOfStates = new int*[totalNbStates];
  for(int s=0;s<totalNbStates;s++)
    TableOfStates[s]= new int[dim+1];
  lengthTabInteractionStates = new int[NbChannels];
  TabInteractionStates = new int **[NbChannels];
   for(int ch=0;ch<NbChannels;ch++)
     {
       TabInteractionStates[ch]=new int*[lengthTabInteractionStates[ch]];
       for(int couple=0;couple<lengthTabInteractionStates[ch];couple++)
	 TabInteractionStates[ch][couple]=new int[2];
     }
   
   // assign the value of Basis into the tables of S  
   for(int s=0;s<totalNbStates;s++)
     for(int qn=0;qn<=dim;qn++)
       TableOfStates[s][qn]=S.TableOfStates[s][qn];
   for(int ch=0;ch<NbChannels;ch++)
     {
       lengthTabInteractionStates[ch]=S.lengthTabInteractionStates[ch];
       for(int couple=0;couple<lengthTabInteractionStates[ch];couple++)
	 for(int state=0;state<2;state++)
	   TabInteractionStates[ch][couple][state]=S.TabInteractionStates[ch][couple][state];
     }
   return *this;
}


// Free dynamic memory
void singleParticleOrbitals:: deallocate ()
{
  delete [] TableOfStates;
  delete [] TabInteractionStates;
  delete [] lengthTabInteractionStates;
  delete [] sizeBlock;
}


// Destructor
singleParticleOrbitals:: ~singleParticleOrbitals ()
{
  deallocate();
}


// read the Dimension, "energy cut" (=max. Energy level) of the basis set, and the index of each quantum number which depend on the dimension of the system
int singleParticleOrbitals:: Dimension () const { return dim; } // a copy of the dimension is returned
int singleParticleOrbitals:: maxEnergylevel () const { return Nmax; } // a copy of the max energy level is returned
int singleParticleOrbitals:: readIndex_l  () const {return (*this).index_l; }
int singleParticleOrbitals:: readIndex_ml () const {return (*this).index_ml; }
int singleParticleOrbitals:: readIndex_ms () const {return (*this).index_ms; }


// read the total number of states of the system
int singleParticleOrbitals:: readTotalNbStates () const { return totalNbStates; } // a copy of the total nb. of states is returned

// read the quantum number of a given state
int singleParticleOrbitals:: readQuantumNumberOfState(const int indexState, const int indexNumber) const
{
#if CHECK2
  if(indexState>totalNbStates-1 || indexNumber>dim)
    if(indexState>totalNbStates-1)
      printf("-----you asked for a state index which does not exist\n");
    else
      printf("-----you asked for a index of quantum number which does not exist\n");
#endif
    return (*this).TableOfStates[indexState][indexNumber];
}

// read the number of channels {M,S} in the system
int singleParticleOrbitals:: readNbChannels () const { return NbChannels; } // a copy of the nb. of channels is returned

// read the total number of couples over all {M,S} in the system
int singleParticleOrbitals:: readTotNbCouples () const { return TotNbCouples; } // a copy of the total nb. of states couples is returned

// total number of states over all energy levels up to Nmax, including spin degeneracy
int singleParticleOrbitals:: update_totalNbStates ()
{
  if(dim == 2) // in 2 DIMENSIONS
    totalNbStates = (Nmax+1)*(Nmax+2); // spin degenracy included
  else if(dim == 3) // in 3 DIMENSIONS
    totalNbStates = (Nmax+1)*(Nmax+2)*(Nmax+3)/3; // spin degenracy included
  return totalNbStates;
}


// Initialize and fill the table of states
void singleParticleOrbitals:: update_Tables ()
{
  int Nlevel=(*this).Nmax;
  //  int index_ml;
  // int index_ms;
  //printf("index_ml=%d\n",index_ml);

  update_totalNbStates();
  int totNbStates=(*this).totalNbStates;
#if COMMENT
  printf("---total number of states:%d\n",totNbStates);
#endif

  // Defining/Filling the brute force table of states
  if( (*this).dim == 2 ) // in 2 DIMENSIONS
    {
      // index_ml=1;
      // index_ms=2;
      (*this).TableOfStates = new int*[totNbStates];
      for(int s=0; s<totNbStates; s++)
	{
	  (*this).TableOfStates[s]=new int[dim+1];// for the 3 quantum numbers n,ml, ms
	}

      // Filling the table
      int state=0; // index on states in the table of states
      
      for(int N=0; N<Nlevel+1; N++)
	{
	  for(int n=0; n<N/2+1; n++)
	    {
	      for(int ml=-N;ml<N+1;ml++)
		{
		  if( (2*n+abs(ml))==N)
		    {
		      for(int spin=0; spin<2; spin++)
			{
			  int ms = spin*2-1;
			  (*this).TableOfStates[state][0]=n;
			  (*this).TableOfStates[state][1]=ml;
			  (*this).TableOfStates[state][2]=ms;
// 			   #if CHECK
// 			  printf("---state #%02d |n=%2d,ml=%2d,ms=%2d>\n",state,n,ml,ms);
// 			  #endif
			  state++;
			}
		    }
		}
	    }	
	}
    }
  else if((*this).dim == 3) // in 3 DIMENSIONS
    {
      //      index_ml=2;
      //index_ms=3;
      (*this).TableOfStates = new int*[totNbStates];
      for(int s=0; s<totNbStates; s++)
	{
	  (*this).TableOfStates[s]=new int[dim+1];// for the 4 quantum numbers n, l, ml, ms
	}

      // Filling the table
      int state=0; // index on states in the table of states
        
      for(int N=0; N<Nlevel+1; N++)
	{
	  for(int n=0; n<N/2+1; n++)
	    {
	      for(int l=0; l<N+1; l++)
		{
		  if( (2*n+l)==N)
		    {
		      for(int ml=-l;ml<l+1;ml++)
			for(int spin=0; spin<2; spin++)
			  {
			    int ms = spin*2-1;
			    (*this).TableOfStates[state][0]=n;
			    (*this).TableOfStates[state][1]=l;
			    (*this).TableOfStates[state][2]=ml;
			    (*this).TableOfStates[state][3]=ms;
// 			    #if CHECK
// 			    printf("---state #%02d |n=%2d,l=%2d,ml=%2d,ms=%2d>\n",state,n,l,ml,ms);
//  			    #endif
			    state++;
			  }
		    }
		}
	    }	
	}
    }

  // sort the Table Of States according to the "ml"  quantum number
#if CHECK
  printf("\nSorting the states by blocks of equal angular momentum and projection\n");
#endif

  if(dim == 2) // sorting according to increasing "ml" values in 2D
    {
      sort_TableOfStates(index_ml);
      update_infoBlocks(index_ml);
    }
  else if(dim == 3) // sorting according to increasing "l", then "ml" values in 3D
    {
      sort_TableOfStates(index_l,index_ml);
      update_infoBlocks(index_l,index_ml);
    }

#if CHECK
  // print the table to screen
  print_TableOfStates();
#endif

  // find the biggest total angular momentum considering 2 particles
  Mmax=0;
  for(int s=0; s<totNbStates ; s++)
    {
      int maximumM=(*this).TableOfStates[s][index_ml];
      if( 2*maximumM>Mmax )
	{
	  Mmax = 2*maximumM;
	}
    }
#if CHECK
  printf("\n---Classification of couple of states with same {M,S}\nMmax=%d\n",Mmax);
#endif

  // New table of coupled states with classification by channel ch.={M,S}, i.e. by couple {total angular momentum M=m1+m2, total spin S=s1+s2}
  // Mapping channel <-> {M,S}: channel = 2*M+S 
  // = stores the couples using the index of TableOfStates

  // definition of the new table of states
  int NbChannel=2*(Mmax+1);
  (*this).NbChannels= NbChannel;
  (*this).TabInteractionStates= new int**[NbChannel];
  lengthTabInteractionStates= new int[NbChannel];

  // for a given M and S, fill in the table with the coupled states corresponding
  int totalM=0;
  int totalS=0;
  int nbCouples;
  (*this).TotNbCouples=0;
  int numChannel;
  
  // Count the number of couples with the same {M,S}
  // fill in the table lengthTabInteractionStates
  for(int M=0; M<Mmax+1; M++)
    for(int S=0; S<2; S++)
      {
	numChannel=2*M+S;
	nbCouples=0;
	for(int s1=0; s1 <totNbStates-1; s1++)
	  for(int s2=s1+1; s2<totNbStates; s2++)
	  {
	    totalM = abs((*this).TableOfStates[s1][index_ml]+(*this).TableOfStates[s2][index_ml]);
	    totalS = abs((*this).TableOfStates[s1][index_ms]+(*this).TableOfStates[s2][index_ms]);

	    if(totalM==M && totalS==2*S)
	      nbCouples++;
	  }
	(*this).lengthTabInteractionStates[numChannel]=nbCouples;
	TotNbCouples+=nbCouples;
      }
  
  // initialize the size of TabInteractionStates
   for(int ch=0; ch<NbChannel; ch++)
     {
       nbCouples=(*this).lengthTabInteractionStates[ch];
       (*this).TabInteractionStates[ch]= new int*[nbCouples];
       for(int indexCouple=0; indexCouple<nbCouples; indexCouple++)
	 (*this).TabInteractionStates[ch][indexCouple]= new int[2];
     }
   
  int indexCouple;
  // fill in the table TabInteractionStates now
  for(int M=0; M<Mmax+1; M++)
    for(int S=0; S<2; S++)
      {
	indexCouple=0;
	numChannel=2*M+S;
	nbCouples=(*this).lengthTabInteractionStates[numChannel];
	while(indexCouple<nbCouples)
	  {
	    for(int s1=0; s1 <totNbStates-1; s1++)
	      for(int s2=s1+1; s2<totNbStates; s2++)
		{
		  totalM = abs((*this).TableOfStates[s1][index_ml]+(*this).TableOfStates[s2][index_ml]);
		  totalS = abs((*this).TableOfStates[s1][index_ms]+(*this).TableOfStates[s2][index_ms]);
		  
		  if(totalM==M && totalS==2*S)
		    {
		      (*this).TabInteractionStates[numChannel][indexCouple][0]=s1;
		      (*this).TabInteractionStates[numChannel][indexCouple][1]=s2;
		      indexCouple++;
		    }
		}
	  }
      } 
 
  // print resulting table on screen
  for(int M=0; M<Mmax+1; M++)
    for(int S=0; S<2; S++)
      {
#if CHECK
	printf("---channel #%d {M=%d S=%d} : ",2*M+S,M,S);
#endif
	int sizebuf=lengthTabInteractionStates[2*M+S];
	 #if CHECK
	for(int i=0; i<sizebuf;i++)
	  printf("(%d,%d) ",TabInteractionStates[2*M+S][i][0],TabInteractionStates[2*M+S][i][1]);
	printf("\n");
	 #endif
      }
  

   #if COMMENT
  printf("---Table of states generated \n");
  #endif
}// end of update_TableOfStates()


// given ch=2M+S, and the index of a couple in this channel, read each Single particle states out of the table.
int singleParticleOrbitals:: mappingCoupleToSingleStates (const int channel, const int indexCouple, const int leftRightState) const
{
  int indexSingleState;
  indexSingleState = (*this).TabInteractionStates[channel][indexCouple][leftRightState]; // leftRightState = {0 or 1}
  return indexSingleState;
}

// read the array lengthTabInteractionStates, to have the number of state couples for each channel.
 int * singleParticleOrbitals:: read_NbCouplesPerChannel () const
 {
   return lengthTabInteractionStates;
 }

// given 2 single particle states, compute the corresponding channel index ch={2M+S}
int singleParticleOrbitals:: read_channelOfCouple (const int s1, const int s2) const
{
  int channel;
  // Couple made of identical single particle states requires special treatment, so need to be identified
if(s1==s2)
    {
      channel = -1;// set the index of such couple to be -1
      //printf("\nERROR: a couple of single particle states build with 2 identical set \n of quantum numbers (including spin) should not be allowed\n");
      //exit(1);
    }
  else
    {
      int totalM = abs(TableOfStates[s1][index_ml]+TableOfStates[s2][index_ml]);
      int totalS = abs(TableOfStates[s1][index_ms]+TableOfStates[s2][index_ms])/2;
      channel = 2*totalM+totalS;
    }
 return channel;
}

// given 2 single particle states, find the corresponding index for the couple.
int singleParticleOrbitals:: mappingSingleStatesToCouple (const int s1, const int s2) const
{
  int channel=read_channelOfCouple(s1,s2);
  int lowerState=min(s1,s2);
  int higherState=max(s1,s2);
  int couple=0;
  while( !(TabInteractionStates[channel][couple][0]==lowerState && TabInteractionStates[channel][couple][1]==higherState) )
    couple++;
    return couple;
}

 // sort the states in the TableOfStates with increasing values of the chosen quantum number
void singleParticleOrbitals:: sort_TableOfStates (const int indexQuantumNumber)
{
  // find smallest/biggest value of the quantum number
  int minQN, maxQN;
  int * listValues= new int[totalNbStates];
  for(int s=0; s<totalNbStates; s++)
    listValues[s]=TableOfStates[s][indexQuantumNumber];
  minQN=*min_element(listValues,listValues+totalNbStates);
  maxQN=*max_element(listValues,listValues+totalNbStates);
  
  // Save a copy of TableOfStates
  int **Old_TableOfStates = new int*[totalNbStates];
  for(int s=0;s<totalNbStates;s++)
    Old_TableOfStates[s]= new int[dim+1];
  for(int s=0;s<totalNbStates;s++)
     for(int qn=0;qn<dim+1;qn++)
       Old_TableOfStates[s][qn]= TableOfStates[s][qn];
  
  // re-arrange the TableOfStates with increasing values of the selected quantum number
  int var_qn=minQN;
  int indexState=0;

  while(indexState <= totalNbStates-1)
    {
      for(int s=0;s<totalNbStates;s++) // loop over the states of one block
	{
	  if(Old_TableOfStates[s][indexQuantumNumber]==var_qn)
	    {
	      // re-assign the states with its corresponding quantum numbers
	      for(int qn=0;qn<dim+1;qn++)
		{
		  TableOfStates[indexState][qn]=Old_TableOfStates[s][qn];
		  //cout << "TableOfStates[" << indexState << "][" << qn<< "]="<<Old_TableOfStates[s][qn]<<endl;
		}
	      indexState++;
	    }
	}// end loop over one block with the same QN
      var_qn++;
    }// end of while()
#if CHECK2
  // Print for test
  printf("\n print test for re-ordering TableOfStates for increasing ml values:\n");
    printf("\nAfter Sorting first QN\n");
  print_TableOfStates();
#endif
  delete [] listValues;
  delete [] Old_TableOfStates; 
}// end of sort_TableOfStates (qn1)

 
// find the number of blocks with the same QN and the nb of states inside each block
void singleParticleOrbitals:: update_infoBlocks (const int indexQuantumNumber)
{
  // find the number of blocks and the number of states per block
  int ValQN= TableOfStates[0][indexQuantumNumber];
  int * MaxSizeBlock = new int[totalNbStates];
  

  for(int s=0; s<totalNbStates; s++) // zeroing MaxSizeBlock
    MaxSizeBlock[s]=0;

  nbBlocks=1;MaxSizeBlock[0]=1;// the first state belongs to the first block and count for 1

  for(int s=1; s<totalNbStates; s++)
    {
      if( TableOfStates[s][indexQuantumNumber] != ValQN)
	{
	  //printf("state #%d not in the same block\n",s);
	  nbBlocks++;
	  MaxSizeBlock[nbBlocks-1]=MaxSizeBlock[nbBlocks-1]+1;
	  ValQN = TableOfStates[s][indexQuantumNumber];

	}
      else
	{
	  //printf("state #%d, in the same block\n",s);
	  MaxSizeBlock[nbBlocks-1]=MaxSizeBlock[nbBlocks-1]+1;
	}
      //printf("state #%d not in the same block, nbBlock=%d, MaxSizeBl[%d]=%d\n",s,nbBlocks,nbBlocks-1,MaxSizeBlock[nbBlocks-1]);
    }

  sizeBlock = new int[nbBlocks];
  for(int block=0; block<nbBlocks; block++)
    sizeBlock[block]=MaxSizeBlock[block];
  
  delete [] MaxSizeBlock;
#if CHECK
  // Print for test
  printf("\nNb blocks first QN=%2d\n",nbBlocks);printf("sizeBlock=[");
  for(int block=0; block<nbBlocks; block++)
    printf("%2d, ",sizeBlock[block]);printf("]\n");
#endif

}//end of update_infoBlocks(QN)



// sort the states in the TableOfStates with increasing values of the chosen quantum number 1 first,
// then for a given value of the qn1, sort with increasing values of the qn2
void singleParticleOrbitals:: sort_TableOfStates (const int indexQn1, const int indexQn2)
{
  // sort according to the first quantum number
  sort_TableOfStates(indexQn1);
  
  // count the blocks and states inside it according to the first QN
  update_infoBlocks(indexQn1);

  // Save a copy of TableOfStates
  int **Old_TableOfStates = new int*[totalNbStates];
  for(int s=0;s<totalNbStates;s++)
    Old_TableOfStates[s]= new int[dim+1];
  for(int s=0;s<totalNbStates;s++)
    for(int qn=0;qn<dim+1;qn++)
      Old_TableOfStates[s][qn]= TableOfStates[s][qn];
  
  // within each block
  int start=0;// index of the state at the beginning of each block at the level of QN1
  for(int bl=0; bl<nbBlocks; bl++)
    {
      // find the smallest value of the quantum number 2 within the block
      int minQN2;
 
      int * listValues2= new int[sizeBlock[bl]];// array with just the values of QN2 for the states in this block
      for(int s=0; s<sizeBlock[bl]; s++)
	listValues2[s]=TableOfStates[start+s][indexQn2];
      minQN2=*min_element(listValues2,listValues2+sizeBlock[bl]);  

      // re-arrange states with increasing QN2
      int var_qn=minQN2;
      int indexState=0;
      while(indexState <= sizeBlock[bl]-1)
	{
	  for(int s=start;s<start+sizeBlock[bl];s++) // loop over the states of one block
	    {
	      if(Old_TableOfStates[s][indexQn2]==var_qn)
		{
		  // re-assign the states with its corresponding quantum numbers
		  for(int qn=0;qn<dim+1;qn++)
		    {
		      TableOfStates[start+indexState][qn]=Old_TableOfStates[s][qn];
		      //cout << "TableOfStates[" << indexState << "][" << qn<< "]="<<Old_TableOfStates[s][qn]<<endl;
		    }
		  indexState++;
		}
	    }// end loop over one block with the same QN
	  var_qn++;
	}// end of while()

      start += sizeBlock[bl];// update the index of the first state in the block
      delete [] listValues2;  
    }// end of loop over the blocks

#if CHECK2
  // Print for test
  printf("\n print test for re-ordering TableOfStates for increasing l and ml values:\n");
    printf("\nAfter Sorting second QN\n");
  print_TableOfStates();
#endif
 
  delete [] Old_TableOfStates; 
}// end of sort_TableOfStates (qn1, qn2)




// find the number of blocks with the same QN1 and QN2 and the size of each block (the nb of states) 
void singleParticleOrbitals:: update_infoBlocks (const int indexQuantumNumber1,const int indexQuantumNumber2)
{
//   // count the blocks and states inside it according to the first QN
//   update_infoBlocks(indexQuantumNumber1);
  
  int nbBlockQN1 = nbBlocks;
  nbBlocks=0; 
  //  printf("nbBlocks=%d",nbBlocks);cin.get();
  int * sizeBlockQN1 = new int [nbBlockQN1];
  for(int i=0;i<nbBlockQN1;i++)
    sizeBlockQN1[i]=sizeBlock[i];
  int * sizeBlockQN2 = new int [totalNbStates];
  // for each block at the level of QN1, find the number of subblocks at the level of QN2
  int startBlock=0;// index of the states at the beginning of each block at the level of QN1
  for(int bl1=0;bl1<nbBlockQN1;bl1++)
    {
      nbBlocks++;// at the beginning of each new block at level QN1, this is counted as a new block/subBlock
      //      printf("\nnbBlocks=%d",nbBlocks);cin.get();
      //      printf("hello2, blockQN1=%d",bl1);cin.get();
      int ValQN= TableOfStates[startBlock][indexQuantumNumber2];
      int * MaxSizeSubBlock = new int[sizeBlockQN1[bl1]];// the maximum nb of subblocks may be the size of the block if the couples of QN is always unique
      for(int s=0; s<sizeBlockQN1[bl1]; s++) // zeroing MaxSizeSubBlock
	MaxSizeSubBlock[s]=0;
      int NbSubBlock=1;MaxSizeSubBlock[0]=1;// the first state belongs to the first subblock and count for 1 state
      
      // increase the nb of states in the block as long as they have the same QNs
      for(int s=1; s<sizeBlockQN1[bl1]; s++)
	{
	  //	  printf("hello3, state %d in the block",s);cin.get();
	  if( TableOfStates[s+startBlock][indexQuantumNumber2] != ValQN)
	    {

	      // not in the same block
	      nbBlocks++;
	      // printf("nbBlocks=%d",nbBlocks);cin.get();
	      NbSubBlock++;
	      MaxSizeSubBlock[NbSubBlock-1]=MaxSizeSubBlock[NbSubBlock-1]+1;
	      ValQN = TableOfStates[s+startBlock][indexQuantumNumber2];
	    }
	  else
	    {
	      // in the same block
	      MaxSizeSubBlock[NbSubBlock-1] += 1;
	    }
	  //printf("state #%d not in the same block, nbBlock=%d, MaxSizeBl[%d]=%d\n",s,nbBlocks,nbBlocks-1,MaxSizeSubBlock[nbBlocks-1]);
	}

      //      printf("hello4, in block[%d] # subBlock=%d",bl1,NbSubBlock);cin.get();
      // store the info about the subblocks in a temp array
      for(int block=0; block<NbSubBlock; block++)
	{
	  int indexBlockTemp=nbBlocks-NbSubBlock+block;
	  //	  printf("hello5, copy nb of block[%d]=%d in temp array",block,MaxSizeSubBlock[block]);cin.get();
	  sizeBlockQN2[indexBlockTemp]=MaxSizeSubBlock[block];
	  //	  printf("sizeBlockQN2[%d]=%d",indexBlockTemp,MaxSizeSubBlock[block]);cin.get();
	}
      // find the number of blocks and the number of states per block
      startBlock += sizeBlockQN1[bl1];
      delete [] MaxSizeSubBlock;
    }// end loop over block at QN1 level
  //  printf("hello6");cin.get();
  //store the info over the total numbers of (sub)blocks and their size in sizeBlock
  sizeBlock = new int[nbBlocks];
  for(int block=0; block<nbBlocks; block++)
    sizeBlock[block]=sizeBlockQN2[block];

#if CHECK
  // Print for test
  printf("\nNb blocks second QN=%2d\n",nbBlocks);printf("sizeBlock=[");
  for(int block=0; block<nbBlocks; block++)
    printf("%2d, ",sizeBlock[block]);printf("]\n");
#endif
  delete [] sizeBlockQN1;
  delete [] sizeBlockQN2;
}//end of update_infoBlocks(QN1,QN2)








// print to screen the list of possible single particle states 
void singleParticleOrbitals:: print_TableOfStates ()
{
  //  printf("\n---total number of states: %d\n",totalNbStates);
  printf("\nTable Of States:\n");
  for(int indexState=0;indexState<totalNbStates;indexState++)
    {
      if(dim==2)
	printf("---state #%02d |n=%2d,ml=%2d,ms=%2d>\n",indexState,readQuantumNumberOfState(indexState,0),readQuantumNumberOfState(indexState,1),readQuantumNumberOfState(indexState,2));
      if(dim==3)
	printf("---state #%02d |n=%2d,l=%2d,ml=%2d,ms=%2d>\n",indexState,readQuantumNumberOfState(indexState,0),readQuantumNumberOfState(indexState,1),readQuantumNumberOfState(indexState,2),readQuantumNumberOfState(indexState,3));
    }
}

// Write to file the content of the sparse Coulomb matrix
void singleParticleOrbitals:: save2file(ostream & mFile)
{
  cout << "Saving tables of states ..."; 
  mFile << endl << endl << "%%%%%%%% Table of states with their quantum numbers, grouped by identical angular momentum %%%%%%%%" << endl;
  mFile << "NbStates= "<< totalNbStates << ";   % Total number of states" <<endl;
  mFile << "states= sparse("<< totalNbStates <<","<< dim+1 <<");" << endl;
  if(dim==2) // in 2D
    {
      for(int indexState=0;indexState<totalNbStates;indexState++)
	for(int qn=0;qn<dim+1;qn++)
	  {
	    mFile << "states(" << indexState+1 << "," << qn+1 << ")=" << readQuantumNumberOfState(indexState,qn)  << "; % ---state #"<< indexState<<" |n=" << readQuantumNumberOfState(indexState,0)<< ",ml="<<readQuantumNumberOfState(indexState,1)<<",ms=" << readQuantumNumberOfState(indexState,2) << ">"<< endl;
	  }
    }
  if(dim==3) // in 3D
    {
      for(int indexState=0;indexState<totalNbStates;indexState++)
	for(int qn=0;qn<dim+1;qn++)
	  mFile << "states(" << indexState+1 << "," << qn+1 << ")=" << readQuantumNumberOfState(indexState,qn) << "; % ---state #"<< indexState<<" |n=" << readQuantumNumberOfState(indexState,0)<< ",l="<<readQuantumNumberOfState(indexState,1)<< ",ml="<<readQuantumNumberOfState(indexState,2)<<",ms=" << readQuantumNumberOfState(indexState,3) << ">"<< endl;
    }
    
    
  mFile << endl << endl << "%%%%%%%% Table of couples of states with the same {M,S} %%%%%%%%" << endl;
  int NbChannels= readNbChannels();
  int * NbCouplesPerCh=read_NbCouplesPerChannel();
  for(int M=0; M<Mmax+1; M++)
    for(int S=0; S<2; S++)
      {    
	int NBCouples=NbCouplesPerCh[2*M+S];// number of couples in the given channel
	int sizebuf=lengthTabInteractionStates[2*M+S];// nb. of couples per channel
	mFile << endl <<"StateCouples_" << setfill ('0') << setw (4) << 2*M+S << "=sparse(2,"<< NBCouples <<");" << endl;
	mFile << "%---channel #"<< 2*M+S<< " {M="<<M<<" S=" <<S<<"}: ";
	for(int i=0; i<sizebuf;i++)
	  mFile << "("<< TabInteractionStates[2*M+S][i][0] <<","<<TabInteractionStates[2*M+S][i][1]<<") ";

	mFile << endl;
	for(int i=0; i<sizebuf;i++)
	  {
	    mFile << "StateCouples_" << setfill ('0') << setw (4) << 2*M+S  << "(1,"<< i+1 <<")= "<<TabInteractionStates[2*M+S][i][0]<< ";";
	    mFile << "StateCouples_" << setfill ('0') << setw (4) << 2*M+S  << "(2,"<< i+1 <<")= "<<TabInteractionStates[2*M+S][i][1]<< ";" << endl;
	  }
      }

  cout << "  finished! :)\n";
}// end of save2file()


#endif // end of singleParticleOrbitals_CPP
