#include <iostream>
#include "singleParticleOrbitals.h"

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


  singleParticleOrbitals Basis(dim, Nlevel);
   
  printf("you choose a basis with dimension %d and the max energy level E=%d, so %d possible states \n", Basis.Dimension(), Basis.maxEnergylevel(), Basis.update_totalNbStates());

  int index_state;
  int n, l, ml, ms;
  
#if 0
  int j=0;
  while(j<1)
    {
      j++;
      // Prompt the user to enter the parameter one by one
      cout << "Enter the index of the state for which you want to know the quantum numbers: \n";
      cin >> index_state;
      
      n =Basis.readQuantumNumberOfState(index_state,0);
      // ml=Basis.readQuantumNumberOfState(index_state,1);
      //ms=Basis.readQuantumNumberOfState(index_state,2);
      ml=Basis.readQuantumNumberOfState(index_state,Basis.readIndex_ml());
      ms=Basis.readQuantumNumberOfState(index_state,Basis.readIndex_ms());
      if(dim==2)
	printf("|%d>={%d,%d,%d}\n",index_state,n,ml,ms);
      else
	{
	  l=Basis.readQuantumNumberOfState(index_state,1);
	printf("|%d>={%d,%d,%d,%d}\n",index_state,n,l,ml,ms);
	}
    }
  
  int * nbCouplePerCh= Basis.read_NbCouplesPerChannel();
  printf("in channel %d, there is %d couples\n",5,nbCouplePerCh[5]);
#endif


#if 1
  printf("\nTest of the function: read_NbCouplesPerChannel() \n");
    
  int * nbCouplePerCh= Basis.read_NbCouplesPerChannel();
  int NbCh =Basis.readNbChannels();
  for(int i=0;i<NbCh;i++)
    printf("in channel %d, there is %d couples\n",i,nbCouplePerCh[i]);

  printf("\nTest of the function: read_channelOfCouple(s1,s2) \n");

  while(true)
    {
      int s1,s2;
      
      // Prompt the user to enter the parameter one by one
      cout << "s1:";  cin >> s1;
      cout << "s2: ";      cin >> s2;
      cout << "read_channelOfCouple(s1,s2), channel=" << Basis.read_channelOfCouple(s1,s2) << endl << endl; 
    }
#endif


#if 0
  int index_couple,index_ch;
  j=0;
  while(j<1)
    {
      j++;
      // Prompt the user to enter the parameter one by one

      cout << "Enter a channel and the index of a couple inside this channel: \n";
      cout << "Enter the channel \n";
      cin >> index_ch;
      cout << "Enter the index of a couple to read out the corresponding single states \n";
      cin >> index_couple;
      for(int k=0;k<2;k++)
	{
	  index_state = Basis.mappingCoupleToSingleStates(index_ch,index_couple,k);
	  n =Basis.readQuantumNumberOfState(index_state,0);
	  // ml=Basis.readQuantumNumberOfState(index_state,1);
	  //ms=Basis.readQuantumNumberOfState(index_state,2);
	  ml=Basis.readQuantumNumberOfState(index_state,Basis.readIndex_ml());
	  ms=Basis.readQuantumNumberOfState(index_state,Basis.readIndex_ms());
	  if(dim==2)
	    printf("|%d>={%d,%d,%d}\n",index_state,n,ml,ms);
	  else
	    {
	      l=Basis.readQuantumNumberOfState(index_state,1);
	      printf("|%d>={%d,%d,%d,%d}\n",index_state,n,l,ml,ms);
	    }
	}
    }
#endif

#if 1
  while(true)
    {
 // Prompt the user to enter the parameter one by one
  int s1,s2;
  cout << "Enter index of state 1 \n";
  cin >> s1;
  cout << "Enter index of state 2 \n";
  cin >> s2;
  int indexCouple=Basis.mappingSingleStatesToCouple(s1,s2);
  printf("couple index of (%d,%d) is %d\n\n",s1,s2,indexCouple);
    }
#endif



  return 0;
}
