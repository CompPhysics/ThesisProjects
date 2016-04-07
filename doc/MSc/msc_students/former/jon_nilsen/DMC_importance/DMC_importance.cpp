#include "DMC_importance.hpp"

using namespace std;

// copy walker k to walker l, make space if no space
void DMC::copy_walker(int k, int l){

  int ni = walkers[0].getNoOfParticles();
  int nj = walkers[0].getDimension();

  static int walker_array_length = no_of_walkers;

  // if new walker needs more space
  if(l==walker_array_length)
    walker_array_length = make_room(ni,nj,walker_array_length);

  // copy X(k) into X(l)
  walkers[l]=walkers[k];

}

int DMC::make_room(int ni, int nj, int nk){
  Walker *w_copy = new Walker[nk*2](ni);
  for(int k=0; k!=nk; k++){
    w_copy[k].setParams(params);
  }
  for(int i=0; i!=nk; i++)
    w_copy[i]=walkers[i];
  
  // note doubling nk so we don't have to this every time..
  delete[] walkers;
  walkers = new Walker[nk*2](ni);
  for(int k=0; k!=nk*2; k++){
    walkers[k].setParams(params);
  }
  
  for(int i=0; i!=ni; i++)
    for(int j=0; j!=nj; j++)
      for(int k=0; k!=nk; k++)
	walkers[k]=w_copy[k];
  
  //static int her=0;
  //std::cerr << "her?" << ++her << std::endl;
  
  delete[] w_copy;
  return nk*2;
}

void DMC:: initialize(string filename){

  double random=0;

  ifstream ifile;

  int no_of_params;

  ifile.open(filename.c_str());

  ifile >> no_of_walkers;
  ifile.ignore(100, '\n'); // to ignore the comments in the infile
  ifile >> desired_walkers;
  ifile.ignore(100, '\n');
  ifile >> max_walkers;
  ifile.ignore(100, '\n');
  ifile >> particles;
  ifile.ignore(100, '\n');
  ifile >> dimensions;
  ifile.ignore(100, '\n');
  ifile >> steps;
  ifile.ignore(100, '\n');
  ifile >> termalization;
  ifile.ignore(100, '\n');
  ifile >> step_length;
  ifile.ignore(100, '\n');
  ifile >> idum;
  ifile.ignore(100, '\n');
  ifile >> tau;
  ifile.ignore(100, '\n');
  ifile >> D;
  ifile.ignore(100, '\n');
  ifile >> e_trial;
  ifile.ignore(100, '\n');
  ifile >> no_of_params;
  ifile.ignore(100, '\n');
  params.redim(no_of_params);
  for(int i=0; i!=no_of_params; i++){
    ifile >> params(i);
    ifile.ignore(100, '\n');
  } 

  ifile.close();

  std_dev = 2*D*tau;

}


void DMC:: output(string filename){

  ofstream ofile;

  ofile.open(filename.c_str());
  QickArray pos(no_of_walkers);

  double max_pos=0;
  for(int i=0; i!=no_of_walkers; i++){
    for(int j=0; j!=dimensions; j++)
      for(int k=0; k!=particles; k++){
	pos(i)+=sqrt(sqr(walkers[i].getParticlePosition(k,j)));
      }
    if(pos(i)>max_pos)
      max_pos=pos(i);
  }
  int res = 1000;
  QickArray ploter(res+1);

  // set up histogram, i.e. wf distribution
  for(int i=0; i!=no_of_walkers; i++)
    ploter((int)(pos(i)/max_pos*res))++;
  ploter /= (double)no_of_walkers;

  // and dump it to file
  for(int i=0; i!=res+1; i++)
    ofile << (double)i/res << " " 
	  << ploter(i)/(i*sqr(i)-(i-1)*sqr(i-1))*no_of_walkers << endl;

  ofile.close();
}

void DMC::oneMonteCarloStep(Random &ran, int k){

  double e_local_x, e_local_y, tmp_e_local_y, wf_x, wf_y, branching, 
    diffusion;
  
  // arrays for quantum force
  QickArray fq_x(particles, dimensions), fq_y(particles, dimensions);
  
  // initial info for this walker
  e_local_x = walkers[k].getLocalEnergy();
  wf_x = walkers[k].getWaveFunction();
  fq_x = walkers[k].getQuantumForce();
  
  // move each electron
  double pos;
  for(int i=0; i!=particles; i++)
    for(int j=0; j!=dimensions; j++){
      pos=walkers[k].getParticlePosition(i,j);
      walkers[k].setParticlePosition(i,j,pos+D*tau*fq_x(i,j)+
				     ran.gran(1,0)*sqrt(tau));
    }
  
  // new info for this walker
  wf_y = walkers[k].getWaveFunction();
  
  diffusion = 0;
  double pos_x, pos_y;
  
  double w = walkers[k].getGreensFunction(D, tau)/
    walkers[k].getNewGreensFunction(D, tau);
  
  trials++;
  // metropolis:
  if(sqr(wf_y/wf_x)*w>ran.ran1()){
    for(int i=0; i!=particles; i++)
      walkers[k].updateParticlePosition(i);
    accept++;
  }
  else{
    for(int i=0; i!=particles; i++)
      walkers[k].resetParticlePosition(i);
    return; // skip branching if not accepted
  }
  
  e_local_y = walkers[k].getLocalEnergy();
  branching = -(.5*(e_local_x+e_local_y)-e_trial)*tau;
  branching = exp(branching);
  
  // random integer with average value equal to the branching
  int MB = int(branching);
  if(branching-MB > ran.ran1())
    ++MB;
  // add MB-1 copies at end of list
  // if MB=0, mark this config as dead
  for(int n=0; n<MB-1 && no_of_walkers < max_walkers; n++){
    // copy this walker to the end of the walker list
    copy_walker(k,no_of_walkers);
    no_of_walkers++;
  }

  if(MB == 0){
    walkers[k].killWalker();
  }

}

void DMC::oneTimeStep(Random &ran){

    int M = no_of_walkers;
    
    // let the walkers march
    for(int k=0; k!=M; k++)
      oneMonteCarloStep(ran, k);
    
    // adjust trial energy
    e_trial += log(desired_walkers/double(no_of_walkers))/10;
    
    static int count=0;
    //cerr << ++count << " " << no_of_walkers << " ";
    for(int i=0; i!=M; i++){
      // if dead walker, put in last walker
      if(walkers[i].isDead()){
	no_of_walkers--;
	copy_walker(no_of_walkers,i);
      }
    }
    //cerr << no_of_walkers << " " << e_trial << endl;
    
    energy  += e_trial;
    energy2 += sqr(e_trial);

}

void DMC::diffMC(){

  // object containing various random generators
  Random ran(idum);

  for(int i=0; i!=no_of_walkers; i++){
    //walkers[i].resetWalker(particles, dimensions);
    walkers[i].setParams(params);
  }

  // set initial position
  for(int i=0; i!=particles; i++)
    for(int j=0; j!=dimensions; j++)
      for(int k=0; k!=no_of_walkers; k++){
	walkers[k].setParticlePosition(i,j,
				       step_length*(ran.ran1()-.5)/params(0));
				       //ran.gran(1,0));
	
	walkers[k].updateParticlePosition(i);
      }
  
  accept=trials=0;
  int adjustInterval = int(0.1 * termalization) + 1;
      
  eav=0, energy=0, energy2=0;
  
  for(int i_step=0; i_step!=termalization; i_step++){
    oneTimeStep(ran);
    if ( (i_step+1) % adjustInterval == 0 ) {
      tau *= accept / (0.9 * trials);
      cerr << accept << " " << trials << endl;
      trials = accept = 0;
    }
  } // end for(i_step)
  cerr << "adjustet time step: " << tau << endl;

  eav=0, energy=0, energy2=0;
  
  for(int i_step=0; i_step!=steps; i_step++)
    oneTimeStep(ran);

  energy  /= steps;
  energy2 /= steps;
  energy2 -= sqr(energy);
  
  cerr << "energy= " << energy << " ±" 
       << sqrt(energy2/steps) << endl << "sigma= " 
       << energy2 << endl;
} // end diffMC()
