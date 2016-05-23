#include "walker.hpp"

using namespace std;

Walker:: Walker(int no_of_particles_){
  dim=3;
  no_of_particles = no_of_particles_;
  particles = new Particle[no_of_particles](dim);
  dead = false;
  old_eq_new = true;
  local_e_ptr = new LocalEnergy();
  if(no_of_particles==1)
    wf_ptr = new Wave1s();
  else if(no_of_particles==2)
    wf_ptr = new Helium();
  else
    cerr << "no wave function implemented for " 
	 << no_of_particles << " particles.";
  q_force_ptr = new QForce();
  greens_ptr = new Greens();
  hamilton_ptr = new Hamilton();
  potential_ptr = new Potential();
}

Walker& Walker::operator=(const Walker& clone_me){
  
  old_eq_new = clone_me.old_eq_new;

  // cloning particle positions
  for(int i=0; i!=no_of_particles; i++){
    for(int j=0; j!=dim; j++)
      particles[i].
	setPosition(j,clone_me.particles[i].getPosition(j));
    if(old_eq_new)
      particles[i].updatePosition();
  }

  params = clone_me.params;
  dead = clone_me.dead;
  dim = clone_me.dim;
  no_of_particles = clone_me.no_of_particles;
  particles_moved = clone_me.particles_moved;
  
}
