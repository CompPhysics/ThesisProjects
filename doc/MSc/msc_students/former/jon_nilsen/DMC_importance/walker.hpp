#ifndef walker_hpp_IS_INCLUDED
#define walker_hpp_IS_INCLUDED

// class containing information with regards to a walker

#include "QickArray/QickArray.hpp"
#include "particle.hpp"
#include "functions.hpp"

class Walker{

private:

  Particle  *particles;
  QickArray params;
  bool dead, old_eq_new;
  int dim, no_of_particles;
  int particles_moved;

  Func *local_e_ptr, *wf_ptr, *q_force_ptr, *greens_ptr, *hamilton_ptr,
    *potential_ptr;
  
  QickArray& posVec();
  QickArray& newPosVec();

  inline double sqr(double x){ return x*x; }
public:

  Walker(int no_of_particles_);

  ~Walker(){
    delete[] particles;
    delete   local_e_ptr;
    delete   wf_ptr;
    delete   q_force_ptr;
    delete   hamilton_ptr;
    delete   potential_ptr;
  }

  double getParticlePosition(int i, int j); // returns j'th coord position 
                                    // of i'th particle
  double getNewParticlePosition(int i, int j); // returns j'th coord position 
                                    // of i'th particle's new pos
  void setParticlePosition(int i, int j, double x); // set's value of 
                                                    //j'th coord,
                                                    // i'th particle
  void updateParticlePosition(int i); // pos=new_pos
  void resetParticlePosition(int i); // new_pos=pos

  inline void setParams(QickArray& params_){ params=params_; }

  double getLocalEnergy();

  double getWaveFunction();

  QickArray& getQuantumForce();
  
  double getGreensFunction(double D, double tau);
  
  double getNewGreensFunction(double D, double tau);
  
  double getLocalEnergy(int i);

  double getWaveFunction(int i);

  QickArray& getQuantumForce(int i);
  
  inline void killWalker(){ dead=true; }

  inline bool isDead(){ return dead; }

  inline int getNoOfParticles(){ return no_of_particles; }

  inline int getDimension(){ return dim; }

  Walker& operator=(const Walker& clone_me);

};

inline double Walker:: getParticlePosition(int i, int j){
  return particles[i].getPosition(j);
}

inline double Walker:: getNewParticlePosition(int i, int j){
  return particles[i].getNewPosition(j);
}

inline void Walker:: setParticlePosition(int i, int j, double x){
  particles[i].setPosition(j,x);
  old_eq_new = false;
}

inline void Walker:: updateParticlePosition(int i){
  particles[i].updatePosition();
  old_eq_new = true;
}

inline void Walker:: resetParticlePosition(int i){
  particles[i].resetPosition();
  old_eq_new = true;
}

inline double Walker:: getLocalEnergy(){
  return (*local_e_ptr)((*wf_ptr), (*hamilton_ptr), 
			(*potential_ptr), newPosVec(), params);

  /*
  // pade-jastrow local energy
  static QickArray pos(no_of_particles, dim);
  double r1=0, r2=0, r12=0;
  pos = newPosVec();
  for(int j=0; j!=dim; j++){
    r1 += sqr(pos(0,j));
    r2 += sqr(pos(1,j));
    r12 += sqr(pos(0,j)-pos(1,j));
  }
  r1  = sqrt(r1);
  r2  = sqrt(r2);
  r12 = sqrt(r12);

  double dot_prod=0;
  for(int j=0; j!=dim; j++)
    dot_prod += (pos(0,j)-pos(1,j))/r12*(pos(0,j)/r1-pos(1,j)/r2);
  double denom = 1/(1+params(1)*r12);
  double denom2 = sqr(denom);
  double denom3 = denom2*denom;
  double denom4 = sqr(denom2);

  return -4 + params(1)*(denom+denom2+denom3) - denom4/4 + dot_prod*denom2;
  */
} // end getLocalEnergy


inline double Walker:: getWaveFunction(){
  return (*wf_ptr)(newPosVec(), params);
  /*
  static QickArray pos(no_of_particles, dim);
  double r1=0, r2=0, r12=0;
  pos = newPosVec();
  for(int j=0; j!=dim; j++){
    r1 += sqr(pos(0,j));
    r2 += sqr(pos(1,j));
    r12 += sqr(pos(0,j)-pos(1,j));
  }
  r1  = sqrt(r1);
  r2  = sqrt(r2);
  r12 = sqrt(r12);
  return exp(-params(0)*(r1+r2)+.5*r12/(1+params(1)*r12));
  */
} // end getWaveFunction


inline QickArray& Walker:: getQuantumForce(){
  static QickArray qf(no_of_particles, dim);
  return (*q_force_ptr)((*wf_ptr), newPosVec(), qf, params);
  /*
  double r1=0, r2=0, r12=0;
  static QickArray pos(no_of_particles, dim);
  pos = newPosVec();
  for(int j=0; j!=dim; j++){
    r1 += sqr(pos(0,j));
    r2 += sqr(pos(1,j));
    r12 += sqr(pos(0,j)-pos(1,j));
  }
  r1  = sqrt(r1);
  r2  = sqrt(r2);
  r12 = sqrt(r12);

  qf *= 0;

  double denom2 = sqr(1/(1+params(1)*r12));
  
  for(int j=0; j!=dim; j++){
    qf(0,j) += -4*pos(0,j)/r1 + denom2*(pos(0,j)-pos(1,j))/r12;
    qf(1,j) += -4*pos(1,j)/r2 + denom2*(pos(1,j)-pos(0,j))/r12;
  }

  return qf;
  */
} // end getQuantumForce

inline double Walker:: getGreensFunction(double D, double tau){
  return (*greens_ptr)((*wf_ptr), (*q_force_ptr), posVec(), newPosVec(),
		       params, D, tau);
  /*
  static QickArray pos(no_of_particles, dim);
  pos = posVec();
  static QickArray new_pos(no_of_particles, dim);
  new_pos = newPosVec();
  static QickArray qf(no_of_particles, dim);
  qf = getQuantumForce();
  double G=0;
  for(int i=0; i!=no_of_particles; i++)
    for(int j=0; j!=dim; j++)
      G += sqr(pos(i,j)-new_pos(i,j)-D*tau*qf(i,j));
  return exp(-D*G/tau);
  */
} // end getQuantumForce

inline double Walker:: getNewGreensFunction(double D, double tau){
  return (*greens_ptr)((*wf_ptr), (*q_force_ptr), newPosVec(), posVec(),
		       params, D, tau);
  /*
  static QickArray pos(no_of_particles, dim);
  pos = posVec();
  static QickArray new_pos(no_of_particles, dim);
  new_pos = newPosVec();
  static QickArray qf(no_of_particles, dim);
  qf = getQuantumForce();
  double G=0;
  for(int i=0; i!=no_of_particles; i++)
    for(int j=0; j!=dim; j++)
      G += sqr(new_pos(i,j)-pos(i,j)-D*tau*qf(i,j));
  return exp(-D*G/tau);
  */
} // end getQuantumForce

// not implemented yet
inline double Walker:: getLocalEnergy(int i){
  return 0;
} // end getLocalEnergy

// not implemented yet
inline double Walker:: getWaveFunction(int i){
  return 0;
} // end getWaveFunction

// not implemented yet
inline QickArray& Walker:: getQuantumForce(int i){
  static QickArray dummy;
  return dummy;
} // end getQuantumForce


inline QickArray& Walker:: posVec(){
  static QickArray pos(no_of_particles, dim);
  for(int i=0; i!=no_of_particles; i++)
    for(int j=0; j!=dim; j++)
      pos(i,j) = getParticlePosition(i,j);
  return pos;
}

inline QickArray& Walker:: newPosVec(){
  static QickArray pos(no_of_particles, dim);
  for(int i=0; i!=no_of_particles; i++)
    for(int j=0; j!=dim; j++)
      pos(i,j) = getNewParticlePosition(i,j);
  return pos;
}

#endif
