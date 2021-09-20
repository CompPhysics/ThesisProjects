#ifndef particle_hpp_IS_INCLUDED
#define particle_hpp_IS_INCLUDED

// class containing information with regards to a particle

//#include "QickArray/QickArray.hpp"

class Particle{
private:

  double *pos, *new_pos;
  int dim;

public:

  Particle(int dim_);
  ~Particle();

  double getPosition(int i); // returns i'th coord position
  double getNewPosition(int i); // returns i'th coord position
  void setPosition(int i, double x); // set's value of i'th coord in new_pos
  void updatePosition(); // pos=new_pos
  void resetPosition(); // new_pos=pos

  double operator()(int i){ return getPosition(i); }
  void   operator()(int i, double x){ setPosition(i,x); }
};

inline double Particle:: getPosition(int i){
  return pos[i];
}

inline double Particle:: getNewPosition(int i){
  return new_pos[i];
}

inline void Particle:: setPosition(int i, double x){
  new_pos[i]=x;
}

inline void Particle:: updatePosition(){
  for(int i=0; i!=dim; i++)
    pos[i]=new_pos[i];
}

inline void Particle:: resetPosition(){
  for(int i=0; i!=dim; i++)
    new_pos[i]=pos[i];
}

#endif
