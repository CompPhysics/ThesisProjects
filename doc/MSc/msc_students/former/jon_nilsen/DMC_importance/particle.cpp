#include "particle.hpp"

using namespace std;

Particle:: Particle(int dim_){
  dim=dim_;
  pos = new double[dim];
  new_pos = new double[dim];
}

Particle:: ~Particle(){
  delete[] pos;
  delete[] new_pos;
}
