#ifndef PROJECT2_RANDOMUNIFORM_H
#define PROJECT2_RANDOMUNIFORM_H
#include "initialstate.h"

class RandomUniform : public InitialState {
public:
    RandomUniform(System* system, int numberOfDimensions, int numberOfParticles, int my_rank);
    void setupInitialState();
};

#endif // PROJECT2_RANDOMUNIFORM_H
