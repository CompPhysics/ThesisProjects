#ifndef BASISFUNCTIONS_H
#define BASISFUNCTIONS_H
#include <vector>

class BasisFunctions
{
public:
    BasisFunctions();
    virtual double evaluate(std::vector<class Particle*> particles, int currentParticle) = 0;
};

#endif // BASISFUNCTIONS_H
