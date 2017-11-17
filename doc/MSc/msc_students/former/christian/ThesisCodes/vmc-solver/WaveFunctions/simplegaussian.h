#ifndef PROJECT2_SIMPLEGAUSSIAN_H
#define PROJECT2_SIMPLEGAUSSIAN_H
#include "wavefunction.h"

class SimpleGaussian : public WaveFunction {
public:
    SimpleGaussian(class System* system, double alpha);
    double evaluate(std::vector<class Particle*> particles);
    double computeDoubleDerivative(std::vector<class Particle*> particles);
    double computeMetropolisRatio(std::vector<Particle *> particles, int currentParticle,
                                  std::vector<double> positionChange);
    std::vector<double> computeDerivative(std::vector<class Particle*> particles);
    std::vector<double> computeDerivativeWrtParameters(std::vector<Particle *> particles);
};

#endif // PROJECT2_SIMPLEGAUSSIAN_H
