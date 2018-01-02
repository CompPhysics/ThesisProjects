#ifndef PROJECT2_TWOELECTRONS_H
#define PROJECT2_TWOELECTRONS_H
#include "wavefunction.h"

class TwoElectrons : public WaveFunction {
public:
    TwoElectrons(class System* system, double alpha, double beta, double omega,
                 double a, double C, bool Jastrow);
    double evaluate(std::vector<class Particle*> particles);
    double computeDoubleDerivative(std::vector<class Particle*> particles);
    double computeMetropolisRatio(std::vector<Particle *> particles, int currentParticle,
                                  std::vector<double> positionChange);
    std::vector<double> computeDerivative(std::vector<class Particle*> particles);
    std::vector<double> computeDerivativeWrtParameters(std::vector<Particle *> particles);

private:
    double m_omega = 0;
    double m_a = 0;
    double m_C = 0;
    bool m_Jastrow = false;
};

#endif // PROJECT2_TWOELECTRONS_H
