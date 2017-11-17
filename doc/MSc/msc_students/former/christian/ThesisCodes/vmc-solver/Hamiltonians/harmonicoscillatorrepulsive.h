#ifndef PROJECT2_HARMONICOSCILLATORREPULSIVE_H
#define PROJECT2_HARMONICOSCILLATORREPULSIVE_H
#include "hamiltonian.h"
#include <vector>

class HarmonicOscillatorRepulsive : public Hamiltonian {
public:
    HarmonicOscillatorRepulsive(System* system, double alpha, double omega, double a, double gamma, bool analyticalKinetic);
    std::vector<double> computeLocalEnergy(std::vector<Particle*> particles);

private:
    double m_omega = 0;
    double m_a = 0;
    double m_gamma = 0;
};

#endif // PROJECT2_HARMONICOSCILLATORREPULSIVE_H
