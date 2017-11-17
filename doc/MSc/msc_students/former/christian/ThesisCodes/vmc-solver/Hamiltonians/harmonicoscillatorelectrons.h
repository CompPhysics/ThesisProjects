#ifndef PROJECT2_HARMONICOSCILLATORELECTRONS_H
#define PROJECT2_HARMONICOSCILLATORELECTRONS_H
#include "hamiltonian.h"
#include <vector>

class HarmonicOscillatorElectrons : public Hamiltonian {
public:
    HarmonicOscillatorElectrons(System* system, double alpha, double omega, bool analyticalKinetic, bool repulsion);
    std::vector<double> computeLocalEnergy(std::vector<Particle*> particles);
    double evaluateSingleParticleWF(vec n, std::vector<double> r, int j);
    std::vector<double> computeSPWFDerivative(vec n, std::vector<double> r, int j);
    double computeSPWFDoubleDerivative(vec n, std::vector<double> r, int j);
    double computeSPWFAlphaDerivative(vec n, std::vector<double> r, int j);
    double computeHermitePolynomial(int nValue, double position);
    double computeHermitePolynomialDerivative(int nValue, double position);
    double computeHermitePolynomialDoubleDerivative(int nValue, double position);
    double computeHermitePolynomialAlphaDerivative(int nValue, double position);

private:
    int m_numberOfDimensions = 0;
    double m_omega = 0;
    bool m_repulsion = false;
};

#endif // PROJECT2_HARMONICOSCILLATORELECTRONS_H
