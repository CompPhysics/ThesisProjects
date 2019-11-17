#ifndef SQUAREWELL_H
#define SQUAREWELL_H
#include "hamiltonian.h"
#include <vector>
#include <armadillo>

using namespace arma;

class SquareWell : public Hamiltonian {
public:
    SquareWell(System* system, double V0, double distToWall, double alpha, double omega, bool analyticalKinetic, bool repulsion);
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
    double m_V0 = 0;
    double m_distToWall = 0;
    double m_omega = 0;
    bool m_repulsion = false;
    mat m_eigvals;
};


#endif // SQUAREWELL_H
