#include "harmonicoscillatorelectrons.h"
#include <cassert>
#include <iostream>
#include <cmath>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

using namespace std;

HarmonicOscillatorElectrons::HarmonicOscillatorElectrons(System* system, double alpha, double omega,
                                                         bool analyticalKinetic, bool repulsion) :
    Hamiltonian(system, analyticalKinetic, alpha, omega) {
    assert(omega > 0);
    m_omega = omega;
    m_repulsion = repulsion;
    m_numberOfDimensions = m_system->getNumberOfDimensions();
}

std::vector<double> HarmonicOscillatorElectrons::computeLocalEnergy(std::vector<Particle*> particles) {

    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();

    double potentialEnergy = 0;
    double repulsiveTerm = 0;

    for (int i=0; i < numberOfParticles; i++){
        double rSquared = 0;
        std::vector<double> r_i = particles[i]->getPosition();
        for (int k=0; k < numberOfDimensions; k++){
            rSquared += r_i[k]*r_i[k];
        }
        potentialEnergy += rSquared;

        for (int j=i+1; j < numberOfParticles; j++){
            double r_ijSquared = 0;
            std::vector<double> r_j = particles[j]->getPosition();
            for (int k=0; k < numberOfDimensions; k++){
                    r_ijSquared += (r_i[k] - r_j[k]) * (r_i[k] - r_j[k]);
            }

            double r_ij = sqrt(r_ijSquared);
            repulsiveTerm += 1./r_ij;
        }
    }
    potentialEnergy *= 0.5*m_omega*m_omega;
    if (m_repulsion) { potentialEnergy += repulsiveTerm; }

    double kineticEnergy = 0;

    if (m_analyticalKinetic == true){
        // Using analytical expression.
        double doubleDerivative = m_system->getWaveFunction()->computeDoubleDerivative(particles);
        kineticEnergy = -0.5*doubleDerivative;
    }
    else{
        // Using numerical diff.
        kineticEnergy = computeKineticEnergy(particles);
    }

    std::vector<double> energies(3);
    energies[0] = kineticEnergy + potentialEnergy;
    energies[1] = kineticEnergy;
    energies[2] = potentialEnergy;

    return energies;
}

double HarmonicOscillatorElectrons::evaluateSingleParticleWF(vec n, std::vector<double> r, int j) {
    // Calculates the single particle wave function.

    //double alpha = m_parameters[0];
    double waveFunction = m_expFactor;

    for (int d = 0; d < m_numberOfDimensions; d++) {
        int nd = n[d];
        waveFunction *= m_hermitePolynomials[nd]->eval(r[d]);//computeHermitePolynomial(n[d], r[d]);
    }

//    double waveFunction = computeHermitePolynomial(nx, x)
//                         *computeHermitePolynomial(ny, y)
//                         *m_expFactor;//exp(-m_omega*alpha*(x*x + y*y)*0.5);

    return waveFunction;
}

std::vector<double> HarmonicOscillatorElectrons::computeSPWFDerivative(vec n, std::vector<double> r, int j) {
    // Calculates the single particle wave function differentiated w.r.t. position.
    std::vector<double> derivative(m_numberOfDimensions);
    //double r2 = x*x + y*y;

    for (int d = 0; d < m_numberOfDimensions; d++) {

        int nd = n[d];

//        derivative[d] = (computeHermitePolynomialDerivative(n[d], r[d]) - m_alpha*m_omega*r[d]*computeHermitePolynomial(n[d], r[d]))
//                       *m_expFactor;

        derivative[d] = (m_hermitePolynomialsDerivative[nd]->eval(r[d]) - m_alpha*m_omega*r[d]*m_hermitePolynomials[nd]->eval(r[d]))
                       *m_expFactor;

        for (int dim = 0; dim < m_numberOfDimensions; dim++) {
            int ndim = n[dim];
            if (d != dim) derivative[d] *= m_hermitePolynomials[ndim]->eval(r[dim]);//computeHermitePolynomial(n[dim], r[dim]);
        }
    }

//    derivative[0] = (computeHermitePolynomialDerivative(nx, x) - alpha*m_omega*x*computeHermitePolynomial(nx, x))
//                   *computeHermitePolynomial(ny, y)*m_expFactor;//exp(-alpha*m_omega*r2*0.5);

//    derivative[1] = (computeHermitePolynomialDerivative(ny, y) - alpha*m_omega*y*computeHermitePolynomial(ny, y))
//                   *computeHermitePolynomial(nx, x)*m_expFactor;//exp(-alpha*m_omega*r2*0.5);

    return derivative;
}

double HarmonicOscillatorElectrons::computeSPWFDoubleDerivative(vec n, std::vector<double> r, int j) {

    // Calculates the single particle wave function twice differentiated w.r.t. position.
    double doubleDerivative = 0;
    //double r2 = x*x + y*y;

    for (int d = 0; d < m_numberOfDimensions; d++) {
        int nd = n[d];

//        double term =  computeHermitePolynomialDoubleDerivative(n[d], r[d])
//                     - m_alpha*m_omega*computeHermitePolynomial(n[d], r[d])
//                     - 2*m_alpha*m_omega*r[d]*computeHermitePolynomialDerivative(n[d], r[d])
//                     + m_alpha*m_omega*m_alpha*m_omega*r[d]*r[d]*computeHermitePolynomial(n[d], r[d]);

        double term =  m_hermitePolynomialsDoubleDerivative[nd]->eval(r[d])
                     - m_alpha*m_omega*m_hermitePolynomials[nd]->eval(r[d])
                     - 2*m_alpha*m_omega*r[d]*m_hermitePolynomialsDerivative[nd]->eval(r[d])
                     + m_alpha*m_omega*m_alpha*m_omega*r[d]*r[d]*m_hermitePolynomials[nd]->eval(r[d]);


        for (int dim = 0; dim < m_numberOfDimensions; dim++) {
            int ndim = n[dim];
            if (d != dim) term *= m_hermitePolynomials[ndim]->eval(r[dim]);//computeHermitePolynomial(n[dim], r[dim]);
        }

        doubleDerivative += term;
    }

    doubleDerivative *= m_expFactor;

//    doubleDerivative += computeHermitePolynomial(ny, y)*m_expFactor//exp(-alpha*m_omega*r2*0.5)
//                       *(computeHermitePolynomialDoubleDerivative(nx, x)
//                         - alpha*m_omega*computeHermitePolynomial(nx, x)
//                         - 2*alpha*m_omega*x*computeHermitePolynomialDerivative(nx, x)
//                         + alpha*m_omega*alpha*m_omega*x*x*computeHermitePolynomial(nx, x));

//    doubleDerivative += computeHermitePolynomial(nx, x)*m_expFactor//exp(-alpha*m_omega*r2*0.5)
//                       *(computeHermitePolynomialDoubleDerivative(ny, y)
//                         - alpha*m_omega*computeHermitePolynomial(ny, y)
//                         - 2*alpha*m_omega*y*computeHermitePolynomialDerivative(ny, y)
//                         + alpha*m_omega*alpha*m_omega*y*y*computeHermitePolynomial(ny, y));

    return doubleDerivative;

}

double HarmonicOscillatorElectrons::computeSPWFAlphaDerivative(vec n, std::vector<double> r, int j) {
    // Calculates the single particle wave function differentiated w.r.t. alpha.
    double derivative = 0;
    //double expFactor = m_system->getWaveFunction()->getExpFactor();
    double alpha = m_system->getWaveFunction()->getParameters()[0];
    //double r2 = x*x + y*y;
    double r2 = 0;
    for (int d = 0; d < m_numberOfDimensions; d++) {
        r2 += r[d]*r[d];
    }

    double term1 = -0.5*m_omega*r2;

    for (int d = 0; d < m_numberOfDimensions; d++) {

        term1 *= computeHermitePolynomial(n[d], r[d]);
        double otherTerms = computeHermitePolynomialAlphaDerivative(n[d], r[d]);

        for (int dim = 0; dim < m_numberOfDimensions; dim++) {
            if (d != dim) otherTerms *= computeHermitePolynomial(n[dim], r[dim]);
        }
        derivative += otherTerms;
    }

    derivative += term1;
    derivative *= exp(-0.5*alpha*m_omega*r2);

//    derivative += (-0.5*m_omega*r2*computeHermitePolynomial(nx, x)*computeHermitePolynomial(ny, y)
//                   +computeHermitePolynomialAlphaDerivative(nx, x)*computeHermitePolynomial(ny, y)
//                   +computeHermitePolynomialAlphaDerivative(ny, y)*computeHermitePolynomial(nx, x))
//                 *exp(-0.5*alpha*m_omega*r2);

    return derivative;

}

double HarmonicOscillatorElectrons::computeHermitePolynomial(int nValue, double position) {
    // Computes Hermite polynomials.
    double alphaSqrt = sqrt(m_alpha);
    double omegaSqrt = sqrt(m_omega);
    double factor = 2*alphaSqrt*omegaSqrt*position;

    double HermitePolynomialPP = 0;                 // H_{n-2}
    double HermitePolynomialP = 1;                  // H_{n-1}
    double HermitePolynomial = HermitePolynomialP;  // H_n

    for (int n=1; n <= nValue; n++) {
        HermitePolynomial = factor*HermitePolynomialP - 2*(n-1)*HermitePolynomialPP;
        HermitePolynomialPP = HermitePolynomialP;
        HermitePolynomialP = HermitePolynomial;
    }

    return HermitePolynomial;

}

double HarmonicOscillatorElectrons::computeHermitePolynomialDerivative(int nValue, double position) {
    // Computes Hermite polynomials differentiated w.r.t. position.
    double alphaSqrt = sqrt(m_alpha);
    double omegaSqrt = sqrt(m_omega);
    double factor1 = 2*alphaSqrt*omegaSqrt;
    double factor2 = 2*alphaSqrt*omegaSqrt*position;

    double HPDerivativePP = 0;              // d/dx H_{n-2}
    double HPDerivativeP = 0;               // d/dx H_{n-1}
    double HPDerivative = HPDerivativeP;    // d/dx H_n

    for (int n=1; n <= nValue; n++) {
        HPDerivative = factor1*computeHermitePolynomial(n-1, position)
                      +factor2*HPDerivativeP
                      -2*(n-1)*HPDerivativePP;
        HPDerivativePP = HPDerivativeP;
        HPDerivativeP = HPDerivative;
    }

    return HPDerivative;

}

double HarmonicOscillatorElectrons::computeHermitePolynomialDoubleDerivative(int nValue, double position) {
    // Computes Hermite polynomials twice differentiated w.r.t. position.
    double alphaSqrt = sqrt(m_alpha);
    double omegaSqrt = sqrt(m_omega);
    double factor1 = 4*alphaSqrt*omegaSqrt;
    double factor2 = 2*alphaSqrt*omegaSqrt*position;

    double HPDoubleDerivativePP = 0;                    // d/dx d/dx H_{n-2}
    double HPDoubleDerivativeP = 0;                     // d/dx d/dx H_{n-1}
    double HPDoubleDerivative = HPDoubleDerivativeP;    // d/dx d/dx H_n

    for (int n=1; n <= nValue; n++) {
        HPDoubleDerivative = factor1*computeHermitePolynomialDerivative(n-1, position)
                            +factor2*HPDoubleDerivativeP
                            -2*(n-1)*HPDoubleDerivativePP;
        HPDoubleDerivativePP = HPDoubleDerivativeP;
        HPDoubleDerivativeP = HPDoubleDerivative;
    }

    return HPDoubleDerivative;

}

double HarmonicOscillatorElectrons::computeHermitePolynomialAlphaDerivative(int nValue, double position) {
    // Computes Hermite polynomials differentiated w.r.t. alpha.
    double alpha = m_system->getWaveFunction()->getParameters()[0];
    double alphaSqrt = sqrt(alpha);
    double omegaSqrt = sqrt(m_omega);
    double factor1 = omegaSqrt/alphaSqrt*position;
    double factor2 = 2*alphaSqrt*omegaSqrt*position;

    double HPDerivativePP = 0;              // d/dα H_{n-2}
    double HPDerivativeP = 0;               // d/dα H_{n-1}
    double HPDerivative = HPDerivativeP;    // d/dα H_n

    for (int n=1; n <= nValue; n++) {
        HPDerivative = factor1*computeHermitePolynomial(n-1, position)
                      +factor2*HPDerivativeP
                      -2*(n-1)*HPDerivativePP;
        HPDerivativePP = HPDerivativeP;
        HPDerivativeP = HPDerivative;
    }

    return HPDerivative;
}
