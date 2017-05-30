#include "doubleharmonicoscillator.h"
#include <cassert>
#include <iostream>
#include <cmath>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

using namespace std;

DoubleHarmonicOscillator::DoubleHarmonicOscillator(System *system, vec L, double alpha, double omega, bool analyticalKinetic, bool repulsion) :
    Hamiltonian(system, analyticalKinetic, alpha, omega) {
    m_alpha = alpha;
    assert(omega > 0);
    m_omega = omega;
    m_repulsion = repulsion;
    m_L = L;
    m_numberOfDimensions = m_system->getNumberOfDimensions();
}

std::vector<double> DoubleHarmonicOscillator::computeLocalEnergy(std::vector<Particle*> particles) {

    int numberOfParticles = m_system->getNumberOfParticles();
    int numberOfDimensions = m_system->getNumberOfDimensions();

    double potentialEnergy = 0;
    double repulsiveTerm = 0;

    for (int i=0; i < numberOfParticles; i++){
        double rSquared = 0;
        double term2 = 0;
        double term3 = 0;

        std::vector<double> r_i = particles[i]->getPosition();
        for (int k=0; k < numberOfDimensions; k++){
            rSquared += r_i[k]*r_i[k];
            term2 += 2.*abs(r_i[k])*m_L(k);
            term3 += m_L(k)*m_L(k);
        }

        potentialEnergy += rSquared - term2 + term3;

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

double DoubleHarmonicOscillator::evaluateSingleParticleWF(vec n, std::vector<double> r, int j) {
    // Calculates the single particle wave function.

    //double alpha = m_parameters[0];
    std::vector<double> r_p(m_numberOfDimensions);
    std::vector<double> r_m(m_numberOfDimensions);
    double r_p2 = 0;
    double r_m2 = 0;

    for (int d = 0; d < m_numberOfDimensions; d++) {
        r_p[d] = r[d] + m_L(d);
        r_m[d] = r[d] - m_L(d);

        r_p2 += r_p[d]*r_p[d];
        r_m2 += r_m[d]*r_m[d];
    }

    double waveFunction_p = exp(-m_alpha*m_omega*r_p2*0.5);

    for (int d = 0; d < m_numberOfDimensions; d++) {
        int nd = n[d];
        waveFunction_p *= m_hermitePolynomials[nd]->eval(r_p[d]);//computeHermitePolynomial(n[d], r_p[d]);
    }

    double waveFunction_m = exp(-m_alpha*m_omega*r_m2*0.5);

    for (int d = 0; d < m_numberOfDimensions; d++) {
        int nd = n[d];
        waveFunction_m *= m_hermitePolynomials[nd]->eval(r_m[d]);//computeHermitePolynomial(n[d], r_m[d]);
    }

    int sign = -2*(j%2)+1;
    //cout << sign << endl;

    if (m_system->getNumberOfParticles() == 1) { waveFunction_p = 0; }
    double waveFunction = waveFunction_p + sign*waveFunction_m;

//    double waveFunction = computeHermitePolynomial(nx, x)
//                         *computeHermitePolynomial(ny, y)
//                         *m_expFactor;//exp(-m_omega*alpha*(x*x + y*y)*0.5);

    return waveFunction;
}

std::vector<double> DoubleHarmonicOscillator::computeSPWFDerivative(vec n, std::vector<double> r, int j) {
    // Calculates the single particle wave function differentiated w.r.t. position.
    std::vector<double> derivative(m_numberOfDimensions);
    //double r2 = x*x + y*y;

    std::vector<double> r_p(m_numberOfDimensions);
    std::vector<double> r_m(m_numberOfDimensions);
    double r_p2 = 0;
    double r_m2 = 0;

    for (int d = 0; d < m_numberOfDimensions; d++) {
        r_p[d] = r[d] + m_L(d);
        r_m[d] = r[d] - m_L(d);

        r_p2 += r_p[d]*r_p[d];
        r_m2 += r_m[d]*r_m[d];
    }

    double expFactor_p = exp(-m_alpha*m_omega*r_p2*0.5);

    std::vector<double> derivative_p(m_numberOfDimensions);

    for (int d = 0; d < m_numberOfDimensions; d++) {
        int nd = n[d];
//        derivative_p[d] = (computeHermitePolynomialDerivative(n[d], r_p[d]) - m_alpha*m_omega*r_p[d]*computeHermitePolynomial(n[d], r_p[d]))
//                       *expFactor_p;
        derivative_p[d] = (m_hermitePolynomialsDerivative[nd]->eval(r_p[d]) - m_alpha*m_omega*r_p[d]*m_hermitePolynomials[nd]->eval(r_p[d]))
                       *expFactor_p;
        for (int dim = 0; dim < m_numberOfDimensions; dim++) {
            int ndim = n[dim];
            if (d != dim) derivative_p[d] *= m_hermitePolynomials[ndim]->eval(r_p[dim]);//computeHermitePolynomial(n[dim], r_p[dim]);
        }
    }

    double expFactor_m = exp(-m_alpha*m_omega*r_m2*0.5);
    std::vector<double> derivative_m(m_numberOfDimensions);

    for (int d = 0; d < m_numberOfDimensions; d++) {
        int nd = n[d];
//        derivative_m[d] = (computeHermitePolynomialDerivative(n[d], r_m[d]) - m_alpha*m_omega*r_m[d]*computeHermitePolynomial(n[d], r_m[d]))
//                       *expFactor_m;
        derivative_m[d] = (m_hermitePolynomialsDerivative[nd]->eval(r_m[d]) - m_alpha*m_omega*r_m[d]*m_hermitePolynomials[nd]->eval(r_m[d]))
                       *expFactor_m;
        for (int dim = 0; dim < m_numberOfDimensions; dim++) {
            int ndim = n[dim];
            if (d != dim) derivative_m[d] *= m_hermitePolynomials[ndim]->eval(r_m[dim]);//computeHermitePolynomial(n[dim], r_m[dim]);
        }
    }

    int sign = -2*(j%2)+1;

    if (m_system->getNumberOfParticles() == 1) {
        for (int d = 0; d < m_numberOfDimensions; d++) {
            derivative_p[d] = 0;
        }
    }

    for (int d = 0; d < m_numberOfDimensions; d++) {
        derivative[d] = derivative_p[d] + sign*derivative_m[d];
    }

//    derivative[0] = (computeHermitePolynomialDerivative(nx, x) - alpha*m_omega*x*computeHermitePolynomial(nx, x))
//                   *computeHermitePolynomial(ny, y)*m_expFactor;//exp(-alpha*m_omega*r2*0.5);

//    derivative[1] = (computeHermitePolynomialDerivative(ny, y) - alpha*m_omega*y*computeHermitePolynomial(ny, y))
//                   *computeHermitePolynomial(nx, x)*m_expFactor;//exp(-alpha*m_omega*r2*0.5);

    return derivative;
}

double DoubleHarmonicOscillator::computeSPWFDoubleDerivative(vec n, std::vector<double> r, int j) {

    // Calculates the single particle wave function twice differentiated w.r.t. position.
    double doubleDerivative = 0;
    //double r2 = x*x + y*y;

    std::vector<double> r_p(m_numberOfDimensions);
    std::vector<double> r_m(m_numberOfDimensions);
    double r_p2 = 0;
    double r_m2 = 0;

    for (int d = 0; d < m_numberOfDimensions; d++) {
        r_p[d] = r[d] + m_L(d);
        r_m[d] = r[d] - m_L(d);

        r_p2 += r_p[d]*r_p[d];
        r_m2 += r_m[d]*r_m[d];
    }

    double doubleDerivative_p = 0;
    double expFactor_p = exp(-m_alpha*m_omega*r_p2*0.5);

    for (int d = 0; d < m_numberOfDimensions; d++) {

        int nd = n[d];

//        double term =  computeHermitePolynomialDoubleDerivative(n[d], r_p[d])
//                     - m_alpha*m_omega*computeHermitePolynomial(n[d], r_p[d])
//                     - 2*m_alpha*m_omega*r_p[d]*computeHermitePolynomialDerivative(n[d], r_p[d])
//                     + m_alpha*m_omega*m_alpha*m_omega*r_p[d]*r_p[d]*computeHermitePolynomial(n[d], r_p[d]);

        double term =  m_hermitePolynomialsDoubleDerivative[nd]->eval(r_p[d])
                     - m_alpha*m_omega*m_hermitePolynomials[nd]->eval(r_p[d])
                     - 2*m_alpha*m_omega*r_p[d]*m_hermitePolynomialsDerivative[nd]->eval(r_p[d])
                     + m_alpha*m_omega*m_alpha*m_omega*r_p[d]*r_p[d]*m_hermitePolynomials[nd]->eval(r_p[d]);

        for (int dim = 0; dim < m_numberOfDimensions; dim++) {
            int ndim = n[dim];
            if (d != dim) term *= m_hermitePolynomials[ndim]->eval(r_p[dim]);//computeHermitePolynomial(n[dim], r_p[dim]);
        }

        doubleDerivative_p += term;
    }

    doubleDerivative_p *= expFactor_p;

    double doubleDerivative_m = 0;
    double expFactor_m = exp(-m_alpha*m_omega*r_m2*0.5);

    for (int d = 0; d < m_numberOfDimensions; d++) {

        int nd = n[d];

//        double term =  computeHermitePolynomialDoubleDerivative(n[d], r_m[d])
//                     - m_alpha*m_omega*computeHermitePolynomial(n[d], r_m[d])
//                     - 2*m_alpha*m_omega*r_m[d]*computeHermitePolynomialDerivative(n[d], r_m[d])
//                     + m_alpha*m_omega*m_alpha*m_omega*r_m[d]*r_m[d]*computeHermitePolynomial(n[d], r_m[d]);

        double term =  m_hermitePolynomialsDoubleDerivative[nd]->eval(r_m[d])
                     - m_alpha*m_omega*m_hermitePolynomials[nd]->eval(r_m[d])
                     - 2*m_alpha*m_omega*r_m[d]*m_hermitePolynomialsDerivative[nd]->eval(r_m[d])
                     + m_alpha*m_omega*m_alpha*m_omega*r_m[d]*r_m[d]*m_hermitePolynomials[nd]->eval(r_m[d]);

        for (int dim = 0; dim < m_numberOfDimensions; dim++) {
            int ndim = n[dim];
            if (d != dim) term *= m_hermitePolynomials[ndim]->eval(r_m[dim]);//computeHermitePolynomial(n[dim], r_m[dim]);
        }

        doubleDerivative_m += term;
    }

    doubleDerivative_m *= expFactor_m;

    int sign = -2*(j%2)+1;

    if (m_system->getNumberOfParticles() == 1) { doubleDerivative_p = 0; }

    doubleDerivative = doubleDerivative_p + sign*doubleDerivative_m;

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

double DoubleHarmonicOscillator::computeSPWFAlphaDerivative(vec n, std::vector<double> r, int j) {
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

double DoubleHarmonicOscillator::computeHermitePolynomial(int nValue, double position) {
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

double DoubleHarmonicOscillator::computeHermitePolynomialDerivative(int nValue, double position) {
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

double DoubleHarmonicOscillator::computeHermitePolynomialDoubleDerivative(int nValue, double position) {
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

double DoubleHarmonicOscillator::computeHermitePolynomialAlphaDerivative(int nValue, double position) {
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
