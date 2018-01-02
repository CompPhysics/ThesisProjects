#include <cassert>
#include <iostream>
#include <cmath>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"
#include "squarewell.h"

using namespace std;

SquareWell::SquareWell(System* system, double V0, double distToWall, double alpha, double omega, bool analyticalKinetic, bool repulsion)
    : Hamiltonian(system, analyticalKinetic, alpha, omega){
    assert(omega > 0);
    m_omega = omega;
    m_repulsion = repulsion;
    m_numberOfDimensions = m_system->getNumberOfDimensions();
    m_distToWall = distToWall;
    m_V0 = V0;
    m_eigvals.load("../diagonalization/PlotAndData/Eigenvalues.dat", arma_ascii);
}

std::vector<double> SquareWell::computeLocalEnergy(std::vector<Particle*> particles) {
    /* Here, you need to compute the kinetic and potential energies. Note that
     * when using numerical differentiation, the computation of the kinetic
     * energy becomes the same for all Hamiltonians, and thus the code for
     * doing this should be moved up to the super-class, Hamiltonian.
     *
     * You may access the wave function currently used through the
     * getWaveFunction method in the m_system object in the super-class, i.e.
     * m_system->getWaveFunction()...
     */

    int numberOfParticles = m_system->getNumberOfParticles();

    double potentialEnergy = 0;
    double repulsiveTerm = 0;

    for (int i=0; i < numberOfParticles; i++){
        std::vector<double> r_i = particles[i]->getPosition();
        bool isOutside = false;

        for (int d = 0; d < m_numberOfDimensions; d++) {
            if (abs(r_i[d]) >= m_distToWall) { //MAKE SURE ABS RETURNS FLOAT!!!
                isOutside = true;
            }
        }
        if (isOutside) potentialEnergy += m_V0;

        for (int j=i+1; j < numberOfParticles; j++){
            double r_ijSquared = 0;
            std::vector<double> r_j = particles[j]->getPosition();
            for (int k=0; k < m_numberOfDimensions; k++){
                    r_ijSquared += (r_i[k] - r_j[k]) * (r_i[k] - r_j[k]);
            }

            double r_ij = sqrt(r_ijSquared);
            repulsiveTerm += 1./r_ij;
        }
    }
    //potentialEnergy *= 0.5*m_omega*m_omega;
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

double SquareWell::evaluateSingleParticleWF(vec n, std::vector<double> r, int j) {
    // Calculates the single particle wave function.

    //double alpha = m_parameters[0];
    double waveFunction = 1;//m_expFactor;

    for (int d = 0; d < m_numberOfDimensions; d++) {
        waveFunction *= sqrt(1./m_distToWall)*sin(((n[d])*M_PI/(2*m_distToWall))*(r[d]+m_distToWall));//computeHermitePolynomial(n[d], r[d]);
    }

//    double waveFunction = computeHermitePolynomial(nx, x)
//                         *computeHermitePolynomial(ny, y)
//                         *m_expFactor;//exp(-m_omega*alpha*(x*x + y*y)*0.5);


    return waveFunction;
}

std::vector<double> SquareWell::computeSPWFDerivative(vec n, std::vector<double> r, int j) {
    // Calculates the single particle wave function differentiated w.r.t. position.
    std::vector<double> derivative(m_numberOfDimensions);
    //double r2 = x*x + y*y;

    for (int d = 0; d < m_numberOfDimensions; d++) {
        derivative[d] = sqrt(1./m_distToWall)*((n[d])*M_PI/(2*m_distToWall))*cos(((n[d])*M_PI/(2*m_distToWall))*(r[d]+m_distToWall));/*(computeHermitePolynomialDerivative(n[d], r[d]) - m_alpha*m_omega*r[d]*computeHermitePolynomial(n[d], r[d]))
                       *m_expFactor;*/
        for (int dim = 0; dim < m_numberOfDimensions; dim++) {
            if (d != dim) derivative[d] *= sqrt(1./m_distToWall)*sin(((n[dim])*M_PI/(2*m_distToWall))*(r[dim]+m_distToWall));//computeHermitePolynomial(n[dim], r[dim]);
        }
    }

//    derivative[0] = (computeHermitePolynomialDerivative(nx, x) - alpha*m_omega*x*computeHermitePolynomial(nx, x))
//                   *computeHermitePolynomial(ny, y)*m_expFactor;//exp(-alpha*m_omega*r2*0.5);

//    derivative[1] = (computeHermitePolynomialDerivative(ny, y) - alpha*m_omega*y*computeHermitePolynomial(ny, y))
//                   *computeHermitePolynomial(nx, x)*m_expFactor;//exp(-alpha*m_omega*r2*0.5);

    return derivative;
}

double SquareWell::computeSPWFDoubleDerivative(vec n, std::vector<double> r, int j) {

    // Calculates the single particle wave function twice differentiated w.r.t. position.
    double doubleDerivative = 0;
    //double r2 = x*x + y*y;

    for (int d = 0; d < m_numberOfDimensions; d++) {

        double term =  -sqrt(1./m_distToWall)
                       *((n[d])*M_PI/(2*m_distToWall))
                       *((n[d])*M_PI/(2*m_distToWall))
                       *sin(((n[d])*M_PI/(2*m_distToWall))*(r[d]+m_distToWall));
                       /*computeHermitePolynomialDoubleDerivative(n[d], r[d])
                     - m_alpha*m_omega*computeHermitePolynomial(n[d], r[d])
                     - 2*m_alpha*m_omega*r[d]*computeHermitePolynomialDerivative(n[d], r[d])
                     + m_alpha*m_omega*m_alpha*m_omega*r[d]*r[d]*computeHermitePolynomial(n[d], r[d]);*/

        for (int dim = 0; dim < m_numberOfDimensions; dim++) {
            if (d != dim) term *= sqrt(1./m_distToWall)*sin(((n[dim])*M_PI/(2*m_distToWall))*(r[dim]+m_distToWall));//computeHermitePolynomial(n[dim], r[dim]);
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

double SquareWell::computeSPWFAlphaDerivative(vec n, std::vector<double> r, int j) {
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

double SquareWell::computeHermitePolynomial(int nValue, double position) {
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

double SquareWell::computeHermitePolynomialDerivative(int nValue, double position) {
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

double SquareWell::computeHermitePolynomialDoubleDerivative(int nValue, double position) {
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

double SquareWell::computeHermitePolynomialAlphaDerivative(int nValue, double position) {
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

//double SquareWell::evaluateSingleParticleWF(vec n, std::vector<double> r, int j) {
//    // Calculates the single particle wave function.

//    //double alpha = m_parameters[0];
//    double waveFunction = 1.;
//    double r2 = 0;

//    for (int d = 0; d < m_numberOfDimensions; d++) {
//        int n_d = n[d];
//        double E = m_eigvals.col(d)[n_d];
//        double k = sqrt(2.*m_eigvals.col(d)[n_d]);
//        double kPrime = sqrt(2*(E-m_V0));
//        double alpha = sqrt(2.*(m_V0-E));

//        double A;
//        double B;
//        double F;
//        double G;
//        double H;
//        double I;

//        if (n_d%2 == 0) {
//            A = 0.;
//            B = 1.;
//            F = 0.;
//            G = B*cos(k*m_distToWall)/exp(-alpha*m_distToWall);
//            H = G;
//            I = 0.;
//        }

//        else {
//            A = 1.;
//            B = 0.;
//            F = 0.;
//            G = -A*sin(k*m_distToWall)/exp(-alpha*m_distToWall);
//            H = -G;
//            I = 0.;
//        }


//        if ( r[d] < -m_distToWall ) {
//            if ( E > m_V0 ) {
//                waveFunction *= sin(kPrime*r[d]) + cos(kPrime*r[d]);
//            }
//            else {
//                waveFunction *= F*exp(-alpha*r[d]) + G*exp(alpha*r[d]);
//            }
//        }
//        else if ( r[d] > m_distToWall ) {
//            if ( E > m_V0 ) {
//                waveFunction *= sin(kPrime*r[d]) + cos(kPrime*r[d]);
//            }
//            else {
//                waveFunction *= H*exp(-alpha*r[d]) + I*exp(alpha*r[d]);
//            }
//        }
//        else {
//            waveFunction *= A*sin(k*r[d]) + B*cos(k*r[d]);
//        }
//    }

////    for (int d = 0; d < m_numberOfDimensions; d++) {
////        r2 += r[d]*r[d];
////    }

////    if (sqrt(r2) > m_distToWall) {
////        for (int d = 0; d < m_numberOfDimensions; d++) {
////            int n_d = n[d];
////            double E = m_eigvals.col(d)[n_d];
////            if (E > m_V0) {
////                double kPrime = sqrt(2*(E-m_V0));
////                waveFunction *= sin(kPrime*r[d]) + cos(kPrime*r[d]);
////            }
////            else if (r[d] < -m_distToWall) {
////                double alpha = sqrt(2.*(m_V0-E));
////                waveFunction *= exp(alpha*r[d]) + exp(-alpha*r[d]);
////            }
////            else if (r[d] > m_distToWall) {

////            }
////        }
////    }
////    else {
////        for (int d = 0; d < m_numberOfDimensions; d++) {
////            int n_d = n[d];
////            double k = sqrt(2.*m_eigvals.col(d)[n_d]);
////            waveFunction *= A*sin(k*r[d]) + B*cos(k*r[d]);
////        }
////    }

//    return waveFunction;
//}

//std::vector<double> SquareWell::computeSPWFDerivative(vec n, std::vector<double> r, int j) {
//    // Calculates the single particle wave function differentiated w.r.t. position.

//    std::vector<double> derivative(m_numberOfDimensions);
//    //double alpha = m_system->getWaveFunction()->getParameters()[0];

//    double r2 = 0;
//    vec alpha(m_numberOfDimensions);
//    vec k(m_numberOfDimensions);
//    vec kPrime(m_numberOfDimensions);
//    vec E(m_numberOfDimensions);

//    double A;
//    double B;
//    double F;
//    double G;
//    double H;
//    double I;

//    for (int d = 0; d < m_numberOfDimensions; d++) {
//        int n_d = n[d];
//        E[d] = m_eigvals.col(d)[n_d];
//        r2 += r[d]*r[d];
//        alpha[d] = sqrt(2.*(m_V0-E[d]));
//        kPrime[d] = sqrt(2.*(E[d]-m_V0));
//        k[d] = sqrt(2.*E[d]);
//    }

//    for (int d = 0; d < m_numberOfDimensions; d++) {

//        int n_d = n[d];

//        if (n_d%2 == 0) {
//            A = 0.;
//            B = 1.;
//            F = 0.;
//            G = B*cos(k[d]*m_distToWall)/exp(-alpha[d]*m_distToWall);
//            H = G;
//            I = 0.;
//        }

//        else {
//            A = 1.;
//            B = 0.;
//            F = 0.;
//            G = -A*sin(k[d]*m_distToWall)/exp(-alpha[d]*m_distToWall);
//            H = -G;
//            I = 0.;
//        }

//        if ( r[d] < -m_distToWall ) {
//            if ( E[d] > m_V0 ) {
//                derivative[d] = kPrime[d]*cos(kPrime[d]*r[d]) - kPrime[d]*sin(kPrime[d]*r[d]);

//                for (int dim = 0; dim < m_numberOfDimensions; dim++) {
//                    if (d != dim) derivative[d] *= sin(kPrime[dim]*r[dim]) + cos(kPrime[dim]*r[dim]);
//                }
//            }
//            else {
//                derivative[d] = alpha[d]*(-F*exp(-alpha[d]*r[d]) + G*exp(alpha[d]*r[d]));

//                for (int dim = 0; dim < m_numberOfDimensions; dim++) {
//                    if (d != dim) derivative[d] *= F*exp(-alpha[dim]*r[dim]) + G*exp(alpha[dim]*r[dim]);
//                }
//            }
//        }
//        else if ( r[d] > m_distToWall ) {
//            if ( E[d] > m_V0 ) {
//                derivative[d] = kPrime[d]*cos(kPrime[d]*r[d]) - kPrime[d]*sin(kPrime[d]*r[d]);

//                for (int dim = 0; dim < m_numberOfDimensions; dim++) {
//                    if (d != dim) derivative[d] *= sin(kPrime[dim]*r[dim]) + cos(kPrime[dim]*r[dim]);
//                }
//            }
//            else {
//                derivative[d] = alpha[d]*(-H*exp(-alpha[d]*r[d]) + I*exp(alpha[d]*r[d]));

//                for (int dim = 0; dim < m_numberOfDimensions; dim++) {
//                    if (d != dim) derivative[d] *= H*exp(-alpha[dim]*r[dim]) + I*exp(alpha[dim]*r[dim]);
//                }
//            }
//        }
//        else {
//            derivative[d] = k[d]*A*cos(k[d]*r[d]) - k[d]*B*sin(k[d]*r[d]);

//            for (int dim = 0; dim < m_numberOfDimensions; dim++) {
//                if (d != dim) derivative[d] *= A*sin(k[dim]*r[dim]) + B*cos(k[dim]*r[dim]);
//            }
//        }
//    }

////    if (sqrt(r2) > m_distToWall) {
////        for (int d = 0; d < m_numberOfDimensions; d++) {
////            if (E[d] > m_V0) {
////                derivative[d] = kPrime[d]*cos(kPrime[d]*r[d]) - kPrime[d]*sin(kPrime[d]*r[d]);

////                for (int dim = 0; dim < m_numberOfDimensions; dim++) {
////                    if (d != dim) derivative[d] *= sin(kPrime[dim]*r[dim]) + cos(kPrime[dim]*r[dim]);
////                }
////            }
////            else {
////                derivative[d] = alpha[d]*(exp(alpha[d]*r[d]) - exp(-alpha[d]*r[d]));

////                for (int dim = 0; dim < m_numberOfDimensions; dim++) {
////                    if (d != dim) derivative[d] *= exp(alpha[dim]*r[dim]) + exp(-alpha[dim]*r[dim]);
////                }
////            }
////        }
////    }

////    else {
////        for (int d = 0; d < m_numberOfDimensions; d++) {
////            derivative[d] = k[d]*cos(k[d]*r[d]) - k[d]*sin(k[d]*r[d]);

////            for (int dim = 0; dim < m_numberOfDimensions; dim++) {
////                if (d != dim) derivative[d] *= sin(k[dim]*r[dim]) + cos(k[dim]*r[dim]);
////            }
////        }
////    }

//    return derivative;
//}

//double SquareWell::computeSPWFDoubleDerivative(vec n, std::vector<double> r, int j) {

//    // Calculates the single particle wave function twice differentiated w.r.t. position.
//    double doubleDerivative = 0;
//    //double alpha = m_system->getWaveFunction()->getParameters()[0];

//    double r2 = 0;
//    vec alpha(m_numberOfDimensions);
//    vec k(m_numberOfDimensions);
//    vec kPrime(m_numberOfDimensions);
//    vec E(m_numberOfDimensions);

//    double A;
//    double B;
//    double F;
//    double G;
//    double H;
//    double I;

//    for (int d = 0; d < m_numberOfDimensions; d++) {
//        int n_d = n[d];
//        E[d] = m_eigvals.col(d)[n_d];
//        r2 += r[d]*r[d];
//        alpha[d] = sqrt(2.*(m_V0-E[d]));
//        kPrime[d] = sqrt(2.*(E[d]-m_V0));
//        k[d] = sqrt(2.*E[d]);
//    }

//    for (int d = 0; d < m_numberOfDimensions; d++) {

//        int n_d = n[d];

//        if (n_d%2 == 0) {
//            A = 0.;
//            B = 1.;
//            F = 0.;
//            G = B*cos(k[d]*m_distToWall)/exp(-alpha[d]*m_distToWall);
//            H = G;
//            I = 0.;
//        }

//        else {
//            A = 1.;
//            B = 0.;
//            F = 0.;
//            G = -A*sin(k[d]*m_distToWall)/exp(-alpha[d]*m_distToWall);
//            H = -G;
//            I = 0.;
//        }

//        double term = 0;
//        if ( r[d] < -m_distToWall ) {
//            if ( E[d] > m_V0 ) {
//                term = -kPrime[d]*kPrime[d]*sin(kPrime[d]*r[d]) - kPrime[d]*kPrime[d]*cos(kPrime[d]*r[d]);

//                for (int dim = 0; dim < m_numberOfDimensions; dim++) {
//                    if (d != dim) term *= sin(kPrime[dim]*r[dim]) + cos(kPrime[dim]*r[dim]);
//                }
//            }
//            else {
//                term = alpha[d]*alpha[d]*(F*exp(-alpha[d]*r[d]) + G*exp(alpha[d]*r[d]));

//                for (int dim = 0; dim < m_numberOfDimensions; dim++) {
//                    if (d != dim) term *= F*exp(-alpha[dim]*r[dim]) + G*exp(alpha[dim]*r[dim]);
//                }
//            }
//        }
//        else if ( r[d] > m_distToWall ) {
//            if ( E[d] > m_V0 ) {
//                term = -kPrime[d]*kPrime[d]*sin(kPrime[d]*r[d]) - kPrime[d]*kPrime[d]*cos(kPrime[d]*r[d]);

//                for (int dim = 0; dim < m_numberOfDimensions; dim++) {
//                    if (d != dim) term *= sin(kPrime[dim]*r[dim]) + cos(kPrime[dim]*r[dim]);
//                }
//            }
//            else {
//                term = alpha[d]*alpha[d]*(H*exp(-alpha[d]*r[d]) + I*exp(alpha[d]*r[d]));

//                for (int dim = 0; dim < m_numberOfDimensions; dim++) {
//                    if (d != dim) term *= H*exp(-alpha[dim]*r[dim]) + I*exp(alpha[dim]*r[dim]);
//                }
//            }
//        }
//        else {
//            double term = -k[d]*k[d]*A*sin(k[d]*r[d]) - k[d]*k[d]*B*cos(k[d]*r[d]);

//            for (int dim = 0; dim < m_numberOfDimensions; dim++) {
//                if (d != dim) term *= A*sin(k[dim]*r[dim]) + B*cos(k[dim]*r[dim]);
//            }
//        }
//        doubleDerivative += term;
//    }

////    if (sqrt(r2) > m_distToWall) {
////        for (int d = 0; d < m_numberOfDimensions; d++) {
////            double term;
////            if (E[d] > m_V0) {
////                term = -kPrime[d]*kPrime[d]*sin(kPrime[d]*r[d]) - kPrime[d]*kPrime[d]*cos(kPrime[d]*r[d]);

////                for (int dim = 0; dim < m_numberOfDimensions; dim++) {
////                    if (d != dim) term *= sin(kPrime[dim]*r[dim]) + cos(kPrime[dim]*r[dim]);
////                }
////            }
////            else {
////                term = alpha[d]*alpha[d]*(exp(alpha[d]*r[d]) + exp(-alpha[d]*r[d]));

////                for (int dim = 0; dim < m_numberOfDimensions; dim++) {
////                    if (d != dim) term *= exp(alpha[dim]*r[dim]) + exp(-alpha[dim]*r[dim]);
////                }
////            }
////            doubleDerivative += term;
////        }
////    }

////    else {
////        for (int d = 0; d < m_numberOfDimensions; d++) {
////            double term = -k[d]*k[d]*sin(k[d]*r[d]) - k[d]*k[d]*cos(k[d]*r[d]);

////            for (int dim = 0; dim < m_numberOfDimensions; dim++) {
////                if (d != dim) term *= sin(k[dim]*r[dim]) + cos(k[dim]*r[dim]);
////            }
////            doubleDerivative += term;
////        }
////    }

//    return doubleDerivative;

//}

//double SquareWell::computeSPWFAlphaDerivative(vec n, std::vector<double> r, int j) {
//    // Calculates the single particle wave function differentiated w.r.t. alpha.
//    //double derivative = 0;
//    //double alpha = m_system->getWaveFunction()->getParameters()[0];
//    //double r2 = x*x + y*y;

//    //return derivative;
//}

